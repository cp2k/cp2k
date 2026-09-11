!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using RI-RS Approximation for molecules
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************

MODULE gw_non_periodic_ri_rs
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_data_size, dbcsr_get_info, &
        dbcsr_get_occupation, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set,&
                                              max_elements_per_block
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_uplo_to_full
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE gw_compute_Z_lP,                 ONLY: compute_Z_lP
   USE gw_fm_utils,                     ONLY: fm_contract_ABA,&
                                              fm_invert,&
                                              fm_sqrt
   USE gw_large_cell_gamma,             ONLY: Fourier_transform_w_to_t,&
                                              G_occ_vir,&
                                              compute_fm_chi_Gamma_freq,&
                                              create_fm_W_MIC_time,&
                                              delete_unnecessary_files,&
                                              fill_fm_Sigma_c_Gamma_time,&
                                              fm_write
   USE gw_optimize_ri_rs_grid,          ONLY: optimize_ri_rs_grid
   USE gw_utils,                        ONLY: analyt_conti_and_print,&
                                              de_init_bs_env,&
                                              time_to_freq
   USE input_constants,                 ONLY: G0W0,&
                                              evGW0,&
                                              rtp_method_bse
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8,&
                                              max_line_length
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE machine,                         ONLY: m_flush,&
                                              m_hostnm,&
                                              m_memory_details,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
   USE post_scf_bandstructure_utils,    ONLY: get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE rirs_grid_utils,                 ONLY: evaluate_ao_on_points,&
                                              get_rirs_grid_filepath,&
                                              initialize_rirs_grid,&
                                              read_rirs_grid_file
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_non_periodic_ri_rs'

   PUBLIC :: gw_calc_non_periodic_ri_rs, ri_rs_grid_assembler, &
             get_basis_offsets, precompute_ri_rs_radii, &
             atomic_basis_at_grid_point, &
             reserve_blocks_within_radius

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for molecules
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE gw_calc_non_periodic_ri_rs(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_non_periodic_ri_rs'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time

      CALL timeset(routineN, handle)

      ! ==============================================================================
      ! 0. Precompute AO and RI radii
      !    Per-atom cutoff radii from the most diffuse Gaussian primitives of
      !    the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets:
      !      α_min,ao = min { ζ_ao | ζ_ao > 10⁻³ },  α_min,ri analogous
      !      r_ao = sqrt( -ln(ε) / α_min,ao )   (radius_ao_per_atom)
      !      r_ri = sqrt( -ln(ε) / α_min,ri )   (radius_ri_per_atom)
      ! ==============================================================================
      CALL precompute_ri_rs_radii(qs_env, bs_env)

      ! ==============================================================================
      ! 1. Grid generation for RI-RS
      !    Modified Lebedev atomic grids (Duchemin & Blase), one per atom,
      !    concatenated into a flat global list:  r_l = R_A + r_l^(A)
      ! ==============================================================================
      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)

      ! ==============================================================================
      ! 2a. Atomic basis evaluation on the grid (grid x AO matrix)
      !     Φ_μl = Φ_μ(r_l)   (mat_phi_mu_l)
      ! ==============================================================================
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)

      ! ==============================================================================
      ! 2b. Print the memory estimate for the RI-RS calculation
      ! ==============================================================================
      CALL print_ri_rs_memory_estimate(qs_env, bs_env)

      ! ==============================================================================
      ! 3. RI-RS fitting coefficients Z_lP (grid x RI matrix)
      !    Per-atom regularized solve, restricted to grid points r_l within a
      !    cutoff distance of atom P:
      !      a. D_ll' = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]²
      !      b. D_lP  = Σ_μν Φ_μ(r_l) Φ_ν(r_l) (μν|P)
      !      c. Jacobi conditioning with d_l = 1/sqrt(D_ll):
      !         D'_ll' = d_l D_ll' d_l' + λδ_ll' ,  D'_lP = d_l D_lP
      !      d. Solve Σ_l' D'_ll' Z'_l'P = D'_lP
      !      e. Rescale Z_lP = d_l Z'_lP                        (mat_Z_lP)
      ! ==============================================================================
      CALL compute_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                        bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! flag the RI-RS grid as built so a subsequent RT-BSE run reuses Z_lP
      ! instead of rebuilding it
      bs_env%ri_rs%grid_built = .TRUE.

      ! ==============================================================================
      ! 4. Polarizability matrix χ on the imaginary-time grid
      !      G^occ_µλ(i|τ|)  = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^vir_µλ(i|τ|)  = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^occ_ll'(i|τ|) = Σ_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')  (G^vir analogous)
      !      χ_ll'(iτ)       = G^occ_ll'(i|τ|) ∘ G^vir_ll'(i|τ|) (element-wise)
      !      χ_PQ(iτ)        = Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q
      ! ==============================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! ==============================================================================
      ! 5. Screened Coulomb interaction W (RI basis)
      !      χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W_PQ(iτ)
      ! ==============================================================================
      CALL compute_W(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      ! ==============================================================================
      ! 6. Exact-exchange self-energy Σ^x
      !      D_µν          = Σ_n^occ C_µn C_νn                (density matrix)
      !      D_ll'         = Σ_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
      !      V^tr_ll'      = Σ_PQ Z_lP V^tr_PQ Z_l'Q          (truncated Coulomb)
      !      Σ^x_ll'       = D_ll' ∘ V^tr_ll'
      !      Σ^x_λσ(k=0)   = -Σ_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
      ! ==============================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

      ! ==============================================================================
      ! 7. Correlation self-energy Σ^c and quasiparticle energies, iterated
      !    until eigenvalue self-consistency if the &EVGW0 section is given and
      !    done in a single pass for G0W0.
      !
      !    (a) W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
      !        Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
      !        Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
      !        Σ^c_λσ(iτ)  =  Σ_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l')
      !    (b) Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
      !        ϵ_n^GW,(i) = ϵ_n^DFT + Σ^c_nn[G^(i-1),W](ϵ_n^GW,(i)) + Σ^x_nn - v^xc_nn
      ! ==============================================================================
      CALL compute_Sigma_c_and_QP_energies(bs_env, fm_W_time, fm_Sigma_x_Gamma)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_non_periodic_ri_rs

! **************************************************************************************************
!> \brief Correlation self-energy and quasiparticle energies, iterated to eigenvalue
!>        self-consistency in G (evGW0).
!>
!>        For G0W0 this runs once with G^(0) built from the DFT eigenvalues. For evGW0 the
!>        Green's function is rebuilt from the quasiparticle energies of all states until the
!>        quasiparticle HOMO, LUMO and HOMO-LUMO gap change by less than EPS_ITER between two
!>        cycles, or MAX_ITER cycles are spent. W, Σ^x and everything computed before this
!>        routine stay frozen: they either do not depend on the eigenvalues at all (Σ^x is G
!>        at τ = 0, a pure density matrix) or are held fixed by construction in GW0. The
!>        imaginary-time grid is fixed as well, since W(iτ) lives on it.
!>
!> \param bs_env ...
!> \param fm_W_time ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_and_QP_energies(bs_env, fm_W_time, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time, fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_and_QP_energies'

      INTEGER                                            :: handle, i_iter, n_iter
      LOGICAL                                            :: converged
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2)                        :: e_fermi_dft
      REAL(KIND=dp), DIMENSION(3, 2)                     :: band_prev
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      converged = .FALSE.

      CALL init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      ! evGW0 self-consistency loop; for G0W0, loop is terminated after one iteration
      DO i_iter = 1, n_iter

         bs_env%ri_rs%evgw0_i_iter = i_iter

         ! W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         ! Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
         ! Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
         ! Σ^c_λσ(iτ)  =  Σ_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l')
         CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                              bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

         ! Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
         ! ϵ_n^GW,(i) = ϵ_n^DFT + Σ^c_nn[G^(i-1),W](ϵ_n^GW,(i)) + Σ^x_nn - v^xc_nn
         CALL compute_QP_energies(bs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

         IF (bs_env%gw_flavour == G0W0) EXIT

         CALL print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

         IF (i_iter == n_iter .OR. converged) EXIT

         ! eigenvalues ϵ_n to be updated in G:
         ! G^occ_µλ(i|τ|) = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         ! G^vir_µλ(i|τ|) = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         CALL update_eigenvalues_G(bs_env)

      END DO

      CALL cp_fm_release(fm_W_time)

      CALL reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      CALL delete_unnecessary_files(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_and_QP_energies

! **************************************************************************************************
!> \brief Sets up the evGW0 eigenvalue self-consistency loop: the cycle count, the eigenvalues
!>        the first Green's function is built from, and the DFT reference that the loop
!>        overwrites. A G0W0 run reduces to a single cycle and needs none of it.
!> \param bs_env ...
!> \param n_iter ...
!> \param band_prev ...
!> \param eigenval_scf_Gamma_dft ...
!> \param e_fermi_dft ...
! **************************************************************************************************
   SUBROUTINE init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(OUT)                               :: n_iter
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(OUT)        :: band_prev
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: e_fermi_dft

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_evGW0'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      n_iter = 1
      band_prev(:, :) = 0.0_dp
      e_fermi_dft(:) = 0.0_dp

      IF (bs_env%gw_flavour == evGW0) THEN
         n_iter = bs_env%ri_rs%evgw0_iter

         ! eigenvalues currently in G; the first cycle starts is a plain G0W0 step
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_scf(:, :, :)
         ! the loop overwrites these; post-GW printing expects the DFT values back
         ALLOCATE (eigenval_scf_Gamma_dft, SOURCE=bs_env%eigenval_scf_Gamma)
         e_fermi_dft(:) = bs_env%e_fermi(:)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_evGW0

! **************************************************************************************************
!> \brief Feeds the quasiparticle energies of the current evGW0 cycle back into the Green's
!>        function used by the next one, and re-centres the Fermi level between the new
!>        quasiparticle HOMO and LUMO.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE update_eigenvalues_G(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'update_eigenvalues_G'

      INTEGER                                            :: handle, i_mo, ispin

      CALL timeset(routineN, handle)

      ! record this cycle's evGW0 result; it is also what the next Green's function is built from
      bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW(:, :, :)

      DO ispin = 1, bs_env%n_spin
         ! Update all physical states; exclude linear-dependency placeholders.
         DO i_mo = 1, bs_env%n_mo_retained
            bs_env%eigenval_scf_Gamma(i_mo, ispin) = bs_env%eigenval_GW(i_mo, 1, ispin)
         END DO
         bs_env%e_fermi(ispin) = &
            0.5_dp*(bs_env%eigenval_GW(bs_env%n_occ(ispin), 1, ispin) + &
                    bs_env%eigenval_GW(bs_env%n_occ(ispin) + 1, 1, ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE update_eigenvalues_G

! **************************************************************************************************
!> \brief Restores the DFT reference that the evGW0 loop overwrote and cleanup
!> \param bs_env ...
!> \param eigenval_scf_Gamma_dft ...
!> \param e_fermi_dft ...
!> \param fm_Sigma_x_Gamma ..
! **************************************************************************************************
   SUBROUTINE reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: e_fermi_dft
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reset_and_clean_bs_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (bs_env%gw_flavour == evGW0) THEN
         ! Update array with final evGW0 result
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW(:, :, :)
         ! put the DFT reference back for the post-GW DOS/band-edge printing
         bs_env%eigenval_scf_Gamma(:, :) = eigenval_scf_Gamma_dft(:, :)
         bs_env%e_fermi(:) = e_fermi_dft(:)
         DEALLOCATE (eigenval_scf_Gamma_dft)
         CALL cp_fm_release(fm_Sigma_x_Gamma)
      END IF

      CALL timestop(handle)

   END SUBROUTINE reset_and_clean_bs_env

! **************************************************************************************************
!> \brief Print the quasiparticle HOMO, LUMO and HOMO-LUMO gap of the current evGW0 cycle
!> \param bs_env ...
!> \param band_prev ...
!> \param i_iter ...
!> \param n_iter ...
!> \param converged ...
! **************************************************************************************************
   SUBROUTINE print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(INOUT)      :: band_prev
      INTEGER, INTENT(IN)                                :: i_iter, n_iter
      LOGICAL, INTENT(OUT)                               :: converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_evGW0_band_edges'

      INTEGER                                            :: handle, homo, ispin, u
      REAL(KIND=dp)                                      :: max_delta
      REAL(KIND=dp), DIMENSION(3)                        :: band

      CALL timeset(routineN, handle)

      u = bs_env%unit_nr

      converged = (i_iter > 1)
      max_delta = 0.0_dp

      IF (u > 0) THEN
         WRITE (u, '(A)') ' '
         WRITE (u, '(T2,A)') REPEAT('-', 79)
         WRITE (u, '(T2,A,I4,A,I4)') 'evGW0 cycle', i_iter, ' /', n_iter
         WRITE (u, '(T2,A)') REPEAT('-', 79)
      END IF

      DO ispin = 1, bs_env%n_spin

         homo = bs_env%n_occ(ispin)
         band(1) = bs_env%eigenval_GW(homo, 1, ispin)
         band(2) = bs_env%eigenval_GW(homo + 1, 1, ispin)
         band(3) = band(2) - band(1)

         IF (i_iter > 1) THEN
            max_delta = MAX(max_delta, MAXVAL(ABS(band(:) - band_prev(:, ispin))))
         END IF

         IF (u > 0) THEN
            IF (bs_env%n_spin == 2) WRITE (u, '(T2,A,I0)') 'Spin ', ispin
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 HOMO (eV)', band(1)*evolt
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 LUMO (eV)', band(2)*evolt
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 HOMO-LUMO gap (eV)', band(3)*evolt
         END IF

         band_prev(:, ispin) = band(:)

      END DO

      IF (i_iter > 1) THEN
         converged = (max_delta < bs_env%ri_rs%evgw0_eps_iter)
         IF (u > 0) WRITE (u, '(T2,A,T61,F20.6)') 'Max. change to previous cycle (eV)', &
            max_delta*evolt
      END IF

      IF (u > 0) WRITE (u, '(T2,A)') REPEAT('-', 79)

      IF (converged) THEN
         IF (u > 0) THEN
            WRITE (u, '(A)') ' '
            WRITE (u, '(T2,A,I4,A)') &
               'evGW0 eigenvalue self-consistency reached in', i_iter, ' cycles.'
            WRITE (u, '(A)') ' '
         END IF
      ELSE IF (i_iter == n_iter) THEN
         CALL cp_warn(__LOCATION__, &
                      "The evGW0 eigenvalue self-consistency cycle did not converge "// &
                      "within MAX_ITER cycles. The reported quasiparticle energies are "// &
                      "those of the last cycle.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_evGW0_band_edges

! **************************************************************************************************
!> \brief Compute per-atom AO and RI basis radii from the most diffuse Gaussian
!>        primitive in the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets.
!>        Stores results in bs_env%ri_rs%radius_ao_per_atom(:) and
!>        bs_env%ri_rs%radius_ri_per_atom(:) and prints a per-atom table.
!>        Radius:  r_kind = sqrt(-log(eps) / alpha_min_kind)
!>        with eps = eps_filter.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE precompute_ri_rs_radii(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precompute_ri_rs_radii'

      INTEGER                                            :: handle, i, iatom, ikind, j, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha_min_ao_kind, alpha_min_ri_kind
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet_ao, zet_ri
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set)
      natom = bs_env%n_atom

      eps = bs_env%eps_filter

      ALLOCATE (alpha_min_ao_kind(nkind), alpha_min_ri_kind(nkind))
      alpha_min_ao_kind = HUGE(1.0_dp)
      alpha_min_ri_kind = HUGE(1.0_dp)

      DO ikind = 1, nkind
         zet_ao => bs_env%basis_set_AO(ikind)%gto_basis_set%zet
         zet_ri => bs_env%basis_set_RI(ikind)%gto_basis_set%zet

         DO i = 1, SIZE(zet_ao, 1)
            DO j = 1, SIZE(zet_ao, 2)
               IF (zet_ao(i, j) > 1.0E-3_dp) THEN
                  alpha_min_ao_kind(ikind) = MIN(alpha_min_ao_kind(ikind), zet_ao(i, j))
               END IF
            END DO
         END DO
         DO i = 1, SIZE(zet_ri, 1)
            DO j = 1, SIZE(zet_ri, 2)
               IF (zet_ri(i, j) > 1.0E-3_dp) THEN
                  alpha_min_ri_kind(ikind) = MIN(alpha_min_ri_kind(ikind), zet_ri(i, j))
               END IF
            END DO
         END DO
      END DO

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      ALLOCATE (bs_env%ri_rs%radius_ao_per_atom(natom))
      ALLOCATE (bs_env%ri_rs%radius_ri_per_atom(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         bs_env%ri_rs%radius_ao_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))
         bs_env%ri_rs%radius_ri_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))
      END DO

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS basis radii (Å):'
         WRITE (bs_env%unit_nr, '(T4,A6,2X,A4,2A14)') 'Kind', 'Elem', 'r_AO (Å)', 'r_RI (Å)'
         DO ikind = 1, nkind
            WRITE (bs_env%unit_nr, '(T4,I6,2X,A4,2F14.4)') &
               ikind, &
               atomic_kind_set(ikind)%element_symbol, &
               SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))*angstrom, &
               SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))*angstrom
         END DO
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      DEALLOCATE (alpha_min_ao_kind, alpha_min_ri_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE precompute_ri_rs_radii

! **************************************************************************************************
!> \brief Spreads the low 21 bits of a into every third bit (bits 0,3,6,...,60): the 1-D helper
!>        for a 3-D Morton (Z-order) code. Standard 64-bit magic-mask implementation.
!> \param a value in [0, 2^21)
!> \param x a with two zero bits inserted between consecutive input bits
! **************************************************************************************************
   SUBROUTINE morton_split3(a, x)
      INTEGER(KIND=int_8), INTENT(IN)                    :: a
      INTEGER(KIND=int_8), INTENT(OUT)                   :: x

      x = IAND(a, INT(z'1FFFFF', int_8))
      x = IAND(IOR(x, ISHFT(x, 32)), INT(z'1F00000000FFFF', int_8))
      x = IAND(IOR(x, ISHFT(x, 16)), INT(z'1F0000FF0000FF', int_8))
      x = IAND(IOR(x, ISHFT(x, 8)), INT(z'100F00F00F00F00F', int_8))
      x = IAND(IOR(x, ISHFT(x, 4)), INT(z'10C30C30C30C30C3', int_8))
      x = IAND(IOR(x, ISHFT(x, 2)), INT(z'1249249249249249', int_8))
   END SUBROUTINE morton_split3

! **************************************************************************************************
!> \brief Returns a permutation of atom indices in Morton (Z-order) space-filling order of their
!>        Cartesian centers, so consecutive atoms are spatial neighbors. The RI-RS grid rows are
!>        laid down in this order, so a contiguous grid panel maps to a compact spatial region and
!>        the CUTOFF_RADIUS_RL_W neighborhood of every panel shrinks. The grid row index
!>        is a summed contraction index, so ANY permutation is result-preserving; this one is
!>        chosen purely to improve locality. Coordinates are normalized to the atom bounding box
!>        and quantized to 21 bits per axis (sub-picometre for any real cell).
!> \param particle_set ...
!> \param order order(i) = atom index placed at layout position i
! **************************************************************************************************
   SUBROUTINE spatial_atom_order(particle_set, order)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: order

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spatial_atom_order'
      INTEGER, PARAMETER                                 :: nbits = 21

      INTEGER                                            :: handle, ia, k, natom
      INTEGER(KIND=int_8)                                :: cmax, ic(3), m1, m2, m3
      INTEGER(KIND=int_8), ALLOCATABLE                   :: mcode(:)
      REAL(KIND=dp)                                      :: hi(3), lo(3), span(3)

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)
      ALLOCATE (order(natom), mcode(natom))
      cmax = ISHFT(1_int_8, nbits) - 1_int_8

      lo(:) = HUGE(1.0_dp)
      hi(:) = -HUGE(1.0_dp)
      DO ia = 1, natom
         DO k = 1, 3
            lo(k) = MIN(lo(k), particle_set(ia)%r(k))
            hi(k) = MAX(hi(k), particle_set(ia)%r(k))
         END DO
      END DO
      span(:) = hi(:) - lo(:)
      DO k = 1, 3
         IF (span(k) <= 0.0_dp) span(k) = 1.0_dp
      END DO

      DO ia = 1, natom
         DO k = 1, 3
            ic(k) = INT(((particle_set(ia)%r(k) - lo(k))/span(k))*REAL(cmax, dp), int_8)
            ic(k) = MIN(cmax, MAX(0_int_8, ic(k)))
         END DO
         CALL morton_split3(ic(1), m1)
         CALL morton_split3(ic(2), m2)
         CALL morton_split3(ic(3), m3)
         mcode(ia) = IOR(IOR(m1, ISHFT(m2, 1)), ISHFT(m3, 2))
      END DO

      ! sort(mcode, natom, order): order(i) = original atom index with the i-th smallest code
      CALL sort(mcode, natom, order)

      DEALLOCATE (mcode)

      CALL timestop(handle)

   END SUBROUTINE spatial_atom_order

! **************************************************************************************************
!> \brief Compute grid points for RI-RS
!>        Right now based on Ivan and Xavier implementation
!>        JCP 150, 174120 (2019), JCTC 17, 2383 (2021)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
! **************************************************************************************************
   SUBROUTINE ri_rs_grid_assembler(qs_env, bs_env, ri_rs_grid_points)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ri_rs_grid_points(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_rs_grid_assembler'

      INTEGER                                            :: atom_idx, end_idx, handle, i_layout, &
                                                            natom, start_idx, total_grid_npts
      INTEGER, ALLOCATABLE                               :: atom_order(:), ri_rs_grid_offsets(:)
      REAL(KIND=dp)                                      :: atomic_center(3)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

       !! Get the information about the atoms in the system
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      natom = bs_env%n_atom

       !! 1. Build the grid cache
      CALL build_grid_cache(bs_env, atomic_kind_set)

      !! 2. Optimize RIRS grid points
      IF (bs_env%ri_rs%grid_opt%enabled) THEN
         CALL optimize_ri_rs_grid(qs_env, bs_env)
      END IF

      !! 3. Calculate grid counts and offsets
      ALLOCATE (ri_rs_grid_offsets(natom + 1))

      ! grid_atom_boundaries(k) = 1-based start of the k-th atom's grid run in LAYOUT order
      ! (the order in which the points are laid down below). Used later to build atom-aligned
      ! DBCSR grid row-blocks in atomic_basis_at_grid_point.
      IF (ALLOCATED(bs_env%ri_rs%grid_atom_boundaries)) DEALLOCATE (bs_env%ri_rs%grid_atom_boundaries)
      ALLOCATE (bs_env%ri_rs%grid_atom_boundaries(natom + 1))

      ! Lay the grid rows down in Morton (spatial) atom order, so a
      ! contiguous grid maps to a compact spatial region.
      CALL spatial_atom_order(particle_set, atom_order)

      total_grid_npts = 0
      DO i_layout = 1, natom
         atom_idx = atom_order(i_layout)
         ri_rs_grid_offsets(atom_idx) = total_grid_npts + 1
         bs_env%ri_rs%grid_atom_boundaries(i_layout) = total_grid_npts + 1
         total_grid_npts = total_grid_npts + bs_env%ri_rs%grid_cache(atom_idx)%npts
      END DO

      ri_rs_grid_offsets(natom + 1) = total_grid_npts + 1
      bs_env%ri_rs%grid_atom_boundaries(natom + 1) = total_grid_npts + 1

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, FMT="(T2,A,T69,I12)") &
            'Total grid points used for RI-RS:', total_grid_npts
         WRITE (bs_env%unit_nr, "(A)") ' '
      END IF

       !! 4. Allocate the global ri_rs_grid arrays
      ALLOCATE (ri_rs_grid_points(3, total_grid_npts))

       !! 5. Parallelize the grid generation loop
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(ri_rs_grid_points, ri_rs_grid_offsets, &
      !$OMP        particle_set, bs_env, natom) &
      !$OMP PRIVATE(atom_idx, atomic_center, start_idx, end_idx) &
      !$OMP SCHEDULE(DYNAMIC, 1)
      DO atom_idx = 1, natom
         atomic_center(:) = particle_set(atom_idx)%r(:)

         start_idx = ri_rs_grid_offsets(atom_idx)
         end_idx = start_idx + bs_env%ri_rs%grid_cache(atom_idx)%npts - 1

         !! Shift the cached origin grid by the atom's center
         ri_rs_grid_points(1, start_idx:end_idx) = bs_env%ri_rs%grid_cache(atom_idx)%raw_points(1, :) + atomic_center(1)
         ri_rs_grid_points(2, start_idx:end_idx) = bs_env%ri_rs%grid_cache(atom_idx)%raw_points(2, :) + atomic_center(2)
         ri_rs_grid_points(3, start_idx:end_idx) = bs_env%ri_rs%grid_cache(atom_idx)%raw_points(3, :) + atomic_center(3)

      END DO
      !$OMP END PARALLEL DO

      !! 6. Cleanup memory
      IF (ALLOCATED(bs_env%ri_rs%grid_cache)) THEN
         DO atom_idx = 1, natom
            IF (ALLOCATED(bs_env%ri_rs%grid_cache(atom_idx)%raw_points)) THEN
               DEALLOCATE (bs_env%ri_rs%grid_cache(atom_idx)%raw_points)
            END IF
         END DO
         DEALLOCATE (bs_env%ri_rs%grid_cache)
      END IF

      DEALLOCATE (atom_order, ri_rs_grid_offsets)

      CALL timestop(handle)

   END SUBROUTINE ri_rs_grid_assembler

! **************************************************************************************************
!> \brief Reads grids from .ion files and stores them in memory based on grid_select
!> \param bs_env ...
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE build_grid_cache(bs_env, atomic_kind_set)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_grid_cache'

      CHARACTER(LEN=default_path_length)                 :: full_path, warning_msg
      CHARACTER(LEN=default_string_length)               :: atom_sym
      INTEGER                                            :: handle, iatom, ikind, &
                                                            n_voronoi_candidates, natom, nkind, &
                                                            requested_npts
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(rirs_grid_type), ALLOCATABLE                  :: kind_grid(:)

      CALL timeset(routineN, handle)

      nkind = SIZE(atomic_kind_set)
      particle_set => bs_env%ri_rs%particle_set
      CPASSERT(ASSOCIATED(particle_set))
      natom = bs_env%n_atom
      ALLOCATE (kind_grid(nkind))

      DO ikind = 1, nkind
         atom_sym = TRIM(atomic_kind_set(ikind)%element_symbol)
         CALL get_rirs_grid_filepath(atom_sym, bs_env%ri_rs%grid_select, &
                                     bs_env%ri_rs%grid_file_suffix, full_path)
         CALL read_rirs_grid_file(TRIM(full_path), kind_grid(ikind)%raw_points)
         kind_grid(ikind)%npts = SIZE(kind_grid(ikind)%raw_points, 2)
      END DO

      ! Make every atom own an independent copy.  The optimizer can therefore adapt symmetry-
      ! equivalent atoms to their individual local environments, and GW consumes those copies
      ! directly without writing intermediate .ion files.
      IF (ALLOCATED(bs_env%ri_rs%grid_cache)) DEALLOCATE (bs_env%ri_rs%grid_cache)
      ALLOCATE (bs_env%ri_rs%grid_cache(natom))
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         bs_env%ri_rs%grid_cache(iatom) = kind_grid(ikind)
         IF (bs_env%ri_rs%grid_opt%enabled) THEN
            atom_sym = TRIM(atomic_kind_set(ikind)%element_symbol)
            requested_npts = bs_env%ri_rs%grid_opt%grid_npoints(ikind)
            CALL initialize_rirs_grid(bs_env%ri_rs%grid_cache(iatom)%raw_points, &
                                      requested_npts, iatom, particle_set, &
                                      n_voronoi_candidates)
            IF (n_voronoi_candidates < requested_npts) THEN
               WRITE (warning_msg, '(A,I0,A,A,A,I0,A,I0,A)') &
                  "Atom ", iatom, " (", TRIM(atom_sym), "): Grid selection requests ", requested_npts, &
                  " points, but only ", n_voronoi_candidates, &
                  " are inside its Voronoi volume. Selecting from the full RI-RS source grid."
               CALL cp_warn(__LOCATION__, TRIM(warning_msg))
            END IF
            bs_env%ri_rs%grid_cache(iatom)%npts = requested_npts
         END IF
      END DO
      DEALLOCATE (kind_grid)

      CALL timestop(handle)

   END SUBROUTINE build_grid_cache

! **************************************************************************************************
!> \brief Evaluates the AO basis on the RI-RS grid and stores it as the sparse DBCSR matrix
!>        Φ_μl = Φ_μ(r_l)  (rows = grid points in atom-aligned blocks of at most
!>        max_elements_per_block points, columns = one block per atom's full AO set).
!>        Grid points outside the reach of an atom's most
!>        diffuse Gaussian (or the CUTOFF_RADIUS_RL_AO) are skipped, and only blocks
!>        with at least one element > eps_filter are stored. This locality is the source of
!>        ALL grid-dimension sparsity used downstream. Also caches the atom centers and the
!>        per-chunk centroids needed by the optional CUTOFF_RADIUS_G_W / CUTOFF_RADIUS_RL_W
!>        operator truncations.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
! **************************************************************************************************
   SUBROUTINE atomic_basis_at_grid_point(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_phi_mu_l

      CHARACTER(LEN=*), PARAMETER :: routineN = 'atomic_basis_at_grid_point'

      INTEGER :: bs_eff, c_size, dimen_ORB, handle, i, i_blk, ia, iatom, natom, npcol, nprow, &
         num_grid_chunks, r_end, r_start, remaining, run, safe_max, total_grid_npts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_row_start, first_sgf
      INTEGER, DIMENSION(:), POINTER                     :: c_blk_sizes, col_dist, r_blk_sizes, &
                                                            row_dist
      REAL(KIND=dp)                                      :: r2_threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_col_buffer
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_ks
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! Extract environment variables
      CALL get_qs_env(qs_env, cell=cell, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      para_env=para_env)

      natom = bs_env%n_atom
      total_grid_npts = SIZE(ri_rs_grid_points, 2)

      ! Map the starting indices of spherical gaussian functions (SGF) for each atom
      ALLOCATE (first_sgf(natom + 1))
      CALL get_basis_offsets(particle_set, qs_kind_set, first_sgf, dimen_ORB)

      ! =========================================================================
      ! 1. SETUP DBCSR MATRIX TOPOLOGY
      ! =========================================================================

      ! A. Define Column Block Sizes (1 Block = 1 Atom's full basis set)
      ALLOCATE (c_blk_sizes(natom))
      DO iatom = 1, natom
         c_blk_sizes(iatom) = first_sgf(iatom + 1) - first_sgf(iatom)
      END DO

      ! B. Define Row Block Sizes: atom-aligned blocks (a block never spans two atoms' grid
      ! runs), each atom's run subdivided into blocks of at most bs_eff points.

      ! Fetch CP2K's default process grid configuration
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist_ks)
      CALL dbcsr_distribution_get(dbcsr_dist_ks, nprows=nprow, npcols=npcol)

      ! Overflow-safe upper bound on the block size (see bs_env%dbcsr_msg_elem_limit).
      safe_max = INT(0.5_dp*REAL(bs_env%dbcsr_msg_elem_limit, dp)* &
                     REAL(MAX(MIN(nprow, npcol), 1), dp)/ &
                     REAL(total_grid_npts, dp))
      safe_max = MAX(1, safe_max)
      ! Block size = CP2K's global max_elements_per_block (GLOBAL/DBCSR input; default 32),
      ! overflow-capped.
      bs_eff = MAX(1, MIN(max_elements_per_block, safe_max))

      ! Count the atom-aligned blocks, then fill r_blk_sizes and each block's starting grid row.
      num_grid_chunks = 0
      DO ia = 1, natom
         run = bs_env%ri_rs%grid_atom_boundaries(ia + 1) - bs_env%ri_rs%grid_atom_boundaries(ia)
         IF (run > 0) num_grid_chunks = num_grid_chunks + (run + bs_eff - 1)/bs_eff
      END DO
      ALLOCATE (r_blk_sizes(num_grid_chunks), blk_row_start(num_grid_chunks))
      i_blk = 0
      r_start = 1
      DO ia = 1, natom
         remaining = bs_env%ri_rs%grid_atom_boundaries(ia + 1) - bs_env%ri_rs%grid_atom_boundaries(ia)
         DO WHILE (remaining > 0)
            i_blk = i_blk + 1
            r_blk_sizes(i_blk) = MIN(bs_eff, remaining)
            blk_row_start(i_blk) = r_start
            r_start = r_start + r_blk_sizes(i_blk)
            remaining = remaining - r_blk_sizes(i_blk)
         END DO
      END DO

      IF (bs_env%unit_nr > 0) THEN
         ! T71 compensates for the two two-byte Greek characters in the label.
         WRITE (bs_env%unit_nr, '(T2,A,T71,I12)') &
            'RI-RS grid row-blocks of Φ_μ(r_l)', num_grid_chunks
         WRITE (bs_env%unit_nr, '(T2,A,T69,I12)') 'RI-RS grid points per block (max)', bs_eff
      END IF

      ! Cache atomic positions: AO and RI blocks are one-block-per-atom, so these are the
      ! block centers used by the optional CUTOFF_RADIUS_G_W atom-pair truncation.
      IF (ALLOCATED(bs_env%ri_rs%atom_centers)) DEALLOCATE (bs_env%ri_rs%atom_centers)
      ALLOCATE (bs_env%ri_rs%atom_centers(3, natom))
      DO iatom = 1, natom
         bs_env%ri_rs%atom_centers(1:3, iatom) = particle_set(iatom)%r(1:3)
      END DO

      ! Cache per-chunk centroids for the optional CUTOFF_RADIUS_RL_W / CUTOFF_RADIUS_W0 block
      ! truncations. Left unallocated otherwise, so PRESENT(centroids) stays .FALSE. at the
      ! contract_grid_panels* calls.
      IF (bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .OR. &
          bs_env%ri_rs%cutoff_radius_w0 > 0.0_dp) THEN
         IF (ALLOCATED(bs_env%ri_rs%chunk_centroids)) DEALLOCATE (bs_env%ri_rs%chunk_centroids)
         ALLOCATE (bs_env%ri_rs%chunk_centroids(3, num_grid_chunks))
         DO i_blk = 1, num_grid_chunks
            r_start = blk_row_start(i_blk)
            r_end = r_start + r_blk_sizes(i_blk) - 1
            bs_env%ri_rs%chunk_centroids(1, i_blk) = &
               SUM(ri_rs_grid_points(1, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
            bs_env%ri_rs%chunk_centroids(2, i_blk) = &
               SUM(ri_rs_grid_points(2, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
            bs_env%ri_rs%chunk_centroids(3, i_blk) = &
               SUM(ri_rs_grid_points(3, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
         END DO
      END IF

      ! C. Build Custom Mappings using Round-Robin across the 2D process grid

      ALLOCATE (row_dist(num_grid_chunks))
      DO i = 1, num_grid_chunks
         row_dist(i) = MOD(i - 1, nprow)
      END DO

      ALLOCATE (col_dist(natom))
      DO i = 1, natom
         col_dist(i) = MOD(i - 1, npcol)
      END DO

      ! E. Create the DBCSR Distribution and Initialize the Matrix
      CALL dbcsr_distribution_new(dist, template=dbcsr_dist_ks, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(mat_phi_mu_l, name="phi_val_sparse", dist=dist, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=r_blk_sizes, col_blk_size=c_blk_sizes)

      ! =========================================================================
      ! 2. STREAM DATA DIRECTLY INTO SPARSE MATRIX
      ! =========================================================================
      ! Iterate over the atoms assigned to this specific MPI rank
      DO iatom = para_env%mepos + 1, natom, para_env%num_pe

         c_size = c_blk_sizes(iatom)

         ! Allocate a temporary dense buffer just for this specific atom
         ALLOCATE (atom_col_buffer(total_grid_npts, c_size))
         atom_col_buffer = 0.0_dp

         ! Evaluate the basis functions on the grid. Skip grid points outside
         ! the spatial extent of the most diffuse AO Gaussian on iatom; beyond
         ! that radius the contribution is guaranteed below eps_filter. A positive
         ! CUTOFF_RADIUS_RL_AO overrides this with a user-defined hard cutoff.
         IF (bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) THEN
            r2_threshold = bs_env%ri_rs%cutoff_radius_ri_ao**2
         ELSE
            r2_threshold = bs_env%ri_rs%radius_ao_per_atom(iatom)**2
         END IF
         CALL fill_phi_for_atom(atom_col_buffer, ri_rs_grid_points, total_grid_npts, &
                                iatom, particle_set, qs_kind_set, cell, r2_threshold)

         ! Slice the dense column into the atom-aligned grid row-blocks and insert into DBCSR
         DO i_blk = 1, num_grid_chunks
            r_start = blk_row_start(i_blk)
            r_end = r_start + r_blk_sizes(i_blk) - 1

            ! Apply dynamic sparsity filtering: Only store blocks with physical significance
            IF (MAXVAL(ABS(atom_col_buffer(r_start:r_end, 1:c_size))) > bs_env%eps_filter) THEN
               CALL dbcsr_put_block(mat_phi_mu_l, row=i_blk, col=iatom, &
                                    block=atom_col_buffer(r_start:r_end, 1:c_size))
            END IF
         END DO

         DEALLOCATE (atom_col_buffer)

      END DO

      CALL dbcsr_finalize(mat_phi_mu_l)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
      CALL print_matrix_occupation(mat_phi_mu_l, 'Φ_μ(r_l)', bs_env)

      ! -------------------------------------------------------------------------
      ! CLEANUP
      ! -------------------------------------------------------------------------
      DEALLOCATE (first_sgf, r_blk_sizes, c_blk_sizes, row_dist, col_dist, blk_row_start)
      CALL dbcsr_distribution_release(dist)

      CALL timestop(handle)

   END SUBROUTINE atomic_basis_at_grid_point

! **************************************************************************************************
!> \brief Computes the AO basis offsets: first_sgf(iatom) is the global index of the first
!>        spherical Gaussian function (SGF) of iatom, first_sgf(natom+1) = total_sgf + 1,
!>        and total_sgf is the total number of AO basis functions.
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param first_sgf ...
!> \param total_sgf ...
! **************************************************************************************************
   SUBROUTINE get_basis_offsets(particle_set, qs_kind_set, first_sgf, total_sgf)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(OUT)                               :: first_sgf(:), total_sgf

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_basis_offsets'

      INTEGER                                            :: handle, iatom, ikind, nsgf

      CALL timeset(routineN, handle)

      total_sgf = 0
      DO iatom = 1, SIZE(particle_set)
         first_sgf(iatom) = total_sgf + 1
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, basis_type="ORB")
         total_sgf = total_sgf + nsgf
      END DO
      first_sgf(SIZE(particle_set) + 1) = total_sgf + 1

      CALL timestop(handle)

   END SUBROUTINE get_basis_offsets

! **************************************************************************************************
!> \brief Computes χ_PQ(iτ) from the occupied and virtual Green's functions.
!> \param bs_env GW calculation environment and imaginary-time grid
!> \param mat_chi_Gamma_tau Response matrices χ_PQ(iτ) to be computed
!> \param mat_phi_mu_l AO values Φ_μ(r_l) on the RI-RS grid
!> \param mat_Z_lP RI-RS fitting matrix Z_lP
! **************************************************************************************************
   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, mat_phi_mu_l, mat_Z_lP)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i_t, ispin, n_panels
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      REAL(KIND=dp)                                      :: grid_occ, t1, tau
      TYPE(dbcsr_type)                                   :: matrix_G_occ_ao, matrix_G_vir_ao

      CALL timeset(routineN, handle)

      ! Panel boundaries for the grid-streaming contraction.
      ! The panels are identical for χ, Σ^x and Σ^c, so the count is reported once here for all three stages.

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

      n_panels = SIZE(pan_first)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T74,I9)') &
            'Number of batches for χ, Σ matrices', n_panels
         WRITE (bs_env%unit_nr, '(A)') ' '
         CALL m_flush(bs_env%unit_nr)
      END IF

      ! =========================================================================
      ! IMAGINARY TIME LOOP
      ! χ_PQ(iτ) = Σ_s g_s · Z^T ( (φ G^occ_s φ^T) ∘ (φ G^vir_s φ^T) ) Z
      ! (g_s = spin degeneracy)
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         t1 = m_walltime()
         tau = bs_env%imag_time_points(i_t)

         DO ispin = 1, bs_env%n_spin

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_ao)
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_ao)

            ! χ_PQ += g_s · Z^T ( (φ G^occ φ^T) ∘ (φ G^vir φ^T) ) Z
            CALL contract_grid_panels(L_A=mat_phi_mu_l, M_A=matrix_G_occ_ao, &
                                      L_B=mat_phi_mu_l, M_B=matrix_G_vir_ao, &
                                      L_out=mat_Z_lP, mat_out=mat_chi_Gamma_tau(i_t)%matrix, &
                                      scale=bs_env%spin_degeneracy, eps=bs_env%eps_filter, &
                                      para_env=bs_env%para_env, &
                                      pan_first=pan_first, pan_last=pan_last, &
                                      lb_eq_la=.TRUE., lout_eq_la=.FALSE., &
                                      zero_out=(ispin == 1), &
                                      keep_sparsity=bs_env%ri_rs%keep_sparsity_rirs, &
                                      centroids=bs_env%ri_rs%chunk_centroids, &
                                      cutoff=bs_env%ri_rs%cutoff_radius_v_w, &
                                      grid_occupation=grid_occ)

            CALL dbcsr_release(matrix_G_occ_ao)
            CALL dbcsr_release(matrix_G_vir_ao)

         END DO ! ispin

         ! Sparsity reports
         IF (i_t == 1) THEN
            CALL print_matrix_occupation(mat_Z_lP, 'Z_lP', bs_env)
            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T73,F7.2,A)') &
                  'Percentage of non-zero matrix elements in G_ll'', χ_ll'', W_ll''', &
                  grid_occ*100.0_dp, ' %'
               CALL m_flush(bs_env%unit_nr)
            END IF
            CALL print_matrix_occupation(mat_chi_Gamma_tau(i_t)%matrix, 'χ_PQ', bs_env, &
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
!> \brief Marks the grid blocks whose centroid lies within cutoff of the bounding box of the
!>        panel [blk0, blk1]'s chunk centroids.
!> \param centroids ...
!> \param blk0 ...
!> \param blk1 ...
!> \param cutoff ...
!> \param used ...
! **************************************************************************************************
   SUBROUTINE mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      INTEGER, INTENT(IN)                                :: blk0, blk1
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mask_grid_blocks_near_panel'

      INTEGER                                            :: c, handle, k
      REAL(KIND=dp)                                      :: cutoff2, d2, dx
      REAL(KIND=dp), DIMENSION(3)                        :: hi, lo

      CALL timeset(routineN, handle)

      cutoff2 = cutoff**2
      lo(:) = MINVAL(centroids(:, blk0:blk1), DIM=2)
      hi(:) = MAXVAL(centroids(:, blk0:blk1), DIM=2)

      ALLOCATE (used(SIZE(centroids, 2)))
      DO c = 1, SIZE(centroids, 2)
         d2 = 0.0_dp
         DO k = 1, 3
            dx = MAX(0.0_dp, lo(k) - centroids(k, c), centroids(k, c) - hi(k))
            d2 = d2 + dx*dx
         END DO
         used(c) = (d2 <= cutoff2)
      END DO

      CALL timestop(handle)

   END SUBROUTINE mask_grid_blocks_near_panel

! **************************************************************************************************
!> \brief Exact allocated-element count of the geo template of panel [blk0, blk1]: the very same
!>        per-block-pair centroid test as build_geo_template_panel, so this is the true DBCSR
!>        data size of A_pan/B_pan/C_pan (DBCSR stores whole blocks).
!> \param r_blk_sizes ...
!> \param centroids ...
!> \param used ...
!> \param blk0 ...
!> \param blk1 ...
!> \param cutoff ...
!> \param nze_tmpl ...
! **************************************************************************************************
   SUBROUTINE panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: used
      INTEGER, INTENT(IN)                                :: blk0, blk1
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER(KIND=int_8), INTENT(OUT)                   :: nze_tmpl

      CHARACTER(LEN=*), PARAMETER :: routineN = 'panel_template_elems'

      INTEGER                                            :: c, handle, ib, n_used
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: used_idx
      REAL(KIND=dp)                                      :: cutoff2

      CALL timeset(routineN, handle)

      ! Compress the near mask once so the pair loop only visits candidate columns.
      n_used = COUNT(used)
      ALLOCATE (used_idx(n_used))
      n_used = 0
      DO c = 1, SIZE(used)
         IF (used(c)) THEN
            n_used = n_used + 1
            used_idx(n_used) = c
         END IF
      END DO

      cutoff2 = cutoff**2
      nze_tmpl = 0_int_8
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(blk0, blk1, n_used, used_idx, centroids, cutoff2, &
      !$OMP r_blk_sizes) PRIVATE(ib, c) REDUCTION(+:nze_tmpl)
      DO ib = blk0, blk1
         DO c = 1, n_used
            IF (SUM((centroids(:, ib) - centroids(:, used_idx(c)))**2) <= cutoff2) THEN
               nze_tmpl = nze_tmpl + INT(r_blk_sizes(ib), int_8)*INT(r_blk_sizes(used_idx(c)), int_8)
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE panel_template_elems

! **************************************************************************************************
!> \brief Per-rank peak memory (GB) of one panel step of the neighborhood-restricted
!>        contractions: three grid x grid panels of the template size (A_pan, B_pan, C_pan)
!>        plus the grid x RI intermediates (tmp2 and the accumulation operand) and the
!>        grid x AO intermediate (tmpA), whose column support is the panel's geometric
!>        neighborhood fraction f_near = width/n_grid. Shared by the panel planner and
!> \param nze_tmpl ...
!> \param pan_rows ...
!> \param width ...
!> \param n_grid_total ...
!> \param n_RI ...
!> \param n_ao ...
!> \param n_procs ...
!> \param mem_GB ...
! **************************************************************************************************
   SUBROUTINE panel_mem_estimate_GB(nze_tmpl, pan_rows, width, n_grid_total, n_RI, n_ao, &
                                    n_procs, mem_GB)

      INTEGER(KIND=int_8), INTENT(IN)                    :: nze_tmpl
      INTEGER, INTENT(IN)                                :: pan_rows, width, n_grid_total, n_RI, &
                                                            n_ao, n_procs
      REAL(KIND=dp), INTENT(OUT)                         :: mem_GB

      REAL(KIND=dp)                                      :: f_near

      f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
      mem_GB = (3.0_dp*REAL(nze_tmpl, dp) + &
                REAL(pan_rows, dp)*f_near*(2.0_dp*REAL(n_RI, dp) + REAL(n_ao, dp)))* &
               8.0_dp/REAL(MAX(n_procs, 1), dp)*1.0E-9_dp

   END SUBROUTINE panel_mem_estimate_GB

! **************************************************************************************************
!> \brief Plans the panel boundaries for the streaming contractions. Panels grow by whole grid
!>        row-blocks towards ~panel_size rows. When the neighborhood restriction is active
!>        (centroids+cutoff), each candidate panel is additionally checked against
!>          (a) the 32-bit message bound with the panel's TRUE occupancy
!>          (b) the per-rank memory budget: panel_mem_estimate_GB <= mem_budget_GB.
!> \param bs_env ...
!> \param r_blk_sizes ...
!> \param panel_size ...
!> \param min_dim ...
!> \param pan_first ...
!> \param pan_last ...
!> \param centroids ...
!> \param cutoff ...
!> \param n_RI ...
!> \param n_ao ...
!> \param n_procs ...
!> \param mem_budget_GB ...
!> \param honor_exact ...
!> \param unsafe ...
! **************************************************************************************************
   SUBROUTINE plan_grid_panels(bs_env, r_blk_sizes, panel_size, min_dim, pan_first, pan_last, &
                               centroids, cutoff, n_RI, n_ao, n_procs, mem_budget_GB, &
                               honor_exact, unsafe)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      INTEGER, INTENT(IN)                                :: panel_size, min_dim
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      INTEGER, INTENT(IN), OPTIONAL                      :: n_RI, n_ao, n_procs
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: mem_budget_GB
      LOGICAL, INTENT(IN), OPTIONAL                      :: honor_exact
      LOGICAL, INTENT(OUT), OPTIONAL                     :: unsafe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'plan_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, ib, n_grid_blocks, &
                                                            n_grid_total, n_panels, rows_acc, &
                                                            TARGET, width
      INTEGER(KIND=int_8)                                :: msg, nze_tmpl, side
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp_first, tmp_last
      LOGICAL                                            :: fits, my_honor_exact, my_unsafe, &
                                                            use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: used
      REAL(KIND=dp)                                      :: f_near, mem_GB

      CALL timeset(routineN, handle)

      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp
      IF (use_cutoff) THEN
         CPASSERT(PRESENT(n_RI) .AND. PRESENT(n_ao) .AND. PRESENT(n_procs))
      END IF

      ! honor_exact: use N_PANELS as requested -- do NOT split a panel further even if it trips
      ! the message-overflow / memory-budget check; instead flag `unsafe` so the caller can warn.
      my_honor_exact = .FALSE.
      IF (PRESENT(honor_exact)) my_honor_exact = honor_exact
      my_unsafe = .FALSE.

      n_grid_blocks = SIZE(r_blk_sizes)
      n_grid_total = SUM(r_blk_sizes)
      ALLOCATE (tmp_first(n_grid_blocks), tmp_last(n_grid_blocks))

      n_panels = 0
      blk0 = 1
      DO WHILE (blk0 <= n_grid_blocks)
         TARGET = panel_size
         DO
            rows_acc = 0
            blk1 = blk0
            DO ib = blk0, n_grid_blocks
               rows_acc = rows_acc + r_blk_sizes(ib)
               blk1 = ib
               IF (rows_acc >= TARGET) EXIT
            END DO
            IF (.NOT. use_cutoff .OR. blk1 == blk0) EXIT
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)
            width = SUM(r_blk_sizes, MASK=used)
            CALL panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)
            f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
            side = INT(REAL(rows_acc, dp)*f_near*REAL(MAX(n_RI, n_ao), dp), int_8)
            msg = MAX(nze_tmpl, side)/INT(MAX(min_dim, 1), int_8)
            fits = (msg <= bs_env%dbcsr_msg_elem_limit/4)
            IF (fits .AND. PRESENT(mem_budget_GB)) THEN
               IF (mem_budget_GB > 0.0_dp) THEN
                  CALL panel_mem_estimate_GB(nze_tmpl, rows_acc, width, n_grid_total, &
                                             n_RI, n_ao, n_procs, mem_GB)
                  fits = (mem_GB <= mem_budget_GB)
               END IF
            END IF
            ! Panel size is bounded only by the message-overflow and memory checks above; there is
            ! no neighborhood-width (f_near) cap. mp_waitall is dominated by the NUMBER of panel
            ! multiplies, so fewer/larger panels are cheaper here -- panel count is driven DOWN by
            ! the N_PANELS keyword (panel_size), not split up by a width heuristic.
            IF (my_honor_exact) THEN
               ! Keep exactly the requested grouping; just record if it exceeds a safety limit.
               IF (.NOT. fits) my_unsafe = .TRUE.
               EXIT
            END IF
            IF (fits) EXIT
            TARGET = MAX(1, MIN(TARGET, rows_acc)/2)
         END DO
         n_panels = n_panels + 1
         tmp_first(n_panels) = blk0
         tmp_last(n_panels) = blk1
         blk0 = blk1 + 1
      END DO

      ALLOCATE (pan_first(n_panels), pan_last(n_panels))
      pan_first(:) = tmp_first(1:n_panels)
      pan_last(:) = tmp_last(1:n_panels)
      DEALLOCATE (tmp_first, tmp_last)

      IF (PRESENT(unsafe)) unsafe = my_unsafe

      CALL timestop(handle)

   END SUBROUTINE plan_grid_panels

! **************************************************************************************************
!> \brief Resolves the panel boundaries for the streaming contractions from the bs_env settings:
!> \param bs_env ...
!> \param mat_phi_mu_l ...
!> \param pan_first ...
!> \param pan_last ...
! **************************************************************************************************
   SUBROUTINE resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last

      CHARACTER(LEN=*), PARAMETER :: routineN = 'resolve_grid_panels'

      CHARACTER(LEN=max_line_length)                     :: msg
      INTEGER                                            :: handle, min_dim, n_grid_total, &
                                                            n_panels_req, npcols, nprows, &
                                                            panel_size, safe_max
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: honor_exact, panels_unsafe, use_cutoff
      REAL(KIND=dp)                                      :: mem_avail_GB, mem_budget_GB
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      IF (ALLOCATED(bs_env%ri_rs%pan_first)) THEN
         ALLOCATE (pan_first, SOURCE=bs_env%ri_rs%pan_first)
         ALLOCATE (pan_last, SOURCE=bs_env%ri_rs%pan_last)
         CALL timestop(handle)
         RETURN
      END IF

      use_cutoff = bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .AND. &
                   ALLOCATED(bs_env%ri_rs%chunk_centroids)

      ! MIN(nprows, npcols) is the divisor that bounds the worst-rank Cannon message: a
      ! P x n_grid panel is replicated into block row strips (P/nprows x n_grid) or column
      ! strips (P x n_grid/npcols) during multiply_cannon, so the largest single-rank
      ! message is ~ P*n_grid / MIN(nprows,npcols) elements.
      CALL dbcsr_get_info(mat_phi_mu_l, nfullrows_total=n_grid_total, row_blk_size=r_blk_sizes, &
                          distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, npcols=npcols)
      min_dim = MAX(MIN(nprows, npcols), 1)

      ! Panel height such that NO per-rank DBCSR message can overflow the 32-bit length field
      ! (see bs_env%dbcsr_msg_elem_limit): requiring the worst-rank message to stay under
      ! 0.5 * HUGE(int_4) gives the safe height P_safe = 0.5 * HUGE(int_4) * min_dim / n_grid.
      IF (use_cutoff) THEN
         safe_max = n_grid_total
      ELSE
         safe_max = INT(0.5_dp*REAL(bs_env%dbcsr_msg_elem_limit, dp)*REAL(min_dim, dp)/ &
                        REAL(n_grid_total, dp))
         safe_max = MAX(1, MIN(safe_max, n_grid_total))
      END IF

      ! A user-set N_PANELS ( > 1 ) is honored EXACTLY: the planner produces that many panels
      ! (up to grid-block granularity) and never force-splits them for the message/memory safety
      ! limits -- if a limit is tripped it warns instead of silently changing the count.
      n_panels_req = bs_env%ri_rs%n_panels
      honor_exact = (n_panels_req > 1)
      panels_unsafe = .FALSE.
      IF (n_panels_req > 1) THEN
         ! ceil(n_grid/n_panels_req) rows per panel => exactly n_panels_req panels. With the
         ! cutoff active safe_max = n_grid_total (no clamp, honored exactly); without it, safe_max
         ! is the int32-overflow ceiling and MUST still bound the panel (the non-cutoff planner
         ! loop has no in-loop message-size check).
         panel_size = MIN((n_grid_total + n_panels_req - 1)/n_panels_req, safe_max)
      ELSE
         ! Default (<= 1): a single whole-grid panel, clamped to the overflow-safe ceiling.
         panel_size = safe_max
      END IF
      panel_size = MAX(1, panel_size)

      IF (use_cutoff) THEN
         ! Half of the measured free memory as panel budget.
         CALL ri_rs_mem_avail_per_proc_GB(bs_env, mem_avail_GB)
         mem_budget_GB = 0.5_dp*mem_avail_GB
         CALL plan_grid_panels(bs_env, r_blk_sizes, panel_size, min_dim, pan_first, pan_last, &
                               centroids=bs_env%ri_rs%chunk_centroids, &
                               cutoff=bs_env%ri_rs%cutoff_radius_v_w, &
                               n_RI=bs_env%n_RI, n_ao=bs_env%n_ao, &
                               n_procs=bs_env%para_env%num_pe, mem_budget_GB=mem_budget_GB, &
                               honor_exact=honor_exact, unsafe=panels_unsafe)
      ELSE
         CALL plan_grid_panels(bs_env, r_blk_sizes, panel_size, min_dim, pan_first, pan_last)
      END IF

      IF (honor_exact .AND. panels_unsafe) THEN
         WRITE (msg, '(A,I0,A)') &
            "N_PANELS = ", n_panels_req, " is used as requested, but one or more panels "// &
            "exceed the DBCSR 32-bit message length or the memory budget. The run may abort "// &
            "or swap; increase N_PANELS if it does."
         CPWARN(TRIM(msg))
      END IF

      ALLOCATE (bs_env%ri_rs%pan_first, SOURCE=pan_first)
      ALLOCATE (bs_env%ri_rs%pan_last, SOURCE=pan_last)

      CALL timestop(handle)

   END SUBROUTINE resolve_grid_panels

! **************************************************************************************************
!> \brief Estimates and prints per-process memory requirements for the RI-RS GW calculation.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_ri_rs_memory_estimate(qs_env, bs_env)

!$    USE OMP_LIB, ONLY: omp_get_max_threads

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_ri_rs_memory_estimate'

      CHARACTER(LEN=max_line_length)                     :: msg
      INTEGER                                            :: handle, iatom, ipan, l, &
                                                            max_n_ao_used, max_n_local_grid, &
                                                            n_ao_used_atom, n_grid_total, &
                                                            n_local_grid, n_loc_ri_max, n_procs, &
                                                            n_procs_per_atom, n_RI, n_threads, &
                                                            natom, pan_rows, pan_width
      INTEGER(KIND=int_8)                                :: nze_tmpl
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used
      REAL(KIND=dp)                                      :: cutoff_ri, mem_avail_GB, mem_D_local_GB, &
                                                            mem_dlp_GB, mem_pan_GB, mem_panels_GB, &
                                                            mem_phi_local_GB, mem_Z_lP_GB, &
                                                            mem_Zlp_peak_GB, pos_P(3)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, nfullrows_total=n_grid_total, &
                          row_blk_size=r_blk_sizes)
      n_RI = bs_env%n_RI
      n_procs = bs_env%para_env%num_pe

      ! Z_lP upper bound: dense n_grid × n_RI, distributed evenly across all ranks.
      ! The actual sparse Z_lP is smaller due to the per-atom locality cutoff.
      mem_Z_lP_GB = REAL(n_grid_total, dp)*REAL(n_RI, dp)*8.0_dp/ &
                    REAL(n_procs, dp)*1.0E-9_dp

      ! Peak panel memory during Σ^c: two G panels (A_occ, A_vir) + one W panel plus the
      ! grid × RI / grid × AO intermediates. With the CUTOFF_RADIUS_RL_W restriction the panel
      ! matrices only allocate the geo-template blocks, so use the same nze-aware model as the
      ! panel planner (panel_mem_estimate_GB); without the cutoff, dense panel_rows × n_grid.
      ! Plus the n_RI × n_RI W_aux matrix. Distributed over n_procs ranks.
      use_cutoff = bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .AND. &
                   ALLOCATED(bs_env%ri_rs%chunk_centroids)
      CALL resolve_grid_panels(bs_env, bs_env%ri_rs%mat_phi_mu_l, pan_first, pan_last)
      mem_panels_GB = 0.0_dp
      DO ipan = 1, SIZE(pan_first)
         pan_rows = SUM(r_blk_sizes(pan_first(ipan):pan_last(ipan)))
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(bs_env%ri_rs%chunk_centroids, pan_first(ipan), &
                                             pan_last(ipan), bs_env%ri_rs%cutoff_radius_v_w, &
                                             grid_used)
            pan_width = SUM(r_blk_sizes, MASK=grid_used)
            CALL panel_template_elems(r_blk_sizes, bs_env%ri_rs%chunk_centroids, &
                                      grid_used, pan_first(ipan), pan_last(ipan), &
                                      bs_env%ri_rs%cutoff_radius_v_w, nze_tmpl)
            CALL panel_mem_estimate_GB(nze_tmpl, pan_rows, pan_width, n_grid_total, &
                                       n_RI, bs_env%n_ao, n_procs, mem_pan_GB)
         ELSE
            pan_width = n_grid_total
            mem_pan_GB = (3.0_dp*REAL(pan_rows, dp)*REAL(pan_width, dp) + &
                          2.0_dp*REAL(pan_rows, dp)*REAL(n_RI, dp))* &
                         8.0_dp/REAL(n_procs, dp)*1.0E-9_dp
         END IF
         mem_panels_GB = MAX(mem_panels_GB, mem_pan_GB)
      END DO
      mem_panels_GB = mem_panels_GB + &
                      REAL(n_RI, dp)*REAL(n_RI, dp)*8.0_dp/REAL(n_procs, dp)*1.0E-9_dp

      ! Z_lP SOLVE peak (compute_Z_lP). For the atom P with the largest integration
      ! sphere, one rank holds simultaneously:
      !   D_local  : n_local_grid x n_local_grid  (dense Gram, BLAS path only; O(n_local_grid^2))
      !   phi_local: n_local_grid x n_ao_used  (AOs reaching into the sphere only)
      !   d_lp     : n_local_grid x n_loc_ri, replicated once + one private copy per OMP thread
      ! n_local_grid = # grid points within cutoff_ri(P) = CUTOFF_RADIUS_RL_RI (if > 0) else
      ! r_c(RI metric) + r_RI(P). This is NOT evenly distributed: n_local_grid depends on the
      ! local density of atoms/grid, so the rank owning the densest atom peaks well above the
      ! average. We report the worst-case (max over atoms) as a per-rank upper bound.
      CALL get_qs_env(qs_env, particle_set=particle_set)
      natom = bs_env%n_atom

      max_n_local_grid = 0
      n_loc_ri_max = 0
      max_n_ao_used = 0
      DO iatom = 1, natom
         IF (bs_env%ri_rs%cutoff_radius_ri_rs > 0.0_dp) THEN
            cutoff_ri = bs_env%ri_rs%cutoff_radius_ri_rs
         ELSE
            cutoff_ri = bs_env%ri_metric%cutoff_radius + bs_env%ri_rs%radius_ri_per_atom(iatom)
         END IF
         pos_P(:) = particle_set(iatom)%r(:)
         n_local_grid = 0
         DO l = 1, n_grid_total
            IF (SUM((bs_env%ri_rs%grid_points(1:3, l) - pos_P(1:3))**2) <= cutoff_ri**2) THEN
               n_local_grid = n_local_grid + 1
            END IF
         END DO
         max_n_local_grid = MAX(max_n_local_grid, n_local_grid)
         CALL get_n_ao_in_sphere(bs_env, particle_set, iatom, cutoff_ri, n_ao_used_atom)
         max_n_ao_used = MAX(max_n_ao_used, n_ao_used_atom)
         n_loc_ri_max = MAX(n_loc_ri_max, &
                            bs_env%i_RI_end_from_atom(iatom) - bs_env%i_RI_start_from_atom(iatom) + 1)
      END DO

      n_procs_per_atom = MIN(MAX(bs_env%ri_rs%n_procs_per_atom_z_lp, 1), n_procs)
      n_threads = 1
!$    n_threads = omp_get_max_threads()

      ! D_local: dense on one rank for the BLAS path; block-cyclic over the subgroup (=> /G) for
      ! the ScaLAPACK path (N_PROCS_PER_ATOM_Z_LP = G > 1). phi_local/d_lp stay per-rank either way.
      IF (n_procs_per_atom > 1) THEN
         mem_D_local_GB = REAL(max_n_local_grid, dp)**2*8.0_dp/REAL(n_procs_per_atom, dp)*1.0E-9_dp
      ELSE
         mem_D_local_GB = REAL(max_n_local_grid, dp)**2*8.0_dp*1.0E-9_dp
      END IF
      mem_phi_local_GB = REAL(max_n_local_grid, dp)*REAL(max_n_ao_used, dp)*8.0_dp*1.0E-9_dp
      mem_dlp_GB = REAL(max_n_local_grid, dp)*REAL(n_loc_ri_max, dp)*8.0_dp* &
                   REAL(1 + n_threads, dp)*1.0E-9_dp
      mem_Zlp_peak_GB = mem_D_local_GB + mem_phi_local_GB + mem_dlp_GB

      ! Available memory per process = node MemLikelyFree / ranks-per-node, min across ranks
      ! (0 on non-Linux => warnings suppressed below).
      CALL ri_rs_mem_avail_per_proc_GB(bs_env, mem_avail_GB)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS memory estimate per MPI process:'
         WRITE (bs_env%unit_nr, '(T4,A,F37.2,A)') &
            'Available memory per process (system)', mem_avail_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F18.2,A)') &
            'Required for Z_lP (dense upper bound; actual is sparser)', mem_Z_lP_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F25.2,A)') &
            'Required for χ, W, Σ panels (peak per panel step)', mem_panels_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F17.2,A)') &
            'Required for Z_lP solve peak (D_local+φ, worst-case atom)', mem_Zlp_peak_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T69,I12)') &
            'Worst-case local-grid number of grid points:', max_n_local_grid
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'Worst-case memory D_local:', mem_D_local_GB, ' GB'

      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Z_lP_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated memory for Z_lP, ", mem_Z_lP_GB, " GB per process, exceeds the "// &
            "available ", mem_avail_GB, " GB. Z_lP (n_grid x n_RI) is distributed across all "// &
            "MPI ranks, so add nodes, use fewer MPI ranks per node, or raise "// &
            "N_PROCS_PER_ATOM_Z_LP to distribute each atom block via ScaLAPACK, which reduces "// &
            "the per-rank memory roughly by the number of ranks per atom."
         CPWARN(TRIM(msg))
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_panels_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the chi/W/Sigma panels, ", mem_panels_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, " GB. Panel memory "// &
            "scales roughly as 3*panel_size*n_grid/n_procs, so add nodes, use fewer MPI ranks "// &
            "per node, or raise N_PANELS for more but smaller panels."
         CPWARN(TRIM(msg))
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Zlp_peak_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the Z_lP solve, ", mem_Zlp_peak_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, &
            " GB. The per-atom Gram matrix D_local(n_local_grid, n_local_grid) dominates and "// &
            "scales as n_local_grid^2, and it is not balanced across ranks: the rank owning "// &
            "the atom with the largest integration sphere peaks well above the average. "// &
            "Either raise N_PROCS_PER_ATOM_Z_LP to distribute D_local block-cyclic via "// &
            "ScaLAPACK, which reduces that term roughly by the number of ranks per atom at no "// &
            "loss of accuracy, or lower CUTOFF_RADIUS_RL_RI, which shrinks D_local as "// &
            "n_local_grid^2 but trades accuracy, or use fewer MPI ranks per node so that each "// &
            "rank has more memory for the peak atom."
         CPWARN(TRIM(msg))
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_ri_rs_memory_estimate

! **************************************************************************************************
!> \brief Creates an empty (panel_chunks x neighborhood_chunks) DBCSR matrix with zero blocks
!>        pre-allocated only where |centroid(panel_row r) - centroid(column c)| <= cutoff.
!>        Used with retain_sparsity=.TRUE. in the subsequent dbcsr_multiply so distant blocks
!>        of the grid-basis panels (φ G φ^T, Z W Z^T, ...) are never computed at all.
!> \param L_pan ...
!> \param L_full ...
!> \param centroids ...
!> \param cutoff ...
!> \param blk0 ...
!> \param A_template ...
!> \param col_map ...
! **************************************************************************************************
   SUBROUTINE build_geo_template_panel(L_pan, L_full, centroids, cutoff, blk0, A_template, col_map)
      TYPE(dbcsr_type), INTENT(IN)                       :: L_pan, L_full
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER, INTENT(IN)                                :: blk0
      TYPE(dbcsr_type), INTENT(OUT)                      :: A_template
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: col_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_geo_template_panel'

      INTEGER                                            :: c, cg, cs, handle, my_pcol, my_prow, &
                                                            n_grid_blks, n_pan_blks, npcols, &
                                                            nprows, r, rs
      INTEGER, DIMENSION(:), POINTER                     :: grid_blk_sizes, pan_blk_sizes
      REAL(KIND=dp)                                      :: cutoff2
      REAL(KIND=dp), ALLOCATABLE                         :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      cutoff2 = cutoff**2
      CALL dbcsr_get_info(L_pan, nblkrows_total=n_pan_blks, row_blk_size=pan_blk_sizes)
      CALL dbcsr_get_info(L_full, nblkrows_total=n_grid_blks, row_blk_size=grid_blk_sizes)

      ! create_product_matrix assigns row r to process MOD(r-1,nprows) and
      ! col c to MOD(c-1,npcols), so we can determine local ownership analytically.
      CALL create_product_matrix(L_pan, L_full, 'N', 'T', A_template)
      CALL dbcsr_get_info(A_template, distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, npcols=npcols, &
                                  myprow=my_prow, mypcol=my_pcol)

      ALLOCATE (zero_blk(MAXVAL(pan_blk_sizes(1:n_pan_blks)), &
                         MAXVAL(grid_blk_sizes(1:n_grid_blks))))
      zero_blk(:, :) = 0.0_dp

      DO r = 1, n_pan_blks
         IF (MOD(r - 1, nprows) /= my_prow) CYCLE
         rs = pan_blk_sizes(r)
         DO c = 1, n_grid_blks
            IF (MOD(c - 1, npcols) /= my_pcol) CYCLE
            cg = c
            IF (PRESENT(col_map)) cg = col_map(c)
            IF ((centroids(1, blk0 + r - 1) - centroids(1, cg))**2 + &
                (centroids(2, blk0 + r - 1) - centroids(2, cg))**2 + &
                (centroids(3, blk0 + r - 1) - centroids(3, cg))**2 <= cutoff2) THEN
               cs = grid_blk_sizes(c)
               CALL dbcsr_put_block(A_template, r, c, zero_blk(1:rs, 1:cs))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(A_template)

      DEALLOCATE (zero_blk)
      CALL timestop(handle)

   END SUBROUTINE build_geo_template_panel

! **************************************************************************************************
!> \brief Slices a contiguous range of grid row-blocks [blk0, blk1] out of a (grid x n) DBCSR
!>        matrix into a new (P x n) panel matrix: iterate the source's local blocks, put the
!>        in-range ones into the panel with a remapped row-block index, then finalize. Row-block
!>        index i of the panel corresponds to source row-block blk0+i-1.
!> \param mat_full ...
!> \param blk0 ...
!> \param blk1 ...
!> \param mat_panel ...
! **************************************************************************************************
   SUBROUTINE extract_grid_panel(mat_full, blk0, blk1, mat_panel)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_full
      INTEGER, INTENT(IN)                                :: blk0, blk1
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_panel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'extract_grid_panel'

      INTEGER                                            :: handle, ib, jb, npb
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_full, col_dist_full, &
                                                            row_blk_full, row_blk_pan, &
                                                            row_dist_full, row_dist_pan
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_pan
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_full, distribution=dist_full, &
                          row_blk_size=row_blk_full, col_blk_size=col_blk_full)
      CALL dbcsr_distribution_get(dist_full, row_dist=row_dist_full, col_dist=col_dist_full)

      npb = blk1 - blk0 + 1
      ALLOCATE (row_dist_pan(npb), row_blk_pan(npb))
      row_dist_pan(:) = row_dist_full(blk0:blk1)
      row_blk_pan(:) = row_blk_full(blk0:blk1)

      CALL dbcsr_distribution_new(dist_pan, template=dist_full, &
                                  row_dist=row_dist_pan, col_dist=col_dist_full)
      CALL dbcsr_create(mat_panel, name="grid_panel", dist=dist_pan, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_pan, col_blk_size=col_blk_full)

      CALL dbcsr_iterator_start(iter, mat_full)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         IF (ib < blk0 .OR. ib > blk1) CYCLE
         CALL dbcsr_put_block(mat_panel, ib - blk0 + 1, jb, blk)
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(mat_panel)

      CALL dbcsr_distribution_release(dist_pan)
      DEALLOCATE (row_dist_pan, row_blk_pan)

      CALL timestop(handle)

   END SUBROUTINE extract_grid_panel

! **************************************************************************************************
!> \brief Marks which column blocks of a DBCSR matrix carry at least one non-zero block anywhere
!>        (global union). Used to restrict the inner index of the panel multiplies to the
!>        AO/RI atoms that actually touch the panel (exact: dropped rows only meet zeros).
!> \param matrix ...
!> \param para_env ...
!> \param used ...
! **************************************************************************************************
   SUBROUTINE collect_used_col_blocks(matrix, para_env, used)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_used_col_blocks'

      INTEGER                                            :: handle, ib, jb, nblkcols
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iused
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nblkcols_total=nblkcols)
      ALLOCATE (iused(nblkcols))
      iused(:) = 0

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         iused(jb) = 1
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL para_env%sum(iused)

      ALLOCATE (used(nblkcols))
      used(:) = (iused(:) > 0)
      DEALLOCATE (iused)

      CALL timestop(handle)

   END SUBROUTINE collect_used_col_blocks

! **************************************************************************************************
!> \brief Copies the flagged block rows (compress_rows=.TRUE.) or block columns (.FALSE.) of a
!>        DBCSR matrix into a compressed matrix. The subset keeps the parent's process assignment
!>        along the compressed dimension, so every block stays on its owning rank: the extraction
!>        is purely local (zero communication), like extract_grid_panel.
!> \param mat_full ...
!> \param used ...
!> \param mat_out ...
!> \param compress_rows ...
!> \param blk_map ...
! **************************************************************************************************
   SUBROUTINE extract_masked_blocks(mat_full, used, mat_out, compress_rows, blk_map)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_full
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: used
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_out
      LOGICAL, INTENT(IN)                                :: compress_rows
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: blk_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'extract_masked_blocks'

      INTEGER                                            :: handle, ib, jb, n_blk, n_sub, r
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inv_map
      INTEGER, DIMENSION(:), POINTER                     :: blk_full, blk_sub, col_blk_full, &
                                                            col_dist_full, dist_full_1d, &
                                                            dist_sub_1d, row_blk_full, &
                                                            row_dist_full
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_sub
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_full, distribution=dist_full, &
                          row_blk_size=row_blk_full, col_blk_size=col_blk_full)
      CALL dbcsr_distribution_get(dist_full, row_dist=row_dist_full, col_dist=col_dist_full)

      IF (compress_rows) THEN
         blk_full => row_blk_full
         dist_full_1d => row_dist_full
      ELSE
         blk_full => col_blk_full
         dist_full_1d => col_dist_full
      END IF
      n_blk = SIZE(blk_full)

      n_sub = COUNT(used)
      CPASSERT(n_sub > 0)
      ALLOCATE (inv_map(n_blk), blk_sub(n_sub), dist_sub_1d(n_sub))
      IF (PRESENT(blk_map)) ALLOCATE (blk_map(n_sub))
      inv_map(:) = 0
      r = 0
      DO ib = 1, n_blk
         IF (used(ib)) THEN
            r = r + 1
            inv_map(ib) = r
            blk_sub(r) = blk_full(ib)
            dist_sub_1d(r) = dist_full_1d(ib)
            IF (PRESENT(blk_map)) blk_map(r) = ib
         END IF
      END DO

      IF (compress_rows) THEN
         CALL dbcsr_distribution_new(dist_sub, template=dist_full, &
                                     row_dist=dist_sub_1d, col_dist=col_dist_full)
         CALL dbcsr_create(mat_out, name="row_subset", dist=dist_sub, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=blk_sub, col_blk_size=col_blk_full)
      ELSE
         CALL dbcsr_distribution_new(dist_sub, template=dist_full, &
                                     row_dist=row_dist_full, col_dist=dist_sub_1d)
         CALL dbcsr_create(mat_out, name="col_subset", dist=dist_sub, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_full, col_blk_size=blk_sub)
      END IF

      CALL dbcsr_iterator_start(iter, mat_full)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         IF (compress_rows) THEN
            IF (inv_map(ib) > 0) CALL dbcsr_put_block(mat_out, inv_map(ib), jb, blk)
         ELSE
            IF (inv_map(jb) > 0) CALL dbcsr_put_block(mat_out, ib, inv_map(jb), blk)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(mat_out)

      CALL dbcsr_distribution_release(dist_sub)
      DEALLOCATE (inv_map, blk_sub, dist_sub_1d)

      CALL timestop(handle)

   END SUBROUTINE extract_masked_blocks

! **************************************************************************************************
!> \brief Pre-seeds a square blocked DBCSR matrix with zero blocks only for block pairs whose
!>        centers lie within radius, for use with copy_fm_to_dbcsr(keep_sparsity=T) or
!>        dbcsr_multiply(retain_sparsity=T). Consumers: the CUTOFF_RADIUS_G_W operator truncation
!>        (atom-blocked, centers = atom_centers) and the RT-BSE CUTOFF_RADIUS_W0 truncation of the
!>        grid-basis W^0 (grid-blocked, centers = chunk_centroids).
!> \param matrix ...
!> \param centers block positions, one column per block row/column of matrix
!> \param radius truncation radius, same units as centers (bohr)
! **************************************************************************************************
   SUBROUTINE reserve_blocks_within_radius(matrix, centers, radius)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centers
      REAL(KIND=dp), INTENT(IN)                          :: radius

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reserve_blocks_within_radius'

      INTEGER                                            :: handle, i, j, my_pcol, my_prow, &
                                                            nblkcols, nblkrows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      REAL(KIND=dp)                                      :: radius2
      REAL(KIND=dp), ALLOCATABLE                         :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows, nblkcols_total=nblkcols, &
                          row_blk_size=row_blk, col_blk_size=col_blk, distribution=dist)
      CALL dbcsr_distribution_get(dist, row_dist=row_dist, col_dist=col_dist, &
                                  myprow=my_prow, mypcol=my_pcol)
      CPASSERT(nblkrows == SIZE(centers, 2))
      CPASSERT(nblkcols == SIZE(centers, 2))

      radius2 = radius**2
      ALLOCATE (zero_blk(MAXVAL(row_blk(1:nblkrows)), MAXVAL(col_blk(1:nblkcols))))
      zero_blk(:, :) = 0.0_dp

      DO i = 1, nblkrows
         IF (row_dist(i) /= my_prow) CYCLE
         DO j = 1, nblkcols
            IF (col_dist(j) /= my_pcol) CYCLE
            IF ((centers(1, i) - centers(1, j))**2 + (centers(2, i) - centers(2, j))**2 + &
                (centers(3, i) - centers(3, j))**2 <= radius2) THEN
               CALL dbcsr_put_block(matrix, i, j, zero_blk(1:row_blk(i), 1:col_blk(j)))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(matrix)

      DEALLOCATE (zero_blk)
      CALL timestop(handle)

   END SUBROUTINE reserve_blocks_within_radius

! **************************************************************************************************
!> \brief Creates the (empty) result matrix of op(mat_left) * op(mat_right) with the correct block
!>        structure and a distribution on the shared process grid, ready to be filled by
!>        dbcsr_multiply. Row structure comes from op(left), column structure from op(right).
!> \param mat_left ...
!> \param mat_right ...
!> \param transa 'N' or 'T' applied to mat_left
!> \param transb 'N' or 'T' applied to mat_right
!> \param mat_out ...
! **************************************************************************************************
   SUBROUTINE create_product_matrix(mat_left, mat_right, transa, transb, mat_out)

      TYPE(dbcsr_type), INTENT(IN)                       :: mat_left, mat_right
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_out

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_product_matrix'

      INTEGER                                            :: handle, i, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_l, col_blk_r, out_col_blk, &
                                                            out_col_dist, out_row_blk, &
                                                            out_row_dist, row_blk_l, row_blk_r
      TYPE(dbcsr_distribution_type)                      :: dist_l, dist_out

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_left, distribution=dist_l, row_blk_size=row_blk_l, col_blk_size=col_blk_l)
      CALL dbcsr_get_info(mat_right, row_blk_size=row_blk_r, col_blk_size=col_blk_r)
      CALL dbcsr_distribution_get(dist_l, nprows=nprows, npcols=npcols)

      ! block SIZES follow op(left)/op(right); DISTRIBUTIONS are freshly round-robined onto the
      ! shared process grid (a transposed operand's row-dist is NOT a valid col-dist on a
      ! non-square grid). dbcsr_multiply redistributes internally, so any valid mapping works.
      IF (transa == 'N') THEN
         out_row_blk => row_blk_l
      ELSE
         out_row_blk => col_blk_l
      END IF
      IF (transb == 'N') THEN
         out_col_blk => col_blk_r
      ELSE
         out_col_blk => row_blk_r
      END IF

      ALLOCATE (out_row_dist(SIZE(out_row_blk)), out_col_dist(SIZE(out_col_blk)))
      DO i = 1, SIZE(out_row_blk)
         out_row_dist(i) = MOD(i - 1, nprows)
      END DO
      DO i = 1, SIZE(out_col_blk)
         out_col_dist(i) = MOD(i - 1, npcols)
      END DO

      CALL dbcsr_distribution_new(dist_out, template=dist_l, &
                                  row_dist=out_row_dist, col_dist=out_col_dist)
      CALL dbcsr_create(mat_out, name="panel_product", dist=dist_out, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=out_row_blk, col_blk_size=out_col_blk)
      CALL dbcsr_distribution_release(dist_out)
      DEALLOCATE (out_row_dist, out_col_dist)

      CALL timestop(handle)

   END SUBROUTINE create_product_matrix

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

      CALL setup_square_topology(template, dist_ao_ao, blk_ao, dist_row_ao)
      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      ! Optional CUTOFF_RADIUS_G_W operator truncation: only atom-pair blocks within the radius
      ! are reserved and filled.
      IF (bs_env%ri_rs%cutoff_radius_g_w > 0.0_dp .AND. &
          ALLOCATED(bs_env%ri_rs%atom_centers)) THEN
         CALL reserve_blocks_within_radius(matrix_G_ao, bs_env%ri_rs%atom_centers, &
                                           bs_env%ri_rs%cutoff_radius_g_w)
         CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.TRUE.)
      ELSE
         CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)
      END IF
      CALL dbcsr_filter(matrix_G_ao, bs_env%eps_filter)

      ! release only the topology; keep matrix_G_ao for the caller
      CALL release_square_topology(dist=dist_ao_ao, mapped_dist=dist_row_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_ao

! **************************************************************************************************
!> \brief Panel-streaming evaluation of  out += scale * L_out^T (A_grid ∘ B_grid) L_out,
!>        where A_grid = L_A M_A L_A^T and B_grid = L_B M_B L_B^T, WITHOUT ever forming the full
!>        grid x grid objects. The grid (row) index is processed in panels of ~panel_size rows; for
!>        each panel only P x grid slabs are built, Hadamard-multiplied, and contracted into the
!>        (small) output. Algebraically identical to L_out^T (A_grid ∘ B_grid) L_out summed over
!>        grid rows, so the result matches the non-streamed path to eps_filter.
!>
!>        Mapping (L in {phi (grid x AO), Z (grid x RI)}, M the AO/RI-space operator):
!>          chi :   L_A=L_B=phi, M_A=G_occ_ao, M_B=G_vir_ao, L_out=Z   -> RI x RI
!>          Sig :   L_A=phi (M_A=D/G), L_B=Z (M_B=V/W),       L_out=phi -> AO x AO
!> \param L_A ...
!> \param M_A ...
!> \param L_B ...
!> \param M_B ...
!> \param L_out ...
!> \param mat_out ...
!> \param scale ...
!> \param eps ...
!> \param para_env ...
!> \param pan_first ...
!> \param pan_last ...
!> \param lb_eq_la ...
!> \param lout_eq_la ...
!> \param zero_out ...
!> \param keep_sparsity ...
!> \param centroids ...
!> \param cutoff ...
!> \param grid_occupation ...
! **************************************************************************************************
   SUBROUTINE contract_grid_panels(L_A, M_A, L_B, M_B, L_out, mat_out, scale, eps, para_env, &
                                   pan_first, pan_last, lb_eq_la, lout_eq_la, zero_out, &
                                   keep_sparsity, centroids, cutoff, grid_occupation)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_A
      TYPE(dbcsr_type), INTENT(INOUT)                    :: M_A
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_B
      TYPE(dbcsr_type), INTENT(INOUT)                    :: M_B
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_out
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_out
      REAL(KIND=dp), INTENT(IN)                          :: scale, eps
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pan_first, pan_last
      LOGICAL, INTENT(IN)                                :: lb_eq_la, lout_eq_la, zero_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: grid_occupation

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, ipan, n_grid_total, &
                                                            ncols_pan, nrows_pan
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap
      LOGICAL                                            :: my_keep_sparsity, use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, usedA, usedB
      TYPE(dbcsr_type)                                   :: A_pan, B_pan, C_pan, LA_pan, LA_panC, &
                                                            LB_pan, LB_panC, Lout_pan, MA_sub, &
                                                            MB_sub, tmp2, tmpA, tmpB
      TYPE(dbcsr_type), POINTER                          :: RB_A, RB_B, RB_out
      TYPE(dbcsr_type), TARGET                           :: LA_near, LB_near, Lout_near

      CALL timeset(routineN, handle)

      my_keep_sparsity = .FALSE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity
      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp
      IF (PRESENT(grid_occupation)) grid_occupation = 0.0_dp

      CALL dbcsr_get_info(L_A, nfullrows_total=n_grid_total)

      IF (zero_out) CALL dbcsr_set(mat_out, 0.0_dp)

      DO ipan = 1, SIZE(pan_first)
         blk0 = pan_first(ipan)
         blk1 = pan_last(ipan)

         ! phi/Z panel slices (P x n)
         CALL extract_grid_panel(L_A, blk0, blk1, LA_pan)
         IF (.NOT. lb_eq_la) CALL extract_grid_panel(L_B, blk0, blk1, LB_pan)
         IF (.NOT. lout_eq_la) CALL extract_grid_panel(L_out, blk0, blk1, Lout_pan)

         ! Which AO/RI atoms (column blocks) actually touch this panel: the inner index of
         ! every multiply below is restricted to them, so only the matching rows of the
         ! system-wide operators M_A/M_B ever enter Cannon (exact: dropped rows meet zeros).
         CALL collect_used_col_blocks(LA_pan, para_env, usedA)
         IF (.NOT. lb_eq_la) THEN
            CALL collect_used_col_blocks(LB_pan, para_env, usedB)
         ELSE
            IF (ALLOCATED(usedB)) DEALLOCATE (usedB)
            ALLOCATE (usedB, SOURCE=usedA)
         END IF
         IF (.NOT. (ANY(usedA) .AND. ANY(usedB))) THEN
            ! empty panel slice: its Hadamard contribution is exactly zero on all ranks
            CALL dbcsr_release(LA_pan)
            IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_pan)
            IF (.NOT. lout_eq_la) CALL dbcsr_release(Lout_pan)
            CYCLE
         END IF

         ! Grid rows within reach of the panel: with the CUTOFF_RADIUS_RL_W truncation only
         ! they can appear as columns of the panel products / inner rows of the L_out multiply,
         ! so the system-wide phi/Z right operands are cut down to this neighborhood slice
         ! (local extraction, zero communication; exact w.r.t. the geo template).
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, grid_used)
            CALL extract_masked_blocks(L_A, grid_used, LA_near, compress_rows=.TRUE., blk_map=gmap)
            IF (.NOT. lb_eq_la) CALL extract_masked_blocks(L_B, grid_used, LB_near, compress_rows=.TRUE.)
            IF (.NOT. lout_eq_la) CALL extract_masked_blocks(L_out, grid_used, Lout_near, compress_rows=.TRUE.)
            RB_A => LA_near
         ELSE
            RB_A => L_A
         END IF
         IF (lb_eq_la) THEN
            RB_B => RB_A
         ELSE IF (use_cutoff) THEN
            RB_B => LB_near
         ELSE
            RB_B => L_B
         END IF
         IF (lout_eq_la) THEN
            RB_out => RB_A
         ELSE IF (use_cutoff) THEN
            RB_out => Lout_near
         ELSE
            RB_out => L_out
         END IF

         ! A_pan = LA_pan * M_A * L_A^T   (P x grid_near).
         ! When cutoff is active, A_pan is pre-seeded with only nearby blocks via
         ! build_geo_template_panel, and the multiply uses retain_sparsity to skip
         ! computing distant blocks entirely (exact: they are zero by locality).
         CALL extract_masked_blocks(LA_pan, usedA, LA_panC, compress_rows=.FALSE.)
         CALL extract_masked_blocks(M_A, usedA, MA_sub, compress_rows=.TRUE.)
         CALL create_product_matrix(LA_panC, MA_sub, 'N', 'N', tmpA)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, LA_panC, MA_sub, 0.0_dp, tmpA, filter_eps=eps)
         CALL dbcsr_release(MA_sub)
         IF (use_cutoff) THEN
            CALL build_geo_template_panel(LA_pan, LA_near, centroids, cutoff, blk0, A_pan, &
                                          col_map=gmap)
         ELSE
            CALL create_product_matrix(tmpA, RB_A, 'N', 'T', A_pan)
         END IF
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_A, 0.0_dp, A_pan, &
                             filter_eps=eps, retain_sparsity=use_cutoff)
         CALL dbcsr_release(tmpA)

         ! Grid-basis occupation of A_pan = φ G φ^T, accumulated over ALL panels into the
         ! occupation of the (never formed) full grid x grid object:
         !   sum_panels nnz(A_pan) / n_grid^2,  with nnz = occ * pan_rows * pan_cols.
         ! Panel-independent by construction -- a single-panel sample would instead report the
         ! local neighbor count of whichever region happens to land in that panel.
         IF (PRESENT(grid_occupation)) THEN
            CALL dbcsr_get_info(A_pan, nfullrows_total=nrows_pan, nfullcols_total=ncols_pan)
            grid_occupation = grid_occupation + dbcsr_get_occupation(A_pan)* &
                              REAL(ncols_pan, dp)*REAL(nrows_pan, dp)/ &
                              (REAL(n_grid_total, dp)*REAL(n_grid_total, dp))
         END IF

         ! B_pan = LB_pan * M_B * L_B^T   (P x grid_near); reuse the L_A slices when L_B == L_A.
         ! With keep_sparsity, B_pan is pre-populated with A_pan's block structure so that
         ! retain_sparsity forces the final multiply to fill only those blocks (exact for ∘).
         IF (lb_eq_la) THEN
            CALL extract_masked_blocks(M_B, usedA, MB_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LA_panC, MB_sub, 'N', 'N', tmpB)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LA_panC, MB_sub, 0.0_dp, tmpB, filter_eps=eps)
         ELSE
            CALL extract_masked_blocks(LB_pan, usedB, LB_panC, compress_rows=.FALSE.)
            CALL extract_masked_blocks(M_B, usedB, MB_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LB_panC, MB_sub, 'N', 'N', tmpB)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LB_panC, MB_sub, 0.0_dp, tmpB, filter_eps=eps)
            CALL dbcsr_release(LB_panC)
         END IF
         CALL dbcsr_release(MB_sub)
         IF (my_keep_sparsity) THEN
            CALL dbcsr_create(B_pan, template=A_pan)
            CALL dbcsr_copy(B_pan, A_pan)
            CALL dbcsr_set(B_pan, 0.0_dp)
            ! The output pattern is already fixed. Omitting redundant on-the-fly filtering also
            ! avoids overflowing DBCSR's single-precision screening norms for conditioned Z fits.
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_B, 0.0_dp, B_pan, &
                                retain_sparsity=.TRUE.)
         ELSE
            CALL create_product_matrix(tmpB, RB_B, 'N', 'T', B_pan)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_B, 0.0_dp, B_pan, filter_eps=eps)
         END IF
         CALL dbcsr_release(tmpB)
         CALL dbcsr_release(LA_panC)

         ! C_pan = scale * (A_pan ∘ B_pan)   (P x grid_near)
         CALL dbcsr_create(C_pan, template=A_pan)
         CALL hadamard_product(A_pan, B_pan, C_pan, scale)
         CALL dbcsr_release(A_pan)
         CALL dbcsr_release(B_pan)

         ! tmp2 = C_pan * L_out   (P x n_out; inner index restricted to the neighborhood)
         CALL create_product_matrix(C_pan, RB_out, 'N', 'N', tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, RB_out, 0.0_dp, tmp2, filter_eps=eps)
         CALL dbcsr_release(C_pan)

         ! mat_out += L_out_pan^T * tmp2   (accumulate: beta = 1)
         IF (lout_eq_la) THEN
            CALL dbcsr_multiply('T', 'N', 1.0_dp, LA_pan, tmp2, 1.0_dp, mat_out, filter_eps=eps)
         ELSE
            CALL dbcsr_multiply('T', 'N', 1.0_dp, Lout_pan, tmp2, 1.0_dp, mat_out, filter_eps=eps)
            CALL dbcsr_release(Lout_pan)
         END IF
         CALL dbcsr_release(tmp2)
         IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_pan)
         CALL dbcsr_release(LA_pan)
         IF (use_cutoff) THEN
            CALL dbcsr_release(LA_near)
            IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_near)
            IF (.NOT. lout_eq_la) CALL dbcsr_release(Lout_near)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE contract_grid_panels

! **************************************************************************************************
!> \brief Σ^c-specific panel loop: computes both the occupied (neg) and virtual (pos) contributions
!>        in a single pass over grid panels, forming W_pan = Z_panel × W_aux × Z^T only ONCE per
!>        panel and reusing it for both the G^occ and G^vir Hadamard contractions.
!>
!>        Computes:
!>          mat_Sigma_neg = φ^T ( (φ G^occ φ^T) ∘ (Z W^MIC Z^T) ) φ
!>          mat_Sigma_pos = φ^T ( (φ G^vir φ^T) ∘ (Z W^MIC Z^T) ) φ
!>
!> \param mat_phi ...
!> \param mat_Z ...
!> \param mat_G_occ_ao ...
!> \param mat_G_vir_ao ...
!> \param mat_W_aux ...
!> \param mat_Sigma_neg ...
!> \param mat_Sigma_pos ...
!> \param eps ...
!> \param para_env ...
!> \param pan_first ...
!> \param pan_last ...
!> \param keep_sparsity ...
!> \param centroids ...
!> \param cutoff ...
! **************************************************************************************************
   SUBROUTINE contract_grid_panels_sigma_c(mat_phi, mat_Z, mat_G_occ_ao, mat_G_vir_ao, &
                                           mat_W_aux, mat_Sigma_neg, mat_Sigma_pos, eps, para_env, &
                                           pan_first, pan_last, keep_sparsity, centroids, cutoff)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: mat_phi, mat_Z
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_G_occ_ao, mat_G_vir_ao, mat_W_aux, &
                                                            mat_Sigma_neg, mat_Sigma_pos
      REAL(KIND=dp), INTENT(IN)                          :: eps
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pan_first, pan_last
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_grid_panels_sigma_c'

      INTEGER                                            :: blk0, blk1, handle, ipan
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap
      LOGICAL                                            :: my_keep_sparsity, use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, used_ao, used_ri
      TYPE(dbcsr_type)                                   :: A_occ, A_vir, C_pan, G_occ_sub, &
                                                            G_vir_sub, phi_pan, phi_panC, tmp2, &
                                                            tmpA, tmpB, W_pan, W_sub, Z_pan, Z_panC
      TYPE(dbcsr_type), POINTER                          :: RB_phi, RB_Z
      TYPE(dbcsr_type), TARGET                           :: phi_near, Z_near

      CALL timeset(routineN, handle)

      my_keep_sparsity = .FALSE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity
      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp

      CALL dbcsr_set(mat_Sigma_neg, 0.0_dp)
      CALL dbcsr_set(mat_Sigma_pos, 0.0_dp)

      DO ipan = 1, SIZE(pan_first)
         blk0 = pan_first(ipan)
         blk1 = pan_last(ipan)

         CALL extract_grid_panel(mat_phi, blk0, blk1, phi_pan)
         CALL extract_grid_panel(mat_Z, blk0, blk1, Z_pan)

         ! AO/RI atoms touching this panel: only the matching rows of G_occ/G_vir/W ever
         ! enter the multiplies below (exact: dropped rows meet zero columns of the panel).
         CALL collect_used_col_blocks(phi_pan, para_env, used_ao)
         CALL collect_used_col_blocks(Z_pan, para_env, used_ri)
         IF (.NOT. (ANY(used_ao) .AND. ANY(used_ri))) THEN
            CALL dbcsr_release(phi_pan)
            CALL dbcsr_release(Z_pan)
            CYCLE
         END IF

         ! Neighborhood slices of phi/Z (grid rows within cutoff of the panel): they replace
         ! the system-wide right operands in every multiply (local extraction, zero comm).
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, grid_used)
            CALL extract_masked_blocks(mat_phi, grid_used, phi_near, compress_rows=.TRUE., blk_map=gmap)
            CALL extract_masked_blocks(mat_Z, grid_used, Z_near, compress_rows=.TRUE.)
            RB_phi => phi_near
            RB_Z => Z_near
         ELSE
            RB_phi => mat_phi
            RB_Z => mat_Z
         END IF

         CALL extract_masked_blocks(phi_pan, used_ao, phi_panC, compress_rows=.FALSE.)
         CALL extract_masked_blocks(Z_pan, used_ri, Z_panC, compress_rows=.FALSE.)
         CALL extract_masked_blocks(mat_G_occ_ao, used_ao, G_occ_sub, compress_rows=.TRUE.)
         CALL extract_masked_blocks(mat_G_vir_ao, used_ao, G_vir_sub, compress_rows=.TRUE.)
         CALL extract_masked_blocks(mat_W_aux, used_ri, W_sub, compress_rows=.TRUE.)

         ! A_occ = phi_pan × G_occ × phi^T  (built first so W_pan can inherit its pattern).
         ! With cutoff active, A_occ is pre-seeded with geo-local blocks so that the
         ! phi^T multiply uses retain_sparsity and never computes distant blocks.
         CALL create_product_matrix(phi_panC, G_occ_sub, 'N', 'N', tmpA)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_occ_sub, 0.0_dp, tmpA, filter_eps=eps)
         IF (use_cutoff) THEN
            CALL build_geo_template_panel(phi_pan, phi_near, centroids, cutoff, blk0, A_occ, &
                                          col_map=gmap)
         ELSE
            CALL create_product_matrix(tmpA, RB_phi, 'N', 'T', A_occ)
         END IF
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_phi, 0.0_dp, A_occ, &
                             filter_eps=eps, retain_sparsity=use_cutoff)
         CALL dbcsr_release(tmpA)
         CALL dbcsr_release(G_occ_sub)

         ! A_vir = phi_pan × G_vir × phi^T  (same pre-screen as A_occ)
         CALL create_product_matrix(phi_panC, G_vir_sub, 'N', 'N', tmpA)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_vir_sub, 0.0_dp, tmpA, filter_eps=eps)
         IF (use_cutoff) THEN
            CALL build_geo_template_panel(phi_pan, phi_near, centroids, cutoff, blk0, A_vir, &
                                          col_map=gmap)
         ELSE
            CALL create_product_matrix(tmpA, RB_phi, 'N', 'T', A_vir)
         END IF
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_phi, 0.0_dp, A_vir, &
                             filter_eps=eps, retain_sparsity=use_cutoff)
         CALL dbcsr_release(tmpA)
         CALL dbcsr_release(G_vir_sub)

         ! W_pan = Z_pan × W_aux × Z^T  (computed once, reused for both Σ^c terms).
         ! With keep_sparsity, W_pan is pre-seeded with the union of A_occ and A_vir block
         ! patterns so that retain_sparsity forces the multiply to fill only those blocks:
         ! exact since W outside G_occ∪G_vir is multiplied by zero in the Hadamard.
         CALL create_product_matrix(Z_panC, W_sub, 'N', 'N', tmpB)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, Z_panC, W_sub, 0.0_dp, tmpB, filter_eps=eps)
         IF (my_keep_sparsity) THEN
            CALL dbcsr_create(W_pan, template=A_occ)
            CALL dbcsr_copy(W_pan, A_occ)
            CALL dbcsr_add(W_pan, A_vir, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(W_pan, 0.0_dp)
            ! The retained pattern is exact for the following Hadamard products; screening it
            ! again is redundant and may overflow DBCSR's single-precision block-norm product.
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_Z, 0.0_dp, W_pan, &
                                retain_sparsity=.TRUE.)
         ELSE
            CALL create_product_matrix(tmpB, RB_Z, 'N', 'T', W_pan)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_Z, 0.0_dp, W_pan, filter_eps=eps)
         END IF
         CALL dbcsr_release(tmpB)
         CALL dbcsr_release(W_sub)
         CALL dbcsr_release(phi_panC)
         CALL dbcsr_release(Z_panC)

         ! Σ^c_neg: φ^T ( A_occ ∘ W_pan ) φ
         CALL dbcsr_create(C_pan, template=A_occ)
         CALL hadamard_product(A_occ, W_pan, C_pan, 1.0_dp)
         CALL dbcsr_release(A_occ)
         CALL create_product_matrix(C_pan, RB_phi, 'N', 'N', tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, RB_phi, 0.0_dp, tmp2, filter_eps=eps)
         CALL dbcsr_release(C_pan)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, phi_pan, tmp2, 1.0_dp, mat_Sigma_neg, filter_eps=eps)
         CALL dbcsr_release(tmp2)

         ! Σ^c_pos: φ^T ( A_vir ∘ W_pan ) φ  — W_pan reused
         CALL dbcsr_create(C_pan, template=A_vir)
         CALL hadamard_product(A_vir, W_pan, C_pan, 1.0_dp)
         CALL dbcsr_release(A_vir)
         CALL create_product_matrix(C_pan, RB_phi, 'N', 'N', tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, RB_phi, 0.0_dp, tmp2, filter_eps=eps)
         CALL dbcsr_release(C_pan)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, phi_pan, tmp2, 1.0_dp, mat_Sigma_pos, filter_eps=eps)
         CALL dbcsr_release(tmp2)

         CALL dbcsr_release(W_pan)
         CALL dbcsr_release(Z_pan)
         CALL dbcsr_release(phi_pan)
         IF (use_cutoff) THEN
            CALL dbcsr_release(phi_near)
            CALL dbcsr_release(Z_near)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE contract_grid_panels_sigma_c

! **************************************************************************************************
!> \brief Computes C = A ◦ B (Element-wise Hadamard product) for sparse DBCSR matrices.
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_C ...
!> \param fac (Scaling factor applied to the product)
! **************************************************************************************************
   SUBROUTINE hadamard_product(matrix_A, matrix_B, matrix_C, fac)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: fac

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hadamard_product'

      INTEGER                                            :: col, handle, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_B, blk_C
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_C, matrix_A)

      CALL dbcsr_iterator_start(iter, matrix_C)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_C)

         CALL dbcsr_get_block_p(matrix_B, row, col, blk_B, found)

         IF (found) THEN
            blk_C(:, :) = fac*blk_C(:, :)*blk_B(:, :)
         ELSE
            ! If B is sparse here, the product is zero
            blk_C(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE hadamard_product

! **************************************************************************************************
!> \brief Computes the screened Coulomb interaction on the imaginary-time grid, entirely in the
!>        RI auxiliary (PQ) basis:
!>          χ_PQ(iω) = Σ_τ w(ω,τ) cos(ωτ) χ_PQ(iτ)                  (cosine transform)
!>          ε(iω)    = Id - V^0.5 M^-1 χ(iω) M^-1 V^0.5             (dielectric function)
!>          W(iω)    = V^0.5 ( ε^-1(iω) - Id ) V^0.5                (correlation part only)
!>          W(iτ)    = Σ_ω w̃(τ,ω) cos(ωτ) W(iω)                    (back transform)
!>          W(iτ)   <- M^-1 W(iτ) M^-1                              (fold in the RI metric)
!>        where V is the bare Coulomb matrix and M the RI metric.
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_time ...
! **************************************************************************************************
   SUBROUTINE compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W'

      INTEGER                                            :: handle, i_t, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type)                                   :: fm_M_inv_V_sqrt, fm_V, fm_V_sqrt

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL create_fm_W_MIC_time(bs_env, fm_W_time)

      ! 1. Allocate V and M matrices
      CALL cp_fm_create(fm_V, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_M_inv_V_sqrt, bs_env%fm_RI_RI%matrix_struct)

      ! Compute V and M^-1 * V^0.5
      CALL compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_M_inv_V_sqrt)

      ! 2. Loop over frequencies
      DO j_w = 1, bs_env%num_time_freq_points
         ! Fourier transformation of χ_PQ(iτ) to χ_PQ(iω_j)
         CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

         ! ε(iω_j) = Id - V^0.5*M^-1*χ(iω_j)*M^-1*V^0.5
         ! W(iω_j) = V^0.5*(ε^-1(iω_j)-Id)*V^0.5
         CALL compute_fm_W_freq(bs_env, bs_env%fm_chi_Gamma_freq, fm_V_sqrt, &
                                fm_M_inv_V_sqrt, bs_env%fm_W_MIC_freq)

         ! Fourier transform from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
         CALL Fourier_transform_w_to_t(bs_env, fm_W_time, bs_env%fm_W_MIC_freq, j_w)
      END DO

      ! M^-1(k=0) W^MIC(iτ) M^-1(k=0) -> fm_W_time
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_W_time)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed W(iτ),', ' Execution time', m_walltime() - t1, ' s'
      END IF

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ! Cleanup
      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_M_inv_V_sqrt)

      ! Marek : Fourier transform W^MIC(itau) back to get it at a specific im.frequency point - iomega = 0
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         t1 = m_walltime()
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         ! Set to zero
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
         ! Sum over all times
         DO i_t = 1, bs_env%num_time_freq_points
            ! Add the relevant structure with correct weight
            CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                     bs_env%imag_time_weights_freq_zero(i_t), fm_W_time(i_t))
         END DO
         ! Done, save to file
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         ! Report calculation
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed W(0),', ' Execution time', m_walltime() - t1, ' s'
         END IF
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 V^0.5 for the RI-RS dielectric function.
!>        The Coulomb matrix V is constructed by the RI-RS k-point path. The inverse metric
!>        M^-1(k=0) is precomputed once in gw_utils and stored in bs_env.
!> \param bs_env GW environment containing M^-1(k=0) and the RI basis
!> \param qs_env electronic-structure environment used to evaluate V(k=0)
!> \param fm_V Coulomb matrix V(k=0)
!> \param fm_V_sqrt symmetric factor V^0.5
!> \param fm_Minv_Vsqrt product M^-1(k=0) V^0.5
! **************************************************************************************************
   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (bs_env%auto_ri%enabled) THEN
         ! -------------------------------------------------------------------
         ! 1a. The optimized AB functions span two atoms, so their previously transformed
         !     Coulomb matrix cannot be rebuilt by the atom-local k-point integral routine.
         ! -------------------------------------------------------------------
         CALL cp_fm_to_fm(bs_env%auto_ri%coulomb_matrix, fm_V)
      ELSE
         ! -------------------------------------------------------------------
         ! 1b. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1.
         ! -------------------------------------------------------------------
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                         qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

         ALLOCATE (mat_V_kp(1:1, 1:2))
         NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)
         CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
         ! The dummy imaginary part is required only by the k-point routine interface.
         CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
         CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                         particle_set, qs_kind_set, atomic_kind_set, &
                                         bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)
         CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
         DEALLOCATE (mat_V_kp)
      END IF

      ! -----------------------------------------------------------------------
      ! 2. V -> V^0.5.
      ! -----------------------------------------------------------------------
      CALL fm_sqrt(fm_V, fm_V_sqrt, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! -----------------------------------------------------------------------
      ! 3. M^-1(k=0) V^0.5.
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI, 1.0_dp, &
                         bs_env%fm_Minv_Gamma, fm_V_sqrt, 0.0_dp, fm_Minv_Vsqrt)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

! **************************************************************************************************
!> \brief Computes the screened interaction at one imaginary frequency:
!>          ε(iω_j)  = Id - (M^-1 V^0.5)^T χ(iω_j) (M^-1 V^0.5)
!>          W(iω_j)  = V^0.5^T ( ε^-1(iω_j) - Id ) V^0.5
!>        ε is inverted via Cholesky; if that fails due to conditioning, via
!>        eigendecomposition (cp_fm_power) with eigenvalue filtering.
!> \param bs_env ...
!> \param fm_chi_freq_j ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
!> \param fm_W_freq_j ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_freq(bs_env, fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt, fm_W_freq_j)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_W_freq_j

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_fm_W_freq'

      INTEGER                                            :: handle, n_RI
      TYPE(cp_fm_type)                                   :: fm_eps_freq_j, fm_work

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_eps_freq_j, fm_chi_freq_j%matrix_struct)
      CALL cp_fm_create(fm_work, fm_chi_freq_j%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. ε(iω_j) = Id - (M^-1 * V^0.5)^T * χ(iω_j) * (M^-1 * V^0.5)
      ! -----------------------------------------------------------------------
      ! work = χ(iω_j) * (M^-1 * V^0.5)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_chi_freq_j, fm_Minv_Vsqrt, 0.0_dp, fm_work)

      ! eps_work = (M^-1 * V^0.5)^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_Minv_Vsqrt, fm_work, 0.0_dp, fm_eps_freq_j)

      ! ε(iω_j) = Id - eps_work  -->  -eps_work + Id
      CALL fm_add_on_diag(fm_eps_freq_j, 1.0_dp)

      ! Force perfect symmetry before Cholesky to avoid info != 0 due to GEMM noise
      CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)

      ! -----------------------------------------------------------------------
      ! 2. W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      ! -----------------------------------------------------------------------

      ! a) Invert ε by Cholesky decomposition or, if that fails, by diagonalization.
      CALL fm_invert(fm_eps_freq_j, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! b) ε^-1(iω_j) - Id
      CALL fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! c) W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      CALL fm_contract_ABA(fm_V_sqrt, fm_eps_freq_j, fm_W_freq_j)

      ! Cleanup
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_eps_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_freq

! **************************************************************************************************
!> \brief Adds a real scalar value to the diagonal of a real full matrix
!> \param fm ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE fm_add_on_diag(fm, alpha)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm
      REAL(KIND=dp), INTENT(IN)                          :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_add_on_diag'

      INTEGER                                            :: handle, i_global, i_row, j_col, &
                                                            j_global, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            i_global = row_indices(i_row)
            IF (j_global == i_global) THEN
               fm%local_data(i_row, j_col) = fm%local_data(i_row, j_col) + alpha
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_add_on_diag

! **************************************************************************************************
!> \brief Computes the exact-exchange part of the GW self-energy:
!>          D_μν        = Σ_n^occ C_μn C_νn                     (density matrix = G^occ at τ=0)
!>          V^tr_PQ     = M^-1 (P|Q)_trunc M^-1                 (truncated Coulomb, RI basis)
!>          Σ^x_λσ(k=0) = -Σ_ll' Φ_λ(r_l) [ (φ D φ^T)_ll' ∘ (Z V^tr Z^T)_ll' ] Φ_σ(r_l')
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
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma, matrix_D_ao, &
                                                            matrix_V_aux

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      CALL dbcsr_create(mat_Sigma_x_Gamma, template=bs_env%mat_ao_ao%matrix)

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

      ! =========================================================================
      ! 1. COMPUTE V^tr_PQ (RI x RI)
      ! =========================================================================
      CALL setup_square_topology(mat_Z_lP, dist_aux_aux, blk_aux, dist_row_aux)

      IF (bs_env%auto_ri%enabled) THEN
         ALLOCATE (fm_Vtr_Gamma(1, 1))
         CALL cp_fm_create(fm_Vtr_Gamma(1, 1), bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_to_fm(bs_env%auto_ri%coulomb_matrix, fm_Vtr_Gamma(1, 1))
      ELSE
         CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                                 bs_env%trunc_coulomb, do_kpoints=.FALSE.)
      END IF

      ! M^-1(k=0) V^tr(τ) M^-1(k=0) -> fm_Vtr_Gamma
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      ! Optional CUTOFF_RADIUS_G_W operator truncation + filter (see build_G_ao)
      IF (bs_env%ri_rs%cutoff_radius_g_w > 0.0_dp .AND. ALLOCATED(bs_env%ri_rs%atom_centers)) THEN
         CALL reserve_blocks_within_radius(matrix_V_aux, bs_env%ri_rs%atom_centers, &
                                           bs_env%ri_rs%cutoff_radius_g_w)
         CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.TRUE.)
      ELSE
         CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)
      END IF
      CALL dbcsr_filter(matrix_V_aux, bs_env%eps_filter)

      ! =========================================================================
      ! 2. SPIN LOOP FOR EXACT EXCHANGE
      ! Σ^x_λσ = -Σ_ll' Φ_λ(r_l) ( D_ll' V^tr_ll' ) Φ_σ(r_l')
      !        = -φ^T ( (φ D φ^T) ∘ (Z V^tr Z^T) ) φ
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! AO-space density matrix D_µν = G^occ at τ = 0
         CALL build_G_ao(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_D_ao)

         CALL contract_grid_panels(L_A=mat_phi_mu_l, M_A=matrix_D_ao, &
                                   L_B=mat_Z_lP, M_B=matrix_V_aux, &
                                   L_out=mat_phi_mu_l, mat_out=mat_Sigma_x_Gamma, &
                                   scale=1.0_dp, eps=bs_env%eps_filter, &
                                   para_env=bs_env%para_env, &
                                   pan_first=pan_first, pan_last=pan_last, &
                                   lb_eq_la=.FALSE., lout_eq_la=.TRUE., zero_out=.TRUE., &
                                   keep_sparsity=bs_env%ri_rs%keep_sparsity_rirs, &
                                   centroids=bs_env%ri_rs%chunk_centroids, &
                                   cutoff=bs_env%ri_rs%cutoff_radius_v_w)
         CALL dbcsr_scale(mat_Sigma_x_Gamma, -1.0_dp)

         CALL dbcsr_release(matrix_D_ao)

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma(ispin))

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
      CALL dbcsr_release(mat_Sigma_x_Gamma)
      CALL release_square_topology(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy on the imaginary-time grid:
!>          Σ^c_λσ(iτ<0) = -Σ_ll' Φ_λ(r_l) [ (φ G^occ φ^T)_ll' ∘ (Z W^MIC Z^T)_ll' ] Φ_σ(r_l')
!>          Σ^c_λσ(iτ>0) = +Σ_ll' Φ_λ(r_l) [ (φ G^vir φ^T)_ll' ∘ (Z W^MIC Z^T)_ll' ] Φ_σ(r_l')
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
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbcsr_type)                                   :: matrix_G_occ_ao, matrix_G_vir_ao, &
                                                            matrix_W_aux

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP AUXILIARY TOPOLOGY AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(mat_Z_lP, dist_aux_aux, blk_aux, dist_row_aux)

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

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
      ! 2. IMAGINARY TIME LOOP
      ! Σ^c_neg_λσ(iτ) = -φ^T ( (φ G^occ φ^T) ∘ (Z W^MIC Z^T) ) φ
      ! Σ^c_pos_λσ(iτ) =  φ^T ( (φ G^vir φ^T) ∘ (Z W^MIC Z^T) ) φ
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%imag_time_points(i_t)

         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, &
                           blk_aux, blk_aux)
         IF (bs_env%ri_rs%cutoff_radius_g_w > 0.0_dp .AND. ALLOCATED(bs_env%ri_rs%atom_centers)) THEN
            CALL reserve_blocks_within_radius(matrix_W_aux, bs_env%ri_rs%atom_centers, &
                                              bs_env%ri_rs%cutoff_radius_g_w)
            CALL copy_fm_to_dbcsr(fm_W_time(i_t), matrix_W_aux, keep_sparsity=.TRUE.)
         ELSE
            CALL copy_fm_to_dbcsr(fm_W_time(i_t), matrix_W_aux, keep_sparsity=.FALSE.)
         END IF
         CALL dbcsr_filter(matrix_W_aux, bs_env%eps_filter)

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_ao)
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_ao)

            ! Σ^c_neg and Σ^c_pos in a single panel loop: W_pan = Z_panel × W × Z^T built once
            CALL contract_grid_panels_sigma_c(mat_phi=mat_phi_mu_l, mat_Z=mat_Z_lP, &
                                              mat_G_occ_ao=matrix_G_occ_ao, &
                                              mat_G_vir_ao=matrix_G_vir_ao, &
                                              mat_W_aux=matrix_W_aux, &
                                              mat_Sigma_neg=mat_Sigma_neg_tau(i_t, ispin)%matrix, &
                                              mat_Sigma_pos=mat_Sigma_pos_tau(i_t, ispin)%matrix, &
                                              eps=bs_env%eps_filter, &
                                              para_env=bs_env%para_env, &
                                              pan_first=pan_first, pan_last=pan_last, &
                                              keep_sparsity=bs_env%ri_rs%keep_sparsity_rirs, &
                                              centroids=bs_env%ri_rs%chunk_centroids, &
                                              cutoff=bs_env%ri_rs%cutoff_radius_v_w)
            CALL dbcsr_scale(mat_Sigma_neg_tau(i_t, ispin)%matrix, -1.0_dp)

            CALL dbcsr_release(matrix_G_occ_ao)
            CALL dbcsr_release(matrix_G_vir_ao)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I15,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         CALL dbcsr_release(matrix_W_aux)

      END DO ! i_t

      ! -------------------------------------------------------------------------
      ! 3. FINALIZE AND CLEANUP
      ! -------------------------------------------------------------------------
      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      ! fm_W_time and the scratch files are released by the caller: in an evGW0 cycle this
      ! routine is entered once per iteration and both have to survive until it is done.
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL release_square_topology(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief Builds the DBCSR distribution.
!> \param matrix_template ...
!> \param square_dist ...
!> \param blk_sizes ...
!> \param mapped_dist ...
! **************************************************************************************************
   SUBROUTINE setup_square_topology(matrix_template, square_dist, blk_sizes, mapped_dist)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: square_dist
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: blk_sizes, mapped_dist

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_square_topology'

      INTEGER                                            :: handle, i, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist
      TYPE(dbcsr_distribution_type)                      :: dist_template

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_template, distribution=dist_template, col_blk_size=col_blk)
      CALL dbcsr_distribution_get(dist_template, col_dist=col_dist, nprows=nprows)

      blk_sizes => col_blk
      ALLOCATE (mapped_dist(SIZE(blk_sizes)))
      DO i = 1, SIZE(blk_sizes)
         mapped_dist(i) = MOD(i - 1, nprows)
      END DO
      CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                  row_dist=mapped_dist, col_dist=col_dist)

      CALL timestop(handle)

   END SUBROUTINE setup_square_topology

! **************************************************************************************************
!> \brief Releases a distribution created by setup_square_topology.
!> \param dist ...
!> \param mapped_dist ...
! **************************************************************************************************
   SUBROUTINE release_square_topology(dist, mapped_dist)

      TYPE(dbcsr_distribution_type), INTENT(INOUT)       :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: mapped_dist

      CALL dbcsr_distribution_release(dist)
      IF (ASSOCIATED(mapped_dist)) THEN
         DEALLOCATE (mapped_dist)
         NULLIFY (mapped_dist)
      END IF

   END SUBROUTINE release_square_topology

! **************************************************************************************************
!> \brief Σ^c_λσ(iτ) -> Σ^c_nn(ϵ) and the quasi-particle levels of the non-periodic RI-RS path,
!>        ϵ_n^GW = ϵ_n^DFT + Σ^c_nn(ϵ_n^GW) + Σ^x_nn - v^xc_nn.
!> \param bs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER                                            :: handle, ispin, j_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_n, V_xc_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_n_freq, Sigma_c_n_time
      TYPE(cp_fm_type)                                   :: fm_ks, fm_mos, fm_s, fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_ks, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_s, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_mos, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_work, bs_env%fm_s_Gamma%matrix_struct)

      ALLOCATE (V_xc_n(bs_env%n_ao), Sigma_x_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_n_freq(bs_env%n_ao, bs_env%num_time_freq_points, 2))

      DO ispin = 1, bs_env%n_spin

         ! 1. Roothaan-Hall H^KS_µν C_νn = S_µν C_νn ϵ_n
         CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(ispin), fm_ks)
         CALL cp_fm_to_fm(bs_env%fm_s_Gamma, fm_s)
         CALL cp_fm_geeig(fm_ks, fm_s, fm_mos, bs_env%eigenval_scf(:, 1, ispin), fm_work)

         ! 2. v^xc_µν -> v^xc_nn  and  Σ^x_µν -> Σ^x_nn
         CALL to_Gamma_and_mo_real(V_xc_n, bs_env%fm_V_xc_Gamma(ispin), fm_mos)
         CALL to_Gamma_and_mo_real(Sigma_x_n, fm_Sigma_x_Gamma(ispin), fm_mos)

         ! 3. Σ^c_µν(+/-i|τ_j|) -> Σ^c_nn(+/-i|τ_j|)
         DO j_t = 1, bs_env%num_time_freq_points
            CALL to_Gamma_and_mo_real(Sigma_c_n_time(:, j_t, 1), &
                                      fm_Sigma_c_Gamma_time(j_t, 1, ispin), fm_mos)
            CALL to_Gamma_and_mo_real(Sigma_c_n_time(:, j_t, 2), &
                                      fm_Sigma_c_Gamma_time(j_t, 2, ispin), fm_mos)
         END DO

         ! 4. Σ^c_nn(iτ) -> Σ^c_nn(iω)
         CALL time_to_freq(bs_env, Sigma_c_n_time, Sigma_c_n_freq, ispin)

         ! 5. Analytic continuation Σ^c_nn(iω) -> Σ^c_nn(ϵ) and the QP levels
         CALL analyt_conti_and_print(bs_env, Sigma_c_n_freq, Sigma_x_n, V_xc_n, &
                                     bs_env%eigenval_scf(:, 1, ispin), 1, ispin)

      END DO ! ispin

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      IF (bs_env%gw_flavour == G0W0) CALL cp_fm_release(fm_Sigma_x_Gamma)
      CALL cp_fm_release(fm_Sigma_c_Gamma_time)

      CALL cp_fm_release(fm_ks)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_mos)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief AO -> MO transform of a Γ-point matrix
!> \param array_n ...
!> \param fm_Gamma ...
!> \param fm_mos ...
! **************************************************************************************************
   SUBROUTINE to_Gamma_and_mo_real(array_n, fm_Gamma, fm_mos)

      REAL(KIND=dp), DIMENSION(:)                        :: array_n
      TYPE(cp_fm_type)                                   :: fm_Gamma, fm_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'to_Gamma_and_mo_real'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fm_mo

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mo, fm_Gamma%matrix_struct)

      ! A_nn' = Σ_μν C_μn A_μν C_νn'
      CALL fm_contract_ABA(fm_mos, fm_Gamma, fm_mo)

      CALL cp_fm_get_diag(fm_mo, array_n)

      CALL cp_fm_release(fm_mo)

      CALL timestop(handle)

   END SUBROUTINE to_Gamma_and_mo_real

! **************************************************************************************************
!> \brief Evaluates the atom-local AO functions Φ_μ(r_l) on the RI-RS grid points.
!> \param phi_val Values Φ_μ(r_l) for the requested atom
!> \param ri_rs_grid Cartesian RI-RS grid points r_l
!> \param npts Number of RI-RS grid points
!> \param iatom Atom whose AO functions are evaluated
!> \param particle_set Atomic positions and kind assignments
!> \param qs_kind_set Atomic-kind basis information
!> \param cell Simulation cell and boundary conditions
!> \param r2_threshold Squared radial cutoff for neglecting AO values
! **************************************************************************************************
   SUBROUTINE fill_phi_for_atom(phi_val, ri_rs_grid, npts, iatom, &
                                particle_set, qs_kind_set, cell, r2_threshold)

      REAL(KIND=dp), INTENT(INOUT)                       :: phi_val(:, :)
      INTEGER, INTENT(IN)                                :: npts
      REAL(KIND=dp), INTENT(IN)                          :: ri_rs_grid(3, npts)
      INTEGER, INTENT(IN)                                :: iatom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: r2_threshold

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fill_phi_for_atom'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL evaluate_ao_on_points(phi_val, ri_rs_grid, iatom, particle_set, qs_kind_set, cell, &
                                 cutoff_squared=r2_threshold)

      CALL timestop(handle)

   END SUBROUTINE fill_phi_for_atom

! **************************************************************************************************
!> \brief Counts the AO functions whose radial support intersects an RI fitting sphere.
!> \param bs_env GW calculation environment and atom-to-AO offsets
!> \param particle_set Atomic positions
!> \param atom_P Atom at the center of the RI fitting sphere
!> \param cutoff_ri Radius of the RI fitting sphere
!> \param n_ao_used Number of intersecting AO functions
! **************************************************************************************************
   SUBROUTINE get_n_ao_in_sphere(bs_env, particle_set, atom_P, cutoff_ri, n_ao_used)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: atom_P
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_ri
      INTEGER, INTENT(OUT)                               :: n_ao_used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_n_ao_in_sphere'

      INTEGER                                            :: handle, ri_atom

      CALL timeset(routineN, handle)

      n_ao_used = 0
      DO ri_atom = 1, bs_env%n_atom
         IF (NORM2(particle_set(ri_atom)%r(:) - particle_set(atom_P)%r(:)) > &
             bs_env%ri_rs%radius_ao_per_atom(ri_atom) + cutoff_ri) CYCLE
         n_ao_used = n_ao_used + bs_env%i_ao_end_from_atom(ri_atom) - &
                     bs_env%i_ao_start_from_atom(ri_atom) + 1
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_n_ao_in_sphere

! **************************************************************************************************
!> \brief Prints the percentage of non-zero elements in a distributed RI-RS matrix.
!> \param matrix Distributed matrix whose occupation is reported
!> \param label Mathematical matrix label used in the output
!> \param bs_env GW calculation environment and output unit
!> \param suffix Optional text appended to the matrix label
! **************************************************************************************************
   SUBROUTINE print_matrix_occupation(matrix, label, bs_env, suffix)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: suffix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_matrix_occupation'

      CHARACTER(LEN=32)                                  :: output_format
      CHARACTER(LEN=max_line_length)                     :: msg, output_label
      INTEGER                                            :: handle, i, unicode_shift
      REAL(KIND=dp)                                      :: frac_2p31, max_loc, occ

      CALL timeset(routineN, handle)

      occ = dbcsr_get_occupation(matrix)
      max_loc = REAL(dbcsr_get_data_size(matrix), dp)
      CALL bs_env%para_env%max(max_loc)

      IF (bs_env%unit_nr > 0) THEN
         frac_2p31 = max_loc/REAL(bs_env%dbcsr_msg_elem_limit, dp)
         output_label = 'Percentage of non-zero matrix elements in '//TRIM(label)
         IF (PRESENT(suffix)) output_label = TRIM(output_label)//TRIM(suffix)
         ! Fortran counts UTF-8 bytes, whereas the terminal displays each Greek letter in one
         ! column. Shift the absolute output tab once for every continuation byte.
         unicode_shift = 0
         DO i = 1, LEN_TRIM(output_label)
            IF (IAND(IACHAR(output_label(i:i)), 192) == 128) unicode_shift = unicode_shift + 1
         END DO
         WRITE (output_format, '(A,I0,A)') '(T2,A,T', 72 + unicode_shift, ',F7.2,A)'
         WRITE (bs_env%unit_nr, output_format) TRIM(output_label), occ*100.0_dp, ' %'
         IF (frac_2p31 > 0.5_dp) THEN
            WRITE (msg, '(3A,F0.2,A)') &
               "The largest per-rank message of ", TRIM(label), " reaches ", frac_2p31, &
               " of the 32-bit limit that DBCSR uses for its message length. Beyond it the "// &
               "length overflows and multiply_cannon fails. Reduce the per-rank block size, "// &
               "for instance with more MPI ranks or a larger N_PANELS."
            CPWARN(TRIM(msg))
         END IF
         CALL m_flush(bs_env%unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_matrix_occupation

! **************************************************************************************************
!> \brief Returns the available memory per MPI process on the most constrained node.
!> \param bs_env GW calculation environment and MPI communicator
!> \param mem_avail_GB Available memory per MPI process in GB
! **************************************************************************************************
   SUBROUTINE ri_rs_mem_avail_per_proc_GB(bs_env, mem_avail_GB)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: mem_avail_GB

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_rs_mem_avail_per_proc_GB'

      CHARACTER(LEN=default_string_length)               :: hostname
      INTEGER                                            :: handle, host_hash, ic, n_procs, &
                                                            ranks_per_node
      INTEGER(KIND=int_8)                                :: h8, mem_buffers, mem_cached, mem_free, &
                                                            mem_likely_free, mem_slab, &
                                                            mem_sreclaimable, mem_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: all_host_hashes

      CALL timeset(routineN, handle)

      n_procs = bs_env%para_env%num_pe
      CALL m_hostnm(hostname)
      h8 = 0_int_8
      DO ic = 1, LEN_TRIM(hostname)
         h8 = MOD(h8*127_int_8 + INT(ICHAR(hostname(ic:ic)), int_8), 2147483647_int_8)
      END DO
      host_hash = INT(h8)
      ALLOCATE (all_host_hashes(n_procs))
      CALL bs_env%para_env%allgather(host_hash, all_host_hashes)
      ranks_per_node = MAX(COUNT(all_host_hashes == host_hash), 1)
      DEALLOCATE (all_host_hashes)

      CALL m_memory_details(MemTotal=mem_total, MemFree=mem_free, Buffers=mem_buffers, &
                            Cached=mem_cached, Slab=mem_slab, SReclaimable=mem_sreclaimable, &
                            MemLikelyFree=mem_likely_free)
      mem_avail_GB = REAL(mem_likely_free, dp)*1.0E-9_dp/REAL(ranks_per_node, dp)
      CALL bs_env%para_env%min(mem_avail_GB)

      CALL timestop(handle)

   END SUBROUTINE ri_rs_mem_avail_per_proc_GB

END MODULE gw_non_periodic_ri_rs
