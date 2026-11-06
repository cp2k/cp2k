!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fully periodic small-cell GW using RI-RS and explicit full k-point sampling.
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE gw_ri_rs_small_cell_full_kp
   USE gw_ri_rs_small_cell_full_kp_types,ONLY: gwkp_type
   USE gw_ri_rs_small_cell_full_kp_utils,ONLY: compute_Sigma_x,&
                                              compute_W0,&
                                              compute_g_and_chi,&
                                              compute_sigma_c_from_w_pq,&
                                              prepare_ri_rs,&
                                              release_gw_work
   USE gw_utils,                        ONLY: analyt_conti_and_print,&
                                              de_init_bs_env
   USE input_constants,                 ONLY: evGW0
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: band_edges_type,&
                                              post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: get_VBM_CBM_bandgaps,&
                                              get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: qs_environment_type
   USE time_frequency_grids,            ONLY: transform_time_to_frequency
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp'

   PUBLIC :: gw_calc_ri_rs_small_cell_full_kp

CONTAINS

! **************************************************************************************************
!> \brief Periodic RI-RS GW from the converged SCF orbitals and eigenvalues.
!>
!> l,l′ label real-space points r_l,r_l′; μ labels AOs, n bands, σ spins, P,Q RI functions.
!> k samples the SCF mesh, q transfers momentum, and R,T are lattice vectors.
!> All phases use Cartesian vectors. ψ_nkσ(r_l) is a sampled Bloch orbital;
!> C_μn(k,σ) is its AO coefficient. σ sums over the two physical spin channels.
!>
!>  (1) ϕ_μ(k,r_l) = Σ_R exp(ik·R) ϕ_μ(r_l-R),
!>      ψ_nkσ(r_l) = Σ_μ ϕ_μ(k,r_l) C_μn(k,σ).
!>
!>  (2) ρ_l(r) = Σ_μR ϕ_μ(r_l-R) ϕ_μ(r-R),
!>      D_ll′(T) = ρ_l(r_l′+T)²,
!>      d_lP(T) = ∫dr ∫dr′ ρ_l(r)² m(r-r′) P(r′-T),
!>      D_ll′(q) = Σ_T exp(iq·T) D_ll′(T),
!>      d_lP(q) = Σ_T exp(iq·T) d_lP(T),
!>      S_l(q) = D_ll(q)^(-1/2),
!>      D̂_ll′(q) = S_l(q) D_ll′(q) S_l′(q),  d̂_lP(q) = S_l(q) d_lP(q),
!>      [D̂(q) + ηI] z_P(q) = d̂_P(q),  Z_lP(q) = S_l(q) z_lP(q).
!>      P(r) is an RI basis function, m(r-r′) the selected RI fitting kernel,
!>      and η = RI_RS/TIKHONOV.
!>
!>  (3) G^occ/vir_ll′(k,iτ;σ) = Σ_n∈occ/vir ψ_nkσ(r_l) ψ_nkσ(r_l′)*
!>                                          × exp(-|ε_nkσ-ε_Fσ| τ),  τ≥0.
!>
!>  (4) G^occ/vir_ll′(R,iτ;σ) = Re[Σ_k w_k exp(-ik·R) G^occ/vir_ll′(k,iτ;σ)].
!>
!>  (5) χ_ll′(R,iτ) = Σ_σ G^occ_ll′(R,iτ;σ) G^vir_ll′(R,iτ;σ),
!>      χ_ll′(q,iτ) = Σ_R exp(iq·R) χ_ll′(R,iτ).
!>
!>  (6) χ_ll′(q,iω_j) = Σ_t C_jt χ_ll′(q,iτ_t),
!>      χ_PQ(q,iω_j) = Σ_ll′ Z_lP(q)* χ_ll′(q,iω_j) Z_l′Q(q).
!>      C_jt supplies the imaginary-time Fourier transformation by cosine quadrature.
!>
!>  (7) K(q) = M(q)⁻¹ V(q)½,
!>      ε_PQ(q,iω) = δ_PQ + Σ_P′Q′ K_P′P(q)* χ_P′Q′(q,iω) K_Q′Q(q),
!>      Wᶜ_PQ(q,iω) = Σ_P′Q′ K_PP′(q) [ε⁻¹_P′Q′(q,iω)-δ_P′Q′] K_QQ′(q)*.
!>      M is the RI fitting metric and V the Coulomb matrix.
!>
!>  (8) Wᶜ_ll′(q,iτ) = Σ_PQ Z_lP(q) Wᶜ_PQ(q,iτ) Z_l′Q(q)*,
!>      Wᶜ_ll′(R,iτ) = Re[Σ_q w_q exp(-iq·R) Wᶜ_ll′(q,iτ)] + ΔWᶜ_Γ,ll′(R,iτ).
!>
!>  (9) V_ll′(q) = Σ_PQ Z_lP(q) [K(q)K(q)†]_PQ Z_l′Q(q)*,
!>      V_ll′(R) = Re[Σ_q w_q exp(-iq·R) V_ll′(q)] + ΔV_Γ,ll′(R),
!>      Σˣ_ll′(R;σ) = -G^occ_ll′(R,0;σ) V_ll′(R).
!>
!> (10) Σˣ_ll′(k;σ) = Σ_R exp(ik·R) Σˣ_ll′(R;σ),
!>      Σˣ_nkσ = Σ_ll′ ψ_nkσ(r_l)* Σˣ_ll′(k;σ) ψ_nkσ(r_l′).
!>
!> (11) Σᶜ_ll′(R,-iτ;σ) = -G^occ_ll′(R,iτ;σ) Wᶜ_ll′(R,iτ),
!>      Σᶜ_ll′(R,+iτ;σ) = +G^vir_ll′(R,iτ;σ) Wᶜ_ll′(R,iτ).
!>
!> (12) Σᶜ_ll′(k,±iτ;σ) = Σ_R exp(ik·R) Σᶜ_ll′(R,±iτ;σ),
!>      Σᶜ_nkσ(±iτ) = Σ_ll′ ψ_nkσ(r_l)* Σᶜ_ll′(k,±iτ;σ) ψ_nkσ(r_l′),
!>      Σᶜ_nkσ(iω) = ∫_0^∞ dτ [exp(iωτ)Σᶜ_nkσ(+iτ)+exp(-iωτ)Σᶜ_nkσ(-iτ)].
!>      Analytic continuation gives Σᶜ at real energies for the quasiparticle equation.
!>
!> (13) ε^GW_nkσ = ε^DFT_nkσ + Σˣ_nkσ + Re Σᶜ_nkσ(ε^GW_nkσ) - vˣᶜ_nkσ.
!>      evGW₀ updates the energies entering G while retaining Wᶜ from the first cycle.
!>
!> Γ-point corrections are included as ΔWᶜ_Γ and ΔV_Γ; their construction is not shown here:
!> compute_G_and_chi prepares ΔWᶜ_Γ through gamma_point_correction.
!> compute_Sigma_c_from_W_PQ adds it to Wᶜ before forming Σᶜ.
!> compute_Sigma_x_from_V computes and applies ΔV_Γ via singular_head_cell_factors.
!>
!> \param qs_env Quickstep environment containing the converged SCF state.
!> \param bs_env Band-structure settings, orbitals, and quasiparticle results.
! **************************************************************************************************
   SUBROUTINE gw_calc_ri_rs_small_cell_full_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_small_cell_full_kp'

      INTEGER                                            :: handle
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      CALL timeset(routineN, handle)

      gwkp_env => bs_env%gwkp_env

      ! Sample AO/Bloch orbitals and fit the auxiliary functions: Eqs. (1)-(2).
      CALL prepare_ri_rs(bs_env, gwkp_env)

      ! G^occ/vir_ll′(R,iτ), χ_ll′(R,iτ), and χ_PQ(q,iω): Eqs. (3)-(6).
      CALL compute_G_and_chi(qs_env, bs_env, gwkp_env)

      ! Screen the Coulomb interaction and transform Wᶜ to imaginary time: Eqs. (7)-(8).
      CALL compute_W0(qs_env, bs_env, gwkp_env)

      ! Σˣ_ll′(R) and its projection onto target Bloch states: Eqs. (9)-(10).
      CALL compute_Sigma_x(bs_env, gwkp_env)

      ! Σᶜ_ll′(R,iτ), band projection, and quasiparticle equation: Eqs. (11)-(13).
      CALL compute_Sigma_c_and_QP(bs_env, gwkp_env)

      CALL release_gw_work(bs_env)
      NULLIFY (gwkp_env)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_small_cell_full_kp

! **************************************************************************************************
!> \brief Iterate Σᶜ and the QP equation, Eqs. (11)-(13), with fixed W₀.
!> G₀W₀ uses one cycle; evGW₀ feeds the QP energies back into G on the SCF mesh.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Temporary periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_and_QP(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_and_QP'

      INTEGER                                            :: handle, i_iter, ispin, n_mo, nkp
      REAL(dp)                                           :: cbm, e_fermi_scf(2), max_delta, vbm
      REAL(dp), ALLOCATABLE                              :: eigenvalues_G_next(:, :, :)

      CALL timeset(routineN, handle)

      nkp = bs_env%nkp_scf_desymm
      n_mo = bs_env%n_mo_retained
      e_fermi_scf = bs_env%e_fermi
      IF (bs_env%gw_flavour == evGW0) ALLOCATE (eigenvalues_G_next(n_mo, nkp, bs_env%n_spin))

      ! G₀W₀ takes one cycle. evGW₀ updates G on the SCF mesh while W₀ stays fixed.
      DO i_iter = 1, bs_env%evgw0_n_iter
         CALL compute_Sigma_c_from_W_PQ(bs_env, gwkp_env)
         CALL compute_QP_energies(bs_env, gwkp_env%sigma_x_bands, &
                                  gwkp_env%sigma_c_bands, i_iter)
         IF (bs_env%evgw0_n_iter == 1) EXIT

         IF (i_iter == 1) CALL print_first_cycle_band_edges(bs_env)
         CALL eigenvalues_for_evGW0_G(bs_env, eigenvalues_G_next)
         max_delta = MAXVAL(ABS(eigenvalues_G_next - &
                                bs_env%eigenval_evGW0(1:n_mo, 1:nkp, 1:bs_env%n_spin)))
         IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,I3,A,I3,A,T66,F12.4,A)') &
            'evGW0 cycle ', i_iter, ' / ', bs_env%evgw0_n_iter, &
            ' Max G eigenvalue change:', max_delta*evolt, ' eV'
         IF (max_delta < bs_env%evgw0_eps_iter) EXIT
         IF (i_iter == bs_env%evgw0_n_iter) THEN
            IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A)') &
               'evGW0 reached MAX_ITER before eigenvalue convergence.'
            EXIT
         END IF

         bs_env%eigenval_evGW0(1:n_mo, 1:nkp, 1:bs_env%n_spin) = eigenvalues_G_next
         gwkp_env%eigenvalues_G(1:n_mo, 1:nkp, 1:bs_env%n_spin) = eigenvalues_G_next
         DO ispin = 1, bs_env%n_spin
            vbm = MAXVAL(eigenvalues_G_next(bs_env%n_occ(ispin), :, ispin))
            cbm = MINVAL(eigenvalues_G_next(bs_env%n_occ(ispin) + 1, :, ispin))
            bs_env%e_fermi(ispin) = 0.5_dp*(vbm + cbm)
         END DO
      END DO

      CALL get_all_VBM_CBM_bandgaps(bs_env)
      IF (bs_env%gw_flavour == evGW0) bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW
      ! G used an updated energy origin; subsequent post-SCF stages retain the SCF reference.
      bs_env%e_fermi = e_fermi_scf

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_and_QP

! **************************************************************************************************
!> \brief Print scalar G0W0 extrema over the SCF mesh and requested band path.
!> The first QP solution is reported before updating G; W remains W0 throughout.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_first_cycle_band_edges(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_first_cycle_band_edges'

      INTEGER                                            :: handle, u
      TYPE(band_edges_type)                              :: band_edges

      CALL timeset(routineN, handle)

      CALL get_VBM_CBM_bandgaps(band_edges, bs_env%eigenval_GW, bs_env)
      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') 'G0W0 from the first evGW0 cycle (without SOC):'
         WRITE (u, '(T2,A,T61,F20.3)') 'G0W0 valence band maximum (eV):', band_edges%VBM*evolt
         WRITE (u, '(T2,A,T61,F20.3)') 'G0W0 conduction band minimum (eV):', band_edges%CBM*evolt
         WRITE (u, '(T2,A,T61,F20.3)') 'G0W0 indirect band gap (eV):', band_edges%IDBG*evolt
         WRITE (u, '(T2,A,T61,F20.3)') 'G0W0 direct band gap (eV):', band_edges%DBG*evolt
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_first_cycle_band_edges

! **************************************************************************************************
!> \brief Build the eigenvalues used by the next Green function from the current QP solution.
!> The 5 eV windows use the DFT VBM/CBM on the full SCF mesh, separately for each spin.
!> Outside them, a mean shift of the window states is applied at each k point and spin.
!> If a window is empty at one k point, its nearest band supplies the shift.
!> This approximation never changes the quasiparticle energies or the final band structure.
!> \param bs_env ...
!> \param eigenvalues_G_next ...
! **************************************************************************************************
   SUBROUTINE eigenvalues_for_evGW0_G(bs_env, eigenvalues_G_next)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: eigenvalues_G_next(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eigenvalues_for_evGW0_G'
      REAL(KIND=dp), PARAMETER                           :: window_ev = 5.0_dp

      INTEGER                                            :: handle, iband, ikp, ispin, n_mo, &
                                                            n_window, nkp, nocc
      REAL(KIND=dp)                                      :: cbm_dft, dft_energy, lower_limit, &
                                                            mean_shift, qp_energy, upper_limit, &
                                                            vbm_dft

      CALL timeset(routineN, handle)

      nkp = bs_env%nkp_scf_desymm
      n_mo = bs_env%n_mo_retained
      CPASSERT(bs_env%nkp_bs_and_DOS >= nkp)
      CPASSERT(n_mo > 0 .AND. n_mo <= bs_env%n_ao)
      CPASSERT(ALL(bs_env%n_occ(1:bs_env%n_spin) > 0))
      CPASSERT(ALL(bs_env%n_occ(1:bs_env%n_spin) < n_mo))
      DO ispin = 1, bs_env%n_spin
         nocc = bs_env%n_occ(ispin)
         vbm_dft = MAXVAL(bs_env%eigenval_scf(nocc, 1:nkp, ispin))
         cbm_dft = MINVAL(bs_env%eigenval_scf(nocc + 1, 1:nkp, ispin))
         lower_limit = vbm_dft - window_ev/evolt
         upper_limit = cbm_dft + window_ev/evolt
         DO ikp = 1, nkp
            mean_shift = 0.0_dp
            n_window = 0
            DO iband = 1, nocc
               dft_energy = bs_env%eigenval_scf(iband, ikp, ispin)
               IF (dft_energy < lower_limit) CYCLE
               qp_energy = bs_env%eigenval_GW(iband, ikp, ispin)
               eigenvalues_G_next(iband, ikp, ispin) = qp_energy
               mean_shift = mean_shift + qp_energy - dft_energy
               n_window = n_window + 1
            END DO
            IF (n_window > 0) THEN
               mean_shift = mean_shift/REAL(n_window, KIND=dp)
            ELSE
               mean_shift = bs_env%eigenval_GW(nocc, ikp, ispin) - &
                            bs_env%eigenval_scf(nocc, ikp, ispin)
            END IF
            DO iband = 1, nocc
               dft_energy = bs_env%eigenval_scf(iband, ikp, ispin)
               IF (dft_energy >= lower_limit) CYCLE
               eigenvalues_G_next(iband, ikp, ispin) = dft_energy + mean_shift
            END DO

            mean_shift = 0.0_dp
            n_window = 0
            DO iband = nocc + 1, n_mo
               dft_energy = bs_env%eigenval_scf(iband, ikp, ispin)
               IF (dft_energy > upper_limit) CYCLE
               qp_energy = bs_env%eigenval_GW(iband, ikp, ispin)
               eigenvalues_G_next(iband, ikp, ispin) = qp_energy
               mean_shift = mean_shift + qp_energy - dft_energy
               n_window = n_window + 1
            END DO
            IF (n_window > 0) THEN
               mean_shift = mean_shift/REAL(n_window, KIND=dp)
            ELSE
               mean_shift = bs_env%eigenval_GW(nocc + 1, ikp, ispin) - &
                            bs_env%eigenval_scf(nocc + 1, ikp, ispin)
            END IF
            DO iband = nocc + 1, n_mo
               dft_energy = bs_env%eigenval_scf(iband, ikp, ispin)
               IF (dft_energy <= upper_limit) CYCLE
               eigenvalues_G_next(iband, ikp, ispin) = dft_energy + mean_shift
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE eigenvalues_for_evGW0_G

! **************************************************************************************************
!> \brief Transform band self-energies and solve the quasiparticle equations.
!> \param bs_env ...
!> \param sigma_x_bands ...
!> \param sigma_c_bands ...
!> \param i_iter Local GW cycle used in the spectrum output.
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env, sigma_x_bands, sigma_c_bands, i_iter)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: sigma_x_bands(:, :, :), &
                                                            sigma_c_bands(:, :, :, :, :)
      INTEGER, INTENT(IN)                                :: i_iter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: Sigma_c_ikp_n_freq(:, :)
      INTEGER                                            :: handle, ikp, ispin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_time

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(sigma_x_bands, 1) == bs_env%n_ao)
      CPASSERT(SIZE(sigma_x_bands, 2) == bs_env%nkp_bs_and_DOS)
      CPASSERT(SIZE(sigma_x_bands, 3) == bs_env%n_spin)
      CPASSERT(SIZE(sigma_c_bands, 1) == bs_env%n_ao)
      CPASSERT(SIZE(sigma_c_bands, 2) == bs_env%nkp_bs_and_DOS)
      CPASSERT(SIZE(sigma_c_bands, 3) == bs_env%num_time_freq_points)
      CPASSERT(SIZE(sigma_c_bands, 4) == 2)
      CPASSERT(SIZE(sigma_c_bands, 5) == bs_env%n_spin)

      ALLOCATE (Sigma_x_ikp_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_ikp_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_ikp_n_freq(bs_env%n_ao, bs_env%num_freq_points_fit))

      DO ispin = 1, bs_env%n_spin

         DO ikp = 1, bs_env%nkp_bs_and_DOS

            ! 1. Σ^x_nn(k) and Σ^c_nn(k,±i|τ_j|), projected in compute_Sigma_x / compute_Sigma_c
            Sigma_x_ikp_n(:) = sigma_x_bands(:, ikp, ispin)
            Sigma_c_ikp_n_time(:, :, :) = sigma_c_bands(:, ikp, :, :, ispin)

            ! 2. Σ^c_nn(k_i,iω) = ∫ from -∞ to ∞ dτ e^-iωτ Σ^c_nn(k_i,iτ)
            CALL transform_time_to_frequency(bs_env%time_frequency_grid, Sigma_c_ikp_n_time, Sigma_c_ikp_n_freq)
            Sigma_c_ikp_n_freq(1:bs_env%n_occ(ispin), :) = CONJG(Sigma_c_ikp_n_freq(1:bs_env%n_occ(ispin), :))

            ! 3. Analytic continuation Σ^c_nn(k_i,iω) -> Σ^c_nn(k_i,ϵ) and
            !    ϵ_nk_i^GW = ϵ_nk_i^DFT + Σ^c_nn(k_i,ϵ) + Σ^x_nn(k_i) - v^xc_nn(k_i)
            CALL analyt_conti_and_print(bs_env, &
                                        Sigma_c_ikp_n_freq, &
                                        Sigma_x_ikp_n, &
                                        bs_env%v_xc_n(:, ikp, ispin), &
                                        bs_env%eigenval_scf(:, ikp, ispin), ikp, ispin, i_iter)

         END DO ! ikp

      END DO ! ispin

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

END MODULE gw_ri_rs_small_cell_full_kp
