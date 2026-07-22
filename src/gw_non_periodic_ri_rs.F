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
   USE atomic_kind_types, ONLY: atomic_kind_type, &
                                get_atomic_kind_set
   USE basis_set_types, ONLY: gto_basis_set_type
   USE cell_types, ONLY: cell_type, &
                         pbc
   USE constants_operator, ONLY: operator_coulomb
   USE cp_blacs_env, ONLY: cp_blacs_env_create, &
                           cp_blacs_env_release, &
                           cp_blacs_env_type
   USE cp_dbcsr_api, ONLY: &
      dbcsr_add, dbcsr_binary_read, dbcsr_binary_write, dbcsr_copy, dbcsr_create, &
      dbcsr_deallocate_matrix, dbcsr_distribution_get, dbcsr_distribution_new, &
      dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, &
      dbcsr_get_block_p, dbcsr_filter, dbcsr_get_data_size, dbcsr_get_info, &
      dbcsr_get_occupation, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
      dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
      dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
      dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib, ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                  copy_fm_to_dbcsr, &
                                  dbcsr_deallocate_matrix_set, &
                                  max_elements_per_block
   USE cp_files, ONLY: close_file, &
                       open_file
   USE cp_fm_basic_linalg, ONLY: cp_fm_scale_and_add, &
                                 cp_fm_uplo_to_full
   USE cp_fm_cholesky, ONLY: cp_fm_cholesky_decompose, &
                             cp_fm_cholesky_invert, &
                             cp_fm_cholesky_solve
   USE cp_fm_diag, ONLY: cp_fm_power
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_get_info, &
                          cp_fm_get_submatrix, &
                          cp_fm_release, &
                          cp_fm_set_all, &
                          cp_fm_set_submatrix, &
                          cp_fm_to_fm, &
                          cp_fm_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_type
   USE cp_output_handling, ONLY: cp_p_file, &
                                 cp_print_key_should_output
   USE gw_integrals, ONLY: build_3c_integral_block_ctx, &
                           gw_3c_ctx_create, &
                           gw_3c_ctx_release, &
                           gw_3c_ctx_type, &
                           gw_3c_ws_create, &
                           gw_3c_ws_release, &
                           gw_3c_ws_type
   USE gw_large_cell_gamma, ONLY: &
      Fourier_transform_w_to_t, G_occ_vir, compute_QP_energies, compute_fm_chi_Gamma_freq, &
      create_fm_W_MIC_time, delete_unnecessary_files, fill_fm_Sigma_c_Gamma_time, fm_write, &
      multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_utils, ONLY: de_init_bs_env
   USE input_constants, ONLY: rtp_method_bse
   USE input_section_types, ONLY: section_vals_type
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp, int_4, int_8
   USE kpoint_coulomb_2c, ONLY: build_2c_coulomb_matrix_kp
   USE machine, ONLY: m_flush, m_hostnm, &
                      m_memory_details, m_walltime
   USE message_passing, ONLY: mp_para_env_type
   USE mp2_ri_2c, ONLY: RI_2c_integral_mat
   USE orbital_pointers, ONLY: indco, &
                               ncoset
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE particle_types, ONLY: particle_type
   USE physcon, ONLY: angstrom
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_kind_types, ONLY: get_qs_kind, &
                            qs_kind_type
   USE util, ONLY: sort

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_non_periodic_ri_rs'

   ! DBCSR ships each ranks matrix panel in one MPI message whose length is stored in a
   ! 32-bit default INTEGER (dbcsr_mpiwrap.F: msglen = SIZE(buffer)). Exceeding this ceiling
   ! overflows the count to a negative/garbage value and multiply_cannon segfaults. This is
   ! that hard limit in elements (= HUGE(int_4) = 2^31-1), the same bound DBCSR uses
   ! internally as mp_max_memory_size. Panel sizing keeps every per-rank message below it.

   INTEGER(KIND=int_8), PARAMETER, PRIVATE :: dbcsr_msg_elem_limit = INT(HUGE(0_int_4), int_8)

   PUBLIC :: gw_calc_non_periodic_ri_rs, ri_rs_grid_assembler, &
             get_basis_offsets, precompute_ri_rs_radii, solve_D_lp_distributed

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
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      ! ========================================================================
      ! 0. Precompute AO and RI radii
      !    Per-atom cutoff radii from the most diffuse Gaussian primitives of
      !    the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets:
      !      α_min,ao = min { ζ_ao | ζ_ao > 10⁻³ },  α_min,ri analogous
      !      r_ao = sqrt( -ln(ε) / α_min,ao )   (radius_ao_per_atom)
      !      r_ri = sqrt( -ln(ε) / α_min,ri )   (radius_ri_per_atom)
      ! ========================================================================
      CALL precompute_ri_rs_radii(qs_env, bs_env)

      ! ========================================================================
      ! 1. Grid generation for RI-RS
      !    Modified Lebedev atomic grids (Duchemin & Blase), one per atom,
      !    concatenated into a flat global list:  r_l = R_A + r_l^(A)
      ! ========================================================================
      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)

      ! ========================================================================
      ! 2a. Atomic basis evaluation on the grid (grid x AO matrix)
      !     Φ_μl = Φ_μ(r_l)   (mat_phi_mu_l)
      ! ========================================================================
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)

      ! ========================================================================
      ! 2b. Print the memory estimate for the RI-RS calculation
      ! ========================================================================
      CALL print_ri_rs_memory_estimate(qs_env, bs_env)

      ! ========================================================================
      ! 3. RI-RS fitting coefficients Z_lP (grid x RI matrix)
      !    Per-atom regularized solve, restricted to grid points r_l within a
      !    cutoff distance of atom P:
      !      a. D_ll' = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]²
      !      b. D_lP  = Σ_μν Φ_μ(r_l) Φ_ν(r_l) (μν|P)
      !      c. Jacobi conditioning with d_l = 1/sqrt(D_ll):
      !         D'_ll' = d_l D_ll' d_l' + λδ_ll' ,  D'_lP = d_l D_lP
      !      d. Solve Σ_l' D'_ll' Z'_l'P = D'_lP
      !      e. Rescale Z_lP = d_l Z'_lP                        (mat_Z_lP)
      ! ========================================================================
      CALL compute_coeff_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                              bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! ========================================================================
      ! 4. Polarizability matrix χ on the imaginary-time grid
      !      G^occ_µλ(i|τ|)  = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^vir_µλ(i|τ|)  = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^occ_ll'(i|τ|) = Σ_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')  (G^vir analogous)
      !      χ_ll'(iτ)       = G^occ_ll'(i|τ|) ∘ G^vir_ll'(i|τ|) (element-wise)
      !      χ_PQ(iτ)        = Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q
      ! ========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! ========================================================================
      ! 5. Screened Coulomb interaction W (RI basis)
      !      χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W_PQ(iτ)
      ! ========================================================================
      CALL compute_W(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      ! ========================================================================
      ! 6. Exact-exchange self-energy Σ^x
      !      D_µν          = Σ_n^occ C_µn C_νn                (density matrix)
      !      D_ll'         = Σ_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
      !      V^tr_ll'      = Σ_PQ Z_lP V^tr_PQ Z_l'Q          (truncated Coulomb)
      !      Σ^x_ll'       = D_ll' ∘ V^tr_ll'
      !      Σ^x_λσ(k=0)   = -Σ_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
      ! ========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

      ! ========================================================================
      ! 7. Correlation self-energy Σ^c on the imaginary-time grid
      !      W_ll'(iτ)     =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
      !      Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
      !      Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
      !      Σ^c_λσ(iτ)    =  Σ_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l')
      ! ========================================================================
      CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

      ! ========================================================================
      ! 8. Quasiparticle energies (analytic continuation iτ -> iω -> real ϵ)
      !      Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
      !      ϵ_n^GW = ϵ_n^DFT + Σ^c_nn(ϵ_n^GW) + Σ^x_nn - v^xc_nn
      ! ========================================================================
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_non_periodic_ri_rs

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
      natom = SIZE(particle_set)

      eps = bs_env%eps_filter

      ALLOCATE (alpha_min_ao_kind(nkind), alpha_min_ri_kind(nkind))
      alpha_min_ao_kind = HUGE(1.0_dp)
      alpha_min_ri_kind = HUGE(1.0_dp)

      DO ikind = 1, nkind
         zet_ao => bs_env%basis_set_AO(ikind)%gto_basis_set%zet
         zet_ri => bs_env%basis_set_RI(ikind)%gto_basis_set%zet

         DO i = 1, SIZE(zet_ao, 1)
            DO j = 1, SIZE(zet_ao, 2)
               IF (zet_ao(i, j) > 1.0E-3_dp) then
                  alpha_min_ao_kind(ikind) = MIN(alpha_min_ao_kind(ikind), zet_ao(i, j))
               END IF
            END DO
         END DO
         DO i = 1, SIZE(zet_ri, 1)
            DO j = 1, SIZE(zet_ri, 2)
               IF (zet_ri(i, j) > 1.0E-3_dp) then
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
         WRITE (bs_env%unit_nr, '(T2,A)') 'Per-kind RI-RS basis radii (Å):'
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
                                                            ikind, j, natom, nkind, start_idx, &
                                                            total_grid_npts
      INTEGER, ALLOCATABLE                               :: atom_order(:), atom_to_kind(:), &
                                                            ri_rs_grid_offsets(:)
      REAL(KIND=dp)                                      :: atomic_center(3)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

       !! Get the information about the atoms in the system
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

       !! 1. Build the grid cache
      CALL build_grid_cache(bs_env, atomic_kind_set)

       !! 2. Calculate grid counts and offsets
      ALLOCATE (ri_rs_grid_offsets(natom + 1))
      ALLOCATE (atom_to_kind(natom))

      ! grid_atom_boundaries(k) = 1-based start of the k-th atom's grid run in LAYOUT order
      ! (the order in which the points are laid down below). Used later to build atom-aligned
      ! DBCSR grid row-blocks in atomic_basis_at_grid_point.
      IF (ALLOCATED(bs_env%ri_rs%grid_atom_boundaries)) DEALLOCATE (bs_env%ri_rs%grid_atom_boundaries)
      ALLOCATE (bs_env%ri_rs%grid_atom_boundaries(natom + 1))

      ! atom -> kind map (needed before the spatial layout loop, which visits atoms in a
      ! kind-independent order).
      DO ikind = 1, nkind
         DO j = 1, SIZE(atomic_kind_set(ikind)%atom_list)
            atom_to_kind(atomic_kind_set(ikind)%atom_list(j)) = ikind
         END DO
      END DO

      ! Lay the grid rows down in Morton (spatial) atom order rather than kind-major, so a
      ! contiguous grid panel maps to a compact spatial region and the CUTOFF_RADIUS_RL_W
      ! neighborhood width of every panel shrinks. The grid row index is summed over in χ/Σ,
      ! so this reordering is result-preserving (see spatial_atom_order).
      CALL spatial_atom_order(particle_set, atom_order)

      total_grid_npts = 0
      DO i_layout = 1, natom
         atom_idx = atom_order(i_layout)
         ikind = atom_to_kind(atom_idx)
         ri_rs_grid_offsets(atom_idx) = total_grid_npts + 1
         bs_env%ri_rs%grid_atom_boundaries(i_layout) = total_grid_npts + 1
         total_grid_npts = total_grid_npts + bs_env%ri_rs%grid_cache(ikind)%npts
      END DO

      ri_rs_grid_offsets(natom + 1) = total_grid_npts + 1
      bs_env%ri_rs%grid_atom_boundaries(natom + 1) = total_grid_npts + 1

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, FMT="(T2,A,T69,I12)") &
            'Total grid points used for RI-RS:', total_grid_npts
         WRITE (bs_env%unit_nr, "(A)") ' '
      END IF

       !! 3. Allocate the global ri_rs_grid arrays
      ALLOCATE (ri_rs_grid_points(3, total_grid_npts))

       !! 4. Parallelize the grid generation loop
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(ri_rs_grid_points, ri_rs_grid_offsets, atom_to_kind, &
      !$OMP        particle_set, bs_env, natom) &
      !$OMP PRIVATE(atom_idx, ikind, atomic_center, start_idx, end_idx) &
      !$OMP SCHEDULE(DYNAMIC, 1)
      DO atom_idx = 1, natom
         ikind = atom_to_kind(atom_idx)
         atomic_center(:) = particle_set(atom_idx)%r(:)

         start_idx = ri_rs_grid_offsets(atom_idx)
         end_idx = start_idx + bs_env%ri_rs%grid_cache(ikind)%npts - 1

         !! Shift the cached origin grid by the atom's center
         ri_rs_grid_points(1, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(1, :) + atomic_center(1)
         ri_rs_grid_points(2, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(2, :) + atomic_center(2)
         ri_rs_grid_points(3, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(3, :) + atomic_center(3)

      END DO
      !$OMP END PARALLEL DO

       !! 5. Cleanup memory
      IF (ALLOCATED(bs_env%ri_rs%grid_cache)) THEN
         DO ikind = 1, nkind
            IF (ALLOCATED(bs_env%ri_rs%grid_cache(ikind)%raw_points)) DEALLOCATE (bs_env%ri_rs%grid_cache(ikind)%raw_points)
         END DO
         DEALLOCATE (bs_env%ri_rs%grid_cache)
      END IF

      DEALLOCATE (atom_order, atom_to_kind, ri_rs_grid_offsets)

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

      CHARACTER(LEN=default_path_length)                 :: filename, full_path, line
      CHARACTER(LEN=default_string_length)               :: atom_sym, suffix
      INTEGER                                            :: colon_idx, handle, i, ierr, ikind, &
                                                            iunit, nkind
      REAL(KIND=dp)                                      :: pt(3)

      CALL timeset(routineN, handle)

      ! Determine the file suffix based on the user's choice
      IF (bs_env%ri_rs%grid_select == 1) THEN
         suffix = "_def2-tzvp-rs.ion"
      ELSE IF (bs_env%ri_rs%grid_select == 2) THEN
         suffix = "_cc-pvtz-rs.ion"
      ELSE
         CPABORT("Unknown grid_select value. Valid options are 1 (def2-TZVPP) or 2 (cc-pVTZ).")
      END IF

      nkind = SIZE(atomic_kind_set)
      IF (.NOT. ALLOCATED(bs_env%ri_rs%grid_cache)) ALLOCATE (bs_env%ri_rs%grid_cache(nkind))

      DO ikind = 1, nkind
         atom_sym = TRIM(atomic_kind_set(ikind)%element_symbol)
         filename = TRIM(atom_sym)//TRIM(suffix)

         full_path = "ri_rs_grid/"//TRIM(filename)

         CALL open_file(file_name=TRIM(full_path), unit_number=iunit, &
                        file_action="READ", file_status="OLD")

         ! Parse the preamble for 'n points'
         bs_env%ri_rs%grid_cache(ikind)%npts = 0
         DO
            READ (iunit, '(A)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            IF (INDEX(line, 'n points') > 0) THEN
               colon_idx = INDEX(line, ':')
               READ (line(colon_idx + 1:), *) bs_env%ri_rs%grid_cache(ikind)%npts
               EXIT
            END IF
         END DO

         ! Allocate the cache array for this specific element
         ALLOCATE (bs_env%ri_rs%grid_cache(ikind)%raw_points(3, bs_env%ri_rs%grid_cache(ikind)%npts))

         REWIND (iunit)
         DO
            READ (iunit, '(A)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            IF (INDEX(line, '<grid_points>') > 0) EXIT
         END DO

         ! Read the raw grid coordinates
         DO i = 1, bs_env%ri_rs%grid_cache(ikind)%npts
            READ (iunit, *, IOSTAT=ierr) pt(1), pt(2), pt(3)
            IF (ierr /= 0) THEN
               CPABORT("Unexpected EOF in grid file ")
            END IF
            bs_env%ri_rs%grid_cache(ikind)%raw_points(:, i) = pt(:)
         END DO

         CALL close_file(unit_number=iunit)
      END DO

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
      INTEGER, DIMENSION(:), POINTER                     :: c_blk_sizes, col_dist, &
                                                            r_blk_sizes, row_dist
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

      natom = SIZE(particle_set)
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

      ! Overflow-safe upper bound on the block size (see dbcsr_msg_elem_limit).
      safe_max = INT(0.5_dp*REAL(dbcsr_msg_elem_limit, dp)*REAL(MAX(MIN(nprow, npcol), 1), dp)/ &
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
         WRITE (bs_env%unit_nr, '(T2,A,T69,I12)') 'RI-RS grid row-blocks of phi(mu,l)', num_grid_chunks
         WRITE (bs_env%unit_nr, '(T2,A,T69,I12)') 'RI-RS grid points per block (max)', bs_eff
      END IF

      ! Cache atomic positions: AO and RI blocks are one-block-per-atom, so these are the
      ! block centers used by the optional CUTOFF_RADIUS_G_W atom-pair truncation.
      IF (ALLOCATED(bs_env%ri_rs%atom_centers)) DEALLOCATE (bs_env%ri_rs%atom_centers)
      ALLOCATE (bs_env%ri_rs%atom_centers(3, natom))
      DO iatom = 1, natom
         bs_env%ri_rs%atom_centers(1:3, iatom) = particle_set(iatom)%r(1:3)
      END DO

      ! Cache per-chunk centroids for the optional CUTOFF_RADIUS_RL_W block truncation.
      IF (bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp) THEN
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

      CALL print_matrix_occupation(mat_phi_mu_l, 'φ(μ,l)', para_env, bs_env%unit_nr)

      ! -------------------------------------------------------------------------
      ! CLEANUP
      ! -------------------------------------------------------------------------
      DEALLOCATE (first_sgf, r_blk_sizes, c_blk_sizes, row_dist, col_dist, blk_row_start)
      CALL dbcsr_distribution_release(dist)

      CALL timestop(handle)

   END SUBROUTINE atomic_basis_at_grid_point

! **************************************************************************************************
!> \brief Evaluates all spherical AO basis functions of one atom on a set of grid points and
!>        ACCUMULATES them into phi_val (+=). For each point within the cutoff radius,
!>          Φ_μ(r) = Σ_pgf Σ_cart sphi(cart,μ) · (x-X_A)^lx (y-Y_A)^ly (z-Z_A)^lz · e^(-ζ_pgf |r-R_A|²)
!>        i.e. contracted Cartesian Gaussians transformed to the spherical basis via the sphi
!>        coefficients. Distances are minimum-image wrapped (pbc); points with
!>        |r - R_A|² > r2_threshold are skipped since every primitive is below eps there.
!>        OMP-parallel over grid points.
!> \param phi_val ...
!> \param ri_rs_grid ...
!> \param npts ...
!> \param iatom ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param r2_threshold ...
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

      INTEGER                                            :: first_sgf, handle, i_pt, ico, iend_co, &
                                                            ikind, ipgf, iset, isgf, ishell, &
                                                            istart_co, l, last_sgf, lx, ly, lz, &
                                                            n_cart_total, row_idx
      REAL(KIND=dp)                                      :: alpha, dist_vec(3), exp_val, poly, r2, &
                                                            r_atom(3), weight
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      CALL timeset(routineN, handle)

      ! Get Atom Info
      ikind = particle_set(iatom)%atomic_kind%kind_number
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
      IF (.NOT. ASSOCIATED(orb_basis_set)) THEN
         CALL timestop(handle)
         RETURN
      END IF

      r_atom = particle_set(iatom)%r

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(phi_val, ri_rs_grid, npts, orb_basis_set, r_atom, cell, ncoset, indco, &
      !$OMP        r2_threshold) &
      !$OMP PRIVATE(i_pt, dist_vec, r2, iset, n_cart_total, ishell, l, istart_co, iend_co, &
      !$OMP         first_sgf, last_sgf, ipgf, alpha, exp_val, isgf, ico, row_idx, weight, &
      !$OMP         lx, ly, lz, poly) &
      !$OMP SCHEDULE(STATIC)
      DO i_pt = 1, npts
         dist_vec = pbc(ri_rs_grid(:, i_pt) - r_atom, cell)
         r2 = DOT_PRODUCT(dist_vec, dist_vec)
         IF (r2 > r2_threshold) CYCLE

         DO iset = 1, orb_basis_set%nset
            n_cart_total = ncoset(orb_basis_set%lmax(iset))

            DO ishell = 1, orb_basis_set%nshell(iset)
               l = orb_basis_set%l(ishell, iset)
               istart_co = ncoset(l - 1) + 1
               iend_co = ncoset(l)

               first_sgf = orb_basis_set%first_sgf(ishell, iset)
               last_sgf = orb_basis_set%last_sgf(ishell, iset)

               DO ipgf = 1, orb_basis_set%npgf(iset)
                  alpha = orb_basis_set%zet(ipgf, iset)
                  exp_val = EXP(-alpha*r2)

                  DO isgf = first_sgf, last_sgf
                     DO ico = istart_co, iend_co
                        row_idx = (ipgf - 1)*n_cart_total + ico
                        weight = orb_basis_set%sphi(row_idx, isgf)
                        lx = indco(1, ico)
                        ly = indco(2, ico)
                        lz = indco(3, ico)
                        poly = (dist_vec(1)**lx)*(dist_vec(2)**ly)*(dist_vec(3)**lz)

                        phi_val(i_pt, isgf) = phi_val(i_pt, isgf) + (weight*poly*exp_val)

                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE fill_phi_for_atom

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
!> \brief Computes the RI-RS fitting coefficients Z_lP by solving, independently for every RI
!>        atom P, a Jacobi-conditioned, Tikhonov-regularized linear system restricted to the
!>        grid points r_l inside P's integration sphere |r_l - R_P| <= cutoff_ri(P):
!>          D_ll'  = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]²             (squared grid Gram matrix, Eq. 13)
!>          D_lP   = Σ_μν Φ_μ(r_l) Φ_ν(r_l) (μν|P)           (grid-RI right-hand side, Eq. 15)
!>          d_l    = 1 / sqrt(D_ll)                          (Jacobi conditioning vector)
!>          D'_ll' = d_l D_ll' d_l' + λ δ_ll'                (λ = TIKHONOV_SIGMA regularization)
!>          Σ_l' D'_ll' Z'_l'P = d_l D_lP                    (Cholesky solve, Eq. 14)
!>          Z_lP   = d_l Z'_l'P                              (undo the conditioning)
!>        Work is distributed over atoms in two phases (planned by classify_z_lp_atoms and
!>        lpt_assign_atoms): Phase A solves "small" atoms with single-rank LAPACK
!>        (dpotrf/dpotrs); Phase B solves "big" atoms, whose dense Gram matrix would exceed one
!>        rank's memory, with ScaLAPACK (pdpotrf/pdpotrs) over rank subgroups of size G.
!>        The solved Z columns are scattered into the sparse global mat_Z_lP.
!>        If a Z_lP restart file exists, it is read instead and the solve is skipped entirely.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************
   SUBROUTINE compute_coeff_Z_lP(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
                                     routineN = 'compute_coeff_Z_lP'

      INTEGER :: atom_j_mepos, atom_j_stride, atom_P, G, handle, handle_dpotrf, handle_dpotrs, &
                 i_blk, idx, info, iphase, j, max_ao_size, my_group, n_ao_total, n_big, n_done, &
                 n_groups, n_loc_ri, n_local_grid, n_my_atoms, n_small, natom, next_pct, &
                 npcol_phi, num_grid_chunks, P_loop_atom, phase_hi
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: big_list, local_grid_idx, &
                                                            my_atoms_A, my_atoms_B, &
                                                            n_local_grid_atom, row_offset, small_list
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_ri, r_blk_sizes, &
                                                            ri_blk_sizes, row_dist_grid
      LOGICAL                                            :: do_scatter, use_dist
      REAL(KIND=dp)                                      :: balance_A, balance_B, cutoff_ri, r_c, t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoff_ri_per_atom, d_vec_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: D_local, d_lp_local, phi_local
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_b, fm_struct_D
      TYPE(cp_fm_type)                                   :: fm_b, fm_D
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist_phi, dist_Z
      TYPE(gw_3c_ctx_type)                               :: ctx_3c
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, input=input, &
                      qs_kind_set=qs_kind_set, cell=cell, atomic_kind_set=atomic_kind_set)

      NULLIFY (para_env_sub, blacs_env_sub)

      natom = SIZE(bs_env%i_RI_start_from_atom)
      n_ao_total = bs_env%i_ao_end_from_atom(natom)

      ! =========================================================================
      ! 1. SETUP DBCSR TOPOLOGY & EXACT OFFSETS
      !    mat_Z_lP inherits the grid row blocking (and row distribution) of
      !    mat_phi_mu_l; its columns are one block per RI atom.
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, row_blk_size=r_blk_sizes, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, npcols=npcol_phi)

      num_grid_chunks = SIZE(r_blk_sizes)

      ALLOCATE (row_offset(num_grid_chunks))
      row_offset(1) = 0
      DO i_blk = 2, num_grid_chunks
         row_offset(i_blk) = row_offset(i_blk - 1) + r_blk_sizes(i_blk - 1)
      END DO

      ALLOCATE (ri_blk_sizes(natom), col_dist_ri(natom))
      DO atom_P = 1, natom
         ri_blk_sizes(atom_P) = bs_env%i_RI_end_from_atom(atom_P) - bs_env%i_RI_start_from_atom(atom_P) + 1
         col_dist_ri(atom_P) = MOD(atom_P - 1, npcol_phi)
      END DO

      CALL dbcsr_distribution_new(dist_Z, template=dist_phi, row_dist=row_dist_grid, col_dist=col_dist_ri)

      IF (bs_env%ri_rs%Z_lP_exists) THEN
         CALL dbcsr_binary_read(filepath=TRIM(bs_env%prefix)//"Z_lP.matrix", &
                                distribution=dist_Z, &
                                matrix_new=mat_Z_lP)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Read Z_lP from file ', ' Execution time', m_walltime() - t1, ' s'
            ! The grid rows are laid out in Morton (spatial) order (spatial_atom_order); a Z_lP.matrix
            ! written by an older build with a different grid ordering would be silently misread into
            ! the current row order. Delete stale Z_lP.matrix files and recompute if in doubt.
            WRITE (bs_env%unit_nr, '(T2,A)') &
               '*** NOTE: Z_lP restart must match the current (spatial) grid row ordering ***'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
      ELSE

         CALL dbcsr_create(mat_Z_lP, name="mat_Z_lP", dist=dist_Z, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=r_blk_sizes, col_blk_size=ri_blk_sizes)

         ! Largest per-atom AO block, needed to size the 3c-integral work buffers.
         max_ao_size = 0
         DO j = 1, SIZE(bs_env%i_ao_start_from_atom)
            max_ao_size = MAX(max_ao_size, bs_env%i_ao_end_from_atom(j) - bs_env%i_ao_start_from_atom(j) + 1)
         END DO

         ! Per-atom RI-RS integration sphere:
         !   cutoff_ri(P) = r_c + r_RI(P)
         ! where r_c is the truncated-Coulomb cutoff of the RI metric and r_RI the radius of
         ! the most diffuse RI auxiliary Gaussian on P. The CUTOFF_RADIUS_RL_RI keyword
         ! (when > 0) overrides the entire cutoff calculation.
         ALLOCATE (cutoff_ri_per_atom(natom))

         IF (bs_env%ri_rs%cutoff_radius_ri_rs > 0.0_dp) THEN
            cutoff_ri_per_atom(:) = bs_env%ri_rs%cutoff_radius_ri_rs
         ELSE
            r_c = bs_env%ri_metric%cutoff_radius
            DO P_loop_atom = 1, natom
               cutoff_ri_per_atom(P_loop_atom) = &
                  r_c + bs_env%ri_rs%radius_ri_per_atom(P_loop_atom)
            END DO
         END IF

         CALL print_sphere_cutoff_table(bs_env, atomic_kind_set, particle_set, cutoff_ri_per_atom)

         ! =========================================================================
         ! 2. PER-ATOM SOLVER CLASSIFICATION
         !    Split the atoms into "small" (single-rank LAPACK, Phase A) and "big"
         !    (distributed ScaLAPACK over subgroups of G ranks, Phase B) by comparing
         !    each atom's estimated solve peak memory against the measured budget.
         ! =========================================================================
         CALL classify_z_lp_atoms(bs_env, para_env, ri_rs_grid_points, particle_set, &
                                  cutoff_ri_per_atom, ri_blk_sizes, n_ao_total, &
                                  n_local_grid_atom, small_list, n_small, big_list, n_big, G)

         ! LPT scheduling: sort the atoms of each phase by estimated solve cost
         ! (n_local_grid^3, Cholesky-dominated) and greedily assign to the least-loaded
         ! rank (Phase A) / subgroup (Phase B).
         CALL lpt_assign_atoms(small_list, n_small, n_local_grid_atom, para_env%num_pe, &
                               para_env%mepos, my_atoms_A, balance_A)
         IF (n_big > 0) THEN
            n_groups = para_env%num_pe/G
            my_group = MIN(para_env%mepos/G, n_groups - 1)
            CALL lpt_assign_atoms(big_list, n_big, n_local_grid_atom, n_groups, my_group, &
                                  my_atoms_B, balance_B)
         ELSE
            ALLOCATE (my_atoms_B(0))
            balance_B = 1.0_dp
         END IF

         ! Atoms this rank will process across both phases for rank-0 progress
         n_my_atoms = SIZE(my_atoms_A) + SIZE(my_atoms_B)
         n_done = 0
         next_pct = 25

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I7,A,I8,A)') &
               'RI-RS Z_lP solver: ', n_small, ' atoms single-rank (BLAS), ', n_big, &
               ' atoms distributed'
            IF (n_small > 0) WRITE (bs_env%unit_nr, '(T4,A,F18.2)') &
               'estimated single-rank load balance (max/mean cost per rank)', balance_A
            IF (n_big > 0) THEN
               WRITE (bs_env%unit_nr, '(T4,A,I44,A)') &
                  'distributed subgroup size G', G, ' ranks'
               WRITE (bs_env%unit_nr, '(T4,A,F17.2)') &
                  'estimated distributed load balance (max/mean cost per group)', balance_B
            END IF
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

         ! Shared context for the three-center integrals (μν|P) of the RHS build
         CALL gw_3c_ctx_create(ctx_3c, qs_env, bs_env%ri_metric, &
                               basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                               basis_i=bs_env%basis_set_RI)

         ! =========================================================================
         ! 3. TWO-PHASE LOOP OVER ATOMS
         !    Phase A processes the "small" atoms with the single-rank BLAS path
         !    Phase B processes the "big" atoms with the distributed ScaLAPACK path
         !    over rank subgroups of size G. phi_local for each atom's cutoff sphere
         !    is built on the fly to avoid replicating a global grid x AO matrix.
         ! =========================================================================
         DO iphase = 1, 2
            IF (iphase == 1) THEN
               use_dist = .FALSE.
               atom_j_mepos = 0
               atom_j_stride = 1
               phase_hi = SIZE(my_atoms_A)
            ELSE
               IF (n_big == 0) CYCLE
               use_dist = .TRUE.
               n_groups = para_env%num_pe/G
               my_group = MIN(para_env%mepos/G, n_groups - 1)
               ALLOCATE (para_env_sub)
               CALL para_env_sub%from_split(para_env, my_group)
               CALL cp_blacs_env_create(blacs_env=blacs_env_sub, para_env=para_env_sub)
               atom_j_mepos = para_env_sub%mepos
               atom_j_stride = para_env_sub%num_pe
               ! All ranks of a subgroup share my_group, hence the identical my_atoms_B list
               ! (the per-atom ScaLAPACK solve is collective over the subgroup).
               phase_hi = SIZE(my_atoms_B)
            END IF

            DO idx = 1, phase_hi
               IF (iphase == 1) THEN
                  atom_P = my_atoms_A(idx)
               ELSE
                  atom_P = my_atoms_B(idx)
               END IF

               n_loc_ri = ri_blk_sizes(atom_P)
               cutoff_ri = cutoff_ri_per_atom(atom_P)

               ! ---------------------------------------------------------------------
               ! A. Sphere-local AO matrix Φ_μ(r_l): select the grid points with
               !    |r_l - R_P| <= cutoff_ri(P), evaluate every AO on them, and drop
               !    points whose largest AO amplitude is below EPS_FILTER.
               ! ---------------------------------------------------------------------
               CALL build_phi_on_sphere(bs_env, particle_set, qs_kind_set, cell, &
                                        ri_rs_grid_points, atom_P, cutoff_ri, n_ao_total, &
                                        local_grid_idx, n_local_grid, phi_local)

               ! ---------------------------------------------------------------------
               ! B. Right-hand side  D_lP = Σ_μν Φ_μ(r_l) Φ_ν(r_l) (μν|P)
               ! ---------------------------------------------------------------------
               ALLOCATE (d_lp_local(n_local_grid, n_loc_ri))
               d_lp_local = 0.0_dp

               CALL compute_d_lp(bs_env, ctx_3c, phi_local, d_lp_local, n_local_grid, &
                                 n_loc_ri, atom_P, max_ao_size, atom_j_mepos, atom_j_stride)

               ! Reduce per-subgroup-rank partials into the replicated d_lp_local.
               ! Skipped for BLAS path: each rank has the full sum locally.
               IF (use_dist) THEN
                  CALL para_env_sub%sum(d_lp_local)
               END IF

               ! ---------------------------------------------------------------------
               ! C. Jacobi conditioning vector d_l = 1/sqrt(D_ll) and, on the BLAS path,
               !    the dense conditioned Gram matrix D'_ll' = d_l D_ll' d_l' + λδ_ll'.
               ! ---------------------------------------------------------------------
               ALLOCATE (d_vec_local(n_local_grid))

               IF (.NOT. use_dist) THEN
                  CALL build_gram_jacobi_blas(phi_local, n_local_grid, n_ao_total, &
                                              bs_env%ri_rs%tikhonov, D_local, d_vec_local)
               ELSE
                  ! ScaLAPACK path: only d_vec is needed here (= 1/||phi(r_l)||^2);
                  ! solve_D_lp_distributed builds its block-cyclic slice of D' internally
                  ! with the squared+scaled values, so no dense D_local on this rank.
                  CALL build_jacobi_diag_from_phi(phi_local, n_local_grid, n_ao_total, &
                                                  d_vec_local)
               END IF

               ! ---------------------------------------------------------------------
               ! D. Pre-scale the RHS:  D'_lP = d_l * D_lP
               ! ---------------------------------------------------------------------
               CALL scale_rows_by_diag(d_lp_local, d_vec_local, n_local_grid, n_loc_ri)

               ! ---------------------------------------------------------------------
               ! E. Cholesky solve  Σ_l' D'_ll' Z'_l'P = D'_lP
               !    (BLAS dpotrf/dpotrs or ScaLAPACK pdpotrf/pdpotrs)
               ! ---------------------------------------------------------------------
               IF (.NOT. use_dist) THEN
                  CALL timeset(routineN//"_dpotrf", handle_dpotrf)
                  CALL dpotrf('L', n_local_grid, D_local, n_local_grid, info)
                  CALL timestop(handle_dpotrf)
                  CALL timeset(routineN//"_dpotrs", handle_dpotrs)
                  CALL dpotrs('L', n_local_grid, n_loc_ri, D_local, n_local_grid, &
                              d_lp_local, n_local_grid, info)
                  CALL timestop(handle_dpotrs)
                  DEALLOCATE (D_local)
               ELSE
                  CALL solve_D_lp_distributed(phi_local, d_vec_local, d_lp_local, &
                                              n_local_grid, n_ao_total, n_loc_ri, &
                                              bs_env%ri_rs%tikhonov, &
                                              para_env_sub, blacs_env_sub, &
                                              fm_struct_D, fm_struct_b, fm_D, fm_b, info)
               END IF

               ! ---------------------------------------------------------------------
               ! F. Undo the conditioning:  Z_lP = d_l * Z'_lP
               ! ---------------------------------------------------------------------
               CALL scale_rows_by_diag(d_lp_local, d_vec_local, n_local_grid, n_loc_ri)

               ! ---------------------------------------------------------------------
               ! G. Scatter the solved Z columns back into the global sparse mat_Z_lP.
               ! ---------------------------------------------------------------------
               do_scatter = .TRUE.
               IF (use_dist) do_scatter = (para_env_sub%mepos == 0)
               IF (do_scatter) THEN
                  CALL scatter_z_columns(mat_Z_lP, d_lp_local, local_grid_idx, n_local_grid, &
                                         n_loc_ri, atom_P, r_blk_sizes, row_offset, &
                                         bs_env%eps_filter)
               END IF

               DEALLOCATE (d_vec_local, d_lp_local)
               DEALLOCATE (local_grid_idx, phi_local)

               ! Progress based on rank 0 only
               n_done = n_done + 1
               IF (bs_env%unit_nr > 0 .AND. n_my_atoms > 0) THEN
                  DO WHILE (next_pct <= 100 .AND. n_done*100 >= next_pct*n_my_atoms)
                     WRITE (bs_env%unit_nr, '(T2,A,I57,A)') 'Computing Z_lP:', next_pct, ' % done'
                     next_pct = next_pct + 25
                  END DO
               END IF
            END DO   ! idx: atoms of this phase owned by this rank / subgroup

            ! Tear down the Phase-B subgroup (all ranks created it collectively).
            IF (iphase == 2) THEN
               CALL cp_blacs_env_release(blacs_env_sub)
               CALL para_env_sub%free()
               DEALLOCATE (para_env_sub)
            END IF
         END DO   ! iphase

         DEALLOCATE (cutoff_ri_per_atom)
         DEALLOCATE (small_list, big_list)

         CALL gw_3c_ctx_release(ctx_3c)

         CALL dbcsr_finalize(mat_Z_lP)

         CALL print_matrix_occupation(mat_Z_lP, 'Z(l,P)', para_env, bs_env%unit_nr)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(A)') ' '
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed Z_lP ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

         logger => cp_get_default_logger()

         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN
            CALL dbcsr_binary_write(matrix=mat_Z_lP, filepath=TRIM(bs_env%prefix)//"Z_lP.matrix")
         END IF

      END IF

      DEALLOCATE (row_offset, ri_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)

      DEALLOCATE (ri_rs_grid_points)

      CALL timestop(handle)

   END SUBROUTINE compute_coeff_Z_lP

! **************************************************************************************************
!> \brief Prints the per-kind maximum RI-RS integration-sphere cutoff table.
!> \param bs_env ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param cutoff_ri_per_atom ...
! **************************************************************************************************
   SUBROUTINE print_sphere_cutoff_table(bs_env, atomic_kind_set, particle_set, cutoff_ri_per_atom)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cutoff_ri_per_atom
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'print_sphere_cutoff_table'

      INTEGER                                            :: handle

      INTEGER                                            :: iatom, ikind, nkind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoff_ri_per_kind

      CALL timeset(routineN, handle)

      IF (bs_env%unit_nr <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      nkind = SIZE(atomic_kind_set)
      ALLOCATE (cutoff_ri_per_kind(nkind))
      cutoff_ri_per_kind(:) = 0.0_dp

      DO iatom = 1, SIZE(particle_set)
         ikind = particle_set(iatom)%atomic_kind%kind_number
         cutoff_ri_per_kind(ikind) = MAX(cutoff_ri_per_kind(ikind), cutoff_ri_per_atom(iatom))
      END DO

      WRITE (bs_env%unit_nr, '(T2,A)') 'Per-kind maximum RI-RS sphere cutoff (Å):'
      WRITE (bs_env%unit_nr, '(T4,A4,A14)') 'Kind', 'cutoff (Å)'
      DO ikind = 1, nkind
         WRITE (bs_env%unit_nr, '(T4,A4,F14.4)') &
            atomic_kind_set(ikind)%element_symbol, &
            cutoff_ri_per_kind(ikind)*angstrom
      END DO
      WRITE (bs_env%unit_nr, '(A)') ' '

      DEALLOCATE (cutoff_ri_per_kind)

      CALL timestop(handle)

   END SUBROUTINE print_sphere_cutoff_table

! **************************************************************************************************
!> \brief Splits the atoms of the Z_lP solve into a single-rank list ("small", Phase A: LAPACK
!>        dpotrf/dpotrs on one rank) and a distributed list ("big", Phase B: ScaLAPACK
!>        pdpotrf/pdpotrs over a rank subgroup of size G), and sizes G.
!>        AUTO mode (N_PROCS_PER_ATOM_Z_LP <= 0, the default): estimate each atom's single-rank
!>        peak memory
!>          peak(P) = 8*n_local_grid(P)^2                (dense Gram matrix D_local)
!>                  + 8*n_local_grid(P)*n_ao             (phi_local)
!>                  + 8*n_local_grid(P)*n_RI(P)*(1+n_threads)   (d_lp + OMP partials)
!>        and send atoms whose peak exceeds mem_safety * available-memory-per-proc to the
!>        distributed path; G is auto-sized so the biggest atom's distributed D_local (/G)
!>        fits alongside the replicated phi_local + d_lp.
!>        MANUAL mode (> 0): 1 forces the single-rank path for every atom; > 1 keeps the
!>        memory-based classification but forces that fixed subgroup size G.
!>        In every mode G is floored by the ScaLAPACK 32-bit index limit (a local block-cyclic
!>        slice of ~n_local_grid^2/G elements must stay below 2^31 or pdpotrf segfaults).
!> \param bs_env ...
!> \param para_env ...
!> \param ri_rs_grid_points ...
!> \param particle_set ...
!> \param cutoff_ri_per_atom ...
!> \param ri_blk_sizes per-atom ...
!> \param n_ao_total ...
!> \param n_local_grid_atom ...
!> \param small_list ...
!> \param n_small ...
!> \param big_list ...
!> \param n_big ...
!> \param G ...
! **************************************************************************************************
   SUBROUTINE classify_z_lp_atoms(bs_env, para_env, ri_rs_grid_points, particle_set, &
                                  cutoff_ri_per_atom, ri_blk_sizes, n_ao_total, &
                                  n_local_grid_atom, small_list, n_small, big_list, n_big, G)

!$    USE OMP_LIB, ONLY: omp_get_max_threads

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cutoff_ri_per_atom
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ri_blk_sizes
      INTEGER, INTENT(IN)                                :: n_ao_total
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: n_local_grid_atom, small_list, big_list
      INTEGER, INTENT(OUT)                               :: n_small, n_big, G
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'classify_z_lp_atoms'

      INTEGER                                            :: handle

      ! Conservative fraction of measured available memory usable per rank for the Z_lP
      REAL(KIND=dp), PARAMETER                           :: mem_safety = 0.8_dp

      ! ScaLAPACK/BLACS index the per-rank local block-cyclic slice (~n_local_grid^2/G
      ! elements) with 32-bit integers; keep it safely below 2^31 or pdpotrf segfaults.
      REAL(KIND=dp), PARAMETER                           :: scalapack_loc_limit = 2.0E9_dp

      INTEGER                                            :: G_atom, G_int32, G_int32_max, l, &
                                                            n_grid_total, n_local_grid, natom, &
                                                            nthreads_cls, P_loop_atom
      LOGICAL                                            :: auto_mode
      REAL(KIND=dp)                                      :: budget_bytes, cutoff_ri, dlp_bytes, &
                                                            mem_avail_GB, ng, nri, peak_bytes, &
                                                            phi_bytes
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)
      n_grid_total = SIZE(ri_rs_grid_points, 2)

      ! Per-atom sphere size: n_local_grid(P) = #{ l : |r_l - R_P| <= cutoff_ri(P) }.
      ! It sets both the memory footprint (D_local is n_local_grid^2) and the solve
      ! cost (~n_local_grid^3), so it drives classification and the LPT load balancing.
      ALLOCATE (n_local_grid_atom(natom), small_list(natom), big_list(natom))
      DO P_loop_atom = 1, natom
         pos_P(:) = particle_set(P_loop_atom)%r(:)
         cutoff_ri = cutoff_ri_per_atom(P_loop_atom)
         n_local_grid = 0
         DO l = 1, n_grid_total
            IF (SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2) <= cutoff_ri**2) then
               n_local_grid = n_local_grid + 1
            end if
         END DO
         n_local_grid_atom(P_loop_atom) = n_local_grid
      END DO

      nthreads_cls = 1
!$    nthreads_cls = omp_get_max_threads()
      ! N_PROCS_PER_ATOM_Z_LP: -1 (default) = AUTO (classify by memory, auto-size G);
      ! 1 = force single-rank BLAS for every atom; >1 = classify by memory but use this
      ! fixed subgroup size G for the big atoms.
      auto_mode = (bs_env%ri_rs%n_procs_per_atom_z_lp <= 0)
      CALL ri_rs_mem_avail_per_proc_GB(bs_env, mem_avail_GB)  ! collective over all ranks
      budget_bytes = mem_safety*mem_avail_GB*1.0E9_dp

      n_small = 0
      n_big = 0
      G = 1
      G_atom = 1              ! max G a big atom needs (memory + ScaLAPACK int32 floor)
      G_int32_max = 1         ! max ScaLAPACK-int32 floor over the distributed atoms
      IF (bs_env%ri_rs%n_procs_per_atom_z_lp == 1) THEN
         ! Force single-rank BLAS for every atom.
         DO P_loop_atom = 1, natom
            n_small = n_small + 1
            small_list(n_small) = P_loop_atom
         END DO
      ELSE IF (mem_avail_GB <= 0.0_dp) THEN
         ! No /proc/meminfo => cannot size by memory.
         IF (auto_mode) THEN
            IF (bs_env%unit_nr > 0) then
               CPWARN("RI-RS Z_lP: no meminfo; single-rank solve for all atoms")
            end if
            DO P_loop_atom = 1, natom
               n_small = n_small + 1
               small_list(n_small) = P_loop_atom
            END DO
         ELSE
            ! Fixed G, no meminfo: distribute all atoms; still floor G by the int32 limit.
            DO P_loop_atom = 1, natom
               ng = REAL(n_local_grid_atom(P_loop_atom), dp)
               G_int32_max = MAX(G_int32_max, CEILING(ng*ng/scalapack_loc_limit))
               n_big = n_big + 1
               big_list(n_big) = P_loop_atom
            END DO
            G = MIN(bs_env%ri_rs%n_procs_per_atom_z_lp, para_env%num_pe)
            IF (G < G_int32_max) THEN
               G = MIN(G_int32_max, para_env%num_pe)
               IF (bs_env%unit_nr > 0) then
                  CPWARN("RI-RS Z_lP: raised G to avoid ScaLAPACK overflow")
               end if
            END IF
         END IF
      ELSE
         ! Classify by memory: peak (D_local + phi_local + d_lp) vs budget. Small -> BLAS,
         ! big -> distributed. Same classification for AUTO and fixed-G modes.
         DO P_loop_atom = 1, natom
            ng = REAL(n_local_grid_atom(P_loop_atom), dp)
            nri = REAL(ri_blk_sizes(P_loop_atom), dp)
            phi_bytes = 8.0_dp*ng*REAL(n_ao_total, dp)
            dlp_bytes = 8.0_dp*ng*nri*REAL(1 + nthreads_cls, dp)
            peak_bytes = 8.0_dp*ng*ng + phi_bytes + dlp_bytes
            IF (peak_bytes <= budget_bytes) THEN
               n_small = n_small + 1
               small_list(n_small) = P_loop_atom
            ELSE
               n_big = n_big + 1
               big_list(n_big) = P_loop_atom
               ! G must satisfy BOTH: (a) memory — distributed D_local (/G) fits next to the
               ! replicated phi_local + d_lp; (b) ScaLAPACK — local ~ng^2/G below the int32 limit.
               G_int32 = CEILING(ng*ng/scalapack_loc_limit)
               G_int32_max = MAX(G_int32_max, G_int32)
               G_atom = MAX(G_atom, G_int32, &
                            CEILING(8.0_dp*ng*ng/MAX(budget_bytes - phi_bytes - dlp_bytes, 1.0_dp)))
            END IF
         END DO
         IF (n_big > 0) THEN
            IF (auto_mode) THEN
               ! Auto-size G from the most demanding big atom.
               IF (G_atom > para_env%num_pe) THEN
                  IF (bs_env%unit_nr > 0) THEN
                     WRITE (bs_env%unit_nr, '(T2,A)') &
                        'RI-RS Z_lP: an atom is too large to fit even distributed over all '// &
                        'ranks. Add nodes, use fewer MPI ranks/node, lower CUTOFF_RADIUS_RL_RI, '// &
                        'or raise EPS_FILTER (more grid screening).'
                  END IF
                  CPABORT("RI-RS Z_lP: atom too large even fully distributed")
               END IF
               G = MIN(MAX(G_atom, 2), para_env%num_pe)
            ELSE
               ! Fixed G from the keyword. Hard-floor by the ScaLAPACK int32 limit (below it
               ! pdpotrf segfaults); warn if it is still below the memory recommendation.
               G = MIN(bs_env%ri_rs%n_procs_per_atom_z_lp, para_env%num_pe)
               IF (G < G_int32_max) THEN
                  G = MIN(G_int32_max, para_env%num_pe)
                  IF (bs_env%unit_nr > 0) then
                     CPWARN("RI-RS Z_lP: raised G to avoid ScaLAPACK overflow")
                  end if
               ELSE IF (G < G_atom .AND. bs_env%unit_nr > 0) THEN
                  CPWARN("RI-RS Z_lP: N_PROCS_PER_ATOM_Z_LP too small for the largest atom")
               END IF
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE classify_z_lp_atoms

! **************************************************************************************************
!> \brief Builds the sphere-local AO matrix phi_local(l, μ) = Φ_μ(r_l) for one RI atom P
!> \param bs_env ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param ri_rs_grid_points ...
!> \param atom_P ...
!> \param cutoff_ri ...
!> \param n_ao_total ...
!> \param local_grid_idx ...
!> \param n_local_grid ...
!> \param phi_local ...
! **************************************************************************************************
   SUBROUTINE build_phi_on_sphere(bs_env, particle_set, qs_kind_set, cell, ri_rs_grid_points, &
                                  atom_P, cutoff_ri, n_ao_total, local_grid_idx, n_local_grid, &
                                  phi_local)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points
      INTEGER, INTENT(IN)                                :: atom_P
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_ri
      INTEGER, INTENT(IN)                                :: n_ao_total
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: local_grid_idx
      INTEGER, INTENT(OUT)                               :: n_local_grid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: phi_local

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_phi_on_sphere'

      INTEGER                                            :: col_end, col_start, handle, j, k, l, &
                                                            loc_idx, n_grid_total, n_keep, &
                                                            source_atom
      REAL(KIND=dp)                                      :: d_sP, dist, r2_threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w_pt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: phi_keep, sphere_grid
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P

      CALL timeset(routineN, handle)

      n_grid_total = SIZE(ri_rs_grid_points, 2)
      pos_P(:) = particle_set(atom_P)%r(:)

      n_local_grid = 0
      DO l = 1, n_grid_total
         dist = SQRT(SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2))
         IF (dist <= cutoff_ri) n_local_grid = n_local_grid + 1
      END DO

      ALLOCATE (local_grid_idx(n_local_grid))

      n_local_grid = 0
      DO l = 1, n_grid_total
         dist = SQRT(SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2))
         IF (dist <= cutoff_ri) THEN
            n_local_grid = n_local_grid + 1
            local_grid_idx(n_local_grid) = l
         END IF
      END DO

      ALLOCATE (sphere_grid(3, n_local_grid))
      DO loc_idx = 1, n_local_grid
         sphere_grid(:, loc_idx) = ri_rs_grid_points(:, local_grid_idx(loc_idx))
      END DO

      ALLOCATE (phi_local(n_local_grid, n_ao_total))
      phi_local = 0.0_dp

      DO source_atom = 1, SIZE(particle_set)
         d_sP = NORM2(particle_set(source_atom)%r(:) - pos_P(:))
         IF (d_sP > bs_env%ri_rs%radius_ao_per_atom(source_atom) + cutoff_ri) CYCLE

         col_start = bs_env%i_ao_start_from_atom(source_atom)
         col_end = bs_env%i_ao_end_from_atom(source_atom)
         ! A positive CUTOFF_RADIUS_RI_AO overrides the per-atom Gaussian radius
         ! with a user-defined hard cutoff.
         IF (bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) THEN
            r2_threshold = bs_env%ri_rs%cutoff_radius_ri_ao**2
         ELSE
            r2_threshold = bs_env%ri_rs%radius_ao_per_atom(source_atom)**2
         END IF

         CALL fill_phi_for_atom(phi_local(:, col_start:col_end), sphere_grid, &
                                n_local_grid, source_atom, particle_set, qs_kind_set, &
                                cell, r2_threshold)
      END DO

      DEALLOCATE (sphere_grid)

      IF (n_local_grid > 0) THEN
         ALLOCATE (w_pt(n_local_grid))
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP SHARED(n_local_grid, n_ao_total, phi_local, w_pt) &
         !$OMP PRIVATE(l, j) SCHEDULE(STATIC)
         DO l = 1, n_local_grid
            w_pt(l) = 0.0_dp
            DO j = 1, n_ao_total
               w_pt(l) = MAX(w_pt(l), ABS(phi_local(l, j)))
            END DO
         END DO
         !$OMP END PARALLEL DO
         n_keep = COUNT(w_pt > bs_env%eps_filter)
         IF (n_keep < n_local_grid) THEN
            ALLOCATE (phi_keep(n_keep, n_ao_total))
            k = 0
            DO l = 1, n_local_grid
               IF (w_pt(l) > bs_env%eps_filter) THEN
                  k = k + 1
                  phi_keep(k, :) = phi_local(l, :)
                  local_grid_idx(k) = local_grid_idx(l)
               END IF
            END DO
            CALL MOVE_ALLOC(phi_keep, phi_local)
            n_local_grid = n_keep
         END IF
         DEALLOCATE (w_pt)
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_phi_on_sphere

! **************************************************************************************************
!> \brief Builds the dense Jacobi-conditioned Gram matrix and the conditioning vector for the
!>        single-rank (BLAS/LAPACK) Z_lP solve:
!>          D_ll'   = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]²     (dsyrk of phi_local, then squared)
!>          d_l     = 1 / sqrt(D_ll)
!>          D'_ll'  = d_l D_ll' d_l' + λ δ_ll'
!>        Only the lower triangle is referenced by the subsequent dpotrf('L')
!> \param phi_local ...
!> \param n_local_grid ...
!> \param n_ao_total ...
!> \param tikhonov ...
!> \param D_local ...
!> \param d_vec_local ...
! **************************************************************************************************
   SUBROUTINE build_gram_jacobi_blas(phi_local, n_local_grid, n_ao_total, tikhonov, D_local, &
                                     d_vec_local)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      INTEGER, INTENT(IN)                                :: n_local_grid, n_ao_total
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: D_local
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: d_vec_local

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_gram_jacobi_blas'

      INTEGER                                            :: handle, handle_dsyrk, i, j

      CALL timeset(routineN, handle)

      ALLOCATE (D_local(n_local_grid, n_local_grid))
      D_local = 0.0_dp

      ! D_ll' = Σ_μ Φ_μ(r_l) Φ_μ(r_l')   (lower triangle only)
      CALL timeset(routineN//"_dsyrk", handle_dsyrk)
      CALL dsyrk("L", "N", n_local_grid, n_ao_total, 1.0_dp, phi_local, &
                 n_local_grid, 0.0_dp, D_local, n_local_grid)
      CALL timestop(handle_dsyrk)

      ! Diagonal: square, derive d_l = 1/sqrt(D_ll), scale, add Tikhonov λ
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_local_grid, D_local, d_vec_local, tikhonov) &
      !$OMP PRIVATE(i) &
      !$OMP SCHEDULE(STATIC)
      DO i = 1, n_local_grid
         D_local(i, i) = D_local(i, i)**2
         d_vec_local(i) = 1.0_dp/SQRT(MAX(D_local(i, i), 1.0E-16_dp))
         D_local(i, i) = (D_local(i, i)*d_vec_local(i)**2) + tikhonov
      END DO
      !$OMP END PARALLEL DO

      ! Off-diagonal: D'_ll' = d_l D_ll'^2 d_l' (mirror to the upper triangle)
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_local_grid, D_local, d_vec_local) &
      !$OMP PRIVATE(j, i) &
      !$OMP SCHEDULE(DYNAMIC)
      DO j = 1, n_local_grid
         DO i = j + 1, n_local_grid
            D_local(i, j) = D_local(i, j)**2
            D_local(i, j) = D_local(i, j)*d_vec_local(i)*d_vec_local(j)
            D_local(j, i) = D_local(i, j)
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE build_gram_jacobi_blas

! **************************************************************************************************
!> \brief Computes the Jacobi conditioning vector directly from phi for the distributed
!>        (ScaLAPACK) Z_lP solve:  d_l = 1 / Σ_μ Φ_μ(r_l)² = 1/sqrt(D_ll), without forming
!>        the Gram matrix.
!> \param phi_local ...
!> \param n_local_grid ...
!> \param n_ao_total ...
!> \param d_vec_local ...
! **************************************************************************************************
   SUBROUTINE build_jacobi_diag_from_phi(phi_local, n_local_grid, n_ao_total, d_vec_local)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      INTEGER, INTENT(IN)                                :: n_local_grid, n_ao_total
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: d_vec_local
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_jacobi_diag_from_phi'

      INTEGER                                            :: handle

      INTEGER                                            :: i, j

      CALL timeset(routineN, handle)

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_local_grid, n_ao_total, phi_local, d_vec_local) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(STATIC)
      DO i = 1, n_local_grid
         d_vec_local(i) = 0.0_dp
         DO j = 1, n_ao_total
            d_vec_local(i) = d_vec_local(i) + phi_local(i, j)*phi_local(i, j)
         END DO
         d_vec_local(i) = 1.0_dp/MAX(d_vec_local(i), 1.0E-16_dp)
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE build_jacobi_diag_from_phi

! **************************************************************************************************
!> \brief Scales every row of a matrix by the corresponding diagonal entry,
!>        A(l, :) <- d_l * A(l, :). Used in Z_lP solve: to pre-scale the RHS
!>        (D'_lP = d_l D_lP) and to undo the conditioning of the solution (Z_lP = d_l Z'_lP).
!> \param d_lp_local ...
!> \param d_vec_local ...
!> \param n_local_grid ...
!> \param n_loc_ri ...
! **************************************************************************************************
   SUBROUTINE scale_rows_by_diag(d_lp_local, d_vec_local, n_local_grid, n_loc_ri)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp_local
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: d_vec_local
      INTEGER, INTENT(IN)                                :: n_local_grid, n_loc_ri
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scale_rows_by_diag'

      INTEGER                                            :: handle

      INTEGER                                            :: i, j_ri

      CALL timeset(routineN, handle)

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_loc_ri, n_local_grid, d_lp_local, d_vec_local) &
      !$OMP PRIVATE(j_ri, i) &
      !$OMP SCHEDULE(STATIC)
      DO j_ri = 1, n_loc_ri
         DO i = 1, n_local_grid
            d_lp_local(i, j_ri) = d_lp_local(i, j_ri)*d_vec_local(i)
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE scale_rows_by_diag

! **************************************************************************************************
!> \brief Scatters the solved Z columns of one atom P from the dense sphere-local solution back
!>        into the global sparse mat_Z_lP
!> \param mat_Z_lP ...
!> \param d_lp_local ...
!> \param local_grid_idx ...
!> \param n_local_grid ...
!> \param n_loc_ri ...
!> \param atom_P ...
!> \param r_blk_sizes ...
!> \param row_offset ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE scatter_z_columns(mat_Z_lP, d_lp_local, local_grid_idx, n_local_grid, n_loc_ri, &
                                atom_P, r_blk_sizes, row_offset, eps_filter)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_Z_lP
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: d_lp_local
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_grid_idx
      INTEGER, INTENT(IN)                                :: n_local_grid, n_loc_ri, atom_P
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes, row_offset
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scatter_z_columns'
      INTEGER                                            :: handle

      INTEGER                                            :: current_chunk_size, g_pt, i_blk, &
                                                            loc_ptr, r_end, r_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Z_blk

      CALL timeset(routineN, handle)

      ALLOCATE (Z_blk(MAXVAL(r_blk_sizes), n_loc_ri))
      loc_ptr = 1

      DO i_blk = 1, SIZE(r_blk_sizes)
         r_start = row_offset(i_blk) + 1
         r_end = row_offset(i_blk) + r_blk_sizes(i_blk)
         current_chunk_size = r_blk_sizes(i_blk)

         Z_blk = 0.0_dp

         ! Copy the sphere points whose global grid index falls inside this block
         DO WHILE (loc_ptr <= n_local_grid)
            g_pt = local_grid_idx(loc_ptr)
            IF (g_pt > r_end) EXIT
            Z_blk(g_pt - r_start + 1, 1:n_loc_ri) = d_lp_local(loc_ptr, 1:n_loc_ri)
            loc_ptr = loc_ptr + 1
         END DO

         IF (MAXVAL(ABS(Z_blk(1:current_chunk_size, 1:n_loc_ri))) > eps_filter) THEN
            CALL dbcsr_put_block(mat_Z_lP, row=i_blk, col=atom_P, &
                                 block=Z_blk(1:current_chunk_size, 1:n_loc_ri))
         END IF
      END DO

      DEALLOCATE (Z_blk)

      CALL timestop(handle)

   END SUBROUTINE scatter_z_columns

! **************************************************************************************************
!> \brief LPT (longest-processing-time) assignment of the Z_lP atoms to workers (MPI ranks in
!>        Phase A, rank subgroups in Phase B): sort by estimated solve cost n_local_grid^3
!>        (the per-atom Cholesky dominates; the n^2 assembly terms order the atoms the same
!>        way) and greedily give each atom to the least-loaded worker.
!> \param atom_list ...
!> \param n_atoms ...
!> \param n_local_grid_atom ...
!> \param n_workers ...
!> \param my_worker ...
!> \param my_atoms ...
!> \param max_over_mean ...
! **************************************************************************************************
   SUBROUTINE lpt_assign_atoms(atom_list, n_atoms, n_local_grid_atom, n_workers, my_worker, &
                               my_atoms, max_over_mean)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_list
      INTEGER, INTENT(IN)                                :: n_atoms
      INTEGER, DIMENSION(:), INTENT(IN)                  :: n_local_grid_atom
      INTEGER, INTENT(IN)                                :: n_workers, my_worker
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: my_atoms
      REAL(KIND=dp), INTENT(OUT)                         :: max_over_mean
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lpt_assign_atoms'
      INTEGER                                            :: handle

      INTEGER                                            :: i, iw, n_mine, w_min
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mine_tmp, perm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cost, load

      CALL timeset(routineN, handle)

      max_over_mean = 1.0_dp
      IF (n_atoms <= 0) THEN
         ALLOCATE (my_atoms(0))
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (cost(n_atoms), perm(n_atoms), mine_tmp(n_atoms), load(n_workers))
      DO i = 1, n_atoms
         cost(i) = REAL(n_local_grid_atom(atom_list(i)), dp)**3
      END DO
      CALL sort(cost, n_atoms, perm)   ! ascending; walk backwards for largest-first

      load(:) = 0.0_dp
      n_mine = 0
      DO i = n_atoms, 1, -1
         w_min = 1
         DO iw = 2, n_workers
            IF (load(iw) < load(w_min)) w_min = iw
         END DO
         load(w_min) = load(w_min) + cost(i)
         IF (w_min - 1 == my_worker) THEN
            n_mine = n_mine + 1
            mine_tmp(n_mine) = atom_list(perm(i))
         END IF
      END DO

      ALLOCATE (my_atoms(n_mine))
      my_atoms(:) = mine_tmp(1:n_mine)
      IF (SUM(load) > 0.0_dp) max_over_mean = MAXVAL(load)*REAL(n_workers, dp)/SUM(load)

      CALL timestop(handle)

   END SUBROUTINE lpt_assign_atoms

! **************************************************************************************************
!> \brief Computes the dense localized RHS d_lp(l,P) = Σ_{μν} Φ_μ(r_l)·Φ_ν(r_l)·(μν|P) for one
!>        RI atom P, OMP-threaded over (atom_j, atom_k) AO-pair blocks: per thread, build the 3c
!>        block, then contract grid-chunked pair densities into a private d_lp partial; partials
!>        are reduced into d_lp at the end.
!>        Pair screening is handled inside build_3c_integral_block_ctx via the `screened` output.
!> \param bs_env ...
!> \param ctx ...
!> \param phi_val ...
!> \param d_lp ...
!> \param n_grid_total ...
!> \param n_loc_ri ...
!> \param atom_P ...
!> \param max_ao_size ...
!> \param atom_j_mepos ...
!> \param atom_j_stride ...
! **************************************************************************************************
   SUBROUTINE compute_d_lp(bs_env, ctx, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, &
                           max_ao_size, atom_j_mepos, atom_j_stride)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_val
      INTEGER, INTENT(IN)                                :: n_grid_total, n_loc_ri
      REAL(KIND=dp), INTENT(INOUT)                       :: d_lp(n_grid_total, n_loc_ri)
      INTEGER, INTENT(IN)                                :: atom_P, max_ao_size, atom_j_mepos, &
                                                            atom_j_stride

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_d_lp'
      INTEGER, PARAMETER                                 :: grid_chunk = 1024

      INTEGER                                            :: atom_j, atom_k, c, handle, &
                                                            handle_dgemm, j, jk_idx, &
                                                            jsize, jstart, k, ksize, kstart, l, &
                                                            l0, ri
      LOGICAL                                            :: screened
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_prv, int_2d_prv, rho_chunk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c_prv
      TYPE(gw_3c_ws_type)                                :: ws

      CALL timeset(routineN, handle)

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, ctx, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, max_ao_size, &
      !$OMP        atom_j_mepos, atom_j_stride) &
      !$OMP PRIVATE(atom_j, atom_k, c, handle_dgemm, j, jk_idx, jsize, jstart, k, ksize, kstart, &
      !$OMP         l, l0, ri, screened, d_lp_prv, int_2d_prv, rho_chunk, int_3c_prv, ws)

      CALL gw_3c_ws_create(ws, ctx)
      ALLOCATE (int_3c_prv(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (int_2d_prv(max_ao_size*max_ao_size, n_loc_ri))
      ALLOCATE (rho_chunk(grid_chunk, max_ao_size*max_ao_size))
      ALLOCATE (d_lp_prv(n_grid_total, n_loc_ri))
      d_lp_prv(:, :) = 0.0_dp

      ! MPI-stride atom_j over the subgroup (atom_j_stride = 1 for the BLAS
      ! path, > 1 for the ScaLAPACK path). The OMP DO parallelizes the inner
      ! atom_k while the outer atom_j carries the MPI stride.
      !$OMP DO SCHEDULE(DYNAMIC)
      DO atom_j = atom_j_mepos + 1, SIZE(bs_env%i_ao_start_from_atom), atom_j_stride
         DO atom_k = 1, SIZE(bs_env%i_ao_start_from_atom)
            jstart = bs_env%i_ao_start_from_atom(atom_j)
            jsize = bs_env%i_ao_end_from_atom(atom_j) - jstart + 1
            kstart = bs_env%i_ao_start_from_atom(atom_k)
            ksize = bs_env%i_ao_end_from_atom(atom_k) - kstart + 1

            int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp

            ! Compute B_{μν,P} = (μν|P); ctx-internal triangle-inequality screening on
            ! kind_radius sets `screened=.TRUE.` for negligible triples.
            CALL build_3c_integral_block_ctx(int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri), &
                                             ctx, ws, atom_j=atom_j, atom_k=atom_k, atom_i=atom_P, &
                                             screened=screened)

            IF (screened) CYCLE

            ! Flatten 3D B_{μν, P} tensor to 2D B_{(μν), P} matrix for BLAS
            DO ri = 1, n_loc_ri
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     int_2d_prv(jk_idx, ri) = int_3c_prv(j, k, ri)
                  END DO
               END DO
            END DO

            ! Pair density ρ(l, μν) = Φ_μ(r_l) Φ_ν(r_l) in grid chunks, contracted on the fly:
            ! d_{l,P} += ρ(l, μν) B_{(μν),P}  (dgemm runs serially inside the parallel region)
            DO l0 = 1, n_grid_total, grid_chunk
               c = MIN(grid_chunk, n_grid_total - l0 + 1)
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     DO l = 1, c
                        rho_chunk(l, jk_idx) = phi_val(l0 + l - 1, jstart + j - 1)* &
                                               phi_val(l0 + l - 1, kstart + k - 1)
                     END DO
                  END DO
               END DO
               CALL timeset(routineN//"_dgemm", handle_dgemm)
               CALL dgemm("N", "N", c, n_loc_ri, jsize*ksize, &
                          1.0_dp, rho_chunk, grid_chunk, &
                          int_2d_prv, max_ao_size*max_ao_size, &
                          1.0_dp, d_lp_prv(l0, 1), n_grid_total)
               CALL timestop(handle_dgemm)
            END DO
         END DO
      END DO
      !$OMP END DO

      !$OMP CRITICAL (compute_d_lp_reduce)
      d_lp(1:n_grid_total, 1:n_loc_ri) = d_lp(1:n_grid_total, 1:n_loc_ri) + &
                                         d_lp_prv(1:n_grid_total, 1:n_loc_ri)
      !$OMP END CRITICAL (compute_d_lp_reduce)

      DEALLOCATE (int_3c_prv, int_2d_prv, rho_chunk, d_lp_prv)
      CALL gw_3c_ws_release(ws)

      !$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE compute_d_lp

! **************************************************************************************************
!> \brief Distributed pdpotrf/pdpotrs solve of D x = b for one atom of the
!>        RI-RS Z_lP build (Phase B, "big" atoms), called with a subgroup of
!>        cooperating ranks and an associated BLACS context.
!>        Each rank in the subgroup holds the (replicated) phi_local and the
!>        (replicated) RHS d_lp; it builds its own block-cyclic slice of the
!>        squared+Jacobi-scaled Gram matrix D via tiled DGEMM, factorizes via
!>        cp_fm_cholesky_decompose (UPLO='U'), and solves with
!>        cp_fm_cholesky_solve. The replicated d_lp is updated in place via
!>        cp_fm_get_submatrix.
!> \param phi_local ...
!> \param d_vec ...
!> \param d_lp ...
!> \param n_loc ...
!> \param n_ao ...
!> \param n_rhs ...
!> \param tikhonov ...
!> \param para_env_sub ...
!> \param blacs_env_sub ...
!> \param fm_struct_D ...
!> \param fm_struct_b ...
!> \param fm_D ...
!> \param fm_b ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE solve_D_lp_distributed(phi_local, d_vec, d_lp, n_loc, n_ao, n_rhs, &
                                     tikhonov, para_env_sub, blacs_env_sub, &
                                     fm_struct_D, fm_struct_b, fm_D, fm_b, info)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: d_vec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      INTEGER, INTENT(IN)                                :: n_loc, n_ao, n_rhs
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_D, fm_struct_b
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_D, fm_b
      INTEGER, INTENT(OUT)                               :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_D_lp_distributed'

      INTEGER                                            :: handle, i_loc, ig, j_loc, jg, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data

      CALL timeset(routineN, handle)
      info = 0

      NULLIFY (fm_struct_D, fm_struct_b)
      CALL cp_fm_struct_create(fm_struct_D, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_loc)
      CALL cp_fm_struct_create(fm_struct_b, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_rhs)
      CALL cp_fm_create(fm_D, fm_struct_D)
      CALL cp_fm_create(fm_b, fm_struct_b)

      ! ---- Build the local block-cyclic slice of fm_D --------------------
      ! Tiled DGEMM build: for each (row_tile x col_tile) sub-block we
      ! gather small phi strips, DGEMM into a small gram tile, then write
      ! the squared & d_vec-scaled result directly into fm_D%local_data.
      ! This avoids materialising the full (nrow_local x n_ao),
      ! (n_ao x ncol_local), and (nrow_local x ncol_local) buffers.
      CALL cp_fm_get_info(fm_D, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data)

      IF (nrow_local > 0 .AND. ncol_local > 0) THEN
         BLOCK
            INTEGER, PARAMETER :: ntile = 1024
            INTEGER :: ib, ie, jb, je, mb, kb, ti, tj, handle_dgemm
            REAL(KIND=dp), ALLOCATABLE :: gram_t(:, :), phi_cols_t(:, :), phi_rows_t(:, :)
            ALLOCATE (phi_rows_t(ntile, n_ao), phi_cols_t(n_ao, ntile), gram_t(ntile, ntile))
            DO ib = 1, nrow_local, ntile
               ie = MIN(ib + ntile - 1, nrow_local)
               mb = ie - ib + 1
               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(mb, n_ao, phi_rows_t, phi_local, row_indices, ib) &
               !$OMP PRIVATE(ti, j_loc) SCHEDULE(STATIC)
               DO j_loc = 1, n_ao
                  DO ti = 1, mb
                     phi_rows_t(ti, j_loc) = phi_local(row_indices(ib + ti - 1), j_loc)
                  END DO
               END DO
               !$OMP END PARALLEL DO
               DO jb = 1, ncol_local, ntile
                  je = MIN(jb + ntile - 1, ncol_local)
                  kb = je - jb + 1
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(kb, n_ao, phi_cols_t, phi_local, col_indices, jb) &
                  !$OMP PRIVATE(tj, i_loc) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     DO i_loc = 1, n_ao
                        phi_cols_t(i_loc, tj) = phi_local(col_indices(jb + tj - 1), i_loc)
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
                  CALL timeset(routineN//"_dgemm", handle_dgemm)
                  CALL dgemm('N', 'N', mb, kb, n_ao, &
                             1.0_dp, phi_rows_t, ntile, phi_cols_t, n_ao, &
                             0.0_dp, gram_t, ntile)
                  CALL timestop(handle_dgemm)
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(mb, kb, gram_t, d_vec, row_indices, col_indices, ib, jb) &
                  !$OMP SHARED(local_data, tikhonov) &
                  !$OMP PRIVATE(ti, tj, ig, jg) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     jg = col_indices(jb + tj - 1)
                     DO ti = 1, mb
                        ig = row_indices(ib + ti - 1)
                        local_data(ib + ti - 1, jb + tj - 1) = &
                           gram_t(ti, tj)*gram_t(ti, tj)*d_vec(ig)*d_vec(jg)
                        IF (ig == jg) THEN
                           local_data(ib + ti - 1, jb + tj - 1) = &
                              local_data(ib + ti - 1, jb + tj - 1) + tikhonov
                        END IF
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
               END DO
            END DO
            DEALLOCATE (phi_rows_t, phi_cols_t, gram_t)
         END BLOCK
      END IF

      ! Load the replicated d_lp into the block-cyclic fm_b
      CALL cp_fm_set_submatrix(fm_b, d_lp)

      ! pdpotrf (Cholesky factorisation; cp_fm_cholesky_decompose
      ! factors with UPLO='U', so pdpotrs must match)
      CALL cp_fm_cholesky_decompose(fm_D, n=n_loc, info_out=info)
      IF (info /= 0) THEN
         CPABORT("pdpotrf failed in solve_D_lp_distributed")
      END IF

      ! pdpotrs/dpotrs (solve in place on fm_b)
      CALL cp_fm_cholesky_solve(fm_D, fm_b, n=n_loc, info_out=info)
      IF (info /= 0) THEN
         CPABORT("pdpotrs failed in solve_D_lp_distributed")
      END IF

      ! Gather distributed solution back into the replicated d_lp
      CALL cp_fm_get_submatrix(fm_b, d_lp)

      CALL cp_fm_release(fm_D)
      CALL cp_fm_release(fm_b)
      CALL cp_fm_struct_release(fm_struct_D)
      CALL cp_fm_struct_release(fm_struct_b)

      CALL timestop(handle)

   END SUBROUTINE solve_D_lp_distributed

! **************************************************************************************************
!> \brief Computes the polarizability matrix in the RI basis for every
!>        imaginary-time point:
!>          G^occ/vir_μν(i|τ|) = Σ_n C_μn e^(-|(ϵ_n-ϵ_F)τ|) C_νn      (AO x AO, build_G_ao)
!>          χ_ll'(iτ) = [Σ_μν Φ_μ(r_l) G^occ_μν Φ_ν(r_l')] ∘ [Σ_μν Φ_μ(r_l) G^vir_μν Φ_ν(r_l')]
!>          χ_PQ(iτ)  = g_s Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q                (g_s = spin degeneracy)
!>        The grid (l) index is streamed in panels (contract_grid_panels).
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
         WRITE (bs_env%unit_nr, '(T2,A,I44)') 'Number of batches for χ, Σ matrices', n_panels
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
         IF (i_t == 1 .AND. bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(A)') ' '
            WRITE (bs_env%unit_nr, '(T2,A,F17.2,A)') &
               'Occupation of non-zero elements of G(l,l''), χ(l,l''), W(l,l'')', &
               grid_occ*100.0_dp, ' %'
            WRITE (bs_env%unit_nr, '(A)') ' '
            CALL m_flush(bs_env%unit_nr)
         END IF
         IF (i_t == 1) then
            CALL print_matrix_occupation(mat_chi_Gamma_tau(i_t)%matrix, 'χ(P,Q)', &
                                         bs_env%para_env, bs_env%unit_nr)
         end if

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
!> \brief Prints the non-zero occupation percentage of a DBCSR matrix on one line.
!> \param matrix ...
!> \param label ...
!> \param para_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_matrix_occupation(matrix, label, para_env, unit_nr)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_matrix_occupation'

      INTEGER                                            :: handle

      REAL(KIND=dp)                                      :: frac_2p31, max_loc, occ

      CALL timeset(routineN, handle)

      occ = dbcsr_get_occupation(matrix)
      max_loc = REAL(dbcsr_get_data_size(matrix), dp)
      CALL para_env%max(max_loc)

      IF (unit_nr > 0) THEN
         frac_2p31 = max_loc/REAL(dbcsr_msg_elem_limit, dp)
         WRITE (unit_nr, '(A)') ' '
         WRITE (unit_nr, '(T2,A,F36.2,A)') &
            'Occupation of non-zero elements of '//TRIM(label), occ*100.0_dp, ' %'
         IF (frac_2p31 > 0.5_dp) then
            WRITE (unit_nr, '(T4,A)') '*** WARNING: max/rank approaching 2^31 -- DBCSR overflow risk ***'
         end if
         CALL m_flush(unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_matrix_occupation

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
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mask_grid_blocks_near_panel'

      INTEGER                                            :: handle
      INTEGER                                            :: c, k
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

      INTEGER                                            :: handle
      INTEGER                                            :: c, ib, n_used
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
            IF (SUM((centroids(:, ib) - centroids(:, used_idx(c)))**2) <= cutoff2) then
               nze_tmpl = nze_tmpl + INT(r_blk_sizes(ib), int_8)*INT(r_blk_sizes(used_idx(c)), int_8)
            end if
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
! **************************************************************************************************
   SUBROUTINE plan_grid_panels(r_blk_sizes, panel_size, min_dim, pan_first, pan_last, &
                               centroids, cutoff, n_RI, n_ao, n_procs, mem_budget_GB, &
                               honor_exact, unsafe)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      INTEGER, INTENT(IN)                                :: panel_size, min_dim
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      INTEGER, INTENT(IN), OPTIONAL                      :: n_RI, n_ao, n_procs
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: mem_budget_GB
      LOGICAL, INTENT(IN), OPTIONAL                      :: honor_exact
      LOGICAL, INTENT(OUT), OPTIONAL                     :: unsafe
      CHARACTER(LEN=*), PARAMETER :: routineN = 'plan_grid_panels'

      INTEGER                                            :: blk0, blk1, ib, n_grid_blocks, &
                                                            n_grid_total, n_panels, rows_acc, &
                                                            target, width
      INTEGER(KIND=int_8)                                :: msg, nze_tmpl, side
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp_first, tmp_last
      LOGICAL                                            :: fits, my_honor_exact, my_unsafe, &
                                                            use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: used
      REAL(KIND=dp)                                      :: f_near, mem_GB
      INTEGER                                            :: handle

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
         target = panel_size
         DO
            rows_acc = 0
            blk1 = blk0
            DO ib = blk0, n_grid_blocks
               rows_acc = rows_acc + r_blk_sizes(ib)
               blk1 = ib
               IF (rows_acc >= target) EXIT
            END DO
            IF (.NOT. use_cutoff .OR. blk1 == blk0) EXIT
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)
            width = SUM(r_blk_sizes, MASK=used)
            CALL panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)
            f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
            side = INT(REAL(rows_acc, dp)*f_near*REAL(MAX(n_RI, n_ao), dp), int_8)
            msg = MAX(nze_tmpl, side)/INT(MAX(min_dim, 1), int_8)
            fits = (msg <= dbcsr_msg_elem_limit/4)
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
            target = MAX(1, MIN(target, rows_acc)/2)
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

      INTEGER                                            :: min_dim, n_grid_total, n_panels_req, &
                                                            npcols, nprows, panel_size, safe_max
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: honor_exact, panels_unsafe, use_cutoff
      REAL(KIND=dp)                                      :: mem_avail_GB, mem_budget_GB
      TYPE(dbcsr_distribution_type)                      :: dist
      INTEGER                                            :: handle

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
      ! (see dbcsr_msg_elem_limit): requiring the worst-rank message to stay under
      ! 0.5 * HUGE(int_4) gives the safe height P_safe = 0.5 * HUGE(int_4) * min_dim / n_grid.
      IF (use_cutoff) THEN
         safe_max = n_grid_total
      ELSE
         safe_max = INT(0.5_dp*REAL(dbcsr_msg_elem_limit, dp)*REAL(min_dim, dp)/ &
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
         CALL plan_grid_panels(r_blk_sizes, panel_size, min_dim, pan_first, pan_last, &
                               centroids=bs_env%ri_rs%chunk_centroids, &
                               cutoff=bs_env%ri_rs%cutoff_radius_v_w, &
                               n_RI=bs_env%n_RI, n_ao=bs_env%n_ao, &
                               n_procs=bs_env%para_env%num_pe, mem_budget_GB=mem_budget_GB, &
                               honor_exact=honor_exact, unsafe=panels_unsafe)
      ELSE
         CALL plan_grid_panels(r_blk_sizes, panel_size, min_dim, pan_first, pan_last)
      END IF

      IF (honor_exact .AND. panels_unsafe .AND. bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,I0,A)') &
            '*** WARNING: N_PANELS = ', n_panels_req, ' used as requested, but one or more '// &
            'panels exceed the DBCSR 32-bit message / memory-budget safety limit. The run may '// &
            'abort or swap; increase N_PANELS if it does. ***'
      END IF

      ALLOCATE (bs_env%ri_rs%pan_first, SOURCE=pan_first)
      ALLOCATE (bs_env%ri_rs%pan_last, SOURCE=pan_last)

      CALL timestop(handle)

   END SUBROUTINE resolve_grid_panels

! **************************************************************************************************
!> \brief Available memory per MPI process (GB). /proc/meminfo reports node-wide memory, so
!>        every rank on a node reads the SAME MemLikelyFree; the per-process share is
!>        node_free / ranks_per_node (ranks grouped by a hostname hash exchanged via allgather).
!>        Returns the MIN across all ranks (most-constrained node)
!> \param bs_env ...
!> \param mem_avail_GB ...
! **************************************************************************************************
   SUBROUTINE ri_rs_mem_avail_per_proc_GB(bs_env, mem_avail_GB)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: mem_avail_GB
      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_rs_mem_avail_per_proc_GB'

      CHARACTER(LEN=default_string_length)               :: hostname
      INTEGER                                            :: host_hash, ic, n_procs, ranks_per_node
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: all_host_hashes
      INTEGER(KIND=int_8)                                :: h8, mem_buffers, mem_cached, mem_free, &
                                                            mem_likely_free, mem_sreclaimable, &
                                                            mem_slab, mem_total
      INTEGER                                            :: handle

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

      INTEGER                                            :: handle, iatom, ipan, l, &
                                                            max_n_local_grid, &
                                                            n_ao_total, n_grid_total, &
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

      ! Z_lP SOLVE peak (compute_coeff_Z_lP). For the atom P with the largest integration
      ! sphere, one rank holds simultaneously:
      !   D_local  : n_local_grid x n_local_grid  (dense Gram, BLAS path only; O(n_local_grid^2))
      !   phi_local: n_local_grid x n_ao_total
      !   d_lp     : n_local_grid x n_loc_ri, replicated once + one private copy per OMP thread
      ! n_local_grid = # grid points within cutoff_ri(P) = CUTOFF_RADIUS_RL_RI (if > 0) else
      ! r_c(RI metric) + r_RI(P). This is NOT evenly distributed: n_local_grid depends on the
      ! local density of atoms/grid, so the rank owning the densest atom peaks well above the
      ! average. We report the worst-case (max over atoms) as a per-rank upper bound.
      CALL get_qs_env(qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      n_ao_total = bs_env%i_ao_end_from_atom(natom)

      max_n_local_grid = 0
      n_loc_ri_max = 0
      DO iatom = 1, natom
         IF (bs_env%ri_rs%cutoff_radius_ri_rs > 0.0_dp) THEN
            cutoff_ri = bs_env%ri_rs%cutoff_radius_ri_rs
         ELSE
            cutoff_ri = bs_env%ri_metric%cutoff_radius + bs_env%ri_rs%radius_ri_per_atom(iatom)
         END IF
         pos_P(:) = particle_set(iatom)%r(:)
         n_local_grid = 0
         DO l = 1, n_grid_total
            IF (SUM((bs_env%ri_rs%grid_points(1:3, l) - pos_P(1:3))**2) <= cutoff_ri**2) then
               n_local_grid = n_local_grid + 1
            end if
         END DO
         max_n_local_grid = MAX(max_n_local_grid, n_local_grid)
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
      mem_phi_local_GB = REAL(max_n_local_grid, dp)*REAL(n_ao_total, dp)*8.0_dp*1.0E-9_dp
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
         WRITE (bs_env%unit_nr, '(T6,A,I21,A,F10.2,A)') &
            'worst-case n_local_grid', max_n_local_grid, ' points  (D_local', mem_D_local_GB, ' GB)'
         WRITE (bs_env%unit_nr, '(A)') ' '

         IF (mem_avail_GB > 0.0_dp .AND. mem_Z_lP_GB > mem_avail_GB) THEN
            WRITE (bs_env%unit_nr, '(T2,A)') &
               '*** WARNING: Estimated Z_lP memory exceeds available memory per process ***'
            WRITE (bs_env%unit_nr, '(T4,A,F6.2,A,F6.2,A)') &
               'Z_lP upper bound: ', mem_Z_lP_GB, ' GB  >  available: ', mem_avail_GB, ' GB'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'Z_lP (n_grid × n_RI) is distributed across all MPI ranks. To reduce per-rank'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'memory: add more nodes, use fewer MPI ranks per node, or increase'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'N_PROCS_PER_ATOM_Z_LP to distribute each atom block via ScaLAPACK (reduces'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'per-rank memory by ~1/G where G = N_PROCS_PER_ATOM_Z_LP).'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
         IF (mem_avail_GB > 0.0_dp .AND. mem_panels_GB > mem_avail_GB) THEN
            WRITE (bs_env%unit_nr, '(T2,A)') &
               '*** WARNING: Estimated χ/W/Σ panel memory exceeds available memory per process ***'
            WRITE (bs_env%unit_nr, '(T4,A,F6.2,A,F6.2,A)') &
               'Panel peak estimate: ', mem_panels_GB, ' GB  >  available: ', mem_avail_GB, ' GB'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'Panel memory scales as ~3×panel_size×n_grid / n_procs. Options:'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '  - More nodes or fewer MPI ranks per node (increases n_procs, reduces share)'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '  - Increase N_PANELS (more, smaller panels → less peak memory per step)'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
         IF (mem_avail_GB > 0.0_dp .AND. mem_Zlp_peak_GB > mem_avail_GB) THEN
            WRITE (bs_env%unit_nr, '(T2,A)') &
               '*** WARNING: Estimated Z_lP solve peak exceeds available memory per process ***'
            WRITE (bs_env%unit_nr, '(T4,A,F8.2,A,F8.2,A)') &
               'Z_lP solve peak: ', mem_Zlp_peak_GB, ' GB  >  available: ', mem_avail_GB, ' GB'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'The per-atom Gram matrix D_local(n_local_grid, n_local_grid) dominates and scales'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'as n_local_grid^2. It is NOT balanced across ranks (the rank owning the atom with'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               'the largest integration sphere peaks well above the average). Options:'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '  - Increase N_PROCS_PER_ATOM_Z_LP=G to distribute D_local block-cyclic via'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '    ScaLAPACK (reduces the D_local term by ~1/G; no accuracy loss)'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '  - Reduce the RI-RS sphere cutoff (CUTOFF_RADIUS_RL_RI): D_local shrinks as'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '    n_local_grid^2, but this trades accuracy'
            WRITE (bs_env%unit_nr, '(T4,A)') &
               '  - Fewer MPI ranks per node (more memory per rank for the peak atom)'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
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
      TYPE(dbcsr_type), INTENT(IN)                      :: L_pan, L_full
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)        :: centroids
      REAL(KIND=dp), INTENT(IN)                         :: cutoff
      INTEGER, INTENT(IN)                               :: blk0
      TYPE(dbcsr_type), INTENT(OUT)                     :: A_template
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL       :: col_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_geo_template_panel'

      INTEGER                                           :: c, cg, cs, handle, my_pcol, my_prow, &
                                                           n_grid_blks, n_pan_blks, npcols, &
                                                           nprows, r, rs
      INTEGER, DIMENSION(:), POINTER                    :: grid_blk_sizes, pan_blk_sizes
      REAL(KIND=dp)                                     :: cutoff2
      REAL(KIND=dp), ALLOCATABLE                        :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                     :: dist

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

      INTEGER                                            :: ib, jb, npb
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_full, col_dist_full, &
                                                            row_blk_full, row_blk_pan, &
                                                            row_dist_full, row_dist_pan
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_pan
      TYPE(dbcsr_iterator_type)                          :: iter
      INTEGER                                            :: handle

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

      INTEGER                                            :: ib, jb, nblkcols
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iused
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter
      INTEGER                                            :: handle

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

      INTEGER                                            :: ib, jb, n_blk, n_sub, r
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inv_map
      INTEGER, DIMENSION(:), POINTER                     :: blk_full, blk_sub, col_blk_full, &
                                                            col_dist_full, dist_full_1d, &
                                                            dist_sub_1d, row_blk_full, &
                                                            row_dist_full
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_sub
      TYPE(dbcsr_iterator_type)                          :: iter
      INTEGER                                            :: handle

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
!> \brief Pre-seeds a square per-atom-blocked DBCSR matrix (G_munu, D_munu, V_PQ, W_PQ) with zero
!>        blocks only for atom pairs within radius, for use with copy_fm_to_dbcsr(keep_sparsity=T):
!>        the CUTOFF_RADIUS_G_W operator truncation.
!> \param matrix ...
!> \param centers ...
!> \param radius ...
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

      INTEGER                                            :: i, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_l, col_blk_r, out_col_blk, &
                                                            out_col_dist, out_row_blk, out_row_dist, &
                                                            row_blk_l, row_blk_r
      TYPE(dbcsr_distribution_type)                      :: dist_l, dist_out
      INTEGER                                            :: handle

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
      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_G_ao'

      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      TYPE(cp_fm_type), POINTER                          :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao
      INTEGER                                            :: handle

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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: grid_occupation

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, ipan, &
                                                            n_grid_total, ncols_pan, nrows_pan
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
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_B, 0.0_dp, B_pan, &
                                filter_eps=eps, retain_sparsity=.TRUE.)
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
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_G_occ_ao, mat_G_vir_ao, mat_W_aux
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_Sigma_neg, mat_Sigma_pos
      REAL(KIND=dp), INTENT(IN)                          :: eps
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pan_first, pan_last
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_grid_panels_sigma_c'

      INTEGER                                            :: blk0, blk1, handle, ipan
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap
      LOGICAL                                            :: my_keep_sparsity, use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, used_ao, used_ri
      TYPE(dbcsr_type)                                   :: A_occ, A_vir, C_pan, G_occ_sub, &
                                                            G_vir_sub, phi_pan, phi_panC, tmp2, &
                                                            tmpA, tmpB, W_pan, W_sub, Z_pan, &
                                                            Z_panC
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
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_Z, 0.0_dp, W_pan, &
                                filter_eps=eps, retain_sparsity=.TRUE.)
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

      ! M^-1*W^MIC(iτ)*M^-1
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_W_time)

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
!> \brief Computes the static RI-basis Coulomb operators entering the dielectric function:
!>        the bare Coulomb matrix V_PQ(k=0), its Cholesky/matrix square root V^0.5, and
!>        M^-1 V^0.5 with M the RI metric (2c integrals of the RI_METRIC operator).
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
! **************************************************************************************************
   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle, info, n_RI, ndep
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      CALL cp_fm_create(fm_work, fm_V%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1
      ! -----------------------------------------------------------------------
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      ALLOCATE (mat_V_kp(1:1, 1:2))
      NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
      ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)

      CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
      CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)

      ! Dummy imaginary part just to satisfy the routine
      CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
      CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                      particle_set, qs_kind_set, atomic_kind_set, &
                                      bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)

      ! Copy real part to fm_V
      CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

      CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
      CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
      DEALLOCATE (mat_V_kp)

      ! -----------------------------------------------------------------------
      ! 2. Get RI-Metric Matrix M(k=0)
      ! -----------------------------------------------------------------------
      CALL RI_2c_integral_mat(qs_env, fm_M, fm_V, n_RI, bs_env%ri_metric, &
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)

      ! -----------------------------------------------------------------------
      ! 3. M -> M^-1
      ! -----------------------------------------------------------------------
      CALL cp_fm_cholesky_decompose(fm_M(1, 1), info_out=info)
      IF (info == 0) THEN
         CALL cp_fm_cholesky_invert(fm_M(1, 1))
         CALL cp_fm_uplo_to_full(fm_M(1, 1), fm_work)
      ELSE
         ! Fallback if Cholesky fails due to conditioning
         CALL cp_fm_power(fm_M(1, 1), fm_work, -1.0_dp, bs_env%eps_eigval_mat_RI, ndep)
         CALL cp_fm_to_fm(fm_work, fm_M(1, 1))
      END IF

      ! -----------------------------------------------------------------------
      ! 4. V -> V^0.5
      ! -----------------------------------------------------------------------
      CALL cp_fm_to_fm(fm_V, fm_V_sqrt)
      CALL cp_fm_cholesky_decompose(fm_V_sqrt, info_out=info)
      IF (info == 0) THEN
         CALL clean_lower_part(fm_V_sqrt)
      ELSE
         CALL cp_fm_power(fm_V, fm_V_sqrt, 0.5_dp, bs_env%eps_eigval_mat_RI, ndep)
      END IF

      ! -----------------------------------------------------------------------
      ! 5. M^-1 * V^0.5
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", n_RI, n_RI, n_RI, 1.0_dp, fm_M(1, 1), fm_V_sqrt, &
                         0.0_dp, fm_Minv_Vsqrt)

      CALL cp_fm_release(fm_M)
      CALL cp_fm_release(fm_work)

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

      INTEGER                                            :: handle, info, n_RI, ndep
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

      ! a) Cholesky decomposition of ε(iω_j)
      CALL cp_fm_cholesky_decompose(fm_eps_freq_j, info_out=info)

      ! b) Inversion
      IF (info == 0) THEN
         CALL cp_fm_cholesky_invert(fm_eps_freq_j)
         CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)
      ELSE
         ! Fallback to expensive diagonalization if Cholesky fails
         CALL cp_fm_power(fm_eps_freq_j, fm_work, -1.0_dp, bs_env%eps_eigval_mat_RI, ndep)
         CALL cp_fm_to_fm(fm_work, fm_eps_freq_j)
      END IF

      ! c) ε^-1(iω_j) - Id
      CALL fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! d) work = (ε^-1(iω_j) - Id) * V^0.5
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_eps_freq_j, fm_V_sqrt, &
                         0.0_dp, fm_work)

      ! e) W(iw) = V^0.5^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_V_sqrt, fm_work, &
                         0.0_dp, fm_W_freq_j)

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
!> \brief Zeroes out the strictly lower triangular part of a real matrix
!> \param fm_mat ...
! **************************************************************************************************
   SUBROUTINE clean_lower_part(fm_mat)
      TYPE(cp_fm_type)                                   :: fm_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'clean_lower_part'

      INTEGER                                            :: handle, i_row, j_col, j_global, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            IF (j_global < row_indices(i_row)) fm_mat%local_data(i_row, j_col) = 0.0_dp
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE clean_lower_part

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

      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)

      ! M^-1 * V^tr * M^-1 directly modifies fm_Vtr_Gamma(:, 1)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      ! Optional CUTOFF_RADIUS_G_W operator truncation + filter (see build_G_ao)
      IF (bs_env%ri_rs%cutoff_radius_g_w > 0.0_dp .AND. &
          ALLOCATED(bs_env%ri_rs%atom_centers)) THEN
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

         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
         IF (bs_env%ri_rs%cutoff_radius_g_w > 0.0_dp .AND. &
             ALLOCATED(bs_env%ri_rs%atom_centers)) THEN
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

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      ! -------------------------------------------------------------------------
      ! 3. FINALIZE AND CLEANUP
      ! -------------------------------------------------------------------------
      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL cp_fm_release(fm_W_time)

      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL release_square_topology(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL delete_unnecessary_files(bs_env)
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'setup_square_topology'

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

END MODULE gw_non_periodic_ri_rs
