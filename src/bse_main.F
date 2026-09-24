!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main routines for GW + Bethe-Salpeter for computing electronic excitations: a GW path
!>        prepares the BSE environment, bse_solve runs the solver on it
!> \par History
!>      04.2024 created [Maximilian Graml]
!>      09.2026 split into prepare_bse_env and bse_solve [Maximilian Graml]
! **************************************************************************************************

MODULE bse_main

   USE bse_full_diag,                   ONLY: create_A_and_B,&
                                              create_hermitian_form_of_ABBA,&
                                              diagonalize_A,&
                                              diagonalize_C
   USE bse_iterative,                   ONLY: solve_bse_iteratively
   USE bse_matvec,                      ONLY: mem_fraction
   USE bse_print,                       ONLY: print_BSE_start_flag
   USE bse_types,                       ONLY: bse_env_type
   USE bse_util,                        ONLY: adapt_BSE_input_params,&
                                              bse_window_from_dft,&
                                              estimate_BSE_resources,&
                                              fulldiag_memory_estimate_GB,&
                                              get_bse_spin_block_layout,&
                                              screen_slabs,&
                                              truncate_BSE_matrices
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE input_constants,                 ONLY: &
        bse_abba, bse_both, bse_fulldiag, bse_iterdiag, bse_memcheck_abort, bse_memcheck_clamp, &
        bse_memcheck_off, bse_memcheck_warn, bse_tda, use_mom_ref_coac
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE mo_window,                       ONLY: mo_window_type
   USE moments_utils,                   ONLY: get_reference_point
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_main'

   PUBLIC :: prepare_bse_env, bse_solve, get_multipoles_ao, bse_fulldiag_memory_check

CONTAINS

! **************************************************************************************************
!> \brief Prepares the BSE environment from the old GW's slabs: run flags, the window on the DFT
!>        axis, the window-sized slabs with both energy sets, [1+Q(0)]^-1 once, the AO multipoles
!>        and the QS pointers. Starts no solver.
!> \param bse_env settings read; state filled here
!> \param fm_S_ia_full B^P_ia of the GW levels per spin, the RPA copy without the energy factor
!> \param fm_S_ij_full B^P_ij per spin
!> \param fm_S_ab_full B^P_ab per spin
!> \param fm_Q the static polarizability Q(0), untouched; the caller releases it afterwards
!> \param Eigenval quasiparticle energies (level, 1, spin)
!> \param Eigenval_scf DFT energies (level, 1, spin)
!> \param homo occupied MOs per spin
!> \param virtual virtual MOs per spin
!> \param dimen_RI RI dimension of the slabs
!> \param dimen_RI_red reduced RI dimension of Q(0)
!> \param gw_corr_lev_occ GW-corrected occupied levels; the window may not exceed them
!> \param bse_lev_virt GW-corrected virtual levels per spin
!> \param qs_env source of the pointers the solver reads and of the AO multipoles
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE prepare_bse_env(bse_env, fm_S_ia_full, fm_S_ij_full, fm_S_ab_full, fm_Q, &
                              Eigenval, Eigenval_scf, homo, virtual, dimen_RI, dimen_RI_red, &
                              gw_corr_lev_occ, bse_lev_virt, qs_env, unit_nr)

      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_S_ia_full, fm_S_ij_full, fm_S_ab_full
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Q
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval, Eigenval_scf
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red, gw_corr_lev_occ
      INTEGER, DIMENSION(:), INTENT(IN)                  :: bse_lev_virt
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prepare_bse_env'

      INTEGER                                            :: handle, ispin, n_window, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_mo, n_ov, offsets
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_dft_spin, eigenval_qp_spin
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_window_type)                               :: window

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      bse_env%nspins = nspins
      bse_env%dimen_RI = dimen_RI
      bse_env%dimen_RI_red = dimen_RI_red
      ! the GW-corrected level counts, the ceiling of the window
      bse_env%gw_corr_lev_occ = gw_corr_lev_occ
      ALLOCATE (bse_env%gw_corr_lev_virt(nspins))
      bse_env%gw_corr_lev_virt(:) = bse_lev_virt(:)

      ! run flags from the settings
      bse_env%do_iterdiag = bse_env%bse_diag_method == bse_iterdiag
      bse_env%do_fulldiag = bse_env%bse_diag_method == bse_fulldiag
      bse_env%do_tda = bse_env%flag_tda == bse_tda .OR. bse_env%flag_tda == bse_both
      bse_env%do_abba = bse_env%flag_tda == bse_abba .OR. bse_env%flag_tda == bse_both

      CALL print_BSE_start_flag(bse_env%do_tda, bse_env%do_abba, unit_nr)

      ! the debug flag follows the debug print level
      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level == debug_print_level) THEN
         bse_env%bse_debug_print = .TRUE.
      END IF

      ! what the solver reads from QS and the GW (pointers, none owned); the para_env is the RPA's,
      ! the one the slabs live on
      bse_env%para_env => fm_S_ia_full(1)%matrix_struct%para_env
      CALL get_qs_env(qs_env, subsys=bse_env%subsys, pw_env=bse_env%pw_env, &
                      dft_control=bse_env%dft_control, mos=bse_env%mos, &
                      exstate_env=bse_env%exstate_env)
      CALL bse_env%para_env%sync()

      IF (nspins > 1) THEN
         IF (bse_env%do_iterdiag) THEN
            CPABORT("Iterative BSE is not yet available for open-shell references")
         END IF
         CALL cp_warn(__LOCATION__, &
                      "Open-shell (UKS/LSD) BSE is a recent addition and has not been "// &
                      "extensively validated. Verify results carefully before using them "// &
                      "for production calculations.")
      END IF

      ! one absolute window on the DFT axis, then the window-sized slabs and both energy sets per spin
      ALLOCATE (n_mo(nspins))
      n_mo(:) = homo(:) + virtual(:)
      CALL bse_window_from_dft(Eigenval_scf(:, 1, :), n_mo, homo, bse_env%bse_cutoff_occ, &
                               bse_env%bse_cutoff_empty, window, warn=.TRUE.)
      DEALLOCATE (n_mo)

      ALLOCATE (bse_env%homo_red(nspins), bse_env%virt_red(nspins), bse_env%homo_full(nspins))
      ALLOCATE (bse_env%fm_S_ia(nspins), bse_env%fm_S_ij(nspins), bse_env%fm_S_ab(nspins))
      bse_env%homo_full(:) = homo(:)
      DO ispin = 1, nspins
         ! the cutoff summary is printed for a closed-shell reference only
         CALL truncate_BSE_matrices(fm_S_ia_full(ispin), fm_S_ij_full(ispin), fm_S_ab_full(ispin), &
                                    Eigenval_scf(:, 1, ispin), Eigenval(:, 1, ispin), &
                                    eigenval_dft_spin, eigenval_qp_spin, homo(ispin), virtual(ispin), &
                                    dimen_RI, unit_nr, ispin, bse_env, window, nspins == 1)
         ! the window is one absolute range, so every spin has the same number of levels
         IF (ispin == 1) THEN
            n_window = SIZE(eigenval_qp_spin)
            ALLOCATE (bse_env%eigenval_dft(n_window, nspins), bse_env%eigenval_qp(n_window, nspins))
         END IF
         bse_env%eigenval_dft(:, ispin) = eigenval_dft_spin(:)
         bse_env%eigenval_qp(:, ispin) = eigenval_qp_spin(:)
         DEALLOCATE (eigenval_dft_spin, eigenval_qp_spin)
      END DO
      ALLOCATE (n_ov(nspins), offsets(nspins))
      CALL get_bse_spin_block_layout(bse_env%homo_red, bse_env%virt_red, n_ov, offsets, bse_env%n_ov_joint)
      DEALLOCATE (n_ov, offsets)

      CALL invert_static_screening(fm_Q, dimen_RI_red, bse_env%fm_eps_inv)

      ! the AO multipoles for the transition dipoles and, when asked for, the exciton descriptors
      CALL get_multipoles_ao(qs_env, 1, bse_env%matrix_dipole_ao, bse_env%multipole_ref_point)
      IF (bse_env%num_print_exc_descr /= 0) THEN
         CALL get_multipoles_ao(qs_env, 2, bse_env%matrix_quadpole_ao, bse_env%multipole_ref_point)
      END IF

      CALL timestop(handle)

   END SUBROUTINE prepare_bse_env

! **************************************************************************************************
!> \brief The inverse dielectric matrix in the RI basis, [1+Q(0)]^-1_PQ on the leading dimen_RI_red
!>        block, P, Q ≤ dimen_RI_red, by Cholesky decomposition and inversion of a copy of Q(0)
!> \param fm_Q the static polarizability Q_PQ(0), untouched
!> \param dimen_RI_red reduced RI dimension of Q(0)
!> \param fm_eps_inv [1+Q(0)]^-1_PQ, created here on fm_Q's structure
! **************************************************************************************************
   SUBROUTINE invert_static_screening(fm_Q, dimen_RI_red, fm_eps_inv)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Q
      INTEGER, INTENT(IN)                                :: dimen_RI_red
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_eps_inv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_static_screening'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      ! a copy: the caller keeps Q(0), and the open-shell path inverts once for both spins
      CALL cp_fm_create(fm_eps_inv, fm_Q%matrix_struct)
      CALL cp_fm_to_fm(fm_Q, fm_eps_inv)
      CALL cp_fm_create(fm_work, fm_Q%matrix_struct)
      CALL cp_fm_set_all(fm_work, 0.0_dp)

      ! δ_PQ + Q_PQ(0) on the leading dimen_RI_red block
      CALL cp_fm_get_info(matrix=fm_eps_inv, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI_red) THEN
               fm_eps_inv%local_data(iiB, jjB) = fm_eps_inv%local_data(iiB, jjB) + 1.0_dp
            END IF
         END DO
      END DO

      ! [1+Q(0)]^-1_PQ by Cholesky, symmetrised from the stored triangle
      CALL cp_fm_cholesky_decompose(matrix=fm_eps_inv, n=dimen_RI_red, info_out=info_chol)
      IF (info_chol /= 0) THEN
         CALL cp_abort(__LOCATION__, 'Cholesky decomposition failed for static polarization in BSE')
      END IF
      CALL cp_fm_cholesky_invert(fm_eps_inv)
      CALL cp_fm_uplo_to_full(fm_eps_inv, fm_work)

      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE invert_static_screening

! **************************************************************************************************
!> \brief Solves the BSE on a prepared environment: input normalisation, the screened slabs from
!>        [1+Q(0)]^-1, then the FULLDIAG or ITERDIAG solver on the energies USE_KS_ENERGIES selects
!> \param bse_env prepared by a GW path; the screened slabs are formed here
!> \param mo_coeff MO coefficients of the reference, per spin, for the properties
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE bse_solve(bse_env, mo_coeff, unit_nr)

      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_solve'

      INTEGER                                            :: handle, ispin, n_window, nspins
      REAL(KIND=dp)                                      :: diag_runtime_est
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_joint
      TYPE(cp_fm_type)                                   :: fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                            fm_inv_sqrt_A_minus_B, &
                                                            fm_sqrt_A_minus_B
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      nspins = bse_env%nspins

      ! normalisation touches no slab; n_ov_joint x 1 is the transition count on either spin path
      CALL adapt_BSE_input_params(bse_env%n_ov_joint, 1, unit_nr, bse_env)

      ! the energies the solver sees, flat joint layout [spin 1 window | spin 2 window]
      n_window = SIZE(bse_env%eigenval_qp, 1)
      ALLOCATE (eigenval_joint(n_window*nspins))
      DO ispin = 1, nspins
         IF (bse_env%use_ks_energies) THEN
            eigenval_joint((ispin - 1)*n_window + 1:ispin*n_window) = bse_env%eigenval_dft(:, ispin)
         ELSE
            eigenval_joint((ispin - 1)*n_window + 1:ispin*n_window) = bse_env%eigenval_qp(:, ispin)
         END IF
      END DO

      ! \bar{B}^P_ij and \bar{B}^P_ia per spin from the one [1+Q(0)]^-1
      ALLOCATE (bse_env%fm_S_bar_ia(nspins), bse_env%fm_S_bar_ij(nspins))
      DO ispin = 1, nspins
         CALL screen_slabs(bse_env%fm_eps_inv, bse_env%fm_S_ij(ispin), bse_env%fm_S_ia(ispin), &
                           bse_env%dimen_RI_red, bse_env%homo_red(ispin), bse_env%virt_red(ispin), &
                           bse_env%fm_S_bar_ij(ispin), bse_env%fm_S_bar_ia(ispin))
      END DO
      CALL cp_fm_release(bse_env%fm_eps_inv)

      IF (bse_env%do_iterdiag) THEN
         CALL solve_bse_iteratively(eigenval_joint, &
                                    bse_env%homo_red, bse_env%virt_red, bse_env%homo_full, &
                                    bse_env%dimen_RI, bse_env%do_tda, bse_env%do_abba, &
                                    bse_env, mo_coeff, unit_nr)
      END IF

      IF (bse_env%do_fulldiag) THEN
         CALL estimate_BSE_resources(bse_env%n_ov_joint, unit_nr, bse_env%do_abba, &
                                     bse_env%para_env, diag_runtime_est)
         CALL create_A_and_B(bse_env%fm_S_ia, bse_env%fm_S_ij, bse_env%fm_S_ab, &
                             bse_env%fm_S_bar_ia, bse_env%fm_S_bar_ij, fm_A_BSE, fm_B_BSE, &
                             bse_env%do_abba, eigenval_joint, unit_nr, bse_env%homo_red, &
                             bse_env%virt_red, bse_env%dimen_RI, bse_env)
         IF (bse_env%do_abba) THEN
            CALL create_hermitian_form_of_ABBA(fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                               fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                               unit_nr, bse_env, diag_runtime_est)
         END IF
         CALL cp_fm_release(fm_B_BSE)
         tddfpt_control => bse_env%dft_control%tddfpt2_control
         ! homo_full: per-spin absolute MO labels for the amplitude table, N_e = SUM(homo_full) for TRK
         IF (bse_env%do_tda .AND. (.NOT. tddfpt_control%do_bse)) THEN
            CALL diagonalize_A(fm_A_BSE, bse_env%homo_red, bse_env%virt_red, bse_env%homo_full, &
                               unit_nr, diag_runtime_est, bse_env, mo_coeff)
         END IF
         CALL cp_fm_release(fm_A_BSE)
         IF (bse_env%do_abba) THEN
            CALL diagonalize_C(fm_C_BSE, bse_env%homo_red, bse_env%virt_red, bse_env%homo_full, &
                               fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                               unit_nr, diag_runtime_est, bse_env, mo_coeff)
         END IF
         CALL cp_fm_release(fm_C_BSE)
      END IF

      DEALLOCATE (eigenval_joint)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A53)') 'BSE|', 'The BSE was successfully calculated. Have a nice day!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE bse_solve

! **************************************************************************************************
!> \brief The AO multipoles M^k_µν = <φ_µ|(r - r_0)^k|φ_ν> of every order up to n_moments, on the
!>        block pattern of the overlap, about r_0, the centre of atomic charges;
!>        k = x, y, z for order 1, then xx, xy, xz, yy, yz, zz for order 2
!> \param qs_env source of the basis, the geometry and the overlap's block pattern
!> \param n_moments highest order, 1 for dipoles, 2 adds quadrupoles
!> \param matrix_multipole M^k_µν, one entry per component k, allocated here; the caller releases it
!> \param rpoint r_0
! **************************************************************************************************
   SUBROUTINE get_multipoles_ao(qs_env, n_moments, matrix_multipole, rpoint)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: n_moments
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_multipole
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: rpoint

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_multipoles_ao'

      INTEGER                                            :: handle, idir, n_multipole
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s, sab_orb, ref_point)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)

      ! components of every order up to n_moments: 3 dipoles, then 6 quadrupoles
      n_multipole = (6 + 11*n_moments + 6*n_moments**2 + n_moments**3)/6 - 1
      NULLIFY (matrix_multipole)
      CALL dbcsr_allocate_matrix_set(matrix_multipole, n_multipole)
      DO idir = 1, n_multipole
         CALL dbcsr_init_p(matrix_multipole(idir)%matrix)
         CALL dbcsr_create(matrix_multipole(idir)%matrix, name="ao_multipole", &
                           template=matrix_s(1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_multipole(idir)%matrix, sab_orb)
         CALL dbcsr_set(matrix_multipole(idir)%matrix, 0._dp)
      END DO

      CALL get_reference_point(rpoint, qs_env=qs_env, reference=use_mom_ref_coac, ref_point=ref_point)
      CALL build_local_moment_matrix(qs_env, matrix_multipole, n_moments, ref_point=rpoint)

      CALL timestop(handle)

   END SUBROUTINE get_multipoles_ao

! **************************************************************************************************
!> \brief Early memory check for the full diagonalization, before any GW work: the window on the
!>        DFT axis gives n_ov, the estimate is estimate_BSE_resources', the budget the one the
!>        iterative solver uses. WARN and CLAMP warn, ABORT aborts, OFF prints the estimate only.
!> \param bse_env settings; the window cutoffs, MEMORY_CHECK and MEMORY_BUDGET_GB
!> \param mos the MO sets the GW will correct, for their DFT eigenvalues and occupations
!> \param para_env the ranks that will hold the matrices
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE bse_fulldiag_memory_check(bse_env, mos, para_env, unit_nr)

      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_fulldiag_memory_check'

      CHARACTER(LEN=16)                                  :: bud_str, mem_str
      CHARACTER(LEN=64)                                  :: remedy
      INTEGER                                            :: handle, ispin, n_ov_joint, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo, n_mo
      LOGICAL                                            :: from_free, over, skipped
      REAL(KIND=dp)                                      :: budget_GB, mem_avail_GB, mem_est_GB
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenval_dft
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(mo_window_type)                               :: window

      CALL timeset(routineN, handle)

      nspins = SIZE(mos)
      ALLOCATE (homo(nspins), n_mo(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), nmo=n_mo(ispin), homo=homo(ispin), eigenvalues=eigenvalues)
         IF (ispin == 1) ALLOCATE (eigenval_dft(n_mo(1), nspins))
         eigenval_dft(:, ispin) = eigenvalues(1:n_mo(ispin))
      END DO
      ! the preparer warns about a combined window; the estimate stays silent
      CALL bse_window_from_dft(eigenval_dft, n_mo, homo, bse_env%bse_cutoff_occ, &
                               bse_env%bse_cutoff_empty, window, warn=.FALSE.)
      n_ov_joint = 0
      DO ispin = 1, nspins
         n_ov_joint = n_ov_joint + (homo(ispin) - window%first_mo + 1)*(window%last_mo - homo(ispin))
      END DO
      DEALLOCATE (homo, n_mo, eigenval_dft)

      mem_est_GB = fulldiag_memory_estimate_GB(n_ov_joint, bse_env%flag_tda /= bse_tda)/ &
                   REAL(para_env%num_pe, KIND=dp)

      ! the budget: as given, or the share mem_fraction of the free memory per rank
      from_free = bse_env%memory_budget_gb < 0.0_dp
      skipped = .FALSE.
      budget_GB = bse_env%memory_budget_gb
      IF (from_free) THEN
         CALL mp_mem_avail_per_rank_GB(para_env, mem_avail_GB)
         budget_GB = mem_fraction*mem_avail_GB
         skipped = mem_avail_GB <= 0.0_dp
      END IF
      over = .NOT. skipped .AND. bse_env%memory_check /= bse_memcheck_off .AND. mem_est_GB > budget_GB

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Early memory check for the full diagonalization'
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Transition pairs in the active window (DFT axis)', &
            n_ov_joint
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Peak memory estimate per MPI rank (GB)', mem_est_GB
         IF (.NOT. skipped) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Memory budget per MPI rank (GB)', budget_GB
         END IF
         IF (bse_env%memory_check == bse_memcheck_off) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: OFF'
         ELSE IF (skipped) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: skipped, free memory not detectable'
         ELSE IF (over) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: the estimate exceeds the budget'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: within the budget'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      IF (over) THEN
         WRITE (mem_str, '(F16.3)') mem_est_GB
         WRITE (bud_str, '(F16.3)') budget_GB
         ! ITERDIAG is closed-shell only, so an open-shell run is not pointed to it
         IF (nspins == 1) THEN
            remedy = "more MPI ranks, a tighter window, ITERDIAG"
         ELSE
            remedy = "more MPI ranks, a tighter window"
         END IF
         SELECT CASE (bse_env%memory_check)
         CASE (bse_memcheck_warn, bse_memcheck_clamp)
            IF (unit_nr > 0) CALL cp_warn(__LOCATION__, &
                                          "BSE full diagonalization: the estimate of "//TRIM(ADJUSTL(mem_str))// &
                                          " GB per MPI rank exceeds the budget of "//TRIM(ADJUSTL(bud_str))// &
                                          " GB. Try one of: "//TRIM(remedy)//".")
         CASE (bse_memcheck_abort)
            CALL cp_abort(__LOCATION__, &
                          "BSE full diagonalization: the estimate of "//TRIM(ADJUSTL(mem_str))// &
                          " GB per MPI rank exceeds the budget of "//TRIM(ADJUSTL(bud_str))// &
                          " GB. Try one of: "//TRIM(remedy)//", MEMORY_CHECK WARN.")
         END SELECT
      END IF

      CALL timestop(handle)

   END SUBROUTINE bse_fulldiag_memory_check

END MODULE bse_main
