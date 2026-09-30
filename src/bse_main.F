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
   USE bse_input,                       ONLY: read_bse_section
   USE bse_iterative,                   ONLY: solve_bse_iteratively
   USE bse_matvec,                      ONLY: mem_fraction
   USE bse_print,                       ONLY: print_BSE_start_flag
   USE bse_types,                       ONLY: bse_env_type
   USE bse_util,                        ONLY: adapt_BSE_input_params,&
                                              bse_window_from_dft,&
                                              bse_window_from_mos,&
                                              estimate_BSE_resources,&
                                              fulldiag_memory_estimate_GB,&
                                              get_bse_spin_block_layout,&
                                              screen_slabs,&
                                              truncate_BSE_matrices
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_add_on_diag,&
                                              cp_fm_scale,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE input_constants,                 ONLY: &
        bse_abba, bse_both, bse_fulldiag, bse_iterdiag, bse_memcheck_abort, bse_memcheck_clamp, &
        bse_memcheck_off, bse_memcheck_warn, bse_screening_w0, bse_tda
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE mo_window,                       ONLY: mo_window_type
   USE pw_poisson_read_input,           ONLY: pw_poisson_read_parameters
   USE pw_poisson_types,                ONLY: pw_poisson_parameter_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: get_multipoles_ao
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_main'

   PUBLIC :: prepare_bse_env_from_ri_rpa_gw, prepare_bse_env_from_bandstructure_gw, bse_solve, &
             bse_fulldiag_memory_check, bse_fulldiag_memory_check_from_mos

CONTAINS

! **************************************************************************************************
!> \brief Entry for the BSE of RI_RPA%GW, which hands over the slabs of all its GW levels and the
!>        static polarizability Q(0): re-reads the settings, cuts the window on the DFT axis out of
!>        the slabs and both energy sets, and inverts 1+Q(0). Starts no solver.
!> \param bse_env settings re-read from its section; state filled here
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
   SUBROUTINE prepare_bse_env_from_ri_rpa_gw(bse_env, fm_S_ia_full, fm_S_ij_full, fm_S_ab_full, fm_Q, &
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_bse_env_from_ri_rpa_gw'

      INTEGER                                            :: handle, ispin, n_window, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_mo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_dft_spin, eigenval_qp_spin
      TYPE(mo_window_type)                               :: window
      TYPE(section_vals_type), POINTER                   :: bse_section

      CALL timeset(routineN, handle)

      ! the settings afresh from the input, since adapt_BSE_input_params normalises them in place and a
      ! later energy evaluation on the same environment runs the BSE again; the local pointer keeps the
      ! section apart from bse_env, whose section pointer the reader sets
      bse_section => bse_env%bse_section
      CALL read_bse_section(bse_section, bse_env)

      nspins = SIZE(homo)
      ! the GW-corrected level counts are the ceiling of the window
      CALL start_bse_env(bse_env, nspins, dimen_RI, dimen_RI_red, gw_corr_lev_occ, bse_lev_virt, unit_nr)

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

      CALL invert_static_screening(fm_Q, dimen_RI_red, bse_env%fm_eps_inv)

      IF (bse_env%wc_scaling /= 1.0_dp) THEN
         IF (bse_env%screening_method /= bse_screening_w0 .OR. bse_env%screening_factor > 0.0_dp) THEN
            CPABORT("SCREENING_IN_W other than W_0 cannot be combined with a scaled W^c (static COHSEX).")
         END IF
         CALL scale_static_screening(bse_env%fm_eps_inv, dimen_RI_red, bse_env%wc_scaling)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A10,T13,A,T71,F10.4)') 'BSE|', &
               'Direct term W = v + alpha W^c(0), alpha =', bse_env%wc_scaling
         END IF
      END IF

      CALL finish_bse_env(bse_env, qs_env, fm_S_ia_full(1)%matrix_struct%para_env)

      CALL timestop(handle)

   END SUBROUTINE prepare_bse_env_from_ri_rpa_gw

! **************************************************************************************************
!> \brief Entry for the BSE of PROPERTIES%BANDSTRUCTURE%GW (PROPERTIES%BSE), which hands over
!>        window-sized slabs and a ready inverse dielectric matrix: nothing is cut or inverted, the
!>        slabs are taken over as they are. Shares start_bse_env and finish_bse_env with
!>        prepare_bse_env_from_ri_rpa_gw. Starts no solver.
!> \param bse_env settings read from its section by the caller; state filled here
!> \param fm_S_ia B^P_ia of the window per spin; moved into bse_env, deallocated on return
!> \param fm_S_ij B^P_ij of the window per spin; moved into bse_env, deallocated on return
!> \param fm_S_ab B^P_ab of the window per spin; moved into bse_env, deallocated on return
!> \param fm_eps_inv the inverse dielectric matrix ε^-1_PQ the core screens the slabs with,
!>        \bar{B}^P_ij = sum_Q ε^-1_PQ B^Q_ij; RI x RI, copied; stays the caller's
!> \param eigenval_qp quasiparticle energies of the window (level, spin)
!> \param eigenval_dft DFT energies of the window (level, spin)
!> \param homo_full occupied MOs per spin of the full reference
!> \param homo_red occupied levels in the window per spin
!> \param virt_red virtual levels in the window per spin
!> \param dimen_RI RI dimension of the slabs and of fm_eps_inv
!> \param qs_env source of the pointers the solver reads and of the AO multipoles
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE prepare_bse_env_from_bandstructure_gw(bse_env, fm_S_ia, fm_S_ij, fm_S_ab, fm_eps_inv, &
                                                    eigenval_qp, eigenval_dft, homo_full, homo_red, virt_red, &
                                                    dimen_RI, qs_env, unit_nr)

      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: fm_S_ia, fm_S_ij, fm_S_ab
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eps_inv
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eigenval_qp, eigenval_dft
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo_full, homo_red, virt_red
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_bse_env_from_bandstructure_gw'

      INTEGER                                            :: handle, ispin, n_window, nspins

      CALL timeset(routineN, handle)

      ! the caller (bse_from_bs_env) read the settings afresh from the input; the entry of RI_RPA%GW reads them itself
      IF (.NOT. ASSOCIATED(bse_env%bse_section)) CPABORT("BSE settings not read: read_bse_section comes first")

      nspins = SIZE(homo_full)
      IF (SIZE(fm_S_ia) /= nspins .OR. SIZE(eigenval_qp, 2) /= nspins) THEN
         CPABORT("BSE needs one set of RI integrals and one set of energies per spin")
      END IF
      ! every window level carries a GW energy here, and the inverse dielectric matrix has the full RI dimension
      CALL start_bse_env(bse_env, nspins, dimen_RI, dimen_RI, MAXVAL(homo_red), virt_red, unit_nr)

      ALLOCATE (bse_env%homo_red(nspins), bse_env%virt_red(nspins), bse_env%homo_full(nspins))
      bse_env%homo_full(:) = homo_full(:)
      bse_env%homo_red(:) = homo_red(:)
      bse_env%virt_red(:) = virt_red(:)
      n_window = SIZE(eigenval_qp, 1)
      DO ispin = 1, nspins
         ! one absolute window for every spin: its occupied and virtual levels fill it
         IF (homo_red(ispin) + virt_red(ispin) /= n_window) THEN
            CPABORT("Occupied and virtual levels do not fill the BSE window")
         END IF
      END DO
      ALLOCATE (bse_env%eigenval_dft(n_window, nspins), bse_env%eigenval_qp(n_window, nspins))
      bse_env%eigenval_dft(:, :) = eigenval_dft(:, :)
      bse_env%eigenval_qp(:, :) = eigenval_qp(:, :)

      ! the slabs are window-sized already: take them over
      CALL MOVE_ALLOC(fm_S_ia, bse_env%fm_S_ia)
      CALL MOVE_ALLOC(fm_S_ij, bse_env%fm_S_ij)
      CALL MOVE_ALLOC(fm_S_ab, bse_env%fm_S_ab)
      CALL cp_fm_create(bse_env%fm_eps_inv, fm_eps_inv%matrix_struct)
      CALL cp_fm_to_fm(fm_eps_inv, bse_env%fm_eps_inv)

      CALL finish_bse_env(bse_env, qs_env, bse_env%fm_S_ia(1)%matrix_struct%para_env)

      CALL timestop(handle)

   END SUBROUTINE prepare_bse_env_from_bandstructure_gw

! **************************************************************************************************
!> \brief The part both entries share before their window: the sizes, the run flags from the
!>        settings, the start banner and the debug flag
!> \param bse_env settings read from its section; sizes and flags set here
!> \param nspins spin channels
!> \param dimen_RI RI dimension of the slabs
!> \param dimen_RI_red RI dimension of the inverse dielectric matrix
!> \param gw_corr_lev_occ GW-corrected occupied levels
!> \param gw_corr_lev_virt GW-corrected virtual levels per spin
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE start_bse_env(bse_env, nspins, dimen_RI, dimen_RI_red, gw_corr_lev_occ, gw_corr_lev_virt, unit_nr)

      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env
      INTEGER, INTENT(IN)                                :: nspins, dimen_RI, dimen_RI_red, &
                                                            gw_corr_lev_occ
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_virt
      INTEGER, INTENT(IN)                                :: unit_nr

      TYPE(cp_logger_type), POINTER                      :: logger

      bse_env%nspins = nspins
      bse_env%dimen_RI = dimen_RI
      bse_env%dimen_RI_red = dimen_RI_red
      bse_env%gw_corr_lev_occ = gw_corr_lev_occ
      ALLOCATE (bse_env%gw_corr_lev_virt(nspins))
      bse_env%gw_corr_lev_virt(:) = gw_corr_lev_virt(:)

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

   END SUBROUTINE start_bse_env

! **************************************************************************************************
!> \brief The part both entries share once the window and the slabs are in place: the
!>        QS pointers, the open-shell checks, the spin-block layout and the AO multipoles
!> \param bse_env state filled here
!> \param qs_env source of the pointers the solver reads and of the AO multipoles
!> \param para_env_slabs the communicator the slabs live on
! **************************************************************************************************
   SUBROUTINE finish_bse_env(bse_env, qs_env, para_env_slabs)

      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env_slabs

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'finish_bse_env'

      INTEGER                                            :: handle, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_ov, offsets
      TYPE(pw_poisson_parameter_type)                    :: poisson_params
      TYPE(section_vals_type), POINTER                   :: input, poisson_section

      CALL timeset(routineN, handle)

      nspins = bse_env%nspins

      ! what the solver reads from QS and the GW (pointers, none owned)
      bse_env%para_env => para_env_slabs
      CALL get_qs_env(qs_env, subsys=bse_env%subsys, pw_env=bse_env%pw_env, &
                      dft_control=bse_env%dft_control, mos=bse_env%mos, &
                      exstate_env=bse_env%exstate_env)
      ! the Poisson periodicity from the input: PROPERTIES%BSE runs after the GW released the grids
      NULLIFY (input, poisson_section)
      CALL get_qs_env(qs_env, input=input)
      poisson_section => section_vals_get_subs_vals(input, "DFT%POISSON")
      CALL pw_poisson_read_parameters(poisson_section, poisson_params)
      bse_env%n_periodic_poisson = COUNT(poisson_params%periodic == 1)
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

      ALLOCATE (n_ov(nspins), offsets(nspins))
      CALL get_bse_spin_block_layout(bse_env%homo_red, bse_env%virt_red, n_ov, offsets, bse_env%n_ov_joint)
      DEALLOCATE (n_ov, offsets)

      ! the AO multipoles for the transition dipoles and, when asked for, the exciton descriptors;
      ! the open-shell path computes no descriptors
      CALL get_multipoles_ao(qs_env, 1, bse_env%matrix_dipole_ao, bse_env%multipole_ref_point)
      IF (bse_env%num_print_exc_descr /= 0 .AND. nspins == 1) THEN
         CALL get_multipoles_ao(qs_env, 2, bse_env%matrix_quadpole_ao, bse_env%multipole_ref_point)
      END IF

      CALL timestop(handle)

   END SUBROUTINE finish_bse_env

! **************************************************************************************************
!> \brief Scale the correlation part of the static screened Coulomb interaction in the RI kernel
!>
!>        The direct terms of the BSE use (pq|W|rs) = Σ_PQ B^P_pq K_PQ B^Q_rs with the RI kernel
!>        K = ε^(-1)(0). With W = v + α W^c(0), W^c(0) = (ε^(-1)(0) - 1) v, the kernel becomes
!>          K_PQ = δ_PQ + α (ε^(-1)_PQ(0) - δ_PQ) = α ε^(-1)_PQ(0) + (1 - α) δ_PQ
!>        such that
!>          A_ia,jb = (ϵ_a - ϵ_i) δ_ij δ_ab + 2 (ia|jb) - (ij|ab) - α (ij|W^c|ab)
!>          B_ia,jb = 2 (ia|jb) - (ib|ja) - α (ib|W^c|ja)
!>        α = 1: standard BSE@W0, α = 0: unscreened (TDHF-like) kernel.
!> \param fm_eps_inv ...
!> \param dimen_RI_red ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE scale_static_screening(fm_eps_inv, dimen_RI_red, alpha)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eps_inv
      INTEGER, INTENT(IN)                                :: dimen_RI_red
      REAL(KIND=dp), INTENT(IN)                          :: alpha

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_static_screening'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! K_PQ = α ε^(-1)_PQ + (1 - α) δ_PQ on the leading dimen_RI_red block
      CALL cp_fm_scale(alpha, fm_eps_inv)
      CALL cp_fm_add_on_diag(fm_eps_inv, 1.0_dp - alpha, dimen_RI_red)

      CALL timestop(handle)

   END SUBROUTINE scale_static_screening

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

      INTEGER                                            :: handle, info_chol
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      ! a copy: the caller keeps Q(0), and the open-shell path inverts once for both spins
      CALL cp_fm_create(fm_eps_inv, fm_Q%matrix_struct)
      CALL cp_fm_to_fm(fm_Q, fm_eps_inv)
      CALL cp_fm_create(fm_work, fm_Q%matrix_struct)
      CALL cp_fm_set_all(fm_work, 0.0_dp)

      ! δ_PQ + Q_PQ(0) on the leading dimen_RI_red block
      CALL cp_fm_add_on_diag(fm_eps_inv, 1.0_dp, dimen_RI_red)

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
!> \brief Early memory check for the full diagonalization, before any GW work: the window on the
!>        DFT axis gives n_ov = sum_spin n_occ n_virt, the estimate is estimate_BSE_resources', the
!>        budget the one the iterative solver uses. WARN and CLAMP warn, ABORT aborts, OFF prints
!>        the estimate only.
!> \param bse_env settings; the TDA flag, MEMORY_CHECK and MEMORY_BUDGET_GB
!> \param window the active MO window on the DFT axis, one absolute range for every spin
!> \param homo occupied MOs per spin
!> \param para_env the ranks that will hold the matrices
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE bse_fulldiag_memory_check(bse_env, window, homo, para_env, unit_nr)

      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_fulldiag_memory_check'

      CHARACTER(LEN=16)                                  :: bud_str, mem_str
      CHARACTER(LEN=64)                                  :: remedy
      INTEGER                                            :: handle, ispin, n_ov_joint, nspins
      LOGICAL                                            :: from_free, over, skipped
      REAL(KIND=dp)                                      :: budget_GB, mem_avail_GB, mem_est_GB

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      n_ov_joint = 0
      DO ispin = 1, nspins
         n_ov_joint = n_ov_joint + (homo(ispin) - window%first_mo + 1)*(window%last_mo - homo(ispin))
      END DO

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

! **************************************************************************************************
!> \brief bse_fulldiag_memory_check for a caller that holds MO sets: the window on the DFT axis
!>        from their eigenvalues, the occupied MOs per spin from their homo
!> \param bse_env settings; the window cutoffs, the TDA flag, MEMORY_CHECK and MEMORY_BUDGET_GB
!> \param mos the MO sets the GW will correct, for their DFT eigenvalues and occupations
!> \param para_env the ranks that will hold the matrices
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE bse_fulldiag_memory_check_from_mos(bse_env, mos, para_env, unit_nr)

      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo
      TYPE(mo_window_type)                               :: window

      ! the preparer warns about a combined window; the estimate stays silent
      CALL bse_window_from_mos(mos, bse_env%bse_cutoff_occ, bse_env%bse_cutoff_empty, window, warn=.FALSE.)
      ALLOCATE (homo(SIZE(mos)))
      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mos(ispin), homo=homo(ispin))
      END DO
      CALL bse_fulldiag_memory_check(bse_env, window, homo, para_env, unit_nr)
      DEALLOCATE (homo)

   END SUBROUTINE bse_fulldiag_memory_check_from_mos

END MODULE bse_main
