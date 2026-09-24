!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The BSE environment: the settings of the &BSE section and the state a GW path prepares
!>        for the solver core (window, energies, window-sized RI slabs, inverse dielectric matrix)
!> \par History
!>      09.2026 created from the bse_type of mp2_types [Maximilian Graml]
! **************************************************************************************************
MODULE bse_types
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE exstates_types,                  ONLY: excited_energy_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_subsys_types,                 ONLY: qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_types'

   PUBLIC :: bse_env_type, bse_env_release

! **************************************************************************************************
!> \brief Settings of the &BSE section (read_bse_section, normalised in place by
!>        adapt_BSE_input_params) and the state prepare_bse_env fills for bse_solve
! **************************************************************************************************
   TYPE bse_env_type
      ! --- settings, one field per keyword of the &BSE section; the reader sets them ---
      INTEGER                                            :: bse_spin_config = 0, &
                                                            bse_diag_method = 0, &
                                                            flag_tda = 0, &
                                                            num_exc_en = 0, &
                                                            num_print_exc = 0, &
                                                            num_print_exc_descr = 0, &
                                                            screening_method = 0, &
                                                            num_buffer_states = 0, &
                                                            block_size = 0, &
                                                            num_guess_transitions = 0, &
                                                            max_subspace_factor = 0, &
                                                            matvec_block_columns = -1, &
                                                            memory_check = 0, &
                                                            max_iter = 0, &
                                                            convergence_criterion = 0, &
                                                            preconditioner = 0, &
                                                            abba_solver = 0
      REAL(KIND=dp)                                      :: eps_res = 0.0_dp, &
                                                            memory_budget_gb = 0.0_dp, &
                                                            eps_energy = 0.0_dp, &
                                                            eps_x = 0.0_dp, &
                                                            screening_factor = 0.0_dp, &
                                                            bse_cutoff_occ = 0.0_dp, &
                                                            bse_cutoff_empty = 0.0_dp
      LOGICAL                                            :: do_bse = .FALSE., &
                                                            bse_debug_print = .FALSE., &
                                                            print_directional_exc_descr = .FALSE., &
                                                            print_directional_crosscorrelation = .FALSE., &
                                                            use_ks_energies = .FALSE.
      !BSE optical spectrum
      REAL(KIND=dp)                                      :: bse_spectrum_freq_step_size = 0.0_dp, &
                                                            bse_spectrum_freq_start = 0.0_dp, &
                                                            bse_spectrum_freq_end = 0.0_dp
      LOGICAL                                            :: bse_print_spectrum = .FALSE.
      REAL(KIND=dp), DIMENSION(:), POINTER               :: bse_eta_spectrum_list => NULL()

      !BSE NTO prints
      LOGICAL                                            :: do_nto_analysis = .FALSE., &
                                                            explicit_nto_list = .FALSE.
      REAL(KIND=dp)                                      :: eps_nto_eigval = 0.0_dp, &
                                                            eps_nto_osc_str = 0.0_dp
      INTEGER                                            :: num_print_exc_ntos = 0
      INTEGER, DIMENSION(:), POINTER                     :: bse_nto_state_list => NULL()
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: bse_nto_state_list_final
      TYPE(section_vals_type), POINTER                   :: bse_section => NULL()
      ! --- state, unset until prepare_bse_env fills it ---
      LOGICAL                                            :: do_tda = .FALSE., &
                                                            do_abba = .FALSE., &
                                                            do_fulldiag = .FALSE., &
                                                            do_iterdiag = .FALSE.
      INTEGER                                            :: nspins = 0, &
                                                            dimen_RI = 0, &
                                                            dimen_RI_red = 0, &
                                                            n_ov_joint = 0, &
                                                            gw_corr_lev_occ = 0
      ! per spin: the window's occupied and virtual counts, the full occupied count (MO labels and
      ! N_e), the OV-pair count, the block offset in the joint matrix, the GW-corrected virtual count
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: homo_red, virt_red, homo_full, n_ov, &
                                                            offsets, gw_corr_lev_virt
      ! windowed energies (level, spin), DFT and quasiparticle; bse_solve picks one by USE_KS_ENERGIES
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: eigenval_dft, eigenval_qp
      ! window-sized RI slabs per spin, [1+Q(0)]^-1, and the screened slabs bse_solve forms
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: fm_S_ia, fm_S_ij, fm_S_ab, &
                                                            fm_S_bar_ia, fm_S_bar_ij
      TYPE(cp_fm_type)                                   :: fm_eps_inv = cp_fm_type()
      ! AO multipoles <φ_µ|(r - r_0)^k|φ_ν> from QS, dipoles always and quadrupoles for the exciton
      ! descriptors, about r_0; the solver transforms them into the window
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_dipole_ao => NULL(), &
                                                            matrix_quadpole_ao => NULL()
      REAL(KIND=dp), DIMENSION(3)                        :: multipole_ref_point = 0.0_dp
      ! QS objects the solver reads, none owned
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env => NULL()
      TYPE(qs_subsys_type), POINTER                      :: subsys => NULL()
      TYPE(pw_env_type), POINTER                         :: pw_env => NULL()
      TYPE(dft_control_type), POINTER                    :: dft_control => NULL()
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s => NULL()
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb => NULL()
      TYPE(excited_energy_type), POINTER                 :: exstate_env => NULL()
   END TYPE bse_env_type

CONTAINS

! **************************************************************************************************
!> \brief Releases the state a GW path prepared and nullifies the QS pointers; the settings stay,
!>        including the section pointer and the two lists into the input tree, since a later energy
!>        evaluation on the same environment runs the BSE again. Safe on an unprepared environment.
!> \param bse_env the environment whose prepared state is released
! **************************************************************************************************
   SUBROUTINE bse_env_release(bse_env)
      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env

      CALL cp_fm_release(bse_env%fm_S_ia)
      CALL cp_fm_release(bse_env%fm_S_ij)
      CALL cp_fm_release(bse_env%fm_S_ab)
      CALL cp_fm_release(bse_env%fm_S_bar_ia)
      CALL cp_fm_release(bse_env%fm_S_bar_ij)
      CALL cp_fm_release(bse_env%fm_eps_inv)
      CALL dbcsr_deallocate_matrix_set(bse_env%matrix_dipole_ao)
      CALL dbcsr_deallocate_matrix_set(bse_env%matrix_quadpole_ao)
      IF (ALLOCATED(bse_env%homo_red)) DEALLOCATE (bse_env%homo_red)
      IF (ALLOCATED(bse_env%virt_red)) DEALLOCATE (bse_env%virt_red)
      IF (ALLOCATED(bse_env%homo_full)) DEALLOCATE (bse_env%homo_full)
      IF (ALLOCATED(bse_env%n_ov)) DEALLOCATE (bse_env%n_ov)
      IF (ALLOCATED(bse_env%offsets)) DEALLOCATE (bse_env%offsets)
      IF (ALLOCATED(bse_env%gw_corr_lev_virt)) DEALLOCATE (bse_env%gw_corr_lev_virt)
      IF (ALLOCATED(bse_env%eigenval_dft)) DEALLOCATE (bse_env%eigenval_dft)
      IF (ALLOCATED(bse_env%eigenval_qp)) DEALLOCATE (bse_env%eigenval_qp)
      IF (ALLOCATED(bse_env%bse_nto_state_list_final)) DEALLOCATE (bse_env%bse_nto_state_list_final)
      NULLIFY (bse_env%para_env, bse_env%blacs_env, bse_env%subsys, bse_env%pw_env, &
               bse_env%dft_control, bse_env%mos, bse_env%matrix_s, bse_env%sab_orb, &
               bse_env%exstate_env)

   END SUBROUTINE bse_env_release

END MODULE bse_types
