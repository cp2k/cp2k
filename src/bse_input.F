!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Reads the &BSE input section into the BSE environment; called by read_mp2_section with
!>        the section under RI_RPA%GW, and by any other GW path with its own section
!> \par History
!>      09.2026 moved out of read_mp2_section [Maximilian Graml]
! **************************************************************************************************
MODULE bse_input
   USE bse_types,                       ONLY: bse_env_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_input'

   PUBLIC :: read_bse_section

CONTAINS

! **************************************************************************************************
!> \brief Reads every keyword of the &BSE section into the settings of the BSE environment
!> \param bse_section the &BSE section of the calling GW path
!> \param bse_env receives every keyword and keeps the section pointer for the NTO print keys
! **************************************************************************************************
   SUBROUTINE read_bse_section(bse_section, bse_env)
      TYPE(section_vals_type), POINTER                   :: bse_section
      TYPE(bse_env_type), INTENT(INOUT)                  :: bse_env

      INTEGER                                            :: print_level_NTO_analysis
      LOGICAL                                            :: explicit

      bse_env%bse_section => bse_section
      CALL section_vals_val_get(bse_section, "_SECTION_PARAMETERS_", &
                                l_val=bse_env%do_bse)
      CALL section_vals_val_get(bse_section, "SPIN_CONFIG", &
                                i_val=bse_env%bse_spin_config)
      CALL section_vals_val_get(bse_section, "ENERGY_CUTOFF_OCC", &
                                r_val=bse_env%bse_cutoff_occ)
      CALL section_vals_val_get(bse_section, "ENERGY_CUTOFF_EMPTY", &
                                r_val=bse_env%bse_cutoff_empty)
      CALL section_vals_val_get(bse_section, "NUM_PRINT_EXC", &
                                i_val=bse_env%num_print_exc)
      CALL section_vals_val_get(bse_section, "NUM_PRINT_EXC_DESCR", &
                                i_val=bse_env%num_print_exc_descr)
      CALL section_vals_val_get(bse_section, "PRINT_DIRECTIONAL_EXC_DESCR", &
                                l_val=bse_env%print_directional_exc_descr)
      CALL section_vals_val_get(bse_section, "PRINT_DIRECTIONAL_CROSSCORRELATION", &
                                l_val=bse_env%print_directional_crosscorrelation, explicit=explicit)
      IF (explicit .AND. bse_env%print_directional_crosscorrelation .AND. &
          .NOT. bse_env%print_directional_exc_descr) THEN
         CALL cp_warn(__LOCATION__, &
                      "PRINT_DIRECTIONAL_CROSSCORRELATION has no effect without PRINT_DIRECTIONAL_EXC_DESCR.")
      END IF
      CALL section_vals_val_get(bse_section, "EPS_X", &
                                r_val=bse_env%eps_x)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%NUM_EXC_EN", &
                                i_val=bse_env%num_exc_en)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%NUM_BUFFER_STATES", &
                                i_val=bse_env%num_buffer_states)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%BLOCK_SIZE", &
                                i_val=bse_env%block_size)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%NUM_GUESS_TRANSITIONS", &
                                i_val=bse_env%num_guess_transitions)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%MAX_SUBSPACE_FACTOR", &
                                i_val=bse_env%max_subspace_factor)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%MATVEC_BLOCK_COLUMNS", &
                                i_val=bse_env%matvec_block_columns)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%MEMORY_CHECK", &
                                i_val=bse_env%memory_check)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%MEMORY_BUDGET_GB", &
                                r_val=bse_env%memory_budget_gb)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%MAX_ITER", &
                                i_val=bse_env%max_iter)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%CONVERGENCE_CRITERION", &
                                i_val=bse_env%convergence_criterion)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%EPS_RES", &
                                r_val=bse_env%eps_res)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%EPS_ENERGY", &
                                r_val=bse_env%eps_energy)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%PRECONDITIONER", &
                                i_val=bse_env%preconditioner)
      CALL section_vals_val_get(bse_section, "BSE_ITERAT%ABBA_SOLVER", &
                                i_val=bse_env%abba_solver)
      CALL section_vals_val_get(bse_section, "BSE_DIAG_METHOD", &
                                i_val=bse_env%bse_diag_method)
      CALL section_vals_val_get(bse_section, "TDA", &
                                i_val=bse_env%flag_tda)
      CALL section_vals_val_get(bse_section, "USE_KS_ENERGIES", &
                                l_val=bse_env%use_ks_energies)
      CALL section_vals_val_get(bse_section, "SCREENING_IN_W%_SECTION_PARAMETERS_", &
                                i_val=bse_env%screening_method)
      CALL section_vals_val_get(bse_section, "SCREENING_IN_W%ALPHA", &
                                r_val=bse_env%screening_factor)
      CALL section_vals_val_get(bse_section, "BSE_DEBUG_PRINT", &
                                l_val=bse_env%bse_debug_print)
      CALL section_vals_val_get(bse_section, "BSE_SPECTRUM%_SECTION_PARAMETERS_", &
                                l_val=bse_env%bse_print_spectrum)
      CALL section_vals_val_get(bse_section, "BSE_SPECTRUM%FREQUENCY_STEP_SIZE", &
                                r_val=bse_env%bse_spectrum_freq_step_size)
      CALL section_vals_val_get(bse_section, "BSE_SPECTRUM%FREQUENCY_STARTING_POINT", &
                                r_val=bse_env%bse_spectrum_freq_start)
      CALL section_vals_val_get(bse_section, "BSE_SPECTRUM%FREQUENCY_END_POINT", &
                                r_val=bse_env%bse_spectrum_freq_end)
      CALL section_vals_val_get(bse_section, "BSE_SPECTRUM%ETA_LIST", &
                                r_vals=bse_env%bse_eta_spectrum_list)
      CALL section_vals_val_get(bse_section, "NTO_ANALYSIS%_SECTION_PARAMETERS_", &
                                i_val=print_level_NTO_analysis)
      ! All other print levels besides "off" should trigger the analysis
      ! cf. input/cp_output_handling.F:cp_print_key_section_create
      IF (print_level_NTO_analysis /= debug_print_level + 1) THEN
         bse_env%do_nto_analysis = .TRUE.
      ELSE
         bse_env%do_nto_analysis = .FALSE.
      END IF
      CALL section_vals_val_get(bse_section, "NTO_ANALYSIS%EPS_NTO_EIGVAL", &
                                r_val=bse_env%eps_nto_eigval)
      CALL section_vals_val_get(bse_section, "NTO_ANALYSIS%EPS_OSC_STR", &
                                r_val=bse_env%eps_nto_osc_str)
      CALL section_vals_val_get(bse_section, "NTO_ANALYSIS%NUM_PRINT_EXC_NTOS", &
                                i_val=bse_env%num_print_exc_ntos)
      CALL section_vals_val_get(bse_section, "NTO_ANALYSIS%STATE_LIST", &
                                EXPLICIT=bse_env%explicit_nto_list)
      IF (bse_env%explicit_nto_list) THEN
         CALL section_vals_val_get(bse_section, "NTO_ANALYSIS%STATE_LIST", &
                                   i_vals=bse_env%bse_nto_state_list)
      END IF

   END SUBROUTINE read_bse_section

END MODULE bse_input
