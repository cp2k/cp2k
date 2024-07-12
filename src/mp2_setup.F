!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! **************************************************************************************************
MODULE mp2_setup
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_init_read_input
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_parser_methods,               ONLY: read_float_object
   USE input_constants,                 ONLY: &
        do_eri_mme, do_potential_short, evGW, evGW0, mp2_method_direct, mp2_method_gpw, &
        mp2_method_none, mp2_ri_optimize_basis, ri_mp2_laplace, ri_mp2_method_gpw, &
        ri_rpa_method_gpw
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              max_line_length
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: erfc_cutoff
   USE mp2_types,                       ONLY: mp2_method_direct,&
                                              mp2_method_gpw,&
                                              mp2_method_none,&
                                              mp2_ri_optimize_basis,&
                                              mp2_type,&
                                              ri_mp2_laplace,&
                                              ri_mp2_method_gpw,&
                                              ri_rpa_method_gpw
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_setup'

   PUBLIC :: read_mp2_section

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param input ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE read_mp2_section(input, mp2_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_mp2_section'

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: string_pointer
      CHARACTER(LEN=max_line_length)                     :: error_message
      INTEGER                                            :: gw_sc, handle, i, i_special_kp, ival, &
                                                            unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: do_mp2, do_opt_ri_basis, do_ri_mp2, &
                                                            do_ri_sos_mp2, do_rpa
      REAL(KIND=dp), DIMENSION(:), POINTER               :: r_vals
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: bse_section, cphf_section, &
                                                            eri_mme_section, gw_section, &
                                                            low_scaling_section, mp2_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      mp2_section => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION")

      mp2_env%method = mp2_method_none

      ! should come from input
      CALL section_vals_val_get(mp2_section, "MEMORY", r_val=mp2_env%mp2_memory)
      CALL section_vals_val_get(mp2_section, "SCALE_S", r_val=mp2_env%scale_S)
      CALL section_vals_val_get(mp2_section, "SCALE_T", r_val=mp2_env%scale_T)
      CALL section_vals_val_get(mp2_section, "GROUP_SIZE", i_val=mp2_env%mp2_num_proc)
      CALL section_vals_val_get(mp2_section, "E_GAP", r_val=mp2_env%e_gap)
      CALL section_vals_val_get(mp2_section, "E_RANGE", r_val=mp2_env%e_range)

      CALL section_vals_val_get(mp2_section, "MP2%_SECTION_PARAMETERS_", l_val=do_mp2)
      CALL section_vals_val_get(mp2_section, "MP2%BIG_SEND", l_val=mp2_env%direct_canonical%big_send)
      IF (do_mp2) THEN
         CALL check_method(mp2_env%method)
         CALL section_vals_val_get(mp2_section, "MP2%METHOD", i_val=mp2_env%method)
      END IF

      CALL section_vals_val_get(mp2_section, "RI_RPA%_SECTION_PARAMETERS_", l_val=do_rpa)
      IF (do_rpa) THEN
         CALL check_method(mp2_env%method)
         mp2_env%method = ri_rpa_method_gpw
      END IF
      CALL section_vals_val_get(mp2_section, "RI_RPA%RPA_NUM_QUAD_POINTS", i_val=mp2_env%ri_rpa%rpa_num_quad_points)
      CALL section_vals_val_get(mp2_section, "RI_RPA%NUM_INTEG_GROUPS", i_val=mp2_env%ri_rpa%rpa_num_integ_groups)
      CALL section_vals_val_get(mp2_section, "RI_RPA%MM_STYLE", i_val=mp2_env%ri_rpa%mm_style)
      CALL section_vals_val_get(mp2_section, "RI_RPA%MINIMAX_QUADRATURE", l_val=mp2_env%ri_rpa%minimax_quad)
      CALL section_vals_val_get(mp2_section, "RI_RPA%ADMM", l_val=mp2_env%ri_rpa%do_admm)
      CALL section_vals_val_get(mp2_section, "RI_RPA%SCALE_RPA", r_val=mp2_env%ri_rpa%scale_rpa)
      mp2_env%ri_rpa%reuse_hfx = .FALSE.

      CALL section_vals_val_get(mp2_section, "RI_RPA%EXCHANGE_CORRECTION%_SECTION_PARAMETERS_", &
                                i_val=mp2_env%ri_rpa%exchange_correction)
      CALL section_vals_val_get(mp2_section, "RI_RPA%EXCHANGE_CORRECTION%BLOCK_SIZE", &
                                i_val=mp2_env%ri_rpa%exchange_block_size)
      CALL section_vals_val_get(mp2_section, "RI_RPA%EXCHANGE_CORRECTION%USE_HFX_IMPLEMENTATION", &
                                l_val=mp2_env%ri_rpa%use_hfx_implementation)

      CALL section_vals_val_get(mp2_section, "RI_RPA%RSE", l_val=mp2_env%ri_rpa%do_rse)
      CALL section_vals_val_get(mp2_section, "RI_RPA%PRINT_DGEMM_INFO", l_val=mp2_env%ri_rpa%print_dgemm_info)

      NULLIFY (gw_section)
      gw_section => section_vals_get_subs_vals(mp2_section, "RI_RPA%GW")
      CALL section_vals_val_get(gw_section, "_SECTION_PARAMETERS_", &
                                l_val=mp2_env%ri_rpa%do_ri_g0w0)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%SELF_CONSISTENCY", &
                                i_val=gw_sc)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%CORR_MOS_OCC", &
                                i_val=mp2_env%ri_g0w0%corr_mos_occ)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%CORR_MOS_VIRT", &
                                i_val=mp2_env%ri_g0w0%corr_mos_virt)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%NUMB_POLES", &
                                i_val=mp2_env%ri_g0w0%num_poles)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%NPARAM_PADE", &
                                i_val=mp2_env%ri_g0w0%nparam_pade)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%ANALYTIC_CONTINUATION", &
                                i_val=mp2_env%ri_g0w0%analytic_continuation)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%OMEGA_MAX_FIT", &
                                r_val=mp2_env%ri_g0w0%omega_max_fit)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%CROSSING_SEARCH", &
                                i_val=mp2_env%ri_g0w0%crossing_search)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%FERMI_LEVEL_OFFSET", &
                                r_val=mp2_env%ri_g0w0%fermi_level_offset)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%HEDIN_SHIFT", &
                                l_val=mp2_env%ri_g0w0%do_hedin_shift)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%EV_GW_ITER", &
                                i_val=mp2_env%ri_g0w0%iter_evGW)
      IF (gw_sc == evGW) mp2_env%ri_g0w0%iter_evGW = 10
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%SC_GW0_ITER", &
                                i_val=mp2_env%ri_g0w0%iter_sc_GW0)
      IF (gw_sc == evGW0) mp2_env%ri_g0w0%iter_sc_GW0 = 10
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%EPS_ITER", &
                                r_val=mp2_env%ri_g0w0%eps_iter)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT_EXX", &
                                i_val=mp2_env%ri_g0w0%print_exx)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT_SELF_ENERGY", &
                                l_val=mp2_env%ri_g0w0%print_self_energy)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%RI_SIGMA_X", &
                                l_val=mp2_env%ri_g0w0%do_ri_Sigma_x)

      NULLIFY (r_vals)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC_CORR_LIST", &
                                r_vals=r_vals)
      ALLOCATE (mp2_env%ri_g0w0%ic_corr_list(1)%array(SIZE(r_vals)))
      mp2_env%ri_g0w0%ic_corr_list(1)%array(:) = r_vals(:)
      NULLIFY (r_vals)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC_CORR_LIST_BETA", &
                                r_vals=r_vals)
      ALLOCATE (mp2_env%ri_g0w0%ic_corr_list(2)%array(SIZE(r_vals)))
      mp2_env%ri_g0w0%ic_corr_list(2)%array(:) = r_vals(:)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%GAMMA_ONLY_SIGMA", &
                                l_val=mp2_env%ri_g0w0%do_gamma_only_sigma)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%UPDATE_XC_ENERGY", &
                                l_val=mp2_env%ri_g0w0%update_xc_energy)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%REGULARIZATION_MINIMAX", &
                                r_val=mp2_env%ri_g0w0%regularization_minimax)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%SOC", &
                                i_val=mp2_env%ri_g0w0%soc_type)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%SOC_ENERGY_WINDOW", &
                                r_val=mp2_env%ri_g0w0%soc_energy_window)

      NULLIFY (bse_section)
      bse_section => section_vals_get_subs_vals(mp2_section, "RI_RPA%GW%BSE")
      CALL section_vals_val_get(bse_section, "_SECTION_PARAMETERS_", &
                                l_val=mp2_env%ri_g0w0%do_bse)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%SPIN_CONFIG", &
                                i_val=mp2_env%ri_g0w0%bse_spin_config)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%ENERGY_CUTOFF_OCC", &
                                r_val=mp2_env%ri_g0w0%bse_cutoff_occ)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%ENERGY_CUTOFF_EMPTY", &
                                r_val=mp2_env%ri_g0w0%bse_cutoff_empty)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%NUM_PRINT_EXC", &
                                i_val=mp2_env%ri_g0w0%num_print_exc)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%EPS_X", &
                                r_val=mp2_env%ri_g0w0%eps_x)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%DAVIDSON_ABORT_COND", &
                                i_val=mp2_env%ri_g0w0%davidson_abort_cond)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%NUM_EXC_EN", &
                                i_val=mp2_env%ri_g0w0%num_exc_en)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%NUM_ADD_START_Z_SPACE", &
                                i_val=mp2_env%ri_g0w0%num_add_start_z_space)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%FAC_MAX_Z_SPACE", &
                                i_val=mp2_env%ri_g0w0%fac_max_z_space)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%NUM_NEW_T", &
                                i_val=mp2_env%ri_g0w0%num_new_t)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%EPS_RES", &
                                r_val=mp2_env%ri_g0w0%eps_res)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%EPS_EXC_EN", &
                                r_val=mp2_env%ri_g0w0%eps_exc_en)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%NUM_DAVIDSON_ITER", &
                                i_val=mp2_env%ri_g0w0%num_davidson_iter)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_ITERAT%Z_SPACE_ENERGY_CUTOFF", &
                                r_val=mp2_env%ri_g0w0%z_space_energy_cutoff)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_DIAG_METHOD", &
                                i_val=mp2_env%ri_g0w0%bse_diag_method)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%TDA", &
                                i_val=mp2_env%ri_g0w0%flag_tda)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%BSE_DEBUG_PRINT", &
                                l_val=mp2_env%ri_g0w0%bse_debug_print)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IMAGE_CHARGE_MODEL", &
                                l_val=mp2_env%ri_g0w0%do_ic_model)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC%PRINT_IC_LIST", &
                                l_val=mp2_env%ri_g0w0%print_ic_values)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC%EPS_DIST", &
                                r_val=mp2_env%ri_g0w0%eps_dist)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION", &
                                l_val=mp2_env%ri_g0w0%do_periodic)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%KPOINTS", &
                                i_vals=mp2_env%ri_g0w0%kp_grid)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%NUM_KP_GRIDS", &
                                i_val=mp2_env%ri_g0w0%num_kp_grids)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%EPS_KPOINT", &
                                r_val=mp2_env%ri_g0w0%eps_kpoint)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%MO_COEFF_GAMMA", &
                                l_val=mp2_env%ri_g0w0%do_mo_coeff_gamma)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%AVERAGE_DEGENERATE_LEVELS", &
                                l_val=mp2_env%ri_g0w0%do_average_deg_levels)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%EPS_EIGENVAL", &
                                r_val=mp2_env%ri_g0w0%eps_eigenval)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%EXTRAPOLATE_KPOINTS", &
                                l_val=mp2_env%ri_g0w0%do_extra_kpoints)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%DO_AUX_BAS_GW", &
                                l_val=mp2_env%ri_g0w0%do_aux_bas_gw)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%FRACTION_AUX_MOS", &
                                r_val=mp2_env%ri_g0w0%frac_aux_mos)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC_CORRECTION%NUM_OMEGA_POINTS", &
                                i_val=mp2_env%ri_g0w0%num_omega_points)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%KPOINT_SET%NPOINTS", &
                                i_val=mp2_env%ri_g0w0%n_kp_in_kp_line)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%KPOINT_SET%SPECIAL_POINT", &
                                n_rep_val=mp2_env%ri_g0w0%n_special_kp)
      ALLOCATE (mp2_env%ri_g0w0%xkp_special_kp(3, mp2_env%ri_g0w0%n_special_kp))
      DO i_special_kp = 1, mp2_env%ri_g0w0%n_special_kp
         CALL section_vals_val_get(mp2_section, "RI_RPA%GW%KPOINT_SET%SPECIAL_POINT", &
                                   i_rep_val=i_special_kp, c_vals=string_pointer)
         CPASSERT(SIZE(string_pointer(:), 1) == 3)
         DO i = 1, 3
            CALL read_float_object(string_pointer(i), &
                                   mp2_env%ri_g0w0%xkp_special_kp(i, i_special_kp), &
                                   error_message)
            IF (LEN_TRIM(error_message) > 0) CPABORT(TRIM(error_message))
         END DO
      END DO

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%KPOINTS_SELF_ENERGY", &
                                i_vals=mp2_env%ri_g0w0%kp_grid_Sigma)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%LOWER_BOUND", &
                                r_val=mp2_env%ri_g0w0%dos_lower)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%UPPER_BOUND", &
                                r_val=mp2_env%ri_g0w0%dos_upper)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%STEP", &
                                r_val=mp2_env%ri_g0w0%dos_prec)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%MIN_LEVEL_SPECTRAL", &
                                i_val=mp2_env%ri_g0w0%dos_min)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%MAX_LEVEL_SPECTRAL", &
                                i_val=mp2_env%ri_g0w0%dos_max)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%MIN_LEVEL_SELF_ENERGY", &
                                i_val=mp2_env%ri_g0w0%min_level_self_energy)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%MAX_LEVEL_SELF_ENERGY", &
                                i_val=mp2_env%ri_g0w0%max_level_self_energy)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%GW_DOS%BROADENING", &
                                r_val=mp2_env%ri_g0w0%dos_eta)

      mp2_env%ri_g0w0%do_kpoints_Sigma = mp2_env%ri_g0w0%n_special_kp > 0 .OR. &
                                         mp2_env%ri_g0w0%kp_grid_Sigma(1)* &
                                         mp2_env%ri_g0w0%kp_grid_Sigma(2)* &
                                         mp2_env%ri_g0w0%kp_grid_Sigma(3) > 0

      mp2_env%ri_g0w0%print_local_bandgap = BTEST(cp_print_key_should_output(logger%iter_info, &
                                                                             mp2_section, "RI_RPA%GW%PRINT%LOCAL_BANDGAP"), &
                                                  cp_p_file)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%LOCAL_BANDGAP%ENERGY_WINDOW", &
                                r_val=mp2_env%ri_g0w0%energy_window_print_loc_bandgap)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%LOCAL_BANDGAP%ENERGY_SPACING", &
                                r_val=mp2_env%ri_g0w0%energy_spacing_print_loc_bandgap)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%LOCAL_BANDGAP%LDOS_THRESHOLD_GAP", &
                                r_val=mp2_env%ri_g0w0%ldos_thresh_print_loc_bandgap)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%LOCAL_BANDGAP%STRIDE", &
                                i_vals=mp2_env%ri_g0w0%stride_loc_bandgap)

      NULLIFY (low_scaling_section)
      low_scaling_section => section_vals_get_subs_vals(mp2_section, "LOW_SCALING")
      CALL section_vals_val_get(low_scaling_section, "_SECTION_PARAMETERS_", &
                                l_val=mp2_env%do_im_time)

      CALL section_vals_val_get(low_scaling_section, "MEMORY_CUT", i_val=mp2_env%ri_rpa_im_time%cut_memory)
      CALL section_vals_val_get(low_scaling_section, "MEMORY_INFO", l_val=mp2_env%ri_rpa_im_time%memory_info)
      CALL section_vals_val_get(low_scaling_section, "EPS_FILTER", r_val=mp2_env%ri_rpa_im_time%eps_filter)
      CALL section_vals_val_get(low_scaling_section, "EPS_STORAGE_SCALING", r_val=mp2_env%ri_rpa_im_time%eps_compress)
      mp2_env%ri_rpa_im_time%eps_compress = mp2_env%ri_rpa_im_time%eps_compress*mp2_env%ri_rpa_im_time%eps_filter
      mp2_env%ri_rpa_im_time%eps_compress = MAX(mp2_env%ri_rpa_im_time%eps_compress, 1.0E-16_dp)
      CALL section_vals_val_get(low_scaling_section, "EPS_FILTER_FACTOR", r_val=mp2_env%ri_rpa_im_time%eps_filter_factor)

      CALL section_vals_val_get(low_scaling_section, "DO_KPOINTS", &
                                l_val=mp2_env%ri_rpa_im_time%do_im_time_kpoints)
      CALL section_vals_val_get(low_scaling_section, "KPOINTS", &
                                i_vals=mp2_env%ri_rpa_im_time%kp_grid)
      mp2_env%ri_rpa_im_time%do_kpoints_from_Gamma = SUM(mp2_env%ri_rpa_im_time%kp_grid) > 0
      IF (mp2_env%ri_rpa_im_time%do_kpoints_from_Gamma) THEN
         CPASSERT(mp2_env%ri_g0w0%do_kpoints_Sigma)
      END IF
      CALL section_vals_val_get(low_scaling_section, "KPOINT_WEIGHTS_W", &
                                i_val=mp2_env%ri_rpa_im_time%kpoint_weights_W_method)
      CALL section_vals_val_get(low_scaling_section, "EXPONENT_TAILORED_WEIGHTS", &
                                r_val=mp2_env%ri_rpa_im_time%exp_tailored_weights)
      CALL section_vals_val_get(low_scaling_section, "REGULARIZATION_RI", &
                                r_val=mp2_env%ri_rpa_im_time%regularization_RI)
      CALL section_vals_val_get(low_scaling_section, "EPS_EIGVAL_S", &
                                r_val=mp2_env%ri_rpa_im_time%eps_eigval_S)
      CALL section_vals_val_get(low_scaling_section, "EPS_EIGVAL_S_GAMMA", &
                                r_val=mp2_env%ri_rpa_im_time%eps_eigval_S_Gamma)
      CALL section_vals_val_get(low_scaling_section, "MAKE_CHI_POS_DEFINITE", &
                                l_val=mp2_env%ri_rpa_im_time%make_chi_pos_definite)
      CALL section_vals_val_get(low_scaling_section, "MAKE_OVERLAP_MAT_AO_POS_DEFINITE", &
                                l_val=mp2_env%ri_rpa_im_time%make_overlap_mat_ao_pos_definite)
      CALL section_vals_val_get(low_scaling_section, "TRUNC_COULOMB_RI_X", &
                                l_val=mp2_env%ri_rpa_im_time%trunc_coulomb_ri_x)
      CALL section_vals_val_get(low_scaling_section, "DO_EXTRAPOLATE_KPOINTS", &
                                l_val=mp2_env%ri_rpa_im_time%do_extrapolate_kpoints)
      CALL section_vals_val_get(low_scaling_section, "REL_CUTOFF_TRUNC_COULOMB_RI_X", &
                                r_val=mp2_env%ri_rpa_im_time%rel_cutoff_trunc_coulomb_ri_x)
      CALL section_vals_val_get(low_scaling_section, "K_MESH_G_FACTOR", &
                                i_val=mp2_env%ri_rpa_im_time%k_mesh_g_factor)

      CALL section_vals_val_get(low_scaling_section, "KEEP_QUADRATURE", &
                                l_val=mp2_env%ri_rpa_im_time%keep_quad)
      NULLIFY (mp2_env%ri_rpa_im_time%tau_tj)
      NULLIFY (mp2_env%ri_rpa_im_time%tau_wj)
      NULLIFY (mp2_env%ri_rpa_im_time%tj)
      NULLIFY (mp2_env%ri_rpa_im_time%wj)
      NULLIFY (mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w)
      NULLIFY (mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t)

      CALL section_vals_val_get(low_scaling_section, "MIN_BLOCK_SIZE", &
                                i_val=mp2_env%ri_rpa_im_time%min_bsize)

      CALL section_vals_val_get(low_scaling_section, "MIN_BLOCK_SIZE_MO", &
                                i_val=mp2_env%ri_rpa_im_time%min_bsize_mo)

      CALL section_vals_val_get(mp2_section, "RI_SOS_MP2%_SECTION_PARAMETERS_", l_val=do_ri_sos_mp2)
      IF (do_ri_sos_mp2) THEN
         CALL check_method(mp2_env%method)
         mp2_env%method = ri_mp2_laplace
      END IF
      CALL section_vals_val_get(mp2_section, "RI_SOS_MP2%QUADRATURE_POINTS", i_val=mp2_env%ri_laplace%n_quadrature)
      CALL section_vals_val_get(mp2_section, "RI_SOS_MP2%NUM_INTEG_GROUPS", i_val=mp2_env%ri_laplace%num_integ_groups)

      CALL section_vals_val_get(mp2_section, "RI_MP2%_SECTION_PARAMETERS_", l_val=do_ri_mp2)
      IF (do_ri_mp2) THEN
         CALL check_method(mp2_env%method)
         mp2_env%method = ri_mp2_method_gpw
      END IF
      CALL section_vals_val_get(mp2_section, "RI_MP2%BLOCK_SIZE", i_val=mp2_env%ri_mp2%block_size)
      CALL section_vals_val_get(mp2_section, "RI_MP2%NUMBER_INTEGRATION_GROUPS", i_val=mp2_env%ri_mp2%number_integration_groups)
      CALL section_vals_val_get(mp2_section, "RI_MP2%PRINT_DGEMM_INFO", l_val=mp2_env%ri_mp2%print_dgemm_info)

      CALL section_vals_val_get(mp2_section, "RI%ROW_BLOCK", i_val=mp2_env%block_size_row)
      CALL section_vals_val_get(mp2_section, "RI%COL_BLOCK", i_val=mp2_env%block_size_col)
      CALL section_vals_val_get(mp2_section, "RI%CALC_COND_NUM", l_val=mp2_env%calc_PQ_cond_num)
      CALL section_vals_val_get(mp2_section, "RI%DO_SVD", l_val=mp2_env%do_svd)
      CALL section_vals_val_get(mp2_section, "RI%ERI_BLKSIZE", i_vals=mp2_env%eri_blksize)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%POTENTIAL_TYPE", i_val=mp2_env%ri_metric%potential_type)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%OMEGA", r_val=mp2_env%ri_metric%omega)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%EPS_RANGE", r_val=mp2_env%eps_range)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%CUTOFF_RADIUS", r_val=mp2_env%ri_metric%cutoff_radius)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%T_C_G_DATA", c_val=mp2_env%ri_metric%filename)
      IF (mp2_env%ri_metric%potential_type == do_potential_short) THEN
         CALL erfc_cutoff(mp2_env%eps_range, mp2_env%ri_metric%omega, mp2_env%ri_metric%cutoff_radius)
      END IF

      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%_SECTION_PARAMETERS_", l_val=do_opt_ri_basis)
      IF (do_opt_ri_basis) THEN
         CALL check_method(mp2_env%method)
         mp2_env%method = mp2_ri_optimize_basis
      END IF
      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%DELTA_I_REL", &
                                r_val=mp2_env%ri_opt_param%DI_rel)
      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%DELTA_RI", &
                                r_val=mp2_env%ri_opt_param%DRI)
      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%EPS_DERIV", &
                                r_val=mp2_env%ri_opt_param%eps_step)
      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%MAX_ITER", &
                                i_val=mp2_env%ri_opt_param%max_num_iter)
      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%BASIS_SIZE", &
                                i_val=mp2_env%ri_opt_param%basis_quality)
      NULLIFY (tmplist)
      CALL section_vals_val_get(mp2_section, "RI%OPT_RI_BASIS%NUM_FUNC", &
                                i_vals=tmplist)
      IF (tmplist(1) > 0) THEN
         ALLOCATE (mp2_env%ri_opt_param%RI_nset_per_l(0:SIZE(tmplist) - 1))
         mp2_env%ri_opt_param%RI_nset_per_l = 0
         DO ival = 1, SIZE(tmplist)
            mp2_env%ri_opt_param%RI_nset_per_l(ival - 1) = tmplist(ival)
         END DO
      END IF

      CALL section_vals_val_get(mp2_section, "INTEGRALS%ERI_METHOD", i_val=mp2_env%eri_method)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%SIZE_LATTICE_SUM", i_val=mp2_env%mp2_gpw%size_lattice_sum)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%WFC_GPW%EPS_FILTER", r_val=mp2_env%mp2_gpw%eps_filter)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%WFC_GPW%EPS_GRID", r_val=mp2_env%mp2_gpw%eps_grid)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%WFC_GPW%CUTOFF", r_val=mp2_env%mp2_gpw%cutoff)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%WFC_GPW%REL_CUTOFF", r_val=mp2_env%mp2_gpw%relative_cutoff)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%WFC_GPW%PRINT_LEVEL", i_val=mp2_env%mp2_gpw%print_level)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%WFC_GPW%EPS_PGF_ORB_S", r_val=mp2_env%mp2_gpw%eps_pgf_orb_S)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%POTENTIAL_TYPE", &
                                i_val=mp2_env%potential_parameter%potential_type)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%TRUNCATION_RADIUS", &
                                r_val=mp2_env%potential_parameter%cutoff_radius)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%POTENTIAL_DATA", &
                                c_val=mp2_env%potential_parameter%filename)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%OMEGA", &
                                r_val=mp2_env%potential_parameter%omega)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%SCALE_COULOMB", &
                                r_val=mp2_env%potential_parameter%scale_coulomb)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%SCALE_LONGRANGE", &
                                r_val=mp2_env%potential_parameter%scale_longrange)

      NULLIFY (mp2_env%eri_mme_param)
      ALLOCATE (mp2_env%eri_mme_param)

      IF (mp2_env%eri_method .EQ. do_eri_mme) THEN
         eri_mme_section => section_vals_get_subs_vals(mp2_section, "INTEGRALS%ERI_MME")
         CALL cp_eri_mme_init_read_input(eri_mme_section, mp2_env%eri_mme_param)
      END IF

      ! Set some parameters in case of P screening
      mp2_env%not_last_hfx = .TRUE.
      mp2_env%p_screen = .TRUE.

      ! Set the CPHF section
      CALL section_vals_val_get(mp2_section, "CANONICAL_GRADIENTS%FREE_HFX_BUFFER", l_val=mp2_env%ri_grad%free_hfx_buffer)
      CALL section_vals_val_get(mp2_section, "CANONICAL_GRADIENTS%EPS_CANONICAL", r_val=mp2_env%ri_grad%eps_canonical)
      CALL section_vals_val_get(mp2_section, "CANONICAL_GRADIENTS%DOT_PRODUCT_BLKSIZE", i_val=mp2_env%ri_grad%dot_blksize)
      CALL section_vals_val_get(mp2_section, "CANONICAL_GRADIENTS%MAX_PARALLEL_COMM", i_val=mp2_env%ri_grad%max_parallel_comm)
      cphf_section => section_vals_get_subs_vals(mp2_section, "CANONICAL_GRADIENTS%CPHF")
      IF (ASSOCIATED(cphf_section)) THEN
         CALL section_vals_val_get(cphf_section, "MAX_ITER", i_val=mp2_env%ri_grad%cphf_max_num_iter)
         CALL section_vals_val_get(cphf_section, "EPS_CONV", r_val=mp2_env%ri_grad%cphf_eps_conv)
         CALL section_vals_val_get(cphf_section, "SCALE_STEP_SIZE", r_val=mp2_env%ri_grad%scale_step_size)
         CALL section_vals_val_get(cphf_section, "SOLVER_METHOD", i_val=mp2_env%ri_grad%z_solver_method)
         CALL section_vals_val_get(cphf_section, "RESTART_EVERY", i_val=mp2_env%ri_grad%cphf_restart)
         CALL section_vals_val_get(cphf_section, "ENFORCE_DECREASE", l_val=mp2_env%ri_grad%enforce_decrease)
         CALL section_vals_val_get(cphf_section, "RECALC_RESIDUAL", l_val=mp2_env%ri_grad%recalc_residual)
         CALL section_vals_val_get(cphf_section, "DO_POLAK_RIBIERE", l_val=mp2_env%ri_grad%polak_ribiere)
      END IF

      ! print some info about the MP2 parameters
      unit_nr = cp_print_key_unit_nr(logger, mp2_section, "PRINT", &
                                     extension=".mp2Log")
      IF ((mp2_env%method .NE. mp2_method_none) .AND. unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') ""
         SELECT CASE (mp2_env%method)
         CASE (mp2_method_direct)
            WRITE (unit_nr, '(T2,A)') "MP2| using direct canonical MP2"
         CASE (mp2_method_gpw)
            WRITE (unit_nr, '(T2,A)') "MP2| using MP2 GPW style"
         CASE (ri_mp2_method_gpw)
            WRITE (unit_nr, '(T2,A)') "MP2| using RI-MP2-GPW"
         CASE (ri_rpa_method_gpw)
            WRITE (unit_nr, '(T2,A)') "RI-RPA| using GPW style"
         CASE (ri_mp2_laplace)
            WRITE (unit_nr, '(T2,A)') "RI-SOS-Laplace-MP2| using GPW style"
         CASE (mp2_ri_optimize_basis)
            WRITE (unit_nr, '(T2,A)') "MP2| Optimize RI auxiliary basis"
         CASE DEFAULT
            CPABORT("")
         END SELECT
         WRITE (unit_nr, '(T2,A)') ""
         CALL m_flush(unit_nr)
      END IF

      CALL cp_print_key_finished_output(unit_nr, logger, mp2_section, &
                                        "PRINT")

      CALL timestop(handle)

   END SUBROUTINE read_mp2_section

! **************************************************************************************************
!> \brief ...
!> \param method ...
! **************************************************************************************************
   SUBROUTINE check_method(method)
      INTEGER, INTENT(IN)                                :: method

      CHARACTER(len=*), PARAMETER                        :: routineN = 'check_method'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (method .NE. mp2_method_none) THEN
         CPABORT("Please use not more than one method to compute the correlation energy.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_method
END MODULE mp2_setup
