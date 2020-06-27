!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: &
        do_eri_mme, do_potential_short, mp2_method_direct, mp2_method_gpw, mp2_method_none, &
        mp2_ri_optimize_basis, ri_mp2_laplace, ri_mp2_method_gpw, ri_rpa_method_gpw
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
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

      CHARACTER(len=*), PARAMETER :: routineN = 'read_mp2_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ival, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: do_mp2, do_opt_ri_basis, do_ri_mp2, &
                                                            do_ri_sos_mp2, do_rpa
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: eri_mme_section, gw_section, &
                                                            im_time_test_section, &
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
      CALL section_vals_val_get(mp2_section, "RI_RPA%SIZE_FREQ_INTEG_GROUP", i_val=mp2_env%ri_rpa%rpa_integ_group_size)
      CALL section_vals_val_get(mp2_section, "RI_RPA%MM_STYLE", i_val=mp2_env%ri_rpa%mm_style)
      CALL section_vals_val_get(mp2_section, "RI_RPA%MINIMAX_QUADRATURE", l_val=mp2_env%ri_rpa%minimax_quad)
      CALL section_vals_val_get(mp2_section, "RI_RPA%ADMM", l_val=mp2_env%ri_rpa%do_admm)
      CALL section_vals_val_get(mp2_section, "RI_RPA%SCALE_RPA", r_val=mp2_env%ri_rpa%scale_rpa)

      CALL section_vals_val_get(mp2_section, "RI_RPA%AXK", l_val=mp2_env%ri_rpa%do_ri_axk)

      CALL section_vals_val_get(mp2_section, "RI_RPA%RSE", l_val=mp2_env%ri_rpa%do_rse)

      NULLIFY (gw_section)
      gw_section => section_vals_get_subs_vals(mp2_section, "RI_RPA%GW")
      CALL section_vals_val_get(gw_section, "_SECTION_PARAMETERS_", &
                                l_val=mp2_env%ri_rpa%do_ri_g0w0)
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
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%EV_GW_ITER", &
                                i_val=mp2_env%ri_g0w0%iter_evGW)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%SC_GW0_ITER", &
                                i_val=mp2_env%ri_g0w0%iter_sc_GW0)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%EPS_ITER", &
                                r_val=mp2_env%ri_g0w0%eps_iter)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT_EXX", &
                                i_val=mp2_env%ri_g0w0%print_exx)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%RI_SIGMA_X", &
                                l_val=mp2_env%ri_g0w0%do_ri_Sigma_x)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC_CORR_LIST", &
                                r_vals=mp2_env%ri_g0w0%ic_corr_list)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC_CORR_LIST_BETA", &
                                r_vals=mp2_env%ri_g0w0%ic_corr_list_beta)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%GAMMA_ONLY_SIGMA", &
                                l_val=mp2_env%ri_g0w0%do_gamma_only_sigma)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%UPDATE_XC_ENERGY", &
                                l_val=mp2_env%ri_g0w0%update_xc_energy)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE", &
                                l_val=mp2_env%ri_g0w0%do_bse)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%NUM_Z_VECTORS", &
                                i_val=mp2_env%ri_g0w0%num_z_vectors)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%THRESHOLD_MIN_TRANS", &
                                r_val=mp2_env%ri_g0w0%eps_min_trans)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%BSE%MAX_ITER", &
                                i_val=mp2_env%ri_g0w0%max_iter_bse)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IMAGE_CHARGE_MODEL", &
                                l_val=mp2_env%ri_g0w0%do_ic_model)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC%PRINT_IC_LIST", &
                                l_val=mp2_env%ri_g0w0%print_ic_values)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%IC%EPS_DIST", &
                                r_val=mp2_env%ri_g0w0%eps_dist)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC", &
                                l_val=mp2_env%ri_g0w0%do_periodic)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%KPOINTS", &
                                i_vals=mp2_env%ri_g0w0%kp_grid)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%NUM_KP_GRIDS", &
                                i_val=mp2_env%ri_g0w0%num_kp_grids)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%EPS_KPOINT", &
                                r_val=mp2_env%ri_g0w0%eps_kpoint)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%MO_COEFF_GAMMA", &
                                l_val=mp2_env%ri_g0w0%do_mo_coeff_gamma)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%AVERAGE_DEGENERATE_LEVELS", &
                                l_val=mp2_env%ri_g0w0%do_average_deg_levels)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%EPS_EIGENVAL", &
                                r_val=mp2_env%ri_g0w0%eps_eigenval)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%EXTRAPOLATE_KPOINTS", &
                                l_val=mp2_env%ri_g0w0%do_extra_kpoints)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%DO_AUX_BAS_GW", &
                                l_val=mp2_env%ri_g0w0%do_aux_bas_gw)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%FRACTION_AUX_MOS", &
                                r_val=mp2_env%ri_g0w0%frac_aux_mos)
      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PERIODIC%NUM_OMEGA_POINTS", &
                                i_val=mp2_env%ri_g0w0%num_omega_points)

      CALL section_vals_val_get(mp2_section, "RI_RPA%GW%PRINT%SELF_ENERGY", &
                                l_val=mp2_env%ri_g0w0%print_self_energy)

      NULLIFY (low_scaling_section)
      low_scaling_section => section_vals_get_subs_vals(mp2_section, "LOW_SCALING")
      CALL section_vals_val_get(low_scaling_section, "_SECTION_PARAMETERS_", &
                                l_val=mp2_env%do_im_time)

      CALL section_vals_val_get(low_scaling_section, "MEMORY_CUT", i_val=mp2_env%ri_rpa_im_time%cut_memory)
      CALL section_vals_val_get(low_scaling_section, "MEMORY_INFO", l_val=mp2_env%ri_rpa_im_time%memory_info)
      CALL section_vals_val_get(low_scaling_section, "EPS_FILTER", r_val=mp2_env%ri_rpa_im_time%eps_filter)
      CALL section_vals_val_get(low_scaling_section, "EPS_FILTER_FACTOR", r_val=mp2_env%ri_rpa_im_time%eps_filter_factor)

      CALL section_vals_val_get(low_scaling_section, "DO_KPOINTS", &
                                l_val=mp2_env%ri_rpa_im_time%do_im_time_kpoints)
      CALL section_vals_val_get(low_scaling_section, "CUTOFF_W", &
                                r_val=mp2_env%ri_rpa_im_time%cutoff)
      CALL section_vals_val_get(low_scaling_section, "KPOINTS", &
                                i_vals=mp2_env%ri_rpa_im_time%kp_grid)
      CALL section_vals_val_get(low_scaling_section, "EXP_KPOINTS", &
                                r_val=mp2_env%ri_rpa_im_time%exp_kpoints)

      CALL section_vals_val_get(low_scaling_section, "MIN_BLOCK_SIZE", &
                                i_val=mp2_env%ri_rpa_im_time%min_bsize)

      NULLIFY (im_time_test_section)
      im_time_test_section => section_vals_get_subs_vals(low_scaling_section, "TEST")

      CALL section_vals_val_get(im_time_test_section, "_SECTION_PARAMETERS_", &
                                l_val=mp2_env%ri_rpa_im_time%group_size_internal)

      mp2_env%ri_rpa_im_time%group_size_internal = .NOT. mp2_env%ri_rpa_im_time%group_size_internal

      IF (.NOT. mp2_env%ri_rpa_im_time%group_size_internal) THEN
         CALL section_vals_val_get(im_time_test_section, "GROUP_SIZE_3c", &
                                   i_val=mp2_env%ri_rpa_im_time%group_size_3c)
         CALL section_vals_val_get(im_time_test_section, "GROUP_SIZE_P", &
                                   i_val=mp2_env%ri_rpa_im_time%group_size_P)
      ENDIF

      CALL section_vals_val_get(mp2_section, "RI_SOS_MP2%_SECTION_PARAMETERS_", l_val=do_ri_sos_mp2)
      IF (do_ri_sos_mp2) THEN
         CALL check_method(mp2_env%method)
         mp2_env%method = ri_mp2_laplace
      END IF
      CALL section_vals_val_get(mp2_section, "RI_SOS_MP2%QUADRATURE_POINTS", i_val=mp2_env%ri_laplace%n_quadrature)
      CALL section_vals_val_get(mp2_section, "RI_SOS_MP2%SIZE_INTEG_GROUP", i_val=mp2_env%ri_laplace%integ_group_size)

      CALL section_vals_val_get(mp2_section, "RI_MP2%_SECTION_PARAMETERS_", l_val=do_ri_mp2)
      IF (do_ri_mp2) THEN
         CALL check_method(mp2_env%method)
         mp2_env%method = ri_mp2_method_gpw
      END IF
      CALL section_vals_val_get(mp2_section, "RI_MP2%BLOCK_SIZE", i_val=mp2_env%ri_mp2%block_size)
      CALL section_vals_val_get(mp2_section, "RI_MP2%EPS_CANONICAL", r_val=mp2_env%ri_mp2%eps_canonical)
      CALL section_vals_val_get(mp2_section, "RI_MP2%FREE_HFX_BUFFER", l_val=mp2_env%ri_mp2%free_hfx_buffer)

      CALL section_vals_val_get(mp2_section, "RI_MP2%CPHF%MAX_ITER", i_val=mp2_env%ri_grad%cphf_max_num_iter)
      CALL section_vals_val_get(mp2_section, "RI_MP2%CPHF%EPS_CONV", r_val=mp2_env%ri_grad%cphf_eps_conv)

      CALL section_vals_val_get(mp2_section, "RI%ROW_BLOCK", i_val=mp2_env%block_size_row)
      CALL section_vals_val_get(mp2_section, "RI%COL_BLOCK", i_val=mp2_env%block_size_col)
      CALL section_vals_val_get(mp2_section, "RI%CALC_COND_NUM", l_val=mp2_env%calc_PQ_cond_num)
      CALL section_vals_val_get(mp2_section, "RI%DO_SVD", l_val=mp2_env%do_svd)
      CALL section_vals_val_get(mp2_section, "RI%EPS_SVD", r_val=mp2_env%eps_svd)
      CALL section_vals_val_get(mp2_section, "RI%ERI_BLKSIZE", i_vals=mp2_env%eri_blksize)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%POTENTIAL_TYPE", i_val=mp2_env%ri_metric%potential_type)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%OMEGA", r_val=mp2_env%ri_metric%omega)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%EPS_RANGE", r_val=mp2_env%eps_range)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%CUTOFF_RADIUS", r_val=mp2_env%ri_metric%cutoff_radius)
      CALL section_vals_val_get(mp2_section, "RI%RI_METRIC%T_C_G_DATA", c_val=mp2_env%ri_metric%filename)
      IF (mp2_env%ri_metric%potential_type == do_potential_short) THEN
         CALL erfc_cutoff(mp2_env%eps_range, mp2_env%ri_metric%omega, mp2_env%ri_metric%cutoff_radius)
      ENDIF

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
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%TShPSC_DATA", &
                                c_val=mp2_env%potential_parameter%filename)
      CALL section_vals_val_get(mp2_section, "INTEGRALS%INTERACTION_POTENTIAL%OMEGA", &
                                r_val=mp2_env%potential_parameter%omega)

      NULLIFY (mp2_env%eri_mme_param)
      ALLOCATE (mp2_env%eri_mme_param)

      IF (mp2_env%eri_method .EQ. do_eri_mme) THEN
         eri_mme_section => section_vals_get_subs_vals(mp2_section, "INTEGRALS%ERI_MME")
         CALL cp_eri_mme_init_read_input(eri_mme_section, mp2_env%eri_mme_param)
      ENDIF

      ! Set some parameters in case of P screening
      mp2_env%not_last_hfx = .TRUE.
      mp2_env%p_screen = .TRUE.

      ! print some info about the MP2 parameters
      unit_nr = cp_print_key_unit_nr(logger, mp2_section, "MP2%MP2_INFO", &
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
      ENDIF

      CALL cp_print_key_finished_output(unit_nr, logger, mp2_section, &
                                        "MP2%MP2_INFO")

      CALL timestop(handle)

   END SUBROUTINE read_mp2_section

! **************************************************************************************************
!> \brief ...
!> \param method ...
! **************************************************************************************************
   SUBROUTINE check_method(method)
      INTEGER, INTENT(IN)                                :: method

      CHARACTER(len=*), PARAMETER :: routineN = 'check_method', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (method .NE. mp2_method_none) THEN
         CPABORT("Please use not more than one method to compute the correlation energy.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_method
END MODULE mp2_setup
