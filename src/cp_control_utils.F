!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities to set up the control types
! **************************************************************************************************
MODULE cp_control_utils
   USE bibliography,                    ONLY: &
        Andreussi2012, Dewar1977, Dewar1985, Elstner1998, Fattebert2002, Grimme2017, Hu2007, &
        Krack2000, Lippert1997, Lippert1999, Porezag1995, Repasky2002, Rocha2006, Schenter2008, &
        Seifert1996, Souza2002, Stengel2009, Stewart1989, Stewart2007, Thiel1992, Umari2002, &
        VanVoorhis2015, VandeVondele2005a, VandeVondele2005b, Yin2017, Zhechkov2005, cite_reference
   USE cp_control_types,                ONLY: &
        admm_control_create, admm_control_type, ddapc_control_create, ddapc_restraint_type, &
        dft_control_create, dft_control_type, efield_type, qs_control_type, sccs_control_create, &
        tddfpt2_control_type, tddfpt_control_create, tddfpt_control_type, xtb_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                              cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        constant_env, custom_env, do_admm_basis_projection, do_admm_blocked_projection, &
        do_admm_blocking_purify_full, do_admm_charge_constrained_projection, &
        do_admm_exch_scaling_merlot, do_admm_purify_mcweeny, do_admm_purify_mo_diag, &
        do_admm_purify_mo_no_diag, do_admm_purify_none, do_admm_purify_none_dm, &
        do_ddapc_constraint, do_ddapc_restraint, do_method_am1, do_method_dftb, do_method_gapw, &
        do_method_gapw_xc, do_method_gpw, do_method_lrigpw, do_method_mndo, do_method_mndod, &
        do_method_ofgpw, do_method_pdg, do_method_pm3, do_method_pm6, do_method_pm6fm, &
        do_method_pnnl, do_method_rigpw, do_method_rm1, do_method_xtb, do_pwgrid_ns_fullspace, &
        do_pwgrid_ns_halfspace, do_pwgrid_spherical, do_s2_constraint, do_s2_restraint, &
        do_se_is_kdso, do_se_is_kdso_d, do_se_is_slater, do_se_lr_ewald, do_se_lr_ewald_gks, &
        do_se_lr_ewald_r3, do_se_lr_none, gaussian_env, numerical, ramp_env, &
        real_time_propagation, sccs_andreussi, sccs_derivative_cd3, sccs_derivative_cd5, &
        sccs_derivative_cd7, sccs_derivative_fft, sccs_fattebert_gygi, sic_ad, sic_eo, &
        sic_list_all, sic_list_unpaired, sic_mauri_spz, sic_mauri_us, sic_none, slater, &
        tddfpt_excitations, use_mom_ref_user, xas_dip_len
   USE input_cp2k_check,                ONLY: xc_functionals_expand
   USE input_cp2k_dft,                  ONLY: create_dft_section
   USE input_enumeration_types,         ONLY: enum_i2c,&
                                              enumeration_type
   USE input_keyword_types,             ONLY: keyword_get,&
                                              keyword_type
   USE input_section_types,             ONLY: &
        section_get_ival, section_get_keyword, section_release, section_type, section_vals_get, &
        section_vals_get_subs_vals, section_vals_type, section_vals_val_get, section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_cdft_utils,                   ONLY: read_cdft_control_section
   USE string_utilities,                ONLY: uppercase
   USE util,                            ONLY: sort
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_input_constants,              ONLY: xc_deriv_collocate
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_write_output,                 ONLY: xc_write
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_control_utils'

   PUBLIC :: read_dft_control, &
             read_mgrid_section, &
             read_qs_section, &
             read_tddfpt_control, &
             read_tddfpt2_control, &
             write_dft_control, &
             write_qs_control, &
             read_ddapc_section
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE read_dft_control(dft_control, dft_section)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_dft_control', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_path_length)                 :: basis_set_file_name, potential_file_name
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: tmpstringlist
      INTEGER                                            :: excitations, irep, isize, method_id, &
                                                            nrep, xc_deriv_method_id
      LOGICAL                                            :: do_ot, do_rtp, explicit, is_present, &
                                                            l_param, not_SE, was_present
      REAL(KIND=dp)                                      :: density_cut, gradient_cut, tau_cut
      REAL(KIND=dp), DIMENSION(:), POINTER               :: pol
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: sccs_section, scf_section, tmp_section, &
                                                            xc_fun_section, xc_section
      TYPE(xc_rho_cflags_type)                           :: needs

      was_present = .FALSE.

      logger => cp_get_default_logger()

      NULLIFY (tmp_section, xc_fun_section, xc_section)
      CALL dft_control_create(dft_control)
      ! determine wheather this is a semiempirical or DFTB run
      ! --> (no XC section needs to be provided)
      not_SE = .TRUE.
      CALL section_vals_val_get(dft_section, "QS%METHOD", i_val=method_id)
      SELECT CASE (method_id)
      CASE (do_method_dftb, do_method_xtb, do_method_mndo, do_method_am1, do_method_pm3, do_method_pnnl, &
            do_method_pm6, do_method_pm6fm, do_method_pdg, do_method_rm1, do_method_mndod)
         not_SE = .FALSE.
      END SELECT
      ! Check for XC section and XC_FUNCTIONAL section
      xc_section => section_vals_get_subs_vals(dft_section, "XC")
      CALL section_vals_get(xc_section, explicit=is_present)
      IF (.NOT. is_present .AND. not_SE) THEN
         CPABORT("XC section missing.")
      END IF
      IF (is_present) THEN
         CALL section_vals_val_get(xc_section, "density_cutoff", r_val=density_cut)
         CALL section_vals_val_get(xc_section, "gradient_cutoff", r_val=gradient_cut)
         CALL section_vals_val_get(xc_section, "tau_cutoff", r_val=tau_cut)
         ! Perform numerical stability checks and possibly correct the issues
         IF (density_cut <= EPSILON(0.0_dp)*100.0_dp) &
            CALL cp_warn(__LOCATION__, &
                         "DENSITY_CUTOFF lower than 100*EPSILON, where EPSILON is the machine precision. "// &
                         "This may lead to numerical problems. Setting up shake_tol to 100*EPSILON! ")
         density_cut = MAX(EPSILON(0.0_dp)*100.0_dp, density_cut)
         IF (gradient_cut <= EPSILON(0.0_dp)*100.0_dp) &
            CALL cp_warn(__LOCATION__, &
                         "GRADIENT_CUTOFF lower than 100*EPSILON, where EPSILON is the machine precision. "// &
                         "This may lead to numerical problems. Setting up shake_tol to 100*EPSILON! ")
         gradient_cut = MAX(EPSILON(0.0_dp)*100.0_dp, gradient_cut)
         IF (tau_cut <= EPSILON(0.0_dp)*100.0_dp) &
            CALL cp_warn(__LOCATION__, &
                         "TAU_CUTOFF lower than 100*EPSILON, where EPSILON is the machine precision. "// &
                         "This may lead to numerical problems. Setting up shake_tol to 100*EPSILON! ")
         tau_cut = MAX(EPSILON(0.0_dp)*100.0_dp, tau_cut)
         CALL section_vals_val_set(xc_section, "density_cutoff", r_val=density_cut)
         CALL section_vals_val_set(xc_section, "gradient_cutoff", r_val=gradient_cut)
         CALL section_vals_val_set(xc_section, "tau_cutoff", r_val=tau_cut)
      END IF
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_get(xc_fun_section, explicit=is_present)
      IF (.NOT. is_present .AND. not_SE) THEN
         CPABORT("XC_FUNCTIONAL section missing.")
      END IF
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      CALL section_vals_val_get(dft_section, "UKS", l_val=dft_control%uks)
      CALL section_vals_val_get(dft_section, "ROKS", l_val=dft_control%roks)
      IF (dft_control%uks .OR. dft_control%roks) THEN
         dft_control%nspins = 2
      ELSE
         dft_control%nspins = 1
      END IF

      dft_control%lsd = (dft_control%nspins > 1)
      needs = xc_functionals_get_needs(xc_fun_section, &
                                       lsd=dft_control%lsd, &
                                       add_basic_components=.TRUE.)
      dft_control%use_kinetic_energy_density = (needs%tau_spin .OR. needs%tau)

      xc_deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      dft_control%drho_by_collocation = (needs%norm_drho .AND. (xc_deriv_method_id == xc_deriv_collocate))
      IF (dft_control%drho_by_collocation) THEN
         CPABORT("derivatives by collocation not implemented")
      END IF

      ! Automatic auxiliary basis set generation
      CALL section_vals_val_get(dft_section, "AUTO_BASIS", n_rep_val=nrep)
      DO irep = 1, nrep
         CALL section_vals_val_get(dft_section, "AUTO_BASIS", i_rep_val=irep, c_vals=tmpstringlist)
         IF (SIZE(tmpstringlist) == 2) THEN
            CALL uppercase(tmpstringlist(2))
            SELECT CASE (tmpstringlist(2))
            CASE ("X")
               isize = -1
            CASE ("SMALL")
               isize = 0
            CASE ("MEDIUM")
               isize = 1
            CASE ("LARGE")
               isize = 2
            CASE ("HUGE")
               isize = 3
            CASE DEFAULT
               CPWARN("Unknown basis size in AUTO_BASIS keyword:"//TRIM(tmpstringlist(1)))
            END SELECT
            !
            SELECT CASE (tmpstringlist(1))
            CASE ("X")
            CASE ("RI_AUX")
               dft_control%auto_basis_ri_aux = isize
            CASE ("AUX_FIT")
               dft_control%auto_basis_aux_fit = isize
            CASE ("LRI_AUX")
               dft_control%auto_basis_lri_aux = isize
            CASE ("RI_HXC")
               dft_control%auto_basis_ri_hxc = isize
            CASE ("RI_XAS")
               dft_control%auto_basis_ri_xas = isize
            CASE ("RI_HFX")
               dft_control%auto_basis_ri_hfx = isize
            CASE DEFAULT
               CPWARN("Unknown basis type in AUTO_BASIS keyword:"//TRIM(tmpstringlist(1)))
            END SELECT
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "AUTO_BASIS keyword in &DFT section has a wrong number of arguments.")
         END IF
      END DO

      !! check if we do wavefunction fitting
      tmp_section => section_vals_get_subs_vals(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD")
      CALL section_vals_get(tmp_section, explicit=is_present)
      dft_control%do_admm = is_present
      dft_control%do_admm_mo = .FALSE.
      dft_control%do_admm_dm = .FALSE.
      IF (is_present) THEN
         do_ot = .FALSE.
         CALL section_vals_val_get(scf_section, "OT%_SECTION_PARAMETERS_", l_val=do_ot)
         CALL admm_control_create(dft_control%admm_control)

         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%EPS_FILTER", &
                                   r_val=dft_control%admm_control%eps_filter)
         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%ADMM_PURIFICATION_METHOD", i_val=method_id)
         dft_control%admm_control%purification_method = method_id

         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%METHOD", i_val=method_id)
         dft_control%admm_control%method = method_id

         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%EXCH_SCALING_MODEL", i_val=method_id)
         dft_control%admm_control%scaling_model = method_id

         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%EXCH_CORRECTION_FUNC", i_val=method_id)
         dft_control%admm_control%aux_exch_func = method_id

         ! parameters for X functional
         dft_control%admm_control%aux_exch_func_param = .FALSE.
         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%OPTX_A1", explicit=explicit, &
                                   r_val=dft_control%admm_control%aux_x_param(1))
         IF (explicit) dft_control%admm_control%aux_exch_func_param = .TRUE.
         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%OPTX_A2", explicit=explicit, &
                                   r_val=dft_control%admm_control%aux_x_param(2))
         IF (explicit) dft_control%admm_control%aux_exch_func_param = .TRUE.
         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%OPTX_GAMMA", explicit=explicit, &
                                   r_val=dft_control%admm_control%aux_x_param(3))
         IF (explicit) dft_control%admm_control%aux_exch_func_param = .TRUE.

         CALL read_admm_block_list(dft_control%admm_control, dft_section)

         !    In the case of charge-constrained projection (e.g. according to Merlot),
         !    there is no purification needed and hence, do_admm_purify_none has to be set.

         IF ((dft_control%admm_control%method == do_admm_blocking_purify_full .OR. &
              dft_control%admm_control%method == do_admm_blocked_projection) &
             .AND. dft_control%admm_control%scaling_model == do_admm_exch_scaling_merlot) THEN
            CPABORT("ADMM: Blocking and Merlot scaling are mutually exclusive.")
         END IF

         IF (dft_control%admm_control%method == do_admm_charge_constrained_projection .AND. &
             dft_control%admm_control%purification_method /= do_admm_purify_none) THEN
            CALL cp_abort(__LOCATION__, &
                          "ADMM: In the case of METHOD=CHARGE_CONSTRAINED_PROJECTION, "// &
                          "ADMM_PURIFICATION_METHOD=NONE has to be set.")
         END IF

         IF (dft_control%admm_control%purification_method == do_admm_purify_mo_diag .OR. &
             dft_control%admm_control%purification_method == do_admm_purify_mo_no_diag) THEN
            IF (dft_control%admm_control%method /= do_admm_basis_projection) &
               CPABORT("ADMM: Chosen purification requires BASIS_PROJECTION")

            IF (.NOT. do_ot) CPABORT("ADMM: MO-based purification requires OT.")
         END IF

         IF (dft_control%admm_control%purification_method == do_admm_purify_none_dm .OR. &
             dft_control%admm_control%purification_method == do_admm_purify_mcweeny) THEN
            dft_control%do_admm_dm = .TRUE.
         ELSE
            dft_control%do_admm_mo = .TRUE.
         ENDIF
      END IF

      ! Set restricted to true, if both OT and ROKS are requested
      !MK in principle dft_control%restricted could be dropped completely like the
      !MK input key by using only dft_control%roks now
      CALL section_vals_val_get(scf_section, "OT%_SECTION_PARAMETERS_", l_val=l_param)
      dft_control%restricted = (dft_control%roks .AND. l_param)

      CALL section_vals_val_get(dft_section, "CHARGE", i_val=dft_control%charge)
      CALL section_vals_val_get(dft_section, "MULTIPLICITY", i_val=dft_control%multiplicity)
      CALL section_vals_val_get(dft_section, "RELAX_MULTIPLICITY", r_val=dft_control%relax_multiplicity)
      IF (dft_control%relax_multiplicity > 0.0_dp) THEN
         IF (.NOT. dft_control%uks) &
            CALL cp_abort(__LOCATION__, "The option RELAX_MULTIPLICITY is only valid for "// &
                          "unrestricted Kohn-Sham (UKS) calculations")
      END IF

      ! check for the presence of the low spin roks section
      tmp_section => section_vals_get_subs_vals(dft_section, "LOW_SPIN_ROKS")
      CALL section_vals_get(tmp_section, explicit=dft_control%low_spin_roks)

      dft_control%sic_method_id = sic_none
      dft_control%sic_scaling_a = 1.0_dp
      dft_control%sic_scaling_b = 1.0_dp

      ! DFT+U
      dft_control%dft_plus_u = .FALSE.
      CALL section_vals_val_get(dft_section, "PLUS_U_METHOD", i_val=method_id)
      dft_control%plus_u_method_id = method_id

      ! Smearing in use
      dft_control%smear = .FALSE.

      ! Surface dipole correction
      dft_control%correct_surf_dip = .FALSE.
      CALL section_vals_val_get(dft_section, "SURFACE_DIPOLE_CORRECTION", l_val=dft_control%correct_surf_dip)
      CALL section_vals_val_get(dft_section, "SURF_DIP_DIR", i_val=dft_control%dir_surf_dip)
      dft_control%pos_dir_surf_dip = -1.0_dp
      CALL section_vals_val_get(dft_section, "SURF_DIP_POS", r_val=dft_control%pos_dir_surf_dip)
! another logical variable, surf_dip_correct_switch, is introduced for
! implementation of "SURF_DIP_SWITCH" [SGh]
      dft_control%switch_surf_dip = .FALSE.
      dft_control%surf_dip_correct_switch = dft_control%correct_surf_dip
      CALL section_vals_val_get(dft_section, "SURF_DIP_SWITCH", l_val=dft_control%switch_surf_dip)
      dft_control%correct_el_density_dip = .FALSE.
      CALL section_vals_val_get(dft_section, "CORE_CORR_DIP", l_val=dft_control%correct_el_density_dip)
      IF (dft_control%correct_el_density_dip) THEN
         IF (dft_control%correct_surf_dip) THEN
            ! Do nothing, move on
         ELSE
            dft_control%correct_el_density_dip = .FALSE.
            CPWARN("CORE_CORR_DIP keyword is activated only if SURFACE_DIPOLE_CORRECTION is TRUE")
         ENDIF
      ENDIF

      CALL section_vals_val_get(dft_section, "BASIS_SET_FILE_NAME", &
                                c_val=basis_set_file_name)
      CALL section_vals_val_get(dft_section, "POTENTIAL_FILE_NAME", &
                                c_val=potential_file_name)

      ! Read the input section
      tmp_section => section_vals_get_subs_vals(dft_section, "sic")
      CALL section_vals_val_get(tmp_section, "SIC_METHOD", &
                                i_val=dft_control%sic_method_id)
      CALL section_vals_val_get(tmp_section, "ORBITAL_SET", &
                                i_val=dft_control%sic_list_id)
      CALL section_vals_val_get(tmp_section, "SIC_SCALING_A", &
                                r_val=dft_control%sic_scaling_a)
      CALL section_vals_val_get(tmp_section, "SIC_SCALING_B", &
                                r_val=dft_control%sic_scaling_b)

      ! Determine if this is a TDDFPT run
      CALL section_vals_val_get(dft_section, "EXCITATIONS", i_val=excitations)
      dft_control%do_tddfpt_calculation = (excitations == tddfpt_excitations)
      IF (dft_control%do_tddfpt_calculation) THEN
         CALL tddfpt_control_create(dft_control%tddfpt_control)
      END IF

      do_rtp = .FALSE.
      tmp_section => section_vals_get_subs_vals(dft_section, "REAL_TIME_PROPAGATION")
      CALL section_vals_get(tmp_section, explicit=is_present)
      IF (is_present) THEN
         CALL read_rtp_section(dft_control, tmp_section)
         do_rtp = .TRUE.
      END IF

      ! Read the input section
      tmp_section => section_vals_get_subs_vals(dft_section, "XAS")
      CALL section_vals_get(tmp_section, explicit=dft_control%do_xas_calculation)
      IF (dft_control%do_xas_calculation) THEN
         ! Override with section parameter
         CALL section_vals_val_get(tmp_section, "_SECTION_PARAMETERS_", &
                                   l_val=dft_control%do_xas_calculation)
      END IF

      tmp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")
      CALL section_vals_get(tmp_section, explicit=dft_control%do_xas_tdp_calculation)
      IF (dft_control%do_xas_tdp_calculation) THEN
         ! Override with section parameter
         CALL section_vals_val_get(tmp_section, "_SECTION_PARAMETERS_", &
                                   l_val=dft_control%do_xas_tdp_calculation)
      END IF

      ! Read the finite field input section
      dft_control%apply_efield = .FALSE.
      dft_control%apply_efield_field = .FALSE. !this is for RTP
      tmp_section => section_vals_get_subs_vals(dft_section, "EFIELD")
      CALL section_vals_get(tmp_section, n_repetition=nrep, explicit=is_present)
      IF (is_present) THEN
         ALLOCATE (dft_control%efield_fields(nrep))
         CALL read_efield_sections(dft_control, tmp_section)
         IF (do_rtp) THEN
            dft_control%apply_efield_field = .TRUE.
         ELSE
            dft_control%apply_efield = .TRUE.
            CPASSERT(nrep == 1)
         END IF
      END IF

      ! Read the finite field input section for periodic fields
      tmp_section => section_vals_get_subs_vals(dft_section, "PERIODIC_EFIELD")
      CALL section_vals_get(tmp_section, explicit=dft_control%apply_period_efield)
      IF (dft_control%apply_period_efield) THEN
         ALLOCATE (dft_control%period_efield)
         CALL section_vals_val_get(tmp_section, "POLARISATION", r_vals=pol)
         dft_control%period_efield%polarisation(1:3) = pol(1:3)
         CALL section_vals_val_get(tmp_section, "D_FILTER", r_vals=pol)
         dft_control%period_efield%d_filter(1:3) = pol(1:3)
         CALL section_vals_val_get(tmp_section, "INTENSITY", &
                                   r_val=dft_control%period_efield%strength)
         dft_control%period_efield%displacement_field = .FALSE.
         CALL section_vals_val_get(tmp_section, "DISPLACEMENT_FIELD", &
                                   l_val=dft_control%period_efield%displacement_field)
         ! periodic fields don't work with RTP
         CPASSERT(.NOT. do_rtp)
         IF (dft_control%period_efield%displacement_field) THEN
            CALL cite_reference(Stengel2009)
         ELSE
            CALL cite_reference(Souza2002)
            CALL cite_reference(Umari2002)
         END IF
      END IF

      ! Read the external potential input section
      tmp_section => section_vals_get_subs_vals(dft_section, "EXTERNAL_POTENTIAL")
      CALL section_vals_get(tmp_section, explicit=dft_control%apply_external_potential)

      ! Read the SCCS input section if present
      sccs_section => section_vals_get_subs_vals(dft_section, "SCCS")
      CALL section_vals_get(sccs_section, explicit=is_present)
      IF (is_present) THEN
         ! Check section parameter if SCCS is activated
         CALL section_vals_val_get(sccs_section, "_SECTION_PARAMETERS_", &
                                   l_val=dft_control%do_sccs)
         IF (dft_control%do_sccs) THEN
            CALL sccs_control_create(dft_control%sccs_control)
            CALL section_vals_val_get(sccs_section, "ALPHA", &
                                      r_val=dft_control%sccs_control%alpha_solvent)
            CALL section_vals_val_get(sccs_section, "BETA", &
                                      r_val=dft_control%sccs_control%beta_solvent)
            CALL section_vals_val_get(sccs_section, "DELTA_RHO", &
                                      r_val=dft_control%sccs_control%delta_rho)
            CALL section_vals_val_get(sccs_section, "DERIVATIVE_METHOD", &
                                      i_val=dft_control%sccs_control%derivative_method)
            CALL section_vals_val_get(sccs_section, "METHOD", &
                                      i_val=dft_control%sccs_control%method_id)
            CALL section_vals_val_get(sccs_section, "DIELECTRIC_CONSTANT", &
                                      r_val=dft_control%sccs_control%epsilon_solvent)
            CALL section_vals_val_get(sccs_section, "EPS_SCCS", &
                                      r_val=dft_control%sccs_control%eps_sccs)
            CALL section_vals_val_get(sccs_section, "EPS_SCF", &
                                      r_val=dft_control%sccs_control%eps_scf)
            CALL section_vals_val_get(sccs_section, "GAMMA", &
                                      r_val=dft_control%sccs_control%gamma_solvent)
            CALL section_vals_val_get(sccs_section, "MAX_ITER", &
                                      i_val=dft_control%sccs_control%max_iter)
            CALL section_vals_val_get(sccs_section, "MIXING", &
                                      r_val=dft_control%sccs_control%mixing)
            SELECT CASE (dft_control%sccs_control%method_id)
            CASE (sccs_andreussi)
               tmp_section => section_vals_get_subs_vals(sccs_section, "ANDREUSSI")
               CALL section_vals_val_get(tmp_section, "RHO_MAX", &
                                         r_val=dft_control%sccs_control%rho_max)
               CALL section_vals_val_get(tmp_section, "RHO_MIN", &
                                         r_val=dft_control%sccs_control%rho_min)
               IF (dft_control%sccs_control%rho_max < dft_control%sccs_control%rho_min) THEN
                  CALL cp_abort(__LOCATION__, &
                                "The SCCS parameter RHO_MAX is smaller than RHO_MIN. "// &
                                "Please, check your input!")
               END IF
               CALL cite_reference(Andreussi2012)
            CASE (sccs_fattebert_gygi)
               tmp_section => section_vals_get_subs_vals(sccs_section, "FATTEBERT-GYGI")
               CALL section_vals_val_get(tmp_section, "BETA", &
                                         r_val=dft_control%sccs_control%beta)
               IF (dft_control%sccs_control%beta < 0.5_dp) THEN
                  CALL cp_abort(__LOCATION__, &
                                "A value smaller than 0.5 for the SCCS parameter beta "// &
                                "causes numerical problems. Please, check your input!")
               END IF
               CALL section_vals_val_get(tmp_section, "RHO_ZERO", &
                                         r_val=dft_control%sccs_control%rho_zero)
               CALL cite_reference(Fattebert2002)
            CASE DEFAULT
               CPABORT("Invalid SCCS model specified. Please, check your input!")
            END SELECT
            CALL cite_reference(Yin2017)
         END IF
      END IF

      ! ZMP added input sections
      ! Read the external density input section
      tmp_section => section_vals_get_subs_vals(dft_section, "EXTERNAL_DENSITY")
      CALL section_vals_get(tmp_section, explicit=dft_control%apply_external_density)

      ! Read the external vxc input section
      tmp_section => section_vals_get_subs_vals(dft_section, "EXTERNAL_VXC")
      CALL section_vals_get(tmp_section, explicit=dft_control%apply_external_vxc)

   END SUBROUTINE read_dft_control

! **************************************************************************************************
!> \brief ...
!> \param qs_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE read_mgrid_section(qs_control, dft_section)

      TYPE(qs_control_type), INTENT(INOUT)               :: qs_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_mgrid_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igrid_level, ngrid_level
      LOGICAL                                            :: explicit, multigrid_set
      REAL(dp)                                           :: cutoff
      REAL(dp), DIMENSION(:), POINTER                    :: cutofflist
      TYPE(section_vals_type), POINTER                   :: mgrid_section

      CALL timeset(routineN, handle)

      NULLIFY (mgrid_section, cutofflist)
      mgrid_section => section_vals_get_subs_vals(dft_section, "MGRID")

      CALL section_vals_val_get(mgrid_section, "NGRIDS", i_val=ngrid_level)
      CALL section_vals_val_get(mgrid_section, "MULTIGRID_SET", l_val=multigrid_set)
      CALL section_vals_val_get(mgrid_section, "CUTOFF", r_val=cutoff)
      CALL section_vals_val_get(mgrid_section, "PROGRESSION_FACTOR", r_val=qs_control%progression_factor)
      CALL section_vals_val_get(mgrid_section, "COMMENSURATE", l_val=qs_control%commensurate_mgrids)
      CALL section_vals_val_get(mgrid_section, "REALSPACE", l_val=qs_control%realspace_mgrids)
      CALL section_vals_val_get(mgrid_section, "REL_CUTOFF", r_val=qs_control%relative_cutoff)
      CALL section_vals_val_get(mgrid_section, "SKIP_LOAD_BALANCE_DISTRIBUTED", &
                                l_val=qs_control%skip_load_balance_distributed)

      ! For SE and DFTB possibly override with new defaults
      IF (qs_control%semi_empirical .OR. qs_control%dftb .OR. qs_control%xtb) THEN
         ngrid_level = 1
         multigrid_set = .FALSE.
         ! Override default cutoff value unless user specified an explicit argument..
         CALL section_vals_val_get(mgrid_section, "CUTOFF", explicit=explicit, r_val=cutoff)
         IF (.NOT. explicit) cutoff = 1.0_dp
      END IF

      ALLOCATE (qs_control%e_cutoff(ngrid_level))
      qs_control%cutoff = cutoff

      IF (multigrid_set) THEN
         ! Read the values from input
         IF (qs_control%commensurate_mgrids) THEN
            CPABORT("Do not specify cutoffs for the commensurate grids (NYI)")
         END IF

         CALL section_vals_val_get(mgrid_section, "MULTIGRID_CUTOFF", r_vals=cutofflist)
         IF (ASSOCIATED(cutofflist)) THEN
            IF (SIZE(cutofflist, 1) /= ngrid_level) THEN
               CPABORT("Number of multi-grids requested and number of cutoff values do not match")
            END IF
            DO igrid_level = 1, ngrid_level
               qs_control%e_cutoff(igrid_level) = cutofflist(igrid_level)
            END DO
         END IF
         ! set cutoff to smallest value in multgrid available with >= cutoff
         DO igrid_level = ngrid_level, 1, -1
            IF (qs_control%cutoff <= qs_control%e_cutoff(igrid_level)) THEN
               qs_control%cutoff = qs_control%e_cutoff(igrid_level)
               EXIT
            END IF
            ! set largest grid value to cutoff
            IF (igrid_level == 1) THEN
               qs_control%cutoff = qs_control%e_cutoff(1)
            END IF
         END DO
      ELSE
         IF (qs_control%commensurate_mgrids) qs_control%progression_factor = 4.0_dp
         qs_control%e_cutoff(1) = qs_control%cutoff
         DO igrid_level = 2, ngrid_level
            qs_control%e_cutoff(igrid_level) = qs_control%e_cutoff(igrid_level - 1)/ &
                                               qs_control%progression_factor
         END DO
      END IF
      ! check that multigrids are ordered
      DO igrid_level = 2, ngrid_level
         IF (qs_control%e_cutoff(igrid_level) > qs_control%e_cutoff(igrid_level - 1)) THEN
            CPABORT("The cutoff values for the multi-grids are not ordered from large to small")
         ELSE IF (qs_control%e_cutoff(igrid_level) == qs_control%e_cutoff(igrid_level - 1)) THEN
            CPABORT("The same cutoff value was specified for two multi-grids")
         END IF
      END DO
      CALL timestop(handle)
   END SUBROUTINE read_mgrid_section

! **************************************************************************************************
!> \brief ...
!> \param qs_control ...
!> \param qs_section ...
! **************************************************************************************************
   SUBROUTINE read_qs_section(qs_control, qs_section)

      TYPE(qs_control_type), INTENT(INOUT)               :: qs_control
      TYPE(section_vals_type), POINTER                   :: qs_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_qs_section', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: clist
      INTEGER                                            :: handle, itmp, j, jj, k, n_rep, n_var, &
                                                            ngauss, nrep
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: explicit, was_present
      REAL(dp)                                           :: tmp, tmpsqrt, value
      REAL(dp), POINTER                                  :: scal(:)
      TYPE(section_vals_type), POINTER :: cdft_control_section, ddapc_restraint_section, &
         dftb_parameter, dftb_section, lri_optbas_section, mull_section, s2_restraint_section, &
         se_section, xtb_parameter, xtb_section

      CALL timeset(routineN, handle)

      was_present = .FALSE.
      NULLIFY (mull_section, ddapc_restraint_section, s2_restraint_section, &
               se_section, dftb_section, xtb_section, dftb_parameter, xtb_parameter, lri_optbas_section, &
               cdft_control_section)

      mull_section => section_vals_get_subs_vals(qs_section, "MULLIKEN_RESTRAINT")
      ddapc_restraint_section => section_vals_get_subs_vals(qs_section, "DDAPC_RESTRAINT")
      s2_restraint_section => section_vals_get_subs_vals(qs_section, "S2_RESTRAINT")
      se_section => section_vals_get_subs_vals(qs_section, "SE")
      dftb_section => section_vals_get_subs_vals(qs_section, "DFTB")
      xtb_section => section_vals_get_subs_vals(qs_section, "xTB")
      dftb_parameter => section_vals_get_subs_vals(dftb_section, "PARAMETER")
      xtb_parameter => section_vals_get_subs_vals(xtb_section, "PARAMETER")
      lri_optbas_section => section_vals_get_subs_vals(qs_section, "OPTIMIZE_LRI_BASIS")
      cdft_control_section => section_vals_get_subs_vals(qs_section, "CDFT")

      ! Setup all defaults values and overwrite input parameters
      ! EPS_DEFAULT should set the target accuracy in the total energy (~per electron) or a closely related value
      CALL section_vals_val_get(qs_section, "EPS_DEFAULT", r_val=value)
      tmpsqrt = SQRT(value) ! a trick to work around a NAG 5.1 optimizer bug

      ! random choice ?
      qs_control%eps_core_charge = value/100.0_dp
      ! correct if all Gaussians would have the same radius (overlap will be smaller than eps_pgf_orb**2).
      ! Can be significantly in error if not... requires fully new screening/pairlist procedures
      qs_control%eps_pgf_orb = tmpsqrt
      qs_control%eps_kg_orb = qs_control%eps_pgf_orb
      ! consistent since also a kind of overlap
      qs_control%eps_ppnl = qs_control%eps_pgf_orb/100.0_dp
      ! accuracy is basically set by the overlap, this sets an empirical shift
      qs_control%eps_ppl = 1.0E-2_dp
      !
      qs_control%gapw_control%eps_cpc = value
      ! expexted error in the density
      qs_control%eps_rho_gspace = value
      qs_control%eps_rho_rspace = value
      ! error in the gradient, can be the sqrt of the error in the energy, ignored if map_consistent
      qs_control%eps_gvg_rspace = tmpsqrt
      !
      CALL section_vals_val_get(qs_section, "EPS_CORE_CHARGE", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_CORE_CHARGE", r_val=qs_control%eps_core_charge)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_GVG_RSPACE", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_GVG_RSPACE", r_val=qs_control%eps_gvg_rspace)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", r_val=qs_control%eps_pgf_orb)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_KG_ORB", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_KG_ORB", r_val=tmp)
         qs_control%eps_kg_orb = SQRT(tmp)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_PPL", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_PPL", r_val=qs_control%eps_ppl)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_PPNL", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_PPNL", r_val=qs_control%eps_ppnl)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_RHO", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_RHO", r_val=qs_control%eps_rho_gspace)
         qs_control%eps_rho_rspace = qs_control%eps_rho_gspace
      END IF
      CALL section_vals_val_get(qs_section, "EPS_RHO_RSPACE", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_RHO_RSPACE", r_val=qs_control%eps_rho_rspace)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_RHO_GSPACE", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_RHO_GSPACE", r_val=qs_control%eps_rho_gspace)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_FILTER_MATRIX", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_FILTER_MATRIX", r_val=qs_control%eps_filter_matrix)
      END IF
      CALL section_vals_val_get(qs_section, "EPS_CPC", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_CPC", r_val=qs_control%gapw_control%eps_cpc)
      END IF

      CALL section_vals_val_get(qs_section, "EPSFIT", r_val=qs_control%gapw_control%eps_fit)
      CALL section_vals_val_get(qs_section, "EPSISO", r_val=qs_control%gapw_control%eps_iso)
      CALL section_vals_val_get(qs_section, "EPSSVD", r_val=qs_control%gapw_control%eps_svd)
      CALL section_vals_val_get(qs_section, "EPSRHO0", r_val=qs_control%gapw_control%eps_Vrho0)
      CALL section_vals_val_get(qs_section, "ALPHA0_HARD", r_val=qs_control%gapw_control%alpha0_hard)
      qs_control%gapw_control%lrho1_eq_lrho0 = .FALSE.
      qs_control%gapw_control%alpha0_hard_from_input = .FALSE.
      IF (qs_control%gapw_control%alpha0_hard /= 0.0_dp) qs_control%gapw_control%alpha0_hard_from_input = .TRUE.
      CALL section_vals_val_get(qs_section, "FORCE_PAW", l_val=qs_control%gapw_control%force_paw)
      CALL section_vals_val_get(qs_section, "MAX_RAD_LOCAL", r_val=qs_control%gapw_control%max_rad_local)

      CALL section_vals_val_get(qs_section, "LS_SCF", l_val=qs_control%do_ls_scf)
      CALL section_vals_val_get(qs_section, "ALMO_SCF", l_val=qs_control%do_almo_scf)
      CALL section_vals_val_get(qs_section, "KG_METHOD", l_val=qs_control%do_kg)

      ! Logicals
      CALL section_vals_val_get(qs_section, "REF_EMBED_SUBSYS", l_val=qs_control%ref_embed_subsys)
      CALL section_vals_val_get(qs_section, "CLUSTER_EMBED_SUBSYS", l_val=qs_control%cluster_embed_subsys)
      CALL section_vals_val_get(qs_section, "HIGH_LEVEL_EMBED_SUBSYS", l_val=qs_control%high_level_embed_subsys)
      CALL section_vals_val_get(qs_section, "DFET_EMBEDDED", l_val=qs_control%dfet_embedded)
      CALL section_vals_val_get(qs_section, "DMFET_EMBEDDED", l_val=qs_control%dmfet_embedded)

      ! Integers gapw
      CALL section_vals_val_get(qs_section, "LMAXN1", i_val=qs_control%gapw_control%lmax_sphere)
      CALL section_vals_val_get(qs_section, "LMAXN0", i_val=qs_control%gapw_control%lmax_rho0)
      CALL section_vals_val_get(qs_section, "LADDN0", i_val=qs_control%gapw_control%ladd_rho0)
      CALL section_vals_val_get(qs_section, "QUADRATURE", i_val=qs_control%gapw_control%quadrature)

      ! Integers grids
      CALL section_vals_val_get(qs_section, "PW_GRID", i_val=itmp)
      SELECT CASE (itmp)
      CASE (do_pwgrid_spherical)
         qs_control%pw_grid_opt%spherical = .TRUE.
         qs_control%pw_grid_opt%fullspace = .FALSE.
      CASE (do_pwgrid_ns_fullspace)
         qs_control%pw_grid_opt%spherical = .FALSE.
         qs_control%pw_grid_opt%fullspace = .TRUE.
      CASE (do_pwgrid_ns_halfspace)
         qs_control%pw_grid_opt%spherical = .FALSE.
         qs_control%pw_grid_opt%fullspace = .FALSE.
      END SELECT

      !   Method for PPL calculation
      CALL section_vals_val_get(qs_section, "CORE_PPL", i_val=itmp)
      qs_control%do_ppl_method = itmp

      CALL section_vals_val_get(qs_section, "PW_GRID_LAYOUT", i_vals=tmplist)
      qs_control%pw_grid_opt%distribution_layout = tmplist
      CALL section_vals_val_get(qs_section, "PW_GRID_BLOCKED", i_val=qs_control%pw_grid_opt%blocked)

      !Integers extrapolation
      CALL section_vals_val_get(qs_section, "EXTRAPOLATION", i_val=qs_control%wf_interpolation_method_nr)
      CALL section_vals_val_get(qs_section, "EXTRAPOLATION_ORDER", i_val=qs_control%wf_extrapolation_order)

      !Method
      CALL section_vals_val_get(qs_section, "METHOD", i_val=qs_control%method_id)
      qs_control%gapw = .FALSE.
      qs_control%gapw_xc = .FALSE.
      qs_control%gpw = .FALSE.
      qs_control%pao = .FALSE.
      qs_control%dftb = .FALSE.
      qs_control%xtb = .FALSE.
      qs_control%semi_empirical = .FALSE.
      qs_control%ofgpw = .FALSE.
      qs_control%lrigpw = .FALSE.
      qs_control%rigpw = .FALSE.
      SELECT CASE (qs_control%method_id)
      CASE (do_method_gapw)
         CALL cite_reference(Lippert1999)
         CALL cite_reference(Krack2000)
         qs_control%gapw = .TRUE.
      CASE (do_method_gapw_xc)
         qs_control%gapw_xc = .TRUE.
      CASE (do_method_gpw)
         CALL cite_reference(Lippert1997)
         CALL cite_reference(VandeVondele2005a)
         qs_control%gpw = .TRUE.
      CASE (do_method_ofgpw)
         qs_control%ofgpw = .TRUE.
      CASE (do_method_lrigpw)
         qs_control%lrigpw = .TRUE.
      CASE (do_method_rigpw)
         qs_control%rigpw = .TRUE.
      CASE (do_method_dftb)
         qs_control%dftb = .TRUE.
         CALL cite_reference(Porezag1995)
         CALL cite_reference(Seifert1996)
      CASE (do_method_xtb)
         qs_control%xtb = .TRUE.
         CALL cite_reference(Grimme2017)
      CASE (do_method_mndo)
         CALL cite_reference(Dewar1977)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_am1)
         CALL cite_reference(Dewar1985)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_pm3)
         CALL cite_reference(Stewart1989)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_pnnl)
         CALL cite_reference(Schenter2008)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_pm6)
         CALL cite_reference(Stewart2007)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_pm6fm)
         CALL cite_reference(VanVoorhis2015)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_pdg)
         CALL cite_reference(Repasky2002)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_rm1)
         CALL cite_reference(Rocha2006)
         qs_control%semi_empirical = .TRUE.
      CASE (do_method_mndod)
         CALL cite_reference(Dewar1977)
         CALL cite_reference(Thiel1992)
         qs_control%semi_empirical = .TRUE.
      END SELECT

      CALL section_vals_get(mull_section, explicit=qs_control%mulliken_restraint)

      IF (qs_control%mulliken_restraint) THEN
         CALL section_vals_val_get(mull_section, "STRENGTH", r_val=qs_control%mulliken_restraint_control%strength)
         CALL section_vals_val_get(mull_section, "TARGET", r_val=qs_control%mulliken_restraint_control%target)
         CALL section_vals_val_get(mull_section, "ATOMS", n_rep_val=n_rep)
         jj = 0
         DO k = 1, n_rep
            CALL section_vals_val_get(mull_section, "ATOMS", i_rep_val=k, i_vals=tmplist)
            jj = jj + SIZE(tmplist)
         END DO
         qs_control%mulliken_restraint_control%natoms = jj
         IF (qs_control%mulliken_restraint_control%natoms < 1) &
            CPABORT("Need at least 1 atom to use mulliken constraints")
         ALLOCATE (qs_control%mulliken_restraint_control%atoms(qs_control%mulliken_restraint_control%natoms))
         jj = 0
         DO k = 1, n_rep
            CALL section_vals_val_get(mull_section, "ATOMS", i_rep_val=k, i_vals=tmplist)
            DO j = 1, SIZE(tmplist)
               jj = jj + 1
               qs_control%mulliken_restraint_control%atoms(jj) = tmplist(j)
            END DO
         END DO
      ENDIF
      CALL section_vals_get(ddapc_restraint_section, n_repetition=nrep, explicit=qs_control%ddapc_restraint)
      IF (qs_control%ddapc_restraint) THEN
         ALLOCATE (qs_control%ddapc_restraint_control(nrep))
         CALL read_ddapc_section(qs_control, qs_section=qs_section)
         qs_control%ddapc_restraint_is_spin = .FALSE.
         qs_control%ddapc_explicit_potential = .FALSE.
      ENDIF

      CALL section_vals_get(s2_restraint_section, explicit=qs_control%s2_restraint)
      IF (qs_control%s2_restraint) THEN
         CALL section_vals_val_get(s2_restraint_section, "STRENGTH", &
                                   r_val=qs_control%s2_restraint_control%strength)
         CALL section_vals_val_get(s2_restraint_section, "TARGET", &
                                   r_val=qs_control%s2_restraint_control%target)
         CALL section_vals_val_get(s2_restraint_section, "FUNCTIONAL_FORM", &
                                   i_val=qs_control%s2_restraint_control%functional_form)
      ENDIF

      CALL section_vals_get(cdft_control_section, explicit=qs_control%cdft)
      IF (qs_control%cdft) THEN
         CALL read_cdft_control_section(qs_control, cdft_control_section)
      ENDIF

      ! Semi-empirical code
      IF (qs_control%semi_empirical) THEN
         CALL section_vals_val_get(se_section, "ORTHOGONAL_BASIS", &
                                   l_val=qs_control%se_control%orthogonal_basis)
         CALL section_vals_val_get(se_section, "DELTA", &
                                   r_val=qs_control%se_control%delta)
         CALL section_vals_val_get(se_section, "ANALYTICAL_GRADIENTS", &
                                   l_val=qs_control%se_control%analytical_gradients)
         CALL section_vals_val_get(se_section, "FORCE_KDSO-D_EXCHANGE", &
                                   l_val=qs_control%se_control%force_kdsod_EX)
         ! Integral Screening
         CALL section_vals_val_get(se_section, "INTEGRAL_SCREENING", &
                                   i_val=qs_control%se_control%integral_screening)
         IF (qs_control%method_id == do_method_pnnl) THEN
            IF (qs_control%se_control%integral_screening /= do_se_IS_slater) &
               CALL cp_warn(__LOCATION__, &
                            "PNNL semi-empirical parameterization supports only the Slater type "// &
                            "integral scheme. Revert to Slater and continue the calculation.")
            qs_control%se_control%integral_screening = do_se_IS_slater
         END IF
         ! Global Arrays variable
         CALL section_vals_val_get(se_section, "GA%NCELLS", &
                                   i_val=qs_control%se_control%ga_ncells)
         ! Long-Range correction
         CALL section_vals_val_get(se_section, "LR_CORRECTION%CUTOFF", &
                                   r_val=qs_control%se_control%cutoff_lrc)
         qs_control%se_control%taper_lrc = qs_control%se_control%cutoff_lrc
         CALL section_vals_val_get(se_section, "LR_CORRECTION%RC_TAPER", &
                                   explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(se_section, "LR_CORRECTION%RC_TAPER", &
                                      r_val=qs_control%se_control%taper_lrc)
         END IF
         CALL section_vals_val_get(se_section, "LR_CORRECTION%RC_RANGE", &
                                   r_val=qs_control%se_control%range_lrc)
         ! Coulomb
         CALL section_vals_val_get(se_section, "COULOMB%CUTOFF", &
                                   r_val=qs_control%se_control%cutoff_cou)
         qs_control%se_control%taper_cou = qs_control%se_control%cutoff_cou
         CALL section_vals_val_get(se_section, "COULOMB%RC_TAPER", &
                                   explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(se_section, "COULOMB%RC_TAPER", &
                                      r_val=qs_control%se_control%taper_cou)
         END IF
         CALL section_vals_val_get(se_section, "COULOMB%RC_RANGE", &
                                   r_val=qs_control%se_control%range_cou)
         ! Exchange
         CALL section_vals_val_get(se_section, "EXCHANGE%CUTOFF", &
                                   r_val=qs_control%se_control%cutoff_exc)
         qs_control%se_control%taper_exc = qs_control%se_control%cutoff_exc
         CALL section_vals_val_get(se_section, "EXCHANGE%RC_TAPER", &
                                   explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(se_section, "EXCHANGE%RC_TAPER", &
                                      r_val=qs_control%se_control%taper_exc)
         END IF
         CALL section_vals_val_get(se_section, "EXCHANGE%RC_RANGE", &
                                   r_val=qs_control%se_control%range_exc)
         ! Screening (only if the integral scheme is of dumped type)
         IF (qs_control%se_control%integral_screening == do_se_IS_kdso_d) THEN
            CALL section_vals_val_get(se_section, "SCREENING%RC_TAPER", &
                                      r_val=qs_control%se_control%taper_scr)
            CALL section_vals_val_get(se_section, "SCREENING%RC_RANGE", &
                                      r_val=qs_control%se_control%range_scr)
         END IF
         ! Periodic Type Calculation
         CALL section_vals_val_get(se_section, "PERIODIC", &
                                   i_val=qs_control%se_control%periodic_type)
         SELECT CASE (qs_control%se_control%periodic_type)
         CASE (do_se_lr_none)
            qs_control%se_control%do_ewald = .FALSE.
            qs_control%se_control%do_ewald_r3 = .FALSE.
            qs_control%se_control%do_ewald_gks = .FALSE.
         CASE (do_se_lr_ewald)
            qs_control%se_control%do_ewald = .TRUE.
            qs_control%se_control%do_ewald_r3 = .FALSE.
            qs_control%se_control%do_ewald_gks = .FALSE.
         CASE (do_se_lr_ewald_gks)
            qs_control%se_control%do_ewald = .FALSE.
            qs_control%se_control%do_ewald_r3 = .FALSE.
            qs_control%se_control%do_ewald_gks = .TRUE.
            IF (qs_control%method_id /= do_method_pnnl) &
               CALL cp_abort(__LOCATION__, &
                             "A periodic semi-empirical calculation was requested with a long-range  "// &
                             "summation on the single integral evaluation. This scheme is supported  "// &
                             "only by the PNNL parameterization.")
         CASE (do_se_lr_ewald_r3)
            qs_control%se_control%do_ewald = .TRUE.
            qs_control%se_control%do_ewald_r3 = .TRUE.
            qs_control%se_control%do_ewald_gks = .FALSE.
            IF (qs_control%se_control%integral_screening /= do_se_IS_kdso) &
               CALL cp_abort(__LOCATION__, &
                             "A periodic semi-empirical calculation was requested with a long-range  "// &
                             "summation for the slowly convergent part 1/R^3, which is not congruent "// &
                             "with the integral screening chosen. The only integral screening supported "// &
                             "by this periodic type calculation is the standard Klopman-Dewar-Sabelli-Ohno.")
         END SELECT

         ! dispersion pair potentials
         CALL section_vals_val_get(se_section, "DISPERSION", &
                                   l_val=qs_control%se_control%dispersion)
         CALL section_vals_val_get(se_section, "DISPERSION_RADIUS", &
                                   r_val=qs_control%se_control%rcdisp)
         CALL section_vals_val_get(se_section, "COORDINATION_CUTOFF", &
                                   r_val=qs_control%se_control%epscn)
         CALL section_vals_val_get(se_section, "D3_SCALING", r_vals=scal)
         qs_control%se_control%sd3(1) = scal(1)
         qs_control%se_control%sd3(2) = scal(2)
         qs_control%se_control%sd3(3) = scal(3)
         CALL section_vals_val_get(se_section, "DISPERSION_PARAMETER_FILE", &
                                   c_val=qs_control%se_control%dispersion_parameter_file)

         ! Stop the execution for non-implemented features
         IF (qs_control%se_control%periodic_type == do_se_lr_ewald_r3) THEN
            CPABORT("EWALD_R3 not implemented yet!")
         END IF

         IF (qs_control%method_id == do_method_mndo .OR. &
             qs_control%method_id == do_method_am1 .OR. &
             qs_control%method_id == do_method_mndod .OR. &
             qs_control%method_id == do_method_pdg .OR. &
             qs_control%method_id == do_method_pm3 .OR. &
             qs_control%method_id == do_method_pm6 .OR. &
             qs_control%method_id == do_method_pm6fm .OR. &
             qs_control%method_id == do_method_pnnl .OR. &
             qs_control%method_id == do_method_rm1) THEN
            qs_control%se_control%orthogonal_basis = .TRUE.
         END IF
      END IF

      ! DFTB code
      IF (qs_control%dftb) THEN
         CALL section_vals_val_get(dftb_section, "ORTHOGONAL_BASIS", &
                                   l_val=qs_control%dftb_control%orthogonal_basis)
         CALL section_vals_val_get(dftb_section, "SELF_CONSISTENT", &
                                   l_val=qs_control%dftb_control%self_consistent)
         CALL section_vals_val_get(dftb_section, "DISPERSION", &
                                   l_val=qs_control%dftb_control%dispersion)
         CALL section_vals_val_get(dftb_section, "DIAGONAL_DFTB3", &
                                   l_val=qs_control%dftb_control%dftb3_diagonal)
         CALL section_vals_val_get(dftb_section, "HB_SR_GAMMA", &
                                   l_val=qs_control%dftb_control%hb_sr_damp)
         CALL section_vals_val_get(dftb_section, "EPS_DISP", &
                                   r_val=qs_control%dftb_control%eps_disp)
         CALL section_vals_val_get(dftb_section, "DO_EWALD", &
                                   l_val=qs_control%dftb_control%do_ewald)
         CALL section_vals_val_get(dftb_parameter, "PARAM_FILE_PATH", &
                                   c_val=qs_control%dftb_control%sk_file_path)
         CALL section_vals_val_get(dftb_parameter, "PARAM_FILE_NAME", &
                                   c_val=qs_control%dftb_control%sk_file_list)
         CALL section_vals_val_get(dftb_parameter, "HB_SR_PARAM", &
                                   r_val=qs_control%dftb_control%hb_sr_para)
         CALL section_vals_val_get(dftb_parameter, "SK_FILE", n_rep_val=n_var)
         ALLOCATE (qs_control%dftb_control%sk_pair_list(3, n_var))
         DO k = 1, n_var
            CALL section_vals_val_get(dftb_parameter, "SK_FILE", i_rep_val=k, &
                                      c_vals=clist)
            qs_control%dftb_control%sk_pair_list(1:3, k) = clist(1:3)
         END DO
         ! Dispersion type
         CALL section_vals_val_get(dftb_parameter, "DISPERSION_TYPE", &
                                   i_val=qs_control%dftb_control%dispersion_type)
         CALL section_vals_val_get(dftb_parameter, "UFF_FORCE_FIELD", &
                                   c_val=qs_control%dftb_control%uff_force_field)
         ! D3 Dispersion
         CALL section_vals_val_get(dftb_parameter, "DISPERSION_RADIUS", &
                                   r_val=qs_control%dftb_control%rcdisp)
         CALL section_vals_val_get(dftb_parameter, "COORDINATION_CUTOFF", &
                                   r_val=qs_control%dftb_control%epscn)
         CALL section_vals_val_get(dftb_parameter, "D3_SCALING", r_vals=scal)
         qs_control%dftb_control%sd3(1) = scal(1)
         qs_control%dftb_control%sd3(2) = scal(2)
         qs_control%dftb_control%sd3(3) = scal(3)
         CALL section_vals_val_get(dftb_parameter, "DISPERSION_PARAMETER_FILE", &
                                   c_val=qs_control%dftb_control%dispersion_parameter_file)

         IF (qs_control%dftb_control%dispersion) CALL cite_reference(Zhechkov2005)
         IF (qs_control%dftb_control%self_consistent) CALL cite_reference(Elstner1998)
         IF (qs_control%dftb_control%hb_sr_damp) CALL cite_reference(Hu2007)
      END IF

      ! xTB code
      IF (qs_control%xtb) THEN
         CALL section_vals_val_get(xtb_section, "DO_EWALD", &
                                   l_val=qs_control%xtb_control%do_ewald)
         CALL section_vals_val_get(xtb_section, "STO_NG", i_val=ngauss)
         qs_control%xtb_control%sto_ng = ngauss
         CALL section_vals_val_get(xtb_section, "HYDROGEN_STO_NG", i_val=ngauss)
         qs_control%xtb_control%h_sto_ng = ngauss
         CALL section_vals_val_get(xtb_parameter, "PARAM_FILE_PATH", &
                                   c_val=qs_control%xtb_control%parameter_file_path)
         CALL section_vals_val_get(xtb_parameter, "PARAM_FILE_NAME", &
                                   c_val=qs_control%xtb_control%parameter_file_name)
         ! D3 Dispersion
         CALL section_vals_val_get(xtb_parameter, "DISPERSION_RADIUS", &
                                   r_val=qs_control%xtb_control%rcdisp)
         CALL section_vals_val_get(xtb_parameter, "COORDINATION_CUTOFF", &
                                   r_val=qs_control%xtb_control%epscn)
         CALL section_vals_val_get(xtb_parameter, "D3BJ_SCALING", r_vals=scal)
         qs_control%xtb_control%s6 = scal(1)
         qs_control%xtb_control%s8 = scal(2)
         CALL section_vals_val_get(xtb_parameter, "D3BJ_PARAM", r_vals=scal)
         qs_control%xtb_control%a1 = scal(1)
         qs_control%xtb_control%a2 = scal(2)
         CALL section_vals_val_get(xtb_parameter, "DISPERSION_PARAMETER_FILE", &
                                   c_val=qs_control%xtb_control%dispersion_parameter_file)
         ! global parameters
         CALL section_vals_val_get(xtb_parameter, "HUCKEL_CONSTANTS", r_vals=scal)
         qs_control%xtb_control%ks = scal(1)
         qs_control%xtb_control%kp = scal(2)
         qs_control%xtb_control%kd = scal(3)
         qs_control%xtb_control%ksp = scal(4)
         qs_control%xtb_control%k2sh = scal(5)
         CALL section_vals_val_get(xtb_parameter, "COULOMB_CONSTANTS", r_vals=scal)
         qs_control%xtb_control%kg = scal(1)
         qs_control%xtb_control%kf = scal(2)
         CALL section_vals_val_get(xtb_parameter, "CN_CONSTANTS", r_vals=scal)
         qs_control%xtb_control%kcns = scal(1)
         qs_control%xtb_control%kcnp = scal(2)
         qs_control%xtb_control%kcnd = scal(3)
         CALL section_vals_val_get(xtb_parameter, "EN_CONSTANT", r_vals=scal)
         qs_control%xtb_control%ken = scal(1)
         ! XB
         CALL section_vals_val_get(xtb_section, "USE_HALOGEN_CORRECTION", &
                                   l_val=qs_control%xtb_control%xb_interaction)
         CALL section_vals_val_get(xtb_parameter, "HALOGEN_BINDING", r_vals=scal)
         qs_control%xtb_control%kxr = scal(1)
         qs_control%xtb_control%kx2 = scal(2)
         ! XB_radius
         CALL section_vals_val_get(xtb_parameter, "XB_RADIUS", r_val=qs_control%xtb_control%xb_radius)
         ! Kab
         CALL section_vals_val_get(xtb_parameter, "KAB_PARAM", n_rep_val=n_rep)
         ! For debug purposes
         CALL section_vals_val_get(xtb_section, "COULOMB_INTERACTION", &
                                   l_val=qs_control%xtb_control%coulomb_interaction)
         CALL section_vals_val_get(xtb_section, "TB3_INTERACTION", &
                                   l_val=qs_control%xtb_control%tb3_interaction)
         ! Check for bad atomic charges
         CALL section_vals_val_get(xtb_section, "CHECK_ATOMIC_CHARGES", &
                                   l_val=qs_control%xtb_control%check_atomic_charges)

         qs_control%xtb_control%kab_nval = n_rep
         IF (n_rep > 0) THEN
            ALLOCATE (qs_control%xtb_control%kab_param(3, n_rep))
            ALLOCATE (qs_control%xtb_control%kab_types(2, n_rep))
            ALLOCATE (qs_control%xtb_control%kab_vals(n_rep))
            DO j = 1, n_rep
               CALL section_vals_val_get(xtb_parameter, "KAB_PARAM", i_rep_val=j, c_vals=clist)
               qs_control%xtb_control%kab_param(1, j) = clist(1)
               CALL get_ptable_info(clist(1), &
                                    ielement=qs_control%xtb_control%kab_types(1, j))
               qs_control%xtb_control%kab_param(2, j) = clist(2)
               CALL get_ptable_info(clist(2), &
                                    ielement=qs_control%xtb_control%kab_types(2, j))
               qs_control%xtb_control%kab_param(3, j) = clist(3)
               READ (clist(3), '(F10.0)') qs_control%xtb_control%kab_vals(j)
            END DO
         END IF
      END IF

      ! Optimize LRI basis set
      CALL section_vals_get(lri_optbas_section, explicit=qs_control%lri_optbas)

      CALL timestop(handle)
   END SUBROUTINE read_qs_section

! **************************************************************************************************
!> \brief ...
!> \param t_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE read_tddfpt_control(t_control, dft_section)
      TYPE(tddfpt_control_type)                          :: t_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_tddfpt_control', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: kenergy_den
      TYPE(section_vals_type), POINTER                   :: sic_section, t_section

      kenergy_den = .FALSE.
      NULLIFY (sic_section, t_section)
      t_section => section_vals_get_subs_vals(dft_section, "TDDFPT")

      CALL section_vals_val_get(t_section, "CONVERGENCE", r_val=t_control%tolerance)
      CALL section_vals_val_get(t_section, "NEV", i_val=t_control%n_ev)
      CALL section_vals_val_get(t_section, "MAX_KV", i_val=t_control%max_kv)
      CALL section_vals_val_get(t_section, "RESTARTS", i_val=t_control%n_restarts)
      CALL section_vals_val_get(t_section, "NREORTHO", i_val=t_control%n_reortho)
      CALL section_vals_val_get(t_section, "RES_ETYPE", i_val=t_control%res_etype)
      CALL section_vals_val_get(t_section, "DIAG_METHOD", i_val=t_control%diag_method)
      CALL section_vals_val_get(t_section, "KERNEL", l_val=t_control%do_kernel)
      CALL section_vals_val_get(t_section, "LSD_SINGLETS", l_val=t_control%lsd_singlets)
      CALL section_vals_val_get(t_section, "INVERT_S", l_val=t_control%invert_S)
      CALL section_vals_val_get(t_section, "PRECOND", l_val=t_control%precond)
      CALL section_vals_val_get(t_section, "OE_CORR", i_val=t_control%oe_corr)

      t_control%use_kinetic_energy_density = .FALSE.
      sic_section => section_vals_get_subs_vals(t_section, "SIC")
      CALL section_vals_val_get(sic_section, "SIC_METHOD", i_val=t_control%sic_method_id)
      CALL section_vals_val_get(sic_section, "ORBITAL_SET", i_val=t_control%sic_list_id)
      CALL section_vals_val_get(sic_section, "SIC_SCALING_A", r_val=t_control%sic_scaling_a)
      CALL section_vals_val_get(sic_section, "SIC_SCALING_B", r_val=t_control%sic_scaling_b)

   END SUBROUTINE read_tddfpt_control

! **************************************************************************************************
!> \brief Read TDDFPT-related input parameters.
!> \param t_control  TDDFPT control parameters
!> \param t_section  TDDFPT input section
!> \param qs_control Quickstep control parameters
! **************************************************************************************************
   SUBROUTINE read_tddfpt2_control(t_control, t_section, qs_control)
      TYPE(tddfpt2_control_type), POINTER                :: t_control
      TYPE(section_vals_type), POINTER                   :: t_section
      TYPE(qs_control_type), POINTER                     :: qs_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_tddfpt2_control', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      LOGICAL                                            :: explicit, multigrid_set
      TYPE(section_vals_type), POINTER                   :: dipole_section, mgrid_section, &
                                                            stda_section, xc_func, xc_section

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(t_section, "_SECTION_PARAMETERS_", l_val=t_control%enabled)

      CALL section_vals_val_get(t_section, "NSTATES", i_val=t_control%nstates)
      CALL section_vals_val_get(t_section, "MAX_ITER", i_val=t_control%niters)
      CALL section_vals_val_get(t_section, "MAX_KV", i_val=t_control%nkvs)
      CALL section_vals_val_get(t_section, "NLUMO", i_val=t_control%nlumo)
      CALL section_vals_val_get(t_section, "NPROC_STATE", i_val=t_control%nprocs)
      CALL section_vals_val_get(t_section, "KERNEL", i_val=t_control%kernel)
      CALL section_vals_val_get(t_section, "OE_CORR", i_val=t_control%oe_corr)
      CALL section_vals_val_get(t_section, "EV_SHIFT", r_val=t_control%ev_shift)
      CALL section_vals_val_get(t_section, "EOS_SHIFT", r_val=t_control%eos_shift)

      CALL section_vals_val_get(t_section, "CONVERGENCE", r_val=t_control%conv)
      CALL section_vals_val_get(t_section, "MIN_AMPLITUDE", r_val=t_control%min_excitation_amplitude)
      CALL section_vals_val_get(t_section, "ORTHOGONAL_EPS", r_val=t_control%orthogonal_eps)

      CALL section_vals_val_get(t_section, "RESTART", l_val=t_control%is_restart)
      CALL section_vals_val_get(t_section, "RKS_TRIPLETS", l_val=t_control%rks_triplets)

      IF (t_control%conv < 0) &
         t_control%conv = ABS(t_control%conv)

      ! DIPOLE_MOMENTS subsection
      dipole_section => section_vals_get_subs_vals(t_section, "DIPOLE_MOMENTS")
      CALL section_vals_val_get(dipole_section, "DIPOLE_FORM", i_val=t_control%dipole_form)
      CALL section_vals_val_get(dipole_section, "REFERENCE", i_val=t_control%dipole_reference)
      CALL section_vals_val_get(dipole_section, "REFERENCE_POINT", explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(dipole_section, "REFERENCE_POINT", r_vals=t_control%dipole_ref_point)
      ELSE
         NULLIFY (t_control%dipole_ref_point)
         IF (t_control%dipole_form == xas_dip_len .AND. t_control%dipole_reference == use_mom_ref_user) THEN
            CPABORT("User-defined reference point should be given explicitly")
         END IF
      END IF

      ! MGRID subsection
      mgrid_section => section_vals_get_subs_vals(t_section, "MGRID")
      CALL section_vals_get(mgrid_section, explicit=t_control%mgrid_is_explicit)

      IF (t_control%mgrid_is_explicit) THEN
         CALL section_vals_val_get(mgrid_section, "NGRIDS", i_val=t_control%mgrid_ngrids, explicit=explicit)
         IF (.NOT. explicit) t_control%mgrid_ngrids = SIZE(qs_control%e_cutoff)

         CALL section_vals_val_get(mgrid_section, "CUTOFF", r_val=t_control%mgrid_cutoff, explicit=explicit)
         IF (.NOT. explicit) t_control%mgrid_cutoff = qs_control%cutoff

         CALL section_vals_val_get(mgrid_section, "PROGRESSION_FACTOR", &
                                   r_val=t_control%mgrid_progression_factor, explicit=explicit)
         IF (explicit) THEN
            IF (t_control%mgrid_progression_factor <= 1.0_dp) &
               CALL cp_abort(__LOCATION__, &
                             "Progression factor should be greater then 1.0 to ensure multi-grid ordering")
         ELSE
            t_control%mgrid_progression_factor = qs_control%progression_factor
         END IF

         CALL section_vals_val_get(mgrid_section, "COMMENSURATE", l_val=t_control%mgrid_commensurate_mgrids, explicit=explicit)
         IF (.NOT. explicit) t_control%mgrid_commensurate_mgrids = qs_control%commensurate_mgrids
         IF (t_control%mgrid_commensurate_mgrids) THEN
            IF (explicit) THEN
               t_control%mgrid_progression_factor = 4.0_dp
            ELSE
               t_control%mgrid_progression_factor = qs_control%progression_factor
            END IF
         END IF

         CALL section_vals_val_get(mgrid_section, "REL_CUTOFF", r_val=t_control%mgrid_relative_cutoff, explicit=explicit)
         IF (.NOT. explicit) t_control%mgrid_relative_cutoff = qs_control%relative_cutoff

         CALL section_vals_val_get(mgrid_section, "MULTIGRID_SET", l_val=multigrid_set, explicit=explicit)
         IF (.NOT. explicit) multigrid_set = .FALSE.
         IF (multigrid_set) THEN
            CALL section_vals_val_get(mgrid_section, "MULTIGRID_CUTOFF", r_vals=t_control%mgrid_e_cutoff)
         ELSE
            NULLIFY (t_control%mgrid_e_cutoff)
         END IF

         CALL section_vals_val_get(mgrid_section, "REALSPACE", l_val=t_control%mgrid_realspace_mgrids, explicit=explicit)
         IF (.NOT. explicit) t_control%mgrid_realspace_mgrids = qs_control%realspace_mgrids

         CALL section_vals_val_get(mgrid_section, "SKIP_LOAD_BALANCE_DISTRIBUTED", &
                                   l_val=t_control%mgrid_skip_load_balance, explicit=explicit)
         IF (.NOT. explicit) t_control%mgrid_skip_load_balance = qs_control%skip_load_balance_distributed

         IF (ASSOCIATED(t_control%mgrid_e_cutoff)) THEN
            IF (SIZE(t_control%mgrid_e_cutoff) /= t_control%mgrid_ngrids) &
               CPABORT("Inconsistent values for number of multi-grids")

            ! sort multi-grids in descending order according to their cutoff values
            t_control%mgrid_e_cutoff = -t_control%mgrid_e_cutoff
            ALLOCATE (inds(t_control%mgrid_ngrids))
            CALL sort(t_control%mgrid_e_cutoff, t_control%mgrid_ngrids, inds)
            DEALLOCATE (inds)
            t_control%mgrid_e_cutoff = -t_control%mgrid_e_cutoff
         END IF
      END IF

      ! expand XC subsection (if given explicitly)
      xc_section => section_vals_get_subs_vals(t_section, "XC")
      xc_func => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_get(xc_func, explicit=explicit)
      IF (explicit) &
         CALL xc_functionals_expand(xc_func, xc_section)

      ! sTDA subsection
      stda_section => section_vals_get_subs_vals(t_section, "STDA")
      CALL section_vals_get(stda_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(stda_section, "HFX_FRACTION", r_val=t_control%stda_control%hfx_fraction)
         CALL section_vals_val_get(stda_section, "EPS_TD_FILTER", r_val=t_control%stda_control%eps_td_filter)
         CALL section_vals_val_get(stda_section, "DO_EWALD", l_val=t_control%stda_control%do_ewald)
      ELSE
         t_control%stda_control%hfx_fraction = 0.0_dp
         t_control%stda_control%eps_td_filter = 1.e-10_dp
         t_control%stda_control%do_ewald = .FALSE.
      END IF

      CALL timestop(handle)
   END SUBROUTINE read_tddfpt2_control

! **************************************************************************************************
!> \brief Write the DFT control parameters to the output unit.
!> \param dft_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE write_dft_control(dft_control, dft_section)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'write_dft_control', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=20)                                  :: tmpStr
      INTEGER                                            :: handle, output_unit
      REAL(kind=dp)                                      :: density_cut, density_smooth_cut_range, &
                                                            gradient_cut, tau_cut
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: section
      TYPE(section_vals_type), POINTER                   :: xc_section

      IF (dft_control%qs_control%semi_empirical) RETURN
      IF (dft_control%qs_control%dftb) RETURN
      IF (dft_control%qs_control%xtb) THEN
         CALL write_xtb_control(dft_control%qs_control%xtb_control, dft_section)
         RETURN
      END IF
      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      output_unit = cp_print_key_unit_nr(logger, dft_section, &
                                         "PRINT%DFT_CONTROL_PARAMETERS", extension=".Log")

      IF (output_unit > 0) THEN

         xc_section => section_vals_get_subs_vals(dft_section, "XC")

         IF (dft_control%uks) THEN
            WRITE (UNIT=output_unit, FMT="(/,T2,A,T78,A)") &
               "DFT| Spin unrestricted (spin-polarized) Kohn-Sham calculation", "UKS"
         ELSE IF (dft_control%roks) THEN
            WRITE (UNIT=output_unit, FMT="(/,T2,A,T77,A)") &
               "DFT| Spin restricted open Kohn-Sham calculation", "ROKS"
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,T2,A,T78,A)") &
               "DFT| Spin restricted Kohn-Sham (RKS) calculation", "RKS"
         END IF

         WRITE (UNIT=output_unit, FMT="(T2,A,T76,I5)") &
            "DFT| Multiplicity", dft_control%multiplicity
         WRITE (UNIT=output_unit, FMT="(T2,A,T76,I5)") &
            "DFT| Number of spin states", dft_control%nspins

         WRITE (UNIT=output_unit, FMT="(T2,A,T76,I5)") &
            "DFT| Charge", dft_control%charge

         IF (dft_control%sic_method_id .NE. sic_none) CALL cite_reference(VandeVondele2005b)
         SELECT CASE (dft_control%sic_method_id)
         CASE (sic_none)
            tmpstr = "NO"
         CASE (sic_mauri_spz)
            tmpstr = "SPZ/MAURI SIC"
         CASE (sic_mauri_us)
            tmpstr = "US/MAURI SIC"
         CASE (sic_ad)
            tmpstr = "AD SIC"
         CASE (sic_eo)
            tmpstr = "Explicit Orbital SIC"
         CASE DEFAULT
            ! fix throughout the cp2k for this option
            CPABORT("SIC option unknown")
         END SELECT

         WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
            "DFT| Self-interaction correction (SIC)", ADJUSTR(TRIM(tmpstr))

         IF (dft_control%sic_method_id /= sic_none) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T66,ES15.6)") &
               "DFT| SIC scaling parameter a", dft_control%sic_scaling_a, &
               "DFT| SIC scaling parameter b", dft_control%sic_scaling_b
         END IF

         IF (dft_control%sic_method_id == sic_eo) THEN
            IF (dft_control%sic_list_id == sic_list_all) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,T66,A)") &
                  "DFT| SIC orbitals", "ALL"
            ENDIF
            IF (dft_control%sic_list_id == sic_list_unpaired) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,T66,A)") &
                  "DFT| SIC orbitals", "UNPAIRED"
            ENDIF
         END IF

         CALL section_vals_val_get(xc_section, "density_cutoff", r_val=density_cut)
         CALL section_vals_val_get(xc_section, "gradient_cutoff", r_val=gradient_cut)
         CALL section_vals_val_get(xc_section, "tau_cutoff", r_val=tau_cut)
         CALL section_vals_val_get(xc_section, "density_smooth_cutoff_range", r_val=density_smooth_cut_range)

         WRITE (UNIT=output_unit, FMT="(T2,A,T66,ES15.6)") &
            "DFT| Cutoffs: density ", density_cut, &
            "DFT|          gradient", gradient_cut, &
            "DFT|          tau     ", tau_cut, &
            "DFT|          cutoff_smoothing_range", density_smooth_cut_range
         CALL section_vals_val_get(xc_section, "XC_GRID%XC_SMOOTH_RHO", &
                                   c_val=tmpStr)
         WRITE (output_unit, '( A, T61, A )') &
            " DFT| XC density smoothing ", ADJUSTR(tmpStr)
         CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", &
                                   c_val=tmpStr)
         WRITE (output_unit, '( A, T61, A )') &
            " DFT| XC derivatives ", ADJUSTR(tmpStr)
         IF (dft_control%dft_plus_u) THEN
            NULLIFY (enum, keyword, section)
            CALL create_dft_section(section)
            keyword => section_get_keyword(section, "PLUS_U_METHOD")
            CALL keyword_get(keyword, enum=enum)
            WRITE (UNIT=output_unit, FMT="(/,T2,A,T41,A40)") &
               "DFT+U| Method", ADJUSTR(TRIM(enum_i2c(enum, dft_control%plus_u_method_id)))
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "DFT+U| Check atomic kind information for details"
            CALL section_release(section)
         END IF

         CALL xc_write(output_unit, xc_section, dft_control%lsd)

         IF (dft_control%do_sccs) THEN
            IF (dft_control%qs_control%gapw) THEN
               CPABORT("SCCS is not yet implemented with GAPW")
            END IF
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
               "SCCS| Self-consistent continuum solvation model"
            SELECT CASE (dft_control%sccs_control%method_id)
            CASE (sccs_andreussi)
               WRITE (UNIT=output_unit, FMT="(T2,A,/,(T2,A,T69,ES12.3))") &
                  "SCCS| Dielectric function proposed by Andreussi et al.", &
                  "SCCS|  rho_max", dft_control%sccs_control%rho_max, &
                  "SCCS|  rho_min", dft_control%sccs_control%rho_min
            CASE (sccs_fattebert_gygi)
               WRITE (UNIT=output_unit, FMT="(T2,A,/,(T2,A,T69,ES12.3))") &
                  "SCCS| Dielectric function proposed by Fattebert and Gygi", &
                  "SCCS|  beta", dft_control%sccs_control%beta, &
                  "SCCS|  rho_zero", dft_control%sccs_control%rho_zero
            CASE DEFAULT
               CPABORT("Invalid SCCS model specified. Please, check your input!")
            END SELECT
            SELECT CASE (dft_control%sccs_control%derivative_method)
            CASE (sccs_derivative_fft)
               WRITE (UNIT=output_unit, FMT="(T2,A,T46,A35)") "SCCS| Numerical derivative calculation", &
                  ADJUSTR("FFT")
            CASE (sccs_derivative_cd3)
               WRITE (UNIT=output_unit, FMT="(T2,A,T46,A35)") "SCCS| Numerical derivative calculation", &
                  ADJUSTR("3-point stencil central differences")
            CASE (sccs_derivative_cd5)
               WRITE (UNIT=output_unit, FMT="(T2,A,T46,A35)") "SCCS| Numerical derivative calculation", &
                  ADJUSTR("5-point stencil central differences")
            CASE (sccs_derivative_cd7)
               WRITE (UNIT=output_unit, FMT="(T2,A,T46,A35)") "SCCS| Numerical derivative calculation", &
                  ADJUSTR("7-point stencil central differences")
            CASE DEFAULT
               CALL cp_abort(__LOCATION__, &
                             "Invalid derivative method specified for SCCS model. "// &
                             "Please, check your input!")
            END SELECT
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,F12.3)") &
               "SCCS| Repulsion parameter alpha [mN/m] = [dyn/cm]", &
               cp_unit_from_cp2k(dft_control%sccs_control%alpha_solvent, "mN/m")
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,F12.3)") &
               "SCCS| Dispersion parameter beta [GPa]", &
               cp_unit_from_cp2k(dft_control%sccs_control%beta_solvent, "GPa")
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,F12.3)") &
               "SCCS| Surface tension gamma [mN/m] = [dyn/cm]", &
               cp_unit_from_cp2k(dft_control%sccs_control%gamma_solvent, "mN/m")
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,F12.3)") &
               "SCCS| Mixing parameter applied during the iteration cycle", &
               dft_control%sccs_control%mixing
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,ES12.3)") &
               "SCCS| Tolerance for the convergence of the SCCS iteration cycle", &
               dft_control%sccs_control%eps_sccs
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,I12)") &
               "SCCS| Maximum number of iteration steps", &
               dft_control%sccs_control%max_iter
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,ES12.3)") &
               "SCCS| SCF convergence threshold for starting the SCCS iteration", &
               dft_control%sccs_control%eps_scf
            WRITE (UNIT=output_unit, FMT="(T2,A,T69,ES12.3)") &
               "SCCS| Numerical increment for the cavity surface calculation", &
               dft_control%sccs_control%delta_rho
         END IF

      END IF

      CALL cp_print_key_finished_output(output_unit, logger, dft_section, &
                                        "PRINT%DFT_CONTROL_PARAMETERS")

      CALL timestop(handle)

   END SUBROUTINE write_dft_control

! **************************************************************************************************
!> \brief Write the xTB control parameters to the output unit.
!> \param xtb_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE write_xtb_control(xtb_control, dft_section)
      TYPE(xtb_control_type), POINTER                    :: xtb_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'write_xtb_control', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, output_unit
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      NULLIFY (logger)
      logger => cp_get_default_logger()

      output_unit = cp_print_key_unit_nr(logger, dft_section, &
                                         "PRINT%DFT_CONTROL_PARAMETERS", extension=".Log")

      IF (output_unit > 0) THEN

         WRITE (UNIT=output_unit, FMT="(/,T2,A,T31,A50)") &
            "xTB| Parameter file", ADJUSTR(TRIM(xtb_control%parameter_file_name))
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
            "xTB| Basis expansion STO-NG", xtb_control%sto_ng
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
            "xTB| Basis expansion STO-NG for Hydrogen", xtb_control%h_sto_ng
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,L10)") &
            "xTB| Halogen interaction potential", xtb_control%xb_interaction
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.3)") &
            "xTB| Halogen interaction potential cutoff radius", xtb_control%xb_radius
         WRITE (UNIT=output_unit, FMT="(T2,A,T31,A50)") &
            "xTB| D3 Dispersion: Parameter file", ADJUSTR(TRIM(xtb_control%dispersion_parameter_file))
         WRITE (UNIT=output_unit, FMT="(T2,A,T51,3F10.3)") &
            "xTB| Huckel constants ks kp kd", xtb_control%ks, xtb_control%kp, xtb_control%kd
         WRITE (UNIT=output_unit, FMT="(T2,A,T61,2F10.3)") &
            "xTB| Huckel constants ksp k2sh", xtb_control%ksp, xtb_control%k2sh
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.3)") &
            "xTB| Mataga-Nishimoto exponent", xtb_control%kg
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.3)") &
            "xTB| Repulsion potential exponent", xtb_control%kf
         WRITE (UNIT=output_unit, FMT="(T2,A,T51,3F10.3)") &
            "xTB| Coordination number scaling kcn(s) kcn(p) kcn(d)", &
            xtb_control%kcns, xtb_control%kcnp, xtb_control%kcnd
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.3)") &
            "xTB| Electronegativity scaling", xtb_control%ken
         WRITE (UNIT=output_unit, FMT="(T2,A,T61,2F10.3)") &
            "xTB| Halogen potential scaling kxr kx2", xtb_control%kxr, xtb_control%kx2
         WRITE (UNIT=output_unit, FMT="(/)")

      END IF

      CALL cp_print_key_finished_output(output_unit, logger, dft_section, &
                                        "PRINT%DFT_CONTROL_PARAMETERS")

      CALL timestop(handle)

   END SUBROUTINE write_xtb_control

! **************************************************************************************************
!> \brief Purpose: Write the QS control parameters to the output unit.
!> \param qs_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE write_qs_control(qs_control, dft_section)
      TYPE(qs_control_type), INTENT(IN)                  :: qs_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'write_qs_control', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=20)                                  :: method, quadrature
      INTEGER                                            :: handle, i, igrid_level, ngrid_level, &
                                                            output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: qs_section
      TYPE(section_vals_type), POINTER                   :: print_section_vals, qs_section_vals

      IF (qs_control%semi_empirical) RETURN
      IF (qs_control%dftb) RETURN
      IF (qs_control%xtb) RETURN
      CALL timeset(routineN, handle)
      NULLIFY (logger, print_section_vals, qs_section, qs_section_vals)
      logger => cp_get_default_logger()
      print_section_vals => section_vals_get_subs_vals(dft_section, "PRINT")
      qs_section_vals => section_vals_get_subs_vals(dft_section, "QS")
      CALL section_vals_get(qs_section_vals, section=qs_section)

      NULLIFY (enum, keyword)
      keyword => section_get_keyword(qs_section, "METHOD")
      CALL keyword_get(keyword, enum=enum)
      method = enum_i2c(enum, qs_control%method_id)

      NULLIFY (enum, keyword)
      keyword => section_get_keyword(qs_section, "QUADRATURE")
      CALL keyword_get(keyword, enum=enum)
      quadrature = enum_i2c(enum, qs_control%gapw_control%quadrature)

      output_unit = cp_print_key_unit_nr(logger, print_section_vals, &
                                         "DFT_CONTROL_PARAMETERS", extension=".Log")
      IF (output_unit > 0) THEN
         ngrid_level = SIZE(qs_control%e_cutoff)
         WRITE (UNIT=output_unit, FMT="(/,T2,A,T61,A20)") &
            "QS| Method:", ADJUSTR(method)
         IF (qs_control%pw_grid_opt%spherical) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T61,A)") &
               "QS| Density plane wave grid type", " SPHERICAL HALFSPACE"
         ELSE IF (qs_control%pw_grid_opt%fullspace) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T57,A)") &
               "QS| Density plane wave grid type", " NON-SPHERICAL FULLSPACE"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T57,A)") &
               "QS| Density plane wave grid type", " NON-SPHERICAL HALFSPACE"
         END IF
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
            "QS| Number of grid levels:", SIZE(qs_control%e_cutoff)
         IF (ngrid_level == 1) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.1)") &
               "QS| Density cutoff [a.u.]:", qs_control%e_cutoff(1)
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.1)") &
               "QS| Density cutoff [a.u.]:", qs_control%cutoff
            IF (qs_control%commensurate_mgrids) &
               WRITE (UNIT=output_unit, FMT="(T2,A)") "QS| Using commensurate multigrids"
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.1)") &
               "QS| Multi grid cutoff [a.u.]: 1) grid level", qs_control%e_cutoff(1)
            WRITE (UNIT=output_unit, FMT="(T2,A,I3,A,T71,F10.1)") &
               ("QS|                         ", igrid_level, ") grid level", &
                qs_control%e_cutoff(igrid_level), &
                igrid_level=2, SIZE(qs_control%e_cutoff))
         END IF
         IF (qs_control%pao) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") "QS| PAO active"
         END IF
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.1)") &
            "QS| Grid level progression factor:", qs_control%progression_factor
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,F10.1)") &
            "QS| Relative density cutoff [a.u.]:", qs_control%relative_cutoff
         WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
            "QS| Interaction thresholds: eps_pgf_orb:", &
            qs_control%eps_pgf_orb, &
            "QS|                         eps_filter_matrix:", &
            qs_control%eps_filter_matrix, &
            "QS|                         eps_core_charge:", &
            qs_control%eps_core_charge, &
            "QS|                         eps_rho_gspace:", &
            qs_control%eps_rho_gspace, &
            "QS|                         eps_rho_rspace:", &
            qs_control%eps_rho_rspace, &
            "QS|                         eps_gvg_rspace:", &
            qs_control%eps_gvg_rspace, &
            "QS|                         eps_ppl:", &
            qs_control%eps_ppl, &
            "QS|                         eps_ppnl:", &
            qs_control%eps_ppnl
         IF (qs_control%gapw) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
               "QS| GAPW|                   eps_fit:", &
               qs_control%gapw_control%eps_fit, &
               "QS| GAPW|                   eps_iso:", &
               qs_control%gapw_control%eps_iso, &
               "QS| GAPW|                   eps_svd:", &
               qs_control%gapw_control%eps_svd, &
               "QS| GAPW|                   eps_cpc:", &
               qs_control%gapw_control%eps_cpc
            WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
               "QS| GAPW|   atom-r-grid: quadrature:", &
               ADJUSTR(quadrature)
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
               "QS| GAPW|      atom-s-grid:  max l :", &
               qs_control%gapw_control%lmax_sphere, &
               "QS| GAPW|      max_l_rho0 :", &
               qs_control%gapw_control%lmax_rho0
            IF (qs_control%gapw_control%non_paw_atoms) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "QS| GAPW|      At least one kind is NOT PAW, i.e. it has only soft AO "
            END IF
            IF (qs_control%gapw_control%nopaw_as_gpw) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "QS| GAPW|      The NOT PAW atoms are treated fully GPW"
            END IF
         END IF
         IF (qs_control%gapw_xc) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
               "QS| GAPW_XC|                eps_fit:", &
               qs_control%gapw_control%eps_fit, &
               "QS| GAPW_XC|                eps_iso:", &
               qs_control%gapw_control%eps_iso, &
               "QS| GAPW_XC|                eps_svd:", &
               qs_control%gapw_control%eps_svd
            WRITE (UNIT=output_unit, FMT="(T2,A,T55,A30)") &
               "QS| GAPW_XC|atom-r-grid: quadrature:", &
               enum_i2c(enum, qs_control%gapw_control%quadrature)
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
               "QS| GAPW_XC|   atom-s-grid:  max l :", &
               qs_control%gapw_control%lmax_sphere
         END IF
         IF (qs_control%mulliken_restraint) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
               "QS| Mulliken restraint target", qs_control%mulliken_restraint_control%target
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
               "QS| Mulliken restraint strength", qs_control%mulliken_restraint_control%strength
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,I8)") &
               "QS| Mulliken restraint atoms: ", qs_control%mulliken_restraint_control%natoms
            WRITE (UNIT=output_unit, FMT="(5I8)") qs_control%mulliken_restraint_control%atoms
         END IF
         IF (qs_control%ddapc_restraint) THEN
            DO i = 1, SIZE(qs_control%ddapc_restraint_control)
               ddapc_restraint_control => qs_control%ddapc_restraint_control(i)%ddapc_restraint_control
               IF (SIZE(qs_control%ddapc_restraint_control) .GT. 1) &
                  WRITE (UNIT=output_unit, FMT="(T2,A,T3,I8)") &
                  "QS| parameters for DDAPC restraint number", i
               WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
                  "QS| ddapc restraint target", ddapc_restraint_control%target
               WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
                  "QS| ddapc restraint strength", ddapc_restraint_control%strength
               WRITE (UNIT=output_unit, FMT="(T2,A,T73,I8)") &
                  "QS| ddapc restraint atoms: ", ddapc_restraint_control%natoms
               WRITE (UNIT=output_unit, FMT="(5I8)") ddapc_restraint_control%atoms
               WRITE (UNIT=output_unit, FMT="(T2,A)") "Coefficients:"
               WRITE (UNIT=output_unit, FMT="(5F6.2)") ddapc_restraint_control%coeff
               SELECT CASE (ddapc_restraint_control%functional_form)
               CASE (do_ddapc_restraint)
                  WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
                     "QS| ddapc restraint functional form :", "RESTRAINT"
               CASE (do_ddapc_constraint)
                  WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
                     "QS| ddapc restraint functional form :", "CONSTRAINT"
               CASE DEFAULT
                  CPABORT("Unknown ddapc restraint")
               END SELECT
            END DO
         END IF
         IF (qs_control%s2_restraint) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
               "QS| s2 restraint target", qs_control%s2_restraint_control%target
            WRITE (UNIT=output_unit, FMT="(T2,A,T73,ES8.1)") &
               "QS| s2 restraint strength", qs_control%s2_restraint_control%strength
            SELECT CASE (qs_control%s2_restraint_control%functional_form)
            CASE (do_s2_restraint)
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
                  "QS| s2 restraint functional form :", "RESTRAINT"
               CPABORT("Not yet implemented")
            CASE (do_s2_constraint)
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
                  "QS| s2 restraint functional form :", "CONSTRAINT"
            CASE DEFAULT
               CPABORT("Unknown ddapc restraint")
            END SELECT
         END IF
      END IF
      CALL cp_print_key_finished_output(output_unit, logger, print_section_vals, &
                                        "DFT_CONTROL_PARAMETERS")

      CALL timestop(handle)

   END SUBROUTINE write_qs_control

! **************************************************************************************************
!> \brief reads the input parameters needed for ddapc.
!> \param qs_control ...
!> \param qs_section ...
!> \param ddapc_restraint_section ...
!> \author fschiff
!> \note
!>      either reads DFT%QS%DDAPC_RESTRAINT or PROPERTIES%ET_coupling
!>      if(qs_section is present the DFT part is read, if ddapc_restraint_section
!>      is present ET_COUPLING is read. Avoid having both!!!
! **************************************************************************************************
   SUBROUTINE read_ddapc_section(qs_control, qs_section, ddapc_restraint_section)

      TYPE(qs_control_type), INTENT(INOUT)               :: qs_control
      TYPE(section_vals_type), OPTIONAL, POINTER         :: qs_section, ddapc_restraint_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_ddapc_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, jj, k, n_rep
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rtmplist
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control
      TYPE(section_vals_type), POINTER                   :: ddapc_section

      IF (PRESENT(ddapc_restraint_section)) THEN
         IF (ASSOCIATED(qs_control%ddapc_restraint_control)) THEN
            IF (SIZE(qs_control%ddapc_restraint_control) .GE. 2) &
               CPABORT("ET_COUPLING cannot be used in combination with a normal restraint")
         ELSE
            ddapc_section => ddapc_restraint_section
            ALLOCATE (qs_control%ddapc_restraint_control(1))
         END IF
      END IF

      IF (PRESENT(qs_section)) THEN
         NULLIFY (ddapc_section)
         ddapc_section => section_vals_get_subs_vals(qs_section, &
                                                     "DDAPC_RESTRAINT")
      END IF

      DO i = 1, SIZE(qs_control%ddapc_restraint_control)

         NULLIFY (qs_control%ddapc_restraint_control(i)%ddapc_restraint_control)
         CALL ddapc_control_create(qs_control%ddapc_restraint_control(i)%ddapc_restraint_control)
         ddapc_restraint_control => qs_control%ddapc_restraint_control(i)%ddapc_restraint_control

         CALL section_vals_val_get(ddapc_section, "STRENGTH", i_rep_section=i, &
                                   r_val=ddapc_restraint_control%strength)
         CALL section_vals_val_get(ddapc_section, "TARGET", i_rep_section=i, &
                                   r_val=ddapc_restraint_control%target)
         CALL section_vals_val_get(ddapc_section, "FUNCTIONAL_FORM", i_rep_section=i, &
                                   i_val=ddapc_restraint_control%functional_form)
         CALL section_vals_val_get(ddapc_section, "ATOMS", i_rep_section=i, &
                                   n_rep_val=n_rep)
         CALL section_vals_val_get(ddapc_section, "TYPE_OF_DENSITY", i_rep_section=i, &
                                   i_val=ddapc_restraint_control%density_type)

         jj = 0
         DO k = 1, n_rep
            CALL section_vals_val_get(ddapc_section, "ATOMS", i_rep_section=i, &
                                      i_rep_val=k, i_vals=tmplist)
            DO j = 1, SIZE(tmplist)
               jj = jj + 1
            END DO
         END DO
         IF (jj < 1) CPABORT("Need at least 1 atom to use ddapc constraints")
         ddapc_restraint_control%natoms = jj
         IF (ASSOCIATED(ddapc_restraint_control%atoms)) &
            DEALLOCATE (ddapc_restraint_control%atoms)
         ALLOCATE (ddapc_restraint_control%atoms(ddapc_restraint_control%natoms))
         jj = 0
         DO k = 1, n_rep
            CALL section_vals_val_get(ddapc_section, "ATOMS", i_rep_section=i, &
                                      i_rep_val=k, i_vals=tmplist)
            DO j = 1, SIZE(tmplist)
               jj = jj + 1
               ddapc_restraint_control%atoms(jj) = tmplist(j)
            END DO
         END DO

         IF (ASSOCIATED(ddapc_restraint_control%coeff)) &
            DEALLOCATE (ddapc_restraint_control%coeff)
         ALLOCATE (ddapc_restraint_control%coeff(ddapc_restraint_control%natoms))
         ddapc_restraint_control%coeff = 1.0_dp

         CALL section_vals_val_get(ddapc_section, "COEFF", i_rep_section=i, &
                                   n_rep_val=n_rep)
         jj = 0
         DO k = 1, n_rep
            CALL section_vals_val_get(ddapc_section, "COEFF", i_rep_section=i, &
                                      i_rep_val=k, r_vals=rtmplist)
            DO j = 1, SIZE(rtmplist)
               jj = jj + 1
               IF (jj > ddapc_restraint_control%natoms) &
                  CPABORT("Need the same number of coeff as there are atoms ")
               ddapc_restraint_control%coeff(jj) = rtmplist(j)
            END DO
         END DO
         IF (jj < ddapc_restraint_control%natoms .AND. jj .NE. 0) &
            CPABORT("Need no or the same number of coeff as there are atoms.")
      END DO
      k = 0
      DO i = 1, SIZE(qs_control%ddapc_restraint_control)
         IF (qs_control%ddapc_restraint_control(i)%ddapc_restraint_control%functional_form == &
             do_ddapc_constraint) k = k + 1
      END DO
      IF (k == 2) CALL cp_abort(__LOCATION__, &
                                "Only a single constraint possible yet, try to use restraints instead ")

   END SUBROUTINE read_ddapc_section

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param efield_section ...
! **************************************************************************************************
   SUBROUTINE read_efield_sections(dft_control, efield_section)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: efield_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_efield_sections', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_path_length)                 :: file_name
      INTEGER                                            :: i, io, j, n, unit_nr
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp_vals
      TYPE(efield_type), POINTER                         :: efield
      TYPE(section_vals_type), POINTER                   :: tmp_section

      DO i = 1, SIZE(dft_control%efield_fields)
         NULLIFY (dft_control%efield_fields(i)%efield)
         ALLOCATE (dft_control%efield_fields(i)%efield)
         efield => dft_control%efield_fields(i)%efield
         NULLIFY (efield%envelop_i_vars, efield%envelop_r_vars)
         CALL section_vals_val_get(efield_section, "INTENSITY", i_rep_section=i, &
                                   r_val=efield%strength)

         CALL section_vals_val_get(efield_section, "POLARISATION", i_rep_section=i, &
                                   r_vals=tmp_vals)
         ALLOCATE (efield%polarisation(SIZE(tmp_vals)))
         efield%polarisation = tmp_vals
         CALL section_vals_val_get(efield_section, "PHASE", i_rep_section=i, &
                                   r_val=efield%phase_offset)
         CALL section_vals_val_get(efield_section, "ENVELOP", i_rep_section=i, &
                                   i_val=efield%envelop_id)
         CALL section_vals_val_get(efield_section, "WAVELENGTH", i_rep_section=i, &
                                   r_val=efield%wavelength)

         IF (efield%envelop_id == constant_env) THEN
            ALLOCATE (efield%envelop_i_vars(2))
            tmp_section => section_vals_get_subs_vals(efield_section, "CONSTANT_ENV", i_rep_section=i)
            CALL section_vals_val_get(tmp_section, "START_STEP", &
                                      i_val=efield%envelop_i_vars(1))
            CALL section_vals_val_get(tmp_section, "END_STEP", &
                                      i_val=efield%envelop_i_vars(2))
         ELSE IF (efield%envelop_id == gaussian_env) THEN
            ALLOCATE (efield%envelop_r_vars(2))
            tmp_section => section_vals_get_subs_vals(efield_section, "GAUSSIAN_ENV", i_rep_section=i)
            CALL section_vals_val_get(tmp_section, "T0", &
                                      r_val=efield%envelop_r_vars(1))
            CALL section_vals_val_get(tmp_section, "SIGMA", &
                                      r_val=efield%envelop_r_vars(2))
         ELSE IF (efield%envelop_id == ramp_env) THEN
            ALLOCATE (efield%envelop_i_vars(4))
            tmp_section => section_vals_get_subs_vals(efield_section, "RAMP_ENV", i_rep_section=i)
            CALL section_vals_val_get(tmp_section, "START_STEP_IN", &
                                      i_val=efield%envelop_i_vars(1))
            CALL section_vals_val_get(tmp_section, "END_STEP_IN", &
                                      i_val=efield%envelop_i_vars(2))
            CALL section_vals_val_get(tmp_section, "START_STEP_OUT", &
                                      i_val=efield%envelop_i_vars(3))
            CALL section_vals_val_get(tmp_section, "END_STEP_OUT", &
                                      i_val=efield%envelop_i_vars(4))
         ELSE IF (efield%envelop_id == custom_env) THEN
            tmp_section => section_vals_get_subs_vals(efield_section, "CUSTOM_ENV", i_rep_section=i)
            CALL section_vals_val_get(tmp_section, "EFIELD_FILE_NAME", c_val=file_name)
            CALL open_file(file_name=TRIM(file_name), file_action="READ", file_status="OLD", unit_number=unit_nr)
            !Determine the number of lines in file
            n = 0
            DO WHILE (.TRUE.)
               READ (unit_nr, *, iostat=io)
               IF (io /= 0) EXIT
               n = n + 1
            END DO
            REWIND (unit_nr)
            ALLOCATE (efield%envelop_r_vars(n + 1))
            !Store the timestep of the list in the first entry of the r_vars
            CALL section_vals_val_get(tmp_section, "TIMESTEP", r_val=efield%envelop_r_vars(1))
            !Read the file
            DO j = 2, n + 1
               READ (unit_nr, *) efield%envelop_r_vars(j)
               efield%envelop_r_vars(j) = cp_unit_to_cp2k(efield%envelop_r_vars(j), "volt/m")
            END DO
            CALL close_file(unit_nr)
         END IF
      END DO
   END SUBROUTINE read_efield_sections

! **************************************************************************************************
!> \brief reads the input parameters needed real time propagation
!> \param dft_control ...
!> \param rtp_section ...
!> \author fschiff
! **************************************************************************************************
   SUBROUTINE read_rtp_section(dft_control, rtp_section)

      TYPE(dft_control_type), INTENT(INOUT)              :: dft_control
      TYPE(section_vals_type), POINTER                   :: rtp_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_rtp_section', &
         routineP = moduleN//':'//routineN

      INTEGER, DIMENSION(:), POINTER                     :: tmp

      ALLOCATE (dft_control%rtp_control)
      CALL section_vals_val_get(rtp_section, "MAX_ITER", &
                                i_val=dft_control%rtp_control%max_iter)
      CALL section_vals_val_get(rtp_section, "MAT_EXP", &
                                i_val=dft_control%rtp_control%mat_exp)
      CALL section_vals_val_get(rtp_section, "ASPC_ORDER", &
                                i_val=dft_control%rtp_control%aspc_order)
      CALL section_vals_val_get(rtp_section, "EXP_ACCURACY", &
                                r_val=dft_control%rtp_control%eps_exp)
      CALL section_vals_val_get(rtp_section, "PROPAGATOR", &
                                i_val=dft_control%rtp_control%propagator)
      CALL section_vals_val_get(rtp_section, "EPS_ITER", &
                                r_val=dft_control%rtp_control%eps_ener)
      CALL section_vals_val_get(rtp_section, "INITIAL_WFN", &
                                i_val=dft_control%rtp_control%initial_wfn)
      CALL section_vals_val_get(rtp_section, "HFX_BALANCE_IN_CORE", &
                                l_val=dft_control%rtp_control%hfx_redistribute)
      CALL section_vals_val_get(rtp_section, "APPLY_DELTA_PULSE", &
                                l_val=dft_control%rtp_control%apply_delta_pulse)
      CALL section_vals_val_get(rtp_section, "PERIODIC", &
                                l_val=dft_control%rtp_control%periodic)
      CALL section_vals_val_get(rtp_section, "DENSITY_PROPAGATION", &
                                l_val=dft_control%rtp_control%linear_scaling)
      CALL section_vals_val_get(rtp_section, "MCWEENY_MAX_ITER", &
                                i_val=dft_control%rtp_control%mcweeny_max_iter)
      CALL section_vals_val_get(rtp_section, "ACCURACY_REFINEMENT", &
                                i_val=dft_control%rtp_control%acc_ref)
      CALL section_vals_val_get(rtp_section, "MCWEENY_EPS", &
                                r_val=dft_control%rtp_control%mcweeny_eps)
      CALL section_vals_val_get(rtp_section, "DELTA_PULSE_SCALE", &
                                r_val=dft_control%rtp_control%delta_pulse_scale)
      CALL section_vals_val_get(rtp_section, "DELTA_PULSE_DIRECTION", &
                                i_vals=tmp)
      dft_control%rtp_control%delta_pulse_direction = tmp
      CALL section_vals_val_get(rtp_section, "SC_CHECK_START", &
                                i_val=dft_control%rtp_control%sc_check_start)

   END SUBROUTINE read_rtp_section

! **************************************************************************************************
!> \brief Parses the BLOCK_LIST keywords from the ADMM section
!> \param admm_control ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE read_admm_block_list(admm_control, dft_section)
      TYPE(admm_control_type), POINTER                   :: admm_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_admm_block_list', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: irep, list_size, n_rep
      INTEGER, DIMENSION(:), POINTER                     :: tmplist

      NULLIFY (tmplist)

      CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%BLOCK_LIST", &
                                n_rep_val=n_rep)

      ALLOCATE (admm_control%blocks(n_rep))

      DO irep = 1, n_rep
         CALL section_vals_val_get(dft_section, "AUXILIARY_DENSITY_MATRIX_METHOD%BLOCK_LIST", &
                                   i_rep_val=irep, i_vals=tmplist)
         list_size = SIZE(tmplist)
         ALLOCATE (admm_control%blocks(irep)%list(list_size))
         admm_control%blocks(irep)%list(:) = tmplist(:)
      END DO

   END SUBROUTINE read_admm_block_list

END MODULE cp_control_utils
