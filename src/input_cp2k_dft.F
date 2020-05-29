!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_dft
   USE basis_set_types,                 ONLY: basis_sort_default,&
                                              basis_sort_zet
   USE bibliography,                    ONLY: &
        Andermatt2016, Andreussi2012, Avezac2005, BaniHashemian2016, Becke1988b, Bengtsson1999, &
        Blochl1995, Brelaz1979, Dewar1977, Dewar1985, Dudarev1997, Dudarev1998, Ehrhardt1985, &
        Fattebert2002, Golze2017a, Golze2017b, Guidon2010, Heinzmann1976, Holmberg2017, &
        Holmberg2018, Hunt2003, Iannuzzi2005, Iannuzzi2006, Iannuzzi2007, Kolafa2004, Krack2000, &
        Krack2002, Kuhne2007, Kunert2003, Lippert1997, Lippert1999, Lu2004, Perdew1981, &
        Repasky2002, Rocha2006, Schenter2008, Schiffmann2015, Stewart1982, Stewart1989, &
        Stewart2007, Thiel1992, Tozer1996, VanVoorhis2015, VandeVondele2003, VandeVondele2005a, &
        VandeVondele2005b, VandeVondele2006, Weber2008, Yin2017, Zhao1994
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level,&
                                              silent_print_level
   USE cp_spline_utils,                 ONLY: pw_interp,&
                                              spline3_nopbc_interp,&
                                              spline3_pbc_interp
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        atomic_guess, becke_cutoff_element, becke_cutoff_global, broyden_type_1, &
        broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, broyden_type_2, &
        broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, casci_canonical, &
        cdft_alpha_constraint, cdft_beta_constraint, cdft_charge_constraint, &
        cdft_magnetization_constraint, cholesky_dbcsr, cholesky_inverse, cholesky_off, &
        cholesky_reduce, cholesky_restore, core_guess, diag_block_davidson, diag_block_krylov, &
        diag_filter_matrix, diag_ot, diag_standard, dmft_model, do_admm_aux_exch_func_bee, &
        do_admm_aux_exch_func_default, do_admm_aux_exch_func_none, do_admm_aux_exch_func_opt, &
        do_admm_aux_exch_func_pbex, do_admm_basis_projection, do_admm_blocked_projection, &
        do_admm_blocking_purify_full, do_admm_charge_constrained_projection, &
        do_admm_exch_scaling_merlot, do_admm_exch_scaling_none, do_admm_purify_cauchy, &
        do_admm_purify_cauchy_subspace, do_admm_purify_mcweeny, do_admm_purify_mo_diag, &
        do_admm_purify_mo_no_diag, do_admm_purify_none, do_admm_purify_none_dm, do_arnoldi, &
        do_bch, do_cn, do_ddapc_constraint, do_ddapc_restraint, do_em, do_etrs, do_full_density, &
        do_gapw_gcs, do_gapw_gct, do_gapw_log, do_loc_both, do_loc_crazy, do_loc_direct, &
        do_loc_homo, do_loc_jacobi, do_loc_l1_norm_sd, do_loc_lumo, do_loc_max, do_loc_min, &
        do_loc_none, do_loc_scdm, do_lri_inv, do_lri_inv_auto, do_lri_opt_all, do_lri_opt_coeff, &
        do_lri_opt_exps, do_lri_pseudoinv_diag, do_lri_pseudoinv_svd, do_method_am1, &
        do_method_dftb, do_method_gapw, do_method_gapw_xc, do_method_gpw, do_method_lrigpw, &
        do_method_mndo, do_method_mndod, do_method_ofgpw, do_method_pdg, do_method_pm3, &
        do_method_pm6, do_method_pm6fm, do_method_pnnl, do_method_rigpw, do_method_rm1, &
        do_method_xtb, do_pade, do_potential_coulomb, do_potential_id, do_potential_short, &
        do_potential_truncated, do_ppl_analytic, do_ppl_grid, do_pwgrid_ns_fullspace, &
        do_pwgrid_ns_halfspace, do_pwgrid_spherical, do_s2_constraint, do_s2_restraint, &
        do_se_is_kdso, do_se_is_kdso_d, do_se_is_slater, do_se_lr_ewald, do_se_lr_ewald_gks, &
        do_se_lr_ewald_r3, do_se_lr_none, do_spin_density, do_taylor, ehrenfest, embed_diff, &
        embed_fa, embed_grid_angstrom, embed_grid_bohr, embed_level_shift, embed_none, &
        embed_quasi_newton, embed_resp, embed_steep_desc, eri_method_full_gpw, eri_method_gpw_ht, &
        eri_operator_coulomb, eri_operator_erf, eri_operator_erfc, eri_operator_gaussian, &
        eri_operator_yukawa, gaussian, general_roks, gto_cartesian, gto_spherical, hf_model, &
        high_spin_roks, history_guess, jacobian_fd1, jacobian_fd1_backward, jacobian_fd1_central, &
        jacobian_fd2, jacobian_fd2_backward, kg_color_dsatur, kg_color_greedy, kg_tnadd_atomic, &
        kg_tnadd_embed, kg_tnadd_embed_ri, kg_tnadd_none, ls_2pnt, ls_3pnt, ls_gold, ls_none, &
        mao_basis_ext, mao_basis_orb, mao_basis_prim, mao_projection, mopac_guess, no_excitations, &
        no_guess, numerical, oe_gllb, oe_lb, oe_none, oe_saop, oe_sic, op_loc_berry, op_loc_boys, &
        op_loc_pipek, orb_dx2, orb_dxy, orb_dy2, orb_dyz, orb_dz2, orb_dzx, orb_px, orb_py, &
        orb_pz, orb_s, ot_algo_irac, ot_algo_taylor_or_diag, ot_chol_irac, ot_lwdn_irac, &
        ot_mini_broyden, ot_mini_cg, ot_mini_diis, ot_mini_sd, ot_poly_irac, ot_precond_full_all, &
        ot_precond_full_kinetic, ot_precond_full_single, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, ot_precond_solver_default, &
        ot_precond_solver_direct, ot_precond_solver_inv_chol, ot_precond_solver_update, &
        outer_scf_basis_center_opt, outer_scf_becke_constraint, outer_scf_cdft_constraint, &
        outer_scf_ddapc_constraint, outer_scf_hirshfeld_constraint, outer_scf_none, &
        outer_scf_optimizer_bisect, outer_scf_optimizer_broyden, outer_scf_optimizer_diis, &
        outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls, outer_scf_optimizer_none, &
        outer_scf_optimizer_sd, outer_scf_optimizer_secant, outer_scf_s2_constraint, &
        plus_u_lowdin, plus_u_mulliken, plus_u_mulliken_charges, radius_covalent, radius_default, &
        radius_single, radius_user, radius_vdw, random_guess, real_time_propagation, &
        ref_charge_atomic, ref_charge_mulliken, rel_dkh, rel_none, rel_pot_erfc, rel_pot_full, &
        rel_sczora_mp, rel_trans_atom, rel_trans_full, rel_trans_molecule, rel_zora, &
        rel_zora_full, rel_zora_mp, restart_guess, rsdft_model, sccs_andreussi, &
        sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7, sccs_derivative_fft, &
        sccs_fattebert_gygi, shape_function_density, shape_function_gaussian, sic_ad, sic_eo, &
        sic_list_all, sic_list_unpaired, sic_mauri_spz, sic_mauri_us, sic_none, slater, &
        smear_energy_window, smear_fermi_dirac, smear_list, sparse_guess, tddfpt_davidson, &
        tddfpt_excitations, tddfpt_lanczos, tddfpt_singlet, tddfpt_spin_cons, tddfpt_spin_flip, &
        tddfpt_triplet, use_coulomb, use_diff, use_mom_ref_coac, use_mom_ref_com, &
        use_mom_ref_user, use_mom_ref_zero, use_no, use_restart_wfn, use_rt_restart, use_scf_wfn, &
        wannier_projection, weight_type_mass, weight_type_unit, wfi_aspc_nr, wfi_frozen_method_nr, &
        wfi_linear_p_method_nr, wfi_linear_ps_method_nr, wfi_linear_wf_method_nr, &
        wfi_ps_method_nr, wfi_use_guess_method_nr, wfi_use_prev_p_method_nr, &
        wfi_use_prev_rho_r_method_nr, wfi_use_prev_wf_method_nr, xas_1s_type, xas_2p_type, &
        xas_2s_type, xas_3d_type, xas_3p_type, xas_3s_type, xas_4d_type, xas_4f_type, xas_4p_type, &
        xas_4s_type, xas_dip_len, xas_dip_vel, xas_dscf, xas_none, xas_not_excited, &
        xas_tdp_by_index, xas_tdp_by_kind, xas_tp_fh, xas_tp_flex, xas_tp_hh, xas_tp_xfh, &
        xas_tp_xhh, xes_tp_val
   USE input_cp2k_almo,                 ONLY: create_almo_scf_section
   USE input_cp2k_distribution,         ONLY: create_distribution_section
   USE input_cp2k_ec,                   ONLY: create_ec_section
   USE input_cp2k_field,                ONLY: create_efield_section,&
                                              create_per_efield_section
   USE input_cp2k_kpoints,              ONLY: create_kpoints_section
   USE input_cp2k_ls,                   ONLY: create_ls_scf_section
   USE input_cp2k_mm,                   ONLY: create_dipoles_section,&
                                              create_neighbor_lists_section
   USE input_cp2k_motion_print,         ONLY: add_format_keyword
   USE input_cp2k_poisson,              ONLY: create_poisson_section
   USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
   USE input_cp2k_tb,                   ONLY: create_dftb_control_section,&
                                              create_xtb_control_section
   USE input_cp2k_transport,            ONLY: create_transport_section
   USE input_cp2k_xc,                   ONLY: create_libxc_section,&
                                              create_xc_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              lchar_t,&
                                              logical_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE pw_grids,                        ONLY: do_pw_grid_blocked_false,&
                                              do_pw_grid_blocked_free,&
                                              do_pw_grid_blocked_true
   USE pw_spline_utils,                 ONLY: no_precond,&
                                              precond_spl3_1,&
                                              precond_spl3_2,&
                                              precond_spl3_3,&
                                              precond_spl3_aint,&
                                              precond_spl3_aint2
   USE qs_density_mixing_types,         ONLY: create_mixing_section
   USE qs_fb_input,                     ONLY: create_filtermatrix_section
   USE qs_mom_types,                    ONLY: create_mom_section
   USE string_utilities,                ONLY: newline,&
                                              s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_dft'

   PUBLIC :: create_dft_section
   PUBLIC :: create_bsse_section, create_qs_section
   PUBLIC :: create_scf_section
   PUBLIC :: create_interp_section, create_localize_section
   PUBLIC :: create_ddapc_restraint_section
   PUBLIC :: create_mgrid_section

CONTAINS

! **************************************************************************************************
!> \brief creates the dft section
!> \param section the section to be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_dft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_dft_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DFT", &
                          description="Parameter needed by LCAO DFT programs", &
                          n_keywords=3, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="BASIS_SET_FILE_NAME", &
                          description="Name of the basis set file, may include a path", &
                          usage="BASIS_SET_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t, repeats=.TRUE., &
                          default_lc_val="BASIS_SET", n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_FILE_NAME", &
                          description="Name of the pseudo potential file, may include a path", &
                          usage="POTENTIAL_FILE_NAME <FILENAME>", &
                          default_lc_val="POTENTIAL")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          variants=(/"RESTART_FILE_NAME"/), &
                          description="Name of the wavefunction restart file, may include a path."// &
                          " If no file is specified, the default is to open the file as generated by the wfn restart print key.", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="UKS", &
                          variants=s2a("UNRESTRICTED_KOHN_SHAM", &
                                       "LSD", &
                                       "SPIN_POLARIZED"), &
                          description="Requests a spin-polarized calculation using alpha "// &
                          "and beta orbitals, i.e. no spin restriction is applied", &
                          usage="LSD", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS", &
                          variants=(/"RESTRICTED_OPEN_KOHN_SHAM"/), &
                          description="Requests a restricted open Kohn-Sham calculation", &
                          usage="ROKS", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="MULTIPLICITY", &
                          variants=(/"MULTIP"/), &
                          description="Two times the total spin plus one. "// &
                          "Specify 3 for a triplet, 4 for a quartet, "// &
                          "and so on. Default is 1 (singlet) for an "// &
                          "even number and 2 (doublet) for an odd number "// &
                          "of electrons.", &
                          usage="MULTIPLICITY 3", &
                          default_i_val=0) ! this default value is just a flag to get the above
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CHARGE", &
                          description="The total charge of the system", &
                          usage="CHARGE -1", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EXCITATIONS", &
                          description="If excitations should be calculated", &
                          usage="EXCITATIONS", &
                          enum_c_vals=s2a("NONE", "TDLR", "TDDFPT"), &
                          enum_i_vals=(/no_excitations, tddfpt_excitations, &
                                        tddfpt_excitations/), &
                          default_i_val=no_excitations)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="PLUS_U_METHOD", &
                          description="Method employed for the calculation of the DFT+U contribution", &
                          repeats=.FALSE., &
                          enum_c_vals=s2a("LOWDIN", "MULLIKEN", "MULLIKEN_CHARGES"), &
                          enum_i_vals=(/plus_u_lowdin, plus_u_mulliken, plus_u_mulliken_charges/), &
                          enum_desc=s2a("Method based on Lowdin population analysis "// &
                                        "(computationally expensive, since the diagonalization of the "// &
                                        "overlap matrix is required, but possibly more robust than Mulliken)", &
                                        "Method based on Mulliken population analysis using the net AO and "// &
                                        "overlap populations (computationally cheap method)", &
                                        "Method based on Mulliken gross orbital populations (GOP)"), &
                          n_var=1, &
                          default_i_val=plus_u_mulliken, &
                          usage="METHOD Lowdin")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RELAX_MULTIPLICITY", &
                          variants=(/"RELAX_MULTIP"/), &
                          description="Tolerance in Hartrees. Do not enforce the occupation "// &
                          "of alpha and beta MOs due to the initially "// &
                          "defined multiplicity, but rather follow the Aufbau principle. "// &
                          "A value greater than zero activates this option. "// &
                          "Larger tolerance values increase the probability for a spin flip. "// &
                          "This option is only valid for unrestricted (i.e. spin polarised) "// &
                          "Kohn-Sham (UKS) calculations.", &
                          usage="RELAX_MULTIPLICITY 0.00001", &
                          repeats=.FALSE., &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUBCELLS", &
                          description="Read the grid size for subcell generation in the construction of "// &
                          "neighbor lists.", usage="SUBCELLS 1.5", &
                          n_var=1, default_r_val=2.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AUTO_BASIS", &
                          description="Specify size of automatically generated auxiliary basis sets: "// &
                          "Options={small,medium,large,huge}", &
                          usage="AUTO_BASIS {basis_type} {basis_size}", &
                          type_of_var=char_t, repeats=.TRUE., n_var=-1, default_c_vals=(/"X", "X"/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURFACE_DIPOLE_CORRECTION", &
                          variants=s2a("SURFACE_DIPOLE", &
                                       "SURF_DIP"), &
                          description="For slab calculations with asymmetric geometries, activate the correction of "// &
                          " the electrostatic potential with "// &
                          "by compensating for the surface dipole. Implemented only for slabs with normal "// &
                          " parallel to one Cartesian axis. The normal direction is given by the keyword SURF_DIP_DIR ", &
                          usage="SURF_DIP", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          citations=(/Bengtsson1999/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURF_DIP_DIR", &
                          description="Cartesian axis parallel to surface normal.", &
                          enum_c_vals=s2a("X", "Y", "Z"), &
                          enum_i_vals=(/1, 2, 3/), &
                          enum_desc=s2a("Along x", "Along y", "Along z"), &
                          n_var=1, &
                          default_i_val=3, &
                          usage="SURF_DIP_DIR Z")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURF_DIP_POS", &
                          description="This keyword assigns an user defined position in Angstroms"// &
                          " in the direction normal to the surface (given by SURF_DIP_DIR). "// &
                          "The default value is -1.0_dp which appplies the correction at a position "// &
                          " that has minimum electron density on the grid.", &
                          usage="SURF_DIP_POS -1.0_dp", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURF_DIP_SWITCH", &
                          description="WARNING: Experimental feature under development that will help the "// &
                          "user to switch parameters to facilitate  SCF convergence. In its current form the "// &
                          "surface dipole correction is switched off if the calculation does not converge in "// &
                          "(0.5*MAX_SCF + 1) outer_scf steps. "// &
                          "The default value is .FALSE.", &
                          usage="SURF_DIP_SWITCH .TRUE.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="CORE_CORR_DIP", &
                          description="If the total CORE_CORRECTION is non-zero and surface dipole "// &
                          "correction is switched on, presence of this keyword will adjust electron "// &
                          "density via MO occupation to reflect the total CORE_CORRECTION. "// &
                          "The default value is .FALSE.", &
                          usage="CORE_CORR_DIP .TRUE.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SORT_BASIS", &
                          description="Sort basis sets according to a certain criterion. ", &
                          enum_c_vals=s2a("DEFAULT", "EXP"), &
                          enum_i_vals=(/basis_sort_default, basis_sort_zet/), &
                          enum_desc=s2a("don't sort", "sort w.r.t. exponent"), &
                          default_i_val=basis_sort_default, &
                          usage="SORT_BASIS EXP")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ls_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_almo_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_kg_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ec_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_admm_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_qs_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_tddfpt_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mgrid_section(subsection, create_subsections=.TRUE.)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_relativistic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_low_spin_roks_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_efield_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_per_efield_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ext_pot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_transport_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! ZMP sections to include the external density or v_xc potential
      CALL create_ext_den_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ext_vxc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_poisson_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_kpoints_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_implicit_solv_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_density_fitting_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xas_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xas_tdp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_rtp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_dft_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sccs_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_dft_section

! **************************************************************************************************
!> \brief Implicit Solvation Model
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_implicit_solv_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_implicit_solv_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, subsection, print_key)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SCRF", &
                          description="Adds an implicit solvation model to the DFT calculation."// &
                          " Know also as Self Consistent Reaction Field.", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_OUT", &
                          description="Value of the dielectric constant outside the sphere", &
                          usage="EPS_OUT <REAL>", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LMAX", &
                          description="Maximum value of L used in the multipole expansion", &
                          usage="LMAX <INTEGER>", &
                          default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_sphere_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing basic info about the method", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_implicit_solv_section

! **************************************************************************************************
!> \brief Create Sphere cavity
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_sphere_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_sphere_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword, subsection)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SPHERE", &
                          description="Treats the implicit solvent environment like a sphere", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="RADIUS", &
                          description="Value of the spherical cavity in the dielectric medium", &
                          usage="RADIUS <REAL>", &
                          unit_str="angstrom", &
                          type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_center_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_sphere_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_center_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_center_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CENTER", &
                          description="Defines the center of the sphere.", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="XYZ", &
                          description="Coordinates of the center of the sphere", &
                          usage="XYZ <REAL> <REAL> <REAL>", &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Defines a list of atoms to define the center of the sphere", &
                          usage="ATOM_LIST <INTEGER> .. <INTEGER>", &
                          type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WEIGHT_TYPE", &
                          description="Defines the weight used to define the center of the sphere"// &
                          " (if ATOM_LIST is provided)", &
                          usage="WEIGHT (UNIT|MASS)", &
                          enum_c_vals=(/"UNIT", "MASS"/), &
                          enum_i_vals=(/weight_type_unit, weight_type_mass/), &
                          default_i_val=weight_type_unit)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIXED", &
                          description="Specify if the center of the sphere should be fixed or"// &
                          " allowed to move", &
                          usage="FIXED <LOGICAL>", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_center_section

! **************************************************************************************************
!> \brief parameters fo the localization of wavefunctions
!> \param section ...
!> \par History
!>      03.2005 created [MI]
! **************************************************************************************************
   SUBROUTINE create_localize_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_localize_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, print_section, subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      NULLIFY (keyword, print_key)
      CALL section_create(section, __LOCATION__, name="LOCALIZE", &
                          description="Use one of the available methods to define the localization "// &
                          " and possibly to optimize it to a minimum or a maximum.", &
                          n_keywords=8, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the MOS localization procedure", &
                          usage="&LOCALIZE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations used for localization methods", &
                          usage="MAX_ITER 2000", default_i_val=10000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MAX_CRAZY_ANGLE", &
         description="Largest allowed angle for the crazy rotations algorithm (smaller is slower but more stable).", &
         usage="MAX_CRAZY_ANGLE 0.1", unit_str="rad", default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CRAZY_SCALE", &
                          description="scale angles", &
                          usage="CRAZY_SCALE 0.9", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CRAZY_USE_DIAG", &
                          description="Use diagonalization (slow) or pade based calculation of matrix exponentials.", &
                          usage="CRAZY_USE_DIAG ", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="USE_HISTORY", &
         description="Generate an improved initial guess based on a history of results, which is useful during MD."// &
         "Will only work if the number of states to be localized remains constant.", &
         usage="USE_HISTORY ", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_OCCUPATION", &
         description="Tolerance in the occupation number to select only fully occupied orbitals for the rotation", &
         usage="EPS_OCCUPATION 1.E-5", default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUT_ITER_EACH", &
                          description="Every how many iterations of the localization algorithm"// &
                          "(Jacobi) the tolerance value is printed out", &
                          usage="OUT_ITER_EACH 100", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LOCALIZATION", &
                          description="Tolerance used in the convergence criterion of the localization methods.", &
                          usage="EPS_LOCALIZATION 1.0E-2", default_r_val=1.0E-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_OR_MAX", &
                          description="Requires the maximization of the spread of the wfn", &
                          usage="MIN_OR_MAX (SPREADMIN|SPREADMAX)", &
                          enum_c_vals=(/"SPREADMIN", "SPREADMAX"/), &
                          enum_i_vals=(/do_loc_min, do_loc_max/), &
                          default_i_val=do_loc_min)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="METHOD", &
         description="Method of optimization if any", &
         usage="METHOD (JACOBI|CRAZY|DIRECT|L1SD|SCDM|NONE)", &
         enum_c_vals=s2a("NONE", "JACOBI", "CRAZY", "L1SD", "DIRECT", "SCDM"), &
         enum_i_vals=(/do_loc_none, &
                       do_loc_jacobi, &
                       do_loc_crazy, &
                       do_loc_l1_norm_sd, &
                       do_loc_direct, do_loc_scdm/), &
         enum_desc=s2a("No localization is applied", &
                       "Using 2 x 2 rotations of the orbitals, slow but robust", &
                       "A new fast method is applied, might be slightly less robust than jacobi, but usually much faster", &
                       "Steepest descent minimization of an approximate l1 norm", &
                       "Using a direct minimisation approacha", "Use QR factorization"), &
         default_i_val=do_loc_jacobi)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBI_FALLBACK", &
                          description="Use Jacobi method in case no convergence was achieved"// &
                          " by using the crazy rotations method.", &
                          usage="JACOBI_FALLBACK", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBI_REFINEMENT", &
                          description="Use Jacobi method to refine the localisation obtained by SCDM", &
                          usage="JACOBI_REFINEMENT", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the localization from a set of orbitals"// &
                          " read from a localization restart file.", &
                          usage="RESTART", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOCHOMO_RESTART_FILE_NAME", &
                          description="File name where to read the MOS from"// &
                          "which to restart the localization procedure for occupied states", &
                          usage="LOCHOMO_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOCLUMO_RESTART_FILE_NAME", &
                          description="File name where to read the MOS from"// &
                          "which to restart the localization procedure for unoccupied states", &
                          usage="LOCLUMO_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR", &
                          description="Type of opertator which defines the spread functional", &
                          usage="OPERATOR (BERRY|BOYS|PIPEK)", &
                          enum_c_vals=s2a("BERRY", "BOYS", "PIPEK"), &
                          enum_i_vals=(/op_loc_berry, op_loc_boys, op_loc_pipek/), &
                          default_i_val=op_loc_berry)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LIST", &
                          description="Indexes of the occupied wfn to be localized"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LIST_UNOCCUPIED", &
                          description="Indexes of the unoccupied states to be localized, "// &
                          "up to now only valid in combination with GPW. "// &
                          "This keyword has to be present if unoccupied states should be localized. "// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATES", &
                          description="Which states to localize, LUMO up to now only available in GPW", &
                          usage="STATES (HOMO|LUMO|ALL)", &
                          enum_c_vals=s2a("OCCUPIED", "UNOCCUPIED", "ALL"), &
                          enum_i_vals=(/do_loc_homo, do_loc_lumo, do_loc_both/), &
                          default_i_val=do_loc_homo)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ENERGY_RANGE", &
         description="Select the orbitals to be localized within the given energy range."// &
         "This type of selection cannot be added on top of the selection through a LIST. It reads to reals that are"// &
         " lower and higher boundaries of the energy range.", &
         usage=" ENERGY_RANGE lower_bound {real}, higher_bound {real}", &
         repeats=.FALSE., &
         n_var=2, default_r_vals=(/0._dp, 0._dp/), unit_str='eV', &
         type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_section)
      CALL section_create(print_section, __LOCATION__, name="PRINT", &
                          description="Collects all printing options related to the Wannier centers and "// &
                          "properties computed with Wannier centers.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing basic info about the method", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)
      ! Add printing of wannier infos
      CALL print_wanniers(print_section)
      NULLIFY (subsection)
      ! Total Dipoles with wannier
      CALL create_dipoles_section(subsection, "TOTAL_DIPOLE", debug_print_level + 1)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Molecular Dipoles with wannier
      CALL create_dipoles_section(subsection, "MOLECULAR_DIPOLES", debug_print_level + 1)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Molecular Mulipole Moments with wannier
      CALL cp_print_key_section_create(subsection, __LOCATION__, name="MOLECULAR_MOMENTS", &
                                       description="Section controlling the calculation of molecular multipole moments.", &
                                       print_level=debug_print_level + 1, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="ORDER", &
                          description="Maximum order of mulitpoles to be calculated.", &
                          usage=" ORDER {integer}", default_i_val=2, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      !
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Molecular States with wannier
      CALL create_molecular_states_section(subsection)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Wannier States with wannier
      CALL create_wannier_states_section(subsection)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      CALL section_add_subsection(section, print_section)
      CALL section_release(print_section)

   END SUBROUTINE create_localize_section

! **************************************************************************************************
!> \brief Controls the printing of the basic info coming from the LOCALIZE
!>        section
!> \param section ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE print_wanniers(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_wanniers', routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(ASSOCIATED(section))
      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_CUBES", &
                                       description="Controls the printing of the wannier functions ", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LU_BOUNDS", &
                          variants=(/"CUBES_LU"/), &
                          description="The lower and upper index of the states to be printed as cube", &
                          usage="CUBES_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/0, -2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LIST", &
                          description="Indexes of the states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="CUBES_LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_CENTERS", &
                                       description="Controls the printing of the wannier functions", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="", &
                                       unit_str="angstrom")

      CALL keyword_create(keyword, __LOCATION__, name="IONS+CENTERS", &
                          description="prints out the wannier centers together with the particles", &
                          usage="IONS+CENTERS", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL add_format_keyword(keyword, print_key, pos=.TRUE., &
                              description="Specifies the format of the output file when IONS+CENTERS is enabled.")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_SPREADS", &
                                       description="Controls the printing of the wannier functions", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "LOC_RESTART", &
                                       description="Controls the printing of restart file for localized MOS", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE print_wanniers

! **************************************************************************************************
!> \brief Create the print dft section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_print_dft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_dft_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, sub_print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options in DFT code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key, keyword, subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_BANNER", &
                                       description="Controls the printing of the banner of the MM program", &
                                       print_level=silent_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "BASIS_SET_FILE", &
                                       description="Controls the printing of a file with all basis sets used.", &
                                       print_level=high_print_level, filename="LOCAL_BASIS_SETS")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "KINETIC_ENERGY", &
                                       description="Controls the printing of the kinetic energy", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "DERIVATIVES", &
                                       description="Print all derivatives after the DFT calculation", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, &
                          name="ndigits", &
                          description="Specify the number of digits used to print derivatives", &
                          default_i_val=6)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="neighbor_lists", &
                                       description="Controls the printing of the neighbor lists", &
                                       print_level=debug_print_level, filename="", unit_str="angstrom")
      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_orb", &
                          description="Activates the printing of the orbital "// &
                          "orbital neighbor lists, "// &
                          "i.e. the overlap neighbor lists", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_aux_fit", &
                          description="Activates the printing of the orbital "// &
                          "orbital neighbor lists wavefunction fitting basis, "// &
                          "i.e. the overlap neighbor lists", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_aux_fit_vs_orb", &
                          description="Activates the printing of the orbital "// &
                          "orbital mixed neighbor lists of wavefunction fitting basis, "// &
                          "and the orbital basis, i.e. the overlap neighbor lists", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_scp", &
                          description="Activates the printing of the vdW SCP "// &
                          "neighbor lists ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_vdw", &
                          description="Activates the printing of the vdW "// &
                          "neighbor lists (from DFT, DFTB, SE), "// &
                          "i.e. the dispersion neighbor lists", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_cn", &
                          description="Activates the printing of the "// &
                          "neighbor lists used for coordination numbers in vdW DFT-D3", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sac_ae", &
                          description="Activates the printing of the orbital "// &
                          "nuclear attraction neighbor lists (erfc potential)", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sac_ppl", &
                          description="Activates the printing of the orbital "// &
                          "GTH-PPL neighbor lists (local part of the "// &
                          "Goedecker-Teter-Hutter pseudo potentials)", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sap_ppnl", &
                          description="Activates the printing of the orbital "// &
                          "GTH-PPNL neighbor lists (non-local part of the"// &
                          "Goedecker-Teter-Hutter pseudo potentials)", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sap_oce", &
                          description="Activates the printing of the orbital "// &
                          "PAW-projector neighbor lists (only GAPW)", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_se", &
                          description="Activates the printing of the two-center "// &
                          "neighbor lists for Coulomb type interactions in NDDO ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_lrc", &
                          description="Activates the printing of the long-range SE correction "// &
                          "neighbor lists (only when doing long-range SE with integral scheme KDSO and KDSO-d)", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_tbe", &
                          description="Activates the printing of the DFTB Ewald "// &
                          "neighbor lists ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_core", &
                          description="Activates the printing of core interaction "// &
                          "neighbor lists ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_xb", &
                          description="Activates the printing of XB interaction from (xTB) "// &
                          "neighbor lists ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="soo_list", &
                          description="Activates the printing of RI orbital-orbital "// &
                          "neighbor lists ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sip_list", &
                          description="Activates the printing of RI basis-projector interaction "// &
                          "neighbor lists ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "SUBCELL", &
                                       description="Activates the printing of the subcells used for the"// &
                                       "generation of neighbor lists.", unit_str="angstrom", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "AO_MATRICES", &
                                       description="Controls the printing of the ao (i.e. contracted gaussian) matrices (debug).", &
                                       print_level=debug_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="OMIT_HEADERS", &
                          description="Print only the matrix data, not the row and column headers", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specify the number of digits used to print the AO matrices", &
                          default_i_val=6)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CORE_HAMILTONIAN", &
                          description="If the printkey is activated controls the printing of the hamiltonian matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DENSITY", &
                          description="If the printkey is activated controls the printing of the density (P) matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KINETIC_ENERGY", &
                          description="If the printkey is activated controls the printing of the kinetic energy matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="KOHN_SHAM_MATRIX", &
                          description="If the printkey is activated controls the printing of the kohn-sham matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="MATRIX_VXC", &
         description="If the printkey is activated compute and print the matrix of the exchange and correlation potential."// &
         "Only the GGA part for GPW is printed", &
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ORTHO", &
                          description="If the printkey is activated controls the printing of the orthogonalization matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OVERLAP", &
                          description="If the printkey is activated controls the printing of the overlap matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="COMMUTATOR_HR", &
                          description="Controls the printing of the [H,r] commutator matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FERMI_CONTACT", &
                          description="If the printkey is activated controls the printing of the Fermi contact matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="PSO", &
         description="If the printkey is activated controls the printing of the paramagnetic spin-orbit matrices", &
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="EFG", &
         description="If the printkey is activated controls the printing of the electric field gradient matrices", &
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_ENERGY", &
                          description="If the printkey is activated controls the printing of the potential energy matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OCE_HARD", &
                          description="If the printkey is activated controls the printing of the OCE HARD matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OCE_SOFT", &
                          description="If the printkey is activated controls the printing of the OCE SOFT matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="W_MATRIX", &
                          description="If the printkey is activated controls the printing of the w matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="W_MATRIX_AUX_FIT", &
                          description="If the printkey is activated controls the printing of the w matrix", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DERIVATIVES", &
                          description="If the printkey is activated controls the printing "// &
                          "of derivatives (for the matrixes that support this)", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "MO", &
         description="Controls the printing of the molecular orbitals. "// &
         "Note that this is only functional with diagonalization based methods, in particular not with OT (see MO_CUBES)", &
         print_level=high_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="Cartesian", &
                          description="If the printkey is activated controls the printing of the mo in the cartesian basis", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EIGENVALUES", variants=s2a("EIGVALS"), &
                          description="If the printkey is activated controls the printing of the eigenvalues of the mos", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EIGENVECTORS", variants=s2a("EIGVECS"), &
                          description="If the printkey is activated controls the printing of the eigenvectors of the mos", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OCCUPATION_NUMBERS", variants=s2a("OCCNUMS"), &
                          description="If the printkey is activated controls the printing of the occupation numbers of the mos", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specify the number of digits used to print the MO eigenvalues and occupation numbers", &
                          default_i_val=6)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="MO_INDEX_RANGE", &
                          variants=s2a("MO_RANGE", "RANGE"), &
                          description="Allows to print only a subset of the MO eigenvectors or eigenvalues. "// &
                          "The indices of the first and the last MO have to be specified", &
                          repeats=.FALSE., &
                          n_var=2, &
                          type_of_var=integer_t, &
                          default_i_vals=(/0, 0/), &
                          usage="MO_INDEX_RANGE 10 15")
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_MOLDEN", &
                                       description="Write the molecular orbitals in Molden file format, for visualisation.", &
                                       print_level=debug_print_level + 1, add_last=add_last_numeric, filename="MOS")
      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specifies the number of significant digits retained. 3 is OK for visualization.", &
                          usage="NDIGITS {int}", &
                          default_i_val=3)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GTO_KIND", &
                          description="Representation of Gaussian-type orbitals", &
                          default_i_val=gto_spherical, &
                          enum_c_vals=s2a("CARTESIAN", "SPHERICAL"), &
                          enum_desc=s2a( &
                          "Cartesian Gaussian orbitals. Use with caution", &
                          "Spherical Gaussian orbitals. Incompatible with VMD"), &
                          enum_i_vals=(/gto_cartesian, gto_spherical/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL create_mo_cubes_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL create_stm_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL create_wfn_mix_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_active_space_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="GAPW", &
                          description="Controls the printing of some gapw related information (debug).", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "projectors", &
                                       description="If the printkey is activated controls if information on"// &
                                       " the projectors is printed.", &
                                       print_level=debug_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "rho0_information", &
                                       description="If the printkey is activated controls if information on rho0 is printed.", &
                                       print_level=debug_print_level, filename="__STD_OUT__", unit_str="angstrom")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "dft_control_parameters", &
                                       description="Controls the printing of dft control parameters.", &
                                       print_level=medium_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "KPOINTS", &
                                       description="Controls the printing of kpoint information.", &
                                       print_level=medium_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (subsection)
      CALL create_bandstructure_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "OVERLAP_CONDITION", &
                                       description="Controls the the checking and printing of and estimate "// &
                                       "of the overlap matrix condition number", &
                                       print_level=debug_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="1-NORM", &
                          description="Calculate an estimate of the 1-norm condition number", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DIAGONALIZATION", &
                          description="Calculate the 1- and 2-norm condition numbers using diagonalization", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ARNOLDI", &
                          description="Calculate the 2-norm condition number using the Arnoldi code (may not be reliable)", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="E_DENSITY_CUBE", &
                                       description="Controls the printing of cube files with "// &
                                       "the electronic density and, for LSD "// &
                                       "calculations, the spin density", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="TOTAL_DENSITY", &
                          description="Print the total electronic density in the case "// &
                          "of a GAPW run. This keyword has only an effect, "// &
                          "if PAW atoms are present. The default is to print "// &
                          "only the soft part of the electronic density and "// &
                          "to ignore the hard part. NOTE: The total density in "// &
                          "real space might exhibit unphysical features like "// &
                          "spikes due to the finite and thus truncated g vector "// &
                          "expansion", &
                          usage="TOTAL_DENSITY {logical}", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XRD_INTERFACE", &
                          description="It activates the print out of exponents and coefficients for the "// &
                          " Gaussian expansion of the core densities, based on atom calculations for each kind."// &
                          " The resulting core dansities are needed to compute the form factors."// &
                          " If GAPW the local densities are also given in terms of a Gaussian expansion,"// &
                          " by fitting the difference between local-fhard and local-soft density for each atom."// &
                          " In this case the keyword TOTAL_DENSITY is set to FALSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NGAUSS", &
                          description="Number of Gaussian functions used in the expansion of atomic (core) density", &
                          usage="NGAUSS 10", n_var=1, default_i_val=12, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "tot_density_cube", &
                                       description="Controls printing of cube files with "// &
                                       "the total density (electrons+atomic core). Note that "// &
                                       "the value of the total density is positive where the "// &
                                       "electron density dominates and negative where the core is.", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "v_hartree_cube", &
                                       description="Controls the printing of a cube file with eletrostatic "// &
                                       " potential generated by the total density (electrons+ions). It is "// &
                                       " valid only for QS with GPW formalism."// &
                                       " Note that by convention the potential has opposite sign than the expected physical one.", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "external_potential_cube", &
                                       description="Controls the printing of a cube file with external "// &
                                       " potential from the DFT%EXTERNAL_POTENTIAL section only.", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

! cube files for data generated by the implicit (generalized) Poisson solver
      CALL create_implicit_psolver_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

! ZMP adding the print section for  the v_xc cube
      CALL cp_print_key_section_create(print_key, __LOCATION__, "v_xc_cube", &
                                       description="Controls the printing of a cube file with xc "// &
                                       " potential generated by the ZMP method (for the moment). It is "// &
                                       " valid only for QS with GPW formalism .", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "efield_cube", &
                                       description="Controls the printing of cube files with electric "// &
                                       " field generated by the total density (electrons+ions). It is "// &
                                       " valid only for QS with GPW formalism .", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "ELF_CUBE", &
                                       description="Controls printing of cube files with "// &
                                       "the electron localization function (ELF). Note that "// &
                                       "the value of ELF is defined between 0 and 1:  Pauli kinetic energy density normalized "// &
                                       " by the kinetic energy density of a uniform el. gas of same density.", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="density_cutoff", &
                          description=" ", &
                          usage="density_cutoff 0.0001", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.0e-10_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "LOCAL_ENERGY_CUBE", &
                                       description="Controls the printing of cube files with the local "// &
                                       " energy. It is valid only for QS with GPW/GAPW formalism ."// &
                                       " Meta and hybrid functionals are not possible.", &
                                       print_level=debug_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "LOCAL_STRESS_CUBE", &
                                       description="Controls the printing of cube files with the local "// &
                                       " stress. It is valid only for QS with GPW/GAPW formalism ."// &
                                       " Meta and hybrid functionals are not possible.", &
                                       print_level=debug_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL create_dos_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL create_pdos_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL create_wannier_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      !Printing of Moments
      CALL create_dipoles_section(print_key, "MOMENTS", high_print_level)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MAX_MOMENT", &
         description="Maximum moment to be calculated. Values higher than 1 not implemented under periodic boundaries.", &
         usage="MAX_MOMENT {integer}", &
         repeats=.FALSE., &
         n_var=1, &
         type_of_var=integer_t, &
         default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="MAGNETIC", &
                          description="Calculate also magnetic moments, only implemented without periodic boundaries", &
                          usage="MAGNETIC yes", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="VEL_REPRS", &
                          description="Calculate expectation values of the el. multipole moments in their velocity"// &
                          "representation during RTP. Implemented up to el. quadrupole moment.", &
                          usage="VEL_REPS yes", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="COM_NL", &
                          description="Include non local commutator for velocity representations."// &
                          "Necessary for origin independent results.", &
                          usage="COM_NL yes", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="SECOND_REFERENCE_POINT", &
                          description="Use second reference point", &
                          usage="SECOND_REFERENCE_POINT .TRUE.", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_2", &
                          variants=s2a("REF_2"), &
                          description="Define a second reference point for the calculation of the electrostatic moment.", &
                          usage="REFERENCE_2 COM", &
                          enum_c_vals=s2a("COM", "COAC", "USER_DEFINED", "ZERO"), &
                          enum_desc=s2a("Use Center of Mass", &
                                        "Use Center of Atomic Charges", &
                                        "Use User Defined Point (Keyword:REF_POINT)", &
                                        "Use Origin of Coordinate System"), &
                          enum_i_vals=(/use_mom_ref_com, &
                                        use_mom_ref_coac, &
                                        use_mom_ref_user, &
                                        use_mom_ref_zero/), &
                          default_i_val=use_mom_ref_zero)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_POINT_2", &
                          variants=s2a("REF_POINT_2"), &
                          description="Fixed second reference point for the calculations of the electrostatic moment.", &
                          usage="REFERENCE_POINT_2 x y z", &
                          repeats=.FALSE., &
                          n_var=3, default_r_vals=(/0._dp, 0._dp, 0._dp/), &
                          type_of_var=real_t, &
                          unit_str='bohr')
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Mulliken population analysis
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MULLIKEN", &
                                       description="Controls the printing of the Mulliken (spin) population analysis", &
                                       print_level=medium_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="PRINT_GOP", &
         description="Print the gross orbital populations (GOP) in addition to the gross atomic populations (GAP) "// &
         "and net charges", &
         usage="PRINT_GOP yes", &
         repeats=.FALSE., &
         n_var=1, &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="PRINT_ALL", &
         description="Print all information including the full net AO and overlap population matrix", &
         usage="PRINT_ALL yes", &
         repeats=.FALSE., &
         n_var=1, &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Lowdin population analysis (fairly expensive to compute, so only at high)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "LOWDIN", &
                                       description="Controls the printing of the Lowdin (spin) population analysis", &
                                       print_level=high_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="PRINT_GOP", &
         description="Print the orbital populations in addition to the atomic populations and net charges", &
         usage="PRINT_GOP yes", &
         repeats=.FALSE., &
         n_var=1, &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="PRINT_ALL", &
         description="Print all information including the full symmetrically orthogonalised density matrix", &
         usage="PRINT_ALL yes", &
         repeats=.FALSE., &
         n_var=1, &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Hirshfeld population analysis
      CALL cp_print_key_section_create(print_key, __LOCATION__, "HIRSHFELD", &
                                       description="Controls the printing of the Hirshfeld (spin) population analysis", &
                                       print_level=medium_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1)
      CALL keyword_create(keyword, __LOCATION__, name="SELF_CONSISTENT", &
                          description="Calculate charges from the Hirscheld-I (self_consistent) method."// &
                          " This scales only the full shape function, not the added charge as in the original scheme.", &
                          usage="SELF_CONSISTENT yes", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SHAPE_FUNCTION", &
                          description="Type of shape function used for Hirshfeld partitioning.", &
                          usage="SHAPE_FUNCTION {Gaussian,Density}", repeats=.FALSE., n_var=1, &
                          default_i_val=shape_function_gaussian, &
                          enum_c_vals=s2a("GAUSSIAN", "DENSITY"), &
                          enum_desc=s2a("Single Gaussian with Colvalent radius", &
                                        "Atomic density expanded in multiple Gaussians"), &
                          enum_i_vals=(/shape_function_gaussian, shape_function_density/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_CHARGE", &
                          description="Charge of atomic partitioning function for Hirshfeld method.", &
                          usage="REFERENCE_CHARGE {Atomic,Mulliken}", repeats=.FALSE., n_var=1, &
                          default_i_val=ref_charge_atomic, &
                          enum_c_vals=s2a("ATOMIC", "MULLIKEN"), &
                          enum_desc=s2a("Use atomic core charges", "Calculate Mulliken charges"), &
                          enum_i_vals=(/ref_charge_atomic, ref_charge_mulliken/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="USER_RADIUS", &
                          description="Use user defined radii to generate Gaussians."// &
                          " These radii are defined by the keyword ATOMIC_RADII", &
                          usage="USER_RADIUS yes", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_RADII", &
                          description="Defines custom radii to setup the spherical Gaussians.", &
                          usage="ATOMIC_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! MAO (modified atomic orbital) analysis
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MAO_ANALYSIS", &
                                       description="Controls the printing of the MAO (modified atomic orbital) analysis", &
                                       print_level=debug_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1, &
                                       citations=(/Heinzmann1976, Ehrhardt1985/))
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold for matrix elements in MAO determination.", &
                          usage="EPS_FILTER reps", repeats=.FALSE., n_var=1, &
                          default_r_val=1.e-8_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_BASIS", &
                          description="Basis set used to construct MAO's.", &
                          usage="REFERENCE_BASIS {ORBITAL,PRIMITIVE,EXTERNAL}", repeats=.FALSE., n_var=1, &
                          default_i_val=mao_basis_orb, &
                          enum_c_vals=s2a("ORBITAL", "PRIMITIVE", "EXTERNAL"), &
                          enum_desc=s2a("Use standard orbital basis set", "Construct basis from primitives of the orbital basis", &
                                        "Read external basis (MAO)"), &
                          enum_i_vals=(/mao_basis_orb, mao_basis_prim, mao_basis_ext/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PRINT_BASIS", &
                          description="Print out MAO reference basis.", &
                          usage="PRINT_BASIS {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_GRAD", &
                          description="Threshold for gradient in MAO optimization.", &
                          usage="EPS_GRAD reps", repeats=.FALSE., n_var=1, &
                          default_r_val=1.e-4_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FUNCTION", &
                          description="Threshold for electron defect in MAO optimization.", &
                          usage="EPS_FUNCTION feps", repeats=.FALSE., n_var=1, &
                          default_r_val=1.e-3_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum allowed iterations for MAO optimization.", &
                          usage="MAX_ITER iter", repeats=.FALSE., n_var=1, &
                          default_i_val=0, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NEGLECT_ABC", &
                          description="Neglect 3 atom terms in MAO analysis.", &
                          usage="NEGLECT_ABC {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="AB_THRESHOLD", &
                          description="Threshold for printing of AB shared electron numbers.", &
                          usage="AB_THRESHOLD thr", repeats=.FALSE., n_var=1, &
                          default_r_val=1.e-2_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ABC_THRESHOLD", &
                          description="Threshold for printing of ABC shared electron numbers.", &
                          usage="ABC_THRESHOLD thr", repeats=.FALSE., n_var=1, &
                          default_r_val=1.e-5_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ANALYZE_UNASSIGNED_CHARGE", &
                          description="Calculate atomic contributions to the unassigned charge.", &
                          usage="ANALYZE_UNASSIGNED_CHARGE {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      !Minimal localized basis analysis
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MINBAS_ANALYSIS", &
                                       description="Controls the printing of the minimal localized basis analysis", &
                                       print_level=debug_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1, &
                                       citations=(/Lu2004/))
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold for matrix elements in basis determination.", &
                          usage="EPS_FILTER reps", repeats=.FALSE., n_var=1, &
                          default_r_val=1.e-8_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FULL_ORTHOGONALIZATION", &
                          description="Orthogonalize the localized minimal basis.", &
                          usage="FULL_ORTHOGONALIZATION {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="BOND_ORDER", &
                          description="Calculate Mayer Bond Orders.", &
                          usage="BOND_ORDER {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      NULLIFY (sub_print_key)
      CALL cp_print_key_section_create(sub_print_key, __LOCATION__, "MINBAS_CUBE", &
                                       description="Write the minimal basis on Cube files.", &
                                       print_level=debug_print_level + 1, add_last=add_last_numeric, filename="MINBAS")
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Indexes of the atoms minimal basis to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="CUBES_LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      NULLIFY (sub_print_key)
      CALL cp_print_key_section_create(sub_print_key, __LOCATION__, "MINBAS_MOLDEN", &
                                       description="Write the minimal basis in Molden file format, for visualisation.", &
                                       print_level=debug_print_level + 1, add_last=add_last_numeric, filename="MINBAS")
      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specifies the number of significant digits retained. 3 is OK for visualization.", &
                          usage="NDIGITS {int}", &
                          default_i_val=3)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GTO_KIND", &
                          description="Representation of Gaussian-type orbitals", &
                          default_i_val=gto_spherical, &
                          enum_c_vals=s2a("CARTESIAN", "SPHERICAL"), &
                          enum_desc=s2a( &
                          "Cartesian Gaussian orbitals. Use with caution", &
                          "Spherical Gaussian orbitals. Incompatible with VMD"), &
                          enum_i_vals=(/gto_cartesian, gto_spherical/))
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      !DOS from density matrix
      CALL cp_print_key_section_create(print_key, __LOCATION__, "ENERGY_WINDOWS", &
                                       description="Controls the printing of the DOS from the density matrix. "// &
                                       "This allows the calculation of the DOS even in density matrix based"// &
                                       "REAL_TIME_PROPAGATION and LS_SCF. "// &
                                       "However, it requires a cubically scaling diagonalization of the Hamiltonian. "// &
                                       "Hartree-Fock NYI, values will be wrong. "// &
                                       "Careful, the orbitals in rtp/emd are not actually eigenstates of the Hamiltonian."// &
                                       "Assumes absence of spin polarization (so far).", &
                                       print_level=high_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("MD"), each_iter_values=(/100/), &
                                       add_last=add_last_numeric, filename="energy-windows")
      CALL keyword_create(keyword, __LOCATION__, name="N_WINDOWS", &
                          description="The number of energy windows.", &
                          usage="N_WINDOWS 200", &
                          default_i_val=100)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Filtering threshold for sparse matrix operations.", &
                          usage="EPS_FILTER 1.0E-6", &
                          default_r_val=1.0E-14_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="RESTRICT_RANGE", &
                          description="Restricts the energy windows to states close to the fermi level", &
                          usage="RESTRICT_RANGE .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="RANGE", &
                          description="If the RESTRICT_RANGE keyword is set, then all energy widnows will"// &
                          " be placed in an interval from from the fermi level minus to the fermi level plus this keyword", &
                          usage="RANGE 1", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PRINT_CUBES", &
                          description="Print the energy windows to cube files", &
                          usage="DENSITY_PROPAGATION .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the energy windows cube files (if enabled) "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Hamiltonian in CSR format
      CALL cp_print_key_section_create(print_key, __LOCATION__, "KS_CSR_WRITE", &
                                       description="Write the KS matrix in CSR format into a file.", &
                                       print_level=debug_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="Threshold", &
                          description="Threshold on the absolute value of the elements to be printed out. "// &
                          "In CP2K all the elements in a (atomic) matrix block are considered non-zero, "// &
                          "if the block contains at least one non-zero element.", &
                          usage="THRESHOLD {real}", &
                          repeats=.FALSE., &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="Upper_triangular", &
                          description="Print only the upper triangular part of the matrix. ", &
                          usage="UPPER_TRIANGULAR {logical}", &
                          repeats=.FALSE., &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="Binary", &
                          description="Whether or not to generate the file in a binary format. ", &
                          usage="BINARY {logical}", &
                          repeats=.FALSE., &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Overlap in CSR format
      CALL cp_print_key_section_create(print_key, __LOCATION__, "S_CSR_WRITE", &
                                       description="Write the overlap matrix in CSR format into a file.", &
                                       print_level=debug_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="Threshold", &
                          description="Threshold on the absolute value of the elements to be printed out. "// &
                          "In CP2K all the elements in a (atomic) matrix block are considered non-zero, "// &
                          "if the block contains at least one non-zero element.", &
                          usage="THRESHOLD {real}", &
                          repeats=.FALSE., &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="Upper_triangular", &
                          description="Print only the upper triangular part of the matrix. ", &
                          usage="UPPER_TRIANGULAR {logical}", &
                          repeats=.FALSE., &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="Binary", &
                          description="Whether or not to generate the file in a binary format. ", &
                          usage="BINARY {logical}", &
                          repeats=.FALSE., &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! interaction adjacency matrix
      CALL cp_print_key_section_create(print_key, __LOCATION__, "ADJMAT_WRITE", &
                                       description="Writes an (upper-triangular) adjacency matrix indicating the "// &
                                       "interaction between atoms (according to overlapping basis functions). The "// &
                                       "columns are: iatom, jatom, ikind, jkind; where iatom and jatom are the atom "// &
                                       "indices (based on the coordinate file), ikind and jkind are the atomic kinds "// &
                                       "(indeces as shown in the ATOMIC KIND INFORMATION section of a CP2K output). ", &
                                       print_level=debug_print_level, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Xray diffraction
      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, name="XRAY_DIFFRACTION_SPECTRUM", &
         description="Calculate and print the coherent X-ray "// &
         "diffraction spectrum", &
         print_level=debug_print_level, &
         filename="", &
         citations=(/Krack2000, Krack2002/))
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="Q_MAX", &
         variants=(/"Q_MAXIMUM"/), &
         description="Maximum Q value calculated for the spectrum", &
         usage="Q_MAX {real}", &
         repeats=.FALSE., &
         n_var=1, &
         type_of_var=real_t, &
         default_r_val=cp_unit_to_cp2k(value=20.0_dp, &
                                       unit_str="angstrom^-1"), &
         unit_str="angstrom^-1")
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="ELECTRIC_FIELD_GRADIENT", &
                                       description="Calculate and print the electric field gradients "// &
                                       "at atomic positions", &
                                       print_level=debug_print_level, &
                                       filename="__STD_OUT__")

      CALL keyword_create(keyword, __LOCATION__, &
                          name="INTERPOLATION", &
                          description="Use interpolation method from real space grid", &
                          usage="INTERPOLATION {logical}", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="GSPACE_SMOOTHING", &
                          description="Use a G-space smoothing function", &
                          usage="GSPACE_SMOOTHING cutoff {real}, width {real}", &
                          repeats=.FALSE., &
                          n_var=2, default_r_vals=(/-1._dp, -1._dp/), &
                          type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DEBUG", &
                          description="Print additional debug output", &
                          usage="DEBUG {logical}", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL create_gspace_interp_section(subsection)
      CALL section_add_subsection(print_key, subsection)
      CALL section_release(subsection)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="BASIS_MOLOPT_QUANTITIES", &
                                       description="Print the two quantities needed in the basis molopt generation:"// &
                                       " total energy and condition number of the overlap matrix (S matrix)", &
                                       print_level=debug_print_level, &
                                       filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="HYPERFINE_COUPLING_TENSOR", &
                                       description="Calculate and print the EPR hyperfine coupling tensor"// &
                                       " at atomic positions", &
                                       print_level=debug_print_level, &
                                       filename="__STD_OUT__")

      CALL keyword_create(keyword, __LOCATION__, &
                          name="INTERACTION_RADIUS", &
                          description="Radius of interaction for EPR hyperfine tensor calculation", &
                          usage="INTERACTION_RADIUS radius {real}", &
                          repeats=.FALSE., &
                          n_var=1, default_r_val=10._dp, &
                          type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="OPTIMIZE_LRI_BASIS", &
                                       description="Optimize the exponents of the LRI basis set", &
                                       print_level=low_print_level, &
                                       filename="OPTIMIZED_LRI_BASIS")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, name="PLUS_U", &
         description="Controls the printing for the DFT+U methods", &
         print_level=high_print_level, &
         filename="__STD_OUT__", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/), &
         citations=(/Dudarev1997, Dudarev1998/))
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, name="SCCS", &
         description="Controls the printing for the SCCS models", &
         print_level=high_print_level, &
         filename="__STD_OUT__", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/), &
         citations=(/Fattebert2002, Andreussi2012, Yin2017/))

      NULLIFY (sub_print_key)

      CALL cp_print_key_section_create( &
         sub_print_key, __LOCATION__, name="DENSITY_GRADIENT", &
         description="Controls the printing of the cube files with "// &
         "the norm of the density gradient |&nabla;&rho;| "// &
         "used by the SCCS model.", &
         print_level=debug_print_level, &
         filename="", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/))
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 "// &
                          "numbers (for X,Y,Z) or 1 number valid for all components.", &
                          n_var=-1, &
                          default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t, &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the cube files when they already exist", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL cp_print_key_section_create( &
         sub_print_key, __LOCATION__, name="DIELECTRIC_FUNCTION", &
         description="Controls the printing of the cube files with "// &
         "the dielectric function used by the SCCS model. "// &
         "This function determines the cavity formed by a solute in "// &
         "a solvent and thus it can be used for the visualisaton of the cavity.", &
         print_level=debug_print_level, &
         filename="", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/), &
         citations=(/Fattebert2002, Andreussi2012, Yin2017/))
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 "// &
                          "numbers (for X,Y,Z) or 1 number valid for all components.", &
                          n_var=-1, &
                          default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t, &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the cube files when they already exist", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL cp_print_key_section_create( &
         sub_print_key, __LOCATION__, name="TOTAL_CHARGE_DENSITY", &
         description="Controls the printing of the cube files with the "// &
         "total charge density &rho;<sup>tot</sup> used by the SCCS model.", &
         print_level=debug_print_level, &
         filename="", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/), &
         citations=(/Fattebert2002, Andreussi2012, Yin2017/))
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 "// &
                          "numbers (for X,Y,Z) or 1 number valid for all components.", &
                          n_var=-1, &
                          default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t, &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the cube files when they already exist", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL cp_print_key_section_create( &
         sub_print_key, __LOCATION__, name="POLARISATION_CHARGE_DENSITY", &
         description="Controls the printing of the cube files with the "// &
         "polarisation charge density &rho;<sup>pol</sup> used by the SCCS model with the "// &
         "total charge density &rho;<sup>tot</sup> = &rho;<sup>sol</sup> + &rho;<sup>pol</sup>", &
         print_level=debug_print_level, &
         filename="", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/), &
         citations=(/Fattebert2002, Andreussi2012, Yin2017/))
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 "// &
                          "numbers (for X,Y,Z) or 1 number valid for all components.", &
                          n_var=-1, &
                          default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t, &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the cube files when they already exist", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL cp_print_key_section_create( &
         sub_print_key, __LOCATION__, name="POLARISATION_POTENTIAL", &
         description="Controls the printing of the cube files with the "// &
         "polarisation potential &phi;<sup>pol</sup> used by the SCCS model with the "// &
         "total potential &phi;<sup>tot</sup> = &phi;<sup>sol</sup> + &phi;<sup>pol</sup>", &
         print_level=debug_print_level, &
         filename="", &
         each_iter_names=s2a("QS_SCF"), &
         each_iter_values=(/0/), &
         citations=(/Fattebert2002, Andreussi2012, Yin2017/))
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 "// &
                          "numbers (for X,Y,Z) or 1 number valid for all components.", &
                          n_var=-1, &
                          default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t, &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the cube files when they already exist", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_dft_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE create_bandstructure_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_bandstructure_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BAND_STRUCTURE", &
                          description="Specifies the k-points used in band structure calculation.", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FILE_NAME", &
                          description="File name used for band structure", &
                          usage="FILE_NAME <filename>", default_c_val="", &
                          type_of_var=char_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ADDED_MOS", &
                          variants=(/"ADDED_BANDS"/), &
                          description="Number of MOs/Bands added to the Band Structure calculation.", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="KPOINT_SET", &
                          description="Specifies a k-point line to be calculated.", &
                          n_keywords=0, n_subsections=0, repeats=.TRUE.)
      ! keywords
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SPECIAL_POINT", &
                          description="Name and coordinates of a special k-point", &
                          usage="SPECIAL_POINT GAMMA 0.0 0.0 0.0", n_var=-1, type_of_var=char_t, repeats=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      !
      CALL keyword_create(keyword, __LOCATION__, name="NPOINTS", &
                          description="Number of k-points along the line.", &
                          default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      !
      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Special k-points are defined either in units"// &
                          " of reciprocal lattice vectors or in Cartesian coordinates in uints of 2Pi/len."// &
                          " B_VECTOR: in multiples of the reciprocal lattice vectors (b)."// &
                          " CART_ANGSTROM: In units of 2*Pi/Angstrom."// &
                          " CART_BOHR: In units of 2*Pi/Bohr.", &
                          usage="UNITS <value>", type_of_var=char_t, default_c_val="B_VECTOR")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      !
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_bandstructure_section

! **************************************************************************************************
!> \brief creates the input section for dealing with homo lumos, including dumping cubes
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_mo_cubes_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_mo_cubes_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_CUBES", &
                                       description="Controls the printing of cubes of the molecular orbitals.", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="write_cube", &
                          description="If the MO cube file should be written. If false, the eigenvalues are still computed."// &
                          " Can also be useful in combination with STM calculations", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="nlumo", &
                          description="If the printkey is activated controls the number of lumos"// &
                          " that are printed and dumped as a cube (-1=all)", &
                          default_i_val=0)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="nhomo", &
         description="If the printkey is activated controls the number of homos that dumped as a cube (-1=all),"// &
         " eigenvalues are always all dumped", &
         default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="homo_list", &
         description="If the printkey is activated controls the index of homos dumped as a cube,"// &
         " eigenvalues are always all dumped. It overrides nhomo.", &
         usage="HOMO_LIST {integer}  {integer} ..  {integer} ", type_of_var=integer_t, &
         n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mo_cubes_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_dos_section(print_key)

      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_dos_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "DOS", &
                                       description="Print Density of States (DOS) (only available states from SCF)", &
                                       print_level=debug_print_level, common_iter_levels=1, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the DOS obtained at different iterations to the output file."// &
                          "By default the file is overwritten", &
                          usage="APPEND", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_E", &
                          description="Histogramm energy spacing.", &
                          usage="DELTA_E 0.0005", type_of_var=real_t, default_r_val=0.001_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_dos_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_pdos_section(print_key)

      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_pdos_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (subsection)
      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PDOS", &
                                       description="Print out the DOS projected per kind and per angular momentum  ", &
                                       print_level=debug_print_level, common_iter_levels=1, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="COMPONENTS", &
                          description="Print out pdos distinguishing all angular momentum components.", &
                          usage="COMPONENTS", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the pdos obtained at different iterations to the pdos  output file."// &
                          "By default the file is overwritten", &
                          usage="APPEND", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="NLUMO", &
         description="Number of virtual orbitals to be added to the MO set (-1=all)."//newline// &
         "CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.", &
         usage="NLUMO integer", default_i_val=0)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OUT_EACH_MO", &
                          description="Output on the status of the calculation every OUT_EACH_MO states. If -1 no output", &
                          usage="OUT_EACH_MO integer", default_i_val=-1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="LDOS", &
                          description="Controls the printing of local PDOS, projected on subsets"// &
                          " of atoms given through lists", &
                          n_keywords=4, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="COMPONENTS", &
                          description="Print out pdos distinguishing all angular momentum components.", &
                          usage="COMPONENTS", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LIST", &
                          description="Specifies a list of indexes of atoms where to project the DOS  ", &
                          usage="LIST {integer}  {integer} ..  {integer} ", type_of_var=integer_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(print_key, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="R_LDOS", &
                          description="Controls the printing of local PDOS, projected on 3D volume in real space,"// &
                          " the volume is defined in terms of position with respect to atoms in the lists", &
                          n_keywords=4, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="LIST", &
                          description="Specifies a list of indexes of atoms used to define the real space volume  ", &
                          usage="LIST {integer}  {integer} ..  {integer} ", type_of_var=integer_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XRANGE", &
                          description="range of positions in Cartesian direction x: all grid points within "// &
                          " this range from at least one atom of the list are considered", &
                          usage="XRANGE -10.0 10.0", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="YRANGE", &
                          description="range of positions in Cartesian direction y: all grid points within "// &
                          " this range from at least one atom of the list are considered", &
                          usage="YRANGE -10.0 10.0", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ZRANGE", &
                          description="range of positions in Cartesian direction z: all grid points within "// &
                          " this range from at least one atom of the list are considered", &
                          usage="ZRANGE -10.0 10.0", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERANGE", &
                          description="only project states with the eigenvalues in the given interval. "// &
                          "Default is all states.", &
                          usage="ERANGE -1.0 1.0", unit_str="hartree", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(print_key, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_pdos_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_wannier_section(print_key)

      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_wannier_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER90", &
                                       description="Interface to Wannier90 code. (EXPERIMENTAL)", &
                                       print_level=debug_print_level, common_iter_levels=1, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="SEED_NAME", &
                          description="The seedname for the Wannier90 calculation (body of filenames).", &
                          usage="SEED_NAME  filename", default_c_val="wannier90", &
                          n_var=1, type_of_var=char_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MP_GRID", &
                          description="The dimensions of the Monkhorst-Pack k-point grid. ", &
                          usage="MP_GRID 6 6 6", n_var=-1, default_i_vals=(/10, 10, 10/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ADDED_MOS", &
                          variants=(/"ADDED_BANDS"/), &
                          description="Number of MOs/Bands added to the Band Structure calculation.", &
                          default_i_val=0)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXCLUDE_BANDS", &
                          description="List of Bands excluded in the Wannier calculation.", &
                          usage="EXCLUDE_BANDS b1 b2 ...", n_var=-1, repeats=.TRUE., &
                          type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WANNIER_FUNCTIONS", &
                          description="Number of Wannier functions to be calculated. ", &
                          usage="WANNIER_FUNCTIONS 6", n_var=1, default_i_val=0, &
                          repeats=.TRUE., type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_wannier_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_stm_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_stm_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "STM", &
                                       description="Controls the printing of cubes for the generation of STM images.", &
                                       print_level=debug_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="nlumo", &
                          description="If the printkey is activated controls the number of additional lumos"// &
                          " that are computed to be able to reproduce STM images obtained"// &
                          "  from positive bias (imaging unoccupied states)", &
                          default_i_val=0)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BIAS", &
                          description="Bias energy for scanning tunneling microscopy (STM) image generation."// &
                          "Orbital densities are summed according to the bias energy."// &
                          "For negative values, states in the range ]EF+bias,EF] are summed,"// &
                          "While positive values sum states in the range [EF,EF+bias[."// &
                          "If positive biases are used, sufficiently many unoccupied stated"// &
                          " (see ADDED_MOS and NLUMO ) should be computed.", &
                          n_var=-1, type_of_var=real_t, default_r_vals=(/0.0_dp/), unit_str='eV')
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TH_TORB", &
                          description="Tip orbital symmetry in Tersoff-Hamann approximation to compute STM images", &
                          repeats=.TRUE., &
                          default_i_val=orb_s, &
                          usage="TH_TORB s dz2", &
                          enum_c_vals=s2a("S", "PX", "PY", "PZ", "DXY", "DYZ", "DZX", "DX2", "DY2", "DZ2"), &
                          enum_i_vals=(/orb_s, orb_px, orb_py, orb_pz, orb_dxy, orb_dyz, orb_dzx, orb_dx2, orb_dy2, orb_dz2/), &
                          enum_desc=s2a("s orbital", "px orbital", "py orbital", "pz orbital", &
                                        "dxy orbital", "dyz orbital", "dzx orbital", "x^2 orbital", "y^2 orbital", "z^2 orbital"))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REF_ENERGY", &
                          description="By default the reference energy is the Fermi energy. In order to compare"// &
                          " with STS experiments, where specific energy ranges are addressed, here"// &
                          "  one can set a different reference energy."// &
                          " The energy range is anyway controlled by the BIAS", &
                          type_of_var=real_t, default_r_val=0.0_dp, unit_str='eV')
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_stm_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_wfn_mix_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_wfn_mix_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (subsection)
      NULLIFY (keyword)

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="WFN_MIX", &
                          description="A section that allows manipulation of the MO coeffs,"// &
                          " e.g. for changing a ground state into an excited state."// &
                          "Starting from a copy of the original MOs, changes can be made"// &
                          "by adding linear combinations of HOMO/LUMO of the original MOs to the result MOs", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="OVERWRITE_MOS", &
                          description="If the keyword is active molecular orbitals in memory will be replaced by the mixed wfn."// &
                          " In combination with RTP or EMD no restart will be required to use the mixed wfn.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="UPDATE", &
                          description="update a result MO with with a linear combination of of original MOs."// &
                          " This section can be repeated to build arbitrary linear combinations using repeatedly y=a*y+b*x.", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_MO_INDEX", &
                          description="Index of the MO (y) to be modified. Counting down in energy with HOMO=1", &
                          usage="RESULT_MO_INDEX 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_MARKED_STATE", &
                          description="Specifies the MO according to "// &
                          "the marks set in MOLECULAR_STATES. The value corresponds to the repetition "// &
                          " of MARK_STATES in MOLECULAR_STATES", &
                          usage="ORIG_MARKED_STATE 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_SPIN_INDEX", &
                          description="Spin of the MO (y) to be modified.", &
                          enum_c_vals=s2a("Alpha", "Beta"), &
                          enum_i_vals=(/1, 2/), & ! direct index in array
                          default_i_val=1, &
                          enum_desc=s2a("Majority spin", "Minority spin"))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_SCALE", &
                          description="Scaling factor of the result variable (a).", &
                          usage="RESULT_SCALE 0.0", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_MO_INDEX", &
                          description="Index of the original MO (x). "// &
                          "Counting down in energy with HOMO=1 or up from LUMO=1, depending on ORIG_IS_VIRTUAL.", &
                          usage="ORIG_MO_INDEX 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_MARKED_STATE", &
                          description="Specifies the MO according to "// &
                          "the marks set in MOLECULAR_STATES. The value corresponds to the repetition "// &
                          " of MARK_STATES in MOLECULAR_STATES", &
                          usage="ORIG_MARKED_STATE 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_SPIN_INDEX", &
                          description="Spin of the MO (x) to be modified.", &
                          enum_c_vals=s2a("Alpha", "Beta"), &
                          enum_i_vals=(/1, 2/), & ! direct index in array
                          default_i_val=1, &
                          enum_desc=s2a("Majority spin", "Minority spin"))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_SCALE", &
                          description="Scaling factor of the result variable (b).", &
                          usage="ORIG_SCALE 0.0", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_IS_VIRTUAL", &
                          description="The original MO (x) is a LUMO.", &
                          usage="ORIG_IS_VIRTUAL", type_of_var=logical_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_wfn_mix_section

! **************************************************************************************************
!> \brief creates the input section for the molecular states
!> \param print_key ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_molecular_states_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_molecular_states_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key2

      CPASSERT(.NOT. ASSOCIATED(print_key))
      NULLIFY (print_key2, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MOLECULAR_STATES", &
                                       description="Controls printing of molecular states ", &
                                       print_level=high_print_level, filename=" ", citations=(/Hunt2003/))

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUBE_EVAL_RANGE", &
         description="only write cubes if the eigenvalues of the corresponding molecular states lie in the given interval. "// &
         "Default is all states.", &
         usage="CUBE_EVAL_RANGE -1.0 1.0", unit_str="hartree", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MARK_STATES", &
                          description="Can be used to mark given molecular states."// &
                          " Sets a mark to both, occupied and unoccupied states. "// &
                          "Occupied states are counted beginning with HOMO=1, "// &
                          "unoccupied states are counted beginning with LUMO=1, "// &
                          "This is only meaningful in combination with WFN_MIX. "// &
                          "First integer specifies the molecule, second integer specifies the state.", &
                          usage="MARK_STATES integer integer", &
                          n_var=2, default_i_vals=(/-1, -1/), type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key2, __LOCATION__, "cubes", &
                                       description="Controls the printing of cube files", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key2, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, print_key2)
      CALL section_release(print_key2)
   END SUBROUTINE create_molecular_states_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_wannier_states_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_wannier_states_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key2

      CPASSERT(.NOT. ASSOCIATED(print_key))
      NULLIFY (print_key2, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_STATES", &
                                       description="Controls printing of molecular states ", &
                                       print_level=high_print_level, filename=" ")

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUBE_EVAL_RANGE", &
         description="only write cubes if the eigenvalues of the corresponding molecular states lie in the given interval. "// &
         "Default is all states.", &
         usage="CUBE_EVAL_RANGE -1.0 1.0", unit_str="hartree", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MARK_STATES", &
                          description="Can be used to mark given molecular states."// &
                          " Sets a mark to both, occupied and unoccupied states. "// &
                          "Occupied states are counted beginning with HOMO=1, "// &
                          "unoccupied states are counted beginning with LUMO=1, "// &
                          "This is only meaningful in combination with WFN_MIX. "// &
                          "First integer specifies the molecule, second integer specifies the state.", &
                          usage="MARK_STATES integer integer", &
                          n_var=2, default_i_vals=(/-1, -1/), type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key2, __LOCATION__, "cubes", &
                                       description="Controls the printing of cube files", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key2, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, print_key2)
      CALL section_release(print_key2)
   END SUBROUTINE create_wannier_states_section

! **************************************************************************************************
!> \brief creates the input section for the qs part
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_qs_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_qs_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="qs", &
                          description="parameters needed to set up the Quickstep framework", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, subsection)

      ! Reals
      CALL keyword_create(keyword, __LOCATION__, name="EPS_DEFAULT", &
                          description="Try setting all EPS_xxx to values leading to an energy correct up to EPS_DEFAULT", &
                          usage="EPS_DEFAULT real", default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CORE_CHARGE", &
                          description="Precision for mapping the core charges.Overrides EPS_DEFAULT/100.0 value", &
                          usage="EPS_CORE_CHARGE real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_GVG_RSPACE", &
         variants=(/"EPS_GVG"/), &
         description="Sets precision of the realspace KS matrix element integration. Overrides SQRT(EPS_DEFAULT) value", &
         usage="EPS_GVG_RSPACE real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_PGF_ORB", &
                          description="Sets precision of the overlap matrix elements. Overrides SQRT(EPS_DEFAULT) value", &
                          usage="EPS_PGF_ORB real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_KG_ORB", &
         description="Sets precision used in coloring the subsets for the Kim-Gordon method. Overrides SQRT(EPS_DEFAULT) value", &
         usage="EPS_KG_ORB 1.0E-8", &
         type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_PPL", &
                          description="Adjusts the precision for the local part of the pseudo potential. ", &
                          usage="EPS_PPL real", type_of_var=real_t, default_r_val=1.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_PPNL", &
         description="Sets precision of the non-local part of the pseudo potential. Overrides sqrt(EPS_DEFAULT) value", &
         usage="EPS_PPNL real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CPC", &
                          description="Sets precision of the GAPW projection. Overrides EPS_DEFAULT value", &
                          usage="EPS_CPC real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RHO", &
                          description="Sets precision of the density mapping on the grids.Overrides EPS_DEFAULT value", &
                          usage="EPS_RHO real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RHO_RSPACE", &
                          description="Sets precision of the density mapping in rspace.Overrides EPS_DEFAULT value."// &
                          ".Overrides EPS_RHO value", &
                          usage="EPS_RHO_RSPACE real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RHO_GSPACE", &
                          description="Sets precision of the density mapping in gspace.Overrides EPS_DEFAULT value."// &
                          ".Overrides EPS_RHO value", &
                          usage="EPS_RHO_GSPACE real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER_MATRIX", &
                          description="Sets the threshold for filtering matrix elements.", &
                          usage="EPS_FILTER_MATRIX 1.0E-6", type_of_var=real_t, default_r_val=0.0E0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPSFIT", &
                          variants=(/"EPS_FIT"/), &
                          description="GAPW: precision to give the extension of a hard gaussian ", &
                          usage="EPSFIT real", default_r_val=1.0E-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPSISO", &
                          variants=(/"EPS_ISO"/), &
                          description="GAPW: precision to determine an isolated projector", &
                          usage="EPSISO real", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPSSVD", &
                          variants=(/"EPS_SVD"/), &
                          description="GAPW: tolerance used in the singular value decomposition of the projector matrix", &
                          usage="EPS_SVD real", default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPSRHO0", &
                          variants=s2a("EPSVRHO0", "EPS_VRHO0"), &
                          description="GAPW : precision to determine the range of V(rho0-rho0soft)", &
                          usage="EPSRHO0 real", default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALPHA0_HARD", &
                          variants=s2a("ALPHA0_H", "ALPHA0"), &
                          description="GAPW: Exponent for hard compensation charge", &
                          usage="ALPHA0_HARD real", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="FORCE_PAW", &
         description="Use the GAPW scheme also for atoms with soft basis sets, i.e. "// &
         " the local densities are computed even if hard and soft should be equal. "// &
         "If this keyword is not set to true, those atoms with soft basis sets are treated by a GPW scheme, i.e. "// &
         "the corresponding density contribution goes on the global grid and is expanded in PW. "// &
         " This option nullifies the effect of the GPW_TYPE in the atomic KIND", &
         usage="FORCE_PAW", &
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_RAD_LOCAL", &
                          description="GAPW : maximum radius of gaussian functions"// &
                          " included in the generation of projectors", &
                          usage="MAX_RAD_LOCAL real", default_r_val=25.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Logicals
      CALL keyword_create(keyword, __LOCATION__, name="LS_SCF", &
                          description="Perform a linear scaling SCF", &
                          usage="LS_SCF", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALMO_SCF", &
                          description="Perform ALMO SCF", &
                          usage="ALMO_SCF", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TRANSPORT", &
                          description="Perform transport calculations (coupling CP2K and OMEN)", &
                          usage="TRANSPORT", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KG_METHOD", &
                          description="Use a Kim-Gordon-like scheme.", &
                          usage="KG_METHOD", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE., citations=(/Iannuzzi2006, Brelaz1979, Andermatt2016/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REF_EMBED_SUBSYS", &
                          description="A total, reference, system in DFT embedding. ", &
                          usage="REF_EMBED_SUBSYS FALSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CLUSTER_EMBED_SUBSYS", &
                          description="A cluster treated with DFT in DFT embedding. ", &
                          usage="CLUSTER_EMBED_SUBSYS FALSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HIGH_LEVEL_EMBED_SUBSYS", &
                          description="A cluster treated with a high-level method in DFT embedding. ", &
                          usage="HIGH_LEVEL_EMBED_SUBSYS FALSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DFET_EMBEDDED", &
                          description="Calculation with DFT-embedding potential. ", &
                          usage="DFET_EMBEDDED FALSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DMFET_EMBEDDED", &
                          description="Calculation with DM embedding potential. ", &
                          usage="DMFET_EMBEDDED FALSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Integers
      CALL keyword_create(keyword, __LOCATION__, name="STO_NG", &
                          description="Order of Gaussian type expansion of Slater orbital basis sets.", &
                          usage="STO_NG", default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LMAXN1", &
                          variants=(/"LMAXRHO1"/), &
                          description="GAPW : max L number for expansion of the atomic densities in spherical gaussians", &
                          usage="LMAXN1 integer", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LMAXN0", &
                          variants=(/"LMAXRHO0"/), &
                          description="GAPW : max L number for the expansion compensation densities in spherical gaussians", &
                          usage="LMAXN0 integer", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LADDN0", &
                          description="GAPW : integer added to the max L of the basis set, used to determine the "// &
                          "maximum value of L for the compensation charge density.", &
                          usage="LADDN0 integer", &
                          default_i_val=99)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Characters
      CALL keyword_create(keyword, __LOCATION__, name="QUADRATURE", &
                          description="GAPW: algorithm to construct the atomic radial grids", &
                          usage="QUADRATURE GC_SIMPLE", &
                          enum_c_vals=s2a("GC_SIMPLE", "GC_TRANSFORMED", "GC_LOG"), &
                          enum_i_vals=(/do_gapw_gcs, do_gapw_gct, do_gapw_log/), &
                          enum_desc=s2a("Gauss-Chebyshev quadrature", &
                                        "Transformed Gauss-Chebyshev quadrature", &
                                        "Logarithmic transformed Gauss-Chebyshev quadrature"), &
                          default_i_val=do_gapw_log)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PW_GRID", &
                          description="What kind of PW_GRID should be employed", &
                          usage="PW_GRID NS-FULLSPACE", &
                          enum_c_vals=s2a("SPHERICAL", "NS-FULLSPACE", "NS-HALFSPACE"), &
                          enum_desc=s2a("- not tested", " tested", " - not tested"), &
                          enum_i_vals=(/do_pwgrid_spherical, do_pwgrid_ns_fullspace, do_pwgrid_ns_halfspace/), &
                          default_i_val=do_pwgrid_ns_fullspace)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PW_GRID_LAYOUT", &
                          description="Force a particular real-space layout for the plane waves grids. "// &
                          "Numbers &le; 0 mean that this dimension is free, incorrect layouts will be ignored. "// &
                          "The default (/-1,-1/) causes CP2K to select a good value, "// &
                          "i.e. plane distributed for large grids, more general distribution for small grids.", &
                          usage="PW_GRID_LAYOUT 4 16", &
                          repeats=.FALSE., n_var=2, &
                          default_i_vals=(/-1, -1/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PW_GRID_BLOCKED", &
                          description="Can be used to set the distribution in g-space for the pw grids and their FFT.", &
                          usage="PW_GRID_BLOCKED FREE", &
                          enum_c_vals=s2a("FREE", "TRUE", "FALSE"), &
                          enum_desc=s2a("CP2K will select an appropriate value", "blocked", "not blocked"), &
                          enum_i_vals=(/do_pw_grid_blocked_free, do_pw_grid_blocked_true, do_pw_grid_blocked_false/), &
                          default_i_val=do_pw_grid_blocked_free)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EXTRAPOLATION", &
         variants=s2a("INTERPOLATION", "WF_INTERPOLATION"), &
         description="Extrapolation strategy for the wavefunction during e.g. MD. "// &
         "Not all options are available for all simulation methods. "// &
         "PS and ASPC are recommended, see also EXTRAPOLATION_ORDER.", &
         citations=(/Kolafa2004, VandeVondele2005a, Kuhne2007/), &
         usage="EXTRAPOLATION PS", &
         enum_c_vals=s2a("USE_GUESS", "USE_PREV_P", "USE_PREV_RHO_R", "LINEAR_WF", &
                         "LINEAR_P", "LINEAR_PS", "USE_PREV_WF", "PS", "FROZEN", "ASPC"), &
         enum_desc=s2a( &
         "Use the method specified with SCF_GUESS, i.e. no extrapolation", &
         "Use the previous density matrix", &
         "Use the previous density in real space", &
         "Linear extrapolation of the wavefunction (not available for K-points)", &
         "Linear extrapolation of the density matrix", &
         "Linear extrapolation of the density matrix times the overlap matrix (not available for K-points)", &
         "Use the previous wavefunction (not available for K-points)", &
         "Higher order extrapolation of the density matrix times the overlap matrix (not available for K-points)", &
         "Frozen ...", &
         "Always stable predictor corrector, similar to PS, but going for MD stability instead of initial guess accuracy. "// &
         "(not available for K-points)"), &
         enum_i_vals=(/ &
         wfi_use_guess_method_nr, &
         wfi_use_prev_p_method_nr, &
         wfi_use_prev_rho_r_method_nr, &
         wfi_linear_wf_method_nr, &
         wfi_linear_p_method_nr, &
         wfi_linear_ps_method_nr, &
         wfi_use_prev_wf_method_nr, &
         wfi_ps_method_nr, &
         wfi_frozen_method_nr, &
         wfi_aspc_nr/), &
         default_i_val=wfi_aspc_nr)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXTRAPOLATION_ORDER", &
                          description="Order for the PS or ASPC extrapolation (typically 2-4). "// &
                          "Higher order might bring more accuracy, but comes, "// &
                          "for large systems, also at some cost. "// &
                          "In some cases, a high order extrapolation is not stable,"// &
                          " and the order needs to be reduced.", &
                          usage="EXTRAPOLATION_ORDER {integer}", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Specifies the electronic structure method that should be employed", &
                          usage="METHOD GAPW", &
                          enum_c_vals=s2a("GAPW", "GAPW_XC", "GPW", "LRIGPW", "RIGPW", &
                                    "MNDO", "MNDOD", "AM1", "PM3", "PM6", "PM6-FM", "PDG", "RM1", "PNNL", "DFTB", "xTB", "OFGPW"), &
                          enum_desc=s2a("Gaussian and augmented plane waves method", &
                                        "Gaussian and augmented plane waves method only for XC", &
                                        "Gaussian and plane waves method", &
                                        "Local resolution of identity method", &
                                        "Resolution of identity method for HXC terms", &
                                        "MNDO semiempirical", "MNDO-d semiempirical", "AM1 semiempirical", &
                                        "PM3 semiempirical", "PM6 semiempirical", "PM6-FM semiempirical", "PDG semiempirical", &
                                        "RM1 semiempirical", &
                                        "PNNL semiempirical", &
                                        "DFTB Density Functional based Tight-Binding", &
                                        "GFN-xTB Extended Tight-Binding", &
                                        "OFGPW Orbital-free GPW method"), &
                          enum_i_vals=(/do_method_gapw, do_method_gapw_xc, do_method_gpw, do_method_lrigpw, do_method_rigpw, &
                                        do_method_mndo, do_method_mndod, do_method_am1, do_method_pm3, &
                                        do_method_pm6, do_method_pm6fm, do_method_pdg, do_method_rm1, &
                                        do_method_pnnl, do_method_dftb, do_method_xtb, do_method_ofgpw/), &
                          citations=(/Lippert1997, Lippert1999, Krack2000, VandeVondele2005a, &
                                      VandeVondele2006, Dewar1977, Dewar1985, Rocha2006, Stewart1989, Thiel1992, &
                                      Repasky2002, Stewart2007, VanVoorhis2015, Schenter2008/), &
                          default_i_val=do_method_gpw)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CORE_PPL", &
                          description="Specifies the method used to calculate the local pseudopotential contribution.", &
                          usage="CORE_PPL ANALYTIC", &
                          enum_c_vals=s2a("ANALYTIC", "GRID"), &
                          enum_desc=s2a("Analytic integration of integrals", &
                                        "Numerical integration on real space grid. Lumped together with core charge"), &
                          enum_i_vals=(/do_ppl_analytic, do_ppl_grid/), &
                          default_i_val=do_ppl_analytic)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_RESTART_FILE_NAME", &
                          description="Root of the file name where to read the embedding "// &
                          "potential guess.", &
                          usage="EMBED_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_CUBE_FILE_NAME", &
                          description="Root of the file name where to read the embedding "// &
                          "potential (guess) as a cube.", &
                          usage="EMBED_CUBE_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_SPIN_CUBE_FILE_NAME", &
                          description="Root of the file name where to read the spin part "// &
                          "of the embedding potential (guess) as a cube.", &
                          usage="EMBED_SPIN_CUBE_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_distribution_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_dftb_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xtb_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mulliken_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ddapc_restraint_section(subsection, "DDAPC_RESTRAINT")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_cdft_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_s2_restraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_lrigpw_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_optimize_lri_basis_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! Embedding subsections: DFET and DMFET
      CALL create_optimize_embed(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_optimize_dmfet(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_qs_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_admm_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_admm_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="AUXILIARY_DENSITY_MATRIX_METHOD", &
                          description="Parameters needed for the ADMM method.", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=(/Guidon2010/))

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ADMM_PURIFICATION_METHOD", &
         description="Method that shall be used for wavefunction fitting. Use MO_DIAG for MD.", &
         enum_c_vals=s2a("NONE", "CAUCHY", "CAUCHY_SUBSPACE", "MO_DIAG", "MO_NO_DIAG", "MCWEENY", "NONE_DM"), &
         enum_i_vals=(/do_admm_purify_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace, &
                       do_admm_purify_mo_diag, do_admm_purify_mo_no_diag, &
                       do_admm_purify_mcweeny, do_admm_purify_none_dm/), &
         enum_desc=s2a("Do not apply any purification", &
                       "Perform purification via general Cauchy representation", &
                       "Perform purification via Cauchy representation in occupied subspace", &
                       "Calculate MO derivatives via Cauchy representation by diagonalization", &
                       "Calculate MO derivatives via Cauchy representation by inversion", &
                       "Perform original McWeeny purification via matrix multiplications", &
                       "Do not apply any purification, works directly with density matrix"), &
         default_i_val=do_admm_purify_mo_diag)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="METHOD", &
         description="Method that shall be used for wavefunction fitting. Use BASIS_PROJECTION for MD.", &
         enum_c_vals=s2a("BASIS_PROJECTION", "BLOCKED_PROJECTION_PURIFY_FULL", "BLOCKED_PROJECTION", &
                         "CHARGE_CONSTRAINED_PROJECTION"), &
         enum_i_vals=(/do_admm_basis_projection, do_admm_blocking_purify_full, do_admm_blocked_projection, &
                       do_admm_charge_constrained_projection/), &
         enum_desc=s2a("Construct auxiliary density matrix from auxiliary basis.", &
                       "Construct auxiliary density from a blocked Fock matrix,"// &
                       " but use the original matrix for purification.", &
                       "Construct auxiliary density from a blocked Fock matrix.", &
                       "Construct auxiliary density from auxiliary basis enforcing charge constrain."), &
         default_i_val=do_admm_basis_projection)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EXCH_SCALING_MODEL", &
         description="Scaling of the exchange correction calculated by the auxiliary density matrix.", &
         enum_c_vals=s2a("NONE", "MERLOT"), &
         enum_i_vals=(/do_admm_exch_scaling_none, do_admm_exch_scaling_merlot/), &
         enum_desc=s2a("No scaling is enabled, refers to methods ADMM1, ADMM2 or ADMMQ.", &
                       "Exchange scaling according to Merlot (2014)"), &
         default_i_val=do_admm_exch_scaling_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EXCH_CORRECTION_FUNC", &
         description="Exchange functional which is used for the ADMM correction.", &
         enum_c_vals=s2a("DEFAULT", "PBEX", "NONE", "OPTX", "BECKE88X"), &
         enum_i_vals=(/do_admm_aux_exch_func_default, do_admm_aux_exch_func_pbex, &
                       do_admm_aux_exch_func_none, do_admm_aux_exch_func_opt, do_admm_aux_exch_func_bee/), &
         enum_desc=s2a("Use PBEX functional for exchange correction, also for truncated and"// &
                       " short Coulomb interaction.", &
                       "Use PBEX functional for exchange correction but always normal form.", &
                       "No correction: X(D)-x(d)-> 0.", &
                       "Use OPTX functional for exchange correction.", &
                       "Use Becke88X functional for exchange correction."), &
         default_i_val=do_admm_aux_exch_func_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="optx_a1", &
                          description="OPTX a1 coefficient", &
                          default_r_val=1.05151_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="optx_a2", &
                          description="OPTX a2 coefficient", &
                          default_r_val=1.43169_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="optx_gamma", &
                          description="OPTX gamma coefficient", &
                          default_r_val=0.006_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BLOCK_LIST", &
                          description="Specifies a list of atoms.", &
                          usage="LIST {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Define accuracy of DBCSR operations", &
                          usage="EPS_FILTER", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_admm_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_se_control_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_se_control_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SE", &
                          description="Parameters needed to set up the Semi-empirical methods", &
                          n_keywords=8, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="ORTHOGONAL_BASIS", &
                          description="Assume orthogonal basis set. This flag is overwritten by "// &
                          "methods with fixed orthogonal/non-orthogonal basis set.", &
                          usage="ORTHOGONAL_BASIS", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STO_NG", &
                          description="Provides the order of the Slater orbital expansion of Gaussian-Type Orbitals.", &
                          usage="STO_NG", default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ANALYTICAL_GRADIENTS", &
                          description="Nuclear Gradients are computed analytically or numerically", &
                          usage="ANALYTICAL_GRADIENTS", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA", &
                          description="Step size in finite difference force calculation", &
                          usage="DELTA {real} ", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTEGRAL_SCREENING", &
                          description="Specifies the functional form for the ", &
                          usage="INTEGRAL_SCREENING (KDSO|KDSO-D|SLATER)", &
                          enum_c_vals=s2a("KDSO", "KDSO-D", "SLATER"), &
                          enum_i_vals=(/do_se_IS_kdso, do_se_IS_kdso_d, do_se_IS_slater/), &
                          enum_desc=s2a("Uses the standard NDDO Klopman-Dewar-Sabelli-Ohno equation "// &
                                        "for the screening of the Coulomb interactions.", &
                                        "Uses a modified Klopman-Dewar-Sabelli-Ohno equation, dumping the screening "// &
                                        "parameter for the Coulomb interactions.", &
                                        "Uses an exponential Slater-type function for modelling the Coulomb interactions."), &
                          default_i_val=do_se_IS_kdso)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC", &
                          description="Specifies the type of treatment for the electrostatic long-range part "// &
                          "in semi-empirical calculations.", &
                          usage="PERIODIC (NONE|EWALD|EWALD_R3|EWALD_GKS)", &
                          enum_c_vals=s2a("NONE", "EWALD", "EWALD_R3", "EWALD_GKS"), &
                          enum_i_vals=(/do_se_lr_none, do_se_lr_ewald, do_se_lr_ewald_r3, do_se_lr_ewald_gks/), &
                          enum_desc=s2a("The long-range part is not explicitly treaten. The interaction "// &
                                        "depends uniquely on the Cutoffs used for the calculation.", &
                                        "Enables the usage of Multipoles Ewald summation schemes. The short-range part "// &
                                        "is tapered according the value of RC_COULOMB.", &
                                        "Enables the usage of Multipoles Ewald summation schemes together with a long-range "// &
                                        "treatment for the 1/R^3 term, deriving from the short-range component. This option "// &
                                        "is active only for K-DSO integral screening type.", &
                                        "Use Ewald directly in Coulomb integral evaluation, works only with Slater screening"), &
                          default_i_val=do_se_lr_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FORCE_KDSO-D_EXCHANGE", &
                          description="This keywords forces the usage of the KDSO-D integral screening "// &
                          "for the Exchange integrals (default is to apply the screening only to the "// &
                          "Coulomb integrals.", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION", &
                          description="Use dispersion correction", &
                          lone_keyword_l_val=.TRUE., &
                          usage="DISPERSION", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_PARAMETER_FILE", &
                          description="Specify file that contains the atomic dispersion parameters", &
                          usage="DISPERSION_PARAMETER_FILE filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_RADIUS", &
                          description="Define radius of dispersion interaction", &
                          usage="DISPERSION_RADIUS", default_r_val=15._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COORDINATION_CUTOFF", &
                          description="Define cutoff for coordination number calculation", &
                          usage="COORDINATION_CUTOFF", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3_SCALING", &
                          description="Scaling parameters (s6,sr6,s8) for the D3 dispersion method,", &
                          usage="D3_SCALING 1.0 1.0 1.0", n_var=3, default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_coulomb_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_exchange_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_screening_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_lr_corr_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_neighbor_lists_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_memory_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_ga_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_se_control_section

! **************************************************************************************************
!> \brief Create the COULOMB se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_lr_corr_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_lr_corr_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="LR_CORRECTION", &
                          description="Setup parameters for the evaluation of the long-range correction term in SE "// &
                          "calculations.", n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="Atomic Cutoff Radius Cutoff for the evaluation of the long-ranbe correction integrals. ", &
                          usage="CUTOFF {real} ", unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(value=6.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering the long-range correction integrals. "// &
                          "If not specified it assumes the same value specified for the CUTOFF.", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2.0*RC_TAPER-20.0*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_lr_corr_section

! **************************************************************************************************
!> \brief Create the COULOMB se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_coulomb_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_coulomb_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="COULOMB", &
                          description="Setup parameters for the evaluation of the COULOMB term in SE "// &
                          "calculations.", n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUTOFF", &
         description="Atomic Cutoff Radius Cutoff for the evaluation of the  Coulomb integrals. "// &
         "For non-periodic calculation the default value is exactly the full cell dimension, in order "// &
         "to evaluate all pair interactions. Instead, for periodic calculations the default numerical value is used.", &
         usage="CUTOFF {real} ", unit_str="angstrom", &
         default_r_val=cp_unit_to_cp2k(value=12.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering Coulomb integrals. "// &
                          "If not specified it assumes the same value specified for the CUTOFF.", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2.0*RC_TAPER-20.0*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_coulomb_section

! **************************************************************************************************
!> \brief Create the EXCHANGE se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_exchange_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_exchange_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXCHANGE", &
                          description="Setup parameters for the evaluation of the EXCHANGE and "// &
                          " core Hamiltonian terms in SE calculations.", n_keywords=0, n_subsections=1, &
                          repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="Atomic Cutoff Radius Cutoff for the evaluation of the Exchange integrals. "// &
                          "For non-periodic calculation the default value is exactly the full cell dimension, in order "// &
                          "to evaluate all pair interactions. Instead, for periodic calculations the default is the "// &
                          "minimum value between 1/4 of the cell dimension and the value specified in input (either"// &
                          " explicitly defined or the default numerical value).", &
                          usage="CUTOFF {real} ", unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(value=12.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering Exchange integrals. "// &
                          "If not specified it assumes the same value specified for the CUTOFF.", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2.0*RC_TAPER-20.0*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_exchange_section

! **************************************************************************************************
!> \brief Create the SCREENING se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_screening_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_screening_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SCREENING", &
                          description="Setup parameters for the tapering of the Coulomb/Exchange Screening in "// &
                          "KDSO-D integral scheme,", n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering the screening term. ", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(value=12.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2*RC_TAPER-20*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_screening_section

! **************************************************************************************************
!> \brief Create the print se section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_se_print_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_se_print_section', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options in SE code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "NEIGHBOR_LISTS", &
                                       description="Activates the printing of the neighbor lists used"// &
                                       " for the periodic SE calculations.", &
                                       print_level=high_print_level, filename="", unit_str="angstrom")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "SUBCELL", &
                                       description="Activates the printing of the subcells used for the"// &
                                       "generation of neighbor lists for periodic SE.", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "EWALD_INFO", &
                                       description="Activates the printing of the information for "// &
                                       "Ewald multipole summation in periodic SE.", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_se_print_section

! **************************************************************************************************
!> \brief creates the input section for use with the GA part of the code
!> \param section the section to create
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE create_se_ga_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_se_ga_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="GA", &
                          description="Sets up memory parameters for the storage of the integrals", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="NCELLS", &
         description="Defines the number of linked cells for the neighbor list. "// &
         "Default value is number of processors", &
         usage="NCELLS 10", &
         default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_se_ga_section

! **************************************************************************************************
!> \brief creates the input section for the se-memory part integral storage
!> \param section the section to create
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE create_se_memory_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_se_memory_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MEMORY", &
                          description="Sets up memory parameters for the storage of the integrals", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EPS_STORAGE", &
         description="Storage threshold for compression is EPS_STORAGE", &
         usage="EPS_STORAGE 1.0E-10", &
         default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MAX_MEMORY", &
         description="Defines the maximum amount of memory [MB] used to store precomputed "// &
         "(possibly compressed) two-electron two-center integrals", &
         usage="MAX_MEMORY 256", &
         default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COMPRESS", &
                          description="Enables the compression of the integrals in memory.", &
                          usage="COMPRESS <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_se_memory_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_mulliken_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_mulliken_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MULLIKEN_RESTRAINT", &
                          description="Use mulliken charges in a restraint (check code for details)", &
                          n_keywords=7, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="STRENGTH", &
                          description="force constant of the restraint", &
                          usage="STRENGTH {real} ", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="target value of the restraint", &
                          usage="TARGET {real} ", default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS", &
                          description="Specifies the list of atoms that is summed in the restraint", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mulliken_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_density_fitting_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_density_fitting_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (keyword, print_key)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DENSITY_FITTING", &
                          description="Setup parameters for density fitting (Bloechl charges or density derived "// &
                          " atomic point charges (DDAPC) charges)", &
                          n_keywords=7, n_subsections=0, repeats=.FALSE., &
                          citations=(/Blochl1995/))

      CALL keyword_create(keyword, __LOCATION__, name="NUM_GAUSS", &
                          description="Specifies the numbers of gaussian used to fit the QM density for each atomic site.", &
                          usage="NUM_GAUSS {integer}", &
                          n_var=1, type_of_var=integer_t, default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PFACTOR", &
                          description="Specifies the progression factor for the gaussian exponent for each atomic site.", &
                          usage="PFACTOR {real}", &
                          n_var=1, type_of_var=real_t, default_r_val=1.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_RADIUS", &
                          description="Specifies the smallest radius of the gaussian used in the fit. All other radius are"// &
                          " obtained with the progression factor.", &
                          usage="MIN_RADIUS {real}", &
                          unit_str="angstrom", n_var=1, type_of_var=real_t, default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RADII", &
                          description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"// &
                          " of this keyword disables all other keywords of this section.", &
                          usage="RADII {real} {real} .. {real}", &
                          unit_str="angstrom", n_var=-1, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GCUT", &
                          description="Cutoff for charge fit in G-space.", &
                          usage="GCUT {real}", &
                          n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing of basic information during the run", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")

      CALL keyword_create(keyword, __LOCATION__, name="CONDITION_NUMBER", &
                          description="Prints information regarding the condition numbers of the A matrix (to be inverted)", &
                          usage="ANALYTICAL_GTERM <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_density_fitting_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \param section_name ...
! **************************************************************************************************
   SUBROUTINE create_ddapc_restraint_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: section_name

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ddapc_restraint_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (keyword, print_key)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name=TRIM(ADJUSTL(section_name)), &
                          description="Use DDAPC charges in a restraint (check code for details)", &
                          n_keywords=7, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE_OF_DENSITY", &
                          description="Specifies the type of density used for the fitting", &
                          usage="TYPE_OF_DENSITY (FULL|SPIN)", &
                          enum_c_vals=s2a("FULL", "SPIN"), &
                          enum_i_vals=(/do_full_density, do_spin_density/), &
                          enum_desc=s2a("Full density", "Spin density"), &
                          default_i_val=do_full_density)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRENGTH", &
                          description="force constant of the restraint", &
                          usage="STRENGTH {real} ", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="target value of the restraint", &
                          usage="TARGET {real} ", default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS", &
                          description="Specifies the list of atoms that is summed in the restraint", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COEFF", &
                          description="Defines the the coefficient of the atom in the atom list (default is one) ", &
                          usage="COEFF 1.0 -1.0", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL_FORM", &
                          description="Specifies the functional form of the term added", &
                          usage="FUNCTIONAL_FORM RESTRAINT", &
                          enum_c_vals=s2a("RESTRAINT", "CONSTRAINT"), &
                          enum_i_vals=(/do_ddapc_restraint, do_ddapc_constraint/), &
                          enum_desc=s2a("Harmonic potential: s*(q-t)**2", "Constraint form: s*(q-t)"), &
                          default_i_val=do_ddapc_restraint)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing basic info about the method", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_ddapc_restraint_section

! **************************************************************************************************
!> \brief Create input section to define CDFT constraint settings specific to Becke weight function.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_becke_constraint_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_becke_constraint_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BECKE_CONSTRAINT", &
                          description="Define settings influencing the construction of the Becke weight function.", &
                          n_keywords=13, repeats=.FALSE., citations=(/Becke1988b/))

      CALL keyword_create(keyword, __LOCATION__, name="ADJUST_SIZE", &
                          description="Adjust Becke cell boundaries with atomic"// &
                          " radii to generate a heteronuclear cutoff profile. These"// &
                          " radii are defined with the keyword ATOMIC_RADII.", &
                          usage="ADJUST_SIZE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_RADII", &
                          description="Defines atomic radii to generate a heteronuclear cutoff profile."// &
                          " Give one value per element in the same order as they"// &
                          " appear in the input coordinates.", &
                          usage="ATOMIC_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SHOULD_SKIP", &
                          description="If grid point is farther than GLOBAL_CUTOFF from all constraint atoms, "// &
                          "move directly to next grid point, thus saving computational resources.", &
                          usage="SHOULD_SKIP", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_CONFINE", &
                          description="Activates Gaussian cavity confinement. The constraint is evaluated only inside "// &
                          "the cavity. The cavity is formed by summing spherical Gaussians centered on the constraint atoms.", &
                          usage="CAVITY_CONFINE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_SHAPE", &
                          description="Specifies the type of Gaussian cavity used.", &
                          usage="CAVITY_SHAPE (SINGLE|VDW|COVALENT|USER)", &
                          enum_c_vals=s2a("DEFAULT", "SINGLE", "VDW", "COVALENT", "USER"), &
                          enum_i_vals=(/radius_default, radius_single, radius_vdw, radius_covalent, radius_user/), &
                          enum_desc=s2a("Use covalent radii (in angstrom) to construct Gaussians, but fixed"// &
                                        " 1.0_dp radius for elements with a radius larger than this value.", &
                                        "Single Gaussian for all atom types with radius given by CAVITY_RADIUS.", &
                                        "Use van der Waals radii to construct Gaussians.", &
                                        "Use covalent radii to construct Gaussians.", &
                                        "Use user defined radii (keyword ATOMIC_RADII) to construct Gaussians."), &
                          default_i_val=radius_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_USE_BOHR", &
                          description="Convert the cavity radius from angstrom to bohr. This results in a larger"// &
                          " confinement cavity than without unit conversion.", &
                          usage="CAVITY_USE_BOHR TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_PRINT", &
                          description="Print cavity in Gaussian cube file format. Currently, printing options"// &
                          " are hardcoded.", &
                          usage="CAVITY_PRINT", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_RADIUS", &
                          description="Radius parameter controlling the creation of Gaussian cavity confinement.", &
                          usage="CAVITY_RADIUS <REAL>", &
                          unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(3.0_dp, "angstrom"), &
                          type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CAVITY", &
                          description="Density threshold for cavity creation. Grid points where the Gaussian"// &
                          " density falls below the threshold are ignored.", &
                          usage="EPS_CAVITY {real} ", default_r_val=1.0e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_TYPE", &
                          description="Specifies the type of cutoff used when building the Becke weight function.", &
                          usage="CUTOFF_TYPE (GLOBAL|ELEMENT)", &
                          enum_c_vals=s2a("GLOBAL", "ELEMENT"), &
                          enum_i_vals=(/becke_cutoff_global, becke_cutoff_element/), &
                          enum_desc=s2a("Use a single value for all elements. Read from GLOBAL_CUTOFF.", &
                                        "Use a different value for all elements. Values read from ELEMENT_CUTOFF."), &
                          default_i_val=becke_cutoff_global)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GLOBAL_CUTOFF", &
                          description="Parameter used to select which atoms contribute to the"// &
                          " weight function at each real space grid point.", &
                          usage="GLOBAL_CUTOFF <REAL>", &
                          unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(3.1750632515_dp, "angstrom"), &
                          type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ELEMENT_CUTOFF", &
                          description="Defines element specific cutoffs to decide which atoms contribute to the"// &
                          " weight function at each real space grid point. Give one value per element in the same"// &
                          " order as they appear in the coordinates.", &
                          usage="ELEMENT_CUTOFF {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IN_MEMORY", &
                          description="Precompute gradients due to Becke constraint during"// &
                          " initial formation of constraint and store them in memory. Useful"// &
                          " in combination with confinement, memory intensive otherwise. Does"// &
                          " nothing if forces are not calculated.", &
                          usage="IN_MEMORY", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_becke_constraint_section

! **************************************************************************************************
!> \brief Creates the input section for defining CDFT constraints.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_cdft_control_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_cdft_control_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: group_section, print_key, subsection

      NULLIFY (keyword, subsection, group_section, print_key)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CDFT", &
                          description="Parameters needed to set up a constrained DFT calculation."// &
                          " Each repetition of the ATOM_GROUP section defines a new constraint."// &
                          " The constraint(s) is (are) converged in a separate external SCF loop with settings"// &
                          " read from the OUTER_SCF section. Supported constraints: Becke and Gaussian"// &
                          " Hirshfeld (partial).", n_keywords=8, n_subsections=2, &
                          repeats=.FALSE., citations=(/Holmberg2017, Holmberg2018/))

      NULLIFY (subsection, keyword)
      CALL create_outer_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_becke_constraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_hirshfeld_constraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE_OF_CONSTRAINT", &
                          description="Specifies the type of constraint used.", &
                          usage="TYPE_OF_CONSTRAINT (NONE|HIRSHFELD|BECKE)", &
                          enum_c_vals=s2a("NONE", "HIRSHFELD", "BECKE"), &
                          enum_i_vals=(/outer_scf_none, outer_scf_hirshfeld_constraint, &
                                        outer_scf_becke_constraint/), &
                          enum_desc=s2a("No constraint (disables section).", &
                                        "Gaussian Hirshfeld constraint. Partial implementation: no forces. "// &
                                        "Requires corresponding section. Not as extensively tested.", &
                                        "Becke constraint. Requires corresponding section."), &
                          citations=(/Becke1988b/), &
                          default_i_val=outer_scf_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRENGTH", &
                          description="Constraint force constants (Lagrange multipliers). "// &
                          "Give one value per constraint group.", &
                          type_of_var=real_t, n_var=-1, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="Constraint target values. Give one value per constraint group. "// &
                          "The target value is the desired number of valence electrons, spin moment, or the number of "// &
                          "alpha or beta electrons on the atoms that define the constraint, suitably multiplied by "// &
                          "atomic coefficients in case a relative constraint between two sets of atoms is employed. "// &
                          "Note that core charges are not subtracted from the target value.", &
                          usage="TARGET {real}", repeats=.FALSE., &
                          type_of_var=real_t, n_var=-1, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_CHARGES", &
                          description="Calculate atomic CDFT charges with selected weight function"// &
                          " (Z = Z_core - Z_CDFT). With fragment based constraints, charges are"// &
                          " relative to the fragment reference state i.e. Z = Z_CDFT - "// &
                          " Z_frag_reference. Note: if the number of atoms is greater than the"// &
                          " default pw_pool max cache, calculation of atomic CDFT charges"// &
                          " will prompt harmless warnings during deallocation of atomic grids.", &
                          usage="ATOMIC_CHARGES", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_A_FILE_NAME", variants=(/"FRAGMENT_A_FILE"/), &
                          description="Name of the reference total electron density cube file for fragment A."// &
                          " May include a path. The reference electron density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_A_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_a.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_B_FILE_NAME", variants=(/"FRAGMENT_B_FILE"/), &
                          description="Name of the reference total electron density cube file for fragment B."// &
                          " May include a path. The reference electron density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_B_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_b.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_A_SPIN_FILE", &
                          variants=(/"FRAGMENT_A_SPIN_FILE_NAME"/), &
                          description="Name of the reference spin density cube file for fragment A."// &
                          " May include a path. The reference spin density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_A_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_a_spin.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_B_SPIN_FILE", &
                          variants=(/"FRAGMENT_B_SPIN_FILE_NAME"/), &
                          description="Name of the reference spin density cube file for fragment B."// &
                          " May include a path. The reference spin density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_B_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_b_spin.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FLIP_FRAGMENT_A", &
                          description="Logical which determines if the reference spin difference density "// &
                          "(rho_alpha-rho_beta) for fragment A should be flipped. With default (off) "// &
                          "value, the fragment is constrained to have more alpha than beta electrons "// &
                          "if the isolated fragment has unpaired electrons. Useful in conjunction with "// &
                          "FLIP_FRAGMENT_B.", &
                          usage="FLIP_FRAGMENT_A", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FLIP_FRAGMENT_B", &
                          description="Logical which determines if the reference spin difference density "// &
                          "(rho_alpha-rho_beta) for fragment B should be flipped. With default (off) "// &
                          "value, the fragment is constrained to have more alpha than beta electrons "// &
                          "if the isolated fragment has unpaired electrons. Useful in conjunction with "// &
                          "FLIP_FRAGMENT_A.", &
                          usage="FLIP_FRAGMENT_B", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of basic info about the method.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")

      CALL section_create(subsection, __LOCATION__, name="WEIGHT_FUNCTION", &
                          description="Controls the printing of cube files with "// &
                          "the CDFT weight function(s). Intended for single-point testing. "// &
                          "In multistep simulations, generated cube files are overwritten each step.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(print_key, subsection)
      CALL section_release(subsection)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL section_create(group_section, __LOCATION__, name="ATOM_GROUP", &
                          description="Define a group of atoms for use in a CDFT constraint. Each repetition of "// &
                          "this section creates a new constraint.", &
                          n_keywords=4, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS", &
                          description="Specifies the list of atoms that are included in the constraint group.", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COEFF", &
                          description="Defines coefficients for the atoms in the list of atoms. Accepts values +/-1.0.", &
                          usage="COEFF 1.0 -1.0", repeats=.TRUE., &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONSTRAINT_TYPE ", &
                          description="Determines what type of constraint to apply. ", &
                          usage="CONSTRAINT_TYPE (CHARGE|MAGNETIZATION|ALPHA|BETA)", &
                          enum_c_vals=s2a("CHARGE", "MAGNETIZATION", "ALPHA", "BETA"), &
                          enum_i_vals=(/cdft_charge_constraint, cdft_magnetization_constraint, &
                                        cdft_alpha_constraint, cdft_beta_constraint/), &
                          enum_desc=s2a("Total charge density constraint (rho_alpha + rho_beta).", &
                                        "Magnetization density constraint (rho_alpha - rho_beta).", &
                                        "Alpha spin density constraint.", &
                                        "Beta spin density constraint."), &
                          default_i_val=cdft_charge_constraint)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_CONSTRAINT", &
                          description="Use a fragment based constraint. "// &
                          "Takes as input the electron densities of two isolated fragments in the "// &
                          "same geometry that they have in the full system. "// &
                          "The isolated fragment densities are read from cube files defined in FRAGMENT_{A,B}_FILE. "// &
                          "For magnetization density constraints, additional files containing the spin difference "// &
                          "densities must be defined with the keywords FRAGMENT_{A,B}_SPIN_FILE. "// &
                          "With this keyword active, the target value of the constraint is calculated from the "// &
                          "the superposition of the isolated fragment densities. Supports only static calculations.", &
                          usage="FRAGMENT_CONSTRAINT", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, group_section)
      CALL section_release(group_section)

      CALL section_create(group_section, __LOCATION__, name="DUMMY_ATOMS", &
                          description="Define an extra group of atoms for which only atomic CDFT charges "// &
                          "should be computed. The section cannot contain any constraint "// &
                          "atoms that were included in section ATOM_GROUP.", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS", &
                          description="Specifies the list of atoms that are included in the DUMMY_ATOMS group.", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, group_section)
      CALL section_release(group_section)

      CALL keyword_create(keyword, __LOCATION__, name="REUSE_PRECOND", &
                          description="Reuse a previously built OT preconditioner between subsequent CDFT SCF iterations "// &
                          "if the inner OT SCF loop converged in PRECOND_FREQ steps or less. Intended mainly for MD "// &
                          "simulations with the FULL_ALL preconditioner to speed up the final iterations of the CDFT SCF loop.", &
                          usage="REUSE_PRECOND yes", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRECOND_FREQ", &
                          description="See REUSE_PRECOND.", &
                          usage="PRECOND_FREQ {int}", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_REUSE", &
                          description="Determines how many times a previously built preconditioner can be reused.", &
                          usage="MAX_REUSE {int}", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PURGE_HISTORY", &
                          description="Purge wavefunction and constraint history to improve SCF convergence during MD."// &
                          " Counts how often the convergence of the first CDFT SCF iteration takes 2 or more outer SCF"// &
                          " iterations and purges the history if the counter exceeds PURGE_FREQ, and PURGE_OFFSET "// &
                          " MD steps have passed since the last purge."// &
                          " The counter is zeroed after each purge.", &
                          usage="PURGE_HISTORY yes", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PURGE_FREQ", &
                          description="See PURGE_HISTORY.", &
                          usage="PURGE_FREQ {int} ", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PURGE_OFFSET", &
                          description="See PURGE_HISTORY.", &
                          usage="PURGE_OFFSET {int} ", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COUNTER", &
                          description="A counter to track the total number of energy evaluations. Needed by"// &
                          " some optimizers to print information. Useful mainly for restarts.", &
                          usage="COUNTER {int} ", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_cdft_control_section

! **************************************************************************************************
!> \brief Creates the input section for defining Gaussian Hirshfeld CDFT constraints.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_hirshfeld_constraint_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_hirshfeld_constraint_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="HIRSHFELD_CONSTRAINT", &
                          description="Parameters for CDFT with a Gaussian Hirshfeld constraint.", &
                          n_keywords=11, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="SHAPE_FUNCTION", &
                          description="Type of shape function used for Hirshfeld partitioning.", &
                          usage="SHAPE_FUNCTION {Gaussian,Density}", repeats=.FALSE., n_var=1, &
                          default_i_val=shape_function_gaussian, &
                          enum_c_vals=s2a("GAUSSIAN", "DENSITY"), &
                          enum_desc=s2a("One Gaussian per atom with radius determined by the keyword GAUSSIAN_SHAPE.", &
                                        "Atomic density expanded in terms of multiple Gaussians."), &
                          enum_i_vals=(/shape_function_gaussian, shape_function_density/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUSSIAN_SHAPE", &
                          description="Specifies the type of Gaussian used for SHAPE_FUNCTION GAUSSIAN.", &
                          usage="GAUSSIAN_SHAPE (SINGLE|VDW|COVALENT|USER)", &
                          enum_c_vals=s2a("DEFAULT", "SINGLE", "VDW", "COVALENT", "USER"), &
                          enum_i_vals=(/radius_default, radius_single, radius_vdw, radius_covalent, radius_user/), &
                          enum_desc=s2a("Use covalent radii (in angstrom) to construct Gaussians, but fixed"// &
                                        " 1.0_dp radius for elements with a radius larger than this value.", &
                                        "Single Gaussian for all atom types with radius given by GAUSSIAN_RADIUS.", &
                                        "Use van der Waals radii to construct Gaussians.", &
                                        "Use covalent radii to construct Gaussians.", &
                                        "Use user defined radii (keyword ATOMIC_RADII) to construct Gaussians."), &
                          default_i_val=radius_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUSSIAN_RADIUS", &
                          description="Radius parameter controlling the creation of Gaussians.", &
                          usage="GAUSSIAN_RADIUS <REAL>", &
                          unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(3.0_dp, "angstrom"), &
                          type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_RADII", &
                          description="Defines custom radii to setup the spherical Gaussians. "// &
                          "Give one value per element in the same order as they "// &
                          "appear in the input coordinates.", &
                          usage="ATOMIC_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_BOHR", &
                          description="Convert the Gaussian radius from angstrom to bohr. This results in a larger "// &
                          "Gaussian than without unit conversion.", &
                          usage="CAVITY_USE_BOHR TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_hirshfeld_constraint_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_s2_restraint_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_s2_restraint_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="S2_RESTRAINT", &
                          description="Use S2 in a re/constraint (OT only)", &
                          n_keywords=7, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="STRENGTH", &
                          description="force constant of the restraint", &
                          usage="STRENGTH {real} ", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="target value of the restraint", &
                          usage="TARGET {real} ", default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL_FORM", &
                          description="Specifies the functional form of the term added", &
                          usage="FUNCTIONAL_FORM RESTRAINT", &
                          enum_c_vals=s2a("RESTRAINT", "CONSTRAINT"), &
                          enum_i_vals=(/do_s2_restraint, do_s2_constraint/), &
                          enum_desc=s2a("Harmonic potential: s*(q-t)**2", "Constraint form: s*(q-t)"), &
                          default_i_val=do_s2_restraint)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_s2_restraint_section

! **************************************************************************************************
!> \brief creates the input section for the tddfpt part
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_tddfpt_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_tddfpt_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="tddfpt", &
                          description="parameters needed to set up the Time Dependent Density Functional PT", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE., &
                          citations=(/Iannuzzi2005/))

      NULLIFY (subsection, keyword)

      ! Integer
      CALL keyword_create(keyword, __LOCATION__, name="MAX_KV", &
                          variants=s2a("MAX_VECTORS"), &
                          description=" maximal number of Krylov space vectors", &
                          usage="MAX_KV someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=60)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTARTS", &
                          variants=s2a("N_RESTARTS"), &
                          description=" maximal number subspace search restarts", &
                          usage="RESTARTS someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NEV", &
                          variants=s2a("N_EV", "EV"), &
                          description=" number of excitations to calculate", &
                          usage="NEV someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NLUMO", &
                          description=" number of additional unoccupied orbitals ", &
                          usage="NLUMO 10", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NREORTHO", &
                          variants=s2a("N_REORTHO", "REORTHO", "REORTHOGONALITAZIONS"), &
                          description=" number of reorthogonalization steps", &
                          usage="NREORTHO someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Logical
      CALL keyword_create(keyword, __LOCATION__, name="KERNEL", &
                          variants=s2a("DO_KERNEL"), &
                          description="compute the kernel (debug purpose only)", &
                          usage="KERNEL logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LSD_SINGLETS", &
                          description="compute singlets using lsd vxc kernel", &
                          usage="LSD_SINGLETS logical_value", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INVERT_S", &
                          variants=s2a("INVERT_OVERLAP"), &
                          description="use the inverse of the overlap matrix", &
                          usage="INVERT_S logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRECONDITIONER", &
                          variants=s2a("PRECOND"), &
                          description="use the preconditioner (only for Davidson)", &
                          usage="PRECONDITIONER logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Character
      CALL keyword_create(keyword, __LOCATION__, name="RES_ETYPE", &
                          variants=s2a("RESTRICTED_EXCITATIONS_TYPE", "RES_E_TYPE"), &
                          description="(singlets/triplets) for restricted calculation", &
                          usage="RES_ETYPE T", &
                          enum_c_vals=s2a("S", "SINGLET", "SINGLETS", "T", "TRIPLET", "TRIPLETS"), &
                          enum_i_vals=(/tddfpt_singlet, tddfpt_singlet, tddfpt_singlet, &
                                        tddfpt_triplet, tddfpt_triplet, tddfpt_triplet/), &
                          default_i_val=tddfpt_singlet)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIAG_METHOD", &
                          variants=s2a("DIAGONALIZATION_METHOD", "METHOD"), &
                          description="Diagonalization method used in tddfpt", &
                          usage="DIAG_METHOD DAVIDSON", &
                          enum_c_vals=s2a("DAVIDSON", "LANCZOS"), &
                          enum_i_vals=(/tddfpt_davidson, tddfpt_lanczos/), &
                          default_i_val=tddfpt_davidson)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OE_CORR", &
                          variants=s2a("ORBITAL_EIGENVALUES_CORRECTION"), &
                          description="Which type of orbital eigenvalue correction to use\n"// &
                          "(to yield better HOMO-LUMO energies)", &
                          usage="OE_CORR SAOP", &
                          enum_c_vals=s2a("NONE", "LB", "LB_ALPHA", "LB94", "GLLB", "GLB", "SAOP", "SIC"), &
                          enum_i_vals=(/oe_none, oe_lb, oe_lb, oe_lb, oe_gllb, oe_gllb, oe_saop, oe_sic/), &
                          default_i_val=oe_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Real
      CALL keyword_create(keyword, __LOCATION__, name="CONVERGENCE", &
                          variants=s2a("CONV"), &
                          description="The convergence of the eigenvalues", &
                          usage="CONVERGENCE 1.0E-6 ", &
                          n_var=1, type_of_var=real_t, &
                          default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_tddfpt_section

! **************************************************************************************************
!> \brief creates the input section for the relativistic part
!> \param section the section to create
!> \author jens
! **************************************************************************************************
   SUBROUTINE create_relativistic_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_relativistic_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="relativistic", &
                          description="parameters needed and setup for relativistic calculations", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="method", &
                          description="type of relativistic correction used", &
                          usage="method (NONE|DKH|ZORA)", default_i_val=rel_none, &
                          enum_c_vals=s2a("NONE", "DKH", "ZORA"), &
                          enum_i_vals=(/rel_none, rel_dkh, rel_zora/), &
                          enum_desc=s2a("Use no relativistic correction", &
                                        "Use Douglas-Kroll-Hess method", &
                                        "Use ZORA method"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DKH_order", &
                          description="The order of the DKH transformation ", &
                          usage="DKH_order 2", default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ZORA_type", &
                          description="Type of ZORA method to be used", &
                          usage="ZORA_type scMP", default_i_val=rel_zora_full, &
                          enum_c_vals=s2a("FULL", "MP", "scMP"), &
                          enum_desc=s2a("Full ZORA method (not implemented)", &
                                        "ZORA with atomic model potential", &
                                        "Scaled ZORA with atomic model potential"), &
                          enum_i_vals=(/rel_zora_full, rel_zora_mp, rel_sczora_mp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="transformation", &
                          description="Type of DKH transformation", &
                          usage="transformation (FULL|MOLECULE|ATOM)", default_i_val=rel_trans_atom, &
                          enum_c_vals=s2a("FULL", "MOLECULE", "ATOM"), &
                          enum_i_vals=(/rel_trans_full, rel_trans_molecule, rel_trans_atom/), &
                          enum_desc=s2a("Use full matrix transformation", &
                                        "Use transformation blocked by molecule", &
                                        "Use atomic blocks"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="z_cutoff", &
                          description="The minimal atomic number considered for atom transformation", &
                          usage="z_cutoff 50", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="potential", &
                          description="External potential used in DKH transformation, full 1/r or erfc(r)/r", &
                          usage="POTENTIAL {FULL,ERFC}", default_i_val=rel_pot_erfc, &
                          enum_c_vals=s2a("FULL", "ERFC"), &
                          enum_i_vals=(/rel_pot_full, rel_pot_erfc/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_relativistic_section

! **************************************************************************************************
!> \brief creates the structure of the section with the DFT SCF parameters
!> \param section will contain the SCF section
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_scf_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_scf_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (print_key)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="scf", &
                          description="Parameters needed to perform an SCF run.", &
                          n_keywords=18, n_subsections=7, repeats=.FALSE.)

      NULLIFY (subsection)

      CALL create_ot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_diagonalization_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_outer_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_smear_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mixing_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mom_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_LUMO", &
                          variants=(/"MAX_ITER_LUMOS"/), &
                          description="The maximum number of iteration for the lumo computation", &
                          usage="MAX_ITER_LUMO 100", default_i_val=299)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LUMO", &
                          variants=(/"EPS_LUMOS"/), &
                          description="target accuracy of the computation of the lumo energy", &
                          usage="EPS_LUMO 1.e-6", default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_SCF", &
                          description="Maximum number of SCF iteration to be performed for one optimization", &
                          usage="MAX_SCF 200", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_SCF_HISTORY", variants=(/"MAX_SCF_HIST"/), &
                          description="Maximum number of SCF iterations after the history pipeline is filled", &
                          usage="MAX_SCF_HISTORY 1", default_i_val=0, lone_keyword_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_DIIS", &
                          variants=(/"MAX_DIIS_BUFFER_SIZE"/), &
                          description="Maximum number of DIIS vectors to be used", &
                          usage="MAX_DIIS 3", default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LEVEL_SHIFT", &
                          variants=(/"LSHIFT"/), &
                          description="Use level shifting to improve convergence", &
                          usage="LEVEL_SHIFT 0.1", default_r_val=0._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCF", &
                          description="Target accuracy for the SCF convergence.", &
                          usage="EPS_SCF 1.e-6", default_r_val=1.e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCF_HISTORY", variants=(/"EPS_SCF_HIST"/), &
                          description="Target accuracy for the SCF convergence after the history pipeline is filled.", &
                          usage="EPS_SCF_HISTORY 1.e-5", default_r_val=0.0_dp, lone_keyword_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHOLESKY", &
                          description="If the cholesky method should be used for computing "// &
                          "the inverse of S, and in this case calling which Lapack routines", &
                          usage="CHOLESKY REDUCE", default_i_val=cholesky_restore, &
                          enum_c_vals=s2a("OFF", "REDUCE", "RESTORE", "INVERSE", "INVERSE_DBCSR"), &
                          enum_desc=s2a("The cholesky algorithm is not used", "Reduce is called", &
                                        "Reduce is replaced by two restore", &
                                        "Restore uses operator multiply by inverse of the triangular matrix", &
                                        "Like inverse, but matrix stored as dbcsr, sparce matrix algebra used when possible"), &
                          enum_i_vals=(/cholesky_off, cholesky_reduce, cholesky_restore, cholesky_inverse, cholesky_dbcsr/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_EIGVAL", &
                          description="Throw away linear combinations of basis functions with a small eigenvalue in S", &
                          usage="EPS_EIGVAL 1.0", default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_DIIS", &
                          description="Threshold on the convergence to start using DIAG/DIIS", &
                          usage="EPS_DIIS 5.0e-2", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="SCF_GUESS", &
         description="Change the initial guess for the wavefunction.", &
         usage="SCF_GUESS RESTART", default_i_val=atomic_guess, &
         enum_c_vals=s2a("ATOMIC", "RESTART", "RANDOM", "CORE", &
                         "HISTORY_RESTART", "MOPAC", "SPARSE", "NONE"), &
         enum_desc=s2a("Generate an atomic density using the atomic code", &
                       "Use the RESTART file as an initial guess (and ATOMIC if not present).", &
                       "Use random wavefunction coefficients.", &
                       "Diagonalize the core hamiltonian for an initial guess.", &
                       "Extrapolated from previous RESTART files.", &
                       "Use same guess as MOPAC for semi-empirical methods or a simple diagonal density matrix for other methods", &
                       "Generate a sparse wavefunction using the atomic code (for OT based methods)", &
                       "Skip initial guess (only for NON-SCC DFTB)."), &
         enum_i_vals=(/atomic_guess, restart_guess, random_guess, core_guess, &
                       history_guess, mopac_guess, sparse_guess, no_guess/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NROW_BLOCK", &
                          description="sets the number of rows in a scalapack block", &
                          usage="NROW_BLOCK 31", default_i_val=32)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NCOL_BLOCK", &
                          description="Sets the number of columns in a scalapack block", &
                          usage="NCOL_BLOCK 31", default_i_val=32)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ADDED_MOS", &
                          description="Number of additional MOS added for each spin", &
                          usage="ADDED_MOS", default_i_val=0, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS_SCHEME", &
                          description="Selects the ROKS scheme when ROKS is applied.", &
                          usage="ROKS_SCHEME HIGH-SPIN", &
                          repeats=.FALSE., &
                          n_var=1, &
                          enum_c_vals=s2a("GENERAL", "HIGH-SPIN"), &
                          enum_i_vals=(/general_roks, high_spin_roks/), &
                          default_i_val=high_spin_roks)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS_F", &
                          variants=(/"F_ROKS"/), &
                          description="Allows to define the parameter f for the "// &
                          "general ROKS scheme.", &
                          usage="ROKS_PARAMETER 1/2", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS_PARAMETERS", &
                          variants=(/"ROKS_PARAMETER"/), &
                          description="Allows to define all parameters for the high-spin "// &
                          "ROKS scheme explicitly. "// &
                          "The full set of 6 parameters has to be specified "// &
                          "in the order acc, bcc, aoo, boo, avv, bvv", &
                          usage="ROKS_PARAMETERS 1/2 1/2 1/2 1/2 1/2 1/2", &
                          repeats=.FALSE., &
                          n_var=6, &
                          type_of_var=real_t, &
                          default_r_vals=(/-0.5_dp, 1.5_dp, 0.5_dp, 0.5_dp, 1.5_dp, -0.5_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NOTCONV_STOPALL", &
                          description="If true, it stops higher level of iteration when SCF does not converge", &
                          usage="NOTCONV_STOPALL logical_value", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="Printing of information during the SCF.", repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of the MO restart file during SCF."// &
                                       "By default keeps a short history of three restarts."// &
                                       "See also RESTART_HISTORY", &
                                       print_level=low_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("QS_SCF"), each_iter_values=(/20/), &
                                       add_last=add_last_numeric, filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum index of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=3)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "RESTART_HISTORY", &
         description="Dumps unique MO restart files during the run keeping all of them.", &
         print_level=low_print_level, common_iter_levels=0, &
         each_iter_names=s2a("__ROOT__", "MD", "GEO_OPT", "ROT_OPT", "NEB", "METADYNAMICS", "QS_SCF"), &
         each_iter_values=(/500, 500, 500, 500, 500, 500, 500/), &
         filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum index of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=3)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "iteration_info", &
                                       description="Controls the printing of basic iteration information during the SCF.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="time_cumul", &
                          description="If the printkey is activated switches the printing of timings"// &
                          " to cumulative (over the SCF).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing of basic information during the SCF.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_ORTHONORMALITY", &
                                       description="Controls the printing relative to the orthonormality of MOs (CT S C).", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_MAGNITUDE", &
                                       description="Prints the min/max eigenvalues of the overlap of the MOs without S (CT C).", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "detailed_energy", &
                                       description="Controls the printing of detailed energy information.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "diis_info", &
                                       description="Controls the printing of DIIS information.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "total_densities", &
                                       description="Controls the printing of total densities.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "Lanczos", &
                                       description="Controls the printing of information on Lanczos refinement iterations.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "DIAG_SUB_SCF", &
         description="Controls the printing of information on subspace diagonalization internal loop. ", &
         print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "Davidson", &
                                       description="Controls the printing of information on Davidson iterations.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FILTER_MATRIX", &
                                       description="Controls the printing of information on Filter Matrix method.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL keyword_create(keyword, __LOCATION__, name="DM_RESTART_WRITE", &
                          description="Write the density matrix into a binary file at the end of the SCF.", &
                          usage="DM_RESTART_WRITE", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_scf_section

! **************************************************************************************************
!> \brief creates the KG section
!> \param section ...
!> \author Martin Haeufel [2012.07]
! **************************************************************************************************
   SUBROUTINE create_kg_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_kg_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="KG_METHOD", &
                          description="Specifies the parameters for a Kim-Gordon-like partitioning"// &
                          " into molecular subunits", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE., &
                          citations=(/Iannuzzi2006, Brelaz1979, Andermatt2016/))

      NULLIFY (keyword, subsection, print_key)

      ! add a XC section
      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! add LRI section
      CALL create_lrigpw_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="COLORING_METHOD", &
                          description="Which algorithm to use for coloring.", &
                          usage="COLORING_METHOD GREEDY", &
                          default_i_val=kg_color_dsatur, &
                          enum_c_vals=s2a("DSATUR", "GREEDY"), &
                          enum_desc=s2a("Maximum degree of saturation, relatively accurate", &
                                        "Greedy, fast coloring, less accurate"), &
                          enum_i_vals=(/kg_color_dsatur, kg_color_greedy/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TNADD_METHOD", &
                          description="Algorithm to use for the calculation of the nonadditive kinetic energy.", &
                          usage="TNADD_METHOD ATOMIC", &
                          default_i_val=kg_tnadd_embed, &
                          enum_c_vals=s2a("EMBEDDING", "RI_EMBEDDING", "ATOMIC", "NONE"), &
                          enum_desc=s2a("Use full embedding potential (see Iannuzzi et al)", &
                                        "Use full embedding potential with RI density fitting", &
                                        "Use sum of atomic model potentials", &
                                        "Do not use kinetic energy embedding"), &
                          enum_i_vals=(/kg_tnadd_embed, kg_tnadd_embed_ri, kg_tnadd_atomic, kg_tnadd_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTEGRATION_GRID", &
                          description="Grid [small,medium,large,huge]to be used for the TNADD integration.", &
                          usage="INTEGRATION_GRID MEDIUM", &
                          default_c_val="MEDIUM")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="Print section", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "NEIGHBOR_LISTS", &
                                       description="Controls the printing of the neighbor lists.", &
                                       print_level=low_print_level, filename="__STD_OUT__", unit_str="angstrom")

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SAB_ORB_FULL", &
                          description="Activates the printing of the full orbital "// &
                          "orbital neighbor lists.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SAB_ORB_MOLECULAR", &
                          description="Activates the printing of the orbital "// &
                          "orbital neighbor lists for molecular subsets.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SAC_KIN", &
                          description="Activates the printing of the orbital "// &
                          "atomic potential neighbor list.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_kg_section

! **************************************************************************************************
!> \brief creates the structure of the section with SCF parameters
!>      controlling an other loop
!> \param section will contain the SCF section
!> \author Joost VandeVondele [2006.03]
! **************************************************************************************************
   SUBROUTINE create_outer_scf_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_outer_scf_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OUTER_SCF", &
                          description="parameters controlling the outer SCF loop", &
                          n_keywords=13, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the outer SCF loop", &
                          usage="&OUTER_SCF ON", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! add CDFT_OPT section
      NULLIFY (subsection)
      CALL create_cdft_opt_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE", &
                          description="Specifies which kind of outer SCF should be employed", &
                          usage="TYPE DDAPC_CONSTRAINT ", &
                          default_i_val=outer_scf_none, &
                          enum_c_vals=s2a("DDAPC_CONSTRAINT", "S2_CONSTRAINT", &
                                          "BASIS_CENTER_OPT", "CDFT_CONSTRAINT", "NONE"), &
                          enum_desc=s2a("Enforce a constraint on the DDAPC, requires the corresponding section", &
                                        "Enforce a constraint on the S2, requires the corresponding section", &
                                        "Optimize positions of basis functions, if atom types FLOATING_BASIS_CENTER "// &
                                        " are defined", &
                                        "Enforce a constraint on a generic CDFT weight population. "// &
                                        "Requires the corresponding section QS&CDFT"// &
                                        " which determines the type of weight used.", &
                                        "Do nothing in the outer loop, useful for resetting the inner loop,"), &
                          enum_i_vals=(/outer_scf_ddapc_constraint, outer_scf_s2_constraint, &
                                        outer_scf_basis_center_opt, outer_scf_cdft_constraint, outer_scf_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPTIMIZER", &
                          description="Method used to bring the outer loop to a stationary point", &
                          usage="OPTIMIZER SD", &
                          default_i_val=outer_scf_optimizer_none, &
                          enum_c_vals=s2a("SD", "DIIS", "NONE", "BISECT", "BROYDEN", "NEWTON", "SECANT", "NEWTON_LS"), &
                          enum_desc=s2a("Takes steps in the direction of the gradient, multiplied by step_size", &
                                        "Uses a Direct Inversion in the Iterative Subspace method", &
                                        "Do nothing, useful only with the none type", &
                                        "Bisection of the gradient, useful for difficult one dimensional cases", &
                                        "Broyden's method. Variant defined in BROYDEN_TYPE.", &
                                        "Newton's method. Only compatible with CDFT constraints.", &
                                        "Secant method. Only for one dimensional cases. See Broyden for "// &
                                        "multidimensional cases.", &
                                        "Newton's method with backtracking line search to find the optimal step size. "// &
                                        "Only compatible with CDFT constraints. Starts from the regular Newton solution "// &
                                        "and successively reduces the step size until the L2 norm of the CDFT gradient "// &
                                        "decreases or MAX_LS steps is reached. Potentially very expensive because "// &
                                        "each iteration performs a full SCF calculation."), &
                          enum_i_vals=(/outer_scf_optimizer_sd, outer_scf_optimizer_diis, outer_scf_optimizer_none, &
                                        outer_scf_optimizer_bisect, outer_scf_optimizer_broyden, &
                                        outer_scf_optimizer_newton, outer_scf_optimizer_secant, &
                                        outer_scf_optimizer_newton_ls/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BISECT_TRUST_COUNT", &
                          description="Maximum number of times the same point will be used in bisection,"// &
                          " a small number guards against the effect of wrongly converged states.", &
                          usage="BISECT_TRUST_COUNT 5", default_i_val=10)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCF", &
                          description="The target gradient of the outer SCF variables. "// &
                          "Notice that the EPS_SCF of the inner loop also determines "// &
                          "the value that can be reached in the outer loop, "// &
                          "typically EPS_SCF of the outer loop must be smaller "// &
                          "than or equal to EPS_SCF of the inner loop.", &
                          usage="EPS_SCF 1.0E-6 ", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIIS_BUFFER_LENGTH", &
                          description="Maximum number of DIIS vectors used ", &
                          usage="DIIS_BUFFER_LENGTH 5", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXTRAPOLATION_ORDER", &
                          description="Number of past states used in the extrapolation of the variables during e.g. MD", &
                          usage="EXTRAPOLATION_ORDER 5", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_SCF", &
                          description="The maximum number of outer loops ", &
                          usage="MAX_SCF 20", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STEP_SIZE", &
                          description="The initial step_size used in the optimizer (currently steepest descent). "// &
                          "Note that in cases where a sadle point is sought for (constrained DFT),"// &
                          " this can be negative. For Newton and Broyden optimizers, use a value less/higher than "// &
                          "the default 1.0 (in absolute value, the sign is not significant) to active an under/overrelaxed "// &
                          "optimizer.", &
                          usage="STEP_SIZE -1.0", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_outer_scf_section

! **************************************************************************************************
!> \brief Create the BSSE section for counterpoise correction
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_bsse_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_bsse_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BSSE", &
                          description="This section is used to set up the BSSE calculation. "// &
                          "It also requires that for each atomic kind X a kind X_ghost is present, "// &
                          "with the GHOST keyword specified, in addition to the other required fields.", &
                          n_keywords=3, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword, subsection)
      ! FRAGMENT SECTION
      CALL section_create(subsection, __LOCATION__, name="FRAGMENT", &
                          description="Specify the atom number belonging to this fragment.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="LIST", &
                          description="Specifies a list of atoms.", &
                          usage="LIST {integer} {integer} .. {integer}", &
                          repeats=.TRUE., n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! CONFIGURATION SECTION
      CALL section_create(subsection, __LOCATION__, name="CONFIGURATION", &
                          description="Specify additional parameters for the combinatorial configurations.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="GLB_CONF", &
                          description="Specifies the global configuration using 1 or 0.", &
                          usage="GLB_CONF {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUB_CONF", &
                          description="Specifies the subconfiguration using 1 or 0 belonging to the global configuration.", &
                          usage="SUB_CONF {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MULTIPLICITY", &
                          variants=(/"MULTIP"/), &
                          description="Specify for each fragment the multiplicity. Two times the total spin plus one. "// &
                          "Specify 3 for a triplet, 4 for a quartet,and so on. Default is 1 (singlet) for an "// &
                          "even number and 2 (doublet) for an odd number of electrons.", &
                          usage="MULTIPLICITY 3", &
                          default_i_val=0) ! this default value is just a flag to get the above
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHARGE", &
                          description="The total charge for each fragment.", &
                          usage="CHARGE -1", &
                          default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="FRAGMENT_ENERGIES", &
                          description="This section contains the energies of the fragments already"// &
                          " computed. It is useful as a summary and specifically for restarting BSSE runs.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="The energy computed for each fragment", repeats=.TRUE., &
                          usage="{REAL}", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_bsse_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_bsse_section

! **************************************************************************************************
!> \brief Create the print bsse section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_print_bsse_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_bsse_section', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options in BSSE code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of information regarding the run.", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of the restart file during BSSE runs."// &
                                       "By default the restart is updated after each configuration calculation is "// &
                                       " completed.", &
                                       print_level=silent_print_level, common_iter_levels=0, &
                                       add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_bsse_section

! **************************************************************************************************
!> \brief creates the interpolation section for the periodic QM/MM
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_gspace_interp_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_gspace_interp_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="interpolator", &
                          description="controls the interpolation for the G-space term", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="aint_precond", &
                          description="the approximate inverse to use to get the starting point"// &
                          " for the linear solver of the spline3 methods", &
                          usage="kind spline3", &
                          default_i_val=precond_spl3_aint, &
                          enum_c_vals=s2a("copy", "spl3_nopbc_aint1", "spl3_nopbc_precond1", &
                                          "spl3_nopbc_aint2", "spl3_nopbc_precond2", "spl3_nopbc_precond3"), &
                          enum_i_vals=(/no_precond, precond_spl3_aint, precond_spl3_1, &
                                        precond_spl3_aint2, precond_spl3_2, precond_spl3_3/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="precond", &
                          description="The preconditioner used"// &
                          " for the linear solver of the spline3 methods", &
                          usage="kind spline3", &
                          default_i_val=precond_spl3_3, &
                          enum_c_vals=s2a("copy", "spl3_nopbc_aint1", "spl3_nopbc_precond1", &
                                          "spl3_nopbc_aint2", "spl3_nopbc_precond2", "spl3_nopbc_precond3"), &
                          enum_i_vals=(/no_precond, precond_spl3_aint, precond_spl3_1, &
                                        precond_spl3_aint2, precond_spl3_2, precond_spl3_3/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_x", &
                          description="accuracy on the solution for spline3 the interpolators", &
                          usage="eps_x 1.e-15", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_r", &
                          description="accuracy on the residual for spline3 the interpolators", &
                          usage="eps_r 1.e-15", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="max_iter", &
                          variants=(/'maxiter'/), &
                          description="the maximum number of iterations", &
                          usage="max_iter 200", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "conv_info", &
                                       description="if convergence information about the linear solver"// &
                                       " of the spline methods should be printed", &
                                       print_level=medium_print_level, each_iter_names=s2a("SPLINE_FIND_COEFFS"), &
                                       each_iter_values=(/10/), filename="__STD_OUT__", &
                                       add_last=add_last_numeric)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_gspace_interp_section

! **************************************************************************************************
!> \brief input section for optional parameters for LRIGPW
!>        LRI: local resolution of identity
!> \param section the section to create
!> \author Dorothea Golze [02.2015]
! **************************************************************************************************
   SUBROUTINE create_lrigpw_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_lrigpw_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="LRIGPW", &
                          description="This section specifies optional parameters for LRIGPW.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE., citations=(/Golze2017b/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LRI_OVERLAP_MATRIX", &
                          description="Specifies whether to calculate the inverse or the "// &
                          "pseudoinverse of the overlap matrix of the auxiliary "// &
                          "basis set. Calculating the pseudoinverse is necessary "// &
                          "for very large auxiliary basis sets, but more expensive. "// &
                          "Using the pseudoinverse, consistent forces are not "// &
                          "guaranteed yet.", &
                          usage="LRI_OVERLAP_MATRIX INVERSE", &
                          enum_c_vals=s2a("INVERSE", "PSEUDO_INVERSE_SVD", "PSEUDO_INVERSE_DIAG", &
                                          "AUTOSELECT"), &
                          enum_desc=s2a("Calculate inverse of the overlap matrix.", &
                                        "Calculate the pseuodinverse of the overlap matrix "// &
                                        "using singular value decomposition.", &
                                        "Calculate the pseudoinverse of the overlap matrix "// &
                                        "by prior diagonalization.", &
                                        "Choose automatically for each pair whether to "// &
                                        "calculate the inverse or pseudoinverse based on the "// &
                                        "condition number of the overlap matrix for each pair. "// &
                                        "Calculating the pseudoinverse is much more expensive."), &
                          enum_i_vals=(/do_lri_inv, do_lri_pseudoinv_svd, &
                                        do_lri_pseudoinv_diag, do_lri_inv_auto/), &
                          default_i_val=do_lri_inv)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_CONDITION_NUM", &
                          description="If AUTOSELECT is chosen for LRI_OVERLAP_MATRIX, this "// &
                          "keyword specifies that the pseudoinverse is calculated "// &
                          "only if the LOG of the condition number of the lri "// &
                          "overlap matrix is larger than the given value.", &
                          usage="MAX_CONDITION_NUM 20.0", default_r_val=20.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_O3_INT", &
                          description="Threshold for ABA and ABB integrals in LRI. "// &
                          "This is used for screening in the KS and "// &
                          "force calculations (tensor contractions).", &
                          usage="EPS_O3_INT 1.e-10", default_r_val=1.0e-14_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DEBUG_LRI_INTEGRALS", &
                          description="Debug the integrals needed for LRIGPW.", &
                          usage="DEBUG_LRI_INTEGRALS TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXACT_1C_TERMS", &
                          description="Don't use LRI for one center densities.", &
                          usage="EXACT_1C_TERMS TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PPL_RI", &
                          description="Use LRI/RI for local pseudopotential.", &
                          usage="PPL_RI TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RI_STATISTIC", &
                          description="Print statistical information on the RI calculation.", &
                          usage="RI_STATISTIC TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISTANT_PAIR_APPROXIMATION", &
                          description="Calculate distant pairs using an independent atom approximation.", &
                          usage="DISTANT_PAIR_APPROXIMATION TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISTANT_PAIR_METHOD", &
                          description="Method used to separate pair density for distant pairs."// &
                          "Options: EW (equal weights); AW (atomic weights); SW (set weights); "// &
                          "LW (shell function weights)", &
                          usage="DISTANT_PAIR_METHOD {method}", &
                          default_c_val="LW")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISTANT_PAIR_RADII", &
                          description="Inner and outer radii used in distant "// &
                          "pair separation. Smooth interpolation between inner and outer "// &
                          "radius is used.", &
                          usage="DISTANT_PAIR_RADII r_inner {real} r_outer {real} ", &
                          n_var=2, default_r_vals=(/8._dp, 12._dp/), unit_str='bohr', &
                          type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SHG_LRI_INTEGRALS", &
                          description="Uses the SHG (solid harmonic Gaussian) integral "// &
                          "scheme instead of Obara-Saika", &
                          usage="SHG_LRI_INTEGRALS TRUE", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., &
                          citations=(/Golze2017a/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RI_SINV", &
                          description="Approximation to be used for the inverse of the "// &
                          "RI overlap matrix. INVF, INVS: exact inverse, apply directly "// &
                          "for solver (F:full matrix, S:sparsematrix). AINV approximate inverse, use with PCG. "// &
                          "NONE: no approximation used with CG solver.", &
                          usage="RI_SINV NONE", default_c_val="INVF")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_lrigpw_section

! **************************************************************************************************
!> \brief input section for optimization of the auxililary basis for LRIGPW
!> \param section the section to create
!> \author Dorothea Golze [05.2014]
! **************************************************************************************************
   SUBROUTINE create_optimize_lri_basis_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_optimize_lri_basis_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPTIMIZE_LRI_BASIS", &
                          description="This section specifies the parameters for optimizing "// &
                          "the lri auxiliary basis sets for LRIGPW. The Powell optimizer is used.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="ACCURACY", &
                          description="Target accuracy for the objective function (RHOEND)", &
                          usage="ACCURACY 5.0E-4", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_FUN", &
                          description="Maximum number of function evaluations", &
                          usage="MAX_FUN 200", default_i_val=4000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STEP_SIZE", &
                          description="Initial step size for search algorithm (RHOBEG)", &
                          usage="STEP_SIZE 1.0E-1", default_r_val=5.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONDITION_WEIGHT", &
                          description="This keyword allows to give different weight "// &
                          "factors to the condition number (LOG(cond) is used).", &
                          usage="CONDITION_WEIGHT 1.0E-4", default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_CONDITION_NUMBER", &
                          description="Determines whether condition number should be part "// &
                          "of optimization or not", &
                          usage="USE_CONDITION_NUMBER", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GEOMETRIC_SEQUENCE", &
                          description="Exponents are assumed to be a geometric sequence. "// &
                          "Only the minimal and maximal exponents of one set are optimized and "// &
                          "the other exponents are obtained by geometric progression.", &
                          usage="GEOMETRIC_SEQUENCE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DEGREES_OF_FREEDOM", &
                          description="Specifies the degrees of freedom in the basis "// &
                          "optimization.", &
                          usage="DEGREES_OF_FREEDOM ALL", &
                          enum_c_vals=s2a("ALL", "COEFFICIENTS", "EXPONENTS"), &
                          enum_desc=s2a("Set all parameters in the basis to be variable.", &
                                        "Set all coefficients in the basis set to be variable.", &
                                        "Set all exponents in the basis to be variable."), &
                          enum_i_vals=(/do_lri_opt_all, do_lri_opt_coeff, do_lri_opt_exps/), &
                          default_i_val=do_lri_opt_exps)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_constrain_exponents_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_optimize_lri_basis_section

! **************************************************************************************************
!> \brief Input for DFT embedding
!> \param section ...
!> \author Vladimir Rybkin [08.2017]
! **************************************************************************************************
   SUBROUTINE create_optimize_embed(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_optimize_embed', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPT_EMBED", &
                          description="This section specifies optional parameters for DFT embedding potential optimization.", &
                          n_keywords=19, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REG_LAMBDA", &
                          description="Parameter for Yang's regularization "// &
                          "involving kinetic matrix.", &
                          usage="REG_LAMBDA 0.0001", default_r_val=0.0001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_ITER", &
                          description="Maximum number of iterations "// &
                          "in the optimization procedure.", &
                          usage="N_ITER 75", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TRUST_RAD", &
                          description="Maximum number of iterations "// &
                          "in the optimization procedure.", &
                          usage="TRUST_RAD 0.5", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENS_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum electron density difference.", &
                          usage="DENS_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENS_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the integrated electron density difference.", &
                          usage="DENS_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_DENS_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum electron density difference.", &
                          usage="DENS_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_DENS_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the integrated electron density difference.", &
                          usage="DENS_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPTIMIZER", &
                          description="Optimize embedding potential.", &
                          usage="OPTIMIZER  LEVEL_SHIFT", &
                          default_i_val=embed_steep_desc, &
                          enum_c_vals=s2a("STEEPEST_DESCENT", "QUASI_NEWTON", "LEVEL_SHIFT"), &
                          enum_desc=s2a("Steepest descent.", "Quasi-Newton.", "Level shift."), &
                          enum_i_vals=(/embed_steep_desc, embed_quasi_newton, embed_level_shift/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID_OPT", &
                          description="Optimize embedding potential on the grid. ", &
                          usage="GRID_OPT .TRUE.", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LEEUWEN-BAERENDS", &
                          description="Van Leeuwen-Baerends iterative update. Alternative to Wu-Yang "// &
                          "optimizer. Use only with ADD_CONTST_POT.", &
                          usage="LEEUWEN-BAERENDS .TRUE.", &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FAB", &
                          description="Finzel-Ayers-Bultinck iterative update. Generally, not reliable. ", &
                          usage="FAB .TRUE.", &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VW_CUTOFF", &
                          description="Cutoff for von Weizsacker potential in "// &
                          "the FAB optimization procedure.", &
                          usage="VW_CUTOFF 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VW_SMOOTH_CUT_RANGE", &
                          description="Smooth cutoff range for von Weizsacker potential in "// &
                          "the FAB optimization procedure.", &
                          usage="VW_CUTOFF 1.0", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POT_GUESS", &
                          description="Specifies the guess of the embedding  "// &
                          "potential. For optimization in finite basis (not grid optimization) "// &
                          "in is a constant part to be added to the one in finite basis. ", &
                          usage="DEGREES_OF_FREEDOM ALL", &
                          enum_c_vals=s2a("NONE", "DIFF", "Fermi_Amaldi", "RESP"), &
                          enum_desc=s2a("Initial guess is zero grid.", &
                                        "Initial density difference. A euristic but working approach.", &
                                        "Fermi-Amaldi potential. More rigorous than DIFF, although less efficient.", &
                                        "Coulomb interaction between the subsystem using RESP charges)"// &
                                        " on the total system."), &
                          enum_i_vals=(/embed_none, embed_diff, embed_fa, embed_resp/), &
                          default_i_val=embed_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHARGE_DISTR_WIDTH", &
                          description="Width of the Gaussian representing "// &
                          "point charges. To be used with ADD_COULOMB_POT.", &
                          usage="CHARGE_DISTR_WIDTH 3.000", default_r_val=1.12490_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="READ_EMBED_POT", &
                          description="Read the embedding potential  "// &
                          "restart vector as a guess.", &
                          usage="READ_EMBED_POT .FALSE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="READ_EMBED_POT_CUBE", &
                          description="Read the embedding potential  "// &
                          "(restart) from the cube file.", &
                          usage="READ_EMBED_POT_CUBE .FALSE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_RESTART_FILE_NAME", &
                          description="Root of the file name where to read the embedding "// &
                          "potential guess.", &
                          usage="EMBED_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_CUBE_FILE_NAME", &
                          description="Root of the file name where to read the embedding "// &
                          "potential (guess) as a cube.", &
                          usage="EMBED_CUBE_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_SPIN_CUBE_FILE_NAME", &
                          description="Root of the file name where to read the spin part "// &
                          "of the embedding potential (guess) as a cube.", &
                          usage="EMBED_SPIN_CUBE_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_print_embed_diff(section)

      CALL create_print_embed_pot_cube(section)

      CALL create_print_embed_restart_vec(section)

      CALL create_print_simple_grid(section)

   END SUBROUTINE create_optimize_embed

! **************************************************************************************************
!> \brief Input for density matrix functional embedding, DMFET
!> \param section ...
!> \author Vladimir Rybkin [08.2018]
! **************************************************************************************************
   SUBROUTINE create_optimize_dmfet(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_optimize_dmfet', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPT_DMFET", &
                          description="This section specifies optional parameters for DMFET matrix potential optimization.", &
                          n_keywords=8, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_ITER", &
                          description="Maximum number of iterations "// &
                          "in the optimization procedure.", &
                          usage="N_ITER 75", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TRUST_RAD", &
                          description="Step length "// &
                          "in the optimization procedure.", &
                          usage="TRUST_RAD 0.5", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DM_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum element of density matrix difference.", &
                          usage="DM_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DM_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the total density matrix difference.", &
                          usage="DM_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BETA_DM_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum element of the beta-spin density "// &
                          "matrix difference.", &
                          usage="DM_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BETA_DM_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the total beta-spin density matrix difference.", &
                          usage="DM_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="READ_DMFET_POT", &
                          description="Read the matrix embedding potential  "// &
                          "(restart) from the cube file.", &
                          usage="READ_DMFET_POT .FALSE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DMFET_RESTART_FILE_NAME", &
                          description="Root of the file name where to read the matrix "// &
                          "potential guess.", &
                          usage="DMFET_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_optimize_dmfet

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_embed_diff(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_embed_diff', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_DENS_DIFF", &
                                       description="Controls the printing of cube files with "// &
                                       " embedding densisty differences", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_diff

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_embed_pot_cube(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_embed_pot_cube', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_POT_CUBE", &
                                       description="Controls the printing of cube files with "// &
                                       "with embedding potential", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 1 1 1", n_var=-1, default_i_vals=(/1, 1, 1/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_pot_cube

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_embed_restart_vec(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_embed_restart_vec', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_POT_VECTOR", &
                                       description="Controls the printing of cube files with "// &
                                       "with embedding potential", &
                                       print_level=silent_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_restart_vec

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_simple_grid(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_simple_grid', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WRITE_SIMPLE_GRID", &
                                       description="Controls the printing of simple grid "// &
                                       "files with embedding potential: X Y Z value", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 1 1 1", n_var=-1, default_i_vals=(/1, 1, 1/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Units of the volumetric file: Angstrom or Bohr.", &
                          usage="UNITS BOHR", &
                          default_i_val=embed_grid_bohr, &
                          enum_c_vals=s2a("BOHR", "ANGSTROM"), &
                          enum_desc=s2a("Atomic units: Bohr", "Metric units: Angstrom."), &
                          enum_i_vals=(/embed_grid_bohr, embed_grid_angstrom/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FOLD_COORD", &
                          description="Activates printing folded coordinates corresponding "// &
                          "to the simple grid. Used as input for external programs.", &
                          usage="FOLD_COORD .TRUE.", n_var=1, type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)

      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_simple_grid

! **************************************************************************************************
!> \brief input section for optional parameters for RIGPW
!> \param section the section to create
!> \author JGH [06.2017]
! **************************************************************************************************
   SUBROUTINE create_rigpw_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_rigpw_section', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RIGPW", &
                          description="This section specifies optional parameters for RIGPW.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

!     CALL keyword_create(keyword, __LOCATION__, name="RI_OVERLAP_MATRIX", &
!                         description="Specifies whether to calculate the inverse or the "// &
!                         "pseudoinverse of the overlap matrix of the auxiliary "// &
!                         "basis set. Calculating the pseudoinverse is necessary "// &
!                         "for very large auxiliary basis sets, but more expensive. "// &
!                         "Using the pseudoinverse, consistent forces are not "// &
!                         "guaranteed yet.", &
!                         usage="RI_OVERLAP_MATRIX INVERSE", &
!                         enum_c_vals=s2a("INVERSE", "PSEUDO_INVERSE_SVD", "PSEUDO_INVERSE_DIAG", &
!                                         "AUTOSELECT"), &
!                         enum_desc=s2a("Calculate inverse of the overlap matrix.", &
!                                       "Calculate the pseuodinverse of the overlap matrix "// &
!                                       "using singular value decomposition.", &
!                                       "Calculate the pseudoinverse of the overlap matrix "// &
!                                       "by prior diagonalization.", &
!                                       "Choose automatically for each pair whether to "// &
!                                       "calculate the inverse or pseudoinverse based on the "// &
!                                       "condition number of the overlap matrix for each pair. "// &
!                                       "Calculating the pseudoinverse is much more expensive."), &
!                         enum_i_vals=(/do_lri_inv, do_lri_pseudoinv_svd, &
!                                       do_lri_pseudoinv_diag, do_lri_inv_auto/), &
!                         default_i_val=do_lri_inv)
!     CALL section_add_keyword(section, keyword)
!     CALL keyword_release(keyword)

   END SUBROUTINE create_rigpw_section

! **************************************************************************************************
!> \brief input section for constraints for auxiliary basis set optimization
!> \param section the section to create
!> \author Dorothea Golze [11.2014]
! **************************************************************************************************
   SUBROUTINE create_constrain_exponents_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_constrain_exponents_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CALL section_create(section, __LOCATION__, name="CONSTRAIN_EXPONENTS", &
                          description="specifies constraints for the exponents of the "// &
                          "lri auxiliary basis sets in the optimization.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          description="Defines the upper and lower boundaries as "// &
                          "(1+scale)*exp and (1-scale)*exp. Fermi-like constraint "// &
                          "function", &
                          usage="SCALE 0.3", default_r_val=0.3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FERMI_EXP", &
                          description="Exponent in the fermi-like constraint function. ", &
                          usage="FERMI_EXP 2.63", default_r_val=2.63391_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_constrain_exponents_section

! **************************************************************************************************
!> \brief creates the multigrid
!> \param section             input section to create
!> \param create_subsections  indicates whether or not subsections INTERPOLATOR and RS_GRID
!>                            should be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_mgrid_section(section, create_subsections)
      TYPE(section_type), POINTER                        :: section
      LOGICAL, INTENT(in)                                :: create_subsections

      CHARACTER(len=*), PARAMETER :: routineN = 'create_mgrid_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="mgrid", &
                          description="multigrid information", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NGRIDS", &
                          description="The number of multigrids to use", &
                          usage="ngrids 1", default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="cutoff", &
                          description="The cutoff of the finest grid level. Default value for "// &
                          "SE or DFTB calculation is 1.0 [Ry].", &
                          usage="cutoff 300", default_r_val=cp_unit_to_cp2k(value=280.0_dp, &
                                                                            unit_str="Ry"), n_var=1, unit_str="Ry")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="progression_factor", &
                          description="Factor used to find the cutoff of the multigrids that"// &
                          " where not given explicitly", &
                          usage="progression_factor <integer>", default_r_val=3._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="commensurate", &
                          description="If the grids should be commensurate. If true overrides "// &
                          "the progression factor and the cutoffs of the sub grids", &
                          usage="commensurate", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="realspace", &
                          description="If both rho and rho_gspace are needed ", &
                          usage="realspace", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REL_CUTOFF", &
                          variants=(/"RELATIVE_CUTOFF"/), &
                          description="Determines the grid at which a Gaussian is mapped,"// &
                          " giving the cutoff used for a gaussian with alpha=1."// &
                          " A value 50+-10Ry might be required for highly accurate results, "// &
                          " Or for simulations with a variable cell."// &
                          " Versions prior to 2.3 used a default of 30Ry.", &
                          usage="RELATIVE_CUTOFF real", default_r_val=20.0_dp, &
                          unit_str="Ry")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIGRID_SET", &
                          description="Activate a manual setting of the multigrids", &
                          usage="MULTIGRID_SET", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SKIP_LOAD_BALANCE_DISTRIBUTED", &
                          description="Skips load balancing on distributed multigrids.  "// &
                          "Memory usage is O(p) so may be used "// &
                          "for all but the very largest runs.", &
                          usage="SKIP_LOAD_BALANCE_DISTRIBUTED", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
!    CALL keyword_create(keyword, __LOCATION__, name="SKIP_LOAD_BALANCE_DISTRIBUTED",&
!         description="Skip load balancing on distributed multigrids, which might be memory intensive."//&
!                     "If not explicitly specified, runs using more than 1024 MPI tasks will default to .TRUE.",&
!         usage="SKIP_LOAD_BALANCE_DISTRIBUTED", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIGRID_CUTOFF", &
                          variants=(/"CUTOFF_LIST"/), &
                          description="List of cutoff values to set up multigrids manually", &
                          usage="MULTIGRID_CUTOFF 200.0 100.0 ", &
                          n_var=-1, &
                          type_of_var=real_t, &
                          unit_str="Ry")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      IF (create_subsections) THEN
         NULLIFY (subsection)
         CALL create_rsgrid_section(subsection)
         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)

         NULLIFY (subsection)
         CALL create_interp_section(subsection)
         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)
      END IF
   END SUBROUTINE create_mgrid_section

! **************************************************************************************************
!> \brief creates the interpolation section
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_interp_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_interp_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="interpolator", &
                          description="kind of interpolation used between the multigrids", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="kind", &
                          description="the interpolator to use", &
                          usage="kind spline3", &
                          default_i_val=pw_interp, &
                          enum_c_vals=s2a("pw", "spline3_nopbc", "spline3"), &
                          enum_i_vals=(/pw_interp, &
                                        spline3_nopbc_interp, spline3_pbc_interp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="safe_computation", &
                          description="if a non unrolled calculation is to be performed in parallel", &
                          usage="safe_computation OFF", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="aint_precond", &
                          description="the approximate inverse to use to get the starting point"// &
                          " for the linear solver of the spline3 methods", &
                          usage="aint_precond copy", &
                          default_i_val=precond_spl3_aint, &
                          enum_c_vals=s2a("copy", "spl3_nopbc_aint1", "spl3_nopbc_aint2", &
                                          "spl3_nopbc_precond1", "spl3_nopbc_precond2", "spl3_nopbc_precond3"), &
                          enum_i_vals=(/no_precond, precond_spl3_aint, precond_spl3_aint2, &
                                        precond_spl3_1, precond_spl3_2, precond_spl3_3/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="precond", &
                          description="The preconditioner used"// &
                          " for the linear solver of the spline3 methods", &
                          usage="PRECOND copy", &
                          default_i_val=precond_spl3_3, &
                          enum_c_vals=s2a("copy", "spl3_nopbc_aint1", "spl3_nopbc_aint2", &
                                          "spl3_nopbc_precond1", "spl3_nopbc_precond2", "spl3_nopbc_precond3"), &
                          enum_i_vals=(/no_precond, precond_spl3_aint, precond_spl3_aint2, &
                                        precond_spl3_1, precond_spl3_2, precond_spl3_3/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_x", &
                          description="accuracy on the solution for spline3 the interpolators", &
                          usage="eps_x 1.e-15", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_r", &
                          description="accuracy on the residual for spline3 the interpolators", &
                          usage="eps_r 1.e-15", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="max_iter", &
                          variants=(/'maxiter'/), &
                          description="the maximum number of iterations", &
                          usage="max_iter 200", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "conv_info", &
                                       description="if convergence information about the linear solver"// &
                                       " of the spline methods should be printed", &
                                       print_level=medium_print_level, each_iter_names=s2a("SPLINE_FIND_COEFFS"), &
                                       each_iter_values=(/10/), filename="__STD_OUT__", &
                                       add_last=add_last_numeric)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_interp_section

! **************************************************************************************************
!> \brief creates the sic (self interaction correction) section
!> \param section ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_sic_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_sic_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="sic", &
                          description="parameters for the self interaction correction", &
                          n_keywords=6, n_subsections=0, repeats=.FALSE., &
                          citations=(/VandeVondele2005b, Perdew1981, Avezac2005/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIC_SCALING_A", &
                          description="Scaling of the coulomb term in sic [experimental]", &
                          usage="SIC_SCALING_A 0.5", &
                          citations=(/VandeVondele2005b/), &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIC_SCALING_B", &
                          description="Scaling of the xc term in sic [experimental]", &
                          usage="SIC_SCALING_B 0.5", &
                          citations=(/VandeVondele2005b/), &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIC_METHOD", &
                          description="Method used to remove the self interaction", &
                          usage="SIC_METHOD MAURI_US", &
                          default_i_val=sic_none, &
                          enum_c_vals=s2a("NONE", "MAURI_US", "MAURI_SPZ", "AD", "EXPLICIT_ORBITALS"), &
                          enum_i_vals=(/sic_none, sic_mauri_us, sic_mauri_spz, sic_ad, sic_eo/), &
                          enum_desc=s2a("Do not apply a sic correction", &
                                        "Employ a (scaled) correction proposed by Mauri and co-workers"// &
                                        " on the spin density / doublet unpaired orbital", &
                                        "Employ a (scaled) Perdew-Zunger expression"// &
                                        " on the spin density / doublet unpaired orbital", &
                                        "The average density correction", &
                                        "(scaled) Perdew-Zunger correction explicitly on a set of orbitals."), &
                          citations=(/VandeVondele2005b, Perdew1981, Avezac2005/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_SET", &
                          description="Type of orbitals treated with the SIC", &
                          usage="ORBITAL_SET ALL", &
                          default_i_val=sic_list_unpaired, &
                          enum_c_vals=s2a("UNPAIRED", "ALL"), &
                          enum_desc=s2a("correction for the unpaired orbitals only, requires a restricted open shell calculation", &
                                        "correction for all orbitals, requires a LSD or ROKS calculation"), &
                          enum_i_vals=(/sic_list_unpaired, sic_list_all/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_sic_section

! **************************************************************************************************
!> \brief creates the low spin roks section
!> \param section ...
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE create_low_spin_roks_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_low_spin_roks_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="LOW_SPIN_ROKS", &
                          description="Specify the details of the low spin ROKS method. "// &
                          "In particular, one can specify various terms added to the energy of the high spin roks configuration"// &
                          " with a energy scaling factor, and a prescription of the spin state.", &
                          n_keywords=6, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_SCALING", &
                          description="The scaling factors for each term added to the total energy. "// &
                          "This list should contain one number for each term added to the total energy.", &
                          usage="ENERGY_SCALING 1.0 -1.0 ", &
                          n_var=-1, type_of_var=real_t, repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="SPIN_CONFIGURATION", &
         description="For each singly occupied orbital, specify if this should be an alpha (=1) or a beta (=2) orbital. "// &
         "This keyword should be repeated, each repetition corresponding to an additional term.", &
         usage="SPIN_CONFIGURATION 1 2", &
         n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_low_spin_roks_section

! **************************************************************************************************
!> \brief Creates the section for applying an electrostatic external potential
!> \param section ...
!> \date 12.2009
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_ext_pot_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_pot_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_POTENTIAL", &
                          description="Section controlling the presence of an electrostatic "// &
                          "external potential dependent on the atomic positions (X,Y,Z). "// &
                          "As the external potential is currently applied via a grid, "// &
                          "it only works with DFT based methods (GPW/GAPW) that already use "// &
                          "a grid based approach to solve the Poisson equation.", &
                          n_keywords=7, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTION", &
                          description="Specifies the functional form in mathematical notation. Variables must be the atomic "// &
                          "coordinates (X,Y,Z) of the grid.", usage="FUNCTION  X^2+Y^2+Z^2+LOG(ABS(X+Y))", &
                          type_of_var=lchar_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAMETERS", &
                          description="Defines the parameters of the functional form", &
                          usage="PARAMETERS a b D", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VALUES", &
                          description="Defines the values of parameter of the functional form", &
                          usage="VALUES ", type_of_var=real_t, &
                          n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Optionally, allows to define valid CP2K unit strings for each parameter value. "// &
                          "It is assumed that the corresponding parameter value is specified in this unit.", &
                          usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATIC", &
                          description="Specifies the external potential as STATIC or time dependent. At the moment "// &
                          "only static potentials are implemented.", &
                          usage="STATIC T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DX", &
                          description="Parameter used for computing the derivative with the Ridders' method.", &
                          usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Checks that the error in computing the derivative is not larger than "// &
                          "the value set; in case error is larger a warning message is printed.", &
                          usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !keyword for reading the external potential from cube file
      CALL keyword_create(keyword, __LOCATION__, name="READ_FROM_CUBE", &
                          description="Switch for reading the external potential from file pot.cube. The values "// &
                          "of the potential must be on the grid points of the realspace grid.", &
                          usage="READ_FROM_CUBE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !keyword for scaling the external potential that is read from file by a constant factor
      CALL keyword_create(keyword, __LOCATION__, name="SCALING_FACTOR", &
                          description="A factor for scaling the the external potential that is read from file."// &
                          "The value of the potential at each grid point is multiplied by this factor.", &
                          usage="SCALING_FACTOR <REAL>", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_ext_pot_section

! **************************************************************************************************
!> \brief ZMP Creates the section for reading user supplied external density
!> \param section ...
!> \date 03.2011
!> \author D. Varsano [daniele.varsano@nano.cnr.it]
! **************************************************************************************************
   SUBROUTINE create_ext_den_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_den_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_DENSITY", &
                          description="Section for the use of the ZMP technique on external densities.", &
                          n_keywords=4, n_subsections=0, repeats=.FALSE., &
                          citations=(/Zhao1994, Tozer1996/))
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FILE_DENSITY", &
                          description="Specifies the filename containing the target density in *.cube format."// &
                          "In the MGRID section it must be imposed NGRID 1, as it works with only"// &
                          "one grid. The number of points in each direction, and the spacing must"// &
                          "be previously defined choosing the plane waves cut-off in section MGRID"// &
                          "keyword CUTOFF, and the cube dimension in section SUBSYS / CELL / keyword ABC", &
                          usage="DENSITY_FILE_NAME <FILENAME>", &
                          type_of_var=char_t, default_c_val="RHO_O.dat", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LAMBDA", &
                          description="Lagrange multiplier defined in the constraint ZMP method. When starting, use"// &
                          "small values when starting from scratch (around 5,10). Then gradually increase"// &
                          "the values depending, restarting from the previous calculation with the smaller"// &
                          "value. To choose the progressive values of LAMBDA look at the convergence of the"// &
                          " eigenvalues.", &
                          usage="DX <REAL>", default_r_val=10.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ZMP_CONSTRAINT", &
                          description="Specify which kind of constraint to solve the ZMP equation. The COULOMB default"// &
                          "option is more stable.", &
                          usage="ZMP_CONSTRAINT <CHAR>", &
                          enum_c_vals=s2a("COULOMB", "DIFF", "NONE"), &
                          enum_i_vals=(/use_coulomb, use_diff, use_no/), &
                          enum_desc=s2a("Coulomb constraint, integral of [rho_0(r)-rho(r)]/|r-r'|", &
                                        "Simple constraint, [rho_0(r)-rho(r)]", &
                                        "No constrain imposed"), &
                          default_i_val=use_coulomb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FERMI_AMALDI", &
                          description="Add the Fermi-Amaldi contribution to the Hartree potential."// &
                          "It leads to a more stable convergence.", &
                          usage="FERMI_AMALDI <LOGICAL>", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ext_den_section

! **************************************************************************************************
!> \brief ZMP Creates the section for creating the external v_xc
!> \param section ...
!> \date 03.2011
!> \author D. Varsano [daniele.varsano@nano.cnr.it]
! **************************************************************************************************
   SUBROUTINE create_ext_vxc_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_vxc_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_VXC", &
                          description="SCF convergence with external v_xc calculated through previous ZMP"// &
                          "calculation", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FILE_VXC", &
                          description="The *.cube filename containing the v_xc potential. This works only"// &
                          "with NGRID 1 imposed in the MGRID section. The number of points in each"// &
                          "direction, and the spacing must equal to those previously used in the ZMP"// &
                          "calculation and defined through the plane wave cut-off and the cube dimension"// &
                          "respectively set in section MGRID / keyword CUTOFF, and in section SUBSYS /"// &
                          "CELL / keyword ABC", &
                          usage="FILE_VXC <FILENAME>", &
                          type_of_var=char_t, default_c_val="VXC_O.dat", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_ext_vxc_section

! **************************************************************************************************
!> \brief makes the orbital transformation section
!> \param section ...
!> \par History
!>      11.2004 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE create_ot_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ot_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OT", &
                          description="Sets the various options for the orbital transformation (OT) method. "// &
                          "Default settings already provide an efficient, yet robust method. "// &
                          "Most systems benefit from using the FULL_ALL preconditioner "// &
                          "combined with a small value (0.001) of ENERGY_GAP. "// &
                          "Well-behaved systems might benefit from using a DIIS minimizer. "// &
                          "<u>Advantages:</u> "// &
                          "It's fast, because no expensive diagonalisation is performed. "// &
                          "If preconditioned correctly, method guaranteed to find minimum. "// &
                          "<u>Disadvantages:</u> "// &
                          "Sensitive to preconditioning. A good preconditioner can be expensive. "// &
                          "No smearing, or advanced SCF mixing possible: POOR convergence for metallic systems.", &
                          n_keywords=27, n_subsections=0, repeats=.FALSE., &
                          citations=(/VandeVondele2003, Weber2008/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the ot method", &
                          usage="&OT T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                          description="Algorithm to be used for OT", &
                          usage="ALGORITHM STRICT", &
                          default_i_val=ot_algo_taylor_or_diag, &
                          enum_c_vals=s2a("STRICT", "IRAC"), &
                          enum_desc=s2a("Strict orthogonality: Taylor or diagonalization based algorithm.", &
                                        "Orbital Transformation based Iterative Refinement "// &
                                        "of the Approximative Congruence transformation (OT/IR)."), &
                          enum_i_vals=(/ot_algo_taylor_or_diag, ot_algo_irac/), &
                          citations=(/VandeVondele2003, VandeVondele2005a, Weber2008/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IRAC_DEGREE", &
                          description="The refinement polynomial degree (2, 3 or 4).", &
                          usage="IRAC_DEGREE 4", &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_IRAC", &
                          description="Maximum allowed refinement iteration.", &
                          usage="MAX_IRAC 5", &
                          default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORTHO_IRAC", &
                          description="The orthogonality method.", &
                          usage="ORTHO_IRAC POLY", &
                          default_i_val=ot_chol_irac, &
                          enum_c_vals=s2a("CHOL", "POLY", "LWDN"), &
                          enum_desc=s2a("Cholesky.", "Polynomial.", "Loewdin."), &
                          enum_i_vals=(/ot_chol_irac, ot_poly_irac, ot_lwdn_irac/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC_FILTER_MATRIX", &
                          description="Sets the threshold for filtering the matrices.", &
                          usage="EPS_IRAC_FILTER_MATRIX 1.0E-5", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC", &
                          description="Targeted accuracy during the refinement iteration.", &
                          usage="EPS_IRAC 1.0E-5", &
                          default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC_QUICK_EXIT", &
                          description="Only one extra refinement iteration is "// &
                          "done when the norm is below this value.", &
                          usage="EPS_IRAC_QUICK_EXIT 1.0E-2", &
                          default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC_SWITCH", &
                          description="The algorithm switches to the polynomial "// &
                          "refinement when the norm is below this value.", &
                          usage="EPS_IRAC_SWITCH 1.0E-3", &
                          default_r_val=1.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ON_THE_FLY_LOC", &
                          description="On the fly localization of the molecular orbitals. "// &
                          "Can only be used with OT/IRAC.", &
                          usage="ON_THE_FLY_LOC T", &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MINIMIZER", &
         description="Minimizer to be used with the OT method", &
         usage="MINIMIZER DIIS", &
         default_i_val=ot_mini_cg, &
         enum_c_vals=s2a("SD", "CG", "DIIS", "BROYDEN"), &
         enum_desc=s2a("Steepest descent: not recommended", "Conjugate Gradients: most reliable, use for difficult systems."// &
                       " The total energy should decrease at every OT CG step if the line search is appropriate.", &
                       "Direct inversion in the iterative subspace: less reliable than CG, but sometimes about 50% faster", &
                       "Broyden mixing approximating the inverse Hessian"), &
         enum_i_vals=(/ot_mini_sd, ot_mini_cg, ot_mini_diis, ot_mini_broyden/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SAFE_DIIS", &
                          variants=(/"SAFER_DIIS"/), &
                          description="Reject DIIS steps if they point away from the"// &
                          " minimum, do SD in that case.", &
                          usage="SAFE_DIIS ON", default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_HISTORY_VEC", &
                          variants=s2a("NDIIS", "N_DIIS", "N_BROYDEN"), &
                          description="Number of history vectors to be used with DIIS or BROYDEN", &
                          usage="N_DIIS 4", &
                          default_i_val=7)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_BETA", &
                          description="Underrelaxation for the broyden mixer", &
                          usage="BROYDEN_BETA 0.9", &
                          default_r_val=0.9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_GAMMA", &
                          description="Backtracking parameter", &
                          usage="BROYDEN_GAMMA 0.5", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_SIGMA", &
                          description="Curvature of energy functional.", &
                          usage="BROYDEN_SIGMA 0.25", &
                          default_r_val=0.25_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_ETA", &
                          description="Dampening of estimated energy curvature.", &
                          usage="BROYDEN_ETA 0.7", &
                          default_r_val=0.7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_OMEGA", &
                          description="Growth limit of curvature.", &
                          usage="BROYDEN_OMEGA 1.1", &
                          default_r_val=1.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_SIGMA_DECREASE", &
                          description="Reduction of curvature on bad approximation.", &
                          usage="BROYDEN_SIGMA_DECREASE 0.7", &
                          default_r_val=0.7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_SIGMA_MIN", &
                          description="Minimum adaptive curvature.", &
                          usage="BROYDEN_SIGMA_MIN 0.05", &
                          default_r_val=0.05_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_FORGET_HISTORY", &
                          description="Forget history on bad approximation", &
                          usage="BROYDEN_FORGET_HISTORY OFF", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_ADAPTIVE_SIGMA", &
                          description="Enable adaptive curvature estimation", &
                          usage="BROYDEN_ADAPTIVE_SIGMA ON", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_ENABLE_FLIP", &
                          description="Ensure positive definite update", &
                          usage="BROYDEN_ENABLE_FLIP ON", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LINESEARCH", &
                          variants=(/"LINE_SEARCH"/), &
                          description="1D line search algorithm to be used with the OT minimizer,"// &
                          " in increasing order of robustness and cost. MINIMIZER CG combined with"// &
                          " LINESEARCH GOLD should always find an electronic minimum. "// &
                          " Whereas the 2PNT minimizer is almost always OK, 3PNT might be needed for systems"// &
                          " in which successive OT CG steps do not decrease the total energy.", &
                          usage="LINESEARCH GOLD", &
                          default_i_val=ls_2pnt, &
                          enum_c_vals=s2a("NONE", "2PNT", "3PNT", "GOLD"), &
                          enum_desc=s2a("take fixed length steps", "extrapolate based on 2 points", &
                                        "... or on 3 points", "perform 1D golden section search of the minimum (very expensive)"), &
                          enum_i_vals=(/ls_none, ls_2pnt, ls_3pnt, ls_gold/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="STEPSIZE", &
         description="Initial stepsize used for the line search, sometimes this parameter can be reduced to stabilize DIIS"// &
         " or to improve the CG behavior in the first few steps."// &
         " The optimal value depends on the quality of the preconditioner."// &
         " A negative values leaves the choice to CP2K depending on the preconditioner.", &
         usage="STEPSIZE 0.4", &
         default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GOLD_TARGET", &
                          description="Target relative uncertainty in the location of the minimum for LINESEARCH GOLD", &
                          usage="GOLD_TARGET 0.1", &
                          default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with all minimization schemes. "// &
         "They differ in effectiveness, cost of construction, cost of application. "// &
         "Properly preconditioned minimization can be orders of magnitude faster than doing nothing.", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=ot_precond_full_kinetic, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", &
                         "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "// &
                       "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "// &
                       "This preconditioner is recommended for almost all systems, except very large systems where "// &
                       "make_preconditioner would dominate the total computational cost.", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, "// &
                       "might be somewhat less robust. Recommended for large systems.", &
                       "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ", &
                       "Cholesky inversion of S and T, fast construction, robust, and relatively good, "// &
                       "use for very large systems.", &
                       "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_full_single, &
                       ot_precond_full_kinetic, ot_precond_s_inverse, ot_precond_none/), &
         citations=(/VandeVondele2003, Weber2008, Schiffmann2015/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHOLESKY", &
                          description="If FULL_ALL the cholesky decomposition of the S matrix is used. "// &
                          "Options on the algorithm to be used.", &
                          usage="CHOLESKY REDUCE", default_i_val=cholesky_reduce, &
                          enum_c_vals=s2a("OFF", "REDUCE", "RESTORE", "INVERSE", "INVERSE_DBCSR"), &
                          enum_desc=s2a("The cholesky algorithm is not used", "Reduce is called", &
                                        "Reduce is replaced by two restore", &
                                        "Restore uses operator multiply by inverse of the triangular matrix", &
                                        "Like inverse, but matrix stored as dbcsr, sparce matrix algebra used when possible"), &
                          enum_i_vals=(/cholesky_off, cholesky_reduce, cholesky_restore, cholesky_inverse, cholesky_dbcsr/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECOND_SOLVER", &
         description="How the preconditioner is applied to the residual.", &
         usage="PRECOND_SOLVER DIRECT", &
         default_i_val=ot_precond_solver_default, &
         enum_c_vals=s2a("DEFAULT", "DIRECT", "INVERSE_CHOLESKY", "INVERSE_UPDATE"), &
         enum_desc=s2a("the default", "Cholesky decomposition followed by triangular solve "// &
                       "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)", &
                       "Cholesky decomposition followed by explicit inversion "// &
                       "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)", &
                       "Performs a Hotelling update of the inverse if a previous preconditioner is present. "// &
                       "Mainly useful for GPU accelerated systems (works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)"), &
         enum_i_vals=(/ot_precond_solver_default, &
                       ot_precond_solver_direct, &
                       ot_precond_solver_inv_chol, &
                       ot_precond_solver_update/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="ENERGY_GAP", &
         description="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, "// &
         "especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate "// &
         "of the gap (can be a small number, e.g. 0.002)."// &
         " FULL_SINGLE_INVERSE takes it as lower bound (values below 0.05 can cause stability issues)."// &
         " In general, higher values will tame the preconditioner in case of poor initial guesses."// &
         " A negative value will leave the choice to CP2K depending on type of preconditioner.", &
         usage="ENERGY_GAP 0.001", &
         default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_TAYLOR", &
         variants=(/"EPSTAYLOR"/), &
         description="Target accuracy of the taylor expansion for the matrix functions, should normally be kept as is.", &
         usage="EPS_TAYLOR 1.0E-15", &
         default_r_val=1.0E-16_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MAX_TAYLOR", &
         description="Maximum order of the Taylor expansion before diagonalisation is preferred, for large parallel runs"// &
         " a slightly higher order could sometimes result in a small speedup.", &
         usage="MAX_TAYLOR 5", &
         default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ROTATION", &
                          description="Introduce additional variables so that rotations of the occupied"// &
                          " subspace are allowed as well, only needed for cases where the energy is not invariant under "// &
                          " a rotation of the occupied subspace such as non-singlet restricted calculations "// &
                          " or fractional occupations.", &
                          usage="ROTATION", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGIES", &
                          description="Optimize orbital energies for use in Fermi-Dirac smearing "// &
                          "(requires ROTATION and FD smearing to be active).", &
                          usage="ENERGIES", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OCCUPATION_PRECONDITIONER", &
                          description="Preconditioner with the occupation numbers (FD smearing)", &
                          usage="OCCUPATION_PRECONDITIONER", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NONDIAG_ENERGY", &
                          description="Add a non-diagonal energy penalty (FD smearing)", &
                          usage="NONDIAG_ENERGY", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NONDIAG_ENERGY_STRENGTH", &
                          description="The prefactor for the non-diagonal energy penalty (FD smearing)", &
                          usage="NONDIAG_ENERGY_STRENGTH", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ot_section

! **************************************************************************************************
!> \brief creates the diagonalization section
!> \param section ...
!> \par History
!>      10.2008 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_diagonalization_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_diagonalization_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DIAGONALIZATION", &
                          description="Set up type and parameters for Kohn-Sham matrix diagonalization.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the diagonalization method", &
                          usage="&DIAGONALIZATION T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                          description="Algorithm to be used for diagonalization", &
                          usage="ALGORITHM STANDARD", &
                          default_i_val=diag_standard, &
                          enum_c_vals=s2a("STANDARD", "OT", "LANCZOS", "DAVIDSON", "FILTER_MATRIX"), &
                          enum_desc=s2a("Standard diagonalization: LAPACK methods or Jacobi.", &
                                        "Iterative diagonalization using OT method", &
                                        "Block Krylov-space approach to self-consistent diagonalisation", &
                                        "Preconditioned blocked Davidson", &
                                        "Filter matrix diagonalization"), &
                          enum_i_vals=(/diag_standard, diag_ot, diag_block_krylov, diag_block_davidson, &
                                        diag_filter_matrix/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBI_THRESHOLD", &
                          description="Controls the accuracy of the pseudo-diagonalization method using Jacobi rotations", &
                          usage="JACOBI_THRESHOLD 1.0E-6", &
                          default_r_val=1.0E-7_dp, &
                          citations=(/Stewart1982/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_JACOBI", &
                          description="Below this threshold value for the SCF convergence the pseudo-diagonalization "// &
                          "method using Jacobi rotations is activated. This method is much faster than a "// &
                          "real diagonalization and it is even speeding up while achieving full convergence."// &
                          "However, it needs a pre-converged wavefunction obtained by at least one real "// &
                          "diagonalization which is further optimized while keeping the original eigenvalue "// &
                          "spectrum. The MO eigenvalues are NOT updated. The method might be useful to speed "// &
                          "up calculations for large systems e.g. using a semi-empirical method.", &
                          usage="EPS_JACOBI 1.0E-5", &
                          default_r_val=0.0_dp, &
                          citations=(/Stewart1982/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ADAPT", &
                          description="Required accuracy in iterative diagonalization as compared to current SCF convergence", &
                          usage="EPS_ADAPT 0.01", &
                          default_r_val=0._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations in iterative diagonalization", &
                          usage="MAX_ITER 20", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Required accuracy in iterative diagonalization", &
                          usage="EPS_ITER 1.e-8", &
                          default_r_val=1.e-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_ot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_krylov_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_diag_subspace_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_davidson_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_filtermatrix_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_diagonalization_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_davidson_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_davidson_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DAVIDSON", &
                          description=" ", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with all minimization schemes. ", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=ot_precond_full_all, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization ", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, might be somewhat less robust. Recommended for large systems.", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_none/), &
         citations=(/VandeVondele2003/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRECOND_SOLVER", &
                          description="How the preconditioner is applied to the residual.", &
                          usage="PRECOND_SOLVER DIRECT", &
                          default_i_val=ot_precond_solver_default, &
                          enum_c_vals=s2a("DEFAULT", "DIRECT", "INVERSE_CHOLESKY"), &
                          enum_desc=s2a("the default", "Cholesky decomposition followed by triangular solve "// &
                                        "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)", &
                                        "Cholesky decomposition followed by explicit inversion "// &
                                        "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)"), &
                          enum_i_vals=(/ot_precond_solver_default, &
                                        ot_precond_solver_direct, &
                                        ot_precond_solver_inv_chol/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="ENERGY_GAP", &
         description="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, "// &
         "especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate "// &
         "of the gap (0.001 doing normally fine). For the other preconditioners, making this value larger (0.2)"// &
         " will tame the preconditioner in case of poor initial guesses.", &
         usage="ENERGY_GAP 0.001", &
         default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NEW_PREC_EACH", &
                          description="Number of SCF iterations after which a new Preconditioner is computed", &
                          usage="NEW_PREC_EACH 10", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIRST_PREC", &
                          description="First SCF iteration at which a Preconditioner is employed", &
                          usage="FIRST_PREC 1", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONV_MOS_PERCENT", &
                          description="Minimal percent of MOS that have to converge within the Davidson loop"// &
                          " before the SCF iteration is completed and a new Hamiltonian is computed", &
                          usage="CONV_MOS_PERCENT 0.8", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPARSE_MOS", &
                          description="Use MOS as sparse matrix and avoid as much as possible multiplications with full matrices", &
                          usage="SPARSE_MOS", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_davidson_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_krylov_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_krylov_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="KRYLOV", &
                          description=" ", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NKRYLOV", &
                          description="Dimension of the Krylov space used for the Lanczos refinement", &
                          usage="NKRYLOV 20", &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NBLOCK", &
                          description="Size of the block of vectors refined simultaneously by the Lanczos procedure", &
                          usage="NBLOCK 1", &
                          default_i_val=32)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_KRYLOV", &
                          description="Convergence criterion for the MOs", &
                          usage="EPS_KRYLOV 0.00001", &
                          default_r_val=0.0000001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_STD_DIAG", &
                          description="Level of convergence to be reached before starting the Lanczos procedure."// &
                          " Above this threshold a standard diagonalization method is used. "// &
                          " If negative Lanczos is started at the first iteration", &
                          usage="EPS_STD_DIAG 0.001", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_MOS_CONV", &
                          description="This requires to check the convergence of MOS also when standard "// &
                          "diagonalization steps are performed, if the block krylov approach is active.", &
                          usage="CHECK_MOS_CONV T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_krylov_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_diag_subspace_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_diag_subspace_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DIAG_SUB_SCF", &
                          description="Activation of self-consistenf subspace refinement by diagonalization "// &
                          "of H by adjusting the occupation but keeping the MOS unchanged.", &
                          n_keywords=2, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of inner SCF loop to refine occupations in MOS subspace", &
                          usage="&DIAG_SUB_SCF T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations for the SCF inner loop", &
                          usage="MAX_ITER 20", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ENE", &
                          description="Required energy accuracy for convergence of subspace diagonalization", &
                          usage="EPS_ENE 1.e-8", &
                          default_r_val=1.e-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ADAPT_SCF", &
                          description="Required density matrix accuracy as compared to current SCF convergence", &
                          usage="EPS_ADAPT_SCF 1.e-1", &
                          default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_SKIP_SUB_DIAG", &
         description="Level of convergence to be reached before starting the internal loop of subspace rotations."// &
         " Above this threshold only the outer diagonalization method is used. "// &
         " If negative the subspace rotation is started at the first iteration", &
         usage="EPS_SKIP_SUB_DIAG 0.001", &
         default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_mixing_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
   END SUBROUTINE create_diag_subspace_section

! **************************************************************************************************
!> \brief makes the input section for core-level spectroscopy simulations
!> \param section ...
!> \par History
!>      03.2005 created [MI]
! **************************************************************************************************
   SUBROUTINE create_xas_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_xas_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xas", &
                          description="Sets the method of choice to calculate core-level excitation spectra. "// &
                          "The occupied states from  which we calculate the "// &
                          "excitation should be specified. "// &
                          "Localization of the orbitals may be useful.", &
                          n_keywords=10, n_subsections=1, repeats=.FALSE., &
                          citations=(/Iannuzzi2007/))

      NULLIFY (keyword, subsection, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of core-level spectroscopy simulations", &
                          usage="&XAS T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          variants=(/"XAS_METHOD"/), &
                          description="Method to be used to calculate core-level excitation spectra", &
                          usage="METHOD TP_HH", &
                          default_i_val=xas_none, &
                          enum_c_vals=s2a("NONE", "TP_HH", "TP_FH", "TP_VAL", "TP_XHH", "TP_XFH", "DSCF", "TP_FLEX"), &
                          enum_desc=s2a( &
                          "No core electron spectroscopy", "Transition potential half-hole", &
                          "Transition potential full-hole", "Hole in homo for X-ray emission only ", &
                          "Transition potential excited half-hole", &
                          "Transition potential excited full-hole ", &
                          "DSCF calculations to compute the first (core)excited state", &
                          "Transition potential with generalized core occupation and total number of electrons"), &
                          enum_i_vals=(/xas_none, xas_tp_hh, xas_tp_fh, xes_tp_val, xas_tp_xhh, &
                                        xas_tp_xfh, xas_dscf, xas_tp_flex/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XAS_CORE", &
                          description="Occupation of the core state in XAS calculation by TP_FLEX.", &
                          usage="XAS_CORE 0.5", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XAS_TOT_EL", &
                          description="Total number of electrons for spin channel alpha, in XAS calculation by TP_FLEX."// &
                          "If it is a negative value, the number of electrons is set to GS number of electrons "// &
                          "minus the amount subtracted from the core state", &
                          usage="XAS_TOT_EL 10", &
                          default_r_val=-1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XES_CORE", &
                          description="Occupation of the core state in XES calculation by TP_VAL. "// &
                          "The HOMO is emptied by the same amount.", &
                          usage="XES_CORE 0.5", &
                          default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XES_EMPTY_HOMO", &
                          description="Set the occupation of the HOMO in XES calculation by TP_VAL. "// &
                          "The HOMO can be emptied or not, if the core is still full.", &
                          usage="XES_EMPTY_HOMO", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_FORM", &
                          variants=(/"DIP_FORM"/), &
                          description="Type of integral to get the oscillator strengths "// &
                          "in the diipole approximation", &
                          usage="DIPOLE_FORM string", &
                          default_i_val=xas_dip_vel, &
                          enum_c_vals=s2a("LENGTH", "VELOCITY"), &
                          enum_desc=s2a("Length form &lang; i | e r | j &rang;", &
                                        "Velocity form &lang; i | d/dr | j &rang;"), &
                          enum_i_vals=(/xas_dip_len, xas_dip_vel/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!   replace the specialized keyword with standard scf section
!    scf_env is added to xas_env

      NULLIFY (subsection)
      CALL create_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_TYPE", &
                          variants=(/"TYPE"/), &
                          description="Type of the orbitals that are excited for the xas spectra calculation", &
                          usage="STATE_TYPE 1S", &
                          default_i_val=xas_1s_type, &
                          enum_c_vals=s2a("1S", "2S", "2P", "3S", "3P", "3D", "4S", "4P", "4D", "4F"), &
                          enum_desc=s2a("1s orbitals", "2s orbitals", "2p orbitals", "3s orbitals", "3p orbitals", &
                                        "3d orbitals", "4s orbitals", "4p orbitals", "4d orbitals", "4f orbitals"), &
                          enum_i_vals=(/xas_1s_type, xas_2s_type, xas_2p_type, xas_3s_type, xas_3p_type, xas_3d_type, &
                                        xas_4s_type, xas_4p_type, xas_4d_type, xas_4f_type/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_SEARCH", &
                          description="# of states where to look for the one to be excited", &
                          usage="STATE_SEARCH 1", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_CHANNEL", &
                          description="# Spin channel of the excited orbital", &
                          usage="SPIN_CHANNEL 1", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_LIST", &
                          variants=(/"AT_LIST"/), &
                          description="Indexes of the atoms to be excited. "// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="ATOMS_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OVERLAP_THRESHOLD", &
                          description="Threshold for including more than one initial core excited state "// &
                          "per atom. The threshold is taken relative to the maximum overlap.", &
                          usage="OVERLAP_THRESHOLD 8.e-1", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_LIST", &
                          variants=(/"ORBITAL_LIST"/), &
                          description="Indices of the localized orbitals to be excited. "// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="ORBITAL_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ADDED_MOS", &
                          description="Number of additional MOS added spin up only", &
                          usage="ADDED_MOS {integer}", default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_ADDED", &
                          description="maximum number of iteration in calculation of added orbitals", &
                          usage="MAX_ITER_ADDED 100", default_i_val=2999)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ADDED", &
                          description="target accuracy incalculation of the added orbitals", &
                          usage="EPS_ADDED 1.e-6", default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NGAUSS", &
                          description="Number of gto's for the expansion of the STO "// &
                          "of the type given by STATE_TYPE", &
                          usage="NGAUSS {integer}", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the excited state if the restart file exists", &
                          usage="RESTART", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          variants=(/"RESTART_FILE_NAME"/), &
                          description="Root of the file names where to read the MOS from "// &
                          "which to restart the calculation of the core level excited states", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="printing of information during the core-level spectroscopy simulation", &
                          repeats=.FALSE.)

      ! Add printing of wannier infos
      CALL print_wanniers(subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "iteration_info", &
                                       description="Controls the printing of basic iteration information during the XAS SCF.", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="time_cumul", &
                          description="If the printkey is activated switches the printing of timings"// &
                          " to cumulative (over the SCF).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing of basic iteration information in CLS", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "XES_SPECTRUM", &
                                       description="Controls the dumping of the CLS output files containing the emission spectra", &
                                       print_level=low_print_level, common_iter_levels=3, filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "XAS_SPECTRUM", &
         description="Controls the dumping of the CLS output files containing the absorption spectra", &
         print_level=low_print_level, common_iter_levels=3, filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL create_pdos_section(print_key)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of MO restart file during the SCF."// &
                                       "of a Core-Level-Spectroscopy calculation. For each new excited atom,"// &
                                       "one different restart file is dumped. These restart files should be"// &
                                       "employed only to restart the same type of CLS calculation, "// &
                                       "i.e. with the same core potential.", &
                                       print_level=low_print_level, common_iter_levels=3, each_iter_names=s2a("XAS_SCF"), &
                                       add_last=add_last_numeric, each_iter_values=(/3/), filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CLS_FUNCTION_CUBES", &
                                       description="Controls the printing of the relaxed orbitals ", &
                                       print_level=high_print_level, common_iter_levels=3, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LU_BOUNDS", &
                          variants=(/"CUBES_LU"/), &
                          description="The lower and upper index of the states to be printed as cube", &
                          usage="CUBES_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/0, -2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LIST", &
                          description="Indexes of the states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="CUBES_LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xas_section

! **************************************************************************************************
!> \brief makes the input section for core-level spectroscopy simulations using
!>        linear response TDDFT
!> \param section ...
!> \par History
!>      11.2017 created [AB]
! **************************************************************************************************
   SUBROUTINE create_xas_tdp_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_xas_tdp_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection, subsubsection, &
                                                            subsubsubsection

      NULLIFY (keyword, print_key, subsection, subsubsection, subsubsubsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="XAS_TDP", &
                          description="XAS simulations using linear-response TDDFT. Excitation from "// &
                          "specified core orbitals are considered one at a time. In case of high "// &
                          "symmetry structures, donor core orbitals should be localized.", &
                          n_keywords=15, n_subsections=3, repeats=.FALSE.)

      NULLIFY (keyword, subsection, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of XAS simulations with linear "// &
                          "response TDDFT.", &
                          usage="&TDP_XAS {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_ONLY", &
                          description="This keyword defines whether the full calculation should "// &
                          "be done or not. If set to .TRUE., only the determination "// &
                          "of donor MOs is conducted. This run option allows for "// &
                          "cheap verification of the input parameters", &
                          usage="CHECK_ONLY {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_FROM_FILE", &
                          variants=s2a("RESTART_FILENAME", "RST_FILENAME", "RESTART_FILE", "RST_FILE"), &
                          description="By providing a RESTART file containing the linear-response "// &
                          "orbitals and excitations energies from a previous calculation, "// &
                          "all computations are skipped except for the corresponding "// &
                          "PDOS and/or CUBE file printing as defined in the PRINT "// &
                          "subsection. Basis sets and geometry need to be consistent.", &
                          usage="RESTART_FROM_FILE <FILENAME>", &
                          type_of_var=char_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXCITATIONS", &
                          variants=(/"EXCITATION"/), &
                          description="Specify the type of excitation to consider. In case of a "// &
                          "resctricted closed-shell ground state calculation, "// &
                          "RCS_SINGLET or/and RCS_TRIPLET can be chosen. In case of a "// &
                          "open-shell ground state calculation (either UKS or ROKS), "// &
                          "standard spin conserving excitation (OS_SPIN_CONS) or/and "// &
                          "spin-flip excitation (OS_SPIN_FLIP) can be chosen.", &
                          usage="EXCITATIONS {string}", &
                          repeats=.TRUE., &
                          default_i_val=tddfpt_singlet, &
                          enum_c_vals=s2a("RCS_SINGLET", "RCS_TRIPLET", "OS_SPIN_CONS", "OS_SPIN_FLIP"), &
                          enum_desc=s2a("Singlet excitation on top of restricted closed-shell ground state", &
                                        "Triplet excitation on top of restricted closed-shell ground state", &
                                        "Spin-conserving excitations on top of open-shell ground state", &
                                        "Spin-flip excitation on top of open-shell ground state"), &
                          enum_i_vals=(/tddfpt_singlet, tddfpt_triplet, tddfpt_spin_cons, tddfpt_spin_flip/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_PGF_XAS", &
                          variants=s2a("EPS_PGF", "EPS_PGF_XAS_TDP"), &
                          description="The threshold used to determine the spacial extent of all "// &
                          "primitive Gaussian functions used for the construction "// &
                          "of neighbor lists in the XAS_TDP method. "// &
                          "By default, takes the value of QS%EPS_PGF_ORB. Useful if "// &
                          "the former value is tiny due to possible ground state HFX "// &
                          "contributions.", &
                          usage="EPS_PGS_XAS {real}", &
                          type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          variants=s2a("EPS_FILTER_MATRIX"), &
                          description="The threshold used for sparse matrix operations", &
                          usage="EPS_FILTER {real}", &
                          type_of_var=real_t, &
                          default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_FORM", &
                          variants=(/"DIP_FORM"/), &
                          description="Type of integral to get the oscillator strengths "// &
                          "in the dipole approximation", &
                          usage="DIPOLE_FORM {string}", &
                          default_i_val=xas_dip_vel, &
                          enum_c_vals=s2a("LENGTH", "VELOCITY"), &
                          enum_desc=s2a("Length form &lang; 0 | e r | n &rang;", &
                                        "Velocity form &lang; 0 | d/dr | n &rang;"), &
                          enum_i_vals=(/xas_dip_len, xas_dip_vel/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="QUADRUPOLE", &
                          variants=s2a("DO_QUADRUPOLE", "DO_QUAD", "QUAD"), &
                          description="Compute the electric quadrupole contribution to the "// &
                          "oscillator strenghts (in the length representation with "// &
                          "the origin set on the relevant excited atom)", &
                          usage="QUADRUPOLE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!  The donor state subsection

      CALL section_create(subsection, __LOCATION__, name="DONOR_STATES", &
                          description="Specifications for the donor states from which core "// &
                          "electrons are excited", &
                          n_keywords=6, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="DEFINE_EXCITED", &
                          description="Whether the atoms to be excited should be defined by "// &
                          "a list of atom indices or by a list of atom kinds.", &
                          usage="DEFINE_EXCITED {string}", &
                          default_i_val=xas_tdp_by_index, &
                          enum_c_vals=s2a("BY_INDEX", "BY_KIND"), &
                          enum_i_vals=(/xas_tdp_by_index, xas_tdp_by_kind/), &
                          enum_desc=s2a("Excited atoms are defined by a list of indices", &
                                        "Excited atoms are defined by a list of atomic kinds"))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          variants=(/"AT_LIST"/), &
                          description="Indices of the atoms to be excited. "// &
                          "Keyword only taken into account if DEFINE_EXCITED = BY_INDEX", &
                          usage="ATOM_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KIND_LIST", &
                          description="Kind of atoms to be excited. "// &
                          "All atoms of the specified kinds are considered. "// &
                          "Keyword only taken into account if DEFINE_EXCITED = BY_KIND", &
                          usage="KIND_LIST {string}  {string} ..  {string} ", &
                          n_var=-1, type_of_var=char_t, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_TYPES", &
                          variants=(/"TYPES"/), &
                          description="Types of orbitals that are excited, for each atom/kind, "// &
                          "in order to do LR-TDDFT driven xas spectra calculation. "// &
                          "This keyword MUST have the same number of entries as the relevant "// &
                          "KIND_LIST or ATOM_LIST. The order of the specified state types must "// &
                          "correspond to the order of the relevant kinds/indices. "// &
                          "This keyword can be repeated, useful when multiple orbital types "// &
                          "should be excited for specific kinds/atoms.", &
                          n_var=-1, default_i_val=xas_not_excited, repeats=.TRUE., &
                          usage="STATE_TYPES {string}  {string} .. {string}", &
                          enum_c_vals=s2a("1S", "2S", "2P", "NE"), &
                          enum_desc=s2a("1s orbital", "2s orbital", "2p orbitals", "not excited"), &
                          enum_i_vals=(/xas_1s_type, xas_2s_type, xas_2p_type, xas_not_excited/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_SEARCH", &
                          description="Number of MOs (per spin) to search to find specified donor core "// &
                          "orbitals, starting from the lowest in energy and upward. By default, "// &
                          "all HOMOs are searched. If the LOCALIZE keyword is used, "// &
                          "then all searched states are first localized.", &
                          usage="N_SEARCH {integer}", &
                          default_i_val=-1)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOCALIZE", &
                          variants=s2a("LOC", "DO_LOC"), &
                          description="Whether the N_SEARCH potential donor states are should be "// &
                          "actively localized. Necessary in case of excited atoms "// &
                          "equivalent under symmetry", &
                          usage="LOCALIZE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of the donor states subsection

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_ORBIT_COUPLING", &
                          variants=(/"SOC"/), &
                          description="Whether spin-orbit coupling should be added. "// &
                          "Note: only applies for spin-restricted calculations with "// &
                          "singlet and triplet excitations OR spin-unrestricted "// &
                          "calculations with both spin-conserving and spin-flip."// &
                          "Only tested with CIS", &
                          usage="SOC {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TAMM_DANCOFF", &
                          variants=(/"TDA"/), &
                          description="Whether the Tamm-Dancoff approximation should be used.", &
                          usage="TAMM_DANCOFF {logical}", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID", &
                          variants=(/"ATOMIC_GRID"/), &
                          description="Specification of the atomic angular and radial grids for "// &
                          "a given atomic kind. This keyword can/should be repeated "// &
                          "for each excited kind. The default grid dimensions are "// &
                          "those set for the GAPW ground state calculation. These "// &
                          "grids are used for the xc-kernel calculations and/or for "// &
                          "spin-orbit coupling. "// &
                          "Usage: GRID < KIND > < LEBEDEV_GRID > < RADIAL_GRID >", &
                          usage="GRID {string} {integer} {integer}", &
                          n_var=3, type_of_var=char_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_EXCITED", &
                          variants=(/"N_ROOTS"/), &
                          description="The number of excited states to compute per donor "// &
                          "molecular orbital. (e.g. if 2p excitations, "// &
                          "3*N_EXCITED excited states are considered). "// &
                          "If N_EXCITED is set to -1, all excitations are considered", &
                          usage="N_EXCITED {integer}", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_RANGE", &
                          variants=s2a("E_RANGE"), &
                          description="The energy range in eV for which excitations are considered. "// &
                          "Only excitated states within the range of: first excitation "// &
                          "energy + ENERGY_RANGE are kept. If ENERGY_RANGE "// &
                          "and N_EXCITED are specified, the former has priority. "// &
                          "Negative values are ignored and N_EXCITED takes over.", &
                          usage="ENERGY_RANGE {real}", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!  The KERNEL subsection
      CALL section_create(subsection, __LOCATION__, name="KERNEL", &
                          description="Defines how the kernel is built in terms of functionals.", &
                          n_keywords=1, &
                          n_subsections=1, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="RI_REGION", &
                          variants=(/"RI_RADIUS"/), &
                          description="The region defined by a sphere of the given radius around "// &
                          "the excited atom defining which RI_XAS basis elements are "// &
                          "considered for the RI projection of the density. Each basis "// &
                          "element which center is in this region is taken. The density "// &
                          "for a given excited atom is expressed as : "// &
                          "sum_ijkl P_ij (ijk) S_kl^-1 xi_l, where P_ij is the density "// &
                          "matrix, i,j span the orbital basis and k,l the RI_XAS basis "// &
                          "in the region. The larger the radius, the more basis "// &
                          "functions to expand the density. However, it is assumed "// &
                          "that it is a small number and the code does not scale well "// &
                          "as the number of RI basis elements gets too large. "// &
                          "Expressed in Angstrom. If the radius is set to 0.0, only "// &
                          "the RI basis elements centered on the excited atom are used.", &
                          usage="RI_REGION {real}", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

!  The XC_FUNCTIONAL subsubsection
      CALL section_create(subsubsection, __LOCATION__, name="XC_FUNCTIONAL", &
                          description="The exchange correlation functional(s) to use for the xc "// &
                          "part of the kernel. Built from any combination of LDA, GGA "// &
                          "and HYB functional from LIBXC.", &
                          n_keywords=0, &
                          n_subsections=1, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the addition of xc functionals to the kernel.", &
                          usage="&XC_FUNCTIONAL {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      ! LIBXC functional
      CALL create_libxc_section(subsubsubsection, "LIBXC", "Uses functionals from LIBXC", &
                                "FUNCTIONAL {name}")
      CALL section_add_subsection(subsubsection, subsubsubsection)
      CALL section_release(subsubsubsection)

      CALL section_add_subsection(subsection, subsubsection)
      CALL section_release(subsubsection)

! The EXACT_EXCHANGE subsubsection
      CALL section_create(subsubsection, __LOCATION__, name="EXACT_EXCHANGE", &
                          description="Whether exact-exchange should be added to the kernel and "// &
                          "if so, with which fraction and operator.", &
                          n_keywords=7, &
                          n_subsections=1, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the addition of exact exchange to the kernel.", &
                          usage="&EXACT_EXCHANGE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_TYPE", &
                          variants=s2a("OP", "OPERATOR", "POTENTIAL"), &
                          description="The type of operator used for exact exchange. The standard "// &
                          "Coulomb operator cannot be used in periodic systems.", &
                          usage="OPERATOR {string}", &
                          repeats=.FALSE., &
                          default_i_val=do_potential_coulomb, &
                          enum_c_vals=s2a("COULOMB", "TRUNCATED", "SHORTRANGE"), &
                          enum_desc=s2a("Standard Coulomb operator: 1/r", &
                                        "Truncated Coulomb operator: 1/r if r < R_c, 0 otherwise ", &
                                        "Short range: erfc(omega*r)/r"), &
                          enum_i_vals=(/do_potential_coulomb, do_potential_truncated, &
                                        do_potential_short/))
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          variants=s2a("R_C", "RC", "RANGE"), &
                          description="The cutoff radius (in Angstrom) for the truncated Coulomb operator.", &
                          usage="CUTOFF_RADIUS {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="T_C_G_DATA", &
                          description="Location of the file t_c_g.dat that contains the data for the "// &
                          "evaluation of the truncated gamma function ", &
                          usage="T_C_G_DATA {string}", &
                          default_c_val="t_c_g.dat")
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OMEGA", &
                          description="The range parameter for the short range operator (in 1/a0).", &
                          usage="OMEGA {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RANGE", &
                          description="The threshold to determine the effective range of the short range "// &
                          "operator: erfc(omega*eff_range)/eff_range = EPS_RANGE", &
                          usage="EPS_RANGE = {double}", &
                          default_r_val=1.0E-6_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCREENING", &
                          variants=s2a("EPS_SCREEN"), &
                          description="A threshold to determine which primitive 3-center integrals "// &
                          "are kept for contraction, as the latter operation can be "// &
                          "expensive (especially for large basis sets )."// &
                          "If |(ab|c)| < EPS_SCREENNING, it is discarded.", &
                          default_r_val=1.0E-8_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          variants=s2a("FRACTION"), &
                          description="Scaling of the exact exchange contribution.", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      !The RI metric subsection
      CALL section_create(subsubsubsection, __LOCATION__, name="RI_METRIC", &
                          description="This subsection allows for the definition of an exchange "// &
                          "RI metric that is different from the main exchange potential. "// &
                          "By default (i.e. if this subsection is ignored), the "// &
                          "exchange kernel is computed in the V approximation: "// &
                          "(ab|ij) = (ab|P) V^-1 (Q|ij), where V = (P|Q). With a RI "// &
                          "metric, we have a 2 step RI involving the metric potential "// &
                          "for the 3-center integrals: "// &
                          "(ab|ij) = (ab!P) (P!Q)^-1 (Q|R) (R!S)^-1 (S!ij), where | "// &
                          "stands for the exchange potential and ! for the metric "// &
                          "potential. This allows for drastic screening of the "// &
                          "3-center integrals by selecting shorter range metric.", &
                          n_keywords=5, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the use of a RI metric.", &
                          usage="&EXACT_EXCHANGE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_TYPE", &
                          variants=s2a("OP", "OPERATOR", "POTENTIAL"), &
                          description="The type of operator used for the metric.", &
                          usage="OPERATOR {string}", &
                          repeats=.FALSE., &
                          default_i_val=do_potential_id, &
                          enum_c_vals=s2a("OVERLAP", "TRUNCATED", "SHORTRANGE"), &
                          enum_desc=s2a("Overlap operator (=identity)", &
                                        "Truncated Coulomb operator: 1/r if r < R_c, 0 otherwise ", &
                                        "Short range: erfc(omega*r)/r"), &
                          enum_i_vals=(/do_potential_id, do_potential_truncated, &
                                        do_potential_short/))
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          variants=s2a("R_C", "RC", "RANGE"), &
                          description="The cutoff radius (in Angstrom) for the truncated Coulomb operator.", &
                          usage="CUTOFF_RADIUS {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="T_C_G_DATA", &
                          description="Location of the file t_c_g.dat that contains the data for the "// &
                          "evaluation of the truncated gamma function ", &
                          usage="T_C_G_DATA {string}", &
                          default_c_val="t_c_g.dat")
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OMEGA", &
                          description="The range parameter for the short range operator (in 1/a0).", &
                          usage="OMEGA {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsubsection, subsubsubsection)
      CALL section_release(subsubsubsection)

      CALL section_add_subsection(subsection, subsubsection)
      CALL section_release(subsubsection)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of Kernel subsection

      CALL section_create(subsection, __LOCATION__, "PRINT", "Controls the printing of information during "// &
                          "XAS TDP calculations", repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="SPECTRUM", &
                                       description="Controles the dumping of the XAS TDP spectrum in output files", &
                                       print_level=low_print_level, filename="", common_iter_levels=3)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL create_pdos_section(print_key)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CUBES", &
                                       description="Controls the printing of the linear-response orbitals "// &
                                       "as *.cube files.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LU_BOUNDS", &
                          variants=(/"CUBES_LU"/), &
                          description="The lower and upper index of the excited states to be printed as cube", &
                          usage="CUBES_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/1, 0/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LIST", &
                          description="Indexes of the excited states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="CUBES_LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of LR-orbitals coefficients "// &
                                       "and corresponding excitation energies such that "// &
                                       "the program can be restarted for PDOS or CUBE "// &
                                       "printing without the heavy computing.", &
                                       print_level=high_print_level, filename="", common_iter_levels=3)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xas_tdp_section

! **************************************************************************************************
!> \brief      Create CP2K input section for the smearing of occupation numbers
!> \param section ...
!> \date       27.08.2008
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE create_smear_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_smear_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, &
                          name="SMEAR", &
                          description="Define the smearing of the MO occupation numbers", &
                          n_keywords=6, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of smearing", &
                          usage="&SMEAR ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="METHOD", &
                          description="Smearing method to be applied", &
                          usage="METHOD Fermi_Dirac", &
                          default_i_val=smear_energy_window, &
                          enum_c_vals=s2a("FERMI_DIRAC", "ENERGY_WINDOW", "LIST"), &
                          enum_i_vals=(/smear_fermi_dirac, smear_energy_window, smear_list/), &
                          enum_desc=s2a("Fermi-Dirac distribution defined by the keyword ELECTRONIC_TEMPERATURE", &
                                        "Energy window defined by the keyword WINDOW_SIZE", &
                                        "Use a fixed list of occupations"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="LIST", &
                          description="A list of fractional occupations to use. Must match the number of states "// &
                          "and sum up to the correct number of electrons", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=real_t, &
                          usage="LIST 2.0 0.6666 0.6666 0.66666 0.0 0.0")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ELECTRONIC_TEMPERATURE", &
                          variants=s2a("ELEC_TEMP", "TELEC"), &
                          description="Electronic temperature in the case of Fermi-Dirac smearing", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=300.0_dp, unit_str="K"), &
                          unit_str="K", &
                          usage="ELECTRONIC_TEMPERATURE [K] 300")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS_FERMI_DIRAC", &
                          description="Accuracy checks on occupation numbers use this as a tolerance", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.0E-10_dp, &
                          usage="EPS_FERMI_DIRAC 1.0E-6")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="WINDOW_SIZE", &
                          description="Size of the energy window centred at the Fermi level", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="au_e", &
                          usage="WINDOW_SIZE [eV] 0.3")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIXED_MAGNETIC_MOMENT", &
                          description="Imposed difference between the numbers of electrons of spin up "// &
                          "and spin down: m = n(up) - n(down). A negative value (default) allows "// &
                          "for a change of the magnetic moment. -1 specifically keeps an integer "// &
                          "number of spin up and spin down electrons.", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=-100.0_dp, &
                          usage="FIXED_MAGNETIC_MOMENT 1.5")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_smear_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_rtp_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_rtp_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, print_section

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="REAL_TIME_PROPAGATION", &
                          description="Parameters needed to set up the real time propagation"// &
                          " for the electron dynamics. This currently works only in the NVE ensemble.", &
                          n_keywords=4, n_subsections=4, repeats=.FALSE., &
                          citations=(/Kunert2003, Andermatt2016/))

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximal number of iterations for the self consistent propagator loop.", &
                          usage="MAX_ITER 10", &
                          default_i_val=10)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Convergence criterion for the self consistent propagator loop.", &
                          usage="EPS_ITER 1.0E-5", &
                          default_r_val=1.0E-7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ASPC_ORDER", &
                          description="Speciefies how many steps will be used for extrapolation. "// &
                          "One will be always used which is means X(t+dt)=X(t)", &
                          usage="ASPC_ORDER 3", &
                          default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAT_EXP", &
                          description="Which method should be used to calculate the exponential"// &
                          " in the propagator. It is recommended to use BCH when employing density_propagation "// &
                          "and ARNOLDI otherwise.", &
                          usage="MAT_EXP TAYLOR", default_i_val=do_arnoldi, &
                          enum_c_vals=s2a("TAYLOR", "PADE", "ARNOLDI", "BCH"), &
                          enum_i_vals=(/do_taylor, do_pade, do_arnoldi, do_bch/), &
                          enum_desc=s2a("exponential is evaluated using scaling and squaring in combination"// &
                                        " with a taylor expansion of the exponential.", &
                                        "uses scaling and squaring together with the pade approximation", &
                                        "uses arnoldi subspace algorithm to compute exp(H)*MO directly, can't be used in "// &
                                        "combination with Crank Nicholson or density propagation", &
                                        "Uses a Baker-Campbell-Hausdorff expansion to propagate the density matrix,"// &
                                        " only works for density propagation"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_PROPAGATION", &
                          description="The density matrix is propagated instead of the molecular orbitals. "// &
                          "This can allow a linear scaling simulation. The density matrix is filtered with "// &
                          "the threshold based on the EPS_FILTER keyword from the LS_SCF section", &
                          usage="DENSITY_PROPAGATION .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SC_CHECK_START", &
                          description="Speciefies how many iteration steps will be done without "// &
                          "a check for self consistency. Can save some time in big calculations.", &
                          usage="SC_CHECK_START 3", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXP_ACCURACY", &
                          description="Accuracy for the taylor and pade approximation. "// &
                          "This is only an upper bound bound since the norm used for the guess "// &
                          "is an upper bound for the needed one.", &
                          usage="EXP_ACCURACY 1.0E-6", &
                          default_r_val=1.0E-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROPAGATOR", &
                          description="Which propagator should be used for the orbitals", &
                          usage="PROPAGATOR ETRS", default_i_val=do_etrs, &
                          enum_c_vals=s2a("ETRS", "CN", "EM"), &
                          enum_i_vals=(/do_etrs, do_cn, do_em/), &
                          enum_desc=s2a("enforced time reversible symmetry", &
                                        "Crank Nicholson propagator", &
                                        "Exponential midpoint propagator"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INITIAL_WFN", &
                          description="Controls the initial WFN used for propagation.", &
                          usage="INITIAL_WFN SCF_WFN", default_i_val=use_scf_wfn, &
                          enum_c_vals=s2a("SCF_WFN", "RESTART_WFN", "RT_RESTART"), &
                          enum_i_vals=(/use_scf_wfn, use_restart_wfn, use_rt_restart/), &
                          enum_desc=s2a("An SCF run is performed to get the initial state.", &
                                        "A wavefunction from a previous SCF is propagated. Especially useful,"// &
                                        " if electronic constraints or restraints are used in the previous calculation, "// &
                                        "since these do not work in the rtp scheme.", &
                                        "use the wavefunction of a real time propagation/ehrenfest run"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_DELTA_PULSE", &
                          description="Applies a delta kick to the initial wfn (only RTP for now - the EMD "// &
                          " case is not yet implemented).", &
                          usage="APPLY_DELTA_PULSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COM_NL", &
                          description="Include non local commutator for periodic delta pulse.", &
                          usage="COM_NL", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC", &
                          description="Apply a delta-kick that is compatible with periodic boundary conditions"// &
                          " for any value of DELTA_PULSE_SCALE. Uses perturbation theory for the preparation of"// &
                          " the initial wfn. Note that the pulse is only applied when INITIAL_WFN is set to SCF_WFN,"// &
                          " and not for restarts (RT_RESTART).", &
                          usage="PERIODIC", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_PULSE_DIRECTION", &
                          description="Direction of the applied electric field. The k vector is given as"// &
                          " 2*Pi*[i,j,k]*inv(h_mat), which for PERIODIC .FALSE. yields exp(ikr) periodic with"// &
                          " the unit cell, only if DELTA_PULSE_SCALE is set to unity. For an orthorhombic cell"// &
                          " [1,0,0] yields [2*Pi/L_x,0,0]. For small cells, this results in a very large kick.", &
                          usage="DELTA_PULSE_DIRECTION 1 1 1", n_var=3, default_i_vals=(/1, 0, 0/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_PULSE_SCALE", &
                          description="Scale the k vector, which for PERIODIC .FALSE. results in exp(ikr) no"// &
                          " longer being periodic with the unit cell. The norm of k is the strength of the"// &
                          " applied electric field in atomic units.", &
                          usage="DELTA_PULSE_SCALE 0.01 ", n_var=1, default_r_val=0.001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HFX_BALANCE_IN_CORE", &
                          description="If HFX is used, this keyword forces a redistribution/recalculation"// &
                          " of the integrals, balanced with respect to the in core steps.", &
                          usage="HFX_BALANCE_IN_CORE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MCWEENY_MAX_ITER", &
                          description="Determines the maximum amount of McWeeny steps used after each converged"// &
                          " step in density propagation", &
                          usage="MCWEENY_MAX_ITER 2", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="ACCURACY_REFINEMENT", &
         description="If using density propagation some parts should be calculated with a higher accuracy than the rest"// &
         " to reduce numerical noise. This factor determines by how much the filtering threshold is"// &
         " reduced for these calculations.", &
         usage="ACCURACY_REFINEMENT", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MCWEENY_EPS", &
                          description="Threshold after which McWeeny is terminated", &
                          usage="MCWEENY_EPS 0.00001", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_section)
      CALL section_create(print_section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options for an RTP runs", &
                          repeats=.FALSE.)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing within real time propagation and Eherenfest dynamics", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of the MO restart file during rtp. "// &
                                       "By default keeps a short history of three restarts."// &
                                       "See also RESTART_HISTORY. In density propagation this controls the printing of P.", &
                                       print_level=low_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("MD"), each_iter_values=(/20/), &
                                       add_last=add_last_numeric, filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum index of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=3)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART_HISTORY", &
                                       description="Dumps unique MO restart files during the run keeping all of them. "// &
                                       "In density propagation it dumps the density matrix instead", &
                                       print_level=low_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=(/500/), &
                                       filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum index of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=3)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CURRENT", &
                                       description="Print the current during an EMD simulation to cube files.", &
                                       print_level=high_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=(/20/), &
                                       filename="current")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum index of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, print_section)
      CALL section_release(print_section)

   END SUBROUTINE create_rtp_section

! **************************************************************************************************
!> \brief      Create CP2K input section for the SCCS model
!> \param section ...
!> \par History:
!>      - Creation (10.10.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE create_sccs_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_sccs_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, &
                          name="SCCS", &
                          description="Define the parameters for self-consistent continuum solvation (SCCS) model", &
                          citations=(/Fattebert2002, Andreussi2012, Yin2017/), &
                          n_keywords=8, &
                          n_subsections=2, &
                          repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the SCCS section", &
                          usage="&SCCS ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ALPHA", &
                          description="Solvent specific tunable parameter for the calculation of "// &
                          "the repulsion term <i>G<sup>rep</sup> = &alpha; S</i> "// &
                          "where <i>S</i> is the (quantum) surface of the cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="mN/m")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="BETA", &
                          description="Solvent specific tunable parameter for the calculation of "// &
                          "the dispersion term <i>G<sup>dis</sup> = &beta; V</i> "// &
                          "where <i>V</i> is the (quantum) volume of the cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="GPa")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DELTA_RHO", &
                          description="Numerical increment for the calculation of the (quantum) "// &
                          "surface of the solute cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=2.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DERIVATIVE_METHOD", &
                          description="Method for the calculation of the numerical derivatives on the real-space grids", &
                          usage="DERIVATIVE_METHOD cd5", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_i_val=sccs_derivative_fft, &
                          enum_c_vals=s2a("FFT", "CD3", "CD5", "CD7"), &
                          enum_i_vals=(/sccs_derivative_fft, &
                                        sccs_derivative_cd3, &
                                        sccs_derivative_cd5, &
                                        sccs_derivative_cd7/), &
                          enum_desc=s2a("Fast Fourier transformation", &
                                        "3-point stencil central differences", &
                                        "5-point stencil central differences", &
                                        "7-point stencil central differences"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DIELECTRIC_CONSTANT", &
                          variants=(/"EPSILON_SOLVENT"/), &
                          description="Dielectric constant of the solvent", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=80.0_dp, &
                          usage="DIELECTRIC_CONSTANT 78.36")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS_SCCS", &
                          variants=s2a("EPS_ITER", "TAU_POL"), &
                          description="Tolerance for the convergence of the polarisation density, "// &
                          "i.e. requested accuracy for the SCCS iteration cycle", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.0E-6_dp, &
                          usage="EPS_ITER 1.0E-7")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS_SCF", &
                          description="The SCCS iteration cycle is activated only if the SCF iteration cycle "// &
                          "is converged to this threshold value", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.5_dp, &
                          usage="EPS_SCF 1.0E-2")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="GAMMA", &
                          variants=s2a("SURFACE_TENSION"), &
                          description="Surface tension of the solvent used for the calculation of "// &
                          "the cavitation term <i>G<sup>cav</sup> = &gamma; S</i> "// &
                          "where <i>S</i> is the (quantum) surface of the cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="mN/m")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MAX_ITER", &
                          description="Maximum number of SCCS iteration steps performed to converge "// &
                          "within the given tolerance", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=integer_t, &
                          default_i_val=100, &
                          usage="MAX_ITER 50")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="METHOD", &
                          description="Method used for the smoothing of the dielectric function", &
                          usage="METHOD Fattebert-Gygi", &
                          default_i_val=sccs_andreussi, &
                          enum_c_vals=s2a("ANDREUSSI", "FATTEBERT-GYGI"), &
                          enum_i_vals=(/sccs_andreussi, sccs_fattebert_gygi/), &
                          enum_desc=s2a("Smoothing function proposed by Andreussi et al.", &
                                        "Smoothing function proposed by Fattebert and Gygi"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MIXING", &
                          variants=(/"ETA"/), &
                          description="Mixing parameter (Hartree damping) employed during the iteration procedure", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.6_dp, &
                          usage="MIXING 0.2")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)

      CALL section_create(subsection, __LOCATION__, &
                          name="ANDREUSSI", &
                          description="Define the parameters of the dielectric smoothing function proposed by "// &
                          "Andreussi et al.", &
                          citations=(/Andreussi2012/), &
                          n_keywords=2, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RHO_MAX", &
                          description="Maximum density value used for the smoothing of the dielectric function", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0035_dp, &
                          usage="RHO_MAX 0.01")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RHO_MIN", &
                          description="Minimum density value used for the smoothing of the dielectric function", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0001_dp, &
                          usage="RHO_MIN 0.0003")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, &
                          name="FATTEBERT-GYGI", &
                          description="Define the parameters of the dielectric smoothing function proposed by "// &
                          "Fattebert and Gygi", &
                          citations=(/Fattebert2002/), &
                          n_keywords=2, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="BETA", &
                          description="Parameter &beta; changes the width of the interface solute-solvent", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.7_dp, &
                          usage="BETA 1.3")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RHO_ZERO", &
                          variants=(/"RHO0"/), &
                          description="Parameter &rho;<sub>0</sub> defines the critical density in the middle "// &
                          "of the interface solute-solvent", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0006_dp, &
                          usage="RHO_ZERO 0.0004")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_sccs_section

! **************************************************************************************************
!> \brief      Create CP2K input section for the calculation of an active space Hamiltonian
!> \param section ...
!> \par History:
!>      - Creation 06.04.2016
!> \author     JHU
! **************************************************************************************************
   SUBROUTINE create_active_space_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_active_space_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="ACTIVE_SPACE", &
                          description="Define parameters and method to calculate an electronic active space", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, subsection, print_key)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the ACTIVE_SPACE section", &
                          usage="&ACTIVE_SPACE ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL", &
                          description="Model used to calculate active space interaction Hamiltonian.", &
                          usage="MODEL HARTREE-FOCK", &
                          default_i_val=hf_model, &
                          enum_c_vals=s2a("HARTREE-FOCK", "RSDFT", "DMFT"), &
                          enum_i_vals=(/hf_model, rsdft_model, dmft_model/), &
                          enum_desc=s2a("Hartree-Fock model for interaction Hamiltonian", &
                                        "Range-separated DFT model for interaction Hamiltonian", &
                                        "DMFT model Hamiltonian"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ELECTRONS", &
                          description="The number of active electrons in the CAS space", &
                          usage="ACTIVE_ELECTRONS 4", n_var=1, default_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INACTIVE_ELECTRONS", &
                          description="The number of inactive electrons in the CAS space", &
                          usage="INACTIVE_ELECTRONS 2 {2}", n_var=-1, default_i_vals=(/-1, -1/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ORBITALS", &
                          description="The number of active orbitals defining the CAS space.", &
                          usage="ACTIVE_ORBITALS 2 {2}", n_var=-1, default_i_vals=(/-1, -1/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INACTIVE_ORBITALS", &
                          description="The number of nactive orbitals.", &
                          usage="INACTIVE_ORBITALS 2 {2}", n_var=-1, default_i_vals=(/-1, -1/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ISOLATED_SYSTEM", &
                          description="System is treated without any periodic boundary conditions.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FCIDUMP", &
                                       description="Controls the writing of a file in FCIDUMP format.", &
                                       print_level=high_print_level, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_SELECTION", &
                          description="Method used to select active space orbitals.", &
                          usage="ORBITAL_SELECTION CANONICAL", &
                          default_i_val=casci_canonical, &
                          enum_c_vals=s2a("CANONICAL", "WANNIER_PROJECTION", "MAO"), &
                          enum_i_vals=(/casci_canonical, wannier_projection, mao_projection/), &
                          enum_desc=s2a("Select orbitals using energy ordering of canoncial orbitals", &
                                        "Select orbitals from projected Wannier functions", &
                                        "Select orbitals from modified atomic orbitals"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUBSPACE_ATOM", &
                          description="Number of atom that defines the subspace to be projected on.", &
                          usage="SUBSPACE_ATOM x", default_i_val=-1, &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUBSPACE_SHELL", &
                          description="Shell definition for subsapce.", &
                          usage="SUBSPACE_SHELL 3d4s", default_c_val="X", &
                          type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_print_orb_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eri_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eri_gpw(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_active_space_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_orb_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_orb_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT_ORBITAL_CUBES", &
                          description="Controls printing of active orbital cube files.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FILENAME", &
                          description="Body of Filename for the cube files.", &
                          usage="FILENAME {name}", default_c_val="ActiveOrbital", &
                          type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALIST", &
                          description="List of alpha orbitals to be printed. -1 defaults to all values", &
                          usage="ALIST {1 2 3 ...}", n_var=-1, default_i_vals=(/-1/), &
                          lone_keyword_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BLIST", &
                          description="List of beta orbitals to be printed. -1 defaults to all values", &
                          usage="ALIST {1 2 3 ...}", n_var=-1, default_i_vals=(/-1/), &
                          lone_keyword_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file"// &
                          "(larger values result in smaller cube files)."// &
                          " You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components", &
                          usage="STRIDE {2 2 2}", n_var=-1, default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_print_orb_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_eri_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ERI", &
                          description="Parameters for the electron repulsion integrals.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Method used in ERI calculation.", &
                          usage="METHOD FULL_GPW", &
                          enum_c_vals=s2a("FULL_GPW", "GPW_HALF_TRANSFORM"), &
                          enum_i_vals=(/eri_method_full_gpw, eri_method_gpw_ht/), &
                          enum_desc=s2a("Use the GPW approach with MOs", &
                                        "Use the GPW approach for half-transformed MO ERIs"), &
                          default_i_val=eri_method_gpw_ht)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR", &
                          description="Operator used in ERI calculation.", &
                          usage="OPERATOR <1/R>", &
                          enum_c_vals=s2a("<1/R>", "<EXP(-A*R)/R>", "<ERF(A*R)/R>", &
                                          "<ERFC(A*R)/R>", "<EXP(-A*R2)/R>"), &
                          enum_i_vals=(/eri_operator_coulomb, eri_operator_yukawa, &
                                        eri_operator_erf, eri_operator_erfc, eri_operator_gaussian/), &
                          enum_desc=s2a("Coulomb operator", &
                                        "Yukawa potential operator", &
                                        "Error function potential operator", &
                                        "Complementary error function potential operator", &
                                        "Gaussian potential operator"), &
                          default_i_val=eri_operator_coulomb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR_PARAMETER", &
                          description="Range parameter for ERI operator.", &
                          usage="OPERATOR_PARAMETER 4.0", type_of_var=real_t, &
                          default_r_val=8.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODICITY", &
                          description="Periodicity used for operators in ERI calclulation.", &
                          usage="PERIODICITY {1 1 1}", n_var=-1, default_i_vals=(/1, 1, 1/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="Radius of operator interactions in non-periodic cases.", &
                          usage="CUTOFF_RADIUS 20", type_of_var=real_t, &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_INTEGRAL", &
         description="Accuracy of ERIs that will be stored.", &
         usage="EPS_FILTER 1.0E-10 ", type_of_var=real_t, &
         default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_eri_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_gpw(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_eri_gpw', routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ERI_GPW", &
                          description="Parameters for the GPW approach to electron repulsion integrals.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_GRID", &
                          description="Determines a threshold for the GPW based integration", &
                          usage="EPS_GRID 1.0E-9 ", type_of_var=real_t, &
                          default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="The cutoff of the finest grid level in the GPW integration.", &
                          usage="CUTOFF 300", type_of_var=real_t, &
                          default_r_val=300.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REL_CUTOFF", &
                          variants=(/"RELATIVE_CUTOFF"/), &
                          description="Determines the grid at which a Gaussian is mapped.", &
                          usage="REL_CUTOFF 50", type_of_var=real_t, &
                          default_r_val=50.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STORE_WFN", &
                          variants=(/"STORE_WAVEFUNCTION"/), &
                          description="Strore wavefunction in real space representation for integration.", &
                          usage="STORE_WFN  T", type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRINT_LEVEL", &
                          variants=(/"IOLEVEL"/), &
                          description="How much output is written by the individual groups.", &
                          usage="PRINT_LEVEL HIGH", &
                          default_i_val=silent_print_level, enum_c_vals= &
                          s2a("SILENT", "LOW", "MEDIUM", "HIGH", "DEBUG"), &
                          enum_desc=s2a("Almost no output", &
                                        "Little output", "Quite some output", "Lots of output", &
                                        "Everything is written out, useful for debugging purposes only"), &
                          enum_i_vals=(/silent_print_level, low_print_level, medium_print_level, &
                                        high_print_level, debug_print_level/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_eri_gpw

! **************************************************************************************************
!> \brief Creates the section for cube files related to the implicit Poisson solver.
!> \param section the section to be created
!> \par History
!>      03.2016 refactored from create_print_dft_section [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE create_implicit_psolver_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_implicit_psolver_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="IMPLICIT_PSOLVER", &
                          description="Controls printing of cube files for data from the implicit "// &
                          "(generalized) Poisson solver.", &
                          citations=(/BaniHashemian2016/), &
                          n_keywords=0, n_subsections=3, repeats=.FALSE.)

      NULLIFY (keyword, print_key)

      ! dielectric constant function
      CALL cp_print_key_section_create(print_key, __LOCATION__, "DIELECTRIC_CUBE", &
                                       description="Controls the printing of a cube file with dielectric constant from "// &
                                       "the implicit (generalized) Poisson solver.", &
                                       print_level=high_print_level, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! dirichlet type constraints
      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "DIRICHLET_BC_CUBE", &
         description="Controls the printing of cube files with unit step functions (constraints) "// &
         "representing Dirichlet-type (boundary) regions defined in the implicit (generalized) Poisson "// &
         "solver section. The regions remain unchanged throughout the calculations. If the Dirichlet "// &
         " regions are relatively large and/or the number of partitions is quite high, in order to save memory,"// &
         " generate the cube files in early steps and perform the rest of the calculations with this keyword"// &
         " switched off.", &
         print_level=high_print_level, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="TILE_CUBES", &
                          description="Print tiles that tessellate the Dirichlet regions into cube files. If TRUE, "// &
                          "generates cube files as many as the total number of tiles.", &
                          usage="TILE_CUBES <logical>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! charge introduced by Lagrange multipliers
      CALL cp_print_key_section_create(print_key, __LOCATION__, "DIRICHLET_CSTR_CHARGE_CUBE", &
                                       description="Controls the printing of cube files with penalty charges induced to "// &
                                       "Dirichlet regions by Lagrange multipliers (implicit Poisson solver).", &
                                       print_level=high_print_level, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_implicit_psolver_section

! **************************************************************************************************
!> \brief creates the input section for parameters related to CDFT specific optimizers
!> \param section the section to be created
!> \par History
!>      03.2018 separated from create_outer_scf_section [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE create_cdft_opt_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_cdft_opt_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CDFT_OPT", &
                          description="Parameters controlling optimization methods that are compatible "// &
                          "only with CDFT based constraints (i.e. CDFT SCF is active). Specifically, "// &
                          "the control parameters for the Broyden and Newton optimizers are defined in this "// &
                          "section.", &
                          n_keywords=10, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_TYPE", &
                          description="Specifies the Broyden optimizer variant to use.", &
                          usage="BROYDEN_TYPE BT1", &
                          default_i_val=broyden_type_1, &
                          enum_c_vals=s2a("BT1", "BT1_EXPLICIT", "BT2", "BT2_EXPLICIT", &
                                          "BT1_LS", "BT1_EXPLICIT_LS", "BT2_LS", "BT2_EXPLICIT_LS"), &
                          enum_desc=s2a("Broyden's first method, also known as the good method. The initial Jacobian"// &
                                        " is built from MD history if available. Otherwise switches to SD for one"// &
                                        " SCF iteration until a Jacobian can be built from the SCF history.", &
                                        "Same as BT1, but computes the explicit Jacobian with finite differences. "// &
                                        "Requires a CDFT SCF procedure to be active.", &
                                        "Same as BT1, but uses Broyden's second method, also known as the bad method.", &
                                        "Same as BT1_EXPLICIT, but using Broyden's second method.", &
                                        "Same as BT1, but uses backtracking line search for optimizing the step size "// &
                                        "(see optimizer NEWTON_LS).", &
                                        "Same as BT1_EXPLICIT, but uses backtracking line search for optimizing the step size.", &
                                        "Same as BT2, but uses backtracking line search for optimizing the step size.", &
                                        "Same as BT2_EXPLICIT, but uses backtracking line search for optimizing the step size."), &
                          enum_i_vals=(/broyden_type_1, broyden_type_1_explicit, broyden_type_2, &
                                        broyden_type_2_explicit, broyden_type_1_ls, broyden_type_1_explicit_ls, &
                                        broyden_type_2_ls, broyden_type_2_explicit_ls/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_TYPE", &
                          description="Finite difference method used to calculate the inverse Jacobian "// &
                          "needed by some optimizers. Compatible only with CDFT constraints.", &
                          usage="JACOBIAN_TYPE FD1", &
                          default_i_val=jacobian_fd1, &
                          enum_c_vals=s2a("FD1", "FD1_BACKWARD", "FD2", "FD2_BACKWARD", "FD1_CENTRAL"), &
                          enum_desc=s2a("First order forward difference (one extra energy evaluation per constraint).", &
                                        "First order backward difference (one extra energy evaluation per constraint).", &
                                        "Second order forward difference (two extra energy evaluations per constraint).", &
                                        "Second order backward difference (two extra energy evaluations per constraint).", &
                                        "First order central difference (two extra energy evaluations per constraint)."), &
                          enum_i_vals=(/jacobian_fd1, jacobian_fd1_backward, jacobian_fd2, &
                                        jacobian_fd2_backward, jacobian_fd1_central/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_STEP", &
                          description="Step size to use in the calculation of the inverse Jacobian with finite differences. "// &
                          "Expects one value for all constraints, or one value per constraint.", &
                          usage="JACOBIAN_STEP 5.0E-3 ", n_var=-1, default_r_val=5.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_FREQ", &
                          description="Defines parameters that control how often the explicit Jacobian is built, "// &
                          "which is needed by some optimizers. Expects two values. The first value "// &
                          "determines how many consecutive CDFT SCF iterations should skip a rebuild,"// &
                          " whereas the latter how many MD steps. The values can be zero (meaning never "// &
                          " rebuild) or positive. Both values cannot be zero.", &
                          usage="JACOBIAN_FREQ 1 1", n_var=2, &
                          default_i_vals=(/1, 1/), type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_RESTART", &
                          description="Restart the inverse Jacobian using the vector defined with keyword JACOBIAN_VECTOR.", &
                          usage="JACOBIAN_RESTART TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_VECTOR", &
                          description="Defines the inverse Jacobian matrix. Useful for restarting calculations. "// &
                          "Expects n^2 values where n is the total number of constraints. "// &
                          "The matrix should be given in row major order.", &
                          usage="JACOBIAN_VECTOR 1.0 0.0", n_var=-1, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_LS", &
                          description="The maximum number of backtracking line search steps to perform.", &
                          usage="MAX_LS 5", default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FACTOR_LS", &
                          description="Control parameter for backtracking line search. The step size is reduced by "// &
                          "this factor on every line search iteration. Value must be between 0 and 1 (exclusive).", &
                          usage="FACTOR_LS 0.5", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONTINUE_LS", &
                          description="Continue backtracking line search until MAX_LS steps are reached or the "// &
                          "norm of the CDFT gradient no longer decreases. Default (false) behavior exits the "// &
                          "line search procedure on the first step that the gradient decreases.", &
                          usage="CONTINUE_LS TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_cdft_opt_section

END MODULE input_cp2k_dft
