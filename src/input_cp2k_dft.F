!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
        Andermatt2016, Andreussi2012, Avezac2005, Bengtsson1999, Blochl1995, Brelaz1979, &
        Dewar1977, Dewar1985, Fattebert2002, Golze2017a, Golze2017b, Guidon2010, Iannuzzi2005, &
        Iannuzzi2006, Iannuzzi2007, Kolafa2004, Krack2000, Kuhne2007, Kunert2003, Lippert1997, &
        Lippert1999, Merlot2014, Perdew1981, Repasky2002, Rocha2006, Schenter2008, Shigeta2001, &
        Stewart1989, Stewart2007, Thiel1992, VanVoorhis2015, VandeVondele2005a, VandeVondele2005b, &
        VandeVondele2006, Yin2017
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
        admm1_type, admm2_type, admmp_type, admmq_type, admms_type, casci_canonical, dmft_model, &
        do_admm_aux_exch_func_bee, do_admm_aux_exch_func_bee_libxc, do_admm_aux_exch_func_default, &
        do_admm_aux_exch_func_default_libxc, do_admm_aux_exch_func_none, &
        do_admm_aux_exch_func_opt, do_admm_aux_exch_func_opt_libxc, do_admm_aux_exch_func_pbex, &
        do_admm_aux_exch_func_pbex_libxc, do_admm_aux_exch_func_sx_libxc, &
        do_admm_basis_projection, do_admm_blocked_projection, do_admm_blocking_purify_full, &
        do_admm_charge_constrained_projection, do_admm_exch_scaling_merlot, &
        do_admm_exch_scaling_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace, &
        do_admm_purify_mcweeny, do_admm_purify_mo_diag, do_admm_purify_mo_no_diag, &
        do_admm_purify_none, do_admm_purify_none_dm, do_arnoldi, do_bch, do_cn, &
        do_ddapc_constraint, do_ddapc_restraint, do_em, do_etrs, do_full_density, do_gapw_gcs, &
        do_gapw_gct, do_gapw_log, do_lri_inv, do_lri_inv_auto, do_lri_opt_all, do_lri_opt_coeff, &
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
        eri_operator_trunc, eri_operator_yukawa, gapw_1c_large, gapw_1c_medium, gapw_1c_orb, &
        gapw_1c_small, gapw_1c_very_large, gaussian, hf_model, kg_color_dsatur, kg_color_greedy, &
        kg_tnadd_atomic, kg_tnadd_embed, kg_tnadd_embed_ri, kg_tnadd_none, manual_selection, &
        mao_projection, no_admm_type, no_excitations, no_solver, numerical, oe_gllb, oe_lb, &
        oe_none, oe_saop, oe_sic, ot_mini_cg, ot_mini_diis, plus_u_lowdin, plus_u_mulliken, &
        plus_u_mulliken_charges, qiskit_solver, real_time_propagation, rel_dkh, rel_none, &
        rel_pot_erfc, rel_pot_full, rel_sczora_mp, rel_trans_atom, rel_trans_full, &
        rel_trans_molecule, rel_zora, rel_zora_full, rel_zora_mp, rsdft_model, sccs_andreussi, &
        sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7, sccs_derivative_fft, &
        sccs_fattebert_gygi, sic_ad, sic_eo, sic_list_all, sic_list_unpaired, sic_mauri_spz, &
        sic_mauri_us, sic_none, slater, tddfpt_davidson, tddfpt_excitations, tddfpt_lanczos, &
        tddfpt_singlet, tddfpt_spin_cons, tddfpt_spin_flip, tddfpt_triplet, use_mom_ref_coac, &
        use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero, use_restart_wfn, use_rt_restart, &
        use_scf_wfn, wannier_projection, weight_type_mass, weight_type_unit, wfi_aspc_nr, &
        wfi_frozen_method_nr, wfi_linear_p_method_nr, wfi_linear_ps_method_nr, &
        wfi_linear_wf_method_nr, wfi_ps_method_nr, wfi_use_guess_method_nr, &
        wfi_use_prev_p_method_nr, wfi_use_prev_rho_r_method_nr, wfi_use_prev_wf_method_nr, &
        xas_1s_type, xas_2p_type, xas_2s_type, xas_3d_type, xas_3p_type, xas_3s_type, xas_4d_type, &
        xas_4f_type, xas_4p_type, xas_4s_type, xas_dip_len, xas_dip_vel, xas_dscf, xas_none, &
        xas_not_excited, xas_tdp_by_index, xas_tdp_by_kind, xas_tp_fh, xas_tp_flex, xas_tp_hh, &
        xas_tp_xfh, xas_tp_xhh, xes_tp_val
   USE input_cp2k_almo,                 ONLY: create_almo_scf_section
   USE input_cp2k_distribution,         ONLY: create_distribution_section
   USE input_cp2k_ec,                   ONLY: create_ec_section
   USE input_cp2k_exstate,              ONLY: create_exstate_section
   USE input_cp2k_external,             ONLY: create_ext_den_section,&
                                              create_ext_pot_section,&
                                              create_ext_vxc_section
   USE input_cp2k_field,                ONLY: create_efield_section,&
                                              create_per_efield_section
   USE input_cp2k_kpoints,              ONLY: create_kpoints_section
   USE input_cp2k_loc,                  ONLY: create_localize_section,&
                                              print_wanniers
   USE input_cp2k_ls,                   ONLY: create_ls_scf_section
   USE input_cp2k_mm,                   ONLY: create_neighbor_lists_section
   USE input_cp2k_poisson,              ONLY: create_poisson_section
   USE input_cp2k_print_dft,            ONLY: create_pdos_section,&
                                              create_print_dft_section
   USE input_cp2k_projection_rtp,       ONLY: create_projection_rtp_section
   USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
   USE input_cp2k_scf,                  ONLY: create_cdft_control_section,&
                                              create_scf_section
   USE input_cp2k_tb,                   ONLY: create_dftb_control_section,&
                                              create_xtb_control_section
   USE input_cp2k_transport,            ONLY: create_transport_section
   USE input_cp2k_xc,                   ONLY: create_xc_fun_section,&
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
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_dft'

   PUBLIC :: create_dft_section
   PUBLIC :: create_bsse_section, create_qs_section
   PUBLIC :: create_interp_section
   PUBLIC :: create_ddapc_restraint_section
   PUBLIC :: create_mgrid_section
   PUBLIC :: create_lrigpw_section

CONTAINS

! **************************************************************************************************
!> \brief creates the dft section
!> \param section the section to be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_dft_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          description="Specify size of automatically generated auxiliary (RI) basis sets: "// &
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
                          "the electrostatic potential with "// &
                          "by compensating for the surface dipole. Implemented only for slabs with normal "// &
                          "parallel to one Cartesian axis. The normal direction is given by the keyword SURF_DIP_DIR", &
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
                          description="This keyword assigns an user defined position in Angstroms "// &
                          "in the direction normal to the surface (given by SURF_DIP_DIR). "// &
                          "The default value is -1.0_dp which appplies the correction at a position "// &
                          "that has minimum electron density on the grid.", &
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

      CALL create_exstate_section(subsection)
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

      CALL create_active_space_section(subsection)
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
!> \brief creates the input section for the qs part
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_qs_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          " Overrides EPS_RHO value", &
                          usage="EPS_RHO_RSPACE real", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RHO_GSPACE", &
                          description="Sets precision of the density mapping in gspace.Overrides EPS_DEFAULT value."// &
                          " Overrides EPS_RHO value", &
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
         "the local densities are computed even if hard and soft should be equal. "// &
         "If this keyword is not set to true, those atoms with soft basis sets are treated by a GPW scheme, i.e. "// &
         "the corresponding density contribution goes on the global grid and is expanded in PW. "// &
         "This option nullifies the effect of the GPW_TYPE in the atomic KIND", &
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

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_1C_BASIS", &
                          description="Specifies how to construct the GAPW one center basis set. "// &
                          "Default is to use the primitives from the orbital basis.", &
                          usage="GAPW_1C_BASIS MEDIUM", &
                          enum_c_vals=s2a("ORB", "EXT_SMALL", "EXT_MEDIUM", "EXT_LARGE", "EXT_VERY_LARGE"), &
                          enum_desc=s2a("Use orbital basis set.", &
                                        "Extension using Small number of primitive Gaussians.", &
                                        "Extension using Medium number of primitive Gaussians.", &
                                        "Extension using Large number of primitive Gaussians.", &
                                        "Extension using Very Large number of primitive Gaussians."), &
                          enum_i_vals=(/gapw_1c_orb, gapw_1c_small, gapw_1c_medium, &
                                        gapw_1c_large, gapw_1c_very_large/), &
                          default_i_val=gapw_1c_orb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_PAIR_LIST_RADIUS", &
                          description="Set the minimum value [Bohr] for the overlap pair list radius."// &
                          " Default is 0.0 Bohr, negative values are changed to the cell size."// &
                          " This allows to control the sparsity of the KS matrix for HFX calculations.", &
                          usage="MIN_PAIR_LIST_RADIUS real", default_r_val=0.0_dp)
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

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="AUXILIARY_DENSITY_MATRIX_METHOD", &
                          description="Parameters needed for the ADMM method.", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=(/Guidon2010/))

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ADMM_TYPE", &
         description="Type of ADMM (sort name) as refered in literature. "// &
         "This sets values for METHOD, ADMM_PURIFICATION_METHOD, and EXCH_SCALING_MODEL", &
         enum_c_vals=s2a("NONE", "ADMM1", "ADMM2", "ADMMS", "ADMMP", "ADMMQ"), &
         enum_desc=s2a("No short name is used, use specific definitions (default)", &
                       "ADMM1 method from Guidon2010", &
                       "ADMM2 method from Guidon2010", &
                       "ADMMS method from Merlot2014", &
                       "ADMMP method from Merlot2014", &
                       "ADMMQ method from Merlot2014"), &
         enum_i_vals=(/no_admm_type, admm1_type, admm2_type, admms_type, admmp_type, admmq_type/), &
         default_i_val=no_admm_type, &
         citations=(/Guidon2010, Merlot2014/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

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
         description="Exchange functional which is used for the ADMM correction. "// &
         "LibXC implementations require linking with LibXC", &
         enum_c_vals=s2a("DEFAULT", "PBEX", "NONE", "OPTX", "BECKE88X", &
                         "PBEX_LIBXC", "BECKE88X_LIBXC", "OPTX_LIBXC", "DEFAULT_LIBXC", "LDA_X_LIBXC"), &
         enum_i_vals=(/do_admm_aux_exch_func_default, do_admm_aux_exch_func_pbex, &
                       do_admm_aux_exch_func_none, do_admm_aux_exch_func_opt, do_admm_aux_exch_func_bee, &
                       do_admm_aux_exch_func_pbex_libxc, do_admm_aux_exch_func_bee_libxc, &
                       do_admm_aux_exch_func_opt_libxc, do_admm_aux_exch_func_default_libxc, &
                       do_admm_aux_exch_func_sx_libxc/), &
         enum_desc=s2a("Use PBE-based corrections according to the chosen interaction operator.", &
                       "Use PBEX functional for exchange correction.", &
                       "No correction: X(D)-x(d)-> 0.", &
                       "Use OPTX functional for exchange correction.", &
                       "Use Becke88X functional for exchange correction.", &
                       "Use PBEX functional (LibXC implementation) for exchange correction.", &
                       "Use Becke88X functional (LibXC implementation) for exchange correction.", &
                       "Use OPTX functional (LibXC implementation) for exchange correction.", &
                       "Use PBE-based corrections (LibXC where possible) to the chosen interaction operator.", &
                       "Use Slater X functional (LibXC where possible) for exchange correction."), &
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

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXCHANGE", &
                          description="Setup parameters for the evaluation of the EXCHANGE and "// &
                          "core Hamiltonian terms in SE calculations.", n_keywords=0, n_subsections=1, &
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
                                       description="Activates the printing of the subcells used for the "// &
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

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (keyword, print_key)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DENSITY_FITTING", &
                          description="Setup parameters for density fitting (Bloechl charges or density derived "// &
                          "atomic point charges (DDAPC) charges)", &
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
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_s2_restraint_section(section)
      TYPE(section_type), POINTER                        :: section

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

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="tddfpt", &
                          description="Old TDDFPT code. Use new version in CP2K_INPUT / FORCE_EVAL / PROPERTIES / TDDFPT", &
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
                          description="Which type of orbital eigenvalue correction to use "// &
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
!> \brief creates the KG section
!> \param section ...
!> \author Martin Haeufel [2012.07]
! **************************************************************************************************
   SUBROUTINE create_kg_section(section)
      TYPE(section_type), POINTER                        :: section

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
!> \brief Create the BSSE section for counterpoise correction
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_bsse_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          description="Specify additional parameters for the combinatorial configurations. "// &
                          "Use this section to manually specify charge and multiplicity of the fragments "// &
                          "and their combinations.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="GLB_CONF", &
                          description="Specifies the global configuration using 1 or 0 for each fragment. "// &
                          "1 specifies the respective fragment as used, 0 as unused.", &
                          usage="GLB_CONF {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUB_CONF", &
                          description="Specifies the subconfiguration using 1 or 0 belonging to the global configuration. "// &
                          "1 specifies the respective fragment as real, 0 as ghost.", &
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
                                       description="Controls the dumping of the restart file during BSSE runs. "// &
                                       "By default the restart is updated after each configuration calculation is "// &
                                       "completed.", &
                                       print_level=silent_print_level, common_iter_levels=0, &
                                       add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_bsse_section

! **************************************************************************************************
!> \brief input section for optional parameters for LRIGPW
!>        LRI: local resolution of identity
!> \param section the section to create
!> \author Dorothea Golze [02.2015]
! **************************************************************************************************
   SUBROUTINE create_lrigpw_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          description="Method used to separate pair density for distant pairs. "// &
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
                          description="Read the matrix embedding potential "// &
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

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_DENS_DIFF", &
                                       description="Controls the printing of cube files with "// &
                                       "embedding densisty differences", &
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
                          " A value 50+-10Ry might be required for highly accurate results,"// &
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
!> \brief makes the input section for core-level spectroscopy simulations
!> \param section ...
!> \par History
!>      03.2005 created [MI]
! **************************************************************************************************
   SUBROUTINE create_xas_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          description="Total number of electrons for spin channel alpha, in XAS calculation by TP_FLEX. "// &
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
                          "This keyword can be repeated several times "// &
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
                          "This keyword can be repeated several times "// &
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
                                       description="Controls the dumping of MO restart file during the SCF. "// &
                                       "of a Core-Level-Spectroscopy calculation. For each new excited atom, "// &
                                       "one different restart file is dumped. These restart files should be "// &
                                       "employed only to restart the same type of CLS calculation, "// &
                                       "i.e. with the same core potential.", &
                                       print_level=low_print_level, common_iter_levels=3, each_iter_names=s2a("XAS_SCF"), &
                                       add_last=add_last_numeric, each_iter_values=(/3/), filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FULL_RESTART", &
                                       description="Controls the dumping of a standard MO restart file "// &
                                       "where coefficients and occupation numbers are those of the TP scheme, "// &
                                       "i.e. with emptied core state.", &
                                       print_level=high_print_level, common_iter_levels=3, each_iter_names=s2a("XAS_SCF"), &
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
                          description="Indexes of the states to be printed as cube files "// &
                          "This keyword can be repeated several times "// &
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

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection, subsubsection, &
                                                            subsubsubsection

      NULLIFY (keyword, print_key, subsection, subsubsection, subsubsubsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="XAS_TDP", &
                          description="XAS simulations using linear-response TDDFT. Excitation from "// &
                          "specified core orbitals are considered one at a time. In case of high "// &
                          "symmetry structures, donor core orbitals should be localized.", &
                          n_keywords=19, n_subsections=4, repeats=.FALSE.)

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

      CALL keyword_create(keyword, __LOCATION__, name="XYZ_DIPOLE", &
                          variants=s2a("DIPOLE_XYZ"), &
                          description="Whether the detailed contributions of the dipole oscillator "// &
                          "strengths along the X,Y,Z directions should be printed.", &
                          usage="XYZ_DIPOLE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!  the GW2X correction subsection
      CALL section_create(subsection, __LOCATION__, name="GW2X", &
                          description="Specifications for the GW2X calculation of core "// &
                          "ionization potentials. Note that donor states need to be actively "// &
                          "localized using the LOCALIZE keyword in DONOR_STATES. N_SEARCH "// &
                          "should be kept to the minimum, such that only core states are localized.", &
                          citations=(/Shigeta2001/), &
                          n_keywords=8, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the GW2X correction of the core ionization potentials", &
                          usage="&GW2X {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XPS_ONLY", &
                          description="If set to .TRUE., only run GW2X calculations for XPS "// &
                          "spectroscopy and ignore all XAS calculations. It is still "// &
                          "required to define the DONOR_STATES and KERNEL%EXACT_EXCHANGE "// &
                          "subsections.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BATCH_SIZE", &
                          description="MOs batch size for batched tensor contraction. Larger "// &
                          "size is faster, but uses more memory. Default should be safe "// &
                          "in most cases.", &
                          default_i_val=64)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_GW2X", &
                          description="Convergence threshold for GW2X iterations (in eV)", &
                          default_r_val=1.E-2_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="C_OS", &
                          description="Opposite-spin scling factor. SCS => 6/5, SOS => 1.3", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="C_SS", &
                          description="Same-spin scling factor. SCS => 1/3, SOS => 0.0", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_GW2X_ITER", &
                          description="Maximum number of iterations for GW2X", &
                          default_i_val=10)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PSEUDO_CANONICAL", &
                          variants=(/"PSEUDO_CANO"/), &
                          description="Whether the pseudo-canonical version of GW2X should be used "// &
                          "(versus only using the diagonal of the generalized Fock matrix)", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

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
                          description="Whether the N_SEARCH potential donor states should be "// &
                          "actively localized. Necessary in case of excited atoms "// &
                          "equivalent under symmetry or GW2X correction.", &
                          usage="LOCALIZE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of the donor states subsection

! The OT solver subsection
      CALL section_create(subsection, __LOCATION__, name="OT_SOLVER", &
                          description="Specifications for the iterative OT solver. Note: only "// &
                          "available within the Tamm-Dancoff approximation. Recommanded if excitations "// &
                          "from multiple donor states take place.", &
                          n_keywords=4, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the usage of the OT iterator solver", &
                          usage="&OT_SOLVER {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations allowed for the OT solver", &
                          usage="MAX_ITER {integer}", &
                          default_i_val=50)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Convergence threshold for the OT solver", &
                          usage="EPS_ITER {double}", &
                          default_r_val=1.0E-4_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MINIMIZER", &
                          description="Minimizer to be used with the OT solver", &
                          usage="MINIMIZER DIIS", &
                          default_i_val=ot_mini_diis, &
                          enum_c_vals=s2a("CG", "DIIS"), &
                          enum_desc=s2a("Conjugated gradient: safer", &
                                        "Direct inversion of the iterative subspace: faster"), &
                          enum_i_vals=(/ot_mini_cg, ot_mini_diis/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of the OT solver subsection

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_ORBIT_COUPLING", &
                          variants=(/"SOC"/), &
                          description="Whether spin-orbit coupling should be added. "// &
                          "Note: only applies for spin-restricted calculations with "// &
                          "singlet and triplet excitations OR spin-unrestricted "// &
                          "calculations with both spin-conserving and spin-flip.", &
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
                          "grids are used for the xc-kernel integration. "// &
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
      CALL create_xc_fun_section(subsubsection)
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
                          "expensive (especially for large basis sets ). "// &
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

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="RESTART_WFN", &
                                       description="Controles the dumping of a MO restart file for a given "// &
                                       "excited state index. Only for K-edge RKS calculations. "// &
                                       "Can be repeated to get multiple *.wfn files at once.", &
                                       print_level=debug_print_level, filename="", common_iter_levels=1)
      CALL keyword_create(keyword, __LOCATION__, name="EXCITED_STATE_INDEX", variants=(/"INDEX"/), &
                          description="The index of the excited state that should be dumped", &
                          usage="INDEX {int}", default_i_val=1, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

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
                          description="Indexes of the excited states to be printed as cube files "// &
                          "This keyword can be repeated several times "// &
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
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_rtp_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          description="Controls the initial WFN used for propagation. "// &
                          "Note that some energy contributions may not be "// &
                          "initialized in the restart cases, for instance "// &
                          "electronic entropy energy in the case of smearing.", &
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

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_WFN_MIX_INIT_RESTART", &
                          description="If set to True and in the case of INITIAL_WFN=RESTART_WFN, call the "// &
                          "DFT%PRINT%WFN_MIX section to mix the read initial wfn. The starting wave-function of the "// &
                          "RTP will be the mixed one. Setting this to True without a defined WFN_MIX section will "// &
                          "not do anything as defining a WFN_MIX section without this keyword for RTP run with "// &
                          "INITIAL_WFN=RESTART_WFN. Note that if INITIAL_WFN=SCF_WFN, this keyword is not needed to "// &
                          "apply the mixing defined in the WFN_MIX section. Default is False.", &
                          usage="APPLY_WFN_MIX_INIT_RESTART", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_DELTA_PULSE", &
                          description="Applies a delta kick to the initial wfn (only RTP for now - the EMD"// &
                          " case is not yet implemented). Only work for INITIAL_WFN=SCF_WFN", &
                          usage="APPLY_DELTA_PULSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_DELTA_PULSE_MAG", &
                          description="Applies a magnetic delta kick to the initial wfn (only RTP for now - the EMD"// &
                          " case is not yet implemented). Only work for INITIAL_WFN=SCF_WFN", &
                          usage="APPLY_DELTA_PULSE_MAG", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VELOCITY_GAUGE", &
                          description="Perform propagation in the velocity gauge using the explicit vector potential"// &
                          " only a constant vector potential as of now (corresonding to a delta-pulse)."// &
                          " uses DELTA_PULSE_SCALE and DELTA_PULSE_DIRECTION to define the vector potential", &
                          usage="VELOCITY_GAUGE T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUGE_ORIG", &
                          description="Define gauge origin for magnetic perturbation", &
                          usage="GAUGE_ORIG COM", &
                          enum_c_vals=s2a("COM", "COAC", "USER_DEFINED", "ZERO"), &
                          enum_desc=s2a("Use Center of Mass", &
                                        "Use Center of Atomic Charges", &
                                        "Use User Defined Point (Keyword:REF_POINT)", &
                                        "Use Origin of Coordinate System"), &
                          enum_i_vals=(/use_mom_ref_com, &
                                        use_mom_ref_coac, &
                                        use_mom_ref_user, &
                                        use_mom_ref_zero/), &
                          default_i_val=use_mom_ref_com)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUGE_ORIG_MANUAL", &
                          description="Manually defined gauge origin for magnetic perturbation [in Bohr!]", &
                          usage="GAUGE_ORIG_MANUAL x y z", &
                          repeats=.FALSE., &
                          n_var=3, default_r_vals=(/0._dp, 0._dp, 0._dp/), &
                          type_of_var=real_t, &
                          unit_str='bohr')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VG_COM_NL", &
                          description="apply gauge transformed non-local potential term"// &
                          " only affects VELOCITY_GAUGE=.TRUE.", &
                          usage="VG_COM_NL T", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COM_NL", &
                          description="Include non-local commutator for periodic delta pulse."// &
                          " only affects PERIODIC=.TRUE.", &
                          usage="COM_NL", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LEN_REP", &
                          description="Use length representation delta pulse (in conjunction with PERIODIC T)."// &
                          " This corresponds to a 1st order perturbation in the length gauge."// &
                          " Note that this is NOT compatible with a periodic calculation!"// &
                          " Uses the reference point defined in DFT%PRINT%MOMENTS ", &
                          usage="LEN_REP T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC", &
                          description="Apply a delta-kick that is compatible with periodic boundary conditions"// &
                          " for any value of DELTA_PULSE_SCALE. Uses perturbation theory for the preparation of"// &
                          " the initial wfn with the velocity operator as perturbation."// &
                          " If LEN_REP is .FALSE. this corresponds to a first order velocity gauge."// &
                          " Note that the pulse is only applied when INITIAL_WFN is set to SCF_WFN,"// &
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
                                       "By default keeps a short history of three restarts. "// &
                                       "See also RESTART_HISTORY. In density propagation this controls the printing of P.", &
                                       print_level=low_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("MD"), each_iter_values=(/20/), &
                                       add_last=add_last_numeric, filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
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
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FIELD", &
                                       description="Print the time-dependent field applied during an EMD simulation in "// &
                                       "atomic unit.", &
                                       print_level=high_print_level, common_iter_levels=-1, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=(/1/), &
                                       filename="applied_field")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL create_projection_rtp_section(print_key)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CURRENT", &
                                       description="Print the current during an EMD simulation to cube files.", &
                                       print_level=high_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=(/20/), &
                                       filename="current")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
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

      CALL cp_print_key_section_create(print_key, __LOCATION__, "E_CONSTITUENTS", &
                                       description="Print the energy constituents (relevant to RTP) which make up "// &
                                       "the Total Energy", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=(/1/), &
                                       filename="rtp")
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
                          "the repulsion term $G^\text{rep} = \alpha S$ "// &
                          "where $S$ is the (quantum) surface of the cavity", &
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
                          "the dispersion term $G^\text{dis} = \beta V$ "// &
                          "where $V$ is the (quantum) volume of the cavity", &
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
                          name="RELATIVE_PERMITTIVITY", &
                          variants=s2a("DIELECTRIC_CONSTANT", "EPSILON_RELATIVE", "EPSILON_SOLVENT"), &
                          description="Relative permittivity (dielectric constant) of the solvent (medium)", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=80.0_dp, &
                          usage="RELATIVE_PERMITTIVITY 78.36")
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
                          "the cavitation term $G^\text{cav} = \gamma S$ "// &
                          "where $S$ is the (quantum) surface of the cavity", &
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
                          description="Parameter $\rho_0$ defines the critical density in the middle "// &
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

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ORBITALS", &
                          description="The number of active orbitals defining the CAS space.", &
                          usage="ACTIVE_ORBITALS 2", n_var=1, default_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ORBITAL_INDICES", &
                          description="The indices of the active orbitals. Requires ORBITAL_SELECTION MANUAL!", &
                          usage="ACTIVE_ORBITAL_INDICES 2 3 {...}", n_var=-1, default_i_vals=(/-1/), &
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
                          enum_c_vals=s2a("CANONICAL", "WANNIER_PROJECTION", "MAO", "MANUAL"), &
                          enum_i_vals=(/casci_canonical, wannier_projection, mao_projection, manual_selection/), &
                          enum_desc=s2a("Select orbitals using energy ordering of canoncial orbitals", &
                                        "Select orbitals from projected Wannier functions", &
                                        "Select orbitals from modified atomic orbitals", &
                                        "Select orbitals manually via ACTIVE_ORBITAL_INDICES"))

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

      CALL keyword_create(keyword, __LOCATION__, name="SCF_EMBEDDING", &
                          description="Whether to turn on the self-consistent embedding scheme", &
                          default_l_val=.FALSE., lone_keyword_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="QCSCHEMA", &
                          description="Name of the QCSchema file, may include a path", &
                          usage="QCSCHEMA <FILENAME>", &
                          type_of_var=lchar_t, repeats=.FALSE., &
                          default_lc_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AS_SOLVER", &
                          description="The external active space solver for the embedding approach", &
                          usage="AS_SOLVER QISKIT", &
                          default_i_val=no_solver, &
                          enum_c_vals=s2a("NONE", "QISKIT"), &
                          enum_i_vals=(/no_solver, qiskit_solver/), &
                          enum_desc=s2a("NO solver, used to produce FCIDUMP/QCSchema files", &
                                        "QISKIT active space solver"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Energy convergence threshold of the DFT embedding scheme.", &
                          usage="EPS_ITER 1.0E-6 ", type_of_var=real_t, &
                          default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Max number of iterations for the DFT embedding scheme.", &
                          usage="MAX_ITER 50", type_of_var=integer_t, &
                          default_i_val=50)
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

      CALL create_socket_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_active_space_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_socket_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SOCKET", &
                          description="Parameters to set up the socket communicating to the external active space solver.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="INET", &
                          description="Use an INET socket rather than a UNIX socket.", &
                          usage="INET <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PORT", &
                          description="Port number for the socket client.", &
                          usage="port <INTEGER>", &
                          default_i_val=12345)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HOST", &
                          description="Host name for the socket client.", &
                          usage="host <HOSTNAME>", &
                          default_c_val="embedding_socket")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_socket_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_orb_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          " (larger values result in smaller cube files)."// &
                          " You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components", &
                          usage="STRIDE {2 2 2}", n_var=-1, default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STOP_AFTER_CUBES", &
                          description="Whether to stop the computation after printing the cubes.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_print_orb_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          default_i_val=eri_method_full_gpw)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR", &
                          description="Operator used in ERI calculation.", &
                          usage="OPERATOR <1/R>", &
                          enum_c_vals=s2a("<1/R>", "<EXP(-A*R)/R>", "<ERF(A*R)/R>", &
                                          "<ERFC(A*R)/R>", "<EXP(-A*R2)/R>", "<H(A-R)/R>"), &
                          enum_i_vals=(/eri_operator_coulomb, eri_operator_yukawa, &
                                        eri_operator_erf, eri_operator_erfc, eri_operator_gaussian, eri_operator_trunc/), &
                          enum_desc=s2a("Coulomb operator", &
                                        "Yukawa potential operator", &
                                        "Error function potential operator", &
                                        "Complementary error function potential operator", &
                                        "Gaussian potential operator", &
                                        "Truncated Coulomb potential"), &
                          default_i_val=eri_operator_coulomb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR_PARAMETER", &
                          description="Range parameter for ERI operator.", &
                          usage="OPERATOR_PARAMETER 4.0", type_of_var=real_t, &
                          default_r_val=0.5_dp)
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
         usage="EPS_INTEGRAL 1.0E-10 ", type_of_var=real_t, &
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

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Determines a threshold for the sparse matrix multiplications if METHOD "// &
                          "GPW_HALF_TRANSFORM is used", &
                          usage="EPS_FILTER 1.0E-9 ", type_of_var=real_t, &
                          default_r_val=1.0E-9_dp)
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
                          description="Store wavefunction in real space representation for integration.", &
                          usage="STORE_WFN  T", type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GROUP_SIZE", &
                          description="Sets the size of a subgroup for ERI calculation, "// &
                          "each of which with a full set of work grids, arrays or orbitals "// &
                          "depending on the method of grids (work grids, arrays, orbitals). "// &
                          "Small numbers reduce communication but increase the memory demands. "// &
                          "A negative number indicates all processes (default).", &
                          usage="GROUP_SIZE 2", type_of_var=integer_t, &
                          default_i_val=-1)
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

END MODULE input_cp2k_dft
