!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the print section of the dft input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      07.2024 moved out of input_cp2k_dft [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_print_dft
   USE basis_set_types, ONLY: basis_sort_default, &
                              basis_sort_zet
   USE bibliography, ONLY: &
      Andermatt2016, Andreussi2012, Avezac2005, BaniHashemian2016, Becke1988b, Bengtsson1999, &
      Blochl1995, Brehm2018, Brelaz1979, Dewar1977, Dewar1985, Dudarev1997, Dudarev1998, &
      Ehrhardt1985, Eriksen2020, Fattebert2002, Golze2017a, Golze2017b, Guidon2010, &
      Heinzmann1976, Holmberg2017, Holmberg2018, Iannuzzi2005, Iannuzzi2006, Iannuzzi2007, &
      Knizia2013, Kolafa2004, Krack2000, Krack2002, Kuhne2007, Kunert2003, Lippert1997, &
      Lippert1999, Lu2004, Merlot2014, Perdew1981, Repasky2002, Rocha2006, Schenter2008, Schiffmann2015, &
      Shigeta2001, Stewart1982, Stewart1989, Stewart2007, Thiel1992, VanVoorhis2015, &
      VandeVondele2003, VandeVondele2005a, VandeVondele2005b, VandeVondele2006, Weber2008, &
      Yin2017, Pracht2019, Caldeweyher2019, Caldeweyher2020
   USE cp_output_handling, ONLY: add_last_numeric, &
                                 cp_print_key_section_create, &
                                 debug_print_level, &
                                 high_print_level, &
                                 low_print_level, &
                                 medium_print_level, &
                                 silent_print_level
   USE cp_spline_utils, ONLY: pw_interp, &
                              spline3_nopbc_interp, &
                              spline3_pbc_interp
   USE cp_units, ONLY: cp_unit_to_cp2k
   USE input_constants, ONLY: &
      atomic_guess, becke_cutoff_element, becke_cutoff_global, bqb_opt_exhaustive, &
      bqb_opt_normal, bqb_opt_off, bqb_opt_patient, bqb_opt_quick, broyden_type_1, &
      broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, broyden_type_2, &
      broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, casci_canonical, &
      cdft_alpha_constraint, cdft_beta_constraint, cdft_charge_constraint, &
      cdft_magnetization_constraint, cholesky_dbcsr, cholesky_inverse, cholesky_off, &
      cholesky_reduce, cholesky_restore, core_guess, diag_block_davidson, diag_block_krylov, &
      diag_filter_matrix, diag_ot, diag_standard, do_admm_aux_exch_func_bee, &
      do_admm_aux_exch_func_bee_libxc, do_admm_aux_exch_func_default, &
      do_admm_aux_exch_func_default_libxc, do_admm_aux_exch_func_none, &
      do_admm_aux_exch_func_opt, do_admm_aux_exch_func_opt_libxc, do_admm_aux_exch_func_pbex, &
      do_admm_aux_exch_func_pbex_libxc, do_admm_aux_exch_func_sx_libxc, &
      do_admm_basis_projection, do_admm_blocked_projection, do_admm_blocking_purify_full, &
      do_admm_charge_constrained_projection, do_admm_exch_scaling_merlot, &
      do_admm_exch_scaling_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace, &
      do_admm_purify_mcweeny, do_admm_purify_mo_diag, do_admm_purify_mo_no_diag, &
      do_admm_purify_none, do_admm_purify_none_dm, do_arnoldi, do_bch, do_cn, &
      do_ddapc_constraint, do_ddapc_restraint, do_em, do_etrs, do_full_density, do_gapw_gcs, &
      do_gapw_gct, do_gapw_log, do_iaoloc_energy, do_iaoloc_enone, do_iaoloc_l1, do_iaoloc_occ, &
      do_iaoloc_pm2, do_iaoloc_pm4, do_lri_inv, do_lri_inv_auto, do_lri_opt_all, &
      do_lri_opt_coeff, do_lri_opt_exps, do_lri_pseudoinv_diag, do_lri_pseudoinv_svd, &
      do_method_am1, do_method_dftb, do_method_gapw, do_method_gapw_xc, do_method_gpw, &
      do_method_lrigpw, do_method_mndo, do_method_mndod, do_method_ofgpw, do_method_pdg, &
      do_method_pm3, do_method_pm6, do_method_pm6fm, do_method_pnnl, do_method_rigpw, &
      do_method_rm1, do_method_xtb, do_pade, do_potential_coulomb, do_potential_id, &
      do_potential_short, do_potential_truncated, do_ppl_analytic, do_ppl_grid, &
      do_pwgrid_ns_fullspace, do_pwgrid_ns_halfspace, do_pwgrid_spherical, do_s2_constraint, &
      do_s2_restraint, do_se_is_kdso, do_se_is_kdso_d, do_se_is_slater, do_se_lr_ewald, &
      do_se_lr_ewald_gks, do_se_lr_ewald_r3, do_se_lr_none, do_spin_density, do_taylor, &
      ehrenfest, embed_diff, embed_fa, embed_grid_angstrom, embed_grid_bohr, embed_level_shift, &
      embed_none, embed_quasi_newton, embed_resp, embed_steep_desc, eri_method_full_gpw, &
      eri_method_gpw_ht, eri_operator_trunc, eri_operator_coulomb, eri_operator_erf, eri_operator_erfc, &
      eri_operator_gaussian, eri_operator_yukawa, gapw_1c_large, gapw_1c_medium, gapw_1c_orb, &
      gapw_1c_small, gapw_1c_very_large, gaussian, general_roks, gto_cartesian, gto_spherical, &
      high_spin_roks, history_guess, jacobian_fd1, jacobian_fd1_backward, &
      jacobian_fd1_central, jacobian_fd2, jacobian_fd2_backward, kg_color_dsatur, &
      kg_color_greedy, kg_tnadd_atomic, kg_tnadd_embed, kg_tnadd_embed_ri, kg_tnadd_none, &
      ls_2pnt, ls_3pnt, ls_gold, ls_none, manual_selection, mao_basis_ext, mao_basis_orb, &
      mao_basis_prim, mao_projection, mopac_guess, no_excitations, no_guess, no_solver, &
      numerical, oe_gllb, oe_lb, oe_none, oe_saop, oe_sic, orb_dx2, orb_dxy, orb_dy2, orb_dyz, &
      orb_dz2, orb_dzx, orb_px, orb_py, orb_pz, orb_s, ot_algo_irac, ot_algo_taylor_or_diag, &
      ot_chol_irac, ot_lwdn_irac, ot_mini_broyden, ot_mini_cg, ot_mini_diis, ot_mini_sd, &
      ot_poly_irac, ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
      ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
      ot_precond_solver_default, ot_precond_solver_direct, ot_precond_solver_inv_chol, &
      ot_precond_solver_update, outer_scf_basis_center_opt, outer_scf_becke_constraint, &
      outer_scf_cdft_constraint, outer_scf_ddapc_constraint, outer_scf_hirshfeld_constraint, &
      outer_scf_none, outer_scf_optimizer_bisect, outer_scf_optimizer_broyden, &
      outer_scf_optimizer_diis, outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls, &
      outer_scf_optimizer_none, outer_scf_optimizer_sd, outer_scf_optimizer_secant, &
      outer_scf_s2_constraint, plus_u_lowdin, plus_u_mulliken, plus_u_mulliken_charges, &
      qiskit_solver, radius_covalent, radius_default, radius_single, radius_user, radius_vdw, &
      random_guess, real_time_propagation, ref_charge_atomic, ref_charge_mulliken, rel_dkh, &
      rel_none, rel_pot_erfc, rel_pot_full, rel_sczora_mp, rel_trans_atom, rel_trans_full, &
      rel_trans_molecule, rel_zora, rel_zora_full, rel_zora_mp, restart_guess, &
      sccs_andreussi, sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7, &
      sccs_derivative_fft, sccs_fattebert_gygi, shape_function_density, shape_function_gaussian, &
      sic_ad, sic_eo, sic_list_all, sic_list_unpaired, sic_mauri_spz, sic_mauri_us, sic_none, &
      slater, smear_energy_window, smear_fermi_dirac, smear_list, sparse_guess, tddfpt_davidson, &
      tddfpt_excitations, tddfpt_lanczos, tddfpt_singlet, tddfpt_spin_cons, tddfpt_spin_flip, &
      tddfpt_triplet, use_mom_ref_coac, use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero, &
      use_restart_wfn, use_rt_restart, use_scf_wfn, wannier_projection, weight_type_mass, &
      weight_type_unit, wfi_aspc_nr, wfi_frozen_method_nr, wfi_linear_p_method_nr, &
      wfi_linear_ps_method_nr, wfi_linear_wf_method_nr, wfi_ps_method_nr, &
      wfi_use_guess_method_nr, wfi_use_prev_p_method_nr, wfi_use_prev_rho_r_method_nr, &
      wfi_use_prev_wf_method_nr, wfn_mix_orig_external, wfn_mix_orig_occ, wfn_mix_orig_virtual, &
      xas_1s_type, xas_2p_type, xas_2s_type, xas_3d_type, xas_3p_type, xas_3s_type, xas_4d_type, &
      xas_4f_type, xas_4p_type, xas_4s_type, xas_dip_len, xas_dip_vel, xas_dscf, xas_none, &
      xas_not_excited, xas_tdp_by_index, xas_tdp_by_kind, xas_tp_fh, xas_tp_flex, xas_tp_hh, &
      xas_tp_xfh, xas_tp_xhh, xes_tp_val, &
      no_admm_type, admm1_type, admm2_type, admms_type, admmp_type, admmq_type, &
      e_dens_total_hard_approx, e_dens_total_density, e_dens_soft_density
   USE input_cp2k_almo, ONLY: create_almo_scf_section
   USE input_cp2k_distribution, ONLY: create_distribution_section
   USE input_cp2k_ec, ONLY: create_ec_section
   USE input_cp2k_exstate, ONLY: create_exstate_section
   USE input_cp2k_external, ONLY: create_ext_den_section, &
                                  create_ext_pot_section, &
                                  create_ext_vxc_section
   USE input_cp2k_field, ONLY: create_efield_section, &
                               create_per_efield_section
   USE input_cp2k_kpoints, ONLY: create_kpoint_set_section
   USE input_cp2k_loc, ONLY: create_localize_section, &
                             print_wanniers
   USE input_cp2k_ls, ONLY: create_ls_scf_section
   USE input_cp2k_mm, ONLY: create_dipoles_section, &
                            create_neighbor_lists_section
   USE input_cp2k_poisson, ONLY: create_poisson_section
   USE input_cp2k_projection_rtp, ONLY: create_projection_rtp_section
   USE input_cp2k_rsgrid, ONLY: create_rsgrid_section
   USE input_cp2k_tb, ONLY: create_dftb_control_section, &
                            create_xtb_control_section
   USE input_cp2k_transport, ONLY: create_transport_section
   USE input_cp2k_voronoi, ONLY: create_print_voronoi_section
   USE input_cp2k_scf, ONLY: create_scf_section, &
                             create_cdft_control_section
   USE input_cp2k_xc, ONLY: create_xc_fun_section, &
                            create_xc_section
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type
   USE input_val_types, ONLY: char_t, &
                              integer_t, &
                              lchar_t, &
                              logical_t, &
                              real_t
   USE kinds, ONLY: dp
   USE pw_grids, ONLY: do_pw_grid_blocked_false, &
                       do_pw_grid_blocked_free, &
                       do_pw_grid_blocked_true
   USE pw_spline_utils, ONLY: no_precond, &
                              precond_spl3_1, &
                              precond_spl3_2, &
                              precond_spl3_3, &
                              precond_spl3_aint, &
                              precond_spl3_aint2
   USE qs_density_mixing_types, ONLY: create_mixing_section
   USE qs_fb_input, ONLY: create_filtermatrix_section
   USE qs_mom_types, ONLY: create_mom_section
   USE string_utilities, ONLY: newline, &
                               s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_print_dft'

   PUBLIC :: create_print_dft_section, create_pdos_section

CONTAINS

! **************************************************************************************************
!> \brief Create the print dft section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_print_dft_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          "GTH-PPNL neighbor lists (non-local part of the "// &
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
                          name="sab_xtbe", &
                          description="Activates the printing of the xTB sr-Coulomb "// &
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
                          name="sab_xtb_nonbond", &
                          description="Activates the printing of nonbonded interaction from (xTB) "// &
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

      CALL keyword_create(keyword, __LOCATION__, &
                          name="sab_cneo", &
                          description="Activates the printing of the nuclear orbital "// &
                          "nuclear repulsion neighbor lists (erfc potential)", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "SUBCELL", &
                                       description="Activates the printing of the subcells used for the "// &
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
         description="If the printkey is activated compute and print the matrix of the exchange and correlation potential. "// &
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
         description="Controls the printing of the molecular orbital (MO) information. The requested MO information "// &
         "is printed for all occupied MOs by default. Use the MO_INDEX_RANGE keyword to restrict the number "// &
         "of the MOs or to print the MO information for unoccupied MOs. With diagonalization, additional MOs "// &
         "have to be made available for printout using the ADDED_MOS keyword in the SCF section. With OT, "// &
         "it is sufficient to specify the desired MO_INDEX_RANGE. The OT eigensolver can be controlled with "// &
         "the EPS_LUMO and MAX_ITER_LUMO keywords in the SCF section.", &
         print_level=high_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="CARTESIAN", &
                          description="Print the MOs in the Cartesian basis instead of the default spherical basis.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ENERGIES", &
                          variants=s2a("EIGENVALUES", "EIGVALS"), &
                          description="Print the MO energies (eigenvalues).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="COEFFICIENTS", &
                          variants=s2a("EIGENVECTORS", "EIGVECS"), &
                          description="Print the MO coefficients (eigenvectors).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OCCUPATION_NUMBERS", &
                          variants=s2a("OCCNUMS"), &
                          description="Print the MO occupation numbers.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OCCUPATION_NUMBERS_STATS", &
                          variants=s2a("OCCNUMSTATS"), &
                          description="Print some stats (max number of occupied MOs, etc.) of the MO occupation numbers."// &
                          " First logical toggles stats printing, first real is the occupied threshold.", &
                          type_of_var=char_t, n_var=-1, &
                          default_c_vals=[".FALSE.", "1.0E-6 "], &
                          lone_keyword_c_vals=[".TRUE."], &
                          usage="OCCUPATION_NUMBERS_STATS {Logical} [{Real}]")
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specify the number of digits used to print the MO information.", &
                          default_i_val=6)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="MO_INDEX_RANGE", &
                          variants=s2a("MO_RANGE", "RANGE"), &
                          description="Print only the requested subset of MOs. The indices of the first and "// &
                          "the last MO have to be specified to define the range. -1 as the last MO index "// &
                          "prints all available orbitals with diagonalisation (ADDED_MOS) and all orbitals with OT.", &
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

      CALL section_create(subsection, __LOCATION__, name="TREXIO", &
                          description="Write a TREXIO file to disk.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="FILENAME", &
                          description="Body of Filename for the trexio file.", &
                          usage="FILENAME {name}", default_c_val="TREXIO", &
                          type_of_var=char_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CARTESIAN", &
                          description="Store the MOs in the Cartesian basis instead of the default spherical basis.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
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
                                       description="Controls the checking and printing of an estimate "// &
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
                                       "the electronic density and, for LSD calculations, the spin density.", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_INCLUDE", &
                          description="Which parts of the density to include. In GAPW the electronic density "// &
                          "is divided into a hard and a soft component, and the default (TOTAL_HARD_APPROX) "// &
                          "is to approximate the hard density as a spherical gaussian and to print the smooth "// &
                          "density accurately. This avoids potential artefacts originating from the hard density. "// &
                          "If the TOTAL_DENSITY keyword is used the hard density will be computed more accurately "// &
                          "but may introduce non-physical features. The SOFT_DENSITY keyword will lead to only the "// &
                          "soft density being printed. In GPW these options have no effect and the cube file will "// &
                          "only contain the valence electron density.", &
                          usage="DENSITY_INCLUDE TOTAL_HARD_APPROX", &
                          enum_c_vals=s2a("TOTAL_HARD_APPROX", "TOTAL_DENSITY", "SOFT_DENSITY"), &
                          enum_desc=s2a("Print (hard+soft) density where the hard components shape is approximated", &
                                        "Print (hard+soft) density. Only has an effect "// &
                                        "if PAW atoms are present. NOTE: The total "// &
                                        "in real space might exhibit unphysical features "// &
                                        "like spikes due to the finite and thus "// &
                                        "truncated g vector", &
                                        "Print only the soft density"), &
                          enum_i_vals=(/e_dens_total_hard_approx, &
                                        e_dens_total_density, &
                                        e_dens_soft_density/), &
                          default_i_val=e_dens_total_hard_approx)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XRD_INTERFACE", &
                          description="It activates the print out of exponents and coefficients for the"// &
                          " Gaussian expansion of the core densities, based on atom calculations for each kind."// &
                          " The resulting core dansities are needed to compute the form factors."// &
                          " If GAPW the local densities are also given in terms of a Gaussian expansion,"// &
                          " by fitting the difference between local-fhard and local-soft density for each atom."// &
                          " In this case the keyword SOFT_DENSITY is enabled.", &
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
                                       "electron density dominates and negative where the core is. "// &
                                       "When GPW is enabled this will simply print the combined density "// &
                                       "of the valence electrons and charge-balanced core. In GAPW the "// &
                                       "electronic density (hard+soft plus a correction term) is printed "// &
                                       "together with the charge-balanced core density to produce a complete "// &
                                       "representation of the total density.", &
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
                                       description="Controls the printing of a cube file with eletrostatic"// &
                                       " potential generated by the total density (electrons+ions). It is"// &
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
                                       description="Controls the printing of a cube file with external"// &
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

      ! Output of BQB volumetric files
      CALL cp_print_key_section_create(print_key, __LOCATION__, name="E_DENSITY_BQB", &
                                       description="Controls the output of the electron density to the losslessly"// &
                                       " compressed BQB file format, see [Brehm2018]"// &
                                       " (via LibBQB see <https://brehm-research.de/bqb>)."// &
                                       " Currently does not work with changing cell vector (NpT ensemble).", &
                                       print_level=debug_print_level + 1, filename="", &
                                       citations=(/Brehm2018/))

      CALL keyword_create(keyword, __LOCATION__, name="SKIP_FIRST", &
                          description="Skips the first step of a MD run (avoids duplicate step if restarted).", &
                          usage="SKIP_FIRST T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STORE_STEP_NUMBER", &
                          description="Stores the step number and simulation time in the comment line of each BQB"// &
                          " frame. Switch it off for binary compatibility with original CP2k CUBE files.", &
                          usage="STORE_STEP_NUMBER F", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK", &
                          description="Performs an on-the-fly decompression of each compressed BQB frame to check"// &
                          " whether the volumetric data exactly matches, and aborts the run if not so.", &
                          usage="CHECK T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OVERWRITE", &
                          description="Specify this keyword to overwrite the output BQB file if"// &
                          " it already exists. By default, the data is appended to an existing file.", &
                          usage="OVERWRITE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HISTORY", &
                          description="Controls how many previous steps are taken into account for extrapolation in"// &
                          " compression. Use a value of 1 to compress the frames independently.", &
                          usage="HISTORY 10", n_var=1, default_i_val=10, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAMETER_KEY", &
                          description="Allows to supply previously optimized compression parameters via a"// &
                          " parameter key (alphanumeric character sequence starting with 'at')."// &
                          " Just leave away the 'at' sign here, because CP2k will otherwise"// &
                          " assume it is a variable name in the input", &
                          usage="PARAMETER_KEY <KEY>", n_var=1, default_c_val="", type_of_var=char_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPTIMIZE", &
                          description="Controls the time spent to optimize the parameters for compression efficiency.", &
                          usage="OPTIMIZE {OFF,QUICK,NORMAL,PATIENT,EXHAUSTIVE}", repeats=.FALSE., n_var=1, &
                          default_i_val=bqb_opt_quick, &
                          enum_c_vals=s2a("OFF", "QUICK", "NORMAL", "PATIENT", "EXHAUSTIVE"), &
                          enum_desc=s2a("No optimization (use defaults)", "Quick optimization", &
                                        "Standard optimization", "Precise optimization", "Exhaustive optimization"), &
                          enum_i_vals=(/bqb_opt_off, bqb_opt_quick, bqb_opt_normal, bqb_opt_patient, bqb_opt_exhaustive/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Voronoi Integration via LibVori
      CALL create_print_voronoi_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! cube files for data generated by the implicit (generalized) Poisson solver
      CALL create_implicit_psolver_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! ZMP adding the print section for  the v_xc cube
      CALL cp_print_key_section_create(print_key, __LOCATION__, "v_xc_cube", &
                                       description="Controls the printing of a cube file with xc"// &
                                       " potential generated by the ZMP method (for the moment). It is"// &
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
                                       description="Controls the printing of cube files with electric"// &
                                       " field generated by the total density (electrons+ions). It is"// &
                                       " valid only for QS with GPW formalism.", &
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
                                       description="Controls printing of cube files with"// &
                                       " the electron localization function (ELF). Note that"// &
                                       " the value of ELF is defined between 0 and 1: Pauli kinetic energy density normalized"// &
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
                                       description="Controls the printing of cube files with the local"// &
                                       " energy. It is valid only for QS with GPW/GAPW formalism."// &
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
                                       description="Controls the printing of cube files with the local"// &
                                       " stress. It is valid only for QS with GPW/GAPW formalism."// &
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
                          description="Calculate expectation values of the el. multipole moments in their velocity "// &
                          "representation during RTP. Implemented up to el. quadrupole moment.", &
                          usage="VEL_REPRS yes", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="COM_NL", &
                          description="Include non local commutator for velocity representations. "// &
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

      ! Print EEQ Charges
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EEQ_CHARGES", &
                                       description="Controls the printing of the EEQ charges", &
                                       print_level=debug_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1, &
                                       citations=(/Pracht2019, Caldeweyher2019, Caldeweyher2020/))
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
                          description="Indexes of the atoms minimal basis to be printed as cube files "// &
                          "This keyword can be repeated several times "// &
                          "(useful if you have to specify many indexes).", &
                          usage="ATOM_LIST 1 2", &
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

      !Energy Decomposition Analysis
      CALL cp_print_key_section_create(print_key, __LOCATION__, "ENERGY_DECOMPOSITION_ANALYSIS", &
                                       description="Controls energy decomposition analysis", &
                                       print_level=debug_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1, &
                                       citations=(/Eriksen2020/))
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_ORB_CANONICAL", &
                          description="Use reference orbitals in canonical form.", &
                          usage="REFERENCE_ORB_CANONICAL {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SKIP_LOCALIZATION", &
                          description="Don't localize the MOs.", &
                          usage="SKIP_LOCALIZATION {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DETAILED_ENERGY", &
                          description="Calculate detailed atomic decomposition energies.", &
                          usage="DETAILED_ENERGY {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EWALD_ALPHA_PARAMETER", &
                          description="Calculate Energy Decomposition for a specific alpha value. "// &
                          "alpha = 1/(2*rc**2), see GTH pseudopotentials.", &
                          usage="EWALD_ALPHA_PARAMETER alpha", repeats=.FALSE., n_var=1, &
                          default_r_val=0.0_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! IAO (Intrinsic atomic orbital) analysis
      CALL cp_print_key_section_create(print_key, __LOCATION__, "IAO_ANALYSIS", &
                                       description="Controls the printing of the IAO (intrinsic atomic orbital) analysis", &
                                       print_level=debug_print_level, filename="__STD_OUT__", &
                                       common_iter_levels=1, &
                                       citations=(/Knizia2013/))
      CALL keyword_create(keyword, __LOCATION__, name="EPS_SVD", &
                          description="Threshold for matrix inversion eigenvalues.", &
                          usage="EPS_SVD reps", repeats=.FALSE., n_var=1, &
                          default_r_val=0.0_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_OCC", &
                          description="Threshold in occupation for vectors included.", &
                          usage="EPS_OCC reps", repeats=.FALSE., n_var=1, &
                          default_r_val=0.0_dp, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_CHARGES", &
                          description="Calculate atomic charges from IAO.", &
                          usage="ATOMIC_CHARGES {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      ! IAO_MOLDEN
      NULLIFY (sub_print_key)
      CALL cp_print_key_section_create(sub_print_key, __LOCATION__, "IAO_MOLDEN", &
                                       description="Write the IAO basis in Molden file format, for visualisation.", &
                                       print_level=debug_print_level + 1, add_last=add_last_numeric, filename="IAOBAS")
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
      ! IAO_CUBES
      NULLIFY (sub_print_key)
      CALL cp_print_key_section_create(sub_print_key, __LOCATION__, "IAO_CUBES", &
                                       description="Controls the printing of the IAO basis "// &
                                       "as *.cube files.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Indices of the atoms to be included in basis CUBE file printing. ", &
                          usage="ATOM_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)
      ! One Center Expansion of IAO
      NULLIFY (sub_print_key)
      CALL cp_print_key_section_create(sub_print_key, __LOCATION__, "ONE_CENTER_EXPANSION", &
                                       description="Calculates single center expansion of IAOs ", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="LMAX", &
                          description="Maximum l quantum number used in the expansion.", &
                          usage="LMAX 2", n_var=1, default_i_val=3, type_of_var=integer_t)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NBAS", &
                          description="Max number of basis functions used in the expansion."// &
                          " Default is determined by the orbital basis set.", &
                          usage="NBAS 10", n_var=1, default_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the OCE basis files when it already exists", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)
      ! Intrinsic Bond orbitals
      NULLIFY (sub_print_key)
      CALL cp_print_key_section_create(sub_print_key, __LOCATION__, "BOND_ORBITALS", &
                                       description="Calculate intrinsic bond orbitals using "// &
                                       "localized MOs in IAO basis.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="LOCALIZATION_OPERATOR", &
                          description="Operator to be optimized for orbital localization", &
                          enum_c_vals=s2a("PIPEK_MEZEY", "PIPEK_MEZEY_4", "L1NORM"), &
                          enum_i_vals=(/do_iaoloc_pm2, do_iaoloc_pm4, do_iaoloc_l1/), &
                          enum_desc=s2a("Use Pipek-Mezey operator (order 2)", &
                                        "Use Pipek-Mezey operator (order 4)", &
                                        "Use L1 norm"), &
                          default_i_val=do_iaoloc_pm2)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_LOCALIZATION_FUNCTION", &
                          description="Function for energy localization: f(e_i), e_i orbital energy", &
                          enum_c_vals=s2a("NONE", "ENERGY", "OCCUPATION"), &
                          enum_i_vals=(/do_iaoloc_enone, do_iaoloc_energy, do_iaoloc_occ/), &
                          enum_desc=s2a("Don't use energy localization.", &
                                        "Use orbital energies for localization.", &
                                        "Use occupation numbers for localization."), &
                          default_i_val=do_iaoloc_enone)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_LOCALIZATION_WEIGHT", &
                          description="Weight given to energy localization, using f(e_i) function", &
                          usage="ENERGY_LOCALIZATION_WEIGHT 0.1", n_var=1, &
                          default_r_val=0.0_dp, type_of_var=real_t)
      CALL section_add_keyword(sub_print_key, keyword)
      CALL keyword_release(keyword)

      ! CHARGE CENTER AND SPREAD
      NULLIFY (subsection)
      CALL cp_print_key_section_create(subsection, __LOCATION__, "CHARGE_CENTER", &
                                       description="Calculation and printing of centers and spreads "// &
                                       "of localized orbitals.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="POSITION_OPERATOR_BERRY", &
                          description="Use Berry phase position operator.", &
                          usage="POSITION_OPERATOR_BERRY T", n_var=1, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(sub_print_key, subsection)
      CALL section_release(subsection)
      ! IBO_MOLDEN
      NULLIFY (subsection)
      CALL cp_print_key_section_create(subsection, __LOCATION__, "IBO_MOLDEN", &
                                       description="Write the IBO orbitals in Molden file format, for visualisation.", &
                                       print_level=debug_print_level + 1, add_last=add_last_numeric, filename="IBOBAS")
      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specifies the number of significant digits retained. 3 is OK for visualization.", &
                          usage="NDIGITS {int}", &
                          default_i_val=3)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GTO_KIND", &
                          description="Representation of Gaussian-type orbitals", &
                          default_i_val=gto_spherical, &
                          enum_c_vals=s2a("CARTESIAN", "SPHERICAL"), &
                          enum_desc=s2a( &
                          "Cartesian Gaussian orbitals. Use with caution", &
                          "Spherical Gaussian orbitals. Incompatible with VMD"), &
                          enum_i_vals=(/gto_cartesian, gto_spherical/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(sub_print_key, subsection)
      CALL section_release(subsection)
      ! IAO_CUBES
      NULLIFY (subsection)
      CALL cp_print_key_section_create(subsection, __LOCATION__, "IBO_CUBES", &
                                       description="Controls the printing of the IBO orbitals "// &
                                       "as *.cube files.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="STATE_LIST", &
                          description="Indices of the orbitals to be included in IBO CUBE file printing. ", &
                          usage="STATE_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(sub_print_key, subsection)
      CALL section_release(subsection)
      CALL section_add_subsection(print_key, sub_print_key)
      CALL section_release(sub_print_key)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)
      ! END OF IAO_ANALYSIS SECTION

      !DOS from density matrix
      CALL cp_print_key_section_create(print_key, __LOCATION__, "ENERGY_WINDOWS", &
                                       description="Controls the printing of the DOS from the density matrix. "// &
                                       "This allows the calculation of the DOS even in density matrix based "// &
                                       "REAL_TIME_PROPAGATION and LS_SCF. "// &
                                       "However, it requires a cubically scaling diagonalization of the Hamiltonian. "// &
                                       "Hartree-Fock NYI, values will be wrong. "// &
                                       "Careful, the orbitals in rtp/emd are not actually eigenstates of the Hamiltonian. "// &
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
                          usage="PRINT_CUBES .TRUE.", &
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
      CALL keyword_create(keyword, __LOCATION__, name="Real_space", &
                          description="Print the KS matrix in real-space instead of k-space.. ", &
                          usage="REAL_SPACE {logical}", &
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
      CALL keyword_create(keyword, __LOCATION__, name="Real_space", &
                          description="Print the overlap matrix in real-space instead of k-space.. ", &
                          usage="REAL_SPACE {logical}", &
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
         print_key, __LOCATION__, name="CHARGEMOL", &
         description="Write .wfx input file for Chargemol", &
         print_level=debug_print_level + 1, &
         filename="CHARGEMOL", &
         add_last=add_last_numeric)
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC", &
                          description="Write information about cell periodicity.", &
                          usage="PERIODIC {LOGICAL}", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
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
         "total charge density $\rho^\text{tot}$ used by the SCCS model.", &
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
         "polarisation charge density $\rho^\text{pol}$ used by the SCCS model with the "// &
         "total charge density $\rho^\text{tot} = \rho^\text{sol} + \rho^\text{pol}", &
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
         "polarisation potential $\phi^\text{pol}$ used by the SCCS model with the "// &
         "total potential $\phi^\text{tot} = \phi^\text{sol} + \phi^\text{pol}$", &
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
      CALL create_kpoint_set_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_bandstructure_section

! **************************************************************************************************
!> \brief creates the input section for dealing with homo lumos, including dumping cubes
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_mo_cubes_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_CUBES", &
                                       description="Controls the printing of the molecular orbitals (MOs) as cube files."// &
                                       " It can be used during a Real Time calculation to print the MOs."// &
                                       " In this case, the density corresponding to the time dependent MO is printed"// &
                                       " instead of the wave-function.", &
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

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "DOS", &
                                       description="Print Density of States (DOS) (only available states from SCF)", &
                                       print_level=debug_print_level, common_iter_levels=1, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the DOS obtained at different iterations to the output file. "// &
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

      CALL keyword_create(keyword, __LOCATION__, name="NDIGITS", &
                          description="Specify the number of digits used to print density and occupation", &
                          default_i_val=4)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MP_GRID", &
                          description="Specify a Monkhorst-Pack grid with which to compute the density of states. "// &
                          "Works only for a k-point calculation", &
                          usage="MP_GRID {integer} {integer} {integer}", default_i_vals=(/-1/), &
                          n_var=3, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_dos_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_pdos_section(print_key)

      TYPE(section_type), POINTER                        :: print_key

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (subsection)
      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PDOS", &
                                       description="Print out the DOS projected per kind and angular momentum  ", &
                                       print_level=debug_print_level, common_iter_levels=1, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="COMPONENTS", &
                          description="Print out PDOS distinguishing all angular momentum components.", &
                          usage="COMPONENTS", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the PDOS obtained at different iterations to the PDOS output file. "// &
                          "By default the file is overwritten", &
                          usage="APPEND", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="NLUMO", &
         description="Number of virtual orbitals to be added to the MO set (-1=all)."//newline// &
         "CAUTION: Setting this value to be higher than the "// &
         "number of states present may cause a Cholesky error."//newline// &
         "This keyword only applies to MO sets, not to XAS_TDP states.", &
         usage="NLUMO integer", default_i_val=0)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OUT_EACH_STATE", &
                          variants=(/"OUT_EACH_MO"/), &
                          description="Output on the status of the calculation every OUT_EACH_MO states. If -1 no output", &
                          usage="OUT_EACH_STATE integer", default_i_val=-1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="LDOS", &
                          description="Controls the printing of local PDOS, projected on subsets"// &
                          " of atoms given through lists", &
                          n_keywords=4, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="COMPONENTS", &
                          description="Print out PDOS distinguishing all angular momentum components.", &
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
                          "this range from at least one atom of the list are considered", &
                          usage="XRANGE -10.0 10.0", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="YRANGE", &
                          description="range of positions in Cartesian direction y: all grid points within "// &
                          "this range from at least one atom of the list are considered", &
                          usage="YRANGE -10.0 10.0", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ZRANGE", &
                          description="range of positions in Cartesian direction z: all grid points within "// &
                          "this range from at least one atom of the list are considered", &
                          usage="ZRANGE -10.0 10.0", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERANGE", &
                          description="Only project states with the energy values in the given interval. "// &
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
                          description="Bias energy for scanning tunneling microscopy (STM) image generation. "// &
                          "Orbital densities are summed according to the bias energy. "// &
                          "For negative values, states in the range ]EF+bias,EF] are summed, "// &
                          "While positive values sum states in the range [EF,EF+bias[. "// &
                          "If positive biases are used, sufficiently many unoccupied stated "// &
                          "(see ADDED_MOS and NLUMO ) should be computed.", &
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

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (subsection)
      NULLIFY (keyword)

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="WFN_MIX", &
                          description="A section that allows manipulation of the MO coeffs, "// &
                          "e.g. for changing a ground state into an excited state. "// &
                          "Starting from a copy of the original MOs, changes can be made "// &
                          "by adding linear combinations of HOMO/LUMO of the original MOs to the result MOs. "// &
                          "This method is called after an SCF optimization or before an RTP run if "// &
                          "INITIAL_WFN=RESTART_WFN. Note that if called after an SCF optimization, a restart file "// &
                          "with the mixed MOs is saved. This is not the case for an RTP with "// &
                          "INITIAL_WFN=RESTART_WFN.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="OVERWRITE_MOS", &
                        description="If set to True, the active molecular orbitals in memory will be replaced by the mixed wfn "// &
                        "at the end of the wfn mixing procedure. For instance, you can then use this new set of MOs to perform "// &
                     "RTP or EMD directly. Note that in the case of an RTP run with INITIAL_WFN=RESTART_WFN, the OVERWRITE_MOS "// &
                          "keyword is not used.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="UPDATE", &
                          description="Update a result MO with with a linear combination of original MOs."// &
                          " This section can be repeated to build arbitrary linear combinations using repeatedly y=a*y+b*x. "// &
                          "RESULT is (y), RESULT_SCALE is (a), ORIG is (x), ORIG_SCALE is (b)", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_MO_INDEX", &
                          description="Index of the MO (y) to be modified. Counting down in energy: "// &
                          "set to 1 for the highest MO, to 3 for the highest MO-2.", &
                          usage="RESULT_MO_INDEX 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_MARKED_STATE", &
                          description="Specifies the MO according to "// &
                          "the marks set in MOLECULAR_STATES. The value corresponds to the repetition "// &
                          "of MARK_STATES in MOLECULAR_STATES", &
                          usage="RESULT_MARKED_STATE 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REVERSE_MO_INDEX", &
                          description="Reverses the index order of the OCCUPIED and EXTERNAL MOs. With this keyword "// &
                          "ORIG_MO_INDEX/RESULT_MO_INDEX 1 point to the lowest energy MO (instead of the highest) "// &
                          "and counts up in energy. The VIRTUAL MOs indexing is unchanged.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
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
                          "If ORIG_TYPE is OCCUPIED, it counts down in energy: set to 1 to point to "// &
                          "the highest MO and to 3 for the highest MO-2. "// &
                          "If ORIG_TYPE is VIRTUAL, it counts up in energy: set to 1 to point to "// &
                          "the lowest virtual MO and to 3 for the lowest MO+2. "// &
                          "If ORIG_TYPE is EXTERNAL, it counts down in energy for the external "// &
                          "set of MOs: set to 1 to point to the highest MO and to 3 for the highest MO-2. "// &
                          "Do not set to zero or negative values.", &
                          usage="ORIG_MO_INDEX 1", type_of_var=integer_t, default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_MARKED_STATE", &
                          description="Specifies the MO according to "// &
                          "the marks set in MOLECULAR_STATES. The value corresponds to the repetition "// &
                          "of MARK_STATES in MOLECULAR_STATES", &
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
                          description="Scaling factor of the original variable (b).", &
                          usage="ORIG_SCALE 0.0", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_TYPE", &
                          description="Type of the original MO. Note that if ADDED_MOS was used in the "// &
                          "SCF construction of the MO matrix, these extra MOs are also treated as OCCUPIED. ", &
                          enum_c_vals=s2a("OCCUPIED", "VIRTUAL", 'EXTERNAL'), &
                          usage="ORIG_TYPE OCCUPIED", &
                          default_i_val=wfn_mix_orig_occ, &
                          enum_desc=s2a("The original MO is the result of the SCF procedure. This can also contain "// &
                                        "unoccupied MOs if the SCF%ADDED_MOS keyword was used.", &
                                        "The original MO is taken from the result of additional MOs calculated a "// &
                                        "posteriori of the SCF by request of the user. E.g. by specifying print%mo_cubes%nlumo. ", &
                                        "The orginal MO is from an external .wfn file. Use the keyword "// &
                                        "ORIG_EXT_FILE_NAME to define its name."), &
                          enum_i_vals=(/wfn_mix_orig_occ, wfn_mix_orig_virtual, wfn_mix_orig_external/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORIG_EXT_FILE_NAME", &
                          description="Name of the wavefunction file to read the original MO from. "// &
                          "For instance, a restart wfn file from SCF calculation or an excited state from XAS_TDP calculation. "// &
                          "If no file is specified, the run will crash. "// &
                          "Currently, a RTP restart file (.rtpwfn) cannot be used as reference. "// &
                          "Currently, this file SHALL have the basis set, number of MO and the same number of spin as the one "// &
                          "from the SCF cycle.", &
                          usage="ORIG_EXT_FILE_NAME <FILENAME>", &
                          default_lc_val="EMPTY")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_wfn_mix_section

! **************************************************************************************************
!> \brief Creates the section for cube files related to the implicit Poisson solver.
!> \param section the section to be created
!> \par History
!>      03.2016 refactored from create_print_dft_section [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE create_implicit_psolver_section(section)
      TYPE(section_type), POINTER                        :: section

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
         description="Controls the printing of cube files with unit step functions (constraints)"// &
         " representing Dirichlet-type (boundary) regions defined in the implicit (generalized) Poisson"// &
         " solver section. The regions remain unchanged throughout the calculations. If the Dirichlet"// &
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
!> \brief creates the interpolation section for the periodic QM/MM
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_gspace_interp_section(section)
      TYPE(section_type), POINTER                        :: section

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
                          usage="aint_precond copy", &
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
                          usage="precond copy", &
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

END MODULE input_cp2k_print_dft
