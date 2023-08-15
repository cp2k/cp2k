#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(fypp-sources)
include(cp2k_utils)
include(GNUInstallDirs) # required to get a proper LIBDIR variable
include(CMakePackageConfigHelpers)
configure_file(base/base_uses.f90 base/base_uses.f90 @ONLY)

# =================================================================================================
# SOURCE FILE LISTS
# =================================================================================================

# do not worry about formatting or sorting. make pretty will fix this for you

# the sources files are split by directory which hopefully make the list a bit
# clearer. Any Fortran new file can be added to CP2K_SRCS_F list. the list
# sorted but it is *not* a requirement. C files should be added preferencially
# to CP2K_SRCS_C etc....
#
# The other option is to simply create a new list with the syntax list(APPEND
# MY_FOO foo.c foo.F foo.cu) and add this list to the list of files to compile.

list(
  APPEND
  CP2K_SRCS_F
  admm_dm_methods.F
  admm_dm_types.F
  admm_methods.F
  admm_types.F
  admm_utils.F
  almo_scf_diis_types.F
  almo_scf_env_methods.F
  almo_scf.F
  almo_scf_lbfgs_types.F
  almo_scf_methods.F
  almo_scf_optimizer.F
  almo_scf_qs.F
  almo_scf_types.F
  al_system_types.F
  atom_admm_methods.F
  atom_basis.F
  atom_electronic_structure.F
  atom_energy.F
  atom.F
  atom_fit.F
  atom_grb.F
  atomic_charges.F
  atom_kind_orbitals.F
  atom_operators.F
  atom_optimization.F
  atom_output.F
  atom_pseudo.F
  atom_sgp.F
  atoms_input.F
  atom_types.F
  atom_upf.F
  atom_utils.F
  atom_xc.F
  auto_basis.F
  bandstructure_methods.F
  basis_set_output.F
  beta_gamma_psi.F
  block_p_types.F
  bse.F
  bsse.F
  cell_methods.F
  colvar_methods.F
  colvar_utils.F
  commutator_rkinetic.F
  commutator_rpnl.F
  constraint_3x3.F
  constraint_4x6.F
  constraint_clv.F
  constraint.F
  constraint_fxd.F
  constraint_util.F
  constraint_vsite.F
  core_ae.F
  core_ppl.F
  core_ppnl.F
  cp2k_debug.F
  cp2k_info.F
  cp_control_types.F
  cp_control_utils.F
  cp_dbcsr_cholesky.F
  cp_dbcsr_cp2k_link.F
  cp_dbcsr_diag.F
  cp_dbcsr_operations.F
  cp_dbcsr_output.F
  cp_ddapc.F
  cp_ddapc_forces.F
  cp_ddapc_methods.F
  cp_ddapc_types.F
  cp_ddapc_util.F
  cp_eri_mme_interface.F
  cp_external_control.F
  parallel_gemm_api.F
  cp_realspace_grid_cube.F
  cp_realspace_grid_init.F
  cp_subsys_methods.F
  cp_symmetry.F
  cryssym.F
  csvr_system_types.F
  csvr_system_utils.F
  ct_methods.F
  ct_types.F
  debug_os_integrals.F
  dft_plus_u.F
  dgemm_counter_types.F
  distribution_2d_types.F
  distribution_methods.F
  dkh_main.F
  dm_ls_chebyshev.F
  dm_ls_scf_create.F
  dm_ls_scf_curvy.F
  dm_ls_scf.F
  dm_ls_scf_methods.F
  dm_ls_scf_qs.F
  dm_ls_scf_types.F
  domain_submatrix_methods.F
  domain_submatrix_types.F
  ec_efield_local.F
  ec_environment.F
  ec_env_types.F
  ec_methods.F
  ec_orth_solver.F
  ed_analysis.F
  efield_tb_methods.F
  efield_utils.F
  eip_environment.F
  eip_environment_types.F
  eip_silicon.F
  embed_environment.F
  embed_main.F
  embed_types.F
  energy_corrections.F
  environment.F
  eri_mme/eri_mme_error_control.F
  eri_mme/eri_mme_gaussian.F
  eri_mme/eri_mme_integrate.F
  eri_mme/eri_mme_lattice_summation.F
  eri_mme/eri_mme_test.F
  eri_mme/eri_mme_types.F
  eri_mme/eri_mme_util.F
  et_coupling.F
  et_coupling_proj.F
  et_coupling_types.F
  ewald_environment_types.F
  ewald_methods_tb.F
  ewald_pw_methods.F
  ewald_pw_types.F
  ewalds.F
  ewalds_multipole.F
  ewald_spline_util.F
  excited_states.F
  exclusion_types.F
  ex_property_calculation.F
  exstates_types.F
  extended_system_types.F
  external_potential_methods.F
  f77_interface.F
  farming_methods.F
  farming_types.F
  fermi_utils.F
  fist_efield_methods.F
  fist_efield_types.F
  fist_energy_types.F
  fist_environment.F
  fist_environment_types.F
  fist_force.F
  fist_intra_force.F
  fist_main.F
  fist_neighbor_list_control.F
  fist_neighbor_lists.F
  fist_neighbor_list_types.F
  fist_nonbond_env_types.F
  fist_nonbond_force.F
  fist_pol_scf.F
  force_env_methods.F
  force_env_types.F
  force_env_utils.F
  force_fields_all.F
  force_fields_ext.F
  force_fields.F
  force_fields_input.F
  force_fields_util.F
  force_field_types.F
  fp_methods.F
  fp_types.F
  free_energy_types.F
  generic_os_integrals.F
  gle_system_types.F
  global_types.F
  graphcon.F
  grid/grid_api.F
  group_dist_types.F
  grrm_utils.F
  gapw_1c_basis_set.F
  gw_communication.F
  gw_methods.F
  gw_types.F
  gw_utils.F
  hartree_local_methods.F
  hartree_local_types.F
  header.F
  hfx_admm_utils.F
  hfx_communication.F
  hfx_compression_methods.F
  hfx_derivatives.F
  hfx_energy_potential.F
  hfx_exx.F
  hfx_helpers.F
  hfx_libint_interface.F
  hfx_load_balance_methods.F
  hfx_pair_list_methods.F
  hfx_pw_methods.F
  hfx_ri.F
  hfx_screening_methods.F
  hfx_types.F
  hirshfeld_methods.F
  hirshfeld_types.F
  input_constants.F
  input_cp2k_almo.F
  input_cp2k_atom.F
  input_cp2k_atprop.F
  input_cp2k_binary_restarts.F
  input_cp2k_check.F
  input_cp2k_colvar.F
  input_cp2k_constraints.F
  input_cp2k_dft.F
  input_cp2k_distribution.F
  input_cp2k_ec.F
  input_cp2k_eip.F
  input_cp2k_embed.F
  input_cp2k_exstate.F
  input_cp2k_external.F
  input_cp2k_field.F
  input_cp2k_force_eval.F
  input_cp2k_free_energy.F
  input_cp2k_global.F
  input_cp2k_hfx.F
  input_cp2k_kpoints.F
  input_cp2k_loc.F
  input_cp2k_ls.F
  input_cp2k_mixed.F
  input_cp2k_mm.F
  input_cp2k_motion_print.F
  input_cp2k_mp2.F
  input_cp2k_neb.F
  input_cp2k_negf.F
  input_cp2k_nnp.F
  input_cp2k_poisson.F
  input_cp2k_projection_rtp.F
  input_cp2k_properties_dft.F
  input_cp2k_pwdft.F
  input_cp2k_qmmm.F
  input_cp2k_read.F
  input_cp2k_resp.F
  input_cp2k_restarts_util.F
  input_cp2k_rsgrid.F
  input_cp2k_subsys.F
  input_cp2k_tb.F
  input_cp2k_thermostats.F
  input_cp2k_transport.F
  input_cp2k_voronoi.F
  input_cp2k_xc.F
  input_optimize_basis.F
  input_optimize_input.F
  input_restart_force_eval.F
  input_restart_rng.F
  integration_grid_types.F
  iao_analysis.F
  iao_types.F
  ipi_driver.F
  iterate_matrix.F
  kg_correction.F
  kg_environment.F
  kg_environment_types.F
  kg_tnadd_mat.F
  kg_vertex_coloring_methods.F
  kpoint_coulomb_2c.F
  kpoint_io.F
  kpoint_methods.F
  kpoint_transitional.F
  kpoint_types.F
  libint_2c_3c.F
  libint_wrapper.F
  library_tests.F
  linesearch.F
  localization_tb.F
  lri_compression.F
  lri_environment_init.F
  lri_environment_methods.F
  lri_environment_types.F
  lri_forces.F
  lri_integrals.F
  lri_ks_methods.F
  lri_optimize_ri_basis.F
  lri_optimize_ri_basis_types.F
  ls_matrix_exp.F
  manybody_eam.F
  manybody_gal21.F
  manybody_gal.F
  manybody_allegro.F
  manybody_nequip.F
  manybody_potential.F
  manybody_quip.F
  manybody_siepmann.F
  manybody_tersoff.F
  mao_basis.F
  mao_methods.F
  mao_optimizer.F
  mao_types.F
  mao_wfn_analysis.F
  matrix_exp.F
  maxwell_solver_interface.F
  mdctrl_types.F
  metadynamics.F
  metadynamics_types.F
  metadynamics_utils.F
  metadyn_tools/graph_methods.F
  metadyn_tools/graph_utils.F
  min_basis_set.F
  minbas_methods.F
  minbas_wfn_analysis.F
  mixed_cdft_methods.F
  mixed_cdft_types.F
  mixed_cdft_utils.F
  mixed_energy_types.F
  mixed_environment.F
  mixed_environment_types.F
  mixed_environment_utils.F
  mixed_main.F
  mm_collocate_potential.F
  mm_mapping_library.F
  mode_selective.F
  molden_utils.F
  molecular_dipoles.F
  molecular_moments.F
  molecular_states.F
  mol_force.F
  molsym.F
  moments_utils.F
  motion_utils.F
  mp2_cphf.F
  mp2_direct_method.F
  mp2_eri.F
  mp2_eri_gpw.F
  mp2.F
  mp2_gpw.F
  mp2_gpw_method.F
  mp2_grids.F
  mp2_integrals.F
  mp2_laplace.F
  mp2_optimize_ri_basis.F
  mp2_ri_2c.F
  mp2_ri_gpw.F
  mp2_ri_grad.F
  mp2_ri_grad_util.F
  mp2_ri_libint.F
  mp2_setup.F
  mp2_types.F
  mscfg_methods.F
  mscfg_types.F
  mulliken.F
  negf_alloc_types.F
  negf_atom_map.F
  negf_control_types.F
  negf_env_types.F
  negf_green_cache.F
  negf_green_methods.F
  negf_integr_cc.F
  negf_integr_simpson.F
  negf_integr_utils.F
  negf_matrix_utils.F
  negf_methods.F
  negf_subgroup_types.F
  negf_vectors.F
  nnp_acsf.F
  nnp_environment.F
  nnp_environment_types.F
  nnp_force.F
  nnp_model.F
  local_gemm_api.F
  offload/offload_api.F
  optbas_fenv_manipulation.F
  optbas_opt_utils.F
  optimize_basis.F
  optimize_basis_types.F
  optimize_basis_utils.F
  optimize_dmfet_potential.F
  optimize_embedding_potential.F
  optimize_input.F
  outer_scf_control_types.F
  paw_proj_set_types.F
  pair_potential_coulomb.F
  pair_potential.F
  pair_potential_types.F
  pair_potential_util.F
  pao_input.F
  pao_io.F
  pao_linpot_full.F
  pao_linpot_rotinv.F
  pao_main.F
  pao_methods.F
  pao_ml_descriptor.F
  pao_ml.F
  pao_ml_gaussprocess.F
  pao_ml_neuralnet.F
  pao_optimizer.F
  pao_param_exp.F
  pao_param.F
  pao_param_equi.F
  pao_param_fock.F
  pao_param_gth.F
  pao_param_linpot.F
  pao_param_methods.F
  pao_potentials.F
  pao_types.F
  particle_methods.F
  paw_proj_set_types.F
  pexsi_interface.F
  pexsi_methods.F
  pexsi_types.F
  pilaenv_hack.F
  pme.F
  pme_tools.F
  population_analyses.F
  preconditioner_apply.F
  preconditioner.F
  preconditioner_makes.F
  preconditioner_solvers.F
  preconditioner_types.F
  qmmm_create.F
  qmmm_elpot.F
  qmmm_ff_fist.F
  qmmm_force.F
  qmmm_gaussian_data.F
  qmmm_gaussian_init.F
  qmmm_gaussian_input.F
  qmmm_gaussian_types.F
  qmmm_gpw_energy.F
  qmmm_gpw_forces.F
  qmmm_image_charge.F
  qmmm_init.F
  qmmm_links_methods.F
  qmmm_per_elpot.F
  qmmm_pw_grid.F
  qmmm_se_energy.F
  qmmm_se_forces.F
  qmmm_tb_coulomb.F
  qmmm_tb_methods.F
  qmmm_topology_util.F
  qmmm_types.F
  qmmm_types_low.F
  qmmm_util.F
  qmmmx_create.F
  qmmmx_force.F
  qmmmx_types.F
  qmmmx_update.F
  qmmmx_util.F
  qs_2nd_kernel_ao.F
  qs_active_space_methods.F
  qs_active_space_types.F
  qs_atomic_block.F
  qs_band_structure.F
  qs_basis_gradient.F
  qs_block_davidson_types.F
  qs_cdft_methods.F
  qs_cdft_opt_types.F
  qs_cdft_scf_utils.F
  qs_cdft_types.F
  qs_cdft_utils.F
  qs_charge_mixing.F
  qs_chargemol.F
  qs_charges_types.F
  qs_collocate_density.F
  qs_commutators.F
  qs_condnum.F
  qs_core_energies.F
  qs_core_hamiltonian.F
  qs_dcdr_ao.F
  qs_dcdr.F
  qs_dcdr_utils.F
  qs_density_matrices.F
  qs_density_mixing_types.F
  qs_dftb3_methods.F
  qs_dftb_coulomb.F
  qs_dftb_dispersion.F
  qs_dftb_matrices.F
  qs_dftb_parameters.F
  qs_dftb_types.F
  qs_dftb_utils.F
  qs_diis.F
  qs_diis_types.F
  qs_dispersion_nonloc.F
  qs_dispersion_pairpot.F
  qs_dispersion_types.F
  qs_dispersion_utils.F
  qs_dos.F
  qs_efield_berry.F
  qs_efield_local.F
  qs_elec_field.F
  qs_electric_field_gradient.F
  qs_elf_methods.F
  qs_energy.F
  qs_energy_init.F
  qs_energy_matrix_w.F
  qs_energy_types.F
  qs_energy_utils.F
  qs_energy_window.F
  qs_environment.F
  qs_environment_methods.F
  qs_environment_types.F
  qs_epr_hyp.F
  qs_external_density.F
  qs_external_potential.F
  qs_fb_atomic_halo_types.F
  qs_fb_atomic_matrix_methods.F
  qs_fb_buffer_types.F
  qs_fb_com_tasks_types.F
  qs_fb_distribution_methods.F
  qs_fb_env_methods.F
  qs_fb_env_types.F
  qs_fb_filter_matrix_methods.F
  qs_fb_hash_table_types.F
  qs_fb_input.F
  qs_fb_matrix_data_types.F
  qs_fb_trial_fns_types.F
  qs_fermi_contact.F
  qs_force.F
  qs_force_types.F
  qs_fxc.F
  qs_gamma2kp.F
  qs_gapw_densities.F
  qs_gcp_method.F
  qs_gcp_types.F
  qs_gcp_utils.F
  qs_grid_atom.F
  qs_gspace_mixing.F
  qs_harmonics_atom.F
  qs_hash_table_functions.F
  qs_initial_guess.F
  qs_integral_utils.F
  qs_integrate_potential.F
  qs_integrate_potential_product.F
  qs_integrate_potential_single.F
  qs_interactions.F
  qs_kernel_methods.F
  qs_kernel_types.F
  qs_kind_types.F
  qs_kinetic.F
  qs_kpp1_env_methods.F
  qs_kpp1_env_types.F
  qs_ks_apply_restraints.F
  qs_ks_atom.F
  qs_ks_methods.F
  qs_ks_qmmm_methods.F
  qs_ks_qmmm_types.F
  qs_ks_types.F
  qs_ks_utils.F
  qs_ks_reference.F
  qs_linres_atom_current.F
  qs_linres_current.F
  qs_linres_current_utils.F
  qs_linres_epr_nablavks.F
  qs_linres_epr_ownutils.F
  qs_linres_epr_utils.F
  qs_linres_issc_utils.F
  qs_linres_kernel.F
  qs_linres_methods.F
  qs_linres_module.F
  qs_linres_nmr_epr_common_utils.F
  qs_linres_nmr_shift.F
  qs_linres_nmr_utils.F
  qs_linres_op.F
  qs_linres_polar_utils.F
  qs_linres_types.F
  qs_localization_methods.F
  qs_local_properties.F
  qs_local_rho_types.F
  qs_loc_dipole.F
  qs_loc_main.F
  qs_loc_methods.F
  qs_loc_molecules.F
  qs_loc_states.F
  qs_loc_types.F
  qs_loc_utils.F
  qs_matrix_pools.F
  qs_matrix_w.F
  qs_mfp.F
  qs_mixing_utils.F
  qs_mo_io.F
  qs_moments.F
  qs_mo_methods.F
  qs_mom_methods.F
  qs_mom_types.F
  qs_mo_occupation.F
  qs_mo_types.F
  qs_neighbor_lists.F
  qs_neighbor_list_types.F
  qs_nl_hash_table_types.F
  qs_o3c_methods.F
  qs_o3c_types.F
  qs_oce_methods.F
  qs_oce_types.F
  qs_operators_ao.F
  qs_ot_eigensolver.F
  qs_ot.F
  qs_ot_minimizer.F
  qs_ot_scf.F
  qs_ot_types.F
  qs_outer_scf.F
  qs_overlap.F
  qs_pdos.F
  qs_p_env_methods.F
  qs_p_env_types.F
  qs_period_efield_types.F
  qs_resp.F
  qs_rho0_ggrid.F
  qs_rho0_methods.F
  qs_rho0_types.F
  qs_rho_atom_methods.F
  qs_rho_atom_types.F
  qs_rho_methods.F
  qs_rho_types.F
  qs_sccs.F
  qs_scf_block_davidson.F
  qs_scf_csr_write.F
  qs_scf_diagonalization.F
  qs_scf.F
  qs_scf_initialization.F
  qs_scf_lanczos.F
  qs_scf_loop_utils.F
  qs_scf_methods.F
  qs_scf_output.F
  qs_scf_post_gpw.F
  qs_scf_post_scf.F
  qs_scf_post_se.F
  qs_scf_post_tb.F
  qs_scf_types.F
  qs_scf_wfn_mix.F
  qs_spin_orbit.F
  qs_subsys_methods.F
  qs_subsys_types.F
  qs_tddfpt2_assign.F
  qs_tddfpt2_densities.F
  qs_tddfpt2_eigensolver.F
  qs_tddfpt2_fhxc.F
  qs_tddfpt2_fhxc_forces.F
  qs_tddfpt2_forces.F
  qs_tddfpt2_lri_utils.F
  qs_tddfpt2_methods.F
  qs_tddfpt2_operators.F
  qs_tddfpt2_properties.F
  qs_tddfpt2_restart.F
  qs_tddfpt2_soc.F
  qs_tddfpt2_soc_utils.F
  qs_tddfpt2_soc_types.F
  qs_tddfpt2_stda_types.F
  qs_tddfpt2_stda_utils.F
  qs_tddfpt2_subgroups.F
  qs_tddfpt2_types.F
  qs_tddfpt2_utils.F
  qs_tddfpt_eigensolver.F
  qs_tddfpt_module.F
  qs_tddfpt_types.F
  qs_tddfpt_utils.F
  qs_tddfpt2_fprint.F
  qs_tensors.F
  qs_tensors_types.F
  qs_update_s_mstruct.F
  qs_vcd.F
  qs_vcd_ao.F
  qs_vcd_utils.F
  qs_vxc_atom.F
  qs_vxc.F
  qs_wannier90.F
  qs_wf_history_methods.F
  qs_wf_history_types.F
  rel_control_types.F
  replica_methods.F
  replica_types.F
  response_solver.F
  restraint.F
  ri_environment_methods.F
  rmsd.F
  rpa_axk.F
  rpa_communication.F
  rpa_grad.F
  rpa_gw.F
  rpa_gw_ic.F
  rpa_gw_im_time_util.F
  rpa_gw_kpoints_util.F
  rpa_gw_sigma_x.F
  rpa_im_time.F
  rpa_main.F
  rpa_rse.F
  rpa_im_time_force_types.F
  rpa_im_time_force_methods.F
  rpa_util.F
  rtp_admm_methods.F
  rt_propagation_forces.F
  rt_propagation_types.F
  rt_propagation_velocity_gauge.F
  sap_kind_types.F
  scf_control_types.F
  scine_utils.F
  se_core_core.F
  se_core_matrix.F
  se_fock_matrix_coulomb.F
  se_fock_matrix_dbg.F
  se_fock_matrix_exchange.F
  se_fock_matrix.F
  se_fock_matrix_integrals.F
  semi_empirical_expns3_methods.F
  semi_empirical_expns3_types.F
  semi_empirical_int3_utils.F
  semi_empirical_int_ana.F
  semi_empirical_int_arrays.F
  semi_empirical_int_debug.F
  semi_empirical_integrals.F
  semi_empirical_int_gks.F
  semi_empirical_int_num.F
  semi_empirical_int_utils.F
  semi_empirical_mpole_methods.F
  semi_empirical_mpole_types.F
  semi_empirical_parameters.F
  semi_empirical_par_utils.F
  semi_empirical_store_int_types.F
  semi_empirical_types.F
  semi_empirical_utils.F
  simpar_types.F
  sirius_interface.F
  spglib_f08.F
  splines_methods.F
  splines_types.F
  spme.F
  s_square_methods.F
  statistical_methods.F
  stm_images.F
  subcell_types.F
  submatrix_dissection.F
  submatrix_types.F
  surface_dipole.F
  taper_types.F
  task_list_methods.F
  task_list_types.F
  tip_scan_methods.F
  tip_scan_types.F
  topology_amber.F
  topology_cif.F
  topology_connectivity_util.F
  topology_constraint_util.F
  topology_coordinate_util.F
  topology_cp2k.F
  topology.F
  topology_generate_util.F
  topology_gromos.F
  topology_input.F
  topology_multiple_unit_cell.F
  topology_pdb.F
  topology_psf.F
  topology_types.F
  topology_util.F
  topology_xtl.F
  topology_xyz.F
  torch_api.F
  transport_env_types.F
  transport.F
  uff_vdw_radii_table.F
  virial_methods.F
  voronoi_interface.F
  wannier90.F
  wannier_states.F
  wannier_states_types.F
  xas_control.F
  xas_env_types.F
  xas_methods.F
  xas_restart.F
  xas_tdp_atom.F
  xas_tdp_correction.F
  xas_tdp_integrals.F
  xas_tdp_kernel.F
  xas_tdp_methods.F
  xas_tdp_types.F
  xas_tdp_utils.F
  xas_tp_scf.F
  xray_diffraction.F
  xtb_coulomb.F
  xtb_ehess.F
  xtb_ehess_force.F
  xtb_matrices.F
  xtb_parameters.F
  xtb_types.F)

list(
  APPEND
  CP2K_SRCS_F
  aobasis/ai_angmom.F
  aobasis/ai_contraction.F
  aobasis/ai_contraction_sphi.F
  aobasis/ai_coulomb.F
  aobasis/ai_coulomb_test.F
  aobasis/ai_derivatives.F
  aobasis/ai_elec_field.F
  aobasis/ai_eri_debug.F
  aobasis/ai_fermi_contact.F
  aobasis/ai_kinetic.F
  aobasis/ai_moments.F
  aobasis/ai_onecenter.F
  aobasis/ai_oneelectron.F
  aobasis/ai_operator_ra2m.F
  aobasis/ai_operators_r12.F
  aobasis/ai_os_rr.F
  aobasis/ai_overlap3_debug.F
  aobasis/ai_overlap3.F
  aobasis/ai_overlap_aabb.F
  aobasis/ai_overlap_debug.F
  aobasis/ai_overlap.F
  aobasis/ai_overlap_ppl.F
  aobasis/ai_spin_orbit.F
  aobasis/ai_verfc.F
  aobasis/ao_util.F
  aobasis/aux_basis_set.F
  aobasis/basis_set_container_types.F
  aobasis/basis_set_types.F
  aobasis/orbital_pointers.F
  aobasis/orbital_symbols.F
  aobasis/orbital_transformation_matrices.F
  aobasis/paw_basis_types.F
  aobasis/soft_basis_set.F
  aobasis/sto_ng.F)

list(
  APPEND
  CP2K_SRCS_F
  arnoldi/arnoldi_api.F
  arnoldi/arnoldi_data_methods.F
  arnoldi/arnoldi_geev.F
  arnoldi/arnoldi_methods.F
  arnoldi/arnoldi_types.F)

list(APPEND CP2K_SRCS_F base/base_hooks.F base/kinds.F base/machine.F)

list(
  APPEND
  CP2K_SRCS_F
  mpiwrap/mp_perf_env.F
  common/bessel_lib.F
  common/bibliography.F
  common/callgraph.F
  common/cg_test.F
  common/cp_array_sort.F
  common/cp_array_utils.F
  common/cp_error_handling.F
  common/cp_files.F
  common/cp_iter_types.F
  common/cp_log_handling.F
  common/cp_min_heap.F
  common/cp_result_methods.F
  common/cp_result_types.F
  common/cp_units.F
  common/d3_poly.F
  common/distribution_1d_types.F
  common/eigenvalueproblems.F
  common/fparser.F
  common/gamma.F
  common/kahan_sum.F
  common/lapack.F
  common/lebedev.F
  common/linear_systems.F
  common/list_callstackentry.F
  common/list.F
  common/list_routinereport.F
  common/list_routinestat.F
  common/list_timerenv.F
  common/mathconstants.F
  common/mathlib.F
  common/memory_utilities.F
  common/parallel_rng_types.F
  common/periodic_table.F
  common/physcon.F
  common/powell.F
  common/print_messages.F
  common/reference_manager.F
  common/routine_map.F
  common/spherical_harmonics.F
  common/splines.F
  common/string_table.F
  common/string_utilities.F
  common/structure_factors.F
  common/structure_factor_types.F
  common/t_c_g0.F
  common/timings_base_type.F
  common/timings.F
  common/timings_report.F
  common/timings_types.F
  common/t_sh_p_s_c.F
  common/util.F
  common/whittaker.F)

list(APPEND CP2K_SRCS_F dbcsrx/dbcsr_vector.F)

list(APPEND CP2K_SRCS_F dbm/dbm_api.F dbm/dbm_tests.F)

list(
  APPEND
  CP2K_SRCS_F
  dbt/dbt_allocate_wrap.F
  dbt/dbt_api.F
  dbt/dbt_array_list_methods.F
  dbt/dbt_block.F
  dbt/dbt_index.F
  dbt/dbt_io.F
  dbt/dbt_methods.F
  dbt/dbt_reshape_ops.F
  dbt/dbt_split.F
  dbt/dbt_test.F
  dbt/dbt_types.F
  dbt/tas/dbt_tas_base.F
  dbt/tas/dbt_tas_global.F
  dbt/tas/dbt_tas_io.F
  dbt/tas/dbt_tas_mm.F
  dbt/tas/dbt_tas_reshape_ops.F
  dbt/tas/dbt_tas_split.F
  dbt/tas/dbt_tas_test.F
  dbt/tas/dbt_tas_types.F
  dbt/tas/dbt_tas_util.F)

list(
  APPEND
  CP2K_SRCS_F
  emd/rt_delta_pulse.F
  emd/rt_hfx_utils.F
  emd/rt_make_propagators.F
  emd/rt_projection_mo_utils.F
  emd/rt_propagation_methods.F
  emd/rt_propagation_output.F
  emd/rt_propagation_utils.F
  emd/rt_propagator_init.F)

list(
  APPEND
  CP2K_SRCS_F
  fm/cp_blacs_env.F
  fm/cp_blacs_types.F
  fm/cp_cfm_basic_linalg.F
  fm/cp_cfm_diag.F
  fm/cp_cfm_basic_linalg.F
  fm/cp_cfm_types.F
  fm/cp_fm_basic_linalg.F
  fm/cp_fm_cholesky.F
  fm/cp_fm_cusolver_api.F
  fm/cp_fm_diag.F
  fm/cp_fm_diag_utils.F
  fm/cp_fm_elpa.F
  fm/cp_fm_pool_types.F
  fm/cp_fm_types.F
  fm/cp_fm_struct.F
  fm/cp_linked_list_fm.F)

list(APPEND CP2K_SRCS_C fm/cp_fm_cusolver.c)

list(APPEND CP2K_SRCS_F hfxbase/hfx_compression_core_methods.F
     hfxbase/hfx_contract_block.F hfxbase/hfx_contraction_methods.F)

list(
  APPEND
  CP2K_SRCS_F
  input/cp_linked_list_input.F
  input/cp_output_handling.F
  input/cp_parser_buffer_types.F
  input/cp_parser_ilist_methods.F
  input/cp_parser_ilist_types.F
  input/cp_parser_inpp_methods.F
  input/cp_parser_inpp_types.F
  input/cp_parser_methods.F
  input/cp_parser_status_types.F
  input/cp_parser_types.F
  input/input_enumeration_types.F
  input/input_keyword_types.F
  input/input_parsing.F
  input/input_section_types.F
  input/input_val_types.F)

list(APPEND CP2K_SRCS_F mpiwrap/mp_perf_env.F mpiwrap/message_passing.F
     mpiwrap/mp_perf_test.F)

list(
  APPEND
  CP2K_SRCS_F
  minimax/minimax_exp.F
  minimax/minimax_exp_gw.F
  minimax/minimax_exp_k15.F
  minimax/minimax_exp_k53.F
  minimax/minimax_rpa.F)

list(
  APPEND
  CP2K_SRCS_F
  motion/averages_types.F
  motion/bfgs_optimizer.F
  motion/cell_opt.F
  motion/cell_opt_types.F
  motion/cell_opt_utils.F
  motion/cg_optimizer.F
  motion/cg_utils.F
  motion/cp_lbfgs.F
  motion/cp_lbfgs_geo.F
  motion/cp_lbfgs_optimizer_gopt.F
  motion/dimer_methods.F
  motion/dimer_types.F
  motion/dimer_utils.F
  motion/free_energy_methods.F
  motion/geo_opt.F
  motion/glbopt_callback.F
  motion/gopt_f77_methods.F
  motion/gopt_f_methods.F
  motion/gopt_f_types.F
  motion/gopt_param_types.F
  motion/helium_common.F
  motion/helium_interactions.F
  motion/helium_io.F
  motion/helium_methods.F
  motion/helium_sampling.F
  motion/helium_types.F
  motion/helium_worm.F
  motion/input_cp2k_md.F
  motion/input_cp2k_restarts.F
  motion/input_cp2k_vib.F
  motion/integrator.F
  motion/integrator_utils.F
  motion/mc/mc_control.F
  motion/mc/mc_coordinates.F
  motion/mc/mc_ensembles.F
  motion/mc/mc_environment_types.F
  motion/mc/mc_ge_moves.F
  motion/mc/mc_misc.F
  motion/mc/mc_move_control.F
  motion/mc/mc_moves.F
  motion/mc/mc_run.F
  motion/mc/mc_types.F
  motion/mc/tamc_run.F
  motion/md_conserved_quantities.F
  motion/mdctrl_methods.F
  motion/md_energies.F
  motion/md_ener_types.F
  motion/md_environment_types.F
  motion/md_run.F
  motion/md_util.F
  motion/md_vel_utils.F
  motion/neb_io.F
  motion/neb_md_utils.F
  motion/neb_methods.F
  motion/neb_opt_utils.F
  motion/neb_types.F
  motion/neb_utils.F
  motion/pint_gle.F
  motion/pint_io.F
  motion/pint_methods.F
  motion/pint_normalmode.F
  motion/pint_piglet.F
  motion/pint_pile.F
  motion/pint_public.F
  motion/pint_qtb.F
  motion/pint_staging.F
  motion/pint_transformations.F
  motion/pint_types.F
  motion/reftraj_types.F
  motion/reftraj_util.F
  motion/rt_propagation.F
  motion/shell_opt.F
  motion/simpar_methods.F
  motion/space_groups.F
  motion/space_groups_types.F
  motion/thermal_region_types.F
  motion/thermal_region_utils.F
  motion/thermostat/al_system_dynamics.F
  motion/thermostat/al_system_init.F
  motion/thermostat/al_system_mapping.F
  motion/thermostat/barostat_types.F
  motion/thermostat/barostat_utils.F
  motion/thermostat/csvr_system_dynamics.F
  motion/thermostat/csvr_system_init.F
  motion/thermostat/csvr_system_mapping.F
  motion/thermostat/extended_system_dynamics.F
  motion/thermostat/extended_system_init.F
  motion/thermostat/extended_system_mapping.F
  motion/thermostat/gle_system_dynamics.F
  motion/thermostat/input_cp2k_barostats.F
  motion/thermostat/thermostat_mapping.F
  motion/thermostat/thermostat_methods.F
  motion/thermostat/thermostat_types.F
  motion/thermostat/thermostat_utils.F
  motion/velocity_verlet_control.F
  motion/vibrational_analysis.F
  motion/wiener_process.F)

list(
  APPEND
  CP2K_SRCS_F
  pw/cp_linked_list_pw.F
  pw/cube_utils.F
  pw/dct.F
  pwdft_environment.F
  pwdft_environment_types.F
  pw/dg_rho0_types.F
  pw/dgs.F
  pw/dg_types.F
  pw/dielectric_methods.F
  pw/dielectric_types.F
  pw/dirichlet_bc_methods.F
  pw/dirichlet_bc_types.F
  pw_env/cp_spline_utils.F
  pw_env/gaussian_gridlevels.F
  pw_env_methods.F
  pw_env/pw_env_types.F
  pw_env/rs_pw_interface.F
  pw/fast.F
  pw/fft/fft_kinds.F
  pw/fft/fft_lib.F
  pw/fft/fft_plan.F
  pw/fft/fftsg_lib.F
  pw/fft/fftw3_lib.F
  pw/fft/mltfftsg_tools.F
  pw/fft_tools.F
  pw/lazy.F
  pw/mt_util.F
  pw_poisson_read_input.F
  pw/ps_implicit_methods.F
  pw/ps_implicit_types.F
  pw/ps_wavelet_base.F
  pw/ps_wavelet_fft3d.F
  pw/ps_wavelet_kernel.F
  pw/ps_wavelet_methods.F
  pw/ps_wavelet_scaling_function.F
  pw/ps_wavelet_types.F
  pw/ps_wavelet_util.F
  pw/pw_copy_all.F
  pw/pw_fpga.F
  pw/pw_grid_info.F
  pw/pw_grids.F
  pw/pw_gpu.F
  pw/pw_grid_types.F
  pw/pw_methods.F
  pw/pw_poisson_methods.F
  pw/pw_poisson_types.F
  pw/pw_pool_types.F
  pw/pw_spline_utils.F
  pw/pw_types.F
  pw/realspace_grid_cube.F
  pw/realspace_grid_types.F
  pw/rs_methods.F)

if(CP2K_USE_PW_GPU)
  list(APPEND CP2K_SRCS_F pw/pw_gpu.F)
endif()

list(
  APPEND
  CP2K_SRCS_F
  shg_int/constants_operator.F
  shg_int/construct_shg.F
  shg_integrals_test.F
  shg_int/generic_shg_integrals.F
  shg_int/generic_shg_integrals_init.F
  shg_int/s_contract_shg.F)

list(
  APPEND
  CP2K_SRCS_F
  start/cp2k_runs.F
  start/cp2k_shell.F
  start/input_cp2k.F
  start/input_cp2k_motion.F
  start/libcp2k.F)

list(
  APPEND
  CP2K_SRCS_F
  subsys/atomic_kind_list_types.F
  subsys/atomic_kind_types.F
  subsys/atprop_types.F
  subsys/cell_types.F
  subsys/colvar_types.F
  subsys/cp_subsys_types.F
  subsys/damping_dipole_types.F
  subsys/external_potential_types.F
  subsys/force_field_kind_types.F
  subsys/molecule_kind_list_types.F
  subsys/molecule_kind_types.F
  subsys/molecule_list_types.F
  subsys/molecule_types.F
  subsys/multipole_types.F
  subsys/particle_list_types.F
  subsys/particle_types.F
  subsys/shell_potential_types.F
  subsys/virial_types.F)

list(
  APPEND
  CP2K_SRCS_F
  swarm/glbopt_history.F
  swarm/glbopt_input.F
  swarm/glbopt_master.F
  swarm/glbopt_mincrawl.F
  swarm/glbopt_minhop.F
  swarm/glbopt_worker.F
  swarm/swarm.F
  swarm/swarm_input.F
  swarm/swarm_master.F
  swarm/swarm_message.F
  swarm/swarm_mpi.F
  swarm/swarm_worker.F)

list(
  APPEND
  CP2K_SRCS_F
  tmc/input_cp2k_tmc.F
  tmc/tmc_analysis.F
  tmc/tmc_analysis_types.F
  tmc/tmc_calculations.F
  tmc/tmc_cancelation.F
  tmc/tmc_dot_tree.F
  tmc/tmc_file_io.F
  tmc/tmc_master.F
  tmc/tmc_messages.F
  tmc/tmc_move_handle.F
  tmc/tmc_moves.F
  tmc/tmc_move_types.F
  tmc/tmc_setup.F
  tmc/tmc_stati.F
  tmc/tmc_tree_acceptance.F
  tmc/tmc_tree_build.F
  tmc/tmc_tree_references.F
  tmc/tmc_tree_search.F
  tmc/tmc_tree_types.F
  tmc/tmc_types.F
  tmc/tmc_utils.F
  tmc/tmc_worker.F)

list(
  APPEND
  CP2K_SRCS_F
  xc_adiabatic_methods.F
  xc_adiabatic_utils.F
  xc/cp_linked_list_xc_deriv.F
  xc_pot_saop.F
  xc_write_output.F
  xc/xc_atom.F
  xc/xc_b97.F
  xc/xc_b97_fxc.F
  xc/xc_cs1.F
  xc/xc_derivative_desc.F
  xc/xc_derivative_set_types.F
  xc/xc_derivatives.F
  xc/xc_derivative_types.F
  xc/xc_exchange_gga.F
  xc/xc.F
  xc/xc_functionals_utilities.F
  xc/xc_fxc_kernel.F
  xc/xc_hcth.F
  xc/xc_input_constants.F
  xc/xc_ke_gga.F
  xc/xc_libxc.F
  xc/xc_libxc_wrap.F
  xc/xc_lyp_adiabatic.F
  xc/xc_lyp.F
  xc/xc_optx.F
  xc/xc_pade.F
  xc/xc_pbe.F
  xc/xc_perdew86.F
  xc/xc_perdew_wang.F
  xc/xc_perdew_zunger.F
  xc/xc_rho_cflags_types.F
  xc/xc_rho_set_types.F
  xc/xc_tfw.F
  xc/xc_thomas_fermi.F
  xc/xc_tpss.F
  xc/xc_util.F
  xc/xc_vwn.F
  xc/xc_xalpha.F
  xc/xc_xbecke88.F
  xc/xc_xbecke88_long_range.F
  xc/xc_xbecke88_lr_adiabatic.F
  xc/xc_xbecke_roussel.F
  xc/xc_xbeef.F
  xc/xc_xbr_pbe_lda_hole_t_c_lr.F
  xc/xc_xlda_hole_t_c_lr.F
  xc/xc_xpbe_hole_t_c_lr.F
  xc/xc_xwpbe.F)

list(APPEND CP2K_SRCS_C sockets.c base/machine_cpuid.c)

set(CP2K_DBM_SRCS_C
    dbm/dbm_distribution.c
    dbm/dbm_library.c
    dbm/dbm_matrix.c
    dbm/dbm_mempool.c
    dbm/dbm_mpi.c
    dbm/dbm_multiply.c
    dbm/dbm_multiply_comm.c
    dbm/dbm_multiply_cpu.c
    dbm/dbm_shard.c)

set(CP2K_GRID_SRCS_C
    grid/common/grid_basis_set.c
    grid/common/grid_library.c
    grid/common/grid_sphere_cache.c
    grid/dgemm/grid_dgemm_coefficients.c
    grid/dgemm/grid_dgemm_collocation_integration.c
    grid/dgemm/grid_dgemm_collocate.c
    grid/dgemm/grid_dgemm_context.c
    grid/dgemm/grid_dgemm_integrate.c
    grid/dgemm/grid_dgemm_prepare_pab.c
    grid/dgemm/grid_dgemm_non_orthorombic_corrections.c
    grid/dgemm/grid_dgemm_tensor_local.c
    grid/dgemm/grid_dgemm_utils.c
    grid/grid_task_list.c
    grid/ref/grid_ref_collocate.c
    grid/ref/grid_ref_integrate.c
    grid/ref/grid_ref_prepare_pab.c
    grid/ref/grid_ref_task_list.c
    grid/cpu/grid_cpu_collocate.c
    grid/cpu/grid_cpu_integrate.c
    grid/cpu/grid_cpu_prepare_pab.c
    grid/cpu/grid_cpu_task_list.c
    grid/grid_replay.c)

set(CP2K_FPGA_SRC_C pw/fpga/fft_fpga.c pw/fpga/opencl_utils.c)
set(CP2K_PW_SRCS_C pw/gpu/pw_gpu_internal.c)
set(CP2K_OFFLOAD_SRCS_C offload/offload_buffer.c offload/offload_library.c)
set(CP2K_SRCS_CPP torch_c_api.cpp)

list(APPEND CP2K_SRCS_C "${CP2K_OFFLOAD_SRCS_C}")

# ##############################################################################
# gpu code                                     #
# ##############################################################################

set(CP2K_SRCS_GPU "")

set(CP2K_GRID_SRCS_GPU
    grid/gpu/grid_gpu_collocate.cu grid/gpu/grid_gpu_integrate.cu
    grid/gpu/grid_gpu_task_list.cu)

set(CP2K_GRID_SRCS_HIP
    grid/hip/grid_hip_collocate.cu grid/hip/grid_hip_integrate.cu
    grid/hip/grid_hip_context.cu)

set(CP2K_DBM_SRCS_GPU dbm/dbm_multiply_gpu_kernel.cu)
set(CP2K_DBM_SRCS_GPU_C dbm/dbm_multiply_gpu.c)

set(CP2K_PW_SRCS_GPU pw/gpu/pw_gpu_kernels.cu)

set(CP2K_PROGS_C
    grid/grid_miniapp.c;grid/grid_unittest.c;dbm/dbm_miniapp.c;start/libcp2k_unittest.c
)
set(CP2K_PROGS_F
    dbt/dbt_unittest.F
    dbt/tas/dbt_tas_unittest.F
    start/cp2k.F
    common/memory_utilities_unittest.F
    common/parallel_rng_types_unittest.F
    metadyn_tools/graph.F
    motion/dumpdcd.F
    motion/xyz2dcd.F
    nequip_unittest.F)

if(CP2K_USE_CUDA OR CP2K_USE_HIP)
  if(NOT CP2K_DISABLE_PW_GPU)
    list(APPEND CP2K_SRCS_GPU "pw/gpu/pw_gpu_kernels.cu")
    list(APPEND CP2K_SRCS_C ${CP2K_PW_SRCS_C})
  endif()

  if(NOT CP2K_DISABLE_DBM_GPU)
    list(APPEND CP2K_DBM_SRCS_C ${CP2K_DBM_SRCS_GPU_C})
    list(APPEND CP2K_SRCS_GPU ${CP2K_DBM_SRCS_GPU})
  endif()

  if(NOT CP2K_DISABLE_GRID_GPU)
    if(CP2K_USE_HIP)
      list(APPEND CP2K_SRCS_GPU ${CP2K_GRID_SRCS_HIP})
    endif()
    list(APPEND CP2K_SRCS_GPU ${CP2K_GRID_SRCS_GPU})
  endif()
  if(CP2K_USE_HIP)
    set_source_files_properties(${CP2K_SRCS_GPU} PROPERTIES LANGUAGE HIP)
  endif()
endif()

list(APPEND CP2K_SRCS_C ${CP2K_DBM_SRCS_C} ${CP2K_GRID_SRCS_C})

# ##############################################################################
#
# ##############################################################################

# check that the files registered in CMakeLists.txt are all present and return
# an error otherwise. cmake only runs this test nothing else below

if(CP2K_ENABLE_CONSISTENCY_CHECKS)
  cp2k_compare_src_with_list(
    "${CP2K_FPGA_SRC_C};${CP2K_SRCS_C};${CP2K_OFFLOAD_SRCS_C};${CP2K_DBM_SRCS_C};${CP2K_GRID_SRCS_C};${CP2K_PW_SRCS_C};${CP2K_DBM_SRCS_GPU_C};${CP2K_PROGS_C}"
    "*.c")
  cp2k_compare_src_with_list(
    "${CP2K_SRCS_F};${CP2K_OFFLOAD_SRCS_F};${CP2K_PROGS_F}" "*.F")
  cp2k_compare_src_with_list(
    "${CP2K_GRID_SRCS_GPU};${CP2K_GRID_SRCS_HIP};${CP2K_DBM_SRCS_GPU};${CP2K_PW_SRCS_GPU}"
    "*.cu")
  cp2k_compare_src_with_list("${CP2K_SRCS_CPP}" "*.cpp")
  return()
endif()

# ##############################################################################
# alias targets encapsulating cuda and hip libraries #
# ##############################################################################

# create some short version of these libraries to ease things afterwards
if(CP2K_USE_CUDA)
  add_library(cp2k_cuda_libs INTERFACE)
  if(CP2K_USE_NVHPC)
    target_link_libraries(cp2k_cuda_libs INTERFACE NVHPC::MATH NVHPC::CUDA)
  else()
    target_link_libraries(
      cp2k_cuda_libs
      INTERFACE $<$<BOOL:${CP2K_USE_CUSOLVER_MP}>:CUDA::cusolver> CUDA::cufft
                CUDA::cufftw CUDA::cublas CUDA::cudart CUDA::cuda_driver)
  endif()
endif()

if(CP2K_USE_HIP)
  add_library(cp2k_hip_libs INTERFACE)
  target_link_libraries(cp2k_hip_libs INTERFACE hip::hipfft roc::hipblas
                                                hip::host)
endif()

add_library(cp2k_linalg_libs INTERFACE)
target_link_libraries(
  cp2k_linalg_libs
  INTERFACE CP2K::LAPACK::lapack
            $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmf>
            $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmext>
            CP2K::BLAS::blas
            $<$<BOOL:${CP2K_USE_CUDA}>:cp2k_cuda_libs>
            $<$<BOOL:${CP2K_USE_HIP}>:cp2k_hip_libs>
            $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_Fortran>
            $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_C>
            $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_CXX>
            OpenMP::OpenMP_Fortran
            OpenMP::OpenMP_C
            OpenMP::OpenMP_CXX)

set(CP2K_GPU_DFLAGS
    $<$<BOOL:${CP2K_USE_HIP}>:__HIP_PLATFORM_AMD__
    __OFFLOAD_HIP>
    $<$<COMPILE_LANGUAGE:HIP>:__HIP_PLATFORM_AMD__
    __OFFLOAD_HIP>
    $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
    $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>)

# ##############################################################################
# set properties for the cp2k target.
# ##############################################################################

add_fypp_sources(CP2K_SRCS ${CP2K_SRCS_F})

# set the __SHORT_FILE__ per file for CP2K sources
foreach(cp2k_src ${CP2K_SRCS})
  # add_fypp_sources returns a path in the current binary dir
  get_filename_component(short_file "${cp2k_src}" NAME)
  set_source_files_properties(
    ${cp2k_src} PROPERTIES COMPILE_DEFINITIONS __SHORT_FILE__="${short_file}")
endforeach()

if(CP2K_USE_LIBTORCH)
  add_library(cp2k
              "${CP2K_SRCS};${CP2K_SRCS_C};${CP2K_SRCS_GPU};${CP2K_SRCS_CPP}")
else()
  add_library(cp2k "${CP2K_SRCS};${CP2K_SRCS_C};${CP2K_SRCS_GPU}")
endif()

target_compile_features(cp2k PUBLIC cxx_std_14)
target_compile_features(cp2k PUBLIC c_std_99)
target_compile_features(cp2k PUBLIC cuda_std_11)
target_include_directories(cp2k PUBLIC $<INSTALL_INTERFACE:include/cp2k>)

set_target_properties(
  cp2k
  PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
             Fortran_MODULE_DIRECTORY mod_files
             POSITION_INDEPENDENT_CODE ON
             VERSION ${cp2k_VERSION}
             SOVERSION ${cp2k_APIVERSION})

if(CP2K_USE_ACCEL MATCHES "CUDA")
  set_property(TARGET cp2k PROPERTY CUDA_ARCHITECTURES
                                    ${CMAKE_CUDA_ARCHITECTURES})
  target_compile_definitions(
    cp2k PUBLIC __OFFLOAD_CUDA $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>)
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
elseif(CP2K_USE_ACCEL MATCHES "HIP")
  set_property(TARGET cp2k PROPERTY HIP_ARCHITECTURES
                                    ${CMAKE_HIP_ARCHITECTURES})
endif()

# mix the target and variables pointing to the include directories.

include_directories(
  ${CMAKE_MPI_INCLUDE_DIRECTORIES}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/base>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbcsrx>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/motion>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbm>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>)

# include CPATH paths to the list of directories to include. gfortran does not
# support CPATH at all
if(DEFINED ENV{CPATH})
  string(REPLACE ":" ";" INCLUDE_LIST $ENV{CPATH})
  include_directories(${INCLUDE_LIST})
endif()

add_library(cp2k_link_libs INTERFACE)

# the order is important so keep it that way
target_link_libraries(
  cp2k_link_libs
  INTERFACE
    $<$<BOOL:${CP2K_USE_SIRIUS}>:sirius::sirius>
    $<$<BOOL:${CP2K_USE_VORI}>:CP2K::VORI::vori>
    $<$<BOOL:${CP2K_USE_PEXSI}>:CP2K::PEXSI::pexsi>
    $<$<BOOL:${CP2K_USE_PEXSI}>:CP2K::ptscotch::ptscotch>
    $<$<BOOL:${CP2K_USE_SPGLIB}>:CP2K::LIBSPG::libspg>
    $<$<BOOL:${CP2K_USE_LIBXC}>:CP2K::Libxc::xc>
    $<$<BOOL:${CP2K_USE_ELPA}>:CP2K::ELPA::elpa>
    $<$<BOOL:${CP2K_USE_FFTW3}>:CP2K::FFTW3::fftw3>
    $<$<BOOL:${CP2K_ENABLE_FFTW3_THREADS_SUPPORT}>:CP2K::FFTW3::fftw3_threads>
    $<$<BOOL:${CP2K_ENABLE_FFTW3_OPENMP_SUPPORT}>:CP2K::FFTW3::fftw3_omp>
    $<$<BOOL:${CP2K_USE_SPLA}>:SPLA::spla>
    $<$<BOOL:${CP2K_USE_LIBINT2}>:CP2K::Libint2::int2>
    $<$<BOOL:${CP2K_USE_TORCH}>:${TORCH_LIBRARIES}>
    $<$<BOOL:${CP2K_USE_COSMA}>:cosma::cosma_prefixed_pxgemm>
    $<$<BOOL:${CP2K_USE_COSMA}>:cosma::cosma>
    DBCSR::dbcsr
    $<$<BOOL:${CP2K_USE_MPI}>:CP2K::SCALAPACK::scalapack>
    cp2k_linalg_libs)

string(TIMESTAMP CP2K_TIMESTAMP "%Y-%m-%d %H:%M:%S")

target_link_libraries(cp2k PUBLIC cp2k_link_libs)
target_compile_definitions(
  cp2k
  PUBLIC $<$<BOOL:${CP2K_USE_MPI}>:__parallel>
         $<$<BOOL:${CP2K_USE_MPI}>:__SCALAPACK>
         $<$<BOOL:${CP2K_ENABLE_F08_MPI}>:__MPI_08>
         __COMPILE_DATE=\"${CP2K_TIMESTAMP}\"
         __COMPILE_HOST=\"${CP2K_HOST_NAME}\"
         __COMPILE_REVISION=\"${CP2K_GIT_HASH}\"
         __DATA_DIR=\"${CMAKE_INSTALL_FULL_DATAROOTDIR}/cp2k/data\"
         __COMPILE_ARCH=\"${CMAKE_SYSTEM_PROCESSOR}\"
         $<$<CONFIG:Release>:NDEBUG>
         $<$<CONFIG:COVERAGE>:__NO_ABORT>
         $<$<CONFIG:DEBUG>:__HAS_IEEE_EXCEPTIONS>
         $<$<CONFIG:DEBUG>:__CHECK_DIAG>
         $<$<BOOL:${CP2K_USE_PEXSI}>:__PEXSI>
         $<$<BOOL:${CP2K_USE_PLUMED2}>:__PLUMED2>
         $<$<BOOL:${CP2K_USE_QUIP}>:__QUIP>
         $<$<BOOL:${CP2K_USE_MAXWELL}>:__LIBMAXWELL>
         $<$<BOOL:${CP2K_USE_VORI}>:__LIBVORI>
         $<$<BOOL:${CP2K_USE_SPGLIB}>:__SPGLIB>
         $<$<BOOL:${CP2K_USE_SIRIUS}>:__SIRIUS>
         $<$<BOOL:${CP2K_USE_SPLA}>:__SPLA>
         $<$<BOOL:${CP2K_USE_SpFFT}>:__SPFFT>
         $<$<BOOL:${CP2K_USE_SPLA_GEMM_OFFLOADING}>:__OFFLOAD_GEMM>
         $<$<BOOL:${CP2K_USE_ELPA}>:__ELPA>
         $<$<BOOL:${CP2K_USE_LIBXC}>:__LIBXC>
         $<$<BOOL:${CP2K_USE_FFTW3}>:__FFTW3>
         $<$<BOOL:${CP2K_USE_LIBINT2}>:__LIBINT>
         $<$<BOOL:${CP2K_USE_PEXSI}>:__LIBPEXSI>
         $<$<BOOL:${CP2K_USE_LIBTORCH}>:__LIBTORCH>
         $<$<BOOL:${CP2K_USE_COSMA}>:__COSMA>
         $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
         $<$<STREQUAL:${CP2K_BLAS_VENDOR},MKL>:__MKL>
         $<$<STREQUAL:${CP2K_BLAS_VENDOR},Apple>:__ACCELERATE>
         $<$<BOOL:${CP2K_USE_CUSOLVER_MP}>:__CUSOLVERMP>
         $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
         $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>
         $<$<BOOL:${CP2K_USE_HIP}>:__HIP_PLATFORM_AMD__
         __OFFLOAD_HIP>
         $<$<COMPILE_LANGUAGE:HIP>:__HIP_PLATFORM_AMD__
         __OFFLOAD_HIP>)

if(CP2K_USE_CUDA OR CP2K_USE_HIP)
  # these checks are only relevant when the debug mode is on

  target_compile_definitions(
    cp2k
    PUBLIC $<$<NOT:$<BOOL:${CP2K_DBCSR_USE_CPU_ONLY}>>:__DBCSR_ACC>
           $<$<BOOL:${CP2K_DISABLE_GRID_GPU}>:__NO_OFFLOAD_GRID>
           $<$<BOOL:${CP2K_DISABLE_PW_GPU}>:__NO_OFFLOAD_PW>
           $<$<BOOL:${CP2K_DISABLE_DBM_GPU}>:__NO_OFFLOAD_DBM>)
endif()

# ##############################################################################
# applications
# ##############################################################################

list(
  APPEND
  __CP2K_APPS
  "memory_utilities_unittest"
  "parallel_rng_types_unittest"
  "graph"
  "dumpdcd"
  "xyz2dcd"
  "libcp2k_unittest"
  "nequip_unittest"
  "dbt_unittest"
  "dbt_tas_unittest")

add_executable(cp2k-bin start/cp2k.F)
# for linking

add_executable(memory_utilities_unittest common/memory_utilities_unittest.F)
add_executable(parallel_rng_types_unittest common/parallel_rng_types_unittest.F)
add_executable(graph metadyn_tools/graph.F)
add_executable(dumpdcd motion/dumpdcd.F)
add_executable(xyz2dcd motion/xyz2dcd.F)
add_executable(libcp2k_unittest start/libcp2k_unittest.c)
add_executable(nequip_unittest nequip_unittest.F)
add_executable(dbt_unittest dbt/dbt_unittest.F)
add_executable(dbt_tas_unittest dbt/tas/dbt_tas_unittest.F)

set_target_properties(
  cp2k-bin
  PROPERTIES CXX_STANDARD 14
             IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
             POSITION_INDEPENDENT_CODE ON
             LINKER_LANGUAGE "CXX"
             OUTPUT_NAME "cp2k.${__cp2k_ext}")

target_link_libraries(cp2k-bin PUBLIC cp2k)

foreach(_app ${__CP2K_APPS})
  set_target_properties(
    ${_app}
    PROPERTIES CXX_STANDARD 14
               IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
               POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE "CXX"
               OUTPUT_NAME "${_app}.${__cp2k_ext}")

  target_compile_features(${_app} PUBLIC c_std_99)
  if(CP2K_USE_ACCEL MATCHES HIP)
    set_property(TARGET ${_app} PROPERTY HIP_ARCHITECTURES
                                         ${CMAKE_HIP_ARCHITECTURES})
  endif()

  target_link_libraries(${_app} PUBLIC cp2k)
endforeach()

# ##############################################################################
# different miniapps
# ##############################################################################

add_executable(dbm_miniapp dbm/dbm_miniapp.c ${CP2K_DBM_SRCS_C}
                           ${CP2K_DBM_SRCS_GPU} ${CP2K_OFFLOAD_SRCS_C})
target_compile_definitions(
  dbm_miniapp
  PRIVATE $<$<CONFIG:Release>:NDEBUG> $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
          $<$<STREQUAL:${CP2K_BLAS_VENDOR},"MKL">:__MKL>)

if(NOT CP2K_DISABLE_DBM_GPU)
  target_compile_definitions(dbm_miniapp PRIVATE __NO_OFFLOAD_DBM)
else()
  target_compile_definitions(dbm_miniapp PRIVATE ${CP2K_GPU_DFLAGS})
endif()

foreach(__app grid_miniapp grid_unittest)
  add_executable(
    ${__app}
    "grid/${__app}.c"
    ${CP2K_GRID_SRCS_C}
    $<$<OR:$<BOOL:${CP2K_USE_CUDA}>,$<BOOL:${CP2K_USE_HIP}>>:${CP2K_GRID_SRCS_GPU}>
    $<$<BOOL:${CP2K_USE_HIP}>:${CP2K_GRID_SRCS_HIP}>
    ${CP2K_OFFLOAD_SRCS_C})
  target_compile_definitions(
    ${__app}
    PRIVATE $<$<CONFIG:Release>:NDEBUG>
            $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
            $<$<STREQUAL:${CP2K_BLAS_VENDOR},"MKL">:__MKL>)

  if(CP2K_DISABLE_GRID_GPU)
    target_compile_definitions(${__app} PUBLIC __NO_OFFLOAD_GRID)
  else()
    target_compile_definitions(${__app} PUBLIC ${CP2K_GPU_DFLAGS})
  endif()
endforeach()

# ##############################################################################
# common to all c executables
# ##############################################################################

# list of all targets that share the same properties. simplify code as most
# cmake function can work on list of targets

foreach(__app grid_miniapp grid_unittest dbm_miniapp)
  set_target_properties(
    ${__app}
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE C
               POSITION_INDEPENDENT_CODE ON
               OUTPUT_NAME "${__app}.${__cp2k_ext}")

  if(CP2K_USE_HIP)
    set_property(TARGET ${__app} PROPERTY HIP_ARCHITECTURES
                                          ${CMAKE_HIP_ARCHITECTURES})
  endif()

  if(CP2K_USE_CUDA)
    set_target_properties(${__app} PROPERTIES CUDA_ARCHITECTURES
                                              ${CMAKE_CUDA_ARCHITECTURES})
  endif()

  target_compile_features(${__app} PUBLIC cuda_std_14 c_std_99 cxx_std_14)

  target_include_directories(
    ${__app}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/base>
           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/common>
           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
           $<$<BOOL:${CP2K_USE_CUDA}>:${CMAKE_CUDA_INCLUDE_DIRECTORIES}>
           $<$<BOOL:${CP2K_USE_HIP}>:${CMAKE_HIP_INCLUDE_DIRECTORIES}>)

  target_link_libraries(${__app} PUBLIC cp2k_linalg_libs m)
endforeach()

# ##############################################################################
# installation
# ##############################################################################
install(
  TARGETS ${__CP2K_APPS} cp2k-bin grid_miniapp grid_unittest dbm_miniapp
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  INCLUDES
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k")

install(DIRECTORY "${CMAKE_SOURCE_DIR}/data"
        DESTINATION "${CMAKE_INSTALL_FULL_DATAROOTDIR}/cp2k")

# ##############################################################################
# cp2kConfig.cmake, etc...
# ##############################################################################

set(cp2k_export_libs cp2k cp2k_link_libs cp2k_linalg_libs)
if(CP2K_USE_CUDA)
  list(APPEND cp2k_export_libs cp2k_cuda_libs)
endif()

if(CP2K_USE_HIP)
  list(APPEND cp2k_export_libs cp2k_hip_libs)
endif()

install(
  TARGETS ${cp2k_export_libs}
  EXPORT cp2k_targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(
  EXPORT cp2k_targets
  FILE cp2kTargets.cmake
  NAMESPACE cp2k::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/cp2k")

install(FILES start/libcp2k.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k")

install(
  DIRECTORY "${PROJECT_BINARY_DIR}/src/mod_files"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k"
  FILES_MATCHING
  PATTERN "*.mod")
