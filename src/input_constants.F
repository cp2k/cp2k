!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collects all constants needed in input so that they can be used without
!>      circular dependencies
!> \par History
!>      02.2005 created [fawzi&Joost]
!>      10.2008 Teodoro Laino [tlaino] - University of Zurich
!>              Cleaned this file of all arrays of strings. Here must be kept
!>              only integer/real global constants
! *****************************************************************************
MODULE input_constants
  USE f77_blas
  USE kinds,                           ONLY: dp

  IMPLICIT NONE
  PRIVATE

  ! Tags for development/release versions
  INTEGER, PARAMETER, PUBLIC               :: id_development_version=0,&
                                              id_release_version=1

  ! Program Names
  INTEGER, PARAMETER, PUBLIC               :: do_atom=1,&
                                              do_farming=2,&
                                              do_test=3,&
                                              do_cp2k=4,&
                                              do_optimize_input=5,&
                                              do_opt_basis=6,&
                                              do_tree_mc=7,&
                                              do_tree_mc_ana=8,&
                                              do_swarm=9

  INTEGER, PARAMETER, PUBLIC               :: opt_force_matching=1

  INTEGER, PARAMETER, PUBLIC               :: do_full_density=0,&
                                              do_spin_density=1

  INTEGER, PARAMETER, PUBLIC               :: do_no_et=0,&
                                              do_et_ddapc=1,&
                                              do_et_becke=2
  ! Method Names
  INTEGER, PARAMETER, PUBLIC               :: do_qs=1,&
                                              do_fist=2,&
                                              do_qmmm=4,&
                                              do_eip=5,&
                                              do_ep=6,&
                                              do_mixed=7

  ! QMMM
  REAL(KIND=dp), PARAMETER, PUBLIC         :: RADIUS_QMMM_DEFAULT = 0.80_dp,& ! Angstrom
                                              ALPHA_IMOMM_DEFAULT = 1.38_dp,&
                                              CHARGE_SCALE_FACTOR = 1.00_dp
  ! BAND
  INTEGER, PARAMETER, PUBLIC               :: band_md_opt=0,&
                                              band_diis_opt=1

  INTEGER, PARAMETER, PUBLIC               :: do_band_cartesian=0,&
                                              do_band_collective=1
  ! Weight_types
  INTEGER, PARAMETER, PUBLIC               :: weight_type_unit=0,&
                                              weight_type_mass=1

  ! diag methods
  INTEGER, PARAMETER, PUBLIC               :: do_diag_syevd = 1,&
                                              do_diag_syevx = 2
  ! initial matrix methods methods
  INTEGER, PARAMETER, PUBLIC               :: do_mat_random = 1,&
                                              do_mat_read = 2

  INTEGER, PARAMETER, PUBLIC               :: do_constr_none=1,&
                                              do_constr_atomic=2,&
                                              do_constr_molec=3

  INTEGER, PARAMETER, PUBLIC               :: do_cell_cp2k=1,&
                                              do_cell_xsc=2

  INTEGER, PARAMETER, PUBLIC               :: rmsd_all=0,&
                                              rmsd_list=1,&
                                              rmsd_weightlist=2

  INTEGER, PARAMETER, PUBLIC               :: do_rep_interleaved=0,&
                                              do_rep_blocked=1

  ! MD Ensemble
  INTEGER, PARAMETER, PUBLIC               :: nve_ensemble=1,&
                                              nvt_ensemble=2,&
                                              npt_i_ensemble=3,&
                                              npt_f_ensemble=4,&
                                              nph_uniaxial_ensemble=5,&
                                              nph_uniaxial_damped_ensemble=6,&
                                              nph_ensemble=7,&
                                              isokin_ensemble=8,&
                                              reftraj_ensemble=9,&
                                              langevin_ensemble=10,&
                                              npe_f_ensemble=11,&
                                              npe_i_ensemble=12,&
                                              nvt_adiabatic_ensemble=13

  INTEGER, PARAMETER, PUBLIC               :: do_thermo_same_as_part=100,&
                                              do_thermo_nose=101,&
                                              do_thermo_csvr=102,&
                                              do_thermo_gle=103, &
                                              do_thermo_al=104


  ! Run Types
  INTEGER, PARAMETER, PUBLIC               :: none_run=1,&
                                              energy_run=2,&
                                              energy_force_run=3,&
                                              mol_dyn_run=4,&
                                              geo_opt_run=5,&
                                              mon_car_run=6,&
                                              electronic_spectra_run=7,&
                                              debug_run=8,&
                                              bsse_run=9,&
                                              linear_response_run=10,&
                                              pint_run=11,&
                                              vib_anal=12,&
                                              do_band=13,&
                                              cell_opt_run=14,&
                                              real_time_propagation=15,&
                                              ehrenfest=16, &
                                              do_tamc=17, &
                                              tree_mc_run=18,&
                                              driver_run=19

  ! Run Types of Atom Code
  INTEGER, PARAMETER, PUBLIC               :: atom_no_run=1,&
                                              atom_energy_run=2,&
                                              atom_basis_run=3,&
                                              atom_pseudo_run=4

  ! Integral types in Atom Code
  INTEGER, PARAMETER, PUBLIC               :: do_analytic=1,&
                                              do_semi_analytic=2,&
                                              do_numeric=3
  ! Basis set type for Atom Code
  INTEGER, PARAMETER, PUBLIC               :: gaussian=1,&
                                              geometrical_gto=2,&
                                              contracted_gto=3,&
                                              slater=4,&
                                              numerical=5
  ! Method type for Atom Code
  INTEGER, PARAMETER, PUBLIC               :: do_rks_atom=1,&
                                              do_uks_atom=2,&
                                              do_rhf_atom=3,&
                                              do_uhf_atom=4,&
                                              do_rohf_atom=5
  ! Relativistic correction for Atom Code
  INTEGER, PARAMETER, PUBLIC               :: do_nonrel_atom=100,&
                                              do_zoramp_atom=101,&
                                              do_sczoramp_atom=102,&
                                              do_dkh0_atom=110,&
                                              do_dkh1_atom=111,&
                                              do_dkh2_atom=112,&
                                              do_dkh3_atom=113,&
                                              do_dkh4_atom=114,&
                                              do_dkh5_atom=115

  ! Pseudopotential type for Atom Code
  INTEGER, PARAMETER, PUBLIC               :: no_pseudo=0,&
                                              gth_pseudo=1
  ! BAND
  INTEGER, PARAMETER, PUBLIC               :: do_b_neb=0,&
                                              do_it_neb=1,&
                                              do_ci_neb=2,&
                                              do_d_neb=3,&
                                              do_sm=4,&
                                              do_eb=5

  INTEGER, PARAMETER, PUBLIC               :: pot_neb_full=0,&
                                              pot_neb_fe=1,&
                                              pot_neb_me=2

  INTEGER, PARAMETER, PUBLIC               :: cholesky_off=0,&
                                              cholesky_reduce=1,&
                                              cholesky_restore=2,&
                                              cholesky_inverse=3,&
                                              cholesky_dbcsr=4

  INTEGER, PARAMETER, PUBLIC               :: atomic_guess=1,&
                                              restart_guess=2,&
                                              random_guess=3,&
                                              core_guess=4,&
                                              densities_guess=5,&
                                              history_guess=6,&
                                              no_guess=7,&
                                              mopac_guess=8,&
                                              sparse_guess=9,&
                                              almo_scf_guess=10

  INTEGER, PARAMETER, PUBLIC               :: pw_interp=1,&
                                              spline3_nopbc_interp=2,&
                                              spline3_pbc_interp=3

  INTEGER, PARAMETER, PUBLIC               :: do_gapw_gcs=1,&
                                              do_gapw_gct=2,&
                                              do_gapw_log=3

  INTEGER, PARAMETER, PUBLIC               :: do_pwgrid_spherical=1,&
                                              do_pwgrid_ns_fullspace=2,&
                                              do_pwgrid_ns_halfspace=3

  INTEGER, PARAMETER, PUBLIC               :: wfi_use_guess_method_nr=0,&
                                              wfi_use_prev_p_method_nr=1,&
                                              wfi_use_prev_rho_r_method_nr=2,&
                                              wfi_linear_wf_method_nr=3,&
                                              wfi_linear_p_method_nr=4,&
                                              wfi_linear_ps_method_nr=5,&
                                              wfi_use_prev_wf_method_nr=6,&
                                              wfi_ps_method_nr=7,&
                                              wfi_frozen_method_nr=8,&
                                              wfi_aspc_nr=9

  INTEGER, PARAMETER, PUBLIC               :: do_method_undef=0,&
                                              do_method_gapw=1,&
                                              do_method_gapw_xc=2,&
                                              do_method_gpw=3,&
                                              do_method_lrigpw=4,&
                                              do_method_mndo=7,&
                                              do_method_am1=8,&
                                              do_method_pm3=9,&
                                              do_method_pm6=10,&
                                              do_method_mndod=11,&
                                              do_method_pdg=12,&
                                              do_method_rm1=13,&
                                              do_method_dftb=14,&
                                              do_method_pchg=15,&
                                              do_method_pnnl=16,&
                                              do_method_ofgpw=17,&
                                              do_method_scptb=18

  INTEGER, PARAMETER, PUBLIC               :: do_multipole_section_on=0,&
                                              do_multipole_section_off=1

  INTEGER, PARAMETER, PUBLIC               :: do_ppl_analytic=0,&
                                              do_ppl_grid=1

  INTEGER, PARAMETER, PUBLIC               :: do_se_IS_kdso=0,&
                                              do_se_IS_kdso_d=1,&
                                              do_se_IS_slater=2

  INTEGER, PARAMETER, PUBLIC               :: do_se_lr_none=0,&
                                              do_se_lr_ewald=1,&
                                              do_se_lr_ewald_r3=2,&
                                              do_se_lr_ewald_gks=3

  INTEGER, PARAMETER, PUBLIC               :: do_par_atom=0,&
                                              do_par_grid=1

  INTEGER, PARAMETER, PUBLIC               :: do_qmmm_center_every_step=0,&
                                              do_qmmm_center_setup_only=1,&
                                              do_qmmm_center_never=2

  INTEGER, PARAMETER, PUBLIC                    :: do_qmmm_center_max_minus_min=0,&
                                              do_qmmm_center_pbc_aware=1
  ! Thermostats
  INTEGER, PARAMETER, PUBLIC               :: do_thermo_no_communication=0,&
                                              do_thermo_communication=1,&
                                              do_thermo_only_master=2

  ! pint coordinate transformations
  INTEGER, PARAMETER, PUBLIC               :: transformation_normal=1,&
                                              transformation_stage=2

  ! helium pbc unit cell shapes
  INTEGER, PARAMETER, PUBLIC               :: helium_cell_shape_none=0,&
                                              helium_cell_shape_cube=1,&
                                              helium_cell_shape_octahedron=2

  ! Free Energy methods
  INTEGER, PARAMETER, PUBLIC               :: do_fe_meta=0,&
                                              do_fe_ui=1,&
                                              do_fe_ac=2

  INTEGER, PARAMETER, PUBLIC               :: do_wall_none=0,&
                                              do_wall_reflective=1,&
                                              do_wall_quadratic=2,&
                                              do_wall_quartic=3,&
                                              do_wall_gaussian=4

  INTEGER, PARAMETER, PUBLIC               :: do_wall_p=0,&
                                              do_wall_m=1
  ! Output formats
  INTEGER, PARAMETER, PUBLIC               :: dump_atomic=1,&
                                              dump_xmol=2,&
                                              dump_dcd=3,&
                                              dump_pdb=4,&
                                              dump_dcd_aligned_cell=5

  INTEGER, PARAMETER, PUBLIC               :: fmt_id_xyz=1,&
                                              fmt_id_pdb=2

  INTEGER, PARAMETER, PUBLIC               :: do_qmmm_none=0,&
                                              do_qmmm_coulomb=1,&
                                              do_qmmm_gauss=2,&
                                              do_qmmm_swave=3,&
                                              do_qmmm_link_imomm=1,&
                                              do_qmmm_link_gho=2,&
                                              do_qmmm_link_pseudo=3,&
                                              do_qmmm_center_none=0,&
                                              do_qmmm_center_grid=1,&
                                              do_qmmm_wall_none=0,&
                                              do_qmmm_wall_reflective=1,&
                                              do_qmmm_wall_quadratic=2,&
                                              do_qmmm_image_calcmatrix=0,&
                                              do_qmmm_image_iter=1

  INTEGER, PARAMETER, PUBLIC               :: do_fm_mom_conserv_none=0, &
                                              do_fm_mom_conserv_equal_f=1, &
                                              do_fm_mom_conserv_equal_a=2

  INTEGER, PARAMETER, PUBLIC               :: do_fm_mom_conserv_core=0, &
                                              do_fm_mom_conserv_qm=1, &
                                              do_fm_mom_conserv_buffer=2

  INTEGER, PARAMETER, PUBLIC               :: do_region_global=0, &
                                              do_region_molecule=1,&
                                              do_region_massive=2,&
                                              do_region_none=3,&
                                              do_region_defined=4,&
                                              do_bondparm_covalent=0,&
                                              do_bondparm_vdw=1,&
                                              do_skip_11=0,&
                                              do_skip_12=1,&
                                              do_skip_13=2,&
                                              do_skip_14=3,&
                                              do_conn_psf=1,&
                                              do_conn_psf_u=2,&
                                              do_conn_generate=3,&
                                              do_conn_off=4,&
                                              do_conn_mol_set=5,&
                                              do_conn_g87=6,&
                                              do_conn_g96=7,&
                                              do_conn_amb7=8,&
                                              do_conn_user=9,&
                                              do_coord_off=0,&
                                              do_coord_pdb=1,&
                                              do_coord_xyz=2,&
                                              do_coord_g96=3,&
                                              do_coord_crd=4,&
                                              do_coord_cif=5,&
                                              do_coord_xtl=6,&
                                              do_coord_cp2k=7

  INTEGER, PARAMETER, PUBLIC               :: do_stress_none=0,&
                                              do_stress_analytical=1,&
                                              do_stress_numerical=2,&
                                              do_stress_diagonal_anal=3,&
                                              do_stress_diagonal_numer=4

  INTEGER, PARAMETER, PUBLIC               :: do_mc_traditional=0,&
                                              do_mc_gemc_nvt=1,&
                                              do_mc_gemc_npt=2,&
                                              do_mc_virial=3

  INTEGER, PARAMETER, PUBLIC               :: do_loc_none=0,&
                                              do_loc_jacobi=1,&
                                              do_loc_crazy=2,&
                                              do_loc_direct=3,&
                                              do_loc_l1_norm_sd=4

  INTEGER, PARAMETER, PUBLIC               :: do_loc_min=0,&
                                              do_loc_max=1,&
                                              op_loc_berry=0,&
                                              op_loc_boys=1,&
                                              op_loc_pipek=2,&
                                              op_loc_l1_norm=3,&
                                              state_loc_all=0,&
                                              state_loc_range=1,&
                                              state_loc_list=2,&
                                              energy_loc_range=3,&
                                              state_loc_none=4

  INTEGER, PARAMETER, PUBLIC               :: do_loc_homo=0,&
                                              do_loc_lumo=1,&
                                              do_loc_both=2

  INTEGER, PARAMETER, PUBLIC               :: orb_s=0,&
                                              orb_px=1,&
                                              orb_py=2,&
                                              orb_pz=3,&
                                              orb_dxy=4,&
                                              orb_dyz=5,&
                                              orb_dzx=6,&
                                              orb_dx2=7,&
                                              orb_dy2=8,&
                                              orb_dz2=9

  INTEGER, PARAMETER, PUBLIC               :: use_mom_ref_com=0,&
                                              use_mom_ref_coac=1,&
                                              use_mom_ref_user=2,&
                                              use_mom_ref_zero=3

  INTEGER, PARAMETER, PUBLIC               :: current_gauge_r=1,&
                                              current_gauge_r_and_step_func=2,&
                                              current_gauge_atom=3,&
                                              current_orb_center_wannier=1,&
                                              current_orb_center_common=2,&
                                              current_orb_center_atom=3,&
                                              current_orb_center_box=4

  INTEGER, PARAMETER, PUBLIC               :: ot_mini_sd=1,ot_mini_cg=2,ot_mini_diis=3,ot_mini_broyden=4
  INTEGER, PARAMETER, PUBLIC               :: ot_algo_taylor_or_diag=1,ot_algo_irac=2
  INTEGER, PARAMETER, PUBLIC               :: ot_chol_irac=1,ot_poly_irac=2,ot_lwdn_irac=3
  INTEGER, PARAMETER, PUBLIC               :: ls_none=1,ls_2pnt=2,ls_3pnt=3,ls_gold=4,ls_fit=5
  INTEGER, PARAMETER, PUBLIC               :: ot_precond_none=0, &
                                              ot_precond_full_single=1, &
                                              ot_precond_full_kinetic=2, &
                                              ot_precond_s_inverse=3, &
                                              ot_precond_full_all=5,&
                                              ot_precond_full_single_inverse=6
  INTEGER, PARAMETER, PUBLIC               :: ot_precond_solver_default=1,&
                                              ot_precond_solver_inv_chol=2,&
                                              ot_precond_solver_direct=3

  INTEGER, PARAMETER, PUBLIC               :: diag_standard=1,&
                                              diag_ot=2, diag_block_krylov=3,&
                                              diag_block_davidson=4

  INTEGER, PUBLIC, PARAMETER               :: default_lbfgs_method_id=1,&
                                              default_bfgs_method_id=2,&
                                              default_cg_method_id=3

  INTEGER, PUBLIC, PARAMETER               :: default_minimization_method_id=1,&
                                              default_ts_method_id=2,&
                                              default_cell_method_id=3,&
                                              default_shellcore_method_id=4

  INTEGER, PUBLIC, PARAMETER               :: none_ts_method_id=0,&
                                              default_dimer_method_id=1

  INTEGER, PUBLIC, PARAMETER               :: do_first_rotation_step=0,&
                                              do_second_rotation_step=1,&
                                              do_third_rotation_step=2

  INTEGER, PUBLIC, PARAMETER               :: default_cell_geo_opt_id=0,&
                                              default_cell_md_id=1,&
                                              default_cell_direct_id=2

  INTEGER, PARAMETER, PUBLIC               :: add_last_no=0,&
                                              add_last_numeric=1,&
                                              add_last_symbolic=2
  INTEGER, PARAMETER, PUBLIC               :: silent_print_level=0,&
                                              low_print_level=1,&
                                              medium_print_level=2,&
                                              high_print_level=3,&
                                              debug_print_level=4

  INTEGER, PARAMETER, PUBLIC               :: xc_funct_no_shortcut=0,&
                                              xc_funct_blyp=1,&
                                              xc_funct_pade=2,&
                                              xc_funct_pbe=3,&
                                              xc_funct_tpss=4,&
                                              xc_funct_hcth120=5,&
                                              xc_funct_olyp=6,&
                                              xc_funct_bp=7,&
                                              xc_none=8,&
                                              xc_funct_xwpbe=9,&
                                              xc_funct_b3lyp=10,&
                                              xc_funct_pbe0=11,&
                                              xc_funct_beefvdw=12
  INTEGER, PARAMETER, PUBLIC               :: sic_none=0,&
                                              sic_mauri_us=1,&
                                              sic_mauri_spz=2,&
                                              sic_ad=3,&
                                              sic_eo=4
  INTEGER, PARAMETER, PUBLIC               :: sic_list_all=1,&
                                              sic_list_unpaired=2
  INTEGER, PARAMETER, PUBLIC               :: tddfpt_singlet=0,&
                                              tddfpt_triplet=1
  INTEGER, PARAMETER, PUBLIC               :: tddfpt_lanczos=0,&
                                              tddfpt_davidson=1
  INTEGER, PARAMETER, PUBLIC               :: oe_none=0,&
                                              oe_lb=1,&
                                              oe_gllb=2,&
                                              oe_saop=3,&
                                              oe_sic=4
  INTEGER, PARAMETER, PUBLIC               :: no_excitations=0,&
                                              tddfpt_excitations=1
  INTEGER, PARAMETER, PUBLIC               :: xc_vdw_fun_none=100,&
                                              xc_vdw_fun_pairpot=101,&
                                              xc_vdw_fun_nonloc=102
  INTEGER, PARAMETER, PUBLIC               :: vdw_pairpot_dftd2=1,&
                                              vdw_pairpot_dftd3=2,&
                                              vdw_pairpot_dftd3bj=3
  INTEGER, PARAMETER, PUBLIC               :: vdw_nl_DRSLL=1,&
                                              vdw_nl_LMKLL=2,&
                                              vdw_nl_RVV10=3

  ! Input constants for relativistic calculations
  INTEGER, PARAMETER, PUBLIC               :: rel_none=0,&
                                              rel_dkh=1,&
                                              rel_zora=2
  INTEGER, PARAMETER, PUBLIC               :: rel_zora_full=100,&
                                              rel_zora_mp=101,&
                                              rel_sczora_mp=102
  INTEGER, PARAMETER, PUBLIC               :: rel_trans_full=0,&
                                              rel_trans_molecule=1,&
                                              rel_trans_atom=2
  INTEGER, PARAMETER, PUBLIC               :: rel_pot_full=0,&
                                              rel_pot_erfc=1

  INTEGER, PARAMETER, PUBLIC               :: do_add=1,&
                                              do_remove=2


  INTEGER, PARAMETER, PUBLIC               :: use_perd_x=0,&
                                              use_perd_y=1,&
                                              use_perd_z=2,&
                                              use_perd_xy=3,&
                                              use_perd_xz=4,&
                                              use_perd_yz=5,&
                                              use_perd_xyz=6,&
                                              use_perd_none=7

  INTEGER, PARAMETER, PUBLIC               :: xc_pot_energy_none=0,&
                                              xc_pot_energy_xc_functional=1,&
                                              xc_pot_energy_sum_eigenvalues=2

  INTEGER, PARAMETER, PUBLIC               :: xc_pot_none = 1000,&
                                              xc_pot_saop = 1001

  ! X-Ray spectra Calculation
  INTEGER, PARAMETER, PUBLIC               :: xas_none=0,&
                                              xas_tp_hh=1,&
                                              xas_tp_fh=2,&
                                              xes_tp_val=3,&
                                              xas_tp_xhh=4,&
                                              xas_tp_xfh=5,&
                                              xas_dscf=6
  INTEGER, PARAMETER, PUBLIC               :: xas_1s_type=1,&
                                              xas_2s_type=2,&
                                              xas_2p_type=3,&
                                              xas_3s_type=4,&
                                              xas_3p_type=5,&
                                              xas_3d_type=6,&
                                              xas_4s_type=7,&
                                              xas_4p_type=8,&
                                              xas_4d_type=9,&
                                              xas_4f_type=10
  INTEGER, PARAMETER, PUBLIC               :: xas_dip_len=1,&
                                              xas_dip_vel=2,&
                                              xas_dip_len2=3
  INTEGER, PARAMETER, PUBLIC               :: xas_scf_default=0,&
                                              xas_scf_general=1

  ! Linear Response for properties
  INTEGER, PARAMETER, PUBLIC               :: lr_none=0,&
                                              lr_chemshift=1,&
                                              lr_gtensor=2,&
                                              lr_current=3

  ! EIP models
  INTEGER, PARAMETER, PUBLIC               :: use_bazant_eip=1, &
                                              use_lenosky_eip=2

  ! ddapc restraint forms
  INTEGER, PARAMETER, PUBLIC               :: do_ddapc_restraint=773,&
                                              do_ddapc_constraint=774

  ! outer scf types
  INTEGER, PARAMETER, PUBLIC               :: outer_scf_ddapc_constraint=123, &
                                              outer_scf_s2_constraint=124,&
                                              outer_scf_becke_constraint=125, &
                                              outer_scf_scp=126, &
                                              outer_scf_none=127

  ! outer scf optimizers
  INTEGER, PARAMETER, PUBLIC               :: outer_scf_optimizer_sd=1001, &
                                              outer_scf_optimizer_diis=1002, &
                                              outer_scf_optimizer_none=1003, &
                                              outer_scf_optimizer_bisect=1004

  ! s2 restraint forms
  INTEGER, PARAMETER, PUBLIC               :: do_s2_restraint=872,&
                                              do_s2_constraint=873

  ! ROKS schemes
  INTEGER, PARAMETER, PUBLIC               :: general_roks   = 1,&
                                              high_spin_roks = 2

  ! mixing of force envs
  INTEGER, PARAMETER, PUBLIC               :: mix_linear_combination=701, &
                                              mix_minimum=702, &
                                              mix_coupled=703, &
                                              mix_restrained=704,&
                                              mix_generic=705,&
                                              mix_mymix=706

  ! cost models for distribution 2d
  INTEGER, PARAMETER, PUBLIC               :: model_block_count=453, &
                                              model_block_surface=454, &
                                              model_block_lmax=455
  ! optimization targets for distribution 2d
  INTEGER, PARAMETER, PUBLIC               :: opt_target_minimum=456, &
                                              opt_target_median=457

  ! clustering methods
  INTEGER, PARAMETER, PUBLIC               :: clustering_none=0, &
                                              clustering_box=1, &
                                              clustering_h_on_heavier=2, &
                                              clustering_small_on_larger=3

  ! HFX potential type
  INTEGER, PARAMETER, PUBLIC               :: do_hfx_potential_coulomb=1,&
                                              do_hfx_potential_short=2,&
                                              do_hfx_potential_long=3,&
                                              do_hfx_potential_mix_cl=4,&
                                              do_hfx_potential_gaussian=5,&
                                              do_hfx_potential_mix_lg=6,&
                                              do_hfx_potential_id=7,&
                                              do_hfx_potential_truncated=8,&
                                              do_hfx_potential_mix_cl_trunc=9,&
                                              do_mp2_potential_TShPSC=1009

  ! HFX periodic number of shells
  INTEGER, PARAMETER, PUBLIC               :: do_hfx_auto_shells = -1

  ! adiabatic hybrid type
  INTEGER, PARAMETER, PUBLIC               :: do_adiabatic_hybrid_mcy3=1

  ! adiabatic model type
  INTEGER, PARAMETER, PUBLIC               :: do_adiabatic_model_pade=1

  ! HFX evaluation type
  INTEGER, PARAMETER, PUBLIC               :: hfx_do_eval_energy=1,&
                                              hfx_do_eval_forces=2

  ! mode selctive vibrational analysis
  INTEGER, PARAMETER, PUBLIC               :: ms_guess_bfgs=1, &
                                              ms_guess_atomic=2, &
                                              ms_guess_restart=3,&
                                              ms_guess_restart_vec=4,&
                                              ms_guess_molden=5

  ! BASIS SET PARAMETER
  ! This parameter identifies basis sets
  INTEGER, PARAMETER, PUBLIC               :: use_orb_basis_set=1,&
                                              use_aux_fit_basis_set=2,&
                                              use_aux_basis_set=3,&
                                              use_scp_basis_set=4,&
                                              use_ri_aux_basis_set=5,&
                                              use_lri_basis_set=6
  ! ADMM PARAMETER
  ! This parameter defines the method to use for the purification
  INTEGER, PARAMETER, PUBLIC               :: do_admm_purify_none=1,&
                                              do_admm_purify_cauchy=2,&
                                              do_admm_purify_cauchy_subspace=3,&
                                              do_admm_purify_mo_diag=4,&
                                              do_admm_purify_mo_no_diag=5,&
                                              do_admm_purify_mcweeny=6
  ! This parameter defines when to use blocking
  INTEGER, PARAMETER, PUBLIC               :: do_admm_basis_projection=0,&
                                              do_admm_blocking_purify_full=1,&
                                              do_admm_blocked_projection=2

  ! Define DFT+U method identifier codes
  INTEGER, PARAMETER, PUBLIC               :: plus_u_lowdin=1,&
                                              plus_u_mulliken=2,&
                                              plus_u_mulliken_charges=3

  ! Define identifier codes for smearing method
  INTEGER, PARAMETER, PUBLIC               :: smear_energy_window=1,&
                                              smear_fermi_dirac=2,&
                                              smear_list=3

  ! Define identifier codes for mixing method

  INTEGER, PARAMETER, PUBLIC               :: no_mix =0,&
                                              direct_p_mix=1,&
                                              kerker_mix =2,&
                                              pulay_mix=3,&
                                              broy_mix=4,&
                                              broy_mix_new=5,&
                                              multisec_mix=6

  ! DIAGONALIZATION library 
  INTEGER, PARAMETER, PUBLIC               :: do_diag_sl=1,&
                                              do_diag_sl2=2,&
                                              do_diag_elpa=3

  ! FFT library
  ! these might need sync with fft_lib.F
  INTEGER, PARAMETER, PUBLIC               :: do_fft_sg=1,&
                                              do_fft_fftw3=3
  INTEGER, PARAMETER, PUBLIC               :: fftw_plan_estimate=1, &
                                              fftw_plan_measure=2, &
                                              fftw_plan_patient=3, &
                                              fftw_plan_exhaustive=4

  ! SCP-DFT
  INTEGER, PARAMETER, PUBLIC               :: copy_coeff=1,  &
                                              copy_fcoeff=2
  !real time propagation
  INTEGER, PARAMETER, PUBLIC               :: do_taylor=1,  &
                                              do_pade=2,&
                                              do_arnoldi=3,&
                                              do_bch=4

  INTEGER, PARAMETER, PUBLIC               :: do_etrs=1,&
                                              do_cn=2,&
                                              do_em=3

  INTEGER, PARAMETER, PUBLIC               :: use_scf_wfn=1,  &
                                              use_restart_wfn=2,&
                                              use_rt_restart=3

  INTEGER, PARAMETER, PUBLIC               :: constant_env=1,  &
                                              gaussian_env=2, &
                                              ramp_env=3

  ! how to solve polarizable force fields
  INTEGER, PARAMETER, PUBLIC               :: do_fist_pol_none = 1, &
                                              do_fist_pol_sc = 2, &
                                              do_fist_pol_cg = 3

  ! some dm_ls_scf parameters
  INTEGER, PARAMETER, PUBLIC               :: ls_s_preconditioner_none=7, &
                                              ls_s_preconditioner_atomic=13,&
                                              ls_s_preconditioner_molecular=17

  INTEGER, PARAMETER, PUBLIC               :: ls_cluster_atomic=111, &
                                              ls_cluster_molecular=112

  INTEGER, PARAMETER, PUBLIC               :: ls_s_inversion_hotelling=3,&
                                              ls_s_inversion_sign_sqrt=5
  
  INTEGER, PARAMETER, PUBLIC               :: ls_scf_ns=17, ls_scf_trs4=18,&
                                              ls_scf_tc2=19

  INTEGER, PARAMETER, PUBLIC               :: ls_scf_line_search_3point=3,&
                                              ls_scf_line_search_3point_2d=6
  ! some ZMP paramenters
  INTEGER, PARAMETER, PUBLIC               :: use_coulomb = 1,&
                                              use_diff  = 2,&
                                              use_no  = 3

  ! parameters for ALMO methods
  INTEGER, PARAMETER, PUBLIC   :: almo_domain_layout_orbital=1, &
                                  almo_domain_layout_atomic=10, &
                                  almo_domain_layout_molecular=100

  INTEGER, PARAMETER, PUBLIC   :: almo_mat_distr_atomic=1, &
                                  almo_mat_distr_molecular=2

  INTEGER, PARAMETER, PUBLIC   :: almo_constraint_block_diagonal=0, &
                                  almo_constraint_distance=1,&
                                  almo_constraint_ao_overlap=2

  INTEGER, PARAMETER, PUBLIC   :: almo_deloc_none=0,&
                                  almo_deloc_x=10,&
                                  almo_deloc_qx=11,&
                                  almo_deloc_xk=12,&
                                  almo_deloc_qscf=13,&
                                  almo_deloc_full_scf=100,&
                                  almo_deloc_x_then_full_scf=110
  
  INTEGER, PARAMETER, PUBLIC   :: almo_scf_dm_sign=1,&
                                  almo_scf_diag=2,&
                                  almo_scf_pcg=3

  ! parameters for CT methods
  INTEGER, PARAMETER, PUBLIC   :: cg_zero=0,&
                                  cg_polak_ribiere=1,&
                                  cg_fletcher_reeves=2,&
                                  cg_hestenes_stiefel=3,&
                                  cg_fletcher=4,&
                                  cg_liu_storey=5,&
                                  cg_dai_yuan=6,& 
                                  cg_hager_zhang=7 
  
  INTEGER, PARAMETER, PUBLIC   :: tensor_orthogonal = 1, & 
                                  tensor_up_down = 2 
  
  INTEGER, PARAMETER, PUBLIC   :: virt_full = 1,&
                                  virt_minimal = 2,&
                                  virt_occ_size = 3,&
                                  virt_number = 4


  ! some MP2 parameters
  INTEGER, PARAMETER, PUBLIC               :: mp2_method_none=0,&
                                              mp2_method_laplace=2,&
                                              mp2_method_direct=4,&
                                              mp2_method_gpw=6,&
                                              ri_mp2_method_gpw=7,&
                                              ri_rpa_method_gpw=8,&
                                              ri_mp2_laplace=9,&
                                              mp2_ri_optimize_basis=10,&
                                              wfc_mm_style_gemm=11,&
                                              wfc_mm_style_syrk=12

  ! periodic RESP parameters
  INTEGER, PARAMETER, PUBLIC               :: do_resp_x_dir=0,&
                                              do_resp_y_dir=1,&
                                              do_resp_z_dir=2,& 
                                              do_resp_minus_x_dir=3,&
                                              do_resp_minus_y_dir=4,& 
                                              do_resp_minus_z_dir=5 

  ! basis optimization parameters
  INTEGER, PARAMETER, PUBLIC               :: do_opt_all=0,&
                                              do_opt_none=1,&
                                              do_opt_coeff=2,&
                                              do_opt_exps=3

  ! callgraph parameters
  INTEGER, PARAMETER, PUBLIC               :: callgraph_none=0,&
                                              callgraph_master=1,&
                                              callgraph_all=2 

  ! kg coloring parameters
  INTEGER, PARAMETER, PUBLIC               :: kg_color_dsatur=0,&
                                              kg_color_greedy=1
  ! kg kinetic energy embedding potential
  INTEGER, PARAMETER, PUBLIC               :: kg_tnadd_embed =100,&
                                              kg_tnadd_atomic=200

  ! swarm parameters
  INTEGER, PARAMETER, PUBLIC               :: swarm_do_glbopt=1

! global opt parameters
  INTEGER, PARAMETER, PUBLIC               :: glbopt_do_minhop  =1,&
                                              glbopt_do_mincrawl=2

! fm matrix multiplication
  INTEGER, PARAMETER, PUBLIC               :: do_pdgemm  =1,&
                                              do_dbcsr   =2

  ! Dispersion DFTB
  INTEGER, PARAMETER, PUBLIC               :: dispersion_uff=100,&
                                              dispersion_d3 =200

END MODULE input_constants
