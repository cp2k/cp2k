!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****h* cp2k/input_constants [1.0] *
!!
!!   NAME
!!     input_constants
!!
!!   FUNCTION
!!     collects all constants needed in input so that they can be used without
!!     circular dependencies
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     02.2005 created [fawzi&Joost]
!!
!!   SOURCE
!****************************************************************************
MODULE input_constants
  USE kinds,                           ONLY: default_string_length,&
                                             dp

  IMPLICIT NONE
  PRIVATE

  ! Tags for development/release versions
  INTEGER, PARAMETER, PUBLIC :: id_development_version=0,&
                                id_release_version=1

  ! Program Names
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(3), PUBLIC :: program_names=&
       (/"FARMING   ",&
         "TEST      ",&
         "CP2K      "/)

  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(3), PUBLIC :: program_name_desc=&
       (/"Runs N independent jobs in a single run",&
         "Do some benchmarking and testing       ",&
         "Runs one of the CP2K package           "/)
         
  INTEGER, PARAMETER, PUBLIC :: do_farming=1,&
                                do_test=2,&
                                do_cp2k=3

  INTEGER, PARAMETER, PUBLIC, DIMENSION(3) :: program_name_index=&
       (/do_farming, do_test, do_cp2k/)

  INTEGER, PARAMETER, PUBLIC :: do_full_density=0,&
                                do_spin_density=1

  INTEGER, PARAMETER, PUBLIC :: do_no_et=0,&
                                do_et_ddapc=1,&
                                do_et_becke=2
  ! Method Names
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: method_names=&
       (/"QS        ",&
         "FIST      ",&
         "KG        ",&
         "QMMM      ",&
         "EIP       ",&
         "QUICKSTEP ",&
         "KIM-GORDON",&
         "KIM_GORDON",&
         "EP        ",&
         "MIXED     "/)

  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: method_name_desc=&
       (/"Electronic structure methods (DFT, ...)",&
         "Molecular Mechanics                    ",&
         "Kim Gordon models (DFT)                ",&
         "Hybrid quantum classical               ",&
         "Empirical Interatomic Potential        ",&
         "Electronic structure methods (DFT, ...)",&    ! Aliases to the previous unique names
         "Kim Gordon models (DFT)                ",&    !
         "Kim Gordon models (DFT)                ",&    !
         "Non-orthogonal perturbation            ",&    !
         "Use a combination of two of the above  "/)    

  INTEGER, PARAMETER, PUBLIC :: do_qs=1,&
                                do_fist=2,&
                                do_kg=3,&
                                do_qmmm=4,&
                                do_eip=5,&
                                do_ep=6,&
                                do_mixed=7

  INTEGER, PARAMETER, PUBLIC, DIMENSION(10) :: method_name_index=&
       (/do_qs, do_fist, do_kg, do_qmmm, do_eip, &
         do_qs, do_kg, do_kg, do_ep, do_mixed/)


  REAL(KIND=dp), PARAMETER, PUBLIC :: RADIUS_QMMM_DEFAULT = 0.80_dp,&
                                      ALPHA_IMOMM_DEFAULT = 1.38_dp,&
                                      CHARGE_SCALE_FACTOR = 1.00_dp
  ! BAND
  INTEGER, PARAMETER, PUBLIC           :: band_md_opt=0,&
                                          band_diis_opt=1
  ! Weight_types
  INTEGER, PARAMETER, PUBLIC :: weight_type_unit=0,&
                                weight_type_mass=1

  ! Blacs type of distribution 
  INTEGER, PARAMETER, PUBLIC :: blacs_square=1, blacs_row=2, blacs_col=3

  ! EWALD methods       
  INTEGER, PARAMETER, PUBLIC :: do_ewald_none=1,&
                                do_ewald_ewald=2,&
                                do_ewald_pme=3,&
                                do_ewald_spme=4

  INTEGER, PARAMETER, PUBLIC :: do_constr_none=1,&
                                do_constr_atomic=2,&
                                do_constr_molec=3
  !
  ! MD Ensemble
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: ensemble_name_types=&
       (/"nve                ",&
         "nvt                ",&
         "npt_i              ",&
         "npt_f              ",&
         "msst               ",&
         "msst_damped        ",&
         "hydrostaticshock   ",&
         "isokin             ",&
         "reftraj            ",&
         "langevin           " /)
  INTEGER, PARAMETER, PUBLIC :: nve_ensemble=1,&
                                nvt_ensemble=2,&
                                npt_i_ensemble=3,&
                                npt_f_ensemble=4,&
                                nph_uniaxial_ensemble=5,&
                                nph_uniaxial_damped_ensemble=6,&
                                nph_ensemble=7,&
                                isokin_ensemble=8,&
                                reftraj_ensemble=9,&
                                langevin_ensemble=10
  INTEGER, PARAMETER, PUBLIC, DIMENSION(10) :: ensemble_name_index=&
       (/nve_ensemble,nvt_ensemble,npt_i_ensemble,&
         npt_f_ensemble,nph_uniaxial_ensemble,&
         nph_uniaxial_damped_ensemble,nph_ensemble,isokin_ensemble,reftraj_ensemble,&
         langevin_ensemble/)
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: ensemble_name_desc=&
   (/"constant energy (microcanonical)                                             ",&
     "constant temperature (canonical, Nose-Hoover)                                ",&
     "constant pressure using an isotropic cell                                    ",&
     "constant pressure using a flexible cell                                      ",&
     "simulate steady shock (uniaxial)                                             ",&
     "simulate steady shock (uniaxial) with extra viscosity                        ",&
     "simulate steady shock with hydrostatic pressure                              ",&
     "constant kinetic energy                                                      ",&
     "reading frames from a file called reftraj.xyz (e.g. for property calculation)",&
     "langevin dynamics (constant temperature)                                     "/)

  INTEGER, PARAMETER, PUBLIC :: do_thermo_same_as_part=100,&
                                do_thermo_nose=101,&
                                do_thermo_csvr=102
  !
  ! Run Types
  INTEGER, PARAMETER, PUBLIC :: none_run=1,&
                                energy_run=2,&
                                energy_force_run=3,&
                                mol_dyn_run=4,&
                                geo_opt_run=5,&
                                mon_car_run=6,&
                                electronic_spectra_run=7,&
                                qs_debug_run=8,&
                                bsse_run=9,&
                                linear_response_run=10,&
                                pint_run=11,&
                                vib_anal=12,&
                                do_band=13
  
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(21), PUBLIC :: run_type_names=&
       (/"NONE                     ",&
         "ENERGY                   ",&
         "ENERGY_FORCE             ",&
         "MD                       ",&
         "GEO_OPT                  ",&
         "MC                       ",&
         "SPECTRA                  ",&
         "DEBUG                    ",&
         "BSSE                     ",&
         "LR                       ",&
         "PINT                     ",&
         "WFN_OPT                  ",&   ! Aliases to the previous unique names
         "WAVEFUNCTION_OPTIMIZATION",&   !
         "MOLECULAR_DYNAMICS       ",&   !
         "GEOMETRY_OPTIMIZATION    ",&   !
         "MONTECARLO               ",&   !
         "ELECTRONIC_SPECTRA       ",&   !
         "LINEAR_RESPONSE          ",&   !
         "BAND                     ",&   
         "VIBRATIONAL_ANALYSIS     ",&   
         "NORMAL_MODES             " /)  !

  INTEGER, PARAMETER, PUBLIC, DIMENSION(21) :: run_type_index=&
       (/ none_run,& 
          energy_run,&
          energy_force_run,&
          mol_dyn_run,&
          geo_opt_run,&
          mon_car_run,&
          electronic_spectra_run,&
          qs_debug_run,&
          bsse_run,&
          linear_response_run,&
          pint_run,&
          energy_run,&
          energy_run,&
          mol_dyn_run,&
          geo_opt_run,&
          mon_car_run,&
          electronic_spectra_run,&
          linear_response_run,&
          do_band,&
          vib_anal,&
          vib_anal/)

  ! BAND
  INTEGER, PARAMETER, PUBLIC :: do_b_neb=0,&
                                do_it_neb=1,&
                                do_ci_neb=2,&
                                do_d_neb=3,&
                                do_sm=4,&
                                do_eb=5

  INTEGER, PARAMETER, PUBLIC :: atomic_guess=1,restart_guess=2,random_guess=3,&
       core_guess=4, densities_guess=5, history_guess=6, no_guess=7
  INTEGER, PARAMETER, PUBLIC :: pw_interp=1,&
       spline3_nopbc_interp=2,spline3_pbc_interp=3 ! to move somewhere else
  INTEGER, PARAMETER, PUBLIC :: do_gapw_gcs=1,&
       do_gapw_gct=2,&
       do_gapw_log=3
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(3), PUBLIC :: quadrature_name=&
       (/"Gauss-Chebyshev         ",&
         "Trans_Gauss-Chebyshev   ",&
         "LogTrans-Gauss-Chebyshev"/)

  INTEGER, PARAMETER, PUBLIC :: do_pwgrid_spherical=1,&
       do_pwgrid_ns_fullspace=2,&
       do_pwgrid_ns_halfspace=3

  INTEGER, PARAMETER, PUBLIC ::&
       do_extrpl_use_guess=0,&
       do_extrpl_use_prev_p=1,&
       do_extrpl_use_prev_rho_r=2,&
       do_extrpl_linear_wf=3,&
       do_extrpl_linear_p=4,&
       do_extrpl_linear_ps=5,&
       do_extrpl_use_prev_wf=6,&
       do_extrpl_ps=7,&
       do_extrpl_frozen=8,&
       do_extrpl_aspc=9

  INTEGER, PARAMETER, PUBLIC ::&
       do_method_gapw=1,&
       do_method_gapw_xc=2,&
       do_method_gpw=3,&
       do_method_kg_pol=4,&
       do_method_kg_nopol=5,&
       do_method_kg_gpw=6,&
       do_method_mndo=7,&
       do_method_am1=8,&
       do_method_pm3=9,&
       do_method_pdg=10,&
       do_method_rm1=11,&
       do_method_dftb=12,&
       do_method_pchg=13
   INTEGER, PARAMETER, PUBLIC :: use_unit_angstrom=0,&
       use_unit_bohr=1,&
       use_unit_scaled_angstrom=2,&
       use_unit_scaled_bohr=3
  INTEGER, PARAMETER, PUBLIC :: rsgrid_distributed=0,&
                                rsgrid_replicated=1,&
                                rsgrid_automatic=2

  INTEGER, PARAMETER, PUBLIC :: do_par_atom=0,&
                                do_par_grid=1
  ! Thermostats
  INTEGER, PARAMETER, PUBLIC :: do_thermo_no_communication=0,&
                                do_thermo_communication=1,&
                                do_thermo_only_master=2
  !
  ! pint coordinate transformations
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(2), PUBLIC :: transformation_name_types=&
  (/"NORMAL", "STAGE "/)
  INTEGER, PARAMETER, PUBLIC :: transformation_normal=1,&
                                transformation_stage=2
  INTEGER, PARAMETER, PUBLIC, DIMENSION(2) :: transformation_name_index=&
       (/transformation_normal,transformation_stage/)
  ! Free Energy methods
  INTEGER, PARAMETER, PUBLIC :: do_fe_meta=0,&
                                do_fe_ui=1,&
                                do_fe_ac=2

  INTEGER, PARAMETER, PUBLIC :: do_wall_none=0,&
                                do_wall_reflective=1,&
                                do_wall_quadratic=2,&
                                do_wall_gaussian=3

  INTEGER, PARAMETER, PUBLIC :: do_wall_p=0,&
                                do_wall_m=1
  !
  ! Output format files
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(3), PUBLIC :: output_name_types=&
  (/"ATOMIC",&
    "XMOL  ",&
    "DCD   "/)
  INTEGER, PARAMETER, PUBLIC :: dump_atomic=1,&
                                dump_xmol=2,&
                                dump_dcd=3
  INTEGER, PARAMETER, PUBLIC, DIMENSION(3) :: output_name_index=&
       (/dump_atomic,dump_xmol,dump_dcd/)

  INTEGER, PARAMETER, PUBLIC :: do_qmmm_none=0,&
                                do_qmmm_coulomb=1,&
                                do_qmmm_gauss=2,&
                                do_qmmm_swave=3,&
                                do_qmmm_link_imomm=1,&
                                do_qmmm_link_gho=2,&
                                do_qmmm_link_pseudo=3,&
                                do_center_none=0,&
                                do_center_grid=1

  INTEGER, PARAMETER, PUBLIC :: do_region_global=0, &
                                do_region_molecule=1,&
                                do_region_massive=2,&
                                do_region_none=3,&
                                do_region_defined=4,&
                                do_bondparm_covalent=0,&
                                do_bondparm_vdw=1,&
                                do_skip_12=0,&
                                do_skip_13=1,&
                                do_skip_14=2,&
                                do_conn_psf=1,&
                                do_conn_psf_u=2,&
                                do_conn_generate=3,&
                                do_conn_off=4,&
                                do_conn_mol_set=5,&
                                do_conn_g87=6,&
                                do_conn_g96=7,&
                                do_coord_pdb=0,&
                                do_coord_xyz=1,&
                                do_coord_g96=2,&
                                do_coord_off=3

  INTEGER, PARAMETER, PUBLIC :: do_ff_undef=0,&
                                do_ff_quartic=1,&
                                do_ff_g96=2,&
                                do_ff_charmm=3,&
                                do_ff_harmonic=4,&
                                do_ff_g87=5
  
  INTEGER, PARAMETER, PUBLIC :: do_stress_none=0,&
                                do_stress_analytical=1,&
                                do_stress_numerical=2,&
                                do_stress_diagonal_anal=3,&
                                do_stress_diagonal_numer=4
  
  INTEGER, PARAMETER, PUBLIC :: do_mc_traditional=0,&
                                do_mc_gemc_nvt=1,&
                                do_mc_gemc_npt=2,&
                                do_mc_virial=3

  INTEGER, PARAMETER, PUBLIC :: do_multipole_charge=0,&
                                do_multipole_dipole=1,&
                                do_multipole_quadrupole=2

  INTEGER, PARAMETER, PUBLIC :: do_loc_none=0,&
                                do_loc_jacobi=1,&
                                do_loc_crazy=2,&
                                do_loc_direct=3,&
                                do_loc_l1_norm_sd=4

  INTEGER, PARAMETER, PUBLIC :: do_loc_min=0,&
                                do_loc_max=1,&
                                op_loc_berry=0,&
                                op_loc_boys=1,&
                                op_loc_pipek=2,&
                                op_loc_l1_norm=3,&
                                state_loc_all=0,&
                                state_loc_range=1,&
                                state_loc_list=2,&
                                state_loc_none=3

  INTEGER, PARAMETER, PUBLIC :: use_mom_ref_com=0,&
                                use_mom_ref_coac=1,&
                                use_mom_ref_user=2,&
                                use_mom_ref_zero=3,&
                                use_mom_ref_nondef=4

  INTEGER, PARAMETER, PUBLIC :: nmr_gauge_r=1,&
                                nmr_gauge_r_and_Bader=2,&
                                nmr_gauge_r_and_step_function=3

  INTEGER, PARAMETER, PUBLIC :: epr_gauge_r=1,&
                                epr_gauge_r_and_Bader=2,&
                                epr_gauge_r_and_step_function=3

  INTEGER, PARAMETER, PUBLIC :: ot_mini_sd=1,ot_mini_cg=2,ot_mini_diis=3
  INTEGER, PARAMETER, PUBLIC :: ot_algo_taylor_or_diag=1,ot_algo_irac=2
  INTEGER, PARAMETER, PUBLIC :: ot_chol_irac=1,ot_poly_irac=2,ot_lwdn_irac=3
  INTEGER, PARAMETER, PUBLIC :: ls_none=1,ls_2pnt=2,ls_3pnt=3,ls_gold=4,ls_fit=5
  INTEGER, PARAMETER, PUBLIC :: ot_precond_none=0, &
                                ot_precond_full_single=1, &
                                ot_precond_full_kinetic=2, &
                                ot_precond_s_inverse=3, &
                                ot_precond_sparse_diag=4,&
                                ot_precond_sparse_firstorder=5,&
                                ot_precond_full_all=6,&
                                ot_precond_full_single_inverse=7
  INTEGER, PARAMETER,PUBLIC :: dyn_coeff_mass_global=1, dyn_coeff_mass_kind=2,&
   dyn_coeff_mass_atom=3, dyn_coeff_mass_coeff=4

  INTEGER, PUBLIC, PARAMETER :: default_lbfgs_method_id=1,&
                                default_bfgs_method_id=2,&
                                default_cg_method_id=3

  INTEGER, PUBLIC, PARAMETER :: default_minimization_method_id=1,&
                                default_ts_method_id=2

  INTEGER, PUBLIC, PARAMETER :: none_ts_method_id=0,&
                                default_dimer_method_id=1
  
  INTEGER, PUBLIC, PARAMETER :: do_first_rotation_step=0,&
                                do_second_rotation_step=1,&
                                do_third_rotation_step=2

  INTEGER, PARAMETER, PUBLIC :: silent_print_level=0, low_print_level=1,&
       medium_print_level=2, high_print_level=3, debug_print_level=4

  INTEGER, PARAMETER, PUBLIC :: xc_rho_no_smooth=0, xc_rho_spline2_smooth=1,&
       xc_rho_spline3_smooth=2, xc_rho_nn10=3, xc_rho_nn50=4
  INTEGER, PARAMETER, PUBLIC :: xc_deriv_pw=0, xc_deriv_spline2=1,&
       xc_deriv_spline3=2, xc_deriv_spline2_smooth=3,&
       xc_deriv_spline3_smooth=4, xc_deriv_nn10_smooth=5,&
       xc_deriv_nn50_smooth=6
  INTEGER, PARAMETER, PUBLIC :: xc_old_f_routine=0, xc_new_f_routine=1,&
       xc_test_lsd_f_routine=3, xc_debug_new_routine=4
  INTEGER, PARAMETER, PUBLIC :: xc_funct_no_shortcut=0,&
       xc_funct_blyp=1, xc_funct_pade=2, xc_funct_pbe=3,&
       xc_funct_tpss=4, xc_funct_hcth120=5, xc_funct_olyp=6,&
       xc_funct_bp=7, xc_none=8, xc_funct_xwpbe=9, xc_funct_b3lyp=10,&
       xc_funct_pbe0=11
  INTEGER, PARAMETER, PUBLIC :: sic_none=0, sic_mauri_us=1, sic_mauri_spz=2, sic_ad=3, sic_eo=4
  INTEGER, PARAMETER, PUBLIC :: sic_list_all=1, sic_list_unpaired=2 
  INTEGER, PARAMETER, PUBLIC :: tddfpt_singlet=0, tddfpt_triplet=1
  INTEGER, PARAMETER, PUBLIC :: tddfpt_lanczos=0, tddfpt_davidson=1
  INTEGER, PARAMETER, PUBLIC :: oe_none=0, oe_lb=1, oe_gllb=2, oe_saop=3, oe_sic=4
  INTEGER, PARAMETER, PUBLIC :: no_excitations=0, tddfpt_excitations=1
  
  INTEGER, PARAMETER, PUBLIC :: do_add=1, do_remove=2

  INTEGER, PARAMETER, PUBLIC :: no_precond=0,precond_spl3_aint=1, precond_spl3_1=2,&
       precond_spl3_aint2=3, precond_spl3_2=4, precond_spl3_3=5

  INTEGER, PARAMETER, PUBLIC :: PERIODIC3D = 1000, WAVELET3D = 1001,&
       ANALYTIC2D = 1001,  ANALYTIC1D = 1002,  ANALYTIC0D = 1003,&
       MT2D = 1101,        MT1D = 1102,        MT0D = 1103,&
       HOCKNEY2D = 1201,   HOCKNEY1D = 1202,   HOCKNEY0D = 1203,&
       MULTIPOLE2D = 1301, MULTIPOLE1D = 1302, MULTIPOLE0D = 1303,&
       WAVELET2D = 1401,   WAVELET1D = 1402,   WAVELET0D = 1303

  INTEGER, PARAMETER, PUBLIC :: use_none=0,&
                                use_periodic = 1,&
                                use_analytic = 2,&
                                use_mt=3,&
                                use_hockney=5,&
                                use_multipole=4,&
                                use_wavelet=6

  INTEGER, PARAMETER, PUBLIC :: use_rs_grid=0,&
                                use_gs_grid = 1

  INTEGER, PARAMETER, PUBLIC :: xgga_b88=101,xgga_b88x=102,&
                                xgga_pw86=103, xgga_pw91=104,&
                                xgga_pbex=105, xgga_revpbe=106, &
                                xgga_optx=107,&
                                xgga_hcth=108, xgga_tpss=109, &
                                xalpha=110, xgga_pbe=111, xgga_opt=112
  INTEGER, PARAMETER, PUBLIC :: cgga_lyp=201, c_pade=202, cgga_hcth93=203, &
                                cgga_hcth120=204, cgga_hcth147=205,&
                                cgga_hcth407=206, cgga_cs1=207, cgga_pw86=208,&
                                c_pw92=209, c_pw92dmc=210, c_pw92vmc=211, &
                                c_pz=212, c_pzdmc=213, c_pzvmc=214, c_vwn=215,&
                                cgga_pbe=216

  INTEGER, PARAMETER, PUBLIC :: ke_ol1=301, ke_ol2=302,ke_llp=303,ke_pw86=304,&
                                ke_pw91=305, ke_lc=306, ke_t92=307, ke_pbe=308, &
                                ke_revpbe=309, ke_tf=310, ke_tfw=311

  INTEGER, PARAMETER, PUBLIC :: pw_orig=209,pw_dmc=210,pw_vmc=211
  
  INTEGER, PARAMETER, PUBLIC :: pz_orig=212, pz_dmc=213, pz_vmc=214

  INTEGER, PARAMETER, PUBLIC :: use_perd_x=0,&
                                use_perd_y=1,&
                                use_perd_z=2,&
                                use_perd_xy=3,&
                                use_perd_xz=4,&
                                use_perd_yz=5,&
                                use_perd_xyz=6,&
                                use_perd_none=7

  INTEGER, PARAMETER, PUBLIC :: do_clv_geo_center=0,&
                                do_clv_fix_point =1,&
                                do_clv_xyz       =0,&
                                do_clv_x         =1,&
                                do_clv_y         =2,&
                                do_clv_z         =3,&
                                do_clv_xy        =4,&
                                do_clv_xz        =5,&
                                do_clv_yz        =6

  INTEGER, PARAMETER, PUBLIC :: xc_pot_energy_none=0, &
                                xc_pot_energy_xc_functional=1, &
                                xc_pot_energy_sum_eigenvalues=2

  INTEGER, PARAMETER, PUBLIC :: xc_pot_none = 1000, xc_pot_saop = 1001
 
  ! X-Ray spectra Calculation
  INTEGER, PARAMETER, PUBLIC :: xas_none=0, xas_tp_hh=1, xas_tp_fh=2,&
                                xes_tp_val=3, xas_tp_xhh=4,xas_tp_xfh=5,&
                                xas_dscf=6 
  INTEGER, PARAMETER, PUBLIC :: xas_1s_type=1, xas_2s_type=2, xas_2p_type=3,&
                                xas_3s_type=4, xas_3p_type=5, xas_3d_type=6,&
                                xas_4s_type=7, xas_4p_type=8, xas_4d_type=9, xas_4f_type=10
  INTEGER, PARAMETER, PUBLIC :: xas_dip_len=1, xas_dip_vel=2, xas_dip_len2=3
  INTEGER, PARAMETER, PUBLIC :: xas_scf_default=0,xas_scf_general=1

  ! Linear Response for properties
  INTEGER, PARAMETER, PUBLIC :: lr_none=0, lr_chemshift=1, lr_gtensor=2
  INTEGER, PARAMETER, PUBLIC :: lr_cg=1, lr_diis=2, lr_sd=3
  INTEGER, PARAMETER, PUBLIC :: lr_precond_none=0, &
                                lr_precond_single=1, &
                                lr_precond_h0_inverse=2, &
                                lr_precond_s_inverse=3,&
                                lr_precond_kinetic=4, &
                                lr_precond_full_all=6

  INTEGER, PARAMETER, PUBLIC :: xc_pbe_orig=11,xc_rev_pbe=12

  ! EIP models
  INTEGER, PARAMETER, PUBLIC :: use_bazant_eip=1, &
                                use_lenosky_eip=2

  ! ddapc restraint forms
  INTEGER, PARAMETER, PUBLIC :: do_ddapc_restraint=773, do_ddapc_constraint=774

  ! outer scf types
  INTEGER, PARAMETER, PUBLIC :: outer_scf_ddapc_constraint=123, &
                                outer_scf_s2_constraint=124,&
                                outer_scf_becke=125, &
                                outer_scf_none=126

  ! outer scf optimizers
  INTEGER, PARAMETER, PUBLIC :: outer_scf_optimizer_sd=1001, &
                                outer_scf_optimizer_diis=1002, &
                                outer_scf_optimizer_none=1003, &
                                outer_scf_optimizer_bisect=1004

  ! s2 restraint forms
  INTEGER, PARAMETER, PUBLIC :: do_s2_restraint=872, do_s2_constraint=873

  ! ROKS schemes
  INTEGER, PARAMETER, PUBLIC :: general_roks   = 1,&
                                high_spin_roks = 2

  ! VWN xc functional flavours
  INTEGER, PARAMETER, PUBLIC :: do_vwn5=1,&
                                do_vwn3=2

  ! mixing of force envs
  INTEGER, PARAMETER, PUBLIC :: mix_linear_combination=701, &
                                mix_minimum=702, &
                                mix_coupled=703, &
                                mix_restrained=704,&
                                mix_generic=705

  ! cost models for distribution 2d
  INTEGER, PARAMETER, PUBLIC :: model_block_count=453, &
                                model_block_surface=454, &
                                model_block_lmax=455


  ! HFX potential type
  INTEGER, PARAMETER, PUBLIC :: do_hf_potential_coulomb=1,&
                                do_hf_potential_short=2,&
                                do_hf_potential_long=3,&
                                do_hf_potential_mix_cl=4

  ! mode selctive vibrational analysis
  INTEGER, PARAMETER, PUBLIC :: ms_guess_bfgs=1, &
                                ms_guess_atomic=2, &
                                ms_guess_restart=3,&
                                no_ms_precond=4, &
                                bfgs_ms_precond=5     

  ! FFT library
  ! these might need sync with fft_lib.F
  INTEGER, PARAMETER, PUBLIC :: do_fft_sg=1,&
                                do_fft_fftw2=2,&
                                do_fft_fftw3=3,&
                                do_fft_essl=4,&
                                do_fft_acml=5,&
                                do_fft_mkl=6,&
                                do_fft_sci=7
  ! fftw will automatically be either FFTW2 or FFTW3 
#if defined ( __FFTW2 )
  INTEGER, PARAMETER, PUBLIC :: do_fft_fftw=do_fft_fftw2
#else
  INTEGER, PARAMETER, PUBLIC :: do_fft_fftw=do_fft_fftw3
#endif

END MODULE input_constants
