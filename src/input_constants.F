!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
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
  USE kinds, ONLY: dp,&
                   default_string_length
  IMPLICIT NONE
  PRIVATE

  ! Program Names
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: program_names=&
       (/"QS        ",&
         "FIST      ",&
         "KG        ",&
         "QMMM      ",&
         "ATOM      ",&
         "FARMING   ",&
         "TEST      ",&
         "QUICKSTEP ",&
         "KIM-GORDON",&
         "KIM_GORDON" /)

  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: program_name_desc=&
       (/"Electronic structure methods (DFT, ...)",&
         "Molecular Mechanics                    ",&
         "Kim Gordon models (DFT)                ",&
         "Hybrid quantum classical               ",&
         "Atomic calculations                    ",&
         "Runs N independent jobs in a single run",&
         "Do some benchmarking and testing       ",&
         "Electronic structure methods (DFT, ...)",&    ! Aliases to the previous unique names
         "Kim Gordon models (DFT)                ",&    !
         "Kim Gordon models (DFT)                " /)   !

  INTEGER, PARAMETER, public :: do_qs=1,&
                                do_fist=2,&
                                do_kg=3,&
                                do_qmmm=4,&
                                do_atom=5,&
                                do_farming=6,&
                                do_test=7

  INTEGER, PARAMETER, PUBLIC, DIMENSION(10) :: program_name_index=&
       (/do_qs, do_fist, do_kg, do_qmmm, do_atom, do_farming, do_test, do_qs, do_kg, do_kg/)

  REAL(KIND=dp), PARAMETER, PUBLIC :: RADIUS_QMMM_DEFAULT = 0.80_dp,&
                                      ALPHA_IMOMM_DEFAULT = 1.38_dp,&
                                      CHARGE_SCALE_FACTOR = 1.00_dp
  ! Blacs type of distribution 
  INTEGER, PARAMETER, PUBLIC :: blacs_square=1, blacs_row=2, blacs_col=3

  ! EWALD methods
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(4), PUBLIC :: ewald_name_types=&
       (/"none     ",&
         "ewald    ",&
         "pme      ",&
         "spme     " /)
  INTEGER, PARAMETER, PUBLIC :: do_ewald_none=1,&
                                do_ewald_ewald=2,&
                                do_ewald_pme=3,&
                                do_ewald_spme=4

  INTEGER, PARAMETER, PUBLIC, DIMENSION(4) :: ewald_name_index=&
       (/do_ewald_none, do_ewald_ewald, do_ewald_pme, do_ewald_spme/)

  !
  ! MD Ensemble
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: ensemble_name_types=&
       (/"nve                ",&
         "nvt                ",&
         "npt_i              ",&
         "npt_f              ",&
         "nph                ",&
         "nph_uniaxial       ",&
         "nph_uniaxial_damped",&
         "isokin             ",&
         "reftraj            ",&
         "langevin           " /)
  INTEGER, PARAMETER, PUBLIC :: nve_ensemble=1,&
                                nvt_ensemble=2,&
                                npt_i_ensemble=3,&
                                npt_f_ensemble=4,&
                                nph_ensemble=5,&
                                nph_uniaxial_ensemble=6,&
                                nph_uniaxial_damped_ensemble=7,&
                                isokin_ensemble=8,&
                                reftraj_ensemble=9,&
                                langevin_ensemble=10
  INTEGER, PARAMETER, PUBLIC, DIMENSION(10) :: ensemble_name_index=&
       (/nve_ensemble,nvt_ensemble,npt_i_ensemble,&
         npt_f_ensemble,nph_ensemble,nph_uniaxial_ensemble,&
         nph_uniaxial_damped_ensemble,isokin_ensemble,reftraj_ensemble,&
         langevin_ensemble/)

  INTEGER, PARAMETER, PUBLIC :: nhc_global=101,nhc_massive=102
  !
  ! MD restart
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(4), PUBLIC :: restart_name_types=&
       (/"init  ",&
         "pos   ",&
         "posvel",&
         "all   "/)
  INTEGER, PARAMETER, PUBLIC :: md_init_start=1,&
                                md_pos_restart=2,&
                                md_pos_vel_restart=3,&
                                md_all_restart=4
  INTEGER, PARAMETER, PUBLIC, DIMENSION(4) :: restart_name_index=&
       (/md_init_start,md_pos_restart, md_pos_vel_restart, md_all_restart/)

  !
  ! Run Types
  INTEGER, PARAMETER, public :: none_run=1,&
                                energy_run=2,&
                                energy_force_run=3,&
                                mol_dyn_run=4,&
                                geo_opt_run=5,&
                                mon_car_run=6,&
                                electronic_spectra_run=7,&
                                qs_debug_run=8,&
                                bsse_run=9,&
                                linear_response_run=10,&
                                ep_run=11,&
                                hmc_run=12
  
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(22), PUBLIC :: run_type_names=&
       (/"NONE                     ",&
         "ENERGY                   ",&
         "ENERGY_FORCE             ",&
         "MD                       ",&
         "GEO_OPT                  ",&
         "MC                       ",&
         "LDMC                     ",&
         "HMC                      ",&
         "SPECTRA                  ",&
         "DEBUG                    ",&
         "BSSE                     ",&
         "LR                       ",&
         "EP                       ",&
         "WFN_OPT                  ",&   ! Aliases to the previous unique names
         "WAVEFUNCTION_OPTIMIZATION",&   !
         "MOLECULAR_DYNAMICS       ",&   !
         "GEOMETRY_OPTIMIZATION    ",&   !
         "MONTECARLO               ",&   !
         "LANGEVIN_MONTECARLO      ",&   !
         "HYBRID_MONTECARLO        ",&   !
         "ELECTRONIC_SPECTRA       ",&   !
         "LINEAR_RESPONSE          " /)  !

  INTEGER, PARAMETER, PUBLIC, DIMENSION(22) :: run_type_index=&
       (/ none_run,& 
          energy_run,&
          energy_force_run,&
          mol_dyn_run,&
          geo_opt_run,&
          mon_car_run,&
          hmc_run,&
          hmc_run,&
          electronic_spectra_run,&
          qs_debug_run,&
          bsse_run,&
          linear_response_run,&
          ep_run,&
          energy_run,&
          energy_run,&
          mol_dyn_run,&
          geo_opt_run,&
          mon_car_run,&
          hmc_run,&
          hmc_run,&
          electronic_spectra_run,&
          linear_response_run /)

  INTEGER, PARAMETER, public :: atomic_guess=1,restart_guess=2,random_guess=3,&
       core_guess=4
  INTEGER, PARAMETER, public :: pw_interp=1,linear_interp=2,&
       spline3_nopbc_interp=3,spline3_pbc_interp=4 ! to move somewhere else
  integer, parameter, public :: do_gapw_gcs=1,&
       do_gapw_gct=2,&
       do_gapw_log=3
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(3), PUBLIC :: quadrature_name=&
       (/"Gauss-Chebyshev         ",&
         "Trans_Gauss-Chebyshev   ",&
         "LogTrans-Gauss-Chebyshev"/)

  integer, parameter, public :: do_pwgrid_spherical=1,&
       do_pwgrid_ns_fullspace=2,&
       do_pwgrid_ns_halfspace=3

  INTEGER, PARAMETER, PUBLIC :: do_extrpl_use_guess=0,&
       do_extrpl_use_prev_p=1,&
       do_extrpl_use_prev_rho_r=2,&
       do_extrpl_linear_wf=3,&
       do_extrpl_linear_p=4,&
       do_extrpl_linear_ps=5,&
       do_extrpl_use_prev_wf=6,&
       do_extrpl_ps=7,&
       do_extrpl_frozen=8
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(10), PUBLIC :: elstr_method_name=&
      (/"GAPW    ", "GAPW_XC ","GPW     ",&
        "KG_POL  ", "KG_NOPOL", "KG_GPW  ", &
        "EHT     ", "MNDO    ", "AM1     ", "PM3     "/)
  INTEGER, PARAMETER, PUBLIC :: do_method_gapw=1,&
       do_method_gapw_xc=2,&
       do_method_gpw=3,&
       do_method_kg_pol=4,&
       do_method_kg_nopol=5,&
       do_method_kg_gpw=6,&
       do_method_eht=7,&
       do_method_mndo=8,&
       do_method_am1=9,&
       do_method_pm3=10
   integer, parameter, public :: use_unit_angstrom=0,&
       use_unit_bohr=1,&
       use_unit_scaled_angstrom=2,&
       use_unit_scaled_bohr=3
  integer, parameter, public :: do_rsgrid_distributed=0,&
       do_rsgrid_replicated=1

  !
  ! Output format files
  !
  CHARACTER (LEN=default_string_length), PARAMETER, DIMENSION(5), PUBLIC :: output_name_types=&
  (/"XYZ   ",&
    "CRD   ",&
    "ATOMIC",&
    "XMOL  ",&
    "DCD   "/)
  INTEGER, PARAMETER, PUBLIC :: dump_xyz=1,&
                                dump_crd=2,&
                                dump_atomic=3,&
                                dump_xmol=4,&
                                dump_dcd=5
  INTEGER, PARAMETER, PUBLIC, DIMENSION(5) :: output_name_index=&
       (/dump_xyz,dump_crd,dump_atomic,dump_xmol,dump_dcd/)

  integer, parameter, public :: do_qmmm_none=0,&
                                do_qmmm_grid=1,&
                                do_qmmm_gauss=2,&
                                do_qmmm_link_imomm=1,&
                                do_qmmm_link_gho=2,&
                                do_qmmm_link_pseudo=3,&
                                do_center_none=0,&
                                do_center_grid=1

  INTEGER, PARAMETER, PUBLIC :: do_nhcopt_global=0, &
                                do_nhcopt_molecule=1,&
                                do_nhcopt_massive=2,&
                                do_nhcopt_none=3,&
                                do_bondparm_cov=0,&
                                do_bondparm_vdw=1,&
                                do_skip_ei_12=0,&
                                do_skip_ei_13=1,&
                                do_skip_ei_14=2,&
                                do_skip_vdw_12=0,&
                                do_skip_vdw_13=1,&
                                do_skip_vdw_14=2,&
                                do_conn_psf=0,&
                                do_conn_gromos=1,&
                                do_conn_generate=2,&
                                do_conn_off=3,&
                                do_coord_pdb=0,&
                                do_coord_chm=1,&
                                do_coord_off=2

  INTEGER, PARAMETER, PUBLIC :: do_mc_traditional=0,&
                                do_mc_gemc=1

  INTEGER, PARAMETER, public :: do_loc_jacobi=0,&
                                do_loc_crazy=1,&
                                do_loc_direct=2,&
                                do_loc_none=3
  INTEGER, PARAMETER, public :: do_loc_min=0,&
                                do_loc_max=1,&
                                op_loc_berry=0,&
                                op_loc_boys=1,&
                                op_loc_pipek=2,&
                                state_loc_all=0,&
                                state_loc_upton=1,&
                                state_loc_fromitoj=2,&
                                state_loc_list=3,&
                                state_loc_none=4

  INTEGER, PARAMETER, public :: ot_mini_sd=1,ot_mini_cg=2,ot_mini_diis=3
  INTEGER, PARAMETER, public :: ot_ls_none=1,ot_ls_2pnt=2,ot_ls_3pnt=3,ot_ls_gold=4
  INTEGER, PARAMETER, public :: ot_precond_none=0, &
                        ot_precond_full_single=1, &
                        ot_precond_full_kinetic=2, &
                        ot_precond_s_inverse=3, &
                        ot_precond_sparse_diag=4,&
                        ot_precond_sparse_firstorder=5
  INTEGER, PARAMETER,PUBLIC :: dyn_coeff_mass_global=1, dyn_coeff_mass_kind=2,&
   dyn_coeff_mass_atom=3, dyn_coeff_mass_coeff=4
  INTEGER, PARAMETER, PUBLIC :: kg_coeff_init_start=1,kg_coeff_pos_restart=2, &
   kg_coeff_pos_vel_restart=3

  INTEGER, PUBLIC, PARAMETER :: default_lbfgs_method_id=1,&
                                default_bfgs_method_id=2,&
                                default_cg_method_id=3


  INTEGER, PARAMETER, PUBLIC :: silent_print_level=0, low_print_level=1,&
       medium_print_level=2, high_print_level=3, debug_print_level=4

  INTEGER, PARAMETER, PUBLIC :: xc_rho_no_smooth=0, xc_rho_spline2_smooth=1,&
       xc_rho_spline3_smooth=2, xc_rho_nn10=3, xc_rho_nn50=4
  INTEGER, PARAMETER, PUBLIC :: xc_deriv_pw=0, xc_deriv_spline2=1,&
       xc_deriv_spline3=2, xc_deriv_spline2_smooth=3,&
       xc_deriv_spline3_smooth=4, xc_deriv_nn10_smooth=5,&
       xc_deriv_fd=6, xc_deriv_nn50_smooth=7
  INTEGER, PARAMETER, PUBLIC :: xc_old_f_routine=0, xc_new_f_routine=1,&
       xc_test_lsd_f_routine=3, xc_debug_new_routine=4
  INTEGER, PARAMETER, PUBLIC :: sic_none=0, sic_mauri_us=1, sic_mauri_spz=2, sic_ad=3
  INTEGER, PARAMETER, PUBLIC :: tddfpt_singlet=0, tddfpt_triplet=1
  INTEGER, PARAMETER, PUBLIC :: tddfpt_lanczos=0, tddfpt_davidson=1
  INTEGER, PARAMETER, PUBLIC :: oe_none=0, oe_lb=1, oe_gllb=2, oe_saop=3, oe_sic=4
  INTEGER, PARAMETER, PUBLIC :: no_excitations=0, tddfpt_excitations=1
  
  
  INTEGER, PARAMETER, PUBLIC :: no_precond=0,precond_spl3_aint=1, precond_spl3_1=2,&
       precond_spl3_aint2=3, precond_spl3_2=4, precond_spl3_3=5

  INTEGER, PARAMETER, PUBLIC :: PERIODIC3D = 1000, ANALYTIC2D = 1001,  &
                        ANALYTIC1D = 1002, ANALYTIC0D = 1003,  &
                        MT2D = 1101, MT1D = 1102, MT0D = 1103,&
                        HOCKNEY2D = 1201, HOCKNEY1D = 1202, HOCKNEY0D = 1203,&
                        MULTIPOLE2D = 1301, MULTIPOLE1D = 1302, MULTIPOLE0D = 1303

  INTEGER, PARAMETER, PUBLIC :: use_periodic = 0,&
                                use_analytic = 1,&
                                use_mt=2,&
                                use_hockney=3,&
                                use_multipole=4

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
 
  ! X-Ray spectra Calculation
  INTEGER, PARAMETER, PUBLIC :: xas_none=0, xas_tp_hh=1, xas_tp_fh=2,  xas_tddft=3, xas_dscf=4 
  INTEGER, PARAMETER, PUBLIC :: xas_1s_type=1, xas_2s_type=2, xas_2p_type=3,&
                                xas_3s_type=4, xas_3p_type=5, xas_3d_type=6,&
                                xas_4s_type=7, xas_4p_type=8, xas_4d_type=9, xas_4f_type=10
  INTEGER, PARAMETER, PUBLIC :: xas_dip_len=1, xas_dip_vel=2

  ! Linear Response for properties
  INTEGER, PARAMETER, PUBLIC :: lr_none=0, lr_chemshift=1
  INTEGER, PARAMETER, PUBLIC :: lr_cg=1, lr_diis=2, lr_sd=3
  INTEGER, PARAMETER, PUBLIC :: lr_precond_none=0, &
                                lr_precond_single=1, &
                                lr_precond_kinetic=2, &
                                lr_precond_s_inverse=3

END MODULE input_constants
