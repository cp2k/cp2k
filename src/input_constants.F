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
  USE kinds, ONLY: dp
  implicit none
  private

  INTEGER, PARAMETER, public :: do_qs=1,do_fist=2,do_kg=3,do_qmmm=4,do_atom=5,&
       do_farming=6,do_test=8
  REAL(KIND=dp), PARAMETER, PUBLIC :: RADIUS_QMMM_DEFAULT = 0.80_dp,&
                                      ALPHA_IMOMM_DEFAULT = 1.38_dp,&
                                      CHARGE_SCALE_FACTOR = 1.00_dp
  INTEGER, PARAMETER, PUBLIC :: nve_ensemble=101,&
                                nvt_ensemble=102,&
                                npt_i_ensemble=103,&
                                npt_f_ensemble=104,&
                                nph_ensemble=105,&
                                nph_uniaxial=106,&
                                nph_uniaxial_damped_ensamble=107,&
                                isokin_ensemble=108
  INTEGER, PARAMETER, PUBLIC :: nhc_global=101,nhc_massive=102
  INTEGER, PARAMETER, PUBLIC :: md_init_start=101,md_pos_restart=102,&
       md_pos_vel_restart=103
  INTEGER, PARAMETER, public :: energy_run=1,energy_force_run=2,&
       md_run=3, geo_opt_run=4,mc_run=5, electronic_spectra_run=6, qs_debug=7
  INTEGER, PARAMETER, public :: atomic_guess=1,restart_guess=2,random_guess=3,&
       core_guess=4
  INTEGER, PARAMETER, public :: pw_interp=1,linear_interp=2,&
       spline3_nopbc_interp=3,spline3_pbc_interp=4 ! to move somewhere else
  integer, parameter, public :: do_gapw_gcs=1,&
       do_gapw_gct=2,&
       do_gapw_log=3
  integer, parameter, public :: do_pwgrid_spherical=1,&
       do_pwgrid_ns_fullspace=2,&
       do_pwgrid_ns_halfspace=3
  INTEGER, PARAMETER, PUBLIC :: do_etrxpl_use_guess=0,&
       do_etrxpl_use_prev_p=1,&
       do_etrxpl_use_prev_rho_r=2,&
       do_etrxpl_linear_wf=3,&
       do_etrxpl_linear_p=4,&
       do_etrxpl_linear_ps=5,&
       do_etrxpl_use_prev_wf=6,&
       do_etrxpl_ps=7
  INTEGER, PARAMETER, PUBLIC :: do_method_gapw=1,&
       do_method_gpw=2,&
       do_method_kg_pol=3,&
       do_method_kg_nopol=4,&
       do_method_kg_gpw=5,&
       do_method_eht=6,&
       do_method_mndo=7,&
       do_method_am1=8,&
       do_method_pm3=9
   integer, parameter, public :: use_unit_angstrom=0,&
       use_unit_bohr=1,&
       use_unit_scaled_angstrom=2,&
       use_unit_scaled_bohr=3
  integer, parameter, public :: do_rsgrid_distributed=0,&
       do_rsgrid_replicated=1

  integer, parameter, public :: do_extrpl_ps=0,&
       do_extrpl_use_guess=1,&
       do_extrpl_use_prev_p=2,&
       do_extrpl_use_prev_wf=3,&
       do_extrpl_use_prev_rho_r=4,&
       do_extrpl_linear_wf=5,&
       do_extrpl_linear_p=6,&
       do_extrpl_linear_ps=7

  INTEGER, PARAMETER, PUBLIC :: dump_xyz=0,&
                                dump_crd=1

  integer, parameter, public :: do_qmmm_none=0,&
                                do_qmmm_grid=1,&
                                do_qmmm_gauss=2,&
                                do_qmmm_link_imomm=1,&
                                do_qmmm_link_gho=2,&
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

  INTEGER, PARAMETER, PUBLIC :: do_ewald=0,&
                                do_pme=1,&
                                do_spme=2

  INTEGER, PARAMETER, public :: do_loc_jacobi=0,&
                                do_loc_crazy=1,&
                                do_loc_direct=2,&
                                do_loc_none=3

  INTEGER, PARAMETER, public :: ot_mini_sd=1,ot_mini_cg=2,ot_mini_diis=3
  INTEGER, PARAMETER, public :: ot_ls_none=1,ot_ls_2pnt=2,ot_ls_3pnt=3,ot_ls_gold=4
  INTEGER, PARAMETER, public :: ot_precond_none=0, &
                        ot_precond_full_single=1, &
                        ot_precond_full_kinetic=2, &
                        ot_precond_s_inverse=3, &
                        ot_precond_sparse_diag=4
  INTEGER, PARAMETER,PUBLIC :: dyn_coeff_mass_global=1, dyn_coeff_mass_kind=2,&
   dyn_coeff_mass_atom=3, dyn_coeff_mass_coeff=4
  INTEGER, PARAMETER, PUBLIC :: kg_coeff_init_start=1,kg_coeff_pos_restart=2, &
   kg_coeff_pos_vel_restart=3

  INTEGER, PUBLIC, PARAMETER :: default_lbfgs_method_id=1,&
                                default_bfgs_method_id=2


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
  INTEGER, PARAMETER, PUBLIC :: oe_none=0, oe_saop=1, oe_sic=2
  INTEGER, PARAMETER, PUBLIC :: no_excitations=0, tddfpt_excitations=1
  
  
  INTEGER, PARAMETER, PUBLIC :: no_precond=0,precond_spl3_aint=1, precond_spl3_1=2

  INTEGER, PARAMETER, PUBLIC :: PERIODIC3D = 1000, ANALYTIC2D = 1001,  &
                        ANALYTIC1D = 1002, ANALYTIC0D = 1003,  &
                        MT2D = 1101, MT1D = 1102, MT0D = 1103,&
                        HOCKNEY2D = 1201, HOCKNEY1D = 1202, HOCKNEY0D = 1203

  INTEGER, PARAMETER, PUBLIC :: use_periodic = 0,&
                                use_analytic = 1,&
                                use_mt=2,&
                                use_hockney=3

  INTEGER, PARAMETER, PUBLIC :: xgga_b88=1, xgga_pw86=2, xgga_pw91=3, xgga_pbe=4,&
       xgga_revpbe=5, xgga_opt=6

  INTEGER, PARAMETER, PUBLIC :: ke_ol1=1, ke_ol2=2,ke_llp=3,ke_pw86=4,&
       ke_pw91=5, ke_lc=6,ke_t92=7

  INTEGER, PARAMETER, PUBLIC :: pw_orig=1,pw_dmc=2,pw_vmc=3
  
  INTEGER, PARAMETER, PUBLIC :: pz_orig=1, pz_dmc=2, pz_vmc=3

  INTEGER, PARAMETER, PUBLIC :: use_perd_x=0,&
                                use_perd_y=1,&
                                use_perd_z=2,&
                                use_perd_xy=3,&
                                use_perd_xz=4,&
                                use_perd_yz=5,&
                                use_perd_xyz=6,&
                                use_perd_none=7
 
  INTEGER, PARAMETER, PUBLIC :: xas_none=0, xas_tp_hh=1, xas_tp_fh=2,  xas_tddft=3 

END MODULE input_constants
