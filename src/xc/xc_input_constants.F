!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief input constants for xc
! *****************************************************************************
MODULE xc_input_constants
  INTEGER, PARAMETER, PUBLIC               :: xc_deriv_pw=0,&
                                              xc_deriv_spline2=1,&
                                              xc_deriv_spline3=2,&
                                              xc_deriv_spline2_smooth=3,&
                                              xc_deriv_spline3_smooth=4,&
                                              xc_deriv_nn10_smooth=5,&
                                              xc_deriv_nn50_smooth=6,&
                                              xc_deriv_collocate=7

  INTEGER, PARAMETER, PUBLIC               :: xc_rho_no_smooth=0,&
                                              xc_rho_spline2_smooth=1,&
                                              xc_rho_spline3_smooth=2,&
                                              xc_rho_nn10=3,&
                                              xc_rho_nn50=4

  !VWN xc functional flavours
  INTEGER, PARAMETER, PUBLIC               :: do_vwn5=1,&
                                              do_vwn3=2


  INTEGER, PARAMETER, PUBLIC               :: xc_pbe_orig=11,&
                                              xc_pbe_rev=12,&
                                              xc_pbe_sol=13

  INTEGER, PARAMETER, PUBLIC               :: xc_old_f_routine=0,&
                                              xc_new_f_routine=1,&
                                              xc_test_lsd_f_routine=3,&
                                              xc_debug_new_routine=4

  INTEGER, PARAMETER, PUBLIC               :: xgga_b88=101,&
                                              xgga_b88x=102,&
                                              xgga_pw86=103,&
                                              xgga_pw91=104,&
                                              xgga_pbex=105,&
                                              xgga_revpbe=106,&
                                              xgga_optx=107,&
                                              xgga_hcth=108,&
                                              xgga_tpss=109,&
                                              xalpha=110,&
                                              xgga_pbe=111,&
                                              xgga_opt=112,&
                                              xgga_ev93=113

  INTEGER, PARAMETER, PUBLIC               :: cgga_lyp=201,&
                                              c_pade=202,&
                                              cgga_hcth93=203,&
                                              cgga_hcth120=204,&
                                              cgga_hcth147=205,&
                                              cgga_hcth407=206,&
                                              cgga_cs1=207,&
                                              cgga_pw86=208,&
                                              c_pw92=209,&
                                              c_pw92dmc=210,&
                                              c_pw92vmc=211,&
                                              c_pz=212,&
                                              c_pzdmc=213,&
                                              c_pzvmc=214,&
                                              c_vwn=215,&
                                              cgga_pbe=216

  INTEGER, PARAMETER, PUBLIC               :: ke_ol1=301,&
                                              ke_ol2=302,&
                                              ke_llp=303,&
                                              ke_pw86=304,&
                                              ke_pw91=305,&
                                              ke_lc=306,&
                                              ke_t92=307,&
                                              ke_pbe=308,&
                                              ke_revpbe=309,&
                                              ke_tf=310,&
                                              ke_tfw=311

  INTEGER, PARAMETER, PUBLIC               :: pw_orig=209,&
                                              pw_dmc=210,&
                                              pw_vmc=211

  INTEGER, PARAMETER, PUBLIC               :: pz_orig=212,&
                                              pz_dmc=213,&
                                              pz_vmc=214

  INTEGER, PARAMETER, PUBLIC               :: xc_b97_orig=401,&
                                              xc_b97_grimme=402,&
                                              xc_b97_mardirossian=403

END MODULE xc_input_constants
