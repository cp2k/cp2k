!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the xc section of the input
!> \par History
!>      10.2009 moved out of input_cp2k_dft [jgh]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_xc
   USE bibliography,                    ONLY: &
        Becke1988, Becke1997, BeckeRoussel1989, Goedecker1996, Grimme2006, Grimme2010, Grimme2011, &
        Heyd2004, Kruse2012, Lee1988, Lehtola2018, Marques2012, Ortiz1994, Paziani2006, &
        Perdew1981, Perdew1996, Perdew2008, Proynov2007, Tao2003, Toulouse2004, Tran2013, &
        Vosko1980, Wellendorff2012, Zhang1998
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              high_print_level
   USE input_constants,                 ONLY: &
        do_adiabatic_hybrid_mcy3, do_adiabatic_model_pade, gaussian, slater, vdw_nl_drsll, &
        vdw_nl_lmkll, vdw_nl_rvv10, vdw_pairpot_dftd2, vdw_pairpot_dftd3, vdw_pairpot_dftd3bj, &
        xc_funct_b3lyp, xc_funct_beefvdw, xc_funct_blyp, xc_funct_bp, xc_funct_hcth120, &
        xc_funct_no_shortcut, xc_funct_olyp, xc_funct_pade, xc_funct_pbe, xc_funct_pbe0, &
        xc_funct_tpss, xc_none, xc_pot_energy_none, xc_pot_energy_sum_eigenvalues, &
        xc_pot_energy_xc_functional, xc_vdw_fun_none, xc_vdw_fun_nonloc, xc_vdw_fun_pairpot
   USE input_cp2k_hfx,                  ONLY: create_hfx_section
   USE input_cp2k_mp2,                  ONLY: create_mp2_section
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
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
   USE xc_input_constants,              ONLY: &
        c_pw92, c_pw92dmc, c_pw92vmc, c_pz, c_pzdmc, c_pzvmc, do_vwn3, do_vwn5, ke_lc, ke_llp, &
        ke_ol1, ke_ol2, ke_pbe, ke_pw86, ke_pw91, ke_t92, pz_orig, xalpha, xc_b97_3c, &
        xc_b97_grimme, xc_b97_mardirossian, xc_b97_orig, xc_debug_new_routine, xc_deriv_collocate, &
        xc_deriv_nn10_smooth, xc_deriv_nn50_smooth, xc_deriv_pw, xc_deriv_spline2, &
        xc_deriv_spline2_smooth, xc_deriv_spline3, xc_deriv_spline3_smooth, xc_new_f_routine, &
        xc_old_f_routine, xc_pbe_orig, xc_pbe_rev, xc_pbe_sol, xc_rho_nn10, xc_rho_nn50, &
        xc_rho_no_smooth, xc_rho_spline2_smooth, xc_rho_spline3_smooth, xc_test_lsd_f_routine, &
        xgga_b88x, xgga_ev93, xgga_opt, xgga_pbex, xgga_pw86, xgga_pw91, xgga_revpbe
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_xc'

   PUBLIC :: create_xc_section, create_xc_fun_section

CONTAINS

! **************************************************************************************************
!> \brief creates the structure of the section needed to select the xc functional
!> \param section the section that will be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_xc_fun_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xc_functional", &
                          description="The xc functional to use", &
                          n_keywords=0, n_subsections=4, repeats=.FALSE., &
                          citations=(/Ortiz1994, Becke1988, Perdew1996, Zhang1998, Lee1988, &
                                      Heyd2004, Vosko1980, Goedecker1996, Perdew1981, &
                                      Tao2003, Wellendorff2012/))

      NULLIFY (subsection, keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
         description="Shortcut for the most common functional combinations.", &
         usage="&xc_functional BLYP", &
         enum_c_vals=s2a("B3LYP", "PBE0", "BLYP", "BP", "PADE", "LDA", "PBE", &
                         "TPSS", "HCTH120", "OLYP", "BEEFVDW", "NO_SHORTCUT", "NONE"), &
         enum_i_vals=(/xc_funct_b3lyp, xc_funct_pbe0, xc_funct_blyp, xc_funct_bp, xc_funct_pade, xc_funct_pade, xc_funct_pbe, &
                       xc_funct_tpss, xc_funct_hcth120, xc_funct_olyp, xc_funct_beefvdw, xc_funct_no_shortcut, xc_none/), &
         enum_desc=s2a("B3LYP", &
                       "PBE0 (see note in section XC/XC_FUNCTIONAL/PBE)", &
                       "BLYP", "BP", "PADE", "Alias for PADE", &
                       "PBE (see note in section XC/XC_FUNCTIONAL/PBE)", &
                       "TPSS (not available with LSD, use LIBXC version instead)", "HCTH120", "OLYP", &
                       "BEEFVDW", "NO_SHORTCUT", "NONE"), &
         default_i_val=xc_funct_no_shortcut, &
         lone_keyword_i_val=xc_funct_no_shortcut)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="BECKE88", &
                          description="Uses the Becke 88 exchange functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Becke1988/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="scale_x", &
         description="scales the exchange part of the functional", &
         default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="LYP_ADIABATIC", &
                          description="Uses the LYP correlation functional in an adiabatic fashion", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Lee1988/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="LAMBDA", &
                          description="Defines the parameter of the adiabatic curve.", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="BECKE88_LR_ADIABATIC", &
                          description="Uses the Becke 88 longrange exchange functional in an adiabatic fashion", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Becke1988/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OMEGA", &
                          description="Potential parameter in erf(omega*r)/r", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="LAMBDA", &
                          description="Defines the parameter of the adiabatic curve", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="BECKE88_LR", &
                          description="Uses the Becke 88 longrange exchange functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Becke1988/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="OMEGA", &
                          description="Potential parameter in erf(omega*r)/r", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="LYP", &
                          description="Uses the LYP functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Lee1988/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PADE", &
                          description="Uses the PADE functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Goedecker1996/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="HCTH", &
                          description="Uses the HCTH class of functionals", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETER_SET", &
                          description="Which version of the parameters should be used", &
                          usage="PARAMETER_SET 407", &
                          enum_c_vals=(/"93 ", "120", "147", "407", "HLE"/), &
                          enum_i_vals=(/93, 120, 147, 407, 408/), &
                          default_i_val=120)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="OPTX", &
                          description="Uses the OPTX functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="a1", &
                          description="OPTX a1 coefficient", &
                          default_r_val=1.05151_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="a2", &
                          description="OPTX a2 coefficient", &
                          default_r_val=1.43169_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="gamma", &
                          description="OPTX gamma coefficient", &
                          default_r_val=0.006_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_libxc_section(subsection, "LIBXC", &
                                "Uses functionals from LIBXC, see also "// &
                                "https://www.tddft.org/programs/libxc/functionals/", &
                                "FUNCTIONAL GGA_X_PBE")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_libxc_section(subsection, "KE_LIBXC", &
                                "To be used for KG runs. Uses kinetic energy functionals from LIBXC, "// &
                                "https://www.tddft.org/programs/libxc/functionals/", &
                                "FUNCTIONAL GGA_K_LLP")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="CS1", &
                          description="Uses the CS1 functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="XGGA", &
                          description="Uses one of the XGGA functionals (optimized versions of "// &
                          "some of these functionals might be available outside this section).", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL", &
                          description="Which one of the XGGA functionals should be used", &
                          usage="FUNCTIONAL PW86X", &
                          enum_c_vals=(/ &
                          "BECKE88X", &
                          "PW86X   ", &
                          "PW91X   ", &
                          "PBEX    ", &
                          "REV_PBEX", &
                          "OPTX    ", &
                          "EV93    "/), &
                          enum_i_vals=(/xgga_b88x, xgga_pw86, xgga_pw91, xgga_pbex, xgga_revpbe, xgga_opt, xgga_ev93/), &
                          default_i_val=xgga_b88x)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="KE_GGA", &
                          description="Uses one of the KE_GGA functionals (optimized versions of "// &
                          "some of these functionals might be available outside this section). "// &
                          "These functionals are needed for the computation of the kinetic "// &
                          "energy in the Kim-Gordon method.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL", &
                          description="Which one of the KE_GGA functionals should be used", &
                          usage="FUNCTIONAL (OL1|OL2|LLP|PW86|PW91|LC|T92|PBE)", &
                          enum_c_vals=(/"OL1 ", "OL2 ", "LLP ", "PW86", "PW91", "LC  ", "T92 ", "PBE "/), &
                          enum_i_vals=(/ke_ol1, ke_ol2, ke_llp, ke_pw86, ke_pw91, ke_lc, ke_t92, ke_pbe/), &
                          enum_desc=s2a("Uses first Ou-Yang and Levy functional, currently not producing correct results", &
                                        "Uses second Ou-Yang and Levy functional, currently not producing correct results", &
                                        "Uses Lee, Lee, and Parr functional", &
                                        "Uses Perdew and Wang's 1986 functional", &
                                        "Uses Perdew and Wang's 1991 functional", &
                                        "Uses Lembarki and Chermette functional", &
                                        "Uses Thakkar functional", &
                                        "Uses the 1996 functional of Perdew, Burke and Ernzerhof"), &
                          default_i_val=ke_llp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="P86C", &
                          description="Uses the P86C functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PW92", &
                          description="Uses the PerdewWang correlation functional.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          description="Scaling of the energy functional", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETRIZATION", &
                          description="Which one of parametrizations should be used", &
                          usage="PARAMETRIZATION DMC", &
                          enum_c_vals=(/ &
                          "ORIGINAL", &
                          "DMC     ", &
                          "VMC     "/), &
                          enum_i_vals=(/c_pw92, c_pw92dmc, c_pw92vmc/), &
                          default_i_val=c_pw92)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PZ81", &
                          description="Uses the PZ functional.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE., &
                          citations=(/Perdew1981, Ortiz1994/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETRIZATION", &
                          description="Which one of parametrizations should be used", &
                          usage="PARAMETRIZATION DMC", &
                          enum_c_vals=(/ &
                          "ORIGINAL", &
                          "DMC     ", &
                          "VMC     "/), &
                          enum_i_vals=(/c_pz, c_pzdmc, c_pzvmc/), &
                          default_i_val=pz_orig)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="TFW", &
                          description="Uses the TFW functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="TF", &
                          description="Uses the TF functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="VWN", &
                          description="Uses the VWN functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Vosko1980/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL_TYPE", &
                          description="Which version of the VWN functional should be used", &
                          usage="FUNCTIONAL_TYPE VWN5", &
                          enum_c_vals=s2a("VWN5", "VWN3"), &
                          enum_i_vals=(/do_vwn5, do_vwn3/), &
                          enum_desc=s2a("This is the recommended (correct) version of the VWN functional", &
                                        "This version is the default in Gaussian, but not recommended."// &
                                        "Notice that it is also employed in Gaussian's default version of B3LYP"), &
                          default_i_val=do_vwn5)

      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="XALPHA", &
                          description="Uses the XALPHA (SLATER) functional.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="XA", &
                          description="Value of the xa parameter (this does not change the exponent, "// &
                          "just the mixing)", &
                          usage="XA 0.7", default_r_val=2._dp/3._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="TPSS", &
                          description="Uses the TPSS functional. Note, that there is no LSD version available. "// &
                          "In such cases, use the LIBXC version instead.", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Tao2003/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PBE", &
                          description="Uses the PBE functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Perdew1996, Zhang1998, Perdew2008/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="parametrization", &
                          description="switches between the different "// &
                          "parametrizations of the functional. "// &
                          "Note: Beta parameters used have only 5 significant digits, "// &
                          "as published. For higher precision and program comparison "// &
                          "use section XC/XC_FUNCTIONAL/LIBXC.", &
                          enum_i_vals=(/xc_pbe_orig, xc_pbe_rev, xc_pbe_sol/), &
                          enum_c_vals=(/"ORIG  ", "revPBE", "PBEsol"/), &
                          enum_desc=(/"original PBE                        ", &
                                      "revised PBE (revPBE)                ", &
                                      "PBE for solids and surfaces (PBEsol)"/), &
                          default_i_val=xc_pbe_orig)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="XWPBE", &
                          description="Uses the short range PBE functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Heyd2004/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x0", &
                          description="scales the exchange part of the original hole PBE-functional", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="omega", &
                          description="screening parameter", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="BECKE97", &
                          description="Uses the Becke 97 exchange correlation functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Becke1997, Grimme2006/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                description="scales the exchange part of the functional, if -1 the default for the given parametrization is used", &
                          default_r_val=-1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="parametrization", &
                          description="switches between the B97 and Grimme parametrization ", &
                          enum_i_vals=(/xc_b97_orig, xc_b97_grimme, xc_b97_grimme, xc_b97_mardirossian, xc_b97_3c/), &
                          enum_c_vals=(/"ORIG      ", "B97GRIMME ", "B97_GRIMME", "wB97X-V   ", "B97-3c    "/), &
                          default_i_val=xc_b97_orig)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="BECKE_ROUSSEL", &
                          description="Becke Roussel exchange hole model. Can be used"// &
                          "as long range correction with a truncated coulomb potential", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/BeckeRoussel1989, Proynov2007/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="Defines the cutoff radius for the truncation. "// &
                          "If put to zero, the standard full range potential will be used", &
                          usage="CUTOFF_RADIUS 2.0", default_r_val=0.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GAMMA", &
                          description="Parameter in the exchange hole. "// &
                          "Usually this is put to 1.0 or 0.8", &
                          usage="GAMMA 0.8", default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="LDA_HOLE_T_C_LR", &
                          description="LDA exchange hole model in truncated coulomb potential", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE_X", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="Defines cutoff for lower integration boundary", &
                          default_r_val=0.0_dp, unit_str="angstrom")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PBE_HOLE_T_C_LR", &
                          description="PBE exchange hole model in trucanted coulomb potential", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE_X", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="Defines cutoff for lower integration boundary", &
                          default_r_val=1.0_dp, unit_str="angstrom")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="GV09", &
                          description="Combination of three different exchange hole models", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE_X", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="Defines cutoff for lower integration boundary", &
                          default_r_val=0.0_dp, unit_str="angstrom")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GAMMA", &
                          description="Parameter for Becke Roussel hole", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="BEEF", & !rk: BEEF Exchange
                          description="Uses the BEEFvdW exchange functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Wellendorff2012/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="SRLDA", &
                          description="Uses the short-range LDA functional", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Paziani2006, Toulouse2004/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the functional", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_x", &
                          description="scales the exchange part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="scale_c", &
                          description="scales the correlation part of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="omega", &
                          description="provide the range-separation parameter of the functional", &
                          default_r_val=1._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETRIZATION", &
                          description="Which one parametrizations of the underlying PW92 functional should be used", &
                          usage="PARAMETRIZATION DMC", &
                          enum_c_vals=(/ &
                          "ORIGINAL", &
                          "DMC     ", &
                          "VMC     "/), &
                          enum_i_vals=(/c_pw92, c_pw92dmc, c_pw92vmc/), &
                          default_i_val=c_pw92)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xc_fun_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \param name ...
!> \param description ...
!> \param usage ...
! **************************************************************************************************
   SUBROUTINE create_libxc_section(section, name, description, usage)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: name, description, usage

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CPASSERT(name == "LIBXC" .OR. name == "KE_LIBXC")
      NULLIFY (keyword)
      CALL section_create(section, __LOCATION__, name, description, &
                          n_keywords=3, n_subsections=0, repeats=.TRUE., &
                          citations=(/Marques2012, Lehtola2018/))
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates functionals from libxc library", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FUNCTIONAL", &
                          description="name of the functional, see also "// &
                          "https://www.tddft.org/programs/libxc/functionals/ "// &
                          "The precise list of available functionals depends on "// &
                          "the version of the linked libxc.", &
                          usage=usage, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          description="scaling factors for the functional", &
                          usage="SCALE 0.8", type_of_var=real_t, &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETERS", &
                          description="list of external parameters for the functional", &
                          usage="PARAMETERS 0.3 [0.5 [1.4]]", type_of_var=real_t, &
                          default_r_vals=(/HUGE(0.0_dp)/), n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE

! **************************************************************************************************
!> \brief creates the structure of the section needed to select an xc potential
!> \param section the section that will be created
!> \author thomas chassaing
! **************************************************************************************************
   SUBROUTINE create_xc_potential_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xc_potential", &
                          description="The xc potential to use (CAREFUL: xc potential here refers "// &
                          "to potentials that are not derived from an xc functional, but rather are "// &
                          "modelled directly. Therefore there is no consistent xc energy available. "// &
                          "To still get an energy expression, see ENERGY below", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE.)

      NULLIFY (subsection, keyword)
      CALL section_create(subsection, __LOCATION__, name="SAOP", &
                          description="Uses the SAOP potential", &
                          n_keywords=3, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="ALPHA", &
                          description="Value of the alpha parameter (default = 1.19).", &
                          usage="ALPHA 1.19", default_r_val=1.19_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="BETA", &
                          description="Value of the beta parameter (default = 0.01).", &
                          usage="BETA 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="K_RHO", &
                          description="Value of the K_rho parameter (default = 0.42).", &
                          usage="ALPHA 0.42", default_r_val=0.42_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY", &
                          description="How to determine the total energy.", &
                          usage="ENERGY [NONE,XC_FUNCTIONAL,SUM_EIGENVALUES", &
                          enum_c_vals=s2a("NONE", "XC_FUNCTIONAL", "FUNCTIONAL", "SUM_EIGENVALUES", "SOE"), &
                          enum_i_vals=(/xc_pot_energy_none, &
                                        xc_pot_energy_xc_functional, &
                                        xc_pot_energy_xc_functional, &
                                        xc_pot_energy_sum_eigenvalues, &
                                        xc_pot_energy_sum_eigenvalues/), &
                          default_i_val=xc_pot_energy_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xc_potential_section

! **************************************************************************************************
!> \brief creates the structure of the section needed for vdW potentials
!> \param section the section that will be created
!> \author jgh
! **************************************************************************************************
   SUBROUTINE create_vdw_potential_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="vdw_potential", &
                          description="This section combines all possible additional dispersion "// &
                          "corrections to the normal XC functionals. This can be more functionals "// &
                          "or simple empirical pair potentials. ", &
                          citations=(/grimme2006, Tran2013/), &
                          n_keywords=1, n_subsections=1, repeats=.FALSE.)

      NULLIFY (subsection, keyword)
      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_TYPE", &
                          variants=s2a("DISPERSION_FUNCTIONAL"), &
                          description="Type of dispersion/vdW functional or potential to use", &
                          usage="POTENTIAL_TYPE (NONE|PAIR_POTENTIAL|NON_LOCAL)", &
                          enum_c_vals=s2a("NONE", "PAIR_POTENTIAL", "NON_LOCAL"), &
                          enum_i_vals=(/xc_vdw_fun_none, xc_vdw_fun_pairpot, xc_vdw_fun_nonloc/), &
                          enum_desc=s2a("No dispersion/van der Waals functional", &
                                        "Pair potential van der Waals density functional", &
                                        "Nonlocal van der Waals density functional"), &
                          default_i_val=xc_vdw_fun_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="PAIR_POTENTIAL", &
                          description="Information on the pair potential to calculate dispersion", &
                          n_keywords=5, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="R_CUTOFF", &
                          description="Range of potential. The cutoff will be 2 times this value", &
                          usage="R_CUTOFF 24.0", default_r_val=20.0_dp, &
                          unit_str="angstrom")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="TYPE", &
                          description="Type of potential", &
                          citations=(/grimme2006, grimme2010, grimme2011/), &
                          usage="TYPE (DFTD2|DFTD3|DFTD3(BJ))", &
                          enum_c_vals=s2a("DFTD2", "DFTD3", "DFTD3(BJ)"), &
                          enum_i_vals=(/vdw_pairpot_dftd2, vdw_pairpot_dftd3, vdw_pairpot_dftd3bj/), &
                          enum_desc=s2a("Grimme D2 method", &
                                        "Grimme D3 method (zero damping)", &
                                        "Grimme D3 method (Becke-Johnson damping)"), &
                          default_i_val=vdw_pairpot_dftd3)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETER_FILE_NAME", &
                          description="Name of the parameter file, may include a path", &
                          usage="PARAMETER_FILE_NAME <FILENAME>", &
                          default_lc_val="DISPERSION_PARAMETERS")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_FUNCTIONAL", &
                          description="Use parameters for this specific density functional. " &
                          //"For available D3 and D3(BJ) parameters see: " &
                          //"https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/dft-d3", &
                          usage="REFERENCE_FUNCTIONAL <functional>", &
                          type_of_var=char_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALING", &
                          description="XC Functional dependent scaling parameter, if set to zero CP2K attempts"// &
                          " to guess the xc functional that is in use and sets the associated scaling parameter.", &
                          usage="SCALING 0.2", default_r_val=0._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EXP_PRE", &
                          description="Prefactor in exponential damping factor (DFT-D2 potential)", &
                          usage="EXP_PRE 20.", default_r_val=20._dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_CN", &
                          description="Cutoff value for coordination number function (DFT-D3 method)", &
                          usage="EPS_CN 1.e-6_dp", default_r_val=1.e-6_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="D3_SCALING", &
                          description="XC Functional dependent scaling parameters (s6,sr6,s8) for the DFT-D3 method,"// &
                          " if set to zero CP2K attempts"// &
                          " to guess the xc functional from REFERENCE_FUNCTIONAL and sets the associated scaling parameter.", &
                          usage="D3_SCALING 1.0 1.0 1.0", n_var=3, &
                          default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_SCALING", &
                          description="XC Functional dependent scaling parameters (s6,a1,s8,a2) for the DFT-D3(BJ) method,"// &
                          " if set to zero CP2K attempts"// &
                          " to guess the xc functional from REFERENCE_FUNCTIONAL and sets the associated scaling parameter.", &
                          usage="D3BJ_SCALING 1.0 1.0 1.0 1.0", n_var=4, &
                          default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CALCULATE_C9_TERM", &
                          description="Calculate C9 terms in DFT-D3 model", &
                          usage="CALCULATE_C9_TERM", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_C9_TERM", &
                          description="Calculate C9 terms in DFT-D3 model using reference coordination numbers", &
                          usage="REFERENCE_C9_TERM", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="LONG_RANGE_CORRECTION", &
                          description="Calculate a long range correction to the DFT-D3 model", &
                          usage="LONG_RANGE_CORRECTION", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SHORT_RANGE_CORRECTION", &
                          description="Calculate a short-range bond correction to the DFT-D3 model", &
                          usage="SHORT_RANGE_CORRECTION", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SHORT_RANGE_CORRECTION_PARAMETERS", &
                          description="Parameters for the short-range bond correction to the DFT-D3 model. "// &
                          " s*(za*zb)^t1*EXP(-g*dr*r0ab^t2), parameters: s, g, t1, t2 "// &
                          " Defaults: s=0.08, g=10.0, t1=0.5, t2=-1.0 ", &
                          usage="SHORT_RANGE_CORRECTION_PARAMETRS", default_r_vals=(/0.08_dp, 10.0_dp, 0.5_dp, -1.0_dp/), &
                          n_var=4, type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      ! KG molecular corrections
      CALL keyword_create(keyword, __LOCATION__, name="MOLECULE_CORRECTION", &
                          description="Calculate a intermolecular correction to the DFT-D3 model", &
                          usage="MOLECULE_CORRECTION", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="MOLECULE_CORRECTION_C8", &
                          description="Calculate a intermolecular correction to the C8 term in the DFT-D3 model", &
                          usage="MOLECULE_CORRECTION_C8  1.0 ", default_r_val=0.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="VERBOSE_OUTPUT", &
                          description="Extensive output for the DFT-D2 and DFT-D3 models", &
                          usage="VERBOSE_OUTPUT", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Set coordination numbers by atom numbers
      CALL keyword_create(keyword, __LOCATION__, name="D3_EXCLUDE_KIND", &
                          description="Specifies the atomic kinds excluded in the DFT-D3 calculation.", &
                          usage="D3_EXCLUDE_KIND kind1 kind2 ... ", repeats=.FALSE., &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Ignore selected pair interactins
      CALL keyword_create(keyword, __LOCATION__, name="D3_EXCLUDE_KIND_PAIR", &
                          description="Specifies the atomic kinds for interactions excluded from the DFT-D3 calculation.", &
                          usage="D3_EXCLUDE_KIND_PAIR kind1 kind2 ", repeats=.TRUE., &
                          n_var=2, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Set coordination numbers by atom kinds
      CALL keyword_create(keyword, __LOCATION__, name="KIND_COORDINATION_NUMBERS", &
                          description="Specifies the coordination number for a kind for the C9 term in DFT-D3.", &
                          usage="KIND_COORDINATION_NUMBERS CN kind ", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      ! Set coordination numbers by atom numbers
      CALL keyword_create(keyword, __LOCATION__, name="ATOM_COORDINATION_NUMBERS", &
                          description="Specifies the coordination number of a set of atoms for the C9 term in DFT-D3.", &
                          usage="ATOM_COORDINATION_NUMBERS CN atom1 atom2 ... ", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! parameter specification atom by atom
      CALL keyword_create(keyword, __LOCATION__, name="ATOMPARM", &
                          description="Specifies parameters for atom types (in atomic units). If "// &
                          "not provided default parameters are used (DFT-D2).", &
                          usage="ATOMPARM <ELEMENT> <C6_parameter> <vdw_radii>", &
                          repeats=.TRUE., n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PRINT_DFTD", &
                                       description="Controls the printing of some info about DFTD contributions", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! nonlocal section
      NULLIFY (subsection, keyword)
      CALL section_create(subsection, __LOCATION__, name="NON_LOCAL", &
                          description="Information on the non local part of dispersion functionals."// &
                          "Correct functionals require a corresponding setting of XC_FUNCTIONAL.", &
                          n_keywords=0, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE", &
                          description="Type of functional (the corresponding kernel data file should be selected)."// &
                          "Allows for common forms such as vdW-DF, vdW-DF2, optB88-vdW, rVV10.", &
                          usage="TYPE DRSLL", &
                          enum_c_vals=s2a("DRSLL", "LMKLL", "RVV10"), &
                          enum_i_vals=(/vdw_nl_DRSLL, vdw_nl_LMKLL, vdw_nl_RVV10/), &
                          enum_desc=s2a("Dion-Rydberg-Schroeder-Langreth-Lundqvist nonlocal van der Waals density functional", &
                                        "Lee-Murray-Kong-Lundqvist-Langreth nonlocal van der Waals density functional", &
                                        "Revised Vydrov-van Voorhis nonlocal van der Waals density functional"), &
                          citations=(/Tran2013/), &
                          default_i_val=vdw_nl_DRSLL)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="VERBOSE_OUTPUT", &
                          description="Extensive output for non local functionals", &
                          usage="VERBOSE_OUTPUT", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="KERNEL_FILE_NAME", &
                          description="Name of the kernel data file, may include a path."// &
                          "vdW_kernel_table.dat is for DRSLL and LMKLL and"// &
                          "rVV10_kernel_table.dat is for rVV10.", &
                          usage="KERNEL_FILE_NAME <FILENAME>", &
                          default_lc_val="vdW_kernel_table.dat")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="The cutoff of the FFT grid used in the calculation "// &
                          "of the nonlocal vdW functional [Ry].", &
                          usage="CUTOFF 300", &
                          default_r_val=-1._dp, unit_str="Ry")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETERS", &
                          description="Parameters b and C of the rVV10 functional", &
                          usage="PARAMETERS 6.3 0.0093", &
                          type_of_var=real_t, default_r_vals=(/6.3_dp, 0.0093_dp/), n_var=2)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          description="Scales the energy contribution of the rVV10 functional", &
                          usage="SCALE 1.0", &
                          type_of_var=real_t, default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_vdw_potential_section

! **************************************************************************************************
!> \brief creates the structure of the section needed for gCP potentials
!> \param section the section that will be created
!> \author jgh
! **************************************************************************************************
   SUBROUTINE create_gcp_potential_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="gcp_potential", &
                          description="This section combines geometrical counterpoise potentials."// &
                          " This is a simple empirical pair potential to correct for BSSE. ", &
                          citations=(/Kruse2012/), &
                          n_keywords=1, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAMETER_FILE_NAME", &
                          description="Name of the parameter file, may include a path", &
                          usage="PARAMETER_FILE_NAME <FILENAME>", &
                          default_lc_val="---")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GLOBAL_PARAMETERS", &
                          description="Global parameters of the gCP method."// &
                          " Parameters are sigma, alpha, beta, eta from the original paper.", &
                          usage="GLOBAL_PARAMETERS 1.0 1.0 1.0 1.0", n_var=4, &
                          default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_ENERGY", &
                          description="Specify the delta energy [Hartree] term for an atom kind", &
                          usage="DELTA_ENERGY  type value", &
                          type_of_var=char_t, repeats=.TRUE., n_var=-1, default_c_vals=(/"XX ", "0.0"/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VERBOSE", &
                          description="Verbose output for gCP calculation", &
                          usage="VERBOSE logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_gcp_potential_section

! **************************************************************************************************
!> \brief creates the input section for the xc part
!> \param section the section to create
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_xc_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xc", &
                          description="Parameters needed for the calculation of the eXchange and Correlation potential", &
                          n_keywords=5, n_subsections=2, repeats=.FALSE.)

      NULLIFY (subsection, keyword)

      CALL keyword_create(keyword, __LOCATION__, name="density_cutoff", &
                          description="The cutoff on the density used by the xc calculation", &
                          usage="density_cutoff 1.e-11", default_r_val=1.0e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="gradient_cutoff", &
                          description="The cutoff on the gradient of the density used by the "// &
                          "xc calculation", &
                          usage="gradient_cutoff 1.e-11", default_r_val=1.0e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_SMOOTH_CUTOFF_RANGE", &
                          description="Parameter for the smoothing procedure in"// &
                          "xc calculation", &
                          usage="gradient_cutoff {real}", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="tau_cutoff", &
                          description="The cutoff on tau used by the xc calculation", &
                          usage="tau_cutoff 1.e-11", default_r_val=1.0e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="FUNCTIONAL_ROUTINE", &
         description="Select the code for xc calculation", &
         usage="FUNCTIONAL_ROUTINE (NEW|OLD|TEST_LSD|DEBUG)", &
         default_i_val=xc_new_f_routine, &
         enum_c_vals=s2a("NEW", "OLD", "TEST_LSD", "DEBUG"), &
         enum_i_vals=(/xc_new_f_routine, xc_old_f_routine, &
                       xc_test_lsd_f_routine, xc_debug_new_routine/), &
         enum_desc=s2a("Use new code for exchange-correlation functional calculation", &
                       "Use old code for exchange-correlation functional calculation", &
                       "Use test local-spin-density approximation code for exchange-correlation functional calculation", &
                       "Use debug new code for exchange-correlation functional calculation"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="xc_grid", & !FM to do
                          description="The xc parameters used when calculating the xc on the grid", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="xc_smooth_rho", &
                          description="The density smoothing used for the xc calculation", &
                          usage="xc_smooth_rho nn10", default_i_val=xc_rho_no_smooth, &
                          enum_c_vals=s2a("NONE", "NN50", "NN10", "SPLINE2", "NN6", "SPLINE3", "NN4"), &
                          enum_i_vals=(/xc_rho_no_smooth, xc_rho_nn50, xc_rho_nn10, &
                                        xc_rho_spline2_smooth, xc_rho_spline2_smooth, xc_rho_spline3_smooth, &
                                        xc_rho_spline3_smooth/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="xc_deriv", &
                          description="The method used to compute the derivatives", &
                          usage="xc_deriv NN10_SMOOTH", default_i_val=xc_deriv_pw, &
                          enum_c_vals=s2a("PW", "SPLINE3", "SPLINE2", "NN50_SMOOTH", "NN10_SMOOTH", &
                                          "SPLINE2_SMOOTH", "NN6_SMOOTH", "SPLINE3_SMOOTH", "NN4_SMOOTH", "COLLOCATE"), &
                          enum_i_vals=(/xc_deriv_pw, xc_deriv_spline3, xc_deriv_spline2, &
                                        xc_deriv_nn50_smooth, xc_deriv_nn10_smooth, xc_deriv_spline2_smooth, &
                                        xc_deriv_spline2_smooth, xc_deriv_spline3_smooth, xc_deriv_spline3_smooth, &
                                        xc_deriv_collocate/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="use_finer_grid", &
                          description="Uses a finer grid only to calculate the xc", &
                          usage="use_finer_grid", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="2ND_DERIV_ANALYTICAL", &
                          description="Use analytical formulas or finite differences for 2nd derivatives of XC", &
                          usage="2ND_DERIV_ANALYTICAL logical", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="3RD_DERIV_ANALYTICAL", &
                          description="Use analytical formulas or finite differences for 3rd derivatives of XC", &
                          usage="3RD_DERIV_ANALYTICAL logical", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xc_fun_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_hfx_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mp2_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_adiabatic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xc_potential_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_vdw_potential_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_gcp_potential_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xc_section

! **************************************************************************************************
!> \brief creates the section for adiabatic hybrid functionals
!> \param section ...
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE create_adiabatic_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ADIABATIC_RESCALING", &
                          description="Parameters for self interaction corrected hybrid functionals", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="FUNCTIONAL_TYPE", &
         description="Which Hybrid functional should be used. "// &
         "(Has to be consistent with the definitions in XC and HF).", &
         usage="FUNCTIONAL_TYPE MCY3", &
         enum_c_vals=s2a("MCY3"), &
         enum_i_vals=(/do_adiabatic_hybrid_mcy3/), &
         enum_desc=s2a("Use MCY3 hybrid functional"), &
         default_i_val=do_adiabatic_hybrid_mcy3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="LAMBDA", &
         description="The point to be used along the adiabatic curve (0 &lt; &lambda; &lt; 1)", &
         usage="LAMBDA 0.71", &
         default_r_val=0.71_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="OMEGA", &
         description="Long-range parameter", &
         usage="OMEGA 0.2", &
         default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="FUNCTIONAL_MODEL", &
         description="Which model for the coupling constant integration should be used. ", &
         usage="FUNCTIONAL_MODEL PADE", &
         enum_c_vals=s2a("PADE"), &
         enum_i_vals=(/do_adiabatic_model_pade/), &
         enum_desc=s2a("Use pade model: W(lambda)=a+(b*lambda)/(1+c*lambda)"), &
         default_i_val=do_adiabatic_model_pade)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_adiabatic_section

END MODULE input_cp2k_xc
