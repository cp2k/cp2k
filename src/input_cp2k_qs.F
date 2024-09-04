!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the QS section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      07.2024 moved out of input_cp2k_dft [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_qs
   USE bibliography,                    ONLY: &
        Andermatt2016, Brelaz1979, Dewar1977, Dewar1985, Golze2017a, Golze2017b, Iannuzzi2006, &
        Kolafa2004, Krack2000, Kuhne2007, Lippert1997, Lippert1999, Repasky2002, Rocha2006, &
        Schenter2008, Stewart1989, Stewart2007, Thiel1992, VanVoorhis2015, VandeVondele2005a, &
        VandeVondele2006
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              low_print_level
   USE input_constants,                 ONLY: &
        do_ddapc_constraint, do_ddapc_restraint, do_full_density, do_gapw_gcs, do_gapw_gct, &
        do_gapw_log, do_lri_inv, do_lri_inv_auto, do_lri_pseudoinv_diag, do_lri_pseudoinv_svd, &
        do_method_am1, do_method_dftb, do_method_gapw, do_method_gapw_xc, do_method_gpw, &
        do_method_lrigpw, do_method_mndo, do_method_mndod, do_method_ofgpw, do_method_pdg, &
        do_method_pm3, do_method_pm6, do_method_pm6fm, do_method_pnnl, do_method_rigpw, &
        do_method_rm1, do_method_xtb, do_ppl_analytic, do_ppl_grid, do_pwgrid_ns_fullspace, &
        do_pwgrid_ns_halfspace, do_pwgrid_spherical, do_s2_constraint, do_s2_restraint, &
        do_spin_density, gapw_1c_large, gapw_1c_medium, gapw_1c_orb, gapw_1c_small, &
        gapw_1c_very_large, gaussian, numerical, slater, wfi_aspc_nr, wfi_frozen_method_nr, &
        wfi_linear_p_method_nr, wfi_linear_ps_method_nr, wfi_linear_wf_method_nr, &
        wfi_ps_method_nr, wfi_use_guess_method_nr, wfi_use_prev_p_method_nr, &
        wfi_use_prev_rho_r_method_nr, wfi_use_prev_wf_method_nr
   USE input_cp2k_distribution,         ONLY: create_distribution_section
   USE input_cp2k_opt,                  ONLY: create_optimize_dmfet,&
                                              create_optimize_embed,&
                                              create_optimize_lri_basis_section
   USE input_cp2k_scf,                  ONLY: create_cdft_control_section
   USE input_cp2k_se,                   ONLY: create_se_control_section
   USE input_cp2k_tb,                   ONLY: create_dftb_control_section,&
                                              create_xtb_control_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t,&
                                              lchar_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE pw_grids,                        ONLY: do_pw_grid_blocked_false,&
                                              do_pw_grid_blocked_free,&
                                              do_pw_grid_blocked_true
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_qs'

   PUBLIC :: create_qs_section, create_lrigpw_section, create_ddapc_restraint_section

CONTAINS

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

END MODULE input_cp2k_qs
