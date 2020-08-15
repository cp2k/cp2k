!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      01.2013 moved out of input_cp2k_dft [MI]
!> \author MI
! **************************************************************************************************
MODULE input_cp2k_properties_dft

   USE bibliography,                    ONLY: Futera2017,&
                                              Iannuzzi2005,&
                                              Kondov2007,&
                                              Luber2014,&
                                              Putrino2000,&
                                              Putrino2002,&
                                              Sebastiani2001,&
                                              Weber2009
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level,&
                                              silent_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        current_gauge_atom, current_gauge_r, current_gauge_r_and_step_func, &
        current_orb_center_atom, current_orb_center_box, current_orb_center_common, &
        current_orb_center_wannier, do_et_ddapc, do_full_density, do_no_et, do_spin_density, &
        gto_cartesian, gto_spherical, oe_gllb, oe_lb, oe_none, oe_saop, oe_shift, &
        ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
        tddfpt_dipole_berry, tddfpt_dipole_length, tddfpt_dipole_velocity, tddfpt_kernel_full, &
        tddfpt_kernel_stda, use_mom_ref_coac, use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero
   USE input_cp2k_atprop,               ONLY: create_atprop_section
   USE input_cp2k_dft,                  ONLY: create_ddapc_restraint_section,&
                                              create_interp_section,&
                                              create_mgrid_section
   USE input_cp2k_loc,                  ONLY: create_localize_section
   USE input_cp2k_resp,                 ONLY: create_resp_section
   USE input_cp2k_xc,                   ONLY: create_xc_section
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
                                              logical_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_properties_dft'

   PUBLIC :: create_properties_section

CONTAINS

! **************************************************************************************************
!> \brief Create the PROPERTIES section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_properties_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PROPERTIES", &
                          description="This section is used to set up the PROPERTIES calculation.", &
                          n_keywords=0, n_subsections=6, repeats=.FALSE.)

      NULLIFY (subsection, keyword)

      CALL create_linres_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_et_coupling_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_resp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_atprop_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL cp_print_key_section_create(subsection, __LOCATION__, name="FIT_CHARGE", &
                                       description="This section is used to print the density derived atomic point charges."// &
                                       "The fit of the charges is controlled through the DENSITY_FITTING section", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="TYPE_OF_DENSITY", &
                          description="Specifies the type of density used for the fitting", &
                          usage="TYPE_OF_DENSITY (FULL|SPIN)", &
                          enum_c_vals=s2a("FULL", "SPIN"), &
                          enum_i_vals=(/do_full_density, do_spin_density/), &
                          enum_desc=s2a("Full density", "Spin density"), &
                          default_i_val=do_full_density)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_tddfpt2_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_properties_section

! **************************************************************************************************
!> \brief creates the input structure used to activate
!>      a linear response calculation
!>      Available properties : none
!> \param section the section to create
!> \author MI
! **************************************************************************************************
   SUBROUTINE create_linres_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, subsection, print_key)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="linres", &
                          description="The linear response is used to calculate one of the "// &
                          " following properties: nmr, epr, raman, ... ", &
                          n_keywords=5, n_subsections=2, repeats=.FALSE., &
                          citations=(/Putrino2000/))

      CALL keyword_create(keyword, __LOCATION__, name="EPS", &
                          description="target accuracy for the convergence of the conjugate gradient.", &
                          usage="EPS 1.e-6", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Filter threshold for response density matrix.", &
                          usage="EPS 1.e-8", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of conjugate gradient iteration to be performed for one optimization.", &
                          usage="MAX_ITER 200", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_EVERY", &
                          description="Restart the conjugate gradient after the specified number of iterations.", &
                          usage="RESTART_EVERY 200", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with all minimization schemes. "// &
         "They differ in effectiveness, cost of construction, cost of application. "// &
         "Properly preconditioned minimization can be orders of magnitude faster than doing nothing.", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=ot_precond_none, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", &
                         "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "// &
                       "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "// &
                       "This preconditioner is recommended for almost all systems, except very large systems where "// &
                       "make_preconditioner would dominate the total computational cost.", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, "// &
                       "might be somewhat less robust. Recommended for large systems.", &
                       "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ", &
                       "Cholesky inversion of S and T, fast construction, robust, and relatively good, "// &
                       "use for very large systems.", &
                       "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_full_single, &
                       ot_precond_full_kinetic, ot_precond_s_inverse, ot_precond_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_GAP", &
                          description="Energy gap estimate [a.u.] for preconditioning", &
                          usage="ENERGY_GAP 0.1", &
                          default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the response calculation if the restart file exists", &
                          usage="RESTART", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          variants=(/"RESTART_FILE_NAME"/), &
                          description="Root of the file names where to read the response functions from"// &
                          "which to restart the calculation of the linear response", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_current_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_nmr_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_spin_spin_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_epr_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_polarizability_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="printing of information during the linear response calculation", &
                          repeats=.FALSE.)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "program_run_info", &
         description="Controls the printing of basic iteration information during the LINRES calculation", &
         print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of restart file of the response wavefunction."// &
                                       "For each set of response functions, i.e. for each perturbation,"// &
                                       "one different restart file is dumped. These restart files should be"// &
                                       "employed only to restart the same type of LINRES calculation, "// &
                                       "i.e. with the same perturbation.", &
                                       print_level=low_print_level, common_iter_levels=3, each_iter_names=s2a("ITER"), &
                                       add_last=add_last_numeric, each_iter_values=(/3/), filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_linres_section

! **************************************************************************************************
!> \brief creates the input structure used to activate
!>      calculation of induced current  DFPT
!>      Available properties : none
!> \param section the section to create
!> \author  MI/VW
! **************************************************************************************************
   SUBROUTINE create_current_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, print_key, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="current", &
                          description="The induced current density is calculated by DFPT.", &
                          n_keywords=4, n_subsections=1, repeats=.FALSE., &
                          citations=(/Sebastiani2001, Weber2009/))

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the induced current calculation", &
                          usage="&CURRENT T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUGE", &
                          description="The gauge used to compute the induced current within GAPW.", &
                          usage="GAUGE R", &
                          default_i_val=current_gauge_r_and_step_func, &
                          enum_c_vals=s2a("R", "R_AND_STEP_FUNCTION", "ATOM"), &
                          enum_desc=s2a("Position gauge (doesnt work well).", &
                                        "Position and step function for the soft and the local parts, respectively.", &
                                        "Atoms."), &
                          enum_i_vals=(/current_gauge_r, current_gauge_r_and_step_func, current_gauge_atom/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUGE_ATOM_RADIUS", &
                          description="Build the gauge=atom using only the atoms within this radius.", &
                          usage="GAUGE_ATOM_RADIUS 10.0", &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=4.0_dp, unit_str="angstrom"), &
                          unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_OLD_GAUGE_ATOM", &
                          description="Use the old way to compute the gauge.", &
                          usage="USE_OLD_GAUGE_ATOM T", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_CENTER", &
                          description="The orbital center.", &
                          usage="ORBITAL_CENTER WANNIER", &
                          default_i_val=current_orb_center_wannier, &
                          enum_c_vals=s2a("WANNIER", "COMMON", "ATOM", "BOX"), &
                          enum_desc=s2a("Use the Wannier centers.", &
                                        "Use a common center (works only for an isolate molecule).", &
                                        "Use the atoms as center.", &
                                        "Boxing."), &
                          enum_i_vals=(/current_orb_center_wannier, current_orb_center_common, &
                                        current_orb_center_atom, current_orb_center_box/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COMMON_CENTER", &
                          description="The common center ", usage="COMMON_CENTER 0.0 1.0 0.0", &
                          n_var=3, default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp/), type_of_var=real_t, &
                          unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NBOX", &
                          description="How many boxes along each directions ", usage="NBOX 6 6 5", &
                          n_var=3, default_i_vals=(/4, 4, 4/), type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHI_PBC", &
                          description="Calculate the succeptibility correction to the shift with PBC", &
                          usage="CHI_PBC T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FORCE_NO_FULL", &
                          description="Avoid the calculation of the state dependent perturbation term, "// &
                          "even if the orbital centers are set at Wannier centers or at Atom centers", &
                          usage="FORCE_NO_FULL T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SELECTED_STATES_ON_ATOM_LIST", &
                          description="Indexes of the atoms for selecting"// &
                          " the states to be used for the response calculations.", &
                          usage="SELECTED_STATES_ON_ATOM_LIST 1 2 10", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SELECTED_STATES_ATOM_RADIUS", &
                          description="Select all the states included in the given radius around each atoms "// &
                          "in SELECTED_STATES_ON_ATOM_LIST.", &
                          usage="SELECTED_STATES_ATOM_RADIUS 2.0", &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=4.0_dp, unit_str="angstrom"), &
                          unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_CURRENT", &
                          description="Restart the induced current density calculation"// &
                          " from a previous run (not working yet).", &
                          usage="RESTART_CURRENT", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="print results of induced current density calculation", &
                          repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CURRENT_CUBES", &
                                       description="Controls the printing of the induced current density (not working yet).", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components (not working yet).", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESPONSE_FUNCTION_CUBES", &
                                       description="Controls the printing of the response functions (not working yet).", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components (not working yet).", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LU_BOUNDS", &
                          variants=(/"CUBES_LU"/), &
                          description="The lower and upper index of the states to be printed as cube (not working yet).", &
                          usage="CUBES_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/0, -2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LIST", &
                          description="Indexes of the states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes) (not working yet).", &
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

      NULLIFY (subsection)
      CALL create_interp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_current_section

! **************************************************************************************************
!> \brief creates the input structure used to activate
!>         calculation of NMR chemical shift using
!>         the induced current obtained from DFPT
!>      Available properties : none
!> \param section the section to create
!> \author  MI/VW
! **************************************************************************************************
   SUBROUTINE create_nmr_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, print_key, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="nmr", &
                          description="The chemical shift is calculated by DFPT.", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE., &
                          citations=(/Weber2009/))

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the nmr calculation", &
                          usage="&NMR T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTERPOLATE_SHIFT", &
                          description="Calculate the soft part of the chemical shift by interpolation ", &
                          usage="INTERPOLATE_SHIFT T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NICS", &
                          description="Calculate the chemical shift in a set of points  "// &
                          " given from an external file", usage="NICS", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NICS_FILE_NAME", &
                          description="Name of the file with the NICS points coordinates", &
                          usage="NICS_FILE_NAME nics_file", &
                          default_lc_val="nics_file")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_NMR", &
                          description="Restart the NMR calculation from a previous run (NOT WORKING YET)", &
                          usage="RESTART_NMR", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SHIFT_GAPW_RADIUS", &
                          description="While computing the local part of the shift (GAPW), "// &
                          "the integration is restricted to nuclei that are within this radius.", &
                          usage="SHIFT_GAPW_RADIUS 20.0", &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=60.0_dp, unit_str="angstrom"), &
                          unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="print results of nmr calculation", &
                          repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESPONSE_FUNCTION_CUBES", &
                                       description="Controls the printing of the response functions ", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
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
                          description="Indexes of the states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
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

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CHI_TENSOR", &
                                       description="Controls the printing of susceptibility", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "SHIELDING_TENSOR", &
                                       description="Controls the printing of the chemical shift", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_LU_BOUNDS", &
                          variants=(/"ATOMS_LU"/), &
                          description="The lower and upper atomic index for which the tensor is printed", &
                          usage="ATOMS_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/0, -2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_LIST", &
                          description="list of atoms for which the shift is printed into a file ", &
                          usage="LIST_ATOMS 1 2", n_var=-1, &
                          type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_interp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_nmr_section

! **************************************************************************************************
!> \brief creates the input structure used to activate
!>      calculation of NMR spin-spin coupling (implementation not operating)
!>      Available properties : none
!> \param section the section to create
!> \author  VW
! **************************************************************************************************
   SUBROUTINE create_spin_spin_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, print_key, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="spinspin", &
                          description="Compute indirect spin-spin coupling constants.", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the nmr calculation", &
                          usage="&SPINSPIN T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_SPINSPIN", &
                          description="Restart the spin-spin calculation from a previous run (NOT WORKING YET)", &
                          usage="RESTART_SPINSPIN", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ISSC_ON_ATOM_LIST", &
                          description="Atoms for which the issc is computed.", &
                          usage="ISSC_ON_ATOM_LIST 1 2 10", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_FC", &
                          description="Compute the Fermi contact contribution", &
                          usage="DO_FC F", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_SD", &
                          description="Compute the spin-dipolar contribution", &
                          usage="DO_SD F", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_PSO", &
                          description="Compute the paramagnetic spin-orbit contribution", &
                          usage="DO_PSO F", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_DSO", &
                          description="Compute the diamagnetic spin-orbit contribution (NOT YET IMPLEMENTED)", &
                          usage="DO_DSO F", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="print results of the indirect spin-spin calculation", &
                          repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "K_MATRIX", &
                                       description="Controls the printing of the indirect spin-spin matrix", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_LIST", &
                          description="list of atoms for which the indirect spin-spin is printed into a file ", &
                          usage="LIST_ATOMS 1 2", n_var=-1, &
                          type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_interp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_spin_spin_section

! **************************************************************************************************
!> \brief creates the input structure used to activate
!>         calculation of EPR using
!>         the induced current obtained from DFPT
!>      Available properties : none
!> \param section the section to create
!> \author  VW
! **************************************************************************************************
   SUBROUTINE create_epr_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection, subsubsection

      NULLIFY (keyword, print_key, subsection, subsubsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EPR", &
                          description="The g tensor is calculated by DFPT ", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE., &
                          citations=(/Weber2009/))

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the epr calculation", &
                          usage="&EPR T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_EPR", &
                          description="Restart the EPR calculation from a previous run (NOT WORKING)", &
                          usage="RESTART_EPR", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="print results of epr calculation", &
                          repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "NABLAVKS_CUBES", &
                                       description="Controls the printing of the components of nabla v_ks ", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
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

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "G_TENSOR", &
                                       description="Controls the printing of the g tensor", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL create_xc_section(subsubsection)
      CALL section_add_subsection(print_key, subsubsection)
      CALL section_release(subsubsection)

      CALL keyword_create(keyword, __LOCATION__, name="GAPW_MAX_ALPHA", &
                          description="Maximum alpha of GTH potentials allowed on the soft grids ", &
                          usage="GAPW_MAX_ALPHA real", default_r_val=5.0_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOO_RHO_HARD", &
                          description="Whether or not to include the atomic parts of the density "// &
                          "in the SOO part of the g tensor", usage="SOO_RHO_HARD", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESPONSE_FUNCTION_CUBES", &
                                       description="Controls the printing of the response functions ", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
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
                          description="Indexes of the states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
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

      NULLIFY (subsection)
      CALL create_interp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_epr_section

! **************************************************************************************************
!> \brief creates the input structure used to activate
!>      calculation of polarizability tensor DFPT
!>      Available properties : none
!> \param section the section to create
!> \author SL
! **************************************************************************************************
   SUBROUTINE create_polarizability_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, print_key, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="POLAR", &
                          description="Compute polarizabilities.", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE., &
                          citations=(/Putrino2002/))

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the polarizability calculation", &
                          usage="&POLAR T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_RAMAN", &
                          description="Compute the electric-dipole--electric-dipole polarizability", &
                          usage="DO_RAMAN F", &
                          citations=(/Luber2014/), &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC_DIPOLE_OPERATOR", &
                          description="Type of dipole operator: Berry phase(T) or Local(F)", &
                          usage="PERIODIC_DIPOLE_OPERATOR T", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="print results of the polarizability calculation", &
                          repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "POLAR_MATRIX", &
                                       description="Controls the printing of the polarizabilities", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="")

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_interp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_polarizability_section

! **************************************************************************************************
!> \brief creates the section for electron transfer coupling
!> \param section ...
!> \author fschiff
! **************************************************************************************************
   SUBROUTINE create_et_coupling_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ET_COUPLING", &
                          description="specifies the two constraints/restraints for extracting ET coupling elements", &
                          n_keywords=1, n_subsections=4, repeats=.FALSE., citations=(/Kondov2007, Futera2017/))

      NULLIFY (subsection)
      CALL create_ddapc_restraint_section(subsection, "DDAPC_RESTRAINT_A")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_ddapc_restraint_section(subsection, "DDAPC_RESTRAINT_B")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_projection(subsection, "PROJECTION")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE_OF_CONSTRAINT", &
                          description="Specifies the type of constraint", &
                          usage="TYPE_OF_CONSTRAINT DDAPC", &
                          enum_c_vals=s2a("NONE", "DDAPC"), &
                          enum_i_vals=(/do_no_et, do_et_ddapc/), &
                          enum_desc=s2a("NONE", "DDAPC Constraint"), &
                          default_i_val=do_no_et)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing basic info about the method", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_et_coupling_section

! **************************************************************************************************
!> \brief defines input sections for specification of Hilbert space partitioning
!>        in projection-operator approach of electronic coupling calulation
!> \param section pointer to the section data structure
!> \param section_name name of the projection section
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE create_projection(section, section_name)

      ! Routine arguments
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: section_name

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, section_block, section_print

! Routine name for dubug purposes

      ! Sanity check
      CPASSERT(.NOT. ASSOCIATED(section))

      ! Initialization
      NULLIFY (keyword)
      NULLIFY (print_key)
      NULLIFY (section_block)
      NULLIFY (section_print)

      ! Input-file section definition
      CALL section_create(section, __LOCATION__, name=TRIM(ADJUSTL(section_name)), &
                          description="Projection-operator approach fo ET coupling calculation", &
                          n_keywords=0, n_subsections=2, repeats=.FALSE.)

      ! Subsection #0: Log printing
      CALL cp_print_key_section_create(print_key, __LOCATION__, 'PROGRAM_RUN_INFO', &
                                       description="Controls printing of data and informations to log file", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Subsection #1: Atomic blocks
      CALL section_create(section_block, __LOCATION__, name='BLOCK', &
                          description="Part of the system (donor, acceptor, bridge,...)", &
                          n_keywords=2, n_subsections=1, repeats=.TRUE.)
      CALL section_add_subsection(section, section_block)

      ! S#1 - Keyword #1: Atom IDs defining a Hilbert space block
      CALL keyword_create(keyword, __LOCATION__, name='ATOMS', &
                          description="Array of atom IDs in the system part", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t, repeats=.FALSE.)
      CALL section_add_keyword(section_block, keyword)
      CALL keyword_release(keyword)

      ! S#1 - Keyword #1: Atom IDs defining a Hilbert space block
      CALL keyword_create(keyword, __LOCATION__, name='NELECTRON', &
                          description="Number of electrons expected in the system part", &
                          usage="NELECTRON {integer}", default_i_val=0)
      CALL section_add_keyword(section_block, keyword)
      CALL keyword_release(keyword)

      ! S#1 - Subsection #1: Printing setting
      CALL section_create(section_print, __LOCATION__, name='PRINT', &
                          description="Possible printing options in ET system part", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL section_add_subsection(section_block, section_print)

      ! S#1 - S#1 - Keyword #1: MO coefficient on specific atom
      CALL keyword_create(keyword, __LOCATION__, name='MO_COEFF_ATOM', &
                          description="Print out MO coeffiecients on given atom", &
                          usage="MO_COEFF_ATOM {integer} {integer} .. {integer}", &
                          type_of_var=integer_t, n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section_print, keyword)
      CALL keyword_release(keyword)

      ! S#1 - S#1 - Keyword #1: MO coefficient of specific state
      CALL keyword_create(keyword, __LOCATION__, name='MO_COEFF_ATOM_STATE', &
                          description="Print out MO coeffiecients of specific state", &
                          usage="MO_COEFF_ATOM_STATE {integer} {integer} .. {integer}", &
                          type_of_var=integer_t, n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section_print, keyword)
      CALL keyword_release(keyword)

      ! S#1 - S#1 - Subsection #1: Saving MOs to CUBE files
      CALL cp_print_key_section_create(print_key, __LOCATION__, 'MO_CUBES', &
                                       description="Controls saving of MO cube files", &
                                       print_level=high_print_level, filename="")

      ! S#1 - S#1 - S#1 - Keyword #1: Stride
      CALL keyword_create(keyword, __LOCATION__, name='STRIDE', &
                          description="The stride (X,Y,Z) used to write the cube file", &
                          usage="STRIDE {integer} {integer} {integer}", n_var=-1, &
                          default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      ! S#1 - S#1 - S#1 - Keyword #2: List of MO IDs
      CALL keyword_create(keyword, __LOCATION__, name='MO_LIST', &
                          description="Indices of molecular orbitals to save", &
                          usage="MO_LIST {integer} {integer} .. {integer}", &
                          type_of_var=integer_t, n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      ! S#1 - S#1 - S#1 - Keyword #2: Number of unoccupied states
      CALL keyword_create(keyword, __LOCATION__, name='NLUMO', &
                          description="Number of unoccupied molecular orbitals to save", &
                          usage="NLUMO {integer}", default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      ! S#1 - S#1 - S#1 - Keyword #3: Number of occupied states
      CALL keyword_create(keyword, __LOCATION__, name='NHOMO', &
                          description="Number of occupied molecular orbitals to save", &
                          usage="NHOMO {integer}", default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section_print, print_key)
      CALL section_release(print_key)

      ! S#1 - S#1 - Clean
      CALL section_release(section_print)

      ! S#1 - Clean
      CALL section_release(section_block)

      ! S#1 - Subsection #1: Printing setting
      CALL section_create(section_print, __LOCATION__, name='PRINT', &
                          description="Possible printing options in ET", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL section_add_subsection(section, section_print)

      ! Print couplings
      CALL cp_print_key_section_create(print_key, __LOCATION__, 'COUPLINGS', &
                                       description="Controls printing couplings onto file", &
                                       print_level=low_print_level, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section_print, print_key)
      CALL section_release(print_key)

      CALL section_release(section_print)

   END SUBROUTINE create_projection

! **************************************************************************************************
!> \brief creates an input section for tddfpt calculation
!> \param section section to create
!> \par History
!>    * 05.2016 forked from create_tddfpt_section [Sergey Chulkov]
!>    * 08.2016 moved from module input_cp2k_dft [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_tddfpt2_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TDDFPT", &
                          description="Parameters needed to set up the Time-Dependent "// &
                          "Density Functional Perturbation Theory. "// &
                          "Current implementation works for hybrid functionals. "// &
                          "Can be used with Gaussian and Plane Waves (GPW) method only.", &
                          n_keywords=12, n_subsections=4, repeats=.FALSE., &
                          citations=(/Iannuzzi2005/))

      NULLIFY (keyword, print_key, subsection)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the TDDFPT procedure", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Integer
      CALL keyword_create(keyword, __LOCATION__, name="NSTATES", &
                          description="Number of excited states to converge.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximal number of iterations to be performed.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_KV", &
                          description="Maximal number of Krylov space vectors. "// &
                          "Davidson iterations will be restarted upon reaching this limit.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=5000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NLUMO", &
                          description="Number of unoccupied orbitals to consider. "// &
                          " Default is to use all unoccupied orbitals (-1).", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NPROC_STATE", &
                          description="Number of MPI processes to be used per excited state. "// &
                          " Default is to use all processors (0).", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! kernel type
      CALL keyword_create(keyword, __LOCATION__, name="KERNEL", &
                          description="Options to compute the kernel", &
                          usage="KERNEL FULL", &
                          enum_c_vals=s2a("FULL", "sTDA"), &
                          enum_i_vals=(/tddfpt_kernel_full, tddfpt_kernel_stda/), &
                          default_i_val=tddfpt_kernel_full)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OE_CORR", &
                          description="Orbital energy correction potential.", &
                          enum_c_vals=s2a("NONE", "LB94", "GLLB", "SAOP", "SHIFT"), &
                          enum_i_vals=(/oe_none, oe_lb, oe_gllb, oe_saop, oe_shift/), &
                          enum_desc=s2a("No orbital correction scheme is used", &
                                        "van Leeuwen and Baerends. PRA, 49:2421, 1994", &
                                        "Gritsenko, van Leeuwen, van Lenthe, Baerends. PRA, 51:1944, 1995", &
                                        "Gritsenko, Schipper, Baerends. Chem. Phys. Lett., 302:199, 1999", &
                                        "Constant shift of virtual and/or open-shell orbitals"), &
                          default_i_val=oe_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! SHIFTS
      CALL keyword_create(keyword, __LOCATION__, name="EV_SHIFT", &
                          variants=s2a("VIRTUAL_SHIFT"), &
                          description="Constant shift of virtual state eigenvalues.", &
                          usage="EV_SHIFT 0.500", &
                          n_var=1, type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      !
      CALL keyword_create(keyword, __LOCATION__, name="EOS_SHIFT", &
                          variants=s2a("OPEN_SHELL_SHIFT"), &
                          description="Constant shift of open shell eigenvalues.", &
                          usage="EOS_SHIFT 0.200", &
                          n_var=1, type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Real
      CALL keyword_create(keyword, __LOCATION__, name="CONVERGENCE", &
                          description="Target accuracy for excited state energies.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_AMPLITUDE", &
                          description="The smallest excitation amplitude to print.", &
                          n_var=1, type_of_var=real_t, &
                          default_r_val=5.0e-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORTHOGONAL_EPS", &
                          description="The largest possible overlap between the ground state and "// &
                          "orthogonalised excited state wave-functions. Davidson iterations "// &
                          "will be restarted when the overlap goes beyond this threshold in "// &
                          "order to prevent numerical instability.", &
                          n_var=1, type_of_var=real_t, &
                          default_r_val=1.0e-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Logical
      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the TDDFPT calculation if a restart file exists", &
                          n_var=1, type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RKS_TRIPLETS", &
                          description="Compute triplet excited states using spin-unpolarised molecular orbitals.", &
                          n_var=1, type_of_var=logical_t, &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Strings
      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          variants=(/"RESTART_FILE_NAME"/), &
                          description="Name of the wave function restart file, may include a path."// &
                          " If no file is specified, the default is to open the file as generated by"// &
                          " the wave function restart print key.", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! DIPOLE subsection
      CALL section_create(subsection, __LOCATION__, name="DIPOLE_MOMENTS", &
                          description="Parameters to compute oscillator strengths in the dipole approximation.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_FORM", &
                          description="Form of dipole transition integrals.", &
                          default_i_val=tddfpt_dipole_velocity, &
                          enum_c_vals=s2a("BERRY", "LENGTH", "VELOCITY"), &
                          enum_desc=s2a("Based on Berry phase formula (valid for fully periodic molecular systems only)", &
                                        "Length form &lang; i | r | j &rang; (valid for non-periodic molecular systems only)", &
                                        "Velocity form &lang; i | d/dr | j &rang;"), &
                          enum_i_vals=(/tddfpt_dipole_berry, tddfpt_dipole_length, tddfpt_dipole_velocity/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE", &
                          description="Reference point to calculate electric "// &
                          "dipole moments using the dipole integrals in the length form.", &
                          enum_c_vals=s2a("COM", "COAC", "USER_DEFINED", "ZERO"), &
                          enum_desc=s2a("Use Center of Mass", &
                                        "Use Center of Atomic Charges", &
                                        "Use User-defined Point", &
                                        "Use Origin of Coordinate System"), &
                          enum_i_vals=(/use_mom_ref_com, &
                                        use_mom_ref_coac, &
                                        use_mom_ref_user, &
                                        use_mom_ref_zero/), &
                          default_i_val=use_mom_ref_com)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_POINT", &
                          description="User-defined reference point.", &
                          usage="REFERENCE_POINT x y z", &
                          repeats=.FALSE., n_var=3, type_of_var=real_t, unit_str='bohr')
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! kernel XC functional
      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! MGRID subsection
      CALL create_mgrid_section(subsection, create_subsections=.FALSE.)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! sTDA subsection
      CALL create_stda_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! PRINT subsection
      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="Printing of information during the TDDFT run.", repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="PROGRAM_BANNER", &
                                       description="Controls the printing of the banner for TDDFPT program", &
                                       print_level=silent_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="GUESS_VECTORS", &
                                       description="Controls the printing of initial guess vectors.", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="ITERATION_INFO", &
                                       description="Controls the printing of basic iteration information "// &
                                       "during the TDDFT run.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="DETAILED_ENERGY", &
                                       description="Controls the printing of detailed energy information "// &
                                       "during the TDDFT run.", &
                                       print_level=medium_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="RESTART", &
                                       description="Controls the dumping of the MO restart file during TDDFPT. "// &
                                       "By default keeps a short history of three restarts.", &
                                       print_level=low_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("TDDFT_SCF"), each_iter_values=(/10/), &
                                       add_last=add_last_numeric, filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum index of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=3)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="NTO_ANALYSIS", &
                                       description="Perform a natural transition orbital analysis.", &
                                       print_level=medium_print_level)
      CALL keyword_create(keyword, __LOCATION__, name="THRESHOLD", &
                          description="Threshold for sum of NTO eigenvalues considered", &
                          usage="Threshold 0.95", &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.975_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CUBE_FILES", &
                          description="Print NTOs on Cube Files", &
                          usage="CUBE_FILES {logical}", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., type_of_var=logical_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). Provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MOS_MOLDEN", &
                                       description="Write the NTO in Molden file format, for visualisation.", &
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
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_tddfpt2_section

! **************************************************************************************************
!> \brief creates the stda input section (simplified Tamm Dancoff Approximation)
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_stda_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="sTDA", &
                          description="parameters needed and setup for sTDA calculations", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRACTION", &
                          variants=(/"HFX_FRACTION"/), &
                          description="The fraction of Hartree-Fock to add to the total energy. "// &
                          "1.0 implies standard Hartree-Fock if used with XC_FUNCTIONAL NONE. "// &
                          "NOTE: In a mixed potential calculation this should be set to 1.0, otherwise "// &
                          "all parts are multiplied with this factor. ", &
                          usage="FRACTION 0.0", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_EWALD", &
                          description="Use Ewald type method for Coulomb interaction", &
                          usage="DO_EWALD", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_TD_FILTER", &
                          description="Threshold for filtering the transition density matrix", &
                          usage="EPS_TD_FILTER", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_stda_section

END MODULE input_cp2k_properties_dft
