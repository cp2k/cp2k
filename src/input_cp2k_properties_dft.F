!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      01.2013 moved out of input_cp2k_dft [MI]
!> \author MI
! *****************************************************************************
MODULE input_cp2k_properties_dft

  USE bibliography,                    ONLY: Putrino2000,&
                                             Putrino2002,&
                                             Sebastiani2001,&
                                             Weber2009
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       current_gauge_atom, current_gauge_r, current_gauge_r_and_step_func, &
       current_orb_center_atom, current_orb_center_box, &
       current_orb_center_common, current_orb_center_wannier, &
       do_ddapc_constraint, do_ddapc_restraint, do_et_becke, do_et_ddapc, &
       do_full_density, do_no_et, do_spin_density, ot_precond_full_all, &
       ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse
  USE input_cp2k_atprop,               ONLY: create_atprop_section
  USE input_cp2k_dft,                  ONLY: create_interp_section,&
                                             create_localize_section
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
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_properties_dft'

  PUBLIC :: create_properties_section

CONTAINS

! *****************************************************************************
!> \brief Create the PROPERTIES section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_properties_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_properties_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PROPERTIES",&
            description="This section is used to set up the PROPERTIES calculation.",&
            n_keywords=4, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY (subsection,keyword)

       CALL create_linres_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_et_coupling_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_resp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_atprop_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,name="FIT_CHARGE",&
            description="This section is used to print the density derived atomic point charges."//&
            "The fit of the charges is controlled through the DENSITY_FITTING section",&
            print_level=high_print_level,filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword, name="TYPE_OF_DENSITY",&
            description="Specifies the type of density used for the fitting",&
            usage="TYPE_OF_DENSITY (FULL|SPIN)",&
            enum_c_vals=s2a("FULL","SPIN"),&
            enum_i_vals=(/ do_full_density, do_spin_density/),&
            enum_desc=s2a("Full density","Spin density"),&
            default_i_val=do_full_density, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_properties_section

! *****************************************************************************
!> \brief creates the input structure used to activate
!>      a linear response calculation
!>      Available properties : none
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MI
! *****************************************************************************
  SUBROUTINE create_linres_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_linres_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,subsection,print_key)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="linres",&
            description="The linear response is used to calculate one of the "//&
            " following properties: nmr, epr, raman, ... ",&
            n_keywords=5, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            citations=(/Putrino2000/),&
            error=error)

       CALL keyword_create(keyword, name="EPS",&
            description="target accuracy for the convergence of the conjugate gradient.",&
            usage="EPS 1.e-6", default_r_val=1.e-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Maximum number of conjugate gradient iteration to be performed for one optimization.",&
            usage="MAX_ITER 200", default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_EVERY",&
            description="Restart the conjugate gradient after the specified number of iterations.",&
            usage="RESTART_EVERY 200", default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            description="Type of preconditioner to be used with all minimization schemes. "//&
            "They differ in effectiveness, cost of construction, cost of application. "//&
            "Properly preconditioned minimization can be orders of magnitude faster than doing nothing.",&
            usage="PRECONDITIONER FULL_ALL",&
            default_i_val=ot_precond_none,&
            enum_c_vals=s2a("FULL_ALL","FULL_SINGLE_INVERSE","FULL_SINGLE","FULL_KINETIC","FULL_S_INVERSE",&
            "NONE"),&
            enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "//&
            "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "//&
            "This preconditioner is recommended for almost all systems, except very large systems where "//&
            "make_preconditioner would dominate the total computational cost.",&
            "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "//&
            "but cheaper to construct, "//&
            "might be somewhat less robust. Recommended for large systems.",&
            "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ",&
            "Cholesky inversion of S and T, fast construction, robust, and relatively good, "//&
            "use for very large systems.",&
            "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.",&
            "skip preconditioning"),&
            enum_i_vals=(/ot_precond_full_all,ot_precond_full_single_inverse,ot_precond_full_single, &
            ot_precond_full_kinetic,ot_precond_s_inverse,ot_precond_none/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_GAP",&
            description="Energy gap estimate [a.u.] for preconditioning",&
            usage="ENERGY_GAP 0.1",&
            default_r_val=0.2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART",&
            description="Restart the response calculation if the restart file exists",&
            usage="RESTART",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WFN_RESTART_FILE_NAME",&
            variants=(/"RESTART_FILE_NAME"/),&
            description="Root of the file names where to read the response functions from"//&
            "which to restart the calculation of the linear response",&
            usage="WFN_RESTART_FILE_NAME <FILENAME>",&
            type_of_var=lchar_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_localize_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_current_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_nmr_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_spin_spin_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_epr_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_polarizability_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,"PRINT","printing of information during the linear response calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic iteration information during the LINRES calculation", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of restart file of the response wavefunction."//&
            "For each set of response functions, i.e. for each perturbation,"//&
            "one different restart file is dumped. These restart files should be"//&
            "employed only to restart the same type of LINRES calculation, "//&
            "i.e. with the same perturbation.", &
            print_level=low_print_level,common_iter_levels=3,each_iter_names=s2a("ITER"),&
            add_last=add_last_numeric,each_iter_values=(/3/),filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_linres_section

! *****************************************************************************
!> \brief creates the input structure used to activate
!>      calculation of induced current  DFPT
!>      Available properties : none
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author  MI/VW
! *****************************************************************************
  SUBROUTINE create_current_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_current_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,print_key, subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="current",&
            description="The induced current density is calculated by DFPT.",&
            n_keywords=4, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Sebastiani2001,Weber2009/),&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the induced current calculation",&
            usage="&CURRENT T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="GAUGE",&
            description="The gauge used to compute the induced current within GAPW.",&
            usage="GAUGE R",&
            default_i_val=current_gauge_r_and_step_func,&
            enum_c_vals=s2a("R", "R_AND_STEP_FUNCTION", "ATOM"),&
            enum_desc=s2a("Position gauge (doesnt work well).",&
                          "Position and step function for the soft and the local parts, respectively.",&
                          "Atoms."),&
            enum_i_vals=(/current_gauge_r,current_gauge_r_and_step_func,current_gauge_atom/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GAUGE_ATOM_RADIUS",&
            description="Build the gauge=atom using only the atoms within this radius.", &
            usage="GAUGE_ATOM_RADIUS 10.0",&
            type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=4.0_dp,unit_str="angstrom",error=error),&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="USE_OLD_GAUGE_ATOM",&
            description="Use the old way to compute the gauge.", &
            usage="USE_OLD_GAUGE_ATOM T",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ORBITAL_CENTER",&
            description="The orbital center.",&
            usage="ORBITAL_CENTER WANNIER",&
            default_i_val=current_orb_center_wannier,&
            enum_c_vals=s2a("WANNIER", "COMMON", "ATOM", "BOX"),&
            enum_desc=s2a("Use the Wannier centers.",&
                          "Use a common center (works only for an isolate molecule).",&
                          "Use the atoms as center.",&
                          "Boxing."),&
            enum_i_vals=(/current_orb_center_wannier,current_orb_center_common,&
                          current_orb_center_atom,current_orb_center_box/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COMMON_CENTER",&
            description="The common center ", usage="COMMON_CENTER 0.0 1.0 0.0",&
            n_var=3,default_r_vals=(/0.0_dp,0.0_dp,0.0_dp/),type_of_var=real_t,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NBOX",&
            description="How many boxes along each directions ", usage="NBOX 6 6 5",&
            n_var=3,default_i_vals=(/4,4,4/),type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHI_PBC",&
            description="Calculate the succeptibility correction to the shift with PBC",&
            usage="CHI_PBC T",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="SELECTED_STATES_ON_ATOM_LIST",&
            description="Indexes of the atoms for selecting"//&
            " the states to be used for the response calculations.",&
            usage="SELECTED_STATES_ON_ATOM_LIST 1 2 10",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="SELECTED_STATES_ATOM_RADIUS",&
            description="Select all the states included in the given radius arround each atoms "//&
            "in SELECTED_STATES_ON_ATOM_LIST.",&
            usage="SELECTED_STATES_ATOM_RADIUS 2.0",&
            type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=4.0_dp,unit_str="angstrom",error=error),&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_CURRENT",&
            description="Restart the induced current density calculation"//&
            " from a previous run (not working yet).",&
            usage="RESTART_CURRENT",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of induced current density calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"CURRENT_CUBES",&
            description="Controls the printing of the induced current density (not working yet).", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components (not working yet).",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESPONSE_FUNCTION_CUBES",&
            description="Controls the printing of the response functions (not working yet).", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components (not working yet).",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
            variants=(/"CUBES_LU"/),&
            description="The lower and upper index of the states to be printed as cube (not working yet).",&
            usage="CUBES_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LIST",&
            description="Indexes of the states to be printed as cube files"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes) (not working yet).",&
            usage="CUBES_LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_current_section

! *****************************************************************************
!> \brief creates the input structure used to activate
!>         calculation of NMR chemical shift using 
!>         the induced current obtained from DFPT
!>      Available properties : none
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author  MI/VW
! *****************************************************************************
  SUBROUTINE create_nmr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_nmr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,print_key, subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nmr",&
            description="The chemical shift is calculated by DFPT.",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Weber2009/),&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the nmr calculation",&
            usage="&NMR T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTERPOLATE_SHIFT",&
            description="Calculate the soft part of the chemical shift by interpolation ",&
            usage="INTERPOLATE_SHIFT T",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NICS",&
            description="Calculate the chemical shift in a set of points  "//&
            " given from an external file", usage="NICS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NICS_FILE_NAME",&
            description="Name of the file with the NICS points coordinates",&
            usage="NICS_FILE_NAME nics_file",&
            default_lc_val="nics_file",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_NMR",&
            description="Restart the NMR calculation from a previous run (NOT WORKING YET)",&
            usage="RESTART_NMR",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SHIFT_GAPW_RADIUS",&
            description="While computing the local part of the shift (GAPW), "//&
            "the integration is restricted to nuclei that are within this radius.", &
            usage="SHIFT_GAPW_RADIUS 20.0",&
            type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=60.0_dp,unit_str="angstrom",error=error),&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of nmr calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"RESPONSE_FUNCTION_CUBES",&
            description="Controls the printing of the response functions ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
            variants=(/"CUBES_LU"/),&
            description="The lower and upper index of the states to be printed as cube",&
            usage="CUBES_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LIST",&
            description="Indexes of the states to be printed as cube files"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="CUBES_LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CHI_TENSOR",&
            description="Controls the printing of susceptibility",&
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SHIELDING_TENSOR",&
            description="Controls the printing of the chemical shift",&
            print_level=low_print_level,add_last=add_last_numeric,filename="",&
            error=error)

       CALL keyword_create(keyword,name="ATOMS_LU_BOUNDS",&
            variants=(/"ATOMS_LU"/),&
            description="The lower and upper atomic index for which the tensor is printed",&
            usage="ATOMS_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_LIST",&
            description="list of atoms for which the shift is printed into a file ",&
            usage="LIST_ATOMS 1 2",n_var=-1,&
            type_of_var=integer_t,repeats=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_nmr_section

! *****************************************************************************
!> \brief creates the input structure used to activate
!>      calculation of NMR spin-spin coupling (implementation not operating)
!>      Available properties : none
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author  VW
! *****************************************************************************
  SUBROUTINE create_spin_spin_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_spin_spin_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,print_key,subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="spinspin",&
            description="Compute indirect spin-spin coupling constants.",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the nmr calculation",&
            usage="&SPINSPIN T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_SPINSPIN",&
            description="Restart the spin-spin calculation from a previous run (NOT WORKING YET)",&
            usage="RESTART_SPINSPIN",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ISSC_ON_ATOM_LIST",&
            description="Atoms for which the issc is computed.",&
            usage="ISSC_ON_ATOM_LIST 1 2 10",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_FC",&
            description="Compute the Fermi contact contribution",&
            usage="DO_FC F",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_SD",&
            description="Compute the spin-dipolar contribution",&
            usage="DO_SD F",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_PSO",&
            description="Compute the paramagnetic spin-orbit contribution",&
            usage="DO_PSO F",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_DSO",&
            description="Compute the diamagnetic spin-orbit contribution (NOT YET IMPLEMENTED)",&
            usage="DO_DSO F",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of the indirect spin-spin calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"K_MATRIX",&
            description="Controls the printing of the indirect spin-spin matrix",&
            print_level=low_print_level,add_last=add_last_numeric,filename="",&
            error=error)

       CALL keyword_create(keyword, name="ATOMS_LIST",&
            description="list of atoms for which the indirect spin-spin is printed into a file ",&
            usage="LIST_ATOMS 1 2",n_var=-1,&
            type_of_var=integer_t,repeats=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_spin_spin_section

! *****************************************************************************
!> \brief creates the input structure used to activate
!>         calculation of EPR using
!>         the induced current obtained from DFPT
!>      Available properties : none
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author  VW
! *****************************************************************************
  SUBROUTINE create_epr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_epr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    NULLIFY(keyword,print_key, subsection, subsubsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EPR",&
            description="The g tensor is calculated by DFPT ",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Weber2009/),&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
           description="controls the activation of the epr calculation",&
           usage="&EPR T",&
           default_l_val=.FALSE.,&
           lone_keyword_l_val=.TRUE.,&
           error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_EPR",&
            description="Restart the EPR calculation from a previous run (NOT WORKING)",&
            usage="RESTART_EPR",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of epr calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"NABLAVKS_CUBES",&
            description="Controls the printing of the components of nabla v_ks ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"G_TENSOR",&
            description="Controls the printing of the g tensor",&
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL create_xc_section(subsubsection,error)
       CALL section_add_subsection(print_key, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL keyword_create(keyword, name="GAPW_MAX_ALPHA",&
            description="Maximum alpha of GTH potentials allowed on the soft grids ",&
            usage="GAPW_MAX_ALPHA real", default_r_val=5.0_dp,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SOO_RHO_HARD",&
            description="Whether or not to include the atomic parts of the density "//&
            "in the SOO part of the g tensor", usage="SOO_RHO_HARD", &
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESPONSE_FUNCTION_CUBES",&
            description="Controls the printing of the response functions ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
            variants=(/"CUBES_LU"/),&
            description="The lower and upper index of the states to be printed as cube",&
            usage="CUBES_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LIST",&
            description="Indexes of the states to be printed as cube files"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="CUBES_LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

     END IF

  END SUBROUTINE create_epr_section

! *****************************************************************************
!> \brief creates the input structure used to activate
!>      calculation of polarizability tensor DFPT
!>      Available properties : none
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author SL
! *****************************************************************************
  SUBROUTINE create_polarizability_section(section,error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_polarizability_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,print_key,subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="POLAR",&
            description="Compute polarizabilities.",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Putrino2002/),&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the polarizability calculation",&
            usage="&POLAR T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="DO_RAMAN",&
            description="Compute the electric-dipole--electric-dipole polarizability",&
            usage="DO_RAMAN F",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of the polarizability calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"POLAR_MATRIX",&
            description="Controls the printing of the polarizabilities",&
            print_level=low_print_level,add_last=add_last_numeric,filename="",&
            error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_polarizability_section

! *****************************************************************************
!> \brief creates the section for electron transfer coupling
!> \param section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_et_coupling_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_et_coupling_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"ET_COUPLING",&
            description="specifies the two constraints/restraints for extracting ET coupling elements",&
            n_keywords=1, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"DDAPC_RESTRAINT_A",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"DDAPC_RESTRAINT_B",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"BECKE_RESTRAINT_A",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"BECKE_RESTRAINT_B",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="TYPE_OF_CONSTRAINT",&
            description="Specifies the type of constraint",&
            usage="TYPE_OF_CONSTRAINT DDAPC",&
            enum_c_vals=s2a("NONE","DDAPC","BECKE"),&
            enum_i_vals=(/ do_no_et, do_et_ddapc, do_et_becke/),&
            enum_desc=s2a("NONE","ddapc_restraint","Sperical potential"),&
            default_i_val=do_no_et,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_et_coupling_section

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param section_name ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_restraint_A(section,section_name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: section_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_restraint_A', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(keyword, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,TRIM(ADJUSTL(section_name)),&
            description="Use DDAPC charges in a restraint (check code for details),"//&
            " section can be repeated, but only one constraint is possible at the moment.",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="force constant of the restraint",&
            usage="STRENGTH {real} ",default_r_val=0.1_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TYPE_OF_DENSITY",&
            description="Specifies the type of density used for the fitting",&
            usage="TYPE_OF_DENSITY (FULL|SPIN)",&
            enum_c_vals=s2a("FULL","SPIN"),&
            enum_i_vals=(/ do_full_density, do_spin_density/),&
            enum_desc=s2a("Full density","Spin density"),&
            default_i_val=do_full_density, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="TARGET",&
            description="target value of the restraint",&
            usage="TARGET {real} ",default_r_val=1._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the list of atoms that is summed in the restraint",&
            usage="ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="COEFF",&
            description="Defines the the coefficient of the atom in the atom list (default is one), currently DDAPC only ",&
            usage="COEFF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_FORM",&
            description="Specifies the functional form of the term added",&
            usage="FUNCTIONAL_FORM RESTRAINT",&
            enum_c_vals=s2a("RESTRAINT","CONSTRAINT"),&
            enum_i_vals=(/ do_ddapc_restraint, do_ddapc_constraint/),&
            enum_desc=s2a("Harmonic potential: s*(q-t)**2","Constraint form: s*(q-t)"),&
            default_i_val=do_ddapc_restraint, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_restraint_A

END MODULE input_cp2k_properties_dft
