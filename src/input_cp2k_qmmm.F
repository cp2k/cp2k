!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief creates the qmmm section of the input
!> \note
!>      moved out of input_cp2k
!> \par History
!>      10.2005 split out of input_cp2k
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_qmmm
  USE bibliography,                    ONLY: Bernstein2009,&
                                             Bernstein2012,&
                                             Golze2013,&
                                             Laino2005,&
                                             Laino2006
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: &
       add_last_numeric, alpha_imomm_default, charge_scale_factor, &
       debug_print_level, do_fm_mom_conserv_buffer, do_fm_mom_conserv_core, &
       do_fm_mom_conserv_equal_a, do_fm_mom_conserv_equal_f, &
       do_fm_mom_conserv_none, do_fm_mom_conserv_qm, &
       do_multipole_section_off, do_multipole_section_on, do_par_atom, &
       do_par_grid, do_qmmm_center_every_step, do_qmmm_center_max_minus_min, &
       do_qmmm_center_never, do_qmmm_center_pbc_aware, &
       do_qmmm_center_setup_only, do_qmmm_coulomb, do_qmmm_gauss, &
       do_qmmm_image_calcmatrix, do_qmmm_image_iter, do_qmmm_link_gho, &
       do_qmmm_link_imomm, do_qmmm_link_pseudo, do_qmmm_none, do_qmmm_swave, &
       do_qmmm_wall_none, do_qmmm_wall_quadratic, do_qmmm_wall_reflective, &
       gaussian, high_print_level, low_print_level, medium_print_level, &
       radius_qmmm_default, silent_print_level, spline3_nopbc_interp, &
       use_perd_none
  USE input_cp2k_mm,                   ONLY: create_GENPOT_section,&
                                             create_Goodwin_section,&
                                             create_LJ_section,&
                                             create_NONBONDED14_section,&
                                             create_Williams_section
  USE input_cp2k_poisson,              ONLY: create_gspace_interp_section
  USE input_cp2k_subsys,               ONLY: create_cell_section
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
                                             lchar_t,&
                                             logical_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE pw_spline_utils,                 ONLY: no_precond,&
                                             precond_spl3_1,&
                                             precond_spl3_2,&
                                             precond_spl3_3,&
                                             precond_spl3_aint,&
                                             precond_spl3_aint2
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_qmmm'

  PUBLIC :: create_qmmm_section

!***
CONTAINS

! *****************************************************************************
!> \brief Creates the QM/MM section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_qmmm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="qmmm",&
            description="Input for QM/MM calculations.",&
            n_keywords=6, n_subsections=3, repeats=.FALSE., required=.TRUE.,&
            citations=(/Laino2005,Laino2006/),&
            error=error)

       NULLIFY(keyword, subsection )
       CALL keyword_create(keyword, name="E_COUPL",&
            variants=s2a("QMMM_COUPLING","ECOUPL"),&
            description="Specifies the type of the QM - MM electrostatic coupling.",&
            usage="E_COUPL GAUSS",&
            enum_c_vals=s2a( "NONE","COULOMB","GAUSS","S-WAVE"),&
            enum_i_vals=(/ do_qmmm_none, do_qmmm_coulomb, do_qmmm_gauss, do_qmmm_swave/),&
            enum_desc=s2a("Mechanical coupling (i.e. classical point charge based)", &
                          "Using analytical 1/r potential (Coulomb) - not available for GPW/GAPW",&
                          "Using fast gaussian expansion of the electrostatic potential (Erf(r/rc)/r)",&
                          "Using fast gaussian expansion of the s-wave electrostatic potential"),&
            default_i_val=do_qmmm_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_POTENTIAL_FILE_NAME",&
            description="Name of the file containing the potential expansion in gaussians. See the "//&
                        "USE_GEEP_LIB keyword.",&
            usage="MM_POTENTIAL_FILE_NAME {filename}",&
            default_lc_val="MM_POTENTIAL",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="use_geep_lib",&
            description=" This keyword enables the use of the internal GEEP library to generate the "//&
            "gaussian expansion of the MM potential. Using this keyword there's no need to provide "//&
            "the MM_POTENTIAL_FILENAME. It expects a number from 2 to 15 (the number of gaussian funtions"//&
            " to be used in the expansion.",&
            usage="use_geep_lib INTEGER",&
            default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="nocompatibility",&
            description="This keyword disables the compatibility of QM/MM "//&
            "potential between CPMD and CP2K implementations. The compatibility"//&
            " is achieved using an MM potential of the form: Erf[x/rc]/x + (1/rc -2/(pi^1/2*rc))*Exp[-(x/rc)^2] ."//&
            "This keyword has effect only selecting GAUSS E_COUPLING type.",&
            usage="nocompatibility LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_mm_rspace",&
            description="Set the threshold for the collocation of the GEEP gaussian functions."//&
            "this keyword affects only the GAUSS E_COUPLING.",&
            usage="eps_mm_rspace real",&
            default_r_val=1.0E-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SPHERICAL_CUTOFF",&
            description="Set the spherical cutoff for the QMMM electrostatic interaction. "//&
            "This acts like a charge multiplicative factor dependent on cutoff. For MM atoms "//&
            "farther than the SPHERICAL_CUTOFF(1) their charge is zero. The switch is performed "//&
            "with a smooth function: 0.5*(1-TANH((r-[SPH_CUT(1)-20*SPH_CUT(2)])/(SPH_CUT(2)))). "//&
            "Two values are required: the first one is the distance cutoff. The second one controls "//&
            "the stiffness of the smoothing.",&
            usage="SPHERICAL_CUTOFF  <REAL>", default_r_vals=(/-1.0_dp,0.0_dp/), n_var=2,&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="parallel_scheme",&
            description="Chooses the parallel_scheme for the long range Potential ",&
            usage="parallel_scheme (ATOM|GRID)",&
            enum_c_vals=s2a( "ATOM","GRID"),&
            enum_desc=s2a("parallelizes on atoms. grids replicated. "//&
            "Replication of the grids can be quite expensive memory wise if running on a system "//&
            "with limited memory per core. The grid option may be preferred in this case.",&
            "parallelizes on grid slices. atoms replicated."),&
            enum_i_vals=(/ do_par_atom, do_par_grid /),&
            default_i_val=do_par_atom, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Centering keywords
       CALL keyword_create(keyword, name="CENTER",&
            description="This keyword sets when the qm system is automatically "//&
            "centered.  Default is EVERY_STEP.",&
            usage="center (EVERY_STEP|SETUP_ONLY|NEVER)",&
            enum_c_vals=s2a("EVERY_STEP","SETUP_ONLY","NEVER"),&
            enum_desc=s2a("re-center every step",&
                          "center at first step only", &
                          "never center"),&
            enum_i_vals=(/ do_qmmm_center_every_step, do_qmmm_center_setup_only, do_qmmm_center_never /),&
            default_i_val=do_qmmm_center_every_step, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CENTER_TYPE",&
            description="How to do the centering",&
            usage="center_type (MAX_MINUS_MIN|PBC_AWARE_MAX_MINUS_MIN)",&
            enum_c_vals=s2a("MAX_MINUS_MIN","PBC_AWARE_MAX_MINUS_MIN"),&
            enum_desc=s2a("center of box defined by maximum coordinate minus minimum coordinate",&
                          "PBC-aware centering.  Useful  for &QMMM&FORCE_MIXING."),&
            enum_i_vals=(/ do_qmmm_center_max_minus_min, do_qmmm_center_pbc_aware /),&
            default_i_val=do_qmmm_center_max_minus_min, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CENTER_GRID",&
            description="This keyword specifies whether the QM system is centered in units of the grid spacing.",&
            usage="grid_center LOGICAL",&
            default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="initial_translation_vector",&
            description="This keyword specify the initial translation vector to be applied to the system.",&
            usage="initial_translation_vector <REAL> <REAL> <REAL>",&
            n_var=3, default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELTA_CHARGE",&
            description="Additional net charge relative to that specified in DFT section.  Used automatically by force mixing",&
            usage="DELTA_CHARGE q",default_i_val=0,&
            n_var=1,type_of_var=integer_t,repeats=.FALSE.,required=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! NB: remember to create these
       CALL create_qmmm_force_mixing_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)    

       CALL create_qmmm_qm_kinds(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_mm_kinds(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_cell_section(subsection,periodic=use_perd_none,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_periodic_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_link_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_forcefield_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_walls_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_image_charge_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_qmmm_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_qmmm_section

! *****************************************************************************
!> \brief Input section to create MM kinds sections
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_mm_kinds(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_mm_kinds', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="MM_KIND",&
         description="Information about the mm kind in the qm/mm scheme",&
         n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
         error=error)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="The MM  kind",usage="O",n_var=1,type_of_var=char_t,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RADIUS",&
         description="Specifies the radius of the atomic kinds",&
         usage="RADIUS real",n_var=1,type_of_var=real_t,unit_str="angstrom",&
         default_r_val=cp_unit_to_cp2k(RADIUS_QMMM_DEFAULT,"angstrom",error=error),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CORR_RADIUS",&
         description="Specifies the correction radius of the atomic kinds"//&
         " The correction radius is connected to the use of the compatibility keyword.",&
         usage="RADIUS real",n_var=1,type_of_var=real_t,unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_qmmm_mm_kinds

! *****************************************************************************
!> \brief Input section to create FORCE_MIXING sections
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author noam
! *****************************************************************************
  SUBROUTINE create_qmmm_force_mixing_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_qmmm_force_mixing_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: link_subsection, print_key, &
                                                qm_kinds_subsection, &
                                                subsection

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="FORCE_MIXING",&
            description="This section enables and defines parameters for force-mixing based QM/MM,"//&
            " which actually does two conventional QM/MM calculations, on a small "//&
            " and a large QM region, and combines the MM forces from one and QM "//&
            " forces from the other to create a complete set of forces.  Energy is "//&
            " not conserved (although the QM/MM energy from the large QM region calculation is reported) "//&
            " so a proper thermostat (i.e. massive, and able to handle dissipation, such as "//&
            " Adaptive Langevin (AD_LANGEVIN)) must be used. For some propagation algorithms "//&
            " (NVT and REFTRAJ MD ensembles) algorithm is adaptive, "//&
            " including molecules hysteretically based on their instantaneous distance from the core region. "//&
            " Information on core/QM/buffer labels can be written in PDB file using "//&
            " MOTION&PRINT&FORCE_MIXING_LABELS.  Will fail if calculation requires a "//&
            " meaningfull stress, or an energy that is consistent with the forces. "//&
            " For GEO_OPT this means "//&
            " only MOTION&GEO_OPT&TYPE CG, MOTION&GEO_OPT&CG&LINE_SEARCH&TYPE 2PNT, and "//&
            " MOTION&GEO_OPT&CG&LINE_SEARCH&2PNT&LINMIN_GRAD_ONLY T",&
         n_keywords=5, n_subsections=3, repeats=.FALSE., required=.FALSE.,&
         citations=(/Bernstein2009,Bernstein2012/),&
         error=error)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_", &
        description="Enables force-mixing", &
        default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MOMENTUM_CONSERVATION_TYPE",&
        description="How to apply force to get momentum conservation", &
        usage="MOMENTUM_CONSERVATION_TYPE <type>",&
        enum_c_vals=s2a( "NONE","EQUAL_F","EQUAL_A"),&
        enum_i_vals=(/ do_fm_mom_conserv_none, do_fm_mom_conserv_equal_f, do_fm_mom_conserv_equal_a/),&
        enum_desc=s2a("No momentum conservation", &
                      "Equal force on each atom",&
                      "Equal acceleration on each atom"),&
        default_i_val=do_fm_mom_conserv_equal_a, required=.FALSE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MOMENTUM_CONSERVATION_REGION",&
        description="Region to apply correction force to for momentum conservation", &
        usage="MOMENTUM_CONSERVATION_REGION <label>",&
        enum_c_vals=s2a( "CORE","QM","BUFFER"),&
        enum_i_vals=(/ do_fm_mom_conserv_core, do_fm_mom_conserv_QM, do_fm_mom_conserv_buffer/),&
        enum_desc=s2a("Apply to QM core region", &
                      "Apply to full QM (dynamics) region",&
                      "Apply to QM+buffer regions"),&
        default_i_val=do_fm_mom_conserv_QM, required=.FALSE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="R_CORE",&
        description="Specify the inner and outer radii of core QM region. "//&
        " All molecules with any atoms within this distance (hysteretically) of any atoms "//&
        " specified as QM in enclosing QM/MM section  will be core QM atoms in the force-mixing calculation.", &
        usage="R_CORE <real> <real>",n_var=2,type_of_var=real_t,&
        default_r_vals=(/ cp_unit_to_cp2k(0.0_dp,"angstrom",error=error),&
                          cp_unit_to_cp2k(0.0_dp,"angstrom",error=error) /),&
        unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="R_QM",&
        description="Specify the inner and outer radii of QM dynamics region. "//&
        " All molecules with atoms within this distance (hysteretically) of any atoms in "//&
        " core will follow QM dynamics in the force-mixing calculation.", &
        usage="R_QM <real> <real>",n_var=2,type_of_var=real_t,&
        default_r_vals=(/ cp_unit_to_cp2k(0.5_dp,"angstrom",error=error),&
                          cp_unit_to_cp2k(1.0_dp,"angstrom",error=error) /),&
        unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="QM_EXTENDED_SEED_IS_ONLY_CORE_LIST",&
        description="Makes the extended QM zone be defined hysterestically "//&
        " by distance from QM core list (i.e. atoms specified explicitly by "//&
        " user) instead of from full QM core region (specified by user + hysteretic "//&
        " selection + unbreakable bonds)",&
        usage="QM_EXTENDED_SEED_IS_ONLY_CORE_LIST <logical>",n_var=1,type_of_var=logical_t,&
        default_l_val=.FALSE., repeats = .FALSE., required = .FALSE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="R_BUF",&
      description="Specify the inner and outer radii of buffer region. "//&
      " All atoms within this distance (hysteretically) of any QM atoms "//&
      " will be buffer atoms in the force-mixing calculation.", &
      usage="R_BUF <real> <real>",n_var=2,type_of_var=real_t,&
      default_r_vals=(/ cp_unit_to_cp2k(0.5_dp,"angstrom",error=error),&
                      cp_unit_to_cp2k(1.0_dp,"angstrom",error=error) /),&
      unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="QM_KIND_ELEMENT_MAPPING",&
         description="Mapping from elements to QM_KINDs for adaptively included atoms.",&
         usage="QM_KIND_ELEMENT_MAPPING {El} {QM_KIND}",&
         n_var=2,type_of_var=char_t,repeats=.TRUE.,required=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_N_QM",&
         description="Maximum number of QM atoms, for detection of runaway adaptive selection.",&
         usage="MAX_N_QM int",default_i_val=300,&
         n_var=1,type_of_var=integer_t,repeats=.FALSE.,required=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ADAPTIVE_EXCLUDE_MOLECULES",&
         description="List of molecule names to exclude from adaptive regions (e.g. big things like proteins)",&
         usage="ADAPTIVE_EXCLUDE_MOLECULES molec1 molec2 ...",&
         n_var=-1,type_of_var=char_t,repeats=.FALSE.,required=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EXTENDED_DELTA_CHARGE",&
         description="Additional net charge in extended region relative to core (core charge is "//&
         " specified in DFT section, as usual for a convetional QM/MM calculation)",&
         usage="EXTENDED_DELTA_CHARGE q",default_i_val=0,&
         n_var=1,type_of_var=integer_t,repeats=.FALSE.,required=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    ! QM_NON_ADAPTIVE subsection
    NULLIFY(subsection)
    CALL section_create(subsection,name="QM_NON_ADAPTIVE",&
         description="List of atoms always in QM region, non-adaptively",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
         error=error)

    NULLIFY(qm_kinds_subsection)
    CALL create_qmmm_qm_kinds(qm_kinds_subsection,error)
    CALL section_add_subsection(subsection,qm_kinds_subsection,error=error)
    CALL section_release(qm_kinds_subsection,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)    

    ! BUFFER_NON_ADAPTIVE subsection
    NULLIFY(subsection)
    CALL section_create(subsection,name="BUFFER_NON_ADAPTIVE",&
         description="List of atoms always in buffer region, non-adaptively, and any needed LINK sections",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
         error=error)

    NULLIFY(qm_kinds_subsection)
    CALL create_qmmm_qm_kinds(qm_kinds_subsection,error)
    CALL section_add_subsection(subsection,qm_kinds_subsection,error=error)
    CALL section_release(qm_kinds_subsection,error=error)
    NULLIFY(link_subsection)
    CALL create_qmmm_link_section(link_subsection,error)
    CALL section_add_subsection(subsection,link_subsection,error=error)
    CALL section_release(link_subsection,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)    


    ![NB] also need a list?
    ![NB] maybe not list+links , but some sort of link template
    ![NB] also, breakable bonds?
    ! BUFFER_LINKS subsection
    NULLIFY(subsection)
    CALL section_create(subsection,name="BUFFER_LINKS",&
         description="Information about possible links for automatic covalent bond breaking for the buffer QM/MM calculation."//&
         "Ignored - need to implement buffer selection by atom and walking of connectivity data.",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
         error=error)

    NULLIFY(link_subsection)
    CALL create_qmmm_link_section(link_subsection,error)
    CALL section_add_subsection(subsection,link_subsection,error=error)
    CALL section_release(link_subsection,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)    

    ! RESTART_INFO subsection
    NULLIFY(subsection)
    CALL section_create(subsection,name="RESTART_INFO",&
            description="This section provides information about old force-mixing indices and labels, "//&
            "for restarts.", &
         n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

      CALL keyword_create(keyword, name="INDICES",&
          description="Indices of atoms in previous step QM regions.", &
          usage="INDICES 1 2 ...", &
          n_var=-1,type_of_var=integer_t,repeats=.TRUE.,required=.TRUE.,error=error)
      CALL section_add_keyword(subsection,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="LABELS",&
          description="Labels of atoms in previous step QM regions.", &
          usage="LABELS 1 1 ...", &
          n_var=-1,type_of_var=integer_t,repeats=.TRUE.,required=.TRUE.,error=error)
      CALL section_add_keyword(subsection,keyword,error=error)
      CALL keyword_release(keyword,error=error)

     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)    

      ! PRINT subsection, with keys for neighbor list
       CALL section_create(subsection,name="print",&
            description="Section of possible print options in FORCE_MIXING.",&
            n_keywords=0, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"SUBCELL",&
            description="Activates the printing of the subcells used for the"//&
            "generation of neighbor lists.", unit_str="angstrom",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Activates the printing of the neighbor lists used"//&
            " for the hysteretic region calculations.", &
            print_level=high_print_level,filename="",unit_str="angstrom",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)    


  END SUBROUTINE create_qmmm_force_mixing_section

! *****************************************************************************
!> \brief Input section to create QM kinds sections
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_qm_kinds(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_qm_kinds', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="QM_KIND",&
         description="Information about the qm kind in the qm/mm scheme",&
         n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
         error=error)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="The qm kind",usage="O",n_var=1,type_of_var=char_t,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MM_INDEX",&
         description="The indexes of the mm atoms that have this kind. This keyword can be"//&
         " repeated several times (useful if you have to specify many indexes).",&
         usage="MM_INDEX 1 2",&
         n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_qmmm_qm_kinds

! *****************************************************************************
!> \brief Input section to set QM/MM periodic boundary conditions
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_walls_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_walls_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="WALLS",&
            description="Enables Walls for the QM box. This can be used to avoid that QM "//&
            " atoms move out of the QM box.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="WALL_SKIN",&
            description="Specify the value of the skin of the Wall in each dimension. "//&
            "The wall's effect is felt when atoms fall within the skin of the Wall.",&
            usage="WALL_SKIN <real> <real> <real>",n_var=3,type_of_var=real_t,&
            default_r_vals=(/ cp_unit_to_cp2k(0.5_dp,"angstrom",error=error),&
                              cp_unit_to_cp2k(0.5_dp,"angstrom",error=error),&
                              cp_unit_to_cp2k(0.5_dp,"angstrom",error=error) /),&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TYPE",&
            description="Specifies the type of wall",&
            usage="TYPE REFLECTIVE",&
            enum_c_vals=s2a( "NONE","REFLECTIVE","QUADRATIC"),&
            enum_i_vals=(/ do_qmmm_wall_none, do_qmmm_wall_reflective, do_qmmm_wall_quadratic/),&
            enum_desc=s2a("No Wall around QM box", &
                          "Reflective Wall around QM box",&
                          "Quadratic Wall around QM box"),&
            default_i_val=do_qmmm_wall_reflective, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Specify the value of the the force constant for the quadratic wall",&
            usage="K <real>",unit_str='internal_cp2k',&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_qmmm_walls_section

! ****************************************************************************
!> \brief Input section for QM/MM image charge calculations
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Dorothea Golze
! ****************************************************************************
  SUBROUTINE create_qmmm_image_charge_section (section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_qmmm_image_charge_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="IMAGE_CHARGE",&
            description="Inclusion of polarization effects within the image charge "//&
                        "approach for systems where QM molecules are physisorbed on e.g. metal "//&
                        "surfaces described by MM. QM box size has to be equal to MM box size.",&
            n_keywords=3, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Golze2013/),&
            error=error)

       CALL keyword_create(keyword, name="MM_ATOM_LIST",&
            description="List of MM atoms carrying an induced Gaussian charge. "//&
                        "If this keyword is not given, all MM atoms will carry an image charge.",&
            usage="MM_ATOM_LIST 1 2 3 or 1..3 ",n_var=-1,type_of_var=integer_t,&
            repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WIDTH",&
            description="Specifies the value of the width of the (induced) Gaussian "//&
                        "charge distribution carried by each MM atom.",&
            usage="WIDTH <real> ",n_var=1,type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=3.0_dp,unit_str="angstrom^-2", error=error), &
            unit_str="angstrom^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXT_POTENTIAL",&
            description="External potential applied to the metal electrode ",&
            usage="EXT_POTENTIAL <real> ",n_var=1,type_of_var=real_t,&
            default_r_val= 0.0_dp,&
            unit_str="volt",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DETERM_COEFF",&
            description="Specifies how the coefficients are determined.",&
            usage="DETERM_COEFF ITERATIVE",&
            enum_c_vals=s2a( "CALC_MATRIX","ITERATIVE"),&
            enum_i_vals=(/ do_qmmm_image_calcmatrix, do_qmmm_image_iter/),&
            enum_desc=s2a("Calculates image matrix and solves linear set of equations", &
                          "Uses an iterative scheme to calculate the coefficients"),&
            default_i_val=do_qmmm_image_calcmatrix, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_IMAGE_MATRIX",&
            description="Restart the image matrix. Useful when "//&
            "calculating coefficients iteratively (the image matrix "//&
            "is used as preconditioner in that case)",&
            usage="RESTART_IMAGE_MATRIX", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="IMAGE_RESTART_FILE_NAME",&
            description="File name where to read the image matrix used "//&
            "as preconditioner in the iterative scheme",&
            usage="IMAGE_RESTART_FILE_NAME <FILENAME>",&
            type_of_var=lchar_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_qmmm_image_charge_section

! *****************************************************************************
!> \brief Input section to set QM/MM periodic boundary conditions
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_periodic_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_periodic_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(keyword, subsection)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PERIODIC",&
            description="Specify parameters for QM/MM periodic boundary conditions calculations",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Laino2006/),&
            error=error)

       CALL keyword_create(keyword, name="GMAX",&
            description="Specifies the maximum value of G in the reciprocal space over which perform the Ewald sum.",&
            usage="GMAX <real>",n_var=1,default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REPLICA",&
            description="Specifies the number of replica to take into consideration for the real part of the "//&
            "calculation. Default is letting the qmmm module decide how many replica you really need.",&
            usage="REPLICA <integer>",n_var=1,default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NGRIDS",&
            description="Specifies the number of grid points used for the Interpolation of the G-space term",&
            usage="NGRIDS <integer> <iteger> <integer> ",n_var=3,default_i_vals=(/50,50,50/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_multipole_qmmm_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_gspace_interp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,"check_spline",&
            description="Controls the checking of the G-space term Spline Interpolation.",&
            print_level=medium_print_level,filename="GSpace-SplInterp",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_periodic_section

! *****************************************************************************
!> \brief Section to set-up parameters for decoupling using the Bloechl scheme
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      Dorothea Golze [04.2014] copied from input_cp2k_poisson.F and
!>      enabled switch-on/off 
!> \author teo
! *****************************************************************************
  SUBROUTINE create_multipole_qmmm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_multipole_qmmm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="MULTIPOLE",&
            description="This section is used to set up the decoupling of QM periodic images with "//&
            "the use of density derived atomic point charges. Switched on by default even if not "//&
            "explicitly given. Can be switched off if e.g. QM and MM box are of the same size.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Defines the usage of the multipole section",&
            usage="ON",&
            enum_c_vals=s2a("ON","OFF"),&
            enum_i_vals=(/do_multipole_section_on, do_multipole_section_off/),&
            enum_desc=s2a("switch on MULTIPOLE section", &
                          "switch off MULTIPOLE section"), &
            default_i_val=do_multipole_section_on,lone_keyword_i_val=do_multipole_section_on,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Real space cutoff for the Ewald sum.",&
            usage="RCUT {real}", n_var=1, type_of_var=real_t,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EWALD_PRECISION",&
            description="Precision achieved in the Ewald sum.",&
            usage="EWALD_PRECISION {real}", n_var=1, type_of_var=real_t,&
            unit_str="hartree",default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANALYTICAL_GTERM",&
            description="Evaluates the Gterm in the Ewald Scheme analytically instead of using Splines.",&
            usage="ANALYTICAL_GTERM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NGRIDS",&
            description="Specifies the number of grid points used for the Interpolation of the G-space term",&
            usage="NGRIDS <integer> <iteger> <integer> ",n_var=3,default_i_vals=(/50,50,50/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_gspace_interp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,"check_spline",&
            description="Controls the checking of the G-space term Spline Interpolation.",&
            print_level=medium_print_level,filename="GSpace-SplInterp",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,"program_run_info",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_multipole_qmmm_section

! *****************************************************************************
!> \brief creates the qm/mm forcefield section to override to the FF specification
!>      given in the FIST input
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_forcefield_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_qmmm_forcefield_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(subsection,keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FORCEFIELD",&
            description="Specify information on the QM/MM forcefield",&
            n_keywords=0, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="MULTIPLE_POTENTIAL",&
            description="Enables the possibility to define NONBONDED and NONBONDED14 as a"//&
            " sum of different kinds of potential. Useful for piecewise defined potentials.",&
            usage="MULTIPLE_POTENTIAL T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_qmmm_ff_nb_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED14_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_forcefield_section

! *****************************************************************************
!> \brief creates the qm/mm forcefield section to override to the FF specification
!>      given in the FIST input - NONBONDED PART
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_ff_nb_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_ff_nb_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="NONBONDED",&
            description="Specify information on the QM/MM non-bonded forcefield",&
            n_keywords=0, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_GENPOT_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_ff_nb_section

! *****************************************************************************
!> \brief creates the qm/mm link section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_link_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_link_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(keyword, subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="LINK",&
            description="Specify information on the QM/MM link treatment",&
            n_keywords=7, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="QM_INDEX",&
            variants=(/ "QM" /),&
            description="Specifies the index of the QM atom involved in the QM/MM link",&
            usage="QM_INDEX integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_KIND",&
            description="Specifies the element of the QM capping atom involved in the QM/MM link",&
            usage="QM_KIND char",n_var=1,type_of_var=char_t,&
            default_c_val="H",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_INDEX",&
            variants=(/ "MM" /),&
            description="Specifies the index of the MM atom involved in the QM/MM link, Default hydrogen.",&
            usage="MM_INDEX integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Overwrite the specification of the radius only for the MM atom involved in the link."//&
            "Default is to use the same radius as for the specified type.",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Overwrite the specification of the correction radius only for the MM atom involved in the link."//&
            "Default is to use the same correction radius as for the specified type.",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

        CALL keyword_create(keyword, name="LINK_TYPE",&
            variants=(/ "LINK ","TYPE ","LTYPE"/),&
            description="Specifies the method to use to treat the defined QM/MM link",&
            usage="LINK_TYPE char",&
            enum_c_vals=s2a( "IMOMM","GHO","PSEUDO"),&
            enum_i_vals=(/ do_qmmm_link_imomm, do_qmmm_link_gho, do_qmmm_link_pseudo/),&
            default_i_val=do_qmmm_link_imomm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA_IMOMM",&
            variants=s2a("ALPHA"),&
            description="Specifies the scaling factor to be used for projecting the forces "//&
            "on the capping hydrogen in the IMOMM QM/MM link scheme to the MM atom of the link. "//&
            "A good guess can be derived from the bond distances of the forcefield: "//&
            "alpha = r_eq(QM-MM) / r_eq(QM-H).",&
            usage="ALPHA_IMOMM real",n_var=1,type_of_var=real_t,&
            default_r_val=ALPHA_IMOMM_DEFAULT,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QMMM_SCALE_FACTOR",&
            variants=(/ "QMMM_CHARGE_SCALE ",&
                        "QMMM_CHARGE_FACTOR",&
                        "QMMM_SCALE_CHARGE " /),&
            description="Specifies the scaling factor for the MM charge involved in the link QM/MM."//&
            " This keyword affects only the QM/MM potential, it doesn't affect the electrostatic in "//&
            " the classical part of the code. "//&
            "Default 1.0 i.e. no charge rescaling of the MM atom of the QM/MM link bond.",&
            usage="SCALE_FACTOR real",n_var=1,type_of_var=real_t,&
            default_r_val=CHARGE_SCALE_FACTOR,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FIST_SCALE_FACTOR",&
            variants=(/ "FIST_CHARGE_SCALE ",&
                        "FIST_CHARGE_FACTOR",&
                        "FIST_SCALE_CHARGE " /),&
            description="Specifies the scaling factor for the MM charge involved in the link QM/MM."//&
            " This keyword modifies the MM charge in FIST. The modified charge will be used then also"//&
            " for the generation of the QM/MM potential. "//&
            "Default 1.0 i.e. no charge rescaling of the MM atom of the QM/MM link bond.",&
            usage="SCALE_FACTOR real",n_var=1,type_of_var=real_t,&
            default_r_val=CHARGE_SCALE_FACTOR,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="MOVE_MM_CHARGE",&
            description="Specify information to move a classical charge before the"//&
            " QM/MM energies and forces evaluation",&
            n_keywords=4, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_1",&
            variants=(/ "MM1" /),&
            description="Specifies the index of the MM atom involved in the QM/MM link to be moved",&
            usage="ATOM_INDEX_1 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_2",&
            variants=(/ "MM2" /),&
            description="Specifies the index of the second atom defining the direction along which "//&
            " the atom will be moved",&
            usage="ATOM_INDEX_2 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the scaling factor that defines the movement along the defined direction",&
            usage="ALPHA real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius used for the QM/MM electrostatic coupling after movement",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,unit_str="angstrom", default_r_val=0.0_dp, &
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius used for the QM/MM electrostatic coupling after movement",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,unit_str="angstrom", default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="ADD_MM_CHARGE",&
            description="Specify information to add a classical charge before the"//&
            " QM/MM energies and forces evaluation",&
            n_keywords=5, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_1",&
            variants=(/ "MM1" /),&
            description="Specifies the index of the first atom defining the direction along which"//&
            " the atom will be added",&
            usage="ATOM_INDEX_1 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_2",&
            variants=(/ "MM2" /),&
            description="Specifies the index of the second atom defining the direction along which "//&
            " the atom will be added",&
            usage="ATOM_INDEX_2 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the scaling factor that defines the movement along the defined direction",&
            usage="ALPHA real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius used for the QM/MM electrostatic coupling for the added source",&
            usage="RADIUS real",n_var=1,unit_str="angstrom",&
            default_r_val=cp_unit_to_cp2k(RADIUS_QMMM_DEFAULT,"angstrom",error=error),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius used for the QM/MM electrostatic coupling for the added source",&
            usage="RADIUS real",n_var=1,unit_str="angstrom",&
            default_r_val=cp_unit_to_cp2k(RADIUS_QMMM_DEFAULT,"angstrom",error=error),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Specifies the charge for the added source of QM/MM potential",&
            usage="CHARGE real",default_r_val=0.0_dp,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_qmmm_link_section

! *****************************************************************************
!> \brief creates the interpolation section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_qmmm_interp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qmmm_interp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="kind of interpolation used between the multigrids",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="kind",&
            description="the interpolator to use",&
            usage="kind spline3",&
            default_i_val=spline3_nopbc_interp,&
            enum_c_vals=s2a("spline3_nopbc"),&
            enum_i_vals=(/spline3_nopbc_interp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="safe_computation",&
            description="if a non unrolled calculation is to be performed in parallel",&
            usage="safe_computation OFF",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=precond_spl3_aint,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=precond_spl3_3,&
            enum_c_vals=s2a("copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_x",&
            description="accuracy on the solution for spline3 the interpolators",&
            usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_r",&
            description="accuracy on the residual for spline3 the interpolators",&
            usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="max_iter",&
            variants=(/'maxiter'/),&
            description="the maximum number of iterations",&
            usage="max_iter 200", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"conv_info",&
            description="if convergence information about the linear solver"//&
            " of the spline methods should be printed", &
            print_level=medium_print_level,each_iter_names=s2a("SPLINE_FIND_COEFFS"),&
            each_iter_values=(/10/),filename="__STD_OUT__",&
            add_last=add_last_numeric,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"spl_coeffs",&
            description="outputs a cube with the coefficents calculated for "//&
            "the spline interpolation", &
            print_level=debug_print_level,&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_qmmm_interp_section

! *****************************************************************************
!> \brief Create the print qmmm section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_qmmm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_qmmm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword, print_key)
       CALL section_create(section,name="print",&
            description="Section of possible print options specific of the QMMM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Controls the printing of the DIPOLE in a QM/MM calculations."//&
            " It requires that the DIPOLE calculations is "//&
            " requested both for the QS  and for the MM  part.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PGF",&
            description="Controls the printing of the gaussian expansion basis set of the"//&
            " electrostatic potential", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the QMMM  potential",&
            print_level=high_print_level,filename="MM_ELPOT_QMMM",&
            common_iter_levels=1,error=error)

       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MM_POTENTIAL",&
            description="Controls the printing of the MM unidimensional potential on file",&
            print_level=high_print_level,filename="MM_ELPOT",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"QMMM_MATRIX",&
            description="Controls the printing of the QMMM 1 electron Hamiltonian Matrix"//&
            " for methods like semiempirical and DFTB",&
            print_level=high_print_level,filename="__STD_OUT__",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PERIODIC_INFO",&
            description="Controls the printing of information regarding the periodic boundary condition.",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"GRID_INFORMATION",&
            description="Controls the printing of information regarding the PW grid structures"//&
            " for PERIODIC qm/mm calculations.",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"derivatives",&
            description="Print all derivatives after QM/MM calculation", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_charges",&
            description="Print all charges generating the QM/MM potential", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_link_info",&
            description="Print all information on QM/MM links", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qs_derivatives",&
            description="Print QM derivatives after QS calculation", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"image_charge_info",&
            description="Prints image charge coefficients and detailed energy info", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
 
       CALL cp_print_key_section_create(print_key,"image_charge_restart",&
            description="Controls the printing of the restart file for "// &
            "the image matrix when using the iterative scheme",&
            print_level=low_print_level,add_last=add_last_numeric,filename="RESTART",&
            common_iter_levels=3, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_qmmm_section

END MODULE input_cp2k_qmmm
