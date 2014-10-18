!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief creates the mm section of the input
!> \note
!>      moved out of input_cp2k
!> \par History
!>      04.2004 created
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k_mm
  USE bibliography,                    ONLY: &
       Devynck2012, Dick1958, Foiles1986, Mitchell1993, QUIP_ref, &
       Siepmann1995, Tersoff1988, Tosi1964a, Tosi1964b, Yamada2000
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE force_field_kind_types,          ONLY: &
       do_ff_amber, do_ff_charmm, do_ff_cubic, do_ff_fues, do_ff_g87, &
       do_ff_g96, do_ff_harmonic, do_ff_mixed_bend_stretch, do_ff_mm2, &
       do_ff_mm3, do_ff_mm4, do_ff_morse, do_ff_quartic, do_ff_undef
  USE input_constants,                 ONLY: &
       debug_print_level, high_print_level, low_print_level, &
       medium_print_level, silent_print_level, use_mom_ref_coac, &
       use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero
  USE input_cp2k_poisson,              ONLY: create_poisson_section
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
                                             real_t
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mm'

  PUBLIC :: create_mm_section, create_dipoles_section
  PUBLIC :: create_NONBONDED14_section, create_LJ_section, create_Williams_section,&
            create_Goodwin_section,&
            create_GENPOT_section, create_neighbor_lists_section
  PUBLIC :: create_CHARGE_section
!***
CONTAINS

! *****************************************************************************
!> \brief Create the input section for FIST.. Come on.. Let's get woohooo
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_mm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mm",&
            description="This section contains all information to run a MM calculation.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_forcefield_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_neighbor_lists_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_poisson_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_mm_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_mm_section

! *****************************************************************************
!> \brief Create the print mm section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_mm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_mm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in MM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key, keyword)

       CALL cp_print_key_section_create(print_key,"DERIVATIVES",&
            description="Controls the printing of derivatives.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"EWALD_INFO",&
            description="Controls the printing of Ewald energy components during the "//&
            "evaluation of the electrostatics.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_dipoles_section(print_key,"DIPOLE",medium_print_level,error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Activates the printing of the neighbor lists.", &
            print_level=high_print_level,filename="",unit_str="angstrom",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ITER_INFO",&
            description="Activates the printing of iteration info during the self-consistent "//&
            "calculation of a polarizable forcefield.", &
            print_level=medium_print_level,filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SUBCELL",&
            description="Activates the printing of the subcells used for the"//&
            "generation of neighbor lists.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
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
            print_level=low_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key, "FF_PARAMETER_FILE", description=&
            "Controls the printing of Force Field parameter file", &
            print_level=debug_print_level+1,filename="",common_iter_levels=2,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key, "FF_INFO", description=&
            "Controls the printing of information in the forcefield settings", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)

       CALL keyword_create(keyword,"spline_info",&
            description="if the printkey is active prints information regarding the splines"//&
            " used in the nonbonded interactions",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,"spline_data",&
            description="if the printkey is active prints on separated files the splined function"//&
            " together with the reference one. Useful to check the spline behavior.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_mm_section

! *****************************************************************************
!> \brief Create the forcefield section. This section is useful to set up the
!>      proper force_field for FIST calculations
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_forcefield_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_forcefield_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FORCEFIELD",&
            description="Section specifying information regarding how to set up properly"// &
            " a force_field for the classical calculations.",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="PARMTYPE",&
            description="Define the kind of torsion potential",&
            usage="PARMTYPE {OFF,CHM,G87,G96}",&
            enum_c_vals=s2a("OFF","CHM","G87","G96","AMBER"),&
            enum_desc=s2a("Provides force field parameters through the input file",&
                          "Provides force field parameters through an external file with CHARMM format",&
                          "Provides force field parameters through an external file with GROMOS 87 format",&
                          "Provides force field parameters through an external file with GROMOS 96 format",&
                          "Provides force field parameters through an external file with AMBER format (from v.8 on)"),&
            enum_i_vals=(/do_ff_undef,&
                          do_ff_charmm,&
                          do_ff_g87,&
                          do_ff_g96,&
                          do_ff_amber/),&
            default_i_val=do_ff_undef, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARM_FILE_NAME",&
            description="Specifies the filename that contains the parameters of the FF.",&
            usage="PARM_FILE_NAME {FILENAME}",type_of_var=lchar_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VDW_SCALE14",&
            description="Scaling factor for the VDW 1-4 ",&
            usage="VDW_SCALE14 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EI_SCALE14",&
            description="Scaling factor for the electrostatics 1-4 ",&
            usage="EI_SCALE14 1.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SHIFT_CUTOFF",&
            description="Add a constant energy shift to the real-space "//&
            "non-bonding interactions (both Van der Waals and "//&
            "electrostatic) such that the energy at the cutoff radius is "//&
            "zero. This makes the non-bonding interactions continuous at "//&
            "the cutoff.",&
            usage="SHIFT_CUTOFF <LOGICAL>", default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_NONBONDED",&
            description="Controls the computation of all the real-sapce "//&
                        "(short-range) nonbonded interactions. This also "//&
                        "includes the real-space corrections for excluded "//&
                        "or scaled 1-2, 1-3 and 1-4 interactions. When set "//&
                        "to F, the neighborlists are not created and all "//&
                        "interactions that depend on them are not computed.",&
            usage="DO_NONBONDED T",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="IGNORE_MISSING_CRITICAL_PARAMS",&
            description="Do not abort when critical force-field parameters "//&
            "are missing. CP2K will run as if the terms containing the "//&
            "missing parameters are zero.",&
            usage="IGNORE_MISSING_BOND_PARAMS T",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MULTIPLE_POTENTIAL",&
            description="Enables the possibility to define NONBONDED and NONBONDED14 as a"//&
            " sum of different kinds of potential. Useful for piecewise defined potentials.",&
            usage="MULTIPLE_POTENTIAL T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       !Universal scattering potential at very short distances
       CALL keyword_create(keyword, name="ZBL_SCATTERING",&
            description="A short range repulsive potential is added, to simulate "//&
            "collisions and scattering.",&
            usage="ZBL_SCATTERING T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !
       ! subsections
       !
       CALL create_SPLINE_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED14_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_CHARGE_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_CHARGES_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_SHELL_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BOND_section(subsection, "BOND", error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BEND_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_TORSION_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_IMPROPER_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_OPBEND_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_DIPOLE_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_QUADRUPOLE_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

   END IF

  END SUBROUTINE create_forcefield_section

! *****************************************************************************
!> \brief This section specifies the parameters for the splines
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_SPLINE_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_SPLINE_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="SPLINE",&
            description="specifies parameters to set up the splines used in the"//&
            " nonboned interactions (both pair body potential and many body potential)",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="R0_NB",&
            description="Specify the minimum value of the distance interval "//&
            " that brackets the value of emax_spline.",&
            usage="R0_NB <REAL>",default_r_val=cp_unit_to_cp2k(value=0.9_dp,&
                                        unit_str="bohr",&
                                        error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT_NB",&
            description="Cutoff radius for nonbonded interactions. This value overrides "//&
            " the value specified in the potential definition and is global for all potentials.",&
            usage="RCUT_NB {real}", default_r_val=cp_unit_to_cp2k(value=-1.0_dp,&
                                        unit_str="angstrom",&
                                        error=error),&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EMAX_SPLINE",&
            description="Specify the maximum value of the potential up to which"//&
            " splines will be constructed",&
            usage="EMAX_SPLINE <REAL>",&
            default_r_val=0.5_dp,unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EMAX_ACCURACY",&
            description="Specify the maximum value of energy used to check the accuracy"//&
            " requested through EPS_SPLINE. Energy values larger than EMAX_ACCURACY"//&
            " generally do not  satisfy the requested accuracy",&
            usage="EMAX_ACCURACY <REAL>", default_r_val=0.02_dp,unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SPLINE",&
            description="Specify the threshold for the choice of the number of"//&
            " points used in the splines (comparing the splined value with the "//&
            " analytically evaluated one)",&
            usage="EPS_SPLINE <REAL>", default_r_val=1.0E-7_dp,unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NPOINTS",&
            description="Override the default search for an accurate spline by specifying a fixed number of spline points.",&
            usage="NPOINTS 1024", default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="UNIQUE_SPLINE",&
            description="For few potentials (Lennard-Jones) one global optimal spline is generated instead"//&
            " of different optimal splines for each kind of potential",&
            usage="UNIQUE_SPLINE <LOGICAL>",lone_keyword_l_val=.TRUE., default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_SPLINE_section

! *****************************************************************************
!> \brief This section specifies the torsion of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_TORSION_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TORSION_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TORSION",&
            description="Specifies the torsion potential of the MM system.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the tors.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3} {KIND4}", required=.TRUE., type_of_var=char_t,&
            n_var=4, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of torsion potential",&
            usage="KIND CHARMM",&
            enum_c_vals=s2a("CHARMM","G87","G96","AMBER"),&
            enum_desc=s2a("Functional Form (CHARMM|G87|G96|AMBER): K * [ 1 +  cos[M*PHI - PHI0]]",&
                          "Functional Form (CHARMM|G87|G96|AMBER): K * [ 1 +  cos[M*PHI - PHI0]]",&
                          "Functional Form (CHARMM|G87|G96|AMBER): K * [ 1 +  cos[M*PHI - PHI0]]",&
                          "Functional Form (CHARMM|G87|G96|AMBER): K * [ 1 +  cos[M*PHI - PHI0]]"),&
            enum_i_vals=(/ do_ff_charmm,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_amber/),&
            default_i_val=do_ff_charmm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PHI0",&
            description="Defines the phase of the potential.",&
            usage="PHI0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="rad", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="M",&
            description="Defines the multiplicity of the potential.",&
            usage="M {integer}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_TORSION_section

! *****************************************************************************
!> \brief This section specifies the improper torsion of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author louis vanduyfhuys
! *****************************************************************************
  SUBROUTINE create_IMPROPER_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_IMPROPER_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="IMPROPER",&
            description="Specifies the improper torsion potential of the MM system.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the improper tors.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3} {KIND4}", required=.TRUE., type_of_var=char_t,&
            n_var=4, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of improper torsion potential",&
            usage="KIND CHARMM",&
            enum_c_vals=s2a("CHARMM","G87","G96","HARMONIC"),&
            enum_desc=s2a("Functional Form (CHARMM): K * [ PHI - PHI0 ]**2",&
                          "Functional Form (G87|G96|HARMONIC): 0.5 * K * [ PHI - PHI0 ]**2",&
                          "Functional Form (G87|G96|HARMONIC): 0.5 * K * [ PHI - PHI0 ]**2",&
                          "Functional Form (G87|G96|HARMONIC): 0.5 * K * [ PHI - PHI0 ]**2"),&
            enum_i_vals=(/ do_ff_charmm,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_harmonic/),&
            default_i_val=do_ff_charmm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*rad^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PHI0",&
            description="Defines the phase of the potential.",&
            usage="PHI0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="rad", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_IMPROPER_section

! *****************************************************************************
!> \brief This section specifies the out of plane bend of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author louis vanduyfhuys
! *****************************************************************************
  SUBROUTINE create_OPBEND_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_OPBEND_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="OPBEND",&
            description="Specifies the out of plane bend potential of the MM system."//&
            "(Only defined for atom quadruples which are also defined as an improper"//&
            " pattern in the topology.)",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the opbend.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3} {KIND4}", required=.TRUE., type_of_var=char_t,&
            n_var=4, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of out of plane bend potential",&
            usage="KIND HARMONIC",&
            enum_c_vals=s2a("HARMONIC","MM2","MM3","MM4"),&
            enum_desc=s2a("Functional Form (HARMONIC): 0.5 * K * [ PHI - PHI0 ]**2",&
                          "Functional Form (MM2|MM3|MM4): K * [ PHI - PHI0 ]**2",&
                          "Functional Form (MM2|MM3|MM4): K * [ PHI - PHI0 ]**2",&
                          "Functional Form (MM2|MM3|MM4): K * [ PHI - PHI0 ]**2"),&
            enum_i_vals=(/ do_ff_harmonic,&
                           do_ff_mm2,&
                           do_ff_mm3,&
                           do_ff_mm4/),&
            default_i_val=do_ff_harmonic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*rad^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PHI0",&
            description="Defines the phase of the potential.",&
            usage="PHI0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="rad", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_OPBEND_section

! *****************************************************************************
!> \brief This section specifies the bend of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_BEND_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BEND_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BEND",&
            description="Specifies the bend potential of the MM system.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword,subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bend.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3}", required=.TRUE., type_of_var=char_t,&
            n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of bend potential",&
            usage="KIND HARMONIC",&
            enum_c_vals=s2a("HARMONIC","CHARMM","AMBER","G87","G96","CUBIC","MIXED_BEND_STRETCH","MM3"),&
            enum_desc=s2a("Functional Form (HARMONIC|G87): 1/2*K*(THETA-THETA0)^2",&
                          "Functional Form (CHARMM|AMBER): K*(THETA-THETA0)^2",&
                          "Functional Form (CHARMM|AMBER): K*(THETA-THETA0)^2",&
                          "Functional Form (HARMONIC|G87): 1/2*K*(THETA-THETA0)^2",&
                          "Functional Form (G96): 1/2*K*(COS(THETA)-THETA0)^2",&
                          "Functional Form (CUBIC): K*(THETA-THETA0)**2*(1+CB*(THETA-THETA0))",&
                          "Functional Form (MIXED_BEND_STRETCH): K*(THETA-THETA0)**2*(1+CB*(THETA-THETA0))+"//&
                             " KSS*(R12-R012)*(R32-R032)+KBS12*(R12-R012)*(THETA-THETA0)+KBS32*(R32-R032)*(THETA-THETA0)",&
                          "Functional Form (MM3): 1/2*K*(THETA-THETA0)**2*(1-0.014*(THETA-THETA0)+5.6E-5*(THETA-THETA0)**2"//&
                             " -7.0E-7*(THETA-THETA0)**3+9.0E-10*(THETA-THETA0)**4)+KBS12*(R12-R012)*(THETA-THETA0)+"//&
                             " KBS32*(R32-R032)*(THETA-THETA0)"),&
            enum_i_vals=(/ do_ff_harmonic,&
                           do_ff_charmm,&
                           do_ff_amber,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_cubic,&
                           do_ff_mixed_bend_stretch,&
                           do_ff_mm3/),&
            default_i_val=do_ff_charmm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*rad^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CB",&
            description="Defines the the cubic force constant of the bend",&
            usage="CB {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="rad^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R012",&
            description="Mixed bend stretch parameter",&
            usage="R012 {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="bohr",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="R032",&
            description="Mixed bend stretch parameter",&
            usage="R032 {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="bohr",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KBS12",&
            description="Mixed bend stretch parameter",&
            usage="KBS12 {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="hartree*bohr^-1*rad^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KBS32",&
            description="Mixed bend stretch parameter",&
            usage="KBS32 {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="hartree*bohr^-1*rad^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KSS",&
            description="Mixed bend stretch parameter",&
            usage="KSS {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="hartree*bohr^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="THETA0",&
            description="Defines the equilibrium angle.",&
            usage="THETA0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='rad',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Create the Urey-Bradley section
       CALL create_BOND_section(subsection,"UB",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_BEND_section

! *****************************************************************************
!> \brief This section specifies the bond of the MM atoms
!> \param section the section to create
!> \param label ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_BOND_section(section,label,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: label
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BOND_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: tag
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword)

       IF (TRIM(label)=="UB") THEN
          tag = " Urey-Bradley "
          CALL section_create(section,name=TRIM(label),&
               description="Specifies the Urey-Bradley potential between the external atoms"//&
               " defining the angle",&
               n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
               error=error)

       ELSE
          tag = " Bond "
          CALL section_create(section,name=TRIM(label),&
               description="Specifies the bond potential",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)

          CALL keyword_create(keyword, name="ATOMS",&
               description="Defines the atomic kinds involved in the bond.",&
               usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
               n_var=2, error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END IF

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of"//TRIM(tag)//"potential.",&
            usage="KIND HARMONIC",&
            enum_c_vals=s2a("HARMONIC","CHARMM","AMBER","G87","G96","QUARTIC",&
                            "MORSE","CUBIC","FUES"),&
            enum_desc=s2a("Functional Form (HARMONIC|G87): 1/2*K*(R-R0)^2",&
                          "Functional Form (CHARMM|AMBER): K*(R-R0)^2",&
                          "Functional Form (CHARMM|AMBER): K*(R-R0)^2",&
                          "Functional Form (HARMONIC|G87): 1/2*K*(R-R0)^2",&
                          "Functional Form (G96): 1/4*K*(R^2-R0^2)^2",&
                          "Functional Form (QUARTIC): (1/2*K1+[1/3*K2+1/4*K3*|R-R0|]*|R-R0|)(R-R0)^2",&
                          "Functional Form (MORSE): K1*[(1-exp(-K2*(R-R0)))^2-1])",&
                          "Functional Form (CUBIC): K*(R-R0)^2*(1+cs*(R-R0)+7/12*(cs^2*(R-R0)^2))",&
                          "Functional Form (FUES): 1/2*K*R0^2*(1+R0/R*(R0/R-2))"),&
            enum_i_vals=(/ do_ff_harmonic,&
                           do_ff_charmm,&
                           do_ff_amber,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_quartic, &
                           do_ff_morse, &
                           do_ff_cubic, &
                           do_ff_fues/),&
            default_i_val=do_ff_charmm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential. "//&
            "For MORSE potentials 2 numbers are expected. "//&
            "For QUARTIC potentials 3 numbers are expected.",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=-1, unit_str="internal_cp2k",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CS",&
            description="Defines the cubic stretch term.",&
            usage="CS {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="bohr^-1", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            description="Defines the equilibrium distance.",&
            usage="R0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="bohr", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BOND_section

! *****************************************************************************
!> \brief This section specifies the charge of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_charges_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_charges_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="charges",&
            description="Allow to specify an array of classical charges, thus avoiding the"//&
            " packing and permitting the usage of different charges for same atomic types.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Value of the charge for the individual atom. Order MUST reflect"//&
            " the one specified for the geometry.", repeats=.TRUE., usage="{Real}", &
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_charges_section

! *****************************************************************************
!> \brief This section specifies the charge of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_charge_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_charge_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="charge",&
            description="This section specifies the charge of the MM atoms",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind of the charge.",&
            usage="ATOM {KIND1}", required=.TRUE., type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Defines the charge of the MM atom in electron charge unit.",&
            usage="CHARGE {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_charge_section

! *****************************************************************************
!> \brief This section specifies the isotropic polarizability of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Marcel Baer
! *****************************************************************************
  SUBROUTINE create_quadrupole_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_quadrupole_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="QUADRUPOLE",&
            description="This section specifies that we will perform an SCF quadrupole calculation of the MM atoms. "//&
                        "Needs KEYWORD POL_SCF in POISSON secton",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind of the scf quadrupole.",&
            usage="ATOM {KIND1}", required=.TRUE., type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CPOL",&
            description="Defines the isotropic polarizability of the MM atom.",&
            usage="CPOL {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='internal_cp2k', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_quadrupole_section

! *****************************************************************************
!> \brief This section specifies the isotropic polarizability of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Marcel Baer
! *****************************************************************************
  SUBROUTINE create_dipole_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dipole_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DIPOLE",&
            description="This section specifies that we will perform an SCF dipole calculation of the MM atoms. "//&
                        "Needs KEYWORD POL_SCF in POISSON secton",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind of the scf dipole.",&
            usage="ATOM {KIND1}", required=.TRUE., type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="APOL",&
            description="Defines the isotropic polarizability of the MM atom.",&
            usage="APOL {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='angstrom^3', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_DAMPING_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_dipole_section

! *****************************************************************************
!> \brief This section specifies the idamping parameters for polarizable atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Rodolphe Vuilleumier
! *****************************************************************************
  SUBROUTINE create_damping_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_damping_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DAMPING",&
            description="This section specifies optional electric field damping for the polarizable atoms. ",&
            n_keywords=4, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind for this damping function.",&
            usage="ATOM {KIND1}", required=.TRUE., type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="TYPE",&
            description="Defines the damping type.",&
            usage="TYPE {string}", required=.TRUE., type_of_var=char_t,&
            n_var=1, default_c_val="TANG-TOENNIES", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORDER",&
            description="Defines the order for this damping.",&
            usage="ORDER {integer}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BIJ",&
            description="Defines the BIJ parameter for this damping.",&
            usage="BIJ {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='angstrom^-1', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CIJ",&
            description="Defines the CIJ parameter for this damping.",&
            usage="CIJ {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_damping_section

! *****************************************************************************
!> \brief This section specifies the charge of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_shell_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_shell_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section, name="SHELL",&
            description="This section specifies the parameters for shell-model potentials",&
            n_keywords=6, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            citations=(/Dick1958,Mitchell1993,Devynck2012/),&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The kind for which the shell potential parameters are given ",&
            usage="H", default_c_val="DEFAULT", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORE_CHARGE",&
            variants=(/"CORE"/),&
            description="Partial charge assigned to the core (electron charge units)",&
            usage="CORE_CHARGE {real}",&
            default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SHELL_CHARGE",&
            variants=(/"SHELL"/),&
            description="Partial charge assigned to the shell (electron charge units)",&
            usage="SHELL_CHARGE {real}",&
            default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS_FRACTION",&
            variants=(/"MASS"/),&
            description="Fraction of the mass of the atom to be assigned to the shell",&
            usage="MASS_FRACTION {real}",&
            default_r_val=0.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K2_SPRING",&
            variants=s2a("K2","SPRING"),&
            description="Force constant k2 of the spring potential 1/2*k2*r^2 + 1/24*k4*r^4 "//&
                        "binding a core-shell pair when a core-shell potential is employed.",&
            required=.TRUE.,&
            repeats=.FALSE.,&
            usage="K2_SPRING {real}",&
            default_r_val=-1.0_dp,&
            unit_str="hartree*bohr^-2",&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K4_SPRING",&
            variants=s2a("K4"),&
            description="Force constant k4 of the spring potential 1/2*k2*r^2 + 1/24*k4*r^4 "//&
                        "binding a core-shell pair when a core-shell potential is employed. "//&
                        "By default a harmonic spring potential is used, i.e. k4 is zero.",&
            required=.FALSE.,&
            repeats=.FALSE.,&
            usage="K4_SPRING {real}",&
            default_r_val=0.0_dp,&
            unit_str="hartree*bohr^-4",&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DISTANCE",&
            description="Assign a maximum elongation of the spring, "//&
            "if negative no limit is imposed",&
            usage="MAX_DISTANCE 0.0",&
            default_r_val=-1.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SHELL_CUTOFF",&
            description="Define a screening function to exclude some neighbors "//&
            " of the shell when electrostatic interaction are considered, "//&
            "if negative no screening is operated",&
            usage="SHELL_CUTOFF -1.0",&
            default_r_val=-1.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_shell_section

! *****************************************************************************
!> \brief This section specifies the input parameters for 1-4 NON-BONDED
!>      Interactions
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_NONBONDED14_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_NONBONDED14_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nonbonded14",&
            description="This section specifies the input parameters for 1-4 NON-BONDED interactions.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
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

  END SUBROUTINE create_NONBONDED14_section

! *****************************************************************************
!> \brief This section specifies the input parameters for 1-4 NON-BONDED
!>      Interactions
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_NONBONDED_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_NONBONDED_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nonbonded",&
            description="This section specifies the input parameters for NON-BONDED interactions.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_EAM_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_QUIP_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_IPBV_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BMHFT_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BMHFTD_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Buck4r_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Buckmorse_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_GENPOT_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Tersoff_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Siepmann_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_NONBONDED_section

! *****************************************************************************
!> \brief This section specifies the input parameters for generation of
!>      neighbor lists
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo [07.2007] - Zurich University
! *****************************************************************************
  SUBROUTINE create_neighbor_lists_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_neighbor_lists_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="neighbor_lists",&
            description="This section specifies the input parameters for the construction of"//&
            " neighbor lists.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="VERLET_SKIN",&
            description="Defines the Verlet Skin for the generation of the neighbor lists",&
            usage="VERLET_SKIN {real}", default_r_val=cp_unit_to_cp2k(value=1.0_dp,&
                                        unit_str="angstrom",&
                                        error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="neighbor_lists_from_scratch",&
            description="This keyword enables the building of the neighbouring list from scratch.",&
            usage="neighbor_lists_from_scratch logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GEO_CHECK",&
            description="This keyword enables the check that two atoms are never below the minimum"//&
            " value used to construct the splines during the construction of the neighbouring list."//&
            " Disabling this keyword avoids CP2K to abort in case two atoms are below the minimum "//&
            " value of the radius used to generate the splines.",&
            usage="GEO_CHECK",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_neighbor_lists_section

! *****************************************************************************
!> \brief This section specifies the input parameters for a generic potential form
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_GENPOT_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_GENPOT_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="GENPOT",&
            description="This section specifies the input parameters for a generic potential type."//&
            "A functional form is specified. Mathematical Operators recognized are +, -, *, /, ** "//&
            "or alternatively ^, whereas symbols for brackets must be (). "//&
            "The function parser recognizes the (single argument) Fortran 90 intrinsic functions "//&
            "abs, exp, log10, log, sqrt, sinh, cosh, tanh, sin, cos, tan, asin, acos, atan. "//&
            "Parsing for INTRINSIC functions is CASE INsensitive.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the generic potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTION",&
            description="Specifies the functional form in mathematical notation.",&
            usage="FUNCTION  a*EXP(-b*x^2)/x+D*log10(x)", required=.TRUE., type_of_var=lchar_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VARIABLES",&
            description="Defines the variable of the functional form.",&
            usage="VARIABLES x", required=.TRUE., type_of_var=char_t,&
            n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARAMETERS",&
            description="Defines the parameters of the functional form",&
            usage="PARAMETERS a b D", required=.TRUE., type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VALUES",&
            description="Defines the values of parameter of the functional form",&
            usage="VALUES ", required=.TRUE., type_of_var=real_t,&
            n_var=-1, repeats=.TRUE., unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="UNITS",&
            description="Optionally, allows to define valid CP2K unit strings for each parameter value. "//&
                        "It is assumed that the corresponding parameter value is specified in this unit.",&
            usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K", required=.FALSE., type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the generic potential",&
            usage="RCUT {real}", default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_GENPOT_section

! *****************************************************************************
!> \brief This section specifies the input parameters for EAM  potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_EAM_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_EAM_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EAM",&
            description="This section specifies the input parameters for EAM potential type.",&
            citations=(/Foiles1986/),n_keywords=1, n_subsections=0, repeats=.TRUE., &
            required=.TRUE., error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARM_FILE_NAME",&
            variants=(/"PARMFILE"/),&
            description="Specifies the filename that contains the tabulated EAM potential. "//&
            "File structure: the first line of the potential file contains a title. "//&
            "The second line contains: atomic number, mass and lattice constant. "//&
            "These information are parsed but not used in CP2K. The third line contains: "//&
            "dr: increment of r for the tabulated values of density and phi (assuming r starts in 0) [angstrom]; "//&
            "drho: increment of density for the tabulated values of the embedding function (assuming rho starts "//&
            "in 0) [au_c]; cutoff: cutoff of the EAM potential; npoints: number of points in tabulated. Follow "//&
            "in order npoints lines for rho [au_c] and its derivative [au_c*angstrom^-1]; npoints lines for "//&
            "PHI [ev] and its derivative [ev*angstrom^-1] and npoint lines for the embedded function [ev] "//&
            "and its derivative [ev*au_c^-1].",&
            usage="PARM_FILE_NAME {FILENAME}",default_lc_val=" ",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_EAM_section

! *****************************************************************************
!> \brief This section specifies the input parameters for QUIP  potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_QUIP_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_QUIP_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="QUIP",&
            description="This section specifies the input parameters for QUIP potential type. " //&
                        "Mainly intended for things like GAP corrections to DFT "//&
                        "to achieve correlated-wavefunction-like accuracy. " //&
                        "Requires linking with quip library from <a href=""http://www.libatoms.org"" "//&
                        "target=""_blank"">http://www.libatoms.org</a> .",&
            citations=(/QUIP_ref/),n_keywords=1, n_subsections=0, repeats=.TRUE., &
            required=.TRUE., error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the QUIP potential. "//&
            "For more than 2 elements, &QUIP section must be repeated until each element "//&
            "has been mentioned at least once. Set IGNORE_MISSING_CRITICAL_PARAMS to T "//&
            "in enclosing &FORCEFIELD section to avoid having to list every pair of elements separately.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARM_FILE_NAME",&
            variants=(/"PARMFILE"/),&
            description="Specifies the filename that contains the QUIP potential.",&
            usage="PARM_FILE_NAME {FILENAME}",default_lc_val="quip_params.xml",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INIT_ARGS",&
            description="Specifies the potential initialization arguments for the QUIP potential. "//&
            "If blank (default) first potential defined in QUIP parameter file will be used.",&
            usage="INIT_ARGS",default_c_vals=(/""/),&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CALC_ARGS",&
            description="Specifies the potential calculation arguments for the QUIP potential.",&
            usage="CALC_ARGS",default_c_vals=(/""/), &
            n_var=-1,type_of_var=char_t, required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_QUIP_section


! *****************************************************************************
!> \brief This section specifies the input parameters for Lennard-Jones potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_LJ_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_LJ_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="lennard-jones",&
            description="This section specifies the input parameters for LENNARD-JONES potential type."//&
            "Functional form: V(r) = 4.0 * EPSILON * [(SIGMA/r)^12-(SIGMA/r)^6].",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSILON",&
            description="Defines the EPSILON parameter of the LJ potential",&
            usage="EPSILON {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIGMA",&
            description="Defines the SIGMA parameter of the LJ potential",&
            usage="SIGMA {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the LJ potential",&
            usage="RCUT {real}", default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_LJ_section

! *****************************************************************************
!> \brief This section specifies the input parameters for Williams potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_Williams_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Williams_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="williams",&
            description="This section specifies the input parameters for WILLIAMS potential type."//&
            "Functional form: V(r) = A*EXP(-B*r) - C / r^6 .",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of the Williams potential",&
            usage="A {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of the Williams potential",&
            usage="B {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of the Williams potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Williams potential",&
            usage="RCUT {real}",default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Williams_section

! *****************************************************************************
!> \brief This section specifies the input parameters for Goodwin potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_Goodwin_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Goodwin_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="goodwin",&
            description="This section specifies the input parameters for GOODWIN potential type."//&
            "Functional form: V(r) = EXP(M*(-(r/DC)**MC+(D/DC)**MC))*VR0*(D/r)**M.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VR0",&
            description="Defines the VR0 parameter of the Goodwin potential",&
            usage="VR0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of the Goodwin potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DC",&
            description="Defines the DC parameter of the Goodwin potential",&
            usage="DC {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="M",&
            description="Defines the M parameter of the Goodwin potential",&
            usage="M {real}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MC",&
            description="Defines the MC parameter of the Goodwin potential",&
            usage="MC {real}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Goodwin potential",&
            usage="RCUT {real}",default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Goodwin_section

! *****************************************************************************
!> \brief This section specifies the input parameters for IPBV potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_ipbv_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ipbv_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ipbv",&
            description="This section specifies the input parameters for IPBV potential type."//&
            "Functional form: Implicit table function.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the IPBV nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the IPBV potential",&
            usage="RCUT {real}",default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_ipbv_section

! *****************************************************************************
!> \brief This section specifies the input parameters for BMHFT potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_BMHFT_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BMHFT_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BMHFT",&
            description="This section specifies the input parameters for BMHFT potential type."//&
            "Functional form: V(r) = A * EXP(-B*r) - C/r^6 - D/r^8."//&
            "Values available inside cp2k only for the Na/Cl pair.",&
            citations=(/Tosi1964a,Tosi1964b/),n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the BMHFT nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAP_ATOMS",&
            description="Defines the kinds for which internally is defined the BMHFT nonbond potential"//&
            " at the moment only Na and Cl.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the BMHFT potential",&
            usage="RCUT {real}", default_r_val=7.8_dp,&
            unit_str="angstrom",  error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of the Fumi-Tosi Potential",&
            usage="A {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of the Fumi-Tosi Potential",&
            usage="B {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of the Fumi-Tosi Potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of the Fumi-Tosi Potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*angstrom^8",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BMHFT_section

! *****************************************************************************
!> \brief This section specifies the input parameters for BMHFTD potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mathieu Salanne 05.2010
! *****************************************************************************
  SUBROUTINE create_BMHFTD_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BMHFTD_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BMHFTD",&
            description="This section specifies the input parameters for BMHFTD potential type."//&
            "Functional form: V(r) = A*exp(-B*r) - f_6*(r)C/r^6 - f_8(r)*D/r^8."//&
            "where f_order(r)=1-exp(-BD * r) * \sum_{k=0}^order (BD * r)^k / k! ."//&
            "(Tang-Toennies damping function)"//&
            "No values available inside cp2k.",&
            citations=(/Tosi1964a,Tosi1964b/),n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the BMHFTD nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAP_ATOMS",&
            description="Defines the kinds for which internally is defined the BMHFTD nonbond potential"//&
            " at the moment no species included.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the BMHFTD potential",&
            usage="RCUT {real}", default_r_val=7.8_dp,&
            unit_str="angstrom",  error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of the dispersion-damped Fumi-Tosi Potential",&
            usage="A {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of the dispersion-damped Fumi-Tosi Potential",&
            usage="B {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of the dispersion-damped Fumi-Tosi Potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of the dispersion-damped Fumi-Tosi Potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*angstrom^8",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BD",&
            description="Defines the BD parameter of the dispersion-damped Fumi-Tosi Potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORDER",&
            description="Defines the order for this damping.",&
            usage="ORDER {integer}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BMHFTD_section

! *****************************************************************************
!> \brief This section specifies the input parameters for Buckingham 4 ranges potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MI
! *****************************************************************************
  SUBROUTINE create_Buck4r_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Buck4r_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BUCK4RANGES",&
            description="This section specifies the input parameters for the Buckingham 4-ranges"//&
            " potential type."//&
            " Functional form:<ul>"//&
            "<li>V(r) = A*EXP(-B*r) for r &lt; r<sub>1</sub></li>"//&
            "<li>V(r) = Sum_n POLY1(n)*r<sup>n</sup> for r<sub>1</sub> &le; r &lt; r<sub>2</sub></li>"//&
            "<li>V(r) = Sum_n POLY2(n)*r<sup>n</sup> for r<sub>2</sub> &le; r &lt; r<sub>3</sub></li>"//&
            "<li>V(r) = -C/r<sup>6</sup> for r &ge; r<sub>3</sub></li></ul>",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of the Buckingham potential",&
            usage="A {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of the Buckingham potential",&
            usage="B {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of the Buckingham  potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R1",&
            description="Defines the upper bound of the first range ",&
            usage="R1 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R2",&
            description="Defines the upper bound of the second range ",&
            usage="R2 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R3",&
            description="Defines the upper bound of the third range ",&
            usage="R3 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="POLY1",&
           description="Coefficients of the polynomial used in the second range"//&
           "This keyword can be repeated several times.",&
           usage="POLY1 C1 C2 C3 ..",&
           n_var=-1,unit_str="K_e",type_of_var=real_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="POLY2",&
           description="Coefficients of the polynomial used in the third range"//&
           "This keyword can be repeated several times.",&
           usage="POLY1 C1 C2 C3 ..",&
           n_var=-1,unit_str="K_e",type_of_var=real_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Buckingham potential",&
            usage="RCUT {real}",default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Buck4r_section

! *****************************************************************************
!> \brief This section specifies the input parameters for Buckingham + Morse potential type
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MI
! *****************************************************************************
  SUBROUTINE create_Buckmorse_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Buckmorse_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BUCKMORSE",&
            description="This section specifies the input parameters for"//&
            " Buckingham plus Morse potential type "//&
            " Functional Form: V(r) = F0*(B1+B2)*EXP([A1+A2-r]/[B1+B2])-C/r^6+D*{EXP[-2*beta*(r-R0)]-2*EXP[-beta*(r-R0)]}.",&
            citations=(/Yamada2000/),n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="F0",&
            description="Defines the f0 parameter of Buckingham+Morse potential",&
            usage="F0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e*angstrom^-1", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A1",&
            description="Defines the A1 parameter of Buckingham+Morse potential",&
            usage="A1 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A2",&
            description="Defines the A2 parameter of Buckingham+Morse potential",&
            usage="A2 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B1",&
            description="Defines the B1 parameter of Buckingham+Morse potential",&
            usage="B1 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B2",&
            description="Defines the B2 parameter of Buckingham+Morse potential",&
            usage="B2 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of Buckingham+Morse  potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the amplitude for the Morse part ",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            description="Defines the equilibrium distance for the Morse part ",&
            usage="R0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Beta",&
            description="Defines the width for the Morse part ",&
            usage="Beta {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Buckingham potential",&
            usage="RCUT {real}",default_r_val=cp_unit_to_cp2k(value=10.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Defines the lower bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMIN {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Defines the upper bound of the potential. If not set the range is the"//&
            " full range generate by the spline", required=.TRUE., usage="RMAX {real}", &
            type_of_var=real_t,  unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Buckmorse_section

! *****************************************************************************
!> \brief This section specifies the input parameters for Tersoff potential type
!>      (Tersoff, J. PRB 39(8), 5566, 1989)
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_Tersoff_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Tersoff_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TERSOFF",&
            description="This section specifies the input parameters for Tersoff potential type.",&
            citations=(/Tersoff1988/),n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of Tersoff potential",&
            usage="A {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=1.8308E3_dp,&
                                           unit_str="eV",&
                                           error=error),&
            n_var=1, unit_str="eV", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of Tersoff potential",&
            usage="B {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=4.7118E2_dp,&
                                           unit_str="eV",&
                                           error=error),&
            n_var=1, unit_str="eV", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="lambda1",&
            description="Defines the lambda1 parameter of Tersoff potential",&
            usage="lambda1 {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=2.4799_dp,&
                                          unit_str="angstrom^-1",&
                                          error=error),&
            n_var=1, unit_str="angstrom^-1", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="lambda2",&
            description="Defines the lambda2 parameter of Tersoff potential",&
            usage="lambda2 {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=1.7322_dp,&
                                          unit_str="angstrom^-1",&
                                          error=error),&
            n_var=1, unit_str="angstrom^-1", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="alpha",&
            description="Defines the alpha parameter of Tersoff potential",&
            usage="alpha {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=0.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="beta",&
            description="Defines the beta parameter of Tersoff potential",&
            usage="beta {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=1.0999E-6_dp,&
            n_var=1, unit_str="", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n",&
            description="Defines the n parameter of Tersoff potential",&
            usage="n {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=7.8734E-1_dp,&
            n_var=1, unit_str="", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="c",&
            description="Defines the c parameter of Tersoff potential",&
            usage="c {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=1.0039E5_dp,&
            n_var=1, unit_str="", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="d",&
            description="Defines the d parameter of Tersoff potential",&
            usage="d {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=1.6218E1_dp,&
            n_var=1, unit_str="", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="h",&
            description="Defines the h parameter of Tersoff potential",&
            usage="h {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=-5.9826E-1_dp,&
            n_var=1, unit_str="", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="lambda3",&
            description="Defines the lambda3 parameter of Tersoff potential",&
            usage="lambda3 {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=1.7322_dp,&
                                          unit_str="angstrom^-1",&
                                          error=error),&
            n_var=1, unit_str="angstrom^-1", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="bigR",&
            description="Defines the bigR parameter of Tersoff potential",&
            usage="bigR {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=2.85_dp,&
                                          unit_str="angstrom",&
                                          error=error),&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="bigD",&
            description="Defines the D parameter of Tersoff potential",&
            usage="bigD {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=0.15_dp,&
                                          unit_str="angstrom",&
                                          error=error),&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the tersoff potential. "//&
            " This parameter is in principle already defined by the values of "//&
            " bigD and bigR. But it is necessary to define it when using the tersoff "//&
            " in conjuction with other potentials (for the same atomic pair) in order to have"//&
            " the same consistent definition of RCUT for all potentials.",&
            usage="RCUT {real}",required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Tersoff_section

! *****************************************************************************
!> \brief This section specifies the input parameters for Siepmann-Sprik
!>        potential type
!>      (Siepmann and Sprik, J. Chem. Phys. 102(1) 511, 1995)  
!> \param section ...
!> \param error ...
! *****************************************************************************
 SUBROUTINE create_Siepmann_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Siepmann_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

   failure=.FALSE.

   CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
   IF (.NOT. failure) THEN
      CALL section_create(section,name="SIEPMANN",&
           description="This section specifies the input parameters for the Siepmann-Sprik"//&
                       " potential type. Consist of 4 terms: T1+T2+T3+T4. The terms T1=A/rij^alpha"//&
                       " and T2=-C/rij^6 have to be given via the GENPOT section. The terms T3+T4"//& 
                       " are obtained from the SIEPMANN section.",&
           citations=(/Siepmann1995/),n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
           error=error)

      NULLIFY(keyword)
 
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of Siepmann potential",&
            usage="B {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=0.6_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of Siepmann potential",&
            usage="D {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=3.688388_dp,&
                                           unit_str="internal_cp2k",&
                                           error=error),&
            n_var=1, unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="E",&
            description="Defines the E parameter of Siepmann potential",&
            usage="E {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=9.069025_dp,&
                                           unit_str="internal_cp2k",&
                                           error=error),&
            n_var=1, unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="F",&
            description="Defines the F parameter of Siepmann potential",&
            usage="B {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=13.3_dp,n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
!
       CALL keyword_create(keyword, name="beta",&
            description="Defines the beta parameter of Siepmann potential",&
            usage="beta {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=10.0_dp, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
!
       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of Siepmann potential",&
            usage="RCUT {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=3.2_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
!
   END IF
 END SUBROUTINE create_Siepmann_section

! *****************************************************************************
!> \brief creates the input section for the qs part
!> \param print_key ...
!> \param label ...
!> \param print_level ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_dipoles_section(print_key,label,print_level,error)
    TYPE(section_type), POINTER              :: print_key
    CHARACTER(LEN=*), INTENT(IN)             :: label
    INTEGER, INTENT(IN)                      :: print_level
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dipoles_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_print_key_section_create(print_key_section=print_key,name=TRIM(label),&
            description="Section controlling the calculation of "//TRIM(label)//"."//&
            " Note that the result in the periodic case might be defined modulo a certain period,"//&
            " determined by the lattice vectors. During MD, this can lead to jumps.",&
            print_level=print_level,filename="__STD_OUT__",error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword=keyword,&
            name="PERIODIC",&
            description="Use Berry phase formula (PERIODIC=T) or simple operator (PERIODIC=F). "//&
                        "The latter normally requires that the CELL is periodic NONE.",&
            usage="PERIODIC {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="REFERENCE",&
            variants=s2a("REF"),&
            description="Define the reference point for the calculation of the electrostatic moment.",&
            usage="REFERENCE COM",&
            enum_c_vals=s2a("COM","COAC","USER_DEFINED","ZERO"),&
            enum_desc=s2a("Use Center of Mass",&
                          "Use Center of Atomic Charges",&
                          "Use User Defined Point (Keyword:REF_POINT)",&
                          "Use Origin of Coordinate System"),&
            enum_i_vals=(/ use_mom_ref_com,&
                           use_mom_ref_coac,&
                           use_mom_ref_user,&
                           use_mom_ref_zero /),&
            default_i_val=use_mom_ref_zero, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="REFERENCE_POINT",&
            variants=s2a("REF_POINT"),&
            description="Fixed reference point for the calculations of the electrostatic moment.",&
            usage="REFERENCE_POINT x y z", &
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=3,default_r_vals=(/0._dp,0._dp,0._dp/),&
            type_of_var=real_t, &
            unit_str='bohr',&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_dipoles_section

END MODULE input_cp2k_mm
