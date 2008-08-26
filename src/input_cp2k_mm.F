!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  USE bibliography,                    ONLY: Tersoff1988,&
                                             Tosi1964a,&
                                             Tosi1964b,&
                                             Yamada2000
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants
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
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mm'

  PUBLIC :: create_mm_section
  PUBLIC :: create_NONBONDED14_section, create_LJ_section, create_Williams_section, &
            create_Goodwin_section, create_Buck4r_section, create_Buckmorse_section,&
            create_GENPOT_section, create_neighbor_lists_section
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

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Controls the printing of dipole information. Computed with"//&
            " the Berry phase formalism. Results in atomic Units.", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)

       CALL keyword_create(keyword=keyword, name="PERIODIC",&
            description="Use Berry phase formula (PERIODIC=T) or simple operator (PERIODIC=F)",&
            usage="PERIODIC {logical}",repeats=.FALSE.,&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="REFERENCE",&
            variants=s2a("REF"),&
            description="Define the reference point for electrical moment calculations.",&
            usage="REFERENCE COM",&
            enum_c_vals=s2a("COM","COAC","USER_DEFINED","ZERO"),&
            enum_desc=s2a("Use Center of Mass",&
                          "Use Center of Atomic Charges",&
                          "Use User Defined Point (Keyword:REF_POINT)",&
                          "Use Origin of Coordinate System"),&
            enum_i_vals=(/ use_mom_ref_com,&
                           use_mom_ref_coac,&
                           use_mom_ref_user,&
                           use_mom_ref_zero/),&
            default_i_val=use_mom_ref_zero, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="REFERENCE_POINT",&
            variants=s2a("REF_POINT"),&
            description="Fixed reference point for dipole moment calculations.",&
            usage="REFERENCE_POINT x y z", repeats=.FALSE., required=.FALSE.,&
            n_var=3, default_r_vals=(/0._dp,0._dp,0._dp/), unit_str='bohr', error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Activates the printing of the neighbor lists.", &
            print_level=high_print_level,filename="",&
            error=error)
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
            print_level=silent_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
    CALL cp_print_key_section_create(print_key, "ff_info", description=&
         "controls the printing of information in the forcefield settings", &
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
            enum_c_vals=s2a("OFF","CHM","G87","G96"),&
            enum_desc=s2a("Provides force field parameters through the input file",&
                          "Provides force field parameters through an external file with CHARMM format",&
                          "Provides force field parameters through an external file with GROMOS 87 format",&
                          "Provides force field parameters through an external file with GROMOS 96 format"),&
            enum_i_vals=(/do_ff_undef,&
                          do_ff_charmm,&
                          do_ff_g87,&
                          do_ff_g96/),&
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

       CALL keyword_create(keyword, name="SCALE_CUTOFF",&
            description="Enables the scaling of the energy at the cutoff radius."//&
            " In this way the energy at the cutoff radius is zero.",&
            usage="SCALE_CUTOFF <LOGICAL>", default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NONBONDED",&
            description="Controls the computation of the nonbonded interactions.",&
            usage="NONBONDED T",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NONBONDED14",&
            description="Controls the computation of the nonbonded-14 interactions.",&
            usage="NONBONDED14 T",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MULTIPLE_POTENTIAL",&
            description="Enables the possibility to define NONBONDED and NONBONDED14 as a"//&
            " sum of different kinds of potential. Useful for piecewise defined potentials.",&
            usage="MULTIPLE_POTENTIAL T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
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

       CALL create_SHELL_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BOND_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BEND_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_TORSION_section(subsection, error=error)
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
            usage="KIND {HARMONIC,CHARMM,G87,G96}",&
            enum_c_vals=s2a("CHARMM","G87","G96"),&
            enum_desc=s2a("Functional Form (CHARMM|G87|G96): K * [ 1 +  cos[M*PHI] * cos[PHI0]]",&
                          "Functional Form (CHARMM|G87|G96): K * [ 1 +  cos[M*PHI] * cos[PHI0]]",&
                          "Functional Form (CHARMM|G87|G96): K * [ 1 +  cos[M*PHI] * cos[PHI0]]"),&
            enum_i_vals=(/ do_ff_charmm,&
                           do_ff_g87,&
                           do_ff_g96/),&
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
            usage="KIND {HARMONIC,CHARMM,G87,G96}",&
            enum_c_vals=s2a("HARMONIC","CHARMM","G87","G96","CUBIC","MIXED_BEND_STRETCH"),&
            enum_desc=s2a("Functional Form (HARMONIC/G87): 1/2*K*(THETA-THETA0)^2",&
                          "Functional Form (CHARMM): K*(THETA-THETA0)^2",&
                          "Functional Form (HARMONIC/G87): 1/2*K*(THETA-THETA0)^2",&
                          "Functional Form (G96): 1/2*K*(COS(THETA)-THETA0)^2",&
                          "Functional Form (Cubic): K*(THETA-THETA0)**2*(1+CB*(THETA-THETA0))",&
                          "Functional Form (Mixed_Bend_Stretch): K*(THETA-THETA0)**2*(1+CB*(THETA-THETA0))+"//&
                             " KSS*(R12-R012)*(R32-R032)+KBS12*(R12-R012)*(THETA-THETA0)+KBS32*(R32-R032)*(THETA-THETA0)"),&
            enum_i_vals=(/ do_ff_harmonic,&
                           do_ff_charmm,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_cubic,&
                           do_ff_mixed_bend_stretch/),&
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
            usage="CB {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="bohr",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="R032",&
            description="Mixed bend stretch parameter",&
            usage="CB {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="bohr",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KBS12",&
            description="Mixed bend stretch parameter",&
            usage="CB {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="hartree*bohr^-1*rad^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KBS32",&
            description="Mixed bend stretch parameter",&
            usage="CB {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="hartree*bohr^-1*rad^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KSS",&
            description="Mixed bend stretch parameter",&
            usage="CB {real}", default_r_val=0.0_dp , type_of_var=real_t,&
            n_var=1, unit_str="hartree*bohr^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="THETA0",&
            description="Defines the equilibrium angle.",&
            usage="THETA0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='rad',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="UB",&
            description="Specifies the bond potential the MM system",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of Urey-Bradley potential.",&
            usage="KIND {HARMONIC,CHARMM,G87,G96,QUARTIC}",&
            enum_c_vals=s2a("HARMONIC","CHARMM","G87","G96","QUARTIC"),&
            enum_desc=s2a("Functional Form (QUARTIC): (1/2*K1+[1/3*K2+1/4*K3*|R-R0|]*|R-R0|)(R-R0)^2",&
                          "Functional Form (HARMONIC|CHARMM|G87): K*(R-R0)^2",&
                          "Functional Form (HARMONIC|CHARMM|G87): K*(R-R0)^2",&
                          "Functional Form (G96): 1/4*K*(R-R0)^2",&
                          "Functional Form (HARMONIC|CHARMM|G87): K*(R-R0)^2"),&
            enum_i_vals=(/ do_ff_harmonic,&
                           do_ff_charmm,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_quartic/),&
            default_i_val=do_ff_charmm, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential. For QUARTIC potentials"//&
            " three numbers are expected.",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=-1, unit_str="hartree*bohr^-2",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            description="Defines the equilibrium distance.",&
            usage="R0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="bohr", error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_BEND_section

! *****************************************************************************
!> \brief This section specifies the bond of the MM atoms
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author teo
! *****************************************************************************
  SUBROUTINE create_BOND_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BOND_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BOND",&
            description="Specifies the bond potential the MM system",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="KIND",&
            description="Define the kind of bond potential.",&
            usage="KIND HARMONIC",&
            enum_c_vals=s2a("HARMONIC","CHARMM","G87","G96","QUARTIC","MORSE","CUBIC"),&
            enum_desc=s2a("Functional Form (HARMONIC|G87): 1/2*K*(R-R0)^2",&
                          "Functional Form (CHARMM): K*(R-R0)^2",&
                          "Functional Form (HARMONIC|G87): 1/2*K*(R-R0)^2",&
                          "Functional Form (G96): 1/4*K*(R-R0)^2",&
                          "Functional Form (QUARTIC): (1/2*K1+[1/3*K2+1/4*K3*|R-R0|]*|R-R0|)(R-R0)^2",&
                          "Functional Form (Morse): K1*[(1-exp(-K2*(R-R0)))**2-1])",&
                          "Functional Form (Cubic): K*(R-R0)**2*(1+cs*(R-R0)+7/12*(cs**2*(R-R0)**2))"),&
            enum_i_vals=(/ do_ff_harmonic,&
                           do_ff_charmm,&
                           do_ff_g87,&
                           do_ff_g96,&
                           do_ff_quartic, &
                           do_ff_morse, &
                           do_ff_cubic/),&
            default_i_val=do_ff_charmm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential. For QUARTIC potentials"//&
            " three numbers are expected.",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=-1, unit_str="hartree*bohr^-2",error=error)
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
       CALL section_create(section,name="shell",&
            description="This section specifies the parameters for shell-model potentials",&
            n_keywords=6, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

       CALL keyword_create(keyword, name="SPRING",&
            description="Force constant of the spring binding core and shell",&
            usage="SPRING {real}",&
            default_r_val=-1._dp,unit_str="hartree*bohr^-2",error=error)
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

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_IPBV_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BMHFT_section(subsection, error=error)
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
            description="Defines the values of  parameter of the functional form",&
            usage="VALUES ", required=.TRUE., type_of_var=real_t,&
            n_var=-1, repeats=.TRUE., unit_str="internal_cp2k", error=error)
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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARM_FILE_NAME",&
            variants=(/"PARMFILE"/),&
            description="Specifies the filename that contains the parameters of the EAM FF.",&
            usage="PARM_FILE_NAME {FILENAME}",default_lc_val=" ",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_EAM_section

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
            n_var=1, unit_str="angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of the Fumi-Tosi Potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^8",error=error)
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
            description="This section specifies the input parameters for Buckingham 4"//&
            " ranges potential type."//&
            "Functional form: (r<r1) V(r) = A* EXP(-B*r); (r1<r<r2) V(r) = SUM_n POLY1(n)*R^n;"//&
            "(r2<r<r3) V(r) = SUM_n POLY2(n)*R^n; (r>r3) V(r) = -C/r^6",&
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

END MODULE input_cp2k_mm
