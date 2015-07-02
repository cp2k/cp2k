!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_free_energy
  USE bibliography,                    ONLY: BarducBus2008,&
                                             VandenCic2006
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       do_fe_ac, do_fe_meta, do_fe_ui, do_wall_gaussian, do_wall_m, &
       do_wall_none, do_wall_p, do_wall_quadratic, do_wall_quartic, &
       do_wall_reflective, gaussian
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
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_free_energy'

PUBLIC :: create_metavar_section,&
          create_fe_section

!***
CONTAINS

! *****************************************************************************
!> \brief creates the free energy section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_fe_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_fe_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY (subsection,keyword,print_key)
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="free_energy",&
            description="Controls the calculation of free energy and free energy derivatives"//&
            " with different possible methods",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="METHOD",&
            description="Defines the method to use to compute free energy.",&
            usage="METHOD (METADYN|UI|AC)",&
            enum_c_vals=s2a( "METADYN","UI","AC"),&
            enum_i_vals=(/do_fe_meta,do_fe_ui,do_fe_ac/),&
            enum_desc=s2a("Metadynamics",&
                          "Umbrella Integration",&
                          "Alchemical Change"),&
            default_i_val=do_fe_meta,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_metadyn_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_ui_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_ac_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(print_key,"free_energy_info",&
            description="Controls the printing of basic and summary information during the"//&
            " Free Energy calculation", &
            print_level=low_print_level,each_iter_names=s2a("MD"),&
            each_iter_values=(/1/),add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_fe_section

! *****************************************************************************
!> \brief creates the metadynamics section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_metadyn_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_metadyn_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="metadyn",&
            description="This section sets parameters to set up a calculation of metadynamics.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., &
            citations=(/VandenCic2006/),error=error)

       NULLIFY(subsection,keyword,print_key)

       CALL keyword_create(keyword, name="USE_PLUMED",&
            description="Specify whether to use plumed as an external metadynamics driver.",&
            usage="USE_PLUMED .FALSE./.TRUE.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PLUMED_INPUT_FILE",&
            description="Specify the file name of the external plumed input file",&
            usage="PLUMED_INPUT_FILE ./FILENAME",&
            default_c_val="./plumed.dat",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MIN_NT_HILLS",&
            description="Specify the minimum MD step interval between spawning "//&
            "two hills. If specified, it must be >= than NT_HILLS. In case MIN_DISP "//&
            "is used, if MIN_DISP is satisfied before MIN_NT_HILLS MD steps have been "//&
            "performed, the MD will continue without any spawning until MIN_NT_HILLS is "//&
            "reached. The default value has the net effect of skipping this check.",&
            usage="MIN_NT_HILLS {integer}",default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NT_HILLS",&
            description="Specify the maximum MD step interval between spawning "//&
            "two hills. When negative, no new hills are spawned and only "//&
            "the hills read from SPAWNED_HILLS_* are in effect. The latter"//&
            "is useful when one wants to add a custom constant bias potential.",&
            usage="NT_HILLS {integer}",default_i_val=30,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="If a Lagrangian scheme is used the temperature for the collective "//&
            "variables is specified. ",usage="TEMPERATURE <REAL>",&
            default_r_val=0.0_dp,unit_str='K',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !RG Adaptive hills
       CALL keyword_create(keyword, name="MIN_DISP",&
            description="Minimum displacement between hills before placing a new hill.",&
            usage="MIN_DISP <REAL>",&
            default_r_val=-1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="OLD_HILL_NUMBER",&
            description="Index of the last hill spawned for this walker.Needed to calculate MIN_DISP",&
            usage="OLD_HILL_NUMBER <INT>",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="OLD_HILL_STEP",&
            description="Timestep of the last hill spawned for this walker.Needed to calculate MIN_DISP",&
            usage="OLD_HILL_STEP <INT>",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       !RG Adaptive hills

       !Hills tail damping
       CALL keyword_create(keyword, name="HILL_TAIL_CUTOFF",&
            description="By setting this variable larger than 0 the tail of the Gaussian hill "//&
            " is damped to zero faster. The Gaussian function is multiplied by a cutoff function"//&
            " that becomes active at |x-X0|>HILL_TAIL_CUTOFF*SCALE, where X0 is the location of "//&
            "the Gaussian and SCALE is the width of the Gaussian. For more than one METAVAR"//&
            " X0 and SCALE are METAVAR-dependent. "//&
            " (1-(|x-X0|/HILL_TAIL_CUTOFF*SCALE)^P_EXP)/(1-(|x-X0|/HILL_TAIL_CUTOFF*SCALE)^Q_EXP)",&
            usage="HILL_TAIL_CUTOFF <REAL>",&
            default_r_val=-1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="P_EXPONENT",&
            description="Exponent at the numerator of the cutoff function to damp the tail of the Gaussian.",&
            usage="P_EXPONENT <INT>",&
            default_i_val=8,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="Q_EXPONENT",&
            description="Exponent at the denominator of the cutoff function to damp the tail of the Gaussian.",&
            usage="Q_EXPONENT <INT>",&
            default_i_val=20,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SLOW_GROWTH",&
            description="Let the last hill grow slowly over NT_HILLS. ",&
            usage="SLOW_GROWTH {logical}",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="TEMP_TOL",&
            description="If a Lagrangian scheme is used the temperature tolerance for the collective "//&
            "variables is specified.",usage="TEMP_TOL <REAL>",&
            unit_str='K',default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LANGEVIN",&
            description="If a Lagrangian scheme is used the eq. motion of the COLVARS are integrated "//&
            "with a LANGEVIN scheme.",&
            usage="LANGEVIN {logical}",&
            citations=(/VandenCic2006/),&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WW",&
            description="Specifies the height of the gaussian to spawn. Default 0.1 .",&
            usage="WW <REAL>",unit_str='hartree',default_r_val=0.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_HILLS",&
            description="This keyword enables the spawning of the hills. Default .FALSE.",&
            usage="DO_HILLS",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WELL_TEMPERED",&
            description="This keyword enables Well-tempered metadynamics. Default .FALSE.",&
            usage="WELL_TEMPERED",citations=(/BarducBus2008/),&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELTA_T",&
            description="If Well-tempered metaD is used, the temperature parameter "//&
            "must be specified.",usage="DELTA_T <REAL>",&
            unit_str='K',default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WTGAMMA",&
            description="If Well-tempered metaD is used, the gamma parameter "//&
            "must be specified if not DELTA_T.",usage="WTGAMMA <REAL>",&
            default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAGRANGE",&
            description="Specifies whether an extended-lagrangian should be used. Default .FALSE.",&
            usage="LAGRANGE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="step_start_val",&
             description="The starting step value for metadynamics",&
             usage="step_start_val <integer>",default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="nhills_start_val",&
             description="The starting value of previously spawned hills",&
             usage="nhills_start_val <integer>",default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="COLVAR_AVG_TEMPERATURE_RESTART",&
             description="COLVAR average temperature. Only for restarting purposes.",&
             usage="COLVAR_AVG_TEMPERATURE_RESTART 0.0",default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="TAMCSteps",&
             description="Number of sampling points for z",&
             usage="TAMCSteps <integer>",default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword=keyword,name="timestep",&
            description="The length of an integration step for colvars (TAMC only)",&
            usage="timestep <real>",default_r_val=cp_unit_to_cp2k(value=0.5_dp,unit_str="fs",error=error),&
            unit_str="fs",error=error)

       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_metavar_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_multiple_walkers_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="print",&
            description="Controls the printing properties during an metadynamics run",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., &
            error=error)
       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic and summary information during"//&
            " metadynamics.", &
            print_level=low_print_level,each_iter_names=s2a("MD","METADYNAMICS"),&
            each_iter_values=(/1,1/),add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"temperature_colvar",&
            description="Controls the printing of the temperature of COLVARS in an "//&
            "extended lagrangian scheme.", &
            print_level=low_print_level,each_iter_names=s2a("MD","METADYNAMICS"),&
            each_iter_values=(/1,1/),add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COLVAR",&
            description="Controls the printing of COLVAR summary information during"//&
            " metadynamics. When an extended Lagrangian use used, the files"//&
            " contain (in order): colvar value of the extended Lagrangian, "//&
            " instantaneous colvar value, force due to the harmonic term of the extended "//&
            " Lagrangian and the force due to the previously spawned hills,"//&
            " the force due to the walls, the velocities in the extended "//&
            " Lagrangian, the potential of the harmonic term of the"//&
            " Lagrangian, the potential energy of the hills, the potential"//&
            " energy of the walls and the temperature of the extended"//&
            " Lagrangian. When the extended Lagrangian is not used, all"//&
            " related fields are omitted.",&
            print_level=low_print_level,each_iter_names=s2a("MD","METADYNAMICS"),&
            each_iter_values=(/1,1/),add_last=add_last_numeric,filename="COLVAR",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"HILLS",&
            description="Controls the printing of HILLS summary information during"//&
            " metadynamics. The file contains: instantaneous colvar value, width of "//&
            " the spawned gaussian and height of the gaussian. According the value of "//&
            " the EACH keyword this file may not be synchronized with the COLVAR file.", &
            print_level=high_print_level,each_iter_names=s2a("MD","METADYNAMICS"),&
            each_iter_values=(/1,1/),add_last=add_last_numeric,filename="HILLS",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_metadyn_history(subsection, section, error)
    END IF
  END SUBROUTINE create_metadyn_section

! *****************************************************************************
!> \brief creates the multiple walker section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teodoro laino [tlaino] 10.2008
! *****************************************************************************
  SUBROUTINE create_multiple_walkers_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_multiple_walkers_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="MULTIPLE_WALKERS",&
            description="Enables and configures the metadynamics using multiple walkers.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection,keyword)
       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Controls the usage of the multiple walkers in a metadynamics run.",&
            usage="&MULTIPLE_WALKERS T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALKER_ID",&
            description="Sets the walker ID for the local metadynamics run.",&
            usage="WALKER_ID <INTEGER>",type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NUMBER_OF_WALKERS",&
            description="Sets the total number of walkers in the metadynamic run.",&
            usage="NUMBER_OF_WALKERS <INTEGER>",type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALKER_COMM_FREQUENCY",&
            description="Sets the frequency (in unit of spawned hills) for the "//&
            "communication between the several walkers, in order to update the "//&
            "local list of hills with the ones coming from the other walkers",&
            usage="WALKER_COMM_FREQUENCY <INTEGER>",default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALKERS_STATUS",&
            description="Stores the status of the several walkers in the local run.",&
            usage="WALKERS_STATUS <INTEGER> .. <INTEGER>",type_of_var=integer_t,n_var=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="WALKERS_FILE_NAME",&
            description="Specify the basename for the NUMBER_OF_WALKERS files used to "//&
            "communicate between the walkers. Absolute path can be input as well "//&
            "together with the filename. One file will be created for each spawned hill.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specified the communication filename for each walker.",repeats=.TRUE.,&
            usage="{String}", type_of_var=lchar_t, n_var=1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_multiple_walkers_section

! *****************************************************************************
!> \brief creates the alchemical section for free energy evaluation
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teodoro laino [tlaino] 04.2007
! *****************************************************************************
  SUBROUTINE create_ac_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ac_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(keyword)
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="ALCHEMICAL_CHANGE",&
            description="Controls the calculation of delta free energies"//&
            " with the alchemical change method.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="PARAMETER",&
            description="Defines the perturbing parameter of the alchemical change tranformation",&
            usage="PARAMETERS k",  type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHTING_FUNCTION",&
            description="Specifies the weighting function (umbrella potential, part of the mixing function)",&
            usage="WEIGHTING_FUNCTION (E1+E2-LOG(E1/E2))", type_of_var=lchar_t,&
            n_var=1, default_lc_val="0",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_CONV",&
            description="Set the relative tolerance for the convergence of the free energy derivative",&
            usage="EPS_CONV <REAL>",&
            default_r_val=1.0E-2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NEQUIL_STEPS",&
            description="Set the number of equilibration steps, skipped to compute averages",&
            usage="NEQUIL_STEPS <INTEGER>",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_ac_section

! *****************************************************************************
!> \brief creates the umbrella integration section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teodoro laino [tlaino] 01.2007
! *****************************************************************************
  SUBROUTINE create_ui_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ui_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="umbrella_integration",&
            description="Controls the calculation of free energy derivatives"//&
            " with the umbrella integration method.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection)
       CALL create_uvar_conv_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_uvar_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_ui_section


! *****************************************************************************
!> \brief Creates the velocity section
!> \param section the section to create
!> \param metadyn_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_metadyn_history(section, metadyn_section, error)
    TYPE(section_type), POINTER              :: section, metadyn_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_metadyn_history', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="SPAWNED_HILLS_POS",&
            description="The position of the spawned hills during metadynamics."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify the spawned hills",repeats=.TRUE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

       CALL section_create(section,name="SPAWNED_HILLS_SCALE",&
            description="The scales of the spawned hills during metadynamics."//&
            "Used for RESTART. When a scale is zero in one or more "//&
            "directions, the Gaussian hill is assumed to be infinitely wide "//&
            "in those directions. The latter can be used to combine spawned "//&
            "hills from multiple 1D metadynamics runs in one multidimensional "//&
            "metadynamics run.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify the spawned hills",repeats=.TRUE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

       CALL section_create(section,name="SPAWNED_HILLS_HEIGHT",&
            description="The height of the spawned hills during metadynamics."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify the spawned hills",repeats=.TRUE.,&
            usage="{Real}", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

       CALL section_create(section,name="SPAWNED_HILLS_INVDT",&
            description="The inverse of the DELTA_T parameter used for Well-Tempered metadynamics."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify the spawned hills",repeats=.TRUE.,&
            usage="{Real}", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)
       !
       ! Extended Lagrangian
       !
       CALL section_create(section,name="EXT_LAGRANGE_SS0",&
            description="Colvar position within an extended Lagrangian formalism."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specified the positions",repeats=.TRUE.,&
            usage="{Real}", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

       CALL section_create(section,name="EXT_LAGRANGE_VVP",&
            description="Colvar velocities within an extended Lagrangian formalism."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specified the velocities",repeats=.TRUE.,&
            usage="{Real}", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

       CALL section_create(section,name="EXT_LAGRANGE_SS",&
            description="Colvar Theta within an extended Lagrangian formalism."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specified the theta",repeats=.TRUE.,&
            usage="{Real}", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

       CALL section_create(section,name="EXT_LAGRANGE_FS",&
            description="Colvar force within an extended Lagrangian formalism."//&
            "Used for RESTART.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specified the theta",repeats=.TRUE.,&
            usage="{Real}", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(metadyn_section, section, error=error)
       CALL section_release(section,error=error)

    END IF
  END SUBROUTINE create_metadyn_history

! *****************************************************************************
!> \brief creates the metavar section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_metavar_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_metavar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection, wall_section

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="METAVAR",&
            description="This section specify the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., &
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specifies the lambda parameter of the collective variable in the"//&
            " extended lagrangian scheme.",&
            usage="LAMBDA <REAL>",unit_str='internal_cp2k',type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS",&
            description="Specifies the mass parameter of the collective variable in the"//&
            " extended lagrangian scheme.",usage="MASS <REAL>",unit_str='amu',type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GAMMA",&
            description="Specifies the friction term in Langevin integration of the collective variable in the"//&
            " extended lagrangian scheme.",&
            citations=(/VandenCic2006/),&
            usage="GAMMA {real}",type_of_var=real_t,unit_str="fs^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE",&
            variants=(/"WIDTH"/),&
            description="Specifies the scale factor for the following collective variable. The history "//&
            "dependent term has the expression: WW * Sum_{j=1}^{nhills} Prod_{k=1}^{ncolvar} "//&
            "[EXP[-0.5*((ss-ss0(k,j))/SCALE(k))^2]], "//&
            "where ncolvar is the number of defined METAVAR and nhills is the number of spawned hills. ",&
            usage="SCALE <REAL>",type_of_var=real_t,unit_str='internal_cp2k',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COLVAR",&
            description="Specifies the colvar on which to apply metadynamics.",&
            usage="COLVAR {integer}", type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Wall section
       NULLIFY(wall_section, subsection)
       CALL section_create(wall_section,name="WALL",&
            description="Controls the activation of walls on COLVAR during a metadynamic run.",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., &
            error=error)

       CALL keyword_create(keyword, name="TYPE",&
            description="Specify the type of wall",&
            usage=" TYPE (REFLECTIVE|QUADRATIC|QUARTIC|GAUSSIAN|NONE)",&
            enum_c_vals=s2a( "REFLECTIVE","QUADRATIC","QUARTIC","GAUSSIAN","NONE"),&
            enum_desc=s2a("Reflective wall. Colvar velocity is inverted when the colvar is beyond the wall position.",&
                          "Applies a quadratic potential at the wall position.",&
                          "Applies a quartic potential at the wall position.",&
                          "Applies a gaussian potential at the wall position.",&
                          "No walls are applied."),&
            enum_i_vals=(/do_wall_reflective,do_wall_quadratic,do_wall_quartic,do_wall_gaussian,do_wall_none/),&
            default_i_val=do_wall_none,error=error)
       CALL section_add_keyword(wall_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POSITION",&
            description="Specify the value of the colvar for the wall position",&
            usage="POSITION <REAL>",unit_str='internal_cp2k',&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(wall_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Reflective wall
       CALL section_create(subsection,name="REFLECTIVE",&
            description="Parameters controlling the reflective wall",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="DIRECTION",&
            description="Specify the direction of the wall.",&
            usage=" TYPE (WALL_PLUS|WALL_MINUS)",&
            enum_c_vals=s2a( "WALL_PLUS","WALL_MINUS"),&
            enum_desc=s2a("Wall extends from the position towards larger values of COLVAR",&
                          "Wall extends from the position towards smaller values of COLVAR"),&
            enum_i_vals=(/do_wall_p,do_wall_m/),default_i_val=do_wall_p,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(wall_section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Quadratic wall
       CALL section_create(subsection,name="QUADRATIC",&
            description="Parameters controlling the quadratic wall",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="DIRECTION",&
            description="Specify the direction of the wall.",&
            usage=" TYPE (WALL_PLUS|WALL_MINUS)",&
            enum_c_vals=s2a( "WALL_PLUS","WALL_MINUS"),&
            enum_desc=s2a("Wall extends from the position towards larger values of COLVAR",&
                          "Wall extends from the position towards smaller values of COLVAR"),&
            enum_i_vals=(/do_wall_p,do_wall_m/),default_i_val=do_wall_p,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Specify the value of the quadratic potential constant: K*(CV-POS)^2",&
            usage="K <REAL>",unit_str='hartree',&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(wall_section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Quartic wall
       CALL section_create(subsection,name="QUARTIC",&
            description="Parameters controlling the quartic wall",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="DIRECTION",&
            description="Specify the direction of the wall.",&
            usage=" TYPE (WALL_PLUS|WALL_MINUS)",&
            enum_c_vals=s2a( "WALL_PLUS","WALL_MINUS"),&
            enum_desc=s2a("Wall extends from the position towards larger values of COLVAR",&
                          "Wall extends from the position towards smaller values of COLVAR"),&
            enum_i_vals=(/do_wall_p,do_wall_m/),default_i_val=do_wall_p,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Specify the value of the quartic potential constant: K*(CV-(POS+/-(1/K^(1/4))))^4",&
            usage="K <REAL>",unit_str='hartree',&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(wall_section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Gaussian wall
       CALL section_create(subsection,name="GAUSSIAN",&
            description="Parameters controlling the gaussian wall.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="WW",&
            description="Specify the height of the gaussian: WW*e^(-((CV-POS)/sigma)^2)",&
            usage="K <REAL>",unit_str='hartree',&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIGMA",&
            description="Specify the width of the gaussian: WW*e^(-((CV-POS)/sigma)^2)",&
            usage="SIGMA <REAL>",unit_str='internal_cp2k',&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(wall_section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_add_subsection(section, wall_section, error=error)
       CALL section_release(wall_section,error=error)

    END IF

  END SUBROUTINE create_metavar_section

! *****************************************************************************
!> \brief creates the uvar section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_uvar_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_uvar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="UVAR",&
            description="This section specify the nature of the collective variables"//&
            " used in computing the free energy.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., &
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="COLVAR",&
            description="Specifies the colvar used to compute free energy",&
            usage="COLVAR {integer}", type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_uvar_section

! *****************************************************************************
!> \brief  creates the section specifying parameters to control the convergence
!>         of the free energy
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teodoro laino [tlaino] 01.2007
! *****************************************************************************
  SUBROUTINE create_uvar_conv_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_uvar_conv_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="CONVERGENCE_CONTROL",&
            description="This section specify parameters controlling the convergence"//&
            " of the free energy.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., &
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="COARSE_GRAINED_WIDTH",&
            variants=(/"CG_WIDTH"/),&
            description="Width of segments in MD steps to generate the set of"//&
            " coarse grained data, providing a correlation independent data set.",&
            usage="COARSE_GRAINED_WIDTH <INTEGER>", default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_COARSE_GRAINED_WIDTH",&
            variants=(/"MAX_CG_WIDTH"/),&
            description="Max Width of segments in MD steps to generate the set of"//&
            " coarse grained data.",&
            usage="MAX_COARSE_GRAINED_WIDTH <INTEGER>", default_i_val=200,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COARSE_GRAINED_POINTS",&
            variants=(/"CG_POINTS"/),&
            description="Set the minimum amount of coarse grained points to collect"//&
            " before starting the statistical analysis",&
            usage="COARSE_GRAINED_POINTS <INTEGER>", default_i_val=30,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_CONV",&
            description="Set the relative tolerance for the convergence of the collective"//&
            " variable averages used to compute the free energy.",&
            usage="EPS_CONV <REAL>",&
            default_r_val=1.0E-2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K_CONFIDENCE_LIMIT",&
            description="Set the confidence limit for the Mann-Kendall trend test.",&
            usage="K_CONFIDENCE_LIMIT <REAL>",&
            default_r_val=0.90_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SW_CONFIDENCE_LIMIT",&
            description="Set the confidence limit for the Shapiro-Wilks normality test.",&
            usage="SW_CONFIDENCE_LIMIT <REAL>",&
            default_r_val=0.90_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VN_CONFIDENCE_LIMIT",&
            description="Set the confidence limit for the Von Neumann serial correlation test.",&
            usage="VN_CONFIDENCE_LIMIT <REAL>",&
            default_r_val=0.90_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_uvar_conv_section

END MODULE input_cp2k_free_energy
