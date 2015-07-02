!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_thermostats
  USE bibliography,                    ONLY: Bussi2007,&
                                             Ceriotti2009,&
                                             Ceriotti2009b,&
                                             Jones2011,&
                                             Nose1984a,&
                                             Nose1984b
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       do_constr_atomic, do_constr_molec, do_constr_none, do_region_defined, &
       do_region_global, do_region_massive, do_region_molecule, &
       do_region_none, do_thermo_al, do_thermo_csvr, do_thermo_gle, &
       do_thermo_nose, do_thermo_same_as_part
  USE input_cp2k_subsys,               ONLY: create_rng_section
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
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_thermostats'

PUBLIC :: create_thermostat_section,&
          create_thermo_fast_section,&
          create_thermo_slow_section,&
          create_coord_section,&
          create_region_section,&
          create_velocity_section,&
          create_mass_section,&
          create_gle_section

!***
CONTAINS
! *****************************************************************************
!> \brief Specifies parameter for thermostat for constant temperature ensembles
!> \param section will contain the coeff section
!> \param coupled_thermostat ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
  SUBROUTINE create_thermo_slow_section(section, coupled_thermostat, error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN), OPTIONAL            :: coupled_thermostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermo_slow_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_coupled_thermostat
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: nose_section, region_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_coupled_thermostat = .FALSE.
       IF (PRESENT(coupled_thermostat)) my_coupled_thermostat = coupled_thermostat
       NULLIFY(nose_section, region_section)

       CALL section_create(section,name="THERMOSTAT_SLOW",&
            description="Specify thermostat type and parameters controlling the thermostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., error=error)
       NULLIFY(keyword)

       IF (.NOT.my_coupled_thermostat) THEN
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_nose,&
               enum_c_vals=s2a("NOSE"),&
               enum_i_vals=(/do_thermo_nose/),&
               enum_desc=s2a("Uses only the Nose-Hoover thermostat."),&
               error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)

          CALL keyword_create(keyword, name="REGION",&
               description="Determines the defined region for slow thermostat",&
               usage="REGION (GLOBAL||MOLECULE||MASSIVE||DEFINED||NONE)", &
               enum_c_vals=s2a( "GLOBAL", "MOLECULE", "MASSIVE", "DEFINED", "NONE"),&
               enum_i_vals=(/do_region_global, do_region_molecule,&
                             do_region_massive, do_region_defined, do_region_none/),&
                             default_i_val=do_region_global,error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)

          CALL create_region_section(region_section,"slow thermostat",error=error)
          CALL section_add_subsection(section, region_section, error=error)
          CALL section_release(region_section,error=error)
       ELSE
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_same_as_part,&
               enum_c_vals=s2a("SAME_AS_PARTICLE","NOSE","CSVR"),&
               enum_i_vals=(/do_thermo_same_as_part,do_thermo_nose,do_thermo_csvr/),&
               enum_desc=s2a("Use the same kind of thermostat used for particles.",&
                             "Uses the Nose-Hoover thermostat.",&
                             "Uses the canonical sampling through velocity rescaling."),&
               error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END IF

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(section, nose_section, error=error)
       CALL section_release(nose_section,error=error)

       ! Print Section
!       CALL create_print_section(subsection, error=error)
!       CALL section_add_subsection(section, subsection, error=error)
!       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_thermo_slow_section

! *****************************************************************************
!> \brief Specifies parameter for thermostat for constant temperature ensembles
!> \param section will contain the coeff section
!> \param coupled_thermostat ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
  SUBROUTINE create_thermo_fast_section(section, coupled_thermostat, error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN), OPTIONAL            :: coupled_thermostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermo_fast_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_coupled_thermostat
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: nose_section, region_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_coupled_thermostat = .FALSE.
       IF (PRESENT(coupled_thermostat)) my_coupled_thermostat = coupled_thermostat
       NULLIFY(nose_section, region_section)

       CALL section_create(section,name="THERMOSTAT_FAST",&
            description="Specify thermostat type and parameters controlling the thermostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., error=error)
       NULLIFY(keyword)

       IF (.NOT.my_coupled_thermostat) THEN
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_nose,&
               enum_c_vals=s2a("NOSE"),&
               enum_i_vals=(/do_thermo_nose/),&
               enum_desc=s2a("Uses only the Nose-Hoover thermostat."),&
               error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)

          CALL keyword_create(keyword, name="REGION",&
               description="Determines the defined region for fast thermostat",&
               usage="REGION (GLOBAL||MOLECULE||MASSIVE||DEFINED||NONE)", &
               enum_c_vals=s2a( "GLOBAL", "MOLECULE", "MASSIVE", "DEFINED", "NONE"),&
               enum_i_vals=(/do_region_global, do_region_molecule,&
                             do_region_massive, do_region_defined, do_region_none/),&
                             default_i_val=do_region_global,error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)

          CALL create_region_section(region_section,"fast thermostat",error=error)
          CALL section_add_subsection(section, region_section, error=error)
          CALL section_release(region_section,error=error)
       ELSE
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_same_as_part,&
               enum_c_vals=s2a("SAME_AS_PARTICLE","NOSE","CSVR"),&
               enum_i_vals=(/do_thermo_same_as_part,do_thermo_nose,do_thermo_csvr/),&
               enum_desc=s2a("Use the same kind of thermostat used for particles.",&
                             "Uses the Nose-Hoover thermostat.",&
                             "Uses the canonical sampling through velocity rescaling."),&
               error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END IF

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(section, nose_section, error=error)
       CALL section_release(nose_section,error=error)

       ! Print Section
!       CALL create_print_section(subsection, error=error)
!       CALL section_add_subsection(section, subsection, error=error)
!       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_thermo_fast_section


! *****************************************************************************
!> \brief Specifies parameter for thermostat for constant temperature ensembles
!> \param section will contain the coeff section
!> \param coupled_thermostat ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
  SUBROUTINE create_thermostat_section(section, coupled_thermostat, error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN), OPTIONAL            :: coupled_thermostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermostat_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_coupled_thermostat
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: al_section, csvr_section, &
                                                gle_section, nose_section, &
                                                region_section, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_coupled_thermostat = .FALSE.
       IF (PRESENT(coupled_thermostat)) my_coupled_thermostat = coupled_thermostat
       NULLIFY(csvr_section, gle_section, al_section, nose_section, subsection, region_section)

       CALL section_create(section,name="THERMOSTAT",&
            description="Specify thermostat type and parameters controlling the thermostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., error=error)
       NULLIFY(keyword)

       IF (.NOT.my_coupled_thermostat) THEN
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_nose,&
               enum_c_vals=s2a("NOSE","CSVR","GLE","AD_LANGEVIN"),&
               enum_i_vals=(/do_thermo_nose,&
               do_thermo_csvr,do_thermo_gle,do_thermo_al/),&
               enum_desc=s2a("Uses the Nose-Hoover thermostat.",&
               "Uses the canonical sampling through velocity rescaling.",&
               "Uses GLE thermostat",&
               "Uses adaptive-Langevin thermostat"),&
               error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)

          CALL keyword_create(keyword, name="REGION",&
               description="Determines the region each thermostat is attached to.",&
               usage="REGION (GLOBAL||MOLECULE||MASSIVE||DEFINED||NONE)", &
               enum_c_vals=s2a( "GLOBAL", "MOLECULE", "MASSIVE", "DEFINED", "NONE"),&
               enum_i_vals=(/do_region_global, do_region_molecule,&
                             do_region_massive, do_region_defined, do_region_none/),&
                             default_i_val=do_region_global,error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)

          CALL create_region_section(region_section,"thermostat",error=error)
          CALL section_add_subsection(section, region_section, error=error)
          CALL section_release(region_section,error=error)
       ELSE
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_same_as_part,&
               enum_c_vals=s2a("SAME_AS_PARTICLE","NOSE","CSVR"),&
               enum_i_vals=(/do_thermo_same_as_part,do_thermo_nose,do_thermo_csvr/),&
               enum_desc=s2a("Use the same kind of thermostat used for particles.",&
                             "Uses the Nose-Hoover thermostat.",&
                             "Uses the canonical sampling through velocity rescaling."),&
               error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END IF

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(section, nose_section, error=error)
       CALL section_release(nose_section,error=error)

       CALL create_csvr_section(csvr_section, error=error)
       CALL section_add_subsection(section, csvr_section, error=error)
       CALL section_release(csvr_section,error=error)

       CALL create_gle_section(gle_section, error=error)
       CALL section_add_subsection(section, gle_section, error=error)
       CALL section_release(gle_section,error=error)

       CALL create_al_section(al_section, error=error)
       CALL section_add_subsection(section, al_section, error=error)
       CALL section_release(al_section,error=error)

       ! Print Section
       CALL create_print_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_thermostat_section

! *****************************************************************************
!> \brief Creates print section for thermostat section
!> \param section ...
!> \param error ...
!> \author teo [tlaino] - University of Zurich - 02.2008
! *****************************************************************************
  SUBROUTINE create_print_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(print_key)
       CALL section_create(section,name="PRINT",&
            description="Collects all print_keys for thermostat",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       CALL cp_print_key_section_create(print_key,"THERMOSTAT_INFO",&
            description="Controls output information of the corresponding thermostat.", &
            print_level=low_print_level, common_iter_levels=1,&
            filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TEMPERATURE",&
            description="Controls the output of the temperatures of the regions corresponding to "//&
                        "the present thermostat", &
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ENERGY",&
            description="Controls the output of kinetic energy, and potential energy "//&
            " of the defined thermostat.", print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_print_section

! *****************************************************************************
!> \brief Creates a section to arbitrary define a region to thermostat
!> \param section will contain the coeff section
!> \param label ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_region_section(section, label, error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: label
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_region_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="DEFINE_REGION",&
            description="This section provides the possibility to define arbitrary region "//&
            " for the "//TRIM(label)//".",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., &
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms to thermostat.",&
            usage="LIST {integer} {integer} .. {integer}", repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecules to thermostat",&
            usage="MOLNAME WAT MEOH",  repeats=.TRUE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_SUBSYS",&
            variants=(/"PROTEIN"/),&
            description="In a QM/MM run all  MM atoms are specified as a whole ensemble to be thermostated",&
            usage="MM_SUBSYS (NONE|ATOMIC|MOLECULAR)",&
            enum_c_vals=s2a( "NONE","ATOMIC","MOLECULAR"),&
            enum_i_vals=(/do_constr_none,do_constr_atomic,do_constr_molec/),&
            enum_desc=s2a("Thermostat nothing",&
                          "Only the MM atoms itself",&
                          "The full molecule/residue that contains a MM atom"),&
            default_i_val=do_constr_none,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_SUBSYS",&
            description="In a QM/MM run all QM atoms are specified as a whole ensemble to be thermostated",&
            usage="QM_SUBSYS (NONE|ATOMIC|MOLECULAR)",&
            enum_c_vals=s2a( "NONE","ATOMIC","MOLECULAR"),&
            enum_desc=s2a("Thermostat nothing",&
                          "Only the QM atoms itself",&
                          "The full molecule/residue that contains a QM atom"),&
            enum_i_vals=(/do_constr_none,do_constr_atomic,do_constr_molec/),&
            default_i_val=do_constr_none,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_region_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the ewald section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author gloria
! *****************************************************************************
  SUBROUTINE create_nose_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_nose_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nose",&
            description="paramameters of the Nose Hoover thermostat chain",&
            citations=(/Nose1984a,Nose1984b/),error=error)

       NULLIFY(keyword,subsection)
       CALL keyword_create(keyword, name="length",&
            description="length of the Nose-Hoover chain", usage="length integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Yoshida",&
            description="order of the yoshida integretor used for the thermostat",&
            usage="Yoshida integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timecon",&
            description="timeconstant of the thermostat chain",&
            usage="timecon <REAL>", &
            default_r_val=cp_unit_to_cp2k(1000.0_dp,"fs",error=error),&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="mts", &
            variants=s2a("multiple_time_steps","mult_t_steps"),&
            description="number of multiple timesteps to be used for the NoseHoover chain",&
            usage="mts integer", &
            default_i_val=2,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_coord_section(subsection,"NOSE HOOVER",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_velocity_section(subsection,"NOSE HOOVER",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_mass_section(subsection,"NOSE HOOVER",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_force_section(subsection,"NOSE HOOVER",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_nose_section


! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
!> \param
!> \author
! *****************************************************************************
  SUBROUTINE create_gle_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gle_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="GLE",&
            description="paramameters of the gle thermostat. This section can be generated "//&
                        " from https://epfl-cosmo.github.io/gle4md/index.html?page=matrix ",&
            citations=(/Ceriotti2009,Ceriotti2009b/),error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="NDIM",&
            description="Size of the gle matrix", usage="NDIM 6", &
            default_i_val=5, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="A_SCALE",&
            description="scaling factor for matrix A (for generic matrix A, depends "//&
            "on the characteristic frequency of the system).", usage="A_SCALE 0.5", &
            default_r_val=cp_unit_to_cp2k(1.0_dp,"ps^-1",error=error), unit_str="ps^-1", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="A_LIST",&
            description="A matrix The defaults give optimal sampling for most "//&
            "cristalline and liquid compounds. Generated with the parameters set kv_4-4.a"//&
            "centered on w_0=40 cm^-1.", usage="A_LIST real real real", &
             type_of_var=real_t, unit_str="internal_cp2k",&
             n_var=-1, repeats=.TRUE.,  &
!             default_r_vals=(/ &
!    1.859575861256e+2_dp,  2.726385349840e-1_dp,  1.152610045461e+1_dp, -3.641457826260e+1_dp,  2.317337581602e+2_dp, &
!   -2.780952471206e-1_dp,  8.595159180871e-5_dp,  7.218904801765e-1_dp, -1.984453934386e-1_dp,  4.240925758342e-1_dp, &
!   -1.482580813121e+1_dp, -7.218904801765e-1_dp,  1.359090212128e+0_dp,  5.149889628035e+0_dp, -9.994926845099e+0_dp, &
!   -1.037218912688e+1_dp,  1.984453934386e-1_dp, -5.149889628035e+0_dp,  2.666191089117e+1_dp,  1.150771549531e+1_dp, &
!    2.180134636042e+2_dp, -4.240925758342e-1_dp,  9.994926845099e+0_dp, -1.150771549531e+1_dp,  3.095839456559e+2_dp /), &
             error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="C_LIST",&
            description="C matrix", usage="C_LIST real real real", &
            unit_str="K_e",&
            type_of_var=real_t, n_var=-1, repeats=.TRUE.,  &
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_thermo_energy_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_rng_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_gles_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

   END IF

  END SUBROUTINE create_gle_section

! *****************************************************************************
!> \brief Creates the gles section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_gles_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gles_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="s",&
            description="The s variable for GLE used for restart",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify s variable for GLE thermostat ",repeats=.FALSE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_gles_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the ewald section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
  SUBROUTINE create_csvr_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_csvr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="csvr",&
            description="Parameters of the canonical sampling through velocity rescaling thermostat.",&
            citations=(/Bussi2007/),error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="timecon",&
            description="Time constant of the CSVR thermostat. A small time "//&
            "constant will result in strong thermostatting (useful for "//&
            "initial equilibrations) and a large time constant would be adequate "//&
            "to get weak thermostatting in production runs.",&
            usage="timecon <REAL>", &
            default_r_val=cp_unit_to_cp2k(1000.0_dp,"fs",error=error),&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_thermo_energy_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_rng_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_csvr_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the adaptive langevin section 
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Noam [bernstei] 
! *****************************************************************************
  SUBROUTINE create_al_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_al_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN  
       CALL section_create(section,name="ad_langevin",&
            description="Parameters of the adaptive-Langevin thermostat."//&
            " Known to work with NVT ensemble, but not tested with"//&
            " other ensembles.  Also tested with FIXED_ATOMS constraints, but"//&
            " may not work with other constraints (restraints should be OK, but"//&
            " haven't been well tested)",&
            citations=(/Jones2011/),error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="timecon_nh",&
            description="Time constant of the Nose-Hoover part of the AD_LANGEVIN thermostat. A small time "//&
            "constant will result in strong thermostatting (useful for "//&
            "initial equilibrations) and a large time constant would be adequate "//&
            "to get weak thermostatting in production runs.",&
            usage="timecon_nh <REAL>", &
            default_r_val=cp_unit_to_cp2k(1000.0_dp,"fs",error=error),&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timecon_langevin",&
            description="Time constant of the Langevin part of the AD_LANGEVIN thermostat. A small time "//&
            "constant will result in strong thermostatting (useful for "//&
            "initial equilibrations) and a large time constant would be adequate "//&
            "to get weak thermostatting in production runs.",&
            usage="timecon_langevin <REAL>", &
            default_r_val=cp_unit_to_cp2k(1000.0_dp,"fs",error=error),&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_thermo_chi_mass_section(subsection,"CHI",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_thermo_chi_mass_section(subsection,"MASS",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_al_section

! *****************************************************************************
!> \brief Creates the thermostat chi restarting section
!> \param section the section to create
!> \param sec_name ...
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author teo
! *****************************************************************************
  SUBROUTINE create_thermo_chi_mass_section(section,sec_name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*)                         :: sec_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_thermo_chi_mass_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name=TRIM(sec_name),&
            description="Information to initialize the Ad-Langevin thermostat DOF "//TRIM(sec_name),&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify an initial thermostat DOF "//TRIM(sec_name)//&
                        " for Ad-Langevin thermostat.",repeats=.TRUE.,&
            unit_str="fs^-1",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_thermo_chi_mass_section

! *****************************************************************************
!> \brief Creates the thermostat energy restarting section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_thermo_energy_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermo_energy_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="THERMOSTAT_ENERGY",&
            description="Information to initialize the CSVR thermostat energy.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify an initial thermostat energy for CSVR thermostat.",&
            repeats=.TRUE., unit_str="internal_cp2k",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_thermo_energy_section

! *****************************************************************************
!> \brief Creates the mass section
!> \param section the section to create
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_force_section(section,name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_force_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="force",&
            description="The forces for "//TRIM(name)//" used for restart",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify masses of the system",repeats=.FALSE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_force_section

! *****************************************************************************
!> \brief Creates the mass section
!> \param section the section to create
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_mass_section(section,name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mass_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mass",&
            description="The masses for "//TRIM(name)//" used for restart",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify masses of the system",repeats=.FALSE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mass_section

! *****************************************************************************
!> \brief Creates the velocity section
!> \param section the section to create
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_velocity_section(section,name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_velocity_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="velocity",&
            description="The velocities for "//TRIM(name)//" used for restart",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify velocities of the system",repeats=.TRUE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_velocity_section

! *****************************************************************************
!> \brief Creates the coord section
!> \param section the section to create
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_coord_section(section,name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_coord_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="coord",&
            description="The positions for "//TRIM(name)//" used for restart",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify positions of the system",repeats=.TRUE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_coord_section

END MODULE input_cp2k_thermostats
