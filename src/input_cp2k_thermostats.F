!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_thermostats *
!!
!!   NAME
!!     input_cp2k_thermostats
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     teo & fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2005 split input_cp2k into smaller modules [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_thermostats
  USE bibliography,                    ONLY: Bussi2007,&
                                             Nose1984a,&
                                             Nose1984b
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_thermostats'

PUBLIC :: create_thermostat_section,&
          create_coord_section,&
          create_velocity_section,&
          create_mass_section

!***
!****************************************************************************
CONTAINS

!***************************************************************************
!!****f* input_cp2k_thermostats/create_thermostat_section *
!!
!!   NAME
!!     create_thermostat_section
!!
!!   FUNCTION
!!     Specifies parameter for thermostat for constant temperature ensembles
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 09.2007
!!
!!*** **********************************************************************
  SUBROUTINE create_thermostat_section(section, coupled_thermostat, error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN), OPTIONAL            :: coupled_thermostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermostat_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_coupled_thermostat
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: csvr_section, nose_section, &
                                                print_key, region_section

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_coupled_thermostat = .FALSE.
       IF (PRESENT(coupled_thermostat)) my_coupled_thermostat = coupled_thermostat
       NULLIFY(csvr_section, nose_section, print_key, region_section)
       
       CALL section_create(section,name="THERMOSTAT",&
            description="Specify thermostat type and parameters controlling the thermostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.NOT.my_coupled_thermostat,&
            error=error)
       NULLIFY(keyword)       

       IF (.NOT.my_coupled_thermostat) THEN
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify the thermostat used for the constant temperature ensembles.",&
               usage="thermostat NOSE",&
               default_i_val=do_thermo_nose,&
               enum_c_vals=s2a("NOSE","CSVR"),&
               enum_i_vals=(/do_thermo_nose,&
               do_thermo_csvr/),&
               enum_desc=s2a("Uses the Nose-Hoover thermostat.",&
               "Uses the canonical sampling through velocity rescaling."),&
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

          CALL create_region_section(region_section,error=error)
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

       CALL create_nose_section(nose_section, my_coupled_thermostat, error=error)
       CALL section_add_subsection(section, nose_section, error=error)
       CALL section_release(nose_section,error=error)
       
       CALL create_csvr_section(csvr_section, my_coupled_thermostat, error=error)
       CALL section_add_subsection(section, csvr_section, error=error)
       CALL section_release(csvr_section,error=error)

       CALL cp_print_key_section_create(print_key,"thermostat_info",&
            description="Controls output information of the corresponding thermostat.", &
            print_level=low_print_level, common_iter_levels=1,&
            filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)               

    END IF
  END SUBROUTINE create_thermostat_section

!***************************************************************************
!!****f* input_cp2k_thermostats/create_region_section *
!!
!!   NAME
!!     create_region_section
!!
!!   FUNCTION
!!     Creates a section to arbitrary define a region to thermostat
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_region_section(section, error)
    TYPE(section_type), POINTER              :: section
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
            " on which applying the thermostats.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms to thermostat.",&
            usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANGE",&
            description="Specifies the starting and the ending index of atoms range you wish to thermostat.",&
            usage="RANGE {integer} {integer}", required=.TRUE.,repeats=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecules to thermostat",&
            usage="MOLNAME WAT MEOH", required=.TRUE., repeats=.TRUE.,&
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
            default_i_val=do_constr_none,repeats=.FALSE.,required=.TRUE.,error=error)
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
            default_i_val=do_constr_none,repeats=.FALSE.,required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
    END IF
  END SUBROUTINE create_region_section

!***************************************************************************
!!****f* input_cp2k_thermostats/create_nose_section *
!!
!!   NAME
!!     create_nose_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the ewald section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     gloria
!!
!!*** **********************************************************************
  SUBROUTINE create_nose_section(section, coupled_thermostat, error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN)                      :: coupled_thermostat
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
            citations=(/Nose1984a,Nose1984b/),required=.NOT. coupled_thermostat,error=error)

       NULLIFY(keyword,subsection)
       CALL keyword_create(keyword, name="length",&
            description="length of the Nose-Hoover chain", usage="length integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Yoshida",&
            description="order of the yoshida integretor usef for the thermostat",&
            usage="Yoshida integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timecon",&
            description="timeconstant of the thermostat chain",&
            usage="timeconstant real", &
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

!***************************************************************************
!!****f* input_cp2k_thermostats/create_csvr_section *
!!
!!   NAME
!!     create_csvr_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the ewald section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 09.2007
!!
!!*** **********************************************************************

  SUBROUTINE create_csvr_section(section, coupled_thermostat, error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN)                      :: coupled_thermostat
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
            description="Paramameters of the canonical sampling through velocity rescaling thermostat.",&
            citations=(/Bussi2007/),required=.NOT.coupled_thermostat,error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="timecon",&
            description="timeconstant of the thermostat.",&
            usage="timeconstant real", &
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
!***************************************************************************

!***************************************************************************
!!****f* input_cp2k_thermostats/create_thermo_energy_section [1.0] *
!!
!!   NAME
!!     create_thermo_energy_section
!!
!!   FUNCTION
!!     Creates the thermostat energy restarting section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify an initial thermostat energy  for CSVR thermostat.",repeats=.TRUE.,&
            unit_str="internal_cp2k",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_thermo_energy_section
!***************************************************************************

!***************************************************************************
!!****f* input_cp2k_thermostats/create_force_section [1.0] *
!!
!!   NAME
!!     create_force_section
!!
!!   FUNCTION
!!     Creates the mass section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify masses of the system",repeats=.FALSE.,&
            usage="<REAL> <REAL> .. <REAL>", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_force_section

!***************************************************************************
!!****f* input_cp2k_thermostats/create_mass_section [1.0] *
!!
!!   NAME
!!     create_mass_section
!!
!!   FUNCTION
!!     Creates the mass section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify masses of the system",repeats=.FALSE.,&
            usage="<REAL> <REAL> .. <REAL>", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mass_section

!***************************************************************************
!!****f* input_cp2k_thermostats/create_velocity_section [1.0] *
!!
!!   NAME
!!     create_velocity_section
!!
!!   FUNCTION
!!     Creates the velocity section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify velocities of the system",repeats=.TRUE.,&
            usage="<REAL> <REAL> .. <REAL>", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_velocity_section

!***************************************************************************
!!****f* input_cp2k_thermostats/create_coord_section [1.0] *
!!
!!   NAME
!!     create_coord_section
!!
!!   FUNCTION
!!     Creates the coord section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify positions of the system",repeats=.TRUE.,&
            usage="<REAL> <REAL> .. <REAL>", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_coord_section

END MODULE input_cp2k_thermostats
