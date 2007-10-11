!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/thermostat_types [1.0] *
!!
!!   NAME
!!     thermostat_types
!!
!!   FUNCTION
!!     Thermostat structure: module containing thermostat available for MD
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 09.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE thermostat_types
  USE csvr_system_types,               ONLY: csvr_dealloc,&
                                             csvr_init,&
                                             csvr_system_type
  USE extended_system_types,           ONLY: lnhc_dealloc,&
                                             lnhc_init,&
                                             lnhc_parameters_type
  USE input_constants,                 ONLY: do_region_global,&
                                             do_thermo_csvr,&
                                             do_thermo_nose
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: thermostats_type,&
            thermostat_type,&
            allocate_thermostats,&
            release_thermostats,&
            retain_thermostats,&
            create_thermostat_type,&
            release_thermostat_type,&
            thermostat_info_type,&
            release_thermostat_info

!!****f* thermostat_types/thermostats_type [1.0] *
!!
!!   NAME
!!    thermostats_type 
!!
!!   FUNCTION
!!    Define thermostat types
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino] - Teodoro Laino - University of Zurich
!!
!!*** **********************************************************************
  TYPE thermostats_type
     INTEGER :: id_nr, ref_count
     TYPE (thermostat_info_type), POINTER :: thermostat_info_part
     TYPE (thermostat_info_type), POINTER :: thermostat_info_shell
     TYPE(thermostat_type), POINTER :: thermostat_part
     TYPE(thermostat_type), POINTER :: thermostat_coef
     TYPE(thermostat_type), POINTER :: thermostat_shell
     TYPE(thermostat_type), POINTER :: thermostat_baro
  END TYPE thermostats_type

  ! Single thermostat_type
  TYPE thermostat_type
     INTEGER :: id_nr, ref_count, type_of_thermostat
     TYPE (lnhc_parameters_type), POINTER :: nhc
     TYPE (csvr_system_type), POINTER     :: csvr
  END TYPE thermostat_type

  ! Global info type
  TYPE thermostat_info_type
     INTEGER                              :: sum_of_thermostats
     INTEGER                              :: number_of_thermostats
     INTEGER                              :: dis_type
  END TYPE thermostat_info_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_types'
  INTEGER, PRIVATE, SAVE :: last_thermostat_id_nr=0,&
                            last_thermostats_id_nr=0

CONTAINS

!!****f* thermostat_types/allocate_thermostats [1.0] *
!!
!!   NAME
!!    allocate_thermostats 
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE allocate_thermostats(thermostats, error )
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_thermostats', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    check = .NOT.ASSOCIATED(thermostats)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)

    ! Preliminary allocation for thermostats
    ALLOCATE(thermostats, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_thermostats_id_nr = last_thermostats_id_nr + 1
    thermostats%id_nr = last_thermostats_id_nr
    thermostats%ref_count =  1

    ! Thermostats Info
    ALLOCATE(thermostats%thermostat_info_part, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(thermostats%thermostat_info_shell, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(thermostats%thermostat_part)
    NULLIFY(thermostats%thermostat_shell)
    NULLIFY(thermostats%thermostat_coef)
    NULLIFY(thermostats%thermostat_baro)

  END SUBROUTINE allocate_thermostats

!!****f* thermostat_types/retain_thermostats [1.0] *
!!
!!   NAME
!!     retain_thermostats
!!
!!   FUNCTION
!!     retains the  full set of thermostats
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   AUTHOR
!!    Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!! 
!!*** **********************************************************************
  SUBROUTINE retain_thermostats(thermostats,error)
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_thermostats', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(thermostats)) THEN
          CPPrecondition(thermostats%ref_count>0,cp_failure_level,routineP,error,failure)
          thermostats%ref_count=thermostats%ref_count+1
       END IF
    END IF

  END SUBROUTINE retain_thermostats

!!****f* thermostat_types/release_thermostats [1.0] *
!!
!!   NAME
!!     release_thermostats
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE release_thermostats(thermostats, error)
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_thermostats', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    check = ASSOCIATED(thermostats)
    IF (check) THEN
       check = thermostats%ref_count>0
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       thermostats%ref_count=thermostats%ref_count-1
       IF (thermostats%ref_count<1) THEN
          CALL release_thermostat_info(thermostats%thermostat_info_part, error)
          CALL release_thermostat_info(thermostats%thermostat_info_shell, error)
          IF (ASSOCIATED(thermostats%thermostat_part)) THEN
             CALL release_thermostat_type(thermostats%thermostat_part, error)
          END IF
          IF (ASSOCIATED(thermostats%thermostat_shell)) THEN
             CALL release_thermostat_type(thermostats%thermostat_shell, error)
          END IF
          IF (ASSOCIATED(thermostats%thermostat_baro)) THEN
             CALL release_thermostat_type(thermostats%thermostat_baro, error)
          END IF
          IF (ASSOCIATED(thermostats%thermostat_coef)) THEN
             CALL release_thermostat_type(thermostats%thermostat_coef, error)
          END IF
          DEALLOCATE(thermostats, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE release_thermostats

!!****f* thermostat_types/create_thermostat_type [1.0] *
!!
!!   NAME
!!    create_thermostat_type 
!!
!!   FUNCTION
!!    Create a thermostat type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   AUTHOR
!!    Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!! 
!!*** **********************************************************************
  SUBROUTINE create_thermostat_type(thermostat, simpar, section, skip_region, error)  
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(section_vals_type), POINTER         :: section
    LOGICAL, INTENT(IN), OPTIONAL            :: skip_region
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: region, stat
    LOGICAL                                  :: failure, skip_region_loc
    TYPE(section_vals_type), POINTER         :: csvr_section, nose_section

    failure = .FALSE.
    skip_region_loc = .FALSE.
    IF (PRESENT(skip_region)) skip_region_loc = skip_region
    ALLOCATE(thermostat, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_thermostat_id_nr = last_thermostat_id_nr + 1
    thermostat%id_nr = last_thermostat_id_nr
    thermostat%ref_count =  1
    region = do_region_global

    NULLIFY(thermostat%nhc)
    NULLIFY(thermostat%csvr)
    CALL section_vals_val_get(section,"TYPE",i_val=thermostat%type_of_thermostat,error=error)
    IF (.NOT.skip_region_loc) CALL section_vals_val_get(section,"REGION",i_val=region,error=error)
    IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
       nose_section => section_vals_get_subs_vals(section,"NOSE",error=error)
       ALLOCATE(thermostat%nhc, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL lnhc_init(thermostat%nhc, nose_section, error=error)
       thermostat%nhc%region = region
    ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN
       csvr_section => section_vals_get_subs_vals(section,"CSVR",error=error)
       ALLOCATE(thermostat%csvr, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL csvr_init(thermostat%csvr, simpar, csvr_section, error=error)
       thermostat%csvr%region=region
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE create_thermostat_type

!!****f* thermostat_types/release_thermostat_info [1.0] *
!!
!!   NAME
!!     release_thermostat_info
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE release_thermostat_info(thermostat_info, error)
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_thermostat_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure =.FALSE.
    IF (ASSOCIATED(thermostat_info)) THEN
       DEALLOCATE(thermostat_info, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE release_thermostat_info
!!*****

!!****f* thermostat_types/release_thermostat_type [1.0] *
!!
!!   NAME
!!     release_thermostat_type
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE release_thermostat_type(thermostat, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_thermostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    check = ASSOCIATED(thermostat)
    IF (check) THEN
       check = thermostat%ref_count>0
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       thermostat%ref_count=thermostat%ref_count-1
       IF (thermostat%ref_count<1) THEN
          IF (ASSOCIATED(thermostat%nhc)) THEN
             CALL lnhc_dealloc(thermostat%nhc,error=error)
          END IF
          IF (ASSOCIATED(thermostat%csvr)) THEN
             CALL csvr_dealloc(thermostat%csvr,error=error)
          END IF
          DEALLOCATE(thermostat, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE release_thermostat_type
!!*****
!******************************************************************************

END MODULE thermostat_types
!******************************************************************************
