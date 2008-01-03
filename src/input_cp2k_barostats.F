!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_barostats *
!!
!!   NAME
!!     input_cp2k_barostats
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
MODULE input_cp2k_barostats
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants
  USE input_cp2k_thermostats,          ONLY: create_mass_section,&
                                             create_thermostat_section,&
                                             create_velocity_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_barostats'

PUBLIC :: create_barostat_section

!***
!****************************************************************************
CONTAINS

!***************************************************************************
!!****f* input_cp2k_barostats/create_barostat_section *
!!
!!   NAME
!!     create_barostat_section
!!
!!   FUNCTION
!!     -
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
  SUBROUTINE create_barostat_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_barostat_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection, thermo_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN  
       CALL section_create(section,name="barostat",&
            description="Parameters of barostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,error=error)

       NULLIFY(keyword,subsection,thermo_section)
       CALL keyword_create(keyword, name="PRESSURE",&
            description="Initial pressure",&
            usage="PRESSURE real",&
            default_r_val=0._dp,unit_str='bar',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="TIMECON",&
            description="Barostat time constant",&
            usage="TIMECON real",&
            default_r_val=cp_unit_to_cp2k(1000.0_dp,"fs",error=error),&
            unit_str='fs',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)        

       CALL create_velocity_section(subsection,"BAROSTAT",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_mass_section(subsection,"BAROSTAT",error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_thermostat_section(thermo_section, coupled_thermostat=.TRUE., error=error)
       CALL section_add_subsection(section, thermo_section, error=error)
       CALL section_release(thermo_section,error=error)

    END IF
  END SUBROUTINE create_barostat_section

END MODULE input_cp2k_barostats
