!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_barostats
  USE barostat_types,                  ONLY: do_clv_x,&
                                             do_clv_xy,&
                                             do_clv_xyz,&
                                             do_clv_xz,&
                                             do_clv_y,&
                                             do_clv_yz,&
                                             do_clv_z
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             high_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
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
  USE input_val_types,                 ONLY: real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "../../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_barostats'

PUBLIC :: create_barostat_section

!***
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param section will contain the coeff section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., error=error)

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

       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="Barostat initial temperature. If not set, the ensemble temperature is used instead.",&
            usage="TEMPERATURE real",type_of_var=real_t,&
            unit_str='K',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMP_TOL",&
            description="Maximum oscillation of the Barostat temperature imposed by rescaling.",&
            usage="TEMP_TOL real",default_r_val=0._dp, &
            unit_str='K',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VIRIAL",&
            description="For NPT_F only: allows the screening of one or more components of the virial in order"//&
            " to relax the cell only along specific cartesian axis",&
            usage="VIRIAL (XYZ | X | Y | Z | XY| XZ | YZ)",&
            enum_c_vals=s2a( "XYZ","X", "Y", "Z", "XY", "XZ", "YZ"),&
            enum_i_vals=(/ do_clv_xyz, do_clv_x, do_clv_y,do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/),&
            default_i_val=do_clv_xyz, error=error)
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

       CALL create_print_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_barostat_section

! *****************************************************************************
!> \brief Creates print section for barostat section
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
            description="Collects all print_keys for barostat",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., error=error)

       CALL cp_print_key_section_create(print_key,"ENERGY",&
            description="Controls the output of kinetic energy, and potential energy "//&
            " of the defined barostat.", print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_print_section

END MODULE input_cp2k_barostats
