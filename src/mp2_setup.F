!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! *****************************************************************************
MODULE mp2_setup
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE mp2_types,                       ONLY: mp2_method_direct,&
                                             mp2_method_gpw,&
                                             mp2_method_none,&
                                             mp2_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_setup'

  PUBLIC :: read_mp2_section

CONTAINS

  SUBROUTINE read_mp2_section(input,mp2_env,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(mp2_type), INTENT(INOUT)            :: mp2_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_mp2_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: mp2_section

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    failure=.FALSE.

    mp2_section => section_vals_get_subs_vals(input,"DFT%XC%MP2",error=error)


    ! should come from input
    CALL section_vals_val_get(mp2_section,"METHOD",i_val=mp2_env%method,error=error)
    CALL section_vals_val_get(mp2_section,"MEMORY",r_val=mp2_env%mp2_memory,error=error)
    CALL section_vals_val_get(mp2_section,"GROUP_SIZE",i_val=mp2_env%mp2_num_proc,error=error) 
    CALL section_vals_val_get(mp2_section,"DIRECT_CANONICAL%BIG_SEND",l_val=mp2_env%direct_canonical%big_send,error=error)
    CALL section_vals_val_get(mp2_section,"INTERACTION_POTENTIAL%POTENTIAL_TYPE",&
                              i_val=mp2_env%potential_parameter%potential_type, error=error)
    CALL section_vals_val_get(mp2_section,"INTERACTION_POTENTIAL%TRUNCATION_RADIUS",&
                              r_val=mp2_env%potential_parameter%truncation_radius, error=error)
    CALL section_vals_val_get(mp2_section,"INTERACTION_POTENTIAL%TShPSC_DATA",&
                              c_val=mp2_env%potential_parameter%filename, error=error)

    CALL section_vals_val_get(mp2_section,"MP2_GPW%EPS_FILTER",&
                              r_val=mp2_env%mp2_gpw%eps_filter, error=error)
    CALL section_vals_val_get(mp2_section,"MP2_GPW%EPS_GRID",&
                              r_val=mp2_env%mp2_gpw%eps_grid, error=error)
    CALL section_vals_val_get(mp2_section,"MP2_GPW%PRINT_LEVEL",&
                              i_val=mp2_env%mp2_gpw%print_level, error=error)

    ! print some info about the MP2 parameters
    unit_nr = cp_print_key_unit_nr(logger,mp2_section,"MP2_INFO",&
         extension=".mp2Log",error=error)
    IF ((mp2_env%method .NE. mp2_method_none) .AND. unit_nr>0) THEN
       WRITE(unit_nr,'(T2,A)') ""
       SELECT CASE(mp2_env%method)
       CASE(mp2_method_direct)
         WRITE(unit_nr,'(T2,A)') "MP2| using direct canonical MP2"
       CASE(mp2_method_gpw)
         WRITE(unit_nr,'(T2,A)') "MP2| using MP2 GPW style"
       CASE DEFAULT
         CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       WRITE(unit_nr,'(T2,A)') ""
    ENDIF
    CALL cp_print_key_finished_output(unit_nr,logger,mp2_section,&
                                     "MP2_INFO", error=error)

    CALL timestop(handle)

  END SUBROUTINE read_mp2_section

END MODULE mp2_setup


