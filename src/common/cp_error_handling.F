!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that contains the routines for error handling
!> \author Ole Schuett
! *****************************************************************************
MODULE cp_error_handling
  USE base_hooks,                      ONLY: cp_abort_hook,&
                                             cp_warn_hook
  USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_abort
  USE termination,                     ONLY: print_message
  USE timings,                         ONLY: print_stack
#include "../base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_error_handling'

  !API public routines
  PUBLIC :: cp_error_handling_setup

  !API (via pointer assignment to hook, PR67982, not meant to be called directly)
  PUBLIC :: cp_abort_handler, cp_warn_handler

CONTAINS

! *****************************************************************************
!> \brief Registers handlers with base_hooks.F
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_error_handling_setup()
    cp_abort_hook => cp_abort_handler
    cp_warn_hook => cp_warn_handler
  END SUBROUTINE cp_error_handling_setup

! *****************************************************************************
!> \brief Abort program with error message
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_abort_handler(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    INTEGER                                  :: unit_nr

    unit_nr = cp_logger_get_default_io_unit()
    IF(unit_nr>0) THEN
       CALL print_message("ABORT in "//TRIM(location)//' :: '//TRIM(ADJUSTL(message)),unit_nr,2,1,1)
       CALL m_flush(unit_nr)
       CALL print_stack(unit_nr)
       CALL m_flush(unit_nr)
    ENDIF

    CALL mp_abort(message)
  END SUBROUTINE cp_abort_handler

! *****************************************************************************
!> \brief Signal a warning
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_warn_handler(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    INTEGER                                  :: unit_nr

    unit_nr = cp_logger_get_default_io_unit()
    IF(unit_nr>0) THEN
       CALL print_message("WARNING in "//TRIM(location)//' :: '//TRIM(ADJUSTL(message)),unit_nr,1,1,1)
       CALL m_flush(unit_nr)
    ENDIF
  END SUBROUTINE cp_warn_handler

END MODULE cp_error_handling
