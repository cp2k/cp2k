!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom
  USE atom_basis,                      ONLY: atom_basis_opt
  USE atom_energy,                     ONLY: atom_energy_opt
  USE atom_pseudo,                     ONLY: atom_pseudo_opt
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE header,                          ONLY: atom_footer,&
                                             atom_header
  USE input_constants,                 ONLY: atom_basis_run,&
                                             atom_energy_run,&
                                             atom_no_run,&
                                             atom_pseudo_run
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: atom_code

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom'

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param root_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_code(root_section,error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_code', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, run_type_id
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: atom_section

    failure=.FALSE.

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    NULLIFY(atom_section)
    atom_section => section_vals_get_subs_vals(root_section,"ATOM",error=error)

    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log",error=error)
    CALL atom_header(iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)

    CALL atom_test(atom_section,iw,error)

    CALL section_vals_val_get(atom_section,"RUN_TYPE",i_val=run_type_id,error=error)
    SELECT CASE (run_type_id)
      CASE (atom_no_run)
        ! do (almost) nothing
      CASE (atom_energy_run)
        CALL atom_energy_opt(atom_section,error)
      CASE (atom_basis_run)
        CALL atom_basis_opt(atom_section,error)
      CASE (atom_pseudo_run)
        CALL atom_pseudo_opt(atom_section,error)
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log",error=error)
    CALL atom_footer(iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)

    CALL timestop(handle)

  END SUBROUTINE atom_code

! *****************************************************************************
!> \brief ...
!> \param atom_section ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_test(atom_section,iw,error)
    TYPE(section_vals_type), POINTER         :: atom_section
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_test', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    INTEGER                                  :: i, ie, z
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL section_vals_val_get(atom_section,"ATOMIC_NUMBER", i_val=z, error=error)
    CALL section_vals_val_get(atom_section,"ELEMENT", c_val=elem, error=error)

    ie = 0
    DO i=1,nelem
      IF ( ptable(i)%symbol == elem ) THEN
        ie = i
        EXIT
      END IF
    END DO
    IF (ie /= z) THEN
      IF ( ie /= 1 .AND. z /= 1 ) THEN
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      END IF
    END IF

  END SUBROUTINE atom_test

! *****************************************************************************

END MODULE atom
