!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
! *****************************************************************************
  SUBROUTINE atom_code(root_section)
    TYPE(section_vals_type), POINTER         :: root_section

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_code', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, run_type_id
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: atom_section

    failure=.FALSE.

    CALL timeset(routineN,handle)

    logger => cp_get_default_logger()
    NULLIFY(atom_section)
    atom_section => section_vals_get_subs_vals(root_section,"ATOM")

    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log")
    CALL atom_header(iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER")

    CALL atom_test(atom_section)

    CALL section_vals_val_get(atom_section,"RUN_TYPE",i_val=run_type_id)
    SELECT CASE (run_type_id)
      CASE (atom_no_run)
        ! do (almost) nothing
      CASE (atom_energy_run)
        CALL atom_energy_opt(atom_section)
      CASE (atom_basis_run)
        CALL atom_basis_opt(atom_section)
      CASE (atom_pseudo_run)
        CALL atom_pseudo_opt(atom_section)
      CASE default
        CPASSERT(.FALSE.)
    END SELECT

    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log")
    CALL atom_footer(iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER")

    CALL timestop(handle)

  END SUBROUTINE atom_code

! *****************************************************************************
!> \brief ...
!> \param atom_section ...
! *****************************************************************************
  SUBROUTINE atom_test(atom_section)
    TYPE(section_vals_type), POINTER         :: atom_section

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_test', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    INTEGER                                  :: i, ie, z
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL section_vals_val_get(atom_section,"ATOMIC_NUMBER", i_val=z)
    CALL section_vals_val_get(atom_section,"ELEMENT", c_val=elem)

    ie = 0
    DO i=1,nelem
      IF ( ptable(i)%symbol == elem ) THEN
        ie = i
        EXIT
      END IF
    END DO
    IF (ie /= z) THEN
      IF ( ie /= 1 .AND. z /= 1 ) THEN
        CPASSERT(.FALSE.)
      END IF
    END IF

  END SUBROUTINE atom_test

! *****************************************************************************

END MODULE atom
