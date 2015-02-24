!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! ****************************************************************************
!> \brief   Machine Hardware support for cp2k
!           Interface to the HWLOC library
!> \author  Christiane Pousa Ribeiro
!> \date    2011-24-10
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-24-10
! ****************************************************************************
MODULE ma_errors

  USE machine,                         ONLY: default_output_unit

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_errors'

  PUBLIC :: ma_error_stop

  INTEGER, PARAMETER, PUBLIC  :: ma_error_allocation = 1
  INTEGER, PARAMETER, PUBLIC  :: ma_error_finalize = 2
  INTEGER, PARAMETER, PUBLIC  :: ma_error_option  = 3
  INTEGER, PARAMETER, PUBLIC  :: ma_error_database = 4
  INTEGER, PARAMETER, PUBLIC  :: ma_error_remove = 5

CONTAINS

! ****************************************************************************
!> \brief To show the library erros
!> \param message         The error message
! *****************************************************************************
  SUBROUTINE ma_error (message)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: message

    INTEGER                                  :: unit_num

      unit_num = default_output_unit

      IF (PRESENT(message)) THEN
       WRITE(unit_num,"(T2,A)") "MACHINE| "//TRIM(message)
       WRITE(unit_num,"(T2,A)") "MACHINE| Abnormal program termination"
      ENDIF
      STOP 1
  END SUBROUTINE ma_error

! ****************************************************************************
!> \brief To show the library erros
!> \param error_type ...
! *****************************************************************************
  SUBROUTINE ma_error_stop (error_type)
    INTEGER                                  :: error_type

      SELECT CASE (error_type)
       CASE (ma_error_allocation)
        CALL ma_error("Memory allocation error")
       CASE (ma_error_finalize)
        CALL ma_error("Finalization of arch_lib can not be performed")
       CASE (ma_error_database)
        CALL ma_error("Database cannot be read")
       CASE (ma_error_option)
        CALL ma_error("Option not valid")
       CASE (ma_error_remove)
        CALL ma_error("Can't remove the file")
       END SELECT
      STOP 1
  END SUBROUTINE ma_error_stop



END MODULE ma_errors
