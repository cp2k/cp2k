!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Split md_ener module from md_environment_type
!>  \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
MODULE md_ener_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE md_ener_type
     INTEGER       :: id_nr, ref_count
     INTEGER       :: nfree, nfree_qm, nfree_shell
     REAL(KIND=dp) :: constant
     REAL(KIND=dp) :: delta_cons, delta_epot
     REAL(KIND=dp) :: epot
     REAL(KIND=dp) :: ekin, ekin_qm
     REAL(KIND=dp) :: temp_part, temp_qm
     REAL(KIND=dp) :: temp_baro
     REAL(KIND=dp) :: ekin_coefs
     REAL(KIND=dp) :: temp_coefs
     REAL(KIND=dp) :: ekin_shell, temp_shell
     REAL(KIND=dp) :: thermostat_part_kin, thermostat_part_pot
     REAL(KIND=dp) :: thermostat_fast_kin, thermostat_fast_pot
     REAL(KIND=dp) :: thermostat_slow_kin, thermostat_slow_pot
     REAL(KIND=dp) :: thermostat_baro_kin, thermostat_baro_pot
     REAL(KIND=dp) :: thermostat_coef_kin, thermostat_coef_pot
     REAL(KIND=dp) :: thermostat_shell_kin, thermostat_shell_pot
     REAL(KIND=dp) :: baro_kin, baro_pot
     REAL(KIND=dp) :: vcom(3), total_mass
     REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: temp_kind
     INTEGER, DIMENSION(:), POINTER       :: nfree_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_shell_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: temp_shell_kind
     INTEGER, DIMENSION(:), POINTER       :: nfree_shell_kind
  END TYPE md_ener_type

! *** Public subroutines and data types ***
  PUBLIC :: create_md_ener,  release_md_ener,retain_md_ener, md_ener_type, zero_md_ener

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_ener_types'
  INTEGER, SAVE, PRIVATE :: last_md_ener_id=0

CONTAINS

! *****************************************************************************
!> \brief retains the given md_ener structure
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE create_md_ener(md_ener, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_md_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE(md_ener,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       last_md_ener_id=last_md_ener_id+1
       md_ener%id_nr=last_md_ener_id

       md_ener%ref_count = 1

       NULLIFY(md_ener%ekin_kind)
       NULLIFY(md_ener%temp_kind)
       NULLIFY(md_ener%nfree_kind)
       NULLIFY(md_ener%ekin_shell_kind)
       NULLIFY(md_ener%temp_shell_kind)
       NULLIFY(md_ener%nfree_shell_kind)
    END IF
  END SUBROUTINE create_md_ener

! *****************************************************************************
!> \brief retains the given md_ener structure
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE retain_md_ener(md_ener, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_md_ener', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(md_ener%ref_count>0,cp_failure_level,routineP,error,failure)
       md_ener%ref_count=md_ener%ref_count+1
    END IF
  END SUBROUTINE retain_md_ener

! *****************************************************************************
!> \brief releases the given md_ener structure
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE release_md_ener(md_ener, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_md_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(md_ener)) THEN
       CPPrecondition(md_ener%ref_count>0,cp_failure_level,routineP,error,failure)
       md_ener%ref_count=md_ener%ref_count-1
       IF (md_ener%ref_count==0) THEN
          IF(ASSOCIATED(md_ener%temp_kind))THEN
             DEALLOCATE(md_ener%temp_kind,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(md_ener%ekin_kind))THEN
             DEALLOCATE(md_ener%ekin_kind,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(md_ener%nfree_kind))THEN
             DEALLOCATE(md_ener%nfree_kind,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(md_ener%temp_shell_kind))THEN
             DEALLOCATE(md_ener%temp_shell_kind,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(md_ener%ekin_shell_kind))THEN
             DEALLOCATE(md_ener%ekin_shell_kind,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(md_ener%nfree_shell_kind))THEN
             DEALLOCATE(md_ener%nfree_shell_kind,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF

          DEALLOCATE(md_ener,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE release_md_ener

! *****************************************************************************
!> \brief initialize to zero energies and temperatures
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE zero_md_ener(md_ener, tkind, tshell, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    LOGICAL, INTENT(IN)                      :: tkind, tshell
    TYPE(cp_error_type), INTENT(inout)       :: error

    md_ener%ekin        = 0.0_dp
    md_ener%temp_part   = 0.0_dp
    md_ener%temp_baro   = 0.0_dp
    md_ener%ekin_coefs  = 0.0_dp
    md_ener%temp_coefs  = 0.0_dp
    md_ener%ekin_qm     = 0.0_dp
    md_ener%temp_qm     = 0.0_dp
    md_ener%ekin_shell  = 0.0_dp
    md_ener%temp_shell  = 0.0_dp
    md_ener%constant    = 0.0_dp
    md_ener%delta_cons  = 0.0_dp
    md_ener%delta_epot  = 0.0_dp
    md_ener%thermostat_part_kin  = 0.0_dp
    md_ener%thermostat_part_pot  = 0.0_dp
    md_ener%thermostat_fast_kin  = 0.0_dp
    md_ener%thermostat_fast_pot  = 0.0_dp
    md_ener%thermostat_slow_kin  = 0.0_dp
    md_ener%thermostat_slow_pot  = 0.0_dp
    md_ener%thermostat_coef_kin  = 0.0_dp
    md_ener%thermostat_coef_pot  = 0.0_dp
    md_ener%thermostat_baro_kin  = 0.0_dp
    md_ener%thermostat_baro_pot  = 0.0_dp
    md_ener%thermostat_shell_kin  = 0.0_dp
    md_ener%thermostat_shell_pot  = 0.0_dp
    md_ener%baro_kin = 0.0_dp
    md_ener%baro_pot = 0.0_dp
    IF(tkind) THEN
       md_ener%temp_kind  = 0.0_dp
       md_ener%ekin_kind  = 0.0_dp

       IF(tshell) THEN
          md_ener%temp_shell_kind  = 0.0_dp
          md_ener%ekin_shell_kind  = 0.0_dp
       END IF
    END IF
    md_ener%vcom(:) = 0.0_dp
    md_ener%total_mass = 0.0_dp
  END SUBROUTINE zero_md_ener

END MODULE md_ener_types
