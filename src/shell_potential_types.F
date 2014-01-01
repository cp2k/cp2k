!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author MI (12.01.2007)
! *****************************************************************************
MODULE shell_potential_types

  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'shell_potential_types'

! *****************************************************************************
!> \brief Define the shell type
! *****************************************************************************
  TYPE shell_kind_type
    INTEGER                                :: ref_count
    REAL(dp)                               :: charge_core, &
                                              charge_shell
    REAL(dp)                               :: mass_core, &
                                              massfrac, &
                                              mass_shell
    REAL(dp)                               :: k2_spring,k4_spring
    REAL(dp)                               :: max_dist
    REAL(dp)                               :: shell_cutoff
  END TYPE shell_kind_type

! *****************************************************************************
  TYPE shell_p_type
    CHARACTER (LEN=default_string_length)   :: atm_name
    TYPE(shell_kind_type), POINTER          :: shell
  END TYPE shell_p_type

! Public subroutines

  PUBLIC :: get_shell, shell_create, shell_p_create, &
            shell_p_release, shell_release, shell_retain

! Public data types

  PUBLIC :: shell_p_type, shell_kind_type

CONTAINS

! *****************************************************************************
  SUBROUTINE get_shell(shell,charge,charge_core,charge_shell,mass_core,&
                       mass_shell,k2_spring,k4_spring,max_dist,shell_cutoff,&
                       error)

    TYPE(shell_kind_type), POINTER           :: shell
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: charge, charge_core, &
      charge_shell, mass_core, mass_shell, k2_spring, k4_spring, max_dist, &
      shell_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_shell', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(shell)) THEN
      IF (PRESENT(charge)) charge = shell%charge_core + shell%charge_shell
      IF (PRESENT(charge_core)) charge_core = shell%charge_core
      IF (PRESENT(charge_shell)) charge_shell = shell%charge_shell
      IF (PRESENT(mass_core)) mass_core = shell%mass_core
      IF (PRESENT(mass_shell)) mass_shell = shell%mass_shell
      IF (PRESENT(k2_spring)) k2_spring = shell%k2_spring
      IF (PRESENT(k4_spring)) k4_spring = shell%k4_spring
      IF (PRESENT(max_dist)) max_dist = shell%max_dist
      IF (PRESENT(shell_cutoff)) shell_cutoff = shell%shell_cutoff
    END IF

  END SUBROUTINE
! *****************************************************************************
  SUBROUTINE shell_create(shell,error)

    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(shell),cp_failure_level,routineP,error,failure)
    ALLOCATE(shell,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    shell%ref_count = 1

  END SUBROUTINE shell_create

! *****************************************************************************
  SUBROUTINE shell_p_create(shell_list,ndim,error)

    TYPE(shell_p_type), DIMENSION(:), &
      POINTER                                :: shell_list
    INTEGER, INTENT(IN)                      :: ndim
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPostcondition(.NOT.ASSOCIATED(shell_list),cp_failure_level,routineP,error,failure)
    ALLOCATE(shell_list(ndim), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i = 1,ndim
      NULLIFY (shell_list(i)%shell)
      CALL shell_create(shell_list(i)%shell,error)
      shell_list(i)%atm_name=''
    END DO

  END SUBROUTINE shell_p_create

! *****************************************************************************
  SUBROUTINE shell_retain(shell,error)

    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(shell),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      CPPrecondition(shell%ref_count>0,cp_failure_level,routineP,error,failure)
      shell%ref_count=shell%ref_count+1
    END IF

  END SUBROUTINE shell_retain

! *****************************************************************************
  SUBROUTINE shell_release(shell,error)

    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF(ASSOCIATED(shell)) THEN
      CPPrecondition(shell%ref_count>0,cp_failure_level,routineP,error,failure)
      shell%ref_count=shell%ref_count-1
      IF(shell%ref_count==0) THEN
        DEALLOCATE(shell,STAT=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
      END IF
    END IF
    NULLIFY(shell)

  END SUBROUTINE shell_release

! *****************************************************************************
  SUBROUTINE shell_p_release(shell_list, error)
    TYPE(shell_p_type), DIMENSION(:), &
      POINTER                                :: shell_list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat

    IF (ASSOCIATED(shell_list)) THEN
      DO i = 1,SIZE(shell_list)
         CALL shell_release(shell_list(i)%shell,error=error)
      END DO
      DEALLOCATE(shell_list,stat=istat)
      CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    END IF

    NULLIFY (shell_list)

  END SUBROUTINE shell_p_release

END MODULE shell_potential_types
