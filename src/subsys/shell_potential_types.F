!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author MI (12.01.2007)
! *****************************************************************************
MODULE shell_potential_types

  
  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "../common/cp_common_uses.f90"

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
!> \brief ...
!> \param shell ...
!> \param charge ...
!> \param charge_core ...
!> \param charge_shell ...
!> \param mass_core ...
!> \param mass_shell ...
!> \param k2_spring ...
!> \param k4_spring ...
!> \param max_dist ...
!> \param shell_cutoff ...
! *****************************************************************************
  SUBROUTINE get_shell(shell,charge,charge_core,charge_shell,mass_core,&
                       mass_shell,k2_spring,k4_spring,max_dist,shell_cutoff)

    TYPE(shell_kind_type), POINTER           :: shell
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: charge, charge_core, &
      charge_shell, mass_core, mass_shell, k2_spring, k4_spring, max_dist, &
      shell_cutoff

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
!> \brief ...
!> \param shell ...
! *****************************************************************************
  SUBROUTINE shell_create(shell)

    TYPE(shell_kind_type), POINTER           :: shell

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(shell))
    ALLOCATE(shell)
    shell%ref_count = 1

  END SUBROUTINE shell_create

! *****************************************************************************
!> \brief ...
!> \param shell_list ...
!> \param ndim ...
! *****************************************************************************
  SUBROUTINE shell_p_create(shell_list,ndim)

    TYPE(shell_p_type), DIMENSION(:), &
      POINTER                                :: shell_list
    INTEGER, INTENT(IN)                      :: ndim

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(.NOT.ASSOCIATED(shell_list))
    ALLOCATE(shell_list(ndim))

    DO i = 1,ndim
      NULLIFY (shell_list(i)%shell)
      CALL shell_create(shell_list(i)%shell)
      shell_list(i)%atm_name=''
    END DO

  END SUBROUTINE shell_p_create

! *****************************************************************************
!> \brief ...
!> \param shell ...
! *****************************************************************************
  SUBROUTINE shell_retain(shell)

    TYPE(shell_kind_type), POINTER           :: shell

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(ASSOCIATED(shell))
    CPASSERT(shell%ref_count>0)
    shell%ref_count=shell%ref_count+1

  END SUBROUTINE shell_retain

! *****************************************************************************
!> \brief ...
!> \param shell ...
! *****************************************************************************
  SUBROUTINE shell_release(shell)

    TYPE(shell_kind_type), POINTER           :: shell

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF(ASSOCIATED(shell)) THEN
      CPASSERT(shell%ref_count>0)
      shell%ref_count=shell%ref_count-1
      IF(shell%ref_count==0) THEN
        DEALLOCATE(shell)
      END IF
    END IF
    NULLIFY(shell)

  END SUBROUTINE shell_release

! *****************************************************************************
!> \brief ...
!> \param shell_list ...
! *****************************************************************************
  SUBROUTINE shell_p_release(shell_list)
    TYPE(shell_p_type), DIMENSION(:), &
      POINTER                                :: shell_list

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (ASSOCIATED(shell_list)) THEN
      DO i = 1,SIZE(shell_list)
         CALL shell_release(shell_list(i)%shell)
      END DO
      DEALLOCATE(shell_list)
    END IF

    NULLIFY (shell_list)

  END SUBROUTINE shell_p_release

END MODULE shell_potential_types
