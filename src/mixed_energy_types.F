!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author fschiff
!> \date   11.06
! *****************************************************************************
MODULE mixed_energy_types

  
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_energy_types'
  PRIVATE

! *****************************************************************************
  TYPE mixed_energy_type
     REAL ( kind = dp ) :: pot
     REAL ( kind = dp ) :: kin
  END TYPE mixed_energy_type

! *****************************************************************************
  TYPE mixed_force_type
     REAL(KIND=dp), DIMENSION(:,:), POINTER         :: forces
  END TYPE mixed_force_type

! Public data types
  PUBLIC :: mixed_energy_type,&
            mixed_force_type

! Public subroutines
  PUBLIC :: allocate_mixed_energy,&
            deallocate_mixed_energy

CONTAINS

! *****************************************************************************
!> \brief   Allocate and/or initialise a mixed energy data structure.
!> \param mixed_energy ...
!> \date    11.06
!> \author  fschiff
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_mixed_energy(mixed_energy)
    TYPE(mixed_energy_type), POINTER         :: mixed_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_mixed_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.ASSOCIATED(mixed_energy)) THEN
      ALLOCATE (mixed_energy,STAT=stat)
      CPASSERT(stat==0)
    END IF
    CALL init_mixed_energy(mixed_energy)
  END SUBROUTINE allocate_mixed_energy

! *****************************************************************************
!> \brief   Deallocate a mixed energy data structure.
!> \param mixed_energy ...
!> \date    11.06
!> \author  fschiff
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_mixed_energy(mixed_energy)
    TYPE(mixed_energy_type), POINTER         :: mixed_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_mixed_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(mixed_energy)) THEN
      DEALLOCATE (mixed_energy,STAT=stat)
      CPASSERT(stat==0)
    END IF
  END SUBROUTINE deallocate_mixed_energy

! *****************************************************************************
!> \brief ...
!> \param mixed_energy ...
! *****************************************************************************
  SUBROUTINE init_mixed_energy(mixed_energy)
    TYPE(mixed_energy_type), POINTER         :: mixed_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'init_mixed_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(mixed_energy)) THEN
      mixed_energy%pot = 0.0_dp
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The mixed_energy pointer is not associated "//&
                        "and cannot be initialised")
    END IF
  END SUBROUTINE init_mixed_energy

END MODULE mixed_energy_types
