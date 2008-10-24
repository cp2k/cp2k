!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE kg_energy_types

  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_energy_types'
  PRIVATE

! *****************************************************************************
  TYPE kg_energy_type
    REAL(KIND = dp) :: bond_corr,&
                core_overlap,& 
                core_self,&
                exc,&
                hartree,&
                pseudo,&
                ppseudo,&
                kin,&
                total
  END TYPE kg_energy_type

! *** Public data types ***

  PUBLIC :: kg_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_kg_energy,&
            deallocate_kg_energy,&
            init_kg_energy

CONTAINS

! *****************************************************************************
!> \brief   Allocate and/or initialise a KG energy data structure.
!> \author  MK
!> \date    13.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_kg_energy(kg_energy)
    TYPE(kg_energy_type), POINTER            :: kg_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_kg_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (.NOT.ASSOCIATED(kg_energy)) THEN
      ALLOCATE (kg_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"kg_energy",8*dp_size)
    END IF

    CALL init_kg_energy(kg_energy)

  END SUBROUTINE allocate_kg_energy

! *****************************************************************************
!> \brief   Deallocate a KG energy data structure.
!> \author  MK
!> \date    13.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_kg_energy(kg_energy)
    TYPE(kg_energy_type), POINTER            :: kg_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_kg_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(kg_energy)) THEN
      DEALLOCATE (kg_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"kg_energy")
    ELSE
      CALL stop_program(routineP,"The kg_energy pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_kg_energy

! *****************************************************************************
!> \brief   Initialise a KG energy data structure.
!> \author  MK
!> \date    13.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_kg_energy(kg_energy)
    TYPE(kg_energy_type), POINTER            :: kg_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'init_kg_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(kg_energy)) THEN
      kg_energy%bond_corr = 0.0_dp
      kg_energy%core_overlap = 0.0_dp
      kg_energy%core_self = 0.0_dp
      kg_energy%exc = 0.0_dp
      kg_energy%hartree = 0.0_dp
      kg_energy%total = 0.0_dp
      kg_energy%pseudo = 0.0_dp
      kg_energy%ppseudo = 0.0_dp
      kg_energy%kin = 0.0_dp
    ELSE
      CALL stop_program(routineP,"The kg_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_kg_energy

END MODULE kg_energy_types
