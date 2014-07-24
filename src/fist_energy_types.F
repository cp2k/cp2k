!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (11.08.2002) exchange and correlation energy now in exc
!> \author MK (13.06.2002)
! *****************************************************************************
MODULE fist_energy_types

  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_energy_types'

  PRIVATE

! *****************************************************************************
  TYPE fist_energy_type
     REAL ( kind = dp ) :: kin, pot, e_gspace, e_self, e_neut, e_bonded, e_induction
     REAL ( kind = dp ) :: kin_shell, harm_shell
  END TYPE fist_energy_type

! *** Public data types ***

  PUBLIC :: fist_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_fist_energy,&
            deallocate_fist_energy

CONTAINS

! *****************************************************************************
!> \brief   Allocate and/or initialise a Fist energy data structure.
!> \param fist_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_fist_energy(fist_energy)
    TYPE(fist_energy_type), POINTER          :: fist_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_fist_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (.NOT.ASSOCIATED(fist_energy)) THEN
      ALLOCATE (fist_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "fist_energy",8*dp_size)
    END IF

    CALL init_fist_energy(fist_energy)

  END SUBROUTINE allocate_fist_energy

! *****************************************************************************
!> \brief   Deallocate a Fist energy data structure.
!> \param fist_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_fist_energy(fist_energy)
    TYPE(fist_energy_type), POINTER          :: fist_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_fist_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(fist_energy)) THEN
      DEALLOCATE (fist_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"fist_energy")
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The fist_energy pointer is not associated "//&
                        "and cannot be deallocated.")
    END IF

  END SUBROUTINE deallocate_fist_energy

! *****************************************************************************
!> \brief   Initialise a Fist energy data structure.
!> \param fist_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_fist_energy(fist_energy)
    TYPE(fist_energy_type), POINTER          :: fist_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'init_fist_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(fist_energy)) THEN
      fist_energy%kin = 0.0_dp
      fist_energy%pot = 0.0_dp
      fist_energy%e_gspace = 0.0_dp
      fist_energy%e_self = 0.0_dp
      fist_energy%e_neut = 0.0_dp
      fist_energy%e_bonded = 0.0_dp
      fist_energy%e_induction = 0.0_dp
      fist_energy%kin_shell = 0.0_dp
      fist_energy%harm_shell = 0.0_dp
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The fist_energy pointer is not associated "//&
                        "and cannot be initialised.")
    END IF

  END SUBROUTINE init_fist_energy

END MODULE fist_energy_types
