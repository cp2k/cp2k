!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE scp_energy_types
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_energy_types'
  PRIVATE

! *****************************************************************************
  TYPE scp_energy_type
    REAL(KIND = dp) :: e_scp_core, e_scp_self, e_scp_total, &
                       e_scp_hartree, e_scp_ks, e_scp_kernel, e_scp_qmmm, &
                       e_scp_ks_self, e_scp_disp
  END TYPE scp_energy_type

! *** Public data types ***

  PUBLIC :: scp_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_scp_energy,&
            deallocate_scp_energy

CONTAINS

! *****************************************************************************
  SUBROUTINE allocate_scp_energy(scp_energy)

    TYPE(scp_energy_type), POINTER           :: scp_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_scp_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (.NOT.ASSOCIATED(scp_energy)) THEN
      ALLOCATE (scp_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "scp_energy",0)
    END IF

    CALL init_scp_energy(scp_energy)

  END SUBROUTINE allocate_scp_energy

! *****************************************************************************
  SUBROUTINE deallocate_scp_energy(scp_energy)

    TYPE(scp_energy_type), POINTER           :: scp_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_scp_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(scp_energy)) THEN
      DEALLOCATE (scp_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"scp_energy")
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The scp_energy pointer is not associated "//&
                        "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_scp_energy

! *****************************************************************************
  SUBROUTINE init_scp_energy(scp_energy)

    TYPE(scp_energy_type), POINTER           :: scp_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'init_scp_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(scp_energy)) THEN
      scp_energy % e_scp_core = 0._dp
      scp_energy % e_scp_self = 0._dp
      scp_energy % e_scp_hartree = 0._dp
      scp_energy % e_scp_ks = 0._dp
      scp_energy % e_scp_ks_self = 0._dp
      scp_energy % e_scp_qmmm = 0._dp
      scp_energy % e_scp_kernel = 0._dp
      scp_energy % e_scp_total = 0._dp
      scp_energy % e_scp_disp = 0._dp
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The scp_energy pointer is not associated "//&
                        "and cannot be initialised")
    END IF

  END SUBROUTINE init_scp_energy

END MODULE scp_energy_types
