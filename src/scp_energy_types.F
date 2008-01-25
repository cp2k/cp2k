!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_energy_types [1.0] *
!!
!!   NAME
!!     scp_energy_types
!!
!!   FUNCTION
!!
!******************************************************************************

MODULE scp_energy_types

  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE scp_energy_type
    REAL(KIND = dp) :: e_scp_core, e_scp_self, e_scp_total, &
                       e_scp_hartree, e_scp_ks, e_scp_kernel, e_scp_qmmm, &
                       e_scp_ks_self
  END TYPE scp_energy_type

! *** Public data types ***

  PUBLIC :: scp_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_scp_energy,&
            deallocate_scp_energy,&
            init_scp_energy

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_scp_energy(scp_energy)

!   Purpose: Allocate and/or initialise a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(scp_energy_type), POINTER           :: scp_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_scp_energy (MODULE scp_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(scp_energy)) THEN
      ALLOCATE (scp_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"scp_energy",8*dp_size)
    END IF

    CALL init_scp_energy(scp_energy)

  END SUBROUTINE allocate_scp_energy

! *****************************************************************************

  SUBROUTINE deallocate_scp_energy(scp_energy)

!   Purpose: Deallocate a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(scp_energy_type), POINTER           :: scp_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_scp_energy (MODULE scp_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(scp_energy)) THEN
      DEALLOCATE (scp_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"scp_energy")
    ELSE
      CALL stop_program(routine,"The scp_energy pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_scp_energy

! *****************************************************************************

  SUBROUTINE init_scp_energy(scp_energy)

!   Purpose: Initialise a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(scp_energy_type), POINTER           :: scp_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_scp_energy (MODULE scp_energy_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(scp_energy)) THEN
      scp_energy % e_scp_core = 0._dp
      scp_energy % e_scp_self = 0._dp
      scp_energy % e_scp_hartree = 0._dp
      scp_energy % e_scp_ks = 0._dp
      scp_energy % e_scp_ks_self = 0._dp
      scp_energy % e_scp_qmmm = 0._dp
      scp_energy % e_scp_kernel = 0._dp
      scp_energy % e_scp_total = 0._dp
    ELSE
      CALL stop_program(routine,"The scp_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_scp_energy

! *****************************************************************************

END MODULE scp_energy_types
