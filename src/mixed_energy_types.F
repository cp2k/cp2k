!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/mixed_energy_types [1.0] *
!!
!!   NAME
!!     mixed_energy_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     fschiff 11.06
!!
!!
!!   SOURCE
!******************************************************************************

MODULE mixed_energy_types

  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE mixed_energy_type
     REAL ( kind = dp ) :: pot
     REAL ( kind = dp ) :: kin     
  END TYPE mixed_energy_type

! *** Public data types ***

  PUBLIC :: mixed_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_mixed_energy,&
            deallocate_mixed_energy,&
            init_mixed_energy

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_mixed_energy(mixed_energy)

!   Purpose: Allocate and/or initialise a mixed energy data structure.

!   History: - Creation (11.06,fschiff)

!   ***************************************************************************

    TYPE(mixed_energy_type), POINTER          :: mixed_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_mixed_energy (MODULE mixed_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(mixed_energy)) THEN
      ALLOCATE (mixed_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"mixed_energy",8*dp_size)
    END IF

    CALL init_mixed_energy(mixed_energy)

  END SUBROUTINE allocate_mixed_energy

! *****************************************************************************

  SUBROUTINE deallocate_mixed_energy(mixed_energy)

!   Purpose: Deallocate a mixed energy data structure.

!   History: - Creation (11.06,fschiff)

!   ***************************************************************************

    TYPE(mixed_energy_type), POINTER          :: mixed_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_mixed_energy (MODULE mixed_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(mixed_energy)) THEN
      DEALLOCATE (mixed_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"mixed_energy")
    ELSE
      CALL stop_program(routine,"The mixed_energy pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_mixed_energy

! *****************************************************************************

  SUBROUTINE init_mixed_energy(mixed_energy)

!   Purpose: Initialise a mixed energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(mixed_energy_type), POINTER          :: mixed_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_mixed_energy (MODULE mixed_energy_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(mixed_energy)) THEN
      mixed_energy%pot = 0.0_dp
    ELSE
      CALL stop_program(routine,"The mixed_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_mixed_energy

! *****************************************************************************

END MODULE mixed_energy_types
