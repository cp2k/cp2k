!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2002 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_energy_types [1.0] *
!!
!!   NAME
!!     qs_energy_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (13.06.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (11.08.2002) exchange and correlation energy now in exc
!!
!!   SOURCE
!******************************************************************************

MODULE qs_energy_types

  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE qs_energy_type
    REAL(KIND = dp) :: core,&
                core_overlap,&
                core_self,&
                exc,&
                hartree,&
                total
  END TYPE qs_energy_type

! *** Public data types ***

  PUBLIC :: qs_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_qs_energy,&
            deallocate_qs_energy,&
            init_qs_energy

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_qs_energy(qs_energy)

!   Purpose: Allocate and/or initialise a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_qs_energy (MODULE qs_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(qs_energy)) THEN
      ALLOCATE (qs_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"qs_energy",8*dp_size)
    END IF

    CALL init_qs_energy(qs_energy)

  END SUBROUTINE allocate_qs_energy

! *****************************************************************************

  SUBROUTINE deallocate_qs_energy(qs_energy)

!   Purpose: Deallocate a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_qs_energy (MODULE qs_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(qs_energy)) THEN
      DEALLOCATE (qs_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"qs_energy")
    ELSE
      CALL stop_program(routine,"The qs_energy pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_qs_energy

! *****************************************************************************

  SUBROUTINE init_qs_energy(qs_energy)

!   Purpose: Initialise a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_qs_energy (MODULE qs_energy_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(qs_energy)) THEN
      qs_energy%core = 0.0_dp
      qs_energy%core_overlap = 0.0_dp
      qs_energy%core_self = 0.0_dp
      qs_energy%exc = 0.0_dp
      qs_energy%hartree = 0.0_dp
      qs_energy%total = 0.0_dp
    ELSE
      CALL stop_program(routine,"The qs_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_qs_energy

! *****************************************************************************

END MODULE qs_energy_types
