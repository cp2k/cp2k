!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2002 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_energy_types [1.0] *
!!
!!   NAME
!!     kg_energy_types
!!
!!   FUNCTION
!!
!******************************************************************************

MODULE kg_energy_types

  USE kinds, ONLY: dp,&
                   dp_size

  USE global_types, ONLY: global_environment_type
  USE termination,  ONLY: stop_memory,&
                          stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE kg_energy_type
    REAL(KIND = dp) :: core_overlap,& 
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

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_kg_energy(kg_energy)

!   Purpose: Allocate and/or initialise a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(kg_energy_type), POINTER :: kg_energy

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_kg_energy (MODULE kg_energy_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(kg_energy)) THEN
      ALLOCATE (kg_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"kg_energy",8*dp_size)
    END IF

    CALL init_kg_energy(kg_energy)

  END SUBROUTINE allocate_kg_energy

! *****************************************************************************

  SUBROUTINE deallocate_kg_energy(kg_energy)

!   Purpose: Deallocate a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(kg_energy_type), POINTER :: kg_energy

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_kg_energy (MODULE kg_energy_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(kg_energy)) THEN
      DEALLOCATE (kg_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"kg_energy")
    ELSE
      CALL stop_program(routine,"The kg_energy pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_kg_energy

! *****************************************************************************

  SUBROUTINE init_kg_energy(kg_energy)

!   Purpose: Initialise a Quickstep energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(kg_energy_type), POINTER :: kg_energy

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_kg_energy (MODULE kg_energy_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(kg_energy)) THEN
      kg_energy%core_overlap = 0.0_dp
      kg_energy%core_self = 0.0_dp
      kg_energy%exc = 0.0_dp
      kg_energy%hartree = 0.0_dp
      kg_energy%total = 0.0_dp
      kg_energy%pseudo = 0.0_dp
      kg_energy%ppseudo = 0.0_dp
      kg_energy%kin = 0.0_dp
    ELSE
      CALL stop_program(routine,"The kg_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_kg_energy

! *****************************************************************************

END MODULE kg_energy_types
