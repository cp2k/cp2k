!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2002 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_energy_types [1.0] *
!!
!!   NAME
!!     fist_energy_types
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

MODULE fist_energy_types

  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE fist_energy_type
     REAL ( kind = dp ) :: kin, pot, e_gspace, e_self, e_neut, e_bonded
  END TYPE fist_energy_type

! *** Public data types ***

  PUBLIC :: fist_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_fist_energy,&
            deallocate_fist_energy,&
            init_fist_energy

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_fist_energy(fist_energy)

!   Purpose: Allocate and/or initialise a Fist energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(fist_energy_type), POINTER            :: fist_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_fist_energy (MODULE fist_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(fist_energy)) THEN
      ALLOCATE (fist_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"fist_energy",8*dp_size)
    END IF

    CALL init_fist_energy(fist_energy)

  END SUBROUTINE allocate_fist_energy

! *****************************************************************************

  SUBROUTINE deallocate_fist_energy(fist_energy)

!   Purpose: Deallocate a Fist energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(fist_energy_type), POINTER            :: fist_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_fist_energy (MODULE fist_energy_types)"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(fist_energy)) THEN
      DEALLOCATE (fist_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"fist_energy")
    ELSE
      CALL stop_program(routine,"The fist_energy pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_fist_energy

! *****************************************************************************

  SUBROUTINE init_fist_energy(fist_energy)

!   Purpose: Initialise a Fist energy data structure.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(fist_energy_type), POINTER            :: fist_energy

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_fist_energy (MODULE fist_energy_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(fist_energy)) THEN
      fist_energy%kin = 0.0_dp
      fist_energy%pot = 0.0_dp
      fist_energy%e_gspace = 0.0_dp
      fist_energy%e_self = 0.0_dp
      fist_energy%e_neut = 0.0_dp
      fist_energy%e_bonded = 0.0_dp
    ELSE
      CALL stop_program(routine,"The fist_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_fist_energy

! *****************************************************************************

END MODULE fist_energy_types
