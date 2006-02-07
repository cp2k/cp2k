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
!!     TL  (25.05.2004) qmmm energy
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_types'

  TYPE qs_energy_type
    REAL(KIND = dp) :: core,&
                core_overlap,&
                core_self,&
                exc,&
                exc1,&
                hartree,&
                hartree_1c,&
                hartree_2c,&
                hartree_3c,&
                kg_gpw_ekin_mol,&
                qmmm_el,&
                qmmm_nu,&
                mulliken,&
                total
  END TYPE qs_energy_type

! *** Public data types ***

  PUBLIC :: qs_energy_type

! *** Public subroutines ***

  PUBLIC :: allocate_qs_energy,&
            deallocate_qs_energy,&
            init_qs_energy,&
            duplicate_qs_energy,&
            copy_qs_energy

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
      qs_energy%qmmm_el = 0.0_dp
      qs_energy%qmmm_nu = 0.0_dp
      qs_energy%exc = 0.0_dp
      qs_energy%exc1 = 0.0_dp
      qs_energy%hartree = 0.0_dp
      qs_energy%hartree_1c = 0.0_dp
      qs_energy%hartree_2c = 0.0_dp
      qs_energy%hartree_3c = 0.0_dp
      qs_energy%kg_gpw_ekin_mol = 0.0_dp
      qs_energy%mulliken = 0.0_dp
      qs_energy%total = 0.0_dp
    ELSE
      CALL stop_program(routine,"The qs_energy pointer is not associated "//&
                                "and cannot be initialised")
    END IF

  END SUBROUTINE init_qs_energy

! *****************************************************************************

!!****f* qs_energy_types/duplicate_qs_energy [1.0] *
!!
!!   NAME
!!     duplicate_qs_energy
!!
!!   FUNCTION
!!     Creates a physical duplicate of an qs_energy_type instance.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_energy_in: The qs_energy_type variable to be duplicated.
!!     - qs_energy_out: The allocated and copied output variable.
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     08.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE duplicate_qs_energy(qs_energy_in, qs_energy_out)

    TYPE(qs_energy_type), POINTER            :: qs_energy_in, qs_energy_out

!   ------------------------------------------------------------------------

    IF (ASSOCIATED(qs_energy_out)) THEN
      CALL deallocate_qs_energy(qs_energy_out)
    ELSE
      NULLIFY(qs_energy_out)
    END IF

    CALL allocate_qs_energy(qs_energy_out)

    CALL copy_qs_energy(qs_energy_in, qs_energy_out)
    
  END SUBROUTINE duplicate_qs_energy

! *****************************************************************************

!!****f* qs_energy_types/copy_qs_energy [1.0] *
!!
!!   NAME
!!     copy_qs_energy
!!
!!   FUNCTION
!!     Copies the whole content of a qs_energy_type variable.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_energy_in: The qs_energy_type varibale to be copied.
!!     - qs_energy_out: The output variable which contains the 
!!                      same data as the original.
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     08.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE copy_qs_energy(qs_energy_in, qs_energy_out)

    TYPE(qs_energy_type), POINTER            :: qs_energy_in, qs_energy_out

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE copy_qs_energy (MODULE qs_energy_types)"

!   ------------------------------------------------------------------------

    IF (ASSOCIATED(qs_energy_in) .AND. ASSOCIATED(qs_energy_out)) THEN
      qs_energy_out%core            = qs_energy_in%core
      qs_energy_out%core_overlap    = qs_energy_in%core_overlap
      qs_energy_out%core_self       = qs_energy_in%core_self
      qs_energy_out%qmmm_el         = qs_energy_in%qmmm_el
      qs_energy_out%qmmm_nu         = qs_energy_in%qmmm_nu
      qs_energy_out%exc             = qs_energy_in%exc
      qs_energy_out%exc1            = qs_energy_in%exc1
      qs_energy_out%hartree         = qs_energy_in%hartree
      qs_energy_out%hartree_1c      = qs_energy_in%hartree_1c
      qs_energy_out%hartree_2c      = qs_energy_in%hartree_2c
      qs_energy_out%hartree_3c      = qs_energy_in%hartree_3c
      qs_energy_out%kg_gpw_ekin_mol = qs_energy_in%kg_gpw_ekin_mol
      qs_energy_out%mulliken        = qs_energy_in%mulliken
      qs_energy_out%total           = qs_energy_in%total
    ELSE
      CALL stop_program(routine,"The qs_energy pointer is not associated ")
    END IF

  END SUBROUTINE copy_qs_energy

! *****************************************************************************

END MODULE qs_energy_types
