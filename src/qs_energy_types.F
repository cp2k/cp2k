!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (11.08.2002) exchange and correlation energy now in exc
!>      TL  (25.05.2004) qmmm energy
!> \author MK (13.06.2002)
! *****************************************************************************
MODULE qs_energy_types

  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_types'

! *****************************************************************************
  TYPE qs_energy_type
    REAL(KIND = dp) :: core,&
                       core_overlap,&
                       core_self,&
                       repulsive,&
                       dispersion,& 
                       ex,&
                       exc,&
                       exc_aux_fit,&
                       exc1,&
                       hartree,&
                       hartree_1c,&
                       hartree_2c,&
                       hartree_3c,&
                       e_hartree,& ! contains the hartree energy of electrons only
                                   ! computed if requested
                       kg_gpw_ekin_mol,&
                       qmmm_el,&
                       qmmm_nu,&
                       mulliken,&
                       becke,&
                       s2_restraint,&
                       dft_plus_u,& ! DFT+U energy contribution
                       ktS,& ! electronic entropic contribution
                       total
    REAL(KIND = dp), DIMENSION(:), POINTER :: ddapc_restraint
  END TYPE qs_energy_type

  ! Public data types

  PUBLIC :: qs_energy_type

  ! Public subroutines

  PUBLIC :: allocate_qs_energy,&
            deallocate_qs_energy,&
            init_qs_energy,&
            duplicate_qs_energy,&
            copy_qs_energy

CONTAINS

! *****************************************************************************
!> \brief   Allocate and/or initialise a Quickstep energy data structure.
!> \author  MK
!> \date    13.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_qs_energy(qs_energy)

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_qs_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (.NOT.ASSOCIATED(qs_energy)) THEN
      ALLOCATE (qs_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"qs_energy",8*dp_size)
      NULLIFY (qs_energy%ddapc_restraint)
    END IF

    CALL init_qs_energy(qs_energy)

  END SUBROUTINE allocate_qs_energy

! *****************************************************************************
!> \brief   Deallocate a Quickstep energy data structure.
!> \author  MK
!> \date    13.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_qs_energy(qs_energy)
    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_qs_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(qs_energy)) THEN
      IF (ASSOCIATED(qs_energy%ddapc_restraint)) THEN
        DEALLOCATE (qs_energy%ddapc_restraint,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineP,"qs_energy%ddapc_restraint")
      END IF
      DEALLOCATE (qs_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"qs_energy")
    ELSE
      CALL stop_program(routineP,"The qs_energy pointer is not associated "//&
                                 "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_qs_energy

! *****************************************************************************
!> \brief   Initialise a Quickstep energy data structure.
!> \author  MK
!> \date    13.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_qs_energy(qs_energy)

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'init_qs_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(qs_energy)) THEN
      qs_energy%core = 0.0_dp
      qs_energy%core_overlap = 0.0_dp
      qs_energy%core_self = 0.0_dp
      qs_energy%repulsive = 0.0_dp
      qs_energy%dispersion = 0.0_dp
      qs_energy%qmmm_el = 0.0_dp
      qs_energy%qmmm_nu = 0.0_dp
      qs_energy%ex = 0.0_dp
      qs_energy%exc = 0.0_dp
      qs_energy%exc_aux_fit = 0.0_dp
      qs_energy%exc1 = 0.0_dp
      qs_energy%e_hartree = 0.0_dp
      qs_energy%hartree = 0.0_dp
      qs_energy%hartree_1c = 0.0_dp
      qs_energy%hartree_2c = 0.0_dp
      qs_energy%hartree_3c = 0.0_dp
      qs_energy%kg_gpw_ekin_mol = 0.0_dp
      qs_energy%mulliken = 0.0_dp
      qs_energy%becke = 0.0_dp
      qs_energy%s2_restraint = 0.0_dp
      qs_energy%dft_plus_u = 0.0_dp
      qs_energy%kTS = 0.0_dp
      qs_energy%total = 0.0_dp
      IF (.NOT.ASSOCIATED(qs_energy%ddapc_restraint)) THEN
        ALLOCATE (qs_energy%ddapc_restraint(1),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"qs_energy%ddapc_restraint",dp_size)
        END IF
        qs_energy%ddapc_restraint(:) = 0.0_dp
      END IF
    ELSE
      CALL stop_program(routineP,"The qs_energy pointer is not associated "//&
                                 "and cannot be initialised")
    END IF

  END SUBROUTINE init_qs_energy

! *****************************************************************************
!> \brief Creates a physical duplicate of an qs_energy_type instance.
!> \param qs_energy_in The qs_energy_type variable to be duplicated.
!> \param qs_energy_out The allocated and copied output variable.
!> \par History
!>      08.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE duplicate_qs_energy(qs_energy_in, qs_energy_out)

    TYPE(qs_energy_type), POINTER            :: qs_energy_in, qs_energy_out

    IF (ASSOCIATED(qs_energy_out)) THEN
      CALL deallocate_qs_energy(qs_energy_out)
    ELSE
      NULLIFY (qs_energy_out)
    END IF

    CALL allocate_qs_energy(qs_energy_out)

    CALL copy_qs_energy(qs_energy_in,qs_energy_out)
    
  END SUBROUTINE duplicate_qs_energy

! *****************************************************************************
!> \brief Copies the whole content of a qs_energy_type variable.
!> \param qs_energy_in The qs_energy_type varibale to be copied.
!> \param qs_energy_out The output variable which contains the 
!>                       same data as the original.
!> \par History
!>      08.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE copy_qs_energy(qs_energy_in,qs_energy_out)

    TYPE(qs_energy_type), POINTER            :: qs_energy_in, qs_energy_out

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_qs_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(qs_energy_in).AND.ASSOCIATED(qs_energy_out)) THEN
      qs_energy_out%core            = qs_energy_in%core
      qs_energy_out%core_overlap    = qs_energy_in%core_overlap
      qs_energy_out%core_self       = qs_energy_in%core_self
      qs_energy_out%repulsive       = qs_energy_in%repulsive
      qs_energy_out%dispersion      = qs_energy_in%dispersion
      qs_energy_out%qmmm_el         = qs_energy_in%qmmm_el
      qs_energy_out%qmmm_nu         = qs_energy_in%qmmm_nu
      qs_energy_out%ex              = qs_energy_in%ex
      qs_energy_out%exc             = qs_energy_in%exc
      qs_energy_out%exc_aux_fit     = qs_energy_in%exc_aux_fit
      qs_energy_out%exc1            = qs_energy_in%exc1
      qs_energy_out%hartree         = qs_energy_in%hartree
      qs_energy_out%hartree_1c      = qs_energy_in%hartree_1c
      qs_energy_out%hartree_2c      = qs_energy_in%hartree_2c
      qs_energy_out%hartree_3c      = qs_energy_in%hartree_3c
      qs_energy_out%kg_gpw_ekin_mol = qs_energy_in%kg_gpw_ekin_mol
      qs_energy_out%mulliken        = qs_energy_in%mulliken
      qs_energy_out%ddapc_restraint = qs_energy_in%ddapc_restraint
      qs_energy_out%becke           = qs_energy_in%becke
      qs_energy_out%s2_restraint    = qs_energy_in%s2_restraint
      qs_energy_out%dft_plus_u      = qs_energy_in%dft_plus_u
      qs_energy_out%ktS             = qs_energy_in%ktS
      qs_energy_out%total           = qs_energy_in%total
      IF (SIZE(qs_energy_in%ddapc_restraint) /= 1) THEN
        DEALLOCATE (qs_energy_out%ddapc_restraint,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineP,"qs_energy%ddapc_restraint")
        ALLOCATE (qs_energy_out%ddapc_restraint(SIZE(qs_energy_in%ddapc_restraint)),&
                  STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"qs_energy%ddapc_restraint",&
                           SIZE(qs_energy_in%ddapc_restraint)*dp_size)
        END IF
      END IF
    ELSE
      CALL stop_program(routineP,"The qs_energy pointer is not associated ")
    END IF

  END SUBROUTINE copy_qs_energy

END MODULE qs_energy_types
