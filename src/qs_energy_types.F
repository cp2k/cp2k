!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (11.08.2002) exchange and correlation energy now in exc
!>      TL  (25.05.2004) qmmm energy
!> \author MK (13.06.2002)
! *****************************************************************************
MODULE qs_energy_types

  
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_types'

  TYPE qs_energy_type
    REAL(KIND = dp) :: core,&
                       core_overlap,&
                       core_overlap0,&
                       core_self,&
                       repulsive,&
                       dispersion,&
                       ex,&
                       exc,&
                       exc_aux_fit,&
                       exc1,&
                       hartree,&
                       hartree_1c,&
                       e_hartree,& ! contains the hartree energy of electrons only
                                   ! computed if requested
                       image_charge,&
                       qmmm_el,&
                       qmmm_nu,&
                       mulliken,&
                       becke,&
                       ee,&
                       ee_core,&
                       efield,&
                       efield_core,&
                       s2_restraint,&
                       dft_plus_u,& ! DFT+U energy contribution
                       sccs_hartree,& ! SCCS Hartree energy of solute + solvent
                       sccs_pol,& ! SCCS polarisation energy
                       sccs_mpc,& ! Makov-Payne energy correction for charged supercells
                       sccs_cav,& ! SCCS cavitation energy
                       sccs_dis,& ! SCCS dispersion free energy
                       sccs_rep,& ! SCCS repulsion free energy
                       ktS,& ! electronic entropic contribution
                       efermi,& ! Fermi energy
                       dftb3,& ! DFTB 3rd order correction
                       mp2,&
                       total,&
                       tot_old,&
                       kinetic,& !total kinetic energy [rk]
                       s_square,&
                       surf_dipole
    REAL(KIND = dp), DIMENSION(:), POINTER :: ddapc_restraint
  END TYPE qs_energy_type

  ! Public data types

  PUBLIC :: qs_energy_type

  ! Public subroutines

  PUBLIC :: allocate_qs_energy,&
            deallocate_qs_energy,&
            init_qs_energy

CONTAINS

! *****************************************************************************
!> \brief   Allocate and/or initialise a Quickstep energy data structure.
!> \param qs_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_qs_energy(qs_energy)

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_qs_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat

    CALL timeset(routineN,handle)
    IF (.NOT.ASSOCIATED(qs_energy)) THEN
      ALLOCATE (qs_energy,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "qs_energy",int_size)
      NULLIFY (qs_energy%ddapc_restraint)
    END IF

    CALL init_qs_energy(qs_energy)
    CALL timestop(handle)

  END SUBROUTINE allocate_qs_energy

! *****************************************************************************
!> \brief   Deallocate a Quickstep energy data structure.
!> \param qs_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_qs_energy(qs_energy)
    TYPE(qs_energy_type), POINTER            :: qs_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_qs_energy', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(qs_energy)) THEN
      IF (ASSOCIATED(qs_energy%ddapc_restraint)) THEN
        DEALLOCATE (qs_energy%ddapc_restraint)
      END IF
      DEALLOCATE (qs_energy)
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The qs_energy pointer is not associated "//&
                        "and cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_qs_energy

! *****************************************************************************
!> \brief   Initialise a Quickstep energy data structure.
!> \param qs_energy ...
!> \date    13.06.2002
!> \author  MK
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
      qs_energy%core_overlap0 = 0.0_dp
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
      qs_energy%image_charge = 0.0_dp
      qs_energy%mulliken = 0.0_dp
      qs_energy%becke = 0.0_dp
      qs_energy%efield = 0.0_dp
      qs_energy%efield_core = 0.0_dp
      qs_energy%ee = 0.0_dp
      qs_energy%ee_core = 0.0_dp
      qs_energy%s2_restraint = 0.0_dp
      qs_energy%dft_plus_u = 0.0_dp
      qs_energy%sccs_hartree = 0.0_dp
      qs_energy%sccs_pol = 0.0_dp
      qs_energy%sccs_mpc = 0.0_dp
      qs_energy%sccs_cav = 0.0_dp
      qs_energy%sccs_dis = 0.0_dp
      qs_energy%sccs_rep = 0.0_dp
      qs_energy%dftb3 = 0.0_dp
      qs_energy%kTS = 0.0_dp
      qs_energy%mp2 = 0.0_dp
      qs_energy%efermi = 0.0_dp
      qs_energy%kinetic = 0.0_dp
      qs_energy%surf_dipole = 0.0_dp
      qs_energy%total = 0.0_dp
      IF (.NOT.ASSOCIATED(qs_energy%ddapc_restraint)) THEN
        ALLOCATE (qs_energy%ddapc_restraint(1),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,&
                           "qs_energy%ddapc_restraint",dp_size)
        END IF
      END IF
      qs_energy%ddapc_restraint(:) = 0.0_dp
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The qs_energy pointer is not associated "//&
                        "and cannot be initialised")
    END IF

  END SUBROUTINE init_qs_energy

END MODULE qs_energy_types
