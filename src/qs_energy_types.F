!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (11.08.2002) exchange and correlation energy now in exc
!>      TL  (25.05.2004) qmmm energy
!> \author MK (13.06.2002)
! **************************************************************************************************
MODULE qs_energy_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_types'

   TYPE qs_energy_type
      REAL(KIND=dp) :: core, &
                       core_overlap, &
                       core_overlap0, &
                       core_self, &
                       repulsive, &
                       dispersion, &
                       gcp, &
                       ex, &
                       exc, &
                       exc_aux_fit, &
                       exc1, &
                       exc1_aux_fit, & !the GAPW contribution to ADMM
                       hartree, &
                       hartree_1c, &
                       e_hartree, & ! contains the hartree energy of electrons only
                       ! computed if requested
                       image_charge, &
                       qmmm_el, &
                       qmmm_nu, &
                       mulliken, &
                       cdft, &
                       ee, &
                       ee_core, &
                       efield, &
                       efield_core, &
                       s2_restraint, &
                       dft_plus_u, & ! DFT+U energy contribution
                       sccs_hartree, & ! SCCS Hartree energy of solute + solvent
                       sccs_pol, & ! SCCS polarisation energy
                       sccs_mpc, & ! Makov-Payne energy correction for charged supercells
                       sccs_cav, & ! SCCS cavitation energy
                       sccs_dis, & ! SCCS dispersion free energy
                       sccs_rep, & ! SCCS repulsion free energy
                       ktS, & ! electronic entropic contribution
                       efermi, & ! Fermi energy
                       dftb3, & ! DFTB 3rd order correction
                       nonscf_correction, & ! e.g. Harris correction
                       mp2, &
                       ! single excitations correction for all
                       ! non-scf orbital(density) corrections
                       ! for example, almo delocalization corrction
                       singles_corr, &
                       ! excitation energy
                       excited_state, &
                       total, &
                       tot_old, &
                       kinetic, & !total kinetic energy [rk]
                       s_square, &
                       surf_dipole, &
                       embed_corr, & ! correction for embedding potential
                       xtb_xb_inter, & ! correction for halogen bonding within GFN1-xTB
                       xtb_nonbonded ! correction for nonbonded interactions within GFN1-xTB
      REAL(KIND=dp), DIMENSION(:), POINTER :: ddapc_restraint
   END TYPE qs_energy_type

   ! Public data types

   PUBLIC :: qs_energy_type

   ! Public subroutines

   PUBLIC :: allocate_qs_energy, &
             deallocate_qs_energy, &
             init_qs_energy

CONTAINS

! **************************************************************************************************
!> \brief   Allocate and/or initialise a Quickstep energy data structure.
!> \param qs_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_qs_energy(qs_energy)

      TYPE(qs_energy_type), POINTER                      :: qs_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_qs_energy'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      IF (.NOT. ASSOCIATED(qs_energy)) THEN
         ALLOCATE (qs_energy)
         NULLIFY (qs_energy%ddapc_restraint)
      END IF

      CALL init_qs_energy(qs_energy)
      CALL timestop(handle)

   END SUBROUTINE allocate_qs_energy

! **************************************************************************************************
!> \brief   Deallocate a Quickstep energy data structure.
!> \param qs_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_qs_energy(qs_energy)
      TYPE(qs_energy_type), POINTER                      :: qs_energy

      IF (ASSOCIATED(qs_energy)) THEN
         IF (ASSOCIATED(qs_energy%ddapc_restraint)) THEN
            DEALLOCATE (qs_energy%ddapc_restraint)
         END IF
         DEALLOCATE (qs_energy)
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The qs_energy pointer is not associated "// &
                       "and cannot be deallocated")
      END IF

   END SUBROUTINE deallocate_qs_energy

! **************************************************************************************************
!> \brief   Initialise a Quickstep energy data structure.
!> \param qs_energy ...
!> \date    13.06.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_qs_energy(qs_energy)

      TYPE(qs_energy_type), POINTER                      :: qs_energy

      IF (ASSOCIATED(qs_energy)) THEN
         qs_energy%core = 0.0_dp
         qs_energy%core_overlap = 0.0_dp
         qs_energy%core_overlap0 = 0.0_dp
         qs_energy%core_self = 0.0_dp
         qs_energy%repulsive = 0.0_dp
         qs_energy%dispersion = 0.0_dp
         qs_energy%gcp = 0.0_dp
         qs_energy%qmmm_el = 0.0_dp
         qs_energy%qmmm_nu = 0.0_dp
         qs_energy%ex = 0.0_dp
         qs_energy%exc = 0.0_dp
         qs_energy%exc_aux_fit = 0.0_dp
         qs_energy%exc1 = 0.0_dp
         qs_energy%exc1_aux_fit = 0.0_dp
         qs_energy%e_hartree = 0.0_dp
         qs_energy%hartree = 0.0_dp
         qs_energy%hartree_1c = 0.0_dp
         qs_energy%image_charge = 0.0_dp
         qs_energy%mulliken = 0.0_dp
         qs_energy%cdft = 0.0_dp
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
         qs_energy%excited_state = 0.0_dp
         qs_energy%total = 0.0_dp
         qs_energy%singles_corr = 0.0_dp
         qs_energy%nonscf_correction = 0.0_dp
         IF (.NOT. ASSOCIATED(qs_energy%ddapc_restraint)) THEN
            ALLOCATE (qs_energy%ddapc_restraint(1))
         END IF
         qs_energy%ddapc_restraint(:) = 0.0_dp
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The qs_energy pointer is not associated "// &
                       "and cannot be initialised")
      END IF

   END SUBROUTINE init_qs_energy

END MODULE qs_energy_types
