!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for Quickstep NON-SCF run.
!> \par History
!>      - initial setup [JGH, 2024]
!> \author JGH (13.05.2024)
! **************************************************************************************************
MODULE qs_nonscf_utils
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_nonscf_utils'

   PUBLIC :: qs_nonscf_print_summary

CONTAINS

! **************************************************************************************************
!> \brief writes a summary of information after diagonalization
!> \param qs_env ...
!> \param tdiag ...
!> \param nelectron_total ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE qs_nonscf_print_summary(qs_env, tdiag, nelectron_total, iounit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: tdiag
      INTEGER, INTENT(IN)                                :: nelectron_total, iounit

      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho

      IF (iounit > 0) THEN
         CALL get_qs_env(qs_env=qs_env, energy=energy, dft_control=dft_control)
         IF (qs_env%harris_method) THEN
            CPASSERT(.NOT. dft_control%qs_control%gapw)
            CPASSERT(.NOT. dft_control%qs_control%gapw_xc)

            CALL get_qs_env(qs_env=qs_env, rho=rho, qs_charges=qs_charges)
            CALL qs_rho_get(rho, tot_rho_r=tot_rho_r)
            WRITE (UNIT=iounit, FMT="(/,(T3,A,T41,2F20.10))") &
               "Electronic density on regular grids: ", &
               SUM(tot_rho_r), &
               SUM(tot_rho_r) + nelectron_total, &
               "Core density on regular grids:", &
               qs_charges%total_rho_core_rspace, &
               qs_charges%total_rho_core_rspace - REAL(nelectron_total + dft_control%charge, dp)
            WRITE (UNIT=iounit, FMT="(T3,A,T41,F20.10)") &
               "Total charge density on r-space grids:     ", &
               SUM(tot_rho_r) + &
               qs_charges%total_rho_core_rspace, &
               "Total charge density g-space grids:     ", &
               qs_charges%total_rho_gspace

            WRITE (UNIT=iounit, FMT="(/,T2,A,T40,A,F10.2,T61,F20.10)") &
               "Diagonalization", "Time:", tdiag, energy%band

            WRITE (UNIT=iounit, FMT="((T3,A,T56,F25.14))") &
               "Core Hamiltonian energy:                       ", energy%core, &
               "Overlap energy of the core charge distribution:", energy%core_overlap, &
               "Self energy of the core charge distribution:   ", energy%core_self, &
               "Hartree energy:                                ", energy%hartree, &
               "Exchange-correlation energy:                   ", energy%exc
            IF (energy%dispersion /= 0.0_dp) &
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "Dispersion energy:                             ", energy%dispersion
            IF (energy%gcp /= 0.0_dp) &
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "gCP energy:                                    ", energy%gcp
            IF (energy%efield /= 0.0_dp) &
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "Electric field interaction energy:          ", energy%efield

         ELSEIF (dft_control%qs_control%semi_empirical) THEN
            CPABORT("NONSCF not available")
         ELSEIF (dft_control%qs_control%dftb) THEN
            CPASSERT(energy%dftb3 == 0.0_dp)
            energy%total = energy%total + energy%band + energy%qmmm_el
            WRITE (UNIT=iounit, FMT="(/,T2,A,T40,A,F10.2,T61,F20.10)") &
               "Diagonalization", "Time:", tdiag, energy%total
            WRITE (UNIT=iounit, FMT="((T3,A,T56,F25.14))") &
               "Core Hamiltonian energy:                       ", energy%core, &
               "Repulsive potential energy:                    ", energy%repulsive, &
               "Dispersion energy:                             ", energy%dispersion
            IF (energy%efield /= 0.0_dp) &
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "Electric field interaction energy:          ", energy%efield
         ELSEIF (dft_control%qs_control%xtb) THEN
            energy%total = energy%total + energy%band + energy%qmmm_el
            WRITE (UNIT=iounit, FMT="(/,T2,A,T40,A,F10.2,T61,F20.10)") &
               "Diagonalization", "Time:", tdiag, energy%total
            CPASSERT(dft_control%qs_control%xtb_control%gfn_type == 0)
            WRITE (UNIT=iounit, FMT="((T3,A,T56,F25.14))") &
               "Core Hamiltonian energy:                       ", energy%core, &
               "Repulsive potential energy:                    ", energy%repulsive, &
               "SRB Correction energy:                         ", energy%srb, &
               "Charge equilibration energy:                   ", energy%eeq, &
               "Dispersion energy:                             ", energy%dispersion
            IF (dft_control%qs_control%xtb_control%do_nonbonded) &
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "Correction for nonbonded interactions:         ", energy%xtb_nonbonded
            IF (energy%efield /= 0.0_dp) &
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "Electric field interaction energy:          ", energy%efield
         ELSE
            CPABORT("NONSCF not available")
         END IF
         IF (dft_control%smear) THEN
            WRITE (UNIT=iounit, FMT="((T3,A,T56,F25.14))") &
               "Electronic entropic energy:", energy%kTS
            WRITE (UNIT=iounit, FMT="((T3,A,T56,F25.14))") &
               "Fermi energy:", energy%efermi
         END IF
         IF (energy%qmmm_el /= 0.0_dp) THEN
            WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
            IF (qs_env%qmmm_env_qm%image_charge) THEN
               WRITE (UNIT=iounit, FMT="(T3,A,T56,F25.14)") &
                  "QM/MM image charge energy:                ", energy%image_charge
            END IF
         END IF

         WRITE (UNIT=iounit, FMT="(/,(T3,A,T56,F25.14))") &
            "Total energy:                                  ", energy%total

         CALL m_flush(iounit)
      END IF

   END SUBROUTINE qs_nonscf_print_summary

END MODULE qs_nonscf_utils
