!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Split md_ener module from md_environment_type
!>  \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! **************************************************************************************************
MODULE md_ener_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
   TYPE md_ener_type
      INTEGER       :: nfree = 0, nfree_shell = 0
      REAL(KIND=dp) :: constant = 0.0_dp
      REAL(KIND=dp) :: delta_cons = 0.0_dp, delta_epot = 0.0_dp
      REAL(KIND=dp) :: epot = 0.0_dp
      REAL(KIND=dp) :: ekin = 0.0_dp, ekin_qm = 0.0_dp
      REAL(KIND=dp) :: temp_part = 0.0_dp, temp_qm = 0.0_dp
      REAL(KIND=dp) :: temp_baro = 0.0_dp
      REAL(KIND=dp) :: ekin_coefs = 0.0_dp
      REAL(KIND=dp) :: temp_coefs = 0.0_dp
      REAL(KIND=dp) :: ekin_shell = 0.0_dp, temp_shell = 0.0_dp
      REAL(KIND=dp) :: thermostat_part_kin = 0.0_dp, thermostat_part_pot = 0.0_dp
      REAL(KIND=dp) :: thermostat_fast_kin = 0.0_dp, thermostat_fast_pot = 0.0_dp
      REAL(KIND=dp) :: thermostat_slow_kin = 0.0_dp, thermostat_slow_pot = 0.0_dp
      REAL(KIND=dp) :: thermostat_baro_kin = 0.0_dp, thermostat_baro_pot = 0.0_dp
      REAL(KIND=dp) :: thermostat_coef_kin = 0.0_dp, thermostat_coef_pot = 0.0_dp
      REAL(KIND=dp) :: thermostat_shell_kin = 0.0_dp, thermostat_shell_pot = 0.0_dp
      REAL(KIND=dp) :: baro_kin = 0.0_dp, baro_pot = 0.0_dp
      REAL(KIND=dp) :: vcom(3) = 0.0_dp, total_mass = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_kind => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: temp_kind => NULL()
      INTEGER, DIMENSION(:), POINTER       :: nfree_kind => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_shell_kind => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: temp_shell_kind => NULL()
      INTEGER, DIMENSION(:), POINTER       :: nfree_shell_kind => NULL()
   END TYPE md_ener_type

! *** Public subroutines and data types ***
   PUBLIC :: create_md_ener, release_md_ener, md_ener_type, zero_md_ener

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_ener_types'

CONTAINS

! **************************************************************************************************
!> \brief retains the given md_ener structure
!> \param md_ener ...
!> \par History
!>      10.2007 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE create_md_ener(md_ener)
      TYPE(md_ener_type), INTENT(OUT)                    :: md_ener

      MARK_USED(md_ener)

   END SUBROUTINE create_md_ener

! **************************************************************************************************
!> \brief releases the given md_ener structure
!> \param md_ener ...
!> \par History
!>      10.2007 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE release_md_ener(md_ener)
      TYPE(md_ener_type), INTENT(INOUT)                  :: md_ener

      IF (ASSOCIATED(md_ener%temp_kind)) THEN
         DEALLOCATE (md_ener%temp_kind)
      END IF
      IF (ASSOCIATED(md_ener%ekin_kind)) THEN
         DEALLOCATE (md_ener%ekin_kind)
      END IF
      IF (ASSOCIATED(md_ener%nfree_kind)) THEN
         DEALLOCATE (md_ener%nfree_kind)
      END IF
      IF (ASSOCIATED(md_ener%temp_shell_kind)) THEN
         DEALLOCATE (md_ener%temp_shell_kind)
      END IF
      IF (ASSOCIATED(md_ener%ekin_shell_kind)) THEN
         DEALLOCATE (md_ener%ekin_shell_kind)
      END IF
      IF (ASSOCIATED(md_ener%nfree_shell_kind)) THEN
         DEALLOCATE (md_ener%nfree_shell_kind)
      END IF

   END SUBROUTINE release_md_ener

! **************************************************************************************************
!> \brief initialize to zero energies and temperatures
!> \param md_ener ...
!> \param tkind ...
!> \param tshell ...
!> \par History
!>      10.2007 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE zero_md_ener(md_ener, tkind, tshell)
      TYPE(md_ener_type), INTENT(INOUT)                  :: md_ener
      LOGICAL, INTENT(IN)                                :: tkind, tshell

      md_ener%ekin = 0.0_dp
      md_ener%temp_part = 0.0_dp
      md_ener%temp_baro = 0.0_dp
      md_ener%ekin_coefs = 0.0_dp
      md_ener%temp_coefs = 0.0_dp
      md_ener%ekin_qm = 0.0_dp
      md_ener%temp_qm = 0.0_dp
      md_ener%ekin_shell = 0.0_dp
      md_ener%temp_shell = 0.0_dp
      md_ener%constant = 0.0_dp
      md_ener%delta_cons = 0.0_dp
      md_ener%delta_epot = 0.0_dp
      md_ener%thermostat_part_kin = 0.0_dp
      md_ener%thermostat_part_pot = 0.0_dp
      md_ener%thermostat_fast_kin = 0.0_dp
      md_ener%thermostat_fast_pot = 0.0_dp
      md_ener%thermostat_slow_kin = 0.0_dp
      md_ener%thermostat_slow_pot = 0.0_dp
      md_ener%thermostat_coef_kin = 0.0_dp
      md_ener%thermostat_coef_pot = 0.0_dp
      md_ener%thermostat_baro_kin = 0.0_dp
      md_ener%thermostat_baro_pot = 0.0_dp
      md_ener%thermostat_shell_kin = 0.0_dp
      md_ener%thermostat_shell_pot = 0.0_dp
      md_ener%baro_kin = 0.0_dp
      md_ener%baro_pot = 0.0_dp
      IF (tkind) THEN
         md_ener%temp_kind = 0.0_dp
         md_ener%ekin_kind = 0.0_dp

         IF (tshell) THEN
            md_ener%temp_shell_kind = 0.0_dp
            md_ener%ekin_shell_kind = 0.0_dp
         END IF
      END IF
      md_ener%vcom(:) = 0.0_dp
      md_ener%total_mass = 0.0_dp
   END SUBROUTINE zero_md_ener

END MODULE md_ener_types
