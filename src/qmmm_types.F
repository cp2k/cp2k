!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Basic container type for QM/MM.
!> \author Ole Schuett
! **************************************************************************************************
MODULE qmmm_types
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE fist_energy_types,               ONLY: fist_energy_type
   USE fist_environment_types,          ONLY: fist_env_get,&
                                              fist_env_release,&
                                              fist_environment_type
   USE kinds,                           ONLY: dp
   USE qmmm_types_low,                  ONLY: qmmm_env_qm_release,&
                                              qmmm_env_qm_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_types'

   PUBLIC :: qmmm_env_type, qmmm_env_release, qmmm_env_get

   TYPE qmmm_env_type
      TYPE(qs_environment_type), POINTER                      :: qs_env => Null()
      TYPE(fist_environment_type), POINTER                    :: fist_env => Null()
      TYPE(qmmm_env_qm_type), POINTER                         :: qm => Null()
   END TYPE qmmm_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
! **************************************************************************************************
   SUBROUTINE qmmm_env_get(qmmm_env, subsys, potential_energy, kinetic_energy)
      TYPE(qmmm_env_type), INTENT(IN)                    :: qmmm_env
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: potential_energy, kinetic_energy

      TYPE(fist_energy_type), POINTER                    :: thermo
      TYPE(qs_energy_type), POINTER                      :: qs_energy

      NULLIFY (qs_energy, thermo)

      IF (PRESENT(kinetic_energy)) THEN
         CALL fist_env_get(qmmm_env%fist_env, thermo=thermo)
         kinetic_energy = thermo%kin
      END IF
      IF (PRESENT(subsys)) THEN
         CALL fist_env_get(qmmm_env%fist_env, subsys=subsys)
      END IF
      IF (PRESENT(potential_energy)) THEN
         ! get the underlying energies from primary subsys.  This is the only subsys
         ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
         CALL fist_env_get(qmmm_env%fist_env, thermo=thermo)
         CALL get_qs_env(qmmm_env%qs_env, energy=qs_energy)
         potential_energy = thermo%pot + qs_energy%total
      END IF
   END SUBROUTINE qmmm_env_get

! **************************************************************************************************
!> \brief releases the given qmmm_env (see doc/ReferenceCounting.html)
!> \param qmmm_env the object to release
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE qmmm_env_release(qmmm_env)
      TYPE(qmmm_env_type), INTENT(INOUT)                 :: qmmm_env

      CALL qs_env_release(qmmm_env%qs_env)
      DEALLOCATE (qmmm_env%qs_env)
      CALL qmmm_env_qm_release(qmmm_env%qm)
      DEALLOCATE (qmmm_env%qm)
      CALL fist_env_release(qmmm_env%fist_env)
      DEALLOCATE (qmmm_env%fist_env)

   END SUBROUTINE qmmm_env_release

END MODULE qmmm_types
