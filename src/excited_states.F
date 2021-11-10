!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for total energy and forces of excited states
!> \par History
!>       01.2020 created
!> \author JGH
! **************************************************************************************************
MODULE excited_states
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE ex_property_calculation,         ONLY: ex_properties
   USE exstates_types,                  ONLY: excited_energy_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              deallocate_qs_force,&
                                              qs_force_type,&
                                              sum_qs_force,&
                                              zero_qs_force
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE response_solver,                 ONLY: response_equation,&
                                              response_force,&
                                              response_force_xtb
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'excited_states'

   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.

   PUBLIC :: excited_state_energy

CONTAINS

! **************************************************************************************************
!> \brief Excited state energy and forces
!>
!> \param qs_env ...
!> \param calculate_forces ...
!> \par History
!>       03.2014 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE excited_state_energy(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'excited_state_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nkind, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: natom_of_kind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: ks_force, lr_force
      TYPE(qs_p_env_type), POINTER                       :: p_env

      CALL timeset(routineN, handle)

      ! Check for energy correction
      IF (qs_env%excited_state) THEN
         logger => cp_get_default_logger()
         IF (logger%para_env%ionode) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         END IF

         CALL get_qs_env(qs_env, exstate_env=ex_env, energy=energy)

         energy%excited_state = ex_env%evalue
         energy%total = energy%total + ex_env%evalue
         IF (calculate_forces) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 27), &
                  " Excited State Forces ", REPEAT("-", 28), "!"
            END IF
            ! prepare force array
            CALL get_qs_env(qs_env, force=ks_force, atomic_kind_set=atomic_kind_set)
            nkind = SIZE(atomic_kind_set)
            ALLOCATE (natom_of_kind(nkind))
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
            NULLIFY (lr_force)
            CALL allocate_qs_force(lr_force, natom_of_kind)
            DEALLOCATE (natom_of_kind)
            CALL zero_qs_force(lr_force)
            CALL set_qs_env(qs_env, force=lr_force)
            !
            NULLIFY (p_env)
            CALL response_equation(qs_env, p_env, ex_env%cpmos, unit_nr)
            !
            CALL get_qs_env(qs_env, dft_control=dft_control)
            IF (dft_control%qs_control%semi_empirical) THEN
               CPABORT("Not available")
            ELSEIF (dft_control%qs_control%dftb) THEN
               CPABORT("Not available")
            ELSEIF (dft_control%qs_control%xtb) THEN
               CALL response_force_xtb(qs_env, p_env, ex_env%matrix_hz, ex_env)
            ELSE
               ! KS-DFT
               CALL response_force(qs_env=qs_env, vh_rspace=ex_env%vh_rspace, &
                                   vxc_rspace=ex_env%vxc_rspace, vtau_rspace=ex_env%vtau_rspace, &
                                   vadmm_rspace=ex_env%vadmm_rspace, matrix_hz=ex_env%matrix_hz, &
                                   matrix_pz=ex_env%matrix_px1, matrix_pz_admm=p_env%p1_admm, &
                                   matrix_wz=p_env%w1, &
                                   p_env=p_env, ex_env=ex_env)
            END IF
            ! add TD and KS forces
            CALL get_qs_env(qs_env, force=lr_force)
            CALL sum_qs_force(ks_force, lr_force)
            CALL set_qs_env(qs_env, force=ks_force)
            CALL deallocate_qs_force(lr_force)
            !
            CALL ex_properties(qs_env, ex_env%matrix_pe, p_env)
            !
            CALL p_env_release(p_env)
            !
         ELSE
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 27), &
                  " Excited State Energy ", REPEAT("-", 28), "!"
               WRITE (unit_nr, '(T2,A,T65,F16.10)') "Excitation Energy [Hartree] ", ex_env%evalue
               WRITE (unit_nr, '(T2,A,T65,F16.10)') "Total Energy [Hartree]", energy%total
            END IF
         END IF

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,A,A)') "!", REPEAT("-", 77), "!"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE excited_state_energy

! **************************************************************************************************

END MODULE excited_states
