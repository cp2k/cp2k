!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Propagation of the electronic degree of freedom in constant-potential MD:
!>        generalized force, momentum kick, position drift and CSVR potentiostat step.
!> \author Nanchen DONGFANG
! **************************************************************************************************
MODULE qs_electrons_utils
   USE csvr_system_utils,               ONLY: rescaling_factor
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: femtoseconds,&
                                              kelvin
   USE qs_electrons_types,              ONLY: qs_electrons_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: qs_electrons_force, &
             qs_electrons_kick, &
             qs_electrons_drift, &
             qs_electrons_thermostat_step

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electrons_utils'
   REAL(KIND=dp), PARAMETER, PRIVATE :: HARTREE_TO_EV = 27.211386245988_dp

CONTAINS

! **************************************************************************************************
!> \brief Generalized electronic force: deviation of the current work function
!>        from the target, driving q toward the requested constant-potential state.
!> \param electrons ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_electrons_force(electrons) RESULT(force)
      TYPE(qs_electrons_type), POINTER                   :: electrons
      REAL(KIND=dp)                                      :: force

      CPASSERT(ASSOCIATED(electrons))
      IF (electrons%control_fermi) THEN
         ! Fermi-level control. The sign is OPPOSITE to work-function control,
         ! and deliberately so: dPhi/dq < 0 while d(e_F)/dq > 0, so the restoring
         ! force has to change sign with the controlled variable:
         !     current_workfunction - target_workfunction   (workfunction mode)
         !     target_ef - current_ef                       (fermi mode)
         force = (electrons%target_fermi_level - electrons%current_fermi_level_eV)/HARTREE_TO_EV
      ELSE
         force = (electrons%current_workfunction - electrons%target_workfunction)/HARTREE_TO_EV
      END IF

   END FUNCTION qs_electrons_force

! **************************************************************************************************
!> \brief Half-kick of the electronic momentum from the current generalized force.
!> \param electrons ...
!> \param dt ...
! **************************************************************************************************
   SUBROUTINE qs_electrons_kick(electrons, dt)
      TYPE(qs_electrons_type), POINTER                   :: electrons
      REAL(KIND=dp), INTENT(IN)                          :: dt

      CPASSERT(ASSOCIATED(electrons))
      ! A fixed mixing weight means q carries no dynamics: freeze it rather than
      ! letting it drift under a force that no longer influences the mixing.
      IF (electrons%lambda_is_fixed) RETURN
      electrons%p = electrons%p + qs_electrons_force(electrons)*dt

   END SUBROUTINE qs_electrons_kick

! **************************************************************************************************
!> \brief Drift of the electronic coordinate (total electron number).
!> \param electrons ...
!> \param dt ...
! **************************************************************************************************
   SUBROUTINE qs_electrons_drift(electrons, dt)
      TYPE(qs_electrons_type), POINTER                   :: electrons
      REAL(KIND=dp), INTENT(IN)                          :: dt

      CPASSERT(ASSOCIATED(electrons))
      IF (electrons%lambda_is_fixed) RETURN
      electrons%q = electrons%q + electrons%p*dt/electrons%mass

   END SUBROUTINE qs_electrons_drift

! **************************************************************************************************
!> \brief Stochastic velocity rescaling ("O" step) of the electronic momentum,
!>        reusing CP2K's own CSVR rescaling_factor (Bussi-Donadio-Parrinello).
!> \param electrons ...
!> \param dt ...
! **************************************************************************************************
   SUBROUTINE qs_electrons_thermostat_step(electrons, dt)
      TYPE(qs_electrons_type), POINTER                   :: electrons
      REAL(KIND=dp), INTENT(IN)                          :: dt

      REAL(KIND=dp)                                      :: alpha, kk, old_kinetic, sigma, taut

      CPASSERT(ASSOCIATED(electrons))

      old_kinetic = 0.5_dp*electrons%p**2/electrons%mass
      kk = old_kinetic
      ! kelvin is the [a.u.] -> [K] conversion factor from physcon, so dividing
      ! a temperature in K by it directly yields k_B*T in Hartree atomic units
      ! -- this is the same idiom CP2K's own thermostat code uses, avoiding any
      ! SI/atomic-unit mismatch with the (atomic-units) kinetic energy above.
      ! sigma is the target kinetic energy ndeg*k_B*T/2 of the single (ndeg = 1)
      ! electronic degree of freedom.
      sigma = 0.5_dp*electrons%temp_potentiostat/kelvin

      ! rescaling_factor expects taut expressed in units of "how often this
      ! routine is called" (see csvr_system_utils.F).  tau_potentiostat is a
      ! plain number in femtoseconds, so the fs -> a.u. conversion is required.
      !
      ! NOTE ON THE ARGUMENT: dt here is the INTERVAL BETWEEN SUCCESSIVE CALLS,
      ! not the MD timestep.  It used to be the MD timestep with a hard-wired
      ! 0.5 factor baked in, which silently assumed the OBABO placement (two
      ! calls per step).  Under BAOAB the potentiostat fires ONCE per step, so
      ! the caller must be free to say so; nvt() now passes dt_md for BAOAB and
      ! langevin() passes 0.5*dt_md for its OBABO placement.  Getting this wrong
      ! changes the damping by a factor of two without any other symptom.
      taut = (electrons%tau_potentiostat/femtoseconds)/dt

      ! ndeg = 1: q is a single genuine degree of freedom.  With ndeg /= 1 the
      ! sumnoises(ndeg-1, ...) term samples toward 2-DOF kinetic-energy
      ! statistics and rescaling_factor also disables the momentum-reversal
      ! branch that gives the correct Langevin limit at ndeg = 1.
      alpha = rescaling_factor(kk, sigma, 1, taut, electrons%gaussian_rng_stream)
      electrons%p = electrons%p*alpha
      electrons%thermostat_energy = electrons%thermostat_energy + &
                                    old_kinetic - 0.5_dp*electrons%p**2/electrons%mass

   END SUBROUTINE qs_electrons_thermostat_step

END MODULE qs_electrons_utils
