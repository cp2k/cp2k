!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Force-based FIRE and FIRE 2.0 steps, independent of force evaluation.
! **************************************************************************************************
MODULE fire_minimizer
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: bohr,&
                                              femtoseconds,&
                                              massunit
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fire_parameters_type, fire_state_type, fire_init, fire_prepare_step, &
             fire_step, fire_accept_step

   TYPE fire_parameters_type
      REAL(KIND=dp) :: timestep = 0.5_dp/femtoseconds
      REAL(KIND=dp) :: max_timestep_factor = 10.0_dp
      REAL(KIND=dp) :: max_step = 0.2_dp*bohr
   END TYPE fire_parameters_type

   TYPE fire_state_type
      PRIVATE
      TYPE(fire_parameters_type) :: params = fire_parameters_type()
      REAL(KIND=dp), ALLOCATABLE :: velocity(:), SCALE(:)
      REAL(KIND=dp) :: dt = 0.0_dp, alpha = 0.0_dp, alpha_start = 0.0_dp, mixing = 0.0_dp
      INTEGER :: ncart = 0, iterations = 0, downhill = 0, uphill = 0, delay = 0
      LOGICAL :: fire2 = .FALSE.
   END TYPE fire_state_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize fictitious dynamics in mass-scaled coordinates.
!> \param state integrator state
!> \param params timestep and displacement bounds
!> \param ndof number of coordinates
!> \param ncart number of Cartesian components at the start of the vector
!> \param fire2 use FIRE 2.0 instead of FIRE
!> \param tail_mass mass in amu for non-Cartesian components, relative to 1 amu for atoms
! **************************************************************************************************
   SUBROUTINE fire_init(state, params, ndof, ncart, fire2, tail_mass)
      TYPE(fire_state_type), INTENT(OUT)                 :: state
      TYPE(fire_parameters_type), INTENT(IN)             :: params
      INTEGER, INTENT(IN)                                :: ndof, ncart
      LOGICAL, INTENT(IN)                                :: fire2
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tail_mass

      CPASSERT(ndof > 0)
      CPASSERT(ncart >= 0 .AND. ncart <= ndof .AND. MOD(ncart, 3) == 0)
      CPASSERT(params%timestep > 0.0_dp)
      CPASSERT(params%max_timestep_factor >= 1.0_dp)
      CPASSERT(params%max_step > 0.0_dp)
      state%params = params
      state%fire2 = fire2
      state%ncart = ncart
      state%iterations = 0
      state%downhill = 0
      state%uphill = 0
      state%alpha_start = MERGE(0.25_dp, 0.1_dp, fire2)
      state%delay = MERGE(20, 5, fire2)
      state%dt = params%timestep
      state%alpha = state%alpha_start
      state%mixing = 0.0_dp
      ALLOCATE (state%velocity(ndof), state%scale(ndof))
      state%velocity = 0.0_dp
      state%scale = 1.0_dp/SQRT(massunit)
      IF (PRESENT(tail_mass)) THEN
         CPASSERT(tail_mass > 0.0_dp)
         state%scale(ncart + 1:) = state%scale(ncart + 1:)/SQRT(tail_mass)
      END IF
   END SUBROUTINE fire_init

! **************************************************************************************************
!> \brief Adapt the timestep and return the FIRE 2.0 half-step correction, if any.
!> \param state integrator state
!> \param force force at the current position
!> \param backstep displacement to apply before taking the next step
!> \param stalled whether the FIRE 2.0 limit on consecutive resets was reached
! **************************************************************************************************
   SUBROUTINE fire_prepare_step(state, force, backstep, stalled)
      TYPE(fire_state_type), INTENT(INOUT)               :: state
      REAL(KIND=dp), INTENT(IN)                          :: force(:)
      REAL(KIND=dp), INTENT(OUT)                         :: backstep(:)
      LOGICAL, INTENT(OUT)                               :: stalled

      REAL(KIND=dp)                                      :: reduced_dt

      stalled = .FALSE.
      backstep = 0.0_dp
      state%mixing = 0.0_dp
      IF (state%iterations == 0) RETURN
      IF (ALL(force == 0.0_dp)) THEN
         state%velocity = 0.0_dp
         RETURN
      END IF
      IF (DOT_PRODUCT(force*state%scale, state%velocity) > 0.0_dp) THEN
         state%uphill = 0
         state%downhill = state%downhill + 1
         state%mixing = state%alpha
         IF (state%downhill > state%delay) THEN
            state%dt = MIN(1.1_dp*state%dt, state%params%max_timestep_factor*state%params%timestep)
            state%alpha = 0.99_dp*state%alpha
         END IF
         IF (state%fire2) state%mixing = state%alpha
      ELSE
         state%downhill = 0
         state%uphill = state%uphill + 1
         IF (state%fire2 .AND. state%uphill > 2000) THEN
            stalled = .TRUE.
            RETURN
         END IF
         IF (.NOT. state%fire2 .OR. state%iterations >= state%delay) THEN
            reduced_dt = 0.5_dp*state%dt
            IF (.NOT. state%fire2 .OR. reduced_dt >= 0.02_dp*state%params%timestep) state%dt = reduced_dt
            state%alpha = state%alpha_start
         END IF
         IF (state%fire2) backstep = -0.5_dp*state%dt*state%scale*state%velocity
         state%velocity = 0.0_dp
      END IF
   END SUBROUTINE fire_prepare_step

! **************************************************************************************************
!> \brief Advance with semi-implicit Euler using forces at the (possibly corrected) position.
!> \param state integrator state
!> \param force force at the position from which to advance
!> \param step proposed displacement, capped per atom and for the cell block
! **************************************************************************************************
   SUBROUTINE fire_step(state, force, step)
      TYPE(fire_state_type), INTENT(INOUT)               :: state
      REAL(KIND=dp), INTENT(IN)                          :: force(:)
      REAL(KIND=dp), INTENT(OUT)                         :: step(:)

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: force_norm, max_move
      REAL(KIND=dp), DIMENSION(SIZE(force))              :: scaled_force

      scaled_force = state%scale*force
      force_norm = NORM2(scaled_force)
      IF (state%fire2) state%velocity(:) = state%velocity + state%dt*scaled_force
      IF (state%mixing > 0.0_dp .AND. force_norm > TINY(1.0_dp)) THEN
         state%velocity(:) = (1.0_dp - state%mixing)*state%velocity + &
                             state%mixing*NORM2(state%velocity)*(scaled_force/force_norm)
      END IF
      IF (.NOT. state%fire2) state%velocity(:) = state%velocity + state%dt*scaled_force

      step = state%dt*state%scale*state%velocity
      max_move = NORM2(step(state%ncart + 1:))
      DO i = 1, state%ncart, 3
         max_move = MAX(max_move, NORM2(step(i:i + 2)))
      END DO
      IF (max_move > state%params%max_step) THEN
         step = step*(state%params%max_step/max_move)
         CALL fire_accept_step(state, step)
      END IF
      state%iterations = state%iterations + 1
   END SUBROUTINE fire_step

! **************************************************************************************************
!> \brief Synchronize the velocity after the caller projects or shortens a step.
!> \param state integrator state
!> \param step actual displacement from the corrected position
! **************************************************************************************************
   SUBROUTINE fire_accept_step(state, step)
      TYPE(fire_state_type), INTENT(INOUT)               :: state
      REAL(KIND=dp), INTENT(IN)                          :: step(:)

      state%velocity(:) = (step/state%scale)/state%dt
   END SUBROUTINE fire_accept_step

END MODULE fire_minimizer
