!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods to apply a simple Lagevin thermostat to PI runs.
!>         v_new = c1*vold + SQRT(kT/m)*c2*random
!> \author Felix Uhl
!> \par History
!>      10.2014 created [Felix Uhl]
! **************************************************************************************************
MODULE pint_pile
   USE input_constants,                 ONLY: propagator_cmd
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                              rng_record_length,&
                                              rng_stream_type,&
                                              rng_stream_type_from_record
   USE pint_types,                      ONLY: normalmode_env_type,&
                                              pile_therm_type,&
                                              pint_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pint_pile_step, &
             pint_pile_init, &
             pint_pile_release, &
             pint_calc_pile_energy

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_pile'

CONTAINS

   ! ***************************************************************************
   !> \brief initializes the data for a pile run
   !> \author Felix Uhl
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pile_therm ...
!> \param pint_env ...
!> \param normalmode_env ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE pint_pile_init(pile_therm, pint_env, normalmode_env, section)
      TYPE(pile_therm_type), POINTER                     :: pile_therm
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(normalmode_env_type), POINTER                 :: normalmode_env
      TYPE(section_vals_type), POINTER                   :: section

      CHARACTER(LEN=rng_record_length)                   :: rng_record
      INTEGER                                            :: i, i_propagator, j, p
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: dti2, ex
      REAL(KIND=dp), DIMENSION(3, 2)                     :: initial_seed
      TYPE(section_vals_type), POINTER                   :: pint_section, rng_section

      pint_env%e_pile = 0.0_dp
      ALLOCATE (pile_therm)
      pile_therm%ref_count = 1
      pile_therm%thermostat_energy = 0.0_dp
      !Get input parameter
      CALL section_vals_val_get(section, "TAU", r_val=pile_therm%tau)
      CALL section_vals_val_get(section, "LAMBDA", r_val=pile_therm%lamb)
      CALL section_vals_val_get(section, "THERMOSTAT_ENERGY", r_val=pile_therm%thermostat_energy)
      pint_section => section_vals_get_subs_vals(pint_env%input, "MOTION%PINT")
      CALL section_vals_val_get(pint_section, "PROPAGATOR", i_val=i_propagator)

      IF (i_propagator == propagator_cmd) THEN
         pile_therm%tau = 0.0_dp
      END IF

      p = pint_env%p
      dti2 = 0.5_dp*pint_env%dt
      ALLOCATE (pile_therm%c1(p))
      ALLOCATE (pile_therm%c2(p))
      ALLOCATE (pile_therm%g_fric(p))
      ALLOCATE (pile_therm%massfact(p, pint_env%ndim))
      !Initialize everything
      ! If tau is negative or zero the thermostat does not act on the centroid
      ! (TRPMD)
      IF (pile_therm%tau <= 0.0_dp) THEN
         pile_therm%g_fric(1) = 0.0_dp
      ELSE
         pile_therm%g_fric(1) = 1.0_dp/pile_therm%tau
      END IF
      DO i = 2, p
         pile_therm%g_fric(i) = 2.0_dp*pile_therm%lamb*SQRT(normalmode_env%lambda(i))
      END DO
      DO i = 1, p
         ex = -dti2*pile_therm%g_fric(i)
         pile_therm%c1(i) = EXP(ex)
         ex = pile_therm%c1(i)*pile_therm%c1(i)
         pile_therm%c2(i) = SQRT(1.0_dp - ex)
      END DO
      DO j = 1, pint_env%ndim
         DO i = 1, pint_env%p
            pile_therm%massfact(i, j) = SQRT(pint_env%kT/pint_env%mass_fict(i, j))
         END DO
      END DO

      !prepare Random number generator
      NULLIFY (rng_section)
      rng_section => section_vals_get_subs_vals(section, &
                                                subsection_name="RNG_INIT")
      CALL section_vals_get(rng_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(rng_section, "_DEFAULT_KEYWORD_", &
                                   i_rep_val=1, c_val=rng_record)

         pile_therm%gaussian_rng_stream = rng_stream_type_from_record(rng_record)
      ELSE
         initial_seed(:, :) = REAL(pint_env%thermostat_rng_seed, dp)
         pile_therm%gaussian_rng_stream = rng_stream_type( &
                                          name="pile_rng_gaussian", distribution_type=GAUSSIAN, &
                                          extended_precision=.TRUE., &
                                          seed=initial_seed)
      END IF

   END SUBROUTINE pint_pile_init

! **************************************************************************************************
!> \brief ...
!> \param vold ...
!> \param vnew ...
!> \param p ...
!> \param ndim ...
!> \param first_mode ...
!> \param masses ...
!> \param pile_therm ...
! **************************************************************************************************
   SUBROUTINE pint_pile_step(vold, vnew, p, ndim, first_mode, masses, pile_therm)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vold, vnew
      INTEGER, INTENT(IN)                                :: p, ndim, first_mode
      REAL(kind=dp), DIMENSION(:, :), INTENT(IN)         :: masses
      TYPE(pile_therm_type), POINTER                     :: pile_therm

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pint_pile_step'

      INTEGER                                            :: handle, ibead, idim
      REAL(KIND=dp)                                      :: delta_ekin

      CALL timeset(routineN, handle)
      delta_ekin = 0.0_dp
      DO idim = 1, ndim
         DO ibead = first_mode, p
            vnew(ibead, idim) = pile_therm%c1(ibead)*vold(ibead, idim) + &
                                pile_therm%massfact(ibead, idim)*pile_therm%c2(ibead)* &
                                pile_therm%gaussian_rng_stream%next()
            delta_ekin = delta_ekin + masses(ibead, idim)*( &
                         vnew(ibead, idim)*vnew(ibead, idim) - &
                         vold(ibead, idim)*vold(ibead, idim))
         END DO
      END DO
      pile_therm%thermostat_energy = pile_therm%thermostat_energy - 0.5_dp*delta_ekin

      CALL timestop(handle)
   END SUBROUTINE pint_pile_step

   ! ***************************************************************************
   !> \brief releases the pile environment
   !> \param pile_therm pile data to be released
   !> \author Felix Uhl
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pile_therm ...
! **************************************************************************************************
   SUBROUTINE pint_pile_release(pile_therm)

      TYPE(pile_therm_type), POINTER                     :: pile_therm

      IF (ASSOCIATED(pile_therm)) THEN
         pile_therm%ref_count = pile_therm%ref_count - 1
         IF (pile_therm%ref_count == 0) THEN
            DEALLOCATE (pile_therm%c1)
            DEALLOCATE (pile_therm%c2)
            DEALLOCATE (pile_therm%g_fric)
            DEALLOCATE (pile_therm%massfact)
            DEALLOCATE (pile_therm)
         END IF
      END IF
      NULLIFY (pile_therm)

      RETURN
   END SUBROUTINE pint_pile_release

   ! ***************************************************************************
   !> \brief returns the pile kinetic energy contribution
   !> \author Felix Uhl
   ! ***************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
! **************************************************************************************************
   SUBROUTINE pint_calc_pile_energy(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      IF (ASSOCIATED(pint_env%pile_therm)) THEN
         pint_env%e_pile = pint_env%pile_therm%thermostat_energy
      END IF

      RETURN

   END SUBROUTINE pint_calc_pile_energy
END MODULE
