!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief all routins needed for a nonperiodic  electric field
! **************************************************************************************************

MODULE efield_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              efield_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE input_constants,                 ONLY: constant_env,&
                                              custom_env,&
                                              gaussian_env,&
                                              ramp_env
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'efield_utils'

! *** Public subroutines ***

   PUBLIC :: efield_potential_lengh_gauge, &
             calculate_ecore_efield, &
             make_field

CONTAINS

! **************************************************************************************************
!> \brief Replace the original implementation of the electric-electronic
!>        interaction in the lenght gauge. This calculation is no longer done in
!>        the grid but using matrices to match the velocity gauge implementation.
!>        Note: The energy is store in energy%core and computed later on.
!> \param qs_env ...
!> \author Guillaume Le Breton (02.23)
! **************************************************************************************************

   SUBROUTINE efield_potential_lengh_gauge(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'efield_potential_lengh_gauge'

      INTEGER                                            :: handle, i, image
      REAL(kind=dp)                                      :: field(3)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_h_kp=matrix_h, &
                      matrix_s=matrix_s)

      NULLIFY (moments)
      CALL dbcsr_allocate_matrix_set(moments, 3)
      DO i = 1, 3
         ALLOCATE (moments(i)%matrix)
         CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
         CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
      END DO

      CALL build_local_moment_matrix(qs_env, moments, 1)

      CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)

      DO i = 1, 3
         DO image = 1, dft_control%nimages
            CALL dbcsr_add(matrix_h(1, image)%matrix, moments(i)%matrix, 1.0_dp, field(i))
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(moments)

      CALL timestop(handle)

   END SUBROUTINE efield_potential_lengh_gauge

! **************************************************************************************************
!> \brief computes the amplitude of the efield within a given envelop
!> \param dft_control ...
!> \param field ...
!> \param sim_step ...
!> \param sim_time ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE make_field(dft_control, field, sim_step, sim_time)
      TYPE(dft_control_type), INTENT(IN)                 :: dft_control
      REAL(dp), INTENT(OUT)                              :: field(3)
      INTEGER, INTENT(IN)                                :: sim_step
      REAL(KIND=dp), INTENT(IN)                          :: sim_time

      INTEGER                                            :: i, lower, nfield, upper
      REAL(dp)                                           :: c, env, nu, pol(3), strength
      REAL(KIND=dp)                                      :: dt
      TYPE(efield_type), POINTER                         :: efield

      c = 137.03599962875_dp
      field = 0._dp
      nfield = SIZE(dft_control%efield_fields)
      DO i = 1, nfield
         efield => dft_control%efield_fields(i)%efield
         IF (.NOT. efield%envelop_id == custom_env) nu = c/(efield%wavelength) !in case of a custom efield we do not need nu
         strength = SQRT(efield%strength/(3.50944_dp*10.0_dp**16))
         IF (DOT_PRODUCT(efield%polarisation, efield%polarisation) == 0) THEN
            pol(:) = 1.0_dp/3.0_dp
         ELSE
            pol(:) = efield%polarisation(:)/(SQRT(DOT_PRODUCT(efield%polarisation, efield%polarisation)))
         END IF
         IF (efield%envelop_id == constant_env) THEN
            IF (sim_step .GE. efield%envelop_i_vars(1) .AND. &
                (sim_step .LE. efield%envelop_i_vars(2) .OR. efield%envelop_i_vars(2) .LT. 0)) THEN
               field = field + strength*COS(sim_time*nu*2.0_dp*pi + &
                                            efield%phase_offset*pi)*pol(:)
            END IF
         ELSE IF (efield%envelop_id == ramp_env) THEN
            IF (sim_step .GE. efield%envelop_i_vars(1) .AND. sim_step .LE. efield%envelop_i_vars(2)) &
               strength = strength*(sim_step - efield%envelop_i_vars(1))/(efield%envelop_i_vars(2) - efield%envelop_i_vars(1))
            IF (sim_step .GE. efield%envelop_i_vars(3) .AND. sim_step .LE. efield%envelop_i_vars(4)) &
               strength = strength*(efield%envelop_i_vars(4) - sim_step)/(efield%envelop_i_vars(4) - efield%envelop_i_vars(3))
            IF (sim_step .GT. efield%envelop_i_vars(4) .AND. efield%envelop_i_vars(4) .GT. 0) strength = 0.0_dp
            IF (sim_step .LE. efield%envelop_i_vars(1)) strength = 0.0_dp
            field = field + strength*COS(sim_time*nu*2.0_dp*pi + &
                                         efield%phase_offset*pi)*pol(:)
         ELSE IF (efield%envelop_id == gaussian_env) THEN
            env = EXP(-0.5_dp*((sim_time - efield%envelop_r_vars(1))/efield%envelop_r_vars(2))**2.0_dp)
            field = field + strength*env*COS(sim_time*nu*2.0_dp*pi + &
                                             efield%phase_offset*pi)*pol(:)
         ELSE IF (efield%envelop_id == custom_env) THEN
            dt = efield%envelop_r_vars(1)
            IF (sim_time .LT. (SIZE(efield%envelop_r_vars) - 2)*dt) THEN
               !make a linear interpolation between the two next points
               lower = FLOOR(sim_time/dt)
               upper = lower + 1
     strength = (efield%envelop_r_vars(lower + 2)*(upper*dt - sim_time) + efield%envelop_r_vars(upper + 2)*(sim_time - lower*dt))/dt
            ELSE
               strength = 0.0_dp
            END IF
            field = field + strength*pol(:)
         END IF
      END DO

   END SUBROUTINE make_field

! **************************************************************************************************
!> \brief Computes the force and the energy due to a efield on the cores
!>        Note: In the velocity gauge, the energy term is not added because
!>        it would lead to an unbalanced energy (center of negative charge not
!>        involved in the electric energy in this gauge).
!> \param qs_env ...
!> \param calculate_forces ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE calculate_ecore_efield(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, OPTIONAL                                  :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ecore_efield'

      INTEGER                                            :: atom_a, handle, iatom, ikind, natom, &
                                                            nkind
      INTEGER, DIMENSION(:), POINTER                     :: list
      LOGICAL                                            :: my_force
      REAL(KIND=dp)                                      :: efield_ener, zeff
      REAL(KIND=dp), DIMENSION(3)                        :: field, r
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (dft_control)
      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%apply_efield_field .OR. dft_control%apply_vector_potential) THEN
         my_force = .FALSE.
         IF (PRESENT(calculate_forces)) my_force = calculate_forces

         CALL get_qs_env(qs_env=qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         energy=energy, &
                         particle_set=particle_set, &
                         cell=cell)
         efield_ener = 0.0_dp
         nkind = SIZE(atomic_kind_set)
         CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)

         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list, natom=natom)
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)

            natom = SIZE(list)
            DO iatom = 1, natom
               IF (dft_control%apply_efield_field) THEN
                  atom_a = list(iatom)
                  r(:) = pbc(particle_set(atom_a)%r(:), cell)
                  efield_ener = efield_ener - zeff*DOT_PRODUCT(r, field)
               END IF
               IF (my_force) THEN
                  CALL get_qs_env(qs_env=qs_env, force=force)
                  force(ikind)%efield(:, iatom) = force(ikind)%efield(:, iatom) - field*zeff
               END IF
            END DO

         END DO
         IF (dft_control%apply_efield_field) energy%efield_core = efield_ener
      END IF
      CALL timestop(handle)
   END SUBROUTINE calculate_ecore_efield
END MODULE efield_utils
