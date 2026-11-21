!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native partial-core XC source and derivatives from the same original UPF Fourier map.
MODULE qs_upf_inverse_core
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_integral_ab,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE upf_core,                        ONLY: upf_core_fourier
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_core_collocate, qs_upf_inverse_core_force

CONTAINS

! **************************************************************************************************
!> \brief Replace the auxiliary Gaussian core fit by the original core's native Fourier image.
!>        Called once for each newly prepared inverse geometry. Real and reciprocal core sources
!>        share the native real-field projection and prescribed spin split.
!> \param ks_env Live inverse KS environment with allocated native partial-core storage.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_core_collocate(ks_env)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env

      COMPLEX(dp), ALLOCATABLE                           :: values(:)
      INTEGER                                            :: atom, kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type), POINTER                      :: core_g
      TYPE(pw_r3d_rs_type), POINTER                      :: core_r
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_ks_env(ks_env, rho_nlcc=core_r, rho_nlcc_g=core_g, dft_control=control, &
                      cell=cell, particle_set=atoms, qs_kind_set=kinds)
      IF (.NOT. ASSOCIATED(core_r)) RETURN
      CPASSERT(ASSOCIATED(core_g))
      CPASSERT(ASSOCIATED(core_r%pw_grid, core_g%pw_grid))
      CPASSERT(ALL(cell%perd == 1))
      ALLOCATE (values(SIZE(core_g%array)))
      core_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(kind)%ghost) CYCLE
         CALL upf_core_fourier(kinds(kind)%sgp_potential%upf_data, core_g%pw_grid%g, atoms(atom)%r, &
                               core_g%pw_grid%vol, values)
         core_g%array = core_g%array + values/REAL(control%nspins, dp)
      END DO
      CALL pw_transfer(core_g, core_r)
      CALL pw_transfer(core_r, core_g)
   END SUBROUTINE qs_upf_inverse_core_collocate

! **************************************************************************************************
!> \brief Apply the adjoint of the original partial-core Fourier source to native XC potentials.
!>        The native XC consumer retains its own explicit functional volume/gradient terms.
!>        These derivatives include the core source's phase, reciprocal vector and 1/volume.
!> \param qs_env Native geometry and force/virial storage.
!> \param potential Folded unweighted XC spin potentials on the auxiliary PW grid.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_core_force(qs_env, potential)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)

      COMPLEX(dp), ALLOCATABLE                           :: density(:), dg(:, :), dr(:, :)
      INTEGER                                            :: atom, axis, index, kind, other, spin
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      REAL(dp)                                           :: gradient, spin_factor
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type)                               :: derivative_g, potential_g
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(virial_type), POINTER                         :: virial

      CALL get_qs_env(qs_env, cell=cell, particle_set=atoms, qs_kind_set=kinds, force=force, &
                      virial=virial, atomic_kind_set=atomic_kinds, para_env=para_env)
      CPASSERT(SIZE(potential) > 0 .AND. ASSOCIATED(force))
      CPASSERT(ALL(cell%perd == 1))
      CALL get_atomic_kind_set(atomic_kinds, atom_of_kind=atom_of_kind, kind_of=kind_of)
      CALL potential_g%create(potential(1)%pw_grid)
      CALL derivative_g%create(potential(1)%pw_grid)
      ALLOCATE (density(SIZE(derivative_g%array)), dr(3, SIZE(derivative_g%array)), &
                dg(3, SIZE(derivative_g%array)))
      spin_factor = 1.0_dp/REAL(SIZE(potential), dp)
      DO atom = 1, SIZE(atoms)
         kind = kind_of(atom)
         index = atom_of_kind(atom)
         IF (kinds(kind)%ghost) CYCLE
         IF (.NOT. kinds(kind)%sgp_potential%upf_data%core_correction) CYCLE
         CALL upf_core_fourier(kinds(kind)%sgp_potential%upf_data, derivative_g%pw_grid%g, atoms(atom)%r, &
                               derivative_g%pw_grid%vol, density, dr, dg)
         DO spin = 1, SIZE(potential)
            CPASSERT(ASSOCIATED(potential(spin)%pw_grid, potential_g%pw_grid))
            CALL pw_transfer(potential(spin), potential_g)
            DO axis = 1, 3
               derivative_g%array = spin_factor*dr(axis, :)
               gradient = pw_integral_ab(potential_g, derivative_g)
               IF (para_env%is_source()) force(kind)%gth_nlcc(axis, index) = &
                  force(kind)%gth_nlcc(axis, index) + gradient
            END DO
            IF (.NOT. virial%pv_availability .OR. virial%pv_numer) CYCLE
            DO axis = 1, 3
               DO other = 1, 3
                  derivative_g%array = atoms(atom)%r(other)*dr(axis, :) - &
                                       derivative_g%pw_grid%g(axis, :)*dg(other, :)
                  IF (axis == other) derivative_g%array = derivative_g%array - density
                  derivative_g%array = spin_factor*derivative_g%array
                  gradient = pw_integral_ab(potential_g, derivative_g)
                  IF (para_env%is_source()) THEN
                     virial%pv_nlcc(axis, other) = virial%pv_nlcc(axis, other) - gradient
                     virial%pv_virial(axis, other) = virial%pv_virial(axis, other) - gradient
                  END IF
               END DO
            END DO
         END DO
      END DO
      CALL potential_g%release()
      CALL derivative_g%release()
   END SUBROUTINE qs_upf_inverse_core_force

END MODULE qs_upf_inverse_core
