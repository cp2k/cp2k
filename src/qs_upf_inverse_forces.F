!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native forces and virial from the stationary inverse USPP density and metric responses.
!>        Ionic/core forces and explicit Poisson/XC virial terms retain their native consumers.
! **************************************************************************************************
MODULE qs_upf_inverse_forces
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_integral_ab
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_upf_inverse_scf,              ONLY: qs_upf_inverse_scf_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_scf_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_force

CONTAINS

! **************************************************************************************************
!> \brief Add the stationary inverse electronic gradients to the native force and virial arrays.
!>        Replicated contractions contribute on one global rank, before the native final reduction.
!>        Remove the folded-potential operator measure before adding native functional volume terms.
!> \param qs_env Native converged environment with prepared P and energy-weighted W.
!> \param potential Folded native spin potentials, including their native dvol factor.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_force(qs_env, potential)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential(:)

      INTEGER                                            :: atom, axis, index, kind, nk, nspin, &
                                                            other, spin
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      REAL(dp)                                           :: field_measure, metric_response, &
                                                            response, strain(3, 3)
      REAL(dp), ALLOCATABLE                              :: atom_response(:, :), positions(:, :), &
                                                            reduced_k_response(:, :)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kinds(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_r3d_rs_type), ALLOCATABLE                  :: vdensity(:), vgradient(:, :), vtau(:)
      TYPE(pw_r3d_rs_type), POINTER                      :: density(:)
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(virial_type), POINTER                         :: virial

      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, force=force, cell=cell, virial=virial, &
                      particle_set=atoms, atomic_kind_set=atomic_kinds, para_env=para_env)
      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state) .AND. ASSOCIATED(force))
      CALL get_atomic_kind_set(atomic_kinds, atom_of_kind=atom_of_kind, kind_of=kind_of)
      nspin = SIZE(potential)
      nk = 1
      ALLOCATE (atom_response(3, SIZE(atoms)), positions(3, SIZE(atoms)), reduced_k_response(3, nk), &
                vdensity(nspin), vgradient(3, nspin), vtau(nspin))
      reduced_k_response = 0.0_dp
      DO atom = 1, SIZE(atoms)
         positions(:, atom) = atoms(atom)%r
      END DO
      DO spin = 1, nspin
         CALL vdensity(spin)%create(potential(spin)%pw_grid)
         vdensity(spin)%array = potential(spin)%array/potential(spin)%pw_grid%dvol
         CALL vtau(spin)%create(potential(spin)%pw_grid)
         vtau(spin)%array = 0.0_dp
         DO axis = 1, 3
            CALL vgradient(axis, spin)%create(potential(spin)%pw_grid)
            vgradient(axis, spin)%array = 0.0_dp
         END DO
      END DO
      DO atom = 1, SIZE(atoms)
         kind = kind_of(atom)
         index = atom_of_kind(atom)
         DO axis = 1, 3
            atom_response = 0.0_dp
            atom_response(axis, atom) = 1.0_dp
            CALL qs_upf_inverse_scf_response(ks_env, vdensity, vgradient, vtau, atom_response, &
                                             0.0_dp*cell%hmat, reduced_k_response, response, metric_response)
            IF (para_env%is_source()) THEN
               force(kind)%all_potential(axis, index) = force(kind)%all_potential(axis, index) + response - metric_response
               force(kind)%overlap(axis, index) = force(kind)%overlap(axis, index) + metric_response
            END IF
         END DO
      END DO
      IF (virial%pv_availability .AND. .NOT. virial%pv_numer) THEN
         IF (ANY(cell%perd == 0)) THEN
            CPABORT("Inverse analytic stress requires periodic electrostatics in all directions")
         END IF
         CALL qs_rho_get(rho, rho_r=density)
         CPASSERT(SIZE(density) == nspin)
         field_measure = 0.0_dp
         DO spin = 1, nspin
            field_measure = field_measure + pw_integral_ab(vdensity(spin), density(spin))
         END DO
         DO axis = 1, 3
            DO other = 1, 3
               strain = 0.0_dp
               strain(axis, other) = 1.0_dp
               atom_response(:, :) = MATMUL(strain, positions)
               CALL qs_upf_inverse_scf_response(ks_env, vdensity, vgradient, vtau, atom_response, &
                                                MATMUL(strain, cell%hmat), reduced_k_response, response, metric_response)
               IF (axis == other) response = response - field_measure
               IF (para_env%is_source()) THEN
                  virial%pv_virial(axis, other) = virial%pv_virial(axis, other) - response
                  virial%pv_overlap(axis, other) = virial%pv_overlap(axis, other) - metric_response
               END IF
            END DO
         END DO
      END IF
      DO spin = 1, nspin
         CALL vdensity(spin)%release()
         CALL vtau(spin)%release()
         DO axis = 1, 3
            CALL vgradient(axis, spin)%release()
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_force

END MODULE qs_upf_inverse_forces
