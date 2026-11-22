!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief PAW one-center functional in the native inverse Gaussian representation.
!>        Fixed kinetic and ionic terms remain in the original PP_DIJ.
! **************************************************************************************************
MODULE qs_upf_inverse_paw
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
   USE paw_onecenter,                   ONLY: paw_onecenter_hartree,&
                                              paw_onecenter_init
   USE paw_types,                       ONLY: paw_kind_type,&
                                              paw_potential_type
   USE paw_xc,                          ONLY: paw_xc_evaluate,&
                                              paw_xc_init
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_response
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_grid_type
   USE upf_projector_integrals,         ONLY: upf_projector_type
   USE upf_types,                       ONLY: atom_upfpot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_paw_init, qs_upf_inverse_paw_evaluate, qs_upf_inverse_paw_response

CONTAINS

! **************************************************************************************************
!> \brief Prepare original PAW radial kernels without changing the supplied operators or waves.
!> \param pot Original scalar PAW dataset.
!> \param projectors Matching radial and magnetic projector mapping.
!> \param nang Requested atomic angular quadrature.
!> \param cache Owned per-kind PAW data.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_init(pot, projectors, nang, cache)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      INTEGER, INTENT(IN)                                :: nang
      TYPE(paw_kind_type), INTENT(OUT)                   :: cache

      CPASSERT(pot%is_paw .AND. .NOT. pot%has_so)
      CPASSERT(ALLOCATED(pot%paw))
      CALL paw_onecenter_init(pot, projectors, cache%onecenter)
      CALL paw_xc_init(pot, nang, cache%xcgrid)
      cache%core_energy = pot%paw%core_energy
      cache%active = .TRUE.
   END SUBROUTINE qs_upf_inverse_paw_init

! **************************************************************************************************
!> \brief Evaluate PAW Hartree, XC and core energy and their Gaussian Hamiltonian derivative.
!> \param kinds Prepared per-kind PAW kernels. Inactive NC and USPP kinds contribute zero.
!> \param frames Native inverse frames for the complete global k mesh, identically ordered on all ranks.
!> \param density Hermitian Gaussian occupation matrices, AO by AO by global k point by spin.
!> \param kweights Normalized global k weights. They enter the projected density once.
!> \param functional Activated XC_FUNCTIONAL section.
!> \param density_cutoff XC density cutoff.
!> \param gradient_cutoff XC gradient cutoff.
!> \param energies One-center Hartree, XC and frozen-core constant, counted once per atom.
!> \param hamiltonian Unweighted Gaussian Hamiltonian correction for each global k point and spin.
!> \param potential Owned per-center projector potentials for the corresponding geometry response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_evaluate(kinds, frames, density, kweights, functional, density_cutoff, &
                                          gradient_cutoff, energies, hamiltonian, potential)
      TYPE(paw_kind_type), INTENT(IN)                    :: kinds(:)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(section_vals_type), POINTER                   :: functional
      REAL(dp), INTENT(IN)                               :: density_cutoff, gradient_cutoff
      REAL(dp), INTENT(OUT)                              :: energies(3)
      COMPLEX(dp), INTENT(OUT)                           :: hamiltonian(:, :, :, :)
      TYPE(paw_potential_type), ALLOCATABLE, INTENT(OUT) :: potential(:)

      COMPLEX(dp), ALLOCATABLE                           :: projected(:, :)
      INTEGER                                            :: center, first, ik, ikind, last, n, nao, &
                                                            nk, nproj, ns, spin
      REAL(dp)                                           :: eh, ex, hermitian_error, hermitian_limit
      REAL(dp), ALLOCATABLE                              :: dm(:, :, :), vh(:, :), vxc(:, :, :)

      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      nao = SIZE(frames(1)%ao%coefficients, 2)
      nproj = SIZE(frames(1)%ao%coefficients, 1)
      ns = SIZE(density, 4)
      CPASSERT(ns == 1 .OR. ns == 2)
      CPASSERT(ALL(SHAPE(density) == [nao, nao, nk, ns]))
      CPASSERT(ALL(SHAPE(hamiltonian) == SHAPE(density)))
      CPASSERT(ALL(ieee_is_finite(kweights)) .AND. ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) < 1.e-12_dp)
      DO ik = 1, nk
         CPASSERT(SIZE(frames(ik)%ao%coefficients, 1) == nproj)
         CPASSERT(SIZE(frames(ik)%ao%coefficients, 2) == nao)
         CPASSERT(ALL(frames(ik)%ao%center_atoms == frames(1)%ao%center_atoms))
         CPASSERT(ALL(frames(ik)%ao%center_kinds == frames(1)%ao%center_kinds))
         IF (.NOT. (ALL(frames(ik)%ao%projector_first == frames(1)%ao%projector_first))) THEN
            CPABORT("Inconsistent inverse PAW frame or density")
         END IF
         CPASSERT(ALL(frames(ik)%ao%projector_last == frames(1)%ao%projector_last))
         DO spin = 1, ns
            CPASSERT(ALL(ieee_is_finite(REAL(density(:, :, ik, spin), dp))))
            CPASSERT(ALL(ieee_is_finite(AIMAG(density(:, :, ik, spin)))))
            hermitian_error = MAXVAL(ABS(density(:, :, ik, spin) - CONJG(TRANSPOSE(density(:, :, ik, spin)))))
            hermitian_limit = 1.e-12_dp*(1.0_dp + MAXVAL(ABS(density(:, :, ik, spin))))
            CPASSERT(hermitian_error <= hermitian_limit)
         END DO
      END DO
      ALLOCATE (potential(SIZE(frames(1)%ao%center_atoms)))
      energies = 0.0_dp
      hamiltonian = 0.0_dp
      DO center = 1, SIZE(potential)
         ikind = frames(1)%ao%center_kinds(center)
         CPASSERT(ikind > 0 .AND. ikind <= SIZE(kinds))
         IF (.NOT. kinds(ikind)%active) CYCLE
         first = frames(1)%ao%projector_first(center)
         last = frames(1)%ao%projector_last(center)
         n = last - first + 1
         CPASSERT(n == kinds(ikind)%onecenter%nproj)
         ALLOCATE (potential(center)%matrix(n, n, ns), SOURCE=0.0_dp)
         ! The global density and inverse coefficients are replicated. One rank owns
         ! each atomic functional, including the core constant, independently of k groups.
         IF (MODULO(center - 1, frames(1)%ao%point_comm%num_pe) == frames(1)%ao%point_comm%mepos) THEN
            ALLOCATE (dm(n, n, ns), projected(n, n), vh(n, n), vxc(n, n, ns))
            dm = 0.0_dp
            DO spin = 1, ns
               DO ik = 1, nk
                  projected(:, :) = MATMUL(frames(ik)%ao%coefficients(first:last, :), &
                                           MATMUL(density(:, :, ik, spin), &
                                                  CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :)))))
                  dm(:, :, spin) = dm(:, :, spin) + kweights(ik)*REAL(projected, dp)
               END DO
               dm(:, :, spin) = 0.5_dp*(dm(:, :, spin) + TRANSPOSE(dm(:, :, spin)))
            END DO
            CALL paw_onecenter_hartree(kinds(ikind)%onecenter, SUM(dm, DIM=3), eh, vh)
            CALL paw_xc_evaluate(kinds(ikind)%xcgrid, kinds(ikind)%onecenter, functional, dm, &
                                 density_cutoff, gradient_cutoff, ex, vxc)
            energies = energies + [eh, ex, kinds(ikind)%core_energy]
            DO spin = 1, ns
               potential(center)%matrix(:, :, spin) = vh + vxc(:, :, spin)
            END DO
            DEALLOCATE (dm, projected, vh, vxc)
         END IF
         CALL frames(1)%ao%point_comm%sum(potential(center)%matrix)
         DO spin = 1, ns
            DO ik = 1, nk
               hamiltonian(:, :, ik, spin) = hamiltonian(:, :, ik, spin) + &
                                             MATMUL(CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :))), &
                                                    MATMUL(CMPLX(potential(center)%matrix(:, :, spin), KIND=dp), &
                                                           frames(ik)%ao%coefficients(first:last, :)))
            END DO
         END DO
      END DO
      CALL frames(1)%ao%point_comm%sum(energies)
   END SUBROUTINE qs_upf_inverse_paw_evaluate

! **************************************************************************************************
!> \brief Differentiate the same PAW functional through the native inverse transformation.
!> \param frames Native frames used for the energy and projector potential.
!> \param density Fixed Gaussian occupation matrices.
!> \param kweights Normalized global k weights.
!> \param potential PAW projector derivative at this density and geometry.
!> \param atom_response Cartesian displacements of the atomic positions.
!> \param cell_response Cell response with fixed fractional grid coordinates.
!> \param reduced_k_response Response of the reduced Bloch vectors, one per global k point.
!> \param derivative Energy response. Rigid atomic quadratures have no additional measure response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_paw_response(frames, density, kweights, potential, atom_response, cell_response, &
                                          reduced_k_response, derivative)
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frames(:)
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :, :)
      REAL(dp), INTENT(IN)                               :: kweights(:)
      TYPE(paw_potential_type), INTENT(IN)               :: potential(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: derivative

      COMPLEX(dp), ALLOCATABLE                           :: dw(:, :), projected(:, :)
      INTEGER                                            :: center, first, ik, info, last, nao, nk, &
                                                            nproj, ns, spin
      REAL(dp)                                           :: inverse_cell(3, 3), strain(3, 3), &
                                                            volume_response
      REAL(dp), ALLOCATABLE                              :: points_response(:, :), &
                                                            weights_response(:)

      nk = SIZE(frames)
      CPASSERT(nk > 0 .AND. SIZE(kweights) == nk)
      nao = SIZE(frames(1)%ao%coefficients, 2)
      nproj = SIZE(frames(1)%ao%coefficients, 1)
      ns = SIZE(density, 4)
      CPASSERT(ALL(SHAPE(density) == [nao, nao, nk, ns]))
      CPASSERT(SIZE(potential) == SIZE(frames(1)%ao%center_atoms))
      CPASSERT(ALL(SHAPE(reduced_k_response) == [3, nk]))
      CPASSERT(ALL(ieee_is_finite(kweights)) .AND. ALL(kweights >= 0.0_dp))
      CPASSERT(ABS(SUM(kweights) - 1.0_dp) < 1.e-12_dp)
      ALLOCATE (dw(nproj, nao))
      derivative = 0.0_dp
      DO ik = 1, nk
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frames(ik)%ao%positions)))
         inverse_cell = frames(ik)%ao%cell
         CALL invmat(inverse_cell, info)
         CPASSERT(info == 0)
         strain = MATMUL(cell_response, inverse_cell)
         volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
         ALLOCATE (points_response(3, SIZE(frames(ik)%ao%weights)), &
                   weights_response(SIZE(frames(ik)%ao%weights)))
         points_response(:, :) = MATMUL(strain, frames(ik)%ao%points)
         weights_response(:) = volume_response*frames(ik)%ao%weights
         CALL qs_upf_inverse_response(frames(ik)%ao, points_response, atom_response, cell_response, &
                                      reduced_k_response(:, ik), weights_response, dw)
         DEALLOCATE (points_response, weights_response)
         DO center = 1, SIZE(potential)
            IF (.NOT. ALLOCATED(potential(center)%matrix)) CYCLE
            first = frames(ik)%ao%projector_first(center)
            last = frames(ik)%ao%projector_last(center)
            CPASSERT(SIZE(potential(center)%matrix, 1) == last - first + 1)
            CPASSERT(SIZE(potential(center)%matrix, 2) == last - first + 1)
            CPASSERT(SIZE(potential(center)%matrix, 3) == ns)
            ALLOCATE (projected(last - first + 1, last - first + 1))
            DO spin = 1, ns
               projected(:, :) = MATMUL(dw(first:last, :), MATMUL(density(:, :, ik, spin), &
                                                                  CONJG(TRANSPOSE(frames(ik)%ao%coefficients(first:last, :)))))
               derivative = derivative + 2.0_dp*kweights(ik)* &
                            SUM(potential(center)%matrix(:, :, spin)*REAL(TRANSPOSE(projected), dp))
            END DO
            DEALLOCATE (projected)
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_paw_response

END MODULE qs_upf_inverse_paw
