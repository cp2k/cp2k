!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Direct reciprocal-space reference for scalar UPF augmentation.
!>        Uses the complete tabulated radial mesh and CP2K real spherical harmonics.
!>        The overlap correction is the zero-frequency limit of this same transform.
! **************************************************************************************************
MODULE uspp_augmentation
   USE atom_upf,                        ONLY: atom_upfpot_type,&
                                              upf_radial_weights
   USE bessel_lib,                      ONLY: bessel0
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE spherical_harmonics,             ONLY: y_lm
   USE upf_projector_integrals,         ONLY: upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'uspp_augmentation'

   TYPE uspp_augmentation_type
      INTEGER :: nproj = 0, lmax = -1
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:), lbeta(:)
      REAL(KIND=dp), ALLOCATABLE :: r(:), weighted_q(:, :, :), gaunt(:, :, :), q(:, :)
   END TYPE uspp_augmentation_type

   PUBLIC :: uspp_augmentation_type, uspp_augmentation_init, uspp_augmentation_fourier, &
             uspp_augmentation_collocate, uspp_augmentation_integrate

CONTAINS

! **************************************************************************************************
!> \brief Prepare radial quadrature and Gaunt coefficients in the native beta convention.
!> \param pot Original scalar USPP with angular-resolved augmentation.
!> \param projectors Prepared native projectors, with Cartesian degree at least max(lbeta).
!> \param augmentation Immutable augmentation data, in the same magnetic-projector ordering.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_init(pot, projectors, augmentation)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_projector_type), INTENT(IN)               :: projectors
      TYPE(uspp_augmentation_type), INTENT(OUT)          :: augmentation

      INTEGER                                            :: i, j, l, m, n, nc, p, pair, q, s
      REAL(KIND=dp), ALLOCATABLE                         :: weights(:)

      CPASSERT(.NOT. pot%has_so .AND. .NOT. pot%is_paw)
      CPASSERT(ALLOCATED(pot%qijl))
      CPASSERT(projectors%degree >= MAXVAL(pot%lbeta))
      CPASSERT(ALL(projectors%lbeta == pot%lbeta))
      n = SIZE(pot%r)
      CPASSERT(n >= 3 .AND. SIZE(pot%rab) == n)
      CPASSERT(ALL(pot%r >= 0.0_dp) .AND. ALL(pot%rab > 0.0_dp))
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)))
      nc = pot%number_of_proj*(pot%number_of_proj + 1)/2
      augmentation%nproj = projectors%nproj
      augmentation%lmax = 2*MAXVAL(pot%lbeta)
      CALL init_orbital_pointers(augmentation%lmax)
      CALL init_spherical_harmonics(augmentation%lmax, -1)
      CPASSERT(ALL(SHAPE(pot%qijl) == [n, nc, augmentation%lmax + 1]))
      CPASSERT(LBOUND(pot%qijl, 3) == 0)
      ALLOCATE (augmentation%r(n), weights(n))
      augmentation%r(:) = pot%r
      weights(:) = upf_radial_weights(pot)
      ALLOCATE (augmentation%weighted_q(n, nc, 0:augmentation%lmax))
      DO l = 0, augmentation%lmax
         DO pair = 1, nc
            augmentation%weighted_q(:, pair, l) = weights*pot%qijl(:, pair, l)
         END DO
      END DO
      n = augmentation%nproj
      ALLOCATE (augmentation%radial_index(n), augmentation%angular_index(n), &
                augmentation%lbeta(pot%number_of_proj))
      augmentation%radial_index(:) = projectors%radial_index
      augmentation%angular_index(:) = projectors%angular_index
      augmentation%lbeta(:) = pot%lbeta
      ALLOCATE (augmentation%gaunt(n, n, (augmentation%lmax + 1)**2), augmentation%q(n, n))
      augmentation%gaunt = 0.0_dp
      augmentation%q = 0.0_dp
      DO q = 1, n
         j = projectors%radial_index(q)
         l = pot%lbeta(j)
         m = projectors%angular_index(q) - l*l
         DO p = 1, n
            i = projectors%radial_index(p)
            ! Contract the existing analytic monomial/Y/Y integrals with Y_j's polynomial.
            DO s = 1, nco(l)
               augmentation%gaunt(p, q, :) = augmentation%gaunt(p, q, :) + &
                                             orbtramat(l)%slm(m, s)*projectors%gaunt(ncoset(l - 1) + s, &
                                                                           projectors%angular_index(p), :(augmentation%lmax + 1)**2)
            END DO
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            IF (projectors%angular_index(p) == projectors%angular_index(q)) THEN
               augmentation%q(p, q) = SUM(augmentation%weighted_q(:, pair, 0))
            END IF
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_init

! **************************************************************************************************
!> \brief Evaluate Q_ij(G) = integral Q_ij(r-R) exp(-i G.r) dr / volume.
!>        Direct radial quadrature is a reference without reciprocal interpolation.
!> \param augmentation Prepared data.
!> \param g Cartesian reciprocal vector in inverse bohr.
!> \param position Augmentation center in bohr.
!> \param volume Cell volume in bohr cubed.
!> \param qg Complex symmetric form-factor matrix. Q(-G) = conjugate(Q(G)).
!> \param strain Optional dQ/d epsilon_ab at fixed fractional center and reciprocal index.
!>        Includes the inverse-volume derivative; the center phase is strain invariant.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_fourier(augmentation, g, position, volume, qg, strain)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: g(3), position(3), volume
      COMPLEX(KIND=dp), INTENT(OUT)                      :: qg(:, :)
      COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL            :: strain(:, :, :, :)

      COMPLEX(KIND=dp)                                   :: phase, phase_l
      INTEGER                                            :: a, axis, i, j, l, lm, m, p, pair, &
                                                            powers(3), q
      REAL(KIND=dp)                                      :: angular, angular_derivative(3), &
                                                            direction(3, 1), monomial, norm_g, &
                                                            value(1)
      REAL(KIND=dp), ALLOCATABLE                         :: harmonics(:), radial(:, :), &
                                                            radial_next(:, :), &
                                                            solid_gradient(:, :), work(:)

      CPASSERT(volume > 0.0_dp)
      CPASSERT(ALL(SHAPE(qg) == augmentation%nproj))
      IF (PRESENT(strain)) THEN
         CPASSERT(ALL(SHAPE(strain) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         strain = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      norm_g = SQRT(SUM(g*g))
      IF (norm_g == 0.0_dp) THEN
         qg(:, :) = CMPLX(augmentation%q/volume, 0.0_dp, KIND=dp)
         IF (PRESENT(strain)) THEN
            DO axis = 1, 3
               strain(:, :, axis, axis) = -qg
            END DO
         END IF
         RETURN
      END IF
      direction(:, 1) = g/norm_g
      ALLOCATE (work(SIZE(augmentation%r)), harmonics((augmentation%lmax + 1)**2))
      ALLOCATE (radial(SIZE(augmentation%weighted_q, 2), 0:augmentation%lmax))
      IF (PRESENT(strain)) THEN
         ALLOCATE (solid_gradient(3, (augmentation%lmax + 1)**2), &
                   radial_next(SIZE(radial, 1), 0:augmentation%lmax))
         solid_gradient = 0.0_dp
      END IF
      DO l = 0, augmentation%lmax
         IF (PRESENT(strain)) THEN
            work(:) = augmentation%r*bessel0(norm_g*augmentation%r, l + 1)
            radial_next(:, l) = MATMUL(work, augmentation%weighted_q(:, :, l))
         END IF
         work(:) = bessel0(norm_g*augmentation%r, l)
         radial(:, l) = MATMUL(work, augmentation%weighted_q(:, :, l))
         DO m = -l, l
            lm = l*l + l + m + 1
            CALL y_lm(direction, value, l, m)
            harmonics(lm) = value(1)
            IF (PRESENT(strain)) THEN
               ! Cartesian solid-harmonic derivatives avoid polar-axis singularities.
               DO a = 1, nco(l)
                  DO axis = 1, 3
                     powers = indco(:, ncoset(l - 1) + a)
                     IF (powers(axis) == 0) CYCLE
                     monomial = REAL(powers(axis), dp)
                     powers(axis) = powers(axis) - 1
                     solid_gradient(axis, lm) = solid_gradient(axis, lm) + &
                                                orbtramat(l)%slm(l + m + 1, a)*monomial*PRODUCT(direction(:, 1)**powers)
                  END DO
               END DO
            END IF
         END DO
      END DO
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, position), KIND=dp))*fourpi/volume
      qg = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO q = 1, augmentation%nproj
         j = augmentation%radial_index(q)
         DO p = 1, q
            i = augmentation%radial_index(p)
            pair = MAX(i, j)*(MAX(i, j) - 1)/2 + MIN(i, j)
            DO l = ABS(augmentation%lbeta(i) - augmentation%lbeta(j)), &
               augmentation%lbeta(i) + augmentation%lbeta(j), 2
               angular = DOT_PRODUCT(augmentation%gaunt(p, q, l*l + 1:(l + 1)**2), &
                                     harmonics(l*l + 1:(l + 1)**2))
               phase_l = CMPLX(0.0_dp, -1.0_dp, KIND=dp)**l
               qg(p, q) = qg(p, q) + phase_l*angular*radial(pair, l)
               IF (PRESENT(strain)) THEN
                  angular_derivative = MATMUL(solid_gradient(:, l*l + 1:(l + 1)**2), &
                                              augmentation%gaunt(p, q, l*l + 1:(l + 1)**2))
                  ! d[j_l(Gr)Y_l(G)]/dG combines the l/G terms analytically.
                  angular_derivative = radial(pair, l)/norm_g*angular_derivative &
                                       - radial_next(pair, l)*angular*direction(:, 1)
                  DO axis = 1, 3
                     strain(p, q, :, axis) = strain(p, q, :, axis) - phase*phase_l*g*angular_derivative(axis)
                  END DO
               END IF
            END DO
            qg(p, q) = phase*qg(p, q)
            qg(q, p) = qg(p, q)
            IF (PRESENT(strain)) THEN
               DO axis = 1, 3
                  strain(p, q, axis, axis) = strain(p, q, axis, axis) - qg(p, q)
               END DO
               strain(q, p, :, :) = strain(p, q, :, :)
            END IF
         END DO
      END DO
   END SUBROUTINE uspp_augmentation_fourier

! **************************************************************************************************
!> \brief Add the augmentation density on a list of reciprocal vectors.
!> \param augmentation Prepared data.
!> \param gvecs Cartesian reciprocal vectors (3,ng).
!> \param position Center in bohr.
!> \param volume Cell volume.
!> \param density_matrix Real symmetric projector density matrix, with occupations included.
!> \param rho Fourier density coefficients (integral divided by volume), accumulated in place.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_collocate(augmentation, gvecs, position, volume, density_matrix, rho)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: gvecs(:, :), position(3), volume, &
                                                            density_matrix(:, :)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: rho(:)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: qg(:, :)
      INTEGER                                            :: ig

      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(rho))
      CPASSERT(ALL(SHAPE(density_matrix) == augmentation%nproj))
      ALLOCATE (qg(augmentation%nproj, augmentation%nproj))
      DO ig = 1, SIZE(rho)
         CALL uspp_augmentation_fourier(augmentation, gvecs(:, ig), position, volume, qg)
         rho(ig) = rho(ig) + SUM(density_matrix*qg)
      END DO
   END SUBROUTINE uspp_augmentation_collocate

! **************************************************************************************************
!> \brief Adjoint of augmentation collocation: integrate a real local potential against Q.
!> \param augmentation Prepared data.
!> \param gvecs Cartesian reciprocal vectors (3,ng).
!> \param position Center in bohr.
!> \param volume Cell volume.
!> \param vg Fourier potential coefficients (integral divided by volume).
!> \param multiplicity One for a complete reciprocal list; two for a retained +/-G pair.
!> \param matrix Real symmetric screened projector contribution, overwritten.
!> \param derivative Optional derivatives of matrix with respect to the augmentation center.
!> \param strain Optional density strain contribution at fixed Fourier potential and volume weight.
! **************************************************************************************************
   SUBROUTINE uspp_augmentation_integrate(augmentation, gvecs, position, volume, vg, multiplicity, matrix, derivative, strain)
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation
      REAL(KIND=dp), INTENT(IN)                          :: gvecs(:, :), position(3), volume
      COMPLEX(KIND=dp), INTENT(IN)                       :: vg(:)
      REAL(KIND=dp), INTENT(IN)                          :: multiplicity(:)
      REAL(KIND=dp), INTENT(OUT)                         :: matrix(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: derivative(:, :, :), strain(:, :, :, :)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: dq(:, :, :, :), qg(:, :)
      INTEGER                                            :: idir, ig, jdir

      CPASSERT(SIZE(gvecs, 1) == 3 .AND. SIZE(gvecs, 2) == SIZE(vg))
      CPASSERT(SIZE(multiplicity) == SIZE(vg))
      CPASSERT(ALL(SHAPE(matrix) == augmentation%nproj))
      ALLOCATE (qg(augmentation%nproj, augmentation%nproj))
      matrix = 0.0_dp
      IF (PRESENT(derivative)) THEN
         CPASSERT(ALL(SHAPE(derivative) == [augmentation%nproj, augmentation%nproj, 3]))
         derivative = 0.0_dp
      END IF
      IF (PRESENT(strain)) THEN
         CPASSERT(ALL(SHAPE(strain) == [augmentation%nproj, augmentation%nproj, 3, 3]))
         ALLOCATE (dq(augmentation%nproj, augmentation%nproj, 3, 3))
         strain = 0.0_dp
      END IF
      DO ig = 1, SIZE(vg)
         CALL uspp_augmentation_fourier(augmentation, gvecs(:, ig), position, volume, qg, dq)
         IF (PRESENT(strain)) THEN
            DO jdir = 1, 3
               DO idir = 1, 3
                  strain(:, :, idir, jdir) = strain(:, :, idir, jdir) + volume*multiplicity(ig)* &
                                             REAL(CONJG(dq(:, :, idir, jdir))*vg(ig), KIND=dp)
               END DO
            END DO
         END IF
         matrix(:, :) = matrix + volume*multiplicity(ig)*REAL(CONJG(qg)*vg(ig), KIND=dp)
         IF (PRESENT(derivative)) THEN
            DO idir = 1, 3
               derivative(:, :, idir) = derivative(:, :, idir) + volume*multiplicity(ig)* &
                                        REAL(CMPLX(0.0_dp, gvecs(idir, ig), dp)*CONJG(qg)*vg(ig), KIND=dp)
            END DO
         END IF
      END DO
   END SUBROUTINE uspp_augmentation_integrate

END MODULE uspp_augmentation
