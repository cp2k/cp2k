!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \brief Fourier projection of the Gaussian ionic charge and its geometry derivatives.
MODULE qs_ionic_charge
   USE cp_control_types,                ONLY: qs_control_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_real_fourier_aliases
   USE pw_types,                        ONLY: pw_c1d_gs_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: uses_augmented_gpw, ionic_charge_fourier, ionic_charge_geometry

CONTAINS

! **************************************************************************************************
!> \brief Identify the GPW grid used for USPP or PAW densities.
!> \param control Quickstep method selection.
!> \param kinds Gaussian kinds in the system.
!> \return Whether GPW uses an ultrasoft or PAW dataset.
! **************************************************************************************************
   LOGICAL FUNCTION uses_augmented_gpw(control, kinds) RESULT(enabled)
      TYPE(qs_control_type), INTENT(IN)                  :: control
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)

      INTEGER                                            :: ikind

      enabled = .FALSE.
      IF (.NOT. control%gpw) RETURN
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(ikind)%sgp_potential)) CYCLE
         IF (.NOT. ALLOCATED(kinds(ikind)%sgp_potential%upf_data)) CYCLE
         ASSOCIATE (pot => kinds(ikind)%sgp_potential%upf_data)
            enabled = pot%is_ultrasoft .OR. pot%is_paw
         END ASSOCIATE
         IF (enabled) RETURN
      END DO
   END FUNCTION uses_augmented_gpw

! **************************************************************************************************
!> \brief Evaluate the physical wave vectors represented by one real FFT coefficient.
!> \param position Cartesian atomic center.
!> \param alpha Gaussian exponent.
!> \param charge Signed integrated Gaussian charge divided by the cell volume.
!> \param vectors Physical wave vectors represented by the current coefficient.
!> \param count Number of physical wave vectors for this coefficient.
!> \param terms Fourier coefficients including the averaging weights.
! **************************************************************************************************
   SUBROUTINE ionic_fourier_terms(position, alpha, charge, vectors, count, terms)
      REAL(dp), INTENT(IN)                               :: position(3), alpha, charge, vectors(3, 8)
      INTEGER, INTENT(IN)                                :: count
      COMPLEX(dp), INTENT(OUT)                           :: terms(8)

      INTEGER                                            :: j
      REAL(dp)                                           :: g(3)

      DO j = 1, count
         g = vectors(:, j)
         terms(j) = charge/REAL(count, dp)*EXP(CMPLX(-DOT_PRODUCT(g, g)/(4.0_dp*alpha), -DOT_PRODUCT(g, position), dp))
      END DO
   END SUBROUTINE ionic_fourier_terms

! **************************************************************************************************
!> \brief Project the negative ionic Gaussian density onto the auxiliary FFT cell.
!>        The Poisson kernel supplies the physical boundary conditions.
!> \param atoms Atomic centers contributing to this field.
!> \param kinds Native Gaussian charge definitions.
!> \param density Output Fourier coefficients.
!> \param axis Optional Cartesian derivative for translating the supplied centers.
! **************************************************************************************************
   SUBROUTINE ionic_charge_fourier(atoms, kinds, density, axis)
      TYPE(particle_type), INTENT(IN)                    :: atoms(:)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: density
      INTEGER, INTENT(IN), OPTIONAL                      :: axis

      COMPLEX(dp)                                        :: terms(8)
      INTEGER                                            :: atom, count, i, kind
      REAL(dp)                                           :: coefficient, vectors(3, 8), &
                                                            charge(SIZE(atoms)), alpha(SIZE(atoms))
      TYPE(pw_grid_type), POINTER                        :: grid

      IF (PRESENT(axis)) THEN
         CPASSERT(axis >= 1 .AND. axis <= 3)
      END IF
      grid => density%pw_grid
      charge = 0.0_dp
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), alpha_core_charge=alpha(atom), ccore_charge=coefficient)
         IF (alpha(atom) == 0.0_dp .OR. coefficient == 0.0_dp) CYCLE
         charge(atom) = -coefficient*(pi/alpha(atom))**1.5_dp/grid%vol
      END DO
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atoms, grid, alpha, charge, axis, density) PRIVATE(atom, i, vectors, terms, count)
      DO i = 1, SIZE(density%array)
         density%array(i) = 0.0_dp
         CALL pw_real_fourier_aliases(grid, i, vectors, count)
         DO atom = 1, SIZE(atoms)
            IF (charge(atom) == 0.0_dp) CYCLE
            CALL ionic_fourier_terms(atoms(atom)%r, alpha(atom), charge(atom), vectors, count, terms)
            IF (PRESENT(axis)) terms(:count) = CMPLX(0.0_dp, -vectors(axis, :count), dp)*terms(:count)
            density%array(i) = density%array(i) + SUM(terms(:count))
         END DO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE ionic_charge_fourier

! **************************************************************************************************
!> \brief Differentiate the ionic field at fixed fractional centers for strain.
!> \param atoms Atomic centers.
!> \param kinds Native Gaussian charge definitions.
!> \param potential Unweighted reciprocal potential conjugate to the ionic field.
!> \param forces Local negative Cartesian energy gradients.
!> \param virial Local negative strain derivative, excluding the integral measure.
!> \param energies Optional local per-atom contractions with the potential.
! **************************************************************************************************
   SUBROUTINE ionic_charge_geometry(atoms, kinds, potential, forces, virial, energies)
      TYPE(particle_type), INTENT(IN)                    :: atoms(:)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(OUT)                              :: forces(:, :), virial(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: energies(:)

      COMPLEX(dp)                                        :: covector, term, terms(8)
      INTEGER                                            :: atom, axis, count, i, j, kind
      REAL(dp)                                           :: alpha, charge, coefficient, g(3), &
                                                            multiplicity, vectors(3, 8)
      TYPE(pw_grid_type), POINTER                        :: grid

      grid => potential%pw_grid
      forces = 0.0_dp
      virial = 0.0_dp
      IF (PRESENT(energies)) energies = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atoms, kinds, potential, grid, forces, energies) &
!$OMP PRIVATE(atom, axis, count, i, j, kind, alpha, charge, coefficient, g, multiplicity, vectors, covector, term, terms) &
!$OMP REDUCTION(+:virial)
      DO atom = 1, SIZE(atoms)
         kind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(kind), alpha_core_charge=alpha, ccore_charge=coefficient)
         IF (alpha == 0.0_dp .OR. coefficient == 0.0_dp) CYCLE
         charge = -coefficient*(pi/alpha)**1.5_dp/grid%vol
         DO i = 1, SIZE(potential%array)
            multiplicity = 1.0_dp
            IF (grid%grid_span == HALFSPACE .AND. .NOT. (grid%have_g0 .AND. i == 1)) multiplicity = 2.0_dp
            covector = grid%vol*multiplicity*CONJG(potential%array(i))
            CALL pw_real_fourier_aliases(grid, i, vectors, count)
            CALL ionic_fourier_terms(atoms(atom)%r, alpha, charge, vectors, count, terms)
            DO j = 1, count
               g = vectors(:, j)
               term = covector*terms(j)
               forces(:, atom) = forces(:, atom) - AIMAG(term)*g
               DO axis = 1, 3
                  virial(:, axis) = virial(:, axis) - REAL(term, dp)*g*g(axis)/(2.0_dp*alpha)
                  virial(axis, axis) = virial(axis, axis) + REAL(term, dp)
               END DO
               IF (PRESENT(energies)) energies(atom) = energies(atom) + REAL(term, dp)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE ionic_charge_geometry

END MODULE qs_ionic_charge
