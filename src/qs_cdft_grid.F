!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Pointwise CDFT partition functions for nonuniform integration grids.
! **************************************************************************************************
MODULE qs_cdft_grid
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: hirshfeld_type
   USE input_constants,                 ONLY: becke_cutoff_element,&
                                              becke_cutoff_global,&
                                              outer_scf_becke_constraint,&
                                              outer_scf_hirshfeld_constraint
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_cdft_types,                   ONLY: becke_constraint_type,&
                                              cdft_control_type,&
                                              hirshfeld_constraint_type
   USE qs_cdft_utils,                   ONLY: hirshfeld_constraint_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_grid'

   TYPE cdft_point_context_type
      INTEGER                                            :: method = -1, natom = 0, ngroup = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: numexp, cavity_numexp
      LOGICAL                                            :: calculate_derivatives = .FALSE., &
                                                            cavity_confine = .FALSE.
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: constraint_atom
      REAL(KIND=dp)                                      :: eps = 0.0_dp, eps_cavity = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoffs, distances, cell_functions
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords, aij, coefficients, &
                                                            alpha, amplitude, cavity_alpha, &
                                                            cavity_amplitude, displacement, &
                                                            datom_numerator, datom_sum, dcell_point, &
                                                            density_atom_derivative
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: pair_vectors, dcell_atom
      TYPE(cell_type), POINTER                           :: cell => NULL()
   END TYPE cdft_point_context_type

   PUBLIC :: cdft_point_context_create, &
             cdft_point_context_release, &
             cdft_point_context_type, &
             cdft_point_weights

CONTAINS

! **************************************************************************************************
!> \brief Initialize reusable data for pointwise CDFT partition evaluation.
!> \param qs_env Quickstep environment
!> \param context pointwise partition context
!> \param calculate_derivatives allocate scratch space for coordinate derivatives
! **************************************************************************************************
   SUBROUTINE cdft_point_context_create(qs_env, context, calculate_derivatives)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cdft_point_context_type), INTENT(OUT)         :: context
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_derivatives

      INTEGER                                            :: atom, iatom, igroup, ikind, jatom, natom
      REAL(KIND=dp)                                      :: chi, ircov, jrcov, uij
      REAL(KIND=dp), DIMENSION(3)                        :: pair_vector
      REAL(KIND=dp), DIMENSION(:), POINTER               :: radii, radii_list
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(becke_constraint_type), POINTER               :: becke_control
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_constraint_type), POINTER           :: hirshfeld_control
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (becke_control, cdft_control, dft_control, hirshfeld_control, &
               hirshfeld_env, atomic_kind_set, particle_set, qs_kind_set, radii, radii_list)
      CALL get_qs_env(qs_env, cell=context%cell, dft_control=dft_control, &
                      natom=natom, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set)
      CPASSERT(ASSOCIATED(context%cell))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      cdft_control => dft_control%qs_control%cdft_control
      CPASSERT(ASSOCIATED(cdft_control))

      context%method = cdft_control%type
      context%natom = natom
      context%ngroup = SIZE(cdft_control%group)
      context%calculate_derivatives = .FALSE.
      IF (PRESENT(calculate_derivatives)) context%calculate_derivatives = calculate_derivatives
      ALLOCATE (context%atom_coords(3, natom), &
                context%coefficients(context%ngroup, natom), &
                context%constraint_atom(natom), &
                context%distances(natom), context%displacement(3, natom), &
                context%cell_functions(natom))
      context%coefficients = 0.0_dp
      context%constraint_atom = .FALSE.
      DO atom = 1, natom
         context%atom_coords(:, atom) = particle_set(atom)%r
      END DO
      DO igroup = 1, context%ngroup
         DO iatom = 1, SIZE(cdft_control%group(igroup)%atoms)
            atom = cdft_control%group(igroup)%atoms(iatom)
            context%coefficients(igroup, atom) = cdft_control%group(igroup)%coeff(iatom)
            context%constraint_atom(atom) = .TRUE.
         END DO
      END DO

      SELECT CASE (context%method)
      CASE (outer_scf_becke_constraint)
         becke_control => cdft_control%becke_control
         CPASSERT(ASSOCIATED(becke_control))
         ALLOCATE (context%cutoffs(natom), context%aij(natom, natom), &
                   context%pair_vectors(3, natom, natom))
         IF (context%calculate_derivatives) THEN
            ALLOCATE (context%dcell_point(3, natom), context%dcell_atom(3, natom, natom), &
                      context%datom_sum(3, natom), context%datom_numerator(3, natom))
         END IF
         IF (ASSOCIATED(becke_control%cutoffs)) THEN
            context%cutoffs(:) = becke_control%cutoffs
         ELSE
            SELECT CASE (becke_control%cutoff_type)
            CASE (becke_cutoff_global)
               context%cutoffs = becke_control%rglobal
            CASE (becke_cutoff_element)
               CPASSERT(ASSOCIATED(becke_control%cutoffs_tmp))
               CPASSERT(SIZE(becke_control%cutoffs_tmp) == SIZE(atomic_kind_set))
               DO atom = 1, natom
                  CALL get_atomic_kind(particle_set(atom)%atomic_kind, kind_number=ikind)
                  context%cutoffs(atom) = becke_control%cutoffs_tmp(ikind)
               END DO
            CASE DEFAULT
               CPABORT("Unknown Becke cutoff type.")
            END SELECT
         END IF
         context%aij = 0.0_dp
         IF (becke_control%adjust) THEN
            IF (ASSOCIATED(becke_control%aij)) THEN
               context%aij(:, :) = becke_control%aij
            ELSE
               IF (ASSOCIATED(becke_control%radii)) THEN
                  radii => becke_control%radii
               ELSE
                  radii => becke_control%radii_tmp
               END IF
               CPASSERT(ASSOCIATED(radii))
               CPASSERT(SIZE(radii) == SIZE(atomic_kind_set))
               DO iatom = 1, natom - 1
                  CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                  ircov = radii(ikind)
                  DO jatom = iatom + 1, natom
                     CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=ikind)
                     jrcov = radii(ikind)
                     IF (ircov /= jrcov) THEN
                        chi = ircov/jrcov
                        uij = (chi - 1.0_dp)/(chi + 1.0_dp)
                        context%aij(iatom, jatom) = MAX(-0.5_dp, MIN(0.5_dp, &
                                                                     uij/(uij**2 - 1.0_dp)))
                        context%aij(jatom, iatom) = -context%aij(iatom, jatom)
                     END IF
                  END DO
               END DO
            END IF
         END IF
         context%pair_vectors = 0.0_dp
         DO iatom = 1, natom - 1
            DO jatom = iatom + 1, natom
               pair_vector = pbc(context%atom_coords(:, jatom), &
                                 context%atom_coords(:, iatom), context%cell)
               context%pair_vectors(:, iatom, jatom) = pair_vector
               context%pair_vectors(:, jatom, iatom) = -pair_vector
            END DO
         END DO
         context%cavity_confine = becke_control%cavity_confine
         context%eps_cavity = becke_control%eps_cavity
         IF (context%cavity_confine) THEN
            hirshfeld_env => becke_control%cavity_env
            CPASSERT(ASSOCIATED(hirshfeld_env))
            IF (.NOT. ASSOCIATED(hirshfeld_env%kind_shape_fn)) THEN
               IF (ASSOCIATED(becke_control%radii)) THEN
                  radii => becke_control%radii
               ELSE IF (ASSOCIATED(becke_control%radii_tmp)) THEN
                  radii => becke_control%radii_tmp
               END IF
               IF (ASSOCIATED(radii)) THEN
                  ALLOCATE (radii_list(SIZE(radii)))
                  DO ikind = 1, SIZE(radii)
                     IF (hirshfeld_env%use_bohr) THEN
                        radii_list(ikind) = radii(ikind)
                     ELSE
                        radii_list(ikind) = cp_unit_from_cp2k(radii(ikind), "angstrom")
                     END IF
                  END DO
               END IF
               CALL create_shape_function(hirshfeld_env, qs_kind_set, atomic_kind_set, &
                                          radius=becke_control%rcavity, radii_list=radii_list)
               IF (ASSOCIATED(radii_list)) DEALLOCATE (radii_list)
            END IF
            CALL store_shape_functions(hirshfeld_env, particle_set, context%cavity_numexp, &
                                       context%cavity_alpha, context%cavity_amplitude, &
                                       include_charge=.FALSE.)
         END IF
      CASE (outer_scf_hirshfeld_constraint)
         hirshfeld_control => cdft_control%hirshfeld_control
         CPASSERT(ASSOCIATED(hirshfeld_control))
         hirshfeld_env => hirshfeld_control%hirshfeld_env
         CPASSERT(ASSOCIATED(hirshfeld_env))
         IF (.NOT. ASSOCIATED(hirshfeld_env%kind_shape_fn) .OR. &
             .NOT. ASSOCIATED(hirshfeld_env%charges)) THEN
            CALL hirshfeld_constraint_init(qs_env)
         END IF
         context%eps = hirshfeld_control%eps_cutoff
         CALL store_shape_functions(hirshfeld_env, particle_set, context%numexp, &
                                    context%alpha, context%amplitude, include_charge=.TRUE.)
         IF (context%calculate_derivatives) THEN
            ALLOCATE (context%density_atom_derivative(3, natom))
         END IF
      CASE DEFAULT
         CPABORT("Unknown CDFT partition type.")
      END SELECT

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param environment ...
!> \param particles ...
!> \param numexp ...
!> \param alpha ...
!> \param amplitude ...
!> \param include_charge ...
! **************************************************************************************************
      SUBROUTINE store_shape_functions(environment, particles, numexp, alpha, amplitude, &
                                       include_charge)
      TYPE(hirshfeld_type), POINTER                      :: environment
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: numexp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: alpha, amplitude
      LOGICAL, INTENT(IN)                                :: include_charge

      INTEGER                                            :: atom, iexp, ikind, maxexp
      REAL(KIND=dp)                                      :: charge

         ALLOCATE (numexp(natom))
         numexp = 0
         maxexp = 0
         DO atom = 1, natom
            CALL get_atomic_kind(particles(atom)%atomic_kind, kind_number=ikind)
            numexp(atom) = environment%kind_shape_fn(ikind)%numexp
            maxexp = MAX(maxexp, numexp(atom))
         END DO
         ALLOCATE (alpha(maxexp, natom), amplitude(maxexp, natom))
         alpha = 0.0_dp
         amplitude = 0.0_dp
         DO atom = 1, natom
            CALL get_atomic_kind(particles(atom)%atomic_kind, kind_number=ikind)
            charge = 1.0_dp
            IF (include_charge) charge = environment%charges(atom)
            DO iexp = 1, numexp(atom)
               alpha(iexp, atom) = environment%kind_shape_fn(ikind)%zet(iexp)
               amplitude(iexp, atom) = charge*environment%kind_shape_fn(ikind)%coef(iexp)
            END DO
         END DO
      END SUBROUTINE store_shape_functions

   END SUBROUTINE cdft_point_context_create

! **************************************************************************************************
!> \brief Release a pointwise CDFT partition context.
!> \param context pointwise partition context
! **************************************************************************************************
   SUBROUTINE cdft_point_context_release(context)
      TYPE(cdft_point_context_type), INTENT(INOUT)       :: context

      IF (ALLOCATED(context%numexp)) DEALLOCATE (context%numexp)
      IF (ALLOCATED(context%cavity_numexp)) DEALLOCATE (context%cavity_numexp)
      IF (ALLOCATED(context%constraint_atom)) DEALLOCATE (context%constraint_atom)
      IF (ALLOCATED(context%cutoffs)) DEALLOCATE (context%cutoffs)
      IF (ALLOCATED(context%distances)) DEALLOCATE (context%distances)
      IF (ALLOCATED(context%cell_functions)) DEALLOCATE (context%cell_functions)
      IF (ALLOCATED(context%datom_sum)) DEALLOCATE (context%datom_sum)
      IF (ALLOCATED(context%atom_coords)) DEALLOCATE (context%atom_coords)
      IF (ALLOCATED(context%aij)) DEALLOCATE (context%aij)
      IF (ALLOCATED(context%coefficients)) DEALLOCATE (context%coefficients)
      IF (ALLOCATED(context%alpha)) DEALLOCATE (context%alpha)
      IF (ALLOCATED(context%amplitude)) DEALLOCATE (context%amplitude)
      IF (ALLOCATED(context%cavity_alpha)) DEALLOCATE (context%cavity_alpha)
      IF (ALLOCATED(context%cavity_amplitude)) DEALLOCATE (context%cavity_amplitude)
      IF (ALLOCATED(context%displacement)) DEALLOCATE (context%displacement)
      IF (ALLOCATED(context%datom_numerator)) DEALLOCATE (context%datom_numerator)
      IF (ALLOCATED(context%dcell_point)) DEALLOCATE (context%dcell_point)
      IF (ALLOCATED(context%density_atom_derivative)) DEALLOCATE (context%density_atom_derivative)
      IF (ALLOCATED(context%pair_vectors)) DEALLOCATE (context%pair_vectors)
      IF (ALLOCATED(context%dcell_atom)) DEALLOCATE (context%dcell_atom)
      NULLIFY (context%cell)
      context%method = -1
      context%natom = 0
      context%ngroup = 0
      context%calculate_derivatives = .FALSE.
   END SUBROUTINE cdft_point_context_release

! **************************************************************************************************
!> \brief Evaluate CDFT weights and coordinate derivatives at one point.
!> \param context pointwise partition context
!> \param point Cartesian point
!> \param weights group weights
!> \param point_derivative derivatives with respect to the Cartesian point
!> \param atom_derivative derivatives with respect to atom positions at fixed point
!> \param atomic_weights optional individual atomic partition weights
! **************************************************************************************************
   SUBROUTINE cdft_point_weights(context, point, weights, point_derivative, atom_derivative, &
                                 atomic_weights)
      TYPE(cdft_point_context_type), INTENT(INOUT)       :: context
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: point_derivative
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: atom_derivative
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: atomic_weights

      CPASSERT(SIZE(weights) == context%ngroup)
      CPASSERT(SIZE(point_derivative, 1) == 3)
      CPASSERT(SIZE(point_derivative, 2) == context%ngroup)
      CPASSERT(SIZE(atom_derivative, 1) == 3)
      CPASSERT(SIZE(atom_derivative, 2) == context%natom)
      CPASSERT(SIZE(atom_derivative, 3) == context%ngroup)
      IF (PRESENT(atomic_weights)) THEN
         CPASSERT(SIZE(atomic_weights) == context%natom)
      END IF
      SELECT CASE (context%method)
      CASE (outer_scf_becke_constraint)
         CALL becke_point_weights(context, point, weights, point_derivative, atom_derivative, &
                                  atomic_weights, context%calculate_derivatives)
      CASE (outer_scf_hirshfeld_constraint)
         CALL hirshfeld_point_weights(context, point, weights, point_derivative, atom_derivative, &
                                      atomic_weights, context%calculate_derivatives)
      CASE DEFAULT
         CPABORT("Unknown CDFT partition type.")
      END SELECT
   END SUBROUTINE cdft_point_weights

! **************************************************************************************************
!> \brief Evaluate Becke weights at one point.
!> \param context ...
!> \param point ...
!> \param weights ...
!> \param point_derivative ...
!> \param atom_derivative ...
!> \param atomic_weights ...
!> \param calculate_derivatives ...
! **************************************************************************************************
   SUBROUTINE becke_point_weights(context, point, weights, point_derivative, atom_derivative, &
                                  atomic_weights, calculate_derivatives)
      TYPE(cdft_point_context_type), INTENT(INOUT)       :: context
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: point_derivative
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: atom_derivative
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: atomic_weights
      LOGICAL, INTENT(IN)                                :: calculate_derivatives

      INTEGER                                            :: atom, iatom, iexp, igroup, jatom
      REAL(KIND=dp)                                      :: adjusted_mu, cavity_density, delta, &
                                                            dmu_factor, f1, f2, f3, mu, numerator, &
                                                            old_cell, pair_distance, s, sum_cell
      REAL(KIND=dp), DIMENSION(3)                        :: dmu_i, dmu_j, dmu_point, &
                                                            dpoint_numerator, ds_i, ds_j, &
                                                            ds_point, dsum_point, unit_i, unit_j

      weights = 0.0_dp
      IF (calculate_derivatives) THEN
         point_derivative = 0.0_dp
         atom_derivative = 0.0_dp
      END IF
      IF (PRESENT(atomic_weights)) atomic_weights = 0.0_dp

      IF (context%cavity_confine) THEN
         cavity_density = 0.0_dp
         DO atom = 1, context%natom
            IF (.NOT. context%constraint_atom(atom)) CYCLE
            context%displacement(:, atom) = pbc(context%atom_coords(:, atom), point, context%cell)
            DO iexp = 1, context%cavity_numexp(atom)
               cavity_density = cavity_density + context%cavity_amplitude(iexp, atom)* &
                                EXP(-context%cavity_alpha(iexp, atom)* &
                                    DOT_PRODUCT(context%displacement(:, atom), context%displacement(:, atom)))
            END DO
         END DO
         IF (cavity_density < context%eps_cavity) RETURN
      END IF

      context%cell_functions = 1.0_dp
      IF (calculate_derivatives) THEN
         context%dcell_point = 0.0_dp
         context%dcell_atom = 0.0_dp
      END IF
      DO atom = 1, context%natom
         context%displacement(:, atom) = pbc(context%atom_coords(:, atom), point, context%cell)
         context%distances(atom) = NORM2(context%displacement(:, atom))
      END DO
      DO iatom = 1, context%natom
         IF (context%distances(iatom) > context%cutoffs(iatom)) THEN
            context%cell_functions(iatom) = 0.0_dp
            CYCLE
         END IF
         IF (calculate_derivatives) THEN
            unit_i = 0.0_dp
            IF (context%distances(iatom) > 1.0E-14_dp) THEN
               unit_i = -context%displacement(:, iatom)/context%distances(iatom)
            END IF
         END IF
         DO jatom = 1, context%natom
            IF (jatom == iatom) CYCLE
            pair_distance = NORM2(context%pair_vectors(:, iatom, jatom))
            IF (pair_distance <= 1.0E-14_dp) CYCLE
            delta = context%distances(iatom) - context%distances(jatom)
            mu = delta/pair_distance
            adjusted_mu = mu + context%aij(iatom, jatom)*(1.0_dp - mu**2)
            f1 = 1.5_dp*adjusted_mu - 0.5_dp*adjusted_mu**3
            f2 = 1.5_dp*f1 - 0.5_dp*f1**3
            f3 = 1.5_dp*f2 - 0.5_dp*f2**3
            s = 0.5_dp*(1.0_dp - f3)
            old_cell = context%cell_functions(iatom)
            IF (calculate_derivatives) THEN
               dmu_factor = 1.0_dp - 2.0_dp*context%aij(iatom, jatom)*mu
               unit_j = 0.0_dp
               IF (context%distances(jatom) > 1.0E-14_dp) THEN
                  unit_j = -context%displacement(:, jatom)/context%distances(jatom)
               END IF
               dmu_i = unit_i/pair_distance - &
                       delta*context%pair_vectors(:, iatom, jatom)/pair_distance**3
               dmu_j = -unit_j/pair_distance + &
                       delta*context%pair_vectors(:, iatom, jatom)/pair_distance**3
               dmu_point = (-unit_i + unit_j)/pair_distance
               dmu_factor = -0.5_dp*dmu_factor*1.5_dp*(1.0_dp - adjusted_mu**2)* &
                            1.5_dp*(1.0_dp - f1**2)*1.5_dp*(1.0_dp - f2**2)
               ds_i = dmu_factor*dmu_i
               ds_j = dmu_factor*dmu_j
               ds_point = dmu_factor*dmu_point
               context%dcell_atom(:, :, iatom) = context%dcell_atom(:, :, iatom)*s
               context%dcell_point(:, iatom) = context%dcell_point(:, iatom)*s
               context%dcell_atom(:, iatom, iatom) = &
                  context%dcell_atom(:, iatom, iatom) + old_cell*ds_i
               context%dcell_atom(:, jatom, iatom) = &
                  context%dcell_atom(:, jatom, iatom) + old_cell*ds_j
               context%dcell_point(:, iatom) = context%dcell_point(:, iatom) + old_cell*ds_point
            END IF
            context%cell_functions(iatom) = old_cell*s
         END DO
      END DO

      sum_cell = SUM(context%cell_functions)
      IF (sum_cell <= 1.0E-6_dp) RETURN
      IF (PRESENT(atomic_weights)) atomic_weights = context%cell_functions/sum_cell
      IF (calculate_derivatives) THEN
         dsum_point = SUM(context%dcell_point, DIM=2)
         context%datom_sum(:, :) = SUM(context%dcell_atom, DIM=3)
      END IF
      DO igroup = 1, context%ngroup
         numerator = DOT_PRODUCT(context%coefficients(igroup, :), context%cell_functions)
         weights(igroup) = numerator/sum_cell
         IF (calculate_derivatives) THEN
            dpoint_numerator = MATMUL(context%dcell_point, context%coefficients(igroup, :))
            point_derivative(:, igroup) = &
               (dpoint_numerator*sum_cell - numerator*dsum_point)/sum_cell**2
            DO atom = 1, context%natom
               context%datom_numerator(:, atom) = &
                  MATMUL(context%dcell_atom(:, atom, :), context%coefficients(igroup, :))
               atom_derivative(:, atom, igroup) = &
                  (context%datom_numerator(:, atom)*sum_cell - &
                   numerator*context%datom_sum(:, atom))/sum_cell**2
            END DO
         END IF
      END DO
   END SUBROUTINE becke_point_weights

! **************************************************************************************************
!> \brief Evaluate Hirshfeld weights at one point.
!> \param context ...
!> \param point ...
!> \param weights ...
!> \param point_derivative ...
!> \param atom_derivative ...
!> \param atomic_weights ...
!> \param calculate_derivatives ...
! **************************************************************************************************
   SUBROUTINE hirshfeld_point_weights(context, point, weights, point_derivative, atom_derivative, &
                                      atomic_weights, calculate_derivatives)
      TYPE(cdft_point_context_type), INTENT(INOUT)       :: context
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: point_derivative
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: atom_derivative
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: atomic_weights
      LOGICAL, INTENT(IN)                                :: calculate_derivatives

      INTEGER                                            :: atom, iexp, igroup
      REAL(KIND=dp)                                      :: exponential, numerator, sum_density
      REAL(KIND=dp), DIMENSION(3)                        :: dpoint_numerator, dsum_point

      weights = 0.0_dp
      IF (calculate_derivatives) THEN
         point_derivative = 0.0_dp
         atom_derivative = 0.0_dp
      END IF
      IF (PRESENT(atomic_weights)) atomic_weights = 0.0_dp
      context%cell_functions = 0.0_dp
      IF (calculate_derivatives) context%density_atom_derivative = 0.0_dp
      DO atom = 1, context%natom
         context%displacement(:, atom) = pbc(context%atom_coords(:, atom), point, context%cell)
         DO iexp = 1, context%numexp(atom)
            exponential = context%amplitude(iexp, atom)* &
                          EXP(-context%alpha(iexp, atom)* &
                              DOT_PRODUCT(context%displacement(:, atom), context%displacement(:, atom)))
            context%cell_functions(atom) = context%cell_functions(atom) + exponential
            IF (calculate_derivatives) THEN
               context%density_atom_derivative(:, atom) = &
                  context%density_atom_derivative(:, atom) + &
                  2.0_dp*context%alpha(iexp, atom)*context%displacement(:, atom)*exponential
            END IF
         END DO
      END DO
      sum_density = SUM(context%cell_functions)
      IF (sum_density <= context%eps) THEN
         RETURN
      END IF
      IF (PRESENT(atomic_weights)) atomic_weights = context%cell_functions/sum_density
      IF (calculate_derivatives) dsum_point = -SUM(context%density_atom_derivative, DIM=2)
      DO igroup = 1, context%ngroup
         numerator = DOT_PRODUCT(context%coefficients(igroup, :), context%cell_functions)
         weights(igroup) = numerator/sum_density
         IF (calculate_derivatives) THEN
            dpoint_numerator(:) = &
               MATMUL(context%density_atom_derivative, context%coefficients(igroup, :))
            dpoint_numerator = -dpoint_numerator
            point_derivative(:, igroup) = &
               (dpoint_numerator*sum_density - numerator*dsum_point)/sum_density**2
            DO atom = 1, context%natom
               atom_derivative(:, atom, igroup) = &
                  (context%coefficients(igroup, atom) - weights(igroup))* &
                  context%density_atom_derivative(:, atom)/sum_density
            END DO
         END IF
      END DO
   END SUBROUTINE hirshfeld_point_weights

END MODULE qs_cdft_grid
