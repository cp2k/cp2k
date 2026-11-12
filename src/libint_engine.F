!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief C-interoperable bindings for the Libint C++ Engine.
!> \note The C ABI, ownership, normalization, and result layout are documented in
!>       libint_engine_c_api.h. Each thread owns its handles and result scratch.
!>       Mixed radial operators are composed behind the existing integral interfaces.
! **************************************************************************************************
MODULE libint_engine
   USE ISO_C_BINDING,                   ONLY: &
        C_ASSOCIATED, C_CHAR, C_DOUBLE, C_F_POINTER, C_INT, C_NULL_CHAR, C_NULL_PTR, C_PTR, &
        C_SIZE_T
   USE input_constants,                 ONLY: &
        do_potential_coulomb, do_potential_gaussian, do_potential_id, do_potential_long, &
        do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, do_potential_short, &
        do_potential_truncated
   USE integral_library_types,          ONLY: compare_potential_types,&
                                              coulomb_operator_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_engine'

!> \brief Thread-owned Engines and scratch for a fixed integral rank and derivative order.
!> \note Pass by reference; intrinsic copies share owning C handles and scratch.
   TYPE, PUBLIC :: cp_libint_engine_context_type
      PRIVATE
      TYPE(C_PTR) :: engines(2) = C_NULL_PTR
      TYPE(coulomb_operator_type) :: potential = coulomb_operator_type()
      INTEGER :: ncenters = 0, derivative_order = 0, max_l = 0, max_nprim = 0, nterms = 0
      REAL(dp) :: weights(2) = 0.0_dp
      REAL(dp), POINTER :: values(:, :) => NULL()
   END TYPE cp_libint_engine_context_type

   INTEGER(C_INT), PARAMETER, PUBLIC :: libint_engine_coulomb = 1, &
                                        libint_engine_erf_coulomb = 2, &
                                        libint_engine_erfc_coulomb = 3, &
                                        libint_engine_truncated_coulomb = 4, &
                                        libint_engine_gaussian = 5, &
                                        libint_engine_delta = 6, &
                                        libint_engine_erfx_coulomb = 7

   PUBLIC :: cp_libint_engine_initialize, cp_libint_engine_finalize, cp_libint_engine_create, &
             cp_libint_engine_set_shell, cp_libint_engine_compute, cp_libint_engine_destroy
   PUBLIC :: cp_libint_engine_init, cp_libint_engine_release, cp_libint_engine_evaluate

   INTERFACE
      !> Called by CP2K's process-level Libint initialization and finalization.
      FUNCTION cp_libint_engine_initialize(error, error_capacity) RESULT(status) &
         BIND(C, NAME="cp2k_libint_engine_initialize")
         IMPORT :: C_CHAR, C_INT, C_SIZE_T
         CHARACTER(KIND=C_CHAR), INTENT(OUT) :: error(*)
         INTEGER(C_SIZE_T), VALUE :: error_capacity
         INTEGER(C_INT) :: status
      END FUNCTION cp_libint_engine_initialize

      FUNCTION cp_libint_engine_finalize(error, error_capacity) RESULT(status) &
         BIND(C, NAME="cp2k_libint_engine_finalize")
         IMPORT :: C_CHAR, C_INT, C_SIZE_T
         CHARACTER(KIND=C_CHAR), INTENT(OUT) :: error(*)
         INTEGER(C_SIZE_T), VALUE :: error_capacity
         INTEGER(C_INT) :: status
      END FUNCTION cp_libint_engine_finalize

      !> handle must initially be C_NULL_PTR; destruction returns it to C_NULL_PTR.
      FUNCTION cp_libint_engine_create(operator_id, ncenters, derivative_order, max_nprim, max_l, &
                                       parameters, target_precision, handle, error, error_capacity) RESULT(status) &
         BIND(C, NAME="cp2k_libint_engine_create")
         IMPORT :: C_CHAR, C_DOUBLE, C_INT, C_PTR, C_SIZE_T
         INTEGER(C_INT), VALUE :: operator_id, ncenters, derivative_order, max_nprim, max_l
         REAL(C_DOUBLE), INTENT(IN) :: parameters(3)
         REAL(C_DOUBLE), VALUE :: target_precision
         TYPE(C_PTR), INTENT(INOUT) :: handle
         CHARACTER(KIND=C_CHAR), INTENT(OUT) :: error(*)
         INTEGER(C_SIZE_T), VALUE :: error_capacity
         INTEGER(C_INT) :: status
      END FUNCTION cp_libint_engine_create

      FUNCTION cp_libint_engine_set_shell(handle, slot, l, nprim, exponents, coefficients, center, &
                                          error, error_capacity) RESULT(status) &
         BIND(C, NAME="cp2k_libint_engine_set_shell")
         IMPORT :: C_CHAR, C_DOUBLE, C_INT, C_PTR, C_SIZE_T
         TYPE(C_PTR), VALUE :: handle
         INTEGER(C_INT), VALUE :: slot, l, nprim
         REAL(C_DOUBLE), INTENT(IN) :: exponents(*), coefficients(*), center(3)
         CHARACTER(KIND=C_CHAR), INTENT(OUT) :: error(*)
         INTEGER(C_SIZE_T), VALUE :: error_capacity
         INTEGER(C_INT) :: status
      END FUNCTION cp_libint_engine_set_shell

      !> Map targets to ncomponents C_PTR entries with C_F_POINTER, then map each
      !> non-null entry to ncart REAL(C_DOUBLE) elements. Treat these views as read-only.
      FUNCTION cp_libint_engine_compute(handle, targets, ncart, ncomponents, error, error_capacity) RESULT(status) &
         BIND(C, NAME="cp2k_libint_engine_compute")
         IMPORT :: C_CHAR, C_INT, C_PTR, C_SIZE_T
         TYPE(C_PTR), VALUE :: handle
         TYPE(C_PTR), INTENT(OUT) :: targets
         INTEGER(C_SIZE_T), INTENT(OUT) :: ncart
         INTEGER(C_INT), INTENT(OUT) :: ncomponents
         CHARACTER(KIND=C_CHAR), INTENT(OUT) :: error(*)
         INTEGER(C_SIZE_T), VALUE :: error_capacity
         INTEGER(C_INT) :: status
      END FUNCTION cp_libint_engine_compute

      SUBROUTINE cp_libint_engine_destroy(handle) BIND(C, NAME="cp2k_libint_engine_destroy")
         IMPORT :: C_PTR
         TYPE(C_PTR), INTENT(INOUT) :: handle
      END SUBROUTINE cp_libint_engine_destroy
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Record the capacities of a fresh, owning Engine context.
!> \param context uninitialized owning context
!> \param ncenters number of physical centers, two through four
!> \param max_l maximum angular momentum
!> \param derivative_order zero for values, one for center derivatives
! **************************************************************************************************
   SUBROUTINE cp_libint_engine_init(context, ncenters, max_l, derivative_order)
      TYPE(cp_libint_engine_context_type), INTENT(OUT)   :: context
      INTEGER, INTENT(IN)                                :: ncenters, max_l, derivative_order

      context = cp_libint_engine_context_type()
      context%ncenters = ncenters
      context%max_l = max_l
      context%derivative_order = derivative_order
   END SUBROUTINE cp_libint_engine_init

! **************************************************************************************************
!> \brief Release the handles and scratch of their original owner.
!> \param context owning context
! **************************************************************************************************
   SUBROUTINE cp_libint_engine_release(context)
      TYPE(cp_libint_engine_context_type), INTENT(INOUT) :: context

      INTEGER                                            :: iterm

      DO iterm = 1, SIZE(context%engines)
         CALL cp_libint_engine_destroy(context%engines(iterm))
      END DO
      IF (ASSOCIATED(context%values)) DEALLOCATE (context%values)
      context%nterms = 0
      context%max_nprim = 0
   END SUBROUTINE cp_libint_engine_release

! **************************************************************************************************
!> \brief Evaluate Cartesian shells, composing mixed operators in reusable scratch.
!> \param context thread-owned Engine context
!> \param angular_momenta physical shell angular momenta in Engine order
!> \param centers physical shell centers
!> \param nprim number of primitives in each shell
!> \param exponents primitive exponents, padded to MAXVAL(nprim)
!> \param coefficients coefficients of normalization-free Cartesian primitives
!> \param potential radial operator; TShPSC uses its separate legacy recurrence
!> \param values borrowed result with last physical shell index fastest; derivatives
!>        are ordered by physical center, then x/y/z. Valid until the next evaluation or release.
!> \note Internal primitive screening is disabled; existing callers retain their screening policy.
! **************************************************************************************************
   SUBROUTINE cp_libint_engine_evaluate(context, angular_momenta, centers, nprim, exponents, coefficients, potential, values)
      TYPE(cp_libint_engine_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN)                                :: angular_momenta(:)
      REAL(dp), INTENT(IN)                               :: centers(:, :)
      INTEGER, INTENT(IN)                                :: nprim(:)
      REAL(dp), INTENT(IN)                               :: exponents(:, :), coefficients(:, :)
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential
      REAL(dp), INTENT(OUT), POINTER                     :: values(:, :)

      CHARACTER(KIND=C_CHAR, LEN=512)                    :: error
      INTEGER                                            :: component, iterm, ncart, ncomponents
      INTEGER(C_INT)                                     :: actual_components, operators(2), status
      INTEGER(C_SIZE_T)                                  :: actual_size
      REAL(C_DOUBLE)                                     :: parameters(3, 2)
      REAL(C_DOUBLE), POINTER                            :: component_values(:)
      TYPE(C_PTR)                                        :: targets
      TYPE(C_PTR), POINTER                               :: component_ptrs(:)

      CPASSERT(SIZE(angular_momenta) == context%ncenters)
      CPASSERT(SIZE(nprim) == context%ncenters)
      IF (context%nterms == 0 .OR. MAXVAL(nprim) > context%max_nprim .OR. &
          MAXVAL(angular_momenta) > context%max_l .OR. &
          .NOT. compare_potential_types(potential, context%potential)) THEN
         DO iterm = 1, SIZE(context%engines)
            CALL cp_libint_engine_destroy(context%engines(iterm))
         END DO
         context%potential = potential
         context%max_l = MAX(context%max_l, MAXVAL(angular_momenta))
         context%max_nprim = MAX(context%max_nprim, MAXVAL(nprim))
         context%nterms = 1
         context%weights = [1.0_dp, 0.0_dp]
         parameters = 0.0_dp
         SELECT CASE (potential%potential_type)
         CASE (do_potential_coulomb)
            operators(1) = libint_engine_coulomb
         CASE (do_potential_truncated)
            operators(1) = libint_engine_truncated_coulomb
            parameters(1, 1) = potential%cutoff_radius
         CASE (do_potential_short, do_potential_long)
            operators(1) = libint_engine_erfc_coulomb
            IF (potential%potential_type == do_potential_long) operators(1) = libint_engine_erf_coulomb
            parameters(1, 1) = potential%omega
         CASE (do_potential_id)
            operators(1) = libint_engine_delta
         CASE (do_potential_gaussian)
            operators(1) = libint_engine_gaussian
            parameters(1, 1) = potential%omega**2
         CASE (do_potential_mix_cl)
            operators(1) = libint_engine_erfx_coulomb
            parameters(:, 1) = [potential%omega, potential%scale_coulomb + potential%scale_longrange, &
                                potential%scale_coulomb]
         CASE (do_potential_mix_cl_trunc)
            context%nterms = 2
            operators = [libint_engine_truncated_coulomb, libint_engine_erfc_coulomb]
            parameters(1, :) = [potential%cutoff_radius, potential%omega]
            context%weights = [potential%scale_coulomb + potential%scale_longrange, -potential%scale_longrange]
         CASE (do_potential_mix_lg)
            context%nterms = 2
            operators = [libint_engine_erf_coulomb, libint_engine_gaussian]
            parameters(1, :) = [potential%omega, potential%omega**2]
            context%weights = [potential%scale_longrange, potential%scale_gaussian]
         CASE DEFAULT
            CPABORT("Operator is unavailable through the Libint Engine interface")
         END SELECT
         DO iterm = 1, context%nterms
            status = cp_libint_engine_create(operators(iterm), INT(context%ncenters, C_INT), &
                                             INT(context%derivative_order, C_INT), INT(context%max_nprim, C_INT), &
                                             INT(context%max_l, C_INT), parameters(:, iterm), 0.0_C_DOUBLE, &
                                             context%engines(iterm), error, INT(LEN(error), C_SIZE_T))
            IF (status /= 0) THEN
               CPABORT("Libint Engine creation failed: "//error(:INDEX(error, C_NULL_CHAR) - 1))
            END IF
         END DO
      END IF

      ncart = PRODUCT((angular_momenta + 1)*(angular_momenta + 2)/2)
      ncomponents = 1
      IF (context%derivative_order == 1) ncomponents = 3*context%ncenters
      IF (ASSOCIATED(context%values)) THEN
         IF (SIZE(context%values, 1) < ncart) DEALLOCATE (context%values)
      END IF
      IF (.NOT. ASSOCIATED(context%values)) ALLOCATE (context%values(ncart, ncomponents))
      values => context%values(1:ncart, :)
      values = 0.0_dp
      DO iterm = 1, context%nterms
         DO component = 1, context%ncenters
            status = cp_libint_engine_set_shell(context%engines(iterm), INT(component, C_INT), &
                                                INT(angular_momenta(component), C_INT), INT(nprim(component), C_INT), &
                                                exponents(:, component), coefficients(:, component), centers(:, component), &
                                                error, INT(LEN(error), C_SIZE_T))
            IF (status /= 0) THEN
               CPABORT("Libint shell setup failed: "//error(:INDEX(error, C_NULL_CHAR) - 1))
            END IF
         END DO
         status = cp_libint_engine_compute(context%engines(iterm), targets, actual_size, actual_components, &
                                           error, INT(LEN(error), C_SIZE_T))
         IF (status /= 0) THEN
            CPABORT("Libint Engine evaluation failed: "//error(:INDEX(error, C_NULL_CHAR) - 1))
         END IF
         CPASSERT(actual_size == INT(ncart, C_SIZE_T))
         CPASSERT(actual_components == ncomponents)
         IF (.NOT. C_ASSOCIATED(targets)) CYCLE
         CALL C_F_POINTER(targets, component_ptrs, [ncomponents])
         DO component = 1, ncomponents
            IF (.NOT. C_ASSOCIATED(component_ptrs(component))) CYCLE
            CALL C_F_POINTER(component_ptrs(component), component_values, [ncart])
            values(:, component) = values(:, component) + context%weights(iterm)*component_values
         END DO
      END DO
   END SUBROUTINE cp_libint_engine_evaluate
END MODULE libint_engine
