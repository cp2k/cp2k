!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief C-interoperable bindings for the Libint C++ Engine.
!> \note The C ABI, ownership, normalization, and result layout are documented in
!>       libint_engine_c_api.h. Each thread owns its handle. Production integral
!>       callers continue to use libint_wrapper until explicitly migrated.
! **************************************************************************************************
MODULE libint_engine
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_PTR,&
                                              C_SIZE_T

   IMPLICIT NONE
   PRIVATE

   INTEGER(C_INT), PARAMETER, PUBLIC :: libint_engine_coulomb = 1, &
                                        libint_engine_erf_coulomb = 2, &
                                        libint_engine_erfc_coulomb = 3, &
                                        libint_engine_truncated_coulomb = 4, &
                                        libint_engine_gaussian = 5, &
                                        libint_engine_delta = 6, &
                                        libint_engine_erfx_coulomb = 7

   PUBLIC :: cp_libint_engine_initialize, cp_libint_engine_finalize, cp_libint_engine_create, &
             cp_libint_engine_set_shell, cp_libint_engine_compute, cp_libint_engine_destroy

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
END MODULE libint_engine
