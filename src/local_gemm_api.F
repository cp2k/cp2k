!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE local_gemm_api
   USE ISO_C_BINDING, ONLY: C_NULL_PTR, &
                            C_PTR
   USE kinds, ONLY: dp
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
   USE input_constants, ONLY: do_dgemm_spla
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_LOC
   USE spla, ONLY: SPLA_OP_NONE, &
                   SPLA_OP_TRANSPOSE, &
                   SPLA_OP_CONJ_TRANSPOSE, &
                   spla_ctx_create, &
                   spla_ctx_destroy, &
                   spla_dgemm, &
                   spla_zgemm, &
                   spla_ctx_set_op_threshold_gpu, &
                   SPLA_SUCCESS
#endif

   USE cp_log_handling, ONLY: cp_to_string
   USE offload_api, ONLY: offload_activate_chosen_device

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'local_gemm_api'

   PUBLIC :: local_gemm_ctxt_type, &
             local_gemm_set_library

   INTEGER, PARAMETER, PUBLIC :: &
      LOCAL_GEMM_PU_HOST = 0, &
      LOCAL_GEMM_PU_GPU = 1

   INTEGER, PRIVATE :: do_dgemm = 1

   TYPE local_gemm_ctxt_type
      TYPE(C_PTR) :: spla_context = C_NULL_PTR
      LOGICAL, PRIVATE :: timing = .TRUE.
   CONTAINS
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: create => local_gemm_create
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: destroy => local_gemm_destroy
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: set_op_threshold_gpu => local_gemm_set_op_threshold_gpu
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE, PRIVATE :: gemm_d => local_dgemm
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE, PRIVATE :: gemm_z => local_zgemm
      GENERIC :: gemm => gemm_d, gemm_z
   END TYPE

CONTAINS

   #:for scalar, prefix in [('REAL', 'd'), ('COMPLEX', 'z')]
! **************************************************************************************************
!> \brief Local GEMM on contiguous arrays, using BLAS or the configured SPLA backend.
!>        Each concurrent caller must own its context. No distributed matrix metadata is used.
!> \param opA operation on A (N/T/C, case insensitive)
!> \param opB operation on B (N/T/C, case insensitive)
!> \param m output rows
!> \param n output columns
!> \param k contraction dimension
!> \param alpha product scale
!> \param A left operand
!> \param lda leading dimension of A
!> \param B right operand
!> \param ldb leading dimension of B
!> \param beta output scale
!> \param C output, must not overlap A or B
!> \param ldc leading dimension of C
!> \param ctx caller-owned context
! **************************************************************************************************
      SUBROUTINE local_${prefix}$gemm(opA, opB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, ctx)
         CHARACTER, INTENT(IN)                              :: opA, opB
         INTEGER, INTENT(IN)                                :: m, n, k, lda, ldb, ldc
         ${scalar}$ (KIND=dp), INTENT(IN)                    :: alpha, beta
         ${scalar}$ (KIND=dp), INTENT(IN), TARGET            :: A(lda, *), B(ldb, *)
         ${scalar}$ (KIND=dp), INTENT(INOUT), TARGET         :: C(ldc, *)
         CLASS(local_gemm_ctxt_type), INTENT(INOUT)        :: ctx

         INTEGER                                            :: handle
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
         INTEGER                                            :: spla_error
#endif
         CHARACTER(LEN=*), PARAMETER :: routineN = 'local_gemm'

         IF (m == 0 .OR. n == 0) RETURN
         IF (ctx%timing) CALL timeset(routineN, handle)
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
         IF (do_dgemm == do_dgemm_spla) THEN
            CPASSERT(C_ASSOCIATED(ctx%spla_context))
            CALL offload_activate_chosen_device()
            spla_error = spla_${prefix}$gemm(spla_operation(opA), spla_operation(opB), m, n, k, alpha, &
                                             C_LOC(A(1, 1)), lda, C_LOC(B(1, 1)), ldb, beta, &
                                             C_LOC(C(1, 1)), ldc, ctx%spla_context)
            IF (spla_error /= SPLA_SUCCESS) &
               CALL cp_abort(__LOCATION__, "spla_${prefix}$gemm failed: "//cp_to_string(spla_error))
         ELSE
#endif
            CALL ${prefix}$gemm(opA, opB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
         END IF
#endif
         IF (ctx%timing) CALL timestop(handle)

      END SUBROUTINE local_${prefix}$gemm
   #:endfor

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
! **************************************************************************************************
!> \brief Translate a BLAS transpose flag for SPLA, including complex conjugation.
!> \param trans BLAS operation
!> \return SPLA operation
! **************************************************************************************************
   FUNCTION spla_operation(trans) RESULT(op)
      CHARACTER, INTENT(IN)                              :: trans
      INTEGER                                            :: op

      SELECT CASE (trans)
      CASE ('N', 'n')
         op = SPLA_OP_NONE
      CASE ('T', 't')
         op = SPLA_OP_TRANSPOSE
      CASE ('C', 'c')
         op = SPLA_OP_CONJ_TRANSPOSE
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "Invalid local GEMM transpose flag.")
      END SELECT
   END FUNCTION spla_operation
#endif

! **************************************************************************************************
!> \brief Create a local GEMM context; destroy an existing context before recreating it.
!> \param ctx newly created context, with timing enabled by default
!> \param pu processing unit for local GEMM
!> \param timing collect per-GEMM timings (default true); disable for timed batches
! **************************************************************************************************
   SUBROUTINE local_gemm_create(ctx, pu, timing)
      CLASS(local_gemm_ctxt_type), INTENT(OUT)           :: ctx
      INTEGER, INTENT(IN)                                :: pu
      LOGICAL, INTENT(IN), OPTIONAL                      :: timing

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER                                            :: spla_error
#endif

      IF (PRESENT(timing)) ctx%timing = timing
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      IF (do_dgemm == do_dgemm_spla) THEN
         CALL offload_activate_chosen_device()

         spla_error = spla_ctx_create(ctx%spla_context, pu)
         IF (spla_error /= SPLA_SUCCESS) &
            CALL cp_abort(__LOCATION__, &
                          "spla_ctx_create failed: "//cp_to_string(spla_error))
      END IF
#else
      MARK_USED(pu)
#endif
   END SUBROUTINE local_gemm_create

! **************************************************************************************************
!> \brief Release an owned SPLA context, independently of the current backend preference.
!> \param ctx context to release; an empty context is allowed
! **************************************************************************************************
   SUBROUTINE local_gemm_destroy(ctx)
      CLASS(local_gemm_ctxt_type), INTENT(INOUT)         :: ctx

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER                                            :: spla_error

      IF (C_ASSOCIATED(ctx%spla_context)) THEN
         CALL offload_activate_chosen_device()

         spla_error = spla_ctx_destroy(ctx%spla_context)
         IF (spla_error /= SPLA_SUCCESS) &
            CALL cp_abort(__LOCATION__, &
                          "spla_ctx_destroy failed: "//cp_to_string(spla_error))
      END IF
#endif
      ctx%spla_context = C_NULL_PTR
   END SUBROUTINE local_gemm_destroy

! **************************************************************************************************
!> \brief Set the SPLA GPU operation threshold; no-op when no SPLA context is allocated.
!> \param ctx local GEMM context
!> \param opThresholdGPU operation-count threshold for GPU offloading
! **************************************************************************************************
   SUBROUTINE local_gemm_set_op_threshold_gpu(ctx, opThresholdGPU)
      CLASS(local_gemm_ctxt_type), INTENT(INOUT)         :: ctx
      INTEGER, INTENT(IN)                                :: opThresholdGPU

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER                                            :: spla_error

      IF (C_ASSOCIATED(ctx%spla_context)) THEN
         CALL offload_activate_chosen_device()

         spla_error = spla_ctx_set_op_threshold_gpu(ctx%spla_context, opThresholdGPU)
         IF (spla_error /= SPLA_SUCCESS) &
            CALL cp_abort(__LOCATION__, &
                          "spla_ctx_set_op_threshold_gpu failed: "//cp_to_string(spla_error))
      END IF
#else
      MARK_USED(ctx)
      MARK_USED(opThresholdGPU)
#endif
   END SUBROUTINE local_gemm_set_op_threshold_gpu

! **************************************************************************************************
!> \brief Select the backend for subsequent local GEMM calls and context creation.
!> \param dgemm_library backend selector from input_constants (SPLA or BLAS)
! **************************************************************************************************
   SUBROUTINE local_gemm_set_library(dgemm_library)
      INTEGER, INTENT(IN)                                :: dgemm_library

      do_dgemm = dgemm_library
   END SUBROUTINE local_gemm_set_library

END MODULE local_gemm_api
