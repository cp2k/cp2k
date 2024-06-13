!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE local_gemm_api
   USE ISO_C_BINDING, ONLY: C_NULL_PTR, &
                            C_PTR
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
   USE input_constants, ONLY: do_dgemm_spla
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_LOC
   USE spla, ONLY: SPLA_PU_HOST, &
                   SPLA_PU_GPU, &
                   SPLA_OP_NONE, &
                   SPLA_OP_TRANSPOSE, &
                   SPLA_OP_CONJ_TRANSPOSE, &
                   spla_ctx_create, &
                   spla_ctx_destroy, &
                   spla_dgemm, &
                   spla_sgemm, &
                   spla_cgemm, &
                   spla_zgemm, &
                   spla_ctx_set_op_threshold_gpu, &
                   SPLA_SUCCESS
#endif

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
   CONTAINS
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: create => local_gemm_create
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: destroy => local_gemm_destroy
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: set_op_threshold_gpu => local_gemm_set_op_threshold_gpu
      PROCEDURE, PASS(ctx), NON_OVERRIDABLE :: gemm => local_gemm
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param opA ...
!> \param opB ...
!> \param m ...
!> \param n ...
!> \param k ...
!> \param alpha ...
!> \param A ...
!> \param lda ...
!> \param B ...
!> \param ldb ...
!> \param beta ...
!> \param C ...
!> \param ldc ...
!> \param ctx ...
! **************************************************************************************************
   SUBROUTINE local_gemm(opA, opB, m, n, k, &
                         alpha, A, lda, B, ldb, &
                         beta, C, ldc, ctx)
      CHARACTER, INTENT(in) :: opA
      CHARACTER, INTENT(in) :: opB
      INTEGER, INTENT(in) :: m
      INTEGER, INTENT(in) :: n
      INTEGER, INTENT(in) :: k
      REAL(8), INTENT(in) :: alpha
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      REAL(8), DIMENSION(*), INTENT(in), TARGET :: A
#else
      REAL(8), DIMENSION(:, :), INTENT(in), TARGET :: A
#endif
      INTEGER, INTENT(in) :: lda
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      REAL(8), DIMENSION(*), INTENT(in), TARGET :: B
#else
      REAL(8), DIMENSION(:, :), INTENT(in), TARGET :: B
#endif

      INTEGER, INTENT(in) :: ldb
      REAL(8), INTENT(in) :: beta
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      REAL(8), DIMENSION(*), INTENT(inout), TARGET ::C
#else
      REAL(8), DIMENSION(:, :), INTENT(inout), TARGET :: C
#endif
      INTEGER, INTENT(in) :: ldc
      CLASS(local_gemm_ctxt_type), INTENT(inout) :: ctx

      INTEGER                                            :: handle
!     no point of using SPLA offloading on CPU ONLY nodes
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER :: spla_op_A, spla_op_B, spla_error
#endif
      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_gemm'
      CALL timeset(routineN, handle)

!     no point of using SPLA offloading on CPU ONLY nodes
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      IF (do_dgemm == do_dgemm_spla) THEN

         IF (opA == 'N') spla_op_A = SPLA_OP_NONE
         IF (opA == 'T') spla_op_A = SPLA_OP_TRANSPOSE

         IF (opB == 'N') spla_op_B = SPLA_OP_NONE
         IF (opB == 'T') spla_op_B = SPLA_OP_TRANSPOSE

#if __GNUC__ >= 9
         CPASSERT(IS_CONTIGUOUS(A))
         CPASSERT(IS_CONTIGUOUS(B))
         CPASSERT(IS_CONTIGUOUS(C))
#endif

         CALL offload_activate_chosen_device()
         spla_error = spla_dgemm(spla_op_A, spla_op_B, &
                                 m, n, k, alpha, &
                                 c_loc(A), lda, &
                                 c_loc(B), ldb, &
                                 beta, c_loc(C), ldc, ctx%spla_context)
         CPASSERT(spla_error == SPLA_SUCCESS)
      ELSE
#endif
         CALL dgemm(opA, opB, m, n, k, alpha, &
                    A, lda, &
                    B, ldb, beta, C, ldc)
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      END IF
#else
      MARK_USED(ctx)
#endif
      CALL timestop(handle)

   END SUBROUTINE local_gemm

! **************************************************************************************************
!> \brief create a context for handling gemm offloading
!> \param ctx newly created context
!> \param pu processing unit to run the (s,d,c,z}dgemm
! **************************************************************************************************
   SUBROUTINE local_gemm_create(ctx, pu)
      CLASS(local_gemm_ctxt_type), INTENT(out) :: ctx
      INTEGER, INTENT(in) :: pu

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER :: error_

      IF (.NOT. C_ASSOCIATED(ctx%spla_context)) THEN
         IF (do_dgemm == do_dgemm_spla) THEN
            CALL offload_activate_chosen_device()

            error_ = spla_ctx_create(ctx%spla_context, pu)
            CPASSERT(error_ == SPLA_SUCCESS)
         ELSE
            ctx%spla_context = C_NULL_PTR
         END IF
      END IF
#else
      MARK_USED(pu)
      ctx%spla_context = C_NULL_PTR
#endif
   END SUBROUTINE local_gemm_create

! **************************************************************************************************
!> \brief release resources associated to a gemm context
!> \param ctx handle
! **************************************************************************************************
   SUBROUTINE local_gemm_destroy(ctx)
      CLASS(local_gemm_ctxt_type), INTENT(inout) :: ctx

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER :: error_

      IF (do_dgemm == do_dgemm_spla) THEN
         CALL offload_activate_chosen_device()

         error_ = spla_ctx_destroy(ctx%spla_context)
         CPASSERT(error_ == SPLA_SUCCESS)
      END IF
#endif
      ctx%spla_context = C_NULL_PTR
   END SUBROUTINE local_gemm_destroy

! **************************************************************************************************
!> \brief ...
!> \param ctx ...
!> \param opThresholdGPU ...
! **************************************************************************************************
   SUBROUTINE local_gemm_set_op_threshold_gpu(ctx, opThresholdGPU)
      CLASS(local_gemm_ctxt_type), INTENT(INOUT)                                        :: ctx
      INTEGER, INTENT(in)                                :: opThresholdGPU

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER                                            :: error__

      CALL offload_activate_chosen_device()
      error__ = spla_ctx_set_op_threshold_gpu(ctx%spla_context, opThresholdGPU)
#else
      MARK_USED(ctx)
      MARK_USED(opThresholdGPU)
#endif
   END SUBROUTINE local_gemm_set_op_threshold_gpu

! **************************************************************************************************
!> \brief ...
!> \param dgemm_library ...
! **************************************************************************************************
   SUBROUTINE local_gemm_set_library(dgemm_library)
      INTEGER, INTENT(IN)                                :: dgemm_library

      do_dgemm = dgemm_library
   END SUBROUTINE local_gemm_set_library

END MODULE local_gemm_api
