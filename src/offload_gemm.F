!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE offload_gemm
   USE ISO_C_BINDING, ONLY: C_LOC, &
                            C_NULL_PTR, &
                            C_PTR
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
   USE input_constants, ONLY: do_dgemm_spla
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'offload_gemm'

   PUBLIC :: offload_dgemm, &
             offload_gemm_create, &
             offload_gemm_destroy, &
             offload_gemm_set_op_threshold_gpu, &
             offload_gemm_set_library

   INTEGER, PARAMETER, PUBLIC :: &
      OFFLOAD_GEMM_PU_HOST = 0, &
      OFFLOAD_GEMM_PU_GPU = 1

   INTEGER, PRIVATE :: do_dgemm = 1

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
   SUBROUTINE offload_dgemm(opA, opB, m, n, k, &
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
      TYPE(C_ptr), OPTIONAL, INTENT(inout) :: ctx

      INTEGER                                            :: handle
!     no point of using SPLA offloading on CPU ONLY nodes
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER :: spla_op_A, spla_op_B, spla_error
#endif
      CHARACTER(LEN=*), PARAMETER :: routineN = 'offload_gemm'
      CALL timeset(routineN, handle)

!     no point of using SPLA offloading on CPU ONLY nodes
#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      IF (PRESENT(ctx) .AND. do_dgemm == do_dgemm_spla) THEN

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
                                 beta, c_loc(C), ldc, ctx)
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

   END SUBROUTINE offload_dgemm

! **************************************************************************************************
!> \brief create a context for handling gemm offloading
!> \param ctx newly created context
!> \param pu processing unit to run the (s,d,c,z}dgemm
! **************************************************************************************************
   SUBROUTINE offload_gemm_create(ctx, pu)
      TYPE(c_ptr), INTENT(out) :: ctx
      INTEGER, INTENT(in) :: pu

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER :: error_

      IF (do_dgemm == do_dgemm_spla) THEN
         CALL offload_activate_chosen_device()

         error_ = spla_ctx_create(ctx, pu)
         CPASSERT(error_ == SPLA_SUCCESS)
      ELSE
         ctx = C_NULL_PTR
      END IF
#else
      MARK_USED(pu)
      MARK_USED(ctx)
      ctx = C_NULL_PTR
#endif
   END SUBROUTINE offload_gemm_create

! **************************************************************************************************
!> \brief release resources associated to a gemm context
!> \param ctx handle
! **************************************************************************************************
   SUBROUTINE offload_gemm_destroy(ctx)
      TYPE(c_ptr), INTENT(inout) :: ctx

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER :: error_

      IF (do_dgemm == do_dgemm_spla) THEN
         CALL offload_activate_chosen_device()

         error_ = spla_ctx_destroy(ctx)
         CPASSERT(error_ == SPLA_SUCCESS)
      END IF
#else
      MARK_USED(ctx)
#endif
      ctx = C_NULL_PTR
   END SUBROUTINE offload_gemm_destroy

! **************************************************************************************************
!> \brief ...
!> \param ctx ...
!> \param opThresholdGPU ...
! **************************************************************************************************
   SUBROUTINE offload_gemm_set_op_threshold_gpu(ctx, opThresholdGPU)
      TYPE(c_ptr)                                        :: ctx
      INTEGER, INTENT(in)                                :: opThresholdGPU

#if defined(__SPLA) && defined(__OFFLOAD_GEMM)
      INTEGER                                            :: error__

      CALL offload_activate_chosen_device()
      error__ = spla_ctx_set_op_threshold_gpu(ctx, opThresholdGPU)
#else
      MARK_USED(ctx)
      MARK_USED(opThresholdGPU)
#endif
   END SUBROUTINE offload_gemm_set_op_threshold_gpu

! **************************************************************************************************
!> \brief ...
!> \param dgemm_library ...
! **************************************************************************************************
   SUBROUTINE offload_gemm_set_library(dgemm_library)
      INTEGER, INTENT(IN)                                :: dgemm_library

      do_dgemm = dgemm_library
   END SUBROUTINE offload_gemm_set_library

END MODULE offload_gemm
