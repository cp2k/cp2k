!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common DBCSR matrix operations used by GW modules.
!> \par History
!>      09.2026 created Jan Wilhelm
! **************************************************************************************************
MODULE gw_utils_dbcsr
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_release, dbcsr_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: dbcsr_contract_ABA, hadamard_product, hadamard_product_inplace

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils_dbcsr'

CONTAINS

! **************************************************************************************************
!> \brief Computes the scaled element-wise product C = factor (A ◦ B) while preserving the
!>        block structure of A. Blocks absent from B are retained in C with zero values.
!> \param matrix_A First factor and source of the block structure
!> \param matrix_B Second factor
!> \param matrix_C Scaled element-wise product
!> \param factor Scaling factor
! **************************************************************************************************
   SUBROUTINE hadamard_product(matrix_A, matrix_B, matrix_C, factor)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: factor

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hadamard_product'

      INTEGER                                            :: handle, icol, irow
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_B, block_C
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_C, matrix_A)
      CALL dbcsr_iterator_start(iterator, matrix_C)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, irow, icol, block_C)
         CALL dbcsr_get_block_p(matrix_B, irow, icol, block_B, found)
         IF (found) THEN
            block_C(:, :) = factor*block_C(:, :)*block_B(:, :)
         ELSE
            block_C(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)

      CALL timestop(handle)

   END SUBROUTINE hadamard_product

! **************************************************************************************************
!> \brief Form A = factor * (A element-wise B) without changing A's block structure.
!> \param matrix_A First factor, overwritten by the product.
!> \param matrix_B Second factor; a missing block represents zero.
!> \param factor Product scale factor.
! **************************************************************************************************
   SUBROUTINE hadamard_product_inplace(matrix_A, matrix_B, factor)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B
      REAL(KIND=dp), INTENT(IN)                          :: factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hadamard_product_inplace'

      INTEGER                                            :: handle, icol, irow
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_A, block_B
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iterator, matrix_A)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, irow, icol, block_A)
         CALL dbcsr_get_block_p(matrix_B, irow, icol, block_B, found)
         IF (found) THEN
            block_A(:, :) = factor*block_A(:, :)*block_B(:, :)
         ELSE
            block_A(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)

      CALL timestop(handle)
   END SUBROUTINE hadamard_product_inplace

! **************************************************************************************************
!> \brief Computes C=A B A^T or C=A^T B A for DBCSR matrices.
!> \param trans_A_left transposition applied to the left occurrence of A
!> \param trans_A_right transposition applied to the right occurrence of A
!> \param matrix_A left and right matrix A
!> \param matrix_B input matrix B
!> \param matrix_C output matrix C
!> \param eps_filter filtering threshold for both matrix multiplications
!> \param retain_sparsity if true, only existing blocks of C are filled
! **************************************************************************************************
   SUBROUTINE dbcsr_contract_ABA(trans_A_left, trans_A_right, matrix_A, matrix_B, matrix_C, &
                                 eps_filter, retain_sparsity)
      CHARACTER(LEN=1), INTENT(IN)                       :: trans_A_left, trans_A_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_contract_ABA'

      INTEGER                                            :: handle
      LOGICAL                                            :: my_retain_sparsity
      TYPE(dbcsr_type)                                   :: work

      CALL timeset(routineN, handle)

      my_retain_sparsity = .FALSE.
      IF (PRESENT(retain_sparsity)) my_retain_sparsity = retain_sparsity

      CALL dbcsr_create(work, template=matrix_A)

      IF (trans_A_left == "N" .AND. trans_A_right == "T") THEN
         ! C = A B A^T
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_A, matrix_B, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", 1.0_dp, work, matrix_A, &
                             0.0_dp, matrix_C, filter_eps=eps_filter, &
                             retain_sparsity=my_retain_sparsity)
      ELSE IF (trans_A_left == "T" .AND. trans_A_right == "N") THEN
         ! C = A^T B A
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_B, matrix_A, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_A, work, &
                             0.0_dp, matrix_C, filter_eps=eps_filter, &
                             retain_sparsity=my_retain_sparsity)
      ELSE
         CPABORT("Unsupported transposition pair in dbcsr_contract_ABA")
      END IF
      CALL dbcsr_release(work)

      CALL timestop(handle)

   END SUBROUTINE dbcsr_contract_ABA

END MODULE gw_utils_dbcsr
