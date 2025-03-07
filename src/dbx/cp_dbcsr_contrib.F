!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dbcsr_contrib
   USE OMP_LIB,                         ONLY: omp_get_num_threads
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_clear, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_data_size, dbcsr_get_info, dbcsr_get_num_blocks, &
        dbcsr_get_occupation, dbcsr_get_readonly_block_p, dbcsr_get_stored_coordinates, &
        dbcsr_has_symmetry, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_readonly_start, dbcsr_iterator_start, dbcsr_iterator_stop, &
        dbcsr_iterator_type, dbcsr_put_block, dbcsr_reserve_blocks, dbcsr_type
   USE dbm_tests,                       ONLY: generate_larnv_seed
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: dbcsr_hadamard_product
   PUBLIC :: dbcsr_maxabs
   PUBLIC :: dbcsr_frobenius_norm
   PUBLIC :: dbcsr_gershgorin_norm
   PUBLIC :: dbcsr_init_random
   PUBLIC :: dbcsr_reserve_diag_blocks
   PUBLIC :: dbcsr_reserve_all_blocks
   PUBLIC :: dbcsr_add_on_diag
   PUBLIC :: dbcsr_dot
   PUBLIC :: dbcsr_trace
   PUBLIC :: dbcsr_get_block_diag
   PUBLIC :: dbcsr_scale_by_vector
   PUBLIC :: dbcsr_get_diag
   PUBLIC :: dbcsr_set_diag
   PUBLIC :: dbcsr_checksum
   PUBLIC :: dbcsr_print

CONTAINS

! **************************************************************************************************
!> \brief Hadamard product: C = A . B (C needs to be different from A and B)
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
! **************************************************************************************************
   SUBROUTINE dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_c

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hadamard_product'

      INTEGER                                            :: col, handle, nblkrows_tot_a, &
                                                            nblkrows_tot_b, nblkrows_tot_c, row
      LOGICAL                                            :: found_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_a, block_b
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_get_info(matrix_a, nblkrows_total=nblkrows_tot_a)
      CALL dbcsr_get_info(matrix_b, nblkrows_total=nblkrows_tot_b)
      CALL dbcsr_get_info(matrix_c, nblkrows_total=nblkrows_tot_c)
      IF (nblkrows_tot_a /= nblkrows_tot_b .OR. nblkrows_tot_a /= nblkrows_tot_c) THEN
         CPABORT("matrices not consistent")
      END IF

      CALL dbcsr_clear(matrix_c)
      CALL dbcsr_iterator_readonly_start(iter, matrix_a)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block_a)
         CALL dbcsr_get_readonly_block_p(matrix_b, row, col, block_b, found_b)
         IF (found_b) THEN
            CALL dbcsr_put_block(matrix_c, row, col, block_a*block_b)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(matrix_c)
      CALL timestop(handle)
   END SUBROUTINE dbcsr_hadamard_product

! **************************************************************************************************
!> \brief Compute the maxabs norm of a dbcsr matrix
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_maxabs(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(dp)                                           :: norm

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_maxabs'

      INTEGER                                            :: handle
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      norm = 0.0_dp
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, block=block)
         norm = MAX(norm, MAXVAL(ABS(block)))
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(matrix, group=group)
      CALL group%max(norm)

      CALL timestop(handle)
   END FUNCTION dbcsr_maxabs

! **************************************************************************************************
!> \brief Compute the frobenius norm of a dbcsr matrix
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_frobenius_norm(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(dp)                                           :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_frobenius_norm'

      INTEGER                                            :: col, handle, row
      LOGICAL                                            :: has_symmetry
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      has_symmetry = dbcsr_has_symmetry(matrix)
      norm = 0.0_dp
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block)
         IF (has_symmetry .AND. row /= col) THEN
            norm = norm + 2.0_dp*SUM(block**2)
         ELSE
            norm = norm + SUM(block**2)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(matrix, group=group)
      CALL group%sum(norm)
      norm = SQRT(norm)

      CALL timestop(handle)
   END FUNCTION dbcsr_frobenius_norm

! **************************************************************************************************
!> \brief Compute the gershgorin norm of a dbcsr matrix
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_gershgorin_norm(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(dp)                                           :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_gershgorin_norm'

      INTEGER                                            :: col, col_offset, handle, i, j, ncol, &
                                                            nrow, row, row_offset
      LOGICAL                                            :: has_symmetry
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      has_symmetry = dbcsr_has_symmetry(matrix)
      CALL dbcsr_get_info(matrix, nfullrows_total=nrow, nfullcols_total=ncol)
      CPASSERT(nrow == ncol)
      ALLOCATE (buffer(nrow))
      buffer = 0.0_dp

      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_offset=row_offset, col_offset=col_offset)
         DO j = 1, SIZE(block, 2)
            DO i = 1, SIZE(block, 1)
               buffer(row_offset + i - 1) = buffer(row_offset + i - 1) + ABS(block(i, j))
               IF (has_symmetry .AND. row /= col) THEN
                  buffer(col_offset + j - 1) = buffer(col_offset + j - 1) + ABS(block(i, j))
               END IF
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(matrix, group=group)
      CALL group%sum(buffer)
      norm = MAXVAL(buffer)
      DEALLOCATE (buffer)

      CALL timestop(handle)
   END FUNCTION dbcsr_gershgorin_norm

! **************************************************************************************************
!> \brief Fills the given matrix with random numbers.
!> \param matrix ...
!> \param keep_sparsity ...
! **************************************************************************************************
   SUBROUTINE dbcsr_init_random(matrix, keep_sparsity)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, OPTIONAL                                  :: keep_sparsity

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_init_random'

      INTEGER                                            :: col, col_size, handle, ncol, nrow, row, &
                                                            row_size
      INTEGER, DIMENSION(4)                              :: iseed
      LOGICAL                                            :: my_keep_sparsity
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      my_keep_sparsity = .FALSE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity
      IF (.NOT. my_keep_sparsity) CALL dbcsr_reserve_all_blocks(matrix)
      CALL dbcsr_get_info(matrix, nblkrows_total=nrow, nblkcols_total=ncol)

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_size=row_size, col_size=col_size)
         ! set the seed for dlarnv, is here to guarantee same value of the random numbers
         ! for all layouts (and block distributions)
         iseed = generate_larnv_seed(row, nrow, col, ncol, 1)
         CALL dlarnv(1, iseed, row_size*col_size, block(1, 1))
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_init_random

! **************************************************************************************************
!> \brief Reserves all diagonal blocks.
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_diag_blocks(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_diag_blocks'

      INTEGER                                            :: handle, i, k, mynode, nblkcols_total, &
                                                            nblkrows_total, owner
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_diag
      INTEGER, DIMENSION(:), POINTER                     :: local_rows
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows_total, nblkcols_total=nblkcols_total)
      CPASSERT(nblkrows_total == nblkcols_total)

      CALL dbcsr_get_info(matrix, local_rows=local_rows, distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)
      ALLOCATE (local_diag(SIZE(local_rows)))

      k = 0
      DO i = 1, SIZE(local_rows)
         CALL dbcsr_get_stored_coordinates(matrix, row=local_rows(i), column=local_rows(i), processor=owner)
         IF (owner == mynode) THEN
            k = k + 1
            local_diag(k) = local_rows(i)
         END IF
      END DO

      CALL dbcsr_reserve_blocks(matrix, rows=local_diag(1:k), cols=local_diag(1:k))
      DEALLOCATE (local_diag)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_reserve_diag_blocks

! **************************************************************************************************
!> \brief Reserves all blocks.
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_all_blocks(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_all_blocks'

      INTEGER                                            :: handle, i, j, k, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cols, rows
      INTEGER, DIMENSION(:), POINTER                     :: local_cols, local_rows

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_get_info(matrix, local_rows=local_rows, local_cols=local_cols)
      n = SIZE(local_rows)*SIZE(local_cols)
      ALLOCATE (rows(n), cols(n))

      k = 0
      DO i = 1, SIZE(local_rows)
      DO j = 1, SIZE(local_cols)
         k = k + 1
         rows(k) = local_rows(i)
         cols(k) = local_cols(j)
      END DO
      END DO

      CALL dbcsr_reserve_blocks(matrix, rows=rows(1:k), cols=cols(1:k))
      DEALLOCATE (rows, cols)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_reserve_all_blocks

! **************************************************************************************************
!> \brief Adds the given scalar to the diagonal of the matrix. Reserves any missing diagonal blocks.
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add_on_diag(matrix, alpha)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(kind=dp), INTENT(IN)                          :: alpha

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_add_on_diag'

      INTEGER                                            :: col, col_size, handle, i, row, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_reserve_diag_blocks(matrix)

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_size=row_size, col_size=col_size)
         IF (row == col) THEN
            CPASSERT(row_size == col_size)
            DO i = 1, row_size
               block(i, i) = block(i, i) + alpha
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_add_on_diag

! **************************************************************************************************
!> \brief Computes the dot product of two matrices, also known as the trace of their matrix product.
!> \param matrix_a ...
!> \param matrix_b ...
!> \param trace ...
! **************************************************************************************************
   SUBROUTINE dbcsr_dot(matrix_a, matrix_b, trace)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(OUT)                         :: trace

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_dot'

      INTEGER                                            :: col, handle, row
      LOGICAL                                            :: found_b, has_symmetry
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_a, block_b
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CPASSERT(dbcsr_has_symmetry(matrix_a) .EQV. dbcsr_has_symmetry(matrix_b))
      has_symmetry = dbcsr_has_symmetry(matrix_a)

      trace = 0.0_dp
      CALL dbcsr_iterator_readonly_start(iter, matrix_a)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block_a)
         IF (SIZE(block_a) == 0) CYCLE ! Skip zero-sized blocks.
         CALL dbcsr_get_readonly_block_p(matrix_b, row, col, block_b, found_b)
         IF (found_b) THEN
            IF (has_symmetry .AND. row /= col) THEN
               trace = trace + 2.0_dp*SUM(block_a*block_b)
            ELSE
               trace = trace + SUM(block_a*block_b)
            END IF
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(matrix_a, group=group)
      CALL group%sum(trace)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_dot

! **************************************************************************************************
!> \brief Computes the trace of the given matrix, also known as the sum of its diagonal elements.
!> \param matrix ...
!> \param trace ...
! **************************************************************************************************
   SUBROUTINE dbcsr_trace(matrix, trace)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp), INTENT(OUT)                         :: trace

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_trace'

      INTEGER                                            :: col, col_size, handle, i, row, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      trace = 0.0_dp
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_size=row_size, col_size=col_size)
         IF (row == col) THEN
            CPASSERT(row_size == col_size)
            DO i = 1, row_size
               trace = trace + block(i, i)
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(matrix, group=group)
      CALL group%sum(trace)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_trace

! **************************************************************************************************
!> \brief Copies the diagonal blocks of matrix into diag.
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_block_diag(matrix, diag)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: diag

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_diag'

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: col, handle, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_get_info(matrix, name=name)
      CALL dbcsr_create(diag, template=matrix, name='diag of '//name)

      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block)
         IF (row == col) THEN
            CALL dbcsr_put_block(diag, row, col, block)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(diag)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_get_block_diag

! **************************************************************************************************
!> \brief Scales the rows/columns of given matrix.
!> \param matrix Matrix to be scaled in-place.
!> \param alpha Vector with scaling factors.
!> \param side Side from which to apply the vector. Allowed values are 'right' and 'left'.
! **************************************************************************************************
   SUBROUTINE dbcsr_scale_by_vector(matrix, alpha, side)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: alpha
      CHARACTER(LEN=*), INTENT(IN)                       :: side

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_by_vector'

      INTEGER                                            :: col_offset, col_size, handle, i, &
                                                            nfullcols_total, nfullrows_total, &
                                                            row_offset, row_size
      LOGICAL                                            :: right
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      IF (side == 'right') THEN
         right = .TRUE.
      ELSE IF (side == 'left') THEN
         right = .FALSE.
      ELSE
         CPABORT("Unknown side: "//TRIM(side))
      END IF

      ! Check that alpha and matrix have matching sizes.
      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      IF (right) THEN
         CPASSERT(nfullcols_total == SIZE(alpha))
      ELSE
         CPASSERT(nfullrows_total == SIZE(alpha))
      END IF

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, block=block, row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         IF (SIZE(block) == 0) CYCLE ! Skip zero-sized blocks.
         IF (right) THEN
            DO i = 1, col_size
               block(:, i) = block(:, i)*alpha(col_offset + i - 1)
            END DO
         ELSE
            DO i = 1, row_size
               block(i, :) = block(i, :)*alpha(row_offset + i - 1)
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_scale_by_vector

! **************************************************************************************************
!> \brief Copies the diagonal elements from the given matrix into the given array.
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_diag(matrix, diag)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: diag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_get_diag'

      INTEGER                                            :: col, col_size, handle, i, &
                                                            nfullcols_total, nfullrows_total, row, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CPASSERT(nfullrows_total == nfullcols_total)
      CPASSERT(nfullrows_total == SIZE(diag))

      diag(:) = 0.0_dp
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_size=row_size, &
                                        col_size=col_size, row_offset=row_offset)
         IF (row == col) THEN
            CPASSERT(row_size == col_size)
            DO i = 1, row_size
               diag(row_offset + i - 1) = block(i, i)
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_get_diag

! **************************************************************************************************
!> \brief Copies the diagonal elements from the given array into the given matrix.
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_set_diag(matrix, diag)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: diag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_set_diag'

      INTEGER                                            :: col, col_size, handle, i, &
                                                            nfullcols_total, nfullrows_total, row, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CPASSERT(nfullrows_total == nfullcols_total)
      CPASSERT(nfullrows_total == SIZE(diag))

      CALL dbcsr_reserve_diag_blocks(matrix)

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_size=row_size, &
                                        col_size=col_size, row_offset=row_offset)
         IF (row == col) THEN
            CPASSERT(row_size == col_size)
            DO i = 1, row_size
               block(i, i) = diag(row_offset + i - 1)
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_set_diag

! **************************************************************************************************
!> \brief Calculates the checksum of a DBCSR matrix.
!> \param matrix ...
!> \param pos Enable position-dependent checksum.
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_checksum(matrix, pos) RESULT(checksum)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: pos
      REAL(KIND=dp)                                      :: checksum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_checksum'

      INTEGER                                            :: col_offset, col_size, handle, i, j, &
                                                            row_offset, row_size
      LOGICAL                                            :: my_pos
      REAL(KIND=dp)                                      :: position_factor
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      my_pos = .FALSE.
      IF (PRESENT(pos)) THEN
         my_pos = pos
      END IF

      checksum = 0.0_dp
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, block=block, row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         IF (my_pos) THEN
            DO i = 1, row_size
            DO j = 1, col_size
               position_factor = LOG(REAL((row_offset + i - 1)*(col_offset + j - 1), KIND=dp))
               checksum = checksum + block(i, j)*position_factor
            END DO
            END DO
         ELSE
            checksum = checksum + SUM(block**2)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(matrix, group=group)
      CALL group%sum(checksum)

      CALL timestop(handle)
   END FUNCTION dbcsr_checksum

! **************************************************************************************************
!> \brief Prints given matrix in matlab format (only present blocks).
!> \param matrix ...
!> \param variable_name ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE dbcsr_print(matrix, variable_name, unit_nr)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: variable_name
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbcsr_print'

      CHARACTER(len=default_string_length)               :: my_variable_name, name
      INTEGER :: col_offset, col_size, handle, i, iw, j, nblkcols_total, nblkrows_total, &
         nfullcols_total, nfullrows_total, row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)
      CPASSERT(omp_get_num_threads() == 1)

      iw = default_output_unit
      IF (PRESENT(unit_nr)) iw = unit_nr

      my_variable_name = 'a'
      IF (PRESENT(variable_name)) my_variable_name = variable_name

      ! Print matrix properties.
      CALL dbcsr_get_info(matrix, name=name, &
                          nblkrows_total=nblkrows_total, nblkcols_total=nblkcols_total, &
                          nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      WRITE (iw, *) "===", routineN, "==="
      WRITE (iw, *) "Name:", name
      WRITE (iw, *) "Symmetry:", dbcsr_has_symmetry(matrix)
      WRITE (iw, *) "Number of blocks:", dbcsr_get_num_blocks(matrix)
      WRITE (iw, *) "Data size:", dbcsr_get_data_size(matrix)
      WRITE (iw, *) "Occupation:", dbcsr_get_occupation(matrix)
      WRITE (iw, *) "Full size:", nfullrows_total, "x", nfullcols_total
      WRITE (iw, *) "Blocked size:", nblkrows_total, "x", nblkcols_total

      ! Print matrix blocks.
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, block=block, row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO i = 1, row_size
         DO j = 1, col_size
            WRITE (iw, '(A,I4,A,I4,A,E23.16,A)') TRIM(my_variable_name)//'(', &
               row_offset + i - 1, ',', col_offset + j - 1, ')=', block(i, j), ';'
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_print

END MODULE cp_dbcsr_contrib
