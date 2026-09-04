!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block-circulant inverse operators for replicated-cell OT preconditioners.
! **************************************************************************************************
MODULE lattice_preconditioner_operator
   USE fft_tools,                       ONLY: BWFFT,&
                                              FFT_RADIX_NEXT,&
                                              FWFFT,&
                                              fft_1d_many,&
                                              fft_alloc,&
                                              fft_dealloc,&
                                              fft_radix_operations
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE kpoint_lattice_fft,              ONLY: cell_to_k_grid_fft
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lattice_preconditioner_operator'

   INTEGER, PARAMETER, PUBLIC           :: lattice_fft_auto_min_cells = 8
   REAL(KIND=dp), PARAMETER, PUBLIC     :: lattice_fft_auto_max_error = 1.0E-3_dp
   REAL(KIND=dp), PARAMETER, PUBLIC     :: lattice_fft_auto_max_cost_ratio = 0.8_dp
   REAL(KIND=dp), PARAMETER, PUBLIC     :: lattice_fft_auto_max_storage_ratio = 4.0_dp

   PUBLIC :: apply_lattice_inverse_dense, &
             apply_lattice_state_inverse_dense, &
             build_lattice_inverse, &
             build_lattice_inverse_from_blocks, &
             build_lattice_local_correction, &
             lattice_fft_cost_ratio, &
             lattice_fft_auto_select, &
             lattice_fft_storage_ratio, &
             lattice_local_cost_ratio, &
             lattice_local_storage_ratio, &
             lattice_grid

CONTAINS

! **************************************************************************************************
!> \brief Enumerate a regular direct/reciprocal lattice with the first direction varying fastest.
!> \param lattice_dims number of cells in every lattice direction
!> \param index_to_cell integer direct-lattice coordinates
!> \param xkp fractional reciprocal-lattice coordinates
! **************************************************************************************************
   SUBROUTINE lattice_grid(lattice_dims, index_to_cell, xkp)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, DIMENSION(:, :), INTENT(OUT)              :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: xkp

      INTEGER                                            :: icell, ix, iy, iz

      IF (ANY(lattice_dims <= 0)) CPABORT("Invalid lattice dimensions")
      IF (SIZE(index_to_cell, 1) < 3 .OR. SIZE(xkp, 1) < 3 .OR. &
          SIZE(index_to_cell, 2) /= PRODUCT(lattice_dims) .OR. &
          SIZE(xkp, 2) /= PRODUCT(lattice_dims)) THEN
         CPABORT("Invalid lattice-grid output shape")
      END IF

      icell = 0
      DO iz = 0, lattice_dims(3) - 1
         DO iy = 0, lattice_dims(2) - 1
            DO ix = 0, lattice_dims(1) - 1
               icell = icell + 1
               index_to_cell(:, icell) = [ix, iy, iz]
               xkp(:, icell) = [REAL(ix, KIND=dp)/REAL(lattice_dims(1), KIND=dp), &
                                REAL(iy, KIND=dp)/REAL(lattice_dims(2), KIND=dp), &
                                REAL(iz, KIND=dp)/REAL(lattice_dims(3), KIND=dp)]
            END DO
         END DO
      END DO

   END SUBROUTINE lattice_grid

! **************************************************************************************************
!> \brief Project a cell-major SPD operator onto lattice translations and invert its Fourier blocks.
!> \param operator_matrix dense SPD matrix with identical cell blocks stored consecutively
!> \param lattice_dims number of replicated cells in every lattice direction
!> \param inverse_k inverse Hermitian operator block at every reciprocal-grid point
!> \param projection_error relative Frobenius error of the block-circulant projection
!> \param used_fft reports whether construction used the FFT path
!> \param info zero on success; otherwise the failing Cholesky factorization index
! **************************************************************************************************
   SUBROUTINE build_lattice_inverse(operator_matrix, lattice_dims, inverse_k, projection_error, &
                                    used_fft, info)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: operator_matrix
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(OUT)                 :: inverse_k
      REAL(KIND=dp), INTENT(OUT)                         :: projection_error
      LOGICAL, INTENT(OUT)                               :: used_fft
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: n, ncell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blocks

      IF (SIZE(operator_matrix, 1) /= SIZE(operator_matrix, 2)) THEN
         CPABORT("Lattice preconditioner operator must be square")
      END IF
      IF (ANY(lattice_dims <= 0)) CPABORT("Lattice dimensions must be positive")
      n = SIZE(operator_matrix, 1)
      ncell = PRODUCT(lattice_dims)
      IF (MODULO(n, ncell) /= 0) THEN
         CPABORT("AO dimension does not match replicated cells")
      END IF
      CALL project_lattice_operator(operator_matrix, lattice_dims, blocks, projection_error)
      CALL build_lattice_inverse_from_blocks(blocks, lattice_dims, inverse_k, used_fft, info)
      DEALLOCATE (blocks)

   END SUBROUTINE build_lattice_inverse

! **************************************************************************************************
!> \brief Transform and invert already projected real-space lattice blocks.
!> \param blocks translation-averaged real-space blocks
!> \param lattice_dims number of replicated cells in every lattice direction
!> \param inverse_k inverse Hermitian operator block at every reciprocal-grid point
!> \param used_fft reports whether construction used the FFT path
!> \param info zero on success; otherwise the failing Cholesky factorization index
! **************************************************************************************************
   SUBROUTINE build_lattice_inverse_from_blocks(blocks, lattice_dims, inverse_k, used_fft, info)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: blocks
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(OUT)                 :: inverse_k
      LOGICAL, INTENT(OUT)                               :: used_fft
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: block_size, ncell
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp

      IF (ANY(lattice_dims <= 0)) CPABORT("Lattice dimensions must be positive")
      ncell = PRODUCT(lattice_dims)
      block_size = SIZE(blocks, 1)
      IF (block_size <= 0 .OR. SIZE(blocks, 2) /= block_size .OR. SIZE(blocks, 3) /= ncell) THEN
         CPABORT("Projected lattice blocks have inconsistent dimensions")
      END IF

      ALLOCATE (index_to_cell(3, ncell), xkp(3, ncell))
      CALL lattice_grid(lattice_dims, index_to_cell, xkp)
      ALLOCATE (inverse_k(block_size, block_size, ncell))
      CALL cell_to_k_grid_fft(blocks, index_to_cell, xkp, lattice_dims, inverse_k, used_fft)
      CALL invert_hermitian_blocks(inverse_k, info)
      DEALLOCATE (index_to_cell, xkp)

   END SUBROUTINE build_lattice_inverse_from_blocks

! **************************************************************************************************
!> \brief Build a balanced local correction for selected non-circulant cell blocks.
!> \param operator_matrix dense SPD operator
!> \param lattice_dims number of replicated cells in every lattice direction
!> \param local_cell_count number of cell blocks retained in the correction
!> \param correction_indices selected global AO indices
!> \param operator_columns columns A U for the selected coordinate basis U
!> \param coarse_inverse inverse of U^T A U
!> \param corrected_projection_error residual projection error after removing selected rows/columns
!> \param info zero on success; otherwise the failing Cholesky factorization index
! **************************************************************************************************
   SUBROUTINE build_lattice_local_correction(operator_matrix, lattice_dims, local_cell_count, &
                                             correction_indices, operator_columns, coarse_inverse, &
                                             corrected_projection_error, info)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: operator_matrix
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: local_cell_count
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: correction_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: operator_columns, coarse_inverse
      REAL(KIND=dp), INTENT(OUT)                         :: corrected_projection_error
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: block_size, cell, chosen, i, index, &
                                                            local_rank, n, ncell, row_first, &
                                                            row_last
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected_cells
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: selected
      REAL(KIND=dp)                                      :: denominator
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: scores
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: projected_operator, residual
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blocks

      IF (SIZE(operator_matrix, 1) /= SIZE(operator_matrix, 2)) THEN
         CPABORT("Lattice local-correction operator must be square")
      END IF
      IF (ANY(lattice_dims <= 0)) CPABORT("Lattice dimensions must be positive")
      n = SIZE(operator_matrix, 1)
      ncell = PRODUCT(lattice_dims)
      IF (MODULO(n, ncell) /= 0) CPABORT("AO dimension does not match replicated cells")
      IF (local_cell_count <= 0 .OR. local_cell_count >= ncell) THEN
         CPABORT("Lattice local correction requires between one and ncell-1 cells")
      END IF
      block_size = n/ncell
      local_rank = local_cell_count*block_size

      CALL project_lattice_operator(operator_matrix, lattice_dims, blocks, &
                                    corrected_projection_error, projected_operator)
      DEALLOCATE (blocks)
      ALLOCATE (residual(n, n))
      residual(:, :) = operator_matrix - projected_operator
      DEALLOCATE (projected_operator)
      ALLOCATE (scores(ncell), selected(ncell), selected_cells(local_cell_count))
      selected = .FALSE.
      DO chosen = 1, local_cell_count
         scores = -1.0_dp
         DO cell = 1, ncell
            IF (selected(cell)) CYCLE
            row_first = (cell - 1)*block_size + 1
            row_last = cell*block_size
            scores(cell) = SUM(residual(row_first:row_last, :)**2) + &
                           SUM(residual(:, row_first:row_last)**2) - &
                           SUM(residual(row_first:row_last, row_first:row_last)**2)
         END DO
         index = MAXLOC(scores, DIM=1, MASK=.NOT. selected)
         selected(index) = .TRUE.
         selected_cells(chosen) = index
         row_first = (index - 1)*block_size + 1
         row_last = index*block_size
         residual(row_first:row_last, :) = 0.0_dp
         residual(:, row_first:row_last) = 0.0_dp
      END DO

      ALLOCATE (correction_indices(local_rank))
      index = 0
      DO chosen = 1, local_cell_count
         cell = selected_cells(chosen)
         DO i = 1, block_size
            index = index + 1
            correction_indices(index) = (cell - 1)*block_size + i
         END DO
      END DO

      denominator = SUM(operator_matrix**2)
      corrected_projection_error = SQRT(SUM(residual**2)/MAX(denominator, TINY(denominator)))

      ALLOCATE (operator_columns(n, local_rank), coarse_inverse(local_rank, local_rank))
      operator_columns(:, :) = operator_matrix(:, correction_indices)
      coarse_inverse(:, :) = operator_matrix(correction_indices, correction_indices)
      CALL dpotrf('U', local_rank, coarse_inverse, local_rank, info)
      IF (info == 0) CALL dpotri('U', local_rank, coarse_inverse, local_rank, info)
      IF (info == 0) THEN
         DO i = 1, local_rank - 1
            coarse_inverse(i + 1:local_rank, i) = coarse_inverse(i, i + 1:local_rank)
         END DO
      END IF

      DEALLOCATE (residual, scores, selected, selected_cells)

   END SUBROUTINE build_lattice_local_correction

! **************************************************************************************************
!> \brief Project a dense cell-major operator onto its block-circulant translation average.
!> \param operator_matrix dense square operator
!> \param lattice_dims number of replicated cells in every lattice direction
!> \param blocks translation-averaged real-space blocks
!> \param projection_error relative Frobenius projection error
!> \param projected_operator optional dense block-circulant projection
! **************************************************************************************************
   SUBROUTINE project_lattice_operator(operator_matrix, lattice_dims, blocks, projection_error, &
                                       projected_operator)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: operator_matrix
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: blocks
      REAL(KIND=dp), INTENT(OUT)                         :: projection_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: projected_operator

      INTEGER                                            :: alpha, beta, block_size, delta_index, &
                                                            icell, jcell, n, ncell
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      INTEGER, DIMENSION(3)                              :: delta
      REAL(KIND=dp)                                      :: denominator, difference, numerator
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp

      IF (SIZE(operator_matrix, 1) /= SIZE(operator_matrix, 2)) THEN
         CPABORT("Lattice preconditioner operator must be square")
      END IF
      IF (ANY(lattice_dims <= 0)) CPABORT("Lattice dimensions must be positive")
      n = SIZE(operator_matrix, 1)
      ncell = PRODUCT(lattice_dims)
      IF (MODULO(n, ncell) /= 0) CPABORT("AO dimension does not match replicated cells")
      block_size = n/ncell

      ALLOCATE (index_to_cell(3, ncell), xkp(3, ncell))
      CALL lattice_grid(lattice_dims, index_to_cell, xkp)
      ALLOCATE (blocks(block_size, block_size, ncell), source=0.0_dp)
      DO jcell = 1, ncell
         DO icell = 1, ncell
            ! With X(k)=sum_R exp(+ikR) X(R), the convolution block is indexed by R-T.
            delta = MODULO(index_to_cell(:, icell) - index_to_cell(:, jcell), lattice_dims)
            delta_index = cell_linear_index(delta, lattice_dims)
            DO beta = 1, block_size
               DO alpha = 1, block_size
                  blocks(alpha, beta, delta_index) = blocks(alpha, beta, delta_index) + &
                                                     operator_matrix((icell - 1)*block_size + alpha, &
                                                                     (jcell - 1)*block_size + beta)/ &
                                                     REAL(ncell, KIND=dp)
               END DO
            END DO
         END DO
      END DO

      IF (PRESENT(projected_operator)) ALLOCATE (projected_operator(n, n))
      numerator = 0.0_dp
      denominator = 0.0_dp
      DO jcell = 1, ncell
         DO icell = 1, ncell
            delta = MODULO(index_to_cell(:, icell) - index_to_cell(:, jcell), lattice_dims)
            delta_index = cell_linear_index(delta, lattice_dims)
            DO beta = 1, block_size
               DO alpha = 1, block_size
                  IF (PRESENT(projected_operator)) THEN
                     projected_operator((icell - 1)*block_size + alpha, &
                                        (jcell - 1)*block_size + beta) = &
                        blocks(alpha, beta, delta_index)
                  END IF
                  difference = operator_matrix((icell - 1)*block_size + alpha, &
                                               (jcell - 1)*block_size + beta) - &
                               blocks(alpha, beta, delta_index)
                  numerator = numerator + difference*difference
                  denominator = denominator + &
                                operator_matrix((icell - 1)*block_size + alpha, &
                                                (jcell - 1)*block_size + beta)**2
               END DO
            END DO
         END DO
      END DO
      projection_error = SQRT(numerator/MAX(denominator, TINY(denominator)))

      DEALLOCATE (index_to_cell, xkp)

   END SUBROUTINE project_lattice_operator

! **************************************************************************************************
!> \brief Estimate lattice-FFT application work relative to a dense inverse application.
!> \param lattice_dims number of explicitly replicated cells
!> \param block_size number of orbitals in the reference cell
!> \param rhs_count number of orbital columns transformed together
!> \param num_pe number of message-passing ranks sharing the operator
!> \return estimated lattice-to-dense work ratio
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION lattice_fft_cost_ratio(lattice_dims, block_size, rhs_count, num_pe) &
      RESULT(ratio)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: block_size, rhs_count, num_pe

      INTEGER                                            :: ncell, nrhs
      REAL(KIND=dp)                                      :: block_size_real, communication_work, &
                                                            dense_work, fft_work, matrix_work, &
                                                            ncell_real, nrhs_real, setup_work

      ncell = MAX(1, PRODUCT(lattice_dims))
      nrhs = MAX(1, rhs_count)
      IF (block_size <= 0) THEN
         ratio = HUGE(ratio)
         RETURN
      END IF
      ncell_real = REAL(ncell, KIND=dp)
      nrhs_real = REAL(nrhs, KIND=dp)
      block_size_real = REAL(block_size, KIND=dp)

      ! This deliberately simple model counts real dense multiply-adds as the reference.
      ! Four scalar operations per FFT butterfly account for the two complex transforms;
      ! the fixed packing term prevents AUTO from selecting tiny batches.  The logarithmic
      ! communication term represents the input/output all-to-all routing.
      dense_work = (ncell_real*block_size_real)**2*nrhs_real
      matrix_work = ncell_real*block_size_real**2*nrhs_real
      fft_work = 4.0_dp*ncell_real*block_size_real*nrhs_real* &
                 LOG(REAL(MAX(2, ncell), KIND=dp))/LOG(2.0_dp)
      setup_work = 8.0_dp*ncell_real*block_size_real
      communication_work = 2.0_dp*ncell_real*block_size_real*nrhs_real* &
                           LOG(REAL(MAX(2, num_pe + 1), KIND=dp))/LOG(2.0_dp)
      ratio = (matrix_work + fft_work + setup_work + communication_work)/MAX(dense_work, TINY(dense_work))

   END FUNCTION lattice_fft_cost_ratio

! **************************************************************************************************
!> \brief Estimate reciprocal-block storage relative to one dense real operator.
!> \param lattice_dims number of explicitly replicated cells
!> \param operator_count number of distinct state-dependent reciprocal operators
!> \return estimated lattice-to-dense storage ratio
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION lattice_fft_storage_ratio(lattice_dims, operator_count) RESULT(ratio)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: operator_count

      INTEGER                                            :: ncell

      ncell = MAX(1, PRODUCT(lattice_dims))
      ! Reciprocal blocks are complex; the block-size dependence cancels against
      ! the reference dense matrix.  A state-independent operator has count one.
      ratio = 2.0_dp*REAL(MAX(1, operator_count), KIND=dp)/REAL(ncell, KIND=dp)

   END FUNCTION lattice_fft_storage_ratio

! **************************************************************************************************
!> \brief Estimate balanced-local-correction work relative to a dense inverse application.
!> \param lattice_dims number of replicated cells
!> \param local_cell_count number of cell blocks retained in the correction
!> \return estimated correction-to-dense work ratio
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION lattice_local_cost_ratio(lattice_dims, local_cell_count) &
      RESULT(ratio)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: local_cell_count

      INTEGER                                            :: ncell
      REAL(KIND=dp)                                      :: local_fraction

      ncell = MAX(1, PRODUCT(lattice_dims))
      IF (local_cell_count <= 0) THEN
         ratio = 0.0_dp
         RETURN
      END IF
      local_fraction = REAL(local_cell_count, KIND=dp)/REAL(ncell, KIND=dp)
      ! Both corrected and dense applications scale linearly with the right-hand-side count.
      ratio = 2.0_dp*(local_fraction + local_fraction**2)

   END FUNCTION lattice_local_cost_ratio

! **************************************************************************************************
!> \brief Estimate balanced-local-correction storage relative to one dense real operator.
!> \param lattice_dims number of replicated cells
!> \param operator_count number of distinct state-dependent local corrections
!> \param local_cell_count number of cell blocks retained in each correction
!> \return estimated correction-to-dense storage ratio
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION lattice_local_storage_ratio(lattice_dims, operator_count, &
                                                           local_cell_count) RESULT(ratio)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: operator_count, local_cell_count

      INTEGER                                            :: ncell
      REAL(KIND=dp)                                      :: local_fraction

      ncell = MAX(1, PRODUCT(lattice_dims))
      IF (local_cell_count <= 0) THEN
         ratio = 0.0_dp
         RETURN
      END IF
      local_fraction = REAL(local_cell_count, KIND=dp)/REAL(ncell, KIND=dp)
      ratio = REAL(MAX(1, operator_count), KIND=dp)*(local_fraction + local_fraction**2)

   END FUNCTION lattice_local_storage_ratio

! **************************************************************************************************
!> \brief Conservative error-, work-, and storage-aware AUTO selector.
!> \param lattice_dims number of explicitly replicated cells
!> \param block_size number of orbitals in the reference cell
!> \param rhs_count number of orbital columns transformed together
!> \param operator_count number of distinct reciprocal operators stored
!> \param num_pe number of message-passing ranks sharing the operator
!> \param projection_error relative Frobenius projection error
!> \param used_fft whether construction used the FFT path
!> \param info zero when all projected Fourier blocks were positive definite
!> \return true when AUTO may use the lattice inverse
! **************************************************************************************************
   PURE LOGICAL FUNCTION lattice_fft_auto_select(lattice_dims, block_size, rhs_count, operator_count, &
                                                 num_pe, projection_error, used_fft, info) &
      RESULT(selected)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: block_size, rhs_count, operator_count, &
                                                            num_pe
      REAL(KIND=dp), INTENT(IN)                          :: projection_error
      LOGICAL, INTENT(IN)                                :: used_fft
      INTEGER, INTENT(IN)                                :: info

      selected = PRODUCT(lattice_dims) >= lattice_fft_auto_min_cells .AND. &
                 projection_error <= lattice_fft_auto_max_error .AND. used_fft .AND. info == 0
      selected = selected .AND. &
                 lattice_fft_cost_ratio(lattice_dims, block_size, rhs_count, num_pe) <= &
                 lattice_fft_auto_max_cost_ratio .AND. &
                 lattice_fft_storage_ratio(lattice_dims, operator_count) <= &
                 lattice_fft_auto_max_storage_ratio

   END FUNCTION lattice_fft_auto_select

! **************************************************************************************************
!> \brief Apply inverse Fourier blocks to a dense cell-major real matrix.
!> \param inverse_k inverse overlap blocks on a complete reciprocal grid
!> \param lattice_dims number of replicated cells in every lattice direction
!> \param matrix_in matrix to precondition
!> \param matrix_out preconditioned matrix
!> \param used_fft reports whether both transforms used the FFT path
!> \param imaginary_residual largest discarded imaginary component
! **************************************************************************************************
   SUBROUTINE apply_lattice_inverse_dense(inverse_k, lattice_dims, matrix_in, matrix_out, &
                                          used_fft, imaginary_residual)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: inverse_k
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix_in
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: matrix_out
      LOGICAL, INTENT(OUT)                               :: used_fft
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: imaginary_residual

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: input_k, input_rs, output_k, output_rs
      INTEGER                                            :: block_size, icell, ncell
      LOGICAL                                            :: forward_fft, inverse_fft

      IF (ANY(lattice_dims <= 0)) CPABORT("Lattice dimensions must be positive")
      ncell = PRODUCT(lattice_dims)
      block_size = SIZE(inverse_k, 1)
      IF (SIZE(inverse_k, 2) /= block_size .OR. SIZE(inverse_k, 3) /= ncell) THEN
         CPABORT("Inconsistent inverse blocks in lattice preconditioner")
      END IF
      IF (SIZE(matrix_in, 1) /= block_size*ncell .OR. &
          ANY(SHAPE(matrix_out) /= SHAPE(matrix_in))) THEN
         CPABORT("Inconsistent dense matrix shape in lattice preconditioner")
      END IF

      ALLOCATE (input_rs(block_size, SIZE(matrix_in, 2), ncell))
      DO icell = 1, ncell
         input_rs(:, :, icell) = CMPLX( &
                                 matrix_in((icell - 1)*block_size + 1:icell*block_size, :), 0.0_dp, KIND=dp)
      END DO

      ALLOCATE (input_k(block_size, SIZE(matrix_in, 2), ncell))
      ALLOCATE (output_k(block_size, SIZE(matrix_in, 2), ncell))
      ALLOCATE (output_rs(block_size, SIZE(matrix_in, 2), ncell))
      CALL batched_lattice_fft(input_rs, lattice_dims, BWFFT, input_k, forward_fft)
      DO icell = 1, ncell
         output_k(:, :, icell) = MATMUL(inverse_k(:, :, icell), input_k(:, :, icell))
      END DO
      CALL batched_lattice_fft(output_k, lattice_dims, FWFFT, output_rs, inverse_fft)
      used_fft = forward_fft .AND. inverse_fft
      IF (PRESENT(imaginary_residual)) imaginary_residual = MAXVAL(ABS(AIMAG(output_rs)))

      DO icell = 1, ncell
         matrix_out((icell - 1)*block_size + 1:icell*block_size, :) = REAL(output_rs(:, :, icell), KIND=dp)
      END DO

      DEALLOCATE (input_k, input_rs, output_k, output_rs)

   END SUBROUTINE apply_lattice_inverse_dense

! **************************************************************************************************
!> \brief Apply one inverse Fourier operator to every column of a dense real matrix.
!> \param inverse_k state-dependent inverse blocks on a complete reciprocal grid
!> \param lattice_dims number of replicated cells in every lattice direction
!> \param matrix_in matrix to precondition, with one state per column
!> \param matrix_out preconditioned matrix
!> \param used_fft reports whether both transforms used the FFT path
!> \param imaginary_residual largest discarded imaginary component
!> \param state_indices optional global operator index for every input column
! **************************************************************************************************
   SUBROUTINE apply_lattice_state_inverse_dense(inverse_k, lattice_dims, matrix_in, matrix_out, &
                                                used_fft, imaginary_residual, state_indices)

      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(IN)                                      :: inverse_k
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix_in
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: matrix_out
      LOGICAL, INTENT(OUT)                               :: used_fft
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: imaginary_residual
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: state_indices

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: input_k, input_rs, output_k, output_rs
      INTEGER                                            :: block_size, icell, j, ncell
      LOGICAL                                            :: forward_fft, inverse_fft

      IF (ANY(lattice_dims <= 0)) CPABORT("Lattice dimensions must be positive")
      ncell = PRODUCT(lattice_dims)
      block_size = SIZE(inverse_k, 1)
      IF (SIZE(inverse_k, 2) /= block_size .OR. SIZE(inverse_k, 3) /= ncell) THEN
         CPABORT("Inconsistent state-dependent inverse blocks in lattice preconditioner")
      END IF
      IF (PRESENT(state_indices)) THEN
         IF (SIZE(state_indices) /= SIZE(matrix_in, 2) .OR. &
             ANY(state_indices < 1) .OR. ANY(state_indices > SIZE(inverse_k, 4))) THEN
            CPABORT("Invalid state indices in lattice preconditioner")
         END IF
      ELSE IF (SIZE(inverse_k, 4) /= SIZE(matrix_in, 2)) THEN
         CPABORT("Inconsistent state-dependent inverse blocks in lattice preconditioner")
      END IF
      IF (SIZE(matrix_in, 1) /= block_size*ncell .OR. &
          ANY(SHAPE(matrix_out) /= SHAPE(matrix_in))) THEN
         CPABORT("Inconsistent dense state matrix shape in lattice preconditioner")
      END IF

      ALLOCATE (input_rs(block_size, SIZE(matrix_in, 2), ncell))
      DO icell = 1, ncell
         input_rs(:, :, icell) = CMPLX( &
                                 matrix_in((icell - 1)*block_size + 1:icell*block_size, :), 0.0_dp, KIND=dp)
      END DO

      ALLOCATE (input_k(block_size, SIZE(matrix_in, 2), ncell))
      ALLOCATE (output_k(block_size, SIZE(matrix_in, 2), ncell))
      ALLOCATE (output_rs(block_size, SIZE(matrix_in, 2), ncell))
      CALL batched_lattice_fft(input_rs, lattice_dims, BWFFT, input_k, forward_fft)
      DO icell = 1, ncell
         DO j = 1, SIZE(matrix_in, 2)
            IF (PRESENT(state_indices)) THEN
               output_k(:, j, icell) = &
                  MATMUL(inverse_k(:, :, icell, state_indices(j)), input_k(:, j, icell))
            ELSE
               output_k(:, j, icell) = MATMUL(inverse_k(:, :, icell, j), input_k(:, j, icell))
            END IF
         END DO
      END DO
      CALL batched_lattice_fft(output_k, lattice_dims, FWFFT, output_rs, inverse_fft)
      used_fft = forward_fft .AND. inverse_fft
      IF (PRESENT(imaginary_residual)) imaginary_residual = MAXVAL(ABS(AIMAG(output_rs)))

      DO icell = 1, ncell
         matrix_out((icell - 1)*block_size + 1:icell*block_size, :) = REAL(output_rs(:, :, icell), KIND=dp)
      END DO

      DEALLOCATE (input_k, input_rs, output_k, output_rs)

   END SUBROUTINE apply_lattice_state_inverse_dense

! **************************************************************************************************
!> \brief Transform a complete replicated-cell batch with three batched one-dimensional FFTs.
!> \param values_in matrices on the complete direct or reciprocal lattice
!> \param lattice_dims logical lattice dimensions
!> \param fsign transform direction, BWFFT for cell-to-k and FWFFT for k-to-cell
!> \param values_out transformed matrices on the complete lattice
!> \param used_fft reports whether the batched FFT path succeeded
! **************************************************************************************************
   SUBROUTINE batched_lattice_fft(values_in, lattice_dims, fsign, values_out, used_fft)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: values_in
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, INTENT(IN)                                :: fsign
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_out
      LOGICAL, INTENT(OUT)                               :: used_fft

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: grid
      INTEGER                                            :: attempt, batch_index, d, i, icell, ix, &
                                                            iy, iz, j, radix_length, stat
      INTEGER, DIMENSION(3)                              :: input_stride, nfft, output_stride, stride
      LOGICAL                                            :: compatible
      REAL(KIND=dp)                                      :: normalization

      IF (ANY(lattice_dims <= 0) .OR. SIZE(values_in, 3) /= PRODUCT(lattice_dims) .OR. &
          ANY(SHAPE(values_out) /= SHAPE(values_in))) THEN
         CPABORT("Invalid shape in batched lattice FFT")
      END IF
      IF (fsign /= FWFFT .AND. fsign /= BWFFT) THEN
         CPABORT("Invalid direction in batched lattice FFT")
      END IF
      nfft = lattice_dims
      compatible = .TRUE.
      DO d = 1, 3
         compatible = .FALSE.
         DO attempt = 0, 15
            IF (nfft(d) >= 3) THEN
               CALL fft_radix_operations(nfft(d), radix_length, FFT_RADIX_NEXT)
               IF (radix_length == nfft(d)) THEN
                  compatible = .TRUE.
                  EXIT
               END IF
            END IF
            nfft(d) = nfft(d) + lattice_dims(d)
         END DO
         IF (.NOT. compatible) EXIT
      END DO
      IF (.NOT. compatible) THEN
         values_out = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         used_fft = .FALSE.
         RETURN
      END IF

      stride = nfft/lattice_dims
      input_stride = 1
      output_stride = 1
      IF (fsign == BWFFT) input_stride = stride
      IF (fsign == FWFFT) output_stride = stride
      ALLOCATE (grid(SIZE(values_in, 1)*SIZE(values_in, 2), nfft(1), nfft(2), nfft(3)))
      grid = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      icell = 0
      DO iz = 0, lattice_dims(3) - 1
         DO iy = 0, lattice_dims(2) - 1
            DO ix = 0, lattice_dims(1) - 1
               icell = icell + 1
               DO j = 1, SIZE(values_in, 2)
                  DO i = 1, SIZE(values_in, 1)
                     batch_index = i + SIZE(values_in, 1)*(j - 1)
                     grid(batch_index, ix*input_stride(1) + 1, iy*input_stride(2) + 1, &
                          iz*input_stride(3) + 1) = values_in(i, j, icell)
                  END DO
               END DO
            END DO
         END DO
      END DO

      stat = 0
      DO d = 1, 3
         CALL transform_grid_dimension(grid, d, fsign, stat)
         IF (stat /= 0) EXIT
      END DO
      IF (stat /= 0) THEN
         values_out = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         used_fft = .FALSE.
         DEALLOCATE (grid)
         RETURN
      END IF

      normalization = 1.0_dp
      IF (fsign == FWFFT) THEN
         normalization = REAL(PRODUCT(nfft), KIND=dp)/REAL(PRODUCT(lattice_dims), KIND=dp)
      END IF
      icell = 0
      DO iz = 0, lattice_dims(3) - 1
         DO iy = 0, lattice_dims(2) - 1
            DO ix = 0, lattice_dims(1) - 1
               icell = icell + 1
               DO j = 1, SIZE(values_out, 2)
                  DO i = 1, SIZE(values_out, 1)
                     batch_index = i + SIZE(values_out, 1)*(j - 1)
                     values_out(i, j, icell) = normalization* &
                                               grid(batch_index, ix*output_stride(1) + 1, iy*output_stride(2) + 1, &
                                                    iz*output_stride(3) + 1)
                  END DO
               END DO
            END DO
         END DO
      END DO
      used_fft = .TRUE.
      DEALLOCATE (grid)

   END SUBROUTINE batched_lattice_fft

! **************************************************************************************************
!> \brief Transform every line of one dimension of a four-dimensional FFT batch.
!> \param grid batched three-dimensional grid
!> \param direction grid dimension to transform
!> \param fsign transform direction
!> \param stat zero on success
! **************************************************************************************************
   SUBROUTINE transform_grid_dimension(grid, direction, fsign, stat)

      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: grid
      INTEGER, INTENT(IN)                                :: direction, fsign
      INTEGER, INTENT(OUT)                               :: stat

      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: line_in, line_out
      INTEGER                                            :: ibatch, iline, ix, iy, iz, m, n
      REAL(KIND=dp)                                      :: scale

      n = SIZE(grid, direction + 1)
      m = SIZE(grid)/n
      CALL fft_alloc(line_in, [n, m])
      CALL fft_alloc(line_out, [n, m])
      iline = 0
      SELECT CASE (direction)
      CASE (1)
         DO iz = 1, SIZE(grid, 4)
            DO iy = 1, SIZE(grid, 3)
               DO ibatch = 1, SIZE(grid, 1)
                  iline = iline + 1
                  line_in(:, iline) = grid(ibatch, :, iy, iz)
               END DO
            END DO
         END DO
      CASE (2)
         DO iz = 1, SIZE(grid, 4)
            DO ix = 1, SIZE(grid, 2)
               DO ibatch = 1, SIZE(grid, 1)
                  iline = iline + 1
                  line_in(:, iline) = grid(ibatch, ix, :, iz)
               END DO
            END DO
         END DO
      CASE (3)
         DO iy = 1, SIZE(grid, 3)
            DO ix = 1, SIZE(grid, 2)
               DO ibatch = 1, SIZE(grid, 1)
                  iline = iline + 1
                  line_in(:, iline) = grid(ibatch, ix, iy, :)
               END DO
            END DO
         END DO
      CASE DEFAULT
         CPABORT("Invalid batched FFT grid direction")
      END SELECT

      scale = 1.0_dp
      IF (fsign == FWFFT) scale = 1.0_dp/REAL(n, KIND=dp)
      CALL fft_1d_many(fsign, n, m, .FALSE., .FALSE., n, n, &
                       line_in, line_out, scale, stat)
      IF (stat == 0) THEN
         iline = 0
         SELECT CASE (direction)
         CASE (1)
            DO iz = 1, SIZE(grid, 4)
               DO iy = 1, SIZE(grid, 3)
                  DO ibatch = 1, SIZE(grid, 1)
                     iline = iline + 1
                     grid(ibatch, :, iy, iz) = line_out(:, iline)
                  END DO
               END DO
            END DO
         CASE (2)
            DO iz = 1, SIZE(grid, 4)
               DO ix = 1, SIZE(grid, 2)
                  DO ibatch = 1, SIZE(grid, 1)
                     iline = iline + 1
                     grid(ibatch, ix, :, iz) = line_out(:, iline)
                  END DO
               END DO
            END DO
         CASE (3)
            DO iy = 1, SIZE(grid, 3)
               DO ix = 1, SIZE(grid, 2)
                  DO ibatch = 1, SIZE(grid, 1)
                     iline = iline + 1
                     grid(ibatch, ix, iy, :) = line_out(:, iline)
                  END DO
               END DO
            END DO
         END SELECT
      END IF
      CALL fft_dealloc(line_in)
      CALL fft_dealloc(line_out)

   END SUBROUTINE transform_grid_dimension

! **************************************************************************************************
!> \brief Convert zero-based lattice coordinates to the cell-major linear index.
!> \param coordinate zero-based lattice coordinate
!> \param lattice_dims number of cells in every direction
!> \return one-based linear index
! **************************************************************************************************
   PURE INTEGER FUNCTION cell_linear_index(coordinate, lattice_dims) RESULT(index)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: coordinate, lattice_dims

      index = coordinate(1) + lattice_dims(1)*(coordinate(2) + &
                                               lattice_dims(2)*coordinate(3)) + 1

   END FUNCTION cell_linear_index

! **************************************************************************************************
!> \brief Hermitize, Cholesky factorize, and invert each small complex block.
!> \param blocks matrix blocks, replaced by their inverse
!> \param info zero on success; otherwise the one-based failing block index
! **************************************************************************************************
   SUBROUTINE invert_hermitian_blocks(blocks, info)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: blocks
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(KIND=dp)                                   :: value
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: solution, work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: factor, inverse
      INTEGER                                            :: i, iblock, j, k, m, n
      REAL(KIND=dp)                                      :: diagonal

      n = SIZE(blocks, 1)
      CPASSERT(SIZE(blocks, 2) == n)
      ALLOCATE (factor(n, n), inverse(n, n), solution(n), work(n))
      info = 0
      DO iblock = 1, SIZE(blocks, 3)
         factor = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         inverse(:, :) = 0.5_dp*(blocks(:, :, iblock) + &
                                 CONJG(TRANSPOSE(blocks(:, :, iblock))))

         ! Small-block upper Cholesky factorization, inverse=U^H U.
         DO j = 1, n
            diagonal = REAL(inverse(j, j), KIND=dp)
            DO m = 1, j - 1
               diagonal = diagonal - ABS(factor(m, j))**2
            END DO
            IF (.NOT. ieee_is_finite(diagonal) .OR. diagonal <= 0.0_dp) THEN
               info = iblock
               DEALLOCATE (factor, inverse, solution, work)
               RETURN
            END IF
            factor(j, j) = CMPLX(SQRT(diagonal), 0.0_dp, KIND=dp)
            DO k = j + 1, n
               value = inverse(j, k)
               DO m = 1, j - 1
                  value = value - CONJG(factor(m, j))*factor(m, k)
               END DO
               factor(j, k) = value/factor(j, j)
            END DO
         END DO

         ! Solve U^H U x_j=e_j for every inverse column. A small explicit
         ! triangular solve keeps the complete small-block kernel self-contained.
         DO j = 1, n
            work = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            solution = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO i = 1, n
               value = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               IF (i == j) value = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
               DO m = 1, i - 1
                  value = value - CONJG(factor(m, i))*work(m)
               END DO
               work(i) = value/CONJG(factor(i, i))
            END DO
            DO i = n, 1, -1
               value = work(i)
               DO m = i + 1, n
                  value = value - factor(i, m)*solution(m)
               END DO
               solution(i) = value/factor(i, i)
            END DO
            inverse(:, j) = solution
         END DO
         blocks(:, :, iblock) = 0.5_dp*(inverse + CONJG(TRANSPOSE(inverse)))
      END DO
      DEALLOCATE (factor, inverse, solution, work)

   END SUBROUTINE invert_hermitian_blocks

END MODULE lattice_preconditioner_operator
