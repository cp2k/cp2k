!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM kpoint_lattice_fft_unittest
   USE fft_tools,                       ONLY: finalize_fft,&
                                              init_fft
   USE kinds,                           ONLY: dp
   USE kpoint_lattice_fft,              ONLY: cell_to_k_grid_fft,&
                                              k_grid_to_cell_fft,&
                                              regular_kpoint_grid
   USE lattice_preconditioner_operator, ONLY: apply_lattice_inverse_dense,&
                                              apply_lattice_state_inverse_dense,&
                                              build_lattice_inverse,&
                                              lattice_fft_auto_select,&
                                              lattice_fft_cost_ratio,&
                                              lattice_fft_storage_ratio,&
                                              lattice_grid
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: base_block = 2, base_cell = 8, base_n = 16, &
                                                          base_ncol = 3, ncell = 31, nkp = 64, &
                                                          nround = 16, nselected = 3
   INTEGER, PARAMETER                                  :: large_block = 1000, large_ncol = 1000, &
                                                          large_num_pe = 16
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)   :: inverse_k
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: state_inverse_k
   COMPLEX(KIND=dp), DIMENSION(2, 3, nkp)              :: direct, transformed
   COMPLEX(KIND=dp), DIMENSION(2, 3, nround)           :: reconstructed, round_k
   COMPLEX(KIND=dp), DIMENSION(2, 3, nselected)         :: selected_transformed
   INTEGER                                             :: d, i, icell, info, j, jcell, ik, ix, iy, iz
   INTEGER, DIMENSION(3)                               :: base_dims, delta, large_dims
   INTEGER, DIMENSION(3, base_cell)                    :: base_index
   INTEGER, DIMENSION(3)                               :: nkp_grid
   INTEGER, DIMENSION(3, ncell)                        :: index_to_cell
   INTEGER, DIMENSION(3, nround)                       :: round_index
   INTEGER, DIMENSION(nselected)                       :: selected_kpoints
   LOGICAL                                             :: used_fft
   REAL(KIND=dp), DIMENSION(3, 3)                      :: hmat
   REAL(KIND=dp)                                       :: cost_ratio, expected_cost_ratio, imaginary_residual, &
                                                          large_block_real, large_ncell_real, &
                                                          large_ncol_real, projection_error
   REAL(KIND=dp), DIMENSION(3, base_cell)              :: base_xkp
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: base_overlap
   REAL(KIND=dp), DIMENSION(base_n, base_ncol)          :: base_input, base_input_rotated, &
                                                           base_output, base_output_rotated, &
                                                           base_state_output
   REAL(KIND=dp), DIMENSION(base_ncol, base_ncol)       :: column_rotation
   REAL(KIND=dp), DIMENSION(3, nkp)                    :: xkp, xkp_bad
   REAL(KIND=dp), DIMENSION(2, 3, ncell)               :: values_rs
   REAL(KIND=dp), DIMENSION(2, 3, nround)              :: round_rs
   TYPE(mp_comm_type)                                  :: mp_comm

   CALL mp_world_init(mp_comm)
   CALL init_fft("FFTSG", alltoall=.FALSE., pool_limit=4, &
                 wisdom_file="", plan_style=1)

   nkp_grid = [4, 4, 4]
   hmat = RESHAPE([4.0_dp, 0.1_dp, 0.2_dp, &
                   0.0_dp, 5.0_dp, 0.3_dp, &
                   0.0_dp, 0.0_dp, 6.0_dp], [3, 3])

   ik = 0
   DO iz = 0, nkp_grid(3) - 1
      DO iy = 0, nkp_grid(2) - 1
         DO ix = 0, nkp_grid(1) - 1
            ik = ik + 1
            xkp(:, ik) = [(REAL(ix, dp) + 0.5_dp)/REAL(nkp_grid(1), dp), &
                          REAL(iy, dp)/REAL(nkp_grid(2), dp), &
                          (REAL(iz, dp) + 0.25_dp)/REAL(nkp_grid(3), dp)]
         END DO
      END DO
   END DO
   IF (.NOT. regular_kpoint_grid(xkp, nkp_grid)) ERROR STOP "Regular k-point grid not recognized"

   DO icell = 1, ncell
      index_to_cell(:, icell) = [MODULO(3*icell + 1, 11) - 5, &
                                 MODULO(5*icell + 2, 9) - 4, &
                                 MODULO(7*icell + 3, 7) - 3]
      DO d = 1, 3
         values_rs(1, d, icell) = SIN(0.13_dp*REAL(icell + 2*d, dp))
         values_rs(2, d, icell) = COS(0.17_dp*REAL(2*icell - d, dp))
      END DO
   END DO

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, transformed, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Complete regular grid did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 2.0E-12_dp) THEN
      ERROR STOP "Shifted lattice FFT differs from direct transform"
   END IF

   selected_kpoints = [2, 31, 64]
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, selected_transformed, &
                           used_fft, selected_kpoints=selected_kpoints)
   IF (.NOT. used_fft) ERROR STOP "Selected distributed output did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 2.0E-12_dp) THEN
      ERROR STOP "Selected lattice-FFT output differs from the full transform"
   END IF

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, transformed, used_fft, 2, hmat)
   IF (.NOT. used_fft) ERROR STOP "Derivative did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct, 2, hmat)
   IF (MAXVAL(ABS(transformed - direct)) > 5.0E-12_dp) THEN
      ERROR STOP "Lattice-FFT derivative differs from direct transform"
   END IF

   nkp_grid = [4, 4, 1]
   ik = 0
   DO iy = 0, nkp_grid(2) - 1
      DO ix = 0, nkp_grid(1) - 1
         ik = ik + 1
         xkp(:, ik) = [(REAL(ix, dp) + 0.5_dp)/REAL(nkp_grid(1), dp), &
                       REAL(iy, dp)/REAL(nkp_grid(2), dp), 0.0_dp]
      END DO
   END DO
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, 1:ik), nkp_grid, &
                           transformed(:, :, 1:ik), used_fft)
   IF (.NOT. used_fft) ERROR STOP "Singleton grid direction did not use padded lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp(:, 1:ik), direct(:, :, 1:ik))
   IF (MAXVAL(ABS(transformed(:, :, 1:ik) - direct(:, :, 1:ik))) > 2.0E-12_dp) THEN
      ERROR STOP "Padded lattice FFT differs from direct transform"
   END IF

   icell = 0
   DO iy = 0, nkp_grid(2) - 1
      DO ix = 0, nkp_grid(1) - 1
         icell = icell + 1
         round_index(:, icell) = [ix, iy, 0]
         round_rs(:, :, icell) = values_rs(:, :, icell)
      END DO
   END DO
   CALL cell_to_k_grid_fft(round_rs, round_index, xkp(:, 1:ik), nkp_grid, round_k, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Round-trip forward transform did not use lattice FFT"
   CALL k_grid_to_cell_fft(round_k, xkp(:, 1:ik), nkp_grid, round_index, reconstructed, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Round-trip inverse transform did not use lattice FFT"
   IF (MAXVAL(ABS(reconstructed - CMPLX(round_rs, KIND=dp))) > 3.0E-12_dp) THEN
      ERROR STOP "Forward/inverse lattice-FFT round trip is inconsistent"
   END IF

   nkp_grid = [4, 4, 4]

   xkp_bad = xkp
   xkp_bad(1, nkp) = xkp_bad(1, nkp) + 0.03125_dp
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, used_fft)
   IF (used_fft) ERROR STOP "Irregular k-point list incorrectly used lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp_bad, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 1.0E-14_dp) THEN
      ERROR STOP "Direct lattice-transform fallback is inconsistent"
   END IF

   ! Verify the block-circulant inverse against direct matrix application and right covariance.
   base_dims = [4, 2, 1]
   CALL lattice_grid(base_dims, base_index, base_xkp)
   base_overlap = 0.0_dp
   DO jcell = 1, base_cell
      DO icell = 1, base_cell
         delta = MODULO(base_index(:, icell) - base_index(:, jcell), base_dims)
         DO i = 1, base_block
            IF (ALL(delta == 0)) THEN
               base_overlap((icell - 1)*base_block + i, (jcell - 1)*base_block + i) = &
                  2.8_dp + 0.2_dp*REAL(i, KIND=dp)
            ELSE IF ((delta(1) == 1 .OR. delta(1) == base_dims(1) - 1) .AND. &
                     delta(2) == 0) THEN
               base_overlap((icell - 1)*base_block + i, (jcell - 1)*base_block + i) = -0.12_dp
            ELSE IF (delta(1) == 0 .AND. delta(2) == 1) THEN
               base_overlap((icell - 1)*base_block + i, (jcell - 1)*base_block + i) = -0.08_dp
            END IF
         END DO
         IF (ALL(delta == 0)) THEN
            base_overlap((icell - 1)*base_block + 1, (jcell - 1)*base_block + 2) = 0.15_dp
            base_overlap((icell - 1)*base_block + 2, (jcell - 1)*base_block + 1) = 0.15_dp
         END IF
      END DO
   END DO
   DO j = 1, base_ncol
      DO i = 1, base_n
         base_input(i, j) = SIN(0.11_dp*REAL(2*i + 3*j, KIND=dp))
      END DO
   END DO
   CALL build_lattice_inverse(base_overlap, base_dims, inverse_k, projection_error, used_fft, info)
   IF (info /= 0) ERROR STOP "Positive lattice overlap failed block Cholesky inversion"
   IF (.NOT. used_fft) ERROR STOP "Block-circulant overlap construction did not use lattice FFT"
   IF (projection_error > 2.0E-15_dp) ERROR STOP "Exact block-circulant overlap was changed by projection"
   IF (lattice_fft_auto_select(base_dims, base_block, base_ncol, 1, 1, &
                               projection_error, used_fft, info)) THEN
      ERROR STOP "Cost-aware lattice-FFT AUTO accepted a tiny eight-cell operator"
   END IF
   IF (.NOT. lattice_fft_auto_select(base_dims, 6, 8, 1, 2, &
                                     projection_error, used_fft, info)) THEN
      ERROR STOP "Cost-aware lattice-FFT AUTO rejected an economical eight-cell operator"
   END IF
   IF (lattice_fft_auto_select([4, 1, 1], 6, 8, 1, 2, &
                               projection_error, used_fft, info)) THEN
      ERROR STOP "Conservative lattice-FFT AUTO accepted fewer than eight cells"
   END IF
   IF (lattice_fft_auto_select(base_dims, 6, 8, 1, 2, 2.0E-3_dp, used_fft, info)) THEN
      ERROR STOP "Conservative lattice-FFT AUTO accepted a large projection error"
   END IF
   IF (lattice_fft_auto_select(base_dims, 6, 8, 17, 2, &
                               projection_error, used_fft, info)) THEN
      ERROR STOP "Storage-aware lattice-FFT AUTO accepted too many state operators"
   END IF
   IF (lattice_fft_cost_ratio(base_dims, base_block, base_ncol, 1) <= 0.8_dp) THEN
      ERROR STOP "Lattice-FFT cost estimate missed the tiny-batch overhead"
   END IF
   large_dims = [100, 100, 1]
   large_ncell_real = REAL(PRODUCT(large_dims), KIND=dp)
   large_block_real = REAL(large_block, KIND=dp)
   large_ncol_real = REAL(large_ncol, KIND=dp)
   expected_cost_ratio = ( &
                         large_ncell_real*large_block_real**2*large_ncol_real + &
                         4.0_dp*large_ncell_real*large_block_real*large_ncol_real* &
                         LOG(large_ncell_real)/LOG(2.0_dp) + &
                         8.0_dp*large_ncell_real*large_block_real + &
                         2.0_dp*large_ncell_real*large_block_real*large_ncol_real* &
                         LOG(REAL(large_num_pe + 1, KIND=dp))/LOG(2.0_dp))/ &
                         ((large_ncell_real*large_block_real)**2*large_ncol_real)
   cost_ratio = lattice_fft_cost_ratio(large_dims, large_block, large_ncol, large_num_pe)
   IF (ABS(cost_ratio - expected_cost_ratio) > 1.0E-14_dp*expected_cost_ratio) THEN
      ERROR STOP "Lattice-FFT cost estimate overflowed for a large problem"
   END IF
   IF (ABS(lattice_fft_storage_ratio(base_dims, 1) - 0.25_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Lattice-FFT storage estimate is inconsistent"
   END IF
   CALL apply_lattice_inverse_dense(inverse_k, base_dims, base_input, base_output, &
                                    used_fft, imaginary_residual)
   IF (.NOT. used_fft) ERROR STOP "Block-circulant inverse application did not use lattice FFT"
   IF (imaginary_residual > 3.0E-13_dp) ERROR STOP "Lattice inverse did not reconstruct a real matrix"
   IF (MAXVAL(ABS(MATMUL(base_overlap, base_output) - base_input)) > 2.0E-11_dp) THEN
      ERROR STOP "FFT block-circulant inverse differs from the dense operator inverse"
   END IF

   ALLOCATE (state_inverse_k(base_block, base_block, base_cell, base_ncol))
   DO j = 1, base_ncol
      state_inverse_k(:, :, :, j) = (0.5_dp + 0.25_dp*REAL(j, KIND=dp))*inverse_k
   END DO
   CALL apply_lattice_state_inverse_dense(state_inverse_k, base_dims, base_input, &
                                          base_state_output, used_fft, imaginary_residual)
   IF (.NOT. used_fft) ERROR STOP "State-dependent inverse application did not use lattice FFT"
   DO j = 1, base_ncol
      IF (MAXVAL(ABS(base_state_output(:, j) - &
                     (0.5_dp + 0.25_dp*REAL(j, KIND=dp))*base_output(:, j))) > 3.0E-12_dp) THEN
         ERROR STOP "State-dependent lattice inverse used the wrong column operator"
      END IF
   END DO
   DEALLOCATE (state_inverse_k)

   column_rotation = 0.0_dp
   column_rotation(1, 1) = COS(0.31_dp)
   column_rotation(1, 2) = -SIN(0.31_dp)
   column_rotation(2, 1) = SIN(0.31_dp)
   column_rotation(2, 2) = COS(0.31_dp)
   column_rotation(3, 3) = 1.0_dp
   base_input_rotated = MATMUL(base_input, column_rotation)
   CALL apply_lattice_inverse_dense(inverse_k, base_dims, base_input_rotated, &
                                    base_output_rotated, used_fft)
   IF (MAXVAL(ABS(base_output_rotated - MATMUL(base_output, column_rotation))) > 3.0E-12_dp) THEN
      ERROR STOP "FFT lattice base is not right-covariant under orbital rotations"
   END IF
   DEALLOCATE (inverse_k)

   CALL finalize_fft(mp_comm, "")
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param values_rs ...
!> \param index_to_cell ...
!> \param xkp ...
!> \param values_k ...
!> \param deriv_direction ...
!> \param hmat ...
! **************************************************************************************************
   SUBROUTINE direct_reference(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_k
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector

      values_k = z_zero
      DO ik = 1, SIZE(xkp, 2)
         DO icell = 1, SIZE(values_rs, 3)
            factor = EXP(gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), KIND=dp)))
            IF (PRESENT(deriv_direction)) THEN
               cell_vector = MATMUL(hmat, REAL(index_to_cell(:, icell), KIND=dp))
               factor = factor*gaussi*cell_vector(deriv_direction)
            END IF
            values_k(:, :, ik) = values_k(:, :, ik) + factor*values_rs(:, :, icell)
         END DO
      END DO

   END SUBROUTINE direct_reference

END PROGRAM kpoint_lattice_fft_unittest
