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
                                              k_grid_to_cell_execute,&
                                              k_grid_to_cell_fft,&
                                              k_grid_to_cell_prepare,&
                                              k_grid_to_cell_release,&
                                              k_grid_to_cell_work_type,&
                                              regular_kpoint_grid
   USE lattice_preconditioner_operator, ONLY: &
        apply_lattice_inverse_dense, apply_lattice_state_inverse_dense, build_lattice_inverse, &
        build_lattice_local_correction, lattice_fft_auto_select, lattice_fft_cost_ratio, &
        lattice_fft_storage_ratio, lattice_grid, lattice_local_cost_ratio, &
        lattice_local_storage_ratio
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
   INTEGER, ALLOCATABLE, DIMENSION(:)                  :: correction_indices
   INTEGER, DIMENSION(3)                               :: base_dims, delta, large_dims
   INTEGER, DIMENSION(3, base_cell)                    :: base_index
   INTEGER, DIMENSION(3)                               :: nkp_grid
   INTEGER, DIMENSION(3, ncell)                        :: index_to_cell
   INTEGER, DIMENSION(3, nround)                       :: round_index
   INTEGER, DIMENSION(nselected)                       :: selected_kpoints
   LOGICAL                                             :: used_fft
   REAL(KIND=dp), DIMENSION(3, 3)                      :: hmat
   REAL(KIND=dp)                                       :: base_residual_norm, corrected_projection_error, &
                                                          corrected_residual_norm, cost_ratio, &
                                                          expected_cost_ratio, imaginary_residual, &
                                                          large_block_real, large_ncell_real, &
                                                          large_ncol_real, projection_error
   REAL(KIND=dp), DIMENSION(3, base_cell)              :: base_xkp
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: base_overlap
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: defect_operator
   REAL(KIND=dp), DIMENSION(base_n, base_ncol)          :: base_input, base_input_rotated, &
                                                           base_output, base_output_rotated, &
                                                           base_state_output, defect_residual
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coarse_inverse, defect_base_output, &
                                                          defect_corrected_output, local_operator, &
                                                          selected_output, selected_rhs
   REAL(KIND=dp), DIMENSION(base_ncol, base_ncol)       :: column_rotation
   REAL(KIND=dp), DIMENSION(3, nkp)                    :: xkp, xkp_bad
   REAL(KIND=dp), DIMENSION(2, 3, ncell)               :: values_rs
   REAL(KIND=dp), DIMENSION(2, 3, nround)              :: round_rs
   TYPE(mp_comm_type)                                  :: mp_comm
   TYPE(k_grid_to_cell_work_type)                      :: inverse_work
   REAL(KIND=dp), DIMENSION(9)                         :: subset_weights
   COMPLEX(KIND=dp), DIMENSION(2, 3, ncell)            :: weighted_reference, weighted_result

   CALL mp_world_init(mp_comm)
   CALL check_image_reconstruction(mp_comm)
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

   ! An incomplete shifted grid uses the supplied weights exactly, including
   ! negative translations. Reuse the preparation for separate matrix batches.
   subset_weights = [(REAL(ik, dp)/45.0_dp, ik=1, 9)]
   weighted_reference = z_zero
   DO icell = 1, ncell
      DO ik = 1, 9
         weighted_reference(:, :, icell) = weighted_reference(:, :, icell) + subset_weights(ik)* &
                                           EXP(-gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), dp)))*round_k(:, :, ik)
      END DO
   END DO
   CALL k_grid_to_cell_prepare(inverse_work, xkp(:, 1:9), nkp_grid, index_to_cell, &
                               weights=subset_weights, allow_incomplete=.TRUE.)
   DO i = 1, 2
      CALL k_grid_to_cell_execute(inverse_work, round_k(i:i, :, 1:9), weighted_result(i:i, :, :), used_fft)
      IF (.NOT. used_fft) ERROR STOP "Prepared incomplete grid did not use FFT"
   END DO
   IF (MAXVAL(ABS(weighted_result - weighted_reference)) > 3.0E-12_dp) THEN
      ERROR STOP "Prepared weighted inverse differs from explicit phase sum"
   END IF
   CALL k_grid_to_cell_fft(round_k(:, :, 1:9), xkp(:, 1:9), nkp_grid, index_to_cell, weighted_result, &
                           used_fft, weights=subset_weights, allow_incomplete=.TRUE.)
   IF (MAXVAL(ABS(weighted_result - weighted_reference)) > 3.0E-12_dp) THEN
      ERROR STOP "One-shot and prepared inverse transforms disagree"
   END IF
   CALL k_grid_to_cell_prepare(inverse_work, xkp(:, 1:9), nkp_grid, index_to_cell, &
                               weights=2.0_dp*subset_weights, allow_incomplete=.TRUE.)
   CALL k_grid_to_cell_execute(inverse_work, round_k(:, :, 1:9), weighted_result, used_fft)
   IF (MAXVAL(ABS(weighted_result - 2.0_dp*weighted_reference)) > 6.0E-12_dp) THEN
      ERROR STOP "Repreparing inverse transform retained stale weights"
   END IF
   CALL k_grid_to_cell_release(inverse_work)
   CALL k_grid_to_cell_release(inverse_work)

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
   IF (ABS(lattice_local_cost_ratio(base_dims, 1) - 0.28125_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Lattice local-correction cost estimate is inconsistent"
   END IF
   IF (ABS(lattice_local_storage_ratio(base_dims, 3, 1) - 0.421875_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Lattice local-correction storage estimate is inconsistent"
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

   ! Add one localized positive defect and verify the balanced cell correction.
   defect_operator = base_overlap
   defect_operator(5, 5) = defect_operator(5, 5) + 0.90_dp
   defect_operator(6, 6) = defect_operator(6, 6) + 0.60_dp
   CALL build_lattice_inverse(defect_operator, base_dims, inverse_k, projection_error, used_fft, info)
   IF (info /= 0 .OR. .NOT. used_fft) ERROR STOP "Defective lattice operator failed construction"
   IF (projection_error <= 1.0E-3_dp) ERROR STOP "Localized defect did not break lattice translations"
   IF (lattice_fft_auto_select(base_dims, 6, 8, 1, 2, projection_error, used_fft, info)) THEN
      ERROR STOP "Lattice-FFT AUTO accepted a localized defect"
   END IF
   CALL build_lattice_local_correction( &
      defect_operator, base_dims, 1, correction_indices, local_operator, coarse_inverse, &
      corrected_projection_error, info)
   IF (info /= 0) ERROR STOP "Local correction coarse inverse failed"
   IF (ANY(correction_indices /= [5, 6])) ERROR STOP "Local correction selected the wrong defect cell"
   IF (corrected_projection_error >= projection_error) THEN
      ERROR STOP "Local correction did not reduce the non-circulant residual"
   END IF

   ALLOCATE (defect_base_output(base_n, base_ncol), defect_corrected_output(base_n, base_ncol))
   CALL apply_lattice_inverse_dense(inverse_k, base_dims, base_input, defect_base_output, used_fft)
   CALL apply_balanced_local_inverse(inverse_k, base_dims, correction_indices, local_operator, &
                                     coarse_inverse, base_input, defect_corrected_output)
   defect_residual = MATMUL(defect_operator, defect_base_output) - base_input
   base_residual_norm = NORM2(defect_residual)
   defect_residual = MATMUL(defect_operator, defect_corrected_output) - base_input
   corrected_residual_norm = NORM2(defect_residual)
   IF (corrected_residual_norm >= base_residual_norm) THEN
      ERROR STOP "Balanced local correction did not improve inverse application"
   END IF
   IF (SUM(base_input*defect_corrected_output) <= 0.0_dp) THEN
      ERROR STOP "Balanced local correction lost positive definiteness"
   END IF

   ALLOCATE (selected_rhs(base_n, SIZE(correction_indices)), &
             selected_output(base_n, SIZE(correction_indices)))
   selected_rhs(:, :) = local_operator
   CALL apply_balanced_local_inverse(inverse_k, base_dims, correction_indices, local_operator, &
                                     coarse_inverse, selected_rhs, selected_output)
   selected_rhs = 0.0_dp
   DO i = 1, SIZE(correction_indices)
      selected_rhs(correction_indices(i), i) = 1.0_dp
   END DO
   IF (MAXVAL(ABS(selected_output - selected_rhs)) > 3.0E-11_dp) THEN
      ERROR STOP "Balanced local correction is not exact in its selected subspace"
   END IF
   DEALLOCATE (coarse_inverse, correction_indices, defect_base_output, defect_corrected_output, &
               local_operator, selected_output, selected_rhs)
   DEALLOCATE (inverse_k)

   CALL finalize_fft(mp_comm, "")
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check triangular image storage against rskp_transform for both neighbor orientations.
!> \param comm MPI communicator
! **************************************************************************************************
   SUBROUTINE check_image_reconstruction(comm)
      USE cp_dbcsr_api, ONLY: &
         dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
         dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_p_type, &
         dbcsr_put_block, dbcsr_release, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_symmetric
      USE cp_log_handling, ONLY: cp_add_default_logger, cp_logger_create, cp_logger_release, &
                                 cp_logger_type, cp_rm_default_logger
      USE kpoint_k_r_trafo_simple, ONLY: replicate_rs_matrices, rs_to_kp
      USE kpoint_methods, ONLY: rskp_transform
      USE kpoint_types, ONLY: kpoint_type
      USE machine, ONLY: default_output_unit
      USE message_passing, ONLY: mp_para_env_release, mp_para_env_type
      USE qs_neighbor_list_types, ONLY: add_neighbor_list, add_neighbor_node, &
                                        allocate_neighbor_list_set, deallocate_neighbor_list_set, neighbor_list_p_type
      TYPE(mp_comm_type), INTENT(IN)                     :: comm

      CHARACTER                                          :: matrix_type, mirror_type
      COMPLEX(dp)                                        :: expected(5, 5), phase, RESULT(5, 5)
      INTEGER                                            :: atom, ik, image(3), img, io, mode, &
                                                            neighbor, orientation, owner
      INTEGER, POINTER                                   :: col_dist(:), images(:, :), &
                                                            process_grid(:, :), row_dist(:), &
                                                            sizes(:)
      LOGICAL                                            :: found
      REAL(dp)                                           :: expanded(1, 5, 5, 3), point(3), sign, &
                                                            values(2, 3)
      REAL(dp), POINTER                                  :: cosine(:, :), sine(:, :)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(dbcsr_type)                                   :: cosine_matrix, sine_matrix
      TYPE(kpoint_type), POINTER                         :: points
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_p_type)                         :: neighbors

      ALLOCATE (para_env)
      CALL para_env%from_dup(comm)
      io = -1
      IF (comm%mepos == 0) io = default_output_unit
      NULLIFY (logger)
      CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io, &
                            close_global_unit_on_dealloc=.FALSE.)
      CALL cp_add_default_logger(logger)
      CALL dbcsr_init_lib(comm%get_handle(), io)
      ALLOCATE (process_grid(0:comm%num_pe - 1, 0:0), row_dist(2), col_dist(2), sizes(2))
      process_grid(:, 0) = [(img, img=0, comm%num_pe - 1)]
      col_dist = 0
      sizes = [2, 3]
      ALLOCATE (points, source(1, 3), images(3, 3))
      ALLOCATE (points%cell_to_index(-1:1, 0:0, 0:0), points%sab_nl(1))
      DO img = 1, 3
         points%cell_to_index(img - 2, 0, 0) = img
         images(:, img) = [img - 2, 0, 0]
         ALLOCATE (source(1, img)%matrix)
      END DO
      values = RESHAPE([0.3_dp, -0.2_dp, 0.7_dp, 0.1_dp, -0.4_dp, 0.6_dp], [2, 3])
      DO mode = 1, 2
         matrix_type = MERGE(dbcsr_type_symmetric, dbcsr_type_antisymmetric, mode == 1)
         mirror_type = MERGE(dbcsr_type_antisymmetric, dbcsr_type_symmetric, mode == 1)
         sign = MERGE(1.0_dp, -1.0_dp, mode == 1)
         DO orientation = 1, 2
            owner = MODULO(orientation - 1, comm%num_pe)
            row_dist = [owner, MODULO(owner + 1, comm%num_pe)]
            CALL dbcsr_distribution_new(distribution, group=comm%get_handle(), pgrid=process_grid, &
                                        row_dist=row_dist, col_dist=col_dist)
            atom = orientation
            neighbor = 3 - orientation
            image = [3 - 2*orientation, 0, 0]
            CALL allocate_neighbor_list_set(points%sab_nl(1)%neighbor_list_set, symmetric=.TRUE.)
            CALL add_neighbor_list(points%sab_nl(1)%neighbor_list_set, atom, neighbors%neighbor_list)
            CALL add_neighbor_node(neighbors%neighbor_list, neighbor, image, REAL(image, dp))
            DO img = 1, 3
               CALL dbcsr_create(source(1, img)%matrix, name="Image reconstruction", dist=distribution, &
                                 matrix_type=matrix_type, row_blk_size=sizes, col_blk_size=sizes)
               IF (comm%mepos == owner .AND. img == image(1) + 2) THEN
                  CALL dbcsr_put_block(source(1, img)%matrix, 1, 2, values)
               END IF
               CALL dbcsr_finalize(source(1, img)%matrix)
            END DO
            CALL dbcsr_create(cosine_matrix, name="Image cosine", dist=distribution, matrix_type=matrix_type, &
                              row_blk_size=sizes, col_blk_size=sizes)
            CALL dbcsr_create(sine_matrix, name="Image sine", dist=distribution, matrix_type=mirror_type, &
                              row_blk_size=sizes, col_blk_size=sizes)
            IF (comm%mepos == owner) THEN
               CALL dbcsr_put_block(cosine_matrix, 1, 2, 0.0_dp*values)
               CALL dbcsr_put_block(sine_matrix, 1, 2, 0.0_dp*values)
            END IF
            CALL dbcsr_finalize(cosine_matrix)
            CALL dbcsr_finalize(sine_matrix)
            expanded = 0.0_dp
            CALL replicate_rs_matrices(source, points, expanded, points%cell_to_index)
            DO ik = 1, 2
               point = [0.137_dp, 0.217_dp, 0.089_dp]*REAL(ik - 1, dp)
               CALL rs_to_kp(expanded(1, :, :, :), RESULT, images, point)
               phase = EXP(gaussi*twopi*point(1))
               expected = z_zero
               expected(1:2, 3:5) = phase*values
               expected(3:5, 1:2) = sign*CONJG(TRANSPOSE(expected(1:2, 3:5)))
               IF (MAXVAL(ABS(RESULT - expected)) > 1.e-12_dp) THEN
                  ERROR STOP "Image reconstruction depends on neighbor orientation"
               END IF
               CALL rskp_transform(cosine_matrix, sine_matrix, source, 1, point, &
                                   points%cell_to_index, points%sab_nl)
               CALL dbcsr_get_block_p(cosine_matrix, 1, 2, cosine, found)
               IF (comm%mepos == owner .AND. .NOT. found) ERROR STOP "Missing image cosine block"
               IF (found) THEN
                  CALL dbcsr_get_block_p(sine_matrix, 1, 2, sine, found)
                  IF (.NOT. found) ERROR STOP "Missing image sine block"
                  IF (MAXVAL(ABS(RESULT(1:2, 3:5) - CMPLX(cosine, sine, KIND=dp))) > 1.e-12_dp) THEN
                     ERROR STOP "Image reconstruction differs from rskp_transform"
                  END IF
                  cosine = 0.0_dp
                  sine = 0.0_dp
               END IF
            END DO
            CALL dbcsr_release(cosine_matrix)
            CALL dbcsr_release(sine_matrix)
            DO img = 1, 3
               CALL dbcsr_release(source(1, img)%matrix)
            END DO
            CALL dbcsr_distribution_release(distribution)
            CALL deallocate_neighbor_list_set(points%sab_nl(1)%neighbor_list_set)
         END DO
      END DO
      DO img = 1, 3
         DEALLOCATE (source(1, img)%matrix)
      END DO
      DEALLOCATE (points%cell_to_index, points%sab_nl)
      DEALLOCATE (points, source, images)
      DEALLOCATE (process_grid, row_dist, col_dist, sizes)
      CALL dbcsr_finalize_lib()
      CALL cp_rm_default_logger()
      CALL cp_logger_release(logger)
      CALL mp_para_env_release(para_env)
   END SUBROUTINE check_image_reconstruction

! **************************************************************************************************
!> \brief Apply the balanced two-level inverse used by the production local correction.
!> \param inverse_k block-circulant base inverse
!> \param lattice_dims replicated-cell dimensions
!> \param correction_indices selected coordinate basis
!> \param operator_columns columns A U
!> \param coarse_inverse inverse of U^T A U
!> \param input right-hand sides
!> \param output corrected inverse action
! **************************************************************************************************
   SUBROUTINE apply_balanced_local_inverse(inverse_k, lattice_dims, correction_indices, &
                                           operator_columns, coarse_inverse, input, output)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: inverse_k
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lattice_dims
      INTEGER, DIMENSION(:), INTENT(IN)                  :: correction_indices
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: operator_columns, coarse_inverse, input
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: output

      INTEGER                                            :: i
      LOGICAL                                            :: used_fft
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: base_input, coarse_input, &
                                                            coarse_projected, coarse_residual, &
                                                            coarse_rhs

      ALLOCATE (base_input(SIZE(input, 1), SIZE(input, 2)), &
                coarse_input(SIZE(correction_indices), SIZE(input, 2)), &
                coarse_projected(SIZE(correction_indices), SIZE(input, 2)), &
                coarse_residual(SIZE(correction_indices), SIZE(input, 2)), &
                coarse_rhs(SIZE(correction_indices), SIZE(input, 2)))
      DO i = 1, SIZE(correction_indices)
         coarse_input(i, :) = input(correction_indices(i), :)
      END DO
      coarse_rhs(:, :) = MATMUL(coarse_inverse, coarse_input)
      base_input(:, :) = input - MATMUL(operator_columns, coarse_rhs)
      CALL apply_lattice_inverse_dense(inverse_k, lattice_dims, base_input, output, used_fft)
      IF (.NOT. used_fft) ERROR STOP "Balanced correction base did not use lattice FFT"
      coarse_projected(:, :) = MATMUL(TRANSPOSE(operator_columns), output)
      coarse_residual(:, :) = MATMUL(coarse_inverse, coarse_projected)
      output(correction_indices, :) = output(correction_indices, :) + coarse_rhs - coarse_residual
      DEALLOCATE (base_input, coarse_input, coarse_projected, coarse_residual, coarse_rhs)

   END SUBROUTINE apply_balanced_local_inverse

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
      REAL(KIND=dp), DIMENSION(3)                        :: cell_index_real, cell_vector

      values_k = z_zero
      DO ik = 1, SIZE(xkp, 2)
         DO icell = 1, SIZE(values_rs, 3)
            factor = EXP(gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), KIND=dp)))
            IF (PRESENT(deriv_direction)) THEN
               cell_index_real = REAL(index_to_cell(:, icell), KIND=dp)
               cell_vector = MATMUL(hmat, cell_index_real)
               factor = factor*gaussi*cell_vector(deriv_direction)
            END IF
            values_k(:, :, ik) = values_k(:, :, ik) + factor*values_rs(:, :, icell)
         END DO
      END DO

   END SUBROUTINE direct_reference

END PROGRAM kpoint_lattice_fft_unittest
