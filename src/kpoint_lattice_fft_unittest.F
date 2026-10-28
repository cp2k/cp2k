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
                                              k_grid_to_cell_prepare,&
                                              k_grid_to_cell_release,&
                                              k_grid_to_cell_work_type
   USE lattice_preconditioner_operator, ONLY: apply_lattice_inverse_dense,&
                                              apply_lattice_state_inverse_dense,&
                                              build_lattice_inverse,&
                                              build_lattice_local_correction,&
                                              lattice_fft_cost_ratio,&
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
                                                          nselected = 3
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)   :: inverse_k
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: state_inverse_k
   COMPLEX(KIND=dp), DIMENSION(2, 3, nkp)              :: direct, transformed
   COMPLEX(KIND=dp), DIMENSION(2, 3, nselected)         :: selected_transformed
   INTEGER                                             :: d, i, icell, info, j, jcell, ik, ix, iy, iz
   INTEGER, ALLOCATABLE, DIMENSION(:)                  :: correction_indices
   INTEGER, DIMENSION(3)                               :: base_dims, delta
   INTEGER, DIMENSION(3, base_cell)                    :: base_index
   INTEGER, DIMENSION(3)                               :: nkp_grid
   INTEGER, DIMENSION(3, ncell)                        :: index_to_cell
   INTEGER, DIMENSION(nselected)                       :: selected_kpoints
   LOGICAL                                             :: used_fft
   REAL(KIND=dp), DIMENSION(3, 3)                      :: hmat
   REAL(KIND=dp)                                       :: corrected_projection_error, cost_ratio, &
                                                          imaginary_residual, projection_error
   REAL(KIND=dp), DIMENSION(3, base_cell)              :: base_xkp
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: base_overlap
   REAL(KIND=dp), DIMENSION(base_n, base_n)            :: defect_operator
   REAL(KIND=dp), DIMENSION(base_n, base_ncol)          :: base_input, base_output, base_state_output
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coarse_inverse, local_operator
   REAL(KIND=dp), DIMENSION(base_block, base_block)    :: coarse_residual
   REAL(KIND=dp), DIMENSION(3, nkp)                    :: xkp, xkp_bad
   REAL(KIND=dp), DIMENSION(2, 3, ncell)               :: values_rs
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

   DO icell = 1, ncell
      index_to_cell(:, icell) = [MODULO(3*icell + 1, 11) - 5, &
                                 MODULO(5*icell + 2, 9) - 4, &
                                 MODULO(2*icell + 3, 7) - 3]
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
                           used_fft, 2, hmat, selected_kpoints=selected_kpoints)
   IF (.NOT. used_fft) ERROR STOP "Selected derivative did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct, 2, hmat)
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 5.0E-12_dp) THEN
      ERROR STOP "Selected lattice-FFT derivative differs from direct transform"
   END IF

   xkp_bad = xkp
   xkp_bad(1, nkp) = xkp_bad(1, nkp) + 0.03125_dp
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, used_fft)
   IF (used_fft) ERROR STOP "Irregular k-point list incorrectly used lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp_bad, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 1.0E-14_dp) THEN
      ERROR STOP "Direct lattice-transform fallback is inconsistent"
   END IF

   ! An incomplete shifted grid must use the supplied, non-normalized weights.
   ! Include negative translations and reuse preparation for separate matrix batches.
   subset_weights = [(REAL(ik, dp)/23.0_dp, ik=1, 9)]
   transformed(:, :, 1:9) = CMPLX(values_rs(:, :, 1:9), values_rs(:, :, 2:10), KIND=dp)
   weighted_reference = z_zero
   DO icell = 1, ncell
      DO ik = 1, 9
         weighted_reference(:, :, icell) = weighted_reference(:, :, icell) + subset_weights(ik)* &
                                           EXP(-gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), dp)))* &
                                           transformed(:, :, ik)
      END DO
   END DO
   CALL k_grid_to_cell_prepare(inverse_work, xkp(:, 1:9), nkp_grid, index_to_cell, &
                               weights=subset_weights, allow_incomplete=.TRUE.)
   DO i = 1, 2
      CALL k_grid_to_cell_execute(inverse_work, transformed(i:i, :, 1:9), weighted_result(i:i, :, :), used_fft)
      IF (.NOT. used_fft) ERROR STOP "Prepared incomplete grid did not use FFT"
   END DO
   IF (MAXVAL(ABS(weighted_result - weighted_reference)) > 3.0E-12_dp) THEN
      ERROR STOP "Prepared weighted inverse differs from explicit phase sum"
   END IF
   CALL k_grid_to_cell_release(inverse_work)

   ! Verify the block-circulant inverse by applying the original dense operator.
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
   ! Check the large-problem scale without duplicating the heuristic's formula.
   cost_ratio = lattice_fft_cost_ratio([100, 100, 1], 1000, 1000, 16)
   IF (.NOT. (cost_ratio > 1.0E-4_dp .AND. cost_ratio < 2.0E-4_dp)) THEN
      ERROR STOP "Lattice-FFT cost estimate overflowed for a large problem"
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

   ! Check the production builder's selected columns and coarse inverse directly.
   defect_operator = base_overlap
   defect_operator(5, 5) = defect_operator(5, 5) + 0.90_dp
   defect_operator(6, 6) = defect_operator(6, 6) + 0.60_dp
   CALL build_lattice_inverse(defect_operator, base_dims, inverse_k, projection_error, used_fft, info)
   IF (info /= 0 .OR. .NOT. used_fft) ERROR STOP "Defective lattice operator failed construction"
   IF (projection_error <= 1.0E-3_dp) ERROR STOP "Localized defect did not break lattice translations"
   CALL build_lattice_local_correction( &
      defect_operator, base_dims, 1, correction_indices, local_operator, coarse_inverse, &
      corrected_projection_error, info)
   IF (info /= 0) ERROR STOP "Local correction coarse inverse failed"
   IF (SIZE(correction_indices) /= base_block) ERROR STOP "Local correction has the wrong rank"
   IF (ANY(correction_indices /= [5, 6])) ERROR STOP "Local correction selected the wrong defect cell"
   IF (corrected_projection_error >= projection_error) THEN
      ERROR STOP "Local correction did not reduce the non-circulant residual"
   END IF

   IF (MAXVAL(ABS(local_operator - defect_operator(:, correction_indices))) > 1.0E-14_dp) THEN
      ERROR STOP "Local correction returned incorrect operator columns"
   END IF
   coarse_residual = MATMUL(defect_operator(correction_indices, correction_indices), coarse_inverse)
   DO i = 1, base_block
      coarse_residual(i, i) = coarse_residual(i, i) - 1.0_dp
   END DO
   IF (MAXVAL(ABS(coarse_residual)) > 3.0E-13_dp) ERROR STOP "Local correction coarse inverse is inconsistent"
   DEALLOCATE (coarse_inverse, correction_indices, local_operator)
   DEALLOCATE (inverse_k)

   CALL finalize_fft(mp_comm, "")
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check triangular image storage and its signed mirror for both neighbor orientations.
!> \param comm MPI communicator
! **************************************************************************************************
   SUBROUTINE check_image_reconstruction(comm)
      USE cp_dbcsr_api, ONLY: &
         dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
         dbcsr_finalize, dbcsr_finalize_lib, dbcsr_init_lib, dbcsr_p_type, dbcsr_put_block, &
         dbcsr_release, dbcsr_type_antisymmetric, dbcsr_type_symmetric
      USE cp_log_handling, ONLY: cp_add_default_logger, cp_logger_create, cp_logger_release, &
                                 cp_logger_type, cp_rm_default_logger
      USE kpoint_k_r_trafo_simple, ONLY: replicate_rs_matrices
      USE kpoint_types, ONLY: kpoint_type
      USE machine, ONLY: default_output_unit
      USE message_passing, ONLY: mp_para_env_release, mp_para_env_type
      USE qs_neighbor_list_types, ONLY: add_neighbor_list, add_neighbor_node, &
                                        allocate_neighbor_list_set, deallocate_neighbor_list_set, neighbor_list_p_type
      TYPE(mp_comm_type), INTENT(IN)                     :: comm

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: atom, image(3), img, io, mode, neighbor, &
                                                            orientation, owner
      INTEGER, POINTER                                   :: col_dist(:), process_grid(:, :), &
                                                            row_dist(:), sizes(:)
      REAL(dp)                                           :: expanded(1, 5, 5, 3), expected(5, 5, 3), &
                                                            sign, values(2, 3)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
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
      ALLOCATE (points, source(1, 3))
      ALLOCATE (points%cell_to_index(-1:1, 0:0, 0:0), points%sab_nl(1))
      DO img = 1, 3
         points%cell_to_index(img - 2, 0, 0) = img
         ALLOCATE (source(1, img)%matrix)
      END DO
      values = RESHAPE([0.3_dp, -0.2_dp, 0.7_dp, 0.1_dp, -0.4_dp, 0.6_dp], [2, 3])
      DO mode = 1, 2
         matrix_type = MERGE(dbcsr_type_symmetric, dbcsr_type_antisymmetric, mode == 1)
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
            expanded = 0.0_dp
            CALL replicate_rs_matrices(source, points, expanded, points%cell_to_index)
            expected = 0.0_dp
            expected(1:2, 3:5, 3) = values
            expected(3:5, 1:2, 1) = sign*TRANSPOSE(values)
            IF (MAXVAL(ABS(expanded(1, :, :, :) - expected)) > 1.e-12_dp) THEN
               ERROR STOP "Image reconstruction depends on neighbor orientation"
            END IF
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
      DEALLOCATE (points, source)
      DEALLOCATE (process_grid, row_dist, col_dist, sizes)
      CALL dbcsr_finalize_lib()
      CALL cp_rm_default_logger()
      CALL cp_logger_release(logger)
      CALL mp_para_env_release(para_env)
   END SUBROUTINE check_image_reconstruction

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
