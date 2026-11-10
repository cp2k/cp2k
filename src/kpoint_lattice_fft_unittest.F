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
   USE kpoint_transform,                ONLY: cell_to_k_transform,&
                                              k_to_real_cell_transform,&
                                              rskp_backend_fft,&
                                              rskp_backend_gemm
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
   INTEGER, PARAMETER                                  :: backends(2) = [rskp_backend_gemm, rskp_backend_fft]
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)   :: inverse_k
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: state_inverse_k
   COMPLEX(KIND=dp), DIMENSION(2, 3, nkp)              :: direct, transformed
   COMPLEX(KIND=dp), DIMENSION(2, 3, nselected)         :: selected_transformed
   INTEGER                                             :: backend, d, i, icell, info, j, jcell, ik, ix, iy, iz
   INTEGER, ALLOCATABLE, DIMENSION(:)                  :: correction_indices
   INTEGER, DIMENSION(3)                               :: base_dims, delta
   INTEGER, DIMENSION(3, base_cell)                    :: base_index
   INTEGER, DIMENSION(3)                               :: nkp_grid
   INTEGER, DIMENSION(3, ncell)                        :: index_to_cell
   INTEGER, DIMENSION(3, 259)                          :: panel_cells
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
   REAL(KIND=dp), DIMENSION(2, 3, ncell)               :: real_transformed
   REAL(KIND=dp), ALLOCATABLE                          :: panel_input(:, :, :), panel_inverse(:, :, :)
   COMPLEX(KIND=dp), ALLOCATABLE                       :: panel_output(:, :, :), panel_reference(:, :, :)

   CALL mp_world_init(mp_comm)
   CALL check_image_reconstruction(mp_comm)
   CALL init_fft("FFTW3", alltoall=.FALSE., pool_limit=4, &
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

   BLOCK
      USE cp_fm_struct, ONLY: cp_fm_struct_type
      USE cp_fm_types, ONLY: cp_fm_p_type, cp_fm_type
      USE kinds, ONLY: int_8
      USE kpoint_transform, ONLY: cell_to_k_backend_evaluate_batch, cell_to_k_backend_prepare, &
                                  cell_to_k_sparse_fft, rskp_route_type, rskp_values_type
      USE kpoint_lattice_fft, ONLY: cell_to_k_grid_release
      USE message_passing, ONLY: mp_para_env_type, mp_para_env_release
      INTEGER :: p, sign
      COMPLEX(dp), TARGET :: dense(4, 3, nkp), expected(4, 3, nkp)
      REAL(dp), TARGET :: split(4, 3, nkp, 2)
      REAL(dp) :: scratch(0)
      TYPE(cp_fm_struct_type), TARGET :: structure
      TYPE(cp_fm_type), TARGET :: matrices(nkp, 2)
      TYPE(cp_fm_p_type) :: destinations(nkp, 2)
      TYPE(rskp_route_type) :: route
      TYPE(rskp_values_type) :: snapshot
      TYPE(mp_para_env_type), POINTER :: env

      ALLOCATE (env)
      CALL env%from_dup(mp_comm)
      structure%ref_count = 1
      structure%local_leading_dimension = 4
      route%matrix_struct => structure
      route%world_target = .TRUE.
      route%ncol = 3
      route%transform%backend = rskp_backend_fft
      ALLOCATE (route%entry(8, ncell), source=0)
      route%entry(5:7, :) = index_to_cell
      ALLOCATE (route%block_order(ncell), source=[(p, p=1, ncell)])
      ALLOCATE (route%block_first(2), source=[1, ncell + 1])
      ALLOCATE (route%element_first(2), source=[1, 4])
      ALLOCATE (route%element_offset(2, 3))
      route%element_offset(:, 1) = [1_int_8, 0_int_8]
      route%element_offset(:, 2) = [5_int_8, 2_int_8]
      route%element_offset(:, 3) = [10_int_8, 7_int_8]
      ALLOCATE (route%transform%packed_first(2), source=[1, 3*ncell + 1])
      ALLOCATE (snapshot%packed(3*ncell), source=RESHAPE(values_rs(1, :, :), [3*ncell]))
      CALL cell_to_k_backend_prepare(route, xkp, nkp_grid, 32_int_8*1024_int_8**2)
      DO p = 1, nkp
         DO sign = 1, 2
            matrices(p, sign)%matrix_struct => structure
            matrices(p, sign)%local_data => split(:, :, p, sign)
            destinations(p, sign)%matrix => matrices(p, sign)
         END DO
      END DO
      DO sign = -1, 1, 2
         snapshot%mirror_sign = sign
         snapshot%scale = CMPLX(0.3_dp, 0.7_dp, dp)
         ALLOCATE (snapshot%fft_values(3, nkp))
         CALL cell_to_k_sparse_fft(snapshot, route, xkp, nkp_grid, RESHAPE([1, nkp], [2, 1]), env, &
                                   1_int_8, cache=snapshot%fft_values)
         IF (MAXVAL(ABS(snapshot%fft_values - transformed(1, :, :))) > 2.0E-12_dp) THEN
            ERROR STOP 'Canonical FFT cache changed raw Fourier values'
         END IF
         expected = z_zero
         DO p = 1, nkp
            expected(1, 1, p) = snapshot%scale*transformed(1, 1, p)
            expected(1, 2, p) = snapshot%scale*transformed(1, 2, p)
            expected(2, 3, p) = snapshot%scale*transformed(1, 3, p)
            expected(2, 1, p) = snapshot%scale*REAL(sign, dp)*CONJG(transformed(1, 2, p))
            expected(3, 2, p) = snapshot%scale*REAL(sign, dp)*CONJG(transformed(1, 3, p))
         END DO
         CALL cell_to_k_backend_evaluate_batch(snapshot, route, xkp, 1, scratch, batch=dense)
         IF (MAXVAL(ABS(dense - expected)) > 2.0E-12_dp) ERROR STOP 'Canonical FFT cache scatter is incorrect'
         CALL cell_to_k_backend_evaluate_batch(snapshot, route, xkp(:, 3:5), 3, scratch, batch=dense(:, :, 3:5))
         IF (MAXVAL(ABS(dense(:, :, 3:5) - expected(:, :, 3:5))) > 2.0E-12_dp) THEN
            ERROR STOP 'Canonical FFT cache lost its local k offset'
         END IF
         CALL cell_to_k_backend_evaluate_batch(snapshot, route, xkp, 1, scratch, output=destinations)
         IF (MAXVAL(ABS(CMPLX(split(:, :, :, 1), split(:, :, :, 2), dp) - expected)) > 2.0E-12_dp) THEN
            ERROR STOP 'Canonical FFT cache split-FM scatter is incorrect'
         END IF
         DEALLOCATE (snapshot%fft_values)
         CALL cell_to_k_sparse_fft(snapshot, route, xkp, nkp_grid, RESHAPE([1, nkp], [2, 1]), env, &
                                   1_int_8, output=destinations)
         IF (MAXVAL(ABS(CMPLX(split(:, :, :, 1), split(:, :, :, 2), dp) - expected)) > 2.0E-12_dp) THEN
            ERROR STOP 'Direct FFT output differs from cached evaluation'
         END IF
      END DO
      CALL cell_to_k_grid_release(route%transform%fft_work)
      NULLIFY (route%matrix_struct)
      CALL mp_para_env_release(env)
   END BLOCK

   selected_kpoints = [64, 2, 31]
   ! A reduced request retains the full mesh spacing and shift, including a single point.
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                           selected_transformed, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Regular-grid subset did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 2.0E-12_dp) THEN
      ERROR STOP "Subset lattice FFT differs from direct transform"
   END IF
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                           selected_transformed(:, :, 1:2), used_fft, selected_kpoints=[3, 1])
   IF (.NOT. used_fft) ERROR STOP "Selection within a subset did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed(:, :, 1:2) - direct(:, :, selected_kpoints([3, 1])))) > 2.0E-12_dp) THEN
      ERROR STOP "Subset selection confused local and full-grid indices"
   END IF
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, 31:31), nkp_grid, selected_transformed(:, :, 1:1), used_fft)
   IF (.NOT. used_fft) ERROR STOP "Single requested grid point did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed(:, :, 1) - direct(:, :, 31))) > 2.0E-12_dp) THEN
      ERROR STOP "Single-point lattice FFT lost the mesh shift"
   END IF

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, selected_transformed, &
                           used_fft, 2, hmat, selected_kpoints=selected_kpoints)
   IF (.NOT. used_fft) ERROR STOP "Selected derivative did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct, 2, hmat)
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 5.0E-12_dp) THEN
      ERROR STOP "Selected lattice-FFT derivative differs from direct transform"
   END IF
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                           selected_transformed, used_fft, 2, hmat)
   IF (.NOT. used_fft) ERROR STOP "Subset derivative did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 5.0E-12_dp) THEN
      ERROR STOP "Subset lattice-FFT derivative differs from direct transform"
   END IF

   xkp_bad = xkp
   xkp_bad(1, nkp) = xkp_bad(1, nkp) + 0.03125_dp
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, used_fft)
   IF (used_fft) ERROR STOP "Irregular k-point list incorrectly used lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp_bad, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 1.0E-14_dp) THEN
      ERROR STOP "Direct lattice-transform fallback is inconsistent"
   END IF
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, used_fft, accumulate=.TRUE.)
   IF (used_fft .OR. MAXVAL(ABS(transformed - 2.0_dp*direct)) > 3.0E-12_dp) THEN
      ERROR STOP 'Direct lattice-transform fallback did not preserve accumulation'
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
   CALL k_grid_to_cell_execute(inverse_work, transformed(:, :, 1:9), values_real=real_transformed, used_fft=used_fft)
   IF (.NOT. used_fft .OR. MAXVAL(ABS(real_transformed - REAL(weighted_reference, dp))) > 3.0E-12_dp) THEN
      ERROR STOP 'Real-output weighted inverse differs from explicit phase sum'
   END IF
   CALL k_grid_to_cell_release(inverse_work)
   CALL k_grid_to_cell_prepare(inverse_work, xkp_bad, nkp_grid, index_to_cell)
   CALL k_grid_to_cell_execute(inverse_work, transformed, weighted_result, used_fft)
   IF (used_fft) ERROR STOP 'Irregular inverse incorrectly used lattice FFT'
   CALL k_grid_to_cell_execute(inverse_work, transformed, values_real=real_transformed, used_fft=used_fft)
   IF (used_fft .OR. MAXVAL(ABS(real_transformed - REAL(weighted_result, dp))) > 3.0E-12_dp) THEN
      ERROR STOP 'Real-output inverse fallback differs from the complex result'
   END IF
   CALL k_grid_to_cell_execute(inverse_work, transformed, values_real=real_transformed, &
                               used_fft=used_fft, allow_direct_fallback=.FALSE.)
   IF (used_fft .OR. ANY(real_transformed /= 0.0_dp)) THEN
      ERROR STOP 'Disabled real-output inverse fallback did not return zero'
   END IF
   CALL k_grid_to_cell_release(inverse_work)

   ! Local RSKP panels cover origin, repeated images, +/-R folding and shuffled subsets.
   index_to_cell(:, 1) = 0
   index_to_cell(:, 2) = index_to_cell(:, 3)
   index_to_cell(:, 4) = -index_to_cell(:, 3)
   CALL direct_reference(values_rs, index_to_cell, xkp, direct)
   weighted_reference = z_zero
   DO icell = 1, ncell
      DO ik = 1, nkp
         weighted_reference(:, :, icell) = weighted_reference(:, :, icell) + direct(:, :, ik)* &
                                           EXP(-gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), dp)))/REAL(nkp, dp)
      END DO
   END DO
   DO backend = 1, SIZE(backends)
      CALL cell_to_k_transform(values_rs, index_to_cell, xkp, nkp_grid, transformed, backends(backend))
      IF (MAXVAL(ABS(transformed - direct)) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP forward transform differs from explicit phase sum'
      END IF
      CALL cell_to_k_transform(values_rs, index_to_cell, xkp(:, selected_kpoints), nkp_grid, &
                               selected_transformed, backends(backend))
      IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP subset changed k-point order'
      END IF
      CALL cell_to_k_transform(values_rs, index_to_cell, -xkp, nkp_grid, transformed, backends(backend))
      IF (MAXVAL(ABS(transformed - CONJG(direct))) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP negative phase sign is inconsistent'
      END IF
      CALL cell_to_k_transform(values_rs(:, :, 1:1), index_to_cell(:, 1:1), xkp, nkp_grid, &
                               transformed, backends(backend))
      DO ik = 1, nkp
         IF (MAXVAL(ABS(transformed(:, :, ik) - values_rs(:, :, 1))) > 3.0E-12_dp) THEN
            ERROR STOP 'Local RSKP origin-only transform is inconsistent'
         END IF
      END DO
      CALL cell_to_k_transform(values_rs(:, :, 1:15), index_to_cell(:, 1:15), xkp, nkp_grid, &
                               transformed, backends(backend))
      CALL cell_to_k_transform(values_rs(:, :, 16:), index_to_cell(:, 16:), xkp, nkp_grid, &
                               transformed, backends(backend), accumulate=.TRUE.)
      IF (MAXVAL(ABS(transformed - direct)) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP image-panel accumulation differs from full transform'
      END IF
      CALL k_to_real_cell_transform(direct, xkp, nkp_grid, index_to_cell, real_transformed, backends(backend))
      IF (MAXVAL(ABS(real_transformed - REAL(weighted_reference, dp))) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP inverse changed the phase sign or normalization'
      END IF
      CALL cell_to_k_transform(values_rs(1:0, :, :), index_to_cell, xkp, nkp_grid, &
                               transformed(1:0, :, :), backends(backend))
      CALL k_to_real_cell_transform(direct(1:0, :, :), xkp, nkp_grid, index_to_cell, &
                                    real_transformed(1:0, :, :), backends(backend))
   END DO
   CALL direct_reference(values_rs, index_to_cell, xkp_bad, direct)
   CALL cell_to_k_transform(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, rskp_backend_gemm)
   IF (MAXVAL(ABS(transformed - direct)) > 3.0E-12_dp) THEN
      ERROR STOP 'Local RSKP irregular-grid GEMM differs from explicit phase sum'
   END IF

   ! A final, partial matrix tile must use its own row count while retaining the packed stride.
   ALLOCATE (panel_input(4097, 1, 3), panel_inverse(4097, 1, 3), &
             panel_output(4097, 1, nkp), panel_reference(4097, 1, nkp))
   DO icell = 1, 3
      panel_input(:, 1, icell) = [(SIN(0.017_dp*REAL(i + icell, dp)), i=1, 4097)]
   END DO
   CALL direct_reference(panel_input, index_to_cell(:, 1:3), xkp, panel_reference)
   DO backend = 1, SIZE(backends)
      CALL cell_to_k_transform(panel_input, index_to_cell(:, 1:3), xkp, nkp_grid, &
                               panel_output, backends(backend))
      IF (MAXVAL(ABS(panel_output - panel_reference)) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP forward partial row tile is inconsistent'
      END IF
      CALL k_to_real_cell_transform(panel_output, xkp, nkp_grid, index_to_cell(:, 1:3), &
                                    panel_inverse, backends(backend))
      ! Images 2 and 3 coincide, so their coefficients add under inverse reconstruction.
      IF (MAXVAL(ABS(panel_inverse(:, :, 1) - panel_input(:, :, 1))) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP inverse partial row tile is inconsistent'
      END IF
      IF (MAXVAL(ABS(panel_inverse(:, :, 2) - SUM(panel_input(:, :, 2:3), DIM=3))) > 3.0E-12_dp) THEN
         ERROR STOP 'Local RSKP inverse lost repeated images'
      END IF
   END DO
   DEALLOCATE (panel_input, panel_inverse, panel_output, panel_reference)

   ! Exercise a final partial cell tile after the 256-column GEMM tile.
   ALLOCATE (panel_inverse(2, 3, SIZE(panel_cells, 2)))
   DO icell = 1, SIZE(panel_cells, 2)
      panel_cells(:, icell) = index_to_cell(:, MODULO(icell - 1, ncell) + 1)
   END DO
   CALL direct_reference(values_rs, index_to_cell, xkp, direct)
   DO backend = 1, SIZE(backends)
      CALL k_to_real_cell_transform(direct, xkp, nkp_grid, panel_cells, panel_inverse, backends(backend))
      DO icell = 1, SIZE(panel_cells, 2)
         IF (MAXVAL(ABS(panel_inverse(:, :, icell) - &
                        REAL(weighted_reference(:, :, MODULO(icell - 1, ncell) + 1), dp))) > 3.0E-12_dp) THEN
            ERROR STOP 'Local RSKP inverse partial cell tile is inconsistent'
         END IF
      END DO
   END DO
   DEALLOCATE (panel_inverse)

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
      USE cp_fm_struct, ONLY: cp_fm_struct_type
      USE cp_fm_types, ONLY: cp_fm_p_type, cp_fm_type
      USE cp_dbcsr_api, ONLY: &
         dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
         dbcsr_finalize, dbcsr_finalize_lib, dbcsr_init_lib, dbcsr_p_type, dbcsr_put_block, &
         dbcsr_release, dbcsr_type_antisymmetric, dbcsr_type_symmetric
      USE cp_log_handling, ONLY: cp_add_default_logger, cp_logger_create, cp_logger_release, &
                                 cp_logger_type, cp_rm_default_logger
      USE kpoint_k_r_trafo_simple, ONLY: replicate_rs_matrices
      USE kinds, ONLY: int_8
      USE kpoint_transform, ONLY: k_grid_to_cell_contract, k_grid_to_cell_phase_prepare, &
                                  k_grid_to_cell_phase_type, k_grid_to_cell_profile_type
      USE kpoint_types, ONLY: kpoint_type
      USE local_gemm_api, ONLY: LOCAL_GEMM_PU_HOST, local_gemm_ctxt_type
      USE machine, ONLY: default_output_unit
      USE message_passing, ONLY: mp_para_env_release, mp_para_env_type
      USE qs_neighbor_list_types, ONLY: add_neighbor_list, add_neighbor_node, &
                                        allocate_neighbor_list_set, deallocate_neighbor_list_set, neighbor_list_p_type
!$    USE omp_lib, ONLY: omp_get_max_threads
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
      BLOCK
         INTEGER, PARAMETER :: nk = 1031, nimages = 11
         INTEGER :: c, e, g, j, k, layout(4, 2, nk, 2, 2), nc, request, s, shift(3), worker, workers
         INTEGER :: image_data(5, nimages), k_index(nk)
         INTEGER(KIND=int_8) :: gathered, expected_gather
         INTEGER, TARGET :: row_map(5), col_map(5)
         LOGICAL :: reverse, transposed
         REAL(dp) :: arg, cosine, sine, re(6), im(6), raw(6, 2), coordinates(3, nk), weight(nk)
         REAL(dp) :: reference(6, nimages, 2)
         REAL(dp), ALLOCATABLE, TARGET :: data(:, :, :, :, :)
         REAL(dp), ALLOCATABLE :: result(:, :, :)
         TYPE(cp_fm_struct_type), TARGET :: structure
         TYPE(cp_fm_type), ALLOCATABLE, TARGET :: matrices(:, :, :)
         TYPE(cp_fm_p_type), ALLOCATABLE :: views(:, :, :)
         TYPE(local_gemm_ctxt_type), ALLOCATABLE :: contexts(:)
         TYPE(k_grid_to_cell_profile_type), ALLOCATABLE :: stats(:)
         TYPE(k_grid_to_cell_phase_type) :: phase_table

         ! Borrowed FM views over test-owned storage; no owning FM release is used.
         ALLOCATE (DATA(5, 5, nk, 2, 2), matrices(nk, 2, 2), views(nk, 2, 2))
         structure%row_indices => row_map
         structure%col_indices => col_map
         structure%local_leading_dimension = 5
         DO s = 1, 2
            DO c = 1, 2
               DO k = 1, nk
                  matrices(k, c, s)%matrix_struct => structure
                  matrices(k, c, s)%local_data => DATA(:, :, k, c, s)
                  views(k, c, s)%matrix => matrices(k, c, s)
                  DO j = 1, 5
                     DO atom = 1, 5
                        DATA(atom, j, k, c, s) = 0.02_dp*REAL(s, dp) + 0.003_dp*REAL(atom - j, dp) + &
                                                 0.1_dp*SIN(0.017_dp*REAL(k, dp) + REAL(c, dp))
                     END DO
                  END DO
               END DO
            END DO
         END DO
         DO k = 1, nk
            k_index(k) = nk + 1 - k
            coordinates(:, k) = REAL(k, dp)*[0.0017_dp, -0.0009_dp, 0.0003_dp]
            weight(k) = (1.0_dp + 0.05_dp*REAL(MOD(k, 5), dp))/REAL(nk, dp)
         END DO
         image_data = RESHAPE([0, 0, 0, 2, 0, 1, 0, 0, 1, 1, -1, 0, 0, 2, -1, &
                               0, 1, 0, 3, 0, 0, -1, 0, 1, 2, 1, 1, 0, 1, 1, -1, -1, 0, 1, -1, &
                               0, 0, 1, 2, 1, 0, 0, -1, 1, 1, 1, 0, 0, 2, 2, 17, -19, 0, 1, -2], [5, nimages])
         DO request = 1, 2
            workers = 1
!$          IF (request == 2) workers = MIN(4, omp_get_max_threads())
            ALLOCATE (contexts(workers))
            ALLOCATE (stats(workers), source=k_grid_to_cell_profile_type())
            DO worker = 1, workers
               contexts(worker) = local_gemm_ctxt_type()
               CALL contexts(worker)%create(LOCAL_GEMM_PU_HOST, timing=.FALSE.)
            END DO
            CALL k_grid_to_cell_phase_prepare(phase_table, image_data(1:3, :), coordinates, stats)
            gathered = 0_int_8
            DO worker = 1, workers
               gathered = gathered + stats(worker)%counts(2)
            END DO
            ! Distinct nonzero magnitudes: x={1,17}, y={1,19}, z={1}.
            IF (gathered /= 10_int_8*INT(nk, int_8)) ERROR STOP 'Axis phase table did not deduplicate coordinates'
            DO mode = 1, 5
               stats(:) = k_grid_to_cell_profile_type()
               nc = MERGE(1, 2, mode == 3)
               reverse = mode /= 1
               transposed = mode == 4
               shift = [1, -2, 1]
               IF (nc == 1 .OR. mode == 5) shift = 0
               row_map = [1, 2, 3, 4, 5]
               col_map = [1, 2, 3, 4, 5]
               DO s = 1, 2
                  DO c = 1, 2
                     DO k = 1, nk
                        layout(:, 1, k, c, s) = [1, 3, 1, 3]
                        layout(:, 2, k, c, s) = [4, 5, 4, 5]
                        IF (mode == 2 .AND. k <= nk/2) layout(2, 1, k, c, s) = 0
                        IF (mode == 2 .AND. s == 2) layout(2, 1, k, c, s) = 0
                        IF (mode == 4) layout(2, 1, k, c, s) = 2
                     END DO
                  END DO
               END DO
               IF (mode == 4) row_map = [1, 3, 4, 5, 6]
               reference = 0.0_dp
               expected_gather = 0_int_8
               DO k = 1, nk
                  g = k_index(k)
                  arg = twopi*SUM(REAL(shift, dp)*coordinates(:, k))
                  cosine = COS(arg)
                  sine = SIN(arg)
                  DO s = 1, 2
                     raw = 0.0_dp
                     DO c = 1, nc
                        DO j = layout(3, 2, g, c, s), layout(4, 2, g, c, s)
                           DO atom = layout(1, 1, g, c, s), layout(2, 1, g, c, s)
                              e = row_map(atom) + 3*(col_map(j) - 4)
                              raw(e, c) = DATA(atom, j, g, c, s)
                              expected_gather = expected_gather + 1_int_8
                           END DO
                        END DO
                     END DO
                     re = cosine*raw(:, 1) + MERGE(-1.0_dp, 1.0_dp, reverse)*sine*raw(:, 2)
                     im = -sine*raw(:, 1) + MERGE(-1.0_dp, 1.0_dp, reverse)*cosine*raw(:, 2)
                     DO img = 1, nimages
                        arg = twopi*SUM(REAL(image_data(1:3, img), dp)*coordinates(:, k))
                        reference(:, img, s) = reference(:, img, s) + weight(k)*( &
                                               REAL(image_data(4, img), dp)*re*COS(arg) + &
                                               MERGE(-1.0_dp, 1.0_dp, transposed)* &
                                               REAL(image_data(5, img), dp)*im*SIN(arg))
                     END DO
                  END DO
               END DO
               CALL k_grid_to_cell_contract(views(:, 1:nc, :), layout(:, :, :, 1:nc, :), k_index, &
                                            coordinates, weight, image_data, [1, 2], [3, 2], [1, 4], &
                                            shift, reverse, transposed, phase_table, RESULT, contexts, stats)
               IF (MAXVAL(ABS(RESULT - reference)) > 2.0E-11_dp) THEN
                  ERROR STOP 'Parallel k-batch density contraction differs from explicit phase sum'
               END IF
               gathered = 0_int_8
               DO worker = 1, workers
                  gathered = gathered + stats(worker)%counts(3)
               END DO
               IF (gathered /= expected_gather) ERROR STOP 'K-batch contraction gathered a source element more than once'
               gathered = 0_int_8
               DO worker = 1, workers
                  gathered = gathered + stats(worker)%counts(2)
               END DO
               IF (ALL(shift == 0)) THEN
                  IF (gathered /= 0_int_8) ERROR STOP 'Zero-shift contraction recomputed trigonometric phases'
               ELSE
                  IF (gathered /= 2_int_8*INT(nk, int_8)) ERROR STOP 'Contraction recomputed image phases'
               END IF
            END DO
            DO worker = 1, workers
               CALL contexts(worker)%destroy()
            END DO
            DEALLOCATE (contexts, stats)
         END DO
      END BLOCK
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
