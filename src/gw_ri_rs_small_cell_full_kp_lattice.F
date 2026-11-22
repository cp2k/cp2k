!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI lattice matrices, real-space cells and data distribution shared by χ, W, Σ and the
!>        Gamma correction of periodic RI-RS GW.
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE gw_ri_rs_small_cell_full_kp_lattice
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_release,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_dist2d_to_dist
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE input_constants,                 ONLY: do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              dp_size,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: pi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp_lattice'

   PUBLIC :: assemble_AA_contraction, &
             block_range, &
             build_full_cell_grid, &
             build_ri_2c_lattice, &
             complex_bytes, &
             gamma_auxiliary_charges, &
             grid_pair_indices, &
             grid_pair_gemm_panel, &
             grid_pair_upper_multiply, &
             grid_pair_shift, &
             grid_pair_stencil_weight, &
             image_weight, &
             owned_mo_coefficients, &
             real_bytes, &
             ri_2c_lattice_to_kp, &
             ri_metric_head_zero

   ! Bytes per stored real and complex number.
   REAL(KIND=dp), PARAMETER :: real_bytes = REAL(dp_size, dp), complex_bytes = 2.0_dp*real_bytes
   INTEGER, PARAMETER, PUBLIC :: grid_pair_stencil_denominator = 4**3
   INTEGER, PARAMETER, PUBLIC :: grid_pair_panel_cols = 32

CONTAINS

! **************************************************************************************************
!> \brief Decode an unordered grid pair p=i+j*(j-1)/2, 1<=i<=j.
!> Integer corrections keep the triangular boundaries exact for large grids.
!> \param ipair packed upper-triangle index
!> \param il row
!> \param ilp column
! **************************************************************************************************
   PURE SUBROUTINE grid_pair_indices(ipair, il, ilp)
      INTEGER, INTENT(IN)                                :: ipair
      INTEGER, INTENT(OUT)                               :: il, ilp

      INTEGER(KIND=int_8)                                :: column, pair

      pair = INT(ipair, int_8)
      column = INT((1.0_dp + SQRT(8.0_dp*REAL(pair, dp) - 7.0_dp))/2.0_dp, int_8)
      DO WHILE (column*(column + 1_int_8)/2_int_8 < pair)
         column = column + 1_int_8
      END DO
      DO WHILE (column*(column - 1_int_8)/2_int_8 >= pair)
         column = column - 1_int_8
      END DO
      ilp = INT(column)
      il = INT(pair - column*(column - 1_int_8)/2_int_8)
   END SUBROUTINE grid_pair_indices

! **************************************************************************************************
!> \brief Nearest lattice image of an unordered grid pair, including both half-cell ties.
!> \param il row of the unordered pair
!> \param ilp column of the unordered pair
!> \param fractional_grid periodic fractional grid coordinates
!> \param pair_shift lower nearest image on a tie
!> \param pair_tie directions with two equidistant images
! **************************************************************************************************
   PURE SUBROUTINE grid_pair_shift(il, ilp, fractional_grid, pair_shift, pair_tie)
      INTEGER, INTENT(IN)                                :: il, ilp
      REAL(KIND=dp), INTENT(IN)                          :: fractional_grid(:, :)
      INTEGER, INTENT(OUT)                               :: pair_shift(3)
      LOGICAL, INTENT(OUT)                               :: pair_tie(3)

      INTEGER                                            :: idir
      REAL(KIND=dp)                                      :: delta(3), half_boundary, tolerance

      delta = fractional_grid(:, il) - fractional_grid(:, ilp)
      pair_shift = NINT(delta)
      pair_tie = .FALSE.
      DO idir = 1, 3
         half_boundary = ANINT(delta(idir) - 0.5_dp) + 0.5_dp
         tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(delta(idir)))
         IF (ABS(delta(idir) - half_boundary) <= tolerance) THEN
            pair_shift(idir) = NINT(half_boundary - 0.5_dp)
            pair_tie(idir) = .TRUE.
         END IF
      END DO
   END SUBROUTINE grid_pair_shift

! **************************************************************************************************
!> \brief Compute a column panel of A B^H through its last diagonal.
!> \param n grid points
!> \param ninner contracted dimension
!> \param first first column of the panel
!> \param width number of columns
!> \param left A
!> \param right B; the caller supplies a Hermitian product
!> \param panel ...
! **************************************************************************************************
   SUBROUTINE grid_pair_gemm_panel(n, ninner, first, width, left, right, panel)
      INTEGER, INTENT(IN)                                :: n, ninner, first, width
      COMPLEX(KIND=dp), INTENT(IN)                       :: left(n, ninner), right(n, ninner)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: panel(n, width)

      INTEGER                                            :: col, row

      CALL zgemm('N', 'C', first + width - 1, width, ninner, z_one, left, n, &
                 right(first, 1), n, z_zero, panel, n)
      DO col = 1, width
         row = first + col - 1
         panel(row, col) = CMPLX(REAL(panel(row, col), dp), 0.0_dp, dp)
      END DO
   END SUBROUTINE grid_pair_gemm_panel

! **************************************************************************************************
!> \brief Y=U Z with U=triu(H,1)+diag(H)/2, so that H=U+U^H.
!> Hermitian projections reconstruct the adjoint after contraction, saving the lower multiply.
!> \param n grid points
!> \param nrhs right-hand sides
!> \param pairs packed upper triangle of H
!> \param Z dense right-hand sides
!> \param Y result
!> \param panel reusable thread-local workspace
! **************************************************************************************************
   SUBROUTINE grid_pair_upper_multiply(n, nrhs, pairs, Z, Y, panel)
      INTEGER, INTENT(IN)                                :: n, nrhs
      COMPLEX(KIND=dp), INTENT(IN)                       :: pairs(:), Z(n, nrhs)
      COMPLEX(KIND=dp), INTENT(OUT) :: Y(n, nrhs), panel(n, MIN(grid_pair_panel_cols, n))

      INTEGER                                            :: col, first, j, offset, width

      CPASSERT(SIZE(pairs, KIND=int_8) == INT(n, int_8)*(n + 1)/2)
      Y = z_zero
      DO first = 1, n, grid_pair_panel_cols
         width = MIN(grid_pair_panel_cols, n - first + 1)
         DO col = 1, width
            j = first + col - 1
            offset = INT(INT(j, int_8)*INT(j - 1, int_8)/2_int_8)
            panel(1:j - 1, col) = pairs(offset + 1:offset + j - 1)
            panel(j, col) = CMPLX(0.5_dp*REAL(pairs(offset + j), dp), 0.0_dp, dp)
            panel(j + 1:first + width - 1, col) = z_zero
         END DO
         CALL zgemm('N', 'N', first + width - 1, nrhs, width, z_one, panel, n, Z(first, 1), n, z_one, Y, n)
      END DO
   END SUBROUTINE grid_pair_upper_multiply

! **************************************************************************************************
!> \brief Finite-supercell image weight, averaged over equidistant pair centers.
!> \param c image relative to the lower pair center
!> \param n_mesh mesh length
!> \param tie two equidistant centers
!> \param boundary weight of a supercell endpoint
!> \return image weight
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION image_weight(c, n_mesh, tie, boundary) RESULT(weight)
      INTEGER, INTENT(IN)                                :: c, n_mesh
      LOGICAL, INTENT(IN)                                :: tie
      REAL(KIND=dp), INTENT(IN)                          :: boundary

      INTEGER                                            :: i_center, n_centers

      n_centers = MERGE(2, 1, tie)
      weight = 0.0_dp
      DO i_center = 0, n_centers - 1
         IF (2*ABS(c - i_center) < n_mesh) THEN
            weight = weight + 1.0_dp/REAL(n_centers, dp)
         ELSE IF (2*ABS(c - i_center) == n_mesh) THEN
            weight = weight + boundary/REAL(n_centers, dp)
         END IF
      END DO
   END FUNCTION image_weight

! **************************************************************************************************
!> \brief Signed stencil numerator on an inversion-closed mesh, with denominator 64.
!> All centered images share exp(2*pi*i*offset*cell); their remaining phases are signs.
!> Each one-axis image weight is a multiple of 1/4, including half-cell averaging.
!> \param cell canonical source cell
!> \param shift lower pair center
!> \param tie half-cell ties
!> \param nmesh Born-von Karman mesh
!> \param half_shift twice the mesh-scaled reciprocal offset
!> \param displacement grid-point displacement in Cartesian coordinates
!> \param hmat lattice matrix
!> \param cutoff_squared squared distance cutoff; nonpositive disables truncation
!> \return signed image-weight numerator
! **************************************************************************************************
   PURE INTEGER FUNCTION grid_pair_stencil_weight(cell, shift, tie, nmesh, half_shift, &
                                                  displacement, hmat, cutoff_squared) RESULT(numerator)
      INTEGER, INTENT(IN)                                :: cell(3), shift(3)
      LOGICAL, INTENT(IN)                                :: tie(3)
      INTEGER, INTENT(IN)                                :: nmesh(3), half_shift(3)
      REAL(KIND=dp), INTENT(IN)                          :: displacement(3), hmat(3, 3), &
                                                            cutoff_squared

      INTEGER                                            :: base(3), ix, iy, iz, translated(3), &
                                                            weight, weights(0:1, 3)
      REAL(KIND=dp)                                      :: boundary(3), distance(3)

      ! Periodic image weights partition unity when every supercell phase is one.
      IF (cutoff_squared <= 0.0_dp .AND. ALL(MODULO(half_shift, 2) == 0)) THEN
         numerator = grid_pair_stencil_denominator
         RETURN
      END IF
      base = MODULO(cell - shift + nmesh/2, nmesh) - nmesh/2
      boundary = 0.5_dp*REAL(1 - MODULO(half_shift, 2), dp)
      DO ix = 1, 3
         weights(:, ix) = NINT(4.0_dp*[image_weight(base(ix), nmesh(ix), tie(ix), boundary(ix)), &
                                       image_weight(base(ix) + nmesh(ix), nmesh(ix), tie(ix), boundary(ix))])
      END DO
      ! Without a distance cutoff, the image sum separates into three one-axis sums.
      IF (cutoff_squared <= 0.0_dp) THEN
         numerator = PRODUCT(weights(0, :) + (1 - 2*MODULO(half_shift, 2))*weights(1, :))
         IF (MODULO(DOT_PRODUCT(half_shift, (shift + base - cell)/nmesh), 2) /= 0) numerator = -numerator
         RETURN
      END IF
      numerator = 0
      DO iz = 0, 1
         IF (weights(iz, 3) == 0) CYCLE
         DO iy = 0, 1
            IF (weights(iy, 2) == 0) CYCLE
            DO ix = 0, 1
               weight = weights(ix, 1)*weights(iy, 2)*weights(iz, 3)
               IF (weight == 0) CYCLE
               translated = shift + base + nmesh*[ix, iy, iz]
               distance = displacement - MATMUL(hmat, REAL(translated, dp))
               IF (SUM(distance**2) > cutoff_squared) CYCLE
               IF (MODULO(DOT_PRODUCT(half_shift, (translated - cell)/nmesh), 2) /= 0) weight = -weight
               numerator = numerator + weight
            END DO
         END DO
      END DO
   END FUNCTION grid_pair_stencil_weight

! **************************************************************************************************
!> \brief Assembles the q-independent, atom-diagonal AA contraction from stored coefficients.
!> \param bs_env ...
!> \param U_Pp Atom-diagonal reference-to-optimized AA contraction.
! **************************************************************************************************
   SUBROUTINE assemble_AA_contraction(bs_env, U_Pp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: U_Pp(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'assemble_AA_contraction'

      INTEGER                                            :: first, first_opt, first_ref, handle, &
                                                            iatom, last, nopt, nref

      CALL timeset(routineN, handle)

      IF (.NOT. bs_env%auto_ri%ready .OR. &
          .NOT. ALLOCATED(bs_env%auto_ri%sizes_ref_RI)) THEN
         CPABORT('Periodic AUTO_RI contraction is not initialized')
      END IF
      ALLOCATE (U_Pp(SUM(bs_env%auto_ri%sizes_ref_RI), &
                     SUM(bs_env%auto_ri%sizes_opt_RI)), source=z_zero)
      first_ref = 1
      first_opt = 1
      DO iatom = 1, bs_env%n_atom
         nref = bs_env%auto_ri%sizes_ref_RI(iatom)
         nopt = bs_env%auto_ri%sizes_opt_RI(iatom)
         first = bs_env%auto_ri%U_Pp_atom_offset(iatom)
         last = first + nref*nopt - 1
         U_Pp(first_ref:first_ref + nref - 1, first_opt:first_opt + nopt - 1) = &
            RESHAPE(bs_env%auto_ri%U_Pp_atom(first:last), [nref, nopt])
         first_ref = first_ref + nref
         first_opt = first_opt + nopt
      END DO
      IF (ANY(SUM(ABS(U_Pp)**2, DIM=1) <= TINY(1.0_dp))) THEN
         CPABORT('Periodic AUTO_RI contraction has an empty optimized direction')
      END IF

      CALL timestop(handle)

   END SUBROUTINE assemble_AA_contraction

! **************************************************************************************************
!> \brief Dense MO coefficients of the cyclically owned SCF k points (owner of k: MODULO(k-1, ranks)).
!> Pack the local BLACS tiles by destination in bounded k-point batches. The common matrix
!> structure supplies the row and column indices for unpacking each sender's tile.
!> \param bs_env ...
!> \param mo_coeff distributed MO coefficients per k point and spin
!> \param coeff coeff(:, :, k_local, spin) for the owned k points
! **************************************************************************************************
   SUBROUTINE owned_mo_coefficients(bs_env, mo_coeff, coeff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), DIMENSION(:, :), INTENT(IN)     :: mo_coeff
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: coeff(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'owned_mo_coefficients'
      INTEGER(KIND=int_8), PARAMETER :: buffer_elements = 64000000_int_8/INT(2*dp_size, int_8)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: received(:), sent(:)
      INTEGER :: batch_size, col, first_k, first_owned, first_owned_k, handle, i_col, i_owned, &
         i_rank, i_row, ikp, ipcol, iprow, ispin, last_k, my_rank, n_ao, n_k, n_mo, n_owned, &
         n_owned_batch, ncol_local, ncol_src, nrow_local, nrow_src, num_pe, offset
      INTEGER(KIND=int_8)                                :: matrix_size, max_tile_size
      INTEGER, ALLOCATABLE                               :: cols(:, :), recv_counts(:), &
                                                            recv_displs(:), rows(:, :), &
                                                            send_counts(:), send_displs(:), &
                                                            tile_sizes(:)
      INTEGER, POINTER                                   :: local_cols(:), local_rows(:)
      LOGICAL                                            :: same_layout
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: layout
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      n_k = SIZE(mo_coeff, 1)
      para_env => bs_env%para_env
      num_pe = para_env%num_pe
      my_rank = para_env%mepos
      CALL cp_cfm_get_info(mo_coeff(1, 1), nrow_global=n_ao, ncol_global=n_mo, nrow_local=nrow_local, &
                           ncol_local=ncol_local, context=context, matrix_struct=layout)
      CPASSERT(layout%para_env%num_pe == num_pe)
      CPASSERT(layout%para_env%mepos == my_rank)
      CALL cp_fm_struct_get(layout, nrow_locals=local_rows, ncol_locals=local_cols)
      matrix_size = INT(n_ao, int_8)*INT(n_mo, int_8)
      max_tile_size = INT(MAXVAL(local_rows), int_8)*INT(MAXVAL(local_cols), int_8)
      CPASSERT(matrix_size > 0_int_8 .AND. max_tile_size > 0_int_8)
      CPASSERT(matrix_size + max_tile_size <= INT(HUGE(0), int_8))
      n_owned = (n_k + num_pe - 1 - my_rank)/num_pe
      batch_size = INT(MIN(INT(n_k, int_8), INT(num_pe, int_8)*MAX(1_int_8, buffer_elements/matrix_size), &
                           MAX(1_int_8, buffer_elements/max_tile_size)))
      ALLOCATE (coeff(n_ao, n_mo, n_owned, SIZE(mo_coeff, 2)), &
                sent(MAX(1, nrow_local*ncol_local*batch_size)), &
                received(MAX(1, INT(matrix_size)*MIN(n_owned, (batch_size + num_pe - 1)/num_pe))))
      ALLOCATE (rows(MAXVAL(local_rows), 0:context%num_pe(1) - 1), &
                cols(MAXVAL(local_cols), 0:context%num_pe(2) - 1), &
                recv_counts(0:num_pe - 1), recv_displs(0:num_pe - 1), &
                send_counts(0:num_pe - 1), send_displs(0:num_pe - 1), tile_sizes(0:num_pe - 1))
      DO iprow = 0, context%num_pe(1) - 1
         DO i_row = 1, local_rows(iprow)
            rows(i_row, iprow) = layout%l2g_row(i_row, iprow)
         END DO
      END DO
      DO ipcol = 0, context%num_pe(2) - 1
         DO i_col = 1, local_cols(ipcol)
            cols(i_col, ipcol) = layout%l2g_col(i_col, ipcol)
         END DO
      END DO
      DO i_rank = 0, num_pe - 1
         iprow = context%mpi2blacs(1, i_rank)
         ipcol = context%mpi2blacs(2, i_rank)
         tile_sizes(i_rank) = local_rows(iprow)*local_cols(ipcol)
      END DO
      CPASSERT(SUM(tile_sizes) == INT(matrix_size))

      DO ispin = 1, SIZE(mo_coeff, 2)
         DO first_k = 1, n_k, batch_size
            last_k = MIN(first_k + batch_size - 1, n_k)
            first_owned_k = first_k + MODULO(my_rank + 1 - first_k, num_pe)
            first_owned = (first_owned_k - 1)/num_pe + 1
            n_owned_batch = (last_k - first_owned_k + num_pe)/num_pe
            recv_counts(:) = n_owned_batch*tile_sizes
            recv_displs(0) = 0
            DO i_rank = 1, num_pe - 1
               recv_displs(i_rank) = recv_displs(i_rank - 1) + recv_counts(i_rank - 1)
            END DO

            offset = 0
            DO i_rank = 0, num_pe - 1
               send_displs(i_rank) = offset
               DO ikp = first_k + MODULO(i_rank + 1 - first_k, num_pe), last_k, num_pe
                  same_layout = ASSOCIATED(mo_coeff(ikp, ispin)%matrix_struct, layout)
                  CPASSERT(same_layout)
                  DO i_col = 1, ncol_local
                     sent(offset + 1:offset + nrow_local) = mo_coeff(ikp, ispin)%local_data(:nrow_local, i_col)
                     offset = offset + nrow_local
                  END DO
               END DO
               send_counts(i_rank) = offset - send_displs(i_rank)
            END DO

            CALL para_env%alltoall(sent, send_counts, send_displs, received, recv_counts, recv_displs)
            DO i_rank = 0, num_pe - 1
               IF (recv_counts(i_rank) == 0) CYCLE
               iprow = context%mpi2blacs(1, i_rank)
               ipcol = context%mpi2blacs(2, i_rank)
               nrow_src = local_rows(iprow)
               ncol_src = local_cols(ipcol)
               offset = recv_displs(i_rank)
               DO i_owned = first_owned, first_owned + n_owned_batch - 1
                  DO i_col = 1, ncol_src
                     col = cols(i_col, ipcol)
                     coeff(rows(:nrow_src, iprow), col, i_owned, ispin) = received(offset + 1:offset + nrow_src)
                     offset = offset + nrow_src
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE owned_mo_coefficients

! **************************************************************************************************
!> \brief Construct a complete cell grid dual to the uniform SCF k-point mesh.
!> \param bs_env ...
!> \param index_to_cell ...
! **************************************************************************************************
   SUBROUTINE build_full_cell_grid(bs_env, index_to_cell)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: index_to_cell(:, :)

      INTEGER                                            :: cell(3), i1, i2, i3, i_cell, &
                                                            i_direction, lower(3), upper(3)

      DO i_direction = 1, 3
         lower(i_direction) = -bs_env%kpoints_scf_desymm%nkp_grid(i_direction)/2
         upper(i_direction) = lower(i_direction) + &
                              bs_env%kpoints_scf_desymm%nkp_grid(i_direction) - 1
      END DO
      ALLOCATE (index_to_cell(3, PRODUCT(bs_env%kpoints_scf_desymm%nkp_grid)))
      i_cell = 0
      DO i1 = lower(1), upper(1)
         DO i2 = lower(2), upper(2)
            DO i3 = lower(3), upper(3)
               i_cell = i_cell + 1
               cell = [i1, i2, i3]
               index_to_cell(:, i_cell) = cell
            END DO
         END DO
      END DO

   END SUBROUTINE build_full_cell_grid

! **************************************************************************************************
!> \brief Zero-wave-vector Fourier value of the selected RI metric kernel.
!> \param bs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION ri_metric_head_zero(bs_env) RESULT(head)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: head

      IF (bs_env%ri_metric%potential_type == do_potential_short) THEN
         IF (bs_env%ri_metric%omega <= 0.0_dp) THEN
            CPABORT('RI attenuation must be positive')
         END IF
         head = pi/bs_env%ri_metric%omega**2
      ELSE
         IF (bs_env%ri_metric%potential_type /= do_potential_truncated) THEN
            CPABORT('Unsupported RI metric')
         END IF
         head = twopi*bs_env%ri_metric%cutoff_radius**2
      END IF
   END FUNCTION ri_metric_head_zero

! **************************************************************************************************
!> \brief Form distributed I_PQ(R) = (φ_P,0 | v | φ_Q,R) over every screened neighbor cell.
!>
!> The finite-range potential v and the Gaussian screening radii determine the cell set.
!> No cell is dropped because it lies outside the SCF Born-von Karman cell. Each real matrix
!> remains distributed over the existing RI matrix communicator; callers release fm_R.
!> The reciprocal interaction is I_PQ(q) = sum_R exp(i q.R) I_PQ(R).
!> \param fm_R ...
!> \param cells ...
!> \param potential ...
!> \param regularization ...
!> \param bs_env ...
!> \param qs_env ...
!> \param reference_basis ...
! **************************************************************************************************
   SUBROUTINE build_ri_2c_lattice(fm_R, cells, potential, regularization, bs_env, qs_env, reference_basis)
      TYPE(cp_fm_type), ALLOCATABLE, INTENT(OUT)         :: fm_R(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells(:, :)
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential
      REAL(KIND=dp), INTENT(IN)                          :: regularization
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: reference_basis

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_ri_2c_lattice'

      INTEGER                                            :: bounds(3), cell_index(3), handle, icell, &
                                                            ix, iy, iz, ncells
      INTEGER, ALLOCATABLE                               :: occupied(:, :, :), occupied_flat(:), &
                                                            sizes(:)
      INTEGER, CONTIGUOUS, POINTER                       :: col_sizes(:), row_sizes(:)
      LOGICAL                                            :: use_reference_basis
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(cp_fm_struct_type), POINTER                   :: reference_matrix_struct
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_type), ALLOCATABLE                      :: matrices(:)
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(kpoint_type), POINTER                         :: kpoints_ptr
      TYPE(kpoint_type), TARGET                          :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para
      TYPE(neighbor_list_iterator_p_type), POINTER       :: iterator(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: neighbors(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL timeset(routineN, handle)

      use_reference_basis = .FALSE.
      IF (PRESENT(reference_basis)) use_reference_basis = reference_basis
      NULLIFY (blacs, para, reference_matrix_struct)
      CALL get_qs_env(qs_env, distribution_2d=dist_2d, particle_set=particles, qs_kind_set=qs_kinds, &
                      blacs_env=blacs, para_env=para)
      ALLOCATE (sizes(bs_env%n_atom))
      CALL get_particle_set(particles, qs_kinds, nsgf=sizes, basis=bs_env%basis_set_RI)
      NULLIFY (neighbors, iterator)
      CALL build_2c_neighbor_lists(neighbors, bs_env%basis_set_RI, bs_env%basis_set_RI, &
                                   potential, "full_2c_RI", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)
      bounds = 0
      CALL neighbor_list_iterator_create(iterator, neighbors)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, cell=cell_index)
         bounds = MAX(bounds, ABS(cell_index))
      END DO
      CALL neighbor_list_iterator_release(iterator)
      CALL bs_env%para_env%max(bounds)
      ALLOCATE (occupied(-bounds(1):bounds(1), -bounds(2):bounds(2), -bounds(3):bounds(3)))
      occupied = 0
      CALL neighbor_list_iterator_create(iterator, neighbors)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, cell=cell_index)
         occupied(cell_index(1), cell_index(2), cell_index(3)) = 1
      END DO
      CALL neighbor_list_iterator_release(iterator)
      ALLOCATE (occupied_flat(SIZE(occupied)))
      occupied_flat(:) = RESHAPE(occupied, [SIZE(occupied)])
      CALL bs_env%para_env%max(occupied_flat)
      occupied(:, :, :) = RESHAPE(occupied_flat, SHAPE(occupied))
      DEALLOCATE (occupied_flat)
      ncells = SUM(occupied)
      CPASSERT(ncells > 0)
      ALLOCATE (cells(3, ncells), fm_R(ncells))
      ALLOCATE (kpoints%cell_to_index(-bounds(1):bounds(1), -bounds(2):bounds(2), -bounds(3):bounds(3)))
      kpoints%cell_to_index = 0
      icell = 0
      DO ix = -bounds(1), bounds(1)
         DO iy = -bounds(2), bounds(2)
            DO iz = -bounds(3), bounds(3)
               IF (occupied(ix, iy, iz) == 0) CYCLE
               icell = icell + 1
               cells(:, icell) = [ix, iy, iz]
               kpoints%cell_to_index(ix, iy, iz) = icell
            END DO
         END DO
      END DO
      DEALLOCATE (occupied)
      kpoints_ptr => kpoints
      ALLOCATE (matrices(ncells))
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, distribution)
      ALLOCATE (row_sizes(SIZE(sizes)), col_sizes(SIZE(sizes)))
      row_sizes = sizes
      col_sizes = sizes
      CALL dbcsr_create(matrices(1), "full RI two-center", distribution, dbcsr_type_no_symmetry, row_sizes, col_sizes)
      DEALLOCATE (row_sizes, col_sizes)
      DO icell = 2, ncells
         CALL dbcsr_create(matrices(icell), template=matrices(1))
      END DO
      CALL build_2c_integrals(matrices, 0.0_dp, qs_env, neighbors, bs_env%basis_set_RI, &
                              bs_env%basis_set_RI, potential, do_kpoints=.TRUE., &
                              ext_kpoints=kpoints_ptr, regularization_RI=regularization)
      DO icell = 1, ncells
         IF (use_reference_basis) THEN
            IF (.NOT. bs_env%auto_ri%enabled) THEN
               CPABORT('Reference RI lattice requires AUTO_RI')
            END IF
            IF (.NOT. ASSOCIATED(reference_matrix_struct)) THEN
               CALL cp_fm_struct_create(reference_matrix_struct, context=blacs, nrow_global=SUM(sizes), &
                                        ncol_global=SUM(sizes), para_env=para)
            END IF
            CALL cp_fm_create(fm_R(icell), reference_matrix_struct)
         ELSE
            CALL cp_fm_create(fm_R(icell), bs_env%fm_RI_RI%matrix_struct)
         END IF
         CALL copy_dbcsr_to_fm(matrices(icell), fm_R(icell))
         CALL dbcsr_release(matrices(icell))
      END DO
      CALL dbcsr_distribution_release(distribution)
      CALL release_neighbor_list_sets(neighbors)
      IF (ASSOCIATED(reference_matrix_struct)) CALL cp_fm_struct_release(reference_matrix_struct)
      DEALLOCATE (matrices, sizes, kpoints%cell_to_index)

      CALL timestop(handle)

   END SUBROUTINE build_ri_2c_lattice

! **************************************************************************************************
!> \brief Evaluate I_PQ(q) = sum_R exp(i q.R) I_PQ(R) collectively into an existing dense matrix.
!>
!> q is in reciprocal fractional coordinates, so q.R includes 2π. The real-space matrices
!> stay distributed; no global RI-by-RI-by-cell array is replicated on the ranks.
!> \param fm_R ...
!> \param cells ...
!> \param qpoint ...
!> \param matrix_q ...
! **************************************************************************************************
   SUBROUTINE ri_2c_lattice_to_kp(fm_R, cells, qpoint, matrix_q)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_R(:)
      INTEGER, INTENT(IN)                                :: cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: qpoint(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: matrix_q(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_2c_lattice_to_kp'

      COMPLEX(KIND=dp)                                   :: phase
      INTEGER                                            :: handle, icell, icol, irow, ncols, nrows
      INTEGER, POINTER                                   :: cols(:), rows(:)

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(fm_R) > 0)
      CPASSERT(SIZE(cells, 1) == 3)
      CPASSERT(SIZE(cells, 2) == SIZE(fm_R))
      CALL cp_fm_get_info(fm_R(1), nrow_local=nrows, ncol_local=ncols, row_indices=rows, col_indices=cols)
      matrix_q = z_zero
      DO icell = 1, SIZE(fm_R)
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(qpoint, REAL(cells(:, icell), dp)), dp))
         DO icol = 1, ncols
            DO irow = 1, nrows
               matrix_q(rows(irow), cols(icol)) = matrix_q(rows(irow), cols(icol)) + &
                                                  phase*fm_R(icell)%local_data(irow, icol)
            END DO
         END DO
      END DO
      CALL fm_R(1)%matrix_struct%para_env%sum(matrix_q)

      CALL timestop(handle)

   END SUBROUTINE ri_2c_lattice_to_kp

! **************************************************************************************************
!> \brief Integrate each auxiliary function for the Gamma charge-mode projection.
!> For a primitive x^(2m) exp(-α x²), I_0=sqrt(π/α), I_2m=(2m-1)I_(2m-2)/(2α).
!> Odd Cartesian moments vanish; the existing spherical contraction supplies normalization.
!> \param qs_env Atomic auxiliary bases.
!> \param charge Integrals of the reference auxiliary functions.
! **************************************************************************************************
   SUBROUTINE gamma_auxiliary_charges(qs_env, charge)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(KIND=dp), INTENT(OUT)                      :: charge(:)

      INTEGER                                            :: first, iatom, ico, il, ipgf, iset, isgf, &
                                                            lmax, ncart, offset
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE                         :: moments(:), primitive(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL get_qs_env(qs_env, particle_set=particles, qs_kind_set=qs_kinds)
      charge = z_zero
      first = 0
      DO iatom = 1, SIZE(particles)
         CALL get_qs_kind(qs_kinds(particles(iatom)%atomic_kind%kind_number), &
                          basis_set=basis, basis_type='RI_AUX')
         DO iset = 1, basis%nset
            lmax = basis%lmax(iset)
            ncart = ncoset(lmax)
            ALLOCATE (moments(0:lmax), primitive(ncart*basis%npgf(iset)))
            primitive = 0.0_dp
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               moments = 0.0_dp
               moments(0) = SQRT(pi/alpha)
               DO il = 2, lmax, 2
                  moments(il) = REAL(il - 1, dp)*moments(il - 2)/(2.0_dp*alpha)
               END DO
               DO ico = ncoset(basis%lmin(iset) - 1) + 1, ncart
                  primitive((ipgf - 1)*ncart + ico) = PRODUCT(moments(indco(:, ico)))
               END DO
            END DO
            offset = basis%first_sgf(1, iset)
            DO isgf = offset, offset + basis%nsgf_set(iset) - 1
               charge(first + isgf) = CMPLX(SUM(basis%sphi(1:SIZE(primitive), isgf)*primitive), 0.0_dp, dp)
            END DO
            DEALLOCATE (moments, primitive)
         END DO
         first = first + basis%nsgf
      END DO
      CPASSERT(first == SIZE(charge))

   END SUBROUTINE gamma_auxiliary_charges

! **************************************************************************************************
!> \brief Return the contiguous part of 1..n owned by one rank of a block distribution.
!> \param n Number of distributed items (grid rows or unordered grid pairs).
!> \param rank ...
!> \param n_rank ...
!> \param first First owned item.
!> \param n_local Number of owned items.
! **************************************************************************************************
   PURE SUBROUTINE block_range(n, rank, n_rank, first, n_local)
      INTEGER, INTENT(IN)                                :: n, rank, n_rank
      INTEGER, INTENT(OUT)                               :: first, n_local

      first = INT(INT(n, int_8)*rank/n_rank) + 1
      n_local = INT(INT(n, int_8)*(rank + 1)/n_rank) - first + 1

   END SUBROUTINE block_range

END MODULE gw_ri_rs_small_cell_full_kp_lattice
