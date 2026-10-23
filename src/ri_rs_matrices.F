!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR collocation, layouts and panel contractions for real-space resolution of the identity.
! **************************************************************************************************
MODULE ri_rs_matrices
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_data_size, &
        dbcsr_get_info, dbcsr_get_occupation, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              max_elements_per_block
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE gw_auto_ri_types,                ONLY: auto_ri_type
   USE gw_utils_dbcsr,                  ONLY: hadamard_product_inplace
   USE gw_utils_fm,                     ONLY: fm_contract_ABA
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              max_line_length
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_types,                     ONLY: ri_rs_env
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: atomic_basis_at_grid_point, setup_square_topology, &
             release_dbcsr_topology_and_matrices, reserve_blocks_within_radius, &
             print_matrix_occupation, &
             contract_grid_panels, mask_grid_blocks_near_panel, &
             panel_template_elems, panel_mem_estimate_GB, resolve_grid_panels, &
             build_ri_rs_bare_aux_interaction

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_matrices'

CONTAINS

! **************************************************************************************************
!> \brief Builds the grid-by-AO collocation matrix with bounded dense evaluation panels.
!>        The molecular path uses atom-aligned row blocks and minimum-image AO values.
!>        The Gamma path uses uniform row blocks and sums neighbouring periodic images.
!> \param ri_rs grid, borrowed geometry references and output collocation matrix
!> \param basis_set_AO AO basis sets, one per atomic kind
!> \param matrix_template template for the process grid and AO column block sizes
!> \param para_env parallel environment
!> \param eps_filter block screening threshold
!> \param msg_elem_limit maximum number of elements per DBCSR message
!> \param unit_nr output unit
!> \param periodic_images select neighbouring-image summation for Gamma
! **************************************************************************************************
   SUBROUTINE atomic_basis_at_grid_point(ri_rs, basis_set_AO, matrix_template, &
                                         para_env, eps_filter, msg_elem_limit, unit_nr, periodic_images)
      TYPE(ri_rs_env), INTENT(INOUT)                     :: ri_rs
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis_set_AO(:)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER(KIND=int_8), INTENT(IN)                    :: msg_elem_limit
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: periodic_images

      CHARACTER(LEN=*), PARAMETER :: routineN = 'atomic_basis_at_grid_point'
      INTEGER, PARAMETER                                 :: grid_panel_size = 4096

      INTEGER :: blocks_per_panel, bs_eff, c_size, first_block, handle, i, i_blk, ia, iatom, &
         ikind, last_block, n_grid, natom, npcol, nprow, num_grid_chunks, panel_size, panel_start, &
         r_end, r_start, remaining, row_first, row_last, run, safe_max
      INTEGER, ALLOCATABLE                               :: blk_row_start(:)
      INTEGER, POINTER                                   :: c_blk_sizes(:), col_dist(:), &
                                                            r_blk_sizes(:), row_dist(:)
      REAL(KIND=dp)                                      :: cutoff_squared
      REAL(KIND=dp), ALLOCATABLE                         :: atom_col_buffer(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist, template_dist
      TYPE(gto_basis_set_type), POINTER                  :: ao_basis

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(ri_rs%cell))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(ri_rs%particle_set))
      CALL dbcsr_get_info(matrix_template, distribution=template_dist, col_blk_size=c_blk_sizes)
      natom = SIZE(c_blk_sizes)
      n_grid = SIZE(ri_rs%grid_points, 2)
      CPASSERT(SIZE(ri_rs%particle_set) == natom)
      CPASSERT(n_grid == ri_rs%n_grid_points)
      CALL dbcsr_distribution_get(template_dist, nprows=nprow, npcols=npcol)

      bs_eff = max_elements_per_block
      IF (.NOT. periodic_images) THEN
         safe_max = INT(msg_elem_limit*MAX(MIN(nprow, npcol), 1)/n_grid/2)
         bs_eff = MAX(1, MIN(bs_eff, MAX(1, safe_max)))
      END IF
      IF (periodic_images) THEN
         num_grid_chunks = n_grid/bs_eff
         IF (MOD(n_grid, bs_eff) /= 0) num_grid_chunks = num_grid_chunks + 1
      ELSE
         num_grid_chunks = 0
         DO ia = 1, natom
            run = ri_rs%grid_atom_boundaries(ia + 1) - ri_rs%grid_atom_boundaries(ia)
            IF (run > 0) num_grid_chunks = num_grid_chunks + (run + bs_eff - 1)/bs_eff
         END DO
      END IF
      ALLOCATE (r_blk_sizes(num_grid_chunks), blk_row_start(num_grid_chunks))
      i_blk = 0
      r_start = 1
      DO ia = 1, MERGE(1, natom, periodic_images)
         IF (periodic_images) THEN
            remaining = n_grid
         ELSE
            remaining = ri_rs%grid_atom_boundaries(ia + 1) - ri_rs%grid_atom_boundaries(ia)
         END IF
         DO WHILE (remaining > 0)
            i_blk = i_blk + 1
            r_blk_sizes(i_blk) = MIN(bs_eff, remaining)
            blk_row_start(i_blk) = r_start
            r_start = r_start + r_blk_sizes(i_blk)
            remaining = remaining - r_blk_sizes(i_blk)
         END DO
      END DO

      IF (.NOT. periodic_images) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,T71,I12)') 'RI-RS grid row-blocks of ϕ_μ(r_l)', num_grid_chunks
            WRITE (unit_nr, '(T2,A,T69,I12)') 'RI-RS grid points per block (max)', bs_eff
         END IF
         IF (ri_rs%cutoff_radius_v_w > 0.0_dp .OR. ri_rs%cutoff_radius_w0 > 0.0_dp) THEN
            IF (ALLOCATED(ri_rs%chunk_centroids)) DEALLOCATE (ri_rs%chunk_centroids)
            ALLOCATE (ri_rs%chunk_centroids(3, num_grid_chunks))
            DO i_blk = 1, num_grid_chunks
               r_start = blk_row_start(i_blk)
               r_end = r_start + r_blk_sizes(i_blk) - 1
               DO i = 1, 3
                  ri_rs%chunk_centroids(i, i_blk) = &
                     SUM(ri_rs%grid_points(i, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
               END DO
            END DO
         END IF
      END IF

      ALLOCATE (row_dist(num_grid_chunks), col_dist(natom))
      DO i = 1, num_grid_chunks
         row_dist(i) = MOD(i - 1, nprow)
      END DO
      DO i = 1, natom
         col_dist(i) = MOD(i - 1, npcol)
      END DO
      CALL dbcsr_distribution_new(dist, template=template_dist, row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(ri_rs%mat_phi_mu_l, name='phi_val_sparse', dist=dist, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=r_blk_sizes, col_blk_size=c_blk_sizes)

      blocks_per_panel = MAX(1, grid_panel_size/MAXVAL(r_blk_sizes))
      DO iatom = para_env%mepos + 1, natom, para_env%num_pe
         c_size = c_blk_sizes(iatom)
         ALLOCATE (atom_col_buffer(MIN(n_grid, blocks_per_panel*MAXVAL(r_blk_sizes)), c_size))
         cutoff_squared = ri_rs%radius_ao_per_atom(iatom)**2
         IF (.NOT. periodic_images .AND. ri_rs%cutoff_radius_ri_ao > 0.0_dp) THEN
            cutoff_squared = ri_rs%cutoff_radius_ri_ao**2
         END IF
         ikind = ri_rs%particle_set(iatom)%atomic_kind%kind_number
         ao_basis => basis_set_AO(ikind)%gto_basis_set
         DO first_block = 1, num_grid_chunks, blocks_per_panel
            last_block = MIN(first_block + blocks_per_panel - 1, num_grid_chunks)
            panel_start = blk_row_start(first_block)
            panel_size = blk_row_start(last_block) + r_blk_sizes(last_block) - panel_start
            atom_col_buffer(:panel_size, :) = 0.0_dp
            CALL evaluate_ao_basis_on_points(atom_col_buffer(:panel_size, :), &
                                             ri_rs%grid_points(:, panel_start:panel_start + panel_size - 1), &
                                             ao_basis, ri_rs%particle_set(iatom)%r, ri_rs%cell, &
                                             cutoff_squared=cutoff_squared, periodic_images=periodic_images)
            DO i_blk = first_block, last_block
               row_first = blk_row_start(i_blk) - panel_start + 1
               row_last = row_first + r_blk_sizes(i_blk) - 1
               IF (MAXVAL(ABS(atom_col_buffer(row_first:row_last, :))) > eps_filter) THEN
                  CALL dbcsr_put_block(ri_rs%mat_phi_mu_l, row=i_blk, col=iatom, &
                                       block=atom_col_buffer(row_first:row_last, :))
               END IF
            END DO
         END DO
         DEALLOCATE (atom_col_buffer)
      END DO
      CALL dbcsr_finalize(ri_rs%mat_phi_mu_l)
      DEALLOCATE (r_blk_sizes, row_dist, col_dist)
      CALL dbcsr_distribution_release(dist)
      CALL timestop(handle)
   END SUBROUTINE atomic_basis_at_grid_point

! **************************************************************************************************
!> \brief Creates a square topology from a selected row or column layout.
!> \param matrix_template source matrix
!> \param dim_type ROW or COL
!> \param square_dist output distribution
!> \param blk_sizes borrowed block sizes
!> \param mapped_dist owned mapping in the other process-grid direction
! **************************************************************************************************
   SUBROUTINE setup_square_topology(matrix_template, dim_type, square_dist, blk_sizes, mapped_dist)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      CHARACTER(LEN=*), INTENT(IN)                       :: dim_type
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: square_dist
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: blk_sizes, mapped_dist

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_square_topology'

      INTEGER                                            :: handle, i, np, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      TYPE(dbcsr_distribution_type)                      :: dist_template

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_template, distribution=dist_template, &
                          row_blk_size=row_blk, col_blk_size=col_blk)
      CALL dbcsr_distribution_get(dist_template, row_dist=row_dist, col_dist=col_dist, &
                                  nprows=nprows, npcols=npcols)

      IF (TRIM(dim_type) == 'ROW') THEN
         ! Creates ROW x ROW (e.g., Grid x Grid from mat_phi_mu_l)
         blk_sizes => row_blk
         np = npcols
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=row_dist, col_dist=mapped_dist)

      ELSE IF (TRIM(dim_type) == 'COL') THEN
         ! Creates COL x COL (e.g., Aux x Aux from mat_Z_lP)
         blk_sizes => col_blk
         np = nprows
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=mapped_dist, col_dist=col_dist)
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_square_topology

! **************************************************************************************************
!> \brief Releases RI-RS topology mappings and temporary DBCSR matrices.
!> \param dist distribution to release
!> \param mapped_dist owned mapping to deallocate
!> \param m1 optional matrix
!> \param m2 optional matrix
!> \param m3 optional matrix
!> \param m4 optional matrix
! **************************************************************************************************
   SUBROUTINE release_dbcsr_topology_and_matrices(dist, mapped_dist, m1, m2, m3, m4)

      TYPE(dbcsr_distribution_type), INTENT(INOUT), &
         OPTIONAL                                        :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: mapped_dist
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: m1, m2, m3, m4

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_dbcsr_topology_and_matrices'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (PRESENT(dist)) CALL dbcsr_distribution_release(dist)
      IF (PRESENT(mapped_dist)) THEN
         IF (ASSOCIATED(mapped_dist)) THEN
            DEALLOCATE (mapped_dist)
            NULLIFY (mapped_dist)
         END IF
      END IF
      IF (PRESENT(m1)) CALL dbcsr_release(m1)
      IF (PRESENT(m2)) CALL dbcsr_release(m2)
      IF (PRESENT(m3)) CALL dbcsr_release(m3)
      IF (PRESENT(m4)) CALL dbcsr_release(m4)

      CALL timestop(handle)

   END SUBROUTINE release_dbcsr_topology_and_matrices

! **************************************************************************************************
!> \brief Reserves grid blocks with block-centre distances within a fixed radius.
!> \param matrix target matrix
!> \param centers block-centre coordinates
!> \param radius truncation radius
! **************************************************************************************************
   SUBROUTINE reserve_blocks_within_radius(matrix, centers, radius)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centers
      REAL(KIND=dp), INTENT(IN)                          :: radius

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reserve_blocks_within_radius'

      INTEGER                                            :: handle, i, j, my_pcol, my_prow, &
                                                            nblkcols, nblkrows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      REAL(KIND=dp)                                      :: radius2
      REAL(KIND=dp), ALLOCATABLE                         :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows, nblkcols_total=nblkcols, &
                          row_blk_size=row_blk, col_blk_size=col_blk, distribution=dist)
      CALL dbcsr_distribution_get(dist, row_dist=row_dist, col_dist=col_dist, &
                                  myprow=my_prow, mypcol=my_pcol)
      CPASSERT(nblkrows == SIZE(centers, 2))
      CPASSERT(nblkcols == SIZE(centers, 2))

      radius2 = radius**2
      ALLOCATE (zero_blk(MAXVAL(row_blk(1:nblkrows)), MAXVAL(col_blk(1:nblkcols))))
      zero_blk(:, :) = 0.0_dp

      DO i = 1, nblkrows
         IF (row_dist(i) /= my_prow) CYCLE
         DO j = 1, nblkcols
            IF (col_dist(j) /= my_pcol) CYCLE
            IF (DOT_PRODUCT(centers(:, i) - centers(:, j), centers(:, i) - centers(:, j)) <= radius2) THEN
               CALL dbcsr_put_block(matrix, i, j, zero_blk(1:row_blk(i), 1:col_blk(j)))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(matrix)

      DEALLOCATE (zero_blk)
      CALL timestop(handle)

   END SUBROUTINE reserve_blocks_within_radius

! **************************************************************************************************
!> \brief Prints the percentage of non-zero elements in a distributed RI-RS matrix.
!> \param matrix Distributed matrix whose occupation is reported
!> \param label Mathematical matrix label used in the output
!> \param para_env communicator shared by all matrix ranks
!> \param unit_nr output unit; nonpositive values suppress printing
!> \param msg_elem_limit maximum number of elements per DBCSR message
!> \param suffix Optional text appended to the matrix label
! **************************************************************************************************
   SUBROUTINE print_matrix_occupation(matrix, label, para_env, unit_nr, msg_elem_limit, suffix)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER(KIND=int_8), INTENT(IN)                    :: msg_elem_limit
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: suffix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_matrix_occupation'

      CHARACTER(LEN=32)                                  :: output_format
      CHARACTER(LEN=max_line_length)                     :: msg, output_label
      INTEGER                                            :: handle, i, unicode_shift
      REAL(KIND=dp)                                      :: frac_2p31, max_loc, occ

      CALL timeset(routineN, handle)

      occ = dbcsr_get_occupation(matrix)
      max_loc = REAL(dbcsr_get_data_size(matrix), dp)
      CALL para_env%max(max_loc)

      IF (unit_nr > 0) THEN
         frac_2p31 = max_loc/REAL(msg_elem_limit, dp)
         output_label = 'Percentage of non-zero matrix elements in '//TRIM(label)
         IF (PRESENT(suffix)) output_label = TRIM(output_label)//TRIM(suffix)
         ! Fortran counts UTF-8 bytes, whereas the terminal displays each Greek letter in one
         ! column. Shift the absolute output tab once for every continuation byte.
         unicode_shift = 0
         DO i = 1, LEN_TRIM(output_label)
            IF (IAND(IACHAR(output_label(i:i)), 192) == 128) unicode_shift = unicode_shift + 1
         END DO
         WRITE (output_format, '(A,I0,A)') '(T2,A,T', 72 + unicode_shift, ',F7.2,A)'
         WRITE (unit_nr, output_format) TRIM(output_label), occ*100.0_dp, ' %'
         IF (frac_2p31 > 0.5_dp) THEN
            WRITE (msg, '(3A,F0.2,A)') &
               "The largest per-rank message of ", TRIM(label), " reaches ", frac_2p31, &
               " of the 32-bit limit that DBCSR uses for its message length. Beyond it the "// &
               "length overflows and multiply_cannon fails. Reduce the per-rank block size, "// &
               "for instance with more MPI ranks."
            CPWARN(TRIM(msg))
         END IF
         CALL m_flush(unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_matrix_occupation

! **************************************************************************************************
!> \brief Marks the grid blocks whose centroid lies within cutoff of the bounding box of the
!>        panel [blk0, blk1]'s chunk centroids.
!> \param centroids ...
!> \param blk0 ...
!> \param blk1 ...
!> \param cutoff ...
!> \param used ...
! **************************************************************************************************
   SUBROUTINE mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      INTEGER, INTENT(IN)                                :: blk0, blk1
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mask_grid_blocks_near_panel'

      INTEGER                                            :: c, handle, k
      REAL(KIND=dp)                                      :: cutoff2, d2, dx
      REAL(KIND=dp), DIMENSION(3)                        :: hi, lo

      CALL timeset(routineN, handle)

      cutoff2 = cutoff**2
      lo(:) = MINVAL(centroids(:, blk0:blk1), DIM=2)
      hi(:) = MAXVAL(centroids(:, blk0:blk1), DIM=2)

      ALLOCATE (used(SIZE(centroids, 2)))
      DO c = 1, SIZE(centroids, 2)
         d2 = 0.0_dp
         DO k = 1, 3
            dx = MAX(0.0_dp, lo(k) - centroids(k, c), centroids(k, c) - hi(k))
            d2 = d2 + dx*dx
         END DO
         used(c) = (d2 <= cutoff2)
      END DO

      CALL timestop(handle)

   END SUBROUTINE mask_grid_blocks_near_panel

! **************************************************************************************************
!> \brief Exact allocated-element count of the geo template of panel [blk0, blk1]: the very same
!>        per-block-pair centroid test as build_geo_template_panel, so this is the true DBCSR
!>        data size of A_pan/B_pan (DBCSR stores whole blocks).
!> \param r_blk_sizes ...
!> \param centroids ...
!> \param used ...
!> \param blk0 ...
!> \param blk1 ...
!> \param cutoff ...
!> \param nze_tmpl ...
! **************************************************************************************************
   SUBROUTINE panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: used
      INTEGER, INTENT(IN)                                :: blk0, blk1
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER(KIND=int_8), INTENT(OUT)                   :: nze_tmpl

      CHARACTER(LEN=*), PARAMETER :: routineN = 'panel_template_elems'

      INTEGER                                            :: c, handle, ib, n_used
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: used_idx
      REAL(KIND=dp)                                      :: cutoff2

      CALL timeset(routineN, handle)

      ! Compress the near mask once so the pair loop only visits candidate columns.
      n_used = COUNT(used)
      ALLOCATE (used_idx(n_used))
      n_used = 0
      DO c = 1, SIZE(used)
         IF (used(c)) THEN
            n_used = n_used + 1
            used_idx(n_used) = c
         END IF
      END DO

      cutoff2 = cutoff**2
      nze_tmpl = 0_int_8
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(blk0, blk1, n_used, used_idx, centroids, cutoff2, &
      !$OMP r_blk_sizes) PRIVATE(ib, c) REDUCTION(+:nze_tmpl)
      DO ib = blk0, blk1
         DO c = 1, n_used
            IF (SUM((centroids(:, ib) - centroids(:, used_idx(c)))**2) <= cutoff2) THEN
               nze_tmpl = nze_tmpl + INT(r_blk_sizes(ib), int_8)*INT(r_blk_sizes(used_idx(c)), int_8)
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE panel_template_elems

! **************************************************************************************************
!> \brief Per-rank peak memory (GB) of one panel step of the neighborhood-restricted
!>        contractions: up to three grid x grid panels (two A operands and the shared B operand)
!>        plus the grid x RI intermediates (tmp2 and the accumulation operand) and the
!>        grid x AO intermediate (tmpA), whose column support is the panel's geometric
!>        neighborhood fraction f_near = width/n_grid. Used by planning and memory reporting.
!> \param nze_tmpl ...
!> \param pan_rows ...
!> \param width ...
!> \param n_grid_total ...
!> \param n_RI ...
!> \param n_ao ...
!> \param n_procs ...
!> \param mem_GB ...
! **************************************************************************************************
   SUBROUTINE panel_mem_estimate_GB(nze_tmpl, pan_rows, width, n_grid_total, n_RI, n_ao, &
                                    n_procs, mem_GB)

      INTEGER(KIND=int_8), INTENT(IN)                    :: nze_tmpl
      INTEGER, INTENT(IN)                                :: pan_rows, width, n_grid_total, n_RI, &
                                                            n_ao, n_procs
      REAL(KIND=dp), INTENT(OUT)                         :: mem_GB

      REAL(KIND=dp)                                      :: f_near

      f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
      mem_GB = (3.0_dp*REAL(nze_tmpl, dp) + &
                REAL(pan_rows, dp)*f_near*(2.0_dp*REAL(n_RI, dp) + REAL(n_ao, dp)))* &
               8.0_dp/REAL(MAX(n_procs, 1), dp)*1.0E-9_dp

   END SUBROUTINE panel_mem_estimate_GB

! **************************************************************************************************
!> \brief Plans whole-grid-block panels using a message-size target without a cutoff.
!>        A positive cutoff additionally checks estimated message and memory costs.
!>        A single indivisible cutoff block exceeding an estimated limit is reported
!>        and accepted, since the panel cannot be split further.
!> \param msg_elem_limit maximum number of elements per DBCSR message
!> \param r_blk_sizes grid row-block sizes
!> \param min_dim smaller dimension of the DBCSR process grid
!> \param pan_first first grid block of each panel
!> \param pan_last last grid block of each panel
!> \param n_RI auxiliary basis dimension
!> \param n_ao AO basis dimension
!> \param n_procs number of MPI processes
!> \param mem_budget_GB per-rank panel memory budget; nonpositive disables the memory check
!> \param centroids optional grid-block centres
!> \param cutoff optional geometric truncation radius
! **************************************************************************************************
   SUBROUTINE plan_grid_panels(msg_elem_limit, r_blk_sizes, min_dim, pan_first, pan_last, &
                               n_RI, n_ao, n_procs, mem_budget_GB, centroids, cutoff)

      INTEGER(KIND=int_8), INTENT(IN)                    :: msg_elem_limit
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      INTEGER, INTENT(IN)                                :: min_dim
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last
      INTEGER, INTENT(IN)                                :: n_RI, n_ao, n_procs
      REAL(KIND=dp), INTENT(IN)                          :: mem_budget_GB
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'plan_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, ib, n_grid_blocks, &
                                                            n_grid_total, n_panels, panel_size, &
                                                            rows_acc, TARGET, width
      INTEGER(KIND=int_8)                                :: msg, msg_limit, nze_tmpl, side
      INTEGER, ALLOCATABLE                               :: tmp_first(:), tmp_last(:)
      LOGICAL                                            :: fits, use_cutoff, warned_memory, &
                                                            warned_message
      LOGICAL, ALLOCATABLE                               :: used(:)
      REAL(KIND=dp)                                      :: f_near, mem_GB

      CALL timeset(routineN, handle)

      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp
      msg_limit = msg_elem_limit/4
      warned_memory = .FALSE.
      warned_message = .FALSE.
      n_grid_blocks = SIZE(r_blk_sizes)
      n_grid_total = SUM(r_blk_sizes)
      panel_size = n_grid_total
      IF (.NOT. use_cutoff) THEN
         ! Preserve the original message-derived target and whole-block growth.
         panel_size = INT(0.5_dp*REAL(msg_elem_limit, dp)*REAL(min_dim, dp)/REAL(n_grid_total, dp))
         panel_size = MAX(1, MIN(panel_size, n_grid_total))
      END IF
      ALLOCATE (tmp_first(n_grid_blocks), tmp_last(n_grid_blocks))

      n_panels = 0
      blk0 = 1
      DO WHILE (blk0 <= n_grid_blocks)
         TARGET = panel_size
         DO
            rows_acc = 0
            DO ib = blk0, n_grid_blocks
               rows_acc = rows_acc + r_blk_sizes(ib)
               blk1 = ib
               IF (rows_acc >= TARGET) EXIT
            END DO
            IF (.NOT. use_cutoff) EXIT
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)
            width = SUM(r_blk_sizes, MASK=used)
            CALL panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)
            f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
            side = INT(REAL(rows_acc, dp)*f_near*REAL(MAX(n_RI, n_ao), dp), int_8)
            msg = MAX(nze_tmpl, side)/INT(MAX(min_dim, 1), int_8)
            CALL panel_mem_estimate_GB(nze_tmpl, rows_acc, width, n_grid_total, &
                                       n_RI, n_ao, n_procs, mem_GB)
            fits = msg <= msg_limit
            IF (mem_budget_GB > 0.0_dp) fits = fits .AND. mem_GB <= mem_budget_GB
            IF (fits) EXIT
            IF (blk1 == blk0) THEN
               IF (msg > msg_limit .AND. .NOT. warned_message) THEN
                  CPWARN("A single RI-RS grid block exceeds the estimated panel message limit.")
                  warned_message = .TRUE.
               END IF
               IF (mem_budget_GB > 0.0_dp .AND. mem_GB > mem_budget_GB .AND. .NOT. warned_memory) THEN
                  CPWARN("A single RI-RS grid block exceeds the estimated panel memory budget.")
                  warned_memory = .TRUE.
               END IF
               EXIT
            END IF
            ! Preserve whole blocks and prefer large panels to reduce the number of multiplies.
            TARGET = MAX(1, MIN(TARGET, rows_acc)/2)
         END DO
         n_panels = n_panels + 1
         tmp_first(n_panels) = blk0
         tmp_last(n_panels) = blk1
         blk0 = blk1 + 1
      END DO

      ALLOCATE (pan_first, SOURCE=tmp_first(:n_panels))
      ALLOCATE (pan_last, SOURCE=tmp_last(:n_panels))
      CALL timestop(handle)

   END SUBROUTINE plan_grid_panels

! **************************************************************************************************
!> \brief Creates an empty (panel_chunks x neighborhood_chunks) DBCSR matrix with zero blocks
!>        pre-allocated only where |centroid(panel_row r) - centroid(column c)| <= cutoff.
!>        Used with retain_sparsity=.TRUE. in the subsequent dbcsr_multiply so distant blocks
!>        of the grid-basis panels (ϕ G ϕ^T, Z W Z^T, ...) are never computed at all.
!> \param L_pan ...
!> \param L_full ...
!> \param centroids ...
!> \param cutoff ...
!> \param blk0 ...
!> \param A_template ...
!> \param col_map ...
! **************************************************************************************************
   SUBROUTINE build_geo_template_panel(L_pan, L_full, centroids, cutoff, blk0, A_template, col_map)
      TYPE(dbcsr_type), INTENT(IN)                       :: L_pan, L_full
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER, INTENT(IN)                                :: blk0
      TYPE(dbcsr_type), INTENT(OUT)                      :: A_template
      INTEGER, DIMENSION(:), INTENT(IN)                  :: col_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_geo_template_panel'

      INTEGER                                            :: c, cg, cs, handle, my_pcol, my_prow, &
                                                            n_grid_blks, n_pan_blks, npcols, &
                                                            nprows, r, rs
      INTEGER, DIMENSION(:), POINTER                     :: grid_blk_sizes, pan_blk_sizes
      REAL(KIND=dp)                                      :: cutoff2
      REAL(KIND=dp), ALLOCATABLE                         :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      cutoff2 = cutoff**2
      CALL dbcsr_get_info(L_pan, nblkrows_total=n_pan_blks, row_blk_size=pan_blk_sizes)
      CALL dbcsr_get_info(L_full, nblkrows_total=n_grid_blks, row_blk_size=grid_blk_sizes)

      ! create_product_matrix assigns row r to process MOD(r-1,nprows) and
      ! col c to MOD(c-1,npcols), so we can determine local ownership analytically.
      CALL create_product_matrix(L_pan, L_full, 'T', A_template)
      CALL dbcsr_get_info(A_template, distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, npcols=npcols, &
                                  myprow=my_prow, mypcol=my_pcol)

      ALLOCATE (zero_blk(MAXVAL(pan_blk_sizes(1:n_pan_blks)), &
                         MAXVAL(grid_blk_sizes(1:n_grid_blks))))
      zero_blk(:, :) = 0.0_dp

      DO r = 1, n_pan_blks
         IF (MOD(r - 1, nprows) /= my_prow) CYCLE
         rs = pan_blk_sizes(r)
         DO c = 1, n_grid_blks
            IF (MOD(c - 1, npcols) /= my_pcol) CYCLE
            cg = col_map(c)
            IF ((centroids(1, blk0 + r - 1) - centroids(1, cg))**2 + &
                (centroids(2, blk0 + r - 1) - centroids(2, cg))**2 + &
                (centroids(3, blk0 + r - 1) - centroids(3, cg))**2 <= cutoff2) THEN
               cs = grid_blk_sizes(c)
               CALL dbcsr_put_block(A_template, r, c, zero_blk(1:rs, 1:cs))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(A_template)

      DEALLOCATE (zero_blk)
      CALL timestop(handle)

   END SUBROUTINE build_geo_template_panel

! **************************************************************************************************
!> \brief Slices a contiguous range of grid row-blocks [blk0, blk1] out of a (grid x n) DBCSR
!>        matrix into a new (P x n) panel matrix: iterate the source's local blocks, put the
!>        in-range ones into the panel with a remapped row-block index, then finalize. Row-block
!>        index i of the panel corresponds to source row-block blk0+i-1.
!> \param mat_full ...
!> \param blk0 ...
!> \param blk1 ...
!> \param mat_panel ...
! **************************************************************************************************
   SUBROUTINE extract_grid_panel(mat_full, blk0, blk1, mat_panel)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_full
      INTEGER, INTENT(IN)                                :: blk0, blk1
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_panel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'extract_grid_panel'

      INTEGER                                            :: handle, ib, jb, npb
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_full, col_dist_full, &
                                                            row_blk_full, row_blk_pan, &
                                                            row_dist_full, row_dist_pan
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_pan
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_full, distribution=dist_full, &
                          row_blk_size=row_blk_full, col_blk_size=col_blk_full)
      CALL dbcsr_distribution_get(dist_full, row_dist=row_dist_full, col_dist=col_dist_full)

      npb = blk1 - blk0 + 1
      ALLOCATE (row_dist_pan(npb), row_blk_pan(npb))
      row_dist_pan(:) = row_dist_full(blk0:blk1)
      row_blk_pan(:) = row_blk_full(blk0:blk1)

      CALL dbcsr_distribution_new(dist_pan, template=dist_full, &
                                  row_dist=row_dist_pan, col_dist=col_dist_full)
      CALL dbcsr_create(mat_panel, name="grid_panel", dist=dist_pan, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_pan, col_blk_size=col_blk_full)

      CALL dbcsr_iterator_start(iter, mat_full)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         IF (ib < blk0 .OR. ib > blk1) CYCLE
         CALL dbcsr_put_block(mat_panel, ib - blk0 + 1, jb, blk)
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(mat_panel)

      CALL dbcsr_distribution_release(dist_pan)
      DEALLOCATE (row_dist_pan, row_blk_pan)

      CALL timestop(handle)

   END SUBROUTINE extract_grid_panel

! **************************************************************************************************
!> \brief Marks which column blocks of a DBCSR matrix carry at least one non-zero block anywhere
!>        (global union). Used to restrict the inner index of the panel multiplies to the
!>        AO/RI atoms that actually touch the panel (exact: dropped rows only meet zeros).
!> \param matrix ...
!> \param para_env ...
!> \param used ...
! **************************************************************************************************
   SUBROUTINE collect_used_col_blocks(matrix, para_env, used)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_used_col_blocks'

      INTEGER                                            :: handle, ib, jb, nblkcols
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iused
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nblkcols_total=nblkcols)
      ALLOCATE (iused(nblkcols))
      iused(:) = 0

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         iused(jb) = 1
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL para_env%sum(iused)

      ALLOCATE (used(nblkcols))
      used(:) = (iused(:) > 0)
      DEALLOCATE (iused)

      CALL timestop(handle)

   END SUBROUTINE collect_used_col_blocks

! **************************************************************************************************
!> \brief Copies the flagged block rows (compress_rows=.TRUE.) or block columns (.FALSE.) of a
!>        DBCSR matrix into a compressed matrix. The subset keeps the parent's process assignment
!>        along the compressed dimension, so every block stays on its owning rank: the extraction
!>        is purely local (zero communication), like extract_grid_panel.
!> \param mat_full ...
!> \param used ...
!> \param mat_out ...
!> \param compress_rows ...
!> \param blk_map ...
! **************************************************************************************************
   SUBROUTINE extract_masked_blocks(mat_full, used, mat_out, compress_rows, blk_map)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_full
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: used
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_out
      LOGICAL, INTENT(IN)                                :: compress_rows
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: blk_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'extract_masked_blocks'

      INTEGER                                            :: handle, ib, jb, n_blk, n_sub, r
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inv_map
      INTEGER, DIMENSION(:), POINTER                     :: blk_full, blk_sub, col_blk_full, &
                                                            col_dist_full, dist_full_1d, &
                                                            dist_sub_1d, row_blk_full, &
                                                            row_dist_full
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_sub
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_full, distribution=dist_full, &
                          row_blk_size=row_blk_full, col_blk_size=col_blk_full)
      CALL dbcsr_distribution_get(dist_full, row_dist=row_dist_full, col_dist=col_dist_full)

      IF (compress_rows) THEN
         blk_full => row_blk_full
         dist_full_1d => row_dist_full
      ELSE
         blk_full => col_blk_full
         dist_full_1d => col_dist_full
      END IF
      n_blk = SIZE(blk_full)

      n_sub = COUNT(used)
      CPASSERT(n_sub > 0)
      ALLOCATE (inv_map(n_blk), blk_sub(n_sub), dist_sub_1d(n_sub))
      IF (PRESENT(blk_map)) ALLOCATE (blk_map(n_sub))
      inv_map(:) = 0
      r = 0
      DO ib = 1, n_blk
         IF (used(ib)) THEN
            r = r + 1
            inv_map(ib) = r
            blk_sub(r) = blk_full(ib)
            dist_sub_1d(r) = dist_full_1d(ib)
            IF (PRESENT(blk_map)) blk_map(r) = ib
         END IF
      END DO

      IF (compress_rows) THEN
         CALL dbcsr_distribution_new(dist_sub, template=dist_full, &
                                     row_dist=dist_sub_1d, col_dist=col_dist_full)
         CALL dbcsr_create(mat_out, name="row_subset", dist=dist_sub, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=blk_sub, col_blk_size=col_blk_full)
      ELSE
         CALL dbcsr_distribution_new(dist_sub, template=dist_full, &
                                     row_dist=row_dist_full, col_dist=dist_sub_1d)
         CALL dbcsr_create(mat_out, name="col_subset", dist=dist_sub, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_full, col_blk_size=blk_sub)
      END IF

      CALL dbcsr_iterator_start(iter, mat_full)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         IF (compress_rows) THEN
            IF (inv_map(ib) > 0) CALL dbcsr_put_block(mat_out, inv_map(ib), jb, blk)
         ELSE
            IF (inv_map(jb) > 0) CALL dbcsr_put_block(mat_out, ib, inv_map(jb), blk)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(mat_out)

      CALL dbcsr_distribution_release(dist_sub)
      DEALLOCATE (inv_map, blk_sub, dist_sub_1d)

      CALL timestop(handle)

   END SUBROUTINE extract_masked_blocks

! **************************************************************************************************
!> \brief Creates the (empty) result matrix of op(mat_left) * op(mat_right) with the correct block
!>        structure and a distribution on the shared process grid, ready to be filled by
!>        dbcsr_multiply. Row structure comes from op(left), column structure from op(right).
!> \param mat_left ...
!> \param mat_right ...
!> \param transb 'N' or 'T' applied to mat_right
!> \param mat_out ...
! **************************************************************************************************
   SUBROUTINE create_product_matrix(mat_left, mat_right, transb, mat_out)

      TYPE(dbcsr_type), INTENT(IN)                       :: mat_left, mat_right
      CHARACTER(LEN=1), INTENT(IN)                       :: transb
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_out

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_product_matrix'

      INTEGER                                            :: handle, i, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_r, out_col_blk, out_col_dist, &
                                                            out_row_blk, out_row_dist, row_blk_l, &
                                                            row_blk_r
      TYPE(dbcsr_distribution_type)                      :: dist_l, dist_out

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_left, distribution=dist_l, row_blk_size=row_blk_l)
      CALL dbcsr_get_info(mat_right, row_blk_size=row_blk_r, col_blk_size=col_blk_r)
      CALL dbcsr_distribution_get(dist_l, nprows=nprows, npcols=npcols)

      ! block SIZES follow op(left)/op(right); DISTRIBUTIONS are freshly round-robined onto the
      ! shared process grid (a transposed operand's row-dist is NOT a valid col-dist on a
      ! non-square grid). dbcsr_multiply redistributes internally, so any valid mapping works.
      out_row_blk => row_blk_l
      IF (transb == 'N') THEN
         out_col_blk => col_blk_r
      ELSE
         out_col_blk => row_blk_r
      END IF

      ALLOCATE (out_row_dist(SIZE(out_row_blk)), out_col_dist(SIZE(out_col_blk)))
      DO i = 1, SIZE(out_row_blk)
         out_row_dist(i) = MOD(i - 1, nprows)
      END DO
      DO i = 1, SIZE(out_col_blk)
         out_col_dist(i) = MOD(i - 1, npcols)
      END DO

      CALL dbcsr_distribution_new(dist_out, template=dist_l, &
                                  row_dist=out_row_dist, col_dist=out_col_dist)
      CALL dbcsr_create(mat_out, name="panel_product", dist=dist_out, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=out_row_blk, col_blk_size=out_col_blk)
      CALL dbcsr_distribution_release(dist_out)
      DEALLOCATE (out_row_dist, out_col_dist)

      CALL timestop(handle)

   END SUBROUTINE create_product_matrix

! **************************************************************************************************
!> \brief Projects two operators onto grid panels, takes their Hadamard product, and
!>        back-projects it: out += scale * L_out^T (A_grid ∘ B_grid) L_out,
!>        where A_grid = L_A M_A L_A^T and B_grid = L_B M_B L_B^T.
!>        Each L has shape (grid, basis), and its M has shape (basis, basis).
!>        The three L matrices must share the ordered grid rows, row blocking and process grid.
!>        Only panel x grid slabs are formed; contributions are accumulated over all panels.
!>        The sequence of EPS_FILTER screening and floating-point reductions may differ from
!>        a full-grid evaluation.
!>
!>        The physical contractions and their signs/spin factors are specified in gw_ri_rs_kernels:
!>          χ:   L_A=L_B=ϕ, M_A=G_occ, M_B=G_vir, L_out=Z      -> RI x RI
!>          Σx:  L_A=ϕ, M_A=D, L_B=Z, M_B=M^-1 V^tr M^-1, L_out=ϕ -> AO x AO
!>          Σc:  L_A=ϕ, M_A=[G_occ,G_vir], L_B=Z, M_B=W, L_out=ϕ -> two AO x AO outputs.
!>        The second Σc output reuses the same W grid panel; its projection and back-projection
!>        otherwise follow the same operations as the first output.
!> \param L_A grid-by-basis collocation or fitting matrix for the first operand
!> \param M_A one or two square operators, with row/column blocking matching the columns of L_A
!> \param L_B grid-by-basis collocation or fitting matrix for the second operand
!> \param M_B second square operator, with row/column blocking matching the columns of L_B
!> \param L_out grid-by-output-basis matrix for back-projection
!> \param mat_out borrowed pointers to pre-created outputs, one per M_A, with the column blocking of L_out
!> \param scale prefactor of the Hadamard product, including the caller's sign/spin convention
!> \param eps DBCSR multiplication screening threshold
!> \param para_env communicator containing all matrix ranks
!> \param pan_first first global grid row-block of each ordered, nonoverlapping panel
!> \param pan_last last global grid row-block of each panel; panels together cover the grid
!> \param lb_eq_la L_B equals L_A in values and layout, allowing their grid/AO slices to be shared
!> \param lout_eq_la L_out equals L_A in values and layout, allowing back-projection slices to be shared
!> \param zero_out clear the output matrices before accumulating, otherwise add to their existing values
!> \param centroids optional global grid-block centres, paired with a positive cutoff for truncation
!> \param cutoff retain grid-block pairs within this centroid distance; nonpositive disables truncation
!> \param grid_occupation allocated-element fraction of the first grid operand before the Hadamard product
! **************************************************************************************************
   SUBROUTINE contract_grid_panels(L_A, M_A, L_B, M_B, L_out, mat_out, scale, eps, para_env, &
                                   pan_first, pan_last, lb_eq_la, lout_eq_la, zero_out, &
                                   centroids, cutoff, grid_occupation)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_A
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: M_A
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_B
      TYPE(dbcsr_type), INTENT(INOUT)                    :: M_B
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_out
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: mat_out
      REAL(KIND=dp), INTENT(IN)                          :: scale, eps
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pan_first, pan_last
      LOGICAL, INTENT(IN)                                :: lb_eq_la, lout_eq_la, zero_out
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: grid_occupation

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, i_A, ipan, &
                                                            n_grid_total, ncols_pan, nrows_pan
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap
      LOGICAL                                            :: use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, usedA, usedB
      TYPE(dbcsr_type)                                   :: A_pan(2), B_pan, LB_pan, LB_panC, &
                                                            MA_sub, MB_sub, tmp2, tmpA, tmpB
      TYPE(dbcsr_type), POINTER                          :: panel_out, RB_A, RB_B, RB_out
      TYPE(dbcsr_type), TARGET                           :: LA_near, LA_pan, LA_panC, LB_near, &
                                                            Lout_near, Lout_pan

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(M_A) >= 1 .AND. SIZE(M_A) <= 2)
      CPASSERT(SIZE(mat_out) == SIZE(M_A))

      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp
      IF (PRESENT(grid_occupation)) grid_occupation = 0.0_dp

      CALL dbcsr_get_info(L_A, nfullrows_total=n_grid_total)

      IF (zero_out) THEN
         DO i_A = 1, SIZE(mat_out)
            CALL dbcsr_set(mat_out(i_A)%matrix, 0.0_dp)
         END DO
      END IF

      DO ipan = 1, SIZE(pan_first)
         blk0 = pan_first(ipan)
         blk1 = pan_last(ipan)

         ! Slice only the current grid rows; compact away AO/RI columns that are exactly absent.
         CALL extract_grid_panel(L_A, blk0, blk1, LA_pan)
         CALL collect_used_col_blocks(LA_pan, para_env, usedA)
         IF (.NOT. lb_eq_la) THEN
            CALL extract_grid_panel(L_B, blk0, blk1, LB_pan)
            CALL collect_used_col_blocks(LB_pan, para_env, usedB)
         ELSE
            IF (ALLOCATED(usedB)) DEALLOCATE (usedB)
            ALLOCATE (usedB, SOURCE=usedA)
         END IF
         IF (.NOT. ANY(usedA) .OR. .NOT. ANY(usedB)) THEN
            CALL dbcsr_release(LA_pan)
            IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_pan)
            CYCLE
         END IF
         IF (lout_eq_la) THEN
            panel_out => LA_pan
         ELSE
            CALL extract_grid_panel(L_out, blk0, blk1, Lout_pan)
            panel_out => Lout_pan
         END IF
         CALL extract_masked_blocks(LA_pan, usedA, LA_panC, compress_rows=.FALSE.)

         ! Grid rows within reach of the panel: with the CUTOFF_RADIUS_RL_W truncation only
         ! they can appear as columns of the panel products / inner rows of the L_out multiply,
         ! so the system-wide phi/Z right operands are cut down to this neighborhood slice
         ! (local extraction, zero communication; exact w.r.t. the geo template).
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, grid_used)
            CALL extract_masked_blocks(L_A, grid_used, LA_near, compress_rows=.TRUE., blk_map=gmap)
            IF (.NOT. lb_eq_la) CALL extract_masked_blocks(L_B, grid_used, LB_near, compress_rows=.TRUE.)
            IF (.NOT. lout_eq_la) CALL extract_masked_blocks(L_out, grid_used, Lout_near, compress_rows=.TRUE.)
            RB_A => LA_near
         ELSE
            RB_A => L_A
         END IF
         IF (lb_eq_la) THEN
            RB_B => RB_A
         ELSE IF (use_cutoff) THEN
            RB_B => LB_near
         ELSE
            RB_B => L_B
         END IF
         IF (lout_eq_la) THEN
            RB_out => RB_A
         ELSE IF (use_cutoff) THEN
            RB_out => Lout_near
         ELSE
            RB_out => L_out
         END IF

         DO i_A = 1, SIZE(M_A)
            ! A_pan(i_A) = LA_pan * M_A(i_A) * L_A^T   (P x grid_near).
            ! When cutoff is active, A_pan(i_A) is pre-seeded with only nearby blocks via
            ! build_geo_template_panel, and the multiply uses retain_sparsity to skip
            ! computing distant blocks entirely within the chosen grid cutoff.
            CALL extract_masked_blocks(M_A(i_A), usedA, MA_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LA_panC, MA_sub, 'N', tmpA)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LA_panC, MA_sub, 0.0_dp, tmpA, filter_eps=eps)
            CALL dbcsr_release(MA_sub)
            IF (use_cutoff) THEN
               CALL build_geo_template_panel(LA_pan, LA_near, centroids, cutoff, blk0, A_pan(i_A), &
                                             col_map=gmap)
            ELSE
               CALL create_product_matrix(tmpA, RB_A, 'T', A_pan(i_A))
            END IF
            CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_A, 0.0_dp, A_pan(i_A), &
                                filter_eps=eps, retain_sparsity=use_cutoff)
            CALL dbcsr_release(tmpA)

            ! Grid-basis occupation of A_pan(i_A) = ϕ G ϕ^T, accumulated over ALL panels into the
            ! occupation of the (never formed) full grid x grid object:
            !   sum_panels nnz(A_pan(i_A)) / n_grid^2,  with nnz = occ * pan_rows * pan_cols.
            ! Panel-independent by construction -- a single-panel sample would instead report the
            ! local neighbor count of whichever region happens to land in that panel.
            IF (PRESENT(grid_occupation) .AND. i_A == 1) THEN
               CALL dbcsr_get_info(A_pan(i_A), nfullrows_total=nrows_pan, nfullcols_total=ncols_pan)
               grid_occupation = grid_occupation + dbcsr_get_occupation(A_pan(i_A))* &
                                 REAL(ncols_pan, dp)*REAL(nrows_pan, dp)/ &
                                 (REAL(n_grid_total, dp)*REAL(n_grid_total, dp))
            END IF
         END DO

         ! B_pan = LB_pan * M_B * L_B^T   (P x grid_near); reuse the L_A slices when L_B == L_A.
         ! B_pan retains the union of the A panels' block structures so that
         ! retain_sparsity forces the final multiply to fill only those blocks (exact for ∘).
         IF (lb_eq_la) THEN
            CALL extract_masked_blocks(M_B, usedA, MB_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LA_panC, MB_sub, 'N', tmpB)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LA_panC, MB_sub, 0.0_dp, tmpB, filter_eps=eps)
         ELSE
            CALL extract_masked_blocks(LB_pan, usedB, LB_panC, compress_rows=.FALSE.)
            CALL extract_masked_blocks(M_B, usedB, MB_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LB_panC, MB_sub, 'N', tmpB)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LB_panC, MB_sub, 0.0_dp, tmpB, filter_eps=eps)
            CALL dbcsr_release(LB_panC)
         END IF
         CALL dbcsr_release(MB_sub)
         CALL dbcsr_create(B_pan, template=A_pan(1))
         CALL dbcsr_copy(B_pan, A_pan(1))
         IF (SIZE(M_A) == 2) CALL dbcsr_add(B_pan, A_pan(2), 1.0_dp, 1.0_dp)
         CALL dbcsr_set(B_pan, 0.0_dp)
         ! The output pattern is already fixed. Omitting redundant on-the-fly filtering also
         ! avoids overflowing DBCSR's single-precision screening norms for conditioned Z fits.
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_B, 0.0_dp, B_pan, &
                             retain_sparsity=.TRUE.)
         CALL dbcsr_release(tmpB)
         CALL dbcsr_release(LA_panC)

         DO i_A = 1, SIZE(M_A)
            ! Reuse A_pan storage for its Hadamard product with the shared B panel.
            CALL hadamard_product_inplace(A_pan(i_A), B_pan, scale)
            CALL create_product_matrix(A_pan(i_A), RB_out, 'N', tmp2)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, A_pan(i_A), RB_out, 0.0_dp, tmp2, filter_eps=eps)
            CALL dbcsr_release(A_pan(i_A))

            ! mat_out += L_out_pan^T * tmp2.
            CALL dbcsr_multiply('T', 'N', 1.0_dp, panel_out, tmp2, 1.0_dp, mat_out(i_A)%matrix, filter_eps=eps)
            CALL dbcsr_release(tmp2)
         END DO
         CALL dbcsr_release(B_pan)
         IF (.NOT. lout_eq_la) CALL dbcsr_release(Lout_pan)
         IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_pan)
         CALL dbcsr_release(LA_pan)
         IF (use_cutoff) THEN
            CALL dbcsr_release(LA_near)
            IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_near)
            IF (.NOT. lout_eq_la) CALL dbcsr_release(Lout_near)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE contract_grid_panels

! **************************************************************************************************
!> \brief Build the unfiltered auxiliary interaction M^-1 V^tr M^-1 in the Z-column layout.
!>        AUTO_RI uses its already transformed Coulomb matrix in the optimized auxiliary basis.
!>        Filtering and cache lifetime are owned by the caller.
!> \param qs_env electronic-structure environment for ordinary RI integrals
!> \param metric_inverse inverse auxiliary metric
!> \param trunc_coulomb Coulomb operator for ordinary RI integrals
!> \param auto_ri optimized auxiliary basis and its Coulomb matrix
!> \param mat_Z_lP template defining the auxiliary blocking and process grid
!> \param matrix_V_aux output bare interaction in the auxiliary basis
! **************************************************************************************************
   SUBROUTINE build_ri_rs_bare_aux_interaction(qs_env, metric_inverse, trunc_coulomb, auto_ri, &
                                               mat_Z_lP, matrix_V_aux)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: metric_inverse
      TYPE(libint_potential_type), INTENT(IN)            :: trunc_coulomb
      TYPE(auto_ri_type), INTENT(IN)                     :: auto_ri
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_Z_lP
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix_V_aux

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_ri_rs_bare_aux_interaction'

      INTEGER                                            :: handle
      INTEGER, POINTER                                   :: blk_aux(:), dist_row_aux(:)
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_Vtr_Gamma(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux

      CALL timeset(routineN, handle)
      IF (auto_ri%enabled) THEN
         CPASSERT(auto_ri%ready)
         ALLOCATE (fm_Vtr_Gamma(1, 1))
         CALL cp_fm_create(fm_Vtr_Gamma(1, 1), metric_inverse%matrix_struct)
         CALL cp_fm_to_fm(auto_ri%V_pq, fm_Vtr_Gamma(1, 1))
      ELSE
         CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, metric_inverse%matrix_struct, &
                                 metric_inverse%matrix_struct%nrow_global, trunc_coulomb)
      END IF
      CALL fm_contract_ABA(metric_inverse, fm_Vtr_Gamma(:, 1))
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)
      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)
      CALL cp_fm_release(fm_Vtr_Gamma)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)
      CALL timestop(handle)
   END SUBROUTINE build_ri_rs_bare_aux_interaction

! **************************************************************************************************
!> \brief Resolve and cache grid-panel boundaries for RI-RS streaming contractions.
!> \param ri_rs collocation and fitting matrices, grid geometry and owned panel cache
!> \param para_env parallel environment
!> \param msg_elem_limit maximum number of elements per DBCSR message
!> \param pan_first ...
!> \param pan_last ...
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
! **************************************************************************************************
   SUBROUTINE resolve_grid_panels(ri_rs, para_env, msg_elem_limit, pan_first, pan_last, cutoff)

      TYPE(ri_rs_env), INTENT(INOUT)                     :: ri_rs
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER(KIND=int_8), INTENT(IN)                    :: msg_elem_limit
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'resolve_grid_panels'

      INTEGER                                            :: handle, min_dim, n_ao, n_grid_total, &
                                                            n_RI, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: use_cutoff
      REAL(KIND=dp)                                      :: effective_cutoff, mem_avail_GB, &
                                                            mem_budget_GB
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      use_cutoff = cutoff > 0.0_dp .AND. &
                   ALLOCATED(ri_rs%chunk_centroids)
      effective_cutoff = 0.0_dp
      IF (use_cutoff) effective_cutoff = cutoff
      IF (ALLOCATED(ri_rs%pan_first)) THEN
         ! This is an exact cache key, not a tolerance-based geometric comparison.
         IF (ri_rs%panel_plan_cutoff == effective_cutoff) THEN
            ALLOCATE (pan_first, SOURCE=ri_rs%pan_first)
            ALLOCATE (pan_last, SOURCE=ri_rs%pan_last)
            CALL timestop(handle)
            RETURN
         END IF
         DEALLOCATE (ri_rs%pan_first, ri_rs%pan_last)
      END IF

      ! MIN(nprows, npcols) is the divisor that bounds the worst-rank Cannon message: a
      ! P x n_grid panel is replicated into block row strips (P/nprows x n_grid) or column
      ! strips (P x n_grid/npcols) during multiply_cannon, so the largest single-rank
      ! message is ~ P*n_grid / MIN(nprows,npcols) elements.
      CALL dbcsr_get_info(ri_rs%mat_phi_mu_l, nfullrows_total=n_grid_total, nfullcols_total=n_ao, &
                          row_blk_size=r_blk_sizes, distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, npcols=npcols)
      min_dim = MAX(MIN(nprows, npcols), 1)

      CALL dbcsr_get_info(ri_rs%mat_Z_lP, nfullcols_total=n_RI)
      ! Only the cutoff path uses the measured free memory to limit panel sizes.
      mem_budget_GB = 0.0_dp
      IF (use_cutoff) THEN
         CALL mp_mem_avail_per_rank_GB(para_env, mem_avail_GB)
         mem_budget_GB = 0.5_dp*mem_avail_GB
      END IF
      CALL plan_grid_panels(msg_elem_limit, r_blk_sizes, min_dim, pan_first, pan_last, &
                            n_RI, n_ao, para_env%num_pe, mem_budget_GB, &
                            centroids=ri_rs%chunk_centroids, cutoff=effective_cutoff)

      ALLOCATE (ri_rs%pan_first, SOURCE=pan_first)
      ALLOCATE (ri_rs%pan_last, SOURCE=pan_last)
      ri_rs%panel_plan_cutoff = effective_cutoff

      CALL timestop(handle)

   END SUBROUTINE resolve_grid_panels

END MODULE ri_rs_matrices
