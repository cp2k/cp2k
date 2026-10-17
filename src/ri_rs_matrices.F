!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR collocation matrices and layouts for real-space resolution of the identity.
! **************************************************************************************************
MODULE ri_rs_matrices
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_data_size, dbcsr_get_info, &
        dbcsr_get_occupation, dbcsr_put_block, dbcsr_release, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: max_elements_per_block
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              max_line_length
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_para_env_type
   USE ri_rs_types,                     ONLY: ri_rs_env
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: atomic_basis_at_grid_point, setup_square_topology, &
             release_dbcsr_topology_and_matrices, reserve_blocks_within_radius, print_matrix_occupation

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

END MODULE ri_rs_matrices
