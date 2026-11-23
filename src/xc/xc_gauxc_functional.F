!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#ifdef __GAUXC
#include "gauxc/gauxc_config.f"
#endif

MODULE xc_gauxc_functional
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE external_potential_types,        ONLY: all_potential_type,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: xc_vdw_fun_nonloc
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              has_nlcc,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE skala_model_selection,           ONLY: find_single_functional_subsection,&
                                              model_key_is_none,&
                                              resolve_skala_execution_target,&
                                              resolve_skala_model_path,&
                                              skala_execution_target_type
   USE string_utilities,                ONLY: uppercase
   USE xc_gauxc_cache,                  ONLY: cp_gauxc_cache_params,&
                                              cp_gauxc_cache_type,&
                                              dense_to_dbcsr_layout_type,&
                                              gauxc_cache_ensure_gradient,&
                                              gauxc_cache_init
   USE xc_gauxc_interface,              ONLY: &
        cp_gauxc_basisset_type, cp_gauxc_grid_type, cp_gauxc_integrator_type, &
        cp_gauxc_molecule_type, cp_gauxc_runtime_type, cp_gauxc_status_type, cp_gauxc_xc_type, &
        gauxc_check_status, gauxc_compute_xc, gauxc_compute_xc_gradient, gauxc_create_basisset, &
        gauxc_create_grid, gauxc_create_integrator, gauxc_create_molecule, &
        gauxc_create_shared_runtime, gauxc_destroy_basisset, gauxc_destroy_grid, &
        gauxc_destroy_integrator, gauxc_destroy_molecule, gauxc_destroy_shared_runtime, &
        gauxc_model_active, gauxc_write_basisset_hdf5, gauxc_write_molecule_hdf5
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: apply_gauxc, gauxc_gapw_all_all_electron, gauxc_gapw_all_pseudopotentials, &
             gauxc_gapw_has_paw_pseudopotentials, skala_info

   LOGICAL, SAVE, PRIVATE :: fd_mixed_space_warned = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Return the GAPW kind-set pseudopotential flags in one pass.
!> \param qs_kind_set ...
!> \param has_pseudopotentials ...
!> \param all_pseudopotentials ...
!> \param has_paw_pseudopotentials ...
!> \param all_all_electron ...
! **************************************************************************************************
   SUBROUTINE gauxc_gapw_kind_flags(qs_kind_set, has_pseudopotentials, &
                                    all_pseudopotentials, has_paw_pseudopotentials, &
                                    all_all_electron)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL, INTENT(OUT)                               :: has_pseudopotentials, &
                                                            all_pseudopotentials, &
                                                            has_paw_pseudopotentials, &
                                                            all_all_electron

      INTEGER                                            :: ikind
      LOGICAL                                            :: is_pseudopotential, paw_atom
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CPASSERT(ASSOCIATED(qs_kind_set))
      has_pseudopotentials = .FALSE.
      all_pseudopotentials = .TRUE.
      has_paw_pseudopotentials = .FALSE.
      all_all_electron = .TRUE.
      DO ikind = 1, SIZE(qs_kind_set)
         NULLIFY (all_potential, gth_potential, sgp_potential)
         CALL get_qs_kind(qs_kind_set(ikind), all_potential=all_potential, &
                          gth_potential=gth_potential, paw_atom=paw_atom, &
                          sgp_potential=sgp_potential)
         is_pseudopotential = ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)
         IF (is_pseudopotential) THEN
            has_pseudopotentials = .TRUE.
            IF (paw_atom) has_paw_pseudopotentials = .TRUE.
            all_all_electron = .FALSE.
         ELSE
            all_pseudopotentials = .FALSE.
         END IF
         IF (.NOT. ASSOCIATED(all_potential)) all_all_electron = .FALSE.
      END DO

   END SUBROUTINE gauxc_gapw_kind_flags

! **************************************************************************************************
!> \brief Return whether every GAPW kind uses an all-electron potential.
!> \param qs_kind_set ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_gapw_all_all_electron(qs_kind_set) RESULT(all_all_electron)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL                                            :: all_all_electron

      LOGICAL                                            :: all_pseudopotentials, &
                                                            has_paw_pseudopotentials, &
                                                            has_pseudopotentials

      CALL gauxc_gapw_kind_flags(qs_kind_set, has_pseudopotentials, all_pseudopotentials, &
                                 has_paw_pseudopotentials, all_all_electron)

   END FUNCTION gauxc_gapw_all_all_electron

! **************************************************************************************************
!> \brief ...
!> \param dense ...
!> \param block ...
!> \param transposed ...
!> \param symmetry_sign ...
!> \param has_symmetry ...
!> \param irow ...
!> \param icol ...
!> \param row_start ...
!> \param col_start ...
!> \param nrow ...
!> \param ncol ...
! **************************************************************************************************
   SUBROUTINE copy_block_to_dense(dense, block, transposed, symmetry_sign, has_symmetry, &
                                  irow, icol, row_start, col_start, nrow, ncol)
      REAL(c_double), DIMENSION(:, :), INTENT(INOUT)     :: dense
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: block
      LOGICAL, INTENT(IN)                                :: transposed
      REAL(KIND=dp), INTENT(IN)                          :: symmetry_sign
      LOGICAL, INTENT(IN)                                :: has_symmetry
      INTEGER, INTENT(IN)                                :: irow, icol, row_start, col_start, nrow, &
                                                            ncol

      INTEGER                                            :: col_end, row_end

      row_end = row_start + nrow - 1
      col_end = col_start + ncol - 1
      IF (transposed) THEN
         dense(row_start:row_end, col_start:col_end) = symmetry_sign*TRANSPOSE(block)
      ELSE
         dense(row_start:row_end, col_start:col_end) = block
      END IF
      IF (has_symmetry .AND. irow /= icol) THEN
         ! Write the mirror directly from the block instead of reading back the
         ! just-written slice: for a transposed block the mirror is sign**2*block
         ! (= block), otherwise it is sign*TRANSPOSE(block).
         IF (transposed) THEN
            dense(col_start:col_end, row_start:row_end) = block
         ELSE
            dense(col_start:col_end, row_start:row_end) = symmetry_sign*TRANSPOSE(block)
         END IF
      END IF

   END SUBROUTINE copy_block_to_dense

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param dense_mat ...
!> \param para_env ...
!> \param replicated ...
!> \param fill_replicated ...
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense(dbcsr_mat, dense_mat, para_env, replicated, fill_replicated)
      USE cp_dbcsr_api, ONLY: dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_get_info, &
                              dbcsr_get_matrix_type, dbcsr_iterator_blocks_left, &
                              dbcsr_iterator_next_block, dbcsr_iterator_readonly_start, &
                              dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_type_antisymmetric, &
                              dbcsr_type_symmetric
      TYPE(dbcsr_p_type), INTENT(IN)                     :: dbcsr_mat
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: dense_mat
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: replicated
      LOGICAL, INTENT(IN), OPTIONAL                      :: fill_replicated

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: col_start, icol, irow, izero, ncol, &
                                                            ncols, nrow, nrows, numnodes, row_start
      LOGICAL                                            :: fill_local, has_symmetry, transposed
      REAL(c_double), POINTER                            :: block(:, :)
      REAL(KIND=dp)                                      :: symmetry_sign
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(dbcsr_mat%matrix, &
                          nfullrows_total=nrows, &
                          nfullcols_total=ncols, &
                          distribution=dist)
      CALL dbcsr_distribution_get(dist, numnodes=numnodes)
      IF (PRESENT(replicated)) replicated = (numnodes == 1)
      fill_local = .FALSE.
      IF (PRESENT(fill_replicated)) fill_local = fill_replicated
      matrix_type = dbcsr_get_matrix_type(dbcsr_mat%matrix)
      has_symmetry = matrix_type == dbcsr_type_symmetric .OR. matrix_type == dbcsr_type_antisymmetric
      symmetry_sign = 1.0_dp
      IF (matrix_type == dbcsr_type_antisymmetric) symmetry_sign = -1.0_dp

      IF (.NOT. ALLOCATED(dense_mat)) THEN
         ALLOCATE (dense_mat(nrows, ncols))
      ELSE IF (.NOT. ALL(SHAPE(dense_mat) == [nrows, ncols])) THEN
         DEALLOCATE (dense_mat)
         ALLOCATE (dense_mat(nrows, ncols))
      ELSE
         CPASSERT(ALL(SHAPE(dense_mat) == [nrows, ncols]))
      END IF
      ! Replicated DBCSR blocks must enter a following MPI sum exactly once. If the
      ! caller requests fill_replicated, every rank fills its own copy locally so the
      ! collective can be skipped entirely.
      IF (numnodes == 1 .AND. para_env%num_pe > 1 .AND. para_env%mepos /= 0 .AND. .NOT. fill_local) THEN
         RETURN
      END IF

      ! Stored blocks and their symmetry mirrors occupy disjoint dense slices.
!$OMP PARALLEL DEFAULT(NONE) SHARED(dbcsr_mat, dense_mat, has_symmetry, symmetry_sign, ncols) &
!$OMP PRIVATE(iter, irow, icol, izero, block, transposed, row_start, col_start, nrow, ncol)
!$OMP DO SCHEDULE(STATIC)
      DO izero = 1, ncols
         dense_mat(:, izero) = 0._dp
      END DO
!$OMP END DO
!$OMP BARRIER
      CALL dbcsr_iterator_readonly_start(iter, dbcsr_mat%matrix, dynamic=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block, transposed=transposed, &
                                        row_size=nrow, col_size=ncol, row_offset=row_start, col_offset=col_start)
         CALL copy_block_to_dense(dense_mat, block, transposed, symmetry_sign, has_symmetry, &
                                  irow, icol, row_start, col_start, nrow, ncol)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

   END SUBROUTINE dbcsr_to_dense

! **************************************************************************************************
!> \brief Fill two dense matrices from two identically distributed DBCSR matrices in one iterator pass.
!> \param dbcsr_a ...
!> \param dense_a ...
!> \param dbcsr_b ...
!> \param dense_b ...
!> \param para_env ...
!> \param replicated ...
!> \param fill_replicated ...
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense_pair(dbcsr_a, dense_a, dbcsr_b, dense_b, para_env, replicated, &
                                  fill_replicated)
      USE cp_dbcsr_api, ONLY: dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_get_info, &
                              dbcsr_get_matrix_type, dbcsr_iterator_blocks_left, &
                              dbcsr_iterator_next_block, dbcsr_iterator_readonly_start, &
                              dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_type_antisymmetric, &
                              dbcsr_type_symmetric
      TYPE(dbcsr_p_type), INTENT(IN)                     :: dbcsr_a
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: dense_a
      TYPE(dbcsr_p_type), INTENT(IN)                     :: dbcsr_b
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: dense_b
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: replicated
      LOGICAL, INTENT(IN), OPTIONAL                      :: fill_replicated

      CHARACTER                                          :: matrix_type_a, matrix_type_b
      INTEGER :: col_start, col_start_b, icol, icol_b, irow, irow_b, izero, ncol, ncol_b, ncols, &
         nrow, nrow_b, nrows, numnodes, row_start, row_start_b
      LOGICAL                                            :: fill_local, has_symmetry_a, &
                                                            has_symmetry_b, transposed_a, &
                                                            transposed_b
      REAL(c_double), POINTER                            :: block_a(:, :), block_b(:, :)
      REAL(KIND=dp)                                      :: symmetry_sign_a, symmetry_sign_b
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_iterator_type)                          :: iter_a, iter_b

      CALL dbcsr_get_info(dbcsr_a%matrix, nfullrows_total=nrows, nfullcols_total=ncols, &
                          distribution=dist)
      CALL dbcsr_distribution_get(dist, numnodes=numnodes)
      IF (PRESENT(replicated)) replicated = (numnodes == 1)
      fill_local = .FALSE.
      IF (PRESENT(fill_replicated)) fill_local = fill_replicated
      matrix_type_a = dbcsr_get_matrix_type(dbcsr_a%matrix)
      matrix_type_b = dbcsr_get_matrix_type(dbcsr_b%matrix)
      has_symmetry_a = matrix_type_a == dbcsr_type_symmetric .OR. &
                       matrix_type_a == dbcsr_type_antisymmetric
      has_symmetry_b = matrix_type_b == dbcsr_type_symmetric .OR. &
                       matrix_type_b == dbcsr_type_antisymmetric
      symmetry_sign_a = MERGE(-1.0_dp, 1.0_dp, matrix_type_a == dbcsr_type_antisymmetric)
      symmetry_sign_b = MERGE(-1.0_dp, 1.0_dp, matrix_type_b == dbcsr_type_antisymmetric)

      IF (.NOT. ALLOCATED(dense_a)) THEN
         ALLOCATE (dense_a(nrows, ncols))
      ELSE IF (.NOT. ALL(SHAPE(dense_a) == [nrows, ncols])) THEN
         DEALLOCATE (dense_a)
         ALLOCATE (dense_a(nrows, ncols))
      END IF
      IF (.NOT. ALLOCATED(dense_b)) THEN
         ALLOCATE (dense_b(nrows, ncols))
      ELSE IF (.NOT. ALL(SHAPE(dense_b) == [nrows, ncols])) THEN
         DEALLOCATE (dense_b)
         ALLOCATE (dense_b(nrows, ncols))
      END IF
      IF (numnodes == 1 .AND. para_env%num_pe > 1 .AND. para_env%mepos /= 0 .AND. .NOT. fill_local) THEN
         RETURN
      END IF

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(dbcsr_a, dbcsr_b, dense_a, dense_b, has_symmetry_a, has_symmetry_b, ncols, &
!$OMP        symmetry_sign_a, symmetry_sign_b) &
!$OMP PRIVATE(iter_a, iter_b, irow, icol, irow_b, icol_b, izero, block_a, block_b, transposed_a, &
!$OMP         transposed_b, row_start, col_start, row_start_b, col_start_b, nrow, ncol, nrow_b, &
!$OMP         ncol_b)
!$OMP DO SCHEDULE(STATIC)
      DO izero = 1, ncols
         dense_a(:, izero) = 0._dp
         dense_b(:, izero) = 0._dp
      END DO
!$OMP END DO
!$OMP BARRIER
      CALL dbcsr_iterator_readonly_start(iter_a, dbcsr_a%matrix, dynamic=.TRUE.)
      CALL dbcsr_iterator_readonly_start(iter_b, dbcsr_b%matrix, dynamic=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter_a))
         CPASSERT(dbcsr_iterator_blocks_left(iter_b))
         CALL dbcsr_iterator_next_block(iter_a, irow, icol, block_a, transposed=transposed_a, &
                                        row_size=nrow, col_size=ncol, row_offset=row_start, &
                                        col_offset=col_start)
         CALL dbcsr_iterator_next_block(iter_b, irow_b, icol_b, block_b, transposed=transposed_b, &
                                        row_size=nrow_b, col_size=ncol_b, row_offset=row_start_b, &
                                        col_offset=col_start_b)
         CPASSERT(irow_b == irow .AND. icol_b == icol)
         CPASSERT(nrow_b == nrow .AND. ncol_b == ncol)
         CALL copy_block_to_dense(dense_a, block_a, transposed_a, symmetry_sign_a, has_symmetry_a, &
                                  irow, icol, row_start, col_start, nrow, ncol)
         CALL copy_block_to_dense(dense_b, block_b, transposed_b, symmetry_sign_b, has_symmetry_b, &
                                  irow, icol, row_start_b, col_start_b, nrow_b, ncol_b)
      END DO
      CPASSERT(.NOT. dbcsr_iterator_blocks_left(iter_b))
      CALL dbcsr_iterator_stop(iter_a)
      CALL dbcsr_iterator_stop(iter_b)
!$OMP END PARALLEL

   END SUBROUTINE dbcsr_to_dense_pair

! **************************************************************************************************
!> \brief Fill a DBCSR matrix from a dense symmetric matrix with full upper block structure.
!>        This creates all upper-triangular blocks, not just those present in a template,
!>        reusing the target matrix when its shape already matches.
!>        This is needed because GauXC computes VXC for the full dense density matrix.
!> \param layout ...
!> \param dense_mat Input dense matrix
!> \param template_dbcsr Template DBCSR matrix for distribution and block sizes
!> \param dbcsr_mat Output DBCSR matrix with full upper block structure
!> \param dense_zeta ...
!> \param z_sign ...
! **************************************************************************************************
   SUBROUTINE dense_to_dbcsr_into(layout, dense_mat, template_dbcsr, dbcsr_mat, dense_zeta, z_sign)
      USE cp_dbcsr_api, ONLY: &
         dbcsr_clear, &
         dbcsr_create, &
         dbcsr_distribution_get, &
         dbcsr_distribution_type, &
         dbcsr_finalize, &
         dbcsr_get_info, &
         dbcsr_get_num_blocks, &
         dbcsr_init_p, &
         dbcsr_iterator_blocks_left, &
         dbcsr_iterator_next_block, &
         dbcsr_iterator_start, &
         dbcsr_iterator_stop, &
         dbcsr_iterator_type, &
         dbcsr_put_block, &
         dbcsr_release, &
         dbcsr_type_symmetric, &
         dbcsr_work_create
      TYPE(dense_to_dbcsr_layout_type), INTENT(INOUT)    :: layout
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: dense_mat
      TYPE(dbcsr_p_type), INTENT(IN)                     :: template_dbcsr
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: dbcsr_mat
      REAL(c_double), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: dense_zeta
      REAL(c_double), INTENT(IN), OPTIONAL               :: z_sign

      INTEGER :: c1, c1t, c2, c2t, ib, icol, irow, mynode, mypcol, myprow, nblkcols_total, &
         nblkrows_total, nc, nc_t, ncols, nr, nr_t, nrows, r1, r1t, r2, r2t
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist
      LOGICAL                                            :: rebuild, transposed
      REAL(c_double)                                     :: zsign
      REAL(c_double), DIMENSION(:, :), POINTER           :: block_ptr
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(template_dbcsr%matrix, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nrows, &
                          nfullcols_total=ncols, &
                          distribution=dist)

      CPASSERT(nrows == SIZE(dense_mat, 1))
      CPASSERT(ncols == SIZE(dense_mat, 2))
      IF (PRESENT(dense_zeta)) THEN
         CPASSERT(SIZE(dense_zeta, 1) == nrows .AND. SIZE(dense_zeta, 2) == ncols)
      END IF

      CALL dbcsr_distribution_get(dist, mynode=mynode, myprow=myprow, mypcol=mypcol, &
                                  row_dist=row_dist, col_dist=col_dist)

      IF (.NOT. dense_layout_matches(layout, nrows, ncols, nblkrows_total, nblkcols_total, &
                                     mynode, myprow, mypcol, row_blk_size, col_blk_size, &
                                     row_dist, col_dist)) THEN
         CALL build_dense_layout(layout, nrows, ncols, nblkrows_total, nblkcols_total, &
                                 mynode, myprow, mypcol, row_blk_size, col_blk_size, &
                                 row_dist, col_dist)
      END IF

      zsign = 1.0_dp
      IF (PRESENT(z_sign)) zsign = z_sign

      ! Reuse the existing upper-triangular block pattern when the target matrix
      ! already matches the cached layout: only the block values need refreshing.
      rebuild = .TRUE.
      IF (ASSOCIATED(dbcsr_mat%matrix)) THEN
         IF (dbcsr_target_matches_layout(dbcsr_mat%matrix, layout)) THEN
            IF (dbcsr_get_num_blocks(dbcsr_mat%matrix) == SIZE(layout%blk_row)) rebuild = .FALSE.
         ELSE
            CALL dbcsr_release(dbcsr_mat%matrix)
            DEALLOCATE (dbcsr_mat%matrix)
         END IF
      END IF

      IF (rebuild) THEN
         IF (ASSOCIATED(dbcsr_mat%matrix)) THEN
            CALL dbcsr_clear(dbcsr_mat%matrix)
         ELSE
            CALL dbcsr_init_p(dbcsr_mat%matrix)
            CALL dbcsr_create(dbcsr_mat%matrix, &
                              template=template_dbcsr%matrix, &
                              name="VXC from GauXC (dense)", &
                              matrix_type=dbcsr_type_symmetric)
         END IF
         CALL dbcsr_work_create(dbcsr_mat%matrix, work_mutable=.TRUE.)
         DO ib = 1, SIZE(layout%blk_row)
            irow = layout%blk_row(ib)
            icol = layout%blk_col(ib)
            nr = layout%row_blk_size(irow)
            nc = layout%col_blk_size(icol)
            nr_t = layout%row_blk_size(icol)
            nc_t = layout%col_blk_size(irow)
            r1 = layout%r_offset(irow)
            r2 = r1 + nr - 1
            c1 = layout%c_offset(icol)
            c2 = c1 + nc - 1
            r1t = layout%r_offset(icol)
            r2t = r1t + nr_t - 1
            c1t = layout%c_offset(irow)
            c2t = c1t + nc_t - 1
            IF (PRESENT(dense_zeta)) THEN
               layout%blk(1:nr, 1:nc) = 0.5_dp*( &
                                        (dense_mat(r1:r2, c1:c2) + zsign*dense_zeta(r1:r2, c1:c2)) + &
                                        TRANSPOSE(dense_mat(r1t:r2t, c1t:c2t) + zsign*dense_zeta(r1t:r2t, c1t:c2t)))
            ELSE
               layout%blk(1:nr, 1:nc) = 0.5_dp*( &
                                        dense_mat(r1:r2, c1:c2) + TRANSPOSE(dense_mat(r1t:r2t, c1t:c2t)))
            END IF
            CALL dbcsr_put_block(dbcsr_mat%matrix, irow, icol, layout%blk(1:nr, 1:nc))
         END DO
         CALL dbcsr_finalize(dbcsr_mat%matrix)
      ELSE
         CALL dbcsr_iterator_start(iter, dbcsr_mat%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, block_ptr, transposed=transposed, &
                                           row_size=nr, col_size=nc, row_offset=r1, col_offset=c1)
            r2 = r1 + nr - 1
            c2 = c1 + nc - 1
            r1t = layout%r_offset(icol)
            r2t = r1t + layout%row_blk_size(icol) - 1
            c1t = layout%c_offset(irow)
            c2t = c1t + layout%col_blk_size(irow) - 1
            IF (transposed) THEN
               IF (PRESENT(dense_zeta)) THEN
                  block_ptr(1:nc, 1:nr) = 0.5_dp*( &
                                          (dense_mat(r1t:r2t, c1t:c2t) + zsign*dense_zeta(r1t:r2t, c1t:c2t)) + &
                                          TRANSPOSE(dense_mat(r1:r2, c1:c2) + zsign*dense_zeta(r1:r2, c1:c2)))
               ELSE
                  block_ptr(1:nc, 1:nr) = 0.5_dp*( &
                                          dense_mat(r1t:r2t, c1t:c2t) + TRANSPOSE(dense_mat(r1:r2, c1:c2)))
               END IF
            ELSE
               IF (PRESENT(dense_zeta)) THEN
                  block_ptr(1:nr, 1:nc) = 0.5_dp*( &
                                          (dense_mat(r1:r2, c1:c2) + zsign*dense_zeta(r1:r2, c1:c2)) + &
                                          TRANSPOSE(dense_mat(r1t:r2t, c1t:c2t) + zsign*dense_zeta(r1t:r2t, c1t:c2t)))
               ELSE
                  block_ptr(1:nr, 1:nc) = 0.5_dp*( &
                                          dense_mat(r1:r2, c1:c2) + TRANSPOSE(dense_mat(r1t:r2t, c1t:c2t)))
               END IF
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF
   END SUBROUTINE dense_to_dbcsr_into

! **************************************************************************************************
!> \brief Whether an existing DBCSR matrix can receive the cached block layout.
!>
!> Reusing a matrix with matching totals but a different matrix type, block
!> sizes or distribution would push blocks to the wrong owners. Recreate from
!> the template unless the target matches the layout exactly.
!> \param matrix ...
!> \param layout ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_target_matches_layout(matrix, layout) RESULT(matches)
      USE cp_dbcsr_api, ONLY: dbcsr_distribution_get, dbcsr_distribution_type, &
                              dbcsr_get_info, dbcsr_type, dbcsr_type_symmetric
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(dense_to_dbcsr_layout_type), INTENT(IN)       :: layout
      LOGICAL                                            :: matches

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: mynode, mypcol, myprow, nblkcols, &
                                                            nblkrows, ncols, nrows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist
      TYPE(dbcsr_distribution_type)                      :: dist

      matches = .FALSE.
      IF (.NOT. ASSOCIATED(matrix)) RETURN

      CALL dbcsr_get_info(matrix, nfullrows_total=nrows, nfullcols_total=ncols, &
                          nblkrows_total=nblkrows, nblkcols_total=nblkcols, &
                          matrix_type=matrix_type, distribution=dist, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      IF (matrix_type /= dbcsr_type_symmetric) RETURN
      IF (nrows /= layout%nrows .OR. ncols /= layout%ncols) RETURN
      IF (nblkrows /= layout%nblkrows .OR. nblkcols /= layout%nblkcols) RETURN
      IF (SIZE(row_blk_size) /= SIZE(layout%row_blk_size)) RETURN
      IF (SIZE(col_blk_size) /= SIZE(layout%col_blk_size)) RETURN
      IF (ANY(row_blk_size /= layout%row_blk_size)) RETURN
      IF (ANY(col_blk_size /= layout%col_blk_size)) RETURN

      CALL dbcsr_distribution_get(dist, mynode=mynode, myprow=myprow, mypcol=mypcol, &
                                  row_dist=row_dist, col_dist=col_dist)
      IF (mynode /= layout%mynode .OR. myprow /= layout%myprow .OR. &
          mypcol /= layout%mypcol) RETURN
      IF (SIZE(row_dist) /= SIZE(layout%row_dist)) RETURN
      IF (SIZE(col_dist) /= SIZE(layout%col_dist)) RETURN
      IF (ANY(row_dist /= layout%row_dist)) RETURN
      IF (ANY(col_dist /= layout%col_dist)) RETURN

      matches = .TRUE.

   END FUNCTION dbcsr_target_matches_layout

! **************************************************************************************************
!> \brief Whether the cached dense->DBCSR block layout still matches the template.
!> \param layout ...
!> \param nrows ...
!> \param ncols ...
!> \param nblkrows ...
!> \param nblkcols ...
!> \param mynode ...
!> \param myprow ...
!> \param mypcol ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param row_dist ...
!> \param col_dist ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_layout_matches(layout, nrows, ncols, nblkrows, nblkcols, mynode, myprow, mypcol, &
                                 row_blk_size, col_blk_size, row_dist, col_dist) RESULT(matches)
      TYPE(dense_to_dbcsr_layout_type), INTENT(IN)       :: layout
      INTEGER, INTENT(IN)                                :: nrows, ncols, nblkrows, nblkcols, &
                                                            mynode, myprow, mypcol
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_blk_size, col_blk_size, row_dist, &
                                                            col_dist
      LOGICAL                                            :: matches

      matches = layout%valid .AND. layout%nrows == nrows .AND. layout%ncols == ncols .AND. &
                layout%nblkrows == nblkrows .AND. layout%nblkcols == nblkcols .AND. &
                layout%mynode == mynode .AND. layout%myprow == myprow .AND. &
                layout%mypcol == mypcol
      IF (matches) THEN
         matches = ALLOCATED(layout%row_blk_size) .AND. ALLOCATED(layout%col_blk_size) .AND. &
                   ALLOCATED(layout%row_dist) .AND. ALLOCATED(layout%col_dist)
      END IF
      IF (matches) THEN
         matches = SIZE(layout%row_blk_size) == SIZE(row_blk_size) .AND. &
                   SIZE(layout%col_blk_size) == SIZE(col_blk_size) .AND. &
                   SIZE(layout%row_dist) == SIZE(row_dist) .AND. &
                   SIZE(layout%col_dist) == SIZE(col_dist)
      END IF
      IF (matches) THEN
         matches = ALL(layout%row_blk_size == row_blk_size) .AND. &
                   ALL(layout%col_blk_size == col_blk_size)
      END IF
      IF (matches) THEN
         matches = ALL(layout%row_dist == row_dist) .AND. ALL(layout%col_dist == col_dist)
      END IF
   END FUNCTION dense_layout_matches

! **************************************************************************************************
!> \brief Build the cached offsets, owned upper-triangular block list and workspace.
!> \param layout ...
!> \param nrows ...
!> \param ncols ...
!> \param nblkrows ...
!> \param nblkcols ...
!> \param mynode ...
!> \param myprow ...
!> \param mypcol ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param row_dist ...
!> \param col_dist ...
! **************************************************************************************************
   SUBROUTINE build_dense_layout(layout, nrows, ncols, nblkrows, nblkcols, mynode, myprow, mypcol, &
                                 row_blk_size, col_blk_size, row_dist, col_dist)
      TYPE(dense_to_dbcsr_layout_type), INTENT(INOUT)    :: layout
      INTEGER, INTENT(IN)                                :: nrows, ncols, nblkrows, nblkcols, &
                                                            mynode, myprow, mypcol
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_blk_size, col_blk_size, row_dist, &
                                                            col_dist

      INTEGER                                            :: a, b, col, icol, irow, nblocks, &
                                                            nlocal_cols, nlocal_rows, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_cols, local_rows

      IF (ALLOCATED(layout%row_blk_size)) THEN
         DEALLOCATE (layout%row_blk_size, layout%col_blk_size, layout%row_dist, layout%col_dist, &
                     layout%r_offset, layout%c_offset, layout%blk_row, layout%blk_col, layout%blk)
      END IF

      ALLOCATE (layout%row_blk_size(nblkrows), layout%col_blk_size(nblkcols), &
                layout%row_dist(nblkrows), layout%col_dist(nblkcols), &
                layout%r_offset(nblkrows), layout%c_offset(nblkcols))
      layout%row_blk_size(:) = row_blk_size
      layout%col_blk_size(:) = col_blk_size
      layout%row_dist(:) = row_dist
      layout%col_dist(:) = col_dist
      layout%r_offset(1) = 1
      DO row = 2, nblkrows
         layout%r_offset(row) = layout%r_offset(row - 1) + row_blk_size(row - 1)
      END DO
      layout%c_offset(1) = 1
      DO col = 2, nblkcols
         layout%c_offset(col) = layout%c_offset(col - 1) + col_blk_size(col - 1)
      END DO

      ALLOCATE (local_rows(nblkrows), local_cols(nblkcols))
      nlocal_rows = 0
      DO row = 1, nblkrows
         IF (row_dist(row) == myprow) THEN
            nlocal_rows = nlocal_rows + 1
            local_rows(nlocal_rows) = row
         END IF
      END DO
      nlocal_cols = 0
      DO col = 1, nblkcols
         IF (col_dist(col) == mypcol) THEN
            nlocal_cols = nlocal_cols + 1
            local_cols(nlocal_cols) = col
         END IF
      END DO

      nblocks = 0
      DO a = 1, nlocal_rows
         irow = local_rows(a)
         DO b = 1, nlocal_cols
            icol = local_cols(b)
            IF (irow <= icol .AND. .NOT. BTEST(irow + icol, 0)) nblocks = nblocks + 1
         END DO
      END DO
      DO b = 1, nlocal_cols
         irow = local_cols(b)
         DO a = 1, nlocal_rows
            icol = local_rows(a)
            IF (irow <= icol .AND. BTEST(irow + icol, 0)) nblocks = nblocks + 1
         END DO
      END DO

      ALLOCATE (layout%blk_row(nblocks), layout%blk_col(nblocks))
      ALLOCATE (layout%blk(MAXVAL(row_blk_size), MAXVAL(col_blk_size)))
      nblocks = 0
      DO a = 1, nlocal_rows
         irow = local_rows(a)
         DO b = 1, nlocal_cols
            icol = local_cols(b)
            IF (irow <= icol .AND. .NOT. BTEST(irow + icol, 0)) THEN
               nblocks = nblocks + 1
               layout%blk_row(nblocks) = irow
               layout%blk_col(nblocks) = icol
            END IF
         END DO
      END DO
      DO b = 1, nlocal_cols
         irow = local_cols(b)
         DO a = 1, nlocal_rows
            icol = local_rows(a)
            IF (irow <= icol .AND. BTEST(irow + icol, 0)) THEN
               nblocks = nblocks + 1
               layout%blk_row(nblocks) = irow
               layout%blk_col(nblocks) = icol
            END IF
         END DO
      END DO
      DEALLOCATE (local_rows, local_cols)

      layout%nrows = nrows
      layout%ncols = ncols
      layout%nblkrows = nblkrows
      layout%nblkcols = nblkcols
      layout%mynode = mynode
      layout%myprow = myprow
      layout%mypcol = mypcol
      layout%valid = .TRUE.
   END SUBROUTINE build_dense_layout

! **************************************************************************************************
!> \brief ...
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_gauxc_functional(xc_section) RESULT(gauxc_functional_section)
      TYPE(section_vals_type), INTENT(in), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: gauxc_functional_section

      INTEGER                                            :: n_subsections

      CALL find_single_functional_subsection(xc_section, "GAUXC", gauxc_functional_section, &
                                             n_subsections)
      IF (n_subsections == 0) THEN
         CPABORT("No XC functional found in XC_FUNCTIONAL section")
      END IF
      IF (n_subsections > 1 .OR. .NOT. ASSOCIATED(gauxc_functional_section)) THEN
         CPABORT("GauXC functionals are mutually exclusive with any other functional.")
      END IF
   END FUNCTION get_gauxc_functional

! **************************************************************************************************
!> \brief Return whether GauXC GAPW mode sees pseudopotential kinds.
!> \param qs_kind_set ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_gapw_has_pseudopotentials(qs_kind_set) RESULT(has_pseudopotentials)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL                                            :: has_pseudopotentials

      LOGICAL                                            :: all_all_electron, all_pseudopotentials, &
                                                            has_paw_pseudopotentials

      CALL gauxc_gapw_kind_flags(qs_kind_set, has_pseudopotentials, all_pseudopotentials, &
                                 has_paw_pseudopotentials, all_all_electron)

   END FUNCTION gauxc_gapw_has_pseudopotentials

! **************************************************************************************************
!> \brief Return whether every GAPW kind uses a pseudopotential.
!> \param qs_kind_set ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_gapw_all_pseudopotentials(qs_kind_set) RESULT(all_pseudopotentials)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL                                            :: all_pseudopotentials

      LOGICAL                                            :: all_all_electron, &
                                                            has_paw_pseudopotentials, &
                                                            has_pseudopotentials

      CALL gauxc_gapw_kind_flags(qs_kind_set, has_pseudopotentials, all_pseudopotentials, &
                                 has_paw_pseudopotentials, all_all_electron)

   END FUNCTION gauxc_gapw_all_pseudopotentials

! **************************************************************************************************
!> \brief Return whether GauXC GAPW mode sees pseudopotential one-center GAPW kinds.
!> \param qs_kind_set ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_gapw_has_paw_pseudopotentials(qs_kind_set) RESULT(has_paw_pseudopotentials)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL                                            :: has_paw_pseudopotentials

      LOGICAL                                            :: all_all_electron, all_pseudopotentials, &
                                                            has_pseudopotentials

      CALL gauxc_gapw_kind_flags(qs_kind_set, has_pseudopotentials, all_pseudopotentials, &
                                 has_paw_pseudopotentials, all_all_electron)

   END FUNCTION gauxc_gapw_has_paw_pseudopotentials

! **************************************************************************************************
!> \brief Check the current periodic scope of the CP2K-GauXC bridge
!> \param dft_control ...
!> \param cell ...
!> \param qs_kind_set ...
!> \param do_kpoints ...
!> \param periodic_reference ...
!> \param is_periodic ...
!> \note This path keeps isolated validation cells usable under PERIODIC XYZ.
!>       It intentionally does not implement compact periodic GauXC quadrature.
! **************************************************************************************************
   SUBROUTINE ensure_gauxc_periodic_reference_scope( &
      dft_control, cell, qs_kind_set, do_kpoints, periodic_reference, is_periodic)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL, INTENT(IN)                                :: do_kpoints, periodic_reference, &
                                                            is_periodic

      INTEGER                                            :: ikind
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(qs_kind_set))

      IF (do_kpoints) THEN
         CALL cp_abort(__LOCATION__, &
                       "GauXC currently supports only Gamma-only density matrices in CP2K. "// &
                       "Periodic k-point density matrices require a dedicated GauXC periodic interface.")
      END IF
      IF (dft_control%nimages /= 1) THEN
         CALL cp_abort(__LOCATION__, &
                       "GauXC currently supports only a single AO image in CP2K. "// &
                       "Periodic neighbour-cell AO blocks require a dedicated GauXC periodic interface.")
      END IF
      IF (.NOT. is_periodic) RETURN

      IF (.NOT. periodic_reference) THEN
         CALL cp_abort(__LOCATION__, &
                       "Periodic GauXC calculations in CP2K require GAUXC%PERIODIC_REFERENCE T. "// &
                       "This opt-in documents that the current path is only an isolated-cell, "// &
                       "Gamma-only, single-image METHOD GPW reference path using GauXC molecular "// &
                       "quadrature, not a dedicated periodic GauXC interface.")
      END IF

      IF (.NOT. ALL(cell%perd == 1)) THEN
         CALL cp_abort(__LOCATION__, &
                       "The current GauXC isolated-cell reference path supports only PERIODIC XYZ. "// &
                       "Partial periodicity requires a dedicated GauXC periodic interface.")
      END IF
      IF (.NOT. dft_control%qs_control%gpw) THEN
         CALL cp_abort(__LOCATION__, &
                       "The current GauXC isolated-cell reference path is limited to METHOD GPW with GTH "// &
                       "pseudopotentials. GAPW, GAPW_XC, and other QS methods are not supported here.")
      END IF

      DO ikind = 1, SIZE(qs_kind_set)
         NULLIFY (gth_potential, sgp_potential)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          gth_potential=gth_potential, &
                          sgp_potential=sgp_potential)
         IF (.NOT. ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
            CALL cp_abort(__LOCATION__, &
                          "The current GauXC isolated-cell reference path is limited to GTH pseudopotentials. "// &
                          "Use non-periodic all-electron GAPW validation for molecular GAPW cases.")
         END IF
      END DO

   END SUBROUTINE ensure_gauxc_periodic_reference_scope

! **************************************************************************************************
!> \brief Verify once that the GauXC energy gradient is replicated across MPI ranks.
!>        The force/virial consumers add the gradient on rank 0 only, which is
!>        correct only for a replicated vector. GauXC's C-API integrators
!>        allreduce exc_grad, but make that assumption explicit and audible.
!> \param exc_grad ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE assert_gauxc_gradient_replicated(exc_grad, para_env)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(mp_para_env_type), POINTER                    :: para_env

      REAL(KIND=dp)                                      :: max_diff, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: reference

      IF (para_env%num_pe <= 1) RETURN
      ALLOCATE (reference(SIZE(exc_grad)))
      reference(:) = exc_grad
      CALL para_env%bcast(reference, 0)
      max_diff = MAXVAL(ABS(exc_grad - reference))
      scale = MAX(1.0_dp, MAXVAL(ABS(exc_grad)), MAXVAL(ABS(reference)))
      DEALLOCATE (reference)
      IF (max_diff > 1.0E-10_dp*scale) THEN
         CALL cp_abort(__LOCATION__, &
                       "GauXC XC gradient is not replicated across MPI ranks, but the "// &
                       "force/virial consumers add it on rank 0 only. This indicates an "// &
                       "unsupported non-replicated integrator.")
      END IF

   END SUBROUTINE assert_gauxc_gradient_replicated

! **************************************************************************************************
!> \brief adds a replicated GauXC energy gradient to the local CP2K force accumulator
!> \param exc_grad ...
!> \param force ...
!> \param atomic_kind_set ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE add_gauxc_gradient_to_force(exc_grad, force, atomic_kind_set, para_env)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

      CPASSERT(ASSOCIATED(force))
      CPASSERT(ASSOCIATED(atomic_kind_set))

      IF (para_env%mepos /= 0) RETURN

      DO ikind = 1, SIZE(atomic_kind_set, 1)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
         DO ia = 1, natom_kind
            iatom = atomic_kind%atom_list(ia)
            force(ikind)%rho_elec(:, ia) = force(ikind)%rho_elec(:, ia) + &
                                           exc_grad(3*iatom - 2:3*iatom)
         END DO
      END DO

   END SUBROUTINE add_gauxc_gradient_to_force

! **************************************************************************************************
!> \brief compute a GauXC XC energy for diagnostic finite differences
!> \param particle_set_eval ...
!> \param qs_kind_set ...
!> \param density_scalar ...
!> \param nspins ...
!> \param model_name ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param lwd_kernel ...
!> \param batch_size ...
!> \param device_runtime_fill_fraction ...
!> \param res ...
!> \param density_zeta ...
!> \param runtime ...
! **************************************************************************************************
   SUBROUTINE gauxc_xc_energy_for_particles( &
      particle_set_eval, qs_kind_set, density_scalar, nspins, model_name, &
      xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
      int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, res, density_zeta, &
      runtime)
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set_eval
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_scalar
      INTEGER, INTENT(IN)                                :: nspins
      CHARACTER(len=*), INTENT(IN) :: model_name, xc_fun_name, grid_type, radial_quadrature, &
         pruning_scheme, lb_exec_space, int_exec_space, lwd_kernel
      INTEGER, INTENT(IN)                                :: batch_size
      REAL(KIND=dp), INTENT(IN)                          :: device_runtime_fill_fraction
      TYPE(cp_gauxc_xc_type), INTENT(INOUT)              :: res
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta
      TYPE(cp_gauxc_runtime_type), INTENT(IN), OPTIONAL  :: runtime

      TYPE(cp_gauxc_basisset_type)                       :: gauxc_basis_fd
      TYPE(cp_gauxc_grid_type)                           :: gauxc_grid_fd
      TYPE(cp_gauxc_integrator_type)                     :: gauxc_integrator_fd
      TYPE(cp_gauxc_molecule_type)                       :: gauxc_mol_fd
      TYPE(cp_gauxc_status_type)                         :: gauxc_status

      gauxc_mol_fd = gauxc_create_molecule(particle_set_eval, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_basis_fd = gauxc_create_basisset(qs_kind_set, particle_set_eval, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      IF (PRESENT(runtime)) THEN
         gauxc_grid_fd = gauxc_create_grid( &
                         gauxc_mol_fd, &
                         gauxc_basis_fd, &
                         grid_type, &
                         radial_quadrature, &
                         pruning_scheme, &
                         lb_exec_space, &
                         batch_size, &
                         device_runtime_fill_fraction, &
                         gauxc_status, &
                         mpi_comm=mp_comm_self%get_handle(), &
                         int_exec_space=int_exec_space, &
                         runtime=runtime)
      ELSE
         gauxc_grid_fd = gauxc_create_grid( &
                         gauxc_mol_fd, &
                         gauxc_basis_fd, &
                         grid_type, &
                         radial_quadrature, &
                         pruning_scheme, &
                         lb_exec_space, &
                         batch_size, &
                         device_runtime_fill_fraction, &
                         gauxc_status, &
                         mpi_comm=mp_comm_self%get_handle(), &
                         force_new_runtime=.TRUE., &
                         int_exec_space=int_exec_space)
      END IF
      CALL gauxc_check_status(gauxc_status)
      gauxc_integrator_fd = gauxc_create_integrator( &
                            TRIM(xc_fun_name), &
                            gauxc_grid_fd, &
                            int_exec_space, &
                            lwd_kernel, &
                            nspins, &
                            gauxc_status)
      CALL gauxc_check_status(gauxc_status)

      IF (nspins == 1) THEN
         CALL gauxc_compute_xc( &
            gauxc_integrator_fd, &
            density_scalar, &
            res, &
            nspins=nspins, &
            status=gauxc_status, &
            model=TRIM(model_name), &
            energy_only=.TRUE.)
      ELSE
         CPASSERT(nspins == 2)
         CPASSERT(PRESENT(density_zeta))
         CALL gauxc_compute_xc( &
            gauxc_integrator_fd, &
            density_scalar, &
            res, &
            density_zeta, &
            nspins, &
            gauxc_status, &
            model=TRIM(model_name), &
            energy_only=.TRUE.)
      END IF
      CALL gauxc_check_status(gauxc_status)

      CALL gauxc_destroy_integrator(gauxc_integrator_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_grid(gauxc_grid_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_basisset(gauxc_basis_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_molecule(gauxc_mol_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)

   END SUBROUTINE gauxc_xc_energy_for_particles

! **************************************************************************************************
!> \brief compute a finite-difference GauXC XC nuclear gradient at fixed density
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param density_scalar ...
!> \param nspins ...
!> \param model_name ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param lwd_kernel ...
!> \param batch_size ...
!> \param device_runtime_fill_fraction ...
!> \param dx ...
!> \param para_env ...
!> \param exc_grad ...
!> \param density_zeta ...
! **************************************************************************************************
   SUBROUTINE gauxc_xc_gradient_fd( &
      particle_set, qs_kind_set, density_scalar, nspins, model_name, &
      xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
      int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, dx, para_env, exc_grad, &
      density_zeta)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_scalar
      INTEGER, INTENT(IN)                                :: nspins
      CHARACTER(len=*), INTENT(IN) :: model_name, xc_fun_name, grid_type, radial_quadrature, &
         pruning_scheme, lb_exec_space, int_exec_space, lwd_kernel
      INTEGER, INTENT(IN)                                :: batch_size
      REAL(KIND=dp), INTENT(IN)                          :: device_runtime_fill_fraction, dx
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: exc_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta

      CHARACTER(len=default_string_length)               :: int_key, lb_key
      INTEGER                                            :: iatom, idir, work_index
      LOGICAL                                            :: distribute_fd, has_work, int_device, &
                                                            lb_device, mixed_spaces, &
                                                            use_onedft_model
      REAL(KIND=dp)                                      :: xc_minus, xc_plus
      TYPE(cp_gauxc_runtime_type)                        :: fd_runtime
      TYPE(cp_gauxc_status_type)                         :: gauxc_status
      TYPE(cp_gauxc_xc_type)                             :: xc_result
      TYPE(particle_type), ALLOCATABLE, DIMENSION(:)     :: particle_set_minus, particle_set_plus

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(dx > 0.0_dp)

      IF (ALLOCATED(exc_grad)) THEN
         IF (SIZE(exc_grad) /= 3*SIZE(particle_set)) DEALLOCATE (exc_grad)
      END IF
      IF (.NOT. ALLOCATED(exc_grad)) ALLOCATE (exc_grad(3*SIZE(particle_set)))
      exc_grad = 0.0_dp

      use_onedft_model = gauxc_model_active(model_name)
      lb_key = ADJUSTL(lb_exec_space)
      int_key = ADJUSTL(int_exec_space)
      CALL uppercase(lb_key)
      CALL uppercase(int_key)
      lb_device = (TRIM(lb_key) == "DEVICE")
      int_device = (TRIM(int_key) == "DEVICE")
      mixed_spaces = (lb_device .NEQV. int_device)
      IF (use_onedft_model .AND. para_env%num_pe > 1 .AND. mixed_spaces .AND. &
          .NOT. fd_mixed_space_warned) THEN
         fd_mixed_space_warned = .TRUE.
         CALL cp_warn(__LOCATION__, &
                      "GauXC finite-difference gradients do not support mixed "// &
                      "LB_EXECUTION_SPACE/INT_EXECUTION_SPACE settings; the FD sweep "// &
                      "will stay on rank 0.")
      END IF
      distribute_fd = (para_env%num_pe > 1) .AND. &
                      (.NOT. use_onedft_model .OR. (.NOT. lb_device .AND. .NOT. int_device))
      has_work = distribute_fd .OR. (para_env%mepos == 0)

      IF (has_work) THEN
         CALL gauxc_create_shared_runtime(fd_runtime, mp_comm_self%get_handle(), lb_exec_space, &
                                          int_exec_space, device_runtime_fill_fraction, gauxc_status)
         CALL gauxc_check_status(gauxc_status)
         ALLOCATE (particle_set_minus(SIZE(particle_set)), particle_set_plus(SIZE(particle_set)))

         DO iatom = 1, SIZE(particle_set)
            DO idir = 1, 3
               work_index = 3*(iatom - 1) + (idir - 1)
               IF (distribute_fd) THEN
                  IF (MOD(work_index, para_env%num_pe) /= para_env%mepos) CYCLE
               END IF

               particle_set_minus = particle_set
               particle_set_plus = particle_set
               particle_set_minus(iatom)%r(idir) = particle_set_minus(iatom)%r(idir) - dx
               particle_set_plus(iatom)%r(idir) = particle_set_plus(iatom)%r(idir) + dx
               IF (PRESENT(density_zeta)) THEN
                  CALL gauxc_xc_energy_for_particles(particle_set_plus, qs_kind_set, density_scalar, &
                                                     nspins, model_name, xc_fun_name, grid_type, &
                                                     radial_quadrature, pruning_scheme, lb_exec_space, &
                                                     int_exec_space, lwd_kernel, batch_size, &
                                                     device_runtime_fill_fraction, xc_result, &
                                                     runtime=fd_runtime, density_zeta=density_zeta)
                  xc_plus = xc_result%exc
                  CALL gauxc_xc_energy_for_particles(particle_set_minus, qs_kind_set, density_scalar, &
                                                     nspins, model_name, xc_fun_name, grid_type, &
                                                     radial_quadrature, pruning_scheme, lb_exec_space, &
                                                     int_exec_space, lwd_kernel, batch_size, &
                                                     device_runtime_fill_fraction, xc_result, &
                                                     runtime=fd_runtime, density_zeta=density_zeta)
                  xc_minus = xc_result%exc
               ELSE
                  CALL gauxc_xc_energy_for_particles(particle_set_plus, qs_kind_set, density_scalar, &
                                                     nspins, model_name, xc_fun_name, grid_type, &
                                                     radial_quadrature, pruning_scheme, lb_exec_space, &
                                                     int_exec_space, lwd_kernel, batch_size, &
                                                     device_runtime_fill_fraction, xc_result, &
                                                     runtime=fd_runtime)
                  xc_plus = xc_result%exc
                  CALL gauxc_xc_energy_for_particles(particle_set_minus, qs_kind_set, density_scalar, &
                                                     nspins, model_name, xc_fun_name, grid_type, &
                                                     radial_quadrature, pruning_scheme, lb_exec_space, &
                                                     int_exec_space, lwd_kernel, batch_size, &
                                                     device_runtime_fill_fraction, xc_result, &
                                                     runtime=fd_runtime)
                  xc_minus = xc_result%exc
               END IF
               exc_grad(3*iatom - 3 + idir) = (xc_plus - xc_minus)/(2.0_dp*dx)
            END DO
         END DO

         DEALLOCATE (particle_set_minus, particle_set_plus)

         CALL gauxc_destroy_shared_runtime(fd_runtime, gauxc_status)
         CALL gauxc_check_status(gauxc_status)
      END IF

      IF (distribute_fd) THEN
         CALL para_env%sum(exc_grad)
      ELSE
         CALL para_env%bcast(exc_grad, 0)
      END IF

   END SUBROUTINE gauxc_xc_gradient_fd

! **************************************************************************************************
!> \brief Arithmetic center of the particle positions.
!> \param particle_set ...
!> \return ...
! **************************************************************************************************
   FUNCTION molecular_center(particle_set) RESULT(center)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(3)                        :: center

      INTEGER                                            :: iatom

      CPASSERT(ASSOCIATED(particle_set))
      center = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         center = center + particle_set(iatom)%r
      END DO
      center = center/REAL(SIZE(particle_set), dp)

   END FUNCTION molecular_center

! **************************************************************************************************
!> \brief finite-difference check of the molecular GauXC XC virial diagnostic
!> \param exc_grad ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param density_scalar ...
!> \param nspins ...
!> \param model_name ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param lwd_kernel ...
!> \param batch_size ...
!> \param device_runtime_fill_fraction ...
!> \param dx ...
!> \param para_env ...
!> \param density_zeta ...
! **************************************************************************************************
   SUBROUTINE debug_gauxc_molecular_virial( &
      exc_grad, particle_set, qs_kind_set, density_scalar, nspins, model_name, &
      xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
      int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, dx, para_env, density_zeta)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_scalar
      INTEGER, INTENT(IN)                                :: nspins
      CHARACTER(len=*), INTENT(IN) :: model_name, xc_fun_name, grid_type, radial_quadrature, &
         pruning_scheme, lb_exec_space, int_exec_space, lwd_kernel
      INTEGER, INTENT(IN)                                :: batch_size
      REAL(KIND=dp), INTENT(IN)                          :: device_runtime_fill_fraction, dx
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta

      INTEGER                                            :: iatom, iw
      REAL(KIND=dp)                                      :: analytic_trace, diff_trace, &
                                                            numerical_trace, xc_minus, xc_plus
      REAL(KIND=dp), DIMENSION(3)                        :: center, displacement, grad
      TYPE(cp_gauxc_xc_type)                             :: xc_result
      TYPE(particle_type), ALLOCATABLE, DIMENSION(:)     :: particle_set_minus, particle_set_plus

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(exc_grad) == 3*SIZE(particle_set))

      IF (para_env%mepos /= 0) RETURN

      center = molecular_center(particle_set)

      ALLOCATE (particle_set_minus(SIZE(particle_set)), particle_set_plus(SIZE(particle_set)))
      particle_set_minus = particle_set
      particle_set_plus = particle_set

      analytic_trace = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         grad = exc_grad(3*iatom - 2:3*iatom)
         displacement = particle_set(iatom)%r - center
         analytic_trace = analytic_trace + DOT_PRODUCT(grad, displacement)
         particle_set_minus(iatom)%r = center + (1.0_dp - dx)*displacement
         particle_set_plus(iatom)%r = center + (1.0_dp + dx)*displacement
      END DO
      analytic_trace = analytic_trace/3.0_dp

      IF (PRESENT(density_zeta)) THEN
         CALL gauxc_xc_energy_for_particles( &
            particle_set_plus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, xc_result, &
            density_zeta=density_zeta)
         xc_plus = xc_result%exc
         CALL gauxc_xc_energy_for_particles( &
            particle_set_minus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, xc_result, &
            density_zeta=density_zeta)
         xc_minus = xc_result%exc
      ELSE
         CALL gauxc_xc_energy_for_particles( &
            particle_set_plus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, xc_result)
         xc_plus = xc_result%exc
         CALL gauxc_xc_energy_for_particles( &
            particle_set_minus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, lwd_kernel, batch_size, device_runtime_fill_fraction, xc_result)
         xc_minus = xc_result%exc
      END IF

      numerical_trace = (xc_plus - xc_minus)/(2.0_dp*dx)/3.0_dp
      diff_trace = analytic_trace - numerical_trace

      iw = cp_logger_get_default_io_unit()
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,T2,A,1X,ES11.4)") &
            "GAUXC| Molecular XC virial finite-difference dx", dx
         WRITE (UNIT=iw, FMT="(T2,A,3(1X,ES19.11))") &
            "GAUXC| Molecular XC virial FD 1/3 Trace", &
            analytic_trace, numerical_trace, diff_trace
      END IF

      DEALLOCATE (particle_set_minus, particle_set_plus)

   END SUBROUTINE debug_gauxc_molecular_virial

! **************************************************************************************************
!> \brief prints a force-based molecular XC virial diagnostic from GauXC gradients
!> \param exc_grad ...
!> \param particle_set ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE print_gauxc_molecular_virial(exc_grad, particle_set, para_env)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=1), DIMENSION(3), PARAMETER          :: label = ["x", "y", "z"]

      INTEGER                                            :: i, iatom, iw, j
      REAL(KIND=dp), DIMENSION(3)                        :: center, displacement, grad, grad_sum
      REAL(KIND=dp), DIMENSION(3, 3)                     :: molecular_virial

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(exc_grad) == 3*SIZE(particle_set))

      IF (para_env%mepos /= 0) RETURN

      center = molecular_center(particle_set)

      grad_sum = 0.0_dp
      molecular_virial = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         grad = exc_grad(3*iatom - 2:3*iatom)
         displacement = particle_set(iatom)%r - center
         grad_sum = grad_sum + grad
         DO i = 1, 3
            DO j = 1, 3
               molecular_virial(i, j) = molecular_virial(i, j) + grad(i)*displacement(j)
            END DO
         END DO
      END DO

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN

      WRITE (UNIT=iw, FMT="(/,T2,A)") &
         "GAUXC| Molecular XC gradient virial diagnostic [a.u.]"
      WRITE (UNIT=iw, FMT="(T2,A,T20,A,T40,A,T60,A)") "GAUXC|", "x", "y", "z"
      DO i = 1, 3
         WRITE (UNIT=iw, FMT="(T2,A,1X,A1,3(1X,ES19.11))") &
            "GAUXC|", label(i), molecular_virial(i, :)
      END DO
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "GAUXC| Molecular XC gradient virial 1/3 Trace", &
         (molecular_virial(1, 1) + molecular_virial(2, 2) + molecular_virial(3, 3))/3.0_dp
      WRITE (UNIT=iw, FMT="(T2,A,3(1X,ES19.11))") &
         "GAUXC| Molecular XC gradient sum", grad_sum
      WRITE (UNIT=iw, FMT="(T2,A)") &
         "GAUXC| Diagnostic only; this is not an analytical periodic stress tensor."

   END SUBROUTINE print_gauxc_molecular_virial

! **************************************************************************************************
!> \brief Return information about the Skala functional
!> \param functional section containing the SKALA subsection
!> \param lsd if you are using lsd or lda
!> \param reference the reference to the article where the functional is explained
!> \param shortform the short definition of the functional
!> \param needs the flags corresponding to the inputs needed by this
!>        functional are set to true (the flags not needed aren't touched)
!> \param max_deriv the maximal derivative available
! **************************************************************************************************
   SUBROUTINE skala_info(functional, lsd, reference, shortform, needs, max_deriv)
      TYPE(section_vals_type), POINTER                   :: functional
      LOGICAL, INTENT(in)                                :: lsd
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      CHARACTER(len=default_path_length)                 :: model_name
      CHARACTER(len=default_string_length)               :: xc_fun_name
      LOGICAL                                            :: model_is_none

      IF (functional%section%name == "SKALA_NATIVE") THEN
         CALL section_vals_val_get(functional, "MODEL", c_val=model_name)
         IF (PRESENT(reference)) reference = "CP2K-native SKALA TorchScript model "//TRIM(model_name)
         IF (PRESENT(shortform)) shortform = "SKALA (native)"
         IF (PRESENT(needs)) THEN
            IF (lsd) THEN
               needs%rho_spin = .TRUE.
               needs%drho_spin = .TRUE.
               needs%tau_spin = .TRUE.
            ELSE
               needs%rho = .TRUE.
               needs%drho = .TRUE.
               needs%tau = .TRUE.
            END IF
         END IF
         IF (PRESENT(max_deriv)) max_deriv = 1
         RETURN
      END IF

      CALL section_vals_val_get(functional, "FUNCTIONAL", c_val=xc_fun_name)
      CALL section_vals_val_get(functional, "MODEL", c_val=model_name)
      model_is_none = model_key_is_none(model_name, xc_fun_name)

      IF (PRESENT(reference)) THEN
         IF (model_is_none) THEN
            reference = "Functional computed by GauXC (underlying: "//TRIM(xc_fun_name)//")"
         ELSE
            reference = "Functional computed by GauXC Skala model "//TRIM(model_name)
         END IF
      END IF
      IF (PRESENT(shortform)) THEN
         IF (model_is_none) THEN
            shortform = "GAUXC ("//TRIM(xc_fun_name)//")"
         ELSE
            shortform = "GAUXC Skala"
         END IF
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho = .TRUE.
         IF (lsd) THEN
            needs%rho_spin = .TRUE.
         END IF
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1

   END SUBROUTINE skala_info

! GauXC uses replicated dense density and VXC matrices. The DBCSR density matrix
! is distributed over MPI ranks, so apply_gauxc allreduces the dense copy before
! passing it to GauXC.

! **************************************************************************************************
!> \brief Evaluate molecular GauXC, optionally for an explicit configuration density.
!> \param qs_env ...
!> \param xc_section ...
!> \param calculate_forces ...
!> \param rho_ao_external Optional spin/image AO density; supply all three external arguments together
!> \param exc_external Configuration XC energy, leaving the baseline energy unchanged
!> \param matrix_vxc_external Initially null; caller owns the returned spin AO potential matrices
!> \note External-density mode preserves baseline energy and KS storage and shares
!>       the geometry/model cache.
! **************************************************************************************************
   SUBROUTINE apply_gauxc(qs_env, xc_section, calculate_forces, &
                          rho_ao_external, exc_external, matrix_vxc_external)
      TYPE(qs_environment_type), INTENT(in), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(in), POINTER       :: xc_section
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL, POINTER                               :: rho_ao_external
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: exc_external
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL, POINTER                               :: matrix_vxc_external

      CHARACTER(len=*), PARAMETER :: nonlocal_vdw_abort_message = &
         "GauXC does not support non-local VDW_POTENTIAL corrections. "// &
         "Use an additive PAIR_POTENTIAL dispersion correction or disable GauXC."

      CHARACTER(len=default_path_length)                 :: model_key, model_name, output_path, &
                                                            resolved_model
      CHARACTER(len=default_string_length) :: gradient_runtime, gradient_runtime_key, grid_key, &
         int_exec_space_key, lb_exec_space_key, lwd_kernel_key, pruning_key, skala_runtime, &
         skala_runtime_key
      INTEGER                                            :: assembly_handle, density_handle, &
                                                            gpu_mode, ikind, img, nimages
      LOGICAL :: backend_has_gpu, do_kpoints, external_density, gapw_method, &
         gapw_paw_pseudopotentials, gapw_pseudopotentials, grid_explicit, have_zeta, &
         int_space_explicit, is_periodic, lb_space_explicit, low_spin_scope, molecular_virial, &
         molecular_virial_debug, need_xc_gradient, periodic_reference, pruning_explicit, &
         replicated, replicated_fill, use_cuda, use_skala_model, write_hdf5_output
      REAL(KIND=dp)                                      :: exc_value, force_numerical_gradient_dx, &
                                                            molecular_virial_debug_dx
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_gauxc_cache_params)                        :: params
      TYPE(cp_gauxc_cache_type), POINTER                 :: cache
      TYPE(cp_gauxc_status_type)                         :: gauxc_status
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vxc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_use, rho_xc
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: gauxc_functional_section
      TYPE(skala_execution_target_type)                  :: execution_target

      external_density = PRESENT(rho_ao_external)
      CPASSERT(PRESENT(exc_external) .EQV. external_density)
      CPASSERT(PRESENT(matrix_vxc_external) .EQV. external_density)
      IF (external_density) THEN
         CPASSERT(ASSOCIATED(rho_ao_external))
         CPASSERT(.NOT. ASSOCIATED(matrix_vxc_external))
      END IF

      NULLIFY ( &
         atomic_kind_set, &
         cell, &
         dft_control, &
         energy, &
         force, &
         ks_env, &
         matrix_vxc, &
         para_env, &
         particle_set, &
         qs_kind_set, &
         rho, &
         rho_use, &
         rho_xc, &
         rho_ao, &
         scf_env)

      CALL get_qs_env( &
         qs_env, &
         cell=cell, &
         dft_control=dft_control, &
         do_kpoints=do_kpoints, &
         energy=energy, &
         ks_env=ks_env, &
         matrix_vxc=matrix_vxc, &
         natom=params%natom, &
         atomic_kind_set=atomic_kind_set, &
         force=force, &
         para_env=para_env, &
         particle_set=particle_set, &
         qs_kind_set=qs_kind_set, &
         rho=rho, &
         rho_xc=rho_xc, &
         scf_env=scf_env)

      gapw_method = dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc
      gapw_pseudopotentials = gapw_method .AND. &
                              gauxc_gapw_has_pseudopotentials(qs_kind_set)
      gapw_paw_pseudopotentials = gapw_method .AND. &
                                  gauxc_gapw_has_paw_pseudopotentials(qs_kind_set)
      IF (dft_control%qs_control%gapw_xc) THEN
         CPASSERT(ASSOCIATED(rho_xc))
         rho_use => rho_xc
      ELSE
         CPASSERT(ASSOCIATED(rho))
         rho_use => rho
      END IF
      CALL qs_rho_get( &
         rho_use, &
         rho_ao_kp=rho_ao)

      IF (external_density) THEN
         rho_ao => rho_ao_external
         ! get_qs_env returned the baseline pointer. Never release or replace it.
         NULLIFY (matrix_vxc)
      END IF
      nimages = dft_control%nimages
      params%nspins = dft_control%nspins
      params%nimages = dft_control%nimages
      params%nkind = SIZE(qs_kind_set)
      IF (ALLOCATED(params%kind_basis_names)) DEALLOCATE (params%kind_basis_names)
      ALLOCATE (params%kind_basis_names(params%nkind))
      DO ikind = 1, params%nkind
         NULLIFY (basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
         IF (ASSOCIATED(basis_set)) THEN
            params%kind_basis_names(ikind) = TRIM(basis_set%name)
         ELSE
            params%kind_basis_names(ikind) = ""
         END IF
      END DO
      IF (external_density) THEN
         CPASSERT(SIZE(rho_ao, 1) == params%nspins)
         CPASSERT(SIZE(rho_ao, 2) == nimages)
      END IF
      low_spin_scope = external_density .OR. (dft_control%low_spin_roks .AND. gapw_method)
      is_periodic = .FALSE.
      IF (ASSOCIATED(cell)) is_periodic = ANY(cell%perd /= 0)

      IF (ASSOCIATED(qs_env%dispersion_env)) THEN
         IF (qs_env%dispersion_env%type == xc_vdw_fun_nonloc) THEN
            CPABORT(nonlocal_vdw_abort_message)
         END IF
      END IF

      gauxc_functional_section => get_gauxc_functional(xc_section)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "FUNCTIONAL", &
         c_val=params%xc_fun_name)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MODEL", &
         c_val=model_name)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "GRID", &
         c_val=params%grid_type, &
         explicit=grid_explicit)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "RADIAL_QUADRATURE", &
         c_val=params%radial_quadrature)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "PRUNING_SCHEME", &
         c_val=params%pruning_scheme, &
         explicit=pruning_explicit)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "BATCH_SIZE", &
         i_val=params%batch_size)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "DEVICE_RUNTIME_FILL_FRACTION", &
         r_val=params%device_runtime_fill_fraction)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "PERIODIC_REFERENCE", &
         l_val=periodic_reference)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MOLECULAR_VIRIAL", &
         l_val=molecular_virial)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MOLECULAR_VIRIAL_DEBUG", &
         l_val=molecular_virial_debug)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MOLECULAR_VIRIAL_DEBUG_DX", &
         r_val=molecular_virial_debug_dx)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "USE_GPU", &
         i_val=gpu_mode)
      backend_has_gpu = .FALSE.
#if defined(GAUXC_HAS_DEVICE) && defined(__LIBTORCH_CUDA)
      backend_has_gpu = .TRUE.
#endif
      CALL resolve_skala_execution_target(gpu_mode, backend_has_gpu, "GAUXC", execution_target)
      use_cuda = execution_target%use_cuda
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "LB_EXECUTION_SPACE", &
         c_val=params%lb_exec_space, &
         explicit=lb_space_explicit)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "INT_EXECUTION_SPACE", &
         c_val=params%int_exec_space, &
         explicit=int_space_explicit)
      lb_exec_space_key = ADJUSTL(params%lb_exec_space)
      int_exec_space_key = ADJUSTL(params%int_exec_space)
      CALL uppercase(lb_exec_space_key)
      CALL uppercase(int_exec_space_key)
      IF (.NOT. use_cuda) THEN
         IF ((lb_space_explicit .AND. TRIM(lb_exec_space_key) == "DEVICE") .OR. &
             (int_space_explicit .AND. TRIM(int_exec_space_key) == "DEVICE")) THEN
            CALL cp_abort(__LOCATION__, &
                          "GAUXC%USE_GPU disables GPU execution, but an explicit "// &
                          "DEVICE execution space was requested. Use USE_GPU GPU "// &
                          "or remove the explicit DEVICE setting.")
         END IF
         params%lb_exec_space = "HOST"
         params%int_exec_space = "HOST"
      ELSE
         IF (.NOT. lb_space_explicit) params%lb_exec_space = "DEVICE"
         IF (.NOT. int_space_explicit) params%int_exec_space = "DEVICE"
      END IF
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "LWD_KERNEL", &
         c_val=params%lwd_kernel)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "SKALA_RUNTIME", &
         c_val=skala_runtime)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MODEL_GRADIENT_RUNTIME", &
         c_val=gradient_runtime)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "FORCE_NUMERICAL_GRADIENT", &
         l_val=params%force_fd_gradient)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "FORCE_NUMERICAL_GRADIENT_DX", &
         r_val=force_numerical_gradient_dx)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "OUTPUT_PATH", &
         c_val=output_path)

      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      skala_runtime_key = ADJUSTL(skala_runtime)
      CALL uppercase(skala_runtime_key)
      gradient_runtime_key = ADJUSTL(gradient_runtime)
      CALL uppercase(gradient_runtime_key)
      lb_exec_space_key = ADJUSTL(params%lb_exec_space)
      CALL uppercase(lb_exec_space_key)
      int_exec_space_key = ADJUSTL(params%int_exec_space)
      CALL uppercase(int_exec_space_key)
      params%use_gauxc_model = .NOT. model_key_is_none(model_name, params%xc_fun_name)
      use_skala_model = (TRIM(model_key) == "SKALA")
      lwd_kernel_key = ADJUSTL(params%lwd_kernel)
      CALL uppercase(lwd_kernel_key)
      IF (TRIM(lwd_kernel_key) == "AUTO") THEN
         params%lwd_kernel = "DEFAULT"
#if defined(GAUXC_HAS_CUTLASS)
         IF (use_skala_model .AND. TRIM(int_exec_space_key) == "DEVICE") THEN
            params%lwd_kernel = "SCHEME1-CUTLASS"
         END IF
#endif
      END IF
      params%model_eval_name = model_name
      IF (.NOT. params%use_gauxc_model) THEN
         ! MODEL NONE and MODEL equal to FUNCTIONAL select conventional GauXC.
         params%model_eval_name = "NONE"
      END IF
      IF (gapw_pseudopotentials .AND. params%use_gauxc_model .AND. .NOT. dft_control%qs_control%gapw_xc .AND. &
          .NOT. gapw_paw_pseudopotentials .AND. para_env%mepos == 0 .AND. ASSOCIATED(scf_env)) THEN
         IF (scf_env%iter_count == 1) THEN
            CALL cp_warn( &
               __LOCATION__, &
               "GauXC Skala with METHOD GAPW and GPW_TYPE pseudopotentials evaluates "// &
               "the XC term directly on the molecular AO/valence density; no GAPW one-center "// &
               "XC correction is used for those regular-grid kinds.")
         END IF
      END IF
      IF (params%device_runtime_fill_fraction <= 0.0_dp .OR. params%device_runtime_fill_fraction > 1.0_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%DEVICE_RUNTIME_FILL_FRACTION must be > 0 and <= 1.")
      END IF
      IF (force_numerical_gradient_dx <= 0.0_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%FORCE_NUMERICAL_GRADIENT_DX must be positive.")
      END IF
      IF (molecular_virial_debug) THEN
         IF (molecular_virial_debug_dx <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "GauXC MOLECULAR_VIRIAL_DEBUG_DX must be positive.")
         END IF
         molecular_virial = .TRUE.
      END IF
      need_xc_gradient = calculate_forces .OR. molecular_virial
      IF (low_spin_scope) THEN
         IF (.NOT. dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
            CPABORT("GauXC low spin ROKS requires METHOD GAPW.")
         END IF
         IF (.NOT. gauxc_gapw_all_all_electron(qs_kind_set)) THEN
            CPABORT("GAPW low spin ROKS GauXC requires all-electron kinds.")
         END IF
         IF (is_periodic .OR. do_kpoints .OR. nimages /= 1) THEN
            CPABORT("GAPW low spin ROKS GauXC requires an isolated Gamma-point calculation.")
         END IF
         IF (need_xc_gradient) THEN
            CALL cp_abort(__LOCATION__, &
                          "GAPW low spin ROKS GauXC nuclear gradients and molecular virial are not implemented.")
         END IF
         IF (TRIM(int_exec_space_key) /= "HOST" .OR. TRIM(lb_exec_space_key) /= "HOST") THEN
            CPABORT("GAPW low spin ROKS GauXC currently requires HOST execution.")
         END IF
         IF (params%use_gauxc_model .AND. .NOT. use_skala_model) THEN
            CPABORT("GAPW low spin ROKS GauXC supports only the SKALA model or conventional XC.")
         END IF
         IF (use_skala_model .AND. TRIM(skala_runtime_key) == "MPI") THEN
            CPABORT("GAPW low spin ROKS Skala requires SKALA_RUNTIME AUTO or SELF.")
         END IF
      END IF
      CALL ensure_gauxc_periodic_reference_scope( &
         dft_control, cell, qs_kind_set, do_kpoints, periodic_reference, is_periodic)
      IF (is_periodic .AND. periodic_reference .AND. para_env%mepos == 0) THEN
         IF (ASSOCIATED(scf_env)) THEN
            IF (scf_env%iter_count == 1) THEN
               CALL cp_warn( &
                  __LOCATION__, &
                  "GAUXC%PERIODIC_REFERENCE uses GauXC molecular quadrature for isolated validation "// &
                  "cells. Compact periodic materials require a dedicated periodic GauXC interface.")
            END IF
         END IF
      END IF
      IF (has_nlcc(qs_kind_set)) THEN
         CALL cp_abort(__LOCATION__, &
                       "GauXC does not support NLCC pseudopotentials. The molecular GauXC "// &
                       "quadrature does not receive the frozen core density. Use the native "// &
                       "SKALA implementation (SKALA_NATIVE) for NLCC calculations.")
      END IF
      IF (params%use_gauxc_model) THEN
         IF (.NOT. grid_explicit) params%grid_type = "SUPERFINE"
         IF (.NOT. pruning_explicit) params%pruning_scheme = "UNPRUNED"

         grid_key = ADJUSTL(params%grid_type)
         pruning_key = ADJUSTL(params%pruning_scheme)
         CALL uppercase(grid_key)
         CALL uppercase(pruning_key)
         IF (use_skala_model .AND. need_xc_gradient .AND. &
             (TRIM(grid_key) /= "SUPERFINE" .OR. TRIM(pruning_key) /= "UNPRUNED")) THEN
            CALL cp_warn( &
               __LOCATION__, &
               "GauXC Skala nuclear gradients are sensitive to the GauXC molecular grid. "// &
               "Use GRID SUPERFINE and PRUNING_SCHEME UNPRUNED for quantitative force checks.")
         END IF
         IF (use_skala_model) THEN
            CALL resolve_skala_model_path(model_name, TRIM(params%int_exec_space) == "DEVICE", &
                                          resolved_model)
            params%model_eval_name = resolved_model
         END IF
      END IF
      SELECT CASE (TRIM(skala_runtime_key))
      CASE ("AUTO")
         params%use_self_runtime = use_skala_model .AND. para_env%num_pe > 1 .AND. params%nspins > 1
      CASE ("MPI")
         params%use_self_runtime = .FALSE.
      CASE ("SELF")
         params%use_self_runtime = use_skala_model .AND. para_env%num_pe > 1
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "Unknown GAUXC%SKALA_RUNTIME value.")
      END SELECT
      IF (.NOT. use_skala_model) params%use_self_runtime = .FALSE.
      params%gradient_runtime = gradient_runtime_key
      params%gapw_method = gapw_method
      params%need_xc_gradient = need_xc_gradient

      IF (.NOT. ASSOCIATED(qs_env%gauxc_cache)) ALLOCATE (qs_env%gauxc_cache)
      cache => qs_env%gauxc_cache
      CALL gauxc_cache_init( &
         cache, &
         params, &
         para_env, &
         particle_set, &
         qs_kind_set, &
         gauxc_status)

      CALL gauxc_cache_ensure_gradient(cache, params, gauxc_status)
      CALL gauxc_check_status(gauxc_status)

      write_hdf5_output = (TRIM(output_path) /= "") .AND. para_env%mepos == 0 .AND. &
                          .NOT. external_density
      IF (write_hdf5_output .AND. ASSOCIATED(scf_env)) THEN
         write_hdf5_output = scf_env%iter_count == 1
      END IF
      IF (write_hdf5_output) THEN
         CALL gauxc_write_molecule_hdf5( &
            cache%molecule, &
            output_path, &
            "molecule.h5", &
            "molecule", &
            gauxc_status)
         CALL gauxc_check_status(gauxc_status)
         CALL gauxc_write_basisset_hdf5( &
            cache%basisset, &
            output_path, &
            "basisset.h5", &
            "basisset", &
            gauxc_status)
         CALL gauxc_check_status(gauxc_status)
      END IF

      IF (qs_env%run_rtp) THEN
         CPABORT("GAUXC XC energy currently does not support real-time propagation")
      END IF

      exc_value = 0.0_dp
      ! For replicated DBCSR matrices every rank can fill its dense copy locally,
      ! which avoids the collective. If that policy is disabled, broadcast the
      ! root-filled copy instead of summing one nonzero contribution plus zeros.
      replicated_fill = .TRUE.

      IF (ASSOCIATED(matrix_vxc)) THEN
         IF (SIZE(matrix_vxc) /= params%nspins) THEN
            CALL dbcsr_deallocate_matrix_set(matrix_vxc)
            CALL dbcsr_allocate_matrix_set(matrix_vxc, params%nspins)
         END IF
      ELSE
         CALL dbcsr_allocate_matrix_set(matrix_vxc, params%nspins)
      END IF

      DO img = 1, nimages
         IF (img > 1) THEN
            CPABORT("UNIMPLEMENTED: Handling nimg>1 in k-point integration")
         END IF
         have_zeta = .FALSE.
         IF (params%nspins == 1) THEN
            CALL timeset("gauxc_density_dense", density_handle)
            CALL dbcsr_to_dense(rho_ao(1, img), cache%density_scalar, para_env, &
                                replicated=replicated, fill_replicated=replicated_fill)
            IF (.NOT. replicated) THEN
               CALL para_env%sum(cache%density_scalar)
            ELSE IF (.NOT. replicated_fill .AND. para_env%num_pe > 1) THEN
               CALL para_env%bcast(cache%density_scalar, 0)
            END IF
            CALL timestop(density_handle)
         ELSE
            CPASSERT(params%nspins == 2)
            ! In here:
            ! scalar <- rho_ao(1, :) + rho_ao(2, :)
            ! zeta   <- rho_ao(1, :) - rho_ao(2, :)
            CALL timeset("gauxc_density_dense", density_handle)
            CALL dbcsr_to_dense_pair(rho_ao(1, img), cache%density_scalar, rho_ao(2, img), &
                                     cache%density_zeta, para_env, replicated=replicated, &
                                     fill_replicated=replicated_fill)
            IF (.NOT. replicated) THEN
               CALL para_env%sum(cache%density_scalar)
               CALL para_env%sum(cache%density_zeta)
            ELSE IF (.NOT. replicated_fill .AND. para_env%num_pe > 1) THEN
               CALL para_env%bcast(cache%density_scalar, 0)
               CALL para_env%bcast(cache%density_zeta, 0)
            END IF
            CALL timestop(density_handle)
            ! Do NOT reorder the following lines!
            cache%density_scalar(:, :) = cache%density_scalar(:, :) + cache%density_zeta(:, :)
            ! Factor two because the next line is evaluated after the above line.
            ! We need to subtract density_zeta once to undo the above line and
            ! a second time because that is what UKS requires.
            ! This style lowers memory footprint.
            cache%density_zeta(:, :) = cache%density_scalar(:, :) - 2.0_dp*cache%density_zeta(:, :)
            have_zeta = .TRUE.
         END IF

         CALL evaluate_gauxc_forces(cache, params, para_env, particle_set, qs_kind_set, &
                                    force, atomic_kind_set, calculate_forces, need_xc_gradient, &
                                    molecular_virial, molecular_virial_debug, have_zeta, &
                                    force_numerical_gradient_dx, molecular_virial_debug_dx, &
                                    gauxc_status)

         exc_value = exc_value + cache%xc_result%exc

         CALL timeset("gauxc_vxc_assembly", assembly_handle)
         CALL assemble_gauxc_vxc(cache, rho_ao, matrix_vxc, img)
         CALL timestop(assembly_handle)
      END DO

      IF (external_density) THEN
         exc_external = exc_value
         matrix_vxc_external => matrix_vxc
      ELSE
         energy%exc = exc_value
         CALL set_ks_env(ks_env, matrix_vxc=matrix_vxc)
      END IF

   END SUBROUTINE apply_gauxc

! **************************************************************************************************
!> \brief Evaluate GauXC XC and, when requested, its gradient and force/virial consumers.
!> \param cache ...
!> \param params ...
!> \param para_env ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param force ...
!> \param atomic_kind_set ...
!> \param calculate_forces ...
!> \param need_xc_gradient ...
!> \param molecular_virial ...
!> \param molecular_virial_debug ...
!> \param have_zeta ...
!> \param force_numerical_gradient_dx ...
!> \param molecular_virial_debug_dx ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE evaluate_gauxc_forces(cache, params, para_env, particle_set, qs_kind_set, &
                                    force, atomic_kind_set, calculate_forces, need_xc_gradient, &
                                    molecular_virial, molecular_virial_debug, have_zeta, &
                                    force_numerical_gradient_dx, molecular_virial_debug_dx, &
                                    status)
      TYPE(cp_gauxc_cache_type), INTENT(INOUT)           :: cache
      TYPE(cp_gauxc_cache_params), INTENT(IN)            :: params
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      LOGICAL, INTENT(IN)                                :: calculate_forces, need_xc_gradient, &
                                                            molecular_virial, &
                                                            molecular_virial_debug, have_zeta
      REAL(KIND=dp), INTENT(IN)                          :: force_numerical_gradient_dx, &
                                                            molecular_virial_debug_dx
      TYPE(cp_gauxc_status_type), INTENT(INOUT)          :: status

      IF (have_zeta) THEN
         CALL gauxc_compute_xc( &
            cache%integrator, cache%density_scalar, cache%xc_result, cache%density_zeta, &
            params%nspins, status, model=TRIM(params%model_eval_name))
      ELSE
         CALL gauxc_compute_xc( &
            cache%integrator, cache%density_scalar, cache%xc_result, nspins=params%nspins, &
            status=status, model=TRIM(params%model_eval_name))
      END IF
      CALL gauxc_check_status(status)

      IF (.NOT. need_xc_gradient) RETURN

      IF (params%use_fd_gradient) THEN
         IF (have_zeta) THEN
            CALL gauxc_xc_gradient_fd( &
               particle_set, qs_kind_set, cache%density_scalar, params%nspins, &
               params%model_eval_name, params%xc_fun_name, params%grid_type, &
               params%radial_quadrature, params%pruning_scheme, params%lb_exec_space, &
               params%int_exec_space, params%lwd_kernel, params%batch_size, &
               params%device_runtime_fill_fraction, force_numerical_gradient_dx, para_env, &
               cache%xc_gradient%exc_grad, density_zeta=cache%density_zeta)
         ELSE
            CALL gauxc_xc_gradient_fd( &
               particle_set, qs_kind_set, cache%density_scalar, params%nspins, &
               params%model_eval_name, params%xc_fun_name, params%grid_type, &
               params%radial_quadrature, params%pruning_scheme, params%lb_exec_space, &
               params%int_exec_space, params%lwd_kernel, params%batch_size, &
               params%device_runtime_fill_fraction, force_numerical_gradient_dx, para_env, &
               cache%xc_gradient%exc_grad)
         END IF
      ELSE IF (params%use_gradient_self_runtime) THEN
         IF (have_zeta) THEN
            CALL gauxc_compute_xc_gradient( &
               cache%gradient_integrator, cache%density_scalar, cache%xc_gradient, &
               cache%density_zeta, params%nspins, params%natom, status, &
               model=TRIM(params%model_eval_name))
         ELSE
            CALL gauxc_compute_xc_gradient( &
               cache%gradient_integrator, cache%density_scalar, cache%xc_gradient, &
               nspins=params%nspins, natom=params%natom, status=status, &
               model=TRIM(params%model_eval_name))
         END IF
      ELSE
         IF (have_zeta) THEN
            CALL gauxc_compute_xc_gradient( &
               cache%integrator, cache%density_scalar, cache%xc_gradient, cache%density_zeta, &
               params%nspins, params%natom, status, model=TRIM(params%model_eval_name))
         ELSE
            CALL gauxc_compute_xc_gradient( &
               cache%integrator, cache%density_scalar, cache%xc_gradient, &
               nspins=params%nspins, natom=params%natom, status=status, &
               model=TRIM(params%model_eval_name))
         END IF
      END IF
      CALL gauxc_check_status(status)
      CALL assert_gauxc_gradient_replicated(cache%xc_gradient%exc_grad, para_env)

      IF (calculate_forces) THEN
         CALL add_gauxc_gradient_to_force( &
            cache%xc_gradient%exc_grad, force, atomic_kind_set, para_env)
      END IF
      IF (molecular_virial) THEN
         CALL print_gauxc_molecular_virial(cache%xc_gradient%exc_grad, particle_set, para_env)
      END IF
      IF (molecular_virial_debug) THEN
         IF (have_zeta) THEN
            CALL debug_gauxc_molecular_virial( &
               cache%xc_gradient%exc_grad, particle_set, qs_kind_set, cache%density_scalar, &
               params%nspins, params%model_eval_name, params%xc_fun_name, params%grid_type, &
               params%radial_quadrature, params%pruning_scheme, params%lb_exec_space, &
               params%int_exec_space, params%lwd_kernel, params%batch_size, &
               params%device_runtime_fill_fraction, molecular_virial_debug_dx, para_env, &
               density_zeta=cache%density_zeta)
         ELSE
            CALL debug_gauxc_molecular_virial( &
               cache%xc_gradient%exc_grad, particle_set, qs_kind_set, cache%density_scalar, &
               params%nspins, params%model_eval_name, params%xc_fun_name, params%grid_type, &
               params%radial_quadrature, params%pruning_scheme, params%lb_exec_space, &
               params%int_exec_space, params%lwd_kernel, params%batch_size, &
               params%device_runtime_fill_fraction, molecular_virial_debug_dx, para_env)
         END IF
      END IF

   END SUBROUTINE evaluate_gauxc_forces

! **************************************************************************************************
!> \brief Convert the evaluated GauXC VXC into the CP2K DBCSR matrix channels.
!> \param cache ...
!> \param rho_ao ...
!> \param matrix_vxc ...
!> \param img ...
! **************************************************************************************************
   SUBROUTINE assemble_gauxc_vxc(cache, rho_ao, matrix_vxc, img)
      TYPE(cp_gauxc_cache_type), INTENT(INOUT)           :: cache
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vxc
      INTEGER, INTENT(IN)                                :: img

      IF (SIZE(matrix_vxc) == 1) THEN
         CALL dense_to_dbcsr_into(cache%dense_layout, cache%xc_result%vxc_scalar, &
                                  rho_ao(1, img), matrix_vxc(1))
      ELSE
         CPASSERT(SIZE(matrix_vxc) == 2)
         CALL dense_to_dbcsr_into(cache%dense_layout, cache%xc_result%vxc_scalar, rho_ao(1, img), &
                                  matrix_vxc(1), dense_zeta=cache%xc_result%vxc_zeta, &
                                  z_sign=1.0_dp)
         CALL dense_to_dbcsr_into(cache%dense_layout, cache%xc_result%vxc_scalar, rho_ao(2, img), &
                                  matrix_vxc(2), dense_zeta=cache%xc_result%vxc_zeta, &
                                  z_sign=-1.0_dp)
      END IF

   END SUBROUTINE assemble_gauxc_vxc

END MODULE xc_gauxc_functional
