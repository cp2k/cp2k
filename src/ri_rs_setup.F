!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Builds the shared RI-RS grid, AO collocation and fitting representation.
!> \par History
!>      09.2026 created Jan Wilhelm
!>      09.2026 moved code by Ritaj Tyagi from gw_ri_rs_non_periodic.F
!>      09.2026 added routine to compute Z_lP for automatic RI optimization
! **************************************************************************************************
MODULE ri_rs_setup
   USE ai_contraction_sphi,             ONLY: abc_contract_xsmm
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_binary_read, dbcsr_binary_write, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_multiply, dbcsr_put_block, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE gw_auto_ri_types,                ONLY: auto_ri_type
   USE gw_utils_compute_integrals,      ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              max_line_length
   USE libint_2c_3c,                    ONLY: eri_3center
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE ri_rs_fitting,                   ONLY: solve_ri_rs_fit,&
                                              store_Z_lP_columns
   USE ri_rs_grid_setup,                ONLY: setup_ri_rs_grid
   USE ri_rs_matrices,                  ONLY: atomic_basis_at_grid_point,&
                                              print_matrix_occupation
   USE ri_rs_types,                     ONLY: ri_rs_env
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points,&
                                              evaluate_ao_on_points,&
                                              precompute_ri_rs_radii
   USE util,                            ONLY: locate,&
                                              sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_setup'

   PUBLIC :: compute_Z_lP, ensure_ri_rs_representation

CONTAINS

! **************************************************************************************************
!> \brief Populate the RI-RS grid, AO collocation and auxiliary fitting coefficients once.
!>        Periodic systems use the same neighbouring-image convention for phi and Z.
!> \param qs_env electronic-structure environment
!> \param bs_env geometry, RI-RS settings and owned output matrices
! **************************************************************************************************
   SUBROUTINE ensure_ri_rs_representation(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ensure_ri_rs_representation'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: build_start

      CALL timeset(routineN, handle)
      IF (.NOT. bs_env%ri_rs%grid_built) THEN
         build_start = m_walltime()
         CALL get_qs_env(qs_env, atomic_kind_set=bs_env%ri_rs%atomic_kind_set, &
                         cell=bs_env%ri_rs%cell, particle_set=bs_env%ri_rs%particle_set)
         ! AO/RI support radii determine which atoms and grid points enter each local fit.
         IF (.NOT. ALLOCATED(bs_env%ri_rs%radius_ao_per_atom)) THEN
            CALL precompute_ri_rs_radii(bs_env%ri_rs, bs_env%basis_set_AO, bs_env%basis_set_RI, &
                                        bs_env%eps_filter, bs_env%unit_nr)
         END IF
         ! Place atom-centered RI-RS grids at r_l = R_A + r_l^(A), then flatten them into one grid.
         CALL setup_ri_rs_grid(bs_env, bs_env%ri_rs%grid_points)
         ! Collocation matrix ϕ_μl = ϕ_μ(r_l), including neighbouring AO images for Gamma.
         CALL atomic_basis_at_grid_point(bs_env%ri_rs, bs_env%basis_set_AO, &
                                         bs_env%mat_ao_ao%matrix, bs_env%para_env, bs_env%eps_filter, &
                                         bs_env%dbcsr_msg_elem_limit, bs_env%unit_nr, &
                                         periodic_images=bs_env%do_periodic)
         CALL print_matrix_occupation(bs_env%ri_rs%mat_phi_mu_l, 'ϕ_μ(r_l)', &
                                      bs_env%para_env, bs_env%unit_nr, bs_env%dbcsr_msg_elem_limit)
         ! Fit Z_lP using the integration domains and RI basis described in compute_Z_lP below.
         CALL compute_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                           bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP, &
                           periodic_images=bs_env%do_periodic)
         bs_env%ri_rs%grid_built = .TRUE.
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
               'Built RI-RS representation,', ' Execution time', m_walltime() - build_start, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE ensure_ri_rs_representation

! **************************************************************************************************
!> \brief Computes Z_lP for (1) a tabulated RI basis set or (2) an on-the-fly generated RI basis set.
!>
!>        In both cases, the Z_lP coefficients are computed by solving the linear system
!>
!>          Σ_l' D_ll' Z_l'P = d_lP,
!>
!>        with D_ll' and d_lP given by
!>
!>          D_ll' = [Σ_μ ϕ_μ(r_l) ϕ_μ(r_l')]²,
!>
!>          d_lP = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P).
!>
!>        For an on-the-fly generated RI basis set, an RI function φ_P may be a contraction of
!>        Gaussians on neighboring atoms, which requires special computation of d_lP.
!>
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param periodic_images use Gamma image sums and minimum-image integration spheres
! **************************************************************************************************
   SUBROUTINE compute_Z_lP(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP, periodic_images)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP
      LOGICAL, INTENT(IN), OPTIONAL                      :: periodic_images

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Z_lP'

      INTEGER                                            :: handle
      LOGICAL                                            :: sum_images
      REAL(KIND=dp)                                      :: fitting_stats(2)

      CALL timeset(routineN, handle)
      sum_images = .FALSE.
      IF (PRESENT(periodic_images)) sum_images = periodic_images

      bs_env%ri_rs%fit_max_grid = 0
      bs_env%ri_rs%fit_peak_bytes = 0.0_dp
      IF (bs_env%auto_ri%enabled) THEN
         CALL compute_Z_lP_auto_ri(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP)
         ! Localized AUTO-RI fits can be owned by different process columns.
         fitting_stats = [REAL(bs_env%ri_rs%fit_max_grid, dp), bs_env%ri_rs%fit_peak_bytes]
         CALL bs_env%para_env%max(fitting_stats)
         bs_env%ri_rs%fit_max_grid = NINT(fitting_stats(1))
         bs_env%ri_rs%fit_peak_bytes = fitting_stats(2)
      ELSE
         CALL compute_Z_lP_standard(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP, sum_images)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_Z_lP

! **************************************************************************************************
!> \brief Computes the RI-RS fitting coefficients Z_lP by solving, independently for every RI
!>        atom P, a Jacobi-conditioned, Tikhonov-regularized linear system restricted to the
!>        grid points r_l inside P's integration sphere |r_l - R_P| <= cutoff_ri(P):
!>
!>        1.  D_ll'  = [ Σ_μ ϕ_μ(r_l) ϕ_μ(r_l') ]²
!>        2.  d_l    = 1 / sqrt(max(D_ll, jacobi_floor)) (Jacobi conditioning vector)
!>        3.  D'_ll' = d_l D_ll' d_l' + λ δ_ll'         (λ = TIKHONOV regularization)
!>        4.  d_lP   = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P)
!>        5.  Σ_l' D'_ll' Z'_l'P = d_l d_lP             (solve linear system for Z'_lP)
!>        6.  Z_lP   = d_l Z'_lP                        (undo the conditioning)
!>
!>        See Duchemin and Blase, arXiv:1902.00066, Sec. II.3 for the least-squares construction
!>        and balancing, and solve_ri_rs_fit for the common floor and regularization convention.
!>
!>        The Gamma path uses image-summed AO values and three-center integrals, with RI
!>        atom P fixed in the reference cell. Molecular and Gamma sphere cutoffs are retained.
!>
!>        Work is distributed over atoms in two phases (planned by classify_z_lp_atoms and
!>        lpt_assign_atoms): Phase A solves "small" atoms with single-rank LAPACK
!>        (dpotrf/dpotrs); Phase B solves "big" atoms, whose dense matrix D'_ll' would exceed one
!>        rank's memory, with ScaLAPACK (pdpotrf/pdpotrs) over rank subgroups of size G.
!>        The solved Z columns are scattered into the sparse global mat_Z_lP.
!>        If a Z_lP restart file exists, it is read instead and the solve is skipped entirely.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param periodic_images use Gamma image sums and minimum-image integration spheres
! **************************************************************************************************
   SUBROUTINE compute_Z_lP_standard(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP, periodic_images)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP
      LOGICAL, INTENT(IN)                                :: periodic_images

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'compute_Z_lP_standard'

      INTEGER :: atom_j_mepos, atom_j_stride, atom_P, G, handle, i_blk, iatom, idx, iphase, j, &
         max_ao_size, my_group, n_big, n_done, n_groups, n_loc_ri, n_local_grid, n_my_atoms, &
         n_small, natom, npcol_phi, num_grid_chunks, phase_hi
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ao_atoms, ao_col_start, big_list, local_grid_idx, &
         my_atoms_A, my_atoms_B, n_local_grid_atom, row_offset, small_list
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_ri, r_blk_sizes, ri_blk_sizes, &
                                                            row_dist_grid
      LOGICAL                                            :: do_scatter, use_dist
      LOGICAL, ALLOCATABLE                               :: ao_nonzero(:, :)
      REAL(KIND=dp)                                      :: cutoff_ri, item_start_time, r_c, t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoff_ri_per_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_local, phi_local
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist_phi, dist_Z
      TYPE(gw_3c_ctx_type)                               :: ctx_3c
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, input=input, &
                      qs_kind_set=qs_kind_set, cell=cell)

      NULLIFY (para_env_sub, blacs_env_sub)

      natom = bs_env%n_atom

      ! =========================================================================
      ! 1. SETUP DBCSR TOPOLOGY & EXACT OFFSETS
      !    mat_Z_lP inherits the grid row blocking (and row distribution) of
      !    mat_phi_mu_l; its columns are one block per RI atom.
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, row_blk_size=r_blk_sizes, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, npcols=npcol_phi)

      num_grid_chunks = SIZE(r_blk_sizes)

      ALLOCATE (row_offset(num_grid_chunks))
      row_offset(1) = 0
      DO i_blk = 2, num_grid_chunks
         row_offset(i_blk) = row_offset(i_blk - 1) + r_blk_sizes(i_blk - 1)
      END DO

      ALLOCATE (ri_blk_sizes(natom), col_dist_ri(natom))
      DO iatom = 1, natom
         ri_blk_sizes(iatom) = &
            bs_env%i_RI_end_from_atom(iatom) - bs_env%i_RI_start_from_atom(iatom) + 1
         col_dist_ri(iatom) = MOD(iatom - 1, npcol_phi)
      END DO

      CALL dbcsr_distribution_new(dist_Z, template=dist_phi, &
                                  row_dist=row_dist_grid, col_dist=col_dist_ri)

      IF (bs_env%ri_rs%Z_lP_exists) THEN
         CALL dbcsr_binary_read(filepath=TRIM(bs_env%prefix)//"Z_lP.matrix", &
                                distribution=dist_Z, &
                                matrix_new=mat_Z_lP)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Read Z_lP from file ', ' Execution time', m_walltime() - t1, ' s'
            ! The grid rows use Morton order (spatial_atom_order). A Z_lP.matrix
            ! written with another grid ordering would be silently read into
            ! the current row order. Delete stale Z_lP.matrix files and recompute if in doubt.
            WRITE (bs_env%unit_nr, '(T2,A)') &
               '*** NOTE: Z_lP restart must match the current (spatial) grid row ordering ***'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
      ELSE

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(A)') ' '
            WRITE (bs_env%unit_nr, '(T2,A)') 'Started computing Z_lP'
            CALL m_flush(bs_env%unit_nr)
         END IF
         CALL dbcsr_create(mat_Z_lP, name="mat_Z_lP", dist=dist_Z, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=r_blk_sizes, col_blk_size=ri_blk_sizes)

         ! Largest per-atom AO block, needed to size the 3c-integral work buffers.
         max_ao_size = 0
         DO j = 1, bs_env%n_atom
            max_ao_size = MAX(max_ao_size, &
                              bs_env%i_ao_end_from_atom(j) - &
                              bs_env%i_ao_start_from_atom(j) + 1)
         END DO

         ! Keep each path's integration domain: r_c + r_RI(P) for molecules,
         ! r_c + r_AO(P) for the Gamma image-summed approximation.
         ALLOCATE (cutoff_ri_per_atom(natom))

         r_c = bs_env%ri_metric%cutoff_radius
         DO iatom = 1, natom
            IF (periodic_images) THEN
               cutoff_ri_per_atom(iatom) = r_c + bs_env%ri_rs%radius_ao_per_atom(iatom)
            ELSE
               cutoff_ri_per_atom(iatom) = r_c + bs_env%ri_rs%radius_ri_per_atom(iatom)
            END IF
         END DO

         CALL print_sphere_cutoff_table(bs_env, cutoff_ri_per_atom)

         ! =========================================================================
         ! 2. PER-ATOM SOLVER CLASSIFICATION
         !    Split the atoms into "small" (single-rank LAPACK, Phase A) and "big"
         !    (distributed ScaLAPACK over subgroups of G ranks, Phase B) by comparing
         !    each atom's estimated solve peak memory against the measured budget.
         ! =========================================================================
         CALL classify_z_lp_atoms(bs_env, ri_rs_grid_points, cutoff_ri_per_atom, ri_blk_sizes, &
                                  n_local_grid_atom, small_list, n_small, big_list, n_big, G, periodic_images)

         ! LPT scheduling: sort the atoms of each phase by estimated solve cost
         ! (n_local_grid^3, Cholesky-dominated) and greedily assign to the least-loaded
         ! rank (Phase A) / subgroup (Phase B).
         CALL lpt_assign_atoms(small_list, n_small, n_local_grid_atom, para_env%num_pe, &
                               para_env%mepos, my_atoms_A)
         IF (n_big > 0) THEN
            n_groups = para_env%num_pe/G
            my_group = MIN(para_env%mepos/G, n_groups - 1)
            CALL lpt_assign_atoms(big_list, n_big, n_local_grid_atom, n_groups, my_group, &
                                  my_atoms_B)
         ELSE
            ALLOCATE (my_atoms_B(0))
         END IF

         ! Atoms this rank will process across both phases for rank-0 progress
         n_my_atoms = SIZE(my_atoms_A) + SIZE(my_atoms_B)
         n_done = 0

         ! Shared context for the three-center integrals (μν|P) of the RHS build
         CALL gw_3c_ctx_create(ctx_3c, bs_env, bs_env%ri_metric, &
                               basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                               basis_i=bs_env%basis_set_RI)

         ! =========================================================================
         ! 3. TWO-PHASE LOOP OVER ATOMS
         !    Phase A processes the "small" atoms with the single-rank BLAS path
         !    Phase B processes the "big" atoms with the distributed ScaLAPACK path
         !    over rank subgroups of size G. phi_local for each atom's cutoff sphere
         !    is built on the fly to avoid replicating a global grid x AO matrix.
         ! =========================================================================
         DO iphase = 1, 2
            IF (iphase == 1) THEN
               use_dist = .FALSE.
               atom_j_mepos = 0
               atom_j_stride = 1
               phase_hi = SIZE(my_atoms_A)
            ELSE
               IF (n_big == 0) CYCLE
               use_dist = .TRUE.
               n_groups = para_env%num_pe/G
               my_group = MIN(para_env%mepos/G, n_groups - 1)
               ALLOCATE (para_env_sub)
               CALL para_env_sub%from_split(para_env, my_group)
               CALL cp_blacs_env_create(blacs_env=blacs_env_sub, para_env=para_env_sub)
               atom_j_mepos = para_env_sub%mepos
               atom_j_stride = para_env_sub%num_pe
               ! All ranks of a subgroup share my_group, hence the identical my_atoms_B list
               ! (the per-atom ScaLAPACK solve is collective over the subgroup).
               phase_hi = SIZE(my_atoms_B)
            END IF

            DO idx = 1, phase_hi
               item_start_time = m_walltime()
               IF (iphase == 1) THEN
                  atom_P = my_atoms_A(idx)
               ELSE
                  atom_P = my_atoms_B(idx)
               END IF

               n_loc_ri = ri_blk_sizes(atom_P)
               cutoff_ri = cutoff_ri_per_atom(atom_P)

               ! ---------------------------------------------------------------------
               ! A. Sphere-local AO matrix ϕ_μ(r_l): select the grid points with
               !    |r_l - R_P| <= cutoff_ri(P), using minimum-image distances for Gamma.
               !    The molecular path drops rows below EPS_FILTER; Gamma keeps all sphere rows.
               ! ---------------------------------------------------------------------
               CALL build_phi_on_sphere(bs_env, qs_kind_set, &
                                        ri_rs_grid_points, atom_P, cutoff_ri, &
                                        local_grid_idx, n_local_grid, phi_local, &
                                        ao_atoms, ao_col_start, periodic_images=periodic_images, ao_nonzero=ao_nonzero)

               ! ---------------------------------------------------------------------
               ! B. Right-hand side  d_lP = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P)
               ! ---------------------------------------------------------------------
               ALLOCATE (d_lp_local(n_local_grid, n_loc_ri))
               d_lp_local = 0.0_dp

               CALL compute_d_lp(bs_env, ctx_3c, phi_local, ao_atoms, ao_col_start, ao_nonzero, d_lp_local, n_local_grid, &
                                 n_loc_ri, atom_P, max_ao_size, atom_j_mepos, atom_j_stride, periodic_images)

               DEALLOCATE (ao_atoms, ao_col_start, ao_nonzero)

               ! Reduce per-subgroup-rank partials into the replicated d_lp_local.
               ! Skipped for BLAS path: each rank has the full sum locally.
               IF (use_dist) THEN
                  CALL para_env_sub%sum(d_lp_local)
               END IF

               ! ---------------------------------------------------------------------
               ! C. Solve the Jacobi-scaled, regularized fit with LAPACK or ScaLAPACK.
               ! ---------------------------------------------------------------------
               IF (use_dist) THEN
                  CALL solve_ri_rs_fit(phi_local, d_lp_local, bs_env%ri_rs%tikhonov, &
                                       para_env_sub, blacs_env_sub)
               ELSE
                  CALL solve_ri_rs_fit(phi_local, d_lp_local, bs_env%ri_rs%tikhonov)
               END IF

               ! ---------------------------------------------------------------------
               ! D. Scatter the solved Z columns back into the global sparse mat_Z_lP.
               ! ---------------------------------------------------------------------
               do_scatter = .TRUE.
               IF (use_dist) do_scatter = (para_env_sub%mepos == 0)
               IF (do_scatter) THEN
                  CALL store_Z_lP_columns(mat_Z_lP, d_lp_local, local_grid_idx, n_local_grid, &
                                          n_loc_ri, atom_P, r_blk_sizes, row_offset, &
                                          bs_env%eps_filter)
               END IF

               DEALLOCATE (d_lp_local)
               DEALLOCATE (local_grid_idx, phi_local)

               ! Report each atom completed by the printing rank. No inter-rank communication is
               ! needed; the final message below is printed only after the global synchronization.
               n_done = n_done + 1
               CALL print_Z_lP_progress(bs_env, n_done, n_my_atoms, &
                                        m_walltime() - item_start_time)
            END DO   ! idx: atoms of this phase owned by this rank / subgroup

            ! Tear down the Phase-B subgroup (all ranks created it collectively).
            IF (iphase == 2) THEN
               CALL cp_blacs_env_release(blacs_env_sub)
               CALL para_env_sub%free()
               DEALLOCATE (para_env_sub)
            END IF
         END DO   ! iphase

         DEALLOCATE (cutoff_ri_per_atom)
         DEALLOCATE (small_list, big_list)

         CALL gw_3c_ctx_release(ctx_3c)

         CALL dbcsr_finalize(mat_Z_lP)

         CALL para_env%sync()
         CALL print_Z_lP_progress(bs_env, natom, natom, m_walltime() - t1, &
                                  all_mpi_ranks=.TRUE.)

         logger => cp_get_default_logger()

         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN
            CALL dbcsr_binary_write(matrix=mat_Z_lP, filepath=TRIM(bs_env%prefix)//"Z_lP.matrix")
         END IF

      END IF

      DEALLOCATE (row_offset, ri_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)

      CALL timestop(handle)

   END SUBROUTINE compute_Z_lP_standard

! **************************************************************************************************
!> \brief Prints the per-kind maximum RI-RS integration-sphere cutoff table.
!> \param bs_env ...
!> \param cutoff_ri_per_atom ...
! **************************************************************************************************
   SUBROUTINE print_sphere_cutoff_table(bs_env, cutoff_ri_per_atom)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cutoff_ri_per_atom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_sphere_cutoff_table'

      INTEGER                                            :: handle, iatom, ikind, nkind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoff_ri_per_kind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      atomic_kind_set => bs_env%ri_rs%atomic_kind_set
      particle_set => bs_env%ri_rs%particle_set

      IF (bs_env%unit_nr <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      nkind = SIZE(atomic_kind_set)
      ALLOCATE (cutoff_ri_per_kind(nkind))
      cutoff_ri_per_kind(:) = 0.0_dp

      DO iatom = 1, bs_env%n_atom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         cutoff_ri_per_kind(ikind) = MAX(cutoff_ri_per_kind(ikind), cutoff_ri_per_atom(iatom))
      END DO

      WRITE (bs_env%unit_nr, '(T2,A)') 'Per-kind maximum RI-RS sphere cutoff (Å):'
      WRITE (bs_env%unit_nr, '(T4,A4,A14)') 'Kind', 'cutoff (Å)'
      DO ikind = 1, nkind
         WRITE (bs_env%unit_nr, '(T4,A4,F14.4)') &
            atomic_kind_set(ikind)%element_symbol, &
            cutoff_ri_per_kind(ikind)*angstrom
      END DO
      WRITE (bs_env%unit_nr, '(A)') ' '

      DEALLOCATE (cutoff_ri_per_kind)

      CALL timestop(handle)

   END SUBROUTINE print_sphere_cutoff_table

! **************************************************************************************************
!> \brief LPT (longest-processing-time) assignment of the Z_lP atoms to workers (MPI ranks in
!>        Phase A, rank subgroups in Phase B): sort by estimated solve cost n_local_grid^3
!>        (the per-atom Cholesky dominates; the n^2 assembly terms order the atoms the same
!>        way) and greedily give each atom to the least-loaded worker.
!> \param atom_list ...
!> \param n_atoms ...
!> \param n_local_grid_atom ...
!> \param n_workers ...
!> \param my_worker ...
!> \param my_atoms ...
! **************************************************************************************************
   SUBROUTINE lpt_assign_atoms(atom_list, n_atoms, n_local_grid_atom, n_workers, my_worker, &
                               my_atoms)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_list
      INTEGER, INTENT(IN)                                :: n_atoms
      INTEGER, DIMENSION(:), INTENT(IN)                  :: n_local_grid_atom
      INTEGER, INTENT(IN)                                :: n_workers, my_worker
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: my_atoms

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lpt_assign_atoms'

      INTEGER                                            :: handle, i, iw, n_mine, w_min
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mine_tmp, perm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cost, load

      CALL timeset(routineN, handle)

      IF (n_atoms <= 0) THEN
         ALLOCATE (my_atoms(0))
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (cost(n_atoms), perm(n_atoms), mine_tmp(n_atoms), load(n_workers))
      DO i = 1, n_atoms
         cost(i) = REAL(n_local_grid_atom(atom_list(i)), dp)**3
      END DO
      CALL sort(cost, n_atoms, perm)   ! ascending; walk backwards for largest-first

      load(:) = 0.0_dp
      n_mine = 0
      DO i = n_atoms, 1, -1
         w_min = 1
         DO iw = 2, n_workers
            IF (load(iw) < load(w_min)) w_min = iw
         END DO
         load(w_min) = load(w_min) + cost(i)
         IF (w_min - 1 == my_worker) THEN
            n_mine = n_mine + 1
            mine_tmp(n_mine) = atom_list(perm(i))
         END IF
      END DO

      ALLOCATE (my_atoms(n_mine))
      my_atoms(:) = mine_tmp(1:n_mine)

      CALL timestop(handle)

   END SUBROUTINE lpt_assign_atoms

! **************************************************************************************************
!> \brief Computes the dense localized right-hand side for one RI atom P,
!>
!>          d_lP = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P).
!>
!>        OMP-threaded over unordered AO-atom pairs: build (mu nu|P) in the reference cell,
!>        or sum (mu_R nu_S|P_0) over nearest images for Gamma. Contract grid-chunked pair
!>        densities into a thread-private d_lp partial, then reduce the partials into d_lp.
!>        The shared integral context caches Libint/spherical-transform tables; each thread
!>        owns its integral workspace. The kind-radius triangle screen in
!>        build_3c_integral_block_ctx skips negligible image triples via the screened output.
!> \param bs_env ...
!> \param ctx shared three-center integral context
!> \param phi_val AO values on local-sphere grid rows and compact active-AO columns
!> \param ao_atoms global indices of the active AO atoms
!> \param ao_col_start local phi column starts, including the terminal column offset
!> \param ao_nonzero exact nonzero AO support on the retained grid rows
!> \param d_lp RHS (n_grid_total × n_loc_ri), zeroed by the caller and accumulated here
!> \param n_grid_total number of retained local-sphere grid rows
!> \param n_loc_ri number of RI functions on atom_P
!> \param atom_P RI atom fixed in the reference cell (0,0,0)
!> \param max_ao_size ...
!> \param atom_j_mepos ...
!> \param atom_j_stride ...
!> \param periodic_images sum nearest AO images with RI atom P fixed in the reference cell
! **************************************************************************************************
   SUBROUTINE compute_d_lp(bs_env, ctx, phi_val, ao_atoms, ao_col_start, ao_nonzero, d_lp, n_grid_total, n_loc_ri, atom_P, &
                           max_ao_size, atom_j_mepos, atom_j_stride, periodic_images)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_val
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_atoms, ao_col_start
      LOGICAL, DIMENSION(:, :), INTENT(IN)               :: ao_nonzero
      INTEGER, INTENT(IN)                                :: n_grid_total, n_loc_ri
      REAL(KIND=dp), INTENT(INOUT)                       :: d_lp(n_grid_total, n_loc_ri)
      INTEGER, INTENT(IN)                                :: atom_P, max_ao_size, atom_j_mepos, &
                                                            atom_j_stride
      LOGICAL, INTENT(IN)                                :: periodic_images

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_d_lp'
      INTEGER, PARAMETER                                 :: grid_chunk = 1024

      INTEGER :: atom_j, atom_k, c, cell_R_vec(3), cell_S_vec(3), handle, handle2, &
         image_extent(3), ix_R, ix_S, iy_R, iy_S, iz_R, iz_S, j, jj, jk_idx, jsize, jstart, k, kk, &
         ksize, kstart, l, l0, n_grid_pair, point, ri
      INTEGER, ALLOCATABLE                               :: grid_index(:)
      LOGICAL                                            :: any_kept, screened
      REAL(KIND=dp)                                      :: pair_factor
      REAL(KIND=dp), ALLOCATABLE                         :: grid_result(:, :), int_3c_prv(:, :, :)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_prv, int_2d_prv, rho_chunk
      TYPE(gw_3c_ws_type)                                :: ws

      CALL timeset(routineN, handle)
      image_extent = 0
      IF (periodic_images) image_extent = bs_env%ri_rs%cell%perd

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, ctx, phi_val, ao_atoms, ao_col_start, ao_nonzero, d_lp, &
      !$OMP        n_grid_total, n_loc_ri, atom_P, max_ao_size, &
      !$OMP        atom_j_mepos, atom_j_stride, image_extent) &
      !$OMP PRIVATE(grid_index, grid_result, n_grid_pair, point, pair_factor, &
      !$OMP         jj, kk, atom_j, atom_k, c, handle2, j, jk_idx, jsize, jstart, k, ksize, kstart, &
      !$OMP         l, l0, ri, screened, d_lp_prv, int_2d_prv, rho_chunk, int_3c_prv, ws, &
      !$OMP         any_kept, ix_R, iy_R, iz_R, ix_S, iy_S, iz_S, cell_R_vec, cell_S_vec)

      ! The integral context is shared; all scratch arrays and the accumulating RHS are thread-private.
      CALL gw_3c_ws_create(ws, ctx)
      ALLOCATE (int_3c_prv(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (int_2d_prv(max_ao_size*max_ao_size, n_loc_ri))
      ALLOCATE (rho_chunk(grid_chunk, max_ao_size*max_ao_size))
      ALLOCATE (d_lp_prv(n_grid_total, n_loc_ri))
      ALLOCATE (grid_index(n_grid_total), grid_result(grid_chunk, n_loc_ri))
      d_lp_prv(:, :) = 0.0_dp

      ! MPI assigns each unordered pair through its first atom; OpenMP divides those atoms.
      ! Skip only exact-zero pair-grid support, without a new screening threshold.
      !$OMP DO SCHEDULE(DYNAMIC)
      DO jj = 1, SIZE(ao_atoms)
         atom_j = ao_atoms(jj)
         IF (MOD(atom_j - 1, atom_j_stride) /= atom_j_mepos) CYCLE
         DO kk = jj, SIZE(ao_atoms)
            atom_k = ao_atoms(kk)
            jstart = ao_col_start(jj)
            kstart = ao_col_start(kk)
            jsize = ao_col_start(jj + 1) - jstart
            ksize = ao_col_start(kk + 1) - kstart

            n_grid_pair = 0
            DO point = 1, n_grid_total
               IF (.NOT. (ao_nonzero(point, jj) .AND. ao_nonzero(point, kk))) CYCLE
               n_grid_pair = n_grid_pair + 1
               grid_index(n_grid_pair) = point
            END DO
            IF (n_grid_pair == 0) CYCLE
            ! (μν|P) = (νμ|P): distinct atom pairs contribute twice.
            pair_factor = 1.0_dp
            IF (atom_j /= atom_k) pair_factor = 2.0_dp

            ! For Gamma, sum (mu_R nu_S|P_0); the molecular path has only R=S=0.
            ! AO exchange symmetry applies to the complete image sum, so the same
            ! unordered atom pairs and pair_factor are used in both paths.
            ! The integral builder accumulates in place, including across image triples.
            ! Its kind-radius triangle screen rejects image triples beyond the integral cutoff.
            int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp
            any_kept = .FALSE.
            DO ix_R = -image_extent(1), image_extent(1)
            DO iy_R = -image_extent(2), image_extent(2)
            DO iz_R = -image_extent(3), image_extent(3)
               cell_R_vec = [ix_R, iy_R, iz_R]
               DO ix_S = -image_extent(1), image_extent(1)
               DO iy_S = -image_extent(2), image_extent(2)
               DO iz_S = -image_extent(3), image_extent(3)
                  cell_S_vec = [ix_S, iy_S, iz_S]
                  CALL build_3c_integral_block_ctx(int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri), &
                                                   ctx, ws, atom_j=atom_j, atom_k=atom_k, atom_i=atom_P, &
                                                   cell_j=cell_R_vec, cell_k=cell_S_vec, screened=screened)
                  IF (screened) CYCLE
                  any_kept = .TRUE.
               END DO
               END DO
               END DO
            END DO
            END DO
            END DO
            IF (.NOT. any_kept) CYCLE

            ! Flatten (μν|P) to a matrix with combined AO index (ν-1)*jsize + μ for BLAS.
            DO ri = 1, n_loc_ri
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     int_2d_prv(jk_idx, ri) = int_3c_prv(j, k, ri)
                  END DO
               END DO
            END DO

            ! Pair density ρ(l, μν) = ϕ_μ(r_l) ϕ_ν(r_l) in grid chunks, contracted on the fly:
            ! d_lP += Σ_μν ρ(l, μν) (μν|P), with the atom-pair multiplicity included in pair_factor.
            DO l0 = 1, n_grid_pair, grid_chunk
               c = MIN(grid_chunk, n_grid_pair - l0 + 1)
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     DO l = 1, c
                        point = grid_index(l0 + l - 1)
                        rho_chunk(l, jk_idx) = phi_val(point, jstart + j - 1)* &
                                               phi_val(point, kstart + k - 1)
                     END DO
                  END DO
               END DO
               CALL timeset(routineN//"_dgemm", handle2)
               CALL dgemm("N", "N", c, n_loc_ri, jsize*ksize, &
                          pair_factor, rho_chunk, grid_chunk, &
                          int_2d_prv, max_ao_size*max_ao_size, &
                          0.0_dp, grid_result, grid_chunk)
               DO ri = 1, n_loc_ri
                  DO l = 1, c
                     point = grid_index(l0 + l - 1)
                     d_lp_prv(point, ri) = d_lp_prv(point, ri) + grid_result(l, ri)
                  END DO
               END DO
               CALL timestop(handle2)
            END DO
         END DO
      END DO
      !$OMP END DO

      ! Reduce once per thread; distributed fits subsequently sum the rank partials in the subgroup.
      !$OMP CRITICAL (compute_d_lp_reduce)
      d_lp(1:n_grid_total, 1:n_loc_ri) = d_lp(1:n_grid_total, 1:n_loc_ri) + &
                                         d_lp_prv(1:n_grid_total, 1:n_loc_ri)
      !$OMP END CRITICAL (compute_d_lp_reduce)

      DEALLOCATE (int_3c_prv, int_2d_prv, rho_chunk, d_lp_prv, grid_index, grid_result)
      CALL gw_3c_ws_release(ws)

      !$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE compute_d_lp

! **************************************************************************************************
!> \brief Computes the RI-RS matrix Z_lP for an automatically optimized RI basis.
!>
!>        The fitting coefficients Z_lP satisfy
!>
!>          Σ_l' D_ll' Z_l'P = d_lP,
!>          D_ll' = [Σ_μ ϕ_μ(r_l) ϕ_μ(r_l')]²,
!>          d_lP = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P).
!>
!>        For the standard RI-RS fit, P is an atom-centered RI function. Here an optimized
!>        function may contain reference functions on neighboring atoms,
!>
!>          φ_p(r) = Σ_A Σ_{P∈A} U_Pp^A φ_P^A(r),
!>
!>        so all atomic contributions to d_lp must be accumulated before solving for Z_lp.
!>        If every optimized column uses the complete RI-RS grid, all columns are solved in one
!>        system. Otherwise, each atom-blocked column set is fitted on its own integration sphere;
!>        grid points outside that sphere are excluded from its fitting equations.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************
   SUBROUTINE compute_Z_lP_auto_ri(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'compute_Z_lP_auto_ri'

      INTEGER :: AB_block, atom_A, atom_B, block_size_b, column_first, first_p_AB, fit_atom, &
         handle, last_atom, max_ao_size, max_nRI_ref, mypcol, myprow, n_active, n_candidate, &
         n_done, n_my_atoms, n_to_a, natom, ncol, ngrid, nri, nRI_ref_A, nRI_ref_B, output_offset, &
         ri_atom, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_atoms, ao_col_start, candidate_order, &
                                                            candidates, local_grid_idx, ri_atoms, &
                                                            row_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_ri, ri_blk_sizes, &
                                                            row_dist_grid, row_size_grid
      LOGICAL                                            :: AB_block_local, common_grid_available, &
                                                            have_fitted_columns, &
                                                            reuse_atomic_integrals
      LOGICAL, ALLOCATABLE                               :: ao_nonzero(:, :)
      REAL(KIND=dp)                                      :: cutoff_ri, item_start_time, &
                                                            mem_avail_GB, previous_fit_peak, r_c, &
                                                            t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_all, d_lp_local, phi_local, U_Pp_AB
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: U_active
      REAL(KIND=dp), DIMENSION(3)                        :: center
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist_Z
      TYPE(dbcsr_type)                                   :: mat_rhs
      TYPE(gw_3c_ctx_type)                               :: ctx_3c
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_col
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      t1 = m_walltime()
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
         WRITE (bs_env%unit_nr, '(T2,A)') 'Started computing Z_lP'
         CALL m_flush(bs_env%unit_nr)
      END IF
      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, input=input, &
                      qs_kind_set=qs_kind_set, cell=cell)
      NULLIFY (para_env_col)

      natom = bs_env%n_atom
      CPASSERT(bs_env%auto_ri%AB_block_count > 0)
      CPASSERT(SIZE(particle_set) == natom)

      CALL prepare_Z_lP_auto_ri(bs_env, mat_phi_mu_l, mat_Z_lP, dist_Z, row_size_grid, &
                                row_dist_grid, col_dist_ri, ri_blk_sizes, row_offset, &
                                myprow, mypcol, max_ao_size)

      CALL gw_3c_ctx_create(ctx_3c, bs_env, bs_env%ri_metric, &
                            basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                            basis_i=bs_env%basis_set_RI)

      ! Check whether one grid can represent d_lp for every optimized atomic column block.
      CALL common_Z_lP_grid_available(bs_env, ri_rs_grid_points, &
                                      common_grid_available, have_fitted_columns)

      IF (common_grid_available .AND. have_fitted_columns) THEN
         CALL mp_mem_avail_per_rank_GB(para_env, mem_avail_GB)
         ! Evaluate ϕ_μ(r_l) once on the common grid.
         CALL build_phi_on_complete_grid(bs_env, qs_kind_set, &
                                         ri_rs_grid_points, local_grid_idx, ngrid, &
                                         phi_local, ao_atoms, ao_col_start, ao_nonzero)
         nri = SUM(ri_blk_sizes)
         CALL record_fit_memory(bs_env%ri_rs, ngrid, SIZE(phi_local, 2), nri)
         CALL print_fit_memory_estimate(bs_env, mem_avail_GB)
         ALLOCATE (d_lp_all(ngrid, nri), source=0.0_dp)
         ! d_lp = Σ_A Σ_{P∈A} d_lP U_Pp for every optimized RI function q.
         CALL compute_d_lp_auto_ri_batch(bs_env, ctx_3c, phi_local, ao_atoms, ao_col_start, ao_nonzero, d_lp_all, ngrid, &
                                         max_ao_size, para_env%mepos, para_env%num_pe)
         DEALLOCATE (ao_atoms, ao_col_start, ao_nonzero)
         CALL para_env%sum(d_lp_all)
         CALL solve_ri_rs_fit(phi_local, d_lp_all, bs_env%ri_rs%tikhonov)

         DO AB_block = 1, bs_env%auto_ri%AB_block_count
            atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
            atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
            AB_block_local = col_dist_ri(atom_A) == mypcol
            IF (atom_B /= atom_A) AB_block_local = AB_block_local .OR. col_dist_ri(atom_B) == mypcol
            IF (.NOT. AB_block_local) CYCLE
            ncol = bs_env%auto_ri%AB_size_opt_RI(AB_block)
            n_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
            IF (n_to_a > 0 .AND. col_dist_ri(atom_A) == mypcol) THEN
               output_offset = SUM(ri_blk_sizes(:atom_A - 1)) + &
                               bs_env%auto_ri%AB_first_p_A(AB_block) - 1
               CALL add_Z_lP_columns(mat_Z_lP, &
                                     d_lp_all(:, output_offset + 1:output_offset + n_to_a), &
                                     local_grid_idx, ngrid, atom_A, &
                                     bs_env%auto_ri%AB_first_p_A(AB_block), &
                                     ri_blk_sizes(atom_A), row_size_grid, row_offset, &
                                     row_dist_grid, &
                                     myprow, bs_env%eps_filter)
            END IF
            IF (n_to_a < ncol) THEN
               CPASSERT(atom_B /= atom_A)
               IF (col_dist_ri(atom_B) == mypcol) THEN
                  block_size_b = ri_blk_sizes(atom_B)
                  output_offset = SUM(ri_blk_sizes(:atom_B - 1)) + &
                                  bs_env%auto_ri%AB_first_p_B(AB_block) - 1
                  CALL add_Z_lP_columns(mat_Z_lP, &
                                        d_lp_all(:, output_offset + 1: &
                                                 output_offset + ncol - n_to_a), &
                                        local_grid_idx, ngrid, atom_B, &
                                        bs_env%auto_ri%AB_first_p_B(AB_block), block_size_b, &
                                        row_size_grid, row_offset, row_dist_grid, myprow, &
                                        bs_env%eps_filter)
               END IF
            END IF
         END DO
         DEALLOCATE (d_lp_all, local_grid_idx, phi_local)
      ELSE
         reuse_atomic_integrals = &
            bs_env%ri_rs%cutoff_radius_ri_ao <= 0.0_dp .OR. &
            bs_env%ri_rs%cutoff_radius_ri_ao <= &
            MINVAL(bs_env%ri_rs%radius_ao_per_atom)
         IF (reuse_atomic_integrals) THEN
            ! Form d_lp = Σ_A Σ_{P∈A} d_lP U_Pp before fitting each atomic column block.
            CALL compute_auto_ri_d_lp(qs_env, bs_env, ctx_3c, &
                                      ri_rs_grid_points, mat_phi_mu_l, mat_rhs, &
                                      max_ao_size)
            CALL dbcsr_release(mat_Z_lP)
            CALL dbcsr_create(mat_Z_lP, name='mat_Z_lP localized AA/AB', dist=dist_Z, &
                              matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_size_grid, &
                              col_blk_size=ri_blk_sizes)
            ! Solve Σ_l' D_ll' Z_l'q = d_lp on each atomic fitting grid.
            CALL fit_auto_ri_z_lp(qs_env, bs_env, ri_rs_grid_points, mat_rhs, mat_Z_lP)
            CALL dbcsr_release(mat_rhs)
         ELSE
            ! All ranks query memory before process columns enter independent fitting loops.
            CALL mp_mem_avail_per_rank_GB(para_env, mem_avail_GB)
            ALLOCATE (para_env_col)
            CALL para_env_col%from_split(para_env, mypcol)
            n_my_atoms = 0
            DO fit_atom = 1, natom
               IF (col_dist_ri(fit_atom) /= mypcol) CYCLE
               IF (ri_blk_sizes(fit_atom) == 0) CYCLE
               n_my_atoms = n_my_atoms + 1
            END DO
            n_done = 0
            ALLOCATE (candidates(2*bs_env%auto_ri%AB_block_count), &
                      candidate_order(2*bs_env%auto_ri%AB_block_count))
            DO fit_atom = 1, natom
               IF (col_dist_ri(fit_atom) /= mypcol) CYCLE
               ncol = ri_blk_sizes(fit_atom)
               IF (ncol == 0) CYCLE
               item_start_time = m_walltime()
               ! Only reference RI atoms contributing to this fitted column block are stored.
               n_candidate = 0
               DO AB_block = 1, bs_env%auto_ri%AB_block_count
                  atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
                  atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
                  n_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
                  IF (.NOT. ((fit_atom == atom_A .AND. n_to_a > 0) .OR. &
                             (fit_atom == atom_B .AND. &
                              n_to_a < bs_env%auto_ri%AB_size_opt_RI(AB_block)))) CYCLE
                  n_candidate = n_candidate + 1
                  candidates(n_candidate) = atom_A
                  IF (atom_B == atom_A) CYCLE
                  n_candidate = n_candidate + 1
                  candidates(n_candidate) = atom_B
               END DO
               CPASSERT(n_candidate > 0)
               CALL sort(candidates, n_candidate, candidate_order)
               n_active = 0
               last_atom = 0
               DO slot = 1, n_candidate
                  ri_atom = candidates(slot)
                  IF (ri_atom == last_atom) CYCLE
                  n_active = n_active + 1
                  candidates(n_active) = ri_atom
                  last_atom = ri_atom
               END DO
               ALLOCATE (ri_atoms(n_active), source=candidates(:n_active))
               max_nRI_ref = 0
               DO slot = 1, n_active
                  max_nRI_ref = MAX(max_nRI_ref, get_ref_RI_size(bs_env, ri_atoms(slot)))
               END DO
               ALLOCATE (U_active(max_nRI_ref, ncol, n_active), source=0.0_dp)

               DO AB_block = 1, bs_env%auto_ri%AB_block_count
                  atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
                  atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
                  n_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
                  IF (fit_atom == atom_A .AND. n_to_a > 0) THEN
                     column_first = bs_env%auto_ri%AB_first_p_A(AB_block)
                     first_p_AB = 1
                     ncol = n_to_a
                  ELSE IF (fit_atom == atom_B .AND. &
                           n_to_a < bs_env%auto_ri%AB_size_opt_RI(AB_block)) THEN
                     column_first = bs_env%auto_ri%AB_first_p_B(AB_block)
                     first_p_AB = n_to_a + 1
                     ncol = bs_env%auto_ri%AB_size_opt_RI(AB_block) - n_to_a
                  ELSE
                     CYCLE
                  END IF
                  CALL get_U_Pp_AB(bs_env%auto_ri, AB_block, U_Pp_AB)

                  slot = locate(ri_atoms, atom_A)
                  nRI_ref_A = get_ref_RI_size(bs_env, atom_A)
                  U_active(1:nRI_ref_A, &
                           column_first:column_first + ncol - 1, slot) = &
                     U_Pp_AB( &
                     1:nRI_ref_A, first_p_AB:first_p_AB + ncol - 1)
                  IF (atom_B /= atom_A) THEN
                     slot = locate(ri_atoms, atom_B)
                     nRI_ref_B = get_ref_RI_size(bs_env, atom_B)
                     U_active(1:nRI_ref_B, &
                              column_first:column_first + ncol - 1, slot) = &
                        U_Pp_AB( &
                        nRI_ref_A + 1:nRI_ref_A + nRI_ref_B, &
                        first_p_AB:first_p_AB + ncol - 1)
                  END IF
                  DEALLOCATE (U_Pp_AB)
               END DO

               center = particle_set(fit_atom)%r
               r_c = bs_env%ri_metric%cutoff_radius
               cutoff_ri = 0.0_dp
               DO slot = 1, n_active
                  ri_atom = ri_atoms(slot)
                  cutoff_ri = MAX(cutoff_ri, &
                                  r_c + bs_env%ri_rs%radius_ri_per_atom(ri_atom) + &
                                  NORM2(center - particle_set(ri_atom)%r))
               END DO
               CALL build_phi_on_sphere(bs_env, qs_kind_set, &
                                        ri_rs_grid_points, fit_atom, cutoff_ri, &
                                        local_grid_idx, ngrid, phi_local, ao_atoms, ao_col_start, &
                                        center=center, ao_nonzero=ao_nonzero)
               ncol = ri_blk_sizes(fit_atom)
               previous_fit_peak = bs_env%ri_rs%fit_peak_bytes
               CALL record_fit_memory(bs_env%ri_rs, ngrid, SIZE(phi_local, 2), ncol)
               IF (bs_env%ri_rs%fit_peak_bytes > previous_fit_peak) THEN
                  CALL print_fit_memory_estimate(bs_env, mem_avail_GB)
               END IF
               ALLOCATE (d_lp_local(ngrid, ncol), source=0.0_dp)
               ! d_lp = Σ_A Σ_{P∈A} d_lP U_Pp on the grid of this atomic column block.
               CALL compute_d_lp_auto_ri(bs_env, ctx_3c, phi_local, ao_atoms, ao_col_start, ao_nonzero, &
                                         d_lp_local, ri_atoms, U_active, max_ao_size, &
                                         para_env_col%mepos, para_env_col%num_pe)
               DEALLOCATE (ao_atoms, ao_col_start, ao_nonzero)
               CALL para_env_col%sum(d_lp_local)
               CALL solve_ri_rs_fit(phi_local, d_lp_local, bs_env%ri_rs%tikhonov)
               CALL add_Z_lP_columns(mat_Z_lP, d_lp_local, local_grid_idx, ngrid, fit_atom, 1, &
                                     ri_blk_sizes(fit_atom), row_size_grid, row_offset, &
                                     row_dist_grid, &
                                     myprow, bs_env%eps_filter)
               DEALLOCATE (d_lp_local, local_grid_idx, phi_local)
               DEALLOCATE (U_active, ri_atoms)
               n_done = n_done + 1
               CALL print_Z_lP_progress(bs_env, n_done, n_my_atoms, &
                                        m_walltime() - item_start_time)
            END DO
            DEALLOCATE (candidates, candidate_order)
            CALL para_env_col%free()
            DEALLOCATE (para_env_col)
         END IF
      END IF
      CALL gw_3c_ctx_release(ctx_3c)

      CALL dbcsr_filter(mat_Z_lP, bs_env%eps_filter)
      CALL dbcsr_finalize(mat_Z_lP)
      CALL para_env%sync()
      CALL print_Z_lP_progress(bs_env, natom, natom, m_walltime() - t1, &
                               all_mpi_ranks=.TRUE.)
      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN
         CALL dbcsr_binary_write(matrix=mat_Z_lP, filepath=TRIM(bs_env%prefix)//'Z_lP.matrix')
      END IF

      DEALLOCATE (row_offset, ri_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)
      CALL timestop(handle)

   END SUBROUTINE compute_Z_lP_auto_ri

! **************************************************************************************************
!> \brief Prints completed Z_lP atoms and execution time from the printing rank.
!> \param bs_env ...
!> \param n_done Number of atoms completed by the printing rank.
!> \param n_total Total number of atoms assigned to the printing rank.
!> \param execution_time Wall-clock time used for the reported work.
!> \param all_mpi_ranks Whether all MPI ranks have completed the calculation.
! **************************************************************************************************
   SUBROUTINE print_Z_lP_progress(bs_env, n_done, n_total, execution_time, all_mpi_ranks)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_done, n_total
      REAL(KIND=dp), INTENT(IN)                          :: execution_time
      LOGICAL, INTENT(IN), OPTIONAL                      :: all_mpi_ranks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_Z_lP_progress'

      INTEGER                                            :: handle
      LOGICAL                                            :: completed

      CALL timeset(routineN, handle)

      IF (bs_env%unit_nr > 0) THEN
         completed = .FALSE.
         IF (PRESENT(all_mpi_ranks)) completed = all_mpi_ranks
         IF (completed) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A,/)') &
               'Computed Z_lP (all MPI ranks) for all atoms,', &
               'Execution time', execution_time, ' s'
         ELSE
            WRITE (bs_env%unit_nr, '(T2,A,I11,A,I3,A,F7.1,A)') &
               'Computed Z_lP (MPI rank 0) for atom', n_done, ' /', n_total, &
               ',    Execution time', execution_time, ' s'
         END IF
         CALL m_flush(bs_env%unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_Z_lP_progress

! **************************************************************************************************
!> \brief Prepares the distributed Z_lP matrix and its atom-blocked column layout. Column block A
!>        contains all optimized functions assigned to atom A. Every local block is reserved once
!>        because several AA/AB contraction blocks may contribute to the same block.
!> \param bs_env ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param dist_Z ...
!> \param row_size_grid ...
!> \param row_dist_grid ...
!> \param col_dist_ri ...
!> \param ri_blk_sizes ...
!> \param row_offset ...
!> \param myprow ...
!> \param mypcol ...
!> \param max_ao_size ...
! **************************************************************************************************
   SUBROUTINE prepare_Z_lP_auto_ri(bs_env, mat_phi_mu_l, mat_Z_lP, dist_Z, row_size_grid, &
                                   row_dist_grid, col_dist_ri, ri_blk_sizes, row_offset, &
                                   myprow, mypcol, max_ao_size)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: dist_Z
      INTEGER, DIMENSION(:), POINTER                     :: row_size_grid, row_dist_grid, &
                                                            col_dist_ri, ri_blk_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: row_offset
      INTEGER, INTENT(OUT)                               :: myprow, mypcol, max_ao_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_Z_lP_auto_ri'

      INTEGER                                            :: handle, i_blk, iatom, npcol
      TYPE(dbcsr_distribution_type)                      :: dist_phi

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_phi_mu_l, row_blk_size=row_size_grid, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, npcols=npcol, &
                                  myprow=myprow, mypcol=mypcol)
      ALLOCATE (row_offset(SIZE(row_size_grid)))
      row_offset(1) = 0
      DO i_blk = 2, SIZE(row_size_grid)
         row_offset(i_blk) = row_offset(i_blk - 1) + row_size_grid(i_blk - 1)
      END DO

      ALLOCATE (ri_blk_sizes(bs_env%n_atom), col_dist_ri(bs_env%n_atom))
      ri_blk_sizes = bs_env%auto_ri%sizes_opt_RI
      DO iatom = 1, bs_env%n_atom
         col_dist_ri(iatom) = MOD(iatom - 1, npcol)
      END DO
      CALL dbcsr_distribution_new(dist_Z, template=dist_phi, row_dist=row_dist_grid, &
                                  col_dist=col_dist_ri)
      CALL dbcsr_create(mat_Z_lP, name='mat_Z_lP localized AA/AB', dist=dist_Z, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_size_grid, &
                        col_blk_size=ri_blk_sizes)
      CALL dbcsr_reserve_all_blocks(mat_Z_lP)
      CALL dbcsr_set(mat_Z_lP, 0.0_dp)

      max_ao_size = 0
      DO iatom = 1, bs_env%n_atom
         max_ao_size = MAX(max_ao_size, bs_env%i_ao_end_from_atom(iatom) - &
                           bs_env%i_ao_start_from_atom(iatom) + 1)
      END DO

      CALL timestop(handle)

   END SUBROUTINE prepare_Z_lP_auto_ri

! **************************************************************************************************
!> \brief Tests whether every fitted atomic block can use the complete RI-RS grid.
!>
!>        For a block centered on A, every grid point must satisfy
!>
!>          |r_l - R_A| <= R_A^fit,
!>
!>        and every atom B whose AOs can contribute must satisfy
!>
!>          |R_B - R_A| <= R_B^AO + R_A^fit.
!>
!>        Pair-midpoint grids do not satisfy this atom-centered criterion in general.
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param common_grid_available ...
!> \param have_fitted_columns ...
! **************************************************************************************************
   SUBROUTINE common_Z_lP_grid_available(bs_env, ri_rs_grid_points, &
                                         common_grid_available, have_fitted_columns)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points
      LOGICAL, INTENT(OUT)                               :: common_grid_available, &
                                                            have_fitted_columns

      CHARACTER(LEN=*), PARAMETER :: routineN = 'common_Z_lP_grid_available'

      INTEGER                                            :: AB_block, atom_A, atom_B, fit_atom, &
                                                            handle, iatom, igrid, n_to_a, natom
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: active_atom
      REAL(KIND=dp)                                      :: cutoff_ri
      REAL(KIND=dp), DIMENSION(3)                        :: center
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      particle_set => bs_env%ri_rs%particle_set
      natom = bs_env%n_atom
      common_grid_available = .TRUE.
      have_fitted_columns = .FALSE.
      ALLOCATE (active_atom(natom))
      DO fit_atom = 1, natom
         IF (bs_env%auto_ri%sizes_opt_RI(fit_atom) == 0) CYCLE
         have_fitted_columns = .TRUE.
         active_atom = .FALSE.
         DO AB_block = 1, bs_env%auto_ri%AB_block_count
            atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
            atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
            n_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
            IF (.NOT. (fit_atom == atom_A .AND. n_to_a > 0) .AND. &
                .NOT. (fit_atom == atom_B .AND. &
                       n_to_a < bs_env%auto_ri%AB_size_opt_RI(AB_block))) CYCLE
            active_atom(atom_A) = .TRUE.
            IF (atom_B /= atom_A) active_atom(atom_B) = .TRUE.
         END DO

         center = particle_set(fit_atom)%r
         cutoff_ri = 0.0_dp
         DO iatom = 1, natom
            IF (.NOT. active_atom(iatom)) CYCLE
            cutoff_ri = MAX(cutoff_ri, bs_env%ri_metric%cutoff_radius + &
                            bs_env%ri_rs%radius_ri_per_atom(iatom) + &
                            NORM2(center - particle_set(iatom)%r))
         END DO

         DO igrid = 1, bs_env%ri_rs%n_grid_points
            IF (NORM2(ri_rs_grid_points(1:3, igrid) - center) > cutoff_ri) THEN
               common_grid_available = .FALSE.
               EXIT
            END IF
         END DO
         IF (.NOT. common_grid_available) EXIT

         DO iatom = 1, natom
            IF (NORM2(particle_set(iatom)%r - center) > &
                bs_env%ri_rs%radius_ao_per_atom(iatom) + cutoff_ri) THEN
               common_grid_available = .FALSE.
               EXIT
            END IF
         END DO
         IF (.NOT. common_grid_available) EXIT
      END DO
      DEALLOCATE (active_atom)

      CALL timestop(handle)

   END SUBROUTINE common_Z_lP_grid_available

! **************************************************************************************************
!> \brief Evaluates the AO collocation matrix on the complete RI-RS grid,
!>
!>          ϕ_lμ = ϕ_μ(r_l),  l = 1, ..., N_grid.
!>
!>        "Complete common grid" means that the same full set of grid points is valid for every
!>        optimized RI column block. The enclosing sphere is only an implementation device passed
!>        to build_phi_on_sphere; it contains every r_l and every AO center, so its chosen center
!>        does not select an AA or AB contraction.
!> \param bs_env ...
!> \param qs_kind_set ...
!> \param ri_rs_grid_points ...
!> \param local_grid_idx ...
!> \param ngrid ...
!> \param phi_local ...
!> \param ao_atoms global indices of the active AO atoms
!> \param ao_col_start local phi column starts, including the terminal column offset
!> \param ao_nonzero exact nonzero AO support on the retained grid rows
! **************************************************************************************************
   SUBROUTINE build_phi_on_complete_grid(bs_env, qs_kind_set, &
                                         ri_rs_grid_points, local_grid_idx, ngrid, &
                                         phi_local, ao_atoms, ao_col_start, ao_nonzero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: local_grid_idx
      INTEGER, INTENT(OUT)                               :: ngrid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: phi_local
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ao_atoms, ao_col_start
      LOGICAL, ALLOCATABLE, INTENT(OUT)                  :: ao_nonzero(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_phi_on_complete_grid'

      INTEGER                                            :: handle, iatom, igrid, reference_atom
      REAL(KIND=dp)                                      :: cutoff_ri
      REAL(KIND=dp), DIMENSION(3)                        :: center
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      particle_set => bs_env%ri_rs%particle_set
      reference_atom = 1
      center = particle_set(reference_atom)%r
      cutoff_ri = 0.0_dp
      DO igrid = 1, bs_env%ri_rs%n_grid_points
         cutoff_ri = MAX(cutoff_ri, NORM2(ri_rs_grid_points(1:3, igrid) - center))
      END DO
      DO iatom = 1, bs_env%n_atom
         cutoff_ri = MAX(cutoff_ri, NORM2(particle_set(iatom)%r - center))
      END DO
      cutoff_ri = cutoff_ri + 1.0_dp

      CALL build_phi_on_sphere(bs_env, qs_kind_set, &
                               ri_rs_grid_points, reference_atom, cutoff_ri, &
                               local_grid_idx, ngrid, phi_local, ao_atoms, ao_col_start, &
                               center=center, ao_nonzero=ao_nonzero)

      CALL timestop(handle)

   END SUBROUTINE build_phi_on_complete_grid

! **************************************************************************************************
!> \brief Writes Z_l,p0+p += z_lp into the process-owned grid-row blocks of one atomic column
!>        block. The routine performs no explicit MPI communication; DBCSR summation combines
!>        contributions when a block receives columns from more than one contraction group.
!> \param mat_Z_lP ...
!> \param z_block ...
!> \param local_grid_idx ...
!> \param n_local_grid ...
!> \param atom_index ...
!> \param first_column ...
!> \param atom_block_size ...
!> \param r_blk_sizes ...
!> \param row_offset ...
!> \param row_dist ...
!> \param myprow ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE add_Z_lP_columns(mat_Z_lP, z_block, local_grid_idx, n_local_grid, atom_index, &
                               first_column, atom_block_size, r_blk_sizes, row_offset, row_dist, &
                               myprow, eps_filter)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_Z_lP
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: z_block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_grid_idx
      INTEGER, INTENT(IN)                                :: n_local_grid, atom_index, first_column, &
                                                            atom_block_size
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes, row_offset, row_dist
      INTEGER, INTENT(IN)                                :: myprow
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_Z_lP_columns'

      INTEGER                                            :: current_chunk_size, g_pt, handle, i_blk, &
                                                            loc_ptr, ncolumn, r_end, r_start
      LOGICAL                                            :: row_owned
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Z_blk

      CALL timeset(routineN, handle)

      ncolumn = SIZE(z_block, 2)
      CPASSERT(first_column > 0)
      CPASSERT(first_column + ncolumn - 1 <= atom_block_size)
      ALLOCATE (Z_blk(MAXVAL(r_blk_sizes), atom_block_size), source=0.0_dp)
      loc_ptr = 1
      DO i_blk = 1, SIZE(r_blk_sizes)
         r_start = row_offset(i_blk) + 1
         r_end = row_offset(i_blk) + r_blk_sizes(i_blk)
         current_chunk_size = r_blk_sizes(i_blk)
         row_owned = row_dist(i_blk) == myprow
         Z_blk = 0.0_dp
         DO WHILE (loc_ptr <= n_local_grid)
            g_pt = local_grid_idx(loc_ptr)
            IF (g_pt > r_end) EXIT
            IF (row_owned) THEN
               Z_blk(g_pt - r_start + 1, first_column:first_column + ncolumn - 1) = &
                  z_block(loc_ptr, 1:ncolumn)
            END IF
            loc_ptr = loc_ptr + 1
         END DO
         IF (row_owned .AND. MAXVAL(ABS(Z_blk(1:current_chunk_size, :))) > eps_filter) THEN
            CALL dbcsr_put_block(mat_Z_lP, row=i_blk, col=atom_index, &
                                 block=Z_blk(1:current_chunk_size, :), summation=.TRUE.)
         END IF
      END DO
      DEALLOCATE (Z_blk)

      CALL timestop(handle)

   END SUBROUTINE add_Z_lP_columns

! **************************************************************************************************
!> \brief Computes d_lp = Σ_{μνP} ϕ_μ(r_l)ϕ_ν(r_l)(μν|P)U_Pp for all AA and AB
!>        columns p of one RI-RS matrix block. P can belong to either atom of an AB contraction.
!> \param bs_env ...
!> \param ctx ...
!> \param phi_val ...
!> \param ao_atoms global indices of the active AO atoms
!> \param ao_col_start local phi column starts, including the terminal column offset
!> \param ao_nonzero exact nonzero AO support on the retained grid rows
!> \param d_lp ...
!> \param ri_atoms global indices of contributing RI atoms in ascending order
!> \param ri_coefficients ...
!> \param max_ao_size ...
!> \param atom_j_mepos ...
!> \param atom_j_stride ...
! **************************************************************************************************
   SUBROUTINE compute_d_lp_auto_ri(bs_env, ctx, phi_val, ao_atoms, ao_col_start, ao_nonzero, d_lp, &
                                   ri_atoms, ri_coefficients, max_ao_size, atom_j_mepos, &
                                   atom_j_stride)
!$    USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(IN)                                      :: phi_val
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_atoms, ao_col_start
      LOGICAL, DIMENSION(:, :), INTENT(IN)               :: ao_nonzero
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), INTENT(IN) :: ri_coefficients
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ri_atoms
      INTEGER, INTENT(IN)                                :: max_ao_size, atom_j_mepos, atom_j_stride

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_d_lp_auto_ri'
      INTEGER, PARAMETER                                 :: grid_chunk = 1024

      INTEGER :: jj, kk, active_column, iatom, jatom, katom, c, handle, i_thread, j, jk_idx, jsize, &
                 jstart, k, ksize, kstart, l, l0, &
                 max_active, nactive, ncol, n_grid, nRI_ref, nthreads, ri, slot, thread_id
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_ncol
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: atom_column
      LOGICAL                                            :: any_integral, screened, single_atom
      REAL(KIND=dp)                                      :: pair_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: int_2d_prv, rho_chunk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: coefficient_compact, d_lp_threads, &
                                                            int_3c_prv, int_3c_atom
      TYPE(gw_3c_ws_type)                                :: ws

      INTEGER, ALLOCATABLE, DIMENSION(:) :: grid_index
      INTEGER :: n_grid_pair, grid_l, point
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: grid_result

      CALL timeset(routineN, handle)
      n_grid = SIZE(d_lp, 1)
      ncol = SIZE(d_lp, 2)
      single_atom = SIZE(ri_atoms) == 1
      nthreads = 1
!$    nthreads = omp_get_max_threads()
      CPASSERT(SIZE(ri_coefficients, 2) == ncol)
      CPASSERT(SIZE(ri_coefficients, 3) == SIZE(ri_atoms))
      CPASSERT(SIZE(ri_atoms) > 0)
      max_active = 0
      IF (.NOT. single_atom) THEN
         ALLOCATE (atom_ncol(SIZE(ri_atoms)), atom_column(ncol, SIZE(ri_atoms)))
         atom_ncol = 0
         atom_column = 0
         DO slot = 1, SIZE(ri_atoms)
            nRI_ref = get_ref_RI_size(bs_env, ri_atoms(slot))
            DO ri = 1, ncol
               IF (.NOT. ANY(ri_coefficients(1:nRI_ref, ri, slot) /= 0.0_dp)) CYCLE
               atom_ncol(slot) = atom_ncol(slot) + 1
               atom_column(atom_ncol(slot), slot) = ri
            END DO
         END DO
         max_active = MAXVAL(atom_ncol)
         CPASSERT(max_active > 0)
         ALLOCATE (coefficient_compact(SIZE(ri_coefficients, 1), max_active, &
                                       SIZE(ri_atoms)), source=0.0_dp)
         DO slot = 1, SIZE(ri_atoms)
            nRI_ref = get_ref_RI_size(bs_env, ri_atoms(slot))
            DO active_column = 1, atom_ncol(slot)
               ri = atom_column(active_column, slot)
               coefficient_compact(1:nRI_ref, active_column, slot) = &
                  ri_coefficients(1:nRI_ref, ri, slot)
            END DO
         END DO
      END IF
      ALLOCATE (d_lp_threads(n_grid, ncol, nthreads), source=0.0_dp)

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, ctx, phi_val, ao_atoms, ao_col_start, ao_nonzero, d_lp, n_grid, ri_coefficients, &
      !$OMP        ri_atoms, single_atom, max_ao_size, atom_j_mepos, atom_j_stride, ncol, &
      !$OMP        d_lp_threads, nthreads, atom_ncol, &
      !$OMP        atom_column, coefficient_compact, max_active) &
      !$OMP PRIVATE(grid_index, n_grid_pair, grid_l, point, grid_result, &
      !$OMP         jj, kk, active_column, iatom, jatom, katom, c, i_thread, j, jk_idx, jsize, jstart, &
      !$OMP         k, ksize, kstart, l, l0, &
      !$OMP         nactive, nRI_ref, ri, slot, any_integral, screened, pair_factor, &
      !$OMP         int_2d_prv, rho_chunk, &
      !$OMP         int_3c_prv, int_3c_atom, ws, thread_id)

      thread_id = 1
!$    thread_id = omp_get_thread_num() + 1

      CALL gw_3c_ws_create(ws, ctx)
      ALLOCATE (int_3c_prv(max_ao_size, max_ao_size, ncol))
      IF (.NOT. single_atom) ALLOCATE (int_3c_atom(max_ao_size, max_ao_size, max_active))
      ALLOCATE (int_2d_prv(max_ao_size*max_ao_size, ncol))
      ALLOCATE (rho_chunk(grid_chunk, max_ao_size*max_ao_size))
      ALLOCATE (grid_index(n_grid), grid_result(grid_chunk, ncol))

      !$OMP DO SCHEDULE(DYNAMIC)
      DO jj = 1, SIZE(ao_atoms)
         jatom = ao_atoms(jj)
         IF (MOD(jatom - 1, atom_j_stride) /= atom_j_mepos) CYCLE
         DO kk = jj, SIZE(ao_atoms)
            katom = ao_atoms(kk)
            jstart = ao_col_start(jj)
            kstart = ao_col_start(kk)
            jsize = ao_col_start(jj + 1) - jstart
            ksize = ao_col_start(kk + 1) - kstart

            n_grid_pair = 0
            DO grid_l = 1, n_grid
               IF (.NOT. (ao_nonzero(grid_l, jj) .AND. ao_nonzero(grid_l, kk))) CYCLE
               n_grid_pair = n_grid_pair + 1
               grid_index(n_grid_pair) = grid_l
            END DO
            IF (n_grid_pair == 0) CYCLE
            int_3c_prv(1:jsize, 1:ksize, 1:ncol) = 0.0_dp
            IF (single_atom) THEN
               ! A single contributor writes the transformed block without packing or scattering.
               CALL build_3c_integral_block_auto_ri_ctx( &
                  int_3c_prv(1:jsize, 1:ksize, 1:ncol), ctx, ws, &
                  atom_j=jatom, atom_k=katom, atom_i=ri_atoms(1), &
                  transform=ri_coefficients(:, :, 1), transform_row=1, screened=screened)
               IF (screened) CYCLE
            ELSE
               any_integral = .FALSE.
               DO slot = 1, SIZE(ri_atoms)
                  iatom = ri_atoms(slot)
                  nactive = atom_ncol(slot)
                  IF (nactive == 0) CYCLE
                  nRI_ref = get_ref_RI_size(bs_env, iatom)
                  int_3c_atom(1:jsize, 1:ksize, 1:nactive) = 0.0_dp
                  CALL build_3c_integral_block_auto_ri_ctx( &
                     int_3c_atom(1:jsize, 1:ksize, 1:nactive), ctx, ws, &
                     atom_j=jatom, atom_k=katom, atom_i=iatom, &
                     transform=coefficient_compact(1:nRI_ref, 1:nactive, slot), &
                     transform_row=1, screened=screened)
                  IF (.NOT. screened) THEN
                     any_integral = .TRUE.
                     DO active_column = 1, nactive
                        ri = atom_column(active_column, slot)
                        int_3c_prv(1:jsize, 1:ksize, ri) = &
                           int_3c_prv(1:jsize, 1:ksize, ri) + &
                           int_3c_atom(1:jsize, 1:ksize, active_column)
                     END DO
                  END IF
               END DO
               IF (.NOT. any_integral) CYCLE
            END IF

            DO ri = 1, ncol
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     int_2d_prv(jk_idx, ri) = int_3c_prv(j, k, ri)
                  END DO
               END DO
            END DO

            pair_factor = 1.0_dp
            IF (jatom /= katom) pair_factor = 2.0_dp
            DO l0 = 1, n_grid_pair, grid_chunk
               c = MIN(grid_chunk, n_grid_pair - l0 + 1)
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     DO l = 1, c
                        point = grid_index(l0 + l - 1)
                        rho_chunk(l, jk_idx) = phi_val(point, jstart + j - 1)* &
                                               phi_val(point, kstart + k - 1)
                     END DO
                  END DO
               END DO
               CALL dgemm('N', 'N', c, ncol, jsize*ksize, pair_factor, rho_chunk, grid_chunk, &
                          int_2d_prv, max_ao_size*max_ao_size, 0.0_dp, &
                          grid_result, grid_chunk)
               DO ri = 1, ncol
                  DO l = 1, c
                     point = grid_index(l0 + l - 1)
                     d_lp_threads(point, ri, thread_id) = &
                        d_lp_threads(point, ri, thread_id) + grid_result(l, ri)
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END DO

      !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
      DO ri = 1, ncol
         DO l = 1, n_grid
            DO i_thread = 1, nthreads
               d_lp(l, ri) = d_lp(l, ri) + d_lp_threads(l, ri, i_thread)
            END DO
         END DO
      END DO
      !$OMP END DO
      IF (.NOT. single_atom) DEALLOCATE (int_3c_atom)
      DEALLOCATE (int_3c_prv, int_2d_prv, rho_chunk)
      DEALLOCATE (grid_index, grid_result)
      CALL gw_3c_ws_release(ws)
      !$OMP END PARALLEL

      IF (.NOT. single_atom) DEALLOCATE (coefficient_compact, atom_column, atom_ncol)
      DEALLOCATE (d_lp_threads)

      CALL timestop(handle)

   END SUBROUTINE compute_d_lp_auto_ri

! **************************************************************************************************
!> \brief Computes d_lp = Σ_A Σ_{P∈A} d_lP U_Pp in rank-sized batches on a common grid.
!> \param bs_env ...
!> \param ctx ...
!> \param phi_val ...
!> \param ao_atoms global indices of the active AO atoms
!> \param ao_col_start local phi column starts, including the terminal column offset
!> \param ao_nonzero exact nonzero AO support on the retained grid rows
!> \param d_lp ...
!> \param n_grid ...
!> \param max_ao_size ...
!> \param mepos ...
!> \param num_pe ...
! **************************************************************************************************
   SUBROUTINE compute_d_lp_auto_ri_batch(bs_env, ctx, phi_val, ao_atoms, ao_col_start, ao_nonzero, d_lp, n_grid, &
                                         max_ao_size, mepos, num_pe)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(IN)                                      :: phi_val
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_atoms, ao_col_start
      LOGICAL, DIMENSION(:, :), INTENT(IN)               :: ao_nonzero
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      INTEGER, INTENT(IN)                                :: n_grid, max_ao_size, mepos, num_pe

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_d_lp_auto_ri_batch'

      INTEGER                                            :: column, handle, iatom, n_done, n_total, &
                                                            ncol_batch
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: global_map
      REAL(KIND=dp)                                      :: item_start_time
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_batch
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: U_Pp

      CALL timeset(routineN, handle)

      n_total = 0
      DO iatom = mepos + 1, bs_env%n_atom, num_pe
         n_total = n_total + 1
      END DO
      n_done = 0

      DO iatom = mepos + 1, bs_env%n_atom, num_pe
         item_start_time = m_walltime()
         CALL collect_auto_ri_columns_for_atom(bs_env, iatom, U_Pp, global_map)
         ncol_batch = SIZE(global_map)
         IF (ncol_batch == 0) THEN
            DEALLOCATE (U_Pp, global_map)
            n_done = n_done + 1
            CALL print_Z_lP_progress(bs_env, n_done, n_total, &
                                     m_walltime() - item_start_time)
            CYCLE
         END IF

         ALLOCATE (d_lp_batch(n_grid, ncol_batch), source=0.0_dp)
         CALL compute_d_lp_auto_ri(bs_env, ctx, phi_val, ao_atoms, ao_col_start, ao_nonzero, d_lp_batch, &
                                   [iatom], U_Pp, max_ao_size, 0, 1)
         DO column = 1, ncol_batch
            d_lp(:, global_map(column)) = d_lp(:, global_map(column)) + d_lp_batch(:, column)
         END DO
         DEALLOCATE (U_Pp, global_map, d_lp_batch)
         n_done = n_done + 1
         CALL print_Z_lP_progress(bs_env, n_done, n_total, &
                                  m_walltime() - item_start_time)
      END DO
      CALL timestop(handle)

   END SUBROUTINE compute_d_lp_auto_ri_batch

! **************************************************************************************************
!> \brief Collects every optimized column p containing reference functions P on atom A. It returns
!>
!>          U_Pp^A
!>
!>        together with the global optimized-column index of q.
!> \param bs_env ...
!> \param iatom ...
!> \param U_Pp ...
!> \param global_map ...
! **************************************************************************************************
   SUBROUTINE collect_auto_ri_columns_for_atom(bs_env, iatom, U_Pp, global_map)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: U_Pp(:, :, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: global_map(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_auto_ri_columns_for_atom'

      INTEGER :: AB_block, atom_A, atom_B, column, global_column, handle, local_column, ncol, &
         ncol_batch, nRI_ref, nRI_ref_A, output_offset, row_first
      REAL(KIND=dp), ALLOCATABLE                         :: U_Pp_AB(:, :)

      CALL timeset(routineN, handle)

      ncol_batch = 0
      DO AB_block = 1, bs_env%auto_ri%AB_block_count
         atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
         atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
         IF (iatom == atom_A .OR. iatom == atom_B) THEN
            ncol_batch = ncol_batch + bs_env%auto_ri%AB_size_opt_RI(AB_block)
         END IF
      END DO

      nRI_ref = get_ref_RI_size(bs_env, iatom)
      ALLOCATE (U_Pp(nRI_ref, ncol_batch, 1), source=0.0_dp)
      ALLOCATE (global_map(ncol_batch))
      column = 0
      DO AB_block = 1, bs_env%auto_ri%AB_block_count
         atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
         atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
         IF (iatom /= atom_A .AND. iatom /= atom_B) CYCLE
         ncol = bs_env%auto_ri%AB_size_opt_RI(AB_block)
         CALL get_U_Pp_AB(bs_env%auto_ri, AB_block, U_Pp_AB)
         row_first = 1
         IF (iatom == atom_B .AND. atom_B /= atom_A) THEN
            nRI_ref_A = get_ref_RI_size(bs_env, atom_A)
            row_first = 1 + nRI_ref_A
         END IF
         U_Pp(:, column + 1:column + ncol, 1) = &
            U_Pp_AB( &
            row_first:row_first + nRI_ref - 1, 1:ncol)

         DO local_column = 1, ncol
            IF (local_column <= bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)) THEN
               output_offset = SUM(bs_env%auto_ri%sizes_opt_RI(:atom_A - 1)) + &
                               bs_env%auto_ri%AB_first_p_A(AB_block) - 1
               global_column = output_offset + local_column
            ELSE
               CPASSERT(atom_B /= atom_A)
               output_offset = SUM(bs_env%auto_ri%sizes_opt_RI(:atom_B - 1)) + &
                               bs_env%auto_ri%AB_first_p_B(AB_block) - 1
               global_column = output_offset + local_column - &
                               bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
            END IF
            global_map(column + local_column) = global_column
         END DO
         column = column + ncol
         DEALLOCATE (U_Pp_AB)
      END DO
      CPASSERT(column == ncol_batch)

      CALL timestop(handle)

   END SUBROUTINE collect_auto_ri_columns_for_atom

! **************************************************************************************************
!> \brief Unpacks one AB contraction matrix U_Pp; B=A denotes an AA block.
!> \param auto_ri ...
!> \param AB_block ...
!> \param U_Pp ...
! **************************************************************************************************
   SUBROUTINE get_U_Pp_AB(auto_ri, AB_block, U_Pp)
      TYPE(auto_ri_type), INTENT(IN)                     :: auto_ri
      INTEGER, INTENT(IN)                                :: AB_block
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: U_Pp(:, :)

      INTEGER                                            :: first, last, ncolumn, nrow

      nrow = auto_ri%AB_size_ref_RI(AB_block)
      ncolumn = auto_ri%AB_size_opt_RI(AB_block)
      first = auto_ri%U_Pp_AB_offset(AB_block)
      last = first + nrow*ncolumn - 1
      ALLOCATE (U_Pp(nrow, ncolumn))
      U_Pp(:, :) = RESHAPE(auto_ri%U_Pp_AB(first:last), [nrow, ncolumn])

   END SUBROUTINE get_U_Pp_AB

! **************************************************************************************************
!> \brief Computes the fitting radius for every optimized atomic column block. If q assigned to A
!>        contains reference functions on B, the required radius is
!>
!>          R_A^fit = max_B [r_c + R_B^RI + |R_A - R_B|].
!> \param bs_env ...
!> \param radius ...
! **************************************************************************************************
   SUBROUTINE compute_auto_ri_grid_radii(bs_env, radius)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: radius(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_auto_ri_grid_radii'

      INTEGER                                            :: a, AB_block, b, handle, n_to_a, ncol
      REAL(KIND=dp)                                      :: cutoff, distance
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      particle_set => bs_env%ri_rs%particle_set
      radius = 0.0_dp
      cutoff = bs_env%ri_metric%cutoff_radius
      DO AB_block = 1, bs_env%auto_ri%AB_block_count
         a = bs_env%auto_ri%AB_atom_A(AB_block)
         b = bs_env%auto_ri%AB_atom_B(AB_block)
         n_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
         ncol = bs_env%auto_ri%AB_size_opt_RI(AB_block)
         IF (n_to_a > 0) THEN
            radius(a) = MAX(radius(a), cutoff + bs_env%ri_rs%radius_ri_per_atom(a))
            IF (b /= a) THEN
               distance = NORM2(particle_set(a)%r - particle_set(b)%r)
               radius(a) = MAX(radius(a), cutoff + bs_env%ri_rs%radius_ri_per_atom(b) + distance)
            END IF
         END IF
         IF (n_to_a < ncol) THEN
            CPASSERT(b /= a)
            distance = NORM2(particle_set(b)%r - particle_set(a)%r)
            radius(b) = MAX(radius(b), cutoff + bs_env%ri_rs%radius_ri_per_atom(b), &
                            cutoff + bs_env%ri_rs%radius_ri_per_atom(a) + distance)
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE compute_auto_ri_grid_radii

! **************************************************************************************************
!> \brief Computes
!>
!>          d_lp = Σ_A Σ_{P∈A} d_lP U_Pp^A,
!>          d_lP = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P),
!>
!>        on the union of the atomic RI-RS grids needed by q. Each atom's three-center integrals are
!>        evaluated once and transformed with U_Pp^A.
!> \param qs_env ...
!> \param bs_env ...
!> \param ctx ...
!> \param grid ...
!> \param mat_phi ...
!> \param mat_rhs ...
!> \param max_ao_size ...
! **************************************************************************************************
   SUBROUTINE compute_auto_ri_d_lp(qs_env, bs_env, ctx, grid, mat_phi, mat_rhs, max_ao_size)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_phi
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_rhs
      INTEGER, INTENT(IN)                                :: max_ao_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_auto_ri_d_lp'

      INTEGER :: AB_block, atom_A, atom_B, first, fit_atom, handle, handle_project, handle_rhs, l, &
         last, n_done, n_first_p_ABs, n_to_a, n_total, n_union, natom, ncol, ngrid, npcol, nprow, &
         ri_atom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_atoms, ao_col_start, global_map, &
                                                            local_index, row_offset, union_index
      INTEGER, DIMENSION(:), POINTER                     :: AB_row_dist, col_dist, &
                                                            first_p_ABs_per_atom, retained_size, &
                                                            row_dist, row_size
      LOGICAL, ALLOCATABLE                               :: ao_nonzero(:, :)
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: needed_fit_atom, union_mask
      REAL(KIND=dp)                                      :: item_start_time, radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fit_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: column_map, phi, rhs, union_grid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: U_Pp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dist_AB, dist_phi, dist_t
      TYPE(dbcsr_type)                                   :: AB_d_lp, transform
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, cell=cell)
      natom = bs_env%n_atom
      CALL dbcsr_get_info(mat_phi, row_blk_size=row_size, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist, nprows=nprow, npcols=npcol)
      ALLOCATE (first_p_ABs_per_atom(natom), col_dist(natom), AB_row_dist(natom), &
                retained_size(natom), row_offset(SIZE(row_size)))
      retained_size(:) = bs_env%auto_ri%sizes_opt_RI
      DO ri_atom = 1, natom
         first_p_ABs_per_atom(ri_atom) = 0
         DO AB_block = 1, bs_env%auto_ri%AB_block_count
            IF (ri_atom == bs_env%auto_ri%AB_atom_A(AB_block) .OR. &
                ri_atom == bs_env%auto_ri%AB_atom_B(AB_block)) THEN
               first_p_ABs_per_atom(ri_atom) = &
                  first_p_ABs_per_atom(ri_atom) + &
                  bs_env%auto_ri%AB_size_opt_RI(AB_block)
            END IF
         END DO
         col_dist(ri_atom) = MOD(ri_atom - 1, npcol)
         AB_row_dist(ri_atom) = MOD(ri_atom - 1, nprow)
      END DO
      row_offset(1) = 0
      DO l = 2, SIZE(row_size)
         row_offset(l) = row_offset(l - 1) + row_size(l - 1)
      END DO
      CALL dbcsr_distribution_new(dist_AB, template=dist_phi, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(AB_d_lp, name='AUTO_RI AA/AB d_lp', dist=dist_AB, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_size, &
                        col_blk_size=first_p_ABs_per_atom)
      CALL dbcsr_distribution_new(dist_t, template=dist_phi, &
                                  row_dist=AB_row_dist, col_dist=col_dist)
      CALL dbcsr_create(transform, name='AUTO_RI U_Pp', dist=dist_t, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=first_p_ABs_per_atom, &
                        col_blk_size=retained_size)
      CALL dbcsr_create(mat_rhs, name='AUTO_RI optimized d_lp', dist=dist_AB, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_size, &
                        col_blk_size=retained_size)
      ALLOCATE (fit_radius(natom))
      CALL compute_auto_ri_grid_radii(bs_env, fit_radius)
      ALLOCATE (needed_fit_atom(natom), union_mask(SIZE(grid, 2)))
      n_total = 0
      DO ri_atom = para_env%mepos + 1, natom, para_env%num_pe
         n_total = n_total + 1
      END DO
      n_done = 0
      CALL timeset(routineN//'_AB_d_lp', handle_rhs)
      DO ri_atom = para_env%mepos + 1, natom, para_env%num_pe
         item_start_time = m_walltime()
         needed_fit_atom = .FALSE.
         DO AB_block = 1, bs_env%auto_ri%AB_block_count
            atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
            atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
            IF (ri_atom /= atom_A .AND. ri_atom /= atom_B) CYCLE
            n_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
            IF (n_to_a > 0) needed_fit_atom(atom_A) = .TRUE.
            IF (n_to_a < bs_env%auto_ri%AB_size_opt_RI(AB_block)) THEN
               needed_fit_atom(atom_B) = .TRUE.
            END IF
         END DO
         IF (.NOT. ANY(needed_fit_atom)) THEN
            n_done = n_done + 1
            CALL print_Z_lP_progress(bs_env, n_done, n_total, &
                                     m_walltime() - item_start_time)
            CYCLE
         END IF
         union_mask = .FALSE.
         radius = 0.0_dp
         DO fit_atom = 1, natom
            IF (.NOT. needed_fit_atom(fit_atom)) CYCLE
            radius = MAX(radius, fit_radius(fit_atom) + &
                         NORM2(particle_set(fit_atom)%r - particle_set(ri_atom)%r))
            DO l = 1, SIZE(grid, 2)
               IF (NORM2(grid(:, l) - particle_set(fit_atom)%r) <= &
                   fit_radius(fit_atom)) union_mask(l) = .TRUE.
            END DO
         END DO
         n_union = COUNT(union_mask)
         ALLOCATE (union_index(n_union), union_grid(3, n_union))
         n_union = 0
         DO l = 1, SIZE(grid, 2)
            IF (.NOT. union_mask(l)) CYCLE
            n_union = n_union + 1
            union_index(n_union) = l
            union_grid(:, n_union) = grid(:, l)
         END DO
         CALL build_phi_on_sphere(bs_env, qs_kind_set, union_grid, ri_atom, &
                                  radius + 1.0_dp, local_index, ngrid, phi, ao_atoms, ao_col_start, &
                                  ao_nonzero=ao_nonzero)
         local_index(1:ngrid) = union_index(local_index(1:ngrid))
         n_first_p_ABs = first_p_ABs_per_atom(ri_atom)
         CALL collect_auto_ri_columns_for_atom(bs_env, ri_atom, U_Pp, global_map)
         CPASSERT(SIZE(global_map) == n_first_p_ABs)
         first = 1
         DO fit_atom = 1, natom
            ncol = retained_size(fit_atom)
            last = first + ncol - 1
            IF (ANY(global_map >= first .AND. global_map <= last)) THEN
               ALLOCATE (column_map(n_first_p_ABs, ncol), source=0.0_dp)
               DO l = 1, n_first_p_ABs
                  IF (global_map(l) < first .OR. global_map(l) > last) CYCLE
                  column_map(l, global_map(l) - first + 1) = 1.0_dp
               END DO
               CALL dbcsr_put_block(transform, ri_atom, fit_atom, column_map)
               DEALLOCATE (column_map)
            END IF
            first = last + 1
         END DO
         ALLOCATE (rhs(ngrid, n_first_p_ABs), source=0.0_dp)
         CALL compute_d_lp_auto_ri(bs_env, ctx, phi, ao_atoms, ao_col_start, ao_nonzero, rhs, &
                                   [ri_atom], U_Pp, max_ao_size, 0, 1)
         DEALLOCATE (ao_atoms, ao_col_start, ao_nonzero)
         CALL store_Z_lP_columns(AB_d_lp, rhs, local_index, ngrid, n_first_p_ABs, ri_atom, &
                                 row_size, row_offset, 0.0_dp)
         DEALLOCATE (rhs, U_Pp, global_map, phi, local_index, &
                     union_grid, union_index)
         n_done = n_done + 1
         CALL print_Z_lP_progress(bs_env, n_done, n_total, &
                                  m_walltime() - item_start_time)
      END DO
      CALL dbcsr_finalize(transform)
      CALL dbcsr_finalize(AB_d_lp)
      CALL timestop(handle_rhs)
      CALL timeset(routineN//'_transform', handle_project)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, AB_d_lp, transform, 0.0_dp, &
                          mat_rhs, filter_eps=0.0_dp)
      CALL timestop(handle_project)
      CALL dbcsr_release(AB_d_lp)
      CALL dbcsr_release(transform)
      CALL dbcsr_distribution_release(dist_AB)
      CALL dbcsr_distribution_release(dist_t)
      DEALLOCATE (first_p_ABs_per_atom, retained_size, col_dist, AB_row_dist, row_offset, &
                  fit_radius, needed_fit_atom, union_mask)
      CALL timestop(handle)
   END SUBROUTINE compute_auto_ri_d_lp

! **************************************************************************************************
!> \brief Extracts rhs(i,q) = d_{local_index(i),q} from one optimized atomic DBCSR column block.
!> \param mat_rhs ...
!> \param atom_index ...
!> \param local_index ...
!> \param row_offset ...
!> \param rhs ...
! **************************************************************************************************
   SUBROUTINE extract_atom_d_lp(mat_rhs, atom_index, local_index, row_offset, rhs)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_rhs
      INTEGER, INTENT(IN)                                :: atom_index
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_index, row_offset
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: rhs

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'extract_atom_d_lp'

      INTEGER                                            :: handle, l, next_row, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block

      CALL timeset(routineN, handle)

      rhs = 0.0_dp
      row = 0
      next_row = 1
      NULLIFY (block)
      DO l = 1, SIZE(rhs, 1)
         DO WHILE (next_row <= SIZE(row_offset))
            IF (row_offset(next_row) >= local_index(l)) EXIT
            row = next_row
            next_row = next_row + 1
            CALL dbcsr_get_block_p(mat_rhs, row, atom_index, block, found)
            IF (.NOT. found) NULLIFY (block)
         END DO
         IF (ASSOCIATED(block)) rhs(l, :) = block(local_index(l) - row_offset(row), :)
      END DO

      CALL timestop(handle)
   END SUBROUTINE extract_atom_d_lp

! **************************************************************************************************
!> \brief Solves Σ_l' D_ll' Z_l'q = d_lp for each optimized atomic column block. Small grids are
!>        gathered within one process column and solved by a local Cholesky factorization.
!> \param qs_env ...
!> \param bs_env ...
!> \param grid ...
!> \param mat_rhs ...
!> \param mat_z ...
! **************************************************************************************************
   SUBROUTINE fit_auto_ri_z_lp(qs_env, bs_env, grid, mat_rhs, mat_z)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: grid(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_rhs, mat_z

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fit_auto_ri_z_lp'

      INTEGER                                            :: atom_index, base, group_size, handle, &
                                                            handle_gather, i, mypcol, n_big, &
                                                            n_small, natom, ncol, ngrid, npcol, &
                                                            slot
      INTEGER, ALLOCATABLE :: all_index(:), ao_atoms(:), ao_col_start(:), big_list(:), &
         grid_size(:), local_index(:), row_offset(:), small_list(:)
      INTEGER, POINTER                                   :: row_size(:)
      LOGICAL, ALLOCATABLE                               :: single_rank(:)
      REAL(KIND=dp), ALLOCATABLE                         :: atom_rhs(:, :), buffer(:, :), &
                                                            local_rhs(:, :), phi(:, :), radius(:)
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(mp_para_env_type), POINTER                    :: column_env, para_env
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, para_env=para_env, qs_kind_set=qs_kinds)
      CALL dbcsr_get_info(mat_rhs, distribution=distribution, row_blk_size=row_size)
      CALL dbcsr_distribution_get(distribution, mypcol=mypcol, npcols=npcol)
      ALLOCATE (column_env)
      CALL column_env%from_split(para_env, mypcol)
      natom = bs_env%n_atom
      ALLOCATE (radius(natom), all_index(SIZE(grid, 2)), row_offset(SIZE(row_size)))
      CALL compute_auto_ri_grid_radii(bs_env, radius)
      CALL classify_z_lp_atoms(bs_env, grid, radius, &
                               bs_env%auto_ri%sizes_opt_RI, grid_size, &
                               small_list, n_small, big_list, n_big, group_size)
      ALLOCATE (single_rank(natom), source=.FALSE.)
      single_rank(small_list(:n_small)) = .TRUE.
      DO i = 1, SIZE(all_index)
         all_index(i) = i
      END DO
      row_offset(1) = 0
      DO i = 2, SIZE(row_size)
         row_offset(i) = row_offset(i - 1) + row_size(i - 1)
      END DO
      DO base = mypcol + 1, natom, npcol*column_env%num_pe
         CALL timeset(routineN//'_gather_rhs', handle_gather)
         DO slot = 0, column_env%num_pe - 1
            atom_index = base + slot*npcol
            IF (atom_index > natom) EXIT
            IF (.NOT. single_rank(atom_index)) CYCLE
            ncol = bs_env%auto_ri%sizes_opt_RI(atom_index)
            IF (ncol == 0) CYCLE
            ALLOCATE (buffer(SIZE(grid, 2), ncol))
            CALL extract_atom_d_lp(mat_rhs, atom_index, all_index, row_offset, buffer)
            CALL column_env%sum(buffer, slot)
            IF (column_env%mepos == slot) THEN
               CALL MOVE_ALLOC(buffer, atom_rhs)
            ELSE
               DEALLOCATE (buffer)
            END IF
         END DO
         CALL timestop(handle_gather)
         atom_index = base + column_env%mepos*npcol
         IF (atom_index > natom) CYCLE
         IF (.NOT. single_rank(atom_index)) CYCLE
         ncol = bs_env%auto_ri%sizes_opt_RI(atom_index)
         IF (ncol == 0) CYCLE
         CALL build_phi_on_sphere(bs_env, qs_kinds, grid, atom_index, &
                                  radius(atom_index), local_index, ngrid, phi, ao_atoms, ao_col_start)
         DEALLOCATE (ao_atoms, ao_col_start)
         ALLOCATE (local_rhs(ngrid, ncol))
         DO i = 1, ngrid
            local_rhs(i, :) = atom_rhs(local_index(i), :)
         END DO
         DEALLOCATE (atom_rhs)
         CALL solve_ri_rs_fit(phi, local_rhs, bs_env%ri_rs%tikhonov)
         CALL store_Z_lP_columns(mat_z, local_rhs, local_index, ngrid, ncol, atom_index, &
                                 row_size, row_offset, 0.0_dp)
         DEALLOCATE (local_rhs, phi, local_index)
      END DO
      IF (n_big > 0) THEN
         CALL fit_distributed_auto_ri_z_lp(qs_env, bs_env, grid, mat_rhs, mat_z, radius, &
                                           big_list(:n_big), group_size, row_size, row_offset, &
                                           all_index)
      END IF
      DEALLOCATE (radius, all_index, row_offset, single_rank, grid_size, small_list, big_list)
      CALL column_env%free()
      DEALLOCATE (column_env)
      CALL dbcsr_finalize(mat_z)
      CALL timestop(handle)
   END SUBROUTINE fit_auto_ri_z_lp

! **************************************************************************************************
!> \brief Solves Σ_l' D_ll' Z_l'q = d_lp for large atomic grids with distributed Cholesky
!>        factorization. Each d_lp block is gathered once and distributed over its assigned rank
!>        group.
!> \param qs_env ...
!> \param bs_env ...
!> \param grid ...
!> \param mat_rhs ...
!> \param mat_z ...
!> \param radius ...
!> \param atom_list ...
!> \param group_size ...
!> \param row_size ...
!> \param row_offset ...
!> \param all_index ...
! **************************************************************************************************
   SUBROUTINE fit_distributed_auto_ri_z_lp(qs_env, bs_env, grid, mat_rhs, mat_z, radius, &
                                           atom_list, group_size, row_size, row_offset, &
                                           all_index)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: grid(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_rhs, mat_z
      REAL(KIND=dp), INTENT(IN)                          :: radius(:)
      INTEGER, INTENT(IN)                                :: atom_list(:), group_size, row_size(:), &
                                                            row_offset(:), all_index(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_distributed_auto_ri_z_lp'

      INTEGER                                            :: atom_index, base, handle, i, my_group, &
                                                            ncol, ngrid, ngroups, root, slot
      INTEGER, ALLOCATABLE                               :: ao_atoms(:), ao_col_start(:), &
                                                            local_index(:)
      REAL(KIND=dp), ALLOCATABLE                         :: atom_rhs(:, :), buffer(:, :), &
                                                            local_rhs(:, :), phi(:, :)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(mp_para_env_type), POINTER                    :: group_env, para_env
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, para_env=para_env, qs_kind_set=qs_kinds)
      ngroups = MAX(1, para_env%num_pe/group_size)
      my_group = MIN(para_env%mepos/group_size, ngroups - 1)
      ALLOCATE (group_env)
      CALL group_env%from_split(para_env, my_group)
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=group_env)
      DO base = 1, SIZE(atom_list), ngroups
         DO slot = 0, ngroups - 1
            IF (base + slot > SIZE(atom_list)) EXIT
            atom_index = atom_list(base + slot)
            ncol = bs_env%auto_ri%sizes_opt_RI(atom_index)
            IF (ncol == 0) CYCLE
            root = slot*group_size
            ALLOCATE (buffer(SIZE(grid, 2), ncol))
            CALL extract_atom_d_lp(mat_rhs, atom_index, all_index, row_offset, buffer)
            CALL para_env%sum(buffer, root)
            IF (para_env%mepos == root) THEN
               CALL MOVE_ALLOC(buffer, atom_rhs)
            ELSE
               DEALLOCATE (buffer)
            END IF
         END DO
         IF (base + my_group > SIZE(atom_list)) CYCLE
         atom_index = atom_list(base + my_group)
         ncol = bs_env%auto_ri%sizes_opt_RI(atom_index)
         IF (ncol == 0) CYCLE
         IF (group_env%mepos /= 0) ALLOCATE (atom_rhs(SIZE(grid, 2), ncol))
         CALL group_env%bcast(atom_rhs, 0)
         CALL build_phi_on_sphere(bs_env, qs_kinds, grid, atom_index, &
                                  radius(atom_index), local_index, ngrid, phi, ao_atoms, ao_col_start)
         DEALLOCATE (ao_atoms, ao_col_start)
         ALLOCATE (local_rhs(ngrid, ncol))
         DO i = 1, ngrid
            local_rhs(i, :) = atom_rhs(local_index(i), :)
         END DO
         DEALLOCATE (atom_rhs)
         CALL solve_ri_rs_fit(phi, local_rhs, bs_env%ri_rs%tikhonov, group_env, blacs_env)
         IF (group_env%mepos == 0) THEN
            CALL store_Z_lP_columns(mat_z, local_rhs, local_index, ngrid, ncol, atom_index, &
                                    row_size, row_offset, 0.0_dp)
         END IF
         DEALLOCATE (local_rhs, phi, local_index)
      END DO
      CALL cp_blacs_env_release(blacs_env)
      CALL group_env%free()
      DEALLOCATE (group_env)
      CALL timestop(handle)
   END SUBROUTINE fit_distributed_auto_ri_z_lp

! **************************************************************************************************
!> \brief Splits the atoms of the Z_lP solve into a single-rank list ("small", Phase A: LAPACK
!>        dpotrf/dpotrs on one rank) and a distributed list ("big", Phase B: ScaLAPACK
!>        pdpotrf/pdpotrs over a rank subgroup of size G), and sizes G.
!>        Each atom's single-rank peak memory is estimated as
!>          peak(P) = 8*n_local_grid(P)^2                (dense conditioned Gram matrix)
!>                  + 8*n_local_grid(P)*n_ao_used(P)     (phi_local)
!>                  + 8*n_local_grid(P)*n_RI(P)*(1+n_threads)   (d_lp + OMP partials)
!>        and send atoms whose peak exceeds mem_safety * available-memory-per-proc to the
!>        distributed path; G is auto-sized so the biggest atom's distributed Gram matrix slice
!>        fits alongside the replicated phi_local + d_lp. G is floored by the ScaLAPACK 32-bit
!>        index limit (a local block-cyclic slice of ~n_local_grid^2/G elements must stay below
!>        2^31 or pdpotrf segfaults).
!>        Records and reports the memory estimates before the standard or AUTO_RI solves.
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param cutoff_ri_per_atom ...
!> \param ri_blk_sizes per-atom ...
!> \param n_local_grid_atom ...
!> \param small_list ...
!> \param n_small ...
!> \param big_list ...
!> \param n_big ...
!> \param G ...
!> \param periodic_images use minimum-image distances for the Gamma integration domain
! **************************************************************************************************
   SUBROUTINE classify_z_lp_atoms(bs_env, ri_rs_grid_points, cutoff_ri_per_atom, ri_blk_sizes, &
                                  n_local_grid_atom, small_list, n_small, big_list, n_big, G, periodic_images)

!$    USE OMP_LIB, ONLY: omp_get_max_threads

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cutoff_ri_per_atom
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ri_blk_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: n_local_grid_atom, small_list, big_list
      INTEGER, INTENT(OUT)                               :: n_small, n_big, G
      LOGICAL, INTENT(IN), OPTIONAL :: periodic_images
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'classify_z_lp_atoms'

      INTEGER                                            :: handle

      ! Conservative fraction of measured available memory usable per rank for the Z_lP
      REAL(KIND=dp), PARAMETER                           :: mem_safety = 0.8_dp

      ! ScaLAPACK/BLACS index the per-rank local block-cyclic slice (~n_local_grid^2/G
      ! elements) with 32-bit integers; keep it safely below 2^31 or pdpotrf segfaults.
      REAL(KIND=dp), PARAMETER                           :: scalapack_loc_limit = 2.0E9_dp

      INTEGER                                            :: G_atom, G_int32, l, n_ao_used_atom, &
                                                            n_grid_total, n_local_grid, natom, &
                                                            nthreads_cls, P_loop_atom
      REAL(KIND=dp)                                      :: budget_bytes, cutoff_ri, dlp_bytes, &
                                                            mem_avail_GB, ng, nri, peak_bytes, &
                                                            phi_bytes
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P, displacement
      LOGICAL :: sum_images
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      sum_images = .FALSE.
      IF (PRESENT(periodic_images)) sum_images = periodic_images

      para_env => bs_env%para_env
      particle_set => bs_env%ri_rs%particle_set
      natom = bs_env%n_atom
      n_grid_total = bs_env%ri_rs%n_grid_points
      CPASSERT(SIZE(ri_rs_grid_points, 2) == n_grid_total)

      ! Per-atom sphere size: n_local_grid(P) = number of l with |r_l - R_P| <= cutoff_ri(P).
      ! It sets both the Gram matrix memory footprint (n_local_grid^2) and the solve
      ! cost (~n_local_grid^3), so it drives classification and the LPT load balancing.
      ALLOCATE (n_local_grid_atom(natom), small_list(natom), big_list(natom))
      DO P_loop_atom = 1, natom
         pos_P(:) = particle_set(P_loop_atom)%r(:)
         cutoff_ri = cutoff_ri_per_atom(P_loop_atom)
         n_local_grid = 0
         DO l = 1, n_grid_total
            displacement = ri_rs_grid_points(:, l) - pos_P
            IF (sum_images) displacement = pbc(ri_rs_grid_points(:, l), pos_P, bs_env%ri_rs%cell)
            IF (NORM2(displacement) <= cutoff_ri) THEN
               n_local_grid = n_local_grid + 1
            END IF
         END DO
         n_local_grid_atom(P_loop_atom) = n_local_grid
      END DO

      nthreads_cls = 1
!$    nthreads_cls = omp_get_max_threads()
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)  ! collective over all ranks
      budget_bytes = mem_safety*mem_avail_GB*1.0E9_dp

      n_small = 0
      n_big = 0
      G = 1
      G_atom = 1              ! max G a big atom needs (memory + ScaLAPACK int32 floor)
      IF (mem_avail_GB <= 0.0_dp .AND. bs_env%unit_nr > 0) THEN
         CPWARN("RI-RS Z_lP: no meminfo; single-rank solve for all atoms")
      END IF
      DO P_loop_atom = 1, natom
         ng = REAL(n_local_grid_atom(P_loop_atom), dp)
         nri = REAL(ri_blk_sizes(P_loop_atom), dp)
         CALL get_n_ao_in_sphere(bs_env, P_loop_atom, &
                                 cutoff_ri_per_atom(P_loop_atom), n_ao_used_atom, sum_images)
         CALL record_fit_memory(bs_env%ri_rs, n_local_grid_atom(P_loop_atom), &
                                n_ao_used_atom, ri_blk_sizes(P_loop_atom), peak_bytes)
         IF (mem_avail_GB <= 0.0_dp .OR. peak_bytes <= budget_bytes) THEN
            n_small = n_small + 1
            small_list(n_small) = P_loop_atom
         ELSE
            n_big = n_big + 1
            big_list(n_big) = P_loop_atom
            phi_bytes = 8.0_dp*ng*REAL(n_ao_used_atom, dp)
            dlp_bytes = 8.0_dp*ng*nri*REAL(1 + nthreads_cls, dp)
            ! Bound both memory and ScaLAPACK's 32-bit local indexing.
            G_int32 = CEILING(ng*ng/scalapack_loc_limit)
            G_atom = MAX(G_atom, G_int32, &
                         CEILING(8.0_dp*ng*ng/MAX(budget_bytes - phi_bytes - dlp_bytes, 1.0_dp)))
         END IF
      END DO
      CALL print_fit_memory_estimate(bs_env, mem_avail_GB)
      IF (n_big > 0) THEN
         IF (G_atom > para_env%num_pe) THEN
            CALL cp_abort(__LOCATION__, &
                          "RI-RS Z_lP: an atom is too large to fit even when "// &
                          "distributed over all ranks. Add nodes, use fewer MPI ranks "// &
                          "per node, lower CUTOFF_RADIUS_RI, or raise EPS_FILTER "// &
                          "for more grid screening.")
         END IF
         G = MIN(MAX(G_atom, 2), para_env%num_pe)
      END IF

      CALL timestop(handle)

   END SUBROUTINE classify_z_lp_atoms

! **************************************************************************************************
!> \brief Number of AO basis functions that can be non-zero inside the RI-RS integration sphere
!> \param bs_env ...
!> \param atom_P ...
!> \param cutoff_ri ...
!> \param n_ao_used ...
!> \param periodic_images use minimum-image distances to the Gamma integration sphere
! **************************************************************************************************
   SUBROUTINE get_n_ao_in_sphere(bs_env, atom_P, cutoff_ri, n_ao_used, periodic_images)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom_P
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_ri
      INTEGER, INTENT(OUT)                               :: n_ao_used
      LOGICAL, INTENT(IN)                                :: periodic_images

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_n_ao_in_sphere'

      INTEGER                                            :: handle, ri_atom
      REAL(KIND=dp)                                      :: displacement(3)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      particle_set => bs_env%ri_rs%particle_set
      n_ao_used = 0
      DO ri_atom = 1, bs_env%n_atom
         displacement = particle_set(ri_atom)%r - particle_set(atom_P)%r
         IF (periodic_images) displacement = pbc(particle_set(ri_atom)%r, particle_set(atom_P)%r, bs_env%ri_rs%cell)
         IF (NORM2(displacement) > &
             bs_env%ri_rs%radius_ao_per_atom(ri_atom) + cutoff_ri) CYCLE
         n_ao_used = n_ao_used + bs_env%i_ao_end_from_atom(ri_atom) - &
                     bs_env%i_ao_start_from_atom(ri_atom) + 1
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_n_ao_in_sphere

! **************************************************************************************************
!> \brief Builds the sphere-local AO matrix phi_local(l, μ) = ϕ_μ(r_l) for one RI atom P
!> \param bs_env ...
!> \param qs_kind_set ...
!> \param ri_rs_grid_points ...
!> \param atom_P ...
!> \param cutoff_ri ...
!> \param local_grid_idx ...
!> \param n_local_grid ...
!> \param phi_local ...
!> \param ao_atoms global indices of the active AO atoms
!> \param ao_col_start local phi column starts, including the terminal column offset
!> \param center ...
!> \param periodic_images select Gamma image-summed AO values and minimum-image sphere distances
!> \param ao_nonzero exact nonzero AO support on the retained grid rows
! **************************************************************************************************
   SUBROUTINE build_phi_on_sphere(bs_env, qs_kind_set, ri_rs_grid_points, &
                                  atom_P, cutoff_ri, local_grid_idx, n_local_grid, &
                                  phi_local, ao_atoms, ao_col_start, center, periodic_images, ao_nonzero)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ri_rs_grid_points
      INTEGER, INTENT(IN)                                :: atom_P
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_ri
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: local_grid_idx
      INTEGER, INTENT(OUT)                               :: n_local_grid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: phi_local
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ao_atoms, ao_col_start
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: center
      LOGICAL, INTENT(IN), OPTIONAL                      :: periodic_images
      LOGICAL, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: ao_nonzero(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_phi_on_sphere'

      INTEGER                                            :: col_end, col_start, handle, iactive, j, &
                                                            k, l, loc_idx, n_active, n_ao_used, &
                                                            n_grid_total, n_keep, ri_atom
      INTEGER, ALLOCATABLE                               :: atom_candidates(:)
      LOGICAL                                            :: sum_images
      REAL(KIND=dp)                                      :: d_sP, dist, r2_threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w_pt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: phi_keep, sphere_grid
      REAL(KIND=dp), DIMENSION(3)                        :: displacement, pos_P
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      sum_images = .FALSE.
      IF (PRESENT(periodic_images)) sum_images = periodic_images

      cell => bs_env%ri_rs%cell
      particle_set => bs_env%ri_rs%particle_set
      ! AUTO_RI may pass only the global grid points needed for the current atom block.
      n_grid_total = SIZE(ri_rs_grid_points, 2)
      IF (PRESENT(center)) THEN
         pos_P(:) = center(:)
      ELSE
         pos_P(:) = particle_set(atom_P)%r(:)
      END IF

      n_local_grid = 0
      DO l = 1, n_grid_total
         displacement = ri_rs_grid_points(:, l) - pos_P
         IF (sum_images) displacement = pbc(ri_rs_grid_points(:, l), pos_P, cell)
         dist = NORM2(displacement)
         IF (dist <= cutoff_ri) n_local_grid = n_local_grid + 1
      END DO

      ALLOCATE (local_grid_idx(n_local_grid))

      n_local_grid = 0
      DO l = 1, n_grid_total
         displacement = ri_rs_grid_points(:, l) - pos_P
         IF (sum_images) displacement = pbc(ri_rs_grid_points(:, l), pos_P, cell)
         dist = NORM2(displacement)
         IF (dist <= cutoff_ri) THEN
            n_local_grid = n_local_grid + 1
            local_grid_idx(n_local_grid) = l
         END IF
      END DO

      ALLOCATE (sphere_grid(3, n_local_grid))
      DO loc_idx = 1, n_local_grid
         sphere_grid(:, loc_idx) = ri_rs_grid_points(:, local_grid_idx(loc_idx))
      END DO

      ! Keep atoms in global order; phi columns are contiguous within each active atom.
      ALLOCATE (atom_candidates(bs_env%n_atom))
      n_active = 0
      DO ri_atom = 1, bs_env%n_atom
         displacement = particle_set(ri_atom)%r - pos_P
         IF (sum_images) displacement = pbc(particle_set(ri_atom)%r, pos_P, cell)
         d_sP = NORM2(displacement)
         IF (d_sP > bs_env%ri_rs%radius_ao_per_atom(ri_atom) + cutoff_ri) CYCLE
         IF (bs_env%i_ao_end_from_atom(ri_atom) < bs_env%i_ao_start_from_atom(ri_atom)) CYCLE
         n_active = n_active + 1
         atom_candidates(n_active) = ri_atom
      END DO
      ALLOCATE (ao_atoms(n_active), SOURCE=atom_candidates(:n_active))
      DEALLOCATE (atom_candidates)
      ALLOCATE (ao_col_start(n_active + 1))
      ao_col_start(1) = 1
      DO iactive = 1, n_active
         ri_atom = ao_atoms(iactive)
         ao_col_start(iactive + 1) = ao_col_start(iactive) + &
                                     bs_env%i_ao_end_from_atom(ri_atom) - bs_env%i_ao_start_from_atom(ri_atom) + 1
      END DO
      n_ao_used = ao_col_start(n_active + 1) - 1
      ALLOCATE (phi_local(n_local_grid, n_ao_used), SOURCE=0.0_dp)

      DO iactive = 1, n_active
         ri_atom = ao_atoms(iactive)
         col_start = ao_col_start(iactive)
         col_end = ao_col_start(iactive + 1) - 1
         ! A positive CUTOFF_RADIUS_RI_AO overrides the per-atom Gaussian radius
         ! with a user-defined hard cutoff.
         IF (.NOT. sum_images .AND. bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) THEN
            r2_threshold = bs_env%ri_rs%cutoff_radius_ri_ao**2
         ELSE
            r2_threshold = bs_env%ri_rs%radius_ao_per_atom(ri_atom)**2
         END IF

         IF (sum_images) THEN
            CALL evaluate_ao_basis_on_points(phi_local(:, col_start:col_end), sphere_grid, &
                                             bs_env%basis_set_AO( &
                                             particle_set(ri_atom)%atomic_kind%kind_number)%gto_basis_set, &
                                             particle_set(ri_atom)%r, cell, &
                                             cutoff_squared=r2_threshold, periodic_images=.TRUE.)
         ELSE
            CALL evaluate_ao_on_points(phi_local(:, col_start:col_end), sphere_grid, &
                                       ri_atom, particle_set, qs_kind_set, cell, &
                                       cutoff_squared=r2_threshold)
         END IF
      END DO

      DEALLOCATE (sphere_grid)

      IF (.NOT. sum_images .AND. n_local_grid > 0) THEN
         ALLOCATE (w_pt(n_local_grid))
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP SHARED(n_local_grid, n_ao_used, phi_local, w_pt) &
         !$OMP PRIVATE(l, j) SCHEDULE(STATIC)
         DO l = 1, n_local_grid
            w_pt(l) = 0.0_dp
            DO j = 1, n_ao_used
               w_pt(l) = MAX(w_pt(l), ABS(phi_local(l, j)))
            END DO
         END DO
         !$OMP END PARALLEL DO
         n_keep = COUNT(w_pt > bs_env%eps_filter)
         IF (n_keep < n_local_grid) THEN
            ALLOCATE (phi_keep(n_keep, n_ao_used))
            k = 0
            DO l = 1, n_local_grid
               IF (w_pt(l) > bs_env%eps_filter) THEN
                  k = k + 1
                  phi_keep(k, :) = phi_local(l, :)
                  local_grid_idx(k) = local_grid_idx(l)
               END IF
            END DO
            CALL MOVE_ALLOC(phi_keep, phi_local)
            n_local_grid = n_keep
         END IF
         DEALLOCATE (w_pt)
      END IF

      IF (PRESENT(ao_nonzero)) THEN
         ! Build support on the final retained rows, once for all RHS contractions using phi.
         ALLOCATE (ao_nonzero(n_local_grid, n_active), SOURCE=.FALSE.)
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(n_active, n_local_grid, ao_col_start, phi_local, ao_nonzero) &
         !$OMP PRIVATE(iactive, j, l)
         DO iactive = 1, n_active
            DO j = ao_col_start(iactive), ao_col_start(iactive + 1) - 1
               DO l = 1, n_local_grid
                  ao_nonzero(l, iactive) = ao_nonzero(l, iactive) .OR. phi_local(l, j) /= 0.0_dp
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_phi_on_sphere

! **************************************************************************************************
!> \brief Computes a three-center integral block directly in AUTO_RI contraction columns.
!>
!>        For each angular momentum l, the primitive RI Gaussians on the atoms contributing
!>        to one optimized function are collected before the contraction
!>
!>          (μν|p) = Σ_P (μν|P) U_Pp .
!>
!>        This avoids constructing and retaining the complete reference-basis block (μν|P).
!> \param int_3c ...
!> \param ctx ...
!> \param ws ...
!> \param atom_j ...
!> \param atom_k ...
!> \param atom_i ...
!> \param transform ...
!> \param transform_row ...
!> \param screened ...
! **************************************************************************************************
   SUBROUTINE build_3c_integral_block_auto_ri_ctx(int_3c, ctx, ws, atom_j, atom_k, atom_i, &
                                                  transform, transform_row, screened)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: int_3c
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(gw_3c_ws_type), INTENT(INOUT)                 :: ws
      INTEGER, INTENT(IN)                                :: atom_j, atom_k, atom_i
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: transform
      INTEGER, INTENT(IN)                                :: transform_row
      LOGICAL, INTENT(OUT), OPTIONAL                     :: screened

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_integral_block_auto_ri_ctx'

      INTEGER :: handle_contract, handle_eri, ikind, iset, jkind, jset, kkind, kset, l, ncoi, &
         ncoj, ncok, ncol, npgf_group, nseti, nsetj, nsetk, primitive_first, sgfi, sgfj, sgfk
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      REAL(KIND=dp)                                      :: dij, dik, djk, group_radius, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k, sijk_ext
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rpgf_group, zet_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: spi_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: sijk, sijk_contr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rik, rj, rjk, rk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, zeti, zetj, zetk

      IF (PRESENT(screened)) screened = .FALSE.
      ri = pbc(ctx%particle_set(atom_i)%r(1:3), ctx%cell)
      rj = pbc(ctx%particle_set(atom_j)%r(1:3), ctx%cell)
      rk = pbc(ctx%particle_set(atom_k)%r(1:3), ctx%cell)
      rjk = rk - rj
      rij = rj - ri
      rik = rk - ri
      djk = NORM2(rjk)
      dij = NORM2(rij)
      dik = NORM2(rik)

      ikind = ctx%kind_of(atom_i)
      jkind = ctx%kind_of(atom_j)
      kkind = ctx%kind_of(atom_k)
      CALL get_gto_basis_set(ctx%basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, &
                             lmax=lmax_i, lmin=lmin_i, npgf=npgfi, nset=nseti, &
                             nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                             zet=zeti, kind_radius=kind_radius_i)
      CALL get_gto_basis_set(ctx%basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, &
                             lmax=lmax_j, lmin=lmin_j, npgf=npgfj, nset=nsetj, &
                             nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                             zet=zetj, kind_radius=kind_radius_j)
      CALL get_gto_basis_set(ctx%basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, &
                             lmax=lmax_k, lmin=lmin_k, npgf=npgfk, nset=nsetk, &
                             nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                             zet=zetk, kind_radius=kind_radius_k)

      IF (kind_radius_j + kind_radius_i + ctx%dr_ij < dij .OR. &
          kind_radius_j + kind_radius_k + ctx%dr_jk < djk .OR. &
          kind_radius_k + kind_radius_i + ctx%dr_ik < dik) THEN
         IF (PRESENT(screened)) screened = .TRUE.
         RETURN
      END IF

      ncol = SIZE(transform, 2)
      CPASSERT(SIZE(int_3c, 3) == ncol)
      DO l = 0, ctx%maxli
         npgf_group = 0
         group_radius = 0.0_dp
         DO iset = 1, nseti
            IF (lmin_i(iset) /= l .OR. lmax_i(iset) /= l) CYCLE
            npgf_group = npgf_group + npgfi(iset)
            group_radius = MAX(group_radius, set_radius_i(iset))
         END DO
         IF (npgf_group == 0) CYCLE
         ncoi = npgf_group*ncoset(l)
         ALLOCATE (zet_group(npgf_group), rpgf_group(npgf_group))
         ALLOCATE (spi_group(ncoi, ncol), source=0.0_dp)
         primitive_first = 1
         DO iset = 1, nseti
            IF (lmin_i(iset) /= l .OR. lmax_i(iset) /= l) CYCLE
            zet_group(primitive_first:primitive_first + npgfi(iset) - 1) = &
               zeti(1:npgfi(iset), iset)
            rpgf_group(primitive_first:primitive_first + npgfi(iset) - 1) = &
               rpgf_i(1:npgfi(iset), iset)
            sgfi = first_sgf_i(1, iset)
            spi_group((primitive_first - 1)*ncoset(l) + 1: &
                      (primitive_first + npgfi(iset) - 1)*ncoset(l), :) = &
               MATMUL(ctx%spi(iset, ikind)%array, &
                      transform(transform_row + sgfi - 1: &
                                transform_row + sgfi + nsgfi(iset) - 2, :))
            primitive_first = primitive_first + npgfi(iset)
         END DO

         DO jset = 1, nsetj
            IF (set_radius_j(jset) + group_radius + ctx%dr_ij < dij) CYCLE
            DO kset = 1, nsetk
               IF (set_radius_j(jset) + set_radius_k(kset) + ctx%dr_jk < djk) CYCLE
               IF (set_radius_k(kset) + group_radius + ctx%dr_ik < dik) CYCLE
               ncoj = npgfj(jset)*ncoset(lmax_j(jset))
               ncok = npgfk(kset)*ncoset(lmax_k(kset))
               sgfj = first_sgf_j(1, jset)
               sgfk = first_sgf_k(1, kset)
               IF (ncoj*ncok*ncoi <= 0) CYCLE
               ALLOCATE (sijk(ncoj, ncok, ncoi), source=0.0_dp)
               CALL timeset(routineN//'_eri', handle_eri)
               CALL eri_3center(sijk, &
                                lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), &
                                rpgf_j(:, jset), rj, &
                                lmin_k(kset), lmax_k(kset), npgfk(kset), zetk(:, kset), &
                                rpgf_k(:, kset), rk, l, l, npgf_group, zet_group, &
                                rpgf_group, ri, djk, dij, dik, ws%lib, ctx%potential_parameter, &
                                int_abc_ext=sijk_ext)
               CALL timestop(handle_eri)
               ALLOCATE (sijk_contr(nsgfj(jset), nsgfk(kset), ncol))
               CALL timeset(routineN//'_contract', handle_contract)
               CALL abc_contract_xsmm(sijk_contr, sijk, ctx%tspj(jset, jkind)%array, &
                                      ctx%spk(kset, kkind)%array, spi_group, ncoj, ncok, ncoi, &
                                      nsgfj(jset), nsgfk(kset), ncol, ws%cpp_buffer, ws%ccp_buffer)
               CALL timestop(handle_contract)
               DEALLOCATE (sijk)
               int_3c(sgfj:sgfj + nsgfj(jset) - 1, sgfk:sgfk + nsgfk(kset) - 1, :) = &
                  int_3c(sgfj:sgfj + nsgfj(jset) - 1, sgfk:sgfk + nsgfk(kset) - 1, :) + &
                  sijk_contr
               DEALLOCATE (sijk_contr)
            END DO
         END DO
         DEALLOCATE (zet_group, rpgf_group, spi_group)
      END DO

   END SUBROUTINE build_3c_integral_block_auto_ri_ctx

! **************************************************************************************************
!> \brief Returns the reference RI basis size of one atom.
!> \param bs_env ...
!> \param iatom ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION get_ref_RI_size(bs_env, iatom) RESULT(n)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom

      INTEGER                                            :: ikind

      ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
      n = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf

   END FUNCTION get_ref_RI_size

! **************************************************************************************************
!> \brief Record the dense Gram, phi and threaded RHS memory model for an RI-RS fit.
!>        The single-rank estimate is an upper bound for a distributed Gram solve.
!> \param ri_rs owned fitting statistics, reset with the representation
!> \param n_grid number of fitting grid rows, before screening when obtained by classification
!> \param n_ao number of AO columns
!> \param n_ri number of RHS columns
!> \param peak_bytes optional single-rank estimate used by solver classification
! **************************************************************************************************
   SUBROUTINE record_fit_memory(ri_rs, n_grid, n_ao, n_ri, peak_bytes)
!$    USE OMP_LIB, ONLY: omp_get_max_threads
      TYPE(ri_rs_env), INTENT(INOUT)                     :: ri_rs
      INTEGER, INTENT(IN)                                :: n_grid, n_ao, n_ri
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: peak_bytes

      INTEGER                                            :: nthreads
      REAL(KIND=dp)                                      :: ng, peak

      nthreads = 1
!$    nthreads = omp_get_max_threads()
      ng = REAL(n_grid, dp)
      peak = 8.0_dp*ng*ng + 8.0_dp*ng*REAL(n_ao, dp) + &
             8.0_dp*ng*REAL(n_ri, dp)*REAL(1 + nthreads, dp)
      ri_rs%fit_max_grid = MAX(ri_rs%fit_max_grid, n_grid)
      ri_rs%fit_peak_bytes = MAX(ri_rs%fit_peak_bytes, peak)
      IF (PRESENT(peak_bytes)) peak_bytes = peak
   END SUBROUTINE record_fit_memory

! **************************************************************************************************
!> \brief Prints the fitting memory bound before allocating the local fitting RHS and Gram matrix.
!> \param bs_env recorded fitting sizes and output unit
!> \param mem_avail_GB available memory per rank, measured collectively by the caller
! **************************************************************************************************
   SUBROUTINE print_fit_memory_estimate(bs_env, mem_avail_GB)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: mem_avail_GB

      CHARACTER(LEN=max_line_length)                     :: msg
      REAL(KIND=dp)                                      :: mem_gram_GB, mem_Zlp_peak_GB

      mem_gram_GB = 8.0_dp*REAL(bs_env%ri_rs%fit_max_grid, dp)**2*1.0E-9_dp
      mem_Zlp_peak_GB = bs_env%ri_rs%fit_peak_bytes*1.0E-9_dp
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS fitting memory estimate per MPI process:'
         WRITE (bs_env%unit_nr, '(T4,A,F37.2,A)') &
            'Available memory per process (system)', mem_avail_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F17.2,A)') &
            'Required for Z_lP fit (Gram+phi+RHS, single-rank bound)', mem_Zlp_peak_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T69,I12)') &
            'Largest fitting-grid number of grid points:', bs_env%ri_rs%fit_max_grid
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'Largest dense Gram matrix:', mem_gram_GB, ' GB'
         CALL m_flush(bs_env%unit_nr)
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Zlp_peak_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the Z_lP solve, ", mem_Zlp_peak_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, &
            " GB. The dense Gram matrix scales as n_local_grid^2; "// &
            "the estimate assumes a single-rank solve. Localized fits may distribute the Gram "// &
            "matrix over a ScaLAPACK subgroup. Use fewer MPI ranks per node if more memory "// &
            "per rank is needed."
         CPWARN(TRIM(msg))
      END IF
   END SUBROUTINE print_fit_memory_estimate

END MODULE ri_rs_setup
