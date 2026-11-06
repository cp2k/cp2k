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
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_binary_read, dbcsr_binary_write, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_info, dbcsr_type, dbcsr_type_no_symmetry
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
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
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
   USE util,                            ONLY: sort
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
!>        With AUTO_RI, the optimized functions of atom P are φ_p = Σ_{P∈A} U_Pp φ_P, so step 4
!>        uses the reference functions of P and d_lp = Σ_P d_lP U_Pp.
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
   SUBROUTINE compute_Z_lP(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP, periodic_images)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP
      LOGICAL, INTENT(IN)                                :: periodic_images

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'compute_Z_lP'

      INTEGER :: atom_j_mepos, atom_j_stride, atom_P, G, handle, i_blk, iatom, idx, iphase, j, &
         max_ao_size, my_group, n_big, n_done, n_groups, n_loc_ri, n_local_grid, n_my_atoms, &
         n_ref_ri, n_small, natom, npcol_phi, num_grid_chunks, phase_hi
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ao_atoms, ao_col_start, big_list, local_grid_idx, &
         my_atoms_A, my_atoms_B, n_local_grid_atom, ref_blk_sizes, row_offset, small_list
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_ri, r_blk_sizes, ri_blk_sizes, &
                                                            row_dist_grid
      LOGICAL                                            :: do_scatter, use_dist
      LOGICAL, ALLOCATABLE                               :: ao_nonzero(:, :)
      REAL(KIND=dp)                                      :: cutoff_ri, item_start_time, r_c, t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoff_ri_per_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_local, d_lp_opt, phi_local, U_Pp
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
      bs_env%ri_rs%fit_max_grid = 0
      bs_env%ri_rs%fit_peak_bytes = 0.0_dp

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
      ! (μν|P) and d_lP use the reference RI functions; AUTO_RI contracts them to ri_blk_sizes.
      ALLOCATE (ref_blk_sizes(natom))
      ref_blk_sizes(:) = ri_blk_sizes
      IF (bs_env%auto_ri%enabled) ref_blk_sizes(:) = bs_env%auto_ri%sizes_ref_RI

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
         CALL classify_z_lp_atoms(bs_env, ri_rs_grid_points, cutoff_ri_per_atom, ref_blk_sizes, &
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
               IF (n_loc_ri == 0) CYCLE
               n_ref_ri = ref_blk_sizes(atom_P)
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
               ALLOCATE (d_lp_local(n_local_grid, n_ref_ri))
               d_lp_local = 0.0_dp

               CALL compute_d_lp(bs_env, ctx_3c, phi_local, ao_atoms, ao_col_start, ao_nonzero, d_lp_local, n_local_grid, &
                                 n_ref_ri, atom_P, max_ao_size, atom_j_mepos, atom_j_stride, periodic_images)

               DEALLOCATE (ao_atoms, ao_col_start, ao_nonzero)

               ! Reduce per-subgroup-rank partials into the replicated d_lp_local.
               ! Skipped for BLAS path: each rank has the full sum locally.
               IF (use_dist) THEN
                  CALL para_env_sub%sum(d_lp_local)
               END IF
               IF (bs_env%auto_ri%enabled) THEN
                  CALL get_U_Pp_AA(bs_env%auto_ri, atom_P, U_Pp)
                  ALLOCATE (d_lp_opt(n_local_grid, n_loc_ri))
                  CALL dgemm('N', 'N', n_local_grid, n_loc_ri, n_ref_ri, 1.0_dp, d_lp_local, &
                             n_local_grid, U_Pp, n_ref_ri, 0.0_dp, d_lp_opt, n_local_grid)
                  CALL MOVE_ALLOC(d_lp_opt, d_lp_local)
                  DEALLOCATE (U_Pp)
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

      DEALLOCATE (row_offset, ri_blk_sizes, ref_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)

      CALL timestop(handle)

   END SUBROUTINE compute_Z_lP

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
!> \brief Unpacks the AA contraction matrix of one atom.
!> \param auto_ri ...
!> \param atom ...
!> \param U_Pp ...
! **************************************************************************************************
   SUBROUTINE get_U_Pp_AA(auto_ri, atom, U_Pp)
      TYPE(auto_ri_type), INTENT(IN)                     :: auto_ri
      INTEGER, INTENT(IN)                                :: atom
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: U_Pp(:, :)

      INTEGER                                            :: first, ncolumn, nrow

      nrow = auto_ri%sizes_ref_RI(atom)
      ncolumn = auto_ri%sizes_opt_RI(atom)
      first = auto_ri%U_Pp_atom_offset(atom)
      ALLOCATE (U_Pp(nrow, ncolumn))
      U_Pp(:, :) = RESHAPE(auto_ri%U_Pp_atom(first:first + nrow*ncolumn - 1), [nrow, ncolumn])

   END SUBROUTINE get_U_Pp_AA

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
!>        Records and reports the memory estimates before the Z_lP solves.
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
