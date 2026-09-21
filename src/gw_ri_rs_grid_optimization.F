!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Optimize automatically initialized atom-centred RI-RS grids.
!> \par History
!>      09.2026 created [Jan Wilhelm]
! **************************************************************************************************
MODULE gw_ri_rs_grid_optimization
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_lbfgs,                        ONLY: setulb
   USE gw_ri_rs_grid_initialization,    ONLY: initialize_ri_rs_grid
   USE gw_ri_rs_utils,                  ONLY: evaluate_ao_basis_on_points,&
                                              filter_grid_to_voronoi,&
                                              get_rirs_cluster_atoms
   USE gw_utils_compute_integrals,      ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: optimize_ri_rs_grid

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_grid_optimization'

   TYPE, PRIVATE :: local_cluster_3c_integrals_type
      INTEGER, ALLOCATABLE                 :: atom_indices(:)
      REAL(KIND=dp), ALLOCATABLE           :: Int_3c(:, :, :)
   END TYPE local_cluster_3c_integrals_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize Lebedev grids and subsequently optimize their grid-point coordinates.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE optimize_ri_rs_grid(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'optimize_ri_rs_grid'

      INTEGER                                            :: handle
      TYPE(local_cluster_3c_integrals_type), ALLOCATABLE :: cluster_3c_int(:)

      CALL timeset(routineN, handle)

      ! setup Lebedev initial grid inside Voronoi volume, select points with Cholesky decomp.
      CALL initialize_ri_rs_grid(bs_env)

      ! Validate the initialized RI-RS state and report the start of the optimization.
      CALL prepare_ri_rs_grid_optimization(bs_env)

      ! C_A = {B: |R_A - R_B| < R_cut}; store (μν|P) for μ, ν, and P centered in C_A.
      CALL build_cluster_3c_int(bs_env, cluster_3c_int)

      ! E_A = Σ_{μνP∈C_A} [(μν|P) - Σ_{l∈G_A} ϕ_μ(r_l) ϕ_ν(r_l) Z_lP^(A)]².
      ! Minimize E_loc(norm) = N_atom^(-1) Σ_A [E_A / Σ_{μνP∈C_A} |(μν|P)|²].
      CALL optimize_grid_coordinates(bs_env, cluster_3c_int)

      CALL timestop(handle)

   END SUBROUTINE optimize_ri_rs_grid

! **************************************************************************************************
!> \brief Validate initialized RI-RS state and announce the grid optimization.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE prepare_ri_rs_grid_optimization(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: iatom, ikind, unit_nr

      CPASSERT(ASSOCIATED(bs_env%ri_rs%cell))
      CPASSERT(ASSOCIATED(bs_env%ri_rs%particle_set))
      CPASSERT(ASSOCIATED(bs_env%para_env))
      CPASSERT(ALLOCATED(bs_env%basis_set_AO))
      CPASSERT(ALLOCATED(bs_env%basis_set_RI))
      CPASSERT(ALLOCATED(bs_env%sizes_AO))
      CPASSERT(ALLOCATED(bs_env%ri_rs%atomic_grids))
      CPASSERT(SIZE(bs_env%sizes_AO) == bs_env%n_atom)
      CPASSERT(SIZE(bs_env%ri_rs%atomic_grids) == bs_env%n_atom)

      DO iatom = 1, bs_env%n_atom
         CPASSERT(ALLOCATED(bs_env%ri_rs%atomic_grids(iatom)%raw_points))
         CPASSERT(SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 1) == 3)
         CPASSERT(SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2) > 0)
         ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(bs_env%basis_set_AO))
         CPASSERT(ikind <= SIZE(bs_env%basis_set_RI))
         CPASSERT(ASSOCIATED(bs_env%basis_set_AO(ikind)%gto_basis_set))
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(ikind)%gto_basis_set))
         IF (bs_env%sizes_AO(iatom) < 1 .OR. &
             bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf < 1) THEN
            CPABORT("Every atom in GRID_OPTIMIZATION needs ORB and RI_AUX functions")
         END IF
      END DO

      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'Started RI-RS grid optimization'
         CALL m_flush(unit_nr)
      END IF
   END SUBROUTINE prepare_ri_rs_grid_optimization

! **************************************************************************************************
!> \brief Evaluate the normalized local objective E_loc and its gradient for one L-BFGS trial vector.
!> \param grid_coordinates Flattened atom-relative grid coordinates.
!> \param atom_coordinate_offsets Starting coordinate offset for each atom.
!> \param bs_env ...
!> \param cluster_3c_int Rank-local cluster three-centre integrals.
!> \param normalized_error Mean normalized squared three-centre-integral error.
!> \param coordinate_gradient Derivative of normalized_error with respect to coordinates.
!> \param max_abs_error Largest absolute three-centre-integral error.
!> \param fit_successful Whether every local-cluster evaluation succeeded.
! **************************************************************************************************
   SUBROUTINE grid_objective(grid_coordinates, atom_coordinate_offsets, bs_env, cluster_3c_int, &
                             normalized_error, coordinate_gradient, &
                             max_abs_error, fit_successful)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: grid_coordinates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_coordinate_offsets
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(local_cluster_3c_integrals_type), &
         DIMENSION(:), INTENT(IN)                        :: cluster_3c_int
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coordinate_gradient
      REAL(KIND=dp), INTENT(OUT)                         :: max_abs_error
      LOGICAL, INTENT(OUT)                               :: fit_successful

      INTEGER                                            :: evaluated_cluster_count, icluster, &
                                                            successful_ranks
      LOGICAL                                            :: cluster_fit_successful
      REAL(KIND=dp)                                      :: cluster_max_abs_error, &
                                                            cluster_normalized_error
      REAL(KIND=dp), ALLOCATABLE :: cluster_coordinate_gradient(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      para_env => bs_env%para_env
      CPASSERT(SIZE(coordinate_gradient) == SIZE(grid_coordinates))
      CALL unpack_lbfgs_grid_coordinates(grid_coordinates, atom_coordinate_offsets, &
                                         bs_env%ri_rs%atomic_grids)
      normalized_error = 0.0_dp
      coordinate_gradient = 0.0_dp
      evaluated_cluster_count = 0
      max_abs_error = 0.0_dp
      fit_successful = .TRUE.
      DO icluster = 1, SIZE(cluster_3c_int)
         CALL evaluate_local_cluster(cluster_3c_int(icluster), bs_env, &
                                     cluster_normalized_error, cluster_coordinate_gradient, &
                                     cluster_max_abs_error, cluster_fit_successful)
         IF (.NOT. cluster_fit_successful) THEN
            fit_successful = .FALSE.
            IF (ALLOCATED(cluster_coordinate_gradient)) DEALLOCATE (cluster_coordinate_gradient)
            EXIT
         END IF
         normalized_error = normalized_error + cluster_normalized_error
         CALL accumulate_atom_gradient(cluster_3c_int(icluster)%atom_indices, &
                                       bs_env%ri_rs%atomic_grids, &
                                       atom_coordinate_offsets, cluster_coordinate_gradient, &
                                       coordinate_gradient)
         evaluated_cluster_count = evaluated_cluster_count + 1
         max_abs_error = &
            MAX(max_abs_error, cluster_max_abs_error)
         DEALLOCATE (cluster_coordinate_gradient)
      END DO
      successful_ranks = MERGE(1, 0, fit_successful)
      CALL para_env%sum(normalized_error)
      CALL para_env%sum(coordinate_gradient)
      CALL para_env%sum(evaluated_cluster_count)
      CALL para_env%sum(successful_ranks)
      CALL para_env%max(max_abs_error)
      fit_successful = evaluated_cluster_count == bs_env%n_atom .AND. &
                       successful_ranks == para_env%num_pe
      IF (.NOT. fit_successful) THEN
         normalized_error = HUGE(normalized_error)
         coordinate_gradient = 0.0_dp
         RETURN
      END IF
      ! E_loc(norm) = (1/N_atom) Σ_A [E_A / Σ_{μνP ∈ C_A} |(μν|P)|²].
      normalized_error = normalized_error/REAL(bs_env%n_atom, dp)
      coordinate_gradient = coordinate_gradient/REAL(bs_env%n_atom, dp)
   END SUBROUTINE grid_objective

! **************************************************************************************************
!> \brief Drive L-BFGS, repeatedly evaluating E_loc and its gradient, and return the best trial.
!> \param bs_env ...
!> \param cluster_3c_int Rank-local cluster three-centre integrals.
! **************************************************************************************************
   SUBROUTINE optimize_grid_coordinates(bs_env, cluster_3c_int)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(local_cluster_3c_integrals_type), &
         DIMENSION(:), INTENT(IN)                        :: cluster_3c_int

      INTEGER, PARAMETER                                 :: lbfgs_history = 7
      REAL(KIND=dp), PARAMETER                           :: lbfgs_factr = 0.0_dp, &
                                                            lbfgs_pgtol = 1.0E-9_dp

      CHARACTER(LEN=60)                                  :: line_search_state, optimizer_task
      INTEGER                                            :: evaluation_count, handle, iatom, &
                                                            molecular_outside, npoints, unit_nr
      INTEGER, ALLOCATABLE                               :: atom_coordinate_offsets(:), &
                                                            bound_types(:), integer_workspace(:)
      INTEGER, DIMENSION(44)                             :: integer_state
      LOGICAL                                            :: evaluation_successful, fit_successful, &
                                                            have_best
      LOGICAL, ALLOCATABLE                               :: inside(:)
      LOGICAL, DIMENSION(4)                              :: logical_state
      REAL(KIND=dp)                                      :: best_normalized_error, max_abs_error, &
                                                            normalized_error, start_time
      REAL(KIND=dp), ALLOCATABLE :: best_grid_coordinates(:), coordinate_gradient(:), &
         grid_coordinates(:), lower_bounds(:), points(:, :), upper_bounds(:), workspace(:)
      REAL(KIND=dp), DIMENSION(29)                       :: real_state

      CALL timeset("rirs_grid_LBFGS", handle)
      start_time = m_walltime()

      CALL pack_lbfgs_grid_coordinates(bs_env%ri_rs%atomic_grids, grid_coordinates, &
                                       atom_coordinate_offsets)
      CPASSERT(SIZE(grid_coordinates) > 0)

      ALLOCATE (best_grid_coordinates(SIZE(grid_coordinates)), &
                coordinate_gradient(SIZE(grid_coordinates)), &
                lower_bounds(SIZE(grid_coordinates)), &
                upper_bounds(SIZE(grid_coordinates)), &
                bound_types(SIZE(grid_coordinates)), &
                integer_workspace(3*SIZE(grid_coordinates)), &
                workspace(2*lbfgs_history*SIZE(grid_coordinates) + &
                          5*SIZE(grid_coordinates) + 11*lbfgs_history**2 + 8*lbfgs_history))
      ! L-BFGS-B ignores lower_bounds and upper_bounds when bound_types is zero.
      lower_bounds = 0.0_dp
      upper_bounds = 0.0_dp
      bound_types = 0
      optimizer_task = 'START'
      line_search_state = ''
      normalized_error = HUGE(normalized_error)
      coordinate_gradient = 0.0_dp
      workspace = 0.0_dp
      integer_workspace = 0
      logical_state = .FALSE.
      integer_state = 0
      real_state = 0.0_dp
      evaluation_count = 0
      have_best = .FALSE.
      best_normalized_error = HUGE(best_normalized_error)
      best_grid_coordinates(:) = grid_coordinates

      DO
         CALL setulb(SIZE(grid_coordinates), lbfgs_history, grid_coordinates, &
                     lower_bounds, upper_bounds, bound_types, &
                     normalized_error, coordinate_gradient, lbfgs_factr, lbfgs_pgtol, &
                     workspace, integer_workspace, optimizer_task, -1, line_search_state, &
                     logical_state, integer_state, real_state, -1.0_dp)
         IF (optimizer_task(1:2) == 'FG') THEN
            IF (evaluation_count >= bs_env%ri_rs%grid_opt%max_iter) EXIT
            CALL grid_objective(grid_coordinates, atom_coordinate_offsets, bs_env, cluster_3c_int, &
                                normalized_error, coordinate_gradient, &
                                max_abs_error, evaluation_successful)
            evaluation_count = evaluation_count + 1
            IF (.NOT. evaluation_successful) EXIT
            IF (evaluation_count == 1 .AND. bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T61,ES20.12)') &
                  'RI-RS grid initial normalized error', normalized_error
               CALL m_flush(bs_env%unit_nr)
            END IF
            IF (.NOT. have_best .OR. normalized_error < best_normalized_error) THEN
               have_best = .TRUE.
               best_normalized_error = normalized_error
               best_grid_coordinates(:) = grid_coordinates
            END IF
         ELSE IF (optimizer_task(1:5) == 'NEW_X') THEN
            CYCLE
         ELSE
            EXIT
         END IF
      END DO
      IF (have_best) grid_coordinates(:) = best_grid_coordinates
      CALL grid_objective(grid_coordinates, atom_coordinate_offsets, bs_env, cluster_3c_int, &
                          normalized_error, coordinate_gradient, &
                          max_abs_error, fit_successful)

      IF (.NOT. fit_successful) CPABORT("RI-RS grid optimization produced no regularized fit")

      DO iatom = 1, SIZE(bs_env%ri_rs%atomic_grids)
         bs_env%ri_rs%atomic_grids(iatom)%npts = &
            SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)
      END DO
      bs_env%ri_rs%Z_lP_exists = .FALSE.

      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         npoints = 0
         molecular_outside = 0
         DO iatom = 1, SIZE(bs_env%ri_rs%atomic_grids)
            ALLOCATE (points(3, bs_env%ri_rs%atomic_grids(iatom)%npts), &
                      inside(bs_env%ri_rs%atomic_grids(iatom)%npts))
            points(:, :) = bs_env%ri_rs%atomic_grids(iatom)%raw_points
            CALL filter_grid_to_voronoi(points, iatom, bs_env%ri_rs%particle_set, mask=inside)
            molecular_outside = molecular_outside + COUNT(.NOT. inside)
            npoints = npoints + SIZE(inside)
            DEALLOCATE (points, inside)
         END DO
         WRITE (unit_nr, '(T2,A,T69,F10.1,A)') &
            'RI-RS grid optimization completed, execution time:', m_walltime() - start_time, ' s'
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Normalized 3C error:', normalized_error
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Maximum absolute 3C error:', max_abs_error
         WRITE (unit_nr, '(T2,A,T69,I12,A,I0)') &
            'Optimized points outside molecular Voronoi cells:', molecular_outside, ' / ', npoints
         FLUSH (unit_nr)
      END IF
      CALL timestop(handle)
   END SUBROUTINE optimize_grid_coordinates

! **************************************************************************************************
!> \brief Build each C_A and store its exact (μν|P); this routine owns the integral context.
!> \param bs_env ...
!> \param cluster_3c_int Rank-local cluster three-centre integrals.
! **************************************************************************************************
   SUBROUTINE build_cluster_3c_int(bs_env, cluster_3c_int)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(local_cluster_3c_integrals_type), &
         ALLOCATABLE, INTENT(OUT)                        :: cluster_3c_int(:)

      INTEGER :: handle, i_cluster_atom_j, i_cluster_atom_k, i_cluster_atom_p, iatom, iatom_j, &
         iatom_k, iatom_p, icenter_atom, icluster_atom, ikind, ilocal_cluster, ip, &
         n_local_clusters, nAO_cluster, nRI_cluster, nRI_nonzero
      INTEGER, ALLOCATABLE                               :: AO_offset(:), RI_offset(:), &
                                                            sizes_ref_RI(:)
      LOGICAL                                            :: screened
      REAL(KIND=dp)                                      :: cluster_radius
      REAL(KIND=dp), ALLOCATABLE                         :: Int_3c_nonzero(:, :, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gw_3c_ctx_type)                               :: integral_context
      TYPE(gw_3c_ws_type)                                :: workspace
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset("rirs_grid_build_clusters", handle)

      cell => bs_env%ri_rs%cell
      para_env => bs_env%para_env
      particle_set => bs_env%ri_rs%particle_set
      cluster_radius = bs_env%ri_rs%grid_opt%cutoff_atomic_cluster
      ALLOCATE (sizes_ref_RI(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         sizes_ref_RI(iatom) = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
      END DO
      n_local_clusters = COUNT([(MOD(icenter_atom - 1, para_env%num_pe) == para_env%mepos, &
                                 icenter_atom=1, bs_env%n_atom)])
      ALLOCATE (cluster_3c_int(n_local_clusters), &
                AO_offset(bs_env%n_atom), RI_offset(bs_env%n_atom))
      CALL gw_3c_ctx_create(integral_context, bs_env, bs_env%ri_metric, &
                            basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                            basis_i=bs_env%basis_set_RI)
      CALL gw_3c_ws_create(workspace, integral_context)
      ilocal_cluster = 0
      DO icenter_atom = 1, bs_env%n_atom
         IF (MOD(icenter_atom - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         ilocal_cluster = ilocal_cluster + 1
         CALL get_rirs_cluster_atoms(particle_set, cell, icenter_atom, cluster_radius, &
                                     cluster_3c_int(ilocal_cluster)%atom_indices)
         AO_offset = 0
         RI_offset = 0
         nAO_cluster = 0
         nRI_cluster = 0
         DO icluster_atom = 1, SIZE(cluster_3c_int(ilocal_cluster)%atom_indices)
            iatom = cluster_3c_int(ilocal_cluster)%atom_indices(icluster_atom)
            AO_offset(iatom) = nAO_cluster
            RI_offset(iatom) = nRI_cluster
            nAO_cluster = nAO_cluster + bs_env%sizes_AO(iatom)
            nRI_cluster = nRI_cluster + sizes_ref_RI(iatom)
         END DO
         ALLOCATE (cluster_3c_int(ilocal_cluster)%Int_3c(nAO_cluster, nAO_cluster, nRI_cluster), &
                   source=0.0_dp)
         DO i_cluster_atom_p = 1, SIZE(cluster_3c_int(ilocal_cluster)%atom_indices)
            iatom_p = cluster_3c_int(ilocal_cluster)%atom_indices(i_cluster_atom_p)
            DO i_cluster_atom_k = 1, SIZE(cluster_3c_int(ilocal_cluster)%atom_indices)
               iatom_k = cluster_3c_int(ilocal_cluster)%atom_indices(i_cluster_atom_k)
               DO i_cluster_atom_j = 1, SIZE(cluster_3c_int(ilocal_cluster)%atom_indices)
                  iatom_j = cluster_3c_int(ilocal_cluster)%atom_indices(i_cluster_atom_j)
                  CALL build_3c_integral_block_ctx(cluster_3c_int(ilocal_cluster)%Int_3c, &
                                                   integral_context, workspace, &
                                                   atom_j=iatom_j, atom_k=iatom_k, atom_i=iatom_p, &
                                                   j_offset=AO_offset(iatom_j), &
                                                   k_offset=AO_offset(iatom_k), &
                                                   i_offset=RI_offset(iatom_p), screened=screened)
               END DO
            END DO
         END DO
         ! Zero auxiliary columns contribute neither to the objective nor its derivatives.
         nRI_nonzero = COUNT([(ANY(cluster_3c_int(ilocal_cluster)%Int_3c(:, :, ip) /= 0.0_dp), &
                               ip=1, nRI_cluster)])
         IF (nRI_nonzero < nRI_cluster) THEN
            ALLOCATE (Int_3c_nonzero(nAO_cluster, nAO_cluster, nRI_nonzero))
            nRI_nonzero = 0
            DO ip = 1, nRI_cluster
               IF (.NOT. ANY(cluster_3c_int(ilocal_cluster)%Int_3c(:, :, ip) /= 0.0_dp)) CYCLE
               nRI_nonzero = nRI_nonzero + 1
               Int_3c_nonzero(:, :, nRI_nonzero) = cluster_3c_int(ilocal_cluster)%Int_3c(:, :, ip)
            END DO
            CALL MOVE_ALLOC(Int_3c_nonzero, cluster_3c_int(ilocal_cluster)%Int_3c)
         END IF
      END DO
      CALL gw_3c_ws_release(workspace)
      CALL gw_3c_ctx_release(integral_context)
      CALL timestop(handle)
   END SUBROUTINE build_cluster_3c_int

! **************************************************************************************************
!> \brief Evaluate ϕ_μ(r_l) on G_A and compute the normalized local error E_A and its gradient.
!> \param cluster_3c_int Local atom indices and exact three-centre integrals.
!> \param bs_env ...
!> \param normalized_error Normalized squared three-centre-integral error.
!> \param coordinate_gradient Derivative with respect to the cluster's physical grid coordinates.
!> \param max_abs_error Largest absolute three-centre-integral error.
!> \param fit_successful Whether the regularized fitting equations were solved.
! **************************************************************************************************
   SUBROUTINE evaluate_local_cluster(cluster_3c_int, bs_env, normalized_error, coordinate_gradient, &
                                     max_abs_error, fit_successful)
      TYPE(local_cluster_3c_integrals_type), INTENT(IN)  :: cluster_3c_int
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: coordinate_gradient(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: max_abs_error
      LOGICAL, INTENT(OUT)                               :: fit_successful

      INTEGER                                            :: AO_offset, grid_point_offset, iatom, &
                                                            icluster_atom, ikind, n_grid_points, &
                                                            nAO, nAO_atom
      REAL(KIND=dp), ALLOCATABLE                         :: dPhi_alpha_l_mu(:, :, :), &
                                                            grid_points(:, :), Phi_l_mu(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: AO_basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      cell => bs_env%ri_rs%cell
      particle_set => bs_env%ri_rs%particle_set
      nAO = SIZE(cluster_3c_int%Int_3c, 1)
      n_grid_points = 0
      DO icluster_atom = 1, SIZE(cluster_3c_int%atom_indices)
         iatom = cluster_3c_int%atom_indices(icluster_atom)
         n_grid_points = n_grid_points + SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)
      END DO
      ALLOCATE (grid_points(3, n_grid_points), Phi_l_mu(n_grid_points, nAO), &
                dPhi_alpha_l_mu(3, n_grid_points, nAO), coordinate_gradient(3, n_grid_points))
      Phi_l_mu = 0.0_dp
      dPhi_alpha_l_mu = 0.0_dp
      grid_point_offset = 0
      DO icluster_atom = 1, SIZE(cluster_3c_int%atom_indices)
         iatom = cluster_3c_int%atom_indices(icluster_atom)
         grid_points(:, grid_point_offset + 1:grid_point_offset + &
                     SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)) = &
            SPREAD(particle_set(iatom)%r, 2, &
                   SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)) + &
            bs_env%ri_rs%atomic_grids(iatom)%raw_points
         grid_point_offset = grid_point_offset + SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)
      END DO
      AO_offset = 0
      DO icluster_atom = 1, SIZE(cluster_3c_int%atom_indices)
         iatom = cluster_3c_int%atom_indices(icluster_atom)
         ikind = particle_set(iatom)%atomic_kind%kind_number
         AO_basis => bs_env%basis_set_AO(ikind)%gto_basis_set
         nAO_atom = AO_basis%nsgf
         CALL evaluate_ao_basis_on_points(Phi_l_mu(:, AO_offset + 1:AO_offset + nAO_atom), &
                                          grid_points, AO_basis, particle_set(iatom)%r, cell, &
                                          dphi=dPhi_alpha_l_mu(:, :, AO_offset + 1:AO_offset + nAO_atom))
         AO_offset = AO_offset + nAO_atom
      END DO
      CALL evaluate_rirs_grid_cluster(Phi_l_mu, dPhi_alpha_l_mu, cluster_3c_int%Int_3c, &
                                      bs_env%ri_rs%tikhonov, normalized_error, coordinate_gradient, &
                                      max_abs_error, fit_successful)
   END SUBROUTINE evaluate_local_cluster

! **************************************************************************************************
!> \brief Solve the cluster-local Z_lP for one C_A and evaluate its normalized error and exact gradient.
!>
!> The solve applies the same column Jacobi scaling and Tikhonov parameter as the production
!> Z_lP construction. Z_prime_lP denotes the coefficients before undoing the Jacobi scaling.
!> \param Phi_l_mu AO values ϕ_μ(r_l), indexed (l, μ).
!> \param dPhi_alpha_l_mu Cartesian derivatives of ϕ_μ(r_l), indexed (α, l, μ).
!> \param Int_3c Exact three-centre integrals, indexed (mu,nu,P).
!> \param tikhonov Tikhonov parameter used by the production RI-RS solve.
!> \param normalized_error Normalized squared residual for this cluster.
!> \param coordinate_gradient Analytic derivative of normalized_error with respect to r(alpha,l).
!> \param max_abs_error Largest absolute error in an unweighted three-centre integral.
!> \param fit_successful False if the regularized normal equations cannot be solved.
! **************************************************************************************************
   SUBROUTINE evaluate_rirs_grid_cluster(Phi_l_mu, dPhi_alpha_l_mu, Int_3c, tikhonov, &
                                         normalized_error, coordinate_gradient, &
                                         max_abs_error, fit_successful)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Phi_l_mu
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: dPhi_alpha_l_mu, Int_3c
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: coordinate_gradient
      REAL(KIND=dp), INTENT(OUT)                         :: max_abs_error
      LOGICAL, INTENT(OUT)                               :: fit_successful

      CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_rirs_grid_cluster'
      REAL(KIND=dp), PARAMETER                           :: jacobi_floor = 1.0E-16_dp

      INTEGER                                            :: handle, i_ao_pair, icartesian, &
                                                            igrid_point, imu, inu, ip, &
                                                            lapack_info, n_AO_pairs, &
                                                            n_grid_points, nAO, nRI, phase_handle
      INTEGER, ALLOCATABLE                               :: mu_of_pair(:), nu_of_pair(:)
      REAL(KIND=dp)                                      :: absolute_error, dPhi_munu_l, &
                                                            Int_3c_norm2, jacobi_scale, &
                                                            Phi_munu_norm2, scaling_projection, &
                                                            symmetry_factor
      REAL(KIND=dp), ALLOCATABLE :: D_inverse_Z_prime_lP(:, :), &
         D_inverse_Z_prime_times_Z_prime_transpose_ll(:, :), D_ll(:, :), d_lP(:, :), &
         dE_dPhi_munu_l(:, :), dPhi_munu_l_scaled(:), Int_3c_munu_P(:, :), jacobi_scaling(:), &
         Phi_munu_l(:, :), Phi_munu_l_scaled(:, :), R_munu_P(:, :), Z_prime_lP(:, :)

      CALL timeset(routineN, handle)

      n_grid_points = SIZE(Phi_l_mu, 1)
      nAO = SIZE(Phi_l_mu, 2)
      nRI = SIZE(Int_3c, 3)
      n_AO_pairs = nAO*(nAO + 1)/2

      CPASSERT(SIZE(dPhi_alpha_l_mu, 1) == 3)
      CPASSERT(SIZE(dPhi_alpha_l_mu, 2) == n_grid_points)
      CPASSERT(SIZE(dPhi_alpha_l_mu, 3) == nAO)
      CPASSERT(SIZE(Int_3c, 1) == nAO)
      CPASSERT(SIZE(Int_3c, 2) == nAO)
      CPASSERT(SIZE(coordinate_gradient, 1) == 3)
      CPASSERT(SIZE(coordinate_gradient, 2) == n_grid_points)
      CPASSERT(tikhonov >= 0.0_dp)

      normalized_error = HUGE(normalized_error)
      coordinate_gradient = 0.0_dp
      max_abs_error = HUGE(max_abs_error)
      fit_successful = .FALSE.
      IF (n_grid_points < 1 .OR. nRI < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (Phi_munu_l(n_AO_pairs, n_grid_points), Phi_munu_l_scaled(n_AO_pairs, n_grid_points), &
                Int_3c_munu_P(n_AO_pairs, nRI), mu_of_pair(n_AO_pairs), nu_of_pair(n_AO_pairs), &
                jacobi_scaling(n_grid_points))
      CALL timeset("rirs_cluster_AO_products", phase_handle)
      i_ao_pair = 0
      DO inu = 1, nAO
         DO imu = 1, inu
            i_ao_pair = i_ao_pair + 1
            mu_of_pair(i_ao_pair) = imu
            nu_of_pair(i_ao_pair) = inu
            symmetry_factor = SQRT(REAL(2 - MERGE(1, 0, imu == inu), dp))
            DO igrid_point = 1, n_grid_points
               Phi_munu_l(i_ao_pair, igrid_point) = &
                  symmetry_factor*Phi_l_mu(igrid_point, imu)*Phi_l_mu(igrid_point, inu)
            END DO
            DO ip = 1, nRI
               Int_3c_munu_P(i_ao_pair, ip) = symmetry_factor*Int_3c(imu, inu, ip)
            END DO
         END DO
      END DO

      Int_3c_norm2 = SUM(Int_3c_munu_P*Int_3c_munu_P)
      CALL timestop(phase_handle)
      IF (Int_3c_norm2 <= TINY(1.0_dp)) THEN
         DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, mu_of_pair, nu_of_pair, jacobi_scaling)
         CALL timestop(handle)
         RETURN
      END IF

      ! D'_ll' = d_l D_ll' d_l' + λδ_ll', with d_l = 1/sqrt(D_ll).
      DO igrid_point = 1, n_grid_points
         Phi_munu_norm2 = SUM(Phi_munu_l(:, igrid_point)*Phi_munu_l(:, igrid_point))
         jacobi_scaling(igrid_point) = 1.0_dp/SQRT(MAX(Phi_munu_norm2, jacobi_floor))
         Phi_munu_l_scaled(:, igrid_point) = &
            jacobi_scaling(igrid_point)*Phi_munu_l(:, igrid_point)
      END DO
      ALLOCATE (D_ll(n_grid_points, n_grid_points), d_lP(n_grid_points, nRI), Z_prime_lP(n_grid_points, nRI))
      CALL timeset("rirs_cluster_D_ll", phase_handle)
      CALL dgemm('T', 'N', n_grid_points, n_grid_points, n_AO_pairs, 1.0_dp, Phi_munu_l_scaled, n_AO_pairs, &
                 Phi_munu_l_scaled, n_AO_pairs, 0.0_dp, D_ll, n_grid_points)
      CALL timestop(phase_handle)
      DO igrid_point = 1, n_grid_points
         D_ll(igrid_point, igrid_point) = D_ll(igrid_point, igrid_point) + tikhonov
      END DO
      CALL timeset("rirs_cluster_d_lP", phase_handle)
      CALL dgemm('T', 'N', n_grid_points, nRI, n_AO_pairs, 1.0_dp, Phi_munu_l_scaled, n_AO_pairs, &
                 Int_3c_munu_P, n_AO_pairs, 0.0_dp, d_lP, n_grid_points)
      CALL timestop(phase_handle)
      CALL dpotrf('L', n_grid_points, D_ll, n_grid_points, lapack_info)
      IF (lapack_info /= 0) THEN
         DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, D_ll, mu_of_pair, nu_of_pair, d_lP, &
                     jacobi_scaling, Z_prime_lP)
         CALL timestop(handle)
         RETURN
      END IF
      Z_prime_lP(:, :) = d_lP
      CALL dpotrs('L', n_grid_points, nRI, D_ll, n_grid_points, Z_prime_lP, n_grid_points, lapack_info)
      DEALLOCATE (d_lP)
      IF (lapack_info /= 0) THEN
         DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, D_ll, mu_of_pair, nu_of_pair, &
                     jacobi_scaling, Z_prime_lP)
         CALL timestop(handle)
         RETURN
      END IF

      CALL timeset("rirs_cluster_residual", phase_handle)
      ALLOCATE (R_munu_P(n_AO_pairs, nRI))
      R_munu_P(:, :) = Int_3c_munu_P
      ! R_μνP = (μν|P) - Σ_l ϕ_μ(r_l) ϕ_ν(r_l) d_l Z'_lP.
      CALL dgemm('N', 'N', n_AO_pairs, nRI, n_grid_points, -1.0_dp, Phi_munu_l_scaled, n_AO_pairs, &
                 Z_prime_lP, n_grid_points, 1.0_dp, R_munu_P, n_AO_pairs)
      normalized_error = SUM(R_munu_P*R_munu_P)/Int_3c_norm2
      CALL timestop(phase_handle)

      ! For λ>0, the residual derivative includes the response of the regularized coefficients:
      ! dE/dϕ = [-2 R Z^T - 2λ(R Y^T - ϕ Y Z^T)] / ||Int_3c||^2.
      ALLOCATE (D_inverse_Z_prime_lP(n_grid_points, nRI), dE_dPhi_munu_l(n_AO_pairs, n_grid_points))
      D_inverse_Z_prime_lP(:, :) = Z_prime_lP
      CALL dpotrs('L', n_grid_points, nRI, D_ll, n_grid_points, D_inverse_Z_prime_lP, n_grid_points, lapack_info)
      IF (lapack_info /= 0) THEN
         DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, dE_dPhi_munu_l, D_ll, mu_of_pair, nu_of_pair, &
                     R_munu_P, jacobi_scaling, D_inverse_Z_prime_lP, Z_prime_lP)
         CALL timestop(handle)
         RETURN
      END IF
      CALL timeset("rirs_cluster_gradient_products", phase_handle)
      CALL dgemm('N', 'T', n_AO_pairs, n_grid_points, nRI, -2.0_dp/Int_3c_norm2, R_munu_P, n_AO_pairs, &
                 Z_prime_lP, n_grid_points, 0.0_dp, dE_dPhi_munu_l, n_AO_pairs)
      IF (tikhonov > 0.0_dp) THEN
         CALL dgemm('N', 'T', n_AO_pairs, n_grid_points, nRI, -2.0_dp*tikhonov/Int_3c_norm2, R_munu_P, n_AO_pairs, &
                    D_inverse_Z_prime_lP, n_grid_points, 1.0_dp, dE_dPhi_munu_l, n_AO_pairs)
         ALLOCATE (D_inverse_Z_prime_times_Z_prime_transpose_ll(n_grid_points, n_grid_points))
         CALL dgemm('N', 'T', n_grid_points, n_grid_points, nRI, 1.0_dp, D_inverse_Z_prime_lP, n_grid_points, &
                    Z_prime_lP, n_grid_points, 0.0_dp, D_inverse_Z_prime_times_Z_prime_transpose_ll, n_grid_points)
         CALL dgemm('N', 'N', n_AO_pairs, n_grid_points, n_grid_points, &
                    2.0_dp*tikhonov/Int_3c_norm2, Phi_munu_l_scaled, n_AO_pairs, &
                    D_inverse_Z_prime_times_Z_prime_transpose_ll, n_grid_points, 1.0_dp, dE_dPhi_munu_l, n_AO_pairs)
         DEALLOCATE (D_inverse_Z_prime_times_Z_prime_transpose_ll)
      END IF

      CALL timestop(phase_handle)
      CALL timeset("rirs_cluster_gradient_coordinates", phase_handle)
      ALLOCATE (dPhi_munu_l_scaled(n_AO_pairs))
      DO igrid_point = 1, n_grid_points
         jacobi_scale = jacobi_scaling(igrid_point)
         Phi_munu_norm2 = SUM(Phi_munu_l(:, igrid_point)*Phi_munu_l(:, igrid_point))
         DO icartesian = 1, 3
            DO i_ao_pair = 1, n_AO_pairs
               imu = mu_of_pair(i_ao_pair)
               inu = nu_of_pair(i_ao_pair)
               symmetry_factor = SQRT(REAL(2 - MERGE(1, 0, imu == inu), dp))
               dPhi_munu_l = symmetry_factor*( &
                             dPhi_alpha_l_mu(icartesian, igrid_point, imu)* &
                             Phi_l_mu(igrid_point, inu) + &
                             Phi_l_mu(igrid_point, imu)* &
                             dPhi_alpha_l_mu(icartesian, igrid_point, inu))
               dPhi_munu_l_scaled(i_ao_pair) = jacobi_scale*dPhi_munu_l
            END DO
            IF (Phi_munu_norm2 > jacobi_floor) THEN
               scaling_projection = &
                  DOT_PRODUCT(Phi_munu_l(:, igrid_point), dPhi_munu_l_scaled)/jacobi_scale
               dPhi_munu_l_scaled(:) = dPhi_munu_l_scaled - &
                                       jacobi_scale**3*Phi_munu_l(:, igrid_point)*scaling_projection
            END IF
            coordinate_gradient(icartesian, igrid_point) = &
               DOT_PRODUCT(dE_dPhi_munu_l(:, igrid_point), dPhi_munu_l_scaled)
         END DO
      END DO

      CALL timestop(phase_handle)
      max_abs_error = 0.0_dp
      DO ip = 1, nRI
         DO i_ao_pair = 1, n_AO_pairs
            symmetry_factor = &
               SQRT(REAL(2 - MERGE(1, 0, mu_of_pair(i_ao_pair) == nu_of_pair(i_ao_pair)), dp))
            absolute_error = ABS(R_munu_P(i_ao_pair, ip))/symmetry_factor
            max_abs_error = MAX(max_abs_error, absolute_error)
         END DO
      END DO

      fit_successful = .TRUE.
      DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, dE_dPhi_munu_l, dPhi_munu_l_scaled, D_ll, &
                  mu_of_pair, nu_of_pair, R_munu_P, jacobi_scaling, D_inverse_Z_prime_lP, Z_prime_lP)
      CALL timestop(handle)
   END SUBROUTINE evaluate_rirs_grid_cluster

! **************************************************************************************************
!> \brief Sum gradient contributions from overlapping C_A into each atom-centred grid coordinate.
!> \param cluster_atoms Atom indices in the local cluster.
!> \param grids Atom-centred RI-RS grids.
!> \param atom_coordinate_offsets Starting coordinate offset for each atom.
!> \param cluster_coordinate_gradient Gradient for the cluster's physical grid points.
!> \param coordinate_gradient Global flattened atom-relative gradient to update.
! **************************************************************************************************
   SUBROUTINE accumulate_atom_gradient(cluster_atoms, grids, atom_coordinate_offsets, &
                                       cluster_coordinate_gradient, coordinate_gradient)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: cluster_atoms
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(IN)     :: grids
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_coordinate_offsets
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cluster_coordinate_gradient
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: coordinate_gradient

      INTEGER                                            :: grid_point_offset, iatom, icartesian, &
                                                            icluster_atom, igrid_point

      grid_point_offset = 0
      DO icluster_atom = 1, SIZE(cluster_atoms)
         iatom = cluster_atoms(icluster_atom)
         DO igrid_point = 1, SIZE(grids(iatom)%raw_points, 2)
            DO icartesian = 1, 3
               coordinate_gradient(atom_coordinate_offsets(iatom) + 3*(igrid_point - 1) + icartesian) = &
                  coordinate_gradient(atom_coordinate_offsets(iatom) + 3*(igrid_point - 1) + icartesian) + &
                  cluster_coordinate_gradient(icartesian, grid_point_offset + igrid_point)
            END DO
         END DO
         grid_point_offset = grid_point_offset + SIZE(grids(iatom)%raw_points, 2)
      END DO
   END SUBROUTINE accumulate_atom_gradient

! **************************************************************************************************
!> \brief Pack all atom-centred grids into the single Cartesian vector required by L-BFGS.
!> \param grids Atom-centred RI-RS grids.
!> \param grid_coordinates Flattened atom-relative grid coordinates.
!> \param atom_coordinate_offsets Starting coordinate offset for each atom.
!>
!> SIZE(grid_coordinates) = 3 Σ_A N_grid,A. Coordinates are ordered Cartesian component first,
!> then grid point, then atom, matching the column-major layout of raw_points(3,N_grid,A).
! **************************************************************************************************
   SUBROUTINE pack_lbfgs_grid_coordinates(grids, grid_coordinates, atom_coordinate_offsets)
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(IN)     :: grids
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: grid_coordinates(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: atom_coordinate_offsets(:)

      INTEGER                                            :: atom_coordinate_count, iatom, &
                                                            n_grid_coordinates

      ALLOCATE (atom_coordinate_offsets(SIZE(grids)))
      n_grid_coordinates = 0
      DO iatom = 1, SIZE(grids)
         atom_coordinate_offsets(iatom) = n_grid_coordinates
         n_grid_coordinates = n_grid_coordinates + SIZE(grids(iatom)%raw_points)
      END DO
      ALLOCATE (grid_coordinates(n_grid_coordinates))

      DO iatom = 1, SIZE(grids)
         atom_coordinate_count = SIZE(grids(iatom)%raw_points)
         grid_coordinates(atom_coordinate_offsets(iatom) + 1: &
                          atom_coordinate_offsets(iatom) + atom_coordinate_count) = &
            RESHAPE(grids(iatom)%raw_points, [atom_coordinate_count])
      END DO
   END SUBROUTINE pack_lbfgs_grid_coordinates

! **************************************************************************************************
!> \brief Restore the optimized Cartesian vector to the persistent atom-centred RI-RS grids.
!> \param grid_coordinates Flattened atom-relative grid coordinates.
!> \param atom_coordinate_offsets Starting coordinate offset for each atom.
!> \param grids Atom-centred RI-RS grids to update.
! **************************************************************************************************
   SUBROUTINE unpack_lbfgs_grid_coordinates(grid_coordinates, atom_coordinate_offsets, grids)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: grid_coordinates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_coordinate_offsets
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(INOUT)  :: grids

      INTEGER                                            :: atom_coordinate_count, iatom

      CPASSERT(SIZE(atom_coordinate_offsets) == SIZE(grids))
      DO iatom = 1, SIZE(grids)
         atom_coordinate_count = SIZE(grids(iatom)%raw_points)
         grids(iatom)%raw_points(:, :) = RESHAPE( &
                                         grid_coordinates(atom_coordinate_offsets(iatom) + 1: &
                                                          atom_coordinate_offsets(iatom) + &
                                                          atom_coordinate_count), &
                                         SHAPE(grids(iatom)%raw_points))
      END DO
   END SUBROUTINE unpack_lbfgs_grid_coordinates

END MODULE gw_ri_rs_grid_optimization
