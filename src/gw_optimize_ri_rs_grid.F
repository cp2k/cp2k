!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local-environment optimizer for atom-centred real-space RI grids.
!> \par History
!>      09.2026 created [Jan Wilhelm]
! **************************************************************************************************
MODULE gw_optimize_ri_rs_grid
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_lbfgs,                        ONLY: setulb
   USE gw_integrals,                    ONLY: build_3c_integral_block_ctx,&
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
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE rirs_grid_utils,                 ONLY: evaluate_ao_basis_on_points
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_optimize_ri_rs_grid'

   TYPE :: packed_atom_grids_type
      INTEGER, ALLOCATABLE                 :: atom_coordinate_offsets(:)
      REAL(KIND=dp), ALLOCATABLE           :: coordinates(:), &
                                              lower_bounds(:), &
                                              upper_bounds(:)
   END TYPE packed_atom_grids_type

   TYPE :: local_cluster_type
      INTEGER, ALLOCATABLE                 :: atom_indices(:)
      REAL(KIND=dp), ALLOCATABLE           :: Int_3c(:, :, :)
   END TYPE local_cluster_type

   PUBLIC :: optimize_ri_rs_grid

CONTAINS

! **************************************************************************************************
!> \brief Coordinate initialization, local training-data construction, optimization, and final storage.
!> \param qs_env Quickstep environment used to build the local-cluster three-centre integrals.
!> \param bs_env GW environment containing all configuration, basis, geometry, parallel, and grid data.
! **************************************************************************************************
   SUBROUTINE optimize_ri_rs_grid(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'optimize_ri_rs_grid'

      INTEGER                                            :: handle
      LOGICAL                                            :: fit_successful
      REAL(KIND=dp)                                      :: max_abs_error, normalized_error, &
                                                            start_time
      TYPE(local_cluster_type), ALLOCATABLE              :: clusters(:)
      TYPE(packed_atom_grids_type)                       :: packed_grids

      CALL timeset(routineN, handle)

      start_time = m_walltime()

      ! Validate the nonperiodic RI-RS setup and report the start of the optimization.
      CALL check_RI_RS_grid_generation(bs_env)

      ! Pack x = {r_(αl)} needed for L-BFGS optimizer, and set numerical L-BFGS bounds.
      CALL pack_atom_grids(bs_env%ri_rs%grid_cache, packed_grids)

      ! C_A = {B: |R_A - R_B| < R_cut}; cache (μν|P) for μ, ν, and P centered in C_A.
      CALL build_local_clusters(qs_env, bs_env, clusters)

      ! E_A = Σ_{μνP∈C_A} [(μν|P) - Σ_{l∈G_A} φ_μ(r_l) φ_ν(r_l) Z_lP^(A)]².
      ! Minimize E_loc(norm) = N_atom^(-1) Σ_A [E_A / Σ_{μνP∈C_A} |(μν|P)|²].
      CALL optimize_grid_coordinates(packed_grids, bs_env, clusters, normalized_error, &
                                     max_abs_error, fit_successful)

      ! Store the optimized grids, invalidate the old Z_lP, and report E_loc(norm).
      CALL finalize_RI_RS_grid_optimization(bs_env, normalized_error, max_abs_error, &
                                            fit_successful, start_time)

      CALL timestop(handle)

   END SUBROUTINE optimize_ri_rs_grid

! **************************************************************************************************
!> \brief Validate the RI-RS grid-optimization environment and announce the calculation.
!> \param bs_env GW environment containing the optimizer input and initialized grid data.
! **************************************************************************************************
   SUBROUTINE check_RI_RS_grid_generation(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: atom, atom_kind, unit_nr

      CPASSERT(ASSOCIATED(bs_env%ri_rs%cell))
      CPASSERT(ASSOCIATED(bs_env%ri_rs%particle_set))
      CPASSERT(ASSOCIATED(bs_env%para_env))
      CPASSERT(ALLOCATED(bs_env%basis_set_AO))
      CPASSERT(ALLOCATED(bs_env%basis_set_RI))
      CPASSERT(ALLOCATED(bs_env%sizes_AO))
      CPASSERT(ALLOCATED(bs_env%ri_rs%grid_cache))
      CPASSERT(SIZE(bs_env%sizes_AO) == bs_env%n_atom)
      CPASSERT(SIZE(bs_env%ri_rs%grid_cache) == bs_env%n_atom)

      IF (bs_env%ri_rs%grid_opt%max_iter < 1) THEN
         CPABORT("GRID_OPTIMIZATION%MAX_ITER must be positive")
      END IF
      IF (bs_env%ri_rs%grid_opt%cutoff_atomic_cluster <= 0.0_dp) THEN
         CPABORT("GRID_OPTIMIZATION%CUTOFF_ATOMIC_CLUSTER must be positive")
      END IF
      IF (bs_env%ri_rs%tikhonov < 0.0_dp) THEN
         CPABORT("RI_RS%TIKHONOV must not be negative")
      END IF
      IF (bs_env%do_periodic) THEN
         CPABORT("GRID_OPTIMIZATION currently supports nonperiodic local environments only")
      END IF

      DO atom = 1, bs_env%n_atom
         CPASSERT(ALLOCATED(bs_env%ri_rs%grid_cache(atom)%raw_points))
         CPASSERT(SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 1) == 3)
         CPASSERT(SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 2) > 0)
         atom_kind = bs_env%ri_rs%particle_set(atom)%atomic_kind%kind_number
         CPASSERT(atom_kind >= 1 .AND. atom_kind <= SIZE(bs_env%basis_set_AO))
         CPASSERT(atom_kind <= SIZE(bs_env%basis_set_RI))
         CPASSERT(ASSOCIATED(bs_env%basis_set_AO(atom_kind)%gto_basis_set))
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(atom_kind)%gto_basis_set))
         IF (bs_env%sizes_AO(atom) < 1 .OR. &
             bs_env%basis_set_RI(atom_kind)%gto_basis_set%nsgf < 1) THEN
            CPABORT("Every atom in GRID_OPTIMIZATION needs ORB and RI_AUX functions")
         END IF
      END DO

      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'Started RI-RS grid optimization'
         CALL m_flush(unit_nr)
      END IF
   END SUBROUTINE check_RI_RS_grid_generation

! **************************************************************************************************
!> \brief Store and print the result of a successful RI-RS grid optimization.
!> \param bs_env GW environment containing the optimized grids.
!> \param normalized_error Mean normalized squared three-center-integral error.
!> \param max_abs_error Largest absolute three-center-integral error.
!> \param fit_successful Whether the final regularized fit was successful.
!> \param start_time Wall-clock time at which the optimization started.
! **************************************************************************************************
   SUBROUTINE finalize_RI_RS_grid_optimization(bs_env, normalized_error, max_abs_error, &
                                               fit_successful, start_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: normalized_error, max_abs_error
      LOGICAL, INTENT(IN)                                :: fit_successful
      REAL(KIND=dp), INTENT(IN)                          :: start_time

      INTEGER                                            :: atom, unit_nr

      IF (.NOT. fit_successful) CPABORT("RI-RS grid optimization produced no regularized fit")

      DO atom = 1, SIZE(bs_env%ri_rs%grid_cache)
         bs_env%ri_rs%grid_cache(atom)%npts = &
            SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 2)
      END DO
      bs_env%ri_rs%Z_lP_exists = .FALSE.

      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T69,F10.1,A)') &
            'RI-RS grid optimization completed, execution time:', m_walltime() - start_time, ' s'
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Normalized 3C error:', normalized_error
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Maximum absolute 3C error:', max_abs_error
         FLUSH (unit_nr)
      END IF
   END SUBROUTINE finalize_RI_RS_grid_optimization

! **************************************************************************************************
!> \brief Evaluate the normalized local objective E_loc and its gradient for one L-BFGS trial vector.
!> \param packed_grids Flattened coordinates and their atom-to-vector mapping.
!> \param bs_env GW environment containing RI-RS configuration and grid data.
!> \param clusters Rank-local clusters and their exact three-centre integrals.
!> \param normalized_error Mean normalized squared three-centre-integral error.
!> \param coordinate_gradient Derivative of normalized_error with respect to coordinates.
!> \param max_abs_error Largest absolute three-centre-integral error.
!> \param fit_successful Whether every local-cluster evaluation succeeded.
! **************************************************************************************************
   SUBROUTINE grid_objective(packed_grids, bs_env, clusters, normalized_error, coordinate_gradient, &
                             max_abs_error, fit_successful)
      TYPE(packed_atom_grids_type), INTENT(IN)           :: packed_grids
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(local_cluster_type), DIMENSION(:), INTENT(IN) :: clusters
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coordinate_gradient
      REAL(KIND=dp), INTENT(OUT)                         :: max_abs_error
      LOGICAL, INTENT(OUT)                               :: fit_successful

      INTEGER                                            :: cluster_index, evaluated_cluster_count, &
                                                            successful_ranks
      LOGICAL                                            :: cluster_fit_successful
      REAL(KIND=dp)                                      :: cluster_max_abs_error, &
                                                            cluster_normalized_error
      REAL(KIND=dp), ALLOCATABLE :: cluster_coordinate_gradient(:, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      para_env => bs_env%para_env
      CPASSERT(SIZE(coordinate_gradient) == SIZE(packed_grids%coordinates))
      CALL unpack_atom_grids(packed_grids, bs_env%ri_rs%grid_cache)
      normalized_error = 0.0_dp
      coordinate_gradient = 0.0_dp
      evaluated_cluster_count = 0
      max_abs_error = 0.0_dp
      fit_successful = .TRUE.
      DO cluster_index = 1, SIZE(clusters)
         CALL evaluate_local_cluster(clusters(cluster_index), bs_env, &
                                     cluster_normalized_error, cluster_coordinate_gradient, &
                                     cluster_max_abs_error, cluster_fit_successful)
         IF (.NOT. cluster_fit_successful) THEN
            fit_successful = .FALSE.
            IF (ALLOCATED(cluster_coordinate_gradient)) DEALLOCATE (cluster_coordinate_gradient)
            EXIT
         END IF
         normalized_error = normalized_error + cluster_normalized_error
         CALL accumulate_atom_gradient(clusters(cluster_index)%atom_indices, bs_env%ri_rs%grid_cache, &
                                       packed_grids%atom_coordinate_offsets, cluster_coordinate_gradient, &
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
!> \param packed_grids Coordinates on entry and optimized coordinates on return.
!> \param bs_env GW environment containing RI-RS configuration and grid data.
!> \param clusters Rank-local clusters and their exact three-centre integrals.
!> \param normalized_error Final mean normalized squared three-centre-integral error.
!> \param max_abs_error Final largest absolute three-centre-integral error.
!> \param fit_successful Whether the final regularized fitting equations were solved.
! **************************************************************************************************
   SUBROUTINE optimize_grid_coordinates(packed_grids, bs_env, clusters, normalized_error, &
                                        max_abs_error, fit_successful)
      TYPE(packed_atom_grids_type), INTENT(INOUT)        :: packed_grids
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(local_cluster_type), DIMENSION(:), INTENT(IN) :: clusters
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error, max_abs_error
      LOGICAL, INTENT(OUT)                               :: fit_successful

      INTEGER, PARAMETER                                 :: history_size = 7
      REAL(KIND=dp), PARAMETER                           :: optimizer_accuracy = 1.0E-9_dp

      CHARACTER(LEN=60)                                  :: line_search_state, optimizer_task
      INTEGER                                            :: evaluation_count
      REAL(KIND=dp), &
         DIMENSION(SIZE(packed_grids%coordinates))       :: best_grid_coordinates, &
                                                            coordinate_gradient
      REAL(KIND=dp), DIMENSION(29)                       :: real_state
      REAL(KIND=dp), DIMENSION(2*history_size*SIZE(&
         packed_grids%coordinates)+5*SIZE(packed_grids%&
         coordinates)+11*history_size**2+8*history_size) :: workspace
      REAL(KIND=dp)                                      :: best_normalized_error
      LOGICAL, DIMENSION(4)                              :: logical_state
      LOGICAL                                            :: evaluation_successful, have_best
      INTEGER, DIMENSION(SIZE(packed_grids%coordinates)) :: bound_types
      INTEGER, DIMENSION(44)                             :: integer_state
      INTEGER, &
         DIMENSION(3*SIZE(packed_grids%coordinates))     :: integer_workspace

      CPASSERT(SIZE(packed_grids%coordinates) > 0)
      CPASSERT(ALL(SHAPE(packed_grids%lower_bounds) == SHAPE(packed_grids%coordinates)))
      CPASSERT(ALL(SHAPE(packed_grids%upper_bounds) == SHAPE(packed_grids%coordinates)))
      CPASSERT(ALL(packed_grids%lower_bounds <= packed_grids%upper_bounds))

      bound_types = 2
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
      best_grid_coordinates = packed_grids%coordinates

      DO
         CALL setulb(SIZE(packed_grids%coordinates), history_size, packed_grids%coordinates, &
                     packed_grids%lower_bounds, packed_grids%upper_bounds, bound_types, &
                     normalized_error, coordinate_gradient, 0.0_dp, optimizer_accuracy, &
                     workspace, integer_workspace, optimizer_task, -1, line_search_state, &
                     logical_state, integer_state, real_state, -1.0_dp)
         IF (optimizer_task(1:2) == 'FG') THEN
            IF (evaluation_count >= bs_env%ri_rs%grid_opt%max_iter) EXIT
            CALL grid_objective(packed_grids, bs_env, clusters, normalized_error, coordinate_gradient, &
                                max_abs_error, evaluation_successful)
            evaluation_count = evaluation_count + 1
            IF (.NOT. evaluation_successful) EXIT
            IF (.NOT. have_best .OR. normalized_error < best_normalized_error) THEN
               have_best = .TRUE.
               best_normalized_error = normalized_error
               best_grid_coordinates = packed_grids%coordinates
            END IF
         ELSE IF (optimizer_task(1:5) == 'NEW_X') THEN
            CYCLE
         ELSE
            EXIT
         END IF
      END DO
      IF (have_best) packed_grids%coordinates(:) = best_grid_coordinates
      CALL grid_objective(packed_grids, bs_env, clusters, normalized_error, coordinate_gradient, &
                          max_abs_error, fit_successful)
   END SUBROUTINE optimize_grid_coordinates

! **************************************************************************************************
!> \brief Build each C_A and cache its exact (μν|P); this routine owns the integral context.
!> \param qs_env Quickstep environment required to create the three-centre-integral context.
!> \param bs_env GW environment containing basis, geometry, parallel, and grid data.
!> \param clusters Rank-local clusters and their exact source-basis three-centre integrals.
! **************************************************************************************************
   SUBROUTINE build_local_clusters(qs_env, bs_env, clusters)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(local_cluster_type), ALLOCATABLE, INTENT(OUT) :: clusters(:)

      INTEGER :: AO_atom_j, AO_atom_j_index, AO_atom_k, AO_atom_k_index, atom, atom_kind, &
         center_atom, cluster_atom, cluster_atom_count, local_cluster, n_local_clusters, &
         nAO_cluster, nRI_cluster, nRI_nonzero, P, P_atom, P_atom_index
      INTEGER, ALLOCATABLE                               :: AO_offset(:), atom_buffer(:), &
                                                            RI_offset(:), sizes_ref_RI(:)
      LOGICAL                                            :: screened
      REAL(KIND=dp)                                      :: cluster_radius
      REAL(KIND=dp), ALLOCATABLE                         :: Int_3c_nonzero(:, :, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gw_3c_ctx_type)                               :: integral_context
      TYPE(gw_3c_ws_type)                                :: workspace
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      cell => bs_env%ri_rs%cell
      para_env => bs_env%para_env
      particle_set => bs_env%ri_rs%particle_set
      cluster_radius = bs_env%ri_rs%grid_opt%cutoff_atomic_cluster
      ALLOCATE (sizes_ref_RI(bs_env%n_atom))
      DO atom = 1, bs_env%n_atom
         atom_kind = particle_set(atom)%atomic_kind%kind_number
         sizes_ref_RI(atom) = bs_env%basis_set_RI(atom_kind)%gto_basis_set%nsgf
      END DO
      n_local_clusters = COUNT([(MOD(center_atom - 1, para_env%num_pe) == para_env%mepos, &
                                 center_atom=1, bs_env%n_atom)])
      ALLOCATE (clusters(n_local_clusters), atom_buffer(bs_env%n_atom), &
                AO_offset(bs_env%n_atom), RI_offset(bs_env%n_atom))
      CALL gw_3c_ctx_create(integral_context, qs_env, bs_env%ri_metric, &
                            bs_env%basis_set_AO, bs_env%basis_set_AO, bs_env%basis_set_RI)
      CALL gw_3c_ws_create(workspace, integral_context)
      local_cluster = 0
      DO center_atom = 1, bs_env%n_atom
         IF (MOD(center_atom - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         local_cluster = local_cluster + 1
         cluster_atom_count = 0
         ! C_A = {B : |R_A - R_B| < R_cut}; every C_A includes its centre atom A.
         DO atom = 1, bs_env%n_atom
            IF (SUM(pbc(particle_set(atom)%r - particle_set(center_atom)%r, cell)**2) <= &
                cluster_radius**2) THEN
               cluster_atom_count = cluster_atom_count + 1
               atom_buffer(cluster_atom_count) = atom
            END IF
         END DO
         ALLOCATE (clusters(local_cluster)%atom_indices(cluster_atom_count), &
                   source=atom_buffer(1:cluster_atom_count))
         AO_offset = 0
         RI_offset = 0
         nAO_cluster = 0
         nRI_cluster = 0
         DO cluster_atom = 1, SIZE(clusters(local_cluster)%atom_indices)
            atom = clusters(local_cluster)%atom_indices(cluster_atom)
            AO_offset(atom) = nAO_cluster
            RI_offset(atom) = nRI_cluster
            nAO_cluster = nAO_cluster + bs_env%sizes_AO(atom)
            nRI_cluster = nRI_cluster + sizes_ref_RI(atom)
         END DO
         ALLOCATE (clusters(local_cluster)%Int_3c(nAO_cluster, nAO_cluster, nRI_cluster), &
                   source=0.0_dp)
         DO P_atom_index = 1, SIZE(clusters(local_cluster)%atom_indices)
            P_atom = clusters(local_cluster)%atom_indices(P_atom_index)
            DO AO_atom_k_index = 1, SIZE(clusters(local_cluster)%atom_indices)
               AO_atom_k = clusters(local_cluster)%atom_indices(AO_atom_k_index)
               DO AO_atom_j_index = 1, SIZE(clusters(local_cluster)%atom_indices)
                  AO_atom_j = clusters(local_cluster)%atom_indices(AO_atom_j_index)
                  CALL build_3c_integral_block_ctx(clusters(local_cluster)%Int_3c, &
                                                   integral_context, workspace, &
                                                   atom_j=AO_atom_j, atom_k=AO_atom_k, atom_i=P_atom, &
                                                   j_offset=AO_offset(AO_atom_j), &
                                                   k_offset=AO_offset(AO_atom_k), &
                                                   i_offset=RI_offset(P_atom), screened=screened)
               END DO
            END DO
         END DO
         ! Zero auxiliary columns contribute neither to the objective nor its derivatives.
         nRI_nonzero = COUNT([(ANY(clusters(local_cluster)%Int_3c(:, :, P) /= 0.0_dp), &
                               P=1, nRI_cluster)])
         IF (nRI_nonzero < nRI_cluster) THEN
            ALLOCATE (Int_3c_nonzero(nAO_cluster, nAO_cluster, nRI_nonzero))
            nRI_nonzero = 0
            DO P = 1, nRI_cluster
               IF (.NOT. ANY(clusters(local_cluster)%Int_3c(:, :, P) /= 0.0_dp)) CYCLE
               nRI_nonzero = nRI_nonzero + 1
               Int_3c_nonzero(:, :, nRI_nonzero) = clusters(local_cluster)%Int_3c(:, :, P)
            END DO
            CALL MOVE_ALLOC(Int_3c_nonzero, clusters(local_cluster)%Int_3c)
         END IF
      END DO
      CALL gw_3c_ws_release(workspace)
      CALL gw_3c_ctx_release(integral_context)
   END SUBROUTINE build_local_clusters

! **************************************************************************************************
!> \brief Build Phi_l_mu on G_A and evaluate the normalized local error E_A and its gradient.
!> \param cluster Local atoms and exact three-centre integrals.
!> \param bs_env GW environment containing orbital basis sets and RI-RS parameters.
!> \param normalized_error Normalized squared three-centre-integral error.
!> \param coordinate_gradient Derivative with respect to the cluster's physical grid coordinates.
!> \param max_abs_error Largest absolute three-centre-integral error.
!> \param fit_successful Whether the regularized fitting equations were solved.
! **************************************************************************************************
   SUBROUTINE evaluate_local_cluster(cluster, bs_env, normalized_error, coordinate_gradient, &
                                     max_abs_error, fit_successful)
      TYPE(local_cluster_type), INTENT(IN)               :: cluster
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: coordinate_gradient(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: max_abs_error
      LOGICAL, INTENT(OUT)                               :: fit_successful

      INTEGER                                            :: AO_offset, atom, atom_kind, &
                                                            cluster_atom, grid_point_offset, &
                                                            n_grid_points, nAO, nAO_atom
      REAL(KIND=dp), ALLOCATABLE                         :: dPhi_alpha_l_mu(:, :, :), &
                                                            grid_points(:, :), Phi_l_mu(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: AO_basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      cell => bs_env%ri_rs%cell
      particle_set => bs_env%ri_rs%particle_set
      nAO = SIZE(cluster%Int_3c, 1)
      n_grid_points = 0
      DO cluster_atom = 1, SIZE(cluster%atom_indices)
         atom = cluster%atom_indices(cluster_atom)
         n_grid_points = n_grid_points + SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 2)
      END DO
      ALLOCATE (grid_points(3, n_grid_points), Phi_l_mu(n_grid_points, nAO), &
                dPhi_alpha_l_mu(3, n_grid_points, nAO), coordinate_gradient(3, n_grid_points))
      Phi_l_mu = 0.0_dp
      dPhi_alpha_l_mu = 0.0_dp
      grid_point_offset = 0
      DO cluster_atom = 1, SIZE(cluster%atom_indices)
         atom = cluster%atom_indices(cluster_atom)
         grid_points(:, grid_point_offset + 1:grid_point_offset + &
                     SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 2)) = &
            SPREAD(particle_set(atom)%r, 2, &
                   SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 2)) + &
            bs_env%ri_rs%grid_cache(atom)%raw_points
         grid_point_offset = grid_point_offset + SIZE(bs_env%ri_rs%grid_cache(atom)%raw_points, 2)
      END DO
      AO_offset = 0
      DO cluster_atom = 1, SIZE(cluster%atom_indices)
         atom = cluster%atom_indices(cluster_atom)
         atom_kind = particle_set(atom)%atomic_kind%kind_number
         AO_basis => bs_env%basis_set_AO(atom_kind)%gto_basis_set
         nAO_atom = AO_basis%nsgf
         CALL evaluate_ao_basis_on_points(Phi_l_mu(:, AO_offset + 1:AO_offset + nAO_atom), &
                                          grid_points, AO_basis, particle_set(atom)%r, cell, &
                                          dphi=dPhi_alpha_l_mu(:, :, AO_offset + 1:AO_offset + nAO_atom))
         AO_offset = AO_offset + nAO_atom
      END DO
      CALL evaluate_rirs_grid_cluster(Phi_l_mu, dPhi_alpha_l_mu, cluster%Int_3c, &
                                      bs_env%ri_rs%tikhonov, normalized_error, coordinate_gradient, &
                                      max_abs_error, fit_successful)
   END SUBROUTINE evaluate_local_cluster

! **************************************************************************************************
!> \brief Solve the cluster-local Z_lP for one C_A and evaluate its normalized error and exact gradient.
!>
!> The solve applies the same column Jacobi scaling and Tikhonov parameter as the production
!> Z_lP construction. Z_prime_lP denotes the coefficients before undoing the Jacobi scaling.
!> \param Phi_l_mu AO values, indexed (l,mu).
!> \param dPhi_alpha_l_mu Cartesian derivatives of AO values, indexed (alpha,l,mu).
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

      INTEGER                                            :: AO_pair, cartesian_component, &
                                                            grid_point, handle, lapack_info, mu, &
                                                            n_AO_pairs, n_grid_points, nAO, nRI, &
                                                            nu, P
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
      AO_pair = 0
      DO nu = 1, nAO
         DO mu = 1, nu
            AO_pair = AO_pair + 1
            mu_of_pair(AO_pair) = mu
            nu_of_pair(AO_pair) = nu
            symmetry_factor = SQRT(REAL(2 - MERGE(1, 0, mu == nu), dp))
            DO grid_point = 1, n_grid_points
               Phi_munu_l(AO_pair, grid_point) = &
                  symmetry_factor*Phi_l_mu(grid_point, mu)*Phi_l_mu(grid_point, nu)
            END DO
            DO P = 1, nRI
               Int_3c_munu_P(AO_pair, P) = symmetry_factor*Int_3c(mu, nu, P)
            END DO
         END DO
      END DO

      Int_3c_norm2 = SUM(Int_3c_munu_P*Int_3c_munu_P)
      IF (Int_3c_norm2 <= TINY(1.0_dp)) THEN
         DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, mu_of_pair, nu_of_pair, jacobi_scaling)
         CALL timestop(handle)
         RETURN
      END IF

      ! D_ll' = d_l(Φ^T Φ)_ll'd_l' + λδ_ll', with d_l = 1/sqrt((Φ^T Φ)_ll).
      DO grid_point = 1, n_grid_points
         Phi_munu_norm2 = SUM(Phi_munu_l(:, grid_point)*Phi_munu_l(:, grid_point))
         jacobi_scaling(grid_point) = 1.0_dp/SQRT(MAX(Phi_munu_norm2, jacobi_floor))
         Phi_munu_l_scaled(:, grid_point) = jacobi_scaling(grid_point)*Phi_munu_l(:, grid_point)
      END DO
      ALLOCATE (D_ll(n_grid_points, n_grid_points), d_lP(n_grid_points, nRI), Z_prime_lP(n_grid_points, nRI))
      CALL dgemm('T', 'N', n_grid_points, n_grid_points, n_AO_pairs, 1.0_dp, Phi_munu_l_scaled, n_AO_pairs, &
                 Phi_munu_l_scaled, n_AO_pairs, 0.0_dp, D_ll, n_grid_points)
      DO grid_point = 1, n_grid_points
         D_ll(grid_point, grid_point) = D_ll(grid_point, grid_point) + tikhonov
      END DO
      CALL dgemm('T', 'N', n_grid_points, nRI, n_AO_pairs, 1.0_dp, Phi_munu_l_scaled, n_AO_pairs, &
                 Int_3c_munu_P, n_AO_pairs, 0.0_dp, d_lP, n_grid_points)
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

      ALLOCATE (R_munu_P(n_AO_pairs, nRI))
      R_munu_P(:, :) = Int_3c_munu_P
      ! R_μνP = (μν|P) - Σ_l Φ_μν,l d_l Z'_lP.
      CALL dgemm('N', 'N', n_AO_pairs, nRI, n_grid_points, -1.0_dp, Phi_munu_l_scaled, n_AO_pairs, &
                 Z_prime_lP, n_grid_points, 1.0_dp, R_munu_P, n_AO_pairs)
      normalized_error = SUM(R_munu_P*R_munu_P)/Int_3c_norm2

      ! For λ>0, the residual derivative includes the response of the regularized coefficients:
      ! dE/dΦ = [-2 R Z^T - 2λ(R Y^T - Φ Y Z^T)] / ||Int_3c||^2.
      ALLOCATE (D_inverse_Z_prime_lP(n_grid_points, nRI), dE_dPhi_munu_l(n_AO_pairs, n_grid_points))
      D_inverse_Z_prime_lP(:, :) = Z_prime_lP
      CALL dpotrs('L', n_grid_points, nRI, D_ll, n_grid_points, D_inverse_Z_prime_lP, n_grid_points, lapack_info)
      IF (lapack_info /= 0) THEN
         DEALLOCATE (Phi_munu_l, Phi_munu_l_scaled, Int_3c_munu_P, dE_dPhi_munu_l, D_ll, mu_of_pair, nu_of_pair, &
                     R_munu_P, jacobi_scaling, D_inverse_Z_prime_lP, Z_prime_lP)
         CALL timestop(handle)
         RETURN
      END IF
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

      ALLOCATE (dPhi_munu_l_scaled(n_AO_pairs))
      DO grid_point = 1, n_grid_points
         jacobi_scale = jacobi_scaling(grid_point)
         Phi_munu_norm2 = SUM(Phi_munu_l(:, grid_point)*Phi_munu_l(:, grid_point))
         DO cartesian_component = 1, 3
            DO AO_pair = 1, n_AO_pairs
               mu = mu_of_pair(AO_pair)
               nu = nu_of_pair(AO_pair)
               symmetry_factor = SQRT(REAL(2 - MERGE(1, 0, mu == nu), dp))
               dPhi_munu_l = symmetry_factor*( &
                             dPhi_alpha_l_mu(cartesian_component, grid_point, mu)* &
                             Phi_l_mu(grid_point, nu) + &
                             Phi_l_mu(grid_point, mu)* &
                             dPhi_alpha_l_mu(cartesian_component, grid_point, nu))
               dPhi_munu_l_scaled(AO_pair) = jacobi_scale*dPhi_munu_l
            END DO
            IF (Phi_munu_norm2 > jacobi_floor) THEN
               scaling_projection = &
                  DOT_PRODUCT(Phi_munu_l(:, grid_point), dPhi_munu_l_scaled)/jacobi_scale
               dPhi_munu_l_scaled(:) = dPhi_munu_l_scaled - &
                                       jacobi_scale**3*Phi_munu_l(:, grid_point)*scaling_projection
            END IF
            coordinate_gradient(cartesian_component, grid_point) = &
               DOT_PRODUCT(dE_dPhi_munu_l(:, grid_point), dPhi_munu_l_scaled)
         END DO
      END DO

      max_abs_error = 0.0_dp
      DO P = 1, nRI
         DO AO_pair = 1, n_AO_pairs
            symmetry_factor = &
               SQRT(REAL(2 - MERGE(1, 0, mu_of_pair(AO_pair) == nu_of_pair(AO_pair)), dp))
            absolute_error = ABS(R_munu_P(AO_pair, P))/symmetry_factor
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

      INTEGER                                            :: atom, cartesian_component, cluster_atom, &
                                                            grid_point, grid_point_offset

      grid_point_offset = 0
      DO cluster_atom = 1, SIZE(cluster_atoms)
         atom = cluster_atoms(cluster_atom)
         DO grid_point = 1, SIZE(grids(atom)%raw_points, 2)
            DO cartesian_component = 1, 3
               coordinate_gradient(atom_coordinate_offsets(atom) + 3*(grid_point - 1) + cartesian_component) = &
                  coordinate_gradient(atom_coordinate_offsets(atom) + 3*(grid_point - 1) + cartesian_component) + &
                  cluster_coordinate_gradient(cartesian_component, grid_point_offset + grid_point)
            END DO
         END DO
         grid_point_offset = grid_point_offset + SIZE(grids(atom)%raw_points, 2)
      END DO
   END SUBROUTINE accumulate_atom_gradient

! **************************************************************************************************
!> \brief Pack all atom-centred grids into the single Cartesian vector required by L-BFGS.
!> \param grids Atom-centred RI-RS grids.
!> \param packed_grids Coordinates, atom offsets, and numerical coordinate bounds.
!>
!> SIZE(coordinates) = 3 Σ_A N_grid,A. Coordinates are ordered Cartesian component first,
!> then grid point, then atom, matching the column-major layout of raw_points(3,N_grid,A).
! **************************************************************************************************
   SUBROUTINE pack_atom_grids(grids, packed_grids)
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(IN)     :: grids
      TYPE(packed_atom_grids_type), INTENT(OUT)          :: packed_grids

      REAL(KIND=dp), PARAMETER :: maximum_grid_displacement = 100.0_dp/angstrom

      INTEGER                                            :: atom, atom_coordinate_count, &
                                                            n_grid_coordinates

      ALLOCATE (packed_grids%atom_coordinate_offsets(SIZE(grids)))
      n_grid_coordinates = 0
      DO atom = 1, SIZE(grids)
         packed_grids%atom_coordinate_offsets(atom) = n_grid_coordinates
         n_grid_coordinates = n_grid_coordinates + SIZE(grids(atom)%raw_points)
      END DO
      ALLOCATE (packed_grids%coordinates(n_grid_coordinates), &
                packed_grids%lower_bounds(n_grid_coordinates), &
                packed_grids%upper_bounds(n_grid_coordinates))

      DO atom = 1, SIZE(grids)
         atom_coordinate_count = SIZE(grids(atom)%raw_points)
         packed_grids%coordinates(packed_grids%atom_coordinate_offsets(atom) + 1: &
                                  packed_grids%atom_coordinate_offsets(atom) + atom_coordinate_count) = &
            RESHAPE(grids(atom)%raw_points, [atom_coordinate_count])
      END DO
      packed_grids%lower_bounds(:) = packed_grids%coordinates - maximum_grid_displacement
      packed_grids%upper_bounds(:) = packed_grids%coordinates + maximum_grid_displacement
   END SUBROUTINE pack_atom_grids

! **************************************************************************************************
!> \brief Restore the optimized Cartesian vector to the persistent atom-centred RI-RS grids.
!> \param packed_grids Flattened coordinates and their atom-to-vector mapping.
!> \param grids Atom-centred RI-RS grids to update.
! **************************************************************************************************
   SUBROUTINE unpack_atom_grids(packed_grids, grids)
      TYPE(packed_atom_grids_type), INTENT(IN)           :: packed_grids
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(INOUT)  :: grids

      INTEGER                                            :: atom, atom_coordinate_count

      CPASSERT(SIZE(packed_grids%atom_coordinate_offsets) == SIZE(grids))
      DO atom = 1, SIZE(grids)
         atom_coordinate_count = SIZE(grids(atom)%raw_points)
         grids(atom)%raw_points(:, :) = RESHAPE( &
                                        packed_grids%coordinates(packed_grids%atom_coordinate_offsets(atom) + 1: &
                                                                 packed_grids%atom_coordinate_offsets(atom) + &
                                                                 atom_coordinate_count), &
                                        SHAPE(grids(atom)%raw_points))
      END DO
   END SUBROUTINE unpack_atom_grids

END MODULE gw_optimize_ri_rs_grid
