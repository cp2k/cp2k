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
MODULE ri_rs_grid_optimization
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_lbfgs,                        ONLY: setulb
   USE gw_utils_compute_integrals,      ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE kinds,                           ONLY: dp,&
                                              dp_size
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE ri_rs_fitting,                   ONLY: evaluate_ri_rs_fit
   USE ri_rs_grid_initialization,       ONLY: initialize_ri_rs_grid
   USE ri_rs_types,                     ONLY: rirs_grid_type
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points,&
                                              get_rirs_cluster_atoms
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: optimize_ri_rs_grid

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_grid_optimization'

   ! Share of the available memory (in GB of 1.0E9 bytes) for the cluster integrals.
   REAL(KIND=dp), PARAMETER :: integral_memory_share = 0.8_dp, bytes_per_GB = 1.0E9_dp

   TYPE, PRIVATE :: local_cluster_3c_integrals_type
      INTEGER, ALLOCATABLE                 :: atom_indices(:)
      REAL(KIND=dp), ALLOCATABLE           :: centers(:, :), Int_3c(:, :, :)
      REAL(KIND=dp)                        :: Int_3c_norm2
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
      TYPE(mp_para_env_type), POINTER                    :: atom_env

      CALL timeset(routineN, handle)

      ! Validate all state and user input before the initialization consumes them.
      CALL validate_ri_rs_grid_optimization_input(bs_env)

      ! setup Lebedev initial grid inside Voronoi volume, select points with Cholesky decomp.
      CALL initialize_ri_rs_grid(bs_env)

      ! Validate the initialized RI-RS state and report the start of the optimization.
      CALL prepare_ri_rs_grid_optimization(bs_env)

      ! C_A contains atom images within R_cut; distribute (μν|P) within its MPI group.
      CALL build_cluster_3c_int(bs_env, atom_env, cluster_3c_int)

      ! E_A = Σ_{μνP∈C_A} [(μν|P) - Σ_{l∈G_A} ϕ_μ(r_l) ϕ_ν(r_l) Z_lP^(A)]².
      ! Minimize E_loc(norm) = N_atom^(-1) Σ_A [E_A / Σ_{μνP∈C_A} |(μν|P)|²].
      CALL optimize_grid_coordinates(bs_env, atom_env, cluster_3c_int)
      CALL mp_para_env_release(atom_env)

      CALL timestop(handle)

   END SUBROUTINE optimize_ri_rs_grid

! **************************************************************************************************
!> \brief Validate state and user input needed to initialize and optimize RI-RS grids.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE validate_ri_rs_grid_optimization_input(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: iatom, ikind

      CPASSERT(ASSOCIATED(bs_env%ri_rs%cell))
      CPASSERT(ASSOCIATED(bs_env%ri_rs%particle_set))
      CPASSERT(ASSOCIATED(bs_env%para_env))
      CPASSERT(ALLOCATED(bs_env%basis_set_AO))
      CPASSERT(ALLOCATED(bs_env%basis_set_RI))
      CPASSERT(ALLOCATED(bs_env%sizes_AO))
      CPASSERT(SIZE(bs_env%ri_rs%particle_set) == bs_env%n_atom)
      CPASSERT(SIZE(bs_env%sizes_AO) == bs_env%n_atom)
      CPASSERT(SIZE(bs_env%basis_set_AO) == SIZE(bs_env%basis_set_RI))

      IF (bs_env%ri_rs%grid_opt%rs_ao_ratio <= 0.0_dp) THEN
         CPABORT("GRID_OPTIMIZATION%RS_AO_RATIO must be positive")
      END IF
      IF (bs_env%ri_rs%grid_opt%cutoff_atomic_cluster == 0.0_dp) THEN
         CPABORT("GRID_OPTIMIZATION%CUTOFF_ATOMIC_CLUSTER must be nonzero")
      END IF
      IF (bs_env%ri_rs%grid_opt%max_iter < 1) THEN
         CPABORT("GRID_OPTIMIZATION%MAX_ITER must be positive")
      END IF
      IF (bs_env%ri_rs%tikhonov < 0.0_dp) THEN
         CPABORT("RI_RS%TIKHONOV must not be negative")
      END IF

      DO iatom = 1, bs_env%n_atom
         ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(bs_env%basis_set_AO))
         CPASSERT(ASSOCIATED(bs_env%basis_set_AO(ikind)%gto_basis_set))
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(ikind)%gto_basis_set))
         IF (bs_env%sizes_AO(iatom) < 1 .OR. &
             bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf < 1) THEN
            CPABORT("Every atom in GRID_OPTIMIZATION needs ORB and RI_AUX functions")
         END IF
      END DO
   END SUBROUTINE validate_ri_rs_grid_optimization_input

! **************************************************************************************************
!> \brief Validate initialized RI-RS state and announce the grid optimization.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE prepare_ri_rs_grid_optimization(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: iatom, unit_nr

      CPASSERT(ALLOCATED(bs_env%ri_rs%atomic_grids))
      CPASSERT(SIZE(bs_env%ri_rs%atomic_grids) == bs_env%n_atom)

      DO iatom = 1, bs_env%n_atom
         CPASSERT(ALLOCATED(bs_env%ri_rs%atomic_grids(iatom)%raw_points))
         CPASSERT(SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 1) == 3)
         CPASSERT(SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2) > 0)
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
!> \param atom_env ...
!> \param cluster_3c_int Rank-local cluster three-centre integrals.
!> \param normalized_error Mean normalized squared three-centre-integral error.
!> \param coordinate_gradient Derivative of normalized_error with respect to coordinates.
!> \param fit_successful Whether every local-cluster evaluation succeeded.
!> \param max_abs_error Optional largest absolute three-centre-integral error.
! **************************************************************************************************
   SUBROUTINE grid_objective(grid_coordinates, atom_coordinate_offsets, bs_env, atom_env, &
                             cluster_3c_int, normalized_error, coordinate_gradient, &
                             fit_successful, max_abs_error)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: grid_coordinates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_coordinate_offsets
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: atom_env
      TYPE(local_cluster_3c_integrals_type), &
         DIMENSION(:), INTENT(IN)                        :: cluster_3c_int
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coordinate_gradient
      LOGICAL, INTENT(OUT)                               :: fit_successful
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: max_abs_error

      INTEGER                                            :: evaluated_cluster_count, iatom, &
                                                            icluster, ilocal, jatom, metadata(2), &
                                                            ngroups, nvalues, root, &
                                                            successful_ranks
      INTEGER, ALLOCATABLE                               :: cluster_atoms(:), gradient_offsets(:)
      LOGICAL                                            :: cluster_fit_successful
      REAL(KIND=dp)                                      :: cluster_max_abs_error, &
                                                            cluster_normalized_error
      REAL(KIND=dp), ALLOCATABLE :: cluster_coordinate_gradient(:, :), local_errors(:), &
         local_gradients(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      para_env => bs_env%para_env
      CPASSERT(SIZE(coordinate_gradient) == SIZE(grid_coordinates))
      CALL unpack_lbfgs_grid_coordinates(grid_coordinates, atom_coordinate_offsets, &
                                         bs_env%ri_rs%atomic_grids)
      normalized_error = 0.0_dp
      coordinate_gradient = 0.0_dp
      evaluated_cluster_count = 0
      IF (PRESENT(max_abs_error)) max_abs_error = 0.0_dp
      fit_successful = .TRUE.
      ALLOCATE (gradient_offsets(SIZE(cluster_3c_int) + 1), local_errors(SIZE(cluster_3c_int)))
      gradient_offsets(1) = 1
      DO icluster = 1, SIZE(cluster_3c_int)
         nvalues = 0
         DO jatom = 1, SIZE(cluster_3c_int(icluster)%atom_indices)
            iatom = cluster_3c_int(icluster)%atom_indices(jatom)
            nvalues = nvalues + 3*bs_env%ri_rs%atomic_grids(iatom)%npts
         END DO
         gradient_offsets(icluster + 1) = gradient_offsets(icluster) + nvalues
      END DO
      ALLOCATE (local_gradients(gradient_offsets(SIZE(cluster_3c_int) + 1) - 1))
      DO icluster = 1, SIZE(cluster_3c_int)
         IF (PRESENT(max_abs_error)) THEN
            CALL evaluate_local_cluster(cluster_3c_int(icluster), bs_env, atom_env, &
                                        cluster_normalized_error, cluster_coordinate_gradient, &
                                        cluster_fit_successful, cluster_max_abs_error)
            max_abs_error = MAX(max_abs_error, cluster_max_abs_error)
         ELSE
            CALL evaluate_local_cluster(cluster_3c_int(icluster), bs_env, atom_env, &
                                        cluster_normalized_error, cluster_coordinate_gradient, &
                                        cluster_fit_successful)
         END IF
         IF (.NOT. cluster_fit_successful) THEN
            fit_successful = .FALSE.
            IF (ALLOCATED(cluster_coordinate_gradient)) DEALLOCATE (cluster_coordinate_gradient)
            EXIT
         END IF
         IF (atom_env%mepos == 0) THEN
            local_errors(icluster) = cluster_normalized_error
            nvalues = SIZE(cluster_coordinate_gradient)
            local_gradients(gradient_offsets(icluster):gradient_offsets(icluster + 1) - 1) = &
               RESHAPE(cluster_coordinate_gradient, [nvalues])
            evaluated_cluster_count = evaluated_cluster_count + 1
         END IF
         DEALLOCATE (cluster_coordinate_gradient)
      END DO
      successful_ranks = MERGE(1, 0, fit_successful)
      CALL para_env%sum(evaluated_cluster_count)
      CALL para_env%sum(successful_ranks)
      IF (PRESENT(max_abs_error)) CALL para_env%max(max_abs_error)
      fit_successful = evaluated_cluster_count == bs_env%n_atom .AND. &
                       successful_ranks == para_env%num_pe
      IF (.NOT. fit_successful) THEN
         normalized_error = HUGE(normalized_error)
         coordinate_gradient = 0.0_dp
         RETURN
      END IF
      ! Add complete cluster contributions in atom order, independent of the rank assignment.
      ngroups = MIN(bs_env%n_atom, para_env%num_pe)
      DO iatom = 1, bs_env%n_atom
         root = MOD(iatom - 1, ngroups)
         ilocal = (iatom - 1)/ngroups + 1
         IF (para_env%mepos == root) THEN
            metadata = [SIZE(cluster_3c_int(ilocal)%atom_indices), &
                        gradient_offsets(ilocal + 1) - gradient_offsets(ilocal)]
            cluster_normalized_error = local_errors(ilocal)
         END IF
         CALL para_env%bcast(metadata, root)
         CALL para_env%bcast(cluster_normalized_error, root)
         ALLOCATE (cluster_atoms(metadata(1)), cluster_coordinate_gradient(3, metadata(2)/3))
         IF (para_env%mepos == root) THEN
            cluster_atoms(:) = cluster_3c_int(ilocal)%atom_indices
            cluster_coordinate_gradient(:, :) = RESHAPE( &
                                                local_gradients(gradient_offsets(ilocal):gradient_offsets(ilocal + 1) - 1), &
                                                SHAPE(cluster_coordinate_gradient))
         END IF
         CALL para_env%bcast(cluster_atoms, root)
         CALL para_env%bcast(cluster_coordinate_gradient, root)
         normalized_error = normalized_error + cluster_normalized_error
         CALL accumulate_atom_gradient(cluster_atoms, bs_env%ri_rs%atomic_grids, &
                                       atom_coordinate_offsets, cluster_coordinate_gradient, coordinate_gradient)
         DEALLOCATE (cluster_atoms, cluster_coordinate_gradient)
      END DO
      ! E_loc(norm) = (1/N_atom) Σ_A [E_A / Σ_{μνP ∈ C_A} |(μν|P)|²].
      normalized_error = normalized_error/REAL(bs_env%n_atom, dp)
      coordinate_gradient = coordinate_gradient/REAL(bs_env%n_atom, dp)
   END SUBROUTINE grid_objective

! **************************************************************************************************
!> \brief Drive L-BFGS, repeatedly evaluating E_loc and its gradient, and return the best trial.
!> \param bs_env ...
!> \param atom_env ...
!> \param cluster_3c_int Rank-local cluster three-centre integrals.
! **************************************************************************************************
   SUBROUTINE optimize_grid_coordinates(bs_env, atom_env, cluster_3c_int)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: atom_env
      TYPE(local_cluster_3c_integrals_type), &
         DIMENSION(:), INTENT(IN)                        :: cluster_3c_int

      CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_grid_coordinates'
      INTEGER, PARAMETER                                 :: lbfgs_history = 7
      REAL(KIND=dp), PARAMETER                           :: lbfgs_factr = 0.0_dp, &
                                                            lbfgs_pgtol = 1.0E-9_dp

      CHARACTER(LEN=60)                                  :: line_search_state, optimizer_task
      INTEGER                                            :: evaluation_count, handle, iatom, &
                                                            npoints, unit_nr
      INTEGER, ALLOCATABLE                               :: atom_coordinate_offsets(:), &
                                                            bound_types(:), integer_workspace(:)
      INTEGER, DIMENSION(44)                             :: integer_state
      LOGICAL                                            :: evaluation_successful, have_best
      LOGICAL, DIMENSION(4)                              :: logical_state
      REAL(KIND=dp)                                      :: best_normalized_error, evaluation_start, &
                                                            evaluation_time, max_abs_error, &
                                                            normalized_error, start_time
      REAL(KIND=dp), ALLOCATABLE :: best_grid_coordinates(:), coordinate_gradient(:), &
         grid_coordinates(:), lower_bounds(:), upper_bounds(:), workspace(:)
      REAL(KIND=dp), DIMENSION(29)                       :: real_state

      CALL timeset(routineN, handle)

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

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T69,I12)') &
            'RI-RS grid maximum evaluations (MAX_ITER):', bs_env%ri_rs%grid_opt%max_iter
         CALL m_flush(bs_env%unit_nr)
      END IF
      DO
         CALL setulb(SIZE(grid_coordinates), lbfgs_history, grid_coordinates, &
                     lower_bounds, upper_bounds, bound_types, &
                     normalized_error, coordinate_gradient, lbfgs_factr, lbfgs_pgtol, &
                     workspace, integer_workspace, optimizer_task, -1, line_search_state, &
                     logical_state, integer_state, real_state, -1.0_dp)
         IF (optimizer_task(1:2) == 'FG') THEN
            IF (evaluation_count >= bs_env%ri_rs%grid_opt%max_iter) EXIT
            evaluation_start = m_walltime()
            CALL grid_objective(grid_coordinates, atom_coordinate_offsets, bs_env, atom_env, &
                                cluster_3c_int, normalized_error, coordinate_gradient, &
                                evaluation_successful)
            evaluation_time = m_walltime() - evaluation_start
            evaluation_count = evaluation_count + 1
            IF (evaluation_count == 1 .AND. bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T61,ES20.12)') &
                  'RI-RS grid initial normalized error', normalized_error
            END IF
            IF (bs_env%unit_nr > 0) THEN
               IF (bs_env%ri_rs%grid_opt%max_iter < 50 .OR. &
                   ANY(evaluation_count == [1, 3, 10, 30, 100, 300, 1000, 3000, &
                                            10000, 30000, 100000, 300000, 1000000, &
                                            3000000, 10000000, 30000000, 100000000, &
                                            300000000, 1000000000]) .OR. &
                   evaluation_count == bs_env%ri_rs%grid_opt%max_iter) THEN
                  CALL write_ri_rs_grid_progress(bs_env%unit_nr, evaluation_count, &
                                                 bs_env%ri_rs%grid_opt%max_iter, &
                                                 evaluation_time, normalized_error, evaluation_successful)
               END IF
            END IF
            IF (.NOT. evaluation_successful) EXIT
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
      IF (.NOT. have_best) CPABORT("RI-RS grid optimization produced no regularized fit")
      ! grid_objective restores the best coordinates to atomic_grids and forms the explicit
      ! residual for the maximum error once.
      CALL grid_objective(best_grid_coordinates, atom_coordinate_offsets, bs_env, atom_env, &
                          cluster_3c_int, normalized_error, coordinate_gradient, &
                          evaluation_successful, max_abs_error)
      IF (.NOT. evaluation_successful) THEN
         CPABORT("RI-RS grid optimization produced no regularized fit")
      END IF

      DO iatom = 1, SIZE(bs_env%ri_rs%atomic_grids)
         bs_env%ri_rs%atomic_grids(iatom)%npts = &
            SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)
      END DO
      bs_env%ri_rs%Z_lP_exists = .FALSE.

      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         npoints = 0
         DO iatom = 1, SIZE(bs_env%ri_rs%atomic_grids)
            npoints = npoints + bs_env%ri_rs%atomic_grids(iatom)%npts
         END DO
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Normalized 3C error:', normalized_error
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Maximum absolute 3C error:', max_abs_error
         WRITE (unit_nr, '(T2,A,T69,I12)') 'Total grid points used for RI-RS:', npoints
         WRITE (unit_nr, '(T2,A,T69,F10.1,A)') &
            'RI-RS grid optimization completed, Exec. time:', m_walltime() - start_time, ' s'
         FLUSH (unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE optimize_grid_coordinates

! **************************************************************************************************
!> \brief Distribute the cluster integrals over the RI columns of each atom's MPI group.
!> \param bs_env ...
!> \param atom_env Communicator shared by the ranks processing the same atom clusters.
!> \param cluster_3c_int Each rank stores both AO triangles of its contiguous range of RI columns.
! **************************************************************************************************
   SUBROUTINE build_cluster_3c_int(bs_env, atom_env, cluster_3c_int)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: atom_env
      TYPE(local_cluster_3c_integrals_type), &
         ALLOCATABLE, INTENT(OUT)                        :: cluster_3c_int(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_cluster_3c_int'

      INTEGER :: first_ri, handle, iatom, icenter, icluster, iimage, ij, ik, ikind, imu, inu, ip, &
         jatom, katom, last_ri, my_group, nao, ngroups, nlocal, nri, p_first, p_last, patom
      INTEGER, ALLOCATABLE                               :: ao_offset(:), cells(:, :), ri_offset(:), &
                                                            sizes_ri(:)
      REAL(KIND=dp)                                      :: available_GB, integral_bytes, &
                                                            resident_bytes
      REAL(KIND=dp), ALLOCATABLE                         :: block(:, :, :)
      TYPE(gw_3c_ctx_type)                               :: integral_context
      TYPE(gw_3c_ws_type)                                :: workspace

      CALL timeset(routineN, handle)

      ngroups = MIN(bs_env%n_atom, bs_env%para_env%num_pe)
      my_group = MOD(bs_env%para_env%mepos, ngroups)
      ALLOCATE (atom_env)
      CALL atom_env%from_split(bs_env%para_env, my_group)
      nlocal = COUNT([(MOD(icenter - 1, ngroups) == my_group, icenter=1, bs_env%n_atom)])
      ALLOCATE (cluster_3c_int(nlocal), sizes_ri(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
         sizes_ri(iatom) = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
      END DO
      CALL gw_3c_ctx_create(integral_context, bs_env, bs_env%ri_metric, &
                            basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                            basis_i=bs_env%basis_set_RI)
      CALL gw_3c_ws_create(workspace, integral_context)
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, available_GB)
      resident_bytes = 0.0_dp
      icluster = 0
      DO icenter = 1, bs_env%n_atom
         IF (MOD(icenter - 1, ngroups) /= my_group) CYCLE
         icluster = icluster + 1
         CALL get_rirs_cluster_atoms(bs_env%ri_rs%particle_set, bs_env%ri_rs%cell, icenter, &
                                     bs_env%ri_rs%grid_opt%cutoff_atomic_cluster, &
                                     cluster_3c_int(icluster)%atom_indices, cells, &
                                     cluster_3c_int(icluster)%centers)
         ALLOCATE (ao_offset(SIZE(cells, 2)), ri_offset(SIZE(cells, 2)))
         nao = 0
         nri = 0
         DO iimage = 1, SIZE(cells, 2)
            iatom = cluster_3c_int(icluster)%atom_indices(iimage)
            ao_offset(iimage) = nao
            ri_offset(iimage) = nri
            nao = nao + bs_env%sizes_AO(iatom)
            nri = nri + sizes_ri(iatom)
         END DO
         first_ri = nri*atom_env%mepos/atom_env%num_pe + 1
         last_ri = nri*(atom_env%mepos + 1)/atom_env%num_pe
         integral_bytes = REAL(dp_size, dp)*REAL(nao, dp)**2*REAL(last_ri - first_ri + 1, dp)
         IF (available_GB > 0.0_dp .AND. &
             resident_bytes + integral_bytes > integral_memory_share*available_GB*bytes_per_GB) THEN
            CPABORT('Distributed RI-RS cluster integrals exceed available memory per rank')
         END IF
         ALLOCATE (cluster_3c_int(icluster)%Int_3c(nao, nao, first_ri:last_ri), SOURCE=0.0_dp)
         DO ip = 1, SIZE(cells, 2)
            patom = cluster_3c_int(icluster)%atom_indices(ip)
            ! Columns of this RI image owned by the rank.
            p_first = MAX(first_ri, ri_offset(ip) + 1)
            p_last = MIN(last_ri, ri_offset(ip) + sizes_ri(patom))
            IF (p_first > p_last) CYCLE
            DO ik = 1, SIZE(cells, 2)
               katom = cluster_3c_int(icluster)%atom_indices(ik)
               DO ij = 1, ik
                  jatom = cluster_3c_int(icluster)%atom_indices(ij)
                  ALLOCATE (block(bs_env%sizes_AO(jatom), bs_env%sizes_AO(katom), sizes_ri(patom)), SOURCE=0.0_dp)
                  CALL build_3c_integral_block_ctx(block, integral_context, workspace, &
                                                   atom_j=jatom, atom_k=katom, atom_i=patom, &
                                                   cell_j=cells(:, ij), cell_k=cells(:, ik), cell_i=cells(:, ip))
                  DO inu = 1, bs_env%sizes_AO(katom)
                     ! Within a diagonal image block keep only the upper triangle mu <= nu.
                     cluster_3c_int(icluster)%Int_3c(ao_offset(ij) + 1: &
                                                     ao_offset(ij) + MERGE(inu, bs_env%sizes_AO(jatom), ij == ik), &
                                                     ao_offset(ik) + inu, p_first:p_last) = &
                        block(1:MERGE(inu, bs_env%sizes_AO(jatom), ij == ik), inu, &
                              p_first - ri_offset(ip):p_last - ri_offset(ip))
                  END DO
                  DEALLOCATE (block)
               END DO
            END DO
         END DO
         ! ||B||² of the upper triangle with off-diagonal pairs counted twice.
         cluster_3c_int(icluster)%Int_3c_norm2 = 0.0_dp
         DO ip = first_ri, last_ri
            DO inu = 1, nao
               cluster_3c_int(icluster)%Int_3c_norm2 = cluster_3c_int(icluster)%Int_3c_norm2 + &
                                                       cluster_3c_int(icluster)%Int_3c(inu, inu, ip)**2 + &
                                                       2.0_dp*SUM(cluster_3c_int(icluster)%Int_3c(1:inu - 1, inu, ip)**2)
            END DO
         END DO
         CALL atom_env%sum(cluster_3c_int(icluster)%Int_3c_norm2)
         ! The evaluator multiplies by the full symmetric B_P (one DGEMM per RI function).
         DO ip = first_ri, last_ri
            DO inu = 1, nao
               DO imu = inu + 1, nao
                  cluster_3c_int(icluster)%Int_3c(imu, inu, ip) = cluster_3c_int(icluster)%Int_3c(inu, imu, ip)
               END DO
            END DO
         END DO
         resident_bytes = resident_bytes + integral_bytes
         DEALLOCATE (ao_offset, ri_offset, cells)
      END DO
      CALL gw_3c_ws_release(workspace)
      CALL gw_3c_ctx_release(integral_context)

      CALL timestop(handle)

   END SUBROUTINE build_cluster_3c_int

! **************************************************************************************************
!> \brief Evaluate ϕ_μ(r_l) on G_A and compute the normalized local error E_A and its gradient.
!> \param cluster_3c_int Local atom indices and exact three-centre integrals.
!> \param bs_env ...
!> \param atom_env ...
!> \param normalized_error Normalized squared three-centre-integral error.
!> \param coordinate_gradient Derivative with respect to the cluster's physical grid coordinates.
!> \param fit_successful Whether the regularized fitting equations were solved.
!> \param max_abs_error Optional largest absolute three-centre-integral error.
! **************************************************************************************************
   SUBROUTINE evaluate_local_cluster(cluster_3c_int, bs_env, atom_env, normalized_error, coordinate_gradient, &
                                     fit_successful, max_abs_error)
      TYPE(local_cluster_3c_integrals_type), INTENT(IN)  :: cluster_3c_int
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: atom_env
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: coordinate_gradient(:, :)
      LOGICAL, INTENT(OUT)                               :: fit_successful
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: max_abs_error

      INTEGER                                            :: AO_offset, grid_point_offset, iatom, &
                                                            icartesian, icluster_atom, igrid, &
                                                            ikind, n_grid_points, nAO, nAO_atom
      REAL(KIND=dp), ALLOCATABLE                         :: ao_gradient(:, :), &
                                                            dPhi_alpha_l_mu(:, :, :), &
                                                            grid_points(:, :), Phi_l_mu(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: AO_basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      cell => bs_env%ri_rs%cell
      particle_set => bs_env%ri_rs%particle_set
      nAO = SUM(bs_env%sizes_AO(cluster_3c_int%atom_indices))
      n_grid_points = 0
      DO icluster_atom = 1, SIZE(cluster_3c_int%atom_indices)
         iatom = cluster_3c_int%atom_indices(icluster_atom)
         n_grid_points = n_grid_points + SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)
      END DO
      ALLOCATE (grid_points(3, n_grid_points), Phi_l_mu(n_grid_points, nAO), &
                dPhi_alpha_l_mu(3, n_grid_points, nAO), ao_gradient(n_grid_points, nAO), &
                coordinate_gradient(3, n_grid_points))
      Phi_l_mu = 0.0_dp
      dPhi_alpha_l_mu = 0.0_dp
      grid_point_offset = 0
      DO icluster_atom = 1, SIZE(cluster_3c_int%atom_indices)
         iatom = cluster_3c_int%atom_indices(icluster_atom)
         grid_points(:, grid_point_offset + 1:grid_point_offset + &
                     SIZE(bs_env%ri_rs%atomic_grids(iatom)%raw_points, 2)) = &
            SPREAD(cluster_3c_int%centers(:, icluster_atom) - cluster_3c_int%centers(:, 1), 2, &
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
                                          grid_points, AO_basis, &
                                          cluster_3c_int%centers(:, icluster_atom) - cluster_3c_int%centers(:, 1), cell, &
                                          dphi=dPhi_alpha_l_mu(:, :, AO_offset + 1:AO_offset + nAO_atom), &
                                          minimum_image=.FALSE.)
         AO_offset = AO_offset + nAO_atom
      END DO
      CALL evaluate_ri_rs_fit(Phi_l_mu, cluster_3c_int%Int_3c, cluster_3c_int%Int_3c_norm2, &
                              bs_env%ri_rs%tikhonov, atom_env, normalized_error, fit_successful, &
                              ao_gradient, max_abs_error)
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_grid_points, coordinate_gradient, ao_gradient, dPhi_alpha_l_mu) &
      !$OMP PRIVATE(igrid, icartesian) SCHEDULE(STATIC)
      DO igrid = 1, n_grid_points
         DO icartesian = 1, 3
            coordinate_gradient(icartesian, igrid) = &
               DOT_PRODUCT(ao_gradient(igrid, :), dPhi_alpha_l_mu(icartesian, igrid, :))
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE evaluate_local_cluster

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

! **************************************************************************************************
!> \brief Print one compact line for a grid-objective evaluation.
!> \param unit_nr ...
!> \param evaluation_count ...
!> \param max_iter ...
!> \param evaluation_time ...
!> \param normalized_error ...
!> \param evaluation_successful ...
! **************************************************************************************************
   SUBROUTINE write_ri_rs_grid_progress(unit_nr, evaluation_count, max_iter, &
                                        evaluation_time, normalized_error, evaluation_successful)
      INTEGER, INTENT(IN)                                :: unit_nr, evaluation_count, max_iter
      REAL(KIND=dp), INTENT(IN)                          :: evaluation_time, normalized_error
      LOGICAL, INTENT(IN)                                :: evaluation_successful

      CHARACTER(LEN=10)                                  :: count_text, error_text, max_iter_text
      CHARACTER(LEN=64)                                  :: progress_line
      CHARACTER(LEN=8)                                   :: evaluation_text
      INTEGER                                            :: iexponent

      WRITE (max_iter_text, '(I0)') max_iter
      WRITE (count_text, '(I10)') evaluation_count
      count_text = count_text(11 - LEN_TRIM(max_iter_text):)
      WRITE (evaluation_text, '(F8.1)') evaluation_time
      IF (INDEX(evaluation_text, '*') > 0) WRITE (evaluation_text, '(ES8.1E3)') evaluation_time
      IF (evaluation_successful) THEN
         WRITE (error_text, '(ES10.1E3)') normalized_error
         error_text = ADJUSTL(error_text)
         iexponent = INDEX(error_text, 'E')
         IF (iexponent > 0) THEN
            IF (error_text(iexponent + 2:iexponent + 2) == '0') THEN
               error_text = error_text(:iexponent + 1)//error_text(iexponent + 3:)
            END IF
         END IF
      ELSE
         error_text = 'no fit'
      END IF
      progress_line = 'RI-RS grid eval '//TRIM(count_text)//' / '//TRIM(max_iter_text)// &
                      ',   3C error: '//TRIM(ADJUSTL(error_text))
      IF (LEN_TRIM(progress_line) > 56) progress_line = 'RI-RS eval '//progress_line(16:)
      WRITE (unit_nr, '(T2,A,T59,A,A,A)') TRIM(progress_line), 'Exec. time: ', evaluation_text, ' s'
      CALL m_flush(unit_nr)
   END SUBROUTINE write_ri_rs_grid_progress

END MODULE ri_rs_grid_optimization
