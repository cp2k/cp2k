!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local-environment optimizer for atom-centred real-space RI grids.
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE gw_optimize_ri_rs_grid
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc,&
                                              use_perd_none
   USE cp_lbfgs,                        ONLY: setulb
   USE gw_integrals,                    ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE input_constants,                 ONLY: do_potential_truncated
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: libint_potential_type
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
   REAL(KIND=dp), PARAMETER, PRIVATE :: optimizer_accuracy = 1.0E-9_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: truncated_coulomb_cutoff = 3.0_dp/angstrom

   TYPE :: local_cluster_type
      INTEGER, ALLOCATABLE                 :: atoms(:)
      REAL(KIND=dp), ALLOCATABLE           :: three_center(:, :, :)
   END TYPE local_cluster_type

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief Interface for objective and gradient evaluations used by the grid optimizer.
!> \param x Flattened atom-relative grid coordinates.
!> \param f Objective value.
!> \param g Objective gradient with respect to x.
!> \param successful Whether the objective evaluation succeeded.
! **************************************************************************************************
      SUBROUTINE grid_objective_callback(x, f, g, successful)
         IMPORT :: dp
         REAL(KIND=dp), DIMENSION(:), INTENT(IN)         :: x
         REAL(KIND=dp), INTENT(OUT)                      :: f
         REAL(KIND=dp), DIMENSION(:), INTENT(OUT)        :: g
         LOGICAL, INTENT(OUT)                            :: successful
      END SUBROUTINE grid_objective_callback
   END INTERFACE

   PUBLIC :: optimize_ri_rs_grid

CONTAINS

! **************************************************************************************************
!> \brief Optimize the selected RI-RS grids stored in a GW band-structure environment.
!> \param qs_env Quickstep environment used by the upstream three-center-integral context.
!> \param bs_env GW environment containing all configuration, basis, geometry, parallel, and grid data.
! **************************************************************************************************
   SUBROUTINE optimize_ri_rs_grid(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'optimize_ri_rs_grid'

      INTEGER :: handle, iatom, icluster, ikind, n_variable, successful_clusters, unit_nr
      INTEGER                                            :: periodic(3)
      INTEGER, ALLOCATABLE                               :: ao_size(:), grid_offsets(:), ri_size(:)
      LOGICAL                                            :: successful
      REAL(KIND=dp)                                      :: f
      REAL(KIND=dp), ALLOCATABLE                         :: g(:), lower(:), upper(:), x(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gw_3c_ctx_type)                               :: integral_context
      TYPE(local_cluster_type), ALLOCATABLE              :: clusters(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp)                                      :: maximum_absolute_error

      CALL timeset(routineN, handle)
      NULLIFY (cell, para_env, particle_set)
      cell => bs_env%ri_rs%cell
      particle_set => bs_env%ri_rs%particle_set
      para_env => bs_env%para_env
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ALLOCATED(bs_env%basis_set_AO))
      CPASSERT(ALLOCATED(bs_env%basis_set_RI))
      CPASSERT(ALLOCATED(bs_env%ri_rs%grid_cache))

      IF (bs_env%ri_rs%grid_opt%max_iter < 1) THEN
         CPABORT("GRID_OPTIMIZATION%MAX_ITER must be positive")
      END IF
      IF (bs_env%ri_rs%grid_opt%cutoff_atomic_cluster <= 0.0_dp) THEN
         CPABORT("GRID_OPTIMIZATION%CUTOFF_ATOMIC_CLUSTER must be positive")
      END IF
      IF (bs_env%ri_rs%tikhonov < 0.0_dp) THEN
         CPABORT("RI_RS%TIKHONOV must not be negative")
      END IF

      CALL get_cell(cell, periodic=periodic)
      IF (ANY(periodic /= use_perd_none)) THEN
         CPABORT("GRID_OPTIMIZATION currently supports nonperiodic local environments only")
      END IF

      unit_nr = bs_env%unit_nr
      ALLOCATE (ao_size(SIZE(particle_set)), ri_size(SIZE(particle_set)))
      DO iatom = 1, SIZE(particle_set)
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CPASSERT(ASSOCIATED(bs_env%basis_set_AO(ikind)%gto_basis_set))
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(ikind)%gto_basis_set))
         ao_size(iatom) = bs_env%basis_set_AO(ikind)%gto_basis_set%nsgf
         ri_size(iatom) = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
      END DO
      IF (ANY(ao_size < 1) .OR. ANY(ri_size < 1)) THEN
         CPABORT("Every atom in GRID_OPTIMIZATION needs ORB and RI_AUX functions")
      END IF

      CPASSERT(SIZE(bs_env%ri_rs%grid_cache) == SIZE(particle_set))
      ALLOCATE (grid_offsets(SIZE(bs_env%ri_rs%grid_cache)))
      n_variable = 0
      DO iatom = 1, SIZE(bs_env%ri_rs%grid_cache)
         grid_offsets(iatom) = n_variable
         n_variable = n_variable + 3*SIZE(bs_env%ri_rs%grid_cache(iatom)%raw_points, 2)
      END DO
      BLOCK
         TYPE(libint_potential_type) :: potential
         potential%potential_type = do_potential_truncated
         potential%cutoff_radius = truncated_coulomb_cutoff
         potential%omega = 0.0_dp
         potential%filename = "t_c_g.dat"
         CALL gw_3c_ctx_create(integral_context, qs_env, potential, &
                               bs_env%basis_set_AO, bs_env%basis_set_AO, bs_env%basis_set_RI)
      END BLOCK
      CALL build_local_clusters(integral_context, particle_set, cell, ao_size, ri_size, &
                                bs_env%ri_rs%grid_opt%cutoff_atomic_cluster, para_env, clusters)

      ALLOCATE (x(n_variable), g(n_variable), lower(n_variable), upper(n_variable))
      CALL pack_atom_grids(bs_env%ri_rs%grid_cache, grid_offsets, x)
      lower(:) = x - 100.0_dp/angstrom
      upper(:) = x + 100.0_dp/angstrom
      CALL optimize_grid_coordinates(x, lower, upper, optimizer_accuracy, &
                                     bs_env%ri_rs%grid_opt%max_iter, grid_objective)
      CALL grid_objective(x, f, g, successful)
      IF (.NOT. successful) CPABORT("RI-RS grid optimization produced no regularized fit")

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') 'RI-RS grid optimization terminated'
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') 'Normalized 3C error:', f
         WRITE (unit_nr, '(T2,A,T72,ES9.1)') &
            'Maximum absolute 3C error:', maximum_absolute_error
      END IF

      DO iatom = 1, SIZE(bs_env%ri_rs%grid_cache)
         bs_env%ri_rs%grid_cache(iatom)%npts = &
            SIZE(bs_env%ri_rs%grid_cache(iatom)%raw_points, 2)
      END DO
      IF (unit_nr > 0) FLUSH (unit_nr)

      bs_env%ri_rs%Z_lP_exists = .FALSE.
      CALL gw_3c_ctx_release(integral_context)
      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief L-BFGS-B callback for the regularized three-centre-integral fitting error.
!> \param coordinates Flattened atom-relative grid coordinates.
!> \param value Mean normalized squared three-centre-integral error.
!> \param gradient Derivative of value with respect to coordinates.
!> \param valid Whether every local-cluster evaluation succeeded.
! **************************************************************************************************
      SUBROUTINE grid_objective(coordinates, value, gradient, valid)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coordinates
      REAL(KIND=dp), INTENT(OUT)                         :: value
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: gradient
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: all_valid
      REAL(KIND=dp)                                      :: cluster_maximum_absolute_error, &
                                                            cluster_value
      REAL(KIND=dp), ALLOCATABLE                         :: cluster_gradient(:, :)

         CALL unpack_atom_grids(coordinates, bs_env%ri_rs%grid_cache, grid_offsets)
         value = 0.0_dp
         gradient = 0.0_dp
         successful_clusters = 0
         maximum_absolute_error = 0.0_dp
         valid = .TRUE.
         DO icluster = 1, SIZE(clusters)
            CALL evaluate_local_cluster(clusters(icluster), bs_env%ri_rs%grid_cache, &
                                        bs_env, particle_set, cell, &
                                        cluster_value, cluster_gradient, &
                                        cluster_maximum_absolute_error, valid)
            IF (.NOT. valid) THEN
               IF (ALLOCATED(cluster_gradient)) DEALLOCATE (cluster_gradient)
               EXIT
            END IF
            value = value + cluster_value
            CALL accumulate_atom_gradient(clusters(icluster)%atoms, bs_env%ri_rs%grid_cache, &
                                          grid_offsets, cluster_gradient, gradient)
            successful_clusters = successful_clusters + 1
            maximum_absolute_error = &
               MAX(maximum_absolute_error, cluster_maximum_absolute_error)
            DEALLOCATE (cluster_gradient)
         END DO
         all_valid = MERGE(1, 0, valid)
         CALL para_env%sum(value)
         CALL para_env%sum(gradient)
         CALL para_env%sum(successful_clusters)
         CALL para_env%sum(all_valid)
         CALL para_env%max(maximum_absolute_error)
         valid = successful_clusters == SIZE(particle_set) .AND. all_valid == para_env%num_pe
         IF (.NOT. valid) THEN
            value = HUGE(value)
            gradient = 0.0_dp
            RETURN
         END IF
         value = value/REAL(SIZE(particle_set), dp)
         gradient = gradient/REAL(SIZE(particle_set), dp)
      END SUBROUTINE grid_objective

   END SUBROUTINE optimize_ri_rs_grid

! **************************************************************************************************
!> \brief Minimize the RI-RS grid objective with CP2K's bound-constrained L-BFGS implementation.
!> \param x Coordinates on entry and best coordinates found on return.
!> \param lower Lower bound for each coordinate.
!> \param upper Upper bound for each coordinate.
!> \param accuracy Projected-gradient convergence threshold.
!> \param max_evaluations Maximum number of objective evaluations.
!> \param objective Objective and gradient callback.
! **************************************************************************************************
   SUBROUTINE optimize_grid_coordinates(x, lower, upper, accuracy, max_evaluations, objective)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: x
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: lower, upper
      REAL(KIND=dp), INTENT(IN)                          :: accuracy
      INTEGER, INTENT(IN)                                :: max_evaluations
      PROCEDURE(grid_objective_callback)                 :: objective

      INTEGER, PARAMETER                                 :: memory = 7
      CHARACTER(LEN=60)                                  :: csave, task
      INTEGER                                            :: evaluations
      INTEGER, DIMENSION(3*SIZE(x))                      :: iwa
      INTEGER, DIMENSION(44)                             :: isave
      INTEGER, DIMENSION(SIZE(x))                        :: bound_type
      LOGICAL                                            :: evaluation_ok, have_best
      LOGICAL, DIMENSION(4)                              :: lsave
      REAL(KIND=dp)                                      :: best_f, f
      REAL(KIND=dp), DIMENSION(29)                       :: dsave
      REAL(KIND=dp), DIMENSION(SIZE(x))                  :: best_x, g
      REAL(KIND=dp), DIMENSION(2*memory*SIZE(x) + 5*SIZE(x) + 11*memory**2 + 8*memory) :: wa

      CPASSERT(SIZE(x) > 0)
      CPASSERT(ALL(SHAPE(lower) == SHAPE(x)) .AND. ALL(SHAPE(upper) == SHAPE(x)))
      CPASSERT(ALL(lower <= upper))
      CPASSERT(accuracy > 0.0_dp .AND. max_evaluations > 0)

      bound_type = 2
      task = 'START'
      csave = ''
      f = HUGE(f)
      g = 0.0_dp
      wa = 0.0_dp
      iwa = 0
      lsave = .FALSE.
      isave = 0
      dsave = 0.0_dp
      evaluations = 0
      have_best = .FALSE.
      best_f = HUGE(best_f)
      best_x = x

      DO
         CALL setulb(SIZE(x), memory, x, lower, upper, bound_type, f, g, &
                     0.0_dp, accuracy, wa, iwa, task, -1, csave, lsave, isave, dsave, -1.0_dp)
         IF (task(1:2) == 'FG') THEN
            IF (evaluations >= max_evaluations) EXIT
            CALL objective(x, f, g, evaluation_ok)
            evaluations = evaluations + 1
            IF (.NOT. evaluation_ok) EXIT
            IF (.NOT. have_best .OR. f < best_f) THEN
               have_best = .TRUE.
               best_f = f
               best_x = x
            END IF
         ELSE IF (task(1:5) == 'NEW_X') THEN
            CYCLE
         ELSE
            EXIT
         END IF
      END DO
      IF (have_best) x = best_x
   END SUBROUTINE optimize_grid_coordinates

! **************************************************************************************************
!> \brief Build local clusters and precompute exact (mu nu|P) values on their owning rank.
!> \param context Three-centre-integral evaluation context.
!> \param particle_set Atomic positions and kinds.
!> \param cell Simulation cell.
!> \param ao_size Number of orbital basis functions on each atom.
!> \param ri_size Number of auxiliary basis functions on each atom.
!> \param cutoff Radius of each atom-centred local cluster.
!> \param para_env MPI environment used to distribute cluster ownership.
!> \param clusters Rank-local clusters and their exact three-centre integrals.
! **************************************************************************************************
   SUBROUTINE build_local_clusters(context, particle_set, cell, ao_size, ri_size, cutoff, &
                                   para_env, clusters)
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: context
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_size, ri_size
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(local_cluster_type), ALLOCATABLE, INTENT(OUT) :: clusters(:)

      INTEGER                                            :: center, ia, iatom, ja, jatom, ka, katom, &
                                                            n_ao, n_cluster, n_ri
      INTEGER, ALLOCATABLE                               :: ao_offset(:), atom_buffer(:), &
                                                            ri_offset(:)
      LOGICAL                                            :: screened
      TYPE(gw_3c_ws_type)                                :: workspace

      n_cluster = COUNT([(MOD(center - 1, para_env%num_pe) == para_env%mepos, &
                          center=1, SIZE(particle_set))])
      ALLOCATE (clusters(n_cluster), atom_buffer(SIZE(particle_set)), &
                ao_offset(SIZE(particle_set)), ri_offset(SIZE(particle_set)))
      CALL gw_3c_ws_create(workspace, context)
      n_cluster = 0
      DO center = 1, SIZE(particle_set)
         IF (MOD(center - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         n_cluster = n_cluster + 1
         ia = 0
         DO iatom = 1, SIZE(particle_set)
            IF (SUM(pbc(particle_set(iatom)%r - particle_set(center)%r, cell)**2) <= cutoff**2) THEN
               ia = ia + 1
               atom_buffer(ia) = iatom
            END IF
         END DO
         ALLOCATE (clusters(n_cluster)%atoms(ia), source=atom_buffer(1:ia))
         ao_offset = 0
         ri_offset = 0
         n_ao = 0
         n_ri = 0
         DO ia = 1, SIZE(clusters(n_cluster)%atoms)
            iatom = clusters(n_cluster)%atoms(ia)
            ao_offset(iatom) = n_ao
            ri_offset(iatom) = n_ri
            n_ao = n_ao + ao_size(iatom)
            n_ri = n_ri + ri_size(iatom)
         END DO
         ALLOCATE (clusters(n_cluster)%three_center(n_ao, n_ao, n_ri), source=0.0_dp)
         DO ia = 1, SIZE(clusters(n_cluster)%atoms)
            iatom = clusters(n_cluster)%atoms(ia)
            DO ka = 1, SIZE(clusters(n_cluster)%atoms)
               katom = clusters(n_cluster)%atoms(ka)
               DO ja = 1, SIZE(clusters(n_cluster)%atoms)
                  jatom = clusters(n_cluster)%atoms(ja)
                  CALL build_3c_integral_block_ctx( &
                     clusters(n_cluster)%three_center, context, workspace, &
                     atom_j=jatom, atom_k=katom, atom_i=iatom, &
                     j_offset=ao_offset(jatom), k_offset=ao_offset(katom), &
                     i_offset=ri_offset(iatom), screened=screened)
               END DO
            END DO
         END DO
      END DO
      CALL gw_3c_ws_release(workspace)
   END SUBROUTINE build_local_clusters

! **************************************************************************************************
!> \brief Evaluate one complete cluster and return derivatives for its physical grid points.
!> \param cluster Local atoms and exact three-centre integrals.
!> \param grids Atom-centred RI-RS grids.
!> \param bs_env GW environment containing orbital basis sets and RI-RS parameters.
!> \param particle_set Atomic positions and kinds.
!> \param cell Simulation cell.
!> \param value Normalized squared three-centre-integral error.
!> \param gradient Derivative with respect to the cluster's physical grid coordinates.
!> \param maximum_absolute_error Largest absolute three-centre-integral error.
!> \param successful Whether the regularized fitting equations were solved.
! **************************************************************************************************
   SUBROUTINE evaluate_local_cluster(cluster, grids, bs_env, particle_set, cell, &
                                     value, gradient, maximum_absolute_error, successful)
      TYPE(local_cluster_type), INTENT(IN)               :: cluster
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(IN)     :: grids
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(OUT)                         :: value
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: gradient(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: maximum_absolute_error
      LOGICAL, INTENT(OUT)                               :: successful

      INTEGER                                            :: ao_offset, ia, iatom, ikind, n_ao, &
                                                            n_atom_ao, n_grid, point_offset
      REAL(KIND=dp), ALLOCATABLE                         :: dphi(:, :, :), phi(:, :), points(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      n_ao = SIZE(cluster%three_center, 1)
      n_grid = 0
      DO ia = 1, SIZE(cluster%atoms)
         iatom = cluster%atoms(ia)
         n_grid = n_grid + SIZE(grids(iatom)%raw_points, 2)
      END DO
      ALLOCATE (points(3, n_grid), phi(n_grid, n_ao), &
                dphi(3, n_grid, n_ao), gradient(3, n_grid))
      phi = 0.0_dp
      dphi = 0.0_dp
      point_offset = 0
      DO ia = 1, SIZE(cluster%atoms)
         iatom = cluster%atoms(ia)
         points(:, point_offset + 1:point_offset + SIZE(grids(iatom)%raw_points, 2)) = &
            SPREAD(particle_set(iatom)%r, 2, SIZE(grids(iatom)%raw_points, 2)) + &
            grids(iatom)%raw_points
         point_offset = point_offset + SIZE(grids(iatom)%raw_points, 2)
      END DO
      ao_offset = 0
      DO ia = 1, SIZE(cluster%atoms)
         iatom = cluster%atoms(ia)
         ikind = particle_set(iatom)%atomic_kind%kind_number
         basis => bs_env%basis_set_AO(ikind)%gto_basis_set
         n_atom_ao = basis%nsgf
         CALL evaluate_ao_basis_on_points( &
            phi(:, ao_offset + 1:ao_offset + n_atom_ao), points, basis, &
            particle_set(iatom)%r, cell, &
            dphi=dphi(:, :, ao_offset + 1:ao_offset + n_atom_ao))
         ao_offset = ao_offset + n_atom_ao
      END DO
      CALL evaluate_rirs_grid_cluster(phi, dphi, cluster%three_center, bs_env%ri_rs%tikhonov, &
                                      value, gradient, maximum_absolute_error, successful)
   END SUBROUTINE evaluate_local_cluster

! **************************************************************************************************
!> \brief Evaluate the regularized local-cluster three-centre fit and its coordinate gradient.
!>
!> The solve applies the same column Jacobi scaling and Tikhonov parameter as the production
!> Z_lP construction. The objective contains only the normalized three-centre-integral residual.
!> \param phi AO values, indexed (l,mu).
!> \param dphi Cartesian derivatives of AO values, indexed (alpha,l,mu).
!> \param three_center Exact three-centre integrals, indexed (mu,nu,P).
!> \param tikhonov Tikhonov parameter used by the production RI-RS solve.
!> \param value Normalized squared residual for this cluster.
!> \param gradient Analytic derivative d value/d r(alpha,l).
!> \param maximum_absolute_error Largest absolute error in an unweighted three-centre integral.
!> \param successful False if the regularized normal equations cannot be solved.
! **************************************************************************************************
   SUBROUTINE evaluate_rirs_grid_cluster(phi, dphi, three_center, tikhonov, value, gradient, &
                                         maximum_absolute_error, successful)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: dphi, three_center
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      REAL(KIND=dp), INTENT(OUT)                         :: value
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient
      REAL(KIND=dp), INTENT(OUT)                         :: maximum_absolute_error
      LOGICAL, INTENT(OUT)                               :: successful

      CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_rirs_grid_cluster'
      REAL(KIND=dp), PARAMETER                           :: jacobi_floor = 1.0E-16_dp

      INTEGER                                            :: alpha, handle, info, ipair, l, mu, n_ao, &
                                                            n_grid, n_pair, n_ri, nu, p
      INTEGER, ALLOCATABLE                               :: pair_mu(:), pair_nu(:)
      REAL(KIND=dp)                                      :: absolute_error, column_norm2, d_a, &
                                                            denom, factor, projection, scale
      REAL(KIND=dp), ALLOCATABLE :: a_matrix(:, :), a_scaled(:, :), b_matrix(:, :), &
         derivative_a(:, :), derivative_scaled_a(:), gram(:, :), residual(:, :), rhs(:, :), &
         scale_columns(:), y_matrix(:, :), yz_matrix(:, :), z_matrix(:, :)

      CALL timeset(routineN, handle)

      n_grid = SIZE(phi, 1)
      n_ao = SIZE(phi, 2)
      n_ri = SIZE(three_center, 3)
      n_pair = n_ao*(n_ao + 1)/2

      CPASSERT(SIZE(dphi, 1) == 3)
      CPASSERT(SIZE(dphi, 2) == n_grid)
      CPASSERT(SIZE(dphi, 3) == n_ao)
      CPASSERT(SIZE(three_center, 1) == n_ao)
      CPASSERT(SIZE(three_center, 2) == n_ao)
      CPASSERT(SIZE(gradient, 1) == 3)
      CPASSERT(SIZE(gradient, 2) == n_grid)
      CPASSERT(tikhonov >= 0.0_dp)

      value = HUGE(value)
      gradient = 0.0_dp
      maximum_absolute_error = HUGE(maximum_absolute_error)
      successful = .FALSE.
      IF (n_grid < 1 .OR. n_ri < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (a_matrix(n_pair, n_grid), a_scaled(n_pair, n_grid), &
                b_matrix(n_pair, n_ri), pair_mu(n_pair), pair_nu(n_pair), &
                scale_columns(n_grid))
      ipair = 0
      DO nu = 1, n_ao
         DO mu = 1, nu
            ipair = ipair + 1
            pair_mu(ipair) = mu
            pair_nu(ipair) = nu
            factor = SQRT(REAL(2 - MERGE(1, 0, mu == nu), dp))
            DO l = 1, n_grid
               a_matrix(ipair, l) = factor*phi(l, mu)*phi(l, nu)
            END DO
            DO p = 1, n_ri
               b_matrix(ipair, p) = factor*three_center(mu, nu, p)
            END DO
         END DO
      END DO

      denom = SUM(b_matrix*b_matrix)
      IF (denom <= TINY(1.0_dp)) THEN
         DEALLOCATE (a_matrix, a_scaled, b_matrix, pair_mu, pair_nu, scale_columns)
         CALL timestop(handle)
         RETURN
      END IF

      ! Production Z_lP uses D'=d(A^T A)d+lambda I with d_l=1/sqrt((A^T A)_ll).
      DO l = 1, n_grid
         column_norm2 = SUM(a_matrix(:, l)*a_matrix(:, l))
         scale_columns(l) = 1.0_dp/SQRT(MAX(column_norm2, jacobi_floor))
         a_scaled(:, l) = scale_columns(l)*a_matrix(:, l)
      END DO
      ALLOCATE (gram(n_grid, n_grid), rhs(n_grid, n_ri), z_matrix(n_grid, n_ri))
      CALL dgemm('T', 'N', n_grid, n_grid, n_pair, 1.0_dp, a_scaled, n_pair, &
                 a_scaled, n_pair, 0.0_dp, gram, n_grid)
      DO l = 1, n_grid
         gram(l, l) = gram(l, l) + tikhonov
      END DO
      CALL dgemm('T', 'N', n_grid, n_ri, n_pair, 1.0_dp, a_scaled, n_pair, &
                 b_matrix, n_pair, 0.0_dp, rhs, n_grid)
      CALL dpotrf('L', n_grid, gram, n_grid, info)
      IF (info /= 0) THEN
         DEALLOCATE (a_matrix, a_scaled, b_matrix, gram, pair_mu, pair_nu, rhs, &
                     scale_columns, z_matrix)
         CALL timestop(handle)
         RETURN
      END IF
      z_matrix(:, :) = rhs
      CALL dpotrs('L', n_grid, n_ri, gram, n_grid, z_matrix, n_grid, info)
      DEALLOCATE (rhs)
      IF (info /= 0) THEN
         DEALLOCATE (a_matrix, a_scaled, b_matrix, gram, pair_mu, pair_nu, &
                     scale_columns, z_matrix)
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (residual(n_pair, n_ri))
      residual(:, :) = b_matrix
      CALL dgemm('N', 'N', n_pair, n_ri, n_grid, -1.0_dp, a_scaled, n_pair, &
                 z_matrix, n_grid, 1.0_dp, residual, n_pair)
      value = SUM(residual*residual)/denom

      ! For lambda > 0, the residual-only derivative includes the response of the regularized
      ! coefficients. With Y=(A^T A+lambda I)^(-1) Z:
      ! dE/dA = [-2 R Z^T - 2 lambda (R Y^T - A Y Z^T)] / ||B||^2.
      ALLOCATE (y_matrix(n_grid, n_ri), derivative_a(n_pair, n_grid))
      y_matrix(:, :) = z_matrix
      CALL dpotrs('L', n_grid, n_ri, gram, n_grid, y_matrix, n_grid, info)
      IF (info /= 0) THEN
         DEALLOCATE (a_matrix, a_scaled, b_matrix, derivative_a, gram, pair_mu, pair_nu, &
                     residual, scale_columns, y_matrix, z_matrix)
         CALL timestop(handle)
         RETURN
      END IF
      CALL dgemm('N', 'T', n_pair, n_grid, n_ri, -2.0_dp/denom, residual, n_pair, &
                 z_matrix, n_grid, 0.0_dp, derivative_a, n_pair)
      IF (tikhonov > 0.0_dp) THEN
         CALL dgemm('N', 'T', n_pair, n_grid, n_ri, -2.0_dp*tikhonov/denom, residual, n_pair, &
                    y_matrix, n_grid, 1.0_dp, derivative_a, n_pair)
         ALLOCATE (yz_matrix(n_grid, n_grid))
         CALL dgemm('N', 'T', n_grid, n_grid, n_ri, 1.0_dp, y_matrix, n_grid, &
                    z_matrix, n_grid, 0.0_dp, yz_matrix, n_grid)
         CALL dgemm('N', 'N', n_pair, n_grid, n_grid, 2.0_dp*tikhonov/denom, a_scaled, n_pair, &
                    yz_matrix, n_grid, 1.0_dp, derivative_a, n_pair)
         DEALLOCATE (yz_matrix)
      END IF

      ALLOCATE (derivative_scaled_a(n_pair))
      DO l = 1, n_grid
         scale = scale_columns(l)
         column_norm2 = SUM(a_matrix(:, l)*a_matrix(:, l))
         DO alpha = 1, 3
            DO ipair = 1, n_pair
               mu = pair_mu(ipair)
               nu = pair_nu(ipair)
               factor = SQRT(REAL(2 - MERGE(1, 0, mu == nu), dp))
               d_a = factor*(dphi(alpha, l, mu)*phi(l, nu) + &
                             phi(l, mu)*dphi(alpha, l, nu))
               derivative_scaled_a(ipair) = scale*d_a
            END DO
            IF (column_norm2 > jacobi_floor) THEN
               projection = DOT_PRODUCT(a_matrix(:, l), derivative_scaled_a)/scale
               derivative_scaled_a(:) = derivative_scaled_a - &
                                        scale**3*a_matrix(:, l)*projection
            END IF
            gradient(alpha, l) = DOT_PRODUCT(derivative_a(:, l), derivative_scaled_a)
         END DO
      END DO

      maximum_absolute_error = 0.0_dp
      DO p = 1, n_ri
         DO ipair = 1, n_pair
            factor = SQRT(REAL(2 - MERGE(1, 0, pair_mu(ipair) == pair_nu(ipair)), dp))
            absolute_error = ABS(residual(ipair, p))/factor
            maximum_absolute_error = MAX(maximum_absolute_error, absolute_error)
         END DO
      END DO

      successful = .TRUE.
      DEALLOCATE (a_matrix, a_scaled, b_matrix, derivative_a, derivative_scaled_a, gram, &
                  pair_mu, pair_nu, residual, scale_columns, y_matrix, z_matrix)
      CALL timestop(handle)
   END SUBROUTINE evaluate_rirs_grid_cluster

! **************************************************************************************************
!> \brief Sum each physical occurrence into its atom-relative coordinate.
!> \param atoms Atom indices in the local cluster.
!> \param grids Atom-centred RI-RS grids.
!> \param grid_offsets Starting coordinate offset for each atom.
!> \param physical_gradient Gradient for the cluster's physical grid points.
!> \param gradient Global flattened atom-relative gradient to update.
! **************************************************************************************************
   SUBROUTINE accumulate_atom_gradient(atoms, grids, grid_offsets, physical_gradient, gradient)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atoms
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(IN)     :: grids
      INTEGER, DIMENSION(:), INTENT(IN)                  :: grid_offsets
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: physical_gradient
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: gradient

      INTEGER                                            :: alpha, ia, iatom, l, point_offset

      point_offset = 0
      DO ia = 1, SIZE(atoms)
         iatom = atoms(ia)
         DO l = 1, SIZE(grids(iatom)%raw_points, 2)
            DO alpha = 1, 3
               gradient(grid_offsets(iatom) + 3*(l - 1) + alpha) = &
                  gradient(grid_offsets(iatom) + 3*(l - 1) + alpha) + &
                  physical_gradient(alpha, point_offset + l)
            END DO
         END DO
         point_offset = point_offset + SIZE(grids(iatom)%raw_points, 2)
      END DO
   END SUBROUTINE accumulate_atom_gradient

! **************************************************************************************************
!> \brief Flatten atom-relative Cartesian coordinates into the optimizer vector.
!> \param grids Atom-centred RI-RS grids.
!> \param grid_offsets Starting coordinate offset for each atom.
!> \param coordinates Flattened atom-relative grid coordinates.
! **************************************************************************************************
   SUBROUTINE pack_atom_grids(grids, grid_offsets, coordinates)
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(IN)     :: grids
      INTEGER, DIMENSION(:), INTENT(IN)                  :: grid_offsets
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coordinates

      INTEGER                                            :: iatom, n_coordinate

      CPASSERT(SIZE(grid_offsets) == SIZE(grids))
      DO iatom = 1, SIZE(grids)
         n_coordinate = SIZE(grids(iatom)%raw_points)
         coordinates(grid_offsets(iatom) + 1:grid_offsets(iatom) + n_coordinate) = &
            RESHAPE(grids(iatom)%raw_points, [n_coordinate])
      END DO
   END SUBROUTINE pack_atom_grids

! **************************************************************************************************
!> \brief Restore atom-relative Cartesian coordinates from the optimizer vector.
!> \param coordinates Flattened atom-relative grid coordinates.
!> \param grids Atom-centred RI-RS grids to update.
!> \param grid_offsets Starting coordinate offset for each atom.
! **************************************************************************************************
   SUBROUTINE unpack_atom_grids(coordinates, grids, grid_offsets)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coordinates
      TYPE(rirs_grid_type), DIMENSION(:), INTENT(INOUT)  :: grids
      INTEGER, DIMENSION(:), INTENT(IN)                  :: grid_offsets

      INTEGER                                            :: iatom, n_coordinate

      CPASSERT(SIZE(grid_offsets) == SIZE(grids))
      DO iatom = 1, SIZE(grids)
         n_coordinate = SIZE(grids(iatom)%raw_points)
         grids(iatom)%raw_points(:, :) = RESHAPE( &
                                         coordinates(grid_offsets(iatom) + 1:grid_offsets(iatom) + n_coordinate), &
                                         SHAPE(grids(iatom)%raw_points))
      END DO
   END SUBROUTINE unpack_atom_grids

END MODULE gw_optimize_ri_rs_grid
