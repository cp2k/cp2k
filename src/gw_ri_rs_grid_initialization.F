!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialize atom-owned RI-RS grids by Cholesky selection from Lebedev grids.
! **************************************************************************************************
MODULE gw_ri_rs_grid_initialization
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE gw_ri_rs_utils,                  ONLY: evaluate_ao_basis_on_points,&
                                              filter_grid_to_voronoi,&
                                              get_rirs_cluster_atoms
   USE input_constants,                 ONLY: do_gapw_log
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              deallocate_grid_atom,&
                                              grid_atom_type
   USE util,                            ONLY: locate
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_grid_initialization'

   PUBLIC :: initialize_ri_rs_grid

CONTAINS

! **************************************************************************************************
!> \brief Construct the initial RI-RS grid of every atom.
!>
!> For every atom B, the Lebedev grid points are
!>
!> (1)  r_l = R_B + ρ_s Ω_a,
!>
!> where ρ_s is a radial quadrature point and Ω_a is a Lebedev direction. Lebedev grids are
!> tabulated according to the angular degree L up to which they integrate exactly. To integrate a
!> three-centre integral (μν|P), the smallest available Lebedev grid satisfying
!>
!> (2)  L >= 2 l_AO,max + l_RI,max + ΔL
!>
!> is used. Here, l_AO,max is the maximum angular momentum of the atomic AO basis functions ϕ_μ,
!> l_RI,max is the maximum angular momentum of the auxiliary basis functions φ_P, and the internal
!> angular buffer is ΔL. Enough radial points ρ_s are used that the Lebedev grid of every atom B
!> contains at least
!>
!> (3)  N_initial^B >= α_initial N_AO^B,
!>
!> points, where N_AO^B is the number of atomic orbitals on atom B. The requested initial RI-RS
!> grid of atom A contains
!>
!> (4)  N_R^A = ceil(RS_AO_RATIO * N_AO^A)
!>
!> points. Around every atom A, a cluster is defined as
!>
!> (5)  C_A = {B : |R_B - R_A| <= CUTOFF_ATOMIC_CLUSTER}.
!>
!> The Cholesky selection uses all Lebedev grid points of every atom B in C_A. The molecular
!> Voronoi cell of atom A is
!>
!> (6)  V_A = {r_l : |r_l - R_A| <= |r_l - R_B| for every atom B}.
!>
!> For cluster C_A, Cholesky selection is performed on
!>
!> (7)  D_ll' = [Σ_(μ in C_A) ϕ_μ(r_l) ϕ_μ(r_l')]^2.
!>
!> The first grid point is the point with the largest diagonal element,
!>
!> (8)  d_l^(0) = D_ll,                 q_1 = arg max_l d_l^(0).
!>
!> For every selected point q_k, the Cholesky column and diagonal are updated according to
!>
!> (9)  L_lk = [D_lq_k - Σ_(j<k) L_lj L_q_kj]/sqrt(d_q_k^(k-1)),
!>
!> (10) d_l^(k) = max(0, d_l^(k-1) - L_lk^2),
!>
!> and the next point q_(k+1) is the point with the largest d_l^(k). Points outside V_A take part
!> in Eqs. (7)-(10), but only selected points inside V_A are placed into the RI-RS grid of atom A.
!> Each cluster C_A is processed independently on one MPI rank; only the completed atom grids are
!> communicated.
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************
   SUBROUTINE initialize_ri_rs_grid(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_ri_rs_grid'
      INTEGER, PARAMETER :: l_additional = 2, radial_quadrature = do_gapw_log
      REAL(KIND=dp), PARAMETER                           :: rs_ao_ratio_initial = 30.0_dp

      INTEGER                                            :: handle, ikind
      TYPE(rirs_grid_type), ALLOCATABLE                  :: radial_lebedev_grids(:)

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%grid_opt%rs_ao_ratio <= rs_ao_ratio_initial)

      ! Build the Lebedev grids of Eqs. (1)-(3) for all elements in the calculation
      ALLOCATE (radial_lebedev_grids(SIZE(bs_env%basis_set_AO)))
      DO ikind = 1, SIZE(bs_env%basis_set_AO)
         CALL build_lebedev_grid(bs_env%basis_set_AO(ikind)%gto_basis_set, &
                                 bs_env%basis_set_RI(ikind)%gto_basis_set, &
                                 rs_ao_ratio_initial, l_additional, radial_quadrature, &
                                 radial_lebedev_grids(ikind))
      END DO

      ! Apply Eqs. (5)-(10) for every atom and keep points inside Voronoi volume, Eq. (6)
      CALL cholesky_selection_and_voronoi_filtering(bs_env, radial_lebedev_grids)

      CALL broadcast_ri_rs_grids(bs_env)

      CALL timestop(handle)

   END SUBROUTINE initialize_ri_rs_grid

! **************************************************************************************************
!> \brief Build one fixed-orientation Lebedev grid according to Eqs. (1)-(3).
!> \param ao ...
!> \param ri ...
!> \param ratio α_initial in Eq. (3).
!> \param l_additional ΔL in Eq. (2).
!> \param radial_quadrature ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE build_lebedev_grid(ao, ri, ratio, l_additional, radial_quadrature, grid)
      TYPE(gto_basis_set_type), POINTER                  :: ao, ri
      REAL(KIND=dp), INTENT(IN)                          :: ratio
      INTEGER, INTENT(IN)                                :: l_additional, radial_quadrature
      TYPE(rirs_grid_type), INTENT(OUT)                  :: grid

      INTEGER                                            :: degree, ir, nang, nrad, offset, rule
      TYPE(grid_atom_type), POINTER                      :: radial_grid

      CPASSERT(ASSOCIATED(ao) .AND. ASSOCIATED(ri))
      degree = 2*MAXVAL(ao%lmax) + MAXVAL(ri%lmax) + l_additional
      rule = get_number_of_lebedev_grid(l=degree)
      nang = lebedev_grid(rule)%n
      nrad = MAX(2, CEILING(ratio*REAL(ao%nsgf, dp)/REAL(nang, dp)))

      NULLIFY (radial_grid)
      CALL allocate_grid_atom(radial_grid)
      CALL create_grid_atom(radial_grid, nrad, nang, 0, rule, radial_quadrature)
      grid%npts = nrad*nang
      ALLOCATE (grid%raw_points(3, grid%npts))
      DO ir = 1, nrad
         offset = (ir - 1)*nang
         grid%raw_points(:, offset + 1:offset + nang) = &
            radial_grid%rad(ir)*lebedev_grid(rule)%r
      END DO
      CALL deallocate_grid_atom(radial_grid)
   END SUBROUTINE build_lebedev_grid

! **************************************************************************************************
!> \brief Apply Eqs. (5)-(10) independently for every atom and keep points inside Eq. (6).
!> \param bs_env ...
!> \param radial_lebedev_grids ...
! **************************************************************************************************
   SUBROUTINE cholesky_selection_and_voronoi_filtering(bs_env, radial_lebedev_grids)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rirs_grid_type), INTENT(IN)                   :: radial_lebedev_grids(:)

      INTEGER                                            :: atom_mpi_rank, iatom, ikind, npoints

      CPASSERT(.NOT. ALLOCATED(bs_env%ri_rs%atomic_grids))
      ALLOCATE (bs_env%ri_rs%atomic_grids(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
         npoints = CEILING(bs_env%ri_rs%grid_opt%rs_ao_ratio* &
                           bs_env%basis_set_AO(ikind)%gto_basis_set%nsgf)
         bs_env%ri_rs%atomic_grids(iatom)%npts = npoints
         ALLOCATE (bs_env%ri_rs%atomic_grids(iatom)%raw_points(3, npoints))

         atom_mpi_rank = MOD(iatom - 1, bs_env%para_env%num_pe)
         IF (atom_mpi_rank == bs_env%para_env%mepos) THEN
            CALL select_cholesky_grid_iatom(bs_env, radial_lebedev_grids, iatom)
         END IF
      END DO
   END SUBROUTINE cholesky_selection_and_voronoi_filtering

! **************************************************************************************************
!> \brief Select one atom grid with Eqs. (5)-(10) and the molecular Voronoi cell in Eq. (6).
!> \param bs_env ...
!> \param radial_lebedev_grids ...
!> \param iatom ...
! **************************************************************************************************
   SUBROUTINE select_cholesky_grid_iatom(bs_env, radial_lebedev_grids, iatom)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rirs_grid_type), INTENT(IN)                   :: radial_lebedev_grids(:)
      INTEGER, INTENT(IN)                                :: iatom

      CHARACTER(LEN=3*default_string_length)             :: failure
      INTEGER                                            :: igrid_point, npoints
      INTEGER, ALLOCATABLE :: ao_point_indices(:), atom_n_ao(:), atom_point_offsets(:), &
         atom_value_offsets(:), cluster_atoms(:), selected_indices(:)
      LOGICAL, ALLOCATABLE                               :: inside_voronoi(:)
      REAL(KIND=dp), ALLOCATABLE                         :: ao_values(:), cluster_points(:, :)

      CALL get_rirs_cluster_atoms(bs_env%ri_rs%particle_set, bs_env%ri_rs%cell, iatom, &
                                  bs_env%ri_rs%grid_opt%cutoff_atomic_cluster, cluster_atoms)
      CALL collect_cluster_lebedev_points(bs_env, radial_lebedev_grids, iatom, &
                                          cluster_atoms, cluster_points)
      npoints = SIZE(cluster_points, 2)
      ALLOCATE (inside_voronoi(npoints))
      CALL filter_grid_to_voronoi(cluster_points, iatom, bs_env%ri_rs%particle_set, &
                                  mask=inside_voronoi)
      IF (COUNT(inside_voronoi) < bs_env%ri_rs%atomic_grids(iatom)%npts) THEN
         WRITE (failure, '(A,I0,A,I0,A,I0,A)') 'Atom ', iatom, ': only ', COUNT(inside_voronoi), &
            ' Voronoi grid points for ', bs_env%ri_rs%atomic_grids(iatom)%npts, &
            ' points; reduce RS_AO_RATIO.'
         CPABORT(TRIM(failure))
      END IF

      CALL evaluate_cluster_ao_values(bs_env, iatom, cluster_atoms, cluster_points, &
                                      ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                      atom_n_ao, ao_values)
      ALLOCATE (selected_indices(bs_env%ri_rs%atomic_grids(iatom)%npts))
      CALL select_cholesky_grid_points(ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                       atom_n_ao, ao_values, npoints, &
                                       inside_voronoi, selected_indices)
      IF (ANY(selected_indices == 0)) THEN
         WRITE (failure, '(A,I0,A,I0,A,I0,A)') 'Atom ', iatom, ': Cholesky selection retained ', &
            COUNT(selected_indices > 0), ' of ', SIZE(selected_indices), &
            ' points before numerical rank exhaustion. Reduce RS_AO_RATIO.'
         CPABORT(TRIM(failure))
      END IF
      DO igrid_point = 1, SIZE(selected_indices)
         bs_env%ri_rs%atomic_grids(iatom)%raw_points(:, igrid_point) = &
            cluster_points(:, selected_indices(igrid_point))
      END DO
   END SUBROUTINE select_cholesky_grid_iatom

! **************************************************************************************************
!> \brief Collect all cluster Lebedev points used in Eqs. (5) and (7).
!> \param bs_env ...
!> \param radial_lebedev_grids ...
!> \param iatom ...
!> \param cluster_atoms ...
!> \param cluster_points ...
! **************************************************************************************************
   SUBROUTINE collect_cluster_lebedev_points(bs_env, radial_lebedev_grids, iatom, &
                                             cluster_atoms, cluster_points)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rirs_grid_type), INTENT(IN)                   :: radial_lebedev_grids(:)
      INTEGER, INTENT(IN)                                :: iatom, cluster_atoms(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: cluster_points(:, :)

      INTEGER                                            :: cluster_iatom, i_cluster_atom, ikind, n, &
                                                            npoints, offset

      npoints = 0
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         ikind = bs_env%ri_rs%particle_set(cluster_atoms(i_cluster_atom))%atomic_kind%kind_number
         npoints = npoints + radial_lebedev_grids(ikind)%npts
      END DO
      ALLOCATE (cluster_points(3, npoints))

      offset = 0
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         cluster_iatom = cluster_atoms(i_cluster_atom)
         ikind = bs_env%ri_rs%particle_set(cluster_iatom)%atomic_kind%kind_number
         n = radial_lebedev_grids(ikind)%npts
         cluster_points(:, offset + 1:offset + n) = radial_lebedev_grids(ikind)%raw_points + &
                                                    SPREAD(bs_env%ri_rs%particle_set(cluster_iatom)%r - &
                                                           bs_env%ri_rs%particle_set(iatom)%r, 2, n)
         offset = offset + n
      END DO
   END SUBROUTINE collect_cluster_lebedev_points

! **************************************************************************************************
!> \brief Evaluate ϕ_μ(r_l) needed for D_ll' in Eq. (7), using AO locality. Atom b uses the
!>        point rows atom_point_offsets(b):atom_point_offsets(b+1)-1. Its column-major block
!>        ϕ_μ(r_l) has shape (number of point rows, atom_n_ao(b)) and starts at
!>        atom_value_offsets(b) in ao_values.
!> \param bs_env ...
!> \param iatom ...
!> \param cluster_atoms ...
!> \param cluster_points ...
!> \param ao_point_indices Cluster-point index for every stored point row.
!> \param atom_point_offsets First stored point row for each atom, followed by the final bound.
!> \param atom_value_offsets First AO value for each atom, followed by the final bound.
!> \param atom_n_ao Number of AO functions for each atom.
!> \param ao_values Compact atom-local values ϕ_μ(r_l).
! **************************************************************************************************
   SUBROUTINE evaluate_cluster_ao_values(bs_env, iatom, cluster_atoms, cluster_points, &
                                         ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                         atom_n_ao, ao_values)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom, cluster_atoms(:)
      REAL(KIND=dp), INTENT(IN)                          :: cluster_points(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: ao_point_indices(:), &
                                                            atom_point_offsets(:), &
                                                            atom_value_offsets(:), atom_n_ao(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ao_values(:)

      INTEGER :: cluster_iatom, first_point, first_value, i_cluster_atom, igrid_point, ikind, &
         last_point, last_value, nactive, npoints
      INTEGER, ALLOCATABLE                               :: point_indices(:)
      LOGICAL, ALLOCATABLE                               :: active(:)
      REAL(KIND=dp), ALLOCATABLE                         :: absolute_points(:, :), &
                                                            atom_points(:, :), atom_values(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: ao

      npoints = SIZE(cluster_points, 2)
      ALLOCATE (absolute_points(3, npoints))
      absolute_points(:, :) = cluster_points + &
                              SPREAD(bs_env%ri_rs%particle_set(iatom)%r, 2, npoints)
      ALLOCATE (point_indices(npoints))
      point_indices(:) = [(igrid_point, igrid_point=1, npoints)]
      ALLOCATE (active(npoints), atom_n_ao(SIZE(cluster_atoms)), &
                atom_point_offsets(SIZE(cluster_atoms) + 1), &
                atom_value_offsets(SIZE(cluster_atoms) + 1))

      atom_point_offsets(1) = 1
      atom_value_offsets(1) = 1
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         cluster_iatom = cluster_atoms(i_cluster_atom)
         ikind = bs_env%ri_rs%particle_set(cluster_iatom)%atomic_kind%kind_number
         ao => bs_env%basis_set_AO(ikind)%gto_basis_set
         CPASSERT(ao%kind_radius > 0.0_dp)
         active(:) = SUM((absolute_points - &
                          SPREAD(bs_env%ri_rs%particle_set(cluster_iatom)%r, 2, npoints))**2, DIM=1) &
                     <= ao%kind_radius**2
         nactive = COUNT(active)
         atom_n_ao(i_cluster_atom) = ao%nsgf
         atom_point_offsets(i_cluster_atom + 1) = atom_point_offsets(i_cluster_atom) + nactive
         atom_value_offsets(i_cluster_atom + 1) = atom_value_offsets(i_cluster_atom) + &
                                                  nactive*atom_n_ao(i_cluster_atom)
      END DO

      ALLOCATE (ao_point_indices(atom_point_offsets(SIZE(cluster_atoms) + 1) - 1), &
                ao_values(atom_value_offsets(SIZE(cluster_atoms) + 1) - 1))
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         cluster_iatom = cluster_atoms(i_cluster_atom)
         ikind = bs_env%ri_rs%particle_set(cluster_iatom)%atomic_kind%kind_number
         ao => bs_env%basis_set_AO(ikind)%gto_basis_set
         active(:) = SUM((absolute_points - &
                          SPREAD(bs_env%ri_rs%particle_set(cluster_iatom)%r, 2, npoints))**2, DIM=1) &
                     <= ao%kind_radius**2
         first_point = atom_point_offsets(i_cluster_atom)
         last_point = atom_point_offsets(i_cluster_atom + 1) - 1
         first_value = atom_value_offsets(i_cluster_atom)
         last_value = atom_value_offsets(i_cluster_atom + 1) - 1
         nactive = last_point - first_point + 1
         ao_point_indices(first_point:last_point) = PACK(point_indices, active)
         ALLOCATE (atom_points(3, nactive), atom_values(nactive, ao%nsgf))
         atom_points(:, :) = RESHAPE(PACK(absolute_points, SPREAD(active, 1, 3)), [3, nactive])
         atom_values(:, :) = 0.0_dp
         CALL evaluate_ao_basis_on_points(atom_values, atom_points, ao, &
                                          bs_env%ri_rs%particle_set(cluster_iatom)%r, bs_env%ri_rs%cell)
         ao_values(first_value:last_value) = RESHAPE(atom_values, [nactive*ao%nsgf])
         DEALLOCATE (atom_points, atom_values)
      END DO
   END SUBROUTINE evaluate_cluster_ao_values

! **************************************************************************************************
!> \brief Apply the Cholesky selection of Eqs. (7)-(10) and retain selected points in V_A.
!> \param ao_point_indices Cluster-point index for every stored point row.
!> \param atom_point_offsets First stored point row for each atom, followed by the final bound.
!> \param atom_value_offsets First AO value for each atom, followed by the final bound.
!> \param atom_n_ao Number of AO functions for each atom.
!> \param ao_values Compact atom-local values ϕ_μ(r_l).
!> \param npoints Number of points in the cluster Lebedev grids.
!> \param inside_voronoi True for points inside the molecular Voronoi cell V_A.
!> \param selected_indices Selected indices in V_A; zero denotes rank exhaustion.
! **************************************************************************************************
   SUBROUTINE select_cholesky_grid_points(ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                          atom_n_ao, ao_values, npoints, inside_voronoi, &
                                          selected_indices)
      INTEGER, INTENT(IN)                                :: ao_point_indices(:), &
                                                            atom_point_offsets(:), &
                                                            atom_value_offsets(:), atom_n_ao(:)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET      :: ao_values(:)
      INTEGER, INTENT(IN)                                :: npoints
      LOGICAL, INTENT(IN)                                :: inside_voronoi(npoints)
      INTEGER, INTENT(OUT)                               :: selected_indices(:)

      INTEGER :: capacity, first_point, first_value, iblock, igrid_point, irow, last_point, &
         last_value, n_ao, new_capacity, nretained, nrows, nselected, selected_point
      REAL(KIND=dp)                                      :: selected_diagonal
      REAL(KIND=dp), ALLOCATABLE :: diagonal(:), factor(:, :), grown(:, :), initial_diagonal(:), &
         selected_ao_values(:), selected_factor_values(:), selection_matrix_column(:), values(:)
      REAL(KIND=dp), POINTER                             :: atom_ao_values(:, :)

      ALLOCATE (diagonal(npoints), initial_diagonal(npoints), &
                selection_matrix_column(npoints), values(npoints))

      ! d_l^(0) = D_ll = [Σ_μ ϕ_μ(r_l)^2]^2.
      diagonal(:) = 0.0_dp
      DO iblock = 1, SIZE(atom_n_ao)
         first_point = atom_point_offsets(iblock)
         last_point = atom_point_offsets(iblock + 1) - 1
         first_value = atom_value_offsets(iblock)
         last_value = atom_value_offsets(iblock + 1) - 1
         nrows = last_point - first_point + 1
         atom_ao_values(1:nrows, 1:atom_n_ao(iblock)) => ao_values(first_value:last_value)
         DO irow = 1, nrows
            igrid_point = ao_point_indices(first_point + irow - 1)
            diagonal(igrid_point) = diagonal(igrid_point) + SUM(atom_ao_values(irow, :)**2)
         END DO
      END DO
      diagonal(:) = diagonal**2
      initial_diagonal(:) = diagonal

      capacity = MIN(32, npoints)
      ALLOCATE (factor(npoints, capacity), selected_ao_values(MAXVAL(atom_n_ao)), &
                selected_factor_values(npoints))
      selected_indices(:) = 0
      nselected = 0
      nretained = 0
      DO WHILE (nretained < SIZE(selected_indices) .AND. nselected < npoints)
         selected_point = MAXLOC(diagonal, DIM=1)
         selected_diagonal = diagonal(selected_point)
         IF (selected_diagonal <= 0.0_dp) EXIT

         IF (nselected == capacity) THEN
            new_capacity = MIN(npoints, capacity + MAX(32, capacity/2))
            ALLOCATE (grown(npoints, new_capacity))
            grown(:, :capacity) = factor(:, :capacity)
            CALL MOVE_ALLOC(grown, factor)
            capacity = new_capacity
         END IF

         ! Form D_lq from atom-local AO products without storing the dense D_ll' matrix.
         selection_matrix_column(:) = 0.0_dp
         DO iblock = 1, SIZE(atom_n_ao)
            first_point = atom_point_offsets(iblock)
            last_point = atom_point_offsets(iblock + 1) - 1
            nrows = last_point - first_point + 1
            IF (nrows == 0) CYCLE
            irow = locate(ao_point_indices(first_point:last_point), selected_point)
            IF (irow == 0) CYCLE
            first_value = atom_value_offsets(iblock)
            last_value = atom_value_offsets(iblock + 1) - 1
            n_ao = atom_n_ao(iblock)
            atom_ao_values(1:nrows, 1:n_ao) => ao_values(first_value:last_value)
            selected_ao_values(:n_ao) = atom_ao_values(irow, :)
            CALL dgemv('N', nrows, n_ao, 1.0_dp, atom_ao_values, nrows, &
                       selected_ao_values, 1, 0.0_dp, values, 1)
            DO irow = 1, nrows
               igrid_point = ao_point_indices(first_point + irow - 1)
               selection_matrix_column(igrid_point) = &
                  selection_matrix_column(igrid_point) + values(irow)
            END DO
         END DO
         selection_matrix_column(:) = selection_matrix_column**2

         ! L_lk = [D_lq_k - sum_(j<k) L_lj L_q_kj]/sqrt(d_q_k^(k-1)).
         IF (nselected > 0) THEN
            selected_factor_values(:nselected) = factor(selected_point, :nselected)
            CALL dgemv('N', npoints, nselected, -1.0_dp, factor, npoints, &
                       selected_factor_values, 1, 1.0_dp, selection_matrix_column, 1)
         END IF
         nselected = nselected + 1
         factor(:, nselected) = selection_matrix_column/SQRT(selected_diagonal)
         diagonal(:) = MAX(0.0_dp, diagonal - factor(:, nselected)**2)
         WHERE (diagonal <= 64.0_dp*EPSILON(1.0_dp)*initial_diagonal) diagonal = 0.0_dp
         diagonal(selected_point) = 0.0_dp

         IF (inside_voronoi(selected_point)) THEN
            nretained = nretained + 1
            selected_indices(nretained) = selected_point
         END IF
      END DO
   END SUBROUTINE select_cholesky_grid_points

! **************************************************************************************************
!> \brief Communicate the independently constructed atom grids after Eqs. (5)-(10).
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE broadcast_ri_rs_grids(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: atom_mpi_rank, iatom

      DO iatom = 1, bs_env%n_atom
         atom_mpi_rank = MOD(iatom - 1, bs_env%para_env%num_pe)
         CALL bs_env%para_env%bcast(bs_env%ri_rs%atomic_grids(iatom)%raw_points, atom_mpi_rank)
      END DO
   END SUBROUTINE broadcast_ri_rs_grids

END MODULE gw_ri_rs_grid_initialization
