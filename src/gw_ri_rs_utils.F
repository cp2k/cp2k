!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common setup operations used by the periodic and non-periodic GW RI-RS implementations.
! **************************************************************************************************
MODULE gw_ri_rs_utils
   USE atomic_kind_types,               ONLY: get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: evaluate_ao_basis_on_points, evaluate_ao_on_points, filter_grid_to_voronoi, &
             get_rirs_cluster_atoms, precompute_ri_rs_radii

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_utils'

CONTAINS

! **************************************************************************************************
!> \brief Compute per-atom AO and RI basis radii from the most diffuse Gaussian
!>        primitive in the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets.
!>        Stores results in bs_env%ri_rs%radius_ao_per_atom(:) and
!>        bs_env%ri_rs%radius_ri_per_atom(:) and prints a per-atom table.
!>        Radius: r_kind = sqrt(-log(eps)/alpha_min_kind), with eps = eps_filter.
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************
   SUBROUTINE precompute_ri_rs_radii(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precompute_ri_rs_radii'

      INTEGER                                            :: handle, i, iatom, ikind, j, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha_min_ao_kind, alpha_min_ri_kind
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet_ao, zet_ri

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(bs_env%ri_rs%atomic_kind_set))
      nkind = SIZE(bs_env%ri_rs%atomic_kind_set)
      natom = bs_env%n_atom
      eps = bs_env%eps_filter

      ALLOCATE (alpha_min_ao_kind(nkind), alpha_min_ri_kind(nkind))
      alpha_min_ao_kind = HUGE(1.0_dp)
      alpha_min_ri_kind = HUGE(1.0_dp)

      DO ikind = 1, nkind
         zet_ao => bs_env%basis_set_AO(ikind)%gto_basis_set%zet
         zet_ri => bs_env%basis_set_RI(ikind)%gto_basis_set%zet
         DO i = 1, SIZE(zet_ao, 1)
            DO j = 1, SIZE(zet_ao, 2)
               IF (zet_ao(i, j) > 1.0E-3_dp) THEN
                  alpha_min_ao_kind(ikind) = MIN(alpha_min_ao_kind(ikind), zet_ao(i, j))
               END IF
            END DO
         END DO
         DO i = 1, SIZE(zet_ri, 1)
            DO j = 1, SIZE(zet_ri, 2)
               IF (zet_ri(i, j) > 1.0E-3_dp) THEN
                  alpha_min_ri_kind(ikind) = MIN(alpha_min_ri_kind(ikind), zet_ri(i, j))
               END IF
            END DO
         END DO
      END DO

      CALL get_atomic_kind_set(atomic_kind_set=bs_env%ri_rs%atomic_kind_set, kind_of=kind_of)

      ALLOCATE (bs_env%ri_rs%radius_ao_per_atom(natom))
      ALLOCATE (bs_env%ri_rs%radius_ri_per_atom(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         bs_env%ri_rs%radius_ao_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))
         bs_env%ri_rs%radius_ri_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))
      END DO

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS basis radii (Å):'
         WRITE (bs_env%unit_nr, '(T4,A6,2X,A4,2A14)') 'Kind', 'Elem', 'r_AO (Å)', 'r_RI (Å)'
         DO ikind = 1, nkind
            WRITE (bs_env%unit_nr, '(T4,I6,2X,A4,2F14.4)') &
               ikind, &
               bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol, &
               SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))*angstrom, &
               SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))*angstrom
         END DO
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      DEALLOCATE (alpha_min_ao_kind, alpha_min_ri_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE precompute_ri_rs_radii

!> \brief Retain source-grid points inside the Voronoi volume of one atom.
!> \param points Source-grid coordinates on entry and Voronoi-filtered coordinates on return.
!> \param icenter_atom Atom whose Voronoi volume is retained.
!> \param particle_set Molecular atom positions.
!> \param mask Optional membership mask; if present, leave points unchanged.
!> \param atom_indices Optional subset of nuclei defining the Voronoi partition.
! **************************************************************************************************
   SUBROUTINE filter_grid_to_voronoi(points, icenter_atom, particle_set, mask, atom_indices)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: points(:, :)
      INTEGER, INTENT(IN)                                :: icenter_atom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      LOGICAL, INTENT(OUT), OPTIONAL                     :: mask(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_indices(:)

      INTEGER                                            :: iatom, iatom_index, ipoint, n_atoms, &
                                                            n_keep
      LOGICAL, ALLOCATABLE                               :: keep(:)
      REAL(KIND=dp)                                      :: displacement(3), distance2, max_radius2, &
                                                            tolerance
      REAL(KIND=dp), ALLOCATABLE                         :: filtered_points(:, :)

      ALLOCATE (keep(SIZE(points, 2)), SOURCE=.TRUE.)
      max_radius2 = MAXVAL(SUM(points**2, DIM=1))
      n_atoms = SIZE(particle_set)
      IF (PRESENT(atom_indices)) n_atoms = SIZE(atom_indices)
      DO iatom_index = 1, n_atoms
         IF (PRESENT(atom_indices)) THEN
            iatom = atom_indices(iatom_index)
            CPASSERT(iatom >= 1 .AND. iatom <= SIZE(particle_set))
         ELSE
            iatom = iatom_index
         END IF
         IF (iatom == icenter_atom) CYCLE
         displacement(:) = particle_set(iatom)%r - particle_set(icenter_atom)%r
         distance2 = SUM(displacement**2)
         ! |R_B-R_A| > 2 max_l|r_l-R_A| cannot cut this finite point set.
         IF (distance2 > 4.0_dp*max_radius2) CYCLE
         tolerance = 32.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, distance2)
         DO ipoint = 1, SIZE(points, 2)
            IF (.NOT. keep(ipoint)) CYCLE
            ! |r_l-R_B|² - |r_l-R_A|² = |R_B-R_A|² - 2(r_l-R_A)·(R_B-R_A).
            IF (2.0_dp*DOT_PRODUCT(points(:, ipoint), displacement) > distance2 + tolerance) THEN
               keep(ipoint) = .FALSE.
            ELSE IF (iatom < icenter_atom) THEN
               IF (ABS(2.0_dp*DOT_PRODUCT(points(:, ipoint), displacement) - distance2) <= tolerance) THEN
                  keep(ipoint) = .FALSE.
               END IF
            END IF
         END DO
      END DO
      IF (PRESENT(mask)) THEN
         CPASSERT(SIZE(mask) == SIZE(keep))
         mask(:) = keep
         RETURN
      END IF
      n_keep = COUNT(keep)
      ALLOCATE (filtered_points(3, n_keep))
      IF (n_keep > 0) filtered_points(:, :) = RESHAPE(PACK(points, SPREAD(keep, 1, 3)), [3, n_keep])
      CALL MOVE_ALLOC(filtered_points, points)
   END SUBROUTINE filter_grid_to_voronoi

! **************************************************************************************************
!> \brief Form C_A from nuclei within the specified radius, in global atom order.
!> \param particle_set Molecular nuclei.
!> \param cell Simulation cell.
!> \param icenter_atom Central atom A.
!> \param radius Cluster radius.
!> \param atom_indices Cluster atom indices including A.
! **************************************************************************************************
   SUBROUTINE get_rirs_cluster_atoms(particle_set, cell, icenter_atom, radius, atom_indices)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: icenter_atom
      REAL(KIND=dp), INTENT(IN)                          :: radius
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: atom_indices(:)

      INTEGER                                            :: iatom, ncluster
      INTEGER, ALLOCATABLE                               :: work(:)

      ALLOCATE (work(SIZE(particle_set)))
      ncluster = 0
      DO iatom = 1, SIZE(particle_set)
         IF (SUM(pbc(particle_set(iatom)%r - particle_set(icenter_atom)%r, cell)**2) <= radius**2) THEN
            ncluster = ncluster + 1
            work(ncluster) = iatom
         END IF
      END DO
      ALLOCATE (atom_indices(ncluster), SOURCE=work(:ncluster))
   END SUBROUTINE get_rirs_cluster_atoms

! **************************************************************************************************
!> \brief Evaluate contracted spherical AOs, and optionally their grid-coordinate derivatives.
!>
!> For atom A and displacement d(alpha)=r_l(alpha)-R_A(alpha), the value is
!>
!>   ϕ_μ(r_l) = Σ_(p,c) S_(pc,μ) d_x^lx d_y^ly d_z^lz exp(-ζ_p |d|^2).
!>
!> The optional derivative is evaluated analytically as
!>
!>   d ϕ_μ(r_l)/d r_(l,α)
!>     = sum_(p,c) S_(pc,mu) exp(-zeta_p |d|^2)
!>       [d polynomial_c/d d_alpha - 2 zeta_p d_alpha polynomial_c].
!>
!> \param phi AO values, accumulated into phi(l,mu).
!> \param grid_points Cartesian grid points, indexed (alpha,l).
!> \param iatom Source atom whose contracted AOs are evaluated.
!> \param particle_set Molecular particles.
!> \param qs_kind_set Quickstep atomic kinds containing the ORB bases.
!> \param cell Simulation cell used for the minimum-image displacement.
!> \param dphi Optional AO derivatives, accumulated into dphi(alpha,l,mu).
!> \param cutoff_squared Optional squared AO cutoff radius.
! **************************************************************************************************
   SUBROUTINE evaluate_ao_on_points(phi, grid_points, iatom, particle_set, qs_kind_set, cell, &
                                    dphi, cutoff_squared)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: phi
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, INTENT(IN)                                :: iatom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dphi
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff_squared

      CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_ao_on_points'

      INTEGER                                            :: handle, ikind
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(grid_points, 1) == 3)
      CPASSERT(SIZE(phi, 1) == SIZE(grid_points, 2))
      IF (PRESENT(dphi)) THEN
         CPASSERT(SIZE(dphi, 1) == 3)
         CPASSERT(SIZE(dphi, 2) == SIZE(phi, 1))
         CPASSERT(SIZE(dphi, 3) == SIZE(phi, 2))
      END IF

      ikind = particle_set(iatom)%atomic_kind%kind_number
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type='ORB')
      IF (.NOT. ASSOCIATED(basis)) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL evaluate_ao_basis_on_points(phi, grid_points, basis, &
                                       particle_set(iatom)%r, cell, dphi, cutoff_squared)
      CALL timestop(handle)
   END SUBROUTINE evaluate_ao_on_points

! **************************************************************************************************
!> \brief Evaluate one explicitly supplied contracted Gaussian basis on Cartesian points.
!> \param phi AO values, accumulated into phi(l,mu).
!> \param grid_points Cartesian grid points, indexed (alpha,l).
!> \param basis Contracted Gaussian basis to evaluate.
!> \param source_position Centre of the basis.
!> \param cell Simulation cell used for the minimum-image displacement.
!> \param dphi Optional AO derivatives, accumulated into dphi(alpha,l,mu).
!> \param cutoff_squared Optional squared AO cutoff radius.
! **************************************************************************************************
   SUBROUTINE evaluate_ao_basis_on_points(phi, grid_points, basis, source_position, cell, &
                                          dphi, cutoff_squared)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: phi
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: source_position
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dphi
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff_squared

      INTEGER                                            :: first_sgf, ialpha, ico, iend_co, ipgf, &
                                                            ipoint, irow, iset, isgf, ishell, &
                                                            istart_co, l, last_sgf, lx, ly, lz, &
                                                            n_cart_total
      REAL(KIND=dp)                                      :: exponent, exponential, polynomial, &
                                                            polynomial_derivative(3), radius2, &
                                                            relative(3), weight

      CPASSERT(ASSOCIATED(basis))
      CPASSERT(SIZE(grid_points, 1) == 3)
      CPASSERT(SIZE(phi, 1) == SIZE(grid_points, 2))
      IF (PRESENT(dphi)) THEN
         CPASSERT(SIZE(dphi, 1) == 3)
         CPASSERT(SIZE(dphi, 2) == SIZE(phi, 1))
         CPASSERT(SIZE(dphi, 3) == SIZE(phi, 2))
      END IF

      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(phi, dphi, grid_points, cell, basis, source_position, cutoff_squared, ncoset, indco) &
      !$OMP PRIVATE(ipoint, relative, radius2, iset, n_cart_total, ishell, l, istart_co, &
      !$OMP         iend_co, first_sgf, last_sgf, ipgf, exponent, exponential, isgf, ico, &
      !$OMP         irow, weight, lx, ly, lz, polynomial, polynomial_derivative, ialpha)
      DO ipoint = 1, SIZE(grid_points, 2)
         relative = pbc(grid_points(:, ipoint) - source_position, cell)
         radius2 = DOT_PRODUCT(relative, relative)
         IF (PRESENT(cutoff_squared)) THEN
            IF (radius2 > cutoff_squared) CYCLE
         END IF

         DO iset = 1, basis%nset
            n_cart_total = ncoset(basis%lmax(iset))
            DO ishell = 1, basis%nshell(iset)
               l = basis%l(ishell, iset)
               istart_co = ncoset(l - 1) + 1
               iend_co = ncoset(l)
               first_sgf = basis%first_sgf(ishell, iset)
               last_sgf = basis%last_sgf(ishell, iset)
               DO ipgf = 1, basis%npgf(iset)
                  exponent = basis%zet(ipgf, iset)
                  exponential = EXP(-exponent*radius2)
                  DO isgf = first_sgf, last_sgf
                     DO ico = istart_co, iend_co
                        irow = (ipgf - 1)*n_cart_total + ico
                        weight = basis%sphi(irow, isgf)
                        lx = indco(1, ico)
                        ly = indco(2, ico)
                        lz = indco(3, ico)
                        polynomial = relative(1)**lx*relative(2)**ly*relative(3)**lz
                        phi(ipoint, isgf) = phi(ipoint, isgf) + weight*polynomial*exponential

                        IF (PRESENT(dphi)) THEN
                           polynomial_derivative = 0.0_dp
                           IF (lx > 0) polynomial_derivative(1) = REAL(lx, dp)*relative(1)**(lx - 1)* &
                                                                  relative(2)**ly*relative(3)**lz
                           IF (ly > 0) polynomial_derivative(2) = REAL(ly, dp)*relative(1)**lx* &
                                                                  relative(2)**(ly - 1)*relative(3)**lz
                           IF (lz > 0) polynomial_derivative(3) = REAL(lz, dp)*relative(1)**lx* &
                                                                  relative(2)**ly*relative(3)**(lz - 1)
                           DO ialpha = 1, 3
                              dphi(ialpha, ipoint, isgf) = dphi(ialpha, ipoint, isgf) + weight*exponential* &
                                                           (polynomial_derivative(ialpha) - &
                                                            2.0_dp*exponent*relative(ialpha)*polynomial)
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE evaluate_ao_basis_on_points

END MODULE gw_ri_rs_utils
