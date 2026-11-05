!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Geometry and AO evaluation for real-space resolution of the identity.
! **************************************************************************************************
MODULE ri_rs_utils
   USE atomic_kind_types,               ONLY: get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc,&
                                              pbc_stable
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE ri_rs_types,                     ONLY: ri_rs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: evaluate_ao_basis_on_points, evaluate_ao_on_points, &
             get_rirs_cluster_atoms, precompute_ri_rs_radii, &
             filter_grid_to_voronoi

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_utils'

   ! Automatic cluster radius in units of the nearest interatomic distance.
   REAL(KIND=dp), PARAMETER :: cluster_radius_per_nearest = 1.4_dp

CONTAINS

! **************************************************************************************************
!> \brief Compute per-atom AO and RI basis radii from the most diffuse Gaussian
!>        primitive in the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets.
!>        Stores per-atom screening radii in ri_rs and prints a per-kind table.
!>        Radius: r_kind = sqrt(-log(eps)/alpha_min_kind), with eps = eps_filter.
!> \param ri_rs RI-RS geometry references and screening radii
!> \param basis_set_AO AO basis sets, one per atomic kind
!> \param basis_set_RI auxiliary basis sets, one per atomic kind
!> \param eps screening threshold
!> \param unit_nr output unit
! **************************************************************************************************
   SUBROUTINE precompute_ri_rs_radii(ri_rs, basis_set_AO, basis_set_RI, eps, unit_nr)
      TYPE(ri_rs_env), INTENT(INOUT)                     :: ri_rs
      TYPE(gto_basis_set_p_type), INTENT(IN)             :: basis_set_AO(:), basis_set_RI(:)
      REAL(KIND=dp), INTENT(IN)                          :: eps
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precompute_ri_rs_radii'
      REAL(KIND=dp), PARAMETER :: min_exponent_for_radius = 1.0E-3_dp

      INTEGER                                            :: handle, i, iatom, ikind, j, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha_min_ao_kind, alpha_min_ri_kind
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet_ao, zet_ri

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ri_rs%atomic_kind_set))
      CPASSERT(ASSOCIATED(ri_rs%particle_set))
      CPASSERT(0.0_dp < eps .AND. eps < 1.0_dp)
      nkind = SIZE(ri_rs%atomic_kind_set)
      natom = SIZE(ri_rs%particle_set)

      ALLOCATE (alpha_min_ao_kind(nkind), alpha_min_ri_kind(nkind))
      alpha_min_ao_kind = HUGE(1.0_dp)
      alpha_min_ri_kind = HUGE(1.0_dp)

      DO ikind = 1, nkind
         zet_ao => basis_set_AO(ikind)%gto_basis_set%zet
         zet_ri => basis_set_RI(ikind)%gto_basis_set%zet
         DO i = 1, SIZE(zet_ao, 1)
            DO j = 1, SIZE(zet_ao, 2)
               IF (zet_ao(i, j) > min_exponent_for_radius) THEN
                  alpha_min_ao_kind(ikind) = MIN(alpha_min_ao_kind(ikind), zet_ao(i, j))
               END IF
            END DO
         END DO
         DO i = 1, SIZE(zet_ri, 1)
            DO j = 1, SIZE(zet_ri, 2)
               IF (zet_ri(i, j) > min_exponent_for_radius) THEN
                  alpha_min_ri_kind(ikind) = MIN(alpha_min_ri_kind(ikind), zet_ri(i, j))
               END IF
            END DO
         END DO
      END DO

      CALL get_atomic_kind_set(atomic_kind_set=ri_rs%atomic_kind_set, kind_of=kind_of)

      ALLOCATE (ri_rs%radius_ao_per_atom(natom))
      ALLOCATE (ri_rs%radius_ri_per_atom(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         ri_rs%radius_ao_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))
         ri_rs%radius_ri_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))
      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'RI-RS basis radii (Å):'
         WRITE (unit_nr, '(T4,A6,2X,A4,2A14)') 'Kind', 'Elem', 'r_AO (Å)', 'r_RI (Å)'
         DO ikind = 1, nkind
            WRITE (unit_nr, '(T4,I6,2X,A4,2F14.4)') &
               ikind, &
               ri_rs%atomic_kind_set(ikind)%element_symbol, &
               SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))*angstrom, &
               SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))*angstrom
         END DO
         WRITE (unit_nr, '(A)') ' '
      END IF

      DEALLOCATE (alpha_min_ao_kind, alpha_min_ri_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE precompute_ri_rs_radii

! **************************************************************************************************
!> \brief Keep atom-relative points inside the Voronoi volume of one atom.
!> \param points Atom-relative source-grid coordinates; unchanged when mask is present.
!> \param icenter_atom Atom in the reference cell whose Voronoi volume is retained.
!> \param particle_set Atoms in the reference cell.
!> \param cell Simulation cell and periodicity.
!> \param mask Optional membership mask.
! **************************************************************************************************
   SUBROUTINE filter_grid_to_voronoi(points, icenter_atom, particle_set, cell, mask)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: points(:, :)
      INTEGER, INTENT(IN)                                :: icenter_atom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(OUT), OPTIONAL                     :: mask(:)

      INTEGER                                            :: iimage, ipoint, n_keep
      INTEGER, ALLOCATABLE                               :: atom_indices(:), cell_indices(:, :)
      LOGICAL, ALLOCATABLE                               :: keep(:)
      REAL(KIND=dp)                                      :: displacement(3), distance2, max_radius, &
                                                            tolerance
      REAL(KIND=dp), ALLOCATABLE                         :: centers(:, :), filtered_points(:, :)

      CPASSERT(SIZE(points, 1) == 3)
      ALLOCATE (keep(SIZE(points, 2)), SOURCE=.TRUE.)
      IF (SIZE(points, 2) > 0) THEN
         max_radius = MAXVAL(NORM2(points, DIM=1))
      ELSE
         max_radius = 0.0_dp
      END IF
      CALL get_rirs_cluster_atoms(particle_set, cell, icenter_atom, 2.0_dp*max_radius, &
                                  atom_indices, cell_indices, centers)
      DO iimage = 1, SIZE(atom_indices)
         IF (atom_indices(iimage) == icenter_atom .AND. &
             ALL(cell_indices(:, iimage) == 0)) CYCLE
         displacement(:) = centers(:, iimage) - pbc(particle_set(icenter_atom)%r, cell)
         distance2 = SUM(displacement**2)
         tolerance = 32.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, distance2)
         DO ipoint = 1, SIZE(points, 2)
            IF (.NOT. keep(ipoint)) CYCLE
            ! |r_l-R_B|² - |r_l-R_A|² = |R_B-R_A|² - 2(r_l-R_A)·(R_B-R_A).
            IF (2.0_dp*DOT_PRODUCT(points(:, ipoint), displacement) > distance2 + tolerance) THEN
               keep(ipoint) = .FALSE.
            ELSE IF (atom_indices(iimage) < icenter_atom) THEN
               IF (ABS(2.0_dp*DOT_PRODUCT(points(:, ipoint), displacement) - distance2) <= &
                   tolerance) keep(ipoint) = .FALSE.
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
!> \brief Form a cluster from all atom images inside a sphere about a reference atom.
!>
!> Nonperiodic directions have only the zero translation.
!> Every image is identified by its reference-cell atom and integer cell translation. Thus, the
!> same atom can occur more than once in a small periodic cell. The central atom is returned first;
!> the remaining images follow in atom-major, relative-image-minor order.
!> \param particle_set Atoms in the reference cell.
!> \param cell Simulation cell and periodicity.
!> \param icenter_atom Central atom in cell (0,0,0).
!> \param radius Cluster radius; a negative value selects 1.4 times the nearest-neighbor distance.
!> \param atom_indices Reference-cell atom index for every image.
!> \param cell_indices Integer cell translation for every image.
!> \param centers Cartesian image centres.
!> \param effective_radius Resolved nonnegative radius, in internal length units.
! **************************************************************************************************
   SUBROUTINE get_rirs_cluster_atoms(particle_set, cell, icenter_atom, radius, atom_indices, &
                                     cell_indices, centers, effective_radius)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: icenter_atom
      REAL(KIND=dp), INTENT(IN)                          :: radius
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: atom_indices(:)
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: cell_indices(:, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT), OPTIONAL  :: centers(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: effective_radius

      INTEGER                                            :: iatom, idim, image(3), image_max(3), &
                                                            image_min(3), image_shift(3), ix, iy, &
                                                            iz, max_images, nimages
      INTEGER, ALLOCATABLE                               :: atom_work(:), cell_work(:, :)
      REAL(KIND=dp)                                      :: center(3), center_atom(3), &
                                                            cluster_radius, delta(3), radius2, &
                                                            tolerance
      REAL(KIND=dp), ALLOCATABLE                         :: center_work(:, :)

      cluster_radius = radius
      IF (radius < 0.0_dp) cluster_radius = automatic_cluster_radius(particle_set, cell, icenter_atom)
      IF (PRESENT(effective_radius)) effective_radius = cluster_radius
      center_atom(:) = pbc(particle_set(icenter_atom)%r, cell)
      DO idim = 1, 3
         IF (cell%perd(idim) == 1) THEN
            image_max(idim) = CEILING(cluster_radius*NORM2(cell%h_inv(idim, :))) + 1
         ELSE
            image_max(idim) = 0
         END IF
      END DO
      image_min(:) = -image_max(:)
      max_images = SIZE(particle_set)*PRODUCT(image_max - image_min + 1)
      ALLOCATE (atom_work(max_images), cell_work(3, max_images), center_work(3, max_images))
      radius2 = cluster_radius**2
      tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, radius2)
      nimages = 1
      atom_work(1) = icenter_atom
      cell_work(:, 1) = 0
      center_work(:, 1) = center_atom
      DO iatom = 1, SIZE(particle_set)
         delta(:) = pbc_stable(particle_set(iatom)%r - particle_set(icenter_atom)%r, cell)
         ! δ_i = r_i^w - r_c^w + H s_i; cell labels remain those of the reference-cell atoms.
         image_shift(:) = NINT(MATMUL(cell%h_inv, delta - (pbc(particle_set(iatom)%r, cell) - center_atom)))
         DO iz = image_min(3), image_max(3)
            DO iy = image_min(2), image_max(2)
               DO ix = image_min(1), image_max(1)
                  image(:) = [ix, iy, iz]
                  IF (iatom == icenter_atom .AND. ALL(image == 0)) CYCLE
                  center(:) = delta + MATMUL(cell%hmat, REAL(image, dp))
                  IF (SUM(center**2) > radius2 + tolerance) CYCLE
                  nimages = nimages + 1
                  atom_work(nimages) = iatom
                  cell_work(:, nimages) = image + image_shift
                  center_work(:, nimages) = center_atom + center
               END DO
            END DO
         END DO
      END DO
      ALLOCATE (atom_indices(nimages), SOURCE=atom_work(:nimages))
      IF (PRESENT(cell_indices)) ALLOCATE (cell_indices(3, nimages), SOURCE=cell_work(:, :nimages))
      IF (PRESENT(centers)) ALLOCATE (centers(3, nimages), SOURCE=center_work(:, :nimages))
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
!> \param cell Simulation cell defining minimum-image or image-summed displacements.
!> \param dphi Optional AO derivatives, accumulated into dphi(alpha,l,mu).
!> \param cutoff_squared Optional squared AO cutoff radius.
!> \param periodic_images Sum images -1:1 along periodic axes and only image 0 along other axes.
!> \param minimum_image Wrap displacements into the cell (default true); disable for explicit images.
! **************************************************************************************************
   SUBROUTINE evaluate_ao_basis_on_points(phi, grid_points, basis, source_position, cell, &
                                          dphi, cutoff_squared, periodic_images, minimum_image)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: phi
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: source_position
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dphi
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff_squared
      LOGICAL, INTENT(IN), OPTIONAL                      :: periodic_images, minimum_image

      INTEGER :: first_sgf, ialpha, ico, iend_co, image_extent(3), ipgf, ipoint, irow, iset, isgf, &
         ishell, istart_co, ix, iy, iz, l, last_sgf, lx, ly, lz, n_cart_total
      LOGICAL                                            :: sum_images, use_minimum_image
      REAL(KIND=dp) :: exponent, exponential, hmat(3, 3), polynomial, polynomial_derivative(3), &
         radius2, relative(3), relative_raw(3), weight

      CPASSERT(ASSOCIATED(basis))
      CPASSERT(SIZE(grid_points, 1) == 3)
      CPASSERT(SIZE(phi, 1) == SIZE(grid_points, 2))
      IF (PRESENT(dphi)) THEN
         CPASSERT(SIZE(dphi, 1) == 3)
         CPASSERT(SIZE(dphi, 2) == SIZE(phi, 1))
         CPASSERT(SIZE(dphi, 3) == SIZE(phi, 2))
      END IF

      use_minimum_image = .TRUE.
      IF (PRESENT(minimum_image)) use_minimum_image = minimum_image
      sum_images = .FALSE.
      IF (PRESENT(periodic_images)) sum_images = periodic_images
      image_extent = 0
      hmat = 0.0_dp
      IF (sum_images) THEN
         image_extent = cell%perd
         CALL get_cell(cell, h=hmat)
      END IF

      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(phi, dphi, grid_points, cell, basis, source_position, cutoff_squared, ncoset, indco) &
      !$OMP SHARED(sum_images, image_extent, hmat, use_minimum_image) &
      !$OMP PRIVATE(ipoint, relative, radius2, iset, n_cart_total, ishell, l, istart_co, &
      !$OMP         iend_co, first_sgf, last_sgf, ipgf, exponent, exponential, isgf, ico, &
      !$OMP         irow, weight, lx, ly, lz, polynomial, polynomial_derivative, ialpha, &
      !$OMP         relative_raw, ix, iy, iz)
      DO ipoint = 1, SIZE(grid_points, 2)
         relative_raw = grid_points(:, ipoint) - source_position
         IF (.NOT. sum_images .AND. use_minimum_image) relative_raw = pbc(relative_raw, cell)
         DO ix = -image_extent(1), image_extent(1)
         DO iy = -image_extent(2), image_extent(2)
         DO iz = -image_extent(3), image_extent(3)
            relative = relative_raw
            IF (sum_images) relative = relative_raw - MATMUL(hmat, REAL([ix, iy, iz], dp))
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
         END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE evaluate_ao_basis_on_points

! **************************************************************************************************
!> \brief Resolve radius_A = 1.4 min_(B,T)!=(A,0) |R_B + T - R_A|, or zero for an isolated atom.
!>
!> A wrapped interatomic displacement or a periodic cell vector supplies an upper bound on the
!> nearest distance. Reciprocal-cell bounds enumerate every image inside that sphere, including
!> for skewed cells. Only the central image is excluded; coincident distinct atoms are invalid.
!> \param particle_set ...
!> \param cell ...
!> \param icenter_atom ...
!> \return ...
! **************************************************************************************************
   FUNCTION automatic_cluster_radius(particle_set, cell, icenter_atom) RESULT(radius)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: icenter_atom
      REAL(KIND=dp)                                      :: radius

      INTEGER                                            :: iatom, idim, image_max(3), ix, iy, iz
      REAL(KIND=dp)                                      :: delta(3), distance2, nearest2, upper2

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(icenter_atom >= 1 .AND. icenter_atom <= SIZE(particle_set))
      IF (SIZE(particle_set) == 1 .AND. ALL(cell%perd == 0)) THEN
         radius = 0.0_dp
         RETURN
      END IF
      upper2 = HUGE(1.0_dp)
      DO idim = 1, 3
         IF (cell%perd(idim) == 1) upper2 = MIN(upper2, SUM(cell%hmat(:, idim)**2))
      END DO
      DO iatom = 1, SIZE(particle_set)
         IF (iatom == icenter_atom) CYCLE
         delta = pbc_stable(particle_set(iatom)%r - particle_set(icenter_atom)%r, cell)
         upper2 = MIN(upper2, SUM(delta**2))
      END DO
      DO idim = 1, 3
         image_max(idim) = 0
         IF (cell%perd(idim) == 1) THEN
            image_max(idim) = CEILING(SQRT(upper2)*NORM2(cell%h_inv(idim, :))) + 1
         END IF
      END DO
      nearest2 = upper2
      DO iatom = 1, SIZE(particle_set)
         delta = pbc_stable(particle_set(iatom)%r - particle_set(icenter_atom)%r, cell)
         DO iz = -image_max(3), image_max(3)
            DO iy = -image_max(2), image_max(2)
               DO ix = -image_max(1), image_max(1)
                  IF (iatom == icenter_atom .AND. ix == 0 .AND. iy == 0 .AND. iz == 0) CYCLE
                  distance2 = SUM((delta + MATMUL(cell%hmat, REAL([ix, iy, iz], dp)))**2)
                  nearest2 = MIN(nearest2, distance2)
               END DO
            END DO
         END DO
      END DO
      IF (nearest2 <= 0.0_dp) CPABORT("Automatic RI-RS cluster radius requires distinct atomic positions")
      radius = cluster_radius_per_nearest*SQRT(nearest2)
   END FUNCTION automatic_cluster_radius

END MODULE ri_rs_utils
