!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared RI-RS grid I/O and contracted-Gaussian evaluation utilities.
! **************************************************************************************************
MODULE rirs_grid_utils
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rirs_grid_utils'

   PUBLIC :: evaluate_ao_basis_on_points, evaluate_ao_on_points, initialize_rirs_grid, &
             read_rirs_grid_file

CONTAINS

! **************************************************************************************************
!> \brief Read Cartesian grid points from the existing CP2K RI-RS .ion format.
!> \param filename Complete input filename; no suffix or element-name construction is performed.
!> \param points Grid points in Bohr, indexed (alpha,l).
! **************************************************************************************************
   SUBROUTINE read_rirs_grid_file(filename, points)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: points

      CHARACTER(len=*), PARAMETER :: routineN = 'read_rirs_grid_file'

      CHARACTER(LEN=default_path_length)                 :: line
      INTEGER                                            :: colon, handle, ierr, iunit, l, npoints
      LOGICAL                                            :: found_points, found_size

      CALL timeset(routineN, handle)
      CALL open_file(file_name=TRIM(filename), unit_number=iunit, &
                     file_action='READ', file_status='OLD')

      found_size = .FALSE.
      npoints = 0
      DO
         READ (iunit, '(A)', IOSTAT=ierr) line
         IF (ierr /= 0) EXIT
         IF (INDEX(line, 'n points') > 0) THEN
            colon = INDEX(line, ':')
            IF (colon > 0) THEN
               READ (line(colon + 1:), *, IOSTAT=ierr) npoints
               found_size = ierr == 0 .AND. npoints > 0
            END IF
            EXIT
         END IF
      END DO
      IF (.NOT. found_size) CPABORT('RI-RS .ion file has no valid n points field: '//TRIM(filename))

      REWIND (iunit)
      found_points = .FALSE.
      DO
         READ (iunit, '(A)', IOSTAT=ierr) line
         IF (ierr /= 0) EXIT
         IF (INDEX(line, '<grid_points>') > 0) THEN
            found_points = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT. found_points) CPABORT('RI-RS .ion file has no <grid_points> block: '//TRIM(filename))

      ALLOCATE (points(3, npoints))
      DO l = 1, npoints
         READ (iunit, *, IOSTAT=ierr) points(:, l)
         IF (ierr /= 0) CPABORT('Invalid grid point in RI-RS .ion file: '//TRIM(filename))
      END DO
      CALL close_file(unit_number=iunit)
      CALL timestop(handle)
   END SUBROUTINE read_rirs_grid_file

! **************************************************************************************************
!> \brief Build one deterministic atom-specific grid from a tabulated source grid.
!> \param points Source-grid coordinates on entry and selected coordinates on return.
!> \param n_select Number of requested points.
!> \param center_atom Atom on which the relative source grid is centred.
!> \param particle_set Molecular atom positions.
!> \param n_voronoi_candidates Number of candidates available in the atom's Voronoi volume.
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE initialize_rirs_grid(points, n_select, center_atom, particle_set, n_voronoi_candidates)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: points(:, :)
      INTEGER, INTENT(IN)                                :: n_select, center_atom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(OUT)                               :: n_voronoi_candidates

      REAL(KIND=dp), ALLOCATABLE                         :: source_points(:, :)

      CPASSERT(SIZE(points, 1) == 3)
      CPASSERT(center_atom >= 1 .AND. center_atom <= SIZE(particle_set))
      ALLOCATE (source_points, source=points)
      CALL filter_grid_to_voronoi(points, center_atom, particle_set)
      n_voronoi_candidates = SIZE(points, 2)
      IF (n_voronoi_candidates < n_select) CALL MOVE_ALLOC(source_points, points)
      CALL select_grid_points(points, n_select)
   END SUBROUTINE initialize_rirs_grid

! **************************************************************************************************
!> \brief Retain source-grid points inside the Voronoi volume of one atom.
!> \param points Source-grid coordinates on entry and Voronoi-filtered coordinates on return.
!> \param center_atom Atom whose Voronoi volume is retained.
!> \param particle_set Molecular atom positions.
! **************************************************************************************************
   SUBROUTINE filter_grid_to_voronoi(points, center_atom, particle_set)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: points(:, :)
      INTEGER, INTENT(IN)                                :: center_atom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      INTEGER                                            :: iatom, ipoint, n_keep
      LOGICAL, ALLOCATABLE                               :: keep(:)
      REAL(KIND=dp)                                      :: other_distance_squared, &
                                                            own_distance_squared, physical_point(3)
      REAL(KIND=dp), ALLOCATABLE                         :: filtered_points(:, :)

      ALLOCATE (keep(SIZE(points, 2)), source=.TRUE.)
      DO ipoint = 1, SIZE(points, 2)
         physical_point = particle_set(center_atom)%r + points(:, ipoint)
         own_distance_squared = SUM(points(:, ipoint)**2)
         DO iatom = 1, SIZE(particle_set)
            IF (iatom == center_atom) CYCLE
            other_distance_squared = SUM((physical_point - particle_set(iatom)%r)**2)
            IF (other_distance_squared < own_distance_squared) THEN
               keep(ipoint) = .FALSE.
               EXIT
            END IF
         END DO
      END DO

      n_keep = COUNT(keep)
      ALLOCATE (filtered_points(3, n_keep))
      IF (n_keep > 0) filtered_points(:, :) = RESHAPE(PACK(points, SPREAD(keep, 1, 3)), [3, n_keep])
      CALL MOVE_ALLOC(filtered_points, points)
   END SUBROUTINE filter_grid_to_voronoi

! **************************************************************************************************
!> \brief Select a deterministic maximin subset from a source grid.
!> \param points Candidate coordinates on entry and selected coordinates on return.
!> \param n_select Number of points to select.
! **************************************************************************************************
   SUBROUTINE select_grid_points(points, n_select)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: points(:, :)
      INTEGER, INTENT(IN)                                :: n_select

      INTEGER                                            :: candidate, i, isel, n_source
      LOGICAL, ALLOCATABLE                               :: available(:)
      REAL(KIND=dp)                                      :: best_distance, distance_squared
      REAL(KIND=dp), ALLOCATABLE                         :: nearest_squared(:), selected_points(:, :)

      n_source = SIZE(points, 2)
      IF (n_select <= 0) CPABORT("GRID_SIZE point counts must be positive.")
      IF (n_select > n_source) THEN
         CPABORT("GRID_SIZE exceeds the number of available RI-RS source-grid points.")
      END IF
      IF (n_select == n_source) RETURN

      ALLOCATE (available(n_source), nearest_squared(n_source), selected_points(3, n_select))
      available = .TRUE.

      candidate = MINLOC(SUM(points**2, DIM=1), DIM=1)
      selected_points(:, 1) = points(:, candidate)
      available(candidate) = .FALSE.
      nearest_squared(:) = SUM((points - SPREAD(selected_points(:, 1), 2, n_source))**2, DIM=1)

      DO isel = 2, n_select
         candidate = 0
         best_distance = -1.0_dp
         DO i = 1, n_source
            IF (available(i) .AND. nearest_squared(i) > best_distance) THEN
               candidate = i
               best_distance = nearest_squared(i)
            END IF
         END DO
         CPASSERT(candidate > 0)
         selected_points(:, isel) = points(:, candidate)
         available(candidate) = .FALSE.
         DO i = 1, n_source
            IF (available(i)) THEN
               distance_squared = SUM((points(:, i) - selected_points(:, isel))**2)
               nearest_squared(i) = MIN(nearest_squared(i), distance_squared)
            END IF
         END DO
      END DO

      CALL MOVE_ALLOC(selected_points, points)
   END SUBROUTINE select_grid_points

! **************************************************************************************************
!> \brief Evaluate contracted spherical AOs, and optionally their grid-coordinate derivatives.
!>
!> For atom A and displacement d(alpha)=r_l(alpha)-R_A(alpha), the value is
!>
!>   phi_mu(r_l) = sum_(p,c) S_(pc,mu) d_x^lx d_y^ly d_z^lz exp(-zeta_p |d|^2).
!>
!> The optional derivative is evaluated analytically as
!>
!>   d phi_mu(r_l)/d r_(l,alpha)
!>     = sum_(p,c) S_(pc,mu) exp(-zeta_p |d|^2)
!>       [d polynomial_c/d d_alpha - 2 zeta_p d_alpha polynomial_c].
!>
!> The routine is shared by production GW RI-RS and grid optimization.  Production callers can
!> pass a cutoff_squared; the optimizer omits it so moving points never cross a discontinuous
!> AO-screening boundary.
!>
!> \param phi AO values, accumulated into phi(l,mu).
!> \param grid_points Cartesian grid points, indexed (alpha,l).
!> \param atom_index Source atom whose contracted AOs are evaluated.
!> \param particle_set Molecular particles.
!> \param qs_kind_set Quickstep atomic kinds containing the ORB bases.
!> \param cell Simulation cell used for the minimum-image displacement.
!> \param dphi Optional AO derivatives, accumulated into dphi(alpha,l,mu).
!> \param cutoff_squared Optional squared AO cutoff radius.
! **************************************************************************************************
   SUBROUTINE evaluate_ao_on_points(phi, grid_points, atom_index, particle_set, qs_kind_set, cell, &
                                    dphi, cutoff_squared)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: phi
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, INTENT(IN)                                :: atom_index
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dphi
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff_squared

      CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_ao_on_points'

      INTEGER                                            :: handle, kind_index
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(grid_points, 1) == 3)
      CPASSERT(SIZE(phi, 1) == SIZE(grid_points, 2))
      IF (PRESENT(dphi)) THEN
         CPASSERT(SIZE(dphi, 1) == 3)
         CPASSERT(SIZE(dphi, 2) == SIZE(phi, 1))
         CPASSERT(SIZE(dphi, 3) == SIZE(phi, 2))
      END IF

      kind_index = particle_set(atom_index)%atomic_kind%kind_number
      CALL get_qs_kind(qs_kind_set(kind_index), basis_set=basis, basis_type='ORB')
      IF (.NOT. ASSOCIATED(basis)) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL evaluate_ao_basis_on_points(phi, grid_points, basis, &
                                       particle_set(atom_index)%r, cell, dphi, cutoff_squared)
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

      INTEGER                                            :: alpha, first_sgf, ico, iend_co, ipgf, &
                                                            iset, isgf, ishell, istart_co, l, &
                                                            last_sgf, lx, ly, lz, n_cart_total, &
                                                            point, row_index
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
      !$OMP PRIVATE(point, relative, radius2, iset, n_cart_total, ishell, l, istart_co, &
      !$OMP         iend_co, first_sgf, last_sgf, ipgf, exponent, exponential, isgf, ico, &
      !$OMP         row_index, weight, lx, ly, lz, polynomial, polynomial_derivative, alpha)
      DO point = 1, SIZE(grid_points, 2)
         relative = pbc(grid_points(:, point) - source_position, cell)
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
                        row_index = (ipgf - 1)*n_cart_total + ico
                        weight = basis%sphi(row_index, isgf)
                        lx = indco(1, ico)
                        ly = indco(2, ico)
                        lz = indco(3, ico)
                        polynomial = relative(1)**lx*relative(2)**ly*relative(3)**lz
                        phi(point, isgf) = phi(point, isgf) + weight*polynomial*exponential

                        IF (PRESENT(dphi)) THEN
                           polynomial_derivative = 0.0_dp
                           IF (lx > 0) polynomial_derivative(1) = REAL(lx, dp)*relative(1)**(lx - 1)* &
                                                                  relative(2)**ly*relative(3)**lz
                           IF (ly > 0) polynomial_derivative(2) = REAL(ly, dp)*relative(1)**lx* &
                                                                  relative(2)**(ly - 1)*relative(3)**lz
                           IF (lz > 0) polynomial_derivative(3) = REAL(lz, dp)*relative(1)**lx* &
                                                                  relative(2)**ly*relative(3)**(lz - 1)
                           DO alpha = 1, 3
                              dphi(alpha, point, isgf) = dphi(alpha, point, isgf) + weight*exponential* &
                                                         (polynomial_derivative(alpha) - 2.0_dp*exponent*relative(alpha)*polynomial)
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

END MODULE rirs_grid_utils
