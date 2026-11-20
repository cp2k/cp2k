!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Integrate Bader basins directly from the Quickstep density.
! **************************************************************************************************
MODULE qs_bader_analysis
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_bader_basins,                 ONLY: bader_integrate_grid
   USE realspace_grid_cube,             ONLY: pw_to_cube
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: print_bader_analysis

CONTAINS

! **************************************************************************************************
!> \brief Collect the spin-summed electron density on the print rank and integrate the basins.
!> \param rho plane-wave real-space density for each spin
!> \param cell simulation cell
!> \param particles atoms in the same order as ionic_charge
!> \param ionic_charge nuclear charge or pseudopotential ionic charge in electron units
!> \param nuclear_center true for physical nuclei, false for basis-only centres
!> \param nominal_electrons electron count from physical nuclear charges and DFT charge
!> \param density_label description of the density used to define and integrate basins
!> \param cutoff density cutoff for vacuum integration
!> \param assign_nearest group maxima by proximity rather than unique nuclear seed
!> \param iw print unit, positive only on the output rank
!> \param cube_unit basin cube unit, positive only on the output rank when requested
!> \param cube_stride stride used only to write the basin cube
! **************************************************************************************************
   SUBROUTINE print_bader_analysis(rho, cell, particles, ionic_charge, nuclear_center, nominal_electrons, &
                                   density_label, cutoff, assign_nearest, iw, cube_unit, cube_stride)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho(:)
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      REAL(KIND=dp), INTENT(IN)                          :: ionic_charge(:)
      LOGICAL, INTENT(IN)                                :: nuclear_center(:)
      REAL(KIND=dp), INTENT(IN)                          :: nominal_electrons
      CHARACTER(LEN=*), INTENT(IN)                       :: density_label
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, INTENT(IN)                                :: assign_nearest
      INTEGER, INTENT(IN)                                :: iw, cube_unit
      INTEGER, POINTER                                   :: cube_stride(:)

      INTEGER, PARAMETER                                 :: max_grid_points = 8000000

      INTEGER                                            :: a, d, i, n(3), n_extra, n_other, ng, &
                                                            nnuc, q(3)
      INTEGER, ALLOCATABLE                               :: atom_indices(:), atomic_numbers(:), &
                                                            seeds(:)
      INTEGER, ALLOCATABLE, TARGET                       :: labels(:, :, :)
      INTEGER, POINTER                                   :: labels_flat(:)
      LOGICAL                                            :: periodic(3)
      LOGICAL, ALLOCATABLE                               :: assigned(:)
      REAL(KIND=dp)                                      :: extra(2), integrated, other(2), &
                                                            SPACING(3), vacuum(2)
      REAL(KIND=dp), ALLOCATABLE                         :: atom_grid(:, :), basins(:, :), &
                                                            positions(:, :)
      REAL(KIND=dp), ALLOCATABLE, TARGET                 :: density(:, :, :)
      REAL(KIND=dp), POINTER                             :: values(:)

      IF (.NOT. cell%orthorhombic) CPABORT("BADER requires a Cartesian orthorhombic cell.")
      nnuc = COUNT(nuclear_center)
      IF (nnuc == 0) CPABORT("BADER requires at least one physical nucleus.")
      periodic = cell%perd /= 0
      IF (.NOT. (cutoff >= 0.0_dp)) CPABORT("BADER requires a nonnegative VACUUM_DENSITY.")
      IF (rho(1)%pw_grid%ngpts > max_grid_points) THEN
         CPABORT("BADER supports at most 8000000 grid points in its current implementation.")
      END IF
      n = rho(1)%pw_grid%npts
      ng = PRODUCT(n)
      CALL collect_density(rho, iw, density)
      IF (iw <= 0) RETURN

      ALLOCATE (seeds(nnuc), assigned(nnuc), basins(2, nnuc))
      IF (assign_nearest) ALLOCATE (atom_grid(3, nnuc))
      spacing = [rho(1)%pw_grid%dh(1, 1), rho(1)%pw_grid%dh(2, 2), rho(1)%pw_grid%dh(3, 3)]
      seeds = 0
      a = 0
      DO i = 1, SIZE(particles)
         IF (.NOT. nuclear_center(i)) CYCLE
         a = a + 1
         IF (assign_nearest) atom_grid(:, a) = particles(i)%r/spacing
         q = NINT(particles(i)%r/spacing)
         DO d = 1, 3
            IF (periodic(d)) q(d) = MODULO(q(d), n(d))
         END DO
         IF (ANY(q < 0) .OR. ANY(q >= n)) CYCLE
         seeds(a) = 1 + q(1) + n(1)*(q(2) + n(2)*q(3))
      END DO
      values(1:ng) => density
      NULLIFY (labels_flat)
      IF (cube_unit > 0) THEN
         IF (SIZE(particles) > 99999 .OR. ANY(n > 99999)) THEN
            CPABORT("BADER basin cube exceeds the cube header limits.")
         END IF
         ALLOCATE (labels(n(1), n(2), n(3)))
         labels_flat(1:ng) => labels
      END IF
      CALL bader_integrate_grid(values, n, spacing, periodic, cutoff, seeds, basins, assigned, other, vacuum, &
                                n_other, atom_grid=atom_grid, n_extra=n_extra, extra=extra, labels=labels_flat)

      IF (cube_unit > 0) THEN
         ALLOCATE (atom_indices(nnuc), atomic_numbers(SIZE(particles)), positions(3, SIZE(particles)))
         a = 0
         DO i = 1, SIZE(particles)
            CALL get_atomic_kind(particles(i)%atomic_kind, z=atomic_numbers(i))
            positions(:, i) = particles(i)%r
            IF (.NOT. nuclear_center(i)) CYCLE
            a = a + 1
            atom_indices(a) = i
         END DO
         DO i = 1, ng
            a = labels_flat(i)
            IF (a > 0) labels_flat(i) = atom_indices(a)
         END DO
         CALL pw_to_cube(rho(1), cube_unit, title="Bader basin labels: atom index, 0 vacuum, -1 no unique owner", &
                         particles_r=positions, particles_z=atomic_numbers, particles_zeff=ionic_charge, &
                         stride=cube_stride, gathered_labels=labels)
      END IF

      WRITE (iw, '(/,T2,A)') "BADER| Weighted basins of "//density_label
      IF (assign_nearest) THEN
         WRITE (iw, '(T2,A)') "BADER| Atom assignment: nearest atom to each maximum"
      ELSE
         WRITE (iw, '(T2,A)') "BADER| Atom assignment: unique nuclear seed only"
      END IF
      WRITE (iw, '(T2,A,I10)') "BADER| Maxima without a unique nuclear seed ", n_extra
      WRITE (iw, '(T2,A,F20.10)') "BADER| Electrons without a unique nuclear seed ", extra(1)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Volume without a unique nuclear seed [bohr^3] ", extra(2)
      WRITE (iw, '(T2,A,ES20.10)') "BADER| Vacuum density [e/bohr^3] ", cutoff
      WRITE (iw, '(T2,A)') "BADER| Atom                 Electrons        Net charge [e]       Volume [bohr^3]"
      a = 0
      DO i = 1, SIZE(particles)
         IF (.NOT. nuclear_center(i)) THEN
            WRITE (iw, '(T2,A,I7,2X,A)') "BADER|", i, "Basis-only centre (excluded)"
            CYCLE
         END IF
         a = a + 1
         IF (assigned(a)) THEN
            WRITE (iw, '(T2,A,I7,3F23.10)') "BADER|", i, basins(1, a), ionic_charge(i) - basins(1, a), basins(2, a)
         ELSE
            WRITE (iw, '(T2,A,I7,2X,A)') "BADER|", i, "No unique basin"
         END IF
      END DO
      WRITE (iw, '(T2,A,I10)') "BADER| Unassigned maxima       ", n_other
      WRITE (iw, '(T2,A,F20.10)') "BADER| Unassigned electrons    ", other(1)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Vacuum electrons        ", vacuum(1)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Unassigned volume [bohr^3] ", other(2)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Vacuum volume [bohr^3]      ", vacuum(2)
      integrated = SUM(basins(1, :)) + other(1) + vacuum(1)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Integrated electrons    ", integrated
      WRITE (iw, '(T2,A,F20.10)') "BADER| Nominal electrons        ", nominal_electrons
      WRITE (iw, '(T2,A,ES20.10)') "BADER| Density minus nominal electrons ", integrated - nominal_electrons
      WRITE (iw, '(T2,A,ES20.10)') "BADER| Electron balance error   ", integrated - SUM(values)*PRODUCT(spacing)
   END SUBROUTINE print_bader_analysis

! **************************************************************************************************
!> \brief Gather distributed grid pencils without allocating the full grid on each rank.
!> \param rho electronic number density for each spin
!> \param iw print unit, positive only on the output rank
!> \param density spin-summed electron density on the print rank
! **************************************************************************************************
   SUBROUTINE collect_density(rho, iw, density)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho(:)
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: density(:, :, :)

      INTEGER                                            :: dest(2), ext(3), hi(3), i, idx, lo(3), &
                                                            local_box(6), local_size, n(3), np, pe
      INTEGER, ALLOCATABLE                               :: boxes(:, :), counts(:), offsets(:)
      REAL(KIND=dp), ALLOCATABLE                         :: collected(:), local_values(:)
      TYPE(mp_comm_type)                                 :: comm

      n = rho(1)%pw_grid%npts
      IF (rho(1)%pw_grid%para%mode == PW_MODE_LOCAL) THEN
         IF (iw <= 0) RETURN
         ALLOCATE (density(n(1), n(2), n(3)))
         density(:, :, :) = rho(1)%array
         DO i = 2, SIZE(rho)
            density(:, :, :) = density + rho(i)%array
         END DO
         RETURN
      END IF

      comm = rho(1)%pw_grid%para%group
      np = comm%num_pe
      pe = comm%mepos
      dest = [iw, pe]
      CALL comm%maxloc(dest)
      CPASSERT(dest(1) > 0)

      local_box(1:3) = rho(1)%pw_grid%bounds_local(1, :) - rho(1)%pw_grid%bounds(1, :) + 1
      local_box(4:6) = local_box(1:3) + rho(1)%pw_grid%npts_local - 1
      ALLOCATE (boxes(6, np), counts(np), offsets(np))
      CALL comm%allgather(local_box, boxes)
      DO idx = 1, np
         counts(idx) = PRODUCT(MAX(0, boxes(4:6, idx) - boxes(1:3, idx) + 1))
         offsets(idx) = SUM(counts(:idx - 1))
      END DO
      CPASSERT(SUM(counts) == PRODUCT(n))
      local_size = counts(pe + 1)
      ALLOCATE (local_values(local_size))
      IF (local_size > 0) THEN
         local_values(:) = RESHAPE(rho(1)%array, [local_size])
         DO i = 2, SIZE(rho)
            local_values(:) = local_values + RESHAPE(rho(i)%array, [local_size])
         END DO
      END IF
      ALLOCATE (collected(MERGE(PRODUCT(n), 0, iw > 0)))
      CALL comm%gatherv(local_values, collected, counts, offsets, dest(2))
      IF (iw <= 0) RETURN
      ALLOCATE (density(n(1), n(2), n(3)))
      DO idx = 1, np
         IF (counts(idx) == 0) CYCLE
         lo = boxes(1:3, idx)
         hi = boxes(4:6, idx)
         ext = hi - lo + 1
         density(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = &
            RESHAPE(collected(offsets(idx) + 1:offsets(idx) + counts(idx)), ext)
      END DO
   END SUBROUTINE collect_density
END MODULE qs_bader_analysis
