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
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_bader_basins,                 ONLY: bader_integrate_grid
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
!> \param density_label description of the density used to define and integrate basins
!> \param cutoff density cutoff for vacuum integration
!> \param assign_nearest group maxima by proximity rather than unique nuclear seed
!> \param iw print unit, positive only on the output rank
! **************************************************************************************************
   SUBROUTINE print_bader_analysis(rho, cell, particles, ionic_charge, density_label, cutoff, assign_nearest, iw)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho(:)
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      REAL(KIND=dp), INTENT(IN)                          :: ionic_charge(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: density_label
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, INTENT(IN)                                :: assign_nearest
      INTEGER, INTENT(IN)                                :: iw

      INTEGER, PARAMETER                                 :: max_grid_points = 8000000

      INTEGER                                            :: d, i, n(3), n_extra, n_other, ng, q(3)
      INTEGER, ALLOCATABLE                               :: seeds(:)
      LOGICAL                                            :: periodic(3)
      LOGICAL, ALLOCATABLE                               :: assigned(:)
      REAL(KIND=dp)                                      :: other(2), SPACING(3), vacuum(2)
      REAL(KIND=dp), ALLOCATABLE                         :: atom_grid(:, :), basins(:, :)
      REAL(KIND=dp), ALLOCATABLE, TARGET                 :: density(:, :, :)
      REAL(KIND=dp), POINTER                             :: values(:)

      IF (.NOT. cell%orthorhombic) CPABORT("BADER requires a Cartesian orthorhombic cell.")
      periodic = cell%perd /= 0
      IF (.NOT. (cutoff >= 0.0_dp)) CPABORT("BADER requires a nonnegative VACUUM_DENSITY.")
      IF (rho(1)%pw_grid%ngpts > max_grid_points) THEN
         CPABORT("BADER supports at most 8000000 grid points in its current implementation.")
      END IF
      n = rho(1)%pw_grid%npts
      ng = PRODUCT(n)
      CALL collect_density(rho, iw, density)
      IF (iw <= 0) RETURN

      ALLOCATE (seeds(SIZE(particles)), assigned(SIZE(particles)), basins(2, SIZE(particles)))
      IF (assign_nearest) ALLOCATE (atom_grid(3, SIZE(particles)))
      spacing = [rho(1)%pw_grid%dh(1, 1), rho(1)%pw_grid%dh(2, 2), rho(1)%pw_grid%dh(3, 3)]
      seeds = 0
      DO i = 1, SIZE(particles)
         IF (assign_nearest) atom_grid(:, i) = particles(i)%r/spacing
         q = NINT(particles(i)%r/spacing)
         DO d = 1, 3
            IF (periodic(d)) q(d) = MODULO(q(d), n(d))
         END DO
         IF (ANY(q < 0) .OR. ANY(q >= n)) CYCLE
         seeds(i) = 1 + q(1) + n(1)*(q(2) + n(2)*q(3))
      END DO
      values(1:ng) => density
      CALL bader_integrate_grid(values, n, spacing, periodic, cutoff, seeds, basins, assigned, other, vacuum, &
                                n_other, atom_grid=atom_grid, n_extra=n_extra)

      WRITE (iw, '(/,T2,A)') "BADER| Weighted basins of "//density_label
      IF (assign_nearest) THEN
         WRITE (iw, '(T2,A)') "BADER| Atom assignment: nearest atom to each maximum"
      ELSE
         WRITE (iw, '(T2,A)') "BADER| Atom assignment: unique nuclear seed only"
      END IF
      WRITE (iw, '(T2,A,I10)') "BADER| Maxima without a unique nuclear seed ", n_extra
      WRITE (iw, '(T2,A,ES20.10)') "BADER| Vacuum density [e/bohr^3] ", cutoff
      WRITE (iw, '(T2,A)') "BADER| Atom                 Electrons        Net charge [e]       Volume [bohr^3]"
      DO i = 1, SIZE(particles)
         IF (assigned(i)) THEN
            WRITE (iw, '(T2,A,I7,3F23.10)') "BADER|", i, basins(1, i), ionic_charge(i) - basins(1, i), basins(2, i)
         ELSE
            WRITE (iw, '(T2,A,I7,2X,A)') "BADER|", i, "No unique basin"
         END IF
      END DO
      WRITE (iw, '(T2,A,I10)') "BADER| Unassigned maxima       ", n_other
      WRITE (iw, '(T2,A,F20.10)') "BADER| Unassigned electrons    ", other(1)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Vacuum electrons        ", vacuum(1)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Unassigned volume [bohr^3] ", other(2)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Vacuum volume [bohr^3]      ", vacuum(2)
      WRITE (iw, '(T2,A,F20.10)') "BADER| Integrated electrons    ", SUM(basins(1, :)) + other(1) + vacuum(1)
      WRITE (iw, '(T2,A,ES20.10)') "BADER| Electron balance error   ", &
         SUM(basins(1, :)) + other(1) + vacuum(1) - SUM(values)*PRODUCT(spacing)
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
