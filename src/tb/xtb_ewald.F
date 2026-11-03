!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Full Klopman-Ohno Ewald summation for native xTB.
!>        For harmonic shell hardness h, gamma = 1/r + R - (a_i+a_j)^2 S3/8,
!>        a_i=1/h_i. R decays as r^-5; S3 is the generalized Ewald sum of r^-3.
!>        Uses the finite-part convention and 40 bohr residual translation box of tblite.
!>        See J. Chem. Theory Comput. 2025, 21, 1615-1625, doi:10.1021/acs.jctc.4c01234.
!>        The ordinary 1/r contribution is supplied by the existing SPME path.
! **************************************************************************************************
MODULE xtb_ewald
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: euler,&
                                              pi
   USE mathlib,                         ONLY: expint
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_ewald'
   REAL(dp), PARAMETER :: ko_cutoff = 40.0_dp
   REAL(dp), PARAMETER :: s3_constant = 0.5_dp*(LOG(pi) - 2.0_dp + euler + 2.0_dp*LOG(2.0_dp))

!> Geometry cache owned by a QS environment, never shared between force evaluations.
!> Real-space blocks are distributed by atom rows; reciprocal vectors by vector index.
!> Derivative tensors and the N_atom x N_G phase matrix are deliberately not cached.
   TYPE xtb_ewald_type
      INTEGER :: rank = -1, nproc = 0
      REAL(dp) :: alpha = -1.0_dp, volume = 0.0_dp, k0 = 0.0_dp
      REAL(dp) :: h(3, 3) = 0.0_dp
      INTEGER, ALLOCATABLE :: nshell(:), rows(:)
      REAL(dp), ALLOCATABLE :: xyz(:, :), hardness(:, :)
      REAL(dp), ALLOCATABLE :: strans(:, :), dtrans(:, :)
      REAL(dp), ALLOCATABLE :: vec(:, :), weight(:), strain(:, :, :)
      REAL(dp), ALLOCATABLE :: blocks(:, :, :, :)
   END TYPE xtb_ewald_type

   PUBLIC :: xtb_ewald_type, update_xtb_ewald, apply_xtb_ewald, release_xtb_ewald

CONTAINS

! **************************************************************************************************
!> \brief Release all geometry-dependent data; also used on QS environment reuse.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE release_xtb_ewald(cache)
      TYPE(xtb_ewald_type), INTENT(INOUT)                :: cache

      IF (ALLOCATED(cache%xyz)) THEN
         DEALLOCATE (cache%xyz, cache%hardness, cache%nshell, cache%rows, cache%blocks)
         DEALLOCATE (cache%strans, cache%dtrans, cache%vec, cache%weight, cache%strain)
      END IF
      cache%rank = -1
      cache%nproc = 0
      cache%alpha = -1.0_dp
   END SUBROUTINE release_xtb_ewald

! **************************************************************************************************
!> \brief Rebuild on changes in geometry, hardness, splitting parameter or MPI partition.
!> \param cache ...
!> \param h Cell vectors in columns
!> \param h_inv Inverse cell matrix
!> \param xyz Cartesian coordinates in the input images, shape (3,natom)
!> \param hardness Shell Hubbard parameters, shape (natom,maxshell)
!> \param nshell Number of active shells per atom
!> \param alpha Ewald splitting parameter
!> \param rank Zero-based rank within the QS communicator
!> \param nproc Number of ranks within the QS communicator
! **************************************************************************************************
   SUBROUTINE update_xtb_ewald(cache, h, h_inv, xyz, hardness, nshell, alpha, rank, nproc)
      TYPE(xtb_ewald_type), INTENT(INOUT)                :: cache
      REAL(dp), INTENT(IN)                               :: h(3, 3), h_inv(3, 3), xyz(:, :), &
                                                            hardness(:, :)
      INTEGER, INTENT(IN)                                :: nshell(:)
      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: rank, nproc

      INTEGER                                            :: b, c, handle, iat, ig, iloc, jat, jg, &
                                                            ms, nat, ng, nr, padding(3)
      REAL(dp)                                           :: cutoff, fmax(3), fmin(3), frac(3), g2, &
                                                            pref, rec(3, 3), rec_inv(3, 3), &
                                                            volume, w, x
      REAL(dp), ALLOCATABLE                              :: trans(:, :)

      CALL timeset('update_xtb_ewald', handle)
      CPASSERT(alpha > 0.0_dp .AND. nproc > 0 .AND. rank >= 0 .AND. rank < nproc)
      nat = SIZE(nshell)
      ms = SIZE(hardness, 2)
      CPASSERT(SIZE(xyz, 1) == 3 .AND. SIZE(xyz, 2) == nat .AND. SIZE(hardness, 1) == nat)
      CPASSERT(ALL(nshell >= 0) .AND. ALL(nshell <= ms))
      IF (ALLOCATED(cache%xyz)) THEN
         IF (SIZE(cache%xyz, 2) == nat .AND. SIZE(cache%hardness, 2) == ms) THEN
            IF (ALL(cache%h == h) .AND. cache%alpha == alpha .AND. &
                cache%rank == rank .AND. cache%nproc == nproc) THEN
               IF (ALL(cache%xyz == xyz) .AND. ALL(cache%hardness == hardness) .AND. &
                   ALL(cache%nshell == nshell)) THEN
                  CALL timestop(handle)
                  RETURN
               END IF
            END IF
         END IF
         CALL release_xtb_ewald(cache)
      END IF
      volume = ABS(h(1, 1)*(h(2, 2)*h(3, 3) - h(2, 3)*h(3, 2)) &
                   - h(1, 2)*(h(2, 1)*h(3, 3) - h(2, 3)*h(3, 1)) &
                   + h(1, 3)*(h(2, 1)*h(3, 2) - h(2, 2)*h(3, 1)))
      CPASSERT(volume > 0.0_dp)
      nr = MAX(0, (nat + nproc - 1 - rank)/nproc)
      ALLOCATE (cache%xyz(3, nat), cache%hardness(nat, ms), cache%nshell(nat), cache%rows(nr))
      ALLOCATE (cache%blocks(MAXVAL(nshell), MAXVAL(nshell), nat, nr))
      cache%xyz(:, :) = xyz
      cache%hardness(:, :) = hardness
      cache%nshell(:) = nshell
      cache%h(:, :) = h
      cache%alpha = alpha
      cache%volume = volume
      cache%rank = rank
      cache%nproc = nproc
      cache%k0 = 4.0_dp*pi/volume*(LOG(alpha/SQRT(pi)) + s3_constant)
      DO iloc = 1, nr
         cache%rows(iloc) = rank + 1 + (iloc - 1)*nproc
      END DO
      CALL lattice_box(h, h_inv, ko_cutoff, cache%strans)
      ! Extend the screened real-space box when input images span more than one cell.
      fmin(:) = HUGE(1.0_dp)
      fmax(:) = -HUGE(1.0_dp)
      DO iat = 1, nat
         frac(:) = MATMUL(h_inv, xyz(:, iat))
         fmin(:) = MIN(fmin, frac)
         fmax(:) = MAX(fmax, frac)
      END DO
      padding(:) = MAX(0, CEILING(fmax - fmin) - 1)
      cutoff = SQRT(-LOG(EPSILON(1.0_dp)))/alpha
      CALL lattice_box(h, h_inv, cutoff, cache%dtrans, padding)
      rec(:, :) = 2.0_dp*pi*TRANSPOSE(h_inv)
      rec_inv(:, :) = TRANSPOSE(h)/(2.0_dp*pi)
      cutoff = 2.0_dp*alpha*SQRT(-LOG(EPSILON(1.0_dp)) + MAX(0.0_dp, LOG(2.0_dp*pi/volume)))
      CALL lattice_box(rec, rec_inv, cutoff, trans)
      ! Origin is the first vector. Both signs are retained, as in tblite.
      ng = MAX(0, (SIZE(trans, 2) - 1 + nproc - 1 - rank)/nproc)
      ALLOCATE (cache%vec(3, ng), cache%weight(ng), cache%strain(3, 3, ng))
      pref = 4.0_dp*pi/volume
      DO ig = 1, ng
         jg = 2 + rank + (ig - 1)*nproc
         cache%vec(:, ig) = trans(:, jg)
         g2 = SUM(cache%vec(:, ig)**2)
         x = g2/(4.0_dp*alpha**2)
         cache%weight(ig) = 0.5_dp*pref*expint(1, x)
         w = pref*EXP(-x)/g2
         DO c = 1, 3
            DO b = 1, 3
               cache%strain(b, c, ig) = w*cache%vec(b, ig)*cache%vec(c, ig)
            END DO
            cache%strain(c, c, ig) = cache%strain(c, c, ig) - cache%weight(ig)
         END DO
      END DO
      DEALLOCATE (trans)
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(cache, nr, nat) PRIVATE(iloc, iat, jat)
      DO iloc = 1, nr
         iat = cache%rows(iloc)
         DO jat = 1, nat
            CALL ko_real_pair(cache, iat, jat, cache%blocks(:, :, jat, iloc))
         END DO
      END DO
      !$OMP END PARALLEL DO
      CALL timestop(handle)
   END SUBROUTINE update_xtb_ewald

! **************************************************************************************************
!> \brief Add the local contribution to potential, negative Cartesian gradient and strain.
!>        Caller reduces potentials over MPI; forces and virial use CP2K's later reduction.
!>        The three charge moments give the exact harmonic-hardness contraction of S3.
!> \param cache ...
!> \param charges Shell charges, shape (natom,maxshell)
!> \param gchrg Potential and, if extent 3 is 4, minus its Cartesian derivatives
!> \param sigma Optional local energy derivative with respect to strain (positive sign)
! **************************************************************************************************
   SUBROUTINE apply_xtb_ewald(cache, charges, gchrg, sigma)
      TYPE(xtb_ewald_type), INTENT(IN)                   :: cache
      REAL(dp), INTENT(IN)                               :: charges(:, :)
      REAL(dp), INTENT(INOUT)                            :: gchrg(:, :, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: sigma(3, 3)

      INTEGER                                            :: b, c, handle, iat, ig, iloc, ish, jat, &
                                                            jsh, nat, ng, ni, nj
      LOGICAL                                            :: deriv, stress
      REAL(dp)                                           :: a, coeff, ct, dv(3), field(3, 3), &
                                                            sig(3, 3), st, theta, total(3), v(3)
      REAL(dp), ALLOCATABLE                              :: moment(:, :), sc(:, :), ss(:, :)
      REAL(dp) :: ds(3, 3, SIZE(charges, 2), SIZE(charges, 2)), &
         dg(3, SIZE(charges, 2), SIZE(charges, 2)), block(SIZE(charges, 2), SIZE(charges, 2))

      CALL timeset('apply_xtb_ewald', handle)
      nat = SIZE(charges, 1)
      ng = SIZE(cache%weight)
      deriv = SIZE(gchrg, 3) == 4
      stress = PRESENT(sigma)
      CPASSERT(SIZE(gchrg, 3) == 1 .OR. deriv)
      CPASSERT(.NOT. stress .OR. deriv)
      ALLOCATE (moment(3, nat), sc(3, ng), ss(3, ng))
      moment(:, :) = 0.0_dp
      DO iat = 1, nat
         DO ish = 1, cache%nshell(iat)
            a = 1.0_dp/cache%hardness(iat, ish)
            moment(1, iat) = moment(1, iat) + charges(iat, ish)
            moment(2, iat) = moment(2, iat) + charges(iat, ish)*a
            moment(3, iat) = moment(3, iat) + charges(iat, ish)*a*a
         END DO
      END DO
      ! CP2K stores the negative potential gradient; rho_elec later stores the energy gradient.
      sig(:, :) = 0.0_dp
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(cache, charges, gchrg, nat, deriv, stress) &
      !$OMP PRIVATE(iloc, iat, jat, ni, nj, ish, jsh, b, block, dg, ds) REDUCTION(+:sig)
      DO iloc = 1, SIZE(cache%rows)
         iat = cache%rows(iloc)
         ni = cache%nshell(iat)
         DO jat = 1, nat
            nj = cache%nshell(jat)
            DO jsh = 1, nj
               gchrg(iat, 1:ni, 1) = gchrg(iat, 1:ni, 1) &
                                     + cache%blocks(1:ni, jsh, jat, iloc)*charges(jat, jsh)
            END DO
            IF (.NOT. deriv) CYCLE
            IF (stress) THEN
               CALL ko_real_pair(cache, iat, jat, block, dg, ds)
            ELSE
               CALL ko_real_pair(cache, iat, jat, block, dg)
            END IF
            DO b = 1, 3
               ! A self-image displacement does not move when the atom moves.
               IF (iat /= jat) THEN
                  DO jsh = 1, nj
                     gchrg(iat, 1:ni, b + 1) = gchrg(iat, 1:ni, b + 1) &
                                               - dg(b, 1:ni, jsh)*charges(jat, jsh)
                  END DO
               END IF
            END DO
            IF (stress) THEN
               DO jsh = 1, nj
                  DO ish = 1, ni
                     sig(:, :) = sig(:, :) + 0.5_dp*charges(iat, ish)*charges(jat, jsh)*ds(:, :, ish, jsh)
                  END DO
               END DO
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO
      sc(:, :) = 0.0_dp
      ss(:, :) = 0.0_dp
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(cache, moment, sc, ss, nat, ng) PRIVATE(ig, iat, theta)
      DO ig = 1, ng
         DO iat = 1, nat
            theta = DOT_PRODUCT(cache%xyz(:, iat), cache%vec(:, ig))
            sc(:, ig) = sc(:, ig) + COS(theta)*moment(:, iat)
            ss(:, ig) = ss(:, ig) + SIN(theta)*moment(:, iat)
         END DO
      END DO
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(cache, gchrg, sc, ss, nat, ng, deriv) &
      !$OMP PRIVATE(iat, ig, ish, b, theta, ct, st, v, dv, field, A)
      DO iat = 1, nat
         v(:) = 0.0_dp
         field(:, :) = 0.0_dp
         DO ig = 1, ng
            theta = DOT_PRODUCT(cache%xyz(:, iat), cache%vec(:, ig))
            ct = COS(theta)
            st = SIN(theta)
            v(:) = v(:) + cache%weight(ig)*(ct*sc(:, ig) + st*ss(:, ig))
            IF (deriv) THEN
               dv(:) = cache%weight(ig)*(-st*sc(:, ig) + ct*ss(:, ig))
               DO b = 1, 3
                  field(:, b) = field(:, b) + cache%vec(b, ig)*dv(:)
               END DO
            END IF
         END DO
         DO ish = 1, cache%nshell(iat)
            a = 1.0_dp/cache%hardness(iat, ish)
            gchrg(iat, ish, 1) = gchrg(iat, ish, 1) - (a*a*v(1) + 2.0_dp*a*v(2) + v(3))/8.0_dp
            IF (deriv) THEN
               DO b = 1, 3
                  gchrg(iat, ish, b + 1) = gchrg(iat, ish, b + 1) &
                                           + (a*a*field(1, b) + 2.0_dp*a*field(2, b) + field(3, b))/8.0_dp
               END DO
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO
      IF (stress) THEN
         DO ig = 1, ng
            coeff = -(sc(1, ig)*sc(3, ig) + ss(1, ig)*ss(3, ig) + sc(2, ig)**2 + ss(2, ig)**2)/8.0_dp
            sig(:, :) = sig(:, :) + coeff*cache%strain(:, :, ig)
         END DO
      END IF
      ! G=0 contributes even for a neutral cell: sum(q/h) need not vanish.
      IF (cache%rank == 0) THEN
         total(:) = SUM(moment, DIM=2)
         DO iat = 1, nat
            DO ish = 1, cache%nshell(iat)
               a = 1.0_dp/cache%hardness(iat, ish)
               gchrg(iat, ish, 1) = gchrg(iat, ish, 1) &
                                    - cache%k0*(a*a*total(1) + 2.0_dp*a*total(2) + total(3))/8.0_dp
            END DO
         END DO
         IF (stress) THEN
            coeff = cache%k0*(total(1)*total(3) + total(2)**2)/8.0_dp
            DO c = 1, 3
               sig(c, c) = sig(c, c) + coeff
            END DO
         END IF
      END IF
      IF (stress) sigma(:, :) = sig
      DEALLOCATE (moment, sc, ss)
      CALL timestop(handle)
   END SUBROUTINE apply_xtb_ewald

! **************************************************************************************************
!> \brief Real-space residual and screened r^-3 term for a single atom pair.
!> \param cache ...
!> \param iat ...
!> \param jat ...
!> \param block ...
!> \param dg Cartesian derivative with respect to r_i-r_j
!> \param ds Strain derivative at fixed alpha (including self images)
! **************************************************************************************************
   SUBROUTINE ko_real_pair(cache, iat, jat, block, dg, ds)
      TYPE(xtb_ewald_type), INTENT(IN)                   :: cache
      INTEGER, INTENT(IN)                                :: iat, jat
      REAL(dp), INTENT(OUT)                              :: block(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: dg(:, :, :), ds(:, :, :, :)

      INTEGER                                            :: b, c, ish, it, jsh
      LOGICAL                                            :: deriv, stress
      REAL(dp)                                           :: alpha, dg3(3), ds3(3, 3), dval, ef, ex, &
                                                            gam2, invr, invsqrtpi, r(3), r2, &
                                                            rij(3), rr, s, s3, val, x

      block(:, :) = 0.0_dp
      deriv = PRESENT(dg)
      stress = PRESENT(ds)
      IF (deriv) dg(:, :, :) = 0.0_dp
      IF (stress) ds(:, :, :, :) = 0.0_dp
      rij(:) = cache%xyz(:, iat) - cache%xyz(:, jat)
      alpha = cache%alpha
      invsqrtpi = 1.0_dp/SQRT(pi)
      s3 = 0.0_dp
      dg3(:) = 0.0_dp
      ds3(:, :) = 0.0_dp
      DO it = 1, SIZE(cache%dtrans, 2)
         r(:) = rij + cache%dtrans(:, it)
         r2 = SUM(r*r)
         IF (r2 < 1.0e-20_dp) CYCLE
         rr = SQRT(r2)
         ex = EXP(-alpha**2*r2)*invsqrtpi
         ef = ERFC(alpha*rr)
         s3 = s3 + ef/(r2*rr) + 2.0_dp*alpha*ex/r2
         IF (deriv .OR. stress) THEN
            dval = -3.0_dp*ef/(r2*r2*rr) - 6.0_dp*alpha*ex/r2**2 - 4.0_dp*alpha**3*ex/r2
            dg3(:) = dg3(:) + dval*r
            IF (stress) THEN
               DO c = 1, 3
                  ds3(:, c) = ds3(:, c) + dval*r*r(c)
               END DO
            END IF
         END IF
      END DO
      IF (iat == jat) s3 = s3 - 4.0_dp*alpha**3*invsqrtpi/3.0_dp
      DO jsh = 1, cache%nshell(jat)
         DO ish = 1, cache%nshell(iat)
            gam2 = (0.5_dp/cache%hardness(iat, ish) + 0.5_dp/cache%hardness(jat, jsh))**2
            block(ish, jsh) = -0.5_dp*gam2*s3
            IF (iat == jat) block(ish, jsh) = block(ish, jsh) + 1.0_dp/SQRT(gam2)
            IF (deriv) dg(:, ish, jsh) = -0.5_dp*gam2*dg3
            IF (stress) ds(:, :, ish, jsh) = -0.5_dp*gam2*ds3
            DO it = 1, SIZE(cache%strans, 2)
               r(:) = rij + cache%strans(:, it)
               r2 = SUM(r*r)
               IF (r2 < 1.0e-20_dp) CYCLE
               invr = 1.0_dp/SQRT(r2)
               x = gam2/r2
               s = SQRT(1.0_dp + x)
               ! Stable forms of (1+x)^(-1/2)-1+x/2 and its radial derivative.
               val = invr*x*x*(s + 2.0_dp)/(2.0_dp*s*(s + 1.0_dp)**2)
               block(ish, jsh) = block(ish, jsh) + val
               IF (deriv .OR. stress) THEN
                  dval = -invr**3*x*x*(3.0_dp*s**3 + 6.0_dp*s*s + 4.0_dp*s + 2.0_dp) &
                         /(2.0_dp*s**3*(s + 1.0_dp)**2)
                  IF (deriv) dg(:, ish, jsh) = dg(:, ish, jsh) + dval*r
                  IF (stress) THEN
                     DO b = 1, 3
                        ds(:, b, ish, jsh) = ds(:, b, ish, jsh) + dval*r*r(b)
                     END DO
                  END IF
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE ko_real_pair

! **************************************************************************************************
!> \brief Symmetric lattice translation box with reciprocal-plane bounds, origin first.
!> \param h ...
!> \param h_inv ...
!> \param cutoff ...
!> \param trans ...
!> \param padding ...
! **************************************************************************************************
   SUBROUTINE lattice_box(h, h_inv, cutoff, trans, padding)
      REAL(dp), INTENT(IN)                               :: h(3, 3), h_inv(3, 3), cutoff
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: trans(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: padding(3)

      INTEGER                                            :: b, i, j, k, n, rep(3)

      DO b = 1, 3
         rep(b) = CEILING(cutoff*SQRT(SUM(h_inv(b, :)**2)))
      END DO
      IF (PRESENT(padding)) rep(:) = rep + padding
      ALLOCATE (trans(3, PRODUCT(2*rep + 1)))
      trans(:, 1) = 0.0_dp
      n = 1
      DO k = -rep(3), rep(3)
         DO j = -rep(2), rep(2)
            DO i = -rep(1), rep(1)
               IF (i == 0 .AND. j == 0 .AND. k == 0) CYCLE
               n = n + 1
               trans(:, n) = h(:, 1)*i + h(:, 2)*j + h(:, 3)*k
            END DO
         END DO
      END DO
   END SUBROUTINE lattice_box

END MODULE xtb_ewald
