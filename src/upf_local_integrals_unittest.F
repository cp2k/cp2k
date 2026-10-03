!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM upf_local_integrals_unittest
   USE ai_overlap_ppl,                  ONLY: ppl_integral,&
                                              ppl_integral_ri
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              ncoset
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init,&
                                              upf_local_matrix,&
                                              upf_local_ri,&
                                              upf_local_type

   IMPLICIT NONE
   LOGICAL :: source_root_argument

   INTEGER, PARAMETER :: la = 4, lb = 4, na = 2, nb = 2
   REAL(KIND=dp), PARAMETER :: alpha(na) = [.1_dp, 3._dp], beta(nb) = [.3_dp, 2.1_dp], gamma(2) = [.2_dp, 4._dp]
   REAL(KIND=dp), PARAMETER :: coeff(1, 2) = RESHAPE([-.7_dp, 2.3_dp], [1, 2])
   INTEGER :: n, grid, geometry, i, j, k, ld, ma, mb
   CHARACTER(len=default_path_length) :: filename
   REAL(KIND=dp) :: rabv(3), rac(3), rbc(3), dx, error, force_a(3), force_b(3), ref_a(3), ref_b(3), scale, threshold
   REAL(KIND=dp), ALLOCATABLE :: vri(:), dri(:, :), rri(:), rdri(:, :), energy_only(:, :)
   REAL(KIND=dp), ALLOCATABLE :: r(:), rab(:), potential(:), v(:, :), da(:, :, :), db(:, :, :), &
                                 reference(:, :), density(:, :), s(:, :, :), fs(:, :, :)
   TYPE(upf_local_type) :: cache
   TYPE(atom_upfpot_type) :: pot
   CALL cp2k_init()
   CALL init_orbital_pointers(2*MAX(la, lb) + 3)
   DO grid = 1, 4
      n = 2501
      IF (grid == 2) n = 2500
      ALLOCATE (r(n), rab(n), potential(n))
      DO i = 1, n
         IF (grid == 3) THEN
            dx = LOG(1 + 12.5_dp/.001_dp)/(n - 1)
            r(i) = .001_dp*(EXP((i - 1)*dx) - 1)
            rab(i) = dx*(r(i) + .001_dp)
         ELSE
            r(i) = (i - 1)*.005_dp
            rab(i) = .005_dp
         END IF
      END DO
      potential(:) = coeff(1, 1)*EXP(-gamma(1)*r*r) + coeff(1, 2)*EXP(-gamma(2)*r*r)
      ALLOCATE (pot%r(n), pot%rab(n), pot%vlocal(n))
      pot%r(:) = r
      pot%rab(:) = rab
      pot%zion = 6._dp
      pot%vlocal(:) = potential
      DO i = 1, n
         IF (r(i) > 0._dp) THEN
            pot%vlocal(i) = potential(i) - pot%zion*erf(SQRT(.8_dp)*r(i))/r(i)
         ELSE
            pot%vlocal(i) = potential(i) - 2*pot%zion*SQRT(.8_dp/pi)
         END IF
      END DO
      threshold = 0._dp
      ! The coarse default EPS_PPL may screen neighbors, never radial samples.
      IF (grid == 4) threshold = 1.e-2_dp
      CALL upf_local_init(pot, .8_dp, la + lb + 2, cache, threshold)
      CALL check_fourier(cache)
      ma = na*ncoset(la)
      mb = nb*ncoset(lb)
      ld = ncoset(2*MAX(la, lb) + 2)
      ALLOCATE (v(ma, mb), da(ma, mb, 3), db(ma, mb, 3), reference(ma, mb), density(ma, mb), s(ld, ld, la + lb + 3), fs(ld, ld, 4))
      ALLOCATE (vri(ma), dri(ma, 3), rri(ma), rdri(ma, 3), energy_only(ma, mb))
      DO j = 1, mb
         DO i = 1, ma
            density(i, j) = SIN(REAL(17*i + 3*j, dp))/REAL(ma*mb, dp)
         END DO
      END DO
      DO geometry = 1, 4
         SELECT CASE (geometry)
         CASE (1)
            rabv = 0; rac = 0
         CASE (2)
            rabv = [.3_dp, -.5_dp, .8_dp]; rac = [.2_dp, .4_dp, -.7_dp]
         CASE (3)
            rabv = [3._dp, 2._dp, -1._dp]; rac = [-1._dp, 2._dp, .3_dp]
         CASE (4)
            rabv = [-7._dp, 1._dp, 2._dp]; rac = [.5_dp, -.2_dp, 1.7_dp]
         END SELECT
         rbc = rac - rabv
         CALL upf_local_matrix(cache, la, alpha, lb, beta, rabv, rac, v, da, db)
         CALL upf_local_matrix(cache, la, alpha, lb, beta, rabv, rac, energy_only)
         IF (MAXVAL(ABS(v - energy_only))/MAX(1._dp, MAXVAL(ABS(v))) > 1.e-12_dp) THEN
            ERROR STOP 'Energy-only and derivative local operators differ'
         END IF
         reference = 0; s = 0; fs = 0
         CALL ppl_integral(la, 0, na, SPREAD(100._dp, 1, na), alpha, lb, 0, nb, SPREAD(100._dp, 1, nb), beta, &
                           2, gamma, [1, 1], coeff, 100._dp, rabv, norm2(rabv), rac, norm2(rac), rbc, norm2(rbc), reference, s, &
                           density, force_a, force_b, fs)
         scale = MAX(1._dp, MAXVAL(ABS(reference)))
         error = MAXVAL(ABS(v - reference))/scale
         DO k = 1, 3
            ref_a(k) = SUM(density*da(:, :, k))
            ref_b(k) = SUM(density*db(:, :, k))
         END DO
         WRITE (*, '(A,2I3,3ES18.9)') 'Native local matrix/derivative errors: ', grid, geometry, &
            error, MAXVAL(ABS(force_a - ref_a)), MAXVAL(ABS(force_b - ref_b))
         IF (error > 1.e-10_dp) ERROR STOP 'Local radial integral differs from analytic three-center reference'
         IF (MAXVAL(ABS(force_a - ref_a)) > 1.e-10_dp .OR. MAXVAL(ABS(force_b - ref_b)) > 1.e-10_dp) THEN
            ERROR STOP 'Local derivative mismatch'
         END IF
         CALL upf_local_ri(cache, la, alpha, rac, vri, dri)
         rri = 0; rdri = 0
         CALL ppl_integral_ri(la, 0, na, SPREAD(100._dp, 1, na), alpha, 2, gamma, [1, 1], coeff, 100._dp, &
                              -rac, norm2(rac), rri, rdri)
         scale = MAX(1._dp, MAXVAL(ABS(rri)), MAXVAL(ABS(rdri)))
         error = MAX(MAXVAL(ABS(vri - rri)), MAXVAL(ABS(dri - rdri)))/scale
         PRINT '(A,2I3,ES18.9)', 'Native local RI errors: ', grid, geometry, error
         IF (error > 1.e-10_dp) ERROR STOP 'Local RI operator differs from analytic Gaussian integrals'
         CALL upf_local_ri(cache, la, alpha, rac, vri)
         IF (MAXVAL(ABS(vri - rri))/scale > 1.e-10_dp) ERROR STOP 'Energy-only RI local operator mismatch'
      END DO
      DEALLOCATE (r, rab, potential, v, da, db, reference, density, s, fs, vri, dri, rri, rdri, energy_only)
      CALL atom_release_upf(pot)
   END DO
   ! Optional real-library audit checks finiteness, not physical accuracy.
   ma = na*ncoset(la)
   mb = nb*ncoset(lb)
   ALLOCATE (v(ma, mb), da(ma, mb, 3), db(ma, mb, 3))
   DO i = 1, command_argument_count()
      CALL get_command_argument(i, filename)
      ! do_regtest.py passes the CP2K source root before any optional UPF audit files.
      IF (i == 1) THEN
         INQUIRE (FILE=TRIM(filename)//"/tests/UNIT_TESTS", EXIST=source_root_argument)
         IF (source_root_argument) CYCLE
      END IF
      CALL atom_read_upf(pot, TRIM(filename))
      CALL upf_local_init(pot, .8_dp, la + lb + 1, cache, 1.e-12_dp)
      rabv = [.3_dp, -.5_dp, .8_dp]; rac = [.2_dp, .4_dp, -.7_dp]
      CALL upf_local_matrix(cache, la, alpha, lb, beta, rabv, rac, v, da, db)
      IF (.NOT. ALL(ieee_is_finite(v)) .OR. .NOT. ALL(ieee_is_finite(da)) .OR. .NOT. ALL(ieee_is_finite(db))) THEN
         ERROR STOP 'Non-finite native local UPF integrals'
      END IF
      PRINT '(A,1X,A,2ES18.9)', 'Native local UPF:', TRIM(filename), cache%radial%radius, MAXVAL(ABS(v))
      CALL atom_release_upf(pot)
   END DO
   DEALLOCATE (v, da, db)
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Independent signed-Gaussian transforms and their center/G derivatives.
!> \param local ...
! **************************************************************************************************
   SUBROUTINE check_fourier(local)
      TYPE(upf_local_type), INTENT(IN)                   :: local

      COMPLEX(KIND=dp)                                   :: analytic(8), dg(3, 8), dr(3, 8), &
                                                            fourier(8), only_value(8), phase, &
                                                            ref_dg(3, 8), ref_dr(3, 8)
      INTEGER                                            :: ig, term
      REAL(KIND=dp)                                      :: g(3, 8), position(3), radial, &
                                                            transform_error, volume

      g(:, 1) = 0.0_dp
      g(:, 2) = [1.e-9_dp, -2.e-9_dp, 3.e-9_dp]
      g(:, 3) = [0.3_dp, -0.4_dp, 0.7_dp]
      g(:, 4) = -g(:, 3)
      g(:, 5) = [0.0_dp, 0.0_dp, 2.0_dp]
      g(:, 6) = [4.0_dp, -3.0_dp, 1.0_dp]
      g(:, 7) = [10.0_dp, 0.0_dp, 0.0_dp]
      g(:, 8) = [20.0_dp, 2.0_dp, -1.0_dp]
      position = [0.4_dp, -0.7_dp, 1.2_dp]
      volume = 125.0_dp
      analytic = 0.0_dp
      ref_dg = 0.0_dp
      DO ig = 1, SIZE(analytic)
         phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g(:, ig), position), dp))/volume
         DO term = 1, SIZE(gamma)
            radial = coeff(1, term)*(pi/gamma(term))**1.5_dp*EXP(-SUM(g(:, ig)**2)/(4.0_dp*gamma(term)))
            analytic(ig) = analytic(ig) + phase*radial
            ref_dg(:, ig) = ref_dg(:, ig) + phase*radial*(-g(:, ig)/(2.0_dp*gamma(term)) + &
                                                          CMPLX(0.0_dp, -position, dp))
         END DO
         ref_dr(:, ig) = CMPLX(0.0_dp, -g(:, ig), dp)*analytic(ig)
      END DO
      CALL upf_local_fourier(local, g, position, volume, fourier, dr, dg)
      CALL upf_local_fourier(local, g, position, volume, only_value)
      transform_error = MAX(MAXVAL(ABS(fourier - analytic)), MAXVAL(ABS(dr - ref_dr)), &
                            MAXVAL(ABS(dg - ref_dg)), MAXVAL(ABS(only_value - fourier)))
      WRITE (*, '(A,ES23.15)') 'Native local Fourier/derivative maximum error: ', transform_error
      IF (transform_error > 1.e-10_dp) ERROR STOP 'Local Fourier potential differs from analytic Gaussian transform'
      IF (ABS(fourier(4) - CONJG(fourier(3))) > 1.e-14_dp) ERROR STOP 'Local Fourier potential lacks conjugation symmetry'
   END SUBROUTINE check_fourier
END PROGRAM upf_local_integrals_unittest
