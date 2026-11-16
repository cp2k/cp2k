!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM upf_projector_integrals_unittest
   USE ai_moments,                      ONLY: moment
   USE ai_overlap,                      ONLY: overlap
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_overlap,&
                                              upf_projector_type

   IMPLICIT NONE
   LOGICAL :: source_root_argument

   INTEGER, PARAMETER :: la = 4, lb = 3, nder = 2, na = 5, nb = 5
   REAL(KIND=dp), PARAMETER :: alpha(na) = [0.04_dp, 0.3_dp, 1.2_dp, 5.0_dp, 24.0_dp]
   REAL(KIND=dp), PARAMETER :: eta(nb) = [0.3_dp, 0.7_dp, 1.1_dp, 2.2_dp, 0.4_dp]
   CHARACTER(LEN=default_path_length) :: filename
   INTEGER :: grid, nr, i, j, p, m, l, geometry, id, ld, ldai, first, shell, q, io_unit
   REAL(KIND=dp) :: dx, rac(3), origin(3), error, scale, expected_d
   REAL(KIND=dp), ALLOCATABLE :: sab(:, :), work(:, :, :), moments(:, :, :), integrals(:, :, :), reference(:, :, :)
   TYPE(atom_upfpot_type) :: pot
   TYPE(upf_projector_type) :: projectors

   CALL cp2k_init()
   io_unit = cp_logger_get_default_io_unit()
   pot%number_of_proj = nb
   ALLOCATE (pot%lbeta(nb))
   pot%lbeta(:) = [0, 1, 2, 3, 1]
   ALLOCATE (pot%dion(nb, nb))
   pot%dion = 0.0_dp
   DO i = 1, nb
      pot%dion(i, i) = 0.3_dp*i - 0.8_dp
   END DO
   pot%dion(2, 5) = 0.17_dp
   pot%dion(5, 2) = 0.17_dp

   DO grid = 1, 3
      nr = 2501
      IF (grid == 2) nr = 2500
      ALLOCATE (pot%r(nr), pot%rab(nr), pot%beta(nr, nb))
      DO i = 1, nr
         IF (grid == 3) THEN
            dx = LOG(1.0_dp + 12.5_dp/0.001_dp)/(nr - 1)
            pot%r(i) = 0.001_dp*(EXP((i - 1)*dx) - 1.0_dp)
            pot%rab(i) = dx*(pot%r(i) + 0.001_dp)
         ELSE
            pot%r(i) = (i - 1)*0.005_dp
            pot%rab(i) = 0.005_dp
         END IF
      END DO
      DO i = 1, nb
         pot%beta(:, i) = pot%r**pot%lbeta(i)*EXP(-eta(i)*pot%r**2)
      END DO
      CALL upf_projector_init(pot, la + nder, projectors)
      ! Independent expansion of D, including non-diagonal coupling of two p projectors.
      p = 0
      DO i = 1, nb
         DO m = -pot%lbeta(i), pot%lbeta(i)
            p = p + 1
            q = 0
            DO j = 1, nb
               DO l = -pot%lbeta(j), pot%lbeta(j)
                  q = q + 1
                  expected_d = 0.0_dp
                  IF (pot%lbeta(i) == pot%lbeta(j) .AND. m == l) expected_d = pot%dion(i, j)
                  IF (ABS(projectors%d(p, q) - expected_d) > 1.e-15_dp) ERROR STOP "UPF D expansion failed"
               END DO
            END DO
         END DO
      END DO
      ld = MAX(na*ncoset(la), nb*ncoset(lb))
      ldai = ncoset(la + nder)
      ALLOCATE (sab(ld, ld*ncoset(nder)), work(ldai, ldai, ncoset(nder)), &
                moments(na*ncoset(la), nb*ncoset(lb), ncoset(nder) - 1), &
                integrals(na*ncoset(la), projectors%nproj, ncoset(nder)), &
                reference(na*ncoset(la), projectors%nproj, ncoset(nder)))
      DO geometry = 1, 4
         SELECT CASE (geometry)
         CASE (1)
            rac = 0.0_dp
         CASE (2)
            rac = [0.2_dp, -0.5_dp, 0.7_dp]
         CASE (3)
            rac = [3.0_dp, 2.0_dp, -1.0_dp]
         CASE (4)
            rac = [-7.0_dp, 1.0_dp, 2.0_dp]
         END SELECT
         sab = 0.0_dp
         CALL overlap(la, 0, na, SPREAD(100.0_dp, 1, na), alpha, &
                      lb, 0, nb, SPREAD(100.0_dp, 1, nb), eta, rac, NORM2(rac), &
                      sab, nder, .TRUE., work, ldai)
         DO id = 1, ncoset(nder)
            p = 0
            DO i = 1, nb
               l = pot%lbeta(i)
               first = (id - 1)*ld + (i - 1)*ncoset(lb) + ncoset(l - 1) + 1
               reference(:, p + 1:p + 2*l + 1, id) = &
                  MATMUL(sab(:na*ncoset(la), first:first + nco(l) - 1), TRANSPOSE(orbtramat(l)%slm))
               p = p + 2*l + 1
            END DO
         END DO
         CALL upf_projector_overlap(projectors, la, alpha, rac, nder, integrals)
         DO id = 1, ncoset(nder)
            scale = MAX(1.0_dp, MAXVAL(ABS(reference(:, :, id))))
            error = MAXVAL(ABS(integrals(:, :, id) - reference(:, :, id)))/scale
            IF (io_unit > 0) WRITE (io_unit, '(A,3I4,ES16.7)') &
               "UPF radial Gaussian overlap/derivative:", grid, geometry, id, error
            IF (error > 1.e-10_dp) ERROR STOP "Numerical UPF overlaps disagree with analytic Gaussian integrals"
         END DO
         ! Independent moment implementation uses the OS recurrence with a fixed reference.
         origin = [0.4_dp, 0.6_dp, -0.3_dp]
         CALL moment(la, na, alpha, SPREAD(100.0_dp, 1, na), 0, &
                     lb, nb, eta, SPREAD(100.0_dp, 1, nb), nder, origin, origin + rac, moments)
         DO id = 2, ncoset(nder)
            p = 0
            DO i = 1, nb
               l = pot%lbeta(i)
               shell = (i - 1)*ncoset(lb) + ncoset(l - 1) + 1
               reference(:, p + 1:p + 2*l + 1, id) = &
                  MATMUL(moments(:, shell:shell + nco(l) - 1, id - 1), TRANSPOSE(orbtramat(l)%slm))
               p = p + 2*l + 1
            END DO
         END DO
         CALL upf_projector_overlap(projectors, la, alpha, rac, nder, integrals, moment_origin=origin)
         scale = MAX(1.0_dp, MAXVAL(ABS(reference)))
         error = MAXVAL(ABS(integrals - reference))/scale
         IF (io_unit > 0) WRITE (io_unit, '(A,2I4,ES16.7)') "UPF radial Gaussian moments:", grid, geometry, error
         IF (error > 1.e-10_dp) ERROR STOP "Numerical UPF moments disagree with analytic Gaussian integrals"
      END DO
      DEALLOCATE (pot%r, pot%rab, pot%beta, sab, work, moments, integrals, reference)
   END DO
   CALL atom_release_upf(pot)
   ! Optional library audit checks real radial meshes, angular channels, and
   ! successive cache replacement. Finiteness is not a physical accuracy test.
   DO i = 1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(i, filename)
      ! do_regtest.py passes the CP2K source root before any optional UPF audit files.
      IF (i == 1) THEN
         INQUIRE (FILE=TRIM(filename)//"/tests/UNIT_TESTS", EXIST=source_root_argument)
         IF (source_root_argument) CYCLE
      END IF
      CALL atom_read_upf(pot, TRIM(filename))
      IF (pot%number_of_proj > 0) THEN
         CALL upf_projector_init(pot, la + nder, projectors)
         IF (.NOT. ALL(ieee_is_finite(projectors%d))) ERROR STOP "Non-finite original UPF D matrix"
         ALLOCATE (integrals(na*ncoset(la), projectors%nproj, ncoset(nder)))
         DO geometry = 1, 2
            rac = 0.0_dp
            IF (geometry == 2) rac = [0.3_dp, -0.5_dp, 1.0_dp]
            CALL upf_projector_overlap(projectors, la, alpha, rac, nder, integrals)
            IF (.NOT. ALL(ieee_is_finite(integrals))) ERROR STOP "Non-finite native UPF projector integrals"
         END DO
         IF (io_unit > 0) WRITE (io_unit, '(A,1X,A,1X,I6,1X,ES16.7)') "UPF native integrals:", TRIM(filename), &
            projectors%nproj, MAXVAL(ABS(integrals))
         DEALLOCATE (integrals)
      END IF
      CALL atom_release_upf(pot)
   END DO
   CALL cp2k_finalize()
END PROGRAM upf_projector_integrals_unittest
