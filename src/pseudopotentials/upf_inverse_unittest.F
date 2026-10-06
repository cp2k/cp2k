!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse reconstruction identities, conditioning and analytic field/operator responses.
! **************************************************************************************************
PROGRAM upf_inverse_unittest
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_evaluate,&
                                              gaussian_basis_fields_periodic
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              ncoset
   USE upf_inverse,                     ONLY: &
        upf_inverse_atomic_residual, upf_inverse_atomic_type, upf_inverse_basis, &
        upf_inverse_basis_response, upf_inverse_init, upf_inverse_operator, &
        upf_inverse_operator_response, upf_inverse_project, upf_inverse_project_response, &
        upf_inverse_type

   IMPLICIT NONE
   REAL(dp), PARAMETER :: step = 1.e-5_dp
   REAL(dp) :: max_error = 0.0_dp
   CALL check_inverse()
   CALL check_atomic_residual()
   CALL check_gaussian_fields()
   WRITE (*, '(A)') 'upf_inverse_unittest PASS'
CONTAINS
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_inverse()
      INTEGER, PARAMETER                                 :: nb = 2, nf = 5, np = 3

      COMPLEX(KIND=dp) :: actual(nb, nb), b(np, nb), bd(np, np), chi(nf, nb), dchi(nf, nb), &
         dd(np, np), ddd(np, np), ddelta(nf, np), ddg(np, nb), delta(nf, np), dg(np, nb), &
         dgd(nb, np), dgg(nb, nb), donsite(np, np), dop(nf, nf), dprojector(np, nf), ds(nf, nb), &
         dw(np, nb), empty_b(0, nb), empty_delta(nf, 0), empty_q(0, 0), empty_w(0, nb), &
         gd(nb, np), gg(nb, nb), images(nf, nb, 2), onsite(np, np), op(nf, nf), projected(np, nb), &
         projector(np, nf), reference(nb, nb), response(nb, nb), sampled(nb, nb, 2), &
         shifted_chi(nf, nb), shifted_delta(nf, np), shifted_op(nf, nf), shifted_p(np, nf)
      COMPLEX(KIND=dp) :: shifted_q(np, np), soft(nf, nb), w(np, nb)
      INTEGER                                            :: i, info, j, phase, side
      REAL(KIND=dp)                                      :: shift
      TYPE(upf_inverse_type)                             :: inverse

      DO phase = 1, 2
         DO j = 1, nb
            DO i = 1, nf
               chi(i, j) = CMPLX(SIN(REAL(i + 2*j, dp)), (phase - 1)*COS(REAL(i*j, dp)), dp)
               dchi(i, j) = 0.07_dp*CMPLX(COS(REAL(i + j, dp)), (phase - 1)*SIN(REAL(2*i + j, dp)), dp)
            END DO
         END DO
         DO j = 1, np
            DO i = 1, nf
               delta(i, j) = 0.12_dp*CMPLX(COS(REAL(i*j, dp)), (phase - 1)*SIN(REAL(i + j, dp)), dp)
               ddelta(i, j) = 0.03_dp*CMPLX(SIN(REAL(i + 2*j, dp)), (phase - 1)*COS(REAL(i*j, dp)), dp)
               projector(j, i) = 0.15_dp*CMPLX(SIN(REAL(i + j, dp)), (phase - 1)*COS(REAL(2*i + j, dp)), dp)
               dprojector(j, i) = 0.02_dp*CMPLX(COS(REAL(i + j, dp)), (phase - 1)*SIN(REAL(i*j, dp)), dp)
            END DO
         END DO
         ! A non-Hermitian operator exposes an accidental conjugation/symmetrization.
         DO j = 1, nf
            DO i = 1, nf
               op(i, j) = CMPLX(0.1_dp*i - 0.05_dp*j, (phase - 1)*0.03_dp*(i + j), dp)
               dop(i, j) = CMPLX(0.02_dp*(i + j), (phase - 1)*0.01_dp*(i - j), dp)
            END DO
         END DO
         onsite = op(1:np, 1:np)*0.3_dp
         donsite = dop(1:np, 1:np)*0.4_dp
         bd = MATMUL(projector, delta)
         b = MATMUL(projector, chi)
         CALL upf_inverse_init(bd, inverse, info)
         IF (info /= 0 .OR. .NOT. inverse%ready) ERROR STOP 'Inverse factorization failed'
         CALL upf_inverse_project(inverse, b, w)
         CALL upf_inverse_basis(chi, delta, w, soft)
         projected = MATMUL(projector, soft)
         CALL check(MAXVAL(ABS(projected - w)), 'inverse projector identity')
         CALL check(MAXVAL(ABS(soft + MATMUL(delta, projected) - chi)), 'full reconstruction round trip')
         CALL upf_inverse_project_response(inverse, w, MATMUL(dprojector, chi) + MATMUL(projector, dchi), &
                                           MATMUL(dprojector, delta) + MATMUL(projector, ddelta), dw)
         CALL upf_inverse_basis_response(dchi, delta, ddelta, w, dw, ds)
         gg = MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(op, chi))
         gd = MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(op, delta))
         dg = MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(op, chi))
         dd = MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(op, delta))
         CALL upf_inverse_operator(gg, gd, dg, dd, onsite, w, actual)
         reference = MATMUL(CONJG(TRANSPOSE(soft)), MATMUL(op, soft)) + &
                     MATMUL(CONJG(TRANSPOSE(projected)), MATMUL(onsite, projected))
         CALL check(MAXVAL(ABS(actual - reference)), 'inverse operator versus direct fields')
         dgg = MATMUL(CONJG(TRANSPOSE(dchi)), MATMUL(op, chi)) + &
               MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(dop, chi) + MATMUL(op, dchi))
         dgd = MATMUL(CONJG(TRANSPOSE(dchi)), MATMUL(op, delta)) + &
               MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(dop, delta) + MATMUL(op, ddelta))
         ddg = MATMUL(CONJG(TRANSPOSE(ddelta)), MATMUL(op, chi)) + &
               MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(dop, chi) + MATMUL(op, dchi))
         ddd = MATMUL(CONJG(TRANSPOSE(ddelta)), MATMUL(op, delta)) + &
               MATMUL(CONJG(TRANSPOSE(delta)), MATMUL(dop, delta) + MATMUL(op, ddelta))
         CALL upf_inverse_operator_response(dgg, gd, dg, dd, onsite, w, dgd, ddg, ddd, donsite, dw, response)
         DO side = 1, 2
            shift = (2*side - 3)*step
            shifted_chi = chi + shift*dchi
            shifted_delta = delta + shift*ddelta
            shifted_p = projector + shift*dprojector
            shifted_op = op + shift*dop
            shifted_q = onsite + shift*donsite
            CALL upf_inverse_init(MATMUL(shifted_p, shifted_delta), inverse, info)
            IF (info /= 0) ERROR STOP 'Displaced inverse factorization failed'
            CALL upf_inverse_project(inverse, MATMUL(shifted_p, shifted_chi), w)
            CALL upf_inverse_basis(shifted_chi, shifted_delta, w, images(:, :, side))
            projected = MATMUL(shifted_p, images(:, :, side))
            sampled(:, :, side) = MATMUL(CONJG(TRANSPOSE(images(:, :, side))), &
                                         MATMUL(shifted_op, images(:, :, side))) + &
                                  MATMUL(CONJG(TRANSPOSE(projected)), MATMUL(shifted_q, projected))
         END DO
         CALL check(MAXVAL(ABS(ds - (images(:, :, 2) - images(:, :, 1))/(2*step))), 'inverse basis response', 1.e-9_dp)
         CALL check(MAXVAL(ABS(response - (sampled(:, :, 2) - sampled(:, :, 1))/(2*step))), &
                    'inverse operator response', 1.e-9_dp)
      END DO

      CALL upf_inverse_init(empty_q, inverse, info)
      IF (info /= 0 .OR. .NOT. inverse%ready) ERROR STOP 'Empty inverse failed'
      CALL upf_inverse_project(inverse, empty_b, empty_w)
      CALL upf_inverse_basis(chi, empty_delta, empty_w, soft)
      CALL check(MAXVAL(ABS(soft - chi)), 'identity without reconstruction projectors')
      bd = CMPLX(0.0_dp, 0.0_dp, dp)
      bd(2, 2) = -1.0_dp
      CALL upf_inverse_init(bd, inverse, info)
      IF (info /= 2 .OR. inverse%ready) ERROR STOP 'Singular inverse accepted'
      bd(2, 2) = -1.0_dp + 1.e-12_dp
      CALL upf_inverse_init(bd, inverse, info, minimum_rcond=1.e-10_dp)
      IF (info /= np + 1 .OR. inverse%ready) ERROR STOP 'Ill-conditioned inverse accepted'
      bd(1, 1) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
      CALL upf_inverse_init(bd, inverse, info)
      IF (info /= np + 2 .OR. inverse%ready) ERROR STOP 'Nonfinite inverse accepted'
   END SUBROUTINE check_inverse

! **************************************************************************************************
!> \brief Verify the one-center reduction against an explicit finite transformed operator.
!>        Nonzero closure residuals and off-diagonal channel couplings must remain present.
! **************************************************************************************************
   SUBROUTINE check_atomic_residual()
      COMPLEX(dp)                                        :: b(2, 2), chi(3, 2), coupling(2, 2), &
                                                            metric(2, 2), soft(3, 2), w(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: rs(:, :), rt(:, :)
      INTEGER                                            :: info
      REAL(dp)                                           :: delta(3, 2), hard(3, 2), ps(3, 2), &
                                                            q(2, 2), t(3, 3)
      TYPE(upf_inverse_atomic_type)                      :: atomic

      ps = RESHAPE([0.2_dp, 0.5_dp, -0.3_dp, 0.1_dp, -0.4_dp, 0.7_dp], [3, 2])
      delta = RESHAPE([0.3_dp, -0.1_dp, 0.2_dp, -0.2_dp, 0.4_dp, 0.1_dp], [3, 2])
      hard = ps + delta
      t = RESHAPE([1.0_dp, 0.2_dp, 0.1_dp, 0.2_dp, 2.0_dp, -0.3_dp, 0.1_dp, -0.3_dp, 1.5_dp], [3, 3])
      q = RESHAPE([0.4_dp, 0.1_dp, 0.1_dp, 0.6_dp], [2, 2])
      chi = CMPLX(ps, delta, dp)
      coupling = CMPLX(MATMUL(TRANSPOSE(ps), delta), 0.0_dp, dp)
      CALL upf_inverse_init(coupling, atomic%inverse, info)
      IF (info /= 0) ERROR STOP 'Atomic inverse failed'
      b = MATMUL(TRANSPOSE(ps), chi)
      CALL upf_inverse_project(atomic%inverse, b, w)
      CALL upf_inverse_basis(chi, CMPLX(delta, 0.0_dp, dp), w, soft)
      atomic%hard%nproj = 2
      ALLOCATE (atomic%kinetic, SOURCE=MATMUL(TRANSPOSE(hard), MATMUL(t, hard)) - MATMUL(TRANSPOSE(ps), MATMUL(t, ps)))
      ALLOCATE (atomic%difference_kinetic, SOURCE=MATMUL(TRANSPOSE(delta), MATMUL(t, delta)))
      ALLOCATE (atomic%wave_overlap, SOURCE=MATMUL(TRANSPOSE(hard), hard) - MATMUL(TRANSPOSE(ps), ps))
      ALLOCATE (atomic%difference_overlap, SOURCE=MATMUL(TRANSPOSE(delta), delta))
      CALL upf_inverse_atomic_residual(atomic, w, MATMUL(TRANSPOSE(delta), chi), MATMUL(TRANSPOSE(delta), MATMUL(t, chi)), rt, rs)
      CALL check(MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(t, chi)) + rt - &
                            MATMUL(CONJG(TRANSPOSE(w)), MATMUL(atomic%kinetic, w)) - &
                            MATMUL(CONJG(TRANSPOSE(soft)), MATMUL(t, soft)))), 'atomic kinetic closure residual')
      metric = MATMUL(CONJG(TRANSPOSE(chi)), chi) + rs + &
               MATMUL(CONJG(TRANSPOSE(w)), MATMUL(q - atomic%wave_overlap, w))
      CALL check(MAXVAL(ABS(metric - MATMUL(CONJG(TRANSPOSE(soft)), soft) - &
                            MATMUL(CONJG(TRANSPOSE(w)), MATMUL(q, w)))), 'atomic metric closure residual')
      IF (MAXVAL(ABS(rt)) < 1.e-3_dp .OR. MAXVAL(ABS(rs)) < 1.e-3_dp) ERROR STOP 'Trivial closure test'
   END SUBROUTINE check_atomic_residual

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_gaussian_fields()
      INTEGER, PARAMETER                                 :: nb = 11, nf = 5

      COMPLEX(dp)                                        :: dg(nf, nb, 3), dv(nf, nb), g(nf, nb, 3), &
                                                            gss(nf, nb, 3, 2), other(nf, nb), &
                                                            v(nf, nb), vs(nf, nb, 2)
      INTEGER                                            :: axis, ico, image, images(3, 3), row, &
                                                            side, translated_images(3, 3)
      REAL(dp) :: bloch(3), cell(3, 3), center(3), dc(3, 3), dcenter(3), dk(3), dpnt(3, nf), &
         gradient(nf, nb, 3), gs(nf, nb, 3, 2), hh(nf, nb, 3, 3), samples(nf, nb, 2), shift, &
         value(nf, nb), x(3, nf), xp(3, nf)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      NULLIFY (basis)
      CALL init_orbital_pointers(3)
      CALL allocate_gto_basis_set(basis)
      basis%nset = 2
      basis%nsgf = nb
      ALLOCATE (basis%lmax(2), basis%npgf(2), basis%nshell(2), basis%l(2, 2), &
                basis%first_sgf(2, 2), basis%last_sgf(2, 2), basis%zet(2, 2), basis%sphi(40, nb))
      basis%lmax = [3, 1]
      basis%npgf = [2, 1]
      basis%nshell = [2, 1]
      basis%l = 0
      basis%l(:, 1) = [0, 3]
      basis%l(1, 2) = 1
      basis%first_sgf = 0
      basis%last_sgf = 0
      basis%first_sgf(:, 1) = [1, 2]
      basis%last_sgf(:, 1) = [1, 8]
      basis%first_sgf(1, 2) = 9
      basis%last_sgf(1, 2) = 11
      basis%zet = 0.0_dp
      basis%zet(:, 1) = [.4_dp, .9_dp]
      basis%zet(1, 2) = .7_dp
      basis%sphi = 0.0_dp
      basis%sphi(1, 1) = .6_dp
      basis%sphi(21, 1) = -.2_dp
      DO row = 2, 8
         DO ico = ncoset(2) + 1, ncoset(3)
            basis%sphi(ico, row) = .03_dp*SIN(REAL(ico*row, dp))
            basis%sphi(ico + 20, row) = -.02_dp*COS(REAL(ico + row, dp))
         END DO
      END DO
      DO row = 9, nb
         DO ico = 2, 4
            basis%sphi(ico, row) = .1_dp*SIN(REAL(ico*row, dp))
         END DO
      END DO
      x(:, 1) = 0.0_dp
      x(:, 2) = [.3_dp, 0.0_dp, 0.0_dp]
      x(:, 3) = [0.0_dp, -.7_dp, 0.0_dp]
      x(:, 4) = [.23_dp, .51_dp, -.4_dp]
      x(:, 5) = [-.6_dp, .2_dp, .8_dp]
      CALL gaussian_basis_fields_evaluate(basis, x, value, gradient, hh)
      CALL check(MAXVAL(ABS(value(:, 1) - (.6_dp*EXP(-.4_dp*SUM(x*x, DIM=1)) &
                                           - .2_dp*EXP(-.9_dp*SUM(x*x, DIM=1))))), 'native signed s contraction')
      DO axis = 1, 3
         DO side = 1, 2
            xp = x
            xp(axis, :) = xp(axis, :) + (2*side - 3)*step
            CALL gaussian_basis_fields_evaluate(basis, xp, samples(:, :, side), gs(:, :, :, side))
         END DO
         CALL check(MAXVAL(ABS(gradient(:, :, axis) - (samples(:, :, 2) - samples(:, :, 1))/(2*step))), &
                    'contracted Gaussian point gradient', 2.e-8_dp)
         CALL check(MAXVAL(ABS(hh(:, :, :, axis) - (gs(:, :, :, 2) - gs(:, :, :, 1))/(2*step))), &
                    'contracted Gaussian origin/axis Hessian', 2.e-7_dp)
      END DO
      cell = RESHAPE([2.3_dp, .1_dp, .2_dp, .2_dp, 2.7_dp, -.1_dp, .1_dp, .3_dp, 2.1_dp], [3, 3])
      center = [.1_dp, -.2_dp, .3_dp]
      images = RESHAPE([0, 0, 0, 1, -1, 0, -1, 0, 1], [3, 3])
      bloch = [.23_dp, -.17_dp, .11_dp]
      dc = .03_dp*cell
      dcenter = [.07_dp, -.02_dp, .01_dp]
      dpnt = .02_dp*x
      dk = [.03_dp, .02_dp, -.01_dp]
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, bloch, v, g, dv, dg, dpnt, dcenter, dc, dk)
      DO side = 1, 2
         shift = (2*side - 3)*step
         CALL gaussian_basis_fields_periodic(basis, x + shift*dpnt, center + shift*dcenter, cell + shift*dc, &
                                             images, bloch + shift*dk, vs(:, :, side), gss(:, :, :, side))
      END DO
      CALL check(MAXVAL(ABS(dv - (vs(:, :, 2) - vs(:, :, 1))/(2*step))), 'native Gaussian Bloch joint response', 2.e-8_dp)
     CALL check(MAXVAL(ABS(dg - (gss(:, :, :, 2) - gss(:, :, :, 1))/(2*step))), 'native Gaussian Bloch gradient response', 2.e-7_dp)
      CALL gaussian_basis_fields_periodic(basis, x, center, cell, images, -bloch, other)
      CALL check(MAXVAL(ABS(other - CONJG(v))), 'native Gaussian Bloch time reversal')
      DO image = 1, 3
         translated_images(:, image) = images(:, image) - [1, 0, 0]
      END DO
      CALL gaussian_basis_fields_periodic(basis, x, center + cell(:, 1), cell, translated_images, bloch, other)
      CALL check(MAXVAL(ABS(other - EXP(CMPLX(0.0_dp, -2*pi*bloch(1), dp))*v)), 'native Gaussian lattice representative')
      CALL deallocate_gto_basis_set(basis)
   END SUBROUTINE check_gaussian_fields

! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param label ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check(value, label, tolerance)
      REAL(KIND=dp), INTENT(IN)                          :: value
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tolerance

      REAL(KIND=dp)                                      :: tol

      tol = 2.e-12_dp
      IF (PRESENT(tolerance)) tol = tolerance
      max_error = MAX(max_error, value)
      IF (value > tol) THEN
         WRITE (*, *) TRIM(label), value, tol
         ERROR STOP 'Inverse primitive mismatch'
      END IF
   END SUBROUTINE check
END PROGRAM upf_inverse_unittest
