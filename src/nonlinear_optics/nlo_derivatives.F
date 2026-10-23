!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \brief Matrix algebra for analytic k derivatives through fifth order.
!> \note Each array stores normalized derivatives X_α=(∂_k^α X)/α!,
!>       with α=(α_x,α_y,α_z), α!=α_x!α_y!α_z! and |α|<=5.
!>       Examples: X_(1,1,0)=∂_kx ∂_ky X and X_(2,0,0)=(∂_kx² X)/2.
!>       Products are formed by the multivariable product rule; no finite differences
!>       are used. The overlap inverse square root X=S^(-1/2) satisfies X S X=I.
!>       Its derivative of each order is obtained from lower orders by a Sylvester solve.
! **************************************************************************************************
MODULE nlo_derivatives
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo_derivatives'

   INTEGER, PARAMETER, PUBLIC                          :: nlo_max_derivatives_order = 5, &
                                                          nlo_n_derivatives = 56
   PUBLIC :: nlo_derivative_indices, nlo_derivative_product, nlo_lowdin_derivatives, &
             nlo_orthonormal_operator_derivatives, &
             nlo_equilibrium_density_derivatives, &
             nlo_derivative_count, derivatives_index

CONTAINS

! **************************************************************************************************
!> \brief Number of Cartesian Taylor coefficients through the given total grade.
!> \param grade Highest total derivative grade, from zero through five.
!> \return Number of Cartesian coefficients through the requested grade.
! **************************************************************************************************
   INTEGER FUNCTION nlo_derivative_count(grade) RESULT(count)
      INTEGER, INTENT(IN)                                :: grade

      count = 0
      IF (grade < 0 .OR. grade > nlo_max_derivatives_order) THEN
         RETURN
      END IF
      count = (grade + 1)*(grade + 2)*(grade + 3)/6

   END FUNCTION nlo_derivative_count

! **************************************************************************************************
!> \brief Enumerate the 56 three-dimensional multi-indices through total order five.
!> \param alpha Cartesian derivative counts (α_x, α_y, α_z) for each Taylor coefficient.
! **************************************************************************************************
   SUBROUTINE nlo_derivative_indices(alpha)
      INTEGER, INTENT(OUT)                               :: alpha(3, nlo_n_derivatives)

      INTEGER                                            :: ia, ib, ij, iorder

      ij = 0
      DO iorder = 0, nlo_max_derivatives_order
         DO ia = iorder, 0, -1
            DO ib = iorder - ia, 0, -1
               ij = ij + 1
               alpha(:, ij) = [ia, ib, iorder - ia - ib]
            END DO
         END DO
      END DO

   END SUBROUTINE nlo_derivative_indices

! **************************************************************************************************
!> \brief Apply the product rule to two arrays of normalized matrix k derivatives.
!> \param alpha Cartesian derivative counts (α_x, α_y, α_z) for each Taylor coefficient.
!> \param left Normalized Taylor coefficients of the left matrix factor.
!> \param right Normalized Taylor coefficients of the right matrix factor.
!> \param product Normalized Taylor coefficients of the matrix product.
!> \param active_grade Highest total derivative grade to compute; omit for all five grades.
! **************************************************************************************************
   SUBROUTINE nlo_derivative_product(alpha, left, right, product, active_grade)
      INTEGER, INTENT(IN)                                :: alpha(3, nlo_n_derivatives)
      COMPLEX(dp), INTENT(IN)                            :: left(:, :, :), right(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: product(:, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: active_grade

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_derivative_product'

      COMPLEX(dp), ALLOCATABLE                           :: term(:, :)
      INTEGER                                            :: handle, ib, ij, ik, nj

      CALL timeset(routineN, handle)

      nj = nlo_n_derivatives
      IF (PRESENT(active_grade)) nj = nlo_derivative_count(active_grade)
      product = 0.0_dp
      ALLOCATE (term(SIZE(left, 1), SIZE(right, 2)))
      DO ij = 1, nj
         DO ib = 1, ij
            ik = derivatives_index(alpha, alpha(:, ij) - alpha(:, ib))
            IF (ik == 0) CYCLE
            term(:, :) = MATMUL(left(:, :, ib), right(:, :, ik))
            PRODUCT(:, :, ij) = PRODUCT(:, :, ij) + term
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE nlo_derivative_product

! **************************************************************************************************
!> \brief Solve X(k) S(k) X(k)=I analytically through fifth Cartesian order.
!> \param overlap Normalized Taylor coefficients of S; leading coefficient Hermitian positive
!> \param inverse_root Normalized Taylor coefficients of S**(-1/2)
!> \param active_grade Highest total derivative grade to compute; omit for all five grades.
! **************************************************************************************************
   SUBROUTINE nlo_lowdin_derivatives(overlap, inverse_root, active_grade)
      COMPLEX(dp), INTENT(IN)                            :: overlap(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: inverse_root(:, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: active_grade

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_lowdin_derivatives'

      COMPLEX(dp), ALLOCATABLE                           :: matrix_work(:, :), metric_u(:, :), &
                                                            rhs(:, :), work(:)
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: matrix_product
      INTEGER                                            :: handle, ib, ic, id, ij, info, n, nj
      INTEGER, DIMENSION(3, nlo_n_derivatives)           :: alpha
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), ALLOCATABLE                         :: metric_e(:), work_e(:)

      CALL timeset(routineN, handle)

      n = SIZE(overlap, 1)
      inverse_root = 0.0_dp
      nj = nlo_n_derivatives
      IF (PRESENT(active_grade)) nj = nlo_derivative_count(active_grade)
      CPASSERT(.NOT. (n < 1 .OR. ANY(SHAPE(overlap) /= [n, n, nj])))
      CPASSERT(.NOT. (ANY(SHAPE(inverse_root) /= [n, n, nj])))
      CPASSERT(nj >= 1)
      DO ij = 1, nj
         scale = MAX(1.0_dp, MAXVAL(ABS(overlap(:, :, ij))))
         IF (.NOT. ALL(ieee_is_finite(REAL(overlap(:, :, ij), KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(overlap(:, :, ij)))) .OR. &
             MAXVAL(ABS(overlap(:, :, ij) - &
                        CONJG(TRANSPOSE(overlap(:, :, ij))))) > 1.e-10_dp*scale) THEN
            CPABORT('NONLINEAR_OPTICS Löwdin derivative check failed.')
         END IF
      END DO
      ALLOCATE (metric_u(n, n), rhs(n, n), matrix_work(n, n), &
                matrix_product(n, n), work(MAX(1, 2*n)), &
                metric_e(n), work_e(MAX(1, 3*n)))
      metric_u(:, :) = overlap(:, :, 1)

      CALL zheev('V', 'U', n, metric_u, n, metric_e, work, SIZE(work), work_e, info)
      IF (info /= 0) CPABORT('NONLINEAR_OPTICS Löwdin overlap eigensolver failed.')

      IF (MAXVAL(metric_e) <= 0.0_dp .OR. MINVAL(metric_e) <= 1.e-10_dp*MAXVAL(metric_e)) THEN
         CPABORT('NONLINEAR_OPTICS Löwdin derivative check failed.')
      END IF
      rhs = 0.0_dp
      DO ib = 1, n
         rhs(ib, ib) = 1.0_dp/SQRT(metric_e(ib))
      END DO
      matrix_work(:, :) = MATMUL(rhs, CONJG(TRANSPOSE(metric_u)))
      inverse_root(:, :, 1) = MATMUL(metric_u, matrix_work)
      CALL nlo_derivative_indices(alpha)

      DO ij = 2, nj
         rhs = 0.0_dp
         DO ib = 1, ij
            DO ic = 1, ij
               id = derivatives_index(alpha, alpha(:, ij) - alpha(:, ib) - alpha(:, ic))
               IF (id == 0) CYCLE
               IF (ib == ij .OR. id == ij) CYCLE
               matrix_work(:, :) = MATMUL(overlap(:, :, ic), inverse_root(:, :, id))
               matrix_product(:, :) = MATMUL(inverse_root(:, :, ib), matrix_work)
               rhs(:, :) = rhs - matrix_product
            END DO
         END DO
         matrix_work(:, :) = MATMUL(rhs, metric_u)
         rhs(:, :) = MATMUL(CONJG(TRANSPOSE(metric_u)), matrix_work)
         DO ib = 1, n
            DO ic = 1, n
               rhs(ib, ic) = rhs(ib, ic)/(SQRT(metric_e(ib)) + SQRT(metric_e(ic)))
            END DO
         END DO
         matrix_work(:, :) = MATMUL(rhs, CONJG(TRANSPOSE(metric_u)))
         inverse_root(:, :, ij) = MATMUL(metric_u, matrix_work)
         scale = MAX(1.0_dp, MAXVAL(ABS(inverse_root(:, :, ij))))
         IF (MAXVAL(ABS(inverse_root(:, :, ij) - &
                        CONJG(TRANSPOSE(inverse_root(:, :, ij))))) > 1.e-9_dp*scale) THEN
            CPABORT('NONLINEAR_OPTICS Löwdin derivative check failed.')
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE nlo_lowdin_derivatives

! **************************************************************************************************
!> \brief Express H and the AO position matrix in an orthonormal AO basis.
!> \param overlap Normalized Taylor coefficients of the AO overlap S(k).
!> \param hamiltonian Normalized Taylor coefficients of the AO Hamiltonian H(k).
!> \param moment Normalized Taylor coefficients of the AO position matrices d_a(k).
!> \param orth_h Normalized Taylor coefficients of the orthogonal Hamiltonian.
!> \param connection Normalized Taylor coefficients of the orthogonal position connection.
!> \param active_grade Highest total derivative grade to compute; omit for all five grades.
!> \note AO basis functions overlap, so S(k) is not the identity. The matrix
!>       X(k)=S(k)^(-1/2) makes their overlap equal to the identity.
!>       This routine differentiates H^orth=X H X and the position matrix
!>       A_a^orth=X d_a X-i (dX/dk_a) S X, where d_a is the AO position matrix.
!>       It does not select occupied bands or build the equilibrium density rho0.
!>       A_a^orth is needed through order four; its order-five slots remain zero.
! **************************************************************************************************
   SUBROUTINE nlo_orthonormal_operator_derivatives(overlap, hamiltonian, moment, orth_h, &
                                                   connection, active_grade)
      COMPLEX(dp), DIMENSION(:, :, :), INTENT(IN)        :: overlap, hamiltonian
      COMPLEX(dp), INTENT(IN)                            :: moment(:, :, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: orth_h(:, :, :), connection(:, :, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: active_grade

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_orthonormal_operator_derivatives'

      COMPLEX(dp), ALLOCATABLE                           :: derivative_x(:, :, :), &
                                                            inverse_root(:, :, :), &
                                                            matrix_derivatives(:, :, :), &
                                                            moment_derivatives(:, :, :)
      INTEGER                                            :: alpha(3, nlo_n_derivatives), grade, &
                                                            handle, ia, ij, ik, n, na, nh, &
                                                            sought(3)
      REAL(KIND=dp)                                      :: scale

      CALL timeset(routineN, handle)

      ! AO overlap and Hamiltonian.
      n = SIZE(overlap, 1)
      orth_h = 0.0_dp
      connection = 0.0_dp
      grade = nlo_max_derivatives_order
      IF (PRESENT(active_grade)) grade = active_grade
      nh = nlo_derivative_count(grade)
      na = nlo_derivative_count(grade - 1)
      CPASSERT(grade >= 2 .AND. nh >= 1)
      CPASSERT(ALL(SHAPE(overlap) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(hamiltonian) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(moment) == [n, n, 3, na]))
      CPASSERT(ALL(SHAPE(orth_h) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(connection) == [n, n, 3, na]))
      IF (.NOT. ALL(ieee_is_finite(REAL(hamiltonian, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(hamiltonian))) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(moment, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(moment)))) THEN
         CPABORT('NONLINEAR_OPTICS AO operator derivatives are inconsistent.')
      END IF
      CALL nlo_derivative_indices(alpha)
      DO ij = 1, nh
         scale = MAX(1.0_dp, MAXVAL(ABS(hamiltonian(:, :, ij))))
         IF (MAXVAL(ABS(hamiltonian(:, :, ij) - &
                        CONJG(TRANSPOSE(hamiltonian(:, :, ij))))) > 1.e-9_dp*scale) THEN
            CPABORT('NONLINEAR_OPTICS AO operator derivatives are inconsistent.')
         END IF
         IF (SUM(alpha(:, ij)) >= grade) CYCLE
         DO ia = 1, 3
            sought = alpha(:, ij)
            sought(ia) = sought(ia) + 1
            ik = derivatives_index(alpha, sought)
            scale = MAX(1.0_dp, MAXVAL(ABS(moment(:, :, ia, ij))), &
                        MAXVAL(ABS(overlap(:, :, ik))))
            IF (MAXVAL(ABS(moment(:, :, ia, ij) - &
                           CONJG(TRANSPOSE(moment(:, :, ia, ij))) + &
                           gaussi*REAL(sought(ia), KIND=dp)*overlap(:, :, ik))) > &
                1.e-9_dp*scale) THEN
               CPABORT('NONLINEAR_OPTICS AO operator derivatives are inconsistent.')
            END IF
         END DO
      END DO
      ALLOCATE (inverse_root(n, n, nh), derivative_x(n, n, na), &
                matrix_derivatives(n, n, nh), moment_derivatives(n, n, na))

      ! X = S^(-1/2), with X S X = 1 at every derivative order.
      CALL nlo_lowdin_derivatives(overlap, inverse_root, active_grade=grade)

      ! H^o = X H X; coefficient products use the multivariable product rule.
      CALL nlo_derivative_product(alpha, inverse_root, hamiltonian, matrix_derivatives, &
                                  active_grade=grade)
      CALL nlo_derivative_product(alpha, matrix_derivatives, inverse_root, orth_h, &
                                  active_grade=grade)

      DO ia = 1, 3
         derivative_x = 0.0_dp
         DO ij = 1, na
            IF (SUM(alpha(:, ij)) >= grade) CYCLE
            sought = alpha(:, ij)
            sought(ia) = sought(ia) + 1
            ik = derivatives_index(alpha, sought)
            derivative_x(:, :, ij) = REAL(alpha(ia, ij) + 1, KIND=dp)*inverse_root(:, :, ik)
         END DO
         CALL nlo_derivative_product(alpha, derivative_x, overlap, matrix_derivatives, &
                                     active_grade=grade - 1)
         CALL nlo_derivative_product(alpha, matrix_derivatives, inverse_root, &
                                     moment_derivatives, active_grade=grade - 1)
         connection(:, :, ia, :) = -gaussi*moment_derivatives
         CALL nlo_derivative_product(alpha, inverse_root, moment(:, :, ia, :), &
                                     matrix_derivatives, active_grade=grade - 1)
         CALL nlo_derivative_product(alpha, matrix_derivatives, inverse_root, &
                                     moment_derivatives, active_grade=grade - 1)
         connection(:, :, ia, :) = connection(:, :, ia, :) + moment_derivatives
         DO ij = 1, na
            scale = MAX(1.0_dp, MAXVAL(ABS(connection(:, :, ia, ij))))
            IF (MAXVAL(ABS(connection(:, :, ia, ij) - &
                           CONJG(TRANSPOSE(connection(:, :, ia, ij))))) > 1.e-9_dp*scale) THEN
               CPABORT('NONLINEAR_OPTICS AO operator derivatives are inconsistent.')
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE nlo_orthonormal_operator_derivatives

! **************************************************************************************************
!> \brief Derivatives of the equilibrium density matrix for filled Kohn-Sham bands.
!> \param energies Band energies ε_m at the central k point.
!> \param occupied Unit occupation flags for the bands at the central k point.
!> \param hamiltonian Normalized Taylor coefficients of the AO Hamiltonian H(k).
!> \param rho0 Normalized Taylor coefficients of the equilibrium density matrix.
!> \param active_grade Highest total derivative grade to compute; omit for all five grades.
!> \note At the central k point, rho0_mn = occupied_m delta_mn with occupied_m=0 or 1.
!>       At nearby k, rho0(k) = sum_(m occupied) |u_m(k)><u_m(k)|.
!>       This is a spectral projector because rho0**2 = rho0; it is also the
!>       equilibrium density matrix for unit band occupation. The actual spin
!>       occupation factor is applied later during Brillouin-zone integration.
!>       The routine obtains analytic k derivatives from H rho0 = rho0 H
!>       and rho0**2 = rho0, without differentiating individual eigenvectors.
!>       It does not calculate the field-induced rho^(1), rho^(2), etc.
! **************************************************************************************************
   SUBROUTINE nlo_equilibrium_density_derivatives(energies, occupied, hamiltonian, rho0, &
                                                  active_grade)
      REAL(KIND=dp), INTENT(IN)                          :: energies(:)
      INTEGER, INTENT(IN)                                :: occupied(:)
      COMPLEX(dp), INTENT(IN)                            :: hamiltonian(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: rho0(:, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: active_grade

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_equilibrium_density_derivatives'

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :), forcing(:, :), &
                                                            matrix_term(:, :)
      INTEGER                                            :: handle, ib, ic, ij, ik, n, nj
      INTEGER, DIMENSION(3, nlo_n_derivatives)           :: alpha
      REAL(KIND=dp)                                      :: gap, scale, symmetry_error

      CALL timeset(routineN, handle)

      n = SIZE(energies)
      rho0 = 0.0_dp
      CPASSERT(.NOT. (n < 2 .OR. SIZE(occupied) /= n))
      nj = nlo_n_derivatives
      IF (PRESENT(active_grade)) nj = nlo_derivative_count(active_grade)
      CPASSERT(ALL(SHAPE(hamiltonian) == [n, n, nj]))
      CPASSERT(ALL(SHAPE(rho0) == [n, n, nj]))
      CPASSERT(.NOT. (ANY(occupied < 0) .OR. ANY(occupied > 1)))
      CPASSERT(.NOT. (ALL(occupied == 0) .OR. ALL(occupied == 1)))
      IF (.NOT. ALL(ieee_is_finite(energies)) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(hamiltonian, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(hamiltonian)))) THEN
         CPABORT('NONLINEAR_OPTICS equilibrium density input is not finite.')
      END IF
      scale = MAX(1.0_dp, MAXVAL(ABS(energies)))
      ! The occupied and empty energies must remain separated at this k point.
      gap = MINVAL(energies, MASK=occupied == 0) - MAXVAL(energies, MASK=occupied == 1)
      IF (gap <= 1.e-10_dp*scale) THEN
         CPABORT('NONLINEAR_OPTICS requires a positive occupied-empty gap.')
      END IF
      ALLOCATE (block(n, n), forcing(n, n), matrix_term(n, n))
      matrix_term(:, :) = hamiltonian(:, :, 1)
      DO ib = 1, n
         matrix_term(ib, ib) = matrix_term(ib, ib) - energies(ib)
      END DO
      IF (MAXVAL(ABS(matrix_term)) > 1.e-9_dp*scale) THEN
         CPABORT('NONLINEAR_OPTICS Hamiltonian is not diagonal in the band basis.')
      END IF
      CALL nlo_derivative_indices(alpha)
      ! Zeroth Taylor coefficient: the 0/1 equilibrium density in the band basis.
      DO ib = 1, n
         rho0(ib, ib, 1) = occupied(ib)
      END DO
      ! Differentiate [H,rho0]=0 and rho0**2=rho0 order by order in k.
      DO ij = 2, nj
         block = 0.0_dp
         forcing = 0.0_dp
         DO ik = 2, ij
            ic = derivatives_index(alpha, alpha(:, ij) - alpha(:, ik))
            IF (ic == 0) CYCLE
            matrix_term(:, :) = MATMUL(hamiltonian(:, :, ik), rho0(:, :, ic))
            forcing(:, :) = forcing + matrix_term
            matrix_term(:, :) = MATMUL(rho0(:, :, ic), hamiltonian(:, :, ik))
            forcing(:, :) = forcing - matrix_term
            IF (ic == 1) CYCLE
            matrix_term(:, :) = MATMUL(rho0(:, :, ik), rho0(:, :, ic))
            block(:, :) = block + matrix_term
         END DO
         DO ic = 1, n
            DO ib = 1, n
               ! Cross-gap blocks use the energy difference; same-group blocks
               ! follow from rho0**2 = rho0.
               IF (occupied(ib) /= occupied(ic)) THEN
                  rho0(ib, ic, ij) = -forcing(ib, ic)/(energies(ib) - energies(ic))
               ELSE IF (occupied(ib) == 1) THEN
                  rho0(ib, ic, ij) = -block(ib, ic)
               ELSE
                  rho0(ib, ic, ij) = block(ib, ic)
               END IF
            END DO
         END DO
      END DO
      symmetry_error = 0.0_dp
      DO ij = 1, nj
         symmetry_error = MAX(symmetry_error, &
                              MAXVAL(ABS(rho0(:, :, ij) - CONJG(TRANSPOSE(rho0(:, :, ij)))))/ &
                              MAX(1.0_dp, MAXVAL(ABS(rho0(:, :, ij)))))
      END DO
      IF (.NOT. ALL(ieee_is_finite(REAL(rho0, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(rho0))) .OR. symmetry_error > 1.e-8_dp) THEN
         CPABORT('NONLINEAR_OPTICS equilibrium density derivatives are invalid.')
      END IF

      CALL timestop(handle)

   END SUBROUTINE nlo_equilibrium_density_derivatives

! **************************************************************************************************
!> \brief Find the Taylor coefficient for the requested (x,y,z) derivative counts.
!> \param alpha Cartesian derivative counts (α_x, α_y, α_z) for each Taylor coefficient.
!> \param sought Requested Cartesian derivative counts; negative counts return zero.
!> \return Its index in alpha, or zero when any count is negative or absent.
! **************************************************************************************************
   INTEGER FUNCTION derivatives_index(alpha, sought) RESULT(index)
      INTEGER, INTENT(IN)                                :: alpha(3, nlo_n_derivatives), sought(3)

      INTEGER                                            :: ij

      index = 0
      IF (ANY(sought < 0)) RETURN
      DO ij = 1, nlo_n_derivatives
         IF (ALL(alpha(:, ij) == sought)) THEN
            index = ij
            RETURN
         END IF
      END DO

   END FUNCTION derivatives_index

END MODULE nlo_derivatives
