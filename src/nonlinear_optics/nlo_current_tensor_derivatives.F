!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \brief Local full-tensor independent-particle current response through fifth field order.
!> \note This is a Bloch-point kernel. It does not integrate the Brillouin zone,
!>       convert current to susceptibility, or build AO derivatives from CP2K data.
! **************************************************************************************************
MODULE nlo_current
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi
   USE nlo_T1_T2_band_transition,       ONLY: nlo_band_transition_correction,&
                                              nlo_band_transition_prepare,&
                                              nlo_band_transition_rate,&
                                              nlo_band_transition_type
   USE nlo_derivatives,                 ONLY: derivatives_index,&
                                              nlo_derivative_count,&
                                              nlo_derivative_indices,&
                                              nlo_n_derivatives
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo_current'

   PUBLIC :: nlo_current_tensor_derivatives

CONTAINS

! **************************************************************************************************
!> \brief Compute intrinsically permutation-symmetric χ^(n) current at one k point.
!> \param order Number of incident fields, 2:5
!> \param omega Incident photon energies in Hartree
!> \param energies Central band energies in Hartree
!> \param rate1 Population damping energy in Hartree
!> \param rate2 Coherence damping energy in Hartree
!> \param transition_width Band-energy width of their smooth transition in Hartree
!> \param hamiltonian Orthogonal Bloch Hamiltonian normalized k derivatives through grade n
!> \param connection Orthogonal Bloch position connection derivatives through grade n-1
!> \param curvature Three position-curvature matrices for this k point and spin
!> \param rho0 Equilibrium density derivatives through grade n
!> \param velocity Physical central-band velocity in bohr/atomic time
!> \param current Complex current tensor, 3 by 3**order
!> \param elements Optional output/input selection, same shape and indexing as current.
!>        Only required subset responses are evaluated. Unselected outputs remain zero.
!> \note Label the n incident fields by j=1,...,n, each with frequency ω_j and
!>       Cartesian direction a_j. For every nonempty labeled subset S, form
!>       Ω_S=Σ_{j∈S}ω_j and the normalized density response
!>       ρ_S = |S|⁻¹ Σ_(j in S) R(Ω_S)L_(a_j)ρ_(S without j),
!>       where R L is the analytic normalized k-derivative step of nlo_density_step_derivatives.
!>       This average equals the sum of all ordered responses divided by n!.
!>       Only k-derivative grades n-|S| are retained for each subset. States
!>       with the same labeled subset and field directions are reused across
!>       all full-tensor components, retaining only two adjacent field orders
!>       in memory. The final current is -Tr[v_c ρ_{1,...,n}]
!>       minus n⁻¹ Σ_j Tr[F_(c,a_j) ρ_(S without j)], in atomic units.
!>       This follows from v_c(E)=i[H+E_a r_a,r_c]=v_c(0)+E_a F_ca,
!>       with [r_c,r_a]=iF_ca. The curvature term is always included. The second
!>       array index stores a_1,...,a_n as 1+Σ_j(a_j-1)3^(j-1).
! **************************************************************************************************
   SUBROUTINE nlo_current_tensor_derivatives(order, omega, energies, rate1, rate2, &
                                             transition_width, hamiltonian, &
                                             connection, curvature, rho0, velocity, current, &
                                             elements)
      INTEGER, INTENT(IN)                                :: order
      REAL(KIND=dp), INTENT(IN)                          :: omega(:), energies(:), rate1, rate2, &
                                                            transition_width
      COMPLEX(dp), INTENT(IN)                            :: hamiltonian(:, :, :), &
                                                            connection(:, :, :, :), &
                                                            curvature(:, :, :), rho0(:, :, :), &
                                                            velocity(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: current(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: elements(:, :)

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_current_tensor_derivatives'

      COMPLEX(dp)                                        :: curvature_term, extra(3, 243)
      COMPLEX(dp), ALLOCATABLE                           :: next_density(:, :, :), &
                                                            previous_density(:, :, :)
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :, :)       :: density
      INTEGER :: code, current_start, direction, field, full_mask, handle, i_code, i_field, &
         i_mask, i_radix, i_subset_code, index, iother, iother_direction, iother_position, &
         ioutput, mask, multiplicity, n, na, ncodes(0:31), ncoeff(0:31), nh, &
         offset_current(0:242, 0:31), offset_previous(0:242, 0:31), order_subset, position, &
         previous_code, previous_end, previous_mask, previous_start, radix, &
         representative_code(0:242, 0:31), total_coefficients
      LOGICAL                                            :: evaluate(0:242, 0:31), &
                                                            required(0:242, 0:31), selected(3, 243)
      REAL(KIND=dp)                                      :: subset_omega, weight

      CALL timeset(routineN, handle)

      current = 0.0_dp
      n = SIZE(energies)
      CPASSERT(.NOT. (order < 2 .OR. order > 5 .OR. n < 2))
      CPASSERT(SIZE(omega) == order)
      nh = nlo_derivative_count(order)
      na = nlo_derivative_count(order - 1)
      CPASSERT(ALL(SHAPE(hamiltonian) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(connection) == [n, n, 3, na]))
      CPASSERT(ALL(SHAPE(rho0) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(curvature) == [n, n, 3]))
      CPASSERT(ALL(SHAPE(velocity) == [n, n, 3]))
      CPASSERT(ALL(SHAPE(current) == [3, 3**order]))
      selected = .TRUE.
      IF (PRESENT(elements)) THEN
         CPASSERT(ALL(SHAPE(elements) == SHAPE(current)))
         IF (.NOT. ANY(elements)) THEN
            CALL timestop(handle)
            RETURN
         END IF
         selected(:, 1:3**order) = elements
      END IF
      ! Restrict every labeled subset to the input directions required by selected elements.
      required = .FALSE.
      required(0, 0) = .TRUE.
      DO i_code = 0, 3**order - 1
         IF (.NOT. ANY(selected(:, i_code + 1))) CYCLE
         DO i_mask = 1, 2**order - 1
            i_subset_code = 0
            i_radix = 1
            DO i_field = 1, order
               IF (.NOT. BTEST(i_mask, i_field - 1)) CYCLE
               i_subset_code = i_subset_code + MOD(i_code/3**(i_field - 1), 3)*i_radix
               i_radix = 3*i_radix
            END DO
            required(i_subset_code, i_mask) = .TRUE.
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(omega)) .OR. .NOT. ALL(ieee_is_finite(energies)) .OR. &
          .NOT. ieee_is_finite(rate1) .OR. .NOT. ieee_is_finite(rate2) .OR. &
          .NOT. ieee_is_finite(transition_width) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(velocity, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(velocity))) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(hamiltonian, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(hamiltonian))) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(connection, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(connection))) .OR. &
          .NOT. ALL(ieee_is_finite(REAL(rho0, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(rho0))) .OR. rate1 <= 0.0_dp .OR. &
          rate2 < rate1 .OR. transition_width <= 0.0_dp &
          .OR. ABS(SUM(omega)) <= 1.e-14_dp) THEN
         CPABORT('NLO current response contains invalid data.')
      END IF
      full_mask = 2**order - 1
      DO mask = 0, full_mask
         ncoeff(mask) = nlo_derivative_count(order - POPCNT(mask))
         ncodes(mask) = 3**POPCNT(mask)
      END DO
      ALLOCATE (previous_density(n, n, ncoeff(0)), &
                next_density(n, n, nlo_derivative_count(order - 1)))
      previous_density(:, :, :) = rho0(:, :, 1:ncoeff(0))
      offset_previous = 0
      offset_previous(0, 0) = 1
      extra = 0.0_dp
      DO order_subset = 1, order
         ! Omega_S = sum_(j in S) omega_j; retain each required field subset S.
         CALL subset_density_offsets(order, order_subset, omega, required, ncoeff, &
                                     offset_current, evaluate, representative_code, &
                                     total_coefficients)
         ALLOCATE (density(n, n, total_coefficients))
         density = 0.0_dp
         DO mask = 1, full_mask
            IF (POPCNT(mask) /= order_subset) CYCLE
            subset_omega = 0.0_dp
            DO field = 1, order
               IF (BTEST(mask, field - 1)) subset_omega = subset_omega + omega(field)
            END DO
            DO code = 0, ncodes(mask) - 1
               IF (.NOT. evaluate(code, mask)) CYCLE
               current_start = offset_current(code, mask)
               DO field = 1, order
                  IF (.NOT. BTEST(mask, field - 1)) CYCLE
                  position = POPCNT(IAND(mask, 2**(field - 1) - 1))
                  radix = 3**position
                  direction = MOD(code/radix, 3) + 1
                  multiplicity = 1
                  DO iother = 1, order
                     IF (iother == field .OR. .NOT. BTEST(mask, iother - 1)) CYCLE
                     IF (ABS(omega(iother) - omega(field)) > 0.0_dp) CYCLE
                     iother_position = POPCNT(IAND(mask, 2**(iother - 1) - 1))
                     iother_direction = MOD(code/3**iother_position, 3) + 1
                     IF (iother_direction /= direction) CYCLE
                     IF (iother < field) EXIT
                     multiplicity = multiplicity + 1
                  END DO
                  IF (iother <= order) CYCLE
                  previous_code = MOD(code, radix) + (code/(3*radix))*radix
                  previous_mask = IBCLR(mask, field - 1)
                  previous_start = offset_previous(previous_code, previous_mask)
                  IF (order_subset == order) THEN
                     DO ioutput = 1, 3
                        weight = REAL(multiplicity, KIND=dp)/REAL(order, KIND=dp)
                        curvature_term = curvature_density_trace(curvature, &
                           & previous_density(:, :, previous_start), ioutput, direction)
                        extra(ioutput, code + 1) = extra(ioutput, code + 1) - &
                                                   curvature_term*weight
                     END DO
                  END IF
                  previous_end = previous_start + ncoeff(previous_mask) - 1
                  ! (ε_m-ε_n-Ω_S-iγ_mn)ρ_S = source from lower-order responses.
                  CALL nlo_density_step_derivatives(order - order_subset, direction, &
                                                    subset_omega, energies, rate1, rate2, &
                                                    transition_width, &
                                                    hamiltonian, connection, &
                                                    previous_density(:, :, &
                                                                     previous_start:previous_end), &
                                                    next_density(:, :, 1:ncoeff(mask)))
                  density(:, :, current_start:current_start + ncoeff(mask) - 1) = &
                     density(:, :, current_start:current_start + ncoeff(mask) - 1) + &
                     next_density(:, :, 1:ncoeff(mask))*REAL(multiplicity, &
                                                             KIND=dp)/REAL(order_subset, KIND=dp)
               END DO
            END DO
         END DO
         CALL MOVE_ALLOC(density, previous_density)
         offset_previous = offset_current
      END DO
      DO index = 1, 3**order
         IF (.NOT. ANY(selected(:, index))) CYCLE
         current_start = offset_previous(index - 1, full_mask)
         DO field = 1, 3
            IF (.NOT. selected(field, index)) CYCLE
            ! J_i^(n) = -Tr[v_i rho^(n)] + curvature contribution.
            current(field, index) = -SUM(velocity(:, :, field)* &
                                         TRANSPOSE(previous_density(:, :, current_start))) + &
                                    extra(field, representative_code(index - 1, full_mask) + 1)
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(REAL(current, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(current)))) THEN
         CPABORT('NLO current response contains invalid data.')
      END IF

      CALL timestop(handle)

   END SUBROUTINE nlo_current_tensor_derivatives

! **************************************************************************************************
!> \brief Store each required density once per multiset of (frequency, direction) pairs.
!> \param order Number of incident optical fields.
!> \param degree Number of fields in each subset represented at this stage.
!> \param omega Photon energies of the labeled incident fields.
!> \param required Mask of subset densities needed by the requested tensor elements.
!> \param ncoeff Number of retained k-derivative coefficients for each field subset.
!> \param offsets Starting positions of the compact density arrays.
!> \param evaluate Marks one representative density per equal-frequency field multiset.
!> \param representative_code Tensor-direction code of the representative density.
!> \param total Total retained matrix-derivative coefficients.
!> \note The normalized subset recursion is invariant under relabeling equal-frequency
!>       fields. Integer multiplicities identify exactly equal multisets, with no frequency
!>       tolerance. Offsets address compact matrix-derivative slices of length ncoeff(mask).
! **************************************************************************************************
   SUBROUTINE subset_density_offsets(order, degree, omega, required, ncoeff, offsets, evaluate, &
                                     representative_code, total)
      INTEGER, INTENT(IN)                                :: order, degree
      REAL(KIND=dp), INTENT(IN)                          :: omega(:)
      LOGICAL, INTENT(IN)                                :: required(0:242, 0:31)
      INTEGER, INTENT(IN)                                :: ncoeff(0:31)
      INTEGER, INTENT(OUT)                               :: offsets(0:242, 0:31)
      LOGICAL, INTENT(OUT)                               :: evaluate(0:242, 0:31)
      INTEGER, INTENT(OUT) :: representative_code(0:242, 0:31), total

      INTEGER :: codes(1024), counts(3, 5), groups(5), icode, idir, ifield, imask, iposition, &
         istate, jfield, nstates, signatures(3, 5, 1024), starts(1024)

      groups = 0
      DO ifield = 1, order
         groups(ifield) = ifield
         DO jfield = 1, ifield - 1
            IF (ABS(omega(ifield) - omega(jfield)) > 0.0_dp) CYCLE
            groups(ifield) = groups(jfield)
            EXIT
         END DO
      END DO
      offsets = 0
      evaluate = .FALSE.
      representative_code = 0
      total = 0
      nstates = 0
      DO imask = 1, 2**order - 1
         IF (POPCNT(imask) /= degree) CYCLE
         DO icode = 0, 3**degree - 1
            IF (.NOT. required(icode, imask)) CYCLE
            counts = 0
            iposition = 0
            DO ifield = 1, order
               IF (.NOT. BTEST(imask, ifield - 1)) CYCLE
               idir = MOD(icode/3**iposition, 3) + 1
               counts(idir, groups(ifield)) = counts(idir, groups(ifield)) + 1
               iposition = iposition + 1
            END DO
            DO istate = 1, nstates
               IF (ALL(signatures(:, :, istate) == counts)) EXIT
            END DO
            IF (istate > nstates) THEN
               nstates = istate
               signatures(:, :, istate) = counts
               starts(istate) = total + 1
               codes(istate) = icode
               total = total + ncoeff(imask)
               evaluate(icode, imask) = .TRUE.
            END IF
            offsets(icode, imask) = starts(istate)
            representative_code(icode, imask) = codes(istate)
         END DO
      END DO

   END SUBROUTINE subset_density_offsets

! **************************************************************************************************
!> \brief Trace a precomputed curvature matrix against one lower-order density.
!> \param curvature Three independent Cartesian position-curvature matrices.
!> \param density Lower-order density matrix for the selected field subset.
!> \param ic Cartesian direction of the output current.
!> \param ia Cartesian direction of the incident field.
!> \return Tr[F_(ic,ia) density], with antisymmetry under ic and ia.
! **************************************************************************************************
   FUNCTION curvature_density_trace(curvature, density, ic, ia) RESULT(value)
      COMPLEX(dp), INTENT(IN)                            :: curvature(:, :, :), density(:, :)
      INTEGER, INTENT(IN)                                :: ic, ia
      COMPLEX(dp)                                        :: value

      INTEGER                                            :: lower, pair, upper

      value = 0.0_dp
      IF (ic == ia) RETURN
      lower = MIN(ic, ia)
      upper = MAX(ic, ia)
      pair = (lower - 1)*(6 - lower)/2 + upper - lower
      value = SUM(curvature(:, :, pair)*TRANSPOSE(density))
      IF (ic > ia) value = -value

   END FUNCTION curvature_density_trace

! **************************************************************************************************
!> \brief Advance one optical field in the normalized Cartesian k-derivative recursion.
!> \param max_grade Required highest k-derivative grade, 0:4
!> \param direction Cartesian component of the new input field, 1:3
!> \param omega Cumulative input photon energy in Hartree
!> \param energies Central band energies in Hartree
!> \param rate1 Population damping energy in Hartree
!> \param rate2 Coherence damping energy in Hartree
!> \param transition_width Band-energy width of their smooth transition in Hartree
!> \param hamiltonian Orthogonal-frame Hamiltonian k derivatives through the active response order
!> \param connection Connection k derivatives through one grade below the response order
!> \param previous Previous-order density normalized k derivatives, at least through grade
!> max_grade+1
!> \param response This order's density normalized k derivatives, at least through grade max_grade
!> \note Write X_α=(∂_k^α X)/α! and L_H(X)=[H,X]. At order p, after adding
!>       the field component j with cumulative photon energy Ω, solve
!>       (L_H0-Ω-iΓ_0)ρ_{p,α} = -i(α_j+1)ρ_{p-1,α+e_j}
!>       -Σ_{β+γ=α}[A_{j,β},ρ_{p-1,γ}]
!>       -Σ_{0<β≤α}[H_β,ρ_{p,α-β}]
!>       +iΣ_{0<β≤α}Γ_β(ρ_{p,α-β}).
!>       For equal times Γ_0=1/T1_T2 and all higher Γ derivatives vanish.
!>       For separate times, Γ_H(X)=γ₂X+(γ₁−γ₂)exp[−(ad_H/w)²]X.
!>       The known positive-grade Γ_H terms enter the right-hand side.
!>       The central H_0 is diagonal in this frame.
!> \note If previous has grade <=g+1, this produces the next response through grade g.
!>       Denominators use exp(-i ω t): ε_m-ε_n-ω-iγ_mn.
!>       The field couples as +E.r for an electron, r=i ∂_k+A.
! **************************************************************************************************
   SUBROUTINE nlo_density_step_derivatives(max_grade, direction, omega, energies, rate1, &
                                           rate2, transition_width, &
                                           hamiltonian, connection, previous, response)
      INTEGER, INTENT(IN)                                :: max_grade, direction
      REAL(KIND=dp), INTENT(IN)                          :: omega, energies(:), rate1, rate2, &
                                                            transition_width
      COMPLEX(dp), INTENT(IN)                            :: hamiltonian(:, :, :), &
                                                            connection(:, :, :, :), &
                                                            previous(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: response(:, :, :)

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_density_step_derivatives'

      COMPLEX(dp), ALLOCATABLE                           :: matrix_term(:, :), rhs(:, :), &
                                                            transition_term(:, :, :)
      INTEGER                                            :: grade, handle, i, ia, ib, j, k, n, &
                                                            ngrade, nh, nlower, np, plus(3), &
                                                            remain(3)
      INTEGER, DIMENSION(3, nlo_n_derivatives)           :: alpha
      LOGICAL                                            :: separate_times
      REAL(KIND=dp)                                      :: pair_rate
      TYPE(nlo_band_transition_type)                     :: transition

      CALL timeset(routineN, handle)

      response = 0.0_dp
      n = SIZE(energies)
      CPASSERT(max_grade >= 0 .AND. max_grade <= 4)
      CPASSERT(direction >= 1 .AND. direction <= 3)
      CPASSERT(n >= 2)
      nh = nlo_derivative_count(max_grade)
      np = nlo_derivative_count(max_grade + 1)
      CPASSERT(SIZE(hamiltonian, 1) == n)
      CPASSERT(SIZE(hamiltonian, 2) == n)
      CPASSERT(SIZE(hamiltonian, 3) >= nh)
      CPASSERT(SIZE(connection, 1) == n)
      CPASSERT(SIZE(connection, 2) == n)
      CPASSERT(SIZE(connection, 3) == 3)
      CPASSERT(SIZE(connection, 4) >= nh)
      CPASSERT(SIZE(previous, 1) == n)
      CPASSERT(SIZE(previous, 2) == n)
      CPASSERT(SIZE(previous, 3) >= np)
      CPASSERT(SIZE(response, 1) == n)
      CPASSERT(SIZE(response, 2) == n)
      CPASSERT(SIZE(response, 3) >= nh)
      IF (.NOT. ieee_is_finite(omega) .OR. .NOT. ALL(ieee_is_finite(energies)) .OR. &
          .NOT. ieee_is_finite(rate1) .OR. .NOT. ieee_is_finite(rate2) .OR. &
          .NOT. ieee_is_finite(transition_width) .OR. rate1 <= 0.0_dp .OR. &
          rate2 < rate1 .OR. transition_width <= 0.0_dp .OR. &
          .NOT. ALL(ieee_is_finite(REAL(previous, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(previous)))) THEN
         CPABORT('NLO density response contains invalid data.')
      END IF
      ALLOCATE (matrix_term(n, n), rhs(n, n))
      separate_times = rate1 < rate2
      IF (separate_times) THEN
         ALLOCATE (transition_term(n, n, (max_grade + 1)*(max_grade + 2)/2))
         CALL nlo_band_transition_prepare(energies, transition_width, rate1, rate2, transition)
      END IF
      CALL nlo_derivative_indices(alpha)
      nlower = 0
      DO ia = 1, nlo_derivative_count(max_grade)
         grade = SUM(alpha(:, ia))
         IF (separate_times .AND. grade > 0) THEN
            IF (ia == nlo_derivative_count(grade - 1) + 1) THEN
               nlower = ia - 1
               ngrade = nlo_derivative_count(grade) - nlower
               CALL nlo_band_transition_correction(transition, hamiltonian, response, grade, &
                                                   transition_term(:, :, :ngrade))
            END IF
         END IF
         plus = alpha(:, ia)
         plus(direction) = plus(direction) + 1
         ib = derivatives_index(alpha, plus)
         CPASSERT(ib /= 0)
         ! D_a rho = d_a rho - i[A_a,rho]; lower H derivatives enter recursively.
         rhs(:, :) = -gaussi*REAL(plus(direction), KIND=dp)*previous(:, :, ib)
         DO ib = 1, ia
            remain = alpha(:, ia) - alpha(:, ib)
            j = derivatives_index(alpha, remain)
            IF (j == 0) CYCLE
            matrix_term(:, :) = MATMUL(connection(:, :, direction, ib), previous(:, :, j))
            rhs(:, :) = rhs - matrix_term
            matrix_term(:, :) = MATMUL(previous(:, :, j), connection(:, :, direction, ib))
            rhs(:, :) = rhs + matrix_term
         END DO
         DO ib = 2, ia
            remain = alpha(:, ia) - alpha(:, ib)
            j = derivatives_index(alpha, remain)
            IF (j == 0) CYCLE
            matrix_term(:, :) = MATMUL(hamiltonian(:, :, ib), response(:, :, j))
            rhs(:, :) = rhs - matrix_term
            matrix_term(:, :) = MATMUL(response(:, :, j), hamiltonian(:, :, ib))
            rhs(:, :) = rhs + matrix_term
         END DO
         IF (separate_times .AND. grade > 0) THEN
            rhs(:, :) = rhs + gaussi*transition_term(:, :, ia - nlower)
         END IF
         DO k = 1, n
            DO i = 1, n
               pair_rate = rate1
               IF (separate_times) pair_rate = nlo_band_transition_rate(transition, i, k)
               response(i, k, ia) = rhs(i, k)/ &
                                    CMPLX(energies(i) - energies(k) - omega, -pair_rate, dp)
            END DO
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(REAL(response, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(response)))) THEN
         CPABORT('NLO density response contains invalid data.')
      END IF

      CALL timestop(handle)

   END SUBROUTINE nlo_density_step_derivatives

END MODULE nlo_current
