!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \author Jan Wilhelm, with help from ChatGPT

MODULE nlo_integration
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE kubo_projected,                  ONLY: projected_eigensystem
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi
   USE nlo_current,                     ONLY: nlo_current_tensor_derivatives
   USE nlo_derivatives,                 ONLY: nlo_derivative_count,&
                                              nlo_derivative_indices,&
                                              nlo_equilibrium_density_derivatives,&
                                              nlo_n_derivatives,&
                                              nlo_orthonormal_operator_derivatives
   USE nlo_read_input,                  ONLY: nlo_frequency_case
   USE nlo_types,                       ONLY: nlo_type
   USE physcon,                         ONLY: evolt
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo_integration'

   PUBLIC :: nlo_bz_integration

CONTAINS

! **************************************************************************************************
!> \brief Integrate the density-response current over the Brillouin zone.
!> \param nlo_env NLO settings and data for the calculation
! **************************************************************************************************
   SUBROUTINE nlo_bz_integration(nlo_env)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_bz_integration'

      COMPLEX(dp), ALLOCATABLE                           :: curvature(:, :, :), matrix_term(:, :)
      INTEGER                                            :: handle, i_case, ia, ic, ik, iw, nocc, &
                                                            pair, spin
      LOGICAL                                            :: zero_output
      REAL(KIND=dp)                                      :: energy_ev(5), k(3), local_projection, &
                                                            occ_factor, omega(5)

      CALL timeset(routineN, handle)

      CALL allocate_nlo_quantities(nlo_env)
      ALLOCATE (curvature(nlo_env%n_bands, nlo_env%n_bands, 3), &
                matrix_term(nlo_env%n_bands, nlo_env%n_bands))

      nlo_env%projection_error = 0.0_dp
      nlo_env%min_direct_gap = HUGE(1.0_dp)
      nlo_env%progress_total = 1 + (nlo_env%n_kpoints - 1)/nlo_env%para_env%num_pe
      IF (nlo_env%progress_total > HUGE(nlo_env%progress_total)/nlo_env%n_spin) THEN
         CPABORT('NLO work count exceeds integer indexing.')
      END IF
      nlo_env%progress_total = nlo_env%n_spin*nlo_env%progress_total
      nlo_env%progress_done = 0
      nlo_env%progress_next = 1
      nlo_env%progress_start = m_walltime()

      DO spin = 1, nlo_env%n_spin
         IF (nlo_env%mos(spin)%maxocc <= 0.0_dp) THEN
            CPABORT('NONLINEAR_OPTICS requires positive band occupation.')
         END IF
         occ_factor = nlo_env%mos(spin)%maxocc
         nocc = NINT(REAL(nlo_env%mos(spin)%nelectron, KIND=dp)/occ_factor)
         IF (nocc < 1 .OR. nocc >= nlo_env%n_bands .OR. &
             ABS(nocc*occ_factor - nlo_env%mos(spin)%nelectron) > 1.e-8_dp) THEN
            CPABORT('NLO needs integer occupied-band count and empty bands.')
         END IF
         nlo_env%occupied = 0
         nlo_env%occupied(1:nocc) = 1
         DO ik = 1, nlo_env%n_kpoints
            IF (MOD(ik - 1, nlo_env%para_env%num_pe) /= nlo_env%para_env%mepos) CYCLE
            k = REAL([MOD(ik - 1, nlo_env%k_mesh(1)), &
                      MOD((ik - 1)/nlo_env%k_mesh(1), nlo_env%k_mesh(2)), &
                      (ik - 1)/(nlo_env%k_mesh(1)*nlo_env%k_mesh(2))], KIND=dp)/ &
                REAL(nlo_env%k_mesh, KIND=dp)
            ! Build H(k), S(k), d(k) and their analytic Cartesian k derivatives.
            CALL ao_derivatives_at_k(nlo_env, k, spin)
            CALL hermitize_derivatives(nlo_env%hamiltonian_derivatives, local_projection)

            nlo_env%projection_error = MAX(nlo_env%projection_error, local_projection)
            IF (local_projection > 1.e-5_dp) THEN
               CPABORT('NLO Hamiltonian images are not Hermitian.')
            END IF
            ! X = S^(-1/2) orthonormalizes the AO basis; U diagonalizes X H X here.
            ! Express all operator derivatives in this same eigenvector basis U.
            CALL nlo_ao_to_band_derivatives( &
               nlo_env%overlap_derivatives, nlo_env%hamiltonian_derivatives, &
               nlo_env%position_derivatives, &
               nlo_env%occupied, nlo_env%energies, nlo_env%band_hamiltonian_derivatives, &
               nlo_env%band_connection_derivatives, nlo_env%rho0_derivatives, &
               nlo_env%velocity, &
               active_grade=nlo_env%order)

            nlo_env%min_direct_gap = MIN( &
                                     nlo_env%min_direct_gap, &
                                     nlo_env%energies(nocc + 1) - nlo_env%energies(nocc))
            ! F_ca = d_c A_a - d_a A_c - i[A_c,A_a], reused for every frequency.
            DO ic = 1, 2
               DO ia = ic + 1, 3
                  pair = (ic - 1)*(6 - ic)/2 + ia - ic
                  curvature(:, :, pair) = &
                     nlo_env%band_connection_derivatives(:, :, ia, ic + 1) - &
                     nlo_env%band_connection_derivatives(:, :, ic, ia + 1)
                  matrix_term(:, :) = MATMUL(nlo_env%band_connection_derivatives(:, :, ic, 1), &
                                             nlo_env%band_connection_derivatives(:, :, ia, 1))
                  curvature(:, :, pair) = curvature(:, :, pair) - gaussi*matrix_term
                  matrix_term(:, :) = MATMUL(nlo_env%band_connection_derivatives(:, :, ia, 1), &
                                             nlo_env%band_connection_derivatives(:, :, ic, 1))
                  curvature(:, :, pair) = curvature(:, :, pair) + gaussi*matrix_term
               END DO
            END DO
            DO i_case = 1, nlo_env%n_frequency_combinations
               CALL nlo_frequency_case(nlo_env%frequency, nlo_env%order, i_case, &
                                       energy_ev, zero_output)
               IF (zero_output) THEN
                  CPABORT('Invalid NONLINEAR_OPTICS frequency case.')
               END IF
               omega = 0.0_dp
               omega(1:nlo_env%order) = energy_ev(1:nlo_env%order)/evolt
               ! Compute J_i^(n)(k) from rho^(n): -Tr[v_i rho^(n)] plus curvature current.
               ! Each rho^(n) follows from lower orders at Omega = sum_j omega_j.
               CALL nlo_current_tensor_derivatives( &
                  nlo_env%order, omega(1:nlo_env%order), nlo_env%energies, &
                  1.0_dp/nlo_env%t1_au, 1.0_dp/nlo_env%t2_au, &
                  nlo_env%transition_width_au, nlo_env%band_hamiltonian_derivatives, &
                  nlo_env%band_connection_derivatives, curvature, &
                  nlo_env%rho0_derivatives, &
                  nlo_env%velocity, &
                  nlo_env%current_flat, elements=nlo_env%elements)

               ! Sum J_i^(n)(k) over spin and k; f_spin is the band occupation factor.
               ! I_i^(n) = sum_(spin,k) f_spin J_i^(n)(k) / N_k.
               nlo_env%current_spectrum(:, :, i_case) = &
                  nlo_env%current_spectrum(:, :, i_case) + &
                  occ_factor*nlo_env%current_flat/REAL(nlo_env%n_kpoints, KIND=dp)
            END DO
            IF (nlo_env%para_env%mepos == 0) CALL report_nlo_progress(nlo_env)

         END DO
      END DO
      CALL nlo_env%para_env%max(nlo_env%projection_error)
      CALL nlo_env%para_env%min(nlo_env%min_direct_gap)
      iw = cp_logger_get_default_io_unit()
      IF (iw > 0) THEN
         WRITE (iw, '(T2,A,T69,F9.2,1X,A)') &
            'NONLINEAR_OPTICS| Minimum occupied-empty gap:', &
            nlo_env%min_direct_gap*evolt, 'eV'
      END IF
      CALL release_nlo_quantities(nlo_env)

      CALL timestop(handle)

   END SUBROUTINE nlo_bz_integration

! **************************************************************************************************
!> \brief Project all analytic Hamiltonian derivatives onto the Hermitian image subspace.
!> \param derivatives AO Hamiltonian derivatives over lattice images
!> \param relative_error Largest relative change from enforcing Hermiticity
!> \note This linear projection commutes with all k derivatives.
! **************************************************************************************************
   SUBROUTINE hermitize_derivatives(derivatives, relative_error)
      COMPLEX(dp), INTENT(INOUT)                         :: derivatives(:, :, :)
      REAL(KIND=dp), INTENT(OUT)                         :: relative_error

      COMPLEX(dp) :: adjoint(SIZE(derivatives, 1), SIZE(derivatives, 2))
      INTEGER                                            :: ij
      REAL(KIND=dp)                                      :: scale

      relative_error = 0.0_dp
      DO ij = 1, SIZE(derivatives, 3)
         adjoint = CONJG(TRANSPOSE(derivatives(:, :, ij)))
         scale = MAX(1.0_dp, MAXVAL(ABS(derivatives(:, :, ij))))
         relative_error = MAX(relative_error, MAXVAL(ABS(derivatives(:, :, ij) - adjoint))/scale)
         derivatives(:, :, ij) = 0.5_dp*(derivatives(:, :, ij) + adjoint)
      END DO

   END SUBROUTINE hermitize_derivatives

! **************************************************************************************************
!> \brief Fourier-transform AO H, S, and position matrices at one optical k point.
!> \param nlo_env NLO settings and data for the calculation
!> \param k Optical k point in fractional coordinates
!> \param spin Spin channel of the converged SCF solution
!> \note For X = H, S, or the AO position matrix d_a, and Cartesian k:
!>       dX(k)/dk_x = sum_R i R_x exp(i k·R) X(R).
!>       R is a real-space lattice vector. The input k uses fractional coordinates:
!>       k·R = 2 pi k_frac·R_int. Higher derivatives add factors i R_a;
!>       each stored Taylor coefficient is divided by its factorial.
! **************************************************************************************************
   SUBROUTINE ao_derivatives_at_k(nlo_env, k, spin)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env
      REAL(KIND=dp), INTENT(IN)                          :: k(3)
      INTEGER, INTENT(IN)                                :: spin

      CHARACTER(len=*), PARAMETER :: routineN = 'ao_derivatives_at_k'

      INTEGER                                            :: alpha(3, nlo_n_derivatives), grade, &
                                                            handle, ia, ij, j, n, na, nh, ni
      REAL(KIND=dp)                                      :: divisor

      CALL timeset(routineN, handle)

      ASSOCIATE (h_image => nlo_env%h_ao_R(spin, :, :, :), &
                 s_image => nlo_env%s_ao_R(1, :, :, :), &
                 moment_image => nlo_env%d_ao_R, &
                 image_cells => nlo_env%images%index_to_cell, &
                 hmat => nlo_env%cell%hmat, &
                 h_derivatives => nlo_env%hamiltonian_derivatives, &
                 s_derivatives => nlo_env%overlap_derivatives, &
                 moment_derivatives => nlo_env%position_derivatives)
         h_derivatives = 0.0_dp
         s_derivatives = 0.0_dp
         moment_derivatives = 0.0_dp
         n = SIZE(h_image, 1)
         ni = SIZE(s_image, 3)
         CPASSERT(n >= 2)
         CPASSERT(SIZE(h_image, 2) == n)
         CPASSERT(SIZE(h_image, 3) >= 1)
         CPASSERT(SIZE(s_image, 1) == n)
         CPASSERT(SIZE(s_image, 2) == n)
         CPASSERT(ni >= 1)
         CPASSERT(ALL(SHAPE(moment_image) == [3, n, n, ni]))
         CPASSERT(SIZE(h_image, 3) == ni)
         CPASSERT(ALL(SHAPE(image_cells) == [3, ni]))
         nh = nlo_derivative_count(nlo_env%order)
         na = nlo_derivative_count(nlo_env%order - 1)
         CPASSERT(ALL(SHAPE(h_derivatives) == [n, n, nh]))
         CPASSERT(ALL(SHAPE(s_derivatives) == [n, n, nh]))
         CPASSERT(ALL(SHAPE(moment_derivatives) == [n, n, 3, na]))
         IF (.NOT. ALL(ieee_is_finite(h_image)) .OR. .NOT. ALL(ieee_is_finite(s_image)) .OR. &
             .NOT. ALL(ieee_is_finite(moment_image)) .OR. &
             .NOT. ALL(ieee_is_finite(k)) .OR. .NOT. ALL(ieee_is_finite(hmat))) THEN
            CPABORT('NONLINEAR_OPTICS numerical derivative construction failed.')
         END IF
         grade = nlo_env%order
         CPASSERT(grade >= 2 .AND. grade <= 5)
         CALL nlo_derivative_indices(alpha)
         DO ij = 1, nlo_derivative_count(grade)
            divisor = 1.0_dp
            DO ia = 1, 3
               DO j = 2, alpha(ia, ij)
                  divisor = divisor*REAL(j, KIND=dp)
               END DO
            END DO
            CALL rs_to_kp_multi_deriv(h_image, h_derivatives(:, :, ij), image_cells, k, &
                                      alpha(:, ij), hmat)
            CALL rs_to_kp_multi_deriv(s_image, s_derivatives(:, :, ij), image_cells, k, &
                                      alpha(:, ij), hmat)
            h_derivatives(:, :, ij) = h_derivatives(:, :, ij)/divisor
            s_derivatives(:, :, ij) = s_derivatives(:, :, ij)/divisor
            IF (SUM(alpha(:, ij)) >= grade) CYCLE
            DO ia = 1, 3
               CALL rs_to_kp_multi_deriv(moment_image(ia, :, :, :), &
                                         moment_derivatives(:, :, ia, ij), image_cells, k, &
                                         alpha(:, ij), hmat)
               moment_derivatives(:, :, ia, ij) = moment_derivatives(:, :, ia, ij)/divisor
            END DO
         END DO
         IF (.NOT. ALL(ieee_is_finite(REAL(h_derivatives, KIND=dp))) &
             .OR. .NOT. ALL(ieee_is_finite(AIMAG(h_derivatives))) .OR. &
             .NOT. ALL(ieee_is_finite(REAL(s_derivatives, KIND=dp))) &
             .OR. .NOT. ALL(ieee_is_finite(AIMAG(s_derivatives))) .OR. &
             .NOT. ALL(ieee_is_finite(REAL(moment_derivatives, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(moment_derivatives)))) THEN
            CPABORT('NONLINEAR_OPTICS numerical derivative construction failed.')
         END IF

      END ASSOCIATE

      CALL timestop(handle)

   END SUBROUTINE ao_derivatives_at_k

! **************************************************************************************************
!> \brief Convert AO operators to band matrices and build the equilibrium density.
!> \param overlap AO overlap and its k derivatives
!> \param hamiltonian AO Hamiltonian and its k derivatives
!> \param moment AO position-matrix derivatives through one grade below the response order
!> \param occupied One for an occupied band, zero for an empty band
!> \param energies Band energies in Hartree at this k point
!> \param band_h Hamiltonian derivatives in the eigenvector basis
!> \param band_a Position connection derivatives in the eigenvector basis
!> \param rho0 Equilibrium density and its k derivatives for unit occupation
!> \param velocity Band velocity at this k point, in bohr per atomic time
!> \param active_grade Highest derivative grade required for the response
!> \note First X=S^(-1/2) makes the AO basis orthonormal. U diagonalizes
!>       X H X at this k point. The same U is applied to every k derivative
!>       of H and the position matrix, so no derivative of U is needed.
!>       Occupied bands then define rho0; field-driven densities are built later.
! **************************************************************************************************
   SUBROUTINE nlo_ao_to_band_derivatives(overlap, hamiltonian, moment, occupied, &
                                         energies, band_h, band_a, rho0, velocity, &
                                         active_grade)
      COMPLEX(dp), DIMENSION(:, :, :), INTENT(IN)        :: overlap, hamiltonian
      COMPLEX(dp), INTENT(IN)                            :: moment(:, :, :, :)
      INTEGER, INTENT(IN)                                :: occupied(:)
      REAL(KIND=dp), INTENT(OUT)                         :: energies(:)
      COMPLEX(dp), INTENT(OUT)                           :: band_h(:, :, :), band_a(:, :, :, :), &
                                                            rho0(:, :, :), velocity(:, :, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: active_grade

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_ao_to_band_derivatives'

      COMPLEX(dp), ALLOCATABLE                           :: identity(:, :), matrix_term(:, :), &
                                                            orth_h(:, :, :), u(:, :)
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :, :, :)    :: connection
      INTEGER                                            :: alpha(3, nlo_n_derivatives), grade, &
                                                            handle, ia, ideriv, ij, info, n, na, nh

      CALL timeset(routineN, handle)

      energies = 0.0_dp
      band_h = 0.0_dp
      band_a = 0.0_dp
      rho0 = 0.0_dp
      velocity = 0.0_dp
      n = SIZE(energies)
      CPASSERT(.NOT. (n < 2 .OR. SIZE(occupied) /= n))
      grade = 5
      IF (PRESENT(active_grade)) grade = active_grade
      CPASSERT(grade >= 2 .AND. grade <= 5)
      nh = nlo_derivative_count(grade)
      na = nlo_derivative_count(grade - 1)
      CPASSERT(ALL(SHAPE(overlap) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(hamiltonian) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(moment) == [n, n, 3, na]))
      CPASSERT(ALL(SHAPE(band_h) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(band_a) == [n, n, 3, na]))
      CPASSERT(ALL(SHAPE(rho0) == [n, n, nh]))
      CPASSERT(ALL(SHAPE(velocity) == [n, n, 3]))
      ALLOCATE (connection(n, n, 3, na), identity(n, n), matrix_term(n, n), &
                orth_h(n, n, nh), u(n, n))

      CALL nlo_orthonormal_operator_derivatives(overlap, hamiltonian, moment, orth_h, &
                                                connection, active_grade=grade)

      identity = 0.0_dp
      DO ij = 1, n
         identity(ij, ij) = 1.0_dp
      END DO

      CALL projected_eigensystem(orth_h(:, :, 1), identity, energies, u, info)
      IF (info /= 0) CPABORT('NONLINEAR_OPTICS band eigensolver failed.')

      DO ij = 1, nh
         matrix_term(:, :) = MATMUL(orth_h(:, :, ij), u)
         band_h(:, :, ij) = MATMUL(CONJG(TRANSPOSE(u)), matrix_term)
         IF (ij > na) CYCLE
         DO ia = 1, 3
            matrix_term(:, :) = MATMUL(connection(:, :, ia, ij), u)
            band_a(:, :, ia, ij) = MATMUL(CONJG(TRANSPOSE(u)), matrix_term)
         END DO
      END DO

      CALL nlo_equilibrium_density_derivatives(energies, occupied, band_h, rho0, &
                                               active_grade=grade)

      CALL nlo_derivative_indices(alpha)
      DO ia = 1, 3
         ideriv = 0
         DO ij = 1, nh
            IF (alpha(ia, ij) == 1 .AND. SUM(alpha(:, ij)) == 1) ideriv = ij
         END DO
         CPASSERT(ideriv /= 0)
         matrix_term(:, :) = MATMUL(band_a(:, :, ia, 1), band_h(:, :, 1))
         velocity(:, :, ia) = band_h(:, :, ideriv) - gaussi*matrix_term
         matrix_term(:, :) = MATMUL(band_h(:, :, 1), band_a(:, :, ia, 1))
         velocity(:, :, ia) = velocity(:, :, ia) + gaussi*matrix_term
      END DO
      IF (.NOT. ALL(ieee_is_finite(REAL(velocity, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(velocity)))) THEN
         CPABORT('NONLINEAR_OPTICS numerical derivative construction failed.')
      END IF

      CALL timestop(handle)

   END SUBROUTINE nlo_ao_to_band_derivatives

! **************************************************************************************************
!> \brief Analytic Cartesian mixed k derivative of a real-space AO matrix.
!> \param rs_real Real-space matrix for every lattice image
!> \param ks_complex Requested derivative in k space
!> \param index_to_cell Integer lattice vector for each image
!> \param xkp k point in reciprocal fractional coordinates
!> \param deriv_orders Number of derivatives in Cartesian x, y, z
!> \param hmat Direct-cell vectors as columns, in the same length units as k^-1
!> \note The Fourier factor is exp(2πi k_frac dot R_int). Each Cartesian derivative
!>       contributes i R_cart. This avoids phase-sensitive finite differences.
! **************************************************************************************************
   SUBROUTINE rs_to_kp_multi_deriv(rs_real, ks_complex, index_to_cell, xkp, deriv_orders, hmat)
      REAL(KIND=dp), INTENT(IN)                          :: rs_real(:, :, :)
      COMPLEX(dp), INTENT(OUT)                           :: ks_complex(:, :)
      INTEGER, INTENT(IN)                                :: index_to_cell(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(3)
      INTEGER, INTENT(IN)                                :: deriv_orders(3)
      REAL(KIND=dp), INTENT(IN)                          :: hmat(3, 3)

      COMPLEX(dp)                                        :: factor
      INTEGER                                            :: i_axis, i_image, n_images
      REAL(KIND=dp)                                      :: cell_index(3), cell_vector(3), phase

      IF (ANY(deriv_orders < 0)) CPABORT('Negative k derivative order requested.')
      n_images = SIZE(rs_real, 3)
      IF (n_images < 1) CPABORT('No real-space lattice images supplied.')
      IF (SIZE(index_to_cell, 1) /= 3 .OR. SIZE(index_to_cell, 2) /= n_images) THEN
         CPABORT('Real-space matrices and lattice-image indices are inconsistent.')
      END IF
      IF (SIZE(ks_complex, 1) /= SIZE(rs_real, 1) .OR. SIZE(ks_complex, 2) /= SIZE(rs_real, 2)) THEN
         CPABORT('Real-space and k-space matrix shapes are inconsistent.')
      END IF
      ks_complex = CMPLX(0.0_dp, 0.0_dp, kind=dp)
      DO i_image = 1, n_images
         cell_index = REAL(index_to_cell(:, i_image), KIND=dp)
         cell_vector = MATMUL(hmat, cell_index)
         phase = twopi*SUM(xkp*cell_index)
         factor = CMPLX(COS(phase), SIN(phase), kind=dp)
         DO i_axis = 1, 3
            factor = factor*CMPLX(0.0_dp, cell_vector(i_axis), kind=dp)**deriv_orders(i_axis)
         END DO
         ks_complex = ks_complex + factor*rs_real(:, :, i_image)
      END DO

   END SUBROUTINE rs_to_kp_multi_deriv

! **************************************************************************************************
!> \brief Allocate the per-k matrices and rank-local frequency accumulators once.
!> \param nlo_env NLO settings and data for the calculation
! **************************************************************************************************
   SUBROUTINE allocate_nlo_quantities(nlo_env)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      INTEGER                                            :: na, nh

      nh = nlo_derivative_count(nlo_env%order)
      na = nlo_derivative_count(nlo_env%order - 1)
      ALLOCATE (nlo_env%energies(nlo_env%n_bands))
      ALLOCATE (nlo_env%occupied(nlo_env%n_bands))
      ALLOCATE (nlo_env%hamiltonian_derivatives(nlo_env%n_bands, nlo_env%n_bands, nh))
      ALLOCATE (nlo_env%overlap_derivatives(nlo_env%n_bands, nlo_env%n_bands, nh))
      ALLOCATE (nlo_env%position_derivatives(nlo_env%n_bands, nlo_env%n_bands, 3, na))
      ALLOCATE (nlo_env%band_hamiltonian_derivatives(nlo_env%n_bands, nlo_env%n_bands, nh))
      ALLOCATE (nlo_env%band_connection_derivatives(nlo_env%n_bands, nlo_env%n_bands, 3, na))
      ALLOCATE (nlo_env%rho0_derivatives(nlo_env%n_bands, nlo_env%n_bands, nh))
      ALLOCATE (nlo_env%velocity(nlo_env%n_bands, nlo_env%n_bands, 3))
      ALLOCATE (nlo_env%current_flat(3, 3**nlo_env%order))
      ALLOCATE (nlo_env%current_spectrum(3, 3**nlo_env%order, nlo_env%n_frequency_combinations), &
                SOURCE=(0.0_dp, 0.0_dp))

   END SUBROUTINE allocate_nlo_quantities

! **************************************************************************************************
!> \brief Print completed spin/k tasks when the next percentage threshold is crossed.
!> \param nlo_env NLO settings and data for the calculation
! **************************************************************************************************
   SUBROUTINE report_nlo_progress(nlo_env)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      REAL(KIND=dp), PARAMETER :: &
         marks(9) = [0.1_dp, 0.3_dp, 1.0_dp, 3.0_dp, 10.0_dp, 30.0_dp, 60.0_dp, 90.0_dp, 100.0_dp]

      INTEGER                                            :: iw
      REAL(KIND=dp)                                      :: percent

      nlo_env%progress_done = nlo_env%progress_done + 1
      percent = 100.0_dp*REAL(nlo_env%progress_done, KIND=dp)/REAL(nlo_env%progress_total, KIND=dp)
      CPASSERT(.NOT. (nlo_env%progress_next > SIZE(marks)))
      IF (percent < marks(nlo_env%progress_next)) THEN
         RETURN
      END IF
      iw = cp_logger_get_default_io_unit()
      IF (iw > 0) THEN
         IF (nlo_env%progress_total < 10000) THEN
            WRITE (iw, '(T2,A,I0,A,I0,7X,F5.1,A,4X,A,T69,F10.2,1X,A)') &
               'NONLINEAR_OPTICS| Progress:  ', nlo_env%progress_done, ' / ', &
               nlo_env%progress_total, percent, '%', 'Total time:', &
               m_walltime() - nlo_env%progress_start, 's'
         ELSE
            WRITE (iw, '(T2,A,I0,A,I0,7X,F5.1,A)') &
               'NONLINEAR_OPTICS| Progress:  ', nlo_env%progress_done, ' / ', &
               nlo_env%progress_total, percent, '%'
            WRITE (iw, '(T2,A,T69,F10.2,1X,A)') &
               'NONLINEAR_OPTICS| Total time:', m_walltime() - nlo_env%progress_start, 's'
         END IF
         CALL m_flush(iw)
      END IF
      DO WHILE (nlo_env%progress_next <= SIZE(marks))
         IF (percent < marks(nlo_env%progress_next)) EXIT
         nlo_env%progress_next = nlo_env%progress_next + 1
      END DO

   END SUBROUTINE report_nlo_progress

! **************************************************************************************************
!> \brief Free the per-k workspace after the MPI reduction and validation.
!> \param nlo_env NLO settings and data for the calculation
! **************************************************************************************************
   SUBROUTINE release_nlo_quantities(nlo_env)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      DEALLOCATE (nlo_env%energies, nlo_env%occupied)
      DEALLOCATE (nlo_env%hamiltonian_derivatives, nlo_env%overlap_derivatives)
      DEALLOCATE (nlo_env%position_derivatives)
      DEALLOCATE (nlo_env%band_hamiltonian_derivatives, nlo_env%band_connection_derivatives)
      DEALLOCATE (nlo_env%rho0_derivatives, nlo_env%velocity, nlo_env%current_flat)

   END SUBROUTINE release_nlo_quantities

END MODULE nlo_integration
