!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse of a specified finite-rank UPF reconstruction, including intercenter coupling.
!>        Delta contains original AE-minus-PS waves, P the original projector bras.
!>        For K=I+P*Delta and B=P*chi, W=K^-1*B gives U*chi=chi-Delta*W and P*U*chi=W.
!>        Atomic hard projectors use the nonoverlapping one-center construction. The general
!>        complex inverse remains available for checking the original functional explicitly.
! **************************************************************************************************
MODULE upf_inverse
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_evaluate,&
                                              upf_radial_fields_init,&
                                              upf_radial_fields_type
   USE upf_types,                       ONLY: upf_radial_weights
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE upf_inverse_type
      INTEGER :: nproj = 0
      LOGICAL :: ready = .FALSE.
      REAL(KIND=dp) :: rcond = 0.0_dp
      COMPLEX(KIND=dp), ALLOCATABLE :: lu(:, :)
      INTEGER, ALLOCATABLE :: pivots(:)
   END TYPE upf_inverse_type

   TYPE upf_inverse_kinetic_type
      TYPE(upf_radial_fields_type) :: delta, pseudo
   END TYPE upf_inverse_kinetic_type

   TYPE upf_inverse_atomic_type
      ! The projector coupling is D-T in the original partial-wave channel order.
      TYPE(upf_projector_type) :: hard, difference
      TYPE(upf_inverse_type) :: inverse
      REAL(dp), ALLOCATABLE :: kinetic(:, :), metric(:, :), wave_overlap(:, :), difference_kinetic(:, :), difference_overlap(:, :)
   END TYPE upf_inverse_atomic_type

   PUBLIC :: upf_inverse_type, upf_inverse_init, upf_inverse_project, upf_inverse_project_response
   PUBLIC :: upf_inverse_basis, upf_inverse_basis_response
   PUBLIC :: upf_inverse_operator, upf_inverse_operator_response
   PUBLIC :: upf_inverse_delta_init, upf_inverse_fields_init
   PUBLIC :: upf_inverse_kinetic_type, upf_inverse_kinetic_init, upf_inverse_kinetic_evaluate
   PUBLIC :: upf_inverse_atomic_type, upf_inverse_atomic_init, upf_inverse_kinetic_matrix
   PUBLIC :: upf_inverse_atomic_residual

CONTAINS

! **************************************************************************************************
!> \brief Local hard projectors and D-T for the nonoverlapping partial-wave construction.
!>        Both factors of a nonlocal term use these transformed projectors. Original Q and
!>        norm data remain independent and are not replaced by a partial-wave identity here.
!>        The hard reference may be an AE or a harder NC wave, as supplied by its generator.
!> \param pot Original scalar NC/USPP and compatible projector-indexed hard/soft reference pairs.
!> \param degree Largest Gaussian Cartesian degree plus derivatives or moments.
!> \param atomic New owned one-center inverse, hard projectors and kinetic correction.
!> \param info Zero on success, otherwise the atomic inverse-factorization status.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE upf_inverse_atomic_init(pot, degree, atomic, info, direct)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: degree
      TYPE(upf_inverse_atomic_type), INTENT(OUT)         :: atomic
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: direct

      COMPLEX(dp), ALLOCATABLE                           :: coupling(:, :), identity(:, :), &
                                                            rotation(:, :)
      INTEGER                                            :: i, j, nb, p, q
      REAL(dp)                                           :: norm_difference
      REAL(dp), ALLOCATABLE                              :: delta_kinetic(:, :), &
                                                            delta_overlap(:, :), &
                                                            radial_kinetic(:, :), &
                                                            radial_weights(:), weights(:)
      TYPE(atom_upfpot_type)                             :: hard, radial
      TYPE(upf_inverse_kinetic_type)                     :: kinetic

      nb = pot%number_of_proj
      CPASSERT(nb > 0)
      info = 0
      IF (pot%is_paw .OR. TRIM(pot%pseudo_type) == "PAW") THEN
         CPABORT("Atom-local inverse PAW fixed-term accounting is not implemented")
      END IF
      CPASSERT(pot%is_ultrasoft)
      IF (PRESENT(direct)) THEN
         IF (direct) THEN
            CALL upf_projector_init(pot, degree, atomic%hard)
            CALL compact_projector(atomic%hard)
            RETURN
         END IF
      END IF
      CALL inverse_radial_data(pot, radial)
      CALL upf_projector_init(radial, degree, atomic%difference)
      CALL compact_projector(atomic%difference)
      CALL upf_inverse_kinetic_init(pot, kinetic)
      CALL upf_inverse_kinetic_matrix(pot, kinetic, radial_kinetic, delta_kinetic, delta_overlap)
      ALLOCATE (radial_weights, SOURCE=upf_radial_weights(pot))
      ALLOCATE (weights, SOURCE=radial_weights*pot%r**2)
      ALLOCATE (coupling(nb, nb), identity(nb, nb), rotation(nb, nb))
      coupling = CMPLX(0.0_dp, 0.0_dp, dp)
      identity = CMPLX(0.0_dp, 0.0_dp, dp)
      DO j = 1, nb
         identity(j, j) = CMPLX(1.0_dp, 0.0_dp, dp)
         DO i = 1, nb
            IF (pot%lbeta(i) /= pot%lbeta(j)) CYCLE
            coupling(i, j) = CMPLX(SUM(weights*pot%beta(:, i)*radial%beta(:, j)), 0.0_dp, dp)
         END DO
      END DO
      CALL upf_inverse_init(coupling, atomic%inverse, info)
      IF (info /= 0) RETURN
      CALL upf_inverse_project(atomic%inverse, identity, rotation)
      ! Bra_i = sum_j K^-1_ij beta_j^*. The stored real ket uses the transpose.
      hard%number_of_proj = nb
      hard%has_so = pot%has_so
      ALLOCATE (hard%r, SOURCE=pot%r)
      ALLOCATE (hard%rab, SOURCE=pot%rab)
      ALLOCATE (hard%lbeta, SOURCE=pot%lbeta)
      ALLOCATE (hard%beta, SOURCE=MATMUL(pot%beta, REAL(TRANSPOSE(rotation), dp)))
      ALLOCATE (hard%dion, SOURCE=pot%dion - radial_kinetic)
      CALL upf_projector_init(hard, degree, atomic%hard)
      CALL compact_projector(atomic%hard)
      ALLOCATE (atomic%kinetic(atomic%hard%nproj, atomic%hard%nproj), SOURCE=0.0_dp)
      ALLOCATE (atomic%metric(atomic%hard%nproj, atomic%hard%nproj), SOURCE=0.0_dp)
      ALLOCATE (atomic%wave_overlap(atomic%hard%nproj, atomic%hard%nproj), SOURCE=0.0_dp)
      ALLOCATE (atomic%difference_kinetic(atomic%hard%nproj, atomic%hard%nproj), SOURCE=0.0_dp)
      ALLOCATE (atomic%difference_overlap(atomic%hard%nproj, atomic%hard%nproj), SOURCE=0.0_dp)
      CPASSERT(ALLOCATED(pot%qij))
      CPASSERT(ALL(SHAPE(pot%qij) == [nb, nb]))
      DO q = 1, atomic%hard%nproj
         j = atomic%hard%radial_index(q)
         DO p = 1, atomic%hard%nproj
            IF (atomic%hard%angular_index(p) /= atomic%hard%angular_index(q)) CYCLE
            i = atomic%hard%radial_index(p)
            atomic%kinetic(p, q) = radial_kinetic(i, j)
            atomic%difference_kinetic(p, q) = delta_kinetic(i, j)
            atomic%difference_overlap(p, q) = delta_overlap(i, j)
            ! The Gaussian norm already contains the hard-wave norm within local closure.
            norm_difference = SUM(radial_weights*( &
                                  pot%partial_waves%ae_wfc(:, i)*pot%partial_waves%ae_wfc(:, j) - &
                                  pot%partial_waves%ps_wfc(:, i)*pot%partial_waves%ps_wfc(:, j)))
            atomic%metric(p, q) = pot%qij(i, j) - norm_difference
            atomic%wave_overlap(p, q) = norm_difference
         END DO
      END DO
   END SUBROUTINE upf_inverse_atomic_init

! **************************************************************************************************
!> \brief Drop exactly zero trailing quadrature rows without changing radial weights or support.
!> \param projector Prepared radial integrals, including their original quadrature weights.
! **************************************************************************************************
   SUBROUTINE compact_projector(projector)
      TYPE(upf_projector_type), INTENT(INOUT)            :: projector

      INTEGER                                            :: last
      REAL(dp), ALLOCATABLE                              :: radial(:), weighted(:, :)

      last = SIZE(projector%r)
      DO WHILE (last > 1)
         IF (ANY(projector%weighted_beta(last, :) /= 0.0_dp)) EXIT
         last = last - 1
      END DO
      IF (last == SIZE(projector%r)) RETURN
      ALLOCATE (radial, SOURCE=projector%r(:last))
      ALLOCATE (weighted, SOURCE=projector%weighted_beta(:last, :))
      CALL MOVE_ALLOC(radial, projector%r)
      CALL MOVE_ALLOC(weighted, projector%weighted_beta)
   END SUBROUTINE compact_projector

! **************************************************************************************************
!> \brief Full radial-channel kinetic matrix from the same interpolated hard/soft reference waves.
!>        The spherical-harmonic addition theorem integrates the angular gradient products
!>        analytically. Off-diagonal channels are retained, with one half in Hartree units.
!>        Only the compact difference contributes, so reference-wave tails cancel exactly.
!> \param pot Original radial mesh, quadrature and angular channel order.
!> \param kinetic Compatible reference-wave interpolants.
!> \param matrix Owned symmetric radial-channel T_hard-T_soft, without magnetic multiplicity.
!> \param difference_kinetic Optional owned Delta kinetic matrix, required with difference_overlap.
!> \param difference_overlap Optional owned Delta norm matrix, required with difference_kinetic.
! **************************************************************************************************
   SUBROUTINE upf_inverse_kinetic_matrix(pot, kinetic, matrix, difference_kinetic, difference_overlap)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: matrix(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL       :: difference_kinetic(:, :), &
                                                            difference_overlap(:, :)

      INTEGER                                            :: i, ir, j, l, m, nb, p, q
      INTEGER, ALLOCATABLE                               :: first(:)
      REAL(dp)                                           :: factor
      REAL(dp), ALLOCATABLE                              :: gradient(:, :), kernel(:, :), value(:), &
                                                            weights(:)

      nb = SIZE(pot%lbeta)
      CPASSERT(kinetic%delta%nproj == SUM(2*pot%lbeta + 1))
      CPASSERT(ALL(kinetic%delta%lbeta == pot%lbeta))
      ALLOCATE (matrix(nb, nb), SOURCE=0.0_dp)
      ALLOCATE (kernel(kinetic%delta%nproj, kinetic%delta%nproj), first(nb))
      ALLOCATE (weights, SOURCE=upf_radial_weights(pot)*pot%r**2)
      CPASSERT(PRESENT(difference_kinetic) .EQV. PRESENT(difference_overlap))
      IF (PRESENT(difference_kinetic)) THEN
         ALLOCATE (difference_kinetic(nb, nb), difference_overlap(nb, nb), &
                   value(kinetic%delta%nproj), gradient(3, kinetic%delta%nproj))
         difference_kinetic = 0.0_dp
         difference_overlap = 0.0_dp
      END IF
      p = 1
      DO i = 1, nb
         first(i) = p
         p = p + 2*pot%lbeta(i) + 1
      END DO
      DO ir = 1, SIZE(pot%r)
         IF (weights(ir) == 0.0_dp) CYCLE
         IF (pot%r(ir) >= kinetic%delta%radius) EXIT
         CALL upf_inverse_kinetic_evaluate(kinetic, [0.0_dp, 0.0_dp, pot%r(ir)], kernel)
         IF (PRESENT(difference_kinetic)) THEN
            CALL upf_radial_fields_evaluate(kinetic%delta, [0.0_dp, 0.0_dp, pot%r(ir)], value, gradient)
         END IF
         DO j = 1, nb
            l = pot%lbeta(j)
            factor = fourpi*weights(ir)/REAL(2*l + 1, dp)
            DO i = 1, j
               IF (pot%lbeta(i) /= l) CYCLE
               DO m = 0, 2*l
                  p = first(i) + m
                  q = first(j) + m
                  matrix(i, j) = matrix(i, j) + factor*kernel(p, q)
                  IF (PRESENT(difference_kinetic)) THEN
                     difference_kinetic(i, j) = difference_kinetic(i, j) + 0.5_dp*factor*DOT_PRODUCT(gradient(:, p), gradient(:, q))
                     difference_overlap(i, j) = difference_overlap(i, j) + factor*value(p)*value(q)
                  END IF
               END DO
            END DO
         END DO
      END DO
      DO j = 1, nb
         DO i = 1, j - 1
            matrix(j, i) = matrix(i, j)
            IF (PRESENT(difference_kinetic)) THEN
               difference_kinetic(j, i) = difference_kinetic(i, j)
               difference_overlap(j, i) = difference_overlap(i, j)
            END IF
         END DO
      END DO
   END SUBROUTINE upf_inverse_kinetic_matrix

! **************************************************************************************************
!> \brief Finite partial-wave closure residuals added to Gaussian T plus hard D-T and Gaussian S.
!>        These vanish only when the local hard-wave expansion satisfies the corresponding
!>        operator identity. Nonoverlap alone does not remove them on an arbitrary AO space.
!> \param atomic Original one-center matrices and transformed hard projectors.
!> \param coefficients Hard projector bras against Gaussian AO columns.
!> \param delta_overlap Delta bras against the same Gaussian columns.
!> \param delta_kinetic Delta bras against the Gaussian kinetic operator.
!> \param kinetic_residual Owned AO correction beyond hard D-T.
!> \param metric_residual Owned AO correction beyond the original q minus wave-norm residual.
! **************************************************************************************************
   SUBROUTINE upf_inverse_atomic_residual(atomic, coefficients, delta_overlap, delta_kinetic, &
                                          kinetic_residual, metric_residual)
      TYPE(upf_inverse_atomic_type), INTENT(IN)          :: atomic
      COMPLEX(dp), INTENT(IN)                            :: coefficients(:, :), delta_overlap(:, :), &
                                                            delta_kinetic(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: kinetic_residual(:, :), &
                                                            metric_residual(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: mixed(:, :)
      INTEGER                                            :: n, np

      np = atomic%hard%nproj
      n = SIZE(coefficients, 2)
      CPASSERT(SIZE(coefficients, 1) == np)
      CPASSERT(ALL(SHAPE(delta_overlap) == [np, n]))
      CPASSERT(ALL(SHAPE(delta_kinetic) == [np, n]))
      ALLOCATE (mixed, SOURCE=MATMUL(CONJG(TRANSPOSE(delta_kinetic)), coefficients))
      ALLOCATE (kinetic_residual, SOURCE=-mixed - CONJG(TRANSPOSE(mixed)) + &
                MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(atomic%difference_kinetic + atomic%kinetic, coefficients)))
      mixed(:, :) = MATMUL(CONJG(TRANSPOSE(delta_overlap)), coefficients)
      ALLOCATE (metric_residual, SOURCE=-mixed - CONJG(TRANSPOSE(mixed)) + &
                MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(atomic%difference_overlap + atomic%wave_overlap, coefficients)))
   END SUBROUTINE upf_inverse_atomic_residual

! **************************************************************************************************
!> \brief Independent valence kinetic augmentation from original projector-indexed AE/PS pairs.
!>        Charge Q is not used. The noncompact PS table is retained over its full supplied mesh.
!>        Only the compact AE-minus-PS difference is evaluated outside the reconstruction region.
!>        This kernel does not supply a frozen-core kinetic density or a pseudization of tau.
!> \param pot Original scalar dataset with compatible full-wave pairs.
!> \param kinetic New radial caches in the original magnetic-projector order.
! **************************************************************************************************
   SUBROUTINE upf_inverse_kinetic_init(pot, kinetic)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_inverse_kinetic_type), INTENT(OUT)        :: kinetic

      INTEGER                                            :: i
      TYPE(atom_upfpot_type)                             :: radial

      CALL inverse_radial_data(pot, radial)
      CALL upf_radial_fields_init(radial%r, radial%beta, radial%lbeta, kinetic%delta)
      DO i = 1, pot%number_of_proj
         IF (pot%r(1) == 0.0_dp) THEN
            radial%beta(2:, i) = pot%partial_waves%ps_wfc(2:, i)/pot%r(2:)
            radial%beta(1, i) = 0.0_dp
         ELSE
            radial%beta(:, i) = pot%partial_waves%ps_wfc(:, i)/pot%r
         END IF
      END DO
      CALL upf_radial_fields_init(radial%r, radial%beta, radial%lbeta, kinetic%pseudo, finite_table=.TRUE.)
   END SUBROUTINE upf_inverse_kinetic_init

! **************************************************************************************************
!> \brief One-half grad(AE_i).grad(AE_j) minus one-half grad(PS_i).grad(PS_j).
!>        Factoring through Delta preserves an exact zero outside its compact support.
!>        Values and spatial derivatives use the same interpolated reference waves.
!> \param kinetic Original reference-wave caches, independent of the supplied charge Q.
!> \param displacement Point minus original center, in bohr.
!> \param matrix Symmetric kinetic augmentation in the native magnetic-projector order.
!> \param gradient Optional Cartesian point gradient. Center motion has the opposite sign.
! **************************************************************************************************
   SUBROUTINE upf_inverse_kinetic_evaluate(kinetic, displacement, matrix, gradient)
      TYPE(upf_inverse_kinetic_type), INTENT(IN)         :: kinetic
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: matrix(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gradient(:, :, :)

      INTEGER                                            :: axis, i, j, n
      REAL(dp), ALLOCATABLE                              :: d(:, :), dd(:, :, :), p(:, :), &
                                                            pd(:, :, :), value(:)

      n = kinetic%delta%nproj
      CPASSERT(kinetic%pseudo%nproj == n .AND. ALL(SHAPE(matrix) == [n, n]))
      matrix = 0.0_dp
      IF (PRESENT(gradient)) THEN
         CPASSERT(ALL(SHAPE(gradient) == [n, n, 3]))
         gradient = 0.0_dp
      END IF
      IF (SUM(displacement**2) >= kinetic%delta%radius**2) RETURN
      ALLOCATE (value(n), d(3, n), p(3, n))
      IF (PRESENT(gradient)) THEN
         ALLOCATE (dd(3, 3, n), pd(3, 3, n))
         CALL upf_radial_fields_evaluate(kinetic%delta, displacement, value, d, dd)
         CALL upf_radial_fields_evaluate(kinetic%pseudo, displacement, value, p, pd)
      ELSE
         CALL upf_radial_fields_evaluate(kinetic%delta, displacement, value, d)
         CALL upf_radial_fields_evaluate(kinetic%pseudo, displacement, value, p)
      END IF
      DO j = 1, n
         DO i = 1, j
            matrix(i, j) = 0.5_dp*(DOT_PRODUCT(d(:, i), d(:, j)) + &
                                   DOT_PRODUCT(d(:, i), p(:, j)) + DOT_PRODUCT(p(:, i), d(:, j)))
            matrix(j, i) = matrix(i, j)
            IF (PRESENT(gradient)) THEN
               DO axis = 1, 3
                  gradient(i, j, axis) = 0.5_dp*(DOT_PRODUCT(dd(:, axis, i), d(:, j) + p(:, j)) + &
                                                 DOT_PRODUCT(d(:, i) + p(:, i), dd(:, axis, j)) + &
                                                 DOT_PRODUCT(pd(:, axis, i), d(:, j)) + &
                                                 DOT_PRODUCT(d(:, i), pd(:, axis, j)))
                  gradient(j, i, axis) = gradient(i, j, axis)
               END DO
            END IF
         END DO
      END DO
   END SUBROUTINE upf_inverse_kinetic_evaluate

! **************************************************************************************************
!> \brief Prepare original AE-minus-PS radial overlaps in the native projector/harmonic order.
!>        The ordinary radial/Gaussian integrator supplies overlaps and AO-center derivatives.
!>        This is an overlap cache, not a replacement for the original nonlocal operator.
!> \param pot Original data with a complete full-wave pair for every beta projector.
!> \param degree Required Cartesian degree, including derivatives or moments.
!> \param delta Original radial quadrature and angular maps for (u_AE-u_PS)/r.
! **************************************************************************************************
   SUBROUTINE upf_inverse_delta_init(pot, degree, delta)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: degree
      TYPE(upf_projector_type), INTENT(OUT)              :: delta

      TYPE(atom_upfpot_type)                             :: radial

      CALL inverse_radial_data(pot, radial)
      CALL upf_projector_init(radial, degree, delta)
   END SUBROUTINE upf_inverse_delta_init

! **************************************************************************************************
!> \brief Prepare original beta and reconstruction fields for common intercenter quadratures.
!>        Native radial channels and real harmonics agree with the Gaussian overlap cache.
!>        NC channels have zero reconstruction but retain their original projector bras.
!> \param pot Original UPF data, with complete full waves for reconstructing channels.
!> \param beta Original projector fields, without weights.
!> \param delta Original reconstruction differences, without weights.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE upf_inverse_fields_init(pot, beta, delta, direct)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(upf_radial_fields_type), INTENT(OUT)          :: beta, delta
      LOGICAL, INTENT(IN), OPTIONAL                      :: direct

      TYPE(atom_upfpot_type)                             :: radial

      CALL upf_radial_fields_init(pot%r, pot%beta, pot%lbeta, beta)
      IF (PRESENT(direct)) THEN
         IF (direct) THEN
            ALLOCATE (radial%beta(SIZE(pot%r), SIZE(pot%lbeta)), SOURCE=0.0_dp)
            CALL upf_radial_fields_init(pot%r, radial%beta, pot%lbeta, delta)
            RETURN
         END IF
      END IF
      IF (TRIM(pot%pseudo_type) == "NC" .AND. .NOT. pot%is_ultrasoft .AND. .NOT. pot%is_paw) THEN
         IF (pot%has_so) CPABORT("Spin-orbit UPF inverse fields are not supported")
         ALLOCATE (radial%beta(SIZE(pot%r), SIZE(pot%lbeta)))
         radial%beta = 0.0_dp
         CALL upf_radial_fields_init(pot%r, radial%beta, pot%lbeta, delta)
      ELSE
         CALL inverse_radial_data(pot, radial)
         CALL upf_radial_fields_init(radial%r, radial%beta, radial%lbeta, delta)
      END IF
   END SUBROUTINE upf_inverse_fields_init

! **************************************************************************************************
!> \brief Assemble original reconstruction differences for overlap and point-field caches.
!> \param pot Original UPF data.
!> \param radial Unweighted original differences in the native radial-channel order.
! **************************************************************************************************
   SUBROUTINE inverse_radial_data(pot, radial)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      TYPE(atom_upfpot_type), INTENT(OUT)                :: radial

      INTEGER                                            :: i, n, nb

      IF (pot%has_so) THEN
         CPABORT("Spin-orbit UPF inverse fields are not supported")
      END IF
      IF (.NOT. ALLOCATED(pot%partial_waves)) THEN
         CPABORT("Inverse UPF requires original AE/PS full waves for every projector")
      END IF
      CPASSERT(ALLOCATED(pot%r) .AND. ALLOCATED(pot%rab))
      CPASSERT(ALLOCATED(pot%lbeta))
      IF (.NOT. (ALLOCATED(pot%partial_waves%ae_wfc) .AND. ALLOCATED(pot%partial_waves%ps_wfc))) THEN
         CPABORT("Incomplete original AE/PS full waves")
      END IF
      CPASSERT(ALLOCATED(pot%partial_waves%angular_momenta))
      n = SIZE(pot%r)
      nb = pot%number_of_proj
      CPASSERT(n >= 3 .AND. nb > 0 .AND. ALL(pot%r >= 0.0_dp))
      CPASSERT(SIZE(pot%rab) == n .AND. SIZE(pot%lbeta) == nb)
      CPASSERT(ALL(SHAPE(pot%partial_waves%ae_wfc) == [n, nb]))
      CPASSERT(ALL(SHAPE(pot%partial_waves%ps_wfc) == [n, nb]))
      CPASSERT(ALL(pot%partial_waves%angular_momenta == pot%lbeta))
      CPASSERT(ALL(ieee_is_finite(pot%partial_waves%ae_wfc)))
      CPASSERT(ALL(ieee_is_finite(pot%partial_waves%ps_wfc)))
      radial%number_of_proj = nb
      radial%has_so = pot%has_so
      ALLOCATE (radial%r, SOURCE=pot%r)
      ALLOCATE (radial%rab, SOURCE=pot%rab)
      ALLOCATE (radial%lbeta, SOURCE=pot%lbeta)
      ALLOCATE (radial%beta(n, nb), radial%dion(nb, nb))
      radial%dion(:, :) = 0.0_dp
      DO i = 1, nb
         radial%beta(:, i) = pot%partial_waves%ae_wfc(:, i) - pot%partial_waves%ps_wfc(:, i)
         IF (pot%r(1) == 0.0_dp) THEN
            radial%beta(2:, i) = radial%beta(2:, i)/pot%r(2:)
            ! Zero quadrature weight at r=0; this storage does not define a field continuation.
            radial%beta(1, i) = 0.0_dp
         ELSE
            radial%beta(:, i) = radial%beta(:, i)/pot%r
         END IF
      END DO
   END SUBROUTINE inverse_radial_data

! **************************************************************************************************
!> \brief Factor K with original projector normalization and all supplied intercenter blocks.
!> \param beta_delta P*Delta, including periodic images and Bloch phases.
!> \param inverse New factors; a failed construction cannot reuse an earlier factorization.
!> \param info Zero on success, LAPACK pivot on singularity, nproj+1 for poor conditioning,
!>             nproj+2 for nonfinite data. Failed factors must not be used.
!> \param minimum_rcond Optional condition limit; no regularization changes the reconstruction.
! **************************************************************************************************
   SUBROUTINE upf_inverse_init(beta_delta, inverse, info, minimum_rcond)
      COMPLEX(KIND=dp), INTENT(IN)                       :: beta_delta(:, :)
      TYPE(upf_inverse_type), INTENT(OUT)                :: inverse
      INTEGER, INTENT(OUT)                               :: info
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: minimum_rcond

      COMPLEX(KIND=dp), ALLOCATABLE                      :: work(:)
      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: anorm, limit
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:)
      REAL(KIND=dp), EXTERNAL                            :: zlange

      EXTERNAL zgecon, zgetrf

      n = SIZE(beta_delta, 1)
      CPASSERT(SIZE(beta_delta, 2) == n)
      limit = 64.0_dp*EPSILON(0.0_dp)
      IF (PRESENT(minimum_rcond)) limit = minimum_rcond
      CPASSERT(ieee_is_finite(limit) .AND. limit >= 0.0_dp .AND. limit <= 1.0_dp)
      inverse%nproj = n
      info = 0
      IF (.NOT. ALL(ieee_is_finite(REAL(beta_delta, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(beta_delta)))) THEN
         info = n + 2
         RETURN
      END IF
      ALLOCATE (inverse%lu(n, n), inverse%pivots(n))
      inverse%lu(:, :) = beta_delta
      DO i = 1, n
         inverse%lu(i, i) = inverse%lu(i, i) + 1.0_dp
      END DO
      IF (n == 0) THEN
         inverse%rcond = 1.0_dp
         inverse%ready = .TRUE.
         RETURN
      END IF
      ALLOCATE (work(2*n), rwork(2*n))
      anorm = zlange('1', n, n, inverse%lu, n, rwork)
      CALL zgetrf(n, n, inverse%lu, n, inverse%pivots, info)
      IF (info /= 0) RETURN
      CALL zgecon('1', n, inverse%lu, n, anorm, inverse%rcond, work, rwork, info)
      IF (info /= 0) RETURN
      IF (.NOT. ieee_is_finite(inverse%rcond)) THEN
         info = n + 2
         RETURN
      END IF
      IF (inverse%rcond < limit) THEN
         info = n + 1
         RETURN
      END IF
      inverse%ready = .TRUE.
   END SUBROUTINE upf_inverse_init

! **************************************************************************************************
!> \brief Solve for W=P*U*chi; do not replace original q by an orthogonalized metric.
!> \param inverse ...
!> \param beta_gaussian B=P*chi, with one Gaussian function per column.
!> \param coefficients W, using the already factored full K.
! **************************************************************************************************
   SUBROUTINE upf_inverse_project(inverse, beta_gaussian, coefficients)
      TYPE(upf_inverse_type), INTENT(IN)                 :: inverse
      COMPLEX(KIND=dp), INTENT(IN)                       :: beta_gaussian(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: coefficients(:, :)

      INTEGER                                            :: info, n, nrhs

      EXTERNAL zgetrs

      CPASSERT(inverse%ready)
      n = inverse%nproj
      nrhs = SIZE(beta_gaussian, 2)
      CPASSERT(SIZE(beta_gaussian, 1) == n)
      CPASSERT(ALL(SHAPE(coefficients) == SHAPE(beta_gaussian)))
      coefficients = beta_gaussian
      IF (n == 0 .OR. nrhs == 0) RETURN
      CALL zgetrs('N', n, nrhs, inverse%lu, n, inverse%pivots, coefficients, n, info)
      CPASSERT(info == 0)
   END SUBROUTINE upf_inverse_project

! **************************************************************************************************
!> \brief Differentiate the same reconstruction: K*dW=dB-dK*W.
!>        Displacement/strain responses must include the original waves, projectors,
!>        their moving integration measure, and all relevant Bloch-image terms.
!> \param inverse ...
!> \param coefficients ...
!> \param beta_gaussian_response ...
!> \param beta_delta_response ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_project_response(inverse, coefficients, beta_gaussian_response, &
                                           beta_delta_response, response)
      TYPE(upf_inverse_type), INTENT(IN)                 :: inverse
      COMPLEX(KIND=dp), INTENT(IN)                       :: coefficients(:, :), &
                                                            beta_gaussian_response(:, :), &
                                                            beta_delta_response(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: response(:, :)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: rhs(:, :)
      INTEGER                                            :: n

      n = inverse%nproj
      CPASSERT(SIZE(coefficients, 1) == n)
      CPASSERT(ALL(SHAPE(beta_delta_response) == [n, n]))
      CPASSERT(ALL(SHAPE(beta_gaussian_response) == SHAPE(coefficients)))
      ALLOCATE (rhs(n, SIZE(coefficients, 2)))
      rhs(:, :) = beta_gaussian_response - MATMUL(beta_delta_response, coefficients)
      CALL upf_inverse_project(inverse, rhs, response)
   END SUBROUTINE upf_inverse_project_response

! **************************************************************************************************
!> \brief Evaluate inverse images on a field batch: U*chi=chi-Delta*W.
!>        Spatial derivatives use the corresponding derivative fields, with W fixed.
!> \param gaussian ...
!> \param delta ...
!> \param coefficients ...
!> \param soft ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_basis(gaussian, delta, coefficients, soft)
      COMPLEX(KIND=dp), INTENT(IN)                       :: gaussian(:, :), delta(:, :), &
                                                            coefficients(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: soft(:, :)

      CPASSERT(SIZE(delta, 1) == SIZE(gaussian, 1))
      CPASSERT(SIZE(delta, 2) == SIZE(coefficients, 1))
      CPASSERT(SIZE(coefficients, 2) == SIZE(gaussian, 2))
      CPASSERT(ALL(SHAPE(soft) == SHAPE(gaussian)))
      soft = gaussian - MATMUL(delta, coefficients)
   END SUBROUTINE upf_inverse_basis

! **************************************************************************************************
!> \brief Joint field response d(U*chi)=dchi-dDelta*W-Delta*dW.
!> \param gaussian_response ...
!> \param delta ...
!> \param delta_response ...
!> \param coefficients ...
!> \param coefficients_response ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_basis_response(gaussian_response, delta, delta_response, &
                                         coefficients, coefficients_response, response)
      COMPLEX(KIND=dp), INTENT(IN) :: gaussian_response(:, :), delta(:, :), delta_response(:, :), &
         coefficients(:, :), coefficients_response(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: response(:, :)

      CPASSERT(ALL(SHAPE(delta_response) == SHAPE(delta)))
      CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(coefficients)))
      CALL upf_inverse_basis(gaussian_response, delta_response, coefficients, response)
      response = response - MATMUL(delta, coefficients_response)
   END SUBROUTINE upf_inverse_basis_response

! **************************************************************************************************
!> \brief Pull back A+P^dagger*augmentation*P into reconstructed Gaussian variables.
!> \param gaussian_operator <chi|A|chi>.
!> \param gaussian_delta <chi|A|Delta>.
!> \param delta_gaussian <Delta|A|chi>; supplied separately for non-Hermitian operators.
!> \param delta_operator <Delta|A|Delta>.
!> \param augmentation Original projector-space contribution; q for the metric,
!>                     Dij+integral(v*Q) for nonlocal/local-augmentation terms.
!> \param coefficients W=K^-1*B. Original operators are retained, including residuals.
!> \param transformed <U*chi|A|U*chi>+W^dagger*augmentation*W.
! **************************************************************************************************
   SUBROUTINE upf_inverse_operator(gaussian_operator, gaussian_delta, delta_gaussian, delta_operator, &
                                   augmentation, coefficients, transformed)
      COMPLEX(KIND=dp), INTENT(IN) :: gaussian_operator(:, :), gaussian_delta(:, :), &
         delta_gaussian(:, :), delta_operator(:, :), augmentation(:, :), coefficients(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: transformed(:, :)

      CALL check_operator_shapes(gaussian_operator, gaussian_delta, delta_gaussian, delta_operator, &
                                 augmentation, coefficients, SHAPE(transformed))
      transformed = gaussian_operator - MATMUL(gaussian_delta, coefficients) - &
                    MATMUL(CONJG(TRANSPOSE(coefficients)), delta_gaussian) + &
                    MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(delta_operator + augmentation, coefficients))
   END SUBROUTINE upf_inverse_operator

! **************************************************************************************************
!> \brief Differentiate the full operator pullback, including both inverse-basis legs.
!>        Matrix responses already include their moving original fields and quadratures.
!>        The stationary energy also requires the transformed-metric Pulay term.
!> \param gaussian_response ...
!> \param gaussian_delta ...
!> \param delta_gaussian ...
!> \param delta_operator ...
!> \param augmentation ...
!> \param coefficients ...
!> \param gaussian_delta_response ...
!> \param delta_gaussian_response ...
!> \param delta_operator_response ...
!> \param augmentation_response ...
!> \param coefficients_response ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE upf_inverse_operator_response(gaussian_response, gaussian_delta, delta_gaussian, &
                                            delta_operator, augmentation, coefficients, &
                                            gaussian_delta_response, delta_gaussian_response, &
                                            delta_operator_response, augmentation_response, &
                                            coefficients_response, response)
      COMPLEX(KIND=dp), INTENT(IN) :: gaussian_response(:, :), gaussian_delta(:, :), &
         delta_gaussian(:, :), delta_operator(:, :), augmentation(:, :), coefficients(:, :), &
         gaussian_delta_response(:, :), delta_gaussian_response(:, :), &
         delta_operator_response(:, :), augmentation_response(:, :), coefficients_response(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: response(:, :)

      CALL check_operator_shapes(gaussian_response, gaussian_delta, delta_gaussian, delta_operator, &
                                 augmentation, coefficients, SHAPE(response))
      CPASSERT(ALL(SHAPE(coefficients_response) == SHAPE(coefficients)))
      CALL upf_inverse_operator(gaussian_response, gaussian_delta_response, delta_gaussian_response, &
                                delta_operator_response, augmentation_response, coefficients, response)
      response = response - MATMUL(gaussian_delta, coefficients_response) - &
                 MATMUL(CONJG(TRANSPOSE(coefficients_response)), delta_gaussian) + &
                 MATMUL(CONJG(TRANSPOSE(coefficients_response)), MATMUL(delta_operator + augmentation, coefficients)) + &
                 MATMUL(CONJG(TRANSPOSE(coefficients)), MATMUL(delta_operator + augmentation, coefficients_response))
   END SUBROUTINE upf_inverse_operator_response

! **************************************************************************************************
!> \brief ...
!> \param gg ...
!> \param gd ...
!> \param dg ...
!> \param dd ...
!> \param augmentation ...
!> \param coefficients ...
!> \param output_shape ...
! **************************************************************************************************
   SUBROUTINE check_operator_shapes(gg, gd, dg, dd, augmentation, coefficients, output_shape)
      COMPLEX(KIND=dp), INTENT(IN)                       :: gg(:, :), gd(:, :), dg(:, :), dd(:, :), &
                                                            augmentation(:, :), coefficients(:, :)
      INTEGER, INTENT(IN)                                :: output_shape(2)

      INTEGER                                            :: m, n

      m = SIZE(coefficients, 1)
      n = SIZE(coefficients, 2)
      CPASSERT(ALL(SHAPE(gg) == [n, n]))
      CPASSERT(ALL(SHAPE(gd) == [n, m]))
      CPASSERT(ALL(SHAPE(dg) == [m, n]))
      CPASSERT(ALL(SHAPE(dd) == [m, m]))
      CPASSERT(ALL(SHAPE(augmentation) == [m, m]))
      CPASSERT(ALL(output_shape == [n, n]))
   END SUBROUTINE check_operator_shapes

END MODULE upf_inverse
