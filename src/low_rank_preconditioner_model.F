!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Small algebraic helpers for the rotationally covariant low-rank OT preconditioner.
! **************************************************************************************************
MODULE low_rank_preconditioner_model
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: apply_low_rank_dense, low_rank_inverse_weight, low_rank_select_rank

CONTAINS

! **************************************************************************************************
!> \brief Spectral inverse weight relative to a common occupied reference level.
!> \param eigenvalue complementary-state eigenvalue
!> \param reference common occupied reference level
!> \param energy_gap lower bound for the denominator
!> \return inverse spectral weight
! **************************************************************************************************
   PURE ELEMENTAL REAL(KIND=dp) FUNCTION low_rank_inverse_weight( &
      eigenvalue, reference, energy_gap) RESULT(weight)

      REAL(KIND=dp), INTENT(IN)                          :: eigenvalue, reference, energy_gap

      REAL(KIND=dp)                                      :: delta

      delta = eigenvalue - reference
      weight = 1.0_dp/MAX(energy_gap, delta)

   END FUNCTION low_rank_inverse_weight

! **************************************************************************************************
!> \brief Select a bounded spectral rank without cutting a degenerate boundary manifold.
!> \param eigenvalues ordered full-space eigenvalues
!> \param nocc number of occupied states at the start of eigenvalues
!> \param max_rank hard upper bound
!> \param degeneracy_tolerance relative degeneracy threshold
!> \return selected rank
! **************************************************************************************************
   PURE INTEGER FUNCTION low_rank_select_rank(eigenvalues, nocc, max_rank, &
                                              degeneracy_tolerance) RESULT(rank_used)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      INTEGER, INTENT(IN)                                :: nocc, max_rank
      REAL(KIND=dp), INTENT(IN)                          :: degeneracy_tolerance

      INTEGER                                            :: ncomplement
      REAL(KIND=dp)                                      :: scale

      rank_used = 0
      ncomplement = SIZE(eigenvalues) - nocc
      IF (nocc < 0 .OR. ncomplement <= 0 .OR. max_rank <= 0) RETURN

      rank_used = MIN(max_rank, ncomplement)

      ! Exclude the complete boundary cluster if the hard rank cap would split it.
      DO WHILE (rank_used > 0 .AND. rank_used < ncomplement)
         scale = MAX(1.0_dp, ABS(eigenvalues(nocc + rank_used)), &
                     ABS(eigenvalues(nocc + rank_used + 1)))
         IF (ABS(eigenvalues(nocc + rank_used + 1) - eigenvalues(nocc + rank_used)) > &
             degeneracy_tolerance*scale) EXIT
         rank_used = rank_used - 1
      END DO

   END FUNCTION low_rank_select_rank

! **************************************************************************************************
!> \brief Dense reference application used to verify right covariance under orbital rotations.
!> \param overlap_inverse inverse overlap matrix
!> \param vectors retained complementary-state vectors
!> \param eigenvalues retained complementary-state eigenvalues
!> \param reference common occupied reference level
!> \param energy_gap lower bound for spectral denominators
!> \param spectral_window overlap-inverse replacement window
!> \param matrix_in input orbital-gradient matrix
!> \param matrix_out preconditioned matrix
! **************************************************************************************************
   PURE SUBROUTINE apply_low_rank_dense(overlap_inverse, vectors, eigenvalues, reference, &
                                        energy_gap, spectral_window, matrix_in, matrix_out)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: overlap_inverse, vectors
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      REAL(KIND=dp), INTENT(IN)                          :: reference, energy_gap, spectral_window
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix_in
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: matrix_out

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: correction

      matrix_out = MATMUL(overlap_inverse, matrix_in)/spectral_window
      DO i = 1, SIZE(eigenvalues)
         correction = low_rank_inverse_weight(eigenvalues(i), reference, energy_gap) - &
                      1.0_dp/spectral_window
         matrix_out = matrix_out + correction*SPREAD(vectors(:, i), 2, SIZE(matrix_in, 2))* &
                      SPREAD(MATMUL(vectors(:, i), matrix_in), 1, SIZE(matrix_in, 1))
      END DO

   END SUBROUTINE apply_low_rank_dense

END MODULE low_rank_preconditioner_model
