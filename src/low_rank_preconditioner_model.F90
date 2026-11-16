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

   PUBLIC :: low_rank_inverse_weight, low_rank_select_rank

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

END MODULE low_rank_preconditioner_model
