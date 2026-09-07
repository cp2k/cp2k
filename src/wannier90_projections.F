!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gauge-covariant selection of fixed trial orbitals for Wannier90.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE wannier90_projections
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier90_projections'
   PUBLIC :: wannier90_select_projections, wannier90_projection_rank

CONTAINS

! **************************************************************************************************
!> \brief Select the same trial orbitals at every k-point using a mesh-wide pivoted QR.
!> \param candidates Band/candidate/k-point projections, stored only on the source rank.
!> \param para_env Global parallel environment.
!> \param amn Selected projection matrices, allocated with the target dimensions on all ranks.
!> \param pivots Selected candidate indices, allocated on all ranks.
!> \param min_svalue Smallest singular value of the selected projections over the mesh.
!> \param anchor_kpoint Optional fixed QR anchor; zero selects the stacked mesh. Without an
!>                     explicit anchor, try individual k-points if the stacked selection loses rank.
! **************************************************************************************************
   SUBROUTINE wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue, anchor_kpoint)
      COMPLEX(KIND=dp), INTENT(IN)                       :: candidates(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), INTENT(OUT)                      :: amn(:, :, :)
      INTEGER, INTENT(OUT)                               :: pivots(:)
      REAL(KIND=dp), INTENT(OUT)                         :: min_svalue
      INTEGER, INTENT(IN), OPTIONAL                      :: anchor_kpoint

      COMPLEX(KIND=dp)                                   :: query(1)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: qr(:, :), tau(:), work(:)
      INTEGER :: anchor, best_defect, defect, first_anchor, ic, ik, info, iw, last_anchor, lwork, &
         nb, nc, nk, nr, nw, previous_defect, saved_pivot
      INTEGER, ALLOCATABLE                               :: best_pivots(:), order(:)
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:)

      nb = SIZE(amn, 1)
      nw = SIZE(amn, 2)
      nk = SIZE(amn, 3)
      CPASSERT(nb >= nw .AND. nw > 0 .AND. SIZE(pivots) == nw)
      first_anchor = 0
      last_anchor = nk
      IF (PRESENT(anchor_kpoint)) THEN
         CPASSERT(anchor_kpoint >= 0 .AND. anchor_kpoint <= nk)
         first_anchor = anchor_kpoint
         last_anchor = anchor_kpoint
      END IF
      min_svalue = 0.0_dp
      IF (para_env%is_source()) THEN
         nc = SIZE(candidates, 2)
         CPASSERT(SIZE(candidates, 1) == nb .AND. SIZE(candidates, 3) == nk .AND. nc >= nw)
         ALLOCATE (best_pivots(nw))
         best_defect = HUGE(0)
         DO anchor = first_anchor, last_anchor
            nr = nb*nk
            IF (anchor > 0) nr = nb
            ALLOCATE (qr(nr, nc), tau(MIN(nr, nc)), order(nc), rwork(2*nc))
            IF (anchor == 0) THEN
            DO ik = 1, nk
               qr((ik - 1)*nb + 1:ik*nb, :) = candidates(:, :, ik)
            END DO
            ELSE
            qr(:, :) = candidates(:, :, anchor)
            END IF
            ! A unitary change of the MO gauge acts on rows, not on the trial orbitals.
            ! The selected columns therefore represent one physical trial set on the full mesh.
            order(:) = 0
            CALL zgeqp3(nr, nc, qr, nr, order, tau, query, -1, rwork, info)
            CPASSERT(info == 0)
            lwork = MAX(1, INT(REAL(query(1), KIND=dp)))
            ALLOCATE (work(lwork))
            CALL zgeqp3(nr, nc, qr, nr, order, tau, work, lwork, rwork, info)
            CPASSERT(info == 0)
            pivots(:) = order(1:nw)
            DO ik = 1, nk
               amn(:, :, ik) = candidates(:, pivots, ik)
            END DO
            DEALLOCATE (qr, tau, order, rwork, work)

            min_svalue = wannier90_projection_rank(amn, defect)
            IF (defect < best_defect) THEN
               best_defect = defect
               best_pivots = pivots
            END IF
            IF (min_svalue >= 1.0e-10_dp) EXIT
         END DO
         IF (best_defect > 0 .AND. .NOT. PRESENT(anchor_kpoint)) THEN
            ! QR optimizes one anchor, not rank everywhere. Each accepted exchange
            ! strictly reduces the number of missing trial directions across the mesh.
            pivots = best_pivots
            DO
               previous_defect = best_defect
               repair: DO iw = 1, nw
                  saved_pivot = pivots(iw)
                  DO ic = 1, nc
                     IF (ANY(pivots == ic)) CYCLE
                     pivots(iw) = ic
                     DO ik = 1, nk
                        amn(:, :, ik) = candidates(:, pivots, ik)
                     END DO
                     min_svalue = wannier90_projection_rank(amn, defect)
                     IF (defect < best_defect) THEN
                        best_defect = defect
                        EXIT repair
                     END IF
                     pivots(iw) = saved_pivot
                  END DO
               END DO repair
               IF (best_defect == 0 .OR. best_defect == previous_defect) EXIT
            END DO
            DO ik = 1, nk
               amn(:, :, ik) = candidates(:, pivots, ik)
            END DO
            min_svalue = wannier90_projection_rank(amn)
         END IF
         DEALLOCATE (best_pivots)
      END IF
      CALL para_env%bcast(amn)
      CALL para_env%bcast(pivots)
      CALL para_env%bcast(min_svalue)

   END SUBROUTINE wannier90_select_projections

! **************************************************************************************************
!> \brief Validate the trial rank independently at every k-point.
!> \param amn Band/trial/k-point projections, possibly rectangular.
!> \param rank_defect Total number of missing trial directions summed over the mesh.
!> \return Smallest singular value over all k-points.
! **************************************************************************************************
   FUNCTION wannier90_projection_rank(amn, rank_defect) RESULT(min_svalue)
      COMPLEX(KIND=dp), INTENT(IN)                       :: amn(:, :, :)
      INTEGER, INTENT(OUT), OPTIONAL                     :: rank_defect
      REAL(KIND=dp)                                      :: min_svalue

      COMPLEX(KIND=dp)                                   :: query(1), unused_u(1, 1), unused_v(1, 1)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: matrix(:, :), work(:)
      INTEGER                                            :: ik, info, lwork, nb, nw
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:), svalues(:)

      nb = SIZE(amn, 1)
      nw = SIZE(amn, 2)
      CPASSERT(nb >= nw .AND. nw > 0 .AND. SIZE(amn, 3) > 0)
      ALLOCATE (matrix(nb, nw), svalues(nw), rwork(5*nw))
      matrix = amn(:, :, 1)
      CALL zgesvd('N', 'N', nb, nw, matrix, nb, svalues, unused_u, 1, unused_v, 1, query, -1, rwork, info)
      CPASSERT(info == 0)
      lwork = MAX(1, INT(REAL(query(1), KIND=dp)))
      ALLOCATE (work(lwork))
      min_svalue = HUGE(0.0_dp)
      IF (PRESENT(rank_defect)) rank_defect = 0
      DO ik = 1, SIZE(amn, 3)
         matrix = amn(:, :, ik)
         CALL zgesvd('N', 'N', nb, nw, matrix, nb, svalues, unused_u, 1, unused_v, 1, work, lwork, rwork, info)
         CPASSERT(info == 0)
         min_svalue = MIN(min_svalue, MINVAL(svalues))
         IF (PRESENT(rank_defect)) rank_defect = rank_defect + COUNT(svalues < 1.0e-10_dp)
      END DO
      DEALLOCATE (matrix, work, rwork, svalues)
   END FUNCTION wannier90_projection_rank

END MODULE wannier90_projections
