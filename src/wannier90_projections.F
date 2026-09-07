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
   PUBLIC :: wannier90_select_projections, wannier90_projection_rank, wannier90_hybrid_candidates

CONTAINS

! **************************************************************************************************
!> \brief Form tetrahedral s/p trials from the full AO bank before selecting projections.
!> \param ao Band/AO/k-point projections on the source rank, empty elsewhere.
!> \param groups AO indices of one s and three p functions for each radial shell pair.
!> \param candidates Hybrid trials and AOs not belonging to any hybrid group.
! **************************************************************************************************
   SUBROUTINE wannier90_hybrid_candidates(ao, groups, candidates)
      COMPLEX(KIND=dp), INTENT(IN)                       :: ao(:, :, :)
      INTEGER, INTENT(IN)                                :: groups(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: candidates(:, :, :)

      REAL(KIND=dp), PARAMETER :: hybrid(4, 4) = 0.5_dp*RESHAPE([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, &
         -1.0_dp, 1.0_dp], [4, 4])

      INTEGER                                            :: iao, ic, ig, ik, nao, nb, nk
      LOGICAL, ALLOCATABLE                               :: grouped(:)

      nb = SIZE(ao, 1)
      nao = SIZE(ao, 2)
      nk = SIZE(ao, 3)
      IF (nb == 0) THEN
         ALLOCATE (candidates(0, 0, 0))
         RETURN
      END IF
      CPASSERT(SIZE(groups, 1) == 4)
      CPASSERT(ALL(groups > 0) .AND. ALL(groups <= nao))
      ALLOCATE (grouped(nao))
      grouped = .FALSE.
      DO ig = 1, SIZE(groups, 2)
         DO ic = 1, 4
            CPASSERT(COUNT(groups(:, ig) == groups(ic, ig)) == 1)
         END DO
         grouped(groups(:, ig)) = .TRUE.
      END DO
      ALLOCATE (candidates(nb, COUNT(.NOT. grouped) + 4*SIZE(groups, 2), nk))
      ic = 0
      DO iao = 1, nao
         IF (grouped(iao)) CYCLE
         ic = ic + 1
         candidates(:, ic, :) = ao(:, iao, :)
      END DO
      ! Include every radial s/p pair. Selecting only already-pivoted AO columns
      ! would exclude physical trials that the original QR happened to discard.
      DO ig = 1, SIZE(groups, 2)
         DO ik = 1, nk
            candidates(:, ic + 1:ic + 4, ik) = MATMUL(ao(:, groups(:, ig), ik), hybrid)
         END DO
         ic = ic + 4
      END DO
      DEALLOCATE (grouped)
   END SUBROUTINE wannier90_hybrid_candidates

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

      COMPLEX(KIND=dp), ALLOCATABLE                      :: qr(:, :)
      INTEGER                                            :: anchor, best_defect, defect, &
                                                            first_anchor, ic, ik, iw, last_anchor, &
                                                            nb, nc, nk, nr, nw, previous_defect, &
                                                            saved_pivot
      INTEGER, ALLOCATABLE                               :: best_pivots(:)

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
         CPASSERT(SIZE(candidates, 1) == nb .AND. SIZE(candidates, 3) == nk)
         CPASSERT(nc >= nw)
         ALLOCATE (best_pivots(nw))
         best_defect = HUGE(0)
         DO anchor = first_anchor, last_anchor
            nr = nb*nk
            IF (anchor > 0) nr = nb
            ALLOCATE (qr(nr, nc))
            IF (anchor == 0) THEN
            DO ik = 1, nk
               qr((ik - 1)*nb + 1:ik*nb, :) = candidates(:, :, ik)
            END DO
            ELSE
            qr(:, :) = candidates(:, :, anchor)
            END IF
            ! A unitary change of the MO gauge acts on rows, not on the trial orbitals.
            ! The selected columns therefore represent one physical trial set on the full mesh.
            CALL select_qr_pivots(qr, pivots)
            DO ik = 1, nk
               amn(:, :, ik) = candidates(:, pivots, ik)
            END DO
            DEALLOCATE (qr)

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
!> \brief Pivoted Householder QR with a reproducible choice among numerically tied columns.
!> \param qr Stacked or single-k-point candidate matrix, overwritten by its factorization.
!> \param pivots Original indices of the selected columns.
! **************************************************************************************************
   SUBROUTINE select_qr_pivots(qr, pivots)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: qr(:, :)
      INTEGER, INTENT(OUT)                               :: pivots(:)

      REAL(KIND=dp), PARAMETER                           :: tie_tolerance = 1.0e-12_dp

      COMPLEX(KIND=dp)                                   :: tau(1)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: column(:), work(:)
      INTEGER                                            :: ic, info, ip, iw, nc, nr, original
      INTEGER, ALLOCATABLE                               :: order(:)
      REAL(KIND=dp)                                      :: largest
      REAL(KIND=dp), ALLOCATABLE                         :: norms(:)

      nr = SIZE(qr, 1)
      nc = SIZE(qr, 2)
      CPASSERT(SIZE(pivots) <= MIN(nr, nc))
      ALLOCATE (column(nr), work(MAX(1, nc)), order(nc), norms(nc))
      order = [(ic, ic=1, nc)]
      DO iw = 1, SIZE(pivots)
         DO ic = iw, nc
            norms(ic) = SUM(ABS(qr(iw:nr, ic))**2)
         END DO
         largest = MAXVAL(norms(iw:nc))
         original = HUGE(0)
         ip = iw
         DO ic = iw, nc
            IF (norms(ic) < largest*(1.0_dp - tie_tolerance)) CYCLE
            IF (order(ic) < original) THEN
               ip = ic
               original = order(ic)
            END IF
         END DO
         column = qr(:, iw)
         qr(:, iw) = qr(:, ip)
         qr(:, ip) = column
         order(ip) = order(iw)
         order(iw) = original
         pivots(iw) = original
         ! LAPACK constructs and applies the reflector; only pivot tie-breaking differs.
         CALL zgeqrf(nr - iw + 1, 1, qr(iw, iw), nr, tau, work, SIZE(work), info)
         CPASSERT(info == 0)
         IF (iw < nc) THEN
            CALL zunmqr('L', 'C', nr - iw + 1, nc - iw, 1, qr(iw, iw), nr, tau, &
                        qr(iw, iw + 1), nr, work, SIZE(work), info)
            CPASSERT(info == 0)
         END IF
      END DO
      DEALLOCATE (column, work, order, norms)
   END SUBROUTINE select_qr_pivots

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
