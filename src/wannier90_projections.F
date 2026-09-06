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
   PUBLIC :: wannier90_select_projections

CONTAINS

! **************************************************************************************************
!> \brief Select the same trial orbitals at every k-point using a mesh-wide pivoted QR.
!> \param candidates Band/candidate/k-point projections, stored only on the source rank.
!> \param para_env Global parallel environment.
!> \param amn Selected projection matrices, allocated with the target dimensions on all ranks.
!> \param pivots Selected candidate indices, allocated on all ranks.
!> \param min_svalue Smallest singular value of the selected projections over the mesh.
! **************************************************************************************************
   SUBROUTINE wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
      COMPLEX(KIND=dp), INTENT(IN)                       :: candidates(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), INTENT(OUT)                      :: amn(:, :, :)
      INTEGER, INTENT(OUT)                               :: pivots(:)
      REAL(KIND=dp), INTENT(OUT)                         :: min_svalue

      COMPLEX(KIND=dp)                                   :: query(1), unused_u(1, 1), unused_v(1, 1)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: matrix(:, :), qr(:, :), tau(:), work(:)
      INTEGER                                            :: ik, info, lwork, nb, nc, nk, nr, nw
      INTEGER, ALLOCATABLE                               :: order(:)
      REAL(KIND=dp), ALLOCATABLE                         :: rwork(:), svalues(:)

      nb = SIZE(amn, 1)
      nw = SIZE(amn, 2)
      nk = SIZE(amn, 3)
      CPASSERT(nb >= nw .AND. nw > 0 .AND. SIZE(pivots) == nw)
      min_svalue = 0.0_dp
      IF (para_env%is_source()) THEN
         nc = SIZE(candidates, 2)
         CPASSERT(SIZE(candidates, 1) == nb .AND. SIZE(candidates, 3) == nk .AND. nc >= nw)
         nr = nb*nk
         ALLOCATE (qr(nr, nc), tau(MIN(nr, nc)), order(nc), rwork(2*nc))
         DO ik = 1, nk
            qr((ik - 1)*nb + 1:ik*nb, :) = candidates(:, :, ik)
         END DO
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

         ! Global column independence alone does not guarantee full rank at every k-point.
         ALLOCATE (matrix(nb, nw), svalues(nw), rwork(5*nw))
         matrix(:, :) = amn(:, :, 1)
         CALL zgesvd('N', 'N', nb, nw, matrix, nb, svalues, unused_u, 1, unused_v, 1, &
                     query, -1, rwork, info)
         CPASSERT(info == 0)
         lwork = MAX(1, INT(REAL(query(1), KIND=dp)))
         ALLOCATE (work(lwork))
         min_svalue = HUGE(0.0_dp)
         DO ik = 1, nk
            matrix(:, :) = amn(:, :, ik)
            CALL zgesvd('N', 'N', nb, nw, matrix, nb, svalues, unused_u, 1, unused_v, 1, &
                        work, lwork, rwork, info)
            CPASSERT(info == 0)
            min_svalue = MIN(min_svalue, MINVAL(svalues))
         END DO
         DEALLOCATE (matrix, svalues, rwork, work)
      END IF
      CALL para_env%bcast(amn)
      CALL para_env%bcast(pivots)
      CALL para_env%bcast(min_svalue)

   END SUBROUTINE wannier90_select_projections

END MODULE wannier90_projections
