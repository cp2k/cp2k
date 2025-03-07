!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief provides a unified interface to lapack geev routines
!> \par History
!>       2014.09 created [Florian Schiffmann]
!>       2023.12 Removed support for single-precision [Ole Schuett]
!>       2024.12 Removed support for complex input matrices [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
MODULE arnoldi_geev
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'arnoldi_geev'

   PUBLIC :: arnoldi_general_local_diag, arnoldi_tridiag_local_diag, arnoldi_symm_local_diag

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
! **************************************************************************************************
   SUBROUTINE arnoldi_symm_local_diag(jobvr, matrix, ndim, evals, revec)
      CHARACTER(1)                                       :: jobvr
      REAL(dp), DIMENSION(:, :)                          :: matrix
      INTEGER                                            :: ndim
      COMPLEX(dp), DIMENSION(:)                          :: evals
      COMPLEX(dp), DIMENSION(:, :)                       :: revec

      INTEGER                                            :: i, info, liwork, lwork, iwork(3 + 5*ndim)
      REAL(dp)                                           :: tmp_array(ndim, ndim), &
                                                            work(1 + 6*ndim + 2*ndim**2)
      REAL(dp), DIMENSION(ndim)                          :: eval

      lwork = 1 + 6*ndim + 2*ndim**2
      liwork = 3 + 5*ndim

      tmp_array(:, :) = matrix(:, :)
      CALL dsyevd(jobvr, "U", ndim, tmp_array, ndim, eval, work, lwork, iwork, liwork, info)

      DO i = 1, ndim
         revec(:, i) = CMPLX(tmp_array(:, i), REAL(0.0, dp), dp)
         evals(i) = CMPLX(eval(i), 0.0, dp)
      END DO

   END SUBROUTINE arnoldi_symm_local_diag

! **************************************************************************************************
!> \brief ...
!> \param jobvl ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
!> \param levec ...
! **************************************************************************************************
   SUBROUTINE arnoldi_tridiag_local_diag(jobvl, jobvr, matrix, ndim, evals, revec, levec)
      CHARACTER(1)                                       :: jobvl, jobvr
      REAL(dp), DIMENSION(:, :)                          :: matrix
      INTEGER                                            :: ndim
      COMPLEX(dp), DIMENSION(:)                          :: evals
      COMPLEX(dp), DIMENSION(:, :)                       :: revec, levec

      INTEGER                                            :: i, info
      REAL(dp)                                           :: work(20*ndim)
      REAL(dp), DIMENSION(ndim)                          :: diag, offdiag
      REAL(dp), DIMENSION(ndim, ndim)                    :: evec_r

      MARK_USED(jobvl) !the argument has to be here for the template to work

      levec(1, 1) = CMPLX(0.0, 0.0, dp)
      info = 0
      diag(ndim) = matrix(ndim, ndim)
      DO i = 1, ndim - 1
         diag(i) = matrix(i, i)
         offdiag(i) = matrix(i + 1, i)

      END DO

      CALL dstev(jobvr, ndim, diag, offdiag, evec_r, ndim, work, info)

      DO i = 1, ndim
         revec(:, i) = CMPLX(evec_r(:, i), REAL(0.0, dp), dp)
         evals(i) = CMPLX(diag(i), 0.0, dp)
      END DO
   END SUBROUTINE arnoldi_tridiag_local_diag

! **************************************************************************************************
!> \brief ...
!> \param jobvl ...
!> \param jobvr ...
!> \param matrix ...
!> \param ndim ...
!> \param evals ...
!> \param revec ...
!> \param levec ...
! **************************************************************************************************
   SUBROUTINE arnoldi_general_local_diag(jobvl, jobvr, matrix, ndim, evals, revec, levec)
      CHARACTER(1)                                       :: jobvl, jobvr
      REAL(dp), DIMENSION(:, :)                          :: matrix
      INTEGER                                            :: ndim
      COMPLEX(dp), DIMENSION(:)                          :: evals
      COMPLEX(dp), DIMENSION(:, :)                       :: revec, levec

      INTEGER                                            :: i, info, lwork
      LOGICAL                                            :: selects(ndim)
      REAL(dp)                                           :: norm, tmp_array(ndim, ndim), &
                                                            work(20*ndim)
      REAL(dp), DIMENSION(ndim)                          :: eval1, eval2
      REAL(dp), DIMENSION(ndim, ndim)                    :: evec_l, evec_r

      MARK_USED(jobvr) !the argument has to be here for the template to work
      MARK_USED(jobvl) !the argument has to be here for the template to work

      eval1 = REAL(0.0, dp); eval2 = REAL(0.0, dp)
      tmp_array(:, :) = matrix(:, :)
      ! ask lapack how much space it would like in the work vector, don't ask me why
      lwork = -1
      CALL dhseqr('S', 'I', ndim, 1, ndim, tmp_array, ndim, eval1, eval2, evec_r, ndim, work, lwork, info)

      lwork = MIN(20*ndim, INT(work(1)))
      CALL dhseqr('S', 'I', ndim, 1, ndim, tmp_array, ndim, eval1, eval2, evec_r, ndim, work, lwork, info)
      CALL dtrevc('R', 'B', selects, ndim, tmp_array, ndim, evec_l, ndim, evec_r, ndim, ndim, ndim, work, info)

      ! compose the eigenvectors, lapacks way of storing them is a pain
      ! if eval is complex, then the complex conj pair of evec can be constructed from the i and i+1st evec
      ! Unfortunately dtrevc computes the ev such that the largest is set to one and not normalized
      i = 1
      DO WHILE (i .LE. ndim)
         IF (ABS(eval2(i)) .LT. EPSILON(REAL(0.0, dp))) THEN
            evec_r(:, i) = evec_r(:, i)/SQRT(DOT_PRODUCT(evec_r(:, i), evec_r(:, i)))
            revec(:, i) = CMPLX(evec_r(:, i), REAL(0.0, dp), dp)
            levec(:, i) = CMPLX(evec_l(:, i), REAL(0.0, dp), dp)
            i = i + 1
         ELSE IF (eval2(i) .GT. EPSILON(REAL(0.0, dp))) THEN
            norm = SQRT(SUM(evec_r(:, i)**2.0_dp) + SUM(evec_r(:, i + 1)**2.0_dp))
            revec(:, i) = CMPLX(evec_r(:, i), evec_r(:, i + 1), dp)/norm
            revec(:, i + 1) = CMPLX(evec_r(:, i), -evec_r(:, i + 1), dp)/norm
            levec(:, i) = CMPLX(evec_l(:, i), evec_l(:, i + 1), dp)
            levec(:, i + 1) = CMPLX(evec_l(:, i), -evec_l(:, i + 1), dp)
            i = i + 2
         ELSE
            CPABORT('something went wrong while sorting the EV in arnoldi_geev')
         END IF
      END DO

      ! this is to keep the interface consistent with complex geev
      DO i = 1, ndim
         evals(i) = CMPLX(eval1(i), eval2(i), dp)
      END DO

   END SUBROUTINE arnoldi_general_local_diag

END MODULE arnoldi_geev
