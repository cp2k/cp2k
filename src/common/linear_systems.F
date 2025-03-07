!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Provides interfaces to LAPACK routines for factorisation and
!>      linear system solving
!> \par History
!>      none
!> \author JGH (30-5-2001)
! **************************************************************************************************
MODULE linear_systems

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'linear_systems'

   PUBLIC :: solve_system

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param mysize ...
!> \param eigenvectors ...
! **************************************************************************************************
   SUBROUTINE solve_system(matrix, mysize, eigenvectors)

      REAL(KIND=dp), INTENT(INOUT)                       :: matrix(:, :)
      INTEGER, INTENT(IN)                                :: mysize
      REAL(KIND=dp), INTENT(INOUT)                       :: eigenvectors(:, :)

      INTEGER                                            :: info, lda, ldb, nrhs, ipiv(mysize)

      lda = SIZE(matrix, 1)
      ldb = SIZE(eigenvectors, 1)
      nrhs = SIZE(eigenvectors, 2)

      CALL dgesv(mysize, nrhs, matrix, lda, ipiv, &
                 eigenvectors, ldb, info)
      IF (info /= 0) THEN
         CPABORT("Error in inversion")
      END IF

   END SUBROUTINE solve_system

END MODULE linear_systems

