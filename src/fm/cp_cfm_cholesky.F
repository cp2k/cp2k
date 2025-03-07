!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief various cholesky decomposition related routines
!> \par History
!>      12.2002 Moved routines from cp_cfm_basic_linalg to this new module [Rocco Meli]
! **************************************************************************************************
MODULE cp_cfm_cholesky
   USE cp_cfm_types, ONLY: cp_cfm_type
   USE kinds, ONLY: dp

#if defined(__DLAF)
   USE cp_cfm_dlaf_api, ONLY: cp_cfm_pzpotrf_dlaf
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_cholesky'

   PUBLIC :: cp_cfm_cholesky_decompose, &
             cp_cfm_cholesky_invert

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Used to replace a symmetric positive definite matrix M with its Cholesky
!>      decomposition U: M = U^T * U, with U upper triangular.
!> \param matrix   the matrix to replace with its Cholesky decomposition
!> \param n        the number of row (and columns) of the matrix &
!>                 (defaults to the min(size(matrix)))
!> \param info_out if present, outputs info from (p)zpotrf
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!>      09.2021 removed CPASSERT(info == 0) since there is already check of info [Jan Wilhelm]
!>      12.2024 Added DLA-Future support [Rocco Meli]
!> \author Joost
! **************************************************************************************************
   SUBROUTINE cp_cfm_cholesky_decompose(matrix, n, info_out)
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix
      INTEGER, INTENT(in), OPTIONAL                      :: n
      INTEGER, INTENT(out), OPTIONAL                     :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_cholesky_decompose'

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: handle, info, my_n
#if defined(__parallel)
      INTEGER, DIMENSION(9)                              :: desca
#else
      INTEGER                                            :: lda
#endif

      CALL timeset(routineN, handle)

      my_n = MIN(matrix%matrix_struct%nrow_global, &
                 matrix%matrix_struct%ncol_global)
      IF (PRESENT(n)) THEN
         CPASSERT(n <= my_n)
         my_n = n
      END IF

      a => matrix%local_data

#if defined(__parallel)
      desca(:) = matrix%matrix_struct%descriptor(:)

#if defined(__DLAF)
      CALL cp_cfm_pzpotrf_dlaf('U', my_n, a, 1, 1, desca, info)
#else
      CALL pzpotrf('U', my_n, a(1, 1), 1, 1, desca, info)
#endif
#else
      lda = SIZE(a, 1)
      CALL zpotrf('U', my_n, a(1, 1), lda, info)
#endif

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CALL cp_abort(__LOCATION__, &
                          "Cholesky decompose failed: matrix is not positive definite or ill-conditioned")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_cholesky_decompose

! **************************************************************************************************
!> \brief Used to replace Cholesky decomposition by the inverse.
!> \param matrix : the matrix to invert (must be an upper triangular matrix),
!>                 and is the output of Cholesky decomposition
!> \param n : size of the matrix to invert (defaults to the min(size(matrix)))
!> \param info_out : if present, outputs info of (p)zpotri
!> \par History
!>      05.2002 created Lianheng Tong, based on cp_fm_cholesky_invert
!> \author Lianheng Tong
! **************************************************************************************************
   SUBROUTINE cp_cfm_cholesky_invert(matrix, n, info_out)
      TYPE(cp_cfm_type), INTENT(IN)              :: matrix
      INTEGER, INTENT(in), OPTIONAL              :: n
      INTEGER, INTENT(out), OPTIONAL             :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_cholesky_invert'
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER :: aa
      INTEGER                                    :: info, handle
      INTEGER                                    :: my_n
#if defined(__parallel)
      INTEGER, DIMENSION(9)                      :: desca
#endif

      CALL timeset(routineN, handle)

      my_n = MIN(matrix%matrix_struct%nrow_global, &
                 matrix%matrix_struct%ncol_global)
      IF (PRESENT(n)) THEN
         CPASSERT(n <= my_n)
         my_n = n
      END IF

      aa => matrix%local_data

#if defined(__parallel)
      desca = matrix%matrix_struct%descriptor
      CALL pzpotri('U', my_n, aa(1, 1), 1, 1, desca, info)
#else
      CALL zpotri('U', my_n, aa(1, 1), SIZE(aa, 1), info)
#endif

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CALL cp_abort(__LOCATION__, &
                          "Cholesky invert failed: the matrix is not positive definite or ill-conditioned.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_cholesky_invert

END MODULE cp_cfm_cholesky
