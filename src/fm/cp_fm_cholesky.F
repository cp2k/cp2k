!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief various cholesky decomposition related routines
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_fm_cholesky
   USE cp_dlaf_utils_api,               ONLY: cp_dlaf_create_grid,&
                                              cp_dlaf_initialize
   USE cp_fm_dlaf_api,                  ONLY: cp_pdpotrf_dlaf,&
                                              cp_pdpotri_dlaf,&
                                              cp_pspotrf_dlaf,&
                                              cp_pspotri_dlaf
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE kinds,                           ONLY: dp,&
                                              sp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_cholesky'

   PUBLIC :: cp_fm_cholesky_decompose, cp_fm_cholesky_invert, &
             cp_fm_cholesky_reduce, cp_fm_cholesky_restore

   ! The following saved variables are Cholesky decomposition global
   ! Stores the default library for Cholesky decomposition
   INTEGER, SAVE, PUBLIC   :: cholesky_type = 0
   ! Minimum matrix size for the use of the DLAF Cholesky decomposition.
   ! ScaLAPACK is used as fallback for all smaller cases.
   INTEGER, SAVE, PUBLIC    :: dlaf_cholesky_n_min = 0
   ! Constants for the diag_type above
   INTEGER, PARAMETER, PUBLIC  :: FM_CHOLESKY_TYPE_SCALAPACK = 101, &
                                  FM_CHOLESKY_TYPE_DLAF = 104
   INTEGER, PARAMETER, PUBLIC :: FM_CHOLESKY_TYPE_DEFAULT = FM_CHOLESKY_TYPE_SCALAPACK

!***
CONTAINS

! **************************************************************************************************
!> \brief used to replace a symmetric positive def. matrix M with its cholesky
!>      decomposition U: M = U^T * U, with U upper triangular
!> \param matrix the matrix to replace with its cholesky decomposition
!> \param n the number of row (and columns) of the matrix &
!>        (defaults to the min(size(matrix)))
!> \param info_out ...
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!> \author Joost
! **************************************************************************************************
   SUBROUTINE cp_fm_cholesky_decompose(matrix, n, info_out)
      TYPE(cp_fm_type), INTENT(IN)             :: matrix
      INTEGER, INTENT(in), OPTIONAL            :: n
      INTEGER, INTENT(out), OPTIONAL           :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_cholesky_decompose'

      INTEGER                                  :: handle, info, my_n
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
#if defined(__parallel)
      INTEGER, DIMENSION(9)                    :: desca
#endif

      CALL timeset(routineN, handle)

      my_n = MIN(matrix%matrix_struct%nrow_global, &
                 matrix%matrix_struct%ncol_global)
      IF (PRESENT(n)) THEN
         CPASSERT(n <= my_n)
         my_n = n
      END IF

      a => matrix%local_data
      a_sp => matrix%local_data_sp

#if defined(__parallel)
      desca(:) = matrix%matrix_struct%descriptor(:)
#if defined(__DLAF)
      IF (cholesky_type == FM_CHOLESKY_TYPE_DLAF .AND. matrix%matrix_struct%nrow_global >= dlaf_cholesky_n_min) THEN
         ! Initialize DLA-Future on-demand; if already initialized, does nothing
         CALL cp_dlaf_initialize()

         ! Create DLAF grid from BLACS context; if already present, does nothing
         CALL cp_dlaf_create_grid(matrix%matrix_struct%context%get_handle())

         IF (matrix%use_sp) THEN
            CALL cp_pspotrf_dlaf('U', my_n, a_sp(:, :), 1, 1, desca, info)
         ELSE
            CALL cp_pdpotrf_dlaf('U', my_n, a(:, :), 1, 1, desca, info)
         END IF
      ELSE
#endif
         IF (matrix%use_sp) THEN
            CALL pspotrf('U', my_n, a_sp(1, 1), 1, 1, desca, info)
         ELSE
            CALL pdpotrf('U', my_n, a(1, 1), 1, 1, desca, info)
         END IF
#if defined(__DLAF)
      END IF
#endif
#else
      IF (matrix%use_sp) THEN
         CALL spotrf('U', my_n, a_sp(1, 1), SIZE(a_sp, 1), info)
      ELSE
         CALL dpotrf('U', my_n, a(1, 1), SIZE(a, 1), info)
      END IF
#endif

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Cholesky decompose failed: the matrix is not positive definite or ill-conditioned.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_cholesky_decompose

! **************************************************************************************************
!> \brief used to replace the cholesky decomposition by the inverse
!> \param matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert (defaults to the min(size(matrix)))
!> \param info_out ...
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_fm_cholesky_invert(matrix, n, info_out)
      TYPE(cp_fm_type), INTENT(IN)              :: matrix
      INTEGER, INTENT(in), OPTIONAL             :: n
      INTEGER, INTENT(OUT), OPTIONAL            :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_cholesky_invert'
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER   :: a_sp
      INTEGER                                   :: info, handle, my_n
#if defined(__parallel)
      INTEGER, DIMENSION(9)                     :: desca
#endif

      CALL timeset(routineN, handle)

      my_n = MIN(matrix%matrix_struct%nrow_global, &
                 matrix%matrix_struct%ncol_global)
      IF (PRESENT(n)) THEN
         CPASSERT(n <= my_n)
         my_n = n
      END IF

      a => matrix%local_data
      a_sp => matrix%local_data_sp

#if defined(__parallel)

      desca(:) = matrix%matrix_struct%descriptor(:)

#if defined(__DLAF)
      IF (cholesky_type == FM_CHOLESKY_TYPE_DLAF .AND. matrix%matrix_struct%nrow_global >= dlaf_cholesky_n_min) THEN
         ! Initialize DLA-Future on-demand; if already initialized, does nothing
         CALL cp_dlaf_initialize()

         ! Create DLAF grid from BLACS context; if already present, does nothing
         CALL cp_dlaf_create_grid(matrix%matrix_struct%context%get_handle())

         IF (matrix%use_sp) THEN
            CALL cp_pspotri_dlaf('U', my_n, a_sp(:, :), 1, 1, desca, info)
         ELSE
            CALL cp_pdpotri_dlaf('U', my_n, a(:, :), 1, 1, desca, info)
         END IF
      ELSE
#endif
         IF (matrix%use_sp) THEN
            CALL pspotri('U', my_n, a_sp(1, 1), 1, 1, desca, info)
         ELSE
            CALL pdpotri('U', my_n, a(1, 1), 1, 1, desca, info)
         END IF
#if defined(__DLAF)
      END IF
#endif

#else

      IF (matrix%use_sp) THEN
         CALL spotri('U', my_n, a_sp(1, 1), SIZE(a_sp, 1), info)
      ELSE
         CALL dpotri('U', my_n, a(1, 1), SIZE(a, 1), info)
      END IF

#endif

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CPABORT("Cholesky invert failed: the matrix is not positive definite or ill-conditioned.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_cholesky_invert

! **************************************************************************************************
!> \brief reduce a matrix pencil A,B to normal form
!>      B has to be cholesky decomposed with  cp_fm_cholesky_decompose
!>      before calling this routine
!>      A,B -> inv(U^T)*A*inv(U),1
!>      (AX=BX -> inv(U^T)*A*inv(U)*U*X=U*X hence evecs U*X)
!> \param matrix the symmetric matrix A
!> \param matrixb the cholesky decomposition of matrix B
!> \param itype ...
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_fm_cholesky_reduce(matrix, matrixb, itype)
      TYPE(cp_fm_type), INTENT(IN)              :: matrix, matrixb
      INTEGER, OPTIONAL                         :: itype

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_cholesky_reduce'
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a, b
      INTEGER                                   :: info, handle
      INTEGER                                   :: n, my_itype
#if defined(__parallel)
      REAL(KIND=dp)                             :: scale
      INTEGER, DIMENSION(9)                     :: desca, descb
#endif

      CALL timeset(routineN, handle)

      n = matrix%matrix_struct%nrow_global

      my_itype = 1
      IF (PRESENT(itype)) my_itype = itype

      a => matrix%local_data
      b => matrixb%local_data

#if defined(__parallel)

      desca(:) = matrix%matrix_struct%descriptor(:)
      descb(:) = matrixb%matrix_struct%descriptor(:)

      CALL pdsygst(my_itype, 'U', n, a(1, 1), 1, 1, desca, b(1, 1), 1, 1, descb, scale, info)

      ! this is supposed to be one in current version of lapack
      ! if not, eigenvalues have to be scaled by this number
      IF (scale /= 1.0_dp) &
         CPABORT("scale not equal 1 (scale="//cp_to_string(scale)//")")
#else

      CALL dsygst(my_itype, 'U', n, a(1, 1), n, b(1, 1), n, info)

#endif

      CPASSERT(info == 0)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_cholesky_reduce

! **************************************************************************************************
!> \brief apply Cholesky decomposition
!>        op can be "SOLVE" (out = U^-1 * in) or "MULTIPLY"  (out = U * in)
!>        pos can be "LEFT" or "RIGHT" (U at the left or at the right)
!> \param fm_matrix ...
!> \param neig ...
!> \param fm_matrixb ...
!> \param fm_matrixout ...
!> \param op ...
!> \param pos ...
!> \param transa ...
! **************************************************************************************************
   SUBROUTINE cp_fm_cholesky_restore(fm_matrix, neig, fm_matrixb, fm_matrixout, op, pos, transa)
      TYPE(cp_fm_type), INTENT(IN)              :: fm_matrix, fm_matrixb, fm_matrixout
      INTEGER, INTENT(IN)                       :: neig
      CHARACTER(LEN=*), INTENT(IN)              :: op
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL    :: pos
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL    :: transa

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_cholesky_restore'
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a, b, outm
      REAL(KIND=sp), DIMENSION(:, :), POINTER   :: a_sp, b_sp, outm_sp
      REAL(KIND=dp)                             :: alpha
      INTEGER                                   :: handle, mdim, n
      CHARACTER                                 :: chol_pos, chol_transa
#if defined(__parallel)
      INTEGER                                   :: i
      INTEGER, DIMENSION(9)                     :: desca, descb, descout
#elif defined(__MKL) && (2025 <= __INTEL_MKL__)
      REAL(KIND=dp)                             :: beta
#endif

      CALL timeset(routineN, handle)

      ! wrong argument op
      CPASSERT((op == "SOLVE") .OR. (op == "MULTIPLY"))
      ! not the same precision
      CPASSERT(fm_matrix%use_sp .EQV. fm_matrixb%use_sp)
      CPASSERT(fm_matrix%use_sp .EQV. fm_matrixout%use_sp)

      chol_pos = 'L'
      IF (PRESENT(pos)) chol_pos = pos(1:1)
      CPASSERT((chol_pos == 'L') .OR. (chol_pos == 'R'))

      chol_transa = 'N'
      IF (PRESENT(transa)) chol_transa = transa

      ! notice b is the cholesky guy
      a => fm_matrix%local_data
      b => fm_matrixb%local_data
      outm => fm_matrixout%local_data
      a_sp => fm_matrix%local_data_sp
      b_sp => fm_matrixb%local_data_sp
      outm_sp => fm_matrixout%local_data_sp
      n = fm_matrix%matrix_struct%nrow_global
      mdim = MERGE(1, 2, 'L' == chol_pos)
      alpha = 1.0_dp

#if defined(__parallel)
      desca(:) = fm_matrix%matrix_struct%descriptor(:)
      descb(:) = fm_matrixb%matrix_struct%descriptor(:)
      descout(:) = fm_matrixout%matrix_struct%descriptor(:)
      DO i = 1, neig
         IF (fm_matrix%use_sp) THEN
            CALL pscopy(n, a_sp(1, 1), 1, i, desca, 1, outm_sp(1, 1), 1, i, descout, 1)
         ELSE
            CALL pdcopy(n, a(1, 1), 1, i, desca, 1, outm(1, 1), 1, i, descout, 1)
         END IF
      END DO
      IF (op == "SOLVE") THEN
         IF (fm_matrix%use_sp) THEN
            CALL pstrsm(chol_pos, 'U', chol_transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), 1, 1, descb, &
                        outm_sp(1, 1), 1, 1, descout)
         ELSE
            CALL pdtrsm(chol_pos, 'U', chol_transa, 'N', n, neig, alpha, b(1, 1), 1, 1, descb, outm(1, 1), 1, 1, descout)
         END IF
      ELSE
         IF (fm_matrix%use_sp) THEN
            CALL pstrmm(chol_pos, 'U', chol_transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), 1, 1, descb, &
                        outm_sp(1, 1), 1, 1, descout)
         ELSE
            CALL pdtrmm(chol_pos, 'U', chol_transa, 'N', n, neig, alpha, b(1, 1), 1, 1, descb, outm(1, 1), 1, 1, descout)
         END IF
      END IF
#elif defined(__MKL) && (2025 <= __INTEL_MKL__)
      beta = 0.0_dp
      IF (op == "SOLVE") THEN
         IF (fm_matrix%use_sp) THEN
            CALL strsm_oop(chol_pos, 'U', chol_transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), SIZE(b_sp, mdim), a_sp(1, 1), n, REAL(beta, sp), outm_sp(1, 1), n)
         ELSE
            CALL dtrsm_oop(chol_pos, 'U', chol_transa, 'N', n, neig, alpha, b(1, 1), SIZE(b, mdim), a(1, 1), n, beta, outm(1, 1), n)
         END IF
      ELSE
         IF (fm_matrix%use_sp) THEN
            CALL strmm_oop(chol_pos, 'U', chol_transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), SIZE(b_sp, mdim), a_sp(1, 1), n, REAL(beta, sp), outm_sp(1, 1), n)
         ELSE
            CALL dtrmm_oop(chol_pos, 'U', chol_transa, 'N', n, neig, alpha, b(1, 1), SIZE(b, mdim), a(1, 1), n, beta, outm(1, 1), n)
         END IF
      END IF
#else
      IF (fm_matrix%use_sp) THEN
         CALL scopy(neig*n, a_sp(1, 1), 1, outm_sp(1, 1), 1)
      ELSE
         CALL dcopy(neig*n, a(1, 1), 1, outm(1, 1), 1)
      END IF
      IF (op == "SOLVE") THEN
         IF (fm_matrix%use_sp) THEN
            CALL strsm(chol_pos, 'U', chol_transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), SIZE(b_sp, mdim), outm_sp(1, 1), n)
         ELSE
            CALL dtrsm(chol_pos, 'U', chol_transa, 'N', n, neig, alpha, b(1, 1), SIZE(b, mdim), outm(1, 1), n)
         END IF
      ELSE
         IF (fm_matrix%use_sp) THEN
            CALL strmm(chol_pos, 'U', chol_transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), SIZE(b_sp, mdim), outm_sp(1, 1), n)
         ELSE
            CALL dtrmm(chol_pos, 'U', chol_transa, 'N', n, neig, alpha, b(1, 1), SIZE(b, mdim), outm(1, 1), n)
         END IF
      END IF
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_cholesky_restore

END MODULE cp_fm_cholesky
