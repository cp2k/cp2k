!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common full-matrix operations used by GW modules.
!> \par History
!>      09.2026 created Jan Wilhelm
! **************************************************************************************************
MODULE gw_utils_fm
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_uplo_to_full
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE mathlib,                         ONLY: diag_complex,&
                                              gemm_square
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils_fm'

   PUBLIC :: cfm_contract_ABA, cfm_Wc_from_chi, cfm_Wc_from_cholesky_eps, &
             fm_contract_ABA, fm_invert, fm_sqrt, local_complex_power, &
             local_Wc_from_chi

   INTERFACE fm_contract_ABA
      MODULE PROCEDURE fm_contract_ABA, fm_contract_ABA_set
   END INTERFACE fm_contract_ABA

CONTAINS

! **************************************************************************************************
!> \brief Computes A^H B A for complex full matrices.
!> \param matrix_A left and right matrix A; H denotes its Hermitian transpose
!> \param matrix_B input matrix B; overwritten by A^H B A if matrix_C is absent
!> \param matrix_C optional output matrix C=A^H B A
! **************************************************************************************************
   SUBROUTINE cfm_contract_ABA(matrix_A, matrix_B, matrix_C)
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix_A
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrix_B
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: matrix_C

      INTEGER                                            :: ncol, nrow
      TYPE(cp_cfm_type)                                  :: work

      CPASSERT(ASSOCIATED(matrix_A%matrix_struct))
      CPASSERT(ASSOCIATED(matrix_B%matrix_struct))
      CALL cp_cfm_get_info(matrix_A, nrow_global=nrow, ncol_global=ncol)
      CPASSERT(nrow == ncol)
      CALL cp_cfm_create(work, matrix_A%matrix_struct)

      CALL parallel_gemm("N", "N", nrow, nrow, nrow, z_one, matrix_B, matrix_A, &
                         z_zero, work)
      IF (PRESENT(matrix_C)) THEN
         CPASSERT(ASSOCIATED(matrix_C%matrix_struct))
         CALL parallel_gemm("C", "N", nrow, nrow, nrow, z_one, matrix_A, work, &
                            z_zero, matrix_C)
      ELSE
         CALL parallel_gemm("C", "N", nrow, nrow, nrow, z_one, matrix_A, work, &
                            z_zero, matrix_B)
      END IF

      CALL cp_cfm_release(work)

   END SUBROUTINE cfm_contract_ABA

! **************************************************************************************************
!> \brief Compute the correlation part of the screened interaction from chi.
!>        A=M^-1 B^H and V=B^H B. On return, chi_Wc contains
!>        B^H[(I+A^H chi A)^-1-I]B.
!> \param chi_Wc input chi(k,iw); overwritten by Wc(k,iw)
!> \param A metric-weighted Coulomb factor M^-1 B^H
!> \param B Coulomb factor satisfying V=B^H B
! **************************************************************************************************
   SUBROUTINE cfm_Wc_from_chi(chi_Wc, A, B)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: chi_Wc
      TYPE(cp_cfm_type), INTENT(IN)                      :: A, B

      INTEGER                                            :: info

      CALL cfm_contract_ABA(A, chi_Wc)
      CALL cfm_add_on_diag(chi_Wc, z_one)
      CALL cp_cfm_cholesky_decompose(chi_Wc, info_out=info)
      CPASSERT(info == 0)
      CALL cfm_Wc_from_cholesky_eps(chi_Wc, B)

   END SUBROUTINE cfm_Wc_from_chi

! **************************************************************************************************
!> \brief Compute the correlation part of the screened interaction for local complex matrices.
!>        The input chi_Wc must already have its negative eigenvalues removed. On return,
!>        chi_Wc contains Wc = B (epsilon^-1 - I) B^H, where epsilon = I + A^H chi A.
!> \param chi_Wc input chi; overwritten by Wc
!> \param A metric-weighted Coulomb factor M^-1 B^H
!> \param B Coulomb factor satisfying V=B^H B
!> \param work local complex workspace for epsilon
! **************************************************************************************************
   SUBROUTINE local_Wc_from_chi(chi_Wc, A, B, work)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: chi_Wc
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: A, B
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: work

      INTEGER                                            :: n

      n = SIZE(chi_Wc, 1)
      CPASSERT(n == SIZE(chi_Wc, 2))
      CPASSERT(ALL(SHAPE(A) == SHAPE(chi_Wc)))
      CPASSERT(ALL(SHAPE(B) == SHAPE(chi_Wc)))
      CPASSERT(ALL(SHAPE(work) == SHAPE(chi_Wc)))

      CALL gemm_square(A, 'C', chi_Wc, 'N', A, 'N', work)
      CALL local_add_on_diag(work, z_one)
      CALL local_complex_power(work, -1.0_dp, 0.0_dp)
      CALL local_add_on_diag(work, -z_one)
      CALL gemm_square(B, 'N', work, 'N', B, 'C', chi_Wc)

   END SUBROUTINE local_Wc_from_chi

! **************************************************************************************************
!> \brief Complete Wc after a caller has formed and Cholesky-factorized I+A^H chi A.
!>        The separate entry point allows RPA to evaluate its log determinant first.
!> \param eps_chol input Cholesky factor; overwritten by eps^-1-I or Wc
!> \param B Coulomb factor satisfying V=B^H B
!> \param Wc optional separate output; if absent, eps_chol is overwritten by Wc
!> \param n_active optional number of physical rows/columns
! **************************************************************************************************
   SUBROUTINE cfm_Wc_from_cholesky_eps(eps_chol, B, Wc, n_active)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: eps_chol
      TYPE(cp_cfm_type), INTENT(IN)                      :: B
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: Wc
      INTEGER, INTENT(IN), OPTIONAL                      :: n_active

      CALL cp_cfm_cholesky_invert(eps_chol)
      CALL cp_cfm_uplo_to_full(eps_chol)
      CALL cfm_add_on_diag(eps_chol, -z_one, n_active)
      CALL cfm_contract_ABA(B, eps_chol, Wc)

   END SUBROUTINE cfm_Wc_from_cholesky_eps

! **************************************************************************************************
!> \brief Add a scalar to the diagonal of a distributed complex full matrix.
!> \param matrix input/output matrix
!> \param alpha scalar added to each diagonal element
!> \param n_active optional number of physical rows/columns
! **************************************************************************************************
   SUBROUTINE cfm_add_on_diag(matrix, alpha, n_active)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrix
      COMPLEX(KIND=dp), INTENT(IN)                       :: alpha
      INTEGER, INTENT(IN), OPTIONAL                      :: n_active

      INTEGER                                            :: i_row, j_col, max_index, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      max_index = HUGE(max_index)
      IF (PRESENT(n_active)) max_index = n_active

      CALL cp_cfm_get_info(matrix=matrix, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      DO j_col = 1, ncol_local
         DO i_row = 1, nrow_local
            IF (row_indices(i_row) == col_indices(j_col) .AND. row_indices(i_row) <= max_index) THEN
               matrix%local_data(i_row, j_col) = matrix%local_data(i_row, j_col) + alpha
            END IF
         END DO
      END DO

   END SUBROUTINE cfm_add_on_diag

! **************************************************************************************************
!> \brief Add a scalar to the diagonal of a local complex matrix.
!> \param matrix input/output matrix
!> \param alpha scalar added to each diagonal element
! **************************************************************************************************
   SUBROUTINE local_add_on_diag(matrix, alpha)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: matrix
      COMPLEX(KIND=dp), INTENT(IN)                       :: alpha

      INTEGER                                            :: i, n

      n = SIZE(matrix, 1)
      CPASSERT(n == SIZE(matrix, 2))

      DO i = 1, n
         matrix(i, i) = matrix(i, i) + alpha
      END DO

   END SUBROUTINE local_add_on_diag

! **************************************************************************************************
!> \brief Compute a spectral power of a local complex Hermitian matrix.
!>        Eigenvalues not larger than eps are discarded.
!> \param matrix input/output matrix
!> \param exponent matrix exponent
!> \param eps eigenvalue threshold
!> \param cond_nr optional condition number
!> \param min_ev optional minimum absolute eigenvalue
!> \param max_ev optional maximum absolute eigenvalue
! **************************************************************************************************
   SUBROUTINE local_complex_power(matrix, exponent, eps, cond_nr, min_ev, max_ev)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: exponent, eps
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: cond_nr, min_ev, max_ev

      CHARACTER(len=*), PARAMETER :: routineN = 'local_complex_power'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenvectors
      INTEGER                                            :: handle, i, n
      REAL(KIND=dp)                                      :: pos_eval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(n == SIZE(matrix, 2))
      matrix(:, :) = 0.5_dp*(matrix(:, :) + CONJG(TRANSPOSE(matrix(:, :))))

      ALLOCATE (eigenvalues(n), eigenvectors(n, n))
      CALL diag_complex(matrix, eigenvectors, eigenvalues)

      IF (PRESENT(cond_nr)) cond_nr = MAXVAL(ABS(eigenvalues))/MINVAL(ABS(eigenvalues))
      IF (PRESENT(min_ev)) min_ev = MINVAL(ABS(eigenvalues))
      IF (PRESENT(max_ev)) max_ev = MAXVAL(ABS(eigenvalues))

      DO i = 1, n
         IF (eps < eigenvalues(i)) THEN
            pos_eval = eigenvalues(i)**(0.5_dp*exponent)
         ELSE
            pos_eval = 0.0_dp
         END IF
         eigenvectors(:, i) = eigenvectors(:, i)*pos_eval
      END DO

      CALL gemm_square(eigenvectors, 'N', eigenvectors, 'C', matrix)

      DEALLOCATE (eigenvalues, eigenvectors)

      CALL timestop(handle)

   END SUBROUTINE local_complex_power

! **************************************************************************************************
!> \brief Inverts a symmetric matrix. First, Cholesky decomposition is tried.
!>        If it fails, the matrix is diagonalized and inverted by taking
!>        inverse eigenvalues.
!> \param matrix_A input A; output inverse or filtered pseudoinverse
!> \param eigenvalue_threshold optional relative eigenvalue threshold
!> \param unit_nr optional output unit for a failed Cholesky decomposition
! **************************************************************************************************
   SUBROUTINE fm_invert(matrix_A, eigenvalue_threshold, unit_nr)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix_A
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eigenvalue_threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      INTEGER                                            :: info, n_dependent
      LOGICAL                                            :: cholesky_successful
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_fm_type)                                   :: work

      threshold = 0.0_dp
      IF (PRESENT(eigenvalue_threshold)) threshold = eigenvalue_threshold
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      ! Preserve A because a failed Cholesky decomposition overwrites it.
      CALL cp_fm_to_fm(matrix_A, work)
      CALL cp_fm_cholesky_decompose(matrix_A, info_out=info)
      cholesky_successful = info == 0
      n_dependent = 0
      IF (cholesky_successful) THEN
         CALL cp_fm_cholesky_invert(matrix_A)
         CALL cp_fm_uplo_to_full(matrix_A, work)
      ELSE
         CALL cp_fm_to_fm(work, matrix_A)
         CALL cp_fm_power(matrix_A, work, -1.0_dp, threshold, n_dependent)
         IF (PRESENT(unit_nr)) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') &
                  'Cholesky decomposition failed; matrix inverted by diagonalization.'
               WRITE (unit_nr, '(T2,A,T72,I9)') 'Discarded eigenmodes:', n_dependent
            END IF
         END IF
      END IF
      CALL cp_fm_release(work)

   END SUBROUTINE fm_invert

! **************************************************************************************************
!> \brief For input A, computes B such that B^T B=A.
!>        First, Cholesky decomposition is tried. If it fails, B is computed
!>        by diagonalizing A.
!> \param matrix_A symmetric input matrix A, retained unchanged
!> \param matrix_B output matrix B
!> \param eigenvalue_threshold optional relative eigenvalue threshold
!> \param unit_nr optional output unit for a failed Cholesky decomposition
! **************************************************************************************************
   SUBROUTINE fm_sqrt(matrix_A, matrix_B, eigenvalue_threshold, unit_nr)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_A
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix_B
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eigenvalue_threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      INTEGER                                            :: i_row, info, j_col, n_dependent, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: cholesky_successful
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_fm_type)                                   :: work

      threshold = 0.0_dp
      IF (PRESENT(eigenvalue_threshold)) threshold = eigenvalue_threshold
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      CALL cp_fm_to_fm(matrix_A, matrix_B)
      CALL cp_fm_cholesky_decompose(matrix_B, info_out=info)
      cholesky_successful = info == 0
      IF (cholesky_successful) THEN
         n_dependent = 0
         CALL cp_fm_get_info(matrix=matrix_B, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices)
         DO j_col = 1, ncol_local
            DO i_row = 1, nrow_local
               IF (row_indices(i_row) > col_indices(j_col)) THEN
                  matrix_B%local_data(i_row, j_col) = 0.0_dp
               END IF
            END DO
         END DO
      ELSE
         CALL cp_fm_to_fm(matrix_A, matrix_B)
         CALL cp_fm_power(matrix_B, work, 0.5_dp, threshold, n_dependent)
         IF (PRESENT(unit_nr)) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') &
                  'Cholesky decomposition failed, using diagonalization.'
               WRITE (unit_nr, '(T2,A,T72,I9)') 'Discarded eigenmodes:', n_dependent
            END IF
         END IF
      END IF
      CALL cp_fm_release(work)

   END SUBROUTINE fm_sqrt

! **************************************************************************************************
!> \brief Computes A^T B A.
!> \param matrix_A left and right matrix A
!> \param matrix_B input matrix B; overwritten by A^T B A if matrix_C is absent
!> \param matrix_C optional output matrix C=A^T B A
! **************************************************************************************************
   SUBROUTINE fm_contract_ABA(matrix_A, matrix_B, matrix_C)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_A
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix_B
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL          :: matrix_C

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_type)                                   :: work

      CPASSERT(ASSOCIATED(matrix_A%matrix_struct))
      CPASSERT(ASSOCIATED(matrix_B%matrix_struct))
      CALL cp_fm_get_info(matrix_A, nrow_global=nrow, ncol_global=ncol)
      CPASSERT(nrow == ncol)
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      CALL parallel_gemm("N", "N", nrow, nrow, nrow, 1.0_dp, matrix_B, matrix_A, &
                         0.0_dp, work)
      IF (PRESENT(matrix_C)) THEN
         CPASSERT(ASSOCIATED(matrix_C%matrix_struct))
         CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                            0.0_dp, matrix_C)
      ELSE
         CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                            0.0_dp, matrix_B)
      END IF

      CALL cp_fm_release(work)

   END SUBROUTINE fm_contract_ABA

! **************************************************************************************************
!> \brief Computes A^T B_i A for a set of matrices.
!> \param matrix_A left and right matrix A
!> \param matrix_B input matrices B_i; overwritten if matrix_C is absent
!> \param matrix_C optional output matrices C_i=A^T B_i A
! **************************************************************************************************
   SUBROUTINE fm_contract_ABA_set(matrix_A, matrix_B, matrix_C)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_A
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: matrix_B
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: matrix_C

      INTEGER                                            :: i, ncol, nrow
      TYPE(cp_fm_type)                                   :: work

      CPASSERT(ASSOCIATED(matrix_A%matrix_struct))
      CALL cp_fm_get_info(matrix_A, nrow_global=nrow, ncol_global=ncol)
      CPASSERT(nrow == ncol)
      IF (PRESENT(matrix_C)) THEN
         CPASSERT(SIZE(matrix_C) == SIZE(matrix_B))
      END IF
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      DO i = 1, SIZE(matrix_B)
         CPASSERT(ASSOCIATED(matrix_B(i)%matrix_struct))
         CALL parallel_gemm("N", "N", nrow, nrow, nrow, 1.0_dp, matrix_B(i), matrix_A, &
                            0.0_dp, work)
         IF (PRESENT(matrix_C)) THEN
            CPASSERT(ASSOCIATED(matrix_C(i)%matrix_struct))
            CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                               0.0_dp, matrix_C(i))
         ELSE
            CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                               0.0_dp, matrix_B(i))
         END IF
      END DO

      CALL cp_fm_release(work)

   END SUBROUTINE fm_contract_ABA_set

END MODULE gw_utils_fm
