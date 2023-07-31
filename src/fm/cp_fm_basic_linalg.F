!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief basic linear algebra operations for full matrices
!> \par History
!>      08.2002 split out of qs_blacs [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_fm_basic_linalg
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_fm_struct, ONLY: cp_fm_struct_equivalent
   USE cp_fm_types, ONLY: &
      cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
      cp_fm_release, cp_fm_set_all, cp_fm_set_element, cp_fm_set_submatrix, cp_fm_to_fm, &
      cp_fm_type
   USE cp_log_handling, ONLY: cp_logger_get_default_unit_nr, &
                              cp_to_string
   USE kahan_sum, ONLY: accurate_dot_product, &
                        accurate_sum
   USE kinds, ONLY: dp, &
                    int_8, &
                    sp
   USE machine, ONLY: m_memory
   USE mathlib, ONLY: get_pseudo_inverse_svd, &
                      invert_matrix
   USE message_passing, ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_basic_linalg'

   PUBLIC :: cp_fm_scale, & ! scale a matrix
             cp_fm_scale_and_add, & ! scale and add two matrices
             cp_fm_geadd, & ! general addition
             cp_fm_column_scale, & ! scale columns of a matrix
             cp_fm_row_scale, & ! scale rows of a matrix
             cp_fm_trace, & ! trace of the transpose(A)*B
             cp_fm_contracted_trace, & ! sum_{i,...,k} Tr [A(i,...,k)^T * B(i,...,k)]
             cp_fm_norm, & ! different norms of A
             cp_fm_schur_product, & ! schur product
             cp_fm_transpose, & ! transpose a matrix
             cp_fm_upper_to_full, & ! symmetrise an upper symmetric matrix
             cp_fm_syrk, & ! rank k update
             cp_fm_triangular_multiply, & ! triangular matrix multiply / solve
             cp_fm_symm, & ! multiply a symmetric with a non-symmetric matrix
             cp_fm_gemm, & ! multiply two matrices
             cp_complex_fm_gemm, & ! multiply two complex matrices, represented by non_complex fm matrices
             cp_fm_invert, & ! computes the inverse and determinant
             cp_fm_frobenius_norm, & ! frobenius norm
             cp_fm_triangular_invert, & ! compute the reciprocal of a triangular matrix
             cp_fm_qr_factorization, & ! compute the QR factorization of a rectangular matrix
             cp_fm_solve, & ! solves the equation  A*B=C A and C are input
             cp_fm_pdgeqpf, & ! compute a QR factorization with column pivoting of a M-by-N distributed matrix
             cp_fm_pdorgqr, & ! generates an M-by-N as first N columns of a product of K elementary reflectors
             cp_fm_potrf, & ! Cholesky decomposition
             cp_fm_potri, & ! Invert triangular matrix
             cp_fm_rot_rows, & ! rotates two rows
             cp_fm_rot_cols, & ! rotates two columns
             cp_fm_cholesky_restore, & ! apply Cholesky decomposition
             cp_fm_Gram_Schmidt_orthonorm ! Gram-Schmidt orthonormalization of columns of a full matrix

   REAL(KIND=dp), EXTERNAL :: dlange, pdlange, pdlatra
   REAL(KIND=sp), EXTERNAL :: slange, pslange, pslatra

   INTERFACE cp_fm_trace
      MODULE PROCEDURE cp_fm_trace_a0b0t0
      MODULE PROCEDURE cp_fm_trace_a1b0t1_a
      MODULE PROCEDURE cp_fm_trace_a1b0t1_p
      MODULE PROCEDURE cp_fm_trace_a1b1t1_aa
      MODULE PROCEDURE cp_fm_trace_a1b1t1_ap
      MODULE PROCEDURE cp_fm_trace_a1b1t1_pa
      MODULE PROCEDURE cp_fm_trace_a1b1t1_pp
   END INTERFACE cp_fm_trace

   INTERFACE cp_fm_contracted_trace
      MODULE PROCEDURE cp_fm_contracted_trace_a2b2t2_aa
      MODULE PROCEDURE cp_fm_contracted_trace_a2b2t2_ap
      MODULE PROCEDURE cp_fm_contracted_trace_a2b2t2_pa
      MODULE PROCEDURE cp_fm_contracted_trace_a2b2t2_pp
   END INTERFACE cp_fm_contracted_trace
CONTAINS

! **************************************************************************************************
!> \brief calc A <- alpha*A + beta*B
!>      optimized for alpha == 1.0 (just add beta*B) and beta == 0.0 (just
!>      scale A)
!> \param alpha ...
!> \param matrix_a ...
!> \param beta ...
!> \param matrix_b ...
! **************************************************************************************************
   SUBROUTINE cp_fm_scale_and_add(alpha, matrix_a, beta, matrix_b)

      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_a
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: beta
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: matrix_b

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_scale_and_add'

      INTEGER                                            :: handle, size_a, size_b
      REAL(KIND=dp)                                      :: my_beta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b
      REAL(KIND=sp), DIMENSION(:, :), POINTER            :: a_sp, b_sp

      CALL timeset(routineN, handle)

      IF (PRESENT(matrix_b)) THEN
         my_beta = 1.0_dp
      ELSE
         my_beta = 0.0_dp
      END IF
      IF (PRESENT(beta)) my_beta = beta
      NULLIFY (a, b)

      IF (PRESENT(beta)) THEN
         CPASSERT(PRESENT(matrix_b))
         IF (ASSOCIATED(matrix_a%local_data, matrix_b%local_data)) THEN
            CPWARN("Bad use of routine. Call cp_fm_scale instead")
            CALL cp_fm_scale(alpha + beta, matrix_a)
            CALL timestop(handle)
            RETURN
         END IF
      END IF

      a => matrix_a%local_data
      a_sp => matrix_a%local_data_sp

      IF (matrix_a%use_sp) THEN
         size_a = SIZE(a_sp, 1)*SIZE(a_sp, 2)
      ELSE
         size_a = SIZE(a, 1)*SIZE(a, 2)
      END IF

      IF (alpha /= 1.0_dp) THEN
         IF (matrix_a%use_sp) THEN
            CALL sscal(size_a, REAL(alpha, sp), a_sp, 1)
         ELSE
            CALL dscal(size_a, alpha, a, 1)
         END IF
      END IF
      IF (my_beta .NE. 0.0_dp) THEN
         IF (matrix_a%matrix_struct%context /= matrix_b%matrix_struct%context) &
            CPABORT("matrixes must be in the same blacs context")

         IF (cp_fm_struct_equivalent(matrix_a%matrix_struct, &
                                     matrix_b%matrix_struct)) THEN

            b => matrix_b%local_data
            b_sp => matrix_b%local_data_sp
            IF (matrix_b%use_sp) THEN
               size_b = SIZE(b_sp, 1)*SIZE(b_sp, 2)
            ELSE
               size_b = SIZE(b, 1)*SIZE(b, 2)
            END IF
            IF (size_a /= size_b) &
               CPABORT("Matrixes must have same locale sizes")

            IF (matrix_a%use_sp .AND. matrix_b%use_sp) THEN
               CALL saxpy(size_a, REAL(my_beta, sp), b_sp, 1, a_sp, 1)
            ELSEIF (matrix_a%use_sp .AND. .NOT. matrix_b%use_sp) THEN
               CALL saxpy(size_a, REAL(my_beta, sp), REAL(b, sp), 1, a_sp, 1)
            ELSEIF (.NOT. matrix_a%use_sp .AND. matrix_b%use_sp) THEN
               CALL daxpy(size_a, my_beta, REAL(b_sp, dp), 1, a, 1)
            ELSE
               CALL daxpy(size_a, my_beta, b, 1, a, 1)
            END IF

         ELSE
#ifdef __SCALAPACK
            CPABORT("to do (pdscal,pdcopy,pdaxpy)")
#else
            CPABORT("")
#endif
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_scale_and_add

! **************************************************************************************************
!> \brief interface to BLACS geadd:
!>                matrix_b = beta*matrix_b + alpha*opt(matrix_a)
!>        where opt(matrix_a) can be either:
!>              'N':  matrix_a
!>              'T':  matrix_a^T
!>              'C':  matrix_a^H (Hermitian conjugate)
!>        note that this is a level three routine, use cp_fm_scale_and_add if that
!>        is sufficient for your needs
!> \param alpha  : complex scalar
!> \param trans  : 'N' normal, 'T' transposed
!> \param matrix_a : input matrix_a
!> \param beta   : complex scalar
!> \param matrix_b : input matrix_b, upon out put the updated matrix_b
!> \author  Lianheng Tong
! **************************************************************************************************
   SUBROUTINE cp_fm_geadd(alpha, trans, matrix_a, beta, matrix_b)
      REAL(KIND=dp), INTENT(IN) :: alpha, beta
      CHARACTER, INTENT(IN) :: trans
      TYPE(cp_fm_type), INTENT(IN) :: matrix_a, matrix_b

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_geadd'

      INTEGER :: nrow_global, ncol_global, handle
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: aa, bb
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: desca, descb
#else
      INTEGER :: ii, jj
#endif

      CALL timeset(routineN, handle)

      nrow_global = matrix_a%matrix_struct%nrow_global
      ncol_global = matrix_a%matrix_struct%ncol_global
      CPASSERT(nrow_global .EQ. matrix_b%matrix_struct%nrow_global)
      CPASSERT(ncol_global .EQ. matrix_b%matrix_struct%ncol_global)

      aa => matrix_a%local_data
      bb => matrix_b%local_data

#if defined(__SCALAPACK)
      desca = matrix_a%matrix_struct%descriptor
      descb = matrix_b%matrix_struct%descriptor
      CALL pdgeadd(trans, &
                   nrow_global, &
                   ncol_global, &
                   alpha, &
                   aa, &
                   1, 1, &
                   desca, &
                   beta, &
                   bb, &
                   1, 1, &
                   descb)
#else
      ! dgeadd is not a standard BLAS function, although is implemented
      ! in some libraries like OpenBLAS, so not going to use it here
      SELECT CASE (trans)
      CASE ('T')
         DO jj = 1, ncol_global
            DO ii = 1, nrow_global
               bb(ii, jj) = beta*bb(ii, jj) + alpha*aa(jj, ii)
            END DO
         END DO
      CASE DEFAULT
         DO jj = 1, ncol_global
            DO ii = 1, nrow_global
               bb(ii, jj) = beta*bb(ii, jj) + alpha*aa(ii, jj)
            END DO
         END DO
      END SELECT
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_geadd

! **************************************************************************************************
!> \brief Computes the LU-decomposition of the matrix, and the determinant of the matrix
!>      IMPORTANT : the sign of the determinant is not defined correctly yet ....
!> \param matrix_a ...
!> \param almost_determinant ...
!> \param correct_sign ...
!> \par History
!>      added correct_sign 02.07 (fschiff)
!> \author Joost VandeVondele
!> \note
!>      - matrix_a is overwritten
!>      - the sign of the determinant might be wrong
!>      - SERIOUS WARNING (KNOWN BUG) : the sign of the determinant depends on ipivot
!>      - one should be able to find out if ipivot is an even or an odd permutation...
!>        if you need the correct sign, just add correct_sign==.TRUE. (fschiff)
!>      - Use cp_fm_get_diag instead of n times cp_fm_get_element (A. Bussy)
! **************************************************************************************************
   SUBROUTINE cp_fm_lu_decompose(matrix_a, almost_determinant, correct_sign)
      TYPE(cp_fm_type), INTENT(IN)          :: matrix_a
      REAL(KIND=dp), INTENT(OUT)               :: almost_determinant
      LOGICAL, INTENT(IN), OPTIONAL            :: correct_sign

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_lu_decompose'

      INTEGER                                  :: handle, i, info, n
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
      REAL(KIND=dp)                            :: determinant
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca
      REAL(KIND=dp), DIMENSION(:), POINTER     :: diag
#else
      INTEGER                                  :: lda
#endif

      CALL timeset(routineN, handle)

      a => matrix_a%local_data
      n = matrix_a%matrix_struct%nrow_global
      ALLOCATE (ipivot(n + matrix_a%matrix_struct%nrow_block))

#if defined(__SCALAPACK)
      MARK_USED(correct_sign)
      desca(:) = matrix_a%matrix_struct%descriptor(:)
      CALL pdgetrf(n, n, a, 1, 1, desca, ipivot, info)

      ALLOCATE (diag(n))
      diag(:) = 0.0_dp
      CALL cp_fm_get_diag(matrix_a, diag)
      determinant = 1.0_dp
      DO i = 1, n
         determinant = determinant*diag(i)
      END DO
      DEALLOCATE (diag)
#else
      lda = SIZE(a, 1)
      CALL dgetrf(n, n, a, lda, ipivot, info)
      determinant = 1.0_dp
      IF (correct_sign) THEN
         DO i = 1, n
            IF (ipivot(i) .NE. i) THEN
               determinant = -determinant*a(i, i)
            ELSE
               determinant = determinant*a(i, i)
            END IF
         END DO
      ELSE
         DO i = 1, n
            determinant = determinant*a(i, i)
         END DO
      END IF
#endif
      ! info is allowed to be zero
      ! this does just signal a zero diagonal element
      DEALLOCATE (ipivot)
      almost_determinant = determinant ! notice that the sign is random
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief computes matrix_c = beta * matrix_c + alpha * ( matrix_a  ** transa ) * ( matrix_b ** transb )
!> \param transa : 'N' -> normal   'T' -> transpose
!>      alpha,beta :: can be 0.0_dp and 1.0_dp
!> \param transb ...
!> \param m ...
!> \param n ...
!> \param k ...
!> \param alpha ...
!> \param matrix_a : m x k matrix ( ! for transa = 'N')
!> \param matrix_b : k x n matrix ( ! for transb = 'N')
!> \param beta ...
!> \param matrix_c : m x n matrix
!> \param a_first_col ...
!> \param a_first_row ...
!> \param b_first_col : the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
!> \author Matthias Krack
!> \note
!>      matrix_c should have no overlap with matrix_a, matrix_b
! **************************************************************************************************
   SUBROUTINE cp_fm_gemm(transa, transb, m, n, k, alpha, matrix_a, matrix_b, beta, &
                         matrix_c, a_first_col, a_first_row, b_first_col, b_first_row, &
                         c_first_col, c_first_row)

      CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
      INTEGER, INTENT(IN)                      :: m, n, k
      REAL(KIND=dp), INTENT(IN)                :: alpha
      TYPE(cp_fm_type), INTENT(IN)             :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(IN)                :: beta
      TYPE(cp_fm_type), INTENT(IN)          :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL            :: a_first_col, a_first_row, &
                                                  b_first_col, b_first_row, &
                                                  c_first_col, c_first_row

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_gemm'

      INTEGER                                  :: handle, i_a, i_b, i_c, j_a, &
                                                  j_b, j_c
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
      REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp, b_sp, c_sp
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb, descc
#else
      INTEGER                                  :: lda, ldb, ldc
#endif

      CALL timeset(routineN, handle)

      !sample peak memory
      CALL m_memory()

      a => matrix_a%local_data
      b => matrix_b%local_data
      c => matrix_c%local_data

      a_sp => matrix_a%local_data_sp
      b_sp => matrix_b%local_data_sp
      c_sp => matrix_c%local_data_sp

      IF (PRESENT(a_first_row)) THEN
         i_a = a_first_row
      ELSE
         i_a = 1
      END IF
      IF (PRESENT(a_first_col)) THEN
         j_a = a_first_col
      ELSE
         j_a = 1
      END IF
      IF (PRESENT(b_first_row)) THEN
         i_b = b_first_row
      ELSE
         i_b = 1
      END IF
      IF (PRESENT(b_first_col)) THEN
         j_b = b_first_col
      ELSE
         j_b = 1
      END IF
      IF (PRESENT(c_first_row)) THEN
         i_c = c_first_row
      ELSE
         i_c = 1
      END IF
      IF (PRESENT(c_first_col)) THEN
         j_c = c_first_col
      ELSE
         j_c = 1
      END IF

#if defined(__SCALAPACK)

      desca(:) = matrix_a%matrix_struct%descriptor(:)
      descb(:) = matrix_b%matrix_struct%descriptor(:)
      descc(:) = matrix_c%matrix_struct%descriptor(:)

      IF (matrix_a%use_sp .AND. matrix_b%use_sp .AND. matrix_c%use_sp) THEN

         CALL psgemm(transa, transb, m, n, k, REAL(alpha, sp), a_sp(1, 1), i_a, j_a, desca, b_sp(1, 1), i_b, j_b, &
                     descb, REAL(beta, sp), c_sp(1, 1), i_c, j_c, descc)

      ELSEIF ((.NOT. matrix_a%use_sp) .AND. (.NOT. matrix_b%use_sp) .AND. (.NOT. matrix_c%use_sp)) THEN

         CALL pdgemm(transa, transb, m, n, k, alpha, a, i_a, j_a, desca, b, i_b, j_b, &
                     descb, beta, c, i_c, j_c, descc)

      ELSE
         CPABORT("Mixed precision gemm NYI")
      END IF
#else

      IF (matrix_a%use_sp .AND. matrix_b%use_sp .AND. matrix_c%use_sp) THEN

         lda = SIZE(a_sp, 1)
         ldb = SIZE(b_sp, 1)
         ldc = SIZE(c_sp, 1)

         CALL sgemm(transa, transb, m, n, k, REAL(alpha, sp), a_sp(i_a, j_a), lda, b_sp(i_b, j_b), ldb, &
                    REAL(beta, sp), c_sp(i_c, j_c), ldc)

      ELSEIF ((.NOT. matrix_a%use_sp) .AND. (.NOT. matrix_b%use_sp) .AND. (.NOT. matrix_c%use_sp)) THEN

         lda = SIZE(a, 1)
         ldb = SIZE(b, 1)
         ldc = SIZE(c, 1)

         CALL dgemm(transa, transb, m, n, k, alpha, a(i_a, j_a), lda, b(i_b, j_b), ldb, beta, c(i_c, j_c), ldc)

      ELSE
         CPABORT("Mixed precision gemm NYI")
      END IF

#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_gemm

! **************************************************************************************************
!> \brief computes matrix_c = beta * matrix_c + alpha *  matrix_a  *  matrix_b
!>      computes matrix_c = beta * matrix_c + alpha *  matrix_b  *  matrix_a
!>      where matrix_a is symmetric
!> \param side : 'L' -> matrix_a is on the left 'R' -> matrix_a is on the right
!>      alpha,beta :: can be 0.0_dp and 1.0_dp
!> \param uplo ...
!> \param m ...
!> \param n ...
!> \param alpha ...
!> \param matrix_a : m x m matrix
!> \param matrix_b : m x n matrix
!> \param beta ...
!> \param matrix_c : m x n matrix
!> \author Matthias Krack
!> \note
!>      matrix_c should have no overlap with matrix_a, matrix_b
!>      all matrices in QS are upper triangular, so uplo should be 'U' always
!>      matrix_a is always an m x m matrix
!>      it is typically slower to do cp_fm_symm than cp_fm_gemm (especially in parallel easily 50 percent !)
! **************************************************************************************************
   SUBROUTINE cp_fm_symm(side, uplo, m, n, alpha, matrix_a, matrix_b, beta, matrix_c)

      CHARACTER(LEN=1), INTENT(IN)             :: side, uplo
      INTEGER, INTENT(IN)                      :: m, n
      REAL(KIND=dp), INTENT(IN)                :: alpha
      TYPE(cp_fm_type), INTENT(IN)                :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(IN)                :: beta
      TYPE(cp_fm_type), INTENT(IN)          :: matrix_c

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_symm'

      INTEGER                                  :: handle
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb, descc
#else
      INTEGER                                  :: lda, ldb, ldc
#endif

      CALL timeset(routineN, handle)

      a => matrix_a%local_data
      b => matrix_b%local_data
      c => matrix_c%local_data

#if defined(__SCALAPACK)

      desca(:) = matrix_a%matrix_struct%descriptor(:)
      descb(:) = matrix_b%matrix_struct%descriptor(:)
      descc(:) = matrix_c%matrix_struct%descriptor(:)

      CALL pdsymm(side, uplo, m, n, alpha, a(1, 1), 1, 1, desca, b(1, 1), 1, 1, descb, beta, c(1, 1), 1, 1, descc)

#else

      lda = matrix_a%matrix_struct%local_leading_dimension
      ldb = matrix_b%matrix_struct%local_leading_dimension
      ldc = matrix_c%matrix_struct%local_leading_dimension

      CALL dsymm(side, uplo, m, n, alpha, a(1, 1), lda, b(1, 1), ldb, beta, c(1, 1), ldc)

#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_symm

! **************************************************************************************************
!> \brief computes the Frobenius norm of matrix_a
!> \brief computes the Frobenius norm of matrix_a
!> \param matrix_a : m x n matrix
!> \return ...
!> \author VW
! **************************************************************************************************
   FUNCTION cp_fm_frobenius_norm(matrix_a) RESULT(norm)
      TYPE(cp_fm_type), INTENT(IN)             :: matrix_a
      REAL(KIND=dp)                            :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_frobenius_norm'

      INTEGER                                  :: handle, size_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      REAL(KIND=dp), EXTERNAL                  :: DDOT
#if defined(__SCALAPACK)
      TYPE(mp_comm_type)                       :: group
#endif

      CALL timeset(routineN, handle)

      norm = 0.0_dp
      a => matrix_a%local_data
      size_a = SIZE(a, 1)*SIZE(a, 2)
      norm = DDOT(size_a, a(1, 1), 1, a(1, 1), 1)
#if defined(__SCALAPACK)
      group = matrix_a%matrix_struct%para_env
      CALL group%sum(norm)
#endif
      norm = SQRT(norm)

      CALL timestop(handle)

   END FUNCTION cp_fm_frobenius_norm

! **************************************************************************************************
!> \brief performs a rank-k update of a symmetric matrix_c
!>         matrix_c = beta * matrix_c + alpha * matrix_a * transpose ( matrix_a )
!> \param uplo : 'U'   ('L')
!> \param trans : 'N'  ('T')
!> \param k : number of cols to use in matrix_a
!>      ia,ja ::  1,1 (could be used for selecting subblock of a)
!> \param alpha ...
!> \param matrix_a ...
!> \param ia ...
!> \param ja ...
!> \param beta ...
!> \param matrix_c ...
!> \author Matthias Krack
!> \note
!>      In QS uplo should 'U' (upper part updated)
! **************************************************************************************************
   SUBROUTINE cp_fm_syrk(uplo, trans, k, alpha, matrix_a, ia, ja, beta, matrix_c)
      CHARACTER(LEN=1), INTENT(IN)             :: uplo, trans
      INTEGER, INTENT(IN)                      :: k
      REAL(KIND=dp), INTENT(IN)                :: alpha
      TYPE(cp_fm_type), INTENT(IN)             :: matrix_a
      INTEGER, INTENT(IN)                      :: ia, ja
      REAL(KIND=dp), INTENT(IN)                :: beta
      TYPE(cp_fm_type), INTENT(IN)          :: matrix_c

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syrk'

      INTEGER                                  :: handle, n
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, c
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descc
#else
      INTEGER                                  :: lda, ldc
#endif

      CALL timeset(routineN, handle)

      n = matrix_c%matrix_struct%nrow_global

      a => matrix_a%local_data
      c => matrix_c%local_data

#if defined(__SCALAPACK)

      desca(:) = matrix_a%matrix_struct%descriptor(:)
      descc(:) = matrix_c%matrix_struct%descriptor(:)

      CALL pdsyrk(uplo, trans, n, k, alpha, a(1, 1), ia, ja, desca, beta, c(1, 1), 1, 1, descc)

#else

      lda = SIZE(a, 1)
      ldc = SIZE(c, 1)

      CALL dsyrk(uplo, trans, n, k, alpha, a(ia, ja), lda, beta, c(1, 1), ldc)

#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_syrk

! **************************************************************************************************
!> \brief computes the schur product of two matrices
!>       c_ij = a_ij * b_ij
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_fm_schur_product(matrix_a, matrix_b, matrix_c)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_a, matrix_b, matrix_c

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_schur_product'

      INTEGER                                            :: handle, icol_local, irow_local, mypcol, &
                                                            myprow, ncol_local, npcol, nprow, &
                                                            nrow_local
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b, c
      TYPE(cp_blacs_env_type), POINTER                   :: context

      CALL timeset(routineN, handle)

      context => matrix_a%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      a => matrix_a%local_data
      b => matrix_b%local_data
      c => matrix_c%local_data

      nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
      ncol_local = matrix_a%matrix_struct%ncol_locals(mypcol)

      DO icol_local = 1, ncol_local
         DO irow_local = 1, nrow_local
            c(irow_local, icol_local) = a(irow_local, icol_local)*b(irow_local, icol_local)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE cp_fm_schur_product

! **************************************************************************************************
!> \brief returns the trace of matrix_a^T matrix_b, i.e
!>      sum_{i,j}(matrix_a(i,j)*matrix_b(i,j))
!> \param matrix_a a matrix
!> \param matrix_b another matrix
!> \param trace ...
!> \par History
!>      11.06.2001 Creation (Matthias Krack)
!>      12.2002 added doc [fawzi]
!> \author Matthias Krack
!> \note
!>      note the transposition of matrix_a!
! **************************************************************************************************
   SUBROUTINE cp_fm_trace_a0b0t0(matrix_a, matrix_b, trace)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(OUT)                         :: trace

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_trace_a0b0t0'

      INTEGER                                            :: handle, mypcol, myprow, ncol_local, &
                                                            npcol, nprow, nrow_local
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b
      REAL(KIND=sp), DIMENSION(:, :), POINTER            :: a_sp, b_sp
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)

      context => matrix_a%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      group = matrix_a%matrix_struct%para_env

      a => matrix_a%local_data
      b => matrix_b%local_data

      a_sp => matrix_a%local_data_sp
      b_sp => matrix_b%local_data_sp

      nrow_local = MIN(matrix_a%matrix_struct%nrow_locals(myprow), matrix_b%matrix_struct%nrow_locals(myprow))
      ncol_local = MIN(matrix_a%matrix_struct%ncol_locals(mypcol), matrix_b%matrix_struct%ncol_locals(mypcol))

      ! cries for an accurate_dot_product
      IF (matrix_a%use_sp .AND. matrix_b%use_sp) THEN
         trace = accurate_sum(REAL(a_sp(1:nrow_local, 1:ncol_local)* &
                                   b_sp(1:nrow_local, 1:ncol_local), dp))
      ELSEIF (matrix_a%use_sp .AND. .NOT. matrix_b%use_sp) THEN
         trace = accurate_sum(REAL(a_sp(1:nrow_local, 1:ncol_local), dp)* &
                              b(1:nrow_local, 1:ncol_local))
      ELSEIF (.NOT. matrix_a%use_sp .AND. matrix_b%use_sp) THEN
         trace = accurate_sum(a(1:nrow_local, 1:ncol_local)* &
                              REAL(b_sp(1:nrow_local, 1:ncol_local), dp))
      ELSE
         trace = accurate_dot_product(a(1:nrow_local, 1:ncol_local), &
                                      b(1:nrow_local, 1:ncol_local))
      END IF

      CALL group%sum(trace)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_trace_a0b0t0

   #:mute
      #:set types = [("cp_fm_type", "a", ""), ("cp_fm_p_type", "p","%matrix")]
   #:endmute

! **************************************************************************************************
!> \brief Compute trace(k) = Tr (matrix_a(k)^T matrix_b) for each pair of matrices A_k and B.
!> \param matrix_a list of A matrices
!> \param matrix_b B matrix
!> \param trace    computed traces
!> \par History
!>    * 08.2018 forked from cp_fm_trace() [Sergey Chulkov]
!> \note \parblock
!>      Computing the trace requires collective communication between involved MPI processes
!>      that implies a synchronisation point between them. The aim of this subroutine is to reduce
!>      the amount of time wasted in such synchronisation by performing one large collective
!>      operation which involves all the matrices in question.
!>
!>      The subroutine's suffix reflects dimensionality of dummy arrays; 'a1b0t1' means that
!>      the dummy variables 'matrix_a' and 'trace' are 1-dimensional arrays, while the variable
!>      'matrix_b' is a single matrix.
!>      \endparblock
! **************************************************************************************************
   #:for longname, shortname, appendix in types
      SUBROUTINE cp_fm_trace_a1b0t1_${shortname}$ (matrix_a, matrix_b, trace)
         TYPE(${longname}$), DIMENSION(:), INTENT(in)       :: matrix_a
         TYPE(cp_fm_type), INTENT(IN)                       :: matrix_b
         REAL(kind=dp), DIMENSION(:), INTENT(out)           :: trace

         CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_trace_a1b0t1_${shortname}$'

         INTEGER                                            :: handle, imatrix, n_matrices, &
                                                               ncols_local, nrows_local
         LOGICAL                                            :: use_sp_a, use_sp_b
         REAL(kind=dp), DIMENSION(:, :), POINTER            :: ldata_a, ldata_b
         REAL(kind=sp), DIMENSION(:, :), POINTER            :: ldata_a_sp, ldata_b_sp
         TYPE(mp_comm_type)                                 :: group

         CALL timeset(routineN, handle)

         n_matrices = SIZE(trace)
         CPASSERT(SIZE(matrix_a) == n_matrices)

         CALL cp_fm_get_info(matrix_b, nrow_local=nrows_local, ncol_local=ncols_local)
         use_sp_b = matrix_b%use_sp

         IF (use_sp_b) THEN
            ldata_b_sp => matrix_b%local_data_sp(1:nrows_local, 1:ncols_local)
         ELSE
            ldata_b => matrix_b%local_data(1:nrows_local, 1:ncols_local)
         END IF

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(imatrix, ldata_a, ldata_a_sp, use_sp_a), &
!$OMP             SHARED(ldata_b, ldata_b_sp, matrix_a, matrix_b), &
!$OMP             SHARED(ncols_local, nrows_local, n_matrices, trace, use_sp_b)

         DO imatrix = 1, n_matrices

            use_sp_a = matrix_a(imatrix) ${appendix}$%use_sp

            ! assume that the matrices A(i) and B have identical shapes and distribution schemes
            IF (use_sp_a .AND. use_sp_b) THEN
               ldata_a_sp => matrix_a(imatrix) ${appendix}$%local_data_sp(1:nrows_local, 1:ncols_local)
               trace(imatrix) = accurate_dot_product(ldata_a_sp, ldata_b_sp)
            ELSE IF (.NOT. use_sp_a .AND. .NOT. use_sp_b) THEN
               ldata_a => matrix_a(imatrix) ${appendix}$%local_data(1:nrows_local, 1:ncols_local)
               trace(imatrix) = accurate_dot_product(ldata_a, ldata_b)
            ELSE
               CPABORT("Matrices A and B are of different types")
            END IF
         END DO
!$OMP END PARALLEL DO

         group = matrix_b%matrix_struct%para_env
         CALL group%sum(trace)

         CALL timestop(handle)
      END SUBROUTINE cp_fm_trace_a1b0t1_${shortname}$
   #:endfor

! **************************************************************************************************
!> \brief Compute trace(k) = Tr (matrix_a(k)^T matrix_b(k)) for each pair of matrices A_k and B_k.
!> \param matrix_a list of A matrices
!> \param matrix_b list of B matrices
!> \param trace    computed traces
!> \param accurate ...
!> \par History
!>    * 11.2016 forked from cp_fm_trace() [Sergey Chulkov]
!> \note \parblock
!>      Computing the trace requires collective communication between involved MPI processes
!>      that implies a synchronisation point between them. The aim of this subroutine is to reduce
!>      the amount of time wasted in such synchronisation by performing one large collective
!>      operation which involves all the matrices in question.
!>
!>      The subroutine's suffix reflects dimensionality of dummy arrays; 'a1b1t1' means that
!>      all dummy variables (matrix_a, matrix_b, and trace) are 1-dimensional arrays.
!>      \endparblock
! **************************************************************************************************
   #:for longname1, shortname1, appendix1 in types
      #:for longname2, shortname2, appendix2 in types
         SUBROUTINE cp_fm_trace_a1b1t1_${shortname1}$${shortname2}$ (matrix_a, matrix_b, trace, accurate)
            TYPE(${longname1}$), DIMENSION(:), INTENT(in)       :: matrix_a
            TYPE(${longname2}$), DIMENSION(:), INTENT(in)       :: matrix_b
            REAL(kind=dp), DIMENSION(:), INTENT(out)           :: trace
            LOGICAL, INTENT(IN), OPTIONAL                      :: accurate

            CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_trace_a1b1t1_${shortname1}$${shortname2}$'

            INTEGER                                            :: handle, imatrix, n_matrices, &
                                                                  ncols_local, nrows_local
            LOGICAL                                            :: use_accurate_sum, use_sp_a, use_sp_b
            REAL(kind=dp), DIMENSION(:, :), POINTER            :: ldata_a, ldata_b
            REAL(kind=sp), DIMENSION(:, :), POINTER            :: ldata_a_sp, ldata_b_sp
            TYPE(mp_comm_type)                                 :: group

            CALL timeset(routineN, handle)

            n_matrices = SIZE(trace)
            CPASSERT(SIZE(matrix_a) == n_matrices)
            CPASSERT(SIZE(matrix_b) == n_matrices)

            use_accurate_sum = .TRUE.
            IF (PRESENT(accurate)) use_accurate_sum = accurate

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(imatrix, ldata_a, ldata_a_sp, ldata_b, ldata_b_sp, ncols_local), &
!$OMP             PRIVATE(nrows_local, use_sp_a, use_sp_b), &
!$OMP             SHARED(matrix_a, matrix_b, n_matrices, trace, use_accurate_sum)
            DO imatrix = 1, n_matrices
               CALL cp_fm_get_info(matrix_a(imatrix) ${appendix1}$, nrow_local=nrows_local, ncol_local=ncols_local)

               use_sp_a = matrix_a(imatrix) ${appendix1}$%use_sp
               use_sp_b = matrix_b(imatrix) ${appendix2}$%use_sp

               ! assume that the matrices A(i) and B(i) have identical shapes and distribution schemes
               IF (use_sp_a .AND. use_sp_b) THEN
                  ldata_a_sp => matrix_a(imatrix) ${appendix1}$%local_data_sp(1:nrows_local, 1:ncols_local)
                  ldata_b_sp => matrix_b(imatrix) ${appendix2}$%local_data_sp(1:nrows_local, 1:ncols_local)
                  IF (use_accurate_sum) THEN
                     trace(imatrix) = accurate_dot_product(ldata_a_sp, ldata_b_sp)
                  ELSE
                     trace(imatrix) = SUM(ldata_a_sp*ldata_b_sp)
                  END IF
               ELSE IF (.NOT. use_sp_a .AND. .NOT. use_sp_b) THEN
                  ldata_a => matrix_a(imatrix) ${appendix1}$%local_data(1:nrows_local, 1:ncols_local)
                  ldata_b => matrix_b(imatrix) ${appendix2}$%local_data(1:nrows_local, 1:ncols_local)
                  IF (use_accurate_sum) THEN
                     trace(imatrix) = accurate_dot_product(ldata_a, ldata_b)
                  ELSE
                     trace(imatrix) = SUM(ldata_a*ldata_b)
                  END IF
               ELSE
                  CPABORT("Matrices A and B are of different types")
               END IF
            END DO
!$OMP END PARALLEL DO

            group = matrix_a(1) ${appendix1}$%matrix_struct%para_env
            CALL group%sum(trace)

            CALL timestop(handle)
         END SUBROUTINE cp_fm_trace_a1b1t1_${shortname1}$${shortname2}$
      #:endfor
   #:endfor

! **************************************************************************************************
!> \brief Compute trace(i,j) = \sum_k Tr (matrix_a(k,i)^T matrix_b(k,j)).
!> \param matrix_a list of A matrices
!> \param matrix_b list of B matrices
!> \param trace    computed traces
!> \param accurate ...
! **************************************************************************************************
   #:for longname1, shortname1, appendix1 in types
      #:for longname2, shortname2, appendix2 in types
         SUBROUTINE cp_fm_contracted_trace_a2b2t2_${shortname1}$${shortname2}$ (matrix_a, matrix_b, trace, accurate)
            TYPE(${longname1}$), DIMENSION(:, :), INTENT(in)       :: matrix_a
            TYPE(${longname2}$), DIMENSION(:, :), INTENT(in)       :: matrix_b
            REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: trace
            LOGICAL, INTENT(IN), OPTIONAL                      :: accurate

            CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_contracted_trace_a2b2t2_${shortname1}$${shortname2}$'

            INTEGER                                            :: handle, ia, ib, iz, na, nb, ncols_local, &
                                                                  nrows_local, nz
            INTEGER(kind=int_8)                                :: ib8, itrace, na8, ntraces
            LOGICAL                                            :: use_accurate_sum, use_sp_a, use_sp_b
            REAL(kind=dp)                                      :: t
            REAL(kind=dp), DIMENSION(:, :), POINTER            :: ldata_a, ldata_b
            REAL(kind=sp), DIMENSION(:, :), POINTER            :: ldata_a_sp, ldata_b_sp
            TYPE(mp_comm_type)                                 :: group

            CALL timeset(routineN, handle)

            nz = SIZE(matrix_a, 1)
            CPASSERT(SIZE(matrix_b, 1) == nz)

            na = SIZE(matrix_a, 2)
            nb = SIZE(matrix_b, 2)
            CPASSERT(SIZE(trace, 1) == na)
            CPASSERT(SIZE(trace, 2) == nb)

            use_accurate_sum = .TRUE.
            IF (PRESENT(accurate)) use_accurate_sum = accurate

            ! here we use one running index (itrace) instead of two (ia, ib) in order to
            ! improve load balance between shared-memory threads
            ntraces = na*nb
            na8 = INT(na, kind=int_8)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(ia, ib, ib8, itrace, iz, ldata_a, ldata_a_sp, ldata_b, ldata_b_sp, ncols_local), &
!$OMP             PRIVATE(nrows_local, t, use_sp_a, use_sp_b), &
!$OMP             SHARED(matrix_a, matrix_b, na, na8, nb, ntraces, nz, trace, use_accurate_sum)
            DO itrace = 1, ntraces
               ib8 = (itrace - 1)/na8
               ia = INT(itrace - ib8*na8)
               ib = INT(ib8) + 1

               t = 0.0_dp
               DO iz = 1, nz
                  CALL cp_fm_get_info(matrix_a(iz, ia) ${appendix1}$, nrow_local=nrows_local, ncol_local=ncols_local)
                  use_sp_a = matrix_a(iz, ia) ${appendix1}$%use_sp
                  use_sp_b = matrix_b(iz, ib) ${appendix2}$%use_sp

                  ! assume that the matrices A(iz, ia) and B(iz, ib) have identical shapes and distribution schemes
                  IF (.NOT. use_sp_a .AND. .NOT. use_sp_b) THEN
                     ldata_a => matrix_a(iz, ia) ${appendix1}$%local_data(1:nrows_local, 1:ncols_local)
                     ldata_b => matrix_b(iz, ib) ${appendix2}$%local_data(1:nrows_local, 1:ncols_local)
                     IF (use_accurate_sum) THEN
                        t = t + accurate_dot_product(ldata_a, ldata_b)
                     ELSE
                        t = t + SUM(ldata_a*ldata_b)
                     END IF
                  ELSE IF (use_sp_a .AND. use_sp_b) THEN
                     ldata_a_sp => matrix_a(iz, ia) ${appendix1}$%local_data_sp(1:nrows_local, 1:ncols_local)
                     ldata_b_sp => matrix_b(iz, ib) ${appendix2}$%local_data_sp(1:nrows_local, 1:ncols_local)
                     IF (use_accurate_sum) THEN
                        t = t + accurate_dot_product(ldata_a_sp, ldata_b_sp)
                     ELSE
                        t = t + SUM(ldata_a_sp*ldata_b_sp)
                     END IF
                  ELSE
                     CPABORT("Matrices A and B are of different types")
                  END IF
               END DO
               trace(ia, ib) = t
            END DO
!$OMP END PARALLEL DO

            group = matrix_a(1, 1) ${appendix1}$%matrix_struct%para_env
            CALL group%sum(trace)

            CALL timestop(handle)
         END SUBROUTINE cp_fm_contracted_trace_a2b2t2_${shortname1}$${shortname2}$
      #:endfor
   #:endfor

! **************************************************************************************************
!> \brief multiplies in place by a triangular matrix:
!>       matrix_b = alpha op(triangular_matrix) matrix_b
!>      or (if side='R')
!>       matrix_b = alpha matrix_b op(triangular_matrix)
!>      op(triangular_matrix) is:
!>       triangular_matrix (if transpose_tr=.false. and invert_tr=.false.)
!>       triangular_matrix^T (if transpose_tr=.true. and invert_tr=.false.)
!>       triangular_matrix^(-1) (if transpose_tr=.false. and invert_tr=.true.)
!>       triangular_matrix^(-T) (if transpose_tr=.true. and invert_tr=.true.)
!> \param triangular_matrix the triangular matrix that multiplies the other
!> \param matrix_b the matrix that gets multiplied and stores the result
!> \param side on which side of matrix_b stays op(triangular_matrix)
!>        (defaults to 'L')
!> \param transpose_tr if the triangular matrix should be transposed
!>        (defaults to false)
!> \param invert_tr if the triangular matrix should be inverted
!>        (defaults to false)
!> \param uplo_tr if triangular_matrix is stored in the upper ('U') or
!>        lower ('L') triangle (defaults to 'U')
!> \param unit_diag_tr if the diagonal elements of triangular_matrix should
!>        be assumed to be 1 (defaults to false)
!> \param n_rows the number of rows of the result (defaults to
!>        size(matrix_b,1))
!> \param n_cols the number of columns of the result (defaults to
!>        size(matrix_b,2))
!> \param alpha ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs an mpi env
! **************************************************************************************************
   SUBROUTINE cp_fm_triangular_multiply(triangular_matrix, matrix_b, side, &
                                        transpose_tr, invert_tr, uplo_tr, unit_diag_tr, n_rows, n_cols, &
                                        alpha)
      TYPE(cp_fm_type), INTENT(IN)                       :: triangular_matrix, matrix_b
      CHARACTER, INTENT(in), OPTIONAL                    :: side
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose_tr, invert_tr
      CHARACTER, INTENT(in), OPTIONAL                    :: uplo_tr
      LOGICAL, INTENT(in), OPTIONAL                      :: unit_diag_tr
      INTEGER, INTENT(in), OPTIONAL                      :: n_rows, n_cols
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: alpha

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_triangular_multiply'

      CHARACTER                                          :: side_char, transa, unit_diag, uplo
      INTEGER                                            :: handle, m, n
      LOGICAL                                            :: invert
      REAL(KIND=dp)                                      :: al

      CALL timeset(routineN, handle)
      side_char = 'L'
      unit_diag = 'N'
      uplo = 'U'
      transa = 'N'
      invert = .FALSE.
      al = 1.0_dp
      CALL cp_fm_get_info(matrix_b, nrow_global=m, ncol_global=n)
      IF (PRESENT(side)) side_char = side
      IF (PRESENT(invert_tr)) invert = invert_tr
      IF (PRESENT(uplo_tr)) uplo = uplo_tr
      IF (PRESENT(unit_diag_tr)) THEN
         IF (unit_diag_tr) THEN
            unit_diag = 'U'
         ELSE
            unit_diag = 'N'
         END IF
      END IF
      IF (PRESENT(transpose_tr)) THEN
         IF (transpose_tr) THEN
            transa = 'T'
         ELSE
            transa = 'N'
         END IF
      END IF
      IF (PRESENT(alpha)) al = alpha
      IF (PRESENT(n_rows)) m = n_rows
      IF (PRESENT(n_cols)) n = n_cols

      IF (invert) THEN

#if defined(__SCALAPACK)
         CALL pdtrsm(side_char, uplo, transa, unit_diag, m, n, al, &
                     triangular_matrix%local_data(1, 1), 1, 1, &
                     triangular_matrix%matrix_struct%descriptor, &
                     matrix_b%local_data(1, 1), 1, 1, &
                     matrix_b%matrix_struct%descriptor(1))
#else
         CALL dtrsm(side_char, uplo, transa, unit_diag, m, n, al, &
                    triangular_matrix%local_data(1, 1), &
                    SIZE(triangular_matrix%local_data, 1), &
                    matrix_b%local_data(1, 1), SIZE(matrix_b%local_data, 1))
#endif

      ELSE

#if defined(__SCALAPACK)
         CALL pdtrmm(side_char, uplo, transa, unit_diag, m, n, al, &
                     triangular_matrix%local_data(1, 1), 1, 1, &
                     triangular_matrix%matrix_struct%descriptor, &
                     matrix_b%local_data(1, 1), 1, 1, &
                     matrix_b%matrix_struct%descriptor(1))
#else
         CALL dtrmm(side_char, uplo, transa, unit_diag, m, n, al, &
                    triangular_matrix%local_data(1, 1), &
                    SIZE(triangular_matrix%local_data, 1), &
                    matrix_b%local_data(1, 1), SIZE(matrix_b%local_data, 1))
#endif

      END IF

      CALL timestop(handle)
   END SUBROUTINE cp_fm_triangular_multiply

! **************************************************************************************************
!> \brief scales a matrix
!>      matrix_a = alpha * matrix_b
!> \param alpha ...
!> \param matrix_a ...
!> \note
!>      use cp_fm_set_all to zero (avoids problems with nan)
! **************************************************************************************************
   SUBROUTINE cp_fm_scale(alpha, matrix_a)
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_a

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_fm_scale'

      INTEGER                                            :: handle, size_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a

      CALL timeset(routineN, handle)

      NULLIFY (a)

      a => matrix_a%local_data
      size_a = SIZE(a, 1)*SIZE(a, 2)

      CALL DSCAL(size_a, alpha, a, 1)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_scale

! **************************************************************************************************
!> \brief transposes a matrix
!>      matrixt = matrix ^ T
!> \param matrix ...
!> \param matrixt ...
!> \note
!>      all matrix elements are transposed (see cp_fm_upper_to_half to symmetrise a matrix)
! **************************************************************************************************
   SUBROUTINE cp_fm_transpose(matrix, matrixt)
      TYPE(cp_fm_type), INTENT(IN)          :: matrix, matrixt

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_transpose'

      INTEGER                                  :: handle, ncol_global, &
                                                  nrow_global, ncol_globalt, nrow_globalt
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, c
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descc
#else
      INTEGER                                  :: i, j
#endif

      nrow_global = matrix%matrix_struct%nrow_global
      ncol_global = matrix%matrix_struct%ncol_global
      nrow_globalt = matrixt%matrix_struct%nrow_global
      ncol_globalt = matrixt%matrix_struct%ncol_global
      CPASSERT(nrow_global == ncol_globalt)
      CPASSERT(nrow_globalt == ncol_global)

      CALL timeset(routineN, handle)

      a => matrix%local_data
      c => matrixt%local_data

#if defined(__SCALAPACK)
      desca(:) = matrix%matrix_struct%descriptor(:)
      descc(:) = matrixt%matrix_struct%descriptor(:)
      CALL pdtran(ncol_global, nrow_global, 1.0_dp, a(1, 1), 1, 1, desca, 0.0_dp, c(1, 1), 1, 1, descc)
#else
      DO j = 1, ncol_global
         DO i = 1, nrow_global
            c(j, i) = a(i, j)
         END DO
      END DO
#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_transpose

! **************************************************************************************************
!> \brief given an upper triangular matrix computes the corresponding full matrix
!> \param matrix the upper triangular matrix as input, the full matrix as output
!> \param work a matrix of the same size as matrix
!> \author Matthias Krack
!> \note
!>       the lower triangular part is irrelevant
! **************************************************************************************************
   SUBROUTINE cp_fm_upper_to_full(matrix, work)

      TYPE(cp_fm_type), INTENT(IN)          :: matrix, work

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_upper_to_full'

      INTEGER                                  :: handle, icol_global, irow_global, &
                                                  mypcol, myprow, ncol_global, &
                                                  npcol, nprow, nrow_global
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER   :: a_sp
      TYPE(cp_blacs_env_type), POINTER         :: context

#if defined(__SCALAPACK)
      INTEGER                                  :: icol_local, irow_local, &
                                                  ncol_block, ncol_local, &
                                                  nrow_block, nrow_local
      INTEGER, DIMENSION(9)                    :: desca, descc
      INTEGER, EXTERNAL                        :: indxl2g
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: c
      REAL(KIND=sp), DIMENSION(:, :), POINTER   :: c_sp
#endif

      nrow_global = matrix%matrix_struct%nrow_global
      ncol_global = matrix%matrix_struct%ncol_global
      CPASSERT(nrow_global == ncol_global)
      nrow_global = work%matrix_struct%nrow_global
      ncol_global = work%matrix_struct%ncol_global
      CPASSERT(nrow_global == ncol_global)
      CPASSERT(matrix%use_sp .EQV. work%use_sp)

      CALL timeset(routineN, handle)

      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

#if defined(__SCALAPACK)

      nrow_block = matrix%matrix_struct%nrow_block
      ncol_block = matrix%matrix_struct%ncol_block

      nrow_local = matrix%matrix_struct%nrow_locals(myprow)
      ncol_local = matrix%matrix_struct%ncol_locals(mypcol)

      a => work%local_data
      a_sp => work%local_data_sp
      desca(:) = work%matrix_struct%descriptor(:)
      c => matrix%local_data
      c_sp => matrix%local_data_sp
      descc(:) = matrix%matrix_struct%descriptor(:)

      DO icol_local = 1, ncol_local
         icol_global = indxl2g(icol_local, ncol_block, mypcol, &
                               matrix%matrix_struct%first_p_pos(2), npcol)
         DO irow_local = 1, nrow_local
            irow_global = indxl2g(irow_local, nrow_block, myprow, &
                                  matrix%matrix_struct%first_p_pos(1), nprow)
            IF (irow_global > icol_global) THEN
               IF (matrix%use_sp) THEN
                  c_sp(irow_local, icol_local) = 0.0_sp
               ELSE
                  c(irow_local, icol_local) = 0.0_dp
               END IF
            ELSE IF (irow_global == icol_global) THEN
               IF (matrix%use_sp) THEN
                  c_sp(irow_local, icol_local) = 0.5_sp*c_sp(irow_local, icol_local)
               ELSE
                  c(irow_local, icol_local) = 0.5_dp*c(irow_local, icol_local)
               END IF
            END IF
         END DO
      END DO

      DO icol_local = 1, ncol_local
      DO irow_local = 1, nrow_local
         IF (matrix%use_sp) THEN
            a_sp(irow_local, icol_local) = c_sp(irow_local, icol_local)
         ELSE
            a(irow_local, icol_local) = c(irow_local, icol_local)
         END IF
      END DO
      END DO

      IF (matrix%use_sp) THEN
         CALL pstran(nrow_global, ncol_global, 1.0_sp, a_sp(1, 1), 1, 1, desca, 1.0_sp, c_sp(1, 1), 1, 1, descc)
      ELSE
         CALL pdtran(nrow_global, ncol_global, 1.0_dp, a(1, 1), 1, 1, desca, 1.0_dp, c(1, 1), 1, 1, descc)
      END IF

#else

      a => matrix%local_data
      a_sp => matrix%local_data_sp
      DO irow_global = 1, nrow_global
         DO icol_global = irow_global + 1, ncol_global
            IF (matrix%use_sp) THEN
               a_sp(icol_global, irow_global) = a_sp(irow_global, icol_global)
            ELSE
               a(icol_global, irow_global) = a(irow_global, icol_global)
            END IF
         END DO
      END DO

#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_upper_to_full

! **************************************************************************************************
!> \brief scales column i of matrix a with scaling(i)
!> \param matrixa ...
!> \param scaling : an array used for scaling the columns,
!>                  SIZE(scaling) determines the number of columns to be scaled
!> \author Joost VandeVondele
!> \note
!>      this is very useful as a first step in the computation of C = sum_i alpha_i A_i transpose (A_i)
!>      that is a rank-k update (cp_fm_syrk , cp_sm_plus_fm_fm_t)
!>      this procedure can be up to 20 times faster than calling cp_fm_syrk n times
!>      where every vector has a different prefactor
! **************************************************************************************************
   SUBROUTINE cp_fm_column_scale(matrixa, scaling)
      TYPE(cp_fm_type), INTENT(IN)          :: matrixa
      REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: scaling

      INTEGER                                  :: k, mypcol, myprow, n, ncol_global, &
                                                  npcol, nprow
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
#if defined(__SCALAPACK)
      INTEGER                                  :: icol_global, icol_local, &
                                                  ipcol, iprow, irow_local
#else
      INTEGER                                  :: i
#endif

      myprow = matrixa%matrix_struct%context%mepos(1)
      mypcol = matrixa%matrix_struct%context%mepos(2)
      nprow = matrixa%matrix_struct%context%num_pe(1)
      npcol = matrixa%matrix_struct%context%num_pe(2)

      ncol_global = matrixa%matrix_struct%ncol_global

      a => matrixa%local_data
      a_sp => matrixa%local_data_sp
      IF (matrixa%use_sp) THEN
         n = SIZE(a_sp, 1)
      ELSE
         n = SIZE(a, 1)
      END IF
      k = MIN(SIZE(scaling), ncol_global)

#if defined(__SCALAPACK)

      DO icol_global = 1, k
         CALL infog2l(1, icol_global, matrixa%matrix_struct%descriptor, &
                      nprow, npcol, myprow, mypcol, &
                      irow_local, icol_local, iprow, ipcol)
         IF ((ipcol == mypcol)) THEN
            IF (matrixa%use_sp) THEN
               CALL SSCAL(n, REAL(scaling(icol_global), sp), a_sp(:, icol_local), 1)
            ELSE
               CALL DSCAL(n, scaling(icol_global), a(:, icol_local), 1)
            END IF
         END IF
      END DO
#else
      DO i = 1, k
         IF (matrixa%use_sp) THEN
            CALL SSCAL(n, REAL(scaling(i), sp), a_sp(:, i), 1)
         ELSE
            CALL DSCAL(n, scaling(i), a(:, i), 1)
         END IF
      END DO
#endif
   END SUBROUTINE cp_fm_column_scale

! **************************************************************************************************
!> \brief scales row i of matrix a with scaling(i)
!> \param matrixa ...
!> \param scaling : an array used for scaling the columns,
!> \author JGH
!> \note
! **************************************************************************************************
   SUBROUTINE cp_fm_row_scale(matrixa, scaling)
      TYPE(cp_fm_type), INTENT(IN)          :: matrixa
      REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: scaling

      INTEGER                                  :: n, m, nrow_global, nrow_local, ncol_local
      INTEGER, DIMENSION(:), POINTER           :: row_indices
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
#if defined(__SCALAPACK)
      INTEGER                                  :: irow_global, icol, irow
#else
      INTEGER                                  :: j
#endif

      CALL cp_fm_get_info(matrixa, row_indices=row_indices, nrow_global=nrow_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local)
      CPASSERT(SIZE(scaling) == nrow_global)

      a => matrixa%local_data
      a_sp => matrixa%local_data_sp
      IF (matrixa%use_sp) THEN
         n = SIZE(a_sp, 1)
         m = SIZE(a_sp, 2)
      ELSE
         n = SIZE(a, 1)
         m = SIZE(a, 2)
      END IF

#if defined(__SCALAPACK)
      DO icol = 1, ncol_local
         IF (matrixa%use_sp) THEN
            DO irow = 1, nrow_local
               irow_global = row_indices(irow)
               a(irow, icol) = REAL(scaling(irow_global), dp)*a(irow, icol)
            END DO
         ELSE
            DO irow = 1, nrow_local
               irow_global = row_indices(irow)
               a(irow, icol) = scaling(irow_global)*a(irow, icol)
            END DO
         END IF
      END DO
#else
      IF (matrixa%use_sp) THEN
         DO j = 1, m
            a_sp(1:n, j) = REAL(scaling(1:n), sp)*a_sp(1:n, j)
         END DO
      ELSE
         DO j = 1, m
            a(1:n, j) = scaling(1:n)*a(1:n, j)
         END DO
      END IF
#endif
   END SUBROUTINE cp_fm_row_scale
! **************************************************************************************************
!> \brief Inverts a cp_fm_type matrix, optionally returning the determinant of the input matrix
!> \param matrix_a the matrix to invert
!> \param matrix_inverse the inverse of matrix_a
!> \param det_a the determinant of matrix_a
!> \param eps_svd optional parameter to active SVD based inversion, singular values below eps_svd
!>                are screened
!> \param eigval optionally return matrix eigenvalues/singular values
!> \par History
!>      note of Jan Wilhelm (12.2015)
!>      - computation of determinant corrected
!>      - determinant only computed if det_a is present
!>      12.2016 added option to use SVD instead of LU [Nico Holmberg]
!>      - Use cp_fm_get diag instead of n times cp_fm_get_element (A. Bussy)
!> \author Florian Schiffmann(02.2007)
! **************************************************************************************************
   SUBROUTINE cp_fm_invert(matrix_a, matrix_inverse, det_a, eps_svd, eigval)

      TYPE(cp_fm_type), INTENT(IN)          :: matrix_a, matrix_inverse
      REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: det_a
      REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eps_svd
      REAL(KIND=dp), DIMENSION(:), POINTER, &
         INTENT(INOUT), OPTIONAL               :: eigval

      INTEGER                                  :: n
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
      REAL(KIND=dp)                            :: determinant, my_eps_svd
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      TYPE(cp_fm_type)                :: matrix_lu

#if defined(__SCALAPACK)
      TYPE(cp_fm_type)                :: u, vt, sigma, inv_sigma_ut
      TYPE(mp_comm_type) :: group
      INTEGER                                  :: i, info, liwork, lwork, exponent_of_minus_one
      INTEGER, DIMENSION(9)                    :: desca
      LOGICAL                                  :: quenched
      REAL(KIND=dp)                            :: alpha, beta
      REAL(KIND=dp), DIMENSION(:), POINTER     :: diag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
#else
      LOGICAL                                  :: sign
      REAL(KIND=dp)                            :: eps1
#endif

      my_eps_svd = 0.0_dp
      IF (PRESENT(eps_svd)) my_eps_svd = eps_svd

      CALL cp_fm_create(matrix=matrix_lu, &
                        matrix_struct=matrix_a%matrix_struct, &
                        name="A_lu"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX")
      CALL cp_fm_to_fm(matrix_a, matrix_lu)

      a => matrix_lu%local_data
      n = matrix_lu%matrix_struct%nrow_global
      ALLOCATE (ipivot(n + matrix_a%matrix_struct%nrow_block))
      ipivot(:) = 0
#if defined(__SCALAPACK)
      IF (my_eps_svd .EQ. 0.0_dp) THEN
         ! Use LU decomposition
         lwork = 3*n
         liwork = 3*n
         desca(:) = matrix_lu%matrix_struct%descriptor(:)
         CALL pdgetrf(n, n, a, 1, 1, desca, ipivot, info)

         IF (PRESENT(det_a) .OR. PRESENT(eigval)) THEN

            ALLOCATE (diag(n))
            diag(:) = 0.0_dp
            CALL cp_fm_get_diag(matrix_lu, diag)

            exponent_of_minus_one = 0
            determinant = 1.0_dp
            DO i = 1, n
               determinant = determinant*diag(i)
               IF (ipivot(i) .NE. i) THEN
                  exponent_of_minus_one = exponent_of_minus_one + 1
               END IF
            END DO
            IF (PRESENT(eigval)) THEN
               CPASSERT(.NOT. ASSOCIATED(eigval))
               ALLOCATE (eigval(n))
               eigval(:) = diag
            END IF
            DEALLOCATE (diag)

            group = matrix_lu%matrix_struct%para_env
            CALL group%sum(exponent_of_minus_one)

            determinant = determinant*(-1.0_dp)**exponent_of_minus_one

         END IF

         alpha = 0.0_dp
         beta = 1.0_dp
         CALL cp_fm_set_all(matrix_inverse, alpha, beta)
         CALL pdgetrs('N', n, n, matrix_lu%local_data, 1, 1, desca, ipivot, matrix_inverse%local_data, 1, 1, desca, info)
      ELSE
         ! Use singular value decomposition
         CALL cp_fm_create(matrix=u, &
                           matrix_struct=matrix_a%matrix_struct, &
                           name="LEFT_SINGULAR_MATRIX")
         CALL cp_fm_set_all(u, alpha=0.0_dp)
         CALL cp_fm_create(matrix=vt, &
                           matrix_struct=matrix_a%matrix_struct, &
                           name="RIGHT_SINGULAR_MATRIX")
         CALL cp_fm_set_all(vt, alpha=0.0_dp)
         ALLOCATE (diag(n))
         diag(:) = 0.0_dp
         desca(:) = matrix_lu%matrix_struct%descriptor(:)
         ALLOCATE (work(1))
         ! Workspace query
         lwork = -1
         CALL pdgesvd('V', 'V', n, n, matrix_lu%local_data, 1, 1, desca, diag, u%local_data, &
                      1, 1, desca, vt%local_data, 1, 1, desca, work, lwork, info)
         lwork = INT(work(1))
         DEALLOCATE (work)
         ALLOCATE (work(lwork))
         ! SVD
         CALL pdgesvd('V', 'V', n, n, matrix_lu%local_data, 1, 1, desca, diag, u%local_data, &
                      1, 1, desca, vt%local_data, 1, 1, desca, work, lwork, info)
         ! info == n+1 implies homogeneity error when the number of procs is large
         ! this likely isnt a problem, but maybe we should handle it separately
         IF (info /= 0 .AND. info /= n + 1) &
            CPABORT("Singular value decomposition of matrix failed.")
         ! (Pseudo)inverse and (pseudo)determinant
         CALL cp_fm_create(matrix=sigma, &
                           matrix_struct=matrix_a%matrix_struct, &
                           name="SINGULAR_VALUE_MATRIX")
         CALL cp_fm_set_all(sigma, alpha=0.0_dp)
         determinant = 1.0_dp
         quenched = .FALSE.
         IF (PRESENT(eigval)) THEN
            CPASSERT(.NOT. ASSOCIATED(eigval))
            ALLOCATE (eigval(n))
            eigval(:) = diag
         END IF
         DO i = 1, n
            IF (diag(i) < my_eps_svd) THEN
               diag(i) = 0.0_dp
               quenched = .TRUE.
            ELSE
               determinant = determinant*diag(i)
               diag(i) = 1.0_dp/diag(i)
            END IF
            CALL cp_fm_set_element(sigma, i, i, diag(i))
         END DO
         DEALLOCATE (diag)
         IF (quenched) &
            CALL cp_warn(__LOCATION__, &
                         "Linear dependencies were detected in the SVD inversion of matrix "//TRIM(ADJUSTL(matrix_a%name))// &
                         ". At least one singular value has been quenched.")
         ! Sigma^-1 * U^T
         CALL cp_fm_create(matrix=inv_sigma_ut, &
                           matrix_struct=matrix_a%matrix_struct, &
                           name="SINGULAR_VALUE_MATRIX")
         CALL cp_fm_set_all(inv_sigma_ut, alpha=0.0_dp)
         CALL pdgemm('N', 'T', n, n, n, 1.0_dp, sigma%local_data, 1, 1, desca, &
                     u%local_data, 1, 1, desca, 0.0_dp, inv_sigma_ut%local_data, 1, 1, desca)
         ! A^-1 = V * (Sigma^-1 * U^T)
         CALL cp_fm_set_all(matrix_inverse, alpha=0.0_dp)
         CALL pdgemm('T', 'N', n, n, n, 1.0_dp, vt%local_data, 1, 1, desca, &
                     inv_sigma_ut%local_data, 1, 1, desca, 0.0_dp, matrix_inverse%local_data, 1, 1, desca)
         ! Clean up
         DEALLOCATE (work)
         CALL cp_fm_release(u)
         CALL cp_fm_release(vt)
         CALL cp_fm_release(sigma)
         CALL cp_fm_release(inv_sigma_ut)
      END IF
#else
      IF (my_eps_svd .EQ. 0.0_dp) THEN
         sign = .TRUE.
         CALL invert_matrix(matrix_a%local_data, matrix_inverse%local_data, &
                            eval_error=eps1)
         CALL cp_fm_lu_decompose(matrix_lu, determinant, correct_sign=sign)
         IF (PRESENT(eigval)) &
            CALL cp_abort(__LOCATION__, &
                          "NYI. Eigenvalues not available for return without SCALAPACK.")
      ELSE
         CALL get_pseudo_inverse_svd(matrix_a%local_data, matrix_inverse%local_data, eps_svd, &
                                     determinant, eigval)
      END IF
#endif
      CALL cp_fm_release(matrix_lu)
      DEALLOCATE (ipivot)
      IF (PRESENT(det_a)) det_a = determinant
   END SUBROUTINE cp_fm_invert

! **************************************************************************************************
!> \brief inverts a triangular matrix
!> \param matrix_a ...
!> \param uplo_tr ...
!> \author MI
! **************************************************************************************************
   SUBROUTINE cp_fm_triangular_invert(matrix_a, uplo_tr)

      TYPE(cp_fm_type), INTENT(IN)          :: matrix_a
      CHARACTER, INTENT(IN), OPTIONAL          :: uplo_tr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_triangular_invert'

      CHARACTER                                :: unit_diag, uplo
      INTEGER                                  :: handle, info, ncol_global
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca
#endif

      CALL timeset(routineN, handle)

      unit_diag = 'N'
      uplo = 'U'
      IF (PRESENT(uplo_tr)) uplo = uplo_tr

      ncol_global = matrix_a%matrix_struct%ncol_global

      a => matrix_a%local_data

#if defined(__SCALAPACK)

      desca(:) = matrix_a%matrix_struct%descriptor(:)

      CALL pdtrtri(uplo, unit_diag, ncol_global, a(1, 1), 1, 1, desca, info)

#else
      CALL dtrtri(uplo, unit_diag, ncol_global, a(1, 1), ncol_global, info)
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_triangular_invert

! **************************************************************************************************
!> \brief  performs a QR factorization of the input rectangular matrix A or of a submatrix of A
!>         the computed upper triangular matrix R is in output in the submatrix sub(A) of size NxN
!>         M and M give the dimension of the submatrix that has to be factorized (MxN) with M>N
!> \param matrix_a ...
!> \param matrix_r ...
!> \param nrow_fact ...
!> \param ncol_fact ...
!> \param first_row ...
!> \param first_col ...
!> \author MI
! **************************************************************************************************
   SUBROUTINE cp_fm_qr_factorization(matrix_a, matrix_r, nrow_fact, ncol_fact, first_row, first_col)
      TYPE(cp_fm_type), INTENT(IN)          :: matrix_a, matrix_r
      INTEGER, INTENT(IN), OPTIONAL            :: nrow_fact, ncol_fact, &
                                                  first_row, first_col

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_qr_factorization'

      INTEGER                                  :: handle, i, icol, info, irow, &
                                                  j, lda, lwork, ncol, &
                                                  ndim, nrow
      REAL(dp), ALLOCATABLE, DIMENSION(:)      :: tau, work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: r_mat
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca
#endif

      CALL timeset(routineN, handle)

      ncol = matrix_a%matrix_struct%ncol_global
      nrow = matrix_a%matrix_struct%nrow_global
      lda = nrow

      a => matrix_a%local_data

      IF (PRESENT(nrow_fact)) nrow = nrow_fact
      IF (PRESENT(ncol_fact)) ncol = ncol_fact
      irow = 1
      IF (PRESENT(first_row)) irow = first_row
      icol = 1
      IF (PRESENT(first_col)) icol = first_col

      CPASSERT(nrow >= ncol)
      ndim = SIZE(a, 2)
!    ALLOCATE(ipiv(ndim))
      ALLOCATE (tau(ndim))

#if defined(__SCALAPACK)

      desca(:) = matrix_a%matrix_struct%descriptor(:)

      lwork = -1
      ALLOCATE (work(2*ndim))
      CALL pdgeqrf(nrow, ncol, a, irow, icol, desca, tau, work, lwork, info)
      lwork = INT(work(1))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))
      CALL pdgeqrf(nrow, ncol, a, irow, icol, desca, tau, work, lwork, info)

#else
      lwork = -1
      ALLOCATE (work(2*ndim))
      CALL dgeqrf(nrow, ncol, a, lda, tau, work, lwork, info)
      lwork = INT(work(1))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))
      CALL dgeqrf(nrow, ncol, a, lda, tau, work, lwork, info)

#endif

      ALLOCATE (r_mat(ncol, ncol))
      CALL cp_fm_get_submatrix(matrix_a, r_mat, 1, 1, ncol, ncol)
      DO i = 1, ncol
         DO j = i + 1, ncol
            r_mat(j, i) = 0.0_dp
         END DO
      END DO
      CALL cp_fm_set_submatrix(matrix_r, r_mat, 1, 1, ncol, ncol)

      DEALLOCATE (tau, work, r_mat)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_qr_factorization

! **************************************************************************************************
!> \brief computes the the solution to A*b=A_general using lu decomposition
!>        pay attention, both matrices are overwritten, a_general contais the result
!> \param matrix_a ...
!> \param general_a ...
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE cp_fm_solve(matrix_a, general_a)
      TYPE(cp_fm_type), INTENT(IN)          :: matrix_a, general_a

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_solve'

      INTEGER                                  :: handle, info, n
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, a_general
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb
#else
      INTEGER                                  :: lda, ldb
#endif

      CALL timeset(routineN, handle)

      a => matrix_a%local_data
      a_general => general_a%local_data
      n = matrix_a%matrix_struct%nrow_global
      ALLOCATE (ipivot(n + matrix_a%matrix_struct%nrow_block))

#if defined(__SCALAPACK)
      desca(:) = matrix_a%matrix_struct%descriptor(:)
      descb(:) = general_a%matrix_struct%descriptor(:)
      CALL pdgetrf(n, n, a, 1, 1, desca, ipivot, info)
      CALL pdgetrs("N", n, n, a, 1, 1, desca, ipivot, a_general, &
                   1, 1, descb, info)

#else
      lda = SIZE(a, 1)
      ldb = SIZE(a_general, 1)
      CALL dgetrf(n, n, a, lda, ipivot, info)
      CALL dgetrs("N", n, n, a, lda, ipivot, a_general, ldb, info)

#endif
      ! info is allowed to be zero
      ! this does just signal a zero diagonal element
      DEALLOCATE (ipivot)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Convenience function. Computes the matrix multiplications needed
!>        for the multiplication of complex matrices.
!>        C = beta * C + alpha * ( A  ** transa ) * ( B ** transb )
!> \param transa : 'N' -> normal   'T' -> transpose
!>      alpha,beta :: can be 0.0_dp and 1.0_dp
!> \param transb ...
!> \param m ...
!> \param n ...
!> \param k ...
!> \param alpha ...
!> \param A_re m x k matrix ( ! for transa = 'N'), real part
!> \param A_im m x k matrix ( ! for transa = 'N'), imaginary part
!> \param B_re k x n matrix ( ! for transa = 'N'), real part
!> \param B_im k x n matrix ( ! for transa = 'N'), imaginary part
!> \param beta ...
!> \param C_re m x n matrix, real part
!> \param C_im m x n matrix, imaginary part
!> \param a_first_col ...
!> \param a_first_row ...
!> \param b_first_col : the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
!> \author Samuel Andermatt
!> \note
!>      C should have no overlap with A, B
! **************************************************************************************************
   SUBROUTINE cp_complex_fm_gemm(transa, transb, m, n, k, alpha, A_re, A_im, B_re, B_im, beta, &
                                 C_re, C_im, a_first_col, a_first_row, b_first_col, b_first_row, c_first_col, &
                                 c_first_row)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      INTEGER, INTENT(IN)                                :: m, n, k
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), INTENT(IN)                       :: A_re, A_im, B_re, B_im
      REAL(KIND=dp), INTENT(IN)                          :: beta
      TYPE(cp_fm_type), INTENT(IN)                       :: C_re, C_im
      INTEGER, INTENT(IN), OPTIONAL                      :: a_first_col, a_first_row, b_first_col, &
                                                            b_first_row, c_first_col, c_first_row

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_complex_fm_gemm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_gemm(transa, transb, m, n, k, alpha, A_re, B_re, beta, C_re, &
                      a_first_col=a_first_col, &
                      a_first_row=a_first_row, &
                      b_first_col=b_first_col, &
                      b_first_row=b_first_row, &
                      c_first_col=c_first_col, &
                      c_first_row=c_first_row)
      CALL cp_fm_gemm(transa, transb, m, n, k, -alpha, A_im, B_im, 1.0_dp, C_re, &
                      a_first_col=a_first_col, &
                      a_first_row=a_first_row, &
                      b_first_col=b_first_col, &
                      b_first_row=b_first_row, &
                      c_first_col=c_first_col, &
                      c_first_row=c_first_row)
      CALL cp_fm_gemm(transa, transb, m, n, k, alpha, A_re, B_im, beta, C_im, &
                      a_first_col=a_first_col, &
                      a_first_row=a_first_row, &
                      b_first_col=b_first_col, &
                      b_first_row=b_first_row, &
                      c_first_col=c_first_col, &
                      c_first_row=c_first_row)
      CALL cp_fm_gemm(transa, transb, m, n, k, alpha, A_im, B_re, 1.0_dp, C_im, &
                      a_first_col=a_first_col, &
                      a_first_row=a_first_row, &
                      b_first_col=b_first_col, &
                      b_first_row=b_first_row, &
                      c_first_col=c_first_col, &
                      c_first_row=c_first_row)

      CALL timestop(handle)

   END SUBROUTINE cp_complex_fm_gemm

! **************************************************************************************************
!> \brief inverts a matrix using LU decomposition
!>        the input matrix will be overwritten
!> \param matrix   : input a general square non-singular matrix, outputs its inverse
!> \param info_out : optional, if present outputs the info from (p)zgetri
!> \author Lianheng Tong
! **************************************************************************************************
   SUBROUTINE cp_fm_lu_invert(matrix, info_out)
      TYPE(cp_fm_type), INTENT(IN)          :: matrix
      INTEGER, INTENT(OUT), OPTIONAL           :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_lu_invert'

      INTEGER :: nrows_global, handle, info, lwork
      INTEGER, DIMENSION(:), ALLOCATABLE       :: ipivot
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat
      REAL(KIND=sp), DIMENSION(:, :), POINTER  :: mat_sp
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
      REAL(KIND=sp), DIMENSION(:), ALLOCATABLE :: work_sp
#if defined(__SCALAPACK)
      INTEGER                                  :: liwork
      INTEGER, DIMENSION(9)                    :: desca
      INTEGER, DIMENSION(:), ALLOCATABLE       :: iwork
#else
      INTEGER                                  :: lda
#endif

      CALL timeset(routineN, handle)

      mat => matrix%local_data
      mat_sp => matrix%local_data_sp
      nrows_global = matrix%matrix_struct%nrow_global
      CPASSERT(nrows_global .EQ. matrix%matrix_struct%ncol_global)
      ALLOCATE (ipivot(nrows_global))
      ! do LU decomposition
#if defined(__SCALAPACK)
      desca = matrix%matrix_struct%descriptor
      IF (matrix%use_sp) THEN
         CALL psgetrf(nrows_global, nrows_global, &
                      mat_sp, 1, 1, desca, ipivot, info)
      ELSE
         CALL pdgetrf(nrows_global, nrows_global, &
                      mat, 1, 1, desca, ipivot, info)
      END IF
#else
      lda = SIZE(mat, 1)
      IF (matrix%use_sp) THEN
         CALL sgetrf(nrows_global, nrows_global, &
                     mat_sp, lda, ipivot, info)
      ELSE
         CALL dgetrf(nrows_global, nrows_global, &
                     mat, lda, ipivot, info)
      END IF
#endif
      IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, "LU decomposition has failed")
      END IF
      ! do inversion
      IF (matrix%use_sp) THEN
         ALLOCATE (work(1))
      ELSE
         ALLOCATE (work_sp(1))
      END IF
#if defined(__SCALAPACK)
      ALLOCATE (iwork(1))
      IF (matrix%use_sp) THEN
         CALL psgetri(nrows_global, mat_sp, 1, 1, desca, &
                      ipivot, work_sp, -1, iwork, -1, info)
         lwork = INT(work_sp(1))
         DEALLOCATE (work_sp)
         ALLOCATE (work_sp(lwork))
      ELSE
         CALL pdgetri(nrows_global, mat, 1, 1, desca, &
                      ipivot, work, -1, iwork, -1, info)
         lwork = INT(work(1))
         DEALLOCATE (work)
         ALLOCATE (work(lwork))
      END IF
      liwork = INT(iwork(1))
      DEALLOCATE (iwork)
      ALLOCATE (iwork(liwork))
      IF (matrix%use_sp) THEN
         CALL psgetri(nrows_global, mat_sp, 1, 1, desca, &
                      ipivot, work_sp, lwork, iwork, liwork, info)
      ELSE
         CALL pdgetri(nrows_global, mat, 1, 1, desca, &
                      ipivot, work, lwork, iwork, liwork, info)
      END IF
      DEALLOCATE (iwork)
#else
      IF (matrix%use_sp) THEN
         CALL sgetri(nrows_global, mat_sp, lda, &
                     ipivot, work_sp, -1, info)
         lwork = INT(work_sp(1))
         DEALLOCATE (work_sp)
         ALLOCATE (work_sp(lwork))
         CALL sgetri(nrows_global, mat_sp, lda, &
                     ipivot, work_sp, lwork, info)
      ELSE
         CALL dgetri(nrows_global, mat, lda, &
                     ipivot, work, -1, info)
         lwork = INT(work(1))
         DEALLOCATE (work)
         ALLOCATE (work(lwork))
         CALL dgetri(nrows_global, mat, lda, &
                     ipivot, work, lwork, info)
      END IF
#endif
      IF (matrix%use_sp) THEN
         DEALLOCATE (work_sp)
      ELSE
         DEALLOCATE (work)
      END IF
      DEALLOCATE (ipivot)

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CALL cp_abort(__LOCATION__, "LU inversion has failed")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_lu_invert

! **************************************************************************************************
!> \brief norm of matrix using (p)dlange
!> \param matrix   : input a general matrix
!> \param mode     : 'M' max abs element value,
!>                   '1' or 'O' one norm, i.e. maximum column sum
!>                   'I' infinity norm, i.e. maximum row sum
!>                   'F' or 'E' Frobenius norm, i.e. sqrt of sum of all squares of elements
!> \return : the norm according to mode
!> \author Lianheng Tong
! **************************************************************************************************
   FUNCTION cp_fm_norm(matrix, mode) RESULT(res)
      TYPE(cp_fm_type), INTENT(IN) :: matrix
      CHARACTER, INTENT(IN) :: mode
      REAL(KIND=dp) :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_norm'

      INTEGER :: nrows, ncols, handle, lwork, nrows_local, ncols_local
      REAL(KIND=sp) :: res_sp
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: aa
      REAL(KIND=sp), DIMENSION(:, :), POINTER :: aa_sp
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
      REAL(KIND=sp), DIMENSION(:), ALLOCATABLE :: work_sp
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: desca
#else
      INTEGER :: lda
#endif

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=matrix, &
                          nrow_global=nrows, &
                          ncol_global=ncols, &
                          nrow_local=nrows_local, &
                          ncol_local=ncols_local)
      aa => matrix%local_data
      aa_sp => matrix%local_data_sp

#if defined(__SCALAPACK)
      desca = matrix%matrix_struct%descriptor
      SELECT CASE (mode)
      CASE ('M', 'm')
         lwork = 1
      CASE ('1', 'O', 'o')
         lwork = ncols_local
      CASE ('I', 'i')
         lwork = nrows_local
      CASE ('F', 'f', 'E', 'e')
         lwork = 1
      CASE DEFAULT
         CPABORT("mode input is not valid")
      END SELECT
      IF (matrix%use_sp) THEN
         ALLOCATE (work_sp(lwork))
         res_sp = pslange(mode, nrows, ncols, aa_sp, 1, 1, desca, work_sp)
         DEALLOCATE (work_sp)
         res = REAL(res_sp, KIND=dp)
      ELSE
         ALLOCATE (work(lwork))
         res = pdlange(mode, nrows, ncols, aa, 1, 1, desca, work)
         DEALLOCATE (work)
      END IF
#else
      SELECT CASE (mode)
      CASE ('M', 'm')
         lwork = 1
      CASE ('1', 'O', 'o')
         lwork = 1
      CASE ('I', 'i')
         lwork = nrows
      CASE ('F', 'f', 'E', 'e')
         lwork = 1
      CASE DEFAULT
         CPABORT("mode input is not valid")
      END SELECT
      IF (matrix%use_sp) THEN
         ALLOCATE (work_sp(lwork))
         lda = SIZE(aa_sp, 1)
         res_sp = slange(mode, nrows, ncols, aa_sp, lda, work_sp)
         DEALLOCATE (work_sp)
         res = REAL(res_sp, KIND=dp)
      ELSE
         ALLOCATE (work(lwork))
         lda = SIZE(aa, 1)
         res = dlange(mode, nrows, ncols, aa, lda, work)
         DEALLOCATE (work)
      END IF
#endif

      CALL timestop(handle)

   END FUNCTION cp_fm_norm

! **************************************************************************************************
!> \brief trace of a matrix using pdlatra
!> \param matrix   : input a square matrix
!> \return : the trace
!> \author Lianheng Tong
! **************************************************************************************************
   FUNCTION cp_fm_latra(matrix) RESULT(res)
      TYPE(cp_fm_type), INTENT(IN) :: matrix
      REAL(KIND=dp) :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_latra'

      INTEGER :: nrows, ncols, handle
      REAL(KIND=sp) :: res_sp
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: aa
      REAL(KIND=sp), DIMENSION(:, :), POINTER :: aa_sp
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: desca
#else
      INTEGER :: ii
#endif

      CALL timeset(routineN, handle)

      nrows = matrix%matrix_struct%nrow_global
      ncols = matrix%matrix_struct%ncol_global
      CPASSERT(nrows .EQ. ncols)
      aa => matrix%local_data
      aa_sp => matrix%local_data_sp

#if defined(__SCALAPACK)
      desca = matrix%matrix_struct%descriptor
      IF (matrix%use_sp) THEN
         res_sp = pslatra(nrows, aa_sp, 1, 1, desca)
         res = REAL(res_sp, KIND=dp)
      ELSE
         res = pdlatra(nrows, aa, 1, 1, desca)
      END IF
#else
      IF (matrix%use_sp) THEN
         res_sp = 0.0_sp
         DO ii = 1, nrows
            res_sp = res_sp + aa_sp(ii, ii)
         END DO
         res = REAL(res_sp, KIND=dp)
      ELSE
         res = 0.0_dp
         DO ii = 1, nrows
            res = res + aa(ii, ii)
         END DO
      END IF
#endif

      CALL timestop(handle)

   END FUNCTION cp_fm_latra

! **************************************************************************************************
!> \brief compute a QR factorization with column pivoting of a M-by-N distributed matrix
!>        sub( A ) = A(IA:IA+M-1,JA:JA+N-1)
!> \param matrix   : input M-by-N distributed matrix sub( A ) which is to be factored
!> \param tau      :  scalar factors TAU of the elementary reflectors. TAU is tied to the distributed matrix A
!> \param nrow ...
!> \param ncol ...
!> \param first_row ...
!> \param first_col ...
!> \author MI
! **************************************************************************************************
   SUBROUTINE cp_fm_pdgeqpf(matrix, tau, nrow, ncol, first_row, first_col)

      TYPE(cp_fm_type), INTENT(IN)                    :: matrix
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tau
      INTEGER, INTENT(IN)                                :: nrow, ncol, first_row, first_col

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_pdgeqpf'

      INTEGER                                            :: handle
      INTEGER                                            :: info, lwork
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a
      REAL(KIND=dp), DIMENSION(:), POINTER               :: work
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: descc
#else
      INTEGER :: ii
#endif

      CALL timeset(routineN, handle)

      a => matrix%local_data
      lwork = -1
      ALLOCATE (work(2*nrow))
      ALLOCATE (ipiv(ncol))
      info = 0
      tau = 0

#if defined(__SCALAPACK)
      descc(:) = matrix%matrix_struct%descriptor(:)
      ! Call SCALAPACK routine to get optimal work dimension
      CALL pdgeqpf(nrow, ncol, a, first_row, first_col, descc, ipiv, tau, work, lwork, info)
      lwork = INT(work(1))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))
      tau = 0.0_dp
      ipiv = 0

      ! Call SCALAPACK routine to get QR decomposition of CTs
      CALL pdgeqpf(nrow, ncol, a, first_row, first_col, descc, ipiv, tau, work, lwork, info)
      CPASSERT(INFO == 0)
#else
      ! NIY
      ii = 1
      tau(first_row) = 0
      tau(first_col) = 0
#endif

      DEALLOCATE (work)
      DEALLOCATE (ipiv)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_pdgeqpf

! **************************************************************************************************
!> \brief generates an M-by-N real distributed matrix Q denoting A(IA:IA+M-1,JA:JA+N-1)
!>         with orthonormal columns, which is defined as the first N columns of a product of K
!>         elementary reflectors of order M
!> \param matrix : On entry, the j-th column must contain the vector which defines the elementary reflector
!>                  as returned from PDGEQRF
!>                 On exit it contains  the M-by-N distributed matrix Q
!> \param tau :   contains the scalar factors TAU of elementary reflectors  as returned by PDGEQRF
!> \param nrow ...
!> \param first_row ...
!> \param first_col ...
!> \author MI
! **************************************************************************************************
   SUBROUTINE cp_fm_pdorgqr(matrix, tau, nrow, first_row, first_col)

      TYPE(cp_fm_type), INTENT(IN)                    :: matrix
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tau
      INTEGER, INTENT(IN)                                :: nrow, first_row, first_col

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_pdorgqr'

      INTEGER                                            :: handle
      INTEGER                                            :: info, lwork
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a
      REAL(KIND=dp), DIMENSION(:), POINTER               :: work
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: descc
#else
      INTEGER :: ii
#endif

      CALL timeset(routineN, handle)

      a => matrix%local_data
      lwork = -1
      ALLOCATE (work(2*nrow))
      info = 0

#if defined(__SCALAPACK)
      descc(:) = matrix%matrix_struct%descriptor(:)

      CALL pdorgqr(nrow, nrow, nrow, a, first_row, first_col, descc, tau, work, lwork, info)
      CPASSERT(info == 0)
      lwork = INT(work(1))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))

      ! Call SCALAPACK routine to get Q
      CALL pdorgqr(nrow, nrow, nrow, a, first_row, first_col, descc, tau, work, lwork, info)
      CPASSERT(INFO == 0)
#else
      ! NIY
      ii = 1
      tau(first_row) = 0
      tau(first_col) = 0
#endif

      DEALLOCATE (work)
      CALL timestop(handle)

   END SUBROUTINE cp_fm_pdorgqr

! **************************************************************************************************
!> \brief Applies a planar rotation defined by cs and sn to the i'th and j'th rows.
!> \param cs cosine of the rotation angle
!> \param sn sinus of the rotation angle
!> \param irow ...
!> \param jrow ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cp_fm_rot_rows(matrix, irow, jrow, cs, sn)
      TYPE(cp_fm_type), INTENT(IN)             :: matrix
      INTEGER, INTENT(IN)                      :: irow, jrow
      REAL(dp), INTENT(IN)                     :: cs, sn

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_rot_rows'
      INTEGER                                  :: handle, nrow, ncol

#if defined(__SCALAPACK)
      INTEGER                                  :: info, lwork
      INTEGER, DIMENSION(9)                    :: desc
      REAL(dp), DIMENSION(:), ALLOCATABLE      :: work
#endif

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(matrix, nrow_global=nrow, ncol_global=ncol)

#if defined(__SCALAPACK)
      lwork = 2*ncol + 1
      ALLOCATE (work(lwork))
      desc(:) = matrix%matrix_struct%descriptor(:)
      CALL pdrot(ncol, &
                 matrix%local_data(1, 1), irow, 1, desc, ncol, &
                 matrix%local_data(1, 1), jrow, 1, desc, ncol, &
                 cs, sn, work, lwork, info)
      CPASSERT(info == 0)
      DEALLOCATE (work)
#else
      CALL drot(ncol, matrix%local_data(irow, 1), ncol, matrix%local_data(jrow, 1), ncol, cs, sn)
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_rot_rows

! **************************************************************************************************
!> \brief Applies a planar rotation defined by cs and sn to the i'th and j'th columnns.
!> \param cs cosine of the rotation angle
!> \param sn sinus of the rotation angle
!> \param icol ...
!> \param jcol ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cp_fm_rot_cols(matrix, icol, jcol, cs, sn)
      TYPE(cp_fm_type), INTENT(IN)             :: matrix
      INTEGER, INTENT(IN)                      :: icol, jcol
      REAL(dp), INTENT(IN)                     :: cs, sn

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_rot_cols'
      INTEGER                                  :: handle, nrow, ncol

#if defined(__SCALAPACK)
      INTEGER                                  :: info, lwork
      INTEGER, DIMENSION(9)                    :: desc
      REAL(dp), DIMENSION(:), ALLOCATABLE      :: work
#endif

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(matrix, nrow_global=nrow, ncol_global=ncol)

#if defined(__SCALAPACK)
      lwork = 2*nrow + 1
      ALLOCATE (work(lwork))
      desc(:) = matrix%matrix_struct%descriptor(:)
      CALL pdrot(nrow, &
                 matrix%local_data(1, 1), 1, icol, desc, 1, &
                 matrix%local_data(1, 1), 1, jcol, desc, 1, &
                 cs, sn, work, lwork, info)
      CPASSERT(info == 0)
      DEALLOCATE (work)
#else
      CALL drot(nrow, matrix%local_data(1, icol), 1, matrix%local_data(1, jcol), 1, cs, sn)
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_rot_cols

! **************************************************************************************************
!> \brief Orthonormalizes selected rows and columns of a full matrix, matrix_a
!> \param matrix_a ...
!> \param B ...
!> \param nrows number of rows of matrix_a, optional, defaults to size(matrix_a,1)
!> \param ncols number of columns of matrix_a, optional, defaults to size(matrix_a, 2)
!> \param start_row starting index of rows, optional, defaults to 1
!> \param start_col starting index of columns, optional, defaults to 1
!> \param do_norm ...
!> \param do_print ...
! **************************************************************************************************
   SUBROUTINE cp_fm_Gram_Schmidt_orthonorm(matrix_a, B, nrows, ncols, start_row, start_col, &
                                           do_norm, do_print)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_a
      REAL(kind=dp), DIMENSION(:, :), INTENT(OUT)        :: B
      INTEGER, INTENT(IN), OPTIONAL                      :: nrows, ncols, start_row, start_col
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_norm, do_print

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_Gram_Schmidt_orthonorm'

      INTEGER :: end_col_global, end_col_local, end_row_global, end_row_local, handle, i, j, &
                 j_col, ncol_global, ncol_local, nrow_global, nrow_local, start_col_global, &
                 start_col_local, start_row_global, start_row_local, this_col, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_norm, my_do_print
      REAL(KIND=dp)                                      :: norm
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: a

      CALL timeset(routineN, handle)

      my_do_norm = .TRUE.
      IF (PRESENT(do_norm)) my_do_norm = do_norm

      my_do_print = .FALSE.
      IF (PRESENT(do_print) .AND. (my_do_norm)) my_do_print = do_print

      unit_nr = -1
      IF (my_do_print) THEN
         unit_nr = cp_logger_get_default_unit_nr()
         IF (unit_nr < 1) my_do_print = .FALSE.
      END IF

      IF (SIZE(B) /= 0) THEN
         IF (PRESENT(nrows)) THEN
            nrow_global = nrows
         ELSE
            nrow_global = SIZE(B, 1)
         END IF

         IF (PRESENT(ncols)) THEN
            ncol_global = ncols
         ELSE
            ncol_global = SIZE(B, 2)
         END IF

         IF (PRESENT(start_row)) THEN
            start_row_global = start_row
         ELSE
            start_row_global = 1
         END IF

         IF (PRESENT(start_col)) THEN
            start_col_global = start_col
         ELSE
            start_col_global = 1
         END IF

         end_row_global = start_row_global + nrow_global - 1
         end_col_global = start_col_global + ncol_global - 1

         CALL cp_fm_get_info(matrix=matrix_a, &
                             nrow_global=nrow_global, ncol_global=ncol_global, &
                             nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices)
         IF (end_row_global > nrow_global) THEN
            end_row_global = nrow_global
         END IF
         IF (end_col_global > ncol_global) THEN
            end_col_global = ncol_global
         END IF

         ! find out row/column indices of locally stored matrix elements that
         ! needs to be copied.
         ! Arrays row_indices and col_indices are assumed to be sorted in
         ! ascending order
         DO start_row_local = 1, nrow_local
            IF (row_indices(start_row_local) >= start_row_global) EXIT
         END DO

         DO end_row_local = start_row_local, nrow_local
            IF (row_indices(end_row_local) > end_row_global) EXIT
         END DO
         end_row_local = end_row_local - 1

         DO start_col_local = 1, ncol_local
            IF (col_indices(start_col_local) >= start_col_global) EXIT
         END DO

         DO end_col_local = start_col_local, ncol_local
            IF (col_indices(end_col_local) > end_col_global) EXIT
         END DO
         end_col_local = end_col_local - 1

         a => matrix_a%local_data

         this_col = col_indices(start_col_local) - start_col_global + 1

         B(:, this_col) = a(:, start_col_local)

         IF (my_do_norm) THEN
            norm = SQRT(accurate_dot_product(B(:, this_col), B(:, this_col)))
            B(:, this_col) = B(:, this_col)/norm
            IF (my_do_print) WRITE (unit_nr, '(I3,F8.3)') this_col, norm
         END IF

         DO i = start_col_local + 1, end_col_local
            this_col = col_indices(i) - start_col_global + 1
            B(:, this_col) = a(:, i)
            DO j = start_col_local, i - 1
               j_col = col_indices(j) - start_col_global + 1
               B(:, this_col) = B(:, this_col) - &
                                accurate_dot_product(B(:, j_col), B(:, this_col))* &
                                B(:, j_col)/accurate_dot_product(B(:, j_col), B(:, j_col))
            END DO

            IF (my_do_norm) THEN
               norm = SQRT(accurate_dot_product(B(:, this_col), B(:, this_col)))
               B(:, this_col) = B(:, this_col)/norm
               IF (my_do_print) WRITE (unit_nr, '(I3,F8.3)') this_col, norm
            END IF

         END DO
         CALL matrix_a%matrix_struct%para_env%sum(B)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_Gram_Schmidt_orthonorm

! **************************************************************************************************
!> \brief Cholesky decomposition
!> \param fm_matrix ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE cp_fm_potrf(fm_matrix, n)
      TYPE(cp_fm_type)                         :: fm_matrix
      INTEGER, INTENT(in)                      :: n

      INTEGER                                  :: info
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca
#endif

      a => fm_matrix%local_data
      a_sp => fm_matrix%local_data_sp
#if defined(__SCALAPACK)
      desca(:) = fm_matrix%matrix_struct%descriptor(:)
      IF (fm_matrix%use_sp) THEN
         CALL pspotrf('U', n, a_sp(1, 1), 1, 1, desca, info)
      ELSE
         CALL pdpotrf('U', n, a(1, 1), 1, 1, desca, info)
      END IF
#else
      IF (fm_matrix%use_sp) THEN
         CALL spotrf('U', n, a_sp(1, 1), SIZE(a_sp, 1), info)
      ELSE
         CALL dpotrf('U', n, a(1, 1), SIZE(a, 1), info)
      END IF
#endif
      IF (info /= 0) &
         CPABORT("Cholesky decomposition failed. Matrix ill conditioned ?")

   END SUBROUTINE cp_fm_potrf

! **************************************************************************************************
!> \brief Invert trianguar matrix
!> \param fm_matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert
! **************************************************************************************************
   SUBROUTINE cp_fm_potri(fm_matrix, n)
      TYPE(cp_fm_type)                          :: fm_matrix
      INTEGER, INTENT(in)                       :: n

      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER   :: a_sp
      INTEGER                                   :: info
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                     :: desca
#endif

      a => fm_matrix%local_data
      a_sp => fm_matrix%local_data_sp
#if defined(__SCALAPACK)
      desca(:) = fm_matrix%matrix_struct%descriptor(:)
      IF (fm_matrix%use_sp) THEN
         CALL pspotri('U', n, a_sp(1, 1), 1, 1, desca, info)
      ELSE
         CALL pdpotri('U', n, a(1, 1), 1, 1, desca, info)
      END IF
#else
      IF (fm_matrix%use_sp) THEN
         CALL spotri('U', n, a_sp(1, 1), SIZE(a_sp, 1), info)
      ELSE
         CALL dpotri('U', n, a(1, 1), SIZE(a, 1), info)
      END IF
#endif
      CPASSERT(info == 0)
   END SUBROUTINE cp_fm_potri

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix ...
!> \param neig ...
!> \param fm_matrixb ...
!> \param fm_matrixout ...
!> \param op ...
!> \param pos ...
!> \param transa ...
! **************************************************************************************************
   SUBROUTINE cp_fm_cholesky_restore(fm_matrix, neig, fm_matrixb, fm_matrixout, op, pos, transa)
      TYPE(cp_fm_type)                               :: fm_matrix
      TYPE(cp_fm_type)                               :: fm_matrixb
      TYPE(cp_fm_type)                               :: fm_matrixout
      INTEGER, INTENT(IN)                            :: neig
      CHARACTER(LEN=*), INTENT(IN)                   :: op
      CHARACTER(LEN=*), INTENT(IN)                   :: pos
      CHARACTER(LEN=*), INTENT(IN)                   :: transa

      REAL(KIND=dp), DIMENSION(:, :), POINTER        :: a, b, outm
      REAL(KIND=sp), DIMENSION(:, :), POINTER        :: a_sp, b_sp, outm_sp
      INTEGER                                        :: n, itype
      REAL(KIND=dp)                                  :: alpha
#if defined(__SCALAPACK)
      INTEGER                                        :: i
      INTEGER, DIMENSION(9)                          :: desca, descb, descout
#endif

      ! notice b is the cholesky guy
      a => fm_matrix%local_data
      b => fm_matrixb%local_data
      outm => fm_matrixout%local_data
      a_sp => fm_matrix%local_data_sp
      b_sp => fm_matrixb%local_data_sp
      outm_sp => fm_matrixout%local_data_sp

      n = fm_matrix%matrix_struct%nrow_global
      itype = 1

#if defined(__SCALAPACK)
      desca(:) = fm_matrix%matrix_struct%descriptor(:)
      descb(:) = fm_matrixb%matrix_struct%descriptor(:)
      descout(:) = fm_matrixout%matrix_struct%descriptor(:)
      alpha = 1.0_dp
      DO i = 1, neig
         IF (fm_matrix%use_sp) THEN
            CALL pscopy(n, a_sp(1, 1), 1, i, desca, 1, outm_sp(1, 1), 1, i, descout, 1)
         ELSE
            CALL pdcopy(n, a(1, 1), 1, i, desca, 1, outm(1, 1), 1, i, descout, 1)
         END IF
      END DO
      IF (op .EQ. "SOLVE") THEN
         IF (fm_matrix%use_sp) THEN
            CALL pstrsm(pos, 'U', transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), 1, 1, descb, &
                        outm_sp(1, 1), 1, 1, descout)
         ELSE
            CALL pdtrsm(pos, 'U', transa, 'N', n, neig, alpha, b(1, 1), 1, 1, descb, outm(1, 1), 1, 1, descout)
         END IF
      ELSE
         IF (fm_matrix%use_sp) THEN
            CALL pstrmm(pos, 'U', transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), 1, 1, descb, &
                        outm_sp(1, 1), 1, 1, descout)
         ELSE
            CALL pdtrmm(pos, 'U', transa, 'N', n, neig, alpha, b(1, 1), 1, 1, descb, outm(1, 1), 1, 1, descout)
         END IF
      END IF
#else
      alpha = 1.0_dp
      IF (fm_matrix%use_sp) THEN
         CALL scopy(neig*n, a_sp(1, 1), 1, outm_sp(1, 1), 1)
      ELSE
         CALL dcopy(neig*n, a(1, 1), 1, outm(1, 1), 1)
      END IF
      IF (op .EQ. "SOLVE") THEN
         IF (fm_matrix%use_sp) THEN
            CALL strsm(pos, 'U', transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), SIZE(b_sp, 1), outm_sp(1, 1), n)
         ELSE
            CALL dtrsm(pos, 'U', transa, 'N', n, neig, alpha, b(1, 1), SIZE(b, 1), outm(1, 1), n)
         END IF
      ELSE
         IF (fm_matrix%use_sp) THEN
            CALL strmm(pos, 'U', transa, 'N', n, neig, REAL(alpha, sp), b_sp(1, 1), n, outm_sp(1, 1), n)
         ELSE
            CALL dtrmm(pos, 'U', transa, 'N', n, neig, alpha, b(1, 1), n, outm(1, 1), n)
         END IF
      END IF
#endif

   END SUBROUTINE cp_fm_cholesky_restore

END MODULE cp_fm_basic_linalg
