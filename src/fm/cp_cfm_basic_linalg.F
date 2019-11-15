!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Basic linear algebra operations for complex full matrices.
!> \note
!>      - not all functionality implemented
!> \par History
!>      Nearly literal copy of Fawzi's routines
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE cp_cfm_basic_linalg
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_prod,&
                                              mp_sum
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_basic_linalg'

   PUBLIC :: cp_cfm_cholesky_decompose, &
             cp_cfm_column_scale, &
             cp_cfm_gemm, &
             cp_cfm_lu_decompose, &
             cp_cfm_lu_invert, &
             cp_cfm_norm, &
             cp_cfm_scale, &
             cp_cfm_scale_and_add, &
             cp_cfm_scale_and_add_fm, &
             cp_cfm_schur_product, &
             cp_cfm_solve, &
             cp_cfm_trace, &
             cp_cfm_transpose, &
             cp_cfm_triangular_invert, &
             cp_cfm_triangular_multiply, &
             cp_cfm_cholesky_invert

   REAL(kind=dp), EXTERNAL :: zlange, pzlange

   INTERFACE cp_cfm_scale
      MODULE PROCEDURE cp_cfm_dscale, cp_cfm_zscale
   END INTERFACE cp_cfm_scale

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Computes the element-wise (Schur) product of two matrices: C = A \circ B .
!> \param matrix_a the first input matrix
!> \param matrix_b the second input matrix
!> \param matrix_c matrix to store the result
! **************************************************************************************************
   SUBROUTINE cp_cfm_schur_product(matrix_a, matrix_b, matrix_c)

      TYPE(cp_cfm_type), POINTER                         :: matrix_a, matrix_b, matrix_c

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_schur_product', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a, b, c
      INTEGER                                            :: handle, icol_local, irow_local, mypcol, &
                                                            myprow, ncol_local, nrow_local

      CALL timeset(routineN, handle)

      myprow = matrix_a%matrix_struct%context%mepos(1)
      mypcol = matrix_a%matrix_struct%context%mepos(2)

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

   END SUBROUTINE cp_cfm_schur_product

! **************************************************************************************************
!> \brief Computes the element-wise (Schur) product of two matrices: C = A \circ conjg(B) .
!> \param matrix_a the first input matrix
!> \param matrix_b the second input matrix
!> \param matrix_c matrix to store the result
! **************************************************************************************************
   SUBROUTINE cp_cfm_schur_product_cc(matrix_a, matrix_b, matrix_c)

      TYPE(cp_cfm_type), POINTER                         :: matrix_a, matrix_b, matrix_c

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_schur_product_cc', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a, b, c
      INTEGER                                            :: handle, icol_local, irow_local, mypcol, &
                                                            myprow, ncol_local, nrow_local

      CALL timeset(routineN, handle)

      myprow = matrix_a%matrix_struct%context%mepos(1)
      mypcol = matrix_a%matrix_struct%context%mepos(2)

      a => matrix_a%local_data
      b => matrix_b%local_data
      c => matrix_c%local_data

      nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
      ncol_local = matrix_a%matrix_struct%ncol_locals(mypcol)

      DO icol_local = 1, ncol_local
         DO irow_local = 1, nrow_local
            c(irow_local, icol_local) = a(irow_local, icol_local)*CONJG(b(irow_local, icol_local))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_schur_product_cc

! **************************************************************************************************
!> \brief Scale and add two BLACS matrices (a = alpha*a + beta*b).
!> \param alpha ...
!> \param matrix_a ...
!> \param beta ...
!> \param matrix_b ...
!> \date    11.06.2001
!> \author  Matthias Krack
!> \version 1.0
!> \note
!>    Use explicit loops to avoid temporary arrays, as a compiler reasonably assumes that arrays
!>    matrix_a%local_data and matrix_b%local_data may overlap (they are referenced by pointers).
!>    In general case (alpha*a + beta*b) explicit loops appears to be up to two times more efficient
!>    than equivalent LAPACK calls (zscale, zaxpy). This is because using LAPACK calls implies
!>    two passes through each array, so data need to be retrieved twice if arrays are large
!>    enough to not fit into the processor's cache.
! **************************************************************************************************
   SUBROUTINE cp_cfm_scale_and_add(alpha, matrix_a, beta, matrix_b)
      COMPLEX(kind=dp), INTENT(in)                       :: alpha
      TYPE(cp_cfm_type), POINTER                         :: matrix_a
      COMPLEX(kind=dp), INTENT(in), OPTIONAL             :: beta
      TYPE(cp_cfm_type), OPTIONAL, POINTER               :: matrix_b

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_scale_and_add', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp)                                   :: my_beta
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a, b
      INTEGER                                            :: handle, icol_local, irow_local, mypcol, &
                                                            myprow, ncol_local, nrow_local

      CALL timeset(routineN, handle)

      my_beta = z_zero
      IF (PRESENT(beta)) my_beta = beta
      NULLIFY (a, b)

      CPASSERT(ASSOCIATED(matrix_a))
      CPASSERT(matrix_a%ref_count > 0)
      ! to do: use dscal,dcopy,daxp
      myprow = matrix_a%matrix_struct%context%mepos(1)
      mypcol = matrix_a%matrix_struct%context%mepos(2)

      nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
      ncol_local = matrix_a%matrix_struct%ncol_locals(mypcol)

      a => matrix_a%local_data

      IF (my_beta == z_zero) THEN

         IF (alpha == z_zero) THEN
            a(:, :) = z_zero
         ELSE IF (alpha == z_one) THEN
            CALL timestop(handle)
            RETURN
         ELSE
            a(:, :) = alpha*a(:, :)
         END IF

      ELSE
         CPASSERT(PRESENT(matrix_b))
         CPASSERT(ASSOCIATED(matrix_b))
         CPASSERT(matrix_b%ref_count > 0)
         IF (matrix_a%matrix_struct%context%group /= matrix_b%matrix_struct%context%group) &
            CPABORT("matrixes must be in the same blacs context")

         IF (cp_fm_struct_equivalent(matrix_a%matrix_struct, &
                                     matrix_b%matrix_struct)) THEN

            b => matrix_b%local_data

            IF (alpha == z_zero) THEN
               IF (my_beta == z_one) THEN
                  !a(:, :) = b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = b(irow_local, icol_local)
                     END DO
                  END DO
               ELSE
                  !a(:, :) = my_beta*b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = my_beta*b(irow_local, icol_local)
                     END DO
                  END DO
               END IF
            ELSE IF (alpha == z_one) THEN
               IF (my_beta == z_one) THEN
                  !a(:, :) = a(:, :)+b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = a(irow_local, icol_local) + b(irow_local, icol_local)
                     END DO
                  END DO
               ELSE
                  !a(:, :) = a(:, :)+my_beta*b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = a(irow_local, icol_local) + my_beta*b(irow_local, icol_local)
                     END DO
                  END DO
               END IF
            ELSE
               !a(:, :) = alpha*a(:, :)+my_beta*b(:, :)
               DO icol_local = 1, ncol_local
                  DO irow_local = 1, nrow_local
                     a(irow_local, icol_local) = alpha*a(irow_local, icol_local) + my_beta*b(irow_local, icol_local)
                  END DO
               END DO
            END IF
         ELSE
#if defined(__SCALAPACK)
            CPABORT("to do (pdscal,pdcopy,pdaxpy)")
#else
            CPABORT("")
#endif
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE cp_cfm_scale_and_add

! **************************************************************************************************
!> \brief Scale and add two BLACS matrices (a = alpha*a + beta*b).
!>        where b is a real matrix (adapted from cp_cfm_scale_and_add).
!> \param alpha ...
!> \param matrix_a ...
!> \param beta ...
!> \param matrix_b ...
!> \date    01.08.2014
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_scale_and_add_fm(alpha, matrix_a, beta, matrix_b)
      COMPLEX(kind=dp), INTENT(in)                       :: alpha
      TYPE(cp_cfm_type), POINTER                         :: matrix_a
      COMPLEX(kind=dp), INTENT(in)                       :: beta
      TYPE(cp_fm_type), POINTER                          :: matrix_b

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_scale_and_add_fm', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: handle, icol_local, irow_local, mypcol, &
                                                            myprow, ncol_local, nrow_local
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: b

      CALL timeset(routineN, handle)

      NULLIFY (a, b)

      CPASSERT(ASSOCIATED(matrix_a))
      myprow = matrix_a%matrix_struct%context%mepos(1)
      mypcol = matrix_a%matrix_struct%context%mepos(2)

      nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
      ncol_local = matrix_a%matrix_struct%ncol_locals(mypcol)

      a => matrix_a%local_data

      IF (beta == z_zero) THEN

         IF (alpha == z_zero) THEN
            a(:, :) = z_zero
         ELSE IF (alpha == z_one) THEN
            CALL timestop(handle)
            RETURN
         ELSE
            a(:, :) = alpha*a(:, :)
         END IF

      ELSE
         CPASSERT(ASSOCIATED(matrix_b))
         IF (matrix_a%matrix_struct%context%group /= matrix_b%matrix_struct%context%group) &
            CPABORT("matrices must be in the same blacs context")

         IF (cp_fm_struct_equivalent(matrix_a%matrix_struct, &
                                     matrix_b%matrix_struct)) THEN

            b => matrix_b%local_data

            IF (alpha == z_zero) THEN
               IF (beta == z_one) THEN
                  !a(:, :) = b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = b(irow_local, icol_local)
                     END DO
                  END DO
               ELSE
                  !a(:, :) = beta*b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = beta*b(irow_local, icol_local)
                     END DO
                  END DO
               END IF
            ELSE IF (alpha == z_one) THEN
               IF (beta == z_one) THEN
                  !a(:, :) = a(:, :)+b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = a(irow_local, icol_local) + b(irow_local, icol_local)
                     END DO
                  END DO
               ELSE
                  !a(:, :) = a(:, :)+beta*b(:, :)
                  DO icol_local = 1, ncol_local
                     DO irow_local = 1, nrow_local
                        a(irow_local, icol_local) = a(irow_local, icol_local) + beta*b(irow_local, icol_local)
                     END DO
                  END DO
               END IF
            ELSE
               !a(:, :) = alpha*a(:, :)+beta*b(:, :)
               DO icol_local = 1, ncol_local
                  DO irow_local = 1, nrow_local
                     a(irow_local, icol_local) = alpha*a(irow_local, icol_local) + beta*b(irow_local, icol_local)
                  END DO
               END DO
            END IF
         ELSE
#if defined(__SCALAPACK)
            CPABORT("to do (pdscal,pdcopy,pdaxpy)")
#else
            CPABORT("")
#endif
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE cp_cfm_scale_and_add_fm

! **************************************************************************************************
!> \brief Computes LU decomposition of a given matrix.
!> \param matrix_a     full matrix
!> \param determinant  determinant
!> \date    11.06.2001
!> \author  Matthias Krack
!> \version 1.0
!> \note
!>    The actual purpose right now is to efficiently compute the determinant of a given matrix.
!>    The original content of the matrix is destroyed.
! **************************************************************************************************
   SUBROUTINE cp_cfm_lu_decompose(matrix_a, determinant)
      TYPE(cp_cfm_type), POINTER                         :: matrix_a
      COMPLEX(kind=dp), INTENT(out)                      :: determinant

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_lu_decompose', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: counter, handle, info, irow, nrow_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipivot

#if defined(__SCALAPACK)
      INTEGER                                            :: icol, ncol_local, nrow_local
      INTEGER, DIMENSION(9)                              :: desca
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
#else
      INTEGER                                            :: lda
#endif

      CALL timeset(routineN, handle)

      nrow_global = matrix_a%matrix_struct%nrow_global
      a => matrix_a%local_data

      ALLOCATE (ipivot(nrow_global))
#if defined(__SCALAPACK)
      CALL cp_cfm_get_info(matrix_a, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      desca(:) = matrix_a%matrix_struct%descriptor(:)
      CALL pzgetrf(nrow_global, nrow_global, a(1, 1), 1, 1, desca, ipivot, info)

      counter = 0
      DO irow = 1, nrow_local
         IF (ipivot(irow) .NE. row_indices(irow)) counter = counter + 1
      END DO

      IF (MOD(counter, 2) == 0) THEN
         determinant = z_one
      ELSE
         determinant = -z_one
      END IF

      ! compute product of diagonal elements
      irow = 1
      icol = 1
      DO WHILE (irow <= nrow_local .AND. icol <= ncol_local)
         IF (row_indices(irow) < col_indices(icol)) THEN
            irow = irow + 1
         ELSE IF (row_indices(irow) > col_indices(icol)) THEN
            icol = icol + 1
         ELSE ! diagonal element
            determinant = determinant*a(irow, icol)
            irow = irow + 1
            icol = icol + 1
         END IF
      END DO
      CALL mp_prod(determinant, matrix_a%matrix_struct%para_env%group)
#else
      lda = SIZE(a, 1)
      CALL zgetrf(nrow_global, nrow_global, a(1, 1), lda, ipivot, info)
      counter = 0
      determinant = z_one
      DO irow = 1, nrow_global
         IF (ipivot(irow) .NE. irow) counter = counter + 1
         determinant = determinant*a(irow, irow)
      ENDDO
      IF (MOD(counter, 2) == 1) determinant = -1.0_dp*determinant
#endif

      ! info is allowed to be zero
      ! this does just signal a zero diagonal element
      DEALLOCATE (ipivot)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Performs one of the matrix-matrix operations:
!>        matrix_c = alpha * op1( matrix_a ) * op2( matrix_b ) + beta*matrix_c.
!> \param transa       form of op1( matrix_a ):
!>                     op1( matrix_a ) = matrix_a,   when transa == 'N' ,
!>                     op1( matrix_a ) = matrix_a^T, when transa == 'T' ,
!>                     op1( matrix_a ) = matrix_a^H, when transa == 'C' ,
!> \param transb       form of op2( matrix_b )
!> \param m            number of rows of the matrix op1( matrix_a )
!> \param n            number of columns of the matrix op2( matrix_b )
!> \param k            number of columns of the matrix op1( matrix_a ) as well as
!>                     number of rows of the matrix op2( matrix_b )
!> \param alpha        scale factor
!> \param matrix_a     matrix A
!> \param matrix_b     matrix B
!> \param beta         scale factor
!> \param matrix_c     matrix C
!> \param a_first_col  (optional) the first column of the matrix_a to multiply
!> \param a_first_row  (optional) the first row of the matrix_a to multiply
!> \param b_first_col  (optional) the first column of the matrix_b to multiply
!> \param b_first_row  (optional) the first row of the matrix_b to multiply
!> \param c_first_col  (optional) the first column of the matrix_c
!> \param c_first_row  (optional) the first row of the matrix_c
!> \date    07.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_gemm(transa, transb, m, n, k, alpha, matrix_a, matrix_b, beta, &
                          matrix_c, a_first_col, a_first_row, b_first_col, b_first_row, c_first_col, &
                          c_first_row)
      CHARACTER(len=1), INTENT(in)                       :: transa, transb
      INTEGER, INTENT(in)                                :: m, n, k
      COMPLEX(kind=dp), INTENT(in)                       :: alpha
      TYPE(cp_cfm_type), POINTER                         :: matrix_a, matrix_b
      COMPLEX(kind=dp), INTENT(in)                       :: beta
      TYPE(cp_cfm_type), POINTER                         :: matrix_c
      INTEGER, INTENT(in), OPTIONAL                      :: a_first_col, a_first_row, b_first_col, &
                                                            b_first_row, c_first_col, c_first_row

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_gemm', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a, b, c
      INTEGER                                            :: handle, i_a, i_b, i_c, j_a, j_b, j_c
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                              :: desca, descb, descc
#else
      INTEGER                                            :: lda, ldb, ldc
#endif

      CALL timeset(routineN, handle)
      a => matrix_a%local_data
      b => matrix_b%local_data
      c => matrix_c%local_data

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

      CALL pzgemm(transa, transb, m, n, k, alpha, a(1, 1), i_a, j_a, desca, &
                  b(1, 1), i_b, j_b, descb, beta, c(1, 1), i_c, j_c, descc)
#else
      lda = SIZE(a, 1)
      ldb = SIZE(b, 1)
      ldc = SIZE(c, 1)

      CALL zgemm(transa, transb, m, n, k, alpha, a(i_a, j_a), &
                 lda, b(i_b, j_b), ldb, beta, c(i_c, j_c), ldc)
#endif
      CALL timestop(handle)
   END SUBROUTINE cp_cfm_gemm

! **************************************************************************************************
!> \brief Scales columns of the full matrix by corresponding factors.
!> \param matrix_a matrix to scale
!> \param scaling  scale factors for every column. The actual number of scaled columns is
!>                 limited by the number of scale factors given or by the actual number of columns
!>                 whichever is smaller.
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_cfm_column_scale(matrix_a, scaling)
      TYPE(cp_cfm_type), POINTER                         :: matrix_a
      COMPLEX(kind=dp), DIMENSION(:), INTENT(in)         :: scaling

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_column_scale', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: handle, icol_local, ncol_local, &
                                                            nrow_local
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
#endif

      CALL timeset(routineN, handle)

      a => matrix_a%local_data

#if defined(__SCALAPACK)
      CALL cp_cfm_get_info(matrix_a, nrow_local=nrow_local, ncol_local=ncol_local, col_indices=col_indices)
      ncol_local = MIN(ncol_local, SIZE(scaling))

      DO icol_local = 1, ncol_local
         CALL zscal(nrow_local, scaling(col_indices(icol_local)), a(1, icol_local), 1)
      END DO
#else
      nrow_local = SIZE(a, 1)
      ncol_local = MIN(SIZE(a, 2), SIZE(scaling))

      DO icol_local = 1, ncol_local
         CALL zscal(nrow_local, scaling(icol_local), a(1, icol_local), 1)
      END DO
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_column_scale

! **************************************************************************************************
!> \brief Scales a complex matrix by a real number.
!>      matrix_a = alpha * matrix_b
!> \param alpha    scale factor
!> \param matrix_a complex matrix to scale
! **************************************************************************************************
   SUBROUTINE cp_cfm_dscale(alpha, matrix_a)
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(cp_cfm_type), POINTER                         :: matrix_a

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_dscale', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY (a)

      CPASSERT(ASSOCIATED(matrix_a))
      CPASSERT(matrix_a%ref_count > 0)

      a => matrix_a%local_data

      CALL zdscal(SIZE(a), alpha, a(1, 1), 1)

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_dscale

! **************************************************************************************************
!> \brief Scales a complex matrix by a complex number.
!>      matrix_a = alpha * matrix_b
!> \param alpha    scale factor
!> \param matrix_a complex matrix to scale
!> \note
!>      use cp_fm_set_all to zero (avoids problems with nan)
! **************************************************************************************************
   SUBROUTINE cp_cfm_zscale(alpha, matrix_a)
      COMPLEX(kind=dp), INTENT(in)                       :: alpha
      TYPE(cp_cfm_type), POINTER                         :: matrix_a

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_zscale', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY (a)

      CPASSERT(ASSOCIATED(matrix_a))
      CPASSERT(matrix_a%ref_count > 0)

      a => matrix_a%local_data

      CALL zscal(SIZE(a), alpha, a(1, 1), 1)

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_zscale

! **************************************************************************************************
!> \brief Solve the system of linear equations A*b=A_general using LU decomposition.
!>        Pay attention that both matrices are overwritten on exit and that
!>        the result is stored into the matrix 'general_a'.
!> \param matrix_a     matrix A (overwritten on exit)
!> \param general_a    (input) matrix A_general, (output) matrix B
!> \param determinant  (optional) determinant
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE cp_cfm_solve(matrix_a, general_a, determinant)
      TYPE(cp_cfm_type), POINTER                         :: matrix_a, general_a
      COMPLEX(kind=dp), OPTIONAL                         :: determinant

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_solve', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a, a_general
      INTEGER                                            :: counter, handle, info, irow, nrow_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipivot

#if defined(__SCALAPACK)
      INTEGER                                            :: icol, ncol_local, nrow_local
      INTEGER, DIMENSION(9)                              :: desca, descb
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
#else
      INTEGER                                            :: lda, ldb
#endif

      CALL timeset(routineN, handle)

      a => matrix_a%local_data
      a_general => general_a%local_data
      nrow_global = matrix_a%matrix_struct%nrow_global
      ALLOCATE (ipivot(nrow_global))

#if defined(__SCALAPACK)
      desca(:) = matrix_a%matrix_struct%descriptor(:)
      descb(:) = general_a%matrix_struct%descriptor(:)
      CALL pzgetrf(nrow_global, nrow_global, a(1, 1), 1, 1, desca, ipivot, info)
      IF (PRESENT(determinant)) THEN
         CALL cp_cfm_get_info(matrix_a, nrow_local=nrow_local, ncol_local=ncol_local, &
                              row_indices=row_indices, col_indices=col_indices)

         counter = 0
         DO irow = 1, nrow_local
            IF (ipivot(irow) .NE. row_indices(irow)) counter = counter + 1
         END DO

         IF (MOD(counter, 2) == 0) THEN
            determinant = z_one
         ELSE
            determinant = -z_one
         END IF

         ! compute product of diagonal elements
         irow = 1
         icol = 1
         DO WHILE (irow <= nrow_local .AND. icol <= ncol_local)
            IF (row_indices(irow) < col_indices(icol)) THEN
               irow = irow + 1
            ELSE IF (row_indices(irow) > col_indices(icol)) THEN
               icol = icol + 1
            ELSE ! diagonal element
               determinant = determinant*a(irow, icol)
               irow = irow + 1
               icol = icol + 1
            END IF
         END DO
         CALL mp_prod(determinant, matrix_a%matrix_struct%para_env%group)
      END IF

      CALL pzgetrs("N", nrow_global, nrow_global, a(1, 1), 1, 1, desca, &
                   ipivot, a_general(1, 1), 1, 1, descb, info)
#else
      lda = SIZE(a, 1)
      ldb = SIZE(a_general, 1)
      CALL zgetrf(nrow_global, nrow_global, a(1, 1), lda, ipivot, info)
      IF (PRESENT(determinant)) THEN
         counter = 0
         determinant = z_one
         DO irow = 1, nrow_global
            IF (ipivot(irow) .NE. irow) counter = counter + 1
            determinant = determinant*a(irow, irow)
         ENDDO
         IF (MOD(counter, 2) == 1) determinant = -1.0_dp*determinant
      END IF
      CALL zgetrs("N", nrow_global, nrow_global, a(1, 1), lda, ipivot, a_general(1, 1), ldb, info)
#endif

      ! info is allowed to be zero
      ! this does just signal a zero diagonal element
      DEALLOCATE (ipivot)
      CALL timestop(handle)

   END SUBROUTINE cp_cfm_solve

! **************************************************************************************************
!> \brief Inverts a matrix using LU decomposition. The input matrix will be overwritten.
!> \param matrix     input a general square non-singular matrix, outputs its inverse
!> \param info_out   optional, if present outputs the info from (p)zgetri
!> \author Lianheng Tong
! **************************************************************************************************
   SUBROUTINE cp_cfm_lu_invert(matrix, info_out)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      INTEGER, INTENT(out), OPTIONAL                     :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_lu_invert', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: work
      COMPLEX(kind=dp), DIMENSION(1)                     :: work1
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: mat
      INTEGER                                            :: handle, info, lwork, nrows_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipivot

#if defined(__SCALAPACK)
      INTEGER                                            :: liwork
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      INTEGER, DIMENSION(1)                              :: iwork1
      INTEGER, DIMENSION(9)                              :: desca
#else
      INTEGER                                            :: lda
#endif

      CALL timeset(routineN, handle)

      mat => matrix%local_data
      nrows_global = matrix%matrix_struct%nrow_global
      CPASSERT(nrows_global .EQ. matrix%matrix_struct%ncol_global)
      ALLOCATE (ipivot(nrows_global))

      ! do LU decomposition
#if defined(__SCALAPACK)
      desca = matrix%matrix_struct%descriptor
      CALL pzgetrf(nrows_global, nrows_global, &
                   mat(1, 1), 1, 1, desca, ipivot, info)
#else
      lda = SIZE(mat, 1)
      CALL zgetrf(nrows_global, nrows_global, &
                  mat(1, 1), lda, ipivot, info)
#endif
      IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, "LU decomposition has failed")
      END IF

      ! do inversion
#if defined(__SCALAPACK)
      CALL pzgetri(nrows_global, mat(1, 1), 1, 1, desca, &
                   ipivot, work1, -1, iwork1, -1, info)
      lwork = INT(work1(1))
      liwork = INT(iwork1(1))
      ALLOCATE (work(lwork))
      ALLOCATE (iwork(liwork))
      CALL pzgetri(nrows_global, mat(1, 1), 1, 1, desca, &
                   ipivot, work, lwork, iwork, liwork, info)
      DEALLOCATE (iwork)
#else
      CALL zgetri(nrows_global, mat(1, 1), lda, ipivot, work1, -1, info)
      lwork = INT(work1(1))
      ALLOCATE (work(lwork))
      CALL zgetri(nrows_global, mat(1, 1), lda, ipivot, work, lwork, info)
#endif
      DEALLOCATE (work)
      DEALLOCATE (ipivot)

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CALL cp_abort(__LOCATION__, "LU inversion has failed")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_lu_invert

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
!> \author Joost
! **************************************************************************************************
   SUBROUTINE cp_cfm_cholesky_decompose(matrix, n, info_out)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      INTEGER, INTENT(in), OPTIONAL                      :: n
      INTEGER, INTENT(out), OPTIONAL                     :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_cholesky_decompose', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a
      INTEGER                                            :: handle, info, my_n
#if defined(__SCALAPACK)
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

#if defined(__SCALAPACK)
      desca(:) = matrix%matrix_struct%descriptor(:)
      CALL pzpotrf('U', my_n, a(1, 1), 1, 1, desca, info)
#else
      lda = SIZE(a, 1)
      CALL zpotrf('U', my_n, a(1, 1), lda, info)
#endif

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CALL cp_abort(__LOCATION__, &
                          "Cholesky decompose failed: matrix is not positive definite  or ill-conditioned")
      END IF

      CPASSERT(info == 0)

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
      TYPE(cp_cfm_type), POINTER                 :: matrix
      INTEGER, INTENT(in), OPTIONAL              :: n
      INTEGER, INTENT(out), OPTIONAL             :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_cholesky_invert', &
                                     routineP = moduleN//':'//routineN
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER  :: aa
      INTEGER                                    :: info, handle
      INTEGER                                    :: my_n
#if defined(__SCALAPACK)
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

#if defined(__SCALAPACK)
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

! **************************************************************************************************
!> \brief Returns the trace of matrix_a^T matrix_b, i.e
!>      sum_{i,j}(matrix_a(i,j)*matrix_b(i,j)) .
!> \param matrix_a a complex matrix
!> \param matrix_b another complex matrix
!> \param trace    value of the trace operator
!> \par History
!>    * 09.2017 created [Sergey Chulkov]
!> \author Sergey Chulkov
!> \note
!>      Based on the subroutine cp_fm_trace(). Note the transposition of matrix_a!
! **************************************************************************************************
   SUBROUTINE cp_cfm_trace(matrix_a, matrix_b, trace)
      TYPE(cp_cfm_type), POINTER                         :: matrix_a, matrix_b
      COMPLEX(kind=dp), INTENT(out)                      :: trace

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_trace', routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, mypcol, myprow, &
                                                            ncol_local, npcol, nprow, nrow_local
      TYPE(cp_blacs_env_type), POINTER                   :: context

      CALL timeset(routineN, handle)

      context => matrix_a%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      group = matrix_a%matrix_struct%para_env%group

      nrow_local = MIN(matrix_a%matrix_struct%nrow_locals(myprow), matrix_b%matrix_struct%nrow_locals(myprow))
      ncol_local = MIN(matrix_a%matrix_struct%ncol_locals(mypcol), matrix_b%matrix_struct%ncol_locals(mypcol))

      ! compute an accurate dot-product
      trace = accurate_dot_product(matrix_a%local_data(1:nrow_local, 1:ncol_local), &
                                   matrix_b%local_data(1:nrow_local, 1:ncol_local))

      CALL mp_sum(trace, group)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_trace

! **************************************************************************************************
!> \brief Multiplies in place by a triangular matrix:
!>       matrix_b = alpha op(triangular_matrix) matrix_b
!>      or (if side='R')
!>       matrix_b = alpha matrix_b op(triangular_matrix)
!>      op(triangular_matrix) is:
!>       triangular_matrix (if transa="N" and invert_tr=.false.)
!>       triangular_matrix^T (if transa="T" and invert_tr=.false.)
!>       triangular_matrix^H (if transa="C" and invert_tr=.false.)
!>       triangular_matrix^(-1) (if transa="N" and invert_tr=.true.)
!>       triangular_matrix^(-T) (if transa="T" and invert_tr=.true.)
!>       triangular_matrix^(-H) (if transa="C" and invert_tr=.true.)
!> \param triangular_matrix the triangular matrix that multiplies the other
!> \param matrix_b the matrix that gets multiplied and stores the result
!> \param side on which side of matrix_b stays op(triangular_matrix)
!>        (defaults to 'L')
!> \param transa_tr ...
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
   SUBROUTINE cp_cfm_triangular_multiply(triangular_matrix, matrix_b, side, &
                                         transa_tr, invert_tr, uplo_tr, unit_diag_tr, n_rows, n_cols, &
                                         alpha)
      TYPE(cp_cfm_type), POINTER                         :: triangular_matrix, matrix_b
      CHARACTER, INTENT(in), OPTIONAL                    :: side, transa_tr
      LOGICAL, INTENT(in), OPTIONAL                      :: invert_tr
      CHARACTER, INTENT(in), OPTIONAL                    :: uplo_tr
      LOGICAL, INTENT(in), OPTIONAL                      :: unit_diag_tr
      INTEGER, INTENT(in), OPTIONAL                      :: n_rows, n_cols
      COMPLEX(kind=dp), INTENT(in), OPTIONAL             :: alpha

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_triangular_multiply', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: side_char, transa, unit_diag, uplo
      COMPLEX(kind=dp)                                   :: al
      INTEGER                                            :: handle, m, n
      LOGICAL                                            :: invert

      CALL timeset(routineN, handle)
      side_char = 'L'
      unit_diag = 'N'
      uplo = 'U'
      transa = 'N'
      invert = .FALSE.
      al = CMPLX(1.0_dp, 0.0_dp, dp)
      CALL cp_cfm_get_info(matrix_b, nrow_global=m, ncol_global=n)
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
      IF (PRESENT(transa_tr)) transa = transa_tr
      IF (PRESENT(alpha)) al = alpha
      IF (PRESENT(n_rows)) m = n_rows
      IF (PRESENT(n_cols)) n = n_cols

      IF (invert) THEN

#if defined(__SCALAPACK)
         CALL pztrsm(side_char, uplo, transa, unit_diag, m, n, al, &
                     triangular_matrix%local_data(1, 1), 1, 1, &
                     triangular_matrix%matrix_struct%descriptor, &
                     matrix_b%local_data(1, 1), 1, 1, &
                     matrix_b%matrix_struct%descriptor(1))
#else
         CALL ztrsm(side_char, uplo, transa, unit_diag, m, n, al, &
                    triangular_matrix%local_data(1, 1), &
                    SIZE(triangular_matrix%local_data, 1), &
                    matrix_b%local_data(1, 1), SIZE(matrix_b%local_data, 1))
#endif

      ELSE

#if defined(__SCALAPACK)
         CALL pztrmm(side_char, uplo, transa, unit_diag, m, n, al, &
                     triangular_matrix%local_data(1, 1), 1, 1, &
                     triangular_matrix%matrix_struct%descriptor, &
                     matrix_b%local_data(1, 1), 1, 1, &
                     matrix_b%matrix_struct%descriptor(1))
#else
         CALL ztrmm(side_char, uplo, transa, unit_diag, m, n, al, &
                    triangular_matrix%local_data(1, 1), &
                    SIZE(triangular_matrix%local_data, 1), &
                    matrix_b%local_data(1, 1), SIZE(matrix_b%local_data, 1))
#endif

      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_triangular_multiply

! **************************************************************************************************
!> \brief Inverts a triangular matrix.
!> \param matrix_a ...
!> \param uplo ...
!> \param info_out ...
!> \author MI
! **************************************************************************************************
   SUBROUTINE cp_cfm_triangular_invert(matrix_a, uplo, info_out)
      TYPE(cp_cfm_type), POINTER               :: matrix_a
      CHARACTER, INTENT(in), OPTIONAL          :: uplo
      INTEGER, INTENT(out), OPTIONAL           :: info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_triangular_invert', &
                                     routineP = moduleN//':'//routineN

      CHARACTER                                :: unit_diag, my_uplo
      INTEGER                                  :: handle, info, ncol_global
      COMPLEX(kind=dp), DIMENSION(:, :), &
         POINTER                                :: a
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca
#endif

      CALL timeset(routineN, handle)

      unit_diag = 'N'
      my_uplo = 'U'
      IF (PRESENT(uplo)) my_uplo = uplo

      ncol_global = matrix_a%matrix_struct%ncol_global

      a => matrix_a%local_data

#if defined(__SCALAPACK)
      desca(:) = matrix_a%matrix_struct%descriptor(:)
      CALL pztrtri(my_uplo, unit_diag, ncol_global, a(1, 1), 1, 1, desca, info)
#else
      CALL ztrtri(my_uplo, unit_diag, ncol_global, a(1, 1), ncol_global, info)
#endif

      IF (PRESENT(info_out)) THEN
         info_out = info
      ELSE
         IF (info /= 0) &
            CALL cp_abort(__LOCATION__, &
                          "triangular invert failed: matrix is not positive definite  or ill-conditioned")
      END IF

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_triangular_invert

! **************************************************************************************************
!> \brief Transposes a BLACS distributed complex matrix.
!> \param matrix    input matrix
!> \param trans     'T' for transpose, 'C' for Hermitian conjugate
!> \param matrixt   output matrix
!> \author Lianheng Tong
! **************************************************************************************************
   SUBROUTINE cp_cfm_transpose(matrix, trans, matrixt)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      CHARACTER, INTENT(in)                              :: trans
      TYPE(cp_cfm_type), POINTER                         :: matrixt

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_transpose', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: aa, cc
      INTEGER                                            :: handle, ncol_global, nrow_global
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                              :: desca, descc
#else
      INTEGER                                            :: ii, jj
#endif

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(ASSOCIATED(matrixt))
      nrow_global = matrix%matrix_struct%nrow_global
      ncol_global = matrix%matrix_struct%ncol_global

      CPASSERT(matrixt%matrix_struct%nrow_global == ncol_global)
      CPASSERT(matrixt%matrix_struct%ncol_global == nrow_global)

      aa => matrix%local_data
      cc => matrixt%local_data

#if defined(__SCALAPACK)
      desca = matrix%matrix_struct%descriptor
      descc = matrixt%matrix_struct%descriptor
      SELECT CASE (trans)
      CASE ('T')
         CALL pztranu(nrow_global, ncol_global, &
                      z_one, aa(1, 1), 1, 1, desca, &
                      z_zero, cc(1, 1), 1, 1, descc)
      CASE ('C')
         CALL pztranc(nrow_global, ncol_global, &
                      z_one, aa(1, 1), 1, 1, desca, &
                      z_zero, cc(1, 1), 1, 1, descc)
      CASE DEFAULT
         CPABORT("trans only accepts 'T' or 'C'")
      END SELECT
#else
      SELECT CASE (trans)
      CASE ('T')
         DO jj = 1, ncol_global
            DO ii = 1, nrow_global
               cc(ii, jj) = aa(jj, ii)
            END DO
         END DO
      CASE ('C')
         DO jj = 1, ncol_global
            DO ii = 1, nrow_global
               cc(ii, jj) = CONJG(aa(jj, ii))
            END DO
         END DO
      CASE DEFAULT
         CPABORT("trans only accepts 'T' or 'C'")
      END SELECT
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_transpose

! **************************************************************************************************
!> \brief Norm of matrix using (p)zlange.
!> \param matrix     input a general matrix
!> \param mode       'M' max abs element value,
!>                   '1' or 'O' one norm, i.e. maximum column sum,
!>                   'I' infinity norm, i.e. maximum row sum,
!>                   'F' or 'E' Frobenius norm, i.e. sqrt of sum of all squares of elements
!> \return the norm according to mode
!> \author Lianheng Tong
! **************************************************************************************************
   FUNCTION cp_cfm_norm(matrix, mode) RESULT(res)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      CHARACTER, INTENT(IN)                              :: mode
      REAL(kind=dp)                                      :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_norm', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: aa
      INTEGER                                            :: handle, lwork, ncols, ncols_local, &
                                                            nrows, nrows_local
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: work

#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                              :: desca
#else
      INTEGER                                            :: lda
#endif

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=matrix, &
                           nrow_global=nrows, &
                           ncol_global=ncols, &
                           nrow_local=nrows_local, &
                           ncol_local=ncols_local)
      aa => matrix%local_data

      SELECT CASE (mode)
      CASE ('M', 'm')
         lwork = 1
      CASE ('1', 'O', 'o')
#if defined(__SCALAPACK)
         lwork = ncols_local
#else
         lwork = 1
#endif
      CASE ('I', 'i')
#if defined(__SCALAPACK)
         lwork = nrows_local
#else
         lwork = nrows
#endif
      CASE ('F', 'f', 'E', 'e')
         lwork = 1
      CASE DEFAULT
         CPABORT("mode input is not valid")
      END SELECT

      ALLOCATE (work(lwork))

#if defined(__SCALAPACK)
      desca = matrix%matrix_struct%descriptor
      res = pzlange(mode, nrows, ncols, aa(1, 1), 1, 1, desca, work)
#else
      lda = SIZE(aa, 1)
      res = zlange(mode, nrows, ncols, aa(1, 1), lda, work)
#endif

      DEALLOCATE (work)
      CALL timestop(handle)
   END FUNCTION cp_cfm_norm
END MODULE cp_cfm_basic_linalg
