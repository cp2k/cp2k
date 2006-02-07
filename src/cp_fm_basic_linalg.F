!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_fm_basic_linalg [1.0] *
!!
!!   NAME
!!     cp_fm_basic_linalg
!!
!!   FUNCTION
!!     basic linear algebra operations for full matrixes
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 splitted out of qs_blacs [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_basic_linalg
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_basic_linalg'

  PUBLIC :: cp_fm_scale, &                 ! scale a matrix
            cp_fm_scale_and_add, &         ! scale and add two matrices
            cp_fm_column_scale, &          ! scale colummns of a matrix
            cp_fm_trace, &                 ! trace of the transpose(A)*B
            cp_fm_schur_product, &         ! schur product 
            cp_fm_transpose, &             ! transpose a matrix
            cp_fm_upper_to_full, &         ! symmetrise an upper symmetric matrix
            cp_fm_syrk, &                  ! rank k update
            cp_fm_triangular_multiply, &   ! triangular matrix multiply / solve
            cp_fm_symm, &                  ! multiply a symmetric with a non-symmetric matrix
            cp_fm_gemm, &                  ! multiply two matrices
            cp_fm_lu_decompose             ! computes determinant (and lu decomp)

!****************************************************************************
CONTAINS

!!****f* cp2k/cp_fm_basic_linalg/cp_fm_scale_and_add *
!!
!!   NAME
!!     cp_fm_scale_and_add
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_scale_and_add(alpha, matrix_a, beta, matrix_b,&
!!         error)
!!       Implicit None
!!       Type(cp_fm_type), Pointer:: matrix_a
!!       Real(KIND = dp), Intent (IN):: alpha
!!       Real(KIND = dp), Intent (IN), Optional:: beta
!!       Type(cp_fm_type), Pointer, Optional:: matrix_b
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_scale_and_add
!!
!!   FUNCTION
!!     calc A <- alpha*A + beta*B
!!     optimized for alpha == 1.0 (just add beta*B) and beta == 0.0 (just
!!     scale A)
!!
!!   ARGUMENTS
!!      TYPE(cp_fm_type), POINTER           :: matrix_a
!!      REAL(KIND = dp), INTENT(IN)                         :: alpha
!!      REAL(KIND = dp), INTENT(in), OPTIONAL               :: beta
!!      TYPE(cp_fm_type), POINTER, OPTIONAL :: matrix_b
!!      TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
!!
!!****
  SUBROUTINE cp_fm_scale_and_add(alpha,matrix_a,beta,matrix_b,error)


    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_b
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_scale_and_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, size_a
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_beta
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    IF (PRESENT(matrix_b)) THEN
       my_beta=1.0_dp
    ELSE
       my_beta=0.0_dp
    ENDIF
    IF(PRESENT(beta)) my_beta=beta
    NULLIFY(a,b)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(beta)) THEN
       CPPrecondition(PRESENT(matrix_b),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(matrix_b),cp_failure_level,routineP,error,failure)
       CPPrecondition(matrix_b%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (matrix_a%id_nr==matrix_b%id_nr) THEN
          CALL cp_assert(matrix_a%id_nr/=matrix_b%id_nr, &
                         cp_warning_level, cp_assertion_failed, &
                         fromWhere=routineP, &
                         message="Bad use of routine. Call cp_fm_scale instead: "// &
                         CPSourceFileRef, &
                         error=error, failure=failure)
          CALL cp_fm_scale(alpha+beta, matrix_a, error=error)
          CALL timestop(0.0_dp,handle)
          RETURN
       END IF
    END IF

    a => matrix_a%local_data
    size_a = SIZE(a,1)*SIZE(a,2)

    IF (alpha /= 1.0_dp) CALL dscal ( size_a, alpha, a, 1)

    IF (my_beta.NE.0.0_dp) THEN
       CALL cp_assert(matrix_a%matrix_struct%context%group==&
            matrix_b%matrix_struct%context%group,cp_failure_level,&
            cp_assertion_failed,fromWhere=routineP,&
            message="matrixes must be in the same blacs context"//&
            CPSourceFileRef,&
            error=error,failure=failure)

       IF (cp_fm_struct_equivalent(matrix_a%matrix_struct,&
            matrix_b%matrix_struct,error=error)) THEN

          b => matrix_b%local_data

          CALL daxpy ( size_a, my_beta, b, 1, a, 1 )

       ELSE
#ifdef __SCALAPACK
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="to do (pdscal,pdcopy,pdaxpy)", error=error)
#else
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
#endif
       END IF

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_scale_and_add

!!****f* cp_fm_basic_linalg/cp_fm_decompose [1.0] *
!!
!!   NAME
!!     cp_fm_lu_decompose
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_lu_decompose(matrix_a, almost_determinant)
!!       Type(cp_fm_type), Pointer:: matrix_a
!!       Real(KIND = dp), Intent (OUT):: almost_determinant
!!     End Subroutine cp_fm_lu_decompose
!!
!!   FUNCTION
!!
!!     Computes the LU-decomposition of the matrix, and the determinant of the matrix
!!     IMPORTANT : the sign of the determinant is not defined correctly yet ....
!!
!!   NOTES
!!     - matrix_a is overwritten
!!     - the sign of the determinant might be wrong
!!     - SERIOUS WARNING (KNOWN BUG) : the sign of the determinant depends on ipivot
!!     - one should be able to find out if ipivot is an even or an odd permutation...
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_lu_decompose(matrix_a,almost_determinant)
    TYPE(cp_fm_type), POINTER                :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: almost_determinant

    CHARACTER(LEN=*), PARAMETER :: routine = "SUBROUTINE cp_fm_lu_decompose"

    INTEGER                                  :: handle, i, info, lda, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    INTEGER, DIMENSION(9)                    :: desca
    REAL(KIND=dp)                            :: determinant
    REAL(KIND=dp), DIMENSION(:), POINTER     :: diag
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a

! *** locals ***

    CALL timeset("cp_fm_lu_decompose","I","",handle)

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE(ipivot(n))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    CALL pdgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    ALLOCATE(diag(n))
    diag(:)=0.0_dp
    DO i=1,n
       CALL cp_fm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    determinant=1.0_dp
    DO i=1,n
       determinant=determinant*diag(i)
    ENDDO
    DEALLOCATE(diag)
#else
    lda=SIZE(a,1)
    CALL dgetrf(n,n,a(1,1),lda,ipivot,info)
    determinant=1.0_dp
    DO i=1,n
       determinant=determinant*a(i,i)
    ENDDO
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    almost_determinant=determinant ! notice that the sign is random
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE


!!****f* cp_basic_linalg/cp_fm_gemm [1.0] *
!!
!!   NAME
!!     cp_fm_gemm
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_gemm(transa, transb, m, n, k, alpha, matrix_a,&
!!         matrix_b, beta, matrix_c, error, b_first_col)
!!       Type(cp_fm_type), Pointer:: matrix_a, matrix_b, matrix_c
!!       Character(Len=1), Intent (IN):: transa, transb
!!       Real(KIND = dp), Intent (IN):: alpha, beta
!!       Integer, Intent (IN):: k, m, n
!!       Integer, Optional, Intent (IN):: b_first_col
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_gemm
!!
!!   FUNCTION
!!     computes matrix_c = beta * matrix_c + alpha * ( matrix_a  ** transa ) * ( matrix_b ** transb )
!!
!!   NOTES
!!     matrix_c should have no overlap with matrix_a, matrix_b
!!
!!   INPUTS
!!     matrix_a :: m x k matrix ( ! for transa = 'N')
!!     matrix_b :: k x n matrix ( ! for transb = 'N')
!!     matrix_c :: m x n matrix 
!!     transa :: 'N' -> normal   'T' -> transpose
!!     alpha,beta :: can be 0.0_dp and 1.0_dp  
!!     b_first_col :: the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!!
!!   AUTHOR
!!      Matthias Krack
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)

    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, INTENT(IN), OPTIONAL            :: b_first_col, a_first_row, &
                                                b_first_row, c_first_col, &
                                                c_first_row

    INTEGER                                  :: handle, i_a, i_b, i_c, j_b, &
                                                j_c, lda, ldb, ldc
    INTEGER, DIMENSION(9)                    :: desca, descb, descc
    REAL(KIND=dp)                            :: flops
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_gemm","I","Gflops",handle)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

    IF (PRESENT(b_first_row)) THEN
      i_a = a_first_row
    ELSE
      i_a = 1
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

    flops=2.0_dp*(1E-9_dp*m*n*k)/(matrix_a%matrix_struct%context%num_pe(1)*&
         matrix_a%matrix_struct%context%num_pe(2))
    CALL pdgemm(transa,transb,m,n,k,alpha,a(1,1),i_a,1,desca,b(1,1),i_b,j_b,&
                descb,beta,c(1,1),i_c,j_c,descc)
#else

    lda = SIZE(a,1)
    ldb = SIZE(b,1)
    ldc = SIZE(c,1)

    flops=2.0_dp*(1E-9_dp*m*n*k)
    CALL dgemm(transa,transb,m,n,k,alpha,a(i_a,1),lda,b(i_b,j_b),ldb,beta,c(i_c,j_c),ldc)

#endif
    CALL timestop(flops,handle)

  END SUBROUTINE cp_fm_gemm

!!****f* cp_fm_basic_linalg/cp_fm_symm [1.0] *
!!
!!   NAME
!!     cp_fm_symm
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_symm(side, uplo, m, n, alpha, matrix_a, matrix_b,&
!!         beta, matrix_c, error)
!!       Type(cp_fm_type), Pointer:: matrix_a, matrix_b, matrix_c
!!       Character(Len=1), Intent (IN):: side, uplo
!!       Real(KIND = dp), Intent (IN):: alpha, beta
!!       Integer, Intent (IN):: m, n
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_symm
!!
!!   FUNCTION
!!     computes matrix_c = beta * matrix_c + alpha *  matrix_a  *  matrix_b 
!!     OR
!!     computes matrix_c = beta * matrix_c + alpha *  matrix_b  *  matrix_a 
!!     where matrix_a is symmetric
!!
!!   NOTES
!!     matrix_c should have no overlap with matrix_a, matrix_b
!!     all matrices in QS are upper triangular, so uplo should be 'U' always
!!     matrix_a is always an m x m matrix
!!     it is typically slower to do cp_fm_symm than cp_fm_gemm (especially in parallel easily 50 percent !)
!!
!!   INPUTS
!!     matrix_a :: m x m matrix 
!!     matrix_b :: m x n matrix 
!!     matrix_c :: m x n matrix
!!     side     :: 'L' -> matrix_a is on the left 'R' -> matrix_a is on the right
!!     alpha,beta :: can be 0.0_dp and 1.0_dp
!!
!!   AUTHOR
!!      Matthias Krack
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_symm(side,uplo,m,n,alpha,matrix_a,matrix_b,beta,matrix_c,&
       error)

    CHARACTER(LEN=1), INTENT(IN)             :: side, uplo
    INTEGER, INTENT(IN)                      :: m, n
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, lda, ldb, ldc, &
                                                mypcol, myprow, npcol, nprow
    INTEGER, DIMENSION(9)                    :: desca, descb, descc
    REAL(KIND=dp)                            :: flops
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
    TYPE(cp_blacs_env_type), POINTER         :: context

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_symm","I","Gflops",handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

    flops=(2.0E-9_dp*m*m*n)/(nprow*npcol)
#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = matrix_b%matrix_struct%descriptor(:)
    descc(:) = matrix_c%matrix_struct%descriptor(:)

    CALL pdsymm(side,uplo,m,n,alpha,a(1,1),1,1,desca,b(1,1),1,1,descb,beta,c(1,1),1,1,descc)

#else

    lda = matrix_a%matrix_struct%local_leading_dimension
    ldb = matrix_b%matrix_struct%local_leading_dimension
    ldc = matrix_c%matrix_struct%local_leading_dimension

    CALL dsymm(side,uplo,m,n,alpha,a(1,1),lda,b(1,1),ldb,beta,c(1,1),ldc)

#endif
    CALL timestop(flops,handle)

  END SUBROUTINE cp_fm_symm


!!****f* cp_fm_basic_linalg/cp_fm_syrk [1.0] *
!!
!!   NAME
!!     cp_fm_syrk
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_syrk(uplo, trans, k, alpha, matrix_a, ia, ja, beta,&
!!         matrix_c, error)
!!       Type(cp_fm_type), Pointer:: matrix_a, matrix_c
!!       Character(Len=1), Intent (IN):: trans, uplo
!!       Real(KIND = dp), Intent (IN):: alpha, beta
!!       Integer, Intent (IN):: ia, ja, k
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_syrk
!!
!!   FUNCTION
!!     performs a rank-k update of a symmetric matrix_c
!!        matrix_c = beta * matrix_c + alpha * matrix_a * transpose ( matrix_a )
!!
!!   NOTES
!!     In QS uplo should 'U' (upper part updated)
!!
!!   INPUTS
!!     uplo :: 'U'   ('L')
!!     trans :: 'N'  ('T')
!!     k :: number of cols to use in matrix_a
!!     ia,ja ::  1,1 (could be used for selecting subblock of a)
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************

  SUBROUTINE cp_fm_syrk(uplo,trans,k,alpha,matrix_a,ia,ja,beta,matrix_c,error)
    CHARACTER(LEN=1), INTENT(IN)             :: uplo, trans
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a
    INTEGER, INTENT(IN)                      :: ia, ja
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, lda, ldc, mypcol, &
                                                myprow, n, npcol, nprow
    INTEGER, DIMENSION(9)                    :: desca, descc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, c
    TYPE(cp_blacs_env_type), POINTER         :: context

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_syrk","I","",handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    n = matrix_a%matrix_struct%nrow_global

    a => matrix_a%local_data
    c => matrix_c%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descc(:) = matrix_c%matrix_struct%descriptor(:)

    CALL pdsyrk(uplo,trans,n,k,alpha,a(1,1),ia,ja,desca,beta,c(1,1),1,1,descc)

#else

    lda = SIZE(a,1)
    ldc = SIZE(c,1)

    CALL dsyrk(uplo,trans,n,k,alpha,a(ia,ja),lda,beta,c(1,1),ldc)

#endif
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_syrk


!!****f* cp_fm_basic_linalg/cp_fm_schur_product [1.0] *
!!
!!   NAME
!!     cp_fm_schur_product
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_schur_product(matrix_a, matrix_b, matrix_c, error)
!!       Type(cp_fm_type), Pointer:: matrix_a, matrix_b, matrix_c
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_schur_product
!!
!!   FUNCTION
!!     computes the schur product of two matrices
!!      c_ij = a_ij * b_ij
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_schur_product(matrix_a,matrix_b,matrix_c,error)

    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b, matrix_c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, icol_local, &
                                                irow_local, mypcol, myprow, &
                                                ncol_local, npcol, nprow, &
                                                nrow_local
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
    TYPE(cp_blacs_env_type), POINTER         :: context

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_schur_product","I","",handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

    nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix_a%matrix_struct%ncol_locals(mypcol)

    DO icol_local=1,ncol_local
       DO irow_local=1,nrow_local
          c(irow_local,icol_local) = a(irow_local,icol_local)*b(irow_local,icol_local)
       END DO
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_schur_product

! *****************************************************************************
!!****f* cp_fm_basic_linalg/cp_fm_trace [1.0] *
!!
!!   NAME
!!     cp_fm_trace
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_trace(matrix_a, matrix_b, trace, error)
!!       Type(cp_fm_type), Pointer:: matrix_a, matrix_b
!!       Real(KIND = dp), Intent (OUT):: trace
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_trace
!!
!!   FUNCTION
!!     returns the trace of matrix_a^T matrix_b, i.e 
!!     sum_{i,j}(matrix_a(i,j)*matrix_b(i,j))
!!
!!   NOTES
!!     note the transposition of matrix_a!
!!
!!   ARGUMENTS
!!     - matrix_a: a matrix
!!     - matrix_b: another matrix
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     11.06.2001 Creation (Matthias Krack)
!!     12.2002 added doc [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_trace(matrix_a,matrix_b,trace,error)

    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(OUT)               :: trace
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: group, handle, icol_local, &
                                                irow_local, mypcol, myprow, &
                                                ncol_local, npcol, nprow, &
                                                nrow_local
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    TYPE(cp_blacs_env_type), POINTER         :: context

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_trace","I","",handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    group = matrix_a%matrix_struct%para_env%group

    a => matrix_a%local_data
    b => matrix_b%local_data

    nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix_b%matrix_struct%ncol_locals(mypcol)

    trace = 0.0_dp

    DO icol_local=1,ncol_local
       DO irow_local=1,nrow_local
          trace = trace + a(irow_local,icol_local)*b(irow_local,icol_local)
       END DO
    END DO

    CALL mp_sum(trace,group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_trace

! *****************************************************************************


!!****f* cp_fm_basic_linalg/cp_fm_triangular_multiply [1.0] *
!!
!!   NAME
!!     cp_fm_triangular_multiply
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_triangular_multiply(triangular_matrix, matrix_b,&
!!         side, transpose_tr, invert_tr, uplo_tr, unit_diag_tr, n_rows,&
!!         n_cols, alpha, error)
!!       Type(cp_fm_type), Pointer:: triangular_matrix, matrix_b
!!       Character, Intent (IN), Optional:: side, uplo_tr
!!       Logical, Intent (IN), Optional:: transpose_tr, invert_tr,&
!!         unit_diag_tr
!!       Integer, Intent (IN), Optional:: n_rows, n_cols
!!       Real(KIND = dp), Intent (IN), Optional:: alpha
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_triangular_multiply
!!
!!   FUNCTION
!!     multiplies in place by a triangular matrix:
!!      matrix_b = alpha op(triangular_matrix) matrix_b
!!     or (if side='R')
!!      matrix_b = alpha matrix_b op(triangular_matrix)
!!     op(triangular_matrix) is:
!!      triangular_matrix (if transpose_tr=.false. and invert_tr=.false.)
!!      triangular_matrix^T (if transpose_tr=.true. and invert_tr=.false.)
!!      triangular_matrix^(-1) (if transpose_tr=.false. and invert_tr=.true.)
!!      triangular_matrix^(-T) (if transpose_tr=.true. and invert_tr=.true.)
!!
!!   NOTES
!!     needs an mpi env
!!
!!   INPUTS
!!     - triangular_matrix: the triangular matrix that multiplies the other
!!     - matrix_b: the matrix that gets multiplied and stores the result
!!     - side: on which side of matrix_b stays op(triangular_matrix)
!!       (defaults to 'L')
!!     - transpose_tr: if the triangular matrix should be transposed
!!       (defaults to false)
!!     - invert_tr: if the triangular matrix should be inverted
!!       (defaults to false)
!!     - uplo_tr: if triangular_matrix is stored in the upper ('U') or
!!       lower ('L') triangle (defaults to 'U')
!!     - unit_diag_tr: if the diagonal elements of triangular_matrix should
!!       be assumed to be 1 (defaults to false)
!!     - n_rows: the number of rows of the result (defaults to 
!!       size(matrix_b,1))
!!     - n_cols: the number of columns of the result (defaults to
!!       size(matrix_b,2))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_fm_triangular_multiply(triangular_matrix,matrix_b,side,&
     transpose_tr, invert_tr, uplo_tr,unit_diag_tr, n_rows, n_cols, &
     alpha,error)
    TYPE(cp_fm_type), POINTER                :: triangular_matrix, matrix_b
    CHARACTER, INTENT(in), OPTIONAL          :: side
    LOGICAL, INTENT(in), OPTIONAL            :: transpose_tr, invert_tr
    CHARACTER, INTENT(in), OPTIONAL          :: uplo_tr
    LOGICAL, INTENT(in), OPTIONAL            :: unit_diag_tr
    INTEGER, INTENT(in), OPTIONAL            :: n_rows, n_cols
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_triangular_multiply', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: side_char, transa, unit_diag, &
                                                uplo
    INTEGER                                  :: handle, m, n
    LOGICAL                                  :: failure, invert
    REAL(KIND=dp)                            :: al

  failure=.FALSE.
  
  CALL timeset(routineN//','//moduleN,"I","",handle)
  side_char='L'
  unit_diag='N'
  uplo='U'
  transa='N'
  invert=.FALSE.
  al=1.0_dp
  CALL cp_fm_get_info(matrix_b, nrow_global=m, ncol_global=n, error=error)
  IF (PRESENT(side)) side_char=side
  IF (PRESENT(invert_tr)) invert=invert_tr
  IF (PRESENT(uplo_tr)) uplo=uplo_tr
  IF (PRESENT(unit_diag_tr)) THEN
     IF (unit_diag_tr) THEN
        unit_diag='U'
     ELSE
        unit_diag='N'
     END IF
  END IF
  IF (PRESENT(transpose_tr)) THEN
     IF (transpose_tr) THEN
        transa='T'
     ELSE
        transa='N'
     END IF
  END IF
  IF (PRESENT(alpha)) al=alpha
  IF (PRESENT(n_rows)) m=n_rows
  IF (PRESENT(n_cols)) n=n_cols

  IF (invert) THEN

#if defined(__SCALAPACK)
     CALL pdtrsm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),1,1,&
          triangular_matrix%matrix_struct%descriptor,&
          matrix_b%local_data(1,1),1,1,&
          matrix_b%matrix_struct%descriptor(1))
#else     
     CALL dtrsm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),&
          SIZE(triangular_matrix%local_data,1),&
          matrix_b%local_data(1,1),SIZE(matrix_b%local_data,1))
#endif

  ELSE

#if defined(__SCALAPACK)
     CALL pdtrmm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),1,1,&
          triangular_matrix%matrix_struct%descriptor,&
          matrix_b%local_data(1,1),1,1,&
          matrix_b%matrix_struct%descriptor(1))
#else     
     CALL dtrmm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),&
          SIZE(triangular_matrix%local_data,1),&
          matrix_b%local_data(1,1),SIZE(matrix_b%local_data,1))
#endif

  END IF

  CALL timestop(0.0_dp,handle)     
  END SUBROUTINE cp_fm_triangular_multiply


!!****f* cp_fm_basic_linalg/cp_fm_scale [1.0] *
!!
!!   NAME
!!     cp_fm_scale
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_scale(alpha, matrix_a, error)
!!       Type(cp_fm_type), Pointer:: matrix_a
!!       Real(KIND = dp), Intent (IN):: alpha
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_scale
!!
!!   FUNCTION
!!     scales a matrix
!!     matrix_a = alpha * matrix_b
!!
!!   NOTES
!!     use cp_fm_set_all to zero (avoids problems with nan)
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************

  SUBROUTINE cp_fm_scale(alpha, matrix_a, error)
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, size_a
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(a)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)

    a => matrix_a%local_data
    size_a = SIZE(a,1)*SIZE(a,2)

    CALL DSCAL(size_a, alpha, a, 1)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_scale

!!****f* cp_fm_basic_linalg/cp_fm_transpose [1.0] *
!!
!!   NAME
!!     cp_fm_transpose
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_transpose(matrix, matrixt, error)
!!       Type(cp_fm_type), Pointer:: matrix, matrixt
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_transpose
!!
!!   FUNCTION
!!     transposes a matrix
!!     matrixt = matrix ^ T
!!
!!   NOTES
!!     all matrix elements are transpose (see cp_fm_upper_to_half to symmetrise a matrix)
!!     all matrix elements are transpose (see cp_fm_upper_to_half to symmetrize a matrix)
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************

  SUBROUTINE cp_fm_transpose(matrix,matrixt,error)
    TYPE(cp_fm_type), POINTER                :: matrix, matrixt
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_transpose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, mypcol, myprow, &
                                                ncol_global, npcol, nprow, &
                                                nrow_global
    INTEGER, DIMENSION(9)                    :: desca, descc
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, c
    TYPE(cp_blacs_env_type), POINTER         :: context

!   ---------------------------------------------------------------------------

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrixt),cp_failure_level,routineP,error,failure)
    IF (failure) RETURN
    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global
    CPPrecondition(nrow_global==ncol_global,cp_failure_level,routineP,error,failure)
    IF (failure) RETURN

    CALL timeset("cp_fm_transpose","I","",handle)

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)


    a => matrix%local_data
    c => matrixt%local_data

#if defined(__SCALAPACK)
    desca(:) = matrix%matrix_struct%descriptor(:)
    descc(:) = matrixt%matrix_struct%descriptor(:)
    CALL pdtran(nrow_global,ncol_global,1.0_dp,a(1,1),1,1,desca,0.0_dp,c(1,1),1,1,descc)
#else
    DO j=1,ncol_global
     DO i=1,nrow_global
        c(j,i)=a(i,j)
     ENDDO
    ENDDO
#endif
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_transpose


!!****f* cp_fm_basic_linalg/cp_fm_upper_to_full [1.0] *
!!
!!   NAME
!!     cp_fm_upper_to_full
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_upper_to_full(matrix, work, error)
!!       Type(cp_fm_type), Pointer:: matrix, work
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_upper_to_full
!!
!!   FUNCTION
!!      given an upper triangular matrix computes the corresponding full matrix
!!
!!   NOTES
!!      the lower triangular part is irrelevant
!!
!!   INPUTS
!!      matrix : the upper triangular matrix
!!      work   : a matrix of the same size as matrix
!!
!!   AUTHOR
!!      Matthias Krack
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_upper_to_full(matrix,work,error)

    TYPE(cp_fm_type), POINTER          :: matrix,work
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
               irow_local,mypcol,myprow,ncol_block,ncol_global,ncol_local,&
               npcol,nprow,nrow_block,nrow_global,nrow_local

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(KIND = dp), DIMENSION(:,:), POINTER :: a,c
    TYPE(cp_blacs_env_type), POINTER :: context
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_upper_to_full',&
         routineP=moduleN//':'//routineN 
    LOGICAL :: failure

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxl2g

#endif
!   ---------------------------------------------------------------------------

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(work),cp_failure_level,routineP,error,failure)
    IF (failure) RETURN
    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global
    CPPrecondition(nrow_global==ncol_global,cp_failure_level,routineP,error,failure)
    nrow_global = work%matrix_struct%nrow_global
    ncol_global = work%matrix_struct%ncol_global
    CPPrecondition(nrow_global==ncol_global,cp_failure_level,routineP,error,failure)
    IF (failure) RETURN

    CALL timeset("cp_fm_upper_to_full","I","",handle)

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

#if defined(__SCALAPACK)


    nrow_block = matrix%matrix_struct%nrow_block
    ncol_block = matrix%matrix_struct%ncol_block

    nrow_local = matrix%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix%matrix_struct%ncol_locals(mypcol)

    a => work%local_data
    desca(:) = work%matrix_struct%descriptor(:)
    c => matrix%local_data
    descc(:) = matrix%matrix_struct%descriptor(:)

    DO icol_local=1,ncol_local
      icol_global = indxl2g(icol_local,ncol_block,mypcol,&
           matrix%matrix_struct%first_p_pos(2),npcol)
      DO irow_local=1,nrow_local
        irow_global = indxl2g(irow_local,nrow_block,myprow,&
             matrix%matrix_struct%first_p_pos(1),nprow)
        IF (irow_global > icol_global) THEN
          c(irow_local,icol_local) = 0.0_dp
        ELSE IF (irow_global == icol_global) THEN
          c(irow_local,icol_local) = 0.5_dp*c(irow_local,icol_local)
        END IF
      END DO
    END DO

    a(:,:) = c(:,:)

    CALL pdtran(nrow_global,ncol_global,1.0_dp,a(1,1),1,1,desca,1.0_dp,c(1,1),1,1,descc)

#else

    a => matrix%local_data
    DO irow_global=1,nrow_global
     DO icol_global=irow_global+1,ncol_global
      a(icol_global,irow_global)=a(irow_global,icol_global)
     ENDDO
    ENDDO

#endif
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_upper_to_full

!!****f* cp_fm_basic_linalg/cp_fm_column_scale [1.0] *
!!
!!   NAME
!!     cp_fm_column_scale
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_column_scale(matrixa, scaling)
!!       Type(cp_fm_type), Pointer:: matrixa
!!       Real(KIND = dp), Dimension(:), Intent (IN):: scaling
!!     End Subroutine cp_fm_column_scale
!!
!!   FUNCTION
!!     scales column i of matrix a with scaling(i)
!!
!!   NOTES
!!     this is very useful as a first step in the computation of C = sum_i alpha_i A_i transpose (A_i)
!!     that is a rank-k update (cp_fm_syrk , cp_sm_plus_fm_fm_t) where every vector has a different prefactor
!!     this procedure can be up to 20 times faster than calling cp_fm_syrk n times
!!
!!   INPUTS
!!     scaling :: an array used for scaling the columns, SIZE(scaling) determines the number of columns to be scaled
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_column_scale(matrixa,scaling)
    TYPE(cp_fm_type), POINTER                :: matrixa
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: scaling

    INTEGER                                  :: i, icol_global, icol_local, &
                                                ipcol, iprow, irow_local, k, &
                                                mypcol, myprow, n, npcol, &
                                                nprow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a

!   ---------------------------------------------------------------------------

    myprow=matrixa%matrix_struct%context%mepos(1)
    mypcol=matrixa%matrix_struct%context%mepos(2)
    nprow=matrixa%matrix_struct%context%num_pe(1)
    npcol=matrixa%matrix_struct%context%num_pe(2)

    a => matrixa%local_data
    n =  SIZE(a,1)
    k =  SIZE(scaling)

#if defined(__SCALAPACK)

    DO icol_global=1,k
       CALL infog2l(1,icol_global,matrixa%matrix_struct%descriptor,&
            nprow,npcol,myprow,mypcol,&
            irow_local,icol_local,iprow,ipcol)
       IF ((ipcol == mypcol)) THEN
           CALL DSCAL(n,scaling(icol_global),a(1,icol_local),1)
       END IF
    ENDDO
#else
    DO i=1,k
      CALL DSCAL(n,scaling(i),a(1,i),1)
    ENDDO
#endif
  END SUBROUTINE cp_fm_column_scale

END MODULE cp_fm_basic_linalg
