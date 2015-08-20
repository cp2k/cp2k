!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basic linear algebra operations for full matrixes
!> \par History
!>      08.2002 splitted out of qs_blacs [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_basic_linalg
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, &
       cp_fm_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp,&
                                             sp
  USE machine,                         ONLY: m_memory
  USE mathlib,                         ONLY: invert_matrix
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_basic_linalg'

  PUBLIC :: cp_fm_scale,&  ! scale a matrix
            cp_fm_scale_and_add,&  ! scale and add two matrices
            cp_fm_column_scale,&  ! scale colummns of a matrix
            cp_fm_trace,&  ! trace of the transpose(A)*B
            cp_fm_schur_product,&  ! schur product
            cp_fm_transpose,&  ! transpose a matrix
            cp_fm_upper_to_full,&  ! symmetrise an upper symmetric matrix
            cp_fm_syrk,&  ! rank k update
            cp_fm_triangular_multiply,&  ! triangular matrix multiply / solve
            cp_fm_symm,&  ! multiply a symmetric with a non-symmetric matrix
            cp_fm_gemm,&  ! multiply two matrices
            cp_complex_fm_gemm,&  ! multiply two complex matrices, represented by non_complex fm matrices
            cp_fm_invert,&  ! computes the inverse and determinant
            cp_fm_frobenius_norm,&  ! frobenius norm
            cp_fm_triangular_invert,&  ! compute the reciprocal of a tirangular matrix
            cp_fm_qr_factorization,&  ! compute the QR factorization of a rectangular matrix
            cp_fm_solve  ! solves the equation  A*B=C A and C are input

CONTAINS

! *****************************************************************************
!> \brief calc A <- alpha*A + beta*B
!>      optimized for alpha == 1.0 (just add beta*B) and beta == 0.0 (just
!>      scale A)
!> \param alpha ...
!> \param matrix_a ...
!> \param beta ...
!> \param matrix_b ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_fm_scale_and_add(alpha,matrix_a,beta,matrix_b,error)

    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_b
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_scale_and_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, size_a
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_beta
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp, b_sp

    CALL timeset(routineN,handle)

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
          CALL timestop(handle)
          RETURN
       END IF
    END IF

    a => matrix_a%local_data
    a_sp => matrix_a%local_data_sp

    IF(matrix_a%use_sp) THEN
       size_a = SIZE(a_sp,1)*SIZE(a_sp,2)
    ELSE
       size_a = SIZE(a,1)*SIZE(a,2)
    ENDIF

    IF (alpha /= 1.0_dp) THEN
       IF(matrix_a%use_sp) THEN
          CALL sscal ( size_a, REAL(alpha,sp), a_sp, 1)
       ELSE
          CALL dscal ( size_a, alpha, a, 1)
       ENDIF
    ENDIF
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
          b_sp => matrix_b%local_data_sp

          IF(matrix_a%use_sp.AND.matrix_b%use_sp) THEN
             CALL saxpy ( size_a, REAL(my_beta,sp), b_sp, 1, a_sp, 1 )
          ELSEIF(matrix_a%use_sp.AND..NOT.matrix_b%use_sp) THEN
             CALL saxpy ( size_a, REAL(my_beta,sp), REAL(b,sp), 1, a_sp, 1 )
          ELSEIF(.NOT.matrix_a%use_sp.AND.matrix_b%use_sp) THEN
             CALL daxpy ( size_a, my_beta, REAL(b_sp,dp), 1, a, 1 )
          ELSE
             CALL daxpy ( size_a, my_beta, b, 1, a, 1 )
          ENDIF

       ELSE
#ifdef __SCALAPACK
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="to do (pdscal,pdcopy,pdaxpy)", error=error)
#else
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
#endif
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_fm_scale_and_add

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE cp_fm_lu_decompose(matrix_a,almost_determinant,correct_sign)
    TYPE(cp_fm_type), POINTER                :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: almost_determinant
    LOGICAL, INTENT(IN), OPTIONAL            :: correct_sign

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_lu_decompose', &
      routineP = moduleN//':'//routineN

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

    CALL timeset(routineN,handle)

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE (ipivot(n+matrix_a%matrix_struct%nrow_block))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    CALL pdgetrf(n,n,a(1,1),1,1,desca,ipivot,info)

    ALLOCATE (diag(n))
    diag(:)=0.0_dp
    DO i=1,n
       CALL cp_fm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    determinant=1.0_dp
    DO i=1,n
       determinant=determinant*diag(i)
    ENDDO
    DEALLOCATE (diag)
#else
    lda=SIZE(a,1)
    CALL dgetrf(n,n,a(1,1),lda,ipivot,info)
    determinant=1.0_dp
    IF(correct_sign)THEN
       DO i=1,n
          IF(ipivot(i).NE.i)THEN
             determinant=-determinant*a(i,i)
          ELSE
             determinant=determinant*a(i,i)
          END IF
       END DO
    ELSE
       DO i=1,n
          determinant=determinant*a(i,i)
       ENDDO
    END IF
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    almost_determinant=determinant ! notice that the sign is random
    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
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
!> \param error ...
!> \param b_first_col : the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!> \param a_first_row ...
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
!> \author Matthias Krack
!> \note
!>      matrix_c should have no overlap with matrix_a, matrix_b
! *****************************************************************************
  SUBROUTINE cp_fm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)

    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: b_first_col, a_first_row, &
                                                b_first_row, c_first_col, &
                                                c_first_row

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_gemm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_a, i_b, i_c, j_b, &
                                                j_c
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp, b_sp, c_sp
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb, descc
#else
    INTEGER                                  :: lda, ldb, ldc
#endif

    CALL timeset(routineN,handle)

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

    IF(matrix_a%use_sp.AND.matrix_b%use_sp.AND.matrix_c%use_sp) THEN

       CALL psgemm(transa,transb,m,n,k,REAL(alpha,sp),a_sp(1,1),i_a,1,desca,b_sp(1,1),i_b,j_b,&
                   descb,REAL(beta,sp),c_sp(1,1),i_c,j_c,descc)

    ELSEIF((.NOT.matrix_a%use_sp).AND.(.NOT.matrix_b%use_sp).AND.(.NOT.matrix_c%use_sp)) THEN

       CALL pdgemm(transa,transb,m,n,k,alpha,a(1,1),i_a,1,desca,b(1,1),i_b,j_b,&
                   descb,beta,c(1,1),i_c,j_c,descc)

    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"Mixed precision gemm NYI")
    ENDIF
#else

    IF(matrix_a%use_sp.AND.matrix_b%use_sp.AND.matrix_c%use_sp) THEN

       lda = SIZE(a_sp,1)
       ldb = SIZE(b_sp,1)
       ldc = SIZE(c_sp,1)

       CALL sgemm(transa,transb,m,n,k,REAL(alpha,sp),a_sp(i_a,1),lda,b_sp(i_b,j_b),ldb,&
                  REAL(beta,sp),c_sp(i_c,j_c),ldc)

    ELSEIF((.NOT.matrix_a%use_sp).AND.(.NOT.matrix_b%use_sp).AND.(.NOT.matrix_c%use_sp)) THEN

       lda = SIZE(a,1)
       ldb = SIZE(b,1)
       ldc = SIZE(c,1)

       CALL dgemm(transa,transb,m,n,k,alpha,a(i_a,1),lda,b(i_b,j_b),ldb,beta,c(i_c,j_c),ldc)

    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"Mixed precision gemm NYI")
    ENDIF

#endif
    CALL timestop(handle)

  END SUBROUTINE cp_fm_gemm

! *****************************************************************************
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
!> \param error ...
!> \author Matthias Krack
!> \note
!>      matrix_c should have no overlap with matrix_a, matrix_b
!>      all matrices in QS are upper triangular, so uplo should be 'U' always
!>      matrix_a is always an m x m matrix
!>      it is typically slower to do cp_fm_symm than cp_fm_gemm (especially in parallel easily 50 percent !)
! *****************************************************************************
  SUBROUTINE cp_fm_symm(side,uplo,m,n,alpha,matrix_a,matrix_b,beta,matrix_c,&
       error)

    CHARACTER(LEN=1), INTENT(IN)             :: side, uplo
    INTEGER, INTENT(IN)                      :: m, n
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_symm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb, descc
#else
    INTEGER                                  :: lda, ldb, ldc
#endif

    CALL timeset(routineN,handle)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

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
    CALL timestop(handle)

  END SUBROUTINE cp_fm_symm

! *****************************************************************************
!> \brief computes the Frobenius norm of matrix_a
!> \param matrix_a : m x n matrix
!> \param norm ...
!> \param error ...
!> \author VW
! *****************************************************************************
  SUBROUTINE cp_fm_frobenius_norm(matrix_a,norm,error)
    TYPE(cp_fm_type), POINTER                :: matrix_a
    REAL(KIND=dp), INTENT(inout)             :: norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_frobenius_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, size_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    REAL(KIND=dp), EXTERNAL                  :: DDOT
#if defined(__SCALAPACK)
    INTEGER                                  :: group
#endif

    CALL timeset(routineN,handle)

    norm = 0.0_dp
    a => matrix_a%local_data
    size_a = SIZE(a,1)*SIZE(a,2)
    norm = DDOT( size_a, a(1,1), 1, a(1,1), 1 )
#if defined(__SCALAPACK)
    group = matrix_a%matrix_struct%para_env%group
    CALL mp_sum(norm,group)
#endif
    norm = SQRT(norm)
 
    CALL timestop(handle)

 END SUBROUTINE cp_fm_frobenius_norm

! *****************************************************************************
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
!> \param error ...
!> \author Matthias Krack
!> \note
!>      In QS uplo should 'U' (upper part updated)
! *****************************************************************************
  SUBROUTINE cp_fm_syrk(uplo,trans,k,alpha,matrix_a,ia,ja,beta,matrix_c,error)
    CHARACTER(LEN=1), INTENT(IN)             :: uplo, trans
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a
    INTEGER, INTENT(IN)                      :: ia, ja
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_syrk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, c
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descc
#else
    INTEGER                                  :: lda, ldc
#endif

    CALL timeset(routineN,handle)

    n = matrix_c%matrix_struct%nrow_global

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
    CALL timestop(handle)

  END SUBROUTINE cp_fm_syrk

! *****************************************************************************
!> \brief computes the schur product of two matrices
!>       c_ij = a_ij * b_ij
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param error ...
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_fm_schur_product(matrix_a,matrix_b,matrix_c,error)

    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b, matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_schur_product', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icol_local, &
                                                irow_local, mypcol, myprow, &
                                                ncol_local, npcol, nprow, &
                                                nrow_local
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b, c
    TYPE(cp_blacs_env_type), POINTER         :: context

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)

  END SUBROUTINE cp_fm_schur_product

! *****************************************************************************
!> \brief returns the trace of matrix_a^T matrix_b, i.e
!>      sum_{i,j}(matrix_a(i,j)*matrix_b(i,j))
!> \param matrix_a a matrix
!> \param matrix_b another matrix
!> \param trace ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.06.2001 Creation (Matthias Krack)
!>      12.2002 added doc [fawzi]
!> \author Matthias Krack
!> \note
!>      note the transposition of matrix_a!
! *****************************************************************************
  SUBROUTINE cp_fm_trace(matrix_a,matrix_b,trace,error)

    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(OUT)               :: trace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_trace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, mypcol, &
                                                myprow, ncol_local, npcol, &
                                                nprow, nrow_local
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp, b_sp
    TYPE(cp_blacs_env_type), POINTER         :: context

    CALL timeset(routineN,handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    group = matrix_a%matrix_struct%para_env%group

    a => matrix_a%local_data
    b => matrix_b%local_data

    a_sp => matrix_a%local_data_sp
    b_sp => matrix_b%local_data_sp

    nrow_local = MIN(matrix_a%matrix_struct%nrow_locals(myprow),matrix_b%matrix_struct%nrow_locals(myprow))
    ncol_local = MIN(matrix_a%matrix_struct%ncol_locals(mypcol),matrix_b%matrix_struct%ncol_locals(mypcol))

    ! cries for an accurate_dot_product
    IF(matrix_a%use_sp.AND.matrix_b%use_sp) THEN
       trace = accurate_sum(REAL(a_sp(1:nrow_local,1:ncol_local) * &
                                 b_sp(1:nrow_local,1:ncol_local), dp) )
    ELSEIF(matrix_a%use_sp.AND..NOT.matrix_b%use_sp) THEN
       trace = accurate_sum(REAL(a_sp(1:nrow_local,1:ncol_local), dp) * &
                                    b(1:nrow_local,1:ncol_local) )
    ELSEIF(.NOT.matrix_a%use_sp.AND.matrix_b%use_sp) THEN
       trace = accurate_sum(        a(1:nrow_local,1:ncol_local) * &
                            REAL(b_sp(1:nrow_local,1:ncol_local), dp) )
    ELSE
       trace = accurate_sum(a(1:nrow_local,1:ncol_local) * &
                            b(1:nrow_local,1:ncol_local) )
    ENDIF

    CALL mp_sum(trace,group)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_trace

! *****************************************************************************
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs an mpi env
! *****************************************************************************
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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_triangular_multiply', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: side_char, transa, unit_diag, &
                                                uplo
    INTEGER                                  :: handle, m, n
    LOGICAL                                  :: failure, invert
    REAL(KIND=dp)                            :: al

  failure=.FALSE.

  CALL timeset(routineN,handle)
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

  CALL timestop(handle)
  END SUBROUTINE cp_fm_triangular_multiply

! *****************************************************************************
!> \brief scales a matrix
!>      matrix_a = alpha * matrix_b
!> \param alpha ...
!> \param matrix_a ...
!> \param error ...
!> \note
!>      use cp_fm_set_all to zero (avoids problems with nan)
! *****************************************************************************
  SUBROUTINE cp_fm_scale(alpha, matrix_a, error)
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, size_a
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(a)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)

    a => matrix_a%local_data
    size_a = SIZE(a,1)*SIZE(a,2)

    CALL DSCAL(size_a, alpha, a, 1)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_scale

! *****************************************************************************
!> \brief transposes a matrix
!>      matrixt = matrix ^ T
!> \param matrix ...
!> \param matrixt ...
!> \param error ...
!> \note
!>      all matrix elements are transpose (see cp_fm_upper_to_half to symmetrise a matrix)
!>      all matrix elements are transpose (see cp_fm_upper_to_half to symmetrize a matrix)
! *****************************************************************************
  SUBROUTINE cp_fm_transpose(matrix,matrixt,error)
    TYPE(cp_fm_type), POINTER                :: matrix, matrixt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_transpose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ncol_global, &
                                                nrow_global
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, c
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descc
#else
    INTEGER                                  :: i, j
#endif

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrixt),cp_failure_level,routineP,error,failure)
    IF (failure) RETURN
    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global
    CPPrecondition(nrow_global==ncol_global,cp_failure_level,routineP,error,failure)
    IF (failure) RETURN

    CALL timeset(routineN,handle)

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
    CALL timestop(handle)

  END SUBROUTINE cp_fm_transpose

! *****************************************************************************
!> \brief given an upper triangular matrix computes the corresponding full matrix
!> \param matrix the upper triangular matrix as input, the full matrix as output
!> \param work a matrix of the same size as matrix
!> \param error ...
!> \author Matthias Krack
!> \note
!>       the lower triangular part is irrelevant
! *****************************************************************************
  SUBROUTINE cp_fm_upper_to_full(matrix,work,error)

    TYPE(cp_fm_type), POINTER                :: matrix,work
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_upper_to_full',&
         routineP=moduleN//':'//routineN

    INTEGER                                  :: handle,icol_global,irow_global,&
                                                mypcol,myprow,ncol_global,&
                                                npcol,nprow,nrow_global
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: a
    REAL(KIND=sp), DIMENSION(:,:), POINTER   :: a_sp
    TYPE(cp_blacs_env_type), POINTER         :: context
    LOGICAL                                  :: failure

#if defined(__SCALAPACK)
    INTEGER                                  :: icol_local,irow_local,&
                                                ncol_block,ncol_local,&
                                                nrow_block,nrow_local
    INTEGER, DIMENSION(9)                    :: desca,descc
    INTEGER, EXTERNAL                        :: indxl2g
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: c
    REAL(KIND=sp), DIMENSION(:,:), POINTER   :: c_sp
#endif

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
    CPPrecondition(matrix%use_sp.EQV.work%use_sp,cp_failure_level,routineP,error,failure)
    IF (failure) RETURN

    CALL timeset(routineN,handle)

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
    a_sp => work%local_data_sp
    desca(:) = work%matrix_struct%descriptor(:)
    c => matrix%local_data
    c_sp => matrix%local_data_sp
    descc(:) = matrix%matrix_struct%descriptor(:)

    DO icol_local=1,ncol_local
      icol_global = indxl2g(icol_local,ncol_block,mypcol,&
           matrix%matrix_struct%first_p_pos(2),npcol)
      DO irow_local=1,nrow_local
        irow_global = indxl2g(irow_local,nrow_block,myprow,&
             matrix%matrix_struct%first_p_pos(1),nprow)
        IF (irow_global > icol_global) THEN
           IF(matrix%use_sp) THEN
              c_sp(irow_local,icol_local) = 0.0_sp
           ELSE
              c(irow_local,icol_local) = 0.0_dp
           ENDIF
        ELSE IF (irow_global == icol_global) THEN
           IF(matrix%use_sp) THEN
              c_sp(irow_local,icol_local) = 0.5_sp*c_sp(irow_local,icol_local)
           ELSE
              c(irow_local,icol_local) = 0.5_dp*c(irow_local,icol_local)
           ENDIF
        END IF
      END DO
    END DO

    DO icol_local=1,ncol_local
    DO irow_local=1,nrow_local
       IF(matrix%use_sp) THEN
          a_sp(irow_local,icol_local) = c_sp(irow_local,icol_local)
       ELSE
          a(irow_local,icol_local) = c(irow_local,icol_local)
       ENDIF
    ENDDO
    ENDDO

    IF(matrix%use_sp) THEN
       CALL pstran(nrow_global,ncol_global,1.0_sp,a_sp(1,1),1,1,desca,1.0_sp,c_sp(1,1),1,1,descc)
    ELSE
       CALL pdtran(nrow_global,ncol_global,1.0_dp,a(1,1),1,1,desca,1.0_dp,c(1,1),1,1,descc)
    ENDIF

#else

    a => matrix%local_data
    a_sp => matrix%local_data_sp
    DO irow_global=1,nrow_global
       DO icol_global=irow_global+1,ncol_global
          IF(matrix%use_sp) THEN
             a_sp(icol_global,irow_global)=a_sp(irow_global,icol_global)
          ELSE
             a(icol_global,irow_global)=a(irow_global,icol_global)
          ENDIF
       ENDDO
    ENDDO

#endif
    CALL timestop(handle)

  END SUBROUTINE cp_fm_upper_to_full

! *****************************************************************************
!> \brief scales column i of matrix a with scaling(i)
!> \param matrixa ...
!> \param scaling : an array used for scaling the columns, SIZE(scaling) determines the number of columns to be scaled
!> \author Joost VandeVondele
!> \note
!>      this is very useful as a first step in the computation of C = sum_i alpha_i A_i transpose (A_i)
!>      that is a rank-k update (cp_fm_syrk , cp_sm_plus_fm_fm_t) where every vector has a different prefactor
!>      this procedure can be up to 20 times faster than calling cp_fm_syrk n times
! *****************************************************************************
  SUBROUTINE cp_fm_column_scale(matrixa,scaling)
    TYPE(cp_fm_type), POINTER                :: matrixa
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: scaling

    INTEGER                                  :: k,mypcol,myprow,n,ncol_global,&
                                                npcol,nprow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
#if defined(__SCALAPACK)
    INTEGER                                  :: icol_global,icol_local,&
                                                ipcol,iprow,irow_local
#else
    INTEGER                                  :: i
#endif

    myprow=matrixa%matrix_struct%context%mepos(1)
    mypcol=matrixa%matrix_struct%context%mepos(2)
    nprow=matrixa%matrix_struct%context%num_pe(1)
    npcol=matrixa%matrix_struct%context%num_pe(2)

    ncol_global =matrixa%matrix_struct%ncol_global

    a => matrixa%local_data
    a_sp => matrixa%local_data_sp
    IF(matrixa%use_sp) THEN
       n =  SIZE(a_sp,1)
    ELSE
       n =  SIZE(a,1)
    ENDIF
    k =  MIN(SIZE(scaling),ncol_global)

#if defined(__SCALAPACK)

    DO icol_global=1,k
       CALL infog2l(1,icol_global,matrixa%matrix_struct%descriptor,&
            nprow,npcol,myprow,mypcol,&
            irow_local,icol_local,iprow,ipcol)
       IF ((ipcol == mypcol)) THEN
          IF(matrixa%use_sp) THEN
             CALL SSCAL(n,REAL(scaling(icol_global),sp),a_sp(1,icol_local),1)
          ELSE
             CALL DSCAL(n,scaling(icol_global),a(1,icol_local),1)
          ENDIF
       END IF
    ENDDO
#else
    DO i=1,k
       IF(matrixa%use_sp) THEN
          CALL SSCAL(n,REAL(scaling(i),sp),a_sp(1,i),1)
       ELSE
          CALL DSCAL(n,scaling(i),a(1,i),1)
       ENDIF
    ENDDO
#endif
  END SUBROUTINE cp_fm_column_scale

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_inverse ...
!> \param det_a ...
!> \param error ...
!> \author Florian Schiffmann(02.2007)
! *****************************************************************************
  SUBROUTINE cp_fm_invert(matrix_a,matrix_inverse,det_a,error)

    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_inverse
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: det_a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_invert', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: determinant
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    TYPE(cp_fm_type), POINTER                :: matrix_B, matrix_lu
#if defined(__SCALAPACK)
    INTEGER                                  :: i, info, liwork, lwork
    INTEGER, DIMENSION(9)                    :: desca
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    REAL(KIND=dp)                            :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: diag
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: berr, ferr, work
#else
    LOGICAL                                  :: sign
    REAL(KIND=dp)                            :: eps1
#endif

    failure=.FALSE.

    CALL cp_fm_create(matrix=matrix_lu,&
         matrix_struct=matrix_a%matrix_struct,&
         name="A_lu"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX",&
         error=error)
    CALL cp_fm_to_fm(matrix_a,matrix_lu,error=error)

    CALL cp_fm_create(matrix=matrix_B,&
         matrix_struct=matrix_a%matrix_struct,&
         name="B_mat"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX",&
         error=error)
    a => matrix_lu%local_data
    n = matrix_lu%matrix_struct%nrow_global
    ALLOCATE(ipivot(n+matrix_a%matrix_struct%nrow_block))
#if defined(__SCALAPACK)
    ALLOCATE(ferr(n))
    ALLOCATE(berr(n))
    ALLOCATE(work(3*n))
    ALLOCATE(iwork(3*N))
    lwork=3*n
    liwork=3*n
    desca(:) = matrix_lu%matrix_struct%descriptor(:)
    CALL pdgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    ALLOCATE(diag(n))
    diag(:)=0.0_dp
    DO i=1,n
       CALL cp_fm_get_element(matrix_lu,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    determinant=1.0_dp
    DO i=1,n
       IF(ipivot(i)==i)THEN
          determinant=determinant*diag(i)
       ELSE
          determinant=determinant*diag(i)*(-1.0_dp)
       END IF
    ENDDO
    DEALLOCATE(diag)

    alpha=0.0_dp
    beta=1.0_dp
   CALL cp_fm_set_all(matrix_inverse,alpha,beta,error)
   CALL pdgetrs('N',n,n,matrix_lu%local_data,1,1,desca,ipivot,matrix_inverse%local_data,1,1,desca,info)
!   CALL cp_fm_set_all(matrix_B,alpha,beta,error)
!    DO iter=1,10
!       CALL pdgerfs('N',n,n,matrix_a%local_data,1,1,desca,matrix_lu%local_data,&
!                      1,1,desca,ipivot,matrix_B%local_data,&
!                     1,1,desca,matrix_inverse%local_data,1,1,&
!                      desca,ferr,berr,work,lwork,iwork,liwork,info)
!       eps1=eps2
!       eps2=MAXVAL(ferr)
!       IF (ABS( eps2 - eps1) <= EPSILON(1.0_dp))THEN
!          EXIT
!       END IF
!    END DO
    DEALLOCATE(ferr)
    DEALLOCATE(berr)
    DEALLOCATE(work)
    DEALLOCATE(iwork)

#else
    sign=.TRUE.
    CALL invert_matrix(matrix_a%local_data,matrix_inverse%local_data,&
                       eval_error=eps1,error=error)
    CALL cp_fm_lu_decompose(matrix_lu,determinant,correct_sign=sign)
#endif
    CALL cp_fm_release(matrix_lu,error=error)
    CALL cp_fm_release(matrix_B,error=error)
    DEALLOCATE(ipivot)
    IF(PRESENT(det_a)) det_a = determinant
  END SUBROUTINE cp_fm_invert

! *****************************************************************************
!> \brief inverts a triangular matrix
!> \param matrix_a ...
!> \param uplo_tr ...
!> \param error ...
!> \author MI
! *****************************************************************************
  SUBROUTINE cp_fm_triangular_invert(matrix_a,uplo_tr,error)

    TYPE(cp_fm_type), POINTER                :: matrix_a
    CHARACTER, INTENT(IN), OPTIONAL          :: uplo_tr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_triangular_invert', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: unit_diag, uplo
    INTEGER                                  :: handle, info, ncol_global
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca
#endif

    failure = .FALSE.

    CALL timeset(routineN,handle)

    unit_diag='N'
    uplo='U'
    IF(PRESENT(uplo_tr)) uplo=uplo_tr

    ncol_global =matrix_a%matrix_struct%ncol_global

    a => matrix_a%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)

    CALL pdtrtri( uplo, unit_diag, ncol_global, a(1,1), 1, 1, desca, info )

#else
    CALL dtrtri( uplo, unit_diag, ncol_global, a(1,1), ncol_global, info )
#endif


    CALL timestop(handle)
  END SUBROUTINE cp_fm_triangular_invert


! *****************************************************************************
!> \brief  perfoms a QR factorization of the input rectangular matrix A or of a submatrix of A
!>         the computed upper triangular matrix R is in output in the submatrix sub(A) of size NxN
!>         M and M give the dimention of the submatrix that has to be factorized (MxN) with M>N
!> \param matrix_a ...
!> \param matrix_r ...
!> \param nrow_fact ...
!> \param ncol_fact ...
!> \param first_row ...
!> \param first_col ...
!> \param error ...
!> \author MI
! *****************************************************************************
  SUBROUTINE cp_fm_qr_factorization(matrix_a, matrix_r, nrow_fact, ncol_fact, first_row, first_col, error)
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_r
    INTEGER, INTENT(IN), OPTIONAL            :: nrow_fact, ncol_fact, &
                                                first_row, first_col
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_qr_factorization', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icol, info, irow, &
                                                istat, j, lda, lwork, ncol, &
                                                ndim, nrow
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: tau, work
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: r_mat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca
#endif

    failure = .FALSE.

    CALL timeset(routineN,handle)

    ncol =matrix_a%matrix_struct%ncol_global
    nrow =matrix_a%matrix_struct%nrow_global
    lda = nrow

    a => matrix_a%local_data

    IF(PRESENT(nrow_fact)) nrow = nrow_fact
    IF(PRESENT(ncol_fact)) ncol = ncol_fact
    irow = 1
    IF(PRESENT(first_row)) irow = first_row
    icol = 1
    IF(PRESENT(first_col)) icol = first_col


    CPPrecondition(nrow>=ncol,cp_failure_level,routineP,error,failure)
    ndim = SIZE(a,2)
!    ALLOCATE(ipiv(ndim),istat=STAT)
    ALLOCATE(tau(ndim),STAT=istat)

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)

    lwork = -1
    ALLOCATE(work(2*ndim),STAT=istat)
    CALL pdgeqrf( nrow, ncol, a, irow, icol, desca,  tau, work, lwork, info )
    lwork = INT(work(1))
    DEALLOCATE(work,STAT=istat)
    ALLOCATE(work(lwork),STAT=istat)
    CALL pdgeqrf( nrow, ncol, a, irow, icol, desca, tau, work, lwork, info )

#else
    lwork = -1
    ALLOCATE(work(2*ndim),STAT=istat)
    CALL dgeqrf( nrow, ncol, a, lda,  tau, work, lwork, info )
    lwork = INT(work(1))
    DEALLOCATE(work,STAT=istat)
    ALLOCATE(work(lwork),STAT=istat)
    CALL dgeqrf( nrow, ncol, a, lda,  tau, work, lwork, info )

#endif

    ALLOCATE(r_mat(ncol,ncol),STAT=istat)
    CALL cp_fm_get_submatrix(matrix_a,r_mat,1,1,ncol,ncol,error=error)
    DO i = 1,ncol
      DO j = i+1,ncol
          r_mat(j,i) = 0.0_dp
      END DO
    END DO
    CALL cp_fm_set_submatrix(matrix_r,r_mat,1,1,ncol,ncol,error=error)


    DEALLOCATE(tau, work, r_mat, STAT=istat)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_qr_factorization

! *****************************************************************************
!> \brief computs the the solution to A*b=A_general using lu decomposition
!>        pay attention, both matrices are overwritten, a_general contais the result
!> \param matrix_a ...
!> \param general_a ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE cp_fm_solve(matrix_a,general_a,error)
    TYPE(cp_fm_type), POINTER                :: matrix_a, general_a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_solve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, a_general
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb
#else
    INTEGER                                  :: lda, ldb
#endif

    CALL timeset(routineN,handle)

    a => matrix_a%local_data
    a_general => general_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE(ipivot(n+matrix_a%matrix_struct%nrow_block))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = general_a%matrix_struct%descriptor(:)
    CALL pdgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    CALL pdgetrs("N" , n , n , a(1,1), 1, 1, desca ,ipivot, a_general(1,1) ,&
                 1, 1, descb, info )

#else
    lda=SIZE(a,1)
    ldb=SIZE(a_general,1)
    CALL dgetrf(n,n,a(1,1),lda,ipivot,info)
    CALL dgetrs("N",n,n,a(1,1),lda,ipivot,a_general,ldb,info)

#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
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
!> \param error ...
!> \param b_first_col : the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!> \param a_first_row ...
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
!> \author Samuel Andermatt
!> \note
!>      C should have no overlap with A, B
! *****************************************************************************
SUBROUTINE cp_complex_fm_gemm(transa,transb,m,n,k,alpha,A_re,A_im,B_re,B_im,beta,&
       C_re,C_im,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: A_re, A_im, B_re, B_im
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: C_re, C_im
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: b_first_col, a_first_row, &
                                                b_first_row, c_first_col, &
                                                c_first_row

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_complex_fm_gemm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL cp_fm_gemm(transa,transb,m,n,k,alpha,A_re,B_re,beta,&
       C_re,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)
    CALL cp_fm_gemm(transa,transb,m,n,k,-alpha,A_im,B_im,1.0_dp,&
       C_re,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)
    CALL cp_fm_gemm(transa,transb,m,n,k,alpha,A_re,B_im,beta,&
       C_im,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)
    CALL cp_fm_gemm(transa,transb,m,n,k,alpha,A_im,B_re,1.0_dp,&
       C_im,error,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)

    CALL timestop(handle)

END SUBROUTINE

END MODULE cp_fm_basic_linalg
