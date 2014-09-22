!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basic linear algebra operations for complex full matrixes
!> \note
!>      - not all functionality implemented
!> \par History
!>      Nearly literal copy of Fawzi's routines
!> \author Joost VandeVondele
! *****************************************************************************
MODULE cp_cfm_basic_linalg
  USE cp_cfm_types,                    ONLY: cp_cfm_get_element,&
                                             cp_cfm_get_info,&
                                             cp_cfm_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_basic_linalg'

  PUBLIC :: cp_cfm_add, cp_cfm_add_fm, cp_cfm_gemm, cp_cfm_lu_decompose,&
            cp_cfm_column_scale, cp_cfm_schur_product,&
            cp_cfm_solve, cp_cfm_cholesky_decompose, cp_cfm_triangular_multiply,&
            cp_cfm_triangular_invert

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_cfm_schur_product(matrix_a,matrix_b,matrix_c,error)

    TYPE(cp_cfm_type), POINTER               :: matrix_a, matrix_b, matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_schur_product', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b, c
    INTEGER                                  :: handle, icol_local, &
                                                irow_local, mypcol, myprow, &
                                                ncol_local, nrow_local

    CALL timeset(routineN,handle)

    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)

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

  END SUBROUTINE cp_cfm_schur_product

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_cfm_schur_product_cc(matrix_a,matrix_b,matrix_c,error)

    TYPE(cp_cfm_type), POINTER               :: matrix_a, matrix_b, matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_schur_product_cc', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b, c
    INTEGER                                  :: handle, icol_local, &
                                                irow_local, mypcol, myprow, &
                                                ncol_local, nrow_local

    CALL timeset(routineN,handle)

    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

    nrow_local = matrix_a%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix_a%matrix_struct%ncol_locals(mypcol)

    DO icol_local=1,ncol_local
       DO irow_local=1,nrow_local
          c(irow_local,icol_local) = a(irow_local,icol_local)*CONJG(b(irow_local,icol_local))
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE cp_cfm_schur_product_cc
! *****************************************************************************
!> \brief   Scale and add two BLACS matrices (a <- alpha*a + beta*b).
!> \param alpha ...
!> \param matrix_a ...
!> \param beta ...
!> \param matrix_b ...
!> \param error ...
!> \date    11.06.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_add(alpha,matrix_a,beta,matrix_b,error)
    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    TYPE(cp_cfm_type), POINTER               :: matrix_a
    COMPLEX(KIND=dp), INTENT(in), OPTIONAL   :: beta
    TYPE(cp_cfm_type), OPTIONAL, POINTER     :: matrix_b
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_add', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: my_beta
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b
    INTEGER                                  :: handle, mypcol, myprow
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_beta=CMPLX(0.0_dp,0.0_dp,dp)
    IF(PRESENT(beta)) my_beta=beta
    NULLIFY(a,b)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)
    ! to do: use dscal,dcopy,daxp
    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)

    a => matrix_a%local_data

    IF (my_beta==0.0_dp) THEN

       IF (alpha==CMPLX(0.0_dp,0.0_dp,dp)) THEN
          a(:,:)=CMPLX(0.0_dp,0.0_dp,dp)
       ELSE IF (alpha == CMPLX(1.0_dp,0.0_dp,dp)) THEN
          RETURN
       ELSE
          a(:,:) = alpha*a(:,:)
       END IF

    ELSE
       CPPrecondition(PRESENT(matrix_b),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(matrix_b),cp_failure_level,routineP,error,failure)
       CPPrecondition(matrix_b%ref_count>0,cp_failure_level,routineP,error,failure)
       CALL cp_assert(matrix_a%matrix_struct%context%group==&
            matrix_b%matrix_struct%context%group,cp_failure_level,&
            cp_assertion_failed,fromWhere=routineP,&
            message="matrixes must be in the same blacs context"//&
CPSourceFileRef,&
            error=error,failure=failure)

       IF (cp_fm_struct_equivalent(matrix_a%matrix_struct,&
            matrix_b%matrix_struct,error=error)) THEN

          b => matrix_b%local_data

          IF (alpha == CMPLX(0.0_dp,0.0_dp,dp)) THEN
             IF (my_beta == CMPLX(1.0_dp,0.0_dp,dp)) THEN
                a(:,:) = b(:,:)
             ELSE
                a(:,:) = my_beta*b(:,:)
             END IF
          ELSE IF (alpha == CMPLX(1.0_dp,0.0_dp,dp)) THEN
             IF (my_beta == CMPLX(1.0_dp,0.0_dp,dp)) THEN
                a(:,:) = a(:,:) + b(:,:)
             ELSE
                a(:,:) = a(:,:) + my_beta*b(:,:)
             END IF
          ELSE
             a(:,:) = alpha*a(:,:) + my_beta*b(:,:)
          END IF
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

  END SUBROUTINE cp_cfm_add
! *****************************************************************************
!> \brief   Scale and add two BLACS matrices (a <- alpha*a + beta*b).
!>          where b is a real matrix (adapted from cp_cfm_add)
!> \param alpha ...
!> \param matrix_a ...
!> \param beta ...
!> \param matrix_b ...
!> \param error ...
!> \date    01.08.2014
!> \author  JGH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_add_fm(alpha,matrix_a,beta,matrix_b,error)
    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    TYPE(cp_cfm_type), POINTER               :: matrix_a
    COMPLEX(KIND=dp), INTENT(IN)             :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_b
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_add_fm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: handle, mypcol, myprow
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: b

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(a,b)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)

    a => matrix_a%local_data

    IF (beta==CMPLX(0.0_dp,0.0_dp,dp)) THEN

       IF (alpha==CMPLX(0.0_dp,0.0_dp,dp)) THEN
          a(:,:)=CMPLX(0.0_dp,0.0_dp,dp)
       ELSE IF (alpha == CMPLX(1.0_dp,0.0_dp,dp)) THEN
          RETURN
       ELSE
          a(:,:) = alpha*a(:,:)
       END IF

    ELSE
       CPPrecondition(ASSOCIATED(matrix_b),cp_failure_level,routineP,error,failure)
       CALL cp_assert(matrix_a%matrix_struct%context%group==&
            matrix_b%matrix_struct%context%group,cp_failure_level,&
            cp_assertion_failed,fromWhere=routineP,&
            message="matrices must be in the same blacs context"//&
CPSourceFileRef,&
            error=error,failure=failure)

       IF (cp_fm_struct_equivalent(matrix_a%matrix_struct,&
            matrix_b%matrix_struct,error=error)) THEN

          b => matrix_b%local_data

          IF (alpha == CMPLX(0.0_dp,0.0_dp,dp)) THEN
             IF (beta == CMPLX(1.0_dp,0.0_dp,dp)) THEN
                a(:,:) = b(:,:)
             ELSE
                a(:,:) = beta*b(:,:)
             END IF
          ELSE IF (alpha == CMPLX(1.0_dp,0.0_dp,dp)) THEN
             IF (beta == CMPLX(1.0_dp,0.0_dp,dp)) THEN
                a(:,:) = a(:,:) + b(:,:)
             ELSE
                a(:,:) = a(:,:) + beta*b(:,:)
             END IF
          ELSE
             a(:,:) = alpha*a(:,:) + beta*b(:,:)
          END IF
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

  END SUBROUTINE cp_cfm_add_fm

! *****************************************************************************
!> \brief   Computes the LU decomposition of a given matrix
!>          the actual purpose right now is to compute the determinant of a given matrix
!>          which is most efficiently done this way, but, indeed, destroys the matrix
!>          SERIOUS WARNING (KNOWN BUG) : the sign of the determinant depends on ipivot
!>          one should be able to find out if ipivot is an even or an odd permutation...
!> \param matrix_a ...
!> \param almost_determinant ...
!> \param error ...
!> \date    11.06.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_lu_decompose(matrix_a,almost_determinant,error)
    TYPE(cp_cfm_type), POINTER               :: matrix_a
    COMPLEX(KIND=dp), INTENT(OUT)            :: almost_determinant
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_lu_decompose', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: determinant
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: counter, handle, i, info, &
                                                n, nrow_local
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    INTEGER, DIMENSION(:), POINTER           :: row_indices
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: diag
#else
    INTEGER                                  :: lda
#endif

    CALL timeset(routineN,handle)

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE(ipivot(n))
    CALL cp_cfm_get_info(matrix_a,row_indices=row_indices,nrow_local=nrow_local,error=error)
#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    CALL pzgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    ALLOCATE(diag(n))
    diag(:)=CMPLX(0.0_dp,0.0_dp,dp)
    counter=0
    DO i=1,nrow_local
       IF(ipivot(i).ne.row_indices(i))counter=counter+1
    END DO
    DO i=1,n
       CALL cp_cfm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    determinant=CMPLX(1.0_dp,0.0_dp,dp)
    DO i=1,n
       determinant=determinant*diag(i)
    ENDDO
    CALL mp_sum(counter,matrix_a%matrix_struct%para_env%group)
    IF(MOD(counter,2)==1)determinant=-1.0_dp*determinant
    DEALLOCATE(diag)
#else
    lda=SIZE(a,1)
    CALL zgetrf(n,n,a(1,1),lda,ipivot,info)
    counter=0
    determinant=CMPLX(1.0_dp,0.0_dp,dp)
    DO i=1,n
       IF(ipivot(i).NE.i)counter=counter+1
       determinant=determinant*a(i,i)
    ENDDO
    IF(MOD(counter,2)==1)determinant=-1.0_dp*determinant
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    almost_determinant=determinant ! notice that the sign is random
    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief   BLACS interface to the BLAS routine dgemm.
!> \param transa ...
!> \param transb ...
!> \param m ...
!> \param n ...
!> \param k ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param error ...
!> \param b_first_col ...
!> \date    07.06.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error,b_first_col)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    TYPE(cp_cfm_type), POINTER               :: matrix_a, matrix_b
    COMPLEX(KIND=dp), INTENT(IN)             :: beta
    TYPE(cp_cfm_type), POINTER               :: matrix_c
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: b_first_col

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_gemm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b, c
    INTEGER                                  :: handle, j
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb, descc
#else
    INTEGER                                  :: lda, ldb, ldc
#endif

    CALL timeset(routineN,handle)
    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

    IF (PRESENT(b_first_col)) THEN
      j = b_first_col
    ELSE
      j = 1
    END IF

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = matrix_b%matrix_struct%descriptor(:)
    descc(:) = matrix_c%matrix_struct%descriptor(:)

    CALL pzgemm(transa,transb,m,n,k,alpha,a(1,1),1,1,desca,b(1,1),1,j,&
                descb,beta,c(1,1),1,j,descc)
#else

    lda = SIZE(a,1)
    ldb = SIZE(b,1)
    ldc = SIZE(c,1)

    CALL zgemm(transa,transb,m,n,k,alpha,a(1,1),lda,b(1,j),ldb,beta,c(1,j),ldc)

#endif
    CALL timestop(handle)

  END SUBROUTINE cp_cfm_gemm

! *****************************************************************************
!> \brief scales column i of matrix a with scaling(i)
!> \param matrixa ...
!> \param scaling : an array used for scaling the columns, SIZE(scaling) determines the number of columns to be scaled
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_cfm_column_scale(matrixa,scaling)
    TYPE(cp_cfm_type), POINTER               :: matrixa
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(in)                             :: scaling

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: k, mypcol, myprow, n, npcol, &
                                                nprow
#if defined(__SCALAPACK)
    INTEGER                                  :: icol_global, icol_local, &
                                                ipcol, iprow, irow_local
#else
    INTEGER                                  :: i
#endif

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
           CALL ZSCAL(n,scaling(icol_global),a(1,icol_local),1)
       END IF
    END DO
#else
    DO i=1,k
      CALL ZSCAL(n,scaling(i),a(1,i),1)
    END DO
#endif

  END SUBROUTINE cp_cfm_column_scale

! *****************************************************************************
!> \brief computs the the solution to A*b=A_general using lu decomposition
!>        pay attention, both matrices are overwritten, a_general contais the result
!> \param matrix_a ...
!> \param general_a ...
!> \param determinant ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE cp_cfm_solve(matrix_a,general_a,determinant,error)
    TYPE(cp_cfm_type), POINTER               :: matrix_a, general_a
    COMPLEX(KIND=dp), OPTIONAL               :: determinant
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_solve', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, a_general
    INTEGER                                  :: counter, handle, i, info, &
                                                n, nrow_local
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    INTEGER, DIMENSION(:), POINTER           :: row_indices
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: diag
#else
    INTEGER                                  :: lda, ldb
#endif

    CALL timeset(routineN,handle)

    a => matrix_a%local_data
    a_general => general_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    CALL cp_cfm_get_info(matrix_a,row_indices=row_indices,nrow_local=nrow_local,error=error)
    ALLOCATE(ipivot(n))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = general_a%matrix_struct%descriptor(:)
    CALL pzgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    IF(PRESENT(determinant))THEN
       ALLOCATE(diag(n))
       diag(:)=CMPLX(0.0_dp,0.0_dp,dp)
       counter=0
       DO i=1,nrow_local
          IF(ipivot(i).ne.row_indices(i))counter=counter+1
       END DO
       DO i=1,n
          CALL cp_cfm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
       ENDDO
       determinant=CMPLX(1.0_dp,0.0_dp,dp)
       DO i=1,n
          determinant=determinant*diag(i)
       ENDDO
       CALL mp_sum(counter,matrix_a%matrix_struct%para_env%group)
       IF(MOD(counter,2)==1)determinant=-1.0_dp*determinant
       DEALLOCATE(diag)
    END IF
    CALL pzgetrs("N" , n , n , a(1,1), 1, 1, desca ,ipivot, a_general(1,1) ,&
                 1, 1, descb, info )
#else
    lda=SIZE(a,1)
    ldb=SIZE(a_general,1)
    CALL zgetrf(n,n,a(1,1),lda,ipivot,info)
    IF(PRESENT(determinant))THEN
       counter=0
       determinant=CMPLX(1.0_dp,0.0_dp,dp)
       DO i=1,n
          IF(ipivot(i).NE.i)counter=counter+1
          determinant=determinant*a(i,i)
       ENDDO
       IF(MOD(counter,2)==1)determinant=-1.0_dp*determinant
    END IF
    CALL zgetrs("N",n,n,a(1,1),lda,ipivot,a_general,ldb,info)
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    CALL timestop(handle)

  END SUBROUTINE cp_cfm_solve

! *****************************************************************************
!> \brief used to replace a symmetric positive def. matrix M with its cholesky
!>      decomposition U: M = U^T * U, with U upper triangular
!> \param matrix the matrix to replace with its cholesky decomposition
!> \param n the number of row (and columns) of the matrix &
!>        (defaults to the min(size(matrix)))
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!> \author Joost
! *****************************************************************************
  SUBROUTINE cp_cfm_cholesky_decompose(matrix,n,error)
    TYPE(cp_cfm_type), POINTER               :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_cholesky_decompose', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: handle, info, my_n
    LOGICAL                                  :: failure
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_n = MIN(matrix%matrix_struct%nrow_global,&
         matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       CPPrecondition(n<=my_n,cp_failure_level,routineP,error,failure)
       my_n=n
    END IF

    a => matrix%local_data

#if defined(__SCALAPACK)
    desca(:) = matrix%matrix_struct%descriptor(:)
    CALL pzpotrf('U',my_n,a(1,1),1,1,desca,info)

#else

    CALL zpotrf('U',my_n,a(1,1),SIZE(a,1),info)

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END  SUBROUTINE cp_cfm_cholesky_decompose

! *****************************************************************************
!> \brief multiplies in place by a triangular matrix:
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs an mpi env
! *****************************************************************************
SUBROUTINE cp_cfm_triangular_multiply(triangular_matrix,matrix_b,side,&
     transa_tr, invert_tr, uplo_tr,unit_diag_tr, n_rows, n_cols, &
     alpha,error)
    TYPE(cp_cfm_type), POINTER               :: triangular_matrix, matrix_b
    CHARACTER, INTENT(in), OPTIONAL          :: side, transa_tr
    LOGICAL, INTENT(in), OPTIONAL            :: invert_tr
    CHARACTER, INTENT(in), OPTIONAL          :: uplo_tr
    LOGICAL, INTENT(in), OPTIONAL            :: unit_diag_tr
    INTEGER, INTENT(in), OPTIONAL            :: n_rows, n_cols
    COMPLEX(KIND=dp), INTENT(in), OPTIONAL   :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_triangular_multiply', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: side_char, transa, unit_diag, &
                                                uplo
    INTEGER                                  :: handle, m, n
    LOGICAL                                  :: failure, invert
    COMPLEX(KIND=dp)                            :: al

  failure=.FALSE.

  CALL timeset(routineN,handle)
  side_char='L'
  unit_diag='N'
  uplo='U'
  transa='N'
  invert=.FALSE.
  al=CMPLX(1.0_dp,0.0_dp,dp)
  CALL cp_cfm_get_info(matrix_b, nrow_global=m, ncol_global=n, error=error)
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
  IF (PRESENT(transa_tr)) transa = transa_tr
  IF (PRESENT(alpha)) al=alpha
  IF (PRESENT(n_rows)) m=n_rows
  IF (PRESENT(n_cols)) n=n_cols

  IF (invert) THEN

#if defined(__SCALAPACK)
     CALL pztrsm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),1,1,&
          triangular_matrix%matrix_struct%descriptor,&
          matrix_b%local_data(1,1),1,1,&
          matrix_b%matrix_struct%descriptor(1))
#else
     CALL ztrsm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),&
          SIZE(triangular_matrix%local_data,1),&
          matrix_b%local_data(1,1),SIZE(matrix_b%local_data,1))
#endif

  ELSE

#if defined(__SCALAPACK)
     CALL pztrmm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),1,1,&
          triangular_matrix%matrix_struct%descriptor,&
          matrix_b%local_data(1,1),1,1,&
          matrix_b%matrix_struct%descriptor(1))
#else
     CALL ztrmm(side_char,uplo,transa,unit_diag,m,n,al,&
          triangular_matrix%local_data(1,1),&
          SIZE(triangular_matrix%local_data,1),&
          matrix_b%local_data(1,1),SIZE(matrix_b%local_data,1))
#endif

  END IF

  CALL timestop(handle)

  END SUBROUTINE cp_cfm_triangular_multiply

! *****************************************************************************
!> \brief inverts a triangular matrix
!> \param matrix_a ...
!> \param uplo_tr ...
!> \param error ...
!> \author MI
! *****************************************************************************
  SUBROUTINE cp_cfm_triangular_invert(matrix_a,uplo_tr,error)

    TYPE(cp_cfm_type), POINTER               :: matrix_a
    CHARACTER, INTENT(IN), OPTIONAL          :: uplo_tr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_triangular_invert', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: unit_diag, uplo
    INTEGER                                  :: handle, info, ncol_global
    LOGICAL                                  :: failure
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
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

    CALL pztrtri( uplo, unit_diag, ncol_global, a(1,1), 1, 1, desca, info )

#else
    CALL ztrtri( uplo, unit_diag, ncol_global, a(1,1), ncol_global, info )
#endif


    CALL timestop(handle)
  END SUBROUTINE cp_cfm_triangular_invert

END MODULE cp_cfm_basic_linalg
