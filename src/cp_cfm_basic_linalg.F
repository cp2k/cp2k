!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_cfm_basic_linalg [1.0] *
!!
!!   NAME
!!     cp_cfm_basic_linalg
!!
!!   FUNCTION
!!     basic linear algebra operations for complex full matrixes
!!
!!   NOTES
!!     - not all functionality implemented
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     Nearly literal copy of Fawzi's routines
!!
!!   SOURCE
!****************************************************************************
MODULE cp_cfm_basic_linalg
  USE cp_cfm_types,                    ONLY: cp_cfm_get_element,&
                                             cp_cfm_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_cfm_basic_linalg'

  PUBLIC :: cp_cfm_add, cp_cfm_gemm, cp_cfm_lu_decompose, cp_cfm_column_scale, cp_cfm_schur_product
!!***
!****************************************************************************
CONTAINS

  SUBROUTINE cp_cfm_schur_product(matrix_a,matrix_b,matrix_c,error)

    TYPE(cp_cfm_type), POINTER               :: matrix_a, matrix_b, matrix_c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b, c
    INTEGER                                  :: handle, icol_local, &
                                                irow_local, mypcol, myprow, &
                                                ncol_local, nrow_local

!   ---------------------------------------------------------------------------

    CALL timeset("cp_cfm_schur_product","I","",handle)

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

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_cfm_schur_product


  SUBROUTINE cp_cfm_add(alpha,matrix_a,beta,matrix_b,error)

    !   Purpose: Scale and add two BLACS matrices (a <- alpha*a + beta*b).

    !   History: - Creation (11.06.2001, Matthias Krack)

    !   ***************************************************************************

    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    TYPE(cp_cfm_type), POINTER               :: matrix_a
    COMPLEX(KIND=dp), INTENT(in), OPTIONAL   :: beta
    TYPE(cp_cfm_type), OPTIONAL, POINTER     :: matrix_b
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_add', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: my_beta
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b
    INTEGER                                  :: handle, mypcol, myprow
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

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
          a(:,:)=CMPLX(0.0_dp,0.0_dp)
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
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_cfm_add

! *****************************************************************************

! *****************************************************************************
! computes the LU decomposition of a given matrix
! the actual purpose right now is to compute the determinant of a given matrix
! which is most efficiently done this way, but, indeed, destroys the matrix
! SERIOUS WARNING (KNOWN BUG) : the sign of the determinant depends on ipivot
! one should be able to find out if ipivot is an even or an odd permutation...
! *****************************************************************************
  SUBROUTINE cp_cfm_lu_decompose(matrix_a,almost_determinant)
    TYPE(cp_cfm_type), POINTER               :: matrix_a
    COMPLEX(KIND=dp), INTENT(OUT)            :: almost_determinant

    CHARACTER(LEN=*), PARAMETER :: routine = "SUBROUTINE cp_cfm_lu_decompose"

    COMPLEX(KIND=dp)                         :: determinant
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: diag
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: handle, i, info, lda, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    INTEGER, DIMENSION(9)                    :: desca

! *** locals ***

    CALL timeset("cp_cfm_lu_decompose","I","",handle)

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE(ipivot(n))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    CALL pzgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    ALLOCATE(diag(n))
    diag(:)=CMPLX(0.0_dp,0.0_dp,dp)
    DO i=1,n
       CALL cp_cfm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    determinant=CMPLX(1.0_dp,0.0_dp,dp)
    do i=1,n
       determinant=determinant*diag(i)
    enddo
    DEALLOCATE(diag)
#else
    lda=size(a,1)
    CALL zgetrf(n,n,a(1,1),lda,ipivot,info)
    determinant=CMPLX(1.0_dp,0.0_dp,dp)
    do i=1,n
       determinant=determinant*a(i,i)
    enddo
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    almost_determinant=determinant ! notice that the sign is random
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE


! *****************************************************************************

  SUBROUTINE cp_cfm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error,b_first_col)

!   Purpose: BLACS interface to the BLAS routine dgemm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    TYPE(cp_cfm_type), POINTER               :: matrix_a, matrix_b
    COMPLEX(KIND=dp), INTENT(IN)             :: beta
    TYPE(cp_cfm_type), POINTER               :: matrix_c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, INTENT(IN), OPTIONAL            :: b_first_col

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b, c
    INTEGER                                  :: handle, j, lda, ldb, ldc
    INTEGER, DIMENSION(9)                    :: desca, descb, descc
    REAL(KIND=dp)                            :: flops

!   ---------------------------------------------------------------------------

    CALL timeset("cp_cfm_gemm","I","Gflops",handle)

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

    flops=4*2.0_dp*(1E-9_dp*m*n*k)/(matrix_a%matrix_struct%context%num_pe(1)*&
         matrix_a%matrix_struct%context%num_pe(2))
    CALL pzgemm(transa,transb,m,n,k,alpha,a(1,1),1,1,desca,b(1,1),1,j,&
                descb,beta,c(1,1),1,j,descc)
#else

    lda = SIZE(a,1)
    ldb = SIZE(b,1)
    ldc = SIZE(c,1)

    flops=4*2.0_dp*(1E-9_dp*m*n*k)
    CALL zgemm(transa,transb,m,n,k,alpha,a(1,1),lda,b(1,j),ldb,beta,c(1,j),ldc)

#endif
    CALL timestop(flops,handle)

  END SUBROUTINE cp_cfm_gemm

!!****f* cp_cfm_basic_linalg/cp_cfm_column_scale [1.0] *
!!
!!   NAME
!!     cp_fm_column_scale
!!
!!   SYNOPSIS
!!     Subroutine cp_cfm_column_scale(matrixa, scaling)
!!       Type(cp_cfm_type), Pointer:: matrixa
!!       Complex(KIND = dp), Dimension(:), Intent (IN):: scaling
!!     End Subroutine cp_fm_column_scale
!!
!!   FUNCTION
!!     scales column i of matrix a with scaling(i)
!!
!!   NOTES
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
  SUBROUTINE cp_cfm_column_scale(matrixa,scaling)
    TYPE(cp_cfm_type), POINTER               :: matrixa
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(in)                             :: scaling

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: i, icol_global, icol_local, &
                                                ipcol, iprow, irow_local, k, &
                                                mypcol, myprow, n, npcol, &
                                                nprow

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
           CALL ZSCAL(n,scaling(icol_global),a(1,icol_local),1)
       END IF
    ENDDO
#else
    DO i=1,k
      CALL ZSCAL(n,scaling(i),a(1,i),1)
    ENDDO
#endif
  END SUBROUTINE cp_cfm_column_scale
END MODULE cp_cfm_basic_linalg
