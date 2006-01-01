!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_fm_cholesky [1.0] *
!!
!!   NAME
!!     cp_fm_cholesky
!!
!!   FUNCTION
!!     various cholesky decomposition related routines
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_cholesky
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_cholesky'

  PUBLIC :: cp_fm_cholesky_decompose, cp_fm_cholesky_invert,&
       cp_fm_cholesky_reduce, cp_fm_cholesky_restore

!***
!****************************************************************************
CONTAINS

!!****f* cp_log_handling/cp_fm_cholesky_decompose [1.0] *
!!
!!   NAME
!!     cp_fm_cholesky_decompose
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_cholesky_decompose(matrix, n, error)
!!       Type(cp_fm_type), Pointer:: matrix
!!       Integer, Intent (IN), Optional:: n
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_cholesky_decompose
!!
!!   FUNCTION
!!     used to replace a symmetric positive def. matrix M with its cholesky
!!     decomposition U: M = U^T * U, with U upper triangular
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to replace with its cholesky decomposition
!!     - n: the number of row (and columns) of the matrix &
!!       (defaults to the min(size(matrix)))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [JVdV]
!!     12.2002 updated, added n optional parm [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_cholesky_decompose(matrix,n,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "cp_fm_cholesky_decompose", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, my_n
    INTEGER, DIMENSION(9)                    :: desca
    LOGICAL                                  :: failure, ionode
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a

    CALL timeset(routineN,"I","",handle)

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
    CALL pdpotrf('U',my_n,a(1,1),1,1,desca,info)

#else

    CALL dpotrf('U',my_n,a(1,1),SIZE(a,1),info)

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END  SUBROUTINE cp_fm_cholesky_decompose
!***************************************************************************
  
!!****f* cp_fm_cholesky/cp_fm_cholesky_invert [1.0] *
!!
!!   NAME
!!     cp_fm_cholesky_invert
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_cholesky_invert(matrix, n, error)
!!       Type(cp_fm_type), Pointer:: matrix
!!       Integer, Intent (IN), Optional:: n
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_cholesky_invert
!!
!!   FUNCTION
!!     used to replace the cholesky decomposition by the inverse
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to invert (must be an upper triangular matrix)
!!     - n: size of the matrix to invert (defaults to the min(size(matrix)))
!!     - error : variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [JVdV]
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_cholesky_invert(matrix,n, error)
    TYPE(cp_fm_type), POINTER           :: matrix
    INTEGER, INTENT(in), OPTIONAL                :: n
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_invert',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: a
    INTEGER                                   :: info,handle
    INTEGER                                   :: my_n
    INTEGER, DIMENSION(9)                     :: desca
    LOGICAL                                   :: failure

    CALL timeset("cp_fm_cholesky_invert","I","",handle)
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
    CALL pdpotri('U',my_n,a(1,1),1,1,desca,info)

#else

    CALL dpotri('U',my_n,a(1,1),SIZE(a,1),info)

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END  SUBROUTINE cp_fm_cholesky_invert
!***************************************************************************

!!****f* cp_fm_cholesky/cp_fm_cholesky_reduce [1.0] *
!!
!!   NAME
!!     cp_fm_cholesky_reduce
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_cholesky_reduce(matrix, matrixb, error)
!!       Type(cp_fm_type), Pointer:: matrix, matrixb
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_cholesky_reduce
!!
!!   FUNCTION
!!     reduce a matrix pencil A,B to normal form
!!     B has to be cholesky decomposed with  cp_fm_cholesky_decompose
!!     before calling this routine
!!     A,B -> inv(U^T)*A*inv(U),1
!!     (AX=BX -> inv(U^T)*A*inv(U)*U*X=U*X hence evecs U*X)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the symmetric matrix A
!!     matrixb: the cholesky decomposition of matrix B
!!     error : variable to control error logging, stopping,... 
!!             see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [JVdV]
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_cholesky_reduce(matrix,matrixb, error)
   TYPE(cp_fm_type), POINTER           :: matrix,matrixb
   TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_reduce',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: a,b
    INTEGER                                   :: info,itype,handle
    INTEGER                                   :: n
    INTEGER                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca,descb
    REAL(KIND = dp)                                  :: scale
    LOGICAL                                   :: failure

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n = matrix%matrix_struct%nrow_global
    itype =1

    a => matrix%local_data
    b => matrixb%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix%matrix_struct%descriptor(:)
    descb(:) = matrixb%matrix_struct%descriptor(:)

    CALL pdsygst(itype,'U',n,a(1,1),1,1,desca,b(1,1),1,1,descb,scale,info)

    ! this is supposed to be one in current version of lapack
    ! if not, eigenvalues have to be scaled by this number
    CALL cp_assert(scale == 1.0_dp, cp_failure_level,cp_assertion_failed,&
         routineP,&
         "scale not equal 1 (scale="//cp_to_string(scale)//")"//&
         CPSourceFileRef,&
         error,failure)
#else

    CALL dsygst(itype,'U',n,a(1,1),n,b(1,1),n,info)

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END  SUBROUTINE cp_fm_cholesky_reduce

!******************************************************************
!
! op can be "SOLVE" (out = U^-1 * in ) or "MULTIPLY"   (out = U * in )
! pos can be "LEFT" or "RIGHT" (U at the left or at the right)
!
! DEPRECATED, see cp_fm_basic_linalg:cp_fm_triangular_multiply
!
!******************************************************************
  SUBROUTINE cp_fm_cholesky_restore(matrix,neig,matrixb,matrixout,op,pos,&
       error)
    TYPE(cp_fm_type), POINTER          :: matrix,matrixb,matrixout
    INTEGER, INTENT(IN)                         :: neig
    CHARACTER ( LEN = * ), INTENT ( IN )        :: op
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: pos
    TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_restore',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: a,b,out
    INTEGER                                   :: info,itype,handle
    INTEGER                                   :: n
    LOGICAL                                   :: ionode, failure
    REAL(KIND = dp)                                  :: alpha
    INTEGER                                   :: myprow, mypcol, i
    INTEGER, DIMENSION(9)                     :: desca,descb,descout
    TYPE(cp_blacs_env_type), POINTER          :: context
    CHARACTER                                 :: chol_pos

    CALL timeset("cp_fm_cholesky_restore","I","",handle)

    failure=.FALSE.
    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = matrix%matrix_struct%nrow_global
    itype = 1
    CALL cp_assert(op == "SOLVE" .OR. op == "MULTIPLY",&
         cp_failure_level,cp_assertion_failed,routineP,&
         "wrong argument op"//&
         CPSourceFileRef,&
         error,failure)

    IF (PRESENT(pos)) THEN 
       SELECT CASE(pos)
       CASE("LEFT")
         chol_pos='L'
       CASE("RIGHT")
         chol_pos='R'
       CASE DEFAULT
          CALL cp_assert(.FALSE.,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "wrong argument pos"//&
               CPSourceFileRef,&
               error,failure)
       END SELECT
    ELSE
       chol_pos='L'
    ENDIF

    ! notice b is the cholesky guy
    a => matrix%local_data
    b => matrixb%local_data
    out => matrixout%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix%matrix_struct%descriptor(:)
    descb(:) = matrixb%matrix_struct%descriptor(:)
    descout(:) = matrixout%matrix_struct%descriptor(:)
    alpha=1.0_dp
    DO i=1,neig
      CALL pdcopy(n,a(1,1),1,i,desca,1,out(1,1),1,i,descout,1)
    ENDDO
    IF (op.EQ."SOLVE") THEN
      CALL pdtrsm(chol_pos,'U','N','N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
    ELSE
      CALL pdtrmm(chol_pos,'U','N','N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
    ENDIF
#else

    alpha=1.0_dp
    CALL dcopy(neig*n,a(1,1),1,out(1,1),1)
    IF (op.EQ."SOLVE") THEN
      CALL dtrsm(chol_pos,'U','N','N',n,neig,alpha,b(1,1),n,out(1,1),n)
    ELSE
      CALL dtrmm(chol_pos,'U','N','N',n,neig,alpha,b(1,1),n,out(1,1),n)
    ENDIF

#endif

    CALL timestop(0.0_dp,handle)

  END  SUBROUTINE cp_fm_cholesky_restore

END MODULE cp_fm_cholesky
