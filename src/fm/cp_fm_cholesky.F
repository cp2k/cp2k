!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief various cholesky decomposition related routines
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_cholesky
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE kinds,                           ONLY: dp,&
                                             sp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_cholesky'

  PUBLIC :: cp_fm_cholesky_decompose, cp_fm_cholesky_invert,&
       cp_fm_cholesky_reduce, cp_fm_cholesky_restore

!***
CONTAINS

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
  SUBROUTINE cp_fm_cholesky_decompose(matrix,n,info_out,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n
    INTEGER, INTENT(out), OPTIONAL           :: info_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_cholesky_decompose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, my_n
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca
#endif

    CALL timeset(routineN,handle)

    failure = .FALSE.

    my_n = MIN(matrix%matrix_struct%nrow_global,&
         matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       CPPrecondition(n<=my_n,cp_failure_level,routineP,error,failure)
       my_n=n
    END IF

    a => matrix%local_data
    a_sp => matrix%local_data_sp

#if defined(__SCALAPACK)
    desca(:) = matrix%matrix_struct%descriptor(:)

    IF(matrix%use_sp) THEN
       CALL pspotrf('U',my_n,a_sp(1,1),1,1,desca,info)
    ELSE
       CALL pdpotrf('U',my_n,a(1,1),1,1,desca,info)
    ENDIF

#else

    IF(matrix%use_sp) THEN
       CALL spotrf('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotrf('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF

#endif

    IF (PRESENT(info_out)) THEN
      info_out = info
    ELSE
      CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_decompose

! *****************************************************************************
!> \brief used to replace the cholesky decomposition by the inverse
!> \param matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert (defaults to the min(size(matrix)))
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_invert(matrix,n, error)
    TYPE(cp_fm_type), POINTER           :: matrix
    INTEGER, INTENT(in), OPTIONAL                :: n
    TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_invert',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a
    REAL(KIND = sp), DIMENSION(:,:), POINTER  :: a_sp
    INTEGER                                   :: info,handle
    INTEGER                                   :: my_n
    LOGICAL                                   :: failure
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                     :: desca
#endif

    CALL timeset(routineN,handle)

    failure = .FALSE.

    my_n = MIN(matrix%matrix_struct%nrow_global,&
         matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       CPPrecondition(n<=my_n,cp_failure_level,routineP,error,failure)
       my_n=n
    END IF

    a => matrix%local_data
    a_sp => matrix%local_data_sp

#if defined(__SCALAPACK)

    desca(:) = matrix%matrix_struct%descriptor(:)

    IF(matrix%use_sp) THEN
       CALL pspotri('U',my_n,a_sp(1,1),1,1,desca,info)
    ELSE
       CALL pdpotri('U',my_n,a(1,1),1,1,desca,info)
    ENDIF

#else

    IF(matrix%use_sp) THEN
       CALL spotri('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotri('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_invert

! *****************************************************************************
!> \brief reduce a matrix pencil A,B to normal form
!>      B has to be cholesky decomposed with  cp_fm_cholesky_decompose
!>      before calling this routine
!>      A,B -> inv(U^T)*A*inv(U),1
!>      (AX=BX -> inv(U^T)*A*inv(U)*U*X=U*X hence evecs U*X)
!> \param matrix the symmetric matrix A
!> \param matrixb the cholesky decomposition of matrix B
!> \param error variable to control error logging, stopping,...
!>              see module cp_error_handling
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_reduce(matrix,matrixb, itype, error)
   TYPE(cp_fm_type), POINTER           :: matrix, matrixb
   INTEGER, OPTIONAL                   :: itype
   TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_reduce',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a, b
    INTEGER                                   :: info, handle
    INTEGER                                   :: n, my_itype
    LOGICAL                                   :: failure
#if defined(__SCALAPACK)
    REAL(KIND=dp)                             :: scale
    INTEGER, DIMENSION(9)                     :: desca, descb
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.
    n = matrix%matrix_struct%nrow_global

    my_itype =1
    IF( PRESENT(itype) ) my_itype = itype

    a => matrix%local_data
    b => matrixb%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix%matrix_struct%descriptor(:)
    descb(:) = matrixb%matrix_struct%descriptor(:)

    CALL pdsygst(my_itype,'U',n,a(1,1),1,1,desca,b(1,1),1,1,descb,scale,info)

    ! this is supposed to be one in current version of lapack
    ! if not, eigenvalues have to be scaled by this number
    CALL cp_assert(scale == 1.0_dp, cp_failure_level,cp_assertion_failed,&
         routineP,&
         "scale not equal 1 (scale="//cp_to_string(scale)//")"//&
CPSourceFileRef,&
         error,failure)
#else

    CALL dsygst(my_itype,'U',n,a(1,1),n,b(1,1),n,info)

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_reduce

!
! op can be "SOLVE" (out = U^-1 * in ) or "MULTIPLY"   (out = U * in )
! pos can be "LEFT" or "RIGHT" (U at the left or at the right)
!
! DEPRECATED, see cp_fm_basic_linalg:cp_fm_triangular_multiply
!
! *****************************************************************************
  SUBROUTINE cp_fm_cholesky_restore(matrix,neig,matrixb,matrixout,op,pos,transa,&
       error)
    TYPE(cp_fm_type), POINTER          :: matrix,matrixb,matrixout
    INTEGER, INTENT(IN)                         :: neig
    CHARACTER ( LEN = * ), INTENT ( IN )        :: op
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: pos
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: transa
    TYPE(cp_error_type), INTENT(inout)  :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_cholesky_restore',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: a,b,out
    REAL(KIND = sp), DIMENSION(:,:), POINTER         :: a_sp,b_sp,out_sp
    INTEGER                                   :: itype, handle
    INTEGER                                   :: n
    LOGICAL                                   :: failure
    REAL(KIND = dp)                           :: alpha
    INTEGER                                   :: myprow, mypcol
    TYPE(cp_blacs_env_type), POINTER          :: context
    CHARACTER                                 :: chol_pos,chol_transa
#if defined(__SCALAPACK)
    INTEGER                                   :: i
    INTEGER, DIMENSION(9)                     :: desca, descb, descout
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.
    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
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

    chol_transa='N'
    IF (PRESENT(transa)) chol_transa=transa

    CALL cp_assert((matrix%use_sp.EQV.matrixb%use_sp).AND.(matrix%use_sp.EQV.matrixout%use_sp), &
         cp_failure_level,cp_assertion_failed,&
         routineP,&
         "not the same precision"//&
CPSourceFileRef,&
         error,failure)

    ! notice b is the cholesky guy
    a => matrix%local_data
    b => matrixb%local_data
    out => matrixout%local_data
    a_sp => matrix%local_data_sp
    b_sp => matrixb%local_data_sp
    out_sp => matrixout%local_data_sp

#if defined(__SCALAPACK)

    desca(:) = matrix%matrix_struct%descriptor(:)
    descb(:) = matrixb%matrix_struct%descriptor(:)
    descout(:) = matrixout%matrix_struct%descriptor(:)
    alpha=1.0_dp
    DO i=1,neig
       IF(matrix%use_sp) THEN
          CALL pscopy(n,a_sp(1,1),1,i,desca,1,out_sp(1,1),1,i,descout,1)
       ELSE
          CALL pdcopy(n,a(1,1),1,i,desca,1,out(1,1),1,i,descout,1)
       ENDIF
    ENDDO
    IF (op.EQ."SOLVE") THEN
       IF(matrix%use_sp) THEN
          CALL pstrsm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),1,1,descb,&
                      out_sp(1,1),1,1,descout)
       ELSE
          CALL pdtrsm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
       ENDIF
    ELSE
       IF(matrix%use_sp) THEN
          CALL pstrmm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),1,1,descb,&
                      out_sp(1,1),1,1,descout)
       ELSE
          CALL pdtrmm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
       ENDIF
    ENDIF
#else

    alpha=1.0_dp
    IF(matrix%use_sp) THEN
       CALL scopy(neig*n,a_sp(1,1),1,out_sp(1,1),1)
    ELSE
       CALL dcopy(neig*n,a(1,1),1,out(1,1),1)
    ENDIF
    IF (op.EQ."SOLVE") THEN
       IF(matrix%use_sp) THEN
          CALL strsm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),SIZE(b_sp,1),out_sp(1,1),n)
       ELSE
          CALL dtrsm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),SIZE(b,1),out(1,1),n)
       ENDIF
    ELSE
       IF(matrix%use_sp) THEN
          CALL strmm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),n,out_sp(1,1),n)
       ELSE
          CALL dtrmm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),n,out(1,1),n)
       ENDIF
    ENDIF

#endif

    CALL timestop(handle)

  END  SUBROUTINE cp_fm_cholesky_restore

END MODULE cp_fm_cholesky
