!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Interface to (sca)lapack for the Cholesky based procedures
!> \author  VW
!> \date    2009-09-08
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-09-08
! *****************************************************************************
MODULE cp_dbcsr_cholesky
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_info,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             sp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_cholesky'

  PUBLIC :: cp_dbcsr_cholesky_decompose, cp_dbcsr_cholesky_invert,&
            cp_dbcsr_cholesky_restore

  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief used to replace a symmetric positive def. matrix M with its cholesky
!>      decomposition U: M = U^T * U, with U upper triangular
!> \param matrix the matrix to replace with its cholesky decomposition
!> \param n the number of row (and columns) of the matrix &
!>        (defaults to the min(size(matrix)))
!> \param para_env ...
!> \param blacs_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!> \author Joost
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_decompose(matrix,n,para_env,blacs_env,error)
    TYPE(cp_dbcsr_type)                      :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: n
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_cholesky_decompose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, my_n, &
                                                nfullcols_total, &
                                                nfullrows_total
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: a_sp
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_matrix
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(fm_matrix, fm_struct)
    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)

    my_n = MIN(fm_matrix%matrix_struct%nrow_global,&
         fm_matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       CPPrecondition(n<=my_n,cp_failure_level,routineP,error,failure)
       my_n=n
    END IF

    a => fm_matrix%local_data
    a_sp => fm_matrix%local_data_sp

#if defined(__SCALAPACK)
    desca(:) = fm_matrix%matrix_struct%descriptor(:)

    IF(fm_matrix%use_sp) THEN
       CALL pspotrf('U',my_n,a_sp(1,1),1,1,desca,info)
    ELSE
       CALL pdpotrf('U',my_n,a(1,1),1,1,desca,info)
    ENDIF

#else

    IF(fm_matrix%use_sp) THEN
       CALL spotrf('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotrf('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    CALL copy_fm_to_dbcsr(fm_matrix, matrix, error=error)

    CALL cp_fm_release(fm_matrix, error)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_decompose

! *****************************************************************************
!> \brief used to replace the cholesky decomposition by the inverse
!> \param matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert (defaults to the min(size(matrix)))
!> \param para_env ...
!> \param blacs_env ...
!> \param upper_to_full ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_invert(matrix,n,para_env,blacs_env,upper_to_full,error)
    TYPE(cp_dbcsr_type)                           :: matrix
    INTEGER, INTENT(in), OPTIONAL             :: n
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(cp_blacs_env_type), POINTER          :: blacs_env
    LOGICAL, INTENT(IN)                       :: upper_to_full
    TYPE(cp_error_type), INTENT(inout)        :: error

    CHARACTER(len=*), PARAMETER :: routineN='dbcsr_cholesky_invert',&
         routineP=moduleN//':'//routineN

    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a
    REAL(KIND = sp), DIMENSION(:,:), POINTER  :: a_sp
    INTEGER                                   :: info,handle
    INTEGER                                   :: my_n, nfullrows_total, nfullcols_total
    LOGICAL                                   :: failure
    TYPE(cp_fm_type), POINTER                 :: fm_matrix, fm_matrix_tmp
    TYPE(cp_fm_struct_type), POINTER          :: fm_struct
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                     :: desca
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(fm_matrix, fm_struct)
    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullrows_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)

    my_n = MIN(fm_matrix%matrix_struct%nrow_global,&
         fm_matrix%matrix_struct%ncol_global)
    IF (PRESENT(n)) THEN
       CPPrecondition(n<=my_n,cp_failure_level,routineP,error,failure)
       my_n=n
    END IF

    a => fm_matrix%local_data
    a_sp => fm_matrix%local_data_sp

#if defined(__SCALAPACK)

    desca(:) = fm_matrix%matrix_struct%descriptor(:)

    IF(fm_matrix%use_sp) THEN
       CALL pspotri('U',my_n,a_sp(1,1),1,1,desca,info)
    ELSE
       CALL pdpotri('U',my_n,a(1,1),1,1,desca,info)
    ENDIF

#else

    IF(fm_matrix%use_sp) THEN
       CALL spotri('U',my_n,a_sp(1,1),SIZE(a_sp,1),info)
    ELSE
       CALL dpotri('U',my_n,a(1,1),SIZE(a,1),info)
    ENDIF

#endif

    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    IF(upper_to_full) THEN
       CALL cp_fm_create(fm_matrix_tmp,fm_matrix%matrix_struct,name="fm_matrix_tmp",error=error)
       CALL cp_fm_upper_to_full(fm_matrix, fm_matrix_tmp, error)
       CALL cp_fm_release(fm_matrix_tmp, error)
    ENDIF

    CALL copy_fm_to_dbcsr(fm_matrix, matrix, error=error)

    CALL cp_fm_release(fm_matrix, error)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_invert

! *****************************************************************************
!> \brief reduce a matrix pencil A,B to normal form
!>      B has to be cholesky decomposed with  cp_dbcsr_cholesky_decompose
!>      before calling this routine
!>      A,B -> inv(U^T)*A*inv(U),1
!>      (AX=BX -> inv(U^T)*A*inv(U)*U*X=U*X hence evecs U*X)
!> \param matrix the symmetric matrix A
!> \param matrixb the cholesky decomposition of matrix B
!> \param para_env ...
!> \param blacs_env ...
!> \param error variable to control error logging, stopping,...
!>              see module cp_error_handling
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_reduce(matrix,matrixb,para_env,blacs_env,error)
    TYPE(cp_dbcsr_type)                           :: matrix,matrixb
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(cp_blacs_env_type), POINTER          :: blacs_env
    TYPE(cp_error_type), INTENT(inout)        :: error

    CHARACTER(len=*), PARAMETER :: routineN='dbcsr_cholesky_reduce',&
         routineP=moduleN//':'//routineN

    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a,b
    INTEGER                                   :: info,itype,handle
    INTEGER                                   :: n, nfullrows_total, nfullcols_total
    LOGICAL                                   :: failure
    TYPE(cp_fm_type), POINTER                 :: fm_matrix, fm_matrixb
    TYPE(cp_fm_struct_type), POINTER          :: fm_struct
#if defined(__SCALAPACK)
    REAL(KIND=dp)                             :: scale
    INTEGER, DIMENSION(9)                     :: desca, descb
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(fm_matrix, fm_struct)
    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)

    CALL cp_dbcsr_get_info(matrixb,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)

    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrixb,fm_struct,name="fm_matrixb",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrixb, fm_matrixb, error=error)

    n = fm_matrix%matrix_struct%nrow_global
    itype =1

    a => fm_matrix%local_data
    b => fm_matrixb%local_data

#if defined(__SCALAPACK)

    desca(:) = fm_matrix%matrix_struct%descriptor(:)
    descb(:) = fm_matrixb%matrix_struct%descriptor(:)

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

    CALL copy_fm_to_dbcsr(fm_matrix, matrix, error=error)

    CALL cp_fm_release(fm_matrix, error)
    CALL cp_fm_release(fm_matrixb, error)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_reduce

!
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param neig ...
!> \param matrixb ...
!> \param matrixout ...
!> \param op ...
!> \param pos ...
!> \param transa ...
!> \param para_env ...
!> \param blacs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_cholesky_restore(matrix,neig,matrixb,matrixout,op,pos,transa,&
       para_env,blacs_env,error)
    TYPE(cp_dbcsr_type)                                :: matrix,matrixb,matrixout
    INTEGER, INTENT(IN)                            :: neig
    CHARACTER ( LEN = * ), INTENT ( IN )           :: op
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: pos
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: transa
    TYPE(cp_para_env_type), POINTER                :: para_env
    TYPE(cp_blacs_env_type), POINTER               :: blacs_env
    TYPE(cp_error_type), INTENT(inout)             :: error

    CHARACTER(len=*), PARAMETER :: routineN='dbcsr_cholesky_restore',&
         routineP=moduleN//':'//routineN

    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: a,b,out
    REAL(KIND = sp), DIMENSION(:,:), POINTER  :: a_sp,b_sp,out_sp
    INTEGER                                   :: itype,handle
    INTEGER                                   :: n
    LOGICAL                                   :: failure
    REAL(KIND = dp)                           :: alpha
    INTEGER                                   :: myprow, mypcol, nfullrows_total, &
         nfullcols_total
    TYPE(cp_blacs_env_type), POINTER          :: context
    CHARACTER                                 :: chol_pos,chol_transa
    TYPE(cp_fm_type), POINTER                 :: fm_matrix,fm_matrixb,fm_matrixout
    TYPE(cp_fm_struct_type), POINTER          :: fm_struct
#if defined(__SCALAPACK)
    INTEGER                                   :: i
    INTEGER, DIMENSION(9)                     :: desca, descb, descout
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(fm_matrix, fm_matrixb, fm_matrixout, fm_struct)

    CALL cp_dbcsr_get_info(matrix,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix,fm_struct,name="fm_matrix",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_dbcsr_get_info(matrixb,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrixb,fm_struct,name="fm_matrixb",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_dbcsr_get_info(matrixout,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrixout,fm_struct,name="fm_matrixout",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix, fm_matrix, error=error)
    CALL copy_dbcsr_to_fm(matrixb, fm_matrixb, error=error)
    !CALL copy_dbcsr_to_fm(matrixout, fm_matrixout)

    context => fm_matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = fm_matrix%matrix_struct%nrow_global
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

    CALL cp_assert((fm_matrix%use_sp.EQV.fm_matrixb%use_sp).AND.(fm_matrix%use_sp.EQV.fm_matrixout%use_sp), &
         cp_failure_level,cp_assertion_failed,&
         routineP,&
         "not the same precision"//&
CPSourceFileRef,&
         error,failure)

    ! notice b is the cholesky guy
    a => fm_matrix%local_data
    b => fm_matrixb%local_data
    out => fm_matrixout%local_data
    a_sp => fm_matrix%local_data_sp
    b_sp => fm_matrixb%local_data_sp
    out_sp => fm_matrixout%local_data_sp

#if defined(__SCALAPACK)

    desca(:) = fm_matrix%matrix_struct%descriptor(:)
    descb(:) = fm_matrixb%matrix_struct%descriptor(:)
    descout(:) = fm_matrixout%matrix_struct%descriptor(:)
    alpha=1.0_dp
    DO i=1,neig
       IF(fm_matrix%use_sp) THEN
          CALL pscopy(n,a_sp(1,1),1,i,desca,1,out_sp(1,1),1,i,descout,1)
       ELSE
          CALL pdcopy(n,a(1,1),1,i,desca,1,out(1,1),1,i,descout,1)
       ENDIF
    ENDDO
    IF (op.EQ."SOLVE") THEN
       IF(fm_matrix%use_sp) THEN
          CALL pstrsm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),1,1,descb,&
                      out_sp(1,1),1,1,descout)
       ELSE
          CALL pdtrsm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
       ENDIF
    ELSE
       IF(fm_matrix%use_sp) THEN
          CALL pstrmm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),1,1,descb,&
                      out_sp(1,1),1,1,descout)
       ELSE
          CALL pdtrmm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
       ENDIF
    ENDIF
#else

    alpha=1.0_dp
    IF(fm_matrix%use_sp) THEN
       CALL scopy(neig*n,a_sp(1,1),1,out_sp(1,1),1)
    ELSE
       CALL dcopy(neig*n,a(1,1),1,out(1,1),1)
    ENDIF
    IF (op.EQ."SOLVE") THEN
       IF(fm_matrix%use_sp) THEN
          CALL strsm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),SIZE(b_sp,1),out_sp(1,1),n)
       ELSE
          CALL dtrsm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),SIZE(b,1),out(1,1),n)
       ENDIF
    ELSE
       IF(fm_matrix%use_sp) THEN
          CALL strmm(chol_pos,'U',chol_transa,'N',n,neig,REAL(alpha,sp),b_sp(1,1),n,out_sp(1,1),n)
       ELSE
          CALL dtrmm(chol_pos,'U',chol_transa,'N',n,neig,alpha,b(1,1),n,out(1,1),n)
       ENDIF
    ENDIF

#endif

    CALL copy_fm_to_dbcsr(fm_matrixout, matrixout, error=error)

    CALL cp_fm_release(fm_matrix, error)
    CALL cp_fm_release(fm_matrixb, error)
    CALL cp_fm_release(fm_matrixout, error)

    CALL timestop(handle)

  END  SUBROUTINE cp_dbcsr_cholesky_restore

END MODULE cp_dbcsr_cholesky

