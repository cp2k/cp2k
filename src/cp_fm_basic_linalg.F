!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_from_globenv,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_equivalent,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_create2,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_set_submatrix,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_dealloc_ref,&
                                             cp_get,&
                                             cp_init,&
                                             cp_matrix_block_iterator,&
                                             cp_next
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parallel,                        ONLY: cp2k_is_parallel
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             sparse_plus_rep_rept,&
                                             sparse_times_replicated
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_basic_linalg'

  PUBLIC :: cp_fm_add, cp_sm_fm_multiply, cp_sm_plus_fm_fm_t,&
       cp_fm_gemm, cp_fm_symm, cp_fm_syrk, cp_fm_schur_product,&
       cp_fm_trace, cp_fm_triangular_multiply
!!***
!****************************************************************************
CONTAINS

  SUBROUTINE cp_fm_add(alpha,matrix_a,beta,matrix_b,error)

!   Purpose: Scale and add two BLACS matrices (a <- alpha*a + beta*b).

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_b
    REAL(wp), INTENT(IN)                      :: alpha,beta
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,mypcol,myprow,npcol,nprow,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_add',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    failure=.FALSE.

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_b),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_b%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL cp_assert(matrix_a%matrix_struct%context%group==&
         matrix_b%matrix_struct%context%group,cp_failure_level,&
         cp_assertion_failed,fromWhere=routineP,&
         message="matrixes must be in the same blacs context"//&
         CPSourceFileRef,&
         error=error,failure=failure)

    IF (cp_fm_struct_equivalent(matrix_a%matrix_struct,matrix_b%matrix_struct,&
         error=error)) THEN
! to do: use dscal,dcopy,daxp
       myprow=matrix_a%matrix_struct%context%mepos(1)
       mypcol=matrix_a%matrix_struct%context%mepos(2)

       a => matrix_a%local_data
       b => matrix_b%local_data

       IF (alpha == 0.0_wp) THEN
          IF (beta == 0.0_wp) THEN
             a(:,:) = 0.0_wp
          ELSE IF (beta == 1.0_wp) THEN
             a(:,:) = b(:,:)
          ELSE
             a(:,:) = beta*b(:,:)
          END IF
       ELSE IF (beta == 0.0_wp) THEN
          IF (alpha == 1.0_wp) THEN
             RETURN
          ELSE
             a(:,:) = alpha*a(:,:)
          END IF
       ELSE IF (alpha == 1.0_wp) THEN
          IF (beta == 1.0_wp) THEN
             a(:,:) = a(:,:) + b(:,:)
          ELSE
             a(:,:) = a(:,:) + beta*b(:,:)
          END IF
       ELSE IF (beta == 1.0_wp) THEN
          a(:,:) = alpha*a(:,:) + b(:,:)
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
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_add

! *****************************************************************************

! *****************************************************************************
  SUBROUTINE cp_sm_fm_multiply(sparse_matrix,v_in,v_out,ncol,para_env, &
       alpha,beta,block_owner,error)
! to do: better interface
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_full_matrix_type) , POINTER :: v_in
    TYPE(cp_full_matrix_type) , POINTER :: v_out
    INTEGER, INTENT(IN)               :: ncol
    TYPE(cp_para_env_type), POINTER :: para_env ! of the sparse matrix
    REAL(kind=wp), INTENT(in), OPTIONAL :: alpha, beta
    INTEGER, INTENT(in), DIMENSION(:,:), OPTIONAL :: block_owner
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***
    INTEGER, PARAMETER          :: max_blocksize = 100
    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE cp_sm_fm_multiply (MODULE qs_blacs)"

!   *** Local variables ***
    LOGICAL :: failure
    REAL(wp), DIMENSION(:,:), POINTER :: rep_v_in,rep_v_out
    INTEGER i,nrow_in,nrow_out,ncol_block,handle,blocksize, stat
!   ---------------------------------------------------------------------------
    CALL timeset("cp_sm_fm_multiply","I","",handle)

    failure=.FALSE.
    nrow_in=v_in%matrix_struct%nrow_global
    nrow_out=v_out%matrix_struct%nrow_global
    blocksize=MIN(ncol,max_blocksize)
    ALLOCATE(rep_v_in(blocksize,nrow_in),rep_v_out(blocksize,nrow_out),&
         stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)

    ! zero the rep_v_in and rep_v_out matrices
    CALL dcopy(blocksize*nrow_in, 0.0_wp, 0, rep_v_in, 1)
    CALL dcopy(blocksize*nrow_out, 0.0_wp, 0, rep_v_out, 1)

    DO i=1,ncol,blocksize
       ncol_block=MIN(ncol-i+1,blocksize)
       CALL cp_fm_get_submatrix(v_in,rep_v_in,start_col=i,n_cols=ncol_block,&
            transpose_target_m=.TRUE.)
       CALL sparse_times_replicated(sparse_matrix,rep_v_in,rep_v_out,&
            ncol_block, blocksize, para_env=para_env, block_owner=block_owner)
       CALL cp_fm_set_submatrix(v_out,rep_v_out,&
            start_col=i,n_cols=ncol_block,transpose_new_val=.TRUE.,&
            alpha=alpha,beta=beta)

    ENDDO

    DEALLOCATE(rep_v_in,rep_v_out,stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_fm_multiply

! *****************************************************************************
! computes sparse=sparse+alpha*v*g^T 
! only v*v^T tested
! *****************************************************************************
  SUBROUTINE cp_sm_plus_fm_fm_t(sparse_matrix,matrix_v,matrix_g,ncol,&
       first_col,alpha,error)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_full_matrix_type) , POINTER :: matrix_v
    TYPE(cp_full_matrix_type) , POINTER, OPTIONAL :: matrix_g
    INTEGER, INTENT(IN)               :: ncol
    REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
    INTEGER, OPTIONAL, INTENT(IN)     :: first_col
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

!   *** Local parameters ***
    INTEGER, PARAMETER          :: max_blocksize = 100
    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE cp_sm_plus_fm_fm_t (MODULE qs_blacs)"

!   *** Local variables ***
    REAL(wp), DIMENSION(:,:), POINTER :: rep_v,rep_g
    REAL(wp) :: prefactor
    INTEGER col1,i,nrow,ncol_block,handle,blocksize
    LOGICAL g_present
!   ---------------------------------------------------------------------------
    CALL timeset("cp_sm_plus_fm_fm_t","I","",handle)

    IF (PRESENT(first_col)) THEN
       col1 = first_col
    ELSE
       col1 = 1
    END IF
    blocksize=MIN(max_blocksize,ncol-col1+1)
    g_present=.FALSE.
    IF (PRESENT(matrix_g)) g_present=.TRUE.

    nrow=matrix_v%matrix_struct%nrow_global
    IF (g_present) THEN
       ALLOCATE(rep_v(blocksize,nrow))
       ALLOCATE(rep_g(blocksize,nrow))
    ELSE
       ALLOCATE(rep_v(blocksize,nrow))
       rep_g=>rep_v
    ENDIF

    IF (PRESENT(alpha)) THEN
       prefactor = alpha
    ELSE
       prefactor = 1.0_wp
    END IF


    DO i=col1,ncol,blocksize
       ncol_block=MIN(ncol-i+1,blocksize)
!CALL blacs_replicated_copy(matrix_v,rep_v,i,ncol_block, &
!                                 blocksize,"B2R")
       CALL cp_fm_get_submatrix(matrix_v,rep_v,start_col=i,n_cols=ncol_block,&
            transpose_target_m=.TRUE., error=error)
       IF (g_present) THEN
!CALL blacs_replicated_copy(matrix_g,rep_g,i,ncol_block, &
!                                 blocksize,"B2R")
          CALL cp_fm_get_submatrix(matrix_g,rep_g,start_col=i,&
               n_cols=ncol_block,&
               transpose_target_m=.TRUE., error=error)
       ENDIF
       CALL sparse_plus_rep_rept(sparse_matrix,rep_v,rep_g,ncol_block, &
            blocksize,alpha=prefactor)
    ENDDO

    IF (g_present) THEN
       DEALLOCATE(rep_v,rep_g)
    ELSE
       DEALLOCATE(rep_v)
    ENDIF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_plus_fm_fm_t
! *****************************************************************************


! *****************************************************************************

  SUBROUTINE cp_fm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error)

!   Purpose: BLACS interface to the BLAS routine dgemm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    CHARACTER(LEN=1), INTENT(IN)              :: transa,transb
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: k,m,n
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,lda,ldb,ldc

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c
    REAL(wp) :: flops

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_gemm","I","Gflops",handle)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = matrix_b%matrix_struct%descriptor(:)
    descc(:) = matrix_c%matrix_struct%descriptor(:)

    flops=2.0_wp*(1E-9_wp*m*n*k)/(matrix_a%matrix_struct%context%num_pe(1)*&
         matrix_a%matrix_struct%context%num_pe(2))
    CALL pdgemm(transa,transb,m,n,k,alpha,a(1,1),1,1,desca,b(1,1),1,1, &
         descb,beta,c(1,1),1,1,&
         descc)
#else

    lda = SIZE(a,1)
    ldb = SIZE(b,1)
    ldc = SIZE(c,1)

    flops=2.0_wp*(1E-9_wp*m*n*k)
    CALL dgemm(transa,transb,m,n,k,alpha,a(1,1),lda,b(1,1),ldb,beta,c(1,1),ldc)

#endif
    CALL timestop(flops,handle)

  END SUBROUTINE cp_fm_gemm

  SUBROUTINE cp_fm_symm(side,uplo,m,n,alpha,matrix_a,matrix_b,beta,matrix_c,&
       error)

!   Purpose: BLACS interface to the BLAS routine dsymm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    CHARACTER(LEN=1), INTENT(IN)              :: side,uplo
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: m,n
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,lda,ldb,ldc,mypcol,myprow,npcol,nprow

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c
    TYPE(cp_blacs_env_type), POINTER    :: context

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_symm","I","",handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    a => matrix_a%local_data
    b => matrix_b%local_data
    c => matrix_c%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix_a%matrix_struct%descriptor(:)
    descb(:) = matrix_b%matrix_struct%descriptor(:)
    descc(:) = matrix_c%matrix_struct%descriptor(:)

    CALL pdsymm(side,uplo,m,n,alpha,a(1,1),1,1,desca,b(1,1),1,1,descb,beta,c(1,1),1,1,descc)

#else

    lda = matrix_a%matrix_struct%nrow_global
    ldb = matrix_b%matrix_struct%nrow_global
    ldc = matrix_c%matrix_struct%nrow_global

    CALL dsymm(side,uplo,m,n,alpha,a(1,1),lda,b(1,1),ldb,beta,c(1,1),ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_symm

! *****************************************************************************

  SUBROUTINE cp_fm_syrk(uplo,trans,k,alpha,matrix_a,ia,ja,beta,matrix_c,&
       error)

!   Purpose: BLACS interface to the BLAS routine dsyrk.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_c
    CHARACTER(LEN=1), INTENT(IN)              :: trans,uplo
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: ia,ja,k
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,lda,ldc,mypcol,myprow,n,npcol,nprow

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c
    TYPE(cp_blacs_env_type), POINTER    :: context

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

    lda = matrix_a%matrix_struct%nrow_global
    ldc = matrix_c%matrix_struct%nrow_global

    CALL dsyrk(uplo,trans,n,k,alpha,a(ia,ja),lda,beta,c(1,1),ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_syrk

! *****************************************************************************

  SUBROUTINE cp_fm_schur_product(matrix_a,matrix_b,matrix_c,error)

!   Purpose: Calculate the schur product of two matrices

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,icol_local,irow_local,mypcol,myprow,&
         ncol_local,npcol,nprow,nrow_local

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c
    TYPE(cp_blacs_env_type), POINTER    :: context

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

    nrow_local = matrix_a%matrix_struct%nrow_local(myprow)
    ncol_local = matrix_a%matrix_struct%ncol_local(mypcol)

    DO icol_local=1,ncol_local
       DO irow_local=1,nrow_local
          c(irow_local,icol_local) = a(irow_local,icol_local)*b(irow_local,icol_local)
       END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_schur_product

! *****************************************************************************
!!****f* cp_fm_basic_linalg/cp_fm_trace [1.0] *
!!
!!   NAME
!!     cp_fm_trace
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_trace(matrix_a, matrix_b, trace, error)
!!       Type(cp_full_matrix_type), Pointer:: matrix_a, matrix_b
!!       Real(Kind=wp), Intent (OUT):: trace
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

    TYPE(cp_full_matrix_type), POINTER           :: matrix_a,matrix_b
    REAL(wp), INTENT(OUT)                        :: trace
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: group,handle,icol_local,irow_local,mypcol,myprow,&
         ncol_local,npcol,nprow,nrow_local

    REAL(wp), DIMENSION(:,:), POINTER :: a,b
    TYPE(cp_blacs_env_type), POINTER    :: context

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

    nrow_local = matrix_a%matrix_struct%nrow_local(myprow)
    ncol_local = matrix_b%matrix_struct%ncol_local(mypcol)

    trace = 0.0_wp

    DO icol_local=1,ncol_local
       DO irow_local=1,nrow_local
          trace = trace + a(irow_local,icol_local)*b(irow_local,icol_local)
       END DO
    END DO

    CALL mp_sum(trace,group)

    CALL timestop(0.0_wp,handle)

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
!!       Type(cp_full_matrix_type), Pointer:: triangular_matrix, matrix_b
!!       Character, Intent (IN), Optional:: side, uplo_tr
!!       Logical, Intent (IN), Optional:: transpose_tr, invert_tr,&
!!         unit_diag_tr
!!       Integer, Intent (IN), Optional:: n_rows, n_cols
!!       Real(Kind=wp), Intent (IN), Optional:: alpha
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
  TYPE(cp_full_matrix_type), POINTER :: triangular_matrix, matrix_b
  CHARACTER, INTENT(in), OPTIONAL :: side, uplo_tr
  LOGICAL, INTENT(in), OPTIONAL :: transpose_tr, invert_tr, unit_diag_tr
  INTEGER, INTENT(in), OPTIONAL :: n_rows, n_cols
  REAL(kind=wp), INTENT(in), OPTIONAL :: alpha
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, invert
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_triangular_multiply',&
        routineP=moduleN//':'//routineN
  CHARACTER :: side_char, unit_diag,uplo, transa
  INTEGER :: handle, m,n
  REAL(kind=wp) :: al
  failure=.FALSE.
  
  CALL timeset(routineN//','//moduleN,"I","",handle)
  side_char='L'
  unit_diag='N'
  uplo='U'
  transa='N'
  invert=.FALSE.
  al=1.0_wp
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

  CALL timestop(0.0_wp,handle)     
  END SUBROUTINE cp_fm_triangular_multiply
!***************************************************************************

END MODULE cp_fm_basic_linalg
