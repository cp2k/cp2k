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
                                             cp_fm_get_element,&
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
                                             mp_sync, &
                                             mp_shift
  USE parallel,                        ONLY: cp2k_is_parallel
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             sparse_plus_loc_loct,&
                                             sparse_times_local
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_basic_linalg'

  PUBLIC :: cp_fm_scale_and_add, cp_sm_fm_multiply, cp_sm_plus_fm_fm_t,&
       cp_fm_gemm, cp_fm_symm, cp_fm_syrk, cp_fm_schur_product,&
       cp_fm_trace, cp_fm_triangular_multiply, cp_fm_lu_decompose, &
       cp_fm_scale
!!***
!****************************************************************************
CONTAINS

  SUBROUTINE cp_fm_scale_and_add(alpha,matrix_a,beta,matrix_b,error)

    !   Purpose: Scale and add two BLACS matrices (a <- alpha*a + beta*b).

    !   History: - Creation (11.06.2001, Matthias Krack)

    !   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER           :: matrix_a
    REAL(wp), INTENT(IN)                         :: alpha
    REAL(wp), INTENT(in), OPTIONAL               :: beta
    TYPE(cp_full_matrix_type), POINTER, OPTIONAL :: matrix_b
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    !   *** Local variables ***

    INTEGER :: handle,mypcol,myprow,npcol,nprow,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_scale_and_add',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    REAL(wp) :: my_beta

    !   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    my_beta=0.0_wp
    IF(PRESENT(beta)) my_beta=beta
    NULLIFY(a,b)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)
    ! to do: use dscal,dcopy,daxp
    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)

    a => matrix_a%local_data

    IF (my_beta==0.0_wp) THEN

       IF (alpha==0.0_wp) THEN
          a(:,:)=0.0_wp
       ELSE IF (alpha == 1.0_wp) THEN
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

          IF (alpha == 0.0_wp) THEN
             IF (my_beta == 1.0_wp) THEN
                a(:,:) = b(:,:)
             ELSE
                a(:,:) = my_beta*b(:,:)
             END IF
          ELSE IF (alpha == 1.0_wp) THEN
             IF (my_beta == 1.0_wp) THEN
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
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_scale_and_add

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

    CALL  cp_sm_fm_multiply_new(sparse_matrix,v_in,v_out,ncol,para_env,alpha,beta,error)
    RETURN

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

       CALL dcopy(nrow_out*blocksize,0.0_wp,0,rep_v_out,1)

       CALL sparse_times_local(sparse_matrix,rep_v_in,rep_v_out,&
            ncol_block, blocksize, para_env=para_env, block_owner=block_owner)


#if defined(__parallel)
       CALL mp_sum(rep_v_out,para_env%group)
#endif

       CALL cp_fm_set_submatrix(v_out,rep_v_out,&
            start_col=i,n_cols=ncol_block,transpose_new_val=.TRUE.,&
            alpha=alpha,beta=beta)

    ENDDO

    DEALLOCATE(rep_v_in,rep_v_out,stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_fm_multiply

!******************************************************************************
!
! new version of the cp_sm_fm_multiply, optimized communication, 
!
!*******************************************************************************
  SUBROUTINE cp_sm_fm_multiply_new(sparse_matrix,v_in,v_out,ncol,para_env, &
       alpha,beta,error)
! to do: better interface
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_full_matrix_type) , POINTER :: v_in
    TYPE(cp_full_matrix_type) , POINTER :: v_out
    INTEGER, INTENT(IN)               :: ncol
    TYPE(cp_para_env_type), POINTER :: para_env ! of the sparse matrix
    REAL(kind=wp), INTENT(in), OPTIONAL :: alpha, beta
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine = "cp_sm_fm_multiply_new"

!   *** Local variables ***
    LOGICAL :: failure
    REAL(wp), DIMENSION(:,:), POINTER :: local_v_in, local_v_out
    INTEGER :: i,j,info,nrow_global,handle,num_pe,mepos,ncol_max,stat,ishift
    INTEGER, DIMENSION(:), POINTER :: ncol_loc
    INTEGER :: ictxt_loc, im,in, myprow,mypcol,nprow,npcol,desc(9),rb,cb
    REAL(wp), POINTER, DIMENSION(:,:) :: newdat
    REAL(wp) :: alpha_l,beta_l
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: NUMROC
#endif
    
    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_sm_fm_multiply_new", &
                                   routineP = moduleN//"/"//routineN

!   ---------------------------------------------------------------------------
    failure=.false.
#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"cp_sm_fm_multiply would be most happy to use scalapack libraries"//&
     CPSourceFileRef,&
     error)
#endif
    CALL timeset("cp_sm_fm_multiply_new","I","",handle)

    CPPrecondition(ASSOCIATED(sparse_matrix),cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(v_in),cp_failure_level,routine,error,failure)
    CPPrecondition(v_in%ref_count>0,cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(v_out),cp_failure_level,routine,error,failure)
    CPPrecondition(v_out%ref_count>0,cp_failure_level,routine,error,failure)

    IF (PRESENT(alpha)) THEN
       alpha_l=alpha
    ELSE
       alpha_l=1.0_wp
    ENDIF
    IF (PRESENT(beta)) THEN
       beta_l=beta
    ELSE
       beta_l=0.0_wp
    ENDIF

    failure=.FALSE.
    NULLIFY(local_v_in,local_v_out,newdat)

    CALL cp_fm_get_info(v_in,nrow_global=nrow_global)
    num_pe=para_env%num_pe
    mepos =para_env%mepos
    rb=nrow_global
    cb=1
    im=nrow_global
    ALLOCATE(ncol_loc(0:num_pe-1))
    ncol_loc=0
#if defined(__SCALAPACK)
    CALL blacs_get(-1,0,ictxt_loc)
    CALL blacs_gridinit(ictxt_loc,'R',1,num_pe)
    CALL blacs_gridinfo(ictxt_loc,nprow,npcol,myprow,mypcol)
    CALL descinit(desc,nrow_global,ncol,rb,cb,0,0,ictxt_loc,nrow_global,info)
    CPPostcondition(info==0,cp_failure_level,routine,error,failure)
    in=numroc(ncol,cb,mypcol,0,npcol)

    allocate(newdat(im,MAX(1,in)))

    ! do the actual scalapack to cols reordering
    CALL pdgemr2d(nrow_global,ncol,v_in%local_data(1,1),1,1,&
                  v_in%matrix_struct%descriptor, &
                  newdat(1,1),1,1,desc,ictxt_loc)
 
    ! obtain data in transposed structures
    ncol_loc(mepos)=in
    CALL mp_sum(ncol_loc,para_env%group)
    ncol_max=MAXVAL(ncol_loc)

    ALLOCATE(local_v_in(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    ALLOCATE(local_v_out(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    DO j=1,in
       DO i=1,im
              local_v_in(j,i)=newdat(i,j)
       END DO
    END DO
#else
    in=ncol
    ncol_loc(mepos)=in
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(local_v_in(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    ALLOCATE(local_v_out(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    DO j=1,in
       DO i=1,im
              local_v_in(j,i)=v_in%local_data(i,j)
       END DO
    END DO
#endif

    IF (beta_l.EQ.0.0_wp) THEN
       CALL dcopy(SIZE(local_v_out,1)*SIZE(local_v_out,2),&
            0.0_wp,0,local_v_out(1,1),1) ! zero output vector
    ELSE ! this rather expensive and should probably not be used. Anyway, not tested
#if defined(__SCALAPACK)
         CALL pdgemr2d(nrow_global,ncol,v_out%local_data(1,1),1,1,&
                       v_out%matrix_struct%descriptor, &
                       newdat(1,1),1,1,desc,ictxt_loc)
         DO j=1,in
            DO i=1,im
                local_v_out(j,i)=beta_l*newdat(i,j)
            END DO
         END DO
#else
       DO j=1,in
          DO i=1,im
                 local_v_out(j,i)=beta_l*v_out%local_data(i,j)
          END DO
       END DO
#endif
    ENDIF

    ! now we shift around the data in a ring, multiplying the v_in with the local sparse matrix
    ! adding the results to v_out, we need to do num_pe shifts

    do ishift=0,num_pe-1
       CALL sparse_times_local(sparse_matrix,local_v_in,local_v_out,&
            ncol_loc(MOD(mepos-ishift+num_pe,num_pe)), ncol_max, para_env=para_env,alpha=alpha_l)
       CALL mp_shift(local_v_in,para_env%group)
       CALL mp_shift(local_v_out,para_env%group)
    enddo

#if defined(__SCALAPACK)
    ! copy the result back into the scalapack data structure
    DO j=1,in
       DO i=1,im
              newdat(i,j)=local_v_out(j,i)
       END DO
    END DO
    ! shuffle around the data
    CALL pdgemr2d(nrow_global,ncol, &
                  newdat(1,1),1,1,desc, &
                  v_out%local_data(1,1),1,1,v_out%matrix_struct%descriptor, &
                  ictxt_loc)
    ! give the grid back
    CALL blacs_gridexit(ictxt_loc)
    DEALLOCATE(newdat)
#else
    DO j=1,in
       DO i=1,im
              v_out%local_data(i,j)=local_v_out(j,i)
       END DO
    END DO
#endif

    DEALLOCATE(local_v_in,local_v_out,ncol_loc)

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_fm_multiply_new

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
    CALL cp_sm_plus_fm_fm_t_new(sparse_matrix,matrix_v,matrix_g,ncol, &
                                first_col,alpha,error)
    RETURN
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
       CALL cp_fm_get_submatrix(matrix_v,rep_v,start_col=i,n_cols=ncol_block,&
            transpose_target_m=.TRUE., error=error)
       IF (g_present) THEN
          CALL cp_fm_get_submatrix(matrix_g,rep_g,start_col=i,&
               n_cols=ncol_block,&
               transpose_target_m=.TRUE., error=error)
       ENDIF
       CALL sparse_plus_loc_loct(sparse_matrix,rep_v,rep_g,ncol_block, &
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
! computes sparse=sparse+alpha*v*g^T 
! only v*v^T tested
! *****************************************************************************
  SUBROUTINE cp_sm_plus_fm_fm_t_new(sparse_matrix,matrix_v,matrix_g,last_col,&
       first_col,alpha,error)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_full_matrix_type) , POINTER :: matrix_v
    TYPE(cp_full_matrix_type) , POINTER, OPTIONAL :: matrix_g
    INTEGER, INTENT(IN)               :: last_col
    REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
    INTEGER, OPTIONAL, INTENT(IN)     :: first_col
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE cp_sm_plus_fm_fm_t_new (MODULE qs_blacs)"

!   *** Local variables ***
    REAL(wp), DIMENSION(:,:), POINTER :: loc_v,loc_g,newdat
    REAL(wp) :: prefactor
    INTEGER :: col1,i,j,handle,nrow_global,num_pe,mepos,ictxt_loc,desc(9)
    INTEGER :: nprow,npcol,myprow,mypcol,rb,cb,ncol_max,ishift,in,im,info,handle2,ncol
    INTEGER, DIMENSION(:), POINTER :: ncol_loc
    LOGICAL g_present,FAILURE

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_sm_plus_fm_fm_t_new", &
                                   routineP = moduleN//"/"//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"cp_sm_plus_fm_fm_t would be most happy to use scalapack libraries"//&
     CPSourceFileRef,&
     error)
#endif
    CALL timeset("cp_sm_plus_fm_fm_t_new","I","",handle)
    NULLIFY(newdat,ncol_loc)

    CPPrecondition(ASSOCIATED(sparse_matrix),cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(matrix_v),cp_failure_level,routine,error,failure)
    CPPrecondition(matrix_v%ref_count>0,cp_failure_level,routine,error,failure)

    IF (PRESENT(first_col)) THEN
       col1 = first_col
    ELSE
       col1 = 1
    END IF

    ncol=last_col-col1+1

    g_present=.FALSE.
    IF (PRESENT(matrix_g)) THEN
       CPPrecondition(ASSOCIATED(matrix_g),cp_failure_level,routine,error,failure)
       CPPrecondition(matrix_g%ref_count>0,cp_failure_level,routine,error,failure)
       g_present=.TRUE. 
    END IF
    ! we start with te easiest case

    CALL cp_fm_get_info(matrix_v,nrow_global=nrow_global)
    num_pe=matrix_v%matrix_struct%para_env%num_pe
    mepos =matrix_v%matrix_struct%para_env%mepos
    rb=nrow_global
    cb=1
    im=nrow_global
    ALLOCATE(ncol_loc(0:num_pe-1))

#if defined(__SCALAPACK)
    CALL blacs_get(-1,0,ictxt_loc)
    CALL blacs_gridinit(ictxt_loc,'R',1,num_pe)
    CALL blacs_gridinfo(ictxt_loc,nprow,npcol,myprow,mypcol)
    CALL descinit(desc,nrow_global,ncol,rb,cb,0,0,ictxt_loc,nrow_global,info)
    CPPostcondition(info==0,cp_failure_level,routine,error,failure)
    in=numroc(ncol,cb,mypcol,0,npcol)
    allocate(newdat(im,MAX(1,in)))
    ncol_loc=0
    ncol_loc(mepos)=in
    CALL mp_sum(ncol_loc,matrix_v%matrix_struct%para_env%group)
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(loc_v(ncol_max,nrow_global))

    ! do the actual scalapack to cols reordering
    CALL timeset("cp_sm_plus_fm_fm_t_pdgemr2d","I","",handle2)
    CALL pdgemr2d(nrow_global,ncol,matrix_v%local_data(1,1),1,col1, &
                  matrix_v%matrix_struct%descriptor, &
                  newdat(1,1),1,1,desc,ictxt_loc)
    CALL timestop(0.0_wp,handle2)
    DO j=1,in
       DO i=1,im
          loc_v(j,i)=newdat(i,j)
       END DO
    END DO

    IF (g_present) THEN
       ALLOCATE(loc_g(ncol_max,nrow_global))
       ! do the actual scalapack to cols reordering
       CALL timeset("cp_sm_plus_fm_fm_t_pdgemr2d","I","",handle2)
       CALL pdgemr2d(nrow_global,ncol,matrix_g%local_data(1,1),1,col1, &
                     matrix_g%matrix_struct%descriptor, &
                     newdat(1,1),1,1,desc,ictxt_loc)
       CALL timestop(0.0_wp,handle2)
       DO j=1,in
          DO i=1,im
             loc_g(j,i)=newdat(i,j)
          END DO
       END DO
    ELSE
       loc_g=>loc_v
    ENDIF
    DEALLOCATE(newdat)
    ! EXIT THE GRID 
    CALL blacs_gridexit(ictxt_loc)
#else
    in=ncol
    ncol_loc=0
    ncol_loc(mepos)=in
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(loc_v(ncol_max,nrow_global))
    DO j=1,in
       DO i=1,im
          loc_v(j,i)=matrix_v%local_data(i,col1-1+j)
       END DO
    END DO
    IF (g_present) THEN
       ALLOCATE(loc_g(ncol_max,nrow_global))
       DO j=1,in
          DO i=1,im
             loc_g(j,i)=matrix_g%local_data(i,col1-1+j)
          END DO
       END DO
    ELSE
       loc_g=>loc_v
    ENDIF
#endif

    IF (PRESENT(alpha)) THEN
       prefactor = alpha
    ELSE
       prefactor = 1.0_wp
    END IF

    DO ishift=0,num_pe-1
       CALL sparse_plus_loc_loct(sparse_matrix,loc_v,loc_g,&
                                 ncol_loc(MOD(mepos-ishift+num_pe,num_pe)), &
                                 ncol_max,alpha=prefactor)
       IF (g_present) THEN
          CALL mp_shift(loc_g,matrix_g%matrix_struct%para_env%group)
          CALL mp_shift(loc_v,matrix_v%matrix_struct%para_env%group)
       ELSE
          CALL mp_shift(loc_v,matrix_v%matrix_struct%para_env%group)
       ENDIF
    ENDDO

    IF (g_present) THEN
       DEALLOCATE(loc_v,loc_g)
    ELSE
       DEALLOCATE(loc_v)
    ENDIF

    DEALLOCATE(ncol_loc)
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_plus_fm_fm_t_new
! *****************************************************************************
! computes the LU decomposition of a given matrix
! the actual purpose right now is to compute the determinant of a given matrix
! which is most efficiently done this way, but, indeed, destroys the matrix
! SERIOUS WARNING (KNOWN BUG) : the sign of the determinant depends on ipivot
! one should be able to find out if ipivot is an even or an odd permutation...
! *****************************************************************************
  SUBROUTINE cp_fm_lu_decompose(matrix_a,almost_determinant)
    TYPE(cp_full_matrix_type), POINTER :: matrix_a
    REAL(wp), INTENT(OUT)              :: almost_determinant

! *** locals ***
    INTEGER :: handle,lda,i,n,istat,info,npe
    INTEGER, DIMENSION(9) :: desca
    REAL(wp), DIMENSION(:,:), POINTER :: a
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipivot
    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE cp_fm_lu_decompose"
    REAL(wp), DIMENSION(:), POINTER :: diag
    REAL(wp) :: determinant

    CALL timeset("cp_fm_lu_decompose","I","",handle)

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
    ALLOCATE(ipivot(n))

#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    CALL pdgetrf(n,n,a(1,1),1,1,desca,ipivot,info)
    ALLOCATE(diag(n))
    diag(:)=0.0_wp
    DO i=1,n
       CALL cp_fm_get_element(matrix_a,i,i,diag(i)) !  not completely optimal in speed i would say
    ENDDO
    CALL mp_sum(diag,matrix_a%matrix_struct%para_env%group)
    determinant=1.0_wp
    do i=1,n
       determinant=determinant*diag(i)
    enddo
    DEALLOCATE(diag)
#else
    lda=size(a,1)
    CALL dgetrf(n,n,a(1,1),lda,ipivot,info)
    determinant=1.0_wp
    do i=1,n
       determinant=determinant*a(i,i)
    enddo
#endif
    ! info is allowed to be zero
    ! this does just signal a zero diagonal element
    DEALLOCATE(ipivot)
    almost_determinant=determinant ! notice that the sign is random

  END SUBROUTINE


! *****************************************************************************

  SUBROUTINE cp_fm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,error,b_first_col)

!   Purpose: BLACS interface to the BLAS routine dgemm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    CHARACTER(LEN=1), INTENT(IN)              :: transa,transb
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: k,m,n
    INTEGER, OPTIONAL, INTENT(IN)             :: b_first_col
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,j,lda,ldb,ldc

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c
    REAL(wp) :: flops

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_gemm","I","Gflops",handle)

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

    flops=2.0_wp*(1E-9_wp*m*n*k)/(matrix_a%matrix_struct%context%num_pe(1)*&
         matrix_a%matrix_struct%context%num_pe(2))
    CALL pdgemm(transa,transb,m,n,k,alpha,a(1,1),1,1,desca,b(1,1),1,j,&
                descb,beta,c(1,1),1,j,descc)
#else

    lda = SIZE(a,1)
    ldb = SIZE(b,1)
    ldc = SIZE(c,1)

    flops=2.0_wp*(1E-9_wp*m*n*k)
    CALL dgemm(transa,transb,m,n,k,alpha,a(1,1),lda,b(1,j),ldb,beta,c(1,j),ldc)

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
    REAL(kind=wp) :: flops

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c
    TYPE(cp_blacs_env_type), POINTER    :: context

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

    flops=(2.0E-9_wp*m*m*n)/(nprow*npcol)
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
    CALL timestop(flops,handle)

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

  SUBROUTINE cp_fm_scale(alpha, matrix_a, error)

    TYPE(cp_full_matrix_type), POINTER           :: matrix_a
    REAL(wp), INTENT(IN)                         :: alpha
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    !   *** Local variables ***

    INTEGER                           :: size_a, handle
    REAL(wp), DIMENSION(:,:), POINTER :: a

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_scale',&
                                   routineP=moduleN//':'//routineN
    LOGICAL :: failure

    !   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(a)

    CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_a%ref_count>0,cp_failure_level,routineP,error,failure)

    a => matrix_a%local_data
    size_a = size(a,1)*size(a,2)

    CALL DSCAL(size_a, alpha, a, 1)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_scale

! *****************************************************************************

END MODULE cp_fm_basic_linalg
