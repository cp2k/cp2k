!******************************************************************************
!!
!!   NAME
!!     preconditioner
!!
!!   FUNCTION
!!     computes preconditioners, and implements methods to apply them
!!     currently used in qs_ot
!!
!!   AUTHOR
!!     Joost VandeVondele (09.2002)
!!
!!   MODIFICATION HISTORY
!!     None
!!
!!   SOURCE
!******************************************************************************

MODULE preconditioner

! *****************************************************************************
  USE kinds, ONLY: wp => dp, wp_size => dp_size
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE cp_para_types, ONLY: cp_para_env_type, cp_blacs_env_type
  USE cp_para_env,       ONLY: cp_para_env_retain, cp_para_env_release,&
       cp_para_env_from_globenv
  USE cp_blacs_env,      ONLY: cp_blacs_env_retain, cp_blacs_env_release,&
       cp_blacs_env_from_globenv

  USE qs_blacs,            ONLY: cp_full_matrix_type,cp_fm_get_info, &
                                 cp_full_matrix_p_type, &
                                 cp_fm_set_all, cp_fm_syrk, &
                                 cp_fm_create2, &
                                 cp_fm_release, &
                                 cp_fm_symm,cp_fm_gemm, &
                                 blacs_make_basis,copy_blacs_to_blacs_matrix, &
                                 cp_fm_trace, blacs_transpose, cp_fm_add, &
                                 cp_fm_schur_product, &
                                 blacs_cholesky_decompose, &
                                 blacs_cholesky_invert, &
                                 blacs_cholesky_reduce, &
                                 blacs_syevx, blacs_syevd, &
                                 blacs_scale_matrix, &
                                 blacs_set_local_block, &
                                 cp_sm_fm_multiply, &
                                 copy_sparse_to_blacs_matrix, &
                                 copy_blacs_to_sparse_matrix, &
                                 blacs_scale_matrix, &
                                 blacs_cholesky_restore

  USE global_types,        ONLY: global_environment_type
  USE sparse_matrix_types, ONLY: real_matrix_type, &
                                 replicate_matrix, &
                                 get_matrix_info, &
                                 get_block_node, &
                                 add_block_node , &
                                 real_block_node_type, &
                                 first_block_node, &
                                 next_block_node, &
                                 deallocate_matrix, &
                                 add_all_real_matrix_blocks

  IMPLICIT NONE

  PUBLIC  :: preconditioner_type
  PUBLIC  :: init_preconditioner
  PUBLIC  :: destroy_preconditioner
  PUBLIC  :: make_preconditioner
  PUBLIC  :: apply_preconditioner

  PRIVATE :: make_sparse_diag, make_full, make_local_block

  TYPE preconditioner_type
       TYPE(real_matrix_type),    POINTER :: sparse_matrix
       TYPE(cp_full_matrix_type), POINTER :: full_matrix
       TYPE(cp_para_env_type),  POINTER   :: para_env
       TYPE(cp_blacs_env_type), POINTER   :: ctxt
       INTEGER :: in_use
  END TYPE preconditioner_type

CONTAINS

SUBROUTINE init_preconditioner(preconditioner_env,globenv)

   TYPE(preconditioner_type)                  :: preconditioner_env
   TYPE(global_environment_type)              :: globenv

   NULLIFY(preconditioner_env%sparse_matrix)
   NULLIFY(preconditioner_env%full_matrix)
   preconditioner_env%in_use=0
   preconditioner_env%para_env => cp_para_env_from_globenv(globenv)
   preconditioner_env%ctxt     => cp_blacs_env_from_globenv(globenv)
   CALL cp_para_env_retain(preconditioner_env%para_env)
   CALL cp_blacs_env_retain(preconditioner_env%ctxt)

END SUBROUTINE init_preconditioner

SUBROUTINE destroy_preconditioner(preconditioner_env)

   TYPE(preconditioner_type) :: preconditioner_env

   IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
      CALL deallocate_matrix(preconditioner_env%sparse_matrix)
      NULLIFY(preconditioner_env%sparse_matrix)
   ENDIF

   IF (ASSOCIATED(preconditioner_env%full_matrix)) THEN
      CALL cp_fm_release(preconditioner_env%full_matrix)
   ENDIF

   CALL cp_para_env_release(preconditioner_env%para_env)
   CALL cp_blacs_env_release(preconditioner_env%ctxt)
 
   preconditioner_env%in_use=0

END SUBROUTINE destroy_preconditioner

!*******************************************************************************
! creates a preconditioner for the system (H-energy_homo S)
! this preconditioner is (must be) symmetric positive definite.
! currently uses a atom-block-diagonal form
! each block will be  ....
!*******************************************************************************

SUBROUTINE make_preconditioner(preconditioner_env, use, matrix_h, matrix_s, &
                                  energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  CHARACTER(LEN=*)                                :: use
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap

!----
  integer handle

  CALL timeset("make_preconditioner","I","",handle)

  SELECT CASE (use)
  CASE ("SPARSE_DIAGONAL")
    preconditioner_env%in_use=1
    CALL make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                              energy_homo, energy_gap )
  CASE ("FULL")
    preconditioner_env%in_use=2
    CALL make_full(preconditioner_env, matrix_h, matrix_s, &
                              energy_homo, energy_gap )

  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT

  CALL timestop(0.0_wp,handle)

END SUBROUTINE make_preconditioner

!*******************************************************************************
! applies a previously created preconditioner to a full matrix
!*******************************************************************************
SUBROUTINE apply_preconditioner(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  integer handle

  CALL timeset("apply_preconditioner","I","",handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
    CALL stop_program("apply preconditioner","No preconditioner in use")
  CASE (1)
    CALL apply_sparse_diag(preconditioner_env, matrix_in, matrix_out)
  CASE (2)
    CALL apply_full(preconditioner_env, matrix_in, matrix_out)
  CASE DEFAULT
    CALL stop_program("apply preconditioner","implemented")
  END SELECT

  CALL timestop(0.0_wp,handle)

END SUBROUTINE apply_preconditioner

!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! a sparse block diagonal approximation
SUBROUTINE apply_sparse_diag(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  integer k
  CALL cp_fm_get_info(matrix_in,ncol_global=k)
  CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_in, &
                         matrix_out,k, para_env=preconditioner_env%para_env)

END SUBROUTINE apply_sparse_diag

SUBROUTINE make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                          energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap

!----
  INTEGER :: iblock_row,iblock_col_s,iblock_col_h,nblocks,n,handle

  TYPE(real_block_node_type), POINTER :: block_node_s
  TYPE(real_block_node_type), POINTER :: block_node_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_s
  REAL(wp), DIMENSION(:,:), POINTER :: block_pre


  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN 
     CALL deallocate_matrix(preconditioner_env%sparse_matrix)
     NULLIFY(preconditioner_env%sparse_matrix)
  ENDIF

  CALL Replicate_matrix(matrix_s,preconditioner_env%sparse_matrix,"PRECONDITIONER",.False.)
  CALL get_matrix_info(matrix_s,nblock_row=nblocks)
  DO iblock_row=1,nblocks
     block_node_s => first_block_node(matrix_s,iblock_row)
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_s) .or. ASSOCIATED(block_node_h))
        IF (.NOT.ASSOCIATED(block_node_s) .or. .NOT.ASSOCIATED(block_node_h)) &
            CALL stop_program("qs preconditioner", ".NOT.ASSOCIATED")
         
        CALL get_block_node(block_node=block_node_s,&
                            block_col=iblock_col_s,&
                            block=block_s)

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            block=block_h)

        IF (iblock_col_h .ne. iblock_col_s) &
           CALL stop_program("qs preconditioner","iblock_col_h .ne. iblock_col_s")

        IF (iblock_col_s .eq. iblock_row) THEN
           n=size(block_s,1)
           ALLOCATE(block_pre(n,n))

           CALL make_local_block(block_h,block_s,block_pre,  &
                                                  energy_homo,energy_gap)

           CALL add_block_node(matrix=preconditioner_env%sparse_matrix,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col_s,&
                                 block=block_pre)

           DEALLOCATE(block_pre)
 
        ENDIF
        block_node_s => next_block_node(block_node_s)
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO
END SUBROUTINE make_sparse_diag
SUBROUTINE make_local_block(block_h,block_s,block_pre, &
                                                energy_homo,energy_gap)

  REAL(wp), DIMENSION(:,:), POINTER :: block_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_s
  REAL(wp), DIMENSION(:,:), POINTER :: block_pre
  REAL(wp)                          :: energy_homo,energy_gap

  REAL(wp), DIMENSION(:,:), POINTER :: block_chol
  REAL(wp), DIMENSION(:,:), POINTER :: block_evec
  REAL(wp), DIMENSION(:,:), POINTER :: block_buf1
  REAL(wp), DIMENSION(:,:), POINTER :: block_buf2
  REAL(wp), DIMENSION(:), POINTER :: evals 
  INTEGER,  DIMENSION(:), POINTER :: iwork
  REAL(wp), DIMENSION(:), POINTER :: work
  INTEGER :: n,i,info,liwork,lwork,istat,j

  n=size(block_s,1)
  lwork=1+6*n+2*n**2+50
  liwork=5*n+3
  ALLOCATE(block_chol(n,n))
  ALLOCATE(block_evec(n,n))
  ALLOCATE(block_buf1(n,n))
  ALLOCATE(evals(n))
  ALLOCATE(work(lwork),STAT=istat)
  IF (istat.NE.0) CALL stop_memory("use preconditioner","work")
  ALLOCATE(iwork(liwork),STAT=istat)
  IF (istat.NE.0) CALL stop_memory("use preconditioner","iwork")
 
  block_pre(:,:)=0.0_wp
  SELECT CASE (3)
  CASE(1)
    DO i=1,n
      block_pre(i,i)=1.0_wp
    ENDDO
  CASE(2)
    DO i=1,n
      block_pre(i,i)=1.0_wp/MAX(energy_gap,block_h(i,i)-energy_homo)
    ENDDO
  CASE(3)
     ! more difficult constuct something like S^-0.5 K^T CASE(2) K^T S^-0.5
     block_chol(:,:)=block_s(:,:)
     block_evec(:,:)=block_h(:,:)
     CALL DPOTRF('U',n,block_chol(1,1),n,info)
     IF (info.ne.0)  CALL stop_program("use preconditioner","Error dpotrf")
     CALL DSYGST(1,'U',n,block_evec(1,1),n,block_chol(1,1),n,info)
     IF (info.ne.0)  CALL stop_program("use preconditioner","Error dsygst")
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     IF (info.NE.0) CALL stop_program("use preconditioner","problems")
     block_pre(:,:)=0.0_wp
     DO i=1,n
        block_pre(i,i)=1.0_wp/MAX(evals(i)-energy_homo,energy_gap)
     ENDDO
     ! K = V E V ^ T
     CALL DGEMM('N','N',n,n,n,1.0_wp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_wp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_wp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_wp,block_pre(1,1),n)
     ! inv(U) K inv(U)^T
     CALL DTRSM('L','U','N','N',n,n,1.0_wp,block_chol(1,1),n,block_pre(1,1),n)
     CALL DTRSM('R','U','T','N',n,n,1.0_wp,block_chol(1,1),n,block_pre(1,1),n)
 
  END SELECT

  DEALLOCATE(iwork)
  DEALLOCATE(work)
  DEALLOCATE(block_chol)
  DEALLOCATE(block_evec)
  DEALLOCATE(block_buf1)
  DEALLOCATE(evals)

END SUBROUTINE make_local_block

!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! a full matrix preconditioner
SUBROUTINE apply_full(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  integer n,k
  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  CALL cp_fm_gemm('N','N',n,k,n,1.0_wp,preconditioner_env%full_matrix, &
                  matrix_in,0.0_wp,matrix_out)
END SUBROUTINE apply_full

SUBROUTINE make_full(preconditioner_env, matrix_h, matrix_s, &
                       energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap

! ---
  TYPE(cp_full_matrix_type), POINTER :: fm_h,fm_s
  INTEGER, PARAMETER :: nrow_block = 32 ! where do I get them from
  INTEGER, PARAMETER :: ncol_block = 32 ! where do I get them from
  REAL(wp), DIMENSION(:), POINTER :: evals
  INTEGER :: i,n
! ---

  IF (ASSOCIATED(preconditioner_env%full_matrix)) THEN 
     CALL cp_fm_release(preconditioner_env%full_matrix)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)

  CALL cp_fm_create2(new_matrix=preconditioner_env%full_matrix,&
                             nrow_global=n,&
                             ncol_global=n,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)


  ALLOCATE(evals(n))
  NULLIFY(fm_h,fm_s)
  CALL cp_fm_create2(new_matrix=fm_h,&
                             nrow_global=n,&
                             ncol_global=n,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h",&                             
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_s,&
                             nrow_global=n,&
                             ncol_global=n,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_s",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL copy_sparse_to_blacs_matrix(matrix_h,fm_h)
  CALL copy_sparse_to_blacs_matrix(matrix_s,fm_s)
  CALL blacs_cholesky_decompose(fm_s)
  CALL blacs_cholesky_reduce(fm_h,fm_s)
  CALL blacs_syevd(fm_h,preconditioner_env%full_matrix,evals)
  CALL blacs_cholesky_restore(preconditioner_env%full_matrix,n,fm_s,fm_h,"SOLVE") 
  do i=1,n
	evals(i)=1.0_wp/MAX(evals(i)-energy_homo,energy_gap)
  enddo
  CALL copy_blacs_to_blacs_matrix(fm_h,preconditioner_env%full_matrix)
  CALL blacs_scale_matrix(preconditioner_env%full_matrix,evals,n)
  CALL cp_fm_gemm('N','T',n,n,n,1.0_wp,preconditioner_env%full_matrix,fm_h,0.0_wp,fm_s)
  CALL copy_blacs_to_blacs_matrix(fm_s,preconditioner_env%full_matrix)

  DEALLOCATE(evals)
  CALL cp_fm_release(fm_h)
  CALL cp_fm_release(fm_s)

END SUBROUTINE make_full

END MODULE preconditioner

