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

  USE cp_full_matrix,     ONLY: cp_fm_init_random

  USE qs_blacs,            ONLY: cp_full_matrix_type,cp_fm_get_info, &
                                 cp_full_matrix_p_type, &
                                 cp_fm_set_all, cp_fm_syrk, &
                                 cp_fm_create2, &
                                 cp_fm_release, &
                                 write_blacs_matrix, &
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
                                 blacs_cholesky_restore, &
                                 replicate_blacs_matrix, &
                                 blacs_column_copy, &
                                 blacs_scale_and_d
                                 

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
                                 add_all_real_matrix_blocks, &
                                 add_matrices

  USE message_passing,     ONLY: mp_sum

  IMPLICIT NONE

  PUBLIC  :: preconditioner_type
  PUBLIC  :: init_preconditioner
  PUBLIC  :: destroy_preconditioner
  PUBLIC  :: make_preconditioner
  PUBLIC  :: apply_preconditioner

  INTERFACE make_preconditioner
    MODULE PROCEDURE make_preconditioner_single,make_preconditioner_all
  END INTERFACE

  PRIVATE :: make_sparse_diag, make_full_single, make_local_block

  TYPE preconditioner_type
    PRIVATE
       TYPE(real_matrix_type),    POINTER :: sparse_matrix
       TYPE(cp_full_matrix_type), POINTER :: full_matrix
       TYPE(cp_para_env_type),  POINTER   :: para_env
       TYPE(cp_blacs_env_type), POINTER   :: ctxt
       INTEGER :: in_use
       REAL(wp), DIMENSION(:), POINTER :: all_evals
       REAL(wp), DIMENSION(:), POINTER :: subspace_evals
       REAL(wp) :: energy_gap
  END TYPE preconditioner_type

CONTAINS

SUBROUTINE init_preconditioner(preconditioner_env,globenv)

   TYPE(preconditioner_type)                  :: preconditioner_env
   TYPE(global_environment_type)              :: globenv

   NULLIFY(preconditioner_env%sparse_matrix)
   NULLIFY(preconditioner_env%full_matrix)
   NULLIFY(preconditioner_env%all_evals)
   NULLIFY(preconditioner_env%subspace_evals)
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
   IF (ASSOCIATED(preconditioner_env%all_evals)) THEN
      deallocate(preconditioner_env%all_evals)
   ENDIF
   IF (ASSOCIATED(preconditioner_env%subspace_evals)) THEN
      deallocate(preconditioner_env%subspace_evals)
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
! might overwrite matrix_h, matrix_t
!*******************************************************************************

SUBROUTINE make_preconditioner_single(preconditioner_env, use, matrix_h, matrix_s, &
                                      matrix_t,energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  CHARACTER(LEN=*)                                :: use
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s,matrix_t
  REAL(wp)                                   :: energy_homo,energy_gap

!----
  integer handle

  CALL timeset("make_preconditioner","I","",handle)

  SELECT CASE (use)
  CASE ("SPARSE_DIAG")
    preconditioner_env%in_use=1
    CALL make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                              energy_homo, energy_gap )
  CASE ("FULL_SINGLE")
    preconditioner_env%in_use=2
    CALL make_full_single(preconditioner_env, matrix_h, matrix_s, &
                                          energy_homo, energy_gap )
  CASE ("FULL_S_INVERSE")
    preconditioner_env%in_use=3
    CALL make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                                             energy_homo, energy_gap )
  CASE ("SPARSE_FIRST_ORDER")
    preconditioner_env%in_use=4
    CALL make_sparse_first_order(preconditioner_env, matrix_h, matrix_s, &
                                                 energy_homo, energy_gap )
  CASE ("FULL_S_HALF")
    preconditioner_env%in_use=5
    CALL make_full_s_half(preconditioner_env, matrix_h, matrix_s, &
                                             energy_homo, energy_gap )
  CASE ("FULL_KINETIC")
    preconditioner_env%in_use=6
    CALL make_full_kinetic(preconditioner_env, matrix_t, matrix_s, &
                                             energy_homo, energy_gap )

  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT

  CALL timestop(0.0_wp,handle)

END SUBROUTINE make_preconditioner_single
!********************************************************************************
! is a bit more complicated. for OT the preconditioner it must be such
! that matrix_c0 is the same one a the one from the constraints,
! additionally, we have to rotate it so that diagonalizes it's subspace of h
! and the evecs part of the preconditioner should be such that
! V=[UC0 UX] (otherwise the computation of the lagrangian multipliers gets more
! difficult than one would like
!********************************************************************************
SUBROUTINE make_preconditioner_all(preconditioner_env, use, matrix_h, matrix_s, &
                                   matrix_c0,matrix_sc0,matrix_hc0, energy_gap )

  TYPE(preconditioner_type)                  :: preconditioner_env
  CHARACTER(LEN=*)                           :: use
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_gap
  TYPE(cp_full_matrix_type), POINTER         :: matrix_c0,matrix_sc0,matrix_hc0
!----
  integer handle

  CALL timeset("make_preconditioner","I","",handle)

  SELECT CASE (use)
  CASE ("FULL_ALL")
    preconditioner_env%in_use=103
    CALL make_full_all(preconditioner_env, matrix_h, matrix_s, &
                   matrix_c0,matrix_sc0,matrix_hc0, energy_gap )

  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT

  CALL timestop(0.0_wp,handle)

END SUBROUTINE make_preconditioner_all

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
  CASE (1,4) ! can share the apply since they are both sparse matrices
    CALL apply_sparse_diag(preconditioner_env, matrix_in, matrix_out)
  CASE (2)
    CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
  CASE (3,6)
    CALL apply_full_s_inverse(preconditioner_env, matrix_in, matrix_out)
  CASE (5)
    CALL apply_full_s_half(preconditioner_env, matrix_in, matrix_out)
  CASE (103)
    CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
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
                                                  energy_homo,energy_gap,3)

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
                                 energy_homo,energy_gap, type)

  REAL(wp), DIMENSION(:,:), POINTER :: block_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_s
  REAL(wp), DIMENSION(:,:), POINTER :: block_pre
  REAL(wp)                          :: energy_homo,energy_gap
  INTEGER                           :: type
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
  SELECT CASE (type)
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
  CASE(4)
     block_chol(:,:)=block_s(:,:)
     CALL DPOTRF('U',n,block_chol(1,1),n,info)
     block_pre(:,:)=0.0_wp
     DO i=1,n
          block_pre(i,i)=1.0_wp
     ENDDO
     CALL DTRSM('L','U','N','N',n,n,1.0_wp,block_chol(1,1),n,block_pre(1,1),n)
     CALL DTRSM('R','U','T','N',n,n,1.0_wp,block_chol(1,1),n,block_pre(1,1),n)
  CASE(5) ! like 3 but using s^-0.5 instead of the cholesky decomposition, and not transforming back
     block_evec(:,:)=block_s(:,:)
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     if (info.ne.0) CALL stop_program("preconditioner","DSYEVD S")
     block_pre(:,:)=0.0_wp
     do i=1,n
	block_pre(i,i)=1.0_wp/sqrt(evals(i))
     enddo
     ! block_pre is s^-0.5
     CALL DGEMM('N','N',n,n,n,1.0_wp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_wp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_wp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_wp,block_pre(1,1),n)
     ! transform H
     block_evec(:,:)=block_h(:,:)
     CALL DGEMM('N','N',n,n,n,1.0_wp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_wp,block_buf1(1,1),n)
     CALL DGEMM('N','N',n,n,n,1.0_wp,block_pre(1,1),n,block_buf1(1,1),n, &
                                                  0.0_wp,block_evec(1,1),n)
     ! get evals and evecs
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     if (info.ne.0) CALL stop_program("preconditioner","DSYEVD H")
     block_pre(:,:)=0.0_wp
     DO i=1,n
          block_pre(i,i)=1.0_wp/MAX(evals(i)-energy_homo,energy_gap)
     ENDDO
     CALL DGEMM('N','N',n,n,n,1.0_wp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_wp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_wp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_wp,block_pre(1,1),n)
  CASE(6) ! like 3 not doing any transformation with s before or after (supposedly done by the caller)
     block_evec(:,:)=block_h(:,:)
     ! get evals and evecs
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     if (info.ne.0) CALL stop_program("preconditioner","DSYEVD H")
     block_pre(:,:)=0.0_wp
     DO i=1,n
          block_pre(i,i)=1.0_wp/MAX(evals(i),energy_gap)
     ENDDO
     CALL DGEMM('N','N',n,n,n,1.0_wp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_wp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_wp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_wp,block_pre(1,1),n)
     
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
! a sparse one based on 2inv(A)-inv(A) F inv(A) where F is S
! turns out that this inverse is not symmetric positive definite unless
! inv(A) is scaled by some number
SUBROUTINE make_sparse_first_order(preconditioner_env, matrix_h, matrix_s, &
                                   energy_homo, energy_gap )

  TYPE(preconditioner_type)                  :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap

!----
  INTEGER :: iblock_row,iblock_col_s,iblock_col_h,nblocks,n,handle,k

  TYPE(real_block_node_type), POINTER :: block_node_s
  TYPE(real_block_node_type), POINTER :: block_node_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_s
  REAL(wp), DIMENSION(:,:), POINTER :: block_pre
  REAL(wp), DIMENSION(:,:), POINTER :: block_tmp
  INTEGER, DIMENSION(:), POINTER :: first_row,last_row
  INTEGER, DIMENSION(:),POINTER :: offset
  REAL(wp), DIMENSION(:), POINTER :: block_buffer


  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
     CALL deallocate_matrix(preconditioner_env%sparse_matrix)
     NULLIFY(preconditioner_env%sparse_matrix)
  ENDIF
  ! first we get the all blocks of the block diagonal preconditioner replicated
  ! on every CPU (store in block_buffer)


  CALL get_matrix_info(matrix_s,nblock_row=nblocks,first_row=first_row,last_row=last_row)

  ALLOCATE(offset(nblocks+1))
  offset(1)=0
  DO iblock_row=1,nblocks
     offset(iblock_row+1)=offset(iblock_row)+(last_row(iblock_row)-first_row(iblock_row)+1)**2
  ENDDO
  ALLOCATE(block_buffer(offset(nblocks+1)))
  CALL DCOPY(offset(nblocks+1),0.0_wp,0,block_buffer(1),1)

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
                                                  energy_homo,energy_gap,4)
           CALL DCOPY(n*n,block_pre(1,1),1,block_buffer(offset(iblock_row)+1),1)

           DEALLOCATE(block_pre)

        ENDIF
        block_node_s => next_block_node(block_node_s)
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO

  CALL mp_sum(block_buffer,preconditioner_env%para_env%group)

  CALL DSCAL(offset(nblocks+1),0.2_wp,block_buffer(1),1)

  ! now every diagonal block is everywhere, form the product
  CALL Replicate_matrix(matrix_s,preconditioner_env%sparse_matrix,"PRECONDITIONER",.false.)

  DO iblock_row=1,nblocks
     block_node_s => first_block_node(matrix_s,iblock_row)
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_s) .or. ASSOCIATED(block_node_h))

        CALL get_block_node(block_node=block_node_s,&
                            block_col=iblock_col_s,&
                            block=block_s)

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            block=block_h)

        n=size(block_s,1)
        k=size(block_s,2)
        ALLOCATE(block_pre(n,k))
        ALLOCATE(block_tmp(n,k))
	block_pre(:,:)=-block_s(:,:)
        CALL DGEMM('N','N',n,k,k,1.0_wp,block_pre(1,1),n, &
                                        block_buffer(offset(iblock_col_s)+1), &
                                        k,0.0_wp,block_tmp(1,1),n)
        CALL DGEMM('N','N',n,k,n,1.0_wp,block_buffer(offset(iblock_row)+1),n, &
                                        block_tmp(1,1), &
                                        n,0.0_wp,block_pre(1,1),n)
        if (iblock_row .eq. iblock_col_s) &
           CALL DAXPY(n*k,2.0_wp,block_buffer(offset(iblock_row)+1),1,block_pre(1,1),1)

        CALL add_block_node(preconditioner_env%sparse_matrix,iblock_row,iblock_col_s,block_pre)

        DEALLOCATE(block_pre)
        DEALLOCATE(block_tmp)
        block_node_s => next_block_node(block_node_s)
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO

  DEALLOCATE(offset)
  DEALLOCATE(block_buffer)

END SUBROUTINE make_sparse_first_order

!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! a full matrix preconditioner
SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  integer n,k
  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  CALL cp_fm_gemm('N','N',n,k,n,1.0_wp,preconditioner_env%full_matrix, &
                  matrix_in,0.0_wp,matrix_out)
END SUBROUTINE apply_full_single

SUBROUTINE make_full_single(preconditioner_env, matrix_h, matrix_s, &
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

END SUBROUTINE make_full_single

!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! only inv(S)
! a full matrix preconditioner
SUBROUTINE apply_full_s_inverse(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  integer n,k
  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  CALL cp_fm_symm('L','U',n,k,1.0_wp,preconditioner_env%full_matrix, &
                        matrix_in, &
                        0.0_wp,matrix_out)

END SUBROUTINE apply_full_s_inverse
SUBROUTINE make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                               energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap

! ---
  INTEGER, PARAMETER :: nrow_block = 32 ! where do I get them from
  INTEGER, PARAMETER :: ncol_block = 32 ! where do I get them from
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


  CALL copy_sparse_to_blacs_matrix(matrix_s,preconditioner_env%full_matrix)
  CALL blacs_cholesky_decompose(preconditioner_env%full_matrix)
  CALL blacs_cholesky_invert(preconditioner_env%full_matrix)

END SUBROUTINE make_full_s_inverse
SUBROUTINE make_full_kinetic(preconditioner_env, matrix_t, matrix_s, &
                               energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_t,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap,shift

! ---
  INTEGER, PARAMETER :: nrow_block = 32 ! where do I get them from
  INTEGER, PARAMETER :: ncol_block = 32 ! where do I get them from
  INTEGER :: i,n
! ---

  IF (ASSOCIATED(preconditioner_env%full_matrix)) THEN
     CALL cp_fm_release(preconditioner_env%full_matrix)
  ENDIF
  CALL get_matrix_info(matrix_t,nrow=n)

  CALL cp_fm_create2(new_matrix=preconditioner_env%full_matrix,&
                             nrow_global=n,&
                             ncol_global=n,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  shift=MAX(-energy_homo,energy_gap)
  CALL add_matrices(matrix_t,1.0_wp,matrix_t,shift,matrix_s)

  CALL copy_sparse_to_blacs_matrix(matrix_t,preconditioner_env%full_matrix)
  CALL blacs_cholesky_decompose(preconditioner_env%full_matrix)
  CALL blacs_cholesky_invert(preconditioner_env%full_matrix)

END SUBROUTINE make_full_kinetic


SUBROUTINE apply_full_s_half(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  integer n,k
  TYPE(cp_full_matrix_type), POINTER :: new_matrix
  NULLIFY(new_matrix)

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  CALL replicate_blacs_matrix(matrix_in,new_matrix,"apply prec dummy")

  CALL cp_fm_symm('L','U',n,k,1.0_wp,preconditioner_env%full_matrix, &
                        matrix_in, &
                        0.0_wp,new_matrix)

  CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,new_matrix, &
                         matrix_out,k, para_env=preconditioner_env%para_env)

  ! CALL copy_blacs_to_blacs_matrix(new_matrix,matrix_out)

  CALL cp_fm_symm('L','U',n,k,1.0_wp,preconditioner_env%full_matrix, &
                        matrix_out, &
                        0.0_wp,new_matrix)

  CALL copy_blacs_to_blacs_matrix(new_matrix,matrix_out)

  CALL cp_fm_release(new_matrix)

END SUBROUTINE apply_full_s_half

SUBROUTINE make_full_s_half(preconditioner_env, matrix_h, matrix_s, &
                               energy_homo, energy_gap )

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_homo,energy_gap

! ---
  INTEGER, PARAMETER :: nrow_block = 32 ! where do I get them from
  INTEGER, PARAMETER :: ncol_block = 32 ! where do I get them from
  INTEGER :: i,n
  TYPE(cp_full_matrix_type), POINTER :: fm_y,fm_yi,fm_z,fm_zi

  INTEGER :: iblock_row,iblock_col_s,iblock_col_h,nblocks

  TYPE(real_block_node_type), POINTER :: block_node_s
  TYPE(real_block_node_type), POINTER :: block_node_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_h
  REAL(wp), DIMENSION(:,:), POINTER :: block_s
  REAL(wp), DIMENSION(:,:), POINTER :: block_pre
! ---

  IF (ASSOCIATED(preconditioner_env%full_matrix)) THEN
     CALL cp_fm_release(preconditioner_env%full_matrix)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)


  CALL cp_fm_create2(new_matrix=fm_y, &
                             nrow_global=n, ncol_global=n,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_y",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_z, &
                             nrow_global=n, ncol_global=n,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_z",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_yi, &
                             nrow_global=n, ncol_global=n,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_yn",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_zi, &
                             nrow_global=n, ncol_global=n,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_zn",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL copy_sparse_to_blacs_matrix(matrix_s,fm_y)
  CALL copy_blacs_to_blacs_matrix(fm_y,fm_z)
  CALL blacs_cholesky_decompose(fm_z)
  CALL blacs_cholesky_invert(fm_z)
  CALL blacs_scale_and_d(fm_z,0.5_wp,0.5_wp)
  CALL blacs_scale_and_d(fm_y,0.5_wp,0.5_wp)
  

  do i=1,4
     CALL copy_blacs_to_blacs_matrix(fm_z,fm_zi)
     CALL blacs_cholesky_decompose(fm_zi)
     CALL blacs_cholesky_invert(fm_zi)

     CALL copy_blacs_to_blacs_matrix(fm_y,fm_yi)
     CALL blacs_cholesky_decompose(fm_yi)
     CALL blacs_cholesky_invert(fm_yi)

     CALL cp_fm_add(0.5_wp,fm_y,0.5_wp,fm_zi)
     CALL cp_fm_add(0.5_wp,fm_z,0.5_wp,fm_yi)
  enddo


  CALL cp_fm_create2(new_matrix=preconditioner_env%full_matrix,&
                             nrow_global=n,&
                             ncol_global=n,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL copy_blacs_to_blacs_matrix(fm_z,preconditioner_env%full_matrix)

  CALL copy_sparse_to_blacs_matrix(matrix_h,fm_y)
  CALL copy_sparse_to_blacs_matrix(matrix_s,fm_zi)
  CALL cp_fm_add(1.0_wp,fm_y,-energy_homo,fm_zi) ! fm_y = H-eS
  CALL blacs_scale_and_d(fm_zi,0.0_wp,1.0_wp) ! ugly way to get the full matrix
  CALL cp_fm_symm('L','U',n,n,1.0_wp,fm_y, &
                        fm_zi, &
                        0.0_wp,fm_yi)
  CALL cp_fm_symm('L','U',n,n,1.0_wp,fm_z, &
                        fm_yi, &
                        0.0_wp,fm_zi)
  CALL cp_fm_symm('R','U',n,n,1.0_wp,fm_z, &
                        fm_zi, &
                        0.0_wp,fm_yi)
  CALL copy_blacs_to_sparse_matrix(fm_yi,matrix_h)

  CALL cp_fm_release(fm_zi)
  CALL cp_fm_release(fm_yi)
  CALL cp_fm_release(fm_z)
  CALL cp_fm_release(fm_y)

  ! and add a diagonal approximation for inv(Shalf)Hinv(Shalf)

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
                                                  energy_homo,energy_gap,6)

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

END SUBROUTINE make_full_s_half


!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! a full matrix preconditioner
SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out)

  TYPE(preconditioner_type)               :: preconditioner_env
  TYPE(cp_full_matrix_type), POINTER            :: matrix_in,matrix_out
!----
  TYPE(cp_full_matrix_type), POINTER            :: new_matrix
  integer n,k,ln,lk,i,j
  REAL(wp) :: dum
  REAL(wp), DIMENSION(:,:), POINTER :: local_data
  INTEGER, DIMENSION(:), POINTER :: row_indices,col_indices
  ! ok need an extra buff here, we don't want to destroy matrix_in
  NULLIFY(new_matrix)
  CALL replicate_blacs_matrix(matrix_in,new_matrix,"apply prec dummy")

  CALL cp_fm_get_info(new_matrix,nrow_global=n,ncol_global=k, &
                      nrow_local=ln,ncol_local=lk,local_data=local_data, &
                      row_indices=row_indices,col_indices=col_indices)

  ! first mult by (inv(U)K)^T
  CALL cp_fm_gemm('T','N',n,k,n,1.0_wp,preconditioner_env%full_matrix, &
                  matrix_in,0.0_wp,new_matrix)
  ! do the right scaling
  do j=1,lk
  do i=1,ln
     dum=1.0_wp/MAX(preconditioner_env%energy_gap, &
             preconditioner_env%all_evals(row_indices(i))- &
             preconditioner_env%subspace_evals(col_indices(j))  )
     local_data(i,j)=local_data(i,j)*dum
  enddo
  enddo
  ! mult back
  CALL cp_fm_gemm('N','N',n,k,n,1.0_wp,preconditioner_env%full_matrix, &
                  new_matrix,0.0_wp,matrix_out)
  CALL cp_fm_release(new_matrix)

END SUBROUTINE apply_full_all

SUBROUTINE make_full_all(preconditioner_env, matrix_h, matrix_s, &
                         matrix_c0,matrix_sc0,matrix_hc0, energy_gap )

  TYPE(preconditioner_type)                  :: preconditioner_env
  TYPE(real_matrix_type), POINTER            :: matrix_h,matrix_s
  REAL(wp)                                   :: energy_gap
  TYPE(cp_full_matrix_type), POINTER         :: matrix_c0,matrix_sc0,matrix_hc0

! ---
  TYPE(cp_full_matrix_type), POINTER :: fm_h,fm_s, fm_wvec, fm_evec, fm_xvec, fm_wxvec, fm_block
  INTEGER, PARAMETER :: nrow_block = 32 ! where do I get them from
  INTEGER, PARAMETER :: ncol_block = 32 ! where do I get them from
  REAL(wp), DIMENSION(:), POINTER :: evals
  INTEGER :: i,n,k,x
! ---

  IF (ASSOCIATED(preconditioner_env%full_matrix)) THEN 
     CALL cp_fm_release(preconditioner_env%full_matrix)
  ENDIF
  IF (ASSOCIATED(preconditioner_env%all_evals)) THEN 
     deallocate(preconditioner_env%all_evals)
  ENDIF
  IF (ASSOCIATED(preconditioner_env%subspace_evals)) THEN 
     deallocate(preconditioner_env%subspace_evals)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)
  CALL cp_fm_get_info(matrix_c0,ncol_global=k)

  allocate(preconditioner_env%all_evals(n))
  allocate(preconditioner_env%subspace_evals(k))
  preconditioner_env%energy_gap=energy_gap
  CALL cp_fm_create2(new_matrix=preconditioner_env%full_matrix,&
                             nrow_global=n,&
                             ncol_global=n,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_h",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  ! first step, make c0,sc0,hc0 ritz vectors
  NULLIFY(fm_wvec,fm_evec,fm_block)
  CALL cp_fm_create2(new_matrix=fm_wvec, nrow_global=n, ncol_global=k,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="make full all prec matrix wvec",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_evec, nrow_global=k, ncol_global=k,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="make full all prec matrix evec",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_block, nrow_global=k, ncol_global=k,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="make full all prec matrix evec",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_gemm('T','N',k,k,n,1.0_wp, &
                    matrix_c0,matrix_hc0,0.0_wp,fm_block)
 
  CALL blacs_syevd(fm_block,fm_evec,preconditioner_env%subspace_evals)

  CALL cp_fm_gemm('N','N',n,k,k,1.0_wp, &
                      matrix_c0,fm_evec,0.0_wp,fm_wvec)
  CALL copy_blacs_to_blacs_matrix(fm_wvec,matrix_c0)

  CALL cp_fm_gemm('N','N',n,k,k,1.0_wp, &
                      matrix_sc0,fm_evec,0.0_wp,fm_wvec)
  CALL copy_blacs_to_blacs_matrix(fm_wvec,matrix_sc0)

  CALL cp_fm_gemm('N','N',n,k,k,1.0_wp, &
                      matrix_hc0,fm_evec,0.0_wp,fm_wvec)
  CALL copy_blacs_to_blacs_matrix(fm_wvec,matrix_hc0)

  CALL cp_fm_release(fm_block) ! fm_wvec still needed
  CALL cp_fm_release(fm_evec)

  ! now lets get things in an orthogonal basis and seperate the c0 space from the other space
  NULLIFY(fm_h,fm_s)
  CALL cp_fm_create2(new_matrix=fm_h, nrow_global=n, ncol_global=n,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="matrix_h",&                             
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_s, nrow_global=n, ncol_global=n,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="matrix_s",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL copy_sparse_to_blacs_matrix(matrix_h,fm_h)
  CALL copy_sparse_to_blacs_matrix(matrix_s,fm_s)
  CALL blacs_cholesky_decompose(fm_s)
  CALL blacs_cholesky_reduce(fm_h,fm_s)
  ! let's find a space of orbitals orthogonal to Uc0
  CALL blacs_cholesky_restore(matrix_c0,k,fm_s,fm_wvec,"MULTIPLY") ! orthogonal c0
  x=n-k ! size of the excited subspace
  NULLIFY(fm_xvec,fm_wxvec)

  CALL cp_fm_create2(new_matrix=fm_xvec, nrow_global=n, ncol_global=x,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="matrix_fm_xvec",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_wxvec, nrow_global=n, ncol_global=x,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="matrix_fm_wxvec",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)

  CALL cp_fm_create2(new_matrix=fm_block, nrow_global=x, ncol_global=k,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_block",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  ! orthogonal to uc0
  CALL cp_fm_init_random(fm_xvec,x)
  CALL cp_fm_gemm('T','N',x,k,n,1.0_wp,fm_xvec,fm_wvec,0.0_wp,fm_block)
  CALL cp_fm_gemm('N','T',n,x,k,-1.0_wp,fm_wvec,fm_block,1.0_wp,fm_xvec)
  CALL cp_fm_release(fm_block)
  ! orthonormalise them
  CALL cp_fm_create2(new_matrix=fm_block, nrow_global=x, ncol_global=x,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_block",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL cp_fm_gemm('T','N',x,x,n,1.0_wp,fm_xvec,fm_xvec,0.0_wp,fm_block)
  CALL blacs_cholesky_decompose(fm_block)
  CALL blacs_cholesky_restore(fm_xvec,x,fm_block,fm_wxvec,"SOLVE","RIGHT") 
  CALL copy_blacs_to_blacs_matrix(fm_wxvec,fm_xvec)
  ! get their H matrix subspace
  CALL cp_fm_gemm('N','N',n,x,n,1.0_wp,fm_h,fm_xvec,0.0_wp,fm_wxvec)
  CALL cp_fm_gemm('T','N',x,x,n,1.0_wp,fm_wxvec,fm_xvec,0.0_wp,fm_block)
  ALLOCATE(evals(x))
  CALL cp_fm_create2(new_matrix=fm_evec, nrow_global=x, ncol_global=x,&
                             nrow_block=nrow_block, ncol_block=ncol_block,&
                             name="fm_block",&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL blacs_syevd(fm_block,fm_evec,evals)
  ! and make them ritz vectors
  CALL cp_fm_gemm('N','N',n,x,x,1.0_wp, &
                      fm_xvec,fm_evec,0.0_wp,fm_wxvec)

  ! now assemble c0 and xvec
  CALL blacs_column_copy(fm_wvec,fm_h,k,1,1)
  CALL blacs_column_copy(fm_wxvec,fm_h,x,1,k+1)

  CALL cp_fm_gemm('T','N',n,n,n,1.0_wp,fm_h,fm_h,0.0_wp,preconditioner_env%full_matrix)
  ! and get inv(U)V
  CALL blacs_cholesky_restore(fm_h,n,fm_s,preconditioner_env%full_matrix,"SOLVE")
  preconditioner_env%all_evals(1:k)=preconditioner_env%subspace_evals(1:k)
  preconditioner_env%all_evals(k+1:n)=evals(1:x)
  DEALLOCATE(evals)

  CALL cp_fm_release(fm_block)
  CALL cp_fm_release(fm_xvec)
  CALL cp_fm_release(fm_wxvec)
  CALL cp_fm_release(fm_evec)
  CALL cp_fm_release(fm_wvec)
  CALL cp_fm_release(fm_h)
  CALL cp_fm_release(fm_s)

END SUBROUTINE make_full_all


END MODULE preconditioner

