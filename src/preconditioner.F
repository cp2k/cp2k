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
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             add_matrices,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             replicate_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PUBLIC  :: preconditioner_type, preconditioner_p_type
  PUBLIC  :: init_preconditioner
  PUBLIC  :: preconditioner_in_use
  PUBLIC  :: destroy_preconditioner
  PUBLIC  :: make_preconditioner
  PUBLIC  :: apply_preconditioner

  INTERFACE make_preconditioner
    MODULE PROCEDURE make_preconditioner_single 
  END INTERFACE

  PRIVATE :: make_sparse_diag, make_full_single, make_local_block

  TYPE preconditioner_type
    PRIVATE
       TYPE(real_matrix_type),    POINTER :: sparse_matrix
       TYPE(cp_fm_type), POINTER :: fm
       TYPE(cp_para_env_type),  POINTER   :: para_env
       TYPE(cp_blacs_env_type), POINTER   :: ctxt
       INTEGER :: in_use
       REAL(KIND = dp), DIMENSION(:), POINTER :: subspace_evals
       REAL(KIND = dp) :: energy_gap
  END TYPE preconditioner_type

  TYPE preconditioner_p_type
     TYPE(preconditioner_type), POINTER :: preconditioner
  END TYPE preconditioner_p_type

CONTAINS

FUNCTION preconditioner_in_use(preconditioner)
    TYPE(preconditioner_type)                :: preconditioner
    LOGICAL                                  :: preconditioner_in_use

  preconditioner_in_use = .NOT. (preconditioner%in_use .EQ. 0)
END FUNCTION

SUBROUTINE init_preconditioner(preconditioner_env,para_env,blacs_env)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env

   NULLIFY(preconditioner_env%sparse_matrix)
   NULLIFY(preconditioner_env%fm)
   NULLIFY(preconditioner_env%subspace_evals)
   preconditioner_env%in_use=0
   preconditioner_env%para_env => para_env
   preconditioner_env%ctxt     => blacs_env
   CALL cp_para_env_retain(preconditioner_env%para_env)
   CALL cp_blacs_env_retain(preconditioner_env%ctxt)

END SUBROUTINE init_preconditioner

SUBROUTINE destroy_preconditioner(preconditioner_env)

    TYPE(preconditioner_type)                :: preconditioner_env

   IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
      CALL deallocate_matrix(preconditioner_env%sparse_matrix)
      NULLIFY(preconditioner_env%sparse_matrix)
   ENDIF

   IF (ASSOCIATED(preconditioner_env%fm)) THEN
      CALL cp_fm_release(preconditioner_env%fm)
   ENDIF
   IF (ASSOCIATED(preconditioner_env%subspace_evals)) THEN
      DEALLOCATE(preconditioner_env%subspace_evals)
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

SUBROUTINE make_preconditioner_single(preconditioner_env, precon_type, matrix_h, matrix_s, &
                                      matrix_t,mo_eigenvectors, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    CHARACTER(LEN=*)                         :: precon_type
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s, matrix_t
    TYPE(cp_fm_type), POINTER                :: mo_eigenvectors
    REAL(KIND=dp)                            :: energy_gap

    INTEGER                                  :: handle, k
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot

!----

  CALL timeset("make_preconditioner","I","",handle)

  CALL cp_fm_get_info(mo_eigenvectors,ncol_global=k) 
  ALLOCATE(eigenvalues_ot(k))
  ! XXXXXXXXXXXXXXXXxx makes MOs eigenstates, dangerous for non-equivalent guys...
  CALL calculate_subspace_eigenvalues(mo_eigenvectors,matrix_h,&
          eigenvalues_ot,para_env=preconditioner_env%para_env, &
                         do_rotation = .FALSE. ) 
  energy_homo=eigenvalues_ot(k)
  DEALLOCATE(eigenvalues_ot)

  SELECT CASE (precon_type)
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
  CASE ("FULL_KINETIC")
    preconditioner_env%in_use=6
    CALL make_full_kinetic(preconditioner_env, matrix_t, matrix_s, &
                                             energy_homo, energy_gap )

  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)

END SUBROUTINE make_preconditioner_single


!*******************************************************************************
! applies a previously created preconditioner to a full matrix
!*******************************************************************************
SUBROUTINE apply_preconditioner(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    INTEGER                                  :: handle

!----

  CALL timeset("apply_preconditioner","I","",handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
    CALL stop_program("apply preconditioner","No preconditioner in use")
  CASE (1,4) ! can share the apply since they are both sparse matrices
    CALL apply_sparse_diag(preconditioner_env, matrix_in, matrix_out)
  CASE (2,3,6)
    CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
  ! CASE (3,6)  ! this would be without the symmetrised form, is about 50 percent slower in parallel
  !   CALL apply_full_s_inverse(preconditioner_env, matrix_in, matrix_out)
  CASE DEFAULT
    CALL stop_program("apply preconditioner","implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)

END SUBROUTINE apply_preconditioner

!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! a sparse block diagonal approximation
SUBROUTINE apply_sparse_diag(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    INTEGER                                  :: k

!----

  CALL cp_fm_get_info(matrix_in,ncol_global=k)
  CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_in, &
                         matrix_out,k, para_env=preconditioner_env%para_env)

END SUBROUTINE apply_sparse_diag

SUBROUTINE make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                          energy_homo, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: iblock_col_h, iblock_col_s, &
                                                iblock_row, n, nblocks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s
    TYPE(real_block_node_type), POINTER      :: block_node_h, block_node_s

!----

  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN 
     CALL deallocate_matrix(preconditioner_env%sparse_matrix)
     NULLIFY(preconditioner_env%sparse_matrix)
  ENDIF

  CALL Replicate_matrix(matrix_s,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.)
  CALL get_matrix_info(matrix_s,nblock_row=nblocks)
  DO iblock_row=1,nblocks
     block_node_s => first_block_node(matrix_s,iblock_row)
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_s) .OR. ASSOCIATED(block_node_h))
        IF (.NOT.ASSOCIATED(block_node_s) .OR. .NOT.ASSOCIATED(block_node_h)) &
            CALL stop_program("qs preconditioner", ".NOT.ASSOCIATED")
         
        CALL get_block_node(block_node=block_node_s,&
                            block_col=iblock_col_s,&
                            BLOCK=block_s)

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            BLOCK=block_h)

        IF (iblock_col_h .NE. iblock_col_s) &
           CALL stop_program("qs preconditioner","iblock_col_h .ne. iblock_col_s")

        IF (iblock_col_s .EQ. iblock_row) THEN
           n=SIZE(block_s,1)
           ALLOCATE(block_pre(n,n))

           CALL make_local_block(block_h,block_s,block_pre,  &
                                                  energy_homo,energy_gap,3)

           CALL add_block_node(matrix=preconditioner_env%sparse_matrix,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col_s,&
                                 BLOCK=block_pre)

           DEALLOCATE(block_pre)
 
        ENDIF
        block_node_s => next_block_node(block_node_s)
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO
END SUBROUTINE make_sparse_diag
SUBROUTINE make_local_block(block_h,block_s,block_pre, &
                                 energy_homo,energy_gap, TYPE)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_s, block_pre
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    INTEGER                                  :: TYPE

    INTEGER                                  :: i, info, istat, liwork, &
                                                lwork, n
    INTEGER, DIMENSION(:), POINTER           :: iwork
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals, work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_buf1, block_chol, &
                                                block_evec

  n=SIZE(block_s,1)
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
 
  block_pre(:,:)=0.0_dp
  SELECT CASE (TYPE)
  CASE(1)
    DO i=1,n
      block_pre(i,i)=1.0_dp
    ENDDO
  CASE(2)
    DO i=1,n
      block_pre(i,i)=1.0_dp/MAX(energy_gap,block_h(i,i)-energy_homo)
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
     block_pre(:,:)=0.0_dp
     DO i=1,n
          block_pre(i,i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
     ENDDO
     ! K = V E V ^ T
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_dp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_dp,block_pre(1,1),n)
     ! inv(U) K inv(U)^T
     CALL DTRSM('L','U','N','N',n,n,1.0_dp,block_chol(1,1),n,block_pre(1,1),n)
     CALL DTRSM('R','U','T','N',n,n,1.0_dp,block_chol(1,1),n,block_pre(1,1),n)
  CASE(4)
     block_chol(:,:)=block_s(:,:)
     CALL DPOTRF('U',n,block_chol(1,1),n,info)
     block_pre(:,:)=0.0_dp
     DO i=1,n
          block_pre(i,i)=1.0_dp
     ENDDO
     CALL DTRSM('L','U','N','N',n,n,1.0_dp,block_chol(1,1),n,block_pre(1,1),n)
     CALL DTRSM('R','U','T','N',n,n,1.0_dp,block_chol(1,1),n,block_pre(1,1),n)
  CASE(5) ! like 3 but using s^-0.5 instead of the cholesky decomposition, and not transforming back
     block_evec(:,:)=block_s(:,:)
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     IF (info.ne.0) CALL stop_program("preconditioner","DSYEVD S")
     block_pre(:,:)=0.0_dp
     DO i=1,n
        block_pre(i,i)=1.0_dp/SQRT(evals(i))
     ENDDO
     ! block_pre is s^-0.5
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_dp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_dp,block_pre(1,1),n)
     ! transform H
     block_evec(:,:)=block_h(:,:)
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_dp,block_buf1(1,1),n)
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_pre(1,1),n,block_buf1(1,1),n, &
                                                  0.0_dp,block_evec(1,1),n)
     ! get evals and evecs
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     IF (info.ne.0) CALL stop_program("preconditioner","DSYEVD H")
     block_pre(:,:)=0.0_dp
     DO i=1,n
          block_pre(i,i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
     ENDDO
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_dp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_dp,block_pre(1,1),n)
  CASE(6) ! like 3 not doing any transformation with s before or after (supposedly done by the caller)
     block_evec(:,:)=block_h(:,:)
     ! get evals and evecs
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     IF (info.ne.0) CALL stop_program("preconditioner","DSYEVD H")
     block_pre(:,:)=0.0_dp
     DO i=1,n
          block_pre(i,i)=1.0_dp/MAX(evals(i),energy_gap)
     ENDDO
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_dp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_dp,block_pre(1,1),n)
     
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

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: iblock_col_h, iblock_col_s, &
                                                iblock_row, k, n, nblocks
    INTEGER, DIMENSION(:), POINTER           :: first_row, last_row, offset
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block_buffer
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s, &
                                                block_tmp
    TYPE(real_block_node_type), POINTER      :: block_node_h, block_node_s

!----

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
  CALL DCOPY(offset(nblocks+1),0.0_dp,0,block_buffer(1),1)

  DO iblock_row=1,nblocks
     block_node_s => first_block_node(matrix_s,iblock_row)
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_s) .OR. ASSOCIATED(block_node_h))
        IF (.NOT.ASSOCIATED(block_node_s) .OR. .NOT.ASSOCIATED(block_node_h)) &
            CALL stop_program("qs preconditioner", ".NOT.ASSOCIATED")

        CALL get_block_node(block_node=block_node_s,&
                            block_col=iblock_col_s,&
                            BLOCK=block_s)

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            BLOCK=block_h)

        IF (iblock_col_h .NE. iblock_col_s) &
           CALL stop_program("qs preconditioner","iblock_col_h .ne. iblock_col_s")

        IF (iblock_col_s .EQ. iblock_row) THEN
           n=SIZE(block_s,1)
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

  CALL DSCAL(offset(nblocks+1),0.2_dp,block_buffer(1),1)

  ! now every diagonal block is everywhere, form the product
  CALL Replicate_matrix(matrix_s,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.)

  DO iblock_row=1,nblocks
     block_node_s => first_block_node(matrix_s,iblock_row)
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_s) .OR. ASSOCIATED(block_node_h))

        CALL get_block_node(block_node=block_node_s,&
                            block_col=iblock_col_s,&
                            BLOCK=block_s)

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            BLOCK=block_h)

        n=SIZE(block_s,1)
        k=SIZE(block_s,2)
        ALLOCATE(block_pre(n,k))
        ALLOCATE(block_tmp(n,k))
        block_pre(:,:)=-block_s(:,:)
        CALL DGEMM('N','N',n,k,k,1.0_dp,block_pre(1,1),n, &
                                        block_buffer(offset(iblock_col_s)+1), &
                                        k,0.0_dp,block_tmp(1,1),n)
        CALL DGEMM('N','N',n,k,n,1.0_dp,block_buffer(offset(iblock_row)+1),n, &
                                        block_tmp(1,1), &
                                        n,0.0_dp,block_pre(1,1),n)
        IF (iblock_row .EQ. iblock_col_s) &
           CALL DAXPY(n*k,2.0_dp,block_buffer(offset(iblock_row)+1),1,block_pre(1,1),1)

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

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    INTEGER                                  :: k, n

!----

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  CALL cp_fm_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_in,0.0_dp,matrix_out)
END SUBROUTINE apply_full_single

SUBROUTINE make_full_single(preconditioner_env, matrix_h, matrix_s, &
                       energy_homo, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: i, n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_s

! ---
! ---

  NULLIFY(fm_h,fm_s,fm_struct_tmp,evals)

  IF (ASSOCIATED(preconditioner_env%fm)) THEN 
     CALL cp_fm_release(preconditioner_env%fm)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)
  ALLOCATE(evals(n))

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp, name="preconditioner")
  CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h")
  CALL cp_fm_create(fm_s,fm_struct_tmp, name="fm_s")
  CALL cp_fm_struct_release(fm_struct_tmp)

  CALL copy_sm_to_fm(matrix_h,fm_h)
  CALL copy_sm_to_fm(matrix_s,fm_s)
  CALL cp_fm_cholesky_decompose(fm_s)
  CALL cp_fm_cholesky_reduce(fm_h,fm_s)
  CALL cp_fm_syevd(fm_h,preconditioner_env%fm,evals)
  CALL cp_fm_cholesky_restore(preconditioner_env%fm,n,fm_s,fm_h,"SOLVE") 
  DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
  ENDDO
  CALL cp_fm_to_fm(fm_h,preconditioner_env%fm)
  CALL cp_fm_column_scale(preconditioner_env%fm,evals)
  CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,preconditioner_env%fm,fm_h,0.0_dp,fm_s)
  CALL cp_fm_to_fm(fm_s,preconditioner_env%fm)

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

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    INTEGER                                  :: k, n

!----

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)
  CALL cp_fm_symm('L','U',n,k,1.0_dp,preconditioner_env%fm, &
                        matrix_in, &
                        0.0_dp,matrix_out)

END SUBROUTINE apply_full_s_inverse
SUBROUTINE make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                               energy_homo, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: n
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work

! ---
! ---

  NULLIFY(fm_struct_tmp,fm_work)

  IF (ASSOCIATED(preconditioner_env%fm)) THEN
     CALL cp_fm_release(preconditioner_env%fm)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp, name="preconditioner")
  CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work")
  CALL cp_fm_struct_release(fm_struct_tmp)

  CALL copy_sm_to_fm(matrix_s,preconditioner_env%fm)
  CALL cp_fm_cholesky_decompose(preconditioner_env%fm)
  CALL cp_fm_cholesky_invert(preconditioner_env%fm)
  CALL cp_fm_upper_to_full(preconditioner_env%fm,fm_work) ! symmetrizing makes the apply_ up to 50 percent faster
  CALL cp_fm_release(fm_work)

END SUBROUTINE make_full_s_inverse



SUBROUTINE make_full_kinetic(preconditioner_env, matrix_t, matrix_s, &
                               energy_homo, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: shift
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work

! ---
! ---

  NULLIFY(fm_struct_tmp,fm_work)

  IF (ASSOCIATED(preconditioner_env%fm)) THEN
     CALL cp_fm_release(preconditioner_env%fm)
  ENDIF
  CALL get_matrix_info(matrix_t,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp, name="preconditioner")
  CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work")
  CALL cp_fm_struct_release(fm_struct_tmp)

  shift=MAX(-energy_homo,energy_gap)
  CALL add_matrices(matrix_t,1.0_dp,matrix_t,shift,matrix_s)

  CALL copy_sm_to_fm(matrix_t,preconditioner_env%fm)
  CALL cp_fm_cholesky_decompose(preconditioner_env%fm)
  CALL cp_fm_cholesky_invert(preconditioner_env%fm)
  CALL cp_fm_upper_to_full(preconditioner_env%fm,fm_work)
  CALL cp_fm_release(fm_work)

END SUBROUTINE make_full_kinetic

END MODULE preconditioner

