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
                                             cp_fm_upper_to_full,&
                                             cp_fm_scale_and_add
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
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type,&
                                             cp_fm_maxabsval
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE input_constants,                 ONLY: ot_precond_none,&
                                             ot_precond_full_all,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_s_inverse,&
                                             ot_precond_sparse_diag,&
                                             ot_precond_sparse_firstorder
  USE input_constants,                 ONLY: lr_precond_full_all,&
                                             lr_precond_s_inverse,&
                                             lr_precond_h0_inverse,&
                                             lr_precond_none
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                             make_basis
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             add_matrices,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             real_matrix_p_type,&
                                             replicate_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "preconditioner"

  PUBLIC  :: preconditioner_type, preconditioner_p_type
  PUBLIC  :: init_preconditioner
  PUBLIC  :: preconditioner_in_use
  PUBLIC  :: destroy_preconditioner
  PUBLIC  :: make_preconditioner
  PUBLIC  :: apply_preconditioner,apply_preconditioner_lr

  INTERFACE make_preconditioner
    MODULE PROCEDURE make_preconditioner_single,make_preconditioner_lr
  END INTERFACE

  PRIVATE :: make_sparse_diag, make_full_single, make_local_block

  TYPE preconditioner_type
!    PRIVATE
       TYPE(real_matrix_type),    POINTER :: sparse_matrix
       TYPE(cp_fm_type), POINTER :: fm
       TYPE(cp_para_env_type),  POINTER   :: para_env
       TYPE(cp_blacs_env_type), POINTER   :: ctxt
       INTEGER :: in_use
       REAL(KIND = dp), DIMENSION(:), POINTER :: occ_evals,full_evals
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
   NULLIFY(preconditioner_env%occ_evals)
   NULLIFY(preconditioner_env%full_evals)
   preconditioner_env%in_use=0
   preconditioner_env%para_env => para_env
   preconditioner_env%ctxt     => blacs_env
   CALL cp_para_env_retain(preconditioner_env%para_env)
   CALL cp_blacs_env_retain(preconditioner_env%ctxt)

END SUBROUTINE init_preconditioner

SUBROUTINE destroy_preconditioner(preconditioner_env)

    TYPE(preconditioner_type)                :: preconditioner_env
    INTEGER :: ip

   IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
      CALL deallocate_matrix(preconditioner_env%sparse_matrix)
      NULLIFY(preconditioner_env%sparse_matrix)
   ENDIF

   IF (ASSOCIATED(preconditioner_env%fm)) THEN
      CALL cp_fm_release(preconditioner_env%fm)
   ENDIF
   IF (ASSOCIATED(preconditioner_env%occ_evals)) THEN
      DEALLOCATE(preconditioner_env%occ_evals)
   ENDIF
   IF (ASSOCIATED(preconditioner_env%full_evals)) THEN
      DEALLOCATE(preconditioner_env%full_evals)
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
                                      matrix_t,mo_coeff, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    INTEGER                                  :: precon_type
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s, matrix_t
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    REAL(KIND=dp)                            :: energy_gap

    INTEGER                                  :: handle, k
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot

!----

  CALL timeset("make_preconditioner","I","",handle)

  CALL cp_fm_get_info(mo_coeff,ncol_global=k) 
  ALLOCATE(eigenvalues_ot(k))


  ! XXXXXXXXXXXXXXXX do not touch the initial MOs, could be harmful for either 
  !                  the case of non-equivalent MOs but also for the derivate
  ! we could already have all eigenvalues e.g. full_all and we could skip this
  ! to be optimised later.
  CALL calculate_subspace_eigenvalues(mo_coeff,matrix_h,&
          eigenvalues_ot,para_env=preconditioner_env%para_env, &
                         do_rotation = .FALSE.)

  energy_homo=eigenvalues_ot(k)

  SELECT CASE (precon_type)
  CASE (ot_precond_sparse_diag)
    preconditioner_env%in_use=1
    CALL make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                              energy_homo, energy_gap )
  CASE (ot_precond_full_single)
    preconditioner_env%in_use=2
    CALL make_full_single(preconditioner_env, preconditioner_env%fm,&
                          matrix_h, matrix_s, energy_homo, energy_gap )
  CASE (ot_precond_s_inverse)
    preconditioner_env%in_use=3
    CALL make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                                             energy_homo, energy_gap )
  CASE (ot_precond_sparse_firstorder)
    preconditioner_env%in_use=4
    CALL make_sparse_first_order(preconditioner_env, matrix_h, matrix_s, &
                                                 energy_homo, energy_gap )
  CASE (ot_precond_full_kinetic)
    preconditioner_env%in_use=6
    CALL make_full_kinetic(preconditioner_env, preconditioner_env%fm,&
                          matrix_t, matrix_s, energy_homo, energy_gap )
  CASE (ot_precond_full_all)
    preconditioner_env%in_use=19
    CALL make_full_all(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                       eigenvalues_ot, energy_gap)
  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT


  DEALLOCATE(eigenvalues_ot)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE make_preconditioner_single

SUBROUTINE make_preconditioner_lr(preconditioner_env, precon_type, ev_h0,&
                                  matrix_h, matrix_s, energy_gap ,error)

    TYPE(preconditioner_type),DIMENSION(:),&
      POINTER                                  :: preconditioner_env
    INTEGER                                    :: precon_type
    TYPE(cp_fm_p_type), DIMENSION(:),&
      POINTER                                  :: ev_h0
    TYPE(real_matrix_type), POINTER            :: matrix_s
    TYPE(real_matrix_p_type), DIMENSION(:),&
      POINTER                                  :: matrix_h
    REAL(KIND=dp)                              :: energy_gap
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: RoutineP = module_name//"make_preconditioner"

    INTEGER                                    :: handle, istat, k, x
    INTEGER                                    :: ispin,n,nspins,i
    LOGICAL                                    :: failure
    REAL(dp), DIMENSION(:), POINTER            :: unocc_evals,evals
    TYPE(cp_fm_struct_type), POINTER           :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                  :: fm_work, fm_h
    TYPE(cp_fm_type), POINTER                  :: matrix_x, matrix_sc0
    TYPE(cp_fm_type), POINTER                  :: matrix_xsc0

!----

  CALL timeset("make_preconditioner","I","",handle)

  NULLIFY(fm_struct_tmp,fm_work,matrix_sc0,matrix_x,matrix_xsc0)

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(preconditioner_env),cp_failure_level,routineP,error,failure)
  nspins = SIZE(preconditioner_env,1)

  DO ispin = 1, nspins
    IF (ASSOCIATED(preconditioner_env(ispin)%fm)) THEN
       CALL cp_fm_release(preconditioner_env(ispin)%fm)
    ENDIF
  END DO

  SELECT CASE (precon_type)
  CASE (lr_precond_s_inverse)

    preconditioner_env(1)%in_use=1
    CALL get_matrix_info(matrix_s,nrow=n)

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env(1)%ctxt, &
                             para_env=preconditioner_env(1)%para_env)
    CALL cp_fm_create(preconditioner_env(1)%fm,fm_struct_tmp, name="preconditioner")
    CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work")
    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL copy_sm_to_fm(matrix_s,preconditioner_env(1)%fm)
    CALL cp_fm_cholesky_decompose(preconditioner_env(1)%fm)
    CALL cp_fm_cholesky_invert(preconditioner_env(1)%fm)
    CALL cp_fm_upper_to_full(preconditioner_env(1)%fm,fm_work) ! symmetrizing makes the apply_ up to 50 percent faster
    CALL cp_fm_release(fm_work)
    IF(nspins==2)THEN
      preconditioner_env(2)%fm => preconditioner_env(1)%fm
      preconditioner_env(2)%in_use=1
    END IF

  CASE (lr_precond_h0_inverse)

    preconditioner_env(1)%in_use=1
    DO ispin = 1,nspins
      CALL get_matrix_info(matrix_h(ispin)%matrix,nrow=n)
      ALLOCATE(evals(n))
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env(ispin)%ctxt, &
                             para_env=preconditioner_env(ispin)%para_env)
      CALL cp_fm_create(preconditioner_env(ispin)%fm,fm_struct_tmp, name="preconditioner")
      CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h")
      CALL cp_fm_create(fm_work,fm_struct_tmp, name="fm_s")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL copy_sm_to_fm(matrix_h(ispin)%matrix,fm_h)
      CALL copy_sm_to_fm(matrix_s,fm_work)
      CALL cp_fm_cholesky_decompose(fm_work)
      CALL cp_fm_cholesky_reduce(fm_h,fm_work)
      CALL cp_fm_syevd(fm_h,preconditioner_env(ispin)%fm,evals)
      CALL cp_fm_cholesky_restore(preconditioner_env(ispin)%fm,n,fm_work,fm_h,"SOLVE")
! do i =1,n
!    write(*,*) i, evals(i)
! end do
! stop 'ev'
      k = SIZE(preconditioner_env(ispin)%occ_evals)
      DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-preconditioner_env(ispin)%occ_evals(k),energy_gap)
      ENDDO
      CALL cp_fm_to_fm(fm_h,preconditioner_env(ispin)%fm)
      CALL cp_fm_column_scale(preconditioner_env(ispin)%fm,evals)
      CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,preconditioner_env(ispin)%fm,fm_h,0.0_dp,fm_work)
      CALL cp_fm_to_fm(fm_work,preconditioner_env(ispin)%fm) 
      DEALLOCATE(evals)
      CALL cp_fm_release(fm_h)
      CALL cp_fm_release(fm_work)

    END DO

  CASE (lr_precond_full_all)

    CPPrecondition(ASSOCIATED(ev_h0),cp_failure_level,routineP,error,failure)
    CPPrecondition(nspins==SIZE(ev_h0,1),cp_failure_level,routineP,error,failure)

    preconditioner_env(1)%in_use=2
    DO ispin = 1,nspins
      CALL cp_fm_get_info(ev_h0(ispin)%matrix,nrow_global=n,ncol_global=k)
      !1) create vectors spanning the subspace orthogonal to ev_h0
      CALL cp_fm_create(matrix_sc0,ev_h0(ispin)%matrix%matrix_struct,name="sc0")
      x=n-k
      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=x, &
                             context=preconditioner_env(ispin)%ctxt, &
                             para_env=preconditioner_env(ispin)%para_env)
      CALL cp_fm_create(matrix_x,fm_struct_tmp,name="matrix_x")
      CALL cp_fm_create(fm_work,fm_struct_tmp,name="fm_work")
      CALL cp_fm_struct_release(fm_struct_tmp) 

      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=x, &
                               context=preconditioner_env(ispin)%ctxt, &
                               para_env=preconditioner_env(ispin)%para_env)
      CALL cp_fm_create(matrix_xsc0,fm_struct_tmp,name="xsc0")
      CALL cp_fm_struct_release(fm_struct_tmp)

      !2) orthogonalize
      CALL cp_fm_init_random(matrix_x)
      CALL cp_sm_fm_multiply(matrix_s,ev_h0(ispin)%matrix,matrix_sc0,k)
      CALL cp_fm_gemm('T','N',k,x,n,1.0_dp,matrix_sc0,matrix_x,0.0_dp,matrix_xsc0)
      CALL cp_fm_gemm('N','N',n,x,k,1.0_dp,ev_h0(ispin)%matrix,matrix_xsc0,0.0_dp,fm_work)
      ! with respect C0
      CALL cp_fm_scale_and_add(1.0_dp,matrix_x,-1.0_dp,fm_work)
      CALL cp_sm_fm_multiply(matrix_s,matrix_x,fm_work,x)
      ! among themselves
      CALL make_basis(matrix_x,x, ortho=fm_work,otype="SV")

      !3) make eigenvectors and store eigenvalues
      ALLOCATE(preconditioner_env(ispin)%full_evals(n),STAT=istat) 
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      preconditioner_env(ispin)%full_evals(1:k)=preconditioner_env(ispin)%occ_evals(1:k)

      ALLOCATE(unocc_evals(x),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL calculate_subspace_eigenvalues(matrix_x,matrix_h(ispin)%matrix, unocc_evals , &
                 para_env=preconditioner_env(ispin)%para_env, do_rotation = .TRUE.)
      preconditioner_env(ispin)%full_evals(k+1:n)=unocc_evals
! do i =1,n
!    write(*,*) i, preconditioner_env(ispin)%full_evals(i)
! end do
! stop 'ev'
      DEALLOCATE(unocc_evals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      preconditioner_env(ispin)%energy_gap=energy_gap

      !4) Clean
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(matrix_xsc0)
      CALL cp_fm_release(matrix_sc0)
      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                               context=preconditioner_env(ispin)%ctxt, &
                               para_env=preconditioner_env(ispin)%para_env)
      CALL cp_fm_create(preconditioner_env(ispin)%fm,fm_struct_tmp,name="preconditioner_env%fm")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_to_fm(ev_h0(ispin)%matrix,preconditioner_env(ispin)%fm,k,1,1)
      CALL cp_fm_to_fm(matrix_x,preconditioner_env(ispin)%fm,x,1,k+1)
      CALL cp_fm_release(matrix_x)

    END DO 


  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)

END SUBROUTINE make_preconditioner_lr

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
  CASE (19)
    CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
  CASE DEFAULT
    CALL stop_program("apply preconditioner","implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)

END SUBROUTINE apply_preconditioner

SUBROUTINE apply_preconditioner_lr(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    INTEGER                                  :: handle

!----

  CALL timeset("apply_preconditioner","I","",handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
    CALL stop_program("apply preconditioner","No preconditioner in use")
  CASE (1)
    CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
  CASE (2)
    CALL apply_full_all(preconditioner_env, matrix_in, matrix_out) 
  CASE DEFAULT
    CALL stop_program("apply preconditioner","implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)



END SUBROUTINE apply_preconditioner_lr
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
                         matrix_out,k)

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

!
! generates a preconditioner for each eigenstate i by keeping the factorized form
! U diag( something i ) U^T
! an important trick is to only precondition in the subspace orthogonal to c0.
! not only is it the only part that matters, it also simplifies the computation of
! the lagrangian multipliers in the OT minimization  (i.e. if the c0 here is different 
! from the c0 used in the OT setup, there will be a bug).
!
! the assumption here is that the c0 are already ritz states of (h,s)
!
SUBROUTINE make_full_all(preconditioner_env, matrix_c0, matrix_h, matrix_s, c0_evals, energy_gap)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_gap
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals

    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    INTEGER                                  :: k,x,n
    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_sc0, matrix_xsc0, matrix_tmp
    REAL(KIND=dp)                            :: a
    REAL(KIND=dp), DIMENSION(:), POINTER     :: unocc_evals

    IF (ASSOCIATED(preconditioner_env%fm)) CALL cp_fm_release(preconditioner_env%fm)
    CALL cp_fm_get_info(matrix_c0,nrow_global=n,ncol_global=k)

    ! 1) create the vectors spanning the subspace orthogonal to c0, these will be the basis of the preconditioner
    ! 1a) create matrix_c0, matrix_x, matrix_xsc0
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0")
    x=n-k
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=x, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env) 
    CALL cp_fm_create(matrix_x,fm_struct_tmp,name="matrix_x")
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp,name="matrix_tmp")
    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=x, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
    CALL cp_fm_create(matrix_xsc0,fm_struct_tmp,name="xsc0")
    CALL cp_fm_struct_release(fm_struct_tmp)

    ! 1b) orthogonalise a random x to c0 
    CALL cp_fm_init_random(matrix_x)
    CALL cp_sm_fm_multiply(matrix_s,matrix_c0,matrix_sc0,k)
    CALL cp_fm_gemm('T','N',k,x,n,1.0_dp,matrix_sc0,matrix_x,0.0_dp,matrix_xsc0) 
    CALL cp_fm_gemm('N','N',n,x,k,1.0_dp,matrix_c0,matrix_xsc0,0.0_dp,matrix_tmp)
    CALL cp_fm_scale_and_add(1.0_dp,matrix_x,-1.0_dp,matrix_tmp)
    CALL cp_sm_fm_multiply(matrix_s,matrix_x,matrix_tmp,x)
    CALL make_basis(matrix_x,x, ortho=matrix_tmp,otype="SV")

    !XXXXXXXX an orthogonality check here
    !CALL cp_fm_gemm('T','N',k,x,n,1.0_dp,matrix_sc0,matrix_x,0.0_dp,matrix_xsc0)
    !CALL cp_fm_maxabsval(matrix_xsc0,a)
    !write(6,*) "cp_fm_maxabsval(matrix_xsc0) after ortho",a

    ! 2) make matrix_x eigenvectors of Hks in their subspace, and store all evals we have
    ALLOCATE(preconditioner_env%full_evals(n))
    ALLOCATE(preconditioner_env%occ_evals(k))
    preconditioner_env%occ_evals=c0_evals
    preconditioner_env%full_evals(1:k)=c0_evals ! notice, this choice causes the preconditioner to be constant when applied to sc0 (see apply_full_all) 

    ALLOCATE(unocc_evals(x)) 
    CALL calculate_subspace_eigenvalues(matrix_x,matrix_h, unocc_evals , &
                 para_env=preconditioner_env%para_env, do_rotation = .TRUE.)
    preconditioner_env%full_evals(k+1:n)=unocc_evals
    DEALLOCATE(unocc_evals)

    preconditioner_env%energy_gap=energy_gap

    ! 3) Clean up and asm the full preconditioner matrix
    CALL cp_fm_release(matrix_tmp)
    CALL cp_fm_release(matrix_xsc0)
    CALL cp_fm_release(matrix_sc0)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env) 
    CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp,name="preconditioner_env%fm")
    CALL cp_fm_struct_release(fm_struct_tmp)
    CALL cp_fm_to_fm(matrix_c0,preconditioner_env%fm,k,1,1)
    CALL cp_fm_to_fm(matrix_x,preconditioner_env%fm,x,1,k+1)
    CALL cp_fm_release(matrix_x)

END SUBROUTINE make_full_all

!
! the corresponding apply_full_all uses the decomposed form to apply the preconditioner
!

SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out

    INTEGER                                  :: i, j, k, n
    TYPE(cp_fm_type), POINTER                :: matrix_tmp

    INTEGER, DIMENSION(:), POINTER :: row_indices
    INTEGER, DIMENSION(:), POINTER :: col_indices
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: local_data
    INTEGER :: nrow_local,ncol_local
    REAL(KIND=dp)                            :: dum

!----

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k)

  CALL cp_fm_create(matrix_tmp,matrix_in%matrix_struct,name="apply_full_all")
  CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data)

  ! 
  CALL cp_fm_gemm('T','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_in,0.0_dp,matrix_tmp)

  ! do the right scaling
  DO j=1,ncol_local
  DO i=1,nrow_local
     dum=1.0_dp/MAX(preconditioner_env%energy_gap, &
             preconditioner_env%full_evals(row_indices(i))-preconditioner_env%occ_evals(col_indices(j)))
     local_data(i,j)=local_data(i,j)*dum
  ENDDO
  ENDDO

  ! mult back
  CALL cp_fm_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_tmp,0.0_dp,matrix_out)

  CALL cp_fm_release(matrix_tmp)

END SUBROUTINE apply_full_all

SUBROUTINE make_full_single(preconditioner_env, fm, matrix_h, matrix_s, &
                       energy_homo, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: i, n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_s

! ---
! ---

  NULLIFY(fm_h,fm_s,fm_struct_tmp,evals)

  IF (ASSOCIATED(fm)) THEN 
     CALL cp_fm_release(fm)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)
  ALLOCATE(evals(n))

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner")
  CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h")
  CALL cp_fm_create(fm_s,fm_struct_tmp, name="fm_s")
  CALL cp_fm_struct_release(fm_struct_tmp)

  CALL copy_sm_to_fm(matrix_h,fm_h)
  CALL copy_sm_to_fm(matrix_s,fm_s)
  CALL cp_fm_cholesky_decompose(fm_s)
  CALL cp_fm_cholesky_reduce(fm_h,fm_s)
  CALL cp_fm_syevd(fm_h,fm,evals)
  CALL cp_fm_cholesky_restore(fm,n,fm_s,fm_h,"SOLVE") 
  DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
  ENDDO
  CALL cp_fm_to_fm(fm_h,fm)
  CALL cp_fm_column_scale(fm,evals)
  CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,fm,fm_h,0.0_dp,fm_s)
  CALL cp_fm_to_fm(fm_s,fm)

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



SUBROUTINE make_full_kinetic(preconditioner_env, fm, matrix_t, matrix_s, &
                               energy_homo, energy_gap )

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), POINTER          :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: shift
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work

! ---
! ---

  NULLIFY(fm_struct_tmp,fm_work)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm)
  ENDIF
  CALL get_matrix_info(matrix_t,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner")
  CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work")
  CALL cp_fm_struct_release(fm_struct_tmp)

  shift=MAX(-energy_homo,energy_gap)
  CALL add_matrices(matrix_t,1.0_dp,matrix_t,shift,matrix_s)

  CALL copy_sm_to_fm(matrix_t,fm)
  CALL cp_fm_cholesky_decompose(fm)
  CALL cp_fm_cholesky_invert(fm)
  CALL cp_fm_upper_to_full(fm,fm_work)
  CALL cp_fm_release(fm_work)

END SUBROUTINE make_full_kinetic

END MODULE preconditioner

