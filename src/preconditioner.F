!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
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
                                             cp_fm_scale_and_add,&
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
                                             cp_fm_get_diag,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE input_constants,                 ONLY: &
       lr_precond_full_all, lr_precond_h0_inverse, lr_precond_s_inverse, &
       ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_s_inverse, &
       ot_precond_sparse_diag, ot_precond_sparse_firstorder
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                             make_basis_sm
  USE sparse_matrix_types,             ONLY: &
       add_block_node, add_matrices, deallocate_matrix, first_block_node, &
       get_block_node, get_matrix_info, next_block_node, &
       real_block_node_type, real_matrix_p_type, real_matrix_type, &
       replicate_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'

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

SUBROUTINE init_preconditioner(preconditioner_env,para_env,blacs_env, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

   NULLIFY(preconditioner_env%sparse_matrix)
   NULLIFY(preconditioner_env%fm)
   NULLIFY(preconditioner_env%occ_evals)
   NULLIFY(preconditioner_env%full_evals)
   preconditioner_env%in_use=0
   preconditioner_env%para_env => para_env
   preconditioner_env%ctxt     => blacs_env
   CALL cp_para_env_retain(preconditioner_env%para_env,error=error)
   CALL cp_blacs_env_retain(preconditioner_env%ctxt,error=error)

END SUBROUTINE init_preconditioner

SUBROUTINE destroy_preconditioner(preconditioner_env, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_error_type), INTENT(inout)       :: error

   IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
      CALL deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
      NULLIFY(preconditioner_env%sparse_matrix)
   ENDIF

   IF (ASSOCIATED(preconditioner_env%fm)) THEN
      CALL cp_fm_release(preconditioner_env%fm,error=error)
   ENDIF
   IF (ASSOCIATED(preconditioner_env%occ_evals)) THEN
      DEALLOCATE(preconditioner_env%occ_evals)
   ENDIF
   IF (ASSOCIATED(preconditioner_env%full_evals)) THEN
      DEALLOCATE(preconditioner_env%full_evals)
   ENDIF

   CALL cp_para_env_release(preconditioner_env%para_env,error=error)
   CALL cp_blacs_env_release(preconditioner_env%ctxt,error=error)

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
                                      matrix_t,mo_coeff, energy_gap ,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    INTEGER                                  :: precon_type
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s, matrix_t
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot

!----

  CALL timeset(routineN,"I","",handle)

  CALL cp_fm_get_info(mo_coeff,ncol_global=k,error=error)

  ALLOCATE(eigenvalues_ot(k))


  ! XXXXXXXXXXXXXXXX do not touch the initial MOs, could be harmful for either
  !                  the case of non-equivalent MOs but also for the derivate
  ! we could already have all eigenvalues e.g. full_all and we could skip this
  ! to be optimised later.
  ! one flaw is that not all SCF methods (i.e. that go over mo_derivs directly)
  ! have a 'valid' matrix_h... (we even don't know what evals are in that case)
  CALL calculate_subspace_eigenvalues(mo_coeff,matrix_h,&
          eigenvalues_ot,para_env=preconditioner_env%para_env, &
                         do_rotation = .FALSE.,error=error)

  ! deal with zero occupied orbitals
  IF (k>0) THEN
     energy_homo=eigenvalues_ot(k)
  ELSE
     energy_homo=0.0_dp
  ENDIF

  SELECT CASE (precon_type)
  CASE (ot_precond_sparse_diag)
    preconditioner_env%in_use=ot_precond_sparse_diag
    CALL make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                              energy_homo, energy_gap ,error=error)
  CASE (ot_precond_full_single)
    preconditioner_env%in_use=ot_precond_full_single
    CALL make_full_single(preconditioner_env, preconditioner_env%fm,&
                          matrix_h, matrix_s, energy_homo, energy_gap ,error=error)
  CASE (ot_precond_s_inverse)
    preconditioner_env%in_use=ot_precond_s_inverse
    CALL make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                                             energy_homo, energy_gap ,error=error)
  CASE (ot_precond_sparse_firstorder)
    preconditioner_env%in_use=ot_precond_sparse_firstorder
    CALL make_sparse_first_order(preconditioner_env, matrix_h, matrix_s, &
                                                 energy_homo, energy_gap ,error=error)
  CASE (ot_precond_full_kinetic)
    preconditioner_env%in_use=ot_precond_full_kinetic
    CALL make_full_kinetic(preconditioner_env, preconditioner_env%fm,&
                          matrix_t, matrix_s, energy_homo, energy_gap ,error=error)
  CASE (ot_precond_full_single_inverse)
    preconditioner_env%in_use=ot_precond_full_single_inverse
    CALL make_full_single_inverse(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                       eigenvalues_ot, energy_gap,error=error)
  CASE (ot_precond_full_all)
    preconditioner_env%in_use=ot_precond_full_all
     CALL make_full_all(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                        eigenvalues_ot, energy_gap,error=error)
  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT


  DEALLOCATE(eigenvalues_ot)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE make_preconditioner_single

SUBROUTINE make_preconditioner_lr(preconditioner_env, precon_type, ev_h0,&
                                  matrix_h, matrix_s, energy_gap ,error)

    TYPE(preconditioner_type), &
      DIMENSION(:), POINTER                  :: preconditioner_env
    INTEGER                                  :: precon_type
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: ev_h0
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h
    TYPE(real_matrix_type), POINTER          :: matrix_s
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner_lr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, istat, k, &
                                                n, nspins, x
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: evals, unocc_evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_work, matrix_sc0, &
                                                matrix_x, matrix_xsc0

!----

  CALL timeset(routineN,"I","",handle)

  NULLIFY(fm_struct_tmp,fm_work,matrix_sc0,matrix_x,matrix_xsc0)

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(preconditioner_env),cp_failure_level,routineP,error,failure)
  nspins = SIZE(preconditioner_env,1)

  DO ispin = 1, nspins
    IF (ASSOCIATED(preconditioner_env(ispin)%fm)) THEN
       CALL cp_fm_release(preconditioner_env(ispin)%fm,error=error)
    ENDIF
  END DO

  SELECT CASE (precon_type)
  CASE (lr_precond_s_inverse)
       preconditioner_env(1)%in_use=lr_precond_s_inverse
       CALL get_matrix_info(matrix_s,nrow=n)
       CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                                context=preconditioner_env(1)%ctxt, &
                                para_env=preconditioner_env(1)%para_env,error=error)
       CALL cp_fm_create(preconditioner_env(1)%fm,fm_struct_tmp, name="preconditioner",error=error)
       CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       CALL copy_sm_to_fm(matrix_s,preconditioner_env(1)%fm,error=error)
       CALL cp_fm_cholesky_decompose(preconditioner_env(1)%fm,error=error)
       CALL cp_fm_cholesky_invert(preconditioner_env(1)%fm,error=error)
       CALL cp_fm_upper_to_full(preconditioner_env(1)%fm,fm_work,error=error) ! symmetrizing makes the apply_ up to 50 percent faster
       CALL cp_fm_release(fm_work,error=error)
       DO ispin = 2,nspins
          preconditioner_env(ispin)%in_use=lr_precond_s_inverse
          CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                                   context=preconditioner_env(ispin)%ctxt, &
                                   para_env=preconditioner_env(ispin)%para_env,error=error)
          CALL cp_fm_create(preconditioner_env(ispin)%fm,fm_struct_tmp, name="preconditioner",error=error)
          CALL cp_fm_struct_release(fm_struct_tmp,error=error)
          CALL cp_fm_to_fm(preconditioner_env(1)%fm,preconditioner_env(ispin)%fm,error=error)
       ENDDO
!    IF(nspins==2)THEN
!      preconditioner_env(2)%fm => preconditioner_env(1)%fm
!      preconditioner_env(2)%in_use=lr_precond_s_inverse
!    END IF
  CASE (lr_precond_h0_inverse)

    DO ispin = 1,nspins
       preconditioner_env(ispin)%in_use=lr_precond_h0_inverse
      CALL get_matrix_info(matrix_h(ispin)%matrix,nrow=n)
      ALLOCATE(evals(n))
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env(ispin)%ctxt, &
                             para_env=preconditioner_env(ispin)%para_env,error=error)
      CALL cp_fm_create(preconditioner_env(ispin)%fm,fm_struct_tmp, name="preconditioner",error=error)
      CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h",error=error)
      CALL cp_fm_create(fm_work,fm_struct_tmp, name="fm_s",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      CALL copy_sm_to_fm(matrix_h(ispin)%matrix,fm_h,error=error)
      CALL copy_sm_to_fm(matrix_s,fm_work,error=error)
      CALL cp_fm_cholesky_decompose(fm_work,error=error)
      CALL cp_fm_cholesky_reduce(fm_h,fm_work,error=error)
      CALL cp_fm_syevd(fm_h,preconditioner_env(ispin)%fm,evals,error=error)
      CALL cp_fm_cholesky_restore(preconditioner_env(ispin)%fm,n,fm_work,fm_h,"SOLVE",error=error)
! do i =1,n
!    write(*,*) i, evals(i)
! end do
! stop 'ev'
      k = SIZE(preconditioner_env(ispin)%occ_evals)
      DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-preconditioner_env(ispin)%occ_evals(k),energy_gap)
      ENDDO
      CALL cp_fm_to_fm(fm_h,preconditioner_env(ispin)%fm,error=error)
      CALL cp_fm_column_scale(preconditioner_env(ispin)%fm,evals)
      CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,preconditioner_env(ispin)%fm,fm_h,0.0_dp,fm_work,error=error)
      CALL cp_fm_to_fm(fm_work,preconditioner_env(ispin)%fm,error=error)
      DEALLOCATE(evals)
      CALL cp_fm_release(fm_h,error=error)
      CALL cp_fm_release(fm_work,error=error)

    END DO

  CASE (lr_precond_full_all)

    CPPrecondition(ASSOCIATED(ev_h0),cp_failure_level,routineP,error,failure)
    CPPrecondition(nspins==SIZE(ev_h0,1),cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspins
       preconditioner_env(ispin)%in_use=lr_precond_full_all
      CALL cp_fm_get_info(ev_h0(ispin)%matrix,nrow_global=n,ncol_global=k,error=error)
      !1) create vectors spanning the subspace orthogonal to ev_h0
      CALL cp_fm_create(matrix_sc0,ev_h0(ispin)%matrix%matrix_struct,name="sc0",error=error)
      x=n-k
      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=x, &
                             context=preconditioner_env(ispin)%ctxt, &
                             para_env=preconditioner_env(ispin)%para_env,error=error)
      CALL cp_fm_create(matrix_x,fm_struct_tmp,name="matrix_x",error=error)
      CALL cp_fm_create(fm_work,fm_struct_tmp,name="fm_work",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=x, &
                               context=preconditioner_env(ispin)%ctxt, &
                               para_env=preconditioner_env(ispin)%para_env,error=error)
      CALL cp_fm_create(matrix_xsc0,fm_struct_tmp,name="xsc0",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      !2) orthogonalize
      CALL cp_fm_init_random(matrix_x,error=error)
      CALL cp_sm_fm_multiply(matrix_s,ev_h0(ispin)%matrix,matrix_sc0,k,error=error)
      CALL cp_fm_gemm('T','N',k,x,n,1.0_dp,matrix_sc0,matrix_x,0.0_dp,matrix_xsc0,error=error)
      CALL cp_fm_release(matrix_sc0,error=error)
      CALL cp_fm_gemm('N','N',n,x,k,1.0_dp,ev_h0(ispin)%matrix,matrix_xsc0,0.0_dp,fm_work,error=error)
      CALL cp_fm_release(matrix_xsc0,error=error)
      ! with respect C0
      CALL cp_fm_scale_and_add(1.0_dp,matrix_x,-1.0_dp,fm_work,error=error)
      CALL cp_fm_release(fm_work,error=error)
      ! among themselves
      CALL make_basis_sm(matrix_x,x, matrix_s,error=error)

      !3) make eigenvectors and store eigenvalues
      ALLOCATE(preconditioner_env(ispin)%full_evals(n),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      preconditioner_env(ispin)%full_evals(1:k)=preconditioner_env(ispin)%occ_evals(1:k)

      ALLOCATE(unocc_evals(x),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL calculate_subspace_eigenvalues(matrix_x,matrix_h(ispin)%matrix, unocc_evals , &
                 para_env=preconditioner_env(ispin)%para_env, do_rotation = .TRUE.,error=error)
      preconditioner_env(ispin)%full_evals(k+1:n)=unocc_evals
! do i =1,n
!    write(*,*) i, preconditioner_env(ispin)%full_evals(i)
! end do
! stop 'ev'
      DEALLOCATE(unocc_evals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      preconditioner_env(ispin)%energy_gap=energy_gap

      !4) Clean
      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                               context=preconditioner_env(ispin)%ctxt, &
                               para_env=preconditioner_env(ispin)%para_env,error=error)
      CALL cp_fm_create(preconditioner_env(ispin)%fm,fm_struct_tmp,name="preconditioner_env%fm",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
      CALL cp_fm_to_fm(ev_h0(ispin)%matrix,preconditioner_env(ispin)%fm,k,1,1)
      CALL cp_fm_to_fm(matrix_x,preconditioner_env(ispin)%fm,x,1,k+1)
      CALL cp_fm_release(matrix_x,error=error)

    END DO


  CASE DEFAULT
    CALL stop_program("make preconditioner","Type not implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)

END SUBROUTINE make_preconditioner_lr

!*******************************************************************************
! applies a previously created preconditioner to a full matrix
!*******************************************************************************
SUBROUTINE apply_preconditioner(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

!----

  CALL timeset(routineN,"I","",handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
    CALL stop_program("apply preconditioner","No preconditioner in use")
  CASE (ot_precond_sparse_diag,ot_precond_sparse_firstorder)
    CALL apply_sparse_diag(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_single,ot_precond_s_inverse,ot_precond_full_kinetic,ot_precond_full_single_inverse)
      CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_all)
    CALL apply_full_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE DEFAULT
    CALL stop_program("apply preconditioner","implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)

END SUBROUTINE apply_preconditioner

SUBROUTINE apply_preconditioner_lr(preconditioner_env, matrix_in, matrix_out,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_lr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

!----

  CALL timeset(routineN,"I","",handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
    CALL stop_program("apply_preconditioner_lr","No preconditioner in use")
  CASE (lr_precond_s_inverse,lr_precond_h0_inverse)
    CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (lr_precond_full_all)
    CALL apply_full_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE DEFAULT
    CALL stop_program("apply_preconditioner_lr","implemented")
  END SELECT

  CALL timestop(0.0_dp,handle)



END SUBROUTINE apply_preconditioner_lr
!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! a sparse block diagonal approximation
SUBROUTINE apply_sparse_diag(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: k

!----

  CALL cp_fm_get_info(matrix_in,ncol_global=k,error=error)
  CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_in, &
                         matrix_out,k,error=error)

END SUBROUTINE apply_sparse_diag

SUBROUTINE make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                          energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iblock_col_h, iblock_col_s, &
                                                iblock_row, n, nblocks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s
    TYPE(real_block_node_type), POINTER      :: block_node_h, block_node_s

!----

  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
     CALL deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
     NULLIFY(preconditioner_env%sparse_matrix)
  ENDIF

  CALL Replicate_matrix(matrix_s,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.,error=error)
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
                                 BLOCK=block_pre,error=error)

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
                                   energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iblock_col_h, iblock_col_s, &
                                                iblock_row, k, n, nblocks
    INTEGER, DIMENSION(:), POINTER           :: first_row, last_row, offset
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block_buffer
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s, &
                                                block_tmp
    TYPE(real_block_node_type), POINTER      :: block_node_h, block_node_s

!----

  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
     CALL deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
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
  CALL Replicate_matrix(matrix_s,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.,error=error)

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

        CALL add_block_node(preconditioner_env%sparse_matrix,iblock_row,iblock_col_s,block_pre,error=error)

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
SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: k, n

!----

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)
  CALL cp_fm_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_in,0.0_dp,matrix_out,error=error)
END SUBROUTINE apply_full_single

!!****f* preconditioner/make_full_single_inverse *
!!
!!   NAME
!!     make_full_single_inverse
!!
!!   FUNCTION
!!     generates a preconditioner by cholesky inverting H-lambda S+(SC)shifts(SC)^T
!!
!!   NOTES
!!     this might fail if the initial guess is bad, or if the system has 'holes'
!!
!!   INPUTS
!!    - matrix_c0 must be already rotated correctly.
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!     10.2006 made more robust [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE make_full_single_inverse(preconditioner_env, matrix_c0, matrix_h, matrix_s, c0_evals, energy_gap,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: eval_shift = 5.0_dp , &
                                                fudge_factor = 2.0_dp 

    INTEGER                                  :: k, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: error_estimate, &
                                                preconditioner_shift
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_shift, diag, evals
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_hc0, matrix_s1, &
                                                matrix_sc0, matrix_shc0, &
                                                matrix_tmp, matrix_tmp2

! arbitrary upshift of the occupied evals
! fudge factor for taking the error estimate into account

    CALL cp_fm_get_info(matrix_c0,nrow_global=n,ncol_global=k,error=error)

    IF (ASSOCIATED(preconditioner_env%fm)) CALL cp_fm_release(preconditioner_env%fm,error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp,name="preconditioner_env%fm",error=error)
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp, name="preconditioner matrix_tmp",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! first try to get a ritz error estimate
    ! 0) cholesky decompose the overlap matrix, if this fails the basis is singular,
    !    more than EPS_DEFAULT
    CALL copy_sm_to_fm(matrix_s,matrix_tmp,error=error)
    CALL cp_fm_cholesky_decompose(matrix_tmp,error=error)

    ! get the error estimate
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0",error=error)
    CALL cp_sm_fm_multiply(matrix_s,matrix_c0,matrix_sc0,k,error=error)
    CALL cp_fm_create(matrix_hc0,matrix_c0%matrix_struct,name="hc0",error=error)
    CALL cp_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

    CALL cp_fm_create(matrix_shc0,matrix_c0%matrix_struct,name="shc0",error=error)
    CALL cp_fm_cholesky_restore(matrix_hc0,k,matrix_tmp,matrix_shc0,"SOLVE",transa="T",error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ! since we only use diagonal elements this is a bit of a waste
    CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,matrix_shc0,matrix_shc0,0.0_dp,matrix_s1,error=error)
    ALLOCATE(diag(k))
    CALL cp_fm_get_diag(matrix_s1,diag,error=error)
    error_estimate=MAXVAL(SQRT(ABS(diag-c0_evals**2)))
    DEALLOCATE(diag)
    CALL cp_fm_release(matrix_s1,error=error)
    CALL cp_fm_release(matrix_shc0,error=error)
    CALL cp_fm_release(matrix_hc0,error=error)

    ! shift up the occupied subspace eigenvalues
    ALLOCATE(c0_shift(k))
    c0_shift=SQRT(-(c0_evals-c0_evals(k))+eval_shift)
    CALL cp_fm_column_scale(matrix_sc0,c0_shift)
    CALL cp_fm_gemm('N','T',n,n,k,1.0_dp,matrix_sc0,matrix_sc0,0.0_dp,preconditioner_env%fm,error=error)
    CALL cp_fm_release(matrix_sc0,error=error)
    DEALLOCATE(c0_shift)

    ! get H added to the shift
    CALL copy_sm_to_fm(matrix_h,matrix_tmp,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,preconditioner_env%fm,1.0_dp,matrix_tmp,error=error)
    CALL copy_sm_to_fm(matrix_s,matrix_tmp,error=error)

    ! preconditioner shift, we target the middle of the occupied spectrum, and taking into account the error_estimate
    ! write(6,*) "Error estimate = ",error_estimate
    preconditioner_shift=-(MINVAL(c0_evals)+ MAXVAL(c0_evals))/2.0_dp + &
                           error_estimate*fudge_factor
    CALL copy_sm_to_fm(matrix_s,matrix_tmp,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,preconditioner_env%fm,preconditioner_shift,matrix_tmp,error=error)

    ! check evals
    IF (.FALSE.) THEN
       CALL cp_fm_to_fm(preconditioner_env%fm,matrix_tmp,error=error)
       CALL cp_fm_create(matrix_tmp2,matrix_tmp%matrix_struct,name="matrix_tmp2",error=error)
       ALLOCATE(evals(n))
       CALL cp_fm_syevd(matrix_tmp,matrix_tmp2,evals,error)
       CALL cp_fm_release(matrix_tmp2,error=error)
       WRITE(6,*) "evals ",evals
       DEALLOCATE(evals)
    ENDIF

    ! the preconditioner ought now to be symmetric positive definite (provided C0 & system are OK)
    ! we do catch the error and try to fix things
    CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
    CALL cp_fm_cholesky_decompose(preconditioner_env%fm,error=sub_error)
    failure = .FALSE. ; CALL cp_error_check(sub_error,failure)
    CALL cp_error_dealloc_ref(sub_error)

    ! check if we need to replace the preconditioner by something simple
    ! just to keep going
    IF (failure) THEN
       ! just the inverse of the overlap matrix right now    
       CALL copy_sm_to_fm(matrix_s,preconditioner_env%fm,error=error)
       CALL cp_fm_cholesky_decompose(preconditioner_env%fm,error=error)
    ENDIF

    ! we need to bring it in the final shape (i.e. full matrix, explicit inverse)
    CALL cp_fm_cholesky_invert(preconditioner_env%fm,error=error)
    CALL cp_fm_upper_to_full(preconditioner_env%fm,matrix_tmp,error=error)
    CALL cp_fm_release(matrix_tmp,error=error)

END SUBROUTINE make_full_single_inverse

!!****f* preconditioner/make_full_all *
!!
!!   NAME
!!     make_full_all
!!
!!   FUNCTION
!!     generates a state by state preconditioner based on the full hamiltonian matrix 
!!
!!   NOTES
!!     includes error estimate on the hamiltonian matrix to result in a stable preconditioner
!!     a preconditioner for each eigenstate i is generated by keeping the factorized form
!!     U diag( something i ) U^T. It is important to only precondition in the subspace orthogonal to c0.
!!     not only is it the only part that matters, it also simplifies the computation of
!!     the lagrangian multipliers in the OT minimization  (i.e. if the c0 here is different
!!     from the c0 used in the OT setup, there will be a bug).
!!
!!   INPUTS
!!     - energy_gap:  should be a slight underestimate of the physical energy gap for almost all systems
!!     the c0 are already ritz states of (h,s)
!!
!!   MODIFICATION HISTORY
!!     10.2006 made more stable [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE make_full_all(preconditioner_env, matrix_c0, matrix_h, matrix_s, c0_evals, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: fudge_factor = 0.25_dp, &
                                                lambda_base = 10.0_dp

    INTEGER                                  :: k, n
    REAL(KIND=dp)                            :: error_estimate, lambda
    REAL(KIND=dp), DIMENSION(:), POINTER     :: diag, norms, shifted_evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER :: matrix_hc0, matrix_left, matrix_pre, &
      matrix_s1, matrix_s2, matrix_sc0, matrix_shc0, matrix_tmp, ortho

    IF (ASSOCIATED(preconditioner_env%fm)) CALL cp_fm_release(preconditioner_env%fm,error)
    CALL cp_fm_get_info(matrix_c0,nrow_global=n,ncol_global=k,error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp,name="preconditioner_env%fm",error=error)
    matrix_pre=>preconditioner_env%fm
    CALL cp_fm_create(ortho,fm_struct_tmp,name="ortho",error=error)
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp,name="matrix_tmp",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ALLOCATE(preconditioner_env%full_evals(n))
    ALLOCATE(preconditioner_env%occ_evals(k))

    ! 0) cholesky decompose the overlap matrix, if this fails the basis is singular,
    !    more than EPS_DEFAULT
    CALL copy_sm_to_fm(matrix_s,ortho,error=error)
    CALL cp_fm_cholesky_decompose(ortho,error=error)

    ! 1) Construct a new H matrix, which has the current C0 as eigenvectors,
    !    possibly shifted by an amount lambda,
    !    and the same spectrum as the original H matrix in the space orthogonal to the C0
    !    with P=C0 C0 ^ T
    !    (1 - PS)^T H (1-PS) + (PS)^T (H - lambda S ) (PS)
    !    we exploit that the C0 are already the ritz states of H
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0",error=error)
    CALL cp_sm_fm_multiply(matrix_s,matrix_c0,matrix_sc0,k,error=error)
    CALL cp_fm_create(matrix_hc0,matrix_c0%matrix_struct,name="hc0",error=error)
    CALL cp_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

       ! An aside, try to estimate the error on the ritz values, we'll need it later on
       CALL cp_fm_create(matrix_shc0,matrix_c0%matrix_struct,name="shc0",error=error)
       CALL cp_fm_cholesky_restore(matrix_hc0,k,ortho,matrix_shc0,"SOLVE",transa="T",error=error)
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
       CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       ! since we only use diagonal elements this is a bit of a waste
       CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,matrix_shc0,matrix_shc0,0.0_dp,matrix_s1,error=error)
       ALLOCATE(diag(k))
       CALL cp_fm_get_diag(matrix_s1,diag,error=error)
       error_estimate=MAXVAL(SQRT(ABS(diag-c0_evals**2)))
       DEALLOCATE(diag)
       CALL cp_fm_release(matrix_s1,error=error)
       CALL cp_fm_release(matrix_shc0,error=error)
       ! we'll only use the energy gap, if our estimate of the error on the eigenvalues
       ! is small enough. A large error combined with a small energy gap would otherwise lead to 
       ! an aggressive but bad preconditioner. Only when the error is small (MD), we can precondition
       ! aggressively
       preconditioner_env%energy_gap= MAX(energy_gap,error_estimate*fudge_factor)

    CALL copy_sm_to_fm(matrix_h,matrix_tmp,error=error)
    CALL cp_fm_upper_to_full(matrix_tmp,matrix_pre,error=error)
    ! tmp = H ( 1 - PS )
    CALL cp_fm_gemm('N','T',n,n,k,-1.0_dp,matrix_hc0,matrix_sc0,1.0_dp,matrix_tmp,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_left,fm_struct_tmp,name="matrix_left",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_gemm('T','N',k,n,n,1.0_dp,matrix_c0,matrix_tmp,0.0_dp,matrix_left,error=error)
    ! tmp = (1 - PS)^T H (1-PS)
    CALL cp_fm_gemm('N','N',n,n,k,-1.0_dp,matrix_sc0,matrix_left,1.0_dp,matrix_tmp,error=error)
    CALL cp_fm_release(matrix_left,error=error)

    ALLOCATE(shifted_evals(k))
    lambda = lambda_base + error_estimate
    shifted_evals=c0_evals - lambda
    CALL cp_fm_to_fm(matrix_sc0,matrix_hc0,error=error)
    CALL cp_fm_column_scale(matrix_hc0,shifted_evals)
    CALL cp_fm_gemm('N','T',n,n,k,1.0_dp,matrix_hc0,matrix_sc0,1.0_dp,matrix_tmp,error=error)

    ! 2) diagonalize this operator
    CALL cp_fm_cholesky_reduce(matrix_tmp,ortho,error=error)
    CALL cp_fm_syevd(matrix_tmp,matrix_pre,preconditioner_env%full_evals,error=error)
    CALL cp_fm_cholesky_restore(matrix_pre,n,ortho,matrix_tmp,"SOLVE",error=error)
    CALL cp_fm_to_fm(matrix_tmp,matrix_pre,error=error)

    ! test that the subspace remained conserved
    IF (.FALSE.) THEN
        CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
        CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
        CALL cp_fm_create(matrix_s2,fm_struct_tmp,name="matrix_s2",error=error)
        CALL cp_fm_struct_release(fm_struct_tmp,error=error)
        ALLOCATE(norms(k))
        CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,matrix_sc0,matrix_tmp,0.0_dp,matrix_s1,error=error)
        CALL cp_fm_syevd(matrix_s1,matrix_s2,norms,error=error)
        WRITE(6,*) "matrix norm deviation (should be close to zero): ", MAXVAL(ABS(ABS(norms)-1.0_dp))
        DEALLOCATE(norms)
        CALL cp_fm_release(matrix_s1,error=error)
        CALL cp_fm_release(matrix_s2,error=error)
    ENDIF

    ! 3) replace the lowest k evals and evecs with what they should be
    preconditioner_env%occ_evals=c0_evals
    ! notice, this choice causes the preconditioner to be constant when applied to sc0 (see apply_full_all)
    preconditioner_env%full_evals(1:k)=c0_evals 
    CALL cp_fm_to_fm(matrix_c0,matrix_pre,k,1,1)

    CALL cp_fm_release(matrix_sc0,error=error)
    CALL cp_fm_release(matrix_hc0,error=error)
    CALL cp_fm_release(ortho,error=error)
    CALL cp_fm_release(matrix_tmp,error=error)
    DEALLOCATE(shifted_evals)

END SUBROUTINE make_full_all

!
! the corresponding apply_full_all uses the decomposed form to apply the preconditioner
!

SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, j, k, n, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_fm_type), POINTER                :: matrix_tmp

!----

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)

  CALL cp_fm_create(matrix_tmp,matrix_in%matrix_struct,name="apply_full_all",error=error)
  CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data,error=error)

  !
  CALL cp_fm_gemm('T','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_in,0.0_dp,matrix_tmp,error=error)

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
                  matrix_tmp,0.0_dp,matrix_out,error=error)

  CALL cp_fm_release(matrix_tmp,error=error)

END SUBROUTINE apply_full_all

SUBROUTINE make_full_single(preconditioner_env, fm, matrix_h, matrix_s, &
                       energy_homo, energy_gap , error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_s

! ---
! ---

  NULLIFY(fm_h,fm_s,fm_struct_tmp,evals)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm,error=error)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)
  ALLOCATE(evals(n))

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h",error=error)
  CALL cp_fm_create(fm_s,fm_struct_tmp, name="fm_s",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  CALL copy_sm_to_fm(matrix_h,fm_h,error=error)
  CALL copy_sm_to_fm(matrix_s,fm_s,error=error)
  CALL cp_fm_cholesky_decompose(fm_s,error=error)
  CALL cp_fm_cholesky_reduce(fm_h,fm_s,error=error)
  CALL cp_fm_syevd(fm_h,fm,evals,error=error)
  CALL cp_fm_cholesky_restore(fm,n,fm_s,fm_h,"SOLVE",error=error)
  DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
  ENDDO
  CALL cp_fm_to_fm(fm_h,fm,error=error)
  CALL cp_fm_column_scale(fm,evals)
  CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,fm,fm_h,0.0_dp,fm_s,error=error)
  CALL cp_fm_to_fm(fm_s,fm,error=error)

  DEALLOCATE(evals)
  CALL cp_fm_release(fm_h,error=error)
  CALL cp_fm_release(fm_s,error=error)

END SUBROUTINE make_full_single

!***************************************************************************
! different types of preconditioner come here
!***************************************************************************
! only inv(S)
! a full matrix preconditioner
SUBROUTINE apply_full_s_inverse(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: k, n

!----

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)
  CALL cp_fm_symm('L','U',n,k,1.0_dp,preconditioner_env%fm, &
                        matrix_in, &
                        0.0_dp,matrix_out,error=error)

END SUBROUTINE apply_full_s_inverse
SUBROUTINE make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                               energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work

! ---
! ---

  NULLIFY(fm_struct_tmp,fm_work)

  IF (ASSOCIATED(preconditioner_env%fm)) THEN
     CALL cp_fm_release(preconditioner_env%fm,error=error)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  CALL copy_sm_to_fm(matrix_s,preconditioner_env%fm,error=error)
  CALL cp_fm_cholesky_decompose(preconditioner_env%fm,error=error)
  CALL cp_fm_upper_to_full(preconditioner_env%fm,fm_work,error=error) ! symmetrizing makes the apply_ up to 50 percent faster
  CALL cp_fm_release(fm_work,error=error)

END SUBROUTINE make_full_s_inverse



SUBROUTINE make_full_kinetic(preconditioner_env, fm, matrix_t, matrix_s, &
                               energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), POINTER          :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: shift
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work

! ---

  failure = .FALSE.
  CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
! ---

  NULLIFY(fm_struct_tmp,fm_work)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm,error=error)
  ENDIF
  CALL get_matrix_info(matrix_t,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_create(fm_work,fm_struct_tmp, name="preconditioner_work",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  shift=MAX(-energy_homo,energy_gap)
  CALL add_matrices(matrix_t,1.0_dp,matrix_t,shift,matrix_s)

  CALL copy_sm_to_fm(matrix_t,fm,error=error)
  CALL cp_fm_cholesky_decompose(fm,error=error)
  CALL cp_fm_cholesky_invert(fm,error=error)
  CALL cp_fm_upper_to_full(fm,fm_work,error=error)
  CALL cp_fm_release(fm_work,error=error)

END SUBROUTINE make_full_kinetic

END MODULE preconditioner

