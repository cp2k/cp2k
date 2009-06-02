!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_release
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_add_dist2d,&
                                             cp_dbcsr_get_dist2d,&
                                             dbcsr_sm_sm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, &
       cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_release
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_slot_blk_p,&
                                             dbcsr_slot_col_i,&
                                             dbcsr_slot_row_p,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: &
       dbcsr_addto_index_array, dbcsr_checksum, dbcsr_create, &
       dbcsr_create_dist2d_r, dbcsr_create_dist2d_r_unrot, dbcsr_destroy, &
       dbcsr_dist2d_to_dist, dbcsr_ensure_data_size, dbcsr_extract_column, &
       dbcsr_finalize, dbcsr_from_fm, dbcsr_from_sm, dbcsr_get_info, &
       dbcsr_init, dbcsr_new_transposed, dbcsr_set_diagonal, &
       dbcsr_verify_matrix, dbcsr_work_create, dbcsr_work_destroy
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
!USE f77_blas_generic,                ONLY: BL_DOT
  USE input_constants,                 ONLY: &
       ot_precond_full_all, ot_precond_full_kinetic, &
       ot_precond_full_kinetic_sainv, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_s_inverse, &
       ot_precond_solver_default, ot_precond_solver_direct, &
       ot_precond_solver_inv_chol, ot_precond_sparse_diag, &
       ot_precond_sparse_firstorder, ot_precond_sparse_kinetic, &
       ot_precond_sparse_kinetic_sainv
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE preconditioner_types,            ONLY: preconditioner_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE sparse_matrix_types,             ONLY: &
       add_block_node, add_matrices, deallocate_matrix, first_block_node, &
       get_block_node, get_matrix_info, next_block_node, &
       real_block_node_type, real_matrix_type, replicate_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'

  PUBLIC  :: make_preconditioner
  PUBLIC  :: apply_preconditioner

  INTERFACE make_preconditioner
    MODULE PROCEDURE make_preconditioner_single
  END INTERFACE

  PRIVATE :: make_sparse_diag, make_full_single, make_local_block

! *****************************************************************************

CONTAINS

! *****************************************************************************

! creates a preconditioner for the system (H-energy_homo S)
! this preconditioner is (must be) symmetric positive definite.
! currently uses a atom-block-diagonal form
! each block will be  ....
! might overwrite matrix_h, matrix_t

! *****************************************************************************
SUBROUTINE make_preconditioner_single(preconditioner_env, precon_type, solver_type, matrix_h, matrix_s, &
                                      matrix_t,mo_coeff, energy_gap ,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    INTEGER, INTENT(IN)                      :: precon_type, solver_type
    TYPE(real_matrix_type), POINTER          :: matrix_h
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_s, matrix_t
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, k, &
                                                my_solver_type
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot

!----

    failure=.FALSE.

  CALL timeset(routineN,handle)

  CALL cp_fm_get_info(mo_coeff,ncol_global=k,error=error)

  ALLOCATE(eigenvalues_ot(k),STAT=istat)
  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

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
  !
  !
  my_solver_type = solver_type
  SELECT CASE (precon_type)
  CASE (ot_precond_sparse_diag)
     preconditioner_env%in_use=ot_precond_sparse_diag
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program("make preconditioner","Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                              energy_homo, energy_gap ,error=error)
     ELSE
        CALL make_sparse_diag_ortho(preconditioner_env, matrix_h, &
                              energy_homo, energy_gap ,error=error)
     END IF
  CASE (ot_precond_full_single)
     preconditioner_env%in_use=ot_precond_full_single
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program("make preconditioner","Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_single(preconditioner_env, preconditioner_env%fm,&
                              matrix_h, matrix_s, energy_homo, energy_gap ,error=error)
     ELSE
        CALL make_full_single_ortho(preconditioner_env, preconditioner_env%fm,&
                              matrix_h, energy_homo, energy_gap ,error=error)
     END IF
  CASE (ot_precond_s_inverse)
     preconditioner_env%in_use=ot_precond_s_inverse
     !IF(my_solver_type.NE.ot_precond_solver_default) THEN
     !   CALL stop_program("make preconditioner","Only PRECOND_SOLVER DEFAULT for the moment")
     !ENDIF
     IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                                 energy_homo, energy_gap ,error=error)
     ELSE
        CALL stop_program("make preconditioner","Type for S=1 not implemented")
     END IF
  CASE (ot_precond_sparse_firstorder)
     preconditioner_env%in_use=ot_precond_sparse_firstorder
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program("make preconditioner","Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_sparse_first_order(preconditioner_env, matrix_h, matrix_s, &
                                     energy_homo, energy_gap ,error=error)
     ELSE
        CALL make_sparse_first_order_ortho(preconditioner_env, matrix_h, &
                                           energy_homo, energy_gap ,error=error)
     END IF
  CASE (ot_precond_full_kinetic)
     preconditioner_env%in_use=ot_precond_full_kinetic
     IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
     IF ( PRESENT(matrix_s) .AND. PRESENT(matrix_t) ) THEN
        CALL make_full_kinetic(preconditioner_env, preconditioner_env%fm,&
             &                 matrix_t, matrix_s, energy_homo, energy_gap, &
             &                 error=error)
     ELSE
        CALL stop_program("make preconditioner","Type for S=1 not implemented")
     ENDIF
  CASE (ot_precond_full_single_inverse)
     preconditioner_env%in_use=ot_precond_full_single_inverse
     IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_single_inverse(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                                      eigenvalues_ot, energy_gap,error=error)
     ELSE
        CALL make_full_single_inverse_ortho(preconditioner_env,mo_coeff,matrix_h, &
                                            eigenvalues_ot, energy_gap,error=error)
     END IF
  CASE (ot_precond_full_all)
     preconditioner_env%in_use=ot_precond_full_all
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program("make preconditioner","Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_all(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                           eigenvalues_ot, energy_gap,error=error)
     ELSE
        CALL make_full_all_ortho(preconditioner_env,mo_coeff,matrix_h, &
                                 eigenvalues_ot, energy_gap,error=error)
     END IF
  CASE(ot_precond_sparse_kinetic)
     preconditioner_env%in_use=ot_precond_sparse_kinetic
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program("make preconditioner","Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     CALL make_sparse_kinetic(preconditioner_env, preconditioner_env%sparse_matrix,&
                              matrix_t, matrix_s, energy_homo, energy_gap ,error=error)
     CALL make_diag_inner_precond(preconditioner_env, preconditioner_env%sparse_matrix,&
                                  preconditioner_env%sparse_matrix_inner, error=error)
  CASE (ot_precond_full_kinetic_sainv)
    preconditioner_env%in_use=ot_precond_full_kinetic_sainv
    IF (my_solver_type .NE. ot_precond_solver_inv_chol&
         .AND. my_solver_type .NE. ot_precond_solver_default) THEN
       CALL stop_program ("make preconditioner",&
            "Only PRECOND_SOLVER INVERSE_CHOLESKY (or DEFAULT) for the moment")
    ENDIF
    CALL make_full_kinetic_sainv(preconditioner_env, preconditioner_env%fm,&
         matrix_t, matrix_s, energy_homo, energy_gap, error=error)
  CASE (ot_precond_sparse_kinetic_sainv)
    preconditioner_env%in_use=ot_precond_sparse_kinetic_sainv
    IF (my_solver_type .NE. ot_precond_solver_inv_chol&
         .AND. my_solver_type .NE. ot_precond_solver_default) THEN
       CALL stop_program ("make preconditioner",&
            "Only PRECOND_SOLVER INVERSE_CHOLESKY (or DEFAULT) for the moment")
    ENDIF
    CALL make_sparse_kinetic_sainv(preconditioner_env,&
         preconditioner_env%dbcsr_matrix, matrix_t, matrix_s,&
         energy_homo, energy_gap ,error=error)
    CALL make_full_kinetic_sainv(preconditioner_env, preconditioner_env%fm,&
         matrix_t, matrix_s, energy_homo, energy_gap, error=error)
  CASE DEFAULT
     CALL stop_program("make preconditioner","Type not implemented")
  END SELECT
  !
  ! here comes the solver
  SELECT CASE(my_solver_type)
  CASE (ot_precond_solver_inv_chol)
     !
     ! compute the full inverse
     preconditioner_env%solver=ot_precond_solver_inv_chol
     CALL make_full_inverse_cholesky(preconditioner_env, preconditioner_env%fm, matrix_s, &
          &                          error=error)
  CASE (ot_precond_solver_direct)
     !
     ! prepare for the direct solver
     preconditioner_env%solver=ot_precond_solver_direct
     CALL make_full_fact_cholesky(preconditioner_env, preconditioner_env%fm, matrix_s, &
          &                       energy_homo, energy_gap, error)
   CASE (ot_precond_solver_default)
     preconditioner_env%solver=ot_precond_solver_default
  CASE DEFAULT
     !
     CALL stop_program("make preconditioner","doesnt know this type of solver")
  END SELECT

  DEALLOCATE(eigenvalues_ot,STAT=istat)
  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  CALL timestop(handle)

END SUBROUTINE make_preconditioner_single


! applies a previously created preconditioner to a full matrix
! *****************************************************************************
SUBROUTINE apply_preconditioner(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

!----

  CALL timeset(routineN,handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
     CALL stop_program("apply preconditioner","No preconditioner in use")
  CASE (ot_precond_sparse_diag,ot_precond_sparse_firstorder)
     CALL apply_sparse_diag(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_single)
     CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_all)
     CALL apply_full_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE(ot_precond_sparse_kinetic)
     CALL apply_solve_lin_system(preconditioner_env, matrix_in, matrix_out, error=error)
  CASE (ot_precond_full_kinetic_sainv)
     SELECT CASE (preconditioner_env%solver)
     CASE (ot_precond_solver_inv_chol, ot_precond_solver_default)
        CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE DEFAULT
        CALL stop_program("apply preconditioner","solver not implemented")
     END SELECT
  CASE (ot_precond_sparse_kinetic_sainv)
     SELECT CASE (preconditioner_env%solver)
     CASE (ot_precond_solver_inv_chol, ot_precond_solver_default)
        CALL apply_sparse_single(preconditioner_env, matrix_in, matrix_out,error=error)
        CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE DEFAULT
        CALL stop_program("apply preconditioner","solver not implemented")
     END SELECT
  CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
     SELECT CASE (preconditioner_env%solver)
     CASE(ot_precond_solver_inv_chol)
        CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE(ot_precond_solver_direct)
        CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out,error=error)        
     CASE DEFAULT
        CALL stop_program("apply preconditioner","solver not implemented")
     END SELECT
  CASE DEFAULT
     CALL stop_program("apply preconditioner","implemented")
  END SELECT

  CALL timestop(handle)

END SUBROUTINE apply_preconditioner

! different types of preconditioner come here
! a sparse block diagonal approximation
! *****************************************************************************
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

! *****************************************************************************
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

  CALL Replicate_matrix(matrix_h,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.,error=error)
  CALL get_matrix_info(matrix_h,nblock_row=nblocks)
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
! *****************************************************************************
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
! different types of preconditioner come here
! a sparse one based on 2inv(A)-inv(A) F inv(A) where F is S
! turns out that this inverse is not symmetric positive definite unless
! inv(A) is scaled by some number
! *****************************************************************************
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

  CALL get_matrix_info(matrix_h,nblock_row=nblocks,first_row=first_row,last_row=last_row)

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
  CALL Replicate_matrix(matrix_h,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.,error=error)

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

! different types of preconditioner come here
! a full matrix preconditioner
! *****************************************************************************
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

SUBROUTINE apply_full_direct(preconditioner_env, matrix_in, matrix_out,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: k, n
    TYPE(cp_fm_type), POINTER                :: work

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)
  CALL cp_fm_create(work,matrix_in%matrix_struct,name="apply_full_single",error=error)
  CALL cp_fm_cholesky_restore(matrix_in,k,preconditioner_env%fm,work,&
       &                      "SOLVE",transa="T",error=error)
  CALL cp_fm_cholesky_restore(work,k,preconditioner_env%fm,matrix_out,&
       &                      "SOLVE",transa="N",error=error)
  CALL cp_fm_release(work,error=error)
    
END SUBROUTINE apply_full_direct

! *****************************************************************************
!> \brief generates a preconditioner by cholesky inverting H-lambda S+(SC)shifts(SC)^T
!> \param matrix_c 0 must be already rotated correctly.
!> \note
!>      this might fail if the initial guess is bad, or if the system has 'holes'
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      10.2006 made more robust [Joost VandeVondele]
! *****************************************************************************
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
    REAL(KIND=dp)                            :: error_estimate, &
                                                preconditioner_shift
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_shift, diag, evals
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

    CALL cp_fm_release(matrix_tmp,error=error)

END SUBROUTINE make_full_single_inverse

! *****************************************************************************
!> \brief generates a state by state preconditioner based on the full hamiltonian matrix 
!> \param energy_gap should be a slight underestimate of the physical energy gap for almost all systems
!>      the c0 are already ritz states of (h,s)
!> \note
!>      includes error estimate on the hamiltonian matrix to result in a stable preconditioner
!>      a preconditioner for each eigenstate i is generated by keeping the factorized form
!>      U diag( something i ) U^T. It is important to only precondition in the subspace orthogonal to c0.
!>      not only is it the only part that matters, it also simplifies the computation of
!>      the lagrangian multipliers in the OT minimization  (i.e. if the c0 here is different
!>      from the c0 used in the OT setup, there will be a bug).
!> \par History
!>      10.2006 made more stable [Joost VandeVondele]
! *****************************************************************************
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

! *****************************************************************************
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

! *****************************************************************************
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

! different types of preconditioner come here
! *****************************************************************************
SUBROUTINE make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, &
                               energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp

! ---

  NULLIFY(fm_struct_tmp)

  IF (ASSOCIATED(preconditioner_env%fm)) THEN
     CALL cp_fm_release(preconditioner_env%fm,error=error)
  ENDIF
  CALL get_matrix_info(matrix_h,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  CALL copy_sm_to_fm(matrix_s,preconditioner_env%fm,error=error)

END SUBROUTINE make_full_s_inverse

! *****************************************************************************
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

! ---

  failure = .FALSE.
  CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
! ---

  NULLIFY(fm_struct_tmp)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm,error=error)
  ENDIF
  CALL get_matrix_info(matrix_t,nrow=n)

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  shift=MAX(-energy_homo,energy_gap)
  CALL add_matrices(matrix_t,1.0_dp,matrix_t,shift,matrix_s)
  CALL copy_sm_to_fm(matrix_t,fm,error=error)

END SUBROUTINE make_full_kinetic

  SUBROUTINE make_full_inverse_cholesky(preconditioner_env, fm, matrix_s, error)
    
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_inverse_cholesky', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_fm_type), POINTER                :: fm_work

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)

    NULLIFY(fm_work)

    CALL cp_fm_create(fm_work,fm%matrix_struct,name="fm_work",error=error)
    !
    ! compute the inverse of SPD matrix fm using the Cholesky factorization
    CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
    CALL cp_fm_cholesky_decompose(fm,error=sub_error)
    failure = .FALSE.
    CALL cp_error_check(sub_error,failure)
    CALL cp_error_dealloc_ref(sub_error)
    !
    ! if fm not SPD we go with the overlap matrix
    IF (failure) THEN
       !
       ! just the overlap matrix
       IF(PRESENT(matrix_s)) THEN
          CALL copy_sm_to_fm(matrix_s,fm,error=error)
          CALL cp_fm_cholesky_decompose(fm,error=error)
       ELSE
          CALL cp_fm_set_all(fm,alpha=0._dp,beta=1._dp,error=error)
       ENDIF
    ENDIF
    CALL cp_fm_cholesky_invert(fm,error=error)
    CALL cp_fm_upper_to_full(fm,fm_work,error=error)
    CALL cp_fm_release(fm_work,error=error)

  END SUBROUTINE make_full_inverse_cholesky

  SUBROUTINE make_full_fact_cholesky(preconditioner_env, fm, matrix_s, &
                                     energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_fact_cholesky', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)
    !
    ! compute the inverse of SPD matrix fm using the Cholesky factorization
    CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
    CALL cp_fm_cholesky_decompose(fm,error=sub_error)
    failure = .FALSE.
    CALL cp_error_check(sub_error,failure)
    CALL cp_error_dealloc_ref(sub_error)
    !
    ! if fm not SPD we go with the overlap matrix
    IF (failure) THEN
       !
       ! just the overlap matrix
       IF(PRESENT(matrix_s)) THEN
          CALL copy_sm_to_fm(matrix_s,fm,error=error)
          CALL cp_fm_cholesky_decompose(fm,error=error)
       ELSE
          CALL cp_fm_set_all(fm,alpha=0._dp,beta=1._dp,error=error)
       ENDIF
    ENDIF

  END SUBROUTINE make_full_fact_cholesky

! *****************************************************************************
! Preconditioners for orthogonal basis sets
! *****************************************************************************
SUBROUTINE make_sparse_diag_ortho(preconditioner_env, matrix_h, energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, iblock_col_h, iblock_row, &
                                                n, nblocks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s
    TYPE(real_block_node_type), POINTER      :: block_node_h

!----

  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
     CALL deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
     NULLIFY(preconditioner_env%sparse_matrix)
  ENDIF

  CALL Replicate_matrix(matrix_h,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.,error=error)
  CALL get_matrix_info(matrix_h,nblock_row=nblocks)
  DO iblock_row=1,nblocks
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_h))
        IF (.NOT.ASSOCIATED(block_node_h)) CALL stop_program("qs preconditioner", ".NOT.ASSOCIATED")

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            BLOCK=block_h)

        IF (iblock_col_h .EQ. iblock_row) THEN
           n=SIZE(block_h,1)
           ALLOCATE(block_pre(n,n),block_s(n,n))
           block_s(:,:)=0._dp
           DO i=1,n
             block_s(i,i)=1._dp
           END DO

           CALL make_local_block(block_h,block_s,block_pre,energy_homo,energy_gap,3)

           CALL add_block_node(matrix=preconditioner_env%sparse_matrix,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col_h,&
                                 BLOCK=block_pre,error=error)

           DEALLOCATE(block_pre,block_s)

        ENDIF
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO
END SUBROUTINE make_sparse_diag_ortho
! *****************************************************************************
! different types of preconditioner come here
! a sparse one based on 2inv(A)-inv(A) F inv(A) where F is S
! turns out that this inverse is not symmetric positive definite unless
! inv(A) is scaled by some number
! *****************************************************************************
SUBROUTINE make_sparse_first_order_ortho(preconditioner_env, matrix_h, &
                                   energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: matrix_h
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, iblock_col_h, &
                                                iblock_col_s, iblock_row, k, &
                                                n, nblocks
    INTEGER, DIMENSION(:), POINTER           :: first_row, last_row, offset
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block_buffer
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s, &
                                                block_tmp
    TYPE(real_block_node_type), POINTER      :: block_node_h

!----

  IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
     CALL deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
     NULLIFY(preconditioner_env%sparse_matrix)
  ENDIF
  ! first we get the all blocks of the block diagonal preconditioner replicated
  ! on every CPU (store in block_buffer)

  CALL get_matrix_info(matrix_h,nblock_row=nblocks,first_row=first_row,last_row=last_row)

  ALLOCATE(offset(nblocks+1))
  offset(1)=0
  DO iblock_row=1,nblocks
     offset(iblock_row+1)=offset(iblock_row)+(last_row(iblock_row)-first_row(iblock_row)+1)**2
  ENDDO
  ALLOCATE(block_buffer(offset(nblocks+1)))
  CALL DCOPY(offset(nblocks+1),0.0_dp,0,block_buffer(1),1)

  DO iblock_row=1,nblocks
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_h))

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            BLOCK=block_h)

        IF (iblock_col_h .EQ. iblock_row) THEN
           n=SIZE(block_h,1)
           ALLOCATE(block_pre(n,n),block_s(n,n))
           block_s = 0._dp
           DO i=1,n
             block_s(i,i) = 1._dp
           END DO

           CALL make_local_block(block_h,block_s,block_pre,  &
                                                  energy_homo,energy_gap,4)
           CALL DCOPY(n*n,block_pre(1,1),1,block_buffer(offset(iblock_row)+1),1)

           DEALLOCATE(block_pre,block_s)

        ENDIF
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO

  CALL mp_sum(block_buffer,preconditioner_env%para_env%group)

  CALL DSCAL(offset(nblocks+1),0.2_dp,block_buffer(1),1)

  ! now every diagonal block is everywhere, form the product
  CALL Replicate_matrix(matrix_h,preconditioner_env%sparse_matrix,"PRECONDITIONER",.FALSE.,error=error)

  DO iblock_row=1,nblocks
     block_node_h => first_block_node(matrix_h,iblock_row)
     DO WHILE (ASSOCIATED(block_node_h))

        CALL get_block_node(block_node=block_node_h,&
                            block_col=iblock_col_h,&
                            BLOCK=block_h)

        n=SIZE(block_h,1)
        k=SIZE(block_h,2)
        ALLOCATE(block_pre(n,k),block_tmp(n,k))
        block_pre(:,:)=0._dp
        IF (iblock_col_h .EQ. iblock_row) THEN
          DO i=1,n
            block_pre(i,i)=-1._dp
          END DO
        END IF
        CALL DGEMM('N','N',n,k,k,1.0_dp,block_pre(1,1),n, &
                                        block_buffer(offset(iblock_col_s)+1), &
                                        k,0.0_dp,block_tmp(1,1),n)
        CALL DGEMM('N','N',n,k,n,1.0_dp,block_buffer(offset(iblock_row)+1),n, &
                                        block_tmp(1,1), &
                                        n,0.0_dp,block_pre(1,1),n)
        IF (iblock_row .EQ. iblock_col_s) &
           CALL DAXPY(n*k,2.0_dp,block_buffer(offset(iblock_row)+1),1,block_pre(1,1),1)

        CALL add_block_node(preconditioner_env%sparse_matrix,iblock_row,iblock_col_s,block_pre,error=error)

        DEALLOCATE(block_pre,block_tmp)
        block_node_h => next_block_node(block_node_h)
     END DO
  ENDDO

  DEALLOCATE(offset,block_buffer)

END SUBROUTINE make_sparse_first_order_ortho
! *****************************************************************************
SUBROUTINE make_full_single_ortho(preconditioner_env, fm, matrix_h, &
                       energy_homo, energy_gap , error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), POINTER          :: matrix_h
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_s

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
  CALL cp_fm_syevd(fm_h,fm,evals,error=error)
  DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
  ENDDO
  CALL cp_fm_to_fm(fm,fm_h,error=error)
  CALL cp_fm_column_scale(fm,evals)
  CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,fm,fm_h,0.0_dp,fm_s,error=error)
  CALL cp_fm_to_fm(fm_s,fm,error=error)

  DEALLOCATE(evals)
  CALL cp_fm_release(fm_h,error=error)
  CALL cp_fm_release(fm_s,error=error)

END SUBROUTINE make_full_single_ortho
! *****************************************************************************
SUBROUTINE make_full_single_inverse_ortho(preconditioner_env, matrix_c0, matrix_h, c0_evals, energy_gap,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(real_matrix_type), POINTER          :: matrix_h
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: eval_shift = 5.0_dp , &
                                                fudge_factor = 2.0_dp 

    INTEGER                                  :: k, n
    REAL(KIND=dp)                            :: error_estimate, &
                                                preconditioner_shift
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_shift, diag, evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_hc0, matrix_s1, &
                                                matrix_sc0, matrix_tmp, &
                                                matrix_tmp2

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

    ! get the error estimate
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0",error=error)
    CALL cp_fm_to_fm(matrix_c0,matrix_sc0,error=error)
    CALL cp_fm_create(matrix_hc0,matrix_c0%matrix_struct,name="hc0",error=error)
    CALL cp_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ! since we only use diagonal elements this is a bit of a waste
    CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,matrix_hc0,matrix_hc0,0.0_dp,matrix_s1,error=error)
    ALLOCATE(diag(k))
    CALL cp_fm_get_diag(matrix_s1,diag,error=error)
    error_estimate=MAXVAL(SQRT(ABS(diag-c0_evals**2)))
    DEALLOCATE(diag)
    CALL cp_fm_release(matrix_s1,error=error)
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

    ! preconditioner shift, we target the middle of the occupied spectrum, and taking into account the error_estimate
    ! write(6,*) "Error estimate = ",error_estimate
    preconditioner_shift=-(MINVAL(c0_evals)+ MAXVAL(c0_evals))/2.0_dp + &
                           error_estimate*fudge_factor
    CALL cp_fm_set_all(matrix_tmp,alpha=0._dp,beta=1._dp,error=error)
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

    CALL cp_fm_release(matrix_tmp,error=error)

END SUBROUTINE make_full_single_inverse_ortho
! *****************************************************************************
SUBROUTINE make_full_all_ortho(preconditioner_env, matrix_c0, matrix_h, c0_evals, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(real_matrix_type), POINTER          :: matrix_h
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: fudge_factor = 0.25_dp, &
                                                lambda_base = 10.0_dp

    INTEGER                                  :: k, n
    REAL(KIND=dp)                            :: error_estimate, lambda
    REAL(KIND=dp), DIMENSION(:), POINTER     :: diag, norms, shifted_evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_hc0, matrix_left, &
                                                matrix_pre, matrix_s1, &
                                                matrix_s2, matrix_sc0, &
                                                matrix_tmp

    IF (ASSOCIATED(preconditioner_env%fm)) CALL cp_fm_release(preconditioner_env%fm,error)
    CALL cp_fm_get_info(matrix_c0,nrow_global=n,ncol_global=k,error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp,name="preconditioner_env%fm",error=error)
    matrix_pre=>preconditioner_env%fm
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp,name="matrix_tmp",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ALLOCATE(preconditioner_env%full_evals(n))
    ALLOCATE(preconditioner_env%occ_evals(k))

    ! 1) Construct a new H matrix, which has the current C0 as eigenvectors,
    !    possibly shifted by an amount lambda,
    !    and the same spectrum as the original H matrix in the space orthogonal to the C0
    !    with P=C0 C0 ^ T
    !    (1 - PS)^T H (1-PS) + (PS)^T (H - lambda S ) (PS)
    !    we exploit that the C0 are already the ritz states of H
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0",error=error)
    CALL cp_fm_to_fm(matrix_c0,matrix_sc0,error=error)
    CALL cp_fm_create(matrix_hc0,matrix_c0%matrix_struct,name="hc0",error=error)
    CALL cp_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

       ! An aside, try to estimate the error on the ritz values, we'll need it later on
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
       CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       ! since we only use diagonal elements this is a bit of a waste
       CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,matrix_hc0,matrix_hc0,0.0_dp,matrix_s1,error=error)
       ALLOCATE(diag(k))
       CALL cp_fm_get_diag(matrix_s1,diag,error=error)
       error_estimate=MAXVAL(SQRT(ABS(diag-c0_evals**2)))
       DEALLOCATE(diag)
       CALL cp_fm_release(matrix_s1,error=error)
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
    CALL cp_fm_syevd(matrix_tmp,matrix_pre,preconditioner_env%full_evals,error=error)

    ! test that the subspace remained conserved
    IF (.FALSE.) THEN
        CALL cp_fm_to_fm(matrix_pre,matrix_tmp,error=error)
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
    CALL cp_fm_release(matrix_tmp,error=error)
    DEALLOCATE(shifted_evals)

END SUBROUTINE make_full_all_ortho
! *****************************************************************************

  SUBROUTINE make_diag_inner_precond(preconditioner_env, sm, sm_inner, error)
    !
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: sm, sm_inner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_diag_inner_precond', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iblock_col_h, iblock_row, &
                                                info, istat, liwork, lwork, &
                                                n, nblocks
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: block_evec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_buf1, block_h, block_pre
    TYPE(real_block_node_type), POINTER      :: block_node_h

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(sm),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(sm_inner)) THEN
       CALL deallocate_matrix(sm_inner,error=error)
       NULLIFY(sm_inner)
    ENDIF

    CALL Replicate_matrix(sm,sm_inner,"PRECONDITIONER_INNER",.FALSE.,error=error)
    CALL get_matrix_info(sm,nblock_row=nblocks)
    DO iblock_row=1,nblocks
       block_node_h => first_block_node(sm,iblock_row)
       DO WHILE (ASSOCIATED(block_node_h))
          IF (.NOT.ASSOCIATED(block_node_h)) CALL stop_program("qs preconditioner", ".NOT.ASSOCIATED")

          CALL get_block_node(block_node=block_node_h,&
                              block_col=iblock_col_h,&
                              BLOCK=block_h)

          IF (iblock_col_h .EQ. iblock_row) THEN
             n = SIZE(block_h,1)
             ALLOCATE(block_pre(n,n),block_evec(n,n),block_buf1(n,n),evals(n))
             lwork=1+6*n+2*n**2+50
             liwork=5*n+3
             ALLOCATE(work(lwork),iwork(liwork),STAT=istat)
             IF (istat.NE.0) CALL stop_memory("use preconditioner","work")

             !
             ! A simple block diagonal preconditoner
             block_evec(:,:) = block_h(:,:)
             CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                         iwork(1), liwork, info)
             IF (info.ne.0) CALL stop_program("preconditioner","DSYEVD H")
             block_pre(:,:)=0.0_dp
             DO i=1,n
                IF(evals(i).lt.1e-10_dp) CALL stop_program("preconditioner","evals(i).lt.1e-10_dp")
                block_pre(i,i)=1.0_dp/evals(i)
             ENDDO
             CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                        0.0_dp,block_buf1(1,1),n)
             CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                        0.0_dp,block_pre(1,1),n)
             !
             ! add the block
             CALL add_block_node(matrix=sm_inner,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col_h,&
                                 BLOCK=block_pre,error=error)

             DEALLOCATE(block_pre,block_evec,block_buf1,evals,work,iwork)

          ENDIF
          block_node_h => next_block_node(block_node_h)
       ENDDO
    ENDDO
  END SUBROUTINE make_diag_inner_precond

  SUBROUTINE apply_solve_lin_system(preconditioner_env, matrix_in, matrix_out, error)
    !
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_solve_lin_system', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_iter = 200, &
                                                restart_every = 50
    REAL(dp), PARAMETER                      :: eps = 1.0e-2_dp

    INTEGER                                  :: iter, m, n, output_unit
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, beta, residual_ot, &
                                                tr_pAp, tr_rr, tr_sr_new, &
                                                tr_sr_old
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: Ap, p, r, s
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(preconditioner_env%sparse_matrix),cp_failure_level,routineP,error,failure)
    !
    CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=m,error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=m,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(r ,fm_struct_tmp,name="solve_r" ,error=error)
    CALL cp_fm_create(p ,fm_struct_tmp,name="solve_p" ,error=error)
    CALL cp_fm_create(s ,fm_struct_tmp,name="solve_s" ,error=error)
    CALL cp_fm_create(Ap,fm_struct_tmp,name="solve_Ap",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(p ,0.0_dp,error=error)
    !
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_unit_nr(logger)
    !
    ! residual
    CALL cp_fm_trace(matrix_in,matrix_in,residual_ot,error=error)
    residual_ot = SQRT(residual_ot/(REAL(n,dp)*REAL(m,dp)))
    !
    ! r = b-A*x
    CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_out,r,m,error=error)
    CALL cp_fm_scale_and_add(-1.0_dp,r,1.0_dp,matrix_in,error=error)
    CALL cp_fm_trace(r,r,tr_rr,error=error)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',0,': ||r_ls||=',&
            & SQRT(tr_rr/(REAL(n,dp)*REAL(m,dp))),' ||r_ot||=',residual_ot
    ENDIF
    !
    ! let's go!
    DO iter = 1,max_iter
       !
       ! s = M * r
       IF(ASSOCIATED(preconditioner_env%sparse_matrix_inner)) THEN
          CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix_inner,r,s,m,error=error)
       ELSE
          CALL cp_fm_to_fm(r,s,error=error)
       ENDIF
       CALL cp_fm_trace(s,r,tr_sr_new,error=error)
       !
       ! beta (this might be better when O(N) beta=(s_k,r_k-r_{k-1})/(s_{k-1},r_{k-1}))
       IF(iter.EQ.1) THEN
          beta = 0.0_dp
       ELSE
          beta = tr_sr_new/tr_sr_old
       ENDIF
       !
       ! p = r + beta * p
       CALL cp_fm_scale_and_add(beta,p,1.0_dp,s,error=error)
       !
       ! Ap = A * p
       CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,p,Ap,m,error=error)
       CALL cp_fm_trace(Ap,p,tr_pAp,error=error)
       !
       ! alpha = Tr(s_k'*r_k)/Tr(p_k'*Ap_k)
       alpha = tr_sr_new/tr_pAp
       !
       ! x = x + alpha * p
       CALL cp_fm_scale_and_add(1.0_dp,matrix_out,alpha,p,error=error)
       !
       ! r = r - alpha * Ap or r = b - A * x
       IF(MOD(iter,restart_every).EQ.0) THEN
          CALL cp_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_out,r,m,error=error)
          CALL cp_fm_scale_and_add(-1.0_dp,r,1.0_dp,matrix_in,error=error)
       ELSE
          CALL cp_fm_scale_and_add(1.0_dp,r,-alpha,Ap,error=error)
       ENDIF
       tr_sr_old = tr_sr_new
       !
       ! printing
       IF(MOD(iter,10).EQ.0.AND.output_unit>0) THEN
          WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',iter,': ||r_ls||=',&
               & SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp))),' ||r_ls||/||r_ot||=',&
               & SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot
       ENDIF
       !
       ! exit
       IF(SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot.LT.eps) EXIT
       !
    ENDDO
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',iter,': ||r_ls||=',&
            & SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp))),' ||r_ls||/||r_ot||=',&
            & SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot
    ENDIF
    IF(iter.GT.max_iter) THEN
       CALL stop_program("apply_solve_lin_system","increase max_iter")
    ENDIF
    !
    CALL cp_fm_release(r ,error=error)
    CALL cp_fm_release(p ,error=error)
    CALL cp_fm_release(s ,error=error)
    CALL cp_fm_release(Ap,error=error)
  END SUBROUTINE apply_solve_lin_system

! *****************************************************************************
!> \brief Kinetic Stable Approximate Inverse Preconditioner
!> \param[in,out] preconditioner_env    The preconditioner environment
!> \param         fm                    The precondition matrix
!> \param         matrix_t              The MO matrix
!> \param         matrix_s              The overlap matrix
!> \param[in]     energy_homo           HOMO energy
!> \param[in]     energy_gap            energy gap
!> \param[in,out] error       cp2k error
!> \author Valery Weber
!> \par History
!>      - Created by VW
!>      - 2009-05-13 UB: Merged and brought up-to-date
! *****************************************************************************
  SUBROUTINE make_full_kinetic_sainv(preconditioner_env,&
       fm, matrix_t, matrix_s, energy_homo, energy_gap, error)
    !
    TYPE(preconditioner_type), INTENT(INOUT) :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(real_matrix_type), POINTER          :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_kinetic_sainv', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: thresh_inner = 1.0E-2_dp, &
                                                thresh_post = 1.0E-5_dp

    INTEGER                                  :: i, j, k, n, ncol_local, &
                                                nrow_local, nz
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: DDOT, mx, shift
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p, v
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: buf1, buf2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work1, fm_work2

!
!
! These threshould should be make specifiable in the input. -UB
!

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    !
    WRITE(*,*)'enter sainv thresh_inner=',thresh_inner,' thresh_post=',thresh_post
    IF(ASSOCIATED(fm)) THEN
       CALL cp_fm_release(fm,error=error)
    ENDIF
    !
    ! some infos
    CALL get_matrix_info(matrix_t,nrow=n)
    !
    ! allocate
    ALLOCATE(buf1(n,1),buf2(n,1),p(n),v(n))
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
         context=preconditioner_env%ctxt, &
         para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner",error=error)
    !
    multiproc: IF (preconditioner_env%para_env%num_pe .EQ. 1) THEN
    !
    CALL cp_fm_create(fm_work1,fm_struct_tmp, name="preconditioner_work1",error=error)
    CALL cp_fm_create(fm_work2,fm_struct_tmp, name="preconditioner_work2",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    !
    ! T = T - epsilon * S
    shift=MAX(-energy_homo,energy_gap)
    CALL add_matrices(matrix_t,1.0_dp,matrix_t,shift,matrix_s)
    CALL copy_sm_to_fm(matrix_t,fm_work1,error=error)
    CALL cp_fm_upper_to_full(fm_work1,fm_work2,error=error)
    !
    ! occupation
    CALL cp_fm_get_info(fm_work1,nrow_local=nrow_local,ncol_local=ncol_local, &
         local_data=local_data,error=error)
    nz=0
    DO j=1,ncol_local
       DO i=1,nrow_local
          IF(ABS(local_data(i,j)).LT.1.0E-10) THEN
             local_data(i,j)=0.0_dp
             nz=nz+1
          ENDIF
       ENDDO
    ENDDO
    WRITE(*,*) 'non0(T-epsilon*S)=',DBLE(n**2-nz)/DBLE(n**2)
    !
    ! max(abs(A))
    CALL cp_fm_get_info(fm_work1,nrow_local=nrow_local,ncol_local=ncol_local, &
         local_data=local_data,error=error)
    mx = MAXVAL(ABS(local_data(1:nrow_local,1:ncol_local)))
    CALL mp_max(mx,preconditioner_env%para_env%group)
    WRITE(*,*) 'max(abs(A))=',mx
    !
    ! cond(A)
    CALL cp_fm_to_fm(fm_work1,fm,error=error)
    CALL cp_fm_syevd(fm,fm_work2,v,error=error)
    WRITE(*,*) 'cond(A)=',MAXVAL(ABS(v))/MAX(MINVAL(ABS(v)),EPSILON(0.0_dp))
    WRITE(*,*) 'min(eig(A))=',MINVAL(v),'max(eig(A))=',MAXVAL(v)
    !
    ! Z = diag(ones(n,1));
    CALL cp_fm_set_all(fm_work2,0.0_dp,1.0_dp,error=error)
    !
    ! Let's ainv
    DO i = 1,n
       !
       ! v=A*Z(:,i);
       CALL cp_fm_get_submatrix(fm_work2,buf2 ,&
            &                   1,i,n,1,transpose=.FALSE.,&
            &                   error=error)     
       CALL dgemv('N',n,n,1.0_dp,fm_work1%local_data(1,1),n,buf2(1,1),1,0.0_dp,v(1),1)
       !
       DO j = i,n
          !
          ! p(j) = v(:)'*Z(:,j);
          CALL cp_fm_get_submatrix(fm_work2,buf2 ,&
               &                   1,j,n,1,transpose=.FALSE.,&
               &                   error=error)     
          !p(j) = BL_DOT(n,v(1:n),1,buf2(1:n,1),1)
          p(j) = DDOT(n,v(1:n),1,buf2(1:n,1),1)
       ENDDO
       !
       IF(i.LT.n) THEN
          CALL cp_fm_get_submatrix(fm_work2,buf1 ,&
               &                   1,i,n,1,transpose=.FALSE.,&
               &                   error=error)
          DO j = i+1,n
             !
             ! Z(:,j) = Z(:,j)-(p(j)/p(i))*Z(:,i);
             CALL cp_fm_get_submatrix(fm_work2,buf2 ,&
                  &                   1,j,n,1,transpose=.FALSE.,&
                  &                   error=error)
             IF(p(i).LT.0.0_dp) THEN
                WRITE(*,*) 'make_kinetic_sainv: i=',i,' p(i)=',p(i)
                CALL stop_program("make_kinetic_sainv", "Ooops pivot negative")
             ENDIF
             IF(p(i).LT.1.0E-12_dp.AND.p(i).GE.0.0_dp) THEN
                WRITE(*,*) 'make_kinetic_sainv: i=',i,' p(i)=',p(i)
                CALL stop_program("make_kinetic_sainv", "Ooops pivot too small")
             ENDIF
             buf2 = buf2 - p(j) / p(i) * buf1
             !
             ! filter...           
             DO k=1,j-1
                IF(ABS(buf2(k,1)).LT.mx*thresh_inner) THEN
                   buf2(k,1)=0.0_dp
                ENDIF
             ENDDO
             CALL cp_fm_set_submatrix(fm_work2,buf2 ,&
                  &                   1,j,n,1,transpose=.FALSE.,&
                  &                   error=error)
          ENDDO
       ENDIF
    ENDDO
    !
    ! M = Z * D^(-1) * Z'
    DO i = 1,n
       p(i) = 1.0_dp / p(i)
    ENDDO
    CALL cp_fm_to_fm(fm_work2,fm_work1,error=error)
    CALL cp_fm_column_scale(fm_work1,p)
    CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,fm_work2,fm_work1,0.0_dp,fm,error=error)
    !
    ! post-filter...
    CALL cp_fm_get_info(fm,nrow_local=nrow_local,ncol_local=ncol_local, &
         local_data=local_data,error=error)
    nz=0
    DO j=1,ncol_local
       DO i=1,nrow_local
          IF(i.NE.j) THEN
             IF(ABS(local_data(i,j)).LT.mx*thresh_post) THEN
                local_data(i,j)=0.0_dp
                nz=nz+1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CALL mp_sum(nz,preconditioner_env%para_env%group)
    WRITE(*,*) 'non0(M)=',DBLE(n**2-nz)/DBLE(n**2)
    mx = MAXVAL(ABS(local_data(1:nrow_local,1:ncol_local)))
    CALL mp_max(mx,preconditioner_env%para_env%group)
    WRITE(*,*) 'max(abs(M))=',mx
    !
    ! eig(M)
    CALL cp_fm_to_fm(fm,fm_work1,error=error)
    CALL cp_fm_syevd(fm_work1,fm_work2,v,error=error)
    WRITE(*,*) 'cond(M)=',MAXVAL(ABS(v))/MAX(MINVAL(ABS(v)),EPSILON(0.0_dp))
    WRITE(*,*) 'min(eig(M))=',MINVAL(v),'max(eig(M))=',MAXVAL(v)
    !
    !
    DEALLOCATE(buf1,buf2,p,v)
    CALL cp_fm_release(fm_work1,error=error)
    CALL cp_fm_release(fm_work2,error=error)
    ELSE ! multiproc
    CALL cp_fm_set_all(fm,0.0_dp,1.0_dp,error=error)
    ENDIF multiproc

  END SUBROUTINE make_full_kinetic_sainv


! *****************************************************************************
!> \brief DBCSR Kinetic Stable Approximate Inverse Preconditioner with
!> \param[in,out] preconditioner_env    The preconditioner environment
!> \param         pc                    The precondition matrix
!> \param         matrix_t              The MO matrix
!> \param         matrix_s              The overlap matrix
!> \param[in]     energy_homo           HOMO energy
!> \param[in]     energy_gap            energy gap
!> \param[in,out] error       cp2k error
!> \author Valery Weber
!> \par History
!>      - Created by VW
!>      - 2009-05-13 Converted to use DBCSR matrices
! *****************************************************************************
  SUBROUTINE make_sparse_kinetic_sainv(preconditioner_env,&
       pc, matrix_t, matrix_s, energy_homo, energy_gap, error)
    !
    TYPE(preconditioner_type), INTENT(INOUT) :: preconditioner_env
    TYPE(dbcsr_type), POINTER                :: pc
    TYPE(real_matrix_type), POINTER          :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sparse_kinetic_sainv', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: thresh_inner = 1.0E-2_dp, &
                                                thresh_post = 1.0E-5_dp

    INTEGER                                  :: i, j, n, nblkcols_total, &
                                                nblkrows_total, &
                                                nfullrows_total, nz, stat
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_col_blk_size
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: shift
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p, v
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: buf1, buf2
    TYPE(array_i1d_obj)                      :: col_blk_size_array, &
                                                col_col_blk_size_array, &
                                                row_blk_size_array
    TYPE(dbcsr_distribution_obj)             :: dist_col, dist_main
    TYPE(dbcsr_type)                         :: col1_m, col2_m, col2_t_m, &
                                                mat_t, p_m, v_m, w1_m, w2_m
    TYPE(distribution_2d_type), POINTER      :: dist2d_col, dist2d_main

!
!
!
! These threshould should be make specifiable in the input. -UB
!

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    CALL cp_assert (.NOT. ASSOCIATED (pc), cp_warning_level,&
         cp_wrong_args_error, routineP,&
         "Precondition matrix should not pre-exist.", error=error)
    ALLOCATE (pc, stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "pc", 1)
    !CALL dbcsr_from_sm (mat_t, matrix_t, error)
    !CALL dbcsr_from_sm (mat_s, matrix_s, error)
    !
    ! some infos
    CALL get_matrix_info(matrix_t,nrow=n)
    !
    ! allocate
    ALLOCATE(buf1(n,1),buf2(n,1),p(n),v(n))
    !CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
    !     context=preconditioner_env%ctxt, &
    !     para_env=preconditioner_env%para_env,error=error)
    !> T = T - epsilon * S
    shift=MAX(-energy_homo,energy_gap)
    CALL add_matrices(matrix_t,1.0_dp,matrix_t,shift,matrix_s)
    !CALL copy_sm_to_fm(matrix_t,fm_work1,error=error)
    CALL dbcsr_init (mat_t, error)
    CALL dbcsr_from_sm (mat_t, matrix_t, error=error)
    !
    NULLIFY (dist2d_main)
    CALL dbcsr_get_info (mat_t, row_blk_size=row_blk_size_array,&
         col_blk_size=col_blk_size_array, nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total, nfullrows_total=nfullrows_total,&
         distribution=dist_main, error=error)
    CALL cp_dbcsr_get_dist2d (mat_t, dist2d_main, error)
    col_blk_size => array_data (col_blk_size_array)
    ! Create the square matrices.
    CALL dbcsr_init (w1_m, error)
    CALL dbcsr_create (w1_m, "w1", dist_main, "N", row_blk_size_array,&
         col_blk_size_array, nze=nz, nblks=nblkrows_total, error=error)
    CALL dbcsr_work_create (w1_m, nblks_guess=nblkrows_total,&
         sizedata_guess=nz, error=error)
    CALL dbcsr_init (w2_m, error)
    CALL dbcsr_create (w2_m, "w2", dist_main, "N", row_blk_size_array,&
         col_blk_size_array, nze=nz, nblks=nblkrows_total, error=error)
    CALL dbcsr_work_create (w1_m, nblks_guess=nblkrows_total,&
         sizedata_guess=nz, error=error)
    ! Create the column matrices
    ! The col_col_blk_size array is modified to reflect the changing column
    ! widths.
    NULLIFY (dist2d_col)
    CALL dbcsr_create_dist2d_r (dist2d_col, dist2d_main, ncolumns=1,&
         error=error)
    CALL dbcsr_dist2d_to_dist (dist2d_col, dist_col, error)
    ALLOCATE (col_col_blk_size(1), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "col_col_blk_size", 1)
    col_col_blk_size(1) = 1
    CALL array_new (col_col_blk_size_array, col_col_blk_size, gift=.TRUE.)
    CALL dbcsr_init (col1_m, error)
    CALL dbcsr_create (col1_m, "col1", dist_col, "N", row_blk_size_array,&
         col_col_blk_size_array, nze=nz, nblks=nblkrows_total, error=error)
    CALL dbcsr_work_create (col1_m, nblks_guess=nblkrows_total,&
         sizedata_guess=nz, error=error)
    CALL dbcsr_addto_index_array(col1_m, dbcsr_slot_row_p,&
         reservation=nblkrows_total+1, error=error)
    CALL dbcsr_addto_index_array(col1_m, dbcsr_slot_col_i,&
         reservation=nblkrows_total, error=error)
    CALL dbcsr_addto_index_array(col1_m, dbcsr_slot_blk_p,&
         reservation=nblkrows_total, error=error)
    CALL dbcsr_ensure_data_size(col1_m, MAXVAL(col_blk_size)*nfullrows_total,&
         nocopy=.TRUE., error=error)
    CALL dbcsr_init (col2_m, error)
    CALL dbcsr_create (col2_m, "col2", dist_col, "N", row_blk_size_array,&
         col_col_blk_size_array, nze=nz, nblks=nblkrows_total, error=error)
    CALL dbcsr_work_create (col2_m, nblks_guess=nblkrows_total,&
         sizedata_guess=nz, error=error)
    CALL dbcsr_addto_index_array(col2_m, dbcsr_slot_row_p,&
         reservation=nblkrows_total+1, extra=nblkrows_total*2, error=error)
    CALL dbcsr_addto_index_array(col2_m, dbcsr_slot_col_i,&
         reservation=nblkrows_total, error=error)
    CALL dbcsr_addto_index_array(col2_m, dbcsr_slot_blk_p,&
         reservation=nblkrows_total, error=error)
    CALL dbcsr_ensure_data_size(col2_m, MAXVAL(col_blk_size)*nfullrows_total,&
         nocopy=.TRUE., error=error)
    CALL dbcsr_init (v_m, error)
    CALL dbcsr_create (v_m, "v_m", dist_col, "T", row_blk_size_array,&
         col_col_blk_size_array, nze=nz, nblks=nblkrows_total, error=error)
    CALL dbcsr_work_create (v_m, nblks_guess=nblkrows_total,&
         sizedata_guess=nz, error=error)
    CALL dbcsr_init (p_m, error)
    CALL dbcsr_create (p_m, "p_m", dist_col, "T", row_blk_size_array,&
         col_col_blk_size_array, nze=nz, nblks=nblkrows_total, error=error)
    CALL dbcsr_work_create (p_m, nblks_guess=nblkrows_total,&
         sizedata_guess=nz, error=error)
    CALL distribution_2d_release (dist2d_col, error)
    CALL dbcsr_distribution_release (dist_col)
    !
    ! Z = diag(ones(n,1));
    !!CALL cp_fm_set_all(fm_work2,0.0_dp,1.0_dp,error=error)
    CALL dbcsr_set_diagonal (w2_m, 1.0_dp, error=error)
    CALL dbcsr_finalize (w2_m, resort=.TRUE., reshuffle=.TRUE.,&
         error=error)
    !
    ! Let's ainv
    DO i = 1,nblkrows_total
       !
       ! v=A*Z(:,i);
       !!CALL cp_fm_get_submatrix(fm_work2,buf2 ,&
       !!     &                   1,i,n,1,transpose=.FALSE.,&
       !!     &                   error=error)     
       failure = .FALSE.
       CPPrecondition(i.GE.1.AND.i.LE.SIZE(col_blk_size),cp_warning_level,routineP,error,failure)
       col_col_blk_size(1) = col_blk_size(i)
       WRITE(*,*)'i',i,'column size',col_blk_size(i)
       CALL dbcsr_extract_column (col2_m, w2_m, i, error)
       !CALL dbcsr_print (col2_m, error)
       CALL dbcsr_verify_matrix (col2_m, error, verbosity=1)
       CALL dbcsr_new_transposed (col2_t_m, col2_m, error=error)
       CALL dbcsr_verify_matrix (col2_t_m, error, verbosity=1)
       !!CALL dgemv('N',n,n,1.0_dp,fm_work1%local_data(1,1),n,buf2(1,1),1,0.0_dp,v(1),1)
       !
       !!DO j = i,n
       !!   !
       !!   ! p(j) = v(:)'*Z(:,j);
       !!   !!CALL cp_fm_get_submatrix(fm_work2,buf2 ,&
       !!   !!     &                   1,j,n,1,transpose=.FALSE.,&
       !!   !!     &                   error=error)     
       !!   p(j) = DDOT(n,v(1:n),1,buf2(1:n,1),1)
       !!ENDDO
       !CALL
       !
       !WRITE (*,*)'Multiplying...'
       !CALL dbcsr_sm_sm_multiply (col2_t_m, w2_m, w1_m, col2_t_m%dist, 0,&
       !     error=error)
       IF(i.LT.n) THEN
          !!CALL cp_fm_get_submatrix(fm_work2,buf1 ,&
          !1     &                   1,i,n,1,transpose=.FALSE.,&
          !!     &                   error=error)
          DO j = i+1,n
             !
             ! Z(:,j) = Z(:,j)-(p(j)/p(i))*Z(:,i);
             !!CALL cp_fm_get_submatrix(fm_work2,buf2 ,&
             !1     &                   1,j,n,1,transpose=.FALSE.,&
             !!     &                   error=error)
             !IF(p(i).LT.0.0_dp) THEN
             !   write(*,*) 'make_kinetic_sainv: i=',i,' p(i)=',p(i)
             !   CALL stop_program("make_kinetic_sainv", "Ooops pivot negative")
             !ENDIF
             !IF(p(i).LT.1.0E-12_dp.AND.p(i).GE.0.0_dp) THEN
             !   write(*,*) 'make_kinetic_sainv: i=',i,' p(i)=',p(i)
             !   CALL stop_program("make_kinetic_sainv", "Ooops pivot too small")
             !ENDIF
             !buf2 = buf2 - p(j) / p(i) * buf1
             !
             ! filter...           
             !DO k=1,j-1
             !   IF(ABS(buf2(k,1)).LT.mx*thresh_inner) THEN
             !      buf2(k,1)=0.0_dp
             !   ENDIF
             !ENDDO
             !!CALL cp_fm_set_submatrix(fm_work2,buf2 ,&
             !!     &                   1,j,n,1,transpose=.FALSE.,&
             !!     &                   error=error)
          ENDDO
       ENDIF
    ENDDO
    !
    ! M = Z * D^(-1) * Z'
    !!DO i = 1,n
    !!   p(i) = 1.0_dp / p(i)
    !!ENDDO
    !CALL cp_fm_to_fm(fm_work2,fm_work1,error=error)
    !CALL cp_fm_column_scale(fm_work1,p)
    !CALL cp_fm_gemm('N','T',n,n,n,1.0_dp,fm_work2,fm_work1,0.0_dp,fm,error=error)
    !!
    !! post-filter...
    !CALL cp_fm_get_info(fm,nrow_local=nrow_local,ncol_local=ncol_local, &
    !     local_data=local_data,error=error)
    !nz=0
    !DO j=1,ncol_local
    !   DO i=1,nrow_local
    !      IF(i.NE.j) THEN
    !         IF(ABS(local_data(i,j)).LT.mx*thresh_post) THEN
    !            local_data(i,j)=0.0_dp
    !            nz=nz+1
    !         ENDIF
    !      ENDIF
    !   ENDDO
    !ENDDO
    !CALL mp_sum(nz,preconditioner_env%para_env%group)
    !write(*,*) 'non0(M)=',DBLE(n**2-nz)/DBLE(n**2)
    !mx = maxval(abs(local_data(1:nrow_local,1:ncol_local)))
    !CALL mp_max(mx,preconditioner_env%para_env%group)
    !write(*,*) 'max(abs(M))=',mx
    !!
    !! eig(M)
    !CALL cp_fm_to_fm(fm,fm_work1,error=error)
    !CALL cp_fm_syevd(fm_work1,fm_work2,v,error=error)
    !write(*,*) 'cond(M)=',MAXVAL(abs(v))/MAX(MINVAL(abs(v)),EPSILON(0.0_dp))
    !write(*,*) 'min(eig(M))=',MINVAL(v),'max(eig(M))=',MAXVAL(v)
    !!
    !!
    !DEALLOCATE(buf1,buf2,p,v)
    !CALL cp_fm_release(fm_work1,error=error)
    !CALL cp_fm_release(fm_work2,error=error)
    !
    CALL dbcsr_destroy (mat_t, error)
    CALL dbcsr_init (pc, error)
    pc = w2_m
    CALL cp_dbcsr_add_dist2d (pc, dist2d_main, error)
    CALL dbcsr_work_destroy (col1_m%w, error)
    CALL dbcsr_work_destroy (col2_m%w, error)
    CALL dbcsr_work_destroy (w1_m%w, error)
    CALL dbcsr_work_destroy (w2_m%w, error)
    CALL dbcsr_destroy (col1_m, error)
    CALL dbcsr_destroy (col2_m, error)
    CALL dbcsr_destroy (col2_t_m, error)
    CALL dbcsr_destroy (p_m, error)
    CALL dbcsr_destroy (v_m, error)
    CALL dbcsr_destroy (w1_m, error)
    ! CALL dbcsr_destroy (w2_m, error) ! pc takes it over
    CALL array_release (col_col_blk_size_array)
    CALL dbcsr_verify_matrix (pc, error)
    !CALL dbcsr_print (pc, error)
  END SUBROUTINE make_sparse_kinetic_sainv


  SUBROUTINE apply_sparse_single(preconditioner_env, matrix_in, matrix_out,&
       error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: k
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size
    TYPE(array_i1d_obj)                      :: precond_row_blk_size_array
    TYPE(dbcsr_distribution_obj)             :: dist_result
    TYPE(dbcsr_type)                         :: matrix_result, matrix_right
    TYPE(distribution_2d_type), POINTER      :: dist2d_precond, dist2d_result

!TYPE(dbcsr_type), INTENT(IN)             :: matrix_in
!TYPE(dbcsr_type), INTENT(OUT)            :: matrix_out
!----

    WRITE(*,*)' APPLYING SPARSE SINGLE KINETIC'
    !CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)
    !CALL cp_fm_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
    !     matrix_in,0.0_dp,matrix_out,error=error)
    CALL cp_fm_get_info(matrix_in, ncol_global=k, error=error)
    NULLIFY (dist2d_precond, dist2d_result, col_blk_size)
    CALL cp_dbcsr_get_dist2d (preconditioner_env%dbcsr_matrix, dist2d_precond,&
         error)
    CALL dbcsr_create_dist2d_r_unrot(dist2d_result, dist2d_precond, k,&
         col_blk_size, error=error)
    CALL dbcsr_dist2d_to_dist (dist2d_result, dist_result, error)
    CALL dbcsr_init (matrix_right, error)
    precond_row_blk_size_array =&
         preconditioner_env%dbcsr_matrix%row_blk_size
    CALL dbcsr_from_fm (matrix_right, matrix_in, 0.0d0, dist2d_result,&
         row_blk_size = array_data (precond_row_blk_size_array),&
         col_blk_size = col_blk_size, error=error)
    WRITE(*,*)'Pre checksum:',dbcsr_checksum (matrix_right, error)
    !CALL dbcsr_get_info (matrix_in, dist2d=dist2d, error=error)
    CALL dbcsr_init (matrix_result, error)
    !CALL dbcsr_print (preconditioner_env%dbcsr_matrix, error)
    !CALL dbcsr_print (matrix_right, error)
    CALL dbcsr_sm_sm_multiply (preconditioner_env%dbcsr_matrix,&
         matrix_right, matrix_result, dist_result, 0, error)
    ! CALL dbcsr_to_fm (matrix_result, matrix_out, error)
    !CALL dbcsr_verify_matrix (matrix_result, error)
    CALL dbcsr_distribution_release (dist_result)
    CALL dbcsr_destroy (matrix_right, error)
    WRITE(*,*)'Post checksum:',dbcsr_checksum (matrix_result, error)
    CALL dbcsr_destroy (matrix_result, error)
    CALL distribution_2d_release (dist2d_result, error)
    WRITE(*,*)' DONE'
  END SUBROUTINE apply_sparse_single



  SUBROUTINE make_sparse_kinetic(preconditioner_env, sm, matrix_t, matrix_s, &
       &                         energy_homo, energy_gap, error)
    !
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(real_matrix_type), POINTER          :: sm, matrix_t, matrix_s
    REAL(KIND=dp), INTENT(IN)                :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sparse_kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: shift

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    ! some infos
    CALL get_matrix_info(matrix_t,nrow=n)
    ! allocate
    CALL Replicate_matrix(matrix_t,sm,"PRECONDITIONER",.TRUE.,error=error)
    ! M = T - epsilon * S
    shift=MAX(-energy_homo,energy_gap)
    CALL add_matrices(sm,1.0_dp,matrix_t,shift,matrix_s)
  END SUBROUTINE make_sparse_kinetic

END MODULE preconditioner

