!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief solves the preconditioner, contains to utility function for 
!>        fm<->dbcsr transfers, should be moved soon
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner_solvers
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_arnoldi_ev, cp_dbcsr_col_block_sizes, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_get_info, cp_dbcsr_init_p, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_setup_arnoldi_data, cp_dbcsr_type, dbcsr_arnoldi_data, &
       dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_default, &
       deallocate_arnoldi_data, get_selected_ritz_val
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: ot_precond_solver_default,&
                                             ot_precond_solver_direct,&
                                             ot_precond_solver_inv_chol,&
                                             ot_precond_solver_update
  USE iterate_matrix,                  ONLY: invert_Hotelling
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_solvers'

  PUBLIC :: solve_preconditioner, transfer_fm_to_dbcsr, transfer_dbcsr_to_fm

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param my_solver_type ...
!> \param preconditioner_env ...
!> \param matrix_s ...
!> \param matrix_h ...
!> \param my_mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE solve_preconditioner(my_solver_type, preconditioner_env, matrix_s, &
                                  matrix_h, my_mixed_precision, error)
    INTEGER                                  :: my_solver_type
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    LOGICAL                                  :: my_mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'solve_preconditioner', &
      routineP = moduleN//':'//routineN

! here comes the solver

    SELECT CASE(my_solver_type)
    CASE (ot_precond_solver_inv_chol)
       !
       ! compute the full inverse
       preconditioner_env%solver=ot_precond_solver_inv_chol
       CALL make_full_inverse_cholesky(preconditioner_env, matrix_s, my_mixed_precision, error=error)
    CASE (ot_precond_solver_direct)
       !
       ! prepare for the direct solver
       preconditioner_env%solver=ot_precond_solver_direct
       CALL make_full_fact_cholesky(preconditioner_env, matrix_s, my_mixed_precision, error)
    CASE (ot_precond_solver_update)
       !
       ! uses an update of the full inverse (needs to be computed the first time)
       ! make sure preconditioner_env is not destroyed in between
       IF(preconditioner_env%solver.NE.ot_precond_solver_update)THEN
          preconditioner_env%solver=ot_precond_solver_update
          CALL make_full_inverse_cholesky(preconditioner_env, matrix_s, my_mixed_precision, error=error)
       ELSE
          preconditioner_env%solver=ot_precond_solver_update
          CALL make_inverse_update(preconditioner_env, matrix_h, my_mixed_precision, matrix_s, error)
       END IF
    CASE (ot_precond_solver_default)
       preconditioner_env%solver=ot_precond_solver_default
    CASE DEFAULT
       !
       CALL stop_program(routineN, moduleN, __LINE__, "Doesn't know this type of solver")
    END SELECT

  END SUBROUTINE solve_preconditioner

! *****************************************************************************
!> \brief Compute the inverse using cholseky factorization
!> \param preconditioner_env ...
!> \param matrix_s ...
!> \param mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_inverse_cholesky(preconditioner_env, matrix_s, mixed_precision, &
       error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_inverse_cholesky', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_fm_type), POINTER                :: fm, fm_work

    failure = .FALSE.

    CALL timeset(routineN, handle)

    ! Maybe we will get a sparse Cholesky at a given point then this can go,
    ! if stuff was stored in fm anyway this simple returns
    CALL transfer_dbcsr_to_fm(preconditioner_env%sparse_matrix, preconditioner_env%fm, &
                              preconditioner_env%para_env, preconditioner_env%ctxt, &
                              mixed_precision, error)
    fm=>preconditioner_env%fm

    NULLIFY(fm_work)

    CALL cp_fm_create(fm_work, fm%matrix_struct, name="fm_work", use_sp=mixed_precision, error=error)
    !
    ! compute the inverse of SPD matrix fm using the Cholesky factorization
    CALL cp_error_init(sub_error, template_error=error, stop_level=cp_fatal_level)

    CALL cp_fm_cholesky_decompose(fm, error=sub_error)

    failure = .FALSE.
    CALL cp_error_check(sub_error, failure)
    CALL cp_error_dealloc_ref(sub_error)
    !
    ! if fm not SPD we go with the overlap matrix
    IF (failure) THEN
       !
       ! just the overlap matrix
       IF(PRESENT(matrix_s)) THEN
          CALL copy_dbcsr_to_fm(matrix_s, fm, error=error)
          CALL cp_fm_cholesky_decompose(fm, error=error)
       ELSE
          CALL cp_fm_set_all(fm, alpha=0._dp, beta=1._dp, error=error)
       ENDIF
    ENDIF
    CALL cp_fm_cholesky_invert(fm, error=error)

    CALL cp_fm_upper_to_full(fm, fm_work, error=error)
    CALL cp_fm_release(fm_work, error=error)

    CALL timestop(handle)

  END SUBROUTINE make_full_inverse_cholesky

! *****************************************************************************
!> \brief Only perform the factorization, can be used later to solve the linear
!>        system on the fly
!> \param preconditioner_env ...
!> \param matrix_s ...
!> \param mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_fact_cholesky(preconditioner_env, matrix_s, mixed_precision, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_fact_cholesky', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_fm_type), POINTER                :: fm

    failure = .FALSE.

    CALL timeset(routineN, handle)

    ! Maybe we will get a sparse Cholesky at a given point then this can go,
    ! if stuff was stored in fm anyway this simple returns
    CALL transfer_dbcsr_to_fm(preconditioner_env%sparse_matrix, preconditioner_env%fm, &
                              preconditioner_env%para_env, preconditioner_env%ctxt, &
                              mixed_precision, error)

    fm=>preconditioner_env%fm
    !
    ! compute the inverse of SPD matrix fm using the Cholesky factorization
    CALL cp_error_init(sub_error, template_error=error, stop_level=cp_fatal_level)
    CALL cp_fm_cholesky_decompose(fm, error=sub_error)
    failure = .FALSE.
    CALL cp_error_check(sub_error, failure)
    CALL cp_error_dealloc_ref(sub_error)
    !
    ! if fm not SPD we go with the overlap matrix
    IF (failure) THEN
       !
       ! just the overlap matrix
       IF(PRESENT(matrix_s)) THEN
          CALL copy_dbcsr_to_fm(matrix_s, fm, error=error)
          CALL cp_fm_cholesky_decompose(fm, error=error)
       ELSE
          CALL cp_fm_set_all(fm, alpha=0._dp, beta=1._dp, error=error)
       ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE make_full_fact_cholesky

! *****************************************************************************
!> \brief computes an approximate inverse using Hotelling iterations
!> \param preconditioner_env ...
!> \param matrix_h as S is not always present this is a safe template for the transfer
!> \param mixed_precision ...
!> \param matrix_s ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_inverse_update(preconditioner_env, matrix_h, mixed_precision, matrix_s, error)
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_inverse_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cond_num, filter_eps

    CALL timeset(routineN, handle)
    !
    ! uses an update of the full inverse (needs to be computed the first time)
    ! make sure preconditioner_env is not destroyed in between

    CPPrecondition(ASSOCIATED(preconditioner_env%dbcsr_matrix), cp_failure_level, routineP, error, failure)
    ! Maybe I gonna add a fm Hotelling, ... for now the same as above make sure we are dbcsr
    CALL transfer_fm_to_dbcsr(preconditioner_env%fm, preconditioner_env%sparse_matrix, &
                              matrix_h, mixed_precision, error)

    ! Try to get a reasonbale guess for the filtering threshold
    IF(PRESENT(matrix_s))CALL estimate_cond_num(matrix_s, cond_num, error)
    filter_eps=1.0_dp/cond_num*0.01_dp
    ! We don't need a high accuracy for the inverse so 0.4 is reasonable for convergence 
    CALL invert_Hotelling(preconditioner_env%dbcsr_matrix, preconditioner_env%sparse_matrix, filter_eps, &
                          use_inv_as_guess=.TRUE., norm_convergence=0.4_dp, filter_eps=filter_eps, error=error)

    CALL timestop(handle)

  END SUBROUTINE make_inverse_update


! *****************************************************************************
!> \brief gives an estimate for the condition number of a matrix, should be moved
!>        to a more general place
!> \param matrix ...
!> \param cond_num ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE estimate_cond_num(matrix, cond_num, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    REAL(KIND=dp)                            :: cond_num
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp)                            :: max_ev, min_ev
    TYPE(cp_dbcsr_p_type), DIMENSION(1)      :: mymat
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi

    mymat(1)%matrix=>matrix
    CALL cp_dbcsr_setup_arnoldi_data(arnoldi, mymat, max_iter=50, threshold=1.0E-4_dp, &
                                     selection_crit=1, nval_request=2)
    CALL cp_dbcsr_arnoldi_ev(mymat, arnoldi, error)
    max_eV=REAL(get_selected_ritz_val(arnoldi, 2), dp)
    min_eV=REAL(get_selected_ritz_val(arnoldi, 1), dp)
    cond_num=max_eV/min_eV
    CALL deallocate_arnoldi_data(arnoldi)

   END SUBROUTINE estimate_cond_num

! *****************************************************************************
!> \brief transfers a full matrix to a dbcsr
!> \param fm_matrix a full matrix gets deallocated in the end
!> \param dbcsr_matrix a dbcsr matrix, gets create from a template
!> \param template_mat the template which is used for the structure
!> \param mixed_precision  whether to use single or double precision
!> \param error ...
! *****************************************************************************
  SUBROUTINE transfer_fm_to_dbcsr(fm_matrix, dbcsr_matrix, template_mat, mixed_precision, error)

    TYPE(cp_fm_type), POINTER                :: fm_matrix
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_matrix, template_mat
    LOGICAL                                  :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transfer_fm_to_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)
    IF(ASSOCIATED(fm_matrix))THEN
       IF(.NOT.ASSOCIATED(dbcsr_matrix)) THEN
          CALL cp_dbcsr_init_p(dbcsr_matrix, error=error)
          IF(mixed_precision) THEN
             CALL cp_dbcsr_create(dbcsr_matrix, "transfered_matrix", &
                  cp_dbcsr_distribution(template_mat), dbcsr_type_no_symmetry, &
                  cp_dbcsr_row_block_sizes(template_mat), cp_dbcsr_col_block_sizes(template_mat), &
                  0, 0, dbcsr_type_real_4, error=error)
          ELSE
             CALL cp_dbcsr_create(dbcsr_matrix, "preconditioner_env%dbcsr_matrix", &
                  cp_dbcsr_distribution(template_mat), dbcsr_type_no_symmetry, &
                  cp_dbcsr_row_block_sizes(template_mat), cp_dbcsr_col_block_sizes(template_mat), &
                  0, 0, dbcsr_type_real_default, error=error)
          ENDIF
       ENDIF
!       CALL cp_fm_create(fm_tmp, matrix_struct=fm_matrix%matrix_struct, error=error)
!       CALL cp_fm_upper_to_full(fm_matrix, fm_tmp, error)
       CALL copy_fm_to_dbcsr(fm_matrix, dbcsr_matrix, error=error)
!       CALL cp_fm_release(fm_tmp, error)
       CALL cp_fm_release(fm_matrix, error)
    END IF
  
    CALL timestop(handle)
  
  END SUBROUTINE transfer_fm_to_dbcsr

! *****************************************************************************
!> \brief transfers a dbcsr to a full matrix
!> \param dbcsr_matrix a dbcsr matrix, gets deallocated at the end
!> \param fm_matrix a full matrix gets created if not yet done
!> \param para_env the para_env
!> \param context the blacs context
!> \param mixed_precision  whether to use single or double precision
!> \param error ...
! *****************************************************************************
  SUBROUTINE transfer_dbcsr_to_fm(dbcsr_matrix, fm_matrix, para_env, context, mixed_precision, error)

    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_matrix
    TYPE(cp_fm_type), POINTER                :: fm_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: context
    LOGICAL                                  :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transfer_dbcsr_to_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp

    CALL timeset(routineN, handle)
    IF(ASSOCIATED(dbcsr_matrix))THEN
       NULLIFY(fm_struct_tmp)

       IF (ASSOCIATED(fm_matrix)) THEN
          CALL cp_fm_release(fm_matrix, error=error)
       ENDIF

       CALL cp_dbcsr_get_info(dbcsr_matrix, nfullrows_total=n)
       CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n, ncol_global=n, &
                                  context=context, para_env=para_env, error=error)
       CALL cp_fm_create(fm_matrix, fm_struct_tmp, use_sp=mixed_precision , error=error)
       CALL cp_fm_struct_release(fm_struct_tmp, error=error)

       CALL copy_dbcsr_to_fm(dbcsr_matrix, fm_matrix, error=error)
       CALL cp_dbcsr_release(dbcsr_matrix, error)
       DEALLOCATE(dbcsr_matrix)
    END IF

    CALL timestop(handle)

  END SUBROUTINE transfer_dbcsr_to_fm

END MODULE preconditioner_solvers
