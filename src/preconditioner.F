!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_distribution, &
       cp_dbcsr_get_info, cp_dbcsr_get_matrix_type, cp_dbcsr_init, &
       cp_dbcsr_init_p, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_release_p, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_default
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: choose_eigv_solver
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_diag,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE input_constants,                 ONLY: &
       ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
       ot_precond_solver_default, ot_precond_solver_direct, &
       ot_precond_solver_inv_chol
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             sp
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             preconditioner_p_type,&
                                             preconditioner_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'

  PUBLIC :: make_preconditioner, restart_preconditioner
  PUBLIC :: apply_preconditioner, prepare_preconditioner

  PRIVATE :: make_full_single

  INTERFACE apply_preconditioner
     MODULE PROCEDURE apply_preconditioner_dbcsr
     MODULE PROCEDURE apply_preconditioner_fm
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

! creates a preconditioner for the system (H-energy_homo S)
! this preconditioner is (must be) symmetric positive definite.
! currently uses a atom-block-diagonal form
! each block will be  ....
! might overwrite matrix_h, matrix_t

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param precon_type ...
!> \param solver_type ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param matrix_t ...
!> \param mo_set ...
!> \param energy_gap ...
!> \param mixed_precision ...
!> \param convert_precond_to_dbcsr ...
!> \param error ...
!> \par History
!>      09.2014 removed some unused or unfinished methods
!>              removed sparse preconditioners and the 
!>              sparse approximate inverse at rev 14341 [Florian Schiffmann]
! *****************************************************************************
SUBROUTINE make_preconditioner(preconditioner_env, precon_type, solver_type, matrix_h, matrix_s, &
     matrix_t, mo_set, energy_gap, mixed_precision, convert_precond_to_dbcsr, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    INTEGER, INTENT(IN)                      :: precon_type, solver_type
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s, matrix_t
    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp)                            :: energy_gap
    LOGICAL, INTENT(IN), OPTIONAL            :: mixed_precision, &
                                                convert_precond_to_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, k, &
                                                my_solver_type
    LOGICAL :: failure, my_convert_precond_to_dbcsr, my_mixed_precision, &
      needs_full_spectrum, needs_homo, use_mo_coeff_b
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
    use_mo_coeff_b = mo_set%use_mo_coeff_b
    CALL cp_fm_get_info(mo_coeff,ncol_global=k,error=error)

    my_convert_precond_to_dbcsr = .FALSE.
    IF(PRESENT(convert_precond_to_dbcsr)) my_convert_precond_to_dbcsr = convert_precond_to_dbcsr

    my_mixed_precision = .FALSE.
    IF(PRESENT(mixed_precision)) my_mixed_precision = mixed_precision
    IF(my_mixed_precision) THEN
       SELECT CASE(precon_type)
       CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse)
          !supported
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "This precond with mixed precision is not supported yet")
       END SELECT
    ENDIF

    needs_full_spectrum=.FALSE.
    needs_homo=.FALSE.

    SELECT CASE(precon_type)
    CASE (ot_precond_full_single_inverse,ot_precond_full_all)
       needs_full_spectrum=.TRUE.
    CASE (ot_precond_full_single)
       needs_homo=.TRUE.
       ! XXXX to be removed if homo estimate only is implemented
       needs_full_spectrum=.TRUE.
    CASE (ot_precond_full_kinetic,ot_precond_s_inverse)
       ! these should be happy without an estimate for the homo energy
       ! preconditioning can  not depend on an absolute eigenvalue, only on eigenvalue differences
    CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "The preconditioner is unknown ...")
    END SELECT

    energy_homo=0.0_dp
    IF (needs_full_spectrum) THEN
       ALLOCATE(eigenvalues_ot(k),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ! XXXXXXXXXXXXXXXX do not touch the initial MOs, could be harmful for either
       !                  the case of non-equivalent MOs but also for the derivate
       ! we could already have all eigenvalues e.g. full_all and we could skip this
       ! to be optimised later.
       ! one flaw is that not all SCF methods (i.e. that go over mo_derivs directly)
       ! have a 'valid' matrix_h... (we even don't know what evals are in that case)
       IF(use_mo_coeff_b) THEN
          CALL calculate_subspace_eigenvalues(mo_coeff_b,matrix_h,&
               eigenvalues_ot, do_rotation = .FALSE.,&
               para_env=mo_coeff%matrix_struct%para_env,&
               blacs_env=mo_coeff%matrix_struct%context,error=error)
       ELSE
          CALL calculate_subspace_eigenvalues(mo_coeff,matrix_h,&
               eigenvalues_ot, do_rotation = .FALSE.,error=error)
       ENDIF
       IF (k>0) energy_homo=eigenvalues_ot(k)
    ELSE
       IF (needs_homo) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Not yet implemented")
       ENDIF
    ENDIF

  my_solver_type = solver_type
  SELECT CASE (precon_type)
  CASE (ot_precond_full_single)
     preconditioner_env%in_use=ot_precond_full_single
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Only PRECOND_SOLVER DEFAULT for the moment")
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
     !   CALL stop_program(routineN,moduleN,__LINE__,"Only PRECOND_SOLVER DEFAULT for the moment")
     !ENDIF
     IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, error=error)
     ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Type for S=1 not implemented")
     END IF
  CASE (ot_precond_full_kinetic)
     preconditioner_env%in_use=ot_precond_full_kinetic
     IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
     IF ( PRESENT(matrix_s) .AND. PRESENT(matrix_t) ) THEN
        CALL make_full_kinetic(preconditioner_env, preconditioner_env%fm,&
                               matrix_t, matrix_s, energy_gap, &
                               my_mixed_precision, error=error)
     ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Type for S=1 not implemented")
     ENDIF
  CASE (ot_precond_full_single_inverse)
     preconditioner_env%in_use=ot_precond_full_single_inverse
     IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
     IF(use_mo_coeff_b) THEN
        CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff,error=error)
     ENDIF
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_single_inverse(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                                      eigenvalues_ot, energy_gap, my_mixed_precision, &
                                      error=error)
     ELSE
        CALL make_full_single_inverse_ortho(preconditioner_env,mo_coeff,matrix_h, &
                                            eigenvalues_ot, energy_gap,error=error)
     END IF
  CASE (ot_precond_full_all)
     preconditioner_env%in_use=ot_precond_full_all
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     IF(use_mo_coeff_b) THEN
        CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff,error=error)
     ENDIF
     IF ( PRESENT(matrix_s) ) THEN
        CALL make_full_all(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                           eigenvalues_ot, energy_gap,error=error)
     ELSE
        CALL make_full_all_ortho(preconditioner_env,mo_coeff,matrix_h, &
                                 eigenvalues_ot, energy_gap,error=error)
     END IF
  CASE DEFAULT
     CALL stop_program(routineN,moduleN,__LINE__,"Type not implemented")
  END SELECT
  !
  ! here comes the solver
  SELECT CASE(my_solver_type)
  CASE (ot_precond_solver_inv_chol)
     !
     ! compute the full inverse
     preconditioner_env%solver=ot_precond_solver_inv_chol
     CALL make_full_inverse_cholesky(preconditioner_env, preconditioner_env%fm, matrix_s, &
                                     my_mixed_precision, error=error)
  CASE (ot_precond_solver_direct)
     !
     ! prepare for the direct solver
     preconditioner_env%solver=ot_precond_solver_direct
     CALL make_full_fact_cholesky(preconditioner_env, preconditioner_env%fm, matrix_s, &
                                  error)
   CASE (ot_precond_solver_default)
     preconditioner_env%solver=ot_precond_solver_default
  CASE DEFAULT
     !
     CALL stop_program(routineN,moduleN,__LINE__,"Doesn't know this type of solver")
  END SELECT


  IF(my_convert_precond_to_dbcsr) THEN
     IF(ASSOCIATED(preconditioner_env%dbcsr_matrix)) THEN
        CALL cp_dbcsr_release_p(preconditioner_env%dbcsr_matrix, error=error)
     ENDIF
     CALL cp_dbcsr_init_p(preconditioner_env%dbcsr_matrix,error=error)
     IF(my_mixed_precision) THEN
        CALL cp_dbcsr_create(preconditioner_env%dbcsr_matrix, "preconditioner_env%dbcsr_matrix", &
             cp_dbcsr_distribution(matrix_h), dbcsr_type_no_symmetry,&
             cp_dbcsr_row_block_sizes(matrix_h), cp_dbcsr_col_block_sizes(matrix_h), &
             0, 0, dbcsr_type_real_4, error=error)
     ELSE
        CALL cp_dbcsr_create(preconditioner_env%dbcsr_matrix, "preconditioner_env%dbcsr_matrix", &
             cp_dbcsr_distribution(matrix_h), dbcsr_type_no_symmetry,&
             cp_dbcsr_row_block_sizes(matrix_h), cp_dbcsr_col_block_sizes(matrix_h), &
             0, 0, dbcsr_type_real_default, error=error)
     ENDIF
     IF(ASSOCIATED(preconditioner_env%fm)) THEN
        CALL copy_fm_to_dbcsr(preconditioner_env%fm,preconditioner_env%dbcsr_matrix, error=error)
     ELSEIF(ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
        CALL cp_dbcsr_copy(preconditioner_env%dbcsr_matrix, preconditioner_env%sparse_matrix, error=error)
     ELSE
        CALL stop_program(routineN,moduleN,__LINE__,"Something is wrong")
     ENDIF
  ENDIF


    IF (needs_full_spectrum) THEN
      DEALLOCATE(eigenvalues_ot,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE make_preconditioner


! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_t ...
!> \param matrix_s ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_sparse_kinetic(preconditioner_env, matrix_t, matrix_s, &
                                 energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_t, matrix_s
    REAL(KIND=dp), INTENT(IN)                :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sparse_kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: shift

    failure = .FALSE.
  CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    ! some infos
    CALL cp_dbcsr_get_info(matrix_t,nfullrows_total=n)
    ! allocate

    IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       CALL cp_dbcsr_deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
       NULLIFY(preconditioner_env%sparse_matrix)
    ENDIF

    ALLOCATE(preconditioner_env%sparse_matrix)
    CALL cp_dbcsr_init(preconditioner_env%sparse_matrix,error=error)
    CALL cp_dbcsr_create(preconditioner_env%sparse_matrix, ' PRECONDITIONER ', &
         cp_dbcsr_distribution (matrix_t),&
         cp_dbcsr_get_matrix_type (matrix_t), cp_dbcsr_row_block_sizes(matrix_t),&
         cp_dbcsr_col_block_sizes(matrix_t), 0, 0, error=error)

    ! M = T - epsilon * S
    shift=MAX(0.0_dp,energy_gap)
    CALL cp_dbcsr_add(preconditioner_env%sparse_matrix,matrix_t,&
         alpha_scalar=0.0_dp,beta_scalar=1.0_dp,error=error)
    CALL cp_dbcsr_add(preconditioner_env%sparse_matrix,matrix_s,&
         alpha_scalar=1.0_dp,beta_scalar=shift,error=error)
  CALL timestop(handle)

  END SUBROUTINE make_sparse_kinetic


! applies a previously created preconditioner to a full matrix
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_preconditioner_fm(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

  CALL timeset(routineN,handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
     CALL stop_program(routineN,moduleN,__LINE__,"No preconditioner in use")
  CASE (ot_precond_full_single)
     CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_all)
     CALL apply_full_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
     SELECT CASE (preconditioner_env%solver)
     CASE(ot_precond_solver_inv_chol)
        CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE(ot_precond_solver_direct)
        CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,"Solver not implemented")
     END SELECT
  CASE DEFAULT
     CALL stop_program(routineN,moduleN,__LINE__,"Unknown preconditioner")
  END SELECT

  CALL timestop(handle)

END SUBROUTINE apply_preconditioner_fm

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_preconditioner_dbcsr(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

  CALL timeset(routineN,handle)

  SELECT CASE (preconditioner_env%in_use)
  CASE (0)
     CALL stop_program(routineN,moduleN,__LINE__,"No preconditioner in use")
  CASE (ot_precond_full_single)
     CALL apply_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_all)
     CALL apply_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
     SELECT CASE (preconditioner_env%solver)
     CASE(ot_precond_solver_inv_chol)
        CALL apply_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE(ot_precond_solver_direct)
        CALL stop_program(routineN,moduleN,__LINE__,"Apply_full_direct not supported with ot")
        !CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,"Wrong solver")
     END SELECT
  CASE DEFAULT
     CALL stop_program(routineN,moduleN,__LINE__,"Wrong preconditioner")
  END SELECT

  CALL timestop(handle)

END SUBROUTINE apply_preconditioner_dbcsr

! different types of preconditioner come here
! a full matrix preconditioner
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n

  CALL timeset(routineN,handle)

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)
  CALL cp_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_in,0.0_dp,matrix_out,error=error)
  CALL timestop(handle)

END SUBROUTINE apply_full_single

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_single(preconditioner_env, matrix_in, matrix_out,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

  CALL timeset(routineN,handle)

  IF (.NOT.ASSOCIATED(preconditioner_env%dbcsr_matrix)) &
     CALL stop_program(routineN,moduleN,__LINE__,&
                       "NOT ASSOCIATED preconditioner_env%dbcsr_matrix")
  CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%dbcsr_matrix,matrix_in,&
       0.0_dp,matrix_out,error=error)

  CALL timestop(handle)

END SUBROUTINE apply_single

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_full_direct(preconditioner_env, matrix_in, matrix_out,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_direct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n
    TYPE(cp_fm_type), POINTER                :: work

  CALL timeset(routineN,handle)

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)
  CALL cp_fm_create(work,matrix_in%matrix_struct,name="apply_full_single",&
                    use_sp=matrix_in%use_sp,error=error)
  CALL cp_fm_cholesky_restore(matrix_in,k,preconditioner_env%fm,work,&
       &                      "SOLVE",transa="T",error=error)
  CALL cp_fm_cholesky_restore(work,k,preconditioner_env%fm,matrix_out,&
       &                      "SOLVE",transa="N",error=error)
  CALL cp_fm_release(work,error=error)

  CALL timestop(handle)

END SUBROUTINE apply_full_direct


! *****************************************************************************
!> \brief generates a preconditioner by cholesky inverting H-lambda S+(SC)shifts(SC)^T
!> \param preconditioner_env ...
!> \param matrix_c0 must be already rotated correctly.
!> \param matrix_h ...
!> \param matrix_s ...
!> \param c0_evals ...
!> \param energy_gap ...
!> \param mixed_precision ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      10.2006 made more robust [Joost VandeVondele]
!> \note
!>      this might fail if the initial guess is bad, or if the system has 'holes'
! *****************************************************************************
SUBROUTINE make_full_single_inverse(preconditioner_env, matrix_c0, matrix_h, matrix_s, c0_evals, &
                                    energy_gap, mixed_precision, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_single_inverse', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: eval_shift = 5.0_dp , &
                                                fudge_factor = 2.0_dp

    INTEGER                                  :: handle, k, n
    REAL(KIND=dp)                            :: error_estimate, &
                                                preconditioner_shift
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_shift, diag, evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER :: matrix_hc0, matrix_hc0_sp, matrix_ptr, &
      matrix_s1, matrix_sc0, matrix_sc0_sp, matrix_shc0, matrix_tmp, &
      matrix_tmp2

  CALL timeset(routineN,handle)

! arbitrary upshift of the occupied evals
! fudge factor for taking the error estimate into account

    CALL cp_fm_get_info(matrix_c0,nrow_global=n,ncol_global=k,error=error)

    IF (ASSOCIATED(preconditioner_env%fm)) CALL cp_fm_release(preconditioner_env%fm,error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=n,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp,name="preconditioner_env%fm",&
                      use_sp=mixed_precision,error=error)
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp, name="preconditioner matrix_tmp",&
                      use_sp=mixed_precision,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! first try to get a ritz error estimate
    ! 0) cholesky decompose the overlap matrix, if this fails the basis is singular,
    !    more than EPS_DEFAULT
    CALL copy_dbcsr_to_fm(matrix_s,matrix_tmp,error=error)
    CALL cp_fm_cholesky_decompose(matrix_tmp,error=error)

    ! get the error estimate
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0",error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_s,matrix_c0,matrix_sc0,k,error=error)
    CALL cp_fm_create(matrix_hc0,matrix_c0%matrix_struct,name="hc0",error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

    IF(mixed_precision) THEN
       CALL cp_fm_create(matrix_sc0_sp,matrix_c0%matrix_struct,name="sc0_sp",&
                         use_sp=mixed_precision,error=error)
       CALL cp_fm_create(matrix_hc0_sp,matrix_c0%matrix_struct,name="hc0_sp",&
                         use_sp=mixed_precision,error=error)
       CALL cp_fm_to_fm(matrix_sc0,matrix_sc0_sp,error=error)
       CALL cp_fm_to_fm(matrix_hc0,matrix_hc0_sp,error=error)
       matrix_ptr => matrix_sc0; matrix_sc0 => matrix_sc0_sp; matrix_sc0_sp => matrix_ptr
       CALL cp_fm_release(matrix_sc0_sp,error=error)
       matrix_ptr => matrix_hc0; matrix_hc0 => matrix_hc0_sp; matrix_hc0_sp => matrix_ptr
       CALL cp_fm_release(matrix_hc0_sp,error=error)
    ENDIF

    CALL cp_fm_create(matrix_shc0,matrix_c0%matrix_struct,name="shc0",&
                      use_sp=mixed_precision,error=error)
    CALL cp_fm_cholesky_restore(matrix_hc0,k,matrix_tmp,matrix_shc0,"SOLVE",transa="T",error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",use_sp=mixed_precision,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ! since we only use diagonal elements this is a bit of a waste
    CALL cp_gemm('T','N',k,k,n,1.0_dp,matrix_shc0,matrix_shc0,0.0_dp,matrix_s1,error=error)
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
    CALL cp_gemm('N','T',n,n,k,1.0_dp,matrix_sc0,matrix_sc0,0.0_dp,preconditioner_env%fm,error=error)
    CALL cp_fm_release(matrix_sc0,error=error)
    DEALLOCATE(c0_shift)

    ! get H added to the shift
    CALL copy_dbcsr_to_fm(matrix_h,matrix_tmp,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,preconditioner_env%fm,1.0_dp,matrix_tmp,error=error)

    ! preconditioner shift, we target the middle of the occupied spectrum, and taking into account the error_estimate
    ! write(*,*) "Error estimate = ",error_estimate
    preconditioner_shift=-(MINVAL(c0_evals)+ MAXVAL(c0_evals))/2.0_dp + &
                           error_estimate*fudge_factor
    CALL copy_dbcsr_to_fm(matrix_s,matrix_tmp,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,preconditioner_env%fm,preconditioner_shift,matrix_tmp,error=error)

    ! check evals
    IF (.FALSE.) THEN
       CALL cp_fm_to_fm(preconditioner_env%fm,matrix_tmp,error=error)
       CALL cp_fm_create(matrix_tmp2,matrix_tmp%matrix_struct,name="matrix_tmp2",error=error)
       ALLOCATE(evals(n))
       CALL choose_eigv_solver(matrix_tmp,matrix_tmp2,evals,error=error)

       CALL cp_fm_release(matrix_tmp2,error=error)
       WRITE(*,*) "evals ",evals
       DEALLOCATE(evals)
    ENDIF

    CALL cp_fm_release(matrix_tmp,error=error)
  CALL timestop(handle)

END SUBROUTINE make_full_single_inverse

! *****************************************************************************
!> \brief generates a state by state preconditioner based on the full hamiltonian matrix
!> \param preconditioner_env ...
!> \param matrix_c0 ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param c0_evals ...
!> \param energy_gap should be a slight underestimate of the physical energy gap for almost all systems
!>      the c0 are already ritz states of (h,s)
!> \param error ...
!> \par History
!>      10.2006 made more stable [Joost VandeVondele]
!> \note
!>      includes error estimate on the hamiltonian matrix to result in a stable preconditioner
!>      a preconditioner for each eigenstate i is generated by keeping the factorized form
!>      U diag( something i ) U^T. It is important to only precondition in the subspace orthogonal to c0.
!>      not only is it the only part that matters, it also simplifies the computation of
!>      the lagrangian multipliers in the OT minimization  (i.e. if the c0 here is different
!>      from the c0 used in the OT setup, there will be a bug).
! *****************************************************************************
SUBROUTINE make_full_all(preconditioner_env, matrix_c0, matrix_h, matrix_s, c0_evals, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_all', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: fudge_factor = 0.25_dp, &
                                                lambda_base = 10.0_dp

    INTEGER                                  :: handle, k, n
    REAL(KIND=dp)                            :: error_estimate, lambda
    REAL(KIND=dp), DIMENSION(:), POINTER     :: diag, norms, shifted_evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER :: matrix_hc0, matrix_left, matrix_pre, &
      matrix_s1, matrix_s2, matrix_sc0, matrix_shc0, matrix_tmp, ortho

  CALL timeset(routineN,handle)

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
    CALL copy_dbcsr_to_fm(matrix_s,ortho,error=error)
    CALL cp_fm_cholesky_decompose(ortho,error=error)

    ! 1) Construct a new H matrix, which has the current C0 as eigenvectors,
    !    possibly shifted by an amount lambda,
    !    and the same spectrum as the original H matrix in the space orthogonal to the C0
    !    with P=C0 C0 ^ T
    !    (1 - PS)^T H (1-PS) + (PS)^T (H - lambda S ) (PS)
    !    we exploit that the C0 are already the ritz states of H
    CALL cp_fm_create(matrix_sc0,matrix_c0%matrix_struct,name="sc0",error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_s,matrix_c0,matrix_sc0,k,error=error)
    CALL cp_fm_create(matrix_hc0,matrix_c0%matrix_struct,name="hc0",error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

       ! An aside, try to estimate the error on the ritz values, we'll need it later on
       CALL cp_fm_create(matrix_shc0,matrix_c0%matrix_struct,name="shc0",error=error)
       CALL cp_fm_cholesky_restore(matrix_hc0,k,ortho,matrix_shc0,"SOLVE",transa="T",error=error)
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
       CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       ! since we only use diagonal elements this is a bit of a waste
       CALL cp_gemm('T','N',k,k,n,1.0_dp,matrix_shc0,matrix_shc0,0.0_dp,matrix_s1,error=error)
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
       CALL copy_dbcsr_to_fm(matrix_h,matrix_tmp,error=error)
       CALL cp_fm_upper_to_full(matrix_tmp,matrix_pre,error=error)
    ! tmp = H ( 1 - PS )
    CALL cp_gemm('N','T',n,n,k,-1.0_dp,matrix_hc0,matrix_sc0,1.0_dp,matrix_tmp,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_left,fm_struct_tmp,name="matrix_left",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_gemm('T','N',k,n,n,1.0_dp,matrix_c0,matrix_tmp,0.0_dp,matrix_left,error=error)
    ! tmp = (1 - PS)^T H (1-PS)
    CALL cp_gemm('N','N',n,n,k,-1.0_dp,matrix_sc0,matrix_left,1.0_dp,matrix_tmp,error=error)
    CALL cp_fm_release(matrix_left,error=error)

    ALLOCATE(shifted_evals(k))
    lambda = lambda_base + error_estimate
    shifted_evals=c0_evals - lambda
    CALL cp_fm_to_fm(matrix_sc0,matrix_hc0,error=error)
    CALL cp_fm_column_scale(matrix_hc0,shifted_evals)
    CALL cp_gemm('N','T',n,n,k,1.0_dp,matrix_hc0,matrix_sc0,1.0_dp,matrix_tmp,error=error)

    ! 2) diagonalize this operator
    CALL cp_fm_cholesky_reduce(matrix_tmp,ortho,error=error)
    CALL choose_eigv_solver(matrix_tmp,matrix_pre,preconditioner_env%full_evals,error=error)
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
        CALL cp_gemm('T','N',k,k,n,1.0_dp,matrix_sc0,matrix_tmp,0.0_dp,matrix_s1,error=error)
        CALL choose_eigv_solver(matrix_s1,matrix_s2,norms,error=error)
        WRITE(*,*) "matrix norm deviation (should be close to zero): ", MAXVAL(ABS(ABS(norms)-1.0_dp))
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
  CALL timestop(handle)

END SUBROUTINE make_full_all

!
! the corresponding apply_full_all uses the decomposed form to apply the preconditioner
!

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, n, &
                                                ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_fm_type), POINTER                :: matrix_tmp

  CALL timeset(routineN,handle)

  CALL cp_fm_get_info(matrix_in,nrow_global=n,ncol_global=k,error=error)

  CALL cp_fm_create(matrix_tmp,matrix_in%matrix_struct,name="apply_full_all",error=error)
  CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data,error=error)

  !
  CALL cp_gemm('T','N',n,k,n,1.0_dp,preconditioner_env%fm, &
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
  CALL cp_gemm('N','N',n,k,n,1.0_dp,preconditioner_env%fm, &
                  matrix_tmp,0.0_dp,matrix_out,error=error)

  CALL cp_fm_release(matrix_tmp,error=error)

  CALL timestop(handle)

END SUBROUTINE apply_full_all

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_all(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_offset, col_size, &
                                                handle, i, j, row, &
                                                row_offset, row_size
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: DATA
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

  CALL timeset(routineN,handle)

  CALL cp_dbcsr_init(matrix_tmp,error=error)
  CALL cp_dbcsr_copy(matrix_tmp,matrix_in,name="apply_full_all",error=error)
  CALL cp_dbcsr_multiply('T','N',1.0_dp,preconditioner_env%dbcsr_matrix, &
                  matrix_in,0.0_dp,matrix_tmp,error=error)
  ! do the right scaling
  CALL cp_dbcsr_iterator_start(iter, matrix_tmp)
  DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
     CALL cp_dbcsr_iterator_next_block(iter, row, col, DATA, &
          row_size=row_size, col_size=col_size, &
          row_offset=row_offset, col_offset=col_offset)
     DO j=1,col_size
     DO i=1,row_size
        dum=1.0_dp/MAX(preconditioner_env%energy_gap, &
             preconditioner_env%full_evals( row_offset+i-1 )&
             -preconditioner_env%occ_evals( col_offset+j-1 ))
        DATA(i,j)=DATA(i,j)*dum
     ENDDO
     ENDDO
  ENDDO
  CALL cp_dbcsr_iterator_stop(iter)
  ! mult back
  CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%dbcsr_matrix, &
                  matrix_tmp,0.0_dp,matrix_out,error=error)
  CALL cp_dbcsr_release(matrix_tmp, error=error)
  CALL timestop(handle)

END SUBROUTINE apply_all

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param fm ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param energy_homo ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
SUBROUTINE make_full_single(preconditioner_env, fm, matrix_h, matrix_s, &
                       energy_homo, energy_gap , error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_s

! ---
! ---

  CALL timeset(routineN,handle)

  NULLIFY(fm_h,fm_s,fm_struct_tmp,evals)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm,error=error)
  ENDIF
  CALL cp_dbcsr_get_info(matrix_h,nfullrows_total=n)
  ALLOCATE(evals(n))

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h",error=error)
  CALL cp_fm_create(fm_s,fm_struct_tmp, name="fm_s",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  CALL copy_dbcsr_to_fm(matrix_h,fm_h,error=error)
  CALL copy_dbcsr_to_fm(matrix_s,fm_s,error=error)
  CALL cp_fm_cholesky_decompose(fm_s,error=error)
  CALL cp_fm_cholesky_reduce(fm_h,fm_s,error=error)
  CALL choose_eigv_solver(fm_h,fm,evals,error=error)
  CALL cp_fm_cholesky_restore(fm,n,fm_s,fm_h,"SOLVE",error=error)
  DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
  ENDDO
  CALL cp_fm_to_fm(fm_h,fm,error=error)
  CALL cp_fm_column_scale(fm,evals)
  CALL cp_gemm('N','T',n,n,n,1.0_dp,fm,fm_h,0.0_dp,fm_s,error=error)
  CALL cp_fm_to_fm(fm_s,fm,error=error)

  DEALLOCATE(evals)
  CALL cp_fm_release(fm_h,error=error)
  CALL cp_fm_release(fm_s,error=error)

  CALL timestop(handle)

END SUBROUTINE make_full_single

! different types of preconditioner come here
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_s_inverse(preconditioner_env, matrix_h, matrix_s, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_s_inverse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp

  CALL timeset(routineN,handle)

  NULLIFY(fm_struct_tmp)

  IF (ASSOCIATED(preconditioner_env%fm)) THEN
     CALL cp_fm_release(preconditioner_env%fm,error=error)
  ENDIF
  CALL cp_dbcsr_get_info(matrix_h,nfullrows_total=n)
  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(preconditioner_env%fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  CALL copy_dbcsr_to_fm(matrix_s,preconditioner_env%fm,error=error)

  CALL timestop(handle)

END SUBROUTINE make_full_s_inverse

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param fm ...
!> \param matrix_t ...
!> \param matrix_s ...
!> \param energy_gap ...
!> \param mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_kinetic(preconditioner_env, fm, matrix_t, matrix_s, &
                               energy_gap, mixed_precision, error)
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), POINTER             :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_gap
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: shift
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp

  failure = .FALSE.
  CALL timeset(routineN,handle)

  CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)

  NULLIFY(fm_struct_tmp)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm,error=error)
  ENDIF
  CALL cp_dbcsr_get_info(matrix_t,nfullrows_total=n)
  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner", use_sp=mixed_precision, error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  shift=MAX(0.0_dp,energy_gap)
  CALL cp_dbcsr_add(matrix_t,matrix_s,alpha_scalar=1.0_dp,beta_scalar=shift,error=error)
  CALL copy_dbcsr_to_fm(matrix_t,fm,error=error)
  CALL timestop(handle)

END SUBROUTINE make_full_kinetic

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param fm ...
!> \param matrix_s ...
!> \param mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_inverse_cholesky(preconditioner_env, fm, matrix_s, mixed_precision, &
       error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_inverse_cholesky', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_fm_type), POINTER                :: fm_work

    failure = .FALSE.

    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)

    NULLIFY(fm_work)

    CALL cp_fm_create(fm_work,fm%matrix_struct,name="fm_work",use_sp=mixed_precision,error=error)
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
          CALL copy_dbcsr_to_fm(matrix_s,fm,error=error)
          CALL cp_fm_cholesky_decompose(fm,error=error)
       ELSE
          CALL cp_fm_set_all(fm,alpha=0._dp,beta=1._dp,error=error)
       ENDIF
    ENDIF
    CALL cp_fm_cholesky_invert(fm,error=error)

    CALL cp_fm_upper_to_full(fm,fm_work,error=error)
    CALL cp_fm_release(fm_work,error=error)

    CALL timestop(handle)

  END SUBROUTINE make_full_inverse_cholesky

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param fm ...
!> \param matrix_s ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_fact_cholesky(preconditioner_env, fm, matrix_s, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_fact_cholesky', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error

    failure = .FALSE.

    CALL timeset(routineN,handle)

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
          CALL copy_dbcsr_to_fm(matrix_s,fm,error=error)
          CALL cp_fm_cholesky_decompose(fm,error=error)
       ELSE
          CALL cp_fm_set_all(fm,alpha=0._dp,beta=1._dp,error=error)
       ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE make_full_fact_cholesky

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param fm ...
!> \param matrix_h ...
!> \param energy_homo ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
SUBROUTINE make_full_single_ortho(preconditioner_env, fm, matrix_h, &
                       energy_homo, energy_gap , error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_single_ortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_h, fm_s

  CALL timeset(routineN,handle)
  NULLIFY(fm_h,fm_s,fm_struct_tmp,evals)

  IF (ASSOCIATED(fm)) THEN
     CALL cp_fm_release(fm,error=error)
  ENDIF
  CALL cp_dbcsr_get_info(matrix_h,nfullrows_total=n)
  ALLOCATE(evals(n))

  CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
  CALL cp_fm_create(fm,fm_struct_tmp, name="preconditioner",error=error)
  CALL cp_fm_create(fm_h,fm_struct_tmp, name="fm_h",error=error)
  CALL cp_fm_create(fm_s,fm_struct_tmp, name="fm_s",error=error)
  CALL cp_fm_struct_release(fm_struct_tmp,error=error)

  CALL copy_dbcsr_to_fm(matrix_h,fm_h,error=error)
  
  CALL choose_eigv_solver(fm_h,fm,evals,error=error)
  DO i=1,n
        evals(i)=1.0_dp/MAX(evals(i)-energy_homo,energy_gap)
  ENDDO
  CALL cp_fm_to_fm(fm,fm_h,error=error)
  CALL cp_fm_column_scale(fm,evals)
  CALL cp_gemm('N','T',n,n,n,1.0_dp,fm,fm_h,0.0_dp,fm_s,error=error)
  CALL cp_fm_to_fm(fm_s,fm,error=error)

  DEALLOCATE(evals)
  CALL cp_fm_release(fm_h,error=error)
  CALL cp_fm_release(fm_s,error=error)

  CALL timestop(handle)

END SUBROUTINE make_full_single_ortho
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_c0 ...
!> \param matrix_h ...
!> \param c0_evals ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
SUBROUTINE make_full_single_inverse_ortho(preconditioner_env, matrix_c0, matrix_h, c0_evals, energy_gap,error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'make_full_single_inverse_ortho', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: eval_shift = 5.0_dp , &
                                                fudge_factor = 2.0_dp

    INTEGER                                  :: handle, k, n
    REAL(KIND=dp)                            :: error_estimate, &
                                                preconditioner_shift
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_shift, diag, evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_hc0, matrix_s1, &
                                                matrix_sc0, matrix_tmp, &
                                                matrix_tmp2

! arbitrary upshift of the occupied evals
! fudge factor for taking the error estimate into account

  CALL timeset(routineN,handle)

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
    CALL cp_dbcsr_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    ! since we only use diagonal elements this is a bit of a waste
    CALL cp_gemm('T','N',k,k,n,1.0_dp,matrix_hc0,matrix_hc0,0.0_dp,matrix_s1,error=error)
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
    CALL cp_gemm('N','T',n,n,k,1.0_dp,matrix_sc0,matrix_sc0,0.0_dp,preconditioner_env%fm,error=error)
    CALL cp_fm_release(matrix_sc0,error=error)
    DEALLOCATE(c0_shift)

    ! get H added to the shift
    CALL copy_dbcsr_to_fm(matrix_h,matrix_tmp,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,preconditioner_env%fm,1.0_dp,matrix_tmp,error=error)

    ! preconditioner shift, we target the middle of the occupied spectrum, and taking into account the error_estimate
    ! write(*,*) "Error estimate = ",error_estimate
    preconditioner_shift=-(MINVAL(c0_evals)+ MAXVAL(c0_evals))/2.0_dp + &
                           error_estimate*fudge_factor
    CALL cp_fm_set_all(matrix_tmp,alpha=0._dp,beta=1._dp,error=error)
    CALL cp_fm_scale_and_add(1.0_dp,preconditioner_env%fm,preconditioner_shift,matrix_tmp,error=error)
    ! check evals
    IF (.FALSE.) THEN
       CALL cp_fm_to_fm(preconditioner_env%fm,matrix_tmp,error=error)
       CALL cp_fm_create(matrix_tmp2,matrix_tmp%matrix_struct,name="matrix_tmp2",error=error)
       ALLOCATE(evals(n))
       CALL choose_eigv_solver(matrix_tmp,matrix_tmp2,evals,error=error)
       CALL cp_fm_release(matrix_tmp2,error=error)
       WRITE(*,*) "evals ",evals
       DEALLOCATE(evals)
    ENDIF

    CALL cp_fm_release(matrix_tmp,error=error)

  CALL timestop(handle)

END SUBROUTINE make_full_single_inverse_ortho
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_c0 ...
!> \param matrix_h ...
!> \param c0_evals ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
SUBROUTINE make_full_all_ortho(preconditioner_env, matrix_c0, matrix_h, c0_evals, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c0_evals
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_all_ortho', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: fudge_factor = 0.25_dp, &
                                                lambda_base = 10.0_dp

    INTEGER                                  :: handle, k, n
    REAL(KIND=dp)                            :: error_estimate, lambda
    REAL(KIND=dp), DIMENSION(:), POINTER     :: diag, norms, shifted_evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_hc0, matrix_left, &
                                                matrix_pre, matrix_s1, &
                                                matrix_s2, matrix_sc0, &
                                                matrix_tmp

  CALL timeset(routineN,handle)

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
    CALL cp_dbcsr_sm_fm_multiply(matrix_h,matrix_c0,matrix_hc0,k,error=error)

       ! An aside, try to estimate the error on the ritz values, we'll need it later on
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=k, &
                                context=preconditioner_env%ctxt, &
                                para_env=preconditioner_env%para_env,error=error)
       CALL cp_fm_create(matrix_s1,fm_struct_tmp,name="matrix_s1",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       ! since we only use diagonal elements this is a bit of a waste
       CALL cp_gemm('T','N',k,k,n,1.0_dp,matrix_hc0,matrix_hc0,0.0_dp,matrix_s1,error=error)
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

    CALL copy_dbcsr_to_fm(matrix_h,matrix_tmp,error=error)
    CALL cp_fm_upper_to_full(matrix_tmp,matrix_pre,error=error)
    ! tmp = H ( 1 - PS )
    CALL cp_gemm('N','T',n,n,k,-1.0_dp,matrix_hc0,matrix_sc0,1.0_dp,matrix_tmp,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=k,ncol_global=n, &
                             context=preconditioner_env%ctxt, &
                             para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(matrix_left,fm_struct_tmp,name="matrix_left",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_gemm('T','N',k,n,n,1.0_dp,matrix_c0,matrix_tmp,0.0_dp,matrix_left,error=error)
    ! tmp = (1 - PS)^T H (1-PS)
    CALL cp_gemm('N','N',n,n,k,-1.0_dp,matrix_sc0,matrix_left,1.0_dp,matrix_tmp,error=error)
    CALL cp_fm_release(matrix_left,error=error)

    ALLOCATE(shifted_evals(k))
    lambda = lambda_base + error_estimate
    shifted_evals=c0_evals - lambda
    CALL cp_fm_to_fm(matrix_sc0,matrix_hc0,error=error)
    CALL cp_fm_column_scale(matrix_hc0,shifted_evals)
    CALL cp_gemm('N','T',n,n,k,1.0_dp,matrix_hc0,matrix_sc0,1.0_dp,matrix_tmp,error=error)

    ! 2) diagonalize this operator
     CALL choose_eigv_solver(matrix_tmp,matrix_pre,preconditioner_env%full_evals,error=error)


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
        CALL cp_gemm('T','N',k,k,n,1.0_dp,matrix_sc0,matrix_tmp,0.0_dp,matrix_s1,error=error)
        CALL choose_eigv_solver(matrix_s1,matrix_s2,norms,error=error)

        WRITE(*,*) "matrix norm deviation (should be close to zero): ", MAXVAL(ABS(ABS(norms)-1.0_dp))
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

  CALL timestop(handle)


END SUBROUTINE make_full_all_ortho


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param preconditioner ...
!> \param prec_type ...
!> \param nspins ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE restart_preconditioner(qs_env,preconditioner,prec_type,nspins,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(preconditioner_p_type), &
      DIMENSION(:), POINTER                  :: preconditioner
    INTEGER, INTENT(IN)                      :: prec_type, nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restart_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: failure

    IF (ASSOCIATED(preconditioner)) THEN
       SELECT CASE(prec_type)
       CASE(ot_precond_full_all,ot_precond_full_single, ot_precond_full_single_inverse) ! these depend on the ks matrix
         DO ispin=1,SIZE(preconditioner)
            CALL destroy_preconditioner(preconditioner(ispin)%preconditioner,error=error)
            DEALLOCATE(preconditioner(ispin)%preconditioner)
         ENDDO
         DEALLOCATE(preconditioner,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CASE(ot_precond_none,ot_precond_full_kinetic,ot_precond_s_inverse) ! these are 'independent'
         ! do nothing
       CASE DEFAULT
           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF

    ! add an OT preconditioner if none is present
    IF (.NOT.ASSOCIATED(preconditioner)) THEN
         SELECT CASE(prec_type)
         CASE(ot_precond_full_all,ot_precond_full_single_inverse)
            ALLOCATE(preconditioner(nspins), stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         CASE DEFAULT
            ALLOCATE(preconditioner(1), stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         END SELECT
         DO ispin=1,SIZE(preconditioner)
            ALLOCATE(preconditioner(ispin)%preconditioner)
            CALL init_preconditioner(preconditioner(ispin)%preconditioner,&
                                     para_env=qs_env%para_env,&
                                     blacs_env=qs_env%blacs_env,error=error)
         ENDDO
    END IF

  END SUBROUTINE restart_preconditioner

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param ot_preconditioner ...
!> \param prec_type ...
!> \param solver_type ...
!> \param energy_gap ...
!> \param nspins ...
!> \param has_unit_metric ...
!> \param mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE prepare_preconditioner(qs_env,mos,matrix_ks,matrix_s,&
                                    ot_preconditioner,prec_type,solver_type,&
                                    energy_gap,nspins,has_unit_metric,mixed_precision,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(preconditioner_p_type), &
      DIMENSION(:), POINTER                  :: ot_preconditioner
    INTEGER, INTENT(IN)                      :: prec_type, solver_type
    REAL(dp), INTENT(IN)                     :: energy_gap
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(IN), OPTIONAL            :: has_unit_metric, &
                                                mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_preconditioner', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: msg
    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: do_co_rotate, failure, &
                                                my_has_unit_metric, &
                                                my_mixed_precision, &
                                                use_mo_coeff_b
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: kinetic
    TYPE(cp_dbcsr_type), POINTER             :: matrix_t, mo_coeff_b
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    failure = .FALSE.
    my_has_unit_metric = .FALSE.
    IF(PRESENT(has_unit_metric)) my_has_unit_metric = has_unit_metric
    my_mixed_precision = .FALSE.
    IF(PRESENT(mixed_precision)) my_mixed_precision = mixed_precision

    NULLIFY(matrix_t, mo_coeff_b, mo_coeff, kinetic)

    IF(qs_env%dft_control%qs_control%semi_empirical .OR. qs_env%dft_control%qs_control%dftb) THEN
      IF(prec_type==ot_precond_full_kinetic) THEN
          msg="Full_kinetic not available for semi-empirical methods"
          CPErrorMessage(cp_failure_level,routineP,TRIM(msg),error)
      END IF
      matrix_t => matrix_s(1)%matrix
    ELSE
      CPPrecondition(.NOT. my_has_unit_metric,cp_failure_level,routineP,error,failure)
      CALL get_qs_env(qs_env, kinetic=kinetic, error=error)
      matrix_t => kinetic(1)%matrix
    END IF

    SELECT CASE(prec_type)
        CASE(ot_precond_none)
           DO ispin = 1,SIZE(ot_preconditioner)
             ot_preconditioner(ispin)%preconditioner%in_use=0
           END DO
        CASE(ot_precond_full_all,ot_precond_full_single_inverse)
           do_co_rotate = ASSOCIATED(qs_env%mo_derivs)

           DO ispin=1,nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff_b=mo_coeff_b,mo_coeff=mo_coeff)
             use_mo_coeff_b =mos(ispin)%mo_set%use_mo_coeff_b
             IF(use_mo_coeff_b .AND. do_co_rotate) THEN
               CALL calculate_subspace_eigenvalues(mo_coeff_b,matrix_ks(ispin)%matrix,&
                    do_rotation = .TRUE., &
                    co_rotate=qs_env%mo_derivs(ispin)%matrix,&
                    para_env=qs_env%para_env,&
                    blacs_env=qs_env%blacs_env,error=error)
             ELSEIF(use_mo_coeff_b) THEN
               CALL calculate_subspace_eigenvalues(mo_coeff_b,matrix_ks(ispin)%matrix,&
                    do_rotation = .TRUE., &
                    para_env=qs_env%para_env,&
                    blacs_env=qs_env%blacs_env,error=error)
             ELSE
               CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                    do_rotation = .TRUE., error=error)
             END IF
             IF(my_has_unit_metric) THEN
               CALL make_preconditioner(ot_preconditioner(ispin)%preconditioner, &
                          prec_type, &
                          solver_type, &
                          matrix_h=matrix_ks(ispin)%matrix, &
                          mo_set=mos(ispin)%mo_set, &
                          energy_gap=energy_gap, &
                          convert_precond_to_dbcsr=.TRUE.,&
                          error=error)
             ELSE
               CALL make_preconditioner(ot_preconditioner(ispin)%preconditioner, &
                        prec_type, &
                        solver_type, &
                        matrix_h=matrix_ks(ispin)%matrix,&
                        matrix_s=matrix_s(1)%matrix,&
                        matrix_t=matrix_t, &
                        mo_set=mos(ispin)%mo_set,&
                        energy_gap=energy_gap,&
                        mixed_precision=my_mixed_precision,&
                        convert_precond_to_dbcsr=.TRUE.,&
                        error=error)
             END IF
           ENDDO

        CASE DEFAULT
          IF(my_has_unit_metric) THEN
            CALL make_preconditioner(ot_preconditioner(1)%preconditioner, &
                     prec_type, &
                     solver_type, &
                     matrix_h=matrix_ks(1)%matrix,&
                     mo_set=mos(1)%mo_set,&
                     energy_gap=energy_gap,&
                     convert_precond_to_dbcsr=.TRUE.,&
                     error=error)
          ELSE
            CALL make_preconditioner(ot_preconditioner(1)%preconditioner, &
                     prec_type, &
                     solver_type, &
                     matrix_h=matrix_ks(1)%matrix,&
                     matrix_s=matrix_s(1)%matrix,&
                     matrix_t=matrix_t, &
                     mo_set=mos(1)%mo_set,&
                     energy_gap=energy_gap,&
                     mixed_precision=my_mixed_precision,&
                     convert_precond_to_dbcsr=.TRUE.,&
                     error=error)
          END IF

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE prepare_preconditioner

END MODULE preconditioner

