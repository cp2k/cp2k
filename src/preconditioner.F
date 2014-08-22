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
       array_data, array_i1d_obj, block_add_on_diag, block_chol_inv, &
       cp_create_bl_distribution, cp_dbcsr_add, cp_dbcsr_add_on_diag, &
       cp_dbcsr_btriu, cp_dbcsr_checksum, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_deallocate_matrix, &
       cp_dbcsr_distribution, cp_dbcsr_distribution_release, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_get_block, cp_dbcsr_get_block_p, &
       cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, cp_dbcsr_get_info, &
       cp_dbcsr_get_matrix_type, cp_dbcsr_get_num_blocks, &
       cp_dbcsr_get_occupation, cp_dbcsr_get_stored_coordinates, &
       cp_dbcsr_init, cp_dbcsr_init_p, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_nblkcols_local, cp_dbcsr_nblkrows_local, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_p_type, &
       cp_dbcsr_print, cp_dbcsr_put_block, cp_dbcsr_release, &
       cp_dbcsr_release_p, cp_dbcsr_row_block_sizes, cp_dbcsr_trace, &
       cp_dbcsr_type, cp_dbcsr_verify_matrix, cp_dbcsr_work_create, &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_error_type, dbcsr_mp_group, &
       dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
       dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_obj, &
       dbcsr_norm_frobenius, dbcsr_scale_mat, dbcsr_type_no_symmetry, &
       dbcsr_type_real_4, dbcsr_type_real_default
  USE cp_dbcsr_operations,             ONLY: &
       copy_dbcsr_to_fm, copy_fm_to_dbcsr, cp_dbcsr_copy_vec, &
       cp_dbcsr_from_fm, cp_dbcsr_multiply_vec, cp_dbcsr_pack_vec, &
       cp_dbcsr_sm_fm_multiply, cp_dbcsr_unpack_vec, packed_vec_bcast, &
       packed_vec_bif_tech2, packed_vec_build_u, packed_vec_ini, &
       packed_vec_scale
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace,&
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
       ot_precond_solver_inv_chol, ot_precond_solver_sainv, &
       ot_precond_sparse_diag, ot_precond_sparse_kinetic, &
       ot_precond_sparse_kinetic_sainv
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size,&
                                             sp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             preconditioner_p_type,&
                                             preconditioner_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'

  PUBLIC :: make_preconditioner, restart_preconditioner
  PUBLIC :: apply_preconditioner, prepare_preconditioner

  PRIVATE :: make_sparse_diag, make_full_single, make_local_block

  INTERFACE apply_preconditioner
     MODULE PROCEDURE apply_preconditioner_dbcsr
     MODULE PROCEDURE apply_preconditioner_fm
  END INTERFACE

  INTERFACE apply_solve_lin_system
     MODULE PROCEDURE apply_solve_lin_system_dbcsr
     MODULE PROCEDURE apply_solve_lin_system_fm
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
    CASE (ot_precond_sparse_diag,ot_precond_full_single)
       needs_homo=.TRUE.
       ! XXXX to be removed if homo estimate only is implemented
       needs_full_spectrum=.TRUE.
    CASE (ot_precond_full_kinetic,ot_precond_s_inverse,ot_precond_sparse_kinetic)
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
  CASE (ot_precond_sparse_diag)
     preconditioner_env%in_use=ot_precond_sparse_diag
     IF (my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Only PRECOND_SOLVER DEFAULT for the moment")
     END IF
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
  CASE(ot_precond_sparse_kinetic)
     preconditioner_env%in_use=ot_precond_sparse_kinetic
     IF(my_solver_type.NE.ot_precond_solver_default) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Only PRECOND_SOLVER DEFAULT for the moment")
     ENDIF
     CALL make_sparse_kinetic(preconditioner_env, matrix_t, matrix_s, energy_gap ,error=error)
     CALL make_diag_inner_precond(preconditioner_env, preconditioner_env%sparse_matrix,&
                                  preconditioner_env%sparse_matrix_inner, error=error)
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
  CASE (ot_precond_solver_sainv)
     !
     preconditioner_env%solver=ot_precond_solver_sainv
     CALL make_sparse_inverse_bif(preconditioner_env, preconditioner_env%fm, matrix_s, &
                                  error=error)
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
  CASE (ot_precond_sparse_diag)
     CALL apply_sparse_diag(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_single)
     CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_all)
     CALL apply_full_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE(ot_precond_sparse_kinetic)
     CALL apply_solve_lin_system(preconditioner_env, matrix_in, matrix_out, error=error)
  CASE (ot_precond_sparse_kinetic_sainv)
     SELECT CASE (preconditioner_env%solver)
     CASE (ot_precond_solver_sainv)
        CALL apply_sparse_single(preconditioner_env, matrix_in, matrix_out,error=error)
        !CALL apply_full_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,"Solver not implemented")
     END SELECT
  CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
     SELECT CASE (preconditioner_env%solver)
     CASE(ot_precond_solver_inv_chol,ot_precond_solver_sainv)
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
  CASE (ot_precond_sparse_diag,ot_precond_full_single)
     CALL apply_single(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE (ot_precond_full_all)
     CALL apply_all(preconditioner_env, matrix_in, matrix_out,error=error)
  CASE(ot_precond_sparse_kinetic)
     CALL apply_solve_lin_system(preconditioner_env, matrix_in, matrix_out, error=error)
  CASE (ot_precond_sparse_kinetic_sainv)
     SELECT CASE (preconditioner_env%solver)
     CASE (ot_precond_solver_sainv)
        CALL apply_single(preconditioner_env, matrix_in, matrix_out,error=error)
     CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,"Wrong solver")
     END SELECT
  CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse,ot_precond_s_inverse)
     SELECT CASE (preconditioner_env%solver)
     CASE(ot_precond_solver_inv_chol,ot_precond_solver_sainv)
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
! a sparse block diagonal approximation
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
SUBROUTINE apply_sparse_diag(preconditioner_env, matrix_in, matrix_out, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: k

  CALL cp_fm_get_info(matrix_in,ncol_global=k,error=error)
  CALL cp_dbcsr_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_in, &
                         matrix_out,k,error=error)

END SUBROUTINE apply_sparse_diag

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param energy_homo ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
SUBROUTINE make_sparse_diag(preconditioner_env, matrix_h, matrix_s, &
                          energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sparse_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iblock_col, &
                                                iblock_row, n, nblocks
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s
    TYPE(cp_dbcsr_iterator)                  :: iter

  CALL timeset(routineN,handle)

    IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       CALL cp_dbcsr_deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
       NULLIFY(preconditioner_env%sparse_matrix)
    ENDIF

    ALLOCATE(preconditioner_env%sparse_matrix)
    CALL cp_dbcsr_init(preconditioner_env%sparse_matrix,error=error)
    CALL cp_dbcsr_create(preconditioner_env%sparse_matrix, ' PRECONDITIONER ', &
         cp_dbcsr_distribution (matrix_h),&
         cp_dbcsr_get_matrix_type (matrix_h), cp_dbcsr_row_block_sizes(matrix_h),&
         cp_dbcsr_col_block_sizes(matrix_h), 0, 0, error=error)

  CALL cp_dbcsr_get_info(matrix_h,nfullrows_total=nblocks)

  CALL cp_dbcsr_iterator_start(iter, matrix_h)
  DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
     CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_h, blk)
     CALL cp_dbcsr_get_block_p(matrix=matrix_s,&
          row=iblock_row,col=iblock_col,BLOCK=block_s,found=found)

     IF(.NOT.ASSOCIATED(block_s))CALL stop_program(routineN,moduleN,__LINE__,&
                                                   ".NOT.ASSOCIATED(block_s)")

     IF (iblock_col .EQ. iblock_row) THEN
        n=SIZE(block_s,1)

        ALLOCATE(block_pre(n,n))

        CALL make_local_block(block_h,block_s,block_pre,  &
                                                  energy_homo,energy_gap,3)

        CALL cp_dbcsr_put_block(matrix=preconditioner_env%sparse_matrix,&
                             row=iblock_row,&
                             col=iblock_col,&
                             block=block_pre)

        DEALLOCATE(block_pre)

     ENDIF
  ENDDO

  CALL cp_dbcsr_iterator_stop(iter)
  CALL cp_dbcsr_finalize(preconditioner_env%sparse_matrix,error=error)

  CALL timestop(handle)

END SUBROUTINE make_sparse_diag
! *****************************************************************************
!> \brief ...
!> \param block_h ...
!> \param block_s ...
!> \param block_pre ...
!> \param energy_homo ...
!> \param energy_gap ...
!> \param TYPE ...
! *****************************************************************************
SUBROUTINE make_local_block(block_h,block_s,block_pre, &
                                 energy_homo,energy_gap, TYPE)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_s, block_pre
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    INTEGER                                  :: TYPE

    CHARACTER(len=*), PARAMETER :: routineN = 'make_local_block', &
      routineP = moduleN//':'//routineN

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
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                   "work",dp_size*lwork)
  ALLOCATE(iwork(liwork),STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                   "iwork",int_size*liwork)
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
     IF (info.ne.0)  CALL stop_program(routineN,moduleN,__LINE__,"Error in dpotrf")
     CALL DSYGST(1,'U',n,block_evec(1,1),n,block_chol(1,1),n,info)
     IF (info.ne.0)  CALL stop_program(routineN,moduleN,__LINE__,"Error in dsygst")
     CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                 iwork(1), liwork, info)
     IF (info.NE.0) CALL stop_program(routineN,moduleN,__LINE__,"Error in dsyevd")
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
     IF (info.ne.0) CALL stop_program(routineN,moduleN,__LINE__,"Error in dsyevd S")
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
     IF (info.ne.0) CALL stop_program(routineN,moduleN,__LINE__,"Error in dsyevd H")
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
     IF (info.ne.0) CALL stop_program(routineN,moduleN,__LINE__,"Error in dsyevd H")
     block_pre(:,:)=0.0_dp
     DO i=1,n
          block_pre(i,i)=1.0_dp/MAX(evals(i),energy_gap)
     ENDDO
     CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                                                 0.0_dp,block_buf1(1,1),n)
     CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                                                  0.0_dp,block_pre(1,1),n)

  END SELECT

  DEALLOCATE (iwork,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iwork")
  DEALLOCATE (work,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")
  DEALLOCATE (block_chol,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"block_chol")
  DEALLOCATE (block_evec,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"block_evec")
  DEALLOCATE (block_buf1,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"block_buf1")
  DEALLOCATE (evals,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"evals")

END SUBROUTINE make_local_block
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
! Preconditioners for orthogonal basis sets
! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_h ...
!> \param energy_homo ...
!> \param energy_gap ...
!> \param error ...
! *****************************************************************************
SUBROUTINE make_sparse_diag_ortho(preconditioner_env, matrix_h, energy_homo, energy_gap, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    REAL(KIND=dp)                            :: energy_homo, energy_gap
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sparse_diag_ortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iblock_col, &
                                                iblock_row, n, nblocks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_h, block_pre, block_s
    TYPE(cp_dbcsr_iterator)                  :: iter

  CALL timeset(routineN,handle)

    IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       CALL cp_dbcsr_deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
       NULLIFY(preconditioner_env%sparse_matrix)
    ENDIF

    ALLOCATE(preconditioner_env%sparse_matrix)
    CALL cp_dbcsr_init(preconditioner_env%sparse_matrix,error=error)
    CALL cp_dbcsr_create(preconditioner_env%sparse_matrix, ' PRECONDITIONER ', &
         cp_dbcsr_distribution (matrix_h), cp_dbcsr_get_matrix_type (matrix_h),&
         cp_dbcsr_row_block_sizes(matrix_h),&
         cp_dbcsr_col_block_sizes(matrix_h), 0, 0, error=error)

  CALL cp_dbcsr_get_info(matrix_h,nfullrows_total=nblocks)

  CALL cp_dbcsr_iterator_start(iter, matrix_h)
  DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
     CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_h, blk)

     IF (iblock_col .EQ. iblock_row) THEN
        n=SIZE(block_h,1)
        ALLOCATE(block_pre(n,n),block_s(n,n))
        block_s(:,:)=0._dp
        DO i=1,n
           block_s(i,i)=1._dp
        END DO

        CALL make_local_block(block_h,block_s,block_pre,energy_homo,energy_gap,3)

        CALL cp_dbcsr_put_block(matrix=preconditioner_env%sparse_matrix,&
                             row=iblock_row,&
                             col=iblock_col,&
                             BLOCK=block_pre)

        DEALLOCATE(block_pre,block_s)

     ENDIF
  ENDDO
  CALL cp_dbcsr_iterator_stop(iter)
  CALL cp_dbcsr_finalize(preconditioner_env%sparse_matrix,error=error)

  CALL timestop(handle)

END SUBROUTINE make_sparse_diag_ortho

! *****************************************************************************

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

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param sm ...
!> \param sm_inner ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_diag_inner_precond(preconditioner_env, sm, sm_inner, error)
    !
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: sm, sm_inner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_diag_inner_precond', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iblock_col, &
                                                iblock_row, info, istat, &
                                                liwork, lwork, n, nblocks
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: block_evec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_buf1, block_h, block_pre
    TYPE(cp_dbcsr_iterator)                  :: iter

  CALL timeset(routineN,handle)
    failure = .FALSE.

    CPPrecondition(ASSOCIATED(sm),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(sm_inner)) THEN
       CALL cp_dbcsr_deallocate_matrix(sm_inner,error=error)
       NULLIFY(sm_inner)
    ENDIF

    ALLOCATE(sm_inner)
    CALL cp_dbcsr_init(sm_inner,error=error)
    CALL cp_dbcsr_create(sm_inner, ' PRECONDITIONER ', &
         cp_dbcsr_distribution (sm), cp_dbcsr_get_matrix_type (sm),&
         cp_dbcsr_row_block_sizes(sm),&
         cp_dbcsr_col_block_sizes(sm), 0, 0, error=error)

    CALL cp_dbcsr_get_info(sm,nfullrows_total=nblocks)

    CALL cp_dbcsr_iterator_start(iter, sm)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_h, blk)

       IF (iblock_col .EQ. iblock_row) THEN
          n = SIZE(block_h,1)
          ALLOCATE(block_pre(n,n),block_evec(n,n),block_buf1(n,n),evals(n))
          lwork=1+6*n+2*n**2+50
          liwork=5*n+3
          ALLOCATE (work(lwork),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "work",dp_size*lwork)
          ALLOCATE (iwork(liwork),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "iwork",int_size*liwork)
          !
          ! A simple block diagonal preconditoner
          block_evec(:,:) = block_h(:,:)
          CALL DSYEVD('V','U', n, block_evec(1,1), n, evals(1), work(1), lwork, &
                         iwork(1), liwork, info)
             IF (info.ne.0) CALL stop_program(routineN,moduleN,__LINE__,"DSYEVD H")
             block_pre(:,:)=0.0_dp
             DO i=1,n
                IF (evals(i).lt.1e-10_dp) CALL stop_program(routineN,moduleN,__LINE__,&
                                                            "evals(i).lt.1e-10_dp")
                block_pre(i,i)=1.0_dp/evals(i)
             ENDDO
             CALL DGEMM('N','N',n,n,n,1.0_dp,block_evec(1,1),n,block_pre(1,1),n, &
                        0.0_dp,block_buf1(1,1),n)
             CALL DGEMM('N','T',n,n,n,1.0_dp,block_buf1(1,1),n,block_evec(1,1),n, &
                        0.0_dp,block_pre(1,1),n)
             !
             ! add the block
             CALL cp_dbcsr_put_block(matrix=sm_inner,&
                                  row=iblock_row,&
                                  col=iblock_col,&
                                  BLOCK=block_pre)

             DEALLOCATE(block_pre,block_evec,block_buf1,evals,work,iwork)

          ENDIF
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
       CALL cp_dbcsr_finalize(sm_inner,error=error)

  CALL timestop(handle)

  END SUBROUTINE make_diag_inner_precond

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE apply_solve_lin_system_fm(preconditioner_env, matrix_in, matrix_out, error)
    !
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_solve_lin_system_fm', &
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
    output_unit= cp_logger_get_default_io_unit(logger)
    !
    ! residual
    CALL cp_fm_trace(matrix_in,matrix_in,residual_ot,error=error)
    residual_ot = SQRT(residual_ot/(REAL(n,dp)*REAL(m,dp)))
    !
    ! r = b-A*x
    CALL cp_dbcsr_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_out,r,m,error=error)
    CALL cp_fm_scale_and_add(-1.0_dp,r,1.0_dp,matrix_in,error=error)
    CALL cp_fm_trace(r,r,tr_rr,error=error)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',0,': ||r_ls||=',&
              SQRT(tr_rr/(REAL(n,dp)*REAL(m,dp))),' ||r_ot||=',residual_ot
    ENDIF
    !
    ! let's go!
    DO iter = 1,max_iter
       !
       ! s = M * r
       IF(ASSOCIATED(preconditioner_env%sparse_matrix_inner)) THEN
          CALL cp_dbcsr_sm_fm_multiply(preconditioner_env%sparse_matrix_inner,r,s,m,error=error)
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
       CALL cp_dbcsr_sm_fm_multiply(preconditioner_env%sparse_matrix,p,Ap,m,error=error)
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
          CALL cp_dbcsr_sm_fm_multiply(preconditioner_env%sparse_matrix,matrix_out,r,m,error=error)
          CALL cp_fm_scale_and_add(-1.0_dp,r,1.0_dp,matrix_in,error=error)
       ELSE
          CALL cp_fm_scale_and_add(1.0_dp,r,-alpha,Ap,error=error)
       ENDIF
       tr_sr_old = tr_sr_new
       !
       ! printing
       IF(MOD(iter,10).EQ.0.AND.output_unit>0) THEN
          WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',iter,': ||r_ls||=',&
                 SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp))),' ||r_ls||/||r_ot||=',&
                 SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot
       ENDIF
       !
       ! exit
       IF(SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot.LT.eps) EXIT
       !
    ENDDO
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',iter,': ||r_ls||=',&
              SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp))),' ||r_ls||/||r_ot||=',&
              SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot
    ENDIF
    IF(iter.GT.max_iter) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Increase max_iter")
    ENDIF
    !
    CALL cp_fm_release(r ,error=error)
    CALL cp_fm_release(p ,error=error)
    CALL cp_fm_release(s ,error=error)
    CALL cp_fm_release(Ap,error=error)
  END SUBROUTINE apply_solve_lin_system_fm

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE apply_solve_lin_system_dbcsr(preconditioner_env, matrix_in, matrix_out, error)
    !
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type)                      :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_solve_lin_system_dbcsr', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_iter = 200, &
                                                restart_every = 50
    REAL(dp), PARAMETER                      :: eps = 1.0e-2_dp

    INTEGER                                  :: iter, m, n, output_unit
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, beta, residual_ot, &
                                                tr_pAp, tr_rr, tr_sr_new, &
                                                tr_sr_old
    TYPE(cp_dbcsr_type)                      :: Ap, p, r, s
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(preconditioner_env%sparse_matrix),cp_failure_level,routineP,error,failure)
    !
    CALL cp_dbcsr_get_info(matrix_in,nfullrows_total=n,nfullcols_total=m)
    CALL cp_dbcsr_init(r,error=error)
    CALL cp_dbcsr_init(p,error=error)
    CALL cp_dbcsr_init(s,error=error)
    CALL cp_dbcsr_init(Ap,error=error)
    CALL cp_dbcsr_copy(r,matrix_in,error=error)
    CALL cp_dbcsr_copy(p,matrix_in,error=error)
    CALL cp_dbcsr_copy(s,matrix_in,error=error)
    CALL cp_dbcsr_copy(Ap,matrix_in,error=error)
    !
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    !
    ! residual
    CALL cp_dbcsr_trace(matrix_in,matrix_in,residual_ot,error=error)
    residual_ot = SQRT(residual_ot/(REAL(n,dp)*REAL(m,dp)))
    !
    ! r = b-A*x
    CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%sparse_matrix,matrix_out,&
         0.0_dp,r,error=error)
    CALL cp_dbcsr_add(r,matrix_in,-1.0_dp,1.0_dp,error=error)
    CALL cp_dbcsr_trace(r,r,tr_rr,error=error)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',0,': ||r_ls||=',&
              SQRT(tr_rr/(REAL(n,dp)*REAL(m,dp))),' ||r_ot||=',residual_ot
    ENDIF
    !
    ! let's go!
    DO iter = 1,max_iter
       !
       ! s = M * r
       IF(ASSOCIATED(preconditioner_env%sparse_matrix_inner)) THEN
          CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%sparse_matrix_inner,r,&
               0.0_dp,s,error=error)
       ELSE
          CALL cp_dbcsr_copy(s,r,error=error)
       ENDIF
       CALL cp_dbcsr_trace(s,r,tr_sr_new,error=error)
       !
       ! beta (this might be better when O(N) beta=(s_k,r_k-r_{k-1})/(s_{k-1},r_{k-1}))
       IF(iter.EQ.1) THEN
          beta = 0.0_dp
       ELSE
          beta = tr_sr_new/tr_sr_old
       ENDIF
       !
       ! p = r + beta * p
       CALL cp_dbcsr_add(p,s,beta,1.0_dp,error=error)
       !
       ! Ap = A * p
       CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%sparse_matrix,p,&
            0.0_dp,Ap,error=error)
       CALL cp_dbcsr_trace(Ap,p,tr_pAp,error=error)
       !
       ! alpha = Tr(s_k'*r_k)/Tr(p_k'*Ap_k)
       alpha = tr_sr_new/tr_pAp
       !
       ! x = x + alpha * p
       CALL cp_dbcsr_add(matrix_out,p,1.0_dp,alpha,error=error)
       !
       ! r = r - alpha * Ap or r = b - A * x
       IF(MOD(iter,restart_every).EQ.0) THEN
          CALL cp_dbcsr_multiply('N','N',1.0_dp,preconditioner_env%sparse_matrix,matrix_out,&
               0.0_dp,r,error=error)
          CALL cp_dbcsr_add(r,matrix_in,-1.0_dp,1.0_dp,error=error)
       ELSE
          CALL cp_dbcsr_add(r,Ap,1.0_dp,-alpha,error=error)
       ENDIF
       tr_sr_old = tr_sr_new
       !
       ! printing
       IF(MOD(iter,10).EQ.0.AND.output_unit>0) THEN
          WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',iter,': ||r_ls||=',&
                 SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp))),' ||r_ls||/||r_ot||=',&
                 SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot
       ENDIF
       !
       ! exit
       IF(SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot.LT.eps) EXIT
       !
    ENDDO
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,I3,2(A,E12.5))') 'apply_solve_lin_system:',iter,': ||r_ls||=',&
              SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp))),' ||r_ls||/||r_ot||=',&
              SQRT(tr_sr_new/(REAL(n,dp)*REAL(m,dp)))/residual_ot
    ENDIF
    IF(iter.GT.max_iter) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"increase max_iter")
    ENDIF
    !
    CALL cp_dbcsr_release(r, error=error)
    CALL cp_dbcsr_release(p, error=error)
    CALL cp_dbcsr_release(s, error=error)
    CALL cp_dbcsr_release(Ap, error=error)
  END SUBROUTINE apply_solve_lin_system_dbcsr

! *****************************************************************************
!> \brief DBCSR BIF based Inverse with
!> \param[in,out] preconditioner_env    The preconditioner environment
!> \param         fm                    in input contains the matrix to
!>                                      approximatly inverse
!> \param         matrix_s              The overlap matrix
!> \param[in,out] error       cp2k error
!> \par History
!>      - Created by VW
!> \author Valery Weber
! *****************************************************************************
  SUBROUTINE make_sparse_inverse_bif(preconditioner_env, fm, matrix_s, error)
    !
    TYPE(preconditioner_type), INTENT(INOUT) :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sparse_inverse_bif', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: drop_v = 1.0E-5_dp, &
                                                thresh_inner = 1.0E-2_dp, &
                                                thresh_post = 1.0E-5_dp

    INTEGER :: d_offset, handle, istat, k, k_s, max_blk_size, &
      max_row_blk_size, max_vec_size, mp_group, mynode, mypcol, myprow, n, &
      nblkcols_local, nblkrows_local, nblkrows_total, npcols, nprows, &
      numnodes, proc_holds_diag_blk
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: dbg_print_chksum, &
                                                dbg_print_matrix, failure, &
                                                found, tr
    REAL(dp)                                 :: chksum, occ_in, occ_out, s, &
                                                t(100), trace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: blk, blk_p, pkd_d, pkd_inv_d, &
                                                pkd_u, pkd_u_fac, pkd_v, &
                                                pkd_v_fac, pkd_v_tmp
    TYPE(cp_dbcsr_type)                      :: mat_li, mat_li_scaled, &
                                                mat_lt, mat_lt_scaled, mat_t, &
                                                mat_tmp, mat_v
    TYPE(cp_dbcsr_type), POINTER             :: prcnd
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_work1, fm_work2
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    t=0.0_dp
    dbg_print_matrix = .FALSE.
    dbg_print_chksum = .TRUE.
!

    failure = .FALSE.

    CALL timeset(routineN,handle)

    CALL cpu_time(t(1))

    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)

    !
    ! some infos
    CALL cp_fm_get_info(fm,nrow_global=n,error=error)
    !
    ! allocate
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n,ncol_global=n,&
         context=preconditioner_env%ctxt, &
         para_env=preconditioner_env%para_env,error=error)
    CALL cp_fm_create(fm_work1,fm_struct_tmp, name="preconditioner_work1",error=error)
    CALL cp_fm_create(fm_work2,fm_struct_tmp, name="preconditioner_work2",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    !
    ! here fm should contain the matrix to inverse
    CALL cp_fm_to_fm(fm,fm_work1,error=error)
    CALL cp_fm_upper_to_full(fm_work1,fm_work2,error=error)
    !
    !
    ! copy matrices
    CALL cp_dbcsr_init (mat_t, error)
    CALL cp_dbcsr_create(mat_t, 'mat t', &
         cp_dbcsr_distribution(matrix_s), dbcsr_type_no_symmetry, cp_dbcsr_row_block_sizes(matrix_s),&
         cp_dbcsr_col_block_sizes(matrix_s), 0, 0,&
         cp_dbcsr_get_data_type(matrix_s), error=error)
    CALL copy_fm_to_dbcsr(fm_work1,mat_t, error=error)
    IF(.TRUE.)CALL cp_dbcsr_verify_matrix(mat_t, error=error)
    IF(dbg_print_matrix)CALL cp_dbcsr_print(mat_t,error=error)

    occ_in = cp_dbcsr_get_occupation(mat_t)
    CALL cp_dbcsr_filter(mat_t,1e-5_dp,error=error)
    IF(.FALSE.)CALL cp_dbcsr_verify_matrix(mat_t, error=error)
    occ_out = cp_dbcsr_get_occupation(mat_t)

    WRITE(*,*) routineN//' occ_in',occ_in,' occ_out',occ_out

    CALL cp_dbcsr_init (mat_li,error=error)
    CALL cp_dbcsr_create(mat_li, 'matrix li', &
         cp_dbcsr_distribution(mat_t), dbcsr_type_no_symmetry, cp_dbcsr_row_block_sizes(mat_t),&
         cp_dbcsr_col_block_sizes(mat_t), 0, 0,&
         cp_dbcsr_get_data_type(mat_t), error=error)
    CALL cp_dbcsr_init (mat_v,error=error)
    CALL cp_dbcsr_create(mat_v, 'matrix v', &
         cp_dbcsr_distribution(mat_t), dbcsr_type_no_symmetry, cp_dbcsr_row_block_sizes(mat_t),&
         cp_dbcsr_col_block_sizes(mat_t), 0, 0,&
         cp_dbcsr_get_data_type(mat_t), error=error)

    CALL cp_dbcsr_work_create(mat_li,work_mutable=.TRUE.,&
         nblks_guess=cp_dbcsr_get_num_blocks(mat_t),sizedata_guess=cp_dbcsr_get_data_size(mat_t),&
         error=error)
    CALL cp_dbcsr_work_create(mat_v,work_mutable=.TRUE.,&
         nblks_guess=cp_dbcsr_get_num_blocks(mat_t),sizedata_guess=cp_dbcsr_get_data_size(mat_t),&
         error=error)
    !
    !
    prcnd => preconditioner_env%sparse_matrix
    IF (ASSOCIATED(prcnd)) THEN
       CALL cp_dbcsr_deallocate_matrix(prcnd,error=error)
       NULLIFY(prcnd)
    ENDIF
    ALLOCATE(prcnd)
    CALL cp_dbcsr_init(prcnd,error=error)
    CALL cp_dbcsr_create(prcnd, ' PRECONDITIONER ', &
         cp_dbcsr_distribution(mat_t), dbcsr_type_no_symmetry, cp_dbcsr_row_block_sizes(mat_t),&
         cp_dbcsr_col_block_sizes(mat_t), cp_dbcsr_get_num_blocks(mat_t), cp_dbcsr_get_data_size(mat_t),&
         cp_dbcsr_get_data_type(mat_t), error=error)

    ! some variables
    s = 1.0_dp ! s > 0

    nblkrows_total = cp_dbcsr_nblkrows_total(mat_t)
    nblkrows_local = cp_dbcsr_nblkrows_local(mat_t)
    nblkcols_local = cp_dbcsr_nblkcols_local(mat_t)
    row_blk_size => array_data(cp_dbcsr_row_block_sizes(mat_t))
    mp_obj = dbcsr_distribution_mp (cp_dbcsr_distribution (mat_t))
    mp_group = dbcsr_mp_group (mp_obj)

    max_row_blk_size = MAXVAL(row_blk_size)
    max_blk_size = max_row_blk_size**2
    max_vec_size = 100 * MAX(nblkrows_local,nblkcols_local)*max_blk_size+nblkrows_total+1 ! need to fix that
    CALL mp_max(max_vec_size,mp_group)

    WRITE(*,*) 'nblkrows_local',nblkrows_local
    WRITE(*,*) 'max_blk_size',max_blk_size
    WRITE(*,*) 'nblkrows_local',nblkrows_local
    WRITE(*,*) 'nblkcols_local',nblkcols_local
    WRITE(*,*) 'nblkrows_total',nblkrows_total
    WRITE(*,*) 'max_vec_size',max_vec_size

    ALLOCATE(pkd_u(max_vec_size), pkd_v(max_vec_size), pkd_v_fac(max_vec_size),&
             pkd_u_fac(max_vec_size), pkd_d(max_vec_size), pkd_inv_d(max_vec_size), &
             blk(max_blk_size), pkd_v_tmp(max_vec_size), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    pkd_u(:) = 0.0_dp
    pkd_v(:) = 0.0_dp
    pkd_d(:) = 0.0_dp
    pkd_v_fac(:) = 0.0_dp
    pkd_u_fac(:) = 0.0_dp
    pkd_inv_d(:) = 0.0_dp
    blk(:) = 0.0_dp

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    !blacs2mpi => dbcsr_mp_pgrid (mp_obj)

    WRITE(*,*) 'numnodes',numnodes
    WRITE(*,*) 'mynode',mynode
    WRITE(*,*) 'myprow',myprow
    WRITE(*,*) 'mypcol',mypcol
    WRITE(*,*) 'npcols',npcols
    WRITE(*,*) 'nprows',nprows

    !
    !CALL dbcsr_get_stored_coordinates (matrix, r, c, tr, p)
    !p .EQ. dbcsr_mp_mynode (mp_obj)
    !

    !
    CALL cpu_time(t(2))
    WRITE(*,*) 't1',t(2)-t(1)
    !
    ! let's bif !
    d_offset = 1
    DO k = 1,nblkrows_total
       !write(*,*) 'enter loop k=',k
       tr = .FALSE.
       k_s = k
       CALL cp_dbcsr_get_stored_coordinates (mat_t, k_s, k_s, tr, proc_holds_diag_blk)
       !write(*,*) 'proc_holds_diag_blk',proc_holds_diag_blk, proc_holds_diag_blk .EQ. dbcsr_mp_mynode (mp_obj)

       CALL cpu_time(t(3))
       !------------------------------------------------------------------------
       ! set v and u
       !
       !v(k,1:n)=A(k,1:n);
       CALL cp_dbcsr_copy_vec(mat_v, mat_t, 'row', k, error=error)
       !chksum = cp_dbcsr_checksum(mat_v,error)
       !write(*,*) 'checksum1 checksum(mat_v)',chksum
       CALL cpu_time(t(4))
       t(50) = t(50) + t(4)-t(3)
       !
       !v(k,k)=v(k,k)-s;
       IF(proc_holds_diag_blk .EQ. dbcsr_mp_mynode (mp_obj)) THEN
          CALL cp_dbcsr_get_block_p(mat_v, k, k, blk_p, found)
          IF (.NOT.found) CALL stop_program(routineN,moduleN,__LINE__,&
                                            "Matrix block not found")
          CALL block_add_on_diag(row_blk_size(k), blk_p, -s)
       ENDIF

       !chksum = cp_dbcsr_checksum(mat_v,error)
       !write(*,*) 'checksum2 checksum(mat_v)',chksum
       !
       !u(k,k)=1;
       !CALL block_set_d(row_blk_size(k), row_blk_size(k), blk, 1.0_dp, 0.0_dp)
       !CALL cp_dbcsr_put_block(mat_li, blk, k, k, error=error)!PARA local put
       !old CALL block_set(row_blk_size(k), row_blk_size(k), blk, 1.0_dp, 0.0_dp)
       !old CALL cp_dbcsr_put_block(mat_li, blk, k, k, error=error)!PARA local put
       !
       !------------------------------------------------------------------------
       ! build the factor inv(d)*u*A
       !
       CALL cpu_time(t(5))
       t(51) = t(51) + t(5)-t(4)
       ! compute u(1:k-1,1:n) * A(1:n,k);
       CALL cp_dbcsr_multiply_vec(mat_li, mat_t, 1, k-1, k, pkd_v_fac, error=error)!PARA pkd_v_fac is hold by
       !                                                                        !PARA the same procs as A(1:n,k)
       !write(*,*) 'pkd_v_fac:',pkd_v_fac(1:int(pkd_v_fac(nblkrows_total+1))-1)
       !chksum = sum(abs(   pkd_v_fac  ( nblkrows_total+2: int(pkd_v_fac (nblkrows_total+1))-1) ))
       !write(*,*) 'checksum3.0',chksum
       !
       CALL cpu_time(t(6))
       t(52) = t(52) + t(6)-t(5)
       !
       ! bcast pkd_v_fac
       !CALL packed_vec_bcast(pkd_v_fac, source, 'columnwise', mp_obj, error=error) !PARA replicate columnwise
       !
       ! compute inv(d(i,i)) * uAk(i,k) / s , i=1,k-1
       !chksum = sum(abs(  pkd_inv_d ))
       !write(*,*) 'checksum3.1',chksum
       CALL packed_vec_scale(1.0_dp/s, pkd_inv_d, pkd_v_fac, k-1, k, row_blk_size, 'left', error=error)!PARA local
       CALL cpu_time(t(7))
       t(53) = t(53) + t(7)-t(6)
       !------------------------------------------------------------------------
       ! build the factor inv(d)*v
       !
       ! pack v(:,k)
       !old CALL cp_dbcsr_pack_vec(mat_v, k, pkd_u_fac, 'column', error=error)
       !
       ! bcast pkd_u_fac
       !CALL packed_vec_bcast(pkd_u_fac, source, 'columnwise', mp_obj, error=error) !PARA replicate columnwise
       !
       ! compute inv(d(i,i)) * v(i,k) / s , i=1,k-1
       !old CALL packed_vec_scale(1.0_dp/s, pkd_inv_d, pkd_u_fac, k-1, k, row_blk_size, 'left', error=error)!PARA local
       !------------------------------------------------------------------------
       ! update v and u
       !
       ! pack u(k,:) and v(k,:)
       CALL packed_vec_ini(pkd_v, nblkrows_total, error)

       !CALL cp_dbcsr_pack_vec(mat_v, k, pkd_v, 'row', error=error)  !PARA local
       CALL cpu_time(t(8))
       t(54) = t(54) + t(8)-t(7)
       !old CALL cp_dbcsr_pack_vec(mat_li, k, pkd_u, 'row', error=error) !PARA local
       !
       ! bcast the u and v
       !call packed_vec_bcast(pkd_v, -1, 'all', .false., nblkrows_total, row_blk_size(k)*row_blk_size, mp_obj, error)
       !chksum = sum(abs(pkd_v(nblkrows_total+2:int(pkd_v(nblkrows_total+1))-1)))
       !write(*,*) 'checksum3.2 pkd_v',chksum
       !chksum = sum(abs(pkd_v_fac(nblkrows_total+2:int(pkd_v_fac(nblkrows_total+1))-1)))
       !write(*,*) 'checksum3.2 pkd_v_fac',chksum


       !CALL packed_vec_bcast(pkd_v, source, 'columnwise', mp_obj, error=error)
       !CALL packed_vec_bcast(pkd_u, source, 'columnwise', mp_obj, error=error)
       !
       ! the core
       !v(k,1:n) = v(k,1:n) - uAk(i,k)' * inv(d(i,i)) * v(i,1:n) / s; i = 1,k-1
       !u(k,1:n) = u(k,1:n) -   v(i,k)' * inv(d(i,i)) * u(i,1:n) / s; i = 1,k-1
       !old CALL packed_vec_bif_tech(mat_v, mat_li, pkd_v_fac, pkd_u_fac, k, pkd_v, pkd_u, error=error)!PARA local
       CALL packed_vec_bif_tech2(mat_v, pkd_v_fac, k, pkd_v, error=error)!PARA local
       CALL cpu_time(t(9))
       t(55) = t(55) + t(9)-t(8)
       !chksum = cp_dbcsr_checksum(mat_v,error)
       !write(*,*) 'checksum4.0 mat_v',chksum
       ! bcast the v


       !write(*,'(A,1000F12.6)') '1pkd_v',pkd_v(nblkrows_total+2:int(pkd_v(nblkrows_total+1))-1)
       CALL packed_vec_bcast(pkd_v, -1, 'all', .TRUE. ,nblkrows_total, &
            row_blk_size(k)*row_blk_size, mp_obj, error)


       CALL cp_dbcsr_unpack_vec(mat_v, k, pkd_v, 'row', do_sum=.TRUE., error=error)  !PARA local
       CALL cp_dbcsr_pack_vec(mat_v, k, pkd_v, 'row', error=error)  !PARA local
       CALL packed_vec_bcast(pkd_v, -1, 'all', .FALSE. ,nblkrows_total, &
            row_blk_size(k)*row_blk_size, mp_obj, error)

       !write(*,*) '2pkd_v',pkd_v( 1 : int(pkd_v(nblkrows_total+1)) - 1 )

       !chksum = sum(abs( pkd_v ( nblkrows_total+2 : int(pkd_v(nblkrows_total+1))-1 ) ))
       !write(*,'(A,1000F12.6)') '2pkd_v',pkd_v(nblkrows_total+2:int(pkd_v(nblkrows_total+1))-1)
       !write(*,*) 'checksum4.1',chksum
       !
       ! compute norm
       !
       ! filter the packed v
       !
       ! build u(k,:) = [-v(k,1:k-1)/s 1 0 ... 0]
       CALL packed_vec_build_u(pkd_u, pkd_v, k, nblkrows_total, s, row_blk_size, error=error)
       CALL cpu_time(t(10))
       t(56) = t(56) + t(10)-t(9)
       !
       ! reduce the u and v
       !CALL packed_vec_reduce(pkd_v, to, 'columnwise', vec_blk_size, n, mp_obj, error=error)
       !CALL packed_vec_reduce(pkd_u, to, 'columnwise', vec_blk_size, n, mp_obj, error=error)
       !
       ! unpack u(k,:) and v(k,:)
       CALL cp_dbcsr_unpack_vec(mat_v, k, pkd_v, 'row', error=error)  !PARA local
       CALL cp_dbcsr_unpack_vec(mat_li, k, pkd_u, 'row', error=error) !PARA local
       !chksum = cp_dbcsr_checksum(mat_v,error)
       !write(*,*) 'checksum5 checksum(mat_v)',chksum
       !chksum = cp_dbcsr_checksum(mat_li,error)
       !write(*,*) 'checksum6.0 checksum(mat_li)',chksum
       !write(*,*) 'checksum6.1',sum(abs(pkd_u))
       !write(*,*) 'pkd_u',pkd_u
       CALL cpu_time(t(11))
       t(57) = t(57) + t(11)-t(10)
       !------------------------------------------------------------------------
       ! build d and inv(d)
       !
       ! d(k,k) = v(k,k) / s + eye(k);
       IF(proc_holds_diag_blk .EQ. dbcsr_mp_mynode (mp_obj)) THEN
          CALL cp_dbcsr_get_block(mat_v, k, k, blk, found)
          IF (.NOT.found) CALL stop_program(routineN,moduleN,__LINE__,&
                                            "Matrix block not found" )
          CALL dscal(row_blk_size(k)**2, 1.0_dp/s, blk, 1)
          CALL block_add_on_diag(row_blk_size(k), blk, 1.0_dp)
          CALL dcopy(row_blk_size(k)**2, blk, 1, pkd_d(d_offset), 1)
          !
          ! inv(d(k,k))
          CALL block_chol_inv(row_blk_size(k), blk)
          CALL dcopy(row_blk_size(k)**2, blk, 1, pkd_inv_d(d_offset), 1)
       ENDIF
       !
       ! bcast d(k,k) and inv(d(k,k))
       CALL mp_bcast(pkd_d(d_offset:d_offset+row_blk_size(k)**2-1), &
            proc_holds_diag_blk, dbcsr_mp_group (mp_obj))
       CALL mp_bcast(pkd_inv_d(d_offset:d_offset+row_blk_size(k)**2-1), &
            proc_holds_diag_blk, dbcsr_mp_group (mp_obj))
       !
       ! update the block pointer
       d_offset = d_offset + row_blk_size(k)**2
       !------------------------------------------------------------------------
       CALL cpu_time(t(12))
       t(58) = t(58) + t(12)-t(11)
    ENDDO
    CALL cp_dbcsr_finalize(mat_v, error=error)
    CALL cp_dbcsr_finalize(mat_li, error=error)

    !WRITE(*,*) 'tin1',t(50)
    !WRITE(*,*) 'tin2',t(51)
    !WRITE(*,*) 'tin3',t(52)
    !WRITE(*,*) 'tin4',t(53)
    !WRITE(*,*) 'tin5',t(54)
    !WRITE(*,*) 'tin6',t(55)
    !WRITE(*,*) 'tin7',t(56)
    !WRITE(*,*) 'tin8',t(57)
    !WRITE(*,*) 'tin9',t(58)
    !
    ! D=d;
    ! Linv=u;
    ! Lt=inv(D)*triu(v+s*I);
    CALL cp_dbcsr_init (mat_lt,error=error)
    CALL cp_dbcsr_btriu(mat_lt, mat_v, error=error)
    !chksum = cp_dbcsr_checksum(mat_lt,error)
    !WRITE(*,*) 'checksum(lt)_1',chksum
    CALL cp_dbcsr_add_on_diag(mat_lt, s, error=error)
    !chksum = cp_dbcsr_checksum(mat_lt,error)
    !WRITE(*,*) 'checksum(lt)_2',chksum,SUM(ABS(pkd_inv_d))
    CALL dbcsr_scale_mat(mat_lt%matrix, alpha_matrix=pkd_inv_d, side='left', error=dbcsr_error)
    !chksum = cp_dbcsr_checksum(mat_lt,error)
    !WRITE(*,*) 'checksum(lt)_3',chksum
    CALL cp_dbcsr_init (mat_li_scaled,error=error)
    CALL cp_dbcsr_copy(mat_li_scaled, mat_li, error=error)
    CALL dbcsr_scale_mat(mat_li_scaled%matrix, alpha_matrix=pkd_inv_d, side='left', error=dbcsr_error)

    CALL cp_dbcsr_init (mat_lt_scaled,error=error)
    CALL cp_dbcsr_copy(mat_lt_scaled, mat_lt, error=error)
    CALL dbcsr_scale_mat(mat_lt_scaled%matrix, alpha_matrix=pkd_d, side='left', error=dbcsr_error)

    occ_in = cp_dbcsr_get_occupation(mat_lt)
    occ_out = cp_dbcsr_get_occupation(mat_li)
    WRITE(*,*) routineN//' occ(Lt)',occ_in,' occ(Li)',occ_out

    chksum = cp_dbcsr_checksum(mat_v,error=error)
    IF(dbg_print_chksum)WRITE(*,*) 'checksum(v)',chksum
    chksum = cp_dbcsr_checksum(mat_li,error=error)
    IF(dbg_print_chksum)WRITE(*,*) 'checksum(li)',chksum
    chksum = cp_dbcsr_checksum(mat_li_scaled,error=error)
    IF(dbg_print_chksum)WRITE(*,*) 'checksum(li_scaled)',chksum
    chksum = cp_dbcsr_checksum(mat_lt,error=error)
    IF(dbg_print_chksum)WRITE(*,*) 'checksum(lt)',chksum


    IF(dbg_print_matrix)CALL cp_dbcsr_print(mat_v,error=error)
    IF(dbg_print_matrix)CALL cp_dbcsr_print(mat_li,error=error,matlab_format=.TRUE.)
    IF(dbg_print_matrix)CALL cp_dbcsr_print(mat_li_scaled,error=error,matlab_format=.TRUE.)
    IF(dbg_print_matrix)CALL cp_dbcsr_print(mat_lt,error=error)

    !
    ! At this point we should have inv(L), L^t, D and inv(D)
    ! compute the precond or apply the factors
    !
    ! C = A*B'
    IF(.FALSE.)CALL cp_dbcsr_verify_matrix(prcnd, error=error)
    IF(.FALSE.)CALL cp_dbcsr_verify_matrix(mat_li, error=error)
    IF(.FALSE.)CALL cp_dbcsr_verify_matrix(mat_t, error=error)
    IF(.FALSE.)CALL cp_dbcsr_verify_matrix(mat_li_scaled, error=error)
    CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mat_li, mat_li_scaled, 0.0_dp, prcnd, error=error)
    !CALL cp_dbcsr_print(prcnd, error, matlab_format=.TRUE.)

    chksum = cp_dbcsr_checksum(prcnd,error=error)
    IF(dbg_print_chksum)WRITE(*,*) 'checksum(prcnd)',chksum


    !
    ! checks
    IF(.TRUE.) THEN
       CALL cp_dbcsr_trace(mat_lt, mat_li, trace, error=error)
       WRITE(*,*) 'trace(Lt^t*Li)/n-1=',trace / REAL(SUM(row_blk_size),dp) - 1.0_dp

       CALL cp_dbcsr_trace(mat_li, mat_lt, trace, error=error)
       WRITE(*,*) 'trace(Li^t*Lt)/n-1=',trace / REAL(SUM(row_blk_size),dp) - 1.0_dp

       CALL cp_dbcsr_init(mat_tmp, error=error)
       CALL cp_dbcsr_create (mat_tmp, "tmp mat_t", cp_dbcsr_distribution (mat_t), dbcsr_type_no_symmetry,&
            cp_dbcsr_row_block_sizes (mat_t), cp_dbcsr_col_block_sizes (mat_t),&
            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_li_scaled, mat_t, 0.0_dp, mat_tmp, error=error)
       CALL cp_dbcsr_trace(mat_li, mat_tmp, trace, error=error)
       WRITE(*,*) 'trace(Li^t*A*Li*D)/n-1=',trace / REAL(SUM(row_blk_size),dp) - 1.0_dp

       !
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, prcnd, mat_t, 0.0_dp, mat_tmp, error=error)
       CALL cp_dbcsr_norm(mat_tmp, which_norm=dbcsr_norm_frobenius, norm_scalar=trace, error=error)
       WRITE(*,*) 'norm(A*Li*D*li^t)/n-1 = ',trace**2 / REAL(SUM(row_blk_size),dp) - 1.0_dp

       CALL cp_dbcsr_multiply("T", "T", 1.0_dp, mat_t, prcnd, 0.0_dp, mat_tmp, error=error)
       CALL cp_dbcsr_norm(mat_tmp, which_norm=dbcsr_norm_frobenius, norm_scalar=trace, error=error)
       WRITE(*,*) 'norm(A*Li*D*li^t)/n-1 = ',trace**2 / REAL(SUM(row_blk_size),dp) - 1.0_dp

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_t, prcnd, 0.0_dp, mat_tmp, error=error)
       CALL cp_dbcsr_norm(mat_tmp, which_norm=dbcsr_norm_frobenius, norm_scalar=trace, error=error)
       WRITE(*,*) 'norm(Li*D*li^t*A)/n-1 = ',trace**2 / REAL(SUM(row_blk_size),dp) - 1.0_dp
       !
       !norm(Lt'*D*Lt-A,'fro')
       CALL cp_dbcsr_copy(mat_tmp, mat_t, error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mat_lt, mat_lt_scaled, -1.0_dp, mat_tmp, error=error)
       CALL cp_dbcsr_norm(mat_tmp, which_norm=dbcsr_norm_frobenius, norm_scalar=trace, error=error)
       WRITE(*,*) 'norm(Lt^t*D*Lt-A) = ',trace
       !
       !norm(Lt*Linv'-I,'fro')
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, mat_lt, mat_li, 0.0_dp, mat_tmp, error=error)
       CALL cp_dbcsr_norm(mat_tmp, which_norm=dbcsr_norm_frobenius, norm_scalar=trace, error=error)
       WRITE(*,*) 'norm(Lt*Linv^t)/n-1 = ',trace**2 / REAL(SUM(row_blk_size),dp) - 1.0_dp

    ENDIF
    !
    ! cleanup
    DEALLOCATE(pkd_u, pkd_v, pkd_v_fac, pkd_u_fac, pkd_d, pkd_inv_d, pkd_v_tmp, blk, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL cp_dbcsr_release(mat_t, error=error)
    CALL cp_dbcsr_release(mat_v, error=error)
    CALL cp_dbcsr_release(mat_lt, error=error)
    CALL cp_dbcsr_release(mat_li, error=error)
    CALL cp_dbcsr_release(mat_li_scaled, error=error)
    CALL cp_dbcsr_release(mat_lt_scaled, error=error)
    CALL cp_dbcsr_release(mat_tmp, error=error)

    CALL cp_fm_release(fm_work1,error=error)
    CALL cp_fm_release(fm_work2,error=error)

    CALL copy_dbcsr_to_fm (prcnd, fm,error=error)!for the moment
    CALL cp_dbcsr_deallocate_matrix(prcnd,error=error)!for the moment

CALL cpu_time(t(13))

WRITE(*,*) 'tend',t(13)-t(12)

    CALL timestop(handle)

  END SUBROUTINE make_sparse_inverse_bif

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE apply_sparse_single(preconditioner_env, matrix_in, matrix_out,&
       error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_in, matrix_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_sparse_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(cp_dbcsr_type)                      :: matrix_result, matrix_right
    TYPE(dbcsr_distribution_obj)             :: dist

  CALL timeset(routineN,handle)

    WRITE(*,*)' APPLYING SPARSE SINGLE KINETIC'
    CALL cp_fm_get_info(matrix_in, ncol_global=k, error=error)
    CALL cp_create_bl_distribution (col_dist, col_blk_size, k, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution( preconditioner_env%dbcsr_matrix ))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution( preconditioner_env%dbcsr_matrix )),&
         dbcsr_distribution_row_dist(cp_dbcsr_distribution( preconditioner_env%dbcsr_matrix )), col_dist)
    CALL cp_dbcsr_init (matrix_right, error)
    CALL cp_dbcsr_from_fm (matrix_right, matrix_in, 0.0_dp, dist,&
         row_blk_size = cp_dbcsr_row_block_sizes(preconditioner_env%dbcsr_matrix),&
         col_blk_size = col_blk_size, error=error)
    WRITE(*,*)'Pre checksum:',cp_dbcsr_checksum (matrix_right, error=error)
    CALL cp_dbcsr_init (matrix_result, error)
    CALL cp_dbcsr_create(matrix_result, "result", dist, 'N',&
         cp_dbcsr_row_block_sizes (preconditioner_env%dbcsr_matrix),&
         cp_dbcsr_col_block_sizes (matrix_right), error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, preconditioner_env%dbcsr_matrix,&
         matrix_right, 0.0_dp, matrix_result, error=error)
    CALL cp_dbcsr_distribution_release (dist)
    CALL cp_dbcsr_release (matrix_right, error=error)
    WRITE(*,*)'Post checksum:',cp_dbcsr_checksum (matrix_result, error=error)
    CALL cp_dbcsr_release (matrix_result, error=error)
    WRITE(*,*)' DONE'

  CALL timestop(handle)

  END SUBROUTINE apply_sparse_single

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
       CASE(ot_precond_none,ot_precond_full_kinetic,ot_precond_s_inverse, &
             ot_precond_sparse_diag,ot_precond_sparse_kinetic) ! these are 'independent'
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
    TYPE(cp_dbcsr_type), POINTER             :: matrix_t, mo_coeff_b
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    failure = .FALSE.
    my_has_unit_metric = .FALSE.
    IF(PRESENT(has_unit_metric)) my_has_unit_metric = has_unit_metric
    my_mixed_precision = .FALSE.
    IF(PRESENT(mixed_precision)) my_mixed_precision = mixed_precision

    NULLIFY(matrix_t, mo_coeff_b, mo_coeff)

    IF(qs_env%dft_control%qs_control%semi_empirical .OR. qs_env%dft_control%qs_control%dftb) THEN
      IF(prec_type==ot_precond_full_kinetic) THEN
          msg="Full_kinetic not available for semi-empirical methods"
          CPErrorMessage(cp_failure_level,routineP,TRIM(msg),error)
      END IF
      matrix_t => matrix_s(1)%matrix
    ELSE
      CPPrecondition(.NOT. my_has_unit_metric,cp_failure_level,routineP,error,failure)
      matrix_t => qs_env%kinetic(1)%matrix
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

