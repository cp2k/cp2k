!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Higher-level operations on DBCSR matrices.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - Created 2009-05-12
! *****************************************************************************

MODULE dbcsr_operations
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE dbcsr_config,                    ONLY: is_initialized
  USE dbcsr_data_methods,              ONLY: dbcsr_scalar
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution
  USE dbcsr_error_handling
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution, dbcsr_distribution_mp, dbcsr_distribution_new, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data_type, dbcsr_get_info, dbcsr_init, dbcsr_mp_npcols, &
       dbcsr_release
  USE dbcsr_mm_cannon,                 ONLY: dbcsr_mm_cannon_clear_mempools,&
                                             dbcsr_mm_cannon_lib_finalize,&
                                             dbcsr_mm_cannon_lib_init,&
                                             dbcsr_mm_cannon_multiply
  USE dbcsr_operations_low,            ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_block_in_limits, dbcsr_btriu, &
       dbcsr_copy, dbcsr_copy_columns, dbcsr_copy_into_existing, &
       dbcsr_copy_submatrix, dbcsr_filter, dbcsr_frobenius_norm, &
       dbcsr_gershgorin_norm, dbcsr_get_block_diag, dbcsr_get_diag, &
       dbcsr_hadamard_product, dbcsr_init_random, dbcsr_maxabs, dbcsr_norm, &
       dbcsr_scale, dbcsr_scale_by_vector, dbcsr_scale_mat, dbcsr_set, &
       dbcsr_set_diag, dbcsr_sum_replicated, dbcsr_symmetrize_block_diag, &
       dbcsr_trace, dbcsr_tril, dbcsr_triu
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize
  USE kinds,                           ONLY: dp,&
                                             int_1_size,&
                                             int_2_size,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_init_lib, dbcsr_finalize_lib,&
            dbcsr_multiply,&
            dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_scale_mat, dbcsr_add, dbcsr_copy,&
            dbcsr_copy_submatrix, dbcsr_copy_into_existing,&
            dbcsr_get_diag, dbcsr_set_diag, &
            dbcsr_get_block_diag, dbcsr_hadamard_product, &
            dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_btriu, dbcsr_triu, dbcsr_tril,&
            dbcsr_symmetrize_block_diag, dbcsr_copy_columns,&
            dbcsr_init_random, dbcsr_lanczos, dbcsr_block_in_limits,&
            dbcsr_sum_replicated, dbcsr_norm, &
            dbcsr_gershgorin_norm, dbcsr_maxabs, dbcsr_frobenius_norm, &
            dbcsr_clear_mempools

  INTERFACE dbcsr_multiply
     MODULE PROCEDURE dbcsr_mm_cannon_multiply
     MODULE PROCEDURE dbcsr_multiply_s, dbcsr_multiply_d,&
                      dbcsr_multiply_c, dbcsr_multiply_z
  END INTERFACE





 CONTAINS


! *****************************************************************************
!> \brief Initialize the DBCSR library
!>
!> Prepares the DBCSR library for use.
!> \param[in,out] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_init_lib (group, error)
    INTEGER, INTENT(IN)                      :: group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!n_stack_buffers,mem_type, n_stack_mem_regions, stack_size, nstacks,&
!INTEGER, DIMENSION(3)                    :: nxstacks
!   ---------------------------------------------------------------------------
!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

    IF (is_initialized) RETURN
    CALL dbcsr_error_set(routineN, error_handle, error)
    !
    CALL dbcsr_assert (int_1_size, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of an 8-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_2_size, "EQ", 2,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 16-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_4_size, "EQ", 4,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 32-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_8_size, "EQ", 8,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 64-bit integer size!",&
         __LINE__, error=error)

    !$omp parallel default(none)  shared(error)
    CALL dbcsr_mm_cannon_lib_init(error)
    !$omp end parallel

    is_initialized = .TRUE.
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_init_lib


! *****************************************************************************
!> \brief Finalize the DBCSR library
!>
!> Cleans up after the DBCSR library.  Used to deallocate persistent objects.
!> \param[in,out] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize_lib (group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

    IF (.NOT. is_initialized) RETURN
    CALL dbcsr_error_set(routineN, error_handle, error)

    IF(output_unit>0)THEN
       WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=output_unit,FMT="(T2,A,T35,A,T80,A)") "-","DBCSR STATISTICS","-"
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)
    END IF

    !$omp parallel default(none) shared(output_unit, group, error)
    CALL dbcsr_mm_cannon_lib_finalize(group, output_unit, error)
    !$omp end parallel
    IF(output_unit>0) WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)

    is_initialized = .FALSE.
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_finalize_lib


! *****************************************************************************
!> \brief  Deallocate memory contained in mempools
! *****************************************************************************
  SUBROUTINE dbcsr_clear_mempools(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    !$omp parallel default(none) shared( error)
    CALL dbcsr_mm_cannon_clear_mempools(error)
    !$omp end parallel
  END SUBROUTINE dbcsr_clear_mempools

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_s(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=real_4), INTENT(IN)            :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(KIND=real_4), INTENT(IN)            :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_mm_cannon_multiply(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_s


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_d(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=real_8), INTENT(IN)            :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(KIND=real_8), INTENT(IN)            :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_d', &
      routineP = moduleN//':'//routineN

    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND.&
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4 .AND.&
       dbcsr_get_data_type(matrix_c) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_mm_cannon_multiply(transa, transb,&
            dbcsr_scalar(REAL(alpha,real_4)), matrix_a, matrix_b, &
            dbcsr_scalar(REAL(beta,real_4)), matrix_c,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, &
            filter_eps=filter_eps,&
            error=error, flop=flop)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND.&
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .AND.&
           dbcsr_get_data_type(matrix_c) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_mm_cannon_multiply(transa, transb,&
            dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, &
            filter_eps=filter_eps,&
            error=error, flop=flop)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_internal_error,&
            routineP, "This combination of data types NYI",__LINE__, error)
    ENDIF
  END SUBROUTINE dbcsr_multiply_d

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_c(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(KIND=real_4), INTENT(IN)         :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(KIND=real_4), INTENT(IN)         :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_mm_cannon_multiply(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_c


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_z(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(KIND=real_8), INTENT(IN)         :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(KIND=real_8), INTENT(IN)         :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_mm_cannon_multiply(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_z
  
  
! *****************************************************************************
!> \brief compute the extremal eigenvalues of a symmetric real matrix
!>        with a (simple) Lanczos approach
!> \param[in] matrix
!> \param[in] max_iter            maximum iteration
!> \param[in] eps                 convergence parameter
!> \param[out] min_eig, max_eig   approximation to the extremal eigenvalues
!> \param[out] approx_norm_2      approximation of the 2 norm (1 < max_iter < 10)
!> \param[out] converged          true if the iteration converged
!> \param[inout] error
!>
! *****************************************************************************
  SUBROUTINE dbcsr_lanczos(matrix, max_iter, eps, min_eig, max_eig, &
       approx_norm_2, converged, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(dp), INTENT(in), OPTIONAL           :: eps
    REAL(dp), INTENT(out), OPTIONAL          :: min_eig, max_eig, &
                                                approx_norm_2
    LOGICAL, INTENT(OUT), OPTIONAL           :: converged
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_lanczos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, i, info, &
                                                iwork, lwork, n
    LOGICAL                                  :: my_converged
    REAL(dp)                                 :: alpha, beta, max_eig_old, &
                                                min_eig_old, my_eps, &
                                                my_max_eig, my_min_eig, &
                                                nrm_f, nrm_v
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: work
    REAL(dp), DIMENSION(max_iter)            :: evals
    REAL(dp), DIMENSION(max_iter, max_iter)  :: T, T0
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_obj)                          :: f, v, v0

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_get_info(matrix,nfullrows_total=n)

    CALL create_bl_distribution (col_dist, col_blk_size, 1, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(dbcsr_distribution(matrix))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dbcsr_distribution_row_dist(dbcsr_distribution(matrix)), col_dist)

    CALL dbcsr_init(v)
    CALL dbcsr_create(v, 'v', dist, dbcsr_type_no_symmetry, matrix%m%row_blk_size,&
         col_blk_size, 0, 0, error=error)
    CALL dbcsr_finalize(v, error=error)

    CALL dbcsr_init(v0)
    CALL dbcsr_create(v0, 'v0', dist, dbcsr_type_no_symmetry, matrix%m%row_blk_size,&
         col_blk_size, 0, 0, error=error)
    CALL dbcsr_finalize(v0, error=error)

    CALL dbcsr_init(f)
    CALL dbcsr_create(f, 'f', dist, dbcsr_type_no_symmetry, matrix%m%row_blk_size,&
         col_blk_size, 0, 0, error=error)
    CALL dbcsr_finalize(f, error=error)

    CALL dbcsr_distribution_release (dist)
    CALL array_release (col_blk_size)
    CALL array_release (col_dist)

    lwork = 1+2*max_iter+100
    ALLOCATE(work(lwork))

    my_eps = 1.0e-1_dp
    IF(PRESENT(eps)) my_eps = eps

    min_eig_old = 0.0_dp
    max_eig_old = 0.0_dp
    T(:,:) = 0.0_dp
    ! v = rand(n,1)
    CALL dbcsr_init_random(v,error=error)
    ! v = v / norm(v)
    nrm_v = dbcsr_frobenius_norm(v)
    CALL dbcsr_scale(v,1.0_dp/nrm_v,error=error)
    ! f = A * v
    CALL dbcsr_multiply('N','N',1.0_dp,matrix,v,0.0_dp,f,error=error)
    ! alpha = f' * v
    CALL dbcsr_trace(f,v,alpha,error=error)
    ! f = f - alpha * v
    CALL dbcsr_add(f,v,1.0_dp,-alpha,error=error)
    T(1,1) = alpha
    my_min_eig = alpha; my_max_eig = alpha
    my_converged = .FALSE.
    DO i = 2,max_iter
       ! beta = norm(f)
       beta = dbcsr_frobenius_norm(f)
       ! v0 = v
       CALL dbcsr_copy(v0,v,error=error)
       ! v = f / beta
       CALL dbcsr_add(v,f,0.0_dp,1.0_dp/beta,error=error)
       ! f = A * v
       CALL dbcsr_multiply('N','N',1.0_dp,matrix,v,0.0_dp,f,error=error)
       ! f = f - beta * v0
       CALL dbcsr_add(f,v0,1.0_dp,-beta,error=error)
       ! alpha = f' * v
       CALL dbcsr_trace(f,v,alpha,error=error)
       ! f = f - alpha * v
       CALL dbcsr_add(f,v,1.0_dp,-alpha,error=error)
       T(i  ,i-1) = beta
       T(i-1,i  ) = beta
       T(i  ,i  ) = alpha
       !
       max_eig_old = my_max_eig; min_eig_old = my_min_eig
       T0(:,:) = T(:,:)
       CALL DSYEVD('N','U',i,T0(1,1),max_iter,evals(1),work(1),lwork,iwork,1,info)
       CALL dbcsr_assert(info.EQ.0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "DSYEVD", __LINE__, error)
       my_max_eig = MAXVAL(evals(1:i)); my_min_eig = MINVAL(evals(1:i))
       !write(*,*) routinen//' i',i,'max_eig',my_max_eig,' min_eig',my_min_eig
       IF(ABS(my_max_eig-max_eig_old).LT.my_eps.AND.ABS(my_min_eig-min_eig_old).LT.my_eps) THEN
          my_converged = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF(PRESENT(approx_norm_2)) THEN
       ! norm(f,2)
       nrm_f = dbcsr_frobenius_norm(f)
       ! norm(T,2)
       T0(:,:) = T(:,:)
       CALL DSYEVD('N','U',max_iter,T0(1,1),max_iter,evals(1),work(1),lwork,iwork,1,info)
       CALL dbcsr_assert(info.EQ.0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "DSYEVD", __LINE__, error)
       ! norm(T,2)+norm(f,2)
       approx_norm_2 = MAXVAL(ABS(evals)) + nrm_f
    ENDIF

    IF(PRESENT(min_eig)) min_eig = my_min_eig
    IF(PRESENT(max_eig)) max_eig = my_max_eig
    IF(PRESENT(converged)) converged = my_converged

    DEALLOCATE(work)
    CALL dbcsr_release(v)
    CALL dbcsr_release(v0)
    CALL dbcsr_release(f)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_lanczos

END MODULE dbcsr_operations
