!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner_makes
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_arnoldi_ev, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_get_info, cp_dbcsr_init, &
       cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_setup_arnoldi_data, cp_dbcsr_type, &
       cp_set_arnoldi_initial_vector, dbcsr_arnoldi_data, &
       dbcsr_type_symmetric, deallocate_arnoldi_data, get_selected_ritz_val, &
       get_selected_ritz_vec
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_m_by_n_from_template,&
                                             cp_dbcsr_sm_fm_multiply,&
                                             cp_fm_to_dbcsr_row_template
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
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
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE input_constants,                 ONLY: ot_precond_full_all,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_full_single_inverse,&
                                             ot_precond_s_inverse,&
                                             ot_precond_solver_default,&
                                             ot_precond_solver_inv_chol
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_makes'

  PUBLIC :: make_preconditioner_matrix  

CONTAINS


! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param matrix_t ...
!> \param mo_coeff ...
!> \param energy_homo ...
!> \param eigenvalues_ot ...
!> \param energy_gap ...
!> \param solver_type ...
!> \param my_mixed_precision ...
!> \param my_solver_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_preconditioner_matrix(preconditioner_env, matrix_h, matrix_s, matrix_t, mo_coeff,&
                           energy_homo, eigenvalues_ot, energy_gap, solver_type, my_mixed_precision,&
                           my_solver_type, error)
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s, matrix_t
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot
    REAL(KIND=dp)                            :: energy_gap
    INTEGER, INTENT(IN)                      :: solver_type
    LOGICAL                                  :: my_mixed_precision
    INTEGER                                  :: my_solver_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: precon_type
    LOGICAL                                  :: failure

    failure=.FALSE.

    precon_type=preconditioner_env%in_use                
    SELECT CASE (precon_type)
    CASE (ot_precond_full_single)
       IF(my_solver_type.NE.ot_precond_solver_default) &
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Only PRECOND_SOLVER DEFAULT for the moment")
       IF ( PRESENT(matrix_s) ) THEN
          CALL make_full_single(preconditioner_env, preconditioner_env%fm,&
                                matrix_h, matrix_s, energy_homo, energy_gap ,error=error)
       ELSE
          CALL make_full_single_ortho(preconditioner_env, preconditioner_env%fm,&
                                matrix_h, energy_homo, energy_gap ,error=error)
       END IF
  
    CASE (ot_precond_s_inverse)
       IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
       IF (.NOT. PRESENT(matrix_s) ) &
          CALL stop_program(routineN,moduleN,__LINE__, "Type for S=1 not implemented")
       CALL make_full_s_inverse(preconditioner_env,matrix_s,error)
  
    CASE (ot_precond_full_kinetic)
       IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
       IF (.NOT.( PRESENT(matrix_s) .AND. PRESENT(matrix_t) )) &
          CALL stop_program(routineN,moduleN,__LINE__,"Type for S=1 not implemented")
       CALL make_full_kinetic(preconditioner_env, matrix_t, matrix_s, energy_gap, &
                             my_mixed_precision, error=error)
    CASE (ot_precond_full_single_inverse)
       IF(my_solver_type.EQ.ot_precond_solver_default) my_solver_type=ot_precond_solver_inv_chol
       CALL make_full_single_inverse(preconditioner_env, mo_coeff, matrix_h, energy_gap, &
                           matrix_s=matrix_s,error=error)
    CASE (ot_precond_full_all)
       IF(my_solver_type.NE.ot_precond_solver_default) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Only PRECOND_SOLVER DEFAULT for the moment")
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

  END SUBROUTINE make_preconditioner_matrix

! *****************************************************************************
!> \brief Simply takes the overlap matrix as preconditioner
!> \param preconditioner_env ...
!> \param matrix_s ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_s_inverse(preconditioner_env, matrix_s, error)
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_s_inverse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL timeset(routineN,handle)
  
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
  
    IF(.NOT.ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       ALLOCATE(preconditioner_env%sparse_matrix)
       CALL cp_dbcsr_init(preconditioner_env%sparse_matrix,error=error)
    END IF
    CALL cp_dbcsr_copy(preconditioner_env%sparse_matrix,matrix_s,name="full_kinetic",error=error)
  
    CALL timestop(handle)
  
  END SUBROUTINE make_full_s_inverse

! *****************************************************************************
!> \brief kinetic matrix+shift*overlap as preconditioner. Cheap but could
!>        be better
!> \param preconditioner_env ...
!> \param matrix_t ...
!> \param matrix_s ...
!> \param energy_gap ...
!> \param mixed_precision ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_full_kinetic(preconditioner_env, matrix_t, matrix_s, &
                               energy_gap, mixed_precision, error)
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_t, matrix_s
    REAL(KIND=dp)                            :: energy_gap
    LOGICAL, INTENT(IN)                      :: mixed_precision
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: shift

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(matrix_t),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)

    IF(.NOT.ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       ALLOCATE(preconditioner_env%sparse_matrix)
       CALL cp_dbcsr_init(preconditioner_env%sparse_matrix,error=error)
    END IF
    CALL cp_dbcsr_copy(preconditioner_env%sparse_matrix,matrix_t,name="full_kinetic",error=error)

    shift=MAX(0.0_dp,energy_gap)

    CALL cp_dbcsr_add(preconditioner_env%sparse_matrix,matrix_s,&
                      alpha_scalar=1.0_dp,beta_scalar=shift,error=error)

    CALL timestop(handle)

  END SUBROUTINE make_full_kinetic

! *****************************************************************************
!> \brief full_single_preconditioner
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

! *****************************************************************************
!> \brief full single in the orthonormal basis
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

! *****************************************************************************
!> \brief full all in the orthonormal basis
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
!> \brief generates a preconditioner matrix H-lambda S+(SC)(2.0*CT*H*C+delta)(SC)^T
!>        for later inversion.
!>        H is the Kohn Sham matrix
!>        lambda*S shifts the spectrum of the generalized form up by lambda
!>        the last term only shifts the occupied space (reversing them in energy order)
!>        This form is implicitely multiplied from both sides by S^0.5
!>        This ensures we precondition the correct quantity
!>        Before this reads S^-0.5 H S^-0.5 + lambda + (S^0.5 C)shifts(S^0.5 C)T
!>        which might be a bit more obvious
!>        Replaced the old full_single_inverse at revision 14616
!> \param preconditioner_env the preconditioner env
!> \param matrix_c0 the MO coefficient matrix (fm)
!> \param matrix_h Kohn-Sham matrix (dbcsr)
!> \param energy_gap an additional shift in lambda=-E_homo+energy_gap
!> \param matrix_s the overlap matrix if not orthonormal (dbcsr, optional)
!> \param error ...
! *****************************************************************************
SUBROUTINE make_full_single_inverse(preconditioner_env, matrix_c0, matrix_h, energy_gap, matrix_s, error)
    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_fm_type), POINTER                :: matrix_c0
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    REAL(KIND=dp)                            :: energy_gap
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_full_single_inverse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n
    REAL(KIND=dp)                            :: max_ev, min_ev, pre_shift
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrices
    TYPE(cp_dbcsr_type), TARGET              :: dbcsr_cThc, dbcsr_hc, &
                                                dbcsr_sc, mo_dbcsr
    TYPE(dbcsr_arnoldi_data)                 :: my_arnoldi

    CALL timeset(routineN,handle)
  
    ! Allocate all working matrices needed 
    CALL cp_fm_get_info(matrix_c0,nrow_global=n,ncol_global=k,error=error)
    CALL cp_dbcsr_init(mo_dbcsr,error)
    ! copy the fm MO's to a sparse matrix, can be solved better if the sparse version is already present
    ! but for the time beeing this will do
    CALL cp_fm_to_dbcsr_row_template(mo_dbcsr,matrix_c0,matrix_h,error)
    CALL cp_dbcsr_init(dbcsr_sc,error)
    CALL cp_dbcsr_create(dbcsr_sc,template=mo_dbcsr,error=error)
    CALL cp_dbcsr_init(dbcsr_hc,error)
    CALL cp_dbcsr_create(dbcsr_hc,template=mo_dbcsr,error=error)
    CALL cp_dbcsr_init(dbcsr_cThc,error)
    CALL cp_dbcsr_m_by_n_from_template(dbcsr_cThc,matrix_h,k,k,sym=dbcsr_type_symmetric,error=error)

    ! Check whether the output matrix was already created, if not do it now
    IF(.NOT.ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       ALLOCATE(preconditioner_env%sparse_matrix)
       CALL cp_dbcsr_init(preconditioner_env%sparse_matrix,error=error)
    END IF

    ! Put the first term of the preconditioner (H) into the output matrix 
    CALL cp_dbcsr_copy(preconditioner_env%sparse_matrix,matrix_h,error=error)

    ! Precompute some matrices
    ! S*C, if orthonormal this will be simply C so a copy will do
    IF(PRESENT(matrix_s))THEN
       CALL cp_dbcsr_multiply("N", "N",1.0_dp,matrix_s,mo_dbcsr,0.0_dp,dbcsr_sc,error=error)
    ELSE
       CALL cp_dbcsr_copy(dbcsr_sc,mo_dbcsr,error=error)
    END IF

!----------------------------compute the occupied subspace and shift it ------------------------------------
    ! cT*H*C which will be used to shift the occupied states to 0
    CALL cp_dbcsr_multiply("N", "N",1.0_dp,matrix_h,mo_dbcsr,0.0_dp,dbcsr_hc,error=error)
    CALL cp_dbcsr_multiply("T", "N",1.0_dp,mo_dbcsr,dbcsr_hc,0.0_dp,dbcsr_cThc,error=error)

    ! Compute the Energy of the HOMO. We will use this as a reference energy 
    ALLOCATE(matrices(1))
    matrices(1)%matrix=>dbcsr_cThc
    CALL cp_dbcsr_setup_arnoldi_data(my_arnoldi,matrices,max_iter=20,threshold=1.0E-3_dp,selection_crit=2,&
                                     nval_request=1, nrestarts=8, generalized_ev=.FALSE.,iram=.FALSE.)
    IF(ASSOCIATED(preconditioner_env%max_ev_vector))&
         CALL cp_set_arnoldi_initial_vector(my_arnoldi,preconditioner_env%max_ev_vector)
    CALL cp_dbcsr_arnoldi_ev(matrices,my_arnoldi,error)
    max_ev=REAL(get_selected_ritz_val(my_arnoldi,1),dp)

    ! save the ev as guess for the next time
    IF(.NOT.ASSOCIATED(preconditioner_env%max_ev_vector))ALLOCATE(preconditioner_env%max_ev_vector)
    CALL get_selected_ritz_vec(my_arnoldi,1,matrices(1)%matrix,preconditioner_env%max_ev_vector,error)
    CALL deallocate_arnoldi_data(my_arnoldi)
    DEALLOCATE(matrices)

    ! Lets shift the occupied states a bit further up, -1.0 because we gonna subtract it from H  
    CALL cp_dbcsr_add_on_diag(dbcsr_cThc,-0.5_dp,error=error)
    ! Get the AO representation of the shift (see above why S is needed), W-matrix like object
    CALL cp_dbcsr_multiply("N", "N",2.0_dp,dbcsr_sc,dbcsr_cThc,0.0_dp,dbcsr_hc,error=error)
    CALL cp_dbcsr_multiply("N", "T",-1.0_dp,dbcsr_hc,dbcsr_sc,1.0_dp,preconditioner_env%sparse_matrix,&
                           error=error)

!-------------------------------------compute eigenvalues of H ----------------------------------------------
    ! Setup the arnoldi procedure to compute the lowest ev. if S is present this has to be the generalized ev
    IF(PRESENT(matrix_s))THEN
       ALLOCATE(matrices(2))
       matrices(1)%matrix=>preconditioner_env%sparse_matrix
       matrices(2)%matrix=>matrix_s
       CALL cp_dbcsr_setup_arnoldi_data(my_arnoldi,matrices,max_iter=20,threshold=2.0E-2_dp,selection_crit=3,&
                                        nval_request=1, nrestarts=21,generalized_ev=.TRUE.,iram=.FALSE.)
    ELSE
       ALLOCATE(matrices(1))
       matrices(1)%matrix=>preconditioner_env%sparse_matrix
       CALL cp_dbcsr_setup_arnoldi_data(my_arnoldi,matrices,max_iter=20,threshold=2.0E-2_dp,selection_crit=3,&
                                        nval_request=1, nrestarts=8, generalized_ev=.FALSE.,iram=.FALSE.)
    END IF
    IF(ASSOCIATED(preconditioner_env%min_ev_vector))&
       CALL cp_set_arnoldi_initial_vector(my_arnoldi,preconditioner_env%min_ev_vector)

    ! compute the LUMO energy
    CALL cp_dbcsr_arnoldi_ev(matrices,my_arnoldi,error)
    min_eV=REAL(get_selected_ritz_val(my_arnoldi,1),dp)

    ! save the lumo vector for restarting in the next step
    IF(.NOT.ASSOCIATED(preconditioner_env%min_ev_vector))ALLOCATE(preconditioner_env%min_ev_vector)
    CALL get_selected_ritz_vec(my_arnoldi,1,matrices(1)%matrix,preconditioner_env%min_ev_vector,error)
    CALL deallocate_arnoldi_data(my_arnoldi)
    DEALLOCATE(matrices)

!-------------------------------------compute eigenvalues of H ----------------------------------------------
    ! Shift the Lumo to the 1.5*the computed energy_gap or the external energy gap value
    ! The factor 1.5 is determined by trying. If the LUMO is positive, enough, just leave it alone
    pre_shift=MAX(1.5_dp*(min_ev-max_ev),energy_gap)
    IF(min_ev.LT.pre_shift)THEN
       pre_shift=pre_shift-min_ev
    ELSE
       pre_shift=0.0_dp
    END IF
    IF(PRESENT(matrix_s))THEN
       CALL cp_dbcsr_add(preconditioner_env%sparse_matrix,matrix_s,1.0_dp,pre_shift,error=error)
    ELSE
       CALL cp_dbcsr_add_on_diag(preconditioner_env%sparse_matrix,pre_shift,error=error)
    END IF

    CALL cp_dbcsr_release(mo_dbcsr,error=error)
    CALL cp_dbcsr_release(dbcsr_hc,error=error)
    CALL cp_dbcsr_release(dbcsr_sc,error=error)
    CALL cp_dbcsr_release(dbcsr_cThc,error=error)

    CALL timestop(handle)
    
END SUBROUTINE make_full_single_inverse

END MODULE preconditioner_makes

