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
MODULE preconditioner_makes
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
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
                                             cp_fm_set_all,&
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
  USE kinds,                           ONLY: dp,&
                                             sp
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
       IF ( PRESENT(matrix_s) ) THEN
          CALL make_full_single_inverse(preconditioner_env,mo_coeff,matrix_h, matrix_s, &
                                        eigenvalues_ot, energy_gap, my_mixed_precision, &
                                        error=error)
       ELSE
          CALL make_full_single_inverse_ortho(preconditioner_env,mo_coeff,matrix_h, &
                                              eigenvalues_ot, energy_gap,error=error)
       END IF
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
!> \brief full_sinlge inverse in the orthonormal basis
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

END MODULE preconditioner_makes

