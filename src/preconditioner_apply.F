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
MODULE preconditioner_apply
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_copy, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_release, cp_dbcsr_type
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE input_constants,                 ONLY: ot_precond_full_all,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_full_single_inverse,&
                                             ot_precond_s_inverse,&
                                             ot_precond_solver_direct,&
                                             ot_precond_solver_inv_chol,&
                                             ot_precond_solver_update
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

  PUBLIC :: apply_preconditioner_fm, apply_preconditioner_dbcsr

CONTAINS

! *****************************************************************************
!> \brief applies a previously created preconditioner to a full matrix
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
       CASE(ot_precond_solver_inv_chol,ot_precond_solver_update)
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
       CASE(ot_precond_solver_inv_chol,ot_precond_solver_update)
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

! *****************************************************************************
!> \brief apply to full matrix, complete inversion has already been done
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
!> \brief apply to dbcsr matrix, complete inversion has already been done
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
!> \brief preconditioner contains the factorization, application done by
!>        solving the linear system
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
!> \brief full all to a full matrix
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
!> \brief full all to a dbcsr matrix
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

END MODULE preconditioner_apply
