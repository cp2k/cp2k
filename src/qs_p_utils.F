!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_p_utils *
!!
!!   NAME
!!     qs_p_utils
!!
!!   FUNCTION
!!     Utility functions for the perturbation calculations.
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     22-08-2002, TCH, started development
!!
!!   SEE ALSO
!!     
!!
!!****
MODULE qs_p_utils

  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  USE global_types,         ONLY : global_environment_type
  USE kinds,                ONLY : dbl
  USE qs_blacs
  USE qs_environment_types, ONLY : qs_environment_type, get_qs_env
  USE qs_p_types,           ONLY : qs_p_env_type
  USE qs_p_build_kernel,    ONLY : kpp1_calculate
  USE qs_p_build_kernel, ONLY: qs_K_P_P1_env_type, kpp1_calculate
  USE kinds, ONLY: wp=>dp
  USE global_types, ONLY: global_environment_type
  USE timings, ONLY: timeset, timestop
  USE qs_blacs, ONLY: allocate_blacs_matrix, deallocate_blacs_matrix,&
       blacs_gemm, blacs_add, sparse_times_blacs, allocate_blacs_matrix_vect,&
       deallocate_blacs_matrix_vect, optimal_blacs_row_block_size, &
       optimal_blacs_col_block_size, blacs_matrix_p_type, blacs_matrix_type,&
       sparse_plus_blacs_blacst, blacs_symm, get_blacs_matrix_info,&
       blacs_scale_and_d, copy_blacs_to_sparse_matrix,&
       blacs_cholesky_decompose, blacs_cholesky_restore, blacs_syrk,&
       copy_blacs_to_blacs_matrix,blacs_triangular_multiply
  USE qs_environment_types, ONLY: qs_environment_type, get_qs_env
  USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind_set
  USE dft_types, ONLY: dft_control_type
  USE cp_block_matrix, ONLY: cp_block_matrix_p_type, cp_block_matrix_type,&
       cp_sp_create, cp_sp_multiplies_blacs, cp_sp_array_dealloc,&
       cp_sp_set_to, cp_sp_array_create
  USE qs_build_ks_matrix, ONLY: qs_ks_env_type, qs_ks_did_change, &
       qs_ks_update_qs_env
  USE sparse_matrix_types, ONLY: real_matrix_p_type, real_matrix_type,&
       get_matrix_info
  USE cp_para_env, ONLY: cp_para_env_type
  USE cp_b_matrix_structure, ONLY: cp_b_matrix_struct_type

  IMPLICIT NONE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_utils'

  PRIVATE
  PUBLIC :: p_op_l1, p_env_eval_l1

CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_p_utils/p_op_l1 *
!!
!!   NAME
!!     p_op_l1
!!
!!   FUNCTION
!!     Evaluates Fv - Sv(epsilon) and stores it in res
!!
!!   ARGUMENTS
!!     p_env    - perturbation calculation environment
!!     v        - the matrix to operate on
!!     res      - the result
!!
!!   MODIFICATION HISTORY
!!
!!   SEE ALSO
!!     
!!
!!****
!!FM  subroutine p_op_l1(p_env, v, res, error)
!!FM    
!!FM    type(blacs_matrix_p_type), dimension(:), pointer   :: res, v
!!FM    type(qs_p_env_type), intent(inout)                 :: p_env
!!FM    
  SUBROUTINE p_op_l1(res, p_env, v, glob_env)
    
    TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER   :: res, v
    TYPE(qs_p_env_type), INTENT(inout)                 :: p_env
    TYPE(global_environment_type), INTENT(in)          :: glob_env

    INTEGER                                            :: i, n_spins, ncol
    TYPE(real_matrix_type), POINTER                    :: k_p

    n_spins = SIZE(v)

    DO i=1, n_spins
       k_p => p_env%qs_env%k(i)%matrix
       CALL get_matrix_info(k_p, ncol=ncol)
       CALL sparse_times_blacs(k_p, v(i)%blacs_matrix, res(i)%blacs_matrix, ncol, glob_env%para_env)       
       CALL blacs_symm('R', 'U', p_env%n_ao, p_env%n_mo, &
            1._dbl, p_env%m_epsilon(i)%blacs_matrix, v(i)%blacs_matrix, &
            0._dbl, p_env%tmp_aomo_1)
       CALL sparse_times_blacs(p_env%qs_env%s(i)%matrix, p_env%tmp_aomo_1, &
            res(i)%blacs_matrix, p_env%n_mo,&
            para_env=glob_env%para_env,alpha=1.0_dbl,beta=1.0_dbl)
    END DO

  END SUBROUTINE p_op_l1

! *****************************************************************************

!!****f* p_op_l2 [1.0] *
!!
!!   NAME
!!     p_op_l2
!!
!!   FUNCTION
!!     evaluates target_m=alpha K_P_P1 C + beta target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - target_m: place where to store the result
!!     - C: values of the coeffs of the correcting orbitals
!!     - alpha: scale factor of the result (defaults to 1.0)
!!     - beta: scale factor of the old values (defaults to 0.0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
!!FM  subroutine p_op_l2(p_env, C, target_m, alpha, beta, error)
!!FM    type(qs_p_env_type), intent(inout) :: p_env
!!FM    type(blacs_matrix_p_type), dimension(:),intent(in) :: C
!!FM    type(blacs_matrix_p_type), dimension(:),intent(inout) :: target_m
!!FM    real(kind=wp), intent(in), optional :: alpha,beta
!!FM    type(cp_error_type), optional, intent(inout) :: error
  SUBROUTINE p_op_l2(p_env, C, target_m, global_env, alpha, beta, error)
    TYPE(qs_p_env_type), INTENT(inout) :: p_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(in) :: C
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(inout) :: target_m
    REAL(kind=wp), INTENT(in), OPTIONAL :: alpha,beta
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='p_op_l2',&
         routineP=moduleN//':'//routineN
    INTEGER :: i
    failure=.FALSE.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i=1,SIZE(p_env%P1)
          CALL cp_sp_set_to(p_env%P1(i)%matrix,value=0.0_wp,error=error)
       END DO
       CALL sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,&
            p_env%c_ortho(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix,&
            ncol=p_env%n_mo)
       CALL sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,&
            p_env%psi0d(1)%blacs_matrix,p_env%c_ortho(1)%blacs_matrix,&
            ncol=p_env%n_mo)
       CALL kpp1_calculate(p_env%K_P_P1_env,p1=p_env%P1,&
            kpp1_matrix=p_env%K_P_P1,&
            global_env=global_env,error=error)

       CALL cp_sp_multiplies_blacs(p_env%K_P_P1(1)%matrix,C(1)%blacs_matrix,&
            target_m(1)%blacs_matrix, alpha=alpha,&
            beta=beta)
    END IF
  END SUBROUTINE p_op_l2
!***************************************************************************

!!****f* qs_p_types/p_env_eval_l1 [1.0] *
!!
!!   NAME
!!     p_env_eval_l1
!!
!!   FUNCTION
!!     evaluates the first part of the linear functional:
!!       r = Hrho C - S C epsilon = Hrho C + S C m_epsilon
!!
!!   NOTES
!!     uses tmp_aomo_1, Hrho must already be valid
!!
!!   INPUTS
!!     - p_env: the environment where to perform the calculation
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE p_env_eval_l1(p_env,global_env,error)
    TYPE(qs_p_env_type), INTENT(inout) :: p_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='p_env_eval_l1',&
         routineP=moduleN//':'//routineN
    TYPE(real_matrix_p_type), DIMENSION(:),POINTER :: s
    failure=.FALSE.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(p_env%qs_env,s=s)
       CALL cp_sp_multiplies_blacs(p_env%Hrho(1)%matrix,&
            p_env%c_ortho(1)%blacs_matrix,p_env%r(1)%blacs_matrix)
       CALL blacs_symm('R','U',p_env%n_ao,p_env%n_mo,1.0_wp,&
            p_env%m_epsilon(1)%blacs_matrix,&
            p_env%c_ortho(1)%blacs_matrix,0.0_wp,p_env%tmp_aomo_1)
       CALL sparse_times_blacs(s(1)%matrix,p_env%tmp_aomo_1, &
            p_env%r(1)%blacs_matrix,p_env%n_mo,&
            para_env=global_env%para_env,alpha=1.0_wp,beta=1.0_wp)
    END IF
  END SUBROUTINE p_env_eval_l1
!***************************************************************************

!!****f* p_env_eval_l2 [1.0] *
!!
!!   NAME
!!     p_env_eval_l2
!!
!!   FUNCTION
!!     evaluates target_m=alpha K_P_P1 C + beta target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - target_m: place where to store the result
!!     - C: values of the coeffs of the correcting orbitals
!!     - alpha: scale factor of the result (defaults to 1.0)
!!     - beta: scale factor of the old values (defaults to 0.0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE p_env_eval_l2(p_env, target_m, C, global_env, alpha, beta, error)
    TYPE(qs_p_env_type), INTENT(inout) :: p_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(in) :: C
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(inout) :: target_m
    REAL(kind=wp), INTENT(in), OPTIONAL :: alpha,beta
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='p_env_eval_l2',&
         routineP=moduleN//':'//routineN
    INTEGER :: i
    failure=.FALSE.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i=1,SIZE(p_env%P1)
          CALL cp_sp_set_to(p_env%P1(i)%matrix,value=0.0_wp,error=error)
       END DO
       CALL sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,&
            p_env%c_ortho(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix,&
            ncol=p_env%n_mo)
       CALL sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,&
            p_env%psi0d(1)%blacs_matrix,p_env%c_ortho(1)%blacs_matrix,&
            ncol=p_env%n_mo)
       CALL kpp1_calculate(p_env%K_P_P1_env,p1=p_env%P1,&
            kpp1_matrix=p_env%K_P_P1,&
            global_env=global_env,error=error)

       CALL cp_sp_multiplies_blacs(p_env%K_P_P1(1)%matrix,C(1)%blacs_matrix,&
            target_m(1)%blacs_matrix, alpha=alpha,&
            beta=beta)
    END IF
  END SUBROUTINE p_env_eval_l2
!***************************************************************************

!!****f* qs_p_types/p_env_preortho [1.0] *
!!
!!   NAME
!!     p_env_preortho
!!
!!   FUNCTION
!!     does a preorthogonalization of the given matrix:
!!     C=(I-PS)C
!!
!!   NOTES
!!     uses tmp_aomo_1
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - C: matrix to orthogonalize
!!     - n_cols: the number of columns of C to multiply (defaults to size(C,2))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE p_env_preortho(p_env, C, global_env, &
       n_cols,error)
    TYPE(qs_p_env_type), INTENT(inout) :: p_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(inout) :: C
    INTEGER, INTENT(in), OPTIONAL :: n_cols
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='p_env_preortho',&
         routineP=moduleN//':'//routineN
    INTEGER :: nrow_global, ncol_global, cols
    failure=.FALSE.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    CALL get_blacs_matrix_info(C(1)%blacs_matrix,nrow_global=nrow_global, &
         ncol_global=ncol_global)
    CPPrecondition(nrow_global>=p_env%n_ao,cp_failure_level,routineP,error,failure)
    cols=ncol_global
    IF (PRESENT(n_cols)) THEN
       CPPrecondition(n_cols<=cols,cp_failure_level,routineP,error,failure)
       cols=n_cols
    END IF
    CPPrecondition(cols<=p_env%n_ao,cp_failure_level,routineP,error,failure)
! check target_m?
    IF (.NOT. failure) THEN
       CALL blacs_gemm('N','T',cols,p_env%n_mo,p_env%n_ao,1.0_wp,&
            C(1)%blacs_matrix,p_env%S_psi0(1)%blacs_matrix,&
            0.0_wp,p_env%tmp_aomo_1)
       CALL blacs_gemm('N','T',p_env%n_ao,cols,p_env%n_ao,-1.0_wp,&
            p_env%psi0d(1)%blacs_matrix,p_env%tmp_aomo_1,&
            1.0_wp,C(1)%blacs_matrix)
    END IF
  END SUBROUTINE p_env_preortho
!***************************************************************************

!!****f* qs_p_types/p_env_postortho [1.0] *
!!
!!   NAME
!!     p_env_postortho
!!
!!   FUNCTION
!!     does a postorthogonalization on the given matrix vector:
!!     target_m=(I-SP)C
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - C: matrix to orthogonalize
!!     - n_cols: the number of columns of C to multiply (defaults to size(C,2))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE p_env_postortho(p_env,target_m, C, global_env,&
       n_cols,error)
    TYPE(qs_p_env_type), INTENT(inout) :: p_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(in) :: C
    TYPE(blacs_matrix_p_type), DIMENSION(:),INTENT(inout) :: target_m
    INTEGER, INTENT(in), OPTIONAL :: n_cols
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='p_env_postortho',&
         routineP=moduleN//':'//routineN
    INTEGER :: nrow_global, ncol_global, cols
    failure=.FALSE.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    CALL get_blacs_matrix_info(C(1)%blacs_matrix,nrow_global=nrow_global,&
         ncol_global=ncol_global)
    CPPrecondition(nrow_global>=p_env%n_ao,cp_failure_level,routineP,error,failure)
    cols=ncol_global
    IF (PRESENT(n_cols)) THEN
       CPPrecondition(n_cols<=cols,cp_failure_level,routineP,error,failure)
       cols=n_cols
    END IF
    CPPrecondition(cols<=p_env%n_ao,cp_failure_level,routineP,error,failure)
    ! check target_m?
    IF (.NOT. failure) THEN     
       CALL blacs_gemm('N','T',cols,p_env%n_mo,p_env%n_ao,1.0_wp,&
            C(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix,&
            0.0_wp,p_env%tmp_aomo_1)
       CALL blacs_gemm('N','T',p_env%n_ao,cols,p_env%n_ao,-1.0_wp,&
            p_env%S_psi0(1)%blacs_matrix,p_env%tmp_aomo_1,&
            1.0_wp,C(1)%blacs_matrix)
    END IF
  END SUBROUTINE p_env_postortho
!***************************************************************************

END MODULE qs_p_utils
