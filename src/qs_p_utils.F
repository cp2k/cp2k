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
!!   NOTES
!!     - use a pool of temporary matrizes instead of
!!       allocating them every time needed
!!     - routines are programmed with spins in mind
!!       but are as of now not able to handle them
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
  USE kinds,                ONLY : dbl
  USE qs_blacs
  USE qs_environment_types, ONLY : qs_environment_type, get_qs_env
  USE qs_p_types,           ONLY : qs_p_env_type
  USE qs_p_build_kernel, ONLY: kpp1_calculate
  USE kinds, ONLY: wp=>dp
  USE global_types, ONLY: global_environment_type
  USE timings, ONLY: timeset, timestop
  USE qs_blacs, ONLY: cp_fm_create2, cp_fm_release,&
       cp_fm_gemm, cp_fm_add, cp_sm_fm_multiply, cp_fm_vect_create2,&
       cp_fm_vect_dealloc, optimal_blacs_row_block_size, &
       optimal_blacs_col_block_size, cp_full_matrix_p_type, cp_full_matrix_type,&
       cp_sm_plus_fm_fm_t, cp_fm_symm, cp_fm_get_info,&
       blacs_scale_and_d, copy_blacs_to_sparse_matrix,&
       blacs_cholesky_decompose, blacs_cholesky_restore, cp_fm_syrk,&
       copy_blacs_to_blacs_matrix,blacs_triangular_multiply, &
       replicate_blacs_matrix
  USE qs_environment_types, ONLY: qs_environment_type
  USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind_set
  USE dft_types, ONLY: dft_control_type
  USE cp_block_matrix, ONLY: cp_block_matrix_p_type, cp_block_matrix_type,&
       cp_sp_create, cp_sp_multiplies_blacs, cp_sp_array_dealloc,&
       cp_sp_set_to, cp_sp_array_create
  USE qs_build_ks_matrix, ONLY: qs_ks_env_type, qs_ks_did_change, &
       qs_ks_update_qs_env
  USE sparse_matrix_types, ONLY: real_matrix_p_type, real_matrix_type,&
       get_matrix_info
  USE cp_para_types, ONLY: cp_para_env_type
  USE cp_b_matrix_structure, ONLY: cp_b_matrix_struct_type

  IMPLICIT NONE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_utils'

  PRIVATE
  PUBLIC :: p_op_l1, p_op_l1_spin, p_op_l2, p_preortho, p_postortho

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
!!   NOTES
!!     -
!!     
!!   ARGUMENTS
!!     - p_env : perturbation calculation environment
!!     - v     : the matrix to operate on
!!     - res   : the result
!!     - glob_env : should eventually vanish
!      - error : error handling object (optional)
!!
!!   AUTHOR
!!     Thomas Chassaing
!!
!!   MODIFICATION HISTORY
!!     10.2002, TCH, extracted single spin calculation
!!
!!****
  SUBROUTINE p_op_l1(p_env, v, res, glob_env, error)
 
    IMPLICIT NONE
   
    ! argument
    TYPE(qs_p_env_type), INTENT(in)                          :: p_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(in)    :: v
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(inout) :: res
    TYPE(global_environment_type), INTENT(in)                :: glob_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)             :: error
                                                           
    ! locals                                               
    LOGICAL                          :: failure
    INTEGER                          :: n_spins, spin, ncol
    TYPE(real_matrix_type), POINTER  :: k_p
    TYPE(cp_full_matrix_type), POINTER :: tmp
    CHARACTER(len=*), PARAMETER      :: routineN = 'p_op_l1',&
                                        routineP = moduleN//':'//routineN

    ! code
    NULLIFY(tmp);
    failure = .FALSE.
    
    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)

    n_spins = p_env%qs_env%dft_control%nspins
    DO spin=1, n_spins

       CALL p_op_l1_spin(p_env, spin, v(spin)%matrix, res(spin)%matrix, glob_env, error)
      
    END DO

  END SUBROUTINE p_op_l1

! *****************************************************************************

!!****f* cp2k/qs_p_utils/p_op_l1 *
!!
!!   NAME
!!     p_op_l1_spin
!!
!!   FUNCTION
!!     Evaluates Fv - Sv(epsilon) and stores it in res for a given spin
!!
!!   NOTES
!!     Same as p_op_l1 but takes a spin as additional argument.
!!     
!!   ARGUMENTS
!!     - p_env : perturbation calculation environment
!!     - spin  : the spin to calculate (1 or 2 normally)
!!     - v     : the matrix to operate on
!!     - res   : the result
!!     - glob_env : should eventually vanish
!      - error : error handling object (optional)
!!
!!   AUTHOR
!!     Thomas Chassaing
!!
!!   MODIFICATION HISTORY
!!     10.2002, TCH, created
!!
!!****
  SUBROUTINE p_op_l1_spin(p_env, spin, v, res, glob_env, error)
 
    IMPLICIT NONE
   
    ! argument
    TYPE(qs_p_env_type), INTENT(in)              :: p_env
    INTEGER, INTENT(IN)                          :: spin
    TYPE(cp_full_matrix_type), POINTER           :: v
    TYPE(cp_full_matrix_type), POINTER           :: res
    TYPE(global_environment_type), INTENT(in)    :: glob_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
                                                           
    ! locals                                               
    LOGICAL                          :: failure
    INTEGER                          :: ncol
    TYPE(real_matrix_type), POINTER  :: k_p
    TYPE(cp_full_matrix_type), POINTER :: tmp
    CHARACTER(len=*), PARAMETER      :: routineN = 'p_op_l1',&
                                        routineP = moduleN//':'//routineN

    ! code
    NULLIFY(tmp);
    failure = .FALSE.
    
    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)
    CPPrecondition(0<spin.AND.spin<=p_env%qs_env%dft_control%nspins, cp_failure_level, routineP, error, failure)

    CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, tmp, &
         routineP//": temporary matrix")
    
    k_p => p_env%qs_env%k(spin)%matrix
    CALL get_matrix_info(k_p, ncol=ncol)
    CALL cp_sm_fm_multiply(k_p, v, res, ncol, glob_env%para_env)
    CALL cp_fm_symm('R', 'U', p_env%n_ao(spin), p_env%n_mo(spin), 1._dbl, &
         p_env%m_epsilon(spin)%matrix, v, 0._dbl, tmp)
    CALL cp_sm_fm_multiply(p_env%qs_env%s(spin)%matrix, tmp, &
         res, p_env%n_mo(spin), para_env=glob_env%para_env, alpha=1.0_dbl,beta=1.0_dbl)    
    CALL cp_fm_release(tmp); NULLIFY(tmp)

  END SUBROUTINE p_op_l1_spin

! *****************************************************************************

!!****f* p_op_l2 [1.0] *
!!
!!   NAME
!!     p_op_l2
!!
!!   FUNCTION
!!     evaluates res = alpha K_P_P1 v + beta res
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - res: place where to store the result
!!     - v: values of the coeffs of the correcting orbitals
!!     - alpha: scale factor of the result (defaults to 1.0)
!!     - beta: scale factor of the old values (defaults to 0.0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     02.09.2002 adapted for new qs_p_env_type (TC)
!!
!!*** **********************************************************************

  SUBROUTINE p_op_l2(p_env, v, res, glob_env, alpha, beta, error)
    
    ! arguments
    TYPE(qs_p_env_type), INTENT(inout)                       :: p_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(in)    :: v
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(inout) :: res
    TYPE(global_environment_type), INTENT(in), TARGET        :: glob_env
    REAL(kind=dbl), INTENT(in), OPTIONAL                     :: alpha, beta
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)             :: error

    ! locals
    LOGICAL                     :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_l2',&
                                   routineP = moduleN//':'//routineN
    INTEGER                     :: spin, n_spins

    failure=.FALSE.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       
!!!       do i=1, size(v)
!!!           call cp_sp_set_to(p_env%P1(spin)%matrix, value=0.0_dbl, error=error)
!!!       end do
!!!       call cp_sm_plus_fm_fm_t(p_env%P1(1)%matrix%sm,&
!!!               p_env%c_ortho(1)%matrix,p_env%psi0d(1)%matrix,&
!!!               ncol=p_env%n_mo)
!!!       call cp_sm_plus_fm_fm_t(p_env%P1(1)%matrix%sm,&
!!!            p_env%psi0d(1)%matrix,p_env%c_ortho(1)%matrix,&
!!!            ncol=p_env%n_mo)
!!!       call kpp1_calculate(p_env%K_P_P1_env,p1=p_env%P1,&
!!!            kpp1_matrix=p_env%K_P_P1,&
!!!            global_env=global_env,error=error)
!!!       call cp_sp_multiplies_blacs(p_env%K_P_P1(1)%matrix,C(1)%matrix,&
!!!            target_m(1)%matrix, alpha=alpha,&
!!!            beta=beta)
       n_spins = p_env%qs_env%dft_control%nspins
       do spin=1, n_spins
          
       end do
    END IF

  END SUBROUTINE p_op_l2
!***************************************************************************

!!****f* qs_p_types/p_env_preortho [1.0] *
!!
!!   NAME
!!     p_preortho
!!
!!   FUNCTION
!!     does a preorthogonalization of the given matrix:
!!     v = (I-PS)v
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - v: matrix to orthogonalize
!!     - n_cols: the number of columns of C to multiply (defaults to size(v,2))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     02.09.2002 adapted for new qs_p_env_type (TC)
!!
!!*** **********************************************************************
  SUBROUTINE p_preortho(p_env, v, glob_env, n_cols, error)
    
    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(in)                        :: p_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER     :: v
    TYPE(global_environment_type), INTENT(in), TARGET      :: glob_env
    INTEGER, INTENT(in), OPTIONAL                          :: n_cols
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)           :: error

    ! locals
    LOGICAL                          :: failure
    CHARACTER(len=*), PARAMETER      :: routineN = 'p_env_preortho',&
                                        routineP = moduleN//':'//routineN
    INTEGER                          :: v_rows, v_cols, cols
    INTEGER                          :: n_spins, spin
    TYPE(cp_full_matrix_type), POINTER :: tmp

    ! code
    failure=.FALSE.

    NULLIFY(tmp)

    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)
    
    n_spins = p_env%qs_env%dft_control%nspins
    DO spin=1, n_spins

       CALL cp_fm_get_info(v(1)%matrix, &
                           nrow_global=v_rows, ncol_global=v_cols)
       CPPrecondition(v_rows >= p_env%n_ao(spin), cp_failure_level, routineP, error, failure)
       cols = v_cols
       IF (PRESENT(n_cols)) THEN
          CPPrecondition(n_cols <= cols, cp_failure_level, routineP, error, failure)
          cols=n_cols
       END IF
       CPPrecondition(cols <= p_env%n_ao(spin), cp_failure_level, routineP, error, failure)

       ! the tmp matrix here is of dimension n_ao * v_cols
       ! but n_mo * v_cols would suffice
       CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, tmp, &
                                   routineP//": temporary matrix")

       IF (.NOT. failure) THEN
          CALL cp_fm_gemm('T', 'N', p_env%n_mo(spin), cols, p_env%n_ao(spin), &
                          1.0_dbl, p_env%S_psi0(spin)%matrix, v(spin)%matrix, &
                          0.0_dbl, tmp)
          CALL cp_fm_gemm('N', 'N', p_env%n_ao(spin), cols, p_env%n_mo(spin), &
                          -1.0_dbl, p_env%psi0d(spin)%matrix, tmp, &
                          1.0_dbl, v(spin)%matrix)
          
       END IF

       CALL cp_fm_release(tmp); NULLIFY(tmp)

    END DO
    
  END SUBROUTINE p_preortho
!***************************************************************************

!!****f* qs_p_types/p_env_postortho [1.0] *
!!
!!   NAME
!!     p_postortho
!!
!!   FUNCTION
!!     does a postorthogonalization on the given matrix vector:
!!     res = (I-SP) v
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - v: matrix to orthogonalize
!!     - res: result
!!     - n_cols: the number of columns of C to multiply (defaults to size(v,2))
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

  SUBROUTINE p_postortho(p_env, v, res, glob_env, n_cols, error)

    ! arguments
    TYPE(qs_p_env_type), INTENT(inout)                     :: p_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(in)    :: v
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(inout) :: res
    TYPE(global_environment_type), INTENT(in), TARGET      :: glob_env
    INTEGER, INTENT(in), OPTIONAL :: n_cols
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL                          :: failure
    CHARACTER(len=*), PARAMETER      :: routineN = 'p_env_postortho',&
                                        routineP = moduleN//':'//routineN
    INTEGER                          :: v_rows, v_cols, cols
    INTEGER                          :: n_spins, spin
    TYPE(cp_full_matrix_type), POINTER :: tmp

    ! code
    failure=.FALSE.

    NULLIFY(tmp)

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)

    n_spins = p_env%qs_env%dft_control%nspins
    DO spin=1, n_spins

       CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, tmp, &
                                   routineP//": temporary matrix")

       CALL cp_fm_get_info(v(1)%matrix, &
                                  nrow_global = v_rows, ncol_global = v_cols)
       CPPrecondition(v_rows >= p_env%n_ao(1), cp_failure_level, routineP, error, failure)

       IF (PRESENT(n_cols)) THEN
          CPPrecondition(n_cols <= v_cols, cp_failure_level, routineP, error, failure)
          cols = n_cols
       ELSE
          cols = v_cols
       END IF
       CPPrecondition(cols <= p_env%n_ao(spin), cp_failure_level, routineP, error, failure)

       ! check result matrix ?
       IF (.NOT. failure) THEN     
!!          CALL cp_fm_gemm('N', 'T', cols, p_env%n_mo(spin), p_env%n_ao(spin), 1.0_dbl,&
!!                          v(spin)%matrix, p_env%psi0d(spin)%matrix, &
!!                          0.0_dbl, tmp)
!!          CALL cp_fm_gemm('N', 'T', p_env%n_ao(spin), cols, p_env%n_ao(spin), -1.0_dbl, &
!!                          p_env%S_psi0(spin)%matrix, tmp, &
!!                          1.0_dbl, v(spin)%matrix)
          call cp_fm_gemm('N', 'T', p_env%n_mo(spin), cols, p_env%n_ao(spin), &
                          1.0_dbl, p_env%psi0d(spin)%matrix, v(spin)%matrix, &
                          0.0_dbl, tmp);
          call cp_fm_gemm('N', 'N', p_env%n_ao(spin), cols, p_env%n_mo(spin), &
                          -1.0_dbl, p_env%S_psi0(spin)%matrix, tmp, &
                          1.0_dbl, v(spin)%matrix)
       END IF
       
       CALL cp_fm_release(tmp); NULLIFY(tmp)

    END DO

  END SUBROUTINE p_postortho

!***************************************************************************

END MODULE qs_p_utils
