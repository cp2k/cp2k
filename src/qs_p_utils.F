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
module qs_p_utils

  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use global_types,         only : global_environment_type
  use kinds,                only : dbl
  use qs_blacs
  use qs_environment_types, only : qs_environment_type, get_qs_env
  use qs_p_types,           only : qs_p_env_type
  use qs_p_build_kernel,    only : kpp1_calculate
  use qs_p_build_kernel, only: qs_K_P_P1_env_type, kpp1_calculate
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use qs_blacs, only: allocate_blacs_matrix, deallocate_blacs_matrix,&
       blacs_gemm, blacs_add, sparse_times_blacs, allocate_blacs_matrix_vect,&
       deallocate_blacs_matrix_vect, optimal_blacs_row_block_size, &
       optimal_blacs_col_block_size, blacs_matrix_p_type, blacs_matrix_type,&
       sparse_plus_blacs_blacst, blacs_symm, get_blacs_matrix_info,&
       blacs_scale_and_d, copy_blacs_to_sparse_matrix,&
       blacs_cholesky_decompose, blacs_cholesky_restore, blacs_syrk,&
       copy_blacs_to_blacs_matrix,blacs_triangular_multiply
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind_set
  use dft_types, only: dft_control_type
  use cp_block_matrix, only: cp_block_matrix_p_type, cp_block_matrix_type,&
       cp_sp_create, cp_sp_multiplies_blacs, cp_sp_array_dealloc,&
       cp_sp_set_to, cp_sp_array_create
  use qs_build_ks_matrix, only: qs_ks_env_type, qs_ks_did_change, &
       qs_ks_update_qs_env
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type,&
       get_matrix_info
  use cp_para_env, only: cp_para_env_type
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type

  implicit none

  character(len=*), private, parameter :: moduleN='qs_p_utils'

  private
  public :: p_op_l1

contains

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
!!     unrestricted orbitals work (maybe)
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
!!     02.09.2002 adapted to new qs_p_env_type (TC)
!!
!!****
!!FM  subroutine p_op_l1(p_env, v, res, error)
!!FM    
!!FM    type(blacs_matrix_p_type), dimension(:), pointer   :: res, v
!!FM    type(qs_p_env_type), intent(inout)                 :: p_env
!!FM    
  subroutine p_op_l1(p_env, v, res, glob_env, error)
    
    ! argument
    type(qs_p_env_type), intent(inout)                     :: p_env
    type(blacs_matrix_p_type), dimension(:), intent(in)    :: v
    type(blacs_matrix_p_type), dimension(:), intent(inout) :: res
    type(global_environment_type), intent(in)              :: glob_env
    type(cp_error_type), optional, intent(inout)           :: error
                                                           
    ! locals                                               
    logical                          :: failure
    integer                          :: spin, n_spins, ncol
    type(real_matrix_type), pointer  :: k_p
    type(blacs_matrix_type), pointer :: tmp
    character(len=*), parameter      :: routineN = 'p_op_l1',&
                                        routineP = moduleN//':'//routineN

    ! code
    failure = .false.
    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)

    n_spins = size(v)
    nullify(tmp)
    !! call get_tmp_from_pool

    do spin=1, n_spins
       k_p => p_env%qs_env%k(spin)%matrix
       call get_matrix_info(k_p, ncol=ncol)
       call sparse_times_blacs(k_p, v(spin)%blacs_matrix, &
                               res(spin)%blacs_matrix, ncol, glob_env%para_env)
       call blacs_symm('R', 'U', p_env%n_ao(spin), p_env%n_mo(spin), 1._dbl, &
                       p_env%m_epsilon(spin)%blacs_matrix, &
                       v(spin)%blacs_matrix, 0._dbl, tmp)
       call sparse_times_blacs(p_env%qs_env%s(spin)%matrix, tmp, &
                               res(spin)%blacs_matrix, p_env%n_mo(spin),&
                               para_env=glob_env%para_env, &
                               alpha=1.0_dbl,beta=1.0_dbl)
    end do

  end subroutine p_op_l1

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
!!     is not yet capable of handling spins !!!
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
!!FM  subroutine p_op_l2(p_env, C, target_m, alpha, beta, error)
!!FM    type(qs_p_env_type), intent(inout) :: p_env
!!FM    type(blacs_matrix_p_type), dimension(:),intent(in) :: C
!!FM    type(blacs_matrix_p_type), dimension(:),intent(inout) :: target_m
!!FM    real(kind=wp), intent(in), optional :: alpha,beta
!!FM    type(cp_error_type), optional, intent(inout) :: error

  subroutine p_op_l2(p_env, v, res, glob_env, alpha, beta, error)
    
    ! arguments
    type(qs_p_env_type), intent(inout)                     :: p_env
    type(blacs_matrix_p_type), dimension(:), intent(in)    :: v
    type(blacs_matrix_p_type), dimension(:), intent(inout) :: res
    type(global_environment_type), intent(in), target      :: glob_env
    real(kind=dbl), intent(in), optional                   :: alpha, beta
    type(cp_error_type), optional, intent(inout)           :: error

    ! locals
    logical                     :: failure
    character(len=*), parameter :: routineN = 'p_op_l2',&
                                   routineP = moduleN//':'//routineN
    integer                     :: n_spins, spin

!!!    integer :: i
    failure=.false.

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)

    if (.not. failure) then
       
!!!       do i=1, size(v)
!!!           call cp_sp_set_to(p_env%P1(spin)%matrix, value=0.0_dbl, error=error)
!!!       end do
!!!       call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,&
!!!               p_env%c_ortho(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix,&
!!!               ncol=p_env%n_mo)
!!!       call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,&
!!!            p_env%psi0d(1)%blacs_matrix,p_env%c_ortho(1)%blacs_matrix,&
!!!            ncol=p_env%n_mo)
!!!       call kpp1_calculate(p_env%K_P_P1_env,p1=p_env%P1,&
!!!            kpp1_matrix=p_env%K_P_P1,&
!!!            global_env=global_env,error=error)
!!!       call cp_sp_multiplies_blacs(p_env%K_P_P1(1)%matrix,C(1)%blacs_matrix,&
!!!            target_m(1)%blacs_matrix, alpha=alpha,&
!!!            beta=beta)
    end if

  end subroutine p_op_l2
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
!!     is not yet capable of handling spins !
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
  subroutine p_preortho(p_env, v, glob_env, n_cols, error)
    
    ! arguments
    type(qs_p_env_type), intent(inout)                     :: p_env
    type(blacs_matrix_p_type), dimension(:), intent(inout) :: v
    type(global_environment_type), intent(in), target      :: glob_env
    integer, intent(in), optional                          :: n_cols
    type(cp_error_type), optional, intent(inout)           :: error

    ! locals
    logical                          :: failure
    character(len=*), parameter      :: routineN = 'p_env_preortho',&
                                        routineP = moduleN//':'//routineN
    integer                          :: v_rows, v_cols, cols
    type(blacs_matrix_type), pointer :: tmp

    ! code
    failure=.false.

    nullify(tmp)
    ! call get_tmp_from_pool

    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)
    call get_blacs_matrix_info(v(1)%blacs_matrix,nrow_global=v_rows, &
                               ncol_global=v_cols)
    CPPrecondition(v_rows >= p_env%n_ao(1), &
                   cp_failure_level, routineP,error,failure)
    cols=v_cols
    if (present(n_cols)) then
       CPPrecondition(n_cols <= cols, &
                      cp_failure_level, routineP, error, failure)
       cols=n_cols
    end if
    CPPrecondition(cols <= p_env%n_ao(1), &
                   cp_failure_level, routineP, error, failure)

    if (.not. failure) then
       call blacs_gemm('N', 'T', cols, p_env%n_mo(1), p_env%n_ao(1) ,1.0_dbl, &
                       v(1)%blacs_matrix, p_env%S_psi0(1)%blacs_matrix, & 
                       0.0_dbl, tmp)
       call blacs_gemm('N', 'T', p_env%n_ao(1), cols, p_env%n_ao(1), -1.0_dbl,&
                       p_env%psi0d(1)%blacs_matrix, tmp, &
                       0.0_dbl, v(1)%blacs_matrix)
    end if
  end subroutine p_preortho
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
!!     
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

  subroutine p_postortho(p_env, v, res, glob_env, n_cols, error)

    ! arguments
    type(qs_p_env_type), intent(inout)                     :: p_env
    type(blacs_matrix_p_type), dimension(:), intent(in)    :: v
    type(blacs_matrix_p_type), dimension(:), intent(inout) :: res
    type(global_environment_type), intent(in), target      :: glob_env
    integer, intent(in), optional :: n_cols
    type(cp_error_type), optional, intent(inout) :: error

    ! locals
    logical                          :: failure
    character(len=*), parameter      :: routineN = 'p_env_postortho',&
                                        routineP = moduleN//':'//routineN
    integer                          :: v_rows, v_cols, cols
    integer                          :: n_spins
    type(blacs_matrix_type), pointer :: tmp

    ! code
    failure=.false.

    nullify(tmp)
    ! call get_tmp_from_pool

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    call get_blacs_matrix_info(v(1)%blacs_matrix, nrow_global=v_rows,&
                                                  ncol_global=v_cols)
    CPPrecondition(v_rows >= p_env%n_ao(1), &
                   cp_failure_level, routineP, error, failure)

    if (present(n_cols)) then
       CPPrecondition(n_cols <= v_cols, &
                      cp_failure_level,routineP,error,failure)
       cols = n_cols
    else
       cols = v_cols
    end if
    CPPrecondition(cols <= p_env%n_ao(1), &
                   cp_failure_level,routineP,error,failure)
    ! check res?
    if (.not. failure) then     
       call blacs_gemm('N', 'T', cols, &
                       p_env%n_mo(1), p_env%n_ao(1), 1.0_dbl,&
                       v(1)%blacs_matrix, p_env%psi0d(1)%blacs_matrix, &
                       0.0_dbl, tmp)
       call blacs_gemm('N', 'T', p_env%n_ao(1) ,cols, p_env%n_ao(1), -1.0_dbl, &
                       p_env%S_psi0(1)%blacs_matrix, tmp,&
                       1.0_dbl, v(1)%blacs_matrix)
    end if
  end subroutine p_postortho

!***************************************************************************

end module qs_p_utils
