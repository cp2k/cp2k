!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_p_types [1.0] *
!!
!!   NAME
!!     qs_p_types
!!
!!   FUNCTION
!!     basis types for the calcualtion of the perturbation of density theory.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
module qs_p_types
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use qs_blacs, only: allocate_blacs_matrix, deallocate_blacs_matrix,&
       blacs_gemm, blacs_add, sparse_times_blacs, allocate_blacs_matrix_vect,&
       deallocate_blacs_matrix_vect, optimal_blacs_row_block_size, &
       optimal_blacs_col_block_size, blacs_matrix_p_type, blacs_matrix_type,&
       sparse_plus_blacs_blacst, blacs_symm, get_blacs_matrix_info
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind_set
  use dft_types, only: dft_control_type
  use cp_block_matrix, only: cp_block_matrix_p_type, cp_block_matrix_type,&
       cp_sp_create, cp_sp_multiplies_blacs, cp_sp_array_dealloc,&
       cp_sp_set_to
  use qs_build_ks_matrix, only: qs_ks_env_type
  use qs_p_build_kernel, only: qs_K_P_P1_env_type, kpp1_calculate
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type,&
       get_matrix_info
  use cp_mpi_env, only: cp_mpi_env_type, cp_mpi_get_from_glob
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_p_types'

!***
!****************************************************************************

  !!****s* qs_p_types/qs_p_env_type [1.0] *
  !!
  !!   NAME
  !!     qs_p_env_type
  !!
  !!   FUNCTION
  !!     Represent a qs system that is perturbed.
  !!     Can calculate the linear operator and the rhs of the system 
  !!     of equations that needs to be solved for the perturbation.
  !!
  !!   NOTES
  !!     for the moment no smearing of the orbitals
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     4.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type qs_p_env_type
     logical :: initialized, s_struct_changed
     type(cp_block_matrix_p_type), dimension(:), pointer :: Hrho, P, P1,&
          K_P_P1
     type(blacs_matrix_p_type), dimension(:), pointer :: I_SP, m_epsilon,&
          rhs, c_ortho, r, psi0d
     type(blacs_matrix_type), pointer :: tmp_aomo_1
     type(qs_ks_env_type), pointer :: ks_env
     type(qs_K_P_P1_env_type), pointer :: K_P_P1_env
     type(qs_environment_type), pointer :: qs_env
     integer n_mo,n_ao
  end type qs_p_env_type
  !!***
  !****************************************************************************

contains

!!****f* qs_p_types/p_env_init [1.0] *
!!
!!   NAME
!!     p_env_init
!!
!!   FUNCTION
!!     initializes the perturbation environment (no setup)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the environment to initialize
!!     - qs_env: the qs_environment for the full system
!!     - ks_env: the environment that builds the ks matrix (and H_rho)
!!     - K_P_P1_env: the environment that builds the second order
!!       perturbation kernel.
!!     - global_env: the global environment
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
subroutine p_env_init(p_env, qs_env, ks_env, K_P_P1_env, global_env,&
     c_ortho_ptr, psi0d_ptr, rhs_ptr, I_SP_ptr, error)
  type(qs_p_env_type), intent(out) :: p_env
  type(qs_environment_type), intent(in), target :: qs_env
  type(qs_ks_env_type), intent(in), target :: ks_env
  type(qs_K_P_P1_env_type), intent(in), target :: K_P_P1_env
  type(global_environment_type), intent(in), target :: global_env
  type(blacs_matrix_p_type), dimension(:),pointer, optional :: c_ortho_ptr,&
       psi0d_ptr, rhs_ptr, I_SP_ptr
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='p_env_init',&
        routineP=moduleN//':'//routineN
  integer :: n_ao, n_mo
  type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
  type(dft_control_type), pointer :: dft_control
  failure=.false.
  nullify(atomic_kind_set, dft_control)
  
  nullify(p_env%psi0d, p_env%c_ortho, p_env%r, p_env%rhs,&
       p_env%tmp_aomo_1)
  nullify(p_env%I_SP, p_env%m_epsilon)
  nullify(p_env%ks_env, p_env%K_P_P1_env,p_env%qs_env)
  nullify(p_env%Hrho, p_env%K_P_P1, p_env%P, p_env%P1)
  p_env%qs_env => qs_env
  p_env%ks_env => ks_env
  p_env%K_P_P1_env => K_P_P1_env
  
  call get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, &
       dft_control=dft_control)
  call get_atomic_kind_set(atomic_kind_set, nsgf=n_ao, nelectron=n_mo)
  n_mo=n_mo/2-dft_control%charge ! lda closed shell only !
  p_env%n_mo=n_mo
  p_env%n_ao=n_ao

  if (present(I_SP_ptr)) then
     p_env%I_SP => I_SP_ptr
  else
     call allocate_blacs_matrix_vect(p_env%I_SP,1,n_ao,n_ao,&
          name="p_env_I_SP",globenv=global_env,error=error)
  end if
  call allocate_blacs_matrix_vect(p_env%m_epsilon,1,n_mo,n_mo,&
       name="p_env_m_epsilon",globenv=global_env,error=error)
  if (present(rhs_ptr)) then
     p_env%rhs => rhs_ptr
  else
     call allocate_blacs_matrix_vect(p_env%rhs,1,n_ao,n_mo,&
          name="p_env_rhs",globenv=global_env,error=error)
  end if
  call allocate_blacs_matrix(p_env%tmp_aomo_1,n_ao,n_mo,&
       optimal_blacs_row_block_size,&
       optimal_blacs_col_block_size,"p_env_aomo_1",global_env)
  call allocate_blacs_matrix_vect(p_env%r,1,n_ao,n_mo,&
       name="p_env_r",globenv=global_env,error=error)
  if (present(c_ortho_ptr)) then
     p_env%c_ortho => c_ortho_ptr
  else
     call allocate_blacs_matrix_vect(p_env%c_ortho,1,n_ao,n_mo,&
          name="p_env_c_ortho",globenv=global_env,error=error)
  end if
  if (present(psi0d_ptr)) then
     p_env%psi0d => psi0d_ptr
  else
     call allocate_blacs_matrix_vect(p_env%psi0d,1,n_mo,&
          optimal_blacs_row_block_size,&
          name="p_env_psi0d",globenv=global_env,error=error)
  end if
  p_env%s_struct_changed=.true.
  p_env%initialized=.true.
end subroutine p_env_init
!***************************************************************************

!!****f* qs_p_types/p_env_dealloc_ref [1.0] *
!!
!!   NAME
!!     p_env_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory that has been allocated by the p_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the environment to deallocate
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
subroutine p_env_dealloc_ref(p_env,error)
  type(qs_p_env_type), intent(inout) :: p_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='p_env_dealloc_ref',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     nullify(p_env%qs_env, p_env%ks_env, p_env%K_P_P1_env)
     call deallocate_blacs_matrix_vect(p_env%I_SP,error=error)
     call deallocate_blacs_matrix_vect(p_env%m_epsilon,error=error)
     call deallocate_blacs_matrix_vect(p_env%c_ortho,error=error)
     call deallocate_blacs_matrix_vect(p_env%r,error=error)
     call deallocate_blacs_matrix_vect(p_env%rhs,error=error)
     call deallocate_blacs_matrix_vect(p_env%psi0d,error=error)
     call deallocate_blacs_matrix(p_env%tmp_aomo_1)
     call cp_sp_array_dealloc(p_env%Hrho,error=error)
     call cp_sp_array_dealloc(p_env%K_P_P1,error=error)
     call cp_sp_array_dealloc(p_env%P,error=error)
     call cp_sp_array_dealloc(p_env%P1,error=error)
     p_env%initialized=.false.
  end if
end subroutine p_env_dealloc_ref
!***************************************************************************

! setup to do

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
subroutine p_env_eval_l1(p_env,global_env,error)
  type(qs_p_env_type), intent(inout) :: p_env
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='p_env_eval_l1',&
        routineP=moduleN//':'//routineN
  type(real_matrix_p_type), dimension(:),pointer :: s
  failure=.false.
  
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call get_qs_env(p_env%qs_env,s=s)
     call cp_sp_multiplies_blacs(p_env%Hrho(1)%matrix,&
          p_env%c_ortho(1)%blacs_matrix,p_env%r(1)%blacs_matrix)
     call blacs_symm('R','U',p_env%n_ao,p_env%n_mo,1.0_wp,&
          p_env%m_epsilon(1)%blacs_matrix,&
          p_env%c_ortho(1)%blacs_matrix,0.0_wp,p_env%tmp_aomo_1,&
          global_env)
     call sparse_times_blacs(s(1)%matrix,p_env%tmp_aomo_1, &
          p_env%r(1)%blacs_matrix,p_env%n_mo,&
          mpi_env=cp_mpi_get_from_glob(global_env),alpha=1.0_wp,beta=1.0_wp)
  end if
end subroutine p_env_eval_l1
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
subroutine p_env_eval_l2(p_env, target_m, C, global_env, alpha, beta, error)
  type(qs_p_env_type), intent(inout) :: p_env
  type(global_environment_type), intent(in), target :: global_env
  type(blacs_matrix_p_type), dimension(:),intent(in) :: C
  type(blacs_matrix_p_type), dimension(:),intent(inout) :: target_m
  real(kind=wp), intent(in), optional :: alpha,beta
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='p_env_eval_l2',&
        routineP=moduleN//':'//routineN
  integer :: i
  failure=.false.
  
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     do i=1,size(p_env%P1)
        call cp_sp_set_to(p_env%P1(i)%matrix,value=0.0_wp,error=error)
     end do
     call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%matrix,&
          p_env%c_ortho(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix,&
          ncol=p_env%n_mo, globenv=global_env)
     call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%matrix,&
          p_env%psi0d(1)%blacs_matrix,p_env%c_ortho(1)%blacs_matrix,&
          ncol=p_env%n_mo, globenv=global_env)
     call kpp1_calculate(p_env%K_P_P1_env,p1=p_env%P1,&
          kpp1_matrix=p_env%K_P_P1,&
          global_env=global_env,error=error)
     
     call cp_sp_multiplies_blacs(p_env%K_P_P1(1)%matrix,C(1)%blacs_matrix,&
          target_m(1)%blacs_matrix, alpha=alpha,&
          beta=beta)
  end if
end subroutine p_env_eval_l2
!***************************************************************************

!!****f* qs_p_types/p_env_preortho [1.0] *
!!
!!   NAME
!!     p_env_preortho
!!
!!   FUNCTION
!!     does a preorthogonalization of the given matrix:
!!     target_m=(I-PS)C
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment
!!     - target_m: place where to store the result
!!     - C: matrix to orthogonalize
!!     - alpha: scale factor of the result (defaults to 1.0)
!!     - beta: scale factor of the old values (defaults to 0.0)
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
subroutine p_env_preortho(p_env,target_m, C, global_env, &
     alpha,beta, n_cols,error)
  type(qs_p_env_type), intent(inout) :: p_env
  type(global_environment_type), intent(in), target :: global_env
  type(blacs_matrix_p_type), dimension(:),intent(in) :: C
  type(blacs_matrix_p_type), dimension(:),intent(inout) :: target_m
  real(kind=wp), intent(in), optional :: alpha,beta
  integer, intent(in), optional :: n_cols
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='p_env_preortho',&
        routineP=moduleN//':'//routineN
  integer :: nrow_global, ncol_global, cols
  real(kind=wp) :: al,be
  failure=.false.
  al=1.0_wp; be=0.0_wp
  
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  call get_blacs_matrix_info(C(1)%blacs_matrix,nrow_global=nrow_global, &
       ncol_global=ncol_global)
  CPPrecondition(nrow_global>=p_env%n_ao,cp_failure_level,routineP,error,failure)
  cols=ncol_global
  if (present(n_cols)) then
     CPPrecondition(n_cols<=cols,cp_failure_level,routineP,error,failure)
     cols=n_cols
  end if
  ! check target_m?
  if (present(alpha)) al=alpha
  if (present(beta)) be=beta
  if (.not. failure) then     
     call blacs_gemm('T','N',p_env%n_ao,cols,p_env%n_ao,al,&
          p_env%I_SP(1)%blacs_matrix,C(1)%blacs_matrix,&
          be,target_m(1)%blacs_matrix,global_env)
  end if
end subroutine p_env_preortho
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
!!     - target_m: place where to store the result
!!     - C: matrix to orthogonalize
!!     - alpha: scale factor of the result (defaults to 1.0)
!!     - beta: scale factor of the old values (defaults to 0.0)
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
subroutine p_env_postortho(p_env,target_m, C, global_env,&
     alpha,beta, n_cols,error)
  type(qs_p_env_type), intent(inout) :: p_env
  type(global_environment_type), intent(in), target :: global_env
  type(blacs_matrix_p_type), dimension(:),intent(in) :: C
  type(blacs_matrix_p_type), dimension(:),intent(inout) :: target_m
  real(kind=wp), intent(in), optional :: alpha,beta
  integer, intent(in), optional :: n_cols
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='p_env_postortho',&
        routineP=moduleN//':'//routineN
  integer :: nrow_global, ncol_global, cols
  real(kind=wp) :: al,be
  failure=.false.
  al=1.0_wp; be=0.0_wp
  
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  call get_blacs_matrix_info(C(1)%blacs_matrix,nrow_global=nrow_global,&
       ncol_global=ncol_global)
  CPPrecondition(nrow_global>=p_env%n_ao,cp_failure_level,routineP,error,failure)
  cols=ncol_global
  if (present(n_cols)) then
     CPPrecondition(n_cols<=cols,cp_failure_level,routineP,error,failure)
     cols=n_cols
  end if
  ! check target_m?
  if (present(alpha)) al=alpha
  if (present(beta)) be=beta
  if (.not. failure) then     
     call blacs_gemm('N','N',p_env%n_ao,cols,p_env%n_ao,al,&
          p_env%I_SP(1)%blacs_matrix,C(1)%blacs_matrix,&
          be,target_m(1)%blacs_matrix,global_env)
  end if
end subroutine p_env_postortho
!***************************************************************************

end module qs_p_types
