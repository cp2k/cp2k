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
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
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
  private
  public :: qs_p_env_type, p_env_init, qs_K_P_P1_env_type


  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_p_types'

!***
!****************************************************************************

  !!****s* qs_p_build_kernel/qs_K_P_P1_env_type [1.0] *
  !!
  !!   NAME
  !!     qs_K_P_P1_env_type
  !!
  !!   FUNCTION
  !!     environement that keeps the informations and temporary
  !!     val to build the K_P_P1 kernel matrix
  !!
  !!   NOTES
  !!     share grids with ks_env?
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     07.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type qs_K_P_P1_env_type
     integer :: ref_count
     type(cp_block_matrix_p_type), dimension(:), pointer :: K_P_P1
  end type qs_K_P_P1_env_type
  !!***
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
!!     for the moment no smearing of the orbitals.
!!     I_SP could be replaced by keeping S_orb (ao x ao) and 
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     7.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  type qs_p_env_type
     logical :: initialized, s_mstruct_changed, orthogonal_orbitals
     type(cp_block_matrix_p_type), dimension(:), pointer :: &
          K_P_P1, Hrho, P, P1
     type(blacs_matrix_p_type), dimension(:), pointer :: m_epsilon,&
          psi0d, S_psi0, c_ortho, r
     type(blacs_matrix_type), pointer :: tmp_aomo_1
     type(qs_ks_env_type), pointer :: ks_env
     type(qs_K_P_P1_env_type), pointer :: K_P_P1_env
     type(qs_environment_type), pointer :: qs_env
     integer n_mo,n_ao
  end type qs_p_env_type
!!FM  type qs_p_env_type
!!FM     logical :: initialized, s_mstruct_changed, orthogonal_orbitals
!!FM     type(cp_block_matrix_p_type), dimension(:), pointer :: &
!!FM          K_P_P1
!!FM     type(blacs_matrix_p_type), dimension(:), pointer :: m_epsilon,&
!!FM          psi0d, S_psi0
!!FM     type(qs_ks_env_type), pointer :: ks_env
!!FM     type(qs_K_P_P1_env_type), pointer :: K_P_P1_env
!!FM     type(qs_environment_type), pointer :: qs_env
!!FM     integer n_mo,n_ao
!!FM  end type qs_p_env_type
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
!!     - orthogonal_orbitals: if the orbitals are 
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
       c_ortho_ptr, psi0d_ptr, orthogonal_orbitals, error)
    type(qs_p_env_type), intent(out) :: p_env
    type(qs_environment_type), intent(in), target :: qs_env
    type(qs_ks_env_type), intent(in), target :: ks_env
    type(qs_K_P_P1_env_type), intent(in), target :: K_P_P1_env
    type(global_environment_type), intent(in), target :: global_env
    type(blacs_matrix_p_type), dimension(:),pointer, optional :: c_ortho_ptr,&
         psi0d_ptr
    logical, intent(in), optional :: orthogonal_orbitals
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='p_env_init',&
         routineP=moduleN//':'//routineN
    integer :: n_ao, n_mo
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    type(dft_control_type), pointer :: dft_control
    failure=.false.
    nullify(atomic_kind_set, dft_control)

    nullify(p_env%psi0d, p_env%c_ortho, p_env%r,&
         p_env%tmp_aomo_1, p_env%S_psi0)
    nullify(p_env%m_epsilon)
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
    p_env%orthogonal_orbitals=.false.
    if (present(orthogonal_orbitals)) p_env%orthogonal_orbitals=orthogonal_orbitals

    call allocate_blacs_matrix_vect(p_env%S_psi0,1,n_ao,n_mo,&
         name="p_env_S_psi0",globenv=global_env,error=error)
    call allocate_blacs_matrix_vect(p_env%m_epsilon,1,n_mo,n_mo,&
         name="p_env_m_epsilon",globenv=global_env,error=error)
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
    p_env%s_mstruct_changed=.true.
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
       call deallocate_blacs_matrix_vect(p_env%S_psi0,error=error)
       call deallocate_blacs_matrix_vect(p_env%m_epsilon,error=error)
       call deallocate_blacs_matrix_vect(p_env%c_ortho,error=error)
       call deallocate_blacs_matrix_vect(p_env%r,error=error)
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

!!****f* qs_p_types/p_env_psi0_changed [1.0] *
!!
!!   NAME
!!     p_env_psi0_changed
!!
!!   FUNCTION
!!     to be called after the value of psi0 has changed.
!!     p_env%r should contain psi0, then psi0d, S_psi0, the filtred P,
!!     Hrho, m_epsilon are calculated.
!!
!!   NOTES
!!     Hrho and P are just wrapper of the k and P matrixes in qs_env
!!     (I want to use the higher level block_matrix)
!!
!!   INPUTS
!!     - p_env: the perturbation environment to set
!!     - global_env: the global environment for the calculations
!!     - Hrho_psi0d: is given, then the partial result Hrho_psi0d is stored in
!!       that vector
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
  subroutine p_env_psi0_changed(p_env, global_env, Hrho_psi0d, error)
    type(qs_p_env_type), intent(inout) :: p_env
    type(global_environment_type), intent(in), target :: global_env
    type(blacs_matrix_p_type), dimension(:), intent(inout), optional :: Hrho_psi0d
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='p_env_set_psi0',&
         routineP=moduleN//':'//routineN
    type(cp_b_matrix_struct_type), pointer :: s_mstruct
    type(real_matrix_p_type), dimension(:), pointer :: S,P,K
    failure=.false.
    nullify(S,P,K,s_mstruct)

    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (p_env%s_mstruct_changed) then
          call cp_sp_array_dealloc(p_env%Hrho, error=error)
          call cp_sp_array_dealloc(p_env%P, error=error)
          call cp_sp_array_dealloc(p_env%P1, error=error)
          call cp_sp_array_dealloc(p_env%K_P_P1,error=error)

          call get_qs_env(p_env%qs_env,s_mstruct=s_mstruct,K=K,P=P)
          allocate(p_env%Hrho(1),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          nullify(p_env%Hrho(1)%matrix)
          call cp_sp_create(p_env%Hrho(1)%matrix, matrix_struct=s_mstruct,&
               u_matrix=k(1)%matrix, should_dealloc_matrix=.false.,&
               error=error)
          allocate(p_env%P(1),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          nullify(p_env%P(1)%matrix)
          call cp_sp_create(p_env%P(1)%matrix, matrix_struct=s_mstruct,&
               u_matrix=P(1)%matrix, should_dealloc_matrix=.false.,&
               error=error)
          call cp_sp_array_create(p_env%P1,1,s_mstruct,&
               allocate_nonsparse_blocks=.true., error=error)
          call cp_sp_array_create(p_env%K_P_P1,1,s_mstruct,&
               allocate_nonsparse_blocks=.true., error=error)
          p_env%s_mstruct_changed=.false.
       end if

       if (p_env%orthogonal_orbitals) then
          call copy_blacs_to_blacs_matrix(p_env%r(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix)
       else
          call blacs_syrk('U','N',p_env%n_mo,1.0_wp,p_env%r(1)%blacs_matrix,1,1,&
               0.0_wp,p_env%m_epsilon(1)%blacs_matrix)
          call blacs_cholesky_decompose(p_env%m_epsilon(1)%blacs_matrix)
          call blacs_triangular_multiply(p_env%m_epsilon(1)%blacs_matrix,&
               p_env%r(1)%blacs_matrix,&
               transpose_tr=.true.,invert_tr=.true.)
          call blacs_triangular_multiply(p_env%m_epsilon(1)%blacs_matrix,&
               p_env%r(1)%blacs_matrix,&
               transpose_tr=.false.,invert_tr=.true.)
       end if
       call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,p_env%r(1)%blacs_matrix,&
            p_env%psi0d(1)%blacs_matrix,p_env%n_mo)

       call qs_ks_did_change(p_env%ks_env,qs_env=p_env%qs_env,&
            global_env=global_env,p_changed=.true.,error=error)

       call get_qs_env(p_env%qs_env, s=S)
       call sparse_times_blacs(S(1)%matrix,p_env%r(1)%blacs_matrix,&
            p_env%S_psi0(1)%blacs_matrix,&
            p_env%n_mo,para_env=global_env%para_env)

       call qs_ks_update_qs_env(p_env%ks_env,qs_env=p_env%qs_env,&
            global_env=global_env)

       call cp_sp_multiplies_blacs(p_env%Hrho(1)%matrix,&
            p_env%psi0d(1)%blacs_matrix,&
            p_env%tmp_aomo_1, error=error)
       call blacs_gemm('T','N',p_env%n_mo,p_env%n_mo,p_env%n_ao,&
            1.0_wp,p_env%psi0d(1)%blacs_matrix,p_env%tmp_aomo_1,&
            0.0_wp,p_env%m_epsilon(1)%blacs_matrix)
       call blacs_scale_and_d(p_env%m_epsilon(1)%blacs_matrix,&
            scale=-1.0_wp)
    end if
  end subroutine p_env_psi0_changed
!***************************************************************************

end module qs_p_types
