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
MODULE qs_p_types
  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
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
  USE qs_mo_types, ONLY : get_mo_set
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
  PRIVATE
  PUBLIC :: qs_p_env_type, p_env_init, qs_K_P_P1_env_type


  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_types'

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
  TYPE qs_K_P_P1_env_type
     INTEGER :: ref_count
     TYPE(cp_block_matrix_p_type), DIMENSION(:), POINTER :: K_P_P1
  END TYPE qs_K_P_P1_env_type
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

!!TC  type qs_p_env_type
!!TC     logical :: initialized, s_mstruct_changed, orthogonal_orbitals
!!TC     type(cp_block_matrix_p_type), dimension(:), pointer :: &
!!TC          K_P_P1, Hrho, P, P1
!!TC     type(blacs_matrix_p_type), dimension(:), pointer :: m_epsilon,&
!!TC          psi0d, S_psi0, c_ortho, r
!!TC     type(blacs_matrix_type), pointer :: tmp_aomo_1
!!TC     type(qs_ks_env_type), pointer :: ks_env
!!TC     type(qs_K_P_P1_env_type), pointer :: K_P_P1_env
!!TC     type(qs_environment_type), pointer :: qs_env
!!TC     integer n_mo,n_ao
!!TC  end type qs_p_env_type

  TYPE qs_p_env_type
     LOGICAL :: initialized, s_mstruct_changed, orthogonal_orbitals
     TYPE(cp_block_matrix_p_type), DIMENSION(:), POINTER :: &
          K_P_P1
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: m_epsilon,&
          psi0d, S_psi0
     TYPE(qs_ks_env_type), POINTER :: ks_env
     TYPE(qs_K_P_P1_env_type), POINTER :: K_P_P1_env
     TYPE(qs_environment_type), POINTER :: qs_env
     INTEGER, DIMENSION(2) :: n_mo, & ! no of molecular orbitals
                              n_ao    ! no of basis functions

  END TYPE qs_p_env_type

!!***
!****************************************************************************

CONTAINS

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
  SUBROUTINE p_env_init(p_env, qs_env, ks_env, K_P_P1_env, global_env,&
                        psi0d, orthogonal_orbitals, error)
    ! arguments
    TYPE(qs_p_env_type), INTENT(out) :: p_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(qs_ks_env_type), INTENT(in), TARGET :: ks_env
    TYPE(qs_K_P_P1_env_type), INTENT(in), TARGET :: K_P_P1_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(blacs_matrix_p_type), DIMENSION(:),POINTER, OPTIONAL :: psi0d
    LOGICAL, INTENT(in), OPTIONAL :: orthogonal_orbitals
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='p_env_init',&
         routineP=moduleN//':'//routineN
    INTEGER :: n_ao, n_mo, n_spins, spin, istat
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(dft_control_type), POINTER :: dft_control
    ! needed because qs_mo_set doesn't provide accessor functions for
    ! its private members
    TYPE(blacs_matrix_type), POINTER :: qs_env_c

    ! code

    ! maybe it would be good to store this in the qs_env structure
    n_spins = SIZE(p_env%qs_env%k) 

    failure=.FALSE.

    NULLIFY(atomic_kind_set, dft_control)
    NULLIFY(p_env%K_P_P1, &
            p_env%m_epsilon, &
            p_env%psi0d, &
            p_env%S_psi0, &
            p_env%ks_env, &
            p_env%K_P_P1_env, &
            p_env%qs_env)

    p_env%qs_env => qs_env
    p_env%ks_env => ks_env
    p_env%K_P_P1_env => K_P_P1_env

    DO spin=1, n_spins

       IF (PRESENT(psi0d)) THEN
          CALL get_blacs_matrix_info(psi0d(spin)%blacs_matrix, &
                                     ncol_global=n_mo, nrow_global=n_ao)
       ELSE
          CALL get_mo_set(p_env%qs_env%c(spin)%mo_set, eigenvectors=qs_env_c)
          CALL get_blacs_matrix_info(qs_env_c, &
                                     ncol_global=n_mo, nrow_global=n_ao)
       END IF
       p_env%n_mo(spin) = n_mo
       p_env%n_ao(spin) = n_ao

!!TC       call get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, &
!!TC                               dft_control=dft_control)
!!TC       call get_atomic_kind_set(atomic_kind_set, nsgf=n_ao, nelectron=n_mo)
!!TC       n_mo=n_mo/2-dft_control%charge ! lda closed shell only !
!!TC       p_env%n_mo=n_mo
!!TC       p_env%n_ao=n_ao

       p_env%orthogonal_orbitals=.FALSE.
       IF (PRESENT(orthogonal_orbitals)) &
            p_env%orthogonal_orbitals=orthogonal_orbitals

    END DO

    CALL allocate_blacs_matrix_vect(p_env%S_psi0, n_spins, n_ao, n_mo,&
                                    name="p_env_S_psi0", &
                                    globenv=global_env, error=error)
    CALL allocate_blacs_matrix_vect(p_env%m_epsilon, n_spins, n_mo, n_mo,&
                                    name="p_env_m_epsilon", &
                                    globenv=global_env, error=error)

    IF (PRESENT(psi0d)) THEN
       p_env%psi0d => psi0d
    ELSE
       CALL allocate_blacs_matrix_vect(p_env%psi0d, n_spins, n_mo, &
                                       optimal_blacs_row_block_size,&
                                       name="p_env_psi0d", &
                                       globenv=global_env, error=error)
    END IF

    ! do we need this anymore?
    p_env%s_mstruct_changed=.TRUE.

    !! call p_env_psi0_changed(p_env, ...)

    p_env%initialized=.TRUE.

  END SUBROUTINE p_env_init
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
!!!!  subroutine p_env_dealloc_ref(p_env,error)
!!!!    type(qs_p_env_type), intent(inout) :: p_env
!!!!    type(cp_error_type), optional, intent(inout) :: error
!!!!
!!!!    logical :: failure
!!!!    character(len=*), parameter :: routineN='p_env_dealloc_ref',&
!!!!         routineP=moduleN//':'//routineN
!!!!    failure=.false.
!!!!
!!!!    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
!!!!    if (.not. failure) then
!!!!       nullify(p_env%qs_env, p_env%ks_env, p_env%K_P_P1_env)
!!!!       call deallocate_blacs_matrix_vect(p_env%S_psi0,error=error)
!!!!       call deallocate_blacs_matrix_vect(p_env%m_epsilon,error=error)
!!!!       call deallocate_blacs_matrix_vect(p_env%c_ortho,error=error)
!!!!       call deallocate_blacs_matrix_vect(p_env%r,error=error)
!!!!       call deallocate_blacs_matrix_vect(p_env%psi0d,error=error)
!!!!       call deallocate_blacs_matrix(p_env%tmp_aomo_1)
!!!!       call cp_sp_array_dealloc(p_env%Hrho,error=error)
!!!!       call cp_sp_array_dealloc(p_env%K_P_P1,error=error)
!!!!       call cp_sp_array_dealloc(p_env%P,error=error)
!!!!       call cp_sp_array_dealloc(p_env%P1,error=error)
!!!!       p_env%initialized=.false.
!!!!    end if
!!!!  end subroutine p_env_dealloc_ref
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
!!!!  subroutine p_env_psi0_changed(p_env, global_env, Hrho_psi0d, error)
!!!!    type(qs_p_env_type), intent(inout) :: p_env
!!!!    type(global_environment_type), intent(in), target :: global_env
!!!!    type(blacs_matrix_p_type), dimension(:), intent(inout), optional :: Hrho_psi0d
!!!!    type(cp_error_type), optional, intent(inout) :: error
!!!!
!!!!    logical :: failure
!!!!    integer :: stat
!!!!    character(len=*), parameter :: routineN='p_env_set_psi0',&
!!!!         routineP=moduleN//':'//routineN
!!!!    type(cp_b_matrix_struct_type), pointer :: s_mstruct
!!!!    type(real_matrix_p_type), dimension(:), pointer :: S,P,K
!!!!    failure=.false.
!!!!    nullify(S,P,K,s_mstruct)
!!!!
!!!!    CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
!!!!    if (.not. failure) then
!!!!       if (p_env%s_mstruct_changed) then
!!!!          call cp_sp_array_dealloc(p_env%Hrho, error=error)
!!!!          call cp_sp_array_dealloc(p_env%P, error=error)
!!!!          call cp_sp_array_dealloc(p_env%P1, error=error)
!!!!          call cp_sp_array_dealloc(p_env%K_P_P1,error=error)
!!!!
!!!!          call get_qs_env(p_env%qs_env,s_mstruct=s_mstruct,K=K,P=P)
!!!!          allocate(p_env%Hrho(1),stat=stat)
!!!!          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!!!          nullify(p_env%Hrho(1)%matrix)
!!!!          call cp_sp_create(p_env%Hrho(1)%matrix, matrix_struct=s_mstruct,&
!!!!               u_matrix=k(1)%matrix, should_dealloc_matrix=.false.,&
!!!!               error=error)
!!!!          allocate(p_env%P(1),stat=stat)
!!!!          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!!!          nullify(p_env%P(1)%matrix)
!!!!          call cp_sp_create(p_env%P(1)%matrix, matrix_struct=s_mstruct,&
!!!!               u_matrix=P(1)%matrix, should_dealloc_matrix=.false.,&
!!!!               error=error)
!!!!          call cp_sp_array_create(p_env%P1,1,s_mstruct,&
!!!!               allocate_nonsparse_blocks=.true., error=error)
!!!!          call cp_sp_array_create(p_env%K_P_P1,1,s_mstruct,&
!!!!               allocate_nonsparse_blocks=.true., error=error)
!!!!          p_env%s_mstruct_changed=.false.
!!!!       end if
!!!!
!!!!       if (p_env%orthogonal_orbitals) then
!!!!          call copy_blacs_to_blacs_matrix(p_env%r(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix)
!!!!       else
!!!!          call blacs_syrk('U','N',p_env%n_mo,1.0_wp,p_env%r(1)%blacs_matrix,1,1,&
!!!!               0.0_wp,p_env%m_epsilon(1)%blacs_matrix)
!!!!          call blacs_cholesky_decompose(p_env%m_epsilon(1)%blacs_matrix)
!!!!          call blacs_triangular_multiply(p_env%m_epsilon(1)%blacs_matrix,&
!!!!               p_env%r(1)%blacs_matrix,&
!!!!               transpose_tr=.true.,invert_tr=.true.)
!!!!          call blacs_triangular_multiply(p_env%m_epsilon(1)%blacs_matrix,&
!!!!               p_env%r(1)%blacs_matrix,&
!!!!               transpose_tr=.false.,invert_tr=.true.)
!!!!       end if
!!!!       call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,p_env%r(1)%blacs_matrix,&
!!!!            p_env%psi0d(1)%blacs_matrix,p_env%n_mo)
!!!!
!!!!       call qs_ks_did_change(p_env%ks_env,qs_env=p_env%qs_env,&
!!!!            global_env=global_env,p_changed=.true.,error=error)
!!!!
!!!!       call get_qs_env(p_env%qs_env, s=S)
!!!!       call sparse_times_blacs(S(1)%matrix,p_env%r(1)%blacs_matrix,&
!!!!            p_env%S_psi0(1)%blacs_matrix,&
!!!!            p_env%n_mo,para_env=global_env%para_env)
!!!!
!!!!       call qs_ks_update_qs_env(p_env%ks_env,qs_env=p_env%qs_env,&
!!!!            global_env=global_env)
!!!!
!!!!       call cp_sp_multiplies_blacs(p_env%Hrho(1)%matrix,&
!!!!            p_env%psi0d(1)%blacs_matrix,&
!!!!            p_env%tmp_aomo_1, error=error)
!!!!       call blacs_gemm('T','N',p_env%n_mo,p_env%n_mo,p_env%n_ao,&
!!!!            1.0_wp,p_env%psi0d(1)%blacs_matrix,p_env%tmp_aomo_1,&
!!!!            0.0_wp,p_env%m_epsilon(1)%blacs_matrix)
!!!!       call blacs_scale_and_d(p_env%m_epsilon(1)%blacs_matrix,&
!!!!            scale=-1.0_wp)
!!!!    end if
!!!!  end subroutine p_env_psi0_changed
!***************************************************************************

END MODULE qs_p_types
