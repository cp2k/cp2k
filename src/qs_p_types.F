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
  USE kinds, ONLY: wp=>dp, dbl
  USE global_types, ONLY: global_environment_type
  USE timings, ONLY: timeset, timestop
  USE qs_blacs, ONLY: allocate_blacs_matrix, deallocate_blacs_matrix,&
       blacs_gemm, blacs_add, sparse_times_blacs, allocate_blacs_matrix_vect,&
       deallocate_blacs_matrix_vect, optimal_blacs_row_block_size, &
       optimal_blacs_col_block_size, blacs_matrix_p_type, blacs_matrix_type,&
       sparse_plus_blacs_blacst, blacs_symm, get_blacs_matrix_info,&
       blacs_scale_and_d, copy_blacs_to_sparse_matrix,&
       blacs_cholesky_decompose, blacs_cholesky_restore, blacs_syrk,&
       copy_blacs_to_blacs_matrix,blacs_triangular_multiply, &
       replicate_blacs_matrix
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
!!     - glob_env: the global environment
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
  SUBROUTINE p_env_init(p_env, qs_env, ks_env, K_P_P1_env, glob_env,&
                        psi0d, orthogonal_orbitals, error)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(out)                          :: p_env
    TYPE(qs_environment_type), INTENT(in), TARGET             :: qs_env
    TYPE(qs_ks_env_type), INTENT(in), TARGET                  :: ks_env
    TYPE(qs_K_P_P1_env_type), INTENT(in), TARGET              :: K_P_P1_env
    TYPE(global_environment_type), INTENT(in), TARGET         :: glob_env
    TYPE(blacs_matrix_p_type), DIMENSION(:),POINTER, OPTIONAL :: psi0d
    LOGICAL, INTENT(in), OPTIONAL                             :: orthogonal_orbitals
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)              :: error

    ! locals
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER      :: routineN='p_env_init',&
                                        routineP=moduleN//':'//routineN
    INTEGER                          :: n_ao, n_mo, n_spins, spin, istat
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(dft_control_type), POINTER  :: dft_control
    TYPE(blacs_matrix_type), POINTER :: qs_env_c

    ! code

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

    n_spins = p_env%qs_env%n_spins
    DO spin=1, qs_env%n_spins

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
                                    globenv=glob_env, error=error)
    CALL allocate_blacs_matrix_vect(p_env%m_epsilon, n_spins, n_mo, n_mo,&
                                    name="p_env_m_epsilon", &
                                    globenv=glob_env, error=error)

    IF (PRESENT(psi0d)) THEN
       p_env%psi0d => psi0d
    ELSE
       ! what is that??
       CALL allocate_blacs_matrix_vect(p_env%psi0d, n_spins, n_mo, &
                                       optimal_blacs_row_block_size,&
                                       name="p_env_psi0d", &
                                       globenv=glob_env, error=error)
    END IF

    ! do still we need this?
    p_env%s_mstruct_changed=.TRUE.

    p_env%initialized=.TRUE.

    IF (PRESENT(psi0d)) CALL p_env_psi0_changed(p_env, glob_env, error=error)
    
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
  SUBROUTINE p_env_dealloc_ref(p_env, error)

    ! arguments
    TYPE(qs_p_env_type), INTENT(inout)           :: p_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL                     :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_dealloc_ref',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.

    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
       NULLIFY(p_env%K_P_P1_env)
       !! call qs_K_P_P1_env_dealloc_ref(p_env%K_P_P1_env, error)
       NULLIFY(p_env%qs_env, p_env%ks_env)
       CALL deallocate_blacs_matrix_vect(p_env%S_psi0, error=error)
       CALL deallocate_blacs_matrix_vect(p_env%m_epsilon, error=error)
       CALL deallocate_blacs_matrix_vect(p_env%psi0d, error=error)
       CALL cp_sp_array_dealloc(p_env%K_P_P1, error=error)
       p_env%initialized = .FALSE.
    END IF
  END SUBROUTINE p_env_dealloc_ref
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
!!     - glob_env: the global environment for the calculations
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

  SUBROUTINE p_env_psi0_changed(p_env, glob_env, Hrho_psi0d, error)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(inout)                :: p_env
    TYPE(global_environment_type), INTENT(in), TARGET :: glob_env
    TYPE(blacs_matrix_p_type), DIMENSION(:), INTENT(inout), OPTIONAL :: Hrho_psi0d
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)      :: error

    ! locals
    LOGICAL                                :: failure
    INTEGER                                :: stat
    CHARACTER(len=*), PARAMETER            :: routineN = 'p_env_psi0_changed', &
                                              routineP = moduleN//':'//routineN
    TYPE(blacs_matrix_type), POINTER       :: tmp
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: S
    INTEGER                                :: n_spins, spin

!    TYPE(cp_b_matrix_struct_type), POINTER :: s_mstruct
!    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: S, P, K

    ! code
    failure=.FALSE.

    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)

!    NULLIFY(S,P,K,s_mstruct)

    NULLIFY(tmp)

    IF (.NOT. failure) THEN
       
       n_spins = p_env%qs_env%n_spins
       DO spin=1, n_spins

          CALL replicate_blacs_matrix(p_env%psi0d(spin)%blacs_matrix, tmp, &
                                      routineP//": temporary matrix")

       ! still needed?
!!TC       if (p_env%s_mstruct_changed) then
!!TC          call cp_sp_array_dealloc(p_env%Hrho, error=error)
!!TC          call cp_sp_array_dealloc(p_env%P, error=error)
!!TC          call cp_sp_array_dealloc(p_env%P1, error=error)
!!TC          call cp_sp_array_dealloc(p_env%K_P_P1,error=error)
!!TC
!!TC          call get_qs_env(p_env%qs_env,s_mstruct=s_mstruct,K=K,P=P)
!!TC          allocate(p_env%Hrho(1),stat=stat)
!!TC          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!TC          nullify(p_env%Hrho(1)%matrix)
!!TC          call cp_sp_create(p_env%Hrho(1)%matrix, matrix_struct=s_mstruct,&
!!TC               u_matrix=k(1)%matrix, should_dealloc_matrix=.false.,&
!!TC               error=error)
!!TC          allocate(p_env%P(1),stat=stat)
!!TC          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!TC          nullify(p_env%P(1)%matrix)
!!TC          call cp_sp_create(p_env%P(1)%matrix, matrix_struct=s_mstruct,&
!!TC               u_matrix=P(1)%matrix, should_dealloc_matrix=.false.,&
!!TC               error=error)
!!TC          call cp_sp_array_create(p_env%P1,1,s_mstruct,&
!!TC               allocate_nonsparse_blocks=.true., error=error)
!!TC          call cp_sp_array_create(p_env%K_P_P1,1,s_mstruct,&
!!TC               allocate_nonsparse_blocks=.true., error=error)
!!TC          p_env%s_mstruct_changed=.false.
!!TC       end if


          ! caluclates 
          ! m_eplsilon = - psi0d^T times K times psi0d
          !            = - [K times psi0d]^T times psi0d (because K is symmetric)
          IF (p_env%orthogonal_orbitals) THEN
             
             ! tmp = k times psi0d
             CALL sparse_times_blacs(p_env%qs_env%k(spin)%matrix, &
                  p_env%psi0d(spin)%blacs_matrix, &
                  tmp, p_env%n_ao(spin), &
                  glob_env%para_env)
             
             ! m_epsilon = -1 times tmp^T times psi0d
             CALL blacs_gemm('T', 'N', &
                  p_env%n_mo(spin), p_env%n_mo(spin), p_env%n_ao(spin), &
                  -1.0_dbl, tmp, p_env%psi0d(spin)%blacs_matrix, &
                  0.0_dbl, p_env%m_epsilon(spin)%blacs_matrix)

             ! @fawzi: I don't know what the following does. Please
             !         check it and change according to your needs.

!!TC          call copy_blacs_to_blacs_matrix(p_env%r(1)%blacs_matrix,p_env%psi0d(1)%blacs_matrix)
!!TC       else
!!TC          call blacs_syrk('U','N',p_env%n_mo,1.0_wp,p_env%r(1)%blacs_matrix,1,1,&
!!TC               0.0_wp,p_env%m_epsilon(1)%blacs_matrix)
!!TC          call blacs_cholesky_decompose(p_env%m_epsilon(1)%blacs_matrix)
!!TC          call blacs_triangular_multiply(p_env%m_epsilon(1)%blacs_matrix,&
!!TC               p_env%r(1)%blacs_matrix,&
!!TC               transpose_tr=.true.,invert_tr=.true.)
!!TC          call blacs_triangular_multiply(p_env%m_epsilon(1)%blacs_matrix,&
!!TC               p_env%r(1)%blacs_matrix,&
!!TC               transpose_tr=.false.,invert_tr=.true.)
          END IF

!!TC       call sparse_plus_blacs_blacst(p_env%P1(1)%matrix%sm,p_env%r(1)%blacs_matrix,&
!!TC            p_env%psi0d(1)%blacs_matrix,p_env%n_mo)
!!TC
!!TC       call qs_ks_did_change(p_env%ks_env,qs_env=p_env%qs_env,&
!!TC            global_env=glob_env,p_changed=.true.,error=error)
!!TC
          CALL get_qs_env(p_env%qs_env, s=S)
          CALL sparse_times_blacs(S(spin)%matrix, &
               p_env%psi0d(spin)%blacs_matrix, &
               p_env%S_psi0(spin)%blacs_matrix, &
               p_env%n_mo(spin), para_env=glob_env%para_env)
       
!!TC          CALL sparse_times_blacs(S(spin)%matrix, p_env%r(1)%blacs_matrix,&
!!TC               p_env%S_psi0(1)%blacs_matrix,&
!!TC               p_env%n_mo,para_env=glob_env%para_env)
!!TC
!!TC       call qs_ks_update_qs_env(p_env%ks_env,qs_env=p_env%qs_env,&
!!TC            global_env=glob_env)
!!TC
!!TC       call cp_sp_multiplies_blacs(p_env%Hrho(1)%matrix,&
!!TC            p_env%psi0d(1)%blacs_matrix,&
!!TC            p_env%tmp_aomo_1, error=error)
!!TC       call blacs_gemm('T','N',p_env%n_mo,p_env%n_mo,p_env%n_ao,&
!!TC            1.0_wp,p_env%psi0d(1)%blacs_matrix,p_env%tmp_aomo_1,&
!!TC            0.0_wp,p_env%m_epsilon(1)%blacs_matrix)
!!TC       call blacs_scale_and_d(p_env%m_epsilon(1)%blacs_matrix,&
!!TC            scale=-1.0_wp)

          CALL deallocate_blacs_matrix(tmp); NULLIFY(tmp)

       END DO
    END IF


  END SUBROUTINE p_env_psi0_changed
!***************************************************************************

END MODULE qs_p_types
