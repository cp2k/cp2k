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
!!     basis types for the calculation of the perturbation of density theory.
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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_type
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_sp_array_create,&
                                             cp_sp_array_dealloc,&
                                             cp_sp_create,&
                                             cp_sp_multiplies_blacs,&
                                             cp_sp_set_to
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_multiply
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_matrix,&
                                             fm_pool_get_mstruct,&
                                             fm_pool_give_back_matrix,&
                                             fm_pools_create_matrix_vect,&
                                             fm_pools_give_back_matrix_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                             cp_fm_struct_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_copy,&
                                             cp_fm_vect_create2,&
                                             cp_fm_vect_dealloc
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             wp=>dp
  USE qs_blacs,                        ONLY: blacs_scale_and_d,&
                                             copy_blacs_to_blacs_matrix,&
                                             copy_blacs_to_sparse_matrix,&
                                             cp_fm_add,&
                                             cp_fm_create2,&
                                             cp_fm_gemm,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_symm,&
                                             cp_fm_syrk,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t,&
                                             optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size,&
                                             replicate_blacs_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE sparse_matrix_types,             ONLY: get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
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
!!     -
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (how may objects are sharing this one)
!!     - v_rspace: potential in r space. This is used to do an update only
!!       of what has changed. Useful if P1 converges to some density
!!       (you spare in the grid-ao conversion).
!!       With LSD contains the various orbitals.
!!     - v_ao: the potential in the ao basis (used togheter with v_rspace
!!       to update only what changed
!!     - rebuild_each: how often a full rebuild should be performed
!!     - updates_since_rebuild: how many updates have been done since the
!!       last rebuild (if -1, forcces the rebuild)
!!     - id_nr: identification number, unique for each kpp1 env
!!     - iter: counter for the iteration number
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!s
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!     11.2002 ref-count, actually used to calculate [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_K_P_P1_env_type
     INTEGER :: ref_count, rebuild_each, updates_since_rebuild, id_nr,&
          iter
     TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace
     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: v_ao
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
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: m_epsilon,&
          psi0d, S_psi0
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
!!   SYNOPSIS
!!     Subroutine p_env_init(p_env, qs_env, k_p_p1_env, glob_env, psi0d,&
!!         orthogonal_orbitals, error)
!!       Implicit None
!!       Type(qs_p_env_type), Intent (OUT):: p_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(qs_k_p_p1_env_type), Intent (IN), Target:: k_p_p1_env
!!       Type(global_environment_type), Intent (IN), Target:: glob_env
!!       Type(cp_full_matrix_p_type), Dimension(:), Pointer, Optional::&
!!         psi0d
!!       Logical, Intent (IN), Optional:: orthogonal_orbitals
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine p_env_init
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
!!     - K_P_P1_env: the environment that builds the second order
!!       perturbation kernel.
!!     - glob_env: the global environment
!!     - orthogonal_orbitals: if the orbitals are orthogonal 
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
  SUBROUTINE p_env_init(p_env, qs_env, K_P_P1_env, glob_env,&
                        psi0d, orthogonal_orbitals, error)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(out)                          :: p_env
    TYPE(qs_environment_type), INTENT(in), TARGET             :: qs_env
    TYPE(qs_K_P_P1_env_type), INTENT(in), TARGET              :: K_P_P1_env
    TYPE(global_environment_type), INTENT(in), TARGET         :: glob_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:),POINTER, OPTIONAL :: psi0d
    LOGICAL, INTENT(in), OPTIONAL                             :: orthogonal_orbitals
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)              :: error

    ! locals
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER      :: routineN='p_env_init',&
                                        routineP=moduleN//':'//routineN
    INTEGER                          :: n_ao, n_mo, n_spins, spin, istat
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(dft_control_type), POINTER  :: dft_control
    TYPE(cp_full_matrix_type), POINTER :: qs_env_c
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: ao_mo_fm_pools
    TYPE(cp_para_env_type), POINTER :: para_env

    ! code

    failure=.FALSE.

    NULLIFY(atomic_kind_set, dft_control, ao_mo_fm_pools, para_env)
    NULLIFY(p_env%K_P_P1, &
            p_env%m_epsilon, &
            p_env%psi0d, &
            p_env%S_psi0, &
            p_env%K_P_P1_env, &
            p_env%qs_env)

    p_env%qs_env => qs_env
    p_env%K_P_P1_env => K_P_P1_env

    CALL get_qs_env(qs_env, ao_mo_fm_pools=ao_mo_fm_pools,para_env=para_env)

    n_spins = p_env%qs_env%dft_control%nspins
    DO spin=1, n_spins

       IF (PRESENT(psi0d)) THEN
          CALL cp_fm_get_info(psi0d(spin)%matrix, &
                                     ncol_global=n_mo, nrow_global=n_ao)
       ELSE
       CALL get_mo_set(p_env%qs_env%c(spin)%mo_set, eigenvectors=qs_env_c)
       CALL cp_fm_get_info(qs_env_c, &
            ncol_global=n_mo, nrow_global=n_ao)
       END IF
       p_env%n_mo(spin) = n_mo
       p_env%n_ao(spin) = n_ao
    END DO

    p_env%orthogonal_orbitals=.FALSE.
    IF (PRESENT(orthogonal_orbitals)) &
         p_env%orthogonal_orbitals=orthogonal_orbitals

    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=p_env%S_psi0,&
                                     name="p_env_S_psi0", &
                                     error=error)
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=p_env%m_epsilon,&
                                     name="p_env_m_epsilon", &
                                     error=error)

    IF (PRESENT(psi0d)) THEN
       IF (ASSOCIATED(psi0d)) THEN
          CALL cp_fm_vect_copy(psi0d,p_env%psi0d,error=error)
       END IF
    ELSE IF (.NOT.p_env%orthogonal_orbitals) THEN
       CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=p_env%psi0d,&
                                        name="p_env_psi0d", &
                                        error=error)
    END IF

    ! do still we need this?
    p_env%s_mstruct_changed=.TRUE.

    p_env%initialized=.TRUE.

    IF (PRESENT(psi0d)) CALL p_env_psi0_changed(p_env, error=error)
    
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
       NULLIFY(p_env%qs_env)
       CALL cp_fm_vect_dealloc(p_env%S_psi0, error=error)
       CALL cp_fm_vect_dealloc(p_env%m_epsilon, error=error)
       CALL cp_fm_vect_dealloc(p_env%psi0d, error=error)
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
!!     To be called after the value of psi0 has changed.
!!     Recalculates the quantities S_psi0 and m_epsilon.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the perturbation environment to set
!!     - psi0: the value of psi0, if not given defaults to the qs_env mos
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

  SUBROUTINE p_env_psi0_changed(p_env, psi0, Hrho_psi0d, error)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(inout)                :: p_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(in), TARGET, OPTIONAL :: psi0
    TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(inout), OPTIONAL :: Hrho_psi0d
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)      :: error

    ! locals
    LOGICAL                                :: failure
    INTEGER                                :: n_spins, spin, stat
    CHARACTER(len=*), PARAMETER            :: routineN = 'p_env_psi0_changed', &
                                              routineP = moduleN//':'//routineN
    TYPE(cp_full_matrix_type), POINTER                 :: tmp
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: my_psi0
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: S
    TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: c
    TYPE(cp_para_env_type), POINTER                    :: para_env

    ! code
    failure=.FALSE.

    CPPrecondition(p_env%initialized, cp_failure_level, routineP, error, failure)

    NULLIFY(tmp, ao_mo_fm_pools,c,my_psi0, para_env)

    n_spins = p_env%qs_env%dft_control%nspins

    IF (.NOT. failure) THEN
       CALL get_qs_env(p_env%qs_env,ao_mo_fm_pools=ao_mo_fm_pools,&
            c=c, S=S, para_env=para_env)

       ! def my_psi0
       IF (PRESENT(psi0)) THEN
          CALL cp_fm_vect_copy(psi0,my_psi0,error=error)
       ELSE
          ALLOCATE(my_psi0(n_spins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO spin=1,n_spins
                NULLIFY(my_psi0(spin)%matrix)
                CALL get_mo_set(c(spin)%mo_set,&
                     eigenvectors=my_psi0(spin)%matrix)
             END DO
          END IF
       END IF

    END IF

    IF (.NOT.failure) THEN
       
       ! def psi0d
       IF (p_env%orthogonal_orbitals) THEN
          
          p_env%psi0d => my_psi0

       ELSE

          DO spin=1,n_spins
             ! m_epsilon=choleski_decomposition(my_psi0^T S my_psi0)^-1
             ! could be optimized by compining next two calls
             CALL cp_sm_fm_multiply(sparse_matrix=S(1)%matrix,&
                  v_in=my_psi0(spin)%matrix,&
                  v_out=p_env%S_psi0(spin)%matrix,&
                  ncol=p_env%n_mo(spin), para_env=para_env, &
                  alpha=-1.0_wp,error=error)
             CALL cp_fm_gemm(transa='T',transb='N',n=p_env%n_mo(spin),&
                  m=p_env%n_mo(spin),k=p_env%n_ao(spin),alpha=1.0_wp,&
                  matrix_a=my_psi0(spin)%matrix, &
                  matrix_b=p_env%S_psi0(spin)%matrix,&
                  beta=0.0_wp,matrix_c=p_env%m_epsilon(spin)%matrix,&
                  error=error)
             CALL cp_fm_cholesky_decompose(p_env%m_epsilon(spin)%matrix,&
                  error=error)
             
             ! psi0d=my_psi0 (my_psi0^T S my_psi0)^-1
             ! faster using cp_fm_cholesky_invert ?
             CALL copy_blacs_to_blacs_matrix(my_psi0(spin)%matrix,&
                  p_env%psi0d(spin)%matrix, error=error)
             CALL cp_fm_triangular_multiply(&
                  triangular_matrix=p_env%m_epsilon(spin)%matrix,&
                  matrix_b=p_env%psi0d(spin)%matrix,side='R',&
                  invert_tr=.TRUE., n_rows=p_env%n_ao(spin),&
                  n_cols=p_env%n_mo(spin),error=error)
             CALL cp_fm_triangular_multiply(&
                  triangular_matrix=p_env%m_epsilon(spin)%matrix,&
                  matrix_b=p_env%psi0d(spin)%matrix,side='R',&
                  transpose_tr=.TRUE.,&
                  invert_tr=.TRUE., n_rows=p_env%n_ao(spin),&
                  n_cols=p_env%n_mo(spin),error=error)
          END DO

       END IF
       
       !-----------------------------------------------------------------------|
       ! calculates                                                            |
       ! m_epsilon = - psi0d^T times K times psi0d                             |
       !           = - [K times psi0d]^T times psi0d (because K is symmetric)  |
       !-----------------------------------------------------------------------|
       DO spin=1, n_spins   
          ! S_psi0 = k times psi0d
          CALL cp_sm_fm_multiply(p_env%qs_env%k(spin)%matrix, &
               p_env%psi0d(spin)%matrix, &
               p_env%S_psi0(spin)%matrix, p_env%n_ao(spin), &
               para_env)            
          ! m_epsilon = -1 times S_psi0^T times psi0d
          CALL cp_fm_gemm('T', 'N', &
               p_env%n_mo(spin), p_env%n_mo(spin), p_env%n_ao(spin), &
               -1.0_dbl, p_env%S_psi0(spin)%matrix, p_env%psi0d(spin)%matrix, &
               0.0_dbl, p_env%m_epsilon(spin)%matrix)
       END DO

       !----------------------------------|
       ! calculates S_psi0 = S * my_psi0  |
       !----------------------------------|
       ! calculating this reduces the mat mult without storing a full aoxao
       ! matrix (for P). If nspin>1 you might consider calculating it on the
       ! fly to spare some memory
       CALL get_qs_env(p_env%qs_env, s=S)
       DO spin=1,n_spins
          CALL cp_sm_fm_multiply(S(1)%matrix, &
               my_psi0(spin)%matrix, &
               p_env%S_psi0(spin)%matrix, &
               p_env%n_mo(spin), para_env=para_env)
       END DO

       ! releases my_psi0
       IF (p_env%orthogonal_orbitals) THEN
          NULLIFY(my_psi0)
       ELSE
          CALL cp_fm_vect_dealloc(my_psi0,error=error)
       END IF

    END IF

  END SUBROUTINE p_env_psi0_changed
!***************************************************************************

END MODULE qs_p_types
