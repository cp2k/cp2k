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
!!     Split and rename the modules qs_p_types,qs_p_build_kernel and
!!     qs_p_utils to qs_p_env_types, qs_p_env_methods and qs_kpp1_types,
!!     qs_kpp1_methods?
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
! ****************************************************************************
MODULE qs_p_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
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
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_matrix,&
                                             fm_pool_get_mstruct,&
                                             fm_pool_give_back_matrix,&
                                             fm_pools_create_matrix_vect,&
                                             fm_pools_give_back_matrix_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type
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
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qs_p_env_type, qs_p_env_p_type, &
       qs_kpp1_env_type, qs_kpp1_env_p_type
  PUBLIC :: p_env_retain, p_env_release, kpp1_release, kpp1_retain

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_types'

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
     LOGICAL :: s_mstruct_changed, orthogonal_orbitals
     INTEGER :: id_nr, ref_count
     TYPE(cp_block_matrix_p_type), DIMENSION(:), POINTER :: &
          kpp1
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: m_epsilon,&
          psi0d, S_psi0
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
!     TYPE(qs_environment_type), POINTER :: qs_env
     INTEGER, DIMENSION(2) :: n_mo, & ! no of molecular orbitals
                              n_ao    ! no of basis functions
  END TYPE qs_p_env_type

!!***
!****************************************************************************

!!****s* qs_p_types/qs_p_env_p_type [1.0] *
!!
!!   NAME
!!     qs_p_env_p_type
!!
!!   FUNCTION
!!     to have arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     p_env: the pointer to the p_env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_p_env_p_type
     TYPE(qs_p_env_type), POINTER :: p_env
  END TYPE qs_p_env_p_type
!!***
  !****************************************************************************

!!****s* qs_p_build_kernel/qs_kpp1_env_type [1.0] *
!!
!!   NAME
!!     qs_kpp1_env_type
!!
!!   FUNCTION
!!     environement that keeps the informations and temporary
!!     val to build the kpp1 kernel matrix
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
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!     11.2002 ref-count, actually used to calculate [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_kpp1_env_type
     INTEGER :: ref_count, rebuild_each, updates_since_rebuild, id_nr,&
          iter
     TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace
     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: v_ao
  END TYPE qs_kpp1_env_type
!!***
  !****************************************************************************

!!****s* qs_p_build_kernel/qs_kpp1_env_p_type [1.0] *
!!
!!   NAME
!!     qs_kpp1_env_p_type
!!
!!   FUNCTION
!!     just to build array of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     kpp1_env: the pointer to the kpp1_env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_kpp1_env_p_type
     TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
  END TYPE qs_kpp1_env_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* qs_p_types/p_env_retain [1.0] *
!!
!!   NAME
!!     p_env_retain
!!
!!   SYNOPSIS
!!     Subroutine p_env_retain(p_env, error)
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine p_env_retain
!!
!!   FUNCTION
!!     retains the given p_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the p_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE p_env_retain(p_env,error)
  TYPE(qs_p_env_type), POINTER :: p_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='p_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
     p_env%ref_count=p_env%ref_count+1
  END IF
END SUBROUTINE p_env_retain
!***************************************************************************

!!****f* qs_p_types/p_env_release [1.0] *
!!
!!   NAME
!!     p_env_release
!!
!!   SYNOPSIS
!!     Subroutine p_env_release(p_env, error)
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine p_env_release
!!
!!   FUNCTION
!!     relases the given p_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - p_env: the environment to release
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
  SUBROUTINE p_env_release(p_env, error)

    ! arguments
    TYPE(qs_p_env_type), POINTER                 :: p_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL                     :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_release',&
                                   routineP = moduleN//':'//routineN

    failure=.FALSE.

    IF (ASSOCIATED(p_env)) THEN
       CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
       p_env%ref_count=p_env%ref_count-1
       IF (p_env%ref_count<1) THEN
          NULLIFY(p_env%kpp1_env)
          CALL kpp1_release(p_env%kpp1_env, error)
          CALL cp_fm_vect_dealloc(p_env%S_psi0, error=error)
          CALL cp_fm_vect_dealloc(p_env%m_epsilon, error=error)
          CALL cp_fm_vect_dealloc(p_env%psi0d, error=error)
          CALL cp_sp_array_dealloc(p_env%kpp1, error=error)
       END IF
    END IF
    NULLIFY(p_env)
  END SUBROUTINE p_env_release
!***************************************************************************

!!****f* qs_p_types/kpp1_release [1.0] *
!!
!!   NAME
!!     kpp1_release
!!
!!   SYNOPSIS
!!     Subroutine kpp1_release(kpp1_env, error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpp1_release
!!
!!   FUNCTION
!!     releases a kpp1_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to release
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
  SUBROUTINE kpp1_release(kpp1_env,error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: ispin, stat
    failure=.FALSE.

    IF (ASSOCIATED(kpp1_env)) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       kpp1_env%ref_count=kpp1_env%ref_count-1
       IF (kpp1_env%ref_count<1) THEN
          IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
             DO ispin=1,SIZE(kpp1_env%v_rspace)
                CALL coeff_deallocate(kpp1_env%v_rspace(ispin))
             END DO
             DEALLOCATE(kpp1_env%v_rspace,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             DO ispin=1,SIZE(kpp1_env%v_ao)
                IF (ASSOCIATED(kpp1_env%v_ao(ispin)%matrix)) THEN
                   CALL deallocate_matrix(kpp1_env%v_ao(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(kpp1_env%v_ao, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
       END IF
    END IF
    NULLIFY(kpp1_env)
  END SUBROUTINE kpp1_release
!***************************************************************************

!!****f* qs_p_types/kpp1_retain [1.0] *
!!
!!   NAME
!!     kpp1_retain
!!
!!   SYNOPSIS
!!     Subroutine kpp1_retain(kpp1_env, error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpp1_retain
!!
!!   FUNCTION
!!     retains a kpp1_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_retain(kpp1_env,error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       kpp1_env%ref_count=kpp1_env%ref_count+1
    END IF
  END SUBROUTINE kpp1_retain
!***************************************************************************

END MODULE qs_p_types
