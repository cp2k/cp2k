!----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations        !
!   Copyright (C) 2001  CP2K developers group                                !
!----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_types [1.0] *
!!
!!   NAME
!!     pao_types
!!
!!   FUNCTION
!!     definition of most pao related types (and get methods for the 
!!     pao_env)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!     02.2002 splitted in different modules (now almost only type defs here)
!!             [fawzi]
!!     05.2002 adapted to the new qs_* structure, updated doc [fawzi]
!!
!!   SOURCE
!*****************************************************************************
MODULE pao_types
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type,&
                                             cp_get
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_sp_create
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pao_bas_projection_types,        ONLY: cp_get,&
                                             pao_bas_projection_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             copy_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             put_block_node,&
                                             real_matrix_type,&
                                             replicate_matrix
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_types'

! =========== public objects ==========

! types
  PUBLIC :: pao_env_type, pao_glob_angles_type

! pao_env procedures
  PUBLIC :: cp_get, cp_valid

! underlying functions
  PUBLIC :: pao_env_get, pao_env_valid

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  INTERFACE cp_get
     MODULE PROCEDURE pao_env_get
  END INTERFACE

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  INTERFACE cp_valid
     MODULE PROCEDURE pao_env_valid
  END INTERFACE

!!***
!****************************************************************************

!!****s* pao_types/pao_glob_angles_type [1.0] *
!!
!!   NAME
!!     pao_glob_angles_type
!!
!!   FUNCTION
!!     this is the global structure that holds all the angles 
!!     (or something similar, depending form the method)
!!     that define the unitary transformation of pao.
!!
!!   NOTES
!!     They are stored in a big one dimensional array because so
!!     they can be directly used with generic optimizers. 
!!     Should be considered almost private, use only the angles
!!     pointer, and that only if you must, i.e. only for the optimizer
!!
!!   ATTRIBUTES
!!     - initialized: internal flag to (maybe) find errors that come
!!       from not having initialized the structure. is true if 
!!       this type has been initialized
!!     - pao_env: a pointer to the pao_env this angles refer to
!!     - angles : a one dimensional array with all the angles that define
!!       the actual unitary transformation
!!     - ut_calculated: which blocks of unitary_t are valid
!!     - unitary_t: the unitary transformation (be careful, only the
!!       blocks flagged in ut_calculated are valid)
!!     - min_overlap_m: overlap in the minimal basis (valid only
!!       if m_s_m_valid is true)
!!     - min_hamiltonian_m: the hamiltonian in the small basis
!!       valid only if m_h_m_valid is true
!!     - m_s_m_valid: true if the overlap matrix is valid
!!     - m_h_m_valid: true if the hamiltonian matrix is valid
!!     - NUi_injection_calculated: which blocks of the NUi_injection
!!       are valid
!!     - NUi_injection: the projection that goes from the minimal
!!       basis to the full basis. N is an atomic orthogonalisation term
!!       U is the unitary transformation, and "i" is the basis injection.
!!       Only the blocks that have true in NUi_injection_calculated are
!!       valid.
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE pao_glob_angles_type
     !private
     LOGICAL :: initialized, m_s_m_valid, m_h_m_valid
     REAL(KIND = dp), DIMENSION(:), POINTER :: angles
     LOGICAL, DIMENSION(:), POINTER :: ut_calculated, NUi_injection_calculated
     TYPE(cp_block_matrix_type), POINTER :: unitary_t, NUi_injection
     TYPE(cp_block_matrix_p_type), DIMENSION(:), POINTER :: min_overlap_m
     TYPE(cp_block_matrix_p_type), DIMENSION(:), POINTER :: min_hamiltonian_m
  END TYPE pao_glob_angles_type
!!***
!****************************************************************************

!!****s* pao_types/pao_env_type [1.0] *
!!
!!   NAME
!!     pao_env_type
!!
!!   FUNCTION
!!     global environement where the info about the various parameters
!!     of the pao method are stored
!!
!!   NOTES
!!     add data about acheived precision?
!!
!!   ATTRIBUTES
!!     - initialized: control flag, just to be sure we have initialized 
!!       everything
!!     - min_d_m_valid: true if the density min_density_m is valid
!!     - wanted_subspace_pos_precision: the wanted precision for the subspace
!!       position
!!     - wanted_in_subspace_precision: the wanted precision for the 
!!       optimization in the subspace itself
!!     - angles_att: the actual value for the angles
!!     - min_density_m: the actual density matrix in the minimal basis
!!     - min_bas_m_struct: the structure of a symmetric matrix in the 
!!       minimal basis
!!     - full_bas_m_struct: the structure of a symmetric of a matrix in the
!!       full basis
!!     - full_min_m_struct: the structure of a non symmetric matrix with
!!       row in the full basis and columns in the minimal basis
!!     - min_full_m_struct: the structure of a non symmetric matrix with
!!       row in the minimal basis and columns in the full basis
!!     - full_diag_struct: the structure of a non symmetric block diagonal
!!       matrix in the full basis
!!     - atomic_ortho_valid: true if the corresponding atomic orthogonalization
!!       matrix is valid
!!     - atomic_ortho: the atomic orthogonalization matrix 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE pao_env_type
     ! private
     LOGICAL :: initialized ! =.false.
     LOGICAL :: min_d_m_valid
     LOGICAL, DIMENSION(:), POINTER :: atomic_ortho_valid
     REAL(KIND = dp) :: wanted_subspace_pos_precision
     REAL(KIND = dp) :: wanted_in_subspace_precision
     INTEGER, DIMENSION(:), POINTER :: angles_begin_at ! private
     INTEGER :: max_min_bas,max_full_bas,max_excl_bas,tot_min_bas,&
          tot_full_bas, tot_excl_bas, max_l_angles
     TYPE(pao_glob_angles_type), POINTER :: angles_att
     TYPE(cp_block_matrix_p_type), DIMENSION(:), POINTER :: min_density_m
     TYPE(cp_b_matrix_struct_type), POINTER :: min_bas_m_struct, &
          full_bas_m_struct, full_min_bas_m_struct, min_full_bas_m_struct,&
          full_bas_diag_m_struct
     TYPE(cp_block_matrix_type), POINTER :: atomic_ortho
  END TYPE pao_env_type
!!***
!***************************************************************************

CONTAINS

! ================== pao_env_type === get methods ===================

!!****f* pao_env_methods/pao_env_get [1.0] *
!!
!!   NAME
!!     pao_env_get
!!
!!   SYNOPSIS
!!     Subroutine pao_env_get(pao_env, wanted_subspace_pos_precision,&
!!         wanted_in_subspace_precision, angles_att, min_density_m,&
!!         min_bas_m_struct, full_bas_m_struct, full_min_bas_m_struct,&
!!         min_full_bas_m_struct, full_bas_diag_m_struct, error)
!!       Type(pao_env_type), Target, Intent (IN):: pao_env
!!       Real(KIND = dp), Optional, Intent (OUT)::&
!!         wanted_subspace_pos_precision
!!       Real(KIND = dp), Optional, Intent (OUT)::&
!!         wanted_in_subspace_precision
!!       Type(pao_glob_angles_type), Optional, Pointer:: angles_att
!!       Type(cp_block_matrix_p_type), Dimension(:), Optional, Pointer::&
!!         min_density_m
!!       Type(cp_b_matrix_struct_type), Pointer, Optional::&
!!         min_bas_m_struct, full_bas_m_struct, full_min_bas_m_struct,&
!!         min_full_bas_m_struct, full_bas_diag_m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_get
!!
!!   FUNCTION
!!     extracting values from the the pao environement, see pao_qs_env_get
!!     in pao_qs_env_methods for an higher leve get that can access more
!!     values.
!!
!!   NOTES
!!     all the pointers returned are guaranteed to remain valid until
!!     either the size of a basis is changed or this object gets deallocated
!!     (but the actual values stored in the array may change)
!!     There are (or will be ;) funtions to directly extract one component
!!     atomic_ortho is not updated on get (see pao_qs_env_methods)
!!
!!   INPUTS
!!     - pao_env: the environement you get the information from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!    The other variable are described in the pao_env_type type
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE pao_env_get(pao_env, wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, angles_att,min_density_m, &
       min_bas_m_struct, full_bas_m_struct, full_min_bas_m_struct, &
       min_full_bas_m_struct, full_bas_diag_m_struct,error)
    TYPE(pao_env_type), TARGET, INTENT(in) :: pao_env
    REAL(KIND = dp), OPTIONAL, INTENT(out) :: wanted_subspace_pos_precision
    REAL(KIND = dp), OPTIONAL, INTENT(out) :: wanted_in_subspace_precision
    TYPE(pao_glob_angles_type), OPTIONAL, POINTER :: angles_att
    TYPE(cp_block_matrix_p_type), DIMENSION(:), OPTIONAL, POINTER ::&
         min_density_m
    TYPE(cp_b_matrix_struct_type), POINTER, OPTIONAL :: min_bas_m_struct, &
         full_bas_m_struct, full_min_bas_m_struct, min_full_bas_m_struct,&
         full_bas_diag_m_struct
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pao_env_get',&
        routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(pao_env%initialized, cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(wanted_subspace_pos_precision)) &
            wanted_subspace_pos_precision=pao_env%wanted_subspace_pos_precision
       IF (PRESENT(wanted_in_subspace_precision)) &
            wanted_in_subspace_precision=pao_env%wanted_in_subspace_precision
       IF (PRESENT(angles_att)) angles_att => pao_env%angles_att
       IF (PRESENT(min_density_m)) min_density_m => pao_env%min_density_m
       IF (PRESENT(min_bas_m_struct)) &
            min_bas_m_struct => pao_env%min_bas_m_struct
       IF (PRESENT(full_bas_m_struct)) &
            full_bas_m_struct => pao_env%full_bas_m_struct
       IF (PRESENT(full_min_bas_m_struct)) &
            full_min_bas_m_struct => pao_env%full_min_bas_m_struct
       IF (PRESENT(min_full_bas_m_struct)) &
            min_full_bas_m_struct => pao_env%min_full_bas_m_struct
       IF (PRESENT(full_bas_diag_m_struct)) &
            full_bas_diag_m_struct => pao_env%full_bas_diag_m_struct
    END IF
  END SUBROUTINE pao_env_get
!***************************************************************************

!!****f* pao_env_methods/pao_env_valid [1.0] *
!!
!!   NAME
!!     pao_env_valid
!!
!!   SYNOPSIS
!!     Function pao_env_valid(pao_env, error) Result(res)
!!       Logical:: res
!!       Type(pao_env_type), Intent (IN):: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_env_valid
!!
!!   FUNCTION
!!     perform minimal validation (initialized no dangling pointers)
!!     writes errors as warnings
!!
!!   NOTES
!!     can be called often
!!
!!   INPUTS
!!     - pao_env: the environement to validate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION pao_env_valid(pao_env,error) RESULT(res)
    LOGICAL ::res
    TYPE(pao_env_type), INTENT(in) :: pao_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pao_env_valid',&
        routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.(ASSOCIATED(pao_env%angles_begin_at).AND.&
            ASSOCIATED(pao_env%full_min_bas_m_struct).AND.&
            ASSOCIATED(pao_env%full_bas_diag_m_struct).AND.&
            ASSOCIATED(pao_env%full_bas_m_struct).AND.&
            ASSOCIATED(pao_env%min_full_bas_m_struct).AND.&
            ASSOCIATED(pao_env%min_bas_m_struct).AND.&
            (.not.ASSOCIATED(pao_env%atomic_ortho).OR.&
            ASSOCIATED(pao_env%atomic_ortho_valid)))) THEN
          CPAssert(ASSOCIATED(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
          CPAssert(ASSOCIATED(pao_env%full_min_bas_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(ASSOCIATED(pao_env%full_bas_diag_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(ASSOCIATED(pao_env%full_bas_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(ASSOCIATED(pao_env%min_full_bas_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(ASSOCIATED(pao_env%min_bas_m_struct),cp_warning_level,routineP,error,failure)
          IF (ASSOCIATED(pao_env%atomic_ortho)) THEN
             CPAssert(ASSOCIATED(pao_env%atomic_ortho_valid),cp_warning_level,routineP,error,failure)
          END IF
       END IF
    END IF
    res=.not.failure
  END FUNCTION pao_env_valid
!***************************************************************************

! ==================== pao_glob_angles_type ====================

! ==================== pao_matrix_type ====================

! ======================= pao_bas_proj =====================

END MODULE pao_types

