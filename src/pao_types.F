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
module pao_types
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level
  use cp_error_handling, only: cp_debug, cp_error_type, cp_assert, &
       cp_error_message
  use kinds, only: wp=>dp
  use sparse_matrix_types, only: real_matrix_type, replicate_matrix, copy_matrix,&
       allocate_matrix, deallocate_matrix, get_matrix_info, get_block_node,&
       add_block_node, put_block_node
  use global_types, only: global_environment_type
  use cp_sparse_matrix, only: cp_sparse_matrix_type,cp_sparse_matrix_p_type,&
       cp_sp_create
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type, cp_get
  use pao_bas_projection_types, only: pao_bas_projection_type, cp_get
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_types'

! =========== public objects ==========

! types
  public :: pao_env_type, pao_glob_angles_type

! pao_env procedures
  public :: cp_get, cp_valid

! underlying functions
  public :: pao_env_get, pao_env_valid

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure pao_env_get
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure pao_env_valid
  end interface

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
  type pao_glob_angles_type
     !private
     logical :: initialized, m_s_m_valid, m_h_m_valid
     real(kind=wp), dimension(:), pointer :: angles
     logical, dimension(:), pointer :: ut_calculated, NUi_injection_calculated
     type(cp_sparse_matrix_type), pointer :: unitary_t, NUi_injection
     type(cp_sparse_matrix_p_type), dimension(:), pointer :: min_overlap_m
     type(cp_sparse_matrix_p_type), dimension(:), pointer :: min_hamiltonian_m
  end type pao_glob_angles_type
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
  type pao_env_type
     ! private
     logical :: initialized ! =.false.
     logical :: min_d_m_valid
     logical, dimension(:), pointer :: atomic_ortho_valid
     real(kind=wp) :: wanted_subspace_pos_precision
     real(kind=wp) :: wanted_in_subspace_precision
     integer, dimension(:), pointer :: angles_begin_at ! private
     integer :: max_min_bas,max_full_bas,max_excl_bas,tot_min_bas,&
          tot_full_bas, tot_excl_bas, max_l_angles
     type(pao_glob_angles_type), pointer :: angles_att
     type(cp_sparse_matrix_p_type), dimension(:), pointer :: min_density_m
     type(cp_b_matrix_struct_type), pointer :: min_bas_m_struct, &
          full_bas_m_struct, full_min_bas_m_struct, min_full_bas_m_struct,&
          full_bas_diag_m_struct
     type(cp_sparse_matrix_type), pointer :: atomic_ortho
  end type pao_env_type
!!***
!***************************************************************************

contains

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
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         wanted_subspace_pos_precision
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         wanted_in_subspace_precision
!!       Type(pao_glob_angles_type), Optional, Pointer:: angles_att
!!       Type(cp_sparse_matrix_p_type), Dimension(:), Optional, Pointer::&
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
  subroutine pao_env_get(pao_env, wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, angles_att,min_density_m, &
       min_bas_m_struct, full_bas_m_struct, full_min_bas_m_struct, &
       min_full_bas_m_struct, full_bas_diag_m_struct,error)
    type(pao_env_type), target, intent(in) :: pao_env
    real(kind=wp), optional, intent(out) :: wanted_subspace_pos_precision
    real(kind=wp), optional, intent(out) :: wanted_in_subspace_precision
    type(pao_glob_angles_type), optional, pointer :: angles_att
    type(cp_sparse_matrix_p_type), dimension(:), optional, pointer ::&
         min_density_m
    type(cp_b_matrix_struct_type), pointer, optional :: min_bas_m_struct, &
         full_bas_m_struct, full_min_bas_m_struct, min_full_bas_m_struct,&
         full_bas_diag_m_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_env_get',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized, cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(wanted_subspace_pos_precision)) &
            wanted_subspace_pos_precision=pao_env%wanted_subspace_pos_precision
       if (present(wanted_in_subspace_precision)) &
            wanted_in_subspace_precision=pao_env%wanted_in_subspace_precision
       if (present(angles_att)) angles_att => pao_env%angles_att
       if (present(min_density_m)) min_density_m => pao_env%min_density_m
       if (present(min_bas_m_struct)) &
            min_bas_m_struct => pao_env%min_bas_m_struct
       if (present(full_bas_m_struct)) &
            full_bas_m_struct => pao_env%full_bas_m_struct
       if (present(full_min_bas_m_struct)) &
            full_min_bas_m_struct => pao_env%full_min_bas_m_struct
       if (present(min_full_bas_m_struct)) &
            min_full_bas_m_struct => pao_env%min_full_bas_m_struct
       if (present(full_bas_diag_m_struct)) &
            full_bas_diag_m_struct => pao_env%full_bas_diag_m_struct
    end if
  end subroutine pao_env_get
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
  function pao_env_valid(pao_env,error) result(res)
    logical ::res
    type(pao_env_type), intent(in) :: pao_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_env_valid',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (.not.(associated(pao_env%angles_begin_at).and.&
            associated(pao_env%full_min_bas_m_struct).and.&
            associated(pao_env%full_bas_diag_m_struct).and.&
            associated(pao_env%full_bas_m_struct).and.&
            associated(pao_env%min_full_bas_m_struct).and.&
            associated(pao_env%min_bas_m_struct).and.&
            (.not.associated(pao_env%atomic_ortho).or.&
            associated(pao_env%atomic_ortho_valid)))) then
          CPAssert(associated(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%full_min_bas_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%full_bas_diag_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%full_bas_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%min_full_bas_m_struct),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%min_bas_m_struct),cp_warning_level,routineP,error,failure)
          if (associated(pao_env%atomic_ortho)) then
             CPAssert(associated(pao_env%atomic_ortho_valid),cp_warning_level,routineP,error,failure)
          end if
       end if
    end if
    res=.not.failure
  end function pao_env_valid
!***************************************************************************

! ==================== pao_glob_angles_type ====================

! ==================== pao_matrix_type ====================

! ======================= pao_bas_proj =====================

end module pao_types

