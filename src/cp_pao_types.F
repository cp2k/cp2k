!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_types [1.0] *
!!
!!   NAME
!!     cp_pao_types
!!
!!   FUNCTION
!!     definition of most pao related types (and get methods for the 
!!     pao_env)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_pao_types
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use matrix_types, only: real_matrix_type, replicate_matrix, copy_matrix,&
       allocate_matrix, deallocate_matrix, get_matrix_info, get_block_node,&
       add_block_node, put_block_node
  use global_types, only: global_environment_type
  use atoms, only: atom_info, atom_info_type
  use atomic_kinds, only: kind_info,nkind
  use cp_pao_utils, only: cp_1d_i_array
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_types'

! =========== public objects ==========

! types
  public :: cp_pao_env_type, cp_pao_glob_angles_type, &
       cp_pao_local_angles_type, cp_pao_matrix_type

! pao_env procedures
  public :: cp_get, cp_valid

! special functions
  public :: cp_get_n_atom,&
       cp_get_min_bas_size, cp_get_full_bas_size, &
       cp_get_p_index, cp_set_p_index,&
       cp_get_non_p_index, cp_get_global_env

! underlying functions
  public :: cp_pao_env_get,&
       cp_pao_env_valid, cp_pao_env_get_n_atom,&
       cp_pao_env_get_min_bas_size, cp_pao_env_get_full_bas_size,&
       cp_pao_env_get_excl_bas_size, cp_pao_env_get_p_index,&
       cp_pao_env_set_p_index, cp_pao_env_get_non_p_index,&
       cp_pao_env_get_global_env

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure cp_pao_env_get
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure cp_pao_env_valid
  end interface

!! returns the number of atoms
  interface cp_get_n_atom
     module procedure cp_pao_env_get_n_atom
  end interface

!! returns the size of the minimal basis
  interface cp_get_min_bas_size
     module procedure cp_pao_env_get_min_bas_size
  end interface

!! returns the size of the full basis
  interface cp_get_full_bas_size
     module procedure cp_pao_env_get_full_bas_size
  end interface

!! returns the size of the excluded basis
  interface cp_get_excl_bas_size
     module procedure cp_pao_env_get_excl_bas_size
  end interface

!! returns the projection index
  interface cp_get_p_index
     module procedure cp_pao_env_get_p_index
  end interface

!! sets the projection index
  interface cp_set_p_index
     module procedure cp_pao_env_set_p_index
  end interface

!! returns the index of the excluded basis
  interface cp_get_non_p_index
     module procedure cp_pao_env_get_non_p_index
  end interface

!! returns the global environement
  interface cp_get_global_env
     module procedure cp_pao_env_get_global_env
  end interface

!!***
!****************************************************************************

!!****s* cp_pao_types/cp_pao_glob_angles_type [1.0] *
!!
!!   NAME
!!     cp_pao_glob_angles_type
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
!!     initialized: internal flag to (maybe) find errors that come
!!          from not having initialized the structure. is true if 
!!          this type has been initialized
!!     pao_env: a pointer to the pao_env this angles refer to
!!     angles : a one dimensional array with all the angles that define
!!          the actual unitary transformation
!!     ut_calculated: which blocks of unitary_t are valid
!!     unitary_t: the unitary transformation (be careful, only the
!!           blocks flagged in ut_calculated are valid)
!!     min_overlap_m: overlap in the minimal basis (valid only
!!           if m_s_m_valid is true)
!!     min_hamiltonian_m: the hamiltonian in the small basis
!!           valid only if m_h_m_valid is true
!!     m_s_m_valid: true if the overlap matrix is valid
!!     m_h_m_valid: true if the hamiltonian matrix is valid
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  type cp_pao_glob_angles_type
     !private
     logical :: initialized, m_s_m_valid, m_h_m_valid
     type(cp_pao_env_type), pointer :: pao_env
     real(kind=wp), dimension(:), pointer :: angles
     logical, dimension(:), pointer :: ut_calculated
     type(cp_pao_matrix_type), pointer :: unitary_t
     type(cp_pao_matrix_type), pointer :: min_overlap_m
     type(cp_pao_matrix_type), pointer :: min_hamiltonian_m
  end type cp_pao_glob_angles_type
!!***
!****************************************************************************

!!****s* cp_pao_types/cp_pao_local_angles_type [1.0] *
!!
!!   NAME
!!     cp_pao_local_angles_type
!!
!!   FUNCTION
!!     This structure hold the angles that regard one atom.
!!     it should be used to access the angles when you want an "atom view"
!!     (for example to build the unitary transformation).
!!
!!   NOTES
!!      Do not use the global structure, the idea is that parallelization
!!      (distributing various atoms to various processors) should be easier
!!      using this structure.
!!
!!   ATTRIBUTES
!!     initialized: internal flag to (maybe) find errors that come
!!          from not having initialized the structure. is true if 
!!          this type has been initialized
!!     atom_nr : the atom this object focalizes on
!!     glob_arg: a pointer to the global angles
!!      
!!     readonly (virtual) attributes:
!!
!!     n_min_bas: the size of the minimal basis on this atom.
!!     n_full_bas: the size of the full basis on this atom.
!!     pao_env: the pao environement this local angles lives in.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  type cp_pao_local_angles_type
     !private
     logical :: initialized
     integer :: atom_nr, begins_at
     type(cp_pao_glob_angles_type), pointer :: glob_angles
  end type cp_pao_local_angles_type
!!***
!****************************************************************************

!!****s* cp_pao_types/cp_pao_matrix_type [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_type
!!
!!   FUNCTION
!!     This objects hold a matrix in the full or minimal basis, 
!!     like (for example) the unitary transformation
!!
!!   NOTES
!!     I liked the idea of having different types that could be
!!     distingued at compile time, but they are largely equal,
!!     so i did just one (and I find a standin pointer type
!!     a little overkill). Did someone say templates ;)
!!
!!   ATTRIBUTES
!!     initialized: internal flag to (maybe) find errors that come
!!          from not having initialized the structure. is true if 
!!          this type has been initialized
!!     pao_env: a pointer to the pao_env this unitary transformation
!!          refers to (SHARED)
!!     matrix: normally COPIES the matrix into this object, but if at the
!!             same time you also set should_dealloc_matrix to false
!!             (i.e. the matrix is shared) it is SHARED.
!!     min_bas_rows: true if the rows are in the min bas (otherwise 
!!          they are in the full basis)
!!     min_bas_cols: true if the columns are in the min bas (otherwise 
!!          they are in the full bas)
!!     should_dealloc_matrix: if true (default) deallocates the matrix
!!          when cp_dealloc_ref is called
!!
!!     Read only (virtual) attributes:
!!
!!     col_size: a vector with the size of the columns of the different
!!           blocks (might disappear)
!!     row_size: a vector with the size of the rows of the different
!!           blocks (might disappear)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  type cp_pao_matrix_type
     !private
     logical :: initialized, min_bas_rows, min_bas_cols, should_dealloc_matrix
     type(cp_pao_env_type), pointer :: pao_env
     type(real_matrix_type), pointer ::  matrix
  end type cp_pao_matrix_type
!!***
!****************************************************************************

!!****s* cp_pao_types/cp_pao_env_type [1.0] *
!!
!!   NAME
!!     cp_pao_env_type
!!
!!   FUNCTION
!!     global environement where the info about the various parameters
!!     of the pao method are stored
!!
!!   NOTES
!!     add data about acheived precision?
!!
!!   ATTRIBUTES
!!     initialized: control flag, just to be sure we have initialized 
!!                  everything
!!     global_env: the global environement, for the future
!!     wanted_subspace_pos_precision: the wanted precision for the subspace
!!        position
!!     wanted_in_subspace_precision: the wanted precision for the 
!!        optimization in the subspace itself
!!     nr_min_bas: size of the minimum basis on each atom
!!     nr_full_bas: size of the full basis on each atom
!!     angles_att: the actual value for the angles
!!     subspace_att: the actual unitary transformation (and thus subspace)
!!     projection: represent a projection (this togheter with a 
!!         unitary transformation changes from the small basis to the 
!!         big basis). It is used because pao is more efficent if it uses
!!         small angles, but should be take through local_angles.
!!         The projection should probabily be associated with the kind of
!!         atom (thus different projections might point to the same 
!!         index array, not directly with the atom: this (ugly) structure 
!!         might change.
!!     density_matrix_att: the actual density matrix in the full basis
!!
!!     The following objects are also available here, but they are actually
!!     shared. At the moment they are golbal objects, but I do not feel that
!!     they should be global (their place is in my opinion in a 
!!     cp_environement), so to shield myself from their possible movement
!!     I have pointers to them.
!!
!!     full_hamiltonian_m: the hamiltonian matrix in the full basis
!!     full_overlap_m: the overlap matrix in the full basis
!!
!!     Read only (virtual) attributes:
!!     atomic_info: the atom_info (see above about why I give here this 
!!           global variable)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  type cp_pao_env_type
     !private
     logical :: initialized ! =.false.
     logical :: min_d_m_valid, dealloc_projections
     type(global_environment_type), pointer :: global_env
     real(kind=wp) :: wanted_subspace_pos_precision
     real(kind=wp) :: wanted_in_subspace_precision
     integer, dimension(:), pointer :: nr_min_bas
     integer, dimension(:), pointer :: nr_full_bas
     integer, dimension(:), pointer :: angles_begin_at ! private
     integer :: n_atoms
     type(cp_1d_i_array), dimension(:), pointer :: projection
     type(cp_pao_glob_angles_type), pointer :: angles_att
     type(cp_pao_matrix_type), pointer :: full_overlap_m
     type(cp_pao_matrix_type), pointer :: full_hamiltonian_m
     type(cp_pao_matrix_type), pointer :: min_density_m
  end type cp_pao_env_type
!!***
!***************************************************************************

contains

! ================== cp_pao_env_type === get methods ===================

!!****f* cp_pao_env_methods/cp_pao_env_get [1.0] *
!!
!!   NAME
!!     cp_pao_env_get
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_get(pao_env, global_env,&
!!         wanted_subspace_pos_precision, wanted_in_subspace_precision,&
!!         nr_min_bas, nr_full_bas, angles_att, subspace_att, n_atom,&
!!         ntot_full_bas, ntot_min_bas, projection, min_density_m,&
!!         full_overlap_m, full_hamiltonian_m, atomic_info, error)
!!       Type(cp_pao_env_type), Target, Intent (IN):: pao_env
!!       Type(global_environment_type), Optional, Pointer:: global_env
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         wanted_subspace_pos_precision
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         wanted_in_subspace_precision
!!       Integer, Optional, Dimension(:), Pointer:: nr_min_bas
!!       Integer, Optional, Dimension(:), Pointer:: nr_full_bas
!!       Type(cp_pao_glob_angles_type), Optional, Pointer:: angles_att
!!       Type(cp_pao_matrix_type), Optional, Pointer:: subspace_att,&
!!         min_density_m, full_hamiltonian_m, full_overlap_m
!!       Integer, Optional, Intent (OUT):: n_atom, ntot_min_bas,&
!!         ntot_full_bas
!!       Type(cp_1d_i_array), Dimension(:), Pointer, Optional:: projection
!!       Type(atom_info_type), Dimension(:), Pointer, Optional::&
!!         atomic_info
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_get
!!
!!   FUNCTION
!!     extracting values from the the pao environement
!!
!!   NOTES
!!     all the pointers returned are guaranteed to remain valid until
!!     either the size of a basis is changed or this object gets deallocated
!!     (but the actual values stored in the array may change)
!!     There are (or will be ;) funtions to direclty extract one component
!!
!!   INPUTS
!!     pao_env: the environement you get the information from
!!     n_atom: the number of atoms
!!     ntot_full_bas: the total size of the full basis
!!     ntot_full_bas: the total size of the minimal basis
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     The other variable are described in the cp_pao_env_type type
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!*** **********************************************************************
  subroutine cp_pao_env_get(pao_env,global_env, wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, nr_min_bas,nr_full_bas, angles_att, &
       subspace_att, n_atom, ntot_full_bas, ntot_min_bas, projection,&
       min_density_m, full_overlap_m, full_hamiltonian_m, atomic_info,&
       error)
    type(cp_pao_env_type), target, intent(in) :: pao_env
    type(global_environment_type), optional, pointer :: global_env
    real(kind=wp), optional, intent(out) :: wanted_subspace_pos_precision
    real(kind=wp), optional, intent(out) :: wanted_in_subspace_precision
    integer, optional, dimension(:), pointer :: nr_min_bas
    integer, optional, dimension(:), pointer :: nr_full_bas
    type(cp_pao_glob_angles_type), optional, pointer :: angles_att
    type(cp_pao_matrix_type), optional, pointer :: subspace_att, &
         min_density_m, full_hamiltonian_m, full_overlap_m
    integer, optional, intent(out) :: n_atom, ntot_min_bas, ntot_full_bas
    type(cp_1d_i_array), dimension(:), pointer, optional :: projection
    type(atom_info_type), dimension(:), pointer, optional :: atomic_info
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(global_env)) &
            global_env=>pao_env%global_env
       if (present(wanted_subspace_pos_precision)) &
            wanted_subspace_pos_precision=pao_env%wanted_subspace_pos_precision
       if (present(wanted_in_subspace_precision)) &
            wanted_in_subspace_precision=pao_env%wanted_in_subspace_precision
       if (present(nr_min_bas)) &
            nr_min_bas => pao_env%nr_min_bas
       if (present(nr_full_bas)) &
            nr_full_bas => pao_env%nr_full_bas
       if (present(angles_att)) &
            angles_att => pao_env%angles_att
       if (present(n_atom)) &
            n_atom=cp_get_n_atom(pao_env, error=error)
       if (present(projection)) &
            projection => pao_env%projection
       if (present(min_density_m)) &
            min_density_m => pao_env%min_density_m
       if (present(full_overlap_m)) &
            full_overlap_m => pao_env%full_overlap_m
       if (present(full_hamiltonian_m)) &
            full_hamiltonian_m => pao_env%full_hamiltonian_m
       if (present(atomic_info)) &
            atomic_info => atom_info
       if (present(ntot_full_bas)) then
          CPPrecondition(associated(nr_full_bas), cp_warning_level,routineP,error,failure)
          if (associated(nr_full_bas)) then
             ntot_full_bas=sum(nr_full_bas)
          else
             ntot_full_bas=-1 ! use 0?
          end if
       end if
       if (present(ntot_min_bas)) then
          CPPrecondition(associated(nr_min_bas), cp_warning_level,routineP,error,failure)
          if (associated(nr_min_bas)) then
             ntot_min_bas=sum(nr_min_bas)
          else
             ntot_min_bas=-1 ! use 0?
          end if
       end if
    end if failureIf
  end subroutine cp_pao_env_get
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_atomic_info [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_atomic_info
!!
!!   FUNCTION
!!     returns the atom_info (you should use this if you do not want problems
!!     when/if the atom_info will not be global anymore)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_atomic_info(pao_env, error) result(res)
    type(cp_pao_env_type), intent(in) :: pao_env
    type(atom_info_type), dimension(:), pointer ::res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_atomic_info',&
        routineP=moduleN//':'//routineN
    failure=.false.
    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    res => atom_info
  end function cp_pao_env_get_atomic_info
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_p_index [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_p_index
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_p_index(pao_env, atom_nr, i, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr, i
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_get_p_index
!!
!!   FUNCTION
!!     returns the i index of the projection on atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao environement you get the data from
!!     atom_nr: the number of the atom
!!     i: the index of the index you want, in 1:nr_min_size(atom_nr)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_p_index(pao_env, atom_nr, i,error) result(res)
    integer ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr, i
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_p_index',&
        routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug .or. debug_this_module) then
       CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%projection),cp_failure_level,routineP,error,failure)
       CPPrecondition(atom_nr>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(i>0,cp_failure_level,routineP,error,failure)
       if (cp_debug .and. debug_this_module) then
          if (.not.failure) then
             CPPrecondition(atom_nr<=size(pao_env%projection),cp_failure_level,routineP,error,failure)
             call cp_assert(associated(pao_env%projection(atom_nr)%array),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "projection array not associated in "//&
                  CPSourceFileRef,&
                  error,failure)
          end if
          if (.not.failure) then
             call cp_assert(i<=size(pao_env%projection(atom_nr)%array),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "i out of bounds"//&
                  CPSourceFileRef,&
                  error,failure)
          end if
       end if
    end if

    failureIf: if (.not. failure) then
       res=pao_env%projection(atom_nr)%array(i)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_env_get_p_index
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_non_p_index [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_non_p_index
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_non_p_index(pao_env, atom_nr, i, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr, i
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_get_non_p_index
!!
!!   FUNCTION
!!     returns the i index of the excluded basis on atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao environement you get the data from
!!     atom_nr: the number of the atom
!!     i: the index of the index you want, in 1:nr_excl_size(atom_nr)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_non_p_index(pao_env, atom_nr, i,error) result(res)
    integer ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr, i
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_p_index',&
        routineP=moduleN//':'//routineN
    integer :: j,k
    failure=.false.

    if (cp_debug .or. debug_this_module) then
       CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%projection),cp_failure_level,routineP,error,failure)
       CPPrecondition(atom_nr>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(i>0,cp_failure_level,routineP,error,failure)
       if (cp_debug .and. debug_this_module) then
          if (.not.failure) then
             CPPrecondition(atom_nr<=size(pao_env%projection),cp_failure_level,routineP,error,failure)
             call cp_assert(associated(pao_env%projection(atom_nr)%array),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "projection array not associated in "//&
                  CPSourceFileRef,&
                  error,failure)
          end if
          if (.not.failure) then
             call cp_assert(i<=size(pao_env%projection(atom_nr)%array),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "i out of bounds"//&
                  CPSourceFileRef,&
                  error,failure)
          end if
       end if
    end if

    failureIf: if (.not. failure) then
! to do
! this is a dirty hack, I should cache everything (but at the
! beginning I did not think about it)
       res=1
       k=1
       do j=1,i
          do
             if (k <= size(pao_env%projection(atom_nr)%array)) then
                if(res==pao_env%projection(atom_nr)%array(k)) then
                   res=res+1
                   k=k+1
                   cycle
                end if
             end if
             exit
          end do
          res=res+1
       end do
       res=res-1
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_env_get_non_p_index
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_set_p_index [1.0] *
!!
!!   NAME
!!     cp_pao_env_set_p_index
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_set_p_index(pao_env, atom_nr, i, value, error)
!!       Type(cp_pao_env_type), Intent (INOUT):: pao_env
!!       Integer, Intent (IN):: atom_nr, i, value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_set_p_index
!!
!!   FUNCTION
!!     sets the index of the i value of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: environement where you set the value
!!     atom_nr: the index of the atom you want to set
!!     i: the index yo want set
!!     value: the value you want the index to take
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_env_set_p_index(pao_env,atom_nr,i, value,error)
    type(cp_pao_env_type), intent(inout) :: pao_env
    integer, intent(in) :: atom_nr, i, value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_set_p_index',&
        routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug .or. debug_this_module) then
       CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%projection),cp_failure_level,routineP,error,failure)
       CPPrecondition(atom_nr>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(i>0,cp_failure_level,routineP,error,failure)
       if (cp_debug .and. debug_this_module) then
          if (.not.failure) then
             CPPrecondition(atom_nr<=size(pao_env%projection),cp_failure_level,routineP,error,failure)
             call cp_assert(associated(pao_env%projection(atom_nr)%array),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "projection array not associated in "//&
                  CPSourceFileRef,&
                  error,failure)
          end if
          if (.not.failure) then
             call cp_assert(i<=size(pao_env%projection(atom_nr)%array),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "i out of bounds"//&
                  CPSourceFileRef,&
                  error,failure)
          end if
       end if
    end if

    failureIf: if (.not. failure) then
       pao_env%projection(atom_nr)%array(i)=value
    end if failureIf
  end subroutine cp_pao_env_set_p_index
!***************************************************************************


!!****f* cp_pao_env_methods/cp_pao_env_get_n_atom [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_n_atom
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_n_atom(pao_env, error) Result(res)
!!       Integer:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_get_n_atom
!!
!!   FUNCTION
!!     returns the numer of atoms
!!
!!   NOTES
!!     aliased in the generic interface cp_get_n_atom
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_n_atom(pao_env, error) result(res)
    integer ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_n_atom',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(pao_env%nr_full_bas),cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       res=size(pao_env%nr_full_bas)
    else failureIf
       res=-1 ! use 0?
    end if failureIf
  end function cp_pao_env_get_n_atom
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_min_bas_size [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_min_bas_size
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_min_bas_size(pao_env, atom_nr, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_get_min_bas_size
!!
!!   FUNCTION
!!     returns the size of the minimal basis on atom atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the environement you get the information from
!!     atom_nr: the number of the atom you are interested in
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_min_bas_size(pao_env, atom_nr, error) result(res)
    integer ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_min_bas_size',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(pao_env%nr_min_bas), cp_failure_level,routineP,error,failure)
    CPPrecondition(atom_nr>0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPrecondition(size(pao_env%nr_min_bas)>= atom_nr,cp_failure_level,routineP,error,failure)
    end if
    failureIf: if (.not. failure) then
       res=pao_env%nr_min_bas(atom_nr)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_env_get_min_bas_size
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_full_bas_size [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_full_bas_size
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_full_bas_size(pao_env, atom_nr, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_get_full_bas_size
!!
!!   FUNCTION
!!     returns the size of the full basis on atom atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the environement you get the information from
!!     atom_nr: the number of the atom you are interested in
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_full_bas_size(pao_env, atom_nr, error) result(res)
    integer ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_full_bas_size',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(atom_nr>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(pao_env%nr_full_bas), cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPrecondition(size(pao_env%nr_full_bas)>= atom_nr, cp_failure_level,routineP,error,failure)
    end if
    failureIf: if (.not. failure) then
       res=pao_env%nr_full_bas(atom_nr)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_env_get_full_bas_size
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_excl_bas_size [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_excl_bas_size
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_excl_bas_size(pao_env, atom_nr, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_get_excl_bas_size
!!
!!   FUNCTION
!!     returns the size of the excluded basis on atom atom_nr
!!     (full_bas-min_bas)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the environement you get the information from
!!     atom_nr: the number of the atom you are interested in
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_get_excl_bas_size(pao_env, atom_nr, error) result(res)
    integer ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_get_excl_bas_size',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(atom_nr>0,cp_failure_level,routineP,error,failure)
    if (cp_debug.and.debug_this_module) then
       CPPrecondition(associated(pao_env%nr_full_bas), cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env%nr_min_bas), cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(size(pao_env%nr_full_bas)>= atom_nr, cp_failure_level,routineP,error,failure)
       end if
    end if
    failureIf: if (.not. failure) then
       res=pao_env%nr_full_bas(atom_nr)-pao_env%nr_min_bas(atom_nr)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_env_get_excl_bas_size
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_get_global_env [1.0] *
!!
!!   NAME
!!     cp_pao_env_get_global_env
!!
!!   SYNOPSIS
!!     Function cp_pao_env_get_global_env(pao_env, error) Result(res)
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Type(global_environment_type), Pointer:: res
!!     End Function cp_pao_env_get_global_env
!!
!!   FUNCTION
!!     returns the global environement
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao environement you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
function cp_pao_env_get_global_env(pao_env,error) result(res)
  type(cp_pao_env_type), intent(in) :: pao_env
  type(cp_error_type), optional, intent(inout) :: error
  type(global_environment_type), pointer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_pao_env_get_global_env',&
        routineP=moduleN//':'//routineN
  failure=.false.
  CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     res => pao_env%global_env
  else
     nullify(res)
  end if
end function cp_pao_env_get_global_env
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_valid [1.0] *
!!
!!   NAME
!!     cp_pao_env_valid
!!
!!   SYNOPSIS
!!     Function cp_pao_env_valid(pao_env, error) Result(res)
!!       Logical:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_valid
!!
!!   FUNCTION
!!     perform minimal validation (initialized no dangling pointers)
!!     writes errors as warnings
!!
!!   NOTES
!!     can be called often
!!
!!   INPUTS
!!     pao_env: the environement to validate
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_valid(pao_env,error) result(res)
    logical ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_pao_env_valid',&
        routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (.not.(associated(pao_env%nr_min_bas).and.&
            associated(pao_env%nr_full_bas).and.&
            associated(pao_env%projection).and.&
            associated(pao_env%angles_begin_at).and.&
            associated(pao_env%full_overlap_m).and.&
            associated(pao_env%full_hamiltonian_m))) then
          CPAssert(associated(pao_env%nr_min_bas),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%nr_full_bas),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%projection),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%full_hamiltonian_m),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_env%full_overlap_m),cp_warning_level,routineP,error,failure)
!          associated(pao_env%angles_att).and.&
!          associated(pao_env%min_density_m).and.&
       end if
    end if
    res=.not.failure
  end function cp_pao_env_valid
!***************************************************************************

! ==================== cp_pao_glob_angles_type ====================

! ==================== cp_pao_local_angles_type ====================

! ==================== cp_pao_matrix_type ====================

! ======================= cp_pao_projection =====================

end module cp_pao_types

