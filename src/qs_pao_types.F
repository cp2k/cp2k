!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_pao_types [1.0] *
!!
!!   NAME
!!     qs_pao_types
!!
!!   FUNCTION
!!     A central place to keep the global things related to pao
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
module qs_pao_types
  use qs_log_handling
  use qs_error_handling
  use kinds, only: wp=>dp
  use matrix_types, only: real_matrix_type, replicate_matrix, copy_matrix,&
       allocate_matrix, deallocate_matrix, get_matrix_info, get_block_node,&
       add_block_node, put_block_node
  use global_types, only: global_environment_type
  use atoms, only: atom_info, atom_info_type
  use atomic_kinds, only: kind_info,nkind
  use qs_pao_utils
  use qs_pao_plain_rot
  implicit none

  logical, private, parameter :: debug_this_module=.true.
  logical, private, parameter :: standalone=.true.

  private

! =========== public objects ==========
  public :: qs_pao_l_compose_u ! to do remove from here and plain rot tests

! types (with private attributes)
  public :: qs_pao_env_type, qs_pao_glob_angles_type, qs_pao_local_angles_type, &
       qs_pao_matrix_type

! core procedures
  public :: qs_init, qs_dealloc_ref, qs_set, qs_get, qs_create, qs_dealloc,&
       qs_did_change

! core functions
  public :: qs_validate, qs_valid

! special procedures
  public :: qs_set_angle, qs_set_to_identity, qs_set_unitary_t

! special functions
  public :: qs_get_n_atom, qs_get_pao_env, qs_get_angle, &
       qs_get_min_bas_size, qs_get_full_bas_size, &
       qs_get_p_index, qs_set_p_index, qs_get_matrix,&
       qs_get_non_p_index, qs_get_unitary_t,qs_get_glob_angles

! projection - injection
  public :: qs_pao_l_atom_projection, qs_pao_l_atom_injection, &
       qs_pao_r_atom_projection, qs_pao_r_atom_injection, &
       qs_pao_atom_m_to_f_matrix, qs_pao_atom_f_to_m_matrix, &
       qs_pao_matrix_transf_to_m

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface qs_init
     module procedure qs_init_pao_env, qs_init_pao_glob_angles, &
          qs_init_pao_matrix, qs_init_pao_local_angles
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface qs_dealloc_ref
     module procedure qs_dealloc_pao_env_ref, qs_dealloc_pao_glob_angles_ref, &
          qs_dealloc_pao_matrix_ref, qs_dealloc_pao_local_angles_ref
  end interface

!! Allocates and initialize a new object returning a pointer to it in 
!! the first argument
  interface qs_create
     module procedure qs_create_pao_matrix, qs_create_pao_glob_angles
  end interface

!! Deallocates the memory used by the first argument, and the object
!! itself
  interface qs_dealloc
     module procedure qs_dealloc_pao_matrix, qs_dealloc_pao_glob_angles
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface qs_validate
     module procedure qs_validate_pao_env, qs_validate_pao_glob_angles, &
          qs_validate_pao_matrix, qs_validate_pao_local_angles
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface qs_valid
     module procedure qs_pao_env_valid, qs_pao_matrix_valid, &
          qs_pao_glob_angles_valid, qs_pao_l_angles_valid
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface qs_set
     module procedure qs_set_pao_env, qs_set_pao_glob_angles, &
          qs_set_pao_local_angles, qs_set_pao_matrix
  end interface

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface qs_get
     module procedure qs_pao_env_get, qs_pao_glob_angles_get, &
          qs_pao_matrix_get, qs_pao_local_angles_get
  end interface

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  interface qs_did_change
     module procedure qs_l_angles_did_change, qs_glob_angles_did_change
  end interface

! +++ end of the core methods +++

!! returns the number of atoms
  interface qs_get_n_atom
     module procedure qs_pao_env_get_n_atom
  end interface

!! returns the pao environement the first object lives in
  interface qs_get_pao_env
     module procedure qs_pao_glob_angles_get_pao_env, &
          qs_pao_matrix_get_pao_env, qs_pao_local_angles_get_pao_env
  end interface

!! returns the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) 
  interface qs_get_angle
     module procedure qs_pao_local_angles_get_angle
  end interface

!! sets the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) to the value in the fourth (value)
  interface qs_set_angle
     module procedure qs_pao_local_angles_set_angle
  end interface

!! returns the size of the minimal basis
  interface qs_get_min_bas_size
     module procedure qs_pao_env_get_min_bas_size, qs_l_angles_get_min_bas_size
  end interface

!! returns the size of the full basis
  interface qs_get_full_bas_size
     module procedure qs_pao_env_get_full_bas_size, &
          qs_l_angles_get_full_bas_size
  end interface

!! returns the size of the excluded basis
  interface qs_get_excl_bas_size
     module procedure qs_pao_env_get_excl_bas_size,&
          qs_l_angles_get_excl_bas_size
  end interface

!! returns the projection index
  interface qs_get_p_index
     module procedure qs_pao_env_get_p_index, qs_pao_local_angles_get_p_index
  end interface

!! sets the projection index
  interface qs_set_p_index
     module procedure qs_pao_env_set_p_index, qs_pao_local_angles_set_p_index
  end interface

!! returns the index of the excluded basis
  interface qs_get_non_p_index
     module procedure qs_pao_env_get_non_p_index,&
          qs_pao_l_angles_get_non_p_index
  end interface

!! sets the matrix to a unit matrix (and thus also guarantee that 
!! the diagonal blocks are present)
  interface qs_set_to_identity
     module procedure qs_pao_matrix_set_to_unity
  end interface

!! gets the underlying real matrix
  interface qs_get_matrix
     module procedure qs_pao_matrix_get_matrix
  end interface

!! gets the unitary transformation
  interface qs_get_unitary_t
     module procedure qs_l_angles_get_unitary_t
  end interface

!! sets the unitary transformation
  interface qs_set_unitary_t
     module procedure qs_l_angles_set_unitary_t
  end interface

!! gets the undelying glob_angles object
  interface qs_get_glob_angles
     module procedure qs_l_angles_get_glob_angles
  end interface
!!***
!****************************************************************************

!!****s* qs_pao_types/qs_pao_glob_angles_type [1.0] *
!!
!!   NAME
!!     qs_pao_glob_angles_type
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
  type qs_pao_glob_angles_type
     private
     logical :: initialized, m_s_m_valid, m_h_m_valid
     type(qs_pao_env_type), pointer :: pao_env
     real(kind=wp), dimension(:), pointer :: angles
     logical, dimension(:), pointer :: ut_calculated
     type(qs_pao_matrix_type), pointer :: unitary_t
     type(qs_pao_matrix_type), pointer :: min_overlap_m
     type(qs_pao_matrix_type), pointer :: min_hamiltonian_m
  end type qs_pao_glob_angles_type
!!***
!****************************************************************************

!!****s* qs_pao_types/qs_pao_local_angles_type [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_type
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
  type qs_pao_local_angles_type
     private
     logical :: initialized
     integer :: atom_nr, begins_at
     type(qs_pao_glob_angles_type), pointer :: glob_angles
  end type qs_pao_local_angles_type
!!***
!****************************************************************************

!!****s* qs_pao_types/qs_pao_matrix_type [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_type
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
!!          when qs_dealloc_ref is called
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
  type qs_pao_matrix_type
     private
     logical :: initialized, min_bas_rows, min_bas_cols, should_dealloc_matrix
     type(qs_pao_env_type), pointer :: pao_env
     type(real_matrix_type), pointer ::  matrix
  end type qs_pao_matrix_type
!!***
!****************************************************************************

!!****s* qs_pao_types/qs_pao_env_type [1.0] *
!!
!!   NAME
!!     qs_pao_env_type
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
!!     qs_environement), so to shield myself from their possible movement
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
  type qs_pao_env_type
     private
     logical :: initialized ! =.false.
     logical :: min_d_m_valid, dealloc_projections
     type(global_environment_type), pointer :: global_env
     real(kind=wp) :: wanted_subspace_pos_precision
     real(kind=wp) :: wanted_in_subspace_precision
     integer, dimension(:), pointer :: nr_min_bas
     integer, dimension(:), pointer :: nr_full_bas
     integer, dimension(:), pointer :: angles_begin_at ! private
     integer :: n_atoms
     type(qs_1d_i_array), dimension(:), pointer :: projection
     type(qs_pao_glob_angles_type), pointer :: angles_att
     type(qs_pao_matrix_type), pointer :: full_overlap_m
     type(qs_pao_matrix_type), pointer :: full_hamiltonian_m
     type(qs_pao_matrix_type), pointer :: min_density_m
  end type qs_pao_env_type
!!***
!***************************************************************************

contains

! ==================== qs_pao_env_type ====================
!begin

! ====== core methods =====

!!****f* qs_pao_types/qs_init_pao_env [1.0] *
!!
!!   NAME
!!     qs_init_pao_env
!!
!!   FUNCTION
!!     initialize the pao environement
!!
!!   NOTES
!!     every environement HAS to be initialized with this routine before 
!!     using.
!!     Unless said all arguments are COPIED if nedeed, not shared.
!!     global_env, is SHARED.
!!
!!   INPUTS
!!     pao_env: the environement that has to be initialized
!!     projection: the projection, if given, the projection subarrays 
!!          are shared
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!     the other arguments are actually a shortcut not to call
!!     qs_set_pao_env, for their description see qs_set_pao_env.
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
  subroutine qs_init_pao_env(pao_env,nr_min_bas, nr_full_bas,global_env,&
       full_overlap_sm, full_hamiltonian_sm,&
       wanted_subspace_pos_precision,&
       wanted_in_subspace_precision,&
       projection, error)
    type(qs_pao_env_type), intent(out) :: pao_env
    type(global_environment_type), intent(in), target :: global_env
    integer, dimension(:), intent(in) :: nr_min_bas, nr_full_bas
    type(real_matrix_type), target, intent(in), optional :: &
         full_hamiltonian_sm, full_overlap_sm
    real(kind=wp), optional, intent(in) :: wanted_in_subspace_precision, &
         wanted_subspace_pos_precision
    type(qs_1d_i_array), dimension(:), optional, intent(in) :: projection
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_init_pao_env'
    type(qs_error_type) :: iError
    integer :: stat, i
    failure=.false.
    call qs_init(iError, template_error=error)

    nullify(pao_env%nr_min_bas, pao_env%global_env, pao_env%nr_full_bas)
    nullify(pao_env%projection, pao_env%angles_begin_at)
    pao_env%min_d_m_valid=.false.
    pao_env%n_atoms=-1
! global_env
    pao_env%global_env => global_env
    call qs_assert(associated(pao_env%global_env),qs_warning_level,&
         qs_assertion_failed,routineN,"glob env not associated in "//&
         QSSourceFileRef(),error)
! min_bas
    call qs_guarantee_size(pao_env%nr_min_bas,size(nr_min_bas), &
         error=iError)
    QSPostcondition(qs_get_error_level(iError)<qs_failure_level,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       pao_env%nr_min_bas=nr_min_bas
    end if
! full_bas
    call qs_guarantee_size(pao_env%nr_full_bas,size(nr_full_bas), &
         error=iError)
    QSPostcondition(qs_get_error_level(iError)<qs_failure_level,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       pao_env%nr_full_bas=nr_full_bas
    end if
    call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
         error,failure,qs_failure_level)
! check min_bas/full_bas
    QSPrecondition(associated(pao_env%nr_min_bas),qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(pao_env%nr_full_bas),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_assert(size(pao_env%nr_full_bas)==size(pao_env%nr_min_bas),&
            qs_failure_level,qs_assertion_failed,routineN,&
            "min_bas and full bas must have the same size in "//&
            QSSourceFileRef(),error,failure)
    end if
! angles_begin_at
    if (.not.failure) then
       pao_env%n_atoms=size(pao_env%nr_full_bas)
       allocate(pao_env%angles_begin_at(pao_env%n_atoms+1),stat=stat)
       QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    end if
    if (.not.failure) then
       if (pao_env%n_atoms>0) pao_env%angles_begin_at(1)=0
       do i=1,pao_env%n_atoms
          pao_env%angles_begin_at(i+1)=pao_env%angles_begin_at(i)+&
               (pao_env%nr_full_bas(i)-pao_env%nr_min_bas(i))*&
               pao_env%nr_min_bas(i)
       end do
    end if
! others
    pao_env%wanted_subspace_pos_precision=1.0e-4_wp
    pao_env%wanted_in_subspace_precision=1.0e-4_wp
    pao_env%initialized=.true.
! projection
    if (present(projection)) then
       allocate(pao_env%projection(size(projection)),stat=stat)
       QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
       if (.not.failure) then
          pao_env%projection=projection
       end if
       pao_env%dealloc_projections=.false.
    else
       pao_env%dealloc_projections=.true.
       call qs_pao_env_init_projection(pao_env,error)
    end if
! n_projection
! call qs_pao_env_init_n_projection(pao_env,error)
! full_h
    allocate(pao_env%full_hamiltonian_m,stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_init(pao_env%full_hamiltonian_m, matrix=full_hamiltonian_sm, &
            min_bas_cols=.false., min_bas_rows=.false., pao_env=pao_env,&
            should_dealloc_matrix=.not.present(full_hamiltonian_sm), error=error)
    end if
! full_s
    allocate(pao_env%full_overlap_m,stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_init(pao_env%full_overlap_m, matrix=full_overlap_sm,&
            min_bas_cols=.false., min_bas_rows=.false., pao_env=pao_env,&
            should_dealloc_matrix=.not.present(full_overlap_sm), error=error)
    end if
! min_d
    allocate(pao_env%min_density_m,stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_init_pao_matrix(pao_env%min_density_m, min_bas_rows=.true.,&
            min_bas_cols=.true., pao_env=pao_env,&
            should_dealloc_matrix=.true., error=error)
    end if
! angles_att
    allocate(pao_env%angles_att,stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_init_pao_glob_angles(pao_env%angles_att, pao_env=pao_env,&
            error=error)
    end if
    call qs_set_pao_env(pao_env,&
         wanted_subspace_pos_precision=wanted_subspace_pos_precision,&
         wanted_in_subspace_precision=wanted_in_subspace_precision,&
         error=error)
    if (qs_debug.and.debug_this_module) then
       call qs_assert(qs_validate(pao_env,test_pos_att=.false.,error=error),&
            qs_failure_level,qs_assertion_failed,routineN,&
            "pao_env is not valid after init in "//&
            QSSourceFileRef(),error,failure)
    end if
    call qs_dealloc_ref(iError,error)
  end subroutine qs_init_pao_env
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_dealloc_pao_env_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_env_ref
!!
!!   FUNCTION
!!     deallocate the memory that the pao environement has allocated,
!!     but NOT the pao environement itself.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the environement to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_pao_env_ref(pao_env, error)
    type(qs_pao_env_type), intent(inout) :: pao_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_dealloc_pao_env_ref'
    integer :: stat,i
    failure=.false.

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (associated(pao_env%angles_att)) then
          call qs_dealloc_pao_glob_angles(pao_env%angles_att, error=error)
       end if
       call qs_dealloc_pao_matrix(pao_env%min_density_m,&
            error=error)
       call qs_dealloc_pao_matrix(pao_env%full_hamiltonian_m,&
            error=error)
       call qs_dealloc_pao_matrix(pao_env%full_overlap_m,&
            error=error)
       nullify(pao_env%global_env)
       if (associated(pao_env%projection)) then
          if (pao_env%dealloc_projections) then
             do i=1,size(pao_env%projection)
                if (associated(pao_env%projection(i)%array)) then
                   deallocate(pao_env%projection(i)%array, stat=stat)
                   QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
                end if
             end do
          end if
          deallocate(pao_env%projection,stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
       if (associated(pao_env%nr_min_bas)) then
          deallocate(pao_env%nr_min_bas, stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure) 
       end if
       if (associated(pao_env%nr_full_bas)) then
          deallocate(pao_env%nr_full_bas, stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure) 
       end if
       if (associated(pao_env%angles_begin_at)) then
          deallocate(pao_env%angles_begin_at,stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
    end if failureIf
  end subroutine qs_dealloc_pao_env_ref
!!***
!***************************************************************************


!!****f* qs_pao_types/qs_validate_pao_env [1.0] *
!!
!!   NAME
!!     qs_validate_pao_env
!!
!!   FUNCTION
!!     returns true if the pao environement is valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     pao_env: the environement to be tested
!!     test_pos_att: if it is false does not test the validity of the
!!          actual angles and unitary transformation. defaults to true
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  function qs_validate_pao_env(pao_env, test_pos_att, error) result(res)
    logical ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    logical ,optional ::test_pos_att
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure, test_pos
    character(len=*), parameter :: routineN='qs_pao_types:qs_validate_pao_env'
    type(qs_logger), pointer :: logger
    integer :: i,m_bas_att,f_bas_att,j
    failure=.false.
    logger => qs_get_error_logger(error)
    if (present(test_pos_att)) then
       test_pos=test_pos_att
    else
       test_pos=.true.
    end if
    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       QSAssert(associated(pao_env%global_env),qs_warning_level,routineN,error,failure)
       QSAssert(pao_env%wanted_subspace_pos_precision<1.0e-2_wp,qs_warning_level,routineN,error,failure)
       QSAssert(pao_env%wanted_in_subspace_precision<1.0e-2_wp,qs_warning_level,routineN,error,failure)
       QSAssert(associated(pao_env%nr_min_bas), qs_warning_level,routineN,error,failure)
       QSAssert(associated(pao_env%nr_full_bas), qs_warning_level,routineN,error,failure)
       QSAssert(size(pao_env%nr_min_bas)>0,qs_warning_level,routineN,error,failure)
       call qs_assert(size(pao_env%nr_min_bas)==size(pao_env%nr_full_bas),&
            qs_warning_level,qs_assertion_failed,routineN,&
            "min_bas and full_s must have the same length in"//&
            QSSourceFileRef(),error,failure)
       QSAssert(size(pao_env%nr_min_bas)==pao_env%n_atoms,qs_warning_level,routineN,error,failure)
       if (.not. standalone) then
          QSAssert(size(pao_env%nr_min_bas)==size(atom_info),qs_warning_level,routineN,error,failure)
       end if
       failureIf2: if (.not.failure .and.(qs_debug.or.test_pos)) then 
! I could use all() and any(), but then the info on the atom nr 
! would be lost...
          do i=1,size(pao_env%nr_min_bas)
             m_bas_att=pao_env%nr_min_bas(i)
             f_bas_att=pao_env%nr_full_bas(i)
             call qs_assert (m_bas_att >= 1, qs_warning_level, &
                  qs_assertion_failed, routineN,&
                  'the atom '//qs_to_string(i)// &
                  ' has a minmal basis that is less than 1', error,failure)
             call qs_assert(m_bas_att<= f_bas_att,qs_warning_level, &
                  qs_assertion_failed, routineN,&
                  ' min bas>number of orbitals for atom'//qs_to_string(i),&
                  error,failure)
             if (.not.standalone) then
                call qs_assert( f_bas_att == &
                     atom_info(i)%last_sgf - atom_info(i)%first_sgf +1, &
                     qs_warning_level,qs_assertion_failed, routineN, &
                     'full bas is different from hamiltonian basis on atom '//&
                     qs_to_string(i),error,failure)
             end if
          end do
       end if failureIf2
       QSAssert(associated(pao_env%projection),qs_warning_level,routineN,error,failure)
       if (.not.failure) then
          do i=1,size(pao_env%projection)
             QSAssert(associated(pao_env%projection(i)%array),qs_warning_level,routineN,error,failure)
             if (.not.failure) then
                call qs_assert(size(pao_env%projection(i)%array)==&
                     pao_env%nr_min_bas(i),qs_warning_level,&
                     qs_assertion_failed,routineN,&
                     "projection has invalid size in "//&
                     QSSourceFileRef(),error,failure)
                QSAssert(all(pao_env%projection(i)%array>0),qs_warning_level,routineN,error,failure)
                do j=2,size(pao_env%projection(i)%array)
                   call qs_assert(pao_env%projection(i)%array(j) > &
                        pao_env%projection(i)%array(j-1),qs_warning_level,&
                        qs_assertion_failed, routineN, &
                        "projection indices are not ordered in "//&
                        QSSourceFileRef(),error,failure)
                end do
             end if
          end do
       end if
! angles_begin_at
       QSPrecondition(associated(pao_env%angles_begin_at),qs_warning_level,routineN,error,failure)
       if (.not.failure) then
          QSPrecondition(pao_env%n_atoms+1==size(pao_env%angles_begin_at),qs_warning_level,routineN,error,failure)
       end if
       if (.not.failure) then
          call qs_assert(all((pao_env%angles_begin_at(2:pao_env%n_atoms+1)-&
               pao_env%angles_begin_at(1:pao_env%n_atoms))==&
               pao_env%nr_min_bas*(pao_env%nr_full_bas-pao_env%nr_min_bas)),&
               qs_warning_level,qs_assertion_failed,routineN,&
               "angles_begin_at is not compatible with nr_min_bas nr_full_bas in "//&
               QSSourceFileRef(),error,failure)
       end if
! pos att
       if (test_pos) then
          if (.not.qs_validate_pao_glob_angles(pao_env%angles_att,&
               error=error)) failure=.true.
          if (.not.qs_validate_pao_matrix(pao_env%min_density_m,&
               error=error)) failure=.true. ! also test idempotency ?
       end if
       if (.not.standalone) then
          if (.not.qs_validate_pao_matrix(pao_env%full_overlap_m,&
               error=error)) failure=.true. ! should be positive definite
          if (.not.qs_validate_pao_matrix(pao_env%full_hamiltonian_m,&
               error=error)) failure=.true. ! should be hermitian
       end if
    end if failureIf
    res=.not.failure
  end function qs_validate_pao_env
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_set_pao_env [1.0] *
!!
!!   NAME
!!     qs_set_pao_env
!!
!!   FUNCTION
!!     changes some of the parameters of the environement
!!
!!   NOTES
!!     arguments, unless otherwise stated are always copied
!!
!!   INPUTS
!!     pao_env: the environement that has to be initialized
!!     global_env: the global environement (for the future), optional,
!!                 it is SHARED
!!     n_atoms: the number of atoms (used to initialize some structures)
!!              optional, it is inferred from nr_min_bas or nr_full_bas
!!              they are present
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!     for the description of the other variables see the attributes 
!!     of qs_pao_env_type
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
  subroutine qs_set_pao_env(pao_env, wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, angles_att, &
       min_density_m,full_overlap_m, full_hamiltonian_m,&
       min_d_valid,error)
    type(qs_pao_env_type), target, intent(inout) :: pao_env
    real(kind=wp), optional, intent(in) :: wanted_in_subspace_precision, &
         wanted_subspace_pos_precision
    type(qs_pao_glob_angles_type),optional, intent(in) :: angles_att
    type(qs_pao_matrix_type),optional, intent(in) :: &
         min_density_m, full_hamiltonian_m, full_overlap_m
    logical, optional :: min_d_valid
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_set_pao_env'
    integer :: stat
    type(qs_error_type) :: iError
    failure=.false.
    call qs_init_error(iError,template_error=error)

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (present(wanted_subspace_pos_precision)) &
         pao_env%wanted_subspace_pos_precision=wanted_subspace_pos_precision
    if (present(wanted_in_subspace_precision)) &
         pao_env%wanted_in_subspace_precision=wanted_in_subspace_precision
    if (present(angles_att)) then
       call qs_set(pao_env%angles_att, source=angles_att, pao_env=pao_env,&
            error=error)
    end if
    if (present(min_density_m)) then
       pao_env%min_d_m_valid=.true.
       call qs_set(pao_env%min_density_m, source=min_density_m,&
            pao_env=pao_env, error=error)
    end if
    if (present(full_hamiltonian_m)) then
       call qs_set(pao_env%full_hamiltonian_m, source=full_hamiltonian_m,&
            pao_env=pao_env, error=error)
    end if
    if (present(full_overlap_m)) then
       call qs_set(pao_env%full_overlap_m, source=full_overlap_m,&
            pao_env=pao_env, error=error)
    end if
    if (present(min_d_valid)) then
       pao_env%min_d_m_valid=min_d_valid
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_set_pao_env
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_get [1.0] *
!!
!!   NAME
!!     qs_pao_env_get
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
!!            see module qs_error_handling 
!!
!!     The other variable are described in the qs_pao_env_type type
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
  subroutine qs_pao_env_get(pao_env,global_env, wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, nr_min_bas,nr_full_bas, angles_att, &
       subspace_att, n_atom, ntot_full_bas, ntot_min_bas, projection,&
       min_density_m, full_overlap_m, full_hamiltonian_m, atomic_info,&
       error)
    type(qs_pao_env_type), target, intent(in) :: pao_env
    type(global_environment_type), optional, pointer :: global_env
    real(kind=wp), optional, intent(out) :: wanted_subspace_pos_precision
    real(kind=wp), optional, intent(out) :: wanted_in_subspace_precision
    integer, optional, dimension(:), pointer :: nr_min_bas
    integer, optional, dimension(:), pointer :: nr_full_bas
    type(qs_pao_glob_angles_type), optional, pointer :: angles_att
    type(qs_pao_matrix_type), optional, pointer :: subspace_att, &
         min_density_m, full_hamiltonian_m, full_overlap_m
    integer, optional, intent(out) :: n_atom, ntot_min_bas, ntot_full_bas
    type(qs_1d_i_array), dimension(:), pointer, optional :: projection
    type(atom_info_type), dimension(:), pointer, optional :: atomic_info
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get'
    failure=.false.

    QSPrecondition(pao_env%initialized, qs_failure_level,routineN,error,failure)
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
            n_atom=qs_get_n_atom(pao_env, error=error)
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
          QSPrecondition(associated(nr_full_bas), qs_warning_level,routineN,error,failure)
          if (associated(nr_full_bas)) then
             ntot_full_bas=sum(nr_full_bas)
          else
             ntot_full_bas=-1 ! use 0?
          end if
       end if
       if (present(ntot_min_bas)) then
          QSPrecondition(associated(nr_min_bas), qs_warning_level,routineN,error,failure)
          if (associated(nr_min_bas)) then
             ntot_min_bas=sum(nr_min_bas)
          else
             ntot_min_bas=-1 ! use 0?
          end if
       end if
    end if failureIf
  end subroutine qs_pao_env_get
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_valid [1.0] *
!!
!!   NAME
!!     qs_pao_env_valid
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_env_valid(pao_env,error) result(res)
    logical ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_valid'
    failure=.false.

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (.not.(associated(pao_env%nr_min_bas).and.&
            associated(pao_env%nr_full_bas).and.&
            associated(pao_env%projection).and.&
            associated(pao_env%angles_begin_at).and.&
            associated(pao_env%full_overlap_m).and.&
            associated(pao_env%full_hamiltonian_m))) then
          QSAssert(associated(pao_env%nr_min_bas),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_env%nr_full_bas),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_env%angles_begin_at),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_env%projection),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_env%full_hamiltonian_m),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_env%full_overlap_m),qs_warning_level,routineN,error,failure)
!          associated(pao_env%angles_att).and.&
!          associated(pao_env%min_density_m).and.&
       end if
    end if
    res=.not.failure
  end function qs_pao_env_valid
!!***
!***************************************************************************

! ====== special methods ======

!!****f* qs_pao_types/qs_pao_env_get_atomic_info [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_atomic_info
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_env_get_atomic_info(pao_env, error) result(res)
    type(qs_pao_env_type), intent(in) :: pao_env
    type(atom_info_type), dimension(:), pointer ::res
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_atomic_info'
    failure=.false.
    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    res => atom_info
  end function qs_pao_env_get_atomic_info
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_get_p_index [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_p_index
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_env_get_p_index(pao_env, atom_nr, i,error) result(res)
    integer ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr, i
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_p_index'
    failure=.false.

    if (qs_debug .or. debug_this_module) then
       QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%projection),qs_failure_level,routineN,error,failure)
       QSPrecondition(atom_nr>0,qs_failure_level,routineN,error,failure)
       QSPrecondition(i>0,qs_failure_level,routineN,error,failure)
       if (qs_debug .and. debug_this_module) then
          if (.not.failure) then
             QSPrecondition(atom_nr<=size(pao_env%projection),qs_failure_level,routineN,error,failure)
             call qs_assert(associated(pao_env%projection(atom_nr)%array),&
                  qs_failure_level,qs_assertion_failed,routineN,&
                  "projection array not associated in "//&
                  QSSourceFileRef(),error,failure)
          end if
          if (.not.failure) then
             call qs_assert(i<=size(pao_env%projection(atom_nr)%array),&
                  qs_failure_level,qs_assertion_failed,routineN,&
                  "i out of bounds"//&
                  QSSourceFileRef(),error,failure)
          end if
       end if
    end if

    failureIf: if (.not. failure) then
       res=pao_env%projection(atom_nr)%array(i)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_env_get_p_index
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_get_non_p_index [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_non_p_index
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_env_get_non_p_index(pao_env, atom_nr, i,error) result(res)
    integer ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr, i
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_p_index'
    integer :: j,k
    failure=.false.

    if (qs_debug .or. debug_this_module) then
       QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%projection),qs_failure_level,routineN,error,failure)
       QSPrecondition(atom_nr>0,qs_failure_level,routineN,error,failure)
       QSPrecondition(i>0,qs_failure_level,routineN,error,failure)
       if (qs_debug .and. debug_this_module) then
          if (.not.failure) then
             QSPrecondition(atom_nr<=size(pao_env%projection),qs_failure_level,routineN,error,failure)
             call qs_assert(associated(pao_env%projection(atom_nr)%array),&
                  qs_failure_level,qs_assertion_failed,routineN,&
                  "projection array not associated in "//&
                  QSSourceFileRef(),error,failure)
          end if
          if (.not.failure) then
             call qs_assert(i<=size(pao_env%projection(atom_nr)%array),&
                  qs_failure_level,qs_assertion_failed,routineN,&
                  "i out of bounds"//&
                  QSSourceFileRef(),error,failure)
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
  end function qs_pao_env_get_non_p_index
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_set_p_index [1.0] *
!!
!!   NAME
!!     qs_pao_env_set_p_index
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_env_set_p_index(pao_env,atom_nr,i, value,error)
    type(qs_pao_env_type), intent(inout) :: pao_env
    integer, intent(in) :: atom_nr, i, value
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_set_p_index'
    failure=.false.

    if (qs_debug .or. debug_this_module) then
       QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%projection),qs_failure_level,routineN,error,failure)
       QSPrecondition(atom_nr>0,qs_failure_level,routineN,error,failure)
       QSPrecondition(i>0,qs_failure_level,routineN,error,failure)
       if (qs_debug .and. debug_this_module) then
          if (.not.failure) then
             QSPrecondition(atom_nr<=size(pao_env%projection),qs_failure_level,routineN,error,failure)
             call qs_assert(associated(pao_env%projection(atom_nr)%array),&
                  qs_failure_level,qs_assertion_failed,routineN,&
                  "projection array not associated in "//&
                  QSSourceFileRef(),error,failure)
          end if
          if (.not.failure) then
             call qs_assert(i<=size(pao_env%projection(atom_nr)%array),&
                  qs_failure_level,qs_assertion_failed,routineN,&
                  "i out of bounds"//&
                  QSSourceFileRef(),error,failure)
          end if
       end if
    end if

    failureIf: if (.not. failure) then
       pao_env%projection(atom_nr)%array(i)=value
    end if failureIf
  end subroutine qs_pao_env_set_p_index
!!***
!***************************************************************************


!!****f* qs_pao_types/qs_pao_env_get_n_atom [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_n_atom
!!
!!   FUNCTION
!!     returns the numer of atoms
!!
!!   NOTES
!!     aliased in the generic interface qs_get_n_atom
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_env_get_n_atom(pao_env, error) result(res)
    integer ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_n_atom'
    failure=.false.

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(pao_env%nr_full_bas),qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       res=size(pao_env%nr_full_bas)
    else failureIf
       res=-1 ! use 0?
    end if failureIf
  end function qs_pao_env_get_n_atom
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_get_min_bas_size [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_min_bas_size
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
!!            see module qs_error_handling 
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
  function qs_pao_env_get_min_bas_size(pao_env, atom_nr, error) result(res)
    integer ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_min_bas_size'
    failure=.false.

    QSPrecondition(pao_env%initialized, qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(pao_env%nr_min_bas), qs_failure_level,routineN,error,failure)
    QSPrecondition(atom_nr>0,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       QSPrecondition(size(pao_env%nr_min_bas)>= atom_nr,qs_failure_level,routineN,error,failure)
    end if
    failureIf: if (.not. failure) then
       res=pao_env%nr_min_bas(atom_nr)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_env_get_min_bas_size
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_get_full_bas_size [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_full_bas_size
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
!!            see module qs_error_handling 
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
  function qs_pao_env_get_full_bas_size(pao_env, atom_nr, error) result(res)
    integer ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_full_bas_size'
    failure=.false.

    QSPrecondition(pao_env%initialized, qs_failure_level,routineN,error,failure)
    QSPrecondition(atom_nr>0,qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(pao_env%nr_full_bas), qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       QSPrecondition(size(pao_env%nr_full_bas)>= atom_nr, qs_failure_level,routineN,error,failure)
    end if
    failureIf: if (.not. failure) then
       res=pao_env%nr_full_bas(atom_nr)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_env_get_full_bas_size
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_env_get_excl_bas_size [1.0] *
!!
!!   NAME
!!     qs_pao_env_get_excl_bas_size
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
!!            see module qs_error_handling 
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
  function qs_pao_env_get_excl_bas_size(pao_env, atom_nr, error) result(res)
    integer ::res
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_get_excl_bas_size'
    failure=.false.

    QSPrecondition(pao_env%initialized, qs_failure_level,routineN,error,failure)
    QSPrecondition(atom_nr>0,qs_failure_level,routineN,error,failure)
    if (qs_debug.and.debug_this_module) then
       QSPrecondition(associated(pao_env%nr_full_bas), qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%nr_min_bas), qs_failure_level,routineN,error,failure)
       if (.not.failure) then
          QSPrecondition(size(pao_env%nr_full_bas)>= atom_nr, qs_failure_level,routineN,error,failure)
       end if
    end if
    failureIf: if (.not. failure) then
       res=pao_env%nr_full_bas(atom_nr)-pao_env%nr_min_bas(atom_nr)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_env_get_excl_bas_size
!!***
!***************************************************************************

! ========= private methods ========

!!****f* qs_pao_types/qs_pao_env_init_projection [1.0] *
!!
!!   NAME
!!     qs_pao_env_init_projection
!!
!!   FUNCTION
!!     checks the pao_env projection and if not valid initializes it
!!
!!   NOTES
!!     private method
!!
!!   INPUTS
!!     pao_env: the pao environement where check the projection
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_env_init_projection(pao_env, error)
    type(qs_pao_env_type), intent(inout) :: pao_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_env_init_projection'
    integer :: natom, i,stat,j
    failure=.false.
    natom=qs_get_n_atom(pao_env)

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(natom>0,qs_warning_level,routineN,error,failure)
    QSPrecondition(associated(pao_env%nr_min_bas),qs_warning_level,routineN,error,failure)
    if (.not.failure) then
       QSPrecondition(size(pao_env%nr_min_bas)==natom,qs_warning_level,routineN,error,failure)
    end if
    failureIf: if (.not. failure) then
       if (.not.associated(pao_env%projection)) then
          allocate(pao_env%projection(natom),stat=stat)
          QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
          if (.not.failure) then
             do i=1,natom
                nullify(pao_env%projection(i)%array)
             end do
          end if
       end if
       QSAssert(size(pao_env%projection)==qs_get_n_atom(pao_env),qs_failure_level,routineN,error,failure)
       if (.not.failure) then
          do i=1,natom
             if (associated(pao_env%projection(i)%array)) then
                if (size(pao_env%projection(i)%array) /= pao_env%nr_min_bas(i)) then
                   call qs_error_message(qs_warning_level,&
                        "minbas size changed on atom "// qs_to_string(i)//&
                        QSSourceFileRef(),routineN,error=error)
                   deallocate(pao_env%projection(i)%array ,stat=stat)
                   QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
                end if
             end if
             if (.not.associated(pao_env%projection(i)%array)) then
                allocate(pao_env%projection(i)%array(pao_env%nr_min_bas(i)),&
                     stat=stat)
                QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
                if (.not.failure) then
                   do j=1,pao_env%nr_min_bas(i)
                      pao_env%projection(i)%array(j)=j
                   end do
                end if
             end if
          end do
       end if
    end if failureIf
  end subroutine qs_pao_env_init_projection
!!***
!***************************************************************************

!end
! ==================== qs_pao_glob_angles_type ====================
!begin


! ======= core methods ======

!!****f* qs_pao_types/qs_init_pao_glob_angles [1.0] *
!!
!!   NAME
!!     qs_init_pao_glob_angles
!!
!!   FUNCTION
!!     initializes the global angles
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the global angles to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_init_pao_glob_angles(glob_angles, pao_env, angles, source,&
       ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    type(qs_pao_glob_angles_type), intent(out) :: glob_angles
    type(qs_pao_env_type), intent(in), target :: pao_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(qs_pao_glob_angles_type), intent(in), optional :: source
    logical, dimension(:), intent(in), optional :: ut_calculated
    type(qs_pao_matrix_type), intent(in), optional :: unitary_t,&
         min_hamiltonian_m, min_overlap_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_init_pao_glob_angles'
    failure=.false.

    nullify(glob_angles%pao_env, glob_angles%angles, glob_angles%unitary_t,&
         glob_angles%ut_calculated, glob_angles%min_overlap_m, &
         glob_angles%min_hamiltonian_m)
    glob_angles%m_s_m_valid=.false.
    glob_angles%m_h_m_valid=.false.
    glob_angles%pao_env => pao_env
    glob_angles%initialized = .true.
    QSPrecondition(associated(glob_angles%pao_env),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_set_pao_glob_angles(glob_angles, angles=angles,&
            source=source, ut_calculated=ut_calculated,&
            unitary_t=unitary_t,min_hamiltonian_m=min_hamiltonian_m,&
            min_overlap_m=min_overlap_m, error=error)
    end if
  end subroutine qs_init_pao_glob_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_dealloc_pao_glob_angles_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_glob_angles_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by a glob_angles object.
!!     does NOT deallocate the object itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_glob_angles: the object to release
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_pao_glob_angles_ref(glob_angles, error)
    type(qs_pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_dealloc_pao_glob_angles_ref'
    integer ::stat
    failure=.false.

    QSPrecondition(glob_angles%initialized, qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       nullify(glob_angles%pao_env)
       if (associated(glob_angles%angles)) then
          deallocate(glob_angles%angles, stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
       if (associated(glob_angles%ut_calculated)) then
          deallocate(glob_angles%ut_calculated,stat=stat)
          QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
       end if
       call qs_dealloc_pao_matrix(glob_angles%unitary_t,error=error)
       call qs_dealloc_pao_matrix(glob_angles%min_hamiltonian_m,error=error)
       call qs_dealloc_pao_matrix(glob_angles%min_overlap_m,error=error)
    end if failureIf
  end subroutine qs_dealloc_pao_glob_angles_ref
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_dealloc_pao_glob_angles [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_glob_angles
!!
!!   FUNCTION
!!     deallocated the memory allocated by this object and the oject itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the object to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_pao_glob_angles(glob_angles,error)
    type(qs_pao_glob_angles_type), pointer :: glob_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_dealloc_pao_glob_angles'
    integer :: stat
    failure=.false.

    if (associated(glob_angles)) then
       QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
       if (.not. failure) then
          call qs_dealloc_pao_glob_angles_ref(glob_angles,error=error)
       end if
       deallocate(glob_angles,stat=stat)
       QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
    end if
  end subroutine qs_dealloc_pao_glob_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_create_pao_glob_angles [1.0] *
!!
!!   NAME
!!     qs_create_pao_glob_angles
!!
!!   FUNCTION
!!     creates a new glob_angles object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: will be a pointer to the newly created object
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!     for the others inputs see qs_init_pao_glob_angles
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_create_pao_glob_angles(glob_angles, pao_env, angles, source,&
       ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    type(qs_pao_glob_angles_type), pointer :: glob_angles
    type(qs_pao_env_type), intent(in), target :: pao_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(qs_pao_glob_angles_type), intent(in), optional :: source
    logical, dimension(:), intent(in) :: ut_calculated
    type(qs_pao_matrix_type), intent(in) :: unitary_t, min_hamiltonian_m,&
         min_overlap_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_create_pao_glob_angles'
    integer :: stat
    failure=.false.

    allocate(glob_angles, stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       call qs_init_pao_glob_angles(glob_angles, pao_env, angles, source,&
            ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    end if
  end subroutine qs_create_pao_glob_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_validate_pao_glob_angles [1.0] *
!!
!!   NAME
!!     qs_validate_pao_glob_angles
!!
!!   FUNCTION
!!     returns true if the angles are valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     glob_angles: the gangles to be tested
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  function qs_validate_pao_glob_angles(glob_angles,error) result(res)
    logical ::res
    type(qs_pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_validate_pao_glob_angles'
    failure=.false.
    QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       QSAssert(associated(glob_angles%angles),qs_warning_level,routineN,error,failure)
       QSAssert(associated(glob_angles%pao_env),qs_warning_level,routineN,error,failure)
       if (.not.failure) then
          call qs_assert(qs_validate_pao_env(glob_angles%pao_env,&
               test_pos_att=.false.,error=error),qs_warning_level,&
               qs_assertion_failed,routineN,&
               "validation of referred globAngles failed in "//&
               QSSourceFileRef(),error,failure)
       end if
       if (.not.failure) then
          call qs_assert(glob_angles%pao_env%angles_begin_at &
               (glob_angles%pao_env%n_atoms+1)== &
               size(glob_angles%angles),qs_warning_level,&
               qs_assertion_failed, routineN,&
               "number of angles different from last begin_at in "//&
               QSSourceFileRef(),error,failure)
          if (dot_product(glob_angles%pao_env%nr_min_bas,&
               glob_angles%pao_env%nr_full_bas-glob_angles%pao_env%nr_min_bas)&
               /= size(glob_angles%angles)) then
             call qs_assert(.false., qs_warning_level, qs_assertion_failed,&
                  routineN,'wrong number of angles in qs_pao_glob_angles_type',&
                  error,failure)
          end if
       end if
       if (associated(glob_angles%ut_calculated)) then
          call qs_assert(size(glob_angles%ut_calculated)==&
               glob_angles%pao_env%n_atoms, qs_warning_level,&
               qs_assertion_failed, routineN,&
               "cached invalid ut_calculated in "//&
               QSSourceFileRef(),error,failure)
       end if
       if (associated(glob_angles%unitary_t)) then
          call qs_assert(qs_validate_pao_matrix(glob_angles%unitary_t,&
               error=error), qs_warning_level,&
               qs_assertion_failed, routineN,&
               "cached invalid unitary_t in "//&
               QSSourceFileRef(),error,failure)
          QSAssert(associated(glob_angles%ut_calculated),qs_warning_level,routineN,error,failure)
       end if
       if (associated(glob_angles%min_hamiltonian_m)) then
          call qs_assert(qs_validate_pao_matrix(glob_angles%min_hamiltonian_m,&
               error=error), qs_warning_level,&
               qs_assertion_failed, routineN,&
               "cached invalid min_hamiltonian_m in "//&
               QSSourceFileRef(),error,failure)
       else
          QSAssert(.not.glob_angles%m_h_m_valid,qs_warning_level,routineN,error,failure)
       end if
       if (associated(glob_angles%min_overlap_m)) then
          call qs_assert(qs_validate_pao_matrix(glob_angles%min_overlap_m,&
               error=error), qs_warning_level,&
               qs_assertion_failed, routineN,&
               "cached invalid min_overlap_m in "//&
               QSSourceFileRef(),error,failure)
       else
          QSAssert(.not.glob_angles%m_s_m_valid,qs_warning_level,routineN,error,failure)
       end if
    end if failureIf
    res=.not.failure
  end function qs_validate_pao_glob_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_set_pao_glob_angles [1.0] *
!!
!!   NAME
!!     qs_set_pao_glob_angles
!!
!!   FUNCTION
!!     can change the value of the various components of pao_glob_angles
!!
!!   NOTES
!!     you should almost never need this function
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     pao_env: the pao env that this glob_angles refers to (SHARED)
!!     angles: the angles themselves as a long 1d array (copied)
!!     source: a glob angles that is copied into this (other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     ut_calculated: flags to see which blocka have been calculated
!!     unitary_t: the cached unitary matrix. Be careful: only the 
!!             blocks that have been calculated are valid.
!!             Every processor calculates oly the blocks he needs.
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_set_pao_glob_angles(glob_angles, pao_env, angles,source,&
       ut_calculated, unitary_t,min_overlap_m, min_hamiltonian_m,&
       m_s_m_valid, m_h_m_valid, error)
    type(qs_pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_pao_env_type), optional, intent(in), target :: pao_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(qs_pao_glob_angles_type), optional, intent(in) :: source
    type(qs_pao_matrix_type), optional, intent(in) :: min_hamiltonian_m,&
         min_overlap_m, unitary_t
    logical, dimension(:), optional, intent(in) :: ut_calculated
    logical, intent(in), optional :: m_h_m_valid, m_s_m_valid
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_set_pao_glob_angles'
    integer :: stat
    type(qs_error_type) :: iError
    call qs_init_error(iError,template_error=error)
    failure=.false.

    QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (present(source)) then
          QSPrecondition(source%initialized,qs_failure_level,routineN,error,failure)
          if (.not.failure) then
             glob_angles%pao_env => source%pao_env
             call qs_set(iError, print_level=max(qs_failure_level,&
                  qs_get_error_print_level(iError)),&
                  stop_level=max(qs_failure_level,&
                  qs_get_error_stop_level(iError)))
             if (associated(source%angles)) then
                call qs_guarantee_size(glob_angles%angles,size(source%angles),&
                     error=iError)
                call qs_assert(qs_get_error_level(iError)<qs_failure_level,&
                     qs_get_error_level(iError),qs_assertion_failed,routineN,&
                     "sub routine qs_guarantee_size error in"//&
                     QSSourceFileRef(),error,failure)
                glob_angles%angles=source%angles
             else if (associated(glob_angles%angles)) then
                deallocate(glob_angles%angles, stat=stat)
                QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
             end if
             call qs_pao_g_angles_set_unitary_t(glob_angles,&
                  ut_calculated=source%ut_calculated,&
                  new_value=source%unitary_t,&
                  error=error)
             call qs_pao_g_angles_set_min_s_m(glob_angles,&
                  new_value=source%min_overlap_m,&
                  error=error)
             call qs_pao_g_angles_set_min_h_m(glob_angles,&
                  new_value=source%min_hamiltonian_m,&
                  error=error)
             glob_angles%m_h_m_valid=source%m_h_m_valid
             glob_angles%m_s_m_valid=source%m_s_m_valid
             call qs_dealloc_ref(iError,error)
             call qs_init(iError, template_error=error)
          end if
       end if
       if (present(pao_env)) &
            glob_angles%pao_env => pao_env
       if (present(angles)) then
          call qs_guarantee_size(glob_angles%angles,size(angles),&
               error=iError)
          call qs_assert(qs_get_error_level(iError)<qs_failure_level,&
               qs_failure_level,qs_assertion_failed,routineN,&
               "error in qs_guarantee_size called in "//&
               QSSourceFileRef(),error,failure)
          if (.not. failure) then
             glob_angles%angles = angles
          end if
       end if
       call qs_pao_g_angles_set_unitary_t(glob_angles,&
            ut_calculated=ut_calculated,&
            new_value=unitary_t,&
            error=error)
       if (present(min_overlap_m)) then
          call qs_pao_g_angles_set_min_s_m(glob_angles,&
               new_value=min_overlap_m,&
               error=error)
       end if
       if (present(min_hamiltonian_m)) then
          call qs_pao_g_angles_set_min_h_m(glob_angles,&
               new_value=min_hamiltonian_m,&
               error=error)
       end if
       if (qs_debug.and.debug_this_module) then
          if (associated(glob_angles%pao_env)) then
             call qs_assert(qs_validate_pao_env(glob_angles%pao_env, &
                  test_pos_att=.false.), qs_failure_level, &
                  qs_precondition_failed, routineN, &
                  "if set the pao_env must be valid in"//&
                  QSSourceFileRef(),error,failure)
          end if
       end if
       if (.not.failure) then
          call qs_guarantee_size(glob_angles%angles,&
               dot_product(glob_angles%pao_env%nr_min_bas, &
               glob_angles%pao_env%nr_full_bas-&
               glob_angles%pao_env%nr_min_bas), &
               error=error)
       end if
    end if failureIf
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_set_pao_glob_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_glob_angles_get [1.0] *
!!
!!   NAME
!!     qs_pao_glob_angles_get
!!
!!   FUNCTION
!!     retrive the various attributes of the glob angles object
!!
!!   NOTES
!!     try not to use the angles attribute
!!
!!   INPUTS
!!     glob_angles: the object you get the info from
!!     pao_env: the pao environement this object lives in SHARED
!!     angles: a big 1D array with all the angles that define the unitary
!!            transformation (try not to use this)
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling
!!     min_hamiltonian_m: the hamiltonian in the small basis. 
!!            If not yet calculated calculates it (allocating it if necessary)
!!     min_overlap_m: the overlap matrix in the small basis
!!            If not yet calculated calculates it (allocating it if necessary)
!!     ut_calculated: which block diagonal blocks of the unitary trasformation
!!            are valid (allocated if necessary)
!!     unitary_t: unitary trasformation. Be careful, only the blocks flagged
!!            in ut_calculated are valid (allocated if necessary)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_glob_angles_get(glob_angles, pao_env, angles,&
       m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t, &
       min_overlap_m, min_hamiltonian_m, error)
    type(qs_pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_pao_env_type), optional, pointer :: pao_env
    real(kind=wp), optional, dimension(:), pointer :: angles
    logical, intent(out), optional :: m_h_m_valid, m_s_m_valid
    logical, dimension(:), pointer, optional :: ut_calculated
    type(qs_pao_matrix_type), pointer, optional :: unitary_t, min_hamiltonian_m,&
         min_overlap_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_glob_angles_get'
    type(qs_error_type) :: iError
    integer :: stat
    call qs_init(iError,template_error=error)
    failure=.false.

    QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (present(pao_env)) pao_env => glob_angles%pao_env
       if (present(angles)) angles => glob_angles%angles
       if (present(min_overlap_m)) then
          if (.not.glob_angles%m_s_m_valid) then
             if (.not.associated(glob_angles%min_overlap_m)) then
                call qs_create(glob_angles%min_overlap_m, &
                     min_bas_cols=.true., min_bas_rows=.true.,&
                     pao_env=glob_angles%pao_env)
             end if
             call qs_pao_matrix_transf_to_m(source_m=&
                  glob_angles%pao_env%full_overlap_m, &
                  target_m=glob_angles%min_overlap_m,&
                  glob_angles=glob_angles,error=error)
             glob_angles%m_s_m_valid=.true.
          end if
          QSAssert(associated(glob_angles%min_overlap_m),qs_failure_level,routineN,error,failure)
          if (glob_angles%m_s_m_valid) then
             min_overlap_m => glob_angles%min_overlap_m
          else
             nullify(min_overlap_m)
          end if
       end if
       if (present(min_hamiltonian_m)) then
          if (.not.glob_angles%m_h_m_valid) then
             if (.not.associated(glob_angles%min_hamiltonian_m)) then
                call qs_create(glob_angles%min_hamiltonian_m, &
                     min_bas_cols=.true., min_bas_rows=.true.,&
                     pao_env=glob_angles%pao_env,error=iError)
                call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
                     error,failure)
             end if
             if (.not.failure) then
                call qs_pao_matrix_transf_to_m(source_m=&
                     glob_angles%pao_env%full_hamiltonian_m, &
                     target_m=glob_angles%min_hamiltonian_m,&
                     glob_angles=glob_angles,error=error)
                glob_angles%m_h_m_valid=.true.
             end if
          end if
          QSAssert(associated(glob_angles%min_overlap_m),qs_failure_level,routineN,error,failure)
          if (glob_angles%m_h_m_valid) then
             min_hamiltonian_m => glob_angles%min_hamiltonian_m
          else
             nullify(min_hamiltonian_m)
          end if
       end if
       if (present(m_s_m_valid)) m_s_m_valid = glob_angles%m_s_m_valid
       if (present(m_h_m_valid)) m_h_m_valid=glob_angles%m_h_m_valid
       if (present(unitary_t)) then
          if (.not.associated(glob_angles%unitary_t)) then
             call qs_create(glob_angles%unitary_t, min_bas_rows=.false.,&
                  min_bas_cols=.true., pao_env=glob_angles%pao_env,&
                  error=error)
          end if
          if (.not.associated(glob_angles%ut_calculated)) then
             allocate(glob_angles%ut_calculated(glob_angles%pao_env%n_atoms),&
                  stat=stat)
             QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
             if (stat==0) glob_angles%ut_calculated=.false.
          end if
       end if
       if (present(ut_calculated).and.&
            .not.associated(glob_angles%ut_calculated)) then
          allocate(glob_angles%ut_calculated(glob_angles%pao_env%n_atoms),&
               stat=stat)
          QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
          if (stat==0) glob_angles%ut_calculated=.false.
       end if
    end if failureIf
  end subroutine qs_pao_glob_angles_get
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_glob_angles_valid [1.0] *
!!
!!   NAME
!!     qs_pao_glob_angles_valid
!!
!!   FUNCTION
!!     minimal validation on the glob angles object (no dangling pointers)
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     glob_angles: the globa angles object to test
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_glob_angles_valid(glob_angles, error) result(res)
    logical ::res
    type(qs_pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_glob_angles_valid'
    failure=.false.

    QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (.not.(associated(glob_angles%pao_env).and.&
            associated(glob_angles%angles).and.&
            (.not.associated(glob_angles%unitary_t).or.&
            associated(glob_angles%ut_calculated)) .and.&
            (.not.glob_angles%m_s_m_valid .or.&
            associated(glob_angles%min_overlap_m)).and.&
            (.not.glob_angles%m_h_m_valid .or.&
            associated(glob_angles%min_hamiltonian_m))&
            )) then
          QSAssert(associated(glob_angles%pao_env),qs_warning_level,routineN,error,failure)
          QSAssert(associated(glob_angles%angles),qs_warning_level,routineN,error,failure)
          if (associated(glob_angles%unitary_t)) then
             QSAssert(associated(glob_angles%ut_calculated),qs_warning_level,routineN,error,failure)
          end if
          if (glob_angles%m_s_m_valid) then
             QSAssert(associated(glob_angles%min_overlap_m),qs_warning_level,routineN,error,failure)
          end if
          if (glob_angles%m_h_m_valid) then
             QSAssert(associated(glob_angles%min_hamiltonian_m),qs_warning_level,routineN,error,failure)
          end if
       else if(qs_debug.and.debug_this_module) then
          if (.not.(associated(glob_angles%pao_env%nr_full_bas).and.&
               associated(glob_angles%pao_env%nr_min_bas))) then
             failure=qs_valid(glob_angles%pao_env, error=error)
          end if
       end if
    end if
    res=.not.failure
  end function qs_pao_glob_angles_valid
!!***
!***************************************************************************

! ===== special methods ====

!!****f* qs_pao_types/qs_pao_g_angles_set_unitary_t [1.0] *
!!
!!   NAME
!!     qs_pao_g_angles_set_unitary_t
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     new_value: the new cached unitary matrix.
!!     ut_calculated: which blocks are valid
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_g_angles_set_unitary_t(glob_angles, new_value, &
       ut_calculated,error)
    type(qs_pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_pao_matrix_type), target, optional, intent(in) :: new_value
    logical, dimension(:), intent(in), optional, target :: ut_calculated
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_g_angles_set_unitary_t'
    logical, dimension(:),pointer :: ut_c
    type(qs_pao_matrix_type), pointer :: n_value
    type (qs_error_type) :: iError
    call qs_init(iError, template_error=error, &
         print_level=max(qs_failure_level, qs_get_error_print_level(error)))
    failure=.false.

    if (present(ut_calculated)) then
       ut_c => ut_calculated
       if (associated(ut_c)) then
          call qs_guarantee_size(glob_angles%ut_calculated,&
               size(ut_calculated),error=iError)
          call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
               error,failure,qs_failure_level)
          glob_angles%ut_calculated=ut_calculated
       else if (associated(glob_angles%ut_calculated)) then
          deallocate(glob_angles%ut_calculated)
       end if
    end if
    if (present(new_value)) then
       n_value => new_value
       if (associated(n_value)) then
          if (.not.associated(glob_angles%unitary_t)) then
             call qs_create(glob_angles%unitary_t, &
                  pao_env=glob_angles%pao_env,&
                  source=n_value,&
                  error=error)
          else
             call qs_set(glob_angles%unitary_t,&
                  source=n_value,&
                  error=error)
          end if
       else if (associated(glob_angles%unitary_t)) then
! do not dealloc ??
          call qs_dealloc(glob_angles%unitary_t,error=error)
       end if
    end if
  end subroutine qs_pao_g_angles_set_unitary_t
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_g_angles_set_min_s_m [1.0] *
!!
!!   NAME
!!     qs_pao_g_angles_set_min_s_m
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     min_overlap_m: the overlap matrix in the minimal basis
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_g_angles_set_min_s_m(glob_angles, new_value,error)
    type(qs_pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_pao_matrix_type), target :: new_value
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_g_angles_set_min_s_m'
    type(qs_pao_matrix_type), pointer :: n_value
    failure=.false.
    n_value => new_value

    if (associated(n_value)) then
       glob_angles%m_s_m_valid=.true.
       if (.not.associated(glob_angles%min_overlap_m)) then
          call qs_create(glob_angles%min_overlap_m, &
               pao_env=glob_angles%pao_env,&
               source=n_value,&
               error=error)
       else
          call qs_set(glob_angles%min_overlap_m,&
               source=n_value,&
               error=error)
       end if
    else if (associated(glob_angles%min_overlap_m)) then
! do not dealloc ??
       call qs_dealloc(glob_angles%min_overlap_m,error=error)
    end if
  end subroutine qs_pao_g_angles_set_min_s_m
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_g_angles_set_min_h_m [1.0] *
!!
!!   NAME
!!     qs_pao_g_angles_set_min_h_m
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     min_hamiltonian_m: the hamiltonian matrix in the minimal basis
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_g_angles_set_min_h_m(glob_angles, new_value,error)
    type(qs_pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_pao_matrix_type), target :: new_value
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_g_angles_set_min_h_m'
    type(qs_pao_matrix_type), pointer :: n_value
    failure=.false.
    n_value => new_value

    if (associated(n_value)) then
       glob_angles%m_h_m_valid=.true.
       if (.not.associated(glob_angles%min_hamiltonian_m)) then
          call qs_create(glob_angles%min_hamiltonian_m, &
               pao_env=glob_angles%pao_env,&
               source=n_value,&
               error=error)
       else
          call qs_set(glob_angles%min_hamiltonian_m,&
               source=n_value,&
               error=error)
       end if
    else if (associated(glob_angles%min_hamiltonian_m)) then
! do not dealloc ??
       call qs_dealloc(glob_angles%min_hamiltonian_m,error=error)
    end if
  end subroutine qs_pao_g_angles_set_min_h_m
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_glob_angles_get_pao_env [1.0] *
!!
!!   NAME
!!     qs_pao_glob_angles_get_pao_env
!!
!!   FUNCTION
!!     returns the global environement the global angles object lives in
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the global angles object you are interested in
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_glob_angles_get_pao_env(glob_angles, error) result(res)
    type(qs_pao_env_type), pointer ::res
    type(qs_pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_glob_angles_get_pao_env'
    failure=.false.
    QSPrecondition(glob_angles%initialized, qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       res => glob_angles%pao_env
    else failureIf
       nullify(res)
    end if failureIf
  end function qs_pao_glob_angles_get_pao_env
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_glob_angles_did_change [1.0] *
!!
!!   NAME
!!     qs_glob_angles_did_change
!!
!!   FUNCTION
!!     Invalidate all the cached values.
!!     Should be called when the value of the angles change.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the globa angles object that should reset his cache
!!     full_reset: if true deallocates all the memory of the cache
!!            nedeed when the matrix structure changes (i.e. when the 
!!            atoms move). defaults to false
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_glob_angles_did_change(glob_angles, full_reset, error)
    type(qs_pao_glob_angles_type), intent(inout) :: glob_angles
    logical, intent(in), optional :: full_reset
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_glob_angles_did_change'
    integer :: stat
    failure=.false.

    QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       glob_angles%m_h_m_valid=.false.
       glob_angles%m_s_m_valid=.false.
       if (associated(glob_angles%ut_calculated)) then
          glob_angles%ut_calculated=.false.
       end if
       if (present(full_reset)) then
          if (full_reset) then
             call qs_dealloc(glob_angles%unitary_t,error=error)
             call qs_dealloc(glob_angles%min_hamiltonian_m,error=error)
             call qs_dealloc(glob_angles%min_overlap_m)
             if (associated(glob_angles%ut_calculated)) then
! do not dealloc??
                deallocate(glob_angles%ut_calculated, stat=stat)
                QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
             end if
          end if
       end if
    end if
  end subroutine qs_glob_angles_did_change
!!***
!***************************************************************************

!end
! ==================== qs_pao_local_angles_type ====================
!begin

! ====== core methods ======

!!****f* qs_pao_types/qs_init_pao_local_angles [1.0] *
!!
!!   NAME
!!     qs_init_pao_local_angles
!!
!!   FUNCTION
!!     Initializes the pao angles structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the pao_local angles to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling
!!
!!     For the other arguments see qs_set_pao_local_angles
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_init_pao_local_angles(local_angles, atom_nr, glob_angles, &
       source, error)
    type(qs_pao_local_angles_type), intent(out):: local_angles
    integer, optional, intent(in) :: atom_nr
    type(qs_pao_glob_angles_type), optional, target :: glob_angles
    type(qs_pao_local_angles_type), intent(in), optional :: source
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_init_pao_local_angles'
    failure=.false.

    local_angles%initialized=.true.
    local_angles%atom_nr=-1
    nullify(local_angles%glob_angles)
    call qs_assert((present(glob_angles).and.present(atom_nr)).or.&
         present(source), qs_warning_level, qs_precondition_failed,&
         routineN,"You did not fully initialize the local_angles",&
         error,failure)
    call qs_set_pao_local_angles(local_angles, atom_nr, glob_angles, &
         source, error=error)
  end subroutine qs_init_pao_local_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_dealloc_pao_local_angles_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_local_angles_ref
!!
!!   FUNCTION
!!     deallocates the memory that the given local angles structure
!!     has allocated. Does not deallocate the local angles structure itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles structure to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_pao_local_angles_ref(local_angles, error)
    real(kind=wp) ::res
    type(qs_pao_local_angles_type), intent(inout) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_dealloc_pao_local_angles_ref'
    failure=.false.

    QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       nullify(local_angles%glob_angles)
    end if failureIf
  end subroutine qs_dealloc_pao_local_angles_ref
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_validate_pao_local_angles [1.0] *
!!
!!   NAME
!!     qs_validate_pao_local_angles
!!
!!   FUNCTION
!!     checks if the local angles object is valid (returns .true. if it is)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object to check
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_validate_pao_local_angles(local_angles, error) result(res)
    logical ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_validate_pao_local_angles'
    type(qs_pao_env_type), pointer :: pao_env
    integer :: i, beginAt
    failure=.false.
    nullify(pao_env)

    QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
    QSAssert(associated(local_angles%glob_angles),qs_warning_level,routineN,error,failure)
    if (.not. failure) then
       call qs_assert(qs_validate_pao_glob_angles(local_angles%glob_angles,&
            error=error),qs_warning_level,qs_assertion_failed,routineN,&
            "referenced pao_glob_angles is invalid, in file "//&
            QSSourceFileRef(),error,failure)
       pao_env => local_angles%glob_angles%pao_env
! pao_env has been checked by glob angles validation
    end if
    if (.not.failure) then
       QSAssert(local_angles%atom_nr>0,qs_warning_level,routineN,error,failure)
       QSAssert(local_angles%atom_nr<=qs_get_n_atom(pao_env),qs_warning_level,routineN,error,failure)
       if (.not.failure) then
          QSAssert(local_angles%begins_at>=0,qs_warning_level,routineN,error,failure)
          call qs_assert(local_angles%begins_at < &
               size(local_angles%glob_angles%angles),&
               qs_warning_level,qs_assertion_failed,routineN,&
               "begins_at out of range in "//QSSourceFileRef(),&
               error,failure)
          beginAt=sum((pao_env%nr_full_bas(1:(local_angles%atom_nr-1))-&
               pao_env%nr_min_bas(1:(local_angles%atom_nr-1)))*&
               pao_env%nr_min_bas(1:(local_angles%atom_nr-1)))
          QSAssert(local_angles%begins_at==beginAt,qs_warning_level,routineN,error,failure)
       end if
    end if
    res=.not.failure
  end function qs_validate_pao_local_angles
!!***
!***************************************************************************


!!****f* qs_pao_types/qs_set_pao_local_angles [1.0] *
!!
!!   NAME
!!     qs_set_pao_local_angles
!!
!!   FUNCTION
!!     sets the values of the various attributes of the local angles
!!     structure.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the object to change
!!     atom_nr: the numbre of the atom this object focalizes on
!!     glob_angles: the global angles object that this object references
!!          SHARED.
!!     source: a local angles object that is copied into this(other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_set_pao_local_angles(local_angles, atom_nr, glob_angles,&
       source,error)
    type(qs_pao_local_angles_type), intent(inout):: local_angles
    integer, optional, intent(in) :: atom_nr
    type(qs_pao_glob_angles_type), optional, target :: glob_angles 
    type(qs_pao_local_angles_type), intent(in), optional :: source
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_set_pao_local_angles'
    type(qs_pao_env_type), pointer :: pao_env
    type(qs_error_type) :: iError
    failure=.false.
    nullify(pao_env)
    call qs_init_error(iError, template_error=error)

    QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (present(source)) then
          QSPrecondition(source%initialized,qs_failure_level,routineN,error,failure)
          if (.not.failure) then
             local_angles%atom_nr=source%atom_nr
             local_angles%begins_at=source%begins_at
             local_angles%glob_angles => source%glob_angles
          end if
       end if
       if (present(glob_angles)) then
          local_angles%glob_angles => glob_angles
       end if
       if (present(atom_nr)) then
          local_angles%atom_nr=atom_nr
          debugIf: if (qs_debug) then
             QSAssert(atom_nr>0,qs_warning_level,routineN,error,failure)
             if (associated(local_angles%glob_angles)) then
                if (associated(local_angles%glob_angles%pao_env)) then
                   call qs_assert(atom_nr <= &
                        qs_get_n_atom(local_angles%glob_angles%pao_env),&
                        qs_warning_level,qs_assertion_failed,routineN,&
                        "atom_nr too big in "//QSSourceFileRef(),&
                        error,failure)
                end if
             end if
          end if debugIf
       end if
! updates begins_at
       if (associated(local_angles%glob_angles) .and. &
            local_angles%atom_nr>0) then
          pao_env => local_angles%glob_angles%pao_env
          if (associated(pao_env)) then
             if (qs_debug) then
                call qs_assert(qs_validate_pao_env(pao_env, &
                     test_pos_att=.false., error=iError),&
                     qs_warning_level,qs_assertion_failed,routineN,&
                     "the environement should be valid when set in "//&
                     QSSourceFileRef(),error,failure)
                call qs_assert(local_angles%atom_nr <= qs_get_n_atom(pao_env),&
                     qs_warning_level,qs_assertion_failed,routineN,&
                     "the atom_nr should be either negative or valid in "//&
                     QSSourceFileRef(),error,failure)
             end if
             if (.not.failure) then
                local_angles%begins_at= &
                     sum(pao_env%nr_min_bas(1:(local_angles%atom_nr-1))*&
                     (pao_env%nr_full_bas(1:(local_angles%atom_nr-1))-&
                     pao_env%nr_min_bas(1:(local_angles%atom_nr-1))))
             end if
          end if
       end if
    end if failureIf
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_set_pao_local_angles
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_local_angles_get [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_get
!!
!!   FUNCTION
!!     gets various attributes of the pao local angles structure
!!
!!   NOTES
!!     I try to keep the glob_angles private.
!!     I have not added al local angles matrix
!!
!!   INPUTS
!!     local_angles: the local angles to get the info from.
!!     atom_nr: the atom this local angles is focalized on.
!!     n_min_bas: the size of the minimal basis on this atom.
!!     n_full_bas: the size of the full basis on this atom.
!!     n_excl_bas: the size of the excluded basis on this atom
!!            (n_full_bas-n_min_bas).
!!     pao_env: the pao environement this local angles lives in.
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_local_angles_get(local_angles, atom_nr, n_min_bas, &
       n_full_bas, pao_env, n_excl_bas,unitary_t,cache_valid, error)
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    integer, optional, intent(out) :: atom_nr, n_min_bas, n_full_bas,&
         n_excl_bas
    type(qs_pao_env_type), optional, pointer :: pao_env
    type(qs_error_type), optional, intent(inout) :: error
    real(kind=wp), dimension(:,:), pointer, optional :: unitary_t
    logical, intent(out), optional :: cache_valid

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_local_angles_get'
    type(qs_pao_env_type), pointer :: my_pao_env
    failure=.false.
    nullify(my_pao_env)

    QSAssert(local_angles%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (present(atom_nr)) then
          atom_nr=local_angles%atom_nr
       end if
       if (associated(local_angles%glob_angles)) then
          my_pao_env => local_angles%glob_angles%pao_env
       end if
       if (present(pao_env)) then
          pao_env => my_pao_env
       end if
       if (present(n_min_bas)) then
          QSAssert(local_angles%atom_nr>0,qs_failure_level,routineN,error,failure)
! also check upper bound?
          n_min_bas=my_pao_env%nr_min_bas(local_angles%atom_nr)
       end if
       if (present(n_full_bas)) then
          QSAssert(local_angles%atom_nr>0,qs_failure_level,routineN,error,failure)
          n_full_bas=my_pao_env%nr_full_bas(local_angles%atom_nr)
       end if
       if (present(n_excl_bas)) then
          QSAssert(local_angles%atom_nr>0,qs_failure_level,routineN,error,failure)
          n_excl_bas=my_pao_env%nr_full_bas(local_angles%atom_nr)-&
               my_pao_env%nr_min_bas(local_angles%atom_nr)
       end if
       if (present(unitary_t)) then
          unitary_t => qs_l_angles_get_unitary_t(local_angles,error=error)
       end if
       if (present(cache_valid)) then
          if (associated(local_angles%glob_angles%ut_calculated).and.&
               local_angles%atom_nr>0) then
             cache_valid=local_angles%glob_angles%ut_calculated(&
                  local_angles%atom_nr)
          else
             cache_valid=.false.
          end if
       end if
       nullify(my_pao_env)
    end if failureIf
  end subroutine qs_pao_local_angles_get
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_l_angles_valid [1.0] *
!!
!!   NAME
!!     qs_pao_l_angles_valid
!!
!!   FUNCTION
!!     performs minimal validation on the object
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     local_angles: the object to validate
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_l_angles_valid(local_angles, error) result(res)
    logical ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_l_angles_valid'
    failure=.false.

    QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
    QSAssert(associated(local_angles%glob_angles),qs_warning_level,routineN,error,failure)
    if (qs_debug.and.debug_this_module) then
       if (.not. failure) then
          if (.not.associated(local_angles%glob_angles%pao_env)) then
             QSAssert(.false.,qs_warning_level,routineN,error,failure)
          else if (.not.(associated(local_angles%glob_angles%pao_env%nr_full_bas)&
               .and.associated(local_angles%glob_angles%pao_env%nr_min_bas))) then
             failure=qs_valid(local_angles%glob_angles%pao_env, error=error)
          end if
       end if
    end if
    res=.not.failure
  end function qs_pao_l_angles_valid
!!***
!***************************************************************************

! ====== special methods =====

!!****f* qs_pao_types/qs_pao_local_angles_get_pao_env [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_get_pao_env
!!
!!   FUNCTION
!!     returns the global environement the local angles object lives in
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you are interested in
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_local_angles_get_pao_env(local_angles, error) result(res)
    type(qs_pao_env_type), pointer ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_local_angles_get_pao_env'
    failure=.false.
    QSPrecondition(local_angles%initialized, qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(local_angles%glob_angles),qs_warning_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       res => local_angles%glob_angles%pao_env
    else failureIf
       nullify(res)
    end if failureIf
  end function qs_pao_local_angles_get_pao_env
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_local_angles_get_p_index [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_get_p_index
!!
!!   FUNCTION
!!     returns the i-te index of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     i: number of the index you want to get
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_local_angles_get_p_index(local_angles, i,error) result(res)
    integer ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: i
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_local_angles_get_p_index'
    type(qs_pao_env_type), pointer :: pao_env
    failure=.false.
    pao_env => qs_get_pao_env(local_angles)
    if (qs_debug .or. debug_this_module) then
       QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    end if

    failureIf: if (.not. failure) then
       res=qs_pao_env_get_p_index(pao_env,atom_nr=local_angles%atom_nr,&
            i=i,error=error)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_local_angles_get_p_index
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_local_angles_set_p_index [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_set_p_index
!!
!!   FUNCTION
!!     set the value of the i-te index of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you set values in
!!     i: number of the index you want to set
!!     value: the new nalue of the i-te index
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_local_angles_set_p_index(local_angles,i,value,error)
    type(qs_pao_local_angles_type), intent(inout) :: local_angles ! actually only pao_env is changed
    integer, intent(in) :: i,value
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_local_angles_set_p_index'
    type(qs_pao_env_type), pointer :: pao_env
    failure=.false.
    pao_env => qs_get_pao_env(local_angles)
    if (qs_debug .or. debug_this_module) then
       QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    end if

    failureIf: if (.not. failure) then
       call qs_pao_env_set_p_index(pao_env,atom_nr=local_angles%atom_nr,&
            i=i,value=value,error=error)
    end if failureIf
  end subroutine qs_pao_local_angles_set_p_index
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_l_angles_get_non_p_index [1.0] *
!!
!!   NAME
!!     qs_pao_l_angles_get_non_p_index
!!
!!   FUNCTION
!!     returns the i-te index of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     i: number of the index you want to get
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_l_angles_get_non_p_index(local_angles, i,error) result(res)
    integer ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: i
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_l_angles_get_non_p_index'
    type(qs_pao_env_type), pointer :: pao_env
    failure=.false.
    pao_env => qs_get_pao_env(local_angles)
    if (qs_debug .or. debug_this_module) then
       QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    end if

    failureIf: if (.not. failure) then
       res=qs_pao_env_get_non_p_index(pao_env,atom_nr=local_angles%atom_nr,&
            i=i,error=error)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_l_angles_get_non_p_index
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_local_angles_get_angle [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_get_angle
!!
!!   FUNCTION
!!     returns the requested angle
!!
!!   NOTES
!!     should be fast, and hopefully inlined
!!
!!   INPUTS
!!     local_angles: the local angles object you want to get the angles from
!!     min_bas: the local index of the minimal basis of which you want 
!!            the angle
!!     excl_bas: the local index of the excluded basis of which you want the
!!            angle
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_local_angles_get_angle(local_angles, min_bas, excl_bas, &
       error) result(res)
    real(kind=wp) ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: min_bas, excl_bas
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_local_angles_get_angle'
    type(qs_pao_glob_angles_type), pointer :: glob_angles
    failure=.false.
    glob_angles => local_angles%glob_angles

    debugIf: if (debug_this_module.and.qs_debug) then
       QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(glob_angles),qs_failure_level,routineN,error,failure)
       if (.not.failure) then
          QSPrecondition(associated(glob_angles%pao_env),qs_failure_level,routineN,error,failure)
       end if
       if (.not.failure) then
          call qs_assert(associated(glob_angles%pao_env%nr_min_bas),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "min bas not initialized in "//QSSourceFileRef(),&
               error,failure)
          call qs_assert(associated(glob_angles%pao_env%nr_full_bas),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "full bas not initialized in "//QSSourceFileRef(),&
               error,failure)
       end if
       if (.not.failure) then
          call qs_assert(min_bas>0.and.min_bas <= &
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "min_bas out of range in "//QSSourceFileRef(),&
               error,failure)
          call qs_assert(excl_bas>0.and.excl_bas <= &
               glob_angles%pao_env%nr_full_bas(local_angles%atom_nr)-&
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "excl_bas out of range in "//QSSourceFileRef(),&
               error,failure)
       end if
! also check begins_at???
    end if debugIf
    failureIf: if (.not. failure) then
       res=glob_angles%angles(local_angles%begins_at + &
            glob_angles%pao_env%nr_min_bas(local_angles%atom_nr)*(min_bas-1) + &
            excl_bas)
    else failureIf
       res=-1
    end if failureIf
  end function qs_pao_local_angles_get_angle
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_local_angles_set_angle [1.0] *
!!
!!   NAME
!!     qs_pao_local_angles_set_angle
!!
!!   FUNCTION
!!     returns the requested angle
!!
!!   NOTES
!!     should be fast, and hopefully inlined
!!
!!   INPUTS
!!     local_angles: the local angles object you want to set the angles from
!!     min_bas: the local index of the minimal basis of which you want 
!!            the angle
!!     full_bas: the local index of the full basis of which you want the
!!            angle
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_local_angles_set_angle(local_angles, min_bas, excl_bas, &
       value,error)
    type(qs_pao_local_angles_type), intent(in) :: local_angles !actually the substructure globalangles is modifed
    integer, intent(in) :: min_bas, excl_bas
    real(kind=wp), intent(in) :: value
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_local_angles_set_angle'
    type(qs_pao_glob_angles_type), pointer :: glob_angles
    failure=.false.
    glob_angles => local_angles%glob_angles

    debugIf: if (debug_this_module.and.qs_debug) then
       QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(glob_angles),qs_failure_level,routineN,error,failure)
       if (.not.failure) then
          QSPrecondition(associated(glob_angles%pao_env),qs_failure_level,routineN,error,failure)
       end if
       if (.not.failure) then
          call qs_assert(associated(local_angles%glob_angles%pao_env%nr_min_bas),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "env min bas not initialized in "//QSSourceFileRef(),&
               error,failure)
          call qs_assert(associated(local_angles%glob_angles%pao_env%nr_full_bas),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "env full bas not initialized in "//QSSourceFileRef(),&
               error,failure)
       end if
       if (.not.failure) then
          call qs_assert(min_bas>0.and.min_bas <= &
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "min_bas out of range in "//QSSourceFileRef(),&
               error,failure)
          call qs_assert(excl_bas>0.and.excl_bas <= &
               glob_angles%pao_env%nr_full_bas(local_angles%atom_nr)-&
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               qs_failure_level,qs_assertion_failed,routineN,&
               "excl_bas out of range in "//QSSourceFileRef(),&
               error,failure)
       end if
! also check begins_at???
    end if debugIf
    failureIf: if (.not. failure) then
       glob_angles%angles(local_angles%begins_at + &
            glob_angles%pao_env%nr_min_bas(local_angles%atom_nr)*(min_bas-1) + &
            excl_bas) = value
    end if failureIf
    nullify(glob_angles)
  end subroutine qs_pao_local_angles_set_angle
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_l_angles_get_excl_bas_size [1.0] *
!!
!!   NAME
!!     qs_l_angles_get_excl_bas_size
!!
!!   FUNCTION
!!     returns the size of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  function qs_l_angles_get_excl_bas_size(local_angles, error) result(res)
    integer ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_get_excl_bas_size'
    type(qs_pao_env_type), pointer :: pao_env
    failure=.false.

    QSPrecondition(local_angles%initialized, qs_failure_level,routineN,error,failure)
    pao_env => qs_get_pao_env(local_angles)
    QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       res=qs_pao_env_get_excl_bas_size(pao_env,atom_nr=local_angles%atom_nr,&
            error=error)
    else
       res=-1
    end if
  end function qs_l_angles_get_excl_bas_size
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_l_angles_get_min_bas_size [1.0] *
!!
!!   NAME
!!     qs_l_angles_get_min_bas_size
!!
!!   FUNCTION
!!     returns the size of the minimal basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  function qs_l_angles_get_min_bas_size(local_angles, error) result(res)
    integer ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_get_min_bas_size'
    type(qs_pao_env_type), pointer :: pao_env
    failure=.false.

    QSPrecondition(local_angles%initialized, qs_failure_level,routineN,error,failure)
    pao_env => qs_get_pao_env(local_angles)
    QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       res=qs_pao_env_get_min_bas_size(pao_env,atom_nr=local_angles%atom_nr,&
            error=error)
    else
       res=-1
    end if
  end function qs_l_angles_get_min_bas_size
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_l_angles_get_full_bas_size [1.0] *
!!
!!   NAME
!!     qs_l_angles_get_full_bas_size
!!
!!   FUNCTION
!!     returns the size of the full basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  function qs_l_angles_get_full_bas_size(local_angles, error) result(res)
    integer ::res
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_get_full_bas_size'
    type(qs_pao_env_type), pointer :: pao_env
    failure=.false.

    QSPrecondition(local_angles%initialized, qs_failure_level,routineN,error,failure)
    pao_env => qs_get_pao_env(local_angles)
    QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       res=qs_pao_env_get_full_bas_size(pao_env,atom_nr=local_angles%atom_nr,&
            error=error)
    else
       res=-1
    end if
  end function qs_l_angles_get_full_bas_size
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_l_angles_get_glob_angles [1.0] *
!!
!!   NAME
!!     qs_l_angles_get_glob_angles
!!
!!   FUNCTION
!!     returns the global angles object that this object accesses
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get the info from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_l_angles_get_glob_angles(local_angles,error) result(res)
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_pao_glob_angles_type), pointer ::res
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_get_glob_angles'
    failure=.false.

    QSPrecondition(qs_valid(local_angles),qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       res => local_angles%glob_angles
    else
       nullify(res)
    end if
  end function qs_l_angles_get_glob_angles
!!***
!***************************************************************************


!!****f* qs_pao_types/qs_l_angles_did_change [1.0] *
!!
!!   NAME
!!     qs_l_angles_did_change
!!
!!   FUNCTION
!!     invalidate the cached unitary transformation of this atom.
!!     to be called when the angles change
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the info from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_l_angles_did_change(local_angles,error)
! changes in local_angles%glob_angles
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_did_change'
    logical, dimension(:), pointer :: ut_calculated
    failure=.false.

    QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(local_angles%atom_nr>0,qs_warning_level,routineN,error,failure)
    if (.not. failure) then
       if (associated(local_angles%glob_angles%ut_calculated)) then
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.false.
       end if
    end if
  end subroutine qs_l_angles_did_change
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_l_angles_set_unitary_t [1.0] *
!!
!!   NAME
!!     qs_l_angles_set_unitary_t
!!
!!   FUNCTION
!!     changes the cached unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you want to change
!!     unitary_t: the new unitary trasformation
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_l_angles_set_unitary_t(local_angles,unitary_t,error)
! changes in glob_angles
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    real(kind=wp), dimension(:,:), intent(in), target :: unitary_t
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_set_unitary_t'
    type(qs_pao_matrix_type), pointer :: full_u_t
    real(kind=wp), dimension(:,:), pointer :: unitary_t_ptr
    failure=.false.

    QSPrecondition(local_angles%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(local_angles%atom_nr>0,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       call qs_pao_glob_angles_get(local_angles%glob_angles,unitary_t=full_u_t,&
            error=error)
       QSPrecondition(associated(full_u_t),qs_warning_level,routineN,error,failure)
       call qs_assert(associated(local_angles%glob_angles%ut_calculated),&
            qs_failure_level,qs_assertion_failed,routineN,&
            "ut_calculated not allocated in "//&
            QSSourceFileRef(),error,failure)
       if (.not.failure) then
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.true.
          unitary_t_ptr => unitary_t
          call put_block_node(qs_get_matrix(full_u_t),&
               block_row=local_angles%atom_nr,&
               block_col=local_angles%atom_nr, block=unitary_t_ptr)
       end if
    end if
  end subroutine qs_l_angles_set_unitary_t
!!***
!***************************************************************************

!!****f* qs_pao_plain_rot/qs_pao_l_compose_u [1.0] *
!!
!!   NAME
!!     qs_pao_l_compose_u
!!
!!   FUNCTION
!!     calculates U*matrix where U is the unitary transformation defined by
!!     the given angles. the sequence is: smallbasis then full basis in 
!!     (local) growing order: u(1,1)*u(2,1)*...*u(size_small_bas,1)*...*
!!      u(size_small_bas,size_full_bas)*matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that will be modified
!!     local_angles: the angles that define U
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_l_compose_u(matrix,local_angles,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_plain_rot:qs_pao_l_compose_u'
    integer :: n_min_bas, n_excl_bas, i_min_bas,i_excl_bas,n_full_bas
    failure=.false.

    call qs_get(local_angles,n_excl_bas=n_excl_bas,n_min_bas=n_min_bas,&
         n_full_bas=n_full_bas)
    QSPrecondition(n_min_bas>0,qs_warning_level,routineN,error,failure)
    QSPrecondition(n_excl_bas>=0,qs_warning_level,routineN,error,failure)
    QSPrecondition(size(matrix,1)==n_full_bas,qs_warning_level,routineN,error,failure)
    QSPrecondition(size(matrix,2)==n_full_bas,qs_warning_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       do i_excl_bas=n_excl_bas,1,-1
          do i_min_bas=n_min_bas,1,-1
             call qs_left_compose_p_rot(matrix=matrix,&
                  angle=qs_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=qs_get_p_index(local_angles,i=i_min_bas,error=error),&
                  j=qs_get_non_p_index(local_angles,i=i_excl_bas,error=error),&
                  error=error)
          end do
       end do
    end if failureIf
  end subroutine qs_pao_l_compose_u
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_l_angles_get_unitary_t [1.0] *
!!
!!   NAME
!!     qs_l_angles_get_unitary_t
!!
!!   FUNCTION
!!     returns the unitary trasformation defined by these angles
!!
!!   NOTES
!!     the trasformation is cached
!!
!!   INPUTS
!!     local_angles: the local angles you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_l_angles_get_unitary_t(local_angles,error) result(res)
    type(qs_pao_local_angles_type), intent(in) :: local_angles
    type(qs_error_type), optional, intent(inout) :: error
    real(kind=wp), dimension(:,:), pointer :: res

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_l_angles_get_unitary_t'
    integer :: i
    type(qs_pao_matrix_type), pointer :: unitary_t
    failure=.false.
    nullify(res)

    QSPrecondition(qs_valid(local_angles),qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       call qs_get(qs_get_glob_angles(local_angles),unitary_t=unitary_t,&
            error=error)
       QSPrecondition(associated(unitary_t),qs_failure_level,routineN,error,failure)
       call qs_assert(associated(local_angles%glob_angles%ut_calculated),&
            qs_failure_level,qs_assertion_failed,routineN,&
            "ut_calculated not associated in "//&
            QSSourceFileRef(),error,failure)
       if (.not.failure) then
          if (local_angles%glob_angles%ut_calculated(local_angles%atom_nr)) then
             call get_block_node(qs_get_matrix(unitary_t),&
                  block=res,block_row=local_angles%atom_nr,&
                  block_col=local_angles%atom_nr)
             QSPostcondition(associated(res),qs_failure_level,routineN,error,failure)
          else
             call get_block_node(qs_get_matrix(unitary_t),&
                  block=res,block_row=local_angles%atom_nr,&
                  block_col=local_angles%atom_nr)
             if (.not.associated(res)) then
                call add_block_node(qs_get_matrix(unitary_t),&
                     local_angles%atom_nr,local_angles%atom_nr,block=res)
             end if
             QSPostcondition(associated(res),qs_failure_level,routineN,error,failure)
             if (.not.failure) then
                QSPrecondition(size(res,1)==size(res,2),qs_failure_level,routineN,error,failure)
             end if
             if (.not.failure) then
                res=0.0_wp
                forall (i=1:size(res,1)) res(i,i)=1.0_wp
                call qs_pao_l_compose_u(res,local_angles,error=error)
                local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.true.
             end if
          end if
       end if
    end if
    if (failure) nullify(res)
  end function qs_l_angles_get_unitary_t
!!***
!***************************************************************************

!end
! ==================== qs_pao_matrix_type ====================
!begin

! ====== core methods ======

!!****f* qs_pao_types/qs_init_pao_matrix [1.0] *
!!
!!   NAME
!!     qs_init_pao_matrix
!!
!!   FUNCTION
!!     initializes a pao matrix. If no matrix is given allocates a block
!!     diagonal matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!     For the others see qs_set_pao_matrix (that is called
!!     by this routine).
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_init_pao_matrix(pao_matrix, pao_env,&
       matrix, source, min_bas_rows, min_bas_cols,&
       should_dealloc_matrix, error)
    type(qs_pao_matrix_type), intent(out) :: pao_matrix
    type(qs_pao_env_type), target :: pao_env
    type(real_matrix_type), optional, target :: matrix
    type(qs_pao_matrix_type), intent(in), optional :: source
    logical, intent(in), optional :: min_bas_cols, min_bas_rows, &
         should_dealloc_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_init_pao_matrix'
    failure=.false.

    pao_matrix%initialized=.true.
    pao_matrix%should_dealloc_matrix=.true.
    nullify(pao_matrix%matrix)
! nullify(pao_matrix%pao_env)
    pao_matrix%pao_env => pao_env
    pao_matrix%min_bas_rows=.true.
    pao_matrix%min_bas_rows=.true.
    call qs_set_pao_matrix(pao_matrix,pao_env,&
         matrix, source, min_bas_rows, min_bas_cols,&
         should_dealloc_matrix, error)
  end subroutine qs_init_pao_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_dealloc_pao_matrix_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_matrix_ref
!!
!!   FUNCTION
!!     Dealloc the memory that the matrix has allocated
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that should dealloc its memory
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_pao_matrix_ref(pao_matrix, error)
    type(qs_pao_matrix_type), intent(inout) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_dealloc_pao_matrix_ref'
    failure=.false.

    QSPrecondition(pao_matrix%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       nullify(pao_matrix%pao_env)
       if (pao_matrix%should_dealloc_matrix) then
          if (associated(pao_matrix%matrix)) then
             call deallocate_matrix(pao_matrix%matrix)
          end if
       else
          nullify(pao_matrix%matrix)
       end if
    end if failureIf
  end subroutine qs_dealloc_pao_matrix_ref
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_create_pao_matrix [1.0] *
!!
!!   NAME
!!     qs_create_pao_matrix
!!
!!   FUNCTION
!!     allocates the matrix and initializes it
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: will point to the newly created matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!     for the others arguments see qs_init_block_pao_matrix
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_create_pao_matrix(pao_matrix, pao_env,&
       matrix, source, min_bas_rows, min_bas_cols,&
       should_dealloc_matrix, error)
    type(qs_pao_matrix_type), pointer :: pao_matrix
    type(qs_pao_env_type), target :: pao_env
    type(real_matrix_type), optional, target :: matrix
    type(qs_pao_matrix_type), intent(in), optional :: source
    logical, intent(in), optional :: min_bas_cols, min_bas_rows, &
         should_dealloc_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_create_pao_matrix'
    integer :: stat
    failure=.false.

    allocate(pao_matrix, stat=stat)
    QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       call qs_init_pao_matrix(pao_matrix, pao_env,&
            matrix, source, min_bas_rows, min_bas_cols,&
            should_dealloc_matrix, error)
    end if
  end subroutine qs_create_pao_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_dealloc_pao_matrix [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_matrix
!!
!!   FUNCTION
!!     deallocates the memory allocated by the matrix ans the matrix
!!     object itself.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_dealloc_pao_matrix(pao_matrix, error)
    type(qs_pao_matrix_type), pointer :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_dealloc_pao_matrix'
    integer :: stat
    failure=.false.

    if (associated(pao_matrix)) then
       call qs_dealloc_pao_matrix_ref(pao_matrix,error=error)
       deallocate(pao_matrix,stat=stat)
       QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
    end if
  end subroutine qs_dealloc_pao_matrix
!!***
!***************************************************************************


!!****f* qs_pao_types/qs_validate_pao_matrix [1.0] *
!!
!!   NAME
!!     qs_validate_pao_matrix
!!
!!   FUNCTION
!!     controls if the matrix  has the right size, returns true
!!     if it is valid, false if it is not.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix to be tested
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_validate_pao_matrix(pao_matrix, error) result(res)
    logical ::res
    type(qs_pao_matrix_type), intent(in) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_validate_pao_matrix'
    integer, dimension(:), pointer :: cols_size, rows_size
    failure=.false.
    nullify(cols_size,rows_size)

    QSPrecondition(pao_matrix%initialized, qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       QSAssert(associated(pao_matrix%matrix),qs_warning_level,routineN,error,failure)
       QSAssert(associated(pao_matrix%pao_env), qs_warning_level,routineN,error,failure)
       if (.not.failure) then
          QSAssert(associated(pao_matrix%pao_env%nr_full_bas),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_matrix%pao_env%nr_min_bas),qs_warning_level,routineN,error,failure)
       end if
       if (.not.failure) then
          call qs_pao_matrix_get(pao_matrix, cols_size=cols_size,&
               rows_size=rows_size)
          call qs_assert(qs_consistent_matrix(matrix=pao_matrix%matrix,&
               rows_size=rows_size, cols_size=cols_size,error=error),&
               qs_warning_level, qs_assertion_failed, routineN,&
               "size in inconsistent with given basis! in "//&
               QSSourceFileRef(),error, failure)
       end if
    end if failureIf
    res = .not.failure
  end function qs_validate_pao_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_set_pao_matrix [1.0] *
!!
!!   NAME
!!     qs_set_pao_matrix
!!
!!   FUNCTION
!!     sets the values of the various elements of the stucture
!!
!!   NOTES
!!     you should rarely nedd this (initialize in init)
!!
!!   INPUTS
!!     pao_matrix: the transformation you want to change
!!     pao_env: the pao environement this transformation refers to.
!!           it is SHARED, and it should be already valid.
!!     matrix: normally COPIES the matrix into this object, but if at the
!!             same time you also set should_dealloc_matrix to false
!!             (i.e. the matrix is shared) it is SHARED.
!!     source: a transformation that is copied into this one(other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     min_bas_rows: if it is true (default) the rows are in the 
!!             minimal basis
!!     min_bas_cols: if it is true (default) the cols are in the
!!             minimal basis
!!     should_dealloc_matrix: if true (default) the matrix will be 
!!             deallocated in qs_dealloc_ref
!!     error: variable to control error logging, stopping,... 
!!           see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_set_pao_matrix(pao_matrix, pao_env,&
       matrix, source, min_bas_rows, min_bas_cols,&
       should_dealloc_matrix, error)
    type(qs_pao_matrix_type), intent(inout) :: pao_matrix
    type(qs_pao_env_type), optional, target :: pao_env
    type(real_matrix_type), optional, target :: matrix
    type(qs_pao_matrix_type), intent(in), optional :: source
    logical, intent(in), optional :: min_bas_cols, min_bas_rows,&
         should_dealloc_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure, should_share_matrix
    character(len=*), parameter :: routineN='qs_pao_types:qs_set_pao_matrix'
    integer :: stat
    character (len=60) :: matrix_name
    type(real_matrix_type), pointer :: tmp_ptr
    failure=.false.
    nullify(tmp_ptr)

    QSPrecondition(pao_matrix%initialized, qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       sourceIf: if (present(source)) then
          QSPrecondition(source%initialized, qs_failure_level,routineN,error,failure)
          if (present(should_dealloc_matrix)) then
             should_share_matrix=.not.should_share_matrix
          else
             should_share_matrix=.false.
          end if
          if (.not.failure) then
             pao_matrix%pao_env => source%pao_env
             pao_matrix%min_bas_rows=source%min_bas_rows
             pao_matrix%min_bas_cols=source%min_bas_cols
             if (associated(source%matrix)) then
                if (.not.present(should_dealloc_matrix)) &
                     should_share_matrix=.not.source%should_dealloc_matrix
                if (.not.should_share_matrix) then
                   if (.not. associated(pao_matrix%matrix)) then
                      call get_matrix_info(source%matrix,matrix_name=matrix_name)
                      call replicate_matrix(source%matrix, &
                           pao_matrix%matrix,&
                           matrix_name)
                   else
                      call copy_matrix(source%matrix, &
                           pao_matrix%matrix)
                   end if
                else ! the matrix is shared
                   pao_matrix%matrix => source%matrix
                end if
             else if (associated(pao_matrix%matrix))then
                if (pao_matrix%should_dealloc_matrix) then
                   deallocate(pao_matrix%matrix, stat=stat)
                   QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
                else
                   nullify(pao_matrix%matrix)
                end if
             end if
          end if
          pao_matrix%should_dealloc_matrix=source%should_dealloc_matrix
       end if sourceIf
       if (present(pao_env)) then
          pao_matrix%pao_env => pao_env
       end if
       matrixIf: if (present(matrix)) then
          if (should_share_matrix) then
             pao_matrix%matrix => matrix
          else
             tmp_ptr => matrix
             if (.not. associated(pao_matrix%matrix)) then
                call get_matrix_info(source%matrix,matrix_name=matrix_name)
                call replicate_matrix(tmp_ptr, &
                     pao_matrix%matrix,&
                     matrix_name)
             else
                call copy_matrix(tmp_ptr, &
                     pao_matrix%matrix)
             end if
          end if
       end if matrixIf
       if (present(min_bas_rows)) &
            pao_matrix%min_bas_rows=min_bas_rows
       if (present(min_bas_cols)) &
            pao_matrix%min_bas_cols=min_bas_cols
       if (present(should_dealloc_matrix)) &
            pao_matrix%should_dealloc_matrix=should_dealloc_matrix
       if (associated(pao_matrix%pao_env).and. &
            .not. associated(pao_matrix%matrix)) then
          call qs_init_block_pao_matrix(pao_matrix, error=error)
       end if
    end if failureIf
  end subroutine qs_set_pao_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_matrix_get [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_get
!!
!!   FUNCTION
!!     Get the attributes of the given matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the pao transformation you want info from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!     For the other arguments see the attributes of qs_pao_matrix_type
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_matrix_get(pao_matrix,pao_env,matrix,min_bas_cols,&
       min_bas_rows, cols_size, rows_size, should_dealloc_matrix, error)
    type(qs_pao_matrix_type), intent(in) :: pao_matrix
    type(qs_pao_env_type), optional, pointer :: pao_env
    type(real_matrix_type), optional, pointer :: matrix
    logical, optional, intent(out) :: min_bas_rows, min_bas_cols, &
         should_dealloc_matrix
    integer, dimension(:), pointer, optional :: cols_size, rows_size
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_matrix_get'
    failure=.false.

    QSPrecondition(pao_matrix%initialized,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (present(pao_env)) pao_env => pao_matrix%pao_env
       if (present(matrix)) matrix => pao_matrix%matrix
       if (present(min_bas_cols)) min_bas_cols=pao_matrix%min_bas_cols
       if (present(min_bas_rows)) min_bas_rows=pao_matrix%min_bas_rows
       if (present(should_dealloc_matrix)) &
            should_dealloc_matrix=pao_matrix%should_dealloc_matrix
       if (present(cols_size)) then
          QSPrecondition(associated(pao_matrix%pao_env),qs_warning_level,routineN,error,failure)
          if (.not.failure) then
             if (pao_matrix%min_bas_cols) then
                cols_size => pao_matrix%pao_env%nr_min_bas
             else
                cols_size => pao_matrix%pao_env%nr_full_bas
             end if
          else
             nullify(cols_size)
          end if
       end if
       if (present(rows_size)) then
          QSPrecondition(associated(pao_matrix%pao_env),qs_warning_level,routineN,error,failure)
          if (.not.failure) then
             if (pao_matrix%min_bas_rows) then
                rows_size => pao_matrix%pao_env%nr_min_bas
             else
                rows_size => pao_matrix%pao_env%nr_full_bas
             end if
          else
             nullify(rows_size)
          end if
       end if
    end if failureIf
  end subroutine qs_pao_matrix_get
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_matrix_valid [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_valid
!!
!!   FUNCTION
!!     performs minimal validadion (no dangling pointers)
!!
!!   NOTES
!!     can be called often
!!
!!   INPUTS
!!     pao_matrix: the object to validate
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_matrix_valid(pao_matrix, error) result(res)
    logical ::res
    type(qs_pao_matrix_type), intent(in) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_matrix_valid'
    failure=.false.

    QSPrecondition(pao_matrix%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (.not.(associated(pao_matrix%matrix).and.&
            associated(pao_matrix%pao_env)) ) then
          QSAssert(associated(pao_matrix%matrix),qs_warning_level,routineN,error,failure)
          QSAssert(associated(pao_matrix%pao_env),qs_warning_level,routineN,error,failure)
       end if
    end if
    res=.not.failure
  end function qs_pao_matrix_valid
!!***
!***************************************************************************

! ===== special methods =====

!!****f* qs_pao_types/qs_pao_matrix_set_to_unity [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_set_to_unity
!!
!!   FUNCTION
!!     sets the matrix to the unit matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_matrix_set_to_unity(pao_matrix, error)
    type(qs_pao_matrix_type), intent(inout) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_matrix_set_to_unity'
    type(qs_error_type) :: iError
    failure=.false.
    call qs_init(iError,template_error=error)

    QSPrecondition(pao_matrix%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       if (.not.associated(pao_matrix%matrix)) then
          call qs_init_block_pao_matrix(pao_matrix, error=iError)
          call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
               error=error,failure=failure)
       end if
    end if
    if (.not. failure) then
       call qs_set_to_identity(pao_matrix%matrix)
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_pao_matrix_set_to_unity
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_matrix_get_pao_env [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_get_pao_env
!!
!!   FUNCTION
!!     returns the global environement the matrix object lives in
!!
!!   NOTES
!!     this might go away (to optimize), so do not rely on it too heavily
!!
!!   INPUTS
!!     pao_matrix: the matrix object you get the info from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_matrix_get_pao_env(pao_matrix, error) result(res)
    type(qs_pao_env_type), pointer ::res
    type(qs_pao_matrix_type), intent(in) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_matrix_get_pao_env'
    failure=.false.
    QSPrecondition(pao_matrix%initialized, qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       res => pao_matrix%pao_env
    else failureIf
       nullify(res)
    end if failureIf
  end function qs_pao_matrix_get_pao_env
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_matrix_get_matrix [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_get_matrix
!!
!!   FUNCTION
!!     returns the undelying matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix you get the info from
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_pao_matrix_get_matrix(pao_matrix, error) result(res)
    type(real_matrix_type), pointer ::res
    type(qs_pao_matrix_type), intent(in) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_matrix_get_matrix'
    failure=.false.
    QSPrecondition(pao_matrix%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       res => pao_matrix%matrix
    else
       nullify(res)
    end if
  end function qs_pao_matrix_get_matrix
!!***
!***************************************************************************


! ===== private utilities =====

!!****f* qs_pao_types/qs_init_block_pao_matrix [1.0] *
!!
!!   NAME
!!     qs_init_block_pao_matrix
!!
!!   FUNCTION
!!     if the matrix is not already right
!!     allocate a blocks sizes of the matrix using the pao_env info
!!     but no blocks are allocated (use qs_set_to_identity to have a 
!!     block diagonal matrix).
!!
!!   NOTES
!!     pao_env should be valid.
!!     private method
!!
!!   INPUTS
!!     pao_matrix: the trasformation that should be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_init_block_pao_matrix(pao_matrix, error)
    type(qs_pao_matrix_type), intent(inout) :: pao_matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_init_block_pao_matrix'
    integer :: n_atom, i, stat
    integer, dimension(:), allocatable :: firstRow, lastRow, firstCol, lastCol
    integer, dimension(:), pointer :: cols_size, rows_size
    failure=.false.
    n_atom=-1

    QSPrecondition(pao_matrix%initialized, qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(pao_matrix%pao_env),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_assert(qs_validate_pao_env(pao_matrix%pao_env, &
            test_pos_att=.false., error=error), qs_failure_level, &
            qs_precondition_failed, routineN,& 
            'pao_env must be valid to initialize a pao matrix in '//&
            QSSourceFileRef(),error,failure) !remove this test?
    end if
    if (.not.failure) then
       call qs_get(pao_matrix, cols_size=cols_size, rows_size=rows_size)
       QSPrecondition(associated(cols_size),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(rows_size),qs_failure_level,routineN,error,failure)
    end if
    failureIf: if (.not. failure) then
       if (associated(pao_matrix%matrix)) then
          if (.not.qs_consistent_matrix(pao_matrix%matrix,&
               cols_size=cols_size,&
               rows_size=rows_size, error=error)) then
             QSErrorMessage(qs_warning_level,routineN,"the matrix was inconsistent",error)
             call deallocate_matrix(pao_matrix%matrix)
          end if
       end if
       if (.not.associated(pao_matrix%matrix)) then
          n_atom=qs_get_n_atom(pao_matrix%pao_env)
          QSPrecondition(n_atom==size(cols_size),qs_failure_level,routineN,error,failure)
          QSPrecondition(n_atom==size(rows_size),qs_failure_level,routineN,error,failure)
          if (.not.failure) then
             allocate(firstRow(n_atom), lastRow(n_atom),&
                  firstCol(n_atom), lastCol(n_atom), stat=stat)
             QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
          end if
          if (.not.failure) then
             call qs_first_last_from_size(cols_size,first_index=firstCol,&
                  last_index=lastCol,error=error)
             call qs_first_last_from_size(rows_size,first_index=firstRow,&
                  last_index=lastRow,error=error)
             call allocate_matrix(matrix=pao_matrix%matrix, &
                  nblock_row=n_atom, nblock_col=n_atom, &
                  nrow=sum(rows_size), ncol=sum(cols_size), &
                  first_row=firstRow, last_row=lastRow, &
                  first_col=firstCol, last_col=lastCol, &
                  matrix_name="pao matrix", &
                  matrix_symmetry="block diagonal")
             if (debug_this_module .and. qs_debug) then
                call qs_assert(&
                     qs_consistent_matrix(pao_matrix%matrix,&
                     cols_size=cols_size, rows_size=rows_size,error=error),&
                     qs_failure_level,qs_assertion_failed,&
                     routineN,"internal error: after its init the pao "//&
                     "transformation is not consistent in "//&
                     QSSourceFileRef(), error, failure)
             end if
             deallocate(firstRow,lastRow, firstCol, lastCol, stat=stat) 
! just to be sure...
             QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
          end if
       end if
    end if failureIf
  end subroutine qs_init_block_pao_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_first_last_from_size [1.0] *
!!
!!   NAME
!!     qs_first_last_from_size
!!
!!   FUNCTION
!!     generates the first index and last index from the sizes of the blocks
!!
!!   NOTES
!!     private method
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_first_last_from_size(block_size, first_index, last_index, error)
    type(qs_error_type), optional, intent(inout) :: error
    integer, dimension(:), intent(out) :: first_index, last_index
    integer, dimension(:), intent(in) :: block_size
    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_first_last_from_size'
    integer :: n_blocks, i
    failure=.false.
    n_blocks=size(block_size)

    QSPrecondition(size(first_index)==n_blocks,qs_failure_level,routineN,error,failure)
    QSPrecondition(size(last_index)==n_blocks,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (n_blocks>0) then
          first_index(1)=1
       end if
       do i=1,n_blocks-1
          last_index(i)=first_index(i)+block_size(i)-1
          first_index(i+1)=last_index(i)+1
       end do
       if (n_blocks>0) then
          last_index(n_blocks)=first_index(n_blocks)+block_size(n_blocks)
       end if
    end if
  end subroutine qs_first_last_from_size
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_consistent_matrix [1.0] *
!!
!!   NAME
!!     qs_consistent_matrix
!!
!!   FUNCTION
!!     Returns true if the matrix size is consistent
!!     with the given sizes
!!
!!   NOTES
!!     private method
!!
!!   INPUTS
!!     matrix: the matrix to check
!!     cols_size: the size of the columns of the blocks
!!     rows_size the size of the rows of the blocks
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  function qs_consistent_matrix(matrix,cols_size,rows_size, error) result(res)
    logical ::res
    type(real_matrix_type), intent(in), target :: matrix
    integer, dimension(:), intent(in) :: cols_size, rows_size
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_consistent_matrix'
    integer :: n_atom, i, stat
    integer, dimension(size(cols_size)) :: firstRow, lastRow,&
         firstCol,lastCol
    integer, dimension(:), pointer :: nr_full_bas, first_col, last_col, &
         first_row, last_row
    type(real_matrix_type), pointer :: tstPtr
    failure=.false.
    nullify(nr_full_bas, first_col, last_col, first_row, last_row)
    tstPtr => matrix

    QSPrecondition(associated(tstPtr),qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       n_atom=size(cols_size)
       QSPrecondition(n_atom==size(rows_size),qs_failure_level,routineN,error,failure)
       if (.not.failure) then
          call qs_first_last_from_size(cols_size,first_index=firstCol,&
               last_index=lastCol)
          call qs_first_last_from_size(rows_size,first_index=firstRow,&
               last_index=lastRow)
          call get_matrix_info(tstPtr, &
               first_col=first_col , last_col=last_col, &
               first_row=first_row , last_row=last_row)
          QSAssert(associated(first_row).and.associated(last_row),qs_warning_level,routineN,error,failure)
          QSAssert(associated(first_col).and.associated(last_row),qs_warning_level,routineN,error,failure)
          if (.not.failure) then
             QSAssert(size(first_row)==n_atom,qs_warning_level,routineN,error,failure)
             QSAssert(size(last_row)==n_atom,qs_warning_level,routineN,error,failure)
             QSAssert(size(first_col)==n_atom,qs_warning_level,routineN,error,failure)
             QSAssert(size(last_col)==n_atom,qs_warning_level,routineN,error,failure)
             if (.not.failure) then
                QSAssert(all(first_row==firstRow),qs_warning_level,routineN,error,failure)
                QSAssert(all(first_col==firstCol),qs_warning_level,routineN,error,failure)
                QSAssert(all(last_row==lastRow),qs_warning_level,routineN,error,failure)
                QSAssert(all(last_col==lastCol),qs_warning_level,routineN,error,failure)
             end if
          end if
          nullify(first_col, last_col, first_row, last_row)
       end if
    end if failureIf
    res= .not. failure
  end function qs_consistent_matrix
!!***
!***************************************************************************

!end
! ======================= qs_pao_projection =====================
!begin

! this routines WILL CHANGE, they are extremly ugly and inefficient
! I will come back to them as soo as the programs works
! to do

!!****f* qs_pao_projection/qs_pao_l_atom_projection [1.0] *
!!
!!   NAME
!!     qs_pao_l_atom_projection
!!
!!   FUNCTION
!!     Left multiplication with the projection in the small basis
!!     in one atom block.
!!     Puts pi*source_m into target_m, where pi is the projection in the
!!     minimal basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_l_atom_projection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_l_atom_projection'
    integer :: size_min_bas, size_full_bas, i
    failure=.false.

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       size_full_bas=qs_get_full_bas_size(pao_env,atom_nr=atom_nr,error=error)
       size_min_bas=qs_get_min_bas_size(pao_env,atom_nr=atom_nr,error=error)
       QSPrecondition(size(source_m,1)==size_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,1)==size_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(source_m,2)==size(target_m,2),qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       do i=1,size_min_bas
          target_m(i,:)=source_m(qs_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error),:)
! reorganization of the operations could probabily make it faster
       end do
    end if
  end subroutine qs_pao_l_atom_projection
!!***
!***************************************************************************

!!****f* qs_pao_projection/qs_pao_l_atom_injection [1.0] *
!!
!!   NAME
!!     qs_pao_l_atom_injection
!!
!!   FUNCTION
!!     Left multiplication with the injection in the full basis
!!     in one atom block.
!!     Puts i*source_m into target_m, where i is the injection in the
!!     full basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_l_atom_injection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_l_atom_injection'
    integer :: size_min_bas, size_full_bas, i
    failure=.false.

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       size_full_bas=qs_get_full_bas_size(pao_env,atom_nr=atom_nr,error=error)
       size_min_bas=qs_get_min_bas_size(pao_env,atom_nr=atom_nr,error=error)
       QSPrecondition(size(source_m,1)==size_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,1)==size_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(source_m,2)==size(target_m,2),qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       target_m=0.0_wp
       do i=1,size_min_bas
          target_m(qs_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error),:)=&
               source_m(i,:)
! reorganization of the operations could probabily make it faster
! and do not set target_m to 0 if it gets overwritten
       end do
    end if
  end subroutine qs_pao_l_atom_injection
!!***
!***************************************************************************

!!****f* qs_pao_projection/qs_pao_r_atom_projection [1.0] *
!!
!!   NAME
!!     qs_pao_r_atom_projection
!!
!!   FUNCTION
!!     Right multiplication with the projection in the small basis
!!     in one atom block.
!!     Puts source_m*pi into target_m, where pi is the projection in the
!!     minimal basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_r_atom_projection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_r_atom_projection'
    integer :: size_min_bas, size_full_bas, i
    failure=.false.

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       size_full_bas=qs_get_full_bas_size(pao_env,atom_nr=atom_nr,error=error)
       size_min_bas=qs_get_min_bas_size(pao_env,atom_nr=atom_nr,error=error)
       QSPrecondition(size(source_m,2)==size_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,2)==size_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(source_m,1)==size(target_m,1),qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       target_m=0.0_wp
       do i=1,size_min_bas
          target_m(:,qs_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error))=&
               source_m(:,i)
! reorganization of the operations could probabily make it faster
       end do
    end if
  end subroutine qs_pao_r_atom_projection
!!***
!***************************************************************************

!!****f* qs_pao_projection/qs_pao_r_atom_injection [1.0] *
!!
!!   NAME
!!     qs_pao_r_atom_injection
!!
!!   FUNCTION
!!     Right multiplication with the injection in the full basis
!!     in one atom block.
!!     Puts source_m*i into target, where i is the injection in the
!!     full basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     atom_nr: the actual atom block
!!     source_m: the source matrix
!!     target_m: the target matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_r_atom_injection(pao_env, atom_nr, source_m,&
       target_m, error)
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_r_atom_injection'
    integer :: size_min_bas, size_full_bas, i
    type(qs_error_type) :: iError
    failure=.false.
    call qs_init(iError,template_error=error)

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       size_full_bas=qs_get_full_bas_size(pao_env,atom_nr=atom_nr,error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_min_bas=qs_get_min_bas_size(pao_env,atom_nr=atom_nr,error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       QSPrecondition(size(source_m,2)==size_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,2)==size_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(source_m,1)==size(target_m,1),qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       do i=1,size_min_bas
          target_m(:,i) = source_m(:,&
               qs_get_p_index(pao_env,atom_nr=atom_nr,i=i,error=error))
! reorganization of the operations could probabily make it faster
       end do
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_pao_r_atom_injection
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_atom_m_to_f_matrix [1.0] *
!!
!!   NAME
!!     qs_pao_atom_m_to_f_matrix
!!
!!   FUNCTION
!!     inject an atomic matrix block into the full basis:
!!     puts i*source_m*pi into target_m, where i the injection, and 
!!     pi(=i^*) is the projection in the small basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     l_atom_nr: the left atom block
!!            (pao_env%nr_min_bas(l_atom_nr)= rows of source_m )
!!     r_atom_nr: the right atom block
!!            (pao_env%nr_min_bas(r_atom_nr)= colss of source_m )
!!     source_m: the source matrix
!!     target_m: the target matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_atom_m_to_f_matrix(pao_env,l_atom_nr,r_atom_nr,&
       source_m, target_m, error)
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: r_atom_nr, l_atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_atom_m_to_f_matrix'
    integer :: size_l_min_bas, size_l_full_bas, size_r_min_bas, size_r_full_bas,&
         i,j
    integer, dimension(:), pointer :: l_index, r_index
    type(qs_error_type) :: iError
    failure=.false.
    nullify(l_index,r_index)
    call qs_init(iError,template_error=error)

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       size_l_full_bas=qs_get_full_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_l_min_bas=qs_get_min_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_r_full_bas=qs_get_full_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_r_min_bas=qs_get_min_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       QSPrecondition(size(source_m,1)==size_l_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(source_m,2)==size_r_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,1)==size_l_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,2)==size_r_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%projection),qs_failure_level,routineN,error,failure)
! check projection more?
    end if
    if (.not. failure) then
       target_m=0.0_wp
       l_index => pao_env%projection(l_atom_nr)%array
       r_index => pao_env%projection(r_atom_nr)%array
       do i=1,size_l_min_bas
          do j=1,size_r_min_bas
             target_m(l_index(i),r_index(j))=source_m(i,j)
          end do
       end do
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_pao_atom_m_to_f_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_atom_f_to_m_matrix [1.0] *
!!
!!   NAME
!!     qs_pao_atom_f_to_m_matrix
!!
!!   FUNCTION
!!     project an atomic matrix block into the minimal basis:
!!     puts pi*source_m*i into target_m, where i the injection, and 
!!     pi(=i^*) is the projection in the small basis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao env where the projection is defined
!!     l_atom_nr: the left atom block 
!!            (pao_env%nr_full_bas(l_atom_nr)= rows of source_m )
!!     r_atom_nr: the right atom block
!!            (pao_env%nr_full_bas(r_atom_nr)= colss of source_m )
!!     source_m: the source matrix
!!     target_m: the target matrix
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_atom_f_to_m_matrix(pao_env,l_atom_nr,r_atom_nr,&
       source_m, target_m, error)
    type(qs_pao_env_type), intent(in) :: pao_env
    integer, intent(in) :: r_atom_nr, l_atom_nr
    real(kind=wp), dimension(:,:), intent(in) :: source_m
    real(kind=wp), dimension(:,:), intent(out) :: target_m
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_atom_f_to_m_matrix'
    integer :: size_l_min_bas, size_l_full_bas, size_r_min_bas, size_r_full_bas,&
         i,j
    integer, dimension(:), pointer :: l_index, r_index
    type(qs_error_type) :: iError
    failure=.false.
    nullify(l_index,r_index)
    call qs_init(iError,template_error=error)

    QSPrecondition(pao_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       size_l_full_bas=qs_get_full_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_l_min_bas=qs_get_min_bas_size(pao_env,atom_nr=l_atom_nr,&
            error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_r_full_bas=qs_get_full_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       size_r_min_bas=qs_get_min_bas_size(pao_env,atom_nr=r_atom_nr,error=iError)
       QSAssert(qs_get_error_level(iError)==0,qs_get_error_level(iError),routineN,error,failure)
       QSPrecondition(size(source_m,1)==size_l_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(source_m,2)==size_r_min_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,1)==size_l_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(size(target_m,2)==size_r_full_bas,qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%projection),qs_failure_level,routineN,error,failure)
! check projection more?
    end if
    if (.not. failure) then
       l_index => pao_env%projection(l_atom_nr)%array
       r_index => pao_env%projection(r_atom_nr)%array
       do i=1,size_l_min_bas
          do j=1,size_r_min_bas
             target_m(i,j)=source_m(l_index(i),r_index(j))
          end do
       end do
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_pao_atom_f_to_m_matrix
!!***
!***************************************************************************

!!****f* qs_pao_types/qs_pao_matrix_transf_to_m [1.0] *
!!
!!   NAME
!!     qs_pao_matrix_transf_to_m
!!
!!   FUNCTION
!!     Transform the given matrix to the minimal basis.
!!     Puts pi U^* source_m U i into the target_m
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: matrix in the full basis
!!     target_m: the result matrix in the small basis
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  subroutine qs_pao_matrix_transf_to_m(source_m, target_m, glob_angles, error)
    type(qs_pao_matrix_type), intent(in) :: source_m
    type(qs_pao_matrix_type), intent(inout) :: target_m
    type(qs_pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types:qs_pao_matrix_transf_to_m'
    type(qs_matrix_block_iterator) :: iterator
    integer :: brow,bcol, stat, max_bsize,i,j
    integer, dimension(:), pointer :: l_index,r_index
    real(kind=wp), dimension(:,:), pointer :: block_val, u_col, u_row
    real(kind=wp),allocatable, dimension(:,:) :: t1_m,t2_m
    type(qs_pao_local_angles_type) :: l_angles
    type(qs_error_type) :: iError
    type(qs_pao_env_type), pointer :: pao_env
    call qs_init(iError,template_error=error)
    nullify(block_val)
    failure=.false.

    QSPrecondition(source_m%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(target_m%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(glob_angles%initialized,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       QSPrecondition(.not.source_m%min_bas_cols,qs_failure_level,routineN,error,failure)
       QSPrecondition(.not.source_m%min_bas_rows,qs_failure_level,routineN,error,failure)
       QSPrecondition(target_m%min_bas_cols,qs_failure_level,routineN,error,failure)
       QSPrecondition(target_m%min_bas_rows,qs_failure_level,routineN,error,failure)
    end if
    if (.not. failure) then
       pao_env => qs_get_pao_env(source_m)
       QSPrecondition(associated(pao_env),qs_failure_level,routineN,error,failure)
    end if
    if (.not.failure) then
       QSPrecondition(associated(pao_env%nr_full_bas),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(pao_env%projection),qs_failure_level,routineN,error,failure)
    end if
    if (.not.failure) then
       call qs_init(l_angles,atom_nr=1, &
            glob_angles=glob_angles, error=error)
       max_bsize=maxval(pao_env%nr_full_bas)
       allocate(t1_m(max_bsize,max_bsize), t2_m(max_bsize,max_bsize), stat=stat)
       QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
       call qs_init(iterator,source_m%matrix,error=iError)
       call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
            error,failure)
       if (.not.failure) then
          blocksDo: do
             if ((.not.qs_next(iterator,error=error)).or.failure) exit
             call qs_get(iterator,block_row=brow,block_col=bcol,&
                  block_val=block_val)
! get u
             call qs_set(l_angles, atom_nr=bcol,error=error)
             call qs_get(l_angles,unitary_t=u_col)
!call get_block_node(unitary_t%matrix, block_col=bcol,&
!     block_row=bcol, block=u_col)
             call qs_set(l_angles, atom_nr=brow,error=error)
             call qs_get(l_angles,unitary_t=u_row)
!call get_block_node(unitary_t%matrix, block_col=brow,&
!     block_row=brow, block=u_row)

! in the following the sizes are put to make sense, i.e. not using
! size(u_row,1)==size(u_row,2)==size(source_m,1)
! size(u_col,1)==size(u_col,2)==size(source_m,2)

! t1_m = source u
             call dgemm('N','N',size(block_val,1),size(block_val,2),&
                  size(u_col,2),1.0_wp,&
                  block_val,size(block_val,1),&
                  u_col,size(u_col,1),&
                  0.0_wp,t1_m,size(t1_m,1))
! t2_m = u^* t1_m
             call dgemm('T','N',size(u_row,1),size(u_row,2),&
                  size(u_col,2),1.0_wp,&
                  u_row,size(u_row,1),&
                  t1_m,size(t1_m,1),&
                  0.0_wp,t2_m,size(t2_m,1))
! get small basis block
             call get_block_node(target_m%matrix, block_col=bcol,&
                  block_row=bcol, block=block_val)
             if (.not.associated(block_val)) then
                call add_block_node(target_m%matrix, block_col=bcol,&
                     block_row=bcol, block=block_val)
                QSPostcondition(associated(block_val),qs_failure_level,routineN,error,failure)
             end if

             l_index => pao_env%projection(bcol)%array
             r_index => pao_env%projection(brow)%array
             QSPrecondition(associated(l_index),qs_failure_level,routineN,error,failure)
             QSPrecondition(associated(r_index),qs_failure_level,routineN,error,failure)
             if (.not.failure) then
                call qs_assert(size(l_index)==size(block_val,1),&
                     qs_failure_level,qs_precondition_failed,routineN,&
                     "size mishmash in "//QSSourceFileRef(),error,failure)
                call qs_assert(size(r_index)==size(block_val,2),&
                     qs_failure_level,qs_precondition_failed,routineN,&
                     "size mishmash in "//QSSourceFileRef(),error,failure)
             end if
             if (.not.failure) then ! put result in the small basis block
                do i=1,size(block_val,1)
                   do j=1,size(block_val,2)
                      block_val(i,j)=t2_m(l_index(i),r_index(j))
                   end do
                end do
             end if

          end do blocksDo
       end if
       call qs_dealloc_ref(iterator,error=error)
       call qs_dealloc_ref(l_angles,error=error)
       deallocate(t1_m,t2_m,stat=stat) ! just to be sure...
       QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_pao_matrix_transf_to_m
!!***
!***************************************************************************

!end
end module qs_pao_types

