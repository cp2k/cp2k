!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_qs_env_methods [1.0] *
!!
!!   NAME
!!     pao_qs_env_methods
!!
!!   FUNCTION
!!     pao related methods that work on the qs_env methods
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
module pao_qs_env_methods
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level
  use cp_error_handling, only: cp_debug, cp_error_type, cp_assert, &
       cp_error_message, cp_error_get_logger
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types, only: pao_env_type, pao_env_get, pao_glob_angles_type
  use pao_bas_projection_types, only: pao_bas_projection_type
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind,&
       get_atomic_kind_set
  use pao_glob_angles_methods, only: pao_g_ang_did_change
  use cp_block_matrix, only: cp_block_matrix_p_type, cp_block_matrix_type,&
       cp_sp_create, cp_sp_release, cp_get_matrix, &
       cp_sp_block_handle_type, cp_sp_block_handle_init, cp_sp_block_handle_release
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type, cp_bmstruct_get,&
       cp_bmstruct_create, cp_bmstruct_release, cp_b_dims_type
  use cp_matrix_utils, only: cp_sm_set_dist_sparsity_from_sm
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type,&
       get_block_node,add_block_node
  use mathlib, only: symmetrize_matrix, power_matrix
  use cp_array_utils, only: cp_2d_r_output
  use particle_types, only: particle_type
  use dft_types, only: dft_control_type
  use cp_para_env, only: cp_para_env_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_qs_env_methods'

  public :: cp_get_b_projection, cp_did_change
  ! underlying routines
  public :: pao_qs_env_get_b_projection, pao_qs_env_get, &
       pao_qs_env_did_change, pao_qs_env_get_atomic_ortho
  public :: pao_bas_projection_type ! for the intel compiler
  interface cp_get_b_projection
     module procedure pao_qs_env_get_b_projection
  end interface
  interface cp_did_change
     module procedure pao_qs_env_did_change
  end interface
!***
!****************************************************************************
contains

!!****f* pao_qs_env_methods/pao_qs_env_get_b_projection [1.0] *
!!
!!   NAME
!!     pao_qs_env_get_b_projection
!!
!!   SYNOPSIS
!!     Function pao_qs_env_get_b_projection(qs_env, atom_nr, error)&
!!         Result(res)
!!       Type(pao_bas_projection_type), Pointer:: res
!!       Type(qs_environment_type), Intent (IN):: qs_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_qs_env_get_b_projection
!!
!!   FUNCTION
!!     returns the base projection projection on atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     qs_env: the qs environement you get the data from
!!     atom_nr: the number of the atom
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 1.2002 created [fawzi]
!!     - 19.04.02 used new projection
!!     - 23.04.02 moved to qs
!!
!!*** **********************************************************************
  function pao_qs_env_get_b_projection(qs_env, atom_nr, error) result(res)
    type(pao_bas_projection_type), pointer ::res
    type(qs_environment_type), intent(in) :: qs_env
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_qs_env_get_b_projection',&
        routineP=moduleN//':'//routineN
    type(particle_type), dimension(:), pointer :: particle_set
    failure=.false.

    call get_qs_env(qs_env, particle_set=particle_set)
    call get_atomic_kind(particle_set(atom_nr)%atomic_kind,pao_bas_proj=res)
  end function pao_qs_env_get_b_projection
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_get [1.0] *
!!
!!   NAME
!!     pao_qs_env_get
!!
!!   SYNOPSIS
!!     Subroutine pao_qs_env_get(qs_env, global_env, pao_env,&
!!         pao_max_min_bas, pao_max_full_bas, pao_max_excl_bas,&
!!         pao_max_l_angles, pao_tot_min_bas, pao_tot_full_bas,&
!!         pao_tot_excl_bas, pao_min_density_m, pao_min_bas_m_struct,&
!!         pao_full_bas_m_struct, pao_min_full_bas_m_struct,&
!!         pao_full_min_bas_m_struct, pao_full_bas_diag_m_struct,&
!!         pao_angles_att, pao_w_subspace_pos_precision,&
!!         pao_w_in_subspace_precision, pao_atomic_ortho,&
!!         pao_atomic_ortho_valid, n_atoms, n_electrons, error)
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(pao_env_type), Pointer, Optional:: pao_env
!!       Integer, Intent (OUT), Optional:: pao_max_min_bas,&
!!         pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
!!         pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas, n_atoms,&
!!         n_electrons
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         pao_w_subspace_pos_precision
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         pao_w_in_subspace_precision
!!       Type(pao_glob_angles_type), Optional, Pointer:: pao_angles_att
!!       Type(cp_block_matrix_p_type), Dimension(:), Optional, Pointer::&
!!         pao_min_density_m
!!       Type(cp_b_matrix_struct_type), Pointer, Optional::&
!!         pao_min_bas_m_struct, pao_full_bas_m_struct,&
!!         pao_min_full_bas_m_struct, pao_full_min_bas_m_struct,&
!!         pao_full_bas_diag_m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical, Dimension(:), Pointer, Optional:: pao_atomic_ortho_valid
!!       Type(cp_block_matrix_type), Pointer, Optional:: pao_atomic_ortho
!!     End Subroutine pao_qs_env_get
!!
!!   FUNCTION
!!     returns pao related info from the qs_environment_type
!!
!!   NOTES
!!     use this function rather than get_atomic_kind_set: some values
!!     might be cached in pao_env in the future
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine pao_qs_env_get(qs_env, global_env, pao_env, pao_max_min_bas,&
     pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
     pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas,&
     pao_min_density_m, pao_min_bas_m_struct, pao_full_bas_m_struct,&
     pao_min_full_bas_m_struct, pao_full_min_bas_m_struct,&
     pao_full_bas_diag_m_struct,&
     pao_angles_att, pao_w_subspace_pos_precision,&
     pao_w_in_subspace_precision,pao_atomic_ortho,&
     pao_atomic_ortho_valid, n_atoms, n_electrons, error)
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  type(pao_env_type), pointer, optional :: pao_env
  integer, intent(out), optional :: pao_max_min_bas,&
       pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
       pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas,&
       n_atoms, n_electrons
  real(kind=wp), optional, intent(out) :: pao_w_subspace_pos_precision
  real(kind=wp), optional, intent(out) :: pao_w_in_subspace_precision
  type(pao_glob_angles_type), optional, pointer :: pao_angles_att
  type(cp_block_matrix_p_type), dimension(:), optional, pointer ::&
       pao_min_density_m
  type(cp_b_matrix_struct_type), pointer, optional :: pao_min_bas_m_struct, &
       pao_full_bas_m_struct,pao_min_full_bas_m_struct, &
       pao_full_min_bas_m_struct, pao_full_bas_diag_m_struct
  type(cp_error_type), optional, intent(inout) :: error
  logical, dimension(:), pointer, optional :: pao_atomic_ortho_valid
  type(cp_block_matrix_type), pointer, optional :: pao_atomic_ortho
  
  logical :: failure
  integer :: natoms
  character(len=*), parameter :: routineN='pao_qs_env_get',&
        routineP=moduleN//':'//routineN
  type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
  type(pao_env_type), pointer :: p_env
  type(dft_control_type), pointer :: dft_control
  natoms = size(qs_env%particle_set)
  failure=.false.
  
  call get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pao_env=p_env)
  CPPrecondition(associated(p_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     if (present(pao_max_min_bas)) pao_max_min_bas=p_env%max_min_bas
     if (present(pao_max_full_bas)) pao_max_full_bas=p_env%max_full_bas
     if (present(pao_max_excl_bas)) pao_max_excl_bas=p_env%max_excl_bas
     if (present(pao_tot_min_bas)) pao_tot_min_bas=p_env%tot_min_bas
     if (present(pao_tot_full_bas)) pao_tot_full_bas=p_env%tot_full_bas
     if (present(pao_tot_excl_bas)) pao_tot_excl_bas=p_env%tot_excl_bas
     if (present(pao_max_l_angles)) pao_max_l_angles=p_env%max_l_angles
     call pao_env_get(p_env,&
          wanted_subspace_pos_precision=pao_w_subspace_pos_precision,&
          wanted_in_subspace_precision=pao_w_in_subspace_precision,&
          angles_att=pao_angles_att,min_density_m=pao_min_density_m,&
          min_bas_m_struct=pao_min_bas_m_struct,&
          min_full_bas_m_struct=pao_min_full_bas_m_struct,&
          full_min_bas_m_struct=pao_full_min_bas_m_struct,&
          full_bas_diag_m_struct=pao_full_bas_diag_m_struct,&
          full_bas_m_struct=pao_full_bas_m_struct, error=error)
     if (present(pao_atomic_ortho_valid)) then
        if (.not.associated(pao_env%atomic_ortho_valid)) then
           allocate(p_env%atomic_ortho_valid(natoms))
           p_env%atomic_ortho_valid=.false.
        end if
        pao_atomic_ortho_valid => p_env%atomic_ortho_valid
     end if
     if (present(pao_atomic_ortho)) then
        if (.not.associated(p_env%atomic_ortho_valid)) then
           allocate(p_env%atomic_ortho_valid(natoms))
           p_env%atomic_ortho_valid=.false.
        end if
        if (.not.associated(p_env%atomic_ortho)) then
           call cp_sp_create(p_env%atomic_ortho,&
                matrix_structure=p_env%full_min_bas_m_struct,&
                error=error)
        end if
        pao_atomic_ortho => p_env%atomic_ortho
     end if
     if (present(n_atoms)) n_atoms = natoms
     if (present(n_electrons)) then
        call get_qs_env(qs_env, atomic_kind_set=atomic_kind_set,&
             dft_control=dft_control)
        CPPrecondition(associated(atomic_kind_set),cp_failure_level,routineP,error,failure)
        CPPrecondition(associated(dft_control),cp_failure_level,routineP,error,failure)
        if (.not.failure) then
           call get_atomic_kind_set(atomic_kind_set,nelectron=n_electrons)
           n_electrons=n_electrons-dft_control%charge
        else
           n_electrons=-1
        end if
     end if
  end if
  if (present(pao_env)) pao_env => p_env
end subroutine pao_qs_env_get
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_did_change [1.0] *
!!
!!   NAME
!!     pao_qs_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_qs_env_did_change(qs_env, global_env, h_changed,&
!!         s_changed, angles_changed, full_reset, error)
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_qs_env_did_change
!!
!!   FUNCTION
!!     invalidates the cache when some data changes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environement that has changed
!!     - s_changed: if the overlap has changed (implies h_changed
!!     - angles_changed: if the angles that describe tha actual subspace 
!!       have changed
!!     - full_reset: if the structures should be regenerated from scratch
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
subroutine pao_qs_env_did_change(qs_env,global_env,&
     h_changed, s_changed,&
     angles_changed, full_reset, error)
  type (qs_environment_type), intent(inout) :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  logical, intent(in), optional :: h_changed, s_changed, &
       angles_changed, full_reset
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_qs_env_did_change',&
       routineP=moduleN//':'//routineN
  logical :: s_c, full_r
  integer :: i
  type(pao_env_type), pointer :: pao_env
  type(real_matrix_p_type), dimension(:), pointer :: overlap_rmpv
  type(cp_b_dims_type), pointer :: b_dims
  type(cp_b_matrix_struct_type), pointer :: old_m_struct
  type(cp_para_env_type), pointer :: para_env
  failure=.false.
  s_c=.false.; full_r=.false.
  
  call get_qs_env(qs_env,pao_env=pao_env)
  CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
  end if
  if (.not. failure) then
     if (present(s_changed)) s_c=s_changed
     if (present(full_reset)) full_r=full_reset
     if (s_c.or.full_r) then
        call get_qs_env(qs_env, s=overlap_rmpv)
        ! invalidate min_density_m
        pao_env%min_d_m_valid=.false.
        do i=1,size(pao_env%min_density_m)
           call cp_sp_release(pao_env%min_density_m(i)%matrix,error=error)
        end do
        ! invalidate matrix structure
        call cp_bmstruct_get(pao_env%min_bas_m_struct,&
            b_dims=b_dims, error=error)
        old_m_struct => pao_env%min_bas_m_struct
        call cp_bmstruct_create(pao_env%min_bas_m_struct,&
               symmetric=.true.,&
               para_env=global_env%para_env,&
               b_dims=b_dims,&
               error=error)
        call cp_sm_set_dist_sparsity_from_sm(pao_env%min_bas_m_struct,&
             para_env=global_env%para_env,&
             sample_matrix=overlap_rmpv(1)%matrix,&
             error=error)
        call cp_bmstruct_release(old_m_struct,error=error)

        call cp_bmstruct_get(pao_env%full_bas_m_struct,&
             para_env=para_env, b_dims=b_dims, error=error)
        old_m_struct => pao_env%full_bas_m_struct
        call cp_bmstruct_create(pao_env%full_bas_m_struct,&
             symmetric=.false.,&
             para_env=para_env,&
             b_dims=b_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        call cp_bmstruct_release(old_m_struct,error=error)

        call cp_bmstruct_get(pao_env%full_min_bas_m_struct,&
             para_env=para_env, b_dims=b_dims, error=error)
        old_m_struct => pao_env%full_min_bas_m_struct
        call cp_bmstruct_create(pao_env%full_min_bas_m_struct,&
             symmetric=.false.,&
             para_env=para_env,&
             b_dims=b_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        call cp_bmstruct_release(old_m_struct,error=error)

        call cp_bmstruct_get(pao_env%min_full_bas_m_struct,&
             para_env=para_env, b_dims=b_dims, error=error)
        old_m_struct => pao_env%min_full_bas_m_struct
        call cp_bmstruct_create(pao_env%min_full_bas_m_struct,&
             symmetric=.false.,&
             para_env=para_env,&
             b_dims=b_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        call cp_bmstruct_release(old_m_struct,error=error)
        
        do i=1,size(pao_env%min_density_m)
           call cp_sp_create(pao_env%min_density_m(i)%matrix,&
                matrix_structure= pao_env%min_bas_m_struct,&
                error=error)
        end do
        ! invalidate atomic orthogonalisation
        if (associated(pao_env%atomic_ortho_valid)) &
             pao_env%atomic_ortho_valid=.false.
        ! in general the atom distribution changes
        call cp_sp_release(pao_env%atomic_ortho,error=error)
        call pao_qs_env_preinit_atom_ortho(qs_env, global_env, error=error)
     end if
     call pao_g_ang_did_change(pao_env%angles_att,qs_env=qs_env,&
          global_env=global_env, h_changed=h_changed,&
          s_changed=s_changed, angles_changed=angles_changed,&
          full_reset=full_reset, error=error)
  end if
end subroutine pao_qs_env_did_change
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_get_atomic_ortho [1.0] *
!!
!!   NAME
!!     pao_qs_env_get_atomic_ortho
!!
!!   SYNOPSIS
!!     Function pao_qs_env_get_atomic_ortho(qs_env, global_env, atom_nr,&
!!         error) Result(res)
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(Kind=wp), Dimension(:,:), Pointer:: res
!!     End Function pao_qs_env_get_atomic_ortho
!!
!!   FUNCTION
!!     returns the atomic orthogonalization matrix, calculating it
!!     if necessary.
!!
!!   NOTES
!!     could cache tmp_m (but submatrix indexing on most compilers is
!!     not optimized)
!!
!!   INPUTS
!!     - qs_env: the qs_env that contains the overlap matrix
!!     - global_env: the global environment (for the communication)
!!     - atom_nr: the atom you want the atomic orthogonalization from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     6.2002 created [fawzi]
!!
!!*** **********************************************************************
function pao_qs_env_get_atomic_ortho(qs_env,global_env,atom_nr,error)&
     result(res)
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  integer, intent(in) :: atom_nr
  type(cp_error_type), optional, intent(inout) :: error
  real(kind=wp), dimension(:,:), pointer :: res

  logical :: failure
  character(len=*), parameter :: routineN='pao_qs_env_get_atomic_ortho',&
       routineP=moduleN//':'//routineN
  integer :: stat, n_atoms
  integer, save :: iter=0
  type(cp_block_matrix_type), pointer :: atomic_ortho
  type(pao_env_type), pointer :: pao_env
  type(real_matrix_p_type), dimension(:), pointer :: s
  real(kind=wp), dimension(:,:), pointer :: block_val
  real(kind=wp), dimension(:,:), allocatable :: tmp_m
  real(kind=wp) :: eps_eigval
  failure=.false.
  nullify(res, atomic_ortho, block_val, pao_env, s)

  call pao_qs_env_get(qs_env,global_env,pao_env=pao_env,n_atoms=n_atoms)
  CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     if (.not.associated(pao_env%atomic_ortho_valid)) then
        allocate(pao_env%atomic_ortho_valid(n_atoms),stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        if (.not.failure) then
           pao_env%atomic_ortho_valid=.false.
        end if
     end if
     if (.not.associated(pao_env%atomic_ortho)) then
        call cp_sp_create(pao_env%atomic_ortho,&
             matrix_structure=pao_env%full_bas_diag_m_struct,&
             error=error)
     end if
     atomic_ortho => pao_env%atomic_ortho
  end if
  CPPrecondition(associated(atomic_ortho),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(pao_env%atomic_ortho_valid),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     if (pao_env%atomic_ortho_valid(atom_nr)) then
        call get_block_node(cp_get_matrix(atomic_ortho),&
             block=res,block_row=atom_nr,&
             block_col=atom_nr)
        CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
     else
        call get_qs_env(qs_env,s=s)
        call get_block_node(cp_get_matrix(atomic_ortho),&
             block=res,block_row=atom_nr,&
             block_col=atom_nr)
        if (.not.associated(res)) then
           call add_block_node(cp_get_matrix(atomic_ortho),&
                block_row=atom_nr,block_col=atom_nr,block=res)
        end if
        CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
        if (.not.failure) then
           CPPrecondition(size(res,1)==size(res,2),cp_failure_level,routineP,error,failure)
        end if
        if (.not.failure) then
           allocate(tmp_m(size(res,1),size(res,2)),stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        end if
        if (.not.failure) then
           call get_block_node(s(1)%matrix,block_row=atom_nr,&
                block_col=atom_nr,block=block_val)
           CPPrecondition(size(block_val,1)==size(res,1),cp_failure_level,routineP,error,failure)
           CPPrecondition(size(block_val,2)==size(res,2),cp_failure_level,routineP,error,failure)
           if (.not.failure) then
              tmp_m=block_val
              eps_eigval= 10 * epsilon(0.0_wp) ! to do: use the scf value
              call power_matrix(tmp_m,res,-0.5_wp,threshold=eps_eigval)
              call symmetrize_matrix(res,"upper_to_lower")
              if (cp_debug.and.debug_this_module) then
                 iter=iter+1
                 call cp_2d_r_output(cp_error_get_logger(error),&
                      outputName="atomic_ortho",&
                      fromWhere=routineP, iter=iter,array=res,&
                      comment="atom_nr="//cp_to_string(atom_nr),&
                      error=error)
              end if
           end if
           deallocate(tmp_m,stat=stat) ! just to be sure...
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
           pao_env%atomic_ortho_valid(atom_nr)=.true.
        end if
     end if
  end if
  if (failure) nullify(res)
end function pao_qs_env_get_atomic_ortho
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_preinit_atom_ortho [1.0] *
!!
!!   NAME
!!     pao_qs_env_preinit_atom_ortho
!!
!!   FUNCTION
!!     If pao is run in parallel shares the ortho
!!
!!   NOTES
!!     could cache tmp_m (but submatrix indexing on most compilers is
!!     not optimized)
!!
!!   INPUTS
!!     - qs_env: the qs_env that contains the overlap matrix
!!     - global_env: the global environment (for the communication)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     6.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine pao_qs_env_preinit_atom_ortho(qs_env,global_env,error)
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='pao_qs_env_preinit_atom_ortho',&
       routineP=moduleN//':'//routineN
  integer :: stat, n_atoms,i
  type(cp_block_matrix_type), pointer :: atomic_ortho
  type(pao_env_type), pointer :: pao_env
  type(cp_sp_block_handle_type) :: block_handle
  
  failure=.false.
  nullify(atomic_ortho, pao_env)

  call pao_qs_env_get(qs_env,global_env,pao_env=pao_env,n_atoms=n_atoms)
  CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     if (.not.associated(pao_env%atomic_ortho_valid)) then
        allocate(pao_env%atomic_ortho_valid(n_atoms),stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        if (.not.failure) then
           pao_env%atomic_ortho_valid=.false.
        end if
     end if
     if (.not.associated(pao_env%atomic_ortho)) then
        call cp_sp_create(pao_env%atomic_ortho,&
             matrix_structure=pao_env%full_bas_diag_m_struct,&
             error=error)
     end if
     atomic_ortho => pao_env%atomic_ortho
  end if
  CPPrecondition(associated(atomic_ortho),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(pao_env%atomic_ortho_valid),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
#ifdef __parallel
!!FM     ! stupid parallelization, every processor has everything
!!FM     do i=1,n_atoms
!!FM        if (atomic_ortho%matrix_structure%distribution%block_owner(i,i) ==&
!!FM             atomic_ortho%global_env%mepos) then
!!FM           call pao_qs_env_get_atomic_ortho(qs_env,&
!!FM                global_env=atomic_ortho%global_env,&
!!FM                atom_nr=i,error=error)
!!FM        end if
!!FM     end do
!!FM     call cp_sp_block_handle_init(block_handle,error=error)
!!FM     do i=1,nblock_rows
!!FM        call cp_sp_get_local_block(matrix=atomic_ortho,&
!!FM             block_row=i, block_col=i,&
!!FM             block_handle=block_handle,alloc_block=.true.,error=error,&
!!FM             owner=owner)
!!FM        call mp_bcast(block_handle%block_val,owner,atomic_ortho%global_env%group)
!!FM        call cp_sp_block_handle_release(block_handle,error=error)
!!FM     end do
#endif
  end if
end subroutine pao_qs_env_preinit_atom_ortho
!***************************************************************************

end module pao_qs_env_methods
