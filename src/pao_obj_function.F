!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_obj_function [1.0] *
!!
!!   NAME
!!     pao_obj_function
!!
!!   FUNCTION
!!     contain the objective function to be optimized from the pao procedure
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module pao_obj_function
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_log, cp_logger_type, cp_note_level,cp_get_default_logger
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_error_dealloc_ref, cp_error_get_logger
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use pao_types, only: pao_env_type, pao_glob_angles_type, cp_get
!  use pao_env_methods
  use cp_block_matrix, only: cp_block_matrix_type, cp_get_matrix,&
       cp_sp_create, cp_sp_release, cp_block_matrix_p_type,&
       cp_sp_alloc_nonsparse_blocks
  use pao_glob_angles_methods, only: pao_g_ang_set, pao_g_ang_did_change
  use pao_l_angles_methods, only: pao_local_angles_type, cp_get_angle,&
       pao_l_ang_init, pao_l_ang_get, cp_get_g_index, cp_set, &
       pao_l_ang_dealloc_ref
  use pao_proj_methods, only: cp_m2f_with_B_M_pi,pao_g_ang_get_min_m
  use qs_blacs, only: cp_full_matrix_type,cp_fm_create2,&
       copy_sparse_to_blacs_matrix,cp_fm_release,&
       cp_fm_trace, add_blacs_to_block_diag_sm,&
       copy_blacs_to_blacs_matrix, cp_fm_gemm, copy_blacs_to_sparse_matrix,&
       cp_fm_syrk,blacs_syevx,cp_fm_symm,cp_fm_get_info,&
       symmetrise_blacs_matrix, power_blacs_matrix
  use sparse_matrix_types, only: symmetrise_diagonal_blocks
  use timings, only: timeset, timestop
  USE qs_mo_types, ONLY: calculate_density_matrix,&
                         get_mo_set,allocate_mo_set, deallocate_mo_set,&
                         mo_set_type, mo_set_p_type
  use message_passing, only: mp_sum
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use pao_unitary_t, only: pao_calc_u_du_type, cp_init, cp_next, &
       cp_dealloc_ref, pao_calc_u_du_get
  use pao_qs_env_methods, only: pao_qs_env_get,pao_qs_env_get_atomic_ortho
  use cp_matrix_utils, only: cp_matrix_block_iterator, cp_sm_b_i_init,&
       cp_next, cp_next, cp_sm_b_i_dealloc_ref,cp_d_struct_create,&
       cp_init_diagonal_blocks, cp_get, cp_init, cp_set_to, cp_sm_output
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type,&
       get_matrix_info
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind_set
  use cp_array_utils, only: cp_1d_r_output, cp_2d_r_output
  use cp_output_handling, only: cp_write_output
  use cp_para_env, only: cp_para_env_type
  implicit none
  private
  
  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN="pao_obj_function"
  integer, parameter :: cp_optimal_blacs_block_size=32
  
  ! type
  public :: pao_gerd_functional_type

  ! generic functions
  public :: cp_init, cp_dealloc_ref, cp_eval_at, cp_valid

  ! underlying routines
  public :: pao_gerd_f_eval_at, pao_gerd_f_valid,&
       pao_gerd_f_init, pao_gerd_f_dealloc_ref
  
  interface cp_init
     module procedure pao_gerd_f_init
  end interface
  interface cp_dealloc_ref
     module procedure pao_gerd_f_dealloc_ref
  end interface
  interface cp_valid
     module procedure pao_gerd_f_valid
  end interface
  interface cp_eval_at
     module procedure pao_gerd_f_eval_at
  end interface

!!***
!****************************************************************************

!!****s* pao_obj_function/pao_gerd_functional_type [1.0] *
!!
!!   NAME
!!     pao_gerd_functional_type
!!
!!   FUNCTION
!!     the environement for the gert functional, basically all the
!!     temporary space, this way it is not allocated everytime,
!!     and it is not a save attribute (so it is multithread safe;)
!!
!!   NOTES
!!     the functional is:
!!       f = Tr[P (N U i)^* H (N U i)] + Tr[V(P - P (N U i)^* S (N U i) P)]
!!     where P is the density matrix in the minimal basis, 
!!     N an atomic ortogonalization term ( S_{at}^{-1/2} )
!!     in the full basis, U a rotation between the full basis and the
!!     min basis, i an injection into the min basis,
!!     H the hamiltonian in the full basis, V langrange multiplicators
!!     for the idempotency of P, S is the the overlap matrix in the 
!!     full basis, and Tr[A] is the trace of a with respect to the
!!     standard scalar product (Tr[A]=sum(A_{i,i})).
!!     This functional was proposed by Gerd Breghold in his dissertation
!!     at the Max-Planck-Institut, and is derived from the functional used 
!!     by Lee and Head-Gordon in Chem. Phys. 107,p 9085, (1997)
!!
!!   ATTRIBUTES
!!     - initialized: true if initialized (bug catcher)
!!     - pao_env: a pointer to the actual pao_env. It is assumed that the
!!       size of the min and full basis, and the pao env do not change
!!       between sucessive evaluations
!!     - min_hamiltonian_m: qs_blacs matrix with the hamiltonian in the 
!!       min basis
!!     - min_overlap_m: qs_blacs matrix with the overlap in the min basis
!!     - density_m: qs_blacs matrix with the density matrix
!!     - work_m: qs_blacs matrix in the small basis
!!     - alpha_mo: molecular orbitals
!!     - scf_work1, scf_work2, scf_work3: qs_blacs matrices in the big basis
!!     - diag_matrix: diagonal sparse matrix
!!     - full_work_m: sparse matrix in the full basis
!!     - min_work_m: sparse matrix in the minimal basis
!!     - d_struct: the distribution structure of the diagonal blocks between
!!       the processors (d_struct(i) gives on which processor is blocks i)
!!     - n_eval: the number of function evaluations
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
  type pao_gerd_functional_type
     logical :: initialized, should_dealloc_mo, rebuild_rho
     integer :: n_eval
     type(qs_environment_type), pointer :: qs_env
     type(cp_full_matrix_type), POINTER :: min_blacs_1, min_blacs_2,&
          min_blacs_3
     type(cp_full_matrix_type), pointer :: scf_work1, scf_work2, scf_work3
     type(mo_set_p_type), dimension(:), pointer :: mo_set
     type(cp_block_matrix_type), pointer :: diag_matrix, full_work_m,&
          min_work_m
     integer, dimension(:), pointer :: d_struct
     type(global_environment_type), pointer :: global_env
  end type pao_gerd_functional_type
!!***
!****************************************************************************
contains

!!****f* pao_obj_function/pao_gerd_f_init [1.0] *
!!
!!   NAME
!!     pao_gerd_f_init
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_init(gf_env, qs_env, global_env, scf_work1,&
!!         scf_work2, scf_work3, mo_orbital, should_dealloc_mo, error)
!!       Type(pao_gerd_functional_type), Intent (OUT):: gf_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_full_matrix_type), Intent (IN), Target:: scf_work1,&
!!         scf_work2, scf_work3
!!       Type(mo_set_p_type), Dimension(:), Target, Intent (IN), Optional::&
!!         mo_orbital
!!       Logical, Intent (IN), Optional:: should_dealloc_mo
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_init
!!
!!   FUNCTION
!!     initializes gerd s functional
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gf_env: the object to be initialized
!!     pao_env: the pao environement of this functional
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling.
!!
!!     see pao_gerd_functional_type for the other variables
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_gerd_f_init(gf_env,qs_env,global_env, scf_work1,scf_work2,&
       scf_work3, mo_orbital, should_dealloc_mo, error)
    type(pao_gerd_functional_type), intent(out) :: gf_env
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in),target :: global_env
    type(cp_full_matrix_type), intent(in), target :: scf_work1, scf_work2,&
         scf_work3
    type(mo_set_p_type), dimension(:), target, intent(in),optional ::&
         mo_orbital
    logical, intent(in), optional :: should_dealloc_mo
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_obj_function:pao_gerd_f_init'
    integer :: ntot_min_bas, ntot_full_bas, i, stat, nelectron, n_atoms
    type(real_matrix_p_type), dimension(:), pointer:: full_overlap_m
    type(real_matrix_type), pointer :: diagMatrix
    type(cp_b_matrix_struct_type), pointer :: min_bas_m_struct, full_bas_m_struct
    type(pao_env_type), pointer :: pao_env
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    nullify(diagMatrix, min_bas_m_struct, full_bas_m_struct, pao_env,&
         atomic_kind_set)
    failure=.false.

    gf_env%qs_env => qs_env
    gf_env%global_env => global_env
    gf_env%n_eval = 0
    if (.not.failure) then
       call get_qs_env(gf_env%qs_env, s=full_overlap_m, pao_env=pao_env,&
            atomic_kind_set=atomic_kind_set)
       call pao_qs_env_get(gf_env%qs_env, global_env=global_env,&
            pao_tot_full_bas=ntot_full_bas,&
            pao_tot_min_bas=ntot_min_bas,&
            pao_min_bas_m_struct=min_bas_m_struct,&
            pao_full_bas_diag_m_struct=full_bas_m_struct, n_atoms=n_atoms,&
            error=error)
       call cp_fm_create2(gf_env%min_blacs_1,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min blacs 1',&
            globenv=global_env)
       call cp_fm_create2(gf_env%min_blacs_2,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min blacs 2',&
            globenv=global_env)
       call cp_fm_create2(gf_env%min_blacs_3,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min blacs 3',&
            globenv=global_env)
       gf_env%scf_work1 => scf_work1
       gf_env%scf_work2 => scf_work2
       gf_env%scf_work3 => scf_work3
       if (present(mo_orbital)) then
          allocate(gf_env%mo_set(size(mo_orbital)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          do i=1,size(mo_orbital)
             gf_env%mo_set(i)%mo_set => mo_orbital(i)%mo_set
          end do
          gf_env%should_dealloc_mo=.false.
       else
          call get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
               nelectron=nelectron)
          allocate(gf_env%mo_set(1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          call allocate_mo_set(gf_env%mo_set(1)%mo_set,nao=ntot_min_bas,&
               nmo=nelectron/2,&
               nelectron=nelectron,maxocc=2.0_wp,nrow_block=n_atoms,&
               ncol_block=n_atoms,&
               globenv=global_env)
          gf_env%should_dealloc_mo=.true.
       end if
       if (present(should_dealloc_mo)) &
            gf_env%should_dealloc_mo=should_dealloc_mo
       call cp_sp_create(gf_env%diag_matrix,&
            matrix_struct=full_bas_m_struct)
       call cp_d_struct_create(gf_env%d_struct,&
            matrix=full_overlap_m(1)%matrix,&
            para_env=global_env%para_env,error=error)
       diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
       call cp_init_diagonal_blocks(diagMatrix,&
            para_env=global_env%para_env,&
            default_value=0.0_wp,&
            d_struct=gf_env%d_struct)
       call cp_sp_create(gf_env%min_work_m,&
            matrix_struct=min_bas_m_struct)
       call cp_sp_alloc_nonsparse_blocks(gf_env%min_work_m,error=error)
       call cp_sp_create(gf_env%full_work_m,&
            matrix_struct=full_bas_m_struct)
       gf_env%rebuild_rho=.true.
       gf_env%initialized=.true.
    end if
  end subroutine pao_gerd_f_init
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_gerd_f_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_dealloc_ref(gf_env, error)
!!       Type(pao_gerd_functional_type), Optional, Intent (INOUT):: gf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory that has been allocated by the functional
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gf_env: the functional to be released
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
  subroutine pao_gerd_f_dealloc_ref(gf_env, error)
    type(pao_gerd_functional_type), optional, intent(inout) :: gf_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_obj_function:pao_gerd_f_dealloc_ref'
    integer :: stat,i
    failure=.false.

    CPPrecondition(gf_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_fm_release(gf_env%min_blacs_1)
       call cp_fm_release(gf_env%min_blacs_2)
       call cp_fm_release(gf_env%min_blacs_3)
       nullify(gf_env%scf_work1, gf_env%scf_work2, gf_env%scf_work3)
       call cp_sp_release(gf_env%min_work_m, error=error)
       call cp_sp_release(gf_env%full_work_m, error=error)
       call cp_sp_release(gf_env%diag_matrix, error=error)
       deallocate(gf_env%d_struct, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       if (associated(gf_env%mo_set)) then
          if (gf_env%should_dealloc_mo) then
             do i=1,size(gf_env%mo_set)
                call deallocate_mo_set(gf_env%mo_set(i)%mo_set)
             end do
          end if
          deallocate(gf_env%mo_set,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       nullify(gf_env%global_env)
       gf_env%initialized=.false.
    end if
  end subroutine pao_gerd_f_dealloc_ref
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_valid [1.0] *
!!
!!   NAME
!!     pao_gerd_f_valid
!!
!!   SYNOPSIS
!!     Function pao_gerd_f_valid(gf_env, error) Result(res)
!!       Logical:: res
!!       Type(pao_gerd_functional_type), Intent (IN):: gf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_gerd_f_valid
!!
!!   FUNCTION
!!     returns true if the functional object is valid (minimal testing)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gf_env: the functional to test
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
  function pao_gerd_f_valid(gf_env, error) result(res)
    logical ::res
    type(pao_gerd_functional_type), intent(in) :: gf_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_obj_function:pao_gerd_f_valid'
    failure=.false.

    CPPrecondition(gf_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not.(associated(gf_env%scf_work1).and.&
         associated(gf_env%scf_work2).and.&
         associated(gf_env%scf_work3).and.&
         associated(gf_env%qs_env).and.&
         associated(gf_env%diag_matrix).and.&
         associated(gf_env%full_work_m).and.&
         associated(gf_env%min_work_m).and.&
         associated(gf_env%min_blacs_1).and.&
         associated(gf_env%min_blacs_2).and.&
         associated(gf_env%min_blacs_3))) then
       CPPrecondition(associated(gf_env%scf_work1),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%scf_work2),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%scf_work3),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%qs_env),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%diag_matrix),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%full_work_m),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%min_work_m),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%min_blacs_1),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%min_blacs_2),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%min_blacs_3),cp_warning_level,routineP,error,failure)
    end if
    res=.not.failure
  end function pao_gerd_f_valid
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_eval_at [1.0] *
!!
!!   NAME
!!     pao_gerd_f_eval_at
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_eval_at(gf_env, x, f, gradient, error)
!!       Type(pao_gerd_functional_type), Intent (INOUT), Target:: gf_env
!!       Real(Kind=wp), Intent (IN), Dimension(:):: x
!!       Real(Kind=wp), Intent (OUT), Optional:: f
!!       Real(Kind=wp), Intent (OUT), Optional, Dimension(:):: gradient
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_eval_at
!!
!!   FUNCTION
!!     evaluates the functional at one point (i.e. for the given angles
!!     and unitary transformation
!!
!!   NOTES
!!     min_p, and the density are twice the P used for the pao functional
!!     (closed shell electrons)
!!
!!   INPUTS
!!     gf_env: the gerd functional environement
!!     angles: angles defining the point where the functional is
!!            evaluated. If not given uses the ones in pao env
!!     f: the value of the functional
!!     df_du: the value of the derivative with respect to u
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
  subroutine pao_gerd_f_eval_at(gf_env,x, f, gradient , error)
    type(pao_gerd_functional_type), intent(inout), target :: gf_env
    real(kind=wp), intent(in), dimension(:) :: x
    real(kind=wp), intent(out), optional :: f
    real(kind=wp), intent(out), optional, dimension(:) :: gradient
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_obj_function:pao_gerd_f_eval_at'
    type(pao_glob_angles_type), pointer :: angles
    type(cp_block_matrix_type), pointer :: h_matrix, s_matrix
    type(cp_block_matrix_p_type), dimension(:), pointer :: h_matrixes, &
         s_matrixes
!         full_hamiltonian_m
    type(cp_error_type) :: iError
    type(global_environment_type), pointer :: g_env
    type(pao_calc_u_du_type) :: calc_u_du
    type(cp_matrix_block_iterator) :: block_iterator
    !type(cp_block_matrix_type) :: full_pao_m, min_pao_m

    type(real_matrix_p_type), dimension(:), pointer :: full_h_rmpv,full_s_rmpv
    type(cp_block_matrix_p_type), dimension(:), pointer :: min_d_smpv
    type(real_matrix_type), pointer :: min_density_m

    type(pao_env_type), pointer :: pao_env

    REAL(wp), DIMENSION(:), POINTER  :: mo_eigenvalues,&
                                        mo_occupation_numbers
    integer :: min_bas_size, full_bas_size, block_row, block_col,&
         min_bas, excl_bas, stat,size_null_space
    real(kind=wp) :: f1,f2,f3
    real(kind=wp), dimension(:,:), pointer :: du, block_val
    type(pao_local_angles_type) :: local_angles
    type(real_matrix_type), pointer :: diagMatrix
    call cp_error_init(iError, template_error=error)
    nullify(diagMatrix, mo_eigenvalues, mo_occupation_numbers,min_density_m,&
         min_d_smpv,full_h_rmpv,full_s_rmpv, h_matrixes, s_matrixes, h_matrix,&
         s_matrix, angles, pao_env)
    g_env => gf_env%global_env
    failure=.false.
    
    CPAssert(cp_valid(gf_env,error=error),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       gf_env%n_eval = gf_env%n_eval+1
       call get_qs_env(gf_env%qs_env,pao_env=pao_env,&
            k=full_h_rmpv, s=full_s_rmpv)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
       if (cp_debug .and. debug_this_module) then
          call cp_1d_r_output(logger=cp_error_get_logger(error),&
               outputName='obj_f_x',fromWhere=routineP, iter=gf_env%n_eval,&
               array=x,error=error)
       end if
    end if
    if (.not.failure) then
       call pao_qs_env_get(gf_env%qs_env,global_env=g_env,pao_angles_att=angles,&
            pao_min_density_m=min_d_smpv,&
            pao_tot_full_bas=full_bas_size,&
            pao_tot_min_bas=min_bas_size, error=error)
       CPPostcondition(associated(angles),cp_failure_level,routineP,error,failure)
       CPPrecondition(size(min_d_smpv)==1,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       min_density_m => min_d_smpv(1)%matrix%sm
       call pao_g_ang_set(angles,qs_env=gf_env%qs_env,global_env=g_env,&
            angles=x,error=error)
       call pao_g_ang_did_change(angles, qs_env=gf_env%qs_env,&
            global_env=g_env,&
            angles_changed=.true., s_changed=.false.,&
            h_changed=.false., full_reset=.false., error=error)
       call pao_g_ang_get_min_m(angles, min_hamiltonian_m=h_matrixes,&
            min_overlap_m=s_matrixes, qs_env=gf_env%qs_env, global_env=g_env,&
            error=error)
       h_matrix => h_matrixes(1)%matrix
       s_matrix => s_matrixes(1)%matrix
       if (cp_debug .and. debug_this_module) then
          call cp_sm_output(logger=cp_error_get_logger(error), &
               outputName="obj_f_min_h", &
               fromWhere=routineP, iter=gf_env%n_eval, matrix=h_matrix%sm,&
               para_env=g_env%para_env,error=error)
          call cp_sm_output(logger=cp_error_get_logger(error), &
               outputName="obj_f_min_s", &
               fromWhere=routineP, iter=gf_env%n_eval, matrix=s_matrix%sm,&
               para_env=g_env%para_env,error=error)
          call cp_sm_output(logger=cp_error_get_logger(error), &
               outputName="obj_f_nui", &
               fromWhere=routineP, iter=gf_env%n_eval, &
               matrix=pao_env%angles_att%NUi_injection%sm,&
               para_env=g_env%para_env,error=error)
       end if


       ! calculates the new density matrix with diagonalization
       ! do not do it every time?
       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(s_matrix),&
            blacs_matrix=gf_env%min_blacs_1)
       call calculate_ortho_matrix(ortho=gf_env%min_blacs_1,&
            work=gf_env%min_blacs_2,ndep=size_null_space, globenv=g_env)
       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(h_matrix),&
            blacs_matrix=gf_env%min_blacs_2)
       call orthogonalise_matrix(ortho=gf_env%min_blacs_1,&
            ks=gf_env%min_blacs_2,&
            work=gf_env%min_blacs_3,&
            globenv=g_env)
       
       ! min_blacs1 = S^{-1/2}, min_blacs2= H_{KS}
       call eigensolver(ks=gf_env%min_blacs_2,&
            mo_set=gf_env%mo_set(1)%mo_set,&
            ortho=gf_env%min_blacs_1, work=gf_env%min_blacs_3,&
            globenv=g_env)

       call calculate_density_matrix(mo_set=gf_env%mo_set(1)%mo_set,&
            density_matrix=gf_env%min_blacs_2,&
            globenv=g_env)
       
       call copy_blacs_to_sparse_matrix(gf_env%min_blacs_2,&
            min_density_m)

       call symmetrise_diagonal_blocks(min_density_m)

       CALL get_mo_set(mo_set=gf_env%mo_set(1)%mo_set,&
                       eigenvalues=mo_eigenvalues,&
                       occupation_numbers=mo_occupation_numbers)

       if (present(f)) then
          f=dot_product(mo_occupation_numbers,mo_eigenvalues)/2.0_wp
          if (cp_debug .and. debug_this_module) then
             call cp_sm_output(logger=cp_error_get_logger(error), &
                  outputName="obj_f_min_p", &
                  fromWhere=routineP, iter=gf_env%n_eval,&
                  matrix=min_density_m,&
                  para_env=g_env%para_env,error=error)
             call cp_write_output(logger=cp_error_get_logger(error),&
                  fromWhere=routineP,&
                  outputName="obj_f_f", iter=gf_env%n_eval,&
                  message="f ="//cp_to_string(f), error=error)
          end if
       end if

       if (present(gradient)) then
          ! uses the filtred P matrix
          ! should change it soon
          
          diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
          ! ** gf_env%diag_matrix = df/du_{i,j}
                    
          ! 2 H B P i^*
          call copy_sparse_to_blacs_matrix(full_h_rmpv(1)%matrix,&
               gf_env%scf_work1)
          call cp_m2f_with_B_M_pi(min_density_m,cp_get_matrix(gf_env%full_work_m),&
               glob_angles=angles, qs_env=gf_env%qs_env, global_env=g_env,&
               error=error)
          call copy_sparse_to_blacs_matrix(cp_get_matrix(gf_env%full_work_m),&
               gf_env%scf_work2)
          call cp_fm_symm(side="L",uplo="U",m=full_bas_size,n=full_bas_size,&
               alpha=1.0_wp,matrix_a=gf_env%scf_work1,&
               matrix_b=gf_env%scf_work2,&
               beta=0.0_wp,matrix_c=gf_env%scf_work3)
          call cp_set_to(diagMatrix,value=0.0_wp)
          call add_blacs_to_block_diag_sm(gf_env%scf_work3,&
               diagMatrix, para_env=g_env%para_env,&
               d_struct=gf_env%d_struct, error=error)
          if (cp_debug .and. debug_this_module) then
             call cp_sm_output(logger=cp_error_get_logger(error), &
               outputName="obj_f_BPpi", &
               fromWhere=routineP, iter=gf_env%n_eval, &
               matrix=cp_get_matrix(gf_env%full_work_m),&
               para_env=g_env%para_env,error=error)
             call cp_sm_output(logger=cp_error_get_logger(error), &
               outputName="obj_f_2HBPpi", &
               fromWhere=routineP, iter=gf_env%n_eval, matrix=diagMatrix,&
               para_env=g_env%para_env,error=error)
          end if
          ! 2 S B P minH P i^*
          call copy_sparse_to_blacs_matrix(h_matrix%sm,&
               gf_env%min_blacs_1)
          call symmetrise_blacs_matrix(gf_env%min_blacs_1,gf_env%min_blacs_2)
          call copy_sparse_to_blacs_matrix(min_density_m,&
               gf_env%min_blacs_2)
          call cp_fm_symm('R','U',min_bas_size,min_bas_size,&
               1.0_wp,gf_env%min_blacs_2,gf_env%min_blacs_1, 0.0_wp,&
               gf_env%min_blacs_3)
          call cp_fm_symm('L','U',min_bas_size,min_bas_size,&
               -0.5_wp,gf_env%min_blacs_2,gf_env%min_blacs_3, 0.0_wp,&
               gf_env%min_blacs_1)
          call copy_blacs_to_sparse_matrix(gf_env%min_blacs_1,&
               cp_get_matrix(gf_env%min_work_m))
          call cp_m2f_with_B_M_pi(cp_get_matrix(gf_env%min_work_m),&
               cp_get_matrix(gf_env%full_work_m),glob_angles=angles,&
               qs_env=gf_env%qs_env, global_env=g_env,&
               error=error)
          call copy_sparse_to_blacs_matrix(cp_get_matrix(gf_env%full_work_m),&
               gf_env%scf_work2)
          call copy_sparse_to_blacs_matrix(full_s_rmpv(1)%matrix,&
               gf_env%scf_work1)
          call cp_fm_symm('L','U',full_bas_size,full_bas_size,&
               1.0_wp,gf_env%scf_work1,gf_env%scf_work2,0.0_wp,&
               gf_env%scf_work3)
          call add_blacs_to_block_diag_sm(gf_env%scf_work3,&
               diagMatrix,para_env=g_env%para_env,&
               d_struct=gf_env%d_struct, error=error)
          ! diag_matrix = N diag_matrix = 2 N H B P i^* - 2 N S B P minH P i^*
          call cp_l_multiply_atomic_ortho(diagMatrix,qs_env=gf_env%qs_env,&
               global_env=g_env,error=error)
          
          if (cp_debug .and. debug_this_module) then
             call cp_sm_output(logger=cp_error_get_logger(error), &
                  outputName="obj_f_df_du", &
                  fromWhere=routineP, iter=gf_env%n_eval, matrix=diagMatrix,&
                  para_env=g_env%para_env,error=error)
          end if
          
          ! ** gradient = df/du_{i,j} du_{i,j}/dtheta_k
          gradient=0.0_wp
          call cp_init(block_iterator,matrix=cp_get_matrix(gf_env%diag_matrix))
          call pao_l_ang_init(local_angles, qs_env=gf_env%qs_env, &
               global_env=g_env, atom_nr=1,&
               glob_angles=angles, error=error)
          do
             if (.not.cp_next(block_iterator)) exit
             call cp_get(block_iterator, block_row=block_row,&
                  block_col=block_col, block_val=block_val)
             if (block_col==block_row) then
                call cp_set(local_angles, atom_nr=block_col)
                call cp_init(calc_u_du,local_angles)
                allocate(du(size(block_val,1),size(block_val,2)),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                if (.not.failure) then
                   do
                      if (.not.cp_next(calc_u_du)) exit
                      call pao_calc_u_du_get(calc_u_du, du=du, &
                           min_bas=min_bas, excl_bas=excl_bas)
                      gradient(cp_get_g_index(local_angles,&
                           min_bas=min_bas, excl_bas=excl_bas))=&
                           sum(du*block_val)
                   end do
                   call cp_dealloc_ref(calc_u_du)
                   deallocate(du,stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                end if
             end if
          end do
          call pao_l_ang_dealloc_ref(local_angles,error=error)
          call cp_sm_b_i_dealloc_ref(block_iterator,error=error)
          call mp_sum(gradient,g_env%group)
          if (cp_debug.and.debug_this_module) then
             call cp_1d_r_output(cp_error_get_logger(error),"obj_f_grad",&
                  routineP,gf_env%n_eval,gradient,error=error)
          end if
       end if
    end if
  end subroutine pao_gerd_f_eval_at
!***************************************************************************

!=========== qs_scf rip off (to avoid circular ref) to do =========

  SUBROUTINE calculate_ortho_matrix(ortho,work,ndep,globenv)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(cp_full_matrix_type), POINTER          :: ortho,work
    INTEGER, INTENT(OUT)                      :: ndep

!   *** Local variables ***

    INTEGER :: handle,nao
    real(kind=wp), parameter :: eps_eigval = 1.0E-5_wp, work_syevx = 1.0_wp

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL power_blacs_matrix(ortho,work,-0.5_wp,eps_eigval,ndep,work_syevx)
    CALL symmetrise_blacs_matrix(ortho,work)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ortho_matrix


  SUBROUTINE orthogonalise_matrix(ortho,ks,work,globenv)

!   Purpose: Orthogonalisation matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("orthogonalise_matrix","I","",handle)

    CALL cp_fm_get_info(matrix=ks,nrow_global=nao)
    CALL cp_fm_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp,work)
    CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp,ks)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE orthogonalise_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work,globenv)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set

!   *** Local variables ***

    INTEGER :: handle,imo,nao,nmo
    integer, save :: iter=0

    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors
    REAL(wp), DIMENSION(:), POINTER  :: mo_eigenvalues
    character(len=*), parameter :: routineN="eigensolver",&
         routineP=moduleN//':'//routineN
!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    eigenvalues=mo_eigenvalues,&
                    eigenvectors=mo_eigenvectors)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,1.0_wp)

    CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
                    mo_eigenvectors)

    if (cp_debug .and. debug_this_module) then
       iter=iter+1
       call cp_1d_r_output(logger=cp_get_default_logger(),&
            outputName="obj_f_mo_eval", &
            fromWhere=routineP, iter=iter, array=mo_eigenvalues)
!       call cp_2d_r_output(logger=cp_get_default_logger(),&
!            outputName="obj_f_mo_evect", &
!            fromWhere=routineP, iter=iter, array=mo_eigenvectors)
    end if
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

!!****f* pao_obj_function/cp_l_compose_atomic_ortho [1.0] *
!!
!!   NAME
!!     cp_l_compose_atomic_ortho
!!
!!   FUNCTION
!!     left multiplies with the atomic orthogonalization
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to be multiplied, on output contains
!!       N matrix (where N is the atomic orthogonalization)
!!     - qs_env: the qs_environment that defines the atomic ortho
!!     - global_env: global environement (to communicate)
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
subroutine cp_l_multiply_atomic_ortho(matrix,qs_env,global_env,error)
  type(real_matrix_type), pointer :: matrix
  ! pao_env might be modified
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_l_multiply_atomic_ortho',&
        routineP=moduleN//':'//routineN
  integer :: block_row, stat, max_full_bas
  real(kind=wp), dimension(:,:), pointer :: block_val, a_ortho
  real(kind=wp), dimension(:,:), allocatable :: t1_m
  type(cp_matrix_block_iterator) :: iterator
  failure=.false.
  nullify(block_val,a_ortho)
  
  CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call pao_qs_env_get(qs_env,global_env=global_env,&
          pao_max_full_bas=max_full_bas,error=error)
     allocate(t1_m(max_full_bas,max_full_bas),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  end if
  if (.not.failure) then
     call cp_sm_b_i_init(iterator,matrix=matrix,error=error)
     do
        if (.not. cp_next(iterator,block_row=block_row,&
             block_val=block_val,error=error)) exit
        a_ortho => pao_qs_env_get_atomic_ortho(qs_env=qs_env,&
             global_env=global_env,atom_nr=block_row,error=error)
        CPPrecondition(associated(block_val),cp_failure_level,routineP,error,failure)
        CPPrecondition(associated(a_ortho),cp_failure_level,routineP,error,failure)
        if (.not.failure) then
           call cp_assert(size(block_val,1)==size(a_ortho,2).and.&
                size(a_ortho,1)==size(a_ortho,2),cp_failure_level,&
                cp_assertion_failed,routineP,"size mismach in "//&
                CPSourceFileRef,&
                error=error,failure=failure)
           call cp_assert(size(t1_m,1)>=size(a_ortho,1).and.&
                size(t1_m,2)>=size(block_val,2),cp_failure_level,&
                cp_assertion_failed,routineP,"t1_m too small in "//&
                CPSourceFileRef,&
                error=error,failure=failure)           
        end if
        if (.not.failure) then
           call dgemm('N','N',size(a_ortho,1),size(block_val,2),size(a_ortho,2),&
                1.0_wp,a_ortho,size(a_ortho,1),block_val,size(block_val,1),&
                0.0_wp,t1_m,size(t1_m,1))
           block_val=t1_m(1:size(block_val,1),1:size(block_val,2))
        end if
     end do
     call cp_sm_b_i_dealloc_ref(iterator,error=error)
     deallocate(t1_m,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  end if
end subroutine cp_l_multiply_atomic_ortho
!***************************************************************************


end module pao_obj_function

