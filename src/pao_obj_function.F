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
       cp_log, cp_logger_type, cp_note_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_error_dealloc_ref, cp_error_get_logger, cp_init
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use pao_types, only: pao_env_type, pao_glob_angles_type, cp_get
!  use pao_env_methods
  use cp_sparse_matrix, only: cp_sparse_matrix_type, cp_get_matrix,&
       cp_create, cp_dealloc, cp_sparse_matrix_p_type
  use pao_glob_angles_methods, only: pao_g_ang_set, pao_g_ang_did_change
  use pao_l_angles_methods, only: pao_local_angles_type, cp_get_angle,&
       pao_l_ang_init, pao_l_ang_get, cp_get_g_index, cp_set
  use pao_proj_methods, only: cp_m2f_with_B_M_pi,pao_g_ang_get_min_m
  use qs_blacs, only: blacs_matrix_type,allocate_blacs_matrix,&
       copy_sparse_to_blacs_matrix,deallocate_blacs_matrix,&
       blacs_trace, add_blacs_to_block_diag_sm, blacs_diag_mult,&
       copy_blacs_to_blacs_matrix, blacs_gemm, copy_blacs_to_sparse_matrix,&
       blacs_syrk,blacs_syevx,blacs_symm,get_blacs_matrix_info
  use sparse_matrix_types, only: symmetrise_diagonal_blocks
  use timings, only: timeset, timestop
  USE qs_mo_types, ONLY: calculate_density_matrix,&
                         get_mo_set,allocate_mo_set, deallocate_mo_set,&
                         mo_set_type, mo_set_p_type
  use message_passing, only: mp_sum
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use pao_unitary_t, only: pao_calc_u_du_type, cp_init, cp_next, &
       cp_dealloc_ref, pao_calc_u_du_get
  use pao_qs_env_methods, only: pao_qs_env_get!, cp_get_b_projection
  use cp_matrix_utils, only: cp_matrix_block_iterator, cp_matrix_b_i_init,&
       cp_next, cp_next, cp_matrix_b_i_dealloc_ref,cp_d_struct_create,&
       cp_init_diagonal_blocks, cp_get, cp_init, cp_set_to
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type,&
       get_matrix_info
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind_set
!  use pao_utils
  implicit none
  private
  
  logical, private, parameter :: debug_this_module=.true.
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
!!     initialized: true if initialized (bug catcher)
!!     pao_env: a pointer to the actual pao_env. It is assumed that the
!!           size of the min and full basis, and the pao env do not change
!!           between sucessive evaluations
!!     min_hamiltonian_m: qs_blacs matrix with the hamiltonian in the min basis
!!     min_overlap_m: qs_blacs matrix with the overlap in the min basis
!!     density_m: qs_blacs matrix with the density matrix
!!     work_m: qs_blacs matrix in the small basis
!!     blacs_context: the qs_blacs context
!!     alpha_mo: molecular orbitals
!!     scf_work1, scf_work2, scf_work3: qs_blacs matrices in the big basis
!!     diag_matrix: diagonal sparse matrix
!!     full_work_m: sparse matrix in the full basis
!!     min_work_m: sparse matrix in the minimal basis
!!     d_struct: the distribution structure of the diagonal blocks between
!!          the processors (d_struct(i) gives on which processor is blocks i)
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
     logical :: initialized, should_dealloc_mo
     type(qs_environment_type), pointer :: qs_env
     type(blacs_matrix_type), POINTER :: min_blacs_1, min_blacs_2,&
          min_blacs_3
     type(blacs_matrix_type), pointer :: scf_work1, scf_work2, scf_work3
     integer :: blacs_context
     type(mo_set_p_type), dimension(:), pointer :: mo_set
     type(cp_sparse_matrix_type), pointer :: diag_matrix, full_work_m,&
          min_work_m
     integer, dimension(:), pointer :: d_struct
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
!!       Type(blacs_matrix_type), Intent (IN), Target:: scf_work1,&
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
    type(blacs_matrix_type), intent(in), target :: scf_work1, scf_work2,&
         scf_work3
    type(mo_set_p_type), dimension(:), target, intent(in),optional ::&
         mo_orbital
    logical, intent(in), optional :: should_dealloc_mo
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_obj_function:pao_gerd_f_init'
    integer :: ntot_min_bas, ntot_full_bas, i, stat, nelectron, n_atoms
    type(real_matrix_p_type), pointer:: full_overlap_m
    type(real_matrix_type), pointer :: diagMatrix
    type(cp_b_matrix_struct_type), pointer :: min_bas_m_struct, full_bas_m_struct
    type(pao_env_type), pointer :: pao_env
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    failure=.false.

    gf_env%qs_env => qs_env
    if (.not.failure) then
#ifdef __parallel
       call blacs_get(0,0,gf_env%blacs_context) ! the default context
#else
       gf_env%blacs_context=0
#endif
       call get_qs_env(gf_env%qs_env, s=full_overlap_m, pao_env=pao_env,&
            atomic_kind_set=atomic_kind_set, natoms=n_atoms)
       call pao_qs_env_get(gf_env%qs_env, global_env=global_env,&
            pao_tot_full_bas=ntot_full_bas,&
            pao_tot_min_bas=ntot_min_bas,&
            pao_min_bas_m_struct=min_bas_m_struct,&
            pao_full_bas_m_struct=full_bas_m_struct, error=error)
       call allocate_blacs_matrix(gf_env%min_blacs_1,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min blacs 1',&
            context=gf_env%blacs_context,&
            globenv=global_env)
       call allocate_blacs_matrix(gf_env%min_blacs_2,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min blacs 2',&
            context=gf_env%blacs_context,&
            globenv=global_env)
       call allocate_blacs_matrix(gf_env%min_blacs_3,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min blacs 3',&
            context=gf_env%blacs_context,&
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
               ncol_block=n_atoms,context=gf_env%blacs_context,&
               globenv=global_env)
          gf_env%should_dealloc_mo=.true.
       end if
       if (present(should_dealloc_mo)) &
            gf_env%should_dealloc_mo=should_dealloc_mo
       call cp_create(gf_env%diag_matrix,global_env=global_env,&
            matrix_structure=min_bas_m_struct)
       call cp_d_struct_create(gf_env%d_struct,&
            matrix=full_overlap_m%matrix,&
            global_env=global_env,error=error)
       diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
       call cp_init_diagonal_blocks(diagMatrix,&
            global_env=global_env,&
            default_value=0.0_wp,&
            d_struct=gf_env%d_struct)
       call cp_create(gf_env%min_work_m,global_env=global_env,&
            matrix_structure=min_bas_m_struct)
       call cp_create(gf_env%full_work_m,global_env=global_env,&
            matrix_structure=full_bas_m_struct)
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
       call deallocate_blacs_matrix(gf_env%min_blacs_1)
       call deallocate_blacs_matrix(gf_env%min_blacs_2)
       call deallocate_blacs_matrix(gf_env%min_blacs_3)
       nullify(gf_env%scf_work1, gf_env%scf_work2, gf_env%scf_work3)
       call cp_dealloc(gf_env%min_work_m, error=error)
       call cp_dealloc(gf_env%full_work_m, error=error)
       call cp_dealloc(gf_env%diag_matrix, error=error)
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
!!     -
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
    type(cp_sparse_matrix_type), pointer :: h_matrix, s_matrix
    type(cp_sparse_matrix_p_type), dimension(:), pointer :: h_matrixes, &
         s_matrixes
!         full_hamiltonian_m
    type(cp_error_type) :: iError
    type(global_environment_type), pointer :: g_env
    type(pao_calc_u_du_type) :: calc_u_du
    type(cp_matrix_block_iterator) :: block_iterator
    !type(cp_sparse_matrix_type) :: full_pao_m, min_pao_m

    type(real_matrix_p_type), pointer :: full_h_rmp, full_s_rmp
    type(cp_sparse_matrix_p_type), dimension(:), pointer :: min_d_rmpv
    type(real_matrix_type), pointer :: d_matrix

    type(pao_env_type), pointer :: pao_env

    REAL(wp), DIMENSION(:), POINTER  :: mo_eigenvalues,&
                                        mo_occupation_numbers
    integer :: n_atom, min_bas_size, full_bas_size, block_row, block_col,&
         min_bas, excl_bas
    real(kind=wp) :: f1,f2,f3
    real(kind=wp), dimension(:,:), pointer :: du, block_val
    type(pao_local_angles_type) :: local_angles
    type(real_matrix_type), pointer :: diagMatrix
    call cp_init(iError, template_error=error)
    nullify(diagMatrix)
    failure=.false.
    
    CPAssert(cp_valid(gf_env,error=error),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       ! get actual min_hamiltonian_m and min_overlap_m
       call get_qs_env(gf_env%qs_env,pao_env=pao_env,&
            h=full_h_rmp, s=full_s_rmp, natoms=n_atom)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       call pao_qs_env_get(gf_env%qs_env, global_env=g_env, pao_angles_att=angles,&
            pao_min_density_m=min_d_rmpv,&
            pao_tot_full_bas=full_bas_size,&
            pao_tot_min_bas=min_bas_size, error=error)
       CPPostcondition(associated(angles),cp_failure_level,routineP,error,failure)
       CPPrecondition(size(min_d_rmpv)==1,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       d_matrix => min_d_rmpv(1)%matrix%matrix
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

       ! calculates the new density matrix with diagonalization
       ! do not do it every time?
       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(s_matrix),&
            blacs_matrix=gf_env%min_blacs_1, &
            context=gf_env%blacs_context, &
            globenv=g_env)

       call orthogonalise_matrix(ortho=gf_env%min_blacs_1,&
            ks=gf_env%min_blacs_2,&
            work=gf_env%min_blacs_3,context=gf_env%blacs_context,&
            globenv=g_env)
       
       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(h_matrix),&
            blacs_matrix=gf_env%min_blacs_2, &
            context=gf_env%blacs_context, &
            globenv=g_env)

       ! min_blacs1 = S^{-1/2}, min_blacs2= H_{KS}
       call eigensolver(ks=gf_env%min_blacs_2,&
            mo_set=gf_env%mo_set(1)%mo_set,&
            ortho=gf_env%min_blacs_1, work=gf_env%min_blacs_3,&
            context=gf_env%blacs_context,globenv=g_env)

       call calculate_density_matrix(mo_set=gf_env%mo_set(1)%mo_set,&
            density_matrix=gf_env%min_blacs_2,&
            context=gf_env%blacs_context,globenv=g_env)
       
       call copy_blacs_to_sparse_matrix(gf_env%min_blacs_2,&
            d_matrix,gf_env%blacs_context,g_env)

       call symmetrise_diagonal_blocks(d_matrix)

       ! could be avoided by having two more qs_blacs matrices
       ! (copy from qs_blacs to qs_blacs does not need communication)
!!FM       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(h_matrix),&
!!FM            blacs_matrix=gf_env%min_blacs_1, &
!!FM            context=gf_env%blacs_context, &
!!FM            globenv=g_env)
!!FM       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(s_matrix),&
!!FM            blacs_matrix=gf_env%min_blacs_1,&
!!FM            context=gf_env%blacs_context,&
!!FM            globenv=g_env)
!!FM       ! min_blacs1=H_{KS} min_blacs2=density

       CALL get_mo_set(mo_set=gf_env%mo_set(1)%mo_set,&
                       eigenvalues=mo_eigenvalues,&
                       occupation_numbers=mo_occupation_numbers)

       if (present(f)) then
          f=dot_product(mo_occupation_numbers,mo_eigenvalues)
! good test: check if equal
!!$          ! Tr[P (N U i)^* H (N U i)]
!!$          call blacs_trace(gf_env%min_blacs_2,gf_env%min_hamiltonian_m,f1,&
!!$               gf_env%blacs_context, g_env)
!!$          ! Tr[V(P (N U i)^* S (N U i) P)]
!!$          call blacs_gemm('N','N',min_bas_size,min_bas_size,min_bas_size,&
!!$               1.0_wp, gf_env%min_blacs_2,gf_env%min_blacs_1,0.0_wp,&
!!$               gf_env%min_blacs_3, gf_env%blacs_context,g_env)
!!$          call blacs_diag_mult(mo_eigenvalues,gf_env%min_blacs_3)
!!$          call blacs_trace(gf_env%min_blacs_2,gf_env%min_blacs_3,f2,&
!!$               gf_env%blacs_context,g_env)
!!$          ! Tr[VP]
!!$          f3=cp_trace(mo_eigenvalues,d_matrix)
!!$          f=f1-f2+f3
       end if

       if (present(gradient)) then
          diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
          ! ** gf_env%diag_matrix = df/du_{i,j}
                    
          ! N H B P i^*
          call copy_sparse_to_blacs_matrix(full_h_rmp%matrix,&
               gf_env%scf_work1,gf_env%blacs_context,g_env)
          call cp_m2f_with_B_M_pi(d_matrix,cp_get_matrix(gf_env%full_work_m),&
               glob_angles=angles, qs_env=gf_env%qs_env, global_env=g_env,&
               error=error)
          call copy_sparse_to_blacs_matrix(cp_get_matrix(gf_env%full_work_m),&
               gf_env%scf_work2, gf_env%blacs_context,g_env)
          call blacs_gemm('N','N',full_bas_size,full_bas_size,full_bas_size,&
               1.0_wp,gf_env%scf_work1,gf_env%scf_work2,0.0_wp,&
               gf_env%scf_work3,gf_env%blacs_context,g_env)
          call cp_set_to(diagMatrix,value=0.0_wp)
          call add_blacs_to_block_diag_sm(gf_env%scf_work3,&
               diagMatrix,&
               global_env=g_env,d_struct=gf_env%d_struct, error=error)
          ! N S B P V P i^*
          call blacs_diag_mult(mo_eigenvalues,&
               gf_env%min_blacs_2,gf_env%min_blacs_3, global_env=g_env)
          call blacs_gemm('N','N',min_bas_size,min_bas_size,min_bas_size,&
               1.0_wp,gf_env%min_blacs_2,gf_env%min_blacs_3, 0.0_wp,&
               gf_env%min_blacs_1, context=gf_env%blacs_context,&
               globenv=g_env)
          call copy_blacs_to_sparse_matrix(gf_env%min_blacs_1,&
               cp_get_matrix(gf_env%min_work_m),gf_env%blacs_context, g_env)
          call cp_m2f_with_B_M_pi(cp_get_matrix(gf_env%min_work_m),&
               cp_get_matrix(gf_env%full_work_m),glob_angles=angles,&
               qs_env=gf_env%qs_env, global_env=g_env,&
               error=error)
          !call cp_dealloc_ref(min_pao_m,error=error)
          !call cp_dealloc_ref(full_pao_m,error=error)

          call copy_sparse_to_blacs_matrix(cp_get_matrix(gf_env%full_work_m),&
               gf_env%scf_work2, gf_env%blacs_context, g_env)
          call copy_sparse_to_blacs_matrix(cp_get_matrix(s_matrix),&
               gf_env%scf_work1, gf_env%blacs_context, g_env)
          call blacs_gemm('N','N',full_bas_size,full_bas_size,full_bas_size,&
               1.0_wp,gf_env%scf_work1,gf_env%scf_work2,0.0_wp,&
               gf_env%scf_work3, gf_env%blacs_context, g_env)
          call add_blacs_to_block_diag_sm(gf_env%scf_work3,&
               diagMatrix,&
               global_env=g_env,d_struct=gf_env%d_struct, error=error)
          
          ! ** gradient = df/du_{i,j} du_{i,j}/dtheta_k
          gradient=0.0_wp
          call cp_init(block_iterator,matrix=cp_get_matrix(gf_env%diag_matrix))
          do
             if (.not.cp_next(block_iterator)) exit
             call cp_get(block_iterator, block_row=block_row,&
                  block_col=block_col, block_val=block_val)
             if (block_col==block_row) then
                call cp_set(local_angles, atom_nr=block_col)
                call cp_init(calc_u_du,local_angles)
                do
                   if (.not.cp_next(calc_u_du)) exit
                   call pao_calc_u_du_get(calc_u_du, du=du, &
                        min_bas=min_bas, excl_bas=excl_bas)
                   gradient(cp_get_g_index(local_angles,&
                         min_bas=min_bas, excl_bas=excl_bas))=&
                         sum(du*block_val)
                end do
                call cp_dealloc_ref(calc_u_du)
             end if
          end do
          call mp_sum(gradient,g_env%group)
       end if
    end if
  end subroutine pao_gerd_f_eval_at
!***************************************************************************


!=========== qs_scf rip off (to avoid circular ref) to do =========

  SUBROUTINE orthogonalise_matrix(ortho,ks,work,context,globenv)

!   Purpose: Orthogonalisation matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("orthogonalise_matrix","I","",handle)

    CALL get_blacs_matrix_info(matrix=ks,nrow_global=nao)
    CALL blacs_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp,work,context,&
                    globenv)
    CALL blacs_gemm("T","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp,ks,context,&
                    globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE orthogonalise_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work,context,globenv)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,imo,nao,nmo

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors
    REAL(wp), DIMENSION(:), POINTER  :: mo_eigenvalues

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    eigenvalues=mo_eigenvalues,&
                    eigenvectors=mo_eigenvectors)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,1.0_wp,context,globenv)

    CALL blacs_gemm("N","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
                    mo_eigenvectors,context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

end module pao_obj_function

