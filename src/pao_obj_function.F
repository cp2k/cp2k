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
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use pao_types
  use pao_env_methods
  use cp_sparse_matrix
  use pao_glob_angles_methods
  use pao_l_angles_methods
  use pao_proj_methods
  use qs_blacs, only: blacs_matrix_type,allocate_blacs_matrix,&
       copy_sparse_to_blacs_matrix,deallocate_blacs_matrix,&
       blacs_trace, add_blacs_to_block_diag_sm, blacs_diag_mult,&
       copy_blacs_to_blacs_matrix, blacs_gemm, copy_blacs_to_sparse_matrix,&
       blacs_syrk,blacs_syevx,blacs_symm,get_blacs_matrix_info
  use sparse_matrix_types, only: symmetrise_diagonal_blocks
  use timings, only: timeset, timestop
  USE qs_mo_types, ONLY: calculate_density_matrix,&
                         get_mo_set,&
                         mo_set_type
  use message_passing, only: mp_sum
  use qs_environment_types
  use pao_unitary_t
  use pao_qs_env_methods
  use cp_matrix_utils
  use cp_b_matrix_structure
  use sparse_matrix_types
  use pao_unitary_t
  use pao_utils
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
     logical :: initialized
     type(qs_environment_type), pointer :: qs_env
     type(blacs_matrix_type), POINTER :: min_hamiltonian_m, min_overlap_m,&
          min_density_m, work_m
     type(blacs_matrix_type), pointer :: scf_work1, scf_work2, scf_work3
     integer :: blacs_context
     type(mo_set_type), pointer :: alpha_mo
     type(cp_sparse_matrix_type), pointer :: diag_matrix, full_work_m, min_work_m
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
!!         scf_work2, scf_work3, mo_orbital, error)
!!       Type(pao_gerd_functional_type), Intent (OUT):: gf_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(blacs_matrix_type), Intent (IN), Target:: scf_work1,&
!!         scf_work2, scf_work3
!!       Type(mo_set_type), Target, Intent (IN):: mo_orbital
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
       scf_work3, mo_orbital, error)
    type(pao_gerd_functional_type), intent(out) :: gf_env
    type(qs_environment_type), intent(inout), target :: qs_env ! modifies pao_env
    type(global_environment_type), intent(in),target :: global_env
    type(blacs_matrix_type), intent(in), target :: scf_work1, scf_work2,&
         scf_work3
    type(mo_set_type), target, intent(in) :: mo_orbital
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_obj_function:pao_gerd_f_init'
    integer :: ntot_min_bas, ntot_full_bas
    type(real_matrix_p_type), pointer:: full_overlap_m
    type(real_matrix_type), pointer :: diagMatrix
    type(cp_b_matrix_struct_type), pointer :: min_bas_m_struct, full_bas_m_struct
    type(pao_env_type), pointer :: pao_env
    failure=.false.

    gf_env%qs_env => qs_env
    if (.not.failure) then
#ifdef __parallel
       call blacs_get(0,0,gf_env%blacs_context) ! the default context
#else
       gf_env%blacs_context=0
#endif
       call get_qs_env(gf_env%qs_env, s=full_overlap_m, pao_env=pao_env)
       call pao_qs_env_get(gf_env%qs_env, global_env=global_env,&
            pao_tot_full_bas=ntot_full_bas,&
            pao_tot_min_bas=ntot_min_bas,&
            error=error)
       call cp_get(pao_env,min_bas_m_struct=min_bas_m_struct,&
            full_bas_m_struct=full_bas_m_struct, error=error)
       call allocate_blacs_matrix(gf_env%min_hamiltonian_m,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min hamiltonian',&
            context=gf_env%blacs_context,&
            globenv=global_env)
       call allocate_blacs_matrix(gf_env%min_overlap_m,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min overlap',&
            context=gf_env%blacs_context,&
            globenv=global_env)
       call allocate_blacs_matrix(gf_env%min_density_m,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min density',&
            context=gf_env%blacs_context,&
            globenv=global_env)
       call allocate_blacs_matrix(gf_env%work_m,&
            nrow_global=ntot_min_bas,&
            ncol_global=ntot_min_bas,&
            nrow_block=cp_optimal_blacs_block_size,&
            ncol_block=cp_optimal_blacs_block_size,&
            name='min work matrix',&
            context=gf_env%blacs_context,&
            globenv=global_env)
       gf_env%scf_work1 => scf_work1
       gf_env%scf_work2 => scf_work2
       gf_env%scf_work3 => scf_work3
       gf_env%alpha_mo => mo_orbital
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
    integer :: stat
    failure=.false.

    CPPrecondition(gf_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call deallocate_blacs_matrix(gf_env%min_overlap_m)
       call deallocate_blacs_matrix(gf_env%min_hamiltonian_m)
       call deallocate_blacs_matrix(gf_env%min_density_m)
       call deallocate_blacs_matrix(gf_env%work_m)
       nullify(gf_env%scf_work1, gf_env%scf_work2, gf_env%scf_work3)
       call cp_dealloc(gf_env%min_work_m, error=error)
       call cp_dealloc(gf_env%full_work_m, error=error)
       call cp_dealloc(gf_env%diag_matrix, error=error)
       deallocate(gf_env%d_struct, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
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
         associated(gf_env%min_work_m))) then
       CPPrecondition(associated(gf_env%scf_work1),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%scf_work2),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%scf_work3),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%qs_env),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%diag_matrix),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%full_work_m),cp_warning_level,routineP,error,failure)
       CPPrecondition(associated(gf_env%min_work_m),cp_warning_level,routineP,error,failure)
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
    type(cp_sparse_matrix_type), pointer :: d_matrix, h_matrix, s_matrix, &
         full_hamiltonian_m
    type(cp_error_type) :: iError
    type(global_environment_type), pointer :: g_env
    type(pao_calc_u_du_type) :: calc_u_du
    type(cp_matrix_block_iterator) :: block_iterator
    type(cp_sparse_matrix_type) :: full_pao_m, min_pao_m
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
!!FM    if (.not.failure) then
!!FM       ! get actual min_hamiltonian_m and min_overlap_m
!!FM       call cp_get(gf_env%qs_env,angles_att=angles, global_env=g_env,&
!!FM            min_density_m=d_matrix, full_hamiltonian_m=full_hamiltonian_m,&
!!FM            n_atom=n_atom, ntot_full_bas=full_bas_size,&
!!FM            ntot_min_bas=min_bas_size, error=error)
!!FM       CPPostcondition(associated(angles),cp_failure_level,routineP,error,failure)
!!FM    end if
!!FM    if (.not.failure) then
!!FM       call cp_set(angles,angles=x,error=error)
!!FM       call cp_did_change(angles, angles_changed=.true., s_changed=.false.,&
!!FM            h_changed=.false., full_reset=.false., error=error)
!!FM       call cp_get(angles, min_hamiltonian_m=h_matrix,&
!!FM            min_overlap_m=s_matrix,&
!!FM            error=error)
!!FM
!!FM       ! calculates the new density matrix with diagonalization
!!FM       ! do not do it every time?
!!FM       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(s_matrix),&
!!FM            blacs_matrix=gf_env%min_overlap_m, &
!!FM            context=gf_env%blacs_context, &
!!FM            globenv=g_env)
!!FM       
!!FM       call orthogonalise_matrix(gf_env%min_overlap_m,gf_env%min_density_m,&
!!FM            gf_env%work_m,gf_env%blacs_context,g_env)
!!FM       
!!FM       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(h_matrix),&
!!FM            blacs_matrix=gf_env%min_density_m, &
!!FM            context=gf_env%blacs_context, &
!!FM            globenv=g_env)
!!FM
!!FM       call eigensolver(gf_env%min_density_m,gf_env%alpha_mo,&
!!FM            gf_env%min_overlap_m,gf_env%work_m,gf_env%blacs_context,g_env)
!!FM       
!!FM       call calculate_density_matrix(gf_env%alpha_mo,gf_env%min_density_m,&
!!FM            gf_env%blacs_context,g_env)
!!FM       
!!FM       call copy_blacs_to_sparse_matrix(gf_env%min_density_m,&
!!FM            cp_get_matrix(d_matrix),gf_env%blacs_context,g_env)
!!FM
!!FM       call symmetrise_diagonal_blocks(cp_get_matrix(d_matrix))
!!FM
!!FM       ! could be avoided by having two more qs_blacs matrices
!!FM       ! (copy from qs_blacs to qs_blacs does not need communication)
!!FM       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(h_matrix),&
!!FM            blacs_matrix=gf_env%min_hamiltonian_m, &
!!FM            context=gf_env%blacs_context, &
!!FM            globenv=g_env)
!!FM       call copy_sparse_to_blacs_matrix(sparse_matrix=cp_get_matrix(s_matrix),&
!!FM            blacs_matrix=gf_env%min_overlap_m, &
!!FM            context=gf_env%blacs_context, &
!!FM            globenv=g_env)
!!FM
!!FM       CALL get_mo_set(mo_set=gf_env%alpha_mo,&
!!FM                       eigenvalues=mo_eigenvalues,&
!!FM                       occupation_numbers=mo_occupation_numbers)
!!FM
!!FM       if (present(f)) then
!!FM          f=dot_product(mo_occupation_numbers,mo_eigenvalues)
!!FM! good test: check if equal
!!FM!!$          ! Tr[P (N U i)^* H (N U i)]
!!FM!!$          call blacs_trace(gf_env%min_density_m,gf_env%min_hamiltonian_m,f1,&
!!FM!!$               gf_env%blacs_context, g_env)
!!FM!!$          ! Tr[V(P (N U i)^* S (N U i) P)]
!!FM!!$          call blacs_gemm('N','N',min_bas_size,min_bas_size,min_bas_size,&
!!FM!!$               1.0_wp, gf_env%min_density_m,gf_env%min_overlap_m,0.0_wp,&
!!FM!!$               gf_env%work_m, gf_env%blacs_context,g_env)
!!FM!!$          call blacs_diag_mult(mo_eigenvalues,gf_env%work_m)
!!FM!!$          call blacs_trace(gf_env%min_density_m,gf_env%work_m,f2,&
!!FM!!$               gf_env%blacs_context,g_env)
!!FM!!$          ! Tr[VP]
!!FM!!$          f3=cp_trace(mo_eigenvalues,d_matrix)
!!FM!!$          f=f1-f2+f3
!!FM       end if
!!FM
!!FM       if (present(gradient)) then
!!FM          diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
!!FM          ! ** gf_env%diag_matrix = df/du_{i,j}
!!FM                    
!!FM          ! N H B P i^*
!!FM          call copy_sparse_to_blacs_matrix(cp_get_matrix(full_hamiltonian_m),&
!!FM               gf_env%scf_work1,gf_env%blacs_context,g_env)
!!FM          call cp_m2f_with_B_M_pi(d_matrix,gf_env%full_work_m,angles,&
!!FM               error=error)
!!FM          call copy_sparse_to_blacs_matrix(cp_get_matrix(gf_env%full_work_m),&
!!FM               gf_env%scf_work2, gf_env%blacs_context,g_env)
!!FM          call blacs_gemm('N','N',full_bas_size,full_bas_size,full_bas_size,&
!!FM               1.0_wp,gf_env%scf_work1,gf_env%scf_work2,0.0_wp,&
!!FM               gf_env%scf_work3,gf_env%blacs_context,g_env)
!!FM          call cp_set_to(diagMatrix,value=0.0_wp)
!!FM          call add_blacs_to_block_diag_sm(gf_env%scf_work3,&
!!FM               diagMatrix,&
!!FM               global_env=g_env,d_struct=gf_env%d_struct, error=error)
!!FM          ! N S B P V P i^*
!!FM          call blacs_diag_mult(mo_eigenvalues,&
!!FM               gf_env%min_density_m,gf_env%work_m, global_env=g_env)
!!FM          call blacs_gemm('N','N',min_bas_size,min_bas_size,min_bas_size,&
!!FM               1.0_wp,gf_env%min_density_m,gf_env%work_m, 0.0_wp,&
!!FM               gf_env%min_hamiltonian_m, context=gf_env%blacs_context,&
!!FM               globenv=g_env)
!!FM          call copy_blacs_to_sparse_matrix(gf_env%min_hamiltonian_m,&
!!FM               cp_get_matrix(gf_env%min_work_m),gf_env%blacs_context, g_env)
!!FM          call cp_m2f_with_B_M_pi(gf_env%min_work_m,&
!!FM               gf_env%full_work_m,glob_angles=angles,&
!!FM               error=error)
!!FM          call cp_dealloc_ref(min_pao_m,error=error)
!!FM          call cp_dealloc_ref(full_pao_m,error=error)
!!FM
!!FM          call copy_sparse_to_blacs_matrix(cp_get_matrix(gf_env%full_work_m),&
!!FM               gf_env%scf_work2, gf_env%blacs_context, g_env)
!!FM          call copy_sparse_to_blacs_matrix(cp_get_matrix(s_matrix),&
!!FM               gf_env%scf_work1, gf_env%blacs_context, g_env)
!!FM          call blacs_gemm('N','N',full_bas_size,full_bas_size,full_bas_size,&
!!FM               1.0_wp,gf_env%scf_work1,gf_env%scf_work2,0.0_wp,&
!!FM               gf_env%scf_work3, gf_env%blacs_context, g_env)
!!FM          call add_blacs_to_block_diag_sm(gf_env%scf_work3,&
!!FM               diagMatrix,&
!!FM               global_env=g_env,d_struct=gf_env%d_struct, error=error)
!!FM          
!!FM          ! ** gradient = df/du_{i,j} du_{i,j}/dtheta_k
!!FM          gradient=0.0_wp
!!FM          call cp_init(block_iterator,matrix=cp_get_matrix(gf_env%diag_matrix))
!!FM          do
!!FM             if (.not.cp_next(block_iterator)) exit
!!FM             call cp_get(block_iterator, block_row=block_row,&
!!FM                  block_col=block_col, block_val=block_val)
!!FM             if (block_col==block_row) then
!!FM                call cp_set(local_angles, atom_nr=block_col)
!!FM                call cp_init(calc_u_du,local_angles)
!!FM                do
!!FM                   if (.not.cp_next(calc_u_du)) exit
!!FM                   call cp_get(calc_u_du, du=du, &
!!FM                        min_bas=min_bas, excl_bas=excl_bas)
!!FM                   gradient(cp_get_g_index(local_angles,&
!!FM                         min_bas=min_bas, excl_bas=excl_bas))=&
!!FM                         sum(du*block_val)
!!FM                end do
!!FM                call cp_dealloc_ref(calc_u_du)
!!FM             end if
!!FM          end do
!!FM          call mp_sum(gradient,g_env%group)
!!FM       end if
!!FM    end if
  end subroutine pao_gerd_f_eval_at
!***************************************************************************


!=========== qs_scf rip off (to avoid circulare ref) to do =========

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

