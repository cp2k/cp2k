!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_pao_obj_function [1.0] *
!!
!!   NAME
!!     qs_pao_obj_function
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
module qs_pao_obj_function
  use qs_log_handling
  use qs_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use qs_pao_types
  implicit none

  logical, private, parameter :: debug_this_module=.true.

  interface qs_valid
     module procedure qs_pao_gerd_f_valid
  end interface

!!***
!****************************************************************************

!!****s* qs_pao_obj_function/qs_pao_gerd_functional [1.0] *
!!
!!   NAME
!!     qs_pao_gerd_functional
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
!!     small basis, i an injection into the small basis,
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
!!           size of the small and full basis, and the pao env do not change
!!           between sucessive evaluations
!!     small_hamiltionian: the hamiltonian in the small basis
!!     small_overlap: the overlap in the small basis
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
  type qs_pao_gerd_functional
     logical :: initialized
     type(qs_pao_env_type), pointer :: pao_env
     type(qs_pao_matrix_type) :: small_hamiltonian, small_overlap
  end type qs_pao_gerd_functional
!!***
  !****************************************************************************
contains

!!****f* qs_pao_obj_function/qs_init_pao_gerd_f [1.0] *
!!
!!   NAME
!!     qs_init_pao_gerd_f
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
  subroutine qs_init_pao_gerd_f(gf_env,pao_env, error)
    type(qs_pao_gerd_functional), intent(out) :: gf_env
    type(qs_pao_env_type), intent(in), target :: pao_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_obj_function:qs_init_pao_gerd_f'
    failure=.false.

    QSPrecondition(qs_valid(pao_env),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       gf_env%pao_env => pao_env
       call qs_init(gf_env%small_hamiltonian,pao_env=pao_env,&
            min_bas_rows=.true., min_bas_cols=.true.)
       call qs_init(gf_env%small_overlap,pao_env=pao_env,&
            min_bas_rows=.true., min_bas_cols=.true.)
       gf_env%initialized=.true.
    end if
  end subroutine qs_init_pao_gerd_f
!!***
  !***************************************************************************

!!****f* qs_pao_obj_function/qs_dealloc_pao_gerd_f_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_pao_gerd_f_ref
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
  subroutine qs_dealloc_pao_gerd_f_ref(gf_env, error)
    type(qs_pao_gerd_functional), optional, intent(inout) :: gf_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_obj_function:qs_dealloc_pao_gerd_f_ref'
    failure=.false.

    QSPrecondition(gf_env%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       call qs_dealloc_ref(gf_env%small_overlap, error=error)
       call qs_dealloc_ref(gf_env%small_hamiltonian, error=error)
    end if
  end subroutine qs_dealloc_pao_gerd_f_ref
!!***
  !***************************************************************************

!!****f* qs_pao_obj_function/qs_pao_gerd_f_valid [1.0] *
!!
!!   NAME
!!     qs_pao_gerd_f_valid
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
  function qs_pao_gerd_f_valid(gf_env, error) result(res)
    logical ::res
    type(qs_pao_gerd_functional), intent(in) :: gf_env
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_obj_function:qs_pao_gerd_f_valid'
    failure=.false.

    QSPrecondition(gf_env%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(gf_env%pao_env),qs_warning_level,routineN,error,failure)
    res=.not.failure
  end function qs_pao_gerd_f_valid
!!***
  !***************************************************************************

!!****f* qs_pao_obj_function/qs_pao_gerd_f_eval_at [1.0] *
!!
!!   NAME
!!     qs_pao_gerd_f_eval_at
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
!!     unitary_t_in: the unitary transformation (if not given recalculates
!!            it)
!!     unitary_t_out: writes the unitary transformation in it
!!     f: the value of the functional
!!     df_du: the value of the derivative with respect to u
!!     density_matrix: the density in the small basis (if not given uses
!!            the one in pao_env)
!!     hamiltonian: the hamiltonian in the full basis (if not given uses 
!!            the one in pao_env)
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
  subroutine qs_pao_gerd_f_eval_at(gf_env, f, df_du, density_matrix,&
       unitary_t, error) ! to do
    type(qs_pao_gerd_functional), intent(in) :: gf_env
    real(kind=wp), intent(out), optional :: f
    type(qs_pao_matrix_type), intent(out), optional :: df_du
    type(qs_pao_matrix_type), optional, intent(in), target :: density_matrix, &
         unitary_t
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_obj_function:qs_pao_gerd_f_eval_at'
    type(qs_pao_matrix_type), pointer :: d_matrix, h_matrix, s_matrix, u_matrix
    type(qs_error_type) :: iError
    integer :: n_atom
    call qs_init(iError, template_error=error)
    failure=.false.

    QSAssert(qs_valid(gf_env,error=error),qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       call qs_get(gf_env%pao_env,min_density_m=d_matrix, &
            full_hamiltonian_m=h_matrix, full_overlap_m=s_matrix, &
            n_atom=n_atom, error=error)
       if (present(density_matrix)) &
            d_matrix => density_matrix
       if (present(unitary_t)) &
            u_matrix => unitary_t
       QSPrecondition(associated(d_matrix),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(u_matrix),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(h_matrix),qs_failure_level,routineN,error,failure)
       QSPrecondition(associated(s_matrix),qs_failure_level,routineN,error,failure)
       !call qs_pao_matrix_transf_to_m(h_matrix, gf_env%small_hamiltonian, &
       !     u_matrix, error=iError)
       !call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
       !     error,failure)
       if (.not.failure) then
          !call qs_pao_matrix_transf_to_m(s_matrix, gf_env%small_overlap, &
          !     u_matrix, error=iError)
          !call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
          !     error,failure)
       end if
       ! to do

    end if
  end subroutine qs_pao_gerd_f_eval_at
!!***
  !***************************************************************************

end module qs_pao_obj_function

