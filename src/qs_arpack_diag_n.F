!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_arpack_diag_n *
!!
!!   NAME
!!     qs_arpack_diag_n
!!
!!   SYNOPSIS
!!     A module to compute the eigenvalues and eigenvectors of a
!!     general (non-symmetric) operator. The operator can be arbitrary
!!     as long as the action of the operator on a vector is specified
!!     (in a separate module like, for example, qs_arpack_full).
!!
!!   MODIFICATION HISTORY
!!     09-2002, TCH, created
!!
!!   SEE ALSO
!!     qs_arpack_diag_sy.F
!!
!!****
module qs_arpack_diag_n

  use global_types, only : global_environment_type
  use kinds,        only : dbl
  use qs_blacs,     only : cp_full_matrix_type, cp_fm_get_info, &
                           blacs_find_vector_distribution
  use termination,  only : stop_program

  ! the operator environments and operators
  use qs_p_types,       only : qs_p_env_type
  use qs_arpack_tddfpt, only : arpack_op_tddfpt

  implicit none

!!!  interface arpack_op
!!!     module procedure arpack_op_tddfpt
!!!  end interface

  character(len=*), parameter :: moduleN = "qs_arpack_diag_n"

  private
  public :: arpack_diag_n

contains

! *****************************************************************************

!!****f* cp2k/qs_arpack_diag_n/arpack_diag_n *
!!
!!   NAME
!!     arpack_diag_n
!!
!!   FUNCTION
!!     Given action of an general (non-symmetric) operator A on a vector v,
!!     computes the eigenvalues and eigenvectors of A.
!!
!!   ARGUMENTS
!!     - n_ev      : number of eigenvalues and eigenvectors to calculate
!!     - eval      : array of eigenvalues
!!     - evec      : on input initial guess - on output eigenvectors (matrix)
!!     - glob_env  : global_environment (needed for communication)
!!     ONE of the following optional arguments must also be specified
!!     - op_tddfpt : operator for the tddfpt calculations
!!
!!   NOTES
!!     - Since the problem at hand is not symmetric (hermitian), there
!!       can be complex eigenvalues. Therefore the array eval must
!!       be two-dimensional (one column for real part and one for the
!!       imaginary part). In the same line: evec must be of size
!!       n * (2*n_ev) because they will be complex (in general)
!!     - The array v can probably be omitted and v be replaced with
!!       evec.
!!
!!   MODIFICATION HISTORY
!!     09-2002, TCH, created
!!
!!   SEE ALSO
!!     arpack_diag_sy
!!
!!****
  subroutine arpack_diag_n(n_ev, eval, evec, glob_env, op_tddfpt)
  
    implicit none

  	! arguments
    integer, intent(in)                          :: n_ev
    real(dbl), dimension(:,:), pointer           :: eval
    type(cp_full_matrix_type), pointer           :: evec
    type(global_environment_type), intent(in)    :: glob_env
    type(qs_p_env_type), intent(inout), optional :: op_tddfpt
  
  	! locals
    integer   :: stat              ! return value of allocation operation
    integer   :: op_present        ! how many optional operators are present
    character(len=*), parameter :: routineN = "arpack_diag_n", &
                                   routineP = moduleN//"/"//routineN
    
    ! arpack specific
    character :: bmat*1, which*2
    logical   :: rvec
    logical, dimension(:), pointer :: ev_select
    integer   :: ido, &
                 n_local, &        ! size of matrix on local node?
                 ldv, &            ! leading dimension ?
                 n_kv, &           ! number of krylov space vectors
                 lworkl, &         ! leading dimension of workl
                 info                 
    real(dbl) :: tol               ! tolerance
    real(dbl) :: sigmar, sigmai    ! unused but must be provided
    integer, dimension(11)             :: param, pntr
    real(dbl), dimension(:), pointer   :: resid, workev, workd, workl, vin, vout
    real(dbl), dimension(:,:), pointer :: v
    
    ! check if an operator has been specified
    op_present = 0
    if (present(op_tddfpt)) then
       op_present = op_present + 1
    end if
    if (op_present /= 1) then
       call stop_program(routineP, "exactly 1 operator must be specified", glob_env)
    end if

  	! code
    call check_input(n_ev, eval, evec, glob_env)
    
    ! set up some variables
    ido   = 0
    call blacs_find_vector_distribution(evec, n_local)
    which = 'SA'
    tol   = 1.0e-15_dbl             ! hardcoded ?
    n_kv  = min(2 * n_ev, n_local)
    if (n_kv < n_ev + 1) &
         call stop_program(routineP, "n_kv < n_ev + 1, n_ev too large", glob_env)
    ldv = n_local
    param = 0
    param(1) = 1    ! exact shift mode
    param(3) = 1500 ! maximum number of iterations
    param(7) = 1    ! simple mode
    pntr     = 0
    lworkl   = n_kv*(n_kv+8)
    info     = 0    ! no starting vector as of now

    ! allocate some arrays
    allocate(resid(n_local), v(ldv, n_kv), &
             workev(3*n_kv), workd(3*n_local), workl(lworkl), &
             ev_select(n_kv), &
             stat=stat)
    if (stat /= 0) &
         call stop_program(routineP, "allocation of temporary variables failed", glob_env)

    reverse_comm: do 
#if defined(__ARPACK)
       call dnaupd(ido, bmat, n_local, which, n_ev, tol, resid, &
                   n_kv, v, ldv, param, pntr, workd, workl, lworkl, info)
       
       if (ido==-1.or.ido==1) then
          vin  => workd(pntr(1):)
          vout => workd(pntr(2):)
          if (present(op_tddfpt)) then
             call arpack_op_tddfpt(op_tddfpt, vin, vout)
          end if
       else
          exit reverse_comm
       end if
       if (info < 0) then
          !! report info
          call stop_program(routineP, "error in dnaupd", glob_env)
       end if
#endif

    end do reverse_comm

    if (info /=0) then
       call stop_program(routineP, "errors in revers communication process", glob_env)
    end if
    
#if defined(__ARPACK)
    rvec = .true.          ! we want the eigenvectors
    ev_select = .true.     ! just to be sure
    sigmar = 0.0_dbl; sigmai = 0.0_dbl
    call dneupd(rvec, 'A', ev_select, eval(:,1), eval(:,2), evec, ldv, &
                sigmar, sigmai, workev, &
                bmat, n_local, which, n_ev, tol, resid, &
                n_kv, v, ldv, param, pntr, workd, workl, lworkl, info)
#endif

    if (info /= 0) then
       !! report info
       call stop_program(routineP, "error in dneupd", glob_env)
    end if

    deallocate(resid, v, workev, workd, workl, ev_select, stat=stat)
    if (stat /= 0) &
         call stop_program(routineP, "deallocation of temporary variables failed", glob_env)

  end subroutine arpack_diag_n

! *****************************************************************************

  subroutine check_input(n_ev, eval, evec, glob_env)
    implicit none

    ! arguments
    integer, intent(in)                       :: n_ev
    real(dbl), dimension(:,:), pointer        :: eval
    type(cp_full_matrix_type), pointer        :: evec
    type(global_environment_type), intent(in) :: glob_env

    ! locals
    integer                     :: evec_cols
    character(len=*), parameter :: routineN = "check_input", &
                                   routineP = moduleN//"/"//routineN

    ! code
    if (n_ev < 1) &
         call stop_program(routineP, "n_ev must be at least 1", glob_env)
    if (size(eval,1) < n_ev) &
         call stop_program(routineP, "size(eval) too small", glob_env)
    call cp_fm_get_info(evec, ncol_global=evec_cols)
    ! not quite sure if it really is necessary to have 2*n_ev cols
    if (evec_cols < 2*n_ev) &
         call stop_program(routineP, "evec needs at least n_ev columns", glob_env)

  end subroutine check_input

! *****************************************************************************

end module qs_arpack_diag_n
