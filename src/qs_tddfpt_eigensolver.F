!-----------------------------------------------------------------------------!
!   CP2K: A GENERAL PROGRAM TO PERFORM MOLECULAR DYNAMICS SIMULATIONS         !
!   COPYRIGHT (C) 2000  CP2K DEVELOPERS GROUP                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****M* CP2K/qs_tddfpt_eigensolver *
!!
!!   NAME
!!     qs_tddfpt_eigensolver
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     
!!     
!!****
module qs_tddfpt_eigensolver
  use cp_control_types,                only: tddfpt_control_type,&
                                             tddfpt_davidson,&
                                             tddfpt_lanczos
  use cp_error_handling,               only: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  use cp_log_handling,                 only: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  use cp_fm_basic_linalg,              only: cp_fm_scale,&
                                             cp_fm_symm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  use cp_fm_diag,                      only: cp_fm_syevd
  use cp_fm_pool_types,                only: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  use cp_fm_struct,                    only: cp_fm_struct_create,&
                                             cp_fm_struct_release, &
                                             cp_fm_struct_type
  use cp_fm_types,                     only: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_get_element,&
                                             cp_fm_set_element,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  use cp_log_handling,                 only: cp_to_string
  use qs_mo_types,                     only: get_mo_set
  use cp_sm_fm_interactions,           only: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  use kinds,                           only: dp
  use physcon,                         only: evolt
  use qs_environment_types,            only: qs_environment_type,&
                                             get_qs_env
  use qs_matrix_pools,                 only: mpools_get
  use qs_p_env_methods,                only: p_op_l1,&
                                             p_op_l2,&
                                             p_postortho,&
                                             p_preortho
  use qs_p_env_types,                  only: qs_p_env_type
  use qs_tddfpt_types
  use qs_tddfpt_utils,                 only: calc_initial_guess,&
                                             calc_norm,&
                                             davidson_initial_guess, &
                                             normalize, &
                                             reorthogonalize, &
                                             dot_prod
  use sparse_matrix_types,             only: set_matrix, &
                                             get_matrix_info, &
                                             get_matrix_diagonal
  use termination,                     only: stop_program
  use timings,                         only: timeset,&
                                             timestop

  implicit none

  type cp_fm_struct_p_type
     type(cp_fm_struct_type), pointer :: struct
  end type cp_fm_struct_p_type

  character(len=*), parameter :: moduleN = "qs_tddfpt_eigensolver"

  private

  public :: davidson, lanczos

contains

! *****************************************************************************

  subroutine davidson(t_env, p_env, qs_env, error)

    type(tddfpt_env_type), intent(INOUT)     :: t_env
    type(qs_p_env_type), pointer             :: p_env
    type(qs_environment_type), pointer       :: qs_env
    type(cp_error_type), &
         intent(INOUT), optional             :: error

    integer                                  :: ispin, nspins, handle
    integer                                  :: n_ev, i, restarts, out_unit
    logical                                  :: do_kernel_save, converged, ionode
    real(KIND=dp), dimension(:)              :: evals
    type(cp_fm_p_type), dimension(:,:), &
      pointer                                :: evecs, iguess
    type(cp_fm_pool_p_type), dimension(:), &
      pointer                                :: ao_mo_fm_pools
    type(cp_logger_type), pointer            :: logger

    character(len=*), parameter :: routineN = "lanczos", &
                                   routineP = moduleN//"/"//routineN

    call timeset(routineN,"I","",handle)

    nullify(ao_mo_fm_pools,evecs,logger)

    logger    => cp_error_get_logger(error)
    ionode    = (logger%para_env%mepos==logger%para_env%source)
    out_unit  = cp_logger_get_default_unit_nr(logger)

    n_ev      = qs_env%dft_control%tddfpt_control%n_ev
    nspins    = qs_env%dft_control%nspins

    do_kernel_save = qs_env%dft_control%tddfpt_control%do_kernel
    qs_env%dft_control%tddfpt_control%do_kernel = .false.

    call mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    allocate(evecs(n_ev,nspins),iguess(1,nspins),evals(n_ev))
    do ispin=1,nspins
       do i=1, n_ev
          call fm_pool_create_fm(ao_mo_fm_pools(ispin)%pool,evecs(i,ispin)%matrix)
       end do
    end do
    
    do restarts=1, 5
       if (restarts==1) then
          do ispin=1, nspins
             call cp_fm_set_all(evecs(1,ispin)%matrix,1.0_dp)
             iguess(1,ispin)%matrix => evecs(1,ispin)%matrix
          end do
          if (iterative_solver('L', evals, iguess, t_env, p_env, qs_env, evals, evecs, error=error)) exit
       else
          if (iterative_solver('L', evals, evecs, t_env, p_env, qs_env, evals, evecs, error=error)) exit
       end if
    end do 

    qs_env%dft_control%tddfpt_control%do_kernel = do_kernel_save

    do restarts=1, 5
       if (iterative_solver('D', evals, evecs, t_env, p_env, qs_env, evals, evecs, error=error)) exit
    end do

    do ispin=1,nspins
       do i=1, n_ev
          call fm_pool_give_back_fm(ao_mo_fm_pools(ispin)%pool, evecs(i,ispin)%matrix)
       end do
    end do
    deallocate(evecs,iguess,evals)

    call timestop(0.0_dp, handle)

  end subroutine davidson

! *****************************************************************************

  subroutine lanczos(t_env, p_env, qs_env, error)

    type(tddfpt_env_type), intent(INOUT)     :: t_env
    type(qs_p_env_type), pointer             :: p_env
    type(qs_environment_type), pointer       :: qs_env
    type(cp_error_type), &
         intent(INOUT), optional             :: error

    integer                                  :: ispin, nspins, handle
    integer                                  :: n_ev, i, restarts
    logical                                  :: do_kernel_save, converged
    real(KIND=dp), dimension(:)              :: evals
    type(cp_fm_p_type), dimension(:,:), &
      pointer                                :: evecs, iguess
    type(cp_fm_pool_p_type), dimension(:), &
      pointer                                :: ao_mo_fm_pools

    character(len=*), parameter :: routineN = "lanczos", &
                                   routineP = moduleN//"/"//routineN

    call timeset(routineN,"I","",handle)

    nullify(ao_mo_fm_pools,evecs)

    n_ev      = qs_env%dft_control%tddfpt_control%n_ev
    nspins    = qs_env%dft_control%nspins

    do_kernel_save = qs_env%dft_control%tddfpt_control%do_kernel
    qs_env%dft_control%tddfpt_control%do_kernel = .false.

    call mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    allocate(evecs(n_ev,nspins),evals(n_ev),iguess(1,nspins))
    do ispin=1,nspins
       do i=1, n_ev
          call fm_pool_create_fm(ao_mo_fm_pools(ispin)%pool,evecs(i,ispin)%matrix)
       end do
    end do
    
    do restarts=1, 5

       if (restarts==1) then
          do ispin=1, nspins
             call cp_fm_set_all(evecs(1,ispin)%matrix,1.0_dp)
             iguess(1,ispin)%matrix => evecs(1,ispin)%matrix
          end do
          if (iterative_solver('L', evals, iguess, t_env, p_env, qs_env, evals, evecs, error=error)) exit
       else
          if (iterative_solver('L', evals, evecs, t_env, p_env, qs_env, evals, evecs, error=error)) exit
       end if

    end do 

    qs_env%dft_control%tddfpt_control%do_kernel = do_kernel_save

    do i=1, 5

       if (iterative_solver('L', evals, evecs, t_env, p_env, qs_env, out_evecs=evecs, error=error)) exit

    end do
    
    do ispin=1,nspins
       do i=1, n_ev
          call fm_pool_give_back_fm(ao_mo_fm_pools(ispin)%pool, evecs(i,ispin)%matrix)
       end do
    end do
    deallocate(evecs,evals,iguess)

    call timestop(0.0_dp, handle)

  end subroutine lanczos

! *****************************************************************************

  ! mode     : 'L' = lanczos, 'D' = davidson
  !
  ! in_evals : intial guess eigenvalues  (ignored in lanczos)
  ! in_evecs : intial guess eigenvectors 
  !            - lanczos: can be any number up to n_ev
  !            - davidson: must be n_ev
  ! iev      : which eigenvalue to converge to (ignored in lanczos)
  !
  ! t_env, p_env, qs_env : environments needed for calculation
  !
  ! out_evals : resulting eigenvalues
  ! out_evecs : resulting eigenvectors

  function iterative_solver(mode, &
                            in_evals, in_evecs, &
                            t_env, p_env, qs_env, &
                            out_evals, out_evecs, error) result(res)

    character, intent(in)                    :: mode
    real(KIND=dp), dimension(:), intent(IN)  :: in_evals
    type(cp_fm_p_type), dimension(:,:), &
      pointer                                :: in_evecs
    type(tddfpt_env_type), intent(INOUT)     :: t_env
    type(qs_p_env_type), pointer             :: p_env
    type(qs_environment_type), pointer       :: qs_env
    real(kind=dp), dimension(:), &
         intent(OUT), optional               :: out_evals
    type(cp_fm_p_type), dimension(:,:), &
      pointer, optional                      :: out_evecs
    type(cp_error_type), &
         intent(INOUT), optional             :: error
    logical                                  :: res

    character(len=*), parameter :: routineN = "iterative_solver", &
                                   routineP = moduleN//"/"//routineN

    integer :: col, handle, i, iter, j, max_iter, &
      max_krylovspace_dim, max_kv, n_ev, n_kv, nspins, &
      row, spin, stat, k
    real(dp)                                 :: Atilde_ij, tmp, convergence
    real(kind=dp), allocatable, &
      dimension(:, :)                        :: evals
    type(cp_fm_p_type), dimension(:), &
      pointer                                :: R, X
    type(cp_fm_p_type), dimension(:, :), &
      pointer                                :: b, Ab, Sb
    type(cp_fm_pool_p_type), dimension(:), &
      pointer                                :: ao_mo_fm_pools
    type(cp_fm_struct_p_type), &
      dimension(:), pointer                  :: kv_fm_struct
    type(cp_fm_p_type)                       :: Atilde, Us
    type(cp_fm_struct_type), pointer         :: tilde_fm_struct, Us_fm_struct
    type(tddfpt_control_type), pointer       :: tddfpt_control

    type(cp_logger_type), pointer            :: logger
    logical                                  :: ionode, failure
    integer                                  :: out_unit

    real(KIND=dp) :: dot

    failure   = .false.
    res       = .false.

    call timeset(routineN,"I","",handle)

    nullify(ao_mo_fm_pools, X, R, b, Ab, Sb, tddfpt_control, logger)

    call get_qs_env(qs_env, tddfpt_control=tddfpt_control)
    logger    => cp_error_get_logger(error)
    ionode    = (logger%para_env%mepos==logger%para_env%source)
    out_unit  = cp_logger_get_default_unit_nr(logger)
    n_ev      = tddfpt_control%n_ev
    nspins    = qs_env%dft_control%nspins

    !-----------------------------------------!
    ! determine the size of the problem       !
    ! and how many krylov space vetors to use ! 
    !-----------------------------------------!
    max_krylovspace_dim = sum(p_env%n_ao(1:nspins)*p_env%n_mo(1:nspins))
    max_kv = tddfpt_control%max_kv
    CPPostcondition(max_krylovspace_dim>=max_kv,cp_failure_level,routineP,error,failure)
    if (mode == 'D') then
       max_iter = (max_kv / n_ev) 
    else if (mode == 'L') then
       max_iter = max_kv - n_ev + 1
    end if

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    call mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    call fm_pools_create_fm_vect(ao_mo_fm_pools, X, name=routineP//":X")
    call fm_pools_create_fm_vect(ao_mo_fm_pools, R, name=routineP//":R")

    allocate(evals(max_kv, 0:max_iter), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    allocate(b(max_kv,nspins), Ab(max_kv,nspins), &
             Sb(max_kv,nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    allocate(kv_fm_struct(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    do spin=1, nspins
       nullify (kv_fm_struct(spin)%struct)
       call cp_fm_struct_create(kv_fm_struct(spin)%struct, qs_env%para_env, &
                                qs_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_mo(spin))
    end do
    do spin=1, nspins
       do i=1, max_kv
          nullify(b(i,spin)%matrix, Ab(i,spin)%matrix, Sb(i,spin)%matrix)
       end do
    end do

    if (ionode) then
       write (out_unit, '(2X,A,T69,A)') &
              "nvec", "Convergence"
       write (out_unit, '(2X,A)') &
              "-----------------------------------------------------------------------------"
    end if

    iter = 1
    k    = 0
    iteration: do
       
       n_kv = n_ev
       if (mode=='L' .and. iter/=1) n_kv = 1
              
       call allocate_krylov_vectors(b,  "b-",  k+1, n_kv, nspins, kv_fm_struct)
       call allocate_krylov_vectors(Ab, "Ab-", k+1, n_kv, nspins, kv_fm_struct)
       call allocate_krylov_vectors(Sb, "Sb-", k+1, n_kv, nspins, kv_fm_struct)

       do i=1, n_kv
          k=k+1

          if (k <= size(in_evecs,1) ) then

             ! take the initial guess
             do spin=1, nspins
                call cp_fm_to_fm(in_evecs(k,spin)%matrix, b(k,spin)%matrix)
             end do

          else

             ! create a new vector
             if (mode == 'L') then
                do spin=1, nspins
                   if (tddfpt_control%invert_S) then
                      call cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                           1.0_dp, t_env%invS(spin)%matrix, Ab(k-1,spin)%matrix, &
                           0.0_dp, b(k,spin)%matrix)
                   else
                      call cp_fm_to_fm(Ab(k-1,spin)%matrix, b(k,spin)%matrix)
                   end if
                end do
             else if (mode == 'D') then

                ! create the new davidson vector
                do spin=1, nspins

                   call cp_fm_set_all(X(spin)%matrix, 0.0_dp)
                   do j=1, k-i
                      call cp_fm_get_element(Us%matrix, j, i, tmp)
                      call cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                               tmp, Ab(j,spin)%matrix)
                      call cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                               -evals(i,iter-1), Sb(j,spin)%matrix)
                   end do
                   do row=1, p_env%n_ao(spin)
                      do col=1, p_env%n_mo(spin)
                         call cp_fm_get_element(X(spin)%matrix, row, col, tmp)
                         if (col<=n_ev) then
                            tmp = tmp / (evals(i,iter-1) - in_evals(col))
                         else
                            tmp = tmp / (evals(i,iter-1) - (in_evals(n_ev)+1.0_dp))
                         end if
                         call cp_fm_set_element(X(spin)%matrix, row, col, tmp)
                      end do
                   end do

                   if (tddfpt_control%invert_S) then
                      call cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                                      1.0_dp, t_env%invS(spin)%matrix, X(spin)%matrix, &
                                      0.0_dp, b(k,spin)%matrix)
                   else
                      call cp_fm_to_fm(X(spin)%matrix, b(k,spin)%matrix)
                   end if
                   
                end do

             else
                if (ionode) write (*,*) "unknown mode"
                CPPostcondition(.false.,cp_failure_level,routineP,error,failure)
             end if

          end if
          
          call p_preortho(p_env, qs_env, b(k,:))
          do j=1, tddfpt_control%n_reortho
             call reorthogonalize(b(k,:), b, Sb, R, k-1) ! R is temp
          end do
          call normalize(b(k,:), R, qs_env%matrix_s) ! R is temp
          do spin=1, nspins
             call cp_fm_to_fm(b(k,spin)%matrix, X(spin)%matrix)
          end do
          call apply_op(X, Ab(k,:), t_env, p_env, qs_env)
          do spin=1, nspins
             call cp_sm_fm_multiply(qs_env%matrix_s(spin)%matrix, &
                  b(k,spin)%matrix, &
                  Sb(k,spin)%matrix, &
                  p_env%n_mo(spin))
          end do
       end do
       
       !--------------------------------------------!
       ! deallocate memory for the reduced matrices !
       !--------------------------------------------!
       if (associated(Atilde%matrix)) call cp_fm_release(Atilde%matrix)
       if (associated(Us%matrix)) call cp_fm_release(Us%matrix)
       if (associated(tilde_fm_struct)) call cp_fm_struct_release(tilde_fm_struct)

       !------------------------------------------!
       ! allocate memory for the reduced matrices !
       !------------------------------------------!
       call cp_fm_struct_create(tilde_fm_struct, &
            qs_env%para_env, &
            qs_env%blacs_env, &
            k, k)
       call cp_fm_create(Atilde%matrix, &
            tilde_fm_struct, &
            routineP//"Atilde")
       call cp_fm_create(Us%matrix, &
            tilde_fm_struct, &
            routineP//"Us")

       !---------------------------------------!
       ! calc the matrix Atilde = transp(b)*Ab !
       !---------------------------------------!
       do i=1, k
          do j=1, k
             Atilde_ij= 0.0_dp
             do spin=1, nspins
                call cp_fm_trace(b(i,spin)%matrix, Ab(j,spin)%matrix, tmp)
                Atilde_ij = Atilde_ij + tmp
             end do
             call cp_fm_set_element(Atilde%matrix, i, j, Atilde_ij)
          end do
       end do

       !--------------------!
       ! diagonalize Atilde !
       !--------------------!
       call cp_fm_syevd(Atilde%matrix, Us%matrix, evals(:,iter))
!!       do i=1, n_ev
!!          if (ionode) write (out_unit,*) evals(i,iter)
!!       end do

       !-------------------!
       ! check convergence !
       !-------------------!
       if (iter/=1) then

          convergence = maxval(abs((evals(1:n_ev,iter)-evals(1:n_ev,iter-1))/evals(1:n_ev,iter-1)))

          if (ionode) write (out_unit, '(2X,I4,T69,ES11.4)') k, convergence

          if (convergence < tddfpt_control%tolerance) then
             res = .true.
             exit iteration
          end if
          if (iter == max_iter) exit iteration

       end if

       iter = iter + 1

    end do iteration

    if (present(out_evals)) then
       out_evals(1:n_ev) = evals(1:n_ev,iter)
    end if

    if (present(out_evecs)) then
       do spin=1, nspins
          do j=1, n_ev
             call cp_fm_set_all(out_evecs(j,spin)%matrix, 0.0_dp)
             do i=1, k
                call cp_fm_get_element(Us%matrix, i, j, tmp)
                call cp_fm_scale_and_add(1.0_dp, out_evecs(j,spin)%matrix, &
                     tmp, b(i,spin)%matrix)
             end do
          end do
       end do
    end if

    do spin=1, nspins
       do i=1, max_kv
          if (associated(b(i,spin)%matrix)) &
               call cp_fm_release(b(i,spin)%matrix)
          if (associated(Ab(i,spin)%matrix)) &
               call cp_fm_release(Ab(i,spin)%matrix)
          if (associated(Sb(i,spin)%matrix)) &
               call cp_fm_release(Sb(i,spin)%matrix)
       end do
    end do

    !-------------------!
    ! Print the results !
    !-------------------!
    if (ionode) then
       write (out_unit, *)
       write (out_unit, *)
       write (out_unit, '(1X,A,T19,A,T45,A)') &
            "Nr.", "Eigenvalue (hartree)", "Eigenvalue (eV)"
       do i=1, min(n_ev,max_kv)
          write (out_unit, '(1X,I3,T20,F19.10,2X,F19.10)') &
               i, evals(i, iter), evals(i, iter)*evolt
       end do
       write (out_unit, *)
    end if

    !----------!
    ! clean up !
    !----------!
    if (associated(Atilde%matrix)) call cp_fm_release(Atilde%matrix)
    if (associated(Us%matrix)) call cp_fm_release(Us%matrix)
    if (associated(tilde_fm_struct)) call cp_fm_struct_release(tilde_fm_struct)
    call fm_pools_give_back_fm_vect(ao_mo_fm_pools, X)
    call fm_pools_give_back_fm_vect(ao_mo_fm_pools, R)
    do spin=1, nspins
       call cp_fm_struct_release(kv_fm_struct(spin)%struct)
       do i=1, max_kv
          if (associated(b(i,spin)%matrix)) &
               call cp_fm_release(b(i,spin)%matrix)
          if (associated(Ab(i,spin)%matrix)) &
               call cp_fm_release(Ab(i,spin)%matrix)
          if (associated(Sb(i,spin)%matrix)) &
               call cp_fm_release(Sb(i,spin)%matrix)
       end do
    end do
    deallocate(b, Ab, Sb, evals, kv_fm_struct, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    nullify(X, R, b, Ab, Sb, kv_fm_struct)

    call timestop(0.0_dp, handle)

  end function iterative_solver

! *****************************************************************************

  ! X        : the vector on which to apply the op
  ! R        : the result
  ! t_env    : td-dft environment (mainly control information)
  ! p_env    : perturbation environment (variables)
  !            both of these carry info for the tddfpt calculation
  ! qs_env   : info about a quickstep ground state calculation

  subroutine apply_op(X, R, t_env, p_env, qs_env)


    type(cp_fm_p_type), dimension(:)         :: X, R
    type(tddfpt_env_type), intent(INOUT)     :: t_env
    type(qs_p_env_type), pointer             :: p_env
    type(qs_environment_type), pointer       :: qs_env

    character(LEN=*), parameter :: routineN = "apply_op", &
                                   routineP = moduleN//"/"//routineN

    integer                                  :: handle, nspins, spin
    integer, save                            :: counter = 0

    call timeset(routineN,"I","",handle)

    counter = counter + 1
    
    nspins = qs_env%dft_control%nspins

    !------------!
    ! R = HX-SXL !
    !------------!
    call p_op_l1(p_env, qs_env, X, R)  ! acts on both spins, result in R

    !-----------------!
    ! calc P1 and     !
    ! R = R + K(P1)*C !
    !-----------------!
    if (qs_env%dft_control%tddfpt_control%do_kernel) then
       do spin=1, nspins
          call set_matrix(p_env%p1(spin)%matrix, value=0.0_dp) ! optimize?
          call cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=p_env%psi0d(spin)%matrix,&
                                  matrix_g=X(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dp) 
          call cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=X(spin)%matrix,&
                                  matrix_g=p_env%psi0d(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dp) 
       end do
       call p_op_l2(p_env, qs_env, p_env%p1, X, &
                    alpha=1.0_dp, beta=0.0_dp) ! X = beta*X + alpha*K(P1)*C
       do spin=1, nspins
          call cp_fm_scale_and_add(1.0_dp, R(spin)%matrix, &
                                   1.0_dp, X(spin)%matrix) ! add X to R
       end do
    end if
    call p_postortho(p_env, qs_env, R)
       
    call timestop(0.0_dp, handle)

  end subroutine apply_op

! *****************************************************************************

  subroutine allocate_krylov_vectors(vectors, vectors_name, &
                                     startv, n_v, nspins, fm_struct)

    type(cp_fm_p_type), dimension(:,:), &
      pointer                                :: vectors
    character(LEN=*), intent(IN)             :: vectors_name
    integer, intent(IN)                      :: startv, n_v, nspins
    type(cp_fm_struct_p_type), &
      dimension(:), pointer                  :: fm_struct

    character(LEN=*), parameter :: routineN = "allocate_krylov_vectors", &
                                   routineP = moduleN//"/"//routineN

    character(LEN=80)                        :: mat_name
    integer                                  :: i, index, spin

    do spin=1, nspins
       do index=startv, startv + n_v - 1
          nullify(vectors(index,spin)%matrix)
          mat_name = routineP//vectors_name//trim(cp_to_string(index))&
                 //","//trim(cp_to_string(spin))
          call cp_fm_create(vectors(index,spin)%matrix,  &
               fm_struct(spin)%struct, mat_name)
          if (.not.associated(vectors(index,spin)%matrix)) &
               call stop_program(routineP, &
               "Could not allocate "//trim(mat_name)//".")
       end do
    end do

  end subroutine allocate_krylov_vectors

! *****************************************************************************

end module qs_tddfpt_eigensolver
