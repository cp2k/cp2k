!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/m *
!!
!!   NAME
!!     m
!!
!!   FUNCTION
!!     
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!
!!****     
module qs_arpack_diag_sy

  use kinds, only: wp => dp, wp_size => dp_size
  use global_types,        only: global_environment_type
  use termination,         only: stop_memory, stop_program
  use timings,             only: timeset, timestop
  use qs_blacs,            only: cp_full_matrix_type,cp_fm_get_info, &
                                 blacs_vector_symv, &
                                 blacs_find_vector_distribution, &
                                 copy_vector_to_blacs, &
                                 copy_blacs_to_vector,blacs_daxpy

  use qs_arpack_full,      only: arpack_full_m_op, full_m_op_type

  implicit none

  character(len=*), parameter :: moduleN = "qs_arpack_diag_sy"

  private
  public :: arpack_diag_sy

contains

! *****************************************************************************

!!****f* cp2k/qs_arpack/arpack_diag_sy *
!!
!!   NAME
!!     arpack_diag_sy
!!
!!   FUNCTION
!!     Given an symmetric operation, computes the lowest
!!     neig eigenvalues and eigenvectors stores these in eval
!!     and vmatrix.
!!     
!!   ARGUMENTS
!!     - vmatrix     : initial guess and eigenvectors output
!!     - eval        : eigenvalues
!!     - nev         : how many eigenvalues to calculate
!!     - ncv         : total number of krylov space vectors
!!     - use_vectors : should initial guess be used?
!!     - globenv     : global environment variables
!!     The following are optional but at least one of them
!!     must be provided. The represent how the object that should be 
!!     diagonalized acts on a vector.
!!     - full_m_op   : represents a simple matrix vector multiply
!!
!!   NOTES
!!     - Uses lanczos with ncv eigenvectors (buffers are allocated
!!     in the routine).
!!     - Setting use_vectors to .true. uses the eigenvectors in vmatrix
!!      as an initial guess
!!     - amatrix has to be upper symmetric and bmatrix should
!!     not be present, otherwise arpack_op_blacs 
!!     has to be modified.
!!     - bmatrix is optional (and not yet implemented to solve for Ax=BxE)
!!
!!   MODIFICATION HISTORY
!!     
!!
!!****
subroutine arpack_diag_sy(vmatrix, &
     eval, nev, ncv, use_vectors,  globenv, &
     full_m_op)

  implicit none

  ! arguments
  type(cp_full_matrix_type), pointer           :: vmatrix
  real(WP), dimension(:), pointer              :: eval
  integer, intent(IN)                          :: nev
  integer, intent(IN)                          :: ncv 
  logical, intent(IN)                          :: use_vectors
  type(global_environment_type), intent(IN)    :: globenv
  type(full_m_op_type), intent(in), optional   :: full_m_op

  ! locals
  logical                           :: op_defined
  integer                           :: context,handle,istat
  real(wp), dimension(:), pointer   :: vin,vout
  character(len=*), parameter       :: routineN = "arpack_diag_sy", &
                                       routineP = moduleN//"/"//routineN
                                    
  ! locals arpack                   
  integer                           :: ploc,nloc,n,ncolsv
  character                         :: bmat*1,which*2
  real(wp), dimension(:), pointer   :: workl,workd,d,resid,ax
  real(wp), dimension(:,:), pointer :: cv
  logical, dimension(:), pointer    :: select
  logical                           :: rvec
  integer                           :: iparam(11), ipntr(11), &
                                       ierr, lworkl, info, ido
  integer                           :: arpack_comm,i
  integer, dimension(9)             :: descv
  real(wp)                          :: tol,sigma,alpha

  call timeset(routineN,'I',' ',handle)

  ! code to check if at least one of the optional operation
  ! types has been specified
  op_defined = .false.
  if (present(full_m_op)) then
     op_defined = .true.
  end if
  if (.not.op_defined) &
       call stop_program(routineP, "An operation must be specified")

  context = globenv%blacs_env%group
  call arpack_confirm(globenv)

  ! get size of problem
  if (present(full_m_op)) then
     n = full_m_op%a_size
  end if

  ! this should go into qs_arpack_full
  call cp_fm_get_info(vmatrix, ncol_global=ncolsv)
  if (ncolsv.lt.2) then
      call stop_program(routineP,'too few cols in (buffer) vmatrix',globenv)
  endif

  if (nev .lt. 1) then
     call stop_program(routineP,"nev too small",globenv)
  endif
  if (ncv .lt. nev+1) then
     call stop_program(routineP,"ncv to small",globenv)
  endif
! limitation of arpack ... no way to get all eigenvectors (well, wouldn't make sence anyway)
  if (ncv .gt. n) then
     call stop_program(routineP,"ncv to large",globenv)
  endif

  ! here we learn how big nloc has to be
  call blacs_find_vector_distribution(vmatrix,nloc)

!
! allocate the scratch space, notice that some of these (cv,workl) are not small..
! if we would force the user to provide sufficient buffer in vmatrix we could
! reuse it, but specifying the right size would be non-trivial
!

  lworkl = ncv*(ncv+8)
  allocate(workl(lworkl),STAT=istat)
  if (istat.ne.0) then
    call stop_memory(routineP,"workl",lworkl*wp_size)
  endif

  allocate(cv(nloc,ncv),STAT=istat)
  if (istat.ne.0) then
    call stop_memory(routineP,"cv",nloc*ncv*wp_size)
  endif

  allocate(select(ncv),STAT=istat)
  if (istat.ne.0) then
    ! size of logical ... ?
    call stop_memory(routineP,"select",ncv*4)
  endif

  allocate(workd(3*nloc),STAT=istat)
  if (istat.ne.0) then
    call stop_memory(routineP,"workd",3*nloc*wp_size)
  endif

  allocate(d(2*ncv),STAT=istat)
  if (istat.ne.0) then
    call stop_memory(routineP,"d",2*ncv*wp_size)
  endif

  allocate(resid(nloc),STAT=istat)
  if (istat.ne.0) then
    call stop_memory(routineP,"resid",nloc*wp_size)
  endif

  allocate(ax(nloc),STAT=istat)
  if (istat.ne.0) then
    call stop_memory(routineP,"ax",nloc*wp_size)
  endif

  arpack_comm = 0

#if defined(__parallel)
  ! arpack_comm = context       ! use BLACS / seems to be a bug on SUN
  arpack_comm = globenv%group   ! use MPI / needs renaming of pxlamch
                                ! conflicts with scalapack
#endif
  
  
  ! arpack settings
  !
  if (present(full_m_op)) then
     if (full_m_op%general) then
        bmat = 'G'   ! solve general problem Ax = BxE
     else
        bmat = 'I'   ! solve simple problem Ax = xE
     end if
  end if
  which = 'SA'       ! compute the nev smallest eigenvalues (algebraic)
  tol  = 1.0E-15_wp  ! accuracy
  info = 0           ! no intial guess for the eigenvectors
                     ! one starting guess can be given in resid if info=1
  ido  = 0           ! initial go
  iparam(:) = 0
  iparam(7) = 1      ! simple mode
  iparam(1) = 1      ! with exact shifts
  iparam(3) = 1500   ! maxiter
  select(:) = .true.


  if (use_vectors) then
    alpha=1.0D0
    ! sum all eigenvectors in the first one
    do i=2, nev ! neig
       call blacs_daxpy(vmatrix,i,alpha,1)
    enddo
    !  and use as an initial state
    info=1
    call copy_blacs_to_vector(vmatrix,1,resid)
  endif 
   
  !  the famous reverse communication let dsaupd 'call' arpack_op_blacs
  do
     !write(6,*) "in the loop",globenv%mepos

#if defined(__ARPACK)
#if defined(__parallel)
     call pdsaupd(arpack_comm, ido, bmat, nloc, which, nev, tol, resid, &
                 ncv, cv, nloc, iparam, ipntr, workd, workl, &
                 lworkl, info )
#else
     call dsaupd ( ido, bmat, nloc, which, nev, tol, resid, &
                 ncv, cv, nloc, iparam, ipntr, workd, workl, &
                 lworkl, info )
#endif
#endif

     if (ido .eq. -1 .or. ido .eq. 1) then
       vin =>workd(ipntr(1):)
       vout=>workd(ipntr(2):)
       call arpack_full_m_op(full_m_op, vin, vout)
     else
       exit
     endif
     if (info.lt.0) then
        call stop_program(routineP,"error in dsaupd",globenv) 
     endif
  enddo

  if (info.lt.0) then
        write(6,*) info
     call stop_program(routineP,"error in dsaupd",globenv) 
  endif


  rvec = .true.
  ierr = 0 

#if defined(__ARPACK)
#if defined(__parallel)
  call pdseupd ( arpack_comm, rvec, 'All', select, d, cv, nloc, sigma, &        
                bmat, nloc, which, nev, tol, resid, ncv, cv, nloc, &
                iparam, ipntr, workd, workl, lworkl, ierr )
#else
  call dseupd ( rvec, 'All', select, d, cv, nloc, sigma, &        
                bmat, nloc, which, nev, tol, resid, ncv, cv, nloc, &
                iparam, ipntr, workd, workl, lworkl, ierr )
#endif
#endif

  if ( ierr .ne. 0) then
     call stop_program(routineP,"dseupd error",globenv)
  endif

  if ( iparam(5) .ne. nev) then
     call stop_program(routineP,"dseupd not all converged",globenv)
  endif

  if ( info .eq. 1) then
     call stop_program(routineP,"max number of iterations reached",globenv)
  endif

  ! copy eigenvalues
  eval(1:nev)=d(1:nev)
  do i=1,nev
     ! is there any clean way to do this ... ?
     workd(1:nloc)=cv(1:nloc,i)
     call copy_vector_to_blacs(workd,vmatrix,i)
  enddo

  deallocate(cv,workl,workd,d,resid,ax,select)
  call timestop(0.0_wp,handle)
  
end subroutine arpack_diag_sy

! *****************************************************************************

! stops the program if arpack is not present
! can also be used to check if the program is present
! in that case use the optional logical
subroutine arpack_confirm(globenv,is_present)

  implicit none

  ! arguments
  type(global_environment_type), intent(IN) :: globenv
  logical, optional, intent(INOUT)          :: is_present


#if defined(__ARPACK)
  if (present(is_present)) then ! hmmmm 
     is_present=.true.
  endif
#else
  if (present(is_present)) then ! hmmmm 
     is_present=.false.
  else
     call stop_program("arpack_confirm",& 
              "(p)arpack library not present but needed",globenv)
  endif
#endif

end subroutine arpack_confirm

! *****************************************************************************

end module

