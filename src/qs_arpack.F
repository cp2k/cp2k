!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_arpack [1.0] *
!!
!!   NAME
!!     qs_arpack
!!
!!   FUNCTION
!!     provide an interface to (p)arpack
!!
!!   AUTHOR
!!     Joost VandeVondele (05.2002)
!!
!!   MODIFICATION HISTORY
!!     None
!!
!!   SOURCE
!******************************************************************************


! given the symmetric matrix amatrix,
! computes the lowest neig eigenvalues and eigenvectors
! stores these in eval,vmatrix
! uses lanczos with ncv eigenvectors (buffers allocated in the routine)
! setting use_vectors to .true. uses the eigenvectors in vmatrix as an initial guess
! bmatrix is optional (and not yet implemented to solve for Ax=BxE)
! seems to be efficient for neig << n as the work is +- neig^2*n
! the prefactor is such that, it seems not to be highly efficient for 
! finding all the occupied orbitals if a modest (normal) basis set is used.
! with small modifications should also be able to compute eigenvectors/eigenvalues
! of different parts of the spectrum / without cholesky decomposition
!
! access of the amatrix is only through multiplications of a vector with amatrix
! full arpack info/code/documentation can be found as
!
! http://www.netlib.org/arpack/
! http://www.caam.rice.edu/software/ARPACK/index.html
! 
!
! this module needs the MPI version of PARPACK to linked 
! notice that the routines PXLAMCH conflict with the scalapack routines, and have to be
! renamed througout the (P)ARPACK code (see the tools directory ....)
! the same is true for the version of SECOND used
! hence any downloaded binary of (P)ARPACK is not likely to work
!
! notice that the code will only work if __ARPACK is defined
! otherwise produces a more or less empty body
! in this way, people can compile more easily the code for other parts of cp2k 
! without going through the effort of compiling a library for themselves
!
! the idea is that this kind of code should be easy to interface to several matrix formats
! since in principle only the vector matrix multiplication is needed
! found below is the blacs matrix interface

MODULE qs_arpack

! *****************************************************************************
  USE kinds, ONLY: wp => dp, wp_size => dp_size
  USE global_types,        ONLY: global_environment_type
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE qs_blacs,            ONLY: blacs_matrix_type,get_blacs_matrix_info, &
                                 blacs_vector_symv,blacs_find_vector_distribution, &
                                 copy_vector_to_blacs,copy_blacs_to_vector,blacs_daxpy

  IMPLICIT NONE

  PRIVATE :: arpack_op
  PUBLIC  :: arpack_diag,arpack_confirm

  PRIVATE :: arpack_diag_blacs,arpack_op_blacs

  interface arpack_diag
     module procedure arpack_diag_blacs
  end interface

  interface arpack_op
     module procedure arpack_op_blacs
  end interface

CONTAINS

! stops the program if arpack is not present
! can also be used to check if the program is present
! in that case use the optional logical
SUBROUTINE arpack_confirm(globenv,is_present)
  TYPE(global_environment_type), INTENT(IN) :: globenv
  LOGICAL, OPTIONAL, INTENT(INOUT) :: is_present
#if defined(__ARPACK)
  if (PRESENT(is_present)) then ! hmmmm 
     is_present=.true.
  endif
#else
  if (PRESENT(is_present)) then ! hmmmm 
     is_present=.false.
  else
     CALL stop_program("arpack_confirm",& 
              "(p)arpack library not present but needed",globenv)
  endif
#endif

END SUBROUTINE

!
! uses blacs matrices. Notice that amatrix has to be upper symmetric / bmatrix should not be present
! otherwise arpack_op_blacs has to be modified
! vmatrix should contain at least two column
! 
SUBROUTINE arpack_diag_blacs(amatrix,vmatrix,eval,neig,ncv,use_vectors,context,globenv,bmatrix)

  ! arguments
  TYPE(blacs_matrix_type), POINTER           :: amatrix
  TYPE(blacs_matrix_type), POINTER           :: vmatrix
  REAL(WP), DIMENSION(:), POINTER            :: eval
  INTEGER, INTENT(IN)                        :: neig,context
  TYPE(blacs_matrix_type), POINTER, OPTIONAL :: bmatrix
  TYPE(global_environment_type), INTENT(IN) :: globenv
  LOGICAL, INTENT(IN)                        :: use_vectors

  ! locals
  INTEGER :: handle,istat
  REAL(wp), DIMENSION(:), POINTER :: vin,vout

  ! locals arpack
  integer :: nev,ploc,nloc,n,ncolsv
  character :: bmat*1,which*2
  REAL(wp), DIMENSION(:), POINTER :: workl,workd,d,resid,ax
  REAL(wp), DIMENSION(:,:), POINTER :: cv
  logical, DIMENSION(:), POINTER :: select
  logical :: rvec
  integer :: iparam(11),ipntr(11),ierr,lworkl,info,ido,ncv
  integer :: arpack_comm,i
  INTEGER, DIMENSION(9) :: descv
  real(wp) :: tol,sigma,alpha

  CALL timeset('arpack_diag_blacs','I',' ',handle)


  CALL arpack_confirm(globenv)
  CALL get_blacs_matrix_info(amatrix,nrow_global=n)
  CALL get_blacs_matrix_info(vmatrix,ncol_global=ncolsv)
  if (ncolsv.lt.2) then
      CALL stop_program('arpack_diag_blacs','too few cols in (buffer) vmatrix',globenv)
  endif
  nev =  neig
  if (nev .lt. 1) then
     CALL stop_program("arpack_diag","nev too small",globenv)
  endif
  if (ncv .lt. nev+1) then
     CALL stop_program("arpack_diag","ncv to small",globenv)
  endif
! limitation of arpack ... no way to get all eigenvectors (well, wouldn't make sence anyway)
  if (ncv .gt. n) then
     CALL stop_program("arpack_diag","ncv to large",globenv)
  endif

  ! here we learn how big nloc has to be
  CALL blacs_find_vector_distribution(vmatrix,nloc,context,globenv)

!
! allocate the scratch space, notice that some of these (cv,workl) are not small...
! if we would force the user to provide sufficient buffer in vmatrix we could
! reuse it, but specifying the right size would be non-trivial
!

  lworkl = ncv*(ncv+8)
  allocate(workl(lworkl),STAT=istat)
  if (istat.ne.0) then
    CALL stop_memory("arpack_diag","workl",lworkl*wp_size)
  endif

  allocate(cv(nloc,ncv),STAT=istat)
  if (istat.ne.0) then
    CALL stop_memory("arpack_diag","cv",nloc*ncv*wp_size)
  endif

  allocate(select(ncv),STAT=istat)
  if (istat.ne.0) then
    ! size of logical ... ?
    CALL stop_memory("arpack_diag","select",ncv*4)
  endif

  allocate(workd(3*nloc),STAT=istat)
  if (istat.ne.0) then
    CALL stop_memory("arpack_diag","workd",3*nloc*wp_size)
  endif

  allocate(d(2*ncv),STAT=istat)
  if (istat.ne.0) then
    CALL stop_memory("arpack_diag","d",2*ncv*wp_size)
  endif

  allocate(resid(nloc),STAT=istat)
  if (istat.ne.0) then
    CALL stop_memory("arpack_diag","resid",nloc*wp_size)
  endif

  allocate(ax(nloc),STAT=istat)
  if (istat.ne.0) then
    CALL stop_memory("arpack_diag","ax",nloc*wp_size)
  endif

  arpack_comm = 0

#if defined(__parallel)
  ! arpack_comm = context       ! to use BLACS / seems to be a bug on SUN
  arpack_comm = globenv%group ! to use MPI / needs renaming of pxlamch
                              ! conflicts with scalapack
#endif
  
  ! some standard arpack settings
  ! solves the eval problem
  ! AX=XE
  ! for the lowest eigenvalues (E) and eigenvectors
  ! in a kind of default mode
  ! many options can be given to arpack to satisfy all your dreams (about this kind of subject, OK ?)
  bmat = 'I'
  which = 'SA'
  tol  = 1.0E-15_wp
  info = 0 ! no intial guess for the eigenvectors
           ! one starting guess can be given in resid if info=1
  ido  = 0 ! initial go
  iparam(:) = 0
  iparam(7) = 1   ! simple mode
  iparam(1) = 1   ! with exact shifts
  iparam(3) = 1500 ! maxiter
  select(:) = .true.


  if (use_vectors) then
    alpha=1.0D0
    ! sum all eigenvectors in the first one
    DO i=2,neig
    CALL blacs_daxpy(vmatrix,i,alpha,1,context,globenv)
    ENDDO
    !  and use as an initial state
    info=1
    CALL copy_blacs_to_vector(vmatrix,1,resid,context,globenv)
  ENDIF 
   
  !  the famous reverse communication let dsaupd 'call' arpack_op_blacs
  DO
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
       ! vmatrix is buff here
       IF (present(bmatrix)) THEN
         call arpack_op(amatrix,vin,vout,vmatrix,context,globenv,bmatrix)
       ELSE
         call arpack_op(amatrix,vin,vout,vmatrix,context,globenv)
       ENDIF
     else
	exit
     endif
     if (info.lt.0) then
        call stop_program("arpack_diag","error in dsaupd",globenv) 
     endif
  ENDDO

  if (info.lt.0) then
        write(6,*) info
     call stop_program("arpack_diag","error in dsaupd",globenv) 
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
     CALL stop_program("arpack_diag","dseupd error",globenv)
  endif

  if ( iparam(5) .ne. neig) then
     CALL stop_program("arpack_diag","dseupd not all converged",globenv)
  endif

  if ( info .eq. 1) then
     CALL stop_program("arpack_diag","max number of iterations reached",globenv)
  endif

  ! copy eigenvalues
  eval(1:neig)=d(1:neig)
  do i=1,neig
     ! is there any clean way to do this ... ?
     workd(1:nloc)=cv(1:nloc,i)
     call copy_vector_to_blacs(workd,vmatrix,i,context,globenv)
  enddo

  deallocate(cv,workl,workd,d,resid,ax,select)
  CALL timestop(0.0_wp,handle)
END SUBROUTINE arpack_diag_blacs

!
! notice that bmatrix can be used to solve the general eigenvalue problem.
! bmatrix has to be the cholesky decomp.
! L^X A L^X  should be calculate as described in the ARPACK manual
! right now, it is not implemented
!
SUBROUTINE arpack_op_blacs(amatrix,vin,vout,matrix_buf,context,globenv,bmatrix)

  ! arguments
  TYPE(blacs_matrix_type), POINTER           :: amatrix,matrix_buf
  TYPE(blacs_matrix_type), POINTER, OPTIONAL :: bmatrix
  REAL(wp), DIMENSION(:), POINTER            :: vin,vout
  integer, intent(IN)                        :: context
  TYPE(global_environment_type), INTENT(IN) :: globenv

  IF (PRESENT(bmatrix)) THEN
     call stop_program("arpack_op_blacs","bmatrix not yet implemented",globenv)
  ENDIF

  CALL  blacs_vector_symv(amatrix,vin,vout,matrix_buf,context,globenv)

END SUBROUTINE arpack_op_blacs


END MODULE qs_arpack
