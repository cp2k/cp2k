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
MODULE qs_arpack_diag_sy

  USE kinds, ONLY: wp => dp, wp_size => dp_size
  USE global_types,        ONLY: global_environment_type
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE qs_arpack_utils,     ONLY: arpack_confirm
  USE qs_blacs,            ONLY: cp_full_matrix_type,cp_fm_get_info, &
                                 blacs_vector_symv, &
                                 blacs_find_vector_distribution, &
                                 copy_vector_to_blacs, &
                                 copy_blacs_to_vector,blacs_daxpy

  ! the operators
  USE qs_arpack_full,      ONLY: arpack_op_symv, op_symv_type

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_arpack_diag_sy"

  PRIVATE
  PUBLIC :: arpack_diag_sy

CONTAINS

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
!!     - op_symv   : represents a simple matrix vector multiply
!!
!!   NOTES
!!     - Uses lanczos with ncv krylov vectors (buffers are allocated
!!       in the routine).
!!     - Setting use_vectors to .true. uses the eigenvectors in vmatrix
!!       as an initial guess
!!     - amatrix has to be upper symmetric and bmatrix should
!!       not be present, otherwise arpack_op_blacs 
!!       has to be modified.
!!     - bmatrix is optional (and not yet implemented to solve for Ax=BxE)
!!
!!   MODIFICATION HISTORY
!!     
!!
!!****
SUBROUTINE arpack_diag_sy(vmatrix, &
                          eval, nev, ncv, use_vectors, globenv, &
                          op_symv)

  IMPLICIT NONE

  ! arguments
  TYPE(cp_full_matrix_type), POINTER           :: vmatrix
  REAL(WP), DIMENSION(:), POINTER              :: eval
  INTEGER, INTENT(IN)                          :: nev
  INTEGER, INTENT(IN)                          :: ncv 
  LOGICAL, INTENT(IN)                          :: use_vectors
  TYPE(global_environment_type), INTENT(IN)    :: globenv
  TYPE(op_symv_type), INTENT(in), OPTIONAL     :: op_symv

  ! locals
  INTEGER                           :: op_present
  INTEGER                           :: context,handle,istat
  REAL(wp), DIMENSION(:), POINTER   :: vin,vout
  CHARACTER(len=*), PARAMETER       :: routineN = "arpack_diag_sy", &
                                       routineP = moduleN//"/"//routineN
                                    
  ! locals arpack                   
  INTEGER                           :: ploc, nloc, n, vmatrix_cols
  CHARACTER                         :: bmat*1,which*2
  REAL(wp), DIMENSION(:), POINTER   :: workl,workd,d,resid,ax
  REAL(wp), DIMENSION(:,:), POINTER :: cv
  LOGICAL, DIMENSION(:), POINTER    :: SELECT
  LOGICAL                           :: rvec
  INTEGER                           :: iparam(11), ipntr(11), &
                                       ierr, lworkl, info, ido
  INTEGER                           :: arpack_comm,i
  INTEGER, DIMENSION(9)             :: descv
  REAL(wp)                          :: tol,sigma,alpha

  CALL timeset(routineN,'I',' ',handle)

  ! code to check if at least one of the optional operation
  ! types has been specified
  op_present = 0
  IF (PRESENT(op_symv)) THEN
     n = op_symv%a_size
     op_present = op_present + 1
  END IF
  IF (op_present /= 1) THEN
     CALL stop_program(routineP, "An operation must be specified", globenv)
  END IF

  context = globenv%blacs_env%group
  CALL arpack_confirm(globenv)

  CALL cp_fm_get_info(vmatrix, ncol_global=vmatrix_cols)
  IF (vmatrix_cols < nev) THEN
     CALL stop_program(routineP, "vmatrix should have at least nev columns", globenv)
  END IF

  IF (nev .LT. 1) THEN
     CALL stop_program(routineP,"nev too small",globenv)
  ENDIF
  IF (ncv .LT. nev+1) THEN
     CALL stop_program(routineP,"ncv to small",globenv)
  ENDIF
! limitation of arpack ... no way to get all eigenvectors (well, wouldn't make sence anyway)
  IF (ncv .GT. n) THEN
     CALL stop_program(routineP,"ncv to large",globenv)
  ENDIF

  ! here we learn how big nloc has to be
  CALL blacs_find_vector_distribution(vmatrix,nloc)

!
! allocate the scratch space, notice that some of these (cv,workl) are not small..
! if we would force the user to provide sufficient buffer in vmatrix we could
! reuse it, but specifying the right size would be non-trivial
!

  lworkl = ncv*(ncv+8)
  ALLOCATE(workl(lworkl),STAT=istat)
  IF (istat.NE.0) THEN
    CALL stop_memory(routineP,"workl",lworkl*wp_size)
  ENDIF

  ALLOCATE(cv(nloc,ncv),STAT=istat)
  IF (istat.NE.0) THEN
    CALL stop_memory(routineP,"cv",nloc*ncv*wp_size)
  ENDIF

  ALLOCATE(SELECT(ncv),STAT=istat)
  IF (istat.NE.0) THEN
    ! size of logical ... ?
    CALL stop_memory(routineP,"select",ncv*4)
  ENDIF

  ALLOCATE(workd(3*nloc),STAT=istat)
  IF (istat.NE.0) THEN
    CALL stop_memory(routineP,"workd",3*nloc*wp_size)
  ENDIF

  ALLOCATE(d(2*ncv),STAT=istat)
  IF (istat.NE.0) THEN
    CALL stop_memory(routineP,"d",2*ncv*wp_size)
  ENDIF

  ALLOCATE(resid(nloc),STAT=istat)
  IF (istat.NE.0) THEN
    CALL stop_memory(routineP,"resid",nloc*wp_size)
  ENDIF

  ALLOCATE(ax(nloc),STAT=istat)
  IF (istat.NE.0) THEN
    CALL stop_memory(routineP,"ax",nloc*wp_size)
  ENDIF

  arpack_comm = 0

#if defined(__parallel)
  ! arpack_comm = context       ! use BLACS / seems to be a bug on SUN
  arpack_comm = globenv%group   ! use MPI / needs renaming of pxlamch
                                ! conflicts with scalapack
#endif
  
  
  ! arpack settings
  !
  IF (PRESENT(op_symv)) THEN
     IF (op_symv%general) THEN
        bmat = 'G'   ! solve general problem Ax = BxE
        CALL stop_program(routineP, "General problem not yet implemented")
     ELSE
        bmat = 'I'   ! solve simple problem Ax = xE
     END IF
  END IF
  which = 'SA'       ! compute the nev smallest eigenvalues (algebraic)
  tol  = 1.0E-15_wp  ! accuracy
  info = 0           ! no intial guess for the eigenvectors
                     ! one starting guess can be given in resid if info=1
  ido  = 0           ! initial go
  iparam(:) = 0
  iparam(7) = 1      ! simple mode
  iparam(1) = 1      ! with exact shifts
  iparam(3) = 1500   ! maxiter
  SELECT(:) = .TRUE.


  IF (use_vectors) THEN
    alpha=1.0D0
    ! sum all eigenvectors in the first one
    DO i=2, nev ! neig
       CALL blacs_daxpy(vmatrix,i,alpha,1)
    ENDDO
    !  and use as an initial state
    info=1
    CALL copy_blacs_to_vector(vmatrix,1,resid)
  ENDIF 
   
  !  the famous reverse communication let dsaupd 'call' arpack_op_blacs
  DO
     !write(6,*) "in the loop",globenv%mepos

#if defined(__ARPACK)
#if defined(__parallel)
     CALL pdsaupd(arpack_comm, ido, bmat, nloc, which, nev, tol, resid, &
                 ncv, cv, nloc, iparam, ipntr, workd, workl, &
                 lworkl, info )
#else
     CALL dsaupd ( ido, bmat, nloc, which, nev, tol, resid, &
                 ncv, cv, nloc, iparam, ipntr, workd, workl, &
                 lworkl, info )
#endif
#endif

     IF (info.LT.0) THEN
        CALL stop_program(routineP,"error in dsaupd",globenv) 
     ENDIF

     IF (ido .EQ. -1 .OR. ido .EQ. 1) THEN
       vin =>workd(ipntr(1):)
       vout=>workd(ipntr(2):)
       IF (PRESENT(op_symv)) THEN
          CALL arpack_op_symv(op_symv, vin, vout, globenv)
       END IF
     ELSE
       EXIT
     ENDIF
  ENDDO

!TC  IF (info.LT.0) THEN
!TC        WRITE(6,*) info
!TC     CALL stop_program(routineP,"error in dsaupd",globenv) 
!TC  ENDIF

  rvec = .TRUE.
  ierr = 0 

#if defined(__ARPACK)
#if defined(__parallel)
  CALL pdseupd ( arpack_comm, rvec, 'All', SELECT, d, cv, nloc, sigma, &        
                bmat, nloc, which, nev, tol, resid, ncv, cv, nloc, &
                iparam, ipntr, workd, workl, lworkl, ierr )
#else
  CALL dseupd ( rvec, 'All', SELECT, d, cv, nloc, sigma, &        
                bmat, nloc, which, nev, tol, resid, ncv, cv, nloc, &
                iparam, ipntr, workd, workl, lworkl, ierr )
#endif
#endif

  IF ( ierr .NE. 0) THEN
     CALL stop_program(routineP,"dseupd error",globenv)
  ENDIF

  IF ( iparam(5) .NE. nev) THEN
     CALL stop_program(routineP,"dseupd not all converged",globenv)
  ENDIF

  IF ( info .EQ. 1) THEN
     CALL stop_program(routineP,"max number of iterations reached",globenv)
  ENDIF

  ! copy eigenvalues
  eval(1:nev)=d(1:nev)
  DO i=1,nev
     ! is there any clean way to do this ... ?
     workd(1:nloc)=cv(1:nloc,i)
     CALL copy_vector_to_blacs(workd,vmatrix,i)
  ENDDO

  DEALLOCATE(cv,workl,workd,d,resid,ax,SELECT)
  CALL timestop(0.0_wp,handle)
  
END SUBROUTINE arpack_diag_sy

! *****************************************************************************

END MODULE

