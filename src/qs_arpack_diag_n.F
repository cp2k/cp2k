!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_arpack_diag_n *
!!
!!   NAME
!!     qs_arpack_diag_n
!!
!!   FUNCTION
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
MODULE qs_arpack_diag_n
  ! the operator environments and operators
  USE cp_error_handling,               ONLY: cp_simple_assert
  USE cp_fm_utils,                     ONLY: copy_rep_vec_to_fm_col
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_arpack_tddfpt,                ONLY: arpack_op_tddfpt,&
                                             op_tddfpt_type
  USE qs_arpack_utils,                 ONLY: arpack_confirm
  USE qs_blacs,                        ONLY: blacs_find_vector_distribution,&
                                             cp_fm_get_info,&
                                             cp_full_matrix_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE termination,                     ONLY: stop_program
  IMPLICIT NONE
  
  LOGICAL, PARAMETER :: DEBUG_THIS_MODULE = .FALSE.

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_arpack_diag_n"

  PRIVATE
  PUBLIC :: arpack_diag_n

CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_arpack_diag_n/arpack_diag_n *
!!
!!   NAME
!!     arpack_diag_n
!!
!!   SYNOPSIS
!!     Subroutine arpack_diag_n(n_ev, eval, evec, glob_env, op_tddfpt)
!!       Implicit None
!!       Integer, Intent (IN):: n_ev
!!       Real(Kind=dbl), Dimension(:,:), Pointer:: eval
!!       Type(cp_full_matrix_type), Pointer:: evec
!!       Type(global_environment_type), Intent (IN):: glob_env
!!       Type(op_tddfpt_type), Intent (INOUT), Optional:: op_tddfpt
!!     End Subroutine arpack_diag_n
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
!!     exactly ONE of the following optional arguments must also be specified
!!     - op_tddfpt : operator for the tddfpt calculations
!!
!!   NOTES
!!     - Since the problem at hand is not symmetric (hermitian), there
!!       can be complex eigenvalues. Therefore the array eval must
!!       be two-dimensional (one column for real part and one for the
!!       imaginary part).
!!
!!   SEE ALSO
!!     arpack_diag_sy
!!
!!****
  SUBROUTINE arpack_diag_n(n_ev, eval, evec, glob_env, op_tddfpt)
  
    IMPLICIT NONE

  	! arguments
    INTEGER, INTENT(in)                           :: n_ev
    REAL(dbl), DIMENSION(:,:), POINTER            :: eval
    TYPE(cp_full_matrix_type), POINTER            :: evec
    TYPE(global_environment_type), INTENT(in)     :: glob_env
    TYPE(op_tddfpt_type), INTENT(inout), OPTIONAL :: op_tddfpt
  
  	! locals
    INTEGER   :: ev, stat
    INTEGER   :: op_present        ! how many optional operators are present
    CHARACTER(len=*), PARAMETER :: routineN = "arpack_diag_n", &
                                   routineP = moduleN//"/"//routineN
    
    ! debugging
    INTEGER   :: arpack_counter
    
    ! arpack specific
    CHARACTER :: bmat*1, which*2
    LOGICAL   :: rvec
    LOGICAL, DIMENSION(:), POINTER :: ev_select
    INTEGER   :: ido, &
                 nrows, &        ! size of matrix on local node?
                 ldv, &            ! leading dimension ?
                 n_kv, &           ! number of krylov space vectors
                 lworkl, &         ! leading dimension of workl
                 info                 
    REAL(dbl) :: tol               ! tolerance
    REAL(dbl) :: sigmar, sigmai    ! unused but must be provided
    INTEGER, DIMENSION(14)             :: param, pntr
    REAL(dbl), DIMENSION(:), POINTER   :: resid, workev, workd, workl, vin, vout
    REAL(dbl), DIMENSION(:,:), POINTER :: v
    
    CALL arpack_confirm(glob_env)
    
    ! check if an operator has been specified
    op_present = 0
    IF (PRESENT(op_tddfpt)) THEN
       op_present = op_present + 1
    END IF
    IF (op_present /= 1) THEN
       CALL stop_program(routineP, "exactly 1 operator must be specified", glob_env)
    END IF

    CALL cp_fm_get_info(matrix=evec, nrow_global=nrows)

    CALL check_input(n_ev, nrows, eval, evec, glob_env)
    
    !--------------------------------------!
    ! calculate the number of krylov space !
    ! vectors (must be at least n_ev+2)    !
    !--------------------------------------!
    n_kv = MIN(2*n_ev, nrows)
    IF (n_kv < n_ev+2) n_kv = n_ev+2

    !-----------------------------------!
    ! set up variables needed by arpack !
    !-----------------------------------!
    ido      = 0
    which    = 'SR'                    ! eigenvalues with the smallest real part
    tol      = 1.0e-15_dbl             ! hardcoded ?
    ldv      = nrows
    param    = 0
    param(1) = 1                       ! exact shift mode
    param(3) = 1500                    ! maximum number of iterations
    param(7) = 1                       ! simple mode (A*x=lambda*x , B=I)
    pntr     = 0
    lworkl   = 3*(n_kv**2) + 6*(n_kv)  ! do larger numbers give better performance?
    info     = 0                       ! no starting vector as of now

    IF (PRESENT(op_tddfpt)) THEN
       bmat = 'I'
    END IF

    ! allocate some arrays
    ALLOCATE(resid(nrows), v(ldv, n_kv), &
             workev(3*n_kv), workd(3*nrows), workl(lworkl), &
             ev_select(n_kv), &
             stat=stat)
    IF (stat /= 0) &
         CALL stop_program(routineP, "allocation of temporary variables failed", glob_env)

    IF (DEBUG_THIS_MODULE) THEN
       WRITE (*,*) routineP
       WRITE (*,*) "nrows : ", nrows, "n_ev : ", n_ev, "ldv : ", ldv, "n_kv : ", n_kv
    END IF

#if defined(__ARPACK)
    arpack_counter = 0
    reverse_comm: DO 

       IF (DEBUG_THIS_MODULE) THEN
          arpack_counter = arpack_counter + 1
          IF (glob_env%ionode) &
               WRITE (glob_env%scr,'(1X,A,I4)') 'iteration ', arpack_counter
       END IF
       
       CALL dnaupd(ido, bmat, nrows, which, n_ev, tol, resid, &
                   n_kv, v, ldv, param, pntr, workd, workl, lworkl, info)
       IF (info < 0) THEN
          IF (glob_env%ionode) WRITE (glob_env%scr, '(A,I5)') 'info = ', info
          CALL stop_program(routineP, "error in dnaupd", glob_env)
       END IF

       IF (ido==-1.OR.ido==1) THEN
          vin  => workd(pntr(1):pntr(1)+nrows-1)
          CPAAssert(SIZE(vin)==nrows)
          vout => workd(pntr(2):pntr(2)+nrows-1)
          CPAAssert(SIZE(vout)==nrows)

          !----------------------------!
          ! call the provided operator !
          !----------------------------!
          IF (PRESENT(op_tddfpt)) THEN
             CALL arpack_op_tddfpt(op_tddfpt, vin, vout, glob_env)
          END IF

       ELSE IF (ido == 99) THEN
          EXIT reverse_comm
       ELSE
          IF (glob_env%ionode) WRITE (glob_env%scr, '(A,I5)') 'ido = ', ido
          CALL stop_program(routineP, "error in reverse communication")
       END IF

    END DO reverse_comm

    rvec = .TRUE.          ! we want the eigenvectors
    ev_select = .TRUE.     ! just to be sure
    sigmar = 0.0_dbl; sigmai = 0.0_dbl
    CALL dneupd(rvec, 'A', ev_select, eval(:,1), eval(:,2), v, ldv, &
                sigmar, sigmai, workev, &
                bmat, nrows, which, n_ev, tol, resid, &
                n_kv, v, ldv, param, pntr, workd, workl, lworkl, info)
#endif

    ! copy the eigenvectors back into the evec array
    DO ev=1, n_ev
       vout => v(:,ev)
       CALL copy_rep_vec_to_fm_col(vout, evec, ev)
    END DO

    IF (info /= 0) THEN
       IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A,I5)') routineP//": info = ", info
       CALL stop_program(routineP, "error in dneupd", glob_env)
    END IF

    DEALLOCATE(resid, v, workev, workd, workl, ev_select, stat=stat)
    IF (stat /= 0) &
         CALL stop_program(routineP, "deallocation of temporary variables failed", glob_env)

  END SUBROUTINE arpack_diag_n

! *****************************************************************************

  SUBROUTINE check_input(n_ev, nrows, eval, evec, glob_env)
    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(in)                       :: n_ev
    INTEGER, INTENT(in)                       :: nrows
    REAL(dbl), DIMENSION(:,:), POINTER        :: eval
    TYPE(cp_full_matrix_type), POINTER        :: evec
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    INTEGER                     :: evec_cols
    CHARACTER(len=*), PARAMETER :: routineN = "check_input", &
                                   routineP = moduleN//"/"//routineN

    ! code
    CALL cp_fm_get_info(evec, ncol_global=evec_cols)

    IF (n_ev < 1) &
         CALL stop_program(routineP, "n_ev must be at least 1", glob_env)

    IF (SIZE(eval,1) < n_ev+1 .OR. SIZE(eval,2) < 2) &
         CALL stop_program(routineP, "array eval too small", glob_env)

    IF (n_ev > nrows - 2) &
       CALL stop_program(routineP, "n_kv >= n_ev+2 and n_kv <= nrows but n_ev > nrows-2")

    IF (evec_cols < n_ev)  &
       CALL stop_program(routineP, "evec needs at least n_ev columns", glob_env)

  END SUBROUTINE check_input

! *****************************************************************************

END MODULE qs_arpack_diag_n
