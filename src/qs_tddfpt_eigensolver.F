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
MODULE qs_tddfpt_eigensolver

  USE cp_error_handling,    ONLY : cp_simple_assert
  USE cp_fm_basic_linalg,   ONLY : cp_fm_gemm, cp_fm_symm, cp_fm_add
  USE cp_fm_cholesky,       ONLY : cp_fm_cholesky_decompose, cp_fm_cholesky_invert
  USE cp_fm_pool_types,     ONLY : cp_fm_pool_p_type,&
                                   cp_fm_pool_type,&
                                   fm_pools_create_matrix_vect,&
                                   fm_pools_give_back_matrix_vect
  USE cp_fm_struct,         ONLY : cp_fm_struct_type, &
                                   cp_fm_struct_create, &
                                   optimal_blacs_col_block_size, &
                                   optimal_blacs_row_block_size
  USE cp_fm_basic_linalg,   ONLY : cp_sm_fm_multiply, &
                                   cp_fm_trace, cp_fm_scale
  USE cp_fm_types,          ONLY : cp_full_matrix_p_type, &
                                   cp_full_matrix_type, &
                                   cp_fm_create, &
                                   cp_fm_release, &
                                   cp_fm_get_element
  USE cp_fm_utils,          ONLY : copy_rep_vec_to_fm, &
                                   copy_fm_to_rep_vec
  USE cp_log_handling,      ONLY : cp_to_string
  USE global_types,         ONLY : global_environment_type
  USE kinds,                ONLY : dbl
  USE qs_arpack_utils,      ONLY : arpack_confirm
  USE qs_blacs,             ONLY : cp_sm_plus_fm_fm_t, &
                                   copy_real_matrix_to_full_matrix, &
                                   cp_sm_fm_multiply, &
                                   copy_blacs_to_blacs_matrix
  USE qs_environment_types, ONLY : qs_environment_type, get_qs_env
  USE qs_p_env_methods,     ONLY : p_op_l1, p_op_l2, p_preortho, p_postortho
  USE qs_p_env_types,       ONLY : qs_p_env_type
  USE qs_tddfpt_types
  USE qs_tddfpt_utils
  USE sparse_matrix_types,  ONLY : set_matrix, &
                                   replicate_matrix_structure
  USE termination,          ONLY : stop_program

  IMPLICIT NONE

  INTERFACE tddfpt_es
     MODULE PROCEDURE tddfpt_es_lanczos
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_eigensolver"

  PRIVATE

  PUBLIC :: tddfpt_es

CONTAINS

! *****************************************************************************

  SUBROUTINE tddfpt_es_arpack(t_env, p_env, qs_env, glob_env)
    
    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(INOUT)      :: t_env
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! L0CALS
    CHARACTER :: bmat*1, which*2
    INTEGER   :: spin, n_spins
    INTEGER   :: index, ev
    INTEGER   :: n_ev              ! number of eigenvalues
    INTEGER   :: ido, &
                 nrows, &          ! size of matrix on local node?
                 ldv, &            ! leading dimension ?
                 n_kv, &           ! number of krylov space vectors
                 lworkl, &         ! leading dimension of workl
                 info                 
    INTEGER                                            :: stat
    INTEGER, DIMENSION(14)                             :: param, pntr
    LOGICAL                                            :: rvec
    LOGICAL, DIMENSION(:), POINTER                     :: ev_select
    REAL(dbl)                                          :: tol               ! tolerance
    REAL(dbl)                                          :: sigmar, sigmai    ! unused but must be provided
    REAL(dbl), DIMENSION(:), POINTER                   :: resid, workev, workd, workl, vin, vout
    REAL(dbl), DIMENSION(:,:), POINTER                 :: v
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X, R
    TYPE(cp_fm_struct_type), POINTER                   :: v_struct
    INTEGER                                            :: v_rows, v_cols

    CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_es_arpack", &
                                   routineP = moduleN//"/"//routineN


    CALL arpack_confirm(glob_env)

    NULLIFY(ao_mo_fm_pools, resid, v, workev, workd, workl, ev_select, X, R)
    n_ev = t_env%control%n_ev
    n_spins = qs_env%dft_control%nspins

    ! muss noch abgeaendert werden fuer parallel ??
    ! CALL cp_fm_get_info(matrix=evec, nrow_global=nrows)
    nrows = 0
    DO spin=1, n_spins
       nrows = nrows + p_env%n_ao(spin)*p_env%n_mo(spin)
    END DO

    !--------------------------------------!
    ! calculate the number of krylov space !
    ! vectors (must be at least n_ev+2)    !
    !--------------------------------------!
    n_kv = MIN(3*n_ev, nrows) ! cannot be bigger than nrows
    IF (n_kv < n_ev+2) n_kv = n_ev+2

    !-----------------------------------!
    ! set up variables needed by arpack !
    !-----------------------------------!
    bmat     = 'I'
    ido      = 0
    which    = 'SR'                    ! eigenvalues with the smallest real part
    tol      = 1.0e-9_dbl             ! hardcoded ?
    ldv      = nrows
    param    = 0
    param(1) = 1                       ! exact shift mode
    param(3) = 10000                     ! maximum number of iterations
    param(7) = 1                       ! simple mode (A*x=lambda*x , B=I)
    pntr     = 0
    lworkl   = 3*(n_kv**2) + 6*(n_kv)  ! do larger numbers give better performance?
    info     = 1                       ! use a starting vector

    !------------------------!
    ! allocate the variables !
    !------------------------!
    ALLOCATE(resid(nrows), v(ldv, n_kv), &
             workev(3*n_kv), workd(3*nrows), workl(lworkl), &
             ev_select(n_kv), &
             X(n_spins), R(n_spins), &
             stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, "allocation of temporary variables failed", glob_env)
    CALL get_qs_env(qs_env, ao_mo_fm_pools=ao_mo_fm_pools)! ,para_env=para_env)
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=X,&
                                     name=routineP//":X")!, error=error)
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=R,&
                                     name=routineP//":R")!, error=error)

    !--------------------------!
    ! fill the starting vector !
    !--------------------------!
    index = 1
    ev = 1 ! use values stored in X(1)
    DO spin=1, n_spins
!!       CALL print_fm_matrix(t_env%X(ev)%spins(spin)%matrix, p_env%n_ao(spin), p_env%n_mo(spin), &
!!            "X.dat")
       vout  => resid(index:index+(p_env%n_ao(spin)*p_env%n_mo(spin))-1)
       CALL copy_fm_to_rep_vec(t_env%X(ev)%spins(spin)%matrix, vout)
       index = index+(p_env%n_ao(spin)*p_env%n_mo(spin))
    END DO

    !-----------------------------------!
    ! start the arpack reverse comm run !
    !-----------------------------------!
#if defined(__ARPACK)
    reverse_comm: DO 

       CALL dnaupd(ido, bmat, nrows, which, n_ev, tol, resid, &
                   n_kv, v, ldv, param, pntr, workd, workl, lworkl, info)
       IF (info < 0) THEN
          IF (glob_env%ionode) WRITE (glob_env%scr, '(A,I5)') 'info = ', info
          CALL stop_program(routineP, "error in dnaupd", glob_env)
       END IF

       IF (ido==-1.OR.ido==1) THEN

          !----------------------------------------------------------!
          ! Call the superoperator on the supervector X, result is R !
          !----------------------------------------------------------!          
          ! vin => X
          index=pntr(1)
          DO spin=1, n_spins
             vin  => workd(index:index+(p_env%n_ao(spin)*p_env%n_mo(spin))-1)
             CALL copy_rep_vec_to_fm(vin, X(spin)%matrix)
             index = index+(p_env%n_ao(spin)*p_env%n_mo(spin))
          END DO
          CALL apply_op(X, R, t_env, p_env, qs_env, glob_env)
          ! X => vout
          index=pntr(2)
          DO spin=1, n_spins
             vout  => workd(index:index+(p_env%n_ao(spin)*p_env%n_mo(spin))-1)
             CALL copy_fm_to_rep_vec(X(spin)%matrix, vout)
             index = index+(p_env%n_ao(spin)*p_env%n_mo(spin))
          END DO
          
       ELSE IF (ido == 99) THEN
          EXIT reverse_comm
       ELSE
          IF (glob_env%ionode) WRITE (glob_env%scr, '(A,I5)') 'ido = ', ido
          CALL stop_program(routineP, "error in reverse communication")
       END IF

    END DO reverse_comm

!!    call tddfpt_es_lanczos()

    rvec = .TRUE.          ! we want the eigenvectors
    ev_select = .TRUE.     ! just to be sure
    sigmar = 0.0_dbl; sigmai = 0.0_dbl
    CALL dneupd(rvec, 'A', ev_select, t_env%ev(:,1), t_env%ev(:,2), v, ldv, &
                sigmar, sigmai, workev, &
                bmat, nrows, which, n_ev, tol, resid, &
                n_kv, v, ldv, param, pntr, workd, workl, lworkl, info)
#endif
    
    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, X)
    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, R)
    NULLIFY(X, R)

!!    ! COPY THE VECTORS FROM V TO T_ENV%X
!!    DO ev=1, n_ev
!!       index=1
!!       DO spin=1, n_spins
!!          vout  => v(index:index+(p_env%n_ao(spin)*p_env%n_mo(spin))-1,ev)
!!          CALL copy_fm_to_rep_vec(t_env%X(ev)%spins(spin)%matrix, vout)
!!          index = index+(p_env%n_ao(spin)*p_env%n_mo(spin))
!!       END DO
!!    END DO

    DEALLOCATE(resid, v, workev, workd, workl, ev_select, stat=stat)

  END SUBROUTINE tddfpt_es_arpack

! *****************************************************************************

  ! X is changed with this routine
  SUBROUTINE apply_op(X, R, t_env, p_env, qs_env, glob_env)

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X     ! on entering the left hand side matrix
                                                                ! on exit the result of the operation
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: R     ! temporary matrix
    TYPE(tddfpt_env_type), INTENT(INOUT)               :: t_env
    TYPE(qs_p_env_type), POINTER                       :: p_env
    TYPE(qs_environment_type)                          :: qs_env
    TYPE(global_environment_type), INTENT(IN)          :: glob_env

    ! LOCALS
    INTEGER                                            :: spin, n_spins
    TYPE(cp_fm_struct_type), POINTER                   :: v_struct
    INTEGER                                            :: v_rows, v_cols
    INTEGER, SAVE                                      :: counter = 0
    REAL(DBL)                                          :: trace

    CHARACTER(LEN=*), PARAMETER                        :: routineN = "apply_op", &
                                                          routineP = moduleN//"/"//routineN

    counter = counter + 1
    
    n_spins = qs_env%dft_control%nspins

    ! print the qs_env%k and qs_env%s matrices
!!    DO spin=1, n_spins
!!       CALL copy_real_matrix_to_full_matrix(qs_env%s(spin)%matrix, inv_U(spin)%matrix)
!!       CALL print_fm_matrix(inv_U(spin)%matrix, p_env%n_ao(spin), p_env%n_ao(spin), &
!!            "S_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat", sym='T')
!!       CALL copy_real_matrix_to_full_matrix(qs_env%k(spin)%matrix, inv_U(spin)%matrix)
!!       CALL print_fm_matrix(inv_U(spin)%matrix, p_env%n_ao(spin), p_env%n_ao(spin), &
!!            "H_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat", sym='T')
!!    END DO

!!    IF (counter == 2) THEN
!!       DO spin=1, n_spins
!!          CALL print_fm_matrix(X(spin)%matrix, p_env%n_ao(spin), p_env%n_mo(spin), &
!!               "X_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat")
!!       END DO
!!    END IF

!!    CALL p_preortho(p_env, qs_env, X) ! make X orthogonal to psi0

!!    ! print X_ortho
!!    DO spin=1, n_spins
!!       CALL print_fm_matrix(X(spin)%matrix, p_env%n_ao(spin), p_env%n_mo(spin), &
!!            "X_ortho_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat")
!!    END DO

    CALL p_op_l1(p_env, qs_env, X, R)  ! acts on both spins

!!    ! print the result of l1(x)
!!    DO spin=1, n_spins
!!       CALL print_fm_matrix(R(spin)%matrix, p_env%n_ao(spin), p_env%n_mo(spin), &
!!            "R1_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat")
!!    END DO

!!    DO spin=1, n_spins                 ! calculate p1(spin)
!!       CALL set_matrix(p_env%p1(spin)%matrix, value=0.0_dbl) ! optimize?
!!       CALL cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
!!            matrix_v=p_env%psi0d(spin)%matrix,&
!!            matrix_g=X(spin)%matrix,&
!!            ncol=p_env%n_mo(spin))            
!!    END DO
!!    CALL p_op_l2(p_env, qs_env, p_env%p1, R, alpha=1.0_dbl, beta=1.0_dbl)
!!    CALL p_op_l2(p_env, qs_env, p_env%p1, R, alpha=2.0_dbl, beta=1.0_dbl)

    CALL p_postortho(p_env, qs_env, R)


!!    DO spin=1, n_spins
!!       CALL print_fm_matrix(R(spin)%matrix, p_env%n_ao(spin), p_env%n_mo(spin), &
!!            "R2_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat")
!!    END DO

!    x=invs*r
    DO spin=1, n_spins
       CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                       1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
                       0.0_dbl, X(spin)%matrix)
    END DO

!!    CALL p_op_l1(p_env, qs_env, X, R) ! acts on both spins
!!
!!!    x=invs*r
!!    DO spin=1, n_spins
!!       CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
!!                       1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
!!                       0.0_dbl, X(spin)%matrix)
!!    END DO
!!
!!!!    STOP
!!
!!    IF (counter==1) THEN
!!    DO spin=1, n_spins
!!       CALL print_fm_matrix(X(spin)%matrix, p_env%n_ao(spin), p_env%n_mo(spin), &
!!            "R3_"//TRIM(ADJUSTL(cp_to_string(spin)))//".dat")
!!    END DO
!!    END IF

!!!    !----------------------!
!!!    ! normalize the result !
!!!    !----------------------!
!!!    DO spin=1, n_spins
!!!       CALL cp_sm_fm_multiply(qs_env%s(spin)%matrix, X(spin)%matrix, R(spin)%matrix, &
!!!            p_env%n_mo(spin), glob_env%para_env, 1.0_dbl, 0.0_dbl)
!!!       CALL cp_fm_trace(X(spin)%matrix, R(spin)%matrix, trace)
!!!       CPAAssert(trace>0)
!!!       trace = 1.2_dbl
!!!       CALL cp_fm_scale(1.0_dbl/sqrt(trace), X(spin)%matrix)
!!!    END DO

  END SUBROUTINE apply_op

! *****************************************************************************

  SUBROUTINE print_fm_matrix(matrix, imax, jmax, name, sym)

    IMPLICIT NONE

    TYPE(cp_full_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)                :: imax, jmax
    CHARACTER(*), INTENT(IN)           :: name
    CHARACTER, OPTIONAL, INTENT(IN)    :: sym

    INTEGER   :: i, j, ioerr
    REAL(DBL) :: element
    CHARACTER :: symmetric

    symmetric = 'F'
    IF (PRESENT(sym)) symmetric=sym

    OPEN(UNIT=8, FILE=name, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)

    DO i=1, imax
       DO j=1, jmax
          IF ((symmetric == 'T') .AND. (i>j)) THEN
             CALL cp_fm_get_element(matrix, j, i, element)
          ELSE
             CALL cp_fm_get_element(matrix, i, j, element)
          END IF
          WRITE (8, '(1X,F20.10)') element
       END DO
    END DO

    CLOSE(8)

  END SUBROUTINE print_fm_matrix

! *****************************************************************************

  SUBROUTINE tddfpt_es_lanczos(t_env, p_env, qs_env, glob_env)
    
    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(INOUT)      :: t_env
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! L0CALS
    INTEGER   :: spin, n_spins, &
                 n_ev, n_kv, n_rows, &
                 n_reortho, &
                 i, j, index, &
                 stat
    REAL(dbl) :: tol, &
                 beta, &
                 tmp
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X, R, kv

    CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_es_lanczos", &
                                   routineP = moduleN//"/"//routineN


    NULLIFY(ao_mo_fm_pools, X, R, kv)
    n_ev      = t_env%control%n_ev
    n_spins   = qs_env%dft_control%nspins
    tol       = 1.0e-10_dbl
    n_reortho = 2 ! muss konfigurierbar sein

    !----------------------------!
    ! determine number of krylov !
    ! space vetors to use        ! 
    !----------------------------!
    n_rows = 0
    DO spin=1, n_spins
       n_rows = n_rows + p_env%n_ao(spin)*p_env%n_mo(spin)
    END DO
    n_kv = MIN(5*n_ev, n_rows) ! cannot be bigger than nrows

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    ALLOCATE(kv(n_spins*n_kv), X(n_spins), R(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
         "allocation of temporary variables failed", glob_env)
    CALL get_qs_env(qs_env, ao_mo_fm_pools=ao_mo_fm_pools)
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools, matrixes=kv,&
                                     name=routineP//":kv")
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=X,&
                                     name=routineP//":X")
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools,matrixes=R,&
                                     name=routineP//":R")

    !------------------------!
    ! provide starting value !
    !------------------------!
    DO spin=1, n_spins
       CALL copy_blacs_to_blacs_matrix(t_env%X(1)%spins(spin)%matrix, kv(spin)%matrix)
    END DO

    DO j=1, n_kv-1

       !---------------------------------------------!
       ! copy current krylov vector to working space !
       !---------------------------------------------!
       index = n_spins*(j-1)
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(kv(index+spin)%matrix, X(spin)%matrix)
       END DO


       CALL apply_op(X, R, t_env, p_env, qs_env, glob_env) ! result in X


       DO i=1, n_reortho
          CALL reorthogonalize(X, kv, j, n_spins)
       END DO


       beta = 0.0_dbl
       DO spin=1, n_spins
          CALL cp_fm_trace(X(spin)%matrix, X(spin)%matrix, tmp)
          beta = beta + tmp
       END DO


       IF (ABS(beta) < tol) EXIT

       !-------------!
       ! normalize X ! 
       !-------------!
       DO spin=1, n_spins
          CALL cp_fm_scale(1._dbl/beta, X(spin)%matrix)
       END DO

       !----------------------------------------------------!
       ! copy X(1:n_spins) to kv(n_spins*j+1:n_spins*(j+1)) !
       !----------------------------------------------------!
       index = n_spins*j
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(X(spin)%matrix, kv(index+spin)%matrix)
       END DO

    END DO

    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, X)
    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, R)
    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, kv)
    NULLIFY(X, R, kv)
    
    DEALLOCATE(X, R, kv)

  END SUBROUTINE tddfpt_es_lanczos

! *****************************************************************************

  SUBROUTINE reorthogonalize(X, kv, j, n_spins)
    IMPLICIT NONE
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X, kv
    INTEGER, INTENT(IN) :: j, n_spins
    INTEGER :: i, index, spin
    REAL(DBL) :: tmp, dot_product

    DO i=1, j-1
       dot_product = 0.0_dbl
       DO spin=1, n_spins
          index = n_spins*(i-1)+spin
          CALL cp_fm_trace(X(spin)%matrix, kv(index)%matrix, tmp)
          dot_product = dot_product + tmp
       END DO
       DO spin=1, n_spins
          index = n_spins*(i-1)+spin
          CALL cp_fm_add(1.0_dbl, X(spin)%matrix, -dot_product, kv(index)%matrix)
       END DO
    END DO
  END SUBROUTINE reorthogonalize

! *****************************************************************************

END MODULE qs_tddfpt_eigensolver
