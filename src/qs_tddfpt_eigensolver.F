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
  USE cp_fm_basic_linalg,   ONLY : cp_fm_gemm, cp_fm_symm, cp_fm_scale_and_add
  USE cp_fm_cholesky,       ONLY : cp_fm_cholesky_decompose, cp_fm_cholesky_invert
  USE cp_fm_pool_types,     ONLY : cp_fm_pool_p_type,&
                                   cp_fm_pool_type,&
                                   fm_pool_create_matrix_vect,&
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
                                   cp_fm_get_element, &
                                   cp_fm_set_element
  USE cp_fm_utils,          ONLY : copy_rep_vec_to_fm, &
                                   copy_fm_to_rep_vec
  USE cp_log_handling,      ONLY : cp_to_string
  USE global_types,         ONLY : global_environment_type
  USE kinds,                ONLY : dbl
  USE qs_arpack_utils,      ONLY : arpack_confirm
  USE qs_blacs,             ONLY : cp_sm_plus_fm_fm_t, &
                                   copy_real_matrix_to_full_matrix, &
                                   copy_blacs_to_blacs_matrix, &
                                   blacs_syevd
  USE qs_environment_types, ONLY : qs_environment_type, get_qs_env
  USE qs_p_env_methods,     ONLY : p_op_l1, p_op_l2, p_preortho, p_postortho
  USE qs_p_env_types,       ONLY : qs_p_env_type
  USE qs_tddfpt_types
  USE qs_tddfpt_utils
  USE sparse_matrix_types,  ONLY : set_matrix, &
                                   replicate_matrix_structure, &
                                   real_matrix_p_type
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
                 row, col, &
                 i, j, index, &
                 stat, n_kv_fac
    REAL(dbl) :: tol, &
                 beta, &
                 tmp, &
                 dot_product
    REAL(dbl), DIMENSION(:), ALLOCATABLE               :: eigenvals
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X, R, kv, T, Tev, opkv
    TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

    CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_es_lanczos", &
                                   routineP = moduleN//"/"//routineN


    NULLIFY(ao_mo_fm_pools, X, R, kv)
    n_ev      = t_env%control%n_ev
    n_spins   = qs_env%dft_control%nspins
    tol       = 1.0e-10_dbl
    n_reortho = t_env%control%n_reortho
    n_kv_fac  = t_env%control%n_kv_fac + n_ev

    !----------------------------!
    ! determine number of krylov !
    ! space vetors to use        ! 
    !----------------------------!
    n_rows = 0
    DO spin=1, n_spins
       n_rows = n_rows + p_env%n_ao(spin)*p_env%n_mo(spin)
    END DO
    n_kv = (n_kv_fac) + (n_ev-1)*(n_kv_fac/2)
    n_kv = MIN(n_kv, n_rows) ! cannot be bigger than nrows


    !----------------------!
    ! allocate the vectors !
    !----------------------!
    ALLOCATE(eigenvals(n_kv))
    CALL get_qs_env(qs_env, ao_mo_fm_pools=ao_mo_fm_pools)
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools, X, &
                                     name=routineP//":X")
    CALL fm_pools_create_matrix_vect(ao_mo_fm_pools, R, &
                                     name=routineP//":R")

    ALLOCATE(kv(n_spins*n_kv), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, "allocation of kv failed", glob_env)
    DO i=1, n_spins*n_kv
       spin = MOD(i, n_spins)+1
       NULLIFY(fm_struct, kv(i)%matrix)
       CALL cp_fm_struct_create(fm_struct, glob_env%para_env, glob_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_mo(spin), &
                                optimal_blacs_row_block_size, optimal_blacs_col_block_size)
       CALL cp_fm_create(kv(i)%matrix, fm_struct, routineP//"kv-"//TRIM(cp_to_string(i)))
    END DO
    ALLOCATE(opkv(n_spins*n_kv), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, "allocation of opkv failed", glob_env)
    DO i=1, n_spins*n_kv
       spin = MOD(i, n_spins)+1
       NULLIFY(fm_struct, opkv(i)%matrix)
       CALL cp_fm_struct_create(fm_struct, glob_env%para_env, glob_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_mo(spin), &
                                optimal_blacs_row_block_size, optimal_blacs_col_block_size)
       CALL cp_fm_create(opkv(i)%matrix, fm_struct, routineP//"opkv-"//TRIM(cp_to_string(i)))
    END DO
    NULLIFY(fm_struct)



    !---------------------------------------!
    ! provide starting vector and normalize !
    !---------------------------------------!
    dot_product = 0.0_dbl
    DO spin=1, n_spins
       CALL copy_blacs_to_blacs_matrix(t_env%X(1)%spins(spin)%matrix, kv(spin)%matrix)
       CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
            kv(spin)%matrix, &
            R(spin)%matrix, &
            p_env%n_mo(spin), &
            glob_env%para_env, 1.0_dbl, 0.0_dbl)
       CALL cp_fm_trace(R(spin)%matrix, kv(spin)%matrix, tmp)
       dot_product = dot_product + tmp
    END DO
    IF (dot_product > tol) THEN
       DO spin=1, n_spins
          CALL cp_fm_scale(SQRT(1.0_dbl/dot_product), kv(spin)%matrix)
       END DO
    ELSE
       CALL stop_program(routineP, "The provided start vector is almost zero.")
    END IF


    beta = 0.0_dbl ! beta of the first vector is 0
    

    krylov: DO j=1, n_kv-1
       

       !---------------------------------------------!
       ! copy current krylov vector to working space !
       !---------------------------------------------!
       index = n_spins*(j-1)
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(kv(index+spin)%matrix, X(spin)%matrix)
       END DO
       
       
       CALL apply_op(X, R, t_env, p_env, qs_env, glob_env) ! result in X


       !-------------------------------!
       ! keep result for later use (T) !
       !-------------------------------!
       index = n_spins*(j-1)
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(X(spin)%matrix, opkv(index+spin)%matrix)
       END DO


       DO i=1, n_reortho
          CALL reorthogonalize(X, kv, qs_env%s, R, j, n_spins, glob_env)
       END DO
       

       !-------------------------!
       ! calc beta=transp(X)*S*X !
       !-------------------------!
       beta = 0.0_dbl
       DO spin=1, n_spins
          CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
               X(spin)%matrix, &
               R(spin)%matrix, &
               p_env%n_mo(spin), &
               glob_env%para_env, 1.0_dbl, 0.0_dbl)
          CALL cp_fm_trace(R(spin)%matrix, X(spin)%matrix, tmp)
          beta = beta + tmp
       END DO


       !-------------!
       ! normalize X ! 
       !-------------!
       IF (ABS(beta) > tol) THEN
          DO spin=1, n_spins
             CALL cp_fm_scale(SQRT(1._dbl/beta), X(spin)%matrix)
          END DO
       ELSE
          ! we can quit n_kv=j
          n_kv = j
          IF (glob_env%ionode) THEN
             WRITE (*,*)
             WRITE (*,*) "Converged with ", n_kv," krylov vectors!"
             WRITE (*,*)
          END IF
          EXIT krylov
       END IF


       !---------------------------!
       ! write current vector back !
       !---------------------------!
       index = n_spins*j
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(X(spin)%matrix, kv(index+spin)%matrix)
       END DO


    END DO krylov

#if 0
    !--------------------------------!
    ! print the krylov space vectors !
    !--------------------------------!
    IF (glob_env%ionode) THEN
       DO i=1, n_kv
          DO spin=1, n_spins
             index = n_spins*(i-1)+spin
             CALL print_fm_matrix(kv(index)%matrix, & 
                  p_env%n_ao(spin), p_env%n_mo(spin), &
                  "kv-"//TRIM(ADJUSTL(cp_to_string(index)))//".dat")
          END DO
       END DO
    END IF
#endif

#if 0
    !-------------------------------------!
    ! check if the vectors are orthogonal !
    !-------------------------------------!
    DO i=1, n_kv
       DO spin=1, n_spins
          index=n_spins*(i-1)+spin
          CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
                  kv(i)%matrix, &
                  R(spin)%matrix, &
                  p_env%n_mo(spin), &
                  glob_env%para_env, 1.0_dbl, 0.0_dbl)
          END DO
          DO j=1, i-1
             dot_product = 0.0_dbl
             DO spin=1, n_spins
                index=n_spins*(j-1)+spin
                CALL cp_fm_trace(kv(index)%matrix, R(spin)%matrix, tmp)
                dot_product = dot_product + tmp
             END DO
             IF (dot_product > 1.0e-9_dbl) &
                  WRITE (*,*) "vectors ",i," and ",j," are not orthogonal:", dot_product
       END DO
    END DO
#endif


    !-----------------------------------------!
    ! calc the matrix T = transp(kv)*S*OP(kv) !
    !-----------------------------------------!
    ALLOCATE(T(1), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, "allocation of T failed", glob_env)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct, glob_env%para_env, glob_env%blacs_env, &
         n_kv, n_kv, &
         optimal_blacs_row_block_size, optimal_blacs_col_block_size)
    CALL cp_fm_create(T(1)%matrix, fm_struct, routineP//"T-"//TRIM(cp_to_string(i)))
    NULLIFY(fm_struct)
    ALLOCATE(Tev(1), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, "allocation of T failed", glob_env)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct, glob_env%para_env, glob_env%blacs_env, &
         n_kv, n_kv, &
         optimal_blacs_row_block_size, optimal_blacs_col_block_size)
    CALL cp_fm_create(Tev(1)%matrix, fm_struct, routineP//"T-"//TRIM(cp_to_string(i)))
    NULLIFY(fm_struct)
    DO col=1, n_kv
       IF (col < n_kv) THEN ! was saved in opkv
          DO spin=1, n_spins
             index = n_spins*(col-1) + spin
             CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
                  opkv(index)%matrix, &
                  R(spin)%matrix, &
                  p_env%n_mo(spin), &
                  glob_env%para_env, 1.0_dbl, 0.0_dbl)
          END DO
       ELSE
          DO spin=1, n_spins
             index=n_spins*(col-1)+spin
             CALL copy_blacs_to_blacs_matrix(kv(index)%matrix, X(spin)%matrix)
          END DO
          CALL apply_op(X, R, t_env, p_env, qs_env, glob_env) ! result in X
          DO spin=1, n_spins
             CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
                  X(spin)%matrix, &
                  R(spin)%matrix, &
                  p_env%n_mo(spin), &
                  glob_env%para_env, 1.0_dbl, 0.0_dbl)
          END DO
       END IF
       DO row=1, n_kv
          dot_product = 0.0_dbl
          DO spin=1, n_spins
             index=n_spins*(row-1)+spin
             CALL cp_fm_trace(kv(index)%matrix, R(spin)%matrix, tmp)
             dot_product = dot_product + tmp
          END DO
          CALL cp_fm_set_element(T(1)%matrix, row, col, dot_product)
       END DO
    END DO


    !---------------------------!
    ! calc the eigenvalues of T !
    !---------------------------!
    CALL blacs_syevd(T(1)%matrix, Tev(1)%matrix, eigenvals)


    !------------------------!
    ! store results in t_env !
    !------------------------!
    WRITE (*,*)
    DO i=1, n_kv
       WRITE (*,'(1X,F20.10)') eigenvals(i)
    END DO
    WRITE (*,*)
    ! remains to be done

    !----------!
    ! clean up !
    !----------!
    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, X)
    CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools, R)
    DO i=1, n_spins*n_kv
       CALL cp_fm_release(kv(i)%matrix)
       CALL cp_fm_release(opkv(i)%matrix)
    END DO
    CALL cp_fm_release(T(1)%matrix)
    CALL cp_fm_release(Tev(1)%matrix)
    DEALLOCATE(X, R, kv, opkv, T, Tev, eigenvals)
    NULLIFY(X, R, kv, opkv, T, Tev)
    
  END SUBROUTINE tddfpt_es_lanczos

! *****************************************************************************
  
  ! X is changed with this routine
  SUBROUTINE apply_op(X, R, t_env, p_env, qs_env, glob_env)

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X 
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: R 
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


!!    CALL p_preortho(p_env, qs_env, X) ! make X orthogonal to psi0

    !------------!
    ! R = HX-SXL !
    !------------!
    CALL p_op_l1(p_env, qs_env, X, R)  ! acts on both spins

#if 1
    !-----------------!
    ! calc P1 and     !
    ! R = R + K(P1)*C !
    !-----------------!
    DO spin=1, n_spins                 
       CALL set_matrix(p_env%p1(spin)%matrix, value=0.0_dbl) ! optimize?
       CALL cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
            matrix_v=p_env%psi0d(spin)%matrix,&
            matrix_g=X(spin)%matrix,&
            ncol=p_env%n_mo(spin))            
    END DO
    CALL p_op_l2(p_env, qs_env, p_env%p1, R, alpha=1.0_dbl, beta=1.0_dbl) ! tam-dancoff
#endif

    !----------------------!
    ! orthogonalize R to C !
    !----------------------!
    CALL p_postortho(p_env, qs_env, R)

    !----------!
    ! X=invS*R !
    !----------!
    DO spin=1, n_spins
       CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                       1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
                       0.0_dbl, X(spin)%matrix)
    END DO

#if 0
    CALL p_op_l1(p_env, qs_env, X, R) ! acts on both spins

    !----------!
    ! X=invS*R !
    !----------!
    DO spin=1, n_spins
       CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                       1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
                       0.0_dbl, X(spin)%matrix)
    END DO
#endif

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

  ! X       : the vector to reorthogonalize
  ! kv      : set of vectors to which to reorthogonalize
  ! invS    : metric of the vectorspace
  ! m_tmp   : temporary matrix
  ! j       : up to which vector of kv to reorthogonalize
  ! n_spins : number of spins to treat

  SUBROUTINE reorthogonalize(X, kv, S, m_tmp, j, n_spins, glob_env)

    IMPLICIT NONE

    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X, kv, m_tmp
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: S
    INTEGER, INTENT(IN) :: j, n_spins
    TYPE(global_environment_type), INTENT(in) :: glob_env

    INTEGER :: i, index, spin
    REAL(DBL) :: tmp, dot_product

    IF (j<1) RETURN
    DO i=1, j
!!    DO i=1, j-1
       dot_product = 0.0_dbl
       DO spin=1, n_spins
          CALL cp_sm_fm_multiply(S(spin)%matrix, X(spin)%matrix, m_tmp(spin)%matrix, &
                                 X(spin)%matrix%matrix_struct%ncol_global, &
                                 glob_env%para_env, 1.0_dbl, 0.0_dbl)
          index = n_spins*(i-1)+spin
          CALL cp_fm_trace(m_tmp(spin)%matrix, kv(index)%matrix, tmp)
          dot_product = dot_product + tmp
       END DO
       DO spin=1, n_spins
          index = n_spins*(i-1)+spin
          CALL cp_fm_scale_and_add(1.0_dbl, X(spin)%matrix, -1.0_dbl*dot_product, kv(index)%matrix)
       END DO
    END DO
  END SUBROUTINE reorthogonalize

! *****************************************************************************

END MODULE qs_tddfpt_eigensolver
