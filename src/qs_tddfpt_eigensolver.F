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

  USE cp_error_handling,               ONLY: cp_simple_assert
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_trace,&
                                             cp_sm_fm_multiply
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_fm_vect,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_type,&
                                             optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_fm_utils,                     ONLY: copy_fm_to_rep_vec,&
                                             copy_rep_vec_to_fm
  USE cp_log_handling,                 ONLY: cp_to_string
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_arpack_utils,                 ONLY: arpack_confirm
  USE qs_blacs,                        ONLY: blacs_syevd,&
                                             copy_blacs_to_blacs_matrix,&
                                             copy_real_matrix_to_full_matrix,&
                                             cp_sm_plus_fm_fm_t
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_p_env_methods,                ONLY: p_op_l1,&
                                             p_op_l2,&
                                             p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_tddfpt_types
  USE qs_tddfpt_utils
  USE sparse_matrix_types,             ONLY: add_all_real_matrix_blocks,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  TYPE cp_fm_struct_p_type
     TYPE(cp_fm_struct_type), POINTER :: struct
  END TYPE cp_fm_struct_p_type

  INTERFACE tddfpt_es
     MODULE PROCEDURE tddfpt_es_lanczos
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_eigensolver"
  INTEGER, PARAMETER          :: MAX_ITER_GLOBAL = 20
  REAL(kind=dbl), PARAMETER   :: TOLERANCE = 1.0E-4

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
                 n_ev, n_rows, &
                 n_kv, n_kv_max, &
                 max_iter, &
                 n_reortho, &
                 row, col, &
                 i, j, index, &
                 stat, n_kv_fac, &
                 iteration
    REAL(dbl) :: tol, &
                 beta, &
                 tmp, &
                 dot_product, &
                 scale, difference
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE               :: eigenvals
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER       :: ao_mo_fm_pools
    TYPE(cp_full_matrix_p_type)                          :: T, Tev           
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER   :: X, R
    TYPE(cp_full_matrix_p_type), DIMENSION(:,:), POINTER :: kv, opkv
    TYPE(cp_fm_struct_p_type), DIMENSION(:), POINTER     :: kv_fm_struct
    TYPE(cp_fm_struct_type), POINTER                     :: T_fm_struct, &
                                                            Tev_fm_struct

    CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_es_lanczos", &
                                   routineP = moduleN//"/"//routineN

    NULLIFY(ao_mo_fm_pools, X, R, kv)
    n_ev      = qs_env%dft_control%tddfpt_control%n_ev
    n_spins   = qs_env%dft_control%nspins
    tol       = 1.0E-10
    n_reortho = qs_env%dft_control%tddfpt_control%n_reortho
    iteration = 1

    !-----------------------------------------!
    ! determine the size of the problem       !
    ! and how many krylov space vetors to use ! 
    !-----------------------------------------!
    n_rows = 0
    DO spin=1, n_spins
       n_rows = n_rows + p_env%n_ao(spin)*p_env%n_mo(spin)
    END DO
    n_kv = MAX(3, n_ev) ! how many vectors per iteration
    max_iter = MIN(( n_rows / n_kv ) * n_kv, MAX_ITER_GLOBAL)
    n_kv_max = max_iter * n_kv

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    ALLOCATE(eigenvals(n_kv_max, max_iter))
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, X, &
                                 name=routineP//":X")
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, R, &
                                 name=routineP//":R")
    ALLOCATE(kv(1:n_kv_max,n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of kv failed", glob_env)
    ALLOCATE(opkv(n_kv_max,n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of opkv failed", glob_env)
    ALLOCATE(kv_fm_struct(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of kv_fm_struct failed", &
                                     glob_env)
    DO spin=1, n_spins
       NULLIFY (kv_fm_struct(spin)%struct)
       CALL cp_fm_struct_create(kv_fm_struct(spin)%struct, glob_env%para_env, &
                                glob_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_mo(spin), &
                                optimal_blacs_row_block_size, &
                                optimal_blacs_col_block_size)
    END DO
    CALL allocate_krylov_vectors(kv, opkv, iteration, n_kv, n_spins, kv_fm_struct)

    DO spin=1, n_spins
       CALL copy_blacs_to_blacs_matrix(t_env%X(1)%spins(spin)%matrix, &
                                       X(spin)%matrix)
    END DO

    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, &
              '(2X,A,2X,A,T69,A)') &
              "Step", "Vectors", "Convergence"
       WRITE (glob_env%scr, &
              '(2X,A)') &
              "-----------------------------------------------------------------------------"
    END IF

    restart: DO 

       krylov: DO j=1, n_kv
       
          index = (iteration-1)*n_kv+j
          
          !---------------------------------------!
          ! apply the operator to working vec and !
          ! store the result for later use (T)    !
          !---------------------------------------!
          IF (j /= 1) THEN
             CALL apply_op(X, R, t_env, p_env, qs_env, glob_env) ! result in X
             DO spin=1, n_spins
                CALL copy_blacs_to_blacs_matrix(X(spin)%matrix, &
                                                opkv(index-1,spin)%matrix)
             END DO
          END IF

          DO i=1, n_reortho
             CALL reorthogonalize(X, kv, qs_env%s, R, index, n_spins, glob_env)
          END DO

          CALL calc_norm(norm=beta, x=X, tmp_vec=R, &
                         qs_env=qs_env, glob_env=glob_env)

          !-------------!
          ! normalize X ! 
          !-------------!
          IF (ABS(beta) > tol) THEN
             DO spin=1, n_spins
                CALL cp_fm_scale(SQRT(1._dbl/beta), X(spin)%matrix)
             END DO
          ELSE 
             IF (index==1) THEN
                CALL stop_program(routineP, "Initial guess to close to zero.")
             ELSE
                ! already converged
                index = index-1
                n_kv_max = index
                IF (glob_env%ionode) THEN
                   WRITE (*,*)
                   WRITE (*,*) "Converged with ", index," krylov vectors!"
                   WRITE (*,*)
                END IF
                EXIT krylov
             END IF
          END IF

          !---------------------------!
          ! write current vector back !
          !---------------------------!
          DO spin=1, n_spins
             CALL copy_blacs_to_blacs_matrix(X(spin)%matrix, &
                                             kv(index,spin)%matrix)
          END DO

       END DO krylov

       !---------------------------------------!
       ! compute the last of the Op*kv vectors !
       !---------------------------------------!
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(kv(index,spin)%matrix, &
                                          X(spin)%matrix)
       END DO
       CALL apply_op(X, R, t_env, p_env, qs_env, glob_env) ! result in X
       DO spin=1, n_spins
          CALL copy_blacs_to_blacs_matrix(X(spin)%matrix, &
                                         opkv(index,spin)%matrix)
       END DO
       

       !--------------------------------------------!
       ! allocate memory for the T and Tev matrices !
       !--------------------------------------------!
       IF (ASSOCIATED(T%matrix)) THEN
          CALL cp_fm_release(T%matrix)
       END IF
       NULLIFY(T_fm_struct)
       CALL cp_fm_struct_create(T_fm_struct, &
                                glob_env%para_env, &
                                glob_env%blacs_env, &
                                index, index, &
                                optimal_blacs_row_block_size, &
                                optimal_blacs_col_block_size)
       CALL cp_fm_create(T%matrix, &
                         T_fm_struct, &
                         routineP//"T")
       IF (ASSOCIATED(Tev%matrix)) THEN
          CALL cp_fm_release(Tev%matrix)
       END IF
       NULLIFY(Tev_fm_struct)
       CALL cp_fm_struct_create(Tev_fm_struct, &
                                glob_env%para_env, &
                                glob_env%blacs_env, &
                                index, index, &
                                optimal_blacs_row_block_size, &
                                optimal_blacs_col_block_size)
       CALL cp_fm_create(Tev%matrix, &
                         Tev_fm_struct, &
                         routineP//"Tev")


       !-----------------------------------------!
       ! calc the matrix T = transp(kv)*S*OP(kv) !
       !-----------------------------------------!
       DO col=1, index
          DO spin=1, n_spins
             CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
                                    opkv(col,spin)%matrix, &
                                    R(spin)%matrix, &
                                    p_env%n_mo(spin), &
                                    glob_env%para_env, 1.0_dbl, 0.0_dbl)
          END DO
          DO row=1, index
             dot_product = 0.0_dbl
             DO spin=1, n_spins
                CALL cp_fm_trace(kv(row,spin)%matrix, R(spin)%matrix, tmp)
                dot_product = dot_product + tmp
             END DO
             CALL cp_fm_set_element(T%matrix, row, col, dot_product)
          END DO
       END DO


       !---------------------------!
       ! calc the eigenvalues of T !
       !---------------------------!
       CALL blacs_syevd(T%matrix, Tev%matrix, eigenvals(:,iteration))


       IF (iteration /= 1) THEN
          difference = ABS(eigenvals(n_ev,iteration)-eigenvals(n_ev,iteration-1))
       ELSE
          difference = 9.9999E99
       END IF
       IF (glob_env%ionode) &
       WRITE (glob_env%scr, &
              '(2X,I4,2X,I7,T69,ES11.4)') &
              iteration, index, difference
       IF (difference < TOLERANCE) THEN
          EXIT restart
       END IF

       IF (index == n_kv_max) THEN
          WRITE (*,*) "MAXIMUM NUMBER OF KRYLOV SPACE VECTORS REACHED"
          EXIT restart
       END IF
       
       iteration = iteration + 1

       !--------------------------------------!
       ! allocate matrizes for next iteration !
       !--------------------------------------!
       CALL allocate_krylov_vectors(kv, opkv, iteration, n_kv, n_spins, kv_fm_struct)

    END DO restart


    !-------------------!
    ! Print the results !
    !-------------------!
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, *)
       WRITE (glob_env%scr, *)
       WRITE (glob_env%scr, '(1X,A,T19,A,T7,A)') &
            "Nr.", "Eigenvalue (hartree)", "Eigenvalue (eV)"
       DO i=1, MIN(n_ev+1,n_kv_max)
          WRITE (glob_env%scr, '(1X,I3,T20,F19.10,2X,F19.10)') &
               i, eigenvals(i, iteration), eigenvals(i, iteration)*27.211_dbl  
          ! this should use the physcons library !
       END DO
       WRITE (glob_env%scr, *)
    END IF

    !----------!
    ! clean up !
    !----------!
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, X)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, R)
    DO spin=1, n_spins
       DO i=1, size(kv,spin)
          CALL cp_fm_release(kv(i,spin)%matrix)
          CALL cp_fm_release(opkv(i,spin)%matrix)
       END DO
    END DO
    IF (ASSOCIATED(T%matrix)) CALL cp_fm_release(T%matrix)
    IF (ASSOCIATED(Tev%matrix)) CALL cp_fm_release(Tev%matrix)
    DEALLOCATE(kv, opkv, eigenvals)
    NULLIFY(X, R, kv, opkv)
    
  END SUBROUTINE tddfpt_es_lanczos

! *****************************************************************************

  ! X        : on input the vector on which to apply the op
  !            on output the result
  ! R        : a temporary vector
  ! t_env    : td-dft environment (mainly control information)
  ! p_env    : perturbation environment (variables)
  !            both of these carry info for the tddfpt calculation
  ! qs_env   : info about a quickstep ground state calculation
  ! glob_env : for printing and communication

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
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: p1
    INTEGER                                            :: nrows, ncols, &
                                                          nblock_rows, nblock_cols
    INTEGER, DIMENSION(:), POINTER                     :: first_row, last_row, &
                                                          first_col, last_col



    CHARACTER(LEN=*), PARAMETER :: routineN = "apply_op", &
                                   routineP = moduleN//"/"//routineN

    counter = counter + 1
    
    n_spins = qs_env%dft_control%nspins

    !------------!
    ! R = HX-SXL !
    !------------!
    CALL p_op_l1(p_env, qs_env, X, R)  ! acts on both spins, result in R

    !-----------------!
    ! calc P1 and     !
    ! R = R + K(P1)*C !
    !-----------------!
    IF (qs_env%dft_control%tddfpt_control%do_kernel) THEN
       DO spin=1, n_spins
          CALL set_matrix(p_env%p1(spin)%matrix, value=0.0_dbl) ! optimize?
          CALL cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=p_env%psi0d(spin)%matrix,&
                                  matrix_g=X(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dbl) 
          CALL cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=X(spin)%matrix,&
                                  matrix_g=p_env%psi0d(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dbl) 
       END DO
       CALL p_op_l2(p_env, qs_env, p_env%p1, R, &
                    alpha=2.0_dbl, beta=1.0_dbl) ! R = beta*R + alpha*K(P1)*C
    END IF

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

  ! X        : the vector to reorthogonalize
  ! kv       : set of vectors to which to reorthogonalize
  ! S        : metric of the vectorspaceA
  ! tmp_mat    : temporary matrix
  ! j        : up to which vector of kv to reorthogonalize
  ! n_spins  : number of spins to treat
  ! glob_env : for communication

  SUBROUTINE reorthogonalize(X, kv, S, tmp_mat, j, n_spins, glob_env)

    IMPLICIT NONE

    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: X, tmp_mat
    TYPE(cp_full_matrix_p_type), DIMENSION(:,:), POINTER :: kv
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: S
    INTEGER, INTENT(IN) :: j, n_spins
    TYPE(global_environment_type), INTENT(in) :: glob_env

    INTEGER :: i, index, spin
    REAL(DBL) :: tmp, dot_product

    IF (j<=1) RETURN
    DO i=1, j
       dot_product = 0.0_dbl
       DO spin=1, n_spins
          CALL cp_sm_fm_multiply(S(spin)%matrix, X(spin)%matrix, &
                                 tmp_mat(spin)%matrix, &
                                 X(spin)%matrix%matrix_struct%ncol_global, &
                                 glob_env%para_env, 1.0_dbl, 0.0_dbl)
          CALL cp_fm_trace(tmp_mat(spin)%matrix, kv(i,spin)%matrix, tmp)
          dot_product = dot_product + tmp
       END DO
       DO spin=1, n_spins
          CALL cp_fm_scale_and_add(1.0_dbl, X(spin)%matrix, &
                                   -1.0_dbl*dot_product, kv(i,spin)%matrix)
       END DO
    END DO
  END SUBROUTINE reorthogonalize

! *****************************************************************************

  SUBROUTINE allocate_krylov_vectors(kv, opkv, &
                                     iteration, n_ev, n_spins, kv_fm_struct)

    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_p_type), DIMENSION(:,:), POINTER :: kv, opkv
    INTEGER, INTENT(IN)                                  :: iteration, &
                                                            n_ev, &
                                                            n_spins
    TYPE(cp_fm_struct_p_type), DIMENSION(:), POINTER     :: kv_fm_struct

    ! locals
    INTEGER                     :: spin, i, index
    CHARACTER(LEN=80)           :: mat_name
    CHARACTER(LEN=*), PARAMETER :: routineN = "allocate_krylov_vectors", &
                                   routineP = moduleN//"/"//routineN

    DO spin=1, n_spins
       DO i=1, n_ev
          index = (iteration-1)*n_ev+i
          NULLIFY(kv(index,spin)%matrix)
          mat_name = routineP//"kv-"//TRIM(cp_to_string(index))&
                 //","//TRIM(cp_to_string(spin))
          CALL cp_fm_create(kv(index,spin)%matrix,  &
               kv_fm_struct(spin)%struct, mat_name)
          IF (.NOT.ASSOCIATED(kv(index,spin)%matrix)) &
               CALL stop_program(routineP, &
               "Could not alloacte "//TRIM(mat_name)//".")
          mat_name = routineP//"opkv-"//TRIM(cp_to_string(index))&
                 //","//TRIM(cp_to_string(spin))
          CALL cp_fm_create(opkv(index,spin)%matrix, &
               kv_fm_struct(spin)%struct, mat_name)
          IF (.NOT.ASSOCIATED(opkv(index,spin)%matrix)) &
               CALL stop_program(routineP, &
               "Could not alloacte "//TRIM(mat_name)//".")
       END DO
    END DO

  END SUBROUTINE allocate_krylov_vectors

! *****************************************************************************

  !---------------------------------------!
  ! x must not be changed in this routine !
  ! tmp_vec may be changed                !
  !---------------------------------------!
  SUBROUTINE calc_norm(norm, x, tmp_vec, qs_env, glob_env)

    IMPLICIT NONE

    REAL(KIND=dbl), INTENT(OUT)                        :: norm
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: x, tmp_vec
    TYPE(qs_environment_type)                          :: qs_env
    TYPE(global_environment_type)                      :: glob_env

    INTEGER        :: spin, n_spins
    REAL(KIND=dbl) :: spin_norm

    n_spins = SIZE(x)

    !-------------------------!
    ! calc beta=transp(X)*S*X !
    !-------------------------!
    norm = 0.0_dbl
    DO spin=1, n_spins
       CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, x(spin)%matrix, &
            tmp_vec(spin)%matrix, &
            x(spin)%matrix%matrix_struct%ncol_global, &
            glob_env%para_env, 1.0_dbl, 0.0_dbl)
       CALL cp_fm_trace(tmp_vec(spin)%matrix, x(spin)%matrix, spin_norm)
       norm = norm + spin_norm
    END DO

  END SUBROUTINE calc_norm

! *****************************************************************************

END MODULE qs_tddfpt_eigensolver
