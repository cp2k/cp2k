!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Front-End for any PAO parametrization
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_param
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_get_block_p, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_type
  USE cp_log_handling,                 ONLY: cp_to_string
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE kinds,                           ONLY: dp
  USE pao_input,                       ONLY: pao_exp_param
  USE pao_param_exp,                   ONLY: pao_calc_U_exp,&
                                             pao_param_finalize_exp,&
                                             pao_param_init_exp
  USE pao_param_linpot,                ONLY: pao_calc_U_linpot,&
                                             pao_param_finalize_linpot,&
                                             pao_param_init_linpot
  USE pao_types,                       ONLY: pao_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param'

  PUBLIC :: pao_param_init, pao_param_finalize, pao_calc_U, pao_calc_grad

CONTAINS


! *****************************************************************************
!> \brief Initialize PAO parametrization
!> \param pao ...
!> \param qs_env ...
!> \param reuse_matrix_X ...
! *****************************************************************************
  SUBROUTINE pao_param_init(pao, qs_env, reuse_matrix_X)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: reuse_matrix_X

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init', &
      routineP = moduleN//':'//routineN

!IF(pao%istep>=0) RETURN ! was there a previous pao-run, eg. during MD ?

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        CALL pao_param_init_exp(pao, qs_env, reuse_matrix_X)
      CASE DEFAULT
        CALL pao_param_init_linpot(pao, qs_env, reuse_matrix_X)
    END SELECT

  END SUBROUTINE pao_param_init


! *****************************************************************************
!> \brief Finalize PAO parametrization
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_param_finalize(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize', &
      routineP = moduleN//':'//routineN

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        CALL pao_param_finalize_exp(pao)
      CASE DEFAULT
        CALL pao_param_finalize_linpot(pao)
    END SELECT

  END SUBROUTINE pao_param_finalize


! *****************************************************************************
!> \brief Calculate new matrix U
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_calc_U(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, iatom, n1, n2
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_precon, block_pX, &
                                                block_U, block_X
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
       IF(arow /= acol) CPABORT("encountered off-diagonal block")
       iatom = arow
       n1 = SIZE(block_X,1)
       n2 = SIZE(block_X,2)

       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=iatom, col=iatom, block=block_U, found=found)
       CPASSERT(ASSOCIATED(block_U))
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_precon, row=iatom, col=iatom, block=block_precon, found=found)
       CPASSERT(ASSOCIATED(block_precon))

       ! apply preconditioner
       ALLOCATE(block_pX(n1,n2))
       block_pX = MATMUL(block_precon, block_X)

       CALL pao_calc_U_low(pao, iatom, block_pX, block_U)

       DEALLOCATE(block_pX)
    END DO

    CALL cp_dbcsr_iterator_stop(iter)

    CALL assert_unitary(pao, pao%matrix_U, pao%matrix_Y)

  END SUBROUTINE pao_calc_U


! *****************************************************************************
!> \brief Helper routien for pao_calc_U
!> \param pao ...
!> \param iatom ...
!> \param block_X ...
!> \param block_U ...
! *****************************************************************************
 SUBROUTINE pao_calc_U_low(pao, iatom, block_X, block_U)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER, INTENT(IN)                      :: iatom
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, block_U

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_low', &
      routineP = moduleN//':'//routineN

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        CALL pao_calc_U_exp(pao, iatom, block_X, block_U)
      CASE DEFAULT
        CALL pao_calc_U_linpot(pao, iatom, block_X, block_U)
    END SELECT

  END SUBROUTINE pao_calc_U_low


! *****************************************************************************
!> \brief Calculate the gradient G = dU/dX
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_calc_grad(qs_env,ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, iatom, &
                                                n1, n2
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_G, block_M, block_pG, &
                                                block_precon, block_pX, &
                                                block_U, block_X
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_M
    TYPE(pao_env_type), POINTER              :: pao

    CALL timeset(routineN,handle)

    pao => ls_scf_env%pao_env

    CALL pao_calc_grad_matrix_M(qs_env,ls_scf_env, matrix_M)

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
       IF(arow /= acol) CPABORT("encountered off-diagonal block")
       iatom = arow
       n1 = SIZE(block_X,1)
       n2 = SIZE(block_X,2)

       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=arow, col=acol, block=block_U, found=found)
       CPASSERT(ASSOCIATED(block_U))
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=arow, col=acol, block=block_G, found=found)
       CPASSERT(ASSOCIATED(block_G))
       CALL cp_dbcsr_get_block_p(matrix=matrix_M, row=arow, col=acol, block=block_M, found=found)
       CPASSERT(ASSOCIATED(block_M))
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_precon, row=iatom, col=iatom, block=block_precon, found=found)
       CPASSERT(ASSOCIATED(block_precon))

       ! apply preconditioner to block_X
       ALLOCATE(block_pX(n1,n2), block_pG(n1,n2))
       block_pX = MATMUL(block_precon, block_X)

       SELECT CASE(pao%parameterization)
         CASE(pao_exp_param)
           CALL pao_calc_U_exp(pao, iatom, block_pX, block_U, block_M, block_pG)
         CASE DEFAULT
           CALL pao_calc_U_linpot(pao, iatom, block_pX, block_U, block_M, block_pG)
       END SELECT

       CALL check_gradient(pao, iatom, block_pX, block_pG, block_M)

       ! apply preconditioner to block_G
       block_G = MATMUL(block_precon, block_pG)

       DEALLOCATE(block_pX, block_pG)
    END DO

    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_release(matrix_M)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_grad


! *****************************************************************************
!> \brief Debugging routien for checking the analytic gradient.
!> \param pao ...
!> \param iatom ...
!> \param block_X ...
!> \param block_G ...
!> \param block_M ...
! *****************************************************************************
  SUBROUTINE check_gradient(pao, iatom, block_X, block_G, block_M)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER, INTENT(IN)                      :: iatom
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, block_G, block_M

    CHARACTER(len=*), PARAMETER :: routineN = 'check_gradient', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, M, N
    REAL(dp)                                 :: delta, eps, Gij_num, symm
    REAL(dp), DIMENSION(:, :), POINTER       :: dU1, dU2, dUdX, dX

    IF(pao%check_gradient_tol<0.0_dp) RETURN ! no checking

    N = SIZE(block_M, 1);  M = SIZE(block_M, 2)
    ALLOCATE(dUdX(N,M), dU1(N,N), dU2(N,N))
    ALLOCATE(dX(SIZE(block_X,1), SIZE(block_X,2)))

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        symm = -1.0_dp ! anti-symmetric
      CASE DEFAULT
        symm = 0.0_dp  ! no symmetry
    END SELECT

    eps = 0.00001_dp

    DO i=1, SIZE(block_X,1)
       DO j=1, SIZE(block_X,2)
          dX = block_X
          dX(i,j) = dX(i,j) + eps
          IF(symm/=0.0) dX(j,i) = dX(j,i) + symm*eps
          CALL pao_calc_U_low(pao, iatom, dX, dU1)

          dX = block_X
          dX(i,j) = dX(i,j) - eps
          IF(symm/=0.0) dX(j,i) = dX(j,i) - symm*eps
          CALL pao_calc_U_low(pao, iatom, dX, dU2)

          dUdX = (dU1 - dU2) / (2.0_dp*eps)
          Gij_num = SUM( block_M * dUdX )
          delta = ABS(Gij_num - block_G(i,j))
          IF(delta > pao%check_gradient_tol) THEN
             CALL cp_abort(__LOCATION__,&
                   "Analytic and numeric gradients differ too much:"//cp_to_string(delta))
          ENDIF
       ENDDO
    ENDDO

    DEALLOCATE(dUdX, dU1, dU2, dX)

  END SUBROUTINE check_gradient


! *****************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M ...
! *****************************************************************************
   SUBROUTINE pao_calc_grad_matrix_M(qs_env,ls_scf_env,matrix_M)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    TYPE(cp_dbcsr_type)                      :: matrix_M

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad_matrix_M', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_dbcsr_type) :: matrix_HB, matrix_HBP, matrix_HP, matrix_M_pre, &
      matrix_PHP, matrix_SB, matrix_SBPHP, matrix_tmp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env

    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control,&
                    rho=rho)
    CALL qs_rho_get(rho, rho_ao=rho_ao)
    nspin = dft_control%nspins

    IF(nspin /= 1) CPABORT("open shell not yet implemented")
    !TODO: handle openshell case properly

    CALL cp_dbcsr_init(matrix_HB)
    CALL cp_dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B, 0.0_dp,&
                           matrix_HB)

    CALL cp_dbcsr_init(matrix_HBP)
    CALL cp_dbcsr_create(matrix_HBP, template=ls_mstruct%matrix_B)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_HB, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_HBP)

    CALL cp_dbcsr_init(matrix_HP)
    CALL cp_dbcsr_create(matrix_HP, template=ls_scf_env%matrix_p(1))
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_ks(1), ls_scf_env%matrix_p(1), 0.0_dp,&
                           matrix_HP)

    CALL cp_dbcsr_init(matrix_PHP)
    CALL cp_dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_p(1))
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), matrix_HP, 0.0_dp,&
                           matrix_PHP)

    CALL cp_dbcsr_init(matrix_SB)
    CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B, 0.0_dp,&
                           matrix_SB)

    CALL cp_dbcsr_init(matrix_SBPHP)
    CALL cp_dbcsr_create(matrix_SBPHP, template=ls_mstruct%matrix_B)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_SB, matrix_PHP, 0.0_dp,&
                           matrix_SBPHP)

    IF(SIZE(ls_scf_env%matrix_p)==1) CALL cp_dbcsr_scale(matrix_SBPHP, 0.5_dp)

    CALL cp_dbcsr_init(matrix_M_pre)
    CALL cp_dbcsr_copy(matrix_M_pre, matrix_HBP)
    CALL cp_dbcsr_add(matrix_M_pre, matrix_SBPHP, 2.0_dp, -2.0_dp)

    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp, template=matrix_M_pre)
    CALL cp_dbcsr_multiply("T", "N", 1.0_dp, pao%matrix_N, matrix_M_pre, 0.0_dp,&
                           matrix_tmp)

    CALL cp_dbcsr_init(matrix_M)
    CALL cp_dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, pao%matrix_Y, 0.0_dp,&
                           matrix_M)

    CALL cp_dbcsr_release(matrix_HB)
    CALL cp_dbcsr_release(matrix_HBP)
    CALL cp_dbcsr_release(matrix_HP)
    CALL cp_dbcsr_release(matrix_PHP)
    CALL cp_dbcsr_release(matrix_SB)
    CALL cp_dbcsr_release(matrix_SBPHP)
    CALL cp_dbcsr_release(matrix_M_pre)
    CALL cp_dbcsr_release(matrix_tmp)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_grad_matrix_M


! *****************************************************************************
!> \brief Debugging routine, check unitaryness of U
!> \param pao ...
!> \param matrix_test ...
!> \param matrix_Y ...
! *****************************************************************************
  SUBROUTINE assert_unitary(pao, matrix_test, matrix_Y)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(cp_dbcsr_type)                      :: matrix_test, matrix_Y

    CHARACTER(len=*), PARAMETER :: routineN = 'assert_unitary', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, i, M, N
    LOGICAL                                  :: found
    REAL(dp)                                 :: delta
    REAL(dp), DIMENSION(:, :), POINTER       :: block_test, block_Y, tmp1, &
                                                tmp2
    TYPE(cp_dbcsr_iterator)                  :: iter

    IF(pao%check_unitary_tol<0.0_dp) RETURN ! no checking

    CALL cp_dbcsr_iterator_start(iter, matrix_test)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_test)
       CALL cp_dbcsr_get_block_p(matrix=matrix_Y, row=arow, col=acol, block=block_Y, found=found)
       CPASSERT(ASSOCIATED(block_Y))
       N = SIZE(block_Y, 1)
       M = SIZE(block_Y, 2)
       ALLOCATE(tmp1(N,M), tmp2(M, M))

       ! we only need the upper left "PAO-corner" to be unitary
       tmp1 = MATMUL(block_test, block_Y)
       tmp2 = MATMUL(TRANSPOSE(tmp1), tmp1)
       DO i=1, M
          tmp2(i,i) = tmp2(i,i) - 1.0_dp
       ENDDO

       delta = MAXVAL(ABS(tmp2))
       IF(delta > pao%check_unitary_tol) THEN
          CALL cp_abort(__LOCATION__,&
             "Found bad unitaryness:"//cp_to_string(delta))
       ENDIF

       DEALLOCATE(tmp1, tmp2)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

  END SUBROUTINE assert_unitary

END MODULE pao_param
