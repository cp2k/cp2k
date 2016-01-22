!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Front-End for any PAO parametrization
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_param
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_get_block_p, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_release, &
       cp_dbcsr_type, dbcsr_distribution_mp, dbcsr_mp_group
  USE cp_log_handling,                 ONLY: cp_to_string
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max,&
                                             mp_min
  USE pao_input,                       ONLY: pao_exp_param
  USE pao_param_exp,                   ONLY: pao_calc_U_exp,&
                                             pao_param_finalize_exp,&
                                             pao_param_init_exp
  USE pao_param_linpot,                ONLY: pao_calc_U_linpot,&
                                             pao_param_finalize_linpot,&
                                             pao_param_init_linpot
  USE pao_types,                       ONLY: pao_env_type
  USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param'

  PUBLIC :: pao_update_AB
  PUBLIC :: pao_param_init, pao_param_finalize, pao_calc_U, pao_calc_grad

CONTAINS

! *****************************************************************************
!> \brief Takes current matrix_X and recalculates derived matrices U, A, and B.
!> \param pao ...
!> \param ls_mstruct ...
! *****************************************************************************
  SUBROUTINE pao_update_AB(pao, ls_mstruct)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(ls_mstruct_type)                    :: ls_mstruct

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_update_AB', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)

    CALL pao_calc_U(pao) !update matrix_U = Function of matrix_X

    !update matrix_A
    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp, template=pao%matrix_U)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N_inv, pao%matrix_U, 0.0_dp,&
                           matrix_tmp)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, pao%matrix_Y,&
                           0.0_dp, ls_mstruct%matrix_A)

    ! update matrix_B
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, pao%matrix_U, 0.0_dp,&
                           matrix_tmp)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, pao%matrix_Y,&
                           0.0_dp, ls_mstruct%matrix_B)

    CALL cp_dbcsr_release(matrix_tmp)

    CALL timestop(handle)
  END SUBROUTINE pao_update_AB


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

    INTEGER                                  :: acol, arow, group, handle, &
                                                iatom, n1, n2
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_U, block_X
    REAL(KIND=dp)                            :: min_gap
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    min_gap = 1.0_dp

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
       IF(arow /= acol) CPABORT("encountered off-diagonal block")
       iatom = arow
       n1 = SIZE(block_X,1)
       n2 = SIZE(block_X,2)

       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=iatom, col=iatom, block=block_U, found=found)
       CPASSERT(ASSOCIATED(block_U))

       CALL pao_calc_U_low(pao, iatom, min_gap, block_X, block_U)
    END DO

    CALL cp_dbcsr_iterator_stop(iter)

    group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(pao%matrix_U)))
    CALL mp_min(min_gap,group)
    IF(pao%iw>0) WRITE(pao%iw,*) "PAO| min_gap:", min_gap

    CALL pao_assert_unitary(pao, pao%matrix_U, pao%matrix_Y)
    CALL timestop(handle)
  END SUBROUTINE pao_calc_U


! *****************************************************************************
!> \brief Helper routien for pao_calc_U
!> \param pao ...
!> \param iatom ...
!> \param min_gap ...
!> \param block_X ...
!> \param block_U ...
! *****************************************************************************
 SUBROUTINE pao_calc_U_low(pao, iatom, min_gap, block_X, block_U)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER, INTENT(IN)                      :: iatom
    REAL(KIND=dp), INTENT(INOUT)             :: min_gap
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, block_U

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_low', &
      routineP = moduleN//':'//routineN

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        CALL pao_calc_U_exp(pao, iatom, block_X, block_U)
      CASE DEFAULT
        CALL pao_calc_U_linpot(pao, iatom, min_gap, block_X, block_U)
    END SELECT

  END SUBROUTINE pao_calc_U_low


! *****************************************************************************
!> \brief Calculate the gradient G = dU/dX
!> \param pao ...
!> \param matrix_M ...
! *****************************************************************************
  SUBROUTINE pao_calc_grad(pao, matrix_M)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(cp_dbcsr_type)                      :: matrix_M

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, group, handle, &
                                                iatom
    LOGICAL                                  :: found
    REAL(dp)                                 :: min_gap
    REAL(dp), DIMENSION(:, :), POINTER       :: block_G, block_M, block_U, &
                                                block_X
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    min_gap = 1.0_dp

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
       IF(arow /= acol) CPABORT("encountered off-diagonal block")
       iatom = arow

       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=arow, col=acol, block=block_U, found=found)
       CPASSERT(ASSOCIATED(block_U))
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=arow, col=acol, block=block_G, found=found)
       CPASSERT(ASSOCIATED(block_G))
       CALL cp_dbcsr_get_block_p(matrix=matrix_M, row=arow, col=acol, block=block_M, found=found)
       CPASSERT(ASSOCIATED(block_M))

       SELECT CASE(pao%parameterization)
         CASE(pao_exp_param)
           CALL pao_calc_U_exp(pao, iatom, block_X, block_U, block_M, block_G)
         CASE DEFAULT
           CALL pao_calc_U_linpot(pao, iatom, min_gap, block_X, block_U, block_M, block_G)
       END SELECT
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_M)))
    CALL mp_min(min_gap,group)
    IF(pao%iw>0) WRITE(pao%iw,*) "PAO| min_gap:", min_gap

    CALL timestop(handle)
  END SUBROUTINE pao_calc_grad


! *****************************************************************************
!> \brief Debugging routine, check unitaryness of U
!> \param pao ...
!> \param matrix_test ...
!> \param matrix_Y ...
! *****************************************************************************
  SUBROUTINE pao_assert_unitary(pao, matrix_test, matrix_Y)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(cp_dbcsr_type)                      :: matrix_test, matrix_Y

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_assert_unitary', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, group, handle, i, &
                                                M, N
    LOGICAL                                  :: found
    REAL(dp)                                 :: delta_max
    REAL(dp), DIMENSION(:, :), POINTER       :: block_test, block_Y, tmp1, &
                                                tmp2
    TYPE(cp_dbcsr_iterator)                  :: iter

    IF(pao%check_unitary_tol<0.0_dp) RETURN ! no checking

    CALL timeset(routineN,handle)
    delta_max = 0.0_dp

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

       delta_max = MAX(delta_max, MAXVAL(ABS(tmp2)))
       DEALLOCATE(tmp1, tmp2)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_test)))
    CALL mp_max(delta_max, group)
    IF(pao%iw>0) WRITE(pao%iw,*) 'PAO| checked unitaryness, max delta:', delta_max
    IF(delta_max > pao%check_unitary_tol)&
       CPABORT("Found bad unitaryness:"//cp_to_string(delta_max))

     CALL timestop(handle)
  END SUBROUTINE pao_assert_unitary

END MODULE pao_param
