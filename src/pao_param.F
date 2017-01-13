!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Front-End for any PAO parametrization
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param
   USE cp_log_handling,                 ONLY: cp_to_string
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, dbcsr_distribution, &
        dbcsr_distribution_mp, dbcsr_frobenius_norm, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_mp_group, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_max
   USE pao_input,                       ONLY: pao_exp_param,&
                                              pao_fock_param,&
                                              pao_gth_param,&
                                              pao_rotinv_param
   USE pao_param_exp,                   ONLY: pao_calc_U_exp,&
                                              pao_param_count_exp,&
                                              pao_param_finalize_exp,&
                                              pao_param_init_exp,&
                                              pao_param_initguess_exp
   USE pao_param_gth,                   ONLY: pao_calc_U_gth,&
                                              pao_param_count_gth,&
                                              pao_param_finalize_gth,&
                                              pao_param_init_gth,&
                                              pao_param_initguess_gth
   USE pao_param_linpot,                ONLY: pao_calc_U_linpot,&
                                              pao_param_count_linpot,&
                                              pao_param_finalize_linpot,&
                                              pao_param_init_linpot,&
                                              pao_param_initguess_linpot
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param'

   PUBLIC :: pao_update_AB, pao_param_count, pao_param_initial_guess
   PUBLIC :: pao_param_init, pao_param_finalize, pao_calc_U

CONTAINS

! **************************************************************************************************
!> \brief Takes current matrix_X and recalculates derived matrices U, A, and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_mstruct ...
!> \param penalty ...
! **************************************************************************************************
   SUBROUTINE pao_update_AB(pao, qs_env, ls_mstruct, penalty)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_mstruct_type)                              :: ls_mstruct
      REAL(dp), INTENT(OUT), OPTIONAL                    :: penalty

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_update_AB'

      INTEGER                                            :: acol, arow, handle, iatom
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_N, block_N_inv, &
                                                            block_U, block_Y
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_U

      CALL timeset(routineN, handle)

      CALL pao_calc_U(pao, qs_env, penalty=penalty) !update matrix_U = Function of matrix_X

      ! pao%matrix_U uses pao%diag_distribution, need to redistribute using distribution of matrix_s
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_create(matrix_U, matrix_type="N", template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U)
      CALL dbcsr_complete_redistribute(pao%matrix_U, matrix_U)

      ! Multiplying diagonal matrices is a local operation.
      ! To take advantage of this we're using an iterator instead of calling dbcsr_multiply().
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,matrix_U) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_U,block_Y,block_A,block_B,block_N,block_N_inv,found)
      CALL dbcsr_iterator_start(iter, matrix_U)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_U)
         iatom = arow; CPASSERT(arow == acol)

         CALL dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=block_Y, found=found)
         CPASSERT(ASSOCIATED(block_Y))

         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_A, row=iatom, col=iatom, block=block_A, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=iatom, col=iatom, block=block_N_inv, found=found)
         CPASSERT(ASSOCIATED(block_A) .AND. ASSOCIATED(block_N_inv))

         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_B, row=iatom, col=iatom, block=block_B, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N, row=iatom, col=iatom, block=block_N, found=found)
         CPASSERT(ASSOCIATED(block_B) .AND. ASSOCIATED(block_N))

         block_A = MATMUL(MATMUL(block_N_inv, block_U), block_Y)
         block_B = MATMUL(MATMUL(block_N, block_U), block_Y)
      ENDDO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_release(matrix_U)

      CALL timestop(handle)
   END SUBROUTINE pao_update_AB

! **************************************************************************************************
!> \brief Initialize PAO parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_init_exp(pao, qs_env)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_init_linpot(pao, qs_env)
      CASE (pao_gth_param)
         CALL pao_param_init_gth(pao, qs_env)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_init

! **************************************************************************************************
!> \brief Finalize PAO parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_finalize_exp(pao)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_finalize_linpot(pao)
      CASE (pao_gth_param)
         CALL pao_param_finalize_gth(pao)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_finalize

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param pao ...
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count(pao, qs_env, ikind, nparams)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_count', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_count_exp(qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_count_linpot(pao, qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_gth_param)
         CALL pao_param_count_gth(qs_env, ikind=ikind, nparams=nparams)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_count

! **************************************************************************************************
!> \brief Fills matrix_X with an initial guess
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_initial_guess(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initial_guess', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(dp)                                           :: norm
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_initguess_exp(pao)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_initguess_linpot(pao, qs_env)
      CASE (pao_gth_param)
         CALL pao_param_initguess_gth(pao)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      norm = dbcsr_frobenius_norm(pao%matrix_X)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Made initial guess for matrix_X with norm:", norm

      IF (pao%precondition) THEN
         !TODO: multiplying a matrix into itself while retaining sparsity seems to be broken
         CALL dbcsr_copy(matrix_tmp, pao%matrix_X)
         CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_precon, matrix_tmp, &
                             0.0_dp, pao%matrix_X, retain_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pao_param_initial_guess

! **************************************************************************************************
!> \brief Calculate new matrix U
!> \param pao ...
!> \param qs_env ...
!> \param matrix_M ...
!> \param matrix_G ...
!> \param penalty ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U(pao, qs_env, matrix_M, matrix_G, penalty, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), OPTIONAL                         :: matrix_M, matrix_G
      REAL(dp), INTENT(OUT), OPTIONAL                    :: penalty
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (PRESENT(penalty)) penalty = 0.0_dp

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_calc_U_exp(pao, matrix_M, matrix_G)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_calc_U_linpot(pao, qs_env, penalty, matrix_M, matrix_G, forces)
      CASE (pao_gth_param)
         CALL pao_calc_U_gth(pao, penalty, matrix_M, matrix_G)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL pao_assert_unitary(pao)
      CALL timestop(handle)
   END SUBROUTINE pao_calc_U

! **************************************************************************************************
!> \brief Debugging routine, check unitaryness of U
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_assert_unitary(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_assert_unitary', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: acol, arow, group, handle, i, iatom, M, N
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      REAL(dp)                                           :: delta_max
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_test, tmp1, tmp2
      TYPE(dbcsr_iterator_type)                          :: iter

      IF (pao%check_unitary_tol < 0.0_dp) RETURN ! no checking

      CALL timeset(routineN, handle)
      delta_max = 0.0_dp

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,blk_sizes_pri,blk_sizes_pao,delta_max) &
!$OMP PRIVATE(iter,arow,acol,iatom,N,M,block_test,tmp1,tmp2)
      CALL dbcsr_iterator_start(iter, pao%matrix_U)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_test)
         iatom = arow; CPASSERT(arow == acol)
         N = blk_sizes_pri(iatom) ! size of primary basis
         M = blk_sizes_pao(iatom) ! size of pao basis

         ! we only need the upper left "PAO-corner" to be unitary
         ALLOCATE (tmp1(N, M), tmp2(M, M))
         tmp1 = block_test(:, 1:M)
         tmp2 = MATMUL(TRANSPOSE(tmp1), tmp1)
         DO i = 1, M
            tmp2(i, i) = tmp2(i, i)-1.0_dp
         ENDDO

!$OMP ATOMIC
         delta_max = MAX(delta_max, MAXVAL(ABS(tmp2)))

         DEALLOCATE (tmp1, tmp2)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      group = dbcsr_mp_group(dbcsr_distribution_mp(dbcsr_distribution(pao%matrix_U)))
      CALL mp_max(delta_max, group)
      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO| checked unitaryness, max delta:', delta_max
      IF (delta_max > pao%check_unitary_tol) &
         CPABORT("Found bad unitaryness:"//cp_to_string(delta_max))

      CALL timestop(handle)
   END SUBROUTINE pao_assert_unitary

END MODULE pao_param
