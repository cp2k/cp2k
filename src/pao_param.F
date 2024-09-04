!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Front-End for any PAO parametrization
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_multiply,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE pao_input,                       ONLY: pao_equi_param,&
                                              pao_exp_param,&
                                              pao_fock_param,&
                                              pao_gth_param,&
                                              pao_rotinv_param
   USE pao_param_equi,                  ONLY: pao_calc_AB_equi,&
                                              pao_param_count_equi,&
                                              pao_param_finalize_equi,&
                                              pao_param_init_equi,&
                                              pao_param_initguess_equi
   USE pao_param_exp,                   ONLY: pao_calc_AB_exp,&
                                              pao_param_count_exp,&
                                              pao_param_finalize_exp,&
                                              pao_param_init_exp,&
                                              pao_param_initguess_exp
   USE pao_param_gth,                   ONLY: pao_calc_AB_gth,&
                                              pao_param_count_gth,&
                                              pao_param_finalize_gth,&
                                              pao_param_init_gth,&
                                              pao_param_initguess_gth
   USE pao_param_linpot,                ONLY: pao_calc_AB_linpot,&
                                              pao_param_count_linpot,&
                                              pao_param_finalize_linpot,&
                                              pao_param_init_linpot,&
                                              pao_param_initguess_linpot
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param'

   PUBLIC :: pao_calc_AB, pao_param_count, pao_param_initial_guess
   PUBLIC :: pao_param_init, pao_param_finalize

CONTAINS

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
!> \param penalty ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB(pao, qs_env, ls_scf_env, gradient, penalty, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient
      REAL(dp), INTENT(OUT), OPTIONAL                    :: penalty
      REAL(dp), DIMENSION(:, :), INTENT(OUT), OPTIONAL   :: forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_AB'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (PRESENT(penalty)) penalty = 0.0_dp
      IF (PRESENT(forces)) forces(:, :) = 0.0_dp

      !calculate matrix_A/B = Function of matrix_X
      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_calc_AB_exp(pao, qs_env, ls_scf_env, gradient)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_calc_AB_linpot(pao, qs_env, ls_scf_env, gradient, penalty, forces)
      CASE (pao_gth_param)
         CALL pao_calc_AB_gth(pao, qs_env, ls_scf_env, gradient, penalty)
      CASE (pao_equi_param)
         CALL pao_calc_AB_equi(pao, qs_env, ls_scf_env, gradient, penalty)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE pao_calc_AB

! **************************************************************************************************
!> \brief Initialize PAO parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_param_init'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_init_exp(pao, qs_env)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_init_linpot(pao, qs_env)
      CASE (pao_gth_param)
         CALL pao_param_init_gth(pao, qs_env)
      CASE (pao_equi_param)
         CALL pao_param_init_equi(pao)
      CASE DEFAULT
         CPABORT("PAO: unknown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_init

! **************************************************************************************************
!> \brief Finalize PAO parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_finalize_exp(pao)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_finalize_linpot(pao)
      CASE (pao_gth_param)
         CALL pao_param_finalize_gth(pao)
      CASE (pao_equi_param)
         CALL pao_param_finalize_equi()
      CASE DEFAULT
         CPABORT("PAO: unknown parametrization")
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

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_param_count'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_count_exp(qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_count_linpot(pao, qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_gth_param)
         CALL pao_param_count_gth(qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_equi_param)
         CALL pao_param_count_equi(qs_env, ikind=ikind, nparams=nparams)
      CASE DEFAULT
         CPABORT("PAO: unknown parametrization")
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

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initial_guess'

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
      CASE (pao_equi_param)
         CALL pao_param_initguess_equi(pao, qs_env)
      CASE DEFAULT
         CPABORT("PAO: unknown parametrization")
      END SELECT

      norm = dbcsr_frobenius_norm(pao%matrix_X)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Made initial guess for matrix_X with norm:", norm

      IF (pao%precondition) THEN
         !TODO: multiplying a matrix into itself while retaining sparsity seems to be broken
         CALL dbcsr_copy(matrix_tmp, pao%matrix_X)
         CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_precon, matrix_tmp, &
                             0.0_dp, pao%matrix_X, retain_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE pao_param_initial_guess

END MODULE pao_param
