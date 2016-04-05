!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Optimizers used by pao_main.F
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_optimizer
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_arnoldi_extremal, cp_dbcsr_copy, &
        cp_dbcsr_create, cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_release, &
        cp_dbcsr_reserve_diag_blocks, cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set, &
        cp_dbcsr_trace, cp_dbcsr_type
   USE kinds,                           ONLY: dp
   USE pao_input,                       ONLY: pao_opt_bfgs,&
                                              pao_opt_cg
   USE pao_types,                       ONLY: pao_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_opt_init, pao_opt_finalize, pao_opt_new_dir

CONTAINS

! **************************************************************************************************
!> \brief Initialize the optimizer
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_opt_init(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL cp_dbcsr_init(pao%matrix_G_prev)
      CALL cp_dbcsr_copy(pao%matrix_G_prev, pao%matrix_G)
      CALL cp_dbcsr_set(pao%matrix_G_prev, 0.0_dp)

      CALL cp_dbcsr_init(pao%matrix_D)
      CALL cp_dbcsr_copy(pao%matrix_D, pao%matrix_G)
      CALL cp_dbcsr_set(pao%matrix_D, 0.0_dp)

      IF (pao%optimizer == pao_opt_bfgs) &
         CALL pao_opt_init_bfgs(pao)

   END SUBROUTINE pao_opt_init

! **************************************************************************************************
!> \brief Initialize the BFGS optimizer
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_opt_init_bfgs(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      INTEGER, DIMENSION(:), POINTER                     :: nparams

      nparams => cp_dbcsr_row_block_sizes(pao%matrix_X)

      CALL cp_dbcsr_init(pao%matrix_BFGS)
      CALL cp_dbcsr_create(pao%matrix_BFGS, &
                           template=pao%matrix_X, &
                           row_blk_size=nparams, &
                           col_blk_size=nparams, &
                           name="PAO matrix_BFGS")

      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_BFGS)
      CALL cp_dbcsr_set(pao%matrix_BFGS, 0.0_dp)
      CALL cp_dbcsr_add_on_diag(pao%matrix_BFGS, 1.0_dp)

   END SUBROUTINE pao_opt_init_bfgs

! **************************************************************************************************
!> \brief Finalize the optimizer
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_opt_finalize(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL cp_dbcsr_release(pao%matrix_G_prev)
      CALL cp_dbcsr_release(pao%matrix_D)

      IF (pao%optimizer == pao_opt_bfgs) &
         CALL cp_dbcsr_release(pao%matrix_BFGS)

   END SUBROUTINE pao_opt_finalize

! **************************************************************************************************
!> \brief Calculates the new search direction.
!> \param pao ...
!> \param icycle ...
! **************************************************************************************************
   SUBROUTINE pao_opt_new_dir(pao, icycle)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: icycle

      SELECT CASE (pao%optimizer)
      CASE (pao_opt_cg)
         CALL pao_opt_newdir_cg(pao, icycle)
      CASE (pao_opt_bfgs)
         CALL pao_opt_newdir_bfgs(pao, icycle)
      CASE DEFAULT
         CPABORT("PAO: unknown optimizer")
      END SELECT

      ! backup gradient
      CALL cp_dbcsr_copy(pao%matrix_G_prev, pao%matrix_G)

   END SUBROUTINE pao_opt_new_dir

! **************************************************************************************************
!> \brief Conjugate Gradient algorithm
!> \param pao ...
!> \param icycle ...
! **************************************************************************************************
   SUBROUTINE pao_opt_newdir_cg(pao, icycle)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: icycle

      REAL(KIND=dp)                                      :: beta, change, trace_D, trace_D_Gnew, &
                                                            trace_G_mix, trace_G_new, trace_G_prev

      ! determine CG mixing factor
      IF (icycle <= pao%cg_init_steps) THEN
         IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| warming up with steepest descent"
         beta = 0.0_dp
      ELSE
         CALL cp_dbcsr_trace(pao%matrix_G, pao%matrix_G, trace_G_new)
         CALL cp_dbcsr_trace(pao%matrix_G_prev, pao%matrix_G_prev, trace_G_prev)
         CALL cp_dbcsr_trace(pao%matrix_G, pao%matrix_G_prev, trace_G_mix)
         CALL cp_dbcsr_trace(pao%matrix_D, pao%matrix_G, trace_D_Gnew)
         CALL cp_dbcsr_trace(pao%matrix_D, pao%matrix_D, trace_D)
         IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| trace_G_new ", trace_G_new
         IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| trace_G_prev ", trace_G_prev
         IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| trace_G_mix ", trace_G_mix
         IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| trace_D ", trace_D
         IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| trace_D_Gnew", trace_D_Gnew

         IF (trace_G_prev /= 0.0_dp) THEN
            beta = (trace_G_new-trace_G_mix)/trace_G_prev !Polak–Ribiere
         ENDIF

         IF (beta < 0.0_dp) THEN
            IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| resetting because beta < 0"
            beta = 0.0_dp
         ENDIF

         change = trace_D_Gnew**2/trace_D*trace_G_new
         IF (change > pao%cg_reset_limit) THEN
            IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| resetting because change > CG_RESET_LIMIT"
            beta = 0.0_dp
         ENDIF

      ENDIF

      IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|CG| beta: ", beta

      ! calculate new CG direction matrix_D
      CALL cp_dbcsr_add(pao%matrix_D, pao%matrix_G, beta, -1.0_dp)

   END SUBROUTINE pao_opt_newdir_cg

! **************************************************************************************************
!> \brief Broyden–Fletcher–Goldfarb–Shanno algorithm
!> \param pao ...
!> \param icycle ...
! **************************************************************************************************
   SUBROUTINE pao_opt_newdir_bfgs(pao, icycle)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: icycle

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_opt_newdir_bfgs'

      INTEGER                                            :: handle
      LOGICAL                                            :: arnoldi_converged
      REAL(dp)                                           :: eval_max, eval_min, theta, trace_ry, &
                                                            trace_sy, trace_yHy, trace_yy
      TYPE(cp_dbcsr_type)                                :: matrix_Hy, matrix_Hyr, matrix_r, &
                                                            matrix_rr, matrix_ryH, matrix_ryHyr, &
                                                            matrix_s, matrix_y, matrix_yr

      CALL timeset(routineN, handle)

      !TODO add filtering?

      ! Notation according to the book from Nocedal and Wright, see chapter 6.
      IF (icycle > 1) THEN
         ! y = G - G_prev
         CALL cp_dbcsr_init(matrix_y)
         CALL cp_dbcsr_copy(matrix_y, pao%matrix_G)
         CALL cp_dbcsr_add(matrix_y, pao%matrix_G_prev, 1.0_dp, -1.0_dp) ! dG

         ! s = X - X_prev
         CALL cp_dbcsr_init(matrix_s)
         CALL cp_dbcsr_copy(matrix_s, pao%matrix_D)
         CALL cp_dbcsr_scale(matrix_s, pao%linesearch%step_size) ! dX

         ! sy = MATMUL(TRANPOSE(s), y)
         CALL cp_dbcsr_trace(matrix_s, matrix_y, trace_sy)

         ! heuristic initialization
         IF (icycle == 2) THEN
            CALL cp_dbcsr_trace(matrix_Y, matrix_Y, trace_yy)
            CALL cp_dbcsr_scale(pao%matrix_BFGS, trace_sy/trace_yy)
            IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|BFGS| Initializing with:", trace_sy/trace_yy
         ENDIF

         ! Hy = MATMUL(H, y)
         CALL cp_dbcsr_init(matrix_Hy)
         CALL cp_dbcsr_create(matrix_Hy, template=pao%matrix_X, matrix_type="N")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_BFGS, matrix_y, 0.0_dp, matrix_Hy)

         ! yHy = MATMUL(TRANPOSE(y), Hy)
         CALL cp_dbcsr_trace(matrix_y, matrix_Hy, trace_yHy)

         ! Use damped BFGS algorithm to ensure H remains positive definite.
         ! See chapter 18 in Nocedal and Wright's book for details.
         ! The formulas were adopted to inverse Hessian algorithm.
         IF (trace_sy < 0.2_dp*trace_yHy) THEN
            theta = 0.8_dp*trace_yHy/(trace_yHy-trace_sy)
            IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|BFGS| Dampening theta:", theta
         ELSE
            theta = 1.0
         ENDIF

         ! r = theta*s + (1-theta)*Hy
         CALL cp_dbcsr_init(matrix_r)
         CALL cp_dbcsr_copy(matrix_r, matrix_s)
         CALL cp_dbcsr_add(matrix_r, matrix_Hy, theta, (1.0_dp-theta))

         ! use t instead of y to update B matrix
         CALL cp_dbcsr_trace(matrix_r, matrix_y, trace_ry)
         CPASSERT(trace_RY > 0.0_dp)

         ! yr = MATMUL(y, TRANSPOSE(r))
         CALL cp_dbcsr_init(matrix_yr)
         CALL cp_dbcsr_create(matrix_yr, template=pao%matrix_BFGS, matrix_type="N")
         CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_y, matrix_r, 0.0_dp, matrix_yr)

         ! Hyr = MATMUL(H, yr)
         CALL cp_dbcsr_init(matrix_Hyr)
         CALL cp_dbcsr_create(matrix_Hyr, template=pao%matrix_BFGS, matrix_type="N")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_BFGS, matrix_yr, 0.0_dp, matrix_Hyr)

         ! ryH = MATMUL(TRANSPOSE(yr), H)
         CALL cp_dbcsr_init(matrix_ryH)
         CALL cp_dbcsr_create(matrix_ryH, template=pao%matrix_BFGS, matrix_type="N")
         CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_yr, pao%matrix_BFGS, 0.0_dp, matrix_ryH)

         ! ryHry = MATMUL(ryH,yr)
         CALL cp_dbcsr_init(matrix_ryHyr)
         CALL cp_dbcsr_create(matrix_ryHyr, template=pao%matrix_BFGS, matrix_type="N")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ryH, matrix_yr, 0.0_dp, matrix_ryHyr)

         ! rr = MATMUL(r,TRANSPOSE(r))
         CALL cp_dbcsr_init(matrix_rr)
         CALL cp_dbcsr_create(matrix_rr, template=pao%matrix_BFGS, matrix_type="N")
         CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_r, matrix_r, 0.0_dp, matrix_rr)

         ! H = H - Hyr/ry - ryH/ry + ryHyr/(ry**2) + rr/ry
         CALL cp_dbcsr_add(pao%matrix_BFGS, matrix_HYR, 1.0_dp, -1.0_dp/trace_ry)
         CALL cp_dbcsr_add(pao%matrix_BFGS, matrix_ryH, 1.0_dp, -1.0_dp/trace_ry)
         CALL cp_dbcsr_add(pao%matrix_BFGS, matrix_ryHyr, 1.0_dp, +1.0_dp/(trace_ry**2))
         CALL cp_dbcsr_add(pao%matrix_BFGS, matrix_rr, 1.0_dp, +1.0_dp/trace_ry)

         ! clean up
         CALL cp_dbcsr_release(matrix_y)
         CALL cp_dbcsr_release(matrix_s)
         CALL cp_dbcsr_release(matrix_r)
         CALL cp_dbcsr_release(matrix_Hy)
         CALL cp_dbcsr_release(matrix_yr)
         CALL cp_dbcsr_release(matrix_Hyr)
         CALL cp_dbcsr_release(matrix_ryH)
         CALL cp_dbcsr_release(matrix_ryHyr)
         CALL cp_dbcsr_release(matrix_rr)
      ENDIF

      ! approximate condition of Hessian
      !TODO: good setting for arnoldi?
      CALL cp_dbcsr_arnoldi_extremal(pao%matrix_BFGS, eval_max, eval_min, max_iter=100, &
                                     threshold=1e-2_dp, converged=arnoldi_converged)
      IF (pao%iw_opt > 0) WRITE (pao%iw_opt, *) "PAO|BFGS| arnoldi(H) ", arnoldi_converged, eval_min, eval_max

      ! calculate new direction
      ! d = MATMUL(H, -g)
      CALL cp_dbcsr_multiply("N", "N", -1.0_dp, pao%matrix_BFGS, pao%matrix_G, 0.0_dp, pao%matrix_D)

      CALL timestop(handle)
   END SUBROUTINE pao_opt_newdir_bfgs

END MODULE pao_optimizer
