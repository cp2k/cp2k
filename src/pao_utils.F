!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Helper routines shared by multiple PAO modules
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_utils
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_frobenius_norm,&
                                             cp_dbcsr_iterator,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop,&
                                             cp_dbcsr_trace,&
                                             cp_dbcsr_type
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE pao_types,                       ONLY: pao_env_type
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_utils'

  PUBLIC :: fill_randomly, pao_calc_CG_beta, pao_test_convegence
  PUBLIC :: pao_print_oneliner, pao_print_step_delimiter

CONTAINS

! *****************************************************************************
!> \brief Print single line summary
!> \param pao ...
!> \param ls_scf_env ...
!> \param energy ...
!> \param mode ...
! *****************************************************************************
  SUBROUTINE pao_print_oneliner(pao, ls_scf_env, energy, mode)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    REAL(KIND=dp), INTENT(IN)                :: energy
    CHARACTER(LEN=*), INTENT(IN)             :: mode

    REAL(KIND=dp)                            :: delta_t, energy_diff, &
                                                loop_eps, norm_G, t2

    norm_G = cp_dbcsr_frobenius_norm(pao%matrix_G)

    t2 = m_walltime()
    delta_t = t2 - pao%step_start_time
    pao%step_start_time = t2


    energy_diff = energy - pao%energy_prev
    loop_eps = ABS(energy_diff) / ls_scf_env%nelectron_total
    IF(pao%iw>0) WRITE (pao%iw,'(A,I6,1X,A8,1X,F20.15,1X,E10.3,1X,F8.3,1X,E10.3,1X,F6.3)')&
                 " PAO| step ",&
                 pao%istep,&
                 TRIM(mode),&
                 energy,&
                 norm_G,&
                 pao%linesearch%step_size,& !prev step, which let to the current energy
                 loop_eps,&
                 delta_t

  END SUBROUTINE pao_print_oneliner


! *****************************************************************************
!> \brief Write iteration delimiter
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_print_step_delimiter(pao)
    TYPE(pao_env_type), POINTER              :: pao

    IF(pao%iw>0) WRITE(pao%iw,"(A,I9,A)") " PAO| ======================= Iteration: ",&
                 pao%istep," ============================="

  END SUBROUTINE pao_print_step_delimiter

! *****************************************************************************
!> \brief Test whether the PAO optimization has reached convergence
!> \param pao ...
!> \param ls_scf_env ...
!> \param threshold ...
!> \param new_energy ...
!> \param is_converged ...
! *****************************************************************************
  SUBROUTINE pao_test_convegence(pao, ls_scf_env, threshold, new_energy, is_converged)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    REAL(KIND=dp), INTENT(IN)                :: threshold, new_energy
    LOGICAL, INTENT(OUT)                     :: is_converged

    REAL(KIND=dp)                            :: energy_diff, loop_eps, norm_G

     is_converged = .FALSE.

     ! print energy improvement
     energy_diff = new_energy - pao%energy_prev
     pao%energy_prev = new_energy
     norm_G = cp_dbcsr_frobenius_norm(pao%matrix_G) ! calc norm before CG mixing
     IF(pao%iw>0) WRITE(pao%iw,*) "PAO| energy:", new_energy
     IF(pao%iw>0) WRITE(pao%iw,*) "PAO| gradient norm:", norm_G

     IF(pao%istep>1) THEN
        IF(pao%iw>0) WRITE(pao%iw,*) "PAO| energy improvement:", energy_diff
        IF(energy_diff>0.0_dp) CPWARN("PAO| energy increased")

        !blk_sizes => cp_dbcsr_col_block_sizes(pao%matrix_G)
        !loop_eps = norm_G/SUM(blk_sizes**2)
        loop_eps = ABS(energy_diff) / ls_scf_env%nelectron_total
        IF(pao%iw>0) WRITE(pao%iw,*) "PAO| loop_eps:", loop_eps

        IF(loop_eps < threshold) THEN
          IF(pao%iw>0) WRITE(pao%iw,*) "PAO| converged after ", pao%istep, " steps :-)"
          is_converged = .TRUE.
        ENDIF
     ENDIF
  END SUBROUTINE pao_test_convegence


! *****************************************************************************
!> \brief Conjugate gradient algorithm
!> \param pao ...
!> \param beta ...
! *****************************************************************************
  SUBROUTINE pao_calc_CG_beta(pao, beta)
    TYPE(pao_env_type), POINTER              :: pao
    REAL(KIND=dp), INTENT(OUT)               :: beta

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_CG_beta', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: trace_G_mix, trace_G_new, &
                                                trace_G_prev

    beta = 0.0_dp !default

    IF(pao%istep > pao%cg_init_steps) THEN
       CALL cp_dbcsr_trace(pao%matrix_G, pao%matrix_G, trace_G_new)
       CALL cp_dbcsr_trace(pao%matrix_G_prev, pao%matrix_G_prev, trace_G_prev)
       CALL cp_dbcsr_trace(pao%matrix_G, pao%matrix_G_prev, trace_G_mix)
       IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| trace_G_new ", trace_G_new
       IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| trace_G_prev ", trace_G_prev
       IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| trace_G_mix ", trace_G_mix

       IF(trace_G_prev/=0.0_dp) THEN
          beta = (trace_G_new - trace_G_mix) / trace_G_prev  !Polak–Ribiere
       ENDIF

       IF(beta < 0.0_dp) THEN
          IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| resetting because beta < 0"
          beta = 0.0_dp
       ENDIF

       IF(pao%max_outer_pao>0 .AND. pao%istep-pao%last_reset > pao%max_outer_pao) THEN
          IF(pao%iw>0) WRITE(pao%iw,*) "PAO|CG| resetting because of MAX_OUTER_PAO"
          pao%last_reset = pao%istep
          beta = 0.0_dp
       ENDIF
    ENDIF

    IF(pao%iw_cg>0) WRITE(pao%iw_cg,*) "PAO|CG| beta: ", beta

  END SUBROUTINE pao_calc_CG_beta

! *****************************************************************************
!> \brief Fills a DBCSR matrix with random noise.
!> \param matrix ...
! *****************************************************************************
  SUBROUTINE fill_randomly(matrix)
    TYPE(cp_dbcsr_type)                      :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_randomly', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, j
    REAL(dp), DIMENSION(:, :), POINTER       :: block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(rng_stream_type), POINTER           :: rng_stream

    NULLIFY(rng_stream)
    CALL timeset(routineN,handle)

    CALL cp_dbcsr_iterator_start(iter, matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block)
       CALL create_rng_stream(rng_stream, name="pao")
       DO i=1, SIZE(block, 1)
          DO j=1, SIZE(block, 2)
             block(i,j) = next_random_number(rng_stream)
          ENDDO
       ENDDO
       CALL delete_rng_stream(rng_stream)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE fill_randomly

END MODULE pao_utils
