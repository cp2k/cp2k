!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief optimizes the parallel distribution trying to minimize a cost function
!> \note
!>      currently optimizes a 2D distribution for a 2D grid of CPUs
!>      see distribute_2d_monte_carlo
!>
!>      The approach currently is based on Monte Carlo like annealing to find the
!>      a good approximation of the global minimum [which appears to work well]
!> \par History
!>      04.2007 Created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE distribution_optimize
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: opt_target_median,&
                                             opt_target_minimum
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             reset_to_next_rng_substream,&
                                             rng_stream_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'distribution_optimize'

  PUBLIC :: distribute_2d_monte_carlo

CONTAINS

! *****************************************************************************
!> \brief minimizes the cost per CPU for a 2D distribution
!>      i.e. the maximum value of the sum of the elements of cost_ele_pair
!>      belong to a given CPU is minimized.
!>      For example, if cost_ele_pair=1 each cpu will have approximately the
!>      same number of elements, and the maximum number of blocks on a given CPU is minimum.
!>      it is up to some cost model to define reasonable costs of each element pair
!>
!>      P x Q : processor grid
!>      nele : number of elements
!>      local_rows[1..nele] : to which CPU row [1..P] belongs this ele
!>      local_cols[1..nele] : to which CPU col [1..Q] belongs this ele
!> \par History
!>      04.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE distribute_2d_monte_carlo(cost_ele_pair,P,Q,nele,local_rows,local_cols,mc_section,para_env,error)
    INTEGER, INTENT(IN)                      :: P, Q, nele, &
                                                cost_ele_pair(nele,nele)
    INTEGER, INTENT(OUT)                     :: local_rows(nele), &
                                                local_cols(nele)
    TYPE(section_vals_type), POINTER         :: mc_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_2d_monte_carlo', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=20)                        :: opt_target_used
    INTEGER :: best_cpu, best_iac, best_max_cost, cost_optimal, handle, i, &
      iac, icpu, iele, iswap, itmp, itries, iw, jswap, max_cost, ncycles, &
      nmoves, nreduction, opt_target, pos_max_cost, random_rowcol
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: max_cost_para
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: cost_cpu
    INTEGER :: local_rows_pos(nele), local_rows_best(nele), &
      local_cols_pos(nele), local_cols_best(nele)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: random_move, &
                                                reduction_factor, swap_prob, &
                                                T, T_in, t_initial, tol
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,mc_section,"PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)
    !
    ! get input options
    !
    CALL section_vals_val_get(mc_section,"MOVES",i_val=nmoves,error=error)
    CALL section_vals_val_get(mc_section,"CYCLES",i_val=ncycles,error=error)
    CALL section_vals_val_get(mc_section,"REDUCTION_STEPS",i_val=nreduction,error=error)
    CALL section_vals_val_get(mc_section,"REDUCTION_FACTOR",r_val=reduction_factor,error=error)
    CALL section_vals_val_get(mc_section,"TEMPERATURE",r_val=T_in,error=error)
    CALL section_vals_val_get(mc_section,"TOLERANCE",r_val=tol,error=error)
    CALL section_vals_val_get(mc_section,"SWAP_PROBABILITY",r_val=swap_prob,error=error)
    CALL section_vals_val_get(mc_section,"OPT_TARGET",i_val=opt_target,error=error)

    SELECT CASE(opt_target)
    CASE(opt_target_minimum)
       opt_target_used = 'SMALLEST MAXIMUM'
    CASE(opt_target_median)
       opt_target_used = 'PERFECT BALANCE'
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    ! the best possible cost is a perfectly equal load for each cpu
    ! POTENTIAL OVERFLOW HERE !
    cost_optimal=(SUM(INT(cost_ele_pair,int_8))+INT(P*Q-1,int_8))/INT(P*Q,int_8)

    IF (iw>0) THEN
       WRITE(iw,'(/,T2,A)') "2D_MC| Monte Carlo annealing to optimize the distribution_2d"
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Number of moves",nmoves
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Number of annealing cycles",ncycles
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Number of reduction steps per cycle",nreduction
       WRITE(iw,'(T2,A,T61,F20.10)') "2D_MC| Reduction factor per step",reduction_factor
       WRITE(iw,'(T2,A,T61,F20.10)') "2D_MC| Termination tolerance",tol
       WRITE(iw,'(T2,A,T61,F20.10)') "2D_MC| Maximum temperature",T_in
       WRITE(iw,'(T2,A,T61,F20.10)') "2D_MC| Swap probability",swap_prob
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Number of processor rows ",P
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Number of processor cols ",Q
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Number of elements ",nele
       WRITE(iw,'(T2,A,T61,A20)') "2D_MC| What do we minimize ",ADJUSTR(TRIM(opt_target_used))
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Cost of optimal distribution ",cost_optimal
    ENDIF

    ! adjust this to be in units of the maximum block cost
    T_in=T_in*MAXVAL(cost_ele_pair)
    ! allocate a cost per CPU structure
    ALLOCATE(cost_cpu(P,Q))

    ! an initial guess distribution, quite trivial, the result should not depend on this
    ! each CPU gets a slightly different initial distribution
    DO i=1,nele
        local_rows(i)=MOD(i-1+para_env%mepos,P)+1
        local_cols(i)=MOD(i-1,Q)+1
    ENDDO
    CALL compute_cost_cpu(cost_cpu,cost_ele_pair,P,Q,nele,local_rows,local_cols)
    SELECT CASE(opt_target)
    CASE(opt_target_minimum)
      max_cost=MAXVAL(cost_cpu)
    CASE(opt_target_median)
      max_cost=deviance_cost_cpu(cost_cpu,P,Q)
    END SELECT


    ! initialize the random number generators
    ! each CPU should have a different stream of random numbers
    ! so the optimization problem is parallellized (Ncpu independent searches)
    NULLIFY(rng_stream)
    CALL create_rng_stream(rng_stream,"distribute_2d_monte_carlo",distribution_type=UNIFORM,error=error)
    DO i=1,para_env%mepos
       CALL reset_to_next_rng_substream(rng_stream,error=error)
    ENDDO

    ! initialize the Monte Carlo
    best_max_cost = max_cost
    local_rows_best = local_rows
    local_cols_best = local_cols
    pos_max_cost = max_cost
    local_rows_pos = local_rows
    local_cols_pos = local_cols
    best_iac=0

    iac=0
    best_iac=0

    CALL CPU_TIME(t_initial)
    DO itries=1,nmoves
       ! do annealing
       IF (MOD(itries,nmoves/ncycles)==1) THEN
          T=T_in
       ENDIF
       IF (MOD(itries,nmoves/(ncycles*nreduction))==1) THEN
          T=T*reduction_factor
       ENDIF

       ! generate trial configuration
       local_rows = local_rows_pos
       local_cols = local_cols_pos

       random_move   = next_random_number(rng_stream,error=error)
       random_rowcol = next_random_number(rng_stream,low=0,high=1,error=error)
       IF (random_move < 1.0_dp-swap_prob) THEN ! assign a col/row to a random CPU with 1.0-swap_prob probability
          IF (random_rowcol==0) THEN
             iele=next_random_number(rng_stream,low=1,high=nele,error=error)
             icpu=next_random_number(rng_stream,low=1,high=P,error=error)
             local_rows(iele)=icpu
          ELSE
             iele=next_random_number(rng_stream,low=1,high=nele,error=error)
             icpu=next_random_number(rng_stream,low=1,high=Q,error=error)
             local_cols(iele)=icpu
          ENDIF
       ELSE ! randomly swap a col/row with swap_prob probability
          IF (random_rowcol==0) THEN
             iswap=next_random_number(rng_stream,low=1,high=nele,error=error)
             jswap=next_random_number(rng_stream,low=1,high=nele,error=error)
             itmp=local_rows(iswap)
             local_rows(iswap)=local_rows(jswap)
             local_rows(jswap)=itmp
          ELSE
             iswap=next_random_number(rng_stream,low=1,high=nele,error=error)
             jswap=next_random_number(rng_stream,low=1,high=nele,error=error)
             itmp=local_cols(iswap)
             local_cols(iswap)=local_cols(jswap)
             local_cols(jswap)=itmp
          ENDIF
       ENDIF

       ! evaluate its cost
       iac = iac + 1
       CALL update_cost_cpu_diff(cost_cpu,cost_ele_pair,P,Q,nele,local_rows,local_cols,local_rows_pos,local_cols_pos)

       SELECT CASE(opt_target)
       CASE(opt_target_minimum)
         max_cost=MAXVAL(cost_cpu)
       CASE(opt_target_median)
         max_cost=deviance_cost_cpu(cost_cpu,P,Q)
       END SELECT

        ! accept / reject
       random_move   = next_random_number(rng_stream,error=error)
       IF (max_cost < pos_max_cost .OR. (EXP( -(max_cost - pos_max_cost)/T) .GT. random_move )) THEN
          pos_max_cost   = max_cost
          local_rows_pos = local_rows
          local_cols_pos = local_cols
       ELSE
          ! restore the old position if rejected
          CALL update_cost_cpu_diff(cost_cpu,cost_ele_pair,P,Q,nele,local_rows_pos,local_cols_pos,local_rows,local_cols)
       ENDIF

       ! if this one happens to be the best so far, store this
       IF (max_cost < best_max_cost) THEN
          best_max_cost = max_cost
          local_rows_best = local_rows
          local_cols_best = local_cols
          best_iac = iac
          !
          ! if we are within tol of the optimal cost, stop optimizing, or cost model is certainly not that good
          !
          IF (REAL(max_cost-cost_optimal,KIND=dp)/cost_optimal.LT.tol) EXIT
          ! write(iw,*) best_iac,best_max_cost, T
       ENDIF

    ENDDO

    ! now get the best results from all different CPUs
    ALLOCATE(max_cost_para(0:para_env%num_pe-1))
    max_cost_para=0
    max_cost_para(para_env%mepos)=best_max_cost
    CALL mp_sum(max_cost_para,para_env%group)
    best_max_cost=max_cost_para(0)
    best_cpu=0
    DO i=1,para_env%num_pe-1
       IF (max_cost_para(i)<best_max_cost) THEN
          best_max_cost=max_cost_para(i)
          best_cpu=i
       ENDIF
    ENDDO
    CALL mp_bcast(local_rows_best,best_cpu,para_env%group)
    CALL mp_bcast(local_cols_best,best_cpu,para_env%group)

    CALL compute_cost_cpu(cost_cpu,cost_ele_pair,P,Q,nele,local_rows_best,local_cols_best)
    IF (iw>0) THEN
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Cost of found distribution ",best_max_cost
       WRITE(iw,'(T2,A,T61,I20)') "2D_MC| Difference in percent",CEILING(100*REAL((best_max_cost-cost_optimal))/cost_optimal)
    ENDIF

    CALL cp_print_key_finished_output(iw,logger,mc_section,&
         "PROGRAM_RUN_INFO", error=error)

    CALL delete_rng_stream(rng_stream,error=error)

    local_rows=local_rows_best
    local_cols=local_cols_best

    CALL timestop(handle)

  END SUBROUTINE distribute_2d_monte_carlo

! *****************************************************************************
!> \brief compute the cost for each CPU given a local_rows/local_cols distribution
!> \note
!>      quadratic in the number of elements and thus expensive.
!>      see update update_cost_cpu_row / update_cost_cpu_col / update_cost_cpu_diff
!> \par History
!>      04.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE compute_cost_cpu(cost_cpu,cost_ele_pair,P,Q,nele,local_rows,local_cols)
    INTEGER, INTENT(IN)                      :: P, Q
    INTEGER, INTENT(out)                     :: cost_cpu(P,Q)
    INTEGER, INTENT(IN)                      :: nele, &
                                                cost_ele_pair(nele,nele), &
                                                local_rows(nele), &
                                                local_cols(nele)

    INTEGER                                  :: i, j

    cost_cpu = 0
    DO j=1,nele
      DO i=1,nele
         cost_cpu(local_rows(i),local_cols(j))=cost_cpu(local_rows(i),local_cols(j))+cost_ele_pair(i,j)
      ENDDO
    ENDDO

  END SUBROUTINE

! *****************************************************************************
!> \brief update the cost_cpu matrix after changing the assignment of a single ele
!> \note
!>      linear scaling in the number of eles
!> \par History
!>      04.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE update_cost_cpu_row(cost_cpu,cost_ele_pair,P,Q,nele,local_rows,local_cols,iele,icpu_old,icpu_new)
    INTEGER, INTENT(IN)                      :: P, Q
    INTEGER, INTENT(INOUT)                   :: cost_cpu(P,Q)
    INTEGER, INTENT(IN) :: nele, cost_ele_pair(nele,nele), local_rows(nele), &
      local_cols(nele), iele, icpu_old, icpu_new

    INTEGER                                  :: i

    DO i=1,nele
       cost_cpu(icpu_old,local_cols(i))=cost_cpu(icpu_old,local_cols(i))-cost_ele_pair(iele,i)
       cost_cpu(icpu_new,local_cols(i))=cost_cpu(icpu_new,local_cols(i))+cost_ele_pair(iele,i)
    ENDDO

  END SUBROUTINE update_cost_cpu_row

! *****************************************************************************
!> \brief update the cost_cpu matrix after changing the assignment of a single ele
!> \note
!>      linear scaling in the number of eles
!> \par History
!>      04.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE update_cost_cpu_col(cost_cpu,cost_ele_pair,P,Q,nele,local_rows,local_cols,iele,icpu_old,icpu_new)
    INTEGER, INTENT(IN)                      :: P, Q
    INTEGER, INTENT(INOUT)                   :: cost_cpu(P,Q)
    INTEGER, INTENT(IN)                      :: nele, cost_ele_pair(nele,nele)
    INTEGER                                  :: local_rows(nele), &
                                                local_cols(nele), iele, &
                                                icpu_old, icpu_new

    INTEGER                                  :: i

    DO i=1,nele
       cost_cpu(local_rows(i),icpu_old)=cost_cpu(local_rows(i),icpu_old)-cost_ele_pair(i,iele)
       cost_cpu(local_rows(i),icpu_new)=cost_cpu(local_rows(i),icpu_new)+cost_ele_pair(i,iele)
    ENDDO

  END SUBROUTINE update_cost_cpu_col

! *****************************************************************************
!> \brief update the cost_cpu matrix after changing the assignment of the local_row or local_col
!> \note
!>      only checked if either the row or the col distribution changes
!>      if both change at the same time one has to check this ...
!>      (i.e. is the cost of the elements that change twice taken correctly into account)
!> \par History
!>      04.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE update_cost_cpu_diff(cost_cpu,cost_ele_pair,P,Q,nele,local_rows_new,local_cols_new,local_rows_old,local_cols_old)
    INTEGER, INTENT(IN)                      :: P, Q
    INTEGER, INTENT(INOUT)                   :: cost_cpu(P,Q)
    INTEGER, INTENT(IN)                      :: nele, cost_ele_pair(nele,nele)
    INTEGER :: local_rows_new(nele), local_cols_new(nele), &
      local_rows_old(nele), local_cols_old(nele)

    INTEGER                                  :: i

    DO i=1,nele
       IF (local_rows_new(i).NE.local_rows_old(i)) THEN
          CALL update_cost_cpu_row(cost_cpu,cost_ele_pair,P,Q,nele, &
                                   local_rows_new,local_cols_new,i,local_rows_old(i),local_rows_new(i))
       ENDIF
       IF (local_cols_new(i).NE.local_cols_old(i)) THEN
          CALL update_cost_cpu_col(cost_cpu,cost_ele_pair,P,Q,nele, &
                                   local_rows_new,local_cols_new,i,local_cols_old(i),local_cols_new(i))
       ENDIF
    ENDDO
  END SUBROUTINE update_cost_cpu_diff

! *****************************************************************************
!> \brief Calculate the deviance (Sqrt(Var)) of the cost_cpu array
!> \par History
!>      12.2009 P. Manninen
! *****************************************************************************
  FUNCTION deviance_cost_cpu(cost_cpu,P,Q) RESULT(idev)
    INTEGER, INTENT(IN)                      :: P, Q, cost_cpu(P,Q)
    INTEGER                                  :: idev

    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: a, dev, mean, rpq

    rpq=P*Q
    mean=0.0_dp
    DO j=1,Q
      DO i=1,P
         mean=mean+cost_cpu(i,j)
      ENDDO
    ENDDO
    mean=mean/rpq
    dev=0.0_dp
    DO j=1,Q
      DO i=1,P
        a=cost_cpu(i,j)-mean
        dev=dev+a*a
      END DO
    END DO
    dev=SQRT(dev/rpq)
    idev=NINT(dev)
  END FUNCTION deviance_cost_cpu
END MODULE distribution_optimize

