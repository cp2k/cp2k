!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Optimizes the load balance for a given cost matrix for a one dimensional
!>      distribution of processors
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_distribution_mc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hfx_types,                       ONLY: hfx_load_balance_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             reset_to_next_rng_substream,&
                                             rng_stream_type
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC     optimize_distribution

  CONTAINS

! *****************************************************************************
!> \brief Initializes the module variables
!> \param para_env parallel environnment
!> \param nbin_in size of the cost_matrix
!> \param ncpu_in number of cpus
!> \param bin_costs_in cost per (nbin*ncpu)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE init(para_env, nbin_in, ncpu_in, bin_costs_in, load_balance_parameter,my_nbin, my_ncpu, &
                  my_bin_costs, my_cost_cpu, my_cost_optimal, my_cost_total, my_distribution, my_best_distribution, &
                  my_best_cost, my_rng_stream, my_moves, my_cycles, my_reduction_steps, my_reduction_factor,&
                  my_initial_temperature, my_tolerance, my_swap_probability, my_process_id, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER(int_8)                           :: nbin_in, ncpu_in
    INTEGER(int_8), DIMENSION(:), POINTER    :: bin_costs_in
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER                                  :: my_nbin, my_ncpu
    INTEGER(int_8), DIMENSION(:), POINTER    :: my_bin_costs, my_cost_cpu
    INTEGER(int_8)                           :: my_cost_optimal, my_cost_total
    INTEGER, DIMENSION(:), POINTER           :: my_distribution, &
                                                my_best_distribution
    INTEGER(int_8)                           :: my_best_cost
    TYPE(rng_stream_type), POINTER           :: my_rng_stream
    INTEGER                                  :: my_moves, my_cycles, &
                                                my_reduction_steps
    REAL(dp)                                 :: my_reduction_factor, &
                                                my_initial_temperature, &
                                                my_tolerance, &
                                                my_swap_probability
    INTEGER                                  :: my_process_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i

    my_nbin                = nbin_in
    my_ncpu                = ncpu_in
    
    ALLOCATE(my_cost_cpu(my_ncpu))
    ALLOCATE(my_bin_costs(my_nbin))
    ALLOCATE(my_distribution(my_nbin))
    ALLOCATE(my_best_distribution(my_nbin))
    
    my_bin_costs           = bin_costs_in
    my_cost_total          = SUM(my_bin_costs)
    my_cost_optimal        = (my_cost_total + (my_ncpu)-1 )/my_ncpu
    my_moves               = load_balance_parameter%moves
    my_cycles              = load_balance_parameter%cycles
    my_reduction_steps     = load_balance_parameter%reduction_steps
    my_reduction_factor    = load_balance_parameter%reduction_factor
    my_initial_temperature = load_balance_parameter%temperature
    my_tolerance           = load_balance_parameter%tolerance
    my_swap_probability    = load_balance_parameter%swap_probability

    ! initialize the random number generators
    ! each CPU should have a different stream of random numbers
    ! so the optimization problem is parallellized (Ncpu independent searches)
    NULLIFY(my_rng_stream)
    CALL create_rng_stream(my_rng_stream,"distribute_1d_monte_carlo",distribution_type=UNIFORM,error=error)
    DO i=1,my_process_id
       CALL reset_to_next_rng_substream(my_rng_stream,error=error)
    ENDDO

  END SUBROUTINE init

! *****************************************************************************
!> \brief frees memory
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE destructor(my_nbin, my_ncpu, my_bin_costs, my_cost_cpu, my_cost_optimal, &
                        my_cost_total, my_distribution, my_best_distribution, my_best_cost, &
                        my_rng_stream, my_moves, my_cycles, my_reduction_steps, &
                        my_reduction_factor, my_initial_temperature, my_tolerance, &
                        my_swap_probability, error)
    INTEGER                                  :: my_nbin, my_ncpu
    INTEGER(int_8), DIMENSION(:), POINTER    :: my_bin_costs, my_cost_cpu
    INTEGER(int_8)                           :: my_cost_optimal, my_cost_total
    INTEGER, DIMENSION(:), POINTER           :: my_distribution, &
                                                my_best_distribution
    INTEGER(int_8)                           :: my_best_cost
    TYPE(rng_stream_type), POINTER           :: my_rng_stream
    INTEGER                                  :: my_moves, my_cycles, &
                                                my_reduction_steps
    REAL(dp)                                 :: my_reduction_factor, &
                                                my_initial_temperature, &
                                                my_tolerance, &
                                                my_swap_probability
    TYPE(cp_error_type), INTENT(inout)       :: error

    DEALLOCATE(my_cost_cpu)
    DEALLOCATE(my_bin_costs)
    DEALLOCATE(my_distribution)
    DEALLOCATE(my_best_distribution)

    CALL delete_rng_stream(my_rng_stream,error=error)
  END SUBROUTINE destructor

! *****************************************************************************
!> \brief Calculate the cost of a distribution
!> \param distribution maps cost_matrix to a cpu
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE compute_cost_cpu(distribution, my_nbin, my_cost_cpu, my_bin_costs)
    INTEGER                                  :: my_nbin, distribution(my_nbin)
    INTEGER(int_8), DIMENSION(:)             :: my_cost_cpu, my_bin_costs

    INTEGER                                  :: i

    my_cost_cpu = 0
    DO i=1,my_nbin
      my_cost_cpu(distribution(i)) = my_cost_cpu(distribution(i))+my_bin_costs(i)
    ENDDO
  END SUBROUTINE compute_cost_cpu

! *****************************************************************************
!> \brief Calculates the maximal cost of a distribution
!> \param max_cost variable that stores the max_cost
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE compute_max_cost(max_cost, average, my_cost_cpu)
    INTEGER(int_8)                           :: max_cost, average
    INTEGER(int_8), DIMENSION(:)             :: my_cost_cpu

    max_cost=MAXVAL(my_cost_cpu)+SQRT(REAL(SUM((my_cost_cpu-average)**2)/SIZE(my_cost_cpu),kind=dp))

  END SUBROUTINE compute_max_cost

! *****************************************************************************
!> \brief Optimizes a 1d distribution according to a given cost_matrix
!> \param Ntries Number of trial steps
!> \param T_in Initial temperature (for acceptance criterion)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      Uses two basic random updates:
!>         a) swap the cost of two distributions
!>         b) assign the cost of i to j
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE distribute_monte_carlo(my_nbin, my_ncpu, my_bin_costs, my_cost_cpu, my_cost_optimal, &
                                    my_cost_total, my_distribution, my_best_distribution, my_best_cost, &
                                    my_rng_stream, my_moves, my_cycles, my_reduction_steps, &
                                    my_reduction_factor, my_initial_temperature, my_tolerance, &
                                    my_swap_probability, my_process_id, error)
    INTEGER                                  :: my_nbin, my_ncpu
    INTEGER(int_8), DIMENSION(:), POINTER    :: my_bin_costs, my_cost_cpu
    INTEGER(int_8)                           :: my_cost_optimal, my_cost_total
    INTEGER, DIMENSION(:), POINTER           :: my_distribution, &
                                                my_best_distribution
    INTEGER(int_8)                           :: my_best_cost
    TYPE(rng_stream_type), POINTER           :: my_rng_stream
    INTEGER                                  :: my_moves, my_cycles, &
                                                my_reduction_steps
    REAL(dp)                                 :: my_reduction_factor, &
                                                my_initial_temperature, &
                                                my_tolerance, &
                                                my_swap_probability
    INTEGER                                  :: my_process_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, iswap, itmp, itries, j, &
                                                jswap, ncpu, nele, nstep, &
                                                old_distribution(my_nbin)
    INTEGER(int_8)                           :: average, best_max_cost, &
                                                max_cost, max_cost_old
    INTEGER(int_8), DIMENSION(:), POINTER    :: tmp_cost, tmp_cpu_cost
    INTEGER, DIMENSION(:), POINTER           :: tmp_cpu_index, tmp_index
    REAL(dp)                                 :: random_accept, random_move, &
                                                T, T_local

    nele=SIZE(my_bin_costs)
    ncpu=SIZE(my_cost_cpu)
    nstep=MAX(1,ncpu/2)

    ALLOCATE(tmp_cost(nele))
    ALLOCATE(tmp_index(nele))
    ALLOCATE(tmp_cpu_cost(ncpu))
    ALLOCATE(tmp_cpu_index(ncpu))
    tmp_cost=my_bin_costs
    CALL sort(tmp_cost,nele,tmp_index)
    my_cost_cpu=0
!
! assign the largest remaining bin to the CPU with the smallest load
! gives near perfect distributions for a sufficient number of bins ...
! doing this in chunks of nstep (where nstep ~ ncpu) makes this n log n and gives
! each cpu a similar number of tasks.
! it also avoids degenerate cases where thousands of zero sized tasks
! are assigned to the same (least loaded) cpu
!
    DO i=my_nbin,1,-nstep
       tmp_cpu_cost=my_cost_cpu
       CALL sort(tmp_cpu_cost,ncpu,tmp_cpu_index)
       DO j=1,MIN(i,nstep)
          itmp=tmp_cpu_index(j)
          my_distribution(tmp_index(i-j+1))=itmp
          my_cost_cpu(itmp)=my_cost_cpu(itmp)+my_bin_costs(tmp_index(i-j+1))
       ENDDO
    ENDDO

    CALL compute_cost_cpu(my_distribution, my_nbin, my_cost_cpu, my_bin_costs)
    average=SUM(my_cost_cpu)/SIZE(my_cost_cpu)
    CALL compute_max_cost(max_cost, average, my_cost_cpu)
    best_max_cost = max_cost
    my_best_distribution = my_distribution

    max_cost_old = max_cost
    old_distribution = my_distribution

    DEALLOCATE(tmp_cost,tmp_index,tmp_cpu_cost,tmp_cpu_index)

    !
    ! A HELL OF HACK: FIXME FIXME FIXME
    !
    RETURN

    ! we base the temperature on a variance of the results obtained so far
    ! too high a temperature would melt the system instead of converging to a minimum
    T_local=my_initial_temperature*SQRT(REAL(SUM((my_cost_cpu-average)**2)/SIZE(my_cost_cpu),kind=dp))

    DO itries=0,my_moves
      !do annealing 
      IF (MOD(itries,my_moves/my_cycles)==0) THEN
        T=T_local
      ENDIF
      IF ( MOD(itries,my_moves/(my_cycles*my_reduction_steps))==0 ) THEN
        T=T*my_reduction_factor
      ENDIF
      ! do some Monte Carlo
      random_move = next_random_number(my_rng_stream,error=error)
      ! generate trial configuration
      my_distribution = old_distribution
      IF (random_move<1.0_dp-my_swap_probability) THEN
          iswap = next_random_number(my_rng_stream,low=1,high=my_nbin,error=error)
          jswap = next_random_number(my_rng_stream,low=1,high=my_ncpu,error=error)
          IF (iswap>my_nbin .OR. jswap>my_ncpu) CYCLE
          my_distribution(iswap)=jswap
      ELSE
          iswap = next_random_number(my_rng_stream,low=1,high=my_nbin,error=error)
          jswap = next_random_number(my_rng_stream,low=1,high=my_nbin,error=error)
          IF (iswap>my_nbin .OR. jswap>my_nbin) CYCLE
          itmp=my_distribution(iswap)
          my_distribution(iswap)=my_distribution(jswap)
          my_distribution(jswap)=itmp
      ENDIF

      ! evaluate its cost
      CALL compute_cost_cpu(my_distribution, my_nbin, my_cost_cpu, my_bin_costs) 
      CALL compute_max_cost(max_cost, average, my_cost_cpu)

      ! accept / reject
      random_accept = next_random_number(my_rng_stream,error=error)
      IF (max_cost < max_cost_old .OR. (EXP( -(max_cost - max_cost_old)/T) .GT. random_accept)) THEN
        max_cost_old = max_cost
        old_distribution = my_distribution
      ELSE
        ! restore the old position
        CALL compute_cost_cpu(old_distribution, my_nbin, my_cost_cpu, my_bin_costs) 
      ENDIF

      ! if this one happens to be the best so far, store this
      IF (max_cost < best_max_cost) THEN
        best_max_cost = max_cost
        my_best_distribution = my_distribution
      ENDIF
      IF (REAL((max_cost-my_cost_optimal),dp)/my_cost_optimal.LE.my_tolerance) EXIT
    ENDDO

    CALL compute_cost_cpu(my_best_distribution, my_nbin, my_cost_cpu, my_bin_costs)
    my_best_cost = best_max_cost

  END SUBROUTINE distribute_monte_carlo

! *****************************************************************************
!> \brief Interface to the module
!> \param para_env parallel environment
!> \param nbin_in size of the cost matrix
!> \param ncpu_in number of processors
!> \param bin_costs_in cost_matrix
!> \param distribution_vector 1d vector that assigns cost matrix to process
!> \param steps number of Monte Carlo steps to be performed
!> \param initial_temperature initial temperature for acceptance criterion
!> \param optimal_cost stores the optimal cost
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE optimize_distribution(para_env, nbin_in, ncpu_in, bin_costs_in, distribution_vector, &
                                        load_balance_parameter, optimal_cost, my_process_id, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER(int_8)                           :: nbin_in, ncpu_in
    INTEGER(int_8), DIMENSION(:), POINTER    :: bin_costs_in
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    INTEGER(int_8), INTENT(INOUT)            :: optimal_cost
    INTEGER                                  :: my_process_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: my_cycles, my_moves, my_nbin, &
                                                my_ncpu, my_reduction_steps
    INTEGER(int_8)                           :: my_best_cost, &
                                                my_cost_optimal, my_cost_total
    INTEGER(int_8), DIMENSION(:), POINTER    :: my_bin_costs, my_cost_cpu
    INTEGER, DIMENSION(:), POINTER           :: my_best_distribution, &
                                                my_distribution
    REAL(dp)                                 :: my_initial_temperature, &
                                                my_reduction_factor, &
                                                my_swap_probability, &
                                                my_tolerance
    TYPE(rng_stream_type), POINTER           :: my_rng_stream

    CALL init(para_env, nbin_in, ncpu_in, bin_costs_in, load_balance_parameter, &
              my_nbin, my_ncpu, my_bin_costs, my_cost_cpu, my_cost_optimal, &
              my_cost_total, my_distribution, my_best_distribution, my_best_cost, my_rng_stream, &
              my_moves, my_cycles, my_reduction_steps, my_reduction_factor, my_initial_temperature, &
              my_tolerance, my_swap_probability, my_process_id, error)
    CALL distribute_monte_carlo(my_nbin, my_ncpu, my_bin_costs, my_cost_cpu, my_cost_optimal, &
                                my_cost_total, my_distribution, my_best_distribution, &
                                my_best_cost, my_rng_stream, my_moves, my_cycles, &
                                my_reduction_steps, my_reduction_factor, my_initial_temperature, &
                                my_tolerance, my_swap_probability,my_process_id, error)
    optimal_cost = my_best_cost
    distribution_vector = my_best_distribution
    CALL destructor(my_nbin, my_ncpu, my_bin_costs, my_cost_cpu, my_cost_optimal, &
                    my_cost_total, my_distribution, my_best_distribution, &
                    my_best_cost, my_rng_stream, my_moves, my_cycles, &
                    my_reduction_steps, my_reduction_factor, my_initial_temperature, &
                    my_tolerance, my_swap_probability, error)
  END SUBROUTINE optimize_distribution

END MODULE hfx_distribution_mc
