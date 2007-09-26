!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/x_hf_distribution_mc[1.0] *
!!
!!   NAME
!!     x_hf_distribution_mc 
!!
!!   FUNCTION
!!     Optimizes the load balance for a given cost matrix for a one dimensional
!!     distribution of processors
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!****************************************************************************

MODULE x_hf_distribution_mc

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             reset_to_next_rng_substream,&
                                             rng_stream_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC     x_hf_optimize_distribution

  !! Module Variables
  INTEGER                                    :: my_nbin
  INTEGER                                    :: my_ncpu
  INTEGER(int_8), ALLOCATABLE, DIMENSION(:)  :: my_bin_costs
  INTEGER(int_8), ALLOCATABLE, DIMENSION(:)  :: my_cost_cpu
  INTEGER(int_8)                             :: my_cost_optimal
  INTEGER(int_8)                             :: my_cost_total
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: my_distribution
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: my_best_distribution
  INTEGER(int_8)                             :: my_best_cost
  TYPE(rng_stream_type), POINTER             :: my_rng_stream

 
  CONTAINS


!!***** x_hf_distribution_mc/init [1.0] *
!!
!!   NAME
!!     init 
!!
!!   FUNCTION
!!     Initializes the module variables
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - para_env: parallel environnment
!!     - nbin_in: size of the cost_matrix
!!     - ncpu_in: number of cpus
!!     - bin_costs_in: cost per (nbin*ncpu)
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE init(para_env, nbin_in,ncpu_in, bin_costs_in, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: nbin_in, ncpu_in
    INTEGER(int_8), DIMENSION(:), POINTER    :: bin_costs_in
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i

    my_nbin = nbin_in
    my_ncpu = ncpu_in

    ALLOCATE(my_cost_cpu(my_ncpu))
    ALLOCATE(my_bin_costs(my_nbin))
    ALLOCATE(my_distribution(my_nbin))
    ALLOCATE(my_best_distribution(my_nbin))

    my_bin_costs = bin_costs_in
    my_cost_total=SUM(my_bin_costs)
    my_cost_optimal=(SUM(my_bin_costs)+(my_ncpu)-1)/my_ncpu

    ! initialize the random number generators
    ! each CPU should have a different stream of random numbers
    ! so the optimization problem is parallellized (Ncpu independent searches)
    NULLIFY(my_rng_stream)
    CALL create_rng_stream(my_rng_stream,"distribute_1d_monte_carlo",distribution_type=UNIFORM,error=error)
    DO i=1,para_env%mepos
       CALL reset_to_next_rng_substream(my_rng_stream,error=error)
    ENDDO

  END SUBROUTINE init
!***************************************************************************

!!***** x_hf_distribution_mc/destructor [1.0] *
!!
!!   NAME
!!     destructor 
!!
!!   FUNCTION
!!     frees memory
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE destructor(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

     DEALLOCATE(my_cost_cpu)
     DEALLOCATE(my_bin_costs)
     DEALLOCATE(my_distribution)
     DEALLOCATE(my_best_distribution)

     CALL delete_rng_stream(my_rng_stream,error=error)
  END SUBROUTINE destructor
!***************************************************************************

!!***** x_hf_distribution_mc/compute_cost_cpu [1.0] *
!!
!!   NAME
!!     compute_cost_cpu 
!!
!!   FUNCTION
!!     Calculate the cost of a distribution
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - distribution: maps cost_matrix to a cpu
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE compute_cost_cpu(distribution)
    INTEGER                                  :: distribution(my_nbin)

    INTEGER                                  :: i

    my_cost_cpu = 0
    DO i=1,my_nbin
      my_cost_cpu(distribution(i)) = my_cost_cpu(distribution(i))+my_bin_costs(i)
    ENDDO
  END SUBROUTINE compute_cost_cpu
!***************************************************************************

!!***** x_hf_distribution_mc/compute_max_cost [1.0] *
!!
!!   NAME
!!     compute_max_cost
!!
!!   FUNCTION
!!     Calculates the maximal cost of a distribution
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - max_cost: variable that stores the max_cost
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE compute_max_cost(max_cost)
    INTEGER(int_8)                           :: max_cost

    INTEGER                                  :: i

    max_cost=-HUGE(max_cost)
    DO i=1,my_ncpu
      IF (my_cost_cpu(i).GT.max_cost) THEN 
        max_cost=my_cost_cpu(i)
      ENDIF
    ENDDO

  END SUBROUTINE compute_max_cost
!***************************************************************************

!!***** x_hf_distribution_mc/distribute_monte_carlo [1.0] *
!!
!!   NAME
!!     distribute_monte_carlo
!!
!!   FUNCTION
!!     Optimizes a 1d distribution according to a given cost_matrix
!!
!!   NOTES
!!     Uses two basic random updates:
!!        a) swap the cost of two distributions
!!        b) assign the cost of i to j
!!
!!   INPUTS
!!     - Ntries: Number of trial steps
!!     - T_in: Initial temperature (for acceptance criterion)
!!     - error:  variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE distribute_monte_carlo(Ntries,T_in, error)
    INTEGER                                  :: Ntries
    REAL(dp)                                 :: T_in
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER :: best_iac, i, iac, iswap, itmp, itries, jswap, &
      my_distribution(my_nbin), old_distribution(my_nbin)
    INTEGER(int_8)                           :: best_max_cost, max_cost, &
                                                max_cost_old
    REAL(dp)                                 :: random_accept, random_move, &
                                                random_move_type, T, t_final

!! Initialize a modulo distribution

    DO i=1,my_nbin
      my_distribution(i)=MOD(i-1,my_ncpu)+1
    ENDDO

    CALL compute_cost_cpu(my_distribution)
    CALL compute_max_cost(max_cost)

    best_max_cost = max_cost
    my_best_distribution = my_distribution
    max_cost_old = max_cost
    old_distribution = my_distribution
    
    iac=0
    best_iac=0

    DO itries=0,Ntries
      !do annealing 
      IF (MOD(itries,Ntries/10)==0) THEN
        T=5*T_in
      ENDIF
      IF (MOD(itries,Ntries/100)==0) THEN
        T=T*0.8
      ENDIF
      ! do some Monte Carlo
      random_move = next_random_number(my_rng_stream,error=error)
      random_move_type = next_random_number(my_rng_stream,error=error)
      ! generate trial configuration
      my_distribution = old_distribution
      IF (random_move<0.1) THEN
        IF (random_move_type<0.5) THEN
          iswap = next_random_number(my_rng_stream,low=1,high=my_nbin,error=error)
          jswap = next_random_number(my_rng_stream,low=1,high=my_ncpu,error=error)
          IF (iswap>my_nbin .OR. jswap>my_ncpu) CYCLE
            my_distribution(iswap)=jswap
          ENDIF
        ELSE
         ! randomly swap a col/row with 0.9 probability
         IF (random_move_type<0.5) THEN
          iswap = next_random_number(my_rng_stream,low=1,high=my_nbin,error=error)
          jswap = next_random_number(my_rng_stream,low=1,high=my_ncpu,error=error)
          IF (iswap>my_nbin .OR. jswap>my_nbin) CYCLE
          itmp=my_distribution(iswap)
          my_distribution(iswap)=my_distribution(jswap)
          my_distribution(jswap)=itmp
        ENDIF
      ENDIF

      ! evaluate its cost
      iac = iac + 1
       
      CALL compute_cost_cpu(my_distribution) 
      CALL compute_max_cost(max_cost)

      ! accept / reject
      IF (max_cost < max_cost_old .OR. (EXP( -(max_cost - max_cost_old)/T) .GT. random_accept)) THEN
        max_cost_old = max_cost
        old_distribution = my_distribution
      ELSE
        ! restore the old position
        CALL compute_cost_cpu(old_distribution) 
      ENDIF

      ! if this one happens to be the best so far, store this
      IF (max_cost < best_max_cost) THEN
        best_max_cost = max_cost
        my_best_distribution = my_distribution
        best_iac = iac
        !
        ! if we are within 0.001 of the optimal cost, stop optimizing, or cost model is certainly not that good
        !
      ENDIF

      IF (REAL((max_cost-my_cost_optimal),dp)/my_cost_optimal.LE.0.001) EXIT
    ENDDO
    CALL CPU_TIME(t_final)

    CALL compute_cost_cpu(my_best_distribution)
    my_best_cost = best_max_cost

  END SUBROUTINE distribute_monte_carlo

!***************************************************************************

!!***** x_hf_distribution_mc/x_hf_optimize_distribution [1.0] *
!!
!!   NAME
!!     x_hf_optimize_distribution 
!!
!!   FUNCTION
!!     Interface to the module
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: parallel environment
!!     - nbin_in: size of the cost matrix
!!     - ncpu_in: number of processors
!!     - bin_costs_in: cost_matrix
!!     - distribution_vector: 1d vector that assigns cost matrix to process
!!     - steps: number of Monte Carlo steps to be performed
!!     - initial_temperature: initial temperature for acceptance criterion
!!     - optimal_cost: stores the optimal cost
!!     - error:  variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_optimize_distribution(para_env, nbin_in, ncpu_in, bin_costs_in, distribution_vector, &
                                steps, initial_temperature, optimal_cost, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: nbin_in, ncpu_in
    INTEGER(int_8), DIMENSION(:), POINTER    :: bin_costs_in
    INTEGER, DIMENSION(:), POINTER           :: distribution_vector
    INTEGER, INTENT(IN)                      :: steps
    REAL(dp), INTENT(INOUT)                  :: initial_temperature
    INTEGER(int_8), INTENT(INOUT)            :: optimal_cost
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL init(para_env, nbin_in, ncpu_in, bin_costs_in, error)
    CALL distribute_monte_carlo(steps, initial_temperature, error)
    
    optimal_cost = my_best_cost
    distribution_vector = my_best_distribution
    
    CALL destructor(error)

  END SUBROUTINE x_hf_optimize_distribution

END MODULE x_hf_distribution_mc
