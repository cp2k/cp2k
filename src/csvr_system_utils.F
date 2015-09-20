!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      15.10.2007 Giovanni Bussi - Implementation validated.
!> \author Teodoro Laino - 09.2007 - University of Zurich [tlaino]
! *****************************************************************************
MODULE csvr_system_utils

  
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
 
  PRIVATE

  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.
  PUBLIC                               :: rescaling_factor
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_utils'

CONTAINS

! *****************************************************************************
!> \brief Stochastic velocity rescale, as described in
!>      Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!>
!>      This subroutine implements Eq.(A7) and returns the new value for the kinetic energy,
!>      which can be used to rescale the velocities.
!>      The procedure can be applied to all atoms or to smaller groups.
!>      If it is applied to intersecting groups in sequence, the kinetic energy
!>      that is given as an input (kk) has to be up-to-date with respect to the previous
!>      rescalings.
!>
!>      When applied to the entire system, and when performing standard molecular dynamics
!>      (fixed c.o.m. (center of mass))
!>      the degrees of freedom of the c.o.m. have to be discarded in the calculation of ndeg,
!>      and the c.o.m. momentum HAS TO BE SET TO ZERO.
!>      When applied to subgroups, one can chose to:
!>      (a) calculate the subgroup kinetic energy in the usual reference frame, and count
!>          the c.o.m. in ndeg
!>      (b) calculate the subgroup kinetic energy with respect to its c.o.m. motion, discard
!>          the c.o.m. in ndeg and apply the rescale factor with respect to the subgroup c.o.m.
!>          velocity.
!>      They should be almost equivalent.
!>      If the subgroups are expected to move one respect to the other, the choice (b)
!>      should be better.
!>
!>      If a null relaxation time is required (taut=0.0), the procedure reduces to an istantaneous
!>      randomization of the kinetic energy, as described in paragraph IIA.
!>
!>      HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
!>      The effective-energy (htilde) drift can be used to check the integrator against
!>      discretization errors.
!>      The easiest recipe is:
!>      htilde = h + conint
!>      where h is the total energy (kinetic + potential)
!>      and conint is a quantity accumulated along the trajectory as minus the sum of all
!>      the increments of kinetic energy due to the thermostat.
!>
!>      Variables:
!>       kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
!>       sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
!>       ndeg  ! number of degrees of freedom of the atoms to be thermalized
!>       taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
!> \param kk ...
!> \param sigma ...
!> \param ndeg ...
!> \param taut ...
!> \param rng_stream ...
!> \retval my_res ...
!> \date 09.2007
!> \author Giovanni Bussi - ETH Zurich, Lugano 10.2007
! *****************************************************************************
  FUNCTION rescaling_factor(kk, sigma, ndeg, taut, rng_stream) RESULT(my_res)
    REAL(KIND=dp), INTENT(IN)                :: kk, sigma
    INTEGER, INTENT(IN)                      :: ndeg
    REAL(KIND=dp), INTENT(IN)                :: taut
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp)                            :: my_res

    CHARACTER(len=*), PARAMETER :: routineN = 'rescaling_factor', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: factor, resample, reverse, rr

    my_res=0.0_dp
    IF (kk > 0.0_dp) THEN
       IF(taut>0.1_dp) THEN
          factor=EXP(-1.0_dp/taut)
       ELSE
          factor=0.0_dp
       END IF
       rr = next_random_number(rng_stream)
       reverse=1.0_dp
       ! reverse of momentum is implemented to have the correct limit to Langevin dynamics for ndeg=1
       IF(rr<-SQRT(ndeg*kk*factor/(sigma*(1.0_dp-factor)))) reverse=-1.0_dp
       ! for ndeg/=1, the reverse of momentum is not necessary. in principles, it should be there.
       ! in practice, it is better to skip it to avoid unnecessary slowing down of the dynamics in the small taut regime
       ! anyway, this should not affect the final ensemble
       IF(ndeg/=1) reverse=1.0_dp
       resample = kk + (1.0_dp-factor)* (sigma*(sumnoises(ndeg-1, rng_stream)+rr**2)/REAL(ndeg,KIND=dp)-kk) &
                     +  2.0_dp*rr*SQRT(kk*sigma/ndeg*(1.0_dp-factor)*factor)

       resample=MAX(0.0_dp,resample)
       my_res = reverse*SQRT(resample/kk)
    END IF
  END FUNCTION rescaling_factor

! *****************************************************************************
!> \brief returns the sum of n independent gaussian noises squared
!>      (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
!> \param nn ...
!> \param rng_stream ...
!> \retval sum_gauss ...
!> \date 09.2007
!> \author Teo - University of Zurich
! *****************************************************************************
  FUNCTION sumnoises(nn, rng_stream) RESULT(sum_gauss)
    INTEGER, INTENT(IN)                      :: nn
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp)                            :: sum_gauss

    CHARACTER(len=*), PARAMETER :: routineN = 'sumnoises', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    sum_gauss = 0.0_dp
    DO i = 1, nn
       sum_gauss = sum_gauss + next_random_number(rng_stream)**2
    END DO

  END FUNCTION sumnoises

END MODULE csvr_system_utils
