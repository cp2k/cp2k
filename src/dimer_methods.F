!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/dimer_methods [1.0] *
!!
!!   NAME
!!     dimer_methods
!!
!!   FUNCTION
!!     Contains types used for a Dimer Method calculations
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE dimer_methods
  USE cg_optimizer
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE dimer_types,                     ONLY: dimer_env_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get
  USE geo_opt
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE particle_list_types,             ONLY: particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_methods'

  PUBLIC :: cp_eval_at_ts

CONTAINS
  
!!****f* dimer_methods/cp_eval_at_ts [1.0] *
!!
!!   NAME
!!     cp_eval_at_ts
!!
!!   FUNCTION
!!     Computes the dimer energy/gradients (including the rotation of the dimer)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE  cp_eval_at_ts (gopt_env, x, f, gradient, calc_force, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x
    REAL(KIND=dp), INTENT(OUT)               :: f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient
    LOGICAL, INTENT(IN)                      :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at_ts', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure, fully_optimize, &
                                                interpolate_gradient
    REAL(KIND=dp)                            :: a0, a1, angle1, angle2, &
                                                angle_tol, b1, curvature0, &
                                                curvature1, curvature2, dCdp, &
                                                f1, norm_theta
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: g0, g1, g1p, nvec, theta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dimer_env_type), POINTER            :: dimer_env

    CALL cp_unimplemented_error(fromWhere=routineP, &
         message="DIMER method not yet working", &
         error=error, error_level=cp_failure_level)

    failure=.FALSE.
    NULLIFY(para_env,dimer_env)
    CALL timeset(routineN,"I","",handle)
    CPPostcondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,error,failure)
    dimer_env => gopt_env%dimer_env
    para_env  => gopt_env%force_env%para_env
    CPPostcondition(ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(nvec(SIZE(x)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(g0(SIZE(x)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(g1(SIZE(x)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(g1p(SIZE(x)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(theta(SIZE(x)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    ! Compute energy, gradients and rotation vector for R0
    CALL cp_eval_at_ts_low(gopt_env, x, 0, dimer_env, calc_force, f,  g0, error)

    ! Possibly rotate Dimer..
    IF (dimer_env%rotate_dimer) THEN
       ! Compute energy, gradients and rotation vector for R1
       CALL cp_eval_at_ts_low(gopt_env, x, 1, dimer_env, calc_force, f1, g1, error)    

       ! This loops optimize the direction of the vector dimer on the sphere
       DO WHILE (.TRUE.)
          ! Determine the theta vector 
          nvec       = dimer_env%nvec
          curvature0 = DOT_PRODUCT(g1-g0,dimer_env%nvec)/dimer_env%dr
          theta      = -2.0_dp*(g1-g0)
          theta      = theta - DOT_PRODUCT(theta, nvec) * nvec
          norm_theta = SQRT(DOT_PRODUCT(theta, theta))
          IF (norm_theta/=0.0_dp) theta = theta / norm_theta

          CALL rotate_dimer(dimer_env%nvec,theta,0.0_dp,error=error)

          ! Estimate of the rotational angle
          dCdp   =  2.0_dp*DOT_PRODUCT(g1-g0,theta)/dimer_env%dr
          b1     =  0.5_dp * dCdp
          angle1 = -0.5_dp*ATAN(dCdp/(2.0_dp*ABS(curvature0)))

          IF (angle1<angle_tol) EXIT
          ! Restore the dimer orientation before rotation
          dimer_env%nvec = nvec

          ! Rotating the dimer of dtheta degrees
          CALL rotate_dimer(dimer_env%nvec,theta,angle1,error=error)

          ! Re-compute energy, gradients and rotation vector for new R1
          CALL cp_eval_at_ts_low(gopt_env, x, 1, dimer_env, calc_force, f1, g1p, error)

          curvature1 = DOT_PRODUCT(g1p-g0,dimer_env%nvec)/dimer_env%dr
          a1     = (curvature0 - curvature1 + b1 * SIN(2.0_dp*angle1))/(1.0_dp - COS(2.0_dp*angle1))
          a0     = 2.0_dp*(curvature0-a1)
          dCdp   = 2.0_dp*DOT_PRODUCT(g1p-g0,theta)/dimer_env%dr
          angle2 = 0.5_dp*ATAN(b1/a1)
          curvature2 = a0/2.0_dp + a1 * COS(2.0_dp*angle2)+b1 * SIN(2.0_dp*angle2)
          IF (curvature2>curvature0) angle2 = angle2 + pi/2.0_dp
          dimer_env%nvec = nvec

          ! Rotating the dimer the optimized (in plane) vector position
          CALL rotate_dimer(dimer_env%nvec,theta,angle2,error=error)
          WRITE(*,*)"CURVATURE",curvature0,angle2,angle1,angle_tol
          IF (.NOT.fully_optimize) EXIT
          
          ! Interpolation of the gradient for the minimum point
          IF (interpolate_gradient) THEN
             ! Interpolate gradients for new minimum position
             g1 = SIN(angle1-angle2)/SIN(angle1)*g1+SIN(angle2)/SIN(angle1)*g1p+&
                  (1.0_dp-COS(angle2)-SIN(angle2)*TAN(angle1/2.0_dp))*g0
          ELSE
             ! If requested compute the value instead of performing an interpolation
             CALL cp_eval_at_ts_low(gopt_env, x, 1, dimer_env, calc_force, f1, g1, error)
          END IF
       END DO
    END IF
    STOP

    ! Correct the gradient
    IF (calc_force) THEN
       IF (curvature1>0) THEN
          gradient =  - DOT_PRODUCT(g0,dimer_env%nvec)*dimer_env%nvec
       ELSE
          gradient =  g0 - 2.0_dp*DOT_PRODUCT(g0,dimer_env%nvec)*dimer_env%nvec 
       END IF
    END IF
    
    STOP
    DEALLOCATE(theta, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(g1p, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(g1, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(g0, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(nvec, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_eval_at_ts

!***************************************************************************
!!****f* dimer_methods/cp_eval_at_ts_low [1.0] *
!!
!!   NAME
!!     cp_eval_at_ts_low
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE  cp_eval_at_ts_low (gopt_env, x,  dimer_index, dimer_env, calc_force,&
       f, gradient, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x
    INTEGER, INTENT(IN)                      :: dimer_index
    TYPE(dimer_env_type), POINTER            :: dimer_env
    LOGICAL, INTENT(IN)                      :: calc_force
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: gradient
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at_ts_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idg, idir, ip
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    idg=0
    CALL force_env_get(gopt_env%force_env, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys, particles=particles, error=error)
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          particles%els(ip)%r(idir)=x(idg)+REAL(dimer_index,KIND=dp)*dimer_env%nvec(idg)*dimer_env%dr
       END DO
    END DO

    !Compute energy and forces
    CALL force_env_calc_energy_force(gopt_env%force_env,calc_force=PRESENT(gradient),error=error)

    ! Possibly take the potential energy
    IF (PRESENT(f)) THEN
       CALL force_env_get(gopt_env%force_env,potential_energy=f,error=error)
    END IF
    
    ! Possibly take the gradients
    IF (PRESENT(gradient)) THEN
       idg=0
       CALL cp_subsys_get(subsys, particles=particles, error=error)
       DO ip=1,particles%n_els
          DO idir=1,3
             idg=idg+1
             CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
             gradient(idg)=-particles%els(ip)%f(idir)
          END DO
       END DO
    END IF

  END SUBROUTINE cp_eval_at_ts_low

!!****f* dimer_methods/rotate_dimer [1.0] *
!!
!!   NAME
!!     rotate_dimer
!!
!!   FUNCTION
!!     Performs a rotation of the unit dimer vector
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE rotate_dimer(nvec,theta,dt,error)
    REAL(KIND=dp), DIMENSION(:)              :: nvec, theta
    REAL(KIND=dp)                            :: dt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotate_dimer', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    nvec  = nvec * COS(dt) + theta * SIN(dt)

  END SUBROUTINE rotate_dimer

END MODULE dimer_methods
