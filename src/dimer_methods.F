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
  USE geo_opt,                         ONLY: cp_rot_opt
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE input_constants,                 ONLY: default_minimization_method_id
  USE kinds,                           ONLY: dp
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
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    REAL(KIND=dp), INTENT(OUT)               :: f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient
    LOGICAL, INTENT(IN)                      :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at_ts', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: f1, norm
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dimer_env_type), POINTER            :: dimer_env

    failure=.FALSE.
    NULLIFY(para_env,dimer_env)
    CALL timeset(routineN,"I","",handle)
    CPPostcondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,error,failure)
    dimer_env => gopt_env%dimer_env
    para_env  => gopt_env%force_env%para_env
    CPPostcondition(ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)

    ! Possibly rotate Dimer or just compute Gradients of point 0 for Translation
    IF (gopt_env%type_id==default_minimization_method_id) THEN

       ! Compute energy, gradients and rotation vector for R1
       CALL cp_eval_at_ts_low(gopt_env, x, 1, dimer_env, calc_force, f1, dimer_env%g1, error)    

       ! Determine the theta vector (i.e. the search direction for line minimization) 
       gradient = -2.0_dp*(dimer_env%g1-dimer_env%g0)
       gradient = gradient  - DOT_PRODUCT(gradient, dimer_env%nvec) * dimer_env%nvec
       norm     = SQRT(DOT_PRODUCT(gradient, gradient))
       f        = norm
       IF (norm/=0.0_dp) gradient = gradient/norm

       ! Compute curvature and derivative of the curvature w.r.t. the rotational angle
       dimer_env%curvature  = DOT_PRODUCT(dimer_env%g1-dimer_env%g0,dimer_env%nvec)/dimer_env%dr
       dimer_env%dCdp       = 2.0_dp*DOT_PRODUCT(dimer_env%g1-dimer_env%g0,gradient)/dimer_env%dr
    ELSE
       ! Compute energy, gradients and rotation vector for R0
       CALL cp_eval_at_ts_low(gopt_env, x, 0, dimer_env, calc_force, f,  dimer_env%g0, error)

       IF (dimer_env%rotate_dimer) THEN
          CALL cp_rot_opt(gopt_env%gopt_dimer_env, x, gopt_env%gopt_dimer_param ,&
               gopt_env%gopt_dimer_env%geo_section, error)
       END IF
       ! Correcting gradients for Translation
       IF (calc_force) THEN
          IF (dimer_env%curvature>0) THEN
             gradient =  - DOT_PRODUCT(dimer_env%g0,dimer_env%nvec)*dimer_env%nvec
          ELSE
             gradient =  dimer_env%g0 - 2.0_dp*DOT_PRODUCT(dimer_env%g0,dimer_env%nvec)*dimer_env%nvec 
          END IF
       END IF
       
       ! Fix the pseudo-energy
       f = HUGE(0.0_dp)
    END IF

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
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
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
