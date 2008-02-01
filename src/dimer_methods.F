!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains types used for a Dimer Method calculations
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
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
  USE input_constants,                 ONLY: do_first_rotation_step,&
                                             do_second_rotation_step,&
                                             do_third_rotation_step
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

! *****************************************************************************
!> \brief Computes the dimer energy/gradients (including the rotation of the dimer)
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE  cp_eval_at_ts (gopt_env, x, f, gradient, calc_force, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    REAL(KIND=dp), INTENT(OUT)               :: f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient
    LOGICAL, INTENT(IN)                      :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at_ts', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: eval_analytical, failure
    REAL(KIND=dp)                            :: angle1, angle2, f1, norm
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
    IF (gopt_env%dimer_rotation) THEN
       SELECT CASE(dimer_env%rot%rotation_step)
       CASE(do_first_rotation_step,do_third_rotation_step)
          eval_analytical = .TRUE.
          IF ((dimer_env%rot%rotation_step==do_third_rotation_step).AND.(dimer_env%rot%interpolate_gradient)) THEN
             eval_analytical = .FALSE.
          END IF
          IF (eval_analytical) THEN
             ! Compute energy, gradients and rotation vector for R1
             CALL cp_eval_at_ts_low(gopt_env, x, 1, dimer_env, calc_force, f1, dimer_env%rot%g1, error)    
          ELSE
             angle1 = dimer_env%rot%angle1
             angle2 = dimer_env%rot%angle2
             dimer_env%rot%g1 = SIN(angle1-angle2)/SIN(angle1)*dimer_env%rot%g1+&
                            SIN(angle2)/SIN(angle1)*dimer_env%rot%g1p+&
                            (1.0_dp-COS(angle2)-SIN(angle2)*TAN(angle1/2.0_dp))*dimer_env%rot%g0             
          END IF

          ! Determine the theta vector (i.e. the search direction for line minimization) 
          gradient = -2.0_dp*(dimer_env%rot%g1-dimer_env%rot%g0)
          gradient = gradient  - DOT_PRODUCT(gradient, dimer_env%nvec) * dimer_env%nvec
          norm     = SQRT(DOT_PRODUCT(gradient, gradient))
          f        = norm
          dimer_env%cg_rot%norm_theta_old = dimer_env%cg_rot%norm_theta
          dimer_env%cg_rot%norm_theta     = norm
          IF (norm/=0.0_dp) gradient = gradient/norm

          ! Compute curvature and derivative of the curvature w.r.t. the rotational angle
          dimer_env%rot%curvature  = DOT_PRODUCT(dimer_env%rot%g1-dimer_env%rot%g0,dimer_env%nvec)/dimer_env%dr
          dimer_env%rot%dCdp       = 2.0_dp*DOT_PRODUCT(dimer_env%rot%g1-dimer_env%rot%g0,gradient)/dimer_env%dr

          dimer_env%rot%rotation_step = do_second_rotation_step
          gradient = -gradient
       CASE(do_second_rotation_step)
          ! Compute energy, gradients and rotation vector for R1
          CALL cp_eval_at_ts_low(gopt_env, x, 1, dimer_env, calc_force, f1, dimer_env%rot%g1p, error)
          dimer_env%rot%curvature  = DOT_PRODUCT(dimer_env%rot%g1p-dimer_env%rot%g0,dimer_env%nvec)/dimer_env%dr
          dimer_env%rot%rotation_step = do_third_rotation_step

          ! Determine the theta vector (i.e. the search direction for line minimization)
          ! This is never used for getting a new theta but may is consistent in order to
          ! give back the right value of f
          gradient = -2.0_dp*(dimer_env%rot%g1p-dimer_env%rot%g0)
          gradient = gradient  - DOT_PRODUCT(gradient, dimer_env%nvec) * dimer_env%nvec
          norm     = SQRT(DOT_PRODUCT(gradient, gradient))
          f        = norm
          IF (norm/=0.0_dp) gradient = gradient/norm
       END SELECT
    ELSE
       ! Compute energy, gradients and rotation vector for R0
       CALL cp_eval_at_ts_low(gopt_env, x, 0, dimer_env, calc_force, f,  dimer_env%rot%g0, error)
       
       ! The dimer is rotated only when we are out of the translation line search
       IF (.NOT.gopt_env%do_line_search) THEN
          CALL cp_rot_opt(gopt_env%gopt_dimer_env, x, gopt_env%gopt_dimer_param ,&
               gopt_env%gopt_dimer_env%geo_section, error)
          dimer_env%rot%rotation_step = do_first_rotation_step
       END IF

       ! Correcting gradients for Translation
       IF (dimer_env%rot%curvature>0) THEN
          gradient =  - DOT_PRODUCT(dimer_env%rot%g0,dimer_env%nvec)*dimer_env%nvec
          CALL remove_translation_component(gopt_env, x, gradient, error)
       ELSE
          gradient =  dimer_env%rot%g0 - 2.0_dp*DOT_PRODUCT(dimer_env%rot%g0,dimer_env%nvec)*dimer_env%nvec 
          CALL remove_translation_component(gopt_env, x, gradient, error)
       END IF

       IF (.NOT.gopt_env%do_line_search) THEN
          f =  SQRT(DOT_PRODUCT(gradient,gradient))
       ELSE
          f = -DOT_PRODUCT(gradient,dimer_env%tsl%tls_vec)
       END IF
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_eval_at_ts

! *****************************************************************************
!> \brief This function removes translational forces after project of the gradient
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  SUBROUTINE  remove_translation_component (gopt_env,x,gradient,error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    REAL(KIND=dp), DIMENSION(:)              :: gradient
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_translation_component', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idg, idir, ip
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: norm_gradient_new, &
                                                norm_gradient_old
    REAL(KIND=dp), DIMENSION(3)              :: tot_force
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    CALL force_env_get(gopt_env%force_env, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys, particles=particles, error=error)
    tot_force = 0.0_dp
    idg = 0
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          tot_force(idir)=tot_force(idir)+gradient(idg)
       END DO
    END DO
    tot_force=tot_force/REAL(particles%n_els,KIND=dp)
    ! Remove spurius translational forces
    norm_gradient_old = SQRT(DOT_PRODUCT(gradient,gradient))
    idg=0
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          gradient(idg)=gradient(idg)-tot_force(idir)
       END DO
    END DO
    norm_gradient_new = SQRT(DOT_PRODUCT(gradient,gradient))
    CPPostcondition(norm_gradient_new>0.0_dp,cp_failure_level,routineP,error,failure)
    gradient = (gradient/norm_gradient_new)*norm_gradient_old

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE remove_translation_component

! *****************************************************************************
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
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

    INTEGER                                  :: handle, idg, idir, ip
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    CALL timeset(routineN,"I","",handle)
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
    CALL force_env_calc_energy_force(gopt_env%force_env,calc_force=calc_force,error=error)

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
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_eval_at_ts_low

END MODULE dimer_methods
