!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Geometry optimization using  Conjugate Gradients
!> \author Teodoro Laino [teo]
!>      10.2005
! *****************************************************************************
MODULE cg_utils
  USE cp_external_control,             ONLY: external_control
  USE dimer_types,                     ONLY: dimer_env_type
  USE dimer_utils,                     ONLY: dimer_thrs,&
                                             rotate_dimer
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_cell_method_id,&
                                             default_minimization_method_id,&
                                             default_shellcore_method_id,&
                                             default_ts_method_id,&
                                             ls_2pnt,&
                                             ls_fit,&
                                             ls_gold
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  PUBLIC :: cg_linmin, get_conjugate_direction
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cg_utils'

CONTAINS

! *****************************************************************************
!> \brief Main driver for line minimization routines for CG
!> \param gopt_env ...
!> \param xvec ...
!> \param xi ...
!> \param g ...
!> \param opt_energy ...
!> \param output_unit ...
!> \param gopt_param ...
!> \param globenv ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE cg_linmin(gopt_env, xvec, xi, g, opt_energy, output_unit, gopt_param,&
       globenv)

    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi, g
    REAL(KIND=dp), INTENT(INOUT)             :: opt_energy
    INTEGER                                  :: output_unit
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_linmin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL timeset(routineN,handle)
    gopt_env%do_line_search = .TRUE.
    SELECT CASE(gopt_env%type_id)
    CASE(default_minimization_method_id)
       SELECT CASE(gopt_param%cg_ls%type_id)
       CASE(ls_2pnt)
          CALL linmin_2pnt(gopt_env, xvec, xi, g, opt_energy, gopt_param, use_only_grad=gopt_param%cg_ls%grad_only, &
               output_unit=output_unit)
       CASE(ls_fit)
          CALL linmin_fit(gopt_env, xvec, xi, opt_energy, gopt_param%cg_ls%brack_limit, &
               gopt_param%cg_ls%initial_step, output_unit, gopt_param, globenv)
       CASE(ls_gold)
          CALL linmin_gold(gopt_env, xvec, xi, opt_energy, gopt_param%cg_ls%brent_tol,&
               gopt_param%cg_ls%brent_max_iter, gopt_param%cg_ls%brack_limit, &
               gopt_param%cg_ls%initial_step, output_unit, globenv)
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Line Search type not yet implemented in CG.", &
               error_level=cp_failure_level)
       END SELECT
    CASE(default_ts_method_id)
       SELECT CASE(gopt_param%cg_ls%type_id)
       CASE(ls_2pnt)
          IF (gopt_env%dimer_rotation) THEN
             CALL rotmin_2pnt(gopt_env, gopt_env%dimer_env, xvec, xi, opt_energy)
          ELSE
             CALL tslmin_2pnt(gopt_env, gopt_env%dimer_env, xvec, xi, opt_energy, gopt_param,&
                  output_unit)
          END IF
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Line Search type not yet implemented in CG for TS search.", &
               error_level=cp_failure_level)
       END SELECT
    CASE(default_cell_method_id)
       SELECT CASE(gopt_param%cg_ls%type_id)
       CASE(ls_2pnt)
          CALL linmin_2pnt(gopt_env, xvec, xi, g, opt_energy, gopt_param, use_only_grad=.TRUE.,&
               output_unit=output_unit)
       CASE(ls_gold)
          CALL linmin_gold(gopt_env, xvec, xi, opt_energy, gopt_param%cg_ls%brent_tol,&
               gopt_param%cg_ls%brent_max_iter, gopt_param%cg_ls%brack_limit, &
               gopt_param%cg_ls%initial_step, output_unit, globenv)
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Line Search type not yet implemented in CG for cell optimization.", &
               error_level=cp_failure_level)
       END SELECT
    CASE(default_shellcore_method_id)
       SELECT CASE(gopt_param%cg_ls%type_id)
       CASE(ls_2pnt)
          CALL linmin_2pnt(gopt_env, xvec, xi, g, opt_energy, gopt_param, use_only_grad=.TRUE.,&
               output_unit=output_unit)
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Line Search type not yet implemented in CG for shellcore optimization.", &
               error_level=cp_failure_level)
       END SELECT

    END SELECT
    gopt_env%do_line_search = .FALSE.
    CALL timestop(handle)

  END SUBROUTINE cg_linmin

! *****************************************************************************
!> \brief Line search subroutine based on 2 points (using gradients and energies
!>        or only gradients)
!> \param gopt_env ...
!> \param x0 ...
!> \param ls_vec ...
!> \param g ...
!> \param opt_energy ...
!> \param gopt_param ...
!> \param use_only_grad ...
!> \param output_unit ...
!> \author Teodoro Laino - created [tlaino] - 03.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE linmin_2pnt(gopt_env, x0, ls_vec, g, opt_energy, gopt_param, use_only_grad,&
       output_unit)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0, ls_vec, g
    REAL(KIND=dp), INTENT(INOUT)             :: opt_energy
    TYPE(gopt_param_type), POINTER           :: gopt_param
    LOGICAL, INTENT(IN), OPTIONAL            :: use_only_grad
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'linmin_2pnt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure, my_use_only_grad, &
                                                save_consistent_energy_force
    REAL(KIND=dp) :: a, b, c, dx, dx_min, dx_min_save, dx_thrs, norm_grad1, &
      norm_grad2, norm_ls_vec, opt_energy2, x_grad_zero
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient2, ls_norm

    failure = .FALSE.
    CALL timeset(routineN,handle)
    norm_ls_vec = SQRT(DOT_PRODUCT(ls_vec,ls_vec))
    my_use_only_grad = .FALSE.
    IF (PRESENT(use_only_grad)) my_use_only_grad = use_only_grad
    IF (norm_ls_vec/=0.0_dp) THEN
       ALLOCATE(ls_norm(SIZE(ls_vec)),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
       ALLOCATE(gradient2(SIZE(ls_vec)),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
       ls_norm = ls_vec / norm_ls_vec
       dx      = norm_ls_vec
       dx_thrs = gopt_param%cg_ls%max_step

       x0 = x0 + dx * ls_norm
       ![NB] don't need consistent energies and forces if using only gradient
       save_consistent_energy_force = gopt_env%require_consistent_energy_force
       gopt_env%require_consistent_energy_force = .NOT. my_use_only_grad
       CALL cp_eval_at(gopt_env, x0, opt_energy2, gradient2, master=gopt_env%force_env%para_env%mepos,&
                       final_evaluation=.FALSE., para_env=gopt_env%force_env%para_env)
       gopt_env%require_consistent_energy_force = save_consistent_energy_force

       norm_grad1 = -DOT_PRODUCT(        g,ls_norm)
       norm_grad2 =  DOT_PRODUCT(gradient2,ls_norm)
       IF (my_use_only_grad) THEN
          ! a*x+b=y
          ! per x=0; b=norm_grad1
          b = norm_grad1
          ! per x=dx; a*dx+b=norm_grad2
          a = (norm_grad2 - b)/dx
          x_grad_zero= -b/a
          dx_min = x_grad_zero
       ELSE
          ! ax**2+b*x+c=y
          ! per x=0 ; c=opt_energy
          c = opt_energy
          ! per x=dx;          a*dx**2 + b*dx + c = opt_energy2
          ! per x=dx;        2*a*dx    + b        = norm_grad2
          !
          !                  - a*dx**2        + c = (opt_energy2-norm_grad2*dx)
          !                    a*dx**2            = c - (opt_energy2-norm_grad2*dx)
          a = (c - (opt_energy2-norm_grad2*dx))/dx**2
          b = norm_grad2 - 2.0_dp*a*dx
          dx_min = 0.0_dp
          IF (a/=0.0_dp) dx_min = -b/(2.0_dp*a)
          opt_energy = opt_energy2
       END IF
       dx_min_save = dx_min
       ! In case the solution is larger than the maximum threshold let's assume the maximum allowed
       ! step length
       IF (ABS(dx_min)>dx_thrs)  dx_min = SIGN(1.0_dp,dx_min)*dx_thrs
       x0 = x0 + (dx_min-dx) * ls_norm

       ! Print out LS info
       IF (output_unit>0) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T31,A,T78,A)")&
               "***","2PNT LINE SEARCH INFO","***"
          WRITE(UNIT=output_unit,FMT="(T2,A,T78,A)")"***","***"
          WRITE(UNIT=output_unit,FMT="(T2,A,3X,A,F12.6,T45,A,F12.6,T78,A)")&
               "***","DX (EVALUATED)=",dx,"DX (THRESHOLD)=",dx_thrs,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A,3X,A,F12.6,T45,A,F12.6,T78,A)")&
               "***","DX (FITTED   )=",dx_min_save,"DX (ACCEPTED )=",dx_min,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       END IF
       DEALLOCATE(ls_norm,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
       DEALLOCATE(gradient2,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    ELSE
       ! Do Nothing, since.. if the effective force is 0 means that we are already
       ! in the saddle point..
    END IF
    CALL timestop(handle)
  END SUBROUTINE linmin_2pnt

! *****************************************************************************
!> \brief Translational minimization for the Dimer Method - 2pnt LS
!> \param gopt_env ...
!> \param dimer_env ...
!> \param x0 ...
!> \param tls_vec ...
!> \param opt_energy ...
!> \param gopt_param ...
!> \param output_unit ...
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  SUBROUTINE tslmin_2pnt(gopt_env, dimer_env, x0, tls_vec, opt_energy, gopt_param, output_unit)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(dimer_env_type), POINTER            :: dimer_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0, tls_vec
    REAL(KIND=dp), INTENT(INOUT)             :: opt_energy
    TYPE(gopt_param_type), POINTER           :: gopt_param
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'tslmin_2pnt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, dx_min, dx_min_acc, &
                                                dx_min_save, dx_thrs, &
                                                norm_tls_vec, opt_energy2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tls_norm

    failure = .FALSE.
    CALL timeset(routineN,handle)
    norm_tls_vec = SQRT(DOT_PRODUCT(tls_vec,tls_vec))
    IF (norm_tls_vec/=0.0_dp) THEN
       ALLOCATE(tls_norm(SIZE(tls_vec)),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)

       tls_norm = tls_vec / norm_tls_vec
       dimer_env%tsl%tls_vec => tls_norm

       dx      = norm_tls_vec
       dx_thrs = gopt_param%cg_ls%max_step
       ! If curvature is positive let's make the largest step allowed
       IF (dimer_env%rot%curvature>0) dx = dx_thrs
       x0 = x0 + dx * tls_norm
       CALL cp_eval_at(gopt_env, x0, opt_energy2, master=gopt_env%force_env%para_env%mepos,&
                       final_evaluation=.FALSE., para_env=gopt_env%force_env%para_env)
       IF (dimer_env%rot%curvature>0) THEN
          dx_min      = 0.0_dp
          dx_min_save = dx
          dx_min_acc  = dx
       ELSE
          ! First let's try to interpolate the minimum
          dx_min = -opt_energy/(opt_energy2-opt_energy)*dx
          ! In case the solution is larger than the maximum threshold let's assume the maximum allowed
          ! step length
          dx_min_save = dx_min
          IF (ABS(dx_min)>dx_thrs)  dx_min = SIGN(1.0_dp,dx_min)*dx_thrs
          dx_min_acc  = dx_min
          dx_min      = dx_min - dx
       END IF
       x0 = x0 + dx_min * tls_norm

       ! Print out LS info
       IF (output_unit>0) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T24,A,T78,A)")&
               "***","2PNT TRANSLATIONAL LINE SEARCH INFO","***"
          WRITE(UNIT=output_unit,FMT="(T2,A,T78,A)")"***","***"
          WRITE(UNIT=output_unit,FMT="(T2,A,3X,A,F12.6,T78,A)")&
               "***","LOCAL CURVATURE =",dimer_env%rot%curvature,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A,3X,A,F12.6,T45,A,F12.6,T78,A)")&
               "***","DX (EVALUATED)=",dx,"DX (THRESHOLD)=",dx_thrs,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A,3X,A,F12.6,T45,A,F12.6,T78,A)")&
               "***","DX (FITTED   )=",dx_min_save,"DX (ACCEPTED )=",dx_min_acc,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       END IF

       ! Here we compute the value of the energy in point zero..
       CALL cp_eval_at(gopt_env, x0, opt_energy, master=gopt_env%force_env%para_env%mepos,&
                       final_evaluation=.FALSE., para_env=gopt_env%force_env%para_env)

       DEALLOCATE(tls_norm,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    ELSE
       ! Do Nothing, since.. if the effective force is 0 means that we are already
       ! in the saddle point..
    END IF
    CALL timestop(handle)

  END SUBROUTINE tslmin_2pnt

! *****************************************************************************
!> \brief Rotational minimization for the Dimer Method - 2 pnt LS
!> \param gopt_env ...
!> \param dimer_env ...
!> \param x0 ...
!> \param theta ...
!> \param opt_energy ...
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  SUBROUTINE rotmin_2pnt(gopt_env, dimer_env, x0, theta, opt_energy)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(dimer_env_type), POINTER            :: dimer_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0, theta
    REAL(KIND=dp), INTENT(INOUT)             :: opt_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'rotmin_2pnt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a0, a1, angle, b1, &
                                                curvature0, curvature1, &
                                                curvature2, dCdp, f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work

    failure = .FALSE.
    CALL timeset(routineN,handle)
    curvature0 =  dimer_env%rot%curvature
    dCdp       =  dimer_env%rot%dCdp
    b1         =  0.5_dp*dCdp
    angle      = -0.5_dp*ATAN(dCdp/(2.0_dp*ABS(curvature0)))
    dimer_env%rot%angle1 = angle
    dimer_env%cg_rot%nvec_old = dimer_env%nvec
    IF (angle>dimer_env%rot%angle_tol) THEN
       ! Rotating the dimer of dtheta degrees
       CALL rotate_dimer(dimer_env%nvec,theta,angle)
       ! Re-compute energy, gradients and rotation vector for new R1
       CALL cp_eval_at(gopt_env, x0, f, master=gopt_env%force_env%para_env%mepos,&
                       final_evaluation=.FALSE., para_env=gopt_env%force_env%para_env)

       curvature1 = dimer_env%rot%curvature
       a1         = (curvature0 - curvature1 + b1 * SIN(2.0_dp*angle))/(1.0_dp - COS(2.0_dp*angle))
       a0         = 2.0_dp*(curvature0-a1)
       angle      = 0.5_dp*ATAN(b1/a1)
       curvature2 = a0/2.0_dp + a1 * COS(2.0_dp*angle)+b1 * SIN(2.0_dp*angle)
       IF (curvature2>curvature0) THEN
          angle = angle + pi/2.0_dp
          curvature2 = a0/2.0_dp + a1 * COS(2.0_dp*angle)+b1 * SIN(2.0_dp*angle)
       END IF
       dimer_env%rot%angle2    = angle
       dimer_env%rot%curvature = curvature2
       ! Rotating the dimer the optimized (in plane) vector position
       dimer_env%nvec = dimer_env%cg_rot%nvec_old
       CALL rotate_dimer(dimer_env%nvec,theta,angle)

       ! Evaluate (by interpolation) the norm of the rotational force in the
       ! minimum of the rotational search (this is for print-out only)
       ALLOCATE(work(SIZE(dimer_env%nvec)),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
       work = dimer_env%rot%g1
       work = SIN(dimer_env%rot%angle1-dimer_env%rot%angle2)/SIN(dimer_env%rot%angle1)*dimer_env%rot%g1+&
              SIN(dimer_env%rot%angle2)/SIN(dimer_env%rot%angle1)*dimer_env%rot%g1p+&
              (1.0_dp-COS(dimer_env%rot%angle2)-SIN(dimer_env%rot%angle2)*TAN(dimer_env%rot%angle1/2.0_dp))*&
              dimer_env%rot%g0
       work = -2.0_dp*(work-dimer_env%rot%g0)
       work = work  - DOT_PRODUCT(work, dimer_env%nvec) * dimer_env%nvec
       opt_energy = SQRT(DOT_PRODUCT(work, work))
       DEALLOCATE(work,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    END IF
    dimer_env%rot%angle2 = angle
    CALL timestop(handle)

  END SUBROUTINE rotmin_2pnt

! *****************************************************************************
!> \brief Line Minimization - Fit
!> \param gopt_env ...
!> \param xvec ...
!> \param xi ...
!> \param opt_energy ...
!> \param brack_limit ...
!> \param step ...
!> \param output_unit ...
!> \param gopt_param ...
!> \param globenv ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
!> \note
!>      Given as input the vector XVEC and XI, finds the scalar
!>      xmin that minimizes the energy XVEC+xmin*XI. Replace step
!>      with the optimal value. Enhanced Version
! *****************************************************************************
  SUBROUTINE linmin_fit(gopt_env, xvec, xi, opt_energy, &
       brack_limit, step, output_unit, gopt_param, globenv)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, brack_limit, step
    INTEGER                                  :: output_unit
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'linmin_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, loc_iter, odim, stat
    LOGICAL                                  :: failure, should_stop
    REAL(KIND=dp)                            :: ax, bx, fprev, rms_dr, &
                                                rms_force, scale, xmin, xx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hist

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(pcom,xicom,hist)
    rms_dr    = gopt_param%rms_dr
    rms_force = gopt_param%rms_force
    ALLOCATE(pcom(SIZE(xvec)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    ALLOCATE(xicom(SIZE(xvec)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)

    pcom  = xvec
    xicom = xi
    xicom = xicom / SQRT(DOT_PRODUCT(xicom,xicom))
    step = step * 0.8_dp ! target a little before the minimum for the first point
    ax = 0.0_dp
    xx = step
    CALL cg_mnbrak(gopt_env, ax, xx, bx, pcom, xicom, brack_limit, output_unit, &
         histpoint=hist,globenv=globenv)
    !
    fprev      = 0.0_dp
    opt_energy = MINVAL(hist(:,2))
    odim =  SIZE(hist,1)
    scale = 0.25_dp
    loc_iter = 0
    DO WHILE ( ABS(hist(odim,3)) > rms_force*scale .OR. ABS(hist(odim,1)-hist(odim-1,1)) > scale*rms_dr)
       CALL external_control(should_stop,"LINFIT",globenv=globenv)
       IF(should_stop) EXIT
       !
       loc_iter = loc_iter + 1
       fprev = opt_energy
       xmin  = FindMin(hist(:,1),hist(:,2),hist(:,3))
       CALL reallocate(hist,1,odim+1,1,3)
       hist(odim+1,1) = xmin
       hist(odim+1,3) = cg_deval1d(gopt_env,xmin,pcom,xicom,opt_energy)
       hist(odim+1,2) = opt_energy
       odim = SIZE(hist,1)
    END DO
    !
    xicom = xmin*xicom
    step = xmin
    xvec  = xvec + xicom
    DEALLOCATE(pcom,  stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(xicom, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(hist, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
            "***","FIT LS  - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
       WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
    END IF
    CALL timestop(handle)

  END SUBROUTINE linmin_fit

! *****************************************************************************
!> \brief Line Minimization routine - GOLD
!> \param gopt_env ...
!> \param xvec ...
!> \param xi ...
!> \param opt_energy ...
!> \param brent_tol ...
!> \param brent_max_iter ...
!> \param brack_limit ...
!> \param step ...
!> \param output_unit ...
!> \param globenv ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
!> \note
!>      Given as input the vector XVEC and XI, finds the scalar
!>      xmin that minimizes the energy XVEC+xmin*XI. Replaces XMIN
!>      with the optimal value
! *****************************************************************************
  SUBROUTINE linmin_gold(gopt_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, &
       brack_limit, step, output_unit, globenv)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, brent_tol
    INTEGER                                  :: brent_max_iter
    REAL(KIND=dp)                            :: brack_limit, step
    INTEGER                                  :: output_unit
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'linmin_gold', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ax, bx, xmin, xx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(pcom,xicom)
    ALLOCATE(pcom(SIZE(xvec)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    ALLOCATE(xicom(SIZE(xvec)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)

    pcom  = xvec
    xicom = xi
    xicom = xicom / SQRT(DOT_PRODUCT(xicom,xicom))
    step = step * 0.8_dp ! target a little before the minimum for the first point
    ax = 0.0_dp
    xx = step
    CALL cg_mnbrak(gopt_env, ax, xx, bx, pcom, xicom, brack_limit, output_unit,&
         globenv=globenv)

    opt_energy = cg_dbrent(gopt_env, ax, xx, bx, brent_tol, brent_max_iter, &
                           xmin, pcom, xicom, output_unit, globenv)
    xicom = xmin*xicom
    step = xmin
    xvec  = xvec + xicom
    DEALLOCATE(pcom,  stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(xicom, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    CALL timestop(handle)
  END SUBROUTINE linmin_gold

! *****************************************************************************
!> \brief Routine for intially bracketing a minimum based on the golden search
!>      minimum
!> \param gopt_env ...
!> \param ax ...
!> \param bx ...
!> \param cx ...
!> \param pcom ...
!> \param xicom ...
!> \param brack_limit ...
!> \param output_unit ...
!> \param histpoint ...
!> \param globenv ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
!> \note
!>      Given two distinct initial points ax and bx this routine searches
!>      in the downhill direction and returns new points ax, bx, cx that
!>      bracket the minimum of the function
! *****************************************************************************
  SUBROUTINE cg_mnbrak(gopt_env, ax, bx, cx, pcom, xicom, brack_limit, output_unit,&
       histpoint, globenv)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: ax, bx, cx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: brack_limit
    INTEGER                                  :: output_unit
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: histpoint
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_mnbrak', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, loc_iter, odim, stat
    LOGICAL                                  :: failure, hist, should_stop
    REAL(KIND=dp)                            :: dum, fa, fb, fc, fu, gold, q, &
                                                r, u, ulim

    CALL timeset(routineN,handle)
    failure = .FALSE.
    hist    = PRESENT(histpoint)
    IF (hist) THEN
       CPPrecondition(.NOT.ASSOCIATED(histpoint),cp_failure_level,routineP,failure)
       ALLOCATE(histpoint(3,3), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    END IF
    gold = (1.0_dp+SQRT(5.0_dp))/2.0_dp
    IF (hist) THEN
       histpoint(1,1) = ax
       histpoint(1,3) = cg_deval1d(gopt_env,ax,pcom,xicom,fa)
       histpoint(1,2) = fa
       histpoint(2,1) = bx
       histpoint(2,3) = cg_deval1d(gopt_env,bx,pcom,xicom,fb)
       histpoint(2,2) = fb
    ELSE
       fa=cg_eval1d(gopt_env,ax,pcom,xicom)
       fb=cg_eval1d(gopt_env,bx,pcom,xicom)
    END IF
    IF(fb.GT.fa)THEN
       dum=ax
       ax=bx
       bx=dum
       dum=fb
       fb=fa
       fa=dum
    ENDIF
    cx=bx+gold*(bx-ax)
    IF (hist) THEN
       histpoint(3,1) = cx
       histpoint(3,3) = cg_deval1d(gopt_env,cx,pcom,xicom,fc)
       histpoint(3,2) = fc
    ELSE
       fc=cg_eval1d(gopt_env,cx,pcom,xicom)
    END IF
    loc_iter = 3
    DO WHILE (fb.GE.fc)
       CALL external_control(should_stop,"MNBRACK",globenv=globenv)
       IF(should_stop) EXIT
       !
       r=(bx-ax)*(fb-fc)
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*SIGN(MAX(ABS(q-r),TINY(0.0_dp)),q-r))
       ulim=bx+brack_limit*(cx-bx)
       IF ((bx-u)*(u-cx).GT.0.0_dp) THEN
          IF (hist) THEN
             odim = SIZE(histpoint,1)
             CALL reallocate(histpoint,1,odim+1,1,3)
             histpoint(odim+1,1) = u
             histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu)
             histpoint(odim+1,2) = fu
          ELSE
             fu=cg_eval1d(gopt_env,u,pcom,xicom)
          END IF
          loc_iter = loc_iter + 1
          IF(fu.LT.fc)THEN
             ax=bx
             fa=fb
             bx=u
             fb=fu
             EXIT
          ELSE IF(fu.GT.fb)THEN
             cx=u
             fc=fu
             EXIT
          ENDIF
          u=cx+gold*(cx-bx)
          IF (hist) THEN
             odim = SIZE(histpoint,1)
             CALL reallocate(histpoint,1,odim+1,1,3)
             histpoint(odim+1,1) = u
             histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu)
             histpoint(odim+1,2) = fu
          ELSE
             fu=cg_eval1d(gopt_env,u,pcom,xicom)
          END IF
          loc_iter = loc_iter + 1
       ELSE IF((cx-u)*(u-ulim).GT.0.)THEN
          IF (hist) THEN
             odim = SIZE(histpoint,1)
             CALL reallocate(histpoint,1,odim+1,1,3)
             histpoint(odim+1,1) = u
             histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu)
             histpoint(odim+1,2) = fu
          ELSE
             fu=cg_eval1d(gopt_env,u,pcom,xicom)
          END IF
          loc_iter = loc_iter + 1
          IF(fu.LT.fc)THEN
             bx=cx
             cx=u
             u=cx+gold*(cx-bx)
             fb=fc
             fc=fu
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(gopt_env,u,pcom,xicom)
             END IF
             loc_iter = loc_iter + 1
          ENDIF
       ELSE IF((u-ulim)*(ulim-cx).GE.0.)THEN
          u=ulim
          IF (hist) THEN
             odim = SIZE(histpoint,1)
             CALL reallocate(histpoint,1,odim+1,1,3)
             histpoint(odim+1,1) = u
             histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu)
             histpoint(odim+1,2) = fu
          ELSE
             fu=cg_eval1d(gopt_env,u,pcom,xicom)
          END IF
          loc_iter = loc_iter + 1
       ELSE
          u=cx+gold*(cx-bx)
          IF (hist) THEN
             odim = SIZE(histpoint,1)
             CALL reallocate(histpoint,1,odim+1,1,3)
             histpoint(odim+1,1) = u
             histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu)
             histpoint(odim+1,2) = fu
          ELSE
             fu=cg_eval1d(gopt_env,u,pcom,xicom)
          END IF
          loc_iter = loc_iter + 1
       ENDIF
       ax=bx
       bx=cx
       cx=u
       fa=fb
       fb=fc
       fc=fu
    END DO
    IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
            "***","MNBRACK - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
       WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
    END IF
    CALL timestop(handle)
  END SUBROUTINE cg_mnbrak

! *****************************************************************************
!> \brief Routine implementing the Brent Method
!>      Brent,R.P. Algorithm for Minimization without Derivatives, Chapt.5
!>      1973
!>      Extension in the use of derivatives
!> \param gopt_env ...
!> \param ax ...
!> \param bx ...
!> \param cx ...
!> \param tol ...
!> \param itmax ...
!> \param xmin ...
!> \param pcom ...
!> \param xicom ...
!> \param output_unit ...
!> \param globenv ...
!> \retval dbrent ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
!> \note
!>      Given a bracketing  triplet of abscissas ax, bx, cx (such that bx
!>      is between ax and cx and energy of bx is less than energy of ax and cx),
!>      this routine isolates the minimum to a precision of about tol using
!>      Brent method. This routine implements the extension of the Brent Method
!>      using derivatives
! *****************************************************************************
  FUNCTION cg_dbrent(gopt_env,ax,bx,cx,tol,itmax,xmin,pcom,xicom,output_unit,&
       globenv) RESULT(dbrent)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: ax, bx, cx, tol
    INTEGER                                  :: itmax
    REAL(KIND=dp)                            :: xmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    INTEGER                                  :: output_unit
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp)                            :: dbrent

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_dbrent', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: zeps = 1.0E-8_dp

    INTEGER                                  :: handle, iter, loc_iter
    LOGICAL                                  :: failure, ok1, ok2, &
                                                should_stop, skip0, skip1
    REAL(KIND=dp)                            :: a, b, d, d1, d2, du, dv, dw, &
                                                dx, e, fu, fv, fw, fx, olde, &
                                                tol1, tol2, u, u1, u2, v, w, &
                                                x, xm

    CALL timeset(routineN,handle)
    failure = .FALSE.
    a  = MIN(ax,cx)
    b  = MAX(ax,cx)
    v  = bx ; w  = v; x  = v
    e  = 0.0_dp
    dx = cg_deval1d(gopt_env,x,pcom,xicom,fx)
    fv = fx
    fw = fx
    dv = dx
    dw = dx
    loc_iter = 1
    DO  iter = 1, itmax
       CALL external_control(should_stop,"BRENT",globenv=globenv)
       IF(should_stop) EXIT
       !
       xm    = 0.5_dp * (a+b)
       tol1  = tol*ABS(x)+zeps
       tol2  = 2.0_dp * tol1
       skip0 = .FALSE.
       skip1 = .FALSE.
       IF(ABS(x-xm).LE.(tol2-0.5_dp*(b-a))) EXIT
       IF(ABS(e).GT.tol1) THEN
          d1 = 2.0_dp*(b-a)
          d2 = d1
          IF(dw.NE.dx) d1=(w-x)*dx/(dx-dw)
          IF(dv.NE.dx) d2=(v-x)*dx/(dx-dv)
          u1   = x+d1
          u2   = x+d2
          ok1  = ((a-u1)*(u1-b).GT.0.0_dp).AND.(dx*d1.LE.0.0_dp)
          ok2  = ((a-u2)*(u2-b).GT.0.0_dp).AND.(dx*d2.LE.0.0_dp)
          olde = e
          e    = d
          IF(.NOT.(ok1.OR.ok2)) THEN
             skip0 = .TRUE.
          ELSE IF (ok1.AND.ok2) THEN
             IF(ABS(d1).LT.ABS(d2)) THEN
                d = d1
             ELSE
                d = d2
             ENDIF
          ELSE IF (ok1) THEN
             d = d1
          ELSE
             d = d2
          ENDIF
          IF (.NOT.skip0) THEN
             IF(ABS(d).GT.ABS(0.5_dp*olde)) skip0=.TRUE.
             IF (.NOT.skip0) THEN
                u=x+d
                IF((u-a).LT.tol2 .OR. (b-u).LT.tol2) d=SIGN(tol1,xm-x)
                skip1=.TRUE.
             END IF
          END IF
       ENDIF
       IF (.NOT.skip1) THEN
          IF(dx.GE.0.0_dp) THEN
             e=a-x
          ELSE
             e=b-x
          ENDIF
          d=0.5_dp*e
       END IF
       IF(ABS(d).GE.tol1) THEN
          u=x+d
          du=cg_deval1d(gopt_env,u,pcom,xicom,fu)
          loc_iter = loc_iter + 1
       ELSE
          u=x+SIGN(tol1,d)
          du=cg_deval1d(gopt_env,u,pcom,xicom,fu)
          loc_iter = loc_iter + 1
          IF(fu.GT.fx) EXIT
       ENDIF
       IF(fu.LE.fx) THEN
          IF(u.GE.x) THEN
             a=x
          ELSE
             b=x
          ENDIF
          v=w;   fv=fw;  dv=dw;  w=x
          fw=fx; dw=dx;  x=u;    fx=fu;  dx=du
       ELSE
          IF(u.LT.x) THEN
             a=u
          ELSE
             b=u
          ENDIF
          IF(fu.LE.fw .OR. w.EQ.x) THEN
             v=w; fv=fw; dv=dw
             w=u; fw=fu; dw=du
          ELSE IF(fu.LE.fv .OR. v.EQ.x .OR. v.EQ.w) THEN
             v=u
             fv=fu
             dv=du
          ENDIF
       ENDIF
    END DO
    IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
       WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
            "***","BRENT   - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
       IF (iter==itmax+1) &
            WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,T78,A)")&
            "***","BRENT - NUMBER OF ITERATIONS EXCEEDED ","***"
       WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
    END IF
    CPPrecondition((iter/=itmax+1),cp_failure_level,routineP,failure)
    xmin=x
    dbrent=fx
    CALL timestop(handle)

  END FUNCTION cg_dbrent

! *****************************************************************************
!> \brief Evaluates energy in one dimensional space defined by the point
!>      pcom and with direction xicom, position x
!> \param gopt_env ...
!> \param x ...
!> \param pcom ...
!> \param xicom ...
!> \retval my_val ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  FUNCTION cg_eval1d(gopt_env,x,pcom,xicom) RESULT(my_val)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: my_val

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_eval1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ALLOCATE(xvec(SIZE(pcom)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    xvec = pcom + x*xicom
    CALL cp_eval_at(gopt_env, xvec, my_val, master=gopt_env%force_env%para_env%mepos,&
                    final_evaluation=.FALSE., para_env=gopt_env%force_env%para_env)
    DEALLOCATE(xvec, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END FUNCTION cg_eval1d

! *****************************************************************************
!> \brief Evaluates derivatives in one dimensional space defined by the point
!>      pcom and with direction xicom, position x
!> \param gopt_env ...
!> \param x ...
!> \param pcom ...
!> \param xicom ...
!> \param fval ...
!> \retval my_val ...
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  FUNCTION cg_deval1d(gopt_env,x,pcom,xicom,fval) RESULT(my_val)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: fval, my_val

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_deval1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: grad, xvec

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ALLOCATE(xvec(SIZE(pcom)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    ALLOCATE(grad(SIZE(pcom)), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    xvec = pcom + x * xicom
    CALL cp_eval_at(gopt_env, xvec, energy, grad, master=gopt_env%force_env%para_env%mepos,&
                    final_evaluation=.FALSE., para_env=gopt_env%force_env%para_env)
    my_val = DOT_PRODUCT(grad,xicom)
    fval   = energy
    DEALLOCATE(xvec, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(grad, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    CALL timestop(handle)

  END FUNCTION cg_deval1d

! *****************************************************************************
!> \brief Find the minimum of a parabolic function obtained with a least square fit
!> \param x ...
!> \param y ...
!> \param dy ...
!> \retval res ...
!> \par History
!>      10.2005 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION FindMin(x,y,dy) RESULT(res)
    REAL(kind=dp), DIMENSION(:)              :: x, y, dy
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'FindMin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, iwork(8*3), lwork, &
                                                min_pos, np, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: diag(3), res1(3), res2(3), &
                                                res3(3), spread, sum_x, &
                                                sum_xx, tmpw(1), vt(3,3)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: work
    REAL(kind=dp), DIMENSION(2*SIZE(x), 3)   :: f
    REAL(kind=dp), DIMENSION(2*SIZE(x))      :: b, w
    REAL(kind=dp)                            :: u(2*SIZE(x),3)

    failure=.FALSE.
    np=SIZE(x)
    CPPrecondition(np>1,cp_failure_level,routineP,failure)
    sum_x=0._dp
    sum_xx=0._dp
    min_pos=1
    DO i=1,np
       sum_xx=sum_xx+x(i)**2
       sum_x=sum_x+x(i)
       IF (y(min_pos)>y(i)) min_pos=i
    END DO
    spread=SQRT(sum_xx/REAL(np,dp)-(sum_x/REAL(np,dp))**2)
    DO i=1,np
       w(i)=EXP(-(REAL(np-i,dp))**2/(REAL(2*9,dp)))
       w(i+np)=2._dp*w(i)
    END DO
    DO i=1,np
       f(i,1)=w(i)
       f(i,2)=x(i)*w(i)
       f(i,3)=x(i)**2*w(i)
       f(i+np,1)=0
       f(i+np,2)=w(i+np)
       f(i+np,3)=2*x(i)*w(i+np)
    END DO
    DO i=1,np
       b(i)=y(i)*w(i)
       b(i+np)=dy(i)*w(i+np)
    END DO
    lwork=-1
    CALL dgesdd('S',SIZE(f,1),SIZE(f,2),f,SIZE(f,1),diag,u,SIZE(u,1),vt,SIZE(vt,1),tmpw,lwork,&
         iwork,info)
    lwork=CEILING(tmpw(1))
    ALLOCATE(work(lwork),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    CALL dgesdd('S',SIZE(f,1),SIZE(f,2),f,SIZE(f,1),diag,u,SIZE(u,1),vt,SIZE(vt,1),work,lwork,&
         iwork,info)
    DEALLOCATE(work, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    CALL dgemv('T',SIZE(u,1),SIZE(u,2),1._dp,u,SIZE(u,1),b,1,0._dp,res1,1)
    DO i=1,3
       res2(i)=res1(i)/diag(i)
    END DO
    CALL dgemv('T',3,3,1._dp,vt,SIZE(vt,1),res2,1,0._dp,res3,1)
    res=-0.5*res3(2)/res3(3)
  END FUNCTION FindMin

! *****************************************************************************
!> \brief Computes the Conjugate direction for the next search
!> \param gopt_env ...
!> \param Fletcher_Reeves ...
!> \param g contains the theta  of the previous step.. (norm 1.0 vector)
!> \param xi contains the -theta of the present step.. (norm 1.0 vector)
!> \param h contains the search direction of the previous step (must be orthogonal
!>            to nvec of the previous step (nvec_old))
!> \par   Info for DIMER method
!> \par History
!>      10.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE get_conjugate_direction(gopt_env, Fletcher_Reeves, g, xi, h)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    LOGICAL, INTENT(IN)                      :: Fletcher_Reeves
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g, xi, h

    CHARACTER(len=*), PARAMETER :: routineN = 'get_conjugate_direction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: check, failure
    REAL(KIND=dp)                            :: dgg, gam, gg, norm, norm_h
    TYPE(dimer_env_type), POINTER            :: dimer_env

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(dimer_env)
    IF (.NOT.gopt_env%dimer_rotation) THEN
       gg  = DOT_PRODUCT(g,g)
       IF (Fletcher_Reeves) THEN
          dgg = DOT_PRODUCT(xi,xi)
       ELSE
          dgg = DOT_PRODUCT((xi+g),xi)
       END IF
       gam = dgg/gg
       g   =   h
       h   = -xi + gam * h
    ELSE
       dimer_env => gopt_env%dimer_env
       check = ABS(DOT_PRODUCT(g,g)-1.0_dp)<MAX(1.0E-9_dp,dimer_thrs)
       CPPrecondition(check,cp_failure_level,routineP,failure)

       check = ABS(DOT_PRODUCT(xi,xi)-1.0_dp)<MAX(1.0E-9_dp,dimer_thrs)
       CPPrecondition(check,cp_failure_level,routineP,failure)

       check = ABS(DOT_PRODUCT(h,dimer_env%cg_rot%nvec_old))<MAX(1.0E-9_dp,dimer_thrs)
       CPPrecondition(check,cp_failure_level,routineP,failure)
       gg  = dimer_env%cg_rot%norm_theta_old**2
       IF (Fletcher_Reeves) THEN
          dgg = dimer_env%cg_rot%norm_theta**2
       ELSE
          norm= dimer_env%cg_rot%norm_theta*dimer_env%cg_rot%norm_theta_old
          dgg = dimer_env%cg_rot%norm_theta**2 + DOT_PRODUCT(g,xi)*norm
       END IF
       ! Compute Theta** and store it in nvec_old
       CALL rotate_dimer(dimer_env%cg_rot%nvec_old,g,dimer_env%rot%angle2+pi/2.0_dp)
       gam = dgg/gg
       g   =    h
       h   = - xi * dimer_env%cg_rot%norm_theta + gam * dimer_env%cg_rot%norm_h * dimer_env%cg_rot%nvec_old
       h   =    h - DOT_PRODUCT(h,dimer_env%nvec)*dimer_env%nvec
       norm_h = SQRT(DOT_PRODUCT(h,h))
       IF (norm_h < EPSILON(0.0_dp)) THEN
          h = 0.0_dp
       ELSE
          h = h/norm_h
       END IF
       dimer_env%cg_rot%norm_h = norm_h
    END IF
    CALL timestop(handle)

  END SUBROUTINE get_conjugate_direction

END MODULE cg_utils
