!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/cg_utils [1.0] *
!!
!!   NAME
!!     cg_utils
!!
!!   FUNCTION
!!     Utilities for Geometry optimization using  Conjugate Gradients
!!
!!   AUTHOR
!!     Teodoro Laino [teo]
!!     10.2005
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************

MODULE cg_utils
  USE dimer_types,                     ONLY: dimer_env_type
  USE dimer_utils,                     ONLY: rotate_dimer
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_minimization_method_id,&
                                             default_ts_method_id,&
                                             do_rotate_dimer,&
                                             do_translate_dimer
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  PUBLIC :: cg_linmin, get_conjugate_direction
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cg_utils'

CONTAINS

!!****f* cg_utils/cg_linmin [1.0] *
!!
!!   NAME
!!     cg_linmin
!!
!!   FUNCTION
!!     Main driver for line minimization routines for CG
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE cg_linmin(gopt_env, xvec, xi, opt_energy, step, output_unit,&
       gopt_param, globenv, error)

    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, step
    INTEGER                                  :: output_unit
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_linmin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: brent_max_iter, handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: brack_limit, brent_tol

    failure = .FALSE.
    CALL timeset(routineN,handle)
    SELECT CASE(gopt_env%type_id)
    CASE(default_minimization_method_id)
       IF (gopt_param%use_fit_line_search) THEN
          CALL linminE(gopt_env, xvec, xi, opt_energy, gopt_param%brent_tol,&
               gopt_param%brent_max_iter, gopt_param%brack_limit, &
               step, output_unit, gopt_param, globenv, error)
       ELSE
          CALL linmin(gopt_env, xvec, xi, opt_energy, gopt_param%brent_tol,&
               gopt_param%brent_max_iter, gopt_param%brack_limit, &
               step, output_unit, globenv, error)
       END IF
    CASE(default_ts_method_id)
       SELECT CASE(gopt_env%dimer_env%status_dimer)
       CASE(do_rotate_dimer)
          CALL rotmin(gopt_env, gopt_env%dimer_env, xvec, xi, error)
       CASE(do_translate_dimer)
          ! To Fix
          STOP
       CASE DEFAULT
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END SELECT
    CALL timestop(handle)

  END SUBROUTINE cg_linmin

!!****f* cg_utils/rotmin [1.0] *
!!
!!   NAME
!!     rotmin
!!
!!   FUNCTION
!!     Rotational minimization for the Dimer MEthod
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE rotmin(gopt_env, dimer_env, x0, theta, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(dimer_env_type), POINTER            :: dimer_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0, theta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotmin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a0, a1, angle, b1, &
                                                curvature0, curvature1, &
                                                curvature2, dCdp, f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: nvec

    failure = .FALSE.
    CALL timeset(routineN,handle)
    ALLOCATE(nvec(SIZE(dimer_env%nvec)),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    nvec       =  dimer_env%nvec
    curvature0 =  dimer_env%curvature
    dCdp       =  dimer_env%dCdp
    b1         =  0.5_dp*dCdp 
    angle      = -0.5_dp*ATAN(dCdp/(2.0_dp*ABS(curvature0)))
    dimer_env%angle1 = angle
    IF (angle>dimer_env%angle_tol) THEN
       ! Rotating the dimer of dtheta degrees
       CALL rotate_dimer(dimer_env%nvec,theta,angle,error=error)
       ! Re-compute energy, gradients and rotation vector for new R1
       CALL cp_eval_at(gopt_env, x0, f, master=gopt_env%force_env%para_env%mepos,&
            para_env=gopt_env%force_env%para_env, error=error)

       curvature1 = dimer_env%curvature
       a1         = (curvature0 - curvature1 + b1 * SIN(2.0_dp*angle))/(1.0_dp - COS(2.0_dp*angle))
       a0         = 2.0_dp*(curvature0-a1)
       angle      = 0.5_dp*ATAN(b1/a1)
       curvature2 = a0/2.0_dp + a1 * COS(2.0_dp*angle)+b1 * SIN(2.0_dp*angle)
       IF (curvature2>curvature0) THEN 
          angle = angle + pi/2.0_dp
          curvature2 = a0/2.0_dp + a1 * COS(2.0_dp*angle)+b1 * SIN(2.0_dp*angle)
       END IF
       dimer_env%angle2    = angle
       dimer_env%curvature = curvature2
       ! Rotating the dimer the optimized (in plane) vector position
       dimer_env%nvec = nvec
       CALL rotate_dimer(dimer_env%nvec,theta,angle,error=error)
    END IF
    DEALLOCATE(nvec,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    dimer_env%angle2 = angle
    CALL timestop(handle)

  END SUBROUTINE rotmin

!!****f* cg_utils/LinMinE [1.0] *
!!
!!   NAME
!!     LinMinE
!!
!!   FUNCTION
!!     Line Minimization routine Enhanced
!!
!!   NOTES
!!     Given as input the vector XVEC and XI, finds the scalar
!!     xmin that minimizes the energy XVEC+xmin*XI. Replace step
!!     with the optimal value. Enhanced Version
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE LinMinE(gopt_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, &
       brack_limit, step, output_unit, gopt_param, globenv, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, brent_tol
    INTEGER                                  :: brent_max_iter
    REAL(KIND=dp)                            :: brack_limit, step
    INTEGER                                  :: output_unit
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'LinMinE', &
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
    IF (.NOT.failure) THEN
       ALLOCATE(pcom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xicom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       pcom  = xvec
       xicom = xi
       xicom = xicom / SQRT(DOT_PRODUCT(xicom,xicom))
       step = step * 0.8_dp ! target a little before the minimum for the first point
       ax = 0.0_dp
       xx = step
       CALL cg_mnbrak(gopt_env, ax, xx, bx, pcom, xicom, brack_limit, output_unit, &
            histpoint=hist,globenv=globenv, error=error)
       !
       fprev      = 0.0_dp
       opt_energy = MINVAL(hist(:,2))
       odim =  SIZE(hist,1)
       scale = 0.25_dp
       loc_iter = 0
       DO WHILE ( ABS(hist(odim,3)) > rms_force*scale .OR. ABS(hist(odim,1)-hist(odim-1,1)) > scale*rms_dr)
          CALL external_control(should_stop,"LINFIT",globenv,error)
          IF(should_stop) EXIT
          !
          loc_iter = loc_iter + 1
          fprev = opt_energy
          xmin  = FindMin(hist(:,1),hist(:,2),hist(:,3),error)
          CALL reallocate(hist,1,odim+1,1,3)
          hist(odim+1,1) = xmin
          hist(odim+1,3) = cg_deval1d(gopt_env,xmin,pcom,xicom,opt_energy,error)
          hist(odim+1,2) = opt_energy
          odim = SIZE(hist,1)
       END DO
       !
       xicom = xmin*xicom
       step = xmin
       xvec  = xvec + xicom
       DEALLOCATE(pcom,  stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(xicom, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(hist, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (output_unit>0) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
          WRITE(UNIT=output_unit,FMT="(T2,A,T22,A,I7,T78,A)")&
               "***","FIT LS  - NUMBER OF ENERGY EVALUATIONS : ",loc_iter,"***"
          WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
       END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE LinMinE

!!****f* cg_utils/LinMin [1.0] *
!!
!!   NAME
!!     LinMin
!!
!!   FUNCTION
!!     Line Minimization routine
!!
!!   NOTES
!!     Given as input the vector XVEC and XI, finds the scalar
!!     xmin that minimizes the energy XVEC+xmin*XI. Replaces XMIN
!!     with the optimal value
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE LinMin(gopt_env, xvec, xi, opt_energy, brent_tol, brent_max_iter, &
       brack_limit, step, output_unit, globenv, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec, xi
    REAL(KIND=dp)                            :: opt_energy, brent_tol
    INTEGER                                  :: brent_max_iter
    REAL(KIND=dp)                            :: brack_limit, step
    INTEGER                                  :: output_unit
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'LinMin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ax, bx, xmin, xx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(pcom,xicom)
    IF (.NOT.failure) THEN
       ALLOCATE(pcom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xicom(SIZE(xvec)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       pcom  = xvec
       xicom = xi
       xicom = xicom / SQRT(DOT_PRODUCT(xicom,xicom))
       step = step * 0.8_dp ! target a little before the minimum for the first point
       ax = 0.0_dp
       xx = step
       CALL cg_mnbrak(gopt_env, ax, xx, bx, pcom, xicom, brack_limit, output_unit,&
            globenv=globenv, error=error)
       opt_energy = cg_dbrent(gopt_env, ax, xx, bx, brent_tol, brent_max_iter, &
                              xmin, pcom, xicom, output_unit, globenv, error)
       xicom = xmin*xicom
       step = xmin
       xvec  = xvec + xicom
       DEALLOCATE(pcom,  stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(xicom, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE LinMin

!!****f* cg_utils/cg_mnbrak [1.0] *
!!
!!   NAME
!!     cg_mnbrak
!!
!!   FUNCTION
!!     Routine for intially bracketing a minimum based on the golden search
!!     minimum
!!
!!   NOTES
!!     Given two distinct initial points ax and bx this routine searches
!!     in the downhill direction and returns new points ax, bx, cx that
!!     bracket the minimum of the function
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE cg_mnbrak(gopt_env, ax, bx, cx, pcom, xicom, brack_limit, output_unit,&
       histpoint, globenv, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: ax, bx, cx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: brack_limit
    INTEGER                                  :: output_unit
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: histpoint
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

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
       CPPrecondition(.NOT.ASSOCIATED(histpoint),cp_failure_level,routineP,error,failure)
       ALLOCATE(histpoint(3,3), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    gold = (1.0_dp+SQRT(5.0_dp))/2.0_dp
    IF (.NOT.failure) THEN
       IF (hist) THEN
          histpoint(1,1) = ax
          histpoint(1,3) = cg_deval1d(gopt_env,ax,pcom,xicom,fa,error)
          histpoint(1,2) = fa
          histpoint(2,1) = bx
          histpoint(2,3) = cg_deval1d(gopt_env,bx,pcom,xicom,fb,error)
          histpoint(2,2) = fb
       ELSE
          fa=cg_eval1d(gopt_env,ax,pcom,xicom,error=error)
          fb=cg_eval1d(gopt_env,bx,pcom,xicom,error=error)
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
          histpoint(3,3) = cg_deval1d(gopt_env,cx,pcom,xicom,fc,error)
          histpoint(3,2) = fc
       ELSE
          fc=cg_eval1d(gopt_env,cx,pcom,xicom,error=error)
       END IF
       loc_iter = 3
       DO WHILE (fb.GE.fc)
          CALL external_control(should_stop,"MNBRACK",globenv,error)
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
                histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(gopt_env,u,pcom,xicom,error=error)
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
                histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(gopt_env,u,pcom,xicom,error=error)
             END IF
             loc_iter = loc_iter + 1
          ELSE IF((cx-u)*(u-ulim).GT.0.)THEN
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(gopt_env,u,pcom,xicom,error=error)
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
                   histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
                   histpoint(odim+1,2) = fu
                ELSE
                   fu=cg_eval1d(gopt_env,u,pcom,xicom,error=error)
                END IF
                loc_iter = loc_iter + 1
             ENDIF
          ELSE IF((u-ulim)*(ulim-cx).GE.0.)THEN
             u=ulim
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(gopt_env,u,pcom,xicom,error=error)
             END IF
             loc_iter = loc_iter + 1
          ELSE
             u=cx+gold*(cx-bx)
             IF (hist) THEN
                odim = SIZE(histpoint,1)
                CALL reallocate(histpoint,1,odim+1,1,3)
                histpoint(odim+1,1) = u
                histpoint(odim+1,3) = cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
                histpoint(odim+1,2) = fu
             ELSE
                fu=cg_eval1d(gopt_env,u,pcom,xicom,error=error)
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
    END IF
    CALL timestop(handle)
  END SUBROUTINE cg_mnbrak

!!****f* cg_utils/cg_dbrent [1.0] *
!!
!!   NAME
!!     cg_dbrent
!!
!!   FUNCTION
!!     Routine implementing the Brent Method
!!     Brent,R.P. Algorithm for Minimization without Derivatives, Chapt.5
!!     1973
!!     Extension in the use of derivatives
!!
!!   NOTES
!!     Given a bracketing  triplet of abscissas ax, bx, cx (such that bx
!!     is between ax and cx and energy of bx is less than energy of ax and cx),
!!     this routine isolates the minimum to a precision of about tol using
!!     Brent method. This routine implements the extension of the Brent Method
!!     using derivatives
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  FUNCTION cg_dbrent(gopt_env,ax,bx,cx,tol,itmax,xmin,pcom,xicom,output_unit,&
       globenv, error) RESULT(dbrent)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: ax, bx, cx, tol
    INTEGER                                  :: itmax
    REAL(KIND=dp)                            :: xmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    INTEGER                                  :: output_unit
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error
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
    IF (.NOT.failure) THEN
       a  = MIN(ax,cx)
       b  = MAX(ax,cx)
       v  = bx ; w  = v; x  = v
       e  = 0.0_dp
       dx = cg_deval1d(gopt_env,x,pcom,xicom,fx,error)
       fv = fx
       fw = fx
       dv = dx
       dw = dx
       loc_iter = 1
       DO  iter = 1, itmax
          CALL external_control(should_stop,"BRENT",globenv,error)
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
             du=cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
             loc_iter = loc_iter + 1
          ELSE
             u=x+SIGN(tol1,d)
             du=cg_deval1d(gopt_env,u,pcom,xicom,fu,error)
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
       CPPrecondition((iter/=itmax+1),cp_failure_level,routineP,error,failure)
       xmin=x
       dbrent=fx
    END IF
    CALL timestop(handle)

  END FUNCTION cg_dbrent

!!****f* cg_utils/cg_eval1d [1.0] *
!!
!!   NAME
!!     cg_eval1d
!!
!!   FUNCTION
!!     Evaluates energy in one dimensional space defined by the point
!!     pcom and with direction xicom, position x
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  FUNCTION cg_eval1d(gopt_env,x,pcom,xicom,error) RESULT(my_val)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: my_val

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_eval1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xvec

    CALL timeset(routineN,handle)

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ALLOCATE(xvec(SIZE(pcom)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       xvec = pcom + x*xicom
       CALL cp_eval_at(gopt_env, xvec, my_val, master=gopt_env%force_env%para_env%mepos,&
            para_env=gopt_env%force_env%para_env, error=error)
       DEALLOCATE(xvec, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END FUNCTION cg_eval1d

!!****f* cg_utils/cg_deval1d [1.0] *
!!
!!   NAME
!!     cg_deval1d
!!
!!   FUNCTION
!!     Evaluates derivatives in one dimensional space defined by the point
!!     pcom and with direction xicom, position x
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  FUNCTION cg_deval1d(gopt_env,x,pcom,xicom,fval,error) RESULT(my_val)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pcom, xicom
    REAL(KIND=dp)                            :: fval
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: my_val

    CHARACTER(len=*), PARAMETER :: routineN = 'cg_deval1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: grad, xvec

    CALL timeset(routineN,handle)

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ALLOCATE(xvec(SIZE(pcom)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(grad(SIZE(pcom)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       xvec = pcom + x * xicom
       CALL cp_eval_at(gopt_env, xvec, energy, grad, master=gopt_env%force_env%para_env%mepos,&
            para_env=gopt_env%force_env%para_env, error=error)
       my_val = DOT_PRODUCT(grad,xicom)
       fval   = energy
       DEALLOCATE(xvec, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(grad, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END FUNCTION cg_deval1d


!!****f* cg_utils/FindFit [1.0] *
!!
!!   NAME
!!     FindMin
!!
!!   FUNCTION
!!     Find the minimum of a parabolic function obtained with a least square fit
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION FindMin(x,y,dy,error) RESULT(res)
    REAL(kind=dp), DIMENSION(:)              :: x, y, dy
    TYPE(cp_error_type), INTENT(inout)       :: error
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
    CPPrecondition(np>1,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
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
       CALL dgesdd('S',SIZE(f,1),SIZE(f,2),f,SIZE(f,1),diag,u,SIZE(u,1),vt,SIZE(vt,1),work,lwork,&
            iwork,info)
       DEALLOCATE(work, stat=stat)
       CALL dgemv('T',SIZE(u,1),SIZE(u,2),1._dp,u,SIZE(u,1),b,1,0._dp,res1,1)
       DO i=1,3
          res2(i)=res1(i)/diag(i)
       END DO
       CALL dgemv('T',3,3,1._dp,vt,SIZE(vt,1),res2,1,0._dp,res3,1)
       res=-0.5*res3(2)/res3(3)
    END IF
  END FUNCTION FindMin

!!****f* cg_utils/get_conjugate_direction [1.0] *
!!
!!   NAME
!!     get_conjugate_direction
!!
!!   FUNCTION
!!     Computes the Conjugate direction for the next search
!!
!!   NOTES
!!
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE get_conjugate_direction(dimer_rotation, Fletcher_Reeves, g, xi, h)
    LOGICAL, INTENT(IN)                      :: dimer_rotation, &
                                                Fletcher_Reeves
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g, xi, h

    CHARACTER(len=*), PARAMETER :: routineN = 'get_conjugate_direction', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dgg, gam, gg

    IF (.NOT.dimer_rotation) THEN
       gg  = DOT_PRODUCT(g,g)
       IF (Fletcher_Reeves) THEN
          dgg = DOT_PRODUCT(xi,xi)
       ELSE
          dgg = DOT_PRODUCT((xi+g),xi)
       END IF
       gam = dgg/gg
       g   = -xi
       h   = g + gam * h
    ELSE
       ! To Fix
       STOP
    END IF

  END SUBROUTINE get_conjugate_direction

END MODULE cg_utils
