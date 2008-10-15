!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for Geometry optimization using BFGS algorithm
! *****************************************************************************
MODULE bfgs_optimizer

  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_methods,                  ONLY: gopt_f_ii,&
                                             gopt_f_io,&
                                             gopt_f_io_finalize,&
                                             gopt_f_io_init,&
                                             print_geo_opt_header,&
                                             print_geo_opt_nc
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_cell_method_id,&
                                             default_ts_method_id
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bfgs_optimizer'
  LOGICAL, PARAMETER                       :: debug_this_module=.TRUE.

  PUBLIC :: geoopt_bfgs

CONTAINS

! *****************************************************************************
!> \brief Main driver for BFGS geometry optimizations 
! *****************************************************************************
  SUBROUTINE geoopt_bfgs(force_env,gopt_param,globenv,&
       geo_section, gopt_env, x0, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: geo_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geoopt_bfgs', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    CHARACTER(LEN=5)                         :: wildcard
    CHARACTER(LEN=default_path_length)       :: hes_filename
    INTEGER                                  :: handle, handle2, hesunit, &
                                                hesunit_read, i, indf, its, &
                                                j, maxiter, ndf, nfree, &
                                                output_unit, stat
    LOGICAL                                  :: conv, failure, hesrest, &
                                                ionode, should_stop
    REAL(KIND=dp)                            :: ediff, emin, eold, etot, &
                                                pred, rad, rat, step
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dg, dr, dx, eigval, gold, &
                                                hes, work, xold
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: eigvec, hes_mat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section

    failure = .FALSE.
    NULLIFY(logger, g)
    logger       => cp_error_get_logger(error)
    para_env     => force_env%para_env
    root_section => force_env%root_section

    CALL timeset(routineN,handle)
    ionode      = para_env%mepos==para_env%source
    maxiter     = gopt_param%max_iter
    conv        = .FALSE.
    rad         = 0.10_dp
    rat         = 0.0_dp
    wildcard    = " BFGS"

    ! Stop if not implemented
    SELECT CASE(gopt_env%type_id)
    CASE(default_ts_method_id,default_cell_method_id)
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="BFGS method not yet working with for DIMER or CELL_OPT", &
            error=error, error_level=cp_failure_level)
    END SELECT

    CALL section_vals_val_get(geo_section,"BFGS%RESTART_HESSIAN",l_val=hesrest,error=error)
    output_unit = cp_print_key_unit_nr(logger,geo_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".geoLog",error=error)

    ndf   = SIZE(x0)
    nfree = gopt_env%nfree
    
    ! Initialize hessian (hes = unitary matrix )
    ALLOCATE (hes(ndf*(ndf+1)/2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO indf=1,ndf
       DO j=1,indf
          i = indf*(indf - 1)/2 + j
          hes(i) = zero
       END DO
       i = indf*(indf-1)/2 + indf
       hes(i) = one
    END DO

    ALLOCATE (xold(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    xold(:) = x0(:)

    ALLOCATE (g(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    g(:) = 0.0_dp

    ALLOCATE (gold(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    gold(:) = 0.0_dp

    ALLOCATE (dx(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dx(:) = 0.0_dp

    ALLOCATE (hes_mat(ndf,ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    hes_mat(:,:) = 0.0_dp

    ALLOCATE (dg(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dg(:) = 0.0_dp

    ALLOCATE (eigval(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    eigval(:) = 0.0_dp

    ALLOCATE (eigvec(ndf,ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    eigvec(:,:) = 0.0_dp

    ALLOCATE (work(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    work(:) = 0.0_dp

    ALLOCATE (dr(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dr(:) = 0.0_dp

    ! Geometry optimization starts now
    CALL print_geo_opt_header(gopt_env, output_unit, wildcard)

    ! Calculate Energy & Gradients
    CALL cp_eval_at(gopt_env, x0, etot, g, gopt_env%force_env%para_env%mepos,&
         gopt_env%force_env%para_env, error)

    ! Print info at time 0
    emin = etot
    CALL gopt_f_io_init(gopt_env, output_unit, etot, wildcard=wildcard, error=error)
    DO its = 1, maxiter
       gopt_env%its = its
       CALL cp_iterate(logger%iter_info,iter_nr=its,last=(its==maxiter),error=error)
       CALL gopt_f_ii(its, output_unit)

       ! Hessian update/restarting
       IF ((its == 1).AND.hesrest) THEN
          stat = 0
          IF(ionode)THEN
             CALL section_vals_val_get(geo_section,"BFGS%RESTART_FILE_NAME",c_val=hes_filename,&
                  error=error)
             CALL open_file(file_name=hes_filename,file_status="OLD",&
                  file_form="UNFORMATTED", file_action="READ", unit_number=hesunit_read)
             READ(UNIT=hesunit_read,IOSTAT=stat) hes
             CALL close_file(unit_number=hesunit_read)
          END IF
          CALL mp_sum(stat,para_env%group)

          IF (stat /= 0) THEN
             IF (output_unit>0) WRITE(UNIT=output_unit,FMT="(/,T2,A)")"**  Error while reading HESSIAN **"
             conv = .TRUE.
          ELSE
             IF (output_unit>0)  WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                  "*** Initial Hessian has been restarted successfully ***"
          END IF
          CALL mp_bcast(hes,para_env%source,para_env%group)
       ELSE IF( its > 1 ) THEN
          DO indf = 1, ndf
             dx(indf) = x0(indf) - xold(indf)
             dg(indf) =  g(indf) - gold(indf)
          END DO

          CALL bfgs(ndf,dx,dg,hes,work)

          !Possibly dump the Hessian file
          CALL timeset(routineN//"hes2disk",handle2)
          hesunit = cp_print_key_unit_nr(logger,geo_section,"BFGS%RESTART",&
               extension=".Hessian",file_form="UNFORMATTED",file_action="WRITE",&
               file_position="REWIND",error=error)
          IF(hesunit>0)THEN
             WRITE (UNIT=hesunit) hes
          END IF
          CALL cp_print_key_finished_output(hesunit,logger,geo_section,"BFGS%RESTART", error=error)
          CALL timestop(handle2)
       END IF

       ! Setting the present positions & gradients as old
       xold = x0
       gold = g

       ! Copying hessian hes to (ndf x ndf) matrix hes_mat for diagonalization
       DO indf = 1, ndf
          DO j = 1, indf
             i = indf*(indf-1)/2 + j
             eigvec(indf,j) = hes(i)
             eigvec(j,indf) = hes(i)
          END DO
       END DO

       CALL diamat_all(eigvec,eigval,dac=.FALSE.,error=error)
       CALL set_hes_eig(ndf,eigval,work)
       dx = eigval 
       CALL rat_fun_opt(ndf,dg,eigval,work,eigvec,g)
       CALL geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)
       CALL trust_radius(ndf,nfree,step,rad,rat,dr,output_unit)

       ! Update the atomic positions
       x0 = x0 + dr

       CALL energy_predict(ndf,work,hes,dr,g,conv,pred)
       eold = etot

       ! Energy & Gradients at new step
       CALL cp_eval_at(gopt_env, x0, etot, g, gopt_env%force_env%para_env%mepos,&
            gopt_env%force_env%para_env, error)
       ediff = etot - eold

       ! check for an external exit command
       CALL external_control(should_stop,"GEO",globenv,error)
       IF(should_stop) EXIT

       ! Some IO and Convergence check
       CALL gopt_f_io(gopt_env, force_env, root_section, its, etot,&
            output_unit, eold, emin, wildcard, gopt_param, ndf, dr, g, conv, ediff, pred,&
            rat, step, rad, error)

       IF (conv.OR.(its==maxiter)) EXIT
       IF (etot < emin) emin = etot
       CALL update_trust_rad(rat,rad,step,ediff)
    END DO

    IF(its == maxiter .AND. (.NOT.conv))THEN
       CALL print_geo_opt_nc(gopt_env, output_unit)
    END IF

    ! Write final  information, if converged
    CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=its,error=error)
    CALL gopt_f_io_finalize(gopt_env, force_env, x0, conv, its, root_section,&
         gopt_env%force_env%para_env, gopt_env%force_env%para_env%mepos,error)

    DEALLOCATE (hes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (xold,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (g,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (gold,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (dx,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (hes_mat,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (dg,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (eigval,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (eigvec,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (dr,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_print_key_finished_output(output_unit,logger,geo_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)
    CALL timestop(handle)

  END SUBROUTINE geoopt_bfgs

! *****************************************************************************
  SUBROUTINE rat_fun_opt(ndf,dg,eigval,work,eigvec,g)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: dg(ndf), eigval(ndf), &
                                                work(ndf), eigvec(ndf,ndf), &
                                                g(ndf)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rat_fun_opt', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp

    INTEGER                                  :: handle, indf, iref, iter, &
                                                maxit
    LOGICAL                                  :: bisec, fail, set
    REAL(KIND=dp)                            :: DDOT, fun, fun1, fun2, fun3, &
                                                fung, lam1, lam2, ln, lp, &
                                                ssize, step, stol

    CALL timeset(routineN,handle)

    stol  = 1.0E-8_dp
    ssize = 0.2_dp
    maxit = 999
    fail  = .FALSE.
    bisec = .FALSE.

    DO indf = 1, ndf
      dg(indf) =  DDOT(ndf,eigvec(1,indf),1,g,1)
    END DO

    set = .FALSE.

    DO

!   calculating Lamda

      lp = 0.0_dp
      iref = 1
      ln = 0.0_dp
      IF(eigval(iref) < 0.0_dp )ln = eigval(iref) - 0.01_dp

      iter = 0
      DO
         iter = iter + 1
         fun = 0.0_dp
         fung = 0.0_dp
         DO indf = 1, ndf
            fun = fun + dg(indf)**2/(ln-eigval(indf))
            fung = fung - dg(indf)**2 /(ln-eigval(indf)**2)
         END DO
         fun = fun - ln
         fung = fung - 1.0_dp
         step = fun/fung
         ln = ln - step
         IF(ABS(step) < stol )GOTO 200
         IF(iter >= maxit)EXIT
      END DO
100   CONTINUE
      bisec = .TRUE.
      iter = 0
      maxit = 9999
      lam1 = 0.0_dp
      IF(eigval(iref) < 0.0_dp )lam1 = eigval(iref) - 0.01_dp
      fun1 = 0.0_dp
      DO indf = 1, ndf
        fun1 = fun1 + dg(indf)**2 / (lam1-eigval(indf))
      END DO
      fun1 = fun1 - lam1
      step = ABS(lam1)/1000.0_dp
      IF(step < ssize) step = ssize
      DO
         iter = iter + 1
         IF(iter > maxit)THEN
           ln = 0.0_dp
           lp = 0.0_dp
           fail = .TRUE.
           GOTO 300
         END IF
         fun2 = 0.0_dp
         lam2 = lam1 - iter * step
         DO indf = 1, ndf
           fun2 = fun2 + eigval(indf)**2 / (lam2 - eigval(indf))
         END DO
         fun2 = fun2 - lam2
         IF(fun2*fun1 < 0.0_dp )THEN
           iter = 0
           DO
             iter = iter + 1
             IF(iter > maxit)THEN
               ln = 0.0_dp
               lp = 0.0_dp
               fail = .TRUE.
               GOTO 300
             END IF
             step = (lam1+lam2)/2
             fun3 = 0.0_dp
             DO indf = 1, ndf
              fun3 = fun3 + dg(indf)**2/(step-eigval(indf))
             END DO
             fun3 = fun3 - step

             IF( ABS(step-lam2) < stol)THEN
               ln = step
               GOTO 200
             END IF

             IF( fun3 * fun1 < stol )THEN
                lam2 = step
             ELSE
                lam1 = step
             END IF
           END DO
         END IF
     END DO

200  CONTINUE
     IF( (ln > eigval(iref)).OR. ( (ln > 0.0_dp) .AND. &
              (eigval(iref) > 0.0_dp) ) )THEN

         IF(.NOT. bisec)GOTO 100
         ln = 0.0_dp
         lp = 0.0_dp
         fail = .TRUE.
     END IF

300  CONTINUE

     IF(fail .AND. .NOT. set )THEN
        set = .TRUE.
        DO indf = 1, ndf
           eigval(indf) = eigval(indf) * work(indf)
        END DO
        GOTO 400
     END IF

     IF(.NOT. set) THEN
        work(1:ndf)=1.0_dp
     ENDIF

     DO indf = 1, ndf
       eigval(indf) = eigval(indf) - ln
     END DO
     EXIT
     END DO
400  CONTINUE

     CALL timestop(handle)

  END SUBROUTINE  rat_fun_opt

! *****************************************************************************
  SUBROUTINE bfgs(ndf,dx,dg,hes,work)
    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: dx(ndf), dg(ndf), &
                                                hes(ndf*(ndf+1)/2), work(ndf)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bfgs', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle, i, indf, j
    REAL(KIND=dp)                            :: DDOT, dxw, gdx

    CALL timeset(routineN,handle)

    CALL DSPMV('U',ndf,one,hes,dx,1,zero,work,1)

    gdx = DDOT(ndf,dg,1,dx,1)
    gdx = one/gdx
    dxw = DDOT(ndf,dx,1,work,1)
    dxw = one/dxw

    DO indf = 1, ndf
      DO j = 1, indf
         i = indf * (indf-1)/2 + j
         hes(i) = hes(i) + gdx * dg(indf) * dg(j) - &
                           dxw * work(indf) * work(j)
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE bfgs

! *****************************************************************************
  SUBROUTINE set_hes_eig(ndf,eigval,work)
    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: eigval(ndf), work(ndf)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_hes_eig', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: max_neg = -0.5_dp, &
                                                max_pos = 5.0_dp, &
                                                min_eig = 0.005_dp, &
                                                one = 1.0_dp

    INTEGER                                  :: handle, indf
    LOGICAL                                  :: neg

    CALL timeset(routineN,handle)

     DO indf = 1, ndf
       IF(eigval(indf) < 0.0_dp )neg = .TRUE.
       IF(eigval(indf) > 1000.0_dp)eigval(indf) = 1000.0_dp
    END DO
    DO indf = 1, ndf
       IF( eigval(indf) < 0.0_dp ) THEN
          IF(eigval(indf) < max_neg)THEN
             eigval(indf) = max_neg
          ELSE IF(eigval(indf) > - min_eig )THEN
             eigval(indf) = - min_eig
          END IF
       ELSE IF(eigval(indf) < 1000.0_dp)THEN
          IF(eigval(indf) < min_eig)THEN
             eigval(indf) = min_eig
          ELSE IF(eigval(indf) > max_pos)THEN
             eigval(indf) = max_pos
          END IF
       END IF
     END DO

     DO indf = 1, ndf
       IF( eigval(indf) < 0.0_dp )THEN
         work(indf) = - one
       ELSE
         work(indf) =  one
       END IF
     END DO

     CALL timestop(handle)

  END SUBROUTINE set_hes_eig

! *****************************************************************************
  SUBROUTINE geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: eigval(ndf), eigvec(ndf,ndf), &
                                                hes_mat(ndf,ndf), dr(ndf), &
                                                g(ndf)

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i, indf, j
    REAL(KIND=dp), ALLOCATABLE               :: tmp(:,:)

     hes_mat(:,:) = eigvec(:,:)

     DO indf=1,ndf
       eigval(indf) = one/eigval(indf)
       CALL dscal(ndf,eigval(indf),hes_mat(1,indf),1)
     END DO

    ALLOCATE(tmp(ndf,ndf))

     CALL dgemm ("N","T",ndf,ndf,ndf,one,hes_mat,ndf,eigvec,ndf,zero,&
                 tmp,ndf)

     DO i=1,ndf
       DO j=1,ndf
         hes_mat(j,i) = 0.5_dp*(tmp(j,i) + tmp(i,j))
         hes_mat(i,j) = hes_mat(j,i)
       END DO
     END DO

     DEALLOCATE(tmp)

!    ** New step **

     CALL DGEMV('N',ndf,ndf,-one,hes_mat,ndf,g,1,zero,dr,1)

  END SUBROUTINE geoopt_get_step

! *****************************************************************************
  SUBROUTINE trust_radius(ndf,nfree,step,rad,rat,dr,output_unit)
    INTEGER, INTENT(IN)                      :: ndf, nfree
    REAL(KIND=dp), INTENT(INOUT)             :: step, rad, rat, dr(ndf)
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'trust_radius', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: DDOT, scal

    CALL timeset(routineN,handle)

    step = DDOT(ndf,dr,1,dr,1)
    scal = MAX(one,LOG10(REAL(nfree,dp)))
    step = SQRT(step/scal)

    IF( step > rad )THEN
       rat = rad / step
       CALL DSCAL(ndf,rat,dr,1)
       step = DDOT(ndf,dr,1,dr,1)
       step = SQRT(step/scal)
       IF(output_unit>0) THEN
         WRITE(unit=output_unit,FMT="(/,T2,A,F8.5)") &
             " Step is scaled; Scaling factor = ", rat
         CALL m_flush(output_unit)
       ENDIF
    END IF
    CALL timestop(handle)

  END SUBROUTINE trust_radius

! *****************************************************************************
  SUBROUTINE energy_predict(ndf,work,hes,dr,g,conv,pred)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: work(ndf), &
                                                hes(ndf*(ndf+1)/2), dr(ndf), &
                                                g(ndf)
    LOGICAL, INTENT(INOUT)                   :: conv
    REAL(KIND=dp), INTENT(INOUT)             :: pred

    CHARACTER(LEN=*), PARAMETER :: routineN = 'energy_predict', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: limit = 1.0E-9_dp , &
                                                one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: DDOT, ener1, ener2

    CALL timeset(routineN,handle)

    ener1 = DDOT(ndf,g,1,dr,1)
    CALL DSPMV('U',ndf,one,hes,dr,1,zero,work,1)
    ener2 = DDOT(ndf,dr,1,work,1)
    pred  = ener1 + 0.5_dp * ener2
    conv = .FALSE.
    CALL timestop(handle)

  END SUBROUTINE energy_predict

! *****************************************************************************
   SUBROUTINE update_trust_rad(rat,rad,step,ediff)

    REAL(KIND=dp), INTENT(INOUT)             :: rat, rad, step, ediff

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_trust_rad', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: max_trust = 1.0_dp, &
                                                min_trust = 0.1_dp

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF( rat > 4.0_dp )THEN
      IF( ediff < 0.0_dp)THEN
        rad = step * 0.5_dp
      ELSE
        rad = step * 0.25_dp
      END IF
    ELSE IF ( rat > 2.0_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 0.75_dp
      ELSE
        rad = step * 0.5_dp
      END IF
    ELSE IF ( rat > 4.0_dp/3.0_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step
      ELSE
        rad = step * 0.75_dp
      END IF
    ELSE IF (rat > 10.0_dp/9.0_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 1.25_dp
      ELSE
        rad = step
      END IF
    ELSE IF( rat > 0.9_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 1.5_dp
      ELSE
        rad = step * 1.25_dp
      END IF
    ELSE IF( rat > 0.75_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 1.25_dp
      ELSE
        rad = step
      END IF
    ELSE IF( rat > 0.5_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step
      ELSE
        rad = step * 0.75_dp
      END IF
    ELSE IF( rat > 0.25_dp )THEN
      IF( ediff < 0.0_dp )THEN
        rad = step * 0.75_dp
      ELSE
        rad = step * 0.5_dp
      END IF
    ELSE IF( ediff < 0.0_dp )THEN
        rad = step * 0.5_dp
    ELSE
        rad = step * 0.25_dp
    END IF

    rad = MAX(rad,min_trust)
    rad = MIN(rad,max_trust)
    CALL timestop(handle)

  END SUBROUTINE update_trust_rad

END MODULE bfgs_optimizer
