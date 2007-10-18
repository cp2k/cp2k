!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/bfgs_optimizer [1.0] *
!!
!!   NAME
!!     bfgs_optimizer
!!
!!   FUNCTION
!!     Routines for Geometry optimization using BFGS algorithm
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************

MODULE bfgs_optimizer

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cell_types,                      ONLY: cell_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: cp_eval_at,&
                                             gopt_f_type,&
                                             gopt_param_type,&
                                             write_cycle_infos,&
                                             check_converg
  USE input_constants,                 ONLY: dump_xmol
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates
  USE termination,                     ONLY: external_control,&
                                             stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = 'bfgs_optimizer'
  LOGICAL, PARAMETER          :: debug_this_module=.TRUE.

  PUBLIC :: geoopt_bfgs,&
            write_cycle_infos,&
            check_converg

CONTAINS

!!****** bfgs_optimizer/geoopt_bfgs [1.0] *
!!
!!   NAME
!!    geoopt_bfgs 
!!
!!   FUNCTION
!!    Main driver for BFGS geometry optimizations 
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
  SUBROUTINE geoopt_bfgs(force_env,gopt_param,root_section,globenv,force_env_section,&
       gopt_env, x0, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'geoopt_bfgs'
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    CHARACTER(LEN=80)                        :: title
    INTEGER :: handle, handle2, hesunit, i, ihandle, indf, istat, it, j, &
      maxiter, natom, ndf, nfree, nkind, output_unit, traj_unit
    LOGICAL                                  :: conv, hesrest, ionode, &
                                                should_stop
    REAL(KIND=dp)                            :: ediff, emin, eold, etot, &
                                                pred, rad, rat, step
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dg, dr, dx, eigval, gold, &
                                                hes, work, xold
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: eigvec, hes_mat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    ! -------------------------------------------------------------------------

    NULLIFY(logger, g)
    logger => cp_error_get_logger(error)
    para_env=>force_env%para_env

    CALL cp_add_iter_level(logger%iter_info,"GEO_OPT",ihandle,error)
    CALL timeset(routineN,"I","",handle)

    ionode      = para_env%mepos==para_env%source
    maxiter     = gopt_param%max_iter
    hesrest     = gopt_param%hesrest
    conv        = .FALSE.
    rad         = 0.10_dp
    rat         = 0.0_dp

    output_unit = cp_print_key_unit_nr(logger,root_section,&
         "MOTION%GEO_OPT%PRINT%PROGRAM_RUN_INFO",&
         extension=".geoLog",error=error)

    CALL force_env_get(force_env,&
         cell=cell,&
         subsys=subsys,error=error)

    CALL cp_subsys_get(subsys=subsys,&
         atomic_kinds=atomic_kinds,&
         particles=particles,error=error)

    nkind = atomic_kinds%n_els

    natom = particles%n_els
    particle_set => particles%els

    ndf = SIZE(x0)
    nfree = SIZE(particle_set)

    ! Initialize hessian (hes = unitary matrix )

    ALLOCATE(hes(ndf*(ndf+1)/2))

    DO indf=1,ndf
      DO j=1,indf
        i = indf*(indf - 1)/2 + j
        hes(i) = zero
      END DO
      i = indf*(indf-1)/2 + indf
      hes(i) = one
    END DO

    ALLOCATE (xold(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "xold",ndf*dp_size)
    xold(:) = x0(:)

    ALLOCATE (g(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "g",ndf*dp_size)
    g(:) = 0.0_dp

    ALLOCATE (gold(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gold",ndf*dp_size)
    gold(:) = 0.0_dp

    ALLOCATE (dx(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "dx",ndf*dp_size)
    dx(:) = 0.0_dp

    ALLOCATE (hes_mat(ndf,ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "hes_mat",ndf*ndf*dp_size)
    hes_mat(:,:) = 0.0_dp

    ALLOCATE (dg(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "dg",ndf*dp_size)
    dg(:) = 0.0_dp

    ALLOCATE (eigval(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "eigval",ndf*dp_size)
    eigval(:) = 0.0_dp

    ALLOCATE (eigvec(ndf,ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "eigvec",ndf*ndf*dp_size)
    eigvec(:,:) = 0.0_dp

    ALLOCATE (work(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",ndf*dp_size)
    work(:) = 0.0_dp

    ALLOCATE (dr(ndf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "dr",ndf*dp_size)
    dr(:) = 0.0_dp

    ! Geometry optimization starts now

    IF (output_unit > 0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
      WRITE (UNIT=output_unit,FMT="(T2,A,T28,A,T78,A)")&
        "***","STARTING GEOMETRY OPTIMIZATION","***"
      WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
      CALL m_flush(output_unit)
    END IF

    IF (ionode) THEN
      CALL open_file(file_name="HESSIAN",&
                     file_status="UNKNOWN",&
                     file_form="UNFORMATTED",&
                     file_action="READWRITE",&
                     unit_number=hesunit)
    END IF

    ! Calculate Energy & Gradients
    CALL cp_eval_at(gopt_env, x0, etot, g, gopt_env%force_env%para_env%mepos,&
            gopt_env%force_env%para_env, force_env_section, error)
    !MK Hm, it seems that here only the root process has the correct gradient for all atoms
    CALL mp_bcast(g,para_env%source,para_env%group)

    ! Print info at time 0
    emin = etot
    it = 0
    IF(output_unit>0) CALL write_cycle_infos &
         (output_unit,it,etot,wildcard=" BFGS")

    DO

       it = it + 1

       CALL cp_iterate(logger%iter_info,iter_nr=it,error=error)

       IF (it > maxiter ) THEN
         IF (output_unit > 0) THEN
           WRITE(UNIT=output_unit,FMT="(/,T2,A)") REPEAT("*",79)
           WRITE (UNIT=output_unit,FMT="(T2,A,T20,A,T78,A)")&
             "***","MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED","***"
           WRITE (UNIT=output_unit,FMT="(T2,A,T27,A,T78,A)")&
             "***","EXITING GEOMETRY OPTIMIZATION","***"
           WRITE(UNIT=output_unit,FMT="(T2,A)") REPEAT("*",79)
           CALL m_flush(output_unit)
         END IF
         EXIT
       END IF

       ! check for an external exit command
       CALL external_control(should_stop,"GEO",globenv,error)
       IF(should_stop) EXIT
       CALL write_restart(force_env=force_env,root_section=root_section,&
            force_env_section=force_env_section,error=error)

       IF (ionode) THEN
         traj_unit = cp_print_key_unit_nr(logger,root_section,&
                                          "MOTION%PRINT%TRAJECTORY",&
                                          extension=".xyz",&
                                          error=error)
         IF (traj_unit > 0) THEN
           CALL force_env_get(force_env,&
                              cell=cell,&
                              subsys=subsys,error=error)
           CALL cp_subsys_get(subsys=subsys,&
                              particles=particles,error=error)
           particle_set => particles%els
           WRITE (UNIT=title,FMT="(A,I8,A,F20.10)") " i =",it,", E =",etot
           CALL write_particle_coordinates(particle_set,traj_unit,&
                                           dump_xmol,"POS",TRIM(title))
           CALL m_flush(traj_unit)
         END IF
         CALL cp_print_key_finished_output(traj_unit,logger,root_section,&
                                           "MOTION%PRINT%TRAJECTORY",error=error)
       END IF

       IF ( output_unit > 0) THEN
          WRITE(UNIT=output_unit,FMT="(/,T2,26('-'))")
          WRITE(UNIT=output_unit,FMT="(T2,A,I6)")&
               "OPTIMIZATION STEP: ",it
          WRITE(UNIT=output_unit,FMT="(T2,26('-'))")
          CALL m_flush(output_unit)
       END IF

       ! Hessian update/restarting
       IF ((it == 1).AND.hesrest) THEN
          IF(ionode)THEN
             READ(UNIT=hesunit,IOSTAT=istat) hes
             REWIND(hesunit)
             IF (istat /= 0) THEN
                WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                     "**  Error while reading HESSIAN **"
                conv = .TRUE.
             ELSE
                WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                     "*** Initial Hessian has been restarted successfully ***"
             END IF

          END IF
          CALL mp_bcast(hes,para_env%source,para_env%group)
       ELSE IF( it > 1 ) THEN
          DO indf = 1, ndf
             dx(indf) = x0(indf) - xold(indf)
             dg(indf) =  g(indf) - gold(indf)
          END DO

          CALL bfgs(ndf,dx,dg,hes,work)

          CALL timeset(routineN//"hes2disk","I","",handle2)
          IF(ionode)THEN
             WRITE (UNIT=hesunit) hes
             CALL m_flush(hesunit)
             REWIND (hesunit)
          END IF
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

       CALL diamat_all(eigvec,eigval,dac=.FALSE.)
       CALL set_hes_eig(ndf,eigval,work)
       CALL DCOPY(ndf,eigval,1,dx,1)
       CALL rat_fun_opt(ndf,dg,eigval,work,eigvec,g)
       CALL geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)
       CALL trust_radius(ndf,nfree,step,rad,rat,dr,output_unit,output_unit>0)

       ! Update the atomic positions
       x0 = x0 + dr

       CALL energy_predict(ndf,work,hes,dr,g,output_unit>0,output_unit,conv,pred)
       IF (conv) EXIT
       eold = etot

       ! Energy & Gradients at new step
       CALL cp_eval_at(gopt_env, x0, etot, g, gopt_env%force_env%para_env%mepos,&
            gopt_env%force_env%para_env, force_env_section, error)
       !MK Hm, it seems that here only the root process has the correct gradient for all atoms
       CALL mp_bcast(g,para_env%source,para_env%group)
       ediff = etot - eold

       IF (output_unit > 0) THEN
         CALL write_cycle_infos(output_unit,it,etot,ediff,pred,rat,step,rad,emin," BFGS")
       END IF

       IF (etot < emin) emin = etot
       CALL update_trust_rad(rat,rad,step,ediff)
       CALL check_converg(ndf,dr,g,output_unit>0,output_unit,conv,gopt_param)
       IF (conv) EXIT

    END DO

    ! Write final particle information, if converged

    IF (conv) THEN
      CALL write_restart(force_env=force_env,root_section=root_section,&
                          force_env_section=force_env_section,error=error)
      IF (ionode) THEN
        traj_unit = cp_print_key_unit_nr(logger,root_section,&
                                         "MOTION%PRINT%TRAJECTORY",&
                                         extension=".xyz",&
                                         error=error)
        IF (traj_unit > 0) THEN
          CALL force_env_get(force_env,&
                             cell=cell,&
                             subsys=subsys,error=error)
          CALL cp_subsys_get(subsys=subsys,&
                             particles=particles,error=error)
          particle_set => particles%els
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)") " i =",it,", E =",etot
          CALL write_particle_coordinates(particle_set,traj_unit,&
                                          dump_xmol,"POS",TRIM(title))
          CALL m_flush(traj_unit)
        END IF
        CALL cp_print_key_finished_output(traj_unit,logger,root_section,&
                                          "MOTION%PRINT%TRAJECTORY",error=error)
      END IF
    END IF

    DEALLOCATE (hes,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"hes")
    DEALLOCATE (xold,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"xold")
    DEALLOCATE (g,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"g")
    DEALLOCATE (gold,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"gold")
    DEALLOCATE (dx,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dx")
    DEALLOCATE (hes_mat,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"hes_mat")
    DEALLOCATE (dg,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dg")
    DEALLOCATE (eigval,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eigval")
    DEALLOCATE (eigvec,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eigvec")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")
    DEALLOCATE (dr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dr")

    IF (ionode) THEN
      CALL close_file(unit_number=hesunit)
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,root_section,&
         "MOTION%GEO_OPT%PRINT%PROGRAM_RUN_INFO", error=error)

    CALL cp_rm_iter_level(logger%iter_info,ihandle,error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE geoopt_bfgs

!!****** bfgs_optimizer/rat_fun_opt [1.0] *
!!
!!   NAME
!!    rat_fun_opt 
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
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
        ! CALL DCOPY(ndf,1.0,0,work,1) ! bug
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

!!****** bfgs_optimizer/bfgs [1.0] *
!!
!!   NAME
!!    bfgs 
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
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

!!****** bfgs_optimizer/set_hes_eig [1.0] *
!!
!!   NAME
!!    set_hes_eig 
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
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

!!****** bfgs_optimizer/geoopt_get_step [1.0] *
!!
!!   NAME
!!    geoopt_get_step 
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
  SUBROUTINE geoopt_get_step(ndf,eigval,eigvec,hes_mat,dr,g)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: eigval(ndf), eigvec(ndf,ndf), &
                                                hes_mat(ndf,ndf), dr(ndf), &
                                                g(ndf)

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i, indf, j
    REAL(KIND=dp)                            :: tmp(ndf,ndf)

     hes_mat(:,:) = eigvec(:,:)

     DO indf=1,ndf
       eigval(indf) = one/eigval(indf)
       CALL dscal(ndf,eigval(indf),hes_mat(1,indf),1)
     END DO

     CALL dgemm ("N","T",ndf,ndf,ndf,one,hes_mat,ndf,eigvec,ndf,zero,&
                 tmp,ndf)

     DO i=1,ndf
       DO j=1,ndf
         hes_mat(j,i) = 0.5_dp*(tmp(j,i) + tmp(i,j))
         hes_mat(i,j) = hes_mat(j,i)
       END DO
     END DO

!    ** New step **

     CALL DGEMV('N',ndf,ndf,-one,hes_mat,ndf,g,1,zero,dr,1)

  END SUBROUTINE geoopt_get_step

!!****** bfgs_optimizer/trust_radius [1.0] *
!!
!!   NAME
!!    trust_radius 
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
  SUBROUTINE trust_radius(ndf,nfree,step,rad,rat,dr,output_unit,ionode)
    INTEGER, INTENT(IN)                      :: ndf, nfree
    REAL(KIND=dp), INTENT(INOUT)             :: step, rad, rat, dr(ndf)
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: ionode

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
       IF(ionode) THEN
         WRITE(unit=output_unit,FMT="(/,T2,A,F8.5)") &
             " Step is scaled; Scaling factor = ", rat
         CALL m_flush(output_unit)
       ENDIF
    END IF
    CALL timestop(handle)

  END SUBROUTINE trust_radius

!!****** bfgs_optimizer/energy_predict [1.0] *
!!
!!   NAME
!!    energy_predict  
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
   SUBROUTINE energy_predict(ndf,work,hes,dr,g,ionode,output_unit,conv,pred)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: work(ndf), &
                                                hes(ndf*(ndf+1)/2), dr(ndf), &
                                                g(ndf)
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
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

!!****** bfgs_optimizer/update_trust_rad [1.0] *
!!
!!   NAME
!!     update_trust_rad 
!!
!!   FUNCTION
!!    
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************
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
