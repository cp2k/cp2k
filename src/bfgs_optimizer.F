!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for Geometry optimization using BFGS algorithm
! *****************************************************************************
MODULE bfgs_optimizer

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_external_control,             ONLY: external_control
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_transpose
  USE cp_fm_diag,                      ONLY: choose_eigv_solver
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_read_unformatted,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type,&
                                             cp_fm_write_unformatted
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_iterate,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_methods,                  ONLY: gopt_f_ii,&
                                             gopt_f_io,&
                                             gopt_f_io_finalize,&
                                             gopt_f_io_init,&
                                             print_geo_opt_header,&
                                             print_geo_opt_nc
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_ts_method_id
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_min,&
                                             mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bfgs_optimizer'
  LOGICAL, PARAMETER                       :: debug_this_module=.TRUE.

  PUBLIC :: geoopt_bfgs

CONTAINS

! *****************************************************************************
!> \brief Main driver for BFGS geometry optimizations
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \param error ...
! *****************************************************************************
  RECURSIVE SUBROUTINE geoopt_bfgs(force_env,gopt_param,globenv,geo_section,gopt_env,x0,&
                         error)

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
    INTEGER                                  :: handle, hesunit_read, indf, &
                                                info, iter_nr, its, maxiter, &
                                                ndf, nfree, output_unit, stat
    LOGICAL                                  :: cell_opt, conv, failure, &
                                                hesrest, ionode, &
                                                shell_present, should_stop, &
                                                use_mod_hes, use_rfo
    REAL(KIND=dp)                            :: ediff, emin, eold, etot, &
                                                pred, rad, rat, step, t_diff, &
                                                t_now, t_old
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dg, dr, dx, eigval, gold, &
                                                work, xold
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_hes
    TYPE(cp_fm_type), POINTER                :: eigvec_mat, hess_mat, hess_tmp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: print_key, root_section

    failure = .FALSE.
    NULLIFY(logger, g,blacs_env)
    logger       => cp_error_get_logger(error)
    para_env     => force_env%para_env
    root_section => force_env%root_section
    t_old=m_walltime()

    CALL timeset(routineN,handle)
    CALL section_vals_val_get(geo_section,"BFGS%TRUST_RADIUS",r_val=rad,error=error)
    print_key    => section_vals_get_subs_vals(geo_section,"BFGS%RESTART", error=error)
    ionode      = para_env%mepos==para_env%source
    maxiter     = gopt_param%max_iter
    cell_opt    = .FALSE.
    conv        = .FALSE.
    rat         = 0.0_dp
    wildcard    = " BFGS"

    ! Stop if not yet implemented
    SELECT CASE (gopt_env%type_id)
    CASE (default_ts_method_id)
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="BFGS method not yet working with DIMER", &
            error=error, error_level=cp_failure_level)
    END SELECT

    CALL section_vals_val_get(geo_section,"BFGS%USE_RAT_FUN_OPT",l_val=use_rfo,error=error)
    CALL section_vals_val_get(geo_section,"BFGS%USE_MODEL_HESSIAN",l_val=use_mod_hes,error=error)
    CALL section_vals_val_get(geo_section,"BFGS%RESTART_HESSIAN",l_val=hesrest,error=error)
    output_unit = cp_print_key_unit_nr(logger,geo_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".geoLog",error=error)

    ndf   = SIZE(x0)
    nfree = gopt_env%nfree
    CALL cp_assert((ndf <= 3000),cp_note_level,cp_assertion_failed,routineP,&
                   "The dimension of the Hessian matrix ("//&
                   TRIM(ADJUSTL(cp_to_string(ndf)))//") is greater than 3000. "//&
                   "The diagonalisation of the full Hessian  matrix needed for BFGS "//&
                   "is computationally expensive. You should consider to use the linear "//&
                   "scaling variant L-BFGS instead.",&
                   only_ionode=.TRUE.)

    ! Initialize hessian (hes = unitary matrix or model hessian )
    CALL cp_blacs_env_create(blacs_env,para_env,globenv%blacs_grid_layout,&
                             globenv%blacs_repeatable,error=error)
    CALL cp_fm_struct_create(fm_struct_hes,para_env=para_env,context=blacs_env, &
         nrow_global=ndf,ncol_global=ndf,error=error)
    CALL cp_fm_create(hess_mat, fm_struct_hes,name="hess_mat",error=error)
    CALL cp_fm_create(hess_tmp, fm_struct_hes,name="hess_tmp",error=error)
    CALL cp_fm_create(eigvec_mat, fm_struct_hes,name="eigvec_mat",error=error)
    ALLOCATE (eigval(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    eigval(:) = 0.0_dp

    CALL force_env_get(force_env=force_env,subsys=subsys,error=error)
    CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,error=error)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present)
    IF (shell_present) THEN
      CALL cp_assert((.NOT.use_mod_hes),cp_warning_level,cp_assertion_failed,routineP,&
                     "No model Hessian available for core-shell models. "//&
                     "An unit matrix is used for the initial Hessian.",&
                      only_ionode=.TRUE.)
      use_mod_hes = .FALSE.
    END IF

    IF (use_mod_hes) THEN
       CALL cp_fm_set_all(hess_mat,alpha=zero,beta=0.00_dp,error=error)
       CALL construct_initial_hess(gopt_env%force_env,hess_mat,error)
       CALL cp_fm_to_fm(hess_mat,hess_tmp,error)
       CALL choose_eigv_solver(hess_tmp,eigvec_mat,eigval,info=info, error=error)
       ! In rare cases the diagonalization of hess_mat fails (bug in scalapack?)
       IF(info/=0) THEN
          CALL cp_fm_set_all(hess_mat,alpha=zero,beta=one,error=error)
          IF(output_unit>0) WRITE(output_unit,*) &
             "BFGS: Matrix diagonalization failed, using unity as model Hessian."
       ELSE
          DO its=1,SIZE(eigval)
            IF(eigval(its).lt.0.1_dp)eigval(its)=0.1_dp
          END DO
          CALL cp_fm_to_fm(eigvec_mat,hess_tmp,error)
          CALL cp_fm_column_scale(eigvec_mat,eigval)
          CALL cp_gemm("N","T",ndf,ndf,ndf,one,hess_tmp,eigvec_mat,zero,hess_mat,error=error)
       ENDIF
    ELSE
       CALL cp_fm_set_all(hess_mat,alpha=zero,beta=one,error=error)
    END IF

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

    ALLOCATE (dg(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dg(:) = 0.0_dp

    ALLOCATE (work(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    work(:) = 0.0_dp

    ALLOCATE (dr(ndf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dr(:) = 0.0_dp

    ! refernce cell needed for smooth cell optimizations
    CALL force_env_get(gopt_env%force_env,cell=cell,error=error)

    ! Geometry optimization starts now
    CALL cp_iterate(logger%iter_info,increment=0,iter_nr_out=iter_nr,error=error)
    CALL print_geo_opt_header(gopt_env, output_unit, wildcard)

    ! Calculate Energy & Gradients
    CALL cp_eval_at(gopt_env, x0, etot, g, gopt_env%force_env%para_env%mepos,&
                    .FALSE., gopt_env%force_env%para_env, error)

    ! Print info at time 0
    emin = etot
    t_now=m_walltime()
    t_diff=t_now-t_old
    t_old=t_now
    CALL gopt_f_io_init(gopt_env, output_unit, etot, wildcard=wildcard, its=iter_nr, used_time=t_diff, error=error)
    DO its = iter_nr+1, maxiter
       CALL cp_iterate(logger%iter_info,last=(its==maxiter),error=error)
       CALL section_vals_val_set(geo_section,"STEP_START_VAL",i_val=its,error=error)
       CALL gopt_f_ii(its, output_unit)

       ! Hessian update/restarting
       IF (((its-iter_nr) == 1).AND.hesrest) THEN
          IF(ionode)THEN
             CALL section_vals_val_get(geo_section,"BFGS%RESTART_FILE_NAME",c_val=hes_filename,&
                  error=error)
             CALL open_file(file_name=hes_filename,file_status="OLD",&
                  file_form="UNFORMATTED", file_action="READ", unit_number=hesunit_read)
          END IF
          CALL cp_fm_read_unformatted(hess_mat,hesunit_read,error)
          IF (ionode) CALL close_file(unit_number=hesunit_read)
       ELSE
          IF( (its-iter_nr) > 1 ) THEN
            DO indf = 1, ndf
               dx(indf) = x0(indf) - xold(indf)
               dg(indf) =  g(indf) - gold(indf)
            END DO

            CALL bfgs(ndf,dx,dg,hess_mat,work,para_env,error)

            !Possibly dump the Hessian file
            IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file))THEN
               CALL write_bfgs_hessian(geo_section,ndf,hess_mat,logger,error)
            ENDIF
          ENDIF
       END IF

       ! Setting the present positions & gradients as old
       xold(:) = x0
       gold(:) = g

       ! Copying hessian hes to (ndf x ndf) matrix hes_mat for diagonalization
       CALL cp_fm_to_fm(hess_mat,hess_tmp,error)

       CALL choose_eigv_solver(hess_tmp,eigvec_mat,eigval,info=info,error=error)

       ! In rare cases the diagonalization of hess_mat fails (bug in scalapack?)
       IF(info/=0) THEN
          IF(output_unit>0) WRITE(output_unit,*) &
             "BFGS: Matrix diagonalization failed, resetting Hessian to unity."
          CALL cp_fm_set_all(hess_mat,alpha=zero,beta=one,error=error)
          CALL cp_fm_to_fm(hess_mat,hess_tmp,error)
          CALL choose_eigv_solver(hess_tmp,eigvec_mat,eigval,error=error)
       END IF

       IF(use_rfo)THEN
          CALL set_hes_eig(ndf,eigval,work)
          dx(:) = eigval
          CALL rat_fun_opt(ndf,dg,eigval,work,eigvec_mat,g,para_env,error)
       END IF
       CALL geoopt_get_step(ndf,eigval,eigvec_mat,hess_tmp,dr,g,para_env,use_rfo,error)
       CALL trust_radius(ndf,nfree,step,rad,rat,dr,output_unit)

       ! Update the atomic positions
       x0 = x0 + dr

       CALL energy_predict(ndf,work,hess_mat,dr,g,conv,pred,para_env,error)
       eold = etot

       ! Energy & Gradients at new step
       CALL cp_eval_at(gopt_env, x0, etot, g, gopt_env%force_env%para_env%mepos,&
                       .FALSE., gopt_env%force_env%para_env, error)
       ediff = etot - eold

       ! check for an external exit command
       CALL external_control(should_stop,"GEO",globenv=globenv,error=error)
       IF(should_stop) EXIT

       ! Some IO and Convergence check
       t_now=m_walltime()
       t_diff=t_now-t_old
       t_old=t_now
       CALL gopt_f_io(gopt_env, force_env, root_section, its, etot, output_unit,&
            eold, emin, wildcard, gopt_param, ndf, dr, g, conv, ediff, pred, rat,&
            step, rad, used_time=t_diff, error=error)

       IF (conv.OR.(its==maxiter)) EXIT
       IF (etot < emin) emin = etot
       IF (use_rfo) CALL update_trust_rad(rat,rad,step,ediff)
    END DO

    IF (its == maxiter .AND. (.NOT.conv))THEN
       CALL print_geo_opt_nc(gopt_env, output_unit)
    END IF

    ! Write final  information, if converged
    CALL cp_iterate(logger%iter_info,last=.TRUE.,increment=0,error=error)
    CALL write_bfgs_hessian(geo_section,ndf,hess_mat,logger,error)
    CALL gopt_f_io_finalize(gopt_env, force_env, x0, conv, its, root_section,&
         gopt_env%force_env%para_env, gopt_env%force_env%para_env%mepos, output_unit, error)

    CALL cp_fm_struct_release(fm_struct_hes,error=error)
    CALL cp_fm_release(hess_mat,error=error)
    CALL cp_fm_release(eigvec_mat,error=error)
    CALL cp_fm_release(hess_tmp,error=error)

    CALL cp_blacs_env_release(blacs_env,error=error)
    DEALLOCATE (xold,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (g,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (gold,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (dx,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (dg,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (eigval,STAT=stat)
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
!> \brief ...
!> \param ndf ...
!> \param dg ...
!> \param eigval ...
!> \param work ...
!> \param eigvec_mat ...
!> \param g ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rat_fun_opt(ndf,dg,eigval,work,eigvec_mat,g,para_env,error)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: dg(ndf), eigval(ndf), &
                                                work(ndf)
    TYPE(cp_fm_type), POINTER                :: eigvec_mat
    REAL(KIND=dp), INTENT(INOUT)             :: g(ndf)
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rat_fun_opt', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp

    INTEGER                                  :: handle, i, indf, iref, iter, &
                                                j, k, l, maxit, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: bisec, fail, set
    REAL(KIND=dp)                            :: fun, fun1, fun2, fun3, fung, &
                                                lam1, lam2, ln, lp, ssize, &
                                                step, stol
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data

    CALL timeset(routineN,handle)

    stol  = 1.0E-8_dp
    ssize = 0.2_dp
    maxit = 999
    fail  = .FALSE.
    bisec = .FALSE.

    dg=0._dp

    CALL cp_fm_get_info(eigvec_mat,row_indices=row_indices,col_indices=col_indices, &
         local_data=local_data,nrow_local=nrow_local,ncol_local=ncol_local,error=error)

    DO i=1,nrow_local
       j=row_indices(i)
       DO k=1,ncol_local
          l=col_indices(k)
          dg(l)=dg(l)+local_data(i,k)*g(j)
       END DO
    END DO
    CALL mp_sum(dg,para_env%group)

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
400  CONTINUE
     END DO

     CALL timestop(handle)

  END SUBROUTINE  rat_fun_opt

! *****************************************************************************
!> \brief ...
!> \param ndf ...
!> \param dx ...
!> \param dg ...
!> \param hess_mat ...
!> \param work ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE bfgs(ndf,dx,dg,hess_mat,work,para_env,error)
    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: dx(ndf), dg(ndf)
    TYPE(cp_fm_type), POINTER                :: hess_mat
    REAL(KIND=dp), INTENT(INOUT)             :: work(ndf)
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bfgs', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle, i, j, k, l, &
                                                ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: DDOT, dxw, gdx
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_hes

    CALL timeset(routineN,handle)

    CALL cp_fm_get_info(hess_mat,row_indices=row_indices,col_indices=col_indices, &
         local_data=local_hes,nrow_local=nrow_local,ncol_local=ncol_local,error=error)

    work=0.0_dp
    DO i=1,nrow_local
       j=row_indices(i)
       DO k=1,ncol_local
          l=col_indices(k)
          work(j)=work(j)+local_hes(i,k)*dx(l)
       END DO
    END DO

    CALL mp_sum(work,para_env%group)

    gdx = DDOT(ndf,dg,1,dx,1)
    gdx = one/gdx
    dxw = DDOT(ndf,dx,1,work,1)
    dxw = one/dxw

    DO i=1,nrow_local
       j=row_indices(i)
       DO k=1,ncol_local
          l=col_indices(k)
          local_hes(i,k)=local_hes(i,k)+gdx*dg(j)*dg(l) -&
                                        dxw*work(j)*work(l)
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE bfgs

! *****************************************************************************
!> \brief ...
!> \param ndf ...
!> \param eigval ...
!> \param work ...
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
!> \brief ...
!> \param ndf ...
!> \param eigval ...
!> \param eigvec_mat ...
!> \param hess_tmp ...
!> \param dr ...
!> \param g ...
!> \param para_env ...
!> \param use_rfo ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE geoopt_get_step(ndf,eigval,eigvec_mat,hess_tmp,dr,g,para_env,use_rfo,error)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: eigval(ndf)
    TYPE(cp_fm_type), POINTER                :: eigvec_mat, hess_tmp
    REAL(KIND=dp), INTENT(INOUT)             :: dr(ndf), g(ndf)
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    LOGICAL                                  :: use_rfo
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: i, indf, j, k, l, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_fm_type), POINTER                :: tmp

    CALL cp_fm_to_fm(eigvec_mat,hess_tmp,error)
    IF(use_rfo)THEN
       DO indf=1,ndf
         eigval(indf) = one/eigval(indf)
       END DO
    ELSE
       DO indf=1,ndf
          eigval(indf) = one/MAX(0.0001_dp,eigval(indf))
       END DO
    END IF

    CALL cp_fm_column_scale(hess_tmp,eigval)
    CALL cp_fm_get_info(eigvec_mat,matrix_struct=matrix_struct,error=error)
    CALL cp_fm_create(tmp, matrix_struct ,name="tmp",error=error)

    CALL cp_gemm("N","T",ndf,ndf,ndf,one,hess_tmp,eigvec_mat,zero,tmp,error=error)

    CALL cp_fm_transpose(tmp,hess_tmp,error)
    CALL cp_fm_release(tmp,error)

!    ** New step **

    CALL cp_fm_get_info(hess_tmp,row_indices=row_indices,col_indices=col_indices, &
         local_data=local_data,nrow_local=nrow_local,ncol_local=ncol_local,error=error)

    dr=0.0_dp
    DO i=1,nrow_local
       j=row_indices(i)
       DO k=1,ncol_local
          l=col_indices(k)
          dr(j)=dr(j)-local_data(i,k)*g(l)
       END DO
    END DO

    CALL mp_sum(dr,para_env%group)

  END SUBROUTINE geoopt_get_step

! *****************************************************************************
!> \brief ...
!> \param ndf ...
!> \param nfree ...
!> \param step ...
!> \param rad ...
!> \param rat ...
!> \param dr ...
!> \param output_unit ...
! *****************************************************************************
  SUBROUTINE trust_radius(ndf,nfree,step,rad,rat,dr,output_unit)
    INTEGER, INTENT(IN)                      :: ndf, nfree
    REAL(KIND=dp), INTENT(INOUT)             :: step, rad, rat, dr(ndf)
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'trust_radius', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: scal

    CALL timeset(routineN,handle)

    step = MAXVAL(ABS(dr))
    scal = MAX(one,rad/step)

    IF( step > rad )THEN
       rat = rad / step
       CALL DSCAL(ndf,rat,dr,1)
       step = rad
       IF(output_unit>0) THEN
         WRITE(unit=output_unit,FMT="(/,T2,A,F8.5)") &
             " Step is scaled; Scaling factor = ", rat
         CALL m_flush(output_unit)
       ENDIF
    END IF
    CALL timestop(handle)

  END SUBROUTINE trust_radius

! *****************************************************************************
!> \brief ...
!> \param ndf ...
!> \param work ...
!> \param hess_mat ...
!> \param dr ...
!> \param g ...
!> \param conv ...
!> \param pred ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE energy_predict(ndf,work,hess_mat,dr,g,conv,pred,para_env,error)

    INTEGER, INTENT(IN)                      :: ndf
    REAL(KIND=dp), INTENT(INOUT)             :: work(ndf)
    TYPE(cp_fm_type), POINTER                :: hess_mat
    REAL(KIND=dp), INTENT(INOUT)             :: dr(ndf), g(ndf)
    LOGICAL, INTENT(INOUT)                   :: conv
    REAL(KIND=dp), INTENT(INOUT)             :: pred
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'energy_predict', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: limit = 1.0E-9_dp , &
                                                one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle, i, j, k, l, &
                                                ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: DDOT, ener1, ener2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data

    CALL timeset(routineN,handle)

    ener1 = DDOT(ndf,g,1,dr,1)

    CALL cp_fm_get_info(hess_mat,row_indices=row_indices,col_indices=col_indices, &
         local_data=local_data,nrow_local=nrow_local,ncol_local=ncol_local,error=error)

    work=0.0_dp
    DO i=1,nrow_local
       j=row_indices(i)
       DO k=1,ncol_local
          l=col_indices(k)
          work(j)=work(j)+local_data(i,k)*dr(l)
       END DO
    END DO

    CALL mp_sum(work,para_env%group)
    ener2 = DDOT(ndf,dr,1,work,1)
    pred  = ener1 + 0.5_dp * ener2
    conv = .FALSE.
    CALL timestop(handle)

  END SUBROUTINE energy_predict

! *****************************************************************************
!> \brief ...
!> \param rat ...
!> \param rad ...
!> \param step ...
!> \param ediff ...
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

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param geo_section ...
!> \param ndf ...
!> \param hess_mat ...
!> \param logger ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE write_bfgs_hessian(geo_section,ndf,hess_mat,logger,error)

    TYPE(section_vals_type), POINTER         :: geo_section
    INTEGER                                  :: ndf
    TYPE(cp_fm_type), POINTER                :: hess_mat
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_bfgs_hessian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, hesunit

    CALL timeset(routineN,handle)

    hesunit = cp_print_key_unit_nr(logger,geo_section,"BFGS%RESTART",&
           extension=".Hessian",file_form="UNFORMATTED",file_action="WRITE",&
           file_position="REWIND",error=error)

    CALL cp_fm_write_unformatted(hess_mat,hesunit,error)

    CALL cp_print_key_finished_output(hesunit,logger,geo_section,"BFGS%RESTART", error=error)

    CALL timestop(handle)

  END SUBROUTINE write_bfgs_hessian
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param hess_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE construct_initial_hess(force_env,hess_mat,error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_fm_type), POINTER                :: hess_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'construct_initial_hess', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, iat_col, iat_row, iglobal, iind, j, jat_row, jglobal, jind, &
      k, natom, ncol_local, nrow_local, stat, z
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: at_row
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: d_ij, rho_ij
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r_ij
    REAL(KIND=dp), DIMENSION(3, 3)           :: alpha, r0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fixed, local_data
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
     CALL force_env_get(force_env=force_env,subsys=subsys,cell=cell,error=error)
     CALL cp_subsys_get(subsys,&
                 particles=particles,&
                 error=error)

     alpha(1,:)=(/1._dp,0.3949_dp,0.3949_dp/)
     alpha(2,:)=(/0.3494_dp,0.2800_dp,0.2800_dp/)
     alpha(3,:)=(/0.3494_dp,0.2800_dp,0.1800_dp/)

     r0(1,:)=(/1.35_dp,2.10_dp,2.53_dp/)
     r0(2,:)=(/2.10_dp,2.87_dp,3.40_dp/)
     r0(3,:)=(/2.53_dp,3.40_dp,3.40_dp/)

     CALL cp_fm_get_info(hess_mat,row_indices=row_indices,col_indices=col_indices, &
            local_data=local_data,nrow_local=nrow_local,ncol_local=ncol_local,error=error)
     natom=particles%n_els
     ALLOCATE(at_row(natom),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(rho_ij(natom,natom),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(d_ij(natom,natom),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(r_ij(natom,natom,3),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(fixed(3,natom),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     fixed=1.0_dp
     CALL fix_atom_control( force_env, error, fixed)
     DO i=1,3
        CALL mp_min(fixed(i,:),hess_mat%matrix_struct%para_env%group)
     END DO
     rho_ij=0
      !XXXX insert proper rows !XXX
     at_row=3
     DO i=1,natom
         CALL get_atomic_kind(atomic_kind=particles%els(i)%atomic_kind,z=z)
         IF(z.LE.10)at_row(i)=2
         IF(z.LE.2)at_row(i)=1
     END DO
    DO i=2,natom
        iat_row=at_row(i)
        DO j=1,i-1
           jat_row=at_row(j)
           !pbc for a distance vector
           r_ij(j,i,:)=pbc(particles%els(i)%r,particles%els(j)%r,cell)
           r_ij(i,j,:)=-r_ij(j,i,:)
           d_ij(j,i)=SQRT(DOT_PRODUCT(r_ij(j,i,:),r_ij(j,i,:)))
           d_ij(i,j)=d_ij(j,i)
           rho_ij(j,i)=EXP(alpha(jat_row,iat_row)*(r0(jat_row,iat_row)**2-d_ij(j,i)**2))
           rho_ij(i,j)=rho_ij(j,i)
        END DO
     END DO
     DO i=1,ncol_local
        iglobal=col_indices(i)
        iind=MOD(iglobal-1,3)+1
        iat_col=(iglobal+2)/3
        IF(iat_col.gt.natom)CYCLE
        DO j=1,nrow_local
           jglobal=row_indices(j)
           jind=MOD(jglobal-1,3)+1
           iat_row=(jglobal+2)/3 
           IF(iat_row.gt.natom)CYCLE
           IF(iat_row.ne.iat_col)THEN
               IF(d_ij(iat_row,iat_col).lt.6.0_dp)&
               local_data(j,i)=local_data(j,i)+&
                     angle_second_deriv(r_ij,d_ij,rho_ij,iind,jind,iat_col,iat_row,natom)
           ELSE
               local_data(j,i)=local_data(j,i)+&
                     angle_second_deriv(r_ij,d_ij,rho_ij,iind,jind,iat_col,iat_row,natom)
           END IF
           IF(iat_col.NE.iat_row)THEN
               IF(d_ij(iat_row,iat_col).lt.6.0_dp) &
               local_data(j,i)=local_data(j,i)-&
                               dist_second_deriv(r_ij(iat_col,iat_row,:),&
                               iind,jind,d_ij(iat_row,iat_col),rho_ij(iat_row,iat_col))
           ELSE
              DO k=1,natom
                 IF(k==iat_col)CYCLE
                 IF(d_ij(iat_row,k).lt.6.0_dp) &
                 local_data(j,i)=local_data(j,i)+&
                               dist_second_deriv(r_ij(iat_col,k,:),&
                               iind,jind,d_ij(iat_row,k),rho_ij(iat_row,k))
              END DO
           END IF
           IF(fixed(jind,iat_row).lt.0.5_dp.OR.fixed(iind,iat_col).lt.0.5_dp)THEN
              local_data(j,i)=0.0_dp
              IF(jind==iind.AND.iat_row==iat_col)local_data(j,i)=1.0_dp
           END IF 
        END DO
     END DO
     DEALLOCATE(fixed,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE(rho_ij,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE(d_ij,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE(r_ij,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE(at_row,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE  construct_initial_hess

! *****************************************************************************
!> \brief ...
!> \param r1 ...
!> \param i ...
!> \param j ...
!> \param d ...
!> \param rho ...
!> \retval deriv ...
! *****************************************************************************
  FUNCTION dist_second_deriv(r1,i,j,d,rho)  RESULT(deriv)
    REAL(KIND=dp), DIMENSION(3)              :: r1
    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: d, rho, deriv

      deriv=0.45_dp*rho*(r1(i)*r1(j))/d**2
  END FUNCTION

! *****************************************************************************
!> \brief ...
!> \param r_ij ...
!> \param d_ij ...
!> \param rho_ij ...
!> \param idir ...
!> \param jdir ...
!> \param iat_der ...
!> \param jat_der ...
!> \param natom ...
!> \retval deriv ...
! *****************************************************************************
  FUNCTION angle_second_deriv(r_ij,d_ij,rho_ij,idir,jdir,iat_der,jat_der,natom) RESULT(deriv)
    REAL(KIND=dp), DIMENSION(:, :, :)        :: r_ij
    REAL(KIND=dp), DIMENSION(:, :)           :: d_ij, rho_ij
    INTEGER                                  :: idir, jdir, iat_der, jat_der, &
                                                natom
    REAL(KIND=dp)                            :: deriv

    INTEGER                                  :: i, iat, idr, j, jat, jdr
    REAL(KIND=dp) :: d12, d23, d31, D_mat(3,2), denom1, denom2, denom3, ka1, &
      ka2, ka3, rho12, rho23, rho31, rsst1, rsst2, rsst3
    REAL(KIND=dp), DIMENSION(3)              :: r12, r23, r31

     deriv=0._dp
     IF(iat_der==jat_der)THEN
        DO i=1,natom-1
           IF(rho_ij(iat_der,i).lt.0.00001)CYCLE
           DO j=i+1,natom
              IF(rho_ij(iat_der,j).lt.0.00001)CYCLE
              IF(i==iat_der.OR.j==iat_der)CYCLE
              IF(iat_der.lt.i.OR.iat_der.gt.j)THEN
                 r12=r_ij(iat_der,i,:); r23=r_ij(i,j,:); r31=r_ij(j,iat_der,:)
                 d12=d_ij(iat_der,i); d23=d_ij(i,j); d31=d_ij(j,iat_der)
                 rho12=rho_ij(iat_der,i); rho23=rho_ij(i,j); rho31=rho_ij(j,iat_der)
              ELSE
                 r12=r_ij(iat_der,j,:); r23=r_ij(j,i,:); r31=r_ij(i,iat_der,:)
                 d12=d_ij(iat_der,j); d23=d_ij(j,i); d31=d_ij(i,iat_der)
                 rho12=rho_ij(iat_der,j); rho23=rho_ij(j,i); rho31=rho_ij(i,iat_der)
              END IF
              ka1=0.15_dp*rho12*rho23; ka2=0.15_dp*rho23*rho31; ka3=0.15_dp*rho31*rho12
              rsst1=DOT_PRODUCT(r12,r23);rsst2=DOT_PRODUCT(r23,r31);rsst3=DOT_PRODUCT(r31,r12)
              denom1=1.0_dp-rsst1**2/(d12**2*d23**2); denom2=1.0_dp-rsst2**2/(d23**2*d31**2)
              denom3=1.0_dp-rsst3**2/(d31**2*d12**2)
              denom1=SIGN(1.0_dp,denom1)*MAX(ABS(denom1),0.01_dp)
              denom2=SIGN(1.0_dp,denom2)*MAX(ABS(denom2),0.01_dp)
              denom3=SIGN(1.0_dp,denom3)*MAX(ABS(denom3),0.01_dp)
              D_mat(1,1)=r23(idir)/(d12*d23)-rsst1*r12(idir)/(d12**3*d23)
              D_mat(1,2)=r23(jdir)/(d12*d23)-rsst1*r12(jdir)/(d12**3*d23)
              D_mat(2,1)=-r23(idir)/(d23*d31)+rsst2*r31(idir)/(d23*d31**3)
              D_mat(2,2)=-r23(jdir)/(d23*d31)+rsst2*r31(jdir)/(d23*d31**3)
              D_mat(3,1)=(r31(idir)-r12(idir))/(d31*d12)+rsst3*r31(idir)/(d31**3*d12)-&
                                                         rsst3*r12(idir)/(d31*d12**3)
              D_mat(3,2)=(r31(jdir)-r12(jdir))/(d31*d12)+rsst3*r31(jdir)/(d31**3*d12)-&
                                                         rsst3*r12(jdir)/(d31*d12**3)
              IF(ABS(denom1).le.0.011_dp)D_mat(1,1)=0.0_dp
              IF(ABS(denom2).le.0.011_dp)D_mat(2,1)=0.0_dp
              IF(ABS(denom3).le.0.011_dp)D_mat(3,1)=0.0_dp
              deriv=deriv+ka1*D_mat(1,1)*D_mat(1,2)/denom1+&
                    ka2*D_mat(2,1)*D_mat(2,2)/denom2+&
                    ka3*D_mat(3,1)*D_mat(3,2)/denom3

           END DO
        END DO
     ELSE
        DO i=1,natom
           IF(i==iat_der.OR.i==jat_der)CYCLE
           IF(jat_der.lt.iat_der)THEN
              iat=jat_der; jat=iat_der; idr=jdir; jdr=idir
           ELSE
              iat=iat_der; jat=jat_der; idr=idir; jdr=jdir
           END IF
           IF(jat.lt.i.OR.iat.gt.i)THEN
              r12=r_ij(iat,jat,:); r23=r_ij(jat,i,:); r31=r_ij(i,iat,:)
              d12=d_ij(iat,jat); d23=d_ij(jat,i); d31=d_ij(i,iat)
              rho12=rho_ij(iat,jat); rho23=rho_ij(jat,i); rho31=rho_ij(i,iat)
           ELSE
              r12=r_ij(iat,i,:); r23=r_ij(i,jat,:); r31=r_ij(jat,iat,:)
              d12=d_ij(iat,i); d23=d_ij(i,jat); d31=d_ij(jat,iat)
              rho12=rho_ij(iat,i); rho23=rho_ij(i,jat); rho31=rho_ij(jat,iat)
           END IF
           ka1=0.15_dp*rho12*rho23; ka2=0.15_dp*rho23*rho31; ka3=0.15_dp*rho31*rho12
           rsst1=DOT_PRODUCT(r12,r23);rsst2=DOT_PRODUCT(r23,r31);rsst3=DOT_PRODUCT(r31,r12)
           denom1=1.0_dp-rsst1**2/(d12**2*d23**2); denom2=1.0_dp-rsst2**2/(d23**2*d31**2)
           denom3=1.0_dp-rsst3**2/(d31**2*d12**2)
           denom1=SIGN(1.0_dp,denom1)*MAX(ABS(denom1),0.01_dp)
           denom2=SIGN(1.0_dp,denom2)*MAX(ABS(denom2),0.01_dp)
           denom3=SIGN(1.0_dp,denom3)*MAX(ABS(denom3),0.01_dp)
           D_mat(1,1)=r23(idr)/(d12*d23)-rsst1*r12(idr)/(d12**3*d23)
           D_mat(2,1)=-r23(idr)/(d23*d31)+rsst2*r31(idr)/(d23*d31**3)
           D_mat(3,1)=(r31(idr)-r12(idr))/(d31*d12)+rsst3*r31(idr)/(d31**3*d12)-&
                                                         rsst3*r12(idr)/(d31*d12**3)
           IF(jat.lt.i.OR.iat.gt.i)THEN
              D_mat(1,2)=(r12(jdr)-r23(jdr))/(d12*d23)+rsst1*r12(jdr)/(d12**3*d23)-&
                                                       rsst1*r23(jdr)/(d12*d23**3)
              D_mat(2,2)=r31(jdr)/(d23*d31)-rsst2*r23(jdr)/(d23**3*d31)
              D_mat(3,2)=-r31(jdr)/(d31*d12)+rsst3*r12(jdr)/(d31*d12**3)
           ELSE
              D_mat(1,2)=-r12(jdr)/(d12*d23)+rsst1*r23(jdr)/(d12*d23**3)
              D_mat(2,2)=(r23(jdr)-r31(jdr))/(d23*d31)+rsst2*r23(jdr)/(d23**3*d31)-&
                                                       rsst2*r31(jdr)/(d23*d31**3)
              D_mat(3,2)=r12(jdr)/(d31*d12)-rsst3*r31(jdr)/(d31**3*d12)
           END IF
           IF(ABS(denom1).le.0.011_dp)D_mat(1,1)=0.0_dp
           IF(ABS(denom2).le.0.011_dp)D_mat(2,1)=0.0_dp
           IF(ABS(denom3).le.0.011_dp)D_mat(3,1)=0.0_dp

           deriv=deriv+ka1*D_mat(1,1)*D_mat(1,2)/denom1+&
                 ka2*D_mat(2,1)*D_mat(2,2)/denom2+&
                 ka3*D_mat(3,1)*D_mat(3,2)/denom3
         END DO
      END IF
      deriv=0.25_dp*deriv

  END FUNCTION angle_second_deriv

END MODULE bfgs_optimizer
