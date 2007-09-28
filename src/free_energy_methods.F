!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/free_energy_methods [1.0] *
!!
!!   NAME
!!     free_energy_methods
!!
!!   FUNCTION
!!     Methods to perform free energy and free energy derivatives calculations
!!
!!   AUTHOR
!!     Teodoro Laino (01.2007) [tlaino]
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE free_energy_methods
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE fparser,                         ONLY: evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE free_energy_types,               ONLY: free_energy_type,&
                                             ui_var_type
  USE input_constants,                 ONLY: do_fe_ac,&
                                             do_fe_ui
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'free_energy_methods'
  PUBLIC :: free_energy_evaluate

CONTAINS

!!****f* free_energy_methods/free_energy_evaluate *
!!
!!   NAME
!!     free_energy_evaluate
!!
!!   FUNCTION
!!     Main driver for free energy calculations
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino (01.2007) [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE free_energy_evaluate(md_env,converged,motion_section,error)
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(OUT)                     :: converged
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'free_energy_evaluate', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par
    INTEGER                                  :: handle, ic, icolvar, &
                                                nforce_eval, output_unit
    INTEGER, POINTER                         :: istep
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: beta, dx, lerr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_p_type), POINTER          :: force_env_p( : )
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(ui_var_type), POINTER               :: cv

    NULLIFY(force_env_p, force_env, istep, subsys, cv, root_section, simpar)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,"I","",handle)
    failure   = .FALSE.
    converged = .FALSE.
    ! Metadynamics is also a free energy calculation but is handled in a different
    ! module.
    ! In this routine specifically biased MD are handled.
    CALL get_md_env(md_env, force_env_p=force_env_p, fe_env=fe_env, simpar=simpar,&
         itimes=istep, error=error)
    force_env => force_env_p(1)%force_env
    IF (.NOT.ASSOCIATED(force_env%meta_env).AND.ASSOCIATED(fe_env)) THEN
       SELECT CASE(fe_env%type)
       CASE(do_fe_ui)
          ! Umbrella Integration..
          CALL force_env_get(force_env, subsys=subsys, error=error)
          DO ic=1,fe_env%ncolvar
             cv => fe_env%uivar(ic)
             icolvar = cv%icolvar
             CALL colvar_eval_glob_f(icolvar,force_env,error=error)
             CALL reallocate(cv%ss,1,istep)
             cv%ss(istep) = subsys%colvar_p(icolvar)%colvar%ss
          END DO
          ! Check Free Energy derivative convergence
          output_unit = cp_print_key_unit_nr(logger,motion_section,&
               "FREE_ENERGY%FREE_ENERGY_INFO",&
               extension=".FreeEnergyLog",log_filename=.FALSE.,error=error)
          CALL ui_check_convergence(fe_env,converged,istep,output_unit,error)
          CALL cp_print_key_finished_output(output_unit, logger,&
               motion_section,"FREE_ENERGY%FREE_ENERGY_INFO",error=error)
       CASE(do_fe_ac)
          CALL initf(2)
          ! Alchemical Changes
          CALL cp_assert(ASSOCIATED(force_env%mixed_env),cp_fatal_level,-300,routineP,&
               'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
               ' Free Energy calculations require the definition of a mixed env!',error)
          my_par => force_env%mixed_env%par
          my_val => force_env%mixed_env%val
          dx     =  force_env%mixed_env%dx
          lerr   =  force_env%mixed_env%lerr
          root_section => force_env%root_section
          coupling_function = force_env%mixed_env%coupling_function
          beta   =  1/simpar%temp_ext
          CALL parsef(1,TRIM(coupling_function),my_par)
          nforce_eval = SIZE(force_env%sub_force_env)
          CALL dump_ac_info(my_val, my_par, dx, lerr, root_section, nforce_eval,&
               fe_env%covmx, istep, beta, error)
          CALL finalizef()
       CASE DEFAULT
          ! Do Nothing
       END SELECT
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE free_energy_evaluate

!!****f* free_energy_methods/ui_check_convergence *
!!
!!   NAME
!!     ui_check_convergence
!!
!!   FUNCTION
!!     Convergence criteria for umbrella integration free energy method
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino (01.2007) [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE ui_check_convergence(fe_env,converged,istep,output_unit,error)        
    TYPE(free_energy_type), POINTER          :: fe_env
    LOGICAL, INTENT(OUT)                     :: converged
    INTEGER, POINTER                         :: istep
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ui_check_convergence', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ic, jc, ncolvar, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps_conv, max_error, &
                                                max_error_avg, max_error_cov
    REAL(KIND=dp), DIMENSION(:), POINTER     :: avgs, wrk
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: covmx, covmx_old
    TYPE(ui_var_type), POINTER               :: cv

    failure   = .FALSE.
    converged =  .FALSE.
    eps_conv  =  fe_env%eps_conv
    covmx     => fe_env%covmx 
    ncolvar   =  fe_env%ncolvar
    CALL timeset(routineN,"I","",handle)
    ! Allocate and store temporarily the old covariance matrix
    ALLOCATE(covmx_old(ncolvar,ncolvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    covmx_old = covmx
    ! Compute new covariance matrix
    ! First compute the average for all collective variables
    ALLOCATE(wrk(ncolvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(avgs(ncolvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ic = 1, fe_env%ncolvar
      avgs(ic) = SUM(fe_env%uivar(ic)%ss)/REAL(istep,KIND=dp)
    END DO 
    covmx = 0.0_dp
    ! Compute the covariance matrix
    DO i =  1, istep
      ! Compute the (X-Xavg)_i 
      DO ic=1,ncolvar
         cv => fe_env%uivar(ic)
         wrk(ic) = cv%ss(i) - avgs(ic)
      END DO
      ! Compute the contribution to the covariance matrix for step i
      DO ic = 1, ncolvar
        DO jc = 1, ic
          covmx(jc,ic) = covmx(jc,ic) + wrk(ic)*wrk(jc)
        END DO
      END DO
    END DO
    ! Normalized the covariance matrix
    covmx = covmx/REAL(istep,KIND=dp)
    ! Symmetrize the covariance matrix
    DO ic = 1, ncolvar
       DO jc = ic+1,ncolvar
          covmx(jc,ic) = covmx(ic,jc)
       END DO
    END DO
    IF (istep>1) THEN
       ! Convergence of the averages
       DO ic = 1, ncolvar
          wrk(ic) = SUM(fe_env%uivar(ic)%ss(:istep-1))/REAL(istep-1,KIND=dp)
       END DO
       max_error_avg = MAXVAL(ABS(avgs-wrk)/avgs)
       ! Convergence of the covariant matrix
       max_error_cov = MAXVAL(ABS(covmx-covmx_old)/covmx)
       ! Total Error and tolerance check
       max_error = MAX(max_error_cov, max_error_avg)
       IF (max_error<=eps_conv) converged = .TRUE.
       ! Print some info results
       IF (output_unit>0) THEN
          WRITE(output_unit,'(80("*"))')
          WRITE(output_unit,'(T30,"FREE ENERGY CALCULATION")')
          WRITE(output_unit,'(T7,"AVERAGES",T50,"COVARIANCE MATRIX")')
          DO ic = 1, ncolvar
             WRITE(output_unit,'(T5,E12.6,T29,"|",4(E12.6,1X))')avgs(ic), covmx(:,ic)
          END DO
          WRITE(output_unit,'(/,T2,"MAX_ERROR AVG:",E12.6,T54,"MAX_ERROR COV:",E12.6)')&
               max_error_avg,max_error_cov
          WRITE(output_unit,'(80("*"),/)')
       END IF
    END IF
    DEALLOCATE(covmx_old,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(avgs,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE ui_check_convergence

!!****f* free_energy_methods/dump_der_par *
!!
!!   NAME
!!     dump_der_par
!!
!!   FUNCTION
!!     Dumps information when performing an alchemical change run
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich [tlaino] - 05.2007
!!
!!*** **********************************************************************
SUBROUTINE dump_ac_info(my_val, my_par, dx, lerr, root_section, nforce_eval, cum_res,&
     istep, beta, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par
    REAL(KIND=dp), INTENT(IN)                :: dx, lerr
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: nforce_eval
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cum_res
    INTEGER, POINTER                         :: istep
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dump_ac_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, par, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: data_function
    INTEGER                                  :: i, iforce_eval, ipar, isize, &
                                                iw, j, NEquilStep, output_unit
    LOGICAL                                  :: failure, ionode
    REAL(KIND=dp) :: avg_BP, avg_DET, avg_DUE, d_ene_w, dedf, ene_w, err, &
      Err_DET, Err_DUE, std_DET, std_DUE, tmp, tmp2, wfac
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: fe_section

    NULLIFY(data_function)
    failure = .FALSE.
    output_unit = -1
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)
    fe_section => section_vals_get_subs_vals(root_section,&
         "MOTION%FREE_ENERGY%ALCHEMICAL_CHANGE",error=error)
    CALL section_vals_val_get(fe_section,"PARAMETER",c_val=par,error=error)
    DO i = 1, SIZE(my_par)
       IF (my_par(i)==par) EXIT
    END DO
    CPPrecondition(i<=SIZE(my_par),cp_failure_level,routineP,error,failure)
    ipar = i
    dedf = evalfd(1,ipar,my_val,dx,err)
    IF ((ABS(err)>lerr).AND.(output_unit>0)) THEN
       WRITE(this_error,"(A,G12.6,A)")"(",err,")"
       WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
       CALL compress(this_error,.TRUE.)
       CALL compress(def_error,.TRUE.)
       CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
            'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
            ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
            TRIM(def_error)//' .',error)
    END IF
    
    ! We must print now the energy of the biased system, the weigthing energy 
    ! and the derivative w.r.t.the coupling parameter of the biased energy
    ! Retrieve the expression of the weighting function:
    CALL section_vals_val_get(fe_section,"WEIGHTING_FUNCTION",&
         c_vals=data_function,error=error)
    ! Clean the function structure
    coupling_function = ""
    DO j = 1, SIZE(data_function)
       coupling_function(LEN_TRIM(coupling_function)+1:) = TRIM(data_function(j))
    END DO
    CALL parsef(2,TRIM(coupling_function),my_par)       
    ene_w   = evalf(2,my_val)
    d_ene_w = evalfd(2,ipar,my_val,dx,err)
    IF ((ABS(err)>lerr).AND.(output_unit>0)) THEN
       WRITE(this_error,"(A,G12.6,A)")"(",err,")"
       WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
       CALL compress(this_error,.TRUE.)
       CALL compress(def_error,.TRUE.)
       CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
            'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
            ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
            TRIM(def_error)//' .',error)
    END IF
    CALL section_vals_val_get(fe_section,"NEQUIL_STEPS",i_val=NEquilStep,error=error)
    ! Store results
    IF (istep>NEquilStep) THEN
       isize = SIZE(cum_res,2)+1
       CALL reallocate(cum_res,1,3,1,isize)
       cum_res(1,isize) = dedf
       cum_res(2,isize) = dedf-d_ene_w
       cum_res(3,isize) = ene_w
       ! Compute derivative of biased and total energy
       ! Total Free Energy
       avg_DET = SUM(cum_res(1,1:isize))/REAL(isize,KIND=dp)
       std_DET = SUM(cum_res(1,1:isize)**2)/REAL(isize,KIND=dp)
       ! Unbiased Free Energy
       avg_BP  = SUM(cum_res(3,1:isize))/REAL(isize,KIND=dp)
       wfac = 0.0_dp
       DO j = 1, isize
          wfac = wfac + EXP(beta*(cum_res(3,j)-avg_BP))
       END DO
       avg_DUE = 0.0_dp
       std_DUE = 0.0_dp
       DO j = 1, isize
          tmp = cum_res(2,j) 
          tmp2= EXP(beta*(cum_res(3,j)-avg_BP))/wfac
          avg_DUE = avg_DUE + tmp    * tmp2
          std_DUE = std_DUE + tmp**2 * tmp2
       END DO
       IF (isize>1) THEN
          Err_DUE = SQRT(std_DUE - avg_DUE**2)/SQRT(REAL(isize-1,KIND=dp))
          Err_DET = SQRT(std_DET - avg_DET**2)/SQRT(REAL(isize-1,KIND=dp))
       END IF
       ! Print info
       iw = cp_print_key_unit_nr(logger,root_section,"MOTION%FREE_ENERGY%FREE_ENERGY_INFO",&
            extension=".free_energy",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T2,79("-"),T37," oOo ")')
          DO iforce_eval = 1, nforce_eval
             WRITE(iw,'(T2,"ALCHEMICAL CHANGE| FORCE_EVAL Nr.",I5,T48,"ENERGY (Hartree)= ",F15.9)')&
                  iforce_eval, my_val(iforce_eval)
          END DO
          WRITE(iw,'(T2,"ALCHEMICAL CHANGE| DERIVATIVE OF TOTAL ENERGY  [ PARAMETER (",A,") ]",T66,F15.9)')&
               TRIM(par),dedf
          WRITE(iw,'(T2,"ALCHEMICAL CHANGE| DERIVATIVE OF BIASED ENERGY [ PARAMETER (",A,") ]",T66,F15.9)')&
               TRIM(par),dedf-d_ene_w
          WRITE(iw,'(T2,"ALCHEMICAL CHANGE| BIASING UMBRELLA POTENTIAL  ",T66,F15.9)')&
               ene_w
          
          IF (isize>1) THEN
             WRITE(iw,'(/,T2,"ALCHEMICAL CHANGE| DERIVATIVE TOTAL FREE ENERGY  ",T50,F15.9,1X,"+/-",1X,F11.9)')&
                  avg_DET,Err_DET
             WRITE(iw,'(T2,"ALCHEMICAL CHANGE| DERIVATIVE UNBIASED FREE ENERGY  ",T50,F15.9,1X,"+/-",1X,F11.9)')&
                  avg_DUE,Err_DUE
          ELSE
             WRITE(iw,'(/,T2,"ALCHEMICAL CHANGE| DERIVATIVE TOTAL FREE ENERGY  ",T50,F15.9,1X,"+/-",1X,T76,A)')&
                  avg_DET,"UNDEF"
             WRITE(iw,'(T2,"ALCHEMICAL CHANGE| DERIVATIVE UNBIASED FREE ENERGY  ",T50,F15.9,1X,"+/-",1X,T76,A)')&
                  avg_DUE,"UNDEF"
          END IF
          WRITE(iw,'(T2,79("-"))')
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,root_section,&
         "MOTION%FREE_ENERGY%FREE_ENERGY_INFO",error=error)

  END SUBROUTINE dump_ac_info

END MODULE free_energy_methods


