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
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE free_energy_types,               ONLY: free_energy_type,&
                                             ui_var_type
  USE input_constants,                 ONLY: do_fe_ui
  USE kinds,                           ONLY: dp
  USE md_environment_types,            ONLY: md_environment_type, get_md_env
  USE memory_utilities,                ONLY: reallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE input_section_types, only: section_vals_type
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

    INTEGER                                  :: handle, ic, isubsys, icolvar
    INTEGER, POINTER                         :: istep
    LOGICAL                                  :: cv_converged, failure
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(force_env_p_type), POINTER          :: force_env_p( : )
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(ui_var_type), POINTER               :: cv
    INTEGER :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(force_env_p, force_env, istep, subsys, cv)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,"I","",handle)
    failure   = .FALSE.
    converged = .FALSE.
    ! Metadynamics is also a free energy calculation but is handled in a different
    ! module.
    ! In this routine specifically biased MD are handled.
    CALL get_md_env(md_env, force_env_p=force_env_p, fe_env=fe_env, itimes=istep)
    force_env => force_env_p(1)%force_env
    IF (.NOT.ASSOCIATED(force_env%meta_env).AND.ASSOCIATED(fe_env)) THEN
       SELECT CASE(fe_env%type)
       CASE(do_fe_ui)
          CALL force_env_get(force_env, subsys=subsys, error=error)
          DO ic=1,fe_env%ncolvar
             cv => fe_env%uivar(ic)
             isubsys = cv%isubsys
             icolvar = cv%icolvar
             CALL colvar_eval_glob_f(isubsys,icolvar,force_env,error=error)
             CALL reallocate(cv%ss,1,istep)
             cv%ss(istep) = subsys(isubsys)%subsys%colvar_p(icolvar)%colvar%ss
          END DO
          ! Check Free Energy derivative convergence
          output_unit = cp_print_key_unit_nr(logger,motion_section,&
               "FREE_ENERGY%FREE_ENERGY_INFO",&
               extension=".FreeEnergyLog",log_filename=.FALSE.,error=error)
          CALL ui_check_convergence(fe_env,converged,istep,output_unit,error)
          CALL cp_print_key_finished_output(output_unit, logger,&
               motion_section,"FREE_ENERGY%FREE_ENERGY_INFO",error=error)    
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
    INTEGER, POINTER :: istep
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(free_energy_type), POINTER          :: fe_env
    LOGICAL, INTENT(OUT)                     :: converged
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ui_check_convergence', &
      routineP = moduleN//':'//routineN

    LOGICAL :: failure
    REAL(KIND=dp)                :: eps_conv, max_error_cov, max_error_avg,&
         max_error
    TYPE(ui_var_type), POINTER               :: cv
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: covmx
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: covmx_old
    REAL(KIND=dp), DIMENSION(:), POINTER     :: avgs, wrk
    INTEGER :: stat, i, ic, jc, ncolvar, handle

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

END MODULE free_energy_methods


