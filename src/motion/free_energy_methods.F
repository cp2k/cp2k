!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods to perform free energy and free energy derivatives calculations
!> \author Teodoro Laino (01.2007) [tlaino]
! *****************************************************************************
MODULE free_energy_methods
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
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
  USE mathlib,                         ONLY: diamat_all
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE simpar_types,                    ONLY: simpar_type
  USE statistical_methods,             ONLY: k_test,&
                                             min_sample_size,&
                                             sw_test,&
                                             vn_test
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'free_energy_methods'
  PUBLIC :: free_energy_evaluate

CONTAINS

! *****************************************************************************
!> \brief Main driver for free energy calculations
!>      In this routine we handle specifically biased MD.
!> \param md_env ...
!> \param converged ...
!> \param fe_section ...
!> \param error ...
!> \par History
!>      Teodoro Laino (01.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE free_energy_evaluate(md_env,converged,fe_section,error)
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(OUT)                     :: converged
    TYPE(section_vals_type), POINTER         :: fe_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'free_energy_evaluate', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par
    INTEGER                                  :: handle, ic, icolvar, &
                                                nforce_eval, output_unit, &
                                                stat_sign_points
    INTEGER, POINTER                         :: istep
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: beta, dx, lerr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(ui_var_type), POINTER               :: cv

    NULLIFY(force_env, istep, subsys, cv, simpar)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)
    failure   = .FALSE.
    converged = .FALSE.
    CALL get_md_env(md_env, force_env=force_env, fe_env=fe_env, simpar=simpar,&
         itimes=istep, error=error)
    ! Metadynamics is also a free energy calculation but is handled in a different
    ! module.
    IF (.NOT.ASSOCIATED(force_env%meta_env).AND.ASSOCIATED(fe_env)) THEN
       SELECT CASE(fe_env%type)
       CASE(do_fe_ui)
          ! Umbrella Integration..
          CALL force_env_get(force_env, subsys=subsys, error=error)
          fe_env%nr_points = fe_env%nr_points + 1
          DO ic=1,fe_env%ncolvar
             cv => fe_env%uivar(ic)
             icolvar = cv%icolvar
             CALL colvar_eval_glob_f(icolvar,force_env,error=error)
             CALL reallocate(cv%ss,1,fe_env%nr_points)
             cv%ss(fe_env%nr_points) = subsys%colvar_p(icolvar)%colvar%ss
             WRITE(*,*)"COLVAR::",cv%ss(fe_env%nr_points)
          END DO
          stat_sign_points = fe_env%nr_points - fe_env%nr_rejected
          WRITE(*,*)fe_env%nr_points, stat_sign_points
          ! Start statistical analysis when enough CG data points have been collected
          IF ( (fe_env%conv_par%cg_width*fe_env%conv_par%cg_points <= stat_sign_points).AND.&
               (MOD(stat_sign_points,fe_env%conv_par%cg_width)==0) ) THEN
             output_unit = cp_print_key_unit_nr(logger,fe_section,"FREE_ENERGY_INFO",&
                  extension=".FreeEnergyLog",log_filename=.FALSE.,error=error)
             CALL print_fe_prolog(output_unit,error)
             ! Trend test..  recomputes the number of statistically significant points..
             CALL ui_check_trend(fe_env, fe_env%conv_par%test_k, stat_sign_points, output_unit, error)
             stat_sign_points = fe_env%nr_points - fe_env%nr_rejected
             ! Normality and serial correlation tests..
             IF (fe_env%conv_par%cg_width*fe_env%conv_par%cg_points <= stat_sign_points.AND.&
                  fe_env%conv_par%test_k) THEN
                ! Statistical tests
                CALL ui_check_convergence(fe_env,converged,stat_sign_points,output_unit,error)
             END IF
             CALL print_fe_epilog(output_unit,error)
             CALL cp_print_key_finished_output(output_unit, logger,fe_section,"FREE_ENERGY_INFO",error=error)
          END IF
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
          coupling_function = force_env%mixed_env%coupling_function
          beta   =  1/simpar%temp_ext
          CALL parsef(1,TRIM(coupling_function),my_par)
          nforce_eval = SIZE(force_env%sub_force_env)
          CALL dump_ac_info(my_val, my_par, dx, lerr, fe_section, nforce_eval,&
               fe_env%covmx, istep, beta, error)
          CALL finalizef()
       CASE DEFAULT
          ! Do Nothing
       END SELECT
    END IF
    CALL timestop(handle)

  END SUBROUTINE free_energy_evaluate

! *****************************************************************************
!> \brief Print prolog of free energy output section
!> \param output_unit ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE print_fe_prolog(output_unit,error)
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_fe_prolog', &
      routineP = moduleN//':'//routineN

    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,79("*"))')
       WRITE(output_unit,'(T30,"FREE ENERGY CALCULATION",/)')
    END IF
  END SUBROUTINE print_fe_prolog

! *****************************************************************************
!> \brief Print epilog of free energy output section
!> \param output_unit ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE print_fe_epilog(output_unit,error)
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_fe_epilog', &
      routineP = moduleN//':'//routineN

    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,79("*"),/)')
    END IF
  END SUBROUTINE print_fe_epilog

! *****************************************************************************
!> \brief Test for trend in coarse grained data set
!> \param fe_env ...
!> \param trend_free ...
!> \param nr_points ...
!> \param output_unit ...
!> \param error ...
!> \par History
!>      Teodoro Laino (01.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE ui_check_trend(fe_env, trend_free, nr_points, output_unit, error)
    TYPE(free_energy_type), POINTER          :: fe_env
    LOGICAL, INTENT(OUT)                     :: trend_free
    INTEGER, INTENT(IN)                      :: nr_points, output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ui_check_trend', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, &
                                                my_reject, ncolvar, &
                                                ng_points, rejected_points
    LOGICAL                                  :: failure, test_avg, test_std
    REAL(KIND=dp)                            :: prob, tau, z
    REAL(KIND=dp), DIMENSION(:), POINTER     :: wrk

    CALL timeset(routineN,handle)
    failure    = .FALSE.
    trend_free = .FALSE.
    test_avg   = .TRUE.
    test_std   = .TRUE.
    ncolvar    = fe_env%ncolvar
    ! Number of coarse grained points
    WRITE(*,*)nr_points, fe_env%conv_par%cg_width
    ng_points = nr_points/fe_env%conv_par%cg_width
    my_reject = 0
    ! Allocate storage
    CALL create_tmp_data(fe_env,wrk,ng_points,ncolvar,error)
    ! Compute the Coarse Grained data set using a reverse cumulative strategy
    CALL create_csg_data(fe_env,ng_points,error)
    ! Test on coarse grained average
    DO j = 1, ncolvar
       ii = 1
       DO i = ng_points, 1, -1
          wrk(ii) = fe_env%cg_data(i)%avg(j)
          ii = ii + 1
       END DO
       DO i = my_reject+1, ng_points
          IF ((ng_points-my_reject).LT.min_sample_size) THEN
             my_reject = MAX(0,my_reject-1)
             test_avg = .FALSE.
             EXIT
          END IF
          CALL k_test(wrk,my_reject+1,ng_points,tau,z,prob,error)
          PRINT *,prob,fe_env%conv_par%k_conf_lm
          IF (prob<fe_env%conv_par%k_conf_lm) EXIT
          my_reject = my_reject + 1
       END DO
       my_reject = MIN(ng_points,my_reject)
    END DO
    rejected_points    = my_reject * fe_env%conv_par%cg_width
    ! Print some info
    IF (output_unit>0) THEN
       WRITE(output_unit,*)"Kendall trend test (Average)",test_avg,&
            "number of points rejected:",rejected_points+fe_env%nr_rejected
       WRITE(output_unit,*)"Reject Nr.",my_reject," coarse grained points testing average"
    END IF
    ! Test on coarse grained covariance matrix
    DO j = 1, ncolvar
       DO k = j, ncolvar
          ii = 1
          DO i = ng_points, 1, -1
             wrk(ii) = fe_env%cg_data(i)%var(j,k)
             ii = ii + 1
          END DO
          DO i = my_reject+1, ng_points
             IF ((ng_points-my_reject).LT.min_sample_size) THEN
                my_reject = MAX(0,my_reject-1)
                test_std  = .FALSE.
                EXIT
             END IF
             CALL k_test(wrk,my_reject+1,ng_points,tau,z,prob,error)
             PRINT *,prob,fe_env%conv_par%k_conf_lm
             IF (prob<fe_env%conv_par%k_conf_lm) EXIT
             my_reject = my_reject + 1
          END DO
          my_reject = MIN(ng_points,my_reject)
       END DO
    END DO
    rejected_points    = my_reject * fe_env%conv_par%cg_width
    fe_env%nr_rejected = fe_env%nr_rejected + rejected_points
    trend_free         = test_avg.AND.test_std
    ! Print some info
    IF (output_unit>0) THEN
       WRITE(output_unit,*)"Kendall trend test (Std. Dev.)",test_std,&
            "number of points rejected:",fe_env%nr_rejected
       WRITE(output_unit,*)"Reject Nr.",my_reject," coarse grained points testing standard dev."
       WRITE(output_unit,*)"Kendall test passed:",trend_free
    END IF
    ! Release storage
    CALL destroy_tmp_data(fe_env,wrk,ng_points,error)
    CALL timestop(handle)
  END SUBROUTINE ui_check_trend

! *****************************************************************************
!> \brief Creates temporary data structures
!> \param fe_env ...
!> \param wrk ...
!> \param ng_points ...
!> \param ncolvar ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE create_tmp_data(fe_env,wrk,ng_points,ncolvar,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: wrk
    INTEGER, INTENT(IN)                      :: ng_points, ncolvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_tmp_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ALLOCATE(fe_env%cg_data(ng_points),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, ng_points
          ALLOCATE(fe_env%cg_data(i)%avg(ncolvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fe_env%cg_data(i)%var(ncolvar,ncolvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       IF (PRESENT(wrk)) THEN
          ALLOCATE(wrk(ng_points),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE create_tmp_data

! *****************************************************************************
!> \brief Destroys temporary data structures
!> \param fe_env ...
!> \param wrk ...
!> \param ng_points ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE destroy_tmp_data(fe_env,wrk,ng_points,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: wrk
    INTEGER, INTENT(IN)                      :: ng_points
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_tmp_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       DO i = 1, ng_points
          DEALLOCATE(fe_env%cg_data(i)%avg,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(fe_env%cg_data(i)%var,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(fe_env%cg_data,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (PRESENT(wrk)) THEN
          DEALLOCATE(wrk,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE destroy_tmp_data

! *****************************************************************************
!> \brief Fills in temporary arrays with coarse grained data
!> \param fe_env ...
!> \param ng_points ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE create_csg_data(fe_env,ng_points,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    INTEGER, INTENT(IN)                      :: ng_points
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_csg_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iend, istart

    DO i = 1, ng_points
       istart = fe_env%nr_points -    (i)*fe_env%conv_par%cg_width + 1
       iend   = fe_env%nr_points -  (i-1)*fe_env%conv_par%cg_width
       WRITE(*,*)istart, iend
       CALL eval_cov_matrix(fe_env,cg_index=i,istart=istart,iend=iend,error=error)
    END DO

  END SUBROUTINE create_csg_data

! *****************************************************************************
!> \brief Checks Normality of the distribution and Serial Correlation of
!>      coarse grained data
!> \param fe_env ...
!> \param test_passed ...
!> \param nr_points ...
!> \param output_unit ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE ui_check_norm_sc(fe_env,test_passed,nr_points,output_unit,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    LOGICAL, INTENT(OUT)                     :: test_passed
    INTEGER, INTENT(IN)                      :: nr_points, output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ui_check_norm_sc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ng_points
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure     = .FALSE.
    test_passed = .FALSE.
    DO WHILE(fe_env%conv_par%cg_width<fe_env%conv_par%max_cg_width)
       ng_points = nr_points/fe_env%conv_par%cg_width
       PRINT *, ng_points
       IF (ng_points<min_sample_size) EXIT
       CALL ui_check_norm_sc_low(fe_env,nr_points,output_unit,error)
       test_passed = fe_env%conv_par%test_vn.AND.fe_env%conv_par%test_sw
       IF (test_passed) EXIT
       fe_env%conv_par%cg_width = fe_env%conv_par%cg_width + 1
       WRITE(*,*)"New coarse grained width:",fe_env%conv_par%cg_width
    END DO
    IF (fe_env%conv_par%cg_width==fe_env%conv_par%max_cg_width.AND.(.NOT.(test_passed))) THEN
       CALL ui_check_norm_sc_low(fe_env,nr_points,output_unit,error)
       test_passed = fe_env%conv_par%test_vn.AND.fe_env%conv_par%test_sw
    END IF
    CALL timestop(handle)
  END SUBROUTINE ui_check_norm_sc

! *****************************************************************************
!> \brief Checks Normality of the distribution and Serial Correlation of
!>      coarse grained data - Low Level routine
!> \param fe_env ...
!> \param nr_points ...
!> \param output_unit ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE ui_check_norm_sc_low(fe_env,nr_points,output_unit,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    INTEGER, INTENT(IN)                      :: nr_points, output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ui_check_norm_sc_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, ncolvar, &
                                                ng_points
    LOGICAL                                  :: avg_test_passed, failure, &
                                                sdv_test_passed
    REAL(KIND=dp)                            :: prob, pw, r, u, w
    REAL(KIND=dp), DIMENSION(:), POINTER     :: wrk

    CALL timeset(routineN,handle)
    failure     = .FALSE.
    ncolvar     = fe_env%ncolvar
    ! Compute the Coarse Grained data set using a reverse cumulative strategy
    fe_env%conv_par%test_sw = .FALSE.
    fe_env%conv_par%test_vn = .FALSE.
    ! Number of coarse grained points
    avg_test_passed = .TRUE.
    sdv_test_passed = .TRUE.
    ng_points = nr_points/fe_env%conv_par%cg_width
    CALL create_tmp_data(fe_env,wrk,ng_points,ncolvar,error)
    CALL create_csg_data(fe_env,ng_points,error)
    ! Testing Averages
    DO j = 1, ncolvar
       DO i = 1, ng_points
          wrk(i) = fe_env%cg_data(i)%avg(j)
       END DO
       ! Test of Shapiro - Wilks for normality
       !                 - Average
       CALL sw_test(wrk, ng_points, w, pw, error)
       PRINT *, 1.0_dp-pw,fe_env%conv_par%sw_conf_lm
       avg_test_passed = (1.0_dp-pw)<=fe_env%conv_par%sw_conf_lm
       fe_env%conv_par%test_sw = avg_test_passed
       IF (output_unit>0) THEN
          WRITE(output_unit,*)"Shapiro-Wilks normality test (Avg)",avg_test_passed
       END IF
       ! Test of von Neumann for serial correlation
       !                 - Average
       CALL vn_test(wrk,ng_points,r,u,prob,error)
       PRINT *, prob, fe_env%conv_par%vn_conf_lm
       avg_test_passed = prob<=fe_env%conv_par%vn_conf_lm
       fe_env%conv_par%test_vn = avg_test_passed
       IF (output_unit>0) THEN
          WRITE(output_unit,*)"von Neumann serial correlation test (Avg)",avg_test_passed
       END IF
    END DO
    ! If tests on average are ok let's proceed with Standard Deviation
    IF (fe_env%conv_par%test_vn.AND.fe_env%conv_par%test_sw) THEN
       ! Testing Standard Deviations
       DO j = 1, ncolvar
          DO k = j, ncolvar
             DO i = 1, ng_points
                wrk(i) = fe_env%cg_data(i)%var(j,k)
             END DO
             ! Test of Shapiro - Wilks for normality
             !                 - Standard Deviation
             CALL sw_test(wrk, ng_points, w, pw, error)
             PRINT *, 1.0_dp-pw,fe_env%conv_par%sw_conf_lm
             sdv_test_passed = (1.0_dp-pw)<=fe_env%conv_par%sw_conf_lm
             fe_env%conv_par%test_sw = fe_env%conv_par%test_sw .AND. sdv_test_passed
             IF (output_unit>0) THEN
                WRITE(output_unit,*)"Shapiro-Wilks normality test (Std. Dev.)",sdv_test_passed
             END IF
             ! Test of von Neumann for serial correlation
             !                 - Standard Deviation
             CALL vn_test(wrk,ng_points,r,u,prob,error)
             PRINT *, prob, fe_env%conv_par%vn_conf_lm
             sdv_test_passed = prob<=fe_env%conv_par%vn_conf_lm
             fe_env%conv_par%test_vn = fe_env%conv_par%test_vn .AND. sdv_test_passed
             IF (output_unit>0) THEN
                WRITE(output_unit,*)"von Neumann serial correlation test (Std. Dev.)",sdv_test_passed
             END IF
          END DO
       END DO
       CALL destroy_tmp_data(fe_env,wrk,ng_points,error)
    ELSE
       CALL destroy_tmp_data(fe_env,wrk,ng_points,error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ui_check_norm_sc_low

! *****************************************************************************
!> \brief Convergence criteria (Error on average and covariance matrix)
!>      for free energy method
!> \param fe_env ...
!> \param converged ...
!> \param nr_points ...
!> \param output_unit ...
!> \param error ...
!> \par History
!>      Teodoro Laino (01.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE ui_check_convergence(fe_env,converged,nr_points,output_unit,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    LOGICAL, INTENT(OUT)                     :: converged
    INTEGER, INTENT(IN)                      :: nr_points, output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ui_check_convergence', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ic, ncolvar, &
                                                ng_points, stat
    LOGICAL                                  :: failure, test_passed
    REAL(KIND=dp)                            :: max_error_avg, max_error_std
    REAL(KIND=dp), DIMENSION(:), POINTER     :: avg_std, avgmx
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cov_std, covmx

    CALL timeset(routineN,handle)
    failure   = .FALSE.
    converged = .FALSE.
    ncolvar   = fe_env%ncolvar
    NULLIFY(avgmx, avg_std, covmx, cov_std)
    CALL ui_check_norm_sc(fe_env,test_passed,nr_points,output_unit,error)
    IF (test_passed) THEN
       ng_points = nr_points/fe_env%conv_par%cg_width
       ! We can finally compute the error on average and covariance matrix
       ! and check if we converged..
       CALL create_tmp_data(fe_env,ng_points=ng_points,ncolvar=ncolvar,error=error)
       CALL create_csg_data(fe_env,ng_points,error)
       ALLOCATE(covmx(ncolvar,ncolvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(avgmx(ncolvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(cov_std(ncolvar*(ncolvar+1)/2,ncolvar*(ncolvar+1)/2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(avg_std(ncolvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       covmx = 0.0_dp
       avgmx = 0.0_dp
       DO i = 1, ng_points
          covmx = covmx + fe_env%cg_data(i)%var
          avgmx = avgmx + fe_env%cg_data(i)%avg
       END DO
       covmx = covmx / REAL(ng_points, KIND=dp)
       avgmx = avgmx / REAL(ng_points, KIND=dp)

       ! Compute errors on average and standard deviation
       CALL compute_avg_std_errors(fe_env, ncolvar, avgmx, covmx, avg_std, cov_std, error)
       WRITE(*,*)"pippo",avgmx, covmx
       WRITE(*,*)"pippo",avg_std, cov_std
       ! Convergence of the averages
       max_error_avg = SQRT(MAXVAL(ABS(avg_std))/REAL(ng_points,KIND=dp))/MINVAL(avgmx)
       max_error_std = SQRT(MAXVAL(ABS(cov_std))/REAL(ng_points,KIND=dp))/MINVAL(covmx)
       IF  (max_error_avg <=fe_env%conv_par%eps_conv.AND.&
            max_error_std <=fe_env%conv_par%eps_conv ) converged = .TRUE.

       IF (output_unit>0) THEN
          WRITE(output_unit,'(/,T2,"CG SAMPLING LENGTH = ",I7,20X,"REQUESTED ACCURACY  = ",E12.6)')ng_points,&
               fe_env%conv_par%eps_conv
          WRITE(output_unit,'(T50,"PRESENT ACCURACY AVG= ",E12.6)')max_error_avg
          WRITE(output_unit,'(T50,"PRESENT ACCURACY STD= ",E12.6)')max_error_std
          WRITE(output_unit,'(T50,"CONVERGED FE-DER = ",L12)')converged

          WRITE(output_unit,'(/,T33, "COVARIANCE MATRIX")')
          WRITE(output_unit,'(T8,'//cp_to_string(ncolvar)//'(3X,I7,6X))')(ic,ic=1,ncolvar)
          DO ic = 1, ncolvar
             WRITE(output_unit,'(T2,I6,'//cp_to_string(ncolvar)//'(3X,E12.6,1X))')ic,covmx(ic,:)
          END DO
          WRITE(output_unit,'(T33, "ERROR OF COVARIANCE MATRIX")')
          WRITE(output_unit,'(T8,'//cp_to_string(ncolvar)//'(3X,I7,6X))')(ic,ic=1,ncolvar)
          DO ic = 1, ncolvar
             WRITE(output_unit,'(T2,I6,'//cp_to_string(ncolvar)//'(3X,E12.6,1X))')ic,cov_std(ic,:)
          END DO

          WRITE(output_unit,'(/,T2,"COLVAR Nr.",18X,13X,"AVERAGE",13X,"STANDARD DEVIATION")')
          WRITE(output_unit,'(T2,"CV",I8,21X,7X,E12.6,14X,E12.6)')&
               (ic,avgmx(ic),SQRT(ABS(avg_std(ic))),ic=1,ncolvar)
       END IF
       CALL destroy_tmp_data(fe_env,ng_points=ng_points,error=error)
       DEALLOCATE(covmx,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(avgmx,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cov_std,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(avg_std,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ui_check_convergence

! *****************************************************************************
!> \brief Computes the errors on averages and standard deviations for a
!>      correlation-independent coarse grained data set
!> \param fe_env ...
!> \param ncolvar ...
!> \param avgmx ...
!> \param covmx ...
!> \param avg_std ...
!> \param cov_std ...
!> \param error ...
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE compute_avg_std_errors(fe_env, ncolvar, avgmx, covmx, avg_std, cov_std, error)
    TYPE(free_energy_type), POINTER          :: fe_env
    INTEGER, INTENT(IN)                      :: ncolvar
    REAL(KIND=dp), DIMENSION(:), POINTER     :: avgmx
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: covmx
    REAL(KIND=dp), DIMENSION(:), POINTER     :: avg_std
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cov_std
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_avg_std_errors', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind, j, k, nvar, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fac
    REAL(KIND=dp), DIMENSION(:), POINTER     :: awrk, eig, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: wrk

    failure = .FALSE.
    ! Averages
    nvar = ncolvar
    ALLOCATE(wrk(nvar,nvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eig(nvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    fac = REAL(SIZE(fe_env%cg_data),KIND=dp)
    wrk = 0.0_dp
    eig = 0.0_dp
    DO k = 1, SIZE(fe_env%cg_data)
       DO j=1,nvar
          DO i = j,nvar
             wrk(i,j) = wrk(i,j) + fe_env%cg_data(k)%avg(i)*fe_env%cg_data(k)%avg(j)
          END DO
       END DO
    END DO
    DO j=1,nvar
       DO i = j,nvar
          wrk(i,j) = wrk(i,j) - avgmx(i)*avgmx(j)*fac
          wrk(j,i) = wrk(i,j)
       END DO
    END DO
    wrk = wrk / (fac-1.0_dp)
    ! Diagonalize the covariance matrix and check for the maximum error
    CALL diamat_all(wrk,eig,error=error)
    DO i = 1, nvar
       avg_std(i) = eig(i)
    END DO
    DEALLOCATE(wrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eig,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Standard Deviations
    nvar = ncolvar*(ncolvar+1)/2
    ALLOCATE(wrk(nvar,nvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eig(nvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(awrk(nvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp(nvar),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    wrk = 0.0_dp
    eig = 0.0_dp
    ind = 0
    DO i= 1, ncolvar
       DO j = i, ncolvar
          ind = ind+1
          awrk(ind) = covmx(i,j)
       END DO
    END DO
    DO k = 1, SIZE(fe_env%cg_data)
       ind = 0
       DO i= 1, ncolvar
          DO j = i, ncolvar
             ind = ind+1
             tmp(ind) = fe_env%cg_data(k)%var(i,j)
          END DO
       END DO
       DO i=1,nvar
          DO j = i,nvar
             wrk(i,j) = wrk(i,j) + tmp(i)*tmp(j) - awrk(i)*awrk(j)
          END DO
       END DO
    END DO
    DO i=1,nvar
       DO j = i,nvar
          wrk(i,j) = wrk(i,j) - fac*awrk(i)*awrk(j)
          wrk(j,i) = wrk(i,j)
       END DO
    END DO
    wrk = wrk / (fac-1.0_dp)
    ! Diagonalize the covariance matrix and check for the maximum error
    CALL diamat_all(wrk,eig,error=error)
    ind = 0
    DO i = 1, ncolvar
       DO j = i, ncolvar
          ind = ind + 1
          cov_std(i,j) = eig(ind)
          cov_std(j,i) = cov_std(i,j)
       END DO
    END DO
    DEALLOCATE(wrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eig,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(awrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE compute_avg_std_errors

! *****************************************************************************
!> \brief Computes the covariance matrix
!> \param fe_env ...
!> \param cg_index ...
!> \param istart ...
!> \param iend ...
!> \param covmx ...
!> \param avgs ...
!> \param error ...
!> \par History
!>      Teodoro Laino (01.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE eval_cov_matrix(fe_env,cg_index,istart,iend,covmx,avgs,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    INTEGER, INTENT(IN)                      :: cg_index, istart, iend
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: covmx
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: avgs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_cov_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ic, jc, jstep, &
                                                ncolvar, nlength
    REAL(KIND=dp)                            :: tmp_ic, tmp_jc
    TYPE(ui_var_type), POINTER               :: cv

    CALL timeset(routineN,handle)
    ncolvar   = fe_env%ncolvar
    nlength   = iend-istart+1
    fe_env%cg_data(cg_index)%avg = 0.0_dp
    fe_env%cg_data(cg_index)%var = 0.0_dp
    IF (nlength>1) THEN
       ! Update the info on averages and variances
       DO jstep = istart, iend
          DO ic=1,ncolvar
             cv => fe_env%uivar(ic)
             tmp_ic = cv%ss(jstep)
             fe_env%cg_data(cg_index)%avg(ic) = fe_env%cg_data(cg_index)%avg(ic) + tmp_ic
          END DO
          DO ic=1,ncolvar
             cv => fe_env%uivar(ic)
             tmp_ic = cv%ss(jstep)
             DO jc = 1, ic
                cv => fe_env%uivar(jc)
                tmp_jc = cv%ss(jstep)
                fe_env%cg_data(cg_index)%var(jc,ic) = fe_env%cg_data(cg_index)%var(jc,ic) + tmp_ic*tmp_jc
             END DO
          END DO
       END DO
       ! Normalized the variances and the averages
       ! Unbiased estimator
       fe_env%cg_data(cg_index)%var=fe_env%cg_data(cg_index)%var/REAL(nlength-1,KIND=dp)
       fe_env%cg_data(cg_index)%avg=fe_env%cg_data(cg_index)%avg/REAL(nlength,KIND=dp)
       ! Compute the covariance matrix
       DO ic=1,ncolvar
          tmp_ic = fe_env%cg_data(cg_index)%avg(ic)
          DO jc = 1, ic
             tmp_jc = fe_env%cg_data(cg_index)%avg(jc) * REAL(nlength,KIND=dp)/REAL(nlength-1,KIND=dp)
             fe_env%cg_data(cg_index)%var(jc,ic) = fe_env%cg_data(cg_index)%var(jc,ic) - tmp_ic*tmp_jc
             fe_env%cg_data(cg_index)%var(ic,jc) = fe_env%cg_data(cg_index)%var(jc,ic)
          END DO
       END DO
       WRITE(*,*)"eval_cov_matrix",istart, iend,  fe_env%cg_data(cg_index)%avg, fe_env%cg_data(cg_index)%var
       IF (PRESENT(covmx)) covmx = fe_env%cg_data(cg_index)%var
       IF (PRESENT(avgs))  avgs  = fe_env%cg_data(cg_index)%avg
    END IF
    CALL timestop(handle)
  END SUBROUTINE eval_cov_matrix

! *****************************************************************************
!> \brief Dumps information when performing an alchemical change run
!> \param my_val ...
!> \param my_par ...
!> \param dx ...
!> \param lerr ...
!> \param fe_section ...
!> \param nforce_eval ...
!> \param cum_res ...
!> \param istep ...
!> \param beta ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino - University of Zurich [tlaino] - 05.2007
! *****************************************************************************
  SUBROUTINE dump_ac_info(my_val, my_par, dx, lerr, fe_section, nforce_eval, cum_res,&
       istep, beta, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par
    REAL(KIND=dp), INTENT(IN)                :: dx, lerr
    TYPE(section_vals_type), POINTER         :: fe_section
    INTEGER, INTENT(IN)                      :: nforce_eval
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cum_res
    INTEGER, POINTER                         :: istep
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dump_ac_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, par, this_error
    INTEGER                                  :: i, iforce_eval, ipar, isize, &
                                                iw, j, NEquilStep
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: avg_BP, avg_DET, avg_DUE, d_ene_w, dedf, ene_w, err, &
      Err_DET, Err_DUE, std_DET, std_DUE, tmp, tmp2, wfac
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: alch_section

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    alch_section => section_vals_get_subs_vals(fe_section,"ALCHEMICAL_CHANGE",error=error)
    CALL section_vals_val_get(alch_section,"PARAMETER",c_val=par,error=error)
    DO i = 1, SIZE(my_par)
       IF (my_par(i)==par) EXIT
    END DO
    CPPrecondition(i<=SIZE(my_par),cp_failure_level,routineP,error,failure)
    ipar = i
    dedf = evalfd(1,ipar,my_val,dx,err)
    IF (ABS(err)>lerr) THEN
       WRITE(this_error,"(A,G12.6,A)")"(",err,")"
       WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
       CALL compress(this_error,.TRUE.)
       CALL compress(def_error,.TRUE.)
       CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
            'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
            ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
            TRIM(def_error)//' .',error=error,only_ionode=.TRUE.)
    END IF

    ! We must print now the energy of the biased system, the weigthing energy
    ! and the derivative w.r.t.the coupling parameter of the biased energy
    ! Retrieve the expression of the weighting function:
    CALL section_vals_val_get(alch_section,"WEIGHTING_FUNCTION",c_val=coupling_function,error=error)
    CALL compress(coupling_function, full=.TRUE.)
    CALL parsef(2,TRIM(coupling_function),my_par)
    ene_w   = evalf(2,my_val)
    d_ene_w = evalfd(2,ipar,my_val,dx,err)
    IF (ABS(err)>lerr) THEN
       WRITE(this_error,"(A,G12.6,A)")"(",err,")"
       WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
       CALL compress(this_error,.TRUE.)
       CALL compress(def_error,.TRUE.)
       CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
            'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
            ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
            TRIM(def_error)//' .',error=error,only_ionode=.TRUE.)
    END IF
    CALL section_vals_val_get(alch_section,"NEQUIL_STEPS",i_val=NEquilStep,error=error)
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
       iw = cp_print_key_unit_nr(logger,fe_section,"FREE_ENERGY_INFO",&
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
    CALL cp_print_key_finished_output(iw,logger,fe_section,"FREE_ENERGY_INFO",error=error)

  END SUBROUTINE dump_ac_info

END MODULE free_energy_methods

