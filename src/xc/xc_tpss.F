!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the tpss functional.
!> \note
!>      The derivation of the formulaes is lengthly, and not fully trivial,
!>      so I have put it in doc/tpss.mw
!> \par History
!>      05.2004 created
!> \author fawzi
! *****************************************************************************
MODULE xc_tpss
  USE bibliography,                    ONLY: Tao2003,&
                                             cite_reference
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_tpss'

  PUBLIC :: tpss_lda_info, tpss_lda_eval, tpss_lsd_info, tpss_lsd_eval

!***
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param tpss_params ...
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv the highest derivative available
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE tpss_lda_info(tpss_params,reference,shortform, needs, max_deriv, error)
    TYPE(section_vals_type), POINTER         :: tpss_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_info', &
      routineP = moduleN//':'//routineN

    REAL(kind=dp)                            :: sc, sx

    CALL section_vals_val_get(tpss_params,"SCALE_C",r_val=sc,error=error)
    CALL section_vals_val_get(tpss_params,"SCALE_X",r_val=sx,error=error)

    IF ( PRESENT ( reference ) ) THEN
       IF (sx==1._dp.AND.sc==1._dp) THEN
          reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LDA version}"
       ELSE
          WRITE (reference,"(a,'sx=',f5.3,'sc=',f5.3,' {LDA version}')")&
               "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003)",&
               sx,sc
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       IF (sx==1._dp.AND.sc==1._dp) THEN
          shortform = "TPSS meta-GGA functional (LDA)"
       ELSE
          WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,' (LDA)')")&
               "TPSS meta-GGA functional",&
               sx,sc
       END IF
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%tau=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE tpss_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param tpss_params ...
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv the highest derivative available
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE tpss_lsd_info(tpss_params,reference,shortform, needs, max_deriv, error)
    TYPE(section_vals_type), POINTER         :: tpss_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_info', &
      routineP = moduleN//':'//routineN

    REAL(kind=dp)                            :: sc, sx

    CALL section_vals_val_get(tpss_params,"SCALE_C",r_val=sc,error=error)
    CALL section_vals_val_get(tpss_params,"SCALE_X",r_val=sx,error=error)

    IF ( PRESENT ( reference ) ) THEN
       IF (sx==1._dp.AND.sc==1._dp) THEN
          reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LSD version}"
       ELSE
          WRITE (reference,"(a,'sx=',f5.3,'sc=',f5.3,' {LSD version}')")&
               "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003)",&
               sx,sc
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       IF (sx==1._dp.AND.sc==1._dp) THEN
          shortform = "TPSS meta-GGA functional (LSD)"
       ELSE
          WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,' (LSD)')")&
               "TPSS meta-GGA functional",&
               sx,sc
       END IF
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%tau_spin=.TRUE.
       needs%norm_drho=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE tpss_lsd_info

! *****************************************************************************
!> \brief evaluates the tpss functional in the spin unpolarized (lda) case
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param tpss_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE tpss_lda_eval(rho_set,deriv_set,grad_deriv,tpss_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: tpss_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ii, non_coer, &
                                                npoints, stat, unit_nr
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, epsilon_tau, &
                                                scale_ec, scale_ex
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrho, e_rho, &
                                                e_tau, norm_drho, rho, tau
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL cite_reference(tao2003)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
            tau=tau,drho_cutoff=epsilon_norm_drho,tau_cutoff=epsilon_tau,&
            error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       ! meaningful default for the arrays we don't need: let us make compiler
       ! and debugger happy...
       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> rho
       END IF

       e_0 => dummy
       e_rho => dummy
       e_ndrho => dummy
       e_tau => dummy

       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_tau,error=error)
       END IF
       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       non_coer = 0
       CALL section_vals_val_get(tpss_params,"SCALE_C",r_val=scale_ec,error=error)
       CALL section_vals_val_get(tpss_params,"SCALE_X",r_val=scale_ex,error=error)

       !$omp parallel default(none) &
       !$omp          shared(rho, tau, norm_drho, e_0, e_rho, e_ndrho, e_tau) &
       !$omp          shared(epsilon_rho, epsilon_tau, npoints, grad_deriv) &
       !$omp          shared(epsilon_norm_drho, scale_ec, scale_ex, error) &
       !$omp          reduction(+: non_coer)

       CALL tpss_lda_calc(rho=rho, norm_drho=norm_drho,&
            tau=tau,e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_tau=e_tau,&
            grad_deriv=grad_deriv, npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_norm_drho=epsilon_norm_drho,epsilon_tau=epsilon_tau,&
            scale_ec=scale_ec,scale_ex=scale_ex,non_coer=non_coer,error=error)

       !$omp end parallel

       logger => cp_error_get_logger(error)

       IF (cp_debug) THEN
          unit_nr=cp_logger_get_default_io_unit(logger)
          CALL mp_sum(non_coer,logger%para_env%group)
          ii=npoints
          CALL mp_sum(ii,logger%para_env%group)
          IF (unit_nr>0) THEN
             WRITE(unit_nr,"('non_coer: ',i10,' of ',i10,' =',f6.3,'%')")&
                  non_coer,2*ii,REAL(non_coer,dp)/REAL(ii,dp)*100._dp
          END IF
       END IF

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE tpss_lda_eval

! *****************************************************************************
!> \brief evaluates the tpss functional in the spin polarized (lsd) case
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param tpss_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE tpss_lsd_eval(rho_set,deriv_set,grad_deriv,tpss_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: tpss_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ii, non_coer, &
                                                npoints, stat, unit_nr
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, epsilon_tau, &
                                                scale_ec, scale_ex
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrhoa, e_ndrhob, e_rhoa, e_rhob, e_tau_a, e_tau_b, norm_drho, &
      norm_drhoa, norm_drhob, rhoa, rhob, tau_a, tau_b
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL cite_reference(tao2003)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,norm_drhoa=norm_drhoa,&
            norm_drhob=norm_drhob,norm_drho=norm_drho,tau_a=tau_a,tau_b=tau_b,&
            local_bounds=bo,rho_cutoff=epsilon_rho,&
            drho_cutoff=epsilon_norm_drho,tau_cutoff=epsilon_tau,&
            error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       ! meaningful default for the arrays we don't need: let us make compiler
       ! and debugger happy...
       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> rhoa
       END IF

       e_0 => dummy
       e_rhoa => dummy
       e_rhob => dummy
       e_ndrhoa => dummy
       e_ndrhob => dummy
       e_ndrho => dummy
       e_tau_a => dummy
       e_tau_b => dummy

       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rhoa,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rhob,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrhob,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(tau_a)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_tau_a,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(tau_b)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_tau_b,error=error)
       END IF
       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       non_coer = 0
       CALL section_vals_val_get(tpss_params,"SCALE_C",r_val=scale_ec,error=error)
       CALL section_vals_val_get(tpss_params,"SCALE_X",r_val=scale_ex,error=error)

       !$omp parallel default(none) &
       !$omp          shared(rhoa, rhob,norm_drhoa, norm_drhob, norm_drho) &
       !$omp          shared(tau_a, tau_b, e_0, e_rhoa, e_rhob, e_ndrhoa) &
       !$omp          shared(e_ndrhob, e_ndrho, e_tau_a, e_tau_b, npoints) &
       !$omp          shared(grad_deriv, epsilon_rho, epsilon_norm_drho) &
       !$omp          shared(epsilon_tau, scale_ec, scale_ex, error) &
       !$omp          reduction(+: non_coer)

       CALL tpss_lsd_calc(rhoa=rhoa,rhob=rhob, norm_drhoa=norm_drhoa, &
            norm_drhob=norm_drhob, norm_drho=norm_drho,&
            tau_a=tau_a,tau_b=tau_b,e_0=e_0,e_rhoa=e_rhoa,e_rhob=e_rhob,&
            e_ndrhoa=e_ndrhoa,e_ndrhob=e_ndrhob,e_ndrho=e_ndrho,&
            e_tau_a=e_tau_a,e_tau_b=e_tau_b,&
            grad_deriv=grad_deriv, npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_norm_drho=epsilon_norm_drho,epsilon_tau=epsilon_tau,&
            scale_ec=scale_ec,scale_ex=scale_ex,non_coer=non_coer,error=error)

       !$omp end parallel

       logger => cp_error_get_logger(error)

       IF (cp_debug) THEN
          unit_nr=cp_logger_get_default_io_unit(logger)
          CALL mp_sum(non_coer,logger%para_env%group)
          ii=npoints
          CALL mp_sum(ii,logger%para_env%group)
          IF (unit_nr>0) THEN
             WRITE(unit_nr,"('non_coer: ',i10,' of ',i10,' =',f6.3,'%')")&
                  non_coer,2*ii,REAL(non_coer,dp)/REAL(2*ii,dp)*100._dp
          END IF
       END IF

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE tpss_lsd_eval

! *****************************************************************************
!> \brief low level calculation routine for the unpolarized (lda) tpss
!> \param rho ...
!> \param norm_drho ...
!> \param tau ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param npoints ...
!> \param grad_deriv ...
!> \param epsilon_rho ...
!> \param epsilon_norm_drho ...
!> \param epsilon_tau ...
!> \param scale_ec ...
!> \param scale_ex ...
!> \param non_coer ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      maple is nice, but if you want the uman readable version of the code
!>      look in doc/tpss.mw
! *****************************************************************************
  SUBROUTINE tpss_lda_calc(rho,norm_drho,tau,e_0,e_rho,e_ndrho,e_tau,&
       npoints,grad_deriv,epsilon_rho,epsilon_norm_drho,epsilon_tau,&
       scale_ec,scale_ex,non_coer,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho, norm_drho, tau
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rho, e_ndrho, e_tau
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, &
                                                epsilon_tau, scale_ec, &
                                                scale_ex
    INTEGER, INTENT(inout)                   :: non_coer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: abs_grad_deriv, ii
    LOGICAL                                  :: failure, t571, t639
    REAL(kind=dp) :: A, A_1, A_2, A_s1, A_s1rho, A_s2, A_s2rho, alpha, &
      alpha_1_1, alpha_1_2, alphanorm_drho, alpharho, alphatau, Arho, b, &
      beta, beta_1_1, beta_1_2, beta_2_1, beta_2_2, beta_3_1, beta_3_2, &
      beta_4_1, beta_4_2, c, d, e_c_u_0, e_c_u_0rho, e_c_u_1_s1, &
      e_c_u_1_s1rho, e_c_u_1_s2, e_c_u_1_s2rho, e_var, epsilon_cGGA, &
      epsilon_cGGA_0_1, epsilon_cGGA_1_0, epsilon_cGGArho, epsilon_cRevPKZB, &
      epsilon_cRevPKZBnorm_drho, epsilon_cRevPKZBrho, epsilon_cRevPKZBtau, &
      ex_unif, Fx, gamma_var, Hnorm_drho, k_f_s1, k_f_s1rho, k_s, k_s_s1, &
      k_s_s2, kappa, m, ma, manorm_drho, marho, mb, mbnorm_drho, mbrho, mu
    REAL(kind=dp) :: my_ndrho, my_rho, my_tau, p, p_1, p_2, p_3, phi_s1, &
      phi_s2, pnorm_drho, prho, rs, rs_s1, rs_s1rho, rs_s2, rs_s2rho, rsrho, &
      t, t1, t100, t101, t111, t12, t13, t138, t14, t140, t143, t145, t146, &
      t147, t151, t152, t16, t161, t168, t177, t186, t187, t189, t19, t190, &
      t191, t193, t194, t196, t197, t198, t199, t2, t20, t201, t202, t204, &
      t205, t208, t209, t21, t211, t212, t213, t215, t216, t218, t219, t22, &
      t220, t221, t223, t224, t226, t227, t230, t231, t233, t234, t235, t238, &
      t239, t241, t242, t243, t245, t246, t248, t249, t252, t253, t254, t256, &
      t26, t260, t263, t264, t265, t267
    REAL(kind=dp) :: t268, t269, t27, t271, t272, t274, t275, t276, t277, &
      t278, t279, t28, t280, t281, t284, t286, t288, t29, t290, t291, t293, &
      t294, t295, t299, t3, t301, t302, t303, t305, t307, t310, t313, t316, &
      t319, t322, t325, t327, t328, t329, t331, t337, t340, t343, t344, t35, &
      t351, t36, t370, t371, t376, t383, t385, t386, t39, t390, t391, t395, &
      t396, t398, t4, t403, t404, t406, t41, t410, t411, t419, t42, t430, &
      t437, t445, t450, t452, t464, t472, t475, t485, t489, t49, t490, t5, &
      t505, t513, t517, t536, t541, t542, t546, t547, t549, t55, t554, t555, &
      t557, t561, t562, t569, t574, t58, t585
    REAL(kind=dp) :: t6, t60, t604, t609, t610, t614, t615, t617, t622, t623, &
      t625, t629, t630, t637, t642, t645, t659, t67, t7, t71, t73, t77, t78, &
      t79, t799, t80, t84, t85, t89, t9, t94, t95, t96, t_s1, t_s1norm_drho, &
      t_s1rho, t_s2, t_s2norm_drho, t_s2rho, tau_w, tau_wnorm_drho, tau_wrho, &
      tildeq_b, tildeq_bnorm_drho, tildeq_brho, tildeq_btau, tnorm_drho, &
      trho, z, znorm_drho, zrho, ztau

    failure=.FALSE.

    IF (.FALSE.) THEN
       ! useful for testing, we just hack in a well defined functional of tau, ndrho and rho
       ! and see that things converge properly with OT.
       !$omp do
       DO ii=1,npoints
          my_tau=tau(ii)
          my_rho=rho(ii)
          my_ndrho=norm_drho(ii)
          IF (grad_deriv>=0) THEN
                   e_0(ii) = e_0(ii)+my_tau*my_ndrho*my_rho
          END IF
          IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
             e_rho(ii)=e_rho(ii)+my_tau*my_ndrho
             e_ndrho(ii)=e_ndrho(ii)+my_tau*my_rho
             e_tau(ii)=e_tau(ii)+my_rho*my_ndrho
          ENDIF
       ENDDO
       !$omp end do
       RETURN
    ENDIF

    IF (.NOT. failure) THEN
       abs_grad_deriv=ABS(grad_deriv)

       kappa = 0.804e0_dp
       beta = 0.66725e-1_dp
       mu = 0.21951e0_dp
       gamma_var = (0.1e1_dp - LOG(0.2e1_dp)) / pi ** 2
       b = 0.4e0_dp
       c = 0.159096e1_dp
       e_var = 0.1537e1_dp
       d = 0.28e1_dp
       p_1 = 0.10e1_dp
       A_1 = 0.31091e-1_dp
       alpha_1_1 = 0.21370e0_dp
       beta_1_1 = 0.75957e1_dp
       beta_2_1 = 0.35876e1_dp
       beta_3_1 = 0.16382e1_dp
       beta_4_1 = 0.49294e0_dp
       p_2 = 0.10e1_dp
       A_2 = 0.15545e-1_dp
       alpha_1_2 = 0.20548e0_dp
       beta_1_2 = 0.141189e2_dp
       beta_2_2 = 0.61977e1_dp
       beta_3_2 = 0.33662e1_dp
       beta_4_2 = 0.62517e0_dp
       p_3 = 0.10e1_dp

       t1 = 3._dp ** (0.1e1_dp / 0.3e1_dp)
       t2 = 4._dp ** (0.1e1_dp / 0.3e1_dp)
       t3 = t2 ** 2
       t4 = t1 * t3
       t5 = 2._dp ** (0.1e1_dp / 0.3e1_dp)
       t6 = 0.1e1_dp / pi
       t12 = t5 ** 2

       !$omp do

       DO ii=1,npoints
          my_tau=tau(ii)
          my_rho=rho(ii)
          IF (my_rho>epsilon_rho.AND.my_tau>epsilon_tau) THEN
             my_ndrho=norm_drho(ii)

             t7 = 0.1e1_dp / my_rho
             t254 = my_ndrho ** 2
             tau_w = t254 * t7 / 0.8e1_dp

             IF (my_tau<tau_w) THEN
                ! enforce z=norm_rho**2/(8._dp*rho*tau) <1
                m=0.5_dp*t254+4.0_dp*my_rho*my_tau
                my_tau=m/8._dp/my_rho
                my_ndrho=SQRT(m)
                t254=m
                non_coer=non_coer+1
             END IF

             t9 = (t6 * t7) ** (0.1e1_dp / 0.3e1_dp)
             rs_s1 = t4 * t5 * t9 / 0.4e1_dp
             phi_s1 = t12 / 0.2e1_dp
             t13 = t1 * t12
             t14 = pi ** 2
             t16 = (t14 * my_rho) ** (0.1e1_dp / 0.3e1_dp)
             k_f_s1 = t13 * t16 / 0.2e1_dp
             t19 = SQRT(k_f_s1 * t6)
             k_s_s1 = 0.2e1_dp * t19
             t20 = 0.1e1_dp / phi_s1
             t21 = my_ndrho * t20
             t22 = 0.1e1_dp / k_s_s1
             t_s1 = t21 * t22 * t7 / 0.2e1_dp
             rs_s2 = rs_s1
             phi_s2 = phi_s1
             t26 = SQRT(k_f_s1 * t6)
             k_s_s2 = 0.2e1_dp * t26
             t27 = 0.1e1_dp / phi_s2
             t28 = my_ndrho * t27
             t29 = 0.1e1_dp / k_s_s2
             t_s2 = t28 * t29 * t7 / 0.2e1_dp
             t35 = 0.1e1_dp / A_1
             t36 = SQRT(rs_s2)
             t39 = t36 * rs_s2
             t41 = p_1 + 0.1e1_dp
             t42 = rs_s2 ** t41
             t49 = LOG(0.1e1_dp + t35 / (beta_1_1 * t36 + beta_2_1 * rs_s2 + &
                  beta_3_1 * t39 + beta_4_1 * t42) / 0.2e1_dp)
             t55 = SQRT(rs_s1)
             t58 = t55 * rs_s1
             t60 = rs_s1 ** t41
             t67 = LOG(0.1e1_dp + t35 / (beta_1_1 * t55 + beta_2_1 * rs_s1 + &
                  beta_3_1 * t58 + beta_4_1 * t60) / 0.2e1_dp)
             t71 = 0.1e1_dp + alpha_1_2 * rs_s2
             t73 = 0.1e1_dp / A_2
             t77 = p_2 + 0.1e1_dp
             t78 = rs_s2 ** t77
             t79 = beta_4_2 * t78
             t80 = beta_1_2 * t36 + beta_2_2 * rs_s2 + beta_3_2 * t39 + t79
             t84 = 0.1e1_dp + t73 / t80 / 0.2e1_dp
             t85 = LOG(t84)
             e_c_u_1_s2 = -0.2e1_dp * A_2 * t71 * t85
             t89 = 0.1e1_dp + alpha_1_2 * rs_s1
             t94 = rs_s1 ** t77
             t95 = beta_4_2 * t94
             t96 = beta_1_2 * t55 + beta_2_2 * rs_s1 + beta_3_2 * t58 + t95
             t100 = 0.1e1_dp + t73 / t96 / 0.2e1_dp
             t101 = LOG(t100)
             e_c_u_1_s1 = -0.2e1_dp * A_2 * t89 * t101
             t111 = p_3 + 1._dp
             rs = t4 * t9 / 0.4e1_dp
             t138 = 0.1e1_dp + alpha_1_1 * rs
             t140 = SQRT(rs)
             t143 = t140 * rs
             t145 = rs ** t41
             t146 = beta_4_1 * t145
             t147 = beta_1_1 * t140 + beta_2_1 * rs + beta_3_1 * t143 + t146
             t151 = 0.1e1_dp + t35 / t147 / 0.2e1_dp
             t152 = LOG(t151)
             e_c_u_0 = -0.2e1_dp * A_1 * t138 * t152
             t161 = rs ** t77
             t168 = LOG(0.1e1_dp + t73 / (beta_1_2 * t140 + beta_2_2 * rs + &
                  beta_3_2 * t143 + beta_4_2 * t161) / 0.2e1_dp)
             t177 = rs ** t111
             t186 = 0.1e1_dp / gamma_var
             t187 = beta * t186
             t189 = phi_s1 ** 2
             t190 = t189 * phi_s1
             t191 = 0.1e1_dp / t190
             t193 = EXP(-e_c_u_1_s1 * t186 * t191)
             t194 = t193 - 0.1e1_dp
             A_s1 = t187 / t194
             t196 = gamma_var * t190
             t197 = t_s1 ** 2
             t198 = A_s1 * t197
             t199 = 0.1e1_dp + t198
             t201 = A_s1 ** 2
             t202 = t197 ** 2
             t204 = 0.1e1_dp + t198 + t201 * t202
             t205 = 0.1e1_dp / t204
             t208 = 0.1e1_dp + t187 * t197 * t199 * t205
             t209 = LOG(t208)
             epsilon_cGGA_1_0 = e_c_u_1_s1 + t196 * t209
             t211 = phi_s2 ** 2
             t212 = t211 * phi_s2
             t213 = 0.1e1_dp / t212
             t215 = EXP(-e_c_u_1_s2 * t186 * t213)
             t216 = t215 - 0.1e1_dp
             A_s2 = t187 / t216
             t218 = gamma_var * t212
             t219 = t_s2 ** 2
             t220 = A_s2 * t219
             t221 = t220 + 0.1e1_dp
             t223 = A_s2 ** 2
             t224 = t219 ** 2
             t226 = 0.1e1_dp + t220 + t223 * t224
             t227 = 0.1e1_dp / t226
             t230 = 0.1e1_dp + t187 * t219 * t221 * t227
             t231 = LOG(t230)
             epsilon_cGGA_0_1 = e_c_u_1_s2 + t218 * t231
             t233 = SQRT(t1 * t16 * t6)
             k_s = 0.2e1_dp * t233
             t234 = 0.1e1_dp / k_s
             t235 = my_ndrho * t234
             t = t235 * t7 / 0.2e1_dp
             t238 = EXP(-e_c_u_0 * t186)
             t239 = -0.1e1_dp + t238
             A = t187 / t239
             t241 = t ** 2
             t242 = A * t241
             t243 = 0.1e1_dp + t242
             t245 = A ** 2
             t246 = t241 ** 2
             t248 = 0.1e1_dp + t242 + t245 * t246
             t249 = 0.1e1_dp / t248
             t252 = 0.1e1_dp + t187 * t241 * t243 * t249
             t253 = LOG(t252)
             epsilon_cGGA = e_c_u_0 + gamma_var * t253
             ma = MAX(epsilon_cGGA_1_0, epsilon_cGGA)
             mb = MAX(epsilon_cGGA_0_1, epsilon_cGGA)
             t256 = tau_w ** 2
             t260 = ma / 0.2e1_dp + mb / 0.2e1_dp
             t263 = 0.53e0_dp * epsilon_cGGA * t256 - 0.153e1_dp * t256 * t260
             t264 = my_tau ** 2
             t265 = 0.1e1_dp / t264
             epsilon_cRevPKZB = epsilon_cGGA + t263 * t265
             t267 = my_rho * epsilon_cRevPKZB
             t268 = d * epsilon_cRevPKZB
             t269 = t256 * tau_w
             t271 = 0.1e1_dp / t264 / my_tau
             t272 = t269 * t271
             t274 = 0.1e1_dp + t268 * t272
             t275 = t254 * t1
             t276 = t14 ** (0.1e1_dp / 0.3e1_dp)
             t277 = t276 ** 2
             t278 = 0.1e1_dp / t277
             t279 = my_rho ** 2
             t280 = my_rho ** (0.1e1_dp / 0.3e1_dp)
             t281 = t280 ** 2
             t284 = t278 / t281 / t279
             p = t275 * t284 / 0.12e2_dp
             t286 = 0.1e1_dp / my_tau
             z = tau_w * t286
             t288 = 0.1e1_dp / z - 0.1e1_dp
             alpha = 0.5e1_dp / 0.3e1_dp * p * t288
             t290 = alpha - 0.1e1_dp
             t291 = b * alpha
             t293 = 0.1e1_dp + t291 * t290
             t294 = SQRT(t293)
             t295 = 0.1e1_dp / t294
             tildeq_b = 0.9e1_dp / 0.20e2_dp * t290 * t295 + 0.2e1_dp / 0.3e1_dp * p
             t299 = z ** 2
             t301 = 0.1e1_dp + t299
             t302 = t301 ** 2
             t303 = 0.1e1_dp / t302
             t305 = 0.10e2_dp / 0.81e2_dp + c * t299 * t303
             t307 = tildeq_b ** 2
             t310 = p ** 2
             t313 = SQRT(0.18e2_dp * t299 + 0.50e2_dp * t310)
             t316 = 0.1e1_dp / kappa
             t319 = SQRT(e_var)
             t322 = e_var * mu
             t325 = t305 * p + 0.146e3_dp / 0.2025e4_dp * t307 - 0.73e2_dp / &
                  0.4050e4_dp * tildeq_b * t313 + 0.100e3_dp / 0.6561e4_dp * t316 * &
                  t310 + 0.4e1_dp / 0.45e2_dp * t319 * t299 + t322 * t310 * p
             t327 = 0.1e1_dp + t319 * p
             t328 = t327 ** 2
             t329 = 0.1e1_dp / t328
             t331 = 0.1e1_dp + t325 * t329 * t316
             Fx = 0.1e1_dp + kappa - kappa / t331
             ex_unif = -0.3e1_dp / 0.4e1_dp * t1 * t16 * t6
             t337 = my_rho * ex_unif

             IF (grad_deriv>=0) THEN
                e_0(ii) = e_0(ii)+&
                     scale_ec * t267 * t274 + scale_ex * t337 * Fx
             END IF

             IF (abs_grad_deriv>0) THEN
                t340 = t9 ** 2
                t343 = 0.1e1_dp / t279
                t344 = 0.1e1_dp / t340 * t6 * t343
                rsrho = -t4 * t344 / 0.12e2_dp
                t351 = t147 ** 2
                e_c_u_0rho = -0.2e1_dp * A_1 * alpha_1_1 * rsrho * t152 + t138 / &
                     t351 * (beta_1_1 / t140 * rsrho / 0.2e1_dp + beta_2_1 * rsrho + &
                     0.3e1_dp / 0.2e1_dp * beta_3_1 * t140 * rsrho + t146 * t41 * rsrho /&
                     rs) / t151
                t370 = t16 ** 2
                t371 = 0.1e1_dp / t370
                t376 = k_s ** 2
                trho = -my_ndrho / t376 * t7 / t233 * t1 * t371 * t14 * t6 &
                     / 0.6e1_dp - t235 * t343 / 0.2e1_dp
                t383 = gamma_var ** 2
                t385 = beta / t383
                t386 = t239 ** 2
                Arho = t385 / t386 * e_c_u_0rho * t238
                t390 = t187 * t
                t391 = t243 * t249
                t395 = Arho * t241
                t396 = A * t
                t398 = 0.2e1_dp * t396 * trho
                t403 = t187 * t241
                t404 = t248 ** 2
                t406 = t243 / t404
                t410 = t241 * t
                t411 = t245 * t410
                t419 = 0.1e1_dp / t252
                epsilon_cGGArho = e_c_u_0rho + gamma_var * (0.2e1_dp * t390 * t391&
                     * trho + t187 * t241 * (t395 + t398) * t249 - t403 * t406 * (t395 +&
                     t398 + 0.2e1_dp * A * t246 * Arho + 0.4e1_dp * t411 * trho)) * t419
                tau_wrho = -t254 * t343 / 0.8e1_dp
                prho = -0.2e1_dp / 0.9e1_dp * t275 * t278 / t281 / t279 / my_rho
                zrho = tau_wrho * t286
                t430 = p / t299
                alpharho = 0.5e1_dp / 0.3e1_dp * prho * t288 - 0.5e1_dp / 0.3e1_dp&
                     * t430 * zrho
                t437 = t290 / t294 / t293
                tildeq_brho = 0.9e1_dp / 0.20e2_dp * alpharho * t295 - 0.9e1_dp / &
                     0.40e2_dp * t437 * (b * alpharho * t290 + t291 * alpharho) + &
                     0.2e1_dp / 0.3e1_dp * prho
                t445 = c * z
                t450 = c * t299 * z
                t452 = 0.1e1_dp / t302 / t301
                t464 = tildeq_b / t313
                t472 = t316 * p
                t475 = t319 * z
                t485 = t325 / t328 / t327
                t489 = t331 ** 2
                t490 = 0.1e1_dp / t489
                rs_s1rho = -t4 * t5 * t344 / 0.12e2_dp
                k_f_s1rho = t13 * t371 * t14 / 0.6e1_dp
                t505 = k_s_s1 ** 2
                t_s1rho = -t21 / t505 * t7 / t19 * k_f_s1rho * t6 / 0.2e1_dp - t21&
                     * t22 * t343 / 0.2e1_dp
                t513 = A_2 * alpha_1_2
                t517 = t96 ** 2
                e_c_u_1_s1rho = -0.2e1_dp * t513 * rs_s1rho * t101 + t89 / t517 * &
                     (beta_1_2 / t55 * rs_s1rho / 0.2e1_dp + beta_2_2 * rs_s1rho + &
                     0.3e1_dp / 0.2e1_dp * beta_3_2 * t55 * rs_s1rho + t95 * t77 * &
                     rs_s1rho / rs_s1) / t100
                t536 = t194 ** 2
                A_s1rho = t385 / t536 * e_c_u_1_s1rho * t191 * t193
                t541 = t187 * t_s1
                t542 = t199 * t205
                t546 = A_s1rho * t197
                t547 = A_s1 * t_s1
                t549 = 0.2e1_dp * t547 * t_s1rho
                t554 = t187 * t197
                t555 = t204 ** 2
                t557 = t199 / t555
                t561 = t197 * t_s1
                t562 = t201 * t561
                t569 = 0.1e1_dp / t208
                t571 = epsilon_cGGA .LT. epsilon_cGGA_1_0
                IF (t571) THEN
                   marho = e_c_u_1_s1rho + t196 * (0.2e1_dp * t541 * t542 &
                        * t_s1rho + t187 * t197 * (t546 + t549) * t205 - t554 * t557 * (t546&
                        + t549 + 0.2e1_dp * A_s1 * t202 * A_s1rho + 0.4e1_dp * t562 * &
                        t_s1rho)) * t569
                ELSE
                   marho=epsilon_cGGArho
                END IF
                rs_s2rho = rs_s1rho
                t574 = k_s_s2 ** 2
                t_s2rho = -t28 / t574 * t7 / t26 * k_f_s1rho * t6 / 0.2e1_dp - t28&
                     * t29 * t343 / 0.2e1_dp
                t585 = t80 ** 2
                e_c_u_1_s2rho = -0.2e1_dp * t513 * rs_s2rho * t85 + t71 / t585 * (&
                     beta_1_2 / t36 * rs_s2rho / 0.2e1_dp + beta_2_2 * rs_s2rho + &
                     0.3e1_dp / 0.2e1_dp * beta_3_2 * t36 * rs_s2rho + t79 * t77 * &
                     rs_s2rho / rs_s2) / t84
                t604 = t216 ** 2
                A_s2rho = t385 / t604 * e_c_u_1_s2rho * t213 * t215
                t609 = t187 * t_s2
                t610 = t221 * t227
                t614 = A_s2rho * t219
                t615 = A_s2 * t_s2
                t617 = 0.2e1_dp * t615 * t_s2rho
                t622 = t187 * t219
                t623 = t226 ** 2
                t625 = t221 / t623
                t629 = t219 * t_s2
                t630 = t223 * t629
                t637 = 0.1e1_dp / t230
                t639 = epsilon_cGGA .LT. epsilon_cGGA_0_1
                IF (t639) THEN
                   mbrho = e_c_u_1_s2rho + t218 * (0.2e1_dp * t609 * t610 &
                        * t_s2rho + t187 * t219 * (t614 + t617) * t227 - t622 * t625 * (t614&
                        + t617 + 0.2e1_dp * A_s2 * t224 * A_s2rho + 0.4e1_dp * t630 * &
                        t_s2rho)) * t637
                ELSE
                   mbrho = epsilon_cGGArho
                END IF
                t642 = epsilon_cGGA * tau_w
                t645 = tau_w * t260
                epsilon_cRevPKZBrho = epsilon_cGGArho + (0.53e0_dp * &
                     epsilon_cGGArho * t256 + 0.106e1_dp * t642 * tau_wrho - 0.306e1_dp *&
                     t645 * tau_wrho - 0.153e1_dp * t256 * (marho / 0.2e1_dp + mbrho / &
                     0.2e1_dp)) * t265
                t659 = t256 * t271

                IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
                   e_rho(ii)=e_rho(ii)+&
                        scale_ec * (epsilon_cRevPKZB * t274 + my_rho * &
                        epsilon_cRevPKZBrho * t274 + t267 * (d * epsilon_cRevPKZBrho * t272 &
                        + 0.3e1_dp * t268 * t659 * tau_wrho)) + scale_ex * (ex_unif * Fx - &
                        my_rho * pi * t1 * t371 * Fx / 0.4e1_dp + t337 * &
                        t490 * (((0.2e1_dp * t445 * t303 * zrho - 0.4e1_dp * t450 * t452 * &
                        zrho) * p + t305 * prho + 0.292e3_dp / 0.2025e4_dp * tildeq_b * &
                        tildeq_brho - 0.73e2_dp / 0.4050e4_dp * tildeq_brho * t313 - &
                        0.73e2_dp / 0.8100e4_dp * t464 * (0.36e2_dp * z * zrho + 0.100e3_dp &
                        * p * prho) + 0.200e3_dp / 0.6561e4_dp * t472 * prho + 0.8e1_dp / &
                        0.45e2_dp * t475 * zrho + 0.3e1_dp * t322 * t310 * prho) * t329 - &
                        0.2e1_dp * t485 * t319 * prho))
                END IF

                tnorm_drho = t234 * t7 / 0.2e1_dp
                Hnorm_drho = gamma_var * (0.2e1_dp * t390 * t391 * tnorm_drho + &
                     0.2e1_dp * t187 * t410 * A * tnorm_drho * t249 - t403 * t406 * (&
                     0.2e1_dp * t396 * tnorm_drho + 0.4e1_dp * t411 * tnorm_drho)) * t419
                tau_wnorm_drho = my_ndrho * t7 / 0.4e1_dp
                pnorm_drho = my_ndrho * t1 * t284 / 0.6e1_dp
                znorm_drho = tau_wnorm_drho * t286
                alphanorm_drho = 0.5e1_dp / 0.3e1_dp * pnorm_drho * t288 - &
                     0.5e1_dp / 0.3e1_dp * t430 * znorm_drho
                tildeq_bnorm_drho = 0.9e1_dp / 0.20e2_dp * alphanorm_drho * t295 -&
                     0.9e1_dp / 0.40e2_dp * t437 * (b * alphanorm_drho * t290 + t291 * &
                     alphanorm_drho) + 0.2e1_dp / 0.3e1_dp * pnorm_drho
                t_s1norm_drho = t20 * t22 * t7 / 0.2e1_dp
                IF (t571) THEN
                   manorm_drho = t196 * (0.2e1_dp * t541 * t542 * &
                        t_s1norm_drho + 0.2e1_dp * t187 * t561 * A_s1 * t_s1norm_drho * t205&
                        - t554 * t557 * (0.2e1_dp * t547 * t_s1norm_drho + 0.4e1_dp * t562 &
                        * t_s1norm_drho)) * t569
                ELSE
                   manorm_drho= Hnorm_drho
                END IF
                t_s2norm_drho = t27 * t29 * t7 / 0.2e1_dp
                IF (t639) THEN
                   mbnorm_drho = t218 * (0.2e1_dp * t609 * t610 * &
                        t_s2norm_drho + 0.2e1_dp * t187 * t629 * A_s2 * t_s2norm_drho * t227&
                        - t622 * t625 * (0.2e1_dp * t615 * t_s2norm_drho + 0.4e1_dp * t630 &
                        * t_s2norm_drho)) * t637
                   mbnorm_drho = Hnorm_drho
                END IF
                epsilon_cRevPKZBnorm_drho = Hnorm_drho + (0.53e0_dp * Hnorm_drho *&
                     t256 + 0.106e1_dp * t642 * tau_wnorm_drho - 0.306e1_dp * t645 * &
                     tau_wnorm_drho - 0.153e1_dp * t256 * (manorm_drho / 0.2e1_dp + &
                     mbnorm_drho / 0.2e1_dp)) * t265

                IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
                   e_ndrho(ii)=e_ndrho(ii)+&
                        scale_ec * (my_rho * epsilon_cRevPKZBnorm_drho * &
                        t274 + t267 * (d * epsilon_cRevPKZBnorm_drho * t272 + 0.3e1_dp * &
                        t268 * t659 * tau_wnorm_drho)) + scale_ex * t337 * t490 * (((&
                        0.2e1_dp * t445 * t303 * znorm_drho - 0.4e1_dp * t450 * t452 * &
                        znorm_drho) * p + t305 * pnorm_drho + 0.292e3_dp / 0.2025e4_dp * &
                        tildeq_b * tildeq_bnorm_drho - 0.73e2_dp / 0.4050e4_dp * &
                        tildeq_bnorm_drho * t313 - 0.73e2_dp / 0.8100e4_dp * t464 * (&
                        0.36e2_dp * z * znorm_drho + 0.100e3_dp * p * pnorm_drho) + &
                        0.200e3_dp / 0.6561e4_dp * t472 * pnorm_drho + 0.8e1_dp / 0.45e2_dp &
                        * t475 * znorm_drho + 0.3e1_dp * t322 * t310 * pnorm_drho) * t329 - &
                        0.2e1_dp * t485 * t319 * pnorm_drho)
                END IF

                epsilon_cRevPKZBtau = -0.2e1_dp * t263 * t271
                t799 = t264 ** 2
                ztau = -tau_w * t265
                alphatau = -0.5e1_dp / 0.3e1_dp * t430 * ztau
                tildeq_btau = 0.9e1_dp / 0.20e2_dp * alphatau * t295 - 0.9e1_dp / &
                     0.40e2_dp * t437 * (b * alphatau * t290 + t291 * alphatau)

                IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
                   e_tau(ii)=e_tau(ii)+&
                        scale_ec * (my_rho * epsilon_cRevPKZBtau * t274 + t267 * &
                        (d * epsilon_cRevPKZBtau * t272 - 0.3e1_dp * t268 * t269 / t799)) + &
                        scale_ex * t337 * t490 * ((0.2e1_dp * t445 * t303 * ztau - 0.4e1_dp &
                        * t450 * t452 * ztau) * p + 0.292e3_dp / 0.2025e4_dp * tildeq_b * &
                        tildeq_btau - 0.73e2_dp / 0.4050e4_dp * tildeq_btau * t313 - &
                        0.73e2_dp / 0.225e3_dp * t464 * z * ztau + 0.8e1_dp / 0.45e2_dp * &
                        t475 * ztau) * t329
                END IF
             END IF
          END IF
       END DO

       !$omp end do

    END IF
  END SUBROUTINE tpss_lda_calc

! *****************************************************************************
!> \brief low level routine that calculates the tpss functional in the
!>      polarized (lsd) case
!> \param rhoa ...
!> \param rhob ...
!> \param norm_drhoa ...
!> \param norm_drhob ...
!> \param norm_drho ...
!> \param tau_a ...
!> \param tau_b ...
!> \param e_0 ...
!> \param e_rhoa ...
!> \param e_rhob ...
!> \param e_ndrhoa ...
!> \param e_ndrhob ...
!> \param e_ndrho ...
!> \param e_tau_a ...
!> \param e_tau_b ...
!> \param npoints ...
!> \param grad_deriv ...
!> \param epsilon_rho ...
!> \param epsilon_norm_drho ...
!> \param epsilon_tau ...
!> \param scale_ec ...
!> \param scale_ex ...
!> \param non_coer ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      maple is nice, but if you want the uman readable version of the code
!>      look in doc/tpss.mw
! *****************************************************************************
  SUBROUTINE tpss_lsd_calc(rhoa, rhob, norm_drhoa, norm_drhob, norm_drho,&
       tau_a,tau_b,e_0,e_rhoa,e_rhob,e_ndrhoa,e_ndrhob,e_ndrho,e_tau_a,e_tau_b,&
       npoints,grad_deriv,epsilon_rho,epsilon_norm_drho,epsilon_tau,&
       scale_ec,scale_ex,non_coer,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drhoa, &
                                                norm_drhob, norm_drho, tau_a, &
                                                tau_b
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rhoa, e_rhob, &
                                                e_ndrhoa, e_ndrhob, e_ndrho, &
                                                e_tau_a, e_tau_b
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, &
                                                epsilon_tau, scale_ec, &
                                                scale_ex
    INTEGER, INTENT(inout)                   :: non_coer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: abs_grad_deriv, ii
    LOGICAL                                  :: failure, t784, t863
    REAL(kind=dp) :: A, A_1, A_2, A_3, A_s1, A_s1rhoa, A_s2, A_s2rhob, &
      alpha_1_1, alpha_1_2, alpha_1_3, alpha_c, alpha_sp1, &
      alpha_sp1norm_drhoa, alpha_sp1rhoa, alpha_sp1tau_a, alpha_sp2, &
      alpha_sp2norm_drhob, alpha_sp2rhob, alpha_sp2tau_b, Arhoa, Arhob, b, &
      beta, beta_1_1, beta_1_2, beta_1_3, beta_2_1, beta_2_2, beta_2_3, &
      beta_3_1, beta_3_2, beta_3_3, beta_4_1, beta_4_2, beta_4_3, c, C_chi, &
      C_chi_eps, C_chi_epsnorm_drho, C_chi_epsnorm_drhoa, &
      C_chi_epsnorm_drhob, C_chi_epsrhoa, C_chi_epsrhob, chi, chirhoa, &
      chirhob, d, e_c_u_0, e_c_u_0rhoa, e_c_u_0rhob, e_c_u_1_s1, &
      e_c_u_1_s1rhoa, e_c_u_1_s2
    REAL(kind=dp) :: e_c_u_1_s2rhob, e_var, eps, epsilon_c_unif, &
      epsilon_c_unifrhoa, epsilon_c_unifrhob, epsilon_cGGA, epsilon_cGGA_0_1, &
      epsilon_cGGA_1_0, epsilon_cGGArhoa, epsilon_cGGArhob, epsilon_cRevPKZB, &
      epsilon_cRevPKZBnorm_drho, epsilon_cRevPKZBnorm_drhoa, &
      epsilon_cRevPKZBnorm_drhob, epsilon_cRevPKZBrhoa, epsilon_cRevPKZBrhob, &
      epsilon_cRevPKZBtau_a, epsilon_cRevPKZBtau_b, err_ndrho, err_ndrhoa, &
      err_ndrhob, err_rhoa, err_rhob, err_taua, err_taub, ex_unif_sp1, &
      ex_unif_sp2, f, f_ii_0, frhoa, frhob, Fx_sp1, Fx_sp2, gamma_var, &
      Hnorm_drho, k_frhoa, k_s, k_s_s1, k_s_s2, kappa, ma, manorm_drho
    REAL(kind=dp) :: manorm_drhoa, marhoa, marhob, mb, mbnorm_drho, &
      mbnorm_drhob, mbrhoa, mbrhob, mean_tmp, mu, my_norm_drho, &
      my_norm_drhoa, my_norm_drhob, my_rhoa, my_rhob, my_tau_a, my_tau_b, &
      p_1, p_2, p_3, p_sp1, p_sp1norm_drhoa, p_sp1rhoa, p_sp2, &
      p_sp2norm_drhob, p_sp2rhob, phi, phi_s1, phi_s2, phirhoa, phirhob, rho, &
      rs, rs_s1, rs_s1rhoa, rs_s2, rs_s2rhob, rsrhoa, rsrhob, t, t1, t101, &
      t102, t103, t1043, t107, t1070, t1072, t109, t110, t111, t1111, t1121, &
      t113, t1132, t115, t1151, t1156, t1157, t1161, t1162, t1164, t1169, &
      t1170, t1172, t1176, t1177, t118, t1184, t12, t121, t13, t130, t1303
    REAL(kind=dp) :: t132, t133, t134, t136, t14, t140, t144, t1457, t146, &
      t147, t148, t1496, t15, t150, t1500, t1506, t152, t154, t155, t156, &
      t157, t16, t161, t164, t165, t166, t167, t173, t174, t177, t179, t18, &
      t180, t187, t19, t193, t196, t198, t2, t20, t205, t209, t21, t211, &
      t215, t216, t217, t218, t22, t222, t223, t227, t232, t233, t234, t238, &
      t239, t245, t249, t25, t250, t257, t265, t27, t272, t274, t275, t277, &
      t280, t282, t285, t287, t288, t289, t293, t294, t298, t3, t30, t303, &
      t304, t305, t309, t310, t314, t319, t32, t320, t321, t325, t326, t328, &
      t329, t33, t330, t331, t332, t333, t337
    REAL(kind=dp) :: t338, t339, t340, t341, t343, t345, t346, t348, t349, &
      t35, t351, t352, t353, t355, t356, t358, t359, t36, t360, t361, t363, &
      t364, t366, t367, t37, t370, t371, t373, t374, t375, t377, t378, t380, &
      t381, t382, t383, t385, t386, t388, t389, t392, t393, t394, t395, t397, &
      t398, t400, t401, t402, t403, t404, t406, t407, t408, t409, t41, t411, &
      t412, t414, t415, t416, t417, t419, t420, t422, t423, t426, t427, t43, &
      t430, t431, t436, t436a, t437, t438, t44, t440, t441, t442, t449, t45, &
      t452, t455, t456, t457, t458, t459, t460, t462, t463, t464, t466, t468, &
      t47, t470, t471, t472, t474
    REAL(kind=dp) :: t475, t476, t478, t479, t481, t49, t492, t498, t5, t505, &
      t513, t518, t52, t520, t532, t540, t543, t55, t553, t557, t558, t559, &
      t560, t561, t570, t571, t577, t58, t581, t583, t584, t585, t589, t593, &
      t597, t6, t600, t604, t606, t607, t61, t611, t618, t621, t625, t627, &
      t628, t632, t639, t64, t650, t651, t660, t663, t665, t668, t67, t672, &
      t675, t678, t680, t684, t686, t687, t69, t690, t691, t699, t7, t70, &
      t703, t704, t708, t709, t71, t711, t716, t717, t719, t720, t723, t724, &
      t73, t731, t735, t736, t739, t746, t754, t758, t76, t766, t767, t769, &
      t77, t771, t773, t775, t785, t79
    REAL(kind=dp) :: t795, t806, t825, t827, t828, t83, t833, t834, t838, &
      t839, t841, t846, t847, t849, t85, t853, t854, t86, t861, t87, t873, &
      t879, t88, t882, t888, t893, t9, t908, t91, t914, t921, t929, t93, &
      t934, t936, t948, t956, t959, t96, t969, t973, t974, t975, t976, t98, &
      t99, t_s1, t_s1norm_drhoa, t_s1rhoa, t_s2, t_s2norm_drhob, t_s2rhob, &
      tau, tau_w, tau_w_sp1, tau_w_sp2, tau_wnorm_drho, tau_wrhoa, tau_wrhob, &
      tildeq_b_sp1, tildeq_b_sp1norm_drhoa, tildeq_b_sp1rhoa, &
      tildeq_b_sp1tau_a, tildeq_b_sp2, tildeq_b_sp2norm_drhob, &
      tildeq_b_sp2rhob, tildeq_b_sp2tau_b, tnorm_drho, trhoa, trhob, z_sp1
    REAL(kind=dp) :: z_sp1norm_drhoa, z_sp1rhoa, z_sp1tau_a, z_sp2, &
      z_sp2norm_drhob, z_sp2rhob, z_sp2tau_b

    err_rhoa=0._dp
    err_rhob=0._dp
    err_ndrho=0._dp
    err_ndrhoa=0._dp
    err_ndrhob=0._dp
    err_taua=0._dp
    err_taub=0._dp

    failure=.FALSE.

    kappa = 0.804e0_dp
    beta = 0.66725e-1_dp
    mu = 0.21951e0_dp
    gamma_var = (0.1e1_dp - LOG(0.2e1_dp)) / pi ** 2
    b = 0.4e0_dp
    c = 0.159096e1_dp
    e_var = 0.1537e1_dp
    d = 0.28e1_dp
    f_ii_0 = 0.8e1_dp / 0.9e1_dp / (2._dp * 2._dp ** (0.1e1_dp / 0.3e1_dp) - 2._dp)
    p_1 = 0.10e1_dp
    A_1 = 0.31091e-1_dp
    alpha_1_1 = 0.21370e0_dp
    beta_1_1 = 0.75957e1_dp
    beta_2_1 = 0.35876e1_dp
    beta_3_1 = 0.16382e1_dp
    beta_4_1 = 0.49294e0_dp
    p_2 = 0.10e1_dp
    A_2 = 0.15545e-1_dp
    alpha_1_2 = 0.20548e0_dp
    beta_1_2 = 0.141189e2_dp
    beta_2_2 = 0.61977e1_dp
    beta_3_2 = 0.33662e1_dp
    beta_4_2 = 0.62517e0_dp
    p_3 = 0.10e1_dp
    A_3 = 0.16887e-1_dp
    alpha_1_3 = 0.11125e0_dp
    beta_1_3 = 0.10357e2_dp
    beta_2_3 = 0.36231e1_dp
    beta_3_3 = 0.88026e0_dp
    beta_4_3 = 0.49671e0_dp

    t12 = 3._dp ** (0.1e1_dp / 0.3e1_dp)
    t14 = pi ** 2
    t15 = t14 ** (0.1e1_dp / 0.3e1_dp)
    t16 = t15 ** 2
    t18 = 2._dp ** (0.1e1_dp / 0.3e1_dp)
    t19 = 1._dp / t16 * t18

    IF (.NOT. failure) THEN

       abs_grad_deriv=ABS(grad_deriv)

       !$omp do

       DO ii=1,npoints
          my_tau_a=MAX(0.0_dp,tau_a(ii))
          my_tau_b=MAX(0.0_dp,tau_b(ii))
          my_rhoa=MAX(0.0_dp,rhoa(ii))
          my_rhob=MAX(0.0_dp,rhob(ii))
          rho=my_rhoa+my_rhob
          tau=my_tau_a+my_tau_b
          IF (rho>epsilon_rho.AND.tau>epsilon_tau) THEN
             my_norm_drhoa=MAX(EPSILON(0.0_dp),norm_drhoa(ii))
             my_norm_drhob=MAX(EPSILON(0.0_dp),norm_drhob(ii))
             my_norm_drho=MAX(EPSILON(0.0_dp),norm_drho(ii))
             my_rhoa=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhoa)
             my_rhob=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhob)
             my_tau_a=MAX(EPSILON(0.0_dp)*1.e4_dp,tau_a(ii))
             my_tau_b=MAX(EPSILON(0.0_dp)*1.e4_dp,tau_b(ii))

             ! enforce 8*rho*tau >= norm_drho^2 (also each spin separately)

             t1 = my_norm_drhoa ** 2
             mean_tmp=8._dp*my_rhoa*my_tau_a
             IF (mean_tmp<t1) THEN
                mean_tmp=0.5_dp*(mean_tmp+t1)
                t1=mean_tmp
                my_norm_drhoa=SQRT(mean_tmp)
                my_tau_a=mean_tmp/(8._dp*my_rhoa)
                non_coer=non_coer+1
             END IF
             t5 = my_norm_drhob ** 2
             mean_tmp=8._dp*my_rhob*my_tau_b
             IF (mean_tmp<t5) THEN
                mean_tmp=0.5_dp*(mean_tmp+t5)
                t5=mean_tmp
                my_norm_drhob=SQRT(mean_tmp)
                my_tau_b=mean_tmp/(8._dp*my_rhob)
                non_coer=non_coer+1
             END IF

             rho=my_rhoa+my_rhob
             t9 = my_tau_a + my_tau_b
             tau = t9
             my_norm_drho=MIN(my_norm_drho,my_norm_drhoa+my_norm_drhob)

             t2 = 1._dp / my_rhoa
             t3 = t1 * t2
             t6 = 1._dp / my_rhob
             t7 = t5 * t6
             t13 = t1 * t12
             t20 = my_rhoa ** 2
             t21 = my_rhoa ** (0.1e1_dp / 0.3e1_dp)
             t22 = t21 ** 2
             t25 = t19 / t22 / t20
             p_sp1 = t13 * t25 / 0.24e2_dp
             tau_w_sp1 = t3 / 0.4e1_dp
             t27 = 0.1e1_dp / my_tau_a
             z_sp1 = tau_w_sp1 * t27 / 0.2e1_dp
             t30 = 0.1e1_dp / z_sp1 - 0.1e1_dp
             alpha_sp1 = 0.5e1_dp / 0.3e1_dp * p_sp1 * t30
             t32 = alpha_sp1 - 0.1e1_dp
             t33 = b * alpha_sp1
             t35 = 0.1e1_dp + t33 * t32
             t36 = SQRT(t35)
             t37 = 0.1e1_dp / t36
             tildeq_b_sp1 = 0.9e1_dp / 0.20e2_dp * t32 * t37 + 0.2e1_dp / &
                  0.3e1_dp * p_sp1
             t41 = z_sp1 ** 2
             t43 = 0.1e1_dp + t41
             t44 = t43 ** 2
             t45 = 0.1e1_dp / t44
             t47 = 0.10e2_dp / 0.81e2_dp + c * t41 * t45
             t49 = tildeq_b_sp1 ** 2
             t52 = p_sp1 ** 2
             t55 = SQRT(0.18e2_dp * t41 + 0.50e2_dp * t52)
             t58 = 0.1e1_dp / kappa
             t61 = SQRT(e_var)
             t64 = e_var * mu
             t67 = t47 * p_sp1 + 0.146e3_dp / 0.2025e4_dp * t49 - 0.73e2_dp / &
                  0.4050e4_dp * tildeq_b_sp1 * t55 + 0.100e3_dp / 0.6561e4_dp * t58 * &
                  t52 + 0.4e1_dp / 0.45e2_dp * t61 * t41 + t64 * t52 * p_sp1
             t69 = 0.1e1_dp + t61 * p_sp1
             t70 = t69 ** 2
             t71 = 0.1e1_dp / t70
             t73 = 0.1e1_dp + t67 * t71 * t58
             Fx_sp1 = 0.1e1_dp + kappa - kappa / t73
             t76 = 1._dp / pi
             t77 = t76 * t12
             t79 = (t14 * my_rhoa) ** (0.1e1_dp / 0.3e1_dp)
             ex_unif_sp1 = -0.3e1_dp / 0.4e1_dp * t77 * t18 * t79
             t83 = my_rhoa * ex_unif_sp1
             t85 = t5 * t12
             t86 = my_rhob ** 2
             t87 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
             t88 = t87 ** 2
             t91 = t19 / t88 / t86
             p_sp2 = t85 * t91 / 0.24e2_dp
             tau_w_sp2 = t7 / 0.4e1_dp
             t93 = 0.1e1_dp / my_tau_b
             z_sp2 = tau_w_sp2 * t93 / 0.2e1_dp
             t96 = 0.1e1_dp / z_sp2 - 0.1e1_dp
             alpha_sp2 = 0.5e1_dp / 0.3e1_dp * p_sp2 * t96
             t98 = alpha_sp2 - 0.1e1_dp
             t99 = b * alpha_sp2
             t101 = 0.1e1_dp + t99 * t98
             t102 = SQRT(t101)
             t103 = 0.1e1_dp / t102
             tildeq_b_sp2 = 0.9e1_dp / 0.20e2_dp * t98 * t103 + 0.2e1_dp / &
                  0.3e1_dp * p_sp2
             t107 = z_sp2 ** 2
             t109 = 0.1e1_dp + t107
             t110 = t109 ** 2
             t111 = 0.1e1_dp / t110
             t113 = 0.10e2_dp / 0.81e2_dp + c * t107 * t111
             t115 = tildeq_b_sp2 ** 2
             t118 = p_sp2 ** 2
             t121 = SQRT(0.18e2_dp * t107 + 0.50e2_dp * t118)
             t130 = t113 * p_sp2 + 0.146e3_dp / 0.2025e4_dp * t115 - 0.73e2_dp &
                  / 0.4050e4_dp * tildeq_b_sp2 * t121 + 0.100e3_dp / 0.6561e4_dp * t58&
                  * t118 + 0.4e1_dp / 0.45e2_dp * t61 * t107 + t64 * t118 * p_sp2
             t132 = 0.1e1_dp + t61 * p_sp2
             t133 = t132 ** 2
             t134 = 0.1e1_dp / t133
             t136 = 0.1e1_dp + t130 * t134 * t58
             Fx_sp2 = 0.1e1_dp + kappa - kappa / t136
             t140 = (t14 * my_rhob) ** (0.1e1_dp / 0.3e1_dp)
             ex_unif_sp2 = -0.3e1_dp / 0.4e1_dp * t77 * t18 * t140
             t144 = my_rhob * ex_unif_sp2
             t146 = 4._dp ** (0.1e1_dp / 0.3e1_dp)
             t147 = t146 ** 2
             t148 = t12 * t147
             t150 = (t76 * t2) ** (0.1e1_dp / 0.3e1_dp)
             rs_s1 = t148 * t150 / 0.4e1_dp
             t152 = t18 ** 2
             phi_s1 = t152 / 0.2e1_dp
             t154 = SQRT(t12 * t79 * t76)
             k_s_s1 = 0.2e1_dp * t154
             t155 = 0.1e1_dp / phi_s1
             t156 = my_norm_drhoa * t155
             t157 = 0.1e1_dp / k_s_s1
             t_s1 = t156 * t157 * t2 / 0.2e1_dp
             t161 = (t76 * t6) ** (0.1e1_dp / 0.3e1_dp)
             rs_s2 = t148 * t161 / 0.4e1_dp
             phi_s2 = phi_s1
             t164 = SQRT(t12 * t140 * t76)
             k_s_s2 = 0.2e1_dp * t164
             t165 = 0.1e1_dp / phi_s2
             t166 = my_norm_drhob * t165
             t167 = 0.1e1_dp / k_s_s2
             t_s2 = t166 * t167 * t6 / 0.2e1_dp
             t173 = 0.1e1_dp / A_1
             t174 = SQRT(rs_s2)
             t177 = t174 * rs_s2
             t179 = p_1 + 0.1e1_dp
             t180 = rs_s2 ** t179
             t187 = LOG(0.1e1_dp + t173 / (beta_1_1 * t174 + beta_2_1 * rs_s2 +&
                  beta_3_1 * t177 + beta_4_1 * t180) / 0.2e1_dp)
             t193 = SQRT(rs_s1)
             t196 = t193 * rs_s1
             t198 = rs_s1 ** t179
             t205 = LOG(0.1e1_dp + t173 / (beta_1_1 * t193 + beta_2_1 * rs_s1 +&
                  beta_3_1 * t196 + beta_4_1 * t198) / 0.2e1_dp)
             t209 = 0.1e1_dp + alpha_1_2 * rs_s2
             t211 = 0.1e1_dp / A_2
             t215 = p_2 + 0.1e1_dp
             t216 = rs_s2 ** t215
             t217 = beta_4_2 * t216
             t218 = beta_1_2 * t174 + beta_2_2 * rs_s2 + beta_3_2 * t177 + t217
             t222 = 0.1e1_dp + t211 / t218 / 0.2e1_dp
             t223 = LOG(t222)
             e_c_u_1_s2 = -0.2e1_dp * A_2 * t209 * t223
             t227 = 0.1e1_dp + alpha_1_2 * rs_s1
             t232 = rs_s1 ** t215
             t233 = beta_4_2 * t232
             t234 = beta_1_2 * t193 + beta_2_2 * rs_s1 + beta_3_2 * t196 + t233
             t238 = 0.1e1_dp + t211 / t234 / 0.2e1_dp
             t239 = LOG(t238)
             e_c_u_1_s1 = -0.2e1_dp * A_2 * t227 * t239
             t245 = 0.1e1_dp / A_3
             t249 = p_3 + 0.1e1_dp
             t250 = rs_s2 ** t249
             t257 = LOG(0.1e1_dp + t245 / (beta_1_3 * t174 + beta_2_3 * rs_s2 +&
                  beta_3_3 * t177 + beta_4_3 * t250) / 0.2e1_dp)
             t265 = rs_s1 ** t249
             t272 = LOG(0.1e1_dp + t245 / (beta_1_3 * t193 + beta_2_3 * rs_s1 +&
                  beta_3_3 * t196 + beta_4_3 * t265) / 0.2e1_dp)
             t274 = my_rhoa - my_rhob
             t275 = 0.1e1_dp / rho
             chi = t274 * t275
             t277 = (t76 * t275) ** (0.1e1_dp / 0.3e1_dp)
             rs = t148 * t277 / 0.4e1_dp
             t280 = 0.1e1_dp + alpha_1_1 * rs
             t282 = SQRT(rs)
             t285 = t282 * rs
             t287 = rs ** t179
             t288 = beta_4_1 * t287
             t289 = beta_1_1 * t282 + beta_2_1 * rs + beta_3_1 * t285 + t288
             t293 = 0.1e1_dp + t173 / t289 / 0.2e1_dp
             t294 = LOG(t293)
             e_c_u_0 = -0.2e1_dp * A_1 * t280 * t294
             t298 = 0.1e1_dp + alpha_1_2 * rs
             t303 = rs ** t215
             t304 = beta_4_2 * t303
             t305 = beta_1_2 * t282 + beta_2_2 * rs + beta_3_2 * t285 + t304
             t309 = 0.1e1_dp + t211 / t305 / 0.2e1_dp
             t310 = LOG(t309)
             t314 = 0.1e1_dp + alpha_1_3 * rs
             t319 = rs ** t249
             t320 = beta_4_3 * t319
             t321 = beta_1_3 * t282 + beta_2_3 * rs + beta_3_3 * t285 + t320
             t325 = 0.1e1_dp + t245 / t321 / 0.2e1_dp
             t326 = LOG(t325)
             alpha_c = 0.2e1_dp * A_3 * t314 * t326
             t328 = 0.1e1_dp + chi
             t329 = t328 ** (0.1e1_dp / 0.3e1_dp)
             t330 = t329 * t328
             t331 = 0.1e1_dp - chi
             t332 = t331 ** (0.1e1_dp / 0.3e1_dp)
             t333 = t332 * t331
             t337 = 1._dp / (2 * t18 - 2._dp)
             f = (t330 + t333 - 0.2e1_dp) * t337
             t338 = alpha_c * f
             t339 = 0.1e1_dp / f_ii_0
             t340 = chi ** 2
             t341 = t340 ** 2
             t343 = t339 * (0.1e1_dp - t341)
             t345 = -0.2e1_dp * A_2 * t298 * t310 - e_c_u_0
             t346 = t345 * f
             epsilon_c_unif = e_c_u_0 + t338 * t343 + t346 * t341
             t348 = 0.1e1_dp / gamma_var
             t349 = beta * t348
             t351 = phi_s1 ** 2
             t352 = t351 * phi_s1
             t353 = 0.1e1_dp / t352
             t355 = EXP(-e_c_u_1_s1 * t348 * t353)
             t356 = t355 - 0.1e1_dp
             A_s1 = t349 / t356
             t358 = gamma_var * t352
             t359 = t_s1 ** 2
             t360 = A_s1 * t359
             t361 = 0.1e1_dp + t360
             t363 = A_s1 ** 2
             t364 = t359 ** 2
             t366 = 0.1e1_dp + t360 + t363 * t364
             t367 = 0.1e1_dp / t366
             t370 = 0.1e1_dp + t349 * t359 * t361 * t367
             t371 = LOG(t370)
             epsilon_cGGA_1_0 = e_c_u_1_s1 + t358 * t371
             t373 = phi_s2 ** 2
             t374 = t373 * phi_s2
             t375 = 0.1e1_dp / t374
             t377 = EXP(-e_c_u_1_s2 * t348 * t375)
             t378 = t377 - 0.1e1_dp
             A_s2 = t349 / t378
             t380 = gamma_var * t374
             t381 = t_s2 ** 2
             t382 = A_s2 * t381
             t383 = 0.1e1_dp + t382
             t385 = A_s2 ** 2
             t386 = t381 ** 2
             t388 = 0.1e1_dp + t382 + t385 * t386
             t389 = 0.1e1_dp / t388
             t392 = 0.1e1_dp + t349 * t381 * t383 * t389
             t393 = LOG(t392)
             epsilon_cGGA_0_1 = e_c_u_1_s2 + t380 * t393
             t394 = t329 ** 2
             t395 = t332 ** 2
             phi = t394 / 0.2e1_dp + t395 / 0.2e1_dp
             t397 = t14 * rho
             t398 = t397 ** (0.1e1_dp / 0.3e1_dp)
             t400 = SQRT(t12 * t398 * t76)
             k_s = 0.2e1_dp * t400
             t401 = 0.1e1_dp / phi
             t402 = my_norm_drho * t401
             t403 = 0.1e1_dp / k_s
             t404 = t403 * t275
             t = t402 * t404 / 0.2e1_dp
             t406 = epsilon_c_unif * t348
             t407 = phi ** 2
             t408 = t407 * phi
             t409 = 0.1e1_dp / t408
             t411 = EXP(-t406 * t409)
             t412 = t411 - 0.1e1_dp
             A = t349 / t412
             t414 = gamma_var * t408
             t415 = t ** 2
             t416 = A * t415
             t417 = 0.1e1_dp + t416
             t419 = A ** 2
             t420 = t415 ** 2
             t422 = 0.1e1_dp + t416 + t419 * t420
             t423 = 0.1e1_dp / t422
             t426 = 0.1e1_dp + t349 * t415 * t417 * t423
             t427 = LOG(t426)
             epsilon_cGGA = epsilon_c_unif + t414 * t427
             t430 = my_rhoa * my_rhob
             t431 = my_norm_drho ** 2
             t436a = t1 * t86 + t5 * t20 - t430 * t431 + t430 *&
                  t1 + t430 * t5
             IF (t436a<0._dp) THEN
!                PRINT *, t436a
                t436=0._dp
             ELSE
                t436 = SQRT(t436a)
             END IF
             t437 = rho ** 2
             t438 = 0.1e1_dp / t437
             t440 = t12 ** 2
             t441 = 0.2e1_dp * t436 * t438 * t440
             t442 = 0.1e1_dp / t398
             eps = t441 * t442 / 0.6e1_dp
             C_chi = 0.53e0_dp + 0.87e0_dp * t340 + 0.5e0_dp * t341 + &
                  0.226e1_dp * t341 * t340
             tau_w = t431 * t275 / 0.8e1_dp
             t449 = eps ** 2
             t452 = 0.1e1_dp / t330 + 0.1e1_dp / t333
             t455 = 0.1e1_dp + t449 * t452 / 0.2e1_dp
             t456 = t455 ** 2
             t457 = t456 ** 2
             t458 = 0.1e1_dp / t457
             C_chi_eps = C_chi * t458
             ma = MAX(epsilon_cGGA_1_0, epsilon_cGGA)
             mb = MAX(epsilon_cGGA_0_1, epsilon_cGGA)
             t459 = epsilon_cGGA * C_chi_eps
             t460 = tau_w ** 2
             t462 = 0.1e1_dp + C_chi_eps
             t463 = t462 * t460
             t464 = my_rhoa * t275
             t466 = my_rhob * t275
             t468 = t464 * ma + t466 * mb
             t470 = t459 * t460 - t463 * t468
             t471 = tau ** 2
             t472 = 0.1e1_dp / t471
             epsilon_cRevPKZB = epsilon_cGGA + t470 * t472
             t474 = rho * epsilon_cRevPKZB
             t475 = d * epsilon_cRevPKZB
             t476 = t460 * tau_w
             t478 = 0.1e1_dp / t471 / tau
             t479 = t476 * t478
             t481 = 0.1e1_dp + t475 * t479

             IF (grad_deriv>=0) THEN
                e_0(ii)=e_0(ii)+&
                     scale_ex * (t83 * Fx_sp1 + t144 * Fx_sp2) + scale_ec * &
                     t474 * t481
             END IF

             IF (abs_grad_deriv>0) THEN
                p_sp1rhoa = -t13 * t19 / t22 / t20 / my_rhoa / 0.9e1_dp
                t492 = 0.1e1_dp / t20
                z_sp1rhoa = -t1 * t492 * t27 / 0.8e1_dp
                t498 = p_sp1 / t41
                alpha_sp1rhoa = 0.5e1_dp / 0.3e1_dp * p_sp1rhoa * t30 - 0.5e1_dp /&
                     0.3e1_dp * t498 * z_sp1rhoa
                t505 = t32 / t36 / t35
                tildeq_b_sp1rhoa = 0.9e1_dp / 0.20e2_dp * alpha_sp1rhoa * t37 - &
                     0.9e1_dp / 0.40e2_dp * t505 * (b * alpha_sp1rhoa * t32 + t33 * &
                     alpha_sp1rhoa) + 0.2e1_dp / 0.3e1_dp * p_sp1rhoa
                t513 = c * z_sp1
                t518 = c * t41 * z_sp1
                t520 = 0.1e1_dp / t44 / t43
                t532 = tildeq_b_sp1 / t55
                t540 = t58 * p_sp1
                t543 = t61 * z_sp1
                t553 = t67 / t70 / t69
                t557 = t73 ** 2
                t558 = 0.1e1_dp / t557
                t559 = pi * t12
                t560 = t79 ** 2
                t561 = 0.1e1_dp / t560
                t570 = t274 * t438
                chirhoa = t275 - t570
                t571 = t277 ** 2
                rsrhoa = -t148 / t571 * t76 * t438 / 0.12e2_dp
                t577 = A_1 * alpha_1_1
                t581 = t289 ** 2
                t583 = t280 / t581
                t584 = 0.1e1_dp / t282
                t585 = beta_1_1 * t584
                t589 = beta_3_1 * t282
                t593 = 0.1e1_dp / rs
                t597 = 0.1e1_dp / t293
                e_c_u_0rhoa = -0.2e1_dp * t577 * rsrhoa * t294 + t583 * (t585 * &
                     rsrhoa / 0.2e1_dp + beta_2_1 * rsrhoa + 0.3e1_dp / 0.2e1_dp * t589 *&
                     rsrhoa + t288 * t179 * rsrhoa * t593) * t597
                t600 = A_2 * alpha_1_2
                t604 = t305 ** 2
                t606 = t298 / t604
                t607 = beta_1_2 * t584
                t611 = beta_3_2 * t282
                t618 = 0.1e1_dp / t309
                t621 = A_3 * alpha_1_3
                t625 = t321 ** 2
                t627 = t314 / t625
                t628 = beta_1_3 * t584
                t632 = beta_3_3 * t282
                t639 = 0.1e1_dp / t325
                frhoa = (0.4e1_dp / 0.3e1_dp * t329 * chirhoa - 0.4e1_dp / &
                     0.3e1_dp * t332 * chirhoa) * t337
                t650 = t340 * chi
                t651 = t339 * t650
                t660 = t650 * chirhoa
                epsilon_c_unifrhoa = e_c_u_0rhoa + (0.2e1_dp * t621 * rsrhoa * &
                     t326 - t627 * (t628 * rsrhoa / 0.2e1_dp + beta_2_3 * rsrhoa + &
                     0.3e1_dp / 0.2e1_dp * t632 * rsrhoa + t320 * t249 * rsrhoa * t593) *&
                     t639) * f * t343 + alpha_c * frhoa * t343 - 0.4e1_dp * t338 * t651 &
                     * chirhoa + (-0.2e1_dp * t600 * rsrhoa * t310 + t606 * (t607 * &
                     rsrhoa / 0.2e1_dp + beta_2_2 * rsrhoa + 0.3e1_dp / 0.2e1_dp * t611 *&
                     rsrhoa + t304 * t215 * rsrhoa * t593) * t618 - e_c_u_0rhoa) * f * &
                     t341 + t345 * frhoa * t341 + 0.4e1_dp * t346 * t660
                t663 = 0.1e1_dp / t329
                t665 = 0.1e1_dp / t332
                phirhoa = t663 * chirhoa / 0.3e1_dp - t665 * chirhoa / 0.3e1_dp
                t668 = t398 ** 2
                k_frhoa = t12 / t668 * t14 / 0.3e1_dp
                t672 = 0.1e1_dp / t400
                t675 = my_norm_drho / t407
                t678 = k_s ** 2
                t680 = 0.1e1_dp / t678 * t275
                t684 = t402 * t403 * t438
                trhoa = -t675 * t404 * phirhoa / 0.2e1_dp - t402 * t680 * t672 * &
                     k_frhoa * t76 / 0.2e1_dp - t684 / 0.2e1_dp
                t686 = t412 ** 2
                t687 = 0.1e1_dp / t686
                t690 = t407 ** 2
                t691 = 0.1e1_dp / t690
                Arhoa = -t349 * t687 * (-epsilon_c_unifrhoa * t348 * t409 + &
                     0.3e1_dp * t406 * t691 * phirhoa) * t411
                t699 = gamma_var * t407
                t703 = t349 * t
                t704 = t417 * t423
                t708 = Arhoa * t415
                t709 = A * t
                t711 = 0.2e1_dp * t709 * trhoa
                t716 = t349 * t415
                t717 = t422 ** 2
                t719 = t417 / t717
                t720 = A * t420
                t723 = t415 * t
                t724 = t419 * t723
                t731 = 0.1e1_dp / t426
                epsilon_cGGArhoa = epsilon_c_unifrhoa + 0.3e1_dp * t699 * t427 * &
                     phirhoa + t414 * (0.2e1_dp * t703 * t704 * trhoa + t349 * t415 * (&
                     t708 + t711) * t423 - t716 * t719 * (t708 + t711 + 0.2e1_dp * t720 *&
                     Arhoa + 0.4e1_dp * t724 * trhoa)) * t731
                t735 = 0.1e1_dp / MAX(EPSILON(0._dp),t436 * t438)
                t736 = t5 * my_rhoa
                t739 = my_rhob * t1
                t746 = 0.4e1_dp * t436 / t437 / rho
                t754 = t441 / t398 / t397 * t14 / 0.18e2_dp
                t758 = t341 * chi
                tau_wrhoa = -t431 * t438 / 0.8e1_dp
                t766 = C_chi / t457 / t455
                t767 = eps * t452
                t769 = t328 ** 2
                t771 = 0.1e1_dp / t329 / t769
                t773 = t331 ** 2
                t775 = 0.1e1_dp / t332 / t773
                C_chi_epsrhoa = (0.174e1_dp * chi * chirhoa + 0.20e1_dp * t660 + &
                     0.1356e2_dp * t758 * chirhoa) * t458 - 0.4e1_dp * t766 * (t767 * ((&
                     t735 * (0.2e1_dp * t736 - my_rhob * t431 + t739 + my_rhob * t5) - &
                     t746) * t440 * t442 / 0.6e1_dp - t754) + t449 * (-0.4e1_dp / &
                     0.3e1_dp * t771 * chirhoa + 0.4e1_dp / 0.3e1_dp * t775 * chirhoa) / &
                     0.2e1_dp)
                t784 = epsilon_cGGA .LT. epsilon_cGGA_0_1
                IF (t784) THEN
                   mbrhoa=0._dp
                ELSE
                   mbrhoa =  epsilon_cGGArhoa
                END IF
                t785 = t150 ** 2
                rs_s1rhoa = -t148 / t785 * t76 * t492 / 0.12e2_dp
                t795 = k_s_s1 ** 2
                t_s1rhoa = -t156 / t795 * t2 / t154 * t12 * t561 * t14&
                     * t76 / 0.6e1_dp - t156 * t157 * t492 / 0.2e1_dp
                t806 = t234 ** 2
                e_c_u_1_s1rhoa = -0.2e1_dp * t600 * rs_s1rhoa * t239 + t227 / t806&
                     * (beta_1_2 / t193 * rs_s1rhoa / 0.2e1_dp + beta_2_2 * rs_s1rhoa + &
                     0.3e1_dp / 0.2e1_dp * beta_3_2 * t193 * rs_s1rhoa + t233 * t215 * &
                     rs_s1rhoa / rs_s1) / t238
                t825 = gamma_var ** 2
                t827 = beta / t825
                t828 = t356 ** 2
                A_s1rhoa = t827 / t828 * e_c_u_1_s1rhoa * t353 * t355
                t833 = t349 * t_s1
                t834 = t361 * t367
                t838 = A_s1rhoa * t359
                t839 = A_s1 * t_s1
                t841 = 0.2e1_dp * t839 * t_s1rhoa
                t846 = t349 * t359
                t847 = t366 ** 2
                t849 = t361 / t847
                t853 = t359 * t_s1
                t854 = t363 * t853
                t861 = 0.1e1_dp / t370
                t863 = epsilon_cGGA .LT. epsilon_cGGA_1_0
                IF (t863) THEN
                   marhoa =  e_c_u_1_s1rhoa + t358 * (0.2e1_dp * t833 * &
                        t834 * t_s1rhoa + t349 * t359 * (t838 + t841) * t367 - t846 * t849 *&
                        (t838 + t841 + 0.2e1_dp * A_s1 * t364 * A_s1rhoa + 0.4e1_dp * t854 &
                        * t_s1rhoa)) * t861
                ELSE
                   marhoa = epsilon_cGGArhoa
                END IF
                t873 = t462 * tau_w
                t879 = my_rhoa * t438 * ma
                t882 = my_rhob * t438 * mb
                epsilon_cRevPKZBrhoa = epsilon_cGGArhoa + (epsilon_cGGArhoa * &
                     C_chi_eps * t460 + epsilon_cGGA * C_chi_epsrhoa * t460 + 0.2e1_dp * &
                     t459 * tau_w * tau_wrhoa - C_chi_epsrhoa * t460 * t468 - 0.2e1_dp * &
                     t873 * t468 * tau_wrhoa - t463 * (t275 * ma - t879 + t464 * marhoa -&
                     t882 + t466 * mbrhoa)) * t472
                t888 = epsilon_cRevPKZB * t481
                t893 = t460 * t478

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_rhoa(ii)=e_rhoa(ii)+&
                        scale_ex * (0.2e1_dp * ex_unif_sp1 * Fx_sp1 - my_rhoa * &
                        t559 * t18 * t561 * Fx_sp1 / 0.2e1_dp + 0.2e1_dp * t83 *&
                        t558 * (((0.2e1_dp * t513 * t45 * z_sp1rhoa - 0.4e1_dp * t518 * &
                        t520 * z_sp1rhoa) * p_sp1 + t47 * p_sp1rhoa + 0.292e3_dp / &
                        0.2025e4_dp * tildeq_b_sp1 * tildeq_b_sp1rhoa - 0.73e2_dp / &
                        0.4050e4_dp * tildeq_b_sp1rhoa * t55 - 0.73e2_dp / 0.8100e4_dp * &
                        t532 * (0.36e2_dp * z_sp1 * z_sp1rhoa + 0.100e3_dp * p_sp1 * &
                        p_sp1rhoa) + 0.200e3_dp / 0.6561e4_dp * t540 * p_sp1rhoa + 0.8e1_dp &
                        / 0.45e2_dp * t543 * z_sp1rhoa + 0.3e1_dp * t64 * t52 * p_sp1rhoa) *&
                        t71 - 0.2e1_dp * t553 * t61 * p_sp1rhoa)) / 0.2e1_dp + scale_ec * (&
                        t888 + rho * epsilon_cRevPKZBrhoa * t481 + t474 * (d * &
                        epsilon_cRevPKZBrhoa * t479 + 0.3e1_dp * t475 * t893 * tau_wrhoa))
                END IF

                p_sp2rhob = -t85 * t19 / t88 / t86 / my_rhob / 0.9e1_dp
                t908 = 0.1e1_dp / t86
                z_sp2rhob = -t5 * t908 * t93 / 0.8e1_dp
                t914 = p_sp2 / t107
                alpha_sp2rhob = 0.5e1_dp / 0.3e1_dp * p_sp2rhob * t96 - 0.5e1_dp /&
                     0.3e1_dp * t914 * z_sp2rhob
                t921 = t98 / t102 / t101
                tildeq_b_sp2rhob = 0.9e1_dp / 0.20e2_dp * alpha_sp2rhob * t103 - &
                     0.9e1_dp / 0.40e2_dp * t921 * (b * alpha_sp2rhob * t98 + t99 * &
                     alpha_sp2rhob) + 0.2e1_dp / 0.3e1_dp * p_sp2rhob
                t929 = c * z_sp2
                t934 = c * t107 * z_sp2
                t936 = 0.1e1_dp / t110 / t109
                t948 = tildeq_b_sp2 / t121
                t956 = t58 * p_sp2
                t959 = t61 * z_sp2
                t969 = t130 / t133 / t132
                t973 = t136 ** 2
                t974 = 0.1e1_dp / t973
                t975 = t140 ** 2
                t976 = 0.1e1_dp / t975
                chirhob = -t275 - t570
                rsrhob = rsrhoa
                e_c_u_0rhob = -0.2e1_dp * t577 * rsrhob * t294 + t583 * (t585 * &
                     rsrhob / 0.2e1_dp + beta_2_1 * rsrhob + 0.3e1_dp / 0.2e1_dp * t589 *&
                     rsrhob + t288 * t179 * rsrhob * t593) * t597
                frhob = (0.4e1_dp / 0.3e1_dp * t329 * chirhob - 0.4e1_dp / &
                     0.3e1_dp * t332 * chirhob) * t337
                t1043 = t650 * chirhob
                epsilon_c_unifrhob = e_c_u_0rhob + (0.2e1_dp * t621 * rsrhob * &
                     t326 - t627 * (t628 * rsrhob / 0.2e1_dp + beta_2_3 * rsrhob + &
                     0.3e1_dp / 0.2e1_dp * t632 * rsrhob + t320 * t249 * rsrhob * t593) *&
                     t639) * f * t343 + alpha_c * frhob * t343 - 0.4e1_dp * t338 * t651 &
                     * chirhob + (-0.2e1_dp * t600 * rsrhob * t310 + t606 * (t607 * &
                     rsrhob / 0.2e1_dp + beta_2_2 * rsrhob + 0.3e1_dp / 0.2e1_dp * t611 *&
                     rsrhob + t304 * t215 * rsrhob * t593) * t618 - e_c_u_0rhob) * f * &
                     t341 + t345 * frhob * t341 + 0.4e1_dp * t346 * t1043
                phirhob = t663 * chirhob / 0.3e1_dp - t665 * chirhob / 0.3e1_dp
                trhob = -t675 * t404 * phirhob / 0.2e1_dp - t402 * t680 * t672 * &
                     k_frhoa * t76 / 0.2e1_dp - t684 / 0.2e1_dp
                Arhob = -t349 * t687 * (-epsilon_c_unifrhob * t348 * t409 + &
                     0.3e1_dp * t406 * t691 * phirhob) * t411
                t1070 = Arhob * t415
                t1072 = 0.2e1_dp * t709 * trhob
                epsilon_cGGArhob = epsilon_c_unifrhob + 0.3e1_dp * t699 * t427 * &
                     phirhob + t414 * (0.2e1_dp * t703 * t704 * trhob + t349 * t415 * (&
                     t1070 + t1072) * t423 - t716 * t719 * (t1070 + t1072 + 0.2e1_dp * &
                     t720 * Arhob + 0.4e1_dp * t724 * trhob)) * t731
                tau_wrhob = tau_wrhoa
                C_chi_epsrhob = (0.174e1_dp * chi * chirhob + 0.20e1_dp * t1043 + &
                     0.1356e2_dp * t758 * chirhob) * t458 - 0.4e1_dp * t766 * (t767 * ((&
                     t735 * (0.2e1_dp * t739 - my_rhoa * t431 + my_rhoa * t1 + t736) - &
                     t746) * t440 * t442 / 0.6e1_dp - t754) + t449 * (-0.4e1_dp / &
                     0.3e1_dp * t771 * chirhob + 0.4e1_dp / 0.3e1_dp * t775 * chirhob) / &
                     0.2e1_dp)
                IF (t863) THEN
                   marhob = 0._dp
                ELSE
                   marhob = epsilon_cGGArhob
                END IF
                t1111 = t161 ** 2
                rs_s2rhob = -t148 / t1111 * t76 * t908 / 0.12e2_dp
                t1121 = k_s_s2 ** 2
                t_s2rhob = -t166 / t1121 * t6 / t164 * t12 * t976 * &
                     t14 * t76 / 0.6e1_dp - t166 * t167 * t908 / 0.2e1_dp
                t1132 = t218 ** 2
                e_c_u_1_s2rhob = -0.2e1_dp * t600 * rs_s2rhob * t223 + t209 / &
                     t1132 * (beta_1_2 / t174 * rs_s2rhob / 0.2e1_dp + beta_2_2 * &
                     rs_s2rhob + 0.3e1_dp / 0.2e1_dp * beta_3_2 * t174 * rs_s2rhob + t217&
                     * t215 * rs_s2rhob / rs_s2) / t222
                t1151 = t378 ** 2
                A_s2rhob = t827 / t1151 * e_c_u_1_s2rhob * t375 * t377
                t1156 = t349 * t_s2
                t1157 = t383 * t389
                t1161 = A_s2rhob * t381
                t1162 = A_s2 * t_s2
                t1164 = 0.2e1_dp * t1162 * t_s2rhob
                t1169 = t349 * t381
                t1170 = t388 ** 2
                t1172 = t383 / t1170
                t1176 = t381 * t_s2
                t1177 = t385 * t1176
                t1184 = 0.1e1_dp / t392
                IF (t784) THEN
                   mbrhob = e_c_u_1_s2rhob + t380 * (0.2e1_dp * t1156 * &
                     t1157 * t_s2rhob + t349 * t381 * (t1161 + t1164) * t389 - t1169 * &
                     t1172 * (t1161 + t1164 + 0.2e1_dp * A_s2 * t386 * A_s2rhob + &
                     0.4e1_dp * t1177 * t_s2rhob)) * t1184
                ELSE
                   mbrhob = epsilon_cGGArhob
                END IF
                epsilon_cRevPKZBrhob = epsilon_cGGArhob + (epsilon_cGGArhob * &
                     C_chi_eps * t460 + epsilon_cGGA * C_chi_epsrhob * t460 + 0.2e1_dp * &
                     t459 * tau_w * tau_wrhob - C_chi_epsrhob * t460 * t468 - 0.2e1_dp * &
                     t873 * t468 * tau_wrhob - t463 * (-t879 + t464 * marhob + t275 * mb &
                     - t882 + t466 * mbrhob)) * t472

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_rhob(ii)=e_rhob(ii)+&
                        scale_ex * (0.2e1_dp * ex_unif_sp2 * Fx_sp2 - my_rhob * &
                        t559 * t18 * t976 * Fx_sp2 / 0.2e1_dp + 0.2e1_dp * t144 &
                        * t974 * (((0.2e1_dp * t929 * t111 * z_sp2rhob - 0.4e1_dp * t934 * &
                        t936 * z_sp2rhob) * p_sp2 + t113 * p_sp2rhob + 0.292e3_dp / &
                        0.2025e4_dp * tildeq_b_sp2 * tildeq_b_sp2rhob - 0.73e2_dp / &
                        0.4050e4_dp * tildeq_b_sp2rhob * t121 - 0.73e2_dp / 0.8100e4_dp * &
                        t948 * (0.36e2_dp * z_sp2 * z_sp2rhob + 0.100e3_dp * p_sp2 * &
                        p_sp2rhob) + 0.200e3_dp / 0.6561e4_dp * t956 * p_sp2rhob + 0.8e1_dp &
                        / 0.45e2_dp * t959 * z_sp2rhob + 0.3e1_dp * t64 * t118 * p_sp2rhob) &
                        * t134 - 0.2e1_dp * t969 * t61 * p_sp2rhob)) / 0.2e1_dp + scale_ec *&
                        (t888 + rho * epsilon_cRevPKZBrhob * t481 + t474 * (d * &
                        epsilon_cRevPKZBrhob * t479 + 0.3e1_dp * t475 * t893 * tau_wrhob))
                END IF

                p_sp1norm_drhoa = my_norm_drhoa * t12 * t25 / 0.12e2_dp
                z_sp1norm_drhoa = my_norm_drhoa * t2 * t27 / 0.4e1_dp
                alpha_sp1norm_drhoa = 0.5e1_dp / 0.3e1_dp * p_sp1norm_drhoa * t30 &
                     - 0.5e1_dp / 0.3e1_dp * t498 * z_sp1norm_drhoa
                tildeq_b_sp1norm_drhoa = 0.9e1_dp / 0.20e2_dp * &
                     alpha_sp1norm_drhoa * t37 - 0.9e1_dp / 0.40e2_dp * t505 * (b * &
                     alpha_sp1norm_drhoa * t32 + t33 * alpha_sp1norm_drhoa) + 0.2e1_dp / &
                     0.3e1_dp * p_sp1norm_drhoa
                t_s1norm_drhoa = t155 * t157 * t2 / 0.2e1_dp
                C_chi_epsnorm_drhoa = -0.2e1_dp / 0.3e1_dp * t766 * t767 * t735 * &
                     (0.2e1_dp * my_norm_drhoa * t86 + 0.2e1_dp * t430 * my_norm_drhoa) * &
                     t440 * t442
                IF (t863) THEN
                   manorm_drhoa = t358 * (0.2e1_dp * t833 * t834 * &
                        t_s1norm_drhoa + 0.2e1_dp * t349 * t853 * A_s1 * t_s1norm_drhoa * &
                        t367 - t846 * t849 * (0.2e1_dp * t839 * t_s1norm_drhoa + 0.4e1_dp * &
                        t854 * t_s1norm_drhoa)) * t861
                ELSE
                   manorm_drhoa= 0._dp
                END IF
                epsilon_cRevPKZBnorm_drhoa = (epsilon_cGGA * C_chi_epsnorm_drhoa *&
                     t460 - C_chi_epsnorm_drhoa * t460 * t468 - t463 * t464 * &
                     manorm_drhoa) * t472
                t1303 = t474 * d

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_ndrhoa(ii)=e_ndrhoa(ii)+&
                        scale_ex * t83 * t558 * (((0.2e1_dp * t513 * &
                        t45 * z_sp1norm_drhoa - 0.4e1_dp * t518 * t520 * z_sp1norm_drhoa) * &
                        p_sp1 + t47 * p_sp1norm_drhoa + 0.292e3_dp / 0.2025e4_dp * &
                        tildeq_b_sp1 * tildeq_b_sp1norm_drhoa - 0.73e2_dp / 0.4050e4_dp * &
                        tildeq_b_sp1norm_drhoa * t55 - 0.73e2_dp / 0.8100e4_dp * t532 * (&
                        0.36e2_dp * z_sp1 * z_sp1norm_drhoa + 0.100e3_dp * p_sp1 * &
                        p_sp1norm_drhoa) + 0.200e3_dp / 0.6561e4_dp * t540 * p_sp1norm_drhoa&
                        + 0.8e1_dp / 0.45e2_dp * t543 * z_sp1norm_drhoa + 0.3e1_dp * t64 * &
                        t52 * p_sp1norm_drhoa) * t71 - 0.2e1_dp * t553 * t61 * &
                        p_sp1norm_drhoa) + scale_ec * (rho * epsilon_cRevPKZBnorm_drhoa * &
                        t481 + t1303 * epsilon_cRevPKZBnorm_drhoa * t476 * t478)
                END IF

                p_sp2norm_drhob = my_norm_drhob * t12 * t91 / 0.12e2_dp
                z_sp2norm_drhob = my_norm_drhob * t6 * t93 / 0.4e1_dp
                alpha_sp2norm_drhob = 0.5e1_dp / 0.3e1_dp * p_sp2norm_drhob * t96 &
                     - 0.5e1_dp / 0.3e1_dp * t914 * z_sp2norm_drhob
                tildeq_b_sp2norm_drhob = 0.9e1_dp / 0.20e2_dp * &
                     alpha_sp2norm_drhob * t103 - 0.9e1_dp / 0.40e2_dp * t921 * (b * &
                     alpha_sp2norm_drhob * t98 + t99 * alpha_sp2norm_drhob) + 0.2e1_dp / &
                     0.3e1_dp * p_sp2norm_drhob
                t_s2norm_drhob = t165 * t167 * t6 / 0.2e1_dp
                C_chi_epsnorm_drhob = -0.2e1_dp / 0.3e1_dp * t766 * t767 * t735 * &
                     (0.2e1_dp * my_norm_drhob * t20 + 0.2e1_dp * t430 * my_norm_drhob) * &
                     t440 * t442
                IF (t784) THEN
                   mbnorm_drhob = t380 * (0.2e1_dp * t1156 * t1157 * &
                        t_s2norm_drhob + 0.2e1_dp * t349 * t1176 * A_s2 * t_s2norm_drhob * &
                        t389 - t1169 * t1172 * (0.2e1_dp * t1162 * t_s2norm_drhob + 0.4e1_dp&
                        * t1177 * t_s2norm_drhob)) * t1184
                ELSE
                   mbnorm_drhob= 0._dp
                END IF
                epsilon_cRevPKZBnorm_drhob = (epsilon_cGGA * C_chi_epsnorm_drhob *&
                     t460 - C_chi_epsnorm_drhob * t460 * t468 - t463 * t466 * &
                     mbnorm_drhob) * t472

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_ndrhob(ii)=e_ndrhob(ii)+&
                        scale_ex * t144 * t974 * (((0.2e1_dp * t929 * &
                        t111 * z_sp2norm_drhob - 0.4e1_dp * t934 * t936 * z_sp2norm_drhob) *&
                        p_sp2 + t113 * p_sp2norm_drhob + 0.292e3_dp / 0.2025e4_dp * &
                        tildeq_b_sp2 * tildeq_b_sp2norm_drhob - 0.73e2_dp / 0.4050e4_dp * &
                        tildeq_b_sp2norm_drhob * t121 - 0.73e2_dp / 0.8100e4_dp * t948 * (&
                        0.36e2_dp * z_sp2 * z_sp2norm_drhob + 0.100e3_dp * p_sp2 * &
                        p_sp2norm_drhob) + 0.200e3_dp / 0.6561e4_dp * t956 * p_sp2norm_drhob&
                        + 0.8e1_dp / 0.45e2_dp * t959 * z_sp2norm_drhob + 0.3e1_dp * t64 * &
                        t118 * p_sp2norm_drhob) * t134 - 0.2e1_dp * t969 * t61 * &
                        p_sp2norm_drhob) + scale_ec * (rho * epsilon_cRevPKZBnorm_drhob * &
                        t481 + t1303 * epsilon_cRevPKZBnorm_drhob * t476 * t478)
                END IF

                tnorm_drho = t401 * t403 * t275 / 0.2e1_dp
                Hnorm_drho = t414 * (0.2e1_dp * t703 * t704 * tnorm_drho + &
                     0.2e1_dp * t349 * t723 * A * tnorm_drho * t423 - t716 * t719 * (&
                     0.2e1_dp * t709 * tnorm_drho + 0.4e1_dp * t724 * tnorm_drho)) * t731
                tau_wnorm_drho = my_norm_drho * t275 / 0.4e1_dp
                C_chi_epsnorm_drho = 0.4e1_dp / 0.3e1_dp * t766 * t767 * t735 * &
                     t430 * my_norm_drho * t440 * t442
                IF (t863) THEN
                   manorm_drho = 0._dp
                ELSE
                   manorm_drho = Hnorm_drho
                END IF
                IF (t784) THEN
                   mbnorm_drho =  0._dp
                ELSE
                   mbnorm_drho = Hnorm_drho
                END IF
                epsilon_cRevPKZBnorm_drho = Hnorm_drho + (Hnorm_drho * C_chi_eps *&
                     t460 + epsilon_cGGA * C_chi_epsnorm_drho * t460 + 0.2e1_dp * t459 *&
                     tau_w * tau_wnorm_drho - C_chi_epsnorm_drho * t460 * t468 - &
                     0.2e1_dp * t873 * t468 * tau_wnorm_drho - t463 * (t464 * manorm_drho&
                     + t466 * mbnorm_drho)) * t472

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_ndrho(ii)=e_ndrho(ii)+&
                        scale_ec * (rho * epsilon_cRevPKZBnorm_drho * &
                        t481 + t474 * (d * epsilon_cRevPKZBnorm_drho * t479 + 0.3e1_dp * &
                        t475 * t893 * tau_wnorm_drho))
                END IF

                t1457 = my_tau_a ** 2
                z_sp1tau_a = -tau_w_sp1 / t1457 / 0.2e1_dp
                alpha_sp1tau_a = -0.5e1_dp / 0.3e1_dp * t498 * z_sp1tau_a
                tildeq_b_sp1tau_a = 0.9e1_dp / 0.20e2_dp * alpha_sp1tau_a * t37 - &
                     0.9e1_dp / 0.40e2_dp * t505 * (b * alpha_sp1tau_a * t32 + t33 * &
                     alpha_sp1tau_a)
                epsilon_cRevPKZBtau_a = -0.2e1_dp * t470 * t478
                t1496 = t471 ** 2
                t1500 = 0.3e1_dp * t475 * t476 / t1496

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_tau_a(ii)=e_tau_a(ii)+&
                        scale_ex * t83 * t558 * ((0.2e1_dp * t513 * t45 * &
                        z_sp1tau_a - 0.4e1_dp * t518 * t520 * z_sp1tau_a) * p_sp1 + &
                        0.292e3_dp / 0.2025e4_dp * tildeq_b_sp1 * tildeq_b_sp1tau_a - &
                        0.73e2_dp / 0.4050e4_dp * tildeq_b_sp1tau_a * t55 - 0.73e2_dp / &
                        0.225e3_dp * t532 * z_sp1 * z_sp1tau_a + 0.8e1_dp / 0.45e2_dp * t543&
                        * z_sp1tau_a) * t71 + scale_ec * (rho * epsilon_cRevPKZBtau_a * &
                        t481 + t474 * (d * epsilon_cRevPKZBtau_a * t479 - t1500))
                END IF

                t1506 = my_tau_b ** 2
                z_sp2tau_b = -tau_w_sp2 / t1506 / 0.2e1_dp
                alpha_sp2tau_b = -0.5e1_dp / 0.3e1_dp * t914 * z_sp2tau_b
                tildeq_b_sp2tau_b = 0.9e1_dp / 0.20e2_dp * alpha_sp2tau_b * t103 -&
                     0.9e1_dp / 0.40e2_dp * t921 * (b * alpha_sp2tau_b * t98 + t99 * &
                     alpha_sp2tau_b)
                epsilon_cRevPKZBtau_b = epsilon_cRevPKZBtau_a
                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_tau_b(ii)=e_tau_b(ii)+&
                        scale_ex * t144 * t974 * ((0.2e1_dp * t929 * t111 * &
                        z_sp2tau_b - 0.4e1_dp * t934 * t936 * z_sp2tau_b) * p_sp2 + &
                        0.292e3_dp / 0.2025e4_dp * tildeq_b_sp2 * tildeq_b_sp2tau_b - &
                        0.73e2_dp / 0.4050e4_dp * tildeq_b_sp2tau_b * t121 - 0.73e2_dp / &
                        0.225e3_dp * t948 * z_sp2 * z_sp2tau_b + 0.8e1_dp / 0.45e2_dp * t959&
                        * z_sp2tau_b) * t134 + scale_ec * (rho * epsilon_cRevPKZBtau_b * &
                        t481 + t474 * (d * epsilon_cRevPKZBtau_b * t479 - t1500))
                END IF
             END IF
          END IF
       END DO

       !$omp end do

    END IF
  END SUBROUTINE tpss_lsd_calc

END MODULE xc_tpss
