!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_tpss [1.0] *
!!
!!   NAME
!!     xc_tpss
!!
!!   FUNCTION
!!     Calculates the tpss functional.
!!
!!   NOTES
!!     The derivation of the formulaes is lengthly, and not fully trivial,
!!     so I have put it in doc/tpss.mw
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     05.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE xc_tpss
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_tpss'
  
  PUBLIC :: tpss_lda_info, tpss_lda_eval, tpss_lsd_info, tpss_lsd_eval

!***
!****************************************************************************
CONTAINS

!!****f* xc_tpss/tpss_lda_info [1.0] *
!!
!!   NAME
!!     tpss_lda_info
!!
!!   FUNCTION
!!     return various information on the functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - max_deriv: the highest derivative available
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lda_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "TPSS meta-GGA functional (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%tau=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE tpss_lda_info
  !***************************************************************************

!!****f* xc_tpss/tpss_lsd_info [1.0] *
!!
!!   NAME
!!     tpss_lsd_info
!!
!!   FUNCTION
!!     return various information on the functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - max_deriv: the highest derivative available
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lsd_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "TPSS meta-GGA functional (LSD)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%tau_spin=.TRUE.
       needs%norm_drho=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE tpss_lsd_info
  !***************************************************************************

!!****f* xc_tpss/tpss_lda_eval [1.0] *
!!
!!   NAME
!!     tpss_lda_eval
!!
!!   FUNCTION
!!     evaluates the tpss functional in the spin unpolarized (lda) case
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the density where you want to evaluate the functional
!!     - deriv_set: place where to store the functional derivatives (they are
!!       added to the derivatives)
!!     - grad_deriv: degree of the derivative that should be evalated,
!!       if positive all the derivatives up to the given degree are evaluated,
!!       if negative only the given degree is calculated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lda_eval(rho_set,deriv_set,grad_deriv,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, epsilon_tau
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrho, e_rho, &
                                                e_tau, norm_drho, rho, tau
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo)

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

       CALL tpss_lda_calc(rho=rho, norm_drho=norm_drho,&
            tau=tau,e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_tau=e_tau,&
            grad_deriv=grad_deriv, npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_norm_drho=epsilon_norm_drho,epsilon_tau=epsilon_tau,&
            error=error)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE tpss_lda_eval
  !***************************************************************************

!!****f* xc_tpss/tpss_lsd_eval [1.0] *
!!
!!   NAME
!!     tpss_lsd_eval
!!
!!   FUNCTION
!!     evaluates the tpss functional in the spin polarized (lsd) case
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the density where you want to evaluate the functional
!!     - deriv_set: place where to store the functional derivatives (they are
!!       added to the derivatives)
!!     - grad_deriv: degree of the derivative that should be evalated,
!!       if positive all the derivatives up to the given degree are evaluated,
!!       if negative only the given degree is calculated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lsd_eval(rho_set,deriv_set,grad_deriv,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, epsilon_tau
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrhoa, e_ndrhob, e_rhoa, e_rhob, e_tau_a, e_tau_b, norm_drho, &
      norm_drhoa, norm_drhob, rhoa, rhob, tau_a, tau_b
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo)

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

       CALL tpss_lsd_calc(rhoa=rhoa,rhob=rhob, norm_drhoa=norm_drhoa, &
            norm_drhob=norm_drhob, norm_drho=norm_drho,&
            tau_a=tau_a,tau_b=tau_b,e_0=e_0,e_rhoa=e_rhoa,e_rhob=e_rhob,&
            e_ndrhoa=e_ndrhoa,e_ndrhob=e_ndrhob,e_ndrho=e_ndrho,&
            e_tau_a=e_tau_a,e_tau_b=e_tau_b,&
            grad_deriv=grad_deriv, npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_norm_drho=epsilon_norm_drho,epsilon_tau=epsilon_tau,&
            error=error)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE tpss_lsd_eval
  !***************************************************************************

!!****f* xc_tpss/tpss_lda_calc *
!!
!!   NAME
!!     tpss_lda_calc
!!
!!   FUNCTION
!!     low level calculation routine for the unpolarized (lda) tpss
!!
!!   NOTES
!!     maple is nice, but if you want the uman readable version of the code
!!     look in doc/tpss.mw
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lda_calc(rho,norm_drho,tau,e_0,e_rho,e_ndrho,e_tau,&
       npoints,grad_deriv,epsilon_rho,epsilon_norm_drho,epsilon_tau,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho, norm_drho, tau
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rho, e_ndrho, e_tau
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, epsilon_tau
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: abs_grad_deriv, ii, non_coer
    LOGICAL                                  :: failure, t535, t603
    REAL(kind=dp) :: A, A_s1, A_s1rho, A_s2, A_s2rho, alpha, alpha_c, &
      alpha_c_s1, alpha_c_s2, alphanorm_drho, alpharho, alphatau, Arho, b, &
      beta, beta_s1, beta_s2, c, C_chi, C_chi_eps, chi_s1, chi_s2, d, &
      e_c_u_0, e_c_u_0_s1, e_c_u_0_s2, e_c_u_0rho, e_c_u_1, e_c_u_1_s1, &
      e_c_u_1_s1rho, e_c_u_1_s2, e_c_u_1_s2rho, e_var, ec, ecnorm_drho, &
      ecrho, ectau, epsilon_c_unif, epsilon_c_unif_s1, epsilon_c_unif_s1rho, &
      epsilon_c_unif_s2, epsilon_c_unif_s2rho, epsilon_c_unifrho, &
      epsilon_cGGA, epsilon_cGGA_0_1, epsilon_cGGA_0_1norm_drho, &
      epsilon_cGGA_0_1rho, epsilon_cGGA_1_0, epsilon_cGGA_1_0norm_drho, &
      epsilon_cGGA_1_0rho, epsilon_cGGAnorm_drho
    REAL(kind=dp) :: epsilon_cGGArho, epsilon_cRevPKZB, &
      epsilon_cRevPKZBnorm_drho, epsilon_cRevPKZBrho, epsilon_cRevPKZBtau, &
      ex_lda, ex_ldanorm_drho, ex_ldarho, ex_ldatau, ex_unif, ex_unifrho, &
      f_ii_0, f_ii_0_s1, f_ii_0_s2, f_s1, f_s2, Fx, Fxnorm_drho, Fxrho, &
      Fxtau, gamma_var, gamma_var_s1, gamma_var_s2, H, H_s1, H_s1norm_drho, &
      H_s1rho, H_s2, H_s2norm_drho, H_s2rho, Hnorm_drho, Hrho, k_f, k_f_s1, &
      k_f_s1rho, k_f_s2, k_f_s2rho, k_frho, k_s, k_s_s1, k_s_s1rho, k_s_s2, &
      k_s_s2rho, k_srho, kappa, m, ma, manorm_drho, marho, mb, mbnorm_drho, &
      mbrho, mu, my_ndrho, my_rho, my_tau, p, phi, phi_s1, phi_s2, &
      pnorm_drho, prho, rs
    REAL(kind=dp) :: rs_s1, rs_s1rho, rs_s2, rs_s2rho, rsrho, t, t1, t102, &
      t104, t110, t114, t119, t122, t123, t136, t139, t14, t140, t141, t142, &
      t145, t147, t148, t149, t151, t152, t155, t156, t157, t158, t159, t160, &
      t161, t163, t164, t165, t168, t169, t17, t172, t173, t176, t178, t180, &
      t183, t184, t19, t197, t2, t210, t213, t214, t215, t217, t219, t220, &
      t223, t224, t225, t226, t227, t228, t230, t231, t232, t235, t236, t237, &
      t239, t240, t241, t242, t244, t247, t25, t250, t251, t252, t254, t257, &
      t258, t259, t260, t261, t262, t263, t264, t267, t269, t271, t273, t276, &
      t277, t278, t282, t283, t284, t285, t288, t29, t290
    REAL(kind=dp) :: t293, t296, t3, t303, t305, t306, t307, t309, t315, &
      t316, t319, t320, t325, t334, t34, t342, t343, t348, t355, t357, t362, &
      t366, t367, t369, t37, t374, t376, t38, t380, t381, t390, t4, t401, &
      t408, t417, t423, t434, t435, t437, t450, t453, t454, t468, t477, t489, &
      t497, t499, t5, t505, t509, t51, t510, t512, t517, t519, t523, t524, &
      t532, t538, t547, t557, t56, t565, t567, t573, t577, t578, t58, t580, &
      t585, t587, t591, t592, t6, t60, t600, t605, t609, t61, t623, t63, t66, &
      t67, t675, t677, t68, t69, t7, t72, t74, t75, t76, t763, t78, t79, &
      t790, t82, t83, t84, t85, t86, t87, t88, t9, t90, t91, t92
    REAL(kind=dp) :: t95, t96, t99, t_s1, t_s1norm_drho, t_s1rho, t_s2, &
      t_s2norm_drho, t_s2rho, tau_w, tau_wnorm_drho, tau_wrho, tildeq_b, &
      tildeq_bnorm_drho, tildeq_brho, tildeq_btau, tnorm_drho, trho, x, &
      xnorm_drho, xrho, xtau, z, znorm_drho, zrho, ztau
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_error_get_logger(error)

    IF (.NOT. failure) THEN
       abs_grad_deriv=ABS(grad_deriv)
       non_coer=0

       !$omp parallel do default(private) &
       !$omp             shared(rho,tau,norm_drho,e_0,e_rho,e_ndrho,e_tau,&
       !$omp                    epsilon_rho,epsilon_tau)
       DO ii=1,npoints
          my_tau=tau(ii)
          my_rho=rho(ii)
          IF (my_rho>epsilon_rho.and.my_tau>epsilon_tau) THEN
             my_ndrho=norm_drho(ii)
             
             chi_s1 = 1.0_dp
             t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
             t2 = 4 ** (0.1e1_dp / 0.3e1_dp)
             t3 = t2 ** 2
             t4 = t1 * t3
             t5 = 2 ** (0.1e1_dp / 0.3e1_dp)
             t6 = 0.1e1_dp / 0.3141592654e1_dp
             t7 = 0.1e1_dp / my_rho
             t237 = my_ndrho ** 2
             tau_w = t237 * t7 / 0.8e1_dp

             IF (my_tau<tau_w) THEN
                ! enforce z=norm_rho**2/(8._dp*rho*tau) <1
                m=0.5_dp*t237+4.0_dp*my_rho*my_tau
                my_tau=m/8._dp/my_rho
                my_ndrho=SQRT(m)
                t237=m
                non_coer=non_coer+1
             END IF

             t9 = (t6 * t7) ** (0.1e1_dp / 0.3e1_dp)
             rs_s1 = t4 * t5 * t9 / 0.4e1_dp
             t14 = SQRT(rs_s1)
             t17 = t14 * rs_s1
             t19 = rs_s1 ** 0.20e1_dp
             t25 = LOG(0.1e1_dp + 0.1608182432e2_dp / (0.75957e1_dp * t14 + &
                  0.35876e1_dp * rs_s1 + 0.16382e1_dp * t17 + 0.49294e0_dp * t19))
             e_c_u_0_s1 = -0.62182e-1_dp * (0.1e1_dp + 0.21370e0_dp * rs_s1) * t25
             t29 = 0.1e1_dp + 0.20548e0_dp * rs_s1
             t34 = 0.141189e2_dp * t14 + 0.61977e1_dp * rs_s1 + 0.33662e1_dp * &
                  t17 + 0.62517e0_dp * t19
             t37 = 0.1e1_dp + 0.3216468318e2_dp / t34
             t38 = LOG(t37)
             e_c_u_1_s1 = -0.31090e-1_dp * t29 * t38
             t51 = LOG(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t14 + &
                  0.36231e1_dp * rs_s1 + 0.88026e0_dp * t17 + 0.49671e0_dp * t19))
             alpha_c_s1 = 0.33774e0_dp * (0.1e1_dp + 0.11125e0_dp * rs_s1) * t51
             f_ii_0_s1 = 0.8e1_dp / 0.9e1_dp / (2 * t5 - 2)
             f_s1 = 1
             epsilon_c_unif_s1 = e_c_u_1_s1
             beta_s1 = 0.66725e-1_dp
             t56 = LOG(0.2e1_dp)
             t58 = 0.3141592654e1_dp ** 2
             gamma_var_s1 = (0.1e1_dp - t56) / t58
             t60 = t5 ** 2
             phi_s1 = t60 / 0.2e1_dp
             t61 = t1 * t60
             t63 = (t58 * my_rho) ** (0.1e1_dp / 0.3e1_dp)
             k_f_s1 = t61 * t63 / 0.2e1_dp
             t66 = SQRT(k_f_s1 * t6)
             k_s_s1 = 0.2e1_dp * t66
             t67 = 0.1e1_dp / phi_s1
             t68 = my_ndrho * t67
             t69 = 0.1e1_dp / k_s_s1
             t_s1 = t68 * t69 * t7 / 0.2e1_dp
             t72 = 0.1e1_dp / gamma_var_s1
             t74 = phi_s1 ** 2
             t75 = t74 * phi_s1
             t76 = 0.1e1_dp / t75
             t78 = EXP(-e_c_u_1_s1 * t72 * t76)
             t79 = t78 - 0.1e1_dp
             A_s1 = 0.66725e-1_dp * t72 / t79
             t82 = gamma_var_s1 * t75
             t83 = t_s1 ** 2
             t84 = t72 * t83
             t85 = A_s1 * t83
             t86 = 0.1e1_dp + t85
             t87 = A_s1 ** 2
             t88 = t83 ** 2
             t90 = 0.1e1_dp + t85 + t87 * t88
             t91 = 0.1e1_dp / t90
             t92 = t86 * t91
             t95 = 0.1e1_dp + 0.66725e-1_dp * t84 * t92
             t96 = LOG(t95)
             H_s1 = t82 * t96
             epsilon_cGGA_1_0 = e_c_u_1_s1 + H_s1
             chi_s2 = -1._dp
             rs_s2 = rs_s1
             t99 = SQRT(rs_s2)
             t102 = t99 * rs_s2
             t104 = rs_s2 ** 0.20e1_dp
             t110 = LOG(0.1e1_dp + 0.1608182432e2_dp / (0.75957e1_dp * t99 + &
                  0.35876e1_dp * rs_s2 + 0.16382e1_dp * t102 + 0.49294e0_dp * t104))
             e_c_u_0_s2 = -0.62182e-1_dp * (0.1e1_dp + 0.21370e0_dp * rs_s2) * t110
             t114 = 0.1e1_dp + 0.20548e0_dp * rs_s2
             t119 = 0.141189e2_dp * t99 + 0.61977e1_dp * rs_s2 + 0.33662e1_dp * &
                  t102 + 0.62517e0_dp * t104
             t122 = 0.1e1_dp + 0.3216468318e2_dp / t119
             t123 = LOG(t122)
             e_c_u_1_s2 = -0.31090e-1_dp * t114 * t123
             t136 = LOG(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t99 + &
                  0.36231e1_dp * rs_s2 + 0.88026e0_dp * t102 + 0.49671e0_dp * t104))
             alpha_c_s2 = 0.33774e0_dp * (0.1e1_dp + 0.11125e0_dp * rs_s2) * t136
             f_ii_0_s2 = f_ii_0_s1
             f_s2 = 1
             epsilon_c_unif_s2 = e_c_u_1_s2
             beta_s2 = 0.66725e-1_dp
             gamma_var_s2 = gamma_var_s1
             phi_s2 = phi_s1
             k_f_s2 = k_f_s1
             t139 = SQRT(k_f_s2 * t6)
             k_s_s2 = 0.2e1_dp * t139
             t140 = 0.1e1_dp / phi_s2
             t141 = my_ndrho * t140
             t142 = 0.1e1_dp / k_s_s2
             t_s2 = t141 * t142 * t7 / 0.2e1_dp
             t145 = 0.1e1_dp / gamma_var_s2
             t147 = phi_s2 ** 2
             t148 = t147 * phi_s2
             t149 = 0.1e1_dp / t148
             t151 = EXP(-e_c_u_1_s2 * t145 * t149)
             t152 = t151 - 0.1e1_dp
             A_s2 = 0.66725e-1_dp * t145 / t152
             t155 = gamma_var_s2 * t148
             t156 = t_s2 ** 2
             t157 = t145 * t156
             t158 = A_s2 * t156
             t159 = 0.1e1_dp + t158
             t160 = A_s2 ** 2
             t161 = t156 ** 2
             t163 = 0.1e1_dp + t158 + t160 * t161
             t164 = 0.1e1_dp / t163
             t165 = t159 * t164
             t168 = 0.1e1_dp + 0.66725e-1_dp * t157 * t165
             t169 = LOG(t168)
             H_s2 = t155 * t169
             epsilon_cGGA_0_1 = e_c_u_1_s2 + H_s2
             rs = t4 * t9 / 0.4e1_dp
             t172 = 0.1e1_dp + 0.21370e0_dp * rs
             t173 = SQRT(rs)
             t176 = t173 * rs
             t178 = rs ** 0.20e1_dp
             t180 = 0.75957e1_dp * t173 + 0.35876e1_dp * rs + 0.16382e1_dp * t176&
                  + 0.49294e0_dp * t178
             t183 = 0.1e1_dp + 0.1608182432e2_dp / t180
             t184 = LOG(t183)
             e_c_u_0 = -0.62182e-1_dp * t172 * t184
             t197 = LOG(0.1e1_dp + 0.3216468318e2_dp / (0.141189e2_dp * t173 + &
                  0.61977e1_dp * rs + 0.33662e1_dp * t176 + 0.62517e0_dp * t178))
             e_c_u_1 = -0.31090e-1_dp * (0.1e1_dp + 0.20548e0_dp * rs) * t197
             t210 = LOG(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t173 + &
                  0.36231e1_dp * rs + 0.88026e0_dp * t176 + 0.49671e0_dp * t178))
             alpha_c = 0.33774e0_dp * (0.1e1_dp + 0.11125e0_dp * rs) * t210
             f_ii_0 = f_ii_0_s2
             epsilon_c_unif = e_c_u_0
             beta = 0.66725e-1_dp
             gamma_var = gamma_var_s2
             phi = 1
             k_f = t1 * t63
             t213 = SQRT(k_f * t6)
             k_s = 0.2e1_dp * t213
             t214 = 0.1e1_dp / k_s
             t215 = my_ndrho * t214
             t = t215 * t7 / 0.2e1_dp
             t217 = 0.1e1_dp / gamma_var
             t219 = EXP(-e_c_u_0 * t217)
             t220 = -0.1e1_dp + t219
             A = 0.66725e-1_dp * t217 / t220
             t223 = t ** 2
             t224 = t217 * t223
             t225 = A * t223
             t226 = 0.1e1_dp + t225
             t227 = A ** 2
             t228 = t223 ** 2
             t230 = 0.1e1_dp + t225 + t227 * t228
             t231 = 0.1e1_dp / t230
             t232 = t226 * t231
             t235 = 0.1e1_dp + 0.66725e-1_dp * t224 * t232
             t236 = LOG(t235)
             H = gamma_var * t236
             epsilon_cGGA = e_c_u_0 + H
             d = 0.28e1_dp
             C_chi = 0.53e0_dp
             C_chi_eps = 0.53e0_dp
             ma = MAX(epsilon_cGGA_1_0, epsilon_cGGA)
             mb = MAX(epsilon_cGGA_0_1, epsilon_cGGA)
             t239 = tau_w ** 2
             t240 = my_tau ** 2
             t241 = 0.1e1_dp / t240
             t242 = t239 * t241
             t244 = 0.1e1_dp + 0.53e0_dp * t242
             t247 = ma / 0.2e1_dp + mb / 0.2e1_dp
             epsilon_cRevPKZB = epsilon_cGGA * t244 - 0.153e1_dp * t242 * t247
             t250 = my_rho * epsilon_cRevPKZB
             t251 = t239 * tau_w
             t252 = epsilon_cRevPKZB * t251
             t254 = 0.1e1_dp / t240 / my_tau
             t257 = 0.1e1_dp + 0.28e1_dp * t252 * t254
             ec = t250 * t257
             kappa = 0.804e0_dp
             mu = 0.21951e0_dp
             b = 0.4e0_dp
             c = 0.159096e1_dp
             e_var = 0.1537e1_dp
             t258 = t237 * t1
             t259 = t58 ** (0.1e1_dp / 0.3e1_dp)
             t260 = t259 ** 2
             t261 = 0.1e1_dp / t260
             t262 = my_rho ** 2
             t263 = my_rho ** (0.1e1_dp / 0.3e1_dp)
             t264 = t263 ** 2
             t267 = t261 / t264 / t262
             p = t258 * t267 / 0.12e2_dp
             t269 = 0.1e1_dp / my_tau
             z = tau_w * t269
             t271 = 0.1e1_dp / z - 0.1e1_dp
             alpha = 0.5e1_dp / 0.3e1_dp * p * t271
             t273 = alpha - 0.1e1_dp
             t276 = 0.1e1_dp + 0.4e0_dp * alpha * t273
             t277 = SQRT(t276)
             t278 = 0.1e1_dp / t277
             tildeq_b = 0.9e1_dp / 0.20e2_dp * t273 * t278 + 0.2e1_dp / 0.3e1_dp * p
             t282 = z ** 2
             t283 = 0.1e1_dp + t282
             t284 = t283 ** 2
             t285 = 0.1e1_dp / t284
             t288 = 0.10e2_dp / 0.81e2_dp + 0.159096e1_dp * t282 * t285
             t290 = tildeq_b ** 2
             t293 = p ** 2
             t296 = SQRT(0.18e2_dp * t282 + 0.50e2_dp * t293)
             t303 = t288 * p + 0.146e3_dp / 0.2025e4_dp * t290 - 0.73e2_dp / &
                  0.4050e4_dp * tildeq_b * t296 + 0.1895718785e-1_dp * t293 + &
                  0.1102007148e0_dp * t282 + 0.33738687e0_dp * t293 * p
             t305 = 0.1e1_dp + 0.1239758041e1_dp * p
             t306 = t305 ** 2
             t307 = 0.1e1_dp / t306
             x = t303 * t307
             t309 = 0.1e1_dp + 0.1243781095e1_dp * x
             Fx = 0.1804e1_dp - 0.804e0_dp / t309
             ex_unif = -0.3e1_dp / 0.4e1_dp * t6 * t1 * t63
             t315 = my_rho * ex_unif
             ex_lda = t315 * Fx

             IF (grad_deriv>=0) THEN
                e_0(ii) = e_0(ii)+&
                     ec + ex_lda
             END IF

             IF (abs_grad_deriv>0) THEN
                t316 = t9 ** 2
                t319 = 0.1e1_dp / t262
                t320 = 0.1e1_dp / t316 * t6 * t319
                rsrho = -t4 * t320 / 0.12e2_dp
                t325 = t180 ** 2
                t334 = rs ** 0.10e1_dp
                e_c_u_0rho = -0.1328829340e-1_dp * rsrho * t184 + 0.9999999999e0_dp &
                     * t172 / t325 * (0.3797850000e1_dp / t173 * rsrho + 0.35876e1_dp * &
                     rsrho + 0.2457300000e1_dp * t173 * rsrho + 0.985880e0_dp * t334 * &
                     rsrho) / t183
                epsilon_c_unifrho = e_c_u_0rho
                t342 = t63 ** 2
                t343 = 0.1e1_dp / t342
                k_frho = t1 * t343 * t58 / 0.3e1_dp
                k_srho = 0.1e1_dp / t213 * k_frho * t6
                t348 = k_s ** 2
                trho = -my_ndrho / t348 * t7 * k_srho / 0.2e1_dp - t215 * t319 / 0.2e1_dp
                t355 = gamma_var ** 2
                t357 = t220 ** 2
                Arho = 0.66725e-1_dp / t355 / t357 * e_c_u_0rho * t219
                t362 = t217 * t
                t366 = Arho * t223
                t367 = A * t
                t369 = 0.2e1_dp * t367 * trho
                t374 = t230 ** 2
                t376 = t226 / t374
                t380 = t223 * t
                t381 = t227 * t380
                t390 = 0.1e1_dp / t235
                Hrho = gamma_var * (0.133450e0_dp * t362 * t232 * trho + &
                     0.66725e-1_dp * t224 * (t366 + t369) * t231 - 0.66725e-1_dp * t224 *&
                     t376 * (t366 + t369 + 0.2e1_dp * A * t228 * Arho + 0.4e1_dp * t381 &
                     * trho)) * t390
                epsilon_cGGArho = e_c_u_0rho + Hrho
                tau_wrho = -t237 * t319 / 0.8e1_dp
                prho = -0.2e1_dp / 0.9e1_dp * t258 * t261 / t264 / t262 / my_rho
                zrho = tau_wrho * t269
                t401 = p / t282
                alpharho = 0.5e1_dp / 0.3e1_dp * prho * t271 - 0.5e1_dp / 0.3e1_dp *&
                     t401 * zrho
                t408 = t273 / t277 / t276
                tildeq_brho = 0.9e1_dp / 0.20e2_dp * alpharho * t278 - 0.9e1_dp / &
                     0.40e2_dp * t408 * (0.4e0_dp * alpharho * t273 + 0.4e0_dp * alpha * &
                     alpharho) + 0.2e1_dp / 0.3e1_dp * prho
                t417 = z * t285
                t423 = t282 * z / t284 / t283
                t434 = tildeq_b / t296
                t435 = z * zrho
                t437 = p * prho
                t450 = t303 / t306 / t305
                xrho = ((0.318192e1_dp * t417 * zrho - 0.636384e1_dp * t423 * zrho) &
                     * p + t288 * prho + 0.292e3_dp / 0.2025e4_dp * tildeq_b * &
                     tildeq_brho - 0.73e2_dp / 0.4050e4_dp * tildeq_brho * t296 - &
                     0.73e2_dp / 0.8100e4_dp * t434 * (0.36e2_dp * t435 + 0.100e3_dp * &
                     t437) + 0.3791437570e-1_dp * t437 + 0.2204014296e0_dp * t435 + &
                     0.101216061e1_dp * t293 * prho) * t307 - 0.2479516082e1_dp * t450 * &
                     prho
                t453 = t309 ** 2
                t454 = 0.1e1_dp / t453
                Fxrho = 0.1000000000e1_dp * t454 * xrho
                ex_unifrho = -0.3141592654e1_dp * t1 * t343 / 0.4e1_dp
                ex_ldarho = ex_unif * Fx + my_rho * ex_unifrho * Fx + t315 * Fxrho
                rs_s1rho = -t4 * t5 * t320 / 0.12e2_dp
                t468 = t34 ** 2
                t477 = rs_s1 ** 0.10e1_dp
                e_c_u_1_s1rho = -0.638837320e-2_dp * rs_s1rho * t38 + &
                     0.1000000000e1_dp * t29 / t468 * (0.7059450000e1_dp / t14 * rs_s1rho&
                     + 0.61977e1_dp * rs_s1rho + 0.5049300000e1_dp * t14 * rs_s1rho + &
                     0.1250340e1_dp * t477 * rs_s1rho) / t37
                epsilon_c_unif_s1rho = e_c_u_1_s1rho
                k_f_s1rho = t61 * t343 * t58 / 0.6e1_dp
                k_s_s1rho = 0.1e1_dp / t66 * k_f_s1rho * t6
                t489 = k_s_s1 ** 2
                t_s1rho = -t68 / t489 * t7 * k_s_s1rho / 0.2e1_dp - t68 * t69 * t319&
                     / 0.2e1_dp
                t497 = gamma_var_s1 ** 2
                t499 = t79 ** 2
                A_s1rho = 0.66725e-1_dp / t497 / t499 * e_c_u_1_s1rho * t76 * t78
                t505 = t72 * t_s1
                t509 = A_s1rho * t83
                t510 = A_s1 * t_s1
                t512 = 0.2e1_dp * t510 * t_s1rho
                t517 = t90 ** 2
                t519 = t86 / t517
                t523 = t83 * t_s1
                t524 = t87 * t523
                t532 = 0.1e1_dp / t95
                H_s1rho = t82 * (0.133450e0_dp * t505 * t92 * t_s1rho + &
                     0.66725e-1_dp * t84 * (t509 + t512) * t91 - 0.66725e-1_dp * t84 * &
                     t519 * (t509 + t512 + 0.2e1_dp * A_s1 * t88 * A_s1rho + 0.4e1_dp * &
                     t524 * t_s1rho)) * t532
                epsilon_cGGA_1_0rho = e_c_u_1_s1rho + H_s1rho
                t535 = epsilon_cGGA - epsilon_cGGA_1_0 .LT. 0.0e0_dp
                IF(t535)THEN
                   marho=epsilon_cGGA_1_0rho
                ELSE
                   marho=epsilon_cGGArho
                END IF
                rs_s2rho = rs_s1rho
                t538 = t119 ** 2
                t547 = rs_s2 ** 0.10e1_dp
                e_c_u_1_s2rho = -0.638837320e-2_dp * rs_s2rho * t123 + &
                     0.1000000000e1_dp * t114 / t538 * (0.7059450000e1_dp / t99 * &
                     rs_s2rho + 0.61977e1_dp * rs_s2rho + 0.5049300000e1_dp * t99 * &
                     rs_s2rho + 0.1250340e1_dp * t547 * rs_s2rho) / t122
                epsilon_c_unif_s2rho = e_c_u_1_s2rho
                k_f_s2rho = k_f_s1rho
                k_s_s2rho = 0.1e1_dp / t139 * k_f_s2rho * t6
                t557 = k_s_s2 ** 2
                t_s2rho = -t141 / t557 * t7 * k_s_s2rho / 0.2e1_dp - t141 * t142 * &
                     t319 / 0.2e1_dp
                t565 = gamma_var_s2 ** 2
                t567 = t152 ** 2
                A_s2rho = 0.66725e-1_dp / t565 / t567 * e_c_u_1_s2rho * t149 * t151
                t573 = t145 * t_s2
                t577 = A_s2rho * t156
                t578 = A_s2 * t_s2
                t580 = 0.2e1_dp * t578 * t_s2rho
                t585 = t163 ** 2
                t587 = t159 / t585
                t591 = t156 * t_s2
                t592 = t160 * t591
                t600 = 0.1e1_dp / t168
                H_s2rho = t155 * (0.133450e0_dp * t573 * t165 * t_s2rho + &
                     0.66725e-1_dp * t157 * (t577 + t580) * t164 - 0.66725e-1_dp * t157 *&
                     t587 * (t577 + t580 + 0.2e1_dp * A_s2 * t161 * A_s2rho + 0.4e1_dp *&
                     t592 * t_s2rho)) * t600
                epsilon_cGGA_0_1rho = e_c_u_1_s2rho + H_s2rho
                t603 = epsilon_cGGA - epsilon_cGGA_0_1 .LT. 0.0e0_dp
                IF(t603) THEN
                   mbrho=epsilon_cGGA_0_1rho
                ELSE
                   mbrho=epsilon_cGGArho
                END IF
                t605 = epsilon_cGGA * tau_w
                t609 = tau_w * t241
                epsilon_cRevPKZBrho = epsilon_cGGArho * t244 + 0.106e1_dp * t605 * &
                     t241 * tau_wrho - 0.306e1_dp * t609 * t247 * tau_wrho - 0.153e1_dp *&
                     t242 * (marho / 0.2e1_dp + mbrho / 0.2e1_dp)
                t623 = epsilon_cRevPKZB * t239
                ecrho = epsilon_cRevPKZB * t257 + my_rho * epsilon_cRevPKZBrho * t257 +&
                     t250 * (0.28e1_dp * epsilon_cRevPKZBrho * t251 * t254 + 0.84e1_dp *&
                     t623 * t254 * tau_wrho)

                IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
                   e_rho(ii)=e_rho(ii)+&
                        ecrho + ex_ldarho
                END IF

                tnorm_drho = t214 * t7 / 0.2e1_dp

                Hnorm_drho = gamma_var * (0.133450e0_dp * t362 * t232 * tnorm_drho +&
                     0.133450e0_dp * t217 * t380 * A * tnorm_drho * t231 - 0.66725e-1_dp&
                     * t224 * t376 * (0.2e1_dp * t367 * tnorm_drho + 0.4e1_dp * t381 * &
                     tnorm_drho)) * t390
                epsilon_cGGAnorm_drho = Hnorm_drho
                tau_wnorm_drho = my_ndrho * t7 / 0.4e1_dp
                pnorm_drho = my_ndrho * t1 * t267 / 0.6e1_dp
                znorm_drho = tau_wnorm_drho * t269
                alphanorm_drho = 0.5e1_dp / 0.3e1_dp * pnorm_drho * t271 - 0.5e1_dp &
                     / 0.3e1_dp * t401 * znorm_drho
                tildeq_bnorm_drho = 0.9e1_dp / 0.20e2_dp * alphanorm_drho * t278 - &
                     0.9e1_dp / 0.40e2_dp * t408 * (0.4e0_dp * alphanorm_drho * t273 + &
                     0.4e0_dp * alpha * alphanorm_drho) + 0.2e1_dp / 0.3e1_dp * &
                     pnorm_drho
                t675 = z * znorm_drho
                t677 = p * pnorm_drho
                xnorm_drho = ((0.318192e1_dp * t417 * znorm_drho - 0.636384e1_dp * &
                     t423 * znorm_drho) * p + t288 * pnorm_drho + 0.292e3_dp / &
                     0.2025e4_dp * tildeq_b * tildeq_bnorm_drho - 0.73e2_dp / 0.4050e4_dp&
                     * tildeq_bnorm_drho * t296 - 0.73e2_dp / 0.8100e4_dp * t434 * (&
                     0.36e2_dp * t675 + 0.100e3_dp * t677) + 0.3791437570e-1_dp * t677 + &
                     0.2204014296e0_dp * t675 + 0.101216061e1_dp * t293 * pnorm_drho) * &
                     t307 - 0.2479516082e1_dp * t450 * pnorm_drho
                Fxnorm_drho = 0.1000000000e1_dp * t454 * xnorm_drho
                ex_ldanorm_drho = t315 * Fxnorm_drho
                t_s1norm_drho = t67 * t69 * t7 / 0.2e1_dp
                H_s1norm_drho = t82 * (0.133450e0_dp * t505 * t92 * t_s1norm_drho + &
                     0.133450e0_dp * t72 * t523 * A_s1 * t_s1norm_drho * t91 - &
                     0.66725e-1_dp * t84 * t519 * (0.2e1_dp * t510 * t_s1norm_drho + &
                     0.4e1_dp * t524 * t_s1norm_drho)) * t532
                epsilon_cGGA_1_0norm_drho = H_s1norm_drho
                IF(t535)THEN
                   manorm_drho=H_s1norm_drho
                ELSE
                   manorm_drho=Hnorm_drho
                END IF
                t_s2norm_drho = t140 * t142 * t7 / 0.2e1_dp
                H_s2norm_drho = t155 * (0.133450e0_dp * t573 * t165 * t_s2norm_drho &
                     + 0.133450e0_dp * t145 * t591 * A_s2 * t_s2norm_drho * t164 - &
                     0.66725e-1_dp * t157 * t587 * (0.2e1_dp * t578 * t_s2norm_drho + &
                     0.4e1_dp * t592 * t_s2norm_drho)) * t600
                epsilon_cGGA_0_1norm_drho = H_s2norm_drho
                IF(t603) THEN
                   mbnorm_drho=H_s2norm_drho
                ELSE
                   mbnorm_drho=Hnorm_drho
                END IF
                epsilon_cRevPKZBnorm_drho = Hnorm_drho * t244 + 0.106e1_dp * t605 * &
                     t241 * tau_wnorm_drho - 0.306e1_dp * t609 * t247 * tau_wnorm_drho - &
                     0.153e1_dp * t242 * (manorm_drho / 0.2e1_dp + mbnorm_drho / 0.2e1_dp&
                     )
                ecnorm_drho = my_rho * epsilon_cRevPKZBnorm_drho * t257 + t250 * (&
                     0.28e1_dp * epsilon_cRevPKZBnorm_drho * t251 * t254 + 0.84e1_dp * &
                     t623 * t254 * tau_wnorm_drho)

                IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
                   e_ndrho(ii)=e_ndrho(ii)+&
                        ecnorm_drho + ex_ldanorm_drho
                END IF

                epsilon_cRevPKZBtau = -0.106e1_dp * epsilon_cGGA * t239 * t254 + &
                     0.306e1_dp * t239 * t254 * t247
                t763 = t240 ** 2
                ectau = my_rho * epsilon_cRevPKZBtau * t257 + t250 * (0.28e1_dp * &
                     epsilon_cRevPKZBtau * t251 * t254 - 0.84e1_dp * t252 / t763)
                ztau = -t609
                alphatau = -0.5e1_dp / 0.3e1_dp * t401 * ztau
                tildeq_btau = 0.9e1_dp / 0.20e2_dp * alphatau * t278 - 0.9e1_dp / &
                     0.40e2_dp * t408 * (0.4e0_dp * alphatau * t273 + 0.4e0_dp * alpha * &
                     alphatau)
                t790 = z * ztau
                xtau = ((0.318192e1_dp * t417 * ztau - 0.636384e1_dp * t423 * ztau) &
                     * p + 0.292e3_dp / 0.2025e4_dp * tildeq_b * tildeq_btau - 0.73e2_dp &
                     / 0.4050e4_dp * tildeq_btau * t296 - 0.73e2_dp / 0.225e3_dp * t434 *&
                     t790 + 0.2204014296e0_dp * t790) * t307
                Fxtau = 0.1000000000e1_dp * t454 * xtau
                ex_ldatau = t315 * Fxtau

                IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
                   e_tau(ii)=e_tau(ii)+&
                        ectau + ex_ldatau
                END IF
             END IF
          END IF
       END DO
       IF (cp_debug) THEN
          WRITE(cp_logger_get_default_unit_nr(logger),&
               "('non_coer=',i6,' of ',i6,'=',f8.3,'%')")&
               non_coer,npoints,REAL(non_coer,dp)/REAL(npoints,dp)*100._dp
       END IF
    END IF

  END SUBROUTINE tpss_lda_calc
!***************************************************************************

!!****f* tpss/tpss_lsd_calc *
!!
!!   NAME
!!     tpss_lsd_calc
!!
!!   FUNCTION
!!     low level routine that calculates the tpss functional in the
!!     polarized (lsd) case
!!
!!   NOTES
!!     maple is nice, but if you want the uman readable version of the code
!!     look in doc/tpss.mw
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lsd_calc(rhoa, rhob, norm_drhoa, norm_drhob, norm_drho,&
       tau_a,tau_b,e_0,e_rhoa,e_rhob,e_ndrhoa,e_ndrhob,e_ndrho,e_tau_a,e_tau_b,&
       npoints,grad_deriv,epsilon_rho,epsilon_norm_drho,epsilon_tau,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drhoa, &
                                                norm_drhob, norm_drho, tau_a, &
                                                tau_b
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rhoa, e_rhob, &
                                                e_ndrhoa, e_ndrhob, e_ndrho, &
                                                e_tau_a, e_tau_b
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, epsilon_tau
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: abs_grad_deriv, ii, non_coer
    LOGICAL                                  :: failure, t731, t806
    REAL(kind=dp) :: A, A_s1, A_s1rhoa, A_s2, A_s2rhob, alpha_c, alpha_sp1, &
      alpha_sp1norm_drhoa, alpha_sp1rhoa, alpha_sp1tau_a, alpha_sp2, &
      alpha_sp2norm_drhob, alpha_sp2rhob, alpha_sp2tau_b, Arhoa, Arhob, &
      C_chi, C_chi_eps, C_chi_epsnorm_drho, C_chi_epsnorm_drhoa, &
      C_chi_epsnorm_drhob, C_chi_epsrhoa, C_chi_epsrhob, chi, chirhoa, &
      chirhob, e_c_u_0, e_c_u_0rhoa, e_c_u_0rhob, e_c_u_1_s1, e_c_u_1_s1rhoa, &
      e_c_u_1_s2, e_c_u_1_s2rhob, eps, epsilon_c_unif, epsilon_c_unifrhoa, &
      epsilon_c_unifrhob, epsilon_cGGA, epsilon_cGGA_0_1, epsilon_cGGA_1_0, &
      epsilon_cGGArhoa, epsilon_cGGArhob, epsilon_cRevPKZB, &
      epsilon_cRevPKZBnorm_drho, epsilon_cRevPKZBnorm_drhoa
    REAL(kind=dp) :: epsilon_cRevPKZBnorm_drhob, epsilon_cRevPKZBrhoa, &
      epsilon_cRevPKZBrhob, epsilon_cRevPKZBtau_a, epsilon_cRevPKZBtau_b, &
      ex_unif_sp1, ex_unif_sp2, f, frhoa, frhob, Fx_sp1, Fx_sp2, gamma_var, &
      gamma_var_s1, gamma_var_s2, Hnorm_drho, k_frhoa, k_s, k_s_s1, k_s_s2, &
      ma, manorm_drho, manorm_drhoa, marhoa, marhob, mb, mbnorm_drho, &
      mbnorm_drhob, mbrhoa, mbrhob, mean_tmp, my_norm_drho, my_norm_drhoa, &
      my_norm_drhob, my_rhoa, my_rhob, my_tau_a, my_tau_b, p_sp1, &
      p_sp1norm_drhoa, p_sp1rhoa, p_sp2, p_sp2norm_drhob, p_sp2rhob, phi, &
      phi_s1, phi_s2, phirhoa, phirhob, rho, rs, rs_s1, rs_s1rhoa, rs_s2, &
      rs_s2rhob, rsrhoa
    REAL(kind=dp) :: rsrhob, t, t1, t102, t103, t1039, t104, t1047, t105, &
      t1056, t1068, t1076, t1078, t108, t1084, t1088, t1089, t1091, t1096, &
      t1098, t110, t1102, t1103, t1111, t113, t116, t1171, t1173, t12, t123, &
      t125, t1259, t126, t1261, t127, t129, t13, t133, t137, t1376, t139, &
      t14, t140, t1401, t141, t1415, t1418, t1422, t143, t1447, t147, t15, &
      t150, t152, t158, t16, t162, t167, t170, t171, t18, t184, t188, t189, &
      t19, t192, t194, t195, t196, t197, t2, t20, t200, t202, t203, t204, &
      t206, t207, t21, t210, t211, t212, t213, t214, t215, t216, t218, t219, &
      t22, t220, t223, t224, t226, t230, t233, t235, t241, t245, t25
    REAL(kind=dp) :: t250, t253, t254, t267, t27, t270, t271, t272, t273, &
      t276, t278, t279, t280, t282, t283, t286, t287, t288, t289, t290, t291, &
      t292, t294, t295, t296, t299, t3, t30, t300, t301, t302, t304, t307, &
      t308, t311, t313, t315, t318, t319, t32, t323, t328, t331, t332, t336, &
      t341, t344, t345, t347, t348, t349, t35, t350, t351, t352, t354, t355, &
      t356, t357, t359, t36, t361, t362, t364, t365, t367, t368, t37, t370, &
      t371, t372, t373, t374, t376, t377, t378, t379, t380, t382, t383, t386, &
      t387, t388, t389, t390, t391, t392, t394, t395, t396, t399, t400, t403, &
      t404, t409, t41, t410, t411, t413, t414, t415, t42
    REAL(kind=dp) :: t422, t425, t428, t429, t43, t430, t431, t432, t433, &
      t435, t436, t437, t439, t44, t441, t443, t444, t445, t447, t448, t449, &
      t451, t454, t463, t469, t47, t476, t485, t49, t491, t5, t502, t503, &
      t505, t518, t52, t521, t522, t524, t525, t526, t535, t536, t544, t546, &
      t547, t548, t55, t551, t553, t554, t557, t563, t565, t571, t577, t579, &
      t585, t597, t598, t6, t607, t610, t612, t615, t619, t62, t622, t625, &
      t627, t631, t633, t635, t638, t639, t64, t647, t65, t651, t655, t656, &
      t658, t66, t663, t665, t666, t669, t670, t678, t68, t682, t683, t686, &
      t693, t7, t701, t705, t71, t713, t714, t716, t718, t72, t720
    REAL(kind=dp) :: t722, t732, t74, t740, t749, t761, t769, t771, t777, &
      t78, t781, t782, t784, t789, t791, t795, t796, t80, t804, t81, t816, &
      t82, t822, t825, t83, t831, t837, t850, t856, t86, t863, t872, t878, &
      t88, t889, t890, t892, t9, t905, t908, t909, t91, t911, t912, t923, &
      t926, t928, t93, t96, t97, t970, t98, t997, t999, t_s1, t_s1norm_drhoa, &
      t_s1rhoa, t_s2, t_s2norm_drhob, t_s2rhob, tau, tau_w, tau_w_sp1, &
      tau_w_sp2, tau_wnorm_drho, tau_wrhoa, tau_wrhob, tildeq_b_sp1, &
      tildeq_b_sp1norm_drhoa, tildeq_b_sp1rhoa, tildeq_b_sp1tau_a, &
      tildeq_b_sp2, tildeq_b_sp2norm_drhob, tildeq_b_sp2rhob, &
      tildeq_b_sp2tau_b, tnorm_drho
    REAL(kind=dp) :: trhoa, trhob, z_sp1, z_sp1norm_drhoa, z_sp1rhoa, &
      z_sp1tau_a, z_sp2, z_sp2norm_drhob, z_sp2rhob, z_sp2tau_b
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_error_get_logger(error)

    IF (.NOT. failure) THEN

       abs_grad_deriv=ABS(grad_deriv)
       non_coer=0

       !$omp parallel do default(private) &
       !$omp             shared(rhoa,rhob,tau_a,tau_b,norm_drho,norm_drhoa,&
       !$omp                    norm_drhob,e_0,e_rhoa,e_rhob,e_ndrho,e_ndrhoa,&
       !$omp                    e_ndrhob,e_tau_a,e_tau_b,epsilon_rho,epsilon_tau)
       DO ii=1,npoints
          my_tau_a=MAX(0.0_dp,tau_a(ii))
          my_tau_b=MAX(0.0_dp,tau_b(ii))
          my_rhoa=MAX(0.0_dp,rhoa(ii))
          my_rhob=MAX(0.0_dp,rhob(ii))
          rho=my_rhoa+my_rhob
          tau=my_tau_a+my_tau_b
          IF (rho>epsilon_rho.AND.tau>epsilon_tau) THEN
             my_norm_drhoa=MAX(0.0_dp,norm_drhoa(ii))
             my_norm_drhob=MAX(0.0_dp,norm_drhob(ii))
             my_norm_drho=MAX(0.0_dp,norm_drho(ii))
             my_rhoa=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhoa)
             my_rhob=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhob)
             my_tau_a=MAX(EPSILON(0.0_dp)*1.e4_dp,tau_a(ii))
             my_tau_b=MAX(EPSILON(0.0_dp)*1.e4_dp,tau_b(ii))

             ! enforce 8*rho*tau >= norm_drho^2 (also each spin separately)

             t1 = my_norm_drhoa ** 2
             mean_tmp=8._dp*my_rhoa*my_tau_b
             IF (mean_tmp<t1) THEN
                mean_tmp=0.5_dp*(mean_tmp+t1)
                t1=mean_tmp
                my_norm_drhoa=SQRT(mean_tmp)
                my_tau_a=mean_tmp/(8._dp*my_rhoa)
                non_coer=non_coer+1
             END IF
             t5 = my_norm_drhob ** 2
             mean_tmp=8._dp*my_rhob*my_tau_b
             IF (mean_tmp<t1) THEN
                mean_tmp=0.5_dp*(mean_tmp+t1)
                t1=mean_tmp
                my_norm_drhob=SQRT(mean_tmp)
                my_tau_b=mean_tmp/(8._dp*my_rhob)
                non_coer=non_coer+1
             END IF

             rho=my_rhoa+my_rhob
             t9 = my_tau_a + my_tau_b
             tau = t9
             my_norm_drho=MIN(my_norm_drho,my_norm_drhoa+my_norm_drhob)

             t2 = 0.1e1_dp / my_rhoa
             t3 = t1 * t2
             t6 = 0.1e1_dp / my_rhob
             t7 = t5 * t6

             t12 = 3 ** (0.1e1_dp / 0.3e1_dp)
             t14 = 0.3141592654e1_dp ** 2
             t13 = t1 * t12
             t15 = t14 ** (0.1e1_dp / 0.3e1_dp)
             t16 = t15 ** 2
             t18 = 2 ** (0.1e1_dp / 0.3e1_dp)
             t19 = 0.1e1_dp / t16 * t18
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
             t35 = 0.1e1_dp + 0.4e0_dp * alpha_sp1 * t32
             t36 = SQRT(t35)
             t37 = 0.1e1_dp / t36
             tildeq_b_sp1 = 0.9e1_dp / 0.20e2_dp * t32 * t37 + 0.2e1_dp / &
                  0.3e1_dp * p_sp1
             t41 = z_sp1 ** 2
             t42 = 0.1e1_dp + t41
             t43 = t42 ** 2
             t44 = 0.1e1_dp / t43
             t47 = 0.10e2_dp / 0.81e2_dp + 0.159096e1_dp * t41 * t44
             t49 = tildeq_b_sp1 ** 2
             t52 = p_sp1 ** 2
             t55 = SQRT(0.18e2_dp * t41 + 0.50e2_dp * t52)
             t62 = t47 * p_sp1 + 0.146e3_dp / 0.2025e4_dp * t49 - 0.73e2_dp / &
                  0.4050e4_dp * tildeq_b_sp1 * t55 + 0.1895718785e-1_dp * t52 + &
                  0.1102007148e0_dp * t41 + 0.33738687e0_dp * t52 * p_sp1
             t64 = 0.1e1_dp + 0.1239758041e1_dp * p_sp1
             t65 = t64 ** 2
             t66 = 0.1e1_dp / t65
             t68 = 0.1e1_dp + 0.1243781095e1_dp * t62 * t66
             Fx_sp1 = 0.1804e1_dp - 0.804e0_dp / t68
             t71 = 0.1e1_dp / 0.3141592654e1_dp
             t72 = t71 * t12
             t74 = (t14 * my_rhoa) ** (0.1e1_dp / 0.3e1_dp)
             ex_unif_sp1 = -0.3e1_dp / 0.4e1_dp * t72 * t18 * t74
             t78 = my_rhoa * ex_unif_sp1
             t80 = t5 * t12
             t81 = my_rhob ** 2
             t82 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
             t83 = t82 ** 2
             t86 = t19 / t83 / t81
             p_sp2 = t80 * t86 / 0.24e2_dp
             tau_w_sp2 = t7 / 0.4e1_dp
             t88 = 0.1e1_dp / my_tau_b
             z_sp2 = tau_w_sp2 * t88 / 0.2e1_dp
             t91 = 0.1e1_dp / z_sp2 - 0.1e1_dp
             alpha_sp2 = 0.5e1_dp / 0.3e1_dp * p_sp2 * t91
             t93 = alpha_sp2 - 0.1e1_dp
             t96 = 0.1e1_dp + 0.4e0_dp * alpha_sp2 * t93
             t97 = SQRT(t96)
             t98 = 0.1e1_dp / t97
             tildeq_b_sp2 = 0.9e1_dp / 0.20e2_dp * t93 * t98 + 0.2e1_dp / &
                  0.3e1_dp * p_sp2
             t102 = z_sp2 ** 2
             t103 = 0.1e1_dp + t102
             t104 = t103 ** 2
             t105 = 0.1e1_dp / t104
             t108 = 0.10e2_dp / 0.81e2_dp + 0.159096e1_dp * t102 * t105
             t110 = tildeq_b_sp2 ** 2
             t113 = p_sp2 ** 2
             t116 = SQRT(0.18e2_dp * t102 + 0.50e2_dp * t113)
             t123 = t108 * p_sp2 + 0.146e3_dp / 0.2025e4_dp * t110 - 0.73e2_dp &
                  / 0.4050e4_dp * tildeq_b_sp2 * t116 + 0.1895718785e-1_dp * t113 + &
                  0.1102007148e0_dp * t102 + 0.33738687e0_dp * t113 * p_sp2
             t125 = 0.1e1_dp + 0.1239758041e1_dp * p_sp2
             t126 = t125 ** 2
             t127 = 0.1e1_dp / t126
             t129 = 0.1e1_dp + 0.1243781095e1_dp * t123 * t127
             Fx_sp2 = 0.1804e1_dp - 0.804e0_dp / t129
             t133 = (t14 * my_rhob) ** (0.1e1_dp / 0.3e1_dp)
             ex_unif_sp2 = -0.3e1_dp / 0.4e1_dp * t72 * t18 * t133
             t137 = my_rhob * ex_unif_sp2
             t139 = 4 ** (0.1e1_dp / 0.3e1_dp)
             t140 = t139 ** 2
             t141 = t12 * t140
             t143 = (t71 * t2) ** (0.1e1_dp / 0.3e1_dp)
             rs_s1 = t141 * t143 / 0.4e1_dp
             t147 = SQRT(rs_s1)
             t150 = t147 * rs_s1
             t152 = rs_s1 ** 0.20e1_dp
             t158 = LOG(0.1e1_dp + 0.1608182432e2_dp / (0.75957e1_dp * t147 + &
                  0.35876e1_dp * rs_s1 + 0.16382e1_dp * t150 + 0.49294e0_dp * t152))
             t162 = 0.1e1_dp + 0.20548e0_dp * rs_s1
             t167 = 0.141189e2_dp * t147 + 0.61977e1_dp * rs_s1 + 0.33662e1_dp &
                  * t150 + 0.62517e0_dp * t152
             t170 = 0.1e1_dp + 0.3216468318e2_dp / t167
             t171 = LOG(t170)
             e_c_u_1_s1 = -0.31090e-1_dp * t162 * t171

             t184 = LOG(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t147 + &
                  0.36231e1_dp * rs_s1 + 0.88026e0_dp * t150 + 0.49671e0_dp * t152))
             t188 = 1 / (2 * t18 - 2)
             t189 = LOG(0.2e1_dp)
             gamma_var_s1 = (0.1e1_dp - t189) / t14
             t192 = t18 ** 2
             phi_s1 = t192 / 0.2e1_dp
             t194 = SQRT(t12 * t74 * t71)
             k_s_s1 = 0.2e1_dp * t194
             t195 = 0.1e1_dp / phi_s1
             t196 = my_norm_drhoa * t195
             t197 = 0.1e1_dp / k_s_s1
             t_s1 = t196 * t197 * t2 / 0.2e1_dp
             t200 = 0.1e1_dp / gamma_var_s1
             t202 = phi_s1 ** 2
             t203 = t202 * phi_s1
             t204 = 0.1e1_dp / t203
             t206 = EXP(-e_c_u_1_s1 * t200 * t204)
             t207 = t206 - 0.1e1_dp
             A_s1 = 0.66725e-1_dp * t200 / t207
             t210 = gamma_var_s1 * t203
             t211 = t_s1 ** 2
             t212 = t200 * t211
             t213 = A_s1 * t211
             t214 = 0.1e1_dp + t213
             t215 = A_s1 ** 2
             t216 = t211 ** 2
             t218 = 0.1e1_dp + t213 + t215 * t216
             t219 = 0.1e1_dp / t218
             t220 = t214 * t219
             t223 = 0.1e1_dp + 0.66725e-1_dp * t212 * t220
             t224 = LOG(t223)
             epsilon_cGGA_1_0 = e_c_u_1_s1 + t210 * t224
             t226 = (t71 * t6) ** (0.1e1_dp / 0.3e1_dp)
             rs_s2 = t141 * t226 / 0.4e1_dp
             t230 = SQRT(rs_s2)
             t233 = t230 * rs_s2
             t235 = rs_s2 ** 0.20e1_dp
             t241 = LOG(0.1e1_dp + 0.1608182432e2_dp / (0.75957e1_dp * t230 + &
                  0.35876e1_dp * rs_s2 + 0.16382e1_dp * t233 + 0.49294e0_dp * t235))
             t245 = 0.1e1_dp + 0.20548e0_dp * rs_s2
             t250 = 0.141189e2_dp * t230 + 0.61977e1_dp * rs_s2 + 0.33662e1_dp &
                  * t233 + 0.62517e0_dp * t235
             t253 = 0.1e1_dp + 0.3216468318e2_dp / t250
             t254 = LOG(t253)
             e_c_u_1_s2 = -0.31090e-1_dp * t245 * t254
             t267 = LOG(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t230 + &
                  0.36231e1_dp * rs_s2 + 0.88026e0_dp * t233 + 0.49671e0_dp * t235))
             gamma_var_s2 = gamma_var_s1
             phi_s2 = phi_s1
             t270 = SQRT(t12 * t133 * t71)
             k_s_s2 = 0.2e1_dp * t270
             t271 = 0.1e1_dp / phi_s2
             t272 = my_norm_drhob * t271
             t273 = 0.1e1_dp / k_s_s2
             t_s2 = t272 * t273 * t6 / 0.2e1_dp
             t276 = 0.1e1_dp / gamma_var_s2
             t278 = phi_s2 ** 2
             t279 = t278 * phi_s2
             t280 = 0.1e1_dp / t279
             t282 = EXP(-e_c_u_1_s2 * t276 * t280)
             t283 = t282 - 0.1e1_dp
             A_s2 = 0.66725e-1_dp * t276 / t283
             t286 = gamma_var_s2 * t279
             t287 = t_s2 ** 2
             t288 = t276 * t287
             t289 = A_s2 * t287
             t290 = 0.1e1_dp + t289
             t291 = A_s2 ** 2
             t292 = t287 ** 2
             t294 = 0.1e1_dp + t289 + t291 * t292
             t295 = 0.1e1_dp / t294
             t296 = t290 * t295
             t299 = 0.1e1_dp + 0.66725e-1_dp * t288 * t296
             t300 = LOG(t299)
             epsilon_cGGA_0_1 = e_c_u_1_s2 + t286 * t300
             t301 = my_rhoa - my_rhob
             t302 = 0.1e1_dp / rho
             chi = t301 * t302
             t304 = (t71 * t302) ** (0.1e1_dp / 0.3e1_dp)
             rs = t141 * t304 / 0.4e1_dp
             t307 = 0.1e1_dp + 0.21370e0_dp * rs
             t308 = SQRT(rs)
             t311 = t308 * rs
             t313 = rs ** 0.20e1_dp
             t315 = 0.75957e1_dp * t308 + 0.35876e1_dp * rs + 0.16382e1_dp * &
                  t311 + 0.49294e0_dp * t313
             t318 = 0.1e1_dp + 0.1608182432e2_dp / t315
             t319 = LOG(t318)
             e_c_u_0 = -0.62182e-1_dp * t307 * t319
             t323 = 0.1e1_dp + 0.20548e0_dp * rs
             t328 = 0.141189e2_dp * t308 + 0.61977e1_dp * rs + 0.33662e1_dp * &
                  t311 + 0.62517e0_dp * t313
             t331 = 0.1e1_dp + 0.3216468318e2_dp / t328
             t332 = LOG(t331)
             t336 = 0.1e1_dp + 0.11125e0_dp * rs
             t341 = 0.10357e2_dp * t308 + 0.36231e1_dp * rs + 0.88026e0_dp * &
                  t311 + 0.49671e0_dp * t313
             t344 = 0.1e1_dp + 0.2960857464e1_dp / t341
             t345 = LOG(t344)
             alpha_c = 0.33774e0_dp * t336 * t345
             t347 = 0.1e1_dp + chi
             t348 = t347 ** (0.1e1_dp / 0.3e1_dp)
             t349 = t348 * t347
             t350 = 0.1e1_dp - chi
             t351 = t350 ** (0.1e1_dp / 0.3e1_dp)
             t352 = t351 * t350
             f = (t349 + t352 - 0.2e1_dp) * t188
             t354 = alpha_c * f
             t355 = 0.9e1_dp / 0.8e1_dp / t188
             t356 = chi ** 2
             t357 = t356 ** 2
             t359 = t355 * (0.1e1_dp - t357)
             t361 = -0.31090e-1_dp * t323 * t332 - e_c_u_0
             t362 = t361 * f
             epsilon_c_unif = e_c_u_0 + t354 * t359 + t362 * t357
             gamma_var = gamma_var_s2
             t364 = t348 ** 2
             t365 = t351 ** 2
             phi = t364 / 0.2e1_dp + t365 / 0.2e1_dp
             t367 = t14 * rho
             t368 = t367 ** (0.1e1_dp / 0.3e1_dp)
             t370 = SQRT(t12 * t368 * t71)
             k_s = 0.2e1_dp * t370
             t371 = 0.1e1_dp / phi
             t372 = my_norm_drho * t371
             t373 = 0.1e1_dp / k_s
             t374 = t373 * t302
             t = t372 * t374 / 0.2e1_dp
             t376 = 0.1e1_dp / gamma_var
             t377 = epsilon_c_unif * t376
             t378 = phi ** 2
             t379 = t378 * phi
             t380 = 0.1e1_dp / t379
             t382 = EXP(-t377 * t380)
             t383 = t382 - 0.1e1_dp
             A = 0.66725e-1_dp * t376 / t383
             t386 = gamma_var * t379
             t387 = t ** 2
             t388 = t376 * t387
             t389 = A * t387
             t390 = 0.1e1_dp + t389
             t391 = A ** 2
             t392 = t387 ** 2
             t394 = 0.1e1_dp + t389 + t391 * t392
             t395 = 0.1e1_dp / t394
             t396 = t390 * t395
             t399 = 0.1e1_dp + 0.66725e-1_dp * t388 * t396
             t400 = LOG(t399)
             epsilon_cGGA = epsilon_c_unif + t386 * t400
             t403 = my_rhoa * my_rhob
             t404 = my_norm_drho ** 2
             t409 = SQRT(MAX(0.0_dp,t1 * t81 + t5 * t20 - t403 * t404 + t403 *&
                  t1 + t403 * t5))
             t410 = rho ** 2
             t411 = 0.1e1_dp / t410
             t413 = t12 ** 2
             t414 = 0.2e1_dp * t409 * t411 * t413
             t415 = 0.1e1_dp / t368
             eps = t414 * t415 / 0.6e1_dp
             C_chi = 0.53e0_dp + 0.87e0_dp * t356 + 0.5e0_dp * t357 + &
                  0.226e1_dp * t357 * t356
             tau_w = t404 * t302 / 0.8e1_dp
             t422 = eps ** 2
             t425 = 0.1e1_dp / t349 + 0.1e1_dp / t352
             t428 = 0.1e1_dp + t422 * t425 / 0.2e1_dp
             t429 = t428 ** 2
             t430 = t429 ** 2
             t431 = 0.1e1_dp / t430
             C_chi_eps = C_chi * t431
             ma = MAX(epsilon_cGGA, epsilon_cGGA_1_0)
             mb = MAX(epsilon_cGGA, epsilon_cGGA_0_1)
             t432 = epsilon_cGGA * C_chi_eps
             t433 = tau_w ** 2
             t435 = 0.1e1_dp + C_chi_eps
             t436 = t435 * t433
             t437 = my_rhoa * t302
             t439 = my_rhob * t302
             t441 = t437 * ma + t439 * mb
             t443 = t432 * t433 - t436 * t441
             t444 = tau ** 2
             t445 = 0.1e1_dp / t444
             epsilon_cRevPKZB = epsilon_cGGA + t443 * t445
             t447 = rho * epsilon_cRevPKZB
             t448 = t433 * tau_w
             t449 = epsilon_cRevPKZB * t448
             t451 = 0.1e1_dp / t444 / tau
             t454 = 0.1e1_dp + 0.28e1_dp * t449 * t451

             IF (grad_deriv>=0) THEN
                e_0(ii)=e_0(ii)&
                     + t78 * Fx_sp1 + t137 * Fx_sp2 + t447 * t454
             END IF

             IF (abs_grad_deriv>0) THEN

                p_sp1rhoa = -t13 * t19 / t22 / t20 / my_rhoa / 0.9e1_dp
                t463 = 0.1e1_dp / t20
                z_sp1rhoa = -t1 * t463 * t27 / 0.8e1_dp
                t469 = p_sp1 / t41
                alpha_sp1rhoa = 0.5e1_dp / 0.3e1_dp * p_sp1rhoa * t30 - 0.5e1_dp /&
                     0.3e1_dp * t469 * z_sp1rhoa
                t476 = t32 / t36 / t35
                tildeq_b_sp1rhoa = 0.9e1_dp / 0.20e2_dp * alpha_sp1rhoa * t37 - &
                     0.9e1_dp / 0.40e2_dp * t476 * (0.4e0_dp * alpha_sp1rhoa * t32 + &
                     0.4e0_dp * alpha_sp1 * alpha_sp1rhoa) + 0.2e1_dp / 0.3e1_dp * &
                     p_sp1rhoa
                t485 = z_sp1 * t44
                t491 = t41 * z_sp1 / t43 / t42
                t502 = tildeq_b_sp1 / t55
                t503 = z_sp1 * z_sp1rhoa
                t505 = p_sp1 * p_sp1rhoa
                t518 = t62 / t65 / t64
                t521 = t68 ** 2
                t522 = 0.1e1_dp / t521
                t524 = 0.3141592654e1_dp * t12
                t525 = t74 ** 2
                t526 = 0.1e1_dp / t525
                t535 = t301 * t411
                chirhoa = t302 - t535
                t536 = t304 ** 2
                rsrhoa = -t141 / t536 * t71 * t411 / 0.12e2_dp
                t544 = t315 ** 2
                t546 = t307 / t544
                t547 = 0.1e1_dp / t308
                t548 = t547 * rsrhoa
                t551 = t308 * rsrhoa
                t553 = rs ** 0.10e1_dp
                t554 = t553 * rsrhoa
                t557 = 0.1e1_dp / t318
                e_c_u_0rhoa = -0.1328829340e-1_dp * rsrhoa * t319 + &
                     0.9999999999e0_dp * t546 * (0.3797850000e1_dp * t548 + 0.35876e1_dp &
                     * rsrhoa + 0.2457300000e1_dp * t551 + 0.985880e0_dp * t554) * t557
                t563 = t328 ** 2
                t565 = t323 / t563
                t571 = 0.1e1_dp / t331
                t577 = t341 ** 2
                t579 = t336 / t577
                t585 = 0.1e1_dp / t344
                frhoa = (0.4e1_dp / 0.3e1_dp * t348 * chirhoa - 0.4e1_dp / &
                     0.3e1_dp * t351 * chirhoa) * t188
                t597 = t356 * chi
                t598 = t355 * t597
                t607 = t597 * chirhoa
                epsilon_c_unifrhoa = e_c_u_0rhoa + (0.375735750e-1_dp * rsrhoa * &
                     t345 - 0.9999999999e0_dp * t579 * (0.5178500000e1_dp * t548 + &
                     0.36231e1_dp * rsrhoa + 0.1320390000e1_dp * t551 + 0.993420e0_dp * &
                     t554) * t585) * f * t359 + alpha_c * frhoa * t359 - 0.4e1_dp * t354 &
                     * t598 * chirhoa + (-0.638837320e-2_dp * rsrhoa * t332 + &
                     0.1000000000e1_dp * t565 * (0.7059450000e1_dp * t548 + 0.61977e1_dp &
                     * rsrhoa + 0.5049300000e1_dp * t551 + 0.1250340e1_dp * t554) * t571 &
                     - e_c_u_0rhoa) * f * t357 + t361 * frhoa * t357 + 0.4e1_dp * t362 * &
                     t607
                t610 = 0.1e1_dp / t348
                t612 = 0.1e1_dp / t351
                phirhoa = t610 * chirhoa / 0.3e1_dp - t612 * chirhoa / 0.3e1_dp
                t615 = t368 ** 2
                k_frhoa = t12 / t615 * t14 / 0.3e1_dp
                t619 = 0.1e1_dp / t370
                t622 = my_norm_drho / t378
                t625 = k_s ** 2
                t627 = 0.1e1_dp / t625 * t302
                t631 = t372 * t373 * t411
                trhoa = -t622 * t374 * phirhoa / 0.2e1_dp - t372 * t627 * t619 * &
                     k_frhoa * t71 / 0.2e1_dp - t631 / 0.2e1_dp
                t633 = t383 ** 2
                t635 = t376 / t633
                t638 = t378 ** 2
                t639 = 0.1e1_dp / t638
                Arhoa = -0.66725e-1_dp * t635 * (-epsilon_c_unifrhoa * t376 * t380&
                     + 0.3e1_dp * t377 * t639 * phirhoa) * t382
                t647 = gamma_var * t378
                t651 = t376 * t
                t655 = Arhoa * t387
                t656 = A * t
                t658 = 0.2e1_dp * t656 * trhoa
                t663 = t394 ** 2
                t665 = t390 / t663
                t666 = A * t392
                t669 = t387 * t
                t670 = t391 * t669
                t678 = 0.1e1_dp / t399
                epsilon_cGGArhoa = epsilon_c_unifrhoa + 0.3e1_dp * t647 * t400 * &
                     phirhoa + t386 * (0.133450e0_dp * t651 * t396 * trhoa + &
                     0.66725e-1_dp * t388 * (t655 + t658) * t395 - 0.66725e-1_dp * t388 *&
                     t665 * (t655 + t658 + 0.2e1_dp * t666 * Arhoa + 0.4e1_dp * t670 * &
                     trhoa)) * t678
                t682 = 0.1e1_dp / MAX(t409,1.e4_dp*EPSILON(0.0_dp)) * t411
                t683 = t5 * my_rhoa
                t686 = my_rhob * t1
                t693 = 0.4e1_dp * t409 / t410 / rho
                t701 = t414 / t368 / t367 * t14 / 0.18e2_dp
                t705 = t357 * chi
                tau_wrhoa = -t404 * t411 / 0.8e1_dp
                t713 = C_chi / t430 / t428
                t714 = eps * t425
                t716 = t347 ** 2
                t718 = 0.1e1_dp / t348 / t716
                t720 = t350 ** 2
                t722 = 0.1e1_dp / t351 / t720
                C_chi_epsrhoa = (0.174e1_dp * chi * chirhoa + 0.20e1_dp * t607 + &
                     0.1356e2_dp * t705 * chirhoa) * t431 - 0.4e1_dp * t713 * (t714 * ((&
                     t682 * (0.2e1_dp * t683 - my_rhob * t404 + t686 + my_rhob * t5) - &
                     t693) * t413 * t415 / 0.6e1_dp - t701) + t422 * (-0.4e1_dp / &
                     0.3e1_dp * t718 * chirhoa + 0.4e1_dp / 0.3e1_dp * t722 * chirhoa) / &
                     0.2e1_dp)
                t731 = epsilon_cGGA .LT. epsilon_cGGA_0_1
                IF (t731) THEN
                   mbrhoa = 0._dp
                ELSE
                   mbrhoa = epsilon_cGGArhoa
                END IF
                t732 = t143 ** 2
                rs_s1rhoa = -t141 / t732 * t71 * t463 / 0.12e2_dp
                t740 = t167 ** 2
                t749 = rs_s1 ** 0.10e1_dp
                e_c_u_1_s1rhoa = -0.638837320e-2_dp * rs_s1rhoa * t171 + &
                     0.1000000000e1_dp * t162 / t740 * (0.7059450000e1_dp / t147 * &
                     rs_s1rhoa + 0.61977e1_dp * rs_s1rhoa + 0.5049300000e1_dp * t147 * &
                     rs_s1rhoa + 0.1250340e1_dp * t749 * rs_s1rhoa) / t170
                t761 = k_s_s1 ** 2
                t_s1rhoa = -t196 / t761 * t2 / t194 * t12 * t526 * t14 * t71&
                     / 0.6e1_dp - t196 * t197 * t463 / 0.2e1_dp
                t769 = gamma_var_s1 ** 2
                t771 = t207 ** 2
                A_s1rhoa = 0.66725e-1_dp / t769 / t771 * e_c_u_1_s1rhoa * t204 * t206
                t777 = t200 * t_s1
                t781 = A_s1rhoa * t211
                t782 = A_s1 * t_s1
                t784 = 0.2e1_dp * t782 * t_s1rhoa
                t789 = t218 ** 2
                t791 = t214 / t789
                t795 = t211 * t_s1
                t796 = t215 * t795
                t804 = 0.1e1_dp / t223
                t806 = epsilon_cGGA .LT. epsilon_cGGA_1_0
                IF (t806) THEN
                   marhoa = e_c_u_1_s1rhoa + t210 * (0.133450e0_dp * t777 &
                        * t220 * t_s1rhoa + 0.66725e-1_dp * t212 * (t781 + t784) * t219 - &
                        0.66725e-1_dp * t212 * t791 * (t781 + t784 + 0.2e1_dp * A_s1 * t216 &
                        * A_s1rhoa + 0.4e1_dp * t796 * t_s1rhoa)) * t804
                ELSE
                   marhoa = epsilon_cGGArhoa
                END IF
                t816 = t435 * tau_w
                t822 = my_rhoa * t411 * ma
                t825 = my_rhob * t411 * mb
                epsilon_cRevPKZBrhoa = epsilon_cGGArhoa + (epsilon_cGGArhoa * &
                     C_chi_eps * t433 + epsilon_cGGA * C_chi_epsrhoa * t433 + 0.2e1_dp * &
                     t432 * tau_w * tau_wrhoa - C_chi_epsrhoa * t433 * t441 - 0.2e1_dp * &
                     t816 * t441 * tau_wrhoa - t436 * (t302 * ma - t822 + t437 * marhoa -&
                     t825 + t439 * mbrhoa)) * t445
                t831 = epsilon_cRevPKZB * t454
                t837 = epsilon_cRevPKZB * t433

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_rhoa(ii)=e_rhoa(ii)&
                        + ex_unif_sp1 * Fx_sp1 - my_rhoa * t524 * t18 * t526&
                        * Fx_sp1 / 0.4e1_dp + 0.1000000000e1_dp * t78 * t522 * (((&
                        0.318192e1_dp * t485 * z_sp1rhoa - 0.636384e1_dp * t491 * z_sp1rhoa)&
                        * p_sp1 + t47 * p_sp1rhoa + 0.292e3_dp / 0.2025e4_dp * tildeq_b_sp1&
                        * tildeq_b_sp1rhoa - 0.73e2_dp / 0.4050e4_dp * tildeq_b_sp1rhoa * &
                        t55 - 0.73e2_dp / 0.8100e4_dp * t502 * (0.36e2_dp * t503 + &
                        0.100e3_dp * t505) + 0.3791437570e-1_dp * t505 + 0.2204014296e0_dp *&
                        t503 + 0.101216061e1_dp * t52 * p_sp1rhoa) * t66 - &
                        0.2479516082e1_dp * t518 * p_sp1rhoa) + t831 + rho * &
                        epsilon_cRevPKZBrhoa * t454 + t447 * (0.28e1_dp * &
                        epsilon_cRevPKZBrhoa * t448 * t451 + 0.84e1_dp * t837 * t451 * &
                        tau_wrhoa)
                END IF

                p_sp2rhob = -t80 * t19 / t83 / t81 / my_rhob / 0.9e1_dp
                t850 = 0.1e1_dp / t81
                z_sp2rhob = -t5 * t850 * t88 / 0.8e1_dp
                t856 = p_sp2 / t102
                alpha_sp2rhob = 0.5e1_dp / 0.3e1_dp * p_sp2rhob * t91 - 0.5e1_dp /&
                     0.3e1_dp * t856 * z_sp2rhob
                t863 = t93 / t97 / t96
                tildeq_b_sp2rhob = 0.9e1_dp / 0.20e2_dp * alpha_sp2rhob * t98 - &
                     0.9e1_dp / 0.40e2_dp * t863 * (0.4e0_dp * alpha_sp2rhob * t93 + &
                     0.4e0_dp * alpha_sp2 * alpha_sp2rhob) + 0.2e1_dp / 0.3e1_dp * &
                     p_sp2rhob
                t872 = z_sp2 * t105
                t878 = t102 * z_sp2 / t104 / t103
                t889 = tildeq_b_sp2 / t116
                t890 = z_sp2 * z_sp2rhob
                t892 = p_sp2 * p_sp2rhob
                t905 = t123 / t126 / t125
                t908 = t129 ** 2
                t909 = 0.1e1_dp / t908
                t911 = t133 ** 2
                t912 = 0.1e1_dp / t911
                chirhob = -t302 - t535
                rsrhob = rsrhoa
                t923 = t547 * rsrhob
                t926 = t308 * rsrhob
                t928 = t553 * rsrhob
                e_c_u_0rhob = -0.1328829340e-1_dp * rsrhob * t319 + &
                     0.9999999999e0_dp * t546 * (0.3797850000e1_dp * t923 + 0.35876e1_dp &
                     * rsrhob + 0.2457300000e1_dp * t926 + 0.985880e0_dp * t928) * t557
                frhob = (0.4e1_dp / 0.3e1_dp * t348 * chirhob - 0.4e1_dp / &
                     0.3e1_dp * t351 * chirhob) * t188
                t970 = t597 * chirhob
                epsilon_c_unifrhob = e_c_u_0rhob + (0.375735750e-1_dp * rsrhob * &
                     t345 - 0.9999999999e0_dp * t579 * (0.5178500000e1_dp * t923 + &
                     0.36231e1_dp * rsrhob + 0.1320390000e1_dp * t926 + 0.993420e0_dp * &
                     t928) * t585) * f * t359 + alpha_c * frhob * t359 - 0.4e1_dp * t354 &
                     * t598 * chirhob + (-0.638837320e-2_dp * rsrhob * t332 + &
                     0.1000000000e1_dp * t565 * (0.7059450000e1_dp * t923 + 0.61977e1_dp &
                     * rsrhob + 0.5049300000e1_dp * t926 + 0.1250340e1_dp * t928) * t571 &
                     - e_c_u_0rhob) * f * t357 + t361 * frhob * t357 + 0.4e1_dp * t362 * &
                     t970
                phirhob = t610 * chirhob / 0.3e1_dp - t612 * chirhob / 0.3e1_dp
                trhob = -t622 * t374 * phirhob / 0.2e1_dp - t372 * t627 * t619 * &
                     k_frhoa * t71 / 0.2e1_dp - t631 / 0.2e1_dp
                Arhob = -0.66725e-1_dp * t635 * (-epsilon_c_unifrhob * t376 * t380&
                     + 0.3e1_dp * t377 * t639 * phirhob) * t382
                t997 = Arhob * t387
                t999 = 0.2e1_dp * t656 * trhob
                epsilon_cGGArhob = epsilon_c_unifrhob + 0.3e1_dp * t647 * t400 * &
                     phirhob + t386 * (0.133450e0_dp * t651 * t396 * trhob + &
                     0.66725e-1_dp * t388 * (t997 + t999) * t395 - 0.66725e-1_dp * t388 *&
                     t665 * (t997 + t999 + 0.2e1_dp * t666 * Arhob + 0.4e1_dp * t670 * &
                     trhob)) * t678
                tau_wrhob = tau_wrhoa
                C_chi_epsrhob = (0.174e1_dp * chi * chirhob + 0.20e1_dp * t970 + &
                     0.1356e2_dp * t705 * chirhob) * t431 - 0.4e1_dp * t713 * (t714 * ((&
                     t682 * (0.2e1_dp * t686 - my_rhoa * t404 + my_rhoa * t1 + t683) - &
                     t693) * t413 * t415 / 0.6e1_dp - t701) + t422 * (-0.4e1_dp / &
                     0.3e1_dp * t718 * chirhob + 0.4e1_dp / 0.3e1_dp * t722 * chirhob) / &
                     0.2e1_dp)
                IF (t806) THEN
                   marhob = 0.0_dp
                ELSE
                   marhob = epsilon_cGGArhob
                END IF
                t1039 = t226 ** 2
                rs_s2rhob = -t141 / t1039 * t71 * t850 / 0.12e2_dp
                t1047 = t250 ** 2
                t1056 = rs_s2 ** 0.10e1_dp
                e_c_u_1_s2rhob = -0.638837320e-2_dp * rs_s2rhob * t254 + &
                     0.1000000000e1_dp * t245 / t1047 * (0.7059450000e1_dp / t230 * &
                     rs_s2rhob + 0.61977e1_dp * rs_s2rhob + 0.5049300000e1_dp * t230 * &
                     rs_s2rhob + 0.1250340e1_dp * t1056 * rs_s2rhob) / t253
                t1068 = k_s_s2 ** 2
                t_s2rhob = -t272 / t1068 * t6 / t270 * t12 * t912 * t14 * &
                     t71 / 0.6e1_dp - t272 * t273 * t850 / 0.2e1_dp
                t1076 = gamma_var_s2 ** 2
                t1078 = t283 ** 2
                A_s2rhob = 0.66725e-1_dp / t1076 / t1078 * e_c_u_1_s2rhob * t280 * t282
                t1084 = t276 * t_s2
                t1088 = A_s2rhob * t287
                t1089 = A_s2 * t_s2
                t1091 = 0.2e1_dp * t1089 * t_s2rhob
                t1096 = t294 ** 2
                t1098 = t290 / t1096
                t1102 = t287 * t_s2
                t1103 = t291 * t1102
                t1111 = 0.1e1_dp / t299
                IF (t731) THEN
                   mbrhob = e_c_u_1_s2rhob + t286 * (0.133450e0_dp * t1084&
                        * t296 * t_s2rhob + 0.66725e-1_dp * t288 * (t1088 + t1091) * t295 -&
                        0.66725e-1_dp * t288 * t1098 * (t1088 + t1091 + 0.2e1_dp * A_s2 * &
                        t292 * A_s2rhob + 0.4e1_dp * t1103 * t_s2rhob)) * t1111
                ELSE
                   mbrhob = epsilon_cGGArhob
                END IF
                epsilon_cRevPKZBrhob = epsilon_cGGArhob + (epsilon_cGGArhob * &
                     C_chi_eps * t433 + epsilon_cGGA * C_chi_epsrhob * t433 + 0.2e1_dp * &
                     t432 * tau_w * tau_wrhob - C_chi_epsrhob * t433 * t441 - 0.2e1_dp * &
                     t816 * t441 * tau_wrhob - t436 * (-t822 + t437 * marhob + t302 * mb &
                     - t825 + t439 * mbrhob)) * t445

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_rhob(ii)=e_rhob(ii)&
                        + ex_unif_sp2 * Fx_sp2 - my_rhob * t524 * t18 * t912&
                        * Fx_sp2 / 0.4e1_dp + 0.1000000000e1_dp * t137 * t909 * (((&
                        0.318192e1_dp * t872 * z_sp2rhob - 0.636384e1_dp * t878 * z_sp2rhob)&
                        * p_sp2 + t108 * p_sp2rhob + 0.292e3_dp / 0.2025e4_dp * &
                        tildeq_b_sp2 * tildeq_b_sp2rhob - 0.73e2_dp / 0.4050e4_dp * &
                        tildeq_b_sp2rhob * t116 - 0.73e2_dp / 0.8100e4_dp * t889 * (&
                        0.36e2_dp * t890 + 0.100e3_dp * t892) + 0.3791437570e-1_dp * t892 + &
                        0.2204014296e0_dp * t890 + 0.101216061e1_dp * t113 * p_sp2rhob) * &
                        t127 - 0.2479516082e1_dp * t905 * p_sp2rhob) + t831 + rho * &
                        epsilon_cRevPKZBrhob * t454 + t447 * (0.28e1_dp * &
                        epsilon_cRevPKZBrhob * t448 * t451 + 0.84e1_dp * t837 * t451 * &
                        tau_wrhob)
                END IF

                p_sp1norm_drhoa = my_norm_drhoa * t12 * t25 / 0.12e2_dp
                z_sp1norm_drhoa = my_norm_drhoa * t2 * t27 / 0.4e1_dp
                alpha_sp1norm_drhoa = 0.5e1_dp / 0.3e1_dp * p_sp1norm_drhoa * t30 &
                     - 0.5e1_dp / 0.3e1_dp * t469 * z_sp1norm_drhoa
                tildeq_b_sp1norm_drhoa = 0.9e1_dp / 0.20e2_dp * &
                     alpha_sp1norm_drhoa * t37 - 0.9e1_dp / 0.40e2_dp * t476 * (0.4e0_dp &
                     * alpha_sp1norm_drhoa * t32 + 0.4e0_dp * alpha_sp1 * &
                     alpha_sp1norm_drhoa) + 0.2e1_dp / 0.3e1_dp * p_sp1norm_drhoa
                t1171 = z_sp1 * z_sp1norm_drhoa
                t1173 = p_sp1 * p_sp1norm_drhoa
                t_s1norm_drhoa = t195 * t197 * t2 / 0.2e1_dp
                C_chi_epsnorm_drhoa = -0.2e1_dp / 0.3e1_dp * t713 * t714 * t682 * &
                     (0.2e1_dp * my_norm_drhoa * t81 + 0.2e1_dp * t403 * my_norm_drhoa) * &
                     t413 * t415
                IF (t806) THEN
                   manorm_drhoa =t210 * (0.133450e0_dp * t777 * t220 * &
                     t_s1norm_drhoa + 0.133450e0_dp * t200 * t795 * A_s1 * t_s1norm_drhoa&
                     * t219 - 0.66725e-1_dp * t212 * t791 * (0.2e1_dp * t782 * &
                     t_s1norm_drhoa + 0.4e1_dp * t796 * t_s1norm_drhoa)) * t804
                ELSE
                   manorm_drhoa=0._dp
                END IF
                epsilon_cRevPKZBnorm_drhoa = (epsilon_cGGA * C_chi_epsnorm_drhoa *&
                     t433 - C_chi_epsnorm_drhoa * t433 * t441 - t436 * t437 * &
                     manorm_drhoa) * t445

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_ndrhoa(ii)=e_ndrhoa(ii)&
                        + 0.1000000000e1_dp * t78 * t522 * (((&
                        0.318192e1_dp * t485 * z_sp1norm_drhoa - 0.636384e1_dp * t491 * &
                        z_sp1norm_drhoa) * p_sp1 + t47 * p_sp1norm_drhoa + 0.292e3_dp / &
                        0.2025e4_dp * tildeq_b_sp1 * tildeq_b_sp1norm_drhoa - 0.73e2_dp / &
                        0.4050e4_dp * tildeq_b_sp1norm_drhoa * t55 - 0.73e2_dp / 0.8100e4_dp&
                        * t502 * (0.36e2_dp * t1171 + 0.100e3_dp * t1173) + &
                        0.3791437570e-1_dp * t1173 + 0.2204014296e0_dp * t1171 + &
                        0.101216061e1_dp * t52 * p_sp1norm_drhoa) * t66 - 0.2479516082e1_dp &
                        * t518 * p_sp1norm_drhoa) + rho * epsilon_cRevPKZBnorm_drhoa * t454 &
                        + 0.28e1_dp * t447 * epsilon_cRevPKZBnorm_drhoa * t448 * t451
                END IF

                p_sp2norm_drhob = my_norm_drhob * t12 * t86 / 0.12e2_dp
                z_sp2norm_drhob = my_norm_drhob * t6 * t88 / 0.4e1_dp
                alpha_sp2norm_drhob = 0.5e1_dp / 0.3e1_dp * p_sp2norm_drhob * t91 &
                     - 0.5e1_dp / 0.3e1_dp * t856 * z_sp2norm_drhob
                tildeq_b_sp2norm_drhob = 0.9e1_dp / 0.20e2_dp * &
                     alpha_sp2norm_drhob * t98 - 0.9e1_dp / 0.40e2_dp * t863 * (0.4e0_dp &
                     * alpha_sp2norm_drhob * t93 + 0.4e0_dp * alpha_sp2 * &
                     alpha_sp2norm_drhob) + 0.2e1_dp / 0.3e1_dp * p_sp2norm_drhob
                t1259 = z_sp2 * z_sp2norm_drhob
                t1261 = p_sp2 * p_sp2norm_drhob
                t_s2norm_drhob = t271 * t273 * t6 / 0.2e1_dp
                C_chi_epsnorm_drhob = -0.2e1_dp / 0.3e1_dp * t713 * t714 * t682 * &
                     (0.2e1_dp * my_norm_drhob * t20 + 0.2e1_dp * t403 * my_norm_drhob) * &
                     t413 * t415
                IF (t731) THEN
                   mbnorm_drhob = t286 * (0.133450e0_dp * t1084 * t296 * &
                     t_s2norm_drhob + 0.133450e0_dp * t276 * t1102 * A_s2 * &
                     t_s2norm_drhob * t295 - 0.66725e-1_dp * t288 * t1098 * (0.2e1_dp * &
                     t1089 * t_s2norm_drhob + 0.4e1_dp * t1103 * t_s2norm_drhob)) * &
                     t1111
                ELSE
                   mbnorm_drhob = 0._dp
                END IF
                epsilon_cRevPKZBnorm_drhob = (epsilon_cGGA * C_chi_epsnorm_drhob *&
                     t433 - C_chi_epsnorm_drhob * t433 * t441 - t436 * t439 * &
                     mbnorm_drhob) * t445

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_ndrhob(ii)=e_ndrhob(ii)&
                        + 0.1000000000e1_dp * t137 * t909 * (((&
                        0.318192e1_dp * t872 * z_sp2norm_drhob - 0.636384e1_dp * t878 * &
                        z_sp2norm_drhob) * p_sp2 + t108 * p_sp2norm_drhob + 0.292e3_dp / &
                        0.2025e4_dp * tildeq_b_sp2 * tildeq_b_sp2norm_drhob - 0.73e2_dp / &
                        0.4050e4_dp * tildeq_b_sp2norm_drhob * t116 - 0.73e2_dp / &
                        0.8100e4_dp * t889 * (0.36e2_dp * t1259 + 0.100e3_dp * t1261) + &
                        0.3791437570e-1_dp * t1261 + 0.2204014296e0_dp * t1259 + &
                        0.101216061e1_dp * t113 * p_sp2norm_drhob) * t127 - &
                        0.2479516082e1_dp * t905 * p_sp2norm_drhob) + rho * &
                        epsilon_cRevPKZBnorm_drhob * t454 + 0.28e1_dp * t447 * &
                        epsilon_cRevPKZBnorm_drhob * t448 * t451
                END IF
                tnorm_drho = t371 * t373 * t302 / 0.2e1_dp
                Hnorm_drho = t386 * (0.133450e0_dp * t651 * t396 * tnorm_drho + &
                     0.133450e0_dp * t376 * t669 * A * tnorm_drho * t395 - 0.66725e-1_dp &
                     * t388 * t665 * (0.2e1_dp * t656 * tnorm_drho + 0.4e1_dp * t670 * &
                     tnorm_drho)) * t678
                tau_wnorm_drho = my_norm_drho * t302 / 0.4e1_dp
                C_chi_epsnorm_drho = 0.4e1_dp / 0.3e1_dp * t713 * t714 * t682 * &
                     t403 * my_norm_drho * t413 * t415
                IF (t806) THEN
                   manorm_drho = 0._dp
                ELSE
                   manorm_drho = Hnorm_drho
                END IF
                IF (t731) THEN
                   mbnorm_drho = 0._dp
                ELSE
                   manorm_drho=Hnorm_drho
                END IF
                epsilon_cRevPKZBnorm_drho = Hnorm_drho + (Hnorm_drho * C_chi_eps *&
                     t433 + epsilon_cGGA * C_chi_epsnorm_drho * t433 + 0.2e1_dp * t432 *&
                     tau_w * tau_wnorm_drho - C_chi_epsnorm_drho * t433 * t441 - &
                     0.2e1_dp * t816 * t441 * tau_wnorm_drho - t436 * (t437 * manorm_drho&
                     + t439 * mbnorm_drho)) * t445

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_ndrho(ii)=e_ndrho(ii)&
                        + rho * epsilon_cRevPKZBnorm_drho * t454 + t447 * &
                        (0.28e1_dp * epsilon_cRevPKZBnorm_drho * t448 * t451 + 0.84e1_dp * &
                        t837 * t451 * tau_wnorm_drho)
                END IF

                t1376 = my_tau_a ** 2
                z_sp1tau_a = -tau_w_sp1 / t1376 / 0.2e1_dp
                alpha_sp1tau_a = -0.5e1_dp / 0.3e1_dp * t469 * z_sp1tau_a

                tildeq_b_sp1tau_a = 0.9e1_dp / 0.20e2_dp * alpha_sp1tau_a * t37 - &
                     0.9e1_dp / 0.40e2_dp * t476 * (0.4e0_dp * alpha_sp1tau_a * t32 + &
                     0.4e0_dp * alpha_sp1 * alpha_sp1tau_a)
                t1401 = z_sp1 * z_sp1tau_a
                epsilon_cRevPKZBtau_a = -0.2e1_dp * t443 * t451
                t1415 = t444 ** 2
                t1418 = 0.84e1_dp * t449 / t1415

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_tau_a(ii)=e_tau_a(ii)&
                        + 0.1000000000e1_dp * t78 * t522 * ((0.318192e1_dp * &
                        t485 * z_sp1tau_a - 0.636384e1_dp * t491 * z_sp1tau_a) * p_sp1 + &
                        0.292e3_dp / 0.2025e4_dp * tildeq_b_sp1 * tildeq_b_sp1tau_a - &
                        0.73e2_dp / 0.4050e4_dp * tildeq_b_sp1tau_a * t55 - 0.73e2_dp / &
                        0.225e3_dp * t502 * t1401 + 0.2204014296e0_dp * t1401) * t66 + rho *&
                        epsilon_cRevPKZBtau_a * t454 + t447 * (0.28e1_dp * &
                        epsilon_cRevPKZBtau_a * t448 * t451 - t1418)
                END IF

                t1422 = my_tau_b ** 2
                z_sp2tau_b = -tau_w_sp2 / t1422 / 0.2e1_dp
                alpha_sp2tau_b = -0.5e1_dp / 0.3e1_dp * t856 * z_sp2tau_b
                tildeq_b_sp2tau_b = 0.9e1_dp / 0.20e2_dp * alpha_sp2tau_b * t98 - &
                     0.9e1_dp / 0.40e2_dp * t863 * (0.4e0_dp * alpha_sp2tau_b * t93 + &
                     0.4e0_dp * alpha_sp2 * alpha_sp2tau_b)
                t1447 = z_sp2 * z_sp2tau_b
                epsilon_cRevPKZBtau_b = epsilon_cRevPKZBtau_a

                IF (grad_deriv==-1.OR.grad_deriv>=1) THEN
                   e_tau_b(ii)=e_tau_b(ii)&
                        + 0.1000000000e1_dp * t137 * t909 * ((0.318192e1_dp * &
                        t872 * z_sp2tau_b - 0.636384e1_dp * t878 * z_sp2tau_b) * p_sp2 + &
                        0.292e3_dp / 0.2025e4_dp * tildeq_b_sp2 * tildeq_b_sp2tau_b - &
                        0.73e2_dp / 0.4050e4_dp * tildeq_b_sp2tau_b * t116 - 0.73e2_dp / &
                        0.225e3_dp * t889 * t1447 + 0.2204014296e0_dp * t1447) * t127 + rho &
                        * epsilon_cRevPKZBtau_b * t454 + t447 * (0.28e1_dp * &
                        epsilon_cRevPKZBtau_b * t448 * t451 - t1418)
                END IF
             END IF
          END IF
       END DO

       IF (cp_debug) THEN
          WRITE(cp_logger_get_default_unit_nr(logger),&
               "('non_coer=',i6,' of ',i6,'=',f8.3,'%')")&
               non_coer,2*npoints,REAL(non_coer,dp)/REAL(2*npoints,dp)*100._dp
       END IF
    END IF
  END SUBROUTINE tpss_lsd_calc
!***************************************************************************

END MODULE xc_tpss
