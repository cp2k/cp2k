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
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='xc_tpss'
  
  PUBLIC :: tpss_lda_info, tpss_lda_eval

!***
!****************************************************************************
contains

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
         INTENT(inout), optional                :: needs
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
       needs%rho=.true.
       needs%tau=.true.
       needs%norm_drho=.true.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE tpss_lda_info
  !***************************************************************************

  !!****f* xc_tpss/tpss_lda_eval [1.0] *
  !!
  !!   NAME
  !!     tpss_lda_eval
  !!
  !!   FUNCTION
  !!     evaluates the tpss functional for lda
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
    TYPE(xc_rho_set_type), pointer :: rho_set
    TYPE(xc_derivative_set_type), pointer :: deriv_set
    INTEGER, INTENT(in) :: grad_deriv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='tpss_lda_eval',&
         routineP=moduleN//':'//routineN
    INTEGER :: npoints, handle
    integer, DIMENSION(:,:), pointer :: bo
    REAL(kind=dp),DIMENSION(:,:,:), POINTER :: dummy,rho, norm_drho,&
         tau,e_0,e_rho,e_ndrho,e_tau
    TYPE(xc_derivative_type), POINTER :: deriv
    REAL(kind=dp) :: epsilon_rho, epsilon_norm_drho
    INTEGER :: stat

    call timeset(routineN,"I","",handle)
    failure=.false.
    nullify(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
            tau=tau,drho_cutoff=epsilon_norm_drho,error=error)
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
          call xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
               allocate_deriv=.TRUE.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_tau,error=error)
       END IF
       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          call cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL tpss_lda_calc(rho=rho, norm_drho=norm_drho,&
            tau=tau,e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_tau=e_tau,&
            grad_deriv=grad_deriv, npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_norm_drho=epsilon_norm_drho, error=error)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          nullify(dummy)
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE tpss_lda_eval
  !***************************************************************************

!!****f* xc_tpss/tpss_lda_calc *
!!
!!   NAME
!!     tpss_lda_calc
!!
!!   FUNCTION
!!     low level calculation routine for tpss
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
       npoints,grad_deriv,epsilon_rho,epsilon_norm_drho,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in) :: rho, norm_drho, tau
    REAL(kind=dp), DIMENSION(*), &
         INTENT(inout)                          :: e_0, e_rho, e_ndrho, e_tau
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), INTENT(in)               :: epsilon_rho, epsilon_norm_drho
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, t603, t535
    character(len=*), parameter :: routineN='tpss_lda_calc',&
         routineP=moduleN//':'//routineN
    REAL(kind=dp) :: my_rho,my_ndrho, my_tau
    REAL(kind=dp) :: chi_s1, t1, t2, t3, t4, t5, t6, t7, t9, rs_s1, t14, t17, &
         t19, t25, e_c_u_0_s1, t29, t34, t37, t38, e_c_u_1_s1, t51, &
         alpha_c_s1, f_ii_0_s1, f_s1, epsilon_c_unif_s1, beta_s1, t56, t58, &
         gamma_var_s1, t60, phi_s1, t61, t63, k_f_s1, t66, k_s_s1, t67, t68, &
         t69, t_s1, t72, t74, t75, t76, t78, t79, A_s1, t82, t83, t84, t85, &
         t86, t87, t88, t90, t91, t92, t95, t96, H_s1, epsilon_cGGA_1_0, &
         chi_s2, rs_s2, t99, t102, t104, t110, e_c_u_0_s2, t114, t119, t122, &
         t123, e_c_u_1_s2, t136, alpha_c_s2, f_ii_0_s2, f_s2, &
         epsilon_c_unif_s2, beta_s2, gamma_var_s2, phi_s2, k_f_s2, t139, &
         k_s_s2, t140, t141, t142, t_s2, t145, t147, t148, t149, t151, t152, &
         A_s2, t155, t156, t157, t158, t159, t160, t161, t163, t164, t165, &
         t168, t169, H_s2, epsilon_cGGA_0_1, rs, t172, t173, t176, t178, &
         t180, t183, t184, e_c_u_0, t197, e_c_u_1, t210, alpha_c, f_ii_0, &
         epsilon_c_unif, beta, gamma_var, phi, k_f, t213, k_s, t214, t215, t,&
         t217, t219, t220, A, t223, t224, t225, t226, t227, t228, t230, &
         t231, t232, t235, t236, H, epsilon_cGGA, d, C_chi, t237, tau_w
    real(kind=dp) :: C_chi_eps, ma, mb, t239, t240, t241, t242, t244, t247, &
         epsilon_cRevPKZB, t250, t251, t252, t254, t257, ec, kappa, mu, b, c,&
         e_var, t258, t259, t260, t261, t262, t263, t264, t267, p, t269, z, &
         t271, alpha, t273, t276, t277, t278, tildeq_b, t282, t283, t284, &
         t285, t288, t290, t293, t296, t303, t305, t306, t307, x, t309, Fx, &
         ex_unif, t315, ex_lda, energy, t316, t319, t320, rsrho, t325, t334, &
         e_c_u_0rho, epsilon_c_unifrho, t342, t343, k_frho, k_srho, t348, &
         trho, t355, t357, Arho, t362, t366, t367, t369, t374, t376, t380, &
         t381, t390, Hrho, epsilon_cGGArho, tau_wrho, prho, zrho, t401, &
         alpharho, t408, tildeq_brho, t417, t423, t434, t435, t437, t450, &
         xrho, t453, t454, Fxrho, ex_unifrho, ex_ldarho, rs_s1rho, t468, &
         t477, e_c_u_1_s1rho, epsilon_c_unif_s1rho, k_f_s1rho, k_s_s1rho, &
         t489, t_s1rho, t497, t499, A_s1rho, t505, t509, t510, t512, t517, &
         t519, t523, t524, t532, H_s1rho, epsilon_cGGA_1_0rho, marho
    real(kind=dp) :: rs_s2rho, t538, t547, e_c_u_1_s2rho, epsilon_c_unif_s2rho, &
         k_f_s2rho, k_s_s2rho, t557, t_s2rho, t565, t567, A_s2rho, t573, &
         t577, t578, t580, t585, t587, t591, t592, t600, H_s2rho, &
         epsilon_cGGA_0_1rho, mbrho, t605, t609, epsilon_cRevPKZBrho, &
         t623, ecrho, deriv_rho, tnorm_drho, Hnorm_drho, &
         epsilon_cGGAnorm_drho, tau_wnorm_drho, pnorm_drho, znorm_drho, &
         alphanorm_drho, tildeq_bnorm_drho, t675, t677, xnorm_drho, &
         Fxnorm_drho, ex_ldanorm_drho, t_s1norm_drho, H_s1norm_drho, &
         epsilon_cGGA_1_0norm_drho, manorm_drho, t_s2norm_drho, &
         H_s2norm_drho, epsilon_cGGA_0_1norm_drho, mbnorm_drho, &
         epsilon_cRevPKZBnorm_drho, ecnorm_drho, deriv_norm_drho, &
         epsilon_cRevPKZBtau, t763, ectau, ztau, alphatau, tildeq_btau, t790,&
         xtau, Fxtau, ex_ldatau, deriv_tau
    INTEGER :: ii, abs_grad_deriv

    failure=.false.

    IF (.NOT. failure) THEN
       abs_grad_deriv=abs(grad_deriv)

       DO ii=1,npoints
          my_rho=rho(ii)
          IF (my_rho>epsilon_rho) THEN
             my_ndrho=norm_drho(ii)
             my_tau=tau(ii)

             chi_s1 = 1.0_dp
             t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
             t2 = 4 ** (0.1e1_dp / 0.3e1_dp)
             t3 = t2 ** 2
             t4 = t1 * t3
             t5 = 2 ** (0.1e1_dp / 0.3e1_dp)
             t6 = 0.1e1_dp / 0.3141592654e1_dp
             t7 = 0.1e1_dp / my_rho
             t9 = (t6 * t7) ** (0.1e1_dp / 0.3e1_dp)
             rs_s1 = t4 * t5 * t9 / 0.4e1_dp
             t14 = sqrt(rs_s1)
             t17 = t14 * rs_s1
             t19 = rs_s1 ** 0.20e1_dp
             t25 = log(0.1e1_dp + 0.1608182432e2_dp / (0.75957e1_dp * t14 + &
                  0.35876e1_dp * rs_s1 + 0.16382e1_dp * t17 + 0.49294e0_dp * t19))
             e_c_u_0_s1 = -0.62182e-1_dp * (0.1e1_dp + 0.21370e0_dp * rs_s1) * t25
             t29 = 0.1e1_dp + 0.20548e0_dp * rs_s1
             t34 = 0.141189e2_dp * t14 + 0.61977e1_dp * rs_s1 + 0.33662e1_dp * &
                  t17 + 0.62517e0_dp * t19
             t37 = 0.1e1_dp + 0.3216468318e2_dp / t34
             t38 = log(t37)
             e_c_u_1_s1 = -0.31090e-1_dp * t29 * t38
             t51 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t14 + &
                  0.36231e1_dp * rs_s1 + 0.88026e0_dp * t17 + 0.49671e0_dp * t19))
             alpha_c_s1 = 0.33774e0_dp * (0.1e1_dp + 0.11125e0_dp * rs_s1) * t51
             f_ii_0_s1 = 0.8e1_dp / 0.9e1_dp / (2 * t5 - 2)
             f_s1 = 1
             epsilon_c_unif_s1 = e_c_u_1_s1
             beta_s1 = 0.66725e-1_dp
             t56 = log(0.2e1_dp)
             t58 = 0.3141592654e1_dp ** 2
             gamma_var_s1 = (0.1e1_dp - t56) / t58
             t60 = t5 ** 2
             phi_s1 = t60 / 0.2e1_dp
             t61 = t1 * t60
             t63 = (t58 * my_rho) ** (0.1e1_dp / 0.3e1_dp)
             k_f_s1 = t61 * t63 / 0.2e1_dp
             t66 = sqrt(k_f_s1 * t6)
             k_s_s1 = 0.2e1_dp * t66
             t67 = 0.1e1_dp / phi_s1
             t68 = my_ndrho * t67
             t69 = 0.1e1_dp / k_s_s1
             t_s1 = t68 * t69 * t7 / 0.2e1_dp
             t72 = 0.1e1_dp / gamma_var_s1
             t74 = phi_s1 ** 2
             t75 = t74 * phi_s1
             t76 = 0.1e1_dp / t75
             t78 = exp(-e_c_u_1_s1 * t72 * t76)
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
             t96 = log(t95)
             H_s1 = t82 * t96
             epsilon_cGGA_1_0 = e_c_u_1_s1 + H_s1
             chi_s2 = -1._dp
             rs_s2 = rs_s1
             t99 = sqrt(rs_s2)
             t102 = t99 * rs_s2
             t104 = rs_s2 ** 0.20e1_dp
             t110 = log(0.1e1_dp + 0.1608182432e2_dp / (0.75957e1_dp * t99 + &
                  0.35876e1_dp * rs_s2 + 0.16382e1_dp * t102 + 0.49294e0_dp * t104))
             e_c_u_0_s2 = -0.62182e-1_dp * (0.1e1_dp + 0.21370e0_dp * rs_s2) * t110
             t114 = 0.1e1_dp + 0.20548e0_dp * rs_s2
             t119 = 0.141189e2_dp * t99 + 0.61977e1_dp * rs_s2 + 0.33662e1_dp * &
                  t102 + 0.62517e0_dp * t104
             t122 = 0.1e1_dp + 0.3216468318e2_dp / t119
             t123 = log(t122)
             e_c_u_1_s2 = -0.31090e-1_dp * t114 * t123
             t136 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t99 + &
                  0.36231e1_dp * rs_s2 + 0.88026e0_dp * t102 + 0.49671e0_dp * t104))
             alpha_c_s2 = 0.33774e0_dp * (0.1e1_dp + 0.11125e0_dp * rs_s2) * t136
             f_ii_0_s2 = f_ii_0_s1
             f_s2 = 1
             epsilon_c_unif_s2 = e_c_u_1_s2
             beta_s2 = 0.66725e-1_dp
             gamma_var_s2 = gamma_var_s1
             phi_s2 = phi_s1
             k_f_s2 = k_f_s1
             t139 = sqrt(k_f_s2 * t6)
             k_s_s2 = 0.2e1_dp * t139
             t140 = 0.1e1_dp / phi_s2
             t141 = my_ndrho * t140
             t142 = 0.1e1_dp / k_s_s2
             t_s2 = t141 * t142 * t7 / 0.2e1_dp
             t145 = 0.1e1_dp / gamma_var_s2
             t147 = phi_s2 ** 2
             t148 = t147 * phi_s2
             t149 = 0.1e1_dp / t148
             t151 = exp(-e_c_u_1_s2 * t145 * t149)
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
             t169 = log(t168)
             H_s2 = t155 * t169
             epsilon_cGGA_0_1 = e_c_u_1_s2 + H_s2
             rs = t4 * t9 / 0.4e1_dp
             t172 = 0.1e1_dp + 0.21370e0_dp * rs
             t173 = sqrt(rs)
             t176 = t173 * rs
             t178 = rs ** 0.20e1_dp
             t180 = 0.75957e1_dp * t173 + 0.35876e1_dp * rs + 0.16382e1_dp * t176&
                  + 0.49294e0_dp * t178
             t183 = 0.1e1_dp + 0.1608182432e2_dp / t180
             t184 = log(t183)
             e_c_u_0 = -0.62182e-1_dp * t172 * t184
             t197 = log(0.1e1_dp + 0.3216468318e2_dp / (0.141189e2_dp * t173 + &
                  0.61977e1_dp * rs + 0.33662e1_dp * t176 + 0.62517e0_dp * t178))
             e_c_u_1 = -0.31090e-1_dp * (0.1e1_dp + 0.20548e0_dp * rs) * t197
             t210 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t173 + &
                  0.36231e1_dp * rs + 0.88026e0_dp * t176 + 0.49671e0_dp * t178))
             alpha_c = 0.33774e0_dp * (0.1e1_dp + 0.11125e0_dp * rs) * t210
             f_ii_0 = f_ii_0_s2
             epsilon_c_unif = e_c_u_0
             beta = 0.66725e-1_dp
             gamma_var = gamma_var_s2
             phi = 1
             k_f = t1 * t63
             t213 = sqrt(k_f * t6)
             k_s = 0.2e1_dp * t213
             t214 = 0.1e1_dp / k_s
             t215 = my_ndrho * t214
             t = t215 * t7 / 0.2e1_dp
             t217 = 0.1e1_dp / gamma_var
             t219 = exp(-e_c_u_0 * t217)
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
             t236 = log(t235)
             H = gamma_var * t236
             epsilon_cGGA = e_c_u_0 + H
             d = 0.28e1_dp
             C_chi = 0.53e0_dp
             t237 = my_ndrho ** 2
             tau_w = t237 * t7 / 0.8e1_dp
             C_chi_eps = 0.53e0_dp
             ma = max(epsilon_cGGA_1_0, epsilon_cGGA)
             mb = max(epsilon_cGGA_0_1, epsilon_cGGA)
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
             t277 = sqrt(t276)
             t278 = 0.1e1_dp / t277
             tildeq_b = 0.9e1_dp / 0.20e2_dp * t273 * t278 + 0.2e1_dp / 0.3e1_dp * p
             t282 = z ** 2
             t283 = 0.1e1_dp + t282
             t284 = t283 ** 2
             t285 = 0.1e1_dp / t284
             t288 = 0.10e2_dp / 0.81e2_dp + 0.159096e1_dp * t282 * t285
             t290 = tildeq_b ** 2
             t293 = p ** 2
             t296 = sqrt(0.18e2_dp * t282 + 0.50e2_dp * t293)
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
                t535 = epsilon_cGGA - epsilon_cGGA_1_0 .lt. 0.0e0_dp
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
                t603 = epsilon_cGGA - epsilon_cGGA_0_1 .lt. 0.0e0_dp
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
    END IF

  end subroutine tpss_lda_calc
!***************************************************************************

END MODULE xc_tpss
