!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates the tpss functional.
!> \note
!>      The derivation of the formulaes is lengthly, and not fully trivial,
!>      so I have put it in doc/tpss.mw
!> \par History
!>      05.2004 created
!> \author fawzi
! **************************************************************************************************
MODULE xc_tpss
   USE bibliography,                    ONLY: Tao2003,&
                                              cite_reference
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_get_derivative
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_tpss'

   PUBLIC :: tpss_lda_info, tpss_lda_eval

!***
CONTAINS

! **************************************************************************************************
!> \brief return various information on the functional
!> \param tpss_params ...
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv the highest derivative available
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE tpss_lda_info(tpss_params, reference, shortform, needs, max_deriv)
      TYPE(section_vals_type), POINTER                   :: tpss_params
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      REAL(kind=dp)                                      :: sc, sx

      CALL section_vals_val_get(tpss_params, "SCALE_C", r_val=sc)
      CALL section_vals_val_get(tpss_params, "SCALE_X", r_val=sx)

      IF (PRESENT(reference)) THEN
         IF (sx == 1._dp .AND. sc == 1._dp) THEN
            reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LDA version}"
         ELSE
            WRITE (reference, "(a,'sx=',f5.3,'sc=',f5.3,' {LDA version}')") &
               "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003)", &
               sx, sc
         END IF
      END IF
      IF (PRESENT(shortform)) THEN
         IF (sx == 1._dp .AND. sc == 1._dp) THEN
            shortform = "TPSS meta-GGA functional (LDA)"
         ELSE
            WRITE (shortform, "(a,'sx=',f5.3,'sc=',f5.3,' (LDA)')") &
               "TPSS meta-GGA functional", &
               sx, sc
         END IF
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho = .TRUE.
         needs%tau = .TRUE.
         needs%norm_drho = .TRUE.
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1

   END SUBROUTINE tpss_lda_info

! **************************************************************************************************
!> \brief evaluates the tpss functional in the spin unpolarized (lda) case
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param tpss_params ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE tpss_lda_eval(rho_set, deriv_set, grad_deriv, tpss_params)
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      INTEGER, INTENT(in)                                :: grad_deriv
      TYPE(section_vals_type), POINTER                   :: tpss_params

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tpss_lda_eval'

      INTEGER                                            :: handle, non_coer, npoints
      INTEGER, DIMENSION(:, :), POINTER                  :: bo
      REAL(kind=dp)                                      :: epsilon_rho, epsilon_tau, scale_ec, &
                                                            scale_ex
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: dummy, e_0, e_ndrho, e_rho, e_tau, &
                                                            norm_drho, rho, tau
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)
      NULLIFY (bo)

      CALL cite_reference(tao2003)

      CPASSERT(ASSOCIATED(rho_set))
      CPASSERT(rho_set%ref_count > 0)
      CPASSERT(ASSOCIATED(deriv_set))
      CPASSERT(deriv_set%ref_count > 0)
      CALL xc_rho_set_get(rho_set, rho=rho, &
                          norm_drho=norm_drho, local_bounds=bo, rho_cutoff=epsilon_rho, &
                          tau=tau, tau_cutoff=epsilon_tau)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      dummy => rho

      e_0 => dummy
      e_rho => dummy
      e_ndrho => dummy
      e_tau => dummy

      IF (grad_deriv >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, "", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
      IF (grad_deriv >= 1 .OR. grad_deriv == -1) THEN
         deriv => xc_dset_get_derivative(deriv_set, "(rho)", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho)
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drho)", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
         deriv => xc_dset_get_derivative(deriv_set, "(tau)", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_tau)
      END IF
      IF (grad_deriv > 1 .OR. grad_deriv < -1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF

      non_coer = 0
      CALL section_vals_val_get(tpss_params, "SCALE_C", r_val=scale_ec)
      CALL section_vals_val_get(tpss_params, "SCALE_X", r_val=scale_ex)

!$OMP     PARALLEL DEFAULT(NONE) &
!$OMP              SHARED(rho, tau, norm_drho, e_0, e_rho, e_ndrho, e_tau) &
!$OMP              SHARED(epsilon_rho, epsilon_tau, npoints, grad_deriv) &
!$OMP              SHARED(scale_ec, scale_ex) &
!$OMP              REDUCTION(+: non_coer)

      CALL tpss_lda_calc(rho=rho, norm_drho=norm_drho, &
                         tau=tau, e_0=e_0, e_rho=e_rho, e_ndrho=e_ndrho, e_tau=e_tau, &
                         grad_deriv=grad_deriv, npoints=npoints, epsilon_rho=epsilon_rho, &
                         epsilon_tau=epsilon_tau, scale_ec=scale_ec, scale_ex=scale_ex, non_coer=non_coer)

!$OMP     END PARALLEL

      logger => cp_get_default_logger()

      CALL timestop(handle)
   END SUBROUTINE tpss_lda_eval

! **************************************************************************************************
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
!> \param epsilon_tau ...
!> \param scale_ec ...
!> \param scale_ex ...
!> \param non_coer ...
!> \author fawzi
!> \note
!>      maple is nice, but if you want the uman readable version of the code
!>      look in doc/tpss.mw
! **************************************************************************************************
   SUBROUTINE tpss_lda_calc(rho, norm_drho, tau, e_0, e_rho, e_ndrho, e_tau, &
                            npoints, grad_deriv, epsilon_rho, epsilon_tau, &
                            scale_ec, scale_ex, non_coer)
      REAL(kind=dp), DIMENSION(*), INTENT(in)            :: rho, norm_drho, tau
      REAL(kind=dp), DIMENSION(*), INTENT(inout)         :: e_0, e_rho, e_ndrho, e_tau
      INTEGER, INTENT(in)                                :: npoints, grad_deriv
      REAL(kind=dp), INTENT(in)                          :: epsilon_rho, epsilon_tau, scale_ec, &
                                                            scale_ex
      INTEGER, INTENT(inout)                             :: non_coer

      INTEGER                                            :: abs_grad_deriv, ii
      LOGICAL                                            :: t571, t639
      REAL(kind=dp) :: A, A_1, A_2, A_s1, A_s1rho, A_s2, A_s2rho, alpha, alpha_1_1, alpha_1_2, &
         alphanorm_drho, alpharho, alphatau, Arho, b, beta, beta_1_1, beta_1_2, beta_2_1, &
         beta_2_2, beta_3_1, beta_3_2, beta_4_1, beta_4_2, c, d, e_c_u_0, e_c_u_0rho, e_c_u_1_s1, &
         e_c_u_1_s1rho, e_c_u_1_s2, e_c_u_1_s2rho, e_var, epsilon_cGGA, epsilon_cGGA_0_1, &
         epsilon_cGGA_1_0, epsilon_cGGArho, epsilon_cRevPKZB, epsilon_cRevPKZBnorm_drho, &
         epsilon_cRevPKZBrho, epsilon_cRevPKZBtau, ex_unif, Fx, gamma_var, Hnorm_drho, k_f_s1, &
         k_f_s1rho, k_s, k_s_s1, k_s_s2, kappa, m, ma, manorm_drho, marho, mb, mbnorm_drho, mbrho
      REAL(kind=dp) :: mu, my_ndrho, my_rho, my_tau, p, p_1, p_2, p_3, phi_s1, phi_s2, pnorm_drho, &
         prho, rs, rs_s1, rs_s1rho, rs_s2, rs_s2rho, rsrho, t, t1, t100, t101, t111, t12, t13, &
         t138, t14, t140, t143, t145, t146, t147, t151, t152, t16, t161, t168, t177, t186, t187, &
         t189, t19, t190, t191, t193, t194, t196, t197, t198, t199, t2, t20, t201, t202, t204, &
         t205, t208, t209, t21, t211, t212, t213, t215, t216, t218, t219, t22, t220, t221, t223, &
         t224, t226, t227, t230, t231, t233, t234, t235, t238, t239, t241, t242, t243, t245, t246, &
         t248, t249, t252, t253, t254, t256, t26, t260, t263, t264, t265
      REAL(kind=dp) :: t267, t268, t269, t27, t271, t272, t274, t275, t276, t277, t278, t279, t28, &
         t280, t281, t284, t286, t288, t29, t290, t291, t293, t294, t295, t299, t3, t301, t302, &
         t303, t305, t307, t310, t313, t316, t319, t322, t325, t327, t328, t329, t331, t337, t340, &
         t343, t344, t35, t351, t36, t370, t371, t376, t383, t385, t386, t39, t390, t391, t395, &
         t396, t398, t4, t403, t404, t406, t41, t410, t411, t419, t42, t430, t437, t445, t450, &
         t452, t464, t472, t475, t485, t489, t49, t490, t5, t505, t513, t517, t536, t541, t542, &
         t546, t547, t549, t55, t554, t555, t557, t561, t562, t569, t574
      REAL(kind=dp) :: t58, t585, t6, t60, t604, t609, t610, t614, t615, t617, t622, t623, t625, &
         t629, t630, t637, t642, t645, t659, t67, t7, t71, t73, t77, t78, t79, t799, t80, t84, &
         t85, t89, t9, t94, t95, t96, t_s1, t_s1norm_drho, t_s1rho, t_s2, t_s2norm_drho, t_s2rho, &
         tau_w, tau_wnorm_drho, tau_wrho, tildeq_b, tildeq_bnorm_drho, tildeq_brho, tildeq_btau, &
         tnorm_drho, trho, z, znorm_drho, zrho, ztau

      IF (.FALSE.) THEN
         ! useful for testing, we just hack in a well defined functional of tau, ndrho and rho
         ! and see that things converge properly with OT.
!$OMP        DO
         DO ii = 1, npoints
            my_tau = tau(ii)
            my_rho = rho(ii)
            my_ndrho = norm_drho(ii)
            IF (grad_deriv >= 0) THEN
               e_0(ii) = e_0(ii) + my_tau*my_ndrho*my_rho
            END IF
            IF (grad_deriv >= 1 .OR. grad_deriv == -1) THEN
               e_rho(ii) = e_rho(ii) + my_tau*my_ndrho
               e_ndrho(ii) = e_ndrho(ii) + my_tau*my_rho
               e_tau(ii) = e_tau(ii) + my_rho*my_ndrho
            END IF
         END DO
!$OMP        END DO
         RETURN
      END IF

      abs_grad_deriv = ABS(grad_deriv)

      kappa = 0.804e0_dp
      beta = 0.66725e-1_dp
      mu = 0.21951e0_dp
      gamma_var = (0.1e1_dp - LOG(0.2e1_dp))/pi**2
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

      t1 = 3._dp**(0.1e1_dp/0.3e1_dp)
      t2 = 4._dp**(0.1e1_dp/0.3e1_dp)
      t3 = t2**2
      t4 = t1*t3
      t5 = 2._dp**(0.1e1_dp/0.3e1_dp)
      t6 = 0.1e1_dp/pi
      t12 = t5**2

!$OMP     DO

      DO ii = 1, npoints
         my_tau = tau(ii)
         my_rho = rho(ii)
         IF (my_rho > epsilon_rho .AND. my_tau > epsilon_tau) THEN
            my_ndrho = norm_drho(ii)

            t7 = 0.1e1_dp/my_rho
            t254 = my_ndrho**2
            tau_w = t254*t7/0.8e1_dp

            IF (my_tau < tau_w) THEN
               ! enforce z=norm_rho**2/(8._dp*rho*tau) <1
               m = 0.5_dp*t254 + 4.0_dp*my_rho*my_tau
               my_tau = m/8._dp/my_rho
               my_ndrho = SQRT(m)
               t254 = m
               non_coer = non_coer + 1
            END IF

            t9 = (t6*t7)**(0.1e1_dp/0.3e1_dp)
            rs_s1 = t4*t5*t9/0.4e1_dp
            phi_s1 = t12/0.2e1_dp
            t13 = t1*t12
            t14 = pi**2
            t16 = (t14*my_rho)**(0.1e1_dp/0.3e1_dp)
            k_f_s1 = t13*t16/0.2e1_dp
            t19 = SQRT(k_f_s1*t6)
            k_s_s1 = 0.2e1_dp*t19
            t20 = 0.1e1_dp/phi_s1
            t21 = my_ndrho*t20
            t22 = 0.1e1_dp/k_s_s1
            t_s1 = t21*t22*t7/0.2e1_dp
            rs_s2 = rs_s1
            phi_s2 = phi_s1
            t26 = SQRT(k_f_s1*t6)
            k_s_s2 = 0.2e1_dp*t26
            t27 = 0.1e1_dp/phi_s2
            t28 = my_ndrho*t27
            t29 = 0.1e1_dp/k_s_s2
            t_s2 = t28*t29*t7/0.2e1_dp
            t35 = 0.1e1_dp/A_1
            t36 = SQRT(rs_s2)
            t39 = t36*rs_s2
            t41 = p_1 + 0.1e1_dp
            t42 = rs_s2**t41
            t49 = LOG(0.1e1_dp + t35/(beta_1_1*t36 + beta_2_1*rs_s2 + &
                                      beta_3_1*t39 + beta_4_1*t42)/0.2e1_dp)
            t55 = SQRT(rs_s1)
            t58 = t55*rs_s1
            t60 = rs_s1**t41
            t67 = LOG(0.1e1_dp + t35/(beta_1_1*t55 + beta_2_1*rs_s1 + &
                                      beta_3_1*t58 + beta_4_1*t60)/0.2e1_dp)
            t71 = 0.1e1_dp + alpha_1_2*rs_s2
            t73 = 0.1e1_dp/A_2
            t77 = p_2 + 0.1e1_dp
            t78 = rs_s2**t77
            t79 = beta_4_2*t78
            t80 = beta_1_2*t36 + beta_2_2*rs_s2 + beta_3_2*t39 + t79
            t84 = 0.1e1_dp + t73/t80/0.2e1_dp
            t85 = LOG(t84)
            e_c_u_1_s2 = -0.2e1_dp*A_2*t71*t85
            t89 = 0.1e1_dp + alpha_1_2*rs_s1
            t94 = rs_s1**t77
            t95 = beta_4_2*t94
            t96 = beta_1_2*t55 + beta_2_2*rs_s1 + beta_3_2*t58 + t95
            t100 = 0.1e1_dp + t73/t96/0.2e1_dp
            t101 = LOG(t100)
            e_c_u_1_s1 = -0.2e1_dp*A_2*t89*t101
            t111 = p_3 + 1._dp
            rs = t4*t9/0.4e1_dp
            t138 = 0.1e1_dp + alpha_1_1*rs
            t140 = SQRT(rs)
            t143 = t140*rs
            t145 = rs**t41
            t146 = beta_4_1*t145
            t147 = beta_1_1*t140 + beta_2_1*rs + beta_3_1*t143 + t146
            t151 = 0.1e1_dp + t35/t147/0.2e1_dp
            t152 = LOG(t151)
            e_c_u_0 = -0.2e1_dp*A_1*t138*t152
            t161 = rs**t77
            t168 = LOG(0.1e1_dp + t73/(beta_1_2*t140 + beta_2_2*rs + &
                                       beta_3_2*t143 + beta_4_2*t161)/0.2e1_dp)
            t177 = rs**t111
            t186 = 0.1e1_dp/gamma_var
            t187 = beta*t186
            t189 = phi_s1**2
            t190 = t189*phi_s1
            t191 = 0.1e1_dp/t190
            t193 = EXP(-e_c_u_1_s1*t186*t191)
            t194 = t193 - 0.1e1_dp
            A_s1 = t187/t194
            t196 = gamma_var*t190
            t197 = t_s1**2
            t198 = A_s1*t197
            t199 = 0.1e1_dp + t198
            t201 = A_s1**2
            t202 = t197**2
            t204 = 0.1e1_dp + t198 + t201*t202
            t205 = 0.1e1_dp/t204
            t208 = 0.1e1_dp + t187*t197*t199*t205
            t209 = LOG(t208)
            epsilon_cGGA_1_0 = e_c_u_1_s1 + t196*t209
            t211 = phi_s2**2
            t212 = t211*phi_s2
            t213 = 0.1e1_dp/t212
            t215 = EXP(-e_c_u_1_s2*t186*t213)
            t216 = t215 - 0.1e1_dp
            A_s2 = t187/t216
            t218 = gamma_var*t212
            t219 = t_s2**2
            t220 = A_s2*t219
            t221 = t220 + 0.1e1_dp
            t223 = A_s2**2
            t224 = t219**2
            t226 = 0.1e1_dp + t220 + t223*t224
            t227 = 0.1e1_dp/t226
            t230 = 0.1e1_dp + t187*t219*t221*t227
            t231 = LOG(t230)
            epsilon_cGGA_0_1 = e_c_u_1_s2 + t218*t231
            t233 = SQRT(t1*t16*t6)
            k_s = 0.2e1_dp*t233
            t234 = 0.1e1_dp/k_s
            t235 = my_ndrho*t234
            t = t235*t7/0.2e1_dp
            t238 = EXP(-e_c_u_0*t186)
            t239 = -0.1e1_dp + t238
            A = t187/t239
            t241 = t**2
            t242 = A*t241
            t243 = 0.1e1_dp + t242
            t245 = A**2
            t246 = t241**2
            t248 = 0.1e1_dp + t242 + t245*t246
            t249 = 0.1e1_dp/t248
            t252 = 0.1e1_dp + t187*t241*t243*t249
            t253 = LOG(t252)
            epsilon_cGGA = e_c_u_0 + gamma_var*t253
            ma = MAX(epsilon_cGGA_1_0, epsilon_cGGA)
            mb = MAX(epsilon_cGGA_0_1, epsilon_cGGA)
            t256 = tau_w**2
            t260 = ma/0.2e1_dp + mb/0.2e1_dp
            t263 = 0.53e0_dp*epsilon_cGGA*t256 - 0.153e1_dp*t256*t260
            t264 = my_tau**2
            t265 = 0.1e1_dp/t264
            epsilon_cRevPKZB = epsilon_cGGA + t263*t265
            t267 = my_rho*epsilon_cRevPKZB
            t268 = d*epsilon_cRevPKZB
            t269 = t256*tau_w
            t271 = 0.1e1_dp/t264/my_tau
            t272 = t269*t271
            t274 = 0.1e1_dp + t268*t272
            t275 = t254*t1
            t276 = t14**(0.1e1_dp/0.3e1_dp)
            t277 = t276**2
            t278 = 0.1e1_dp/t277
            t279 = my_rho**2
            t280 = my_rho**(0.1e1_dp/0.3e1_dp)
            t281 = t280**2
            t284 = t278/t281/t279
            p = t275*t284/0.12e2_dp
            t286 = 0.1e1_dp/my_tau
            z = tau_w*t286
            t288 = 0.1e1_dp/z - 0.1e1_dp
            alpha = 0.5e1_dp/0.3e1_dp*p*t288
            t290 = alpha - 0.1e1_dp
            t291 = b*alpha
            t293 = 0.1e1_dp + t291*t290
            t294 = SQRT(t293)
            t295 = 0.1e1_dp/t294
            tildeq_b = 0.9e1_dp/0.20e2_dp*t290*t295 + 0.2e1_dp/0.3e1_dp*p
            t299 = z**2
            t301 = 0.1e1_dp + t299
            t302 = t301**2
            t303 = 0.1e1_dp/t302
            t305 = 0.10e2_dp/0.81e2_dp + c*t299*t303
            t307 = tildeq_b**2
            t310 = p**2
            t313 = SQRT(0.18e2_dp*t299 + 0.50e2_dp*t310)
            t316 = 0.1e1_dp/kappa
            t319 = SQRT(e_var)
            t322 = e_var*mu
            t325 = t305*p + 0.146e3_dp/0.2025e4_dp*t307 - 0.73e2_dp/ &
                   0.4050e4_dp*tildeq_b*t313 + 0.100e3_dp/0.6561e4_dp*t316* &
                   t310 + 0.4e1_dp/0.45e2_dp*t319*t299 + t322*t310*p
            t327 = 0.1e1_dp + t319*p
            t328 = t327**2
            t329 = 0.1e1_dp/t328
            t331 = 0.1e1_dp + t325*t329*t316
            Fx = 0.1e1_dp + kappa - kappa/t331
            ex_unif = -0.3e1_dp/0.4e1_dp*t1*t16*t6
            t337 = my_rho*ex_unif

            IF (grad_deriv >= 0) THEN
               e_0(ii) = e_0(ii) + &
                         scale_ec*t267*t274 + scale_ex*t337*Fx
            END IF

            IF (abs_grad_deriv > 0) THEN
               t340 = t9**2
               t343 = 0.1e1_dp/t279
               t344 = 0.1e1_dp/t340*t6*t343
               rsrho = -t4*t344/0.12e2_dp
               t351 = t147**2
               e_c_u_0rho = -0.2e1_dp*A_1*alpha_1_1*rsrho*t152 + t138/ &
                            t351*(beta_1_1/t140*rsrho/0.2e1_dp + beta_2_1*rsrho + &
                                  0.3e1_dp/0.2e1_dp*beta_3_1*t140*rsrho + t146*t41*rsrho/ &
                                  rs)/t151
               t370 = t16**2
               t371 = 0.1e1_dp/t370
               t376 = k_s**2
               trho = -my_ndrho/t376*t7/t233*t1*t371*t14*t6 &
                      /0.6e1_dp - t235*t343/0.2e1_dp
               t383 = gamma_var**2
               t385 = beta/t383
               t386 = t239**2
               Arho = t385/t386*e_c_u_0rho*t238
               t390 = t187*t
               t391 = t243*t249
               t395 = Arho*t241
               t396 = A*t
               t398 = 0.2e1_dp*t396*trho
               t403 = t187*t241
               t404 = t248**2
               t406 = t243/t404
               t410 = t241*t
               t411 = t245*t410
               t419 = 0.1e1_dp/t252
               epsilon_cGGArho = e_c_u_0rho + gamma_var*(0.2e1_dp*t390*t391 &
                                                         *trho + t187*t241*(t395 + t398)*t249 - t403*t406*(t395 + &
                                                                             t398 + 0.2e1_dp*A*t246*Arho + 0.4e1_dp*t411*trho))*t419
               tau_wrho = -t254*t343/0.8e1_dp
               prho = -0.2e1_dp/0.9e1_dp*t275*t278/t281/t279/my_rho
               zrho = tau_wrho*t286
               t430 = p/t299
               alpharho = 0.5e1_dp/0.3e1_dp*prho*t288 - 0.5e1_dp/0.3e1_dp &
                          *t430*zrho
               t437 = t290/t294/t293
               tildeq_brho = 0.9e1_dp/0.20e2_dp*alpharho*t295 - 0.9e1_dp/ &
                             0.40e2_dp*t437*(b*alpharho*t290 + t291*alpharho) + &
                             0.2e1_dp/0.3e1_dp*prho
               t445 = c*z
               t450 = c*t299*z
               t452 = 0.1e1_dp/t302/t301
               t464 = tildeq_b/t313
               t472 = t316*p
               t475 = t319*z
               t485 = t325/t328/t327
               t489 = t331**2
               t490 = 0.1e1_dp/t489
               rs_s1rho = -t4*t5*t344/0.12e2_dp
               k_f_s1rho = t13*t371*t14/0.6e1_dp
               t505 = k_s_s1**2
               t_s1rho = -t21/t505*t7/t19*k_f_s1rho*t6/0.2e1_dp - t21 &
                         *t22*t343/0.2e1_dp
               t513 = A_2*alpha_1_2
               t517 = t96**2
               e_c_u_1_s1rho = -0.2e1_dp*t513*rs_s1rho*t101 + t89/t517* &
                               (beta_1_2/t55*rs_s1rho/0.2e1_dp + beta_2_2*rs_s1rho + &
                                0.3e1_dp/0.2e1_dp*beta_3_2*t55*rs_s1rho + t95*t77* &
                                rs_s1rho/rs_s1)/t100
               t536 = t194**2
               A_s1rho = t385/t536*e_c_u_1_s1rho*t191*t193
               t541 = t187*t_s1
               t542 = t199*t205
               t546 = A_s1rho*t197
               t547 = A_s1*t_s1
               t549 = 0.2e1_dp*t547*t_s1rho
               t554 = t187*t197
               t555 = t204**2
               t557 = t199/t555
               t561 = t197*t_s1
               t562 = t201*t561
               t569 = 0.1e1_dp/t208
               t571 = epsilon_cGGA .LT. epsilon_cGGA_1_0
               IF (t571) THEN
                  marho = e_c_u_1_s1rho + t196*(0.2e1_dp*t541*t542 &
                                                *t_s1rho + t187*t197*(t546 + t549)*t205 - t554*t557*(t546 &
                                                                              + t549 + 0.2e1_dp*A_s1*t202*A_s1rho + 0.4e1_dp*t562* &
                                                                                                     t_s1rho))*t569
               ELSE
                  marho = epsilon_cGGArho
               END IF
               rs_s2rho = rs_s1rho
               t574 = k_s_s2**2
               t_s2rho = -t28/t574*t7/t26*k_f_s1rho*t6/0.2e1_dp - t28 &
                         *t29*t343/0.2e1_dp
               t585 = t80**2
               e_c_u_1_s2rho = -0.2e1_dp*t513*rs_s2rho*t85 + t71/t585*( &
                               beta_1_2/t36*rs_s2rho/0.2e1_dp + beta_2_2*rs_s2rho + &
                               0.3e1_dp/0.2e1_dp*beta_3_2*t36*rs_s2rho + t79*t77* &
                               rs_s2rho/rs_s2)/t84
               t604 = t216**2
               A_s2rho = t385/t604*e_c_u_1_s2rho*t213*t215
               t609 = t187*t_s2
               t610 = t221*t227
               t614 = A_s2rho*t219
               t615 = A_s2*t_s2
               t617 = 0.2e1_dp*t615*t_s2rho
               t622 = t187*t219
               t623 = t226**2
               t625 = t221/t623
               t629 = t219*t_s2
               t630 = t223*t629
               t637 = 0.1e1_dp/t230
               t639 = epsilon_cGGA .LT. epsilon_cGGA_0_1
               IF (t639) THEN
                  mbrho = e_c_u_1_s2rho + t218*(0.2e1_dp*t609*t610 &
                                                *t_s2rho + t187*t219*(t614 + t617)*t227 - t622*t625*(t614 &
                                                                              + t617 + 0.2e1_dp*A_s2*t224*A_s2rho + 0.4e1_dp*t630* &
                                                                                                     t_s2rho))*t637
               ELSE
                  mbrho = epsilon_cGGArho
               END IF
               t642 = epsilon_cGGA*tau_w
               t645 = tau_w*t260
               epsilon_cRevPKZBrho = epsilon_cGGArho + (0.53e0_dp* &
                                                        epsilon_cGGArho*t256 + 0.106e1_dp*t642*tau_wrho - 0.306e1_dp* &
                                                        t645*tau_wrho - 0.153e1_dp*t256*(marho/0.2e1_dp + mbrho/ &
                                                                                         0.2e1_dp))*t265
               t659 = t256*t271

               IF (grad_deriv >= 1 .OR. grad_deriv == -1) THEN
                  e_rho(ii) = e_rho(ii) + &
                              scale_ec*(epsilon_cRevPKZB*t274 + my_rho* &
                                        epsilon_cRevPKZBrho*t274 + t267*(d*epsilon_cRevPKZBrho*t272 &
                                                                         + 0.3e1_dp*t268*t659*tau_wrho)) + scale_ex*(ex_unif*Fx - &
                                                                                             my_rho*pi*t1*t371*Fx/0.4e1_dp + t337* &
                                                                             t490*(((0.2e1_dp*t445*t303*zrho - 0.4e1_dp*t450*t452* &
                                                                            zrho)*p + t305*prho + 0.292e3_dp/0.2025e4_dp*tildeq_b* &
                                                                            tildeq_brho - 0.73e2_dp/0.4050e4_dp*tildeq_brho*t313 - &
                                                                         0.73e2_dp/0.8100e4_dp*t464*(0.36e2_dp*z*zrho + 0.100e3_dp &
                                                                           *p*prho) + 0.200e3_dp/0.6561e4_dp*t472*prho + 0.8e1_dp/ &
                                                                             0.45e2_dp*t475*zrho + 0.3e1_dp*t322*t310*prho)*t329 - &
                                                                                                           0.2e1_dp*t485*t319*prho))
               END IF

               tnorm_drho = t234*t7/0.2e1_dp
               Hnorm_drho = gamma_var*(0.2e1_dp*t390*t391*tnorm_drho + &
                                       0.2e1_dp*t187*t410*A*tnorm_drho*t249 - t403*t406*( &
                                       0.2e1_dp*t396*tnorm_drho + 0.4e1_dp*t411*tnorm_drho))*t419
               tau_wnorm_drho = my_ndrho*t7/0.4e1_dp
               pnorm_drho = my_ndrho*t1*t284/0.6e1_dp
               znorm_drho = tau_wnorm_drho*t286
               alphanorm_drho = 0.5e1_dp/0.3e1_dp*pnorm_drho*t288 - &
                                0.5e1_dp/0.3e1_dp*t430*znorm_drho
               tildeq_bnorm_drho = 0.9e1_dp/0.20e2_dp*alphanorm_drho*t295 - &
                                   0.9e1_dp/0.40e2_dp*t437*(b*alphanorm_drho*t290 + t291* &
                                                            alphanorm_drho) + 0.2e1_dp/0.3e1_dp*pnorm_drho
               t_s1norm_drho = t20*t22*t7/0.2e1_dp
               IF (t571) THEN
                  manorm_drho = t196*(0.2e1_dp*t541*t542* &
                                      t_s1norm_drho + 0.2e1_dp*t187*t561*A_s1*t_s1norm_drho*t205 &
                                      - t554*t557*(0.2e1_dp*t547*t_s1norm_drho + 0.4e1_dp*t562 &
                                                   *t_s1norm_drho))*t569
               ELSE
                  manorm_drho = Hnorm_drho
               END IF
               t_s2norm_drho = t27*t29*t7/0.2e1_dp
               IF (t639) THEN
                  mbnorm_drho = t218*(0.2e1_dp*t609*t610* &
                                      t_s2norm_drho + 0.2e1_dp*t187*t629*A_s2*t_s2norm_drho*t227 &
                                      - t622*t625*(0.2e1_dp*t615*t_s2norm_drho + 0.4e1_dp*t630 &
                                                   *t_s2norm_drho))*t637
               ELSE
                  mbnorm_drho = Hnorm_drho
               END IF
               epsilon_cRevPKZBnorm_drho = Hnorm_drho + (0.53e0_dp*Hnorm_drho* &
                                                         t256 + 0.106e1_dp*t642*tau_wnorm_drho - 0.306e1_dp*t645* &
                                                         tau_wnorm_drho - 0.153e1_dp*t256*(manorm_drho/0.2e1_dp + &
                                                                                           mbnorm_drho/0.2e1_dp))*t265

               IF (grad_deriv >= 1 .OR. grad_deriv == -1) THEN
                  e_ndrho(ii) = e_ndrho(ii) + &
                                scale_ec*(my_rho*epsilon_cRevPKZBnorm_drho* &
                                          t274 + t267*(d*epsilon_cRevPKZBnorm_drho*t272 + 0.3e1_dp* &
                                                       t268*t659*tau_wnorm_drho)) + scale_ex*t337*t490*((( &
                                                                               0.2e1_dp*t445*t303*znorm_drho - 0.4e1_dp*t450*t452* &
                                                                         znorm_drho)*p + t305*pnorm_drho + 0.292e3_dp/0.2025e4_dp* &
                                                                               tildeq_b*tildeq_bnorm_drho - 0.73e2_dp/0.4050e4_dp* &
                                                                             tildeq_bnorm_drho*t313 - 0.73e2_dp/0.8100e4_dp*t464*( &
                                                                               0.36e2_dp*z*znorm_drho + 0.100e3_dp*p*pnorm_drho) + &
                                                                       0.200e3_dp/0.6561e4_dp*t472*pnorm_drho + 0.8e1_dp/0.45e2_dp &
                                                                          *t475*znorm_drho + 0.3e1_dp*t322*t310*pnorm_drho)*t329 - &
                                                                                                      0.2e1_dp*t485*t319*pnorm_drho)
               END IF

               epsilon_cRevPKZBtau = -0.2e1_dp*t263*t271
               t799 = t264**2
               ztau = -tau_w*t265
               alphatau = -0.5e1_dp/0.3e1_dp*t430*ztau
               tildeq_btau = 0.9e1_dp/0.20e2_dp*alphatau*t295 - 0.9e1_dp/ &
                             0.40e2_dp*t437*(b*alphatau*t290 + t291*alphatau)

               IF (grad_deriv >= 1 .OR. grad_deriv == -1) THEN
                  e_tau(ii) = e_tau(ii) + &
                              scale_ec*(my_rho*epsilon_cRevPKZBtau*t274 + t267* &
                                        (d*epsilon_cRevPKZBtau*t272 - 0.3e1_dp*t268*t269/t799)) + &
                              scale_ex*t337*t490*((0.2e1_dp*t445*t303*ztau - 0.4e1_dp &
                                                   *t450*t452*ztau)*p + 0.292e3_dp/0.2025e4_dp*tildeq_b* &
                                                  tildeq_btau - 0.73e2_dp/0.4050e4_dp*tildeq_btau*t313 - &
                                                  0.73e2_dp/0.225e3_dp*t464*z*ztau + 0.8e1_dp/0.45e2_dp* &
                                                  t475*ztau)*t329
               END IF
            END IF
         END IF
      END DO

!$OMP     END DO

   END SUBROUTINE tpss_lda_calc

END MODULE xc_tpss
