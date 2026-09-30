!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates short range exchange part for wPBE functional and averaged
!>      PBE exchange-hole functional (omega = 0.0 )
!> \par History
!>      Manuel Guidon (05.2007)  : initial version
!> \author Manuel Guidon (05.2007)
! **************************************************************************************************
MODULE xc_xwpbe
   USE bibliography,                    ONLY: Heyd2004,&
                                              cite_reference
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              rootpi
   USE mathlib,                         ONLY: expint
   USE xc_derivative_desc,              ONLY: deriv_norm_drho,&
                                              deriv_norm_drhoa,&
                                              deriv_norm_drhob,&
                                              deriv_rho,&
                                              deriv_rhoa,&
                                              deriv_rhob
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

! *** Global parameters ***

   PUBLIC :: xwpbe_lda_info, xwpbe_lda_eval, xwpbe_lsd_info, &
             xwpbe_lsd_eval

   REAL(KIND=dp), PARAMETER :: alpha1 = -1.128223946706117_dp, &
                               alpha2 = 1.452736265762971_dp, &
                               alpha3 = -1.243162299390327_dp, &
                               alpha4 = 0.971824836115601_dp, &
                               alpha5 = -0.568861079687373_dp, &
                               alpha6 = 0.246880514820192_dp, &
                               alpha7 = -0.065032363850763_dp, &
                               alpha8 = 0.008401793031216_dp
   REAL(KIND=dp), PARAMETER :: beta = 1.455915450052607_dp, &
                               beta2 = 2.0_dp
   REAL(KIND=dp), PARAMETER :: a1 = 0.00979681_dp, &
                               a2 = 0.04108340_dp, &
                               a3 = 0.18744000_dp, &
                               a4 = 0.00120824_dp, &
                               a5 = 0.0347188_dp
   REAL(KIND=dp), PARAMETER :: A = 1.0161144_dp, &
                               B = -0.37170836_dp, &
                               C = -0.077215461_dp, &
                               DD = 0.57786348_dp, &
                               E = -0.051955731_dp, &
                               F1 = 0.47965830_dp, &
                               F2 = 6.4753871_dp, &
                               clda = -0.73855876638202240588423_dp
   REAL(KIND=dp), PARAMETER :: expcutoff = 700.0_dp, &
                               exei1 = 4.0364_dp, &
                               exei2 = 1.15198_dp, &
                               exei3 = 5.03627_dp, &
                               exei4 = 4.19160_dp
   REAL(KIND=dp), PARAMETER :: smax = 8.572844_dp, &
                               sconst = 18.79622316_dp, &
                               scutoff = 8.3_dp
   REAL(KIND=dp), PARAMETER :: gcutoff = 0.08_dp, &
                               g1 = -0.02628417880_dp/E, &
                               g2 = -0.07117647788_dp/E, &
                               g3 = 0.08534541323_dp/E, &
                               g4 = 0.0_dp
   REAL(KIND=dp), PARAMETER :: wcutoff = 14.0_dp
   REAL(KIND=dp), PARAMETER :: f12 = 0.5_dp, f14 = 0.25_dp, f158 = 15.0_dp/8.0_dp, &
                               f1516 = 15.0_dp/16.0_dp, f24364 = 243.0_dp/64.0_dp, &
                               f2716 = 27.0_dp/16.0_dp, f2732 = 27.0_dp/32.0_dp, &
                               f34 = 0.75_dp, f32 = 1.5_dp, f38 = 0.375_dp, f68 = 0.75_dp, &
                               f6561512 = 6561.0_dp/512.0_dp, f8132 = 81.0_dp/32.0_dp, &
                               f8164 = 81.0_dp/64.0_dp, f729128 = 729.0_dp/128.0_dp, &
                               f52 = 2.5_dp, f94 = 9.0_dp/4.0_dp, f916 = 9.0_dp/16.0_dp, &
                               f89 = 8.0_dp/9.0_dp, f2187256 = 2187.0_dp/256.0_dp, &
                               r1 = 1.0_dp, f98 = 9.0_dp/8.0_dp, r15 = 15.0_dp, &
                               r3 = 3.0_dp, r4 = 4.0_dp, r16 = 16.0_dp, r8 = 8.0_dp, &
                               r6 = 6.0_dp, r2 = 2.0_dp

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xwpbe'

CONTAINS

! **************************************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_info(reference, shortform, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "Jochen Heyd and Gustavo E. Scuseria, J. Chem. Phys., 120, 7274 {LDA version}"
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "shortrange part of PBE exchange {LDA}"
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho = .TRUE.
         needs%norm_drho = .TRUE.
      END IF
      ! deriv > 1 are not correct
      IF (PRESENT(max_deriv)) max_deriv = 1
   END SUBROUTINE xwpbe_lda_info

! **************************************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param xwpbe_params input parameters (scaling,omega)
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      The current version provides code for derivatives up to second order.
!>      Using the maple sheet in cp2k/doc it is straightforward to produce routines
!>      for higher derivatives.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_eval(rho_set, deriv_set, order, xwpbe_params)

      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(IN)                                :: order
      TYPE(section_vals_type), POINTER                   :: xwpbe_params

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xwpbe_lda_eval'

      INTEGER                                            :: handle, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp)                                      :: epsilon_norm_drho, epsilon_rho, omega, &
                                                            sx, sx0
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: dummy, e_0, e_ndrho, e_rho, norm_drho, &
                                                            rho
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)

      CALL cite_reference(Heyd2004)

      CALL section_vals_val_get(xwpbe_params, "SCALE_X", r_val=sx)
      CALL section_vals_val_get(xwpbe_params, "SCALE_X0", r_val=sx0)
      CALL section_vals_val_get(xwpbe_params, "OMEGA", r_val=omega)

      CALL xc_rho_set_get(rho_set, rho=rho, &
                          norm_drho=norm_drho, local_bounds=bo, rho_cutoff=epsilon_rho, &
                          drho_cutoff=epsilon_norm_drho)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      dummy => rho

      e_0 => dummy
      e_rho => dummy
      e_ndrho => dummy

      IF (order >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
      IF (order >= 1) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho)
         deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
      END IF
      IF (order > 1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF

!$OMP     PARALLEL DEFAULT(NONE) &
!$OMP     SHARED(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho) &
!$OMP     SHARED(epsilon_rho, sx, sx0, omega)

      CALL xwpbe_lda_calc(npoints, order, rho=rho, norm_drho=norm_drho, &
                          e_0=e_0, e_rho=e_rho, e_ndrho=e_ndrho, &
                          epsilon_rho=epsilon_rho, sx=sx, sx0=sx0, omega=omega)

!$OMP     END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE xwpbe_lda_eval

! **************************************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lda
!> \param npoints ...
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param rho , ndrho: density and norm of the density gradient
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param epsilon_rho ...
!> \param sx , sx0: scaling factor for omega!=0 and omega=0
!> \param sx0 ...
!> \param omega screening parameter
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      In order to avoid numerical instabilities, this routine calls different
!>      subroutines. There are 4 routines for the case omega!=0 and 2 routines
!>      for omega=0.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho, &
                             epsilon_rho, sx, sx0, omega)

      INTEGER, INTENT(in)                                :: npoints, order
      REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: rho, norm_drho, e_0, e_rho, e_ndrho
      REAL(kind=dp), INTENT(in)                          :: epsilon_rho, sx, sx0, omega

      INTEGER                                            :: ip
      REAL(dp)                                           :: my_ndrho, my_rho
      REAL(KIND=dp)                                      :: ss, ss2, sscale, t1, t2, t3, t4, t5, t6, &
                                                            t7, t8, ww

!$OMP     DO

      DO ip = 1, npoints
         my_rho = MAX(rho(ip), 0.0_dp)
         IF (my_rho > epsilon_rho) THEN
            my_ndrho = MAX(norm_drho(ip), 0.0_dp)

            !Do some precalculation in order to catch the correct branch afterwards
            sscale = 1.0_dp
            t1 = pi**2
            t2 = t1*my_rho
            t3 = t2**(0.1e1_dp/0.3e1_dp)
            t4 = 0.1e1_dp/t3
            t5 = omega*t4
            ww = 0.6933612743506347048433524e0_dp*t5
            t6 = my_ndrho*t4
            t7 = 0.1e1_dp/my_rho
            t8 = t7*sscale
            ss = 0.3466806371753173524216762e0_dp*t6*t8
            IF (ss > scutoff) THEN
               ss2 = ss*ss
               sscale = (smax*ss2 - sconst)/(ss2*ss)
            END IF

            IF (sx0 /= 0.0_dp) THEN
               !original PBE hole
               IF (ss*sscale > gcutoff) THEN
                  CALL xwpbe_lda_calc_0(e_0(ip), e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, sscale, sx0, order)
               ELSE
                  CALL xwpbe_lda_calc_01(e_0(ip), e_rho(ip), e_ndrho(ip), my_rho, &
                                         my_ndrho, sscale, sx0, order)
               END IF
            END IF

            IF (sx /= 0.0_dp) THEN
               IF (ww < wcutoff .AND. ss*sscale > gcutoff) THEN
                  CALL xwpbe_lda_calc_1(e_0(ip), e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               ELSE IF (ww < wcutoff .AND. ss*sscale <= gcutoff) THEN
                  CALL xwpbe_lda_calc_2(e_0(ip), e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               ELSE IF (ww >= wcutoff .AND. ss*sscale > gcutoff) THEN
                  CALL xwpbe_lda_calc_3(e_0(ip), e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               ELSE
                  CALL xwpbe_lda_calc_4(e_0(ip), e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               END IF
            END IF
         END IF
      END DO

!$OMP     END DO

   END SUBROUTINE xwpbe_lda_calc

! **************************************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx0 scaling factor
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the exact functional for omega=0.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc_0(e_0, e_rho, e_ndrho, rho, ndrho, sscale, sx0, order)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_0, e_rho, e_ndrho
      REAL(KIND=dp), INTENT(IN)                          :: rho, ndrho, sscale, sx0
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp) :: dQndrho, dQrho, Q, t1, t10, t100, t103, t105, t106, t109, t11, t110, t111, &
         t112, t114, t116, t117, t118, t119, t12, t120, t121, t122, t123, t125, t126, t127, t13, &
         t130, t134, t135, t136, t138, t139, t14, t142, t143, t145, t146, t147, t148, t149, t15, &
         t150, t152, t156, t158, t159, t16, t160, t163, t164, t165, t166, t169, t170, t172, t173, &
         t174, t175, t176, t18, t180, t183, t184, t185, t187, t188, t19, t190, t191, t192, t193, &
         t199, t2, t20, t200, t201, t202, t203, t207, t209, t21, t215, t218, t219, t22, t220, &
         t222, t223, t224, t225, t227, t228, t231, t233, t236, t237, t24
      REAL(KIND=dp) :: t240, t241, t242, t245, t246, t247, t249, t25, t250, t253, t254, t258, t26, &
         t261, t262, t263, t265, t266, t269, t272, t274, t275, t278, t28, t281, t282, t284, t285, &
         t287, t288, t29, t290, t291, t294, t297, t3, t300, t301, t303, t304, t306, t307, t31, &
         t314, t32, t321, t323, t326, t327, t33, t330, t334, t335, t336, t337, t338, t339, t34, &
         t340, t342, t343, t344, t345, t347, t355, t356, t357, t358, t359, t36, t361, t362, t363, &
         t367, t368, t37, t372, t374, t376, t377, t381, t383, t384, t387, t388, t39, t390, t391, &
         t394, t397, t398, t4, t40, t400, t401, t402, t404, t405, t406
      REAL(KIND=dp) :: t408, t409, t410, t411, t412, t413, t414, t415, t416, t417, t418, t42, &
         t423, t425, t426, t429, t430, t432, t433, t435, t44, t440, t443, t444, t446, t447, t449, &
         t46, t463, t465, t47, t471, t472, t475, t476, t478, t479, t480, t484, t485, t489, t491, &
         t495, t497, t5, t500, t501, t504, t505, t507, t508, t510, t511, t513, t514, t515, t52, &
         t520, t53, t54, t56, t58, t59, t6, t61, t63, t65, t66, t67, t68, t69, t7, t70, t71, t72, &
         t73, t74, t76, t77, t78, t79, t8, t80, t81, t83, t84, t86, t88, t89, t9, t91, t92, t93, &
         t94, t95, t97, t98

      t1 = ndrho**2
      t2 = a1*t1
      t3 = r2**2
      t4 = 0.1e1_dp/t3
      t5 = t2*t4
      t6 = pi**2
      t7 = r3*t6
      t8 = t7*rho
      t9 = t8**(0.1e1_dp/0.3e1_dp)
      t10 = t9**2
      t11 = 0.1e1_dp/t10
      t12 = rho**2
      t13 = 0.1e1_dp/t12
      t14 = t11*t13
      t15 = sscale**2
      t16 = t14*t15
      t18 = t1**2
      t19 = a2*t18
      t20 = t3**2
      t21 = 0.1e1_dp/t20
      t22 = t19*t21
      t24 = 0.1e1_dp/t9/t8
      t25 = t12**2
      t26 = 0.1e1_dp/t25
      t28 = t15**2
      t29 = t24*t26*t28
      t31 = t5*t16 + t22*t29
      t32 = f94*t31
      t33 = a3*t18
      t34 = t33*t21
      t36 = t18*ndrho
      t37 = a4*t36
      t39 = 0.1e1_dp/t20/r2
      t40 = t37*t39
      t42 = 0.1e1_dp/t10/t8
      t44 = 0.1e1_dp/t25/rho
      t46 = t28*sscale
      t47 = t42*t44*t46
      t52 = 0.1e1_dp/t20/t3
      t53 = a5*t18*t1*t52
      t54 = r3**2
      t56 = t6**2
      t58 = 0.1e1_dp/t54/t56
      t59 = t25**2
      t61 = t28*t15
      t63 = t58/t59*t61
      t65 = r1 + t34*t29 + t40*t47 + t53*t63
      t66 = 0.1e1_dp/t65
      t67 = t66*t1
      t68 = t32*t67
      t69 = t4*t11
      t70 = t13*t15
      t71 = 0.1e1_dp/A
      t72 = t70*t71
      t73 = t69*t72
      Q = t68*t73
      t74 = rho**(0.1e1_dp/0.3e1_dp)
      t76 = t74*rho*f89
      t77 = B*f12
      t78 = t1*t4
      t79 = t78*t11
      t80 = t31*t66
      t81 = t70*t80
      t83 = t79*t81 + DD
      t84 = 0.1e1_dp/t83
      t86 = F2*t31
      t88 = F1 + t86*t66
      t89 = t70*t88
      t91 = t79*t89 + r1
      t92 = f12*t91
      t93 = t83**2
      t94 = 0.1e1_dp/t93
      t95 = C*t94
      t97 = f34*pi
      t98 = rootpi
      t100 = r6*C
      t103 = r4*B
      t105 = r8*A
      t106 = t93*t83
      t109 = t98*(r15*E + t100*t91*t83 + t103*t93 + t105*t106)
      t110 = 0.1e1_dp/r16
      t111 = SQRT(t83)
      t112 = t111*t106
      t114 = t110/t112
      t116 = SQRT(A)
      t117 = EXP(Q)
      t118 = t116*t117
      t119 = f32*ndrho
      t120 = 0.1e1_dp/r2
      t121 = t119*t120
      t122 = 0.1e1_dp/t9
      t123 = 0.1e1_dp/rho
      t125 = t80*t71
      t126 = SQRT(t125)
      t127 = sscale*t126
      t130 = erfc(t121*t122*t123*t127)
      t134 = 0.1e1_dp/f1516
      t135 = (t97 + t109*t114 - t97*t118*t130)*t134
      t136 = 0.1e1_dp/t98
      t138 = 0.1e1_dp/E
      t139 = t136*t112*t138
      t142 = (-t135*t139 + r1)*E
      t143 = 0.1e1_dp/t106
      t145 = f12*A
      t146 = exei(Q)
      t147 = t78*t14
      t148 = t15*t31
      t149 = t66*t84
      t150 = t148*t149
      t152 = LOG(t147*t150)
      t156 = (t77*t84 + t92*t95 + t142*t143 + t145*(t146 + t152)) &
             *Clda
      e_0 = e_0 + (-t76*t156)*sx0
      IF (order >= 1) THEN
         t158 = t4*t42
         t159 = t2*t158
         t160 = t70*t7
         t163 = t12*rho
         t164 = 0.1e1_dp/t163
         t165 = t11*t164
         t166 = t165*t15
         t169 = t54*t56
         t170 = t169*t12
         t172 = 0.1e1_dp/t9/t170
         t173 = t21*t172
         t174 = t19*t173
         t175 = t26*t28
         t176 = t175*t7
         t180 = t24*t44*t28
         t183 = -0.2e1_dp/0.3e1_dp*t159*t160 - (2._dp*t5*t166) - 0.4e1_dp/ &
                0.3e1_dp*t174*t176 - (4._dp*t22*t180)
         t184 = f94*t183
         t185 = t184*t67
         t187 = t65**2
         t188 = 0.1e1_dp/t187
         t190 = t188*t1*t4
         t191 = t32*t190
         t192 = t15*t71
         t193 = t33*t173
         t199 = 0.1e1_dp/t10/t170
         t200 = t39*t199
         t201 = t37*t200
         t202 = t44*t46
         t203 = t202*t7
         t207 = 0.1e1_dp/t25/t12
         t209 = t42*t207*t46
         t215 = t58/t59/rho*t61
         t218 = -0.4e1_dp/0.3e1_dp*t193*t176 - (4._dp*t34*t180) - 0.5e1_dp &
                /0.3e1_dp*t201*t203 - (5._dp*t40*t209) - (8._dp*t53*t215)
         t219 = t192*t218
         t220 = t14*t219
         t222 = t67*t4
         t223 = t32*t222
         t224 = t42*t13
         t225 = t224*t15
         t227 = t71*r3*t6
         t228 = t225*t227
         t231 = t164*t15
         t233 = t69*t231*t71
         dQrho = t185*t73 - t191*t220 - 0.2e1_dp/0.3e1_dp*t223*t228 - (2._dp &
                                                                       *t68*t233)
         t236 = a1*ndrho
         t237 = t236*t4
         t240 = t1*ndrho
         t241 = a2*t240
         t242 = t241*t21
         t245 = 2._dp*t237*t16 + 4._dp*t242*t29
         t246 = f94*t245
         t247 = t246*t67
         t249 = a3*t240
         t250 = t249*t21
         t253 = a4*t18
         t254 = t253*t39
         t258 = a5*t36*t52
         t261 = 4._dp*t250*t29 + 5._dp*t254*t47 + 6._dp*t258*t63
         t262 = t192*t261
         t263 = t14*t262
         t265 = t66*ndrho
         t266 = t32*t265
         dQndrho = t247*t73 - t191*t263 + 2._dp*t266*t73
         t269 = t74*f89
         t272 = t78*t224
         t274 = t66*r3*t6
         t275 = t148*t274
         t278 = t231*t80
         t281 = t183*t66
         t282 = t70*t281
         t284 = t188*t218
         t285 = t148*t284
         t287 = -0.2e1_dp/0.3e1_dp*t272*t275 - (2._dp*t79*t278) + (t79 &
                                                                   *t282) - t147*t285
         t288 = t94*t287
         t290 = t15*t88
         t291 = t290*t7
         t294 = t231*t88
         t297 = F2*t183
         t300 = t297*t66 - t86*t284
         t301 = t70*t300
         t303 = -0.2e1_dp/0.3e1_dp*t272*t291 - (2._dp*t79*t294) + (t79 &
                                                                   *t301)
         t304 = f12*t303
         t306 = C*t143
         t307 = t306*t287
         t314 = t83*t287
         t321 = t98*(t100*t303*t83 + t100*t91*t287 + 2._dp*t103*t314 &
                     + 3._dp*t105*t93*t287)
         t323 = t93**2
         t326 = t110/t111/t323
         t327 = t326*t287
         t330 = t97*t116
         t334 = rootpi
         t335 = 0.1e1_dp/t334
         t336 = t117*t335
         t337 = f32**2
         t338 = t337*t1
         t339 = t338*t69
         t340 = t70*t125
         t342 = EXP(-t339*t340)
         t343 = t120*t24
         t344 = t119*t343
         t345 = t123*sscale
         t347 = t126*r3*t6
         t355 = t119*t120*t122
         t356 = 0.1e1_dp/t126
         t357 = t281*t71
         t358 = t31*t188
         t359 = t71*t218
         t361 = t357 - t358*t359
         t362 = t356*t361
         t363 = t345*t362
         t367 = t342*(-t344*t345*t347/0.3e1_dp - t121*t122*t13*t127 &
                      + t355*t363/0.2e1_dp)
         t368 = t336*t367
         t372 = (t321*t114 - 0.7e1_dp/0.2e1_dp*t109*t327 - (t330*dQrho &
                                                            *t117*t130) + (2._dp*t330*t368))*t134
         t374 = t135*t136
         t376 = t111*t93*t138
         t377 = t376*t287
         t381 = (-t372*t139 - 0.7e1_dp/0.2e1_dp*t374*t377)*E
         t383 = 0.1e1_dp/t323
         t384 = t383*t287
         t387 = dexeirho(Q, dQrho)
         t388 = t78*t225
         t390 = t84*r3*t6
         t391 = t80*t390
         t394 = t78*t165
         t397 = t15*t183
         t398 = t397*t149
         t400 = t188*t84
         t401 = t400*t218
         t402 = t148*t401
         t404 = t66*t94
         t405 = t404*t287
         t406 = t148*t405
         t408 = -0.2e1_dp/0.3e1_dp*t388*t391 - (2._dp*t394*t150) + t147 &
                *t398 - t147*t402 - t147*t406
         t409 = 0.1e1_dp/t1
         t410 = t408*t409
         t411 = t3*t10
         t412 = t410*t411
         t413 = 0.1e1_dp/t15
         t414 = t12*t413
         t415 = 0.1e1_dp/t31
         t416 = t415*t65
         t417 = t416*t83
         t418 = t414*t417
         t423 = (-t77*t288 + t304*t95 - 2._dp*t92*t307 + t381*t143 - 3._dp &
                 *t142*t384 + t145*(t387 + t412*t418))*Clda
         e_rho = e_rho + (-0.4e1_dp/0.3e1_dp*t269*t156 - t76*t423)*sx0
         t425 = ndrho*t4
         t426 = t425*t11
         t429 = t245*t66
         t430 = t70*t429
         t432 = t188*t261
         t433 = t148*t432
         t435 = 2._dp*t426*t81 + t79*t430 - t147*t433
         t440 = F2*t245
         t443 = t440*t66 - t86*t432
         t444 = t70*t443
         t446 = 2._dp*t426*t89 + t79*t444
         t447 = f12*t446
         t449 = t306*t435
         t463 = t98*(t100*t446*t83 + t100*t91*t435 + 2._dp*t103*t83 &
                     *t435 + 3._dp*t105*t93*t435)
         t465 = t326*t435
         t471 = f32*t120
         t472 = t471*t122
         t475 = t429*t71
         t476 = t71*t261
         t478 = t475 - t358*t476
         t479 = t356*t478
         t480 = t345*t479
         t484 = t342*(t472*t345*t126 + t355*t480/0.2e1_dp)
         t485 = t336*t484
         t489 = (t463*t114 - 0.7e1_dp/0.2e1_dp*t109*t465 - (t330*dQndrho &
                                                            *t117*t130) + (2._dp*t330*t485))*t134
         t491 = t376*t435
         t495 = (-t489*t139 - 0.7e1_dp/0.2e1_dp*t374*t491)*E
         t497 = t383*t435
         t500 = dexeindrho(Q, dQndrho)
         t501 = t425*t14
         t504 = t15*t245
         t505 = t504*t149
         t507 = t400*t261
         t508 = t148*t507
         t510 = t404*t435
         t511 = t148*t510
         t513 = 2._dp*t501*t150 + t147*t505 - t147*t508 - t147*t511
         t514 = t513*t409
         t515 = t514*t411
         t520 = (-t77*t94*t435 + t447*t95 - 2._dp*t92*t449 + t495*t143 &
                 - 3._dp*t142*t497 + t145*(t500 + t515*t418))*Clda
         e_ndrho = e_ndrho + (-t76*t520)*sx0
      END IF

   END SUBROUTINE xwpbe_lda_calc_0

! **************************************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx0 scaling factor
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega=0 using a taylor
!>      expansion for the parameter G.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc_01(e_0, e_rho, e_ndrho, rho, ndrho, sscale, sx0, order)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_0, e_rho, e_ndrho
      REAL(KIND=dp), INTENT(IN)                          :: rho, ndrho, sscale, sx0
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp) :: dQndrho, dQrho, Q, t1, t10, t100, t101, t103, t104, t107, t109, t11, t111, &
         t112, t113, t114, t115, t116, t118, t12, t122, t124, t125, t126, t129, t13, t130, t131, &
         t132, t135, t136, t139, t14, t140, t141, t142, t146, t149, t15, t150, t151, t153, t154, &
         t156, t157, t158, t159, t16, t165, t166, t167, t168, t169, t173, t175, t18, t181, t184, &
         t185, t186, t188, t189, t19, t190, t191, t193, t194, t197, t199, t2, t20, t202, t203, &
         t206, t207, t208, t21, t211, t212, t213, t215, t216, t219, t22, t220, t224, t227, t228, &
         t229, t231, t232, t235, t238, t24, t240, t241, t244, t247, t248
      REAL(KIND=dp) :: t25, t250, t251, t253, t254, t256, t257, t26, t260, t263, t266, t267, t270, &
         t272, t273, t276, t277, t28, t280, t283, t288, t29, t293, t294, t297, t299, t3, t300, &
         t301, t304, t305, t307, t308, t31, t311, t314, t315, t317, t318, t319, t32, t321, t322, &
         t323, t325, t326, t327, t328, t329, t33, t330, t331, t332, t333, t334, t335, t34, t340, &
         t342, t343, t346, t347, t349, t350, t352, t357, t36, t360, t361, t364, t366, t37, t371, &
         t372, t375, t376, t379, t380, t383, t385, t388, t389, t39, t392, t393, t395, t396, t398, &
         t399, t4, t40, t401, t402, t403, t408, t42, t44, t46, t47, t5
      REAL(KIND=dp) :: t52, t53, t54, t56, t58, t59, t6, t61, t63, t65, t66, t67, t68, t69, t7, &
         t70, t71, t72, t73, t74, t76, t77, t78, t79, t8, t80, t81, t83, t84, t86, t88, t89, t9, &
         t92, t93, t94, t95, t97, t98

      t1 = ndrho**2
      t2 = a1*t1
      t3 = r2**2
      t4 = 0.1e1_dp/t3
      t5 = t2*t4
      t6 = pi**2
      t7 = r3*t6
      t8 = t7*rho
      t9 = t8**(0.1e1_dp/0.3e1_dp)
      t10 = t9**2
      t11 = 0.1e1_dp/t10
      t12 = rho**2
      t13 = 0.1e1_dp/t12
      t14 = t11*t13
      t15 = sscale**2
      t16 = t14*t15
      t18 = t1**2
      t19 = a2*t18
      t20 = t3**2
      t21 = 0.1e1_dp/t20
      t22 = t19*t21
      t24 = 0.1e1_dp/t9/t8
      t25 = t12**2
      t26 = 0.1e1_dp/t25
      t28 = t15**2
      t29 = t24*t26*t28
      t31 = t5*t16 + t22*t29
      t32 = f94*t31
      t33 = a3*t18
      t34 = t33*t21
      t36 = t18*ndrho
      t37 = a4*t36
      t39 = 0.1e1_dp/t20/r2
      t40 = t37*t39
      t42 = 0.1e1_dp/t10/t8
      t44 = 0.1e1_dp/t25/rho
      t46 = t28*sscale
      t47 = t42*t44*t46
      t52 = 0.1e1_dp/t20/t3
      t53 = a5*t18*t1*t52
      t54 = r3**2
      t56 = t6**2
      t58 = 0.1e1_dp/t54/t56
      t59 = t25**2
      t61 = t28*t15
      t63 = t58/t59*t61
      t65 = r1 + t34*t29 + t40*t47 + t53*t63
      t66 = 0.1e1_dp/t65
      t67 = t66*t1
      t68 = t32*t67
      t69 = t4*t11
      t70 = t13*t15
      t71 = 0.1e1_dp/A
      t72 = t70*t71
      t73 = t69*t72
      Q = t68*t73
      t74 = rho**(0.1e1_dp/0.3e1_dp)
      t76 = t74*rho*f89
      t77 = B*f12
      t78 = t1*t4
      t79 = t78*t11
      t80 = t31*t66
      t81 = t70*t80
      t83 = t79*t81 + DD
      t84 = 0.1e1_dp/t83
      t86 = F2*t31
      t88 = F1 + t86*t66
      t89 = t70*t88
      t92 = f12*(t79*t89 + r1)
      t93 = t83**2
      t94 = 0.1e1_dp/t93
      t95 = C*t94
      t97 = g2*t1
      t98 = t97*t4
      t100 = g3*t18
      t101 = t100*t21
      t103 = g1 + t98*t16 + t101*t29
      t104 = t70*t103
      t107 = (t79*t104 + r1)*E
      t109 = 0.1e1_dp/t93/t83
      t111 = f12*A
      t112 = exei(Q)
      t113 = t78*t14
      t114 = t15*t31
      t115 = t66*t84
      t116 = t114*t115
      t118 = LOG(t113*t116)
      t122 = (t77*t84 + t92*t95 + t107*t109 + t111*(t112 + t118)) &
             *Clda
      e_0 = e_0 + (-t76*t122)*sx0
      IF (order >= 1) THEN
         t124 = t4*t42
         t125 = t2*t124
         t126 = t70*t7
         t129 = t12*rho
         t130 = 0.1e1_dp/t129
         t131 = t11*t130
         t132 = t131*t15
         t135 = t54*t56
         t136 = t135*t12
         t139 = t21/t9/t136
         t140 = t19*t139
         t141 = t26*t28
         t142 = t141*t7
         t146 = t24*t44*t28
         t149 = -0.2e1_dp/0.3e1_dp*t125*t126 - (2._dp*t5*t132) - 0.4e1_dp/ &
                0.3e1_dp*t140*t142 - (4._dp*t22*t146)
         t150 = f94*t149
         t151 = t150*t67
         t153 = t65**2
         t154 = 0.1e1_dp/t153
         t156 = t154*t1*t4
         t157 = t32*t156
         t158 = t15*t71
         t159 = t33*t139
         t165 = 0.1e1_dp/t10/t136
         t166 = t39*t165
         t167 = t37*t166
         t168 = t44*t46
         t169 = t168*t7
         t173 = 0.1e1_dp/t25/t12
         t175 = t42*t173*t46
         t181 = t58/t59/rho*t61
         t184 = -0.4e1_dp/0.3e1_dp*t159*t142 - (4._dp*t34*t146) - 0.5e1_dp &
                /0.3e1_dp*t167*t169 - (5._dp*t40*t175) - (8._dp*t53*t181)
         t185 = t158*t184
         t186 = t14*t185
         t188 = t67*t4
         t189 = t32*t188
         t190 = t42*t13
         t191 = t190*t15
         t193 = t71*r3*t6
         t194 = t191*t193
         t197 = t130*t15
         t199 = t69*t197*t71
         dQrho = t151*t73 - t157*t186 - 0.2e1_dp/0.3e1_dp*t189*t194 - (2._dp &
                                                                       *t68*t199)
         t202 = a1*ndrho
         t203 = t202*t4
         t206 = t1*ndrho
         t207 = a2*t206
         t208 = t207*t21
         t211 = 2._dp*t203*t16 + 4._dp*t208*t29
         t212 = f94*t211
         t213 = t212*t67
         t215 = a3*t206
         t216 = t215*t21
         t219 = a4*t18
         t220 = t219*t39
         t224 = a5*t36*t52
         t227 = 4._dp*t216*t29 + 5._dp*t220*t47 + 6._dp*t224*t63
         t228 = t158*t227
         t229 = t14*t228
         t231 = t66*ndrho
         t232 = t32*t231
         dQndrho = t213*t73 - t157*t229 + 2._dp*t232*t73
         t235 = t74*f89
         t238 = t78*t190
         t240 = t66*r3*t6
         t241 = t114*t240
         t244 = t197*t80
         t247 = t149*t66
         t248 = t70*t247
         t250 = t154*t184
         t251 = t114*t250
         t253 = -0.2e1_dp/0.3e1_dp*t238*t241 - (2._dp*t79*t244) + (t79 &
                                                                   *t248) - t113*t251
         t254 = t94*t253
         t256 = t15*t88
         t257 = t256*t7
         t260 = t197*t88
         t263 = F2*t149
         t266 = t263*t66 - t86*t250
         t267 = t70*t266
         t270 = f12*(-0.2e1_dp/0.3e1_dp*t238*t257 - (2._dp*t79*t260) + &
                     (t79*t267))
         t272 = C*t109
         t273 = t272*t253
         t276 = t15*t103
         t277 = t276*t7
         t280 = t197*t103
         t283 = t97*t124
         t288 = t100*t139
         t293 = -0.2e1_dp/0.3e1_dp*t283*t126 - (2._dp*t98*t132) - 0.4e1_dp &
                /0.3e1_dp*t288*t142 - (4._dp*t101*t146)
         t294 = t70*t293
         t297 = (-0.2e1_dp/0.3e1_dp*t238*t277 - (2._dp*t79*t280) + (t79 &
                                                                    *t294))*E
         t299 = t93**2
         t300 = 0.1e1_dp/t299
         t301 = t300*t253
         t304 = dexeirho(Q, dQrho)
         t305 = t78*t191
         t307 = t84*r3*t6
         t308 = t80*t307
         t311 = t78*t131
         t314 = t15*t149
         t315 = t314*t115
         t317 = t154*t84
         t318 = t317*t184
         t319 = t114*t318
         t321 = t66*t94
         t322 = t321*t253
         t323 = t114*t322
         t325 = -0.2e1_dp/0.3e1_dp*t305*t308 - (2._dp*t311*t116) + t113 &
                *t315 - t113*t319 - t113*t323
         t326 = 0.1e1_dp/t1
         t327 = t325*t326
         t328 = t3*t10
         t329 = t327*t328
         t330 = 0.1e1_dp/t15
         t331 = t12*t330
         t332 = 0.1e1_dp/t31
         t333 = t332*t65
         t334 = t333*t83
         t335 = t331*t334
         t340 = (-t77*t254 + t270*t95 - 2._dp*t92*t273 + t297*t109 - 3._dp &
                 *t107*t301 + t111*(t304 + t329*t335))*Clda
         e_rho = e_rho + (-0.4e1_dp/0.3e1_dp*t235*t122 - t76*t340)*sx0
         t342 = ndrho*t4
         t343 = t342*t11
         t346 = t211*t66
         t347 = t70*t346
         t349 = t154*t227
         t350 = t114*t349
         t352 = 2._dp*t343*t81 + t79*t347 - t113*t350
         t357 = F2*t211
         t360 = t357*t66 - t86*t349
         t361 = t70*t360
         t364 = f12*(2._dp*t343*t89 + t79*t361)
         t366 = t272*t352
         t371 = g2*ndrho
         t372 = t371*t4
         t375 = g3*t206
         t376 = t375*t21
         t379 = 2._dp*t372*t16 + 4._dp*t376*t29
         t380 = t70*t379
         t383 = (2._dp*t343*t104 + t79*t380)*E
         t385 = t300*t352
         t388 = dexeindrho(Q, dQndrho)
         t389 = t342*t14
         t392 = t15*t211
         t393 = t392*t115
         t395 = t317*t227
         t396 = t114*t395
         t398 = t321*t352
         t399 = t114*t398
         t401 = 2._dp*t389*t116 + t113*t393 - t113*t396 - t113*t399
         t402 = t401*t326
         t403 = t402*t328
         t408 = (-t77*t94*t352 + t364*t95 - 2._dp*t92*t366 + t383*t109 &
                 - 3._dp*t107*t385 + t111*(t388 + t403*t335))*Clda
         e_ndrho = e_ndrho + (-t76*t408)*sx0
      END IF

   END SUBROUTINE xwpbe_lda_calc_01

! **************************************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega scaling factor
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the exact functional for omega!=0.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc_1(e_0, e_rho, e_ndrho, rho, ndrho, omega, sscale, sx, order)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_0, e_rho, e_ndrho
      REAL(KIND=dp), INTENT(IN)                          :: rho, ndrho, omega, sscale, sx
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp) :: dQndrho, dQrho, Q, t1, t10, t1001, t1009, t1011, t102, t1024, t1035, t104, &
         t105, t1052, t1066, t1068, t1069, t108, t109, t11, t110, t111, t113, t115, t116, t117, &
         t118, t119, t12, t120, t121, t123, t124, t125, t126, t127, t128, t129, t13, t131, t132, &
         t133, t136, t14, t140, t141, t142, t144, t145, t147, t148, t149, t15, t150, t151, t152, &
         t153, t154, t155, t156, t157, t158, t159, t16, t162, t163, t166, t167, t168, t169, t17, &
         t170, t174, t176, t177, t18, t180, t181, t185, t186, t189, t19, t190, t192, t193, t194, &
         t195, t197, t198, t199, t2, t200, t202, t205, t209, t21, t210
      REAL(KIND=dp) :: t213, t214, t219, t22, t221, t222, t223, t224, t227, t228, t23, t233, t237, &
         t239, t24, t241, t243, t246, t247, t25, t250, t251, t252, t253, t255, t256, t257, t258, &
         t263, t264, t265, t267, t27, t270, t271, t273, t275, t276, t28, t280, t284, t285, t286, &
         t288, t289, t29, t290, t293, t294, t295, t298, t3, t303, t305, t307, t31, t310, t311, &
         t312, t313, t315, t316, t318, t319, t32, t321, t326, t329, t330, t332, t333, t335, t336, &
         t338, t339, t34, t340, t342, t343, t345, t346, t347, t348, t349, t35, t351, t352, t353, &
         t354, t357, t358, t36, t361, t362, t363, t364, t368, t371, t372
      REAL(KIND=dp) :: t373, t374, t375, t376, t377, t378, t38, t383, t384, t385, t386, t39, t390, &
         t392, t398, t4, t401, t402, t403, t404, t406, t409, t41, t411, t412, t415, t416, t419, &
         t42, t420, t421, t424, t425, t426, t428, t429, t432, t433, t437, t44, t440, t441, t442, &
         t444, t446, t449, t452, t453, t454, t457, t46, t461, t463, t465, t468, t469, t472, t475, &
         t478, t479, t48, t481, t486, t49, t493, t495, t498, t499, t5, t502, t503, t504, t507, &
         t508, t509, t510, t511, t513, t514, t515, t517, t518, t522, t525, t529, t530, t531, t532, &
         t533, t534, t535, t537, t538, t539, t54, t540, t542, t55, t550
      REAL(KIND=dp) :: t551, t552, t553, t554, t556, t557, t558, t56, t562, t563, t567, t569, &
         t571, t572, t575, t576, t577, t578, t58, t581, t582, t583, t587, t588, t589, t592, t596, &
         t597, t6, t60, t602, t603, t608, t61, t613, t617, t620, t621, t624, t626, t629, t63, &
         t631, t634, t638, t639, t644, t645, t65, t653, t655, t656, t657, t658, t659, t663, t666, &
         t669, t67, t673, t679, t68, t681, t685, t689, t69, t690, t691, t697, t698, t7, t70, t701, &
         t71, t711, t713, t717, t718, t72, t721, t728, t73, t735, t736, t739, t74, t740, t746, &
         t747, t748, t75, t752, t753, t754, t755, t759, t761, t763, t765, t769
      REAL(KIND=dp) :: t77, t771, t772, t779, t78, t780, t781, t783, t784, t787, t791, t792, t799, &
         t8, t80, t800, t801, t803, t804, t808, t81, t810, t812, t815, t816, t817, t82, t820, &
         t823, t826, t83, t832, t834, t836, t84, t842, t845, t846, t848, t849, t85, t851, t865, &
         t867, t87, t870, t871, t873, t874, t876, t877, t88, t880, t884, t885, t888, t889, t891, &
         t892, t893, t897, t898, t9, t90, t902, t904, t907, t908, t909, t91, t910, t916, t918, &
         t919, t92, t93, t930, t932, t933, t937, t94, t942, t945, t95, t951, t954, t958, t96, &
         t962, t965, t968, t97, t971, t972, t979, t982, t988, t99

      t1 = ndrho**2
      t2 = r2**2
      t3 = 0.1e1_dp/t2
      t4 = t1*t3
      t5 = pi**2
      t6 = r3*t5
      t7 = t6*rho
      t8 = t7**(0.1e1_dp/0.3e1_dp)
      t9 = t8**2
      t10 = 0.1e1_dp/t9
      t11 = t4*t10
      t12 = rho**2
      t13 = 0.1e1_dp/t12
      t14 = sscale**2
      t15 = t13*t14
      t16 = a1*t1
      t17 = t16*t3
      t18 = t10*t13
      t19 = t18*t14
      t21 = t1**2
      t22 = a2*t21
      t23 = t2**2
      t24 = 0.1e1_dp/t23
      t25 = t22*t24
      t27 = 0.1e1_dp/t8/t7
      t28 = t12**2
      t29 = 0.1e1_dp/t28
      t31 = t14**2
      t32 = t27*t29*t31
      t34 = t17*t19 + t25*t32
      t35 = a3*t21
      t36 = t35*t24
      t38 = t21*ndrho
      t39 = a4*t38
      t41 = 0.1e1_dp/t23/r2
      t42 = t39*t41
      t44 = 0.1e1_dp/t9/t7
      t46 = 0.1e1_dp/t28/rho
      t48 = t31*sscale
      t49 = t44*t46*t48
      t54 = 0.1e1_dp/t23/t2
      t55 = a5*t21*t1*t54
      t56 = r3**2
      t58 = t5**2
      t60 = 0.1e1_dp/t56/t58
      t61 = t28**2
      t63 = t31*t14
      t65 = t60/t61*t63
      t67 = r1 + t36*t32 + t42*t49 + t55*t65
      t68 = 0.1e1_dp/t67
      t69 = t34*t68
      t70 = t15*t69
      t71 = t11*t70
      t72 = omega**2
      t73 = beta*t72
      t74 = t73*t10
      t75 = t71 + t74
      t77 = 0.1e1_dp/A
      Q = f94*t75*t77
      t78 = rho**(0.1e1_dp/0.3e1_dp)
      t80 = t78*rho*f89
      t81 = B*f12
      t82 = t71 + DD
      t83 = 0.1e1_dp/t82
      t84 = t81*t83
      t85 = F2*t34
      t87 = F1 + t85*t68
      t88 = t15*t87
      t90 = t11*t88 + r1
      t91 = f12*t90
      t92 = t82**2
      t93 = 0.1e1_dp/t92
      t94 = C*t93
      t95 = t91*t94
      t96 = f34*pi
      t97 = rootpi
      t99 = r6*C
      t102 = r4*B
      t104 = r8*A
      t105 = t92*t82
      t108 = t97*(r15*E + t99*t90*t82 + t102*t92 + t104*t105)
      t109 = 0.1e1_dp/r16
      t110 = SQRT(t82)
      t111 = t110*t105
      t113 = t109/t111
      t115 = SQRT(A)
      t116 = f94*t34
      t117 = t68*t1
      t118 = t116*t117
      t119 = t3*t10
      t120 = t15*t77
      t121 = t119*t120
      t123 = EXP(t118*t121)
      t124 = t115*t123
      t125 = f32*ndrho
      t126 = 0.1e1_dp/r2
      t127 = t125*t126
      t128 = 0.1e1_dp/t8
      t129 = 0.1e1_dp/rho
      t131 = t69*t77
      t132 = SQRT(t131)
      t133 = sscale*t132
      t136 = erfc(t127*t128*t129*t133)
      t140 = 0.1e1_dp/f1516
      t141 = (t96 + t108*t113 - t96*t124*t136)*t140
      t142 = 0.1e1_dp/t97
      t144 = 0.1e1_dp/E
      t145 = t142*t111*t144
      t147 = -t141*t145 + r1
      t148 = t147*E
      t149 = 0.1e1_dp/t105
      t150 = t148*t149
      t151 = f158*E
      t152 = t147*t83
      t153 = t72*t10
      t154 = t71 + DD + t153
      t155 = t154**2
      t156 = t155**2
      t157 = t156*t154
      t158 = SQRT(t157)
      t159 = 0.1e1_dp/t158
      t162 = SQRT(t154)
      t163 = 0.1e1_dp/t162
      t166 = f68*C
      t167 = t90*t83
      t168 = t155*t154
      t169 = SQRT(t168)
      t170 = 0.1e1_dp/t169
      t174 = (-t151*t152*t159 - t81*t83*t163 - t166*t167*t170) &
             *omega
      t176 = f52*E
      t177 = t147*t93
      t180 = f12*C
      t181 = t90*t93
      t185 = t72*omega
      t186 = (-t176*t177*t159 - t180*t181*t170)*t185
      t189 = 0.1e1_dp/r3/t5
      t190 = t189*t129
      t192 = t72**2
      t193 = t192*omega
      t194 = t159*t193
      t195 = t194*t44
      t197 = f12*A
      t198 = exei(Q)
      t199 = t71 + DD + t74
      t200 = 0.1e1_dp/t199
      t202 = LOG(t75*t200)
      t205 = SQRT(t199)
      t209 = t115*f34
      t210 = exer(Q)
      t213 = (t197*t97/t205 - t209*t210)*alpha1
      t214 = omega*t128
      t219 = (t197*t200 - f98*t198)*alpha2
      t221 = A*f14
      t222 = t199**2
      t223 = t222*t199
      t224 = SQRT(t223)
      t227 = SQRT(t75)
      t228 = 0.1e1_dp/t227
      t233 = 0.1e1_dp/t115
      t237 = (t97*(t221/t224 - f98*t228) + f2716*t210*t233)*alpha3 &
             *t185
      t239 = 0.1e1_dp/t75
      t241 = 0.1e1_dp/t222
      t243 = f8132*t77
      t246 = (-f98*t239 + t197*t241 + t243*t198)*alpha4
      t247 = t192*t27
      t250 = t75**2
      t251 = t250*t75
      t252 = SQRT(t251)
      t253 = 0.1e1_dp/t252
      t255 = f38*A
      t256 = t222**2
      t257 = t256*t199
      t258 = SQRT(t257)
      t263 = A**2
      t264 = t263*A
      t265 = SQRT(t264)
      t267 = f24364/t265
      t270 = (t97*(t243*t228 - f916*t253 + t255/t258) - t267*t210) &
             *alpha5
      t271 = t193*t44
      t273 = 0.1e1_dp/t223
      t275 = 0.1e1_dp/t250
      t276 = f98*t275
      t280 = f729128/t263
      t284 = t192*t72
      t285 = (A*t273 - t276 + t243*r1*t239 - t280*t198)*alpha6 &
             *t284
      t286 = t60*t13
      t288 = f1516*A
      t289 = t256*t223
      t290 = SQRT(t289)
      t293 = t250**2
      t294 = t293*t75
      t295 = SQRT(t294)
      t298 = f8164*t77
      t303 = t263**2
      t305 = SQRT(t303*A)
      t307 = f2187256/t305
      t310 = (t97*(t288/t290 - f2732/t295 + t298*t253 - t280*t228) &
              + t307*t210)*alpha7
      t311 = t192*t185
      t312 = t56*t58
      t313 = t312*t12
      t315 = 0.1e1_dp/t8/t313
      t316 = t311*t315
      t318 = r3*A
      t319 = 0.1e1_dp/t256
      t321 = 0.1e1_dp/t251
      t326 = f6561512/t264
      t329 = (t318*t319 - f94*t321 + t243*t275 - t280*t239 + t326 &
              *t198)*alpha8
      t330 = t192**2
      t332 = 0.1e1_dp/t9/t313
      t333 = t330*t332
      t335 = t84 + t95 + t150 + t174*t128 + t186*t190 - t150*t195 + &
             t197*(t198 + t202) + t213*t214 + t219*t153 + t237*t190 + &
             t246*t247 + t270*t271 + t285*t286 + t310*t316 + t329*t333
      t336 = t335*Clda
      e_0 = e_0 + (-t80*t336)*sx
      IF (order >= 1) THEN
         t338 = t44*t13
         t339 = t4*t338
         t340 = t14*t34
         t342 = t68*r3*t5
         t343 = t340*t342
         t345 = 0.2e1_dp/0.3e1_dp*t339*t343
         t346 = t12*rho
         t347 = 0.1e1_dp/t346
         t348 = t347*t14
         t349 = t348*t69
         t351 = 2._dp*t11*t349
         t352 = t3*t44
         t353 = t16*t352
         t354 = t15*t6
         t357 = t10*t347
         t358 = t357*t14
         t361 = t24*t315
         t362 = t22*t361
         t363 = t29*t31
         t364 = t363*t6
         t368 = t27*t46*t31
         t371 = -0.2e1_dp/0.3e1_dp*t353*t354 - (2._dp*t17*t358) - 0.4e1_dp &
                /0.3e1_dp*t362*t364 - (4._dp*t25*t368)
         t372 = t371*t68
         t373 = t15*t372
         t374 = t11*t373
         t375 = t4*t18
         t376 = t67**2
         t377 = 0.1e1_dp/t376
         t378 = t35*t361
         t383 = t41*t332
         t384 = t39*t383
         t385 = t46*t48
         t386 = t385*t6
         t390 = 0.1e1_dp/t28/t12
         t392 = t44*t390*t48
         t398 = t60/t61/rho*t63
         t401 = -0.4e1_dp/0.3e1_dp*t378*t364 - (4._dp*t36*t368) - 0.5e1_dp &
                /0.3e1_dp*t384*t386 - (5._dp*t42*t392) - (8._dp*t55*t398)
         t402 = t377*t401
         t403 = t340*t402
         t404 = t375*t403
         t406 = t44*r3*t5
         t409 = -t345 - t351 + t374 - t404 - 0.2e1_dp/0.3e1_dp*t73*t406
         dQrho = f94*t409*t77
         t411 = ndrho*t3
         t412 = t411*t10
         t415 = a1*ndrho
         t416 = t415*t3
         t419 = t1*ndrho
         t420 = a2*t419
         t421 = t420*t24
         t424 = 2._dp*t416*t19 + 4._dp*t421*t32
         t425 = t424*t68
         t426 = t15*t425
         t428 = a3*t419
         t429 = t428*t24
         t432 = a4*t21
         t433 = t432*t41
         t437 = a5*t38*t54
         t440 = 4._dp*t429*t32 + 5._dp*t433*t49 + 6._dp*t437*t65
         t441 = t377*t440
         t442 = t340*t441
         t444 = 2._dp*t412*t70 + t11*t426 - t375*t442
         dQndrho = f94*t444*t77
         t446 = t78*f89
         t449 = t60*t347
         t452 = C*t149
         t453 = -t345 - t351 + t374 - t404
         t454 = t452*t453
         t457 = t329*t330
         t461 = t56*r3*t58*t5*t346
         t463 = 0.1e1_dp/t9/t461
         t465 = t463*r3*t5
         t468 = t14*t87
         t469 = t468*t6
         t472 = t348*t87
         t475 = F2*t371
         t478 = t475*t68 - t85*t402
         t479 = t15*t478
         t481 = -0.2e1_dp/0.3e1_dp*t339*t469 - (2._dp*t11*t472) + (t11 &
                                                                   *t479)
         t486 = t82*t453
         t493 = t97*(t99*t481*t82 + t99*t90*t453 + 2._dp*t102*t486 &
                     + 3._dp*t104*t92*t453)
         t495 = t92**2
         t498 = t109/t110/t495
         t499 = t498*t453
         t502 = t96*t115
         t503 = f94*t371
         t504 = t503*t117
         t507 = t377*t1*t3
         t508 = t116*t507
         t509 = t14*t77
         t510 = t509*t401
         t511 = t18*t510
         t513 = t117*t3
         t514 = t116*t513
         t515 = t338*t14
         t517 = t77*r3*t5
         t518 = t515*t517
         t522 = t119*t348*t77
         t525 = t504*t121 - t508*t511 - 0.2e1_dp/0.3e1_dp*t514*t518 - (2._dp &
                                                                       *t118*t522)
         t529 = rootpi
         t530 = 0.1e1_dp/t529
         t531 = t123*t530
         t532 = f32**2
         t533 = t532*t1
         t534 = t533*t119
         t535 = t15*t131
         t537 = EXP(-t534*t535)
         t538 = t126*t27
         t539 = t125*t538
         t540 = t129*sscale
         t542 = t132*r3*t5
         t550 = t125*t126*t128
         t551 = 0.1e1_dp/t132
         t552 = t372*t77
         t553 = t34*t377
         t554 = t77*t401
         t556 = t552 - t553*t554
         t557 = t551*t556
         t558 = t540*t557
         t562 = t537*(-t539*t540*t542/0.3e1_dp - t127*t128*t13*t133 &
                      + t550*t558/0.2e1_dp)
         t563 = t531*t562
         t567 = (t493*t113 - 0.7e1_dp/0.2e1_dp*t108*t499 - (t502*t525 &
                                                            *t123*t136) + (2._dp*t502*t563))*t140
         t569 = t141*t142
         t571 = t110*t92*t144
         t572 = t571*t453
         t575 = -t567*t145 - 0.7e1_dp/0.2e1_dp*t569*t572
         t576 = t575*E
         t577 = t576*t149
         t578 = t189*t13
         t581 = 0.1e1_dp/t158/t157
         t582 = t149*t581
         t583 = t148*t582
         t587 = -t345 - t351 + t374 - t404 - 0.2e1_dp/0.3e1_dp*t72*t44*t6
         t588 = t156*t587
         t589 = t271*t588
         t592 = t219*t72
         t596 = 0.1e1_dp/t224/t223
         t597 = t596*t222
         t602 = 0.1e1_dp/t227/t75
         t603 = f98*t602
         t608 = dexerrho(Q, dQrho)
         t613 = (t97*(-0.3e1_dp/0.2e1_dp*t221*t597*t409 + t603*t409/ &
                      0.2e1_dp) + f2716*t608*t233)*alpha3*t185
         t617 = f12*t481
         t620 = 0.1e1_dp/t495
         t621 = t620*t453
         t624 = t213*omega
         t626 = t27*r3*t5
         t629 = t246*t192
         t631 = t315*r3*t5
         t634 = t602*t409
         t638 = 0.1e1_dp/t252/t251
         t639 = f916*t638
         t644 = 0.1e1_dp/t258/t257
         t645 = t644*t256
         t653 = (t97*(-t243*t634/0.2e1_dp + 0.3e1_dp/0.2e1_dp*t639*t250* &
                      t409 - 0.5e1_dp/0.2e1_dp*t255*t645*t409) - t267*t608)*alpha5
         t655 = -(2._dp*t285*t449) - (2._dp*t91*t454) - 0.8e1_dp/0.3e1_dp &
                *t457*t465 + t577 - t186*t578 + 0.5e1_dp/0.2e1_dp*t583*t589 &
                - 0.2e1_dp/0.3e1_dp*t592*t406 + t613*t190 - t81*t93*t453 + &
                t617*t94 - t237*t578 - (3._dp*t148*t621) - t624*t626/0.3e1_dp &
                - 0.4e1_dp/0.3e1_dp*t629*t631 + t653*t271
         t656 = t149*t159
         t657 = t148*t656
         t658 = t193*t332
         t659 = t658*t6
         t663 = t273*t409
         t666 = dexeirho(Q, dQrho)
         t669 = (t276*t409 - 2._dp*t197*t663 + t243*t666)*alpha4
         t673 = t97/t205/t199
         t679 = (-t197*t673*t409/0.2e1_dp - t209*t608)*alpha1
         t681 = t241*t409
         t685 = (-t197*t681 - f98*t666)*alpha2
         t689 = 0.1e1_dp/t290/t289
         t690 = t256*t222
         t691 = t689*t690
         t697 = f2732/t295/t294
         t698 = t293*t409
         t701 = t638*t250
         t711 = (t97*(-0.7e1_dp/0.2e1_dp*t288*t691*t409 + 0.5e1_dp/0.2e1_dp &
                      *t697*t698 - 0.3e1_dp/0.2e1_dp*t298*t701*t409 + t280*t634/ &
                      0.2e1_dp) + t307*t608)*alpha7
         t713 = 0.1e1_dp/t257
         t717 = 0.1e1_dp/t293
         t718 = f94*t717
         t721 = t321*t409
         t728 = (-4._dp*t318*t713*t409 + 3._dp*t718*t409 - 2._dp*t243*t721 &
                 + t280*t275*t409 + t326*t666)*alpha8
         t735 = t176*t147
         t736 = t656*t453
         t739 = t93*t581
         t740 = t739*t588
         t746 = t180*t90
         t747 = t149*t170
         t748 = t747*t453
         t752 = 0.1e1_dp/t169/t168
         t753 = t93*t752
         t754 = t155*t587
         t755 = t753*t754
         t759 = (-t176*t575*t93*t159 + (2._dp*t735*t736) + 0.5e1_dp/ &
                 0.2e1_dp*(t735)*(t740) - t180*t481*t93*t170 + (2._dp &
                                                                *t746*t748) + 0.3e1_dp/0.2e1_dp*(t746)*(t755))*t185
         t761 = t310*t311
         t763 = 0.1e1_dp/t8/t461
         t765 = t763*r3*t5
         t769 = t75*t241
         t771 = t409*t200 - t769*t409
         t772 = t771*t239
         t779 = t151*t147
         t780 = t93*t159
         t781 = t780*t453
         t783 = t83*t581
         t784 = t783*t588
         t787 = t93*t163
         t791 = 0.1e1_dp/t162/t154
         t792 = t83*t791
         t799 = t166*t90
         t800 = t93*t170
         t801 = t800*t453
         t803 = t83*t752
         t804 = t803*t754
         t808 = (-t151*t575*t83*t159 + t779*t781 + 0.5e1_dp/0.2e1_dp*t779 &
                 *t784 + t81*t787*t453 + t81*t792*t587/0.2e1_dp - t166 &
                 *t481*t83*t170 + t799*t801 + 0.3e1_dp/0.2e1_dp*t799*t804)* &
                omega
         t810 = t148*t620
         t812 = t194*t44*t453
         t815 = t270*t193
         t816 = t332*r3
         t817 = t816*t5
         t820 = A*t319
         t823 = f98*t321
         t826 = r1*t275
         t832 = (-3._dp*t820*t409 + 2._dp*t823*t409 - t243*t826*t409 - t280 &
                 *t666)*alpha6*t284
         t834 = 0.5e1_dp/0.3e1_dp*t657*t659 + t669*t247 + t679*t214 + t685 &
                *t153 - t577*t195 + t711*t316 + t728*t333 - t174*t626 &
                /0.3e1_dp + t759*t190 - 0.7e1_dp/0.3e1_dp*t761*t765 + t197*(t666 &
                                                                  + t772*t199) + t808*t128 + (3._dp*t810*t812) - 0.5e1_dp/0.3e1_dp &
                *t815*t817 + t832*t286
         t836 = (t655 + t834)*Clda
         e_rho = e_rho + (-0.4e1_dp/0.3e1_dp*t446*t336 - t80*t836)*sx
         t842 = F2*t424
         t845 = t842*t68 - t85*t441
         t846 = t15*t845
         t848 = 2._dp*t412*t88 + t11*t846
         t849 = f12*t848
         t851 = t452*t444
         t865 = t97*(t99*t848*t82 + t99*t90*t444 + 2._dp*t102*t82 &
                     *t444 + 3._dp*t104*t92*t444)
         t867 = t498*t444
         t870 = f94*t424
         t871 = t870*t117
         t873 = t509*t440
         t874 = t18*t873
         t876 = t68*ndrho
         t877 = t116*t876
         t880 = t871*t121 - t508*t874 + 2._dp*t877*t121
         t884 = f32*t126
         t885 = t884*t128
         t888 = t425*t77
         t889 = t77*t440
         t891 = t888 - t553*t889
         t892 = t551*t891
         t893 = t540*t892
         t897 = t537*(t885*t540*t132 + t550*t893/0.2e1_dp)
         t898 = t531*t897
         t902 = (t865*t113 - 0.7e1_dp/0.2e1_dp*t108*t867 - (t502*t880 &
                                                            *t123*t136) + (2._dp*t502*t898))*t140
         t904 = t571*t444
         t907 = -t902*t145 - 0.7e1_dp/0.2e1_dp*t569*t904
         t908 = t907*E
         t909 = t908*t149
         t910 = t620*t444
         t916 = t780*t444
         t918 = t156*t444
         t919 = t783*t918
         t930 = t800*t444
         t932 = t155*t444
         t933 = t803*t932
         t937 = (-t151*t907*t83*t159 + t779*t916 + 0.5e1_dp/0.2e1_dp*t779 &
                 *t919 + t81*t787*t444 + t81*t792*t444/0.2e1_dp - t166 &
                 *t848*t83*t170 + t799*t930 + 0.3e1_dp/0.2e1_dp*t799*t933)* &
                omega
         t942 = t656*t444
         t945 = t739*t918
         t951 = t747*t444
         t954 = t753*t932
         t958 = (-t176*t907*t93*t159 + (2._dp*t735*t942) + 0.5e1_dp/ &
                 0.2e1_dp*(t735)*(t945) - t180*t848*t93*t170 + (2._dp &
                                                                *t746*t951) + 0.3e1_dp/0.2e1_dp*(t746)*(t954))*t185
         t962 = t194*t44*t444
         t965 = t271*t918
         t968 = dexeindrho(Q, dQndrho)
         t971 = t444*t200 - t769*t444
         t972 = t971*t239
         t979 = dexerndrho(Q, dQndrho)
         t982 = (-t197*t673*t444/0.2e1_dp - t209*t979)*alpha1
         t988 = (-t197*t241*t444 - f98*t968)*alpha2
         t1001 = (t97*(-0.3e1_dp/0.2e1_dp*t221*t597*t444 + t603*t444/ &
                       0.2e1_dp) + f2716*t979*t233)*alpha3*t185
         t1009 = (t276*t444 - 2._dp*t197*t273*t444 + t243*t968)*alpha4
         t1011 = t602*t444
         t1024 = (t97*(-t243*t1011/0.2e1_dp + 0.3e1_dp/0.2e1_dp*t639*t250 &
                       *t444 - 0.5e1_dp/0.2e1_dp*t255*t645*t444) - t267*t979)*alpha5
         t1035 = (-3._dp*t820*t444 + 2._dp*t823*t444 - t243*t826*t444 - &
                  t280*t968)*alpha6*t284
         t1052 = (t97*(-0.7e1_dp/0.2e1_dp*t288*t691*t444 + 0.5e1_dp/0.2e1_dp &
                       *t697*t293*t444 - 0.3e1_dp/0.2e1_dp*t298*t701*t444 + t280 &
                       *t1011/0.2e1_dp) + t307*t979)*alpha7
         t1066 = (-4._dp*t318*t713*t444 + 3._dp*t718*t444 - 2._dp*t243*t321 &
                  *t444 + t280*t275*t444 + t326*t968)*alpha8
         t1068 = -t81*t93*t444 + t849*t94 - (2._dp*t91*t851) + t909 &
                 - (3._dp*t148*t910) + t937*t128 + t958*t190 - t909*t195 &
                 + (3._dp*t810*t962) + 0.5e1_dp/0.2e1_dp*t583*t965 + t197*(t968 &
                                                                         + t972*t199) + t982*t214 + t988*t153 + t1001*t190 + t1009 &
                 *t247 + t1024*t271 + t1035*t286 + t1052*t316 + t1066* &
                 t333
         t1069 = t1068*Clda
         e_ndrho = e_ndrho + (-t80*t1069)*sx
      END IF

   END SUBROUTINE xwpbe_lda_calc_1

! **************************************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega scaling factor
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega!=0 using a taylor
!>      expansion for the parameter G.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc_2(e_0, e_rho, e_ndrho, rho, ndrho, omega, sscale, sx, order)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_0, e_rho, e_ndrho
      REAL(KIND=dp), INTENT(IN)                          :: rho, ndrho, omega, sscale, sx
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp) :: dQndrho, dQrho, Q, t1, t10, t100, t102, t103, t105, t106, t108, t109, t11, &
         t110, t111, t112, t113, t114, t115, t116, t117, t118, t12, t121, t122, t125, t126, t127, &
         t128, t129, t13, t133, t134, t136, t137, t14, t140, t141, t145, t146, t149, t15, t151, &
         t153, t154, t155, t156, t158, t159, t16, t160, t161, t163, t166, t167, t17, t171, t172, &
         t173, t176, t177, t18, t182, t184, t185, t186, t187, t19, t190, t191, t196, t2, t200, &
         t202, t204, t206, t209, t21, t210, t213, t214, t215, t216, t218, t219, t22, t220, t221, &
         t226, t227, t228, t23, t230, t233, t234, t236, t238, t239, t24
      REAL(KIND=dp) :: t243, t247, t248, t249, t25, t251, t252, t253, t256, t257, t258, t261, &
         t266, t268, t27, t270, t273, t274, t275, t276, t278, t279, t28, t281, t282, t284, t289, &
         t29, t292, t293, t295, t296, t298, t299, t3, t301, t302, t303, t305, t306, t308, t309, &
         t31, t310, t311, t312, t314, t315, t316, t317, t32, t320, t321, t324, t325, t326, t327, &
         t331, t334, t335, t336, t337, t338, t339, t34, t340, t341, t346, t347, t348, t349, t35, &
         t353, t355, t36, t361, t364, t365, t366, t367, t369, t372, t374, t375, t378, t379, t38, &
         t382, t383, t384, t387, t388, t389, t39, t391, t392, t395, t396, t4
      REAL(KIND=dp) :: t400, t403, t404, t405, t407, t409, t41, t413, t416, t417, t42, t420, t423, &
         t428, t433, t434, t436, t44, t440, t441, t442, t443, t447, t448, t452, t453, t454, t457, &
         t458, t46, t461, t464, t467, t468, t470, t474, t475, t476, t48, t480, t481, t482, t483, &
         t487, t489, t49, t491, t496, t5, t500, t503, t505, t506, t507, t510, t512, t513, t514, &
         t516, t519, t522, t523, t524, t530, t531, t535, t536, t54, t541, t542, t549, t55, t551, &
         t552, t558, t559, t56, t561, t565, t566, t569, t574, t577, t579, t58, t583, t584, t585, &
         t587, t588, t591, t595, t596, t6, t60, t603, t604, t605, t607
      REAL(KIND=dp) :: t608, t61, t612, t615, t620, t622, t626, t628, t629, t63, t633, t634, t635, &
         t638, t641, t644, t65, t650, t652, t656, t658, t659, t660, t665, t67, t670, t671, t679, &
         t68, t681, t685, t687, t689, t69, t692, t695, t697, t699, t7, t70, t702, t705, t709, t71, &
         t710, t714, t72, t723, t725, t727, t73, t733, t736, t737, t739, t74, t740, t742, t747, &
         t748, t75, t751, t752, t755, t756, t758, t759, t760, t761, t767, t769, t77, t770, t78, &
         t781, t783, t784, t788, t793, t796, t8, t80, t802, t805, t809, t81, t813, t816, t819, &
         t82, t822, t823, t83, t830, t833, t839, t84, t85, t852, t860
      REAL(KIND=dp) :: t862, t87, t875, t88, t886, t9, t90, t903, t91, t917, t919, t92, t920, t93, &
         t94, t95, t96, t97, t99

      t1 = ndrho**2
      t2 = r2**2
      t3 = 0.1e1_dp/t2
      t4 = t1*t3
      t5 = pi**2
      t6 = r3*t5
      t7 = t6*rho
      t8 = t7**(0.1e1_dp/0.3e1_dp)
      t9 = t8**2
      t10 = 0.1e1_dp/t9
      t11 = t4*t10
      t12 = rho**2
      t13 = 0.1e1_dp/t12
      t14 = sscale**2
      t15 = t13*t14
      t16 = a1*t1
      t17 = t16*t3
      t18 = t10*t13
      t19 = t18*t14
      t21 = t1**2
      t22 = a2*t21
      t23 = t2**2
      t24 = 0.1e1_dp/t23
      t25 = t22*t24
      t27 = 0.1e1_dp/t8/t7
      t28 = t12**2
      t29 = 0.1e1_dp/t28
      t31 = t14**2
      t32 = t27*t29*t31
      t34 = t17*t19 + t25*t32
      t35 = a3*t21
      t36 = t35*t24
      t38 = t21*ndrho
      t39 = a4*t38
      t41 = 0.1e1_dp/t23/r2
      t42 = t39*t41
      t44 = 0.1e1_dp/t9/t7
      t46 = 0.1e1_dp/t28/rho
      t48 = t31*sscale
      t49 = t44*t46*t48
      t54 = 0.1e1_dp/t23/t2
      t55 = a5*t21*t1*t54
      t56 = r3**2
      t58 = t5**2
      t60 = 0.1e1_dp/t56/t58
      t61 = t28**2
      t63 = t31*t14
      t65 = t60/t61*t63
      t67 = r1 + t36*t32 + t42*t49 + t55*t65
      t68 = 0.1e1_dp/t67
      t69 = t34*t68
      t70 = t15*t69
      t71 = t11*t70
      t72 = omega**2
      t73 = beta*t72
      t74 = t73*t10
      t75 = t71 + t74
      t77 = 0.1e1_dp/A
      Q = f94*t75*t77
      t78 = rho**(0.1e1_dp/0.3e1_dp)
      t80 = t78*rho*f89
      t81 = B*f12
      t82 = t71 + DD
      t83 = 0.1e1_dp/t82
      t84 = t81*t83
      t85 = F2*t34
      t87 = F1 + t85*t68
      t88 = t15*t87
      t90 = t11*t88 + r1
      t91 = f12*t90
      t92 = t82**2
      t93 = 0.1e1_dp/t92
      t94 = C*t93
      t95 = t91*t94
      t96 = g2*t1
      t97 = t96*t3
      t99 = g3*t21
      t100 = t99*t24
      t102 = g1 + t97*t19 + t100*t32
      t103 = t15*t102
      t105 = t11*t103 + r1
      t106 = t105*E
      t108 = 0.1e1_dp/t92/t82
      t109 = t106*t108
      t110 = f158*E
      t111 = t105*t83
      t112 = t72*t10
      t113 = t71 + DD + t112
      t114 = t113**2
      t115 = t114**2
      t116 = t115*t113
      t117 = SQRT(t116)
      t118 = 0.1e1_dp/t117
      t121 = SQRT(t113)
      t122 = 0.1e1_dp/t121
      t125 = f68*C
      t126 = t90*t83
      t127 = t114*t113
      t128 = SQRT(t127)
      t129 = 0.1e1_dp/t128
      t133 = (-t110*t111*t118 - t81*t83*t122 - t125*t126*t129) &
             *omega
      t134 = 0.1e1_dp/t8
      t136 = f52*E
      t137 = t105*t93
      t140 = f12*C
      t141 = t90*t93
      t145 = t72*omega
      t146 = (-t136*t137*t118 - t140*t141*t129)*t145
      t149 = 0.1e1_dp/r3/t5
      t151 = t149/rho
      t153 = t72**2
      t154 = t153*omega
      t155 = t118*t154
      t156 = t155*t44
      t158 = f12*A
      t159 = exei(Q)
      t160 = t71 + DD + t74
      t161 = 0.1e1_dp/t160
      t163 = LOG(t75*t161)
      t166 = rootpi
      t167 = SQRT(t160)
      t171 = SQRT(A)
      t172 = t171*f34
      t173 = exer(Q)
      t176 = (t158*t166/t167 - t172*t173)*alpha1
      t177 = omega*t134
      t182 = (t158*t161 - f98*t159)*alpha2
      t184 = A*f14
      t185 = t160**2
      t186 = t185*t160
      t187 = SQRT(t186)
      t190 = SQRT(t75)
      t191 = 0.1e1_dp/t190
      t196 = 0.1e1_dp/t171
      t200 = (t166*(t184/t187 - f98*t191) + f2716*t173*t196)* &
             alpha3*t145
      t202 = 0.1e1_dp/t75
      t204 = 0.1e1_dp/t185
      t206 = f8132*t77
      t209 = (-f98*t202 + t158*t204 + t206*t159)*alpha4
      t210 = t153*t27
      t213 = t75**2
      t214 = t213*t75
      t215 = SQRT(t214)
      t216 = 0.1e1_dp/t215
      t218 = f38*A
      t219 = t185**2
      t220 = t219*t160
      t221 = SQRT(t220)
      t226 = A**2
      t227 = t226*A
      t228 = SQRT(t227)
      t230 = f24364/t228
      t233 = (t166*(t206*t191 - f916*t216 + t218/t221) - t230*t173) &
             *alpha5
      t234 = t154*t44
      t236 = 0.1e1_dp/t186
      t238 = 0.1e1_dp/t213
      t239 = f98*t238
      t243 = f729128/t226
      t247 = t153*t72
      t248 = (A*t236 - t239 + t206*r1*t202 - t243*t159)*alpha6 &
             *t247
      t249 = t60*t13
      t251 = f1516*A
      t252 = t219*t186
      t253 = SQRT(t252)
      t256 = t213**2
      t257 = t256*t75
      t258 = SQRT(t257)
      t261 = f8164*t77
      t266 = t226**2
      t268 = SQRT(t266*A)
      t270 = f2187256/t268
      t273 = (t166*(t251/t253 - f2732/t258 + t261*t216 - t243*t191) &
              + t270*t173)*alpha7
      t274 = t153*t145
      t275 = t56*t58
      t276 = t275*t12
      t278 = 0.1e1_dp/t8/t276
      t279 = t274*t278
      t281 = r3*A
      t282 = 0.1e1_dp/t219
      t284 = 0.1e1_dp/t214
      t289 = f6561512/t227
      t292 = (t281*t282 - f94*t284 + t206*t238 - t243*t202 + t289 &
              *t159)*alpha8
      t293 = t153**2
      t295 = 0.1e1_dp/t9/t276
      t296 = t293*t295
      t298 = t84 + t95 + t109 + t133*t134 + t146*t151 - t109*t156 + &
             t158*(t159 + t163) + t176*t177 + t182*t112 + t200*t151 + &
             t209*t210 + t233*t234 + t248*t249 + t273*t279 + t292*t296
      t299 = t298*Clda
      e_0 = e_0 + (-t80*t299)*sx
      IF (order >= 1) THEN
         t301 = t44*t13
         t302 = t4*t301
         t303 = t14*t34
         t305 = t68*r3*t5
         t306 = t303*t305
         t308 = 0.2e1_dp/0.3e1_dp*t302*t306
         t309 = t12*rho
         t310 = 0.1e1_dp/t309
         t311 = t310*t14
         t312 = t311*t69
         t314 = 2._dp*t11*t312
         t315 = t3*t44
         t316 = t16*t315
         t317 = t15*t6
         t320 = t10*t310
         t321 = t320*t14
         t324 = t24*t278
         t325 = t22*t324
         t326 = t29*t31
         t327 = t326*t6
         t331 = t27*t46*t31
         t334 = -0.2e1_dp/0.3e1_dp*t316*t317 - (2._dp*t17*t321) - 0.4e1_dp &
                /0.3e1_dp*t325*t327 - (4._dp*t25*t331)
         t335 = t334*t68
         t336 = t15*t335
         t337 = t11*t336
         t338 = t4*t18
         t339 = t67**2
         t340 = 0.1e1_dp/t339
         t341 = t35*t324
         t346 = t41*t295
         t347 = t39*t346
         t348 = t46*t48
         t349 = t348*t6
         t353 = 0.1e1_dp/t28/t12
         t355 = t44*t353*t48
         t361 = t60/t61/rho*t63
         t364 = -0.4e1_dp/0.3e1_dp*t341*t327 - (4._dp*t36*t331) - 0.5e1_dp &
                /0.3e1_dp*t347*t349 - (5._dp*t42*t355) - (8._dp*t55*t361)
         t365 = t340*t364
         t366 = t303*t365
         t367 = t338*t366
         t369 = t44*r3*t5
         t372 = -t308 - t314 + t337 - t367 - 0.2e1_dp/0.3e1_dp*t73*t369
         dQrho = f94*t372*t77
         t374 = ndrho*t3
         t375 = t374*t10
         t378 = a1*ndrho
         t379 = t378*t3
         t382 = t1*ndrho
         t383 = a2*t382
         t384 = t383*t24
         t387 = 2._dp*t379*t19 + 4._dp*t384*t32
         t388 = t387*t68
         t389 = t15*t388
         t391 = a3*t382
         t392 = t391*t24
         t395 = a4*t21
         t396 = t395*t41
         t400 = a5*t38*t54
         t403 = 4._dp*t392*t32 + 5._dp*t396*t49 + 6._dp*t400*t65
         t404 = t340*t403
         t405 = t303*t404
         t407 = 2._dp*t375*t70 + t11*t389 - t338*t405
         dQndrho = f94*t407*t77
         t409 = t78*f89
         t413 = t27*r3*t5
         t416 = t14*t102
         t417 = t416*t6
         t420 = t311*t102
         t423 = t96*t315
         t428 = t99*t324
         t433 = -0.2e1_dp/0.3e1_dp*t423*t317 - (2._dp*t97*t321) - 0.4e1_dp &
                /0.3e1_dp*t428*t327 - (4._dp*t100*t331)
         t434 = t15*t433
         t436 = -0.2e1_dp/0.3e1_dp*t302*t417 - (2._dp*t11*t420) + (t11 &
                                                                   *t434)
         t440 = t136*t105
         t441 = t108*t118
         t442 = -t308 - t314 + t337 - t367
         t443 = t441*t442
         t447 = 0.1e1_dp/t117/t116
         t448 = t93*t447
         t452 = -t308 - t314 + t337 - t367 - 0.2e1_dp/0.3e1_dp*t72*t44*t6
         t453 = t115*t452
         t454 = t448*t453
         t457 = t14*t87
         t458 = t457*t6
         t461 = t311*t87
         t464 = F2*t334
         t467 = t464*t68 - t85*t365
         t468 = t15*t467
         t470 = -0.2e1_dp/0.3e1_dp*t302*t458 - (2._dp*t11*t461) + (t11 &
                                                                   *t468)
         t474 = t140*t90
         t475 = t108*t129
         t476 = t475*t442
         t480 = 0.1e1_dp/t128/t127
         t481 = t93*t480
         t482 = t114*t452
         t483 = t481*t482
         t487 = (-t136*t436*t93*t118 + (2._dp*t440*t443) + 0.5e1_dp/ &
                 0.2e1_dp*(t440)*(t454) - t140*t470*t93*t129 + (2._dp &
                                                                *t474*t476) + 0.3e1_dp/0.2e1_dp*(t474)*(t483))*t145
         t489 = t209*t153
         t491 = t278*r3*t5
         t496 = t166/t167/t160
         t500 = dexerrho(Q, dQrho)
         t503 = (-t158*t496*t372/0.2e1_dp - t172*t500)*alpha1
         t505 = t106*t441
         t506 = t154*t295
         t507 = t506*t6
         t510 = f12*t470
         t512 = t92**2
         t513 = 0.1e1_dp/t512
         t514 = t106*t513
         t516 = t155*t44*t442
         t519 = t149*t13
         t522 = 0.1e1_dp/t253/t252
         t523 = t219*t185
         t524 = t522*t523
         t530 = f2732/t258/t257
         t531 = t256*t372
         t535 = 0.1e1_dp/t215/t214
         t536 = t535*t213
         t541 = 0.1e1_dp/t190/t75
         t542 = t541*t372
         t549 = (t166*(-0.7e1_dp/0.2e1_dp*t251*t524*t372 + 0.5e1_dp/0.2e1_dp &
                       *t530*t531 - 0.3e1_dp/0.2e1_dp*t261*t536*t372 + t243*t542 &
                       /0.2e1_dp) + t270*t500)*alpha7
         t551 = C*t108
         t552 = t551*t442
         t558 = t436*E
         t559 = t558*t108
         t561 = 0.1e1_dp/t220
         t565 = 0.1e1_dp/t256
         t566 = f94*t565
         t569 = t284*t372
         t574 = dexeirho(Q, dQrho)
         t577 = (-4._dp*t281*t561*t372 + 3._dp*t566*t372 - 2._dp*t206*t569 &
                 + t243*t238*t372 + t289*t574)*alpha8
         t579 = -t133*t413/0.3e1_dp + t487*t151 - 0.4e1_dp/0.3e1_dp*t489* &
                t491 + t503*t177 + 0.5e1_dp/0.3e1_dp*t505*t507 + t510*t94 + (3._dp &
                                                                            *t514*t516) - t146*t519 + t549*t279 - (2._dp*t91*t552) &
                - t81*t93*t442 - t200*t519 - t559*t156 + t577*t296 &
                + t559
         t583 = t110*t105
         t584 = t93*t118
         t585 = t584*t442
         t587 = t83*t447
         t588 = t587*t453
         t591 = t93*t122
         t595 = 0.1e1_dp/t121/t113
         t596 = t83*t595
         t603 = t125*t90
         t604 = t93*t129
         t605 = t604*t442
         t607 = t83*t480
         t608 = t607*t482
         t612 = (-t110*t436*t83*t118 + t583*t585 + 0.5e1_dp/0.2e1_dp*t583 &
                 *t588 + t81*t591*t442 + t81*t596*t452/0.2e1_dp - t125 &
                 *t470*t83*t129 + t603*t605 + 0.3e1_dp/0.2e1_dp*t603*t608)* &
                omega
         t615 = t236*t372
         t620 = (t239*t372 - 2._dp*t158*t615 + t206*t574)*alpha4
         t622 = t513*t442
         t626 = t75*t204
         t628 = t372*t161 - t626*t372
         t629 = t628*t202
         t633 = t233*t154
         t634 = t295*r3
         t635 = t634*t5
         t638 = A*t282
         t641 = f98*t284
         t644 = r1*t238
         t650 = (-3._dp*t638*t372 + 2._dp*t641*t372 - t206*t644*t372 - t243 &
                 *t574)*alpha6*t247
         t652 = t204*t372
         t656 = (-t158*t652 - f98*t574)*alpha2
         t658 = t108*t447
         t659 = t106*t658
         t660 = t234*t453
         t665 = f916*t535
         t670 = 0.1e1_dp/t221/t220
         t671 = t670*t219
         t679 = (t166*(-t206*t542/0.2e1_dp + 0.3e1_dp/0.2e1_dp*t665*t213 &
                       *t372 - 0.5e1_dp/0.2e1_dp*t218*t671*t372) - t230*t500)*alpha5
         t681 = t292*t293
         t685 = t56*r3*t58*t5*t309
         t687 = 0.1e1_dp/t9/t685
         t689 = t687*r3*t5
         t692 = t60*t310
         t695 = t273*t274
         t697 = 0.1e1_dp/t8/t685
         t699 = t697*r3*t5
         t702 = t176*omega
         t705 = t182*t72
         t709 = 0.1e1_dp/t187/t186
         t710 = t709*t185
         t714 = f98*t541
         t723 = (t166*(-0.3e1_dp/0.2e1_dp*t184*t710*t372 + t714*t372/ &
                       0.2e1_dp) + f2716*t500*t196)*alpha3*t145
         t725 = t612*t134 + t620*t210 - (3._dp*t106*t622) + t158*(t574 &
                                                                  + t629*t160) - 0.5e1_dp/0.3e1_dp*t633*t635 + t650*t249 + &
                t656*t112 + 0.5e1_dp/0.2e1_dp*t659*t660 + t679*t234 - 0.8e1_dp/ &
                0.3e1_dp*t681*t689 - (2._dp*t248*t692) - 0.7e1_dp/0.3e1_dp*t695 &
                *t699 - t702*t413/0.3e1_dp - 0.2e1_dp/0.3e1_dp*t705*t369 + t723 &
                *t151
         t727 = (t579 + t725)*Clda
         e_rho = e_rho + (-0.4e1_dp/0.3e1_dp*t409*t299 - t80*t727)*sx
         t733 = F2*t387
         t736 = t733*t68 - t85*t404
         t737 = t15*t736
         t739 = 2._dp*t375*t88 + t11*t737
         t740 = f12*t739
         t742 = t551*t407
         t747 = g2*ndrho
         t748 = t747*t3
         t751 = g3*t382
         t752 = t751*t24
         t755 = 2._dp*t748*t19 + 4._dp*t752*t32
         t756 = t15*t755
         t758 = 2._dp*t375*t103 + t11*t756
         t759 = t758*E
         t760 = t759*t108
         t761 = t513*t407
         t767 = t584*t407
         t769 = t115*t407
         t770 = t587*t769
         t781 = t604*t407
         t783 = t114*t407
         t784 = t607*t783
         t788 = (-t110*t758*t83*t118 + t583*t767 + 0.5e1_dp/0.2e1_dp*t583 &
                 *t770 + t81*t591*t407 + t81*t596*t407/0.2e1_dp - t125 &
                 *t739*t83*t129 + t603*t781 + 0.3e1_dp/0.2e1_dp*t603*t784)* &
                omega
         t793 = t441*t407
         t796 = t448*t769
         t802 = t475*t407
         t805 = t481*t783
         t809 = (-t136*t758*t93*t118 + (2._dp*t440*t793) + 0.5e1_dp/ &
                 0.2e1_dp*(t440)*(t796) - t140*t739*t93*t129 + (2._dp &
                                                                *t474*t802) + 0.3e1_dp/0.2e1_dp*(t474)*(t805))*t145
         t813 = t155*t44*t407
         t816 = t234*t769
         t819 = dexeindrho(Q, dQndrho)
         t822 = t407*t161 - t626*t407
         t823 = t822*t202
         t830 = dexerndrho(Q, dQndrho)
         t833 = (-t158*t496*t407/0.2e1_dp - t172*t830)*alpha1
         t839 = (-t158*t204*t407 - f98*t819)*alpha2
         t852 = (t166*(-0.3e1_dp/0.2e1_dp*t184*t710*t407 + t714*t407/ &
                       0.2e1_dp) + f2716*t830*t196)*alpha3*t145
         t860 = (t239*t407 - 2._dp*t158*t236*t407 + t206*t819)*alpha4
         t862 = t541*t407
         t875 = (t166*(-t206*t862/0.2e1_dp + 0.3e1_dp/0.2e1_dp*t665*t213 &
                       *t407 - 0.5e1_dp/0.2e1_dp*t218*t671*t407) - t230*t830)*alpha5
         t886 = (-3._dp*t638*t407 + 2._dp*t641*t407 - t206*t644*t407 - t243 &
                 *t819)*alpha6*t247
         t903 = (t166*(-0.7e1_dp/0.2e1_dp*t251*t524*t407 + 0.5e1_dp/0.2e1_dp &
                       *t530*t256*t407 - 0.3e1_dp/0.2e1_dp*t261*t536*t407 + t243 &
                       *t862/0.2e1_dp) + t270*t830)*alpha7
         t917 = (-4._dp*t281*t561*t407 + 3._dp*t566*t407 - 2._dp*t206*t284 &
                 *t407 + t243*t238*t407 + t289*t819)*alpha8
         t919 = -t81*t93*t407 + t740*t94 - (2._dp*t91*t742) + t760 &
                - (3._dp*t106*t761) + t788*t134 + t809*t151 - t760*t156 &
                + (3._dp*t514*t813) + 0.5e1_dp/0.2e1_dp*t659*t816 + t158*(t819 &
                                                                          + t823*t160) + t833*t177 + t839*t112 + t852*t151 + t860 &
                *t210 + t875*t234 + t886*t249 + t903*t279 + t917*t296
         t920 = t919*Clda
         e_ndrho = e_ndrho + (-t80*t920)*sx
      END IF

   END SUBROUTINE xwpbe_lda_calc_2

! **************************************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega screening parameter
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega!=0 using a simple
!>      gaussian expansion for large ww.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc_3(e_0, e_rho, e_ndrho, rho, ndrho, omega, sscale, sx, order)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_0, e_rho, e_ndrho
      REAL(KIND=dp), INTENT(IN)                          :: rho, ndrho, omega, sscale, sx
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp) :: dQndrho, dQrho, Q, t1, t10, t102, t104, t105, t108, t109, t11, t110, t111, &
         t113, t115, t116, t117, t118, t119, t12, t120, t121, t123, t124, t125, t126, t127, t128, &
         t129, t13, t131, t132, t133, t136, t14, t140, t141, t142, t144, t145, t147, t148, t149, &
         t15, t150, t151, t152, t154, t155, t156, t157, t158, t159, t16, t162, t163, t166, t167, &
         t168, t169, t17, t170, t174, t176, t177, t18, t180, t181, t185, t186, t189, t19, t190, &
         t192, t193, t194, t195, t197, t198, t199, t2, t200, t202, t206, t208, t209, t21, t210, &
         t212, t213, t215, t216, t217, t218, t219, t22, t221, t222, t223
      REAL(KIND=dp) :: t224, t227, t228, t23, t231, t232, t234, t235, t236, t237, t238, t24, t242, &
         t245, t246, t247, t248, t249, t25, t250, t251, t252, t258, t259, t260, t261, t262, t266, &
         t268, t27, t274, t277, t278, t279, t28, t280, t285, t287, t288, t29, t291, t292, t295, &
         t296, t297, t3, t300, t301, t302, t304, t305, t308, t309, t31, t313, t316, t317, t318, &
         t32, t320, t322, t325, t328, t329, t332, t335, t338, t339, t34, t341, t342, t344, t345, &
         t35, t352, t359, t36, t361, t364, t365, t368, t369, t370, t373, t374, t375, t376, t377, &
         t379, t38, t380, t381, t383, t384, t388, t39, t391, t395, t396
      REAL(KIND=dp) :: t397, t398, t399, t4, t400, t401, t403, t404, t405, t406, t408, t41, t416, &
         t417, t418, t419, t42, t420, t422, t423, t424, t428, t429, t433, t435, t437, t438, t44, &
         t441, t442, t443, t444, t445, t451, t452, t453, t456, t457, t46, t461, t462, t463, t466, &
         t470, t471, t478, t479, t48, t480, t483, t484, t485, t486, t49, t490, t493, t499, t5, &
         t500, t501, t504, t505, t511, t512, t513, t516, t517, t521, t523, t526, t528, t531, t532, &
         t533, t534, t537, t538, t539, t54, t542, t544, t545, t546, t548, t549, t55, t550, t554, &
         t555, t56, t561, t564, t565, t567, t568, t570, t58, t584, t586
      REAL(KIND=dp) :: t589, t590, t592, t593, t595, t596, t599, t6, t60, t603, t604, t607, t608, &
         t61, t610, t611, t612, t616, t617, t621, t623, t626, t627, t628, t629, t63, t635, t637, &
         t638, t649, t65, t651, t652, t656, t661, t664, t67, t670, t673, t677, t68, t681, t684, &
         t687, t69, t690, t691, t695, t696, t7, t70, t71, t72, t73, t74, t75, t77, t78, t8, t80, &
         t81, t82, t83, t84, t85, t87, t88, t9, t90, t91, t92, t93, t94, t95, t96, t97, t99

      t1 = ndrho**2
      t2 = r2**2
      t3 = 0.1e1_dp/t2
      t4 = t1*t3
      t5 = pi**2
      t6 = r3*t5
      t7 = t6*rho
      t8 = t7**(0.1e1_dp/0.3e1_dp)
      t9 = t8**2
      t10 = 0.1e1_dp/t9
      t11 = t4*t10
      t12 = rho**2
      t13 = 0.1e1_dp/t12
      t14 = sscale**2
      t15 = t13*t14
      t16 = a1*t1
      t17 = t16*t3
      t18 = t10*t13
      t19 = t18*t14
      t21 = t1**2
      t22 = a2*t21
      t23 = t2**2
      t24 = 0.1e1_dp/t23
      t25 = t22*t24
      t27 = 0.1e1_dp/t8/t7
      t28 = t12**2
      t29 = 0.1e1_dp/t28
      t31 = t14**2
      t32 = t27*t29*t31
      t34 = t17*t19 + t25*t32
      t35 = a3*t21
      t36 = t35*t24
      t38 = t21*ndrho
      t39 = a4*t38
      t41 = 0.1e1_dp/t23/r2
      t42 = t39*t41
      t44 = 0.1e1_dp/t9/t7
      t46 = 0.1e1_dp/t28/rho
      t48 = t31*sscale
      t49 = t44*t46*t48
      t54 = 0.1e1_dp/t23/t2
      t55 = a5*t21*t1*t54
      t56 = r3**2
      t58 = t5**2
      t60 = 0.1e1_dp/t56/t58
      t61 = t28**2
      t63 = t31*t14
      t65 = t60/t61*t63
      t67 = r1 + t36*t32 + t42*t49 + t55*t65
      t68 = 0.1e1_dp/t67
      t69 = t34*t68
      t70 = t15*t69
      t71 = t11*t70
      t72 = omega**2
      t73 = beta2*t72
      t74 = t73*t10
      t75 = t71 + t74
      t77 = 0.1e1_dp/A
      Q = f94*t75*t77
      t78 = rho**(0.1e1_dp/0.3e1_dp)
      t80 = t78*rho*f89
      t81 = B*f12
      t82 = t71 + DD
      t83 = 0.1e1_dp/t82
      t84 = t81*t83
      t85 = F2*t34
      t87 = F1 + t85*t68
      t88 = t15*t87
      t90 = t11*t88 + r1
      t91 = f12*t90
      t92 = t82**2
      t93 = 0.1e1_dp/t92
      t94 = C*t93
      t95 = t91*t94
      t96 = f34*pi
      t97 = rootpi
      t99 = r6*C
      t102 = r4*B
      t104 = r8*A
      t105 = t92*t82
      t108 = t97*(r15*E + t99*t90*t82 + t102*t92 + t104*t105)
      t109 = 0.1e1_dp/r16
      t110 = SQRT(t82)
      t111 = t110*t105
      t113 = t109/t111
      t115 = SQRT(A)
      t116 = f94*t34
      t117 = t68*t1
      t118 = t116*t117
      t119 = t3*t10
      t120 = t15*t77
      t121 = t119*t120
      t123 = EXP(t118*t121)
      t124 = t115*t123
      t125 = f32*ndrho
      t126 = 0.1e1_dp/r2
      t127 = t125*t126
      t128 = 0.1e1_dp/t8
      t129 = 0.1e1_dp/rho
      t131 = t69*t77
      t132 = SQRT(t131)
      t133 = sscale*t132
      t136 = erfc(t127*t128*t129*t133)
      t140 = 0.1e1_dp/f1516
      t141 = (t96 + t108*t113 - t96*t124*t136)*t140
      t142 = 0.1e1_dp/t97
      t144 = 0.1e1_dp/E
      t145 = t142*t111*t144
      t147 = -t141*t145 + r1
      t148 = t147*E
      t149 = 0.1e1_dp/t105
      t150 = t148*t149
      t151 = f158*E
      t152 = t147*t83
      t154 = t71 + DD + t72*t10
      t155 = t154**2
      t156 = t155**2
      t157 = t156*t154
      t158 = SQRT(t157)
      t159 = 0.1e1_dp/t158
      t162 = SQRT(t154)
      t163 = 0.1e1_dp/t162
      t166 = f68*C
      t167 = t90*t83
      t168 = t155*t154
      t169 = SQRT(t168)
      t170 = 0.1e1_dp/t169
      t174 = (-t151*t152*t159 - t81*t83*t163 - t166*t167*t170) &
             *omega
      t176 = f52*E
      t177 = t147*t93
      t180 = f12*C
      t181 = t90*t93
      t185 = t72*omega
      t186 = (-t176*t177*t159 - t180*t181*t170)*t185
      t189 = 0.1e1_dp/r3/t5
      t190 = t189*t129
      t192 = t72**2
      t193 = t192*omega
      t194 = t159*t193
      t195 = t194*t44
      t197 = f12*A
      t198 = exei(Q)
      t199 = t71 + DD + t74
      t200 = 0.1e1_dp/t199
      t202 = LOG(t75*t200)
      t206 = (t84 + t95 + t150 + t174*t128 + t186*t190 - t150*t195 &
              + t197*(t198 + t202))*Clda
      e_0 = e_0 + (-t80*t206)*sx
      IF (order >= 1) THEN
         t208 = t44*t13
         t209 = t4*t208
         t210 = t14*t34
         t212 = t68*r3*t5
         t213 = t210*t212
         t215 = 0.2e1_dp/0.3e1_dp*t209*t213
         t216 = t12*rho
         t217 = 0.1e1_dp/t216
         t218 = t217*t14
         t219 = t218*t69
         t221 = 2._dp*t11*t219
         t222 = t3*t44
         t223 = t16*t222
         t224 = t15*t6
         t227 = t10*t217
         t228 = t227*t14
         t231 = t56*t58
         t232 = t231*t12
         t234 = 0.1e1_dp/t8/t232
         t235 = t24*t234
         t236 = t22*t235
         t237 = t29*t31
         t238 = t237*t6
         t242 = t27*t46*t31
         t245 = -0.2e1_dp/0.3e1_dp*t223*t224 - (2._dp*t17*t228) - 0.4e1_dp &
                /0.3e1_dp*t236*t238 - (4._dp*t25*t242)
         t246 = t245*t68
         t247 = t15*t246
         t248 = t11*t247
         t249 = t4*t18
         t250 = t67**2
         t251 = 0.1e1_dp/t250
         t252 = t35*t235
         t258 = 0.1e1_dp/t9/t232
         t259 = t41*t258
         t260 = t39*t259
         t261 = t46*t48
         t262 = t261*t6
         t266 = 0.1e1_dp/t28/t12
         t268 = t44*t266*t48
         t274 = t60/t61/rho*t63
         t277 = -0.4e1_dp/0.3e1_dp*t252*t238 - (4._dp*t36*t242) - 0.5e1_dp &
                /0.3e1_dp*t260*t262 - (5._dp*t42*t268) - (8._dp*t55*t274)
         t278 = t251*t277
         t279 = t210*t278
         t280 = t249*t279
         t285 = -t215 - t221 + t248 - t280 - 0.2e1_dp/0.3e1_dp*t73*t44*r3 &
                *t5
         dQrho = f94*t285*t77
         t287 = ndrho*t3
         t288 = t287*t10
         t291 = a1*ndrho
         t292 = t291*t3
         t295 = t1*ndrho
         t296 = a2*t295
         t297 = t296*t24
         t300 = 2._dp*t292*t19 + 4._dp*t297*t32
         t301 = t300*t68
         t302 = t15*t301
         t304 = a3*t295
         t305 = t304*t24
         t308 = a4*t21
         t309 = t308*t41
         t313 = a5*t38*t54
         t316 = 4._dp*t305*t32 + 5._dp*t309*t49 + 6._dp*t313*t65
         t317 = t251*t316
         t318 = t210*t317
         t320 = 2._dp*t288*t70 + t11*t302 - t249*t318
         dQndrho = f94*t320*t77
         t322 = t78*f89
         t325 = -t215 - t221 + t248 - t280
         t328 = t14*t87
         t329 = t328*t6
         t332 = t218*t87
         t335 = F2*t245
         t338 = t335*t68 - t85*t278
         t339 = t15*t338
         t341 = -0.2e1_dp/0.3e1_dp*t209*t329 - (2._dp*t11*t332) + (t11 &
                                                                   *t339)
         t342 = f12*t341
         t344 = C*t149
         t345 = t344*t325
         t352 = t82*t325
         t359 = t97*(t99*t341*t82 + t99*t90*t325 + 2._dp*t102*t352 &
                     + 3._dp*t104*t92*t325)
         t361 = t92**2
         t364 = t109/t110/t361
         t365 = t364*t325
         t368 = t96*t115
         t369 = f94*t245
         t370 = t369*t117
         t373 = t251*t1*t3
         t374 = t116*t373
         t375 = t14*t77
         t376 = t375*t277
         t377 = t18*t376
         t379 = t117*t3
         t380 = t116*t379
         t381 = t208*t14
         t383 = t77*r3*t5
         t384 = t381*t383
         t388 = t119*t218*t77
         t391 = t370*t121 - t374*t377 - 0.2e1_dp/0.3e1_dp*t380*t384 - (2._dp &
                                                                       *t118*t388)
         t395 = rootpi
         t396 = 0.1e1_dp/t395
         t397 = t123*t396
         t398 = f32**2
         t399 = t398*t1
         t400 = t399*t119
         t401 = t15*t131
         t403 = EXP(-t400*t401)
         t404 = t126*t27
         t405 = t125*t404
         t406 = t129*sscale
         t408 = t132*r3*t5
         t416 = t125*t126*t128
         t417 = 0.1e1_dp/t132
         t418 = t246*t77
         t419 = t34*t251
         t420 = t77*t277
         t422 = t418 - t419*t420
         t423 = t417*t422
         t424 = t406*t423
         t428 = t403*(-t405*t406*t408/0.3e1_dp - t127*t128*t13*t133 &
                      + t416*t424/0.2e1_dp)
         t429 = t397*t428
         t433 = (t359*t113 - 0.7e1_dp/0.2e1_dp*t108*t365 - (t368*t391 &
                                                            *t123*t136) + (2._dp*t368*t429))*t140
         t435 = t141*t142
         t437 = t110*t92*t144
         t438 = t437*t325
         t441 = -t433*t145 - 0.7e1_dp/0.2e1_dp*t435*t438
         t442 = t441*E
         t443 = t442*t149
         t444 = 0.1e1_dp/t361
         t445 = t444*t325
         t451 = t151*t147
         t452 = t93*t159
         t453 = t452*t325
         t456 = 0.1e1_dp/t158/t157
         t457 = t83*t456
         t461 = -t215 - t221 + t248 - t280 - 0.2e1_dp/0.3e1_dp*t72*t44*t6
         t462 = t156*t461
         t463 = t457*t462
         t466 = t93*t163
         t470 = 0.1e1_dp/t162/t154
         t471 = t83*t470
         t478 = t166*t90
         t479 = t93*t170
         t480 = t479*t325
         t483 = 0.1e1_dp/t169/t168
         t484 = t83*t483
         t485 = t155*t461
         t486 = t484*t485
         t490 = (-t151*t441*t83*t159 + t451*t453 + 0.5e1_dp/0.2e1_dp*t451 &
                 *t463 + t81*t466*t325 + t81*t471*t461/0.2e1_dp - t166 &
                 *t341*t83*t170 + t478*t480 + 0.3e1_dp/0.2e1_dp*t478*t486)* &
                omega
         t493 = t27*r3*t5
         t499 = t176*t147
         t500 = t149*t159
         t501 = t500*t325
         t504 = t93*t456
         t505 = t504*t462
         t511 = t180*t90
         t512 = t149*t170
         t513 = t512*t325
         t516 = t93*t483
         t517 = t516*t485
         t521 = (-t176*t441*t93*t159 + (2._dp*t499*t501) + 0.5e1_dp/ &
                 0.2e1_dp*(t499)*(t505) - t180*t341*t93*t170 + (2._dp &
                                                                *t511*t513) + 0.3e1_dp/0.2e1_dp*(t511)*(t517))*t185
         t523 = t189*t13
         t526 = t148*t444
         t528 = t194*t44*t325
         t531 = t149*t456
         t532 = t148*t531
         t533 = t193*t44
         t534 = t533*t462
         t537 = t148*t500
         t538 = t193*t258
         t539 = t538*t6
         t542 = dexeirho(Q, dQrho)
         t544 = t199**2
         t545 = 0.1e1_dp/t544
         t546 = t75*t545
         t548 = t285*t200 - t546*t285
         t549 = 0.1e1_dp/t75
         t550 = t548*t549
         t554 = -t81*t93*t325 + t342*t94 - (2._dp*t91*t345) + t443 &
                - (3._dp*t148*t445) + t490*t128 - t174*t493/0.3e1_dp + t521 &
                *t190 - t186*t523 - t443*t195 + (3._dp*t526*t528) + 0.5e1_dp &
                /0.2e1_dp*t532*t534 + 0.5e1_dp/0.3e1_dp*t537*t539 + t197*(t542 &
                                                                          + t550*t199)
         t555 = t554*Clda
         e_rho = e_rho + (-0.4e1_dp/0.3e1_dp*t322*t206 - t80*t555)*sx
         t561 = F2*t300
         t564 = t561*t68 - t85*t317
         t565 = t15*t564
         t567 = 2._dp*t288*t88 + t11*t565
         t568 = f12*t567
         t570 = t344*t320
         t584 = t97*(t99*t567*t82 + t99*t90*t320 + 2._dp*t102*t82 &
                     *t320 + 3._dp*t104*t92*t320)
         t586 = t364*t320
         t589 = f94*t300
         t590 = t589*t117
         t592 = t375*t316
         t593 = t18*t592
         t595 = t68*ndrho
         t596 = t116*t595
         t599 = t590*t121 - t374*t593 + 2._dp*t596*t121
         t603 = f32*t126
         t604 = t603*t128
         t607 = t301*t77
         t608 = t77*t316
         t610 = t607 - t419*t608
         t611 = t417*t610
         t612 = t406*t611
         t616 = t403*(t604*t406*t132 + t416*t612/0.2e1_dp)
         t617 = t397*t616
         t621 = (t584*t113 - 0.7e1_dp/0.2e1_dp*t108*t586 - (t368*t599 &
                                                            *t123*t136) + (2._dp*t368*t617))*t140
         t623 = t437*t320
         t626 = -t621*t145 - 0.7e1_dp/0.2e1_dp*t435*t623
         t627 = t626*E
         t628 = t627*t149
         t629 = t444*t320
         t635 = t452*t320
         t637 = t156*t320
         t638 = t457*t637
         t649 = t479*t320
         t651 = t155*t320
         t652 = t484*t651
         t656 = (-t151*t626*t83*t159 + t451*t635 + 0.5e1_dp/0.2e1_dp*t451 &
                 *t638 + t81*t466*t320 + t81*t471*t320/0.2e1_dp - t166 &
                 *t567*t83*t170 + t478*t649 + 0.3e1_dp/0.2e1_dp*t478*t652)* &
                omega
         t661 = t500*t320
         t664 = t504*t637
         t670 = t512*t320
         t673 = t516*t651
         t677 = (-t176*t626*t93*t159 + (2._dp*t499*t661) + 0.5e1_dp/ &
                 0.2e1_dp*(t499)*(t664) - t180*t567*t93*t170 + (2._dp &
                                                                *t511*t670) + 0.3e1_dp/0.2e1_dp*(t511)*(t673))*t185
         t681 = t194*t44*t320
         t684 = t533*t637
         t687 = dexeindrho(Q, dQndrho)
         t690 = t320*t200 - t546*t320
         t691 = t690*t549
         t695 = -t81*t93*t320 + t568*t94 - (2._dp*t91*t570) + t628 &
                - (3._dp*t148*t629) + t656*t128 + t677*t190 - t628*t195 &
                + (3._dp*t526*t681) + 0.5e1_dp/0.2e1_dp*t532*t684 + t197*(t687 &
                                                                          + t691*t199)
         t696 = t695*Clda
         e_ndrho = e_ndrho + (-t80*t696)*sx
      END IF

   END SUBROUTINE xwpbe_lda_calc_3

! **************************************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega screening parameter
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega!=0 using a simple
!>      gaussian expansion for large ww and a taylor expansion for the
!>      parameter G.
! **************************************************************************************************
   SUBROUTINE xwpbe_lda_calc_4(e_0, e_rho, e_ndrho, rho, ndrho, omega, sscale, sx, order)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_0, e_rho, e_ndrho
      REAL(KIND=dp), INTENT(IN)                          :: rho, ndrho, omega, sscale, sx
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp) :: dQndrho, dQrho, Q, t1, t10, t100, t102, t103, t105, t106, t108, t109, t11, &
         t110, t111, t113, t114, t115, t116, t117, t118, t12, t121, t122, t125, t126, t127, t128, &
         t129, t13, t133, t134, t136, t137, t14, t140, t141, t145, t146, t149, t15, t151, t153, &
         t154, t155, t156, t158, t159, t16, t160, t161, t163, t167, t169, t17, t170, t171, t173, &
         t174, t176, t177, t178, t179, t18, t180, t182, t183, t184, t185, t188, t189, t19, t192, &
         t193, t195, t196, t197, t198, t199, t2, t203, t206, t207, t208, t209, t21, t210, t211, &
         t212, t213, t219, t22, t220, t221, t222, t223, t227, t229, t23
      REAL(KIND=dp) :: t235, t238, t239, t24, t240, t241, t246, t248, t249, t25, t252, t253, t256, &
         t257, t258, t261, t262, t263, t265, t266, t269, t27, t270, t274, t277, t278, t279, t28, &
         t281, t283, t286, t289, t29, t290, t293, t296, t299, t3, t300, t302, t303, t305, t306, &
         t309, t31, t310, t313, t316, t32, t321, t326, t327, t329, t330, t331, t332, t333, t334, &
         t34, t340, t341, t342, t345, t346, t35, t350, t351, t352, t355, t359, t36, t360, t367, &
         t368, t369, t372, t373, t374, t375, t379, t38, t382, t388, t389, t39, t390, t393, t394, &
         t4, t400, t401, t402, t405, t406, t41, t410, t412, t415, t417
      REAL(KIND=dp) :: t42, t420, t421, t422, t423, t426, t427, t428, t431, t433, t434, t435, &
         t437, t438, t439, t44, t443, t444, t450, t453, t454, t456, t457, t459, t46, t464, t465, &
         t468, t469, t472, t473, t475, t476, t477, t478, t48, t484, t486, t487, t49, t498, t5, &
         t500, t501, t505, t510, t513, t519, t522, t526, t530, t533, t536, t539, t54, t540, t544, &
         t545, t55, t56, t58, t6, t60, t61, t63, t65, t67, t68, t69, t7, t70, t71, t72, t73, t74, &
         t75, t77, t78, t8, t80, t81, t82, t83, t84, t85, t87, t88, t9, t90, t91, t92, t93, t94, &
         t95, t96, t97, t99

      t1 = ndrho**2
      t2 = r2**2
      t3 = 0.1e1_dp/t2
      t4 = t1*t3
      t5 = pi**2
      t6 = r3*t5
      t7 = t6*rho
      t8 = t7**(0.1e1_dp/0.3e1_dp)
      t9 = t8**2
      t10 = 0.1e1_dp/t9
      t11 = t4*t10
      t12 = rho**2
      t13 = 0.1e1_dp/t12
      t14 = sscale**2
      t15 = t13*t14
      t16 = a1*t1
      t17 = t16*t3
      t18 = t10*t13
      t19 = t18*t14
      t21 = t1**2
      t22 = a2*t21
      t23 = t2**2
      t24 = 0.1e1_dp/t23
      t25 = t22*t24
      t27 = 0.1e1_dp/t8/t7
      t28 = t12**2
      t29 = 0.1e1_dp/t28
      t31 = t14**2
      t32 = t27*t29*t31
      t34 = t17*t19 + t25*t32
      t35 = a3*t21
      t36 = t35*t24
      t38 = t21*ndrho
      t39 = a4*t38
      t41 = 0.1e1_dp/t23/r2
      t42 = t39*t41
      t44 = 0.1e1_dp/t9/t7
      t46 = 0.1e1_dp/t28/rho
      t48 = t31*sscale
      t49 = t44*t46*t48
      t54 = 0.1e1_dp/t23/t2
      t55 = a5*t21*t1*t54
      t56 = r3**2
      t58 = t5**2
      t60 = 0.1e1_dp/t56/t58
      t61 = t28**2
      t63 = t31*t14
      t65 = t60/t61*t63
      t67 = r1 + t36*t32 + t42*t49 + t55*t65
      t68 = 0.1e1_dp/t67
      t69 = t34*t68
      t70 = t15*t69
      t71 = t11*t70
      t72 = omega**2
      t73 = beta2*t72
      t74 = t73*t10
      t75 = t71 + t74
      t77 = 0.1e1_dp/A
      Q = f94*t75*t77
      t78 = rho**(0.1e1_dp/0.3e1_dp)
      t80 = t78*rho*f89
      t81 = B*f12
      t82 = t71 + DD
      t83 = 0.1e1_dp/t82
      t84 = t81*t83
      t85 = F2*t34
      t87 = F1 + t85*t68
      t88 = t15*t87
      t90 = t11*t88 + r1
      t91 = f12*t90
      t92 = t82**2
      t93 = 0.1e1_dp/t92
      t94 = C*t93
      t95 = t91*t94
      t96 = g2*t1
      t97 = t96*t3
      t99 = g3*t21
      t100 = t99*t24
      t102 = g1 + t97*t19 + t100*t32
      t103 = t15*t102
      t105 = t11*t103 + r1
      t106 = t105*E
      t108 = 0.1e1_dp/t92/t82
      t109 = t106*t108
      t110 = f158*E
      t111 = t105*t83
      t113 = t71 + DD + t72*t10
      t114 = t113**2
      t115 = t114**2
      t116 = t115*t113
      t117 = SQRT(t116)
      t118 = 0.1e1_dp/t117
      t121 = SQRT(t113)
      t122 = 0.1e1_dp/t121
      t125 = f68*C
      t126 = t90*t83
      t127 = t114*t113
      t128 = SQRT(t127)
      t129 = 0.1e1_dp/t128
      t133 = (-t110*t111*t118 - t81*t83*t122 - t125*t126*t129) &
             *omega
      t134 = 0.1e1_dp/t8
      t136 = f52*E
      t137 = t105*t93
      t140 = f12*C
      t141 = t90*t93
      t145 = t72*omega
      t146 = (-t136*t137*t118 - t140*t141*t129)*t145
      t149 = 0.1e1_dp/r3/t5
      t151 = t149/rho
      t153 = t72**2
      t154 = t153*omega
      t155 = t118*t154
      t156 = t155*t44
      t158 = f12*A
      t159 = exei(Q)
      t160 = t71 + DD + t74
      t161 = 0.1e1_dp/t160
      t163 = LOG(t75*t161)
      t167 = (t84 + t95 + t109 + t133*t134 + t146*t151 - t109*t156 &
              + t158*(t159 + t163))*Clda
      e_0 = e_0 + (-t80*t167)*sx
      IF (order >= 1) THEN
         t169 = t44*t13
         t170 = t4*t169
         t171 = t14*t34
         t173 = t68*r3*t5
         t174 = t171*t173
         t176 = 0.2e1_dp/0.3e1_dp*t170*t174
         t177 = t12*rho
         t178 = 0.1e1_dp/t177
         t179 = t178*t14
         t180 = t179*t69
         t182 = 2._dp*t11*t180
         t183 = t3*t44
         t184 = t16*t183
         t185 = t15*t6
         t188 = t10*t178
         t189 = t188*t14
         t192 = t56*t58
         t193 = t192*t12
         t195 = 0.1e1_dp/t8/t193
         t196 = t24*t195
         t197 = t22*t196
         t198 = t29*t31
         t199 = t198*t6
         t203 = t27*t46*t31
         t206 = -0.2e1_dp/0.3e1_dp*t184*t185 - (2._dp*t17*t189) - 0.4e1_dp &
                /0.3e1_dp*t197*t199 - (4._dp*t25*t203)
         t207 = t206*t68
         t208 = t15*t207
         t209 = t11*t208
         t210 = t4*t18
         t211 = t67**2
         t212 = 0.1e1_dp/t211
         t213 = t35*t196
         t219 = 0.1e1_dp/t9/t193
         t220 = t41*t219
         t221 = t39*t220
         t222 = t46*t48
         t223 = t222*t6
         t227 = 0.1e1_dp/t28/t12
         t229 = t44*t227*t48
         t235 = t60/t61/rho*t63
         t238 = -0.4e1_dp/0.3e1_dp*t213*t199 - (4._dp*t36*t203) - 0.5e1_dp &
                /0.3e1_dp*t221*t223 - (5._dp*t42*t229) - (8._dp*t55*t235)
         t239 = t212*t238
         t240 = t171*t239
         t241 = t210*t240
         t246 = -t176 - t182 + t209 - t241 - 0.2e1_dp/0.3e1_dp*t73*t44*r3 &
                *t5
         dQrho = f94*t246*t77
         t248 = ndrho*t3
         t249 = t248*t10
         t252 = a1*ndrho
         t253 = t252*t3
         t256 = t1*ndrho
         t257 = a2*t256
         t258 = t257*t24
         t261 = 2._dp*t253*t19 + 4._dp*t258*t32
         t262 = t261*t68
         t263 = t15*t262
         t265 = a3*t256
         t266 = t265*t24
         t269 = a4*t21
         t270 = t269*t41
         t274 = a5*t38*t54
         t277 = 4._dp*t266*t32 + 5._dp*t270*t49 + 6._dp*t274*t65
         t278 = t212*t277
         t279 = t171*t278
         t281 = 2._dp*t249*t70 + t11*t263 - t210*t279
         dQndrho = f94*t281*t77
         t283 = t78*f89
         t286 = -t176 - t182 + t209 - t241
         t289 = t14*t87
         t290 = t289*t6
         t293 = t179*t87
         t296 = F2*t206
         t299 = t296*t68 - t85*t239
         t300 = t15*t299
         t302 = -0.2e1_dp/0.3e1_dp*t170*t290 - (2._dp*t11*t293) + (t11 &
                                                                   *t300)
         t303 = f12*t302
         t305 = C*t108
         t306 = t305*t286
         t309 = t14*t102
         t310 = t309*t6
         t313 = t179*t102
         t316 = t96*t183
         t321 = t99*t196
         t326 = -0.2e1_dp/0.3e1_dp*t316*t185 - (2._dp*t97*t189) - 0.4e1_dp &
                /0.3e1_dp*t321*t199 - (4._dp*t100*t203)
         t327 = t15*t326
         t329 = -0.2e1_dp/0.3e1_dp*t170*t310 - (2._dp*t11*t313) + (t11 &
                                                                   *t327)
         t330 = t329*E
         t331 = t330*t108
         t332 = t92**2
         t333 = 0.1e1_dp/t332
         t334 = t333*t286
         t340 = t110*t105
         t341 = t93*t118
         t342 = t341*t286
         t345 = 0.1e1_dp/t117/t116
         t346 = t83*t345
         t350 = -t176 - t182 + t209 - t241 - 0.2e1_dp/0.3e1_dp*t72*t44*t6
         t351 = t115*t350
         t352 = t346*t351
         t355 = t93*t122
         t359 = 0.1e1_dp/t121/t113
         t360 = t83*t359
         t367 = t125*t90
         t368 = t93*t129
         t369 = t368*t286
         t372 = 0.1e1_dp/t128/t127
         t373 = t83*t372
         t374 = t114*t350
         t375 = t373*t374
         t379 = (-t110*t329*t83*t118 + t340*t342 + 0.5e1_dp/0.2e1_dp*t340 &
                 *t352 + t81*t355*t286 + t81*t360*t350/0.2e1_dp - t125 &
                 *t302*t83*t129 + t367*t369 + 0.3e1_dp/0.2e1_dp*t367*t375)* &
                omega
         t382 = t27*r3*t5
         t388 = t136*t105
         t389 = t108*t118
         t390 = t389*t286
         t393 = t93*t345
         t394 = t393*t351
         t400 = t140*t90
         t401 = t108*t129
         t402 = t401*t286
         t405 = t93*t372
         t406 = t405*t374
         t410 = (-t136*t329*t93*t118 + (2._dp*t388*t390) + 0.5e1_dp/ &
                 0.2e1_dp*(t388)*(t394) - t140*t302*t93*t129 + (2._dp &
                                                                *t400*t402) + 0.3e1_dp/0.2e1_dp*(t400)*(t406))*t145
         t412 = t149*t13
         t415 = t106*t333
         t417 = t155*t44*t286
         t420 = t108*t345
         t421 = t106*t420
         t422 = t154*t44
         t423 = t422*t351
         t426 = t106*t389
         t427 = t154*t219
         t428 = t427*t6
         t431 = dexeirho(Q, dQrho)
         t433 = t160**2
         t434 = 0.1e1_dp/t433
         t435 = t75*t434
         t437 = t246*t161 - t435*t246
         t438 = 0.1e1_dp/t75
         t439 = t437*t438
         t443 = -t81*t93*t286 + t303*t94 - (2._dp*t91*t306) + t331 &
                - (3._dp*t106*t334) + t379*t134 - t133*t382/0.3e1_dp + t410 &
                *t151 - t146*t412 - t331*t156 + (3._dp*t415*t417) + 0.5e1_dp &
                /0.2e1_dp*t421*t423 + 0.5e1_dp/0.3e1_dp*t426*t428 + t158*(t431 &
                                                                          + t439*t160)
         t444 = t443*Clda
         e_rho = e_rho + (-0.4e1_dp/0.3e1_dp*t283*t167 - t80*t444)*sx
         t450 = F2*t261
         t453 = t450*t68 - t85*t278
         t454 = t15*t453
         t456 = 2._dp*t249*t88 + t11*t454
         t457 = f12*t456
         t459 = t305*t281
         t464 = g2*ndrho
         t465 = t464*t3
         t468 = g3*t256
         t469 = t468*t24
         t472 = 2._dp*t465*t19 + 4._dp*t469*t32
         t473 = t15*t472
         t475 = 2._dp*t249*t103 + t11*t473
         t476 = t475*E
         t477 = t476*t108
         t478 = t333*t281
         t484 = t341*t281
         t486 = t115*t281
         t487 = t346*t486
         t498 = t368*t281
         t500 = t114*t281
         t501 = t373*t500
         t505 = (-t110*t475*t83*t118 + t340*t484 + 0.5e1_dp/0.2e1_dp*t340 &
                 *t487 + t81*t355*t281 + t81*t360*t281/0.2e1_dp - t125 &
                 *t456*t83*t129 + t367*t498 + 0.3e1_dp/0.2e1_dp*t367*t501)* &
                omega
         t510 = t389*t281
         t513 = t393*t486
         t519 = t401*t281
         t522 = t405*t500
         t526 = (-t136*t475*t93*t118 + (2._dp*t388*t510) + 0.5e1_dp/ &
                 0.2e1_dp*(t388)*(t513) - t140*t456*t93*t129 + (2._dp &
                                                                *t400*t519) + 0.3e1_dp/0.2e1_dp*(t400)*(t522))*t145
         t530 = t155*t44*t281
         t533 = t422*t486
         t536 = dexeindrho(Q, dQndrho)
         t539 = t281*t161 - t435*t281
         t540 = t539*t438
         t544 = -t81*t93*t281 + t457*t94 - (2._dp*t91*t459) + t477 &
                - (3._dp*t106*t478) + t505*t134 + t526*t151 - t477*t156 &
                + (3._dp*t415*t530) + 0.5e1_dp/0.2e1_dp*t421*t533 + t158*(t536 &
                                                                          + t540*t160)
         t545 = t544*Clda
         e_ndrho = e_ndrho + (-t80*t545)*sx
      END IF

   END SUBROUTINE xwpbe_lda_calc_4

! **************************************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE xwpbe_lsd_info(reference, shortform, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "Jochen Heyd and Gustavo E. Scuseria, J. Chem. Phys., 120, 7274 {LSD version}"
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "shortrange part of PBE exchange {LSD}"
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho_spin = .TRUE.
         needs%norm_drho_spin = .TRUE.
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1
   END SUBROUTINE xwpbe_lsd_info

! **************************************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lsd
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param xwpbe_params input parameters (scaling,omega)
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      The current version provides code for derivatives up to second order.
!>      Using the maple sheet in cp2k/doc it is straightforward to produce routines
!>      for higher derivatives.
! **************************************************************************************************
   SUBROUTINE xwpbe_lsd_eval(rho_set, deriv_set, order, xwpbe_params)

      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(IN)                                :: order
      TYPE(section_vals_type), POINTER                   :: xwpbe_params

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xwpbe_lsd_eval'

      INTEGER                                            :: handle, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp)                                      :: epsilon_norm_drho, epsilon_rho, omega, &
                                                            sx, sx0
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: dummy, e_0, e_ndrhoa, e_ndrhob, e_rhoa, &
                                                            e_rhob, norm_drhoa, norm_drhob, rhoa, &
                                                            rhob
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)

      CALL cite_reference(Heyd2004)

      CALL section_vals_val_get(xwpbe_params, "SCALE_X", r_val=sx)
      CALL section_vals_val_get(xwpbe_params, "SCALE_X0", r_val=sx0)
      CALL section_vals_val_get(xwpbe_params, "OMEGA", r_val=omega)

      CALL xc_rho_set_get(rho_set, rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
                          norm_drhob=norm_drhob, local_bounds=bo, rho_cutoff=epsilon_rho, &
                          drho_cutoff=epsilon_norm_drho)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      dummy => rhoa

      e_0 => dummy
      e_rhoa => dummy
      e_rhob => dummy
      e_ndrhoa => dummy
      e_ndrhob => dummy

      IF (order >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
      IF (order >= 1 .OR. order == -1) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rhoa)
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rhob)
         deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa)
         deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_ndrhob)
      END IF
      IF (order > 1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF

!$OMP     PARALLEL DEFAULT(NONE) &
!$OMP     SHARED(npoints, order, rhoa, norm_drhoa, e_0, e_rhoa, e_ndrhoa) &
!$OMP     SHARED(epsilon_rho, sx, sx0, omega) &
!$OMP     SHARED(rhob, norm_drhob, e_rhob, e_ndrhob)

      !Call lsd_calc for alpha - and beta - spins

      CALL xwpbe_lsd_calc(npoints, order, rho=rhoa, norm_drho=norm_drhoa, &
                          e_0=e_0, e_rho=e_rhoa, e_ndrho=e_ndrhoa, &
                          epsilon_rho=epsilon_rho, &
                          sx=sx, sx0=sx0, omega=omega)
      CALL xwpbe_lsd_calc(npoints, order, rho=rhob, norm_drho=norm_drhob, &
                          e_0=e_0, e_rho=e_rhob, e_ndrho=e_ndrhob, &
                          epsilon_rho=epsilon_rho, &
                          sx=sx, sx0=sx0, omega=omega)

!$OMP     END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE xwpbe_lsd_eval

! **************************************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lsd.
!> \param npoints ...
!> \param order degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param rho , ndrho: density and norm of the density gradient
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param epsilon_rho ...
!> \param sx , sx0: scaling factor for omega!=0 and omega=0
!> \param sx0 ...
!> \param omega screening parameter
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - The lsd part is calculated using the spin-scaling relations for the
!>        exchange energy:
!>
!>               Ex[na,nb] = 0.5 * Ex[2*na] + 0.5 * Ex[2*nb].
!>
!>      - In order to avoid numerical instabilities, this routine calls different
!>        subroutines. There are 4 routines for the case omega!=0 and 2 routines
!>        for omega=0.
! **************************************************************************************************
   SUBROUTINE xwpbe_lsd_calc(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho, &
                             epsilon_rho, sx, sx0, omega)

      INTEGER, INTENT(in)                                :: npoints, order
      REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: rho, norm_drho, e_0, e_rho, e_ndrho
      REAL(kind=dp), INTENT(in)                          :: epsilon_rho, sx, sx0, omega

      INTEGER                                            :: ip
      REAL(dp)                                           :: e_0_temp, my_ndrho, my_rho
      REAL(KIND=dp)                                      :: ss, ss2, sscale, t1, t2, t3, t4, t5, t6, &
                                                            t7, t8, ww

!$OMP     DO

      DO ip = 1, npoints
         !According to spin-scaling relation, we need twice the density and its gradient
         my_rho = 2.0_dp*MAX(rho(ip), 0.0_dp)
         IF (my_rho > epsilon_rho) THEN
            my_ndrho = 2.0_dp*MAX(norm_drho(ip), 0.0_dp)

            !Do some precalculation in order to catch the correct branch afterwards
            sscale = 1.0_dp
            t1 = pi**2
            t2 = t1*my_rho
            t3 = t2**(0.1e1_dp/0.3e1_dp)
            t4 = 0.1e1_dp/t3
            t5 = omega*t4
            ww = 0.6933612743506347048433524e0_dp*t5
            t6 = my_ndrho*t4
            t7 = 0.1e1_dp/my_rho
            t8 = t7*sscale
            ss = 0.3466806371753173524216762e0_dp*t6*t8
            IF (ss > scutoff) THEN
               ss2 = ss*ss
               sscale = ((smax)*ss2 - (sconst))/(ss2*ss)
            END IF
            e_0_temp = 0.0_dp
            IF (sx0 /= 0.0_dp) THEN
               !original PBE hole
               IF (ss*sscale > gcutoff) THEN
                  CALL xwpbe_lda_calc_0(e_0_temp, e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, sscale, sx0, order)
               ELSE
                  CALL xwpbe_lda_calc_01(e_0_temp, e_rho(ip), e_ndrho(ip), my_rho, &
                                         my_ndrho, sscale, sx0, order)
               END IF
               !According to spin-scaling relation, we need only half of the energy
               e_0(ip) = e_0(ip) + 0.5_dp*e_0_temp
            END IF
            e_0_temp = 0.0_dp
            IF (sx /= 0.0_dp) THEN
               IF (ww < wcutoff .AND. ss*sscale > gcutoff) THEN
                  CALL xwpbe_lda_calc_1(e_0_temp, e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               ELSE IF (ww < wcutoff .AND. ss*sscale <= gcutoff) THEN
                  CALL xwpbe_lda_calc_2(e_0_temp, e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               ELSE IF (ww >= wcutoff .AND. ss*sscale > gcutoff) THEN
                  CALL xwpbe_lda_calc_3(e_0_temp, e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               ELSE
                  CALL xwpbe_lda_calc_4(e_0_temp, e_rho(ip), e_ndrho(ip), my_rho, &
                                        my_ndrho, omega, sscale, sx, order)
               END IF
               !According to spin-scaling relation, we need only half of the energy
            END IF
            e_0(ip) = e_0(ip) + 0.5_dp*e_0_temp
         END IF
      END DO

!$OMP     END DO

   END SUBROUTINE xwpbe_lsd_calc

! **************************************************************************************************
!> \brief These functions evaluate products exp(x)*Ei(x) and pi*exp(x)*erfc(sqrt(x)),
!>      as well as their derivatives with respect to various combinations of
!>      rho and norm_drho.
!> \param Q , dQrho, dQndrho, d2Qrhondrho :
!>        Argument Q and derivatives with respect to various combinations of
!>        rho and norm_drho
!> \return ...
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>     - In order to avoid numerical instabilities, these routines use Taylor-
!>       expansions for the above core-products for large arguments.
!>     - When adapting this module for higher order derivatives, appropriate
!>       functions have to be provided!
! **************************************************************************************************
   FUNCTION exei(Q)
      REAL(dp), INTENT(IN)                               :: Q
      REAL(dp)                                           :: exei

      exei = 0.0_dp
      IF (Q < expcutoff) THEN
         !Use exact product
         exei = EXP(Q)*expint(1, Q)
      ELSE
         !Use approximation
         exei = (1._dp/Q)*(Q*Q + exei1*Q + exei2)/(Q*Q + exei3*Q + exei4)
      END IF
   END FUNCTION exei

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \return ...
! **************************************************************************************************
   FUNCTION exer(Q)
      REAL(dp), INTENT(IN)                               :: Q
      REAL(dp)                                           :: exer

      REAL(dp)                                           :: Q3, Q5

      exer = 0.0_dp
      IF (Q < expcutoff) THEN
         !Use exact expression
         exer = pi*EXP(Q)*erfc(SQRT(Q))
      ELSE
         !Use approximation
         Q3 = Q*Q*Q
         Q5 = Q3*Q*Q
         exer = pi*(1.0_dp/SQRT(Q*pi) - 1.0_dp/(2.0_dp*SQRT(pi*Q3)) + 3.0_dp/(4.0_dp*(SQRT(pi*Q5))))
      END IF
   END FUNCTION exer

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION dexeirho(Q, dQrho)
      REAL(dp), INTENT(IN)                               :: Q, dQrho
      REAL(dp)                                           :: dexeirho

      dexeirho = dQrho*(exei(Q) - 1.0_dp/Q)
   END FUNCTION dexeirho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION dexeindrho(Q, dQndrho)
      REAL(dp), INTENT(IN)                               :: Q, dQndrho
      REAL(dp)                                           :: dexeindrho

      dexeindrho = dQndrho*(exei(Q) - 1.0_dp/Q)
   END FUNCTION dexeindrho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION dexerrho(Q, dQrho)
      REAL(dp), INTENT(IN)                               :: Q, dQrho
      REAL(dp)                                           :: dexerrho

      dexerrho = dQrho*exer(Q) - dQrho*rootpi/SQRT(Q)
   END FUNCTION dexerrho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION dexerndrho(Q, dQndrho)
      REAL(dp), INTENT(IN)                               :: Q, dQndrho
      REAL(dp)                                           :: dexerndrho

      dexerndrho = dQndrho*exer(Q) - dQndrho*rootpi/SQRT(Q)
   END FUNCTION dexerndrho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param d2Qrhorho ...
!> \return ...
! **************************************************************************************************
   FUNCTION d2exeirhorho(Q, dQrho, d2Qrhorho)
      REAL(dp), INTENT(IN)                               :: Q, dQrho, d2Qrhorho
      REAL(dp)                                           :: d2exeirhorho

      d2exeirhorho = exei(Q)*(d2Qrhorho + dQrho*dQrho) + &
                     1.0_dp/(Q*Q)*(-Q*dQrho*dQrho - Q*d2Qrhorho + dQrho*dQrho)
   END FUNCTION d2exeirhorho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param d2Qrhorho ...
!> \return ...
! **************************************************************************************************
   FUNCTION d2exerrhorho(Q, dQrho, d2Qrhorho)
      REAL(dp), INTENT(IN)                               :: Q, dQrho, d2Qrhorho
      REAL(dp)                                           :: d2exerrhorho

      REAL(dp)                                           :: pi12, Q12

      Q12 = SQRT(Q)
      pi12 = rootpi

      d2exerrhorho = exer(Q)*(d2Qrhorho + dQrho*dQrho) - dQrho*dQrho/(pi12*Q12) + &
                     0.5_dp*dQrho*dQrho/(pi12*Q*Q12) - d2Qrhorho/(pi12*Q12)
   END FUNCTION d2exerrhorho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param dQndrho ...
!> \param d2Qrhondrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION d2exeirhondrho(Q, dQrho, dQndrho, d2Qrhondrho)
      REAL(dp), INTENT(IN)                               :: Q, dQrho, dQndrho, d2Qrhondrho
      REAL(dp)                                           :: d2exeirhondrho

      d2exeirhondrho = exei(Q)*(d2Qrhondrho + dQrho*dQndrho) - &
                       1.0_dp/Q*(dQrho*dQndrho + d2Qrhondrho) + 1.0_dp/(Q*Q)*dQrho*dQndrho
   END FUNCTION d2exeirhondrho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param dQndrho ...
!> \param d2Qrhondrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION d2exerrhondrho(Q, dQrho, dQndrho, d2Qrhondrho)
      REAL(dp), INTENT(IN)                               :: Q, dQrho, dQndrho, d2Qrhondrho
      REAL(dp)                                           :: d2exerrhondrho

      REAL(dp)                                           :: pi12, Q12

      Q12 = SQRT(Q)
      pi12 = rootpi

      d2exerrhondrho = exer(Q)*(d2Qrhondrho + dQrho*dQndrho) - 1.0_dp/(pi12*Q12)*dQrho*dQndrho &
                       + 0.5_dp/(pi12*Q12*Q)*dQrho*dQndrho - 1.0_dp/(pi12*Q12)*d2Qrhondrho
   END FUNCTION d2exerrhondrho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \param d2Qndrhondrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION d2exeindrhondrho(Q, dQndrho, d2Qndrhondrho)
      REAL(dp), INTENT(IN)                               :: Q, dQndrho, d2Qndrhondrho
      REAL(dp)                                           :: d2exeindrhondrho

      d2exeindrhondrho = exei(Q)*(d2Qndrhondrho + dQndrho*dQndrho) + &
                         1.0_dp/(Q*Q)*(-Q*dQndrho*dQndrho - Q*d2Qndrhondrho + dQndrho*dQndrho)
   END FUNCTION d2exeindrhondrho

! **************************************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \param d2Qndrhondrho ...
!> \return ...
! **************************************************************************************************
   FUNCTION d2exerndrhondrho(Q, dQndrho, d2Qndrhondrho)
      REAL(dp), INTENT(IN)                               :: Q, dQndrho, d2Qndrhondrho
      REAL(dp)                                           :: d2exerndrhondrho

      REAL(dp)                                           :: pi12, Q12

      Q12 = SQRT(Q)
      pi12 = rootpi

      d2exerndrhondrho = exer(Q)*(d2Qndrhondrho + dQndrho*dQndrho) - dQndrho*dQndrho/(pi12*Q12) &
                         + 0.5_dp*dQndrho*dQndrho/(pi12*Q*Q12) - d2Qndrhondrho/(pi12*Q12)
   END FUNCTION d2exerndrhondrho

END MODULE xc_xwpbe

