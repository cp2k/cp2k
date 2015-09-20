!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates short range exchange part for wPBE functional and averaged
!>      PBE exchange-hole functional (omega = 0.0 )
!> \par History
!>      Manuel Guidon (05.2007)  : initial version
!> \author Manuel Guidon (05.2007)
! *****************************************************************************
MODULE xc_xwpbe
  USE bibliography,                    ONLY: Heyd2004,&
                                             cite_reference
  USE erf_fn,                          ONLY: erfc
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             rootpi
  USE mathlib,                         ONLY: expint
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

! *** Global parameters ***

  PUBLIC :: xwpbe_lda_info, xwpbe_lda_eval, xwpbe_lsd_info, &
            xwpbe_lsd_eval

  REAL(KIND=dp), PARAMETER :: alpha1 = -1.128223946706117_dp, &
                              alpha2 =  1.452736265762971_dp,&
                              alpha3 = -1.243162299390327_dp, &
                              alpha4 =  0.971824836115601_dp, &
                              alpha5 = -0.568861079687373_dp, &
                              alpha6 =  0.246880514820192_dp, &
                              alpha7 = -0.065032363850763_dp, &
                              alpha8 =  0.008401793031216_dp
  REAL(KIND=dp), PARAMETER :: beta  = 1.455915450052607_dp, &
                              beta2 = 2.0_dp
  REAL(KIND=dp), PARAMETER :: a1 = 0.00979681_dp, &
                              a2 = 0.04108340_dp, &
                              a3 = 0.18744000_dp, &
                              a4 = 0.00120824_dp, &
                              a5 = 0.0347188_dp
  REAL(KIND=dp), PARAMETER :: A  = 1.0161144_dp, &
                              B  = -0.37170836_dp, &
                              C  = -0.077215461_dp, &
                              DD = 0.57786348_dp, &
                              E  = -0.051955731_dp, &
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
                              f6561512 = 6561.0_dp/512.0_dp, f8132 = 81.0_dp/32.0_dp,&
                              f8164 = 81.0_dp/64.0_dp, f729128 = 729.0_dp/128.0_dp,&
                              f52 = 2.5_dp, f94 = 9.0_dp/4.0_dp, f916 = 9.0_dp/16.0_dp,&
                              f89 = 8.0_dp/9.0_dp, f2187256 = 2187.0_dp/256.0_dp, &
                              r1 = 1.0_dp, f98 = 9.0_dp/8.0_dp, r15 = 15.0_dp,&
                              r3 = 3.0_dp, r4 = 4.0_dp, r16 = 16.0_dp, r8 = 8.0_dp, &
                              r6 = 6.0_dp, r2=2.0_dp

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xwpbe'

  CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xwpbe_lda_info ( reference, shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    IF ( PRESENT ( reference ) ) THEN
      reference = "Jochen Heyd and Gustavo E. Scuseria, J. Chem. Phys., 120, 7274 {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "shortrange part of PBE exchange {LDA}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=2
  END SUBROUTINE xwpbe_lda_info

! *****************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evalated,
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
! *****************************************************************************
  SUBROUTINE xwpbe_lda_eval ( rho_set, deriv_set, order, xwpbe_params)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: xwpbe_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xwpbe_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, omega, sx, sx0
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrho, &
                                                e_ndrho_ndrho, e_ndrho_rho, &
                                                e_rho, e_rho_rho, norm_drho, &
                                                rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(bo)
    CALL cite_reference(Heyd2004)

    CALL section_vals_val_get(xwpbe_params,"SCALE_X",r_val=sx)
    CALL section_vals_val_get(xwpbe_params,"SCALE_X0",r_val=sx0)
    CALL section_vals_val_get(xwpbe_params,"OMEGA",r_val=omega)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)


    CALL xc_rho_set_get(rho_set,rho=rho,&
         norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
         drho_cutoff=epsilon_norm_drho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPASSERT(stat==0)
    ELSE
      dummy=> rho
    END IF

    e_0 => dummy
    e_rho => dummy
    e_ndrho => dummy
    e_rho_rho => dummy
    e_ndrho_rho => dummy
    e_ndrho_ndrho => dummy

    IF (order>=0) THEN
      deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (order>=1.OR.order==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rho)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
    END IF
    IF (order>=2.OR.order==-2) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rho_rho)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho)
      deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)
    END IF
    IF (order>2.OR.order<-2) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 2 not implemented", &
            error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp shared(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho) &
    !$omp shared(e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, epsilon_rho) &
    !$omp shared(sx, sx0, omega)

    CALL xwpbe_lda_calc(npoints,order,rho=rho, norm_drho=norm_drho,&
         e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
         e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
         epsilon_rho=epsilon_rho,sx=sx,sx0=sx0,omega=omega)

    !$omp end parallel

    IF (cp_debug) THEN
      DEALLOCATE(dummy,stat=stat)
      CPASSERT(stat==0)
    ELSE
      NULLIFY(dummy)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xwpbe_lda_eval

! *****************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lda
!> \param npoints ...
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param rho , ndrho: density and norm of the density gradient
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
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
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho,&
                            e_rho_rho, e_ndrho_rho, e_ndrho_ndrho,&
                            epsilon_rho, sx, sx0, omega)

    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: rho, norm_drho, e_0, e_rho, &
                                                e_ndrho, e_rho_rho, &
                                                e_ndrho_rho, e_ndrho_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, sx0, omega

    CHARACTER(len=*), PARAMETER :: routineN = 'xwpbe_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: my_ndrho, my_rho
    REAL(KIND=dp)                            :: ss, ss2, sscale, t1, t2, t3, &
                                                t4, t5, t6, t7, t8, ww

    !$omp do

    DO ip =1,npoints
      my_rho = MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho = MAX(norm_drho(ip),0.0_dp)

        !Do some precalculation in order to catch the correct branch afterwards
        sscale = 1.0_dp
        t1 = pi ** 2
        t2 = t1 * my_rho
        t3 = t2 ** (0.1e1_dp / 0.3e1_dp)
        t4 = 0.1e1_dp / t3
        t5 = omega * t4
        ww = 0.6933612743506347048433524e0_dp * t5
        t6 = my_ndrho * t4
        t7 = 0.1e1_dp / my_rho
        t8 = t7 * sscale
        ss = 0.3466806371753173524216762e0_dp * t6 * t8
        IF( ss > scutoff) THEN
          ss2 = ss*ss
          sscale = (smax*ss2-sconst)/(ss2*ss)
        END IF

        IF(sx0/=0.0_dp) THEN
          !original PBE hole
           IF(ss*sscale>gcutoff) THEN
             CALL xwpbe_lda_calc_0(e_0(ip), e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                   e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                   my_ndrho, sscale, sx0, order)
           ELSE
             CALL xwpbe_lda_calc_01(e_0(ip), e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                   e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                   my_ndrho, sscale, sx0, order)
           END IF
        END IF

        IF( sx /= 0.0_dp ) THEN
          IF(ww<wcutoff .AND. ss*sscale>gcutoff) THEN
            CALL xwpbe_lda_calc_1(e_0(ip), e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          ELSE IF(ww<wcutoff .AND. ss*sscale<=gcutoff) THEN
            CALL xwpbe_lda_calc_2(e_0(ip), e_rho(ip), e_ndrho(ip), e_rho_rho(ip), &
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip),my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          ELSE IF(ww>=wcutoff .AND. ss*sscale>gcutoff) THEN
            CALL xwpbe_lda_calc_3(e_0(ip), e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          ELSE
            CALL xwpbe_lda_calc_4(e_0(ip), e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          END IF
        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE xwpbe_lda_calc

! *****************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx0 scaling factor
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the exact functional for omega=0.
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc_0(e_0, e_rho, e_ndrho, e_rho_rho, e_ndrho_rho, &
                              e_ndrho_ndrho, rho, ndrho, sscale, sx0, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho, &
                                                e_rho_rho, e_ndrho_rho, &
                                                e_ndrho_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, sscale, sx0
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: d2Qndrhondrho, d2Qrhondrho, d2Qrhorho, dQndrho, dQrho, &
      Q, t1, t10, t100, t1003, t1007, t1018, t1023, t1024, t103, t1037, t105, &
      t1055, t1056, t1057, t1059, t106, t1060, t1073, t1079, t1082, t109, &
      t11, t110, t1100, t111, t1110, t1114, t1117, t112, t114, t1143, t116, &
      t1165, t117, t118, t119, t12, t120, t1202, t121, t1215, t122, t123, &
      t1242, t125, t1258, t126, t1263, t127, t1286, t13, t130, t1316, t134, &
      t1347, t135, t1352, t1358, t136, t1362, t1365, t1369, t1372, t138, &
      t1382, t1388, t139, t1392, t1395, t14, t1412, t142, t143, t145, t1455, &
      t146, t1465, t147, t148, t149, t15, t150
    REAL(KIND=dp) :: t1510, t152, t1547, t156, t1561, t158, t159, t16, t160, &
      t163, t164, t165, t166, t169, t170, t172, t173, t174, t175, t176, t18, &
      t180, t183, t184, t185, t187, t188, t19, t190, t191, t192, t193, t199, &
      t2, t20, t200, t201, t202, t203, t207, t209, t21, t215, t218, t219, &
      t22, t220, t222, t223, t224, t225, t227, t228, t231, t233, t236, t237, &
      t24, t240, t241, t242, t245, t246, t247, t249, t25, t250, t253, t254, &
      t258, t26, t261, t262, t263, t265, t266, t269, t272, t274, t275, t278, &
      t28, t281, t282, t284, t285, t287, t288, t29, t290, t291, t294, t297, &
      t3, t300, t301, t303, t304, t306
    REAL(KIND=dp) :: t307, t31, t314, t32, t321, t323, t326, t327, t33, t330, &
      t334, t335, t336, t337, t338, t339, t34, t340, t342, t343, t344, t345, &
      t347, t355, t356, t357, t358, t359, t36, t361, t362, t363, t367, t368, &
      t37, t372, t374, t376, t377, t381, t383, t384, t387, t388, t39, t390, &
      t391, t394, t397, t398, t4, t40, t400, t401, t402, t404, t405, t406, &
      t408, t409, t410, t411, t412, t413, t414, t415, t416, t417, t418, t42, &
      t423, t425, t426, t429, t430, t432, t433, t435, t44, t440, t443, t444, &
      t446, t447, t449, t46, t463, t465, t47, t471, t472, t475, t476, t478, &
      t479, t480, t484, t485, t489
    REAL(KIND=dp) :: t491, t495, t497, t5, t500, t501, t504, t505, t507, &
      t508, t510, t511, t513, t514, t515, t52, t520, t53, t530, t537, t54, &
      t540, t542, t546, t550, t553, t557, t56, t566, t569, t570, t577, t579, &
      t58, t59, t6, t61, t616, t620, t621, t627, t628, t63, t632, t647, t65, &
      t655, t657, t66, t663, t67, t678, t68, t685, t69, t7, t70, t707, t71, &
      t716, t72, t73, t735, t74, t744, t751, t755, t76, t761, t77, t778, t78, &
      t784, t788, t79, t791, t8, t80, t81, t819, t824, t83, t84, t854, t856, &
      t857, t86, t867, t872, t875, t878, t88, t887, t888, t889, t89, t9, &
      t905, t91, t910, t911, t92, t923, t924
    REAL(KIND=dp) :: t93, t930, t933, t94, t95, t952, t956, t968, t97, t975, &
      t98, t983

    IF( order >= 0 ) THEN
      t1 = ndrho ** 2
      t2 = a1 * t1
      t3 = r2 ** 2
      t4 = 0.1e1_dp / t3
      t5 = t2 * t4
      t6 = pi ** 2
      t7 = r3 * t6
      t8 = t7 * rho
      t9 = t8 ** (0.1e1_dp / 0.3e1_dp)
      t10 = t9 ** 2
      t11 = 0.1e1_dp / t10
      t12 = rho ** 2
      t13 = 0.1e1_dp / t12
      t14 = t11 * t13
      t15 = sscale ** 2
      t16 = t14 * t15
      t18 = t1 ** 2
      t19 = a2 * t18
      t20 = t3 ** 2
      t21 = 0.1e1_dp / t20
      t22 = t19 * t21
      t24 = 0.1e1_dp / t9 / t8
      t25 = t12 ** 2
      t26 = 0.1e1_dp / t25
      t28 = t15 ** 2
      t29 = t24 * t26 * t28
      t31 = t5 * t16 + t22 * t29
      t32 = f94 * t31
      t33 = a3 * t18
      t34 = t33 * t21
      t36 = t18 * ndrho
      t37 = a4 * t36
      t39 = 0.1e1_dp / t20 / r2
      t40 = t37 * t39
      t42 = 0.1e1_dp / t10 / t8
      t44 = 0.1e1_dp / t25 / rho
      t46 = t28 * sscale
      t47 = t42 * t44 * t46
      t52 = 0.1e1_dp / t20 / t3
      t53 = a5 * t18 * t1 * t52
      t54 = r3 ** 2
      t56 = t6 ** 2
      t58 = 0.1e1_dp / t54 / t56
      t59 = t25 ** 2
      t61 = t28 * t15
      t63 = t58 / t59 * t61
      t65 = r1 + t34 * t29 + t40 * t47 + t53 * t63
      t66 = 0.1e1_dp / t65
      t67 = t66 * t1
      t68 = t32 * t67
      t69 = t4 * t11
      t70 = t13 * t15
      t71 = 0.1e1_dp / A
      t72 = t70 * t71
      t73 = t69 * t72
      Q = t68 * t73
      t74 = rho ** (0.1e1_dp / 0.3e1_dp)
      t76 = t74 * rho * f89
      t77 = B * f12
      t78 = t1 * t4
      t79 = t78 * t11
      t80 = t31 * t66
      t81 = t70 * t80
      t83 = t79 * t81 + DD
      t84 = 0.1e1_dp / t83
      t86 = F2 * t31
      t88 = F1 + t86 * t66
      t89 = t70 * t88
      t91 = t79 * t89 + r1
      t92 = f12 * t91
      t93 = t83 ** 2
      t94 = 0.1e1_dp / t93
      t95 = C * t94
      t97 = f34 * pi
      t98 = rootpi
      t100 = r6 * C
      t103 = r4 * B
      t105 = r8 * A
      t106 = t93 * t83
      t109 = t98 * (r15 * E + t100 * t91 * t83 + t103 * t93 + t105 * t106)
      t110 = 0.1e1_dp / r16
      t111 = SQRT(t83)
      t112 = t111 * t106
      t114 = t110 / t112
      t116 = SQRT(A)
      t117 = EXP(Q)
      t118 = t116 * t117
      t119 = f32 * ndrho
      t120 = 0.1e1_dp / r2
      t121 = t119 * t120
      t122 = 0.1e1_dp / t9
      t123 = 0.1e1_dp / rho
      t125 = t80 * t71
      t126 = SQRT(t125)
      t127 = sscale * t126
      t130 = erfc(t121 * t122 * t123 * t127)
      t134 = 0.1e1_dp / f1516
      t135 = (t97 + t109 * t114 - t97 * t118 * t130) * t134
      t136 = 0.1e1_dp / t98
      t138 = 0.1e1_dp / E
      t139 = t136 * t112 * t138
      t142 = (-t135 * t139 + r1) * E
      t143 = 0.1e1_dp / t106
      t145 = f12 * A
      t146 = exei(Q)
      t147 = t78 * t14
      t148 = t15 * t31
      t149 = t66 * t84
      t150 = t148 * t149
      t152 = LOG(t147 * t150)
      t156 = (t77 * t84 + t92 * t95 + t142 * t143 + t145 * (t146 + t152))  &
      * Clda
      e_0 = e_0 + ( -t76 * t156 ) * sx0
    END IF
    IF(order >= 1 .OR. order == -1 ) THEN
      t158 = t4 * t42
      t159 = t2 * t158
      t160 = t70 * t7
      t163 = t12 * rho
      t164 = 0.1e1_dp / t163
      t165 = t11 * t164
      t166 = t165 * t15
      t169 = t54 * t56
      t170 = t169 * t12
      t172 = 0.1e1_dp / t9 / t170
      t173 = t21 * t172
      t174 = t19 * t173
      t175 = t26 * t28
      t176 = t175 * t7
      t180 = t24 * t44 * t28
      t183 = -0.2e1_dp / 0.3e1_dp * t159 * t160 - (2._dp * t5 * t166) - 0.4e1_dp /  &
      0.3e1_dp * t174 * t176 - (4._dp * t22 * t180)
      t184 = f94 * t183
      t185 = t184 * t67
      t187 = t65 ** 2
      t188 = 0.1e1_dp / t187
      t190 = t188 * t1 * t4
      t191 = t32 * t190
      t192 = t15 * t71
      t193 = t33 * t173
      t199 = 0.1e1_dp / t10 / t170
      t200 = t39 * t199
      t201 = t37 * t200
      t202 = t44 * t46
      t203 = t202 * t7
      t207 = 0.1e1_dp / t25 / t12
      t209 = t42 * t207 * t46
      t215 = t58 / t59 / rho * t61
      t218 = -0.4e1_dp / 0.3e1_dp * t193 * t176 - (4._dp * t34 * t180) - 0.5e1_dp  &
      / 0.3e1_dp * t201 * t203 - (5._dp * t40 * t209) - (8._dp * t53 * t215)
      t219 = t192 * t218
      t220 = t14 * t219
      t222 = t67 * t4
      t223 = t32 * t222
      t224 = t42 * t13
      t225 = t224 * t15
      t227 = t71 * r3 * t6
      t228 = t225 * t227
      t231 = t164 * t15
      t233 = t69 * t231 * t71
      dQrho = t185 * t73 - t191 * t220 - 0.2e1_dp / 0.3e1_dp * t223 * t228 - (2._dp  &
      * t68 * t233)
      t236 = a1 * ndrho
      t237 = t236 * t4
      t240 = t1 * ndrho
      t241 = a2 * t240
      t242 = t241 * t21
      t245 = 2._dp * t237 * t16 + 4._dp * t242 * t29
      t246 = f94 * t245
      t247 = t246 * t67
      t249 = a3 * t240
      t250 = t249 * t21
      t253 = a4 * t18
      t254 = t253 * t39
      t258 = a5 * t36 * t52
      t261 = 4._dp * t250 * t29 + 5._dp * t254 * t47 + 6._dp * t258 * t63
      t262 = t192 * t261
      t263 = t14 * t262
      t265 = t66 * ndrho
      t266 = t32 * t265
      dQndrho = t247 * t73 - t191 * t263 + 2._dp * t266 * t73
      t269 = t74 * f89
      t272 = t78 * t224
      t274 = t66 * r3 * t6
      t275 = t148 * t274
      t278 = t231 * t80
      t281 = t183 * t66
      t282 = t70 * t281
      t284 = t188 * t218
      t285 = t148 * t284
      t287 = -0.2e1_dp / 0.3e1_dp * t272 * t275 - (2._dp * t79 * t278) + (t79  &
      * t282) - t147 * t285
      t288 = t94 * t287
      t290 = t15 * t88
      t291 = t290 * t7
      t294 = t231 * t88
      t297 = F2 * t183
      t300 = t297 * t66 - t86 * t284
      t301 = t70 * t300
      t303 = -0.2e1_dp / 0.3e1_dp * t272 * t291 - (2._dp * t79 * t294) + (t79  &
      * t301)
      t304 = f12 * t303
      t306 = C * t143
      t307 = t306 * t287
      t314 = t83 * t287
      t321 = t98 * (t100 * t303 * t83 + t100 * t91 * t287 + 2._dp * t103 * t314  &
      + 3._dp * t105 * t93 * t287)
      t323 = t93 ** 2
      t326 = t110 / t111 / t323
      t327 = t326 * t287
      t330 = t97 * t116
      t334 = rootpi
      t335 = 0.1e1_dp / t334
      t336 = t117 * t335
      t337 = f32 ** 2
      t338 = t337 * t1
      t339 = t338 * t69
      t340 = t70 * t125
      t342 = EXP(-t339 * t340)
      t343 = t120 * t24
      t344 = t119 * t343
      t345 = t123 * sscale
      t347 = t126 * r3 * t6
      t355 = t119 * t120 * t122
      t356 = 0.1e1_dp / t126
      t357 = t281 * t71
      t358 = t31 * t188
      t359 = t71 * t218
      t361 = t357 - t358 * t359
      t362 = t356 * t361
      t363 = t345 * t362
      t367 = t342 * (-t344 * t345 * t347 / 0.3e1_dp - t121 * t122 * t13 * t127  &
      + t355 * t363 / 0.2e1_dp)
      t368 = t336 * t367
      t372 = (t321 * t114 - 0.7e1_dp / 0.2e1_dp * t109 * t327 - (t330 * dQrho  &
      * t117 * t130) + (2._dp * t330 * t368)) * t134
      t374 = t135 * t136
      t376 = t111 * t93 * t138
      t377 = t376 * t287
      t381 = (-t372 * t139 - 0.7e1_dp / 0.2e1_dp * t374 * t377) * E
      t383 = 0.1e1_dp / t323
      t384 = t383 * t287
      t387 = dexeirho(Q,dQrho)
      t388 = t78 * t225
      t390 = t84 * r3 * t6
      t391 = t80 * t390
      t394 = t78 * t165
      t397 = t15 * t183
      t398 = t397 * t149
      t400 = t188 * t84
      t401 = t400 * t218
      t402 = t148 * t401
      t404 = t66 * t94
      t405 = t404 * t287
      t406 = t148 * t405
      t408 = -0.2e1_dp / 0.3e1_dp * t388 * t391 - (2._dp * t394 * t150) + t147  &
      * t398 - t147 * t402 - t147 * t406
      t409 = 0.1e1_dp / t1
      t410 = t408 * t409
      t411 = t3 * t10
      t412 = t410 * t411
      t413 = 0.1e1_dp / t15
      t414 = t12 * t413
      t415 = 0.1e1_dp / t31
      t416 = t415 * t65
      t417 = t416 * t83
      t418 = t414 * t417
      t423 = (-t77 * t288 + t304 * t95 - 2._dp * t92 * t307 + t381 * t143 - 3._dp  &
      * t142 * t384 + t145 * (t387 + t412 * t418)) * Clda
      e_rho = e_rho + ( -0.4e1_dp / 0.3e1_dp * t269 * t156 - t76 * t423 ) * sx0
      t425 = ndrho * t4
      t426 = t425 * t11
      t429 = t245 * t66
      t430 = t70 * t429
      t432 = t188 * t261
      t433 = t148 * t432
      t435 = 2._dp * t426 * t81 + t79 * t430 - t147 * t433
      t440 = F2 * t245
      t443 = t440 * t66 - t86 * t432
      t444 = t70 * t443
      t446 = 2._dp * t426 * t89 + t79 * t444
      t447 = f12 * t446
      t449 = t306 * t435
      t463 = t98 * (t100 * t446 * t83 + t100 * t91 * t435 + 2._dp * t103 * t83  &
      * t435 + 3._dp * t105 * t93 * t435)
      t465 = t326 * t435
      t471 = f32 * t120
      t472 = t471 * t122
      t475 = t429 * t71
      t476 = t71 * t261
      t478 = t475 - t358 * t476
      t479 = t356 * t478
      t480 = t345 * t479
      t484 = t342 * (t472 * t345 * t126 + t355 * t480 / 0.2e1_dp)
      t485 = t336 * t484
      t489 = (t463 * t114 - 0.7e1_dp / 0.2e1_dp * t109 * t465 - (t330 * dQndrho  &
      * t117 * t130) + (2._dp * t330 * t485)) * t134
      t491 = t376 * t435
      t495 = (-t489 * t139 - 0.7e1_dp / 0.2e1_dp * t374 * t491) * E
      t497 = t383 * t435
      t500 = dexeindrho(Q,dQndrho)
      t501 = t425 * t14
      t504 = t15 * t245
      t505 = t504 * t149
      t507 = t400 * t261
      t508 = t148 * t507
      t510 = t404 * t435
      t511 = t148 * t510
      t513 = 2._dp * t501 * t150 + t147 * t505 - t147 * t508 - t147 * t511
      t514 = t513 * t409
      t515 = t514 * t411
      t520 = (-t77 * t94 * t435 + t447 * t95 - 2._dp * t92 * t449 + t495 * t143  &
      - 3._dp * t142 * t497 + t145 * (t500 + t515 * t418)) * Clda
      e_ndrho = e_ndrho + ( -t76 * t520 ) * sx0
    END IF
    IF( order >= 2 .OR. order == -2 ) THEN
      t530 = t11 * t26
      t537 = t54 * r3 * t56 * t6 * t163
      t540 = t21 / t9 / t537
      t542 = t175 * t169
      t546 = t44 * t28 * t7
      t550 = t24 * t207 * t28
      t553 = 0.10e2_dp / 0.9e1_dp * t2 * t4 * t199 * t70 * t169 + 0.8e1_dp / 0.3e1_dp  &
      * t159 * t231 * t7 + (6._dp * t5 * t530 * t15) + 0.28e2_dp / 0.9e1_dp *  &
      t19 * t540 * t542 + 0.32e2_dp / 0.3e1_dp * t174 * t546 + (20._dp * t22  &
      * t550)
      t557 = t184 * t190
      t566 = 0.1e1_dp / t187 / t65
      t569 = t32 * t566 * t1 * t4
      t570 = t218 ** 2
      t577 = t32 * t188 * t78 * t42
      t579 = t218 * r3 * t6
      t616 = 0.28e2_dp / 0.9e1_dp * t33 * t540 * t542 + 0.32e2_dp / 0.3e1_dp * t193 *  &
      t546 + (20._dp * t34 * t550) + 0.40e2_dp / 0.9e1_dp * t37 * t39 / t10 /  &
      t537 * t202 * t169 + 0.50e2_dp / 0.3e1_dp * t201 * t207 * t46 * t7 + 0.30e2_dp  &
      * t40 * t42 / t25 / t163 * t46 + (72._dp * t53 * t58 / t59 /  &
      t12 * t61)
      t620 = t199 * t13
      t621 = t620 * t15
      t627 = t42 * t164
      t628 = t627 * t15
      t632 = t26 * t15
      d2Qrhorho = f94 * t553 * t67 * t73 - (2._dp * t557 * t220) - 0.4e1_dp / 0.3e1_dp  &
      * t184 * t222 * t228 - (4._dp * t185 * t233) + (2._dp * t569 * t14  &
      * t192 * t570) + 0.4e1_dp / 0.3e1_dp * t577 * t72 * t579 + (4._dp * t191  &
      * t165 * t219) - (t191 * t14 * t192 * t616) + 0.10e2_dp / 0.9e1_dp  &
      * t223 * t621 * t71 * t54 * t56 + 0.8e1_dp / 0.3e1_dp * t223 * t628 *  &
      t227 + 0.6e1_dp * t68 * t69 * t632 * t71
      t647 = -0.4e1_dp / 0.3e1_dp * t236 * t158 * t160 - (4._dp * t237 * t166)  &
      - 0.16e2_dp / 0.3e1_dp * t241 * t173 * t176 - (16._dp * t242 * t180)
      t655 = t246 * t190
      t657 = t359 * t261
      t663 = t32 * t188 * ndrho * t4
      t678 = -0.16e2_dp / 0.3e1_dp * t249 * t173 * t176 - (16._dp * t250 * t180)  &
      - 0.25e2_dp / 0.3e1_dp * t253 * t200 * t203 - (25._dp * t254 * t209) -  &
      (48._dp * t258 * t215)
      t685 = t7 * t261
      d2Qrhondrho = (f94 * t647 * t67 * t73) - t557 * t263 + (2._dp * t184 *  &
      t265 * t73) - (t655 * t220) + (2._dp * t569 * t16 * t657) - (2._dp  &
      * t663 * t220) - (t191 * t14 * t192 * t678) - 0.2e1_dp / 0.3e1_dp  &
      * t246 * t222 * t228 + 0.2e1_dp / 0.3e1_dp * t577 * t72 * t685 - 0.4e1_dp  &
      / 0.3e1_dp * t32 * (t265) * t4 * t228 - (2._dp * t247 * t233) +  &
      (2._dp * t191 * t165 * t262) - (4._dp * t266 * t233)
      t707 = 2._dp * a1 * t4 * t16 + 12._dp * a2 * t1 * t21 * t29
      t716 = t261 ** 2
      t735 = 12._dp * a3 * t1 * t21 * t29 + 20._dp * a4 * t240 * t39 * t47 + 30._dp *  &
      a5 * t18 * t52 * t63
      d2Qndrhondrho = f94 * t707 * t67 * t73 - 2._dp * t655 * t263 + 4._dp * t246 * t265 *  &
      t73 + 2._dp * t569 * t14 * t192 * t716 - 4._dp * t663 * t263 - t191 * t14  &
      * t192 * t735 + 2._dp * t32 * t66 * t4 * t14 * t192
      t744 = t74 ** 2
      t751 = t287 ** 2
      t755 = t78 * t620
      t761 = t78 * t627
      t778 = t553 * t66
      t784 = t566 * t570
      t788 = t188 * t616
      t791 = 0.10e2_dp / 0.9e1_dp * t755 * t148 * t66 * t54 * t56 + 0.8e1_dp / 0.3e1_dp  &
      * t761 * t275 - 0.4e1_dp / 0.3e1_dp * t272 * t397 * t274 + 0.4e1_dp / 0.3e1_dp  &
      * t388 * t358 * t579 + (6._dp * t79 * t632 * t80) - (4._dp * t79  &
      * t231 * t281) + (4._dp * t394 * t285) + (t79 * t70 * t778)  &
      - 0.2e1_dp * t147 * t397 * t284 + 0.2e1_dp * t147 * t148 * t784 - t147  &
      * t148 * t788
      t819 = 0.10e2_dp / 0.9e1_dp * t755 * t290 * t169 + 0.8e1_dp / 0.3e1_dp * t761 *  &
      t291 - 0.4e1_dp / 0.3e1_dp * t272 * t15 * t300 * t7 + (6._dp * t79 * t632  &
      * t88) - 0.4e1_dp * (t79) * t231 * t300 + (t79 * t70 * (F2  &
      * t553 * t66 - 2._dp * t297 * t284 + 2._dp * t86 * t784 - t86 * t788))
      t824 = C * t383
      t854 = t323 * t83
      t856 = 0.1e1_dp / t111 / t854
      t857 = t110 * t856
      t867 = dQrho ** 2
      t872 = t97 * t116 * dQrho
      t875 = t97 * t118
      t878 = t148 * t66
      t887 = t69 * t13
      t888 = t338 * t887
      t889 = t188 * t71
      t905 = t13 * sscale
      t910 = t119 * t343 * t123
      t911 = sscale * t356
      t923 = 0.1e1_dp / t126 / t125
      t924 = t361 ** 2
      t930 = t183 * t188
      t933 = t31 * t566
      t952 = t372 * t136
      t956 = t111 * t83 * t138
      t968 = 0.1e1_dp / t854
      t975 = d2exeirhorho(Q,dQrho,d2Qrhorho)
      t983 = t66 * t143
      t1003 = t358 * t84
      t1007 = t80 * t94
      t1018 = t566 * t84
      t1023 = t78 * t16
      t1024 = t94 * t218
      t1037 = (6._dp * t78 * t530 * t150) - (4._dp * t394 * t398) + (4._dp  &
      * t394 * t402) + (2._dp * t147 * t148 * t983 * t751) - (t147  &
      * t148 * t404 * t791) + (t147 * t15 * t553 * t149) - (2._dp *  &
      t147 * t397 * t401) - (2._dp * t147 * t397 * t405) - 0.4e1_dp / 0.3e1_dp  &
      * t388 * t281 * t390 + 0.4e1_dp / 0.3e1_dp * t388 * t1003 * t579 + 0.4e1_dp  &
      / 0.3e1_dp * t388 * t1007 * t7 * t287 + 0.10e2_dp / 0.9e1_dp * (t78)  &
      * (t621) * (t80) * (t84) * (t54) * (t56) + (2._dp  &
      * t147 * t148 * t1018 * t570) + 0.2e1_dp * t1023 * t358 * t1024  &
      * t287 - (t147 * t148 * t400 * t616) + 0.8e1_dp / 0.3e1_dp * (t78)  &
      * (t628) * (t391) + (4._dp * t394 * t406)
      t1055 = t411 * t12
      t1056 = t410 * t1055
      t1057 = t31 ** 2
      t1059 = t413 / t1057
      t1060 = t65 * t83
      t1073 = (2._dp * t77 * t143 * t751) - (t77 * t94 * t791) + (f12  &
      * t819 * t95) - (4._dp * t304 * t307) + (6._dp * t92 * t824 *  &
      t751) - (2._dp * t92 * t306 * t791) + (-((t98 * (t100 * t819  &
      * t83 + 2._dp * t100 * t303 * t287 + t100 * t91 * t791 + 2._dp * t103 * t751  &
      + 2._dp * t103 * t83 * t791 + 6._dp * t105 * t83 * t751 + 3._dp * t105 * t93  &
      * t791) * t114) - (7._dp * t321 * t327) + 0.63e2_dp / 0.4e1_dp * (t109)  &
      * (t857) * (t751) - 0.7e1_dp / 0.2e1_dp * (t109) * (t326)  &
      * (t791) - t330 * d2Qrhorho * t117 * t130 - t330 * t867 * t117  &
      * t130 + (4._dp * t872 * t368) + 0.2e1_dp * t875 * t335 * (0.2e1_dp /  &
      0.3e1_dp * t338 * t158 * t13 * t878 * t227 + (2._dp * t339 * t231 *  &
      t125) - (t339 * t70 * t357) + t888 * t148 * t889 * t218) * t367  &
      + 0.2e1_dp * t330 * t336 * t342 * (0.4e1_dp / 0.9e1_dp * t119 * t120 * t172  &
      * t345 * t126 * t54 * t56 + 0.2e1_dp / 0.3e1_dp * t344 * t905 * t347  &
      - t910 * t911 * t7 * t361 / 0.3e1_dp + (2._dp * t121 * t122 * t164 *  &
      t127) - t355 * t905 * t362 - t355 * t345 * t923 * t924 / 0.4e1_dp + t355  &
      * t345 * t356 * (t778 * t71 - 2._dp * t930 * t359 + 2._dp * t933 *  &
      t71 * t570 - t358 * t71 * t616) / 0.2e1_dp)) * t134 * t139 - (7._dp  &
      * t952 * t377) - 0.35e2_dp / 0.4e1_dp * (t374) * (t956) * (t751)  &
      - 0.7e1_dp / 0.2e1_dp * (t374) * (t376) * (t791)) * E *  &
      (t143) - (6._dp * t381 * t384) + (12._dp * t142 * t968 * t751)  &
      - (3._dp * t142 * t383 * t791) + t145 * (t975 + t1037 * t409 * t411  &
      * t418 + 0.2e1_dp / 0.3e1_dp * (t410) * (t3) * (t122) *  &
      (t12) * (t413) * (t415) * (t65) * (t83) * (r3)  &
      * (t6) + (2._dp * t412 * rho * t413 * t417) - t1056 * t1059  &
      * t1060 * t183 + (t412) * t414 * (t415) * t218 * (t83)  &
      + (t412) * t414 * t416 * (t287))
      e_rho_rho = e_rho_rho + ( -0.4e1_dp / 0.9e1_dp / t744 * f89 * t156 - 0.8e1_dp / 0.3e1_dp * t269 * t423  &
      - t76 * t1073 * Clda ) * sx0
      t1079 = t143 * t287 * t435
      t1082 = t425 * t224
      t1100 = t647 * t66
      t1110 = t566 * t218 * t261
      t1114 = t188 * t678
      t1117 = -0.4e1_dp / 0.3e1_dp * t1082 * t275 - 0.2e1_dp / 0.3e1_dp * t272 * t504  &
      * t274 + 0.2e1_dp / 0.3e1_dp * t388 * t358 * t685 - (4._dp * t426 * t278)  &
      - (2._dp * t79 * t231 * t429) + (2._dp * t394 * t433) + (2._dp  &
      * t426 * t282) + (t79 * t70 * t1100) - t147 * t397 * t432 - (2._dp  &
      * t501 * t285) - t147 * t504 * t284 + 0.2e1_dp * t147 * t148 * t1110  &
      - t147 * t148 * t1114
      t1143 = -0.4e1_dp / 0.3e1_dp * t1082 * t291 - 0.2e1_dp / 0.3e1_dp * t272 * t15  &
      * t443 * t7 - (4._dp * t426 * t294) - 0.2e1_dp * t79 * t231 * t443 +  &
      (2._dp * t426 * t301) + t79 * t70 * (F2 * t647 * t66 - t297 *  &
      t432 - t440 * t284 + 2._dp * t86 * t1110 - t86 * t1114)
      t1165 = t435 * t287
      t1202 = t97 * t116 * dQndrho
      t1215 = t335 * (-2._dp * t337 * ndrho * t69 * t340 - t339 * t70 * t475  &
      + t888 * t148 * t889 * t261)
      t1242 = t245 * t188
      t1258 = (t98 * (t100 * t1143 * t83 + t100 * t303 * t435 + t100  &
      * t446 * t287 + t100 * t91 * t1117 + 2._dp * t103 * t1165 + 2._dp * t103 *  &
      t83 * t1117 + 6._dp * t105 * t314 * t435 + 3._dp * t105 * t93 * t1117) *  &
      t114) - 0.7e1_dp / 0.2e1_dp * t321 * t465 - 0.7e1_dp / 0.2e1_dp * t463 * t327  &
      + 0.63e2_dp / 0.4e1_dp * (t109) * (t110) * (t856) * (t287)  &
      * (t435) - 0.7e1_dp / 0.2e1_dp * (t109) * (t326) * (t1117)  &
      - t330 * d2Qrhondrho * t117 * t130 - t330 * dQrho * dQndrho * t117 * t130  &
      + (2._dp * t872 * t485) + (2._dp * t1202 * t368) + (2._dp * t875  &
      * t1215 * t367) + 0.2e1_dp * t330 * t336 * t342 * (-t471 * t24 * t123  &
      * t127 * t7 / 0.3e1_dp - t910 * t911 * t7 * t478 / 0.6e1_dp - t472 * t905  &
      * t126 - t355 * t905 * t479 / 0.2e1_dp + t472 * t363 / 0.2e1_dp - t355  &
      * t345 * t923 * t361 * t478 / 0.4e1_dp + t355 * t345 * t356 * (t1100  &
      * t71 - t930 * t476 - t1242 * t359 + 2._dp * t933 * t657 - t358  &
      * t71 * t678) / 0.2e1_dp)
      t1263 = t489 * t136
      t1286 = d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t1316 = -0.4e1_dp / 0.3e1_dp * t425 * t225 * t391 - 0.2e1_dp / 0.3e1_dp * t388  &
      * t429 * t390 + 0.2e1_dp / 0.3e1_dp * t388 * t1003 * t685 + 0.2e1_dp / 0.3e1_dp  &
      * t388 * t1007 * t7 * t435 - 0.4e1_dp * t425 * t165 * t150 - (2._dp  &
      * t394 * t505) + (2._dp * t394 * t508) + (2._dp * t394 * t511) +  &
      (2._dp * t501 * t398) + t147 * t15 * t647 * t149 - t147 * t397 * t507
      t1347 = -t147 * t397 * t510 - 2._dp * t501 * t402 - t147 * t504 * t401  &
      + 2._dp * t1023 * t933 * t84 * t218 * t261 + t1023 * t358 * t1024 * t435  &
      - t147 * t148 * t400 * t678 - 2._dp * t501 * t406 - t147 * t504 * t405  &
      + t1023 * t358 * t288 * t261 + 2._dp * t1023 * t80 * t1079 - t147  &
      * t148 * t404 * t1117
      t1352 = 0.1e1_dp / t240
      t1358 = t1059 * t1060 * t245
      t1362 = t414 * t415 * t261 * t83
      t1365 = t414 * t416 * t435
      t1369 = (2._dp * t77 * t1079) - (t77 * t94 * t1117) + f12 * t1143  &
      * t95 - (2._dp * t304 * t449) - (2._dp * t447 * t307) + (6._dp  &
      * t92 * C * t384 * t435) - (2._dp * t92 * t306 * t1117) + (-t1258  &
      * t134 * t139 - 0.7e1_dp / 0.2e1_dp * t952 * t491 - 0.7e1_dp / 0.2e1_dp * t1263  &
      * t377 - 0.35e2_dp / 0.4e1_dp * t374 * t956 * t1165 - 0.7e1_dp / 0.2e1_dp *  &
      t374 * t376 * (t1117)) * E * t143 - (3._dp * t381 * t497) - (3._dp  &
      * t495 * t384) + (12._dp * t142 * t968 * t287 * t435) - (3._dp  &
      * t142 * t383 * t1117) + (t145 * (t1286 + (t1316 + t1347) *  &
      t409 * t411 * t418 - 2._dp * t408 * t1352 * t411 * t418 - t1056 * t1358  &
      + t412 * t1362 + t412 * t1365))
      e_ndrho_rho = e_ndrho_rho +(  -0.4e1_dp / 0.3e1_dp * t269 * t520 - t76 * t1369 * Clda ) * sx0
      t1372 = t435 ** 2
      t1382 = t707 * t66
      t1388 = t566 * t716
      t1392 = t188 * t735
      t1395 = 2._dp * t887 * t878 + 4._dp * t426 * t430 - 4._dp * t501 * t433 + t79 *  &
      t70 * t1382 - 2._dp * t147 * t504 * t432 + 2._dp * t147 * t148 * t1388 -  &
      t147 * t148 * t1392
      t1412 = 2._dp * t69 * t89 + 4._dp * t426 * t444 + t79 * t70 * (F2 * t707 *  &
      t66 - 2._dp * t440 * t432 + 2._dp * t86 * t1388 - t86 * t1392)
      t1455 = dQndrho ** 2
      t1465 = t478 ** 2
      t1510 = d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
      t1547 = 2._dp * t887 * t150 + 4._dp * t501 * t505 - 4._dp * t501 * t508 - 4._dp * t501  &
      * t511 + t147 * t15 * t707 * t149 - 2._dp * t147 * t504 * t507 - 2._dp  &
      * t147 * t504 * t510 + 2._dp * t147 * t148 * t1018 * t716 + 2._dp * t1023  &
      * t358 * t94 * t261 * t435 - t147 * t148 * t400 * t735 + 2._dp * t147  &
      * t148 * t983 * t1372 - t147 * t148 * t404 * t1395
      t1561 = (2._dp * t77 * t143 * t1372) - (t77 * t94 * t1395) + (f12  &
      * t1412 * t95) - (4._dp * t447 * t449) + (6._dp * t92 * t824  &
      * t1372) - (2._dp * t92 * t306 * t1395) + (-((t98 * (t100 *  &
      t1412 * t83 + 2._dp * t100 * t446 * t435 + t100 * t91 * t1395 + 2._dp * t103  &
      * t1372 + 2._dp * t103 * t83 * t1395 + 6._dp * t105 * t83 * t1372 + 3._dp *  &
      t105 * t93 * t1395) * t114) - (7._dp * t463 * t465) + 0.63e2_dp / 0.4e1_dp  &
      * (t109) * (t857) * (t1372) - 0.7e1_dp / 0.2e1_dp * (t109)  &
      * (t326) * (t1395) - t330 * d2Qndrhondrho * t117 * t130 - t330  &
      * t1455 * t117 * t130 + (4._dp * t1202 * t485) + (2._dp * t875 *  &
      t1215 * t484) + 0.2e1_dp * t330 * t336 * t342 * (t472 * t480 - t355  &
      * t345 * t923 * t1465 / 0.4e1_dp + t355 * t345 * t356 * (t1382 *  &
      t71 - 2._dp * t1242 * t476 + 2._dp * t933 * t71 * t716 - t358 * t71 * t735)  &
      / 0.2e1_dp)) * t134 * t139 - (7._dp * t1263 * t491) - 0.35e2_dp / 0.4e1_dp  &
      * (t374) * (t956) * (t1372) - 0.7e1_dp / 0.2e1_dp * (t374)  &
      * (t376) * (t1395)) * E * (t143) - (6._dp * t495  &
      * t497) + (12._dp * t142 * t968 * t1372) - (3._dp * t142 * t383 *  &
      t1395) + (t145 * (t1510 + t1547 * t409 * t411 * t418 - 2._dp * t513  &
      * t1352 * t411 * t418 - t514 * t1055 * t1358 + t515 * t1362 + t515  &
      * t1365))
      e_ndrho_ndrho = e_ndrho_ndrho + ( -t76 * t1561 * Clda ) * sx0
    END IF

  END SUBROUTINE xwpbe_lda_calc_0

! *****************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx0 scaling factor
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega=0 using a taylor
!>      expansion for the parameter G.
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc_01(e_0, e_rho, e_ndrho, e_rho_rho, e_ndrho_rho, &
                              e_ndrho_ndrho, rho, ndrho, sscale, sx0, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho, &
                                                e_rho_rho, e_ndrho_rho, &
                                                e_ndrho_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, sscale, sx0
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: d2Qndrhondrho, d2Qrhondrho, d2Qrhorho, dQndrho, dQrho, &
      Q, t1, t10, t100, t101, t1019, t103, t104, t1051, t1056, t1062, t1066, &
      t1069, t107, t1073, t1076, t1080, t109, t1094, t1098, t11, t1101, t111, &
      t112, t113, t114, t115, t1154, t116, t118, t1191, t12, t1205, t122, &
      t124, t125, t126, t129, t13, t130, t131, t132, t135, t136, t139, t14, &
      t140, t141, t142, t146, t149, t15, t150, t151, t153, t154, t156, t157, &
      t158, t159, t16, t165, t166, t167, t168, t169, t173, t175, t18, t181, &
      t184, t185, t186, t188, t189, t19, t190, t191, t193, t194, t197, t199, &
      t2, t20, t202, t203, t206, t207
    REAL(KIND=dp) :: t208, t21, t211, t212, t213, t215, t216, t219, t22, &
      t220, t224, t227, t228, t229, t231, t232, t235, t238, t24, t240, t241, &
      t244, t247, t248, t25, t250, t251, t253, t254, t256, t257, t26, t260, &
      t263, t266, t267, t270, t272, t273, t276, t277, t28, t280, t283, t288, &
      t29, t293, t294, t297, t299, t3, t300, t301, t304, t305, t307, t308, &
      t31, t311, t314, t315, t317, t318, t319, t32, t321, t322, t323, t325, &
      t326, t327, t328, t329, t33, t330, t331, t332, t333, t334, t335, t34, &
      t340, t342, t343, t346, t347, t349, t350, t352, t357, t36, t360, t361, &
      t364, t366, t37, t371, t372, t375
    REAL(KIND=dp) :: t376, t379, t380, t383, t385, t388, t389, t39, t392, &
      t393, t395, t396, t398, t399, t4, t40, t401, t402, t403, t408, t410, &
      t412, t415, t418, t419, t42, t425, t428, t430, t434, t438, t44, t441, &
      t445, t454, t457, t458, t46, t465, t467, t47, t5, t504, t508, t509, &
      t515, t516, t52, t520, t53, t535, t54, t543, t552, t56, t567, t574, &
      t58, t59, t596, t6, t605, t61, t624, t63, t633, t640, t644, t65, t650, &
      t656, t66, t67, t674, t678, t68, t681, t69, t7, t70, t71, t714, t72, &
      t73, t74, t759, t76, t766, t77, t773, t777, t78, t79, t8, t80, t806, &
      t81, t811, t812, t820, t828, t83, t84, t847
    REAL(KIND=dp) :: t848, t849, t851, t852, t86, t865, t871, t874, t88, t89, &
      t9, t902, t906, t909, t92, t93, t94, t95, t97, t98, t989

    IF( order >= 0 ) THEN
      t1 = ndrho ** 2
      t2 = a1 * t1
      t3 = r2 ** 2
      t4 = 0.1e1_dp / t3
      t5 = t2 * t4
      t6 = pi ** 2
      t7 = r3 * t6
      t8 = t7 * rho
      t9 = t8 ** (0.1e1_dp / 0.3e1_dp)
      t10 = t9 ** 2
      t11 = 0.1e1_dp / t10
      t12 = rho ** 2
      t13 = 0.1e1_dp / t12
      t14 = t11 * t13
      t15 = sscale ** 2
      t16 = t14 * t15
      t18 = t1 ** 2
      t19 = a2 * t18
      t20 = t3 ** 2
      t21 = 0.1e1_dp / t20
      t22 = t19 * t21
      t24 = 0.1e1_dp / t9 / t8
      t25 = t12 ** 2
      t26 = 0.1e1_dp / t25
      t28 = t15 ** 2
      t29 = t24 * t26 * t28
      t31 = t5 * t16 + t22 * t29
      t32 = f94 * t31
      t33 = a3 * t18
      t34 = t33 * t21
      t36 = t18 * ndrho
      t37 = a4 * t36
      t39 = 0.1e1_dp / t20 / r2
      t40 = t37 * t39
      t42 = 0.1e1_dp / t10 / t8
      t44 = 0.1e1_dp / t25 / rho
      t46 = t28 * sscale
      t47 = t42 * t44 * t46
      t52 = 0.1e1_dp / t20 / t3
      t53 = a5 * t18 * t1 * t52
      t54 = r3 ** 2
      t56 = t6 ** 2
      t58 = 0.1e1_dp / t54 / t56
      t59 = t25 ** 2
      t61 = t28 * t15
      t63 = t58 / t59 * t61
      t65 = r1 + t34 * t29 + t40 * t47 + t53 * t63
      t66 = 0.1e1_dp / t65
      t67 = t66 * t1
      t68 = t32 * t67
      t69 = t4 * t11
      t70 = t13 * t15
      t71 = 0.1e1_dp / A
      t72 = t70 * t71
      t73 = t69 * t72
      Q = t68 * t73
      t74 = rho ** (0.1e1_dp / 0.3e1_dp)
      t76 = t74 * rho * f89
      t77 = B * f12
      t78 = t1 * t4
      t79 = t78 * t11
      t80 = t31 * t66
      t81 = t70 * t80
      t83 = t79 * t81 + DD
      t84 = 0.1e1_dp / t83
      t86 = F2 * t31
      t88 = F1 + t86 * t66
      t89 = t70 * t88
      t92 = f12 * (t79 * t89 + r1)
      t93 = t83 ** 2
      t94 = 0.1e1_dp / t93
      t95 = C * t94
      t97 = g2 * t1
      t98 = t97 * t4
      t100 = g3 * t18
      t101 = t100 * t21
      t103 = g1 + t98 * t16 + t101 * t29
      t104 = t70 * t103
      t107 = (t79 * t104 + r1) * E
      t109 = 0.1e1_dp / t93 / t83
      t111 = f12 * A
      t112 = exei(Q)
      t113 = t78 * t14
      t114 = t15 * t31
      t115 = t66 * t84
      t116 = t114 * t115
      t118 = LOG(t113 * t116)
      t122 = (t77 * t84 + t92 * t95 + t107 * t109 + t111 * (t112 + t118))  &
      * Clda
      e_0 = e_0 + ( -t76 * t122 ) * sx0
    END IF
    IF( order >=1 .OR. order == -1 ) THEN
      t124 = t4 * t42
      t125 = t2 * t124
      t126 = t70 * t7
      t129 = t12 * rho
      t130 = 0.1e1_dp / t129
      t131 = t11 * t130
      t132 = t131 * t15
      t135 = t54 * t56
      t136 = t135 * t12
      t139 = t21 / t9 / t136
      t140 = t19 * t139
      t141 = t26 * t28
      t142 = t141 * t7
      t146 = t24 * t44 * t28
      t149 = -0.2e1_dp / 0.3e1_dp * t125 * t126 - (2._dp * t5 * t132) - 0.4e1_dp /  &
      0.3e1_dp * t140 * t142 - (4._dp * t22 * t146)
      t150 = f94 * t149
      t151 = t150 * t67
      t153 = t65 ** 2
      t154 = 0.1e1_dp / t153
      t156 = t154 * t1 * t4
      t157 = t32 * t156
      t158 = t15 * t71
      t159 = t33 * t139
      t165 = 0.1e1_dp / t10 / t136
      t166 = t39 * t165
      t167 = t37 * t166
      t168 = t44 * t46
      t169 = t168 * t7
      t173 = 0.1e1_dp / t25 / t12
      t175 = t42 * t173 * t46
      t181 = t58 / t59 / rho * t61
      t184 = -0.4e1_dp / 0.3e1_dp * t159 * t142 - (4._dp * t34 * t146) - 0.5e1_dp  &
      / 0.3e1_dp * t167 * t169 - (5._dp * t40 * t175) - (8._dp * t53 * t181)
      t185 = t158 * t184
      t186 = t14 * t185
      t188 = t67 * t4
      t189 = t32 * t188
      t190 = t42 * t13
      t191 = t190 * t15
      t193 = t71 * r3 * t6
      t194 = t191 * t193
      t197 = t130 * t15
      t199 = t69 * t197 * t71
      dQrho = t151 * t73 - t157 * t186 - 0.2e1_dp / 0.3e1_dp * t189 * t194 - (2._dp  &
      * t68 * t199)
      t202 = a1 * ndrho
      t203 = t202 * t4
      t206 = t1 * ndrho
      t207 = a2 * t206
      t208 = t207 * t21
      t211 = 2._dp * t203 * t16 + 4._dp * t208 * t29
      t212 = f94 * t211
      t213 = t212 * t67
      t215 = a3 * t206
      t216 = t215 * t21
      t219 = a4 * t18
      t220 = t219 * t39
      t224 = a5 * t36 * t52
      t227 = 4._dp * t216 * t29 + 5._dp * t220 * t47 + 6._dp * t224 * t63
      t228 = t158 * t227
      t229 = t14 * t228
      t231 = t66 * ndrho
      t232 = t32 * t231
      dQndrho = t213 * t73 - t157 * t229 + 2._dp * t232 * t73
      t235 = t74 * f89
      t238 = t78 * t190
      t240 = t66 * r3 * t6
      t241 = t114 * t240
      t244 = t197 * t80
      t247 = t149 * t66
      t248 = t70 * t247
      t250 = t154 * t184
      t251 = t114 * t250
      t253 = -0.2e1_dp / 0.3e1_dp * t238 * t241 - (2._dp * t79 * t244) + (t79  &
      * t248) - t113 * t251
      t254 = t94 * t253
      t256 = t15 * t88
      t257 = t256 * t7
      t260 = t197 * t88
      t263 = F2 * t149
      t266 = t263 * t66 - t86 * t250
      t267 = t70 * t266
      t270 = f12 * (-0.2e1_dp / 0.3e1_dp * t238 * t257 - (2._dp * t79 * t260) +  &
      (t79 * t267))
      t272 = C * t109
      t273 = t272 * t253
      t276 = t15 * t103
      t277 = t276 * t7
      t280 = t197 * t103
      t283 = t97 * t124
      t288 = t100 * t139
      t293 = -0.2e1_dp / 0.3e1_dp * t283 * t126 - (2._dp * t98 * t132) - 0.4e1_dp  &
      / 0.3e1_dp * t288 * t142 - (4._dp * t101 * t146)
      t294 = t70 * t293
      t297 = (-0.2e1_dp / 0.3e1_dp * t238 * t277 - (2._dp * t79 * t280) + (t79  &
      * t294)) * E
      t299 = t93 ** 2
      t300 = 0.1e1_dp / t299
      t301 = t300 * t253
      t304 = dexeirho(Q,dQrho)
      t305 = t78 * t191
      t307 = t84 * r3 * t6
      t308 = t80 * t307
      t311 = t78 * t131
      t314 = t15 * t149
      t315 = t314 * t115
      t317 = t154 * t84
      t318 = t317 * t184
      t319 = t114 * t318
      t321 = t66 * t94
      t322 = t321 * t253
      t323 = t114 * t322
      t325 = -0.2e1_dp / 0.3e1_dp * t305 * t308 - (2._dp * t311 * t116) + t113  &
      * t315 - t113 * t319 - t113 * t323
      t326 = 0.1e1_dp / t1
      t327 = t325 * t326
      t328 = t3 * t10
      t329 = t327 * t328
      t330 = 0.1e1_dp / t15
      t331 = t12 * t330
      t332 = 0.1e1_dp / t31
      t333 = t332 * t65
      t334 = t333 * t83
      t335 = t331 * t334
      t340 = (-t77 * t254 + t270 * t95 - 2._dp * t92 * t273 + t297 * t109 - 3._dp  &
      * t107 * t301 + t111 * (t304 + t329 * t335)) * Clda
      e_rho = e_rho + ( -0.4e1_dp / 0.3e1_dp * t235 * t122 - t76 * t340 ) * sx0
      t342 = ndrho * t4
      t343 = t342 * t11
      t346 = t211 * t66
      t347 = t70 * t346
      t349 = t154 * t227
      t350 = t114 * t349
      t352 = 2._dp * t343 * t81 + t79 * t347 - t113 * t350
      t357 = F2 * t211
      t360 = t357 * t66 - t86 * t349
      t361 = t70 * t360
      t364 = f12 * (2._dp * t343 * t89 + t79 * t361)
      t366 = t272 * t352
      t371 = g2 * ndrho
      t372 = t371 * t4
      t375 = g3 * t206
      t376 = t375 * t21
      t379 = 2._dp * t372 * t16 + 4._dp * t376 * t29
      t380 = t70 * t379
      t383 = (2._dp * t343 * t104 + t79 * t380) * E
      t385 = t300 * t352
      t388 = dexeindrho(Q,dQndrho)
      t389 = t342 * t14
      t392 = t15 * t211
      t393 = t392 * t115
      t395 = t317 * t227
      t396 = t114 * t395
      t398 = t321 * t352
      t399 = t114 * t398
      t401 = 2._dp * t389 * t116 + t113 * t393 - t113 * t396 - t113 * t399
      t402 = t401 * t326
      t403 = t402 * t328
      t408 = (-t77 * t94 * t352 + t364 * t95 - 2._dp * t92 * t366 + t383 * t109  &
      - 3._dp * t107 * t385 + t111 * (t388 + t403 * t335)) * Clda
      e_ndrho = e_ndrho + ( -t76 * t408 ) * sx0
    END IF
    IF( order >= 2 .OR. order == -2 ) THEN
      t410 = t4 * t165
      t412 = t70 * t135
      t415 = t197 * t7
      t418 = t11 * t26
      t419 = t418 * t15
      t425 = t54 * r3 * t56 * t6 * t129
      t428 = t21 / t9 / t425
      t430 = t141 * t135
      t434 = t44 * t28 * t7
      t438 = t24 * t173 * t28
      t441 = 0.10e2_dp / 0.9e1_dp * t2 * t410 * t412 + 0.8e1_dp / 0.3e1_dp * t125 * t415  &
      + (6._dp * t5 * t419) + 0.28e2_dp / 0.9e1_dp * t19 * t428 * t430 + 0.32e2_dp  &
      / 0.3e1_dp * t140 * t434 + (20._dp * t22 * t438)
      t445 = t150 * t156
      t454 = 0.1e1_dp / t153 / t65
      t457 = t32 * t454 * t1 * t4
      t458 = t184 ** 2
      t465 = t32 * t154 * t78 * t42
      t467 = t184 * r3 * t6
      t504 = 0.28e2_dp / 0.9e1_dp * t33 * t428 * t430 + 0.32e2_dp / 0.3e1_dp * t159 *  &
      t434 + (20._dp * t34 * t438) + 0.40e2_dp / 0.9e1_dp * t37 * t39 / t10 /  &
      t425 * t168 * t135 + 0.50e2_dp / 0.3e1_dp * t167 * t173 * t46 * t7 + 0.30e2_dp  &
      * t40 * t42 / t25 / t129 * t46 + (72._dp * t53 * t58 / t59 /  &
      t12 * t61)
      t508 = t165 * t13
      t509 = t508 * t15
      t515 = t42 * t130
      t516 = t515 * t15
      t520 = t26 * t15
      d2Qrhorho = f94 * t441 * t67 * t73 - (2._dp * t445 * t186) - 0.4e1_dp / 0.3e1_dp  &
      * t150 * t188 * t194 - (4._dp * t151 * t199) + (2._dp * t457 * t14  &
      * t158 * t458) + 0.4e1_dp / 0.3e1_dp * t465 * t72 * t467 + (4._dp * t157  &
      * t131 * t185) - (t157 * t14 * t158 * t504) + 0.10e2_dp / 0.9e1_dp  &
      * t189 * t509 * t71 * t54 * t56 + 0.8e1_dp / 0.3e1_dp * t189 * t516 *  &
      t193 + 0.6e1_dp * t68 * t69 * t520 * t71
      t535 = -0.4e1_dp / 0.3e1_dp * t202 * t124 * t126 - (4._dp * t203 * t132)  &
      - 0.16e2_dp / 0.3e1_dp * t207 * t139 * t142 - (16._dp * t208 * t146)
      t543 = t212 * t156
      t552 = t32 * t154 * ndrho * t4
      t567 = -0.16e2_dp / 0.3e1_dp * t215 * t139 * t142 - (16._dp * t216 * t146)  &
      - 0.25e2_dp / 0.3e1_dp * t219 * t166 * t169 - (25._dp * t220 * t175) -  &
      (48._dp * t224 * t181)
      t574 = t7 * t227
      d2Qrhondrho = (f94 * t535 * t67 * t73) - t445 * t229 + (2._dp * t150 *  &
      t231 * t73) - (t543 * t186) + (2._dp * t457 * t16 * t71 * t184  &
      * t227) - (2._dp * t552 * t186) - (t157 * t14 * t158 * t567)  &
      - 0.2e1_dp / 0.3e1_dp * t212 * t188 * t194 + 0.2e1_dp / 0.3e1_dp * t465 * t72  &
      * t574 - 0.4e1_dp / 0.3e1_dp * t32 * (t231) * t4 * t194 - (2._dp * t213  &
      * t199) + (2._dp * t157 * t131 * t228) - (4._dp * t232 * t199)
      t596 = 2._dp * a1 * t4 * t16 + 12._dp * a2 * t1 * t21 * t29
      t605 = t227 ** 2
      t624 = 12._dp * a3 * t1 * t21 * t29 + 20._dp * a4 * t206 * t39 * t47 + 30._dp *  &
      a5 * t18 * t52 * t63
      d2Qndrhondrho = f94 * t596 * t67 * t73 - 2._dp * t543 * t229 + 4._dp * t212 * t231 *  &
      t73 + 2._dp * t457 * t14 * t158 * t605 - 4._dp * t552 * t229 - t157 * t14  &
      * t158 * t624 + 2._dp * t32 * t66 * t4 * t14 * t158
      t633 = t74 ** 2
      t640 = t253 ** 2
      t644 = t78 * t508
      t650 = t78 * t515
      t656 = t31 * t154
      t674 = t454 * t458
      t678 = t154 * t504
      t681 = 0.10e2_dp / 0.9e1_dp * t644 * t114 * t66 * t54 * t56 + 0.8e1_dp / 0.3e1_dp  &
      * t650 * t241 - 0.4e1_dp / 0.3e1_dp * t238 * t314 * t240 + 0.4e1_dp / 0.3e1_dp  &
      * t305 * t656 * t467 + (6._dp * t79 * t520 * t80) - (4._dp * t79  &
      * t197 * t247) + (4._dp * t311 * t251) + (t79) * t70 * t441  &
      * t66 - 0.2e1_dp * t113 * t314 * t250 + 0.2e1_dp * t113 * t114 * t674 -  &
      t113 * t114 * t678
      t714 = C * t300
      t759 = 0.1e1_dp / t299 / t83
      t766 = d2exeirhorho(Q,dQrho,d2Qrhorho)
      t773 = t656 * t84
      t777 = t80 * t94
      t806 = t454 * t84
      t811 = t78 * t16
      t812 = t94 * t184
      t820 = t66 * t109
      t828 = 0.8e1_dp / 0.3e1_dp * t78 * t516 * t308 - 0.4e1_dp / 0.3e1_dp * t305 * t247  &
      * t307 + 0.4e1_dp / 0.3e1_dp * t305 * t773 * t467 + 0.4e1_dp / 0.3e1_dp *  &
      t305 * t777 * t7 * t253 + 0.10e2_dp / 0.9e1_dp * t78 * t509 * t80 * t84  &
      * t54 * t56 + 0.6e1_dp * t78 * t418 * t116 - (4._dp * t311 * t315) +  &
      (4._dp * t311 * t319) + (4._dp * t311 * t323) + (t113 * t15 *  &
      t441 * t115) - (2._dp * t113 * t314 * t318) - (2._dp * t113 * t314  &
      * t322) + (2._dp * t113 * t114 * t806 * t458) + 0.2e1_dp * t811 * t656  &
      * t812 * t253 - (t113 * t114 * t317 * t504) + (2._dp * t113  &
      * t114 * t820 * t640) - (t113 * t114 * t321 * t681)
      t847 = t328 * t12
      t848 = t327 * t847
      t849 = t31 ** 2
      t851 = t330 / t849
      t852 = t65 * t83
      t865 = (2._dp * t77 * t109 * t640) - (t77 * t94 * t681) + f12 *  &
      (0.10e2_dp / 0.9e1_dp * t644 * t256 * t135 + 0.8e1_dp / 0.3e1_dp * t650 * t257  &
      - 0.4e1_dp / 0.3e1_dp * t238 * t15 * t266 * t7 + (6._dp * t79 * t520 *  &
      t88) - 0.4e1_dp * (t79) * t197 * t266 + (t79 * t70 * (F2 * t441  &
      * t66 - 2._dp * t263 * t250 + 2._dp * t86 * t674 - t86 * t678))) * t95  &
      - (4._dp * t270 * t273) + (6._dp * t92 * t714 * t640) - (2._dp *  &
      t92 * t272 * t681) + (0.10e2_dp / 0.9e1_dp * t644 * t276 * t135 + 0.8e1_dp  &
      / 0.3e1_dp * t650 * t277 - 0.4e1_dp / 0.3e1_dp * t238 * t15 * t293 * t7 + (6._dp  &
      * t79 * t520 * t103) - 0.4e1_dp * (t79) * t197 * t293 + (t79)  &
      * (t70) * (0.10e2_dp / 0.9e1_dp * t97 * t410 * t412 + 0.8e1_dp /  &
      0.3e1_dp * t283 * t415 + (6._dp * t98 * t419) + 0.28e2_dp / 0.9e1_dp * t100  &
      * t428 * t430 + 0.32e2_dp / 0.3e1_dp * t288 * t434 + (20._dp * t101 *  &
      t438))) * E * (t109) - (6._dp * t297 * t301) + (12._dp * t107  &
      * t759 * t640) - (3._dp * t107 * t300 * t681) + t111 * (t766 + t828  &
      * t326 * t328 * t335 + 0.2e1_dp / 0.3e1_dp * t327 * t3 / t9 * t12 * t330  &
      * t332 * t65 * t83 * r3 * t6 + 0.2e1_dp * t329 * rho * t330 * t334  &
      - t848 * t851 * t852 * t149 + t329 * t331 * t332 * t184 * t83 + t329  &
      * t331 * t333 * t253)
      e_rho_rho = e_rho_rho + ( -0.4e1_dp / 0.9e1_dp / t633 * f89 * t122 - 0.8e1_dp / 0.3e1_dp * t235 * t340  &
      - t76 * t865 * Clda ) * sx0
      t871 = t109 * t253 * t352
      t874 = t342 * t190
      t902 = t454 * t184 * t227
      t906 = t154 * t567
      t909 = -0.4e1_dp / 0.3e1_dp * t874 * t241 - 0.2e1_dp / 0.3e1_dp * t238 * t392 *  &
      t240 + 0.2e1_dp / 0.3e1_dp * t305 * t656 * t574 - (4._dp * t343 * t244)  &
      - (2._dp * t79 * t197 * t346) + (2._dp * t311 * t350) + (2._dp *  &
      t343 * t248) + (t79 * t70 * t535 * t66) - t113 * t314 * t349 -  &
      (2._dp * t389 * t251) - t113 * t392 * t250 + 0.2e1_dp * t113 * t114  &
      * t902 - t113 * t114 * t906
      t989 = d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t1019 = -0.4e1_dp / 0.3e1_dp * t342 * t191 * t308 - 0.2e1_dp / 0.3e1_dp * t305  &
      * t346 * t307 + 0.2e1_dp / 0.3e1_dp * t305 * t773 * t574 + 0.2e1_dp / 0.3e1_dp  &
      * t305 * t777 * t7 * t352 - 0.4e1_dp * t342 * t131 * t116 - (2._dp *  &
      t311 * t393) + (2._dp * t311 * t396) + (2._dp * t311 * t399) + (2._dp  &
      * t389 * t315) + t113 * t15 * t535 * t115 - t113 * t314 * t395
      t1051 = -t113 * t314 * t398 - 2._dp * t389 * t319 - t113 * t392 * t318  &
      + 2._dp * t811 * t31 * t454 * t84 * t184 * t227 + t811 * t656 * t812 *  &
      t352 - t113 * t114 * t317 * t567 - 2._dp * t389 * t323 - t113 * t392  &
      * t322 + t811 * t656 * t254 * t227 + 2._dp * t811 * t80 * t871 - t113  &
      * t114 * t321 * t909
      t1056 = 0.1e1_dp / t206
      t1062 = t851 * t852 * t211
      t1066 = t331 * t332 * t227 * t83
      t1069 = t331 * t333 * t352
      t1073 = (2._dp * t77 * t871) - (t77 * t94 * t909) + f12 * (-0.4e1_dp  &
      / 0.3e1_dp * t874 * t257 - 0.2e1_dp / 0.3e1_dp * t238 * t15 * t360 * t7  &
      - (4._dp * t343 * t260) - 0.2e1_dp * t79 * t197 * t360 + (2._dp * t343  &
      * t267) + t79 * t70 * (F2 * t535 * t66 - t263 * t349 - t357  &
      * t250 + 2._dp * t86 * t902 - t86 * t906)) * t95 - (2._dp * t270 * t366)  &
      - (2._dp * t364 * t273) + (6._dp * t92 * C * t301 * t352) - (2._dp  &
      * t92 * t272 * t909) + (-0.4e1_dp / 0.3e1_dp * t874 * t277 - 0.2e1_dp /  &
      0.3e1_dp * t238 * t15 * t379 * t7 - (4._dp * t343 * t280) - 0.2e1_dp *  &
      t79 * t197 * t379 + (2._dp * t343 * t294) + t79 * t70 * (-0.4e1_dp /  &
      0.3e1_dp * t371 * t124 * t126 - (4._dp * t372 * t132) - 0.16e2_dp / 0.3e1_dp  &
      * t375 * t139 * t142 - (16._dp * t376 * t146))) * E * t109 - (3._dp  &
      * t297 * t385) - (3._dp * t383 * t301) + (12._dp * t107 * t759  &
      * t253 * t352) - (3._dp * t107 * t300 * t909) + (t111 * (t989  &
      + (t1019 + t1051) * t326 * t328 * t335 - 2._dp * t325 * t1056 * t328  &
      * t335 - t848 * t1062 + t329 * t1066 + t329 * t1069))
      e_ndrho_rho = e_ndrho_rho + ( -0.4e1_dp / 0.3e1_dp * t235 * t408 - t76 * t1073 * Clda ) * sx0
      t1076 = t352 ** 2
      t1080 = t69 * t13
      t1094 = t454 * t605
      t1098 = t154 * t624
      t1101 = 2._dp * t1080 * t114 * t66 + 4._dp * t343 * t347 - 4._dp * t389 * t350  &
      + t79 * t70 * t596 * t66 - 2._dp * t113 * t392 * t349 + 2._dp * t113 * t114  &
      * t1094 - t113 * t114 * t1098
      t1154 = d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
      t1191 = 2._dp * t1080 * t116 + 4._dp * t389 * t393 - 4._dp * t389 * t396 - 4._dp *  &
      t389 * t399 + t113 * t15 * t596 * t115 - 2._dp * t113 * t392 * t395 -  &
      2._dp * t113 * t392 * t398 + 2._dp * t113 * t114 * t806 * t605 + 2._dp * t811  &
      * t656 * t94 * t227 * t352 - t113 * t114 * t317 * t624 + 2._dp * t113  &
      * t114 * t820 * t1076 - t113 * t114 * t321 * t1101
      t1205 = 2._dp * t77 * t109 * t1076 - t77 * t94 * t1101 + f12 * (2._dp * t69  &
      * t89 + 4._dp * t343 * t361 + t79 * t70 * (F2 * t596 * t66 - 2._dp * t357  &
      * t349 + 2._dp * t86 * t1094 - t86 * t1098)) * t95 - 4._dp * t364 * t366  &
      + 6._dp * t92 * t714 * t1076 - 2._dp * t92 * t272 * t1101 + (2._dp * t69 * t104  &
      + 4._dp * t343 * t380 + t79 * t70 * (2._dp * g2 * t4 * t16 + 12._dp * g3 * t1  &
      * t21 * t29)) * E * t109 - 6._dp * t383 * t385 + 12._dp * t107 * t759 *  &
      t1076 - 3._dp * t107 * t300 * t1101 + t111 * (t1154 + t1191 * t326 * t328  &
      * t335 - 2._dp * t401 * t1056 * t328 * t335 - t402 * t847 * t1062  &
      + t403 * t1066 + t403 * t1069)
      e_ndrho_ndrho = e_ndrho_ndrho + ( -t76 * t1205 * Clda ) * sx0
    END IF

  END SUBROUTINE xwpbe_lda_calc_01

! *****************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega scaling factor
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the exact functional for omega!=0.
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc_1(e_0, e_rho, e_ndrho, e_rho_rho, e_ndrho_rho, &
                              e_ndrho_ndrho, rho, ndrho, omega, sscale, sx, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho, &
                                                e_rho_rho, e_ndrho_rho, &
                                                e_ndrho_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, omega, sscale, sx
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: d2Qndrhondrho, d2Qrhondrho, d2Qrhorho, dQndrho, dQrho, &
      Q, t1, t10, t1001, t1009, t1011, t102, t1024, t1035, t104, t105, t1052, &
      t1066, t1068, t1069, t1071, t1072, t1077, t1078, t1079, t108, t1081, &
      t1082, t1085, t1086, t1087, t109, t1090, t1091, t1094, t1097, t1098, &
      t11, t110, t1100, t111, t1113, t1115, t1119, t1123, t1126, t1127, &
      t1129, t113, t1132, t1134, t1135, t1136, t1139, t115, t116, t1168, &
      t1169, t117, t1171, t1173, t1176, t1178, t118, t1181, t1185, t119, t12, &
      t120, t1208, t1209, t121, t1214, t1220, t123, t1236, t1237, t124, &
      t1240, t1242, t1243, t125, t1257, t1258, t126
    REAL(KIND=dp) :: t1264, t1265, t127, t128, t1281, t1282, t1285, t1287, &
      t129, t1294, t1297, t13, t1301, t1304, t1305, t1306, t1309, t131, &
      t1310, t1317, t1318, t132, t1324, t1325, t1326, t133, t1355, t136, &
      t1379, t1381, t1382, t1392, t14, t140, t1402, t1409, t141, t142, t1433, &
      t1437, t144, t1442, t1445, t145, t1456, t1457, t147, t1473, t1478, &
      t1479, t148, t149, t1491, t1492, t1498, t15, t150, t1501, t151, t152, &
      t1520, t1524, t153, t1531, t1535, t154, t1543, t1546, t1547, t1548, &
      t155, t1551, t1555, t1559, t156, t1562, t1563, t1564, t1568, t157, &
      t1572, t1573, t1579, t158, t1580, t159, t1592
    REAL(KIND=dp) :: t1599, t16, t1602, t1609, t1615, t162, t1624, t1625, &
      t163, t1633, t1649, t1658, t166, t1663, t167, t1677, t168, t1680, &
      t1681, t1682, t1686, t169, t1693, t1697, t17, t170, t1709, t1710, &
      t1714, t1722, t1723, t1726, t1729, t1731, t174, t1740, t176, t1767, &
      t1768, t177, t1770, t1781, t1786, t1791, t18, t180, t1802, t1804, &
      t1805, t181, t1816, t1827, t1835, t1839, t1846, t185, t1853, t1854, &
      t186, t1860, t1868, t1880, t1881, t1885, t189, t1898, t19, t190, t1906, &
      t1907, t192, t1923, t193, t1932, t1938, t194, t195, t1951, t1956, &
      t1961, t1967, t197, t1977, t198, t1983, t1984, t1988
    REAL(KIND=dp) :: t199, t1995, t2, t200, t202, t2028, t205, t2059, t2072, &
      t2073, t209, t2099, t21, t210, t2108, t213, t214, t2142, t2144, t2150, &
      t2173, t2184, t219, t2197, t22, t221, t222, t2224, t223, t224, t2240, &
      t2245, t2254, t2258, t2267, t2269, t227, t2271, t2274, t228, t2280, &
      t2282, t2285, t2291, t2297, t23, t2305, t2311, t2316, t2323, t2329, &
      t233, t2348, t2363, t237, t2371, t2379, t239, t24, t2401, t2405, t241, &
      t2410, t243, t2437, t2442, t2449, t2452, t2455, t2457, t246, t247, &
      t2473, t2484, t2499, t25, t250, t251, t2512, t252, t2529, t253, t2543, &
      t255, t256, t257, t2573, t258, t263
    REAL(KIND=dp) :: t264, t265, t267, t2688, t27, t270, t2707, t271, t2715, &
      t273, t275, t276, t2764, t2774, t28, t280, t2808, t2810, t2838, t284, &
      t2841, t2844, t2846, t285, t286, t2875, t288, t2880, t2884, t289, t29, &
      t290, t2927, t293, t294, t295, t2963, t2971, t2975, t2979, t298, t2994, &
      t3, t3001, t303, t3033, t305, t307, t31, t310, t311, t312, t313, t3139, &
      t315, t316, t3167, t318, t319, t32, t321, t326, t329, t330, t332, t333, &
      t335, t336, t338, t339, t34, t340, t342, t343, t345, t346, t347, t348, &
      t349, t35, t351, t352, t353, t354, t357, t358, t36, t361, t362, t363, &
      t364, t368, t371, t372
    REAL(KIND=dp) :: t373, t374, t375, t376, t377, t378, t38, t383, t384, &
      t385, t386, t39, t390, t392, t398, t4, t401, t402, t403, t404, t406, &
      t409, t41, t411, t412, t415, t416, t419, t42, t420, t421, t424, t425, &
      t426, t428, t429, t432, t433, t437, t44, t440, t441, t442, t444, t446, &
      t449, t452, t453, t454, t457, t46, t461, t463, t465, t468, t469, t472, &
      t475, t478, t479, t48, t481, t486, t49, t493, t495, t498, t499, t5, &
      t502, t503, t504, t507, t508, t509, t510, t511, t513, t514, t515, t517, &
      t518, t522, t525, t529, t530, t531, t532, t533, t534, t535, t537, t538, &
      t539, t54, t540, t542, t55, t550
    REAL(KIND=dp) :: t551, t552, t553, t554, t556, t557, t558, t56, t562, &
      t563, t567, t569, t571, t572, t575, t576, t577, t578, t58, t581, t582, &
      t583, t587, t588, t589, t592, t596, t597, t6, t60, t602, t603, t608, &
      t61, t613, t617, t620, t621, t624, t626, t629, t63, t631, t634, t638, &
      t639, t644, t645, t65, t653, t655, t656, t657, t658, t659, t663, t666, &
      t669, t67, t673, t679, t68, t681, t685, t689, t69, t690, t691, t697, &
      t698, t7, t70, t701, t71, t711, t713, t717, t718, t72, t721, t728, t73, &
      t735, t736, t739, t74, t740, t746, t747, t748, t75, t752, t753, t754, &
      t755, t759, t761, t763, t765, t769
    REAL(KIND=dp) :: t77, t771, t772, t779, t78, t780, t781, t783, t784, &
      t787, t791, t792, t799, t8, t80, t800, t801, t803, t804, t808, t81, &
      t810, t812, t815, t816, t817, t82, t820, t823, t826, t83, t832, t834, &
      t836, t84, t842, t845, t846, t848, t849, t85, t851, t865, t867, t87, &
      t870, t871, t873, t874, t876, t877, t88, t880, t884, t885, t888, t889, &
      t891, t892, t893, t897, t898, t9, t90, t902, t904, t907, t908, t909, &
      t91, t910, t916, t918, t919, t92, t93, t930, t932, t933, t937, t94, &
      t942, t945, t95, t951, t954, t958, t96, t962, t965, t968, t97, t971, &
      t972, t979, t982, t988, t99

    IF( order >= 0 ) THEN
      t1 = ndrho ** 2
      t2 = r2 ** 2
      t3 = 0.1e1_dp / t2
      t4 = t1 * t3
      t5 = pi ** 2
      t6 = r3 * t5
      t7 = t6 * rho
      t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
      t9 = t8 ** 2
      t10 = 0.1e1_dp / t9
      t11 = t4 * t10
      t12 = rho ** 2
      t13 = 0.1e1_dp / t12
      t14 = sscale ** 2
      t15 = t13 * t14
      t16 = a1 * t1
      t17 = t16 * t3
      t18 = t10 * t13
      t19 = t18 * t14
      t21 = t1 ** 2
      t22 = a2 * t21
      t23 = t2 ** 2
      t24 = 0.1e1_dp / t23
      t25 = t22 * t24
      t27 = 0.1e1_dp / t8 / t7
      t28 = t12 ** 2
      t29 = 0.1e1_dp / t28
      t31 = t14 ** 2
      t32 = t27 * t29 * t31
      t34 = t17 * t19 + t25 * t32
      t35 = a3 * t21
      t36 = t35 * t24
      t38 = t21 * ndrho
      t39 = a4 * t38
      t41 = 0.1e1_dp / t23 / r2
      t42 = t39 * t41
      t44 = 0.1e1_dp / t9 / t7
      t46 = 0.1e1_dp / t28 / rho
      t48 = t31 * sscale
      t49 = t44 * t46 * t48
      t54 = 0.1e1_dp / t23 / t2
      t55 = a5 * t21 * t1 * t54
      t56 = r3 ** 2
      t58 = t5 ** 2
      t60 = 0.1e1_dp / t56 / t58
      t61 = t28 ** 2
      t63 = t31 * t14
      t65 = t60 / t61 * t63
      t67 = r1 + t36 * t32 + t42 * t49 + t55 * t65
      t68 = 0.1e1_dp / t67
      t69 = t34 * t68
      t70 = t15 * t69
      t71 = t11 * t70
      t72 = omega ** 2
      t73 = beta * t72
      t74 = t73 * t10
      t75 = t71 + t74
      t77 = 0.1e1_dp / A
      Q = f94 * t75 * t77
      t78 = rho ** (0.1e1_dp / 0.3e1_dp)
      t80 = t78 * rho * f89
      t81 = B * f12
      t82 = t71 + DD
      t83 = 0.1e1_dp / t82
      t84 = t81 * t83
      t85 = F2 * t34
      t87 = F1 + t85 * t68
      t88 = t15 * t87
      t90 = t11 * t88 + r1
      t91 = f12 * t90
      t92 = t82 ** 2
      t93 = 0.1e1_dp / t92
      t94 = C * t93
      t95 = t91 * t94
      t96 = f34 * pi
      t97 = rootpi
      t99 = r6 * C
      t102 = r4 * B
      t104 = r8 * A
      t105 = t92 * t82
      t108 = t97 * (r15 * E + t99 * t90 * t82 + t102 * t92 + t104 * t105)
      t109 = 0.1e1_dp / r16
      t110 = SQRT(t82)
      t111 = t110 * t105
      t113 = t109 / t111
      t115 = SQRT(A)
      t116 = f94 * t34
      t117 = t68 * t1
      t118 = t116 * t117
      t119 = t3 * t10
      t120 = t15 * t77
      t121 = t119 * t120
      t123 = EXP(t118 * t121)
      t124 = t115 * t123
      t125 = f32 * ndrho
      t126 = 0.1e1_dp / r2
      t127 = t125 * t126
      t128 = 0.1e1_dp / t8
      t129 = 0.1e1_dp / rho
      t131 = t69 * t77
      t132 = SQRT(t131)
      t133 = sscale * t132
      t136 = erfc(t127 * t128 * t129 * t133)
      t140 = 0.1e1_dp / f1516
      t141 = (t96 + t108 * t113 - t96 * t124 * t136) * t140
      t142 = 0.1e1_dp / t97
      t144 = 0.1e1_dp / E
      t145 = t142 * t111 * t144
      t147 = -t141 * t145 + r1
      t148 = t147 * E
      t149 = 0.1e1_dp / t105
      t150 = t148 * t149
      t151 = f158 * E
      t152 = t147 * t83
      t153 = t72 * t10
      t154 = t71 + DD + t153
      t155 = t154 ** 2
      t156 = t155 ** 2
      t157 = t156 * t154
      t158 = SQRT(t157)
      t159 = 0.1e1_dp / t158
      t162 = SQRT(t154)
      t163 = 0.1e1_dp / t162
      t166 = f68 * C
      t167 = t90 * t83
      t168 = t155 * t154
      t169 = SQRT(t168)
      t170 = 0.1e1_dp / t169
      t174 = (-t151 * t152 * t159 - t81 * t83 * t163 - t166 * t167 * t170)  &
      * omega
      t176 = f52 * E
      t177 = t147 * t93
      t180 = f12 * C
      t181 = t90 * t93
      t185 = t72 * omega
      t186 = (-t176 * t177 * t159 - t180 * t181 * t170) * t185
      t189 = 0.1e1_dp / r3 / t5
      t190 = t189 * t129
      t192 = t72 ** 2
      t193 = t192 * omega
      t194 = t159 * t193
      t195 = t194 * t44
      t197 = f12 * A
      t198 = exei(Q)
      t199 = t71 + DD + t74
      t200 = 0.1e1_dp / t199
      t202 = LOG(t75 * t200)
      t205 = SQRT(t199)
      t209 = t115 * f34
      t210 = exer(Q)
      t213 = (t197 * t97 / t205 - t209 * t210) * alpha1
      t214 = omega * t128
      t219 = (t197 * t200 - f98 * t198) * alpha2
      t221 = A * f14
      t222 = t199 ** 2
      t223 = t222 * t199
      t224 = SQRT(t223)
      t227 = SQRT(t75)
      t228 = 0.1e1_dp / t227
      t233 = 0.1e1_dp / t115
      t237 = (t97 * (t221 / t224 - f98 * t228) + f2716 * t210 * t233) * alpha3  &
      * t185
      t239 = 0.1e1_dp / t75
      t241 = 0.1e1_dp / t222
      t243 = f8132 * t77
      t246 = (-f98 * t239 + t197 * t241 + t243 * t198) * alpha4
      t247 = t192 * t27
      t250 = t75 ** 2
      t251 = t250 * t75
      t252 = SQRT(t251)
      t253 = 0.1e1_dp / t252
      t255 = f38 * A
      t256 = t222 ** 2
      t257 = t256 * t199
      t258 = SQRT(t257)
      t263 = A ** 2
      t264 = t263 * A
      t265 = SQRT(t264)
      t267 = f24364 / t265
      t270 = (t97 * (t243 * t228 - f916 * t253 + t255 / t258) - t267 * t210)  &
      * alpha5
      t271 = t193 * t44
      t273 = 0.1e1_dp / t223
      t275 = 0.1e1_dp / t250
      t276 = f98 * t275
      t280 = f729128 / t263
      t284 = t192 * t72
      t285 = (A * t273 - t276 + t243 * r1 * t239 - t280 * t198) * alpha6  &
      * t284
      t286 = t60 * t13
      t288 = f1516 * A
      t289 = t256 * t223
      t290 = SQRT(t289)
      t293 = t250 ** 2
      t294 = t293 * t75
      t295 = SQRT(t294)
      t298 = f8164 * t77
      t303 = t263 ** 2
      t305 = SQRT(t303 * A)
      t307 = f2187256 / t305
      t310 = (t97 * (t288 / t290 - f2732 / t295 + t298 * t253 - t280 * t228)  &
      + t307 * t210) * alpha7
      t311 = t192 * t185
      t312 = t56 * t58
      t313 = t312 * t12
      t315 = 0.1e1_dp / t8 / t313
      t316 = t311 * t315
      t318 = r3 * A
      t319 = 0.1e1_dp / t256
      t321 = 0.1e1_dp / t251
      t326 = f6561512 / t264
      t329 = (t318 * t319 - f94 * t321 + t243 * t275 - t280 * t239 + t326  &
      * t198) * alpha8
      t330 = t192 ** 2
      t332 = 0.1e1_dp / t9 / t313
      t333 = t330 * t332
      t335 = t84 + t95 + t150 + t174 * t128 + t186 * t190 - t150 * t195 +  &
      t197 * (t198 + t202) + t213 * t214 + t219 * t153 + t237 * t190 +  &
      t246 * t247 + t270 * t271 + t285 * t286 + t310 * t316 + t329 * t333
      t336 = t335 * Clda
      e_0 = e_0  +(  -t80 * t336 ) * sx
    END IF
    IF( order >= 1 .OR. order == -1 ) THEN
      t338 = t44 * t13
      t339 = t4 * t338
      t340 = t14 * t34
      t342 = t68 * r3 * t5
      t343 = t340 * t342
      t345 = 0.2e1_dp / 0.3e1_dp * t339 * t343
      t346 = t12 * rho
      t347 = 0.1e1_dp / t346
      t348 = t347 * t14
      t349 = t348 * t69
      t351 = 2._dp * t11 * t349
      t352 = t3 * t44
      t353 = t16 * t352
      t354 = t15 * t6
      t357 = t10 * t347
      t358 = t357 * t14
      t361 = t24 * t315
      t362 = t22 * t361
      t363 = t29 * t31
      t364 = t363 * t6
      t368 = t27 * t46 * t31
      t371 = -0.2e1_dp / 0.3e1_dp * t353 * t354 - (2._dp * t17 * t358) - 0.4e1_dp  &
      / 0.3e1_dp * t362 * t364 - (4._dp * t25 * t368)
      t372 = t371 * t68
      t373 = t15 * t372
      t374 = t11 * t373
      t375 = t4 * t18
      t376 = t67 ** 2
      t377 = 0.1e1_dp / t376
      t378 = t35 * t361
      t383 = t41 * t332
      t384 = t39 * t383
      t385 = t46 * t48
      t386 = t385 * t6
      t390 = 0.1e1_dp / t28 / t12
      t392 = t44 * t390 * t48
      t398 = t60 / t61 / rho * t63
      t401 = -0.4e1_dp / 0.3e1_dp * t378 * t364 - (4._dp * t36 * t368) - 0.5e1_dp  &
      / 0.3e1_dp * t384 * t386 - (5._dp * t42 * t392) - (8._dp * t55 * t398)
      t402 = t377 * t401
      t403 = t340 * t402
      t404 = t375 * t403
      t406 = t44 * r3 * t5
      t409 = -t345 - t351 + t374 - t404 - 0.2e1_dp / 0.3e1_dp * t73 * t406
      dQrho = f94 * t409 * t77
      t411 = ndrho * t3
      t412 = t411 * t10
      t415 = a1 * ndrho
      t416 = t415 * t3
      t419 = t1 * ndrho
      t420 = a2 * t419
      t421 = t420 * t24
      t424 = 2._dp * t416 * t19 + 4._dp * t421 * t32
      t425 = t424 * t68
      t426 = t15 * t425
      t428 = a3 * t419
      t429 = t428 * t24
      t432 = a4 * t21
      t433 = t432 * t41
      t437 = a5 * t38 * t54
      t440 = 4._dp * t429 * t32 + 5._dp * t433 * t49 + 6._dp * t437 * t65
      t441 = t377 * t440
      t442 = t340 * t441
      t444 = 2._dp * t412 * t70 + t11 * t426 - t375 * t442
      dQndrho = f94 * t444 * t77
      t446 = t78 * f89
      t449 = t60 * t347
      t452 = C * t149
      t453 = -t345 - t351 + t374 - t404
      t454 = t452 * t453
      t457 = t329 * t330
      t461 = t56 * r3 * t58 * t5 * t346
      t463 = 0.1e1_dp / t9 / t461
      t465 = t463 * r3 * t5
      t468 = t14 * t87
      t469 = t468 * t6
      t472 = t348 * t87
      t475 = F2 * t371
      t478 = t475 * t68 - t85 * t402
      t479 = t15 * t478
      t481 = -0.2e1_dp / 0.3e1_dp * t339 * t469 - (2._dp * t11 * t472) + (t11  &
      * t479)
      t486 = t82 * t453
      t493 = t97 * (t99 * t481 * t82 + t99 * t90 * t453 + 2._dp * t102 * t486  &
      + 3._dp * t104 * t92 * t453)
      t495 = t92 ** 2
      t498 = t109 / t110 / t495
      t499 = t498 * t453
      t502 = t96 * t115
      t503 = f94 * t371
      t504 = t503 * t117
      t507 = t377 * t1 * t3
      t508 = t116 * t507
      t509 = t14 * t77
      t510 = t509 * t401
      t511 = t18 * t510
      t513 = t117 * t3
      t514 = t116 * t513
      t515 = t338 * t14
      t517 = t77 * r3 * t5
      t518 = t515 * t517
      t522 = t119 * t348 * t77
      t525 = t504 * t121 - t508 * t511 - 0.2e1_dp / 0.3e1_dp * t514 * t518 - (2._dp  &
      * t118 * t522)
      t529 = rootpi
      t530 = 0.1e1_dp / t529
      t531 = t123 * t530
      t532 = f32 ** 2
      t533 = t532 * t1
      t534 = t533 * t119
      t535 = t15 * t131
      t537 = EXP(-t534 * t535)
      t538 = t126 * t27
      t539 = t125 * t538
      t540 = t129 * sscale
      t542 = t132 * r3 * t5
      t550 = t125 * t126 * t128
      t551 = 0.1e1_dp / t132
      t552 = t372 * t77
      t553 = t34 * t377
      t554 = t77 * t401
      t556 = t552 - t553 * t554
      t557 = t551 * t556
      t558 = t540 * t557
      t562 = t537 * (-t539 * t540 * t542 / 0.3e1_dp - t127 * t128 * t13 * t133  &
      + t550 * t558 / 0.2e1_dp)
      t563 = t531 * t562
      t567 = (t493 * t113 - 0.7e1_dp / 0.2e1_dp * t108 * t499 - (t502 * t525  &
      * t123 * t136) + (2._dp * t502 * t563)) * t140
      t569 = t141 * t142
      t571 = t110 * t92 * t144
      t572 = t571 * t453
      t575 = -t567 * t145 - 0.7e1_dp / 0.2e1_dp * t569 * t572
      t576 = t575 * E
      t577 = t576 * t149
      t578 = t189 * t13
      t581 = 0.1e1_dp / t158 / t157
      t582 = t149 * t581
      t583 = t148 * t582
      t587 = -t345 - t351 + t374 - t404 - 0.2e1_dp / 0.3e1_dp * t72 * t44 * t6
      t588 = t156 * t587
      t589 = t271 * t588
      t592 = t219 * t72
      t596 = 0.1e1_dp / t224 / t223
      t597 = t596 * t222
      t602 = 0.1e1_dp / t227 / t75
      t603 = f98 * t602
      t608 = dexerrho(Q,dQrho)
      t613 = (t97 * (-0.3e1_dp / 0.2e1_dp * t221 * t597 * t409 + t603 * t409 /  &
      0.2e1_dp) + f2716 * t608 * t233) * alpha3 * t185
      t617 = f12 * t481
      t620 = 0.1e1_dp / t495
      t621 = t620 * t453
      t624 = t213 * omega
      t626 = t27 * r3 * t5
      t629 = t246 * t192
      t631 = t315 * r3 * t5
      t634 = t602 * t409
      t638 = 0.1e1_dp / t252 / t251
      t639 = f916 * t638
      t644 = 0.1e1_dp / t258 / t257
      t645 = t644 * t256
      t653 = (t97 * (-t243 * t634 / 0.2e1_dp + 0.3e1_dp / 0.2e1_dp * t639 * t250 *  &
      t409 - 0.5e1_dp / 0.2e1_dp * t255 * t645 * t409) - t267 * t608) * alpha5
      t655 = -(2._dp * t285 * t449) - (2._dp * t91 * t454) - 0.8e1_dp / 0.3e1_dp  &
      * t457 * t465 + t577 - t186 * t578 + 0.5e1_dp / 0.2e1_dp * t583 * t589  &
      - 0.2e1_dp / 0.3e1_dp * t592 * t406 + t613 * t190 - t81 * t93 * t453 +  &
      t617 * t94 - t237 * t578 - (3._dp * t148 * t621) - t624 * t626 / 0.3e1_dp  &
      - 0.4e1_dp / 0.3e1_dp * t629 * t631 + t653 * t271
      t656 = t149 * t159
      t657 = t148 * t656
      t658 = t193 * t332
      t659 = t658 * t6
      t663 = t273 * t409
      t666 = dexeirho(Q,dQrho)
      t669 = (t276 * t409 - 2._dp * t197 * t663 + t243 * t666) * alpha4
      t673 = t97 / t205 / t199
      t679 = (-t197 * t673 * t409 / 0.2e1_dp - t209 * t608) * alpha1
      t681 = t241 * t409
      t685 = (-t197 * t681 - f98 * t666) * alpha2
      t689 = 0.1e1_dp / t290 / t289
      t690 = t256 * t222
      t691 = t689 * t690
      t697 = f2732 / t295 / t294
      t698 = t293 * t409
      t701 = t638 * t250
      t711 = (t97 * (-0.7e1_dp / 0.2e1_dp * t288 * t691 * t409 + 0.5e1_dp / 0.2e1_dp  &
      * t697 * t698 - 0.3e1_dp / 0.2e1_dp * t298 * t701 * t409 + t280 * t634 /  &
      0.2e1_dp) + t307 * t608) * alpha7
      t713 = 0.1e1_dp / t257
      t717 = 0.1e1_dp / t293
      t718 = f94 * t717
      t721 = t321 * t409
      t728 = (-4._dp * t318 * t713 * t409 + 3._dp * t718 * t409 - 2._dp * t243 * t721  &
      + t280 * t275 * t409 + t326 * t666) * alpha8
      t735 = t176 * t147
      t736 = t656 * t453
      t739 = t93 * t581
      t740 = t739 * t588
      t746 = t180 * t90
      t747 = t149 * t170
      t748 = t747 * t453
      t752 = 0.1e1_dp / t169 / t168
      t753 = t93 * t752
      t754 = t155 * t587
      t755 = t753 * t754
      t759 = (-t176 * t575 * t93 * t159 + (2._dp * t735 * t736) + 0.5e1_dp /  &
      0.2e1_dp * (t735) * (t740) - t180 * t481 * t93 * t170 + (2._dp  &
      * t746 * t748) + 0.3e1_dp / 0.2e1_dp * (t746) * (t755)) * t185
      t761 = t310 * t311
      t763 = 0.1e1_dp / t8 / t461
      t765 = t763 * r3 * t5
      t769 = t75 * t241
      t771 = t409 * t200 - t769 * t409
      t772 = t771 * t239
      t779 = t151 * t147
      t780 = t93 * t159
      t781 = t780 * t453
      t783 = t83 * t581
      t784 = t783 * t588
      t787 = t93 * t163
      t791 = 0.1e1_dp / t162 / t154
      t792 = t83 * t791
      t799 = t166 * t90
      t800 = t93 * t170
      t801 = t800 * t453
      t803 = t83 * t752
      t804 = t803 * t754
      t808 = (-t151 * t575 * t83 * t159 + t779 * t781 + 0.5e1_dp / 0.2e1_dp * t779  &
      * t784 + t81 * t787 * t453 + t81 * t792 * t587 / 0.2e1_dp - t166  &
      * t481 * t83 * t170 + t799 * t801 + 0.3e1_dp / 0.2e1_dp * t799 * t804) *  &
      omega
      t810 = t148 * t620
      t812 = t194 * t44 * t453
      t815 = t270 * t193
      t816 = t332 * r3
      t817 = t816 * t5
      t820 = A * t319
      t823 = f98 * t321
      t826 = r1 * t275
      t832 = (-3._dp * t820 * t409 + 2._dp * t823 * t409 - t243 * t826 * t409 - t280  &
      * t666) * alpha6 * t284
      t834 = 0.5e1_dp / 0.3e1_dp * t657 * t659 + t669 * t247 + t679 * t214 + t685  &
      * t153 - t577 * t195 + t711 * t316 + t728 * t333 - t174 * t626  &
      / 0.3e1_dp + t759 * t190 - 0.7e1_dp / 0.3e1_dp * t761 * t765 + t197 * (t666  &
      + t772 * t199) + t808 * t128 + (3._dp * t810 * t812) - 0.5e1_dp / 0.3e1_dp  &
      * t815 * t817 + t832 * t286
      t836 = (t655 + t834) * Clda
      e_rho = e_rho + ( -0.4e1_dp / 0.3e1_dp * t446 * t336 - t80 * t836 ) * sx
      t842 = F2 * t424
      t845 = t842 * t68 - t85 * t441
      t846 = t15 * t845
      t848 = 2._dp * t412 * t88 + t11 * t846
      t849 = f12 * t848
      t851 = t452 * t444
      t865 = t97 * (t99 * t848 * t82 + t99 * t90 * t444 + 2._dp * t102 * t82  &
      * t444 + 3._dp * t104 * t92 * t444)
      t867 = t498 * t444
      t870 = f94 * t424
      t871 = t870 * t117
      t873 = t509 * t440
      t874 = t18 * t873
      t876 = t68 * ndrho
      t877 = t116 * t876
      t880 = t871 * t121 - t508 * t874 + 2._dp * t877 * t121
      t884 = f32 * t126
      t885 = t884 * t128
      t888 = t425 * t77
      t889 = t77 * t440
      t891 = t888 - t553 * t889
      t892 = t551 * t891
      t893 = t540 * t892
      t897 = t537 * (t885 * t540 * t132 + t550 * t893 / 0.2e1_dp)
      t898 = t531 * t897
      t902 = (t865 * t113 - 0.7e1_dp / 0.2e1_dp * t108 * t867 - (t502 * t880  &
      * t123 * t136) + (2._dp * t502 * t898)) * t140
      t904 = t571 * t444
      t907 = -t902 * t145 - 0.7e1_dp / 0.2e1_dp * t569 * t904
      t908 = t907 * E
      t909 = t908 * t149
      t910 = t620 * t444
      t916 = t780 * t444
      t918 = t156 * t444
      t919 = t783 * t918
      t930 = t800 * t444
      t932 = t155 * t444
      t933 = t803 * t932
      t937 = (-t151 * t907 * t83 * t159 + t779 * t916 + 0.5e1_dp / 0.2e1_dp * t779  &
      * t919 + t81 * t787 * t444 + t81 * t792 * t444 / 0.2e1_dp - t166  &
      * t848 * t83 * t170 + t799 * t930 + 0.3e1_dp / 0.2e1_dp * t799 * t933) *  &
      omega
      t942 = t656 * t444
      t945 = t739 * t918
      t951 = t747 * t444
      t954 = t753 * t932
      t958 = (-t176 * t907 * t93 * t159 + (2._dp * t735 * t942) + 0.5e1_dp /  &
      0.2e1_dp * (t735) * (t945) - t180 * t848 * t93 * t170 + (2._dp  &
      * t746 * t951) + 0.3e1_dp / 0.2e1_dp * (t746) * (t954)) * t185
      t962 = t194 * t44 * t444
      t965 = t271 * t918
      t968 = dexeindrho(Q,dQndrho)
      t971 = t444 * t200 - t769 * t444
      t972 = t971 * t239
      t979 = dexerndrho(Q,dQndrho)
      t982 = (-t197 * t673 * t444 / 0.2e1_dp - t209 * t979) * alpha1
      t988 = (-t197 * t241 * t444 - f98 * t968) * alpha2
      t1001 = (t97 * (-0.3e1_dp / 0.2e1_dp * t221 * t597 * t444 + t603 * t444 /  &
      0.2e1_dp) + f2716 * t979 * t233) * alpha3 * t185
      t1009 = (t276 * t444 - 2._dp * t197 * t273 * t444 + t243 * t968) * alpha4
      t1011 = t602 * t444
      t1024 = (t97 * (-t243 * t1011 / 0.2e1_dp + 0.3e1_dp / 0.2e1_dp * t639 * t250  &
      * t444 - 0.5e1_dp / 0.2e1_dp * t255 * t645 * t444) - t267 * t979) * alpha5
      t1035 = (-3._dp * t820 * t444 + 2._dp * t823 * t444 - t243 * t826 * t444 -  &
      t280 * t968) * alpha6 * t284
      t1052 = (t97 * (-0.7e1_dp / 0.2e1_dp * t288 * t691 * t444 + 0.5e1_dp / 0.2e1_dp  &
      * t697 * t293 * t444 - 0.3e1_dp / 0.2e1_dp * t298 * t701 * t444 + t280  &
      * t1011 / 0.2e1_dp) + t307 * t979) * alpha7
      t1066 = (-4._dp * t318 * t713 * t444 + 3._dp * t718 * t444 - 2._dp * t243 * t321  &
      * t444 + t280 * t275 * t444 + t326 * t968) * alpha8
      t1068 = -t81 * t93 * t444 + t849 * t94 - (2._dp * t91 * t851) + t909  &
      - (3._dp * t148 * t910) + t937 * t128 + t958 * t190 - t909 * t195  &
      + (3._dp * t810 * t962) + 0.5e1_dp / 0.2e1_dp * t583 * t965 + t197 * (t968  &
      + t972 * t199) + t982 * t214 + t988 * t153 + t1001 * t190 + t1009  &
      * t247 + t1024 * t271 + t1035 * t286 + t1052 * t316 + t1066 *  &
      t333
      t1069 = t1068 * Clda
      e_ndrho = e_ndrho +(  -t80 * t1069) * sx
    END IF
    IF( order >= 2 .OR. order == -2 ) THEN
      t1071 = t332 * t13
      t1072 = t4 * t1071
      t1077 = 0.10e2_dp / 0.9e1_dp * t1072 * t340 * t68 * t56 * t58
      t1078 = t44 * t347
      t1079 = t4 * t1078
      t1081 = 0.8e1_dp / 0.3e1_dp * t1079 * t343
      t1082 = t14 * t371
      t1085 = 0.4e1_dp / 0.3e1_dp * t339 * t1082 * t342
      t1086 = t4 * t515
      t1087 = t6 * t401
      t1090 = 0.4e1_dp / 0.3e1_dp * t1086 * t553 * t1087
      t1091 = t29 * t14
      t1094 = 6._dp * t11 * t1091 * t69
      t1097 = 4._dp * t11 * t348 * t372
      t1098 = t4 * t357
      t1100 = 4._dp * t1098 * t403
      t1113 = t24 * t763
      t1115 = t363 * t312
      t1119 = t46 * t31 * t6
      t1123 = t27 * t390 * t31
      t1126 = 0.10e2_dp / 0.9e1_dp * t16 * t3 * t332 * t15 * t312 + 0.8e1_dp / 0.3e1_dp  &
      * t353 * t348 * t6 + (6._dp * t17 * t10 * t29 * t14) + 0.28e2_dp /  &
      0.9e1_dp * t22 * t1113 * t1115 + 0.32e2_dp / 0.3e1_dp * t362 * t1119 + (20._dp  &
      * t25 * t1123)
      t1127 = t1126 * t68
      t1129 = t11 * t15 * t1127
      t1132 = 2._dp * t375 * t1082 * t402
      t1134 = 0.1e1_dp / t376 / t67
      t1135 = t401 ** 2
      t1136 = t1134 * t1135
      t1139 = 2._dp * t375 * t340 * t1136
      t1168 = 0.28e2_dp / 0.9e1_dp * t35 * t1113 * t1115 + 0.32e2_dp / 0.3e1_dp * t378  &
      * t1119 + (20._dp * t36 * t1123) + 0.40e2_dp / 0.9e1_dp * t39 * t41 *  &
      t463 * t385 * t312 + 0.50e2_dp / 0.3e1_dp * t384 * t390 * t48 * t6 + 0.30e2_dp  &
      * t42 * t44 / t28 / t346 * t48 + (72._dp * t55 * t60 / t61 / t12  &
      * t63)
      t1169 = t377 * t1168
      t1171 = t375 * t340 * t1169
      t1173 = t332 * t56 * t58
      t1176 = t1077 + t1081 - t1085 + t1090 + t1094 - t1097 + t1100 + t1129  &
      - t1132 + t1139 - t1171 + 0.10e2_dp / 0.9e1_dp * t73 * t1173
      d2Qrhorho = f94 * t1176 * t77
      t1178 = t411 * t338
      t1181 = t14 * t424
      t1185 = t6 * t440
      t1208 = -0.4e1_dp / 0.3e1_dp * t415 * t352 * t354 - (4._dp * t416 * t358)  &
      - 0.16e2_dp / 0.3e1_dp * t420 * t361 * t364 - (16._dp * t421 * t368)
      t1209 = t1208 * t68
      t1214 = t411 * t18
      t1220 = t1134 * t401 * t440
      t1236 = -0.16e2_dp / 0.3e1_dp * t428 * t361 * t364 - (16._dp * t429 * t368)  &
      - 0.25e2_dp / 0.3e1_dp * t432 * t383 * t386 - (25._dp * t433 * t392)  &
      - (48._dp * t437 * t398)
      t1237 = t377 * t1236
      t1240 = -0.4e1_dp / 0.3e1_dp * t1178 * t343 - 0.2e1_dp / 0.3e1_dp * t339 * t1181  &
      * t342 + 0.2e1_dp / 0.3e1_dp * t1086 * t553 * t1185 - (4._dp * t412 *  &
      t349) - (2._dp * t11 * t348 * t425) + (2._dp * t1098 * t442) + (2._dp  &
      * t412 * t373) + (t11 * t15 * t1209) - t375 * t1082 * t441  &
      - (2._dp * t1214 * t403) - t375 * t1181 * t402 + 0.2e1_dp * t375 * t340  &
      * t1220 - t375 * t340 * t1237
      d2Qrhondrho = f94 * t1240 * t77
      t1242 = t119 * t13
      t1243 = t340 * t68
      t1257 = 2._dp * a1 * t3 * t19 + 12._dp * a2 * t1 * t24 * t32
      t1258 = t1257 * t68
      t1264 = t440 ** 2
      t1265 = t1134 * t1264
      t1281 = 12._dp * a3 * t1 * t24 * t32 + 20._dp * a4 * t419 * t41 * t49 + 30._dp  &
      * a5 * t21 * t54 * t65
      t1282 = t377 * t1281
      t1285 = 2._dp * t1242 * t1243 + 4._dp * t412 * t426 - 4._dp * t1214 * t442 + t11  &
      * t15 * t1258 - 2._dp * t375 * t1181 * t441 + 2._dp * t375 * t340 * t1265  &
      - t375 * t340 * t1282
      d2Qndrhondrho = f94 * t1285 * t77
      t1287 = t78 ** 2
      t1294 = t166 * t481
      t1297 = t453 ** 2
      t1301 = t1077 + t1081 - t1085 + t1090 + t1094 - t1097 + t1100 + t1129  &
      - t1132 + t1139 - t1171
      t1304 = t166 * t181
      t1305 = t752 * t453
      t1306 = t1305 * t754
      t1309 = t587 ** 2
      t1310 = t154 * t1309
      t1317 = t1077 + t1081 - t1085 + t1090 + t1094 - t1097 + t1100 + t1129  &
      - t1132 + t1139 - t1171 + 0.10e2_dp / 0.9e1_dp * t72 * t332 * t312
      t1318 = t155 * t1317
      t1324 = 0.1e1_dp / t169 / t156 / t155
      t1325 = t83 * t1324
      t1326 = t156 * t1309
      t1355 = 0.10e2_dp / 0.9e1_dp * t1072 * t468 * t312 + 0.8e1_dp / 0.3e1_dp * t1079  &
      * t469 - 0.4e1_dp / 0.3e1_dp * t339 * t14 * t478 * t6 + (6._dp * t11 *  &
      t1091 * t87) - 0.4e1_dp * (t11) * t348 * t478 + (t11 * t15 *  &
      (F2 * t1126 * t68 - 2._dp * t475 * t402 + 2._dp * t85 * t1136 - t85 * t1169))
      t1379 = t495 * t82
      t1381 = 0.1e1_dp / t110 / t1379
      t1382 = t109 * t1381
      t1392 = t503 * t507
      t1402 = t116 * t1134 * t1 * t3
      t1409 = t116 * t377 * t4 * t44
      t1433 = f94 * t1126 * t117 * t121 - (2._dp * t1392 * t511) - 0.4e1_dp  &
      / 0.3e1_dp * t503 * t513 * t518 - (4._dp * t504 * t522) + (2._dp * t1402  &
      * t18 * t509 * t1135) + 0.4e1_dp / 0.3e1_dp * t1409 * t120 * t1087 +  &
      (4._dp * t508 * t357 * t510) - (t508 * t18 * t509 * t1168) +  &
      0.10e2_dp / 0.9e1_dp * t514 * t1071 * t14 * t77 * t56 * t58 + 0.8e1_dp / 0.3e1_dp  &
      * t514 * t1078 * t14 * t517 + 0.6e1_dp * t118 * t119 * t1091 * t77
      t1437 = t525 ** 2
      t1442 = t96 * t115 * t525
      t1445 = t96 * t124
      t1456 = t533 * t1242
      t1457 = t377 * t77
      t1473 = t13 * sscale
      t1478 = t125 * t538 * t129
      t1479 = sscale * t551
      t1491 = 0.1e1_dp / t132 / t131
      t1492 = t556 ** 2
      t1498 = t371 * t377
      t1501 = t34 * t1134
      t1520 = t567 * t142
      t1524 = t110 * t82 * t144
      t1531 = -((t97 * (t99 * t1355 * t82 + 2._dp * t99 * t481 * t453 + t99  &
      * t90 * t1301 + 2._dp * t102 * t1297 + 2._dp * t102 * t82 * t1301 + 6._dp *  &
      t104 * t82 * t1297 + 3._dp * t104 * t92 * t1301) * t113) - (7._dp * t493  &
      * t499) + 0.63e2_dp / 0.4e1_dp * (t108) * (t1382) * (t1297)  &
      - 0.7e1_dp / 0.2e1_dp * (t108) * (t498) * (t1301) - t502  &
      * t1433 * t123 * t136 - t502 * t1437 * t123 * t136 + (4._dp * t1442  &
      * t563) + 0.2e1_dp * t1445 * t530 * (0.2e1_dp / 0.3e1_dp * t533 * t352 *  &
      t13 * t1243 * t517 + (2._dp * t534 * t348 * t131) - (t534 * t15  &
      * t552) + t1456 * t340 * t1457 * t401) * t562 + 0.2e1_dp * t502 *  &
      t531 * t537 * (0.4e1_dp / 0.9e1_dp * t125 * t126 * t315 * t540 * t132 *  &
      t56 * t58 + 0.2e1_dp / 0.3e1_dp * t539 * t1473 * t542 - t1478 * t1479 *  &
      t6 * t556 / 0.3e1_dp + (2._dp * t127 * t128 * t347 * t133) - t550 * t1473  &
      * t557 - t550 * t540 * t1491 * t1492 / 0.4e1_dp + t550 * t540 *  &
      t551 * (t1127 * t77 - 2._dp * t1498 * t554 + 2._dp * t1501 * t77 * t1135  &
      - t553 * t77 * t1168) / 0.2e1_dp)) * t140 * t145 - (7._dp * t1520  &
      * t572) - 0.35e2_dp / 0.4e1_dp * (t569) * (t1524) * (t1297)  &
      - 0.7e1_dp / 0.2e1_dp * (t569) * (t571) * (t1301)
      t1535 = t151 * t575
      t1543 = (3._dp * t1294 * t804) - (2._dp * t799 * t747 * t1297) + (t799  &
      * t800 * t1301) - (3._dp * t1304 * t1306) + (3._dp * t799  &
      * t803 * t1310) + 0.3e1_dp / 0.2e1_dp * (t799) * (t803) * (t1318)  &
      - 0.27e2_dp / 0.4e1_dp * (t799) * (t1325) * (t1326) -  &
      t151 * t1531 * t83 * t159 + (2._dp * t1535 * t781) + (5._dp * t1535  &
      * t784) - (2._dp * t779 * t656 * t1297)
      t1546 = t151 * t177
      t1547 = t581 * t453
      t1548 = t1547 * t588
      t1551 = t156 * t1317
      t1555 = t168 * t1309
      t1559 = t156 ** 2
      t1562 = 0.1e1_dp / t158 / t1559 / t155
      t1563 = t83 * t1562
      t1564 = t1559 * t1309
      t1568 = t149 * t163
      t1572 = t81 * t93
      t1573 = t791 * t453
      t1579 = 0.1e1_dp / t162 / t155
      t1580 = t83 * t1579
      t1592 = t779 * t780 * (t1301) - (5._dp * t1546 * t1548) + 0.5e1_dp  &
      / 0.2e1_dp * t779 * t783 * t1551 + 0.10e2_dp * t779 * t783 * t1555 - 0.75e2_dp  &
      / 0.4e1_dp * t779 * t1563 * t1564 - (2._dp * t81 * t1568 * t1297)  &
      - t1572 * t1573 * t587 + (t81 * t787 * t1301) - 0.3e1_dp / 0.4e1_dp  &
      * (t81) * (t1580) * (t1309) + (t81 * t792 * t1317)  &
      / 0.2e1_dp - t166 * t1355 * t83 * t170 + (2._dp * t1294 * t801)
      t1599 = t576 * t582
      t1602 = 0.1e1_dp / t1379
      t1609 = t315 * t56 * t58
      t1615 = t189 * t347
      t1624 = 0.1e1_dp / t690
      t1625 = t409 ** 2
      t1633 = f94 / t294
      t1649 = d2exeirhorho(Q,dQrho,d2Qrhorho)
      t1658 = t148 * t620 * t581
      t1663 = (t1543 + t1592) * omega * t128 + (10._dp * t583 * t271 * t1555)  &
      + (5._dp * t1599 * t589) + (12._dp * t148 * t1602 * t1297) -  &
      (2._dp * t613 * t578) + 0.4e1_dp / 0.9e1_dp * t624 * t1609 - 0.8e1_dp / 0.3e1_dp  &
      * t669 * t192 * t631 + (2._dp * t186 * t1615) + 0.10e2_dp / 0.3e1_dp  &
      * t576 * t656 * t659 + 0.5e1_dp / 0.2e1_dp * (t583) * (t271) * (t1551)  &
      + ((20._dp * t318 * t1624 * t1625 - 4._dp * t318 * t713 * t1176  &
      - 12._dp * t1633 * t1625 + 3._dp * t718 * t1176 + 6._dp * t243 * t717 * t1625  &
      - 2._dp * t243 * t321 * t1176 - 2._dp * t280 * t321 * t1625 + t280 *  &
      t275 * t1176 + t326 * t1649) * alpha8 * t333) - (3._dp * t148 * t620  &
      * t1301) - (15._dp * t1658 * t271 * t588 * t453)
      t1677 = t256 ** 2
      t1680 = 0.1e1_dp / t290 / t1677 / t690
      t1681 = t1677 * t256
      t1682 = t1680 * t1681
      t1686 = t689 * t257
      t1693 = t293 ** 2
      t1697 = f2732 / t295 / t1693 / t250
      t1709 = 0.1e1_dp / t252 / t293 / t250
      t1710 = t1709 * t293
      t1714 = t638 * t75
      t1722 = 0.1e1_dp / t227 / t250
      t1723 = t1722 * t1625
      t1726 = t602 * t1176
      t1729 = 0.147e3_dp / 0.4e1_dp * t288 * t1682 * t1625 - 0.21e2_dp * t288 * t1686  &
      * t1625 - 0.7e1_dp / 0.2e1_dp * t288 * t691 * t1176 - 0.75e2_dp / 0.4e1_dp  &
      * t1697 * t1693 * t1625 + 0.10e2_dp * t697 * t251 * t1625 + 0.5e1_dp /  &
      0.2e1_dp * t697 * t293 * t1176 + 0.27e2_dp / 0.4e1_dp * t298 * t1710 * t1625  &
      - 0.3e1_dp * t298 * t1714 * t1625 - 0.3e1_dp / 0.2e1_dp * t298 * t701 * t1176  &
      - 0.3e1_dp / 0.4e1_dp * t280 * t1723 + t280 * t1726 / 0.2e1_dp
      t1731 = d2exerrhorho(Q,dQrho,d2Qrhorho)
      t1740 = t148 * t1602
      t1767 = t56 ** 2
      t1768 = t58 ** 2
      t1770 = t1767 * t1768 * t28
      t1781 = A * t713
      t1786 = f98 * t717
      t1791 = r1 * t321
      t1802 = ((-2._dp * t823 * t1625 + t276 * t1176 + 6._dp * t197 * t319 *  &
      t1625 - 2._dp * t197 * t273 * t1176 + t243 * t1649) * alpha4 * t247) +  &
      (t97 * t1729 + t307 * t1731) * alpha7 * t316 - 0.40e2_dp / 0.9e1_dp * t657  &
      * t193 * t463 * t312 - (12._dp * t1740 * t194 * t44 * t1297) +  &
      (6._dp * t285 * t60 * t29) + ((2._dp * t197 * t273 * t1625 - t197  &
      * t241 * t1176 - f98 * t1649) * alpha2 * t153) - (4._dp * t832 *  &
      t449) - (2._dp * t759 * t578) + (2._dp * t237 * t1615) - (6._dp *  &
      t576 * t621) + f12 * t1355 * t94 + 0.70e2_dp / 0.9e1_dp * t761 / t8 / t1770  &
      * t56 * t58 + 0.40e2_dp / 0.9e1_dp * t815 * t463 * t56 * t58 + ((12._dp  &
      * t1781 * t1625 - 3._dp * t820 * t1176 - 6._dp * t1786 * t1625 + 2._dp * t823  &
      * t1176 + 2._dp * t243 * t1791 * t1625 - t243 * t826 * t1176 - t280  &
      * t1649) * alpha6 * t284 * t286)
      t1804 = t620 * t159
      t1805 = t148 * t1804
      t1816 = t576 * t620
      t1827 = t148 * t582 * t193
      t1835 = t148 * t149 * t1562
      t1839 = C * t620
      t1846 = t75 * t273
      t1853 = t771 * t275
      t1854 = t199 * t409
      t1860 = t1531 * E * t149
      t1868 = f916 * t1709
      t1880 = 0.1e1_dp / t258 / t1677 / t222
      t1881 = t1880 * t1677
      t1885 = t644 * t223
      t1898 = -(10._dp * t1805 * t658 * t6 * t453) - 0.14e2_dp / 0.3e1_dp * t711  &
      * t311 * t765 - 0.16e2_dp / 0.3e1_dp * t728 * t330 * t465 + (6._dp * t1816  &
      * t812) + (2._dp * t81 * t149 * t1297) + 0.28e2_dp / 0.9e1_dp * t629  &
      * t763 * t56 * t58 - 0.25e2_dp / 0.3e1_dp * t1827 * t332 * t156 * t587  &
      * r3 * t5 - 0.75e2_dp / 0.4e1_dp * t1835 * t271 * t1564 + (6._dp * t91  &
      * t1839 * t1297) + (t197 * (t1649 + (t1176 * t200 - 2._dp * t1625  &
      * t241 + 2._dp * t1846 * t1625 - t769 * t1176) * t239 * t199 - t1853 *  &
      t1854 + t772 * t409)) - t1860 * t195 - (t81 * t93 * t1301) +  &
      (t97 * (0.3e1_dp / 0.4e1_dp * t243 * t1723 - t243 * t1726 / 0.2e1_dp - 0.27e2_dp  &
      / 0.4e1_dp * (t1868) * (t293) * (t1625) + (3._dp * t639  &
      * t75 * t1625) + 0.3e1_dp / 0.2e1_dp * (t639) * (t250) * (t1176)  &
      + 0.75e2_dp / 0.4e1_dp * (t255) * (t1881) * (t1625) -  &
      (10._dp * t255 * t1885 * t1625) - 0.5e1_dp / 0.2e1_dp * (t255) * (t645)  &
      * (t1176)) - t267 * t1731) * alpha5 * t271
      t1906 = 0.1e1_dp / t205 / t222
      t1907 = t97 * t1906
      t1923 = t176 * t575
      t1932 = t176 * t147 * t149
      t1938 = t93 * t1562
      t1951 = t180 * t481
      t1956 = t620 * t170
      t1961 = t180 * t90 * t149
      t1967 = t93 * t1324
      t1977 = -t176 * t1531 * t93 * t159 + (4._dp * t1923 * t736) + (5._dp  &
      * t1923 * t740) - (6._dp * t735 * t1804 * t1297) - (10._dp * t1932  &
      * t1548) + (2._dp * t735 * t656 * t1301) - 0.75e2_dp / 0.4e1_dp * (t735)  &
      * (t1938) * (t1564) + (10._dp * t735 * t739 * t1555)  &
      + 0.5e1_dp / 0.2e1_dp * (t735) * (t739) * (t1551) - t180  &
      * t1355 * t93 * t170 + (4._dp * t1951 * t748) + (3._dp * t1951 * t755)  &
      - (6._dp * t746 * t1956 * t1297) - (6._dp * t1961 * t1306) +  &
      (2._dp * t746 * t747 * t1301) - 0.27e2_dp / 0.4e1_dp * (t746) * (t1967)  &
      * (t1326) + (3._dp * t746 * t753 * t1310) + 0.3e1_dp / 0.2e1_dp  &
      * (t746) * (t753) * (t1318)
      t1983 = 0.1e1_dp / t224 / t690
      t1984 = t1983 * t256
      t1988 = t596 * t199
      t1995 = f98 * t1722
      t2028 = -0.4e1_dp / 0.3e1_dp * t685 * t72 * t406 - 0.2e1_dp / 0.3e1_dp * t679 *  &
      omega * t626 + (0.3e1_dp / 0.4e1_dp * t197 * t1907 * t1625 - t197 * t673  &
      * t1176 / 0.2e1_dp - t209 * t1731) * alpha1 * t214 + 0.4e1_dp / 0.9e1_dp  &
      * t174 * t1609 + t1977 * t185 * t190 + 0.10e2_dp / 0.9e1_dp * t592 * t1173  &
      + (t97 * (0.27e2_dp / 0.4e1_dp * t221 * t1984 * t1625 - 0.3e1_dp * t221  &
      * t1988 * t1625 - 0.3e1_dp / 0.2e1_dp * t221 * t597 * t1176 - 0.3e1_dp / 0.4e1_dp  &
      * t1995 * t1625 + t603 * t1176 / 0.2e1_dp) + f2716 * t1731 * t233)  &
      * alpha3 * t185 * t190 + (3._dp * t810 * t194 * t44 * t1301) - (2._dp  &
      * t91 * t452 * t1301) - 0.2e1_dp / 0.3e1_dp * t808 * t626 - 0.10e2_dp  &
      / 0.3e1_dp * t653 * t193 * t817 + t1860 - (4._dp * t617 * t454) + 0.88e2_dp  &
      / 0.9e1_dp * t457 / t9 / t1770 * t56 * t58
      e_rho_rho = e_rho_rho   -0.4e1_dp / 0.9e1_dp / t1287 * f89 * t336 - 0.8e1_dp / 0.3e1_dp * t446 *  &
      t836 - t80 * (t1663 + t1802 + t1898 + t2028) * Clda
      t2059 = t156 * t1240
      t2072 = t587 * t444
      t2073 = t581 * t156 * t2072
      t2099 = -0.4e1_dp / 0.3e1_dp * t1178 * t469 - 0.2e1_dp / 0.3e1_dp * t339 * t14  &
      * t845 * t6 - (4._dp * t412 * t472) - 0.2e1_dp * t11 * t348 * t845 +  &
      (2._dp * t412 * t479) + t11 * t15 * (F2 * t1208 * t68 - t475 *  &
      t441 - t842 * t402 + 2._dp * t85 * t1220 - t85 * t1237)
      t2108 = t444 * t453
      t2142 = t870 * t507
      t2144 = t554 * t440
      t2150 = t116 * t377 * ndrho * t3
      t2173 = (f94 * t1208 * t117 * t121) - t1392 * t874 + (2._dp * t503  &
      * t876 * t121) - (t2142 * t511) + (2._dp * t1402 * t19 * t2144)  &
      - (2._dp * t2150 * t511) - (t508 * t18 * t509 * t1236) -  &
      0.2e1_dp / 0.3e1_dp * t870 * t513 * t518 + 0.2e1_dp / 0.3e1_dp * t1409 * t120  &
      * t1185 - 0.4e1_dp / 0.3e1_dp * t116 * (t876) * t3 * t518 - (2._dp  &
      * t871 * t522) + (2._dp * t508 * t357 * t873) - (4._dp * t877 * t522)
      t2184 = t96 * t115 * t880
      t2197 = t530 * (-2._dp * t532 * ndrho * t119 * t535 - t534 * t15 * t888  &
      + t1456 * t340 * t1457 * t440)
      t2224 = t424 * t377
      t2240 = (t97 * (t99 * t2099 * t82 + t99 * t481 * t444 + t99 * t848  &
      * t453 + t99 * t90 * t1240 + 2._dp * t102 * t2108 + 2._dp * t102 * t82  &
      * t1240 + 6._dp * t104 * t486 * t444 + 3._dp * t104 * t92 * t1240) * t113)  &
      - 0.7e1_dp / 0.2e1_dp * t493 * t867 - 0.7e1_dp / 0.2e1_dp * t865 * t499 + 0.63e2_dp  &
      / 0.4e1_dp * (t108) * (t109) * (t1381) * (t453)  &
      * (t444) - 0.7e1_dp / 0.2e1_dp * (t108) * (t498) * (t1240)  &
      - t502 * t2173 * t123 * t136 - t502 * t525 * t880 * t123 * t136  &
      + (2._dp * t1442 * t898) + (2._dp * t2184 * t563) + (2._dp * t1445  &
      * t2197 * t562) + 0.2e1_dp * t502 * t531 * t537 * (-t884 * t27 * t129  &
      * t133 * t6 / 0.3e1_dp - t1478 * t1479 * t6 * t891 / 0.6e1_dp - t885  &
      * t1473 * t132 - t550 * t1473 * t892 / 0.2e1_dp + t885 * t558 / 0.2e1_dp  &
      - t550 * t540 * t1491 * t556 * t891 / 0.4e1_dp + t550 * t540 * t551  &
      * (t1209 * t77 - t1498 * t889 - t2224 * t554 + 2._dp * t1501 * t2144  &
      - t553 * t77 * t1236) / 0.2e1_dp)
      t2245 = t902 * t142
      t2254 = -t2240 * t140 * t145 - 0.7e1_dp / 0.2e1_dp * t1520 * t904 - 0.7e1_dp  &
      / 0.2e1_dp * t2245 * t572 - 0.35e2_dp / 0.4e1_dp * t569 * t1524 * t2108 -  &
      0.7e1_dp / 0.2e1_dp * t569 * t571 * t1240
      t2258 = t1547 * t918
      t2267 = t151 * t907
      t2269 = t166 * t167
      t2271 = t752 * t154 * t2072
      t2274 = t81 * t792 * t1240 / 0.2e1_dp - t1572 * t1573 * t444 / 0.2e1_dp +  &
      t81 * t787 * t1240 + 0.5e1_dp / 0.2e1_dp * t779 * t783 * t2059 - 0.2e1_dp  &
      * t81 * t149 * t163 * t453 * t444 + t1294 * t930 + 0.5e1_dp / 0.2e1_dp *  &
      t1535 * t919 - 0.5e1_dp / 0.2e1_dp * t1546 * t2073 - t151 * t2254 * t83  &
      * t159 - 0.5e1_dp / 0.2e1_dp * t1546 * t2258 + t1535 * t916 - 0.2e1_dp * t779  &
      * t656 * t2108 + t779 * t780 * t1240 + t2267 * t781 + (3._dp *  &
      t2269 * t2271)
      t2280 = t151 * t152
      t2282 = t1562 * t1559 * t2072
      t2285 = t166 * t848
      t2291 = t752 * t155 * t2072
      t2297 = t1305 * t932
      t2305 = t581 * t168 * t2072
      t2311 = t155 * t1240
      t2316 = t1324 * t156 * t2072
      t2323 = -(2._dp * t799 * t747 * t2108) + (t799 * t800 * t1240)  &
      - 0.75e2_dp / 0.4e1_dp * t2280 * t2282 + 0.3e1_dp / 0.2e1_dp * t2285 * t804 +  &
      0.5e1_dp / 0.2e1_dp * t2267 * t784 - 0.3e1_dp / 0.2e1_dp * t1304 * t2291 - t166  &
      * t2099 * t83 * t170 - 0.3e1_dp / 0.2e1_dp * t1304 * t2297 - t1572 * t791  &
      * t587 * t444 / 0.2e1_dp + 0.10e2_dp * t2280 * t2305 + 0.3e1_dp / 0.2e1_dp  &
      * t1294 * t933 + t2285 * t801 + 0.3e1_dp / 0.2e1_dp * (t799) * (t803)  &
      * (t2311) - 0.27e2_dp / 0.4e1_dp * t2269 * t2316 - 0.3e1_dp / 0.4e1_dp  &
      * t84 * t1579 * t587 * t444
      t2329 = t908 * t582
      t2348 = t176 * t907
      t2363 = t176 * t177
      t2371 = -t176 * t2254 * t93 * t159 + (2._dp * t1923 * t942) + 0.5e1_dp  &
      / 0.2e1_dp * (t1923) * (t945) + (2._dp * t2348 * t736) - (6._dp  &
      * t735 * t1804 * t2108) - (5._dp * t1932 * t2258) + (2._dp * t735  &
      * t656 * t1240) + 0.5e1_dp / 0.2e1_dp * (t2348) * (t740) - (5._dp  &
      * t1932 * t2073) - 0.75e2_dp / 0.4e1_dp * t2363 * t2282 + 0.10e2_dp *  &
      t2363 * t2305 + 0.5e1_dp / 0.2e1_dp * (t735) * (t739) * (t2059)
      t2379 = t180 * t848
      t2401 = -t180 * t2099 * t93 * t170 + (2._dp * t1951 * t951) + 0.3e1_dp  &
      / 0.2e1_dp * (t1951) * (t954) + (2._dp * t2379 * t748) - (6._dp  &
      * t746 * t1956 * t2108) - (3._dp * t1961 * t2297) + (2._dp * t746  &
      * t747 * t1240) + 0.3e1_dp / 0.2e1_dp * (t2379) * (t755) - (3._dp  &
      * t1961 * t2291) - 0.27e2_dp / 0.4e1_dp * t95 * t2316 + 0.3e1_dp * t95  &
      * t2271 + 0.3e1_dp / 0.2e1_dp * (t746) * (t753) * (t2311)
      t2405 = -0.25e2_dp / 0.6e1_dp * t1827 * t816 * t5 * t156 * t444 + 0.12e2_dp  &
      * t148 * t1602 * t453 * t444 - 0.5e1_dp * t1805 * t658 * t6 * t444 +  &
      0.6e1_dp * t746 * t621 * t444 + (t2274 + t2323) * omega * t128 - (2._dp  &
      * t1035 * t449) + 0.5e1_dp / 0.2e1_dp * t2329 * t589 + 0.5e1_dp / 0.2e1_dp  &
      * t1599 * t965 - t81 * t93 * t1240 - 0.7e1_dp / 0.3e1_dp * t1052 * t311  &
      * t765 - t937 * t626 / 0.3e1_dp + (t2371 + t2401) * t185 * t190
      t2410 = t2254 * E * t149
      t2437 = t698 * t444
      t2442 = t75 * t409 * t444
      t2449 = t1722 * t409 * t444
      t2452 = t602 * t1240
      t2455 = 0.147e3_dp / 0.4e1_dp * t288 * t1680 * t1681 * t409 * t444 - 0.21e2_dp  &
      * t288 * t689 * t257 * t409 * t444 - 0.7e1_dp / 0.2e1_dp * t288 * t691  &
      * t1240 - 0.75e2_dp / 0.4e1_dp * t1697 * t1693 * t409 * t444 + 0.10e2_dp  &
      * t697 * t251 * t409 * t444 + 0.5e1_dp / 0.2e1_dp * t697 * t293 * t1240  &
      + 0.27e2_dp / 0.4e1_dp * t298 * t1709 * t2437 - 0.3e1_dp * t298 * t638 * t2442  &
      - 0.3e1_dp / 0.2e1_dp * t298 * t701 * t1240 - 0.3e1_dp / 0.4e1_dp * t280 *  &
      t2449 + t280 * t2452 / 0.2e1_dp
      t2457 = d2exerrhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t2473 = t409 * t444
      t2484 = d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t2499 = -(3._dp * t908 * t621) - t1001 * t578 - t2410 * t195 + t2410  &
      + (t97 * t2455 + t307 * t2457) * alpha7 * t316 + 0.5e1_dp / 0.2e1_dp *  &
      t583 * t271 * t2059 - 0.12e2_dp * t148 * t1602 * t159 * t271 * t2108  &
      - t982 * omega * t626 / 0.3e1_dp + ((-2._dp * t823 * t2473 + t276 *  &
      t1240 + 6._dp * t197 * t319 * t409 * t444 - 2._dp * t197 * t273 * t1240 +  &
      t243 * t2484) * alpha4 * t247) - 0.15e2_dp / 0.2e1_dp * t1658 * t271 * t588  &
      * (t444) - 0.4e1_dp / 0.3e1_dp * t1009 * t192 * t631 - 0.3e1_dp * t148  &
      * t620 * (t1240)
      t2512 = t199 * t444
      t2529 = t721 * t444
      t2543 = t908 * t620
      t2573 = (3._dp * t1816 * t962) + (t197 * (t2484 + (t1240 * t200  &
      - 2._dp * t681 * t444 + 2._dp * t1846 * t2473 - t769 * t1240) * t239 * t199  &
      - t1853 * t2512 + t772 * t444)) - 0.5e1_dp / 0.3e1_dp * t1024 * t193  &
      * t817 + ((12._dp * t1781 * t2473 - 3._dp * t820 * t1240 - 6._dp * t1786 *  &
      t2473 + 2._dp * t823 * t1240 + 2._dp * t243 * r1 * t2529 - t243 * t826 *  &
      t1240 - t280 * t2484) * alpha6 * t284 * t286) + (2._dp * t81 * t149  &
      * t453 * t444) + (3._dp * t2543 * t812) + f12 * t2099 * t94 + (10._dp  &
      * t583 * t271 * t168 * t587 * t444) + ((2._dp * t197 * t663  &
      * t444 - t197 * t241 * t1240 - f98 * t2484) * alpha2 * t153) - 0.8e1_dp  &
      / 0.3e1_dp * t1066 * t330 * t465 - 0.2e1_dp / 0.3e1_dp * t988 * t72 * t406  &
      - 0.75e2_dp / 0.4e1_dp * (t1835) * (t271) * (t1559) * (t587)  &
      * (t444)
      t2688 = -(2._dp * t91 * t452 * t1240) - (2._dp * t617 * t851) + (t97  &
      * (0.3e1_dp / 0.4e1_dp * t243 * t2449 - t243 * t2452 / 0.2e1_dp - 0.27e2_dp  &
      / 0.4e1_dp * t1868 * t2437 + (3._dp * t639 * t2442) + 0.3e1_dp / 0.2e1_dp *  &
      (t639) * (t250) * (t1240) + 0.75e2_dp / 0.4e1_dp * t255 * t1880  &
      * t1677 * t409 * t444 - 0.10e2_dp * t255 * t644 * t223 * t409 *  &
      t444 - 0.5e1_dp / 0.2e1_dp * t255 * t645 * (t1240)) - t267 * t2457)  &
      * alpha5 * t271 - (3._dp * t576 * t910) - (2._dp * t849 * t454) -  &
      t958 * t578 + (0.20e2_dp * t318 * t1624 * t409 * t444 - 0.4e1_dp * t318  &
      * t713 * (t1240) - (12._dp * t1633 * t2473) + (3._dp * t718 *  &
      t1240) + 0.6e1_dp * t243 * t717 * t409 * t444 - 0.2e1_dp * t243 * t321 *  &
      (t1240) - (2._dp * t280 * t2529) + (t280 * t275 * t1240)  &
      + t326 * t2484) * alpha8 * t333 + (3._dp * t810 * t194 * t44 * t1240)  &
      + (0.3e1_dp / 0.4e1_dp * t197 * t97 * t1906 * t409 * t444 - t197 * t673  &
      * (t1240) / 0.2e1_dp - t209 * t2457) * alpha1 * t214 + 0.5e1_dp  &
      / 0.3e1_dp * t908 * t656 * t659 + (t97 * (0.27e2_dp / 0.4e1_dp * t221 * t1983  &
      * t256 * t409 * t444 - 0.3e1_dp * t221 * t596 * t1854 * t444 - 0.3e1_dp  &
      / 0.2e1_dp * t221 * t597 * (t1240) - 0.3e1_dp / 0.4e1_dp * (t1995)  &
      * (t2473) + (t603 * t1240) / 0.2e1_dp) + f2716 * t2457 * t233)  &
      * alpha3 * t185 * t190 - 0.15e2_dp / 0.2e1_dp * t1658 * t271 * t453  &
      * t156 * t444
      e_ndrho_rho = e_ndrho_rho   -0.4e1_dp / 0.3e1_dp * t446 * t1069 - t80 * (t2405 + t2499 + t2573  &
      + t2688) * Clda
      t2707 = 2._dp * t119 * t88 + 4._dp * t412 * t846 + t11 * t15 * (F2 * t1257  &
      * t68 - 2._dp * t842 * t441 + 2._dp * t85 * t1265 - t85 * t1282)
      t2715 = t444 ** 2
      t2764 = t880 ** 2
      t2774 = t891 ** 2
      t2808 = -((t97 * (t99 * t2707 * t82 + 2._dp * t99 * t848 * t444 + t99  &
      * t90 * t1285 + 2._dp * t102 * t2715 + 2._dp * t102 * t82 * t1285 + 6._dp *  &
      t104 * t82 * t2715 + 3._dp * t104 * t92 * t1285) * t113) - (7._dp * t865  &
      * t867) + 0.63e2_dp / 0.4e1_dp * (t108) * (t1382) * (t2715)  &
      - 0.7e1_dp / 0.2e1_dp * (t108) * (t498) * (t1285) - (t502  &
      * (f94 * t1257 * t117 * t121 - 2._dp * t2142 * t874 + 4._dp * t870 *  &
      t876 * t121 + 2._dp * t1402 * t18 * t509 * t1264 - 4._dp * t2150 * t874 -  &
      t508 * t18 * t509 * t1281 + 2._dp * t116 * t68 * t3 * t18 * t509) * t123  &
      * t136) - (t502 * t2764 * t123 * t136) + (4._dp * t2184 *  &
      t898) + (2._dp * t1445 * t2197 * t897) + 0.2e1_dp * (t502) * t531  &
      * t537 * (t885 * t893 - t550 * t540 * t1491 * t2774 / 0.4e1_dp + t550  &
      * t540 * t551 * (t1258 * t77 - 2._dp * t2224 * t889 + 2._dp * t1501  &
      * t77 * t1264 - t553 * t77 * t1281) / 0.2e1_dp)) * t140 * t145 - (7._dp  &
      * t2245 * t904) - 0.35e2_dp / 0.4e1_dp * (t569) * (t1524) * (t2715)  &
      - 0.7e1_dp / 0.2e1_dp * (t569) * (t571) * (t1285)
      t2810 = t2808 * E * t149
      t2838 = t1722 * t2715
      t2841 = t602 * t1285
      t2844 = 0.147e3_dp / 0.4e1_dp * t288 * t1682 * t2715 - 0.21e2_dp * t288 * t1686  &
      * t2715 - 0.7e1_dp / 0.2e1_dp * t288 * t691 * t1285 - 0.75e2_dp / 0.4e1_dp  &
      * t1697 * t1693 * t2715 + 0.10e2_dp * t697 * t251 * t2715 + 0.5e1_dp /  &
      0.2e1_dp * t697 * t293 * t1285 + 0.27e2_dp / 0.4e1_dp * t298 * t1710 * t2715  &
      - 0.3e1_dp * t298 * t1714 * t2715 - 0.3e1_dp / 0.2e1_dp * t298 * t701 * t1285  &
      - 0.3e1_dp / 0.4e1_dp * t280 * t2838 + t280 * t2841 / 0.2e1_dp
      t2846 = d2exerndrhondrho(Q,dQndrho,d2Qndrhondrho)
      t2875 = d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
      t2880 = t2715 * t156
      t2884 = t1559 * t2715
      t2927 = t156 * t1285
      t2963 = t2810 + (t97 * t2844 + t307 * t2846) * alpha7 * t316 + (3._dp  &
      * t810 * t194 * t44 * t1285) - (6._dp * t908 * t910) - (3._dp *  &
      t148 * t620 * t1285) + (0.3e1_dp / 0.4e1_dp * t197 * t1907 * t2715 - t197  &
      * t673 * (t1285) / 0.2e1_dp - t209 * t2846) * alpha1 * t214 +  &
      (0.2e1_dp * t197 * t273 * t2715 - t197 * t241 * (t1285) - f98 * t2875)  &
      * alpha2 * t153 - (15._dp * t1658 * t271 * t2880) - 0.75e2_dp /  &
      0.4e1_dp * (t1835) * (t271) * (t2884) + (0.20e2_dp * t318 *  &
      t1624 * t2715 - 0.4e1_dp * t318 * t713 * (t1285) - 0.12e2_dp * t1633  &
      * t2715 + (3._dp * t718 * t1285) + 0.6e1_dp * t243 * t717 * t2715 -  &
      0.2e1_dp * t243 * t321 * (t1285) - 0.2e1_dp * t280 * t321 * t2715 +  &
      t280 * t275 * (t1285) + t326 * t2875) * alpha8 * t333 + t197  &
      * (t2875 + ((t1285 * t200) - 0.2e1_dp * t2715 * t241 + 0.2e1_dp * t1846  &
      * t2715 - (t769 * t1285)) * t239 * t199 - t971 * t275 * t2512  &
      + t972 * t444) + 0.5e1_dp / 0.2e1_dp * (t583) * (t271) * (t2927)  &
      + (t97 * (0.3e1_dp / 0.4e1_dp * t243 * t2838 - t243 * t2841 / 0.2e1_dp  &
      - 0.27e2_dp / 0.4e1_dp * t1868 * t293 * t2715 + 0.3e1_dp * t639 * t75  &
      * t2715 + 0.3e1_dp / 0.2e1_dp * t639 * t250 * (t1285) + 0.75e2_dp / 0.4e1_dp  &
      * t255 * t1881 * t2715 - 0.10e2_dp * t255 * t1885 * t2715 - 0.5e1_dp  &
      / 0.2e1_dp * t255 * t645 * (t1285)) - t267 * t2846) * alpha5 * (t271)  &
      + (5._dp * t2329 * t965) - (t81 * t93 * t1285)
      t2971 = t2715 * t155
      t2975 = t154 * t2715
      t2979 = t155 * t1285
      t2994 = t168 * t2715
      t3001 = -(2._dp * t799 * t747 * t2715) + (t799 * t800 * t1285)  &
      - (3._dp * t799 * t753 * t2971) + (3._dp * t799 * t803 * t2975) +  &
      0.3e1_dp / 0.2e1_dp * (t799) * (t803) * (t2979) - 0.27e2_dp / 0.4e1_dp  &
      * (t799) * (t1325) * (t2880) - (2._dp * t779 * t656  &
      * t2715) + (t779 * t780 * t1285) - (5._dp * t779 * t739 * t2880)  &
      + (10._dp * t779 * t783 * t2994) + 0.5e1_dp / 0.2e1_dp * (t779)  &
      * (t783) * (t2927)
      t3033 = -0.75e2_dp / 0.4e1_dp * t779 * t1563 * t2884 - (2._dp * t81 * t1568  &
      * t2715) - (t81 * t93 * t791 * t2715) + (t81 * t787 * t1285)  &
      - 0.3e1_dp / 0.4e1_dp * (t81) * (t1580) * (t2715) + (t81  &
      * t792 * t1285) / 0.2e1_dp - t166 * t2707 * t83 * t170 + (2._dp  &
      * t2285 * t930) + (3._dp * t2285 * t933) - t151 * t2808 * t83 *  &
      t159 + (2._dp * t2267 * t916) + (5._dp * t2267 * t919)
      t3139 = -t176 * t2808 * t93 * t159 + (4._dp * t2348 * t942) + (5._dp  &
      * t2348 * t945) - (6._dp * t735 * t1804 * t2715) - (10._dp * t735  &
      * t582 * t2880) + (2._dp * t735 * t656 * t1285) - 0.75e2_dp / 0.4e1_dp  &
      * (t735) * (t1938) * (t2884) + (10._dp * t735 * t739  &
      * t2994) + 0.5e1_dp / 0.2e1_dp * (t735) * (t739) * (t2927) -  &
      t180 * t2707 * t93 * t170 + (4._dp * t2379 * t951) + (3._dp * t2379  &
      * t954) - (6._dp * t746 * t1956 * t2715) - (6._dp * t746 * t149  &
      * t752 * t2971) + (2._dp * t746 * t747 * t1285) - 0.27e2_dp / 0.4e1_dp  &
      * (t746) * (t1967) * (t2880) + (3._dp * t746 * t753 *  &
      t2975) + 0.3e1_dp / 0.2e1_dp * (t746) * (t753) * (t2979)
      t3167 = f12 * t2707 * t94 + (t3001 + t3033) * omega * t128 + (2._dp  &
      * t81 * t149 * t2715) + (6._dp * t2543 * t962) - t2810 * t195 + (10._dp  &
      * t583 * t271 * t2994) + (12._dp * t148 * t1602 * t2715) +  &
      (t97 * (0.27e2_dp / 0.4e1_dp * (t221) * (t1984) * (t2715) -  &
      (3._dp * t221 * t1988 * t2715) - 0.3e1_dp / 0.2e1_dp * (t221) * (t597)  &
      * (t1285) - 0.3e1_dp / 0.4e1_dp * (t1995) * (t2715) +  &
      (t603 * t1285) / 0.2e1_dp) + f2716 * t2846 * t233) * alpha3 * t185  &
      * t190 + ((-2._dp * t823 * t2715 + t276 * t1285 + 6._dp * t197 * t319  &
      * t2715 - 2._dp * t197 * t273 * t1285 + t243 * t2875) * alpha4 * t247)  &
      - (2._dp * t91 * t452 * t1285) - (4._dp * t849 * t851) + t3139  &
      * t185 * t190 + (6._dp * t91 * t1839 * t2715) + ((12._dp * t1781  &
      * t2715 - 3._dp * t820 * t1285 - 6._dp * t1786 * t2715 + 2._dp * t823 * t1285  &
      + 2._dp * t243 * t1791 * t2715 - t243 * t826 * t1285 - t280 * t2875) *  &
      alpha6 * t284 * t286) - (12._dp * t1740 * t194 * t44 * t2715)
      e_ndrho_ndrho = e_ndrho_ndrho   -t80 * (t2963 + t3167) * Clda
    END IF

  END SUBROUTINE xwpbe_lda_calc_1

! *****************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega scaling factor
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega!=0 using a taylor
!>      expansion for the parameter G.
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc_2(e_0, e_rho, e_ndrho, e_rho_rho, e_ndrho_rho, &
                              e_ndrho_ndrho, rho, ndrho, omega, sscale, sx, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho, &
                                                e_rho_rho, e_ndrho_rho, &
                                                e_ndrho_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, omega, sscale, sx
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: d2Qndrhondrho, d2Qrhondrho, d2Qrhorho, dQndrho, dQrho, &
      Q, t1, t10, t100, t102, t1022, t1024, t1026, t1029, t103, t1031, t1034, &
      t105, t106, t1061, t1067, t1073, t108, t109, t1090, t1093, t1095, t11, &
      t110, t111, t1111, t1118, t1119, t112, t113, t1136, t1139, t114, t1141, &
      t1149, t115, t1150, t1151, t1157, t1158, t1159, t116, t1162, t1167, &
      t1168, t117, t1178, t118, t1181, t1182, t1186, t1195, t12, t121, t122, &
      t1239, t1241, t1243, t125, t1251, t1256, t126, t1261, t1262, t1263, &
      t1266, t127, t1270, t1273, t1274, t1275, t1276, t128, t1280, t1287, &
      t1288, t129, t13, t1317, t1321
    REAL(KIND=dp) :: t1326, t133, t1331, t1332, t1333, t134, t1341, t1342, &
      t1343, t1347, t1351, t1355, t136, t1362, t137, t1374, t1379, t1382, &
      t1383, t1385, t1392, t1397, t14, t140, t1407, t141, t1417, t1426, &
      t1430, t1434, t1435, t1438, t1442, t1443, t145, t1453, t146, t1469, &
      t1473, t1476, t149, t15, t151, t1517, t1519, t153, t1538, t154, t1545, &
      t1546, t155, t1552, t1553, t156, t1577, t158, t1584, t1589, t159, &
      t1594, t16, t160, t161, t1613, t1614, t1615, t1619, t1626, t163, t1630, &
      t1640, t1644, t1655, t166, t1661, t1666, t1667, t167, t1671, t1678, &
      t1691, t1699, t17, t171, t172, t1726, t173
    REAL(KIND=dp) :: t1732, t1736, t1737, t1742, t1756, t176, t1768, t177, &
      t1773, t1785, t1788, t1791, t1795, t18, t1816, t182, t1824, t184, &
      t1841, t185, t1850, t186, t187, t1885, t19, t190, t191, t1938, t1946, &
      t196, t1976, t1978, t2, t200, t2018, t202, t204, t2053, t2056, t206, &
      t2060, t2066, t2069, t2071, t2072, t2076, t2084, t2086, t209, t2090, &
      t2099, t21, t210, t2105, t2111, t2115, t213, t2133, t2136, t214, t215, &
      t2155, t2158, t216, t218, t2180, t219, t2195, t22, t220, t2203, t221, &
      t2211, t2233, t226, t227, t2274, t2279, t228, t2280, t2283, t23, t230, &
      t2306, t2316, t2323, t233, t2336, t234
    REAL(KIND=dp) :: t236, t2365, t238, t239, t2391, t24, t243, t2432, t2452, &
      t2454, t247, t2473, t248, t2486, t249, t25, t2501, t251, t2511, t2515, &
      t2519, t252, t253, t256, t257, t2571, t258, t2604, t261, t266, t2668, &
      t268, t27, t270, t273, t274, t275, t276, t278, t279, t28, t281, t282, &
      t284, t289, t29, t292, t293, t295, t296, t298, t299, t3, t301, t302, &
      t303, t305, t306, t308, t309, t31, t310, t311, t312, t314, t315, t316, &
      t317, t32, t320, t321, t324, t325, t326, t327, t331, t334, t335, t336, &
      t337, t338, t339, t34, t340, t341, t346, t347, t348, t349, t35, t353, &
      t355, t36, t361, t364, t365
    REAL(KIND=dp) :: t366, t367, t369, t372, t374, t375, t378, t379, t38, &
      t382, t383, t384, t387, t388, t389, t39, t391, t392, t395, t396, t4, &
      t400, t403, t404, t405, t407, t409, t41, t413, t416, t417, t42, t420, &
      t423, t428, t433, t434, t436, t44, t440, t441, t442, t443, t447, t448, &
      t452, t453, t454, t457, t458, t46, t461, t464, t467, t468, t470, t474, &
      t475, t476, t48, t480, t481, t482, t483, t487, t489, t49, t491, t496, &
      t5, t500, t503, t505, t506, t507, t510, t512, t513, t514, t516, t519, &
      t522, t523, t524, t530, t531, t535, t536, t54, t541, t542, t549, t55, &
      t551, t552, t558, t559, t56, t561
    REAL(KIND=dp) :: t565, t566, t569, t574, t577, t579, t58, t583, t584, &
      t585, t587, t588, t591, t595, t596, t6, t60, t603, t604, t605, t607, &
      t608, t61, t612, t615, t620, t622, t626, t628, t629, t63, t633, t634, &
      t635, t638, t641, t644, t65, t650, t652, t656, t658, t659, t660, t665, &
      t67, t670, t671, t679, t68, t681, t685, t687, t689, t69, t692, t695, &
      t697, t699, t7, t70, t702, t705, t709, t71, t710, t714, t72, t723, &
      t725, t727, t73, t733, t736, t737, t739, t74, t740, t742, t747, t748, &
      t75, t751, t752, t755, t756, t758, t759, t760, t761, t767, t769, t77, &
      t770, t78, t781, t783, t784, t788, t793
    REAL(KIND=dp) :: t796, t8, t80, t802, t805, t809, t81, t813, t816, t819, &
      t82, t822, t823, t83, t830, t833, t839, t84, t85, t852, t860, t862, &
      t87, t875, t88, t886, t9, t90, t903, t91, t917, t919, t92, t920, t923, &
      t928, t93, t930, t932, t933, t936, t938, t939, t94, t943, t944, t947, &
      t95, t950, t951, t953, t954, t956, t959, t96, t963, t966, t968, t97, &
      t972, t976, t979, t982, t985, t987, t988, t989, t99, t992

    IF( order >= 0 ) THEN
      t1 = ndrho ** 2
      t2 = r2 ** 2
      t3 = 0.1e1_dp / t2
      t4 = t1 * t3
      t5 = pi ** 2
      t6 = r3 * t5
      t7 = t6 * rho
      t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
      t9 = t8 ** 2
      t10 = 0.1e1_dp / t9
      t11 = t4 * t10
      t12 = rho ** 2
      t13 = 0.1e1_dp / t12
      t14 = sscale ** 2
      t15 = t13 * t14
      t16 = a1 * t1
      t17 = t16 * t3
      t18 = t10 * t13
      t19 = t18 * t14
      t21 = t1 ** 2
      t22 = a2 * t21
      t23 = t2 ** 2
      t24 = 0.1e1_dp / t23
      t25 = t22 * t24
      t27 = 0.1e1_dp / t8 / t7
      t28 = t12 ** 2
      t29 = 0.1e1_dp / t28
      t31 = t14 ** 2
      t32 = t27 * t29 * t31
      t34 = t17 * t19 + t25 * t32
      t35 = a3 * t21
      t36 = t35 * t24
      t38 = t21 * ndrho
      t39 = a4 * t38
      t41 = 0.1e1_dp / t23 / r2
      t42 = t39 * t41
      t44 = 0.1e1_dp / t9 / t7
      t46 = 0.1e1_dp / t28 / rho
      t48 = t31 * sscale
      t49 = t44 * t46 * t48
      t54 = 0.1e1_dp / t23 / t2
      t55 = a5 * t21 * t1 * t54
      t56 = r3 ** 2
      t58 = t5 ** 2
      t60 = 0.1e1_dp / t56 / t58
      t61 = t28 ** 2
      t63 = t31 * t14
      t65 = t60 / t61 * t63
      t67 = r1 + t36 * t32 + t42 * t49 + t55 * t65
      t68 = 0.1e1_dp / t67
      t69 = t34 * t68
      t70 = t15 * t69
      t71 = t11 * t70
      t72 = omega ** 2
      t73 = beta * t72
      t74 = t73 * t10
      t75 = t71 + t74
      t77 = 0.1e1_dp / A
      Q = f94 * t75 * t77
      t78 = rho ** (0.1e1_dp / 0.3e1_dp)
      t80 = t78 * rho * f89
      t81 = B * f12
      t82 = t71 + DD
      t83 = 0.1e1_dp / t82
      t84 = t81 * t83
      t85 = F2 * t34
      t87 = F1 + t85 * t68
      t88 = t15 * t87
      t90 = t11 * t88 + r1
      t91 = f12 * t90
      t92 = t82 ** 2
      t93 = 0.1e1_dp / t92
      t94 = C * t93
      t95 = t91 * t94
      t96 = g2 * t1
      t97 = t96 * t3
      t99 = g3 * t21
      t100 = t99 * t24
      t102 = g1 + t97 * t19 + t100 * t32
      t103 = t15 * t102
      t105 = t11 * t103 + r1
      t106 = t105 * E
      t108 = 0.1e1_dp / t92 / t82
      t109 = t106 * t108
      t110 = f158 * E
      t111 = t105 * t83
      t112 = t72 * t10
      t113 = t71 + DD + t112
      t114 = t113 ** 2
      t115 = t114 ** 2
      t116 = t115 * t113
      t117 = SQRT(t116)
      t118 = 0.1e1_dp / t117
      t121 = SQRT(t113)
      t122 = 0.1e1_dp / t121
      t125 = f68 * C
      t126 = t90 * t83
      t127 = t114 * t113
      t128 = SQRT(t127)
      t129 = 0.1e1_dp / t128
      t133 = (-t110 * t111 * t118 - t81 * t83 * t122 - t125 * t126 * t129)  &
      * omega
      t134 = 0.1e1_dp / t8
      t136 = f52 * E
      t137 = t105 * t93
      t140 = f12 * C
      t141 = t90 * t93
      t145 = t72 * omega
      t146 = (-t136 * t137 * t118 - t140 * t141 * t129) * t145
      t149 = 0.1e1_dp / r3 / t5
      t151 = t149 / rho
      t153 = t72 ** 2
      t154 = t153 * omega
      t155 = t118 * t154
      t156 = t155 * t44
      t158 = f12 * A
      t159 = exei(Q)
      t160 = t71 + DD + t74
      t161 = 0.1e1_dp / t160
      t163 = LOG(t75 * t161)
      t166 = rootpi
      t167 = SQRT(t160)
      t171 = SQRT(A)
      t172 = t171 * f34
      t173 = exer(Q)
      t176 = (t158 * t166 / t167 - t172 * t173) * alpha1
      t177 = omega * t134
      t182 = (t158 * t161 - f98 * t159) * alpha2
      t184 = A * f14
      t185 = t160 ** 2
      t186 = t185 * t160
      t187 = SQRT(t186)
      t190 = SQRT(t75)
      t191 = 0.1e1_dp / t190
      t196 = 0.1e1_dp / t171
      t200 = (t166 * (t184 / t187 - f98 * t191) + f2716 * t173 * t196) *  &
      alpha3 * t145
      t202 = 0.1e1_dp / t75
      t204 = 0.1e1_dp / t185
      t206 = f8132 * t77
      t209 = (-f98 * t202 + t158 * t204 + t206 * t159) * alpha4
      t210 = t153 * t27
      t213 = t75 ** 2
      t214 = t213 * t75
      t215 = SQRT(t214)
      t216 = 0.1e1_dp / t215
      t218 = f38 * A
      t219 = t185 ** 2
      t220 = t219 * t160
      t221 = SQRT(t220)
      t226 = A ** 2
      t227 = t226 * A
      t228 = SQRT(t227)
      t230 = f24364 / t228
      t233 = (t166 * (t206 * t191 - f916 * t216 + t218 / t221) - t230 * t173)  &
      * alpha5
      t234 = t154 * t44
      t236 = 0.1e1_dp / t186
      t238 = 0.1e1_dp / t213
      t239 = f98 * t238
      t243 = f729128 / t226
      t247 = t153 * t72
      t248 = (A * t236 - t239 + t206 * r1 * t202 - t243 * t159) * alpha6  &
      * t247
      t249 = t60 * t13
      t251 = f1516 * A
      t252 = t219 * t186
      t253 = SQRT(t252)
      t256 = t213 ** 2
      t257 = t256 * t75
      t258 = SQRT(t257)
      t261 = f8164 * t77
      t266 = t226 ** 2
      t268 = SQRT(t266 * A)
      t270 = f2187256 / t268
      t273 = (t166 * (t251 / t253 - f2732 / t258 + t261 * t216 - t243 * t191)  &
      + t270 * t173) * alpha7
      t274 = t153 * t145
      t275 = t56 * t58
      t276 = t275 * t12
      t278 = 0.1e1_dp / t8 / t276
      t279 = t274 * t278
      t281 = r3 * A
      t282 = 0.1e1_dp / t219
      t284 = 0.1e1_dp / t214
      t289 = f6561512 / t227
      t292 = (t281 * t282 - f94 * t284 + t206 * t238 - t243 * t202 + t289  &
      * t159) * alpha8
      t293 = t153 ** 2
      t295 = 0.1e1_dp / t9 / t276
      t296 = t293 * t295
      t298 = t84 + t95 + t109 + t133 * t134 + t146 * t151 - t109 * t156 +  &
      t158 * (t159 + t163) + t176 * t177 + t182 * t112 + t200 * t151 +  &
      t209 * t210 + t233 * t234 + t248 * t249 + t273 * t279 + t292 * t296
      t299 = t298 * Clda
      e_0 = e_0 + ( -t80 * t299 ) * sx
    END IF
    IF( order >= 1 .OR. order == -1 ) THEN
      t301 = t44 * t13
      t302 = t4 * t301
      t303 = t14 * t34
      t305 = t68 * r3 * t5
      t306 = t303 * t305
      t308 = 0.2e1_dp / 0.3e1_dp * t302 * t306
      t309 = t12 * rho
      t310 = 0.1e1_dp / t309
      t311 = t310 * t14
      t312 = t311 * t69
      t314 = 2._dp * t11 * t312
      t315 = t3 * t44
      t316 = t16 * t315
      t317 = t15 * t6
      t320 = t10 * t310
      t321 = t320 * t14
      t324 = t24 * t278
      t325 = t22 * t324
      t326 = t29 * t31
      t327 = t326 * t6
      t331 = t27 * t46 * t31
      t334 = -0.2e1_dp / 0.3e1_dp * t316 * t317 - (2._dp * t17 * t321) - 0.4e1_dp  &
      / 0.3e1_dp * t325 * t327 - (4._dp * t25 * t331)
      t335 = t334 * t68
      t336 = t15 * t335
      t337 = t11 * t336
      t338 = t4 * t18
      t339 = t67 ** 2
      t340 = 0.1e1_dp / t339
      t341 = t35 * t324
      t346 = t41 * t295
      t347 = t39 * t346
      t348 = t46 * t48
      t349 = t348 * t6
      t353 = 0.1e1_dp / t28 / t12
      t355 = t44 * t353 * t48
      t361 = t60 / t61 / rho * t63
      t364 = -0.4e1_dp / 0.3e1_dp * t341 * t327 - (4._dp * t36 * t331) - 0.5e1_dp  &
      / 0.3e1_dp * t347 * t349 - (5._dp * t42 * t355) - (8._dp * t55 * t361)
      t365 = t340 * t364
      t366 = t303 * t365
      t367 = t338 * t366
      t369 = t44 * r3 * t5
      t372 = -t308 - t314 + t337 - t367 - 0.2e1_dp / 0.3e1_dp * t73 * t369
      dQrho = f94 * t372 * t77
      t374 = ndrho * t3
      t375 = t374 * t10
      t378 = a1 * ndrho
      t379 = t378 * t3
      t382 = t1 * ndrho
      t383 = a2 * t382
      t384 = t383 * t24
      t387 = 2._dp * t379 * t19 + 4._dp * t384 * t32
      t388 = t387 * t68
      t389 = t15 * t388
      t391 = a3 * t382
      t392 = t391 * t24
      t395 = a4 * t21
      t396 = t395 * t41
      t400 = a5 * t38 * t54
      t403 = 4._dp * t392 * t32 + 5._dp * t396 * t49 + 6._dp * t400 * t65
      t404 = t340 * t403
      t405 = t303 * t404
      t407 = 2._dp * t375 * t70 + t11 * t389 - t338 * t405
      dQndrho = f94 * t407 * t77
      t409 = t78 * f89
      t413 = t27 * r3 * t5
      t416 = t14 * t102
      t417 = t416 * t6
      t420 = t311 * t102
      t423 = t96 * t315
      t428 = t99 * t324
      t433 = -0.2e1_dp / 0.3e1_dp * t423 * t317 - (2._dp * t97 * t321) - 0.4e1_dp  &
      / 0.3e1_dp * t428 * t327 - (4._dp * t100 * t331)
      t434 = t15 * t433
      t436 = -0.2e1_dp / 0.3e1_dp * t302 * t417 - (2._dp * t11 * t420) + (t11  &
      * t434)
      t440 = t136 * t105
      t441 = t108 * t118
      t442 = -t308 - t314 + t337 - t367
      t443 = t441 * t442
      t447 = 0.1e1_dp / t117 / t116
      t448 = t93 * t447
      t452 = -t308 - t314 + t337 - t367 - 0.2e1_dp / 0.3e1_dp * t72 * t44 * t6
      t453 = t115 * t452
      t454 = t448 * t453
      t457 = t14 * t87
      t458 = t457 * t6
      t461 = t311 * t87
      t464 = F2 * t334
      t467 = t464 * t68 - t85 * t365
      t468 = t15 * t467
      t470 = -0.2e1_dp / 0.3e1_dp * t302 * t458 - (2._dp * t11 * t461) + (t11  &
      * t468)
      t474 = t140 * t90
      t475 = t108 * t129
      t476 = t475 * t442
      t480 = 0.1e1_dp / t128 / t127
      t481 = t93 * t480
      t482 = t114 * t452
      t483 = t481 * t482
      t487 = (-t136 * t436 * t93 * t118 + (2._dp * t440 * t443) + 0.5e1_dp /  &
      0.2e1_dp * (t440) * (t454) - t140 * t470 * t93 * t129 + (2._dp  &
      * t474 * t476) + 0.3e1_dp / 0.2e1_dp * (t474) * (t483)) * t145
      t489 = t209 * t153
      t491 = t278 * r3 * t5
      t496 = t166 / t167 / t160
      t500 = dexerrho(Q,dQrho)
      t503 = (-t158 * t496 * t372 / 0.2e1_dp - t172 * t500) * alpha1
      t505 = t106 * t441
      t506 = t154 * t295
      t507 = t506 * t6
      t510 = f12 * t470
      t512 = t92 ** 2
      t513 = 0.1e1_dp / t512
      t514 = t106 * t513
      t516 = t155 * t44 * t442
      t519 = t149 * t13
      t522 = 0.1e1_dp / t253 / t252
      t523 = t219 * t185
      t524 = t522 * t523
      t530 = f2732 / t258 / t257
      t531 = t256 * t372
      t535 = 0.1e1_dp / t215 / t214
      t536 = t535 * t213
      t541 = 0.1e1_dp / t190 / t75
      t542 = t541 * t372
      t549 = (t166 * (-0.7e1_dp / 0.2e1_dp * t251 * t524 * t372 + 0.5e1_dp / 0.2e1_dp  &
      * t530 * t531 - 0.3e1_dp / 0.2e1_dp * t261 * t536 * t372 + t243 * t542  &
      / 0.2e1_dp) + t270 * t500) * alpha7
      t551 = C * t108
      t552 = t551 * t442
      t558 = t436 * E
      t559 = t558 * t108
      t561 = 0.1e1_dp / t220
      t565 = 0.1e1_dp / t256
      t566 = f94 * t565
      t569 = t284 * t372
      t574 = dexeirho(Q,dQrho)
      t577 = (-4._dp * t281 * t561 * t372 + 3._dp * t566 * t372 - 2._dp * t206 * t569  &
      + t243 * t238 * t372 + t289 * t574) * alpha8
      t579 = -t133 * t413 / 0.3e1_dp + t487 * t151 - 0.4e1_dp / 0.3e1_dp * t489 *  &
      t491 + t503 * t177 + 0.5e1_dp / 0.3e1_dp * t505 * t507 + t510 * t94 + (3._dp  &
      * t514 * t516) - t146 * t519 + t549 * t279 - (2._dp * t91 * t552)  &
      - t81 * t93 * t442 - t200 * t519 - t559 * t156 + t577 * t296  &
      + t559
      t583 = t110 * t105
      t584 = t93 * t118
      t585 = t584 * t442
      t587 = t83 * t447
      t588 = t587 * t453
      t591 = t93 * t122
      t595 = 0.1e1_dp / t121 / t113
      t596 = t83 * t595
      t603 = t125 * t90
      t604 = t93 * t129
      t605 = t604 * t442
      t607 = t83 * t480
      t608 = t607 * t482
      t612 = (-t110 * t436 * t83 * t118 + t583 * t585 + 0.5e1_dp / 0.2e1_dp * t583  &
      * t588 + t81 * t591 * t442 + t81 * t596 * t452 / 0.2e1_dp - t125  &
      * t470 * t83 * t129 + t603 * t605 + 0.3e1_dp / 0.2e1_dp * t603 * t608) *  &
      omega
      t615 = t236 * t372
      t620 = (t239 * t372 - 2._dp * t158 * t615 + t206 * t574) * alpha4
      t622 = t513 * t442
      t626 = t75 * t204
      t628 = t372 * t161 - t626 * t372
      t629 = t628 * t202
      t633 = t233 * t154
      t634 = t295 * r3
      t635 = t634 * t5
      t638 = A * t282
      t641 = f98 * t284
      t644 = r1 * t238
      t650 = (-3._dp * t638 * t372 + 2._dp * t641 * t372 - t206 * t644 * t372 - t243  &
      * t574) * alpha6 * t247
      t652 = t204 * t372
      t656 = (-t158 * t652 - f98 * t574) * alpha2
      t658 = t108 * t447
      t659 = t106 * t658
      t660 = t234 * t453
      t665 = f916 * t535
      t670 = 0.1e1_dp / t221 / t220
      t671 = t670 * t219
      t679 = (t166 * (-t206 * t542 / 0.2e1_dp + 0.3e1_dp / 0.2e1_dp * t665 * t213  &
      * t372 - 0.5e1_dp / 0.2e1_dp * t218 * t671 * t372) - t230 * t500) * alpha5
      t681 = t292 * t293
      t685 = t56 * r3 * t58 * t5 * t309
      t687 = 0.1e1_dp / t9 / t685
      t689 = t687 * r3 * t5
      t692 = t60 * t310
      t695 = t273 * t274
      t697 = 0.1e1_dp / t8 / t685
      t699 = t697 * r3 * t5
      t702 = t176 * omega
      t705 = t182 * t72
      t709 = 0.1e1_dp / t187 / t186
      t710 = t709 * t185
      t714 = f98 * t541
      t723 = (t166 * (-0.3e1_dp / 0.2e1_dp * t184 * t710 * t372 + t714 * t372 /  &
      0.2e1_dp) + f2716 * t500 * t196) * alpha3 * t145
      t725 = t612 * t134 + t620 * t210 - (3._dp * t106 * t622) + t158 * (t574  &
      + t629 * t160) - 0.5e1_dp / 0.3e1_dp * t633 * t635 + t650 * t249 +  &
      t656 * t112 + 0.5e1_dp / 0.2e1_dp * t659 * t660 + t679 * t234 - 0.8e1_dp /  &
      0.3e1_dp * t681 * t689 - (2._dp * t248 * t692) - 0.7e1_dp / 0.3e1_dp * t695  &
      * t699 - t702 * t413 / 0.3e1_dp - 0.2e1_dp / 0.3e1_dp * t705 * t369 + t723  &
      * t151
      t727 = (t579 + t725) * Clda
      e_rho = e_rho + ( -0.4e1_dp / 0.3e1_dp * t409 * t299 - t80 * t727 ) * sx
      t733 = F2 * t387
      t736 = t733 * t68 - t85 * t404
      t737 = t15 * t736
      t739 = 2._dp * t375 * t88 + t11 * t737
      t740 = f12 * t739
      t742 = t551 * t407
      t747 = g2 * ndrho
      t748 = t747 * t3
      t751 = g3 * t382
      t752 = t751 * t24
      t755 = 2._dp * t748 * t19 + 4._dp * t752 * t32
      t756 = t15 * t755
      t758 = 2._dp * t375 * t103 + t11 * t756
      t759 = t758 * E
      t760 = t759 * t108
      t761 = t513 * t407
      t767 = t584 * t407
      t769 = t115 * t407
      t770 = t587 * t769
      t781 = t604 * t407
      t783 = t114 * t407
      t784 = t607 * t783
      t788 = (-t110 * t758 * t83 * t118 + t583 * t767 + 0.5e1_dp / 0.2e1_dp * t583  &
      * t770 + t81 * t591 * t407 + t81 * t596 * t407 / 0.2e1_dp - t125  &
      * t739 * t83 * t129 + t603 * t781 + 0.3e1_dp / 0.2e1_dp * t603 * t784) *  &
      omega
      t793 = t441 * t407
      t796 = t448 * t769
      t802 = t475 * t407
      t805 = t481 * t783
      t809 = (-t136 * t758 * t93 * t118 + (2._dp * t440 * t793) + 0.5e1_dp /  &
      0.2e1_dp * (t440) * (t796) - t140 * t739 * t93 * t129 + (2._dp  &
      * t474 * t802) + 0.3e1_dp / 0.2e1_dp * (t474) * (t805)) * t145
      t813 = t155 * t44 * t407
      t816 = t234 * t769
      t819 = dexeindrho(Q,dQndrho)
      t822 = t407 * t161 - t626 * t407
      t823 = t822 * t202
      t830 = dexerndrho(Q,dQndrho)
      t833 = (-t158 * t496 * t407 / 0.2e1_dp - t172 * t830) * alpha1
      t839 = (-t158 * t204 * t407 - f98 * t819) * alpha2
      t852 = (t166 * (-0.3e1_dp / 0.2e1_dp * t184 * t710 * t407 + t714 * t407 /  &
      0.2e1_dp) + f2716 * t830 * t196) * alpha3 * t145
      t860 = (t239 * t407 - 2._dp * t158 * t236 * t407 + t206 * t819) * alpha4
      t862 = t541 * t407
      t875 = (t166 * (-t206 * t862 / 0.2e1_dp + 0.3e1_dp / 0.2e1_dp * t665 * t213  &
      * t407 - 0.5e1_dp / 0.2e1_dp * t218 * t671 * t407) - t230 * t830) * alpha5
      t886 = (-3._dp * t638 * t407 + 2._dp * t641 * t407 - t206 * t644 * t407 - t243  &
      * t819) * alpha6 * t247
      t903 = (t166 * (-0.7e1_dp / 0.2e1_dp * t251 * t524 * t407 + 0.5e1_dp / 0.2e1_dp  &
      * t530 * t256 * t407 - 0.3e1_dp / 0.2e1_dp * t261 * t536 * t407 + t243  &
      * t862 / 0.2e1_dp) + t270 * t830) * alpha7
      t917 = (-4._dp * t281 * t561 * t407 + 3._dp * t566 * t407 - 2._dp * t206 * t284  &
      * t407 + t243 * t238 * t407 + t289 * t819) * alpha8
      t919 = -t81 * t93 * t407 + t740 * t94 - (2._dp * t91 * t742) + t760  &
      - (3._dp * t106 * t761) + t788 * t134 + t809 * t151 - t760 * t156  &
      + (3._dp * t514 * t813) + 0.5e1_dp / 0.2e1_dp * t659 * t816 + t158 * (t819  &
      + t823 * t160) + t833 * t177 + t839 * t112 + t852 * t151 + t860  &
      * t210 + t875 * t234 + t886 * t249 + t903 * t279 + t917 * t296
      t920 = t919 * Clda
      e_ndrho = e_ndrho + ( -t80 * t920 ) * sx
    END IF
    IF( order >= 2 .OR. order == -2 ) THEN
      t923 = t4 * t295 * t13
      t928 = 0.10e2_dp / 0.9e1_dp * t923 * t303 * t68 * t56 * t58
      t930 = t4 * t44 * t310
      t932 = 0.8e1_dp / 0.3e1_dp * t930 * t306
      t933 = t14 * t334
      t936 = 0.4e1_dp / 0.3e1_dp * t302 * t933 * t305
      t938 = t4 * t301 * t14
      t939 = t34 * t340
      t943 = 0.4e1_dp / 0.3e1_dp * t938 * t939 * t6 * t364
      t944 = t29 * t14
      t947 = 6._dp * t11 * t944 * t69
      t950 = 4._dp * t11 * t311 * t335
      t951 = t4 * t320
      t953 = 4._dp * t951 * t366
      t954 = t3 * t295
      t956 = t15 * t275
      t959 = t311 * t6
      t963 = t10 * t29 * t14
      t966 = t24 * t697
      t968 = t326 * t275
      t972 = t46 * t31 * t6
      t976 = t27 * t353 * t31
      t979 = 0.10e2_dp / 0.9e1_dp * t16 * t954 * t956 + 0.8e1_dp / 0.3e1_dp * t316 *  &
      t959 + (6._dp * t17 * t963) + 0.28e2_dp / 0.9e1_dp * t22 * t966 * t968 +  &
      0.32e2_dp / 0.3e1_dp * t325 * t972 + (20._dp * t25 * t976)
      t982 = t11 * t15 * t979 * t68
      t985 = 2._dp * t338 * t933 * t365
      t987 = 0.1e1_dp / t339 / t67
      t988 = t364 ** 2
      t989 = t987 * t988
      t992 = 2._dp * t338 * t303 * t989
      t1022 = t340 * (0.28e2_dp / 0.9e1_dp * t35 * t966 * t968 + 0.32e2_dp / 0.3e1_dp  &
      * t341 * t972 + (20._dp * t36 * t976) + 0.40e2_dp / 0.9e1_dp * t39 * t41  &
      * t687 * t348 * t275 + 0.50e2_dp / 0.3e1_dp * t347 * t353 * t48 * t6 +  &
      0.30e2_dp * t42 * t44 / t28 / t309 * t48 + (72._dp * t55 * t60 / t61  &
      / t12 * t63))
      t1024 = t338 * t303 * t1022
      t1026 = t295 * t56 * t58
      t1029 = t928 + t932 - t936 + t943 + t947 - t950 + t953 + t982 - t985  &
      + t992 - t1024 + 0.10e2_dp / 0.9e1_dp * t73 * t1026
      d2Qrhorho = f94 * t1029 * t77
      t1031 = t374 * t301
      t1034 = t14 * t387
      t1061 = -0.4e1_dp / 0.3e1_dp * t378 * t315 * t317 - (4._dp * t379 * t321)  &
      - 0.16e2_dp / 0.3e1_dp * t383 * t324 * t327 - (16._dp * t384 * t331)
      t1067 = t374 * t18
      t1073 = t987 * t364 * t403
      t1090 = t340 * (-0.16e2_dp / 0.3e1_dp * t391 * t324 * t327 - (16._dp * t392  &
      * t331) - 0.25e2_dp / 0.3e1_dp * t395 * t346 * t349 - (25._dp * t396  &
      * t355) - (48._dp * t400 * t361))
      t1093 = -0.4e1_dp / 0.3e1_dp * t1031 * t306 - 0.2e1_dp / 0.3e1_dp * t302 * t1034  &
      * t305 + 0.2e1_dp / 0.3e1_dp * t938 * t939 * t6 * t403 - (4._dp * t375  &
      * t312) - (2._dp * t11 * t311 * t388) + (2._dp * t951 * t405) + (2._dp  &
      * t375 * t336) + (t11 * t15 * t1061 * t68) - t338 * t933  &
      * t404 - (2._dp * t1067 * t366) - t338 * t1034 * t365 + 0.2e1_dp * t338  &
      * t303 * t1073 - t338 * t303 * t1090
      d2Qrhondrho = f94 * t1093 * t77
      t1095 = t3 * t10
      t1111 = 2._dp * a1 * t3 * t19 + 12._dp * a2 * t1 * t24 * t32
      t1118 = t403 ** 2
      t1119 = t987 * t1118
      t1136 = t340 * (12._dp * a3 * t1 * t24 * t32 + 20._dp * a4 * t382 * t41 * t49  &
      + 30._dp * a5 * t21 * t54 * t65)
      t1139 = 2._dp * t1095 * t13 * t303 * t68 + 4._dp * t375 * t389 - 4._dp * t1067  &
      * t405 + t11 * t15 * t1111 * t68 - 2._dp * t338 * t1034 * t404 + 2._dp * t338  &
      * t303 * t1119 - t338 * t303 * t1136
      d2Qndrhondrho = f94 * t1139 * t77
      t1141 = t78 ** 2
      t1149 = 0.1e1_dp / t512 / t82
      t1150 = t106 * t1149
      t1151 = t442 ** 2
      t1157 = 0.1e1_dp / t190 / t213
      t1158 = t372 ** 2
      t1159 = t1157 * t1158
      t1162 = t541 * t1029
      t1167 = 0.1e1_dp / t215 / t256 / t213
      t1168 = f916 * t1167
      t1178 = t219 ** 2
      t1181 = 0.1e1_dp / t221 / t1178 / t185
      t1182 = t1181 * t1178
      t1186 = t670 * t186
      t1195 = d2exerrhorho(Q,dQrho,d2Qrhorho)
      t1239 = 0.10e2_dp / 0.9e1_dp * t923 * t416 * t275 + 0.8e1_dp / 0.3e1_dp * t930  &
      * t417 - 0.4e1_dp / 0.3e1_dp * t302 * t14 * t433 * t6 + (6._dp * t11 * t944  &
      * t102) - 0.4e1_dp * (t11) * t311 * t433 + (t11) * t15 *  &
      (0.10e2_dp / 0.9e1_dp * t96 * t954 * t956 + 0.8e1_dp / 0.3e1_dp * t423 * t959  &
      + (6._dp * t97 * t963) + 0.28e2_dp / 0.9e1_dp * t99 * t966 * t968 + 0.32e2_dp  &
      / 0.3e1_dp * t428 * t972 + (20._dp * t100 * t976))
      t1241 = t1239 * E * t108
      t1243 = t278 * t56 * t58
      t1251 = t136 * t436
      t1256 = t513 * t118
      t1261 = t136 * t105 * t108
      t1262 = t447 * t442
      t1263 = t1262 * t453
      t1266 = t928 + t932 - t936 + t943 + t947 - t950 + t953 + t982 - t985  &
      + t992 - t1024
      t1270 = t115 ** 2
      t1273 = 0.1e1_dp / t117 / t1270 / t114
      t1274 = t93 * t1273
      t1275 = t452 ** 2
      t1276 = t1270 * t1275
      t1280 = t127 * t1275
      t1287 = t928 + t932 - t936 + t943 + t947 - t950 + t953 + t982 - t985  &
      + t992 - t1024 + 0.10e2_dp / 0.9e1_dp * t72 * t295 * t275
      t1288 = t115 * t1287
      t1317 = 0.10e2_dp / 0.9e1_dp * t923 * t457 * t275 + 0.8e1_dp / 0.3e1_dp * t930  &
      * t458 - 0.4e1_dp / 0.3e1_dp * t302 * t14 * t467 * t6 + (6._dp * t11 * t944  &
      * t87) - 0.4e1_dp * (t11) * t311 * t467 + (t11 * t15 * (F2  &
      * t979 * t68 - 2._dp * t464 * t365 + 2._dp * t85 * t989 - t85 * t1022))
      t1321 = t140 * t470
      t1326 = t513 * t129
      t1331 = t140 * t90 * t108
      t1332 = t480 * t442
      t1333 = t1332 * t482
      t1341 = 0.1e1_dp / t128 / t115 / t114
      t1342 = t93 * t1341
      t1343 = t115 * t1275
      t1347 = t113 * t1275
      t1351 = t114 * t1287
      t1355 = -t136 * t1239 * t93 * t118 + (4._dp * t1251 * t443) + (5._dp  &
      * t1251 * t454) - (6._dp * t440 * t1256 * t1151) - (10._dp * t1261  &
      * t1263) + (2._dp * t440 * t441 * t1266) - 0.75e2_dp / 0.4e1_dp * (t440)  &
      * (t1274) * (t1276) + (10._dp * t440 * t448 * t1280)  &
      + 0.5e1_dp / 0.2e1_dp * (t440) * (t448) * (t1288) - t140  &
      * t1317 * t93 * t129 + (4._dp * t1321 * t476) + (3._dp * t1321 * t483)  &
      - (6._dp * t474 * t1326 * t1151) - (6._dp * t1331 * t1333) +  &
      (2._dp * t474 * t475 * t1266) - 0.27e2_dp / 0.4e1_dp * (t474) * (t1342)  &
      * (t1343) + (3._dp * t474 * t481 * t1347) + 0.3e1_dp / 0.2e1_dp  &
      * (t474) * (t481) * (t1351)
      t1362 = t106 * t658 * t154
      t1374 = d2exeirhorho(Q,dQrho,d2Qrhorho)
      t1379 = t558 * t658
      t1382 = t56 ** 2
      t1383 = t58 ** 2
      t1385 = t1382 * t1383 * t28
      t1392 = -(12._dp * t1150 * t155 * t44 * t1151) + (t166 * (0.3e1_dp / 0.4e1_dp  &
      * t206 * t1159 - t206 * t1162 / 0.2e1_dp - 0.27e2_dp / 0.4e1_dp * t1168  &
      * t256 * t1158 + 0.3e1_dp * t665 * t75 * t1158 + 0.3e1_dp / 0.2e1_dp * t665  &
      * t213 * t1029 + 0.75e2_dp / 0.4e1_dp * t218 * t1182 * t1158 - 0.10e2_dp  &
      * t218 * t1186 * t1158 - 0.5e1_dp / 0.2e1_dp * t218 * t671 * t1029) - t230  &
      * t1195) * alpha5 * t234 + 0.28e2_dp / 0.9e1_dp * t489 * t697 * t56  &
      * t58 + 0.10e2_dp / 0.3e1_dp * t558 * t441 * t507 + t1241 + 0.4e1_dp / 0.9e1_dp  &
      * t702 * t1243 + 0.4e1_dp / 0.9e1_dp * t133 * t1243 + t1355 * t145 * t151  &
      - (2._dp * t91 * t551 * t1266) - 0.25e2_dp / 0.3e1_dp * t1362 * t295  &
      * t115 * t452 * r3 * t5 + (0.2e1_dp * t158 * t236 * t1158 - t158 * t204  &
      * t1029 - f98 * t1374) * alpha2 * t112 + (5._dp * t1379 * t660)  &
      + 0.70e2_dp / 0.9e1_dp * t695 / t8 / t1385 * t56 * t58
      t1397 = t106 * t108 * t1273
      t1407 = t110 * t436
      t1417 = t110 * t137
      t1426 = t83 * t1273
      t1430 = t108 * t122
      t1434 = t81 * t93
      t1435 = t595 * t442
      t1438 = -t110 * t1239 * t83 * t118 + (2._dp * t1407 * t585) + (5._dp  &
      * t1407 * t588) - (2._dp * t583 * t441 * t1151) + (t583 * t584  &
      * t1266) - (5._dp * t1417 * t1263) + (10._dp * t583 * t587 * t1280)  &
      + 0.5e1_dp / 0.2e1_dp * (t583) * (t587) * (t1288) - 0.75e2_dp  &
      / 0.4e1_dp * (t583) * (t1426) * (t1276) - (2._dp * t81  &
      * t1430 * t1151) - t1434 * t1435 * t452
      t1442 = 0.1e1_dp / t121 / t114
      t1443 = t83 * t1442
      t1453 = t125 * t470
      t1469 = t83 * t1341
      t1473 = t125 * t141
      t1476 = t81 * t591 * (t1266) - 0.3e1_dp / 0.4e1_dp * t81 * t1443 * t1275  &
      + t81 * t596 * t1287 / 0.2e1_dp - t125 * t1317 * t83 * t129 + (2._dp  &
      * t1453 * t605) + (3._dp * t1453 * t608) - (2._dp * t603 * t475  &
      * t1151) + (t603 * t604 * t1266) + (3._dp * t603 * t607 * t1347)  &
      + 0.3e1_dp / 0.2e1_dp * (t603) * (t607) * (t1351) - 0.27e2_dp  &
      / 0.4e1_dp * (t603) * (t1469) * (t1343) - (3._dp * t1473  &
      * t1333)
      t1517 = -0.16e2_dp / 0.3e1_dp * t577 * t293 * t689 - 0.75e2_dp / 0.4e1_dp * t1397  &
      * t234 * t1276 + 0.5e1_dp / 0.2e1_dp * t659 * t234 * t1288 + (t1438 +  &
      t1476) * omega * t134 + f12 * t1317 * t94 - (4._dp * t650 * t692)  &
      - 0.14e2_dp / 0.3e1_dp * t549 * t274 * t699 + (12._dp * t106 * t1149 *  &
      t1151) - (3._dp * t106 * t513 * t1266) - 0.40e2_dp / 0.9e1_dp * t505 * t154  &
      * t687 * t275 + ((-2._dp * t641 * t1158 + t239 * t1029 + 6._dp * t158  &
      * t282 * t1158 - 2._dp * t158 * t236 * t1029 + t206 * t1374) * alpha4  &
      * t210) - (6._dp * t558 * t622) + (6._dp * t248 * t60 * t29) -  &
      (4._dp * t510 * t552)
      t1519 = t149 * t310
      t1538 = t75 * t236
      t1545 = t628 * t238
      t1546 = t160 * t372
      t1552 = 0.1e1_dp / t167 / t185
      t1553 = t166 * t1552
      t1577 = (2._dp * t146 * t1519) - (2._dp * t723 * t519) + 0.10e2_dp / 0.9e1_dp  &
      * t705 * t1026 - t1241 * t156 - 0.10e2_dp / 0.3e1_dp * t679 * t154  &
      * t635 - 0.2e1_dp / 0.3e1_dp * t503 * omega * t413 + (2._dp * t200 * t1519)  &
      + (t158 * (t1374 + (t1029 * t161 - 2._dp * t1158 * t204 + 2._dp *  &
      t1538 * t1158 - t626 * t1029) * t202 * t160 - t1545 * t1546 + t629  &
      * t372)) + (0.3e1_dp / 0.4e1_dp * (t158) * (t1553) * (t1158)  &
      - (t158 * t496 * t1029) / 0.2e1_dp - t172 * t1195) * alpha1 * t177  &
      - 0.8e1_dp / 0.3e1_dp * t620 * t153 * t491 - t81 * t93 * t1266 - 0.2e1_dp  &
      / 0.3e1_dp * t612 * t413 + 0.88e2_dp / 0.9e1_dp * t681 / t9 / t1385 * t56  &
      * t58
      t1584 = A * t561
      t1589 = f98 * t565
      t1594 = r1 * t284
      t1613 = 0.1e1_dp / t253 / t1178 / t523
      t1614 = t1178 * t219
      t1615 = t1613 * t1614
      t1619 = t522 * t220
      t1626 = t256 ** 2
      t1630 = f2732 / t258 / t1626 / t213
      t1640 = t1167 * t256
      t1644 = t535 * t75
      t1655 = 0.147e3_dp / 0.4e1_dp * t251 * t1615 * t1158 - 0.21e2_dp * t251 * t1619  &
      * t1158 - 0.7e1_dp / 0.2e1_dp * t251 * t524 * t1029 - 0.75e2_dp / 0.4e1_dp  &
      * t1630 * t1626 * t1158 + 0.10e2_dp * t530 * t214 * t1158 + 0.5e1_dp /  &
      0.2e1_dp * t530 * t256 * t1029 + 0.27e2_dp / 0.4e1_dp * t261 * t1640 * t1158  &
      - 0.3e1_dp * t261 * t1644 * t1158 - 0.3e1_dp / 0.2e1_dp * t261 * t536 * t1029  &
      - 0.3e1_dp / 0.4e1_dp * t243 * t1159 + t243 * t1162 / 0.2e1_dp
      t1661 = C * t513
      t1666 = 0.1e1_dp / t187 / t523
      t1667 = t1666 * t219
      t1671 = t709 * t160
      t1678 = f98 * t1157
      t1691 = 0.1e1_dp / t523
      t1699 = f94 / t257
      t1726 = t106 * t1256
      t1732 = t558 * t513
      t1736 = t106 * t513 * t447
      t1737 = t442 * t115
      t1742 = -(2._dp * t487 * t519) + 0.40e2_dp / 0.9e1_dp * t633 * t687 * t56  &
      * t58 + ((12._dp * t1584 * t1158 - 3._dp * t638 * t1029 - 6._dp * t1589 *  &
      t1158 + 2._dp * t641 * t1029 + 2._dp * t206 * t1594 * t1158 - t206 * t644  &
      * t1029 - t243 * t1374) * alpha6 * t247 * t249) - 0.4e1_dp / 0.3e1_dp *  &
      t656 * t72 * t369 + (2._dp * t81 * t108 * t1151) + (t166 * t1655  &
      + t270 * t1195) * alpha7 * t279 + (6._dp * t91 * t1661 * t1151) +  &
      (t166 * (0.27e2_dp / 0.4e1_dp * (t184) * (t1667) * (t1158) -  &
      (3._dp * t184 * t1671 * t1158) - 0.3e1_dp / 0.2e1_dp * (t184) * (t710)  &
      * (t1029) - 0.3e1_dp / 0.4e1_dp * (t1678) * (t1158)  &
      + (t714 * t1029) / 0.2e1_dp) + f2716 * t1195 * t196) * alpha3 * t145  &
      * t151 + ((20._dp * t281 * t1691 * t1158 - 4._dp * t281 * t561 * t1029  &
      - 12._dp * t1699 * t1158 + 3._dp * t566 * t1029 + 6._dp * t206 * t565 * t1158  &
      - 2._dp * t206 * t284 * t1029 - 2._dp * t243 * t284 * t1158 + t243 *  &
      t238 * t1029 + t289 * t1374) * alpha8 * t296) + (10._dp * t659 * t234  &
      * t1280) + (3._dp * t514 * t155 * t44 * t1266) - (10._dp * t1726  &
      * t506 * t442 * r3 * t5) + (6._dp * t1732 * t516) - (15._dp * t1736  &
      * t234 * t1737 * t452)
      e_rho_rho = e_rho_rho + ( -0.4e1_dp / 0.9e1_dp / t1141 * f89 * t299 - 0.8e1_dp / 0.3e1_dp * t409 *  &
      t727 - t80 * (t1392 + t1517 + t1577 + t1742) * Clda ) * sx
      t1756 = t372 * t407
      t1768 = t569 * t407
      t1773 = d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t1785 = t1157 * t372 * t407
      t1788 = t541 * t1093
      t1791 = t531 * t407
      t1795 = t75 * t372 * t407
      t1816 = d2exerrhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t1824 = t759 * t658
      t1841 = t442 * t407
      t1850 = ((20._dp * t281 * t1691 * t372 * t407 - 4._dp * t281 * t561 * t1093  &
      - 12._dp * t1699 * t1756 + 3._dp * t566 * t1093 + 6._dp * t206 * t565 * t372  &
      * t407 - 2._dp * t206 * t284 * t1093 - 2._dp * t243 * t1768 + t243 * t238  &
      * t1093 + t289 * t1773) * alpha8 * t296) - t833 * omega * t413  &
      / 0.3e1_dp - 0.2e1_dp / 0.3e1_dp * t839 * t72 * t369 + (t166 * (0.3e1_dp / 0.4e1_dp  &
      * (t206) * (t1785) - (t206 * t1788) / 0.2e1_dp - 0.27e2_dp  &
      / 0.4e1_dp * t1168 * t1791 + (3._dp * t665 * t1795) + 0.3e1_dp / 0.2e1_dp  &
      * (t665) * (t213) * (t1093) + 0.75e2_dp / 0.4e1_dp * (t218)  &
      * (t1181) * (t1178) * (t372) * (t407) - (10._dp  &
      * t218 * t670 * t186 * t372 * t407) - 0.5e1_dp / 0.2e1_dp * (t218)  &
      * (t671) * (t1093)) - t230 * t1816) * alpha5 * t234 - (3._dp  &
      * t106 * t513 * t1093) + 0.5e1_dp / 0.2e1_dp * t1824 * t660 - (2._dp  &
      * t740 * t552) - 0.8e1_dp / 0.3e1_dp * t917 * t293 * t689 - 0.7e1_dp / 0.3e1_dp  &
      * t903 * t274 * t699 - 0.15e2_dp / 0.2e1_dp * t1736 * t234 * t453 *  &
      (t407) - 0.12e2_dp * (t106) * t1149 * t118 * t234 * t1841 +  &
      0.10e2_dp * t659 * t234 * t127 * t452 * (t407)
      t1885 = t160 * t407
      t1938 = -0.5e1_dp / 0.3e1_dp * t875 * t154 * t635 + ((12._dp * t1584 * t1756  &
      - 3._dp * t638 * t1093 - 6._dp * t1589 * t1756 + 2._dp * t641 * t1093 + 2._dp  &
      * t206 * r1 * t1768 - t206 * t644 * t1093 - t243 * t1773) * alpha6  &
      * t247 * t249) - (2._dp * t510 * t742) - (t81 * t93 * t1093)  &
      + (t158 * (t1773 + (t1093 * t161 - 2._dp * t652 * t407 + 2._dp * t1538  &
      * t1756 - t626 * t1093) * t202 * t160 - t1545 * t1885 + t629 * t407))  &
      + (0.3e1_dp / 0.4e1_dp * (t158) * (t166) * (t1552) * (t372)  &
      * (t407) - (t158 * t496 * t1093) / 0.2e1_dp - t172 *  &
      t1816) * alpha1 * t177 - (3._dp * t759 * t622) - 0.15e2_dp / 0.2e1_dp *  &
      (t1736) * (t234) * (t1737) * (t407) + (3._dp * t514  &
      * t155 * t44 * t1093) + ((2._dp * t158 * t615 * t407 - t158 * t204  &
      * t1093 - f98 * t1773) * alpha2 * t112) + ((-2._dp * t641 * t1756  &
      + t239 * t1093 + 6._dp * t158 * t282 * t372 * t407 - 2._dp * t158 * t236  &
      * t1093 + t206 * t1773) * alpha4 * t210) + (6._dp * t474 * t622  &
      * t407)
      t1946 = t115 * t1093
      t1976 = -0.4e1_dp / 0.3e1_dp * t1031 * t417 - 0.2e1_dp / 0.3e1_dp * t302 * t14  &
      * t755 * t6 - (4._dp * t375 * t420) - 0.2e1_dp * t11 * t311 * t755 +  &
      (2._dp * t375 * t434) + t11 * t15 * (-0.4e1_dp / 0.3e1_dp * t747 * t315  &
      * t317 - (4._dp * t748 * t321) - 0.16e2_dp / 0.3e1_dp * t751 * t324 * t327  &
      - (16._dp * t752 * t331))
      t1978 = t1976 * E * t108
      t2018 = 0.147e3_dp / 0.4e1_dp * t251 * t1613 * t1614 * t372 * t407 - 0.21e2_dp  &
      * t251 * t522 * t220 * t372 * t407 - 0.7e1_dp / 0.2e1_dp * t251 * t524  &
      * t1093 - 0.75e2_dp / 0.4e1_dp * t1630 * t1626 * t372 * t407 + 0.10e2_dp  &
      * t530 * t214 * t372 * t407 + 0.5e1_dp / 0.2e1_dp * t530 * t256 * t1093  &
      + 0.27e2_dp / 0.4e1_dp * t261 * t1167 * t1791 - 0.3e1_dp * t261 * t535 * t1795  &
      - 0.3e1_dp / 0.2e1_dp * t261 * t536 * t1093 - 0.3e1_dp / 0.4e1_dp * t243 *  &
      t1785 + t243 * t1788 / 0.2e1_dp
      t2053 = -0.4e1_dp / 0.3e1_dp * t1031 * t458 - 0.2e1_dp / 0.3e1_dp * t302 * t14  &
      * t736 * t6 - (4._dp * t375 * t461) - 0.2e1_dp * t11 * t311 * t736 +  &
      (2._dp * t375 * t468) + t11 * t15 * (F2 * t1061 * t68 - t464 *  &
      t404 - t733 * t365 + 2._dp * t85 * t1073 - t85 * t1090)
      t2056 = -(3._dp * t558 * t761) - (2._dp * t91 * t551 * t1093) - t809  &
      * t519 + 0.5e1_dp / 0.2e1_dp * t659 * t234 * t1946 + t1978 + 0.5e1_dp / 0.2e1_dp  &
      * t1379 * t816 + (t166 * t2018 + t270 * t1816) * alpha7 * t279  &
      - (2._dp * t886 * t692) + (3._dp * t1732 * t813) - t1978 * t156  &
      - t852 * t519 + f12 * t2053 * t94
      t2060 = t759 * t513
      t2066 = t1262 * t769
      t2069 = t125 * t126
      t2071 = t452 * t407
      t2072 = t480 * t113 * t2071
      t2076 = t447 * t115 * t2071
      t2084 = t110 * t111
      t2086 = t1273 * t1270 * t2071
      t2090 = t480 * t114 * t2071
      t2099 = t447 * t127 * t2071
      t2105 = t125 * t739
      t2111 = t114 * t1093
      t2115 = -0.5e1_dp / 0.2e1_dp * t1417 * t2066 + (3._dp * t2069 * t2072) -  &
      0.5e1_dp / 0.2e1_dp * t1417 * t2076 - (2._dp * t603 * t475 * t1841) + (t603  &
      * t604 * t1093) - 0.75e2_dp / 0.4e1_dp * t2084 * t2086 - 0.3e1_dp /  &
      0.2e1_dp * t1473 * t2090 - (2._dp * t583 * t441 * t1841) + (t583  &
      * t584 * t1093) + 0.10e2_dp * t2084 * t2099 + 0.5e1_dp / 0.2e1_dp * (t583)  &
      * (t587) * (t1946) + 0.3e1_dp / 0.2e1_dp * t2105 * t608 + 0.3e1_dp  &
      / 0.2e1_dp * t1453 * t784 + t2105 * t605 + 0.3e1_dp / 0.2e1_dp * (t603)  &
      * (t607) * (t2111)
      t2133 = t110 * t758
      t2136 = t1332 * t783
      t2155 = t1341 * t115 * t2071
      t2158 = -t1434 * t595 * t452 * t407 / 0.2e1_dp - t110 * t1976 * t83 *  &
      t118 + t1407 * t767 - 0.2e1_dp * t81 * t108 * t122 * t442 * t407 - t125  &
      * t2053 * t83 * t129 + t1453 * t781 + 0.5e1_dp / 0.2e1_dp * t2133 * t588  &
      - 0.3e1_dp / 0.2e1_dp * t1473 * t2136 - t1434 * t1435 * t407 / 0.2e1_dp  &
      + t81 * t591 * t1093 + t81 * t596 * t1093 / 0.2e1_dp - 0.3e1_dp / 0.4e1_dp  &
      * t84 * t1442 * t452 * t407 + 0.5e1_dp / 0.2e1_dp * t1407 * t770 + t2133  &
      * t585 - 0.27e2_dp / 0.4e1_dp * t2069 * t2155
      t2180 = t136 * t758
      t2195 = t136 * t137
      t2203 = -t136 * t1976 * t93 * t118 + (2._dp * t1251 * t793) + 0.5e1_dp  &
      / 0.2e1_dp * (t1251) * (t796) + (2._dp * t2180 * t443) - (6._dp  &
      * t440 * t1256 * t1841) - (5._dp * t1261 * t2066) + (2._dp * t440  &
      * t441 * t1093) + 0.5e1_dp / 0.2e1_dp * (t2180) * (t454) - (5._dp  &
      * t1261 * t2076) - 0.75e2_dp / 0.4e1_dp * t2195 * t2086 + 0.10e2_dp *  &
      t2195 * t2099 + 0.5e1_dp / 0.2e1_dp * (t440) * (t448) * (t1946)
      t2211 = t140 * t739
      t2233 = -t140 * t2053 * t93 * t129 + (2._dp * t1321 * t802) + 0.3e1_dp  &
      / 0.2e1_dp * (t1321) * (t805) + (2._dp * t2211 * t476) - (6._dp  &
      * t474 * t1326 * t1841) - (3._dp * t1331 * t2136) + (2._dp * t474  &
      * t475 * t1093) + 0.3e1_dp / 0.2e1_dp * (t2211) * (t483) - (3._dp  &
      * t1331 * t2090) - 0.27e2_dp / 0.4e1_dp * t95 * t2155 + 0.3e1_dp * t95  &
      * t2072 + 0.3e1_dp / 0.2e1_dp * (t474) * (t481) * (t2111)
      t2274 = -0.4e1_dp / 0.3e1_dp * t860 * t153 * t491 + (3._dp * t2060 * t516)  &
      + 0.5e1_dp / 0.3e1_dp * t759 * t441 * t507 + (t2115 + t2158) * omega *  &
      t134 + (12._dp * t106 * t1149 * t442 * t407) - 0.75e2_dp / 0.4e1_dp * (t1397)  &
      * (t234) * (t1270) * (t452) * (t407) -  &
      t788 * t413 / 0.3e1_dp + (t2203 + t2233) * t145 * t151 + (t166 * (0.27e2_dp  &
      / 0.4e1_dp * (t184) * (t1666) * (t219) * (t372) *  &
      (t407) - (3._dp * t184 * t709 * t1546 * t407) - 0.3e1_dp / 0.2e1_dp  &
      * (t184) * (t710) * (t1093) - 0.3e1_dp / 0.4e1_dp * t1678 *  &
      t1756 + (t714 * t1093) / 0.2e1_dp) + f2716 * t1816 * t196) * alpha3  &
      * t145 * t151 - 0.25e2_dp / 0.6e1_dp * (t1362) * (t634) * (t5)  &
      * (t115) * (t407) + (2._dp * t81 * t108 * t442 * t407)  &
      - (5._dp * t1726 * t506 * t6 * t407)
      e_ndrho_rho = e_ndrho_rho + ( -0.4e1_dp / 0.3e1_dp * t409 * t920 - t80 * (t1850 + t1938 + t2056 +  &
      t2274) * Clda ) * sx
      t2279 = t407 ** 2
      t2280 = t1157 * t2279
      t2283 = t541 * t1139
      t2306 = d2exerndrhondrho(Q,dQndrho,d2Qndrhondrho)
      t2316 = t127 * t2279
      t2323 = t1270 * t2279
      t2336 = d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
      t2365 = 2._dp * t1095 * t88 + 4._dp * t375 * t737 + t11 * t15 * (F2 * t1111  &
      * t68 - 2._dp * t733 * t404 + 2._dp * t85 * t1119 - t85 * t1136)
      t2391 = (t166 * (0.3e1_dp / 0.4e1_dp * t206 * t2280 - t206 * t2283 / 0.2e1_dp  &
      - 0.27e2_dp / 0.4e1_dp * t1168 * t256 * t2279 + 0.3e1_dp * t665 * t75 * t2279  &
      + 0.3e1_dp / 0.2e1_dp * t665 * t213 * t1139 + 0.75e2_dp / 0.4e1_dp * t218  &
      * t1182 * t2279 - 0.10e2_dp * t218 * t1186 * t2279 - 0.5e1_dp / 0.2e1_dp *  &
      t218 * t671 * t1139) - t230 * t2306) * alpha5 * t234 + (6._dp * t2060  &
      * t813) - 0.2e1_dp * t91 * t551 * t1139 + 0.10e2_dp * t659 * t234 *  &
      t2316 + 0.6e1_dp * t91 * t1661 * t2279 - 0.75e2_dp / 0.4e1_dp * t1397 * t234  &
      * t2323 + (-0.2e1_dp * t641 * t2279 + t239 * t1139 + 0.6e1_dp * t158  &
      * t282 * t2279 - 0.2e1_dp * t158 * t236 * t1139 + t206 * t2336) * alpha4  &
      * t210 + 0.3e1_dp * t514 * t155 * t44 * t1139 + (5._dp * t1824 *  &
      t816) - (4._dp * t740 * t742) - t81 * t93 * t1139 + f12 * t2365 *  &
      t94 - (6._dp * t759 * t761) + 0.2e1_dp * t81 * t108 * t2279 + (0.12e2_dp  &
      * t1584 * t2279 - 0.3e1_dp * t638 * t1139 - 0.6e1_dp * t1589 * t2279 +  &
      0.2e1_dp * t641 * t1139 + 0.2e1_dp * t206 * t1594 * t2279 - t206 * t644  &
      * t1139 - t243 * t2336) * alpha6 * t247 * t249
      t2432 = 0.147e3_dp / 0.4e1_dp * t251 * t1615 * t2279 - 0.21e2_dp * t251 * t1619  &
      * t2279 - 0.7e1_dp / 0.2e1_dp * t251 * t524 * t1139 - 0.75e2_dp / 0.4e1_dp  &
      * t1630 * t1626 * t2279 + 0.10e2_dp * t530 * t214 * t2279 + 0.5e1_dp /  &
      0.2e1_dp * t530 * t256 * t1139 + 0.27e2_dp / 0.4e1_dp * t261 * t1640 * t2279  &
      - 0.3e1_dp * t261 * t1644 * t2279 - 0.3e1_dp / 0.2e1_dp * t261 * t536 * t1139  &
      - 0.3e1_dp / 0.4e1_dp * t243 * t2280 + t243 * t2283 / 0.2e1_dp
      t2452 = 2._dp * t1095 * t103 + 4._dp * t375 * t756 + t11 * t15 * (2._dp * g2 *  &
      t3 * t19 + 12._dp * g3 * t1 * t24 * t32)
      t2454 = t2452 * E * t108
      t2473 = t2279 * t115
      t2486 = t115 * t1139
      t2501 = t2279 * t114
      t2511 = t113 * t2279
      t2515 = t114 * t1139
      t2519 = -t136 * t2452 * t93 * t118 + (4._dp * t2180 * t793) + (5._dp  &
      * t2180 * t796) - (6._dp * t440 * t1256 * t2279) - (10._dp * t440  &
      * t658 * t2473) + (2._dp * t440 * t441 * t1139) - 0.75e2_dp / 0.4e1_dp  &
      * (t440) * (t1274) * (t2323) + (10._dp * t440 * t448  &
      * t2316) + 0.5e1_dp / 0.2e1_dp * (t440) * (t448) * (t2486) -  &
      t140 * t2365 * t93 * t129 + (4._dp * t2211 * t802) + (3._dp * t2211  &
      * t805) - (6._dp * t474 * t1326 * t2279) - (6._dp * t474 * t108  &
      * t480 * t2501) + (2._dp * t474 * t475 * t1139) - 0.27e2_dp / 0.4e1_dp  &
      * (t474) * (t1342) * (t2473) + (3._dp * t474 * t481 *  &
      t2511) + 0.3e1_dp / 0.2e1_dp * (t474) * (t481) * (t2515)
      t2571 = -t110 * t2452 * t83 * t118 + (2._dp * t2133 * t767) + (5._dp  &
      * t2133 * t770) - (2._dp * t583 * t441 * t2279) + (t583 * t584  &
      * t1139) - (5._dp * t583 * t448 * t2473) + (10._dp * t583 * t587  &
      * t2316) + 0.5e1_dp / 0.2e1_dp * (t583) * (t587) * (t2486)  &
      - 0.75e2_dp / 0.4e1_dp * (t583) * (t1426) * (t2323) - (2._dp  &
      * t81 * t1430 * t2279) - (t81 * t93 * t595 * t2279)
      t2604 = t81 * t591 * t1139 - 0.3e1_dp / 0.4e1_dp * t81 * t1443 * t2279 +  &
      t81 * t596 * t1139 / 0.2e1_dp - t125 * t2365 * t83 * t129 + (2._dp *  &
      t2105 * t781) + (3._dp * t2105 * t784) - 0.2e1_dp * t603 * t475 * t2279  &
      + t603 * t604 * t1139 - 0.3e1_dp * t603 * t481 * t2501 + 0.3e1_dp * t603  &
      * t607 * t2511 + 0.3e1_dp / 0.2e1_dp * t603 * t607 * t2515 - 0.27e2_dp  &
      / 0.4e1_dp * t603 * t1469 * t2473
      t2668 = ((2._dp * t158 * t236 * t2279 - t158 * t204 * t1139 - f98 *  &
      t2336) * alpha2 * t112) + (t166 * t2432 + t270 * t2306) * alpha7  &
      * t279 - t2454 * t156 - (12._dp * t1150 * t155 * t44 * t2279) + (12._dp  &
      * t106 * t1149 * t2279) + t2519 * t145 * t151 + 0.5e1_dp / 0.2e1_dp  &
      * t659 * t234 * t2486 - 0.15e2_dp * t1736 * t234 * t2473 - (3._dp *  &
      t106 * t513 * t1139) + (0.3e1_dp / 0.4e1_dp * (t158) * (t1553)  &
      * (t2279) - (t158 * t496 * t1139) / 0.2e1_dp - t172 * t2306)  &
      * alpha1 * t177 + t2454 + (t2571 + t2604) * omega * t134 + (t166 *  &
      (0.27e2_dp / 0.4e1_dp * (t184) * (t1667) * (t2279) - (3._dp  &
      * t184 * t1671 * t2279) - 0.3e1_dp / 0.2e1_dp * (t184) * (t710)  &
      * (t1139) - 0.3e1_dp / 0.4e1_dp * (t1678) * (t2279) + (t714  &
      * t1139) / 0.2e1_dp) + f2716 * t2306 * t196) * alpha3 * t145 * t151  &
      + (t158 * (t2336 + (t1139 * t161 - 2._dp * t2279 * t204 + 2._dp *  &
      t1538 * t2279 - t626 * t1139) * t202 * t160 - t822 * t238 * t1885  &
      + t823 * t407)) + ((20._dp * t281 * t1691 * t2279 - 4._dp * t281 * t561  &
      * t1139 - 12._dp * t1699 * t2279 + 3._dp * t566 * t1139 + 6._dp * t206 * t565  &
      * t2279 - 2._dp * t206 * t284 * t1139 - 2._dp * t243 * t284 * t2279 + t243  &
      * t238 * t1139 + t289 * t2336) * alpha8 * t296)
      e_ndrho_ndrho = e_ndrho_ndrho + ( -t80 * (t2391 + t2668) * Clda ) * sx
    END IF

  END SUBROUTINE xwpbe_lda_calc_2

! *****************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega screening parameter
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega!=0 using a simple
!>      gaussian expansion for large ww.
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc_3(e_0, e_rho, e_ndrho, e_rho_rho, e_ndrho_rho, &
                              e_ndrho_ndrho, rho, ndrho, omega, sscale, sx, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho, &
                                                e_rho_rho, e_ndrho_rho, &
                                                e_ndrho_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, omega, sscale, sx
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: d2Qndrhondrho, d2Qrhondrho, d2Qrhorho, dQndrho, dQrho, &
      Q, t1, t10, t1004, t1005, t1015, t102, t1025, t1032, t104, t105, t1056, &
      t1060, t1065, t1068, t1079, t108, t1080, t109, t1096, t11, t110, t1101, &
      t1102, t111, t1114, t1115, t1121, t1124, t113, t1143, t1147, t115, &
      t1154, t1156, t116, t1162, t1169, t117, t1170, t1178, t1179, t118, &
      t1189, t119, t1193, t12, t120, t1202, t1203, t1204, t121, t1210, t1213, &
      t1214, t1215, t1216, t1220, t123, t1230, t1235, t124, t1240, t1241, &
      t1242, t125, t1250, t1251, t1252, t1256, t126, t1260, t1264, t127, &
      t1273, t128, t1288, t129, t1295, t13
    REAL(KIND=dp) :: t1300, t1304, t1308, t1309, t131, t1315, t1316, t132, &
      t1326, t133, t1331, t1337, t1346, t1350, t136, t1363, t1372, t1382, &
      t1386, t1388, t1393, t14, t140, t1400, t1401, t1402, t1408, t141, t142, &
      t1437, t144, t1446, t145, t147, t148, t1480, t1482, t1488, t149, t15, &
      t150, t151, t1511, t152, t1522, t1535, t154, t155, t156, t1562, t157, &
      t1578, t158, t1583, t159, t1592, t1594, t16, t1608, t162, t1620, t1627, &
      t163, t1632, t1645, t1652, t166, t167, t1675, t1678, t168, t1685, &
      t1687, t1688, t169, t1692, t1695, t17, t170, t1712, t1719, t1731, &
      t1737, t1739, t174, t1747, t1753, t176
    REAL(KIND=dp) :: t1762, t1765, t177, t1781, t1796, t18, t180, t1804, &
      t181, t1812, t1834, t185, t186, t1860, t1878, t1886, t189, t19, t190, &
      t192, t193, t1935, t194, t1945, t195, t197, t1979, t198, t1981, t1989, &
      t199, t1999, t2, t200, t2003, t2007, t2013, t202, t2050, t206, t2060, &
      t2064, t2068, t208, t209, t21, t210, t2107, t212, t213, t2140, t215, &
      t2159, t216, t217, t2174, t218, t219, t22, t221, t222, t223, t224, &
      t227, t228, t23, t231, t232, t234, t235, t236, t237, t238, t24, t242, &
      t245, t246, t247, t248, t249, t25, t250, t251, t252, t258, t259, t260, &
      t261, t262, t266, t268, t27, t274, t277
    REAL(KIND=dp) :: t278, t279, t28, t280, t285, t287, t288, t29, t291, &
      t292, t295, t296, t297, t3, t300, t301, t302, t304, t305, t308, t309, &
      t31, t313, t316, t317, t318, t32, t320, t322, t325, t328, t329, t332, &
      t335, t338, t339, t34, t341, t342, t344, t345, t35, t352, t359, t36, &
      t361, t364, t365, t368, t369, t370, t373, t374, t375, t376, t377, t379, &
      t38, t380, t381, t383, t384, t388, t39, t391, t395, t396, t397, t398, &
      t399, t4, t400, t401, t403, t404, t405, t406, t408, t41, t416, t417, &
      t418, t419, t42, t420, t422, t423, t424, t428, t429, t433, t435, t437, &
      t438, t44, t441, t442, t443, t444
    REAL(KIND=dp) :: t445, t451, t452, t453, t456, t457, t46, t461, t462, &
      t463, t466, t470, t471, t478, t479, t48, t480, t483, t484, t485, t486, &
      t49, t490, t493, t499, t5, t500, t501, t504, t505, t511, t512, t513, &
      t516, t517, t521, t523, t526, t528, t531, t532, t533, t534, t537, t538, &
      t539, t54, t542, t544, t545, t546, t548, t549, t55, t550, t554, t555, &
      t56, t561, t564, t565, t567, t568, t570, t58, t584, t586, t589, t590, &
      t592, t593, t595, t596, t599, t6, t60, t603, t604, t607, t608, t61, &
      t610, t611, t612, t616, t617, t621, t623, t626, t627, t628, t629, t63, &
      t635, t637, t638, t649, t65, t651
    REAL(KIND=dp) :: t652, t656, t661, t664, t67, t670, t673, t677, t68, &
      t681, t684, t687, t69, t690, t691, t695, t696, t698, t699, t7, t70, &
      t704, t705, t706, t708, t709, t71, t712, t713, t714, t717, t718, t72, &
      t721, t724, t725, t727, t73, t74, t743, t746, t748, t75, t752, t756, &
      t759, t760, t762, t765, t767, t768, t769, t77, t772, t78, t781, t8, &
      t80, t803, t804, t806, t81, t811, t813, t816, t82, t820, t83, t84, &
      t843, t844, t849, t85, t855, t87, t871, t872, t875, t877, t878, t88, &
      t892, t893, t899, t9, t90, t900, t91, t916, t917, t92, t920, t922, t93, &
      t932, t933, t94, t95, t96, t964, t967, t968
    REAL(KIND=dp) :: t969, t97, t973, t976, t99

    IF ( order >= 0 ) THEN
      t1 = ndrho ** 2
      t2 = r2 ** 2
      t3 = 0.1e1_dp / t2
      t4 = t1 * t3
      t5 = pi ** 2
      t6 = r3 * t5
      t7 = t6 * rho
      t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
      t9 = t8 ** 2
      t10 = 0.1e1_dp / t9
      t11 = t4 * t10
      t12 = rho ** 2
      t13 = 0.1e1_dp / t12
      t14 = sscale ** 2
      t15 = t13 * t14
      t16 = a1 * t1
      t17 = t16 * t3
      t18 = t10 * t13
      t19 = t18 * t14
      t21 = t1 ** 2
      t22 = a2 * t21
      t23 = t2 ** 2
      t24 = 0.1e1_dp / t23
      t25 = t22 * t24
      t27 = 0.1e1_dp / t8 / t7
      t28 = t12 ** 2
      t29 = 0.1e1_dp / t28
      t31 = t14 ** 2
      t32 = t27 * t29 * t31
      t34 = t17 * t19 + t25 * t32
      t35 = a3 * t21
      t36 = t35 * t24
      t38 = t21 * ndrho
      t39 = a4 * t38
      t41 = 0.1e1_dp / t23 / r2
      t42 = t39 * t41
      t44 = 0.1e1_dp / t9 / t7
      t46 = 0.1e1_dp / t28 / rho
      t48 = t31 * sscale
      t49 = t44 * t46 * t48
      t54 = 0.1e1_dp / t23 / t2
      t55 = a5 * t21 * t1 * t54
      t56 = r3 ** 2
      t58 = t5 ** 2
      t60 = 0.1e1_dp / t56 / t58
      t61 = t28 ** 2
      t63 = t31 * t14
      t65 = t60 / t61 * t63
      t67 = r1 + t36 * t32 + t42 * t49 + t55 * t65
      t68 = 0.1e1_dp / t67
      t69 = t34 * t68
      t70 = t15 * t69
      t71 = t11 * t70
      t72 = omega ** 2
      t73 = beta2 * t72
      t74 = t73 * t10
      t75 = t71 + t74
      t77 = 0.1e1_dp / A
      Q = f94 * t75 * t77
      t78 = rho ** (0.1e1_dp / 0.3e1_dp)
      t80 = t78 * rho * f89
      t81 = B * f12
      t82 = t71 + DD
      t83 = 0.1e1_dp / t82
      t84 = t81 * t83
      t85 = F2 * t34
      t87 = F1 + t85 * t68
      t88 = t15 * t87
      t90 = t11 * t88 + r1
      t91 = f12 * t90
      t92 = t82 ** 2
      t93 = 0.1e1_dp / t92
      t94 = C * t93
      t95 = t91 * t94
      t96 = f34 * pi
      t97 = rootpi
      t99 = r6 * C
      t102 = r4 * B
      t104 = r8 * A
      t105 = t92 * t82
      t108 = t97 * (r15 * E + t99 * t90 * t82 + t102 * t92 + t104 * t105)
      t109 = 0.1e1_dp / r16
      t110 = SQRT(t82)
      t111 = t110 * t105
      t113 = t109 / t111
      t115 = SQRT(A)
      t116 = f94 * t34
      t117 = t68 * t1
      t118 = t116 * t117
      t119 = t3 * t10
      t120 = t15 * t77
      t121 = t119 * t120
      t123 = EXP(t118 * t121)
      t124 = t115 * t123
      t125 = f32 * ndrho
      t126 = 0.1e1_dp / r2
      t127 = t125 * t126
      t128 = 0.1e1_dp / t8
      t129 = 0.1e1_dp / rho
      t131 = t69 * t77
      t132 = SQRT(t131)
      t133 = sscale * t132
      t136 = erfc(t127 * t128 * t129 * t133)
      t140 = 0.1e1_dp / f1516
      t141 = (t96 + t108 * t113 - t96 * t124 * t136) * t140
      t142 = 0.1e1_dp / t97
      t144 = 0.1e1_dp / E
      t145 = t142 * t111 * t144
      t147 = -t141 * t145 + r1
      t148 = t147 * E
      t149 = 0.1e1_dp / t105
      t150 = t148 * t149
      t151 = f158 * E
      t152 = t147 * t83
      t154 = t71 + DD + t72 * t10
      t155 = t154 ** 2
      t156 = t155 ** 2
      t157 = t156 * t154
      t158 = SQRT(t157)
      t159 = 0.1e1_dp / t158
      t162 = SQRT(t154)
      t163 = 0.1e1_dp / t162
      t166 = f68 * C
      t167 = t90 * t83
      t168 = t155 * t154
      t169 = SQRT(t168)
      t170 = 0.1e1_dp / t169
      t174 = (-t151 * t152 * t159 - t81 * t83 * t163 - t166 * t167 * t170)  &
      * omega
      t176 = f52 * E
      t177 = t147 * t93
      t180 = f12 * C
      t181 = t90 * t93
      t185 = t72 * omega
      t186 = (-t176 * t177 * t159 - t180 * t181 * t170) * t185
      t189 = 0.1e1_dp / r3 / t5
      t190 = t189 * t129
      t192 = t72 ** 2
      t193 = t192 * omega
      t194 = t159 * t193
      t195 = t194 * t44
      t197 = f12 * A
      t198 = exei(Q)
      t199 = t71 + DD + t74
      t200 = 0.1e1_dp / t199
      t202 = LOG(t75 * t200)
      t206 = (t84 + t95 + t150 + t174 * t128 + t186 * t190 - t150 * t195  &
      + t197 * (t198 + t202)) * Clda
      e_0 = e_0 + ( -t80 * t206 ) * sx
    END IF
    IF( order >= 1 .OR. order == -1 ) THEN
      t208 = t44 * t13
      t209 = t4 * t208
      t210 = t14 * t34
      t212 = t68 * r3 * t5
      t213 = t210 * t212
      t215 = 0.2e1_dp / 0.3e1_dp * t209 * t213
      t216 = t12 * rho
      t217 = 0.1e1_dp / t216
      t218 = t217 * t14
      t219 = t218 * t69
      t221 = 2._dp * t11 * t219
      t222 = t3 * t44
      t223 = t16 * t222
      t224 = t15 * t6
      t227 = t10 * t217
      t228 = t227 * t14
      t231 = t56 * t58
      t232 = t231 * t12
      t234 = 0.1e1_dp / t8 / t232
      t235 = t24 * t234
      t236 = t22 * t235
      t237 = t29 * t31
      t238 = t237 * t6
      t242 = t27 * t46 * t31
      t245 = -0.2e1_dp / 0.3e1_dp * t223 * t224 - (2._dp * t17 * t228) - 0.4e1_dp  &
      / 0.3e1_dp * t236 * t238 - (4._dp * t25 * t242)
      t246 = t245 * t68
      t247 = t15 * t246
      t248 = t11 * t247
      t249 = t4 * t18
      t250 = t67 ** 2
      t251 = 0.1e1_dp / t250
      t252 = t35 * t235
      t258 = 0.1e1_dp / t9 / t232
      t259 = t41 * t258
      t260 = t39 * t259
      t261 = t46 * t48
      t262 = t261 * t6
      t266 = 0.1e1_dp / t28 / t12
      t268 = t44 * t266 * t48
      t274 = t60 / t61 / rho * t63
      t277 = -0.4e1_dp / 0.3e1_dp * t252 * t238 - (4._dp * t36 * t242) - 0.5e1_dp  &
      / 0.3e1_dp * t260 * t262 - (5._dp * t42 * t268) - (8._dp * t55 * t274)
      t278 = t251 * t277
      t279 = t210 * t278
      t280 = t249 * t279
      t285 = -t215 - t221 + t248 - t280 - 0.2e1_dp / 0.3e1_dp * t73 * t44 * r3  &
      * t5
      dQrho = f94 * t285 * t77
      t287 = ndrho * t3
      t288 = t287 * t10
      t291 = a1 * ndrho
      t292 = t291 * t3
      t295 = t1 * ndrho
      t296 = a2 * t295
      t297 = t296 * t24
      t300 = 2._dp * t292 * t19 + 4._dp * t297 * t32
      t301 = t300 * t68
      t302 = t15 * t301
      t304 = a3 * t295
      t305 = t304 * t24
      t308 = a4 * t21
      t309 = t308 * t41
      t313 = a5 * t38 * t54
      t316 = 4._dp * t305 * t32 + 5._dp * t309 * t49 + 6._dp * t313 * t65
      t317 = t251 * t316
      t318 = t210 * t317
      t320 = 2._dp * t288 * t70 + t11 * t302 - t249 * t318
      dQndrho = f94 * t320 * t77
      t322 = t78 * f89
      t325 = -t215 - t221 + t248 - t280
      t328 = t14 * t87
      t329 = t328 * t6
      t332 = t218 * t87
      t335 = F2 * t245
      t338 = t335 * t68 - t85 * t278
      t339 = t15 * t338
      t341 = -0.2e1_dp / 0.3e1_dp * t209 * t329 - (2._dp * t11 * t332) + (t11  &
      * t339)
      t342 = f12 * t341
      t344 = C * t149
      t345 = t344 * t325
      t352 = t82 * t325
      t359 = t97 * (t99 * t341 * t82 + t99 * t90 * t325 + 2._dp * t102 * t352  &
      + 3._dp * t104 * t92 * t325)
      t361 = t92 ** 2
      t364 = t109 / t110 / t361
      t365 = t364 * t325
      t368 = t96 * t115
      t369 = f94 * t245
      t370 = t369 * t117
      t373 = t251 * t1 * t3
      t374 = t116 * t373
      t375 = t14 * t77
      t376 = t375 * t277
      t377 = t18 * t376
      t379 = t117 * t3
      t380 = t116 * t379
      t381 = t208 * t14
      t383 = t77 * r3 * t5
      t384 = t381 * t383
      t388 = t119 * t218 * t77
      t391 = t370 * t121 - t374 * t377 - 0.2e1_dp / 0.3e1_dp * t380 * t384 - (2._dp  &
      * t118 * t388)
      t395 = rootpi
      t396 = 0.1e1_dp / t395
      t397 = t123 * t396
      t398 = f32 ** 2
      t399 = t398 * t1
      t400 = t399 * t119
      t401 = t15 * t131
      t403 = EXP(-t400 * t401)
      t404 = t126 * t27
      t405 = t125 * t404
      t406 = t129 * sscale
      t408 = t132 * r3 * t5
      t416 = t125 * t126 * t128
      t417 = 0.1e1_dp / t132
      t418 = t246 * t77
      t419 = t34 * t251
      t420 = t77 * t277
      t422 = t418 - t419 * t420
      t423 = t417 * t422
      t424 = t406 * t423
      t428 = t403 * (-t405 * t406 * t408 / 0.3e1_dp - t127 * t128 * t13 * t133  &
      + t416 * t424 / 0.2e1_dp)
      t429 = t397 * t428
      t433 = (t359 * t113 - 0.7e1_dp / 0.2e1_dp * t108 * t365 - (t368 * t391  &
      * t123 * t136) + (2._dp * t368 * t429)) * t140
      t435 = t141 * t142
      t437 = t110 * t92 * t144
      t438 = t437 * t325
      t441 = -t433 * t145 - 0.7e1_dp / 0.2e1_dp * t435 * t438
      t442 = t441 * E
      t443 = t442 * t149
      t444 = 0.1e1_dp / t361
      t445 = t444 * t325
      t451 = t151 * t147
      t452 = t93 * t159
      t453 = t452 * t325
      t456 = 0.1e1_dp / t158 / t157
      t457 = t83 * t456
      t461 = -t215 - t221 + t248 - t280 - 0.2e1_dp / 0.3e1_dp * t72 * t44 * t6
      t462 = t156 * t461
      t463 = t457 * t462
      t466 = t93 * t163
      t470 = 0.1e1_dp / t162 / t154
      t471 = t83 * t470
      t478 = t166 * t90
      t479 = t93 * t170
      t480 = t479 * t325
      t483 = 0.1e1_dp / t169 / t168
      t484 = t83 * t483
      t485 = t155 * t461
      t486 = t484 * t485
      t490 = (-t151 * t441 * t83 * t159 + t451 * t453 + 0.5e1_dp / 0.2e1_dp * t451  &
      * t463 + t81 * t466 * t325 + t81 * t471 * t461 / 0.2e1_dp - t166  &
      * t341 * t83 * t170 + t478 * t480 + 0.3e1_dp / 0.2e1_dp * t478 * t486) *  &
      omega
      t493 = t27 * r3 * t5
      t499 = t176 * t147
      t500 = t149 * t159
      t501 = t500 * t325
      t504 = t93 * t456
      t505 = t504 * t462
      t511 = t180 * t90
      t512 = t149 * t170
      t513 = t512 * t325
      t516 = t93 * t483
      t517 = t516 * t485
      t521 = (-t176 * t441 * t93 * t159 + (2._dp * t499 * t501) + 0.5e1_dp /  &
      0.2e1_dp * (t499) * (t505) - t180 * t341 * t93 * t170 + (2._dp  &
      * t511 * t513) + 0.3e1_dp / 0.2e1_dp * (t511) * (t517)) * t185
      t523 = t189 * t13
      t526 = t148 * t444
      t528 = t194 * t44 * t325
      t531 = t149 * t456
      t532 = t148 * t531
      t533 = t193 * t44
      t534 = t533 * t462
      t537 = t148 * t500
      t538 = t193 * t258
      t539 = t538 * t6
      t542 = dexeirho(Q,dQrho)
      t544 = t199 ** 2
      t545 = 0.1e1_dp / t544
      t546 = t75 * t545
      t548 = t285 * t200 - t546 * t285
      t549 = 0.1e1_dp / t75
      t550 = t548 * t549
      t554 = -t81 * t93 * t325 + t342 * t94 - (2._dp * t91 * t345) + t443  &
      - (3._dp * t148 * t445) + t490 * t128 - t174 * t493 / 0.3e1_dp + t521  &
      * t190 - t186 * t523 - t443 * t195 + (3._dp * t526 * t528) + 0.5e1_dp  &
      / 0.2e1_dp * t532 * t534 + 0.5e1_dp / 0.3e1_dp * t537 * t539 + t197 * (t542  &
      + t550 * t199)
      t555 = t554 * Clda
      e_rho = e_rho + ( -0.4e1_dp / 0.3e1_dp * t322 * t206 - t80 * t555 ) * sx
      t561 = F2 * t300
      t564 = t561 * t68 - t85 * t317
      t565 = t15 * t564
      t567 = 2._dp * t288 * t88 + t11 * t565
      t568 = f12 * t567
      t570 = t344 * t320
      t584 = t97 * (t99 * t567 * t82 + t99 * t90 * t320 + 2._dp * t102 * t82  &
      * t320 + 3._dp * t104 * t92 * t320)
      t586 = t364 * t320
      t589 = f94 * t300
      t590 = t589 * t117
      t592 = t375 * t316
      t593 = t18 * t592
      t595 = t68 * ndrho
      t596 = t116 * t595
      t599 = t590 * t121 - t374 * t593 + 2._dp * t596 * t121
      t603 = f32 * t126
      t604 = t603 * t128
      t607 = t301 * t77
      t608 = t77 * t316
      t610 = t607 - t419 * t608
      t611 = t417 * t610
      t612 = t406 * t611
      t616 = t403 * (t604 * t406 * t132 + t416 * t612 / 0.2e1_dp)
      t617 = t397 * t616
      t621 = (t584 * t113 - 0.7e1_dp / 0.2e1_dp * t108 * t586 - (t368 * t599  &
      * t123 * t136) + (2._dp * t368 * t617)) * t140
      t623 = t437 * t320
      t626 = -t621 * t145 - 0.7e1_dp / 0.2e1_dp * t435 * t623
      t627 = t626 * E
      t628 = t627 * t149
      t629 = t444 * t320
      t635 = t452 * t320
      t637 = t156 * t320
      t638 = t457 * t637
      t649 = t479 * t320
      t651 = t155 * t320
      t652 = t484 * t651
      t656 = (-t151 * t626 * t83 * t159 + t451 * t635 + 0.5e1_dp / 0.2e1_dp * t451  &
      * t638 + t81 * t466 * t320 + t81 * t471 * t320 / 0.2e1_dp - t166  &
      * t567 * t83 * t170 + t478 * t649 + 0.3e1_dp / 0.2e1_dp * t478 * t652) *  &
      omega
      t661 = t500 * t320
      t664 = t504 * t637
      t670 = t512 * t320
      t673 = t516 * t651
      t677 = (-t176 * t626 * t93 * t159 + (2._dp * t499 * t661) + 0.5e1_dp /  &
      0.2e1_dp * (t499) * (t664) - t180 * t567 * t93 * t170 + (2._dp  &
      * t511 * t670) + 0.3e1_dp / 0.2e1_dp * (t511) * (t673)) * t185
      t681 = t194 * t44 * t320
      t684 = t533 * t637
      t687 = dexeindrho(Q,dQndrho)
      t690 = t320 * t200 - t546 * t320
      t691 = t690 * t549
      t695 = -t81 * t93 * t320 + t568 * t94 - (2._dp * t91 * t570) + t628  &
      - (3._dp * t148 * t629) + t656 * t128 + t677 * t190 - t628 * t195  &
      + (3._dp * t526 * t681) + 0.5e1_dp / 0.2e1_dp * t532 * t684 + t197 * (t687  &
      + t691 * t199)
      t696 = t695 * Clda
      e_ndrho = e_ndrho + ( -t80 * t696 ) * sx
    END IF
    IF( order >= 2 .OR. order == -2 ) THEN
      t698 = t258 * t13
      t699 = t4 * t698
      t704 = 0.10e2_dp / 0.9e1_dp * t699 * t210 * t68 * t56 * t58
      t705 = t44 * t217
      t706 = t4 * t705
      t708 = 0.8e1_dp / 0.3e1_dp * t706 * t213
      t709 = t14 * t245
      t712 = 0.4e1_dp / 0.3e1_dp * t209 * t709 * t212
      t713 = t4 * t381
      t714 = t6 * t277
      t717 = 0.4e1_dp / 0.3e1_dp * t713 * t419 * t714
      t718 = t29 * t14
      t721 = 6._dp * t11 * t718 * t69
      t724 = 4._dp * t11 * t218 * t246
      t725 = t4 * t227
      t727 = 4._dp * t725 * t279
      t743 = t56 * r3 * t58 * t5 * t216
      t746 = t24 / t8 / t743
      t748 = t237 * t231
      t752 = t46 * t31 * t6
      t756 = t27 * t266 * t31
      t759 = 0.10e2_dp / 0.9e1_dp * t16 * t3 * t258 * t15 * t231 + 0.8e1_dp / 0.3e1_dp  &
      * t223 * t218 * t6 + (6._dp * t17 * t10 * t29 * t14) + 0.28e2_dp /  &
      0.9e1_dp * t22 * t746 * t748 + 0.32e2_dp / 0.3e1_dp * t236 * t752 + (20._dp  &
      * t25 * t756)
      t760 = t759 * t68
      t762 = t11 * t15 * t760
      t765 = 2._dp * t249 * t709 * t278
      t767 = 0.1e1_dp / t250 / t67
      t768 = t277 ** 2
      t769 = t767 * t768
      t772 = 2._dp * t249 * t210 * t769
      t781 = 0.1e1_dp / t9 / t743
      t803 = 0.28e2_dp / 0.9e1_dp * t35 * t746 * t748 + 0.32e2_dp / 0.3e1_dp * t252 *  &
      t752 + (20._dp * t36 * t756) + 0.40e2_dp / 0.9e1_dp * t39 * t41 * t781  &
      * t261 * t231 + 0.50e2_dp / 0.3e1_dp * t260 * t266 * t48 * t6 + 0.30e2_dp *  &
      t42 * t44 / t28 / t216 * t48 + (72._dp * t55 * t60 / t61 / t12 *  &
      t63)
      t804 = t251 * t803
      t806 = t249 * t210 * t804
      t811 = t704 + t708 - t712 + t717 + t721 - t724 + t727 + t762 - t765  &
      + t772 - t806 + 0.10e2_dp / 0.9e1_dp * t73 * t258 * t56 * t58
      d2Qrhorho = f94 * t811 * t77
      t813 = t287 * t208
      t816 = t14 * t300
      t820 = t6 * t316
      t843 = -0.4e1_dp / 0.3e1_dp * t291 * t222 * t224 - (4._dp * t292 * t228)  &
      - 0.16e2_dp / 0.3e1_dp * t296 * t235 * t238 - (16._dp * t297 * t242)
      t844 = t843 * t68
      t849 = t287 * t18
      t855 = t767 * t277 * t316
      t871 = -0.16e2_dp / 0.3e1_dp * t304 * t235 * t238 - (16._dp * t305 * t242)  &
      - 0.25e2_dp / 0.3e1_dp * t308 * t259 * t262 - (25._dp * t309 * t268) -  &
      (48._dp * t313 * t274)
      t872 = t251 * t871
      t875 = -0.4e1_dp / 0.3e1_dp * t813 * t213 - 0.2e1_dp / 0.3e1_dp * t209 * t816 *  &
      t212 + 0.2e1_dp / 0.3e1_dp * t713 * t419 * t820 - (4._dp * t288 * t219)  &
      - (2._dp * t11 * t218 * t301) + (2._dp * t725 * t318) + (2._dp *  &
      t288 * t247) + (t11 * t15 * t844) - t249 * t709 * t317 - (2._dp  &
      * t849 * t279) - t249 * t816 * t278 + 0.2e1_dp * t249 * t210 * t855  &
      - t249 * t210 * t872
      d2Qrhondrho = f94 * t875 * t77
      t877 = t119 * t13
      t878 = t210 * t68
      t892 = 2._dp * a1 * t3 * t19 + 12._dp * a2 * t1 * t24 * t32
      t893 = t892 * t68
      t899 = t316 ** 2
      t900 = t767 * t899
      t916 = 12._dp * a3 * t1 * t24 * t32 + 20._dp * a4 * t295 * t41 * t49 + 30._dp *  &
      a5 * t21 * t54 * t65
      t917 = t251 * t916
      t920 = 2._dp * t877 * t878 + 4._dp * t288 * t302 - 4._dp * t849 * t318 + t11 *  &
      t15 * t893 - 2._dp * t249 * t816 * t317 + 2._dp * t249 * t210 * t900 - t249  &
      * t210 * t917
      d2Qndrhondrho = f94 * t920 * t77
      t922 = t78 ** 2
      t932 = t704 + t708 - t712 + t717 + t721 - t724 + t727 + t762 - t765  &
      + t772 - t806 + 0.10e2_dp / 0.9e1_dp * t72 * t258 * t231
      t933 = t156 * t932
      t964 = 0.10e2_dp / 0.9e1_dp * t699 * t328 * t231 + 0.8e1_dp / 0.3e1_dp * t706 *  &
      t329 - 0.4e1_dp / 0.3e1_dp * t209 * t14 * t338 * t6 + (6._dp * t11 * t718  &
      * t87) - 0.4e1_dp * (t11) * t218 * t338 + (t11 * t15 * (F2  &
      * t759 * t68 - 2._dp * t335 * t278 + 2._dp * t85 * t769 - t85 * t804))
      t967 = t361 * t82
      t968 = 0.1e1_dp / t967
      t969 = t325 ** 2
      t973 = t442 * t444
      t976 = t704 + t708 - t712 + t717 + t721 - t724 + t727 + t762 - t765  &
      + t772 - t806
      t1004 = 0.1e1_dp / t110 / t967
      t1005 = t109 * t1004
      t1015 = t369 * t373
      t1025 = t116 * t767 * t1 * t3
      t1032 = t116 * t251 * t4 * t44
      t1056 = f94 * t759 * t117 * t121 - (2._dp * t1015 * t377) - 0.4e1_dp /  &
      0.3e1_dp * t369 * t379 * t384 - (4._dp * t370 * t388) + (2._dp * t1025  &
      * t18 * t375 * t768) + 0.4e1_dp / 0.3e1_dp * t1032 * t120 * t714 + (4._dp  &
      * t374 * t227 * t376) - (t374 * t18 * t375 * t803) + 0.10e2_dp  &
      / 0.9e1_dp * t380 * t698 * t14 * t77 * t56 * t58 + 0.8e1_dp / 0.3e1_dp *  &
      t380 * t705 * t14 * t383 + 0.6e1_dp * t118 * t119 * t718 * t77
      t1060 = t391 ** 2
      t1065 = t96 * t115 * t391
      t1068 = t96 * t124
      t1079 = t399 * t877
      t1080 = t251 * t77
      t1096 = t13 * sscale
      t1101 = t125 * t404 * t129
      t1102 = sscale * t417
      t1114 = 0.1e1_dp / t132 / t131
      t1115 = t422 ** 2
      t1121 = t245 * t251
      t1124 = t34 * t767
      t1143 = t433 * t142
      t1147 = t110 * t82 * t144
      t1154 = -((t97 * (t99 * t964 * t82 + 2._dp * t99 * t341 * t325 + t99  &
      * t90 * t976 + 2._dp * t102 * t969 + 2._dp * t102 * t82 * t976 + 6._dp * t104  &
      * t82 * t969 + 3._dp * t104 * t92 * t976) * t113) - (7._dp * t359 *  &
      t365) + 0.63e2_dp / 0.4e1_dp * (t108) * (t1005) * (t969) - 0.7e1_dp  &
      / 0.2e1_dp * (t108) * (t364) * (t976) - t368 * t1056  &
      * t123 * t136 - t368 * t1060 * t123 * t136 + (4._dp * t1065 * t429)  &
      + 0.2e1_dp * t1068 * t396 * (0.2e1_dp / 0.3e1_dp * t399 * t222 * t13 * t878  &
      * t383 + (2._dp * t400 * t218 * t131) - (t400 * t15 * t418)  &
      + t1079 * t210 * t1080 * t277) * t428 + 0.2e1_dp * t368 * t397 * t403  &
      * (0.4e1_dp / 0.9e1_dp * t125 * t126 * t234 * t406 * t132 * t56 * t58  &
      + 0.2e1_dp / 0.3e1_dp * t405 * t1096 * t408 - t1101 * t1102 * t6 * t422  &
      / 0.3e1_dp + (2._dp * t127 * t128 * t217 * t133) - t416 * t1096 * t423  &
      - t416 * t406 * t1114 * t1115 / 0.4e1_dp + t416 * t406 * t417 * (t760  &
      * t77 - 2._dp * t1121 * t420 + 2._dp * t1124 * t77 * t768 - t419 *  &
      t77 * t803) / 0.2e1_dp)) * t140 * t145 - (7._dp * t1143 * t438) - 0.35e2_dp  &
      / 0.4e1_dp * (t435) * (t1147) * (t969) - 0.7e1_dp / 0.2e1_dp  &
      * (t435) * (t437) * (t976)
      t1156 = t1154 * E * t149
      t1162 = t442 * t531
      t1169 = t148 * t444 * t456
      t1170 = t325 * t156
      t1178 = t444 * t159
      t1179 = t148 * t1178
      t1189 = 0.5e1_dp / 0.2e1_dp * t532 * t533 * t933 - (6._dp * t442 * t445)  &
      + f12 * t964 * t94 + (12._dp * t148 * t968 * t969) + (6._dp * t973  &
      * t528) - (3._dp * t148 * t444 * t976) + t1156 - (t81 * t93 *  &
      t976) + (2._dp * t81 * t149 * t969) + (5._dp * t1162 * t534) + 0.10e2_dp  &
      / 0.3e1_dp * (t442) * (t500) * (t539) - 0.15e2_dp * t1169  &
      * t533 * t1170 * t461 - (2._dp * t91 * t344 * t976) - (10._dp *  &
      t1179 * t538 * t325 * r3 * t5) + 0.4e1_dp / 0.9e1_dp * t174 * t234 * t56  &
      * t58
      t1193 = t176 * t441
      t1202 = t176 * t147 * t149
      t1203 = t456 * t325
      t1204 = t1203 * t462
      t1210 = t156 ** 2
      t1213 = 0.1e1_dp / t158 / t1210 / t155
      t1214 = t93 * t1213
      t1215 = t461 ** 2
      t1216 = t1210 * t1215
      t1220 = t168 * t1215
      t1230 = t180 * t341
      t1235 = t444 * t170
      t1240 = t180 * t90 * t149
      t1241 = t483 * t325
      t1242 = t1241 * t485
      t1250 = 0.1e1_dp / t169 / t156 / t155
      t1251 = t93 * t1250
      t1252 = t156 * t1215
      t1256 = t154 * t1215
      t1260 = t155 * t932
      t1264 = -t176 * t1154 * t93 * t159 + (4._dp * t1193 * t501) + (5._dp  &
      * t1193 * t505) - (6._dp * t499 * t1178 * t969) - (10._dp * t1202  &
      * t1204) + (2._dp * t499 * t500 * t976) - 0.75e2_dp / 0.4e1_dp * (t499)  &
      * (t1214) * (t1216) + (10._dp * t499 * t504 * t1220)  &
      + 0.5e1_dp / 0.2e1_dp * (t499) * (t504) * (t933) - t180 * t964  &
      * t93 * t170 + (4._dp * t1230 * t513) + (3._dp * t1230 * t517)  &
      - (6._dp * t511 * t1235 * t969) - (6._dp * t1240 * t1242) + (2._dp  &
      * t511 * t512 * t976) - 0.27e2_dp / 0.4e1_dp * (t511) * (t1251)  &
      * (t1252) + (3._dp * t511 * t516 * t1256) + 0.3e1_dp / 0.2e1_dp *  &
      (t511) * (t516) * (t1260)
      t1273 = t148 * t149 * t1213
      t1288 = t148 * t531 * t193
      t1295 = t148 * t968
      t1300 = t83 * t1213
      t1304 = t149 * t163
      t1308 = t81 * t93
      t1309 = t470 * t325
      t1315 = 0.1e1_dp / t162 / t155
      t1316 = t83 * t1315
      t1326 = t166 * t341
      t1331 = t166 * t181
      t1337 = -0.75e2_dp / 0.4e1_dp * t451 * t1300 * t1216 - (2._dp * t81 * t1304  &
      * t969) - t1308 * t1309 * t461 + (t81 * t466 * t976) - 0.3e1_dp  &
      / 0.4e1_dp * (t81) * (t1316) * (t1215) + (t81 * t471  &
      * t932) / 0.2e1_dp - t166 * t964 * t83 * t170 + (2._dp * t1326 * t480)  &
      + (3._dp * t1326 * t486) - (3._dp * t1331 * t1242) - (2._dp *  &
      t478 * t512 * t969)
      t1346 = t83 * t1250
      t1350 = t151 * t441
      t1363 = t151 * t177
      t1372 = (t478 * t479 * t976) + (3._dp * t478 * t484 * t1256) +  &
      0.3e1_dp / 0.2e1_dp * (t478) * (t484) * (t1260) - 0.27e2_dp / 0.4e1_dp  &
      * (t478) * (t1346) * (t1252) + (2._dp * t1350 * t453)  &
      - t151 * t1154 * t83 * t159 + (5._dp * t1350 * t463) - (2._dp  &
      * t451 * t500 * t969) + (t451 * t452 * t976) - (5._dp * t1363  &
      * t1204) + (10._dp * t451 * t457 * t1220) + 0.5e1_dp / 0.2e1_dp * (t451)  &
      * (t457) * (t933)
      t1382 = C * t444
      t1386 = d2exeirhorho(Q,dQrho,d2Qrhorho)
      t1388 = t285 ** 2
      t1393 = t75 / t544 / t199
      t1400 = t75 ** 2
      t1401 = 0.1e1_dp / t1400
      t1402 = t548 * t1401
      t1408 = t1264 * t185 * t190 + (10._dp * t532 * t533 * t1220) - (2._dp  &
      * t521 * t523) - 0.75e2_dp / 0.4e1_dp * (t1273) * (t533) * (t1216)  &
      + (3._dp * t526 * t194 * t44 * t976) - t1156 * t195 - 0.2e1_dp  &
      / 0.3e1_dp * t490 * t493 + (2._dp * t186 * t189 * t217) - 0.25e2_dp  &
      / 0.3e1_dp * t1288 * t258 * t156 * t461 * r3 * t5 - (12._dp * t1295 *  &
      t194 * t44 * t969) + (t1337 + t1372) * omega * t128 - 0.40e2_dp / 0.9e1_dp  &
      * t537 * t193 * t781 * t231 - (4._dp * t342 * t345) + (6._dp *  &
      t91 * t1382 * t969) + (t197 * (t1386 + (t811 * t200 - 2._dp * t1388  &
      * t545 + 2._dp * t1393 * t1388 - t546 * t811) * t549 * t199 - t1402  &
      * t199 * t285 + t550 * t285))
      e_rho_rho = e_rho_rho + ( -0.4e1_dp / 0.9e1_dp / t922 * f89 * t206 - 0.8e1_dp / 0.3e1_dp * t322 * t555  &
      - t80 * (t1189 + t1408) * Clda ) * sx
      t1437 = -0.4e1_dp / 0.3e1_dp * t813 * t329 - 0.2e1_dp / 0.3e1_dp * t209 * t14 *  &
      t564 * t6 - (4._dp * t288 * t332) - 0.2e1_dp * t11 * t218 * t564 + (2._dp  &
      * t288 * t339) + t11 * t15 * (F2 * t843 * t68 - t335 * t317  &
      - t561 * t278 + 2._dp * t85 * t855 - t85 * t872)
      t1446 = t320 * t325
      t1480 = t589 * t373
      t1482 = t420 * t316
      t1488 = t116 * t251 * ndrho * t3
      t1511 = (f94 * t843 * t117 * t121) - t1015 * t593 + (2._dp * t369  &
      * t595 * t121) - (t1480 * t377) + (2._dp * t1025 * t19 * t1482)  &
      - (2._dp * t1488 * t377) - (t374 * t18 * t375 * t871) - 0.2e1_dp  &
      / 0.3e1_dp * t589 * t379 * t384 + 0.2e1_dp / 0.3e1_dp * t1032 * t120 *  &
      t820 - 0.4e1_dp / 0.3e1_dp * t116 * (t595) * t3 * t384 - (2._dp * t590  &
      * t388) + (2._dp * t374 * t227 * t592) - (4._dp * t596 * t388)
      t1522 = t96 * t115 * t599
      t1535 = t396 * (-2._dp * t398 * ndrho * t119 * t401 - t400 * t15 * t607  &
      + t1079 * t210 * t1080 * t316)
      t1562 = t300 * t251
      t1578 = (t97 * (t99 * t1437 * t82 + t99 * t341 * t320 + t99 * t567  &
      * t325 + t99 * t90 * t875 + 2._dp * t102 * t1446 + 2._dp * t102 * t82  &
      * t875 + 6._dp * t104 * t352 * t320 + 3._dp * t104 * t92 * t875) * t113) -  &
      0.7e1_dp / 0.2e1_dp * t359 * t586 - 0.7e1_dp / 0.2e1_dp * t584 * t365 + 0.63e2_dp  &
      / 0.4e1_dp * (t108) * (t109) * (t1004) * (t325) * (t320)  &
      - 0.7e1_dp / 0.2e1_dp * (t108) * (t364) * (t875) -  &
      t368 * t1511 * t123 * t136 - t368 * t391 * t599 * t123 * t136 + (2._dp  &
      * t1065 * t617) + (2._dp * t1522 * t429) + (2._dp * t1068 *  &
      t1535 * t428) + 0.2e1_dp * t368 * t397 * t403 * (-t603 * t27 * t129 *  &
      t133 * t6 / 0.3e1_dp - t1101 * t1102 * t6 * t610 / 0.6e1_dp - t604 * t1096  &
      * t132 - t416 * t1096 * t611 / 0.2e1_dp + t604 * t424 / 0.2e1_dp - t416  &
      * t406 * t1114 * t422 * t610 / 0.4e1_dp + t416 * t406 * t417 * (t844  &
      * t77 - t1121 * t608 - t1562 * t420 + 2._dp * t1124 * t1482 -  &
      t419 * t77 * t871) / 0.2e1_dp)
      t1583 = t621 * t142
      t1592 = -t1578 * t140 * t145 - 0.7e1_dp / 0.2e1_dp * t1143 * t623 - 0.7e1_dp  &
      / 0.2e1_dp * t1583 * t438 - 0.35e2_dp / 0.4e1_dp * t435 * t1147 * t1446 -  &
      0.7e1_dp / 0.2e1_dp * t435 * t437 * t875
      t1594 = t1592 * E * t149
      t1608 = d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t1620 = t199 * t320
      t1627 = t627 * t444
      t1632 = t156 * t875
      t1645 = t627 * t531
      t1652 = -t1594 * t195 - t656 * t493 / 0.3e1_dp - (3._dp * t442 * t629)  &
      - 0.25e2_dp / 0.6e1_dp * t1288 * t258 * r3 * t5 * t156 * t320 - (3._dp  &
      * t627 * t445) + t197 * (t1608 + (t875 * t200 - 0.2e1_dp * t285 * t545  &
      * t320 + 0.2e1_dp * t1393 * t285 * t320 - t546 * t875) * t549 * t199  &
      - t1402 * t1620 + t550 * t320) + 0.5e1_dp / 0.2e1_dp * t1162 * t684 +  &
      (3._dp * t1627 * t528) + (3._dp * t973 * t681) + t1594 + 0.5e1_dp / 0.2e1_dp  &
      * t532 * t533 * t1632 - t677 * t523 + 0.5e1_dp / 0.3e1_dp * (t627)  &
      * (t500) * (t539) - 0.12e2_dp * t148 * t968 * t159 * t533  &
      * t1446 + 0.5e1_dp / 0.2e1_dp * t1645 * t534 + 0.2e1_dp * t81 * t149 * t325  &
      * t320
      t1675 = t1241 * t651
      t1678 = t151 * t626
      t1685 = t151 * t152
      t1687 = t461 * t320
      t1688 = t1213 * t1210 * t1687
      t1692 = t483 * t155 * t1687
      t1695 = t1203 * t637
      t1712 = t155 * t875
      t1719 = 0.3e1_dp / 0.2e1_dp * t1326 * t652 - 0.3e1_dp / 0.2e1_dp * t1331 * t1675  &
      + t1678 * t453 - (2._dp * t451 * t500 * t1446) + (t451 * t452  &
      * t875) - 0.75e2_dp / 0.4e1_dp * t1685 * t1688 - 0.3e1_dp / 0.2e1_dp * t1331  &
      * t1692 - 0.5e1_dp / 0.2e1_dp * t1363 * t1695 + t1350 * t635 - 0.3e1_dp /  &
      0.4e1_dp * t84 * t1315 * t461 * t320 + 0.5e1_dp / 0.2e1_dp * t1678 * t463 -  &
      t1308 * t470 * t461 * t320 / 0.2e1_dp - t1308 * t1309 * t320 / 0.2e1_dp  &
      + 0.3e1_dp / 0.2e1_dp * t478 * t484 * t1712 - 0.2e1_dp * t478 * t512 * (t1446)
      t1731 = t166 * t567
      t1737 = t166 * t167
      t1739 = t483 * t154 * t1687
      t1747 = t1250 * t156 * t1687
      t1753 = t456 * t156 * t1687
      t1762 = t456 * t168 * t1687
      t1765 = (t478 * t479 * t875) - t166 * t1437 * t83 * t170 + t1326  &
      * t649 - (2._dp * t81 * t149 * t163 * t325 * t320) + 0.3e1_dp / 0.2e1_dp  &
      * t1731 * t486 + 0.5e1_dp / 0.2e1_dp * t451 * t457 * t1632 + (3._dp *  &
      t1737 * t1739) - t151 * t1592 * t83 * t159 + t1731 * t480 - 0.27e2_dp  &
      / 0.4e1_dp * (t1737) * (t1747) + (t81 * t466 * t875) -  &
      0.5e1_dp / 0.2e1_dp * t1363 * t1753 + 0.5e1_dp / 0.2e1_dp * t1350 * t638 + (t81  &
      * t471 * t875) / 0.2e1_dp + (10._dp * t1685 * t1762)
      t1781 = t176 * t626
      t1796 = t176 * t177
      t1804 = -t176 * t1592 * t93 * t159 + (2._dp * t1193 * t661) + 0.5e1_dp  &
      / 0.2e1_dp * (t1193) * (t664) + (2._dp * t1781 * t501) - (6._dp  &
      * t499 * t1178 * t1446) - (5._dp * t1202 * t1695) + (2._dp * t499  &
      * t500 * t875) + 0.5e1_dp / 0.2e1_dp * (t1781) * (t505) - (5._dp  &
      * t1202 * t1753) - 0.75e2_dp / 0.4e1_dp * t1796 * t1688 + 0.10e2_dp *  &
      t1796 * t1762 + 0.5e1_dp / 0.2e1_dp * (t499) * (t504) * (t1632)
      t1812 = t180 * t567
      t1834 = -t180 * t1437 * t93 * t170 + (2._dp * t1230 * t670) + 0.3e1_dp  &
      / 0.2e1_dp * (t1230) * (t673) + (2._dp * t1812 * t513) - (6._dp  &
      * t511 * t1235 * t1446) - (3._dp * t1240 * t1675) + (2._dp * t511  &
      * t512 * t875) + 0.3e1_dp / 0.2e1_dp * (t1812) * (t517) - (3._dp  &
      * t1240 * t1692) - 0.27e2_dp / 0.4e1_dp * t95 * t1747 + 0.3e1_dp * t95  &
      * t1739 + 0.3e1_dp / 0.2e1_dp * (t511) * (t516) * (t1712)
      t1860 = (6._dp * t511 * t445 * t320) + (3._dp * t526 * t194 * t44 *  &
      t875) - (t81 * t93 * t875) - (3._dp * t148 * t444 * t875) - (2._dp  &
      * t568 * t345) + f12 * t1437 * t94 - (5._dp * t1179 * t538 *  &
      t6 * t320) + (t1719 + t1765) * omega * t128 - 0.75e2_dp / 0.4e1_dp * (t1273)  &
      * (t533) * (t1210) * (t461) * (t320) + (t1804  &
      + t1834) * t185 * t190 - 0.15e2_dp / 0.2e1_dp * (t1169) * (t533)  &
      * (t462) * (t320) + (12._dp * t148 * t968 * t325 *  &
      t320) - (2._dp * t342 * t570) - 0.15e2_dp / 0.2e1_dp * (t1169) * (t533)  &
      * (t1170) * (t320) + (10._dp * t532 * t533 * t168  &
      * t461 * t320) - (2._dp * t91 * t344 * t875)
      e_ndrho_rho = e_ndrho_rho + ( -0.4e1_dp / 0.3e1_dp * t322 * t696 - t80 * (t1652 + t1860) * Clda ) * sx
      t1878 = 2._dp * t119 * t88 + 4._dp * t288 * t565 + t11 * t15 * (F2 * t892 *  &
      t68 - 2._dp * t561 * t317 + 2._dp * t85 * t900 - t85 * t917)
      t1886 = t320 ** 2
      t1935 = t599 ** 2
      t1945 = t610 ** 2
      t1979 = -((t97 * (t99 * t1878 * t82 + 2._dp * t99 * t567 * t320 + t99  &
      * t90 * t920 + 2._dp * t102 * t1886 + 2._dp * t102 * t82 * t920 + 6._dp * t104  &
      * t82 * t1886 + 3._dp * t104 * t92 * t920) * t113) - (7._dp * t584  &
      * t586) + 0.63e2_dp / 0.4e1_dp * (t108) * (t1005) * (t1886)  &
      - 0.7e1_dp / 0.2e1_dp * (t108) * (t364) * (t920) - (t368  &
      * (f94 * t892 * t117 * t121 - 2._dp * t1480 * t593 + 4._dp * t589 * t595  &
      * t121 + 2._dp * t1025 * t18 * t375 * t899 - 4._dp * t1488 * t593 - t374  &
      * t18 * t375 * t916 + 2._dp * t116 * t68 * t3 * t18 * t375) * t123 * t136)  &
      - (t368 * t1935 * t123 * t136) + (4._dp * t1522 * t617) +  &
      (2._dp * t1068 * t1535 * t616) + 0.2e1_dp * (t368) * t397 * t403  &
      * (t604 * t612 - t416 * t406 * t1114 * t1945 / 0.4e1_dp + t416 * t406  &
      * t417 * (t893 * t77 - 2._dp * t1562 * t608 + 2._dp * t1124 * t77 *  &
      t899 - t419 * t77 * t916) / 0.2e1_dp)) * t140 * t145 - (7._dp * t1583  &
      * t623) - 0.35e2_dp / 0.4e1_dp * (t435) * (t1147) * (t1886)  &
      - 0.7e1_dp / 0.2e1_dp * (t435) * (t437) * (t920)
      t1981 = t1979 * E * t149
      t1989 = t1886 * t156
      t1999 = t168 * t1886
      t2003 = t156 * t920
      t2007 = t1210 * t1886
      t2013 = -t1981 * t195 + t1981 + (6._dp * t1627 * t681) + (3._dp * t526  &
      * t194 * t44 * t920) - (15._dp * t1169 * t533 * t1989) + (5._dp  &
      * t1645 * t684) - (12._dp * t1295 * t194 * t44 * t1886) + (10._dp  &
      * t532 * t533 * t1999) + 0.5e1_dp / 0.2e1_dp * (t532) * (t533)  &
      * (t2003) - 0.75e2_dp / 0.4e1_dp * (t1273) * (t533) * (t2007)  &
      - (4._dp * t568 * t570)
      t2050 = t1886 * t155
      t2060 = t154 * t1886
      t2064 = t155 * t920
      t2068 = -t176 * t1979 * t93 * t159 + (4._dp * t1781 * t661) + (5._dp  &
      * t1781 * t664) - (6._dp * t499 * t1178 * t1886) - (10._dp * t499  &
      * t531 * t1989) + (2._dp * t499 * t500 * t920) - 0.75e2_dp / 0.4e1_dp  &
      * (t499) * (t1214) * (t2007) + (10._dp * t499 * t504 *  &
      t1999) + 0.5e1_dp / 0.2e1_dp * (t499) * (t504) * (t2003) -  &
      t180 * t1878 * t93 * t170 + (4._dp * t1812 * t670) + (3._dp * t1812  &
      * t673) - (6._dp * t511 * t1235 * t1886) - (6._dp * t511 * t149  &
      * t483 * t2050) + (2._dp * t511 * t512 * t920) - 0.27e2_dp / 0.4e1_dp *  &
      (t511) * (t1251) * (t1989) + (3._dp * t511 * t516 * t2060)  &
      + 0.3e1_dp / 0.2e1_dp * (t511) * (t516) * (t2064)
      t2107 = -(2._dp * t451 * t500 * t1886) + (t451 * t452 * t920) -  &
      (5._dp * t451 * t504 * t1989) + (5._dp * t1678 * t638) + (10._dp  &
      * t451 * t457 * t1999) + 0.5e1_dp / 0.2e1_dp * (t451) * (t457) *  &
      (t2003) - 0.75e2_dp / 0.4e1_dp * (t451) * (t1300) * (t2007)  &
      - (2._dp * t81 * t1304 * t1886) - (t81 * t93 * t470 * t1886)  &
      + (t81 * t466 * t920) - 0.3e1_dp / 0.4e1_dp * (t81) * (t1316)  &
      * (t1886)
      t2140 = t81 * t471 * t920 / 0.2e1_dp - t166 * t1878 * t83 * t170 + (2._dp  &
      * t1731 * t649) + (3._dp * t1731 * t652) - (2._dp * t478 * t512  &
      * t1886) + (t478) * t479 * t920 - (3._dp * t478 * t516 * t2050)  &
      + (3._dp * t478 * t484 * t2060) + 0.3e1_dp / 0.2e1_dp * (t478) *  &
      (t484) * (t2064) - 0.27e2_dp / 0.4e1_dp * (t478) * (t1346)  &
      * (t1989) - t151 * t1979 * t83 * t159 + (2._dp * t1678 * t635)
      t2159 = d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
      t2174 = t2068 * t185 * t190 + 6._dp * t91 * t1382 * t1886 - 2._dp * t91 * t344  &
      * t920 + (t2107 + t2140) * omega * t128 + f12 * t1878 * t94 -  &
      t81 * t93 * t920 - 6._dp * t627 * t629 - 3._dp * t148 * t444 * t920 + 2._dp *  &
      t81 * t149 * t1886 + 12._dp * t148 * t968 * t1886 + t197 * (t2159 + (t920  &
      * t200 - 2._dp * t1886 * t545 + 2._dp * t1393 * t1886 - t546 * t920) *  &
      t549 * t199 - t690 * t1401 * t1620 + t691 * t320)
      e_ndrho_ndrho = e_ndrho_ndrho + ( -t80 * (t2013 + t2174) * Clda ) * sx
    END IF

  END SUBROUTINE xwpbe_lda_calc_3

! *****************************************************************************
!> \brief Evaluates the screened hole averaged PBE exchange functional for lda.
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param rho , ndrho: density and norm of the density gradient
!> \param ndrho ...
!> \param omega screening parameter
!> \param sscale scaling factor to enforce Lieb-Oxford bound
!> \param sx scaling factor
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine evaluates the functional for omega!=0 using a simple
!>      gaussian expansion for large ww and a taylor expansion for the
!>      parameter G.
! *****************************************************************************
  SUBROUTINE xwpbe_lda_calc_4(e_0, e_rho, e_ndrho, e_rho_rho, e_ndrho_rho, &
                              e_ndrho_ndrho, rho, ndrho, omega, sscale, sx, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho, &
                                                e_rho_rho, e_ndrho_rho, &
                                                e_ndrho_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, omega, sscale, sx
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: d2Qndrhondrho, d2Qrhondrho, d2Qrhorho, dQndrho, dQrho, &
      Q, t1, t10, t100, t1001, t1011, t1017, t102, t1026, t103, t1033, t1035, &
      t1040, t1047, t1048, t1049, t105, t106, t1065, t1066, t1071, t1074, &
      t108, t1082, t1089, t109, t1098, t11, t110, t111, t1111, t1118, t113, &
      t114, t115, t1155, t1157, t116, t117, t1174, t118, t1181, t1184, t1189, &
      t1190, t1198, t12, t1205, t1208, t121, t1210, t1218, t122, t1224, &
      t1231, t125, t1255, t126, t1261, t1264, t1266, t127, t1270, t1277, &
      t128, t1288, t129, t1299, t13, t1319, t1324, t133, t1336, t134, t1351, &
      t136, t137, t1382, t1397, t14, t140
    REAL(KIND=dp) :: t1405, t141, t1413, t1435, t1443, t1447, t1448, t145, &
      t1452, t146, t1481, t1483, t149, t15, t1500, t151, t1529, t153, t1533, &
      t1537, t154, t155, t1552, t156, t1562, t1566, t1570, t1576, t158, t159, &
      t16, t160, t161, t1618, t163, t1652, t167, t1672, t169, t17, t170, &
      t171, t173, t174, t176, t177, t178, t179, t18, t180, t182, t183, t184, &
      t185, t188, t189, t19, t192, t193, t195, t196, t197, t198, t199, t2, &
      t203, t206, t207, t208, t209, t21, t210, t211, t212, t213, t219, t22, &
      t220, t221, t222, t223, t227, t229, t23, t235, t238, t239, t24, t240, &
      t241, t246, t248, t249, t25, t252
    REAL(KIND=dp) :: t253, t256, t257, t258, t261, t262, t263, t265, t266, &
      t269, t27, t270, t274, t277, t278, t279, t28, t281, t283, t286, t289, &
      t29, t290, t293, t296, t299, t3, t300, t302, t303, t305, t306, t309, &
      t31, t310, t313, t316, t32, t321, t326, t327, t329, t330, t331, t332, &
      t333, t334, t34, t340, t341, t342, t345, t346, t35, t350, t351, t352, &
      t355, t359, t36, t360, t367, t368, t369, t372, t373, t374, t375, t379, &
      t38, t382, t388, t389, t39, t390, t393, t394, t4, t400, t401, t402, &
      t405, t406, t41, t410, t412, t415, t417, t42, t420, t421, t422, t423, &
      t426, t427, t428, t431, t433, t434
    REAL(KIND=dp) :: t435, t437, t438, t439, t44, t443, t444, t450, t453, &
      t454, t456, t457, t459, t46, t464, t465, t468, t469, t472, t473, t475, &
      t476, t477, t478, t48, t484, t486, t487, t49, t498, t5, t500, t501, &
      t505, t510, t513, t519, t522, t526, t530, t533, t536, t539, t54, t540, &
      t544, t545, t548, t55, t553, t555, t557, t558, t56, t561, t563, t564, &
      t568, t569, t572, t575, t576, t578, t579, t58, t581, t584, t588, t594, &
      t597, t599, t6, t60, t603, t607, t61, t610, t613, t616, t618, t619, &
      t620, t623, t63, t632, t65, t655, t657, t662, t664, t667, t67, t68, &
      t69, t694, t7, t70, t700, t706, t71, t72
    REAL(KIND=dp) :: t723, t726, t728, t73, t74, t744, t75, t751, t752, t769, &
      t77, t772, t774, t78, t782, t783, t784, t789, t792, t793, t794, t795, &
      t799, t8, t80, t803, t804, t807, t81, t811, t812, t819, t82, t83, t84, &
      t848, t85, t852, t862, t863, t864, t865, t868, t87, t872, t878, t879, &
      t88, t880, t9, t90, t91, t916, t92, t920, t93, t930, t931, t932, t935, &
      t939, t94, t943, t95, t956, t96, t961, t966, t97, t972, t985, t99, &
      t990, t995

    IF( order >= 0 ) THEN
      t1 = ndrho ** 2
      t2 = r2 ** 2
      t3 = 0.1e1_dp / t2
      t4 = t1 * t3
      t5 = pi ** 2
      t6 = r3 * t5
      t7 = t6 * rho
      t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
      t9 = t8 ** 2
      t10 = 0.1e1_dp / t9
      t11 = t4 * t10
      t12 = rho ** 2
      t13 = 0.1e1_dp / t12
      t14 = sscale ** 2
      t15 = t13 * t14
      t16 = a1 * t1
      t17 = t16 * t3
      t18 = t10 * t13
      t19 = t18 * t14
      t21 = t1 ** 2
      t22 = a2 * t21
      t23 = t2 ** 2
      t24 = 0.1e1_dp / t23
      t25 = t22 * t24
      t27 = 0.1e1_dp / t8 / t7
      t28 = t12 ** 2
      t29 = 0.1e1_dp / t28
      t31 = t14 ** 2
      t32 = t27 * t29 * t31
      t34 = t17 * t19 + t25 * t32
      t35 = a3 * t21
      t36 = t35 * t24
      t38 = t21 * ndrho
      t39 = a4 * t38
      t41 = 0.1e1_dp / t23 / r2
      t42 = t39 * t41
      t44 = 0.1e1_dp / t9 / t7
      t46 = 0.1e1_dp / t28 / rho
      t48 = t31 * sscale
      t49 = t44 * t46 * t48
      t54 = 0.1e1_dp / t23 / t2
      t55 = a5 * t21 * t1 * t54
      t56 = r3 ** 2
      t58 = t5 ** 2
      t60 = 0.1e1_dp / t56 / t58
      t61 = t28 ** 2
      t63 = t31 * t14
      t65 = t60 / t61 * t63
      t67 = r1 + t36 * t32 + t42 * t49 + t55 * t65
      t68 = 0.1e1_dp / t67
      t69 = t34 * t68
      t70 = t15 * t69
      t71 = t11 * t70
      t72 = omega ** 2
      t73 = beta2 * t72
      t74 = t73 * t10
      t75 = t71 + t74
      t77 = 0.1e1_dp / A
      Q = f94 * t75 * t77
      t78 = rho ** (0.1e1_dp / 0.3e1_dp)
      t80 = t78 * rho * f89
      t81 = B * f12
      t82 = t71 + DD
      t83 = 0.1e1_dp / t82
      t84 = t81 * t83
      t85 = F2 * t34
      t87 = F1 + t85 * t68
      t88 = t15 * t87
      t90 = t11 * t88 + r1
      t91 = f12 * t90
      t92 = t82 ** 2
      t93 = 0.1e1_dp / t92
      t94 = C * t93
      t95 = t91 * t94
      t96 = g2 * t1
      t97 = t96 * t3
      t99 = g3 * t21
      t100 = t99 * t24
      t102 = g1 + t97 * t19 + t100 * t32
      t103 = t15 * t102
      t105 = t11 * t103 + r1
      t106 = t105 * E
      t108 = 0.1e1_dp / t92 / t82
      t109 = t106 * t108
      t110 = f158 * E
      t111 = t105 * t83
      t113 = t71 + DD + t72 * t10
      t114 = t113 ** 2
      t115 = t114 ** 2
      t116 = t115 * t113
      t117 = SQRT(t116)
      t118 = 0.1e1_dp / t117
      t121 = SQRT(t113)
      t122 = 0.1e1_dp / t121
      t125 = f68 * C
      t126 = t90 * t83
      t127 = t114 * t113
      t128 = SQRT(t127)
      t129 = 0.1e1_dp / t128
      t133 = (-t110 * t111 * t118 - t81 * t83 * t122 - t125 * t126 * t129)  &
      * omega
      t134 = 0.1e1_dp / t8
      t136 = f52 * E
      t137 = t105 * t93
      t140 = f12 * C
      t141 = t90 * t93
      t145 = t72 * omega
      t146 = (-t136 * t137 * t118 - t140 * t141 * t129) * t145
      t149 = 0.1e1_dp / r3 / t5
      t151 = t149 / rho
      t153 = t72 ** 2
      t154 = t153 * omega
      t155 = t118 * t154
      t156 = t155 * t44
      t158 = f12 * A
      t159 = exei(Q)
      t160 = t71 + DD + t74
      t161 = 0.1e1_dp / t160
      t163 = LOG(t75 * t161)
      t167 = (t84 + t95 + t109 + t133 * t134 + t146 * t151 - t109 * t156  &
      + t158 * (t159 + t163)) * Clda
      e_0 = e_0 + ( -t80 * t167 )  * sx
    END IF
    IF( order >= 1 .OR. order == -1 ) THEN
      t169 = t44 * t13
      t170 = t4 * t169
      t171 = t14 * t34
      t173 = t68 * r3 * t5
      t174 = t171 * t173
      t176 = 0.2e1_dp / 0.3e1_dp * t170 * t174
      t177 = t12 * rho
      t178 = 0.1e1_dp / t177
      t179 = t178 * t14
      t180 = t179 * t69
      t182 = 2._dp * t11 * t180
      t183 = t3 * t44
      t184 = t16 * t183
      t185 = t15 * t6
      t188 = t10 * t178
      t189 = t188 * t14
      t192 = t56 * t58
      t193 = t192 * t12
      t195 = 0.1e1_dp / t8 / t193
      t196 = t24 * t195
      t197 = t22 * t196
      t198 = t29 * t31
      t199 = t198 * t6
      t203 = t27 * t46 * t31
      t206 = -0.2e1_dp / 0.3e1_dp * t184 * t185 - (2._dp * t17 * t189) - 0.4e1_dp  &
      / 0.3e1_dp * t197 * t199 - (4._dp * t25 * t203)
      t207 = t206 * t68
      t208 = t15 * t207
      t209 = t11 * t208
      t210 = t4 * t18
      t211 = t67 ** 2
      t212 = 0.1e1_dp / t211
      t213 = t35 * t196
      t219 = 0.1e1_dp / t9 / t193
      t220 = t41 * t219
      t221 = t39 * t220
      t222 = t46 * t48
      t223 = t222 * t6
      t227 = 0.1e1_dp / t28 / t12
      t229 = t44 * t227 * t48
      t235 = t60 / t61 / rho * t63
      t238 = -0.4e1_dp / 0.3e1_dp * t213 * t199 - (4._dp * t36 * t203) - 0.5e1_dp  &
      / 0.3e1_dp * t221 * t223 - (5._dp * t42 * t229) - (8._dp * t55 * t235)
      t239 = t212 * t238
      t240 = t171 * t239
      t241 = t210 * t240
      t246 = -t176 - t182 + t209 - t241 - 0.2e1_dp / 0.3e1_dp * t73 * t44 * r3  &
      * t5
      dQrho = f94 * t246 * t77
      t248 = ndrho * t3
      t249 = t248 * t10
      t252 = a1 * ndrho
      t253 = t252 * t3
      t256 = t1 * ndrho
      t257 = a2 * t256
      t258 = t257 * t24
      t261 = 2._dp * t253 * t19 + 4._dp * t258 * t32
      t262 = t261 * t68
      t263 = t15 * t262
      t265 = a3 * t256
      t266 = t265 * t24
      t269 = a4 * t21
      t270 = t269 * t41
      t274 = a5 * t38 * t54
      t277 = 4._dp * t266 * t32 + 5._dp * t270 * t49 + 6._dp * t274 * t65
      t278 = t212 * t277
      t279 = t171 * t278
      t281 = 2._dp * t249 * t70 + t11 * t263 - t210 * t279
      dQndrho = f94 * t281 * t77
      t283 = t78 * f89
      t286 = -t176 - t182 + t209 - t241
      t289 = t14 * t87
      t290 = t289 * t6
      t293 = t179 * t87
      t296 = F2 * t206
      t299 = t296 * t68 - t85 * t239
      t300 = t15 * t299
      t302 = -0.2e1_dp / 0.3e1_dp * t170 * t290 - (2._dp * t11 * t293) + (t11  &
      * t300)
      t303 = f12 * t302
      t305 = C * t108
      t306 = t305 * t286
      t309 = t14 * t102
      t310 = t309 * t6
      t313 = t179 * t102
      t316 = t96 * t183
      t321 = t99 * t196
      t326 = -0.2e1_dp / 0.3e1_dp * t316 * t185 - (2._dp * t97 * t189) - 0.4e1_dp  &
      / 0.3e1_dp * t321 * t199 - (4._dp * t100 * t203)
      t327 = t15 * t326
      t329 = -0.2e1_dp / 0.3e1_dp * t170 * t310 - (2._dp * t11 * t313) + (t11  &
      * t327)
      t330 = t329 * E
      t331 = t330 * t108
      t332 = t92 ** 2
      t333 = 0.1e1_dp / t332
      t334 = t333 * t286
      t340 = t110 * t105
      t341 = t93 * t118
      t342 = t341 * t286
      t345 = 0.1e1_dp / t117 / t116
      t346 = t83 * t345
      t350 = -t176 - t182 + t209 - t241 - 0.2e1_dp / 0.3e1_dp * t72 * t44 * t6
      t351 = t115 * t350
      t352 = t346 * t351
      t355 = t93 * t122
      t359 = 0.1e1_dp / t121 / t113
      t360 = t83 * t359
      t367 = t125 * t90
      t368 = t93 * t129
      t369 = t368 * t286
      t372 = 0.1e1_dp / t128 / t127
      t373 = t83 * t372
      t374 = t114 * t350
      t375 = t373 * t374
      t379 = (-t110 * t329 * t83 * t118 + t340 * t342 + 0.5e1_dp / 0.2e1_dp * t340  &
      * t352 + t81 * t355 * t286 + t81 * t360 * t350 / 0.2e1_dp - t125  &
      * t302 * t83 * t129 + t367 * t369 + 0.3e1_dp / 0.2e1_dp * t367 * t375) *  &
      omega
      t382 = t27 * r3 * t5
      t388 = t136 * t105
      t389 = t108 * t118
      t390 = t389 * t286
      t393 = t93 * t345
      t394 = t393 * t351
      t400 = t140 * t90
      t401 = t108 * t129
      t402 = t401 * t286
      t405 = t93 * t372
      t406 = t405 * t374
      t410 = (-t136 * t329 * t93 * t118 + (2._dp * t388 * t390) + 0.5e1_dp /  &
      0.2e1_dp * (t388) * (t394) - t140 * t302 * t93 * t129 + (2._dp  &
      * t400 * t402) + 0.3e1_dp / 0.2e1_dp * (t400) * (t406)) * t145
      t412 = t149 * t13
      t415 = t106 * t333
      t417 = t155 * t44 * t286
      t420 = t108 * t345
      t421 = t106 * t420
      t422 = t154 * t44
      t423 = t422 * t351
      t426 = t106 * t389
      t427 = t154 * t219
      t428 = t427 * t6
      t431 = dexeirho(Q,dQrho)
      t433 = t160 ** 2
      t434 = 0.1e1_dp / t433
      t435 = t75 * t434
      t437 = t246 * t161 - t435 * t246
      t438 = 0.1e1_dp / t75
      t439 = t437 * t438
      t443 = -t81 * t93 * t286 + t303 * t94 - (2._dp * t91 * t306) + t331  &
      - (3._dp * t106 * t334) + t379 * t134 - t133 * t382 / 0.3e1_dp + t410  &
      * t151 - t146 * t412 - t331 * t156 + (3._dp * t415 * t417) + 0.5e1_dp  &
      / 0.2e1_dp * t421 * t423 + 0.5e1_dp / 0.3e1_dp * t426 * t428 + t158 * (t431  &
      + t439 * t160)
      t444 = t443 * Clda
      e_rho = e_rho + ( -0.4e1_dp / 0.3e1_dp * t283 * t167 - t80 * t444 ) * sx
      t450 = F2 * t261
      t453 = t450 * t68 - t85 * t278
      t454 = t15 * t453
      t456 = 2._dp * t249 * t88 + t11 * t454
      t457 = f12 * t456
      t459 = t305 * t281
      t464 = g2 * ndrho
      t465 = t464 * t3
      t468 = g3 * t256
      t469 = t468 * t24
      t472 = 2._dp * t465 * t19 + 4._dp * t469 * t32
      t473 = t15 * t472
      t475 = 2._dp * t249 * t103 + t11 * t473
      t476 = t475 * E
      t477 = t476 * t108
      t478 = t333 * t281
      t484 = t341 * t281
      t486 = t115 * t281
      t487 = t346 * t486
      t498 = t368 * t281
      t500 = t114 * t281
      t501 = t373 * t500
      t505 = (-t110 * t475 * t83 * t118 + t340 * t484 + 0.5e1_dp / 0.2e1_dp * t340  &
      * t487 + t81 * t355 * t281 + t81 * t360 * t281 / 0.2e1_dp - t125  &
      * t456 * t83 * t129 + t367 * t498 + 0.3e1_dp / 0.2e1_dp * t367 * t501) *  &
      omega
      t510 = t389 * t281
      t513 = t393 * t486
      t519 = t401 * t281
      t522 = t405 * t500
      t526 = (-t136 * t475 * t93 * t118 + (2._dp * t388 * t510) + 0.5e1_dp /  &
      0.2e1_dp * (t388) * (t513) - t140 * t456 * t93 * t129 + (2._dp  &
      * t400 * t519) + 0.3e1_dp / 0.2e1_dp * (t400) * (t522)) * t145
      t530 = t155 * t44 * t281
      t533 = t422 * t486
      t536 = dexeindrho(Q,dQndrho)
      t539 = t281 * t161 - t435 * t281
      t540 = t539 * t438
      t544 = -t81 * t93 * t281 + t457 * t94 - (2._dp * t91 * t459) + t477  &
      - (3._dp * t106 * t478) + t505 * t134 + t526 * t151 - t477 * t156  &
      + (3._dp * t415 * t530) + 0.5e1_dp / 0.2e1_dp * t421 * t533 + t158 * (t536  &
      + t540 * t160)
      t545 = t544 * Clda
      e_ndrho = e_ndrho + ( -t80 * t545 ) * sx
    END IF
    IF( order >= 2 .OR. order == -2 ) THEN
      t548 = t4 * t219 * t13
      t553 = 0.10e2_dp / 0.9e1_dp * t548 * t171 * t68 * t56 * t58
      t555 = t4 * t44 * t178
      t557 = 0.8e1_dp / 0.3e1_dp * t555 * t174
      t558 = t14 * t206
      t561 = 0.4e1_dp / 0.3e1_dp * t170 * t558 * t173
      t563 = t4 * t169 * t14
      t564 = t34 * t212
      t568 = 0.4e1_dp / 0.3e1_dp * t563 * t564 * t6 * t238
      t569 = t29 * t14
      t572 = 6._dp * t11 * t569 * t69
      t575 = 4._dp * t11 * t179 * t207
      t576 = t4 * t188
      t578 = 4._dp * t576 * t240
      t579 = t3 * t219
      t581 = t15 * t192
      t584 = t179 * t6
      t588 = t10 * t29 * t14
      t594 = t56 * r3 * t58 * t5 * t177
      t597 = t24 / t8 / t594
      t599 = t198 * t192
      t603 = t46 * t31 * t6
      t607 = t27 * t227 * t31
      t610 = 0.10e2_dp / 0.9e1_dp * t16 * t579 * t581 + 0.8e1_dp / 0.3e1_dp * t184 *  &
      t584 + (6._dp * t17 * t588) + 0.28e2_dp / 0.9e1_dp * t22 * t597 * t599 +  &
      0.32e2_dp / 0.3e1_dp * t197 * t603 + (20._dp * t25 * t607)
      t613 = t11 * t15 * t610 * t68
      t616 = 2._dp * t210 * t558 * t239
      t618 = 0.1e1_dp / t211 / t67
      t619 = t238 ** 2
      t620 = t618 * t619
      t623 = 2._dp * t210 * t171 * t620
      t632 = 0.1e1_dp / t9 / t594
      t655 = t212 * (0.28e2_dp / 0.9e1_dp * t35 * t597 * t599 + 0.32e2_dp / 0.3e1_dp  &
      * t213 * t603 + (20._dp * t36 * t607) + 0.40e2_dp / 0.9e1_dp * t39 * t41  &
      * t632 * t222 * t192 + 0.50e2_dp / 0.3e1_dp * t221 * t227 * t48 * t6 +  &
      0.30e2_dp * t42 * t44 / t28 / t177 * t48 + (72._dp * t55 * t60 / t61  &
      / t12 * t63))
      t657 = t210 * t171 * t655
      t662 = t553 + t557 - t561 + t568 + t572 - t575 + t578 + t613 - t616  &
      + t623 - t657 + 0.10e2_dp / 0.9e1_dp * t73 * t219 * t56 * t58
      d2Qrhorho = f94 * t662 * t77
      t664 = t248 * t169
      t667 = t14 * t261
      t694 = -0.4e1_dp / 0.3e1_dp * t252 * t183 * t185 - (4._dp * t253 * t189)  &
      - 0.16e2_dp / 0.3e1_dp * t257 * t196 * t199 - (16._dp * t258 * t203)
      t700 = t248 * t18
      t706 = t618 * t238 * t277
      t723 = t212 * (-0.16e2_dp / 0.3e1_dp * t265 * t196 * t199 - (16._dp * t266  &
      * t203) - 0.25e2_dp / 0.3e1_dp * t269 * t220 * t223 - (25._dp * t270 *  &
      t229) - (48._dp * t274 * t235))
      t726 = -0.4e1_dp / 0.3e1_dp * t664 * t174 - 0.2e1_dp / 0.3e1_dp * t170 * t667 *  &
      t173 + 0.2e1_dp / 0.3e1_dp * t563 * t564 * t6 * t277 - (4._dp * t249 *  &
      t180) - (2._dp * t11 * t179 * t262) + (2._dp * t576 * t279) + (2._dp  &
      * t249 * t208) + (t11 * t15 * t694 * t68) - t210 * t558 * t278  &
      - (2._dp * t700 * t240) - t210 * t667 * t239 + 0.2e1_dp * t210 *  &
      t171 * t706 - t210 * t171 * t723
      d2Qrhondrho = f94 * t726 * t77
      t728 = t3 * t10
      t744 = 2._dp * a1 * t3 * t19 + 12._dp * a2 * t1 * t24 * t32
      t751 = t277 ** 2
      t752 = t618 * t751
      t769 = t212 * (12._dp * a3 * t1 * t24 * t32 + 20._dp * a4 * t256 * t41 * t49  &
      + 30._dp * a5 * t21 * t54 * t65)
      t772 = 2._dp * t728 * t13 * t171 * t68 + 4._dp * t249 * t263 - 4._dp * t700 * t279  &
      + t11 * t15 * t744 * t68 - 2._dp * t210 * t667 * t278 + 2._dp * t210 *  &
      t171 * t752 - t210 * t171 * t769
      d2Qndrhondrho = f94 * t772 * t77
      t774 = t78 ** 2
      t782 = 0.1e1_dp / t332 / t82
      t783 = t106 * t782
      t784 = t286 ** 2
      t789 = t115 ** 2
      t792 = 0.1e1_dp / t117 / t789 / t114
      t793 = t83 * t792
      t794 = t350 ** 2
      t795 = t789 * t794
      t799 = t108 * t122
      t803 = t81 * t93
      t804 = t359 * t286
      t807 = t553 + t557 - t561 + t568 + t572 - t575 + t578 + t613 - t616  &
      + t623 - t657
      t811 = 0.1e1_dp / t121 / t114
      t812 = t83 * t811
      t819 = t553 + t557 - t561 + t568 + t572 - t575 + t578 + t613 - t616  &
      + t623 - t657 + 0.10e2_dp / 0.9e1_dp * t72 * t219 * t192
      t848 = 0.10e2_dp / 0.9e1_dp * t548 * t289 * t192 + 0.8e1_dp / 0.3e1_dp * t555 *  &
      t290 - 0.4e1_dp / 0.3e1_dp * t170 * t14 * t299 * t6 + (6._dp * t11 * t569  &
      * t87) - 0.4e1_dp * (t11) * t179 * t299 + (t11 * t15 * (F2  &
      * t610 * t68 - 2._dp * t296 * t239 + 2._dp * t85 * t620 - t85 * t655))
      t852 = t125 * t302
      t862 = -0.75e2_dp / 0.4e1_dp * t340 * t793 * t795 - (2._dp * t81 * t799 *  &
      t784) - t803 * t804 * t350 + (t81 * t355 * t807) - 0.3e1_dp / 0.4e1_dp  &
      * (t81) * (t812) * (t794) + (t81 * t360 * t819)  &
      / 0.2e1_dp - t125 * t848 * t83 * t129 + (2._dp * t852 * t369) + (3._dp  &
      * t852 * t375) - (2._dp * t367 * t401 * t784) + (t367 * t368  &
      * t807)
      t863 = t125 * t141
      t864 = t372 * t286
      t865 = t864 * t374
      t868 = t113 * t794
      t872 = t114 * t819
      t878 = 0.1e1_dp / t128 / t115 / t114
      t879 = t83 * t878
      t880 = t115 * t794
      t916 = 0.10e2_dp / 0.9e1_dp * t548 * t309 * t192 + 0.8e1_dp / 0.3e1_dp * t555 *  &
      t310 - 0.4e1_dp / 0.3e1_dp * t170 * t14 * t326 * t6 + (6._dp * t11 * t569  &
      * t102) - 0.4e1_dp * (t11) * t179 * t326 + (t11) * t15 * (0.10e2_dp  &
      / 0.9e1_dp * t96 * t579 * t581 + 0.8e1_dp / 0.3e1_dp * t316 * t584 +  &
      (6._dp * t97 * t588) + 0.28e2_dp / 0.9e1_dp * t99 * t597 * t599 + 0.32e2_dp  &
      / 0.3e1_dp * t321 * t603 + (20._dp * t100 * t607))
      t920 = t110 * t329
      t930 = t110 * t137
      t931 = t345 * t286
      t932 = t931 * t351
      t935 = t127 * t794
      t939 = t115 * t819
      t943 = -(3._dp * t863 * t865) + (3._dp * t367 * t373 * t868) + 0.3e1_dp  &
      / 0.2e1_dp * (t367) * (t373) * (t872) - 0.27e2_dp / 0.4e1_dp  &
      * (t367) * (t879) * (t880) - t110 * t916 * t83 * t118  &
      + (2._dp * t920 * t342) + (5._dp * t920 * t352) - (2._dp * t340 *  &
      t389 * t784) + (t340 * t341 * t807) - (5._dp * t930 * t932) +  &
      (10._dp * t340 * t346 * t935) + 0.5e1_dp / 0.2e1_dp * (t340) * (t346)  &
      * (t939)
      t956 = t136 * t329
      t961 = t333 * t118
      t966 = t136 * t105 * t108
      t972 = t93 * t792
      t985 = t140 * t302
      t990 = t333 * t129
      t995 = t140 * t90 * t108
      t1001 = t93 * t878
      t1011 = -t136 * t916 * t93 * t118 + (4._dp * t956 * t390) + (5._dp  &
      * t956 * t394) - (6._dp * t388 * t961 * t784) - (10._dp * t966 * t932)  &
      + (2._dp * t388 * t389 * t807) - 0.75e2_dp / 0.4e1_dp * (t388)  &
      * (t972) * (t795) + (10._dp * t388 * t393 * t935) + 0.5e1_dp  &
      / 0.2e1_dp * (t388) * (t393) * (t939) - t140 * t848 * t93  &
      * t129 + (4._dp * t985 * t402) + (3._dp * t985 * t406) - (6._dp *  &
      t400 * t990 * t784) - (6._dp * t995 * t865) + (2._dp * t400 * t401  &
      * t807) - 0.27e2_dp / 0.4e1_dp * (t400) * (t1001) * (t880)  &
      + (3._dp * t400 * t405 * t868) + 0.3e1_dp / 0.2e1_dp * (t400) * (t405)  &
      * (t872)
      t1017 = t106 * t961
      t1026 = t106 * t420 * t154
      t1033 = d2exeirhorho(Q,dQrho,d2Qrhorho)
      t1035 = t246 ** 2
      t1040 = t75 / t433 / t160
      t1047 = t75 ** 2
      t1048 = 0.1e1_dp / t1047
      t1049 = t437 * t1048
      t1065 = t106 * t333 * t345
      t1066 = t286 * t115
      t1071 = t330 * t420
      t1074 = -(12._dp * t783 * t155 * t44 * t784) + (t862 + t943) * omega  &
      * t134 + f12 * t848 * t94 + 0.4e1_dp / 0.9e1_dp * t133 * t195 * t56 *  &
      t58 + t1011 * t145 * t151 + (12._dp * t106 * t782 * t784) - (10._dp  &
      * t1017 * t427 * t286 * r3 * t5) - (t81 * t93 * t807) - 0.25e2_dp  &
      / 0.3e1_dp * (t1026) * (t219) * (t115) * (t350) *  &
      (r3) * (t5) + (t158 * (t1033 + (t662 * t161 - 2._dp * t1035  &
      * t434 + 2._dp * t1040 * t1035 - t435 * t662) * t438 * t160 - t1049  &
      * t160 * t246 + t439 * t246)) + (3._dp * t415 * t155 * t44 * t807)  &
      + (2._dp * t81 * t108 * t784) - (4._dp * t303 * t306) - (15._dp *  &
      t1065 * t422 * t1066 * t350) + (5._dp * t1071 * t423)
      t1082 = t106 * t108 * t792
      t1089 = t330 * t333
      t1098 = t916 * E * t108
      t1111 = C * t333
      t1118 = 0.5e1_dp / 0.2e1_dp * t421 * t422 * t939 + (2._dp * t146 * t149 *  &
      t178) - 0.75e2_dp / 0.4e1_dp * t1082 * t422 * t795 - (3._dp * t106 * t333  &
      * t807) + (6._dp * t1089 * t417) + 0.10e2_dp * t421 * t422 * t935  &
      - 0.2e1_dp / 0.3e1_dp * t379 * t382 + t1098 - 0.40e2_dp / 0.9e1_dp * t426 * t154  &
      * t632 * t192 - (2._dp * t410 * t412) - (2._dp * t91 * t305 * t807)  &
      - (6._dp * t330 * t334) - t1098 * t156 + (6._dp * t91 * t1111  &
      * t784) + 0.10e2_dp / 0.3e1_dp * (t330) * (t389) * (t428)
      e_rho_rho = e_rho_rho + ( -0.4e1_dp / 0.9e1_dp / t774 * f89 * t167 - 0.8e1_dp / 0.3e1_dp * t283 * t444  &
      - t80 * (t1074 + t1118) * Clda ) * sx
      t1155 = -0.4e1_dp / 0.3e1_dp * t664 * t310 - 0.2e1_dp / 0.3e1_dp * t170 * t14 *  &
      t472 * t6 - (4._dp * t249 * t313) - 0.2e1_dp * t11 * t179 * t472 + (2._dp  &
      * t249 * t327) + t11 * t15 * (-0.4e1_dp / 0.3e1_dp * t464 * t183 *  &
      t185 - (4._dp * t465 * t189) - 0.16e2_dp / 0.3e1_dp * t468 * t196 * t199  &
      - (16._dp * t469 * t203))
      t1157 = t1155 * E * t108
      t1174 = t476 * t333
      t1181 = t931 * t486
      t1184 = t114 * t726
      t1189 = t350 * t281
      t1190 = t372 * t114 * t1189
      t1198 = t115 * t726
      t1205 = t110 * t475
      t1208 = t110 * t111
      t1210 = t792 * t789 * t1189
      t1218 = t286 * t281
      t1224 = t125 * t456
      t1231 = -0.5e1_dp / 0.2e1_dp * t930 * t1181 + 0.3e1_dp / 0.2e1_dp * t367 * t373  &
      * t1184 - 0.3e1_dp / 0.2e1_dp * t863 * t1190 - t803 * t804 * t281 / 0.2e1_dp  &
      + t81 * t355 * t726 + 0.5e1_dp / 0.2e1_dp * t340 * t346 * t1198 + t81  &
      * t360 * t726 / 0.2e1_dp + 0.5e1_dp / 0.2e1_dp * t1205 * t352 - 0.75e2_dp / 0.4e1_dp  &
      * t1208 * t1210 - 0.2e1_dp * t81 * t108 * t122 * t286 * t281 - 0.2e1_dp  &
      * t340 * t389 * t1218 + t340 * t341 * t726 + 0.3e1_dp / 0.2e1_dp *  &
      t1224 * t375 - t110 * t1155 * t83 * t118 + t920 * t484
      t1255 = -0.4e1_dp / 0.3e1_dp * t664 * t290 - 0.2e1_dp / 0.3e1_dp * t170 * t14 *  &
      t453 * t6 - (4._dp * t249 * t293) - 0.2e1_dp * t11 * t179 * t453 + (2._dp  &
      * t249 * t300) + t11 * t15 * (F2 * t694 * t68 - t296 * t278  &
      - t450 * t239 + 2._dp * t85 * t706 - t85 * t723)
      t1261 = t345 * t115 * t1189
      t1264 = t125 * t126
      t1266 = t878 * t115 * t1189
      t1270 = t345 * t127 * t1189
      t1277 = t372 * t113 * t1189
      t1288 = t864 * t500
      t1299 = -t125 * t1255 * t83 * t129 + t852 * t498 - 0.5e1_dp / 0.2e1_dp *  &
      t930 * t1261 - 0.27e2_dp / 0.4e1_dp * t1264 * t1266 + (10._dp * t1208 *  &
      t1270) + 0.3e1_dp / 0.2e1_dp * t852 * t501 + t1224 * t369 + 0.3e1_dp * t1264  &
      * t1277 - 0.3e1_dp / 0.4e1_dp * t84 * t811 * t350 * t281 - t803 * t359  &
      * t350 * t281 / 0.2e1_dp - 0.3e1_dp / 0.2e1_dp * t863 * t1288 - (2._dp * t367  &
      * t401 * t1218) + (t367 * t368 * t726) + 0.5e1_dp / 0.2e1_dp * t920  &
      * t487 + t1205 * t342
      t1319 = -0.75e2_dp / 0.4e1_dp * t1082 * t422 * t789 * t350 * t281 - t1157  &
      * t156 - 0.15e2_dp / 0.2e1_dp * t1065 * t422 * t1066 * t281 + t1157 - (2._dp  &
      * t303 * t459) + 0.5e1_dp / 0.3e1_dp * t476 * t389 * t428 - 0.25e2_dp  &
      / 0.6e1_dp * t1026 * t219 * r3 * t5 * t115 * t281 + (3._dp * t1174 *  &
      t417) - t526 * t412 - (2._dp * t91 * t305 * t726) + (t1231 + t1299)  &
      * omega * t134 + 0.6e1_dp * t400 * t334 * t281 - 0.5e1_dp * t1017 * t427  &
      * t6 * t281 + 0.10e2_dp * t421 * t422 * t127 * t350 * t281 - (2._dp  &
      * t457 * t306) + 0.5e1_dp / 0.2e1_dp * t1071 * t533
      t1324 = d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
      t1336 = t160 * t281
      t1351 = t476 * t420
      t1382 = t136 * t475
      t1397 = t136 * t137
      t1405 = -t136 * t1155 * t93 * t118 + (2._dp * t956 * t510) + 0.5e1_dp  &
      / 0.2e1_dp * (t956) * (t513) + (2._dp * t1382 * t390) - (6._dp  &
      * t388 * t961 * t1218) - (5._dp * t966 * t1181) + (2._dp * t388  &
      * t389 * t726) + 0.5e1_dp / 0.2e1_dp * (t1382) * (t394) - (5._dp  &
      * t966 * t1261) - 0.75e2_dp / 0.4e1_dp * t1397 * t1210 + 0.10e2_dp * t1397  &
      * t1270 + 0.5e1_dp / 0.2e1_dp * (t388) * (t393) * (t1198)
      t1413 = t140 * t456
      t1435 = -t140 * t1255 * t93 * t129 + (2._dp * t985 * t519) + 0.3e1_dp  &
      / 0.2e1_dp * (t985) * (t522) + (2._dp * t1413 * t402) - (6._dp  &
      * t400 * t990 * t1218) - (3._dp * t995 * t1288) + (2._dp * t400  &
      * t401 * t726) + 0.3e1_dp / 0.2e1_dp * (t1413) * (t406) - (3._dp  &
      * t995 * t1190) - 0.27e2_dp / 0.4e1_dp * t95 * t1266 + 0.3e1_dp * t95 * t1277  &
      + 0.3e1_dp / 0.2e1_dp * (t400) * (t405) * (t1184)
      t1443 = -0.15e2_dp / 0.2e1_dp * t1065 * t422 * t351 * t281 + t158 * (t1324  &
      + (t726 * t161 - 0.2e1_dp * t246 * t434 * t281 + 0.2e1_dp * t1040 * t246  &
      * t281 - t435 * t726) * t438 * t160 - t1049 * t1336 + t439 * t281)  &
      - 0.12e2_dp * t106 * t782 * t118 * t422 * t1218 + 0.5e1_dp / 0.2e1_dp *  &
      t421 * t422 * t1198 - (3._dp * t330 * t478) + 0.5e1_dp / 0.2e1_dp * t1351  &
      * t423 + 0.12e2_dp * t106 * t782 * t286 * t281 - 0.3e1_dp * t106 * t333  &
      * t726 - t81 * t93 * t726 + f12 * t1255 * t94 + (3._dp * t1089  &
      * t530) - (3._dp * t476 * t334) + 0.2e1_dp * t81 * t108 * t286 * t281  &
      - t505 * t382 / 0.3e1_dp + (t1405 + t1435) * t145 * t151 + 0.3e1_dp * t415  &
      * t155 * t44 * t726
      e_ndrho_rho = e_ndrho_rho + ( -0.4e1_dp / 0.3e1_dp * t283 * t545 - t80 * (t1319 + t1443) * Clda ) * sx
      t1447 = t281 ** 2
      t1448 = t1447 * t115
      t1452 = d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
      t1481 = 2._dp * t728 * t103 + 4._dp * t249 * t473 + t11 * t15 * (2._dp * g2 * t3  &
      * t19 + 12._dp * g3 * t1 * t24 * t32)
      t1483 = t1481 * E * t108
      t1500 = 2._dp * t728 * t88 + 4._dp * t249 * t454 + t11 * t15 * (F2 * t744 *  &
      t68 - 2._dp * t450 * t278 + 2._dp * t85 * t752 - t85 * t769)
      t1529 = t789 * t1447
      t1533 = t127 * t1447
      t1537 = t115 * t772
      t1552 = t1447 * t114
      t1562 = t113 * t1447
      t1566 = t114 * t772
      t1570 = -t136 * t1481 * t93 * t118 + (4._dp * t1382 * t510) + (5._dp  &
      * t1382 * t513) - (6._dp * t388 * t961 * t1447) - (10._dp * t388  &
      * t420 * t1448) + (2._dp * t388 * t389 * t772) - 0.75e2_dp / 0.4e1_dp *  &
      (t388) * (t972) * (t1529) + (10._dp * t388 * t393 * t1533)  &
      + 0.5e1_dp / 0.2e1_dp * (t388) * (t393) * (t1537) - t140  &
      * t1500 * t93 * t129 + (4._dp * t1413 * t519) + (3._dp * t1413  &
      * t522) - (6._dp * t400 * t990 * t1447) - (6._dp * t400 * t108 * t372  &
      * t1552) + (2._dp * t400 * t401 * t772) - 0.27e2_dp / 0.4e1_dp * (t400)  &
      * (t1001) * (t1448) + (3._dp * t400 * t405 * t1562)  &
      + 0.3e1_dp / 0.2e1_dp * (t400) * (t405) * (t1566)
      t1576 = -15._dp * t1065 * t422 * t1448 + t158 * (t1452 + (t772 * t161 -  &
      2._dp * t1447 * t434 + 2._dp * t1040 * t1447 - t435 * t772) * t438 * t160  &
      - t539 * t1048 * t1336 + t540 * t281) + t1483 - t81 * t93 * t772  &
      + f12 * t1500 * t94 + 2._dp * t81 * t108 * t1447 - 6._dp * t476 * t478 - 3._dp  &
      * t106 * t333 * t772 - 4._dp * t457 * t459 + t1570 * t145 * t151 + 10._dp  &
      * t421 * t422 * t1533
      t1618 = -(2._dp * t81 * t799 * t1447) - (t81 * t93 * t359 * t1447)  &
      + (t81 * t355 * t772) - 0.3e1_dp / 0.4e1_dp * (t81) * (t812)  &
      * (t1447) + (t81 * t360 * t772) / 0.2e1_dp - t125 * t1500  &
      * t83 * t129 + (2._dp * t1224 * t498) + (3._dp * t1224 * t501) -  &
      (2._dp * t367 * t401 * t1447) + (t367 * t368 * t772) - (3._dp  &
      * t367 * t405 * t1552)
      t1652 = (3._dp * t367 * t373 * t1562) + 0.3e1_dp / 0.2e1_dp * (t367)  &
      * (t373) * (t1566) - 0.27e2_dp / 0.4e1_dp * (t367) * (t879)  &
      * (t1448) - t110 * t1481 * t83 * t118 + (2._dp * t1205 * t484)  &
      + (5._dp * t1205 * t487) - (2._dp * t340 * t389 * t1447) + (t340  &
      * t341 * t772) - (5._dp * t340 * t393 * t1448) + (10._dp *  &
      t340 * t346 * t1533) + 0.5e1_dp / 0.2e1_dp * (t340) * (t346) *  &
      (t1537) - 0.75e2_dp / 0.4e1_dp * (t340) * (t793) * (t1529)
      t1672 = 0.5e1_dp / 0.2e1_dp * t421 * t422 * t1537 - 0.75e2_dp / 0.4e1_dp * t1082  &
      * t422 * t1529 + (6._dp * t91 * t1111 * t1447) - (2._dp * t91 *  &
      t305 * t772) + (t1618 + t1652) * omega * t134 + (12._dp * t106 * t782  &
      * t1447) - t1483 * t156 + (6._dp * t1174 * t530) + (5._dp * t1351  &
      * t533) - (12._dp * t783 * t155 * t44 * t1447) + (3._dp * t415  &
      * t155 * t44 * t772)
      e_ndrho_ndrho = e_ndrho_ndrho + ( -t80 * (t1576 + t1672) * Clda ) * sx
    END IF

  END SUBROUTINE xwpbe_lda_calc_4

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xwpbe_lsd_info ( reference, shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    IF ( PRESENT ( reference ) ) THEN
      reference = "Jochen Heyd and Gustavo E. Scuseria, J. Chem. Phys., 120, 7274 {LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "shortrange part of PBE exchange {LSD}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho_spin=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=2
  END SUBROUTINE xwpbe_lsd_info

! *****************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lsd
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evalated,
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
! *****************************************************************************
  SUBROUTINE xwpbe_lsd_eval ( rho_set, deriv_set, order, xwpbe_params)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: xwpbe_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xwpbe_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, omega, sx, sx0
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrhoa, &
      e_ndrhoa_ndrhoa, e_ndrhoa_rhoa, e_ndrhob, e_ndrhob_ndrhob, &
      e_ndrhob_rhob, e_rhoa, e_rhoa_rhoa, e_rhob, e_rhob_rhob, norm_drhoa, &
      norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(bo)

    CALL cite_reference(Heyd2004)

    CALL section_vals_val_get(xwpbe_params,"SCALE_X",r_val=sx)
    CALL section_vals_val_get(xwpbe_params,"SCALE_X0",r_val=sx0)
    CALL section_vals_val_get(xwpbe_params,"OMEGA",r_val=omega)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)

    CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,norm_drhoa=norm_drhoa,&
         norm_drhob=norm_drhob,local_bounds=bo,rho_cutoff=epsilon_rho,&
         drho_cutoff=epsilon_norm_drho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPASSERT(stat==0)
    ELSE
      dummy=> rhoa
    END IF

    e_0 => dummy
    e_rhoa => dummy
    e_rhob => dummy
    e_ndrhoa => dummy
    e_ndrhob => dummy
    e_rhoa_rhoa => dummy
    e_rhob_rhob => dummy
    e_ndrhoa_rhoa => dummy
    e_ndrhob_rhob => dummy
    e_ndrhoa_ndrhoa => dummy
    e_ndrhob_ndrhob => dummy

    IF (order>=0) THEN
      deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (order>=1.OR.order==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rhoa)
      deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rhob)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrhob)
    END IF
   IF (order>=2.OR.order==-2) THEN
     deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
           allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa)
     deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
           allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
           allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhoa)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
           allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhob)
     deriv => xc_dset_get_derivative(deriv_set,&
           "(norm_drhoa)(norm_drha)", allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhoa)
     deriv => xc_dset_get_derivative(deriv_set,&
           "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_ndrhob)
   END IF
    IF (order>2.OR.order<-2) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 2 not implemented", &
            error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp shared(npoints, order, rhoa, norm_drhoa, e_0, e_rhoa, e_ndrhoa) &
    !$omp shared(e_rhoa_rhoa, e_ndrhoa_rhoa, e_ndrhoa_ndrhoa, epsilon_rho) &
    !$omp shared(sx, sx0, omega) &
    !$omp shared(rhob, norm_drhob, e_rhob, e_ndrhob, e_rhob_rhob) &
    !$omp shared(e_ndrhob_rhob, e_ndrhob_ndrhob)

    !Call lsd_calc for alpha - and beta - spins

    CALL xwpbe_lsd_calc(npoints,order, rho=rhoa, norm_drho=norm_drhoa, &
         e_0=e_0,e_rho=e_rhoa,e_ndrho=e_ndrhoa,e_rho_rho=e_rhoa_rhoa, &
         e_ndrho_rho=e_ndrhoa_rhoa, e_ndrho_ndrho=e_ndrhoa_ndrhoa, &
         epsilon_rho=epsilon_rho, &
         sx=sx,sx0=sx0,omega=omega)
    CALL xwpbe_lsd_calc(npoints,order, rho=rhob, norm_drho=norm_drhob, &
         e_0=e_0,e_rho=e_rhob,e_ndrho=e_ndrhob,e_rho_rho=e_rhob_rhob, &
         e_ndrho_rho=e_ndrhob_rhob, e_ndrho_ndrho=e_ndrhob_ndrhob, &
         epsilon_rho=epsilon_rho,&
         sx=sx,sx0=sx0,omega=omega)

    !$omp end parallel

    IF (cp_debug) THEN
      DEALLOCATE(dummy,stat=stat)
      CPASSERT(stat==0)
    ELSE
      NULLIFY(dummy)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xwpbe_lsd_eval

! *****************************************************************************
!> \brief evaluates the screened hole averaged PBE exchange functional for lsd.
!> \param npoints ...
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param rho , ndrho: density and norm of the density gradient
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
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
! *****************************************************************************
  SUBROUTINE xwpbe_lsd_calc(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho, &
                            e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, &
                            epsilon_rho,sx, sx0, omega)

    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: rho, norm_drho, e_0, e_rho, &
                                                e_ndrho, e_rho_rho, &
                                                e_ndrho_rho, e_ndrho_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, sx0, omega

    CHARACTER(len=*), PARAMETER :: routineN = 'xwpbe_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: e_0_temp, my_ndrho, my_rho
    REAL(KIND=dp)                            :: ss, ss2, sscale, t1, t2, t3, &
                                                t4, t5, t6, t7, t8, ww

    !$omp do

    DO ip =1,npoints
      !According to spin-scaling relation, we need twice the density and its gradient
      my_rho = 2.0_dp*MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho = 2.0_dp*MAX(norm_drho(ip),0.0_dp)

        !Do some precalculation in order to catch the correct branch afterwards
        sscale = 1.0_dp
        t1 = pi ** 2
        t2 = t1 * my_rho
        t3 = t2 ** (0.1e1_dp / 0.3e1_dp)
        t4 = 0.1e1_dp / t3
        t5 = omega * t4
        ww = 0.6933612743506347048433524e0_dp * t5
        t6 = my_ndrho * t4
        t7 = 0.1e1_dp / my_rho
        t8 = t7 * sscale
        ss = 0.3466806371753173524216762e0_dp * t6 * t8
        IF( ss > scutoff) THEN
          ss2 = ss*ss
          sscale = ((smax)*ss2-(sconst))/(ss2*ss)
        END IF
        e_0_temp = 0.0_dp
        IF(sx0/=0.0_dp) THEN
          !original PBE hole
           IF(ss*sscale>gcutoff) THEN
             CALL xwpbe_lda_calc_0(e_0_temp, e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                   e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                   my_ndrho, sscale, sx0, order)
           ELSE
             CALL xwpbe_lda_calc_01(e_0_temp, e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                   e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                   my_ndrho, sscale, sx0, order)
           END IF
           !According to spin-scaling relation, we need only half of the energy
           e_0(ip) = e_0(ip) + 0.5_dp * e_0_temp
        END IF
        e_0_temp = 0.0_dp
        IF(sx/=0.0_dp) THEN
          IF(ww<wcutoff .AND. ss*sscale>gcutoff) THEN
            CALL xwpbe_lda_calc_1(e_0_temp, e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          ELSE IF(ww<wcutoff .AND. ss*sscale<=gcutoff) THEN
            CALL xwpbe_lda_calc_2(e_0_temp, e_rho(ip), e_ndrho(ip), e_rho_rho(ip), &
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip),my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          ELSE IF(ww>=wcutoff .AND. ss*sscale>gcutoff) THEN
            CALL xwpbe_lda_calc_3(e_0_temp, e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          ELSE
            CALL xwpbe_lda_calc_4(e_0_temp, e_rho(ip), e_ndrho(ip), e_rho_rho(ip),&
                                  e_ndrho_rho(ip), e_ndrho_ndrho(ip), my_rho,&
                                  my_ndrho, omega, sscale, sx, order)
          END IF
         !According to spin-scaling relation, we need only half of the energy
        END IF
        e_0(ip) = e_0(ip) + 0.5_dp * e_0_temp
      END IF
    END DO

    !$omp end do

  END SUBROUTINE xwpbe_lsd_calc

! *****************************************************************************
!> \brief These functions evaluate products exp(x)*Ei(x) and pi*exp(x)*erfc(sqrt(x)),
!>      as well as their derivatives with respect to various combinations of
!>      rho and norm_drho.
!> \param Q , dQrho, dQndrho, d2Qrhondrho.... :
!>        Argument Q and derivatives with respect to various combinations of
!>        rho and norm_drho
!> \retval exei ...
!> \par History
!>      05.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>     - In order to avoid numerical instabilities, these routines use Taylor-
!>       expansions for the above core-products for large arguments.
!>     - When adapting this module for higher order derivatives, appropriate
!>       functions have to be provided!
! *****************************************************************************
  FUNCTION exei(Q)
    REAL(dp), INTENT(IN)                     :: Q
    REAL(dp)                                 :: exei

    exei = 0.0_dp
    IF( Q<expcutoff ) THEN
      !Use exact product
      exei = EXP(Q)*expint(1,Q)
    ELSE
      !Use approximation
      exei = (1._dp/Q)*(Q*Q+exei1*Q+exei2)/(Q*Q+exei3*Q+exei4)
    END IF
  END FUNCTION exei

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \retval exer ...
! *****************************************************************************
  FUNCTION exer(Q)
    REAL(dp), INTENT(IN)                     :: Q
    REAL(dp)                                 :: exer

    REAL(dp)                                 :: Q3, Q5

    exer = 0.0_dp
    IF( Q<expcutoff ) THEN
      !Use exact expression
      exer = pi*EXP(Q)*erfc(SQRT(Q))
    ELSE
      !Use approximation
      Q3 = Q*Q*Q
      Q5 = Q3*Q*Q
      exer = pi* (1.0_dp/SQRT(Q*pi) - 1.0_dp/(2.0_dp*SQRT(pi*Q3)) + 3.0_dp/(4.0_dp*(SQRT(pi*Q5))))
    END IF
  END FUNCTION exer

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \retval dexeirho ...
! *****************************************************************************
  FUNCTION dexeirho(Q,dQrho)
    REAL(dp), INTENT(IN)                     :: Q, dQrho
    REAL(dp)                                 :: dexeirho

    dexeirho = dQrho*(exei(Q)-1.0_dp/Q)
  END FUNCTION dexeirho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \retval dexeindrho ...
! *****************************************************************************
  FUNCTION dexeindrho(Q,dQndrho)
    REAL(dp), INTENT(IN)                     :: Q, dQndrho
    REAL(dp)                                 :: dexeindrho

    dexeindrho = dQndrho*(exei(Q)-1.0_dp/Q)
  END FUNCTION dexeindrho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \retval dexerrho ...
! *****************************************************************************
  FUNCTION dexerrho(Q,dQrho)
    REAL(dp), INTENT(IN)                     :: Q, dQrho
    REAL(dp)                                 :: dexerrho

    dexerrho = dQrho*exer(Q)-dQrho*rootpi/SQRT(Q)
  END FUNCTION dexerrho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \retval dexerndrho ...
! *****************************************************************************
  FUNCTION dexerndrho(Q,dQndrho)
    REAL(dp), INTENT(IN)                     :: Q, dQndrho
    REAL(dp)                                 :: dexerndrho

    dexerndrho = dQndrho*exer(Q)-dQndrho*rootpi/SQRT(Q)
  END FUNCTION dexerndrho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param d2Qrhorho ...
!> \retval d2exeirhorho ...
! *****************************************************************************
  FUNCTION d2exeirhorho(Q,dQrho,d2Qrhorho)
    REAL(dp), INTENT(IN)                     :: Q, dQrho, d2Qrhorho
    REAL(dp)                                 :: d2exeirhorho

   d2exeirhorho = exei(Q)*(d2Qrhorho+dQrho*dQrho) + &
                  1.0_dp/(Q*Q)*(-Q*dQrho*dQrho-Q*d2Qrhorho+dQrho*dQrho)
  END FUNCTION d2exeirhorho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param d2Qrhorho ...
!> \retval d2exerrhorho ...
! *****************************************************************************
  FUNCTION d2exerrhorho(Q,dQrho,d2Qrhorho)
    REAL(dp), INTENT(IN)                     :: Q, dQrho, d2Qrhorho
    REAL(dp)                                 :: d2exerrhorho

    REAL(dp)                                 :: pi12, Q12

    Q12 = SQRT(Q)
    pi12 = rootpi

    d2exerrhorho = exer(Q)*(d2Qrhorho+dQrho*dQrho) - dQrho*dQrho/(pi12*Q12) + &
                   0.5_dp*dQrho*dQrho/(pi12*Q*Q12)-d2Qrhorho/(pi12*Q12)
  END FUNCTION d2exerrhorho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param dQndrho ...
!> \param d2Qrhondrho ...
!> \retval d2exeirhondrho ...
! *****************************************************************************
  FUNCTION d2exeirhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
    REAL(dp), INTENT(IN)                     :: Q, dQrho, dQndrho, d2Qrhondrho
    REAL(dp)                                 :: d2exeirhondrho

    d2exeirhondrho = exei(Q)*(d2Qrhondrho+dQrho*dQndrho) - &
                     1.0_dp/Q*(dQrho*dQndrho+d2Qrhondrho) + 1.0_dp/(Q*Q)*dQrho*dQndrho
  END FUNCTION d2exeirhondrho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQrho ...
!> \param dQndrho ...
!> \param d2Qrhondrho ...
!> \retval d2exerrhondrho ...
! *****************************************************************************
  FUNCTION d2exerrhondrho(Q,dQrho,dQndrho,d2Qrhondrho)
    REAL(dp), INTENT(IN)                     :: Q, dQrho, dQndrho, d2Qrhondrho
    REAL(dp)                                 :: d2exerrhondrho

    REAL(dp)                                 :: pi12, Q12

    Q12 = SQRT(Q)
    pi12 = rootpi

    d2exerrhondrho = exer(Q)*(d2Qrhondrho+dQrho*dQndrho) - 1.0_dp/(pi12*Q12)*dQrho*dQndrho &
                     + 0.5_dp/(pi12*Q12*Q)*dQrho*dQndrho - 1.0_dp/(pi12*Q12)*d2Qrhondrho
  END FUNCTION d2exerrhondrho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \param d2Qndrhondrho ...
!> \retval d2exeindrhondrho ...
! *****************************************************************************
  FUNCTION d2exeindrhondrho(Q,dQndrho,d2Qndrhondrho)
    REAL(dp), INTENT(IN)                     :: Q, dQndrho, d2Qndrhondrho
    REAL(dp)                                 :: d2exeindrhondrho

    d2exeindrhondrho = exei(Q)*(d2Qndrhondrho+dQndrho*dQndrho) + &
                       1.0_dp/(Q*Q)*(-Q*dQndrho*dQndrho-Q*d2Qndrhondrho+dQndrho*dQndrho)
  END FUNCTION d2exeindrhondrho

! *****************************************************************************
!> \brief ...
!> \param Q ...
!> \param dQndrho ...
!> \param d2Qndrhondrho ...
!> \retval d2exerndrhondrho ...
! *****************************************************************************
  FUNCTION d2exerndrhondrho(Q,dQndrho,d2Qndrhondrho)
    REAL(dp), INTENT(IN)                     :: Q, dQndrho, d2Qndrhondrho
    REAL(dp)                                 :: d2exerndrhondrho

    REAL(dp)                                 :: pi12, Q12

    Q12 = SQRT(Q)
    pi12 = rootpi

    d2exerndrhondrho = exer(Q)*(d2Qndrhondrho+dQndrho*dQndrho) - dQndrho*dQndrho/(pi12*Q12) &
                       + 0.5_dp*dQndrho*dQndrho/(pi12*Q*Q12)-d2Qndrhondrho/(pi12*Q12)
  END FUNCTION d2exerndrhondrho

END MODULE xc_xwpbe

