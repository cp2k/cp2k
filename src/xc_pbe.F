!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/xc_pbe [1.0] *
!!
!!   NAME
!!     xc_pbe
!!
!!   FUNCTION
!!     calculates the pbe correlation functional
!!
!!   NOTES
!!     This was generated with the help of a maple worksheet that you can
!!     find under doc/pbe.mw .
!!     I did not add 3. derivatives in the polazied (lsd) case because it
!!     would have added 2500 lines of code. If you need them the worksheet
!!     is already prepared for them, and by uncommenting a couple of lines
!!     you should be able to generate them.
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_pbe
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
  USE input_section_types, only: section_vals_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_pbe'
  REAL(kind=dp), PARAMETER, PRIVATE :: a=0.04918_dp, b=0.132_dp,&
       c=0.2533_dp,d=0.349_dp

  PUBLIC :: pbe_lda_info, pbe_lsd_info, pbe_lda_eval, pbe_lsd_eval
!!***
!****************************************************************************
CONTAINS

!!****f* xc_pbe/pbe_lda_info [1.0] *
!!
!!   NAME
!!     pbe_lda_info
!!
!!   SYNOPSIS
!!     Subroutine pbe_lda_info(reference, shortform, needs, max_deriv,&
!!         error)
!!       Character(Len=*), Intent (OUT), Optional:: reference, shortform
!!       Type(xc_rho_cflags_type), Intent (INOUT), Optional:: needs
!!       Integer, Intent (OUT), Optional:: max_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pbe_lda_info
!!
!!   FUNCTION
!!     return various information on the functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pbe_params: section selecting the varius parameters for the functional
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE pbe_lda_info(pbe_params,reference,shortform, needs, max_deriv,&
       error)
    TYPE(section_vals_type), POINTER :: pbe_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
         INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lda_info', &
         routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "J.P.Perdew, K.Burke, M.Ernzerhof, "//&
            "Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868, (1996)"//&
            "{LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Perdew-Burke-Ernzerhof exchange-correlation energy functional (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%rho_1_3=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE pbe_lda_info
!***************************************************************************

!!****f* xc_pbe/pbe_lsd_info [1.0] *
!!
!!   NAME
!!     pbe_lsd_info
!!
!!   SYNOPSIS
!!     Subroutine pbe_lsd_info(reference, shortform, needs, max_deriv,&
!!         error)
!!       Character(Len=*), Intent (OUT), Optional:: reference, shortform
!!       Type(xc_rho_cflags_type), Intent (INOUT), Optional:: needs
!!       Integer, Intent (OUT), Optional:: max_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pbe_lsd_info
!!
!!   FUNCTION
!!     return various information on the functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pbe_params: section selecting the varius parameters for the functional
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pbe_lsd_info(pbe_params,reference,shortform, needs, max_deriv,&
     error)
    TYPE(section_vals_type), POINTER :: pbe_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_info', &
      routineP = moduleN//':'//routineN

  IF ( PRESENT ( reference ) ) THEN
     reference = "J.P.Perdew, K.Burke, M.Ernzerhof, "//&
          "Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868, (1996)"//&
          "{LSD version}"
  END IF
  IF ( PRESENT ( shortform ) ) THEN
     shortform = "Perdew-Burke-Ernzerhof exchange-correlation energy functional (LSD)"
  END IF
  IF (PRESENT(needs)) THEN
     needs%rho_spin=.TRUE.
     needs%norm_drho_spin=.TRUE.
     needs%norm_drho=.TRUE.
  END IF
  IF (PRESENT(max_deriv)) max_deriv=2

END SUBROUTINE pbe_lsd_info
!***************************************************************************

!!****f* xc_pbe/pbe_lda_eval [1.0] *
!!
!!   NAME
!!     pbe_lda_eval
!!
!!   SYNOPSIS
!!     Subroutine pbe_lda_eval(rho_set, deriv_set, grad_deriv, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: grad_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pbe_lda_eval
!!
!!   FUNCTION
!!     evaluates the pbe correlation functional for lda
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
SUBROUTINE pbe_lda_eval(rho_set,deriv_set,grad_deriv,pbe_params,error)
    TYPE(section_vals_type), POINTER :: pbe_params
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, e_ndrho_rho_rho, e_rho, &
      e_rho_rho, e_rho_rho_rho, norm_drho, rho, rho_1_3, e_ndrho_ndrho_ndrho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(bo)

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
          norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
          drho_cutoff=epsilon_norm_drho,error=error)
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
     e_rho_rho => dummy
     e_ndrho_rho => dummy
     e_ndrho_ndrho => dummy
     e_rho_rho_rho => dummy
     e_ndrho_rho_rho => dummy
     e_ndrho_ndrho_rho => dummy
     e_ndrho_ndrho_ndrho => dummy

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
     END IF
     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
     END IF
     IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(rho)(rho)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)(rho)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,&
             error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)
     END IF
     IF (grad_deriv>3.OR.grad_deriv<-3) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 3 not implemented", &
             error=error, error_level=cp_failure_level)
     END IF

     CALL pbe_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
          e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
          e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
          e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
          e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,&
          grad_deriv=grad_deriv,&
          npoints=npoints,epsilon_rho=epsilon_rho,epsilon_norm_drho=epsilon_norm_drho,&
          pbe_params=pbe_params,error=error)

     IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     ELSE
        NULLIFY(dummy)
     END IF
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE pbe_lda_eval
!***************************************************************************

!!****f* xc_pbe/pbe_lda_calc [1.0] *
!!
!!   NAME
!!     pbe_lda_eval
!!
!!   SYNOPSIS
!!     Subroutine pbe_lda_eval(rho_set, deriv_set, grad_deriv, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: grad_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pbe_lda_eval
!!
!!   FUNCTION
!!     evaluates the pbe correlation functional for lda
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
SUBROUTINE pbe_lda_calc(rho, rho_1_3, norm_drho,&
     e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
     e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
     e_ndrho_ndrho_ndrho,grad_deriv,npoints,epsilon_rho, epsilon_norm_drho,&
     pbe_params,error)
  INTEGER, INTENT(in)                      :: npoints, grad_deriv
  REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: e_ndrho_ndrho_rho, &
       e_ndrho_rho_rho, e_rho_rho_rho, e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, &
       e_ndrho, e_rho, e_0, e_ndrho_ndrho_ndrho
  REAL(kind=dp), DIMENSION(1:npoints), &
       INTENT(in)                             :: norm_drho, rho_1_3, rho
  REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_norm_drho
  TYPE(section_vals_type), POINTER :: pbe_params
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lda_calc', &
       routineP = moduleN//':'//routineN

  INTEGER                                  :: ii
  real(kind=dp) :: my_rho, my_norm_drho, t1, t2, t3, t4, t5, t6, t7, t8, rs,&
       t11, t12, t15, t17, t19, t22, t23, e_c_u_0, t36, t49, t51, t55, &
       t57, t58, gamma_var, t59, t60, k_f, t61, t62, k_s, t63, t64, t, t66,&
       t68, t69, A, t72, t73, t74, t75, t76, t77, t79, t80, t81, t84, t85,&
       epsilon_cGGA, mu, kf, ex_unif, t88, t89, s, t91, t94, Fx, t97,&
       t98, t100, t101, t102, rsrho, t108, t109, t110, t111, t112, t117, &
       t120, t121, t122, e_c_u_0rho, t125, k_frho, t129, k_srho, t131, &
       t132, t133, t134, t136, trho, t138, t140, t141, t142, Arho, t145, &
       t146, t149, t150, t152, t153, t154, t157, t158, t159, t160, t163, &
       t164, t167, t168, t171, t172, t173, epsilon_cGGArho, kfrho, &
       ex_unifrho, t177, t178, t179, t180, t182, srho, t184, t186, Fxrho, &
       t190, tnorm_drho, t194, t197, t206, t210, t211, Hnorm_drho,&
       snorm_drho, Fxnorm_drho, t217, t218, t219, t224
  real(kind=dp) :: rsrhorho, t229, t230, t234, t235, t236, t237, t240, t241, t242, &
       t244, t245, t246, t250, t253, t254, t255, t258, t259, t260, t261, &
       t262, t263, e_c_u_0rhorho, e_c_u_01rho, t269, k_frhorho, t273, t274,&
       t279, k_srhorho, k_s1rho, t281, t282, t285, t287, t288, t291, t293,&
       t294, trhorho, t295, t1rho, t299, t301, t302, t303, t304, t311, &
       Arhorho, A1rho, t318, t321, t323, t324, t325, t326, t329, t330, &
       t335, t336, t339, t342, t345, t346, t348, t349, t351, t352, t354, &
       t356, t357, t358, t361, t365, t366, t369, t373, t375, t376, t377, &
       t380, t383, t384, t389, t392, t393, t398, t402, t403, t405, t406, &
       t407, t415, t416, t418, epsilon_cGGArhorho, kfrhorho, ex_unifrhorho,&
       ex_unif1rho, t423, t424, t425, t426, t428, t429, t430, t432, t433, &
       srhorho, s1rho, t436, t437, t438, t442, Fxrhorho, Fx1rho, t453
  real(kind=dp) :: t457, t460, t462, tnorm_drhorho, t464, t467, t470, &
       t471, t474, t477, t478, t479, t482, t486, t490, t493, t494, t500, &
       t501, t506, t513, t518, t522, t523, t525, Hnorm_drhorho, t528, &
       snorm_drhorho, t532, t536, Fxnorm_drhorho, t545, &
       t546, t549, t553, t559, t566, t571, t573, t574, Hnorm_drhonorm_drho,&
       t576, Fxnorm_drhonorm_drho, t586, &
       rsrhorhorho, rs2rho, t603, t611, t612, t618, t628, t629, t634, t652,&
       t658, t664, t686, e_c_u_0rhorhorho, e_c_u_0rho1rho, e_c_u_01rhorho,&
       e_c_u_02rho, k_frhorhorho, k_f2rho, t730, t738, t741, k_srho1rho
  real(kind=dp) :: k_s1rhorho, k_s2rho, t750, t754, t760, trhorhorho, t797, trho1rho, &
       t1rhorho, t2rho, t806, t807, t808, t824, t852, Arhorhorho, Arho1rho,&
       A1rhorho, A2rho, t874, t876, t877, t878, t879, t882, t890, t891, &
       t898, t900, t903, t905, t906, t908, t909, t911, t913, t914, t915, &
       t921, t927, t930, t943, t951, t956, t958, t960, t962, t964, t967, &
       t972, t980, t984, t992, t996, t1004, t1006, t1019, t1021, t1049, &
       t1050, t1053, t1054, t1058, t1061, t1068, t1071, t1074, t1076, &
       t1079, t1080, t1081, t1083, t1086, t1087, t1089, t1090, t1091, &
       t1094, t1096, t1098, t1100, t1102, t1105, t1111, t1112, t1130, &
       t1133, t1159, t1166, t1167, t1172, t1176, t1180, t1187, t1194, &
       t1207, t1208, t1236, t1237, t1240, t1273, kfrhorhorho, kf2rho
  real(kind=dp) :: ex_unifrho1rho, ex_unif1rhorho, ex_unif2rho, t1286, srho1rho, &
       s1rhorho, s2rho, t1324, t1327, t1328, t1329, t1330, t1333, &
       Fxrho1rho, Fx1rhorho, Fx2rho, ex_ldarhorhorho, &
       t1406, tnorm_drhorhorho, tnorm_drho1rho, t1425, t1433, t1438, t1439,&
       t1442, t1443, t1446, t1451, t1455, t1458, t1459, t1477, t1480, &
       t1484, t1497, t1501, t1506, t1514, t1517, t1520, t1533, t1552, &
       t1563, t1564, t1571, t1578, t1584, t1588, t1627, t1636, t1648, &
       t1683, t1684, t1686, snorm_drho1rho, t1702, t1720, Fxnorm_drho1rho, &
       t1747, t1759, t1766, t1786, t1794, t1806
  real(kind=dp) :: t1811, t1815, t1829, t1833, t1837, t1842, t1848, t1852, &
       t2norm_drho, t1865, t1871, t1876, &
       t1883, t1885, t1891, t1910, t1913, t1924, t1932, t1943, t1951, &
       t1956, s2norm_drho

  t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
  t2 = 4 ** (0.1e1_dp / 0.3e1_dp)
  t3 = t2 ** 2
  t4 = t1 * t3
  t5 = 0.1e1_dp / 0.3141592654e1_dp
  SELECT CASE(grad_deriv)
  CASE (0)
     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_norm_drho = norm_drho(ii)

           t6 = 0.1e1_dp / my_rho
           t7 = t5 * t6
           t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
           rs = t4 * t8 / 0.4e1_dp
           t11 = 0.1e1_dp + 0.21370e0_dp * rs
           t12 = sqrt(rs)
           t15 = t12 * rs
           t17 = rs ** 0.20e1_dp
           t19 = 0.75957e1_dp * t12 + 0.35876e1_dp * rs + 0.16382e1_dp * t15 &
                + 0.49294e0_dp * t17
           t22 = 0.1e1_dp + 0.1608182432e2_dp / t19
           t23 = log(t22)
           e_c_u_0 = -0.62182e-1_dp * t11 * t23
           t36 = log(0.1e1_dp + 0.3216468318e2_dp / (0.141189e2_dp * t12 + &
                0.61977e1_dp * rs + 0.33662e1_dp * t15 + 0.62517e0_dp * t17))
           t49 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t12 + &
                0.36231e1_dp * rs + 0.88026e0_dp * t15 + 0.49671e0_dp * t17))
           t51 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t55 = log(0.2e1_dp)
           t57 = 0.3141592654e1_dp ** 2
           t58 = 0.1e1_dp / t57
           gamma_var = (0.1e1_dp - t55) * t58
           t59 = t57 * my_rho
           t60 = t59 ** (0.1e1_dp / 0.3e1_dp)
           k_f = t1 * t60
           t61 = k_f * t5
           t62 = sqrt(t61)
           k_s = 0.2e1_dp * t62
           t63 = 0.1e1_dp / k_s
           t64 = my_norm_drho * t63
           t = t64 * t6 / 0.2e1_dp
           t66 = 0.1e1_dp / gamma_var
           t68 = exp(-e_c_u_0 * t66)
           t69 = -0.1e1_dp + t68
           A = 0.66725e-1_dp * t66 / t69
           t72 = t ** 2
           t73 = t66 * t72
           t74 = A * t72
           t75 = 0.1e1_dp + t74
           t76 = A ** 2
           t77 = t72 ** 2
           t79 = 0.1e1_dp + t74 + t76 * t77
           t80 = 0.1e1_dp / t79
           t81 = t75 * t80
           t84 = 0.1e1_dp + 0.66725e-1_dp * t73 * t81
           t85 = log(t84)
           epsilon_cGGA = e_c_u_0 + gamma_var * t85
           mu = 0.2224166667e-1_dp * t57
           kf = k_f
           ex_unif = -0.3e1_dp / 0.4e1_dp * t5 * kf
           t88 = 0.1e1_dp / kf
           t89 = my_norm_drho * t88
           s = t89 * t6 / 0.2e1_dp
           t91 = s ** 2
           t94 = 0.1e1_dp + 0.1243781095e1_dp * mu * t91
           Fx = 0.1804e1_dp - 0.804e0_dp / t94
           t97 = my_rho * ex_unif

           e_0(ii) = e_0(ii)&
                + t97 * Fx + my_rho * epsilon_cGGA
        END IF
     END DO
  CASE(1)

     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_norm_drho = norm_drho(ii)

           t6 = 0.1e1_dp / my_rho
           t7 = t5 * t6
           t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
           rs = t4 * t8 / 0.4e1_dp
           t11 = 0.1e1_dp + 0.21370e0_dp * rs
           t12 = sqrt(rs)
           t15 = t12 * rs
           t17 = rs ** 0.20e1_dp
           t19 = 0.75957e1_dp * t12 + 0.35876e1_dp * rs + 0.16382e1_dp * t15 &
                + 0.49294e0_dp * t17
           t22 = 0.1e1_dp + 0.1608182432e2_dp / t19
           t23 = log(t22)
           e_c_u_0 = -0.62182e-1_dp * t11 * t23
           t36 = log(0.1e1_dp + 0.3216468318e2_dp / (0.141189e2_dp * t12 + &
                0.61977e1_dp * rs + 0.33662e1_dp * t15 + 0.62517e0_dp * t17))
           t49 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t12 + &
                0.36231e1_dp * rs + 0.88026e0_dp * t15 + 0.49671e0_dp * t17))
           t51 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t55 = log(0.2e1_dp)
           t57 = 0.3141592654e1_dp ** 2
           t58 = 0.1e1_dp / t57
           gamma_var = (0.1e1_dp - t55) * t58
           t59 = t57 * my_rho
           t60 = t59 ** (0.1e1_dp / 0.3e1_dp)
           k_f = t1 * t60
           t61 = k_f * t5
           t62 = sqrt(t61)
           k_s = 0.2e1_dp * t62
           t63 = 0.1e1_dp / k_s
           t64 = my_norm_drho * t63
           t = t64 * t6 / 0.2e1_dp
           t66 = 0.1e1_dp / gamma_var
           t68 = exp(-e_c_u_0 * t66)
           t69 = -0.1e1_dp + t68
           A = 0.66725e-1_dp * t66 / t69
           t72 = t ** 2
           t73 = t66 * t72
           t74 = A * t72
           t75 = 0.1e1_dp + t74
           t76 = A ** 2
           t77 = t72 ** 2
           t79 = 0.1e1_dp + t74 + t76 * t77
           t80 = 0.1e1_dp / t79
           t81 = t75 * t80
           t84 = 0.1e1_dp + 0.66725e-1_dp * t73 * t81
           t85 = log(t84)
           epsilon_cGGA = e_c_u_0 + gamma_var * t85
           mu = 0.2224166667e-1_dp * t57
           kf = k_f
           ex_unif = -0.3e1_dp / 0.4e1_dp * t5 * kf
           t88 = 0.1e1_dp / kf
           t89 = my_norm_drho * t88
           s = t89 * t6 / 0.2e1_dp
           t91 = s ** 2
           t94 = 0.1e1_dp + 0.1243781095e1_dp * mu * t91
           Fx = 0.1804e1_dp - 0.804e0_dp / t94
           t97 = my_rho * ex_unif

           e_0(ii) = e_0(ii)&
                + t97 * Fx + my_rho * epsilon_cGGA

           t98 = t8 ** 2
           t100 = 0.1e1_dp / t98 * t5
           t101 = my_rho ** 2
           t102 = 0.1e1_dp / t101
           rsrho = -t4 * t100 * t102 / 0.12e2_dp
           t108 = t19 ** 2
           t109 = 0.1e1_dp / t108
           t110 = t11 * t109
           t111 = 0.1e1_dp / t12
           t112 = t111 * rsrho
           t117 = rs ** 0.10e1_dp
           t120 = 0.3797850000e1_dp * t112 + 0.35876e1_dp * rsrho + &
                0.2457300000e1_dp * t12 * rsrho + 0.985880e0_dp * t117 * rsrho
           t121 = 0.1e1_dp / t22
           t122 = t120 * t121
           e_c_u_0rho = -0.1328829340e-1_dp * rsrho * t23 + 1.e0_dp&
                * t110 * t122
           t125 = t60 ** 2
           k_frho = t1 / t125 * t57 / 0.3e1_dp
           t129 = 0.1e1_dp / t62
           k_srho = t129 * k_frho * t5
           t131 = k_s ** 2
           t132 = 0.1e1_dp / t131
           t133 = my_norm_drho * t132
           t134 = t6 * k_srho
           t136 = t64 * t102
           trho = -t133 * t134 / 0.2e1_dp - t136 / 0.2e1_dp
           t138 = gamma_var ** 2
           t140 = t69 ** 2
           t141 = 0.1e1_dp / t140
           t142 = 0.1e1_dp / t138 * t141
           Arho = 0.66725e-1_dp * t142 * e_c_u_0rho * t68
           t145 = t66 * t
           t146 = t81 * trho
           t149 = Arho * t72
           t150 = A * t
           t152 = 0.2e1_dp * t150 * trho
           t153 = t149 + t152
           t154 = t153 * t80
           t157 = t79 ** 2
           t158 = 0.1e1_dp / t157
           t159 = t75 * t158
           t160 = A * t77
           t163 = t72 * t
           t164 = t76 * t163
           t167 = t149 + t152 + 0.2e1_dp * t160 * Arho + 0.4e1_dp * t164 * trho
           t168 = t159 * t167
           t171 = 0.133450e0_dp * t145 * t146 + 0.66725e-1_dp * t73 * t154 - &
                0.66725e-1_dp * t73 * t168
           t172 = gamma_var * t171
           t173 = 0.1e1_dp / t84
           epsilon_cGGArho = e_c_u_0rho + t172 * t173
           kfrho = k_frho
           ex_unifrho = -0.3e1_dp / 0.4e1_dp * t5 * kfrho
           t177 = kf ** 2
           t178 = 0.1e1_dp / t177
           t179 = my_norm_drho * t178
           t180 = t6 * kfrho
           t182 = t89 * t102
           srho = -t179 * t180 / 0.2e1_dp - t182 / 0.2e1_dp
           t184 = t94 ** 2
           t186 = 0.1e1_dp / t184 * mu
           Fxrho = 0.2e1_dp * t186 * s * srho
           t190 = my_rho * ex_unifrho

           e_rho(ii) = e_rho(ii)&
                +ex_unif * Fx + t190 * Fx + t97 * Fxrho + epsilon_cGGA + &
                my_rho * epsilon_cGGArho

           tnorm_drho = t63 * t6 / 0.2e1_dp
           t194 = t81 * tnorm_drho
           t197 = t66 * t163
           t206 = 0.2e1_dp * t150 * tnorm_drho + 0.4e1_dp * t164 * tnorm_drho
           t210 = 0.133450e0_dp * t145 * t194 + 0.133450e0_dp * t197 * A * &
                tnorm_drho * t80 - 0.66725e-1_dp * t73 * t159 * t206
           t211 = gamma_var * t210
           Hnorm_drho = t211 * t173
           snorm_drho = t88 * t6 / 0.2e1_dp
           Fxnorm_drho = 0.2e1_dp * t186 * s * snorm_drho


           e_ndrho(ii) = e_ndrho(ii)&
                +t97 * Fxnorm_drho + my_rho * Hnorm_drho
        END IF
     END DO
  CASE(2)
     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_norm_drho = norm_drho(ii)

           t6 = 0.1e1_dp / my_rho
           t7 = t5 * t6
           t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
           rs = t4 * t8 / 0.4e1_dp
           t11 = 0.1e1_dp + 0.21370e0_dp * rs
           t12 = sqrt(rs)
           t15 = t12 * rs
           t17 = rs ** 0.20e1_dp
           t19 = 0.75957e1_dp * t12 + 0.35876e1_dp * rs + 0.16382e1_dp * t15 &
                + 0.49294e0_dp * t17
           t22 = 0.1e1_dp + 0.1608182432e2_dp / t19
           t23 = log(t22)
           e_c_u_0 = -0.62182e-1_dp * t11 * t23
           t36 = log(0.1e1_dp + 0.3216468318e2_dp / (0.141189e2_dp * t12 + &
                0.61977e1_dp * rs + 0.33662e1_dp * t15 + 0.62517e0_dp * t17))
           t49 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t12 + &
                0.36231e1_dp * rs + 0.88026e0_dp * t15 + 0.49671e0_dp * t17))
           t51 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t55 = log(0.2e1_dp)
           t57 = 0.3141592654e1_dp ** 2
           t58 = 0.1e1_dp / t57
           gamma_var = (0.1e1_dp - t55) * t58
           t59 = t57 * my_rho
           t60 = t59 ** (0.1e1_dp / 0.3e1_dp)
           k_f = t1 * t60
           t61 = k_f * t5
           t62 = sqrt(t61)
           k_s = 0.2e1_dp * t62
           t63 = 0.1e1_dp / k_s
           t64 = my_norm_drho * t63
           t = t64 * t6 / 0.2e1_dp
           t66 = 0.1e1_dp / gamma_var
           t68 = exp(-e_c_u_0 * t66)
           t69 = -0.1e1_dp + t68
           A = 0.66725e-1_dp * t66 / t69
           t72 = t ** 2
           t73 = t66 * t72
           t74 = A * t72
           t75 = 0.1e1_dp + t74
           t76 = A ** 2
           t77 = t72 ** 2
           t79 = 0.1e1_dp + t74 + t76 * t77
           t80 = 0.1e1_dp / t79
           t81 = t75 * t80
           t84 = 0.1e1_dp + 0.66725e-1_dp * t73 * t81
           t85 = log(t84)
           epsilon_cGGA = e_c_u_0 + gamma_var * t85
           mu = 0.2224166667e-1_dp * t57
           kf = k_f
           ex_unif = -0.3e1_dp / 0.4e1_dp * t5 * kf
           t88 = 0.1e1_dp / kf
           t89 = my_norm_drho * t88
           s = t89 * t6 / 0.2e1_dp
           t91 = s ** 2
           t94 = 0.1e1_dp + 0.1243781095e1_dp * mu * t91
           Fx = 0.1804e1_dp - 0.804e0_dp / t94
           t97 = my_rho * ex_unif

           e_0(ii) = e_0(ii)&
                + t97 * Fx + my_rho * epsilon_cGGA

           t98 = t8 ** 2
           t100 = 0.1e1_dp / t98 * t5
           t101 = my_rho ** 2
           t102 = 0.1e1_dp / t101
           rsrho = -t4 * t100 * t102 / 0.12e2_dp
           t108 = t19 ** 2
           t109 = 0.1e1_dp / t108
           t110 = t11 * t109
           t111 = 0.1e1_dp / t12
           t112 = t111 * rsrho
           t117 = rs ** 0.10e1_dp
           t120 = 0.3797850000e1_dp * t112 + 0.35876e1_dp * rsrho + &
                0.2457300000e1_dp * t12 * rsrho + 0.985880e0_dp * t117 * rsrho
           t121 = 0.1e1_dp / t22
           t122 = t120 * t121
           e_c_u_0rho = -0.1328829340e-1_dp * rsrho * t23 + 1.e0_dp&
                * t110 * t122
           t125 = t60 ** 2
           k_frho = t1 / t125 * t57 / 0.3e1_dp
           t129 = 0.1e1_dp / t62
           k_srho = t129 * k_frho * t5
           t131 = k_s ** 2
           t132 = 0.1e1_dp / t131
           t133 = my_norm_drho * t132
           t134 = t6 * k_srho
           t136 = t64 * t102
           trho = -t133 * t134 / 0.2e1_dp - t136 / 0.2e1_dp
           t138 = gamma_var ** 2
           t140 = t69 ** 2
           t141 = 0.1e1_dp / t140
           t142 = 0.1e1_dp / t138 * t141
           Arho = 0.66725e-1_dp * t142 * e_c_u_0rho * t68
           t145 = t66 * t
           t146 = t81 * trho
           t149 = Arho * t72
           t150 = A * t
           t152 = 0.2e1_dp * t150 * trho
           t153 = t149 + t152
           t154 = t153 * t80
           t157 = t79 ** 2
           t158 = 0.1e1_dp / t157
           t159 = t75 * t158
           t160 = A * t77
           t163 = t72 * t
           t164 = t76 * t163
           t167 = t149 + t152 + 0.2e1_dp * t160 * Arho + 0.4e1_dp * t164 * trho
           t168 = t159 * t167
           t171 = 0.133450e0_dp * t145 * t146 + 0.66725e-1_dp * t73 * t154 - &
                0.66725e-1_dp * t73 * t168
           t172 = gamma_var * t171
           t173 = 0.1e1_dp / t84
           epsilon_cGGArho = e_c_u_0rho + t172 * t173
           kfrho = k_frho
           ex_unifrho = -0.3e1_dp / 0.4e1_dp * t5 * kfrho
           t177 = kf ** 2
           t178 = 0.1e1_dp / t177
           t179 = my_norm_drho * t178
           t180 = t6 * kfrho
           t182 = t89 * t102
           srho = -t179 * t180 / 0.2e1_dp - t182 / 0.2e1_dp
           t184 = t94 ** 2
           t186 = 0.1e1_dp / t184 * mu
           Fxrho = 0.2e1_dp * t186 * s * srho
           t190 = my_rho * ex_unifrho

           e_rho(ii) = e_rho(ii)&
                +ex_unif * Fx + t190 * Fx + t97 * Fxrho + epsilon_cGGA + &
                my_rho * epsilon_cGGArho

           tnorm_drho = t63 * t6 / 0.2e1_dp
           t194 = t81 * tnorm_drho
           t197 = t66 * t163
           t206 = 0.2e1_dp * t150 * tnorm_drho + 0.4e1_dp * t164 * tnorm_drho
           t210 = 0.133450e0_dp * t145 * t194 + 0.133450e0_dp * t197 * A * &
                tnorm_drho * t80 - 0.66725e-1_dp * t73 * t159 * t206
           t211 = gamma_var * t210
           Hnorm_drho = t211 * t173
           snorm_drho = t88 * t6 / 0.2e1_dp
           Fxnorm_drho = 0.2e1_dp * t186 * s * snorm_drho


           e_ndrho(ii) = e_ndrho(ii)&
                +t97 * Fxnorm_drho + my_rho * Hnorm_drho

           t217 = 0.1e1_dp / t98 / t7 * t58
           t218 = t101 ** 2
           t219 = 0.1e1_dp / t218
           t224 = 0.1e1_dp / t101 / my_rho
           rsrhorho = -t4 * t217 * t219 / 0.18e2_dp + t4 * t100 *&
                t224 / 0.6e1_dp
           t229 = 0.1328829340e-1_dp * rsrhorho * t23
           t230 = rsrho * t109
           t234 = 0.1e1_dp / t108 / t19
           t235 = t11 * t234
           t236 = t120 ** 2
           t237 = t236 * t121
           t240 = 0.1e1_dp / t15
           t241 = rsrho ** 2
           t242 = t240 * t241
           t244 = t111 * rsrhorho
           t245 = 0.3797850000e1_dp * t244
           t246 = 0.35876e1_dp * rsrhorho
           t250 = 0.2457300000e1_dp * t12 * rsrhorho
           t253 = 0.985880e0_dp * t117 * rsrhorho
           t254 = -0.1898925000e1_dp * t242 + t245 + t246 + 0.1228650000e1_dp&
                * t111 * t241 + t250 + 0.9858800e0_dp * t241 + t253
           t255 = t254 * t121
           t258 = t108 ** 2
           t259 = 0.1e1_dp / t258
           t260 = t11 * t259
           t261 = t22 ** 2
           t262 = 0.1e1_dp / t261
           t263 = t236 * t262
           e_c_u_0rhorho = -t229 + 0.4274000000e0_dp * t230 * t122 - &
                0.2000000000e1_dp * t235 * t237 + 1.e0_dp * t110 * t255 + &
                0.1608182432e2_dp * t260 * t263
           e_c_u_01rho = e_c_u_0rho
           t269 = t57 ** 2
           k_frhorho = -0.2e1_dp / 0.9e1_dp * t1 / t125 / t59 * t269
           t273 = 0.1e1_dp / t62 / t61
           t274 = k_frho ** 2
           t279 = t129 * k_frhorho * t5
           k_srhorho = -t273 * t274 * t58 / 0.2e1_dp + t279
           k_s1rho = k_srho
           t281 = 0.1e1_dp / t131 / k_s
           t282 = my_norm_drho * t281
           t285 = t102 * k_srho
           t287 = t133 * t285 / 0.2e1_dp
           t288 = t6 * k_srhorho
           t291 = t102 * k_s1rho
           t293 = t133 * t291 / 0.2e1_dp
           t294 = t64 * t224
           trhorho = t282 * t134 * k_s1rho + t287 - t133 * t288 / 0.2e1_dp + &
                t293 + t294
           t295 = t6 * k_s1rho
           t1rho = -t133 * t295 / 0.2e1_dp - t136 / 0.2e1_dp
           t299 = 0.1e1_dp / t138 / gamma_var
           t301 = 0.1e1_dp / t140 / t69
           t302 = t299 * t301
           t303 = t68 ** 2
           t304 = e_c_u_0rho * t303
           t311 = t299 * t141
           Arhorho = 0.133450e0_dp * t302 * t304 * e_c_u_01rho + &
                0.66725e-1_dp * t142 * e_c_u_0rhorho * t68 - 0.66725e-1_dp * t311 * &
                e_c_u_0rho * e_c_u_01rho * t68
           A1rho = 0.66725e-1_dp * t142 * e_c_u_01rho * t68
           t318 = t66 * t1rho
           t321 = A1rho * t72
           t323 = 0.2e1_dp * t150 * t1rho
           t324 = t321 + t323
           t325 = t324 * t80
           t326 = t325 * trho
           t329 = t145 * t75
           t330 = t158 * trho
           t335 = t321 + t323 + 0.2e1_dp * t160 * A1rho + 0.4e1_dp * t164 * t1rho
           t336 = t330 * t335
           t339 = t81 * trhorho
           t342 = t154 * t1rho
           t345 = Arhorho * t72
           t346 = Arho * t
           t348 = 0.2e1_dp * t346 * t1rho
           t349 = A1rho * t
           t351 = 0.2e1_dp * t349 * trho
           t352 = A * t1rho
           t354 = 0.2e1_dp * t352 * trho
           t356 = 0.2e1_dp * t150 * trhorho
           t357 = t345 + t348 + t351 + t354 + t356
           t358 = t357 * t80
           t361 = t153 * t158
           t365 = t158 * t167
           t366 = t365 * t1rho
           t369 = t324 * t158
           t373 = t73 * t75
           t375 = 0.1e1_dp / t157 / t79
           t376 = t375 * t167
           t377 = t376 * t335
           t380 = A1rho * t77
           t383 = A * t163
           t384 = Arho * t1rho
           t389 = trho * A1rho
           t392 = t76 * t72
           t393 = trho * t1rho
           t398 = t345 + t348 + t351 + t354 + t356 + 0.2e1_dp * t380 * Arho +&
                0.8e1_dp * t383 * t384 + 0.2e1_dp * t160 * Arhorho + 0.8e1_dp * &
                t383 * t389 + 0.12e2_dp * t392 * t393 + 0.4e1_dp * t164 * trhorho
           t402 = 0.133450e0_dp * t318 * t146 + 0.133450e0_dp * t145 * t326 -&
                0.133450e0_dp * t329 * t336 + 0.133450e0_dp * t145 * t339 + &
                0.133450e0_dp * t145 * t342 + 0.66725e-1_dp * t73 * t358 - &
                0.66725e-1_dp * t73 * t361 * t335 - 0.133450e0_dp * t329 * t366 - &
                0.66725e-1_dp * t73 * t369 * t167 + 0.133450e0_dp * t373 * t377 - &
                0.66725e-1_dp * t73 * t159 * t398
           t403 = gamma_var * t402
           t405 = t84 ** 2
           t406 = 0.1e1_dp / t405
           t407 = t81 * t1rho
           t415 = 0.133450e0_dp * t145 * t407 + 0.66725e-1_dp * t73 * t325 - &
                0.66725e-1_dp * t73 * t159 * t335
           t416 = t406 * t415
           t418 = gamma_var * t415
           epsilon_cGGArhorho = e_c_u_0rhorho + t403 * t173 - t172 * t416
           kfrhorho = k_frhorho
           ex_unifrhorho = -0.3e1_dp / 0.4e1_dp * t5 * kfrhorho
           ex_unif1rho = ex_unifrho
           t423 = 0.1e1_dp / t177 / kf
           t424 = my_norm_drho * t423
           t425 = kfrho ** 2
           t426 = t6 * t425
           t428 = t102 * kfrho
           t429 = t179 * t428
           t430 = t6 * kfrhorho
           t432 = t179 * t430 / 0.2e1_dp
           t433 = t89 * t224
           srhorho = t424 * t426 + t429 - t432 + t433
           s1rho = srho
           t436 = mu ** 2
           t437 = 0.1e1_dp / t184 / t94 * t436
           t438 = t91 * srho
           t442 = s1rho * srho
           Fxrhorho = -0.9950248765e1_dp * t437 * t438 * s1rho + &
                0.2e1_dp * t186 * t442 + 0.2e1_dp * t186 * s * &
                srhorho
           Fx1rho = 0.2e1_dp * t186 * s * s1rho
           t453 = my_rho * ex_unifrhorho
           t457 = my_rho * ex_unif1rho

           e_rho_rho(ii) = e_rho_rho(ii)&
                +ex_unif1rho * Fx + ex_unif * Fx1rho + ex_unifrho * &
                Fx + t453 * Fx + t190 * Fx1rho + ex_unif * Fxrho + t457 * Fxrho + &
                t97 * Fxrhorho + e_c_u_01rho + t418 * t173 + epsilon_cGGArho + my_rho *&
                epsilon_cGGArhorho

           t460 = t132 * t6
           t462 = t63 * t102
           tnorm_drhorho = -t460 * k_srho / 0.2e1_dp - t462 / 0.2e1_dp
           t464 = t66 * trho
           t467 = t154 * tnorm_drho
           t470 = t158 * tnorm_drho
           t471 = t470 * t167
           t474 = t81 * tnorm_drhorho
           t477 = t73 * A
           t478 = tnorm_drho * t80
           t479 = t478 * trho
           t482 = Arho * tnorm_drho
           t486 = A * tnorm_drhorho
           t490 = t197 * A
           t493 = t158 * t206
           t494 = t493 * trho
           t500 = t375 * t206
           t501 = t500 * t167
           t506 = A * trho
           t513 = tnorm_drho * trho
           t518 = 0.2e1_dp * t346 * tnorm_drho + 0.2e1_dp * t506 * tnorm_drho&
                + 0.2e1_dp * t150 * tnorm_drhorho + 0.8e1_dp * t383 * t482 + &
                0.12e2_dp * t392 * t513 + 0.4e1_dp * t164 * tnorm_drhorho
           t522 = 0.133450e0_dp * t464 * t194 + 0.133450e0_dp * t145 * t467 -&
                0.133450e0_dp * t329 * t471 + 0.133450e0_dp * t145 * t474 + &
                0.400350e0_dp * t477 * t479 + 0.133450e0_dp * t197 * t482 * t80 + &
                0.133450e0_dp * t197 * t486 * t80 - 0.133450e0_dp * t490 * t471 - &
                0.133450e0_dp * t329 * t494 - 0.66725e-1_dp * t73 * t361 * t206 + &
                0.133450e0_dp * t373 * t501 - 0.66725e-1_dp * t73 * t159 * t518
           t523 = gamma_var * t522
           t525 = t406 * t171
           Hnorm_drhorho = t523 * t173 - t211 * t525
           t528 = t178 * t6
           snorm_drhorho = -t528 * kfrho / 0.2e1_dp - t88 * t102 / 0.2e1_dp
           t532 = t91 * snorm_drho
           t536 = srho * snorm_drho
           Fxnorm_drhorho = -0.9950248765e1_dp * t437 * t532 * srho + &
                0.2e1_dp * t186 * t536 + 0.2e1_dp * t186 * s * &
                snorm_drhorho

           e_ndrho_rho(ii) = e_ndrho_rho(ii)&
                +ex_unif * Fxnorm_drho + t190 * Fxnorm_drho + &
                t97 * Fxnorm_drhorho + Hnorm_drho + my_rho * Hnorm_drhorho

           t545 = tnorm_drho ** 2
           t546 = t66 * t545
           t549 = A * t545
           t553 = t470 * t206
           t559 = t206 ** 2
           t566 = 0.2e1_dp * t549 + 0.12e2_dp * t392 * t545
           t571 = gamma_var * (0.133450e0_dp * t546 * t81 + 0.667250e0_dp * &
                t73 * t549 * t80 - 0.266900e0_dp * t329 * t553 - 0.266900e0_dp * &
                t490 * t553 + 0.133450e0_dp * t73 * t75 * t375 * t559 - &
                0.66725e-1_dp * t73 * t159 * t566)
           t573 = t210 ** 2
           t574 = gamma_var * t573
           Hnorm_drhonorm_drho = t571 * t173 - t574 * t406
           t576 = snorm_drho ** 2
           Fxnorm_drhonorm_drho = -0.9950248765e1_dp * t437 * t91 * t576 + &
                0.2e1_dp * t186 * t576

           e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
                +t97 * Fxnorm_drhonorm_drho + my_rho * &
                Hnorm_drhonorm_drho
        END IF
     END DO
  case default
     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_norm_drho = norm_drho(ii)

           t6 = 0.1e1_dp / my_rho
           t7 = t5 * t6
           t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
           rs = t4 * t8 / 0.4e1_dp
           t11 = 0.1e1_dp + 0.21370e0_dp * rs
           t12 = sqrt(rs)
           t15 = t12 * rs
           t17 = rs ** 0.20e1_dp
           t19 = 0.75957e1_dp * t12 + 0.35876e1_dp * rs + 0.16382e1_dp * t15 &
                + 0.49294e0_dp * t17
           t22 = 0.1e1_dp + 0.1608182432e2_dp / t19
           t23 = log(t22)
           e_c_u_0 = -0.62182e-1_dp * t11 * t23
           t36 = log(0.1e1_dp + 0.3216468318e2_dp / (0.141189e2_dp * t12 + &
                0.61977e1_dp * rs + 0.33662e1_dp * t15 + 0.62517e0_dp * t17))
           t49 = log(0.1e1_dp + 0.2960857464e1_dp / (0.10357e2_dp * t12 + &
                0.36231e1_dp * rs + 0.88026e0_dp * t15 + 0.49671e0_dp * t17))
           t51 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t55 = log(0.2e1_dp)
           t57 = 0.3141592654e1_dp ** 2
           t58 = 0.1e1_dp / t57
           gamma_var = (0.1e1_dp - t55) * t58
           t59 = t57 * my_rho
           t60 = t59 ** (0.1e1_dp / 0.3e1_dp)
           k_f = t1 * t60
           t61 = k_f * t5
           t62 = sqrt(t61)
           k_s = 0.2e1_dp * t62
           t63 = 0.1e1_dp / k_s
           t64 = my_norm_drho * t63
           t = t64 * t6 / 0.2e1_dp
           t66 = 0.1e1_dp / gamma_var
           t68 = exp(-e_c_u_0 * t66)
           t69 = -0.1e1_dp + t68
           A = 0.66725e-1_dp * t66 / t69
           t72 = t ** 2
           t73 = t66 * t72
           t74 = A * t72
           t75 = 0.1e1_dp + t74
           t76 = A ** 2
           t77 = t72 ** 2
           t79 = 0.1e1_dp + t74 + t76 * t77
           t80 = 0.1e1_dp / t79
           t81 = t75 * t80
           t84 = 0.1e1_dp + 0.66725e-1_dp * t73 * t81
           t85 = log(t84)
           epsilon_cGGA = e_c_u_0 + gamma_var * t85
           mu = 0.2224166667e-1_dp * t57
           kf = k_f
           ex_unif = -0.3e1_dp / 0.4e1_dp * t5 * kf
           t88 = 0.1e1_dp / kf
           t89 = my_norm_drho * t88
           s = t89 * t6 / 0.2e1_dp
           t91 = s ** 2
           t94 = 0.1e1_dp + 0.1243781095e1_dp * mu * t91
           Fx = 0.1804e1_dp - 0.804e0_dp / t94
           t97 = my_rho * ex_unif

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)&
                   + t97 * Fx + my_rho * epsilon_cGGA
           END IF

           t98 = t8 ** 2
           t100 = 0.1e1_dp / t98 * t5
           t101 = my_rho ** 2
           t102 = 0.1e1_dp / t101
           rsrho = -t4 * t100 * t102 / 0.12e2_dp
           t108 = t19 ** 2
           t109 = 0.1e1_dp / t108
           t110 = t11 * t109
           t111 = 0.1e1_dp / t12
           t112 = t111 * rsrho
           t117 = rs ** 0.10e1_dp
           t120 = 0.3797850000e1_dp * t112 + 0.35876e1_dp * rsrho + &
                0.2457300000e1_dp * t12 * rsrho + 0.985880e0_dp * t117 * rsrho
           t121 = 0.1e1_dp / t22
           t122 = t120 * t121
           e_c_u_0rho = -0.1328829340e-1_dp * rsrho * t23 + 1.e0_dp&
                * t110 * t122
           t125 = t60 ** 2
           k_frho = t1 / t125 * t57 / 0.3e1_dp
           t129 = 0.1e1_dp / t62
           k_srho = t129 * k_frho * t5
           t131 = k_s ** 2
           t132 = 0.1e1_dp / t131
           t133 = my_norm_drho * t132
           t134 = t6 * k_srho
           t136 = t64 * t102
           trho = -t133 * t134 / 0.2e1_dp - t136 / 0.2e1_dp
           t138 = gamma_var ** 2
           t140 = t69 ** 2
           t141 = 0.1e1_dp / t140
           t142 = 0.1e1_dp / t138 * t141
           Arho = 0.66725e-1_dp * t142 * e_c_u_0rho * t68
           t145 = t66 * t
           t146 = t81 * trho
           t149 = Arho * t72
           t150 = A * t
           t152 = 0.2e1_dp * t150 * trho
           t153 = t149 + t152
           t154 = t153 * t80
           t157 = t79 ** 2
           t158 = 0.1e1_dp / t157
           t159 = t75 * t158
           t160 = A * t77
           t163 = t72 * t
           t164 = t76 * t163
           t167 = t149 + t152 + 0.2e1_dp * t160 * Arho + 0.4e1_dp * t164 * trho
           t168 = t159 * t167
           t171 = 0.133450e0_dp * t145 * t146 + 0.66725e-1_dp * t73 * t154 - &
                0.66725e-1_dp * t73 * t168
           t172 = gamma_var * t171
           t173 = 0.1e1_dp / t84
           epsilon_cGGArho = e_c_u_0rho + t172 * t173
           kfrho = k_frho
           ex_unifrho = -0.3e1_dp / 0.4e1_dp * t5 * kfrho
           t177 = kf ** 2
           t178 = 0.1e1_dp / t177
           t179 = my_norm_drho * t178
           t180 = t6 * kfrho
           t182 = t89 * t102
           srho = -t179 * t180 / 0.2e1_dp - t182 / 0.2e1_dp
           t184 = t94 ** 2
           t186 = 0.1e1_dp / t184 * mu
           Fxrho = 0.2e1_dp * t186 * s * srho
           t190 = my_rho * ex_unifrho

           IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
              e_rho(ii) = e_rho(ii)&
                   +ex_unif * Fx + t190 * Fx + t97 * Fxrho + epsilon_cGGA + &
                   my_rho * epsilon_cGGArho
           END IF

           tnorm_drho = t63 * t6 / 0.2e1_dp
           t194 = t81 * tnorm_drho
           t197 = t66 * t163
           t206 = 0.2e1_dp * t150 * tnorm_drho + 0.4e1_dp * t164 * tnorm_drho
           t210 = 0.133450e0_dp * t145 * t194 + 0.133450e0_dp * t197 * A * &
                tnorm_drho * t80 - 0.66725e-1_dp * t73 * t159 * t206
           t211 = gamma_var * t210
           Hnorm_drho = t211 * t173
           snorm_drho = t88 * t6 / 0.2e1_dp
           Fxnorm_drho = 0.2e1_dp * t186 * s * snorm_drho


           IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
              e_ndrho(ii) = e_ndrho(ii)&
                   +t97 * Fxnorm_drho + my_rho * Hnorm_drho
           END IF

           t217 = 0.1e1_dp / t98 / t7 * t58
           t218 = t101 ** 2
           t219 = 0.1e1_dp / t218
           t224 = 0.1e1_dp / t101 / my_rho
           rsrhorho = -t4 * t217 * t219 / 0.18e2_dp + t4 * t100 *&
                t224 / 0.6e1_dp
           t229 = 0.1328829340e-1_dp * rsrhorho * t23
           t230 = rsrho * t109
           t234 = 0.1e1_dp / t108 / t19
           t235 = t11 * t234
           t236 = t120 ** 2
           t237 = t236 * t121
           t240 = 0.1e1_dp / t15
           t241 = rsrho ** 2
           t242 = t240 * t241
           t244 = t111 * rsrhorho
           t245 = 0.3797850000e1_dp * t244
           t246 = 0.35876e1_dp * rsrhorho
           t250 = 0.2457300000e1_dp * t12 * rsrhorho
           t253 = 0.985880e0_dp * t117 * rsrhorho
           t254 = -0.1898925000e1_dp * t242 + t245 + t246 + 0.1228650000e1_dp&
                * t111 * t241 + t250 + 0.9858800e0_dp * t241 + t253
           t255 = t254 * t121
           t258 = t108 ** 2
           t259 = 0.1e1_dp / t258
           t260 = t11 * t259
           t261 = t22 ** 2
           t262 = 0.1e1_dp / t261
           t263 = t236 * t262
           e_c_u_0rhorho = -t229 + 0.4274000000e0_dp * t230 * t122 - &
                0.2000000000e1_dp * t235 * t237 + 1.e0_dp * t110 * t255 + &
                0.1608182432e2_dp * t260 * t263
           e_c_u_01rho = e_c_u_0rho
           t269 = t57 ** 2
           k_frhorho = -0.2e1_dp / 0.9e1_dp * t1 / t125 / t59 * t269
           t273 = 0.1e1_dp / t62 / t61
           t274 = k_frho ** 2
           t279 = t129 * k_frhorho * t5
           k_srhorho = -t273 * t274 * t58 / 0.2e1_dp + t279
           k_s1rho = k_srho
           t281 = 0.1e1_dp / t131 / k_s
           t282 = my_norm_drho * t281
           t285 = t102 * k_srho
           t287 = t133 * t285 / 0.2e1_dp
           t288 = t6 * k_srhorho
           t291 = t102 * k_s1rho
           t293 = t133 * t291 / 0.2e1_dp
           t294 = t64 * t224
           trhorho = t282 * t134 * k_s1rho + t287 - t133 * t288 / 0.2e1_dp + &
                t293 + t294
           t295 = t6 * k_s1rho
           t1rho = -t133 * t295 / 0.2e1_dp - t136 / 0.2e1_dp
           t299 = 0.1e1_dp / t138 / gamma_var
           t301 = 0.1e1_dp / t140 / t69
           t302 = t299 * t301
           t303 = t68 ** 2
           t304 = e_c_u_0rho * t303
           t311 = t299 * t141
           Arhorho = 0.133450e0_dp * t302 * t304 * e_c_u_01rho + &
                0.66725e-1_dp * t142 * e_c_u_0rhorho * t68 - 0.66725e-1_dp * t311 * &
                e_c_u_0rho * e_c_u_01rho * t68
           A1rho = 0.66725e-1_dp * t142 * e_c_u_01rho * t68
           t318 = t66 * t1rho
           t321 = A1rho * t72
           t323 = 0.2e1_dp * t150 * t1rho
           t324 = t321 + t323
           t325 = t324 * t80
           t326 = t325 * trho
           t329 = t145 * t75
           t330 = t158 * trho
           t335 = t321 + t323 + 0.2e1_dp * t160 * A1rho + 0.4e1_dp * t164 * t1rho
           t336 = t330 * t335
           t339 = t81 * trhorho
           t342 = t154 * t1rho
           t345 = Arhorho * t72
           t346 = Arho * t
           t348 = 0.2e1_dp * t346 * t1rho
           t349 = A1rho * t
           t351 = 0.2e1_dp * t349 * trho
           t352 = A * t1rho
           t354 = 0.2e1_dp * t352 * trho
           t356 = 0.2e1_dp * t150 * trhorho
           t357 = t345 + t348 + t351 + t354 + t356
           t358 = t357 * t80
           t361 = t153 * t158
           t365 = t158 * t167
           t366 = t365 * t1rho
           t369 = t324 * t158
           t373 = t73 * t75
           t375 = 0.1e1_dp / t157 / t79
           t376 = t375 * t167
           t377 = t376 * t335
           t380 = A1rho * t77
           t383 = A * t163
           t384 = Arho * t1rho
           t389 = trho * A1rho
           t392 = t76 * t72
           t393 = trho * t1rho
           t398 = t345 + t348 + t351 + t354 + t356 + 0.2e1_dp * t380 * Arho +&
                0.8e1_dp * t383 * t384 + 0.2e1_dp * t160 * Arhorho + 0.8e1_dp * &
                t383 * t389 + 0.12e2_dp * t392 * t393 + 0.4e1_dp * t164 * trhorho
           t402 = 0.133450e0_dp * t318 * t146 + 0.133450e0_dp * t145 * t326 -&
                0.133450e0_dp * t329 * t336 + 0.133450e0_dp * t145 * t339 + &
                0.133450e0_dp * t145 * t342 + 0.66725e-1_dp * t73 * t358 - &
                0.66725e-1_dp * t73 * t361 * t335 - 0.133450e0_dp * t329 * t366 - &
                0.66725e-1_dp * t73 * t369 * t167 + 0.133450e0_dp * t373 * t377 - &
                0.66725e-1_dp * t73 * t159 * t398
           t403 = gamma_var * t402
           t405 = t84 ** 2
           t406 = 0.1e1_dp / t405
           t407 = t81 * t1rho
           t415 = 0.133450e0_dp * t145 * t407 + 0.66725e-1_dp * t73 * t325 - &
                0.66725e-1_dp * t73 * t159 * t335
           t416 = t406 * t415
           t418 = gamma_var * t415
           epsilon_cGGArhorho = e_c_u_0rhorho + t403 * t173 - t172 * t416
           kfrhorho = k_frhorho
           ex_unifrhorho = -0.3e1_dp / 0.4e1_dp * t5 * kfrhorho
           ex_unif1rho = ex_unifrho
           t423 = 0.1e1_dp / t177 / kf
           t424 = my_norm_drho * t423
           t425 = kfrho ** 2
           t426 = t6 * t425
           t428 = t102 * kfrho
           t429 = t179 * t428
           t430 = t6 * kfrhorho
           t432 = t179 * t430 / 0.2e1_dp
           t433 = t89 * t224
           srhorho = t424 * t426 + t429 - t432 + t433
           s1rho = srho
           t436 = mu ** 2
           t437 = 0.1e1_dp / t184 / t94 * t436
           t438 = t91 * srho
           t442 = s1rho * srho
           Fxrhorho = -0.9950248765e1_dp * t437 * t438 * s1rho + &
                0.2e1_dp * t186 * t442 + 0.2e1_dp * t186 * s * &
                srhorho
           Fx1rho = 0.2e1_dp * t186 * s * s1rho
           t453 = my_rho * ex_unifrhorho
           t457 = my_rho * ex_unif1rho

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              e_rho_rho(ii) = e_rho_rho(ii)&
                   +ex_unif1rho * Fx + ex_unif * Fx1rho + ex_unifrho * &
                   Fx + t453 * Fx + t190 * Fx1rho + ex_unif * Fxrho + t457 * Fxrho + &
                   t97 * Fxrhorho + e_c_u_01rho + t418 * t173 + epsilon_cGGArho + my_rho *&
                   epsilon_cGGArhorho
           END IF

           t460 = t132 * t6
           t462 = t63 * t102
           tnorm_drhorho = -t460 * k_srho / 0.2e1_dp - t462 / 0.2e1_dp
           t464 = t66 * trho
           t467 = t154 * tnorm_drho
           t470 = t158 * tnorm_drho
           t471 = t470 * t167
           t474 = t81 * tnorm_drhorho
           t477 = t73 * A
           t478 = tnorm_drho * t80
           t479 = t478 * trho
           t482 = Arho * tnorm_drho
           t486 = A * tnorm_drhorho
           t490 = t197 * A
           t493 = t158 * t206
           t494 = t493 * trho
           t500 = t375 * t206
           t501 = t500 * t167
           t506 = A * trho
           t513 = tnorm_drho * trho
           t518 = 0.2e1_dp * t346 * tnorm_drho + 0.2e1_dp * t506 * tnorm_drho&
                + 0.2e1_dp * t150 * tnorm_drhorho + 0.8e1_dp * t383 * t482 + &
                0.12e2_dp * t392 * t513 + 0.4e1_dp * t164 * tnorm_drhorho
           t522 = 0.133450e0_dp * t464 * t194 + 0.133450e0_dp * t145 * t467 -&
                0.133450e0_dp * t329 * t471 + 0.133450e0_dp * t145 * t474 + &
                0.400350e0_dp * t477 * t479 + 0.133450e0_dp * t197 * t482 * t80 + &
                0.133450e0_dp * t197 * t486 * t80 - 0.133450e0_dp * t490 * t471 - &
                0.133450e0_dp * t329 * t494 - 0.66725e-1_dp * t73 * t361 * t206 + &
                0.133450e0_dp * t373 * t501 - 0.66725e-1_dp * t73 * t159 * t518
           t523 = gamma_var * t522
           t525 = t406 * t171
           Hnorm_drhorho = t523 * t173 - t211 * t525
           t528 = t178 * t6
           snorm_drhorho = -t528 * kfrho / 0.2e1_dp - t88 * t102 / 0.2e1_dp
           t532 = t91 * snorm_drho
           t536 = srho * snorm_drho
           Fxnorm_drhorho = -0.9950248765e1_dp * t437 * t532 * srho + &
                0.2e1_dp * t186 * t536 + 0.2e1_dp * t186 * s * &
                snorm_drhorho

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              e_ndrho_rho(ii) = e_ndrho_rho(ii)&
                   +ex_unif * Fxnorm_drho + t190 * Fxnorm_drho + &
                   t97 * Fxnorm_drhorho + Hnorm_drho + my_rho * Hnorm_drhorho
           END IF

           t545 = tnorm_drho ** 2
           t546 = t66 * t545
           t549 = A * t545
           t553 = t470 * t206
           t559 = t206 ** 2
           t566 = 0.2e1_dp * t549 + 0.12e2_dp * t392 * t545
           t571 = gamma_var * (0.133450e0_dp * t546 * t81 + 0.667250e0_dp * &
                t73 * t549 * t80 - 0.266900e0_dp * t329 * t553 - 0.266900e0_dp * &
                t490 * t553 + 0.133450e0_dp * t73 * t75 * t375 * t559 - &
                0.66725e-1_dp * t73 * t159 * t566)
           t573 = t210 ** 2
           t574 = gamma_var * t573
           Hnorm_drhonorm_drho = t571 * t173 - t574 * t406
           t576 = snorm_drho ** 2
           Fxnorm_drhonorm_drho = -0.9950248765e1_dp * t437 * t91 * t576 + &
                0.2e1_dp * t186 * t576

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
                   +t97 * Fxnorm_drhonorm_drho + my_rho * &
                   Hnorm_drhonorm_drho
           END IF

           IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
              t586 = 0.1e1_dp / t57 / 0.3141592654e1_dp
              rsrhorhorho = -0.5e1_dp / 0.54e2_dp * t4 / t98 / t58 / t102 &
                   * t586 / t218 / t101 + t4 * t217 / t218 / my_rho / 0.3e1_dp - &
                   t4 * t100 * t219 / 0.2e1_dp
              rs2rho = rsrho
              t603 = rsrhorho * t109
              t611 = 0.3797850000e1_dp * t111 * rs2rho + 0.35876e1_dp * rs2rho +&
                   0.2457300000e1_dp * t12 * rs2rho + 0.985880e0_dp * t117 * rs2rho
              t612 = t611 * t121
              t618 = t122 * t611
              t628 = -0.1898925000e1_dp * t240 * rsrho * rs2rho + t245 + t246 + &
                   0.1228650000e1_dp * t112 * rs2rho + t250 + 0.9858800e0_dp * rsrho * &
                   rs2rho + t253
              t629 = t628 * t121
              t634 = t120 * t262 * t611
              t652 = rs2rho * t109
              t658 = rs ** 2
              t664 = t240 * rsrhorho
              t686 = 0.2848387500e1_dp / t12 / t658 * t241 * rs2rho - &
                   0.3797850000e1_dp * t664 * rsrho - 0.1898925000e1_dp * t664 * rs2rho&
                   + 0.3797850000e1_dp * t111 * rsrhorhorho + 0.35876e1_dp * &
                   rsrhorhorho - 0.6143250000e0_dp * t242 * rs2rho + 0.2457300000e1_dp &
                   * t244 * rsrho + 0.1228650000e1_dp * t244 * rs2rho + &
                   0.2457300000e1_dp * t12 * rsrhorhorho + 0.19717600e1_dp * rsrhorho *&
                   rsrho + 0.9858800e0_dp * rsrhorho * rs2rho + 0.985880e0_dp * t117 *&
                   rsrhorhorho
              e_c_u_0rhorhorho = -0.1328829340e-1_dp * rsrhorhorho * t23 + &
                   0.2137000000e0_dp * t603 * t612 + 0.4274000000e0_dp * t603 * t122 - &
                   0.8548000000e0_dp * rsrho * t234 * t618 + 0.4274000000e0_dp * t230 *&
                   t629 + 0.6873371714e1_dp * rsrho * t259 * t634 - 0.4274000000e0_dp &
                   * rs2rho * t234 * t237 + 0.6000000000e1_dp * t260 * t237 * t611 - &
                   0.4000000000e1_dp * t235 * t629 * t120 - 0.9649094592e2_dp * t11 / &
                   t258 / t19 * t263 * t611 + 0.2137000000e0_dp * t652 * t255 - &
                   0.2000000000e1_dp * t235 * t255 * t611 + 1.e0_dp * t110 * &
                   t686 * t121 + 0.1608182432e2_dp * t260 * t254 * t262 * t611 + &
                   0.3436685857e1_dp * rs2rho * t259 * t263 + 0.3216364864e2_dp * t260 &
                   * t628 * t262 * t120 + 0.5172501469e3_dp * t11 / t258 / t108 * t236 &
                   / t261 / t22 * t611
              e_c_u_0rho1rho = -t229 + 0.2137000000e0_dp * t230 * t612 + &
                   0.2137000000e0_dp * t652 * t122 - 0.2000000000e1_dp * t235 * t618 + &
                   1.e0_dp * t110 * t629 + 0.1608182432e2_dp * t260 * t634
              e_c_u_01rhorho = e_c_u_0rho1rho
              e_c_u_02rho = -0.1328829340e-1_dp * rs2rho * t23 + &
                   1.e0_dp * t110 * t612
              k_frhorhorho = 0.10e2_dp / 0.27e2_dp * t1 / t125 / t101 * t57
              k_f2rho = kfrho
              t730 = k_f ** 2
              t738 = t273 * k_frhorho
              t741 = t58 * k_f2rho
              k_srho1rho = -t273 * k_frho * t741 / 0.2e1_dp + t279
              k_s1rhorho = k_srho1rho
              k_s2rho = t129 * k_f2rho * t5
              t750 = t131 ** 2
              t754 = k_srho * k_s1rho
              t760 = t6 * k_srho1rho
              trhorhorho = -0.3e1_dp * my_norm_drho / t750 * t6 * t754 * k_s2rho - &
                   t282 * t285 * k_s1rho + t282 * t760 * k_s1rho + t282 * t134 * &
                   k_s1rhorho - t282 * t285 * k_s2rho - t133 * t224 * k_srho + t133 * &
                   t102 * k_srho1rho / 0.2e1_dp + t282 * t288 * k_s2rho + t133 * t102 *&
                   k_srhorho / 0.2e1_dp - t133 * t6 * (0.3e1_dp / 0.4e1_dp / t62 / &
                   t730 / t58 * t274 * t586 * k_f2rho - t738 * t58 * k_frho - t738 * &
                   t741 / 0.2e1_dp + t129 * k_frhorhorho * t5) / 0.2e1_dp - t282 * t291&
                   * k_s2rho - t133 * t224 * k_s1rho + t133 * t102 * k_s1rhorho / &
                   0.2e1_dp - t133 * t224 * k_s2rho - 0.3e1_dp * t64 * t219
              t797 = t133 * t102 * k_s2rho / 0.2e1_dp
              trho1rho = t282 * t134 * k_s2rho + t287 - t133 * t760 / 0.2e1_dp +&
                   t797 + t294
              t1rhorho = t282 * t295 * k_s2rho + t293 - t133 * t6 * k_s1rhorho /&
                   0.2e1_dp + t797 + t294
              t2rho = -t133 * t6 * k_s2rho / 0.2e1_dp - t136 / 0.2e1_dp
              t806 = t138 ** 2
              t807 = 0.1e1_dp / t806
              t808 = t140 ** 2
              t824 = t303 * e_c_u_01rho * e_c_u_02rho
              t852 = e_c_u_01rho * e_c_u_02rho * t68
              Arhorhorho = 0.400350e0_dp * t807 / t808 * e_c_u_0rho * t303 * t68&
                   * e_c_u_01rho * e_c_u_02rho + 0.133450e0_dp * t302 * e_c_u_0rho1rho&
                   * t303 * e_c_u_01rho - 0.400350e0_dp * t807 * t301 * e_c_u_0rho * &
                   t824 + 0.133450e0_dp * t302 * t304 * e_c_u_01rhorho + 0.133450e0_dp &
                   * t302 * e_c_u_0rhorho * t303 * e_c_u_02rho + 0.66725e-1_dp * t142 *&
                   e_c_u_0rhorhorho * t68 - 0.66725e-1_dp * t311 * e_c_u_0rhorho * &
                   e_c_u_02rho * t68 - 0.66725e-1_dp * t311 * e_c_u_0rho1rho * &
                   e_c_u_01rho * t68 - 0.66725e-1_dp * t311 * e_c_u_0rho * &
                   e_c_u_01rhorho * t68 + 0.66725e-1_dp * t807 * t141 * e_c_u_0rho * &
                   t852
              Arho1rho = 0.133450e0_dp * t302 * t304 * e_c_u_02rho + &
                   0.66725e-1_dp * t142 * e_c_u_0rho1rho * t68 - 0.66725e-1_dp * t311 *&
                   e_c_u_0rho * e_c_u_02rho * t68
              A1rhorho = 0.133450e0_dp * t302 * t824 + 0.66725e-1_dp * t142 * &
                   e_c_u_01rhorho * t68 - 0.66725e-1_dp * t311 * t852
              A2rho = 0.66725e-1_dp * t142 * e_c_u_02rho * t68
              t874 = A2rho * t72
              t876 = 0.2e1_dp * t150 * t2rho
              t877 = t874 + t876
              t878 = t877 * t80
              t879 = t878 * trho
              t882 = t66 * t2rho
              t890 = t874 + t876 + 0.2e1_dp * t160 * A2rho + 0.4e1_dp * t164 * t2rho
              t891 = t335 * t890
              t898 = t73 * t153
              t900 = t375 * t335 * t890
              t903 = Arho1rho * t72
              t905 = 0.2e1_dp * t346 * t2rho
              t906 = A2rho * t
              t908 = 0.2e1_dp * t906 * trho
              t909 = A * t2rho
              t911 = 0.2e1_dp * t909 * trho
              t913 = 0.2e1_dp * t150 * trho1rho
              t914 = t903 + t905 + t908 + t911 + t913
              t915 = t914 * t80
              t921 = t882 * t75
              t927 = A2rho * t77
              t930 = Arho * t2rho
              t943 = t903 + t905 + t908 + t911 + t913 + 0.2e1_dp * t927 * Arho +&
                   0.8e1_dp * t383 * t930 + 0.2e1_dp * t160 * Arho1rho + 0.8e1_dp * &
                   t383 * trho * A2rho + 0.12e2_dp * t392 * trho * t2rho + 0.4e1_dp * &
                   t164 * trho1rho
              t951 = 0.133450e0_dp * t318 * t879 + 0.133450e0_dp * t882 * t339 +&
                   0.266900e0_dp * t329 * t375 * trho * t891 + 0.133450e0_dp * t145 * &
                   t325 * trho1rho + 0.133450e0_dp * t898 * t900 + 0.133450e0_dp * t145&
                   * t915 * t1rho + 0.133450e0_dp * t882 * t326 - 0.133450e0_dp * t921&
                   * t336 + 0.133450e0_dp * t145 * t154 * t1rhorho + 0.133450e0_dp * &
                   t373 * t375 * t943 * t335 + 0.133450e0_dp * t73 * t877 * t377
              t956 = A1rhorho * t72
              t958 = 0.2e1_dp * t349 * t2rho
              t960 = 0.2e1_dp * t906 * t1rho
              t962 = 0.2e1_dp * t909 * t1rho
              t964 = 0.2e1_dp * t150 * t1rhorho
              t967 = A1rho * t2rho
              t972 = t1rho * A2rho
              t980 = t956 + t958 + t960 + t962 + t964 + 0.2e1_dp * t927 * A1rho &
                   + 0.8e1_dp * t383 * t967 + 0.2e1_dp * t160 * A1rhorho + 0.8e1_dp * &
                   t383 * t972 + 0.12e2_dp * t392 * t1rho * t2rho + 0.4e1_dp * t164 * &
                   t1rhorho
              t984 = t956 + t958 + t960 + t962 + t964
              t992 = t357 * t158
              t996 = t145 * t877
              t1004 = t145 * t153
              t1006 = t158 * t335 * t2rho
              t1019 = -0.133450e0_dp * t329 * t158 * trho1rho * t335 - &
                   0.66725e-1_dp * t73 * t361 * t980 - 0.66725e-1_dp * t73 * t984 * &
                   t158 * t167 + 0.133450e0_dp * t145 * t358 * t2rho - 0.66725e-1_dp * &
                   t73 * t992 * t890 - 0.133450e0_dp * t996 * t336 - 0.133450e0_dp * &
                   t996 * t366 - 0.133450e0_dp * t329 * t330 * t980 - 0.133450e0_dp * &
                   t1004 * t1006 + 0.266900e0_dp * t329 * t376 * t1rho * t890 - &
                   0.133450e0_dp * t329 * t365 * t1rhorho - 0.66725e-1_dp * t73 * t369 &
                   * t943
              t1021 = t984 * t80
              t1049 = t73 * t324
              t1050 = t376 * t890
              t1053 = t318 * t75
              t1054 = t330 * t890
              t1058 = t158 * t1rho * t890
              t1061 = 0.133450e0_dp * t145 * t1021 * trho - 0.133450e0_dp * t921&
                   * t366 + 0.266900e0_dp * t329 * t376 * t335 * t2rho + 0.133450e0_dp&
                   * t373 * t375 * t398 * t890 - 0.66725e-1_dp * t73 * t914 * t158 * &
                   t335 + 0.133450e0_dp * t145 * t878 * trhorho - 0.133450e0_dp * t329 &
                   * t158 * trhorho * t890 + 0.133450e0_dp * t66 * t1rhorho * t146 + &
                   0.133450e0_dp * t1049 * t1050 - 0.133450e0_dp * t1053 * t1054 - &
                   0.133450e0_dp * t1004 * t1058
              t1068 = t145 * t324
              t1071 = t967 * trho
              t1074 = Arhorho * t
              t1076 = 0.2e1_dp * t1074 * t2rho
              t1079 = 0.2e1_dp * Arho1rho * t * t1rho
              t1080 = t930 * t1rho
              t1081 = 0.2e1_dp * t1080
              t1083 = 0.2e1_dp * t346 * t1rhorho
              t1086 = 0.2e1_dp * A1rhorho * t * trho
              t1087 = 0.2e1_dp * t1071
              t1089 = 0.2e1_dp * t349 * trho1rho
              t1090 = t972 * trho
              t1091 = 0.2e1_dp * t1090
              t1094 = 0.2e1_dp * A * t1rhorho * trho
              t1096 = 0.2e1_dp * t352 * trho1rho
              t1098 = 0.2e1_dp * t906 * trhorho
              t1100 = 0.2e1_dp * t909 * trhorho
              t1102 = 0.2e1_dp * t150 * trhorhorho
              t1105 = t76 * t
              t1111 = 0.24e2_dp * t74 * t1071 + t1076 + t1079 + t1081 + t1083 + &
                   t1086 + t1087 + t1089 + t1091 + t1094 + t1096 + t1098 + t1100 + &
                   t1102 + 0.24e2_dp * t74 * t1090 + 0.24e2_dp * t1105 * t393 * t2rho +&
                   0.24e2_dp * t74 * t1080
              t1112 = Arhorhorho * t72
              t1130 = A2rho * t163
              t1133 = A1rho * t163
              t1159 = t1112 + 0.2e1_dp * t380 * Arho1rho + 0.2e1_dp * t927 * &
                   Arhorho + 0.4e1_dp * t164 * trhorhorho + 0.2e1_dp * t160 * &
                   Arhorhorho + 0.8e1_dp * t383 * trho * A1rhorho + 0.12e2_dp * t392 * &
                   trho * t1rhorho + 0.8e1_dp * t383 * trho1rho * A1rho + 0.8e1_dp * &
                   t1130 * t389 + 0.8e1_dp * t1133 * t930 + 0.12e2_dp * t392 * trho1rho&
                   * t1rho + 0.8e1_dp * t383 * Arho1rho * t1rho + 0.2e1_dp * A1rhorho &
                   * t77 * Arho + 0.8e1_dp * t1130 * t384 + 0.12e2_dp * t392 * trhorho &
                   * t2rho + 0.8e1_dp * t383 * Arho * t1rhorho + 0.8e1_dp * t383 * &
                   Arhorho * t2rho + 0.8e1_dp * t383 * trhorho * A2rho
              t1166 = t157 ** 2
              t1167 = 0.1e1_dp / t1166
              t1172 = t877 * t158
              t1176 = t1112 + t1076 + t1079 + t1081 + t1083 + t1086 + t1087 + &
                   t1089 + t1091 + t1094 + t1096 + t1098 + t1100 + t1102
              t1180 = t365 * t2rho
              t1187 = t81 * trho1rho
              t1194 = 0.133450e0_dp * t145 * t81 * trhorhorho + 0.133450e0_dp * &
                   t373 * t376 * t980 - 0.133450e0_dp * t1068 * t1054 - 0.66725e-1_dp *&
                   t73 * t159 * (t1111 + t1159) + 0.133450e0_dp * t882 * t342 - &
                   0.400350e0_dp * t373 * t1167 * t167 * t891 - 0.66725e-1_dp * t73 * &
                   t1172 * t398 + 0.66725e-1_dp * t73 * t1176 * t80 - 0.133450e0_dp * &
                   t1068 * t1180 - 0.133450e0_dp * t329 * t158 * t398 * t2rho + &
                   0.133450e0_dp * t318 * t1187 - 0.133450e0_dp * t329 * t158 * t943 * &
                   t1rho
              t1207 = 0.133450e0_dp * t145 * t81 * t2rho + 0.66725e-1_dp * t73 *&
                   t878 - 0.66725e-1_dp * t73 * t159 * t890
              t1208 = t406 * t1207
              t1236 = 0.133450e0_dp * t882 * t146 + 0.133450e0_dp * t145 * t879 &
                   - 0.133450e0_dp * t329 * t1054 + 0.133450e0_dp * t145 * t1187 + &
                   0.133450e0_dp * t145 * t154 * t2rho + 0.66725e-1_dp * t73 * t915 - &
                   0.66725e-1_dp * t73 * t361 * t890 - 0.133450e0_dp * t329 * t1180 - &
                   0.66725e-1_dp * t73 * t1172 * t167 + 0.133450e0_dp * t373 * t1050 - &
                   0.66725e-1_dp * t73 * t159 * t943
              t1237 = gamma_var * t1236
              t1240 = 0.1e1_dp / t405 / t84
              t1273 = 0.133450e0_dp * t882 * t407 + 0.133450e0_dp * t145 * t878 &
                   * t1rho - 0.133450e0_dp * t329 * t1058 + 0.133450e0_dp * t145 * t81 &
                   * t1rhorho + 0.133450e0_dp * t145 * t325 * t2rho + 0.66725e-1_dp * &
                   t73 * t1021 - 0.66725e-1_dp * t73 * t369 * t890 - 0.133450e0_dp * &
                   t329 * t1006 - 0.66725e-1_dp * t73 * t1172 * t335 + 0.133450e0_dp * &
                   t373 * t900 - 0.66725e-1_dp * t73 * t159 * t980
              kfrhorhorho = k_frhorhorho
              kf2rho = k_f2rho
              ex_unifrho1rho = ex_unifrhorho
              ex_unif1rhorho = ex_unifrho1rho
              ex_unif2rho = -0.3e1_dp / 0.4e1_dp * t5 * kf2rho
              t1286 = t177 ** 2
              srho1rho = t424 * t180 * kf2rho + t429 / 0.2e1_dp - t432 + t179 * &
                   t102 * kf2rho / 0.2e1_dp + t433
              s1rhorho = srho1rho
              s2rho = -t179 * t6 * kf2rho / 0.2e1_dp - t182 / 0.2e1_dp
              t1324 = t184 ** 2
              t1327 = 0.1e1_dp / t1324 * t436 * mu
              t1328 = t91 * s
              t1329 = t1327 * t1328
              t1330 = t442 * s2rho
              t1333 = t437 * s
              Fxrho1rho = -0.9950248765e1_dp * t437 * t438 * s2rho + &
                   0.2e1_dp * t186 * s2rho * srho + 0.2e1_dp * t186 *&
                   s * srho1rho
              Fx1rhorho = -0.9950248765e1_dp * t437 * t91 * s1rho * s2rho + &
                   0.2e1_dp * t186 * s2rho * s1rho + 0.2e1_dp * t186 &
                   * s * s1rhorho
              Fx2rho = 0.2e1_dp * t186 * s * s2rho
              ex_ldarhorhorho = ex_unif1rhorho * Fx + ex_unif1rho * Fx2rho + &
                   ex_unif2rho * Fx1rho + ex_unif * Fx1rhorho + ex_unifrho1rho * Fx + &
                   ex_unifrho * Fx2rho + ex_unifrhorho * Fx - 0.3e1_dp / 0.4e1_dp * my_rho&
                   * t5 * kfrhorhorho * Fx + t453 * Fx2rho + ex_unifrho * Fx1rho + my_rho&
                   * ex_unifrho1rho * Fx1rho + t190 * Fx1rhorho + ex_unif2rho * Fxrho &
                   + ex_unif * Fxrho1rho + ex_unif1rho * Fxrho + my_rho * ex_unif1rhorho *&
                   Fxrho + t457 * Fxrho1rho + ex_unif * Fxrhorho + my_rho * ex_unif2rho *&
                   Fxrhorho + t97 * (0.7425558783e2_dp * t1329 * t1330 - &
                   0.2985074630e2_dp * t1333 * t1330 - 0.9950248765e1_dp * t437 * t91 *&
                   srho1rho * s1rho - 0.9950248765e1_dp * t437 * t438 * s1rhorho + &
                   0.2e1_dp * t186 * s1rhorho * srho + 0.2e1_dp * &
                   t186 * s1rho * srho1rho - 0.9950248765e1_dp * t437 * t91 * srhorho *&
                   s2rho + 0.2e1_dp * t186 * s2rho * srhorho + &
                   0.2e1_dp * t186 * s * (-0.3e1_dp * my_norm_drho / t1286 * t426&
                   * kf2rho - t424 * t102 * t425 + 0.2e1_dp * t424 * t430 * kfrho - &
                   0.2e1_dp * t424 * t428 * kf2rho - 0.2e1_dp * t179 * t224 * kfrho + &
                   0.3e1_dp / 0.2e1_dp * t179 * t102 * kfrhorho + t424 * t430 * kf2rho &
                   - t179 * t6 * kfrhorhorho / 0.2e1_dp - t179 * t224 * kf2rho - &
                   0.3e1_dp * t89 * t219))

              e_rho_rho_rho(ii) = e_rho_rho_rho(ii)&
                   +ex_ldarhorhorho + e_c_u_01rhorho + gamma_var * &
                   t1273 * t173 - t418 * t1208 + e_c_u_0rho1rho + t1237 * t173 - t172 *&
                   t1208 + epsilon_cGGArhorho + my_rho * (e_c_u_0rhorhorho + gamma_var * &
                   (t951 + t1019 + t1061 + t1194) * t173 - t403 * t1208 - t1237 * t416 &
                   + 0.2e1_dp * t172 * t1240 * t415 * t1207 - t172 * t406 * t1273)

              t1406 = t132 * t102
              tnorm_drhorhorho = t281 * t6 * t754 + t1406 * k_srho / 0.2e1_dp - &
                   t460 * k_srhorho / 0.2e1_dp + t1406 * k_s1rho / 0.2e1_dp + t63 * &
                   t224
              tnorm_drho1rho = -t460 * k_s1rho / 0.2e1_dp - t462 / 0.2e1_dp
              t1425 = tnorm_drho * A1rho
              t1433 = 0.2e1_dp * t349 * tnorm_drho + 0.2e1_dp * t352 * &
                   tnorm_drho + 0.2e1_dp * t150 * tnorm_drho1rho + 0.8e1_dp * t383 * &
                   t1425 + 0.12e2_dp * t392 * tnorm_drho * t1rho + 0.4e1_dp * t164 * &
                   tnorm_drho1rho
              t1438 = t197 * Arho
              t1439 = t470 * t335
              t1442 = tnorm_drhorho * t158
              t1443 = t1442 * t335
              t1446 = Arho * tnorm_drho1rho
              t1451 = t167 * t335
              t1455 = t464 * t75
              t1458 = t375 * tnorm_drho
              t1459 = t1458 * t1451
              t1477 = 0.133450e0_dp * t66 * trhorho * t194 + 0.133450e0_dp * &
                   t373 * t375 * t1433 * t167 - 0.133450e0_dp * t1438 * t1439 - &
                   0.133450e0_dp * t490 * t1443 + 0.133450e0_dp * t197 * t1446 * t80 - &
                   0.400350e0_dp * t373 * t1167 * t206 * t1451 - 0.133450e0_dp * t1455 &
                   * t1439 + 0.266900e0_dp * t329 * t1459 - 0.66725e-1_dp * t73 * t369 &
                   * t518 + 0.133450e0_dp * t373 * t500 * t398 + 0.133450e0_dp * t145 *&
                   t325 * tnorm_drhorho + 0.133450e0_dp * t318 * t474 - 0.133450e0_dp &
                   * t329 * t158 * t518 * t1rho
              t1480 = t500 * t335
              t1484 = t478 * t1rho
              t1497 = trho * t335
              t1501 = t81 * tnorm_drho1rho
              t1506 = t145 * A
              t1514 = t158 * tnorm_drho1rho * t167
              t1517 = t470 * t398
              t1520 = -0.133450e0_dp * t329 * t1443 + 0.133450e0_dp * t898 * &
                   t1480 + 0.400350e0_dp * t73 * Arho * t1484 - 0.66725e-1_dp * t73 * &
                   t361 * t1433 + 0.400350e0_dp * t477 * tnorm_drhorho * t80 * t1rho + &
                   0.133450e0_dp * t145 * t358 * tnorm_drho - 0.400350e0_dp * t477 * &
                   t470 * t1497 + 0.133450e0_dp * t464 * t1501 - 0.133450e0_dp * t1068 &
                   * t494 + 0.800700e0_dp * t1506 * t478 * t393 + 0.400350e0_dp * t73 *&
                   A1rho * t479 - 0.133450e0_dp * t329 * t1514 - 0.133450e0_dp * t490 &
                   * t1517
              t1533 = t375 * t518
              t1552 = Arhorho * tnorm_drho
              t1563 = 0.133450e0_dp * t197 * A * tnorm_drhorhorho * t80 - &
                   0.133450e0_dp * t197 * A1rho * t471 + 0.400350e0_dp * t477 * &
                   tnorm_drho1rho * t80 * trho + 0.133450e0_dp * t373 * t1533 * t335 - &
                   0.133450e0_dp * t490 * t1514 - 0.133450e0_dp * t1004 * t1439 + &
                   0.133450e0_dp * t145 * t81 * tnorm_drhorhorho + 0.266900e0_dp * t490&
                   * t1459 + 0.133450e0_dp * t145 * t154 * tnorm_drho1rho + &
                   0.400350e0_dp * t477 * t478 * trhorho + 0.133450e0_dp * t197 * t1552&
                   * t80 - 0.133450e0_dp * t329 * t158 * t1433 * trho - 0.133450e0_dp &
                   * t329 * t493 * trhorho
              t1564 = t167 * t1rho
              t1571 = t325 * tnorm_drho
              t1578 = A1rho * tnorm_drhorho
              t1584 = t384 * tnorm_drho
              t1588 = t389 * tnorm_drho
              t1627 = 0.2e1_dp * t1074 * tnorm_drho + 0.2e1_dp * t1584 + &
                   0.2e1_dp * t346 * tnorm_drho1rho + 0.2e1_dp * t1588 + 0.2e1_dp * A *&
                   trhorho * tnorm_drho + 0.2e1_dp * t506 * tnorm_drho1rho + 0.2e1_dp &
                   * t349 * tnorm_drhorho + 0.2e1_dp * t352 * tnorm_drhorho + 0.2e1_dp &
                   * t150 * tnorm_drhorhorho + 0.8e1_dp * t1133 * t482 + 0.24e2_dp * &
                   t74 * t1584 + 0.8e1_dp * t383 * t1446 + 0.8e1_dp * t383 * t1552 + &
                   0.24e2_dp * t74 * t1588 + 0.24e2_dp * t1105 * t513 * t1rho + &
                   0.12e2_dp * t392 * tnorm_drho1rho * trho + 0.12e2_dp * t392 * &
                   tnorm_drho * trhorho + 0.8e1_dp * t383 * t1578 + 0.12e2_dp * t392 * &
                   tnorm_drhorho * t1rho + 0.4e1_dp * t164 * tnorm_drhorhorho
              t1636 = t493 * t1rho
              t1648 = -0.400350e0_dp * t477 * t470 * t1564 + 0.266900e0_dp * &
                   t329 * t500 * t1497 + 0.133450e0_dp * t464 * t1571 - 0.133450e0_dp *&
                   t1053 * t471 - 0.133450e0_dp * t1053 * t494 + 0.133450e0_dp * t197 &
                   * t1578 * t80 - 0.66725e-1_dp * t73 * t159 * t1627 - 0.66725e-1_dp *&
                   t73 * t992 * t206 + 0.133450e0_dp * t1049 * t501 - 0.133450e0_dp * &
                   t1004 * t1636 - 0.133450e0_dp * t1068 * t471 + 0.133450e0_dp * t318 &
                   * t467 + 0.266900e0_dp * t329 * t500 * t1564 - 0.133450e0_dp * t329 &
                   * t1517
              t1683 = 0.133450e0_dp * t318 * t194 + 0.133450e0_dp * t145 * t1571&
                   - 0.133450e0_dp * t329 * t1439 + 0.133450e0_dp * t145 * t1501 + &
                   0.400350e0_dp * t477 * t1484 + 0.133450e0_dp * t197 * t1425 * t80 + &
                   0.133450e0_dp * t197 * A * tnorm_drho1rho * t80 - 0.133450e0_dp * &
                   t490 * t1439 - 0.133450e0_dp * t329 * t1636 - 0.66725e-1_dp * t73 * &
                   t369 * t206 + 0.133450e0_dp * t373 * t1480 - 0.66725e-1_dp * t73 * &
                   t159 * t1433
              t1684 = gamma_var * t1683
              t1686 = t1240 * t171
              snorm_drho1rho = snorm_drhorho
              t1702 = t536 * s1rho
              t1720 = t91 * snorm_drhorho
              Fxnorm_drho1rho = -0.9950248765e1_dp * t437 * t532 * s1rho + &
                   0.2e1_dp * t186 * s1rho * snorm_drho + 0.2e1_dp * &
                   t186 * s * snorm_drho1rho

              e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii)&
                   +ex_unif1rho * Fxnorm_drho + ex_unif * &
                   Fxnorm_drho1rho + ex_unifrho * Fxnorm_drho + t453 * Fxnorm_drho + &
                   t190 * Fxnorm_drho1rho + ex_unif * Fxnorm_drhorho + t457 * &
                   Fxnorm_drhorho + t97 * (0.7425558783e2_dp * t1329 * t1702 - &
                   0.2985074630e2_dp * t1333 * t1702 - 0.9950248765e1_dp * t437 * t91 *&
                   snorm_drho1rho * srho - 0.9950248765e1_dp * t437 * t532 * srhorho +&
                   0.2e1_dp * t186 * srhorho * snorm_drho + 0.2e1_dp&
                   * t186 * srho * snorm_drho1rho - 0.9950248765e1_dp * t437 * t1720 *&
                   s1rho + 0.2e1_dp * t186 * s1rho * snorm_drhorho + &
                   0.2e1_dp * t186 * s * (t423 * t6 * t425 + t178 * t102 * &
                   kfrho - t528 * kfrhorho / 0.2e1_dp + t88 * t224)) + t1684 * t173 - &
                   t211 * t416 + Hnorm_drhorho + my_rho * (gamma_var * (t1477 + t1520 + &
                   t1563 + t1648) * t173 - t523 * t416 - t1684 * t525 + 0.2e1_dp * t211&
                   * t1686 * t415 - t211 * t406 * t402)

              t1747 = t470 * t518
              t1759 = t375 * t559
              t1766 = Arho * t545
              t1786 = t1442 * t206
              t1794 = -0.266900e0_dp * t490 * t1747 - 0.800700e0_dp * t477 * &
                   t470 * t206 * trho + 0.1334500e1_dp * t477 * t478 * tnorm_drhorho - &
                   0.266900e0_dp * t1004 * t553 + 0.266900e0_dp * t329 * t1759 * trho +&
                   0.266900e0_dp * t373 * t1533 * t206 - 0.66725e-1_dp * t73 * t159 * &
                   (0.2e1_dp * t1766 + 0.4e1_dp * t486 * tnorm_drho + 0.24e2_dp * t74 *&
                   t1766 + 0.24e2_dp * t1105 * t545 * trho + 0.24e2_dp * t392 * &
                   tnorm_drhorho * tnorm_drho) - 0.266900e0_dp * t1455 * t553 - &
                   0.133450e0_dp * t546 * t168 - 0.266900e0_dp * t490 * t1786 - &
                   0.266900e0_dp * t329 * t1786 + 0.266900e0_dp * t66 * tnorm_drhorho *&
                   t194
              t1806 = t500 * tnorm_drho * t167
              t1811 = t545 * t80
              t1815 = t375 * t566
              t1829 = t545 * t158
              t1833 = t1167 * t559
              t1837 = 0.133450e0_dp * t73 * t153 * t375 * t559 - 0.133450e0_dp *&
                   t329 * t158 * t566 * trho - 0.266900e0_dp * t1438 * t553 + &
                   0.533800e0_dp * t329 * t1806 - 0.266900e0_dp * t329 * t1747 + &
                   0.1334500e1_dp * t1506 * t1811 * trho + 0.133450e0_dp * t373 * t1815&
                   * t167 + 0.133450e0_dp * t546 * t154 + 0.533800e0_dp * t490 * t1806&
                   - 0.66725e-1_dp * t73 * t361 * t566 + 0.667250e0_dp * t73 * t1766 *&
                   t80 - 0.667250e0_dp * t477 * t1829 * t167 - 0.400350e0_dp * t373 * &
                   t1833 * t167
              t1842 = t406 * t210
              t1848 = t1328 * t576
              t1852 = s * t576

              e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii)&
                   +ex_unif * Fxnorm_drhonorm_drho + &
                   t190 * Fxnorm_drhonorm_drho + t97 * (0.7425558783e2_dp * t1327 * &
                   t1848 * srho - 0.2985074630e2_dp * t437 * t1852 * srho - &
                   0.1990049753e2_dp * t437 * t1720 * snorm_drho + 0.4000000002e1_dp * &
                   t186 * snorm_drhorho * snorm_drho) + Hnorm_drhonorm_drho + my_rho * (&
                   gamma_var * (t1794 + t1837) * t173 - t571 * t525 - 0.2e1_dp * t523 *&
                   t1842 + 0.2e1_dp * t574 * t1686)

              t2norm_drho = tnorm_drho
              t1865 = A * t2norm_drho
              t1871 = 0.2e1_dp * t1865 * tnorm_drho + 0.12e2_dp * t392 * &
                   tnorm_drho * t2norm_drho
              t1876 = t2norm_drho * t158 * t566
              t1883 = 0.2e1_dp * t150 * t2norm_drho + 0.4e1_dp * t164 * t2norm_drho
              t1885 = t1458 * t206 * t1883
              t1891 = t158 * t1871 * tnorm_drho
              t1910 = t1759 * t2norm_drho
              t1913 = t66 * t2norm_drho
              t1924 = t159 * t1883
              t1932 = 0.266900e0_dp * t373 * t500 * t1871 - 0.133450e0_dp * t490&
                   * t1876 + 0.533800e0_dp * t490 * t1885 + 0.533800e0_dp * t329 * &
                   t1885 - 0.266900e0_dp * t329 * t1891 - 0.266900e0_dp * t490 * t1891 &
                   + 0.133450e0_dp * t373 * t1815 * t1883 - 0.1601400e1_dp * t197 * t75&
                   * t158 * t76 * t545 * t2norm_drho - 0.133450e0_dp * t329 * t1876 + &
                   0.1601400e1_dp * t1506 * t1811 * t2norm_drho + 0.266900e0_dp * t329 &
                   * t1910 - 0.266900e0_dp * t1913 * t75 * t553 - 0.1334500e1_dp * t477&
                   * t470 * t206 * t2norm_drho - 0.667250e0_dp * t477 * t1829 * t1883 &
                   - 0.133450e0_dp * t546 * t1924 + 0.266900e0_dp * t490 * t1910 - &
                   0.400350e0_dp * t373 * t1833 * t1883
              t1943 = 0.133450e0_dp * t145 * t81 * t2norm_drho + 0.133450e0_dp *&
                   t197 * t1865 * t80 - 0.66725e-1_dp * t73 * t1924
              t1951 = t470 * t1883
              t1956 = t493 * t2norm_drho
              s2norm_drho = snorm_drho

              e_ndrho_ndrho_ndrho(ii) = e_ndrho_ndrho_ndrho(ii)&
                   +t97 * (0.7425558783e2_dp * &
                   t1327 * t1848 * s2norm_drho - 0.2985074630e2_dp * t437 * t1852 * &
                   s2norm_drho) + my_rho * (gamma_var * t1932 * t173 - t571 * t406 * t1943&
                   - 0.2e1_dp * gamma_var * (0.133450e0_dp * t1913 * t194 + &
                   0.667250e0_dp * t477 * t478 * t2norm_drho - 0.133450e0_dp * t329 * &
                   t1951 - 0.133450e0_dp * t490 * t1951 - 0.133450e0_dp * t329 * t1956 &
                   - 0.133450e0_dp * t490 * t1956 + 0.133450e0_dp * t373 * t500 * t1883&
                   - 0.66725e-1_dp * t73 * t159 * t1871) * t1842 + 0.2e1_dp * t574 * &
                   t1240 * t1943)
           END IF

        END IF
     END DO
  END SELECT
END SUBROUTINE pbe_lda_calc
!***************************************************************************

!!****f* xc_pbe/pbe_lsd_eval [1.0] *
!!
!!   NAME
!!     pbe_lsd_eval
!!
!!   SYNOPSIS
!!     Subroutine pbe_lsd_eval(rho_set, deriv_set, grad_deriv, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: grad_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pbe_lsd_eval
!!
!!   FUNCTION
!!     evaluates the becke 88 exchange functional for lsd
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pbe_lsd_eval(rho_set,deriv_set,grad_deriv,pbe_param,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), pointer         :: pbe_param
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndr, &
      e_ndr_ndr, e_ndr_ra, e_ndr_rb, e_ndra, e_ndra_ndra, e_ndra_ra, &
      e_ndrb, e_ndrb_ndrb, e_ndrb_rb, e_ra, e_ra_ra, &
      e_ra_rb, e_rb, e_rb_rb, norm_drho, norm_drhoa, norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  NULLIFY(deriv, bo)

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,&
          rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
          norm_drhob=norm_drhob, norm_drho=norm_drho, &
          rho_cutoff=epsilon_rho,&
          drho_cutoff=epsilon_drho, local_bounds=bo, error=error)
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
     e_ra => dummy
     e_rb => dummy
     e_ndra_ra => dummy
     e_ndrb_rb => dummy
     e_ndr_ndr => dummy
     e_ndra_ndra => dummy
     e_ndrb_ndrb => dummy
     e_ndr => dummy
     e_ndra => dummy
     e_ndrb => dummy
     e_ra_ra => dummy
     e_ra_rb => dummy
     e_rb_rb => dummy
     e_ndr_ra => dummy
     e_ndr_rb => dummy

     IF (grad_deriv>=0) THEN
        deriv => xc_dset_get_derivative(deriv_set,"",&
             allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)
     END IF
     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrb,error=error)
     END IF
     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ra_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ra_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rb_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndra_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrb_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr_ndr,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndra_ndra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrb_ndrb,error=error)
     END IF
     IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
        ! to do
     END IF
     CALL pbe_lsd_calc(&
          rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, norm_drhoa=norm_drhoa,&
          norm_drhob=norm_drhob, e_0=e_0, e_ra=e_ra, e_rb=e_rb,&
          e_ra_ndra=e_ndra_ra,&
          e_rb_ndrb=e_ndrb_rb, e_ndr_ndr=e_ndr_ndr,&
          e_ndra_ndra=e_ndra_ndra, e_ndrb_ndrb=e_ndrb_ndrb, e_ndr=e_ndr,&
          e_ndra=e_ndra, e_ndrb=e_ndrb, e_ra_ra=e_ra_ra, &
          e_ra_rb=e_ra_rb, e_rb_rb=e_rb_rb, e_ra_ndr=e_ndr_ra,&
          e_rb_ndr=e_ndr_rb,&
          grad_deriv=grad_deriv, npoints=npoints, &
          epsilon_rho=epsilon_rho,epsilon_drho=epsilon_drho,&
          pbe_param=pbe_param,error=error)

     IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     ELSE
        NULLIFY(dummy)
     END IF
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE pbe_lsd_eval
!***************************************************************************
 
!!****f* xc_pbe/pbe_lsd_calc [1.0] *
!!
!!   NAME
!!     pbe_lsd_calc
!!
!!   FUNCTION
!!     low level calculation of the pbe exchange-correlation functional for lsd
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rhoa,rhob: alpha or beta spin density
!!     - norm_drhoa,norm_drhob,norm_drho: || grad rhoa |||,| grad rhoa ||,
!!       || grad (rhoa+rhob) ||
!!     - e_0: adds to it the local value of the functional
!!     - e_*: derivative of the functional wrt. to the variables
!!       named where the * is.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pbe_lsd_calc(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
     e_0, e_ra, e_rb, e_ra_ndra, e_rb_ndrb, e_ndr_ndr,&
     e_ndra_ndra, e_ndrb_ndrb, e_ndr,&
     e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr,&
     grad_deriv,npoints,epsilon_rho,epsilon_drho,pbe_param,error)
  REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drho, &
       norm_drhoa, norm_drhob
  REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_ra, e_rb, e_ra_ndra, &
       e_rb_ndrb, e_ndr_ndr, e_ndra_ndra, e_ndrb_ndrb,&
       e_ndr, e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr
  INTEGER, INTENT(in)                      :: grad_deriv, npoints
  REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
  TYPE(section_vals_type), pointer :: pbe_param
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_calc', &
       routineP = moduleN//':'//routineN
  REAL(kind=dp), PARAMETER :: small = 1.0e-20, t_1_9 = 1.0_dp/9.0_dp, &
       t_2_3 = 2.0_dp/3.0_dp, t_4_27 = 4.0_dp/27.0_dp, t_4_3 = 4.0_dp/3.0_dp, &
       t_4_9 = 4.0_dp/9.0_dp, t_8_27 = 8.0_dp/27.0_dp

  INTEGER                                  :: ii
  LOGICAL                                  :: failure

  real(kind=dp) :: my_rhoa, my_rhob, my_norm_drho, my_norm_drhoa, &
       my_norm_drhob, kappa, beta, mu, gamma_var, p_1, A_1, alpha_1_1, &
       beta_1_1, beta_2_1, beta_3_1, beta_4_1, p_2, A_2, alpha_1_2, &
       beta_1_2, beta_2_2, beta_3_2, beta_4_2, p_3, A_3, alpha_1_3, &
       beta_1_3, beta_2_3, beta_3_3, beta_4_3, scale_ec, scale_ex, my_rho, t1,&
       t2, chi, t3, t4, t5, t6, t7, t8, t9, rs, t12, t14, t15, t18, t20, &
       t21, t22, t23, t27, t28, e_c_u_0, t32, t34, t38, t39, t40, t41, t45,&
       t46, t50, t52, t56, t57, t58, t59, t63, t64, alpha_c, t66, t69, &
       t70, t71, t72, t73, t74, t75, f, t77, t78, t79, t80, t82, t84, t85, &
       epsilon_c_unif, t87, t88, phi, t90, t91, t92, t93, t94, k_s, t95, &
       t96, t97, t98, t, t100, t101, t102, t103, t104, t105, t107, t108, A,&
       t110, t111, t112, t113, t115, t116, t118, t119, t122, t123
  real(kind=dp) :: epsilon_cGGA, t124, t125, t126, kf_a, ex_unif_a, t129, t130, t131, &
       s_a, t133, t135, t137, Fx_a, t140, t142, t143, kf_b, ex_unif_b, &
       t146, t147, t148, s_b, t150, t153, Fx_b, t156, exc, t162, t163, &
       t164, chirhoa, t165, t167, rsrhoa, t171, t175, t176, t177, t178, &
       t179, t183, t187, t190, t191, t192, e_c_u_0rhoa, t194, t198, t199, &
       t200, t201, t205, t211, t212, t213, e_c_u_1rhoa, t215, t219, t220, &
       t221, t222, t226, t232, t233, t234, alpha_crhoa, frhoa, t240, t242, &
       t244, t245, t246, t248, t249, t250, t252, t254, t256
  real(kind=dp) :: epsilon_c_unifrhoa, t257, t259, phirhoa, t262, k_frhoa, t266, &
       k_srhoa, t268, t269, t272, t273, t274, t277, t278, trhoa, t280, &
       t281, t282, t284, t285, t286, t289, Arhoa, t293, t294, t297, t298, &
       t299, t302, t303, t305, t306, t310, t311, t312, t313, t314, t317, &
       t318, t321, t324, t325, t326, epsilon_cGGArhoa, t329, kf_arhoa, &
       ex_unif_arhoa, t335, t336, t337, t340, t341, s_arhoa, t344, t346, &
       Fx_arhoa, t350, exc_rhoa, chirhob, rsrhob, t368, e_c_u_0rhob, t382, &
       e_c_u_1rhob, t396, alpha_crhob, frhob, t403, t405, t407, t409, t410,&
       t411, t413, t415, t417, epsilon_c_unifrhob, phirhob, k_frhob, &
       k_srhob, trhob, t427, t429, t432, Arhob, t436, t439, t442, t444
  real(kind=dp) :: t445, t453, t456, t457, epsilon_cGGArhob, t460, kf_brhob, &
       ex_unif_brhob, t466, t467, t468, t471, t472, s_brhob, t475, t477, &
       Fx_brhob, t481, exc_rhob, t488, tnorm_drho, t493, t494, t502, t505, &
       t506, Hnorm_drho, exc_norm_drho, s_anorm_drhoa, Fx_anorm_drhoa, &
       exc_norm_drhoa, s_bnorm_drhob, Fx_bnorm_drhob, exc_norm_drhob, t518,&
       t519, chirhoarhoa, t523, t525, rsrhoarhoa, t536, t538, t543, t544, &
       t548, t549, t550, t556, t561, t563, t564, t576, t578, t579, t580, &
       e_c_u_0rhoarhoa, e_c_u_01rhoa, t588, t590, t595, t596, t600, t606
  real(kind=dp) :: t611, t624, t626, t627, t628, t636, t638, t643, t644, t648, t654, &
       t659, t672, t674, t675, t676, alpha_c1rhoa, t681, t682, t687, &
       frhoarhoa, f1rhoa, t705, t708, t711, t726, t733, t736, t745, &
       epsilon_c_unifrhoarhoa, epsilon_c_unif1rhoa, t750, t755, &
       phirhoarhoa, phi1rhoa, t763, k_frhoarhoa, t767, t768, k_s1rhoa, &
       t775, t776, t779, t785, t789, t795, t798, t801, t812, trhoarhoa, &
       t1rhoa, t820, t821, t822, t823, t825, t828, t839, t840, t851, &
       Arhoarhoa, A1rhoa, t858, t865, t867, t868, t876, t879, t880, t904, &
       t908, t909, t911, t914, t917, t919, t924, t936, t944, t953, t959, &
       t962, t965, t966, t967, kf_arhoarhoa, ex_unif_a1rhoa, t985, &
       s_a1rhoa, t998, t999, t1000, t1001, Fx_a1rhoa, exc_rhoa_rhoa
  real(kind=dp) :: chirhoarhob, rsrhoarhob, t1031, t1033, t1038, t1050, &
       e_c_u_0rhoarhob, t1069, t1071, t1104, t1106, frhoarhob, t1164, &
       t1193, epsilon_c_unifrhoarhob, phirhoarhob, k_frhoarhob, t1228, &
       t1231, trhoarhob, Arhoarhob, t1269, t1283, t1285, t1286, t1288, &
       t1291, t1293, t1304, t1327, t1330, exc_rhoa_rhob, chirhobrhob, &
       rsrhobrhob, t1342, t1346, e_c_u_0rhobrhob, e_c_u_01rhob, t1377, &
       t1411, alpha_c1rhob, t1440, frhobrhob, f1rhob, t1462, t1465, t1482, &
       t1489, t1492, t1501, epsilon_c_unifrhobrhob, epsilon_c_unif1rhob
  real(kind=dp) :: phirhobrhob, phi1rhob, t1514, k_s1rhob, t1520, t1529, trhobrhob, &
       t1rhob, t1550, t1552, t1555, Arhobrhob, A1rhob, t1588, t1590, t1591,&
       t1599, t1602, t1603, t1630, t1632, t1635, t1638, t1640, t1674, &
       t1677, t1680, kf_brhobrhob, ex_unif_b1rhob, t1698, s_b1rhob, t1711, &
       t1712, t1713, Fx_b1rhob, exc_rhob_rhob, t1739, t1741, t1743, &
       trhoanorm_drho, t1748, t1765, t1766, t1767, t1801
  real(kind=dp) :: exc_rhoa_norm_drho, trhobnorm_drho, t1829, t1830, t1831, t1865, &
       exc_rhob_norm_drho, t1871, t1876, t1888, t1901, &
       exc_norm_drho_norm_drho, exc_rhoa_norm_drhoa, t1922, &
       exc_norm_drhoa_norm_drhoa, exc_rhob_norm_drhob, t1949, &
       exc_norm_drhob_norm_drhob, r_eqs_lsd4

  failure=.FALSE.

!FM  CALL pbe_lsd_calc_old(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
!FM       e_0, e_ra, e_rb, e_ra_ndra, e_rb_ndrb, e_ndr_ndr,&
!FM       e_ndra_ndra, e_ndrb_ndrb, e_ndr,&
!FM       e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr,&
!FM       grad_deriv,npoints,epsilon_rho,epsilon_drho,error)

  kappa = 0.804e0_dp
  beta = 0.66725e-1_dp
  mu = beta * 0.3141592654e1_dp ** 2 / 0.3e1_dp
  gamma_var = (0.1e1_dp - log(0.2e1_dp)) / 0.3141592654e1_dp ** 2
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
  A_3 = 0.16887e-1_dp ! correct
  !A_3 = 0.16887e0_dp ! wrong
  alpha_1_3 = 0.11125e0_dp
  beta_1_3 = 0.10357e2_dp
  beta_2_3 = 0.36231e1_dp
  beta_3_3 = 0.88026e0_dp
  beta_4_3 = 0.49671e0_dp

  scale_ec = 1._dp ! -0.1e1_dp
  scale_ex = 1._dp ! -0.1e1_dp
  SELECT CASE(grad_deriv)
  CASE default
     DO ii=1,npoints
        my_rhoa=MAX(rhoa(ii),0.0_dp)
        my_rhob=MAX(rhob(ii),0.0_dp)
        my_rho=my_rhoa+my_rhob
        IF (my_rho>epsilon_rho) THEN
           my_rhoa=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhoa)
           my_rhob=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhob)
           my_rho = my_rhoa + my_rhob
           my_norm_drho=norm_drho(ii)
           my_norm_drhoa=norm_drhoa(ii)
           my_norm_drhob=norm_drhob(ii)


           t1 = my_rhoa - my_rhob
           t2 = 0.1e1_dp / my_rho
           chi = t1 * t2
           t3 = 3 ** (0.1e1_dp / 0.3e1_dp)
           t4 = 4 ** (0.1e1_dp / 0.3e1_dp)
           t5 = t4 ** 2
           t6 = t3 * t5
           t7 = 0.1e1_dp / 0.3141592654e1_dp
           t8 = t7 * t2
           t9 = t8 ** (0.1e1_dp / 0.3e1_dp)
           rs = t6 * t9 / 0.4e1_dp
           t12 = 0.1e1_dp + alpha_1_1 * rs
           t14 = 0.1e1_dp / A_1
           t15 = sqrt(rs)
           t18 = t15 * rs
           t20 = p_1 + 0.1e1_dp
           t21 = rs ** t20
           t22 = beta_4_1 * t21
           t23 = beta_1_1 * t15 + beta_2_1 * rs + beta_3_1 * t18 + t22
           t27 = 0.1e1_dp + t14 / t23 / 0.2e1_dp
           t28 = log(t27)
           e_c_u_0 = -0.2e1_dp * A_1 * t12 * t28
           t32 = 0.1e1_dp + alpha_1_2 * rs
           t34 = 0.1e1_dp / A_2
           t38 = p_2 + 0.1e1_dp
           t39 = rs ** t38
           t40 = beta_4_2 * t39
           t41 = beta_1_2 * t15 + beta_2_2 * rs + beta_3_2 * t18 + t40
           t45 = 0.1e1_dp + t34 / t41 / 0.2e1_dp
           t46 = log(t45)
           t50 = 0.1e1_dp + alpha_1_3 * rs
           t52 = 0.1e1_dp / A_3
           t56 = p_3 + 0.1e1_dp
           t57 = rs ** t56
           t58 = beta_4_3 * t57
           t59 = beta_1_3 * t15 + beta_2_3 * rs + beta_3_3 * t18 + t58
           t63 = 0.1e1_dp + t52 / t59 / 0.2e1_dp
           t64 = log(t63)
           alpha_c = 0.2e1_dp * A_3 * t50 * t64
           t66 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t69 = 1 / (2 * t66 - 2)
           t70 = 0.1e1_dp + chi
           t71 = t70 ** (0.1e1_dp / 0.3e1_dp)
           t72 = t71 * t70
           t73 = 0.1e1_dp - chi
           t74 = t73 ** (0.1e1_dp / 0.3e1_dp)
           t75 = t74 * t73
           f = (t72 + t75 - 0.2e1_dp) * t69
           t77 = alpha_c * f
           t78 = 0.9e1_dp / 0.8e1_dp / t69
           t79 = chi ** 2
           t80 = t79 ** 2
           t82 = t78 * (0.1e1_dp - t80)
           t84 = -0.2e1_dp * A_2 * t32 * t46 - e_c_u_0
           t85 = t84 * f
           epsilon_c_unif = e_c_u_0 + t77 * t82 + t85 * t80
           t87 = t71 ** 2
           t88 = t74 ** 2
           phi = t87 / 0.2e1_dp + t88 / 0.2e1_dp
           t90 = 0.3141592654e1_dp ** 2
           t91 = t90 * my_rho
           t92 = t91 ** (0.1e1_dp / 0.3e1_dp)
           t93 = t3 * t92 * t7
           t94 = sqrt(t93)
           k_s = 0.2e1_dp * t94
           t95 = 0.1e1_dp / phi
           t96 = my_norm_drho * t95
           t97 = 0.1e1_dp / k_s
           t98 = t97 * t2
           t = t96 * t98 / 0.2e1_dp
           t100 = 0.1e1_dp / gamma_var
           t101 = beta * t100
           t102 = epsilon_c_unif * t100
           t103 = phi ** 2
           t104 = t103 * phi
           t105 = 0.1e1_dp / t104
           t107 = exp(-t102 * t105)
           t108 = t107 - 0.1e1_dp
           A = t101 / t108
           t110 = gamma_var * t104
           t111 = t ** 2
           t112 = A * t111
           t113 = 0.1e1_dp + t112
           t115 = A ** 2
           t116 = t111 ** 2
           t118 = 0.1e1_dp + t112 + t115 * t116
           t119 = 0.1e1_dp / t118
           t122 = 0.1e1_dp + t101 * t111 * t113 * t119
           t123 = log(t122)
           epsilon_cGGA = epsilon_c_unif + t110 * t123
           t124 = t3 * t66
           t125 = t90 * my_rhoa
           t126 = t125 ** (0.1e1_dp / 0.3e1_dp)
           kf_a = t124 * t126
           ex_unif_a = -0.3e1_dp / 0.4e1_dp * t7 * kf_a
           t129 = 0.1e1_dp / kf_a
           t130 = my_norm_drhoa * t129
           t131 = 0.1e1_dp / my_rhoa
           s_a = t130 * t131 / 0.2e1_dp
           t133 = s_a ** 2
           t135 = 0.1e1_dp / kappa
           t137 = 0.1e1_dp + mu * t133 * t135
           Fx_a = 0.1e1_dp + kappa - kappa / t137
           t140 = my_rhoa * ex_unif_a
           t142 = t90 * my_rhob
           t143 = t142 ** (0.1e1_dp / 0.3e1_dp)
           kf_b = t124 * t143
           ex_unif_b = -0.3e1_dp / 0.4e1_dp * t7 * kf_b
           t146 = 0.1e1_dp / kf_b
           t147 = my_norm_drhob * t146
           t148 = 0.1e1_dp / my_rhob
           s_b = t147 * t148 / 0.2e1_dp
           t150 = s_b ** 2
           t153 = 0.1e1_dp + mu * t150 * t135
           Fx_b = 0.1e1_dp + kappa - kappa / t153
           t156 = my_rhob * ex_unif_b

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)+&
                   scale_ex * (0.2e1_dp * t140 * Fx_a + 0.2e1_dp * t156 * Fx_b)&
                   / 0.2e1_dp + scale_ec * my_rho * epsilon_cGGA
           END IF

           t162 = my_rho ** 2
           t163 = 0.1e1_dp / t162
           t164 = t1 * t163
           chirhoa = t2 - t164
           t165 = t9 ** 2
           t167 = 0.1e1_dp / t165 * t7
           rsrhoa = -t6 * t167 * t163 / 0.12e2_dp
           t171 = A_1 * alpha_1_1
           t175 = t23 ** 2
           t176 = 0.1e1_dp / t175
           t177 = t12 * t176
           t178 = 0.1e1_dp / t15
           t179 = beta_1_1 * t178
           t183 = beta_3_1 * t15
           t187 = 0.1e1_dp / rs
           t190 = t179 * rsrhoa / 0.2e1_dp + beta_2_1 * rsrhoa + 0.3e1_dp / &
                0.2e1_dp * t183 * rsrhoa + t22 * t20 * rsrhoa * t187
           t191 = 0.1e1_dp / t27
           t192 = t190 * t191
           e_c_u_0rhoa = -0.2e1_dp * t171 * rsrhoa * t28 + t177 * t192
           t194 = A_2 * alpha_1_2
           t198 = t41 ** 2
           t199 = 0.1e1_dp / t198
           t200 = t32 * t199
           t201 = beta_1_2 * t178
           t205 = beta_3_2 * t15
           t211 = t201 * rsrhoa / 0.2e1_dp + beta_2_2 * rsrhoa + 0.3e1_dp / &
                0.2e1_dp * t205 * rsrhoa + t40 * t38 * rsrhoa * t187
           t212 = 0.1e1_dp / t45
           t213 = t211 * t212
           e_c_u_1rhoa = -0.2e1_dp * t194 * rsrhoa * t46 + t200 * t213
           t215 = A_3 * alpha_1_3
           t219 = t59 ** 2
           t220 = 0.1e1_dp / t219
           t221 = t50 * t220
           t222 = beta_1_3 * t178
           t226 = beta_3_3 * t15
           t232 = t222 * rsrhoa / 0.2e1_dp + beta_2_3 * rsrhoa + 0.3e1_dp / &
                0.2e1_dp * t226 * rsrhoa + t58 * t56 * rsrhoa * t187
           t233 = 0.1e1_dp / t63
           t234 = t232 * t233
           alpha_crhoa = 0.2e1_dp * t215 * rsrhoa * t64 - t221 * t234
           frhoa = (0.4e1_dp / 0.3e1_dp * t71 * chirhoa - 0.4e1_dp / 0.3e1_dp&
                * t74 * chirhoa) * t69
           t240 = alpha_crhoa * f
           t242 = alpha_c * frhoa
           t244 = t79 * chi
           t245 = t78 * t244
           t246 = t245 * chirhoa
           t248 = 0.4e1_dp * t77 * t246
           t249 = e_c_u_1rhoa - e_c_u_0rhoa
           t250 = t249 * f
           t252 = t84 * frhoa
           t254 = t244 * chirhoa
           t256 = 0.4e1_dp * t85 * t254
           epsilon_c_unifrhoa = e_c_u_0rhoa + t240 * t82 + t242 * t82 - t248 &
                + t250 * t80 + t252 * t80 + t256
           t257 = 0.1e1_dp / t71
           t259 = 0.1e1_dp / t74
           phirhoa = t257 * chirhoa / 0.3e1_dp - t259 * chirhoa / 0.3e1_dp
           t262 = t92 ** 2
           k_frhoa = t3 / t262 * t90 / 0.3e1_dp
           t266 = 0.1e1_dp / t94
           k_srhoa = t266 * k_frhoa * t7
           t268 = 0.1e1_dp / t103
           t269 = my_norm_drho * t268
           t272 = k_s ** 2
           t273 = 0.1e1_dp / t272
           t274 = t273 * t2
           t277 = t97 * t163
           t278 = t96 * t277
           trhoa = -t269 * t98 * phirhoa / 0.2e1_dp - t96 * t274 * k_srhoa / &
                0.2e1_dp - t278 / 0.2e1_dp
           t280 = t108 ** 2
           t281 = 0.1e1_dp / t280
           t282 = epsilon_c_unifrhoa * t100
           t284 = t103 ** 2
           t285 = 0.1e1_dp / t284
           t286 = t285 * phirhoa
           t289 = -t282 * t105 + 0.3e1_dp * t102 * t286
           Arhoa = -t101 * t281 * t289 * t107
           t293 = gamma_var * t103
           t294 = t123 * phirhoa
           t297 = t101 * t
           t298 = t113 * t119
           t299 = t298 * trhoa
           t302 = Arhoa * t111
           t303 = A * t
           t305 = 0.2e1_dp * t303 * trhoa
           t306 = t302 + t305
           t310 = t101 * t111
           t311 = t118 ** 2
           t312 = 0.1e1_dp / t311
           t313 = t113 * t312
           t314 = A * t116
           t317 = t111 * t
           t318 = t115 * t317
           t321 = t302 + t305 + 0.2e1_dp * t314 * Arhoa + 0.4e1_dp * t318 * trhoa
           t324 = 0.2e1_dp * t297 * t299 + t101 * t111 * t306 * t119 - t310 *&
                t313 * t321
           t325 = 0.1e1_dp / t122
           t326 = t324 * t325
           epsilon_cGGArhoa = epsilon_c_unifrhoa + 0.3e1_dp * t293 * t294 + &
                t110 * t326
           t329 = t126 ** 2
           kf_arhoa = t124 / t329 * t90 / 0.3e1_dp
           ex_unif_arhoa = -0.3e1_dp / 0.4e1_dp * t7 * kf_arhoa
           t335 = kf_a ** 2
           t336 = 0.1e1_dp / t335
           t337 = my_norm_drhoa * t336
           t340 = my_rhoa ** 2
           t341 = 0.1e1_dp / t340
           s_arhoa = -t337 * t131 * kf_arhoa / 0.2e1_dp - t130 * t341 / 0.2e1_dp
           t344 = t137 ** 2
           t346 = 0.1e1_dp / t344 * mu
           Fx_arhoa = 0.2e1_dp * t346 * s_a * s_arhoa
           t350 = my_rhoa * ex_unif_arhoa

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ra(ii) = e_ra(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_a * Fx_a + 0.2e1_dp * &
                   t350 * Fx_a + 0.2e1_dp * t140 * Fx_arhoa) / 0.2e1_dp + scale_ec * (&
                   epsilon_cGGA + my_rho * epsilon_cGGArhoa)
           END IF

           chirhob = -t2 - t164
           rsrhob = rsrhoa
           t368 = t179 * rsrhob / 0.2e1_dp + beta_2_1 * rsrhob + 0.3e1_dp / &
                0.2e1_dp * t183 * rsrhob + t22 * t20 * rsrhob * t187
           e_c_u_0rhob = -0.2e1_dp * t171 * rsrhob * t28 + t177 * t368 * t191
           t382 = t201 * rsrhob / 0.2e1_dp + beta_2_2 * rsrhob + 0.3e1_dp / &
                0.2e1_dp * t205 * rsrhob + t40 * t38 * rsrhob * t187
           e_c_u_1rhob = -0.2e1_dp * t194 * rsrhob * t46 + t200 * t382 * t212
           t396 = t222 * rsrhob / 0.2e1_dp + beta_2_3 * rsrhob + 0.3e1_dp / &
                0.2e1_dp * t226 * rsrhob + t58 * t56 * rsrhob * t187
           alpha_crhob = 0.2e1_dp * t215 * rsrhob * t64 - t221 * t396 * t233
           frhob = (0.4e1_dp / 0.3e1_dp * t71 * chirhob - 0.4e1_dp / 0.3e1_dp&
                * t74 * chirhob) * t69
           t403 = alpha_crhob * f
           t405 = alpha_c * frhob
           t407 = t245 * chirhob
           t409 = 0.4e1_dp * t77 * t407
           t410 = e_c_u_1rhob - e_c_u_0rhob
           t411 = t410 * f
           t413 = t84 * frhob
           t415 = t244 * chirhob
           t417 = 0.4e1_dp * t85 * t415
           epsilon_c_unifrhob = e_c_u_0rhob + t403 * t82 + t405 * t82 - t409 &
                + t411 * t80 + t413 * t80 + t417
           phirhob = t257 * chirhob / 0.3e1_dp - t259 * chirhob / 0.3e1_dp
           k_frhob = k_frhoa
           k_srhob = t266 * k_frhob * t7
           trhob = -t269 * t98 * phirhob / 0.2e1_dp - t96 * t274 * k_srhob / &
                0.2e1_dp - t278 / 0.2e1_dp
           t427 = epsilon_c_unifrhob * t100
           t429 = t285 * phirhob
           t432 = -t427 * t105 + 0.3e1_dp * t102 * t429
           Arhob = -t101 * t281 * t432 * t107
           t436 = t123 * phirhob
           t439 = t298 * trhob
           t442 = Arhob * t111
           t444 = 0.2e1_dp * t303 * trhob
           t445 = t442 + t444
           t453 = t442 + t444 + 0.2e1_dp * t314 * Arhob + 0.4e1_dp * t318 * trhob
           t456 = 0.2e1_dp * t297 * t439 + t101 * t111 * t445 * t119 - t310 *&
                t313 * t453
           t457 = t456 * t325
           epsilon_cGGArhob = epsilon_c_unifrhob + 0.3e1_dp * t293 * t436 + &
                t110 * t457
           t460 = t143 ** 2
           kf_brhob = t124 / t460 * t90 / 0.3e1_dp
           ex_unif_brhob = -0.3e1_dp / 0.4e1_dp * t7 * kf_brhob
           t466 = kf_b ** 2
           t467 = 0.1e1_dp / t466
           t468 = my_norm_drhob * t467
           t471 = my_rhob ** 2
           t472 = 0.1e1_dp / t471
           s_brhob = -t468 * t148 * kf_brhob / 0.2e1_dp - t147 * t472 / 0.2e1_dp
           t475 = t153 ** 2
           t477 = 0.1e1_dp / t475 * mu
           Fx_brhob = 0.2e1_dp * t477 * s_b * s_brhob
           t481 = my_rhob * ex_unif_brhob

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_rb(ii) = e_rb(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_b * Fx_b + 0.2e1_dp * &
                   t481 * Fx_b + 0.2e1_dp * t156 * Fx_brhob) / 0.2e1_dp + scale_ec * (&
                   epsilon_cGGA + my_rho * epsilon_cGGArhob)
           END IF

           t488 = t95 * t97
           tnorm_drho = t488 * t2 / 0.2e1_dp
           t493 = t101 * t317
           t494 = A * tnorm_drho
           t502 = 0.2e1_dp * t303 * tnorm_drho + 0.4e1_dp * t318 * tnorm_drho
           t505 = 0.2e1_dp * t297 * t298 * tnorm_drho + 0.2e1_dp * t493 * &
                t494 * t119 - t310 * t313 * t502
           t506 = t505 * t325
           Hnorm_drho = t110 * t506

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndr(ii) = e_ndr(ii)+&
                   scale_ec * my_rho * Hnorm_drho
           END IF

           s_anorm_drhoa = t129 * t131 / 0.2e1_dp
           Fx_anorm_drhoa = 0.2e1_dp * t346 * s_a * s_anorm_drhoa

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndra(ii) = e_ndra(ii)+&
                   scale_ex * t140 * Fx_anorm_drhoa
           END IF

           s_bnorm_drhob = t146 * t148 / 0.2e1_dp
           Fx_bnorm_drhob = 0.2e1_dp * t477 * s_b * s_bnorm_drhob

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndrb(ii) = e_ndrb(ii)+&
                   scale_ex * t156 * Fx_bnorm_drhob
           END IF

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              t518 = 0.1e1_dp / t162 / my_rho
              t519 = t1 * t518
              chirhoarhoa = -0.2e1_dp * t163 + 0.2e1_dp * t519
              t523 = 0.1e1_dp / t90
              t525 = t162 ** 2
              rsrhoarhoa = -t6 / t165 / t8 * t523 / t525 / 0.18e2_dp + &
                   t6 * t167 * t518 / 0.6e1_dp
              t536 = alpha_1_1 * rsrhoa
              t538 = t176 * t190 * t191
              t543 = t12 / t175 / t23
              t544 = t190 ** 2
              t548 = 0.1e1_dp / t18
              t549 = beta_1_1 * t548
              t550 = rsrhoa ** 2
              t556 = beta_3_1 * t178
              t561 = t20 ** 2
              t563 = rs ** 2
              t564 = 0.1e1_dp / t563
              t576 = t175 ** 2
              t578 = t12 / t576
              t579 = t27 ** 2
              t580 = 0.1e1_dp / t579
              e_c_u_0rhoarhoa = -0.2e1_dp * t171 * rsrhoarhoa * t28 + 0.2e1_dp *&
                   t536 * t538 - 0.2e1_dp * t543 * t544 * t191 + t177 * (-t549 * t550 &
                   / 0.4e1_dp + t179 * rsrhoarhoa / 0.2e1_dp + beta_2_1 * rsrhoarhoa + &
                   0.3e1_dp / 0.4e1_dp * t556 * t550 + 0.3e1_dp / 0.2e1_dp * t183 * &
                   rsrhoarhoa + t22 * t561 * t550 * t564 + t22 * t20 * rsrhoarhoa * &
                   t187 - t22 * t20 * t550 * t564) * t191 + t578 * t544 * t580 * t14 / &
                   0.2e1_dp
              e_c_u_01rhoa = e_c_u_0rhoa
              t588 = alpha_1_2 * rsrhoa
              t590 = t199 * t211 * t212
              t595 = t32 / t198 / t41
              t596 = t211 ** 2
              t600 = beta_1_2 * t548
              t606 = beta_3_2 * t178
              t611 = t38 ** 2
              t624 = t198 ** 2
              t626 = t32 / t624
              t627 = t45 ** 2
              t628 = 0.1e1_dp / t627
              t636 = alpha_1_3 * rsrhoa
              t638 = t220 * t232 * t233
              t643 = t50 / t219 / t59
              t644 = t232 ** 2
              t648 = beta_1_3 * t548
              t654 = beta_3_3 * t178
              t659 = t56 ** 2
              t672 = t219 ** 2
              t674 = t50 / t672
              t675 = t63 ** 2
              t676 = 0.1e1_dp / t675
              alpha_c1rhoa = alpha_crhoa
              t681 = 0.1e1_dp / t87
              t682 = chirhoa ** 2
              t687 = 0.1e1_dp / t88
              frhoarhoa = (0.4e1_dp / 0.9e1_dp * t681 * t682 + 0.4e1_dp / &
                   0.3e1_dp * t71 * chirhoarhoa + 0.4e1_dp / 0.9e1_dp * t687 * t682 - &
                   0.4e1_dp / 0.3e1_dp * t74 * chirhoarhoa) * t69
              f1rhoa = frhoa
              t705 = alpha_c1rhoa * f
              t708 = alpha_c * f1rhoa
              t711 = t78 * t79
              t726 = e_c_u_1rhoa - e_c_u_01rhoa
              t733 = t726 * f
              t736 = t84 * f1rhoa
              t745 = -0.4e1_dp * t77 * t245 * chirhoarhoa + (-0.2e1_dp * t194 * &
                   rsrhoarhoa * t46 + 0.2e1_dp * t588 * t590 - 0.2e1_dp * t595 * t596 *&
                   t212 + t200 * (-t600 * t550 / 0.4e1_dp + t201 * rsrhoarhoa / &
                   0.2e1_dp + beta_2_2 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t606 * t550&
                   + 0.3e1_dp / 0.2e1_dp * t205 * rsrhoarhoa + t40 * t611 * t550 * &
                   t564 + t40 * t38 * rsrhoarhoa * t187 - t40 * t38 * t550 * t564) * &
                   t212 + t626 * t596 * t628 * t34 / 0.2e1_dp - e_c_u_0rhoarhoa) * f * &
                   t80 + t249 * f1rhoa * t80 + 0.4e1_dp * t250 * t254 + t726 * frhoa * &
                   t80 + t84 * frhoarhoa * t80 + 0.4e1_dp * t252 * t254 + 0.4e1_dp * &
                   t733 * t254 + 0.4e1_dp * t736 * t254 + 0.12e2_dp * t85 * t79 * t682 &
                   + 0.4e1_dp * t85 * t244 * chirhoarhoa
              epsilon_c_unifrhoarhoa = e_c_u_0rhoarhoa + (0.2e1_dp * t215 * &
                   rsrhoarhoa * t64 - 0.2e1_dp * t636 * t638 + 0.2e1_dp * t643 * t644 *&
                   t233 - t221 * (-t648 * t550 / 0.4e1_dp + t222 * rsrhoarhoa / &
                   0.2e1_dp + beta_2_3 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t654 * t550&
                   + 0.3e1_dp / 0.2e1_dp * t226 * rsrhoarhoa + t58 * t659 * t550 * &
                   t564 + t58 * t56 * rsrhoarhoa * t187 - t58 * t56 * t550 * t564) * &
                   t233 - t674 * t644 * t676 * t52 / 0.2e1_dp) * f * t82 + alpha_crhoa &
                   * f1rhoa * t82 - 0.4e1_dp * t240 * t246 + alpha_c1rhoa * frhoa * t82&
                   + alpha_c * frhoarhoa * t82 - 0.4e1_dp * t242 * t246 - 0.4e1_dp * &
                   t705 * t246 - 0.4e1_dp * t708 * t246 - 0.12e2_dp * t77 * t711 * t682&
                   + t745
              epsilon_c_unif1rhoa = e_c_u_01rhoa + t705 * t82 + t708 * t82 - &
                   t248 + t733 * t80 + t736 * t80 + t256
              t750 = 0.1e1_dp / t72
              t755 = 0.1e1_dp / t75
              phirhoarhoa = -t750 * t682 / 0.9e1_dp + t257 * chirhoarhoa / &
                   0.3e1_dp - t755 * t682 / 0.9e1_dp - t259 * chirhoarhoa / 0.3e1_dp
              phi1rhoa = phirhoa
              t763 = t90 ** 2
              k_frhoarhoa = -0.2e1_dp / 0.9e1_dp * t3 / t262 / t91 * t763
              t767 = 0.1e1_dp / t94 / t93
              t768 = k_frhoa ** 2
              k_s1rhoa = k_srhoa
              t775 = my_norm_drho * t105 * t97
              t776 = t2 * phirhoa
              t779 = t269 * t273
              t785 = t269 * t277 * phirhoa / 0.2e1_dp
              t789 = t2 * k_srhoa
              t795 = t96 / t272 / k_s
              t798 = t273 * t163
              t801 = t96 * t798 * k_srhoa / 0.2e1_dp
              t812 = t96 * t97 * t518
              trhoarhoa = t775 * t776 * phi1rhoa + t779 * t776 * k_s1rhoa / &
                   0.2e1_dp + t785 - t269 * t98 * phirhoarhoa / 0.2e1_dp + t779 * t789 &
                   * phi1rhoa / 0.2e1_dp + t795 * t789 * k_s1rhoa + t801 - t96 * t274 *&
                   (-t767 * t768 * t523 / 0.2e1_dp + t266 * k_frhoarhoa * t7) / &
                   0.2e1_dp + t269 * t277 * phi1rhoa / 0.2e1_dp + t96 * t798 * k_s1rhoa&
                   / 0.2e1_dp + t812
              t1rhoa = -t269 * t98 * phi1rhoa / 0.2e1_dp - t96 * t274 * k_s1rhoa&
                   / 0.2e1_dp - t278 / 0.2e1_dp
              t820 = t101 / t280 / t108
              t821 = t107 ** 2
              t822 = t289 * t821
              t823 = epsilon_c_unif1rhoa * t100
              t825 = t285 * phi1rhoa
              t828 = -t823 * t105 + 0.3e1_dp * t102 * t825
              t839 = 0.1e1_dp / t284 / phi
              t840 = t839 * phirhoa
              t851 = t101 * t281
              Arhoarhoa = 0.2e1_dp * t820 * t822 * t828 - t101 * t281 * (&
                   -epsilon_c_unifrhoarhoa * t100 * t105 + 0.3e1_dp * t282 * t825 + &
                   0.3e1_dp * t823 * t286 - 0.12e2_dp * t102 * t840 * phi1rhoa + &
                   0.3e1_dp * t102 * t285 * phirhoarhoa) * t107 - t851 * t289 * t828 * &
                   t107
              A1rhoa = -t101 * t281 * t828 * t107
              t858 = gamma_var * phi
              t865 = A1rhoa * t111
              t867 = 0.2e1_dp * t303 * t1rhoa
              t868 = t865 + t867
              t876 = t865 + t867 + 0.2e1_dp * t314 * A1rhoa + 0.4e1_dp * t318 * t1rhoa
              t879 = 0.2e1_dp * t297 * t298 * t1rhoa + t101 * t111 * t868 * t119&
                   - t310 * t313 * t876
              t880 = t879 * t325
              t904 = t306 * t119
              t908 = Arhoarhoa * t111
              t909 = Arhoa * t
              t911 = 0.2e1_dp * t909 * t1rhoa
              t914 = 0.2e1_dp * A1rhoa * t * trhoa
              t917 = 0.2e1_dp * A * t1rhoa * trhoa
              t919 = 0.2e1_dp * t303 * trhoarhoa
              t924 = t306 * t312
              t936 = t113 / t311 / t118
              t944 = A * t317
              t953 = t115 * t111
              t959 = t908 + t911 + t914 + t917 + t919 + 0.2e1_dp * A1rhoa * t116&
                   * Arhoa + 0.8e1_dp * t944 * Arhoa * t1rhoa + 0.2e1_dp * t314 * &
                   Arhoarhoa + 0.8e1_dp * t944 * trhoa * A1rhoa + 0.12e2_dp * t953 * &
                   trhoa * t1rhoa + 0.4e1_dp * t318 * trhoarhoa
              t962 = 0.2e1_dp * t101 * t1rhoa * t299 + 0.2e1_dp * t297 * t868 * &
                   t119 * trhoa - 0.2e1_dp * t297 * t313 * trhoa * t876 + 0.2e1_dp * &
                   t297 * t298 * trhoarhoa + 0.2e1_dp * t297 * t904 * t1rhoa + t101 * &
                   t111 * (t908 + t911 + t914 + t917 + t919) * t119 - t310 * t924 * &
                   t876 - 0.2e1_dp * t297 * t313 * t321 * t1rhoa - t310 * t868 * t312 *&
                   t321 + 0.2e1_dp * t310 * t936 * t321 * t876 - t310 * t313 * t959
              t965 = t122 ** 2
              t966 = 0.1e1_dp / t965
              t967 = t324 * t966
              kf_arhoarhoa = -0.2e1_dp / 0.9e1_dp * t124 / t329 / t125 * t763
              ex_unif_a1rhoa = ex_unif_arhoa
              t985 = kf_arhoa ** 2
              s_a1rhoa = s_arhoa
              t998 = mu ** 2
              t999 = 0.1e1_dp / t344 / t137 * t998
              t1000 = t999 * t133
              t1001 = s_arhoa * t135
              Fx_a1rhoa = 0.2e1_dp * t346 * s_a * s_a1rhoa

              e_ra_ra(ii) = e_ra_ra(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_a1rhoa * Fx_a + &
                   0.2e1_dp * ex_unif_a * Fx_a1rhoa + 0.2e1_dp * ex_unif_arhoa * Fx_a -&
                   0.3e1_dp / 0.2e1_dp * my_rhoa * t7 * kf_arhoarhoa * Fx_a + 0.2e1_dp * &
                   t350 * Fx_a1rhoa + 0.2e1_dp * ex_unif_a * Fx_arhoa + 0.2e1_dp * my_rhoa&
                   * ex_unif_a1rhoa * Fx_arhoa + 0.2e1_dp * t140 * (-0.8e1_dp * t1000 &
                   * t1001 * s_a1rhoa + 0.2e1_dp * t346 * s_a1rhoa * s_arhoa + 0.2e1_dp&
                   * t346 * s_a * (my_norm_drhoa / t335 / kf_a * t131 * t985 + t337 * &
                   t341 * kf_arhoa - t337 * t131 * kf_arhoarhoa / 0.2e1_dp + t130 / &
                   t340 / my_rhoa))) / 0.2e1_dp + scale_ec * (epsilon_c_unif1rhoa + &
                   0.3e1_dp * t293 * t123 * phi1rhoa + t110 * t880 + epsilon_cGGArhoa +&
                   my_rho * (epsilon_c_unifrhoarhoa + 0.6e1_dp * t858 * t294 * phi1rhoa +&
                   0.3e1_dp * t293 * t880 * phirhoa + 0.3e1_dp * t293 * t123 * &
                   phirhoarhoa + 0.3e1_dp * t293 * t326 * phi1rhoa + t110 * t962 * t325&
                   - t110 * t967 * t879))

              chirhoarhob = 0.2e1_dp * t519
              rsrhoarhob = rsrhoarhoa
              t1031 = t176 * t368 * t191
              t1033 = alpha_1_1 * rsrhob
              t1038 = rsrhoa * rsrhob
              t1050 = rsrhob * t564 * rsrhoa
              e_c_u_0rhoarhob = -0.2e1_dp * t171 * rsrhoarhob * t28 + t536 * &
                   t1031 + t1033 * t538 - 0.2e1_dp * t543 * t192 * t368 + t177 * (-t549&
                   * t1038 / 0.4e1_dp + t179 * rsrhoarhob / 0.2e1_dp + beta_2_1 * &
                   rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t556 * t1038 + 0.3e1_dp / &
                   0.2e1_dp * t183 * rsrhoarhob + t22 * t561 * t1050 + t22 * t20 * &
                   rsrhoarhob * t187 - t22 * t20 * t1050) * t191 + t578 * t190 * t580 *&
                   t14 * t368 / 0.2e1_dp
              t1069 = t199 * t382 * t212
              t1071 = alpha_1_2 * rsrhob
              t1104 = t220 * t396 * t233
              t1106 = alpha_1_3 * rsrhob
              frhoarhob = (0.4e1_dp / 0.9e1_dp * t681 * chirhoa * chirhob + &
                   0.4e1_dp / 0.3e1_dp * t71 * chirhoarhob + 0.4e1_dp / 0.9e1_dp * t687&
                   * chirhoa * chirhob - 0.4e1_dp / 0.3e1_dp * t74 * chirhoarhob) * &
                   t69
              t1164 = t79 * chirhoa * chirhob
              t1193 = -0.4e1_dp * t77 * t245 * chirhoarhob + (-0.2e1_dp * t194 *&
                   rsrhoarhob * t46 + t588 * t1069 + t1071 * t590 - 0.2e1_dp * t595 * &
                   t213 * t382 + t200 * (-t600 * t1038 / 0.4e1_dp + t201 * rsrhoarhob /&
                   0.2e1_dp + beta_2_2 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t606 * &
                   t1038 + 0.3e1_dp / 0.2e1_dp * t205 * rsrhoarhob + t40 * t611 * t1050&
                   + t40 * t38 * rsrhoarhob * t187 - t40 * t38 * t1050) * t212 + t626 &
                   * t211 * t628 * t34 * t382 / 0.2e1_dp - e_c_u_0rhoarhob) * f * t80 +&
                   t249 * frhob * t80 + 0.4e1_dp * t250 * t415 + t410 * frhoa * t80 + &
                   t84 * frhoarhob * t80 + 0.4e1_dp * t252 * t415 + 0.4e1_dp * t411 * &
                   t254 + 0.4e1_dp * t413 * t254 + 0.12e2_dp * t85 * t1164 + 0.4e1_dp *&
                   t85 * t244 * chirhoarhob
              epsilon_c_unifrhoarhob = e_c_u_0rhoarhob + (0.2e1_dp * t215 * &
                   rsrhoarhob * t64 - t636 * t1104 - t1106 * t638 + 0.2e1_dp * t643 * &
                   t234 * t396 - t221 * (-t648 * t1038 / 0.4e1_dp + t222 * rsrhoarhob /&
                   0.2e1_dp + beta_2_3 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t654 * &
                   t1038 + 0.3e1_dp / 0.2e1_dp * t226 * rsrhoarhob + t58 * t659 * t1050&
                   + t58 * t56 * rsrhoarhob * t187 - t58 * t56 * t1050) * t233 - t674 &
                   * t232 * t676 * t52 * t396 / 0.2e1_dp) * f * t82 + alpha_crhoa * &
                   frhob * t82 - 0.4e1_dp * t240 * t407 + alpha_crhob * frhoa * t82 + &
                   alpha_c * frhoarhob * t82 - 0.4e1_dp * t242 * t407 - 0.4e1_dp * t403&
                   * t246 - 0.4e1_dp * t405 * t246 - 0.12e2_dp * t77 * t78 * t1164 + &
                   t1193
              phirhoarhob = -t750 * chirhoa * chirhob / 0.9e1_dp + t257 * &
                   chirhoarhob / 0.3e1_dp - t755 * chirhoa * chirhob / 0.9e1_dp - t259 &
                   * chirhoarhob / 0.3e1_dp
              k_frhoarhob = k_frhoarhoa
              t1228 = t269 * t277 * phirhob / 0.2e1_dp
              t1231 = t96 * t798 * k_srhob / 0.2e1_dp
              trhoarhob = t775 * t776 * phirhob + t779 * t776 * k_srhob / &
                   0.2e1_dp + t785 - t269 * t98 * phirhoarhob / 0.2e1_dp + t779 * t789 &
                   * phirhob / 0.2e1_dp + t795 * t789 * k_srhob + t801 - t96 * t274 * (&
                   -t767 * k_frhoa * t523 * k_frhob / 0.2e1_dp + t266 * k_frhoarhob * &
                   t7) / 0.2e1_dp + t1228 + t1231 + t812
              Arhoarhob = 0.2e1_dp * t820 * t822 * t432 - t101 * t281 * (&
                   -epsilon_c_unifrhoarhob * t100 * t105 + 0.3e1_dp * t282 * t429 + &
                   0.3e1_dp * t427 * t286 - 0.12e2_dp * t102 * t840 * phirhob + &
                   0.3e1_dp * t102 * t285 * phirhoarhob) * t107 - t851 * t289 * t432 * &
                   t107
              t1269 = t445 * t119
              t1283 = Arhoarhob * t111
              t1285 = 0.2e1_dp * t909 * trhob
              t1286 = Arhob * t
              t1288 = 0.2e1_dp * t1286 * trhoa
              t1291 = 0.2e1_dp * A * trhob * trhoa
              t1293 = 0.2e1_dp * t303 * trhoarhob
              t1304 = t445 * t312
              t1327 = t1283 + t1285 + t1288 + t1291 + t1293 + 0.2e1_dp * Arhob *&
                   t116 * Arhoa + 0.8e1_dp * t944 * Arhoa * trhob + 0.2e1_dp * t314 * &
                   Arhoarhob + 0.8e1_dp * t944 * trhoa * Arhob + 0.12e2_dp * t953 * &
                   trhoa * trhob + 0.4e1_dp * t318 * trhoarhob
              t1330 = 0.2e1_dp * t101 * trhob * t299 + 0.2e1_dp * t297 * t1269 *&
                   trhoa - 0.2e1_dp * t297 * t313 * trhoa * t453 + 0.2e1_dp * t297 * &
                   t298 * trhoarhob + 0.2e1_dp * t297 * t904 * trhob + t101 * t111 * (&
                   t1283 + t1285 + t1288 + t1291 + t1293) * t119 - t310 * t924 * t453 -&
                   0.2e1_dp * t297 * t313 * t321 * trhob - t310 * t1304 * t321 + &
                   0.2e1_dp * t310 * t936 * t321 * t453 - t310 * t313 * t1327

              e_ra_rb(ii) = e_ra_rb(ii)+&
                   scale_ec * (epsilon_cGGArhob + epsilon_cGGArhoa + &
                   my_rho * (epsilon_c_unifrhoarhob + 0.6e1_dp * t858 * t294 * phirhob + &
                   0.3e1_dp * t293 * t457 * phirhoa + 0.3e1_dp * t293 * t123 * &
                   phirhoarhob + 0.3e1_dp * t293 * t326 * phirhob + t110 * t1330 * t325&
                   - t110 * t967 * t456))

              chirhobrhob = 0.2e1_dp * t163 + 0.2e1_dp * t519
              rsrhobrhob = rsrhoarhob
              t1342 = t368 ** 2
              t1346 = rsrhob ** 2
              e_c_u_0rhobrhob = -0.2e1_dp * t171 * rsrhobrhob * t28 + 0.2e1_dp *&
                   t1033 * t1031 - 0.2e1_dp * t543 * t1342 * t191 + t177 * (-t549 * &
                   t1346 / 0.4e1_dp + t179 * rsrhobrhob / 0.2e1_dp + beta_2_1 * &
                   rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t556 * t1346 + 0.3e1_dp / &
                   0.2e1_dp * t183 * rsrhobrhob + t22 * t561 * t1346 * t564 + t22 * t20&
                   * rsrhobrhob * t187 - t22 * t20 * t1346 * t564) * t191 + t578 * &
                   t1342 * t580 * t14 / 0.2e1_dp
              e_c_u_01rhob = e_c_u_0rhob
              t1377 = t382 ** 2
              t1411 = t396 ** 2
              alpha_c1rhob = alpha_crhob
              t1440 = chirhob ** 2
              frhobrhob = (0.4e1_dp / 0.9e1_dp * t681 * t1440 + 0.4e1_dp / &
                   0.3e1_dp * t71 * chirhobrhob + 0.4e1_dp / 0.9e1_dp * t687 * t1440 - &
                   0.4e1_dp / 0.3e1_dp * t74 * chirhobrhob) * t69
              f1rhob = frhob
              t1462 = alpha_c1rhob * f
              t1465 = alpha_c * f1rhob
              t1482 = e_c_u_1rhob - e_c_u_01rhob
              t1489 = t1482 * f
              t1492 = t84 * f1rhob
              t1501 = -0.4e1_dp * t77 * t245 * chirhobrhob + (-0.2e1_dp * t194 *&
                   rsrhobrhob * t46 + 0.2e1_dp * t1071 * t1069 - 0.2e1_dp * t595 * &
                   t1377 * t212 + t200 * (-t600 * t1346 / 0.4e1_dp + t201 * rsrhobrhob &
                   / 0.2e1_dp + beta_2_2 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t606 * &
                   t1346 + 0.3e1_dp / 0.2e1_dp * t205 * rsrhobrhob + t40 * t611 * t1346&
                   * t564 + t40 * t38 * rsrhobrhob * t187 - t40 * t38 * t1346 * t564) &
                   * t212 + t626 * t1377 * t628 * t34 / 0.2e1_dp - e_c_u_0rhobrhob) * f&
                   * t80 + t410 * f1rhob * t80 + 0.4e1_dp * t411 * t415 + t1482 * &
                   frhob * t80 + t84 * frhobrhob * t80 + 0.4e1_dp * t413 * t415 + &
                   0.4e1_dp * t1489 * t415 + 0.4e1_dp * t1492 * t415 + 0.12e2_dp * t85 &
                   * t79 * t1440 + 0.4e1_dp * t85 * t244 * chirhobrhob
              epsilon_c_unifrhobrhob = e_c_u_0rhobrhob + (0.2e1_dp * t215 * &
                   rsrhobrhob * t64 - 0.2e1_dp * t1106 * t1104 + 0.2e1_dp * t643 * &
                   t1411 * t233 - t221 * (-t648 * t1346 / 0.4e1_dp + t222 * rsrhobrhob &
                   / 0.2e1_dp + beta_2_3 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t654 * &
                   t1346 + 0.3e1_dp / 0.2e1_dp * t226 * rsrhobrhob + t58 * t659 * t1346&
                   * t564 + t58 * t56 * rsrhobrhob * t187 - t58 * t56 * t1346 * t564) &
                   * t233 - t674 * t1411 * t676 * t52 / 0.2e1_dp) * f * t82 + &
                   alpha_crhob * f1rhob * t82 - 0.4e1_dp * t403 * t407 + alpha_c1rhob *&
                   frhob * t82 + alpha_c * frhobrhob * t82 - 0.4e1_dp * t405 * t407 - &
                   0.4e1_dp * t1462 * t407 - 0.4e1_dp * t1465 * t407 - 0.12e2_dp * t77 &
                   * t711 * t1440 + t1501
              epsilon_c_unif1rhob = e_c_u_01rhob + t1462 * t82 + t1465 * t82 - &
                   t409 + t1489 * t80 + t1492 * t80 + t417
              phirhobrhob = -t750 * t1440 / 0.9e1_dp + t257 * chirhobrhob / &
                   0.3e1_dp - t755 * t1440 / 0.9e1_dp - t259 * chirhobrhob / 0.3e1_dp
              phi1rhob = phirhob
              t1514 = k_frhob ** 2
              k_s1rhob = k_srhob
              t1520 = t2 * phirhob
              t1529 = t2 * k_srhob
              trhobrhob = t775 * t1520 * phi1rhob + t779 * t1520 * k_s1rhob / &
                   0.2e1_dp + t1228 - t269 * t98 * phirhobrhob / 0.2e1_dp + t779 * &
                   t1529 * phi1rhob / 0.2e1_dp + t795 * t1529 * k_s1rhob + t1231 - t96 &
                   * t274 * (-t767 * t1514 * t523 / 0.2e1_dp + t266 * k_frhoarhob * t7)&
                   / 0.2e1_dp + t269 * t277 * phi1rhob / 0.2e1_dp + t96 * t798 * &
                   k_s1rhob / 0.2e1_dp + t812
              t1rhob = -t269 * t98 * phi1rhob / 0.2e1_dp - t96 * t274 * k_s1rhob&
                   / 0.2e1_dp - t278 / 0.2e1_dp
              t1550 = epsilon_c_unif1rhob * t100
              t1552 = t285 * phi1rhob
              t1555 = -t1550 * t105 + 0.3e1_dp * t102 * t1552
              Arhobrhob = 0.2e1_dp * t820 * t432 * t821 * t1555 - t101 * t281 * &
                   (-epsilon_c_unifrhobrhob * t100 * t105 + 0.3e1_dp * t427 * t1552 + &
                   0.3e1_dp * t1550 * t429 - 0.12e2_dp * t102 * t839 * phirhob * &
                   phi1rhob + 0.3e1_dp * t102 * t285 * phirhobrhob) * t107 - t851 * &
                   t432 * t1555 * t107
              A1rhob = -t101 * t281 * t1555 * t107
              t1588 = A1rhob * t111
              t1590 = 0.2e1_dp * t303 * t1rhob
              t1591 = t1588 + t1590
              t1599 = t1588 + t1590 + 0.2e1_dp * t314 * A1rhob + 0.4e1_dp * t318&
                   * t1rhob
              t1602 = 0.2e1_dp * t297 * t298 * t1rhob + t101 * t111 * t1591 * &
                   t119 - t310 * t313 * t1599
              t1603 = t1602 * t325
              t1630 = Arhobrhob * t111
              t1632 = 0.2e1_dp * t1286 * t1rhob
              t1635 = 0.2e1_dp * A1rhob * t * trhob
              t1638 = 0.2e1_dp * A * t1rhob * trhob
              t1640 = 0.2e1_dp * t303 * trhobrhob
              t1674 = t1630 + t1632 + t1635 + t1638 + t1640 + 0.2e1_dp * A1rhob &
                   * t116 * Arhob + 0.8e1_dp * t944 * Arhob * t1rhob + 0.2e1_dp * t314 &
                   * Arhobrhob + 0.8e1_dp * t944 * trhob * A1rhob + 0.12e2_dp * t953 * &
                   trhob * t1rhob + 0.4e1_dp * t318 * trhobrhob
              t1677 = 0.2e1_dp * t101 * t1rhob * t439 + 0.2e1_dp * t297 * t1591 &
                   * t119 * trhob - 0.2e1_dp * t297 * t313 * trhob * t1599 + 0.2e1_dp *&
                   t297 * t298 * trhobrhob + 0.2e1_dp * t297 * t1269 * t1rhob + t101 *&
                   t111 * (t1630 + t1632 + t1635 + t1638 + t1640) * t119 - t310 * &
                   t1304 * t1599 - 0.2e1_dp * t297 * t313 * t453 * t1rhob - t310 * &
                   t1591 * t312 * t453 + 0.2e1_dp * t310 * t936 * t453 * t1599 - t310 *&
                   t313 * t1674
              t1680 = t456 * t966
              kf_brhobrhob = -0.2e1_dp / 0.9e1_dp * t124 / t460 / t142 * t763
              ex_unif_b1rhob = ex_unif_brhob
              t1698 = kf_brhob ** 2
              s_b1rhob = s_brhob
              t1711 = 0.1e1_dp / t475 / t153 * t998
              t1712 = t1711 * t150
              t1713 = s_brhob * t135
              Fx_b1rhob = 0.2e1_dp * t477 * s_b * s_b1rhob

              e_rb_rb(ii) = e_rb_rb(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_b1rhob * Fx_b + &
                   0.2e1_dp * ex_unif_b * Fx_b1rhob + 0.2e1_dp * ex_unif_brhob * Fx_b -&
                   0.3e1_dp / 0.2e1_dp * my_rhob * t7 * kf_brhobrhob * Fx_b + 0.2e1_dp * &
                   t481 * Fx_b1rhob + 0.2e1_dp * ex_unif_b * Fx_brhob + 0.2e1_dp * my_rhob&
                   * ex_unif_b1rhob * Fx_brhob + 0.2e1_dp * t156 * (-0.8e1_dp * t1712 &
                   * t1713 * s_b1rhob + 0.2e1_dp * t477 * s_b1rhob * s_brhob + 0.2e1_dp&
                   * t477 * s_b * (my_norm_drhob / t466 / kf_b * t148 * t1698 + t468 * &
                   t472 * kf_brhob - t468 * t148 * kf_brhobrhob / 0.2e1_dp + t147 / &
                   t471 / my_rhob))) / 0.2e1_dp + scale_ec * (epsilon_c_unif1rhob + &
                   0.3e1_dp * t293 * t123 * phi1rhob + t110 * t1603 + epsilon_cGGArhob &
                   + my_rho * (epsilon_c_unifrhobrhob + 0.6e1_dp * t858 * t436 * phi1rhob &
                   + 0.3e1_dp * t293 * t1603 * phirhob + 0.3e1_dp * t293 * t123 * &
                   phirhobrhob + 0.3e1_dp * t293 * t457 * phi1rhob + t110 * t1677 * &
                   t325 - t110 * t1680 * t1602))
              t1739 = t268 * t97
              t1741 = t95 * t273
              t1743 = t488 * t163
              trhoanorm_drho = -t1739 * t776 / 0.2e1_dp - t1741 * t789 / &
                   0.2e1_dp - t1743 / 0.2e1_dp
              t1748 = t101 * tnorm_drho
              t1765 = t909 * tnorm_drho
              t1766 = t494 * trhoa
              t1767 = t303 * trhoanorm_drho
              t1801 = 0.2e1_dp * t1748 * t299 + 0.4e1_dp * t310 * t494 * t119 * &
                   trhoa - 0.2e1_dp * t297 * t313 * trhoa * t502 + 0.2e1_dp * t297 * &
                   t298 * trhoanorm_drho + 0.2e1_dp * t297 * t904 * tnorm_drho + t101 *&
                   t111 * (0.2e1_dp * t1765 + 0.2e1_dp * t1766 + 0.2e1_dp * t1767) * &
                   t119 - t310 * t924 * t502 - 0.2e1_dp * t297 * t313 * t321 * &
                   tnorm_drho - 0.2e1_dp * t493 * t494 * t312 * t321 + 0.2e1_dp * t310 &
                   * t936 * t321 * t502 - t310 * t313 * (0.2e1_dp * t1765 + 0.2e1_dp * &
                   t1766 + 0.2e1_dp * t1767 + 0.8e1_dp * t944 * Arhoa * tnorm_drho + &
                   0.12e2_dp * t953 * trhoa * tnorm_drho + 0.4e1_dp * t318 * &
                   trhoanorm_drho)

              e_ra_ndr(ii) = e_ra_ndr(ii)+&
                   scale_ec * (Hnorm_drho + my_rho * (0.3e1_dp * &
                   t293 * t506 * phirhoa + t110 * t1801 * t325 - t110 * t967 * t505))

              trhobnorm_drho = -t1739 * t1520 / 0.2e1_dp - t1741 * t1529 / &
                   0.2e1_dp - t1743 / 0.2e1_dp
              t1829 = t1286 * tnorm_drho
              t1830 = t494 * trhob
              t1831 = t303 * trhobnorm_drho
              t1865 = 0.2e1_dp * t1748 * t439 + 0.4e1_dp * t310 * t494 * t119 * &
                   trhob - 0.2e1_dp * t297 * t313 * trhob * t502 + 0.2e1_dp * t297 * &
                   t298 * trhobnorm_drho + 0.2e1_dp * t297 * t1269 * tnorm_drho + t101 &
                   * t111 * (0.2e1_dp * t1829 + 0.2e1_dp * t1830 + 0.2e1_dp * t1831) * &
                   t119 - t310 * t1304 * t502 - 0.2e1_dp * t297 * t313 * t453 * &
                   tnorm_drho - 0.2e1_dp * t493 * t494 * t312 * t453 + 0.2e1_dp * t310 &
                   * t936 * t453 * t502 - t310 * t313 * (0.2e1_dp * t1829 + 0.2e1_dp * &
                   t1830 + 0.2e1_dp * t1831 + 0.8e1_dp * t944 * Arhob * tnorm_drho + &
                   0.12e2_dp * t953 * trhob * tnorm_drho + 0.4e1_dp * t318 * &
                   trhobnorm_drho)

              e_rb_ndr(ii) = e_rb_ndr(ii)+&
                   scale_ec * (Hnorm_drho + my_rho * (0.3e1_dp * &
                   t293 * t506 * phirhob + t110 * t1865 * t325 - t110 * t1680 * t505))

              t1871 = tnorm_drho ** 2
              t1876 = A * t1871
              t1888 = t502 ** 2
              t1901 = t505 ** 2

              e_ndr_ndr(ii) = e_ndr_ndr(ii)+&
                   scale_ec * my_rho * (t110 * (0.2e1_dp * &
                   t101 * t1871 * t113 * t119 + 0.10e2_dp * t310 * t1876 * t119 - &
                   0.4e1_dp * t297 * t313 * tnorm_drho * t502 - 0.4e1_dp * t493 * t494 &
                   * t312 * t502 + 0.2e1_dp * t310 * t936 * t1888 - t310 * t313 * (&
                   0.2e1_dp * t1876 + 0.12e2_dp * t953 * t1871)) * t325 - t110 * t1901 &
                   * t966)

              e_ra_ndra(ii) = e_ra_ndra(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_a * &
                   Fx_anorm_drhoa + 0.2e1_dp * t350 * Fx_anorm_drhoa + 0.2e1_dp * t140 &
                   * (-0.8e1_dp * t1000 * t1001 * s_anorm_drhoa + 0.2e1_dp * t346 * &
                   s_anorm_drhoa * s_arhoa + 0.2e1_dp * t346 * s_a * (-t336 * t131 * &
                   kf_arhoa / 0.2e1_dp - t129 * t341 / 0.2e1_dp))) / 0.2e1_dp

              t1922 = s_anorm_drhoa ** 2

              e_ndra_ndra(ii) = e_ndra_ndra(ii)+&
                   scale_ex * t140 * (-0.8e1_dp * t999 * &
                   t133 * t1922 * t135 + 0.2e1_dp * t346 * t1922)

              e_rb_ndrb(ii) = e_rb_ndrb(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_b * &
                   Fx_bnorm_drhob + 0.2e1_dp * t481 * Fx_bnorm_drhob + 0.2e1_dp * t156 &
                   * (-0.8e1_dp * t1712 * t1713 * s_bnorm_drhob + 0.2e1_dp * t477 * &
                   s_bnorm_drhob * s_brhob + 0.2e1_dp * t477 * s_b * (-t467 * t148 * &
                   kf_brhob / 0.2e1_dp - t146 * t472 / 0.2e1_dp))) / 0.2e1_dp


              t1949 = s_bnorm_drhob ** 2
              e_ndrb_ndrb(ii) = e_ndrb_ndrb(ii)+&
                   scale_ex * t156 * (-0.8e1_dp * t1711 *&
                   t150 * t1949 * t135 + 0.2e1_dp * t477 * t1949)
           END IF
        END IF
     END DO
  END SELECT
END SUBROUTINE pbe_lsd_calc


!!****f* xc_pbe/pbe_lsd_calc [1.0] *
!!
!!   NAME
!!     pbe_lsd_calc
!!
!!   FUNCTION
!!     low level calculation of the pbe exchange-correlation functional for lsd
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rhoa,rhob: alpha or beta spin density
!!     - norm_drhoa,norm_drhob,norm_drho: || grad rhoa |||,| grad rhoa ||,
!!       || grad (rhoa+rhob) ||
!!     - e_0: adds to it the local value of the functional
!!     - e_*: derivative of the functional wrt. to the variables
!!       named where the * is.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pbe_lsd_calc_old(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
     e_0, e_ra, e_rb, e_ra_ndra, e_rb_ndrb, e_ndr_ndr,&
     e_ndra_ndra, e_ndrb_ndrb, e_ndr,&
     e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr,&
     grad_deriv,npoints,epsilon_rho,epsilon_drho,error)
  REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drho, &
       norm_drhoa, norm_drhob
  REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_ra, e_rb, e_ra_ndra, &
       e_rb_ndrb, e_ndr_ndr, e_ndra_ndra, e_ndrb_ndrb,&
       e_ndr, e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr
  INTEGER, INTENT(in)                      :: grad_deriv, npoints
  REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_calc_old', &
       routineP = moduleN//':'//routineN
  REAL(kind=dp), PARAMETER :: small = 1.0e-20, t_1_9 = 1.0_dp/9.0_dp, &
       t_2_3 = 2.0_dp/3.0_dp, t_4_27 = 4.0_dp/27.0_dp, t_4_3 = 4.0_dp/3.0_dp, &
       t_4_9 = 4.0_dp/9.0_dp, t_8_27 = 8.0_dp/27.0_dp

  INTEGER                                  :: ii
  LOGICAL                                  :: failure

  real(kind=dp) :: my_rhoa, my_rhob, my_norm_drho, my_norm_drhoa, &
       my_norm_drhob, my_rho, t1, t2, chi, t3, t4, t5, t6, t7, t8, t9, rs, &
       t12, t13, t16, t18, t20, t23, t24, e_c_u_0, t28, t33, t36, t37, t41,&
       t46, t49, t50, alpha_c, t52, t55, t56, t57, t58, t59, t60, t61, f, &
       t63, t64, t65, t66, t68, t70, t71, epsilon_c_unif, t73, t75, t76, &
       gamma_var, t77, t78, phi, t80, t81, t82, t83, k_s, t84, t85, t86, &
       t87, t, t89, t90, t91, t92, t93, t95, t96, A, t99, t100, t101, t102,&
       t103, t104, t105, t107, t108, t109, t112, t113, epsilon_cGGA, mu, &
       t114, t115, t116, kf_a, ex_unif_a, t119, t120, t121, s_a, t123, &
       t126, Fx_a, t129, t131, t132, kf_b, ex_unif_b, t135, t136, t137, &
       s_b, t139, t142, Fx_b, t145, t149, t150, t151, chirhoa, t152, &
       t154, rsrhoa, t160, t161, t162, t163, t164, t167, t169, t170, t172, &
       t173, t174, e_c_u_0rhoa, t179, t180, t181, t186, t187, t188, &
       e_c_u_1rhoa, t193, t194, t195, t200, t201, t202, alpha_crhoa, frhoa,&
       t209, t211, t213, t214, t215, t217, t218, t219, t221, t223, t225
  real(kind=dp) :: &
       epsilon_c_unifrhoa, t226, t228, phirhoa, t231, k_frhoa, t235, &
       k_srhoa, t237, t238, t241, t242, t243, t246, t247, trhoa, t249, &
       t251, t252, t254, t255, t256, t259, Arhoa, t263, t264, t267, t268, &
       t271, t272, t274, t275, t276, t279, t280, t281, t282, t285, t286, &
       t289, t293, t294, t295, epsilon_cGGArhoa, t298, kf_arhoa, &
       ex_unif_arhoa, t304, t305, t306, t309, t310, s_arhoa, t313, t315, &
       Fx_arhoa, t319, chirhob, rsrhob, t326, t329, t331, t333, &
       t334, e_c_u_0rhob, t343, t344, e_c_u_1rhob, t353, t354, alpha_crhob,&
       frhob, t361, t363, t365, t367, t368, t369, t371, t373, t375, &
       epsilon_c_unifrhob, phirhob, k_frhob, k_srhob, trhob, t385, t387, &
       t390, Arhob, t394, t397, t400, t402, t403, t404, t411, t415, t416, &
       epsilon_cGGArhob, t419, kf_brhob, ex_unif_brhob, t425, t426, t427
  real(kind=dp) :: &
       t430, t431, s_brhob, t434, t436, Fx_brhob, t440, t445, &
       tnorm_drho, t450, t451, t459, t463, t464, Hnorm_drho,&
       s_anorm_drhoa, Fx_anorm_drhoa, s_bnorm_drhob, &
       Fx_bnorm_drhob, t474, t475, chirhoarhoa, t480, &
       rsrhoarhoa, t490, t495, t496, t500, t501, t502, t504, t507, t509, &
       t512, t518, t520, t521, t522, e_c_u_0rhoarhoa, e_c_u_01rhoa, t528, &
       t533, t534, t549, t551, t552, t553, t559, t564, t565, t580, t582, &
       t583, t584, alpha_c1rhoa, t588, t589, t594, frhoarhoa, f1rhoa, t612
  real(kind=dp) :: &
       t615, t618, t633, t640, t643, t652, epsilon_c_unifrhoarhoa, &
       epsilon_c_unif1rhoa, t657, t662, phirhoarhoa, phi1rhoa, t670, &
       k_frhoarhoa, t674, t675, k_s1rhoa, t682, t683, t686, t692, t696, &
       t702, t705, t708, t719, trhoarhoa, t1rhoa, t727, t728, t729, t730, &
       t732, t735, t746, t747, Arhoarhoa, A1rhoa, t765, t772, t774, t775, &
       t776, t783, t787, t788, t804, t805, t815, t816, t818, t821, t824, &
       t826, t831, t835, t843, t845, t846, t853, t862, t868, t872, t875, &
       t876, t877, kf_arhoarhoa, ex_unif_a1rhoa, t895, s_a1rhoa, t908, &
       t909, t910, Fx_a1rhoa, chirhoarhob, rsrhoarhob, t938,&
       t945, t947, t950, t952, t954, t956, e_c_u_0rhoarhob, t970, t995, &
       frhoarhob, t1045, t1074, epsilon_c_unifrhoarhob, phirhoarhob, &
       k_frhoarhob, t1109, t1112, trhoarhob, Arhoarhob, t1163, t1165
  real(kind=dp) :: &
       t1166, t1168, t1171, t1173, t1184, t1207, t1211, &
       chirhobrhob, rsrhobrhob, t1222, t1226, t1227, t1229, t1232, t1234, &
       t1237, e_c_u_0rhobrhob, e_c_u_01rhob, t1250, t1272, alpha_c1rhob, &
       t1290, frhobrhob, f1rhob, t1312, t1315, t1332, t1339, t1342, t1351, &
       epsilon_c_unifrhobrhob, epsilon_c_unif1rhob, phirhobrhob, phi1rhob, &
       t1364, k_s1rhob, t1370, t1379, trhobrhob, t1rhob, t1400, t1402, &
       t1405, Arhobrhob, A1rhob, t1439, t1441, t1442, t1443, t1450, t1454, &
       t1455, t1471, t1481, t1483, t1486, t1489, t1491, t1499, t1507, &
       t1527, t1531, t1534, kf_brhobrhob, ex_unif_b1rhob, t1552, s_b1rhob, &
       t1565, t1566, Fx_b1rhob, t1590, t1592, t1594, &
       trhoanorm_drho, t1599, t1602, t1603, t1616, t1617, t1618, t1627, &
       t1630, t1651, trhobnorm_drho, t1677, t1678, &
       t1679, t1688, t1711, t1717, t1721, t1726, t1732, t1746, t1766, &
       t1790

  failure=.FALSE.

  SELECT CASE(grad_deriv)
  CASE default
     DO ii=1,npoints
        my_rhoa=MAX(rhoa(ii),0.0_dp)
        my_rhob=MAX(rhob(ii),0.0_dp)
        my_rho=my_rhoa+my_rhob
        IF (my_rho>epsilon_rho) THEN
           my_rhoa=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhoa)
           my_rhob=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhob)
           my_rho=my_rhoa+my_rhob
           my_norm_drho=norm_drho(ii)
           my_norm_drhoa=norm_drhoa(ii)
           my_norm_drhob=norm_drhob(ii)

           t1 = my_rhoa - my_rhob
           t2 = 0.1e1_dp / my_rho
           chi = t1 * t2
           t3 = 3 ** (0.1e1_dp / 0.3e1_dp)
           t4 = 4 ** (0.1e1_dp / 0.3e1_dp)
           t5 = t4 ** 2
           t6 = t3 * t5
           t7 = 0.1e1_dp / 0.3141592654e1_dp
           t8 = t7 * t2
           t9 = t8 ** (0.1e1_dp / 0.3e1_dp)
           rs = t6 * t9 / 0.4e1_dp
           t12 = 0.1e1_dp + 0.21370e0_dp * rs
           t13 = sqrt(rs)
           t16 = t13 * rs
           t18 = rs ** 0.20e1_dp
           t20 = 0.75957e1_dp * t13 + 0.35876e1_dp * rs + 0.16382e1_dp * t16 &
                + 0.49294e0_dp * t18
           t23 = 0.1e1_dp + 0.1608182432e2_dp / t20
           t24 = log(t23)
           e_c_u_0 = -0.62182e-1_dp * t12 * t24
           t28 = 0.1e1_dp + 0.20548e0_dp * rs
           t33 = 0.141189e2_dp * t13 + 0.61977e1_dp * rs + 0.33662e1_dp * t16&
                + 0.62517e0_dp * t18
           t36 = 0.1e1_dp + 0.3216468318e2_dp / t33
           t37 = log(t36)
           t41 = 0.1e1_dp + 0.11125e0_dp * rs
           t46 = 0.10357e2_dp * t13 + 0.36231e1_dp * rs + 0.88026e0_dp * t16 &
                + 0.49671e0_dp * t18
           t49 = 0.1e1_dp + 0.2960857464e1_dp / t46
           t50 = log(t49)
           alpha_c = 0.33774e0_dp * t41 * t50
           t52 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t55 = 1 / (2 * t52 - 2)
           t56 = 0.1e1_dp + chi
           t57 = t56 ** (0.1e1_dp / 0.3e1_dp)
           t58 = t57 * t56
           t59 = 0.1e1_dp - chi
           t60 = t59 ** (0.1e1_dp / 0.3e1_dp)
           t61 = t60 * t59
           f = (t58 + t61 - 0.2e1_dp) * t55
           t63 = alpha_c * f
           t64 = 0.9e1_dp / 0.8e1_dp / t55
           t65 = chi ** 2
           t66 = t65 ** 2
           t68 = t64 * (0.1e1_dp - t66)
           t70 = -0.31090e-1_dp * t28 * t37 - e_c_u_0
           t71 = t70 * f
           epsilon_c_unif = e_c_u_0 + t63 * t68 + t71 * t66
           t73 = log(0.2e1_dp)
           t75 = 0.3141592654e1_dp ** 2
           t76 = 0.1e1_dp / t75
           gamma_var = (0.1e1_dp - t73) * t76
           t77 = t57 ** 2
           t78 = t60 ** 2
           phi = t77 / 0.2e1_dp + t78 / 0.2e1_dp
           t80 = t75 * my_rho
           t81 = t80 ** (0.1e1_dp / 0.3e1_dp)
           t82 = t3 * t81 * t7
           t83 = sqrt(t82)
           k_s = 0.2e1_dp * t83
           t84 = 0.1e1_dp / phi
           t85 = my_norm_drho * t84
           t86 = 0.1e1_dp / k_s
           t87 = t86 * t2
           t = t85 * t87 / 0.2e1_dp
           t89 = 0.1e1_dp / gamma_var
           t90 = epsilon_c_unif * t89
           t91 = phi ** 2
           t92 = t91 * phi
           t93 = 0.1e1_dp / t92
           t95 = exp(-t90 * t93)
           t96 = t95 - 0.1e1_dp
           A = 0.66725e-1_dp * t89 / t96
           t99 = gamma_var * t92
           t100 = t ** 2
           t101 = t89 * t100
           t102 = A * t100
           t103 = 0.1e1_dp + t102
           t104 = A ** 2
           t105 = t100 ** 2
           t107 = 0.1e1_dp + t102 + t104 * t105
           t108 = 0.1e1_dp / t107
           t109 = t103 * t108
           t112 = 0.1e1_dp + 0.66725e-1_dp * t101 * t109
           t113 = log(t112)
           epsilon_cGGA = epsilon_c_unif + t99 * t113
           mu = 0.2224166667e-1_dp * t75
           t114 = t3 * t52
           t115 = t75 * my_rhoa
           t116 = t115 ** (0.1e1_dp / 0.3e1_dp)
           kf_a = t114 * t116
           ex_unif_a = -0.3e1_dp / 0.4e1_dp * t7 * kf_a
           t119 = 0.1e1_dp / kf_a
           t120 = my_norm_drhoa * t119
           t121 = 0.1e1_dp / my_rhoa
           s_a = t120 * t121 / 0.2e1_dp
           t123 = s_a ** 2
           t126 = 0.1e1_dp + 0.1243781095e1_dp * mu * t123
           Fx_a = 0.1804e1_dp - 0.804e0_dp / t126
           t129 = my_rhoa * ex_unif_a
           t131 = t75 * my_rhob
           t132 = t131 ** (0.1e1_dp / 0.3e1_dp)
           kf_b = t114 * t132
           ex_unif_b = -0.3e1_dp / 0.4e1_dp * t7 * kf_b
           t135 = 0.1e1_dp / kf_b
           t136 = my_norm_drhob * t135
           t137 = 0.1e1_dp / my_rhob
           s_b = t136 * t137 / 0.2e1_dp
           t139 = s_b ** 2
           t142 = 0.1e1_dp + 0.1243781095e1_dp * mu * t139
           Fx_b = 0.1804e1_dp - 0.804e0_dp / t142
           t145 = my_rhob * ex_unif_b

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)&
                   +t129 * Fx_a + t145 * Fx_b + my_rho * epsilon_cGGA
           END IF

           t149 = my_rho ** 2
           t150 = 0.1e1_dp / t149
           t151 = t1 * t150
           chirhoa = t2 - t151
           t152 = t9 ** 2
           t154 = 0.1e1_dp / t152 * t7
           rsrhoa = -t6 * t154 * t150 / 0.12e2_dp
           t160 = t20 ** 2
           t161 = 0.1e1_dp / t160
           t162 = t12 * t161
           t163 = 0.1e1_dp / t13
           t164 = t163 * rsrhoa
           t167 = t13 * rsrhoa
           t169 = rs ** 0.10e1_dp
           t170 = t169 * rsrhoa
           t172 = 0.3797850000e1_dp * t164 + 0.35876e1_dp * rsrhoa + &
                0.2457300000e1_dp * t167 + 0.985880e0_dp * t170
           t173 = 0.1e1_dp / t23
           t174 = t172 * t173
           e_c_u_0rhoa = -0.1328829340e-1_dp * rsrhoa * t24 + &
                1.e0_dp * t162 * t174
           t179 = t33 ** 2
           t180 = 0.1e1_dp / t179
           t181 = t28 * t180
           t186 = 0.7059450000e1_dp * t164 + 0.61977e1_dp * rsrhoa + &
                0.5049300000e1_dp * t167 + 0.1250340e1_dp * t170
           t187 = 0.1e1_dp / t36
           t188 = t186 * t187
           e_c_u_1rhoa = -0.638837320e-2_dp * rsrhoa * t37 + &
                0.1000000000e1_dp * t181 * t188
           t193 = t46 ** 2
           t194 = 0.1e1_dp / t193
           t195 = t41 * t194
           t200 = 0.5178500000e1_dp * t164 + 0.36231e1_dp * rsrhoa + &
                0.1320390000e1_dp * t167 + 0.993420e0_dp * t170
           t201 = 0.1e1_dp / t49
           t202 = t200 * t201
           alpha_crhoa = 0.375735750e-1_dp * rsrhoa * t50 - 1.e0_dp&
                * t195 * t202
           frhoa = (0.4e1_dp / 0.3e1_dp * t57 * chirhoa - 0.4e1_dp / 0.3e1_dp&
                * t60 * chirhoa) * t55
           t209 = alpha_crhoa * f
           t211 = alpha_c * frhoa
           t213 = t65 * chi
           t214 = t64 * t213
           t215 = t214 * chirhoa
           t217 = 0.4e1_dp * t63 * t215
           t218 = e_c_u_1rhoa - e_c_u_0rhoa
           t219 = t218 * f
           t221 = t70 * frhoa
           t223 = t213 * chirhoa
           t225 = 0.4e1_dp * t71 * t223
           epsilon_c_unifrhoa = e_c_u_0rhoa + t209 * t68 + t211 * t68 - t217 &
                + t219 * t66 + t221 * t66 + t225
           t226 = 0.1e1_dp / t57
           t228 = 0.1e1_dp / t60
           phirhoa = t226 * chirhoa / 0.3e1_dp - t228 * chirhoa / 0.3e1_dp
           t231 = t81 ** 2
           k_frhoa = t3 / t231 * t75 / 0.3e1_dp
           t235 = 0.1e1_dp / t83
           k_srhoa = t235 * k_frhoa * t7
           t237 = 0.1e1_dp / t91
           t238 = my_norm_drho * t237
           t241 = k_s ** 2
           t242 = 0.1e1_dp / t241
           t243 = t242 * t2
           t246 = t86 * t150
           t247 = t85 * t246
           trhoa = -t238 * t87 * phirhoa / 0.2e1_dp - t85 * t243 * k_srhoa / &
                0.2e1_dp - t247 / 0.2e1_dp
           t249 = t96 ** 2
           t251 = t89 / t249
           t252 = epsilon_c_unifrhoa * t89
           t254 = t91 ** 2
           t255 = 0.1e1_dp / t254
           t256 = t255 * phirhoa
           t259 = -t252 * t93 + 0.3e1_dp * t90 * t256
           Arhoa = -0.66725e-1_dp * t251 * t259 * t95
           t263 = gamma_var * t91
           t264 = t113 * phirhoa
           t267 = t89 * t
           t268 = t109 * trhoa
           t271 = Arhoa * t100
           t272 = A * t
           t274 = 0.2e1_dp * t272 * trhoa
           t275 = t271 + t274
           t276 = t275 * t108
           t279 = t107 ** 2
           t280 = 0.1e1_dp / t279
           t281 = t103 * t280
           t282 = A * t105
           t285 = t100 * t
           t286 = t104 * t285
           t289 = t271 + t274 + 0.2e1_dp * t282 * Arhoa + 0.4e1_dp * t286 * trhoa
           t293 = 0.133450e0_dp * t267 * t268 + 0.66725e-1_dp * t101 * t276 -&
                0.66725e-1_dp * t101 * t281 * t289
           t294 = 0.1e1_dp / t112
           t295 = t293 * t294
           epsilon_cGGArhoa = epsilon_c_unifrhoa + 0.3e1_dp * t263 * t264 + &
                t99 * t295
           t298 = t116 ** 2
           kf_arhoa = t114 / t298 * t75 / 0.3e1_dp
           ex_unif_arhoa = -0.3e1_dp / 0.4e1_dp * t7 * kf_arhoa
           t304 = kf_a ** 2
           t305 = 0.1e1_dp / t304
           t306 = my_norm_drhoa * t305
           t309 = my_rhoa ** 2
           t310 = 0.1e1_dp / t309
           s_arhoa = -t306 * t121 * kf_arhoa / 0.2e1_dp - t120 * t310 / 0.2e1_dp
           t313 = t126 ** 2
           t315 = 0.1e1_dp / t313 * mu
           Fx_arhoa = 0.2e1_dp * t315 * s_a * s_arhoa
           t319 = my_rhoa * ex_unif_arhoa

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ra(ii) = e_ra(ii)&
                   +ex_unif_a * Fx_a + t319 * Fx_a + t129 * Fx_arhoa + &
                   epsilon_cGGA + my_rho * epsilon_cGGArhoa
           END IF

           chirhob = -t2 - t151
           rsrhob = rsrhoa
           t326 = t163 * rsrhob
           t329 = t13 * rsrhob
           t331 = t169 * rsrhob
           t333 = 0.3797850000e1_dp * t326 + 0.35876e1_dp * rsrhob + &
                0.2457300000e1_dp * t329 + 0.985880e0_dp * t331
           t334 = t333 * t173
           e_c_u_0rhob = -0.1328829340e-1_dp * rsrhob * t24 + &
                1.e0_dp * t162 * t334
           t343 = 0.7059450000e1_dp * t326 + 0.61977e1_dp * rsrhob + &
                0.5049300000e1_dp * t329 + 0.1250340e1_dp * t331
           t344 = t343 * t187
           e_c_u_1rhob = -0.638837320e-2_dp * rsrhob * t37 + &
                0.1000000000e1_dp * t181 * t344
           t353 = 0.5178500000e1_dp * t326 + 0.36231e1_dp * rsrhob + &
                0.1320390000e1_dp * t329 + 0.993420e0_dp * t331
           t354 = t353 * t201
           alpha_crhob = 0.375735750e-1_dp * rsrhob * t50 - 1.e0_dp&
                * t195 * t354
           frhob = (0.4e1_dp / 0.3e1_dp * t57 * chirhob - 0.4e1_dp / 0.3e1_dp&
                * t60 * chirhob) * t55
           t361 = alpha_crhob * f
           t363 = alpha_c * frhob
           t365 = t214 * chirhob
           t367 = 0.4e1_dp * t63 * t365
           t368 = e_c_u_1rhob - e_c_u_0rhob
           t369 = t368 * f
           t371 = t70 * frhob
           t373 = t213 * chirhob
           t375 = 0.4e1_dp * t71 * t373
           epsilon_c_unifrhob = e_c_u_0rhob + t361 * t68 + t363 * t68 - t367 &
                + t369 * t66 + t371 * t66 + t375
           phirhob = t226 * chirhob / 0.3e1_dp - t228 * chirhob / 0.3e1_dp
           k_frhob = k_frhoa
           k_srhob = t235 * k_frhob * t7
           trhob = -t238 * t87 * phirhob / 0.2e1_dp - t85 * t243 * k_srhob / &
                0.2e1_dp - t247 / 0.2e1_dp
           t385 = epsilon_c_unifrhob * t89
           t387 = t255 * phirhob
           t390 = -t385 * t93 + 0.3e1_dp * t90 * t387
           Arhob = -0.66725e-1_dp * t251 * t390 * t95
           t394 = t113 * phirhob
           t397 = t109 * trhob
           t400 = Arhob * t100
           t402 = 0.2e1_dp * t272 * trhob
           t403 = t400 + t402
           t404 = t403 * t108
           t411 = t400 + t402 + 0.2e1_dp * t282 * Arhob + 0.4e1_dp * t286 * trhob
           t415 = 0.133450e0_dp * t267 * t397 + 0.66725e-1_dp * t101 * t404 -&
                0.66725e-1_dp * t101 * t281 * t411
           t416 = t415 * t294
           epsilon_cGGArhob = epsilon_c_unifrhob + 0.3e1_dp * t263 * t394 + &
                t99 * t416
           t419 = t132 ** 2
           kf_brhob = t114 / t419 * t75 / 0.3e1_dp
           ex_unif_brhob = -0.3e1_dp / 0.4e1_dp * t7 * kf_brhob
           t425 = kf_b ** 2
           t426 = 0.1e1_dp / t425
           t427 = my_norm_drhob * t426
           t430 = my_rhob ** 2
           t431 = 0.1e1_dp / t430
           s_brhob = -t427 * t137 * kf_brhob / 0.2e1_dp - t136 * t431 / 0.2e1_dp
           t434 = t142 ** 2
           t436 = 0.1e1_dp / t434 * mu
           Fx_brhob = 0.2e1_dp * t436 * s_b * s_brhob
           t440 = my_rhob * ex_unif_brhob

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_rb(ii) = e_rb(ii)&
                   +ex_unif_b * Fx_b + t440 * Fx_b + t145 * Fx_brhob + &
                   epsilon_cGGA + my_rho * epsilon_cGGArhob
           END IF

           t445 = t84 * t86
           tnorm_drho = t445 * t2 / 0.2e1_dp
           t450 = t89 * t285
           t451 = A * tnorm_drho
           t459 = 0.2e1_dp * t272 * tnorm_drho + 0.4e1_dp * t286 * tnorm_drho
           t463 = 0.133450e0_dp * t267 * t109 * tnorm_drho + 0.133450e0_dp * &
                t450 * t451 * t108 - 0.66725e-1_dp * t101 * t281 * t459
           t464 = t463 * t294
           Hnorm_drho = t99 * t464

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndr(ii) = e_ndr(ii)&
                   +my_rho * Hnorm_drho
           END IF

           s_anorm_drhoa = t119 * t121 / 0.2e1_dp
           Fx_anorm_drhoa = 0.2e1_dp * t315 * s_a * s_anorm_drhoa

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndra(ii) = e_ndra(ii)&
                   +t129 * Fx_anorm_drhoa
           END IF

           s_bnorm_drhob = t135 * t137 / 0.2e1_dp
           Fx_bnorm_drhob = 0.2e1_dp * t436 * s_b * s_bnorm_drhob

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndrb(ii) = e_ndrb(ii)&
                   +t145 * Fx_bnorm_drhob
           END IF

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              t474 = 0.1e1_dp / t149 / my_rho
              t475 = t1 * t474
              chirhoarhoa = -0.2e1_dp * t150 + 0.2e1_dp * t475
              t480 = t149 ** 2
              rsrhoarhoa = -t6 / t152 / t8 * t76 / t480 / 0.18e2_dp + &
                            t6 * t154 * t474 / 0.6e1_dp
              t490 = rsrhoa * t161
              t495 = t12 / t160 / t20
              t496 = t172 ** 2
              t500 = 0.1e1_dp / t16
              t501 = rsrhoa ** 2
              t502 = t500 * t501
              t504 = t163 * rsrhoarhoa
              t507 = t163 * t501
              t509 = t13 * rsrhoarhoa
              t512 = t169 * rsrhoarhoa
              t518 = t160 ** 2
              t520 = t12 / t518
              t521 = t23 ** 2
              t522 = 0.1e1_dp / t521
              e_c_u_0rhoarhoa = -0.1328829340e-1_dp * rsrhoarhoa * t24 + &
                   0.4274000000e0_dp * t490 * t174 - 0.2000000000e1_dp * t495 * t496 * &
                   t173 + 1.e0_dp * t162 * (-0.1898925000e1_dp * t502 + &
                   0.3797850000e1_dp * t504 + 0.35876e1_dp * rsrhoarhoa + &
                   0.1228650000e1_dp * t507 + 0.2457300000e1_dp * t509 + 0.9858800e0_dp&
                   * t501 + 0.985880e0_dp * t512) * t173 + 0.1608182432e2_dp * t520 * &
                   t496 * t522
              e_c_u_01rhoa = e_c_u_0rhoa
              t528 = rsrhoa * t180
              t533 = t28 / t179 / t33
              t534 = t186 ** 2
              t549 = t179 ** 2
              t551 = t28 / t549
              t552 = t36 ** 2
              t553 = 0.1e1_dp / t552
              t559 = rsrhoa * t194
              t564 = t41 / t193 / t46
              t565 = t200 ** 2
              t580 = t193 ** 2
              t582 = t41 / t580
              t583 = t49 ** 2
              t584 = 0.1e1_dp / t583
              alpha_c1rhoa = alpha_crhoa
              t588 = 0.1e1_dp / t77
              t589 = chirhoa ** 2
              t594 = 0.1e1_dp / t78
              frhoarhoa = (0.4e1_dp / 0.9e1_dp * t588 * t589 + 0.4e1_dp / &
                   0.3e1_dp * t57 * chirhoarhoa + 0.4e1_dp / 0.9e1_dp * t594 * t589 - &
                   0.4e1_dp / 0.3e1_dp * t60 * chirhoarhoa) * t55
              f1rhoa = frhoa
              t612 = alpha_c1rhoa * f
              t615 = alpha_c * f1rhoa
              t618 = t64 * t65
              t633 = e_c_u_1rhoa - e_c_u_01rhoa
              t640 = t633 * f
              t643 = t70 * f1rhoa
              t652 = -0.4e1_dp * t63 * t214 * chirhoarhoa + (-0.638837320e-2_dp &
                   * rsrhoarhoa * t37 + 0.4109600000e0_dp * t528 * t188 - &
                   0.2000000000e1_dp * t533 * t534 * t187 + 0.1000000000e1_dp * t181 * &
                   (-0.3529725000e1_dp * t502 + 0.7059450000e1_dp * t504 + 0.61977e1_dp&
                   * rsrhoarhoa + 0.2524650000e1_dp * t507 + 0.5049300000e1_dp * t509 &
                   + 0.12503400e1_dp * t501 + 0.1250340e1_dp * t512) * t187 + &
                   0.3216468318e2_dp * t551 * t534 * t553 - e_c_u_0rhoarhoa) * f * t66 &
                   + t218 * f1rhoa * t66 + 0.4e1_dp * t219 * t223 + t633 * frhoa * t66 &
                   + t70 * frhoarhoa * t66 + 0.4e1_dp * t221 * t223 + 0.4e1_dp * t640 *&
                   t223 + 0.4e1_dp * t643 * t223 + 0.12e2_dp * t71 * t65 * t589 + &
                   0.4e1_dp * t71 * t213 * chirhoarhoa
              epsilon_c_unifrhoarhoa = e_c_u_0rhoarhoa + (0.375735750e-1_dp * &
                   rsrhoarhoa * t50 - 0.2225000000e0_dp * t559 * t202 + &
                   0.2000000000e1_dp * t564 * t565 * t201 - 1.e0_dp * t195 * &
                   (-0.2589250000e1_dp * t502 + 0.5178500000e1_dp * t504 + 0.36231e1_dp&
                   * rsrhoarhoa + 0.6601950000e0_dp * t507 + 0.1320390000e1_dp * t509 &
                   + 0.9934200e0_dp * t501 + 0.993420e0_dp * t512) * t201 - &
                   0.2960857464e1_dp * t582 * t565 * t584) * f * t68 + alpha_crhoa * &
                   f1rhoa * t68 - 0.4e1_dp * t209 * t215 + alpha_c1rhoa * frhoa * t68 +&
                   alpha_c * frhoarhoa * t68 - 0.4e1_dp * t211 * t215 - 0.4e1_dp * &
                   t612 * t215 - 0.4e1_dp * t615 * t215 - 0.12e2_dp * t63 * t618 * t589&
                   + t652
              epsilon_c_unif1rhoa = e_c_u_01rhoa + t612 * t68 + t615 * t68 - &
                   t217 + t640 * t66 + t643 * t66 + t225
              t657 = 0.1e1_dp / t58
              t662 = 0.1e1_dp / t61
              phirhoarhoa = -t657 * t589 / 0.9e1_dp + t226 * chirhoarhoa / &
                   0.3e1_dp - t662 * t589 / 0.9e1_dp - t228 * chirhoarhoa / 0.3e1_dp
              phi1rhoa = phirhoa
              t670 = t75 ** 2
              k_frhoarhoa = -0.2e1_dp / 0.9e1_dp * t3 / t231 / t80 * t670
              t674 = 0.1e1_dp / t83 / t82
              t675 = k_frhoa ** 2
              k_s1rhoa = k_srhoa
              t682 = my_norm_drho * t93 * t86
              t683 = t2 * phirhoa
              t686 = t238 * t242
              t692 = t238 * t246 * phirhoa / 0.2e1_dp
              t696 = t2 * k_srhoa
              t702 = t85 / t241 / k_s
              t705 = t242 * t150
              t708 = t85 * t705 * k_srhoa / 0.2e1_dp
              t719 = t85 * t86 * t474
              trhoarhoa = t682 * t683 * phi1rhoa + t686 * t683 * k_s1rhoa / &
                   0.2e1_dp + t692 - t238 * t87 * phirhoarhoa / 0.2e1_dp + t686 * t696 &
                   * phi1rhoa / 0.2e1_dp + t702 * t696 * k_s1rhoa + t708 - t85 * t243 *&
                   (-t674 * t675 * t76 / 0.2e1_dp + t235 * k_frhoarhoa * t7) / &
                   0.2e1_dp + t238 * t246 * phi1rhoa / 0.2e1_dp + t85 * t705 * k_s1rhoa&
                   / 0.2e1_dp + t719
              t1rhoa = -t238 * t87 * phi1rhoa / 0.2e1_dp - t85 * t243 * k_s1rhoa&
                   / 0.2e1_dp - t247 / 0.2e1_dp
              t727 = t89 / t249 / t96
              t728 = t95 ** 2
              t729 = t259 * t728
              t730 = epsilon_c_unif1rhoa * t89
              t732 = t255 * phi1rhoa
              t735 = -t730 * t93 + 0.3e1_dp * t90 * t732
              t746 = 0.1e1_dp / t254 / phi
              t747 = t746 * phirhoa
              Arhoarhoa = 0.133450e0_dp * t727 * t729 * t735 - 0.66725e-1_dp * &
                   t251 * (-epsilon_c_unifrhoarhoa * t89 * t93 + 0.3e1_dp * t252 * t732&
                   + 0.3e1_dp * t730 * t256 - 0.12e2_dp * t90 * t747 * phi1rhoa + &
                   0.3e1_dp * t90 * t255 * phirhoarhoa) * t95 - 0.66725e-1_dp * t251 * &
                   t259 * t735 * t95
              A1rhoa = -0.66725e-1_dp * t251 * t735 * t95
              t765 = gamma_var * phi
              t772 = A1rhoa * t100
              t774 = 0.2e1_dp * t272 * t1rhoa
              t775 = t772 + t774
              t776 = t775 * t108
              t783 = t772 + t774 + 0.2e1_dp * t282 * A1rhoa + 0.4e1_dp * t286 * t1rhoa
              t787 = 0.133450e0_dp * t267 * t109 * t1rhoa + 0.66725e-1_dp * t101&
                   * t776 - 0.66725e-1_dp * t101 * t281 * t783
              t788 = t787 * t294
              t804 = t267 * t103
              t805 = t280 * trhoa
              t815 = Arhoarhoa * t100
              t816 = Arhoa * t
              t818 = 0.2e1_dp * t816 * t1rhoa
              t821 = 0.2e1_dp * A1rhoa * t * trhoa
              t824 = 0.2e1_dp * A * t1rhoa * trhoa
              t826 = 0.2e1_dp * t272 * trhoarhoa
              t831 = t275 * t280
              t835 = t280 * t289
              t843 = t101 * t103
              t845 = 0.1e1_dp / t279 / t107
              t846 = t845 * t289
              t853 = A * t285
              t862 = t104 * t100
              t868 = t815 + t818 + t821 + t824 + t826 + 0.2e1_dp * A1rhoa * t105&
                   * Arhoa + 0.8e1_dp * t853 * Arhoa * t1rhoa + 0.2e1_dp * t282 * &
                   Arhoarhoa + 0.8e1_dp * t853 * trhoa * A1rhoa + 0.12e2_dp * t862 * &
                   trhoa * t1rhoa + 0.4e1_dp * t286 * trhoarhoa
              t872 = 0.133450e0_dp * t89 * t1rhoa * t268 + 0.133450e0_dp * t267 &
                   * t776 * trhoa - 0.133450e0_dp * t804 * t805 * t783 + 0.133450e0_dp &
                   * t267 * t109 * trhoarhoa + 0.133450e0_dp * t267 * t276 * t1rhoa + &
                   0.66725e-1_dp * t101 * (t815 + t818 + t821 + t824 + t826) * t108 - &
                   0.66725e-1_dp * t101 * t831 * t783 - 0.133450e0_dp * t804 * t835 * &
                   t1rhoa - 0.66725e-1_dp * t101 * t775 * t280 * t289 + 0.133450e0_dp *&
                   t843 * t846 * t783 - 0.66725e-1_dp * t101 * t281 * t868
              t875 = t112 ** 2
              t876 = 0.1e1_dp / t875
              t877 = t293 * t876
              kf_arhoarhoa = -0.2e1_dp / 0.9e1_dp * t114 / t298 / t115 * t670
              ex_unif_a1rhoa = ex_unif_arhoa
              t895 = kf_arhoa ** 2
              s_a1rhoa = s_arhoa
              t908 = mu ** 2
              t909 = 0.1e1_dp / t313 / t126 * t908
              t910 = t123 * s_arhoa
              Fx_a1rhoa = 0.2e1_dp * t315 * s_a * s_a1rhoa

              e_ra_ra(ii) = e_ra_ra(ii)&
                   +ex_unif_a1rhoa * Fx_a + ex_unif_a * Fx_a1rhoa + &
                   ex_unif_arhoa * Fx_a - 0.3e1_dp / 0.4e1_dp * my_rhoa * t7 * &
                   kf_arhoarhoa * Fx_a + t319 * Fx_a1rhoa + ex_unif_a * Fx_arhoa + my_rhoa&
                   * ex_unif_a1rhoa * Fx_arhoa + t129 * (-0.9950248765e1_dp * t909 * &
                   t910 * s_a1rhoa + 0.2e1_dp * t315 * s_a1rhoa * s_arhoa + &
                   0.2e1_dp * t315 * s_a * (my_norm_drhoa / t304 / kf_a * t121 * &
                   t895 + t306 * t310 * kf_arhoa - t306 * t121 * kf_arhoarhoa / &
                   0.2e1_dp + t120 / t309 / my_rhoa)) + epsilon_c_unif1rhoa + 0.3e1_dp * &
                   t263 * t113 * phi1rhoa + t99 * t788 + epsilon_cGGArhoa + my_rho * (&
                   epsilon_c_unifrhoarhoa + 0.6e1_dp * t765 * t264 * phi1rhoa + &
                   0.3e1_dp * t263 * t788 * phirhoa + 0.3e1_dp * t263 * t113 * &
                   phirhoarhoa + 0.3e1_dp * t263 * t295 * phi1rhoa + t99 * t872 * t294 &
                   - t99 * t877 * t787)

              chirhoarhob = 0.2e1_dp * t475
              rsrhoarhob = rsrhoarhoa
              t938 = rsrhob * t161
              t945 = t500 * rsrhoa * rsrhob
              t947 = t163 * rsrhoarhob
              t950 = t164 * rsrhob
              t952 = t13 * rsrhoarhob
              t954 = rsrhoa * rsrhob
              t956 = t169 * rsrhoarhob
              e_c_u_0rhoarhob = -0.1328829340e-1_dp * rsrhoarhob * t24 + &
                   0.2137000000e0_dp * t490 * t334 + 0.2137000000e0_dp * t938 * t174 - &
                   0.2000000000e1_dp * t495 * t174 * t333 + 1.e0_dp * t162 * &
                   (-0.1898925000e1_dp * t945 + 0.3797850000e1_dp * t947 + 0.35876e1_dp&
                   * rsrhoarhob + 0.1228650000e1_dp * t950 + 0.2457300000e1_dp * t952 &
                   + 0.9858800e0_dp * t954 + 0.985880e0_dp * t956) * t173 + &
                   0.1608182432e2_dp * t520 * t172 * t522 * t333
              t970 = rsrhob * t180
              t995 = rsrhob * t194
              frhoarhob = (0.4e1_dp / 0.9e1_dp * t588 * chirhoa * chirhob + &
                   0.4e1_dp / 0.3e1_dp * t57 * chirhoarhob + 0.4e1_dp / 0.9e1_dp * t594&
                   * chirhoa * chirhob - 0.4e1_dp / 0.3e1_dp * t60 * chirhoarhob) * &
                   t55
              t1045 = t65 * chirhoa * chirhob
              t1074 = -0.4e1_dp * t63 * t214 * chirhoarhob + (-0.638837320e-2_dp&
                   * rsrhoarhob * t37 + 0.2054800000e0_dp * t528 * t344 + &
                   0.2054800000e0_dp * t970 * t188 - 0.2000000000e1_dp * t533 * t188 * &
                   t343 + 0.1000000000e1_dp * t181 * (-0.3529725000e1_dp * t945 + &
                   0.7059450000e1_dp * t947 + 0.61977e1_dp * rsrhoarhob + &
                   0.2524650000e1_dp * t950 + 0.5049300000e1_dp * t952 + &
                   0.12503400e1_dp * t954 + 0.1250340e1_dp * t956) * t187 + &
                   0.3216468318e2_dp * t551 * t186 * t553 * t343 - e_c_u_0rhoarhob) * f&
                   * t66 + t218 * frhob * t66 + 0.4e1_dp * t219 * t373 + t368 * frhoa &
                   * t66 + t70 * frhoarhob * t66 + 0.4e1_dp * t221 * t373 + 0.4e1_dp * &
                   t369 * t223 + 0.4e1_dp * t371 * t223 + 0.12e2_dp * t71 * t1045 + &
                   0.4e1_dp * t71 * t213 * chirhoarhob
              epsilon_c_unifrhoarhob = e_c_u_0rhoarhob + (0.375735750e-1_dp * &
                   rsrhoarhob * t50 - 0.1112500000e0_dp * t559 * t354 - &
                   0.1112500000e0_dp * t995 * t202 + 0.2000000000e1_dp * t564 * t202 * &
                   t353 - 1.e0_dp * t195 * (-0.2589250000e1_dp * t945 + &
                   0.5178500000e1_dp * t947 + 0.36231e1_dp * rsrhoarhob + &
                   0.6601950000e0_dp * t950 + 0.1320390000e1_dp * t952 + 0.9934200e0_dp&
                   * t954 + 0.993420e0_dp * t956) * t201 - 0.2960857464e1_dp * t582 * &
                   t200 * t584 * t353) * f * t68 + alpha_crhoa * frhob * t68 - 0.4e1_dp&
                   * t209 * t365 + alpha_crhob * frhoa * t68 + alpha_c * frhoarhob * &
                   t68 - 0.4e1_dp * t211 * t365 - 0.4e1_dp * t361 * t215 - 0.4e1_dp * &
                   t363 * t215 - 0.12e2_dp * t63 * t64 * t1045 + t1074
              phirhoarhob = -t657 * chirhoa * chirhob / 0.9e1_dp + t226 * &
                   chirhoarhob / 0.3e1_dp - t662 * chirhoa * chirhob / 0.9e1_dp - t228 &
                   * chirhoarhob / 0.3e1_dp
              k_frhoarhob = k_frhoarhoa
              t1109 = t238 * t246 * phirhob / 0.2e1_dp
              t1112 = t85 * t705 * k_srhob / 0.2e1_dp
              trhoarhob = t682 * t683 * phirhob + t686 * t683 * k_srhob / &
                   0.2e1_dp + t692 - t238 * t87 * phirhoarhob / 0.2e1_dp + t686 * t696 &
                   * phirhob / 0.2e1_dp + t702 * t696 * k_srhob + t708 - t85 * t243 * (&
                   -t674 * k_frhoa * t76 * k_frhob / 0.2e1_dp + t235 * k_frhoarhob * t7&
                   ) / 0.2e1_dp + t1109 + t1112 + t719
              Arhoarhob = 0.133450e0_dp * t727 * t729 * t390 - 0.66725e-1_dp * &
                   t251 * (-epsilon_c_unifrhoarhob * t89 * t93 + 0.3e1_dp * t252 * t387&
                   + 0.3e1_dp * t385 * t256 - 0.12e2_dp * t90 * t747 * phirhob + &
                   0.3e1_dp * t90 * t255 * phirhoarhob) * t95 - 0.66725e-1_dp * t251 * &
                   t259 * t390 * t95
              t1163 = Arhoarhob * t100
              t1165 = 0.2e1_dp * t816 * trhob
              t1166 = Arhob * t
              t1168 = 0.2e1_dp * t1166 * trhoa
              t1171 = 0.2e1_dp * A * trhob * trhoa
              t1173 = 0.2e1_dp * t272 * trhoarhob
              t1184 = t403 * t280
              t1207 = t1163 + t1165 + t1168 + t1171 + t1173 + 0.2e1_dp * Arhob *&
                   t105 * Arhoa + 0.8e1_dp * t853 * Arhoa * trhob + 0.2e1_dp * t282 * &
                   Arhoarhob + 0.8e1_dp * t853 * trhoa * Arhob + 0.12e2_dp * t862 * &
                   trhoa * trhob + 0.4e1_dp * t286 * trhoarhob
              t1211 = 0.133450e0_dp * t89 * trhob * t268 + 0.133450e0_dp * t267 &
                   * t404 * trhoa - 0.133450e0_dp * t804 * t805 * t411 + 0.133450e0_dp &
                   * t267 * t109 * trhoarhob + 0.133450e0_dp * t267 * t276 * trhob + &
                   0.66725e-1_dp * t101 * (t1163 + t1165 + t1168 + t1171 + t1173) * &
                   t108 - 0.66725e-1_dp * t101 * t831 * t411 - 0.133450e0_dp * t804 * &
                   t835 * trhob - 0.66725e-1_dp * t101 * t1184 * t289 + 0.133450e0_dp *&
                   t843 * t846 * t411 - 0.66725e-1_dp * t101 * t281 * t1207

              e_ra_rb(ii) = e_ra_rb(ii)&
                   +epsilon_cGGArhob + epsilon_cGGArhoa + my_rho * (&
                   epsilon_c_unifrhoarhob + 0.6e1_dp * t765 * t264 * phirhob + 0.3e1_dp&
                   * t263 * t416 * phirhoa + 0.3e1_dp * t263 * t113 * phirhoarhob + &
                   0.3e1_dp * t263 * t295 * phirhob + t99 * t1211 * t294 - t99 * t877 *&
                   t415)

              chirhobrhob = 0.2e1_dp * t150 + 0.2e1_dp * t475
              rsrhobrhob = rsrhoarhob
              t1222 = t333 ** 2
              t1226 = rsrhob ** 2
              t1227 = t500 * t1226
              t1229 = t163 * rsrhobrhob
              t1232 = t163 * t1226
              t1234 = t13 * rsrhobrhob
              t1237 = t169 * rsrhobrhob
              e_c_u_0rhobrhob = -0.1328829340e-1_dp * rsrhobrhob * t24 + &
                   0.4274000000e0_dp * t938 * t334 - 0.2000000000e1_dp * t495 * t1222 *&
                   t173 + 1.e0_dp * t162 * (-0.1898925000e1_dp * t1227 + &
                   0.3797850000e1_dp * t1229 + 0.35876e1_dp * rsrhobrhob + &
                   0.1228650000e1_dp * t1232 + 0.2457300000e1_dp * t1234 + &
                   0.9858800e0_dp * t1226 + 0.985880e0_dp * t1237) * t173 + &
                   0.1608182432e2_dp * t520 * t1222 * t522
              e_c_u_01rhob = e_c_u_0rhob
              t1250 = t343 ** 2
              t1272 = t353 ** 2
              alpha_c1rhob = alpha_crhob
              t1290 = chirhob ** 2
              frhobrhob = (0.4e1_dp / 0.9e1_dp * t588 * t1290 + 0.4e1_dp / &
                   0.3e1_dp * t57 * chirhobrhob + 0.4e1_dp / 0.9e1_dp * t594 * t1290 - &
                   0.4e1_dp / 0.3e1_dp * t60 * chirhobrhob) * t55
              f1rhob = frhob
              t1312 = alpha_c1rhob * f
              t1315 = alpha_c * f1rhob
              t1332 = e_c_u_1rhob - e_c_u_01rhob
              t1339 = t1332 * f
              t1342 = t70 * f1rhob
              t1351 = -0.4e1_dp * t63 * t214 * chirhobrhob + (-0.638837320e-2_dp&
                   * rsrhobrhob * t37 + 0.4109600000e0_dp * t970 * t344 - &
                   0.2000000000e1_dp * t533 * t1250 * t187 + 0.1000000000e1_dp * t181 *&
                   (-0.3529725000e1_dp * t1227 + 0.7059450000e1_dp * t1229 + &
                   0.61977e1_dp * rsrhobrhob + 0.2524650000e1_dp * t1232 + &
                   0.5049300000e1_dp * t1234 + 0.12503400e1_dp * t1226 + 0.1250340e1_dp&
                   * t1237) * t187 + 0.3216468318e2_dp * t551 * t1250 * t553 - &
                   e_c_u_0rhobrhob) * f * t66 + t368 * f1rhob * t66 + 0.4e1_dp * t369 *&
                   t373 + t1332 * frhob * t66 + t70 * frhobrhob * t66 + 0.4e1_dp * &
                   t371 * t373 + 0.4e1_dp * t1339 * t373 + 0.4e1_dp * t1342 * t373 + &
                   0.12e2_dp * t71 * t65 * t1290 + 0.4e1_dp * t71 * t213 * chirhobrhob
              epsilon_c_unifrhobrhob = e_c_u_0rhobrhob + (0.375735750e-1_dp * &
                   rsrhobrhob * t50 - 0.2225000000e0_dp * t995 * t354 + &
                   0.2000000000e1_dp * t564 * t1272 * t201 - 1.e0_dp * t195 *&
                   (-0.2589250000e1_dp * t1227 + 0.5178500000e1_dp * t1229 + &
                   0.36231e1_dp * rsrhobrhob + 0.6601950000e0_dp * t1232 + &
                   0.1320390000e1_dp * t1234 + 0.9934200e0_dp * t1226 + 0.993420e0_dp *&
                   t1237) * t201 - 0.2960857464e1_dp * t582 * t1272 * t584) * f * t68 &
                   + alpha_crhob * f1rhob * t68 - 0.4e1_dp * t361 * t365 + alpha_c1rhob&
                   * frhob * t68 + alpha_c * frhobrhob * t68 - 0.4e1_dp * t363 * t365 &
                   - 0.4e1_dp * t1312 * t365 - 0.4e1_dp * t1315 * t365 - 0.12e2_dp * &
                   t63 * t618 * t1290 + t1351
              epsilon_c_unif1rhob = e_c_u_01rhob + t1312 * t68 + t1315 * t68 - &
                   t367 + t1339 * t66 + t1342 * t66 + t375
              phirhobrhob = -t657 * t1290 / 0.9e1_dp + t226 * chirhobrhob / &
                   0.3e1_dp - t662 * t1290 / 0.9e1_dp - t228 * chirhobrhob / 0.3e1_dp
              phi1rhob = phirhob
              t1364 = k_frhob ** 2
              k_s1rhob = k_srhob
              t1370 = t2 * phirhob
              t1379 = t2 * k_srhob
              trhobrhob = t682 * t1370 * phi1rhob + t686 * t1370 * k_s1rhob / &
                   0.2e1_dp + t1109 - t238 * t87 * phirhobrhob / 0.2e1_dp + t686 * &
                   t1379 * phi1rhob / 0.2e1_dp + t702 * t1379 * k_s1rhob + t1112 - t85 &
                   * t243 * (-t674 * t1364 * t76 / 0.2e1_dp + t235 * k_frhoarhob * t7) &
                   / 0.2e1_dp + t238 * t246 * phi1rhob / 0.2e1_dp + t85 * t705 * &
                   k_s1rhob / 0.2e1_dp + t719
              t1rhob = -t238 * t87 * phi1rhob / 0.2e1_dp - t85 * t243 * k_s1rhob&
                   / 0.2e1_dp - t247 / 0.2e1_dp
              t1400 = epsilon_c_unif1rhob * t89
              t1402 = t255 * phi1rhob
              t1405 = -t1400 * t93 + 0.3e1_dp * t90 * t1402
              Arhobrhob = 0.133450e0_dp * t727 * t390 * t728 * t1405 - &
                   0.66725e-1_dp * t251 * (-epsilon_c_unifrhobrhob * t89 * t93 + &
                   0.3e1_dp * t385 * t1402 + 0.3e1_dp * t1400 * t387 - 0.12e2_dp * t90 &
                   * t746 * phirhob * phi1rhob + 0.3e1_dp * t90 * t255 * phirhobrhob) *&
                   t95 - 0.66725e-1_dp * t251 * t390 * t1405 * t95
              A1rhob = -0.66725e-1_dp * t251 * t1405 * t95
              t1439 = A1rhob * t100
              t1441 = 0.2e1_dp * t272 * t1rhob
              t1442 = t1439 + t1441
              t1443 = t1442 * t108
              t1450 = t1439 + t1441 + 0.2e1_dp * t282 * A1rhob + 0.4e1_dp * t286&
                   * t1rhob
              t1454 = 0.133450e0_dp * t267 * t109 * t1rhob + 0.66725e-1_dp * &
                   t101 * t1443 - 0.66725e-1_dp * t101 * t281 * t1450
              t1455 = t1454 * t294
              t1471 = t280 * trhob
              t1481 = Arhobrhob * t100
              t1483 = 0.2e1_dp * t1166 * t1rhob
              t1486 = 0.2e1_dp * A1rhob * t * trhob
              t1489 = 0.2e1_dp * A * t1rhob * trhob
              t1491 = 0.2e1_dp * t272 * trhobrhob
              t1499 = t280 * t411
              t1507 = t845 * t411
              t1527 = t1481 + t1483 + t1486 + t1489 + t1491 + 0.2e1_dp * A1rhob &
                   * t105 * Arhob + 0.8e1_dp * t853 * Arhob * t1rhob + 0.2e1_dp * t282 &
                   * Arhobrhob + 0.8e1_dp * t853 * trhob * A1rhob + 0.12e2_dp * t862 * &
                   trhob * t1rhob + 0.4e1_dp * t286 * trhobrhob
              t1531 = 0.133450e0_dp * t89 * t1rhob * t397 + 0.133450e0_dp * t267&
                   * t1443 * trhob - 0.133450e0_dp * t804 * t1471 * t1450 + &
                   0.133450e0_dp * t267 * t109 * trhobrhob + 0.133450e0_dp * t267 * &
                   t404 * t1rhob + 0.66725e-1_dp * t101 * (t1481 + t1483 + t1486 + &
                   t1489 + t1491) * t108 - 0.66725e-1_dp * t101 * t1184 * t1450 - &
                   0.133450e0_dp * t804 * t1499 * t1rhob - 0.66725e-1_dp * t101 * t1442&
                   * t280 * t411 + 0.133450e0_dp * t843 * t1507 * t1450 - &
                   0.66725e-1_dp * t101 * t281 * t1527
              t1534 = t415 * t876
              kf_brhobrhob = -0.2e1_dp / 0.9e1_dp * t114 / t419 / t131 * t670
              ex_unif_b1rhob = ex_unif_brhob
              t1552 = kf_brhob ** 2
              s_b1rhob = s_brhob
              t1565 = 0.1e1_dp / t434 / t142 * t908
              t1566 = t139 * s_brhob
              Fx_b1rhob = 0.2e1_dp * t436 * s_b * s_b1rhob

              e_rb_rb(ii) = e_rb_rb(ii)&
                   +ex_unif_b1rhob * Fx_b + ex_unif_b * Fx_b1rhob + &
                   ex_unif_brhob * Fx_b - 0.3e1_dp / 0.4e1_dp * my_rhob * t7 * &
                   kf_brhobrhob * Fx_b + t440 * Fx_b1rhob + ex_unif_b * Fx_brhob + my_rhob&
                   * ex_unif_b1rhob * Fx_brhob + t145 * (-0.9950248765e1_dp * t1565 * &
                   t1566 * s_b1rhob + 0.2e1_dp * t436 * s_b1rhob * s_brhob + &
                   0.2e1_dp * t436 * s_b * (my_norm_drhob / t425 / kf_b * t137 * &
                   t1552 + t427 * t431 * kf_brhob - t427 * t137 * kf_brhobrhob / &
                   0.2e1_dp + t136 / t430 / my_rhob)) + epsilon_c_unif1rhob + 0.3e1_dp * &
                   t263 * t113 * phi1rhob + t99 * t1455 + epsilon_cGGArhob + my_rho * (&
                   epsilon_c_unifrhobrhob + 0.6e1_dp * t765 * t394 * phi1rhob + &
                   0.3e1_dp * t263 * t1455 * phirhob + 0.3e1_dp * t263 * t113 * &
                   phirhobrhob + 0.3e1_dp * t263 * t416 * phi1rhob + t99 * t1531 * t294&
                   - t99 * t1534 * t1454)

              t1590 = t237 * t86
              t1592 = t84 * t242
              t1594 = t445 * t150
              trhoanorm_drho = -t1590 * t683 / 0.2e1_dp - t1592 * t696 / &
                   0.2e1_dp - t1594 / 0.2e1_dp
              t1599 = t89 * tnorm_drho
              t1602 = t101 * A
              t1603 = tnorm_drho * t108
              t1616 = t816 * tnorm_drho
              t1617 = t451 * trhoa
              t1618 = t272 * trhoanorm_drho
              t1627 = t835 * tnorm_drho
              t1630 = t450 * A
              t1651 = 0.133450e0_dp * t1599 * t268 + 0.266900e0_dp * t1602 * &
                   t1603 * trhoa - 0.133450e0_dp * t804 * t805 * t459 + 0.133450e0_dp *&
                   t267 * t109 * trhoanorm_drho + 0.133450e0_dp * t267 * t276 * &
                   tnorm_drho + 0.66725e-1_dp * t101 * (0.2e1_dp * t1616 + 0.2e1_dp * &
                   t1617 + 0.2e1_dp * t1618) * t108 - 0.66725e-1_dp * t101 * t831 * &
                   t459 - 0.133450e0_dp * t804 * t1627 - 0.133450e0_dp * t1630 * t1627 &
                   + 0.133450e0_dp * t843 * t846 * t459 - 0.66725e-1_dp * t101 * t281 *&
                   (0.2e1_dp * t1616 + 0.2e1_dp * t1617 + 0.2e1_dp * t1618 + 0.8e1_dp &
                   * t853 * Arhoa * tnorm_drho + 0.12e2_dp * t862 * trhoa * tnorm_drho &
                   + 0.4e1_dp * t286 * trhoanorm_drho)

              e_ra_ndr(ii) = e_ra_ndr(ii)&
                   +Hnorm_drho + my_rho * (0.3e1_dp * t263 * t464 * &
                   phirhoa + t99 * t1651 * t294 - t99 * t877 * t463)

              trhobnorm_drho = -t1590 * t1370 / 0.2e1_dp - t1592 * t1379 / &
                   0.2e1_dp - t1594 / 0.2e1_dp
              t1677 = t1166 * tnorm_drho
              t1678 = t451 * trhob
              t1679 = t272 * trhobnorm_drho
              t1688 = t1499 * tnorm_drho
              t1711 = 0.133450e0_dp * t1599 * t397 + 0.266900e0_dp * t1602 * &
                   t1603 * trhob - 0.133450e0_dp * t804 * t1471 * t459 + 0.133450e0_dp &
                   * t267 * t109 * trhobnorm_drho + 0.133450e0_dp * t267 * t404 * &
                   tnorm_drho + 0.66725e-1_dp * t101 * (0.2e1_dp * t1677 + 0.2e1_dp * &
                   t1678 + 0.2e1_dp * t1679) * t108 - 0.66725e-1_dp * t101 * t1184 * &
                   t459 - 0.133450e0_dp * t804 * t1688 - 0.133450e0_dp * t1630 * t1688 &
                   + 0.133450e0_dp * t843 * t1507 * t459 - 0.66725e-1_dp * t101 * t281 &
                   * (0.2e1_dp * t1677 + 0.2e1_dp * t1678 + 0.2e1_dp * t1679 + 0.8e1_dp&
                   * t853 * Arhob * tnorm_drho + 0.12e2_dp * t862 * trhob * tnorm_drho&
                   + 0.4e1_dp * t286 * trhobnorm_drho)

              e_rb_ndr(ii) = e_rb_ndr(ii)&
                   +Hnorm_drho + my_rho * (0.3e1_dp * t263 * t464 * &
                   phirhob + t99 * t1711 * t294 - t99 * t1534 * t463)

              t1717 = tnorm_drho ** 2
              t1721 = A * t1717
              t1726 = t280 * tnorm_drho * t459
              t1732 = t459 ** 2
              t1746 = t463 ** 2

              e_ndr_ndr(ii) = e_ndr_ndr(ii)&
                   +my_rho * (t99 * (0.133450e0_dp * t89 * &
                   t1717 * t109 + 0.667250e0_dp * t101 * t1721 * t108 - 0.266900e0_dp *&
                   t804 * t1726 - 0.266900e0_dp * t1630 * t1726 + 0.133450e0_dp * t101&
                   * t103 * t845 * t1732 - 0.66725e-1_dp * t101 * t281 * (0.2e1_dp * &
                   t1721 + 0.12e2_dp * t862 * t1717)) * t294 - t99 * t1746 * t876)

              e_ra_ndra(ii) = e_ra_ndra(ii)&
                   +ex_unif_a * Fx_anorm_drhoa + t319 * &
                   Fx_anorm_drhoa + t129 * (-0.9950248765e1_dp * t909 * t910 * &
                   s_anorm_drhoa + 0.2e1_dp * t315 * s_anorm_drhoa * s_arhoa +&
                   0.2e1_dp * t315 * s_a * (-t305 * t121 * kf_arhoa / &
                   0.2e1_dp - t119 * t310 / 0.2e1_dp))

              t1766 = s_anorm_drhoa ** 2

              e_ndra_ndra(ii) = e_ndra_ndra(ii)&
                   +t129 * (-0.9950248765e1_dp * t909 * &
                   t123 * t1766 + 0.2e1_dp * t315 * t1766)

              e_rb_ndrb(ii) = e_rb_ndrb(ii)&
                   +ex_unif_b * Fx_bnorm_drhob + t440 * &
                   Fx_bnorm_drhob + t145 * (-0.9950248765e1_dp * t1565 * t1566 * &
                   s_bnorm_drhob + 0.2e1_dp * t436 * s_bnorm_drhob * s_brhob +&
                   0.2e1_dp * t436 * s_b * (-t426 * t137 * kf_brhob / &
                   0.2e1_dp - t135 * t431 / 0.2e1_dp))

              t1790 = s_bnorm_drhob ** 2

              e_ndrb_ndrb(ii) = e_ndrb_ndrb(ii)&
                   +t145 * (-0.9950248765e1_dp * t1565 * &
                   t139 * t1790 + 0.2e1_dp * t436 * t1790)
           END IF
        END IF
     END DO
  END SELECT

END SUBROUTINE pbe_lsd_calc_old
!***************************************************************************


END MODULE xc_pbe
