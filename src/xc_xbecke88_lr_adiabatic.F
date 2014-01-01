!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the density_scaled BECKE88 long-range functional
!>      when used in adiabatic hybrids.
!>      The energy is given as
!>
!>        Ec = 2*lambda*Ex(rho/lambda) + lambda^2*d/dlambda(Ex(rho/lambda)),
!>
!>      where rho/lambda is the scaled density
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE xc_xbecke88_lr_adiabatic
  USE bibliography,                    ONLY: Becke1988,&
                                             cite_reference
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE erf_fn,                          ONLY: erf
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             rootpi
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xbecke88_lr_adiabatic'
  REAL(kind=dp), PARAMETER :: beta=0.0042_dp

  PUBLIC :: xb88_lr_adiabatic_lda_info, xb88_lr_adiabatic_lsd_info, &
            xb88_lr_adiabatic_lda_eval, xb88_lr_adiabatic_lsd_eval
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_adiabatic_lda_info(reference,shortform, needs, max_deriv,error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_adiabatic_lda_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "A. Becke, Phys. Rev. A 38, 3098 (1988) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Becke 1988 Exchange Functional (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xb88_lr_adiabatic_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_adiabatic_lsd_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_adiabatic_lsd_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "A. Becke, Phys. Rev. A 38, 3098 (1988) {LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Becke 1988 Exchange Functional (LSD)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%rho_spin_1_3=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE xb88_lr_adiabatic_lsd_info

! *****************************************************************************
!> \brief evaluates the becke 88 longrange exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param xb 88_params: input parameters (scaling)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_adiabatic_lda_eval(rho_set,deriv_set,grad_deriv,xb88_lr_ad_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xb88_lr_ad_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_adiabatic_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, lambda, omega, sx
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, &
      e_ndrho_rho_rho, e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(xb88_lr_ad_params,"SCALE_X",r_val=sx,error=error)
    CALL section_vals_val_get(xb88_lr_ad_params, "OMEGA",r_val=omega, error=error)
    CALL section_vals_val_get(xb88_lr_ad_params, "LAMBDA",r_val=lambda, error=error)

    CALL cite_reference(Becke1988)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
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
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)
       END IF
       IF (grad_deriv>3.OR.grad_deriv<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       !$omp parallel default(none) &
       !$omp          shared(rho, norm_drho, e_0, e_rho, e_ndrho, e_rho_rho) &
       !$omp          shared(e_ndrho_rho, e_ndrho_ndrho, e_rho_rho_rho) &
       !$omp          shared(e_ndrho_rho_rho, e_ndrho_ndrho_rho) &
       !$omp          shared(e_ndrho_ndrho_ndrho, grad_deriv, npoints) &
       !$omp          shared(epsilon_rho, sx, omega, lambda, error)

       CALL xb88_lr_adiabatic_lda_calc(rho=rho, norm_drho=norm_drho,&
            e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
            e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
            e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
            e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,&
            e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,grad_deriv=grad_deriv,&
            npoints=npoints,epsilon_rho=epsilon_rho,&
            sx=sx,omega=omega,lambda=lambda,error=error)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE xb88_lr_adiabatic_lda_eval

! *****************************************************************************
!> \brief low level calculation of the becke 88 exchange functional for lda
!> \param rho alpha or beta spin density
!> \param norm_drho || grad rho ||
!> \param e_ 0: adds to it the local value of the functional
!> \param e_ *_: derivative of the functional wrt. to the variables
!>        named where the * is. Everything wrt. to the spin of the arguments.
!> \param sx scaling-parameter for exchange
!> \param omega parameter for erfc
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Just took the lsd code and scaled rho and ndrho by 1/2 (e_0 with 2.0)
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_adiabatic_lda_calc(rho, norm_drho,&
       e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
       e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
       e_ndrho_ndrho_ndrho,grad_deriv,npoints,epsilon_rho, sx,&
       omega, lambda, error)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: &
      e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho_rho, e_rho_rho_rho, &
      e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, omega, lambda
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_adiabatic_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: Cx, epsilon_rho43, my_epsilon_rho, my_ndrho, my_rho, t1, &
      t10, t1000, t1002, t1004, t1008, t1012, t1014, t1019, t102, t1021, &
      t1024, t1026, t103, t1030, t1033, t1034, t1036, t104, t1040, t1044, &
      t1046, t105, t1050, t1054, t1055, t1056, t1058, t1059, t106, t1060, &
      t1063, t1064, t1067, t1068, t107, t1071, t1076, t1079, t108, t1082, &
      t1085, t1089, t1093, t1097, t1099, t11, t1100, t1101, t1106, t1109, &
      t1111, t1119, t112, t1122, t1125, t1126, t1129, t1134, t1137, t114, &
      t1140, t1141, t1142, t1143, t115, t1150, t1157, t116, t1161, t1168, &
      t1172, t1176, t1189, t119, t1191, t1192, t1193, t1194, t12, t1200, t1201
    REAL(KIND=dp) :: t1206, t1207, t1210, t1213, t1216, t122, t1221, t1224, &
      t123, t1235, t1238, t1239, t1246, t1247, t1254, t126, t1264, t1268, &
      t1275, t1283, t1284, t1288, t1292, t13, t130, t1300, t1301, t1302, &
      t1309, t131, t1316, t132, t1320, t133, t1333, t1340, t1347, t1348, &
      t136, t1360, t1367, t1368, t137, t1372, t1373, t1377, t1387, t1388, &
      t1392, t1396, t1397, t14, t140, t1401, t1406, t1408, t1409, t141, &
      t1410, t1411, t1413, t1417, t1420, t1425, t1426, t1429, t1436, t1438, &
      t1439, t144, t1440, t1443, t1444, t1447, t145, t1451, t1455, t1459, &
      t146, t1466, t147, t1470, t1471, t1476, t1479, t1482, t1483, t1485, &
      t1489
    REAL(KIND=dp) :: t149, t1492, t1493, t1496, t1498, t15, t150, t1500, &
      t1502, t1504, t1505, t1506, t1507, t1509, t151, t1510, t1512, t1513, &
      t1515, t1517, t1518, t1519, t1523, t1525, t1528, t1530, t1531, t1535, &
      t1536, t1541, t1543, t1545, t1546, t1547, t1549, t155, t1552, t1553, &
      t1556, t1557, t156, t1560, t1563, t1564, t1565, t1568, t1570, t1571, &
      t1574, t1575, t1577, t1579, t1581, t1582, t1584, t1586, t1588, t1589, &
      t1591, t1592, t1594, t1596, t1597, t1598, t1599, t16, t160, t1601, &
      t1602, t1603, t1605, t1607, t1608, t161, t1610, t1611, t1613, t1615, &
      t1616, t1617, t1618, t1620, t1621, t1622, t1624, t1625, t1626, t1628
    REAL(KIND=dp) :: t1630, t1631, t1632, t1634, t1636, t1638, t164, t1640, &
      t1642, t1644, t1645, t1647, t1649, t165, t1652, t1654, t1656, t1658, &
      t166, t1660, t1662, t1664, t1667, t1669, t1671, t1675, t1677, t1679, &
      t1680, t1681, t1684, t1685, t1688, t1689, t1693, t1694, t1696, t1699, &
      t17, t170, t1700, t1702, t1703, t1705, t1706, t1708, t1709, t171, &
      t1712, t1714, t1715, t1717, t1718, t1719, t172, t1720, t1722, t1724, &
      t1726, t1728, t1730, t1731, t1735, t1737, t1740, t1741, t1743, t1745, &
      t1748, t1752, t1754, t1755, t1757, t1758, t176, t1760, t1762, t1763, &
      t1765, t1766, t1768, t1769, t177, t1773, t1774, t1776, t1777, t1779
    REAL(KIND=dp) :: t178, t1782, t1784, t1785, t1787, t1789, t1790, t1794, &
      t1795, t1797, t18, t1800, t1801, t1804, t1807, t1810, t1811, t1814, &
      t1817, t182, t1821, t1822, t1824, t1825, t1828, t1829, t183, t1832, &
      t1835, t1836, t1837, t1838, t1842, t1843, t1849, t1852, t1854, t1858, &
      t1861, t1864, t1865, t1868, t187, t1871, t1874, t1875, t1878, t1881, &
      t1884, t1885, t1886, t189, t1891, t1892, t1895, t1898, t19, t1901, &
      t1902, t1905, t1908, t1912, t1914, t1917, t192, t1920, t1923, t1924, &
      t1925, t1928, t1932, t1938, t1939, t194, t1940, t1943, t1946, t195, &
      t1952, t1954, t1958, t196, t1961, t1962, t1965, t1966, t1969, t1971
    REAL(KIND=dp) :: t1972, t1973, t1974, t1975, t1977, t1980, t1981, t1984, &
      t1987, t1988, t199, t1991, t1994, t1997, t20, t200, t2000, t2001, &
      t2002, t2006, t2009, t2013, t2014, t2018, t2022, t2027, t2028, t203, &
      t2031, t2032, t2035, t2042, t2046, t2050, t2054, t2057, t2058, t2061, &
      t207, t2073, t2074, t2077, t2081, t2082, t2086, t2088, t2095, t2098, &
      t2099, t2104, t2107, t2113, t2114, t2117, t2118, t2124, t2126, t2130, &
      t2131, t2138, t2139, t214, t2143, t2147, t2160, t2169, t217, t218, &
      t2182, t2183, t2187, t219, t2190, t2193, t2196, t2197, t22, t220, &
      t2200, t2201, t2204, t2205, t221, t2220, t2226, t2227, t2230, t224
    REAL(KIND=dp) :: t225, t2257, t2261, t2264, t2268, t2272, t2274, t2276, &
      t2279, t2280, t2281, t2282, t2285, t2287, t2292, t2295, t2296, t2299, &
      t23, t230, t2300, t2304, t2307, t2308, t2311, t2312, t2319, t2326, &
      t233, t234, t2351, t2368, t237, t2378, t2382, t2392, t24, t2404, t2418, &
      t2424, t2437, t2439, t245, t2451, t2455, t2468, t2480, t249, t2492, &
      t25, t250, t2537, t2552, t2565, t2569, t2576, t2585, t26, t260, t2611, &
      t262, t2621, t2622, t2631, t2635, t2636, t264, t2645, t2655, t266, &
      t267, t268, t27, t2707, t271, t2720, t2724, t2732, t274, t275, t2759, &
      t2763, t2779, t278, t2783, t2808, t2816, t282, t2821, t2823
    REAL(KIND=dp) :: t283, t2833, t2843, t2853, t2869, t2873, t2877, t2885, &
      t289, t290, t2935, t2937, t294, t295, t298, t299, t2992, t2994, t2998, &
      t302, t3026, t3030, t306, t3069, t307, t3070, t308, t3085, t309, t3093, &
      t3097, t3098, t31, t311, t3129, t313, t3141, t316, t317, t318, t3181, &
      t321, t322, t3232, t3244, t3248, t326, t3262, t3264, t3268, t327, &
      t3278, t328, t3288, t329, t3291, t3305, t3306, t3314, t332, t3329, &
      t3339, t3344, t3347, t3352, t3353, t3356, t3357, t336, t3362, t3363, &
      t3366, t3378, t3379, t3382, t339, t341, t342, t343, t3449, t3453, t347, &
      t3472, t3491, t3494, t3497, t35, t350, t3501, t3505, t3507
    REAL(KIND=dp) :: t3512, t3513, t3514, t352, t3525, t3533, t3543, t3546, &
      t3547, t3552, t3553, t356, t3561, t3569, t357, t3575, t3576, t3580, &
      t3581, t36, t360, t3608, t3618, t362, t365, t366, t3668, t367, t3673, &
      t3674, t3677, t3679, t368, t3682, t37, t371, t3711, t372, t3726, t375, &
      t3750, t376, t3763, t379, t3794, t3796, t380, t3826, t3827, t383, t384, &
      t3842, t3849, t3853, t3861, t3865, t3869, t387, t3871, t3875, t3892, &
      t3918, t392, t3925, t393, t3933, t3939, t394, t3940, t3965, t3969, &
      t3975, t3979, t398, t3985, t399, t4, t40, t400, t401, t4010, t4018, &
      t4025, t4032, t405, t406, t4070, t4079, t41, t4101, t4108, t413
    REAL(KIND=dp) :: t4140, t4156, t417, t42, t4201, t4214, t4230, t424, &
      t425, t4255, t4256, t4261, t4272, t428, t429, t4292, t43, t4317, t433, &
      t4345, t4352, t437, t4378, t4392, t44, t4428, t444, t445, t446, t4466, &
      t4498, t45, t450, t4505, t451, t4512, t4545, t455, t4553, t4575, t4585, &
      t46, t4607, t462, t4627, t4657, t466, t4696, t47, t470, t4700, t4701, &
      t4704, t4707, t4708, t471, t4713, t4714, t4717, t472, t4721, t4725, &
      t4726, t4729, t473, t4732, t4733, t4734, t4736, t4738, t4739, t474, &
      t4743, t4744, t4748, t4760, t4764, t4768, t4769, t4770, t4771, t4772, &
      t4776, t478, t4782, t4786, t4789, t479, t4790, t4793, t4796
    REAL(KIND=dp) :: t4797, t48, t4800, t4801, t4804, t4808, t4812, t4817, &
      t4820, t4829, t4830, t4833, t4834, t4837, t4838, t4849, t4850, t4853, &
      t4857, t4861, t4864, t4869, t4872, t4874, t4878, t4881, t4884, t4887, &
      t4889, t489, t4893, t4895, t4897, t4898, t49, t4901, t4904, t4905, &
      t4907, t4908, t4909, t4911, t4915, t4917, t4918, t4920, t4922, t4925, &
      t4929, t4932, t4935, t4940, t4941, t4943, t4945, t4946, t4951, t4952, &
      t4955, t4958, t496, t498, t4994, t4995, t50, t500, t5001, t5002, t5005, &
      t5008, t5009, t5011, t5013, t5016, t5019, t5022, t5025, t5028, t5029, &
      t503, t5032, t5035, t5036, t5040, t5043, t5044, t5047, t5050
    REAL(KIND=dp) :: t5051, t5060, t5065, t5069, t5079, t508, t5083, t5087, &
      t5106, t5108, t511, t5115, t5118, t512, t5122, t5123, t5126, t5139, &
      t5143, t5147, t515, t516, t5163, t517, t5170, t5174, t520, t5226, &
      t5228, t5229, t5237, t5244, t5247, t5251, t5256, t526, t5262, t527, &
      t5270, t5271, t5276, t5280, t5281, t5284, t5291, t5294, t5295, t5298, &
      t53, t530, t5301, t5303, t5305, t5308, t531, t5311, t5314, t5317, &
      t5320, t5323, t5326, t5329, t5330, t5333, t5335, t5338, t5339, t534, &
      t5342, t5343, t5345, t5347, t5351, t5354, t5357, t5360, t5363, t5367, &
      t537, t5370, t5371, t5373, t5376, t5379, t5383, t5386, t5389, t5397
    REAL(KIND=dp) :: t5400, t5402, t5406, t5409, t541, t5412, t5415, t5418, &
      t5421, t5423, t5426, t5429, t5430, t5432, t5436, t5440, t5443, t5446, &
      t5449, t5451, t5454, t5456, t5459, t5460, t5462, t5466, t5468, t5469, &
      t5472, t5474, t5478, t5480, t5482, t5484, t5487, t5489, t5492, t5494, &
      t5496, t5499, t55, t550, t5503, t5507, t5509, t551, t5511, t5516, &
      t5523, t5525, t5528, t5531, t5533, t5535, t5537, t5540, t5541, t5545, &
      t5548, t5551, t5553, t5558, t5560, t5562, t5564, t5565, t5568, t557, &
      t5572, t5574, t5575, t5579, t558, t5581, t5585, t5590, t5593, t5595, &
      t5597, t56, t5600, t5602, t5604, t5607, t561, t5610, t5611, t5613
    REAL(KIND=dp) :: t5619, t5622, t5625, t5628, t5631, t5632, t5635, t5638, &
      t5641, t5643, t5647, t5650, t5653, t5657, t566, t5660, t5663, t5666, &
      t5667, t567, t5670, t5671, t5673, t5676, t5679, t5682, t5686, t5689, &
      t569, t5692, t5695, t5696, t5698, t5699, t57, t5701, t5705, t5713, &
      t5716, t5719, t5722, t5724, t5728, t5731, t5734, t5738, t5741, t5742, &
      t5745, t5749, t575, t5752, t5756, t5760, t5767, t5770, t5773, t5776, &
      t5780, t5783, t5786, t5789, t579, t5792, t5795, t5797, t58, t5801, &
      t5804, t5807, t5810, t5813, t5816, t5819, t5822, t5825, t5828, t5830, &
      t5831, t5833, t5835, t5838, t5841, t5845, t5848, t585, t5850, t5854
    REAL(KIND=dp) :: t5857, t5861, t5864, t5868, t5869, t5875, t5880, t5881, &
      t5885, t5888, t589, t5891, t5898, t59, t590, t5902, t5905, t5908, &
      t5911, t5913, t5918, t5922, t5925, t5926, t5929, t593, t5932, t5938, &
      t5940, t5942, t5944, t5949, t5952, t5954, t5956, t5958, t5960, t5963, &
      t5964, t5966, t597, t5971, t5973, t5975, t5977, t5978, t598, t5985, &
      t5990, t5994, t5996, t6, t60, t6000, t6002, t6004, t6006, t6009, t6017, &
      t6019, t6021, t6023, t6025, t6029, t6031, t6033, t6036, t6040, t6042, &
      t6044, t6047, t6049, t605, t6052, t6053, t6056, t6058, t6060, t6063, &
      t6066, t6074, t6075, t6078, t6079, t608, t6087, t609, t6090
    REAL(KIND=dp) :: t6091, t6094, t6099, t61, t6102, t6106, t6109, t6112, &
      t6116, t6119, t6123, t6126, t6129, t6132, t6136, t6138, t614, t6141, &
      t6144, t615, t6155, t6158, t6160, t6162, t6170, t6173, t6175, t6179, &
      t6181, t6185, t6186, t6189, t6192, t6194, t6198, t62, t6201, t6207, &
      t6210, t6213, t6216, t6219, t6222, t6223, t6227, t624, t6247, t6249, &
      t6252, t6253, t6256, t6257, t6260, t6263, t6267, t6268, t627, t6271, &
      t628, t6281, t6283, t6290, t6291, t63, t6308, t6338, t6339, t6344, &
      t635, t6352, t6363, t6367, t637, t6370, t6373, t6375, t6380, t6389, &
      t640, t6400, t6415, t6428, t6433, t6442, t645, t6451, t646, t6461, t6479
    REAL(KIND=dp) :: t649, t6490, t6499, t65, t650, t6509, t6529, t653, &
      t6539, t6554, t6557, t6569, t657, t658, t6583, t6596, t66, t6603, &
      t6611, t6620, t6621, t6630, t6640, t6658, t6665, t6669, t668, t67, &
      t672, t673, t679, t68, t689, t69, t7, t70, t705, t707, t708, t709, t71, &
      t710, t713, t716, t719, t722, t723, t728, t729, t732, t735, t739, t74, &
      t747, t753, t754, t764, t765, t766, t769, t772, t773, t776, t777, t78, &
      t782, t785, t786, t787, t79, t790, t792, t795, t796, t797, t8, t800, &
      t801, t804, t815, t822, t826, t827, t83, t833, t84, t85, t851, t852, &
      t855, t857, t859, t86, t862, t863, t864, t865, t866, t867, t87, t874
    REAL(KIND=dp) :: t877, t878, t88, t881, t882, t886, t889, t893, t894, &
      t897, t9, t905, t906, t907, t914, t918, t92, t922, t923, t926, t929, &
      t93, t930, t94, t940, t941, t942, t944, t946, t949, t95, t950, t953, &
      t954, t955, t958, t959, t962, t963, t966, t969, t97, t972, t973, t974, &
      t975, t976, t979, t98, t980, t983, t984, t987, t99, t991, t995, t997, xx

    failure=.FALSE.
    my_epsilon_rho=epsilon_rho
    epsilon_rho43 = my_epsilon_rho**(4.0_dp/3.0_dp)
    Cx = 1.5_dp*(3.0_dp/4.0_dp/pi)**(1.0_dp/3.0_dp)

    IF (.NOT. failure) THEN
      !$omp do
      DO ii=1,npoints
        !! scale densities with 0.5 because code comes from LSD
        my_rho=rho(ii) * 0.5_dp
        my_ndrho=norm_drho(ii) * 0.5_dp
        IF (my_rho>my_epsilon_rho) THEN
          IF( grad_deriv >=0 ) THEN
            t1 = my_rho ** (0.1e1_dp / 0.3e1_dp)
            xx = my_ndrho * MAX(1.0_dp/t1/my_rho,epsilon_rho43)
            t4 = lambda * my_rho
            t6 = my_ndrho ** 2
            t7 = beta * t6
            t8 = lambda ** 2
            t9 = t8 ** 2
            t10 = t9 ** 2
            t11 = 0.1e1_dp / t10
            t12 = my_rho ** 2
            t13 = t9 * t8
            t14 = 0.1e1_dp / t13
            t15 = t12 * t14
            t16 = t8 * lambda
            t17 = 0.1e1_dp / t16
            t18 = my_rho * t17
            t19 = t18 ** (0.1e1_dp / 0.3e1_dp)
            t20 = t19 ** 2
            t22 = 0.1e1_dp / t20 / t15
            t23 = t11 * t22
            t24 = beta * my_ndrho
            t25 = 0.1e1_dp / t9
            t26 = t19 * t18
            t27 = 0.1e1_dp / t26
            t31 = LOG(my_ndrho * t25 * t27 + SQRT((my_ndrho * t25 * t27) ** 0.2e1_dp +  &
            0.1e1_dp))
            t35 = 0.10e1_dp + 0.60e1_dp * t24 * t25 * t27 * t31
            t36 = 0.1e1_dp / t35
            t37 = t23 * t36
            t40 = 0.20e1_dp * Cx + 0.20e1_dp * t7 * t37
            t41 = SQRT(t40)
            t42 = t41 * t40
            t43 = t4 * t42
            t44 = rootpi
            t45 = 0.1e1_dp / t44
            t46 = omega * t45
            t47 = 0.1e1_dp / omega
            t48 = 0.1e1_dp / t41
            t49 = t47 * t48
            t50 = t44 * t19
            t53 = erf(0.3000000000e1_dp * t49 * t50)
            t55 = omega * t41
            t56 = 0.1e1_dp / t19
            t57 = t45 * t56
            t58 = omega ** 2
            t59 = 0.1e1_dp / t58
            t60 = 0.1e1_dp / t40
            t61 = t59 * t60
            t62 = pi * t20
            t63 = t61 * t62
            t65 = EXP(-0.8999999998e1_dp * t63)
            t66 = t58 * t40
            t67 = 0.1e1_dp / pi
            t68 = 0.1e1_dp / t20
            t69 = t67 * t68
            t70 = t65 - 0.10e1_dp
            t71 = t69 * t70
            t74 = t65 - 0.15e1_dp - 0.5555555558e-1_dp * t66 * t71
            t78 = t44 * t53 + 0.3333333334e0_dp * t55 * t57 * t74
            t79 = t46 * t78
            t83 = my_rho * t41 * omega
            t84 = t45 * t78
            t85 = t10 * lambda
            t86 = 0.1e1_dp / t85
            t87 = t86 * t22
            t88 = t87 * t36
            t92 = 0.1e1_dp / t10 / t9
            t93 = t7 * t92
            t94 = t12 * my_rho
            t95 = t94 * t86
            t97 = 0.1e1_dp / t20 / t95
            t98 = t97 * t36
            t99 = t98 * my_rho
            t102 = t7 * t11
            t103 = t35 ** 2
            t104 = 0.1e1_dp / t103
            t105 = t22 * t104
            t106 = t9 * lambda
            t107 = 0.1e1_dp / t106
            t108 = t107 * t27
            t112 = t24 * t11
            t114 = 0.1e1_dp / t19 / t15
            t115 = t114 * t31
            t116 = t115 * my_rho
            t119 = t24 * t25
            t122 = my_ndrho * t11
            t123 = t114 * my_rho
            t126 = -4 * my_ndrho * t107 * t27 + 4 * t122 * t123
            t130 = 1 + t6 * t11 * t22
            t131 = SQRT(t130)
            t132 = 0.1e1_dp / t131
            t133 = t27 * t126 * t132
            t136 = -0.240e2_dp * t24 * t108 * t31 + 0.240e2_dp * t112 * t116 + 0.60e1_dp  &
            * t119 * t133
            t137 = t105 * t136
            t140 = -0.160e2_dp * t7 * t88 + 0.160e2_dp * t93 * t99 - 0.20e1_dp * t102 *  &
            t137
            t141 = t84 * t140
            t144 = my_rho * t42
            t145 = SQRT(0.3141592654e1_dp)
            t146 = 0.1e1_dp / t145
            t147 = t44 * t146
            t149 = EXP(-0.9000000000e1_dp * t63)
            t150 = 0.1e1_dp / t42
            t151 = t47 * t150
            t155 = t49 * t44
            t156 = t68 * my_rho
            t160 = -0.1500000000e1_dp * t151 * t50 * t140 - 0.3000000000e1_dp * t155  &
            * t156 * t25
            t161 = t149 * t160
            t164 = omega * t48
            t165 = t164 * t45
            t166 = t56 * t74
            t170 = t55 * t45
            t171 = t27 * t74
            t172 = my_rho * t25
            t176 = t40 ** 2
            t177 = 0.1e1_dp / t176
            t178 = t59 * t177
            t182 = t61 * pi
            t183 = t56 * my_rho
            t187 = 0.8999999998e1_dp * t178 * t62 * t140 + 0.1800000000e2_dp * t182 *  &
            t183 * t25
            t189 = t58 * t140
            t192 = t66 * t67
            t194 = 0.1e1_dp / t20 / t18
            t195 = t194 * t70
            t196 = t195 * t172
            t199 = t68 * t187
            t200 = t199 * t65
            t203 = t187 * t65 - 0.5555555558e-1_dp * t189 * t71 - 0.1111111112e0_dp *  &
            t192 * t196 - 0.5555555558e-1_dp * t192 * t200
            t207 = (2 * t147 * t161) + 0.1666666667e0_dp * t165 * t166 * t140  &
            + 0.3333333334e0_dp * t170 * t171 * t172 + 0.3333333334e0_dp * t55 * t57  &
            * t203
            !! Multiply with 2.0 because Code comes from LSD

            e_0(ii) = e_0(ii) +  (-0.4444444448e0_dp * t43 * t79 + t8 * (-0.3333333336e0_dp * t83 * t141  &
            - 0.2222222224e0_dp * t144 * t46 * t207)) * sx * 2.0_dp
          END IF
          IF( grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
            t214 = lambda * t42
            t217 = t4 * t41
            t218 = t10 * t16
            t219 = 0.1e1_dp / t218
            t220 = t219 * t97
            t221 = t220 * t36
            t224 = t9 * t16
            t225 = 0.1e1_dp / t224
            t230 = t220 * t132
            t233 = -0.8000000000e1_dp * t24 * t225 * t114 * t31 - 0.8000000000e1_dp *  &
            t7 * t230
            t234 = t105 * t233
            t237 = -0.5333333333e1_dp * t7 * t221 - 0.20e1_dp * t102 * t234
            t245 = t44 * t68
            t249 = -0.1500000000e1_dp * t151 * t50 * t237 + 0.1000000000e1_dp * t49 *  &
            t245 * t17
            t250 = t149 * t249
            t260 = t178 * t62 * t237
            t262 = pi * t56
            t264 = t61 * t262 * t17
            t266 = 0.8999999998e1_dp * t260 - 0.5999999999e1_dp * t264
            t267 = t266 * t65
            t268 = t58 * t237
            t271 = t195 * t17
            t274 = t68 * t266
            t275 = t274 * t65
            t278 = t267 - 0.5555555558e-1_dp * t268 * t71 + 0.3703703705e-1_dp * t192  &
            * t271 - 0.5555555558e-1_dp * t192 * t275
            t282 = (2 * t147 * t250) + 0.1666666667e0_dp * t165 * t166 * t237  &
            - 0.1111111111e0_dp * t170 * t171 * t17 + 0.3333333334e0_dp * t55 * t57  &
            * t278
            t283 = t46 * t282
            t289 = my_rho * t48 * omega
            t290 = t140 * t237
            t294 = t45 * t282
            t295 = t294 * t140
            t298 = t92 * t97
            t299 = t298 * t36
            t302 = t7 * t86
            t306 = 0.1e1_dp / t10 / t224
            t307 = t7 * t306
            t308 = t12 ** 2
            t309 = t308 * t92
            t311 = 0.1e1_dp / t20 / t309
            t313 = t311 * t36 * my_rho
            t316 = t97 * t104
            t317 = my_rho * t233
            t318 = t316 * t317
            t321 = t7 * t219
            t322 = t316 * t136
            t326 = 0.1e1_dp / t103 / t35
            t327 = t22 * t326
            t328 = t136 * t233
            t329 = t327 * t328
            t332 = t11 * t114
            t336 = t298 * t132
            t339 = t24 * t219
            t341 = 0.1e1_dp / t19 / t95
            t342 = t341 * t31
            t343 = t342 * my_rho
            t347 = t311 * t132 * my_rho
            t350 = t24 * t225
            t352 = t114 * t126 * t132
            t356 = my_ndrho * t219
            t357 = t341 * my_rho
            t360 = 0.28e2_dp / 0.3e1_dp * t122 * t114 - 0.28e2_dp / 0.3e1_dp * t356 * t357
            t362 = t27 * t360 * t132
            t365 = t6 * my_ndrho
            t366 = beta * t365
            t367 = t308 * my_rho
            t368 = 0.1e1_dp / t367
            t371 = 0.1e1_dp / t131 / t130
            t372 = t368 * t126 * t371
            t375 = 0.5600000000e2_dp * t24 * t332 * t31 + 0.3200000000e2_dp * t7 * t336  &
            - 0.5600000000e2_dp * t339 * t343 - 0.3200000000e2_dp * t307 * t347 -  &
            0.8000000000e1_dp * t350 * t352 + 0.60e1_dp * t119 * t362 + 0.8000000000e1_dp  &
            * t366 * t372
            t376 = t105 * t375
            t379 = 0.5866666667e2_dp * t7 * t299 + 0.160e2_dp * t302 * t234 - 0.5866666667e2_dp  &
            * t307 * t313 - 0.160e2_dp * t93 * t318 + 0.5333333333e1_dp * t321  &
            * t322 + 0.40e1_dp * t102 * t329 - 0.20e1_dp * t102 * t376
            t380 = t84 * t379
            t383 = t42 * omega
            t384 = t45 * t207
            t387 = t384 * t237
            t392 = 0.9000000000e1_dp * t260 - 0.6000000000e1_dp * t264
            t393 = t392 * t149
            t394 = t393 * t160
            t398 = 0.1e1_dp / t41 / t176
            t399 = t47 * t398
            t400 = t399 * t44
            t401 = t19 * t140
            t405 = t151 * t44
            t406 = t68 * t140
            t413 = t25 * t237
            t417 = t194 * my_rho
            t424 = 0.2250000000e1_dp * t400 * t401 * t237 - 0.5000000000e0_dp * t405  &
            * t406 * t17 - 0.1500000000e1_dp * t151 * t50 * t379 + 0.1500000000e1_dp  &
            * t405 * t156 * t413 + 0.2000000000e1_dp * t155 * t417 * t225 - 0.3000000000e1_dp  &
            * t49 * t245 * t25
            t425 = t149 * t424
            t428 = omega * t150
            t429 = t428 * t45
            t433 = t140 * t17
            t437 = t56 * t278
            t444 = t45 * t27
            t445 = t164 * t444
            t446 = t74 * my_rho
            t450 = t114 * t74
            t451 = my_rho * t225
            t455 = t27 * t278
            t462 = t56 * t203
            t466 = t27 * t203
            t470 = t176 * t40
            t471 = 0.1e1_dp / t470
            t472 = t59 * t471
            t473 = t472 * pi
            t474 = t20 * t140
            t478 = t178 * pi
            t479 = t56 * t140
            t489 = t27 * my_rho
            t496 = -0.1800000000e2_dp * t473 * t474 * t237 + 0.5999999999e1_dp * t478  &
            * t479 * t17 + 0.8999999998e1_dp * t178 * t62 * t379 - 0.1800000000e2_dp  &
            * t478 * t183 * t413 - 0.6000000000e1_dp * t182 * t489 * t225 + 0.1800000000e2_dp  &
            * t61 * t262 * t25
            t498 = t187 * t266
            t500 = t58 * t379
            t503 = t189 * t67
            t508 = t268 * t67
            t511 = t22 * t70
            t512 = t511 * t451
            t515 = t67 * t194
            t516 = t66 * t515
            t517 = t267 * t172
            t520 = t195 * t25
            t526 = t65 * t17
            t527 = t194 * t187 * t526
            t530 = t68 * t496
            t531 = t530 * t65
            t534 = t199 * t267
            t537 = t496 * t65 + t498 * t65 - 0.5555555558e-1_dp * t500 * t71 + 0.3703703705e-1_dp  &
            * t503 * t271 - 0.5555555558e-1_dp * t503 * t275 - 0.1111111112e0_dp  &
            * t508 * t196 + 0.1851851853e0_dp * t192 * t512 - 0.1111111112e0_dp  &
            * t516 * t517 - 0.1111111112e0_dp * t192 * t520 - 0.5555555558e-1_dp  &
            * t508 * t200 + 0.3703703705e-1_dp * t192 * t527 - 0.5555555558e-1_dp  &
            * t192 * t531 - 0.5555555558e-1_dp * t192 * t534
            t541 = (2 * t147 * t394) + (2 * t147 * t425) - 0.8333333335e-1_dp  &
            * t429 * t166 * t290 - 0.5555555557e-1_dp * t165 * t171 * t433 +  &
            0.1666666667e0_dp * t165 * t437 * t140 + 0.1666666667e0_dp * t165 * t166  &
            * t379 + 0.1666666667e0_dp * t445 * t446 * t413 - 0.4444444445e0_dp * t170  &
            * t450 * t451 + 0.3333333334e0_dp * t170 * t455 * t172 + 0.3333333334e0_dp  &
            * t170 * t171 * t25 + 0.1666666667e0_dp * t165 * t462 * t237 -  &
            0.1111111111e0_dp * t170 * t466 * t17 + 0.3333333334e0_dp * t55 * t57 *  &
            t537
            e_rho(ii) = e_rho(ii) +  (-0.4444444448e0_dp * t214 * t79 - 0.6666666672e0_dp * t217 * t46  &
            * t78 * t237 - 0.4444444448e0_dp * t43 * t283 + t8 * (-0.3333333336e0_dp  &
            * t55 * t141 - 0.1666666668e0_dp * t289 * t84 * t290 - 0.3333333336e0_dp  &
            * t83 * t295 - 0.3333333336e0_dp * t83 * t380 - 0.2222222224e0_dp * t383  &
            * t384 - 0.3333333336e0_dp * t83 * t387 - 0.2222222224e0_dp * t144 *  &
            t46 * t541)) * sx
            t550 = beta * t25
            t551 = t27 * t31
            t557 = 0.60e1_dp * t550 * t551 + 0.60e1_dp * t24 * t23 * t132
            t558 = t105 * t557
            t561 = 0.40e1_dp * t24 * t37 - 0.20e1_dp * t102 * t558
            t566 = pi * t146
            t567 = t566 * t149
            t569 = t151 * t19 * t561
            t575 = t20 * t561
            t579 = t58 * t561
            t585 = 0.8999999998e1_dp * t478 * t575 * t65 - 0.5555555558e-1_dp * t579  &
            * t71 - 0.5000000001e0_dp * t60 * t561 * t65
            t589 = -0.3000000000e1_dp * t567 * t569 + 0.1666666667e0_dp * t165 * t166  &
            * t561 + 0.3333333334e0_dp * t55 * t57 * t585
            t590 = t46 * t589
            t593 = t140 * t561
            t597 = t45 * t589
            t598 = t597 * t140
            t605 = t24 * t92
            t608 = my_rho * t557
            t609 = t316 * t608
            t614 = t136 * t557
            t615 = t327 * t614
            t624 = beta * t11
            t627 = t97 * t132
            t628 = t627 * my_rho
            t635 = -4 * t108 + 4 * t332 * my_rho
            t637 = t27 * t635 * t132
            t640 = 0.1e1_dp / t308
            t645 = -0.240e2_dp * beta * t107 * t551 - 0.240e2_dp * t24 * t87 * t132 +  &
            0.240e2_dp * t624 * t116 + 0.240e2_dp * t605 * t628 + 0.60e1_dp * t550 * t133  &
            + 0.60e1_dp * t119 * t637 - 0.60e1_dp * t7 * t640 * t126 * t371
            t646 = t105 * t645
            t649 = -0.320e2_dp * t24 * t88 + 0.160e2_dp * t302 * t558 + 0.320e2_dp * t605  &
            * t99 - 0.160e2_dp * t93 * t609 - 0.40e1_dp * t112 * t137 + 0.40e1_dp * t102  &
            * t615 - 0.20e1_dp * t102 * t646
            t650 = t84 * t649
            t653 = t384 * t561
            t657 = t44 * pi * t146
            t658 = t657 * t178
            t668 = t25 * t561
            t672 = 0.2250000000e1_dp * t400 * t401 * t561 - 0.1500000000e1_dp * t151  &
            * t50 * t649 + 0.1500000000e1_dp * t405 * t156 * t668
            t673 = t149 * t672
            t679 = t56 * t585
            t689 = t27 * t585
            t705 = -0.1800000000e2_dp * t473 * t474 * t561 + 0.8999999998e1_dp * t178  &
            * t62 * t649 - 0.1800000000e2_dp * t478 * t183 * t668
            t707 = t187 * t59
            t708 = t707 * t177
            t709 = t561 * t65
            t710 = t62 * t709
            t713 = t58 * t649
            t716 = t140 * t177
            t719 = t579 * t67
            t722 = 0.1e1_dp / lambda
            t723 = t60 * t722
            t728 = t68 * t705
            t729 = t728 * t65
            t732 = t60 * t187
            t735 = t705 * t65 + 0.8999999998e1_dp * t708 * t710 - 0.5555555558e-1_dp  &
            * t713 * t71 - 0.5000000001e0_dp * t716 * t709 - 0.1111111112e0_dp * t719  &
            * t196 - 0.1000000001e1_dp * t723 * t709 - 0.5555555558e-1_dp * t719 *  &
            t200 - 0.5555555558e-1_dp * t192 * t729 - 0.5000000001e0_dp * t732 * t709
            t739 = 0.1800000000e2_dp * t658 * t575 * t161 + (2 * t147 * t673)  &
            - 0.8333333335e-1_dp * t429 * t166 * t593 + 0.1666666667e0_dp * t165 * t679  &
            * t140 + 0.1666666667e0_dp * t165 * t166 * t649 + 0.1666666667e0_dp  &
            * t445 * t446 * t668 + 0.3333333334e0_dp * t170 * t689 * t172 + 0.1666666667e0_dp  &
            * t165 * t462 * t561 + 0.3333333334e0_dp * t55 * t57 * t735
            e_ndrho(ii) = e_ndrho(ii) +  (-0.6666666672e0_dp * t217 * t46 * t78 * t561 - 0.4444444448e0_dp  &
            * t43 * t590 + t8 * (-0.1666666668e0_dp * t289 * t84 * t593 - 0.3333333336e0_dp  &
            * t83 * t598 - 0.3333333336e0_dp * t83 * t650 - 0.3333333336e0_dp  &
            * t83 * t653 - 0.2222222224e0_dp * t144 * t46 * t739)) * sx
          END IF

          IF( grad_deriv >= 2 .OR. grad_deriv == -2 ) THEN
            t747 = lambda * t41 * omega
            t753 = t4 * t48
            t754 = t237 ** 2
            t764 = 0.1e1_dp / t10 / t13
            t765 = t764 * t311
            t766 = t765 * t36
            t769 = t316 * t233
            t772 = t233 ** 2
            t773 = t327 * t772
            t776 = t10 * t8
            t777 = 0.1e1_dp / t776
            t782 = t765 * t132
            t785 = t6 ** 2
            t786 = beta * t785
            t787 = t10 ** 2
            t790 = t308 * t94
            t792 = 0.1e1_dp / t787 / t106
            t795 = 0.1e1_dp / t19 / t790 / t792
            t796 = 0.1e1_dp / t787 / t13 * t795
            t797 = t796 * t371
            t800 = 0.1866666667e2_dp * t24 * t777 * t341 * t31 + 0.4000000000e2_dp *  &
            t7 * t782 - 0.1066666667e2_dp * t786 * t797
            t801 = t105 * t800
            t804 = 0.1955555555e2_dp * t7 * t766 + 0.1066666667e2_dp * t321 * t769 +  &
            0.40e1_dp * t102 * t773 - 0.20e1_dp * t102 * t801
            t815 = t68 * t237
            t822 = t44 * t194
            t826 = 0.2250000000e1_dp * t399 * t50 * t754 - 0.1000000000e1_dp * t405 *  &
            t815 * t17 - 0.1500000000e1_dp * t151 * t50 * t804 - 0.6666666667e0_dp  &
            * t49 * t822 * t14
            t827 = t149 * t826
            t833 = t237 * t17
            t851 = 0.1800000000e2_dp * t472 * t62 * t754
            t852 = t56 * t237
            t855 = 0.1200000000e2_dp * t478 * t852 * t17
            t857 = t178 * t62 * t804
            t859 = pi * t27
            t862 = 0.2000000000e1_dp * t61 * t859 * t14
            t863 = -t851 + t855 + 0.8999999998e1_dp * t857 + t862
            t864 = t863 * t65
            t865 = t266 ** 2
            t866 = t865 * t65
            t867 = t58 * t804
            t874 = t511 * t14
            t877 = t194 * t266
            t878 = t877 * t526
            t881 = t68 * t863
            t882 = t881 * t65
            t886 = t68 * t865 * t65
            t889 = t864 + t866 - 0.5555555558e-1_dp * t867 * t71 + 0.7407407410e-1_dp  &
            * t508 * t271 - 0.1111111112e0_dp * t508 * t275 - 0.6172839508e-1_dp *  &
            t192 * t874 + 0.7407407410e-1_dp * t192 * t878 - 0.5555555558e-1_dp * t192  &
            * t882 - 0.5555555558e-1_dp * t192 * t886
            t893 = (2 * t147 * t393 * t249) + (2 * t147 * t827) - 0.8333333335e-1_dp  &
            * t429 * t166 * t754 - 0.1111111111e0_dp * t165 * t171 * t833  &
            + 0.3333333334e0_dp * t165 * t437 * t237 + 0.1666666667e0_dp * t165  &
            * t166 * t804 + 0.1481481481e0_dp * t170 * t450 * t14 - 0.2222222222e0_dp  &
            * t170 * t455 * t17 + 0.3333333334e0_dp * t55 * t57 * t889
            t894 = t46 * t893
            t897 = t78 * t140
            t905 = my_rho * t150
            t906 = t905 * omega
            t907 = t140 * t754
            t914 = t379 * t237
            t918 = t140 * t804
            t922 = t45 * t893
            t923 = t922 * t140
            t926 = t294 * t379
            t929 = t306 * t311
            t930 = t929 * t36
            t940 = 0.1e1_dp / t787 / t8
            t941 = t7 * t940
            t942 = t367 * t306
            t944 = 0.1e1_dp / t20 / t942
            t946 = t944 * t36 * my_rho
            t949 = t311 * t104
            t950 = t949 * t317
            t953 = t97 * t326
            t954 = my_rho * t772
            t955 = t953 * t954
            t958 = my_rho * t800
            t959 = t316 * t958
            t962 = t7 * t764
            t963 = t949 * t136
            t966 = t953 * t328
            t969 = t316 * t375
            t972 = t103 ** 2
            t973 = 0.1e1_dp / t972
            t974 = t22 * t973
            t975 = t136 * t772
            t976 = t974 * t975
            t979 = t375 * t233
            t980 = t327 * t979
            t983 = t136 * t800
            t984 = t327 * t983
            t987 = t219 * t341
            t991 = t929 * t132
            t995 = 0.1e1_dp / t787 / t224
            t997 = t995 * t795 * t371
            t1000 = t24 * t764
            t1002 = 0.1e1_dp / t19 / t309
            t1004 = t1002 * t31 * my_rho
            t1008 = t944 * t132 * my_rho
            t1012 = 0.1e1_dp / t787 / t776
            t1014 = t308 ** 2
            t1019 = 0.1e1_dp / t19 / t1014 * t787 * t10
            t1021 = t1019 * t371 * my_rho
            t1024 = t24 * t777
            t1026 = t341 * t126 * t132
            t1030 = t114 * t360 * t132
            t1033 = t308 * t12
            t1034 = 0.1e1_dp / t1033
            t1036 = t1034 * t126 * t371
            t1040 = my_ndrho * t764
            t1044 = -0.280e3_dp / 0.9e1_dp * t356 * t341 + 0.280e3_dp / 0.9e1_dp * t1040 *  &
            t1002 * my_rho
            t1046 = t27 * t1044 * t132
            t1050 = t368 * t360 * t371
            t1054 = beta * t785 * my_ndrho
            t1055 = t1054 * t368
            t1056 = t130 ** 2
            t1058 = 0.1e1_dp / t131 / t1056
            t1059 = t126 * t1058
            t1060 = t1059 * t220
            t1063 = -0.1866666667e3_dp * t24 * t987 * t31 - 0.2240000000e3_dp * t7 *  &
            t991 + 0.4266666667e2_dp * t786 * t997 + 0.1866666667e3_dp * t1000 * t1004  &
            + 0.2240000000e3_dp * t941 * t1008 - 0.4266666667e2_dp * t786 * t1012  &
            * t1021 + 0.1866666667e2_dp * t1024 * t1026 - 0.1600000000e2_dp * t350  &
            * t1030 - 0.5066666667e2_dp * t366 * t1036 + 0.60e1_dp * t119 * t1046 +  &
            0.1600000000e2_dp * t366 * t1050 + 0.3200000000e2_dp * t1055 * t1060
            t1064 = t105 * t1063
            t1067 = -0.2737777778e3_dp * t7 * t930 - 0.1173333333e3_dp * t93 * t769 -  &
            0.320e2_dp * t302 * t773 + 0.160e2_dp * t302 * t801 + 0.2737777778e3_dp *  &
            t941 * t946 + 0.1173333333e3_dp * t307 * t950 + 0.320e2_dp * t93 * t955  &
            - 0.160e2_dp * t93 * t959 - 0.1955555555e2_dp * t962 * t963 - 0.2133333334e2_dp  &
            * t321 * t966 + 0.1066666667e2_dp * t321 * t969 - 0.120e2_dp * t102  &
            * t976 + 0.80e1_dp * t102 * t980 + 0.40e1_dp * t102 * t984 - 0.20e1_dp * t102  &
            * t1064
            t1068 = t84 * t1067
            t1071 = t45 * t541
            t1076 = t384 * t754
            t1079 = t1071 * t237
            t1082 = t384 * t804
            t1085 = t379 * t17
            t1089 = t114 * t278
            t1093 = t25 * t804
            t1097 = t187 * t865
            t1099 = t97 * t70
            t1100 = my_rho * t777
            t1101 = t1099 * t1100
            t1106 = t500 * t67
            t1109 = t187 * t863
            t1111 = t496 * t266
            t1119 = t194 * t496 * t526
            t1122 = t199 * t866
            t1125 = t65 * t25
            t1126 = t877 * t1125
            t1129 = t867 * t67
            t1134 = t866 * t172
            t1137 = t268 * t515
            t1140 = t176 ** 2
            t1141 = 0.1e1_dp / t1140
            t1142 = t59 * t1141
            t1143 = t1142 * pi
            t1150 = t20 * t379
            t1157 = t27 * t140
            t1161 = t56 * t379
            t1168 = t25 * t754
            t1172 = t225 * t237
            t1176 = t56 * t25
            t1189 = 0.5400000000e2_dp * t1143 * t474 * t754 - 0.2400000000e2_dp * t473  &
            * t479 * t833 - 0.3600000000e2_dp * t473 * t1150 * t237 - 0.1800000000e2_dp  &
            * t473 * t474 * t804 - 0.2000000000e1_dp * t478 * t1157 * t14 +  &
            0.1200000000e2_dp * t478 * t1161 * t17 + 0.8999999998e1_dp * t178 * t62  &
            * t1067 + 0.3600000000e2_dp * t473 * t183 * t1168 + 0.1200000000e2_dp *  &
            t478 * t489 * t1172 - 0.3600000000e2_dp * t478 * t1176 * t237 - 0.1800000000e2_dp  &
            * t478 * t183 * t1093 + 0.8000000000e1_dp * t182 * t123 *  &
            t777 - 0.1200000000e2_dp * t61 * t859 * t225
            t1191 = (t1097 * t65) - 0.4938271608e0_dp * t192 * t1101 - 0.6172839508e-1_dp  &
            * t503 * t874 - 0.1111111112e0_dp * t1106 * t275 + (t1109  &
            * t65) + (2 * t1111 * t65) + 0.3703703706e0_dp * t508 * t512 +  &
            0.7407407410e-1_dp * t503 * t878 + 0.7407407410e-1_dp * t192 * t1119 - 0.5555555558e-1_dp  &
            * t192 * t1122 - 0.2222222224e0_dp * t192 * t1126 - 0.1111111112e0_dp  &
            * t1129 * t196 - 0.5555555558e-1_dp * t503 * t882 - 0.1111111112e0_dp  &
            * t516 * t1134 - 0.2222222224e0_dp * t1137 * t517 + (t1189  &
            * t65)
            t1192 = t67 * t22
            t1193 = t66 * t1192
            t1194 = t267 * t451
            t1200 = t65 * t14
            t1201 = t22 * t187 * t1200
            t1206 = t68 * t1189
            t1207 = t1206 * t65
            t1210 = t58 * t1067
            t1213 = t511 * t225
            t1216 = t864 * t172
            t1221 = t498 * t526
            t1224 = t530 * t267
            t1235 = t199 * t864
            t1238 = 0.3703703706e0_dp * t1193 * t1194 - 0.5555555558e-1_dp * t503 * t886  &
            - 0.6172839508e-1_dp * t192 * t1201 - 0.5555555558e-1_dp * t1129 * t200  &
            - 0.5555555558e-1_dp * t192 * t1207 - 0.5555555558e-1_dp * t1210 * t71  &
            + 0.3703703706e0_dp * t192 * t1213 - 0.1111111112e0_dp * t516 * t1216  &
            - 0.1111111112e0_dp * t508 * t534 + 0.7407407410e-1_dp * t516 * t1221 -  &
            0.1111111112e0_dp * t192 * t1224 + 0.7407407410e-1_dp * t508 * t527 + 0.7407407410e-1_dp  &
            * t1106 * t271 - 0.1111111112e0_dp * t508 * t531 - 0.2222222224e0_dp  &
            * t508 * t520 - 0.5555555558e-1_dp * t192 * t1235
            t1239 = t1191 + t1238
            t1246 = t428 * t444
            t1247 = t74 * t140
            t1254 = t27 * t889
            t1264 = t27 * t537
            t1268 = t341 * t74
            t1275 = t56 * t889
            t1283 = -t851 + t855 + 0.9000000000e1_dp * t857 + t862
            t1284 = t1283 * t149
            t1288 = -0.1111111111e0_dp * t165 * t171 * t1085 - 0.8888888890e0_dp * t170  &
            * t1089 * t451 + 0.1666666667e0_dp * t445 * t446 * t1093 + 0.3333333334e0_dp  &
            * t55 * t57 * t1239 - 0.8333333335e-1_dp * t429 * t462 * t754  &
            + 0.5555555556e-1_dp * t1246 * t1247 * t833 + (4 * t147 * t393 *  &
            t424) + 0.3333333334e0_dp * t170 * t1254 * t172 - 0.1111111111e0_dp * t165  &
            * t466 * t833 + 0.3333333334e0_dp * t165 * t437 * t379 - 0.2222222222e0_dp  &
            * t170 * t1264 * t17 + 0.1037037037e1_dp * t170 * t1268 * t1100  &
            + 0.1666666667e0_dp * t165 * t166 * t1067 + 0.1666666667e0_dp * t165 *  &
            t1275 * t140 - 0.8333333335e-1_dp * t1246 * t446 * t1168 + (2 *  &
            t147 * t1284 * t160)
            t1292 = t278 * my_rho
            t1300 = 0.1e1_dp / t41 / t470
            t1301 = t47 * t1300
            t1302 = t1301 * t44
            t1309 = t19 * t379
            t1316 = t194 * t140
            t1320 = t68 * t379
            t1333 = t68 * t25
            t1340 = t22 * my_rho
            t1347 = -0.5625000000e1_dp * t1302 * t401 * t754 + 0.1500000000e1_dp * t400  &
            * t406 * t833 + 0.4500000000e1_dp * t400 * t1309 * t237 + 0.2250000000e1_dp  &
            * t400 * t401 * t804 + 0.3333333333e0_dp * t405 * t1316 * t14  &
            - 0.1000000000e1_dp * t405 * t1320 * t17 - 0.1500000000e1_dp * t151 * t50  &
            * t1067 - 0.2250000000e1_dp * t400 * t156 * t1168 - 0.2000000000e1_dp  &
            * t405 * t417 * t1172 + 0.3000000000e1_dp * t405 * t1333 * t237 + 0.1500000000e1_dp  &
            * t405 * t156 * t1093 - 0.3333333333e1_dp * t155 * t1340  &
            * t777 + 0.4000000000e1_dp * t49 * t822 * t225
            t1348 = t149 * t1347
            t1360 = t140 * t14
            t1367 = omega * t398
            t1368 = t1367 * t45
            t1372 = t45 * t114
            t1373 = t164 * t1372
            t1377 = t114 * t203
            t1387 = t392 ** 2
            t1388 = t1387 * t149
            t1392 = t56 * t537
            t1396 = 0.6666666668e0_dp * t170 * t455 * t25 + 0.3333333334e0_dp * t445  &
            * t1292 * t413 + 0.1666666667e0_dp * t165 * t462 * t804 + (2 * t147  &
            * t1348) - 0.1111111111e0_dp * t165 * t455 * t433 - 0.8333333335e-1_dp  &
            * t429 * t166 * t918 + 0.3333333334e0_dp * t165 * t171 * t413 + 0.7407407409e-1_dp  &
            * t165 * t450 * t1360 - 0.1666666667e0_dp * t429 * t437  &
            * t290 + 0.1250000000e0_dp * t1368 * t166 * t907 - 0.4444444445e0_dp * t1373  &
            * t446 * t1172 + 0.1481481481e0_dp * t170 * t1377 * t14 - 0.8888888890e0_dp  &
            * t170 * t450 * t225 - 0.1666666667e0_dp * t429 * t166 * t914  &
            + (2 * t147 * t1388 * t160) + 0.3333333334e0_dp * t165 * t1392  &
            * t237
            t1397 = t1288 + t1396
            t1401 = -0.3333333336e0_dp * t165 * t897 * t237 - 0.6666666672e0_dp * t55  &
            * t295 - 0.6666666672e0_dp * t55 * t380 + 0.8333333340e-1_dp * t906 * t84  &
            * t907 - 0.3333333336e0_dp * t289 * t294 * t290 - 0.3333333336e0_dp *  &
            t289 * t84 * t914 - 0.1666666668e0_dp * t289 * t84 * t918 - 0.3333333336e0_dp  &
            * t83 * t923 - 0.6666666672e0_dp * t83 * t926 - 0.3333333336e0_dp  &
            * t83 * t1068 - 0.4444444448e0_dp * t383 * t1071 - 0.6666666672e0_dp *  &
            t55 * t387 - 0.1666666668e0_dp * t289 * t1076 - 0.6666666672e0_dp * t83  &
            * t1079 - 0.3333333336e0_dp * t83 * t1082 - 0.2222222224e0_dp * t144 * t46  &
            * t1397
            e_rho_rho(ii) = e_rho_rho(ii) +  (-0.1333333334e1_dp * t747 * t84 * t237 - 0.8888888896e0_dp * t214  &
            * t283 - 0.3333333336e0_dp * t753 * t46 * t78 * t754 - 0.1333333334e1_dp  &
            * t217 * t46 * t282 * t237 - 0.6666666672e0_dp * t217 * t46 * t78 *  &
            t804 - 0.4444444448e0_dp * t43 * t894 + t8 * t1401) * sx
            t1406 = 0.6666666672e0_dp * t747 * t84 * t561
            t1408 = 0.4444444448e0_dp * t214 * t590
            t1409 = t4 * t164
            t1410 = t237 * t561
            t1411 = t84 * t1410
            t1413 = 0.3333333336e0_dp * t1409 * t1411
            t1417 = 0.6666666672e0_dp * t217 * t46 * t589 * t237
            t1420 = t316 * t557
            t1425 = t233 * t557
            t1426 = t327 * t1425
            t1429 = beta * t225
            t1436 = t1033 * t940
            t1438 = 0.1e1_dp / t19 / t1436
            t1439 = 0.1e1_dp / t787 / t16 * t1438
            t1440 = t1439 * t371
            t1443 = -0.8000000000e1_dp * t1429 * t115 - 0.2400000000e2_dp * t24 * t230  &
            + 0.8000000000e1_dp * t366 * t1440
            t1444 = t105 * t1443
            t1447 = -0.1066666667e2_dp * t24 * t221 + 0.5333333333e1_dp * t321 * t1420  &
            - 0.40e1_dp * t112 * t234 + 0.40e1_dp * t102 * t1426 - 0.20e1_dp * t102 *  &
            t1444
            t1451 = 0.6666666672e0_dp * t217 * t46 * t78 * t1447
            t1455 = 0.6666666672e0_dp * t217 * t46 * t282 * t561
            t1459 = t19 * t237
            t1466 = t68 * t17
            t1470 = 0.2250000000e1_dp * t400 * t1459 * t561 - 0.1500000000e1_dp * t151  &
            * t50 * t1447 - 0.5000000000e0_dp * t405 * t1466 * t561
            t1471 = t149 * t1470
            t1476 = 0.8333333335e-1_dp * t429 * t166 * t1410
            t1479 = 0.1666666667e0_dp * t165 * t679 * t237
            t1482 = 0.1666666667e0_dp * t165 * t166 * t1447
            t1483 = t17 * t561
            t1485 = t165 * t171 * t1483
            t1489 = 0.1111111111e0_dp * t170 * t689 * t17
            t1492 = 0.1666666667e0_dp * t165 * t437 * t561
            t1493 = t20 * t237
            t1496 = 0.1800000000e2_dp * t473 * t1493 * t561
            t1498 = t178 * t62 * t1447
            t1500 = t56 * t17
            t1502 = t478 * t1500 * t561
            t1504 = -t1496 + 0.8999999998e1_dp * t1498 + 0.5999999999e1_dp * t1502
            t1505 = t1504 * t65
            t1506 = t266 * t59
            t1507 = t1506 * t177
            t1509 = 0.8999999998e1_dp * t1507 * t710
            t1510 = t58 * t1447
            t1512 = 0.5555555558e-1_dp * t1510 * t71
            t1513 = t237 * t177
            t1515 = 0.5000000001e0_dp * t1513 * t709
            t1517 = 0.3703703705e-1_dp * t719 * t271
            t1518 = 0.1e1_dp / my_rho
            t1519 = t60 * t1518
            t1523 = 0.5555555558e-1_dp * t719 * t275
            t1525 = t68 * t1504 * t65
            t1528 = t60 * t266
            t1530 = 0.5000000001e0_dp * t1528 * t709
            t1531 = t1505 + t1509 - t1512 - t1515 + t1517 + 0.3333333334e0_dp * t1519  &
            * t709 - t1523 - 0.5555555558e-1_dp * t192 * t1525 - t1530
            t1535 = 0.1800000000e2_dp * t658 * t575 * t250 + (2 * t147 * t1471)  &
            - t1476 + t1479 + t1482 - 0.5555555555e-1_dp * t1485 - t1489 + t1492  &
            + 0.3333333334e0_dp * t55 * t57 * t1531
            t1536 = t46 * t1535
            t1541 = 0.1666666668e0_dp * t165 * t897 * t561
            t1543 = 0.3333333336e0_dp * t55 * t598
            t1545 = 0.3333333336e0_dp * t55 * t650
            t1546 = t905 * t46
            t1547 = t897 * t1410
            t1549 = 0.8333333340e-1_dp * t1546 * t1547
            t1552 = 0.1666666668e0_dp * t289 * t597 * t290
            t1553 = t649 * t237
            t1556 = 0.1666666668e0_dp * t289 * t84 * t1553
            t1557 = t140 * t1447
            t1560 = 0.1666666668e0_dp * t289 * t84 * t1557
            t1563 = 0.1666666668e0_dp * t289 * t294 * t593
            t1564 = t45 * t1535
            t1565 = t1564 * t140
            t1568 = t294 * t649
            t1570 = 0.3333333336e0_dp * t83 * t1568
            t1571 = t379 * t561
            t1574 = 0.1666666668e0_dp * t289 * t84 * t1571
            t1575 = t597 * t379
            t1577 = 0.3333333336e0_dp * t83 * t1575
            t1579 = 0.1173333333e3_dp * t24 * t299
            t1581 = 0.5866666667e2_dp * t93 * t1420
            t1582 = t24 * t86
            t1584 = 0.320e2_dp * t1582 * t234
            t1586 = 0.320e2_dp * t302 * t1426
            t1588 = 0.160e2_dp * t302 * t1444
            t1589 = t24 * t306
            t1591 = 0.1173333333e3_dp * t1589 * t313
            t1592 = t949 * t608
            t1594 = 0.5866666667e2_dp * t307 * t1592
            t1596 = 0.320e2_dp * t605 * t318
            t1597 = t7 * t298
            t1598 = t326 * my_rho
            t1599 = t1598 * t1425
            t1601 = 0.320e2_dp * t1597 * t1599
            t1602 = my_rho * t1443
            t1603 = t316 * t1602
            t1605 = 0.160e2_dp * t93 * t1603
            t1607 = 0.1066666667e2_dp * t339 * t322
            t1608 = t953 * t614
            t1610 = 0.1066666667e2_dp * t321 * t1608
            t1611 = t316 * t645
            t1613 = 0.5333333333e1_dp * t321 * t1611
            t1615 = 0.80e1_dp * t112 * t329
            t1616 = t7 * t23
            t1617 = t973 * t136
            t1618 = t1617 * t1425
            t1620 = 0.120e2_dp * t1616 * t1618
            t1621 = t645 * t233
            t1622 = t327 * t1621
            t1624 = 0.40e1_dp * t102 * t1622
            t1625 = t136 * t1443
            t1626 = t327 * t1625
            t1628 = 0.40e1_dp * t102 * t1626
            t1630 = 0.40e1_dp * t112 * t376
            t1631 = t375 * t557
            t1632 = t327 * t1631
            t1634 = 0.40e1_dp * t102 * t1632
            t1636 = 0.5600000000e2_dp * t624 * t115
            t1638 = 0.1200000000e3_dp * t24 * t336
            t1640 = 0.1e1_dp / t787 / t9
            t1642 = t1640 * t1438 * t371
            t1644 = 0.3200000000e2_dp * t366 * t1642
            t1645 = beta * t219
            t1647 = 0.5600000000e2_dp * t1645 * t343
            t1649 = 0.1200000000e3_dp * t1589 * t347
            t1652 = t795 * t371 * my_rho
            t1654 = 0.3200000000e2_dp * t366 * t995 * t1652
            t1656 = 0.8000000000e1_dp * t1429 * t352
            t1658 = t114 * t635 * t132
            t1660 = 0.8000000000e1_dp * t350 * t1658
            t1662 = 0.3200000000e2_dp * t7 * t372
            t1664 = 0.60e1_dp * t550 * t362
            t1667 = 0.28e2_dp / 0.3e1_dp * t332 - 0.28e2_dp / 0.3e1_dp * t987 * my_rho
            t1669 = t27 * t1667 * t132
            t1671 = 0.60e1_dp * t119 * t1669
            t1675 = 0.60e1_dp * t7 * t640 * t360 * t371
            t1677 = t368 * t635 * t371
            t1679 = 0.8000000000e1_dp * t366 * t1677
            t1680 = t786 * t368
            t1681 = t1059 * t23
            t1684 = t1636 + t1638 - t1644 - t1647 - t1649 + t1654 - t1656 - t1660  &
            + t1662 + t1664 + t1671 - t1675 + t1679 - 0.2400000000e2_dp * t1680  &
            * t1681
            t1685 = t105 * t1684
            t1688 = t1579 - t1581 + t1584 - t1586 + t1588 - t1591 + t1594 - t1596  &
            + t1601 - t1605 + t1607 - t1610 + t1613 + t1615 - t1620 + t1624  &
            + t1628 - t1630 + t1634 - 0.20e1_dp * t102 * t1685
            t1689 = t84 * t1688
            t1693 = 0.3333333336e0_dp * t55 * t653
            t1694 = t45 * t739
            t1696 = 0.2222222224e0_dp * t383 * t1694
            t1699 = 0.1666666668e0_dp * t289 * t384 * t1410
            t1700 = t1694 * t237
            t1702 = 0.3333333336e0_dp * t83 * t1700
            t1703 = t384 * t1447
            t1705 = 0.3333333336e0_dp * t83 * t1703
            t1706 = t1071 * t561
            t1708 = 0.3333333336e0_dp * t83 * t1706
            t1709 = t187 * t1504
            t1712 = t266 * t561 * t65
            t1714 = 0.5000000001e0_dp * t732 * t1712
            t1715 = t713 * t67
            t1717 = 0.3703703705e-1_dp * t1715 * t271
            t1718 = t705 * t266
            t1719 = t1718 * t65
            t1720 = t140 * t266
            t1722 = t177 * t561 * t65
            t1724 = 0.5000000001e0_dp * t1720 * t1722
            t1726 = t194 * t705 * t526
            t1728 = 0.3703703705e-1_dp * t192 * t1726
            t1730 = 0.1111111112e0_dp * t719 * t520
            t1731 = t58 * t1688
            t1735 = 0.3703703705e-1_dp * t719 * t527
            t1737 = t722 * t561 * t65
            t1740 = t496 * t59
            t1741 = t1740 * t177
            t1743 = 0.8999999998e1_dp * t1741 * t710
            t1745 = t187 * t561 * t65
            t1748 = t140 * t1518
            t1752 = 0.1851851853e0_dp * t719 * t512
            t1754 = 0.5555555558e-1_dp * t1715 * t275
            t1755 = t1510 * t67
            t1757 = 0.1111111112e0_dp * t1755 * t196
            t1758 = t1709 * t65 - t1714 + t1717 + t1719 - t1724 + t1728 - t1730  &
            - 0.5555555558e-1_dp * t1731 * t71 + t1735 + 0.666666666e0_dp * t1519 *  &
            t1737 + t1743 + 0.3333333334e0_dp * t1519 * t1745 + 0.3333333334e0_dp *  &
            t1748 * t1722 + t1752 - t1754 - t1757
            t1760 = 0.5555555558e-1_dp * t719 * t534
            t1762 = 0.5555555558e-1_dp * t1755 * t200
            t1763 = t237 * t722
            t1765 = 0.1000000001e1_dp * t1763 * t1722
            t1766 = t237 * t187
            t1768 = 0.5000000001e0_dp * t1766 * t1722
            t1769 = t199 * t1505
            t1773 = 0.5555555558e-1_dp * t719 * t531
            t1774 = t60 * t496
            t1776 = 0.5000000001e0_dp * t1774 * t709
            t1777 = t379 * t177
            t1779 = 0.5000000001e0_dp * t1777 * t709
            t1782 = t728 * t267
            t1784 = 0.5555555558e-1_dp * t192 * t1782
            t1785 = t579 * t515
            t1787 = 0.1111111112e0_dp * t1785 * t517
            t1789 = 0.5555555558e-1_dp * t508 * t729
            t1790 = t1505 * t172
            t1794 = 0.1000000001e1_dp * t723 * t1712
            t1795 = t498 * t178
            t1797 = 0.8999999998e1_dp * t1795 * t710
            t1800 = 0.5400000000e2_dp * t1143 * t474 * t1410
            t1801 = t20 * t649
            t1804 = 0.1800000000e2_dp * t473 * t1801 * t237
            t1807 = 0.1800000000e2_dp * t473 * t474 * t1447
            t1810 = 0.1200000000e2_dp * t473 * t479 * t1483
            t1811 = t56 * t649
            t1814 = 0.5999999999e1_dp * t478 * t1811 * t17
            t1817 = 0.1800000000e2_dp * t473 * t1150 * t561
            t1821 = t472 * t262
            t1822 = t172 * t1410
            t1824 = 0.3600000000e2_dp * t1821 * t1822
            t1825 = t25 * t1447
            t1828 = 0.1800000000e2_dp * t478 * t183 * t1825
            t1829 = t225 * t561
            t1832 = 0.6000000000e1_dp * t478 * t489 * t1829
            t1835 = 0.1800000000e2_dp * t478 * t1176 * t561
            t1836 = t1800 - t1804 - t1807 - t1810 + t1814 - t1817 + 0.8999999998e1_dp  &
            * t178 * t62 * t1688 + t1824 - t1828 + t1832 - t1835
            t1837 = t68 * t1836
            t1838 = t1837 * t65
            t1842 = -t1760 - t1762 - t1765 - t1768 - 0.5555555558e-1_dp * t192 * t1769  &
            - t1773 - t1776 - t1779 - 0.5555555558e-1_dp * t503 * t1525 - t1784  &
            - t1787 - t1789 - 0.1111111112e0_dp * t516 * t1790 - t1794 + t1797  &
            - 0.5555555558e-1_dp * t192 * t1838 + t1836 * t65
            t1843 = t1758 + t1842
            t1849 = 0.8333333335e-1_dp * t429 * t166 * t1571
            t1852 = 0.8333333335e-1_dp * t429 * t679 * t290
            t1854 = t165 * t466 * t1483
            t1858 = 0.1666666667e0_dp * t445 * t446 * t1825
            t1861 = 0.1800000000e2_dp * t658 * t575 * t425
            t1864 = 0.1666666667e0_dp * t445 * t1292 * t668
            t1865 = t413 * t561
            t1868 = 0.8333333335e-1_dp * t1246 * t446 * t1865
            t1871 = 2 * t147 * t393 * t672
            t1874 = 0.3333333334e0_dp * t170 * t689 * t25
            t1875 = t114 * t585
            t1878 = 0.4444444445e0_dp * t170 * t1875 * t451
            t1881 = 0.1666666667e0_dp * t165 * t437 * t649
            t1884 = -t1496 + 0.9000000000e1_dp * t1498 + 0.6000000000e1_dp * t1502
            t1885 = t1884 * t149
            t1886 = t1885 * t160
            t1891 = 0.5625000000e1_dp * t1302 * t401 * t1410
            t1892 = t19 * t649
            t1895 = 0.2250000000e1_dp * t400 * t1892 * t237
            t1898 = 0.2250000000e1_dp * t400 * t401 * t1447
            t1901 = 0.7500000000e0_dp * t400 * t406 * t1483
            t1902 = t68 * t649
            t1905 = 0.5000000000e0_dp * t405 * t1902 * t17
            t1908 = 0.2250000000e1_dp * t400 * t1309 * t561
            t1912 = t399 * t245
            t1914 = 0.2250000000e1_dp * t1912 * t1822
            t1917 = 0.1500000000e1_dp * t405 * t156 * t1825
            t1920 = 0.1000000000e1_dp * t405 * t417 * t1829
            t1923 = 0.1500000000e1_dp * t405 * t1333 * t561
            t1924 = -t1891 + t1895 + t1898 + t1901 - t1905 + t1908 - 0.1500000000e1_dp  &
            * t151 * t50 * t1688 - t1914 + t1917 - t1920 + t1923
            t1925 = t149 * t1924
            t1928 = t56 * t1531
            t1932 = t27 * t1531
            t1938 = 0.1666666667e0_dp * t165 * t462 * t1447
            t1939 = 0.3333333334e0_dp * t55 * t57 * t1843 - t1849 - t1852 - 0.5555555555e-1_dp  &
            * t1854 + t1858 + t1861 + t1864 - t1868 + t1871 + t1874  &
            - t1878 + t1881 + (2 * t147 * t1886) + (2 * t147 * t1925)  &
            + 0.1666666667e0_dp * t165 * t1928 * t140 + 0.3333333334e0_dp * t170 * t1932  &
            * t172 + t1938
            t1940 = t27 * t735
            t1943 = 0.1111111111e0_dp * t170 * t1940 * t17
            t1946 = 0.1666666667e0_dp * t165 * t1392 * t561
            t1952 = 0.5555555557e-1_dp * t165 * t689 * t433
            t1954 = t1373 * t446 * t1829
            t1958 = 0.8333333335e-1_dp * t429 * t166 * t1553
            t1961 = 0.1666666667e0_dp * t165 * t679 * t379
            t1962 = t585 * my_rho
            t1965 = 0.1666666667e0_dp * t445 * t1962 * t413
            t1966 = t56 * t735
            t1969 = 0.1666666667e0_dp * t165 * t1966 * t237
            t1971 = t657 * t392 * t59
            t1972 = t177 * t20
            t1973 = t561 * t149
            t1974 = t1973 * t160
            t1975 = t1972 * t1974
            t1977 = 0.1800000000e2_dp * t1971 * t1975
            t1980 = 0.2777777778e-1_dp * t1246 * t1247 * t1483
            t1981 = t1367 * t57
            t1984 = 0.1250000000e0_dp * t1981 * t1247 * t1410
            t1987 = 0.8333333335e-1_dp * t429 * t437 * t593
            t1988 = t649 * t17
            t1991 = 0.5555555557e-1_dp * t165 * t171 * t1988
            t1994 = 0.8333333335e-1_dp * t429 * t166 * t1557
            t1997 = 0.1666666667e0_dp * t165 * t171 * t668
            t2000 = 0.8333333335e-1_dp * t429 * t462 * t1410
            t2001 = -t1943 + t1946 + 0.1666666667e0_dp * t165 * t166 * t1688 - t1952  &
            - 0.2222222222e0_dp * t1954 - t1958 + t1961 + t1965 + t1969 + t1977  &
            + t1980 + t1984 - t1987 - t1991 - t1994 + t1997 - t2000
            t2002 = t1939 + t2001
            t2006 = -t1541 - t1543 - t1545 + t1549 - t1552 - t1556 - t1560 - t1563  &
            - 0.3333333336e0_dp * t83 * t1565 - t1570 - t1574 - t1577 - 0.3333333336e0_dp  &
            * t83 * t1689 - t1693 - t1696 - t1699 - t1702 - t1705 -  &
            t1708 - 0.2222222224e0_dp * t144 * t46 * t2002
            e_ndrho_rho(ii) = e_ndrho_rho(ii) +  (-t1406 - t1408 - t1413 - t1417 - t1451 - t1455 - 0.4444444448e0_dp  &
            * t43 * t1536 + t8 * t2006) * sx
            t2009 = t566 * t393
            t2013 = t566 * t149 * t47
            t2014 = t398 * t19
            t2018 = t150 * t68
            t2022 = t19 * t1447
            t2027 = t65 * t237
            t2028 = t575 * t2027
            t2031 = t56 * t561
            t2032 = t2031 * t526
            t2035 = t20 * t1447
            t2042 = -0.1800000000e2_dp * t473 * t2028 + 0.5999999999e1_dp * t478 * t2032  &
            + 0.8999999998e1_dp * t478 * t2035 * t65 + t1509 - t1512 + t1517  &
            - t1523 + t1515 - 0.5000000001e0_dp * t60 * t1447 * t65 - t1530
            t2046 = -0.3000000000e1_dp * t2009 * t569 + 0.4500000000e1_dp * t2013 * t2014  &
            * t1410 - 0.1000000000e1_dp * t2013 * t2018 * t1483 - 0.3000000000e1_dp  &
            * t567 * t151 * t2022 - t1476 - 0.5555555557e-1_dp * t1485 + t1492  &
            + t1482 + t1479 - t1489 + 0.3333333334e0_dp * t55 * t57 * t2042
            t2050 = t45 * t2046
            t2054 = t786 * t640
            t2057 = t1636 + t1638 - t1644 - t1647 - t1649 + t1654 - t1656 + t1664  &
            + t1662 - t1660 + t1671 + t1679 - t1675 - 0.240e2_dp * t2054 * t1060
            t2058 = t105 * t2057
            t2061 = t1579 + t1584 - t1581 - t1586 + t1588 - t1591 - t1596 + t1594  &
            + t1601 - t1605 + t1607 + t1615 - t1630 - t1610 - t1620 + t1634  &
            + t1628 + t1613 + t1624 - 0.20e1_dp * t102 * t2058
            t2073 = -t1891 + t1901 + t1908 + t1898 + t1895 - t1905 - 0.1500000000e1_dp  &
            * t151 * t50 * t2061 - t1914 - t1920 + t1923 + t1917
            t2074 = t149 * t2073
            t2077 = t161 * t17
            t2081 = t657 * t472
            t2082 = t161 * t237
            t2086 = -t1849 - t1852 - 0.5555555557e-1_dp * t1854 + t1858 + t1861 +  &
            t1864 - t1868 + t1871 + t1874 + 0.3333333334e0_dp * t170 * t27 * t2042  &
            * t172 - t1878 + (2 * t147 * t2074) + t1881 + t1938 - t1943  &
            + 0.1200000000e2_dp * t658 * t2031 * t2077 + t1946 - 0.3600000000e2_dp *  &
            t2081 * t575 * t2082
            t2088 = t56 * t2042
            t2095 = t1447 * t65
            t2098 = t140 * t471
            t2099 = t709 * t237
            t2104 = t62 * t2095
            t2107 = -t1714 + t1717 + t1719 - t1724 + t1728 - 0.5000000001e0_dp * t716  &
            * t2095 + 0.1000000000e1_dp * t2098 * t2099 - 0.1000000001e1_dp * t723  &
            * t2095 - t1730 + t1735 + t1743 + 0.8999999998e1_dp * t708 * t2104  &
            + t1752 - t1754 - t1757 - t1760
            t2113 = t177 * pi
            t2114 = t707 * t2113
            t2117 = t471 * pi
            t2118 = t707 * t2117
            t2124 = t1800 - t1810 - t1817 - t1807 - t1804 + t1814 + 0.8999999998e1_dp  &
            * t178 * t62 * t2061 + t1824 + t1832 - t1835 - t1828
            t2126 = t68 * t2124 * t65
            t2130 = -0.5000000001e0_dp * t732 * t2095 - 0.5555555558e-1_dp * t58 * t2061  &
            * t71 - t1762 + t1765 + t1768 + 0.5999999999e1_dp * t2114 * t2032  &
            - t1773 - t1776 - 0.1800000000e2_dp * t2118 * t2028 - 0.5555555558e-1_dp  &
            * t192 * t2126 - t1779 - t1784 - t1787 - t1789 - t1794 + t1797 +  &
            t2124 * t65
            t2131 = t2107 + t2130
            t2138 = -t1952 - 0.2222222223e0_dp * t1954 - t1958 + t1961 + t1965 + t1969  &
            + t1977 + t1980 + 0.1666666667e0_dp * t165 * t2088 * t140 + 0.1800000000e2_dp  &
            * t658 * t2035 * t161 + 0.3333333334e0_dp * t55 * t57 * t2131  &
            + t1984 + 0.1666666667e0_dp * t165 * t166 * t2061 - t1987 - t1991  &
            - t1994 + t1997 - t2000
            t2139 = t2086 + t2138
            t2143 = -t1541 + t1549 - t1563 - t1574 - t1560 - t1543 - t1552 - 0.3333333336e0_dp  &
            * t83 * t2050 * t140 - t1577 - t1545 - t1556 - t1570  &
            - 0.3333333336e0_dp * t83 * t84 * t2061 - t1693 - t1699 - t1708 - t1705  &
            - t1696 - t1702 - 0.2222222224e0_dp * t144 * t46 * t2139
            e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) +  (-t1406 - t1413 - t1455 - t1451 - t1408 - t1417 - 0.4444444448e0_dp  &
            * t43 * t46 * t2046 + t8 * t2143) * sx
            t2147 = lambda * t48 * omega
            t2160 = t754 * t237
            t2169 = t237 * t804
            t2182 = 0.1e1_dp / t787 / lambda
            t2183 = t2182 * t944
            t2187 = t949 * t233
            t2190 = t953 * t772
            t2193 = t316 * t800
            t2196 = t772 * t233
            t2197 = t974 * t2196
            t2200 = t233 * t800
            t2201 = t327 * t2200
            t2204 = t10 * t106
            t2205 = 0.1e1_dp / t2204
            t2220 = beta * t785 * t6
            t2226 = -0.6222222223e2_dp * t24 * t2205 * t1002 * t31 - 0.2115555556e3_dp  &
            * t7 * t2183 * t132 + 0.1315555556e3_dp * t786 / t787 / t85 * t1019  &
            * t371 - 0.4266666668e2_dp * t2220 / t1014 / t94 * t1058
            t2227 = t105 * t2226
            t2230 = -0.9125925923e2_dp * t7 * t2183 * t36 - 0.5866666667e2_dp * t962  &
            * t2187 - 0.3200000001e2_dp * t321 * t2190 + 0.1600000000e2_dp * t321 *  &
            t2193 - 0.120e2_dp * t102 * t2197 + 0.120e2_dp * t102 * t2201 - 0.20e1_dp *  &
            t102 * t2227
            t2257 = 0.5400000000e2_dp * t1142 * t62 * t2160
            t2261 = 0.3600000000e2_dp * t473 * t56 * t754 * t17
            t2264 = 0.5400000000e2_dp * t473 * t1493 * t804
            t2268 = 0.6000000000e1_dp * t478 * t27 * t237 * t14
            t2272 = 0.1800000000e2_dp * t478 * t56 * t804 * t17
            t2274 = t178 * t62 * t2230
            t2276 = pi * t114
            t2279 = 0.2666666667e1_dp * t61 * t2276 * t86
            t2280 = t2257 - t2261 - t2264 - t2268 + t2272 + 0.8999999998e1_dp * t2274  &
            - t2279
            t2281 = t2280 * t65
            t2282 = t881 * t267
            t2285 = t865 * t266
            t2287 = t68 * t2285 * t65
            t2292 = t863 * t266
            t2295 = t2285 * t65
            t2296 = t58 * t2230
            t2299 = t194 * t865
            t2300 = t2299 * t526
            t2304 = t68 * t2280 * t65
            t2307 = t22 * t266
            t2308 = t2307 * t1200
            t2311 = t194 * t863
            t2312 = t2311 * t526
            t2319 = t1099 * t86
            t2326 = -0.1851851853e0_dp * t508 * t874 + t2281 - 0.1666666668e0_dp * t192  &
            * t2282 - 0.5555555558e-1_dp * t192 * t2287 - 0.1666666668e0_dp * t508  &
            * t886 + (3 * t2292 * t65) + t2295 - 0.5555555558e-1_dp * t2296  &
            * t71 + 0.1111111112e0_dp * t192 * t2300 - 0.5555555558e-1_dp * t192 *  &
            t2304 - 0.1851851853e0_dp * t192 * t2308 + 0.1111111112e0_dp * t192 * t2312  &
            + 0.2222222223e0_dp * t508 * t878 - 0.1666666668e0_dp * t508 * t882  &
            + 0.1646090535e0_dp * t192 * t2319 - 0.1666666668e0_dp * t1129 * t275 +  &
            0.1111111112e0_dp * t1129 * t271
            t2351 = t44 * t22
            t2368 = t237 * t14
            t2378 = t804 * t17
            t2382 = t754 * t17
            t2392 = (2 * t147 * t1284 * t249) + (2 * t147 * t1388 * t249)  &
            - 0.3456790122e0_dp * t170 * t1268 * t86 + 0.4444444444e0_dp * t170 *  &
            t1089 * t14 + 0.1666666667e0_dp * t165 * t166 * t2230 + 0.5000000001e0_dp  &
            * t165 * t437 * t804 + 0.3333333334e0_dp * t55 * t57 * t2326 + 0.2e1_dp  &
            * (t147) * t149 * (-0.5625000000e1_dp * t1301 * t50 * t2160 +  &
            0.2250000000e1_dp * t400 * t68 * t754 * t17 + 0.6750000000e1_dp * t400 *  &
            t1459 * t804 + 0.1000000000e1_dp * t405 * t194 * t237 * t14 - 0.1500000000e1_dp  &
            * t405 * t68 * t804 * t17 - 0.1500000000e1_dp * t151 * t50 *  &
            t2230 + 0.1111111111e1_dp * t49 * t2351 * t86) - 0.3333333333e0_dp * t170  &
            * t1254 * t17 + (4 * t147 * t393 * t826) + 0.5000000001e0_dp *  &
            t165 * t1275 * t237 + 0.2222222222e0_dp * t165 * t450 * t2368 - 0.3333333333e0_dp  &
            * t165 * t455 * t833 - 0.2500000000e0_dp * t429 * t166 * t2169  &
            - 0.1666666667e0_dp * t165 * t171 * t2378 + 0.8333333333e-1_dp * t429  &
            * t171 * t2382 - 0.2500000000e0_dp * t429 * t437 * t754 + 0.1250000000e0_dp  &
            * t1368 * t166 * t2160
            t2404 = t45 * t1397
            t2418 = t379 * t754
            t2424 = t1067 * t237
            t2437 = t7 * t792
            t2439 = 0.1e1_dp / t20 / t1436
            t2451 = t940 * t944
            t2455 = t97 * t973
            t2468 = t944 * t104
            t2480 = 0.1e1_dp / t19 / t942
            t2492 = t1014 * my_rho
            t2537 = t764 * t1002
            t2552 = t360 * t1058
            t2565 = t126 / t131 / t1056 / t130
            t2569 = -0.8088888890e3_dp * t24 * t2182 * t2480 * t31 * my_rho - 0.1518222222e4_dp  &
            * t2437 * t2439 * t132 * my_rho + 0.6542222223e3_dp * t786 / t2204  &
            / t19 / t2492 * t218 * t371 * my_rho - 0.6222222223e2_dp * t24 * t2205  &
            * t1002 * t126 * t132 + 0.5600000000e2_dp * t1024 * t341 * t360  &
            * t132 + 0.3288888889e3_dp * t366 / t790 * t126 * t371 - 0.2400000000e2_dp  &
            * t350 * t114 * t1044 * t132 - 0.1520000000e3_dp * t366 * t1034 *  &
            t360 * t371 - 0.3626666667e3_dp * t1054 * t1034 * t1060 + 0.60e1_dp * t119  &
            * t27 * (0.3640e4_dp / 0.27e2_dp * t1040 * t1002 - 0.3640e4_dp / 0.27e2_dp  &
            * my_ndrho * t2182 * t2480 * my_rho) * t132 + 0.8088888890e3_dp * t24 * t2537  &
            * t31 + 0.1518222222e4_dp * t7 * t2451 * t132 - 0.6542222223e3_dp *  &
            t786 * t1012 * t1019 * t371 + 0.2400000000e2_dp * t366 * t368 * t1044  &
            * t371 + 0.9600000000e2_dp * t1055 * t2552 * t220 - 0.1173333333e3_dp  &
            * t1055 * t1059 * t765 + 0.2133333333e3_dp * beta * t785 * t365 * t368  &
            * t2565 * t796
            t2576 = 0.40e1_dp * t102 * t327 * t136 * t2226 + 0.3519999999e3_dp * t93  &
            * t2190 - 0.1551407408e4_dp * t2437 * t2439 * t36 * my_rho + 0.960e2_dp *  &
            t1597 * t1598 * t2200 + 0.120e2_dp * t102 * t327 * t375 * t800 + 0.1551407408e4_dp  &
            * t7 * t2451 * t36 - 0.960e2_dp * t93 * t2455 * my_rho * t2196  &
            + 0.120e2_dp * t102 * t327 * t1063 * t233 + 0.1600000000e2_dp * t321  &
            * t316 * t1063 + 0.9125925923e2_dp * t7 * t2182 * t2468 * t136 + 0.1760000000e3_dp  &
            * t307 * t949 * t958 - 0.360e2_dp * t1616 * t1617 * t2200  &
            - 0.20e1_dp * t102 * t105 * t2569 - 0.6400000001e2_dp * t321 * t953 * t979
            t2585 = 0.1e1_dp / t972 / t35
            t2611 = t311 * t326
            t2621 = -0.3200000001e2_dp * t321 * t953 * t983 - 0.160e2_dp * t93 * t316  &
            * my_rho * t2226 + 0.480e2_dp * t102 * t22 * t2585 * t136 * t2196 - 0.1760000000e3_dp  &
            * t93 * t2193 + 0.960e2_dp * t302 * t2197 + 0.160e2_dp * t302  &
            * t2227 - 0.960e2_dp * t302 * t2201 - 0.8213333332e3_dp * t941 * t2468  &
            * t317 + 0.8213333332e3_dp * t307 * t2187 + 0.9600000002e2_dp * t321 *  &
            t2455 * t975 - 0.360e2_dp * t102 * t974 * t375 * t772 - 0.3519999999e3_dp  &
            * t307 * t2611 * t954 + 0.1173333334e3_dp * t962 * t2611 * t328 -  &
            0.5866666667e2_dp * t962 * t949 * t375
            t2622 = t2576 + t2621
            t2631 = t282 * t140
            t2635 = -0.1000000001e1_dp * t55 * t923 - 0.3333333336e0_dp * t83 * t384  &
            * t2230 - 0.5000000004e0_dp * t165 * t897 * t804 - 0.6666666672e0_dp * t383  &
            * t2404 - 0.1000000001e1_dp * t83 * t1071 * t804 - 0.5000000004e0_dp  &
            * t164 * t1076 - 0.1000000001e1_dp * t289 * t294 * t914 - 0.5000000004e0_dp  &
            * t289 * t294 * t918 + 0.2500000002e0_dp * t906 * t84 * t2418 -  &
            0.2000000001e1_dp * t55 * t1079 - 0.5000000004e0_dp * t289 * t84 * t2424  &
            + 0.2500000002e0_dp * t1546 * t897 * t2169 - 0.3333333336e0_dp * t83 *  &
            t84 * t2622 + 0.2500000002e0_dp * t906 * t294 * t907 - 0.1000000001e1_dp  &
            * t55 * t1082 - 0.1000000001e1_dp * t165 * t2631 * t237
            t2636 = t78 * t379
            t2645 = t379 * t804
            t2655 = t140 * t2230
            t2707 = t74 * t379
            t2720 = t56 * t1239
            t2724 = 0.5000000001e0_dp * t165 * t1392 * t804 - 0.1333333334e1_dp * t165  &
            * t450 * t1172 - 0.2666666667e1_dp * t170 * t1089 * t225 + 0.8333333333e-1_dp  &
            * t429 * t466 * t2382 + 0.1000000000e1_dp * t170 * t1254 * t25  &
            + 0.1666666667e0_dp * t165 * t166 * t2622 + 0.1666666667e0_dp * t165 *  &
            t56 * t2326 * t140 - 0.2500000000e0_dp * t429 * t166 * t2424 - 0.3333333333e0_dp  &
            * t165 * t455 * t1085 - 0.1728395062e0_dp * t165 * t1268 *  &
            t140 * t86 - 0.2500000000e0_dp * t1246 * t446 * t413 * t804 + 0.1666666666e0_dp  &
            * t1246 * t2707 * t833 + 0.1250000000e0_dp * t1368 * t462 * t2160  &
            - 0.1666666667e0_dp * t165 * t466 * t2378 + (6 * t147 * t1283  &
            * t394) + 0.5000000001e0_dp * t165 * t2720 * t237
            t2732 = my_rho * t2205
            t2759 = t25 * t2230
            t2763 = t1367 * t444
            t2779 = t25 * t2160
            t2783 = -0.2500000000e0_dp * t429 * t1392 * t754 + (6 * t147 * t393  &
            * t1347) - 0.3456790123e1_dp * t170 * t1002 * t74 * t2732 + 0.5000000001e0_dp  &
            * t165 * t1275 * t379 - 0.1333333334e1_dp * t170 * t114 * t889  &
            * t451 + 0.3333333334e0_dp * t170 * t27 * t2326 * t172 - 0.2500000000e0_dp  &
            * t429 * t171 * t1168 + 0.3750000000e0_dp * t1981 * t1247 * t2169  &
            + 0.5000000001e0_dp * t165 * t171 * t1093 + 0.8333333334e-1_dp * t1246  &
            * t1247 * t2378 + 0.1666666667e0_dp * t445 * t446 * t2759 - 0.1250000000e0_dp  &
            * t2763 * t1247 * t2382 + 0.5000000001e0_dp * t165 * t437 * t1067  &
            + 0.2222222222e0_dp * t165 * t1089 * t1360 + 0.1000000000e1_dp * t165  &
            * t455 * t413 - 0.2500000000e0_dp * t429 * t462 * t2169 + 0.1250000000e0_dp  &
            * t2763 * t446 * t2779
            t2808 = t19 * t1067
            t2816 = t225 * t754
            t2821 = 0.1e1_dp / t41 / t1140
            t2823 = t47 * t2821 * t44
            t2833 = -0.1687500000e2_dp * t1302 * t401 * t2169 + 0.6750000000e1_dp * t400  &
            * t1309 * t804 + 0.4500000000e1_dp * t400 * t1320 * t833 - 0.1500000000e1_dp  &
            * t400 * t1316 * t2368 + 0.1000000000e1_dp * t405 * t194 * t379  &
            * t14 + 0.8888888888e1_dp * t155 * t97 * my_rho * t2205 + 0.6750000000e1_dp  &
            * t400 * t2808 * t237 - 0.1500000000e1_dp * t405 * t68 * t1067  &
            * t17 + 0.4500000000e1_dp * t400 * t417 * t2816 + 0.1968750000e2_dp * t2823  &
            * t401 * t2160 + 0.1500000000e1_dp * t405 * t156 * t2759 + 0.2250000000e1_dp  &
            * t400 * t401 * t2230
            t2843 = t777 * t237
            t2853 = t172 * t2169
            t2869 = t225 * t804
            t2873 = t194 * t225
            t2877 = -0.1500000000e1_dp * t151 * t50 * t2622 + 0.5625000000e1_dp * t1302  &
            * t156 * t2779 - 0.6750000000e1_dp * t400 * t1333 * t754 + 0.4999999999e1_dp  &
            * t405 * t1340 * t2843 - 0.1000000000e2_dp * t49 * t2351 * t777  &
            - 0.5625000000e1_dp * t1302 * t406 * t2382 - 0.6750000000e1_dp * t1912  &
            * t2853 - 0.5555555555e0_dp * t405 * t22 * t140 * t86 + 0.4500000000e1_dp  &
            * t405 * t1333 * t804 + 0.2250000000e1_dp * t400 * t406 * t2378 -  &
            0.1687500000e2_dp * t1302 * t1309 * t754 - 0.3000000000e1_dp * t405 * t417  &
            * t2869 - 0.6000000000e1_dp * t405 * t2873 * t237
            t2885 = t428 * t1372
            t2935 = -0.2500000000e0_dp * t429 * t166 * t2645 + (2 * t147 * t149  &
            * (t2833 + t2877)) + 0.5000000001e0_dp * t445 * t1292 * t1093 + 0.3333333333e0_dp  &
            * t2885 * t446 * t2816 + 0.2222222222e0_dp * t165 * t1377  &
            * t2368 - 0.3333333333e0_dp * t165 * t1264 * t833 + 0.2e1_dp * (t147)  &
            * (t2257 - t2261 - t2264 - t2268 + t2272 + 0.9000000000e1_dp * t2274  &
            - t2279) * (t149) * t160 - 0.3333333333e0_dp * t170 * t27 * t1239  &
            * t17 - 0.8333333335e-1_dp * t429 * t166 * t2655 + (6 * t147  &
            * t1284 * t424) + 0.2e1_dp * (t147) * t1387 * t392 * (t149) *  &
            t160 - 0.2500000000e0_dp * t1246 * t1292 * t1168 - 0.1666666667e0_dp *  &
            t165 * t1254 * t433 + (6 * t147 * t1388 * t424) + 0.1666666667e0_dp  &
            * t165 * t462 * t2230 + 0.2222222222e0_dp * t165 * t450 * t379 * t14  &
            - 0.5000000001e0_dp * t429 * t437 * t914
            t2937 = t164 * t45 * t341
            t2992 = 0.1e1_dp / t1140 / t40
            t2994 = t59 * t2992 * pi
            t2998 = t27 * t225
            t3026 = 0.1620000000e3_dp * t1143 * t474 * t2169 + 0.1080000000e3_dp * t1143  &
            * t479 * t2382 - 0.2160000000e3_dp * t2994 * t474 * t2160 + 0.3600000000e2_dp  &
            * t478 * t2998 * t237 + 0.1620000000e3_dp * t1143 * t1150 *  &
            t754 + 0.2666666667e1_dp * t478 * t114 * t140 * t86 - 0.3600000000e2_dp  &
            * t473 * t479 * t2378 + 0.1080000000e3_dp * t1821 * t2853 - 0.7200000000e2_dp  &
            * t473 * t1161 * t833 - 0.1800000000e2_dp * t473 * t474 * t2230  &
            + 0.8999999998e1_dp * t178 * t62 * t2622 - 0.5400000000e2_dp * t478 *  &
            t1176 * t804
            t3030 = t20 * t1067
            t3069 = -0.1800000000e2_dp * t478 * t183 * t2759 - 0.5400000000e2_dp * t473  &
            * t3030 * t237 - 0.6000000000e1_dp * t478 * t27 * t379 * t14 - 0.5400000000e2_dp  &
            * t473 * t1150 * t804 + 0.1800000000e2_dp * t478 * t489 *  &
            t2869 + 0.1080000000e3_dp * t473 * t1176 * t754 + 0.1800000000e2_dp * t478  &
            * t56 * t1067 * t17 - 0.3600000000e2_dp * t473 * t489 * t2816 - 0.1866666667e2_dp  &
            * t182 * t357 * t2205 + 0.2400000000e2_dp * t61 * t2276  &
            * t777 - 0.2400000000e2_dp * t478 * t123 * t2843 - 0.1080000000e3_dp *  &
            t1143 * t183 * t2779 + 0.1200000000e2_dp * t473 * t1157 * t2368
            t3070 = t3026 + t3069
            t3085 = t1210 * t67
            t3093 = (3 * t496 * t865 * t65) + (3 * t1189 * t266 * t65)  &
            + 0.1111111112e1_dp * (t192) * (t2307) * (t65) * (t225)  &
            + (t187 * t2280 * t65) + (3 * t496 * t863 * t65) + 0.2222222223e0_dp  &
            * t516 * t1111 * t526 + (t187 * t2285 * t65) + (t3070  &
            * t65) - 0.1851851853e0_dp * t1193 * t498 * t1200 + 0.2222222223e0_dp  &
            * t508 * t1119 + 0.1111111112e0_dp * (t192) * t194 * (t1189)  &
            * t526 + 0.2222222223e0_dp * t1137 * t1221 - 0.1666666668e0_dp * t1129  &
            * t534 + 0.1111111112e0_dp * t3085 * t271 - 0.5555555558e-1_dp * (t192)  &
            * (t199) * (t2281) - 0.1851851853e0_dp * t508 * t1201
            t3097 = t66 * t69
            t3098 = t1109 * t267
            t3129 = t2296 * t67
            t3141 = -0.3333333336e0_dp * t192 * t2299 * t1125 - 0.1666666668e0_dp * t3097  &
            * t3098 + 0.1111111112e0_dp * t503 * t2300 + 0.1111111112e0_dp * t1129  &
            * t527 - 0.1666666668e0_dp * t1129 * t531 - 0.5555555558e-1_dp * t192  &
            * t68 * t3070 * t65 + 0.1646090535e0_dp * t192 * t97 * t187 * t65 *  &
            t86 + 0.5555555559e0_dp * t1129 * t512 - 0.1666666668e0_dp * t192 * t530  &
            * t864 - 0.1666666668e0_dp * t1106 * t882 - 0.1666666668e0_dp * t192 *  &
            t530 * t866 - 0.1481481482e1_dp * t192 * t1099 * t777 - 0.1111111112e0_dp  &
            * t3129 * t196 + 0.5555555559e0_dp * t1193 * t866 * t451 - 0.1666666668e0_dp  &
            * t503 * t2282 - 0.1666666668e0_dp * t3085 * t275 + 0.1111111112e0_dp  &
            * t503 * t2312
            t3181 = -0.1851851853e0_dp * t1106 * t874 - 0.1666666668e0_dp * t508 * t1235  &
            - 0.3333333336e0_dp * t508 * t1224 - 0.1111111112e0_dp * t516 * t2281  &
            * t172 + 0.1111111112e0_dp * t516 * t1097 * t526 + 0.2222222223e0_dp *  &
            t1106 * t878 - 0.1851851853e0_dp * t503 * t2308 - 0.3333333336e0_dp * t1137  &
            * t1216 + (3 * t3098) + 0.1111111112e0_dp * t516 * t1109 * t526  &
            - 0.5555555558e-1_dp * t503 * t2287 - 0.3333333336e0_dp * t192 * t2311  &
            * t1125 - 0.6666666672e0_dp * t508 * t1126 + 0.1646090535e0_dp * t503  &
            * t2319 - 0.5555555558e-1_dp * t3129 * t200 - 0.5555555558e-1_dp * t58  &
            * t2622 * t71 - 0.1666666668e0_dp * t1106 * t886
            t3232 = -0.1666666668e0_dp * t192 * t1206 * t267 - 0.5555555558e-1_dp * t192  &
            * t199 * t2295 - 0.1666666668e0_dp * t508 * t1122 - 0.3333333336e0_dp  &
            * t1129 * t520 - 0.1111111112e0_dp * t516 * t2295 * t172 - 0.1481481482e1_dp  &
            * t66 * t67 * t97 * t267 * t1100 - 0.1851851853e0_dp * t192 *  &
            t22 * t496 * t1200 - 0.3333333336e0_dp * t516 * t2292 * t65 * my_rho *  &
            t25 + 0.1810699590e1_dp * t192 * t311 * t70 * t2732 - 0.5555555558e-1_dp  &
            * t503 * t2304 - 0.1481481482e1_dp * t508 * t1101 - 0.3333333336e0_dp *  &
            t867 * t515 * t517 + 0.1111111112e1_dp * t268 * t1192 * t1194 + 0.5555555559e0_dp  &
            * t1193 * t864 * t451 + 0.1111111112e1_dp * t508 * t1213 -  &
            0.1666666668e0_dp * t508 * t1207 - 0.3333333336e0_dp * t1137 * t1134
            t3244 = t278 * t140
            t3248 = t889 * my_rho
            t3262 = omega * t1300
            t3264 = t140 * t2160
            t3268 = 0.1555555556e1_dp * t2937 * t446 * t2843 + 0.3111111111e1_dp * t170  &
            * t1268 * t777 - 0.2500000000e0_dp * t429 * t1275 * t290 - 0.1333333334e1_dp  &
            * t1373 * t1292 * t1172 + 0.3111111111e1_dp * t170 * t341 * t278  &
            * t1100 - 0.3456790122e0_dp * t170 * t341 * t203 * t86 + 0.3750000000e0_dp  &
            * t1368 * t166 * t2418 + 0.4444444444e0_dp * t170 * t114 * t537  &
            * t14 + 0.3333333334e0_dp * t55 * t57 * (t3093 + t3141 + t3181 + t3232)  &
            - 0.1111111111e0_dp * t2885 * t1247 * t2368 - 0.6666666668e0_dp * t1373  &
            * t446 * t2869 + 0.1666666667e0_dp * t1246 * t3244 * t833 + 0.5000000001e0_dp  &
            * t445 * t3248 * t413 - 0.1666666667e0_dp * t165 * t171 *  &
            t1067 * t17 - 0.2500000000e0_dp * t429 * t437 * t918 + 0.3750000000e0_dp  &
            * t1368 * t437 * t907 - 0.3125000000e0_dp * t3262 * t45 * t166 * t3264
            t3278 = my_rho * t398
            t3288 = -0.1000000001e1_dp * t165 * t2636 * t237 - 0.2000000001e1_dp * t55  &
            * t926 - 0.1000000001e1_dp * t83 * t922 * t379 - 0.5000000004e0_dp * t289  &
            * t84 * t2645 + 0.2500000002e0_dp * t429 * t897 * t754 - 0.1000000001e1_dp  &
            * t83 * t294 * t1067 - 0.1666666668e0_dp * t289 * t84 * t2655  &
            - 0.1000000001e1_dp * t83 * t2404 * t237 - 0.5000000004e0_dp * t289 * t922  &
            * t290 - 0.5000000004e0_dp * t289 * t384 * t2169 - 0.5000000004e0_dp  &
            * t289 * t1071 * t754 - 0.2222222224e0_dp * t144 * t46 * (t2724 + t2783  &
            + t2935 + t3268) - 0.3333333336e0_dp * t83 * t45 * t2392 * t140 -  &
            0.1250000001e0_dp * t3278 * omega * t84 * t3264 - 0.1000000001e1_dp * t55  &
            * t1068 + 0.8333333340e-1_dp * t906 * t384 * t2160
            t3291 = -0.1000000001e1_dp * t2147 * t84 * t754 - 0.4000000002e1_dp * t747  &
            * t294 * t237 - 0.2000000001e1_dp * t747 * t84 * t804 - 0.1333333334e1_dp  &
            * t214 * t894 + 0.1666666668e0_dp * t4 * t150 * t46 * t78 * t2160  &
            - 0.1000000001e1_dp * t753 * t46 * t282 * t754 - 0.1000000001e1_dp * t1409  &
            * t84 * t2169 - 0.2000000001e1_dp * t217 * t46 * t893 * t237 - 0.2000000001e1_dp  &
            * t217 * t46 * t282 * t804 - 0.6666666672e0_dp * t217 *  &
            t46 * t78 * t2230 - 0.4444444448e0_dp * t43 * t46 * t2392 + t8 * (t2635  &
            + t3288)
          END IF
          IF( grad_deriv >= 3 .OR. grad_deriv == -3 ) THEN
            e_rho_rho_rho(ii) = e_rho_rho_rho(ii) +  t3291 * sx
            t3305 = t4 * t428
            t3306 = t754 * t561
            t3314 = t237 * t1447
            t3329 = t804 * t561
            t3339 = t949 * t557
            t3344 = t953 * t1425
            t3347 = t316 * t1443
            t3352 = t772 * t557
            t3353 = t974 * t3352
            t3356 = t233 * t1443
            t3357 = t327 * t3356
            t3362 = t800 * t557
            t3363 = t327 * t3362
            t3366 = beta * t777
            t3378 = 0.1866666667e2_dp * t3366 * t342 + 0.9866666667e2_dp * t24 * t782  &
            - 0.8266666668e2_dp * t366 * t797 + 0.3200000001e2_dp * t1054 / t1014 /  &
            t12 * t1058
            t3379 = t105 * t3378
            t3382 = 0.3911111110e2_dp * t24 * t766 - 0.1955555555e2_dp * t962 * t3339  &
            + 0.2133333334e2_dp * t339 * t769 - 0.2133333334e2_dp * t321 * t3344 +  &
            0.1066666667e2_dp * t321 * t3347 + 0.80e1_dp * t112 * t773 - 0.120e2_dp * t102  &
            * t3353 + 0.80e1_dp * t102 * t3357 - 0.40e1_dp * t112 * t801 + 0.40e1_dp  &
            * t102 * t3363 - 0.20e1_dp * t102 * t3379
            t3449 = t1447 * t17
            t3453 = (2 * t147 * t1885 * t249) + 0.1800000000e2_dp * (t1971)  &
            * (t1972) * (t1973) * (t249) + (2 * t147 * t393  &
            * t1470) + 0.1800000000e2_dp * t658 * t575 * t827 + 0.2e1_dp * (t147)  &
            * t149 * (-0.5625000000e1_dp * t1302 * t19 * t754 * t561 + 0.4500000000e1_dp  &
            * t400 * t1459 * t1447 + 0.1500000000e1_dp * t400 * t815 * t1483  &
            - 0.1000000000e1_dp * t405 * t68 * t1447 * t17 + 0.2250000000e1_dp *  &
            t400 * t19 * t804 * t561 - 0.1500000000e1_dp * t151 * t50 * t3382 + 0.3333333334e0_dp  &
            * t405 * t194 * t14 * t561) + 0.1250000000e0_dp * t1368  &
            * t166 * t3306 - 0.8333333335e-1_dp * t429 * t679 * t754 - 0.1666666667e0_dp  &
            * t429 * t166 * t3314 + 0.5555555555e-1_dp * t1246 * t74 * t237  &
            * t1483 - 0.1111111111e0_dp * t165 * t689 * t833 - 0.1111111111e0_dp *  &
            t165 * t171 * t3449
            t3472 = t14 * t561
            t3491 = 0.5400000000e2_dp * t1143 * t20 * t754 * t561
            t3494 = 0.3600000000e2_dp * t473 * t1493 * t1447
            t3497 = 0.2400000000e2_dp * t473 * t852 * t1483
            t3501 = 0.1200000000e2_dp * t478 * t56 * t1447 * t17
            t3505 = 0.1800000000e2_dp * t473 * t20 * t804 * t561
            t3507 = t178 * t62 * t3382
            t3512 = 0.2000000000e1_dp * t478 * t27 * t14 * t561
            t3513 = t3491 - t3494 - t3497 + t3501 - t3505 + 0.8999999998e1_dp * t3507  &
            - t3512
            t3514 = t3513 * t65
            t3525 = t58 * t3382
            t3533 = t237 * t1518
            t3543 = t3514 + 0.8999999998e1_dp * t863 * t59 * t177 * t710 + (2  &
            * t267 * t1504) + 0.8999999998e1_dp * t865 * t59 * t177 * t710 - 0.5555555558e-1_dp  &
            * t3525 * t71 - 0.5000000001e0_dp * t804 * t177 * t709 +  &
            0.7407407410e-1_dp * t1755 * t271 + 0.6666666668e0_dp * t3533 * t1722 -  &
            0.1111111112e0_dp * t1755 * t275 - 0.1111111112e0_dp * t508 * t1525 - 0.1000000001e1_dp  &
            * t237 * t266 * t1722
            t3546 = 0.1e1_dp / t12
            t3547 = t60 * t3546
            t3552 = t194 * t1504
            t3553 = t3552 * t526
            t3561 = t68 * t3513 * t65
            t3569 = t274 * t1505
            t3575 = -0.6172839508e-1_dp * t719 * t874 - 0.5555555556e0_dp * t3547 * t709  &
            + 0.7407407410e-1_dp * t719 * t878 + 0.7407407410e-1_dp * t192 * t3553  &
            + 0.6666666668e0_dp * t1519 * t1712 - 0.5555555558e-1_dp * t719 * t882  &
            - 0.5555555558e-1_dp * t192 * t3561 - 0.5000000001e0_dp * t60 * t863 *  &
            t709 - 0.5555555558e-1_dp * t719 * t886 - 0.1111111112e0_dp * t192 * t3569  &
            - 0.5000000001e0_dp * t60 * t865 * t709
            t3576 = t3543 + t3575
            t3580 = -0.1666666667e0_dp * t429 * t437 * t1410 + 0.3333333334e0_dp * t165  &
            * t1928 * t237 + 0.3333333334e0_dp * t165 * t437 * t1447 - 0.8333333335e-1_dp  &
            * t429 * t166 * t3329 + 0.1666666667e0_dp * t165 * t679 * t804  &
            + 0.1666666667e0_dp * t165 * t166 * t3382 + 0.7407407405e-1_dp * t165  &
            * t450 * t3472 + 0.1481481481e0_dp * t170 * t1875 * t14 - 0.1111111111e0_dp  &
            * t165 * t455 * t1483 - 0.2222222222e0_dp * t170 * t1932 * t17 +  &
            0.1666666667e0_dp * t165 * t1275 * t561 + 0.3333333334e0_dp * t55 * t57  &
            * t3576
            t3581 = t3453 + t3580
            t3608 = t973 * my_rho
            t3618 = 0.640e2_dp * t1597 * t1598 * t3356 - 0.640e2_dp * t1582 * t773 -  &
            0.1173333333e3_dp * t93 * t3347 - 0.3911111110e2_dp * t1000 * t963 - 0.320e2_dp  &
            * t605 * t959 + 0.1173333333e3_dp * t307 * t949 * t1602 + 0.40e1_dp  &
            * t102 * t327 * t136 * t3378 - 0.960e2_dp * t1597 * t3608 * t3352 -  &
            0.40e1_dp * t112 * t1064 + 0.80e1_dp * t112 * t984 + 0.2737777778e3_dp * t307  &
            * t3339
            t3668 = t2552 * t23
            t3673 = t635 * t1058
            t3674 = t3673 * t220
            t3677 = -0.1866666667e3_dp * t1645 * t342 - 0.6346666667e3_dp * t24 * t991  &
            + 0.3946666667e3_dp * t366 * t997 + 0.1866666667e3_dp * beta * t764 *  &
            t1004 + 0.60e1_dp * t550 * t1046 + 0.60e1_dp * t119 * t27 * (-0.280e3_dp /  &
            0.9e1_dp * t987 + 0.280e3_dp / 0.9e1_dp * t2537 * my_rho) * t132 - 0.60e1_dp * t7  &
            * t640 * t1044 * t371 + 0.1600000000e2_dp * t366 * t368 * t1667 * t371  &
            - 0.4800000000e2_dp * t1680 * t3668 + 0.1600000000e3_dp * t1680 * t1060  &
            + 0.3200000000e2_dp * t1055 * t3674
            t3679 = t2565 * t1439
            t3682 = t24 * t940
            t3711 = -0.1600000000e3_dp * t2220 * t368 * t3679 + 0.6346666667e3_dp * t3682  &
            * t1008 - 0.3946666667e3_dp * t366 * t1012 * t1021 + 0.1866666667e2_dp  &
            * t3366 * t1026 + 0.1866666667e2_dp * t1024 * t341 * t635 * t132  &
            - 0.1706666667e3_dp * t7 * t1036 - 0.1600000000e2_dp * t1429 * t1030 - 0.1600000000e2_dp  &
            * t350 * t114 * t1667 * t132 + 0.6400000000e2_dp * t7 *  &
            t1050 - 0.5066666667e2_dp * t366 * t1034 * t635 * t371 + 0.1520000000e3_dp  &
            * t786 * t1034 * t1681
            t3726 = -0.240e2_dp * t1616 * t1617 * t3356 + 0.160e2_dp * t302 * t3379 +  &
            0.40e1_dp * t102 * t327 * t1063 * t557 + 0.960e2_dp * t302 * t3353 - 0.320e2_dp  &
            * t302 * t3363 - 0.2737777778e3_dp * t941 * t2468 * t608 - 0.120e2_dp  &
            * t102 * t974 * t645 * t772 - 0.2133333334e2_dp * t321 * t953 * t1621  &
            - 0.20e1_dp * t102 * t105 * (t3677 + t3711) + 0.80e1_dp * t102 * t327  &
            * t375 * t1443 - 0.120e2_dp * t1616 * t1617 * t3362 - 0.2133333334e2_dp  &
            * t321 * t953 * t1625
            t3750 = t2585 * t136
            t3763 = 0.2346666666e3_dp * t1589 * t950 + 0.80e1_dp * t102 * t327 * t1684  &
            * t233 + 0.6400000002e2_dp * t7 * t220 * t1618 - 0.240e2_dp * t112 * t976  &
            - 0.160e2_dp * t93 * t316 * my_rho * t3378 + 0.2346666666e3_dp * t93 *  &
            t3344 + 0.3911111110e2_dp * t962 * t2611 * t614 - 0.5475555556e3_dp * t24  &
            * t930 + 0.480e2_dp * t1616 * t3750 * t3352 - 0.1955555555e2_dp * t962  &
            * t949 * t645 - 0.240e2_dp * t1616 * t973 * t375 * t1425 + 0.2133333334e2_dp  &
            * t339 * t969
            t3794 = -0.640e2_dp * t302 * t3357 - 0.4266666668e2_dp * t339 * t966 + 0.640e2_dp  &
            * t605 * t955 + 0.40e1_dp * t102 * t327 * t645 * t800 + 0.5475555556e3_dp  &
            * t3682 * t946 + 0.1066666667e2_dp * t321 * t316 * t1684 + 0.160e2_dp  &
            * t112 * t980 - 0.2346666666e3_dp * t7 * t929 * t1599 + 0.320e2_dp  &
            * t1597 * t1598 * t3362 - 0.2346666666e3_dp * t605 * t769 + 0.320e2_dp  &
            * t1582 * t801 - 0.2133333334e2_dp * t321 * t953 * t1631
            t3796 = t3618 + t3726 + t3763 + t3794
            t3826 = 0.8333333340e-1_dp * t906 * t597 * t907 + 0.1666666668e0_dp * t429  &
            * t1547 - 0.3333333336e0_dp * t83 * t84 * t3796 - 0.3333333336e0_dp *  &
            t83 * t922 * t649 - 0.3333333336e0_dp * t289 * t294 * t1571 - 0.3333333336e0_dp  &
            * t289 * t1071 * t1410 - 0.6666666672e0_dp * t83 * t1564 * t379  &
            - 0.6666666672e0_dp * t55 * t1700 - 0.3333333336e0_dp * t165 * t207 *  &
            t237 * t561 - 0.6666666672e0_dp * t55 * t1689 - 0.3333333336e0_dp * t83  &
            * t1694 * t804 + 0.1666666668e0_dp * t1546 * t897 * t3314
            t3827 = t649 * t804
            t3842 = t3278 * t46
            t3849 = t140 * t3382
            t3853 = t649 * t754
            t3861 = t78 * t649
            t3865 = t589 * t140
            t3869 = -0.1666666668e0_dp * t289 * t84 * t3827 - 0.1666666668e0_dp * t289  &
            * t922 * t593 + 0.1666666668e0_dp * t1546 * t2636 * t1410 - 0.6666666672e0_dp  &
            * t55 * t1703 - 0.3333333336e0_dp * t83 * t597 * t1067 - 0.1250000001e0_dp  &
            * t3842 * t897 * t3306 - 0.3333333336e0_dp * t83 * t384 *  &
            t3382 - 0.1666666668e0_dp * t289 * t84 * t3849 + 0.8333333340e-1_dp * t906  &
            * t84 * t3853 - 0.6666666672e0_dp * t55 * t1568 - 0.6666666672e0_dp *  &
            t55 * t1565 - 0.3333333336e0_dp * t165 * t3861 * t237 - 0.3333333336e0_dp  &
            * t165 * t3865 * t237
            t3871 = t379 * t1447
            t3875 = t45 * t2002
            t3892 = t225 * t1447
            t3918 = t25 * t3382
            t3925 = t777 * t561
            t3933 = t451 * t1410
            t3939 = t1301 * t245
            t3940 = t172 * t3306
            t3965 = 0.1500000000e1_dp * t405 * t156 * t3918 + 0.1968750000e2_dp * t2823  &
            * t401 * t3306 + 0.1666666666e1_dp * t405 * t1340 * t3925 + 0.3000000000e1_dp  &
            * t405 * t1333 * t1447 + 0.3000000000e1_dp * t399 * t822 * t3933  &
            + 0.2250000000e1_dp * t400 * t401 * t3382 + 0.5625000000e1_dp * t3939  &
            * t3940 + 0.3333333333e0_dp * t405 * t194 * t649 * t14 - 0.1125000000e2_dp  &
            * t1302 * t1309 * t1410 + 0.2250000000e1_dp * t400 * t1892 * t804  &
            - 0.5625000000e1_dp * t1302 * t401 * t3329 - 0.5625000000e1_dp * t1302  &
            * t1892 * t754 + 0.2250000000e1_dp * t400 * t2808 * t561 - 0.4500000000e1_dp  &
            * t400 * t1333 * t1410
            t3969 = t172 * t3314
            t3975 = t19 * t1688
            t3979 = t290 * t1483
            t3985 = t172 * t3329
            t4010 = -0.1125000000e2_dp * t1302 * t401 * t3314 - 0.4500000000e1_dp * t1912  &
            * t3969 + 0.4500000000e1_dp * t400 * t1309 * t1447 + 0.4500000000e1_dp  &
            * t400 * t3975 * t237 - 0.3750000000e1_dp * t3939 * t3979 - 0.2000000000e1_dp  &
            * t405 * t417 * t3892 - 0.2250000000e1_dp * t1912 * t3985 -  &
            0.2000000000e1_dp * t405 * t2873 * t561 - 0.5000000000e0_dp * t400 * t1316  &
            * t3472 - 0.1000000000e1_dp * t405 * t68 * t1688 * t17 + 0.1500000000e1_dp  &
            * t400 * t406 * t3449 - 0.1500000000e1_dp * t151 * t50 * t3796  &
            + 0.1500000000e1_dp * t400 * t1902 * t833 + 0.1500000000e1_dp * t400 *  &
            t1320 * t1483
            t4018 = 0.2222222222e0_dp * t2885 * t446 * t1172 * t561 - 0.1666666667e0_dp  &
            * t1246 * t1292 * t1865 - 0.4444444445e0_dp * t1373 * t446 * t3892  &
            - 0.1111111111e0_dp * t165 * t171 * t1688 * t17 + 0.1666666667e0_dp * t165  &
            * t56 * t3576 * t140 + 0.1250000000e0_dp * t2763 * t446 * t1168 *  &
            t561 + 0.3333333334e0_dp * t165 * t689 * t413 + 0.1666666667e0_dp * t165  &
            * t679 * t1067 - 0.1666666667e0_dp * t1246 * t446 * t413 * t1447 +  &
            (2 * t147 * t149 * (t3965 + t4010)) - 0.4444444445e0_dp * t1373 *  &
            t1292 * t1829
            t4025 = t3262 * t57
            t4032 = t74 * t649
            t4070 = t3525 * t67
            t4079 = -0.5555555558e-1_dp * t192 * t199 * t3514 + (t705 * t863 *  &
            t65) + (2 * t496 * t1504 * t65) + (t187 * t3513 * t65) -  &
            0.5555555558e-1_dp * t192 * t728 * t864 - 0.5555555556e0_dp * t140 * t3546  &
            * t1722 + (t705 * t865 * t65) - 0.1111111112e0_dp * t1106 * t1525  &
            - 0.1111111112e0_dp * t4070 * t196 + 0.6666666668e0_dp * t1748 * t266  &
            * t1722 + 0.3703703706e0_dp * t1193 * t1505 * t451
            t4101 = t498 * t1505
            t4108 = -0.5000000001e0_dp * t60 * t1189 * t709 + 0.3703703706e0_dp * t1755  &
            * t512 - 0.5555555558e-1_dp * t503 * t3561 - 0.2222222224e0_dp * t516  &
            * t267 * t172 * t1504 - 0.5555555558e-1_dp * t1715 * t886 - 0.5555555558e-1_dp  &
            * t1715 * t882 - 0.1111111112e0_dp * t516 * t3514 * t172 - 0.5000000001e0_dp  &
            * t1067 * t177 * t709 + (2 * t4101) - 0.1111111112e0_dp  &
            * t508 * t1769 - 0.1111111112e0_dp * t192 * t1837 * t267
            t4140 = -0.5555555558e-1_dp * t719 * t1207 - 0.2000000001e1_dp * t1763 *  &
            t266 * t1722 - 0.1111111112e0_dp * t192 * t530 * t1505 + 0.1333333334e1_dp  &
            * t1519 * t722 * t1712 + 0.7407407410e-1_dp * t516 * t1718 * t526  &
            + 0.6666666668e0_dp * t379 * t1518 * t1722 - 0.5555555558e-1_dp * t192 *  &
            t728 * t866 - 0.2222222224e0_dp * t719 * t1126 + 0.1800000000e2_dp * t1111  &
            * t178 * t710 - 0.1111111112e0_dp * t508 * t1782 - 0.5555555558e-1_dp  &
            * t58 * t3796 * t71
            t4156 = t865 * t561 * t65
            t4201 = t20 * t1688
            t4214 = -0.8000000000e1_dp * t478 * t123 * t3925 - 0.2400000000e2_dp * t473  &
            * t479 * t3449 + 0.1200000000e2_dp * t478 * t489 * t3892 + 0.1200000000e2_dp  &
            * t478 * t2998 * t561 - 0.2000000000e1_dp * t478 * t27 * t649  &
            * t14 + 0.1200000000e2_dp * t478 * t56 * t1688 * t17 - 0.2400000000e2_dp  &
            * t473 * t1161 * t1483 + 0.1080000000e3_dp * t1143 * t474 * t3314 -  &
            0.2160000000e3_dp * t2994 * t474 * t3306 + 0.7200000000e2_dp * t1821 *  &
            t3969 - 0.3600000000e2_dp * t473 * t4201 * t237 - 0.1800000000e2_dp * t473  &
            * t3030 * t561 - 0.2400000000e2_dp * t473 * t1811 * t833 - 0.1800000000e2_dp  &
            * t478 * t183 * t3918
            t4230 = t1142 * t262
            t4255 = 0.5400000000e2_dp * t1143 * t474 * t3329 + 0.5400000000e2_dp * t1143  &
            * t1801 * t754 + 0.7200000000e2_dp * t473 * t1176 * t1410 + 0.8999999998e1_dp  &
            * t178 * t62 * t3796 - 0.2400000000e2_dp * t472 * t859 * t3933  &
            - 0.1080000000e3_dp * t4230 * t3940 + 0.4000000000e1_dp * t473 * t1157  &
            * t3472 + 0.1080000000e3_dp * t1143 * t1150 * t1410 + 0.3600000000e2_dp  &
            * t1821 * t3985 + 0.7200000000e2_dp * t4230 * t3979 - 0.1800000000e2_dp  &
            * t473 * t1801 * t804 - 0.3600000000e2_dp * t473 * t1150 * t1447 -  &
            0.3600000000e2_dp * t478 * t1176 * t1447 - 0.1800000000e2_dp * t473 * t474  &
            * t3382
            t4256 = t4214 + t4255
            t4261 = -0.2222222224e0_dp * t1510 * t515 * t517 + 0.3703703706e0_dp * t719  &
            * t1213 - 0.1000000001e1_dp * t1766 * t266 * t1722 - 0.5555555558e-1_dp  &
            * t4070 * t200 + 0.6666666668e0_dp * t1519 * t496 * t561 * t65 - 0.5000000001e0_dp  &
            * t732 * t4156 - 0.5000000001e0_dp * t140 * t865 * t1722  &
            - 0.5000000001e0_dp * t804 * t187 * t1722 - 0.4938271608e0_dp * t719 *  &
            t1101 + 0.8999999998e1_dp * t1097 * t178 * t710 - 0.5555555558e-1_dp *  &
            t192 * t68 * t4256 * t65
            t4272 = t863 * t561 * t65
            t4292 = (2 * t1836 * t266 * t65) - 0.1111111112e0_dp * t719 * t1224  &
            - 0.2222222224e0_dp * t1137 * t1790 - 0.1000000001e1_dp * t723 * t4272  &
            - 0.1111111112e0_dp * t1785 * t1134 - 0.6172839508e-1_dp * t192 * t22 *  &
            t705 * t1200 - 0.5000000001e0_dp * t140 * t863 * t1722 - 0.1111111112e0_dp  &
            * t1755 * t531 - 0.6172839508e-1_dp * t1715 * t874 - 0.5555555558e-1_dp  &
            * t719 * t1122 + 0.7407407410e-1_dp * t508 * t1726
            t4317 = 0.7407407410e-1_dp * t1755 * t527 - 0.6172839508e-1_dp * t719 * t1201  &
            - 0.1111111112e0_dp * t508 * t1838 - 0.2222222224e0_dp * t192 * t3552  &
            * t1125 + 0.7407407410e-1_dp * t1785 * t1221 - 0.1111111111e1_dp * t3547  &
            * t1737 + 0.6666666668e0_dp * t1519 * t187 * t1712 + 0.7407407410e-1_dp  &
            * t503 * t3553 - 0.1111111112e0_dp * t1755 * t534 - 0.1000000001e1_dp  &
            * t1774 * t1712 + 0.7407407410e-1_dp * t719 * t1119
            t4345 = -0.1000000001e1_dp * t379 * t266 * t1722 - 0.1000000001e1_dp * t804  &
            * t722 * t1722 + 0.8999999998e1_dp * t1109 * t178 * t710 - 0.5555555558e-1_dp  &
            * t1129 * t729 + t4256 * t65 - 0.2222222224e0_dp * t1755 * t520  &
            - 0.1111111112e0_dp * t503 * t3569 - 0.1111111112e0_dp * t1785 * t1216  &
            - 0.1111111112e0_dp * t3097 * t4101 - 0.5000000001e0_dp * t732 * t4272  &
            + 0.7407407410e-1_dp * t192 * t194 * t1836 * t526
            t4352 = t1731 * t67
            t4378 = 0.7407407410e-1_dp * t1715 * t878 + 0.8999999998e1_dp * t1189 * t59  &
            * t177 * t710 + 0.7407407410e-1_dp * t4352 * t271 - 0.1111111112e0_dp  &
            * t4352 * t275 - 0.5555555556e0_dp * t3547 * t1745 + 0.7407407410e-1_dp  &
            * t516 * t1709 * t526 + 0.6666666668e0_dp * t3533 * t187 * t1722 - 0.1000000001e1_dp  &
            * t237 * t496 * t1722 + 0.3703703706e0_dp * t579 * t1192  &
            * t1194 - 0.1000000001e1_dp * t723 * t4156 + 0.1333333334e1_dp * t3533  &
            * t722 * t1722 - 0.5555555558e-1_dp * t719 * t1235
            t4392 = 0.3333333334e0_dp * t165 * t171 * t1825 + 0.1250000000e0_dp * t1981  &
            * t1247 * t3329 - 0.3125000000e0_dp * t4025 * t1247 * t3306 + 0.3333333334e0_dp  &
            * t165 * t1392 * t1447 + 0.5555555556e-1_dp * t1246 * t4032  &
            * t833 - 0.1666666667e0_dp * t1246 * t74 * t25 * t1410 + 0.5555555556e-1_dp  &
            * t1246 * t1247 * t3449 + 0.1800000000e2_dp * t657 * t1283 * t59  &
            * t1975 - 0.1666666667e0_dp * t429 * t166 * t3871 + 0.3333333334e0_dp  &
            * t55 * t57 * (t4079 + t4108 + t4140 + t4261 + t4292 + t4317 + t4345  &
            + t4378) + 0.1481481481e0_dp * t170 * t114 * t735 * t14 - 0.1666666667e0_dp  &
            * t429 * t1928 * t290
            t4428 = -0.1111111111e0_dp * t165 * t1264 * t1483 + 0.1800000000e2_dp * t658  &
            * t575 * t1348 - 0.8333333335e-1_dp * t429 * t679 * t918 - 0.1111111111e0_dp  &
            * t165 * t455 * t1988 + 0.3333333334e0_dp * t170 * t27 * t3576  &
            * t172 - 0.8333333335e-1_dp * t429 * t1966 * t754 + 0.1250000000e0_dp  &
            * t1368 * t679 * t907 + 0.1666666667e0_dp * t445 * t446 * t3918 - 0.4444444445e0_dp  &
            * t165 * t450 * t1829 + 0.2500000000e0_dp * t1981 * t3244  &
            * t1410 - 0.1666666667e0_dp * t429 * t1392 * t1410
            t4466 = -0.1666666667e0_dp * t429 * t437 * t1553 - 0.1111111111e0_dp * t165  &
            * t689 * t1085 + 0.1666666667e0_dp * t165 * t1275 * t649 - 0.2222222222e0_dp  &
            * t170 * t27 * t1843 * t17 - 0.8333333335e-1_dp * t1246 * t1962  &
            * t1168 + 0.1666666667e0_dp * t165 * t2720 * t561 + 0.6666666668e0_dp  &
            * t170 * t1932 * t25 + 0.3333333334e0_dp * t445 * t1292 * t1825 + 0.1250000000e0_dp  &
            * t1368 * t166 * t3853 - 0.4444444445e0_dp * t1373 * t1962  &
            * t1172 + 0.7407407405e-1_dp * t165 * t1377 * t3472 + 0.5555555555e-1_dp  &
            * t1246 * t3244 * t1483
            t4498 = t56 * t1843
            t4505 = 0.1666666667e0_dp * t165 * t166 * t3796 - 0.8333333334e-1_dp * t2763  &
            * t1247 * t833 * t561 + (4 * t147 * t392 * t1886) + (2  &
            * t147 * t1388 * t672) - 0.1111111111e0_dp * t165 * t1940 * t833 + 0.3333333334e0_dp  &
            * t165 * t1928 * t379 + 0.7407407409e-1_dp * t165 * t450  &
            * t649 * t14 - 0.1666666667e0_dp * t429 * t462 * t3314 - 0.3703703704e-1_dp  &
            * t2885 * t1247 * t3472 + 0.3333333334e0_dp * t165 * t4498 * t237  &
            + (2 * t147 * t1284 * t672)
            t4512 = t1531 * my_rho
            t4545 = (4 * t147 * t1885 * t424) + 0.1666666667e0_dp * t165 * t462  &
            * t3382 + 0.3333333334e0_dp * t445 * t4512 * t413 + 0.5185185185e0_dp  &
            * t2937 * t446 * t3925 - 0.8333333335e-1_dp * t429 * t166 * t3827 - 0.1666666667e0_dp  &
            * t429 * t437 * t1571 + 0.5555555555e-1_dp * t1246 * t203  &
            * t237 * t1483 + 0.5555555555e-1_dp * t1246 * t2707 * t1483 + 0.7407407409e-1_dp  &
            * t165 * t1875 * t1360 + (4 * t147 * t393 * t1924)  &
            + 0.3600000000e2_dp * (t1971) * (t1972) * (t1973) * (t424)  &
            + 0.3333333334e0_dp * t165 * t437 * t1688
            t4553 = t585 * t140
            t4575 = t1067 * t561
            t4585 = 0.1666666667e0_dp * t445 * t3248 * t668 - 0.1111111111e0_dp * t165  &
            * t1932 * t433 + 0.5555555556e-1_dp * t1246 * t4553 * t833 - 0.8333333335e-1_dp  &
            * t429 * t1275 * t593 - 0.8333333335e-1_dp * t429 * t462 *  &
            t3329 + 0.2500000000e0_dp * t1981 * t1247 * t3314 + 0.1666666667e0_dp *  &
            t445 * t1962 * t1093 + 0.3333333334e0_dp * t165 * t455 * t668 - 0.8333333335e-1_dp  &
            * t429 * t166 * t3849 - 0.8333333335e-1_dp * t429 * t166 *  &
            t4575 + 0.2e1_dp * t147 * (t3491 - t3494 - t3497 + t3501 - t3505 + 0.9000000000e1_dp  &
            * t3507 - t3512) * t149 * t160
            t4607 = t1688 * t237
            t4627 = 0.1666666667e0_dp * t165 * t1966 * t804 - 0.1666666667e0_dp * t429  &
            * t437 * t1557 + 0.1800000000e2_dp * t657 * t1387 * t59 * t1975 + 0.1037037037e1_dp  &
            * t170 * t341 * t585 * t1100 - 0.8888888890e0_dp * t170  &
            * t114 * t1531 * t451 - 0.8888888890e0_dp * t170 * t1875 * t225 - 0.1666666667e0_dp  &
            * t429 * t166 * t4607 - 0.1111111111e0_dp * t165 * t466  &
            * t3449 - 0.1666666667e0_dp * t429 * t679 * t914 - 0.8333333335e-1_dp *  &
            t1246 * t446 * t1093 * t561 + 0.1250000000e0_dp * t1368 * t462 * t3306  &
            + 0.2500000000e0_dp * t1981 * t2707 * t1410
            t4657 = -0.3333333336e0_dp * t289 * t84 * t3871 - 0.6666666672e0_dp * t83  &
            * t3875 * t237 + 0.8333333340e-1_dp * t906 * t384 * t3306 - 0.1666666668e0_dp  &
            * t289 * t597 * t918 - 0.2222222224e0_dp * t144 * t46 * (t4018  &
            + t4392 + t4428 + t4466 + t4505 + t4545 + t4585 + t4627) - 0.1666666668e0_dp  &
            * t289 * t384 * t3329 - 0.3333333336e0_dp * t165 * t2636 * t561  &
            - 0.3333333336e0_dp * t165 * t897 * t1447 + 0.8333333340e-1_dp * t1546  &
            * t897 * t3329 - 0.6666666672e0_dp * t83 * t1071 * t1447 - 0.3333333336e0_dp  &
            * t83 * t2404 * t561 - 0.4444444448e0_dp * t383 * t3875 - 0.3333333336e0_dp  &
            * t289 * t294 * t1553
            t4696 = -0.6666666672e0_dp * t55 * t1575 - 0.3333333336e0_dp * t289 * t597  &
            * t914 - 0.3333333336e0_dp * t165 * t2631 * t561 - 0.6666666672e0_dp *  &
            t55 * t1706 - 0.3333333336e0_dp * t289 * t384 * t3314 - 0.6666666672e0_dp  &
            * t83 * t294 * t1688 + 0.1666666668e0_dp * t1546 * t2631 * t1410 -  &
            0.3333333336e0_dp * t289 * t1564 * t290 - 0.1666666668e0_dp * t289 * t84  &
            * t4575 - 0.3333333336e0_dp * t289 * t294 * t1557 - 0.1666666668e0_dp  &
            * t289 * t1694 * t754 - 0.3333333336e0_dp * t289 * t84 * t4607 - 0.3333333336e0_dp  &
            * t83 * t45 * t3581 * t140
            t4700 = -0.6666666670e0_dp * t2147 * t1411 - 0.1333333334e1_dp * t747 * t597  &
            * t237 - 0.1333333334e1_dp * t747 * t84 * t1447 - 0.1333333334e1_dp  &
            * t747 * t294 * t561 - 0.8888888896e0_dp * t214 * t1536 + 0.1666666668e0_dp  &
            * t3305 * t84 * t3306 - 0.3333333336e0_dp * t753 * t46 * t589 * t754  &
            - 0.6666666672e0_dp * t1409 * t84 * t3314 - 0.6666666670e0_dp * t1409  &
            * t294 * t1410 - 0.1333333334e1_dp * t217 * t46 * t1535 * t237 - 0.1333333334e1_dp  &
            * t217 * t46 * t282 * t1447 - 0.3333333336e0_dp * t1409  &
            * t84 * t3329 - 0.6666666672e0_dp * t217 * t46 * t589 * t804 - 0.6666666672e0_dp  &
            * t217 * t46 * t78 * t3382 - 0.6666666672e0_dp * t217 * t46  &
            * t893 * t561 - 0.4444444448e0_dp * t43 * t46 * t3581 + t8 * (t3826 +  &
            t3869 + t4657 + t4696)
            e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii) +  t4700 * sx
            t4701 = t561 ** 2
            t4704 = 0.3333333336e0_dp * t2147 * t84 * t4701
            t4707 = 0.1333333334e1_dp * t747 * t597 * t561
            t4708 = t22 * t36
            t4713 = t557 ** 2
            t4714 = t327 * t4713
            t4717 = t22 * t132
            t4721 = 0.1e1_dp / t787 * t2480
            t4725 = 0.120e2_dp * t624 * t4717 - 0.60e1_dp * t7 * t4721 * t371
            t4726 = t105 * t4725
            t4729 = 0.40e1_dp * t624 * t4708 - 0.80e1_dp * t112 * t558 + 0.40e1_dp * t102  &
            * t4714 - 0.20e1_dp * t102 * t4726
            t4732 = 0.6666666672e0_dp * t747 * t84 * t4729
            t4733 = pi ** 2
            t4734 = t4733 * t146
            t4736 = 0.1e1_dp / t58 / omega
            t4738 = t4734 * t4736 * t1300
            t4739 = t4701 * t149
            t4743 = t19 * t4701
            t4744 = t399 * t4743
            t4748 = t151 * t19 * t4729
            t4760 = t20 * t4701
            t4764 = t20 * t4729
            t4768 = t58 ** 2
            t4769 = 0.1e1_dp / t4768
            t4770 = t4769 * t1141
            t4771 = t4770 * t4733
            t4772 = t26 * t4701
            t4776 = t58 * t4729
            t4782 = -0.2250000000e2_dp * t473 * t4760 * t65 + 0.8999999998e1_dp * t478  &
            * t4764 * t65 + 0.8099999996e2_dp * t4771 * t4772 * t65 - 0.5555555558e-1_dp  &
            * t4776 * t71 - 0.5000000001e0_dp * t60 * t4729 * t65
            t4786 = -0.2700000000e2_dp * t4738 * t18 * t4739 + 0.4500000000e1_dp * t567  &
            * t4744 - 0.3000000000e1_dp * t567 * t4748 - 0.8333333335e-1_dp * t429  &
            * t166 * t4701 + 0.3333333334e0_dp * t165 * t679 * t561 + 0.1666666667e0_dp  &
            * t165 * t166 * t4729 + 0.3333333334e0_dp * t55 * t57 * t4782
            t4789 = 0.4444444448e0_dp * t214 * t46 * t4786
            t4790 = t237 * t4701
            t4793 = 0.1666666668e0_dp * t3305 * t84 * t4790
            t4796 = 0.6666666672e0_dp * t1409 * t597 * t1410
            t4797 = t1447 * t561
            t4800 = 0.6666666672e0_dp * t1409 * t84 * t4797
            t4801 = t237 * t4729
            t4804 = 0.3333333336e0_dp * t1409 * t84 * t4801
            t4808 = 0.6666666672e0_dp * t217 * t46 * t4786 * t237
            t4812 = 0.1333333334e1_dp * t217 * t46 * t589 * t1447
            t4817 = t953 * t4713
            t4820 = t316 * t4725
            t4829 = t233 * t4713
            t4830 = t974 * t4829
            t4833 = t1443 * t557
            t4834 = t327 * t4833
            t4837 = t233 * t4725
            t4838 = t327 * t4837
            t4849 = -0.3200000000e2_dp * t1645 * t627 + 0.4800000000e2_dp * t7 * t1440  &
            - 0.2400000000e2_dp * t786 / t2492 * t1058
            t4850 = t105 * t4849
            t4853 = -0.1066666667e2_dp * t1645 * t98 + 0.2133333334e2_dp * t339 * t1420  &
            - 0.1066666667e2_dp * t321 * t4817 + 0.5333333333e1_dp * t321 * t4820  &
            - 0.40e1_dp * t624 * t234 + 0.160e2_dp * t112 * t1426 - 0.80e1_dp * t112 *  &
            t1444 - 0.120e2_dp * t102 * t4830 + 0.80e1_dp * t102 * t4834 + 0.40e1_dp *  &
            t102 * t4838 - 0.20e1_dp * t102 * t4850
            t4857 = 0.6666666672e0_dp * t217 * t46 * t78 * t4853
            t4861 = 0.3333333336e0_dp * t753 * t46 * t282 * t4701
            t4864 = t217 * t46 * t1535 * t561
            t4869 = 0.6666666672e0_dp * t217 * t46 * t282 * t4729
            t4872 = 0.8333333335e-1_dp * t429 * t437 * t4701
            t4874 = t165 * t1928 * t561
            t4878 = 0.1666666667e0_dp * t165 * t437 * t4729
            t4881 = 0.5400000000e2_dp * t1143 * t1493 * t4701
            t4884 = 0.3600000000e2_dp * t473 * t2035 * t561
            t4887 = 0.1800000000e2_dp * t473 * t1493 * t4729
            t4889 = t178 * t62 * t4853
            t4893 = 0.1200000000e2_dp * t473 * t1500 * t4701
            t4895 = t478 * t1500 * t4729
            t4897 = t4881 - t4884 - t4887 + 0.8999999998e1_dp * t4889 - t4893 + 0.5999999999e1_dp  &
            * t4895
            t4898 = t4897 * t65
            t4901 = t1504 * t59 * t177 * t710
            t4904 = t4701 * t65
            t4905 = t62 * t4904
            t4907 = 0.2250000000e2_dp * t1506 * t471 * t4905
            t4908 = t4729 * t65
            t4909 = t62 * t4908
            t4911 = 0.8999999998e1_dp * t1507 * t4909
            t4915 = t4733 * t26 * t4904
            t4917 = 0.8099999996e2_dp * t266 * t4769 * t1141 * t4915
            t4918 = t58 * t4853
            t4920 = 0.5555555558e-1_dp * t4918 * t71
            t4922 = t1447 * t177 * t709
            t4925 = 0.5000000001e0_dp * t1513 * t4908
            t4929 = t60 * t1504 * t709
            t4932 = 0.5000000001e0_dp * t1528 * t4908
            t4935 = 0.1000000000e1_dp * t237 * t471 * t4904
            t4940 = t59 * pi * t20 * t65
            t4941 = t237 * t1141 * t4701 * t4940
            t4943 = t4776 * t67
            t4945 = 0.3703703705e-1_dp * t4943 * t271
            t4946 = t471 * t1518
            t4951 = 0.5555555558e-1_dp * t4943 * t275
            t4952 = t719 * t1525
            t4955 = t68 * t4897 * t65
            t4958 = t4898 + 0.1800000000e2_dp * t4901 - t4907 + t4911 + t4917 - t4920  &
            - 0.1000000000e1_dp * t4922 - t4925 + 0.3333333334e0_dp * t1519 * t4908  &
            - 0.1000000000e1_dp * t4929 - t4932 + t4935 - 0.4500000000e1_dp * t4941  &
            + t4945 + 0.3000000000e1_dp * t4946 * t4701 * t4940 - t4951 - 0.1111111112e0_dp  &
            * t4952 - 0.5555555558e-1_dp * t192 * t4955
            t4994 = t44 * t4733 * t146
            t4995 = t4994 * t4770
            t5001 = 0.1666666667e0_dp * t165 * t166 * t4853
            t5002 = t17 * t4701
            t5005 = 0.2777777778e-1_dp * t429 * t171 * t5002
            t5008 = 0.1111111111e0_dp * t165 * t689 * t1483
            t5009 = t17 * t4729
            t5011 = t165 * t171 * t5009
            t5013 = t27 * t4782
            t5016 = 0.1111111111e0_dp * t170 * t5013 * t17
            t5019 = 0.1250000000e0_dp * t1368 * t166 * t4790
            t5022 = 0.1666666667e0_dp * t429 * t679 * t1410
            t5025 = 0.1666666667e0_dp * t429 * t166 * t4797
            t5028 = 0.8333333335e-1_dp * t429 * t166 * t4801
            t5029 = t56 * t4782
            t5032 = 0.1666666667e0_dp * t165 * t5029 * t237
            t5035 = 0.3333333334e0_dp * t165 * t679 * t1447
            t5036 = -t4872 + 0.3333333334e0_dp * t4874 + t4878 + 0.3333333334e0_dp *  &
            t55 * t57 * t4958 + 0.2e1_dp * t147 * t149 * (-0.5625000000e1_dp * t1302  &
            * t1459 * t4701 + 0.4500000000e1_dp * t400 * t2022 * t561 + 0.2250000000e1_dp  &
            * t400 * t1459 * t4729 - 0.1500000000e1_dp * t151 * t50 * t4853  &
            + 0.7500000000e0_dp * t400 * t1466 * t4701 - 0.5000000000e0_dp * t405  &
            * t1466 * t4729) - 0.3600000000e2_dp * t2081 * t4760 * t250 + 0.1800000000e2_dp  &
            * t658 * t4764 * t250 + 0.3600000000e2_dp * t658 * t575 * t1471  &
            + 0.1620000000e3_dp * t4995 * t4772 * t250 + t5001 + t5005 - t5008  &
            - 0.5555555555e-1_dp * t5011 - t5016 + t5019 - t5022 - t5025 - t5028  &
            + t5032 + t5035
            t5040 = t379 * t4729
            t5043 = 0.1666666668e0_dp * t289 * t84 * t5040
            t5044 = t140 * t4853
            t5047 = 0.1666666668e0_dp * t289 * t84 * t5044
            t5050 = 0.1666666668e0_dp * t165 * t897 * t4729
            t5051 = beta * t86
            t5060 = beta * t92
            t5065 = my_rho * t4713
            t5069 = my_rho * t4725
            t5079 = t136 * t4713
            t5083 = t645 * t557
            t5087 = t136 * t4725
            t5106 = beta * t640
            t5108 = t126 * t371 * my_ndrho
            t5115 = t366 * t640
            t5118 = -0.480e2_dp * t5051 * t4717 + 0.240e2_dp * t7 * t2182 * t2480 * t371  &
            + 0.480e2_dp * t5060 * t628 - 0.240e2_dp * t7 * t1640 * t1438 * t371  &
            * my_rho + 0.120e2_dp * t550 * t637 - 0.180e2_dp * t5106 * t5108 - 0.120e2_dp  &
            * t7 * t640 * t635 * t371 + 0.180e2_dp * t5115 * t1681
            t5122 = -0.320e2_dp * t5051 * t4708 + 0.640e2_dp * t1582 * t558 - 0.320e2_dp  &
            * t302 * t4714 + 0.160e2_dp * t302 * t4726 + 0.320e2_dp * t5060 * t99 -  &
            0.640e2_dp * t605 * t609 + 0.320e2_dp * t93 * t953 * t5065 - 0.160e2_dp *  &
            t93 * t316 * t5069 - 0.40e1_dp * t624 * t137 + 0.160e2_dp * t112 * t615  &
            - 0.80e1_dp * t112 * t646 - 0.120e2_dp * t102 * t974 * t5079 + 0.80e1_dp *  &
            t102 * t327 * t5083 + 0.40e1_dp * t102 * t327 * t5087 - 0.20e1_dp * t102  &
            * t105 * t5118
            t5123 = t5122 * t237
            t5126 = 0.1666666668e0_dp * t289 * t84 * t5123
            t5139 = t25 * t4701
            t5143 = t25 * t4729
            t5147 = -0.5625000000e1_dp * t1302 * t401 * t4701 + 0.4500000000e1_dp * t400  &
            * t1892 * t561 + 0.2250000000e1_dp * t400 * t401 * t4729 - 0.1500000000e1_dp  &
            * t151 * t50 * t5122 - 0.2250000000e1_dp * t400 * t156 * t5139  &
            + 0.1500000000e1_dp * t405 * t156 * t5143
            t5163 = t140 * t4701
            t5170 = t649 * t561
            t5174 = t140 * t4729
            t5226 = 0.5400000000e2_dp * t1143 * t474 * t4701 - 0.3600000000e2_dp * t473  &
            * t1801 * t561 - 0.1800000000e2_dp * t473 * t474 * t4729 + 0.8999999998e1_dp  &
            * t178 * t62 * t5122 + 0.3600000000e2_dp * t473 * t183 * t5139  &
            - 0.1800000000e2_dp * t478 * t183 * t5143
            t5228 = t705 * t59
            t5229 = t5228 * t177
            t5237 = t649 * t177
            t5244 = t60 * t705
            t5247 = t187 * t4769
            t5251 = t58 * t5122
            t5256 = t140 * t1141
            t5262 = t471 * t722
            t5270 = t68 * t5226
            t5271 = t5270 * t65
            t5276 = t5226 * t65 + 0.1800000000e2_dp * t5229 * t710 - 0.2250000000e2_dp  &
            * t707 * t471 * t4905 + 0.8999999998e1_dp * t708 * t4909 - 0.1000000000e1_dp  &
            * t5237 * t709 - 0.5000000001e0_dp * t716 * t4908 - 0.1000000001e1_dp  &
            * t723 * t4908 - 0.1000000000e1_dp * t5244 * t709 + 0.8099999996e2_dp  &
            * t5247 * t1141 * t4915 - 0.5555555558e-1_dp * t5251 * t71 + 0.1000000000e1_dp  &
            * t2098 * t4904 - 0.4500000000e1_dp * t5256 * t4701 * t4940  &
            - 0.1111111112e0_dp * t4943 * t196 - 0.9000000007e1_dp * t5262 * t4701  &
            * t4940 - 0.5555555558e-1_dp * t4943 * t200 - 0.1111111112e0_dp * t719 *  &
            t729 - 0.5555555558e-1_dp * t192 * t5271 - 0.5000000001e0_dp * t732 * t4908
            t5280 = (2 * t147 * t149 * t5147) - 0.3600000000e2_dp * t2081 * t4760  &
            * t161 + 0.1800000000e2_dp * t658 * t4764 * t161 + 0.3600000000e2_dp  &
            * t658 * t575 * t673 + 0.1620000000e3_dp * t4995 * t4772 * t161 + 0.1250000000e0_dp  &
            * t1368 * t166 * t5163 - 0.1666666667e0_dp * t429 * t679  &
            * t593 - 0.1666666667e0_dp * t429 * t166 * t5170 - 0.8333333335e-1_dp *  &
            t429 * t166 * t5174 + 0.1666666667e0_dp * t165 * t5029 * t140 + 0.3333333334e0_dp  &
            * t165 * t679 * t649 + 0.1666666667e0_dp * t165 * t166 * t5122  &
            - 0.8333333335e-1_dp * t1246 * t446 * t5139 + 0.3333333334e0_dp * t445  &
            * t1962 * t668 + 0.1666666667e0_dp * t445 * t446 * t5143 + 0.3333333334e0_dp  &
            * t170 * t5013 * t172 - 0.8333333335e-1_dp * t429 * t462 * t4701  &
            + 0.3333333334e0_dp * t165 * t1966 * t561 + 0.1666666667e0_dp * t165  &
            * t462 * t4729 + 0.3333333334e0_dp * t55 * t57 * t5276
            t5281 = t45 * t5280
            t5284 = 0.3333333336e0_dp * t83 * t5281 * t237
            t5291 = 0.1666666668e0_dp * t289 * t384 * t4801
            t5294 = 0.1666666668e0_dp * t1546 * t3865 * t1410
            t5295 = t45 * t4786
            t5298 = 0.3333333336e0_dp * t83 * t5295 * t379
            t5301 = 0.3333333336e0_dp * t55 * t384 * t4729
            t5303 = t83 * t1564 * t649
            t5305 = -t5043 - t5047 - t5050 - t5126 - t5284 - 0.3333333336e0_dp * t83  &
            * t45 * t5036 * t140 - t5291 + t5294 - t5298 - t5301 - 0.6666666672e0_dp  &
            * t5303
            t5308 = 0.3333333336e0_dp * t289 * t1694 * t1410
            t5311 = 0.1666666668e0_dp * t289 * t1071 * t4701
            t5314 = 0.1666666668e0_dp * t289 * t5295 * t290
            t5317 = 0.3333333336e0_dp * t165 * t3861 * t561
            t5320 = 0.3333333336e0_dp * t83 * t1071 * t4729
            t5323 = 0.1666666668e0_dp * t1546 * t3861 * t1410
            t5326 = 0.1250000001e0_dp * t3842 * t897 * t4790
            t5329 = 0.6666666672e0_dp * t83 * t1694 * t1447
            t5330 = t379 * t4701
            t5333 = 0.8333333340e-1_dp * t906 * t84 * t5330
            t5335 = 0.2222222224e0_dp * t383 * t5281
            t5338 = 0.6666666672e0_dp * t55 * t1694 * t561
            t5339 = t649 * t1447
            t5342 = 0.3333333336e0_dp * t289 * t84 * t5339
            t5343 = -t5308 - t5311 - t5314 - t5317 - t5320 + t5323 - t5326 - t5329  &
            + t5333 - t5335 - t5338 - t5342
            t5345 = t1688 * t561
            t5347 = t289 * t84 * t5345
            t5351 = 0.1666666668e0_dp * t164 * t384 * t4701
            t5354 = 0.3333333336e0_dp * t83 * t384 * t4853
            t5357 = 0.6666666672e0_dp * t55 * t597 * t649
            t5360 = 0.3333333336e0_dp * t289 * t384 * t4797
            t5363 = 0.8333333335e-1_dp * t429 * t166 * t5040
            t5367 = 0.1666666667e0_dp * t445 * t4782 * my_rho * t413
            t5370 = 0.1666666667e0_dp * t429 * t679 * t1557
            t5371 = t225 * t4701
            t5373 = t2885 * t446 * t5371
            t5376 = t165 * t466 * t5009
            t5379 = t165 * t1928 * t649
            t5383 = 0.4166666667e-1_dp * t2763 * t1247 * t5002
            t5386 = 0.8333333335e-1_dp * t429 * t5029 * t290
            t5389 = 0.8333333335e-1_dp * t429 * t166 * t5044
            t5397 = 0.8333333335e-1_dp * t1246 * t446 * t413 * t4729
            t5400 = 0.3333333334e0_dp * t445 * t1962 * t1825
            t5402 = t165 * t679 * t1688
            t5406 = 0.8333333335e-1_dp * t429 * t166 * t5123
            t5409 = 0.8333333335e-1_dp * t429 * t1392 * t4701
            t5412 = 0.8333333335e-1_dp * t429 * t171 * t5139
            t5415 = 0.1666666667e0_dp * t429 * t462 * t4797
            t5418 = 2 * t147 * t393 * t5147
            t5421 = t657 * t1884 * t59 * t1975
            t5423 = -t5363 + t5367 - t5370 + 0.1111111111e0_dp * t5373 - 0.5555555555e-1_dp  &
            * t5376 + 0.3333333334e0_dp * t5379 - t5383 - t5386 - t5389 +  &
            0.3333333334e0_dp * t170 * t27 * t4958 * t172 - t5397 + t5400 + 0.3333333334e0_dp  &
            * t5402 - t5406 - t5409 - t5412 - t5415 + t5418 + 0.3600000000e2_dp  &
            * t5421
            t5426 = 0.8333333335e-1_dp * t1246 * t1292 * t5139
            t5429 = 0.1666666667e0_dp * t445 * t1292 * t5143
            t5430 = t225 * t4729
            t5432 = t1373 * t446 * t5430
            t5436 = 0.1666666667e0_dp * t165 * t5029 * t379
            t5440 = 0.40e1_dp * t102 * t327 * t375 * t4725
            t5443 = 0.960e2_dp * t1597 * t3608 * t4829
            t5446 = 0.480e2_dp * t24 * t23 * t1618
            t5449 = 0.1173333333e3_dp * t307 * t2611 * t5065
            t5451 = 0.640e2_dp * t605 * t1603
            t5454 = 0.2133333334e2_dp * t321 * t953 * t5083
            t5456 = 0.1280e3_dp * t1582 * t1426
            t5459 = 0.5333333333e1_dp * t321 * t316 * t5118
            t5460 = beta * t306
            t5462 = 0.1173333333e3_dp * t5460 * t313
            t5466 = 0.40e1_dp * t102 * t327 * t136 * t4849
            t5468 = 0.2346666666e3_dp * t605 * t1420
            t5469 = t5440 - t5443 - t5446 - t5449 - t5451 - t5454 - t5456 + t5459  &
            - t5462 + t5466 - t5468
            t5472 = 0.3200000001e2_dp * t321 * t2455 * t5079
            t5474 = 0.4266666668e2_dp * t339 * t1608
            t5478 = 0.240e2_dp * t1616 * t973 * t645 * t1425
            t5480 = 0.1760000000e3_dp * t5060 * t627
            t5482 = 0.2160000000e3_dp * t7 * t1642
            t5484 = 0.1760000000e3_dp * t5460 * t347
            t5487 = 0.2160000000e3_dp * t7 * t995 * t1652
            t5489 = 0.1600000000e2_dp * t1429 * t1658
            t5492 = 0.7200000000e2_dp * beta * t368 * t5108
            t5494 = 0.6400000000e2_dp * t7 * t1677
            t5496 = t366 * t368 * t1681
            t5499 = 0.120e2_dp * t550 * t1669
            t5503 = 0.180e2_dp * t5106 * t360 * t371 * my_ndrho
            t5507 = 0.120e2_dp * t7 * t640 * t1667 * t371
            t5509 = 0.180e2_dp * t5115 * t3668
            t5511 = t1680 * t3673 * t23
            t5516 = t5480 - t5482 - t5484 + t5487 - t5489 + t5492 + t5494 - 0.1920000000e3_dp  &
            * t5496 + t5499 - t5503 - t5507 + t5509 - 0.4800000000e2_dp  &
            * t5511 + 0.1200000000e3_dp * t1055 * t2565 * t4721
            t5523 = 0.120e2_dp * t102 * t974 * t375 * t4713
            t5525 = 0.2133333334e2_dp * t339 * t1611
            t5528 = t102 * t327 * t1684 * t557
            t5531 = 0.160e2_dp * t112 * t1626
            t5533 = 0.160e2_dp * t112 * t1632
            t5535 = 0.640e2_dp * t1582 * t1444
            t5537 = 0.160e2_dp * t112 * t1622
            t5540 = 0.320e2_dp * t1597 * t1598 * t4837
            t5541 = t5472 - t5474 - t5478 - 0.20e1_dp * t102 * t105 * t5516 - t5523  &
            + t5525 + 0.80e1_dp * t5528 + t5531 + t5533 + t5535 + t5537 + t5540
            t5545 = 0.480e2_dp * t1616 * t3750 * t4829
            t5548 = 0.120e2_dp * t1616 * t1617 * t4837
            t5551 = 0.1280e3_dp * t24 * t298 * t1599
            t5553 = 0.160e2_dp * t302 * t4850
            t5558 = 0.320e2_dp * t5060 * t97 * t104 * my_rho * t233
            t5560 = 0.2346666666e3_dp * t1589 * t1592
            t5562 = 0.5866666667e2_dp * t93 * t4820
            t5564 = 0.1173333333e3_dp * t5060 * t98
            t5565 = t112 * t1685
            t5568 = 0.1066666667e2_dp * t1645 * t322
            t5572 = 0.40e1_dp * t102 * t327 * t5118 * t233
            t5574 = 0.1173333333e3_dp * t93 * t4817
            t5575 = t5545 - t5548 + t5551 + t5553 - t5558 + t5560 - t5562 + t5564  &
            - 0.80e1_dp * t5565 + t5568 + t5572 + t5574
            t5579 = 0.80e1_dp * t102 * t327 * t645 * t1443
            t5581 = 0.320e2_dp * t302 * t4838
            t5585 = 0.160e2_dp * t93 * t316 * my_rho * t4849
            t5590 = 0.80e1_dp * t624 * t22 * t326 * t136 * t233
            t5593 = 0.1066666667e2_dp * t321 * t953 * t5087
            t5595 = 0.640e2_dp * t302 * t4834
            t5597 = 0.960e2_dp * t302 * t4830
            t5600 = 0.240e2_dp * t1616 * t1617 * t4833
            t5602 = 0.40e1_dp * t624 * t376
            t5604 = 0.320e2_dp * t5051 * t234
            t5607 = 0.5866666667e2_dp * t307 * t949 * t5069
            t5610 = 0.640e2_dp * t1597 * t1598 * t4833
            t5611 = t5579 - t5581 - t5585 + t5590 - t5593 - t5595 + t5597 - t5600  &
            - t5602 + t5604 + t5607 + t5610
            t5613 = t5469 + t5541 + t5575 + t5611
            t5619 = 0.3333333334e0_dp * t170 * t5013 * t25
            t5622 = 0.1666666667e0_dp * t165 * t1392 * t4729
            t5625 = 0.1666666667e0_dp * t429 * t166 * t5339
            t5628 = 0.3600000000e2_dp * t2081 * t4760 * t425
            t5631 = 0.5555555556e-1_dp * t1246 * t4032 * t1483
            t5632 = t25 * t4853
            t5635 = 0.1500000000e1_dp * t405 * t156 * t5632
            t5638 = 0.2250000000e1_dp * t400 * t1333 * t4701
            t5641 = 0.1500000000e1_dp * t405 * t1333 * t4729
            t5643 = t400 * t3975 * t561
            t5647 = 0.5625000000e1_dp * t1302 * t1309 * t4701
            t5650 = 0.5625000000e1_dp * t1302 * t401 * t4801
            t5653 = 0.4500000000e1_dp * t400 * t1892 * t1447
            t5657 = 0.5000000000e0_dp * t405 * t68 * t5122 * t17
            t5660 = 0.1125000000e2_dp * t1302 * t401 * t4797
            t5663 = 0.1500000000e1_dp * t400 * t417 * t5371
            t5666 = 0.1125000000e2_dp * t1302 * t1892 * t1410
            t5667 = t5635 - t5638 + t5641 + 0.4500000000e1_dp * t5643 - t5647 - t5650  &
            + t5653 - t5657 - t5660 + t5663 - t5666
            t5670 = 0.2250000000e1_dp * t400 * t1309 * t4729
            t5671 = t172 * t4801
            t5673 = 0.2250000000e1_dp * t1912 * t5671
            t5676 = 0.2250000000e1_dp * t400 * t401 * t4853
            t5679 = 0.1875000000e1_dp * t1302 * t406 * t5002
            t5682 = 0.1500000000e1_dp * t400 * t1902 * t1483
            t5686 = 0.2250000000e1_dp * t400 * t19 * t5122 * t237
            t5689 = 0.1000000000e1_dp * t405 * t417 * t5430
            t5692 = 0.1968750000e2_dp * t2823 * t401 * t4790
            t5695 = 0.7500000000e0_dp * t400 * t406 * t5009
            t5696 = t172 * t4790
            t5698 = 0.5625000000e1_dp * t3939 * t5696
            t5699 = t172 * t4797
            t5701 = 0.4500000000e1_dp * t1912 * t5699
            t5705 = t5670 - t5673 + t5676 - t5679 + t5682 + t5686 - t5689 + t5692  &
            + t5695 + t5698 - t5701 - 0.1500000000e1_dp * t151 * t50 * t5613
            t5713 = t4739 * t160
            t5716 = 0.1620000000e3_dp * t4994 * t392 * t4769 * t1141 * t26 * t5713
            t5719 = 0.3333333334e0_dp * t165 * t1966 * t1447
            t5722 = 0.1111111111e0_dp * t165 * t689 * t1988
            t5724 = t429 * t166 * t5345
            t5728 = 0.2500000000e0_dp * t1981 * t1247 * t4797
            t5731 = 0.1666666667e0_dp * t445 * t446 * t5632
            t5734 = 0.1666666667e0_dp * t429 * t679 * t1553
            t5738 = 0.4444444445e0_dp * t170 * t114 * t4782 * t451
            t5741 = 0.1666666667e0_dp * t165 * t171 * t5143
            t5742 = -t5426 + t5429 - 0.2222222222e0_dp * t5432 + t5436 + 0.1666666667e0_dp  &
            * t165 * t166 * t5613 + t5619 + t5622 - t5625 - t5628 + t5631  &
            + (2 * t147 * t149 * (t5667 + t5705)) + t5716 + t5719 - t5722  &
            - 0.1666666667e0_dp * t5724 + t5728 + t5731 - t5734 - t5738 + t5741
            t5745 = t165 * t4498 * t561
            t5749 = 0.1666666667e0_dp * t429 * t437 * t5170
            t5752 = 0.8333333335e-1_dp * t429 * t462 * t4801
            t5756 = 0.1111111111e0_dp * t170 * t27 * t5276 * t17
            t5760 = 0.1250000000e0_dp * t2763 * t446 * t413 * t4701
            t5767 = 0.1666666667e0_dp * t165 * t437 * t5122
            t5770 = 0.1620000000e3_dp * t4995 * t4772 * t425
            t5773 = 0.1250000000e0_dp * t1368 * t437 * t5163
            t5776 = 0.2777777778e-1_dp * t1246 * t1247 * t5009
            t5780 = 0.1666666667e0_dp * t165 * t56 * t5276 * t237
            t5783 = 0.5555555557e-1_dp * t165 * t5013 * t433
            t5786 = 0.3333333334e0_dp * t165 * t689 * t668
            t5789 = 0.1250000000e0_dp * t1368 * t166 * t5330
            t5792 = 0.2777777778e-1_dp * t429 * t466 * t5002
            t5795 = 0.1111111111e0_dp * t165 * t1940 * t1483
            t5797 = t1373 * t1962 * t1829
            t5801 = 0.3600000000e2_dp * t1143 * t479 * t5002
            t5804 = 0.2400000000e2_dp * t473 * t1811 * t1483
            t5807 = 0.2160000000e3_dp * t2994 * t474 * t4790
            t5810 = 0.1080000000e3_dp * t1143 * t1801 * t1410
            t5813 = 0.1800000000e2_dp * t478 * t183 * t5632
            t5816 = 0.3600000000e2_dp * t473 * t1176 * t4701
            t5819 = 0.1800000000e2_dp * t478 * t1176 * t4729
            t5822 = 0.1800000000e2_dp * t473 * t474 * t4853
            t5825 = 0.1200000000e2_dp * t473 * t479 * t5009
            t5828 = 0.3600000000e2_dp * t473 * t1801 * t1447
            t5830 = 0.1080000000e3_dp * t4230 * t5696
            t5831 = t5801 - t5804 - t5807 + t5810 - t5813 + t5816 - t5819 - t5822  &
            - t5825 - t5828 - t5830
            t5833 = 0.7200000000e2_dp * t1821 * t5699
            t5835 = 0.3600000000e2_dp * t1821 * t5671
            t5838 = 0.1080000000e3_dp * t1143 * t474 * t4797
            t5841 = 0.5400000000e2_dp * t1143 * t474 * t4801
            t5845 = 0.1800000000e2_dp * t473 * t20 * t5122 * t237
            t5848 = 0.5400000000e2_dp * t1143 * t1150 * t4701
            t5850 = t473 * t4201 * t561
            t5854 = 0.1200000000e2_dp * t473 * t489 * t5371
            t5857 = 0.6000000000e1_dp * t478 * t489 * t5430
            t5861 = 0.5999999999e1_dp * t478 * t56 * t5122 * t17
            t5864 = 0.1800000000e2_dp * t473 * t1150 * t4729
            t5868 = t5833 + t5835 + t5838 + t5841 - t5845 + t5848 - 0.3600000000e2_dp  &
            * t5850 - t5854 + t5857 + t5861 - t5864 + 0.8999999998e1_dp * t178  &
            * t62 * t5613
            t5869 = t5831 + t5868
            t5875 = t1447 * t187 * t1722
            t5880 = 0.3703703705e-1_dp * t192 * t194 * t5226 * t526
            t5881 = t1141 * t4701
            t5885 = t719 * t1769
            t5888 = t60 * t1836 * t709
            t5891 = 0.5555555558e-1_dp * t4943 * t534
            t5898 = t705 * t1504 * t65
            t5902 = 0.1800000000e2_dp * t1718 * t178 * t710
            t5905 = 0.5555555558e-1_dp * t192 * t5270 * t267
            t5908 = t1836 * t59 * t177 * t710
            t5911 = t266 * t4729 * t65
            t5913 = 0.1000000001e1_dp * t723 * t5911
            t5918 = t140 * t1504 * t1722
            t5922 = 0.1111111112e0_dp * t4776 * t515 * t517
            t5925 = 0.8099999996e2_dp * t498 * t4770 * t4915
            t5926 = -0.5555555558e-1_dp * t192 * t68 * t5869 * t65 - 0.1000000000e1_dp  &
            * t5875 + t5880 + 0.3000000000e1_dp * t1748 * t5881 * t4940 - 0.1111111112e0_dp  &
            * t5885 - 0.1000000000e1_dp * t5888 - t5891 + t5869 * t65 +  &
            0.666666666e0_dp * t1519 * t722 * t4729 * t65 + (2 * t5898) + t5902  &
            - t5905 + 0.1800000000e2_dp * t5908 - t5913 - 0.1111111112e0_dp * t516  &
            * t4898 * t172 - 0.1000000000e1_dp * t5918 - t5922 + t5925
            t5929 = 0.2250000000e2_dp * t1740 * t471 * t4905
            t5932 = 0.1000000000e1_dp * t649 * t266 * t1722
            t5938 = t192 * t728 * t1505
            t5940 = t4918 * t67
            t5942 = 0.1111111112e0_dp * t5940 * t196
            t5944 = t1447 * t722 * t1722
            t5949 = 0.8999999998e1_dp * t1795 * t4909
            t5952 = 0.9000000007e1_dp * t1763 * t5881 * t4940
            t5954 = 0.8999999998e1_dp * t1741 * t4909
            t5956 = t471 * t4701 * t65
            t5958 = 0.1000000000e1_dp * t1766 * t5956
            t5960 = 0.5555555558e-1_dp * t5940 * t200
            t5963 = 0.1000000000e1_dp * t379 * t471 * t4904
            t5964 = t1715 * t1525
            t5966 = t1785 * t1790
            t5971 = 0.4500000000e1_dp * t379 * t1141 * t4701 * t4940
            t5973 = 0.5555555558e-1_dp * t4943 * t531
            t5975 = t5226 * t266 * t65
            t5977 = 0.1111111112e0_dp * t1755 * t729
            t5978 = -t5929 - t5932 + 0.3333333334e0_dp * t1519 * t187 * t4729 * t65  &
            - 0.1111111112e0_dp * t5938 - t5942 - 0.2000000002e1_dp * t5944 - 0.5555555558e-1_dp  &
            * t503 * t4955 + t5949 - t5952 + t5954 + t5958 - t5960  &
            + t5963 - 0.1111111112e0_dp * t5964 - 0.2222222224e0_dp * t5966 - t5971  &
            - t5973 + t5975 - t5977
            t5985 = 0.5000000001e0_dp * t1774 * t4908
            t5990 = t177 * t4729 * t65
            t5994 = 0.5000000001e0_dp * t1766 * t5990
            t5996 = 0.3703703705e-1_dp * t4943 * t527
            t6000 = 0.9000000007e1_dp * t5262 * t266 * t4701 * t4940
            t6002 = 0.1851851853e0_dp * t4943 * t512
            t6004 = t1709 * t178 * t710
            t6006 = t719 * t1838
            t6009 = 0.5000000001e0_dp * t1777 * t4908
            t6017 = 0.5000000001e0_dp * t732 * t5911
            t6019 = 0.7407407410e-1_dp * t719 * t1726
            t6021 = t649 * t1518 * t1722
            t6023 = t5251 * t67
            t6025 = 0.3703703705e-1_dp * t6023 * t271
            t6029 = 0.5999999993e1_dp * t4946 * t722 * t4701 * t4940 - t5985 - 0.5555555558e-1_dp  &
            * t192 * t199 * t4898 + 0.3333333334e0_dp * t1748 * t5990  &
            - t5994 + t5996 - t6000 + t6002 + 0.1800000000e2_dp * t6004 - 0.1111111112e0_dp  &
            * t6006 - t6009 + 0.3000000000e1_dp * t4946 * t187 * t4701  &
            * t4940 - 0.6666666668e0_dp * t1748 * t5956 - t6017 + t6019 + 0.6666666668e0_dp  &
            * t6021 + t6025 - 0.5555555558e-1_dp * t58 * t5613 * t71
            t6031 = 0.5000000001e0_dp * t1720 * t5990
            t6033 = 0.1000000000e1_dp * t5244 * t1712
            t6036 = 0.2250000000e2_dp * t498 * t472 * t4905
            t6040 = 0.5555555558e-1_dp * t6023 * t275
            t6042 = 0.1000000001e1_dp * t1763 * t5990
            t6044 = 0.2000000002e1_dp * t1763 * t5956
            t6047 = 0.4500000000e1_dp * t1720 * t5881 * t4940
            t6049 = t1688 * t177 * t709
            t6052 = t1504 * t561 * t65
            t6053 = t732 * t6052
            t6056 = 0.5555555558e-1_dp * t508 * t5271
            t6058 = 0.1111111112e0_dp * t4943 * t520
            t6060 = t1766 * t5881 * t4940
            t6063 = 0.1000000000e1_dp * t1720 * t5956
            t6066 = t1519 * t705 * t561 * t65
            t6074 = 0.8099999996e2_dp * t496 * t4769 * t1141 * t4915
            t6075 = t723 * t6052
            t6078 = 0.1111111112e0_dp * t719 * t1782
            t6079 = -t6031 - t6033 - t6036 + t187 * t4897 * t65 - t6040 - t6042  &
            + t6044 - t6047 - 0.1000000000e1_dp * t6049 - 0.1000000000e1_dp * t6053  &
            - t6056 - t6058 - 0.4500000000e1_dp * t6060 + t6063 + 0.6666666668e0_dp  &
            * t6066 - 0.1000000000e1_dp * t237 * t705 * t1722 + t6074 - 0.2000000002e1_dp  &
            * t6075 - t6078
            t6087 = 0.1666666667e0_dp * t429 * t1966 * t1410
            t6090 = 0.3125000000e0_dp * t4025 * t1247 * t4790
            t6091 = 0.3333333334e0_dp * t5745 - t5749 - t5752 - t5756 + t5760 + 0.1666666667e0_dp  &
            * t165 * t56 * t4958 * t140 + t5767 + t5770 + t5773 +  &
            t5776 + t5780 - t5783 + t5786 + t5789 + t5792 - t5795 - 0.4444444444e0_dp  &
            * t5797 + 0.3333333334e0_dp * t55 * t57 * (t5926 + t5978 + t6029  &
            + t6079) - t6087 - t6090
            t6094 = 0.1250000000e0_dp * t1981 * t1247 * t4801
            t6099 = 0.1800000000e2_dp * t1971 * t1972 * t4729 * t149 * t160
            t6102 = 0.2500000000e0_dp * t1981 * t4553 * t1410
            t6106 = 0.3600000000e2_dp * t1971 * t1972 * t1973 * t672
            t6109 = 0.1666666667e0_dp * t1246 * t1962 * t1865
            t6112 = 0.5555555556e-1_dp * t1246 * t4553 * t1483
            t6116 = 0.3600000000e2_dp * t1971 * t471 * t20 * t5713
            t6119 = 0.1250000000e0_dp * t1368 * t462 * t4790
            t6123 = 0.1666666667e0_dp * t1246 * t446 * t1825 * t561
            t6126 = 0.1800000000e2_dp * t658 * t4764 * t425
            t6129 = 0.8333333335e-1_dp * t429 * t437 * t5174
            t6132 = 0.1666666667e0_dp * t165 * t462 * t4853
            t6136 = 0.5555555557e-1_dp * t165 * t171 * t5122 * t17
            t6138 = t429 * t1928 * t593
            t6141 = t658 * t575 * t1925
            t6144 = t445 * t4512 * t668
            t6155 = 0.2500000000e0_dp * t1981 * t4032 * t1410
            t6158 = 0.1666666667e0_dp * t429 * t679 * t1571
            t6160 = t147 * t1885 * t672
            t6162 = t6094 + t6099 + t6102 + t6106 - t6109 + t6112 - t6116 + t6119  &
            - t6123 + t6126 - t6129 + t6132 - t6136 - 0.1666666667e0_dp * t6138  &
            + 0.3600000000e2_dp * t6141 + 0.3333333334e0_dp * t6144 + 0.2e1_dp * t147  &
            * (t4881 - t4884 - t4887 + 0.9000000000e1_dp * t4889 - t4893 + 0.6000000000e1_dp  &
            * t4895) * t149 * t160 + t6155 - t6158 + (4 * t6160)
            t6170 = 0.3333333336e0_dp * t289 * t597 * t1571
            t6173 = 0.3333333336e0_dp * t55 * t5295 * t140
            t6175 = t83 * t597 * t1688
            t6179 = 0.3333333336e0_dp * t165 * t3865 * t561
            t6181 = t289 * t1564 * t593
            t6185 = 0.3333333336e0_dp * t289 * t597 * t1553
            t6186 = -0.3333333336e0_dp * t5347 - t5351 - t5354 - t5357 - t5360 - 0.2222222224e0_dp  &
            * t144 * t46 * (t5423 + t5742 + t6091 + t6162) - t6170  &
            - t6173 - 0.6666666672e0_dp * t6175 - t6179 - 0.3333333336e0_dp * t6181  &
            - t6185
            t6189 = 0.8333333340e-1_dp * t906 * t294 * t5163
            t6192 = 0.3333333336e0_dp * t83 * t294 * t5122
            t6194 = t83 * t3875 * t561
            t6198 = 0.3333333336e0_dp * t289 * t294 * t5170
            t6201 = 0.3333333336e0_dp * t289 * t597 * t1557
            t6207 = 0.1666666668e0_dp * t1546 * t897 * t4797
            t6210 = 0.1666666668e0_dp * t289 * t294 * t5174
            t6213 = 0.3333333336e0_dp * t55 * t84 * t5122
            t6216 = 0.8333333340e-1_dp * t429 * t897 * t4701
            t6219 = 0.8333333340e-1_dp * t906 * t384 * t4790
            t6222 = 0.8333333340e-1_dp * t1546 * t897 * t4801
            t6223 = t6189 - t6192 - 0.6666666672e0_dp * t6194 - t6198 - t6201 - 0.3333333336e0_dp  &
            * t83 * t84 * t5613 + t6207 - t6210 - t6213 + t6216 +  &
            t6219 + t6222
            t6227 = -t4704 - t4707 - t4732 - t4789 + t4793 - t4796 - t4800 - t4804  &
            - t4808 - t4812 - t4857 - t4861 - 0.1333333334e1_dp * t4864 - t4869  &
            - 0.4444444448e0_dp * t43 * t46 * t5036 + t8 * (t5305 + t5343 + t6186  &
            + t6223)
            e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii) +  t6227 * sx
            t6247 = t4769 * t2992
            t6249 = t4772 * t2027
            t6252 = t56 * t4701
            t6253 = t6252 * t526
            t6256 = t56 * t4729
            t6257 = t6256 * t526
            t6260 = t4743 * t526
            t6263 = t20 * t4853
            t6267 = t26 * t1447
            t6268 = t6267 * t709
            t6271 = t4764 * t2027
            t6281 = -0.1620000000e3_dp * t6247 * t4733 * t6249 - 0.1200000000e2_dp *  &
            t473 * t6253 + 0.5999999999e1_dp * t478 * t6257 + 0.5399999998e2_dp * t4771  &
            * t6260 + 0.8999999998e1_dp * t478 * t6263 * t65 + 0.8099999996e2_dp  &
            * t4771 * t6268 - 0.1800000000e2_dp * t473 * t6271 - 0.4050000000e2_dp  &
            * t473 * t575 * t2095 - 0.5000000001e0_dp * t60 * t4853 * t65 + 0.8999999998e1_dp  &
            * t4901 - t4907 + t4911
            t6283 = t4701 * t1518
            t6290 = t4917 - t4920 + 0.5000000001e0_dp * t4922 + t4925 + 0.3333333334e0_dp  &
            * t6283 * t177 * t65 - 0.5000000001e0_dp * t4929 - t4932 - t4935  &
            + 0.5850000000e2_dp * t4941 + t4945 - t4951 - 0.5555555558e-1_dp * t4952
            t6291 = t6281 + t6290
            t6308 = 0.4050000000e2_dp * t4734 * t4736 * t2821 * t18 * t4739 * t237  &
            - 0.1125000000e2_dp * t2013 * t1300 * t19 * t4790 - t5028 - t5022 +  &
            0.1500000000e1_dp * t2013 * t398 * t68 * t5002 + 0.3333333334e0_dp * t55  &
            * t57 * t6291 + t5019 - 0.3000000000e1_dp * t566 * t1885 * t569 - t5016  &
            - 0.9000000000e1_dp * t4734 * t4736 * t1300 * t4701 * t149 * t17  &
            + t5032 + 0.1666666667e0_dp * t165 * t2088 * t561 - 0.5555555557e-1_dp *  &
            t5011 + t5001
            t6338 = 0.1666666667e0_dp * t4874 - t5025 - 0.2700000000e2_dp * t4734 * t392  &
            * t4736 * t1300 * my_rho * t5002 * t149 - 0.3000000000e1_dp * t2009  &
            * t4748 - 0.2700000000e2_dp * t4738 * t18 * t1973 * t1447 - t5008 -  &
            0.3000000000e1_dp * t567 * t151 * t19 * t4853 + t5035 - 0.1000000000e1_dp  &
            * t2013 * t2018 * t5009 + 0.9000000000e1_dp * t2013 * t2014 * t4797  &
            - t4872 + 0.4500000000e1_dp * t2013 * t2014 * t4801 + t4878 + t5005  &
            + 0.4500000000e1_dp * t2009 * t4744
            t6339 = t6308 + t6338
            t6344 = -t5043 - t5047 - t5050 - t5126 - t5284 - t5291 + t5294 - t5298  &
            - t5301 - 0.3333333336e0_dp * t5303 - t5308 - t5311 - t5314
            t6352 = t5440 - t5443 - t5446 - t5449 - t5451 - t5454 + 0.40e1_dp * t102  &
            * t327 * t2057 * t557 - t5456 + t5459 - t5462 + t5466 - t5468
            t6363 = t5480 - t5482 - t5484 + t5487 - t5489 + t5492 + t5494 - 0.9600000000e2_dp  &
            * t5496 + t5499 - t5503 - t5507 + t5509 - 0.2400000000e2_dp  &
            * t5511 - 0.960e2_dp * t5115 * t1060 - 0.240e2_dp * t2054 * t3674 + 0.1200e3_dp  &
            * t1054 * t640 * t3679
            t6367 = t5472 - t5474 - t5478 - t5523 + t5525 + 0.40e1_dp * t5528 + t5531  &
            + t5533 + t5535 + t5537 - 0.20e1_dp * t102 * t105 * t6363 + t5540
            t6370 = t5545 - t5548 + t5551 + t5553 - t5558 + t5560 - t5562 + t5564  &
            - 0.40e1_dp * t5565 + t5568 + t5572 + t5574
            t6373 = t5579 - t5581 - 0.40e1_dp * t112 * t2058 - t5585 + t5590 - t5593  &
            - t5595 + t5597 - t5600 - t5602 + t5604 + t5607 + t5610
            t6375 = t6352 + t6367 + t6370 + t6373
            t6380 = -t5317 - t5320 + t5323 - t5326 - t5329 - 0.3333333336e0_dp * t83  &
            * t597 * t2061 - 0.3333333336e0_dp * t83 * t84 * t6375 + t5333 - t5335  &
            - t5338 - t5342 - 0.1666666668e0_dp * t5347 - t5351
            t6389 = -t5354 - t5357 - 0.3333333336e0_dp * t83 * t45 * t2139 * t561  &
            - t5360 - t6170 - t6173 - 0.3333333336e0_dp * t6175 - t6179 - 0.1666666668e0_dp  &
            * t6181 - t6185 + t6189 - t6192 - 0.3333333336e0_dp * t6194
            t6400 = t2061 * t561
            t6415 = -t5363 + t5367 - t5370 + 0.1111111112e0_dp * t5373 - 0.5555555557e-1_dp  &
            * t5376 + 0.1666666667e0_dp * t5379 - t5383 - t5386 - t5389 +  &
            0.1200000000e2_dp * t658 * t2031 * t673 * t17 + 0.3333333334e0_dp * t170  &
            * t27 * t6291 * t172 - t5397
            t6428 = t5635 - t5638 + t5641 + 0.2250000000e1_dp * t5643 - t5647 - 0.1500000000e1_dp  &
            * t151 * t50 * t6375 - t5650 + t5653 - t5657 - t5660  &
            + t5663 - t5666
            t6433 = t5670 - t5673 + t5676 - t5679 + t5682 + t5686 - t5689 + t5692  &
            + t5695 + 0.2250000000e1_dp * t400 * t19 * t2061 * t561 + t5698 -  &
            t5701
            t6442 = t5400 + 0.1666666667e0_dp * t5402 - t5406 + 0.1666666667e0_dp * t165  &
            * t56 * t2131 * t561 + 0.1666666667e0_dp * t165 * t679 * t2061 +  &
            (2 * t147 * t149 * (t6428 + t6433)) - t5409 - t5412 - t5415 +  &
            t5418 + 0.1800000000e2_dp * t658 * t575 * t2074 + 0.1800000000e2_dp * t5421
            t6451 = -t5426 + t5429 - 0.2222222223e0_dp * t5432 + t5436 + 0.1620000000e3_dp  &
            * t4995 * t6267 * t1974 + 0.1666666667e0_dp * t165 * t2088 * t649  &
            + t5619 + t5622 - t5625 - t5628 + t5631 + t5716
            t6461 = t5719 - t5722 - 0.8333333335e-1_dp * t5724 + 0.1200000000e2_dp *  &
            t658 * t6256 * t2077 - 0.3240000000e3_dp * t4994 * t6247 * t4772 * t2082  &
            + t5728 + t5731 - t5734 - t5738 + t5741 + 0.1666666667e0_dp * t5745  &
            - t5749
            t6479 = -t5752 + 0.1666666667e0_dp * t445 * t2042 * my_rho * t668 - t5756  &
            - 0.3600000000e2_dp * t2081 * t575 * t673 * t237 + 0.1800000000e2_dp *  &
            t658 * t6263 * t161 + t5760 + t5767 + t5770 - 0.7200000000e2_dp * t2081  &
            * t575 * t161 * t1447 + t5773 + t5776 + t5780
            t6490 = -t5783 - 0.3600000000e2_dp * t2081 * t4764 * t2082 - 0.2400000000e2_dp  &
            * t2081 * t6252 * t2077 + t5786 + t5789 + t5792 - t5795 - 0.4444444445e0_dp  &
            * t5797 - t6087 - 0.8333333335e-1_dp * t429 * t166 * t6400  &
            - t6090 + t6094
            t6499 = t6099 + 0.1666666667e0_dp * t165 * t166 * t6375 + t6102 + 0.1666666667e0_dp  &
            * t165 * t56 * t6291 * t140 + t6106 - t6109 + t6112 - t6116  &
            + t6119 - t6123 + t6126 - t6129
            t6509 = t4853 * t65
            t6529 = 0.4999999999e0_dp * t5875 - 0.5000000001e0_dp * t5244 * t2095 - 0.5000000001e0_dp  &
            * t732 * t6509 + t5880 + 0.1000000000e1_dp * t649 * t471  &
            * t2099 - 0.1620000000e3_dp * t5247 * t2992 * t4733 * t6249 - 0.5555555558e-1_dp  &
            * t5885 + 0.8999999998e1_dp * t140 * t2992 * t4701 * t59 *  &
            t62 * t2027 - 0.5000000001e0_dp * t5888 - t5891 - 0.1800000000e2_dp * t2118  &
            * t6271
            t6539 = t5898 + t5902 - t5905 + 0.8999999998e1_dp * t5908 - t5913 + 0.2000000000e1_dp  &
            * t2098 * t2095 * t561 - 0.5000000001e0_dp * t5918 - t5922  &
            + 0.3333333334e0_dp * t6283 * t187 * t177 * t65 + t5925 - t5929
            t6554 = -t5932 - 0.5555555558e-1_dp * t5938 + 0.666666666e0_dp * t6283 *  &
            t722 * t177 * t65 - t5942 + 0.8999999998e1_dp * t2124 * t59 * t177 *  &
            t710 + 0.1000000001e1_dp * t5944 - 0.5000000001e0_dp * t2061 * t177 * t709  &
            + t5949 + t5952 + t5954 - t5958
            t6557 = t1447 * t59
            t6569 = -t5960 - 0.5000000001e0_dp * t5237 * t2095 + t5963 - 0.4500000000e1_dp  &
            * t5256 * t6557 * t710 - 0.5555555558e-1_dp * t5964 - 0.1111111112e0_dp  &
            * t5966 - t5971 - 0.1200000000e2_dp * t2118 * t6253 - 0.5555555558e-1_dp  &
            * t719 * t2126 - t5973 + t5975 + 0.5999999999e1_dp * t2114 * t6257
            t6583 = -t5977 - t5985 + t5994 - 0.5000000001e0_dp * t716 * t6509 + t5996  &
            - t6000 + t6002 + 0.8999999998e1_dp * t6004 + 0.1000000000e1_dp * t2098  &
            * t4908 * t237 - 0.1000000001e1_dp * t723 * t6509 - 0.1800000000e2_dp  &
            * t5228 * t2117 * t2028
            t6596 = t5247 * t1141 * t4733
            t6603 = -0.5555555558e-1_dp * t6006 - t6009 - t6017 - 0.5555555558e-1_dp  &
            * t58 * t6375 * t71 + t6019 - 0.9000000007e1_dp * t5262 * t6557 * t710  &
            + 0.8999999998e1_dp * t708 * t62 * t6509 + 0.3333333334e0_dp * t6021 +  &
            t6025 + 0.8099999996e2_dp * t6596 * t6268 + 0.5399999998e2_dp * t6596 *  &
            t6260 + 0.8999999998e1_dp * t5229 * t2104
            t6611 = t5801 - t5804 - t5807 + t5810 - t5813 + t5816 - t5819 - t5822  &
            - t5825 - t5828 - t5830 + t5833
            t6620 = t5835 + t5838 + t5841 - t5845 + t5848 - 0.1800000000e2_dp * t5850  &
            - t5854 + t5857 + t5861 - t5864 - 0.1800000000e2_dp * t473 * t20  &
            * t2061 * t561 + 0.8999999998e1_dp * t178 * t62 * t6375
            t6621 = t6611 + t6620
            t6630 = -0.3000000000e1_dp * t5256 * t4904 * t237 + 0.5999999999e1_dp * t5228  &
            * t2113 * t2032 + t6621 * t65 - t6031 - 0.4050000000e2_dp * t471  &
            * t187 * t6557 * t710 - t6033 - t6036 - t6040 + t6042 - 0.5000000001e0_dp  &
            * t60 * t2124 * t709 - t6044
            t6640 = -t6047 - 0.5000000001e0_dp * t6049 - 0.5555555558e-1_dp * t192 *  &
            t68 * t6621 * t65 - 0.5000000001e0_dp * t6053 - t6056 - t6058 + 0.5850000000e2_dp  &
            * t6060 + t6063 + 0.3333333334e0_dp * t6066 + t6074 - 0.1000000001e1_dp  &
            * t6075 - t6078
            t6658 = 0.1800000000e2_dp * t658 * t2035 * t673 - 0.8333333335e-1_dp * t429  &
            * t2088 * t593 + t6132 + 0.3333333334e0_dp * t55 * t57 * (t6529 +  &
            t6539 + t6554 + t6569 + t6583 + t6603 + t6630 + t6640) - t6136 - 0.8333333335e-1_dp  &
            * t6138 + 0.1800000000e2_dp * t6141 + 0.1080000000e3_dp *  &
            t4995 * t4743 * t2077 + 0.1080000000e3_dp * t657 * t1142 * t4760 * t2082  &
            + 0.1666666667e0_dp * t6144 + t6155 - t6158 + (2 * t6160)
            t6665 = -0.3333333336e0_dp * t83 * t45 * t6339 * t140 - t6198 - 0.3333333336e0_dp  &
            * t83 * t2050 * t649 - t6201 + t6207 - t6210 - t6213 + t6216  &
            + t6219 - 0.1666666668e0_dp * t289 * t2050 * t593 - 0.1666666668e0_dp  &
            * t289 * t84 * t6400 + t6222 - 0.2222222224e0_dp * t144 * t46 * (t6415  &
            + t6442 + t6451 + t6461 + t6479 + t6490 + t6499 + t6658)
            t6669 = -t4704 - t4707 - t4732 + t4793 - t4796 - t4800 - t4804 - t4861  &
            - 0.6666666672e0_dp * t4864 - t4869 - t4812 - t4857 - t4789 - t4808  &
            - 0.6666666672e0_dp * t217 * t46 * t2046 * t561 - 0.4444444448e0_dp *  &
            t43 * t46 * t6339 + t8 * (t6344 + t6380 + t6389 + t6665)
            e_ndrho_ndrho_ndrho(ii) = e_ndrho_ndrho_ndrho(ii) +  t6669 * sx
          END IF

        END IF
      END DO
      !$omp end do
    END IF
  END SUBROUTINE xb88_lr_adiabatic_lda_calc

! *****************************************************************************
!> \brief evaluates the becke 88 longrange exchange functional for lsd
!> \param sx scaling-parameter for exchange
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_adiabatic_lsd_eval(rho_set,deriv_set,grad_deriv,xb88_lr_ad_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xb88_lr_ad_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_adiabatic_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, npoints, &
                                                stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, lambda, omega, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0
    TYPE(cp_3d_r_p_type), DIMENSION(2) :: e_ndrho, e_ndrho_ndrho, &
      e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, e_ndrho_rho_rho, &
      e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    CALL cite_reference(Becke1988)

    failure=.FALSE.
    NULLIFY(deriv, bo)
    DO i=1,2
       NULLIFY(norm_drho(i)%array, rho(i)%array)
    END DO

    CALL section_vals_val_get(xb88_lr_ad_params,"SCALE_X",r_val=sx,error=error)
    CALL section_vals_val_get(xb88_lr_ad_params, "OMEGA",r_val=omega, error=error)
    CALL section_vals_val_get(xb88_lr_ad_params, "LAMBDA",r_val=lambda, error=error)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,&
            rhoa=rho(1)%array,&
            rhob=rho(2)%array,norm_drhoa=norm_drho(1)%array, &
            norm_drhob=norm_drho(2)%array,rho_cutoff=epsilon_rho,&
            local_bounds=bo, error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       ! meaningful default for the arrays we don't need: let us make compiler
       ! and debugger happy...
       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> rho(1)%array
       END IF

       e_0=>dummy
       DO i=1,2
          e_rho(i)%array => dummy
          e_ndrho(i)%array => dummy
          e_rho_rho(i)%array => dummy
          e_ndrho_rho(i)%array => dummy
          e_ndrho_ndrho(i)%array => dummy
          e_rho_rho_rho(i)%array => dummy
          e_ndrho_rho_rho(i)%array => dummy
          e_ndrho_ndrho_rho(i)%array => dummy
          e_ndrho_ndrho_ndrho(i)%array => dummy
       END DO

       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho(2)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho(2)%array,error=error)
       END IF
       IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho(2)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho(2)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho(2)%array,error=error)
       END IF
       IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho(2)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhoa)(rhoa)(rhoa)",allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhob)(rhob)(rhob)",allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho(2)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhoa)(norm_drhoa)(rhoa)",allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhob)(norm_drhob)(rhob)",allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho(2)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhoa)(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,&
               error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho(1)%array,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhob)(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,&
               error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho(2)%array,error=error)
       END IF
       IF (grad_deriv>3.OR.grad_deriv<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       DO ispin=1,2

          !$omp parallel default(none) &
          !$omp          shared(rho, ispin, norm_drho, e_0, e_rho, e_ndrho) &
          !$omp          shared(e_rho_rho, e_ndrho_rho, e_ndrho_ndrho) &
          !$omp          shared(e_rho_rho_rho, e_ndrho_rho_rho) &
          !$omp          shared(e_ndrho_ndrho_rho, e_ndrho_ndrho_ndrho) &
          !$omp          shared(grad_deriv, npoints, epsilon_rho, sx, omega) &
          !$omp          shared(lambda, error)

          CALL xb88_lr_adiabatic_lsd_calc(&
               rho_spin=rho(ispin)%array,&
               norm_drho_spin=norm_drho(ispin)%array,&
               e_0=e_0,&
               e_rho_spin=e_rho(ispin)%array,&
               e_ndrho_spin=e_ndrho(ispin)%array,&
               e_rho_rho_spin=e_rho_rho(ispin)%array,&
               e_ndrho_rho_spin=e_ndrho_rho(ispin)%array,&
               e_ndrho_ndrho_spin=e_ndrho_ndrho(ispin)%array,&
               e_rho_rho_rho_spin=e_rho_rho_rho(ispin)%array,&
               e_ndrho_rho_rho_spin=e_ndrho_rho_rho(ispin)%array,&
               e_ndrho_ndrho_rho_spin=e_ndrho_ndrho_rho(ispin)%array,&
               e_ndrho_ndrho_ndrho_spin=e_ndrho_ndrho_ndrho(ispin)%array,&
               grad_deriv=grad_deriv, npoints=npoints,&
               epsilon_rho=epsilon_rho, &
               sx=sx, omega=omega, lambda=lambda, &
               error=error)

          !$omp end parallel

       END DO

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE xb88_lr_adiabatic_lsd_eval

! *****************************************************************************
!> \brief low level calculation of the becke 88 exchange functional for lsd
!> \param rho_spin alpha or beta spin density
!> \param norm_drho_spin || grad rho_spin ||
!> \param e_ 0: adds to it the local value of the functional
!> \param e_ *_spin: derivative of the functional wrt. to the variables
!>        named where the * is. Everything wrt. to the spin of the arguments.
!> \param sx scaling-parameter for exchange
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_adiabatic_lsd_calc(rho_spin, norm_drho_spin,e_0,&
       e_rho_spin,e_ndrho_spin,e_rho_rho_spin,e_ndrho_rho_spin,&
       e_ndrho_ndrho_spin, e_rho_rho_rho_spin, e_ndrho_rho_rho_spin,&
       e_ndrho_ndrho_rho_spin,&
       e_ndrho_ndrho_ndrho_spin,grad_deriv,npoints,epsilon_rho, sx,&
       omega, lambda, error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho_spin, norm_drho_spin
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_rho_spin, &
      e_ndrho_spin, e_rho_rho_spin, e_ndrho_rho_spin, e_ndrho_ndrho_spin, &
      e_rho_rho_rho_spin, e_ndrho_rho_rho_spin, e_ndrho_ndrho_rho_spin, &
      e_ndrho_ndrho_ndrho_spin
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, omega, lambda
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_adiabatic_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: Cx, epsilon_rho43, my_epsilon_rho, ndrho, rho, t1, t10, &
      t1000, t1002, t1004, t1008, t1012, t1014, t1019, t102, t1021, t1024, &
      t1026, t103, t1030, t1033, t1034, t1036, t104, t1040, t1044, t1046, &
      t105, t1050, t1054, t1055, t1056, t1058, t1059, t106, t1060, t1063, &
      t1064, t1067, t1068, t107, t1071, t1076, t1079, t108, t1082, t1085, &
      t1089, t1093, t1097, t1099, t11, t1100, t1101, t1106, t1109, t1111, &
      t1119, t112, t1122, t1125, t1126, t1129, t1134, t1137, t114, t1140, &
      t1141, t1142, t1143, t115, t1150, t1157, t116, t1161, t1168, t1172, &
      t1176, t1189, t119, t1191, t1192, t1193, t1194, t12, t1200, t1201, t1206
    REAL(kind=dp) :: t1207, t1210, t1213, t1216, t122, t1221, t1224, t123, &
      t1235, t1238, t1239, t1246, t1247, t1254, t126, t1264, t1268, t1275, &
      t1283, t1284, t1288, t1292, t13, t130, t1300, t1301, t1302, t1309, &
      t131, t1316, t132, t1320, t133, t1333, t1340, t1347, t1348, t136, &
      t1360, t1367, t1368, t137, t1372, t1373, t1377, t1387, t1388, t1392, &
      t1396, t1397, t14, t140, t1401, t1406, t1408, t1409, t141, t1410, &
      t1411, t1413, t1417, t1420, t1425, t1426, t1429, t1436, t1438, t1439, &
      t144, t1440, t1443, t1444, t1447, t145, t1451, t1455, t1459, t146, &
      t1466, t147, t1470, t1471, t1476, t1479, t1482, t1483, t1485, t1489, &
      t149
    REAL(kind=dp) :: t1492, t1493, t1496, t1498, t15, t150, t1500, t1502, &
      t1504, t1505, t1506, t1507, t1509, t151, t1510, t1512, t1513, t1515, &
      t1517, t1518, t1519, t1523, t1525, t1528, t1530, t1531, t1535, t1536, &
      t1541, t1543, t1545, t1546, t1547, t1549, t155, t1552, t1553, t1556, &
      t1557, t156, t1560, t1563, t1564, t1565, t1568, t1570, t1571, t1574, &
      t1575, t1577, t1579, t1581, t1582, t1584, t1586, t1588, t1589, t1591, &
      t1592, t1594, t1596, t1597, t1598, t1599, t16, t160, t1601, t1602, &
      t1603, t1605, t1607, t1608, t161, t1610, t1611, t1613, t1615, t1616, &
      t1617, t1618, t1620, t1621, t1622, t1624, t1625, t1626, t1628, t1630
    REAL(kind=dp) :: t1631, t1632, t1634, t1636, t1638, t164, t1640, t1642, &
      t1644, t1645, t1647, t1649, t165, t1652, t1654, t1656, t1658, t166, &
      t1660, t1662, t1664, t1667, t1669, t1671, t1675, t1677, t1679, t1680, &
      t1681, t1684, t1685, t1688, t1689, t1693, t1694, t1696, t1699, t17, &
      t170, t1700, t1702, t1703, t1705, t1706, t1708, t1709, t171, t1712, &
      t1714, t1715, t1717, t1718, t1719, t172, t1720, t1722, t1724, t1726, &
      t1728, t1730, t1731, t1735, t1737, t1740, t1741, t1743, t1745, t1748, &
      t1752, t1754, t1755, t1757, t1758, t176, t1760, t1762, t1763, t1765, &
      t1766, t1768, t1769, t177, t1773, t1774, t1776, t1777, t1779, t178
    REAL(kind=dp) :: t1782, t1784, t1785, t1787, t1789, t1790, t1794, t1795, &
      t1797, t18, t1800, t1801, t1804, t1807, t1810, t1811, t1814, t1817, &
      t182, t1821, t1822, t1824, t1825, t1828, t1829, t183, t1832, t1835, &
      t1836, t1837, t1838, t1842, t1843, t1849, t1852, t1854, t1858, t1861, &
      t1864, t1865, t1868, t187, t1871, t1874, t1875, t1878, t1881, t1884, &
      t1885, t1886, t189, t1891, t1892, t1895, t1898, t19, t1901, t1902, &
      t1905, t1908, t1912, t1914, t1917, t192, t1920, t1923, t1924, t1925, &
      t1928, t1932, t1938, t1939, t194, t1940, t1943, t1946, t195, t1952, &
      t1954, t1958, t196, t1961, t1962, t1965, t1966, t1969, t1971, t1972
    REAL(kind=dp) :: t1973, t1974, t1975, t1977, t1980, t1981, t1984, t1987, &
      t1988, t199, t1991, t1994, t1997, t20, t200, t2000, t2001, t2002, &
      t2006, t2009, t2013, t2014, t2018, t2022, t2027, t2028, t203, t2031, &
      t2032, t2035, t2042, t2046, t2050, t2054, t2057, t2058, t2061, t207, &
      t2073, t2074, t2077, t2081, t2082, t2086, t2088, t2095, t2098, t2099, &
      t2104, t2107, t2113, t2114, t2117, t2118, t2124, t2126, t2130, t2131, &
      t2138, t2139, t214, t2143, t2147, t2160, t2169, t217, t218, t2182, &
      t2183, t2187, t219, t2190, t2193, t2196, t2197, t22, t220, t2200, &
      t2201, t2204, t2205, t221, t2220, t2226, t2227, t2230, t224, t225, t2257
    REAL(kind=dp) :: t2261, t2264, t2268, t2272, t2274, t2276, t2279, t2280, &
      t2281, t2282, t2285, t2287, t2292, t2295, t2296, t2299, t23, t230, &
      t2300, t2304, t2307, t2308, t2311, t2312, t2319, t2326, t233, t234, &
      t2351, t2368, t237, t2378, t2382, t2392, t24, t2404, t2418, t2424, &
      t2437, t2439, t245, t2451, t2455, t2468, t2480, t249, t2492, t25, t250, &
      t2537, t2552, t2565, t2569, t2576, t2585, t26, t260, t2611, t262, &
      t2621, t2622, t2631, t2635, t2636, t264, t2645, t2655, t266, t267, &
      t268, t27, t2707, t271, t2720, t2724, t2732, t274, t275, t2759, t2763, &
      t2779, t278, t2783, t2808, t2816, t282, t2821, t2823, t283, t2833
    REAL(kind=dp) :: t2843, t2853, t2869, t2873, t2877, t2885, t289, t290, &
      t2935, t2937, t294, t295, t298, t299, t2992, t2994, t2998, t302, t3026, &
      t3030, t306, t3069, t307, t3070, t308, t3085, t309, t3093, t3097, &
      t3098, t31, t311, t3129, t313, t3141, t316, t317, t318, t3181, t321, &
      t322, t3232, t3244, t3248, t326, t3262, t3264, t3268, t327, t3278, &
      t328, t3288, t329, t3291, t3305, t3306, t3314, t332, t3329, t3339, &
      t3344, t3347, t3352, t3353, t3356, t3357, t336, t3362, t3363, t3366, &
      t3378, t3379, t3382, t339, t341, t342, t343, t3449, t3453, t347, t3472, &
      t3491, t3494, t3497, t35, t350, t3501, t3505, t3507, t3512, t3513
    REAL(kind=dp) :: t3514, t352, t3525, t3533, t3543, t3546, t3547, t3552, &
      t3553, t356, t3561, t3569, t357, t3575, t3576, t3580, t3581, t36, t360, &
      t3608, t3618, t362, t365, t366, t3668, t367, t3673, t3674, t3677, &
      t3679, t368, t3682, t37, t371, t3711, t372, t3726, t375, t3750, t376, &
      t3763, t379, t3794, t3796, t380, t3826, t3827, t383, t384, t3842, &
      t3849, t3853, t3861, t3865, t3869, t387, t3871, t3875, t3892, t3918, &
      t392, t3925, t393, t3933, t3939, t394, t3940, t3965, t3969, t3975, &
      t3979, t398, t3985, t399, t4, t40, t400, t401, t4010, t4018, t4025, &
      t4032, t405, t406, t4070, t4079, t41, t4101, t4108, t413, t4140, t4156
    REAL(kind=dp) :: t417, t42, t4201, t4214, t4230, t424, t425, t4255, &
      t4256, t4261, t4272, t428, t429, t4292, t43, t4317, t433, t4345, t4352, &
      t437, t4378, t4392, t44, t4428, t444, t445, t446, t4466, t4498, t45, &
      t450, t4505, t451, t4512, t4545, t455, t4553, t4575, t4585, t46, t4607, &
      t462, t4627, t4657, t466, t4696, t47, t470, t4700, t4701, t4704, t4707, &
      t4708, t471, t4713, t4714, t4717, t472, t4721, t4725, t4726, t4729, &
      t473, t4732, t4733, t4734, t4736, t4738, t4739, t474, t4743, t4744, &
      t4748, t4760, t4764, t4768, t4769, t4770, t4771, t4772, t4776, t478, &
      t4782, t4786, t4789, t479, t4790, t4793, t4796, t4797, t48
    REAL(kind=dp) :: t4800, t4801, t4804, t4808, t4812, t4817, t4820, t4829, &
      t4830, t4833, t4834, t4837, t4838, t4849, t4850, t4853, t4857, t4861, &
      t4864, t4869, t4872, t4874, t4878, t4881, t4884, t4887, t4889, t489, &
      t4893, t4895, t4897, t4898, t49, t4901, t4904, t4905, t4907, t4908, &
      t4909, t4911, t4915, t4917, t4918, t4920, t4922, t4925, t4929, t4932, &
      t4935, t4940, t4941, t4943, t4945, t4946, t4951, t4952, t4955, t4958, &
      t496, t498, t4994, t4995, t50, t500, t5001, t5002, t5005, t5008, t5009, &
      t5011, t5013, t5016, t5019, t5022, t5025, t5028, t5029, t503, t5032, &
      t5035, t5036, t5040, t5043, t5044, t5047, t5050, t5051, t5060
    REAL(kind=dp) :: t5065, t5069, t5079, t508, t5083, t5087, t5106, t5108, &
      t511, t5115, t5118, t512, t5122, t5123, t5126, t5139, t5143, t5147, &
      t515, t516, t5163, t517, t5170, t5174, t520, t5226, t5228, t5229, &
      t5237, t5244, t5247, t5251, t5256, t526, t5262, t527, t5270, t5271, &
      t5276, t5280, t5281, t5284, t5291, t5294, t5295, t5298, t53, t530, &
      t5301, t5303, t5305, t5308, t531, t5311, t5314, t5317, t5320, t5323, &
      t5326, t5329, t5330, t5333, t5335, t5338, t5339, t534, t5342, t5343, &
      t5345, t5347, t5351, t5354, t5357, t5360, t5363, t5367, t537, t5370, &
      t5371, t5373, t5376, t5379, t5383, t5386, t5389, t5397, t5400, t5402
    REAL(kind=dp) :: t5406, t5409, t541, t5412, t5415, t5418, t5421, t5423, &
      t5426, t5429, t5430, t5432, t5436, t5440, t5443, t5446, t5449, t5451, &
      t5454, t5456, t5459, t5460, t5462, t5466, t5468, t5469, t5472, t5474, &
      t5478, t5480, t5482, t5484, t5487, t5489, t5492, t5494, t5496, t5499, &
      t55, t550, t5503, t5507, t5509, t551, t5511, t5516, t5523, t5525, &
      t5528, t5531, t5533, t5535, t5537, t5540, t5541, t5545, t5548, t5551, &
      t5553, t5558, t5560, t5562, t5564, t5565, t5568, t557, t5572, t5574, &
      t5575, t5579, t558, t5581, t5585, t5590, t5593, t5595, t5597, t56, &
      t5600, t5602, t5604, t5607, t561, t5610, t5611, t5613, t5619, t5622
    REAL(kind=dp) :: t5625, t5628, t5631, t5632, t5635, t5638, t5641, t5643, &
      t5647, t5650, t5653, t5657, t566, t5660, t5663, t5666, t5667, t567, &
      t5670, t5671, t5673, t5676, t5679, t5682, t5686, t5689, t569, t5692, &
      t5695, t5696, t5698, t5699, t57, t5701, t5705, t5713, t5716, t5719, &
      t5722, t5724, t5728, t5731, t5734, t5738, t5741, t5742, t5745, t5749, &
      t575, t5752, t5756, t5760, t5767, t5770, t5773, t5776, t5780, t5783, &
      t5786, t5789, t579, t5792, t5795, t5797, t58, t5801, t5804, t5807, &
      t5810, t5813, t5816, t5819, t5822, t5825, t5828, t5830, t5831, t5833, &
      t5835, t5838, t5841, t5845, t5848, t585, t5850, t5854, t5857, t5861
    REAL(kind=dp) :: t5864, t5868, t5869, t5875, t5880, t5881, t5885, t5888, &
      t589, t5891, t5898, t59, t590, t5902, t5905, t5908, t5911, t5913, &
      t5918, t5922, t5925, t5926, t5929, t593, t5932, t5938, t5940, t5942, &
      t5944, t5949, t5952, t5954, t5956, t5958, t5960, t5963, t5964, t5966, &
      t597, t5971, t5973, t5975, t5977, t5978, t598, t5985, t5990, t5994, &
      t5996, t6, t60, t6000, t6002, t6004, t6006, t6009, t6017, t6019, t6021, &
      t6023, t6025, t6029, t6031, t6033, t6036, t6040, t6042, t6044, t6047, &
      t6049, t605, t6052, t6053, t6056, t6058, t6060, t6063, t6066, t6074, &
      t6075, t6078, t6079, t608, t6087, t609, t6090, t6091, t6094
    REAL(kind=dp) :: t6099, t61, t6102, t6106, t6109, t6112, t6116, t6119, &
      t6123, t6126, t6129, t6132, t6136, t6138, t614, t6141, t6144, t615, &
      t6155, t6158, t6160, t6162, t6170, t6173, t6175, t6179, t6181, t6185, &
      t6186, t6189, t6192, t6194, t6198, t62, t6201, t6207, t6210, t6213, &
      t6216, t6219, t6222, t6223, t6227, t624, t6247, t6249, t6252, t6253, &
      t6256, t6257, t6260, t6263, t6267, t6268, t627, t6271, t628, t6281, &
      t6283, t6290, t6291, t63, t6308, t6338, t6339, t6344, t635, t6352, &
      t6363, t6367, t637, t6370, t6373, t6375, t6380, t6389, t640, t6400, &
      t6415, t6428, t6433, t6442, t645, t6451, t646, t6461, t6479, t649, t6490
    REAL(kind=dp) :: t6499, t65, t650, t6509, t6529, t653, t6539, t6554, &
      t6557, t6569, t657, t658, t6583, t6596, t66, t6603, t6611, t6620, &
      t6621, t6630, t6640, t6658, t6665, t6669, t668, t67, t672, t673, t679, &
      t68, t689, t69, t7, t70, t705, t707, t708, t709, t71, t710, t713, t716, &
      t719, t722, t723, t728, t729, t732, t735, t739, t74, t747, t753, t754, &
      t764, t765, t766, t769, t772, t773, t776, t777, t78, t782, t785, t786, &
      t787, t79, t790, t792, t795, t796, t797, t8, t800, t801, t804, t815, &
      t822, t826, t827, t83, t833, t84, t85, t851, t852, t855, t857, t859, &
      t86, t862, t863, t864, t865, t866, t867, t87, t874, t877, t878
    REAL(kind=dp) :: t88, t881, t882, t886, t889, t893, t894, t897, t9, t905, &
      t906, t907, t914, t918, t92, t922, t923, t926, t929, t93, t930, t94, &
      t940, t941, t942, t944, t946, t949, t95, t950, t953, t954, t955, t958, &
      t959, t962, t963, t966, t969, t97, t972, t973, t974, t975, t976, t979, &
      t98, t980, t983, t984, t987, t99, t991, t995, t997, xx

    failure=.FALSE.
    my_epsilon_rho=0.5_dp*epsilon_rho
    epsilon_rho43 = my_epsilon_rho**(4.0_dp/3.0_dp)
    Cx = 1.5_dp*(3.0_dp/4.0_dp/pi)**(1.0_dp/3.0_dp)

    IF (.NOT. failure) THEN
      !$omp do
      DO ii=1,npoints
        rho=rho_spin(ii)
        ndrho=norm_drho_spin(ii)
        IF (rho>my_epsilon_rho) THEN
          IF( grad_deriv >=0 ) THEN
            t1 = rho ** (0.1e1_dp / 0.3e1_dp)
            xx = ndrho * MAX(1.0_dp/t1/rho,epsilon_rho43)
            t4 = lambda * rho
            t6 = ndrho ** 2
            t7 = beta * t6
            t8 = lambda ** 2
            t9 = t8 ** 2
            t10 = t9 ** 2
            t11 = 0.1e1_dp / t10
            t12 = rho ** 2
            t13 = t9 * t8
            t14 = 0.1e1_dp / t13
            t15 = t12 * t14
            t16 = t8 * lambda
            t17 = 0.1e1_dp / t16
            t18 = rho * t17
            t19 = t18 ** (0.1e1_dp / 0.3e1_dp)
            t20 = t19 ** 2
            t22 = 0.1e1_dp / t20 / t15
            t23 = t11 * t22
            t24 = beta * ndrho
            t25 = 0.1e1_dp / t9
            t26 = t19 * t18
            t27 = 0.1e1_dp / t26
            t31 = LOG(ndrho * t25 * t27 + SQRT((ndrho * t25 * t27) ** 0.2e1_dp +  &
            0.1e1_dp))
            t35 = 0.10e1_dp + 0.60e1_dp * t24 * t25 * t27 * t31
            t36 = 0.1e1_dp / t35
            t37 = t23 * t36
            t40 = 0.20e1_dp * Cx + 0.20e1_dp * t7 * t37
            t41 = SQRT(t40)
            t42 = t41 * t40
            t43 = t4 * t42
            t44 = rootpi
            t45 = 0.1e1_dp / t44
            t46 = omega * t45
            t47 = 0.1e1_dp / omega
            t48 = 0.1e1_dp / t41
            t49 = t47 * t48
            t50 = t44 * t19
            t53 = erf(0.3000000000e1_dp * t49 * t50)
            t55 = omega * t41
            t56 = 0.1e1_dp / t19
            t57 = t45 * t56
            t58 = omega ** 2
            t59 = 0.1e1_dp / t58
            t60 = 0.1e1_dp / t40
            t61 = t59 * t60
            t62 = pi * t20
            t63 = t61 * t62
            t65 = EXP(-0.8999999998e1_dp * t63)
            t66 = t58 * t40
            t67 = 0.1e1_dp / pi
            t68 = 0.1e1_dp / t20
            t69 = t67 * t68
            t70 = t65 - 0.10e1_dp
            t71 = t69 * t70
            t74 = t65 - 0.15e1_dp - 0.5555555558e-1_dp * t66 * t71
            t78 = t44 * t53 + 0.3333333334e0_dp * t55 * t57 * t74
            t79 = t46 * t78
            t83 = rho * t41 * omega
            t84 = t45 * t78
            t85 = t10 * lambda
            t86 = 0.1e1_dp / t85
            t87 = t86 * t22
            t88 = t87 * t36
            t92 = 0.1e1_dp / t10 / t9
            t93 = t7 * t92
            t94 = t12 * rho
            t95 = t94 * t86
            t97 = 0.1e1_dp / t20 / t95
            t98 = t97 * t36
            t99 = t98 * rho
            t102 = t7 * t11
            t103 = t35 ** 2
            t104 = 0.1e1_dp / t103
            t105 = t22 * t104
            t106 = t9 * lambda
            t107 = 0.1e1_dp / t106
            t108 = t107 * t27
            t112 = t24 * t11
            t114 = 0.1e1_dp / t19 / t15
            t115 = t114 * t31
            t116 = t115 * rho
            t119 = t24 * t25
            t122 = ndrho * t11
            t123 = t114 * rho
            t126 = -4 * ndrho * t107 * t27 + 4 * t122 * t123
            t130 = 1 + t6 * t11 * t22
            t131 = SQRT(t130)
            t132 = 0.1e1_dp / t131
            t133 = t27 * t126 * t132
            t136 = -0.240e2_dp * t24 * t108 * t31 + 0.240e2_dp * t112 * t116 + 0.60e1_dp  &
            * t119 * t133
            t137 = t105 * t136
            t140 = -0.160e2_dp * t7 * t88 + 0.160e2_dp * t93 * t99 - 0.20e1_dp * t102 *  &
            t137
            t141 = t84 * t140
            t144 = rho * t42
            t145 = SQRT(0.3141592654e1_dp)
            t146 = 0.1e1_dp / t145
            t147 = t44 * t146
            t149 = EXP(-0.9000000000e1_dp * t63)
            t150 = 0.1e1_dp / t42
            t151 = t47 * t150
            t155 = t49 * t44
            t156 = t68 * rho
            t160 = -0.1500000000e1_dp * t151 * t50 * t140 - 0.3000000000e1_dp * t155  &
            * t156 * t25
            t161 = t149 * t160
            t164 = omega * t48
            t165 = t164 * t45
            t166 = t56 * t74
            t170 = t55 * t45
            t171 = t27 * t74
            t172 = rho * t25
            t176 = t40 ** 2
            t177 = 0.1e1_dp / t176
            t178 = t59 * t177
            t182 = t61 * pi
            t183 = t56 * rho
            t187 = 0.8999999998e1_dp * t178 * t62 * t140 + 0.1800000000e2_dp * t182 *  &
            t183 * t25
            t189 = t58 * t140
            t192 = t66 * t67
            t194 = 0.1e1_dp / t20 / t18
            t195 = t194 * t70
            t196 = t195 * t172
            t199 = t68 * t187
            t200 = t199 * t65
            t203 = t187 * t65 - 0.5555555558e-1_dp * t189 * t71 - 0.1111111112e0_dp *  &
            t192 * t196 - 0.5555555558e-1_dp * t192 * t200
            t207 = (2 * t147 * t161) + 0.1666666667e0_dp * t165 * t166 * t140  &
            + 0.3333333334e0_dp * t170 * t171 * t172 + 0.3333333334e0_dp * t55 * t57  &
            * t203
            e_0(ii) = e_0(ii) +  (-0.4444444448e0_dp * t43 * t79 + t8 * (-0.3333333336e0_dp * t83 * t141  &
            - 0.2222222224e0_dp * t144 * t46 * t207)) * sx
          END IF
          IF( grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
            t214 = lambda * t42
            t217 = t4 * t41
            t218 = t10 * t16
            t219 = 0.1e1_dp / t218
            t220 = t219 * t97
            t221 = t220 * t36
            t224 = t9 * t16
            t225 = 0.1e1_dp / t224
            t230 = t220 * t132
            t233 = -0.8000000000e1_dp * t24 * t225 * t114 * t31 - 0.8000000000e1_dp *  &
            t7 * t230
            t234 = t105 * t233
            t237 = -0.5333333333e1_dp * t7 * t221 - 0.20e1_dp * t102 * t234
            t245 = t44 * t68
            t249 = -0.1500000000e1_dp * t151 * t50 * t237 + 0.1000000000e1_dp * t49 *  &
            t245 * t17
            t250 = t149 * t249
            t260 = t178 * t62 * t237
            t262 = pi * t56
            t264 = t61 * t262 * t17
            t266 = 0.8999999998e1_dp * t260 - 0.5999999999e1_dp * t264
            t267 = t266 * t65
            t268 = t58 * t237
            t271 = t195 * t17
            t274 = t68 * t266
            t275 = t274 * t65
            t278 = t267 - 0.5555555558e-1_dp * t268 * t71 + 0.3703703705e-1_dp * t192  &
            * t271 - 0.5555555558e-1_dp * t192 * t275
            t282 = (2 * t147 * t250) + 0.1666666667e0_dp * t165 * t166 * t237  &
            - 0.1111111111e0_dp * t170 * t171 * t17 + 0.3333333334e0_dp * t55 * t57  &
            * t278
            t283 = t46 * t282
            t289 = rho * t48 * omega
            t290 = t140 * t237
            t294 = t45 * t282
            t295 = t294 * t140
            t298 = t92 * t97
            t299 = t298 * t36
            t302 = t7 * t86
            t306 = 0.1e1_dp / t10 / t224
            t307 = t7 * t306
            t308 = t12 ** 2
            t309 = t308 * t92
            t311 = 0.1e1_dp / t20 / t309
            t313 = t311 * t36 * rho
            t316 = t97 * t104
            t317 = rho * t233
            t318 = t316 * t317
            t321 = t7 * t219
            t322 = t316 * t136
            t326 = 0.1e1_dp / t103 / t35
            t327 = t22 * t326
            t328 = t136 * t233
            t329 = t327 * t328
            t332 = t11 * t114
            t336 = t298 * t132
            t339 = t24 * t219
            t341 = 0.1e1_dp / t19 / t95
            t342 = t341 * t31
            t343 = t342 * rho
            t347 = t311 * t132 * rho
            t350 = t24 * t225
            t352 = t114 * t126 * t132
            t356 = ndrho * t219
            t357 = t341 * rho
            t360 = 0.28e2_dp / 0.3e1_dp * t122 * t114 - 0.28e2_dp / 0.3e1_dp * t356 * t357
            t362 = t27 * t360 * t132
            t365 = t6 * ndrho
            t366 = beta * t365
            t367 = t308 * rho
            t368 = 0.1e1_dp / t367
            t371 = 0.1e1_dp / t131 / t130
            t372 = t368 * t126 * t371
            t375 = 0.5600000000e2_dp * t24 * t332 * t31 + 0.3200000000e2_dp * t7 * t336  &
            - 0.5600000000e2_dp * t339 * t343 - 0.3200000000e2_dp * t307 * t347 -  &
            0.8000000000e1_dp * t350 * t352 + 0.60e1_dp * t119 * t362 + 0.8000000000e1_dp  &
            * t366 * t372
            t376 = t105 * t375
            t379 = 0.5866666667e2_dp * t7 * t299 + 0.160e2_dp * t302 * t234 - 0.5866666667e2_dp  &
            * t307 * t313 - 0.160e2_dp * t93 * t318 + 0.5333333333e1_dp * t321  &
            * t322 + 0.40e1_dp * t102 * t329 - 0.20e1_dp * t102 * t376
            t380 = t84 * t379
            t383 = t42 * omega
            t384 = t45 * t207
            t387 = t384 * t237
            t392 = 0.9000000000e1_dp * t260 - 0.6000000000e1_dp * t264
            t393 = t392 * t149
            t394 = t393 * t160
            t398 = 0.1e1_dp / t41 / t176
            t399 = t47 * t398
            t400 = t399 * t44
            t401 = t19 * t140
            t405 = t151 * t44
            t406 = t68 * t140
            t413 = t25 * t237
            t417 = t194 * rho
            t424 = 0.2250000000e1_dp * t400 * t401 * t237 - 0.5000000000e0_dp * t405  &
            * t406 * t17 - 0.1500000000e1_dp * t151 * t50 * t379 + 0.1500000000e1_dp  &
            * t405 * t156 * t413 + 0.2000000000e1_dp * t155 * t417 * t225 - 0.3000000000e1_dp  &
            * t49 * t245 * t25
            t425 = t149 * t424
            t428 = omega * t150
            t429 = t428 * t45
            t433 = t140 * t17
            t437 = t56 * t278
            t444 = t45 * t27
            t445 = t164 * t444
            t446 = t74 * rho
            t450 = t114 * t74
            t451 = rho * t225
            t455 = t27 * t278
            t462 = t56 * t203
            t466 = t27 * t203
            t470 = t176 * t40
            t471 = 0.1e1_dp / t470
            t472 = t59 * t471
            t473 = t472 * pi
            t474 = t20 * t140
            t478 = t178 * pi
            t479 = t56 * t140
            t489 = t27 * rho
            t496 = -0.1800000000e2_dp * t473 * t474 * t237 + 0.5999999999e1_dp * t478  &
            * t479 * t17 + 0.8999999998e1_dp * t178 * t62 * t379 - 0.1800000000e2_dp  &
            * t478 * t183 * t413 - 0.6000000000e1_dp * t182 * t489 * t225 + 0.1800000000e2_dp  &
            * t61 * t262 * t25
            t498 = t187 * t266
            t500 = t58 * t379
            t503 = t189 * t67
            t508 = t268 * t67
            t511 = t22 * t70
            t512 = t511 * t451
            t515 = t67 * t194
            t516 = t66 * t515
            t517 = t267 * t172
            t520 = t195 * t25
            t526 = t65 * t17
            t527 = t194 * t187 * t526
            t530 = t68 * t496
            t531 = t530 * t65
            t534 = t199 * t267
            t537 = t496 * t65 + t498 * t65 - 0.5555555558e-1_dp * t500 * t71 + 0.3703703705e-1_dp  &
            * t503 * t271 - 0.5555555558e-1_dp * t503 * t275 - 0.1111111112e0_dp  &
            * t508 * t196 + 0.1851851853e0_dp * t192 * t512 - 0.1111111112e0_dp  &
            * t516 * t517 - 0.1111111112e0_dp * t192 * t520 - 0.5555555558e-1_dp  &
            * t508 * t200 + 0.3703703705e-1_dp * t192 * t527 - 0.5555555558e-1_dp  &
            * t192 * t531 - 0.5555555558e-1_dp * t192 * t534
            t541 = (2 * t147 * t394) + (2 * t147 * t425) - 0.8333333335e-1_dp  &
            * t429 * t166 * t290 - 0.5555555557e-1_dp * t165 * t171 * t433 +  &
            0.1666666667e0_dp * t165 * t437 * t140 + 0.1666666667e0_dp * t165 * t166  &
            * t379 + 0.1666666667e0_dp * t445 * t446 * t413 - 0.4444444445e0_dp * t170  &
            * t450 * t451 + 0.3333333334e0_dp * t170 * t455 * t172 + 0.3333333334e0_dp  &
            * t170 * t171 * t25 + 0.1666666667e0_dp * t165 * t462 * t237 -  &
            0.1111111111e0_dp * t170 * t466 * t17 + 0.3333333334e0_dp * t55 * t57 *  &
            t537
            e_rho_spin(ii) = e_rho_spin(ii) +  (-0.4444444448e0_dp * t214 * t79 - 0.6666666672e0_dp * t217 * t46  &
            * t78 * t237 - 0.4444444448e0_dp * t43 * t283 + t8 * (-0.3333333336e0_dp  &
            * t55 * t141 - 0.1666666668e0_dp * t289 * t84 * t290 - 0.3333333336e0_dp  &
            * t83 * t295 - 0.3333333336e0_dp * t83 * t380 - 0.2222222224e0_dp * t383  &
            * t384 - 0.3333333336e0_dp * t83 * t387 - 0.2222222224e0_dp * t144 *  &
            t46 * t541)) * sx
            t550 = beta * t25
            t551 = t27 * t31
            t557 = 0.60e1_dp * t550 * t551 + 0.60e1_dp * t24 * t23 * t132
            t558 = t105 * t557
            t561 = 0.40e1_dp * t24 * t37 - 0.20e1_dp * t102 * t558
            t566 = pi * t146
            t567 = t566 * t149
            t569 = t151 * t19 * t561
            t575 = t20 * t561
            t579 = t58 * t561
            t585 = 0.8999999998e1_dp * t478 * t575 * t65 - 0.5555555558e-1_dp * t579  &
            * t71 - 0.5000000001e0_dp * t60 * t561 * t65
            t589 = -0.3000000000e1_dp * t567 * t569 + 0.1666666667e0_dp * t165 * t166  &
            * t561 + 0.3333333334e0_dp * t55 * t57 * t585
            t590 = t46 * t589
            t593 = t140 * t561
            t597 = t45 * t589
            t598 = t597 * t140
            t605 = t24 * t92
            t608 = rho * t557
            t609 = t316 * t608
            t614 = t136 * t557
            t615 = t327 * t614
            t624 = beta * t11
            t627 = t97 * t132
            t628 = t627 * rho
            t635 = -4 * t108 + 4 * t332 * rho
            t637 = t27 * t635 * t132
            t640 = 0.1e1_dp / t308
            t645 = -0.240e2_dp * beta * t107 * t551 - 0.240e2_dp * t24 * t87 * t132 +  &
            0.240e2_dp * t624 * t116 + 0.240e2_dp * t605 * t628 + 0.60e1_dp * t550 * t133  &
            + 0.60e1_dp * t119 * t637 - 0.60e1_dp * t7 * t640 * t126 * t371
            t646 = t105 * t645
            t649 = -0.320e2_dp * t24 * t88 + 0.160e2_dp * t302 * t558 + 0.320e2_dp * t605  &
            * t99 - 0.160e2_dp * t93 * t609 - 0.40e1_dp * t112 * t137 + 0.40e1_dp * t102  &
            * t615 - 0.20e1_dp * t102 * t646
            t650 = t84 * t649
            t653 = t384 * t561
            t657 = t44 * pi * t146
            t658 = t657 * t178
            t668 = t25 * t561
            t672 = 0.2250000000e1_dp * t400 * t401 * t561 - 0.1500000000e1_dp * t151  &
            * t50 * t649 + 0.1500000000e1_dp * t405 * t156 * t668
            t673 = t149 * t672
            t679 = t56 * t585
            t689 = t27 * t585
            t705 = -0.1800000000e2_dp * t473 * t474 * t561 + 0.8999999998e1_dp * t178  &
            * t62 * t649 - 0.1800000000e2_dp * t478 * t183 * t668
            t707 = t187 * t59
            t708 = t707 * t177
            t709 = t561 * t65
            t710 = t62 * t709
            t713 = t58 * t649
            t716 = t140 * t177
            t719 = t579 * t67
            t722 = 0.1e1_dp / lambda
            t723 = t60 * t722
            t728 = t68 * t705
            t729 = t728 * t65
            t732 = t60 * t187
            t735 = t705 * t65 + 0.8999999998e1_dp * t708 * t710 - 0.5555555558e-1_dp  &
            * t713 * t71 - 0.5000000001e0_dp * t716 * t709 - 0.1111111112e0_dp * t719  &
            * t196 - 0.1000000001e1_dp * t723 * t709 - 0.5555555558e-1_dp * t719 *  &
            t200 - 0.5555555558e-1_dp * t192 * t729 - 0.5000000001e0_dp * t732 * t709
            t739 = 0.1800000000e2_dp * t658 * t575 * t161 + (2 * t147 * t673)  &
            - 0.8333333335e-1_dp * t429 * t166 * t593 + 0.1666666667e0_dp * t165 * t679  &
            * t140 + 0.1666666667e0_dp * t165 * t166 * t649 + 0.1666666667e0_dp  &
            * t445 * t446 * t668 + 0.3333333334e0_dp * t170 * t689 * t172 + 0.1666666667e0_dp  &
            * t165 * t462 * t561 + 0.3333333334e0_dp * t55 * t57 * t735
            e_ndrho_spin(ii) = e_ndrho_spin(ii) +  (-0.6666666672e0_dp * t217 * t46 * t78 * t561 - 0.4444444448e0_dp  &
            * t43 * t590 + t8 * (-0.1666666668e0_dp * t289 * t84 * t593 - 0.3333333336e0_dp  &
            * t83 * t598 - 0.3333333336e0_dp * t83 * t650 - 0.3333333336e0_dp  &
            * t83 * t653 - 0.2222222224e0_dp * t144 * t46 * t739)) * sx
          END IF

          IF( grad_deriv >= 2 .OR. grad_deriv == -2 ) THEN
            t747 = lambda * t41 * omega
            t753 = t4 * t48
            t754 = t237 ** 2
            t764 = 0.1e1_dp / t10 / t13
            t765 = t764 * t311
            t766 = t765 * t36
            t769 = t316 * t233
            t772 = t233 ** 2
            t773 = t327 * t772
            t776 = t10 * t8
            t777 = 0.1e1_dp / t776
            t782 = t765 * t132
            t785 = t6 ** 2
            t786 = beta * t785
            t787 = t10 ** 2
            t790 = t308 * t94
            t792 = 0.1e1_dp / t787 / t106
            t795 = 0.1e1_dp / t19 / t790 / t792
            t796 = 0.1e1_dp / t787 / t13 * t795
            t797 = t796 * t371
            t800 = 0.1866666667e2_dp * t24 * t777 * t341 * t31 + 0.4000000000e2_dp *  &
            t7 * t782 - 0.1066666667e2_dp * t786 * t797
            t801 = t105 * t800
            t804 = 0.1955555555e2_dp * t7 * t766 + 0.1066666667e2_dp * t321 * t769 +  &
            0.40e1_dp * t102 * t773 - 0.20e1_dp * t102 * t801
            t815 = t68 * t237
            t822 = t44 * t194
            t826 = 0.2250000000e1_dp * t399 * t50 * t754 - 0.1000000000e1_dp * t405 *  &
            t815 * t17 - 0.1500000000e1_dp * t151 * t50 * t804 - 0.6666666667e0_dp  &
            * t49 * t822 * t14
            t827 = t149 * t826
            t833 = t237 * t17
            t851 = 0.1800000000e2_dp * t472 * t62 * t754
            t852 = t56 * t237
            t855 = 0.1200000000e2_dp * t478 * t852 * t17
            t857 = t178 * t62 * t804
            t859 = pi * t27
            t862 = 0.2000000000e1_dp * t61 * t859 * t14
            t863 = -t851 + t855 + 0.8999999998e1_dp * t857 + t862
            t864 = t863 * t65
            t865 = t266 ** 2
            t866 = t865 * t65
            t867 = t58 * t804
            t874 = t511 * t14
            t877 = t194 * t266
            t878 = t877 * t526
            t881 = t68 * t863
            t882 = t881 * t65
            t886 = t68 * t865 * t65
            t889 = t864 + t866 - 0.5555555558e-1_dp * t867 * t71 + 0.7407407410e-1_dp  &
            * t508 * t271 - 0.1111111112e0_dp * t508 * t275 - 0.6172839508e-1_dp *  &
            t192 * t874 + 0.7407407410e-1_dp * t192 * t878 - 0.5555555558e-1_dp * t192  &
            * t882 - 0.5555555558e-1_dp * t192 * t886
            t893 = (2 * t147 * t393 * t249) + (2 * t147 * t827) - 0.8333333335e-1_dp  &
            * t429 * t166 * t754 - 0.1111111111e0_dp * t165 * t171 * t833  &
            + 0.3333333334e0_dp * t165 * t437 * t237 + 0.1666666667e0_dp * t165  &
            * t166 * t804 + 0.1481481481e0_dp * t170 * t450 * t14 - 0.2222222222e0_dp  &
            * t170 * t455 * t17 + 0.3333333334e0_dp * t55 * t57 * t889
            t894 = t46 * t893
            t897 = t78 * t140
            t905 = rho * t150
            t906 = t905 * omega
            t907 = t140 * t754
            t914 = t379 * t237
            t918 = t140 * t804
            t922 = t45 * t893
            t923 = t922 * t140
            t926 = t294 * t379
            t929 = t306 * t311
            t930 = t929 * t36
            t940 = 0.1e1_dp / t787 / t8
            t941 = t7 * t940
            t942 = t367 * t306
            t944 = 0.1e1_dp / t20 / t942
            t946 = t944 * t36 * rho
            t949 = t311 * t104
            t950 = t949 * t317
            t953 = t97 * t326
            t954 = rho * t772
            t955 = t953 * t954
            t958 = rho * t800
            t959 = t316 * t958
            t962 = t7 * t764
            t963 = t949 * t136
            t966 = t953 * t328
            t969 = t316 * t375
            t972 = t103 ** 2
            t973 = 0.1e1_dp / t972
            t974 = t22 * t973
            t975 = t136 * t772
            t976 = t974 * t975
            t979 = t375 * t233
            t980 = t327 * t979
            t983 = t136 * t800
            t984 = t327 * t983
            t987 = t219 * t341
            t991 = t929 * t132
            t995 = 0.1e1_dp / t787 / t224
            t997 = t995 * t795 * t371
            t1000 = t24 * t764
            t1002 = 0.1e1_dp / t19 / t309
            t1004 = t1002 * t31 * rho
            t1008 = t944 * t132 * rho
            t1012 = 0.1e1_dp / t787 / t776
            t1014 = t308 ** 2
            t1019 = 0.1e1_dp / t19 / t1014 * t787 * t10
            t1021 = t1019 * t371 * rho
            t1024 = t24 * t777
            t1026 = t341 * t126 * t132
            t1030 = t114 * t360 * t132
            t1033 = t308 * t12
            t1034 = 0.1e1_dp / t1033
            t1036 = t1034 * t126 * t371
            t1040 = ndrho * t764
            t1044 = -0.280e3_dp / 0.9e1_dp * t356 * t341 + 0.280e3_dp / 0.9e1_dp * t1040 *  &
            t1002 * rho
            t1046 = t27 * t1044 * t132
            t1050 = t368 * t360 * t371
            t1054 = beta * t785 * ndrho
            t1055 = t1054 * t368
            t1056 = t130 ** 2
            t1058 = 0.1e1_dp / t131 / t1056
            t1059 = t126 * t1058
            t1060 = t1059 * t220
            t1063 = -0.1866666667e3_dp * t24 * t987 * t31 - 0.2240000000e3_dp * t7 *  &
            t991 + 0.4266666667e2_dp * t786 * t997 + 0.1866666667e3_dp * t1000 * t1004  &
            + 0.2240000000e3_dp * t941 * t1008 - 0.4266666667e2_dp * t786 * t1012  &
            * t1021 + 0.1866666667e2_dp * t1024 * t1026 - 0.1600000000e2_dp * t350  &
            * t1030 - 0.5066666667e2_dp * t366 * t1036 + 0.60e1_dp * t119 * t1046 +  &
            0.1600000000e2_dp * t366 * t1050 + 0.3200000000e2_dp * t1055 * t1060
            t1064 = t105 * t1063
            t1067 = -0.2737777778e3_dp * t7 * t930 - 0.1173333333e3_dp * t93 * t769 -  &
            0.320e2_dp * t302 * t773 + 0.160e2_dp * t302 * t801 + 0.2737777778e3_dp *  &
            t941 * t946 + 0.1173333333e3_dp * t307 * t950 + 0.320e2_dp * t93 * t955  &
            - 0.160e2_dp * t93 * t959 - 0.1955555555e2_dp * t962 * t963 - 0.2133333334e2_dp  &
            * t321 * t966 + 0.1066666667e2_dp * t321 * t969 - 0.120e2_dp * t102  &
            * t976 + 0.80e1_dp * t102 * t980 + 0.40e1_dp * t102 * t984 - 0.20e1_dp * t102  &
            * t1064
            t1068 = t84 * t1067
            t1071 = t45 * t541
            t1076 = t384 * t754
            t1079 = t1071 * t237
            t1082 = t384 * t804
            t1085 = t379 * t17
            t1089 = t114 * t278
            t1093 = t25 * t804
            t1097 = t187 * t865
            t1099 = t97 * t70
            t1100 = rho * t777
            t1101 = t1099 * t1100
            t1106 = t500 * t67
            t1109 = t187 * t863
            t1111 = t496 * t266
            t1119 = t194 * t496 * t526
            t1122 = t199 * t866
            t1125 = t65 * t25
            t1126 = t877 * t1125
            t1129 = t867 * t67
            t1134 = t866 * t172
            t1137 = t268 * t515
            t1140 = t176 ** 2
            t1141 = 0.1e1_dp / t1140
            t1142 = t59 * t1141
            t1143 = t1142 * pi
            t1150 = t20 * t379
            t1157 = t27 * t140
            t1161 = t56 * t379
            t1168 = t25 * t754
            t1172 = t225 * t237
            t1176 = t56 * t25
            t1189 = 0.5400000000e2_dp * t1143 * t474 * t754 - 0.2400000000e2_dp * t473  &
            * t479 * t833 - 0.3600000000e2_dp * t473 * t1150 * t237 - 0.1800000000e2_dp  &
            * t473 * t474 * t804 - 0.2000000000e1_dp * t478 * t1157 * t14 +  &
            0.1200000000e2_dp * t478 * t1161 * t17 + 0.8999999998e1_dp * t178 * t62  &
            * t1067 + 0.3600000000e2_dp * t473 * t183 * t1168 + 0.1200000000e2_dp *  &
            t478 * t489 * t1172 - 0.3600000000e2_dp * t478 * t1176 * t237 - 0.1800000000e2_dp  &
            * t478 * t183 * t1093 + 0.8000000000e1_dp * t182 * t123 *  &
            t777 - 0.1200000000e2_dp * t61 * t859 * t225
            t1191 = (t1097 * t65) - 0.4938271608e0_dp * t192 * t1101 - 0.6172839508e-1_dp  &
            * t503 * t874 - 0.1111111112e0_dp * t1106 * t275 + (t1109  &
            * t65) + (2 * t1111 * t65) + 0.3703703706e0_dp * t508 * t512 +  &
            0.7407407410e-1_dp * t503 * t878 + 0.7407407410e-1_dp * t192 * t1119 - 0.5555555558e-1_dp  &
            * t192 * t1122 - 0.2222222224e0_dp * t192 * t1126 - 0.1111111112e0_dp  &
            * t1129 * t196 - 0.5555555558e-1_dp * t503 * t882 - 0.1111111112e0_dp  &
            * t516 * t1134 - 0.2222222224e0_dp * t1137 * t517 + (t1189  &
            * t65)
            t1192 = t67 * t22
            t1193 = t66 * t1192
            t1194 = t267 * t451
            t1200 = t65 * t14
            t1201 = t22 * t187 * t1200
            t1206 = t68 * t1189
            t1207 = t1206 * t65
            t1210 = t58 * t1067
            t1213 = t511 * t225
            t1216 = t864 * t172
            t1221 = t498 * t526
            t1224 = t530 * t267
            t1235 = t199 * t864
            t1238 = 0.3703703706e0_dp * t1193 * t1194 - 0.5555555558e-1_dp * t503 * t886  &
            - 0.6172839508e-1_dp * t192 * t1201 - 0.5555555558e-1_dp * t1129 * t200  &
            - 0.5555555558e-1_dp * t192 * t1207 - 0.5555555558e-1_dp * t1210 * t71  &
            + 0.3703703706e0_dp * t192 * t1213 - 0.1111111112e0_dp * t516 * t1216  &
            - 0.1111111112e0_dp * t508 * t534 + 0.7407407410e-1_dp * t516 * t1221 -  &
            0.1111111112e0_dp * t192 * t1224 + 0.7407407410e-1_dp * t508 * t527 + 0.7407407410e-1_dp  &
            * t1106 * t271 - 0.1111111112e0_dp * t508 * t531 - 0.2222222224e0_dp  &
            * t508 * t520 - 0.5555555558e-1_dp * t192 * t1235
            t1239 = t1191 + t1238
            t1246 = t428 * t444
            t1247 = t74 * t140
            t1254 = t27 * t889
            t1264 = t27 * t537
            t1268 = t341 * t74
            t1275 = t56 * t889
            t1283 = -t851 + t855 + 0.9000000000e1_dp * t857 + t862
            t1284 = t1283 * t149
            t1288 = -0.1111111111e0_dp * t165 * t171 * t1085 - 0.8888888890e0_dp * t170  &
            * t1089 * t451 + 0.1666666667e0_dp * t445 * t446 * t1093 + 0.3333333334e0_dp  &
            * t55 * t57 * t1239 - 0.8333333335e-1_dp * t429 * t462 * t754  &
            + 0.5555555556e-1_dp * t1246 * t1247 * t833 + (4 * t147 * t393 *  &
            t424) + 0.3333333334e0_dp * t170 * t1254 * t172 - 0.1111111111e0_dp * t165  &
            * t466 * t833 + 0.3333333334e0_dp * t165 * t437 * t379 - 0.2222222222e0_dp  &
            * t170 * t1264 * t17 + 0.1037037037e1_dp * t170 * t1268 * t1100  &
            + 0.1666666667e0_dp * t165 * t166 * t1067 + 0.1666666667e0_dp * t165 *  &
            t1275 * t140 - 0.8333333335e-1_dp * t1246 * t446 * t1168 + (2 *  &
            t147 * t1284 * t160)
            t1292 = t278 * rho
            t1300 = 0.1e1_dp / t41 / t470
            t1301 = t47 * t1300
            t1302 = t1301 * t44
            t1309 = t19 * t379
            t1316 = t194 * t140
            t1320 = t68 * t379
            t1333 = t68 * t25
            t1340 = t22 * rho
            t1347 = -0.5625000000e1_dp * t1302 * t401 * t754 + 0.1500000000e1_dp * t400  &
            * t406 * t833 + 0.4500000000e1_dp * t400 * t1309 * t237 + 0.2250000000e1_dp  &
            * t400 * t401 * t804 + 0.3333333333e0_dp * t405 * t1316 * t14  &
            - 0.1000000000e1_dp * t405 * t1320 * t17 - 0.1500000000e1_dp * t151 * t50  &
            * t1067 - 0.2250000000e1_dp * t400 * t156 * t1168 - 0.2000000000e1_dp  &
            * t405 * t417 * t1172 + 0.3000000000e1_dp * t405 * t1333 * t237 + 0.1500000000e1_dp  &
            * t405 * t156 * t1093 - 0.3333333333e1_dp * t155 * t1340  &
            * t777 + 0.4000000000e1_dp * t49 * t822 * t225
            t1348 = t149 * t1347
            t1360 = t140 * t14
            t1367 = omega * t398
            t1368 = t1367 * t45
            t1372 = t45 * t114
            t1373 = t164 * t1372
            t1377 = t114 * t203
            t1387 = t392 ** 2
            t1388 = t1387 * t149
            t1392 = t56 * t537
            t1396 = 0.6666666668e0_dp * t170 * t455 * t25 + 0.3333333334e0_dp * t445  &
            * t1292 * t413 + 0.1666666667e0_dp * t165 * t462 * t804 + (2 * t147  &
            * t1348) - 0.1111111111e0_dp * t165 * t455 * t433 - 0.8333333335e-1_dp  &
            * t429 * t166 * t918 + 0.3333333334e0_dp * t165 * t171 * t413 + 0.7407407409e-1_dp  &
            * t165 * t450 * t1360 - 0.1666666667e0_dp * t429 * t437  &
            * t290 + 0.1250000000e0_dp * t1368 * t166 * t907 - 0.4444444445e0_dp * t1373  &
            * t446 * t1172 + 0.1481481481e0_dp * t170 * t1377 * t14 - 0.8888888890e0_dp  &
            * t170 * t450 * t225 - 0.1666666667e0_dp * t429 * t166 * t914  &
            + (2 * t147 * t1388 * t160) + 0.3333333334e0_dp * t165 * t1392  &
            * t237
            t1397 = t1288 + t1396
            t1401 = -0.3333333336e0_dp * t165 * t897 * t237 - 0.6666666672e0_dp * t55  &
            * t295 - 0.6666666672e0_dp * t55 * t380 + 0.8333333340e-1_dp * t906 * t84  &
            * t907 - 0.3333333336e0_dp * t289 * t294 * t290 - 0.3333333336e0_dp *  &
            t289 * t84 * t914 - 0.1666666668e0_dp * t289 * t84 * t918 - 0.3333333336e0_dp  &
            * t83 * t923 - 0.6666666672e0_dp * t83 * t926 - 0.3333333336e0_dp  &
            * t83 * t1068 - 0.4444444448e0_dp * t383 * t1071 - 0.6666666672e0_dp *  &
            t55 * t387 - 0.1666666668e0_dp * t289 * t1076 - 0.6666666672e0_dp * t83  &
            * t1079 - 0.3333333336e0_dp * t83 * t1082 - 0.2222222224e0_dp * t144 * t46  &
            * t1397
            e_rho_rho_spin(ii) = e_rho_rho_spin(ii) +  (-0.1333333334e1_dp * t747 * t84 * t237 - 0.8888888896e0_dp * t214  &
            * t283 - 0.3333333336e0_dp * t753 * t46 * t78 * t754 - 0.1333333334e1_dp  &
            * t217 * t46 * t282 * t237 - 0.6666666672e0_dp * t217 * t46 * t78 *  &
            t804 - 0.4444444448e0_dp * t43 * t894 + t8 * t1401) * sx
            t1406 = 0.6666666672e0_dp * t747 * t84 * t561
            t1408 = 0.4444444448e0_dp * t214 * t590
            t1409 = t4 * t164
            t1410 = t237 * t561
            t1411 = t84 * t1410
            t1413 = 0.3333333336e0_dp * t1409 * t1411
            t1417 = 0.6666666672e0_dp * t217 * t46 * t589 * t237
            t1420 = t316 * t557
            t1425 = t233 * t557
            t1426 = t327 * t1425
            t1429 = beta * t225
            t1436 = t1033 * t940
            t1438 = 0.1e1_dp / t19 / t1436
            t1439 = 0.1e1_dp / t787 / t16 * t1438
            t1440 = t1439 * t371
            t1443 = -0.8000000000e1_dp * t1429 * t115 - 0.2400000000e2_dp * t24 * t230  &
            + 0.8000000000e1_dp * t366 * t1440
            t1444 = t105 * t1443
            t1447 = -0.1066666667e2_dp * t24 * t221 + 0.5333333333e1_dp * t321 * t1420  &
            - 0.40e1_dp * t112 * t234 + 0.40e1_dp * t102 * t1426 - 0.20e1_dp * t102 *  &
            t1444
            t1451 = 0.6666666672e0_dp * t217 * t46 * t78 * t1447
            t1455 = 0.6666666672e0_dp * t217 * t46 * t282 * t561
            t1459 = t19 * t237
            t1466 = t68 * t17
            t1470 = 0.2250000000e1_dp * t400 * t1459 * t561 - 0.1500000000e1_dp * t151  &
            * t50 * t1447 - 0.5000000000e0_dp * t405 * t1466 * t561
            t1471 = t149 * t1470
            t1476 = 0.8333333335e-1_dp * t429 * t166 * t1410
            t1479 = 0.1666666667e0_dp * t165 * t679 * t237
            t1482 = 0.1666666667e0_dp * t165 * t166 * t1447
            t1483 = t17 * t561
            t1485 = t165 * t171 * t1483
            t1489 = 0.1111111111e0_dp * t170 * t689 * t17
            t1492 = 0.1666666667e0_dp * t165 * t437 * t561
            t1493 = t20 * t237
            t1496 = 0.1800000000e2_dp * t473 * t1493 * t561
            t1498 = t178 * t62 * t1447
            t1500 = t56 * t17
            t1502 = t478 * t1500 * t561
            t1504 = -t1496 + 0.8999999998e1_dp * t1498 + 0.5999999999e1_dp * t1502
            t1505 = t1504 * t65
            t1506 = t266 * t59
            t1507 = t1506 * t177
            t1509 = 0.8999999998e1_dp * t1507 * t710
            t1510 = t58 * t1447
            t1512 = 0.5555555558e-1_dp * t1510 * t71
            t1513 = t237 * t177
            t1515 = 0.5000000001e0_dp * t1513 * t709
            t1517 = 0.3703703705e-1_dp * t719 * t271
            t1518 = 0.1e1_dp / rho
            t1519 = t60 * t1518
            t1523 = 0.5555555558e-1_dp * t719 * t275
            t1525 = t68 * t1504 * t65
            t1528 = t60 * t266
            t1530 = 0.5000000001e0_dp * t1528 * t709
            t1531 = t1505 + t1509 - t1512 - t1515 + t1517 + 0.3333333334e0_dp * t1519  &
            * t709 - t1523 - 0.5555555558e-1_dp * t192 * t1525 - t1530
            t1535 = 0.1800000000e2_dp * t658 * t575 * t250 + (2 * t147 * t1471)  &
            - t1476 + t1479 + t1482 - 0.5555555555e-1_dp * t1485 - t1489 + t1492  &
            + 0.3333333334e0_dp * t55 * t57 * t1531
            t1536 = t46 * t1535
            t1541 = 0.1666666668e0_dp * t165 * t897 * t561
            t1543 = 0.3333333336e0_dp * t55 * t598
            t1545 = 0.3333333336e0_dp * t55 * t650
            t1546 = t905 * t46
            t1547 = t897 * t1410
            t1549 = 0.8333333340e-1_dp * t1546 * t1547
            t1552 = 0.1666666668e0_dp * t289 * t597 * t290
            t1553 = t649 * t237
            t1556 = 0.1666666668e0_dp * t289 * t84 * t1553
            t1557 = t140 * t1447
            t1560 = 0.1666666668e0_dp * t289 * t84 * t1557
            t1563 = 0.1666666668e0_dp * t289 * t294 * t593
            t1564 = t45 * t1535
            t1565 = t1564 * t140
            t1568 = t294 * t649
            t1570 = 0.3333333336e0_dp * t83 * t1568
            t1571 = t379 * t561
            t1574 = 0.1666666668e0_dp * t289 * t84 * t1571
            t1575 = t597 * t379
            t1577 = 0.3333333336e0_dp * t83 * t1575
            t1579 = 0.1173333333e3_dp * t24 * t299
            t1581 = 0.5866666667e2_dp * t93 * t1420
            t1582 = t24 * t86
            t1584 = 0.320e2_dp * t1582 * t234
            t1586 = 0.320e2_dp * t302 * t1426
            t1588 = 0.160e2_dp * t302 * t1444
            t1589 = t24 * t306
            t1591 = 0.1173333333e3_dp * t1589 * t313
            t1592 = t949 * t608
            t1594 = 0.5866666667e2_dp * t307 * t1592
            t1596 = 0.320e2_dp * t605 * t318
            t1597 = t7 * t298
            t1598 = t326 * rho
            t1599 = t1598 * t1425
            t1601 = 0.320e2_dp * t1597 * t1599
            t1602 = rho * t1443
            t1603 = t316 * t1602
            t1605 = 0.160e2_dp * t93 * t1603
            t1607 = 0.1066666667e2_dp * t339 * t322
            t1608 = t953 * t614
            t1610 = 0.1066666667e2_dp * t321 * t1608
            t1611 = t316 * t645
            t1613 = 0.5333333333e1_dp * t321 * t1611
            t1615 = 0.80e1_dp * t112 * t329
            t1616 = t7 * t23
            t1617 = t973 * t136
            t1618 = t1617 * t1425
            t1620 = 0.120e2_dp * t1616 * t1618
            t1621 = t645 * t233
            t1622 = t327 * t1621
            t1624 = 0.40e1_dp * t102 * t1622
            t1625 = t136 * t1443
            t1626 = t327 * t1625
            t1628 = 0.40e1_dp * t102 * t1626
            t1630 = 0.40e1_dp * t112 * t376
            t1631 = t375 * t557
            t1632 = t327 * t1631
            t1634 = 0.40e1_dp * t102 * t1632
            t1636 = 0.5600000000e2_dp * t624 * t115
            t1638 = 0.1200000000e3_dp * t24 * t336
            t1640 = 0.1e1_dp / t787 / t9
            t1642 = t1640 * t1438 * t371
            t1644 = 0.3200000000e2_dp * t366 * t1642
            t1645 = beta * t219
            t1647 = 0.5600000000e2_dp * t1645 * t343
            t1649 = 0.1200000000e3_dp * t1589 * t347
            t1652 = t795 * t371 * rho
            t1654 = 0.3200000000e2_dp * t366 * t995 * t1652
            t1656 = 0.8000000000e1_dp * t1429 * t352
            t1658 = t114 * t635 * t132
            t1660 = 0.8000000000e1_dp * t350 * t1658
            t1662 = 0.3200000000e2_dp * t7 * t372
            t1664 = 0.60e1_dp * t550 * t362
            t1667 = 0.28e2_dp / 0.3e1_dp * t332 - 0.28e2_dp / 0.3e1_dp * t987 * rho
            t1669 = t27 * t1667 * t132
            t1671 = 0.60e1_dp * t119 * t1669
            t1675 = 0.60e1_dp * t7 * t640 * t360 * t371
            t1677 = t368 * t635 * t371
            t1679 = 0.8000000000e1_dp * t366 * t1677
            t1680 = t786 * t368
            t1681 = t1059 * t23
            t1684 = t1636 + t1638 - t1644 - t1647 - t1649 + t1654 - t1656 - t1660  &
            + t1662 + t1664 + t1671 - t1675 + t1679 - 0.2400000000e2_dp * t1680  &
            * t1681
            t1685 = t105 * t1684
            t1688 = t1579 - t1581 + t1584 - t1586 + t1588 - t1591 + t1594 - t1596  &
            + t1601 - t1605 + t1607 - t1610 + t1613 + t1615 - t1620 + t1624  &
            + t1628 - t1630 + t1634 - 0.20e1_dp * t102 * t1685
            t1689 = t84 * t1688
            t1693 = 0.3333333336e0_dp * t55 * t653
            t1694 = t45 * t739
            t1696 = 0.2222222224e0_dp * t383 * t1694
            t1699 = 0.1666666668e0_dp * t289 * t384 * t1410
            t1700 = t1694 * t237
            t1702 = 0.3333333336e0_dp * t83 * t1700
            t1703 = t384 * t1447
            t1705 = 0.3333333336e0_dp * t83 * t1703
            t1706 = t1071 * t561
            t1708 = 0.3333333336e0_dp * t83 * t1706
            t1709 = t187 * t1504
            t1712 = t266 * t561 * t65
            t1714 = 0.5000000001e0_dp * t732 * t1712
            t1715 = t713 * t67
            t1717 = 0.3703703705e-1_dp * t1715 * t271
            t1718 = t705 * t266
            t1719 = t1718 * t65
            t1720 = t140 * t266
            t1722 = t177 * t561 * t65
            t1724 = 0.5000000001e0_dp * t1720 * t1722
            t1726 = t194 * t705 * t526
            t1728 = 0.3703703705e-1_dp * t192 * t1726
            t1730 = 0.1111111112e0_dp * t719 * t520
            t1731 = t58 * t1688
            t1735 = 0.3703703705e-1_dp * t719 * t527
            t1737 = t722 * t561 * t65
            t1740 = t496 * t59
            t1741 = t1740 * t177
            t1743 = 0.8999999998e1_dp * t1741 * t710
            t1745 = t187 * t561 * t65
            t1748 = t140 * t1518
            t1752 = 0.1851851853e0_dp * t719 * t512
            t1754 = 0.5555555558e-1_dp * t1715 * t275
            t1755 = t1510 * t67
            t1757 = 0.1111111112e0_dp * t1755 * t196
            t1758 = t1709 * t65 - t1714 + t1717 + t1719 - t1724 + t1728 - t1730  &
            - 0.5555555558e-1_dp * t1731 * t71 + t1735 + 0.666666666e0_dp * t1519 *  &
            t1737 + t1743 + 0.3333333334e0_dp * t1519 * t1745 + 0.3333333334e0_dp *  &
            t1748 * t1722 + t1752 - t1754 - t1757
            t1760 = 0.5555555558e-1_dp * t719 * t534
            t1762 = 0.5555555558e-1_dp * t1755 * t200
            t1763 = t237 * t722
            t1765 = 0.1000000001e1_dp * t1763 * t1722
            t1766 = t237 * t187
            t1768 = 0.5000000001e0_dp * t1766 * t1722
            t1769 = t199 * t1505
            t1773 = 0.5555555558e-1_dp * t719 * t531
            t1774 = t60 * t496
            t1776 = 0.5000000001e0_dp * t1774 * t709
            t1777 = t379 * t177
            t1779 = 0.5000000001e0_dp * t1777 * t709
            t1782 = t728 * t267
            t1784 = 0.5555555558e-1_dp * t192 * t1782
            t1785 = t579 * t515
            t1787 = 0.1111111112e0_dp * t1785 * t517
            t1789 = 0.5555555558e-1_dp * t508 * t729
            t1790 = t1505 * t172
            t1794 = 0.1000000001e1_dp * t723 * t1712
            t1795 = t498 * t178
            t1797 = 0.8999999998e1_dp * t1795 * t710
            t1800 = 0.5400000000e2_dp * t1143 * t474 * t1410
            t1801 = t20 * t649
            t1804 = 0.1800000000e2_dp * t473 * t1801 * t237
            t1807 = 0.1800000000e2_dp * t473 * t474 * t1447
            t1810 = 0.1200000000e2_dp * t473 * t479 * t1483
            t1811 = t56 * t649
            t1814 = 0.5999999999e1_dp * t478 * t1811 * t17
            t1817 = 0.1800000000e2_dp * t473 * t1150 * t561
            t1821 = t472 * t262
            t1822 = t172 * t1410
            t1824 = 0.3600000000e2_dp * t1821 * t1822
            t1825 = t25 * t1447
            t1828 = 0.1800000000e2_dp * t478 * t183 * t1825
            t1829 = t225 * t561
            t1832 = 0.6000000000e1_dp * t478 * t489 * t1829
            t1835 = 0.1800000000e2_dp * t478 * t1176 * t561
            t1836 = t1800 - t1804 - t1807 - t1810 + t1814 - t1817 + 0.8999999998e1_dp  &
            * t178 * t62 * t1688 + t1824 - t1828 + t1832 - t1835
            t1837 = t68 * t1836
            t1838 = t1837 * t65
            t1842 = -t1760 - t1762 - t1765 - t1768 - 0.5555555558e-1_dp * t192 * t1769  &
            - t1773 - t1776 - t1779 - 0.5555555558e-1_dp * t503 * t1525 - t1784  &
            - t1787 - t1789 - 0.1111111112e0_dp * t516 * t1790 - t1794 + t1797  &
            - 0.5555555558e-1_dp * t192 * t1838 + t1836 * t65
            t1843 = t1758 + t1842
            t1849 = 0.8333333335e-1_dp * t429 * t166 * t1571
            t1852 = 0.8333333335e-1_dp * t429 * t679 * t290
            t1854 = t165 * t466 * t1483
            t1858 = 0.1666666667e0_dp * t445 * t446 * t1825
            t1861 = 0.1800000000e2_dp * t658 * t575 * t425
            t1864 = 0.1666666667e0_dp * t445 * t1292 * t668
            t1865 = t413 * t561
            t1868 = 0.8333333335e-1_dp * t1246 * t446 * t1865
            t1871 = 2 * t147 * t393 * t672
            t1874 = 0.3333333334e0_dp * t170 * t689 * t25
            t1875 = t114 * t585
            t1878 = 0.4444444445e0_dp * t170 * t1875 * t451
            t1881 = 0.1666666667e0_dp * t165 * t437 * t649
            t1884 = -t1496 + 0.9000000000e1_dp * t1498 + 0.6000000000e1_dp * t1502
            t1885 = t1884 * t149
            t1886 = t1885 * t160
            t1891 = 0.5625000000e1_dp * t1302 * t401 * t1410
            t1892 = t19 * t649
            t1895 = 0.2250000000e1_dp * t400 * t1892 * t237
            t1898 = 0.2250000000e1_dp * t400 * t401 * t1447
            t1901 = 0.7500000000e0_dp * t400 * t406 * t1483
            t1902 = t68 * t649
            t1905 = 0.5000000000e0_dp * t405 * t1902 * t17
            t1908 = 0.2250000000e1_dp * t400 * t1309 * t561
            t1912 = t399 * t245
            t1914 = 0.2250000000e1_dp * t1912 * t1822
            t1917 = 0.1500000000e1_dp * t405 * t156 * t1825
            t1920 = 0.1000000000e1_dp * t405 * t417 * t1829
            t1923 = 0.1500000000e1_dp * t405 * t1333 * t561
            t1924 = -t1891 + t1895 + t1898 + t1901 - t1905 + t1908 - 0.1500000000e1_dp  &
            * t151 * t50 * t1688 - t1914 + t1917 - t1920 + t1923
            t1925 = t149 * t1924
            t1928 = t56 * t1531
            t1932 = t27 * t1531
            t1938 = 0.1666666667e0_dp * t165 * t462 * t1447
            t1939 = 0.3333333334e0_dp * t55 * t57 * t1843 - t1849 - t1852 - 0.5555555555e-1_dp  &
            * t1854 + t1858 + t1861 + t1864 - t1868 + t1871 + t1874  &
            - t1878 + t1881 + (2 * t147 * t1886) + (2 * t147 * t1925)  &
            + 0.1666666667e0_dp * t165 * t1928 * t140 + 0.3333333334e0_dp * t170 * t1932  &
            * t172 + t1938
            t1940 = t27 * t735
            t1943 = 0.1111111111e0_dp * t170 * t1940 * t17
            t1946 = 0.1666666667e0_dp * t165 * t1392 * t561
            t1952 = 0.5555555557e-1_dp * t165 * t689 * t433
            t1954 = t1373 * t446 * t1829
            t1958 = 0.8333333335e-1_dp * t429 * t166 * t1553
            t1961 = 0.1666666667e0_dp * t165 * t679 * t379
            t1962 = t585 * rho
            t1965 = 0.1666666667e0_dp * t445 * t1962 * t413
            t1966 = t56 * t735
            t1969 = 0.1666666667e0_dp * t165 * t1966 * t237
            t1971 = t657 * t392 * t59
            t1972 = t177 * t20
            t1973 = t561 * t149
            t1974 = t1973 * t160
            t1975 = t1972 * t1974
            t1977 = 0.1800000000e2_dp * t1971 * t1975
            t1980 = 0.2777777778e-1_dp * t1246 * t1247 * t1483
            t1981 = t1367 * t57
            t1984 = 0.1250000000e0_dp * t1981 * t1247 * t1410
            t1987 = 0.8333333335e-1_dp * t429 * t437 * t593
            t1988 = t649 * t17
            t1991 = 0.5555555557e-1_dp * t165 * t171 * t1988
            t1994 = 0.8333333335e-1_dp * t429 * t166 * t1557
            t1997 = 0.1666666667e0_dp * t165 * t171 * t668
            t2000 = 0.8333333335e-1_dp * t429 * t462 * t1410
            t2001 = -t1943 + t1946 + 0.1666666667e0_dp * t165 * t166 * t1688 - t1952  &
            - 0.2222222222e0_dp * t1954 - t1958 + t1961 + t1965 + t1969 + t1977  &
            + t1980 + t1984 - t1987 - t1991 - t1994 + t1997 - t2000
            t2002 = t1939 + t2001
            t2006 = -t1541 - t1543 - t1545 + t1549 - t1552 - t1556 - t1560 - t1563  &
            - 0.3333333336e0_dp * t83 * t1565 - t1570 - t1574 - t1577 - 0.3333333336e0_dp  &
            * t83 * t1689 - t1693 - t1696 - t1699 - t1702 - t1705 -  &
            t1708 - 0.2222222224e0_dp * t144 * t46 * t2002
            e_ndrho_rho_spin(ii) = e_ndrho_rho_spin(ii) +  (-t1406 - t1408 - t1413 - t1417 - t1451 - t1455 - 0.4444444448e0_dp  &
            * t43 * t1536 + t8 * t2006) * sx
            t2009 = t566 * t393
            t2013 = t566 * t149 * t47
            t2014 = t398 * t19
            t2018 = t150 * t68
            t2022 = t19 * t1447
            t2027 = t65 * t237
            t2028 = t575 * t2027
            t2031 = t56 * t561
            t2032 = t2031 * t526
            t2035 = t20 * t1447
            t2042 = -0.1800000000e2_dp * t473 * t2028 + 0.5999999999e1_dp * t478 * t2032  &
            + 0.8999999998e1_dp * t478 * t2035 * t65 + t1509 - t1512 + t1517  &
            - t1523 + t1515 - 0.5000000001e0_dp * t60 * t1447 * t65 - t1530
            t2046 = -0.3000000000e1_dp * t2009 * t569 + 0.4500000000e1_dp * t2013 * t2014  &
            * t1410 - 0.1000000000e1_dp * t2013 * t2018 * t1483 - 0.3000000000e1_dp  &
            * t567 * t151 * t2022 - t1476 - 0.5555555557e-1_dp * t1485 + t1492  &
            + t1482 + t1479 - t1489 + 0.3333333334e0_dp * t55 * t57 * t2042
            t2050 = t45 * t2046
            t2054 = t786 * t640
            t2057 = t1636 + t1638 - t1644 - t1647 - t1649 + t1654 - t1656 + t1664  &
            + t1662 - t1660 + t1671 + t1679 - t1675 - 0.240e2_dp * t2054 * t1060
            t2058 = t105 * t2057
            t2061 = t1579 + t1584 - t1581 - t1586 + t1588 - t1591 - t1596 + t1594  &
            + t1601 - t1605 + t1607 + t1615 - t1630 - t1610 - t1620 + t1634  &
            + t1628 + t1613 + t1624 - 0.20e1_dp * t102 * t2058
            t2073 = -t1891 + t1901 + t1908 + t1898 + t1895 - t1905 - 0.1500000000e1_dp  &
            * t151 * t50 * t2061 - t1914 - t1920 + t1923 + t1917
            t2074 = t149 * t2073
            t2077 = t161 * t17
            t2081 = t657 * t472
            t2082 = t161 * t237
            t2086 = -t1849 - t1852 - 0.5555555557e-1_dp * t1854 + t1858 + t1861 +  &
            t1864 - t1868 + t1871 + t1874 + 0.3333333334e0_dp * t170 * t27 * t2042  &
            * t172 - t1878 + (2 * t147 * t2074) + t1881 + t1938 - t1943  &
            + 0.1200000000e2_dp * t658 * t2031 * t2077 + t1946 - 0.3600000000e2_dp *  &
            t2081 * t575 * t2082
            t2088 = t56 * t2042
            t2095 = t1447 * t65
            t2098 = t140 * t471
            t2099 = t709 * t237
            t2104 = t62 * t2095
            t2107 = -t1714 + t1717 + t1719 - t1724 + t1728 - 0.5000000001e0_dp * t716  &
            * t2095 + 0.1000000000e1_dp * t2098 * t2099 - 0.1000000001e1_dp * t723  &
            * t2095 - t1730 + t1735 + t1743 + 0.8999999998e1_dp * t708 * t2104  &
            + t1752 - t1754 - t1757 - t1760
            t2113 = t177 * pi
            t2114 = t707 * t2113
            t2117 = t471 * pi
            t2118 = t707 * t2117
            t2124 = t1800 - t1810 - t1817 - t1807 - t1804 + t1814 + 0.8999999998e1_dp  &
            * t178 * t62 * t2061 + t1824 + t1832 - t1835 - t1828
            t2126 = t68 * t2124 * t65
            t2130 = -0.5000000001e0_dp * t732 * t2095 - 0.5555555558e-1_dp * t58 * t2061  &
            * t71 - t1762 + t1765 + t1768 + 0.5999999999e1_dp * t2114 * t2032  &
            - t1773 - t1776 - 0.1800000000e2_dp * t2118 * t2028 - 0.5555555558e-1_dp  &
            * t192 * t2126 - t1779 - t1784 - t1787 - t1789 - t1794 + t1797 +  &
            t2124 * t65
            t2131 = t2107 + t2130
            t2138 = -t1952 - 0.2222222223e0_dp * t1954 - t1958 + t1961 + t1965 + t1969  &
            + t1977 + t1980 + 0.1666666667e0_dp * t165 * t2088 * t140 + 0.1800000000e2_dp  &
            * t658 * t2035 * t161 + 0.3333333334e0_dp * t55 * t57 * t2131  &
            + t1984 + 0.1666666667e0_dp * t165 * t166 * t2061 - t1987 - t1991  &
            - t1994 + t1997 - t2000
            t2139 = t2086 + t2138
            t2143 = -t1541 + t1549 - t1563 - t1574 - t1560 - t1543 - t1552 - 0.3333333336e0_dp  &
            * t83 * t2050 * t140 - t1577 - t1545 - t1556 - t1570  &
            - 0.3333333336e0_dp * t83 * t84 * t2061 - t1693 - t1699 - t1708 - t1705  &
            - t1696 - t1702 - 0.2222222224e0_dp * t144 * t46 * t2139
            e_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_spin(ii) +  (-t1406 - t1413 - t1455 - t1451 - t1408 - t1417 -&
            0.4444444448e0_dp * t43 * t46 * t2046 + t8 * t2143) * sx
            t2147 = lambda * t48 * omega
            t2160 = t754 * t237
            t2169 = t237 * t804
            t2182 = 0.1e1_dp / t787 / lambda
            t2183 = t2182 * t944
            t2187 = t949 * t233
            t2190 = t953 * t772
            t2193 = t316 * t800
            t2196 = t772 * t233
            t2197 = t974 * t2196
            t2200 = t233 * t800
            t2201 = t327 * t2200
            t2204 = t10 * t106
            t2205 = 0.1e1_dp / t2204
            t2220 = beta * t785 * t6
            t2226 = -0.6222222223e2_dp * t24 * t2205 * t1002 * t31 - 0.2115555556e3_dp  &
            * t7 * t2183 * t132 + 0.1315555556e3_dp * t786 / t787 / t85 * t1019  &
            * t371 - 0.4266666668e2_dp * t2220 / t1014 / t94 * t1058
            t2227 = t105 * t2226
            t2230 = -0.9125925923e2_dp * t7 * t2183 * t36 - 0.5866666667e2_dp * t962  &
            * t2187 - 0.3200000001e2_dp * t321 * t2190 + 0.1600000000e2_dp * t321 *  &
            t2193 - 0.120e2_dp * t102 * t2197 + 0.120e2_dp * t102 * t2201 - 0.20e1_dp *  &
            t102 * t2227
            t2257 = 0.5400000000e2_dp * t1142 * t62 * t2160
            t2261 = 0.3600000000e2_dp * t473 * t56 * t754 * t17
            t2264 = 0.5400000000e2_dp * t473 * t1493 * t804
            t2268 = 0.6000000000e1_dp * t478 * t27 * t237 * t14
            t2272 = 0.1800000000e2_dp * t478 * t56 * t804 * t17
            t2274 = t178 * t62 * t2230
            t2276 = pi * t114
            t2279 = 0.2666666667e1_dp * t61 * t2276 * t86
            t2280 = t2257 - t2261 - t2264 - t2268 + t2272 + 0.8999999998e1_dp * t2274  &
            - t2279
            t2281 = t2280 * t65
            t2282 = t881 * t267
            t2285 = t865 * t266
            t2287 = t68 * t2285 * t65
            t2292 = t863 * t266
            t2295 = t2285 * t65
            t2296 = t58 * t2230
            t2299 = t194 * t865
            t2300 = t2299 * t526
            t2304 = t68 * t2280 * t65
            t2307 = t22 * t266
            t2308 = t2307 * t1200
            t2311 = t194 * t863
            t2312 = t2311 * t526
            t2319 = t1099 * t86
            t2326 = -0.1851851853e0_dp * t508 * t874 + t2281 - 0.1666666668e0_dp * t192  &
            * t2282 - 0.5555555558e-1_dp * t192 * t2287 - 0.1666666668e0_dp * t508  &
            * t886 + (3 * t2292 * t65) + t2295 - 0.5555555558e-1_dp * t2296  &
            * t71 + 0.1111111112e0_dp * t192 * t2300 - 0.5555555558e-1_dp * t192 *  &
            t2304 - 0.1851851853e0_dp * t192 * t2308 + 0.1111111112e0_dp * t192 * t2312  &
            + 0.2222222223e0_dp * t508 * t878 - 0.1666666668e0_dp * t508 * t882  &
            + 0.1646090535e0_dp * t192 * t2319 - 0.1666666668e0_dp * t1129 * t275 +  &
            0.1111111112e0_dp * t1129 * t271
            t2351 = t44 * t22
            t2368 = t237 * t14
            t2378 = t804 * t17
            t2382 = t754 * t17
            t2392 = (2 * t147 * t1284 * t249) + (2 * t147 * t1388 * t249)  &
            - 0.3456790122e0_dp * t170 * t1268 * t86 + 0.4444444444e0_dp * t170 *  &
            t1089 * t14 + 0.1666666667e0_dp * t165 * t166 * t2230 + 0.5000000001e0_dp  &
            * t165 * t437 * t804 + 0.3333333334e0_dp * t55 * t57 * t2326 + 0.2e1_dp  &
            * (t147) * t149 * (-0.5625000000e1_dp * t1301 * t50 * t2160 +  &
            0.2250000000e1_dp * t400 * t68 * t754 * t17 + 0.6750000000e1_dp * t400 *  &
            t1459 * t804 + 0.1000000000e1_dp * t405 * t194 * t237 * t14 - 0.1500000000e1_dp  &
            * t405 * t68 * t804 * t17 - 0.1500000000e1_dp * t151 * t50 *  &
            t2230 + 0.1111111111e1_dp * t49 * t2351 * t86) - 0.3333333333e0_dp * t170  &
            * t1254 * t17 + (4 * t147 * t393 * t826) + 0.5000000001e0_dp *  &
            t165 * t1275 * t237 + 0.2222222222e0_dp * t165 * t450 * t2368 - 0.3333333333e0_dp  &
            * t165 * t455 * t833 - 0.2500000000e0_dp * t429 * t166 * t2169  &
            - 0.1666666667e0_dp * t165 * t171 * t2378 + 0.8333333333e-1_dp * t429  &
            * t171 * t2382 - 0.2500000000e0_dp * t429 * t437 * t754 + 0.1250000000e0_dp  &
            * t1368 * t166 * t2160
            t2404 = t45 * t1397
            t2418 = t379 * t754
            t2424 = t1067 * t237
            t2437 = t7 * t792
            t2439 = 0.1e1_dp / t20 / t1436
            t2451 = t940 * t944
            t2455 = t97 * t973
            t2468 = t944 * t104
            t2480 = 0.1e1_dp / t19 / t942
            t2492 = t1014 * rho
            t2537 = t764 * t1002
            t2552 = t360 * t1058
            t2565 = t126 / t131 / t1056 / t130
            t2569 = -0.8088888890e3_dp * t24 * t2182 * t2480 * t31 * rho - 0.1518222222e4_dp  &
            * t2437 * t2439 * t132 * rho + 0.6542222223e3_dp * t786 / t2204  &
            / t19 / t2492 * t218 * t371 * rho - 0.6222222223e2_dp * t24 * t2205  &
            * t1002 * t126 * t132 + 0.5600000000e2_dp * t1024 * t341 * t360  &
            * t132 + 0.3288888889e3_dp * t366 / t790 * t126 * t371 - 0.2400000000e2_dp  &
            * t350 * t114 * t1044 * t132 - 0.1520000000e3_dp * t366 * t1034 *  &
            t360 * t371 - 0.3626666667e3_dp * t1054 * t1034 * t1060 + 0.60e1_dp * t119  &
            * t27 * (0.3640e4_dp / 0.27e2_dp * t1040 * t1002 - 0.3640e4_dp / 0.27e2_dp  &
            * ndrho * t2182 * t2480 * rho) * t132 + 0.8088888890e3_dp * t24 * t2537  &
            * t31 + 0.1518222222e4_dp * t7 * t2451 * t132 - 0.6542222223e3_dp *  &
            t786 * t1012 * t1019 * t371 + 0.2400000000e2_dp * t366 * t368 * t1044  &
            * t371 + 0.9600000000e2_dp * t1055 * t2552 * t220 - 0.1173333333e3_dp  &
            * t1055 * t1059 * t765 + 0.2133333333e3_dp * beta * t785 * t365 * t368  &
            * t2565 * t796
            t2576 = 0.40e1_dp * t102 * t327 * t136 * t2226 + 0.3519999999e3_dp * t93  &
            * t2190 - 0.1551407408e4_dp * t2437 * t2439 * t36 * rho + 0.960e2_dp *  &
            t1597 * t1598 * t2200 + 0.120e2_dp * t102 * t327 * t375 * t800 + 0.1551407408e4_dp  &
            * t7 * t2451 * t36 - 0.960e2_dp * t93 * t2455 * rho * t2196  &
            + 0.120e2_dp * t102 * t327 * t1063 * t233 + 0.1600000000e2_dp * t321  &
            * t316 * t1063 + 0.9125925923e2_dp * t7 * t2182 * t2468 * t136 + 0.1760000000e3_dp  &
            * t307 * t949 * t958 - 0.360e2_dp * t1616 * t1617 * t2200  &
            - 0.20e1_dp * t102 * t105 * t2569 - 0.6400000001e2_dp * t321 * t953 * t979
            t2585 = 0.1e1_dp / t972 / t35
            t2611 = t311 * t326
            t2621 = -0.3200000001e2_dp * t321 * t953 * t983 - 0.160e2_dp * t93 * t316  &
            * rho * t2226 + 0.480e2_dp * t102 * t22 * t2585 * t136 * t2196 - 0.1760000000e3_dp  &
            * t93 * t2193 + 0.960e2_dp * t302 * t2197 + 0.160e2_dp * t302  &
            * t2227 - 0.960e2_dp * t302 * t2201 - 0.8213333332e3_dp * t941 * t2468  &
            * t317 + 0.8213333332e3_dp * t307 * t2187 + 0.9600000002e2_dp * t321 *  &
            t2455 * t975 - 0.360e2_dp * t102 * t974 * t375 * t772 - 0.3519999999e3_dp  &
            * t307 * t2611 * t954 + 0.1173333334e3_dp * t962 * t2611 * t328 -  &
            0.5866666667e2_dp * t962 * t949 * t375
            t2622 = t2576 + t2621
            t2631 = t282 * t140
            t2635 = -0.1000000001e1_dp * t55 * t923 - 0.3333333336e0_dp * t83 * t384  &
            * t2230 - 0.5000000004e0_dp * t165 * t897 * t804 - 0.6666666672e0_dp * t383  &
            * t2404 - 0.1000000001e1_dp * t83 * t1071 * t804 - 0.5000000004e0_dp  &
            * t164 * t1076 - 0.1000000001e1_dp * t289 * t294 * t914 - 0.5000000004e0_dp  &
            * t289 * t294 * t918 + 0.2500000002e0_dp * t906 * t84 * t2418 -  &
            0.2000000001e1_dp * t55 * t1079 - 0.5000000004e0_dp * t289 * t84 * t2424  &
            + 0.2500000002e0_dp * t1546 * t897 * t2169 - 0.3333333336e0_dp * t83 *  &
            t84 * t2622 + 0.2500000002e0_dp * t906 * t294 * t907 - 0.1000000001e1_dp  &
            * t55 * t1082 - 0.1000000001e1_dp * t165 * t2631 * t237
            t2636 = t78 * t379
            t2645 = t379 * t804
            t2655 = t140 * t2230
            t2707 = t74 * t379
            t2720 = t56 * t1239
            t2724 = 0.5000000001e0_dp * t165 * t1392 * t804 - 0.1333333334e1_dp * t165  &
            * t450 * t1172 - 0.2666666667e1_dp * t170 * t1089 * t225 + 0.8333333333e-1_dp  &
            * t429 * t466 * t2382 + 0.1000000000e1_dp * t170 * t1254 * t25  &
            + 0.1666666667e0_dp * t165 * t166 * t2622 + 0.1666666667e0_dp * t165 *  &
            t56 * t2326 * t140 - 0.2500000000e0_dp * t429 * t166 * t2424 - 0.3333333333e0_dp  &
            * t165 * t455 * t1085 - 0.1728395062e0_dp * t165 * t1268 *  &
            t140 * t86 - 0.2500000000e0_dp * t1246 * t446 * t413 * t804 + 0.1666666666e0_dp  &
            * t1246 * t2707 * t833 + 0.1250000000e0_dp * t1368 * t462 * t2160  &
            - 0.1666666667e0_dp * t165 * t466 * t2378 + (6 * t147 * t1283  &
            * t394) + 0.5000000001e0_dp * t165 * t2720 * t237
            t2732 = rho * t2205
            t2759 = t25 * t2230
            t2763 = t1367 * t444
            t2779 = t25 * t2160
            t2783 = -0.2500000000e0_dp * t429 * t1392 * t754 + (6 * t147 * t393  &
            * t1347) - 0.3456790123e1_dp * t170 * t1002 * t74 * t2732 + 0.5000000001e0_dp  &
            * t165 * t1275 * t379 - 0.1333333334e1_dp * t170 * t114 * t889  &
            * t451 + 0.3333333334e0_dp * t170 * t27 * t2326 * t172 - 0.2500000000e0_dp  &
            * t429 * t171 * t1168 + 0.3750000000e0_dp * t1981 * t1247 * t2169  &
            + 0.5000000001e0_dp * t165 * t171 * t1093 + 0.8333333334e-1_dp * t1246  &
            * t1247 * t2378 + 0.1666666667e0_dp * t445 * t446 * t2759 - 0.1250000000e0_dp  &
            * t2763 * t1247 * t2382 + 0.5000000001e0_dp * t165 * t437 * t1067  &
            + 0.2222222222e0_dp * t165 * t1089 * t1360 + 0.1000000000e1_dp * t165  &
            * t455 * t413 - 0.2500000000e0_dp * t429 * t462 * t2169 + 0.1250000000e0_dp  &
            * t2763 * t446 * t2779
            t2808 = t19 * t1067
            t2816 = t225 * t754
            t2821 = 0.1e1_dp / t41 / t1140
            t2823 = t47 * t2821 * t44
            t2833 = -0.1687500000e2_dp * t1302 * t401 * t2169 + 0.6750000000e1_dp * t400  &
            * t1309 * t804 + 0.4500000000e1_dp * t400 * t1320 * t833 - 0.1500000000e1_dp  &
            * t400 * t1316 * t2368 + 0.1000000000e1_dp * t405 * t194 * t379  &
            * t14 + 0.8888888888e1_dp * t155 * t97 * rho * t2205 + 0.6750000000e1_dp  &
            * t400 * t2808 * t237 - 0.1500000000e1_dp * t405 * t68 * t1067  &
            * t17 + 0.4500000000e1_dp * t400 * t417 * t2816 + 0.1968750000e2_dp * t2823  &
            * t401 * t2160 + 0.1500000000e1_dp * t405 * t156 * t2759 + 0.2250000000e1_dp  &
            * t400 * t401 * t2230
            t2843 = t777 * t237
            t2853 = t172 * t2169
            t2869 = t225 * t804
            t2873 = t194 * t225
            t2877 = -0.1500000000e1_dp * t151 * t50 * t2622 + 0.5625000000e1_dp * t1302  &
            * t156 * t2779 - 0.6750000000e1_dp * t400 * t1333 * t754 + 0.4999999999e1_dp  &
            * t405 * t1340 * t2843 - 0.1000000000e2_dp * t49 * t2351 * t777  &
            - 0.5625000000e1_dp * t1302 * t406 * t2382 - 0.6750000000e1_dp * t1912  &
            * t2853 - 0.5555555555e0_dp * t405 * t22 * t140 * t86 + 0.4500000000e1_dp  &
            * t405 * t1333 * t804 + 0.2250000000e1_dp * t400 * t406 * t2378 -  &
            0.1687500000e2_dp * t1302 * t1309 * t754 - 0.3000000000e1_dp * t405 * t417  &
            * t2869 - 0.6000000000e1_dp * t405 * t2873 * t237
            t2885 = t428 * t1372
            t2935 = -0.2500000000e0_dp * t429 * t166 * t2645 + (2 * t147 * t149  &
            * (t2833 + t2877)) + 0.5000000001e0_dp * t445 * t1292 * t1093 + 0.3333333333e0_dp  &
            * t2885 * t446 * t2816 + 0.2222222222e0_dp * t165 * t1377  &
            * t2368 - 0.3333333333e0_dp * t165 * t1264 * t833 + 0.2e1_dp * (t147)  &
            * (t2257 - t2261 - t2264 - t2268 + t2272 + 0.9000000000e1_dp * t2274  &
            - t2279) * (t149) * t160 - 0.3333333333e0_dp * t170 * t27 * t1239  &
            * t17 - 0.8333333335e-1_dp * t429 * t166 * t2655 + (6 * t147  &
            * t1284 * t424) + 0.2e1_dp * (t147) * t1387 * t392 * (t149) *  &
            t160 - 0.2500000000e0_dp * t1246 * t1292 * t1168 - 0.1666666667e0_dp *  &
            t165 * t1254 * t433 + (6 * t147 * t1388 * t424) + 0.1666666667e0_dp  &
            * t165 * t462 * t2230 + 0.2222222222e0_dp * t165 * t450 * t379 * t14  &
            - 0.5000000001e0_dp * t429 * t437 * t914
            t2937 = t164 * t45 * t341
            t2992 = 0.1e1_dp / t1140 / t40
            t2994 = t59 * t2992 * pi
            t2998 = t27 * t225
            t3026 = 0.1620000000e3_dp * t1143 * t474 * t2169 + 0.1080000000e3_dp * t1143  &
            * t479 * t2382 - 0.2160000000e3_dp * t2994 * t474 * t2160 + 0.3600000000e2_dp  &
            * t478 * t2998 * t237 + 0.1620000000e3_dp * t1143 * t1150 *  &
            t754 + 0.2666666667e1_dp * t478 * t114 * t140 * t86 - 0.3600000000e2_dp  &
            * t473 * t479 * t2378 + 0.1080000000e3_dp * t1821 * t2853 - 0.7200000000e2_dp  &
            * t473 * t1161 * t833 - 0.1800000000e2_dp * t473 * t474 * t2230  &
            + 0.8999999998e1_dp * t178 * t62 * t2622 - 0.5400000000e2_dp * t478 *  &
            t1176 * t804
            t3030 = t20 * t1067
            t3069 = -0.1800000000e2_dp * t478 * t183 * t2759 - 0.5400000000e2_dp * t473  &
            * t3030 * t237 - 0.6000000000e1_dp * t478 * t27 * t379 * t14 - 0.5400000000e2_dp  &
            * t473 * t1150 * t804 + 0.1800000000e2_dp * t478 * t489 *  &
            t2869 + 0.1080000000e3_dp * t473 * t1176 * t754 + 0.1800000000e2_dp * t478  &
            * t56 * t1067 * t17 - 0.3600000000e2_dp * t473 * t489 * t2816 - 0.1866666667e2_dp  &
            * t182 * t357 * t2205 + 0.2400000000e2_dp * t61 * t2276  &
            * t777 - 0.2400000000e2_dp * t478 * t123 * t2843 - 0.1080000000e3_dp *  &
            t1143 * t183 * t2779 + 0.1200000000e2_dp * t473 * t1157 * t2368
            t3070 = t3026 + t3069
            t3085 = t1210 * t67
            t3093 = (3 * t496 * t865 * t65) + (3 * t1189 * t266 * t65)  &
            + 0.1111111112e1_dp * (t192) * (t2307) * (t65) * (t225)  &
            + (t187 * t2280 * t65) + (3 * t496 * t863 * t65) + 0.2222222223e0_dp  &
            * t516 * t1111 * t526 + (t187 * t2285 * t65) + (t3070  &
            * t65) - 0.1851851853e0_dp * t1193 * t498 * t1200 + 0.2222222223e0_dp  &
            * t508 * t1119 + 0.1111111112e0_dp * (t192) * t194 * (t1189)  &
            * t526 + 0.2222222223e0_dp * t1137 * t1221 - 0.1666666668e0_dp * t1129  &
            * t534 + 0.1111111112e0_dp * t3085 * t271 - 0.5555555558e-1_dp * (t192)  &
            * (t199) * (t2281) - 0.1851851853e0_dp * t508 * t1201
            t3097 = t66 * t69
            t3098 = t1109 * t267
            t3129 = t2296 * t67
            t3141 = -0.3333333336e0_dp * t192 * t2299 * t1125 - 0.1666666668e0_dp * t3097  &
            * t3098 + 0.1111111112e0_dp * t503 * t2300 + 0.1111111112e0_dp * t1129  &
            * t527 - 0.1666666668e0_dp * t1129 * t531 - 0.5555555558e-1_dp * t192  &
            * t68 * t3070 * t65 + 0.1646090535e0_dp * t192 * t97 * t187 * t65 *  &
            t86 + 0.5555555559e0_dp * t1129 * t512 - 0.1666666668e0_dp * t192 * t530  &
            * t864 - 0.1666666668e0_dp * t1106 * t882 - 0.1666666668e0_dp * t192 *  &
            t530 * t866 - 0.1481481482e1_dp * t192 * t1099 * t777 - 0.1111111112e0_dp  &
            * t3129 * t196 + 0.5555555559e0_dp * t1193 * t866 * t451 - 0.1666666668e0_dp  &
            * t503 * t2282 - 0.1666666668e0_dp * t3085 * t275 + 0.1111111112e0_dp  &
            * t503 * t2312
            t3181 = -0.1851851853e0_dp * t1106 * t874 - 0.1666666668e0_dp * t508 * t1235  &
            - 0.3333333336e0_dp * t508 * t1224 - 0.1111111112e0_dp * t516 * t2281  &
            * t172 + 0.1111111112e0_dp * t516 * t1097 * t526 + 0.2222222223e0_dp *  &
            t1106 * t878 - 0.1851851853e0_dp * t503 * t2308 - 0.3333333336e0_dp * t1137  &
            * t1216 + (3 * t3098) + 0.1111111112e0_dp * t516 * t1109 * t526  &
            - 0.5555555558e-1_dp * t503 * t2287 - 0.3333333336e0_dp * t192 * t2311  &
            * t1125 - 0.6666666672e0_dp * t508 * t1126 + 0.1646090535e0_dp * t503  &
            * t2319 - 0.5555555558e-1_dp * t3129 * t200 - 0.5555555558e-1_dp * t58  &
            * t2622 * t71 - 0.1666666668e0_dp * t1106 * t886
            t3232 = -0.1666666668e0_dp * t192 * t1206 * t267 - 0.5555555558e-1_dp * t192  &
            * t199 * t2295 - 0.1666666668e0_dp * t508 * t1122 - 0.3333333336e0_dp  &
            * t1129 * t520 - 0.1111111112e0_dp * t516 * t2295 * t172 - 0.1481481482e1_dp  &
            * t66 * t67 * t97 * t267 * t1100 - 0.1851851853e0_dp * t192 *  &
            t22 * t496 * t1200 - 0.3333333336e0_dp * t516 * t2292 * t65 * rho *  &
            t25 + 0.1810699590e1_dp * t192 * t311 * t70 * t2732 - 0.5555555558e-1_dp  &
            * t503 * t2304 - 0.1481481482e1_dp * t508 * t1101 - 0.3333333336e0_dp *  &
            t867 * t515 * t517 + 0.1111111112e1_dp * t268 * t1192 * t1194 + 0.5555555559e0_dp  &
            * t1193 * t864 * t451 + 0.1111111112e1_dp * t508 * t1213 -  &
            0.1666666668e0_dp * t508 * t1207 - 0.3333333336e0_dp * t1137 * t1134
            t3244 = t278 * t140
            t3248 = t889 * rho
            t3262 = omega * t1300
            t3264 = t140 * t2160
            t3268 = 0.1555555556e1_dp * t2937 * t446 * t2843 + 0.3111111111e1_dp * t170  &
            * t1268 * t777 - 0.2500000000e0_dp * t429 * t1275 * t290 - 0.1333333334e1_dp  &
            * t1373 * t1292 * t1172 + 0.3111111111e1_dp * t170 * t341 * t278  &
            * t1100 - 0.3456790122e0_dp * t170 * t341 * t203 * t86 + 0.3750000000e0_dp  &
            * t1368 * t166 * t2418 + 0.4444444444e0_dp * t170 * t114 * t537  &
            * t14 + 0.3333333334e0_dp * t55 * t57 * (t3093 + t3141 + t3181 + t3232)  &
            - 0.1111111111e0_dp * t2885 * t1247 * t2368 - 0.6666666668e0_dp * t1373  &
            * t446 * t2869 + 0.1666666667e0_dp * t1246 * t3244 * t833 + 0.5000000001e0_dp  &
            * t445 * t3248 * t413 - 0.1666666667e0_dp * t165 * t171 *  &
            t1067 * t17 - 0.2500000000e0_dp * t429 * t437 * t918 + 0.3750000000e0_dp  &
            * t1368 * t437 * t907 - 0.3125000000e0_dp * t3262 * t45 * t166 * t3264
            t3278 = rho * t398
            t3288 = -0.1000000001e1_dp * t165 * t2636 * t237 - 0.2000000001e1_dp * t55  &
            * t926 - 0.1000000001e1_dp * t83 * t922 * t379 - 0.5000000004e0_dp * t289  &
            * t84 * t2645 + 0.2500000002e0_dp * t429 * t897 * t754 - 0.1000000001e1_dp  &
            * t83 * t294 * t1067 - 0.1666666668e0_dp * t289 * t84 * t2655  &
            - 0.1000000001e1_dp * t83 * t2404 * t237 - 0.5000000004e0_dp * t289 * t922  &
            * t290 - 0.5000000004e0_dp * t289 * t384 * t2169 - 0.5000000004e0_dp  &
            * t289 * t1071 * t754 - 0.2222222224e0_dp * t144 * t46 * (t2724 + t2783  &
            + t2935 + t3268) - 0.3333333336e0_dp * t83 * t45 * t2392 * t140 -  &
            0.1250000001e0_dp * t3278 * omega * t84 * t3264 - 0.1000000001e1_dp * t55  &
            * t1068 + 0.8333333340e-1_dp * t906 * t384 * t2160
            t3291 = -0.1000000001e1_dp * t2147 * t84 * t754 - 0.4000000002e1_dp * t747  &
            * t294 * t237 - 0.2000000001e1_dp * t747 * t84 * t804 - 0.1333333334e1_dp  &
            * t214 * t894 + 0.1666666668e0_dp * t4 * t150 * t46 * t78 * t2160  &
            - 0.1000000001e1_dp * t753 * t46 * t282 * t754 - 0.1000000001e1_dp * t1409  &
            * t84 * t2169 - 0.2000000001e1_dp * t217 * t46 * t893 * t237 - 0.2000000001e1_dp  &
            * t217 * t46 * t282 * t804 - 0.6666666672e0_dp * t217 *  &
            t46 * t78 * t2230 - 0.4444444448e0_dp * t43 * t46 * t2392 + t8 * (t2635  &
            + t3288)
          END IF
          IF( grad_deriv >= 3 .OR. grad_deriv == -3 ) THEN
            e_rho_rho_rho_spin(ii) = e_rho_rho_rho_spin(ii) +  t3291 * sx
            t3305 = t4 * t428
            t3306 = t754 * t561
            t3314 = t237 * t1447
            t3329 = t804 * t561
            t3339 = t949 * t557
            t3344 = t953 * t1425
            t3347 = t316 * t1443
            t3352 = t772 * t557
            t3353 = t974 * t3352
            t3356 = t233 * t1443
            t3357 = t327 * t3356
            t3362 = t800 * t557
            t3363 = t327 * t3362
            t3366 = beta * t777
            t3378 = 0.1866666667e2_dp * t3366 * t342 + 0.9866666667e2_dp * t24 * t782  &
            - 0.8266666668e2_dp * t366 * t797 + 0.3200000001e2_dp * t1054 / t1014 /  &
            t12 * t1058
            t3379 = t105 * t3378
            t3382 = 0.3911111110e2_dp * t24 * t766 - 0.1955555555e2_dp * t962 * t3339  &
            + 0.2133333334e2_dp * t339 * t769 - 0.2133333334e2_dp * t321 * t3344 +  &
            0.1066666667e2_dp * t321 * t3347 + 0.80e1_dp * t112 * t773 - 0.120e2_dp * t102  &
            * t3353 + 0.80e1_dp * t102 * t3357 - 0.40e1_dp * t112 * t801 + 0.40e1_dp  &
            * t102 * t3363 - 0.20e1_dp * t102 * t3379
            t3449 = t1447 * t17
            t3453 = (2 * t147 * t1885 * t249) + 0.1800000000e2_dp * (t1971)  &
            * (t1972) * (t1973) * (t249) + (2 * t147 * t393  &
            * t1470) + 0.1800000000e2_dp * t658 * t575 * t827 + 0.2e1_dp * (t147)  &
            * t149 * (-0.5625000000e1_dp * t1302 * t19 * t754 * t561 + 0.4500000000e1_dp  &
            * t400 * t1459 * t1447 + 0.1500000000e1_dp * t400 * t815 * t1483  &
            - 0.1000000000e1_dp * t405 * t68 * t1447 * t17 + 0.2250000000e1_dp *  &
            t400 * t19 * t804 * t561 - 0.1500000000e1_dp * t151 * t50 * t3382 + 0.3333333334e0_dp  &
            * t405 * t194 * t14 * t561) + 0.1250000000e0_dp * t1368  &
            * t166 * t3306 - 0.8333333335e-1_dp * t429 * t679 * t754 - 0.1666666667e0_dp  &
            * t429 * t166 * t3314 + 0.5555555555e-1_dp * t1246 * t74 * t237  &
            * t1483 - 0.1111111111e0_dp * t165 * t689 * t833 - 0.1111111111e0_dp *  &
            t165 * t171 * t3449
            t3472 = t14 * t561
            t3491 = 0.5400000000e2_dp * t1143 * t20 * t754 * t561
            t3494 = 0.3600000000e2_dp * t473 * t1493 * t1447
            t3497 = 0.2400000000e2_dp * t473 * t852 * t1483
            t3501 = 0.1200000000e2_dp * t478 * t56 * t1447 * t17
            t3505 = 0.1800000000e2_dp * t473 * t20 * t804 * t561
            t3507 = t178 * t62 * t3382
            t3512 = 0.2000000000e1_dp * t478 * t27 * t14 * t561
            t3513 = t3491 - t3494 - t3497 + t3501 - t3505 + 0.8999999998e1_dp * t3507  &
            - t3512
            t3514 = t3513 * t65
            t3525 = t58 * t3382
            t3533 = t237 * t1518
            t3543 = t3514 + 0.8999999998e1_dp * t863 * t59 * t177 * t710 + (2  &
            * t267 * t1504) + 0.8999999998e1_dp * t865 * t59 * t177 * t710 - 0.5555555558e-1_dp  &
            * t3525 * t71 - 0.5000000001e0_dp * t804 * t177 * t709 +  &
            0.7407407410e-1_dp * t1755 * t271 + 0.6666666668e0_dp * t3533 * t1722 -  &
            0.1111111112e0_dp * t1755 * t275 - 0.1111111112e0_dp * t508 * t1525 - 0.1000000001e1_dp  &
            * t237 * t266 * t1722
            t3546 = 0.1e1_dp / t12
            t3547 = t60 * t3546
            t3552 = t194 * t1504
            t3553 = t3552 * t526
            t3561 = t68 * t3513 * t65
            t3569 = t274 * t1505
            t3575 = -0.6172839508e-1_dp * t719 * t874 - 0.5555555556e0_dp * t3547 * t709  &
            + 0.7407407410e-1_dp * t719 * t878 + 0.7407407410e-1_dp * t192 * t3553  &
            + 0.6666666668e0_dp * t1519 * t1712 - 0.5555555558e-1_dp * t719 * t882  &
            - 0.5555555558e-1_dp * t192 * t3561 - 0.5000000001e0_dp * t60 * t863 *  &
            t709 - 0.5555555558e-1_dp * t719 * t886 - 0.1111111112e0_dp * t192 * t3569  &
            - 0.5000000001e0_dp * t60 * t865 * t709
            t3576 = t3543 + t3575
            t3580 = -0.1666666667e0_dp * t429 * t437 * t1410 + 0.3333333334e0_dp * t165  &
            * t1928 * t237 + 0.3333333334e0_dp * t165 * t437 * t1447 - 0.8333333335e-1_dp  &
            * t429 * t166 * t3329 + 0.1666666667e0_dp * t165 * t679 * t804  &
            + 0.1666666667e0_dp * t165 * t166 * t3382 + 0.7407407405e-1_dp * t165  &
            * t450 * t3472 + 0.1481481481e0_dp * t170 * t1875 * t14 - 0.1111111111e0_dp  &
            * t165 * t455 * t1483 - 0.2222222222e0_dp * t170 * t1932 * t17 +  &
            0.1666666667e0_dp * t165 * t1275 * t561 + 0.3333333334e0_dp * t55 * t57  &
            * t3576
            t3581 = t3453 + t3580
            t3608 = t973 * rho
            t3618 = 0.640e2_dp * t1597 * t1598 * t3356 - 0.640e2_dp * t1582 * t773 -  &
            0.1173333333e3_dp * t93 * t3347 - 0.3911111110e2_dp * t1000 * t963 - 0.320e2_dp  &
            * t605 * t959 + 0.1173333333e3_dp * t307 * t949 * t1602 + 0.40e1_dp  &
            * t102 * t327 * t136 * t3378 - 0.960e2_dp * t1597 * t3608 * t3352 -  &
            0.40e1_dp * t112 * t1064 + 0.80e1_dp * t112 * t984 + 0.2737777778e3_dp * t307  &
            * t3339
            t3668 = t2552 * t23
            t3673 = t635 * t1058
            t3674 = t3673 * t220
            t3677 = -0.1866666667e3_dp * t1645 * t342 - 0.6346666667e3_dp * t24 * t991  &
            + 0.3946666667e3_dp * t366 * t997 + 0.1866666667e3_dp * beta * t764 *  &
            t1004 + 0.60e1_dp * t550 * t1046 + 0.60e1_dp * t119 * t27 * (-0.280e3_dp /  &
            0.9e1_dp * t987 + 0.280e3_dp / 0.9e1_dp * t2537 * rho) * t132 - 0.60e1_dp * t7  &
            * t640 * t1044 * t371 + 0.1600000000e2_dp * t366 * t368 * t1667 * t371  &
            - 0.4800000000e2_dp * t1680 * t3668 + 0.1600000000e3_dp * t1680 * t1060  &
            + 0.3200000000e2_dp * t1055 * t3674
            t3679 = t2565 * t1439
            t3682 = t24 * t940
            t3711 = -0.1600000000e3_dp * t2220 * t368 * t3679 + 0.6346666667e3_dp * t3682  &
            * t1008 - 0.3946666667e3_dp * t366 * t1012 * t1021 + 0.1866666667e2_dp  &
            * t3366 * t1026 + 0.1866666667e2_dp * t1024 * t341 * t635 * t132  &
            - 0.1706666667e3_dp * t7 * t1036 - 0.1600000000e2_dp * t1429 * t1030 - 0.1600000000e2_dp  &
            * t350 * t114 * t1667 * t132 + 0.6400000000e2_dp * t7 *  &
            t1050 - 0.5066666667e2_dp * t366 * t1034 * t635 * t371 + 0.1520000000e3_dp  &
            * t786 * t1034 * t1681
            t3726 = -0.240e2_dp * t1616 * t1617 * t3356 + 0.160e2_dp * t302 * t3379 +  &
            0.40e1_dp * t102 * t327 * t1063 * t557 + 0.960e2_dp * t302 * t3353 - 0.320e2_dp  &
            * t302 * t3363 - 0.2737777778e3_dp * t941 * t2468 * t608 - 0.120e2_dp  &
            * t102 * t974 * t645 * t772 - 0.2133333334e2_dp * t321 * t953 * t1621  &
            - 0.20e1_dp * t102 * t105 * (t3677 + t3711) + 0.80e1_dp * t102 * t327  &
            * t375 * t1443 - 0.120e2_dp * t1616 * t1617 * t3362 - 0.2133333334e2_dp  &
            * t321 * t953 * t1625
            t3750 = t2585 * t136
            t3763 = 0.2346666666e3_dp * t1589 * t950 + 0.80e1_dp * t102 * t327 * t1684  &
            * t233 + 0.6400000002e2_dp * t7 * t220 * t1618 - 0.240e2_dp * t112 * t976  &
            - 0.160e2_dp * t93 * t316 * rho * t3378 + 0.2346666666e3_dp * t93 *  &
            t3344 + 0.3911111110e2_dp * t962 * t2611 * t614 - 0.5475555556e3_dp * t24  &
            * t930 + 0.480e2_dp * t1616 * t3750 * t3352 - 0.1955555555e2_dp * t962  &
            * t949 * t645 - 0.240e2_dp * t1616 * t973 * t375 * t1425 + 0.2133333334e2_dp  &
            * t339 * t969
            t3794 = -0.640e2_dp * t302 * t3357 - 0.4266666668e2_dp * t339 * t966 + 0.640e2_dp  &
            * t605 * t955 + 0.40e1_dp * t102 * t327 * t645 * t800 + 0.5475555556e3_dp  &
            * t3682 * t946 + 0.1066666667e2_dp * t321 * t316 * t1684 + 0.160e2_dp  &
            * t112 * t980 - 0.2346666666e3_dp * t7 * t929 * t1599 + 0.320e2_dp  &
            * t1597 * t1598 * t3362 - 0.2346666666e3_dp * t605 * t769 + 0.320e2_dp  &
            * t1582 * t801 - 0.2133333334e2_dp * t321 * t953 * t1631
            t3796 = t3618 + t3726 + t3763 + t3794
            t3826 = 0.8333333340e-1_dp * t906 * t597 * t907 + 0.1666666668e0_dp * t429  &
            * t1547 - 0.3333333336e0_dp * t83 * t84 * t3796 - 0.3333333336e0_dp *  &
            t83 * t922 * t649 - 0.3333333336e0_dp * t289 * t294 * t1571 - 0.3333333336e0_dp  &
            * t289 * t1071 * t1410 - 0.6666666672e0_dp * t83 * t1564 * t379  &
            - 0.6666666672e0_dp * t55 * t1700 - 0.3333333336e0_dp * t165 * t207 *  &
            t237 * t561 - 0.6666666672e0_dp * t55 * t1689 - 0.3333333336e0_dp * t83  &
            * t1694 * t804 + 0.1666666668e0_dp * t1546 * t897 * t3314
            t3827 = t649 * t804
            t3842 = t3278 * t46
            t3849 = t140 * t3382
            t3853 = t649 * t754
            t3861 = t78 * t649
            t3865 = t589 * t140
            t3869 = -0.1666666668e0_dp * t289 * t84 * t3827 - 0.1666666668e0_dp * t289  &
            * t922 * t593 + 0.1666666668e0_dp * t1546 * t2636 * t1410 - 0.6666666672e0_dp  &
            * t55 * t1703 - 0.3333333336e0_dp * t83 * t597 * t1067 - 0.1250000001e0_dp  &
            * t3842 * t897 * t3306 - 0.3333333336e0_dp * t83 * t384 *  &
            t3382 - 0.1666666668e0_dp * t289 * t84 * t3849 + 0.8333333340e-1_dp * t906  &
            * t84 * t3853 - 0.6666666672e0_dp * t55 * t1568 - 0.6666666672e0_dp *  &
            t55 * t1565 - 0.3333333336e0_dp * t165 * t3861 * t237 - 0.3333333336e0_dp  &
            * t165 * t3865 * t237
            t3871 = t379 * t1447
            t3875 = t45 * t2002
            t3892 = t225 * t1447
            t3918 = t25 * t3382
            t3925 = t777 * t561
            t3933 = t451 * t1410
            t3939 = t1301 * t245
            t3940 = t172 * t3306
            t3965 = 0.1500000000e1_dp * t405 * t156 * t3918 + 0.1968750000e2_dp * t2823  &
            * t401 * t3306 + 0.1666666666e1_dp * t405 * t1340 * t3925 + 0.3000000000e1_dp  &
            * t405 * t1333 * t1447 + 0.3000000000e1_dp * t399 * t822 * t3933  &
            + 0.2250000000e1_dp * t400 * t401 * t3382 + 0.5625000000e1_dp * t3939  &
            * t3940 + 0.3333333333e0_dp * t405 * t194 * t649 * t14 - 0.1125000000e2_dp  &
            * t1302 * t1309 * t1410 + 0.2250000000e1_dp * t400 * t1892 * t804  &
            - 0.5625000000e1_dp * t1302 * t401 * t3329 - 0.5625000000e1_dp * t1302  &
            * t1892 * t754 + 0.2250000000e1_dp * t400 * t2808 * t561 - 0.4500000000e1_dp  &
            * t400 * t1333 * t1410
            t3969 = t172 * t3314
            t3975 = t19 * t1688
            t3979 = t290 * t1483
            t3985 = t172 * t3329
            t4010 = -0.1125000000e2_dp * t1302 * t401 * t3314 - 0.4500000000e1_dp * t1912  &
            * t3969 + 0.4500000000e1_dp * t400 * t1309 * t1447 + 0.4500000000e1_dp  &
            * t400 * t3975 * t237 - 0.3750000000e1_dp * t3939 * t3979 - 0.2000000000e1_dp  &
            * t405 * t417 * t3892 - 0.2250000000e1_dp * t1912 * t3985 -  &
            0.2000000000e1_dp * t405 * t2873 * t561 - 0.5000000000e0_dp * t400 * t1316  &
            * t3472 - 0.1000000000e1_dp * t405 * t68 * t1688 * t17 + 0.1500000000e1_dp  &
            * t400 * t406 * t3449 - 0.1500000000e1_dp * t151 * t50 * t3796  &
            + 0.1500000000e1_dp * t400 * t1902 * t833 + 0.1500000000e1_dp * t400 *  &
            t1320 * t1483
            t4018 = 0.2222222222e0_dp * t2885 * t446 * t1172 * t561 - 0.1666666667e0_dp  &
            * t1246 * t1292 * t1865 - 0.4444444445e0_dp * t1373 * t446 * t3892  &
            - 0.1111111111e0_dp * t165 * t171 * t1688 * t17 + 0.1666666667e0_dp * t165  &
            * t56 * t3576 * t140 + 0.1250000000e0_dp * t2763 * t446 * t1168 *  &
            t561 + 0.3333333334e0_dp * t165 * t689 * t413 + 0.1666666667e0_dp * t165  &
            * t679 * t1067 - 0.1666666667e0_dp * t1246 * t446 * t413 * t1447 +  &
            (2 * t147 * t149 * (t3965 + t4010)) - 0.4444444445e0_dp * t1373 *  &
            t1292 * t1829
            t4025 = t3262 * t57
            t4032 = t74 * t649
            t4070 = t3525 * t67
            t4079 = -0.5555555558e-1_dp * t192 * t199 * t3514 + (t705 * t863 *  &
            t65) + (2 * t496 * t1504 * t65) + (t187 * t3513 * t65) -  &
            0.5555555558e-1_dp * t192 * t728 * t864 - 0.5555555556e0_dp * t140 * t3546  &
            * t1722 + (t705 * t865 * t65) - 0.1111111112e0_dp * t1106 * t1525  &
            - 0.1111111112e0_dp * t4070 * t196 + 0.6666666668e0_dp * t1748 * t266  &
            * t1722 + 0.3703703706e0_dp * t1193 * t1505 * t451
            t4101 = t498 * t1505
            t4108 = -0.5000000001e0_dp * t60 * t1189 * t709 + 0.3703703706e0_dp * t1755  &
            * t512 - 0.5555555558e-1_dp * t503 * t3561 - 0.2222222224e0_dp * t516  &
            * t267 * t172 * t1504 - 0.5555555558e-1_dp * t1715 * t886 - 0.5555555558e-1_dp  &
            * t1715 * t882 - 0.1111111112e0_dp * t516 * t3514 * t172 - 0.5000000001e0_dp  &
            * t1067 * t177 * t709 + (2 * t4101) - 0.1111111112e0_dp  &
            * t508 * t1769 - 0.1111111112e0_dp * t192 * t1837 * t267
            t4140 = -0.5555555558e-1_dp * t719 * t1207 - 0.2000000001e1_dp * t1763 *  &
            t266 * t1722 - 0.1111111112e0_dp * t192 * t530 * t1505 + 0.1333333334e1_dp  &
            * t1519 * t722 * t1712 + 0.7407407410e-1_dp * t516 * t1718 * t526  &
            + 0.6666666668e0_dp * t379 * t1518 * t1722 - 0.5555555558e-1_dp * t192 *  &
            t728 * t866 - 0.2222222224e0_dp * t719 * t1126 + 0.1800000000e2_dp * t1111  &
            * t178 * t710 - 0.1111111112e0_dp * t508 * t1782 - 0.5555555558e-1_dp  &
            * t58 * t3796 * t71
            t4156 = t865 * t561 * t65
            t4201 = t20 * t1688
            t4214 = -0.8000000000e1_dp * t478 * t123 * t3925 - 0.2400000000e2_dp * t473  &
            * t479 * t3449 + 0.1200000000e2_dp * t478 * t489 * t3892 + 0.1200000000e2_dp  &
            * t478 * t2998 * t561 - 0.2000000000e1_dp * t478 * t27 * t649  &
            * t14 + 0.1200000000e2_dp * t478 * t56 * t1688 * t17 - 0.2400000000e2_dp  &
            * t473 * t1161 * t1483 + 0.1080000000e3_dp * t1143 * t474 * t3314 -  &
            0.2160000000e3_dp * t2994 * t474 * t3306 + 0.7200000000e2_dp * t1821 *  &
            t3969 - 0.3600000000e2_dp * t473 * t4201 * t237 - 0.1800000000e2_dp * t473  &
            * t3030 * t561 - 0.2400000000e2_dp * t473 * t1811 * t833 - 0.1800000000e2_dp  &
            * t478 * t183 * t3918
            t4230 = t1142 * t262
            t4255 = 0.5400000000e2_dp * t1143 * t474 * t3329 + 0.5400000000e2_dp * t1143  &
            * t1801 * t754 + 0.7200000000e2_dp * t473 * t1176 * t1410 + 0.8999999998e1_dp  &
            * t178 * t62 * t3796 - 0.2400000000e2_dp * t472 * t859 * t3933  &
            - 0.1080000000e3_dp * t4230 * t3940 + 0.4000000000e1_dp * t473 * t1157  &
            * t3472 + 0.1080000000e3_dp * t1143 * t1150 * t1410 + 0.3600000000e2_dp  &
            * t1821 * t3985 + 0.7200000000e2_dp * t4230 * t3979 - 0.1800000000e2_dp  &
            * t473 * t1801 * t804 - 0.3600000000e2_dp * t473 * t1150 * t1447 -  &
            0.3600000000e2_dp * t478 * t1176 * t1447 - 0.1800000000e2_dp * t473 * t474  &
            * t3382
            t4256 = t4214 + t4255
            t4261 = -0.2222222224e0_dp * t1510 * t515 * t517 + 0.3703703706e0_dp * t719  &
            * t1213 - 0.1000000001e1_dp * t1766 * t266 * t1722 - 0.5555555558e-1_dp  &
            * t4070 * t200 + 0.6666666668e0_dp * t1519 * t496 * t561 * t65 - 0.5000000001e0_dp  &
            * t732 * t4156 - 0.5000000001e0_dp * t140 * t865 * t1722  &
            - 0.5000000001e0_dp * t804 * t187 * t1722 - 0.4938271608e0_dp * t719 *  &
            t1101 + 0.8999999998e1_dp * t1097 * t178 * t710 - 0.5555555558e-1_dp *  &
            t192 * t68 * t4256 * t65
            t4272 = t863 * t561 * t65
            t4292 = (2 * t1836 * t266 * t65) - 0.1111111112e0_dp * t719 * t1224  &
            - 0.2222222224e0_dp * t1137 * t1790 - 0.1000000001e1_dp * t723 * t4272  &
            - 0.1111111112e0_dp * t1785 * t1134 - 0.6172839508e-1_dp * t192 * t22 *  &
            t705 * t1200 - 0.5000000001e0_dp * t140 * t863 * t1722 - 0.1111111112e0_dp  &
            * t1755 * t531 - 0.6172839508e-1_dp * t1715 * t874 - 0.5555555558e-1_dp  &
            * t719 * t1122 + 0.7407407410e-1_dp * t508 * t1726
            t4317 = 0.7407407410e-1_dp * t1755 * t527 - 0.6172839508e-1_dp * t719 * t1201  &
            - 0.1111111112e0_dp * t508 * t1838 - 0.2222222224e0_dp * t192 * t3552  &
            * t1125 + 0.7407407410e-1_dp * t1785 * t1221 - 0.1111111111e1_dp * t3547  &
            * t1737 + 0.6666666668e0_dp * t1519 * t187 * t1712 + 0.7407407410e-1_dp  &
            * t503 * t3553 - 0.1111111112e0_dp * t1755 * t534 - 0.1000000001e1_dp  &
            * t1774 * t1712 + 0.7407407410e-1_dp * t719 * t1119
            t4345 = -0.1000000001e1_dp * t379 * t266 * t1722 - 0.1000000001e1_dp * t804  &
            * t722 * t1722 + 0.8999999998e1_dp * t1109 * t178 * t710 - 0.5555555558e-1_dp  &
            * t1129 * t729 + t4256 * t65 - 0.2222222224e0_dp * t1755 * t520  &
            - 0.1111111112e0_dp * t503 * t3569 - 0.1111111112e0_dp * t1785 * t1216  &
            - 0.1111111112e0_dp * t3097 * t4101 - 0.5000000001e0_dp * t732 * t4272  &
            + 0.7407407410e-1_dp * t192 * t194 * t1836 * t526
            t4352 = t1731 * t67
            t4378 = 0.7407407410e-1_dp * t1715 * t878 + 0.8999999998e1_dp * t1189 * t59  &
            * t177 * t710 + 0.7407407410e-1_dp * t4352 * t271 - 0.1111111112e0_dp  &
            * t4352 * t275 - 0.5555555556e0_dp * t3547 * t1745 + 0.7407407410e-1_dp  &
            * t516 * t1709 * t526 + 0.6666666668e0_dp * t3533 * t187 * t1722 - 0.1000000001e1_dp  &
            * t237 * t496 * t1722 + 0.3703703706e0_dp * t579 * t1192  &
            * t1194 - 0.1000000001e1_dp * t723 * t4156 + 0.1333333334e1_dp * t3533  &
            * t722 * t1722 - 0.5555555558e-1_dp * t719 * t1235
            t4392 = 0.3333333334e0_dp * t165 * t171 * t1825 + 0.1250000000e0_dp * t1981  &
            * t1247 * t3329 - 0.3125000000e0_dp * t4025 * t1247 * t3306 + 0.3333333334e0_dp  &
            * t165 * t1392 * t1447 + 0.5555555556e-1_dp * t1246 * t4032  &
            * t833 - 0.1666666667e0_dp * t1246 * t74 * t25 * t1410 + 0.5555555556e-1_dp  &
            * t1246 * t1247 * t3449 + 0.1800000000e2_dp * t657 * t1283 * t59  &
            * t1975 - 0.1666666667e0_dp * t429 * t166 * t3871 + 0.3333333334e0_dp  &
            * t55 * t57 * (t4079 + t4108 + t4140 + t4261 + t4292 + t4317 + t4345  &
            + t4378) + 0.1481481481e0_dp * t170 * t114 * t735 * t14 - 0.1666666667e0_dp  &
            * t429 * t1928 * t290
            t4428 = -0.1111111111e0_dp * t165 * t1264 * t1483 + 0.1800000000e2_dp * t658  &
            * t575 * t1348 - 0.8333333335e-1_dp * t429 * t679 * t918 - 0.1111111111e0_dp  &
            * t165 * t455 * t1988 + 0.3333333334e0_dp * t170 * t27 * t3576  &
            * t172 - 0.8333333335e-1_dp * t429 * t1966 * t754 + 0.1250000000e0_dp  &
            * t1368 * t679 * t907 + 0.1666666667e0_dp * t445 * t446 * t3918 - 0.4444444445e0_dp  &
            * t165 * t450 * t1829 + 0.2500000000e0_dp * t1981 * t3244  &
            * t1410 - 0.1666666667e0_dp * t429 * t1392 * t1410
            t4466 = -0.1666666667e0_dp * t429 * t437 * t1553 - 0.1111111111e0_dp * t165  &
            * t689 * t1085 + 0.1666666667e0_dp * t165 * t1275 * t649 - 0.2222222222e0_dp  &
            * t170 * t27 * t1843 * t17 - 0.8333333335e-1_dp * t1246 * t1962  &
            * t1168 + 0.1666666667e0_dp * t165 * t2720 * t561 + 0.6666666668e0_dp  &
            * t170 * t1932 * t25 + 0.3333333334e0_dp * t445 * t1292 * t1825 + 0.1250000000e0_dp  &
            * t1368 * t166 * t3853 - 0.4444444445e0_dp * t1373 * t1962  &
            * t1172 + 0.7407407405e-1_dp * t165 * t1377 * t3472 + 0.5555555555e-1_dp  &
            * t1246 * t3244 * t1483
            t4498 = t56 * t1843
            t4505 = 0.1666666667e0_dp * t165 * t166 * t3796 - 0.8333333334e-1_dp * t2763  &
            * t1247 * t833 * t561 + (4 * t147 * t392 * t1886) + (2  &
            * t147 * t1388 * t672) - 0.1111111111e0_dp * t165 * t1940 * t833 + 0.3333333334e0_dp  &
            * t165 * t1928 * t379 + 0.7407407409e-1_dp * t165 * t450  &
            * t649 * t14 - 0.1666666667e0_dp * t429 * t462 * t3314 - 0.3703703704e-1_dp  &
            * t2885 * t1247 * t3472 + 0.3333333334e0_dp * t165 * t4498 * t237  &
            + (2 * t147 * t1284 * t672)
            t4512 = t1531 * rho
            t4545 = (4 * t147 * t1885 * t424) + 0.1666666667e0_dp * t165 * t462  &
            * t3382 + 0.3333333334e0_dp * t445 * t4512 * t413 + 0.5185185185e0_dp  &
            * t2937 * t446 * t3925 - 0.8333333335e-1_dp * t429 * t166 * t3827 - 0.1666666667e0_dp  &
            * t429 * t437 * t1571 + 0.5555555555e-1_dp * t1246 * t203  &
            * t237 * t1483 + 0.5555555555e-1_dp * t1246 * t2707 * t1483 + 0.7407407409e-1_dp  &
            * t165 * t1875 * t1360 + (4 * t147 * t393 * t1924)  &
            + 0.3600000000e2_dp * (t1971) * (t1972) * (t1973) * (t424)  &
            + 0.3333333334e0_dp * t165 * t437 * t1688
            t4553 = t585 * t140
            t4575 = t1067 * t561
            t4585 = 0.1666666667e0_dp * t445 * t3248 * t668 - 0.1111111111e0_dp * t165  &
            * t1932 * t433 + 0.5555555556e-1_dp * t1246 * t4553 * t833 - 0.8333333335e-1_dp  &
            * t429 * t1275 * t593 - 0.8333333335e-1_dp * t429 * t462 *  &
            t3329 + 0.2500000000e0_dp * t1981 * t1247 * t3314 + 0.1666666667e0_dp *  &
            t445 * t1962 * t1093 + 0.3333333334e0_dp * t165 * t455 * t668 - 0.8333333335e-1_dp  &
            * t429 * t166 * t3849 - 0.8333333335e-1_dp * t429 * t166 *  &
            t4575 + 0.2e1_dp * t147 * (t3491 - t3494 - t3497 + t3501 - t3505 + 0.9000000000e1_dp  &
            * t3507 - t3512) * t149 * t160
            t4607 = t1688 * t237
            t4627 = 0.1666666667e0_dp * t165 * t1966 * t804 - 0.1666666667e0_dp * t429  &
            * t437 * t1557 + 0.1800000000e2_dp * t657 * t1387 * t59 * t1975 + 0.1037037037e1_dp  &
            * t170 * t341 * t585 * t1100 - 0.8888888890e0_dp * t170  &
            * t114 * t1531 * t451 - 0.8888888890e0_dp * t170 * t1875 * t225 - 0.1666666667e0_dp  &
            * t429 * t166 * t4607 - 0.1111111111e0_dp * t165 * t466  &
            * t3449 - 0.1666666667e0_dp * t429 * t679 * t914 - 0.8333333335e-1_dp *  &
            t1246 * t446 * t1093 * t561 + 0.1250000000e0_dp * t1368 * t462 * t3306  &
            + 0.2500000000e0_dp * t1981 * t2707 * t1410
            t4657 = -0.3333333336e0_dp * t289 * t84 * t3871 - 0.6666666672e0_dp * t83  &
            * t3875 * t237 + 0.8333333340e-1_dp * t906 * t384 * t3306 - 0.1666666668e0_dp  &
            * t289 * t597 * t918 - 0.2222222224e0_dp * t144 * t46 * (t4018  &
            + t4392 + t4428 + t4466 + t4505 + t4545 + t4585 + t4627) - 0.1666666668e0_dp  &
            * t289 * t384 * t3329 - 0.3333333336e0_dp * t165 * t2636 * t561  &
            - 0.3333333336e0_dp * t165 * t897 * t1447 + 0.8333333340e-1_dp * t1546  &
            * t897 * t3329 - 0.6666666672e0_dp * t83 * t1071 * t1447 - 0.3333333336e0_dp  &
            * t83 * t2404 * t561 - 0.4444444448e0_dp * t383 * t3875 - 0.3333333336e0_dp  &
            * t289 * t294 * t1553
            t4696 = -0.6666666672e0_dp * t55 * t1575 - 0.3333333336e0_dp * t289 * t597  &
            * t914 - 0.3333333336e0_dp * t165 * t2631 * t561 - 0.6666666672e0_dp *  &
            t55 * t1706 - 0.3333333336e0_dp * t289 * t384 * t3314 - 0.6666666672e0_dp  &
            * t83 * t294 * t1688 + 0.1666666668e0_dp * t1546 * t2631 * t1410 -  &
            0.3333333336e0_dp * t289 * t1564 * t290 - 0.1666666668e0_dp * t289 * t84  &
            * t4575 - 0.3333333336e0_dp * t289 * t294 * t1557 - 0.1666666668e0_dp  &
            * t289 * t1694 * t754 - 0.3333333336e0_dp * t289 * t84 * t4607 - 0.3333333336e0_dp  &
            * t83 * t45 * t3581 * t140
            t4700 = -0.6666666670e0_dp * t2147 * t1411 - 0.1333333334e1_dp * t747 * t597  &
            * t237 - 0.1333333334e1_dp * t747 * t84 * t1447 - 0.1333333334e1_dp  &
            * t747 * t294 * t561 - 0.8888888896e0_dp * t214 * t1536 + 0.1666666668e0_dp  &
            * t3305 * t84 * t3306 - 0.3333333336e0_dp * t753 * t46 * t589 * t754  &
            - 0.6666666672e0_dp * t1409 * t84 * t3314 - 0.6666666670e0_dp * t1409  &
            * t294 * t1410 - 0.1333333334e1_dp * t217 * t46 * t1535 * t237 - 0.1333333334e1_dp  &
            * t217 * t46 * t282 * t1447 - 0.3333333336e0_dp * t1409  &
            * t84 * t3329 - 0.6666666672e0_dp * t217 * t46 * t589 * t804 - 0.6666666672e0_dp  &
            * t217 * t46 * t78 * t3382 - 0.6666666672e0_dp * t217 * t46  &
            * t893 * t561 - 0.4444444448e0_dp * t43 * t46 * t3581 + t8 * (t3826 +  &
            t3869 + t4657 + t4696)
            e_ndrho_rho_rho_spin(ii) = e_ndrho_rho_rho_spin(ii) +  t4700 * sx
            t4701 = t561 ** 2
            t4704 = 0.3333333336e0_dp * t2147 * t84 * t4701
            t4707 = 0.1333333334e1_dp * t747 * t597 * t561
            t4708 = t22 * t36
            t4713 = t557 ** 2
            t4714 = t327 * t4713
            t4717 = t22 * t132
            t4721 = 0.1e1_dp / t787 * t2480
            t4725 = 0.120e2_dp * t624 * t4717 - 0.60e1_dp * t7 * t4721 * t371
            t4726 = t105 * t4725
            t4729 = 0.40e1_dp * t624 * t4708 - 0.80e1_dp * t112 * t558 + 0.40e1_dp * t102  &
            * t4714 - 0.20e1_dp * t102 * t4726
            t4732 = 0.6666666672e0_dp * t747 * t84 * t4729
            t4733 = pi ** 2
            t4734 = t4733 * t146
            t4736 = 0.1e1_dp / t58 / omega
            t4738 = t4734 * t4736 * t1300
            t4739 = t4701 * t149
            t4743 = t19 * t4701
            t4744 = t399 * t4743
            t4748 = t151 * t19 * t4729
            t4760 = t20 * t4701
            t4764 = t20 * t4729
            t4768 = t58 ** 2
            t4769 = 0.1e1_dp / t4768
            t4770 = t4769 * t1141
            t4771 = t4770 * t4733
            t4772 = t26 * t4701
            t4776 = t58 * t4729
            t4782 = -0.2250000000e2_dp * t473 * t4760 * t65 + 0.8999999998e1_dp * t478  &
            * t4764 * t65 + 0.8099999996e2_dp * t4771 * t4772 * t65 - 0.5555555558e-1_dp  &
            * t4776 * t71 - 0.5000000001e0_dp * t60 * t4729 * t65
            t4786 = -0.2700000000e2_dp * t4738 * t18 * t4739 + 0.4500000000e1_dp * t567  &
            * t4744 - 0.3000000000e1_dp * t567 * t4748 - 0.8333333335e-1_dp * t429  &
            * t166 * t4701 + 0.3333333334e0_dp * t165 * t679 * t561 + 0.1666666667e0_dp  &
            * t165 * t166 * t4729 + 0.3333333334e0_dp * t55 * t57 * t4782
            t4789 = 0.4444444448e0_dp * t214 * t46 * t4786
            t4790 = t237 * t4701
            t4793 = 0.1666666668e0_dp * t3305 * t84 * t4790
            t4796 = 0.6666666672e0_dp * t1409 * t597 * t1410
            t4797 = t1447 * t561
            t4800 = 0.6666666672e0_dp * t1409 * t84 * t4797
            t4801 = t237 * t4729
            t4804 = 0.3333333336e0_dp * t1409 * t84 * t4801
            t4808 = 0.6666666672e0_dp * t217 * t46 * t4786 * t237
            t4812 = 0.1333333334e1_dp * t217 * t46 * t589 * t1447
            t4817 = t953 * t4713
            t4820 = t316 * t4725
            t4829 = t233 * t4713
            t4830 = t974 * t4829
            t4833 = t1443 * t557
            t4834 = t327 * t4833
            t4837 = t233 * t4725
            t4838 = t327 * t4837
            t4849 = -0.3200000000e2_dp * t1645 * t627 + 0.4800000000e2_dp * t7 * t1440  &
            - 0.2400000000e2_dp * t786 / t2492 * t1058
            t4850 = t105 * t4849
            t4853 = -0.1066666667e2_dp * t1645 * t98 + 0.2133333334e2_dp * t339 * t1420  &
            - 0.1066666667e2_dp * t321 * t4817 + 0.5333333333e1_dp * t321 * t4820  &
            - 0.40e1_dp * t624 * t234 + 0.160e2_dp * t112 * t1426 - 0.80e1_dp * t112 *  &
            t1444 - 0.120e2_dp * t102 * t4830 + 0.80e1_dp * t102 * t4834 + 0.40e1_dp *  &
            t102 * t4838 - 0.20e1_dp * t102 * t4850
            t4857 = 0.6666666672e0_dp * t217 * t46 * t78 * t4853
            t4861 = 0.3333333336e0_dp * t753 * t46 * t282 * t4701
            t4864 = t217 * t46 * t1535 * t561
            t4869 = 0.6666666672e0_dp * t217 * t46 * t282 * t4729
            t4872 = 0.8333333335e-1_dp * t429 * t437 * t4701
            t4874 = t165 * t1928 * t561
            t4878 = 0.1666666667e0_dp * t165 * t437 * t4729
            t4881 = 0.5400000000e2_dp * t1143 * t1493 * t4701
            t4884 = 0.3600000000e2_dp * t473 * t2035 * t561
            t4887 = 0.1800000000e2_dp * t473 * t1493 * t4729
            t4889 = t178 * t62 * t4853
            t4893 = 0.1200000000e2_dp * t473 * t1500 * t4701
            t4895 = t478 * t1500 * t4729
            t4897 = t4881 - t4884 - t4887 + 0.8999999998e1_dp * t4889 - t4893 + 0.5999999999e1_dp  &
            * t4895
            t4898 = t4897 * t65
            t4901 = t1504 * t59 * t177 * t710
            t4904 = t4701 * t65
            t4905 = t62 * t4904
            t4907 = 0.2250000000e2_dp * t1506 * t471 * t4905
            t4908 = t4729 * t65
            t4909 = t62 * t4908
            t4911 = 0.8999999998e1_dp * t1507 * t4909
            t4915 = t4733 * t26 * t4904
            t4917 = 0.8099999996e2_dp * t266 * t4769 * t1141 * t4915
            t4918 = t58 * t4853
            t4920 = 0.5555555558e-1_dp * t4918 * t71
            t4922 = t1447 * t177 * t709
            t4925 = 0.5000000001e0_dp * t1513 * t4908
            t4929 = t60 * t1504 * t709
            t4932 = 0.5000000001e0_dp * t1528 * t4908
            t4935 = 0.1000000000e1_dp * t237 * t471 * t4904
            t4940 = t59 * pi * t20 * t65
            t4941 = t237 * t1141 * t4701 * t4940
            t4943 = t4776 * t67
            t4945 = 0.3703703705e-1_dp * t4943 * t271
            t4946 = t471 * t1518
            t4951 = 0.5555555558e-1_dp * t4943 * t275
            t4952 = t719 * t1525
            t4955 = t68 * t4897 * t65
            t4958 = t4898 + 0.1800000000e2_dp * t4901 - t4907 + t4911 + t4917 - t4920  &
            - 0.1000000000e1_dp * t4922 - t4925 + 0.3333333334e0_dp * t1519 * t4908  &
            - 0.1000000000e1_dp * t4929 - t4932 + t4935 - 0.4500000000e1_dp * t4941  &
            + t4945 + 0.3000000000e1_dp * t4946 * t4701 * t4940 - t4951 - 0.1111111112e0_dp  &
            * t4952 - 0.5555555558e-1_dp * t192 * t4955
            t4994 = t44 * t4733 * t146
            t4995 = t4994 * t4770
            t5001 = 0.1666666667e0_dp * t165 * t166 * t4853
            t5002 = t17 * t4701
            t5005 = 0.2777777778e-1_dp * t429 * t171 * t5002
            t5008 = 0.1111111111e0_dp * t165 * t689 * t1483
            t5009 = t17 * t4729
            t5011 = t165 * t171 * t5009
            t5013 = t27 * t4782
            t5016 = 0.1111111111e0_dp * t170 * t5013 * t17
            t5019 = 0.1250000000e0_dp * t1368 * t166 * t4790
            t5022 = 0.1666666667e0_dp * t429 * t679 * t1410
            t5025 = 0.1666666667e0_dp * t429 * t166 * t4797
            t5028 = 0.8333333335e-1_dp * t429 * t166 * t4801
            t5029 = t56 * t4782
            t5032 = 0.1666666667e0_dp * t165 * t5029 * t237
            t5035 = 0.3333333334e0_dp * t165 * t679 * t1447
            t5036 = -t4872 + 0.3333333334e0_dp * t4874 + t4878 + 0.3333333334e0_dp *  &
            t55 * t57 * t4958 + 0.2e1_dp * t147 * t149 * (-0.5625000000e1_dp * t1302  &
            * t1459 * t4701 + 0.4500000000e1_dp * t400 * t2022 * t561 + 0.2250000000e1_dp  &
            * t400 * t1459 * t4729 - 0.1500000000e1_dp * t151 * t50 * t4853  &
            + 0.7500000000e0_dp * t400 * t1466 * t4701 - 0.5000000000e0_dp * t405  &
            * t1466 * t4729) - 0.3600000000e2_dp * t2081 * t4760 * t250 + 0.1800000000e2_dp  &
            * t658 * t4764 * t250 + 0.3600000000e2_dp * t658 * t575 * t1471  &
            + 0.1620000000e3_dp * t4995 * t4772 * t250 + t5001 + t5005 - t5008  &
            - 0.5555555555e-1_dp * t5011 - t5016 + t5019 - t5022 - t5025 - t5028  &
            + t5032 + t5035
            t5040 = t379 * t4729
            t5043 = 0.1666666668e0_dp * t289 * t84 * t5040
            t5044 = t140 * t4853
            t5047 = 0.1666666668e0_dp * t289 * t84 * t5044
            t5050 = 0.1666666668e0_dp * t165 * t897 * t4729
            t5051 = beta * t86
            t5060 = beta * t92
            t5065 = rho * t4713
            t5069 = rho * t4725
            t5079 = t136 * t4713
            t5083 = t645 * t557
            t5087 = t136 * t4725
            t5106 = beta * t640
            t5108 = t126 * t371 * ndrho
            t5115 = t366 * t640
            t5118 = -0.480e2_dp * t5051 * t4717 + 0.240e2_dp * t7 * t2182 * t2480 * t371  &
            + 0.480e2_dp * t5060 * t628 - 0.240e2_dp * t7 * t1640 * t1438 * t371  &
            * rho + 0.120e2_dp * t550 * t637 - 0.180e2_dp * t5106 * t5108 - 0.120e2_dp  &
            * t7 * t640 * t635 * t371 + 0.180e2_dp * t5115 * t1681
            t5122 = -0.320e2_dp * t5051 * t4708 + 0.640e2_dp * t1582 * t558 - 0.320e2_dp  &
            * t302 * t4714 + 0.160e2_dp * t302 * t4726 + 0.320e2_dp * t5060 * t99 -  &
            0.640e2_dp * t605 * t609 + 0.320e2_dp * t93 * t953 * t5065 - 0.160e2_dp *  &
            t93 * t316 * t5069 - 0.40e1_dp * t624 * t137 + 0.160e2_dp * t112 * t615  &
            - 0.80e1_dp * t112 * t646 - 0.120e2_dp * t102 * t974 * t5079 + 0.80e1_dp *  &
            t102 * t327 * t5083 + 0.40e1_dp * t102 * t327 * t5087 - 0.20e1_dp * t102  &
            * t105 * t5118
            t5123 = t5122 * t237
            t5126 = 0.1666666668e0_dp * t289 * t84 * t5123
            t5139 = t25 * t4701
            t5143 = t25 * t4729
            t5147 = -0.5625000000e1_dp * t1302 * t401 * t4701 + 0.4500000000e1_dp * t400  &
            * t1892 * t561 + 0.2250000000e1_dp * t400 * t401 * t4729 - 0.1500000000e1_dp  &
            * t151 * t50 * t5122 - 0.2250000000e1_dp * t400 * t156 * t5139  &
            + 0.1500000000e1_dp * t405 * t156 * t5143
            t5163 = t140 * t4701
            t5170 = t649 * t561
            t5174 = t140 * t4729
            t5226 = 0.5400000000e2_dp * t1143 * t474 * t4701 - 0.3600000000e2_dp * t473  &
            * t1801 * t561 - 0.1800000000e2_dp * t473 * t474 * t4729 + 0.8999999998e1_dp  &
            * t178 * t62 * t5122 + 0.3600000000e2_dp * t473 * t183 * t5139  &
            - 0.1800000000e2_dp * t478 * t183 * t5143
            t5228 = t705 * t59
            t5229 = t5228 * t177
            t5237 = t649 * t177
            t5244 = t60 * t705
            t5247 = t187 * t4769
            t5251 = t58 * t5122
            t5256 = t140 * t1141
            t5262 = t471 * t722
            t5270 = t68 * t5226
            t5271 = t5270 * t65
            t5276 = t5226 * t65 + 0.1800000000e2_dp * t5229 * t710 - 0.2250000000e2_dp  &
            * t707 * t471 * t4905 + 0.8999999998e1_dp * t708 * t4909 - 0.1000000000e1_dp  &
            * t5237 * t709 - 0.5000000001e0_dp * t716 * t4908 - 0.1000000001e1_dp  &
            * t723 * t4908 - 0.1000000000e1_dp * t5244 * t709 + 0.8099999996e2_dp  &
            * t5247 * t1141 * t4915 - 0.5555555558e-1_dp * t5251 * t71 + 0.1000000000e1_dp  &
            * t2098 * t4904 - 0.4500000000e1_dp * t5256 * t4701 * t4940  &
            - 0.1111111112e0_dp * t4943 * t196 - 0.9000000007e1_dp * t5262 * t4701  &
            * t4940 - 0.5555555558e-1_dp * t4943 * t200 - 0.1111111112e0_dp * t719 *  &
            t729 - 0.5555555558e-1_dp * t192 * t5271 - 0.5000000001e0_dp * t732 * t4908
            t5280 = (2 * t147 * t149 * t5147) - 0.3600000000e2_dp * t2081 * t4760  &
            * t161 + 0.1800000000e2_dp * t658 * t4764 * t161 + 0.3600000000e2_dp  &
            * t658 * t575 * t673 + 0.1620000000e3_dp * t4995 * t4772 * t161 + 0.1250000000e0_dp  &
            * t1368 * t166 * t5163 - 0.1666666667e0_dp * t429 * t679  &
            * t593 - 0.1666666667e0_dp * t429 * t166 * t5170 - 0.8333333335e-1_dp *  &
            t429 * t166 * t5174 + 0.1666666667e0_dp * t165 * t5029 * t140 + 0.3333333334e0_dp  &
            * t165 * t679 * t649 + 0.1666666667e0_dp * t165 * t166 * t5122  &
            - 0.8333333335e-1_dp * t1246 * t446 * t5139 + 0.3333333334e0_dp * t445  &
            * t1962 * t668 + 0.1666666667e0_dp * t445 * t446 * t5143 + 0.3333333334e0_dp  &
            * t170 * t5013 * t172 - 0.8333333335e-1_dp * t429 * t462 * t4701  &
            + 0.3333333334e0_dp * t165 * t1966 * t561 + 0.1666666667e0_dp * t165  &
            * t462 * t4729 + 0.3333333334e0_dp * t55 * t57 * t5276
            t5281 = t45 * t5280
            t5284 = 0.3333333336e0_dp * t83 * t5281 * t237
            t5291 = 0.1666666668e0_dp * t289 * t384 * t4801
            t5294 = 0.1666666668e0_dp * t1546 * t3865 * t1410
            t5295 = t45 * t4786
            t5298 = 0.3333333336e0_dp * t83 * t5295 * t379
            t5301 = 0.3333333336e0_dp * t55 * t384 * t4729
            t5303 = t83 * t1564 * t649
            t5305 = -t5043 - t5047 - t5050 - t5126 - t5284 - 0.3333333336e0_dp * t83  &
            * t45 * t5036 * t140 - t5291 + t5294 - t5298 - t5301 - 0.6666666672e0_dp  &
            * t5303
            t5308 = 0.3333333336e0_dp * t289 * t1694 * t1410
            t5311 = 0.1666666668e0_dp * t289 * t1071 * t4701
            t5314 = 0.1666666668e0_dp * t289 * t5295 * t290
            t5317 = 0.3333333336e0_dp * t165 * t3861 * t561
            t5320 = 0.3333333336e0_dp * t83 * t1071 * t4729
            t5323 = 0.1666666668e0_dp * t1546 * t3861 * t1410
            t5326 = 0.1250000001e0_dp * t3842 * t897 * t4790
            t5329 = 0.6666666672e0_dp * t83 * t1694 * t1447
            t5330 = t379 * t4701
            t5333 = 0.8333333340e-1_dp * t906 * t84 * t5330
            t5335 = 0.2222222224e0_dp * t383 * t5281
            t5338 = 0.6666666672e0_dp * t55 * t1694 * t561
            t5339 = t649 * t1447
            t5342 = 0.3333333336e0_dp * t289 * t84 * t5339
            t5343 = -t5308 - t5311 - t5314 - t5317 - t5320 + t5323 - t5326 - t5329  &
            + t5333 - t5335 - t5338 - t5342
            t5345 = t1688 * t561
            t5347 = t289 * t84 * t5345
            t5351 = 0.1666666668e0_dp * t164 * t384 * t4701
            t5354 = 0.3333333336e0_dp * t83 * t384 * t4853
            t5357 = 0.6666666672e0_dp * t55 * t597 * t649
            t5360 = 0.3333333336e0_dp * t289 * t384 * t4797
            t5363 = 0.8333333335e-1_dp * t429 * t166 * t5040
            t5367 = 0.1666666667e0_dp * t445 * t4782 * rho * t413
            t5370 = 0.1666666667e0_dp * t429 * t679 * t1557
            t5371 = t225 * t4701
            t5373 = t2885 * t446 * t5371
            t5376 = t165 * t466 * t5009
            t5379 = t165 * t1928 * t649
            t5383 = 0.4166666667e-1_dp * t2763 * t1247 * t5002
            t5386 = 0.8333333335e-1_dp * t429 * t5029 * t290
            t5389 = 0.8333333335e-1_dp * t429 * t166 * t5044
            t5397 = 0.8333333335e-1_dp * t1246 * t446 * t413 * t4729
            t5400 = 0.3333333334e0_dp * t445 * t1962 * t1825
            t5402 = t165 * t679 * t1688
            t5406 = 0.8333333335e-1_dp * t429 * t166 * t5123
            t5409 = 0.8333333335e-1_dp * t429 * t1392 * t4701
            t5412 = 0.8333333335e-1_dp * t429 * t171 * t5139
            t5415 = 0.1666666667e0_dp * t429 * t462 * t4797
            t5418 = 2 * t147 * t393 * t5147
            t5421 = t657 * t1884 * t59 * t1975
            t5423 = -t5363 + t5367 - t5370 + 0.1111111111e0_dp * t5373 - 0.5555555555e-1_dp  &
            * t5376 + 0.3333333334e0_dp * t5379 - t5383 - t5386 - t5389 +  &
            0.3333333334e0_dp * t170 * t27 * t4958 * t172 - t5397 + t5400 + 0.3333333334e0_dp  &
            * t5402 - t5406 - t5409 - t5412 - t5415 + t5418 + 0.3600000000e2_dp  &
            * t5421
            t5426 = 0.8333333335e-1_dp * t1246 * t1292 * t5139
            t5429 = 0.1666666667e0_dp * t445 * t1292 * t5143
            t5430 = t225 * t4729
            t5432 = t1373 * t446 * t5430
            t5436 = 0.1666666667e0_dp * t165 * t5029 * t379
            t5440 = 0.40e1_dp * t102 * t327 * t375 * t4725
            t5443 = 0.960e2_dp * t1597 * t3608 * t4829
            t5446 = 0.480e2_dp * t24 * t23 * t1618
            t5449 = 0.1173333333e3_dp * t307 * t2611 * t5065
            t5451 = 0.640e2_dp * t605 * t1603
            t5454 = 0.2133333334e2_dp * t321 * t953 * t5083
            t5456 = 0.1280e3_dp * t1582 * t1426
            t5459 = 0.5333333333e1_dp * t321 * t316 * t5118
            t5460 = beta * t306
            t5462 = 0.1173333333e3_dp * t5460 * t313
            t5466 = 0.40e1_dp * t102 * t327 * t136 * t4849
            t5468 = 0.2346666666e3_dp * t605 * t1420
            t5469 = t5440 - t5443 - t5446 - t5449 - t5451 - t5454 - t5456 + t5459  &
            - t5462 + t5466 - t5468
            t5472 = 0.3200000001e2_dp * t321 * t2455 * t5079
            t5474 = 0.4266666668e2_dp * t339 * t1608
            t5478 = 0.240e2_dp * t1616 * t973 * t645 * t1425
            t5480 = 0.1760000000e3_dp * t5060 * t627
            t5482 = 0.2160000000e3_dp * t7 * t1642
            t5484 = 0.1760000000e3_dp * t5460 * t347
            t5487 = 0.2160000000e3_dp * t7 * t995 * t1652
            t5489 = 0.1600000000e2_dp * t1429 * t1658
            t5492 = 0.7200000000e2_dp * beta * t368 * t5108
            t5494 = 0.6400000000e2_dp * t7 * t1677
            t5496 = t366 * t368 * t1681
            t5499 = 0.120e2_dp * t550 * t1669
            t5503 = 0.180e2_dp * t5106 * t360 * t371 * ndrho
            t5507 = 0.120e2_dp * t7 * t640 * t1667 * t371
            t5509 = 0.180e2_dp * t5115 * t3668
            t5511 = t1680 * t3673 * t23
            t5516 = t5480 - t5482 - t5484 + t5487 - t5489 + t5492 + t5494 - 0.1920000000e3_dp  &
            * t5496 + t5499 - t5503 - t5507 + t5509 - 0.4800000000e2_dp  &
            * t5511 + 0.1200000000e3_dp * t1055 * t2565 * t4721
            t5523 = 0.120e2_dp * t102 * t974 * t375 * t4713
            t5525 = 0.2133333334e2_dp * t339 * t1611
            t5528 = t102 * t327 * t1684 * t557
            t5531 = 0.160e2_dp * t112 * t1626
            t5533 = 0.160e2_dp * t112 * t1632
            t5535 = 0.640e2_dp * t1582 * t1444
            t5537 = 0.160e2_dp * t112 * t1622
            t5540 = 0.320e2_dp * t1597 * t1598 * t4837
            t5541 = t5472 - t5474 - t5478 - 0.20e1_dp * t102 * t105 * t5516 - t5523  &
            + t5525 + 0.80e1_dp * t5528 + t5531 + t5533 + t5535 + t5537 + t5540
            t5545 = 0.480e2_dp * t1616 * t3750 * t4829
            t5548 = 0.120e2_dp * t1616 * t1617 * t4837
            t5551 = 0.1280e3_dp * t24 * t298 * t1599
            t5553 = 0.160e2_dp * t302 * t4850
            t5558 = 0.320e2_dp * t5060 * t97 * t104 * rho * t233
            t5560 = 0.2346666666e3_dp * t1589 * t1592
            t5562 = 0.5866666667e2_dp * t93 * t4820
            t5564 = 0.1173333333e3_dp * t5060 * t98
            t5565 = t112 * t1685
            t5568 = 0.1066666667e2_dp * t1645 * t322
            t5572 = 0.40e1_dp * t102 * t327 * t5118 * t233
            t5574 = 0.1173333333e3_dp * t93 * t4817
            t5575 = t5545 - t5548 + t5551 + t5553 - t5558 + t5560 - t5562 + t5564  &
            - 0.80e1_dp * t5565 + t5568 + t5572 + t5574
            t5579 = 0.80e1_dp * t102 * t327 * t645 * t1443
            t5581 = 0.320e2_dp * t302 * t4838
            t5585 = 0.160e2_dp * t93 * t316 * rho * t4849
            t5590 = 0.80e1_dp * t624 * t22 * t326 * t136 * t233
            t5593 = 0.1066666667e2_dp * t321 * t953 * t5087
            t5595 = 0.640e2_dp * t302 * t4834
            t5597 = 0.960e2_dp * t302 * t4830
            t5600 = 0.240e2_dp * t1616 * t1617 * t4833
            t5602 = 0.40e1_dp * t624 * t376
            t5604 = 0.320e2_dp * t5051 * t234
            t5607 = 0.5866666667e2_dp * t307 * t949 * t5069
            t5610 = 0.640e2_dp * t1597 * t1598 * t4833
            t5611 = t5579 - t5581 - t5585 + t5590 - t5593 - t5595 + t5597 - t5600  &
            - t5602 + t5604 + t5607 + t5610
            t5613 = t5469 + t5541 + t5575 + t5611
            t5619 = 0.3333333334e0_dp * t170 * t5013 * t25
            t5622 = 0.1666666667e0_dp * t165 * t1392 * t4729
            t5625 = 0.1666666667e0_dp * t429 * t166 * t5339
            t5628 = 0.3600000000e2_dp * t2081 * t4760 * t425
            t5631 = 0.5555555556e-1_dp * t1246 * t4032 * t1483
            t5632 = t25 * t4853
            t5635 = 0.1500000000e1_dp * t405 * t156 * t5632
            t5638 = 0.2250000000e1_dp * t400 * t1333 * t4701
            t5641 = 0.1500000000e1_dp * t405 * t1333 * t4729
            t5643 = t400 * t3975 * t561
            t5647 = 0.5625000000e1_dp * t1302 * t1309 * t4701
            t5650 = 0.5625000000e1_dp * t1302 * t401 * t4801
            t5653 = 0.4500000000e1_dp * t400 * t1892 * t1447
            t5657 = 0.5000000000e0_dp * t405 * t68 * t5122 * t17
            t5660 = 0.1125000000e2_dp * t1302 * t401 * t4797
            t5663 = 0.1500000000e1_dp * t400 * t417 * t5371
            t5666 = 0.1125000000e2_dp * t1302 * t1892 * t1410
            t5667 = t5635 - t5638 + t5641 + 0.4500000000e1_dp * t5643 - t5647 - t5650  &
            + t5653 - t5657 - t5660 + t5663 - t5666
            t5670 = 0.2250000000e1_dp * t400 * t1309 * t4729
            t5671 = t172 * t4801
            t5673 = 0.2250000000e1_dp * t1912 * t5671
            t5676 = 0.2250000000e1_dp * t400 * t401 * t4853
            t5679 = 0.1875000000e1_dp * t1302 * t406 * t5002
            t5682 = 0.1500000000e1_dp * t400 * t1902 * t1483
            t5686 = 0.2250000000e1_dp * t400 * t19 * t5122 * t237
            t5689 = 0.1000000000e1_dp * t405 * t417 * t5430
            t5692 = 0.1968750000e2_dp * t2823 * t401 * t4790
            t5695 = 0.7500000000e0_dp * t400 * t406 * t5009
            t5696 = t172 * t4790
            t5698 = 0.5625000000e1_dp * t3939 * t5696
            t5699 = t172 * t4797
            t5701 = 0.4500000000e1_dp * t1912 * t5699
            t5705 = t5670 - t5673 + t5676 - t5679 + t5682 + t5686 - t5689 + t5692  &
            + t5695 + t5698 - t5701 - 0.1500000000e1_dp * t151 * t50 * t5613
            t5713 = t4739 * t160
            t5716 = 0.1620000000e3_dp * t4994 * t392 * t4769 * t1141 * t26 * t5713
            t5719 = 0.3333333334e0_dp * t165 * t1966 * t1447
            t5722 = 0.1111111111e0_dp * t165 * t689 * t1988
            t5724 = t429 * t166 * t5345
            t5728 = 0.2500000000e0_dp * t1981 * t1247 * t4797
            t5731 = 0.1666666667e0_dp * t445 * t446 * t5632
            t5734 = 0.1666666667e0_dp * t429 * t679 * t1553
            t5738 = 0.4444444445e0_dp * t170 * t114 * t4782 * t451
            t5741 = 0.1666666667e0_dp * t165 * t171 * t5143
            t5742 = -t5426 + t5429 - 0.2222222222e0_dp * t5432 + t5436 + 0.1666666667e0_dp  &
            * t165 * t166 * t5613 + t5619 + t5622 - t5625 - t5628 + t5631  &
            + (2 * t147 * t149 * (t5667 + t5705)) + t5716 + t5719 - t5722  &
            - 0.1666666667e0_dp * t5724 + t5728 + t5731 - t5734 - t5738 + t5741
            t5745 = t165 * t4498 * t561
            t5749 = 0.1666666667e0_dp * t429 * t437 * t5170
            t5752 = 0.8333333335e-1_dp * t429 * t462 * t4801
            t5756 = 0.1111111111e0_dp * t170 * t27 * t5276 * t17
            t5760 = 0.1250000000e0_dp * t2763 * t446 * t413 * t4701
            t5767 = 0.1666666667e0_dp * t165 * t437 * t5122
            t5770 = 0.1620000000e3_dp * t4995 * t4772 * t425
            t5773 = 0.1250000000e0_dp * t1368 * t437 * t5163
            t5776 = 0.2777777778e-1_dp * t1246 * t1247 * t5009
            t5780 = 0.1666666667e0_dp * t165 * t56 * t5276 * t237
            t5783 = 0.5555555557e-1_dp * t165 * t5013 * t433
            t5786 = 0.3333333334e0_dp * t165 * t689 * t668
            t5789 = 0.1250000000e0_dp * t1368 * t166 * t5330
            t5792 = 0.2777777778e-1_dp * t429 * t466 * t5002
            t5795 = 0.1111111111e0_dp * t165 * t1940 * t1483
            t5797 = t1373 * t1962 * t1829
            t5801 = 0.3600000000e2_dp * t1143 * t479 * t5002
            t5804 = 0.2400000000e2_dp * t473 * t1811 * t1483
            t5807 = 0.2160000000e3_dp * t2994 * t474 * t4790
            t5810 = 0.1080000000e3_dp * t1143 * t1801 * t1410
            t5813 = 0.1800000000e2_dp * t478 * t183 * t5632
            t5816 = 0.3600000000e2_dp * t473 * t1176 * t4701
            t5819 = 0.1800000000e2_dp * t478 * t1176 * t4729
            t5822 = 0.1800000000e2_dp * t473 * t474 * t4853
            t5825 = 0.1200000000e2_dp * t473 * t479 * t5009
            t5828 = 0.3600000000e2_dp * t473 * t1801 * t1447
            t5830 = 0.1080000000e3_dp * t4230 * t5696
            t5831 = t5801 - t5804 - t5807 + t5810 - t5813 + t5816 - t5819 - t5822  &
            - t5825 - t5828 - t5830
            t5833 = 0.7200000000e2_dp * t1821 * t5699
            t5835 = 0.3600000000e2_dp * t1821 * t5671
            t5838 = 0.1080000000e3_dp * t1143 * t474 * t4797
            t5841 = 0.5400000000e2_dp * t1143 * t474 * t4801
            t5845 = 0.1800000000e2_dp * t473 * t20 * t5122 * t237
            t5848 = 0.5400000000e2_dp * t1143 * t1150 * t4701
            t5850 = t473 * t4201 * t561
            t5854 = 0.1200000000e2_dp * t473 * t489 * t5371
            t5857 = 0.6000000000e1_dp * t478 * t489 * t5430
            t5861 = 0.5999999999e1_dp * t478 * t56 * t5122 * t17
            t5864 = 0.1800000000e2_dp * t473 * t1150 * t4729
            t5868 = t5833 + t5835 + t5838 + t5841 - t5845 + t5848 - 0.3600000000e2_dp  &
            * t5850 - t5854 + t5857 + t5861 - t5864 + 0.8999999998e1_dp * t178  &
            * t62 * t5613
            t5869 = t5831 + t5868
            t5875 = t1447 * t187 * t1722
            t5880 = 0.3703703705e-1_dp * t192 * t194 * t5226 * t526
            t5881 = t1141 * t4701
            t5885 = t719 * t1769
            t5888 = t60 * t1836 * t709
            t5891 = 0.5555555558e-1_dp * t4943 * t534
            t5898 = t705 * t1504 * t65
            t5902 = 0.1800000000e2_dp * t1718 * t178 * t710
            t5905 = 0.5555555558e-1_dp * t192 * t5270 * t267
            t5908 = t1836 * t59 * t177 * t710
            t5911 = t266 * t4729 * t65
            t5913 = 0.1000000001e1_dp * t723 * t5911
            t5918 = t140 * t1504 * t1722
            t5922 = 0.1111111112e0_dp * t4776 * t515 * t517
            t5925 = 0.8099999996e2_dp * t498 * t4770 * t4915
            t5926 = -0.5555555558e-1_dp * t192 * t68 * t5869 * t65 - 0.1000000000e1_dp  &
            * t5875 + t5880 + 0.3000000000e1_dp * t1748 * t5881 * t4940 - 0.1111111112e0_dp  &
            * t5885 - 0.1000000000e1_dp * t5888 - t5891 + t5869 * t65 +  &
            0.666666666e0_dp * t1519 * t722 * t4729 * t65 + (2 * t5898) + t5902  &
            - t5905 + 0.1800000000e2_dp * t5908 - t5913 - 0.1111111112e0_dp * t516  &
            * t4898 * t172 - 0.1000000000e1_dp * t5918 - t5922 + t5925
            t5929 = 0.2250000000e2_dp * t1740 * t471 * t4905
            t5932 = 0.1000000000e1_dp * t649 * t266 * t1722
            t5938 = t192 * t728 * t1505
            t5940 = t4918 * t67
            t5942 = 0.1111111112e0_dp * t5940 * t196
            t5944 = t1447 * t722 * t1722
            t5949 = 0.8999999998e1_dp * t1795 * t4909
            t5952 = 0.9000000007e1_dp * t1763 * t5881 * t4940
            t5954 = 0.8999999998e1_dp * t1741 * t4909
            t5956 = t471 * t4701 * t65
            t5958 = 0.1000000000e1_dp * t1766 * t5956
            t5960 = 0.5555555558e-1_dp * t5940 * t200
            t5963 = 0.1000000000e1_dp * t379 * t471 * t4904
            t5964 = t1715 * t1525
            t5966 = t1785 * t1790
            t5971 = 0.4500000000e1_dp * t379 * t1141 * t4701 * t4940
            t5973 = 0.5555555558e-1_dp * t4943 * t531
            t5975 = t5226 * t266 * t65
            t5977 = 0.1111111112e0_dp * t1755 * t729
            t5978 = -t5929 - t5932 + 0.3333333334e0_dp * t1519 * t187 * t4729 * t65  &
            - 0.1111111112e0_dp * t5938 - t5942 - 0.2000000002e1_dp * t5944 - 0.5555555558e-1_dp  &
            * t503 * t4955 + t5949 - t5952 + t5954 + t5958 - t5960  &
            + t5963 - 0.1111111112e0_dp * t5964 - 0.2222222224e0_dp * t5966 - t5971  &
            - t5973 + t5975 - t5977
            t5985 = 0.5000000001e0_dp * t1774 * t4908
            t5990 = t177 * t4729 * t65
            t5994 = 0.5000000001e0_dp * t1766 * t5990
            t5996 = 0.3703703705e-1_dp * t4943 * t527
            t6000 = 0.9000000007e1_dp * t5262 * t266 * t4701 * t4940
            t6002 = 0.1851851853e0_dp * t4943 * t512
            t6004 = t1709 * t178 * t710
            t6006 = t719 * t1838
            t6009 = 0.5000000001e0_dp * t1777 * t4908
            t6017 = 0.5000000001e0_dp * t732 * t5911
            t6019 = 0.7407407410e-1_dp * t719 * t1726
            t6021 = t649 * t1518 * t1722
            t6023 = t5251 * t67
            t6025 = 0.3703703705e-1_dp * t6023 * t271
            t6029 = 0.5999999993e1_dp * t4946 * t722 * t4701 * t4940 - t5985 - 0.5555555558e-1_dp  &
            * t192 * t199 * t4898 + 0.3333333334e0_dp * t1748 * t5990  &
            - t5994 + t5996 - t6000 + t6002 + 0.1800000000e2_dp * t6004 - 0.1111111112e0_dp  &
            * t6006 - t6009 + 0.3000000000e1_dp * t4946 * t187 * t4701  &
            * t4940 - 0.6666666668e0_dp * t1748 * t5956 - t6017 + t6019 + 0.6666666668e0_dp  &
            * t6021 + t6025 - 0.5555555558e-1_dp * t58 * t5613 * t71
            t6031 = 0.5000000001e0_dp * t1720 * t5990
            t6033 = 0.1000000000e1_dp * t5244 * t1712
            t6036 = 0.2250000000e2_dp * t498 * t472 * t4905
            t6040 = 0.5555555558e-1_dp * t6023 * t275
            t6042 = 0.1000000001e1_dp * t1763 * t5990
            t6044 = 0.2000000002e1_dp * t1763 * t5956
            t6047 = 0.4500000000e1_dp * t1720 * t5881 * t4940
            t6049 = t1688 * t177 * t709
            t6052 = t1504 * t561 * t65
            t6053 = t732 * t6052
            t6056 = 0.5555555558e-1_dp * t508 * t5271
            t6058 = 0.1111111112e0_dp * t4943 * t520
            t6060 = t1766 * t5881 * t4940
            t6063 = 0.1000000000e1_dp * t1720 * t5956
            t6066 = t1519 * t705 * t561 * t65
            t6074 = 0.8099999996e2_dp * t496 * t4769 * t1141 * t4915
            t6075 = t723 * t6052
            t6078 = 0.1111111112e0_dp * t719 * t1782
            t6079 = -t6031 - t6033 - t6036 + t187 * t4897 * t65 - t6040 - t6042  &
            + t6044 - t6047 - 0.1000000000e1_dp * t6049 - 0.1000000000e1_dp * t6053  &
            - t6056 - t6058 - 0.4500000000e1_dp * t6060 + t6063 + 0.6666666668e0_dp  &
            * t6066 - 0.1000000000e1_dp * t237 * t705 * t1722 + t6074 - 0.2000000002e1_dp  &
            * t6075 - t6078
            t6087 = 0.1666666667e0_dp * t429 * t1966 * t1410
            t6090 = 0.3125000000e0_dp * t4025 * t1247 * t4790
            t6091 = 0.3333333334e0_dp * t5745 - t5749 - t5752 - t5756 + t5760 + 0.1666666667e0_dp  &
            * t165 * t56 * t4958 * t140 + t5767 + t5770 + t5773 +  &
            t5776 + t5780 - t5783 + t5786 + t5789 + t5792 - t5795 - 0.4444444444e0_dp  &
            * t5797 + 0.3333333334e0_dp * t55 * t57 * (t5926 + t5978 + t6029  &
            + t6079) - t6087 - t6090
            t6094 = 0.1250000000e0_dp * t1981 * t1247 * t4801
            t6099 = 0.1800000000e2_dp * t1971 * t1972 * t4729 * t149 * t160
            t6102 = 0.2500000000e0_dp * t1981 * t4553 * t1410
            t6106 = 0.3600000000e2_dp * t1971 * t1972 * t1973 * t672
            t6109 = 0.1666666667e0_dp * t1246 * t1962 * t1865
            t6112 = 0.5555555556e-1_dp * t1246 * t4553 * t1483
            t6116 = 0.3600000000e2_dp * t1971 * t471 * t20 * t5713
            t6119 = 0.1250000000e0_dp * t1368 * t462 * t4790
            t6123 = 0.1666666667e0_dp * t1246 * t446 * t1825 * t561
            t6126 = 0.1800000000e2_dp * t658 * t4764 * t425
            t6129 = 0.8333333335e-1_dp * t429 * t437 * t5174
            t6132 = 0.1666666667e0_dp * t165 * t462 * t4853
            t6136 = 0.5555555557e-1_dp * t165 * t171 * t5122 * t17
            t6138 = t429 * t1928 * t593
            t6141 = t658 * t575 * t1925
            t6144 = t445 * t4512 * t668
            t6155 = 0.2500000000e0_dp * t1981 * t4032 * t1410
            t6158 = 0.1666666667e0_dp * t429 * t679 * t1571
            t6160 = t147 * t1885 * t672
            t6162 = t6094 + t6099 + t6102 + t6106 - t6109 + t6112 - t6116 + t6119  &
            - t6123 + t6126 - t6129 + t6132 - t6136 - 0.1666666667e0_dp * t6138  &
            + 0.3600000000e2_dp * t6141 + 0.3333333334e0_dp * t6144 + 0.2e1_dp * t147  &
            * (t4881 - t4884 - t4887 + 0.9000000000e1_dp * t4889 - t4893 + 0.6000000000e1_dp  &
            * t4895) * t149 * t160 + t6155 - t6158 + (4 * t6160)
            t6170 = 0.3333333336e0_dp * t289 * t597 * t1571
            t6173 = 0.3333333336e0_dp * t55 * t5295 * t140
            t6175 = t83 * t597 * t1688
            t6179 = 0.3333333336e0_dp * t165 * t3865 * t561
            t6181 = t289 * t1564 * t593
            t6185 = 0.3333333336e0_dp * t289 * t597 * t1553
            t6186 = -0.3333333336e0_dp * t5347 - t5351 - t5354 - t5357 - t5360 - 0.2222222224e0_dp  &
            * t144 * t46 * (t5423 + t5742 + t6091 + t6162) - t6170  &
            - t6173 - 0.6666666672e0_dp * t6175 - t6179 - 0.3333333336e0_dp * t6181  &
            - t6185
            t6189 = 0.8333333340e-1_dp * t906 * t294 * t5163
            t6192 = 0.3333333336e0_dp * t83 * t294 * t5122
            t6194 = t83 * t3875 * t561
            t6198 = 0.3333333336e0_dp * t289 * t294 * t5170
            t6201 = 0.3333333336e0_dp * t289 * t597 * t1557
            t6207 = 0.1666666668e0_dp * t1546 * t897 * t4797
            t6210 = 0.1666666668e0_dp * t289 * t294 * t5174
            t6213 = 0.3333333336e0_dp * t55 * t84 * t5122
            t6216 = 0.8333333340e-1_dp * t429 * t897 * t4701
            t6219 = 0.8333333340e-1_dp * t906 * t384 * t4790
            t6222 = 0.8333333340e-1_dp * t1546 * t897 * t4801
            t6223 = t6189 - t6192 - 0.6666666672e0_dp * t6194 - t6198 - t6201 - 0.3333333336e0_dp  &
            * t83 * t84 * t5613 + t6207 - t6210 - t6213 + t6216 +  &
            t6219 + t6222
            t6227 = -t4704 - t4707 - t4732 - t4789 + t4793 - t4796 - t4800 - t4804  &
            - t4808 - t4812 - t4857 - t4861 - 0.1333333334e1_dp * t4864 - t4869  &
            - 0.4444444448e0_dp * t43 * t46 * t5036 + t8 * (t5305 + t5343 + t6186  &
            + t6223)
            e_ndrho_ndrho_rho_spin(ii) = e_ndrho_ndrho_rho_spin(ii) +  t6227 * sx
            t6247 = t4769 * t2992
            t6249 = t4772 * t2027
            t6252 = t56 * t4701
            t6253 = t6252 * t526
            t6256 = t56 * t4729
            t6257 = t6256 * t526
            t6260 = t4743 * t526
            t6263 = t20 * t4853
            t6267 = t26 * t1447
            t6268 = t6267 * t709
            t6271 = t4764 * t2027
            t6281 = -0.1620000000e3_dp * t6247 * t4733 * t6249 - 0.1200000000e2_dp *  &
            t473 * t6253 + 0.5999999999e1_dp * t478 * t6257 + 0.5399999998e2_dp * t4771  &
            * t6260 + 0.8999999998e1_dp * t478 * t6263 * t65 + 0.8099999996e2_dp  &
            * t4771 * t6268 - 0.1800000000e2_dp * t473 * t6271 - 0.4050000000e2_dp  &
            * t473 * t575 * t2095 - 0.5000000001e0_dp * t60 * t4853 * t65 + 0.8999999998e1_dp  &
            * t4901 - t4907 + t4911
            t6283 = t4701 * t1518
            t6290 = t4917 - t4920 + 0.5000000001e0_dp * t4922 + t4925 + 0.3333333334e0_dp  &
            * t6283 * t177 * t65 - 0.5000000001e0_dp * t4929 - t4932 - t4935  &
            + 0.5850000000e2_dp * t4941 + t4945 - t4951 - 0.5555555558e-1_dp * t4952
            t6291 = t6281 + t6290
            t6308 = 0.4050000000e2_dp * t4734 * t4736 * t2821 * t18 * t4739 * t237  &
            - 0.1125000000e2_dp * t2013 * t1300 * t19 * t4790 - t5028 - t5022 +  &
            0.1500000000e1_dp * t2013 * t398 * t68 * t5002 + 0.3333333334e0_dp * t55  &
            * t57 * t6291 + t5019 - 0.3000000000e1_dp * t566 * t1885 * t569 - t5016  &
            - 0.9000000000e1_dp * t4734 * t4736 * t1300 * t4701 * t149 * t17  &
            + t5032 + 0.1666666667e0_dp * t165 * t2088 * t561 - 0.5555555557e-1_dp *  &
            t5011 + t5001
            t6338 = 0.1666666667e0_dp * t4874 - t5025 - 0.2700000000e2_dp * t4734 * t392  &
            * t4736 * t1300 * rho * t5002 * t149 - 0.3000000000e1_dp * t2009  &
            * t4748 - 0.2700000000e2_dp * t4738 * t18 * t1973 * t1447 - t5008 -  &
            0.3000000000e1_dp * t567 * t151 * t19 * t4853 + t5035 - 0.1000000000e1_dp  &
            * t2013 * t2018 * t5009 + 0.9000000000e1_dp * t2013 * t2014 * t4797  &
            - t4872 + 0.4500000000e1_dp * t2013 * t2014 * t4801 + t4878 + t5005  &
            + 0.4500000000e1_dp * t2009 * t4744
            t6339 = t6308 + t6338
            t6344 = -t5043 - t5047 - t5050 - t5126 - t5284 - t5291 + t5294 - t5298  &
            - t5301 - 0.3333333336e0_dp * t5303 - t5308 - t5311 - t5314
            t6352 = t5440 - t5443 - t5446 - t5449 - t5451 - t5454 + 0.40e1_dp * t102  &
            * t327 * t2057 * t557 - t5456 + t5459 - t5462 + t5466 - t5468
            t6363 = t5480 - t5482 - t5484 + t5487 - t5489 + t5492 + t5494 - 0.9600000000e2_dp  &
            * t5496 + t5499 - t5503 - t5507 + t5509 - 0.2400000000e2_dp  &
            * t5511 - 0.960e2_dp * t5115 * t1060 - 0.240e2_dp * t2054 * t3674 + 0.1200e3_dp  &
            * t1054 * t640 * t3679
            t6367 = t5472 - t5474 - t5478 - t5523 + t5525 + 0.40e1_dp * t5528 + t5531  &
            + t5533 + t5535 + t5537 - 0.20e1_dp * t102 * t105 * t6363 + t5540
            t6370 = t5545 - t5548 + t5551 + t5553 - t5558 + t5560 - t5562 + t5564  &
            - 0.40e1_dp * t5565 + t5568 + t5572 + t5574
            t6373 = t5579 - t5581 - 0.40e1_dp * t112 * t2058 - t5585 + t5590 - t5593  &
            - t5595 + t5597 - t5600 - t5602 + t5604 + t5607 + t5610
            t6375 = t6352 + t6367 + t6370 + t6373
            t6380 = -t5317 - t5320 + t5323 - t5326 - t5329 - 0.3333333336e0_dp * t83  &
            * t597 * t2061 - 0.3333333336e0_dp * t83 * t84 * t6375 + t5333 - t5335  &
            - t5338 - t5342 - 0.1666666668e0_dp * t5347 - t5351
            t6389 = -t5354 - t5357 - 0.3333333336e0_dp * t83 * t45 * t2139 * t561  &
            - t5360 - t6170 - t6173 - 0.3333333336e0_dp * t6175 - t6179 - 0.1666666668e0_dp  &
            * t6181 - t6185 + t6189 - t6192 - 0.3333333336e0_dp * t6194
            t6400 = t2061 * t561
            t6415 = -t5363 + t5367 - t5370 + 0.1111111112e0_dp * t5373 - 0.5555555557e-1_dp  &
            * t5376 + 0.1666666667e0_dp * t5379 - t5383 - t5386 - t5389 +  &
            0.1200000000e2_dp * t658 * t2031 * t673 * t17 + 0.3333333334e0_dp * t170  &
            * t27 * t6291 * t172 - t5397
            t6428 = t5635 - t5638 + t5641 + 0.2250000000e1_dp * t5643 - t5647 - 0.1500000000e1_dp  &
            * t151 * t50 * t6375 - t5650 + t5653 - t5657 - t5660  &
            + t5663 - t5666
            t6433 = t5670 - t5673 + t5676 - t5679 + t5682 + t5686 - t5689 + t5692  &
            + t5695 + 0.2250000000e1_dp * t400 * t19 * t2061 * t561 + t5698 -  &
            t5701
            t6442 = t5400 + 0.1666666667e0_dp * t5402 - t5406 + 0.1666666667e0_dp * t165  &
            * t56 * t2131 * t561 + 0.1666666667e0_dp * t165 * t679 * t2061 +  &
            (2 * t147 * t149 * (t6428 + t6433)) - t5409 - t5412 - t5415 +  &
            t5418 + 0.1800000000e2_dp * t658 * t575 * t2074 + 0.1800000000e2_dp * t5421
            t6451 = -t5426 + t5429 - 0.2222222223e0_dp * t5432 + t5436 + 0.1620000000e3_dp  &
            * t4995 * t6267 * t1974 + 0.1666666667e0_dp * t165 * t2088 * t649  &
            + t5619 + t5622 - t5625 - t5628 + t5631 + t5716
            t6461 = t5719 - t5722 - 0.8333333335e-1_dp * t5724 + 0.1200000000e2_dp *  &
            t658 * t6256 * t2077 - 0.3240000000e3_dp * t4994 * t6247 * t4772 * t2082  &
            + t5728 + t5731 - t5734 - t5738 + t5741 + 0.1666666667e0_dp * t5745  &
            - t5749
            t6479 = -t5752 + 0.1666666667e0_dp * t445 * t2042 * rho * t668 - t5756  &
            - 0.3600000000e2_dp * t2081 * t575 * t673 * t237 + 0.1800000000e2_dp *  &
            t658 * t6263 * t161 + t5760 + t5767 + t5770 - 0.7200000000e2_dp * t2081  &
            * t575 * t161 * t1447 + t5773 + t5776 + t5780
            t6490 = -t5783 - 0.3600000000e2_dp * t2081 * t4764 * t2082 - 0.2400000000e2_dp  &
            * t2081 * t6252 * t2077 + t5786 + t5789 + t5792 - t5795 - 0.4444444445e0_dp  &
            * t5797 - t6087 - 0.8333333335e-1_dp * t429 * t166 * t6400  &
            - t6090 + t6094
            t6499 = t6099 + 0.1666666667e0_dp * t165 * t166 * t6375 + t6102 + 0.1666666667e0_dp  &
            * t165 * t56 * t6291 * t140 + t6106 - t6109 + t6112 - t6116  &
            + t6119 - t6123 + t6126 - t6129
            t6509 = t4853 * t65
            t6529 = 0.4999999999e0_dp * t5875 - 0.5000000001e0_dp * t5244 * t2095 - 0.5000000001e0_dp  &
            * t732 * t6509 + t5880 + 0.1000000000e1_dp * t649 * t471  &
            * t2099 - 0.1620000000e3_dp * t5247 * t2992 * t4733 * t6249 - 0.5555555558e-1_dp  &
            * t5885 + 0.8999999998e1_dp * t140 * t2992 * t4701 * t59 *  &
            t62 * t2027 - 0.5000000001e0_dp * t5888 - t5891 - 0.1800000000e2_dp * t2118  &
            * t6271
            t6539 = t5898 + t5902 - t5905 + 0.8999999998e1_dp * t5908 - t5913 + 0.2000000000e1_dp  &
            * t2098 * t2095 * t561 - 0.5000000001e0_dp * t5918 - t5922  &
            + 0.3333333334e0_dp * t6283 * t187 * t177 * t65 + t5925 - t5929
            t6554 = -t5932 - 0.5555555558e-1_dp * t5938 + 0.666666666e0_dp * t6283 *  &
            t722 * t177 * t65 - t5942 + 0.8999999998e1_dp * t2124 * t59 * t177 *  &
            t710 + 0.1000000001e1_dp * t5944 - 0.5000000001e0_dp * t2061 * t177 * t709  &
            + t5949 + t5952 + t5954 - t5958
            t6557 = t1447 * t59
            t6569 = -t5960 - 0.5000000001e0_dp * t5237 * t2095 + t5963 - 0.4500000000e1_dp  &
            * t5256 * t6557 * t710 - 0.5555555558e-1_dp * t5964 - 0.1111111112e0_dp  &
            * t5966 - t5971 - 0.1200000000e2_dp * t2118 * t6253 - 0.5555555558e-1_dp  &
            * t719 * t2126 - t5973 + t5975 + 0.5999999999e1_dp * t2114 * t6257
            t6583 = -t5977 - t5985 + t5994 - 0.5000000001e0_dp * t716 * t6509 + t5996  &
            - t6000 + t6002 + 0.8999999998e1_dp * t6004 + 0.1000000000e1_dp * t2098  &
            * t4908 * t237 - 0.1000000001e1_dp * t723 * t6509 - 0.1800000000e2_dp  &
            * t5228 * t2117 * t2028
            t6596 = t5247 * t1141 * t4733
            t6603 = -0.5555555558e-1_dp * t6006 - t6009 - t6017 - 0.5555555558e-1_dp  &
            * t58 * t6375 * t71 + t6019 - 0.9000000007e1_dp * t5262 * t6557 * t710  &
            + 0.8999999998e1_dp * t708 * t62 * t6509 + 0.3333333334e0_dp * t6021 +  &
            t6025 + 0.8099999996e2_dp * t6596 * t6268 + 0.5399999998e2_dp * t6596 *  &
            t6260 + 0.8999999998e1_dp * t5229 * t2104
            t6611 = t5801 - t5804 - t5807 + t5810 - t5813 + t5816 - t5819 - t5822  &
            - t5825 - t5828 - t5830 + t5833
            t6620 = t5835 + t5838 + t5841 - t5845 + t5848 - 0.1800000000e2_dp * t5850  &
            - t5854 + t5857 + t5861 - t5864 - 0.1800000000e2_dp * t473 * t20  &
            * t2061 * t561 + 0.8999999998e1_dp * t178 * t62 * t6375
            t6621 = t6611 + t6620
            t6630 = -0.3000000000e1_dp * t5256 * t4904 * t237 + 0.5999999999e1_dp * t5228  &
            * t2113 * t2032 + t6621 * t65 - t6031 - 0.4050000000e2_dp * t471  &
            * t187 * t6557 * t710 - t6033 - t6036 - t6040 + t6042 - 0.5000000001e0_dp  &
            * t60 * t2124 * t709 - t6044
            t6640 = -t6047 - 0.5000000001e0_dp * t6049 - 0.5555555558e-1_dp * t192 *  &
            t68 * t6621 * t65 - 0.5000000001e0_dp * t6053 - t6056 - t6058 + 0.5850000000e2_dp  &
            * t6060 + t6063 + 0.3333333334e0_dp * t6066 + t6074 - 0.1000000001e1_dp  &
            * t6075 - t6078
            t6658 = 0.1800000000e2_dp * t658 * t2035 * t673 - 0.8333333335e-1_dp * t429  &
            * t2088 * t593 + t6132 + 0.3333333334e0_dp * t55 * t57 * (t6529 +  &
            t6539 + t6554 + t6569 + t6583 + t6603 + t6630 + t6640) - t6136 - 0.8333333335e-1_dp  &
            * t6138 + 0.1800000000e2_dp * t6141 + 0.1080000000e3_dp *  &
            t4995 * t4743 * t2077 + 0.1080000000e3_dp * t657 * t1142 * t4760 * t2082  &
            + 0.1666666667e0_dp * t6144 + t6155 - t6158 + (2 * t6160)
            t6665 = -0.3333333336e0_dp * t83 * t45 * t6339 * t140 - t6198 - 0.3333333336e0_dp  &
            * t83 * t2050 * t649 - t6201 + t6207 - t6210 - t6213 + t6216  &
            + t6219 - 0.1666666668e0_dp * t289 * t2050 * t593 - 0.1666666668e0_dp  &
            * t289 * t84 * t6400 + t6222 - 0.2222222224e0_dp * t144 * t46 * (t6415  &
            + t6442 + t6451 + t6461 + t6479 + t6490 + t6499 + t6658)
            t6669 = -t4704 - t4707 - t4732 + t4793 - t4796 - t4800 - t4804 - t4861  &
            - 0.6666666672e0_dp * t4864 - t4869 - t4812 - t4857 - t4789 - t4808  &
            - 0.6666666672e0_dp * t217 * t46 * t2046 * t561 - 0.4444444448e0_dp *  &
            t43 * t46 * t6339 + t8 * (t6344 + t6380 + t6389 + t6665)
            e_ndrho_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_ndrho_spin(ii) +  t6669 * sx
          END IF
        END IF
      END DO
      !$omp end do
    END IF

  END SUBROUTINE xb88_lr_adiabatic_lsd_calc

END MODULE xc_xbecke88_lr_adiabatic
