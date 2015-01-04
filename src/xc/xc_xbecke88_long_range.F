!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the longrange part of Becke 88 exchange functional
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE xc_xbecke88_long_range
  USE bibliography,                    ONLY: Becke1988,&
                                             cite_reference
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE erf_fn,                          ONLY: erf
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
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xbecke88_long_range'
  REAL(kind=dp), PARAMETER :: beta=0.0042_dp

  PUBLIC :: xb88_lr_lda_info, xb88_lr_lsd_info, xb88_lr_lda_eval, xb88_lr_lsd_eval
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_lda_info(reference,shortform, needs, max_deriv,error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_lda_info', &
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
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE xb88_lr_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_lsd_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_lsd_info', &
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

  END SUBROUTINE xb88_lr_lsd_info

! *****************************************************************************
!> \brief evaluates the becke 88 longrange exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param xb88_lr_params input parameters (scaling)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_lda_eval(rho_set,deriv_set,grad_deriv,xb88_lr_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xb88_lr_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, omega, sx
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, &
      e_ndrho_rho_rho, e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(xb88_lr_params,"scale_x",r_val=sx,error=error)
    CALL section_vals_val_get(xb88_lr_params, "omega",r_val=omega, error=error)

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
       !$omp          shared(epsilon_rho, sx, omega, error)

       CALL xb88_lr_lda_calc(rho=rho, norm_drho=norm_drho,&
            e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
            e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
            e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
            e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,&
            e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,grad_deriv=grad_deriv,&
            npoints=npoints,epsilon_rho=epsilon_rho,&
            sx=sx,omega=omega,error=error)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE xb88_lr_lda_eval

! *****************************************************************************
!> \brief low level calculation of the becke 88 exchange functional for lda
!> \param rho alpha or beta spin density
!> \param norm_drho || grad rho ||
!> \param e_0 adds to it the local value of the functional
!> \param e_rho derivative of the functional wrt. to the variables
!>        named where the * is.
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param e_rho_rho_rho ...
!> \param e_ndrho_rho_rho ...
!> \param e_ndrho_ndrho_rho ...
!> \param e_ndrho_ndrho_ndrho ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sx scaling-parameter for exchange
!> \param omega parameter for erfc
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
!> \note
!>      - Just took the lsd code and scaled rho and ndrho by 1/2 (e_0 with 2.0)
!>      - Derivatives higher than 1 not tested
! *****************************************************************************
  SUBROUTINE xb88_lr_lda_calc(rho, norm_drho,&
       e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
       e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
       e_ndrho_ndrho_ndrho,grad_deriv,npoints,epsilon_rho, sx,&
       omega, error)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: &
      e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho_rho, e_rho_rho_rho, &
      e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, omega
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: Cx, epsilon_rho43, my_epsilon_rho, my_ndrho, my_rho, t1, &
      t10, t100, t1002, t1009, t101, t1013, t103, t104, t1044, t105, t1069, &
      t109, t1091, t11, t1102, t112, t113, t1136, t1141, t1146, t1153, t1156, &
      t116, t1163, t1167, t117, t1177, t118, t12, t122, t1220, t126, t127, &
      t128, t1284, t130, t132, t133, t1334, t1341, t1345, t135, t136, t137, &
      t1370, t1396, t140, t1400, t1404, t1405, t141, t1412, t143, t1449, &
      t1456, t146, t1467, t147, t1472, t148, t151, t155, t1553, t16, t168, &
      t169, t17, t172, t173, t176, t177, t18, t185, t186, t190, t196, t2, &
      t200, t207, t21, t211, t212, t213, t216
    REAL(kind=dp) :: t219, t22, t221, t222, t225, t226, t23, t230, t231, &
      t232, t233, t237, t24, t241, t244, t245, t250, t251, t254, t255, t258, &
      t259, t26, t264, t265, t27, t270, t271, t28, t281, t284, t285, t289, &
      t29, t293, t297, t3, t30, t304, t31, t311, t312, t313, t316, t319, &
      t321, t323, t325, t326, t328, t330, t335, t339, t34, t343, t346, t347, &
      t351, t354, t358, t36, t365, t368, t37, t372, t377, t38, t382, t383, &
      t384, t39, t393, t397, t40, t400, t401, t404, t405, t408, t41, t412, &
      t413, t417, t418, t42, t428, t429, t43, t435, t44, t451, t452, t455, &
      t457, t46, t460, t462, t463, t464, t465
    REAL(kind=dp) :: t466, t467, t47, t470, t473, t478, t479, t48, t482, &
      t486, t489, t49, t492, t496, t5, t501, t502, t505, t508, t51, t513, &
      t514, t519, t52, t521, t525, t526, t529, t530, t533, t534, t536, t537, &
      t539, t55, t562, t566, t570, t571, t572, t573, t574, t575, t576, t580, &
      t586, t59, t590, t6, t60, t601, t605, t606, t613, t624, t628, t632, &
      t639, t64, t640, t641, t657, t667, t669, t677, t68, t687, t69, t7, t71, &
      t716, t72, t722, t735, t739, t746, t75, t76, t769, t77, t79, t791, &
      t793, t8, t838, t84, t842, t846, t85, t857, t86, t860, t867, t87, t880, &
      t90, t91, t933, t94, t95, t961, t98, t99
    REAL(kind=dp) :: xx

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
            t2 = t1 * my_rho
            t3 = 0.1e1_dp / t2
            xx = my_ndrho * MAX(t3,epsilon_rho43)
            t5 = my_ndrho ** 2
            t6 = beta * t5
            t7 = my_rho ** 2
            t8 = t1 ** 2
            t10 = 0.1e1_dp / t8 / t7
            t11 = beta * my_ndrho
            t12 = LOG(xx + SQRT(xx ** 0.2e1_dp + 0.1e1_dp))
            t16 = 0.10e1_dp + 0.60e1_dp * t11 * t3 * t12
            t17 = 0.1e1_dp / t16
            t18 = t10 * t17
            t21 = 0.20e1_dp * Cx + 0.20e1_dp * t6 * t18
            t22 = SQRT(t21)
            t23 = t22 * t21
            t24 = my_rho * t23
            t26 = rootpi
            t27 = 0.1e1_dp / t26
            t28 = 0.1e1_dp / omega
            t29 = 0.1e1_dp / t22
            t30 = t28 * t29
            t31 = t26 * t1
            t34 = erf(0.3000000000e1_dp * t30 * t31)
            t36 = omega * t22
            t37 = 0.1e1_dp / t1
            t38 = t27 * t37
            t39 = omega ** 2
            t40 = 0.1e1_dp / t39
            t41 = 0.1e1_dp / t21
            t42 = t40 * t41
            t43 = pi * t8
            t44 = t42 * t43
            t46 = EXP(-0.8999999998e1_dp * t44)
            t47 = t39 * t21
            t48 = 0.1e1_dp / pi
            t49 = 0.1e1_dp / t8
            t51 = t46 - 0.10e1_dp
            t52 = t48 * t49 * t51
            t55 = t46 - 0.15e1_dp - 0.5555555558e-1_dp * t47 * t52
            t59 = t26 * t34 + 0.3333333334e0_dp * t36 * t38 * t55
            t60 = t27 * t59

            !! Multiply with 2.0 because it code comes from LSD
            e_0(ii) = e_0(ii) -0.2222222224e0_dp * t24 * omega * t60 * sx * 2.0_dp

          END IF
          IF( grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
            t64 = t23 * omega
            t68 = my_rho * t22 * omega
            t69 = t7 * my_rho
            t71 = 0.1e1_dp / t8 / t69
            t72 = t71 * t17
            t75 = t16 ** 2
            t76 = 0.1e1_dp / t75
            t77 = t10 * t76
            t79 = 0.1e1_dp / t1 / t7
            t84 = 1 + t5 * t10
            t85 = SQRT(t84)
            t86 = 0.1e1_dp / t85
            t87 = t71 * t86
            t90 = -0.8000000000e1_dp * t11 * t79 * t12 - 0.8000000000e1_dp * t6 * t87
            t91 = t77 * t90
            t94 = -0.5333333333e1_dp * t6 * t72 - 0.20e1_dp * t6 * t91
            t95 = t60 * t94
            t98 = omega * t27
            t99 = SQRT(0.3141592654e1_dp)
            t100 = 0.1e1_dp / t99
            t101 = t26 * t100
            t103 = EXP(-0.9000000000e1_dp * t44)
            t104 = 0.1e1_dp / t23
            t105 = t28 * t104
            t109 = t26 * t49
            t112 = -0.1500000000e1_dp * t105 * t31 * t94 + 0.1000000000e1_dp * t30 *  &
            t109
            t113 = t103 * t112
            t116 = omega * t29
            t117 = t116 * t27
            t118 = t37 * t55
            t122 = t27 * t3
            t126 = t21 ** 2
            t127 = 0.1e1_dp / t126
            t128 = t40 * t127
            t130 = t128 * t43 * t94
            t132 = pi * t37
            t133 = t42 * t132
            t135 = 0.8999999998e1_dp * t130 - 0.5999999999e1_dp * t133
            t136 = t135 * t46
            t137 = t39 * t94
            t140 = t8 * my_rho
            t141 = 0.1e1_dp / t140
            t143 = t48 * t141 * t51
            t146 = t47 * t48
            t147 = t49 * t135
            t148 = t147 * t46
            t151 = t136 - 0.5555555558e-1_dp * t137 * t52 + 0.3703703705e-1_dp * t47  &
            * t143 - 0.5555555558e-1_dp * t146 * t148
            t155 = REAL(2 * t101 * t113,KIND=dp) + 0.1666666667e0_dp * t117 * t118 * t94 -  &
            0.1111111111e0_dp * t36 * t122 * t55 + 0.3333333334e0_dp * t36 * t38 *  &
            t151

            e_rho(ii) = e_rho(ii) +  (-0.2222222224e0_dp * t64 * t60 - 0.3333333336e0_dp * t68 * t95 -  &
            0.2222222224e0_dp * t24 * t98 * t155) * sx

            t168 = 0.60e1_dp * beta * t3 * t12 + 0.60e1_dp * t11 * t10 * t86
            t169 = t77 * t168
            t172 = 0.40e1_dp * t11 * t18 - 0.20e1_dp * t6 * t169
            t173 = t60 * t172
            t176 = pi * t100
            t177 = t176 * t103
            t185 = t128 * pi
            t186 = t8 * t172
            t190 = t39 * t172
            t196 = 0.8999999998e1_dp * t185 * t186 * t46 - 0.5555555558e-1_dp * t190  &
            * t52 - 0.5000000001e0_dp * t41 * t172 * t46
            t200 = -0.3000000000e1_dp * t177 * t105 * t1 * t172 + 0.1666666667e0_dp *  &
            t117 * t118 * t172 + 0.3333333334e0_dp * t36 * t38 * t196

            e_ndrho(ii) = e_ndrho(ii) +  (-0.3333333336e0_dp * t68 * t173 - 0.2222222224e0_dp * t24 * t98  &
            * t200) * sx

          END IF
          IF( grad_deriv >=2 .OR. grad_deriv == -2 ) THEN
            t207 = t27 * t155
            t211 = my_rho * t29 * omega
            t212 = t94 ** 2
            t213 = t60 * t212
            t216 = t207 * t94
            t219 = t7 ** 2
            t221 = 0.1e1_dp / t8 / t219
            t222 = t221 * t17
            t225 = t71 * t76
            t226 = t225 * t90
            t230 = 0.1e1_dp / t75 / t16
            t231 = t10 * t230
            t232 = t90 ** 2
            t233 = t231 * t232
            t237 = 0.1e1_dp / t1 / t69
            t241 = t221 * t86
            t244 = t5 ** 2
            t245 = beta * t244
            t250 = 0.1e1_dp / t85 / t84
            t251 = 0.1e1_dp / t1 / t219 / t69 * t250
            t254 = 0.1866666667e2_dp * t11 * t237 * t12 + 0.4000000000e2_dp * t6 * t241  &
            - 0.1066666667e2_dp * t245 * t251
            t255 = t77 * t254
            t258 = 0.1955555555e2_dp * t6 * t222 + 0.1066666667e2_dp * t6 * t226 + 0.40e1_dp  &
            * t6 * t233 - 0.20e1_dp * t6 * t255
            t259 = t60 * t258
            t264 = 0.9000000000e1_dp * t130 - 0.6000000000e1_dp * t133
            t265 = t264 * t103
            t270 = 0.1e1_dp / t22 / t126
            t271 = t28 * t270
            t281 = t26 * t141
            t284 = 0.2250000000e1_dp * t271 * t31 * t212 - 0.1000000000e1_dp * t105 *  &
            t109 * t94 - 0.1500000000e1_dp * t105 * t31 * t258 - 0.6666666667e0_dp  &
            * t30 * t281
            t285 = t103 * t284
            t289 = omega * t104 * t27
            t293 = t3 * t55
            t297 = t37 * t151
            t304 = t27 * t79
            t311 = t126 * t21
            t312 = 0.1e1_dp / t311
            t313 = t40 * t312
            t316 = 0.1800000000e2_dp * t313 * t43 * t212
            t319 = 0.1200000000e2_dp * t128 * t132 * t94
            t321 = t128 * t43 * t258
            t323 = pi * t3
            t325 = 0.2000000000e1_dp * t42 * t323
            t326 = -t316 + t319 + 0.8999999998e1_dp * t321 + t325
            t328 = t135 ** 2
            t330 = t39 * t258
            t335 = t137 * t48
            t339 = t48 * t10 * t51
            t343 = t141 * t135 * t46
            t346 = t49 * t326
            t347 = t346 * t46
            t351 = t49 * t328 * t46
            t354 = t326 * t46 + t328 * t46 - 0.5555555558e-1_dp * t330 * t52 + 0.7407407410e-1_dp  &
            * t137 * t143 - 0.1111111112e0_dp * t335 * t148 - 0.6172839508e-1_dp  &
            * t47 * t339 + 0.7407407410e-1_dp * t146 * t343 - 0.5555555558e-1_dp  &
            * t146 * t347 - 0.5555555558e-1_dp * t146 * t351
            t358 = REAL(2 * t101 * t265 * t112,KIND=dp) + REAL(2 * t101 * t285,KIND=dp) - 0.8333333335e-1_dp  &
            * t289 * t118 * t212 - 0.1111111111e0_dp * t117 * t293 * t94  &
            + 0.3333333334e0_dp * t117 * t297 * t94 + 0.1666666667e0_dp * t117 *  &
            t118 * t258 + 0.1481481481e0_dp * t36 * t304 * t55 - 0.2222222222e0_dp *  &
            t36 * t122 * t151 + 0.3333333334e0_dp * t36 * t38 * t354

            e_rho_rho(ii) = e_rho_rho(ii) +  (-0.6666666672e0_dp * t36 * t95 - 0.4444444448e0_dp * t64 * t207 -  &
            0.1666666668e0_dp * t211 * t213 - 0.6666666672e0_dp * t68 * t216 - 0.3333333336e0_dp  &
            * t68 * t259 - 0.2222222224e0_dp * t24 * t98 * t358) * sx

            t365 = t27 * t200
            t368 = t94 * t172
            t372 = t365 * t94
            t377 = t225 * t168
            t382 = t6 * t10
            t383 = t230 * t90
            t384 = t383 * t168
            t393 = beta * t5 * my_ndrho
            t397 = 0.1e1_dp / t1 / t219 / t7 * t250
            t400 = -0.8000000000e1_dp * beta * t79 * t12 - 0.2400000000e2_dp * t11 *  &
            t87 + 0.8000000000e1_dp * t393 * t397
            t401 = t77 * t400
            t404 = -0.1066666667e2_dp * t11 * t72 + 0.5333333333e1_dp * t6 * t377 - 0.40e1_dp  &
            * t11 * t91 + 0.40e1_dp * t382 * t384 - 0.20e1_dp * t6 * t401
            t405 = t60 * t404
            t408 = t207 * t172
            t412 = t26 * pi * t100
            t413 = t412 * t128
            t417 = t271 * t26
            t418 = t1 * t94
            t428 = 0.2250000000e1_dp * t417 * t418 * t172 - 0.1500000000e1_dp * t105  &
            * t31 * t404 - 0.5000000000e0_dp * t105 * t109 * t172
            t429 = t103 * t428
            t435 = t37 * t196
            t451 = t313 * pi
            t452 = t8 * t94
            t455 = 0.1800000000e2_dp * t451 * t452 * t172
            t457 = t128 * t43 * t404
            t460 = t128 * t132 * t172
            t462 = -t455 + 0.8999999998e1_dp * t457 + 0.5999999999e1_dp * t460
            t463 = t462 * t46
            t464 = t135 * t40
            t465 = t464 * t127
            t466 = t172 * t46
            t467 = t43 * t466
            t470 = t39 * t404
            t473 = t94 * t127
            t478 = 0.1e1_dp / my_rho
            t479 = t41 * t478
            t482 = t190 * t48
            t486 = t49 * t462 * t46
            t489 = t41 * t135
            t492 = t463 + 0.8999999998e1_dp * t465 * t467 - 0.5555555558e-1_dp * t470  &
            * t52 - 0.5000000001e0_dp * t473 * t466 + 0.3703703705e-1_dp * t190 * t143  &
            + 0.3333333334e0_dp * t479 * t466 - 0.5555555558e-1_dp * t482 * t148  &
            - 0.5555555558e-1_dp * t146 * t486 - 0.5000000001e0_dp * t489 * t466
            t496 = 0.1800000000e2_dp * t413 * t186 * t113 + REAL(2 * t101 * t429,KIND=dp)  &
            - 0.8333333335e-1_dp * t289 * t118 * t368 + 0.1666666667e0_dp * t117 * t435  &
            * t94 + 0.1666666667e0_dp * t117 * t118 * t404 - 0.5555555555e-1_dp  &
            * t117 * t293 * t172 - 0.1111111111e0_dp * t36 * t122 * t196 + 0.1666666667e0_dp  &
            * t117 * t297 * t172 + 0.3333333334e0_dp * t36 * t38 * t492

            e_ndrho_rho(ii) = e_ndrho_rho(ii) +  (-0.3333333336e0_dp * t36 * t173 - 0.2222222224e0_dp * t64 * t365 -  &
            0.1666666668e0_dp * t211 * t60 * t368 - 0.3333333336e0_dp * t68 * t372  &
            - 0.3333333336e0_dp * t68 * t405 - 0.3333333336e0_dp * t68 * t408 - 0.2222222224e0_dp  &
            * t24 * t98 * t496) * sx

            t501 = t172 ** 2
            t502 = t60 * t501
            t505 = t365 * t172
            t508 = beta * t10
            t513 = t168 ** 2
            t514 = t231 * t513
            t519 = t219 * my_rho
            t521 = 0.1e1_dp / t1 / t519
            t525 = 0.120e2_dp * t508 * t86 - 0.60e1_dp * t6 * t521 * t250
            t526 = t77 * t525
            t529 = 0.40e1_dp * t508 * t17 - 0.80e1_dp * t11 * t169 + 0.40e1_dp * t6 * t514  &
            - 0.20e1_dp * t6 * t526
            t530 = t60 * t529
            t533 = pi ** 2
            t534 = t533 * t100
            t536 = 0.1e1_dp / t39 / omega
            t537 = t534 * t536
            t539 = 0.1e1_dp / t22 / t311
            t562 = t8 * t501
            t566 = t8 * t529
            t570 = t39 ** 2
            t571 = 0.1e1_dp / t570
            t572 = t126 ** 2
            t573 = 0.1e1_dp / t572
            t574 = t571 * t573
            t575 = t574 * t533
            t576 = t2 * t501
            t580 = t39 * t529
            t586 = -0.2250000000e2_dp * t451 * t562 * t46 + 0.8999999998e1_dp * t185  &
            * t566 * t46 + 0.8099999996e2_dp * t575 * t576 * t46 - 0.5555555558e-1_dp  &
            * t580 * t52 - 0.5000000001e0_dp * t41 * t529 * t46
            t590 = -0.2700000000e2_dp * t537 * t539 * my_rho * t501 * t103 + 0.4500000000e1_dp  &
            * t177 * t271 * t1 * t501 - 0.3000000000e1_dp * t177 * t105 *  &
            t1 * t529 - 0.8333333335e-1_dp * t289 * t118 * t501 + 0.3333333334e0_dp  &
            * t117 * t435 * t172 + 0.1666666667e0_dp * t117 * t118 * t529 + 0.3333333334e0_dp  &
            * t36 * t38 * t586

            e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) +  (-0.1666666668e0_dp * t211 * t502 - 0.6666666672e0_dp * t68 * t505  &
            - 0.3333333336e0_dp * t68 * t530 - 0.2222222224e0_dp * t24 * t98 * t590)  &
            * sx

          END IF
          IF( grad_deriv >=3 .OR. grad_deriv == -3 ) THEN
            t601 = t27 * t358
            t605 = my_rho * t104 * omega
            t606 = t212 * t94
            t613 = t94 * t258
            t624 = 0.1e1_dp / t8 / t519
            t628 = t221 * t76
            t632 = t71 * t230
            t639 = t75 ** 2
            t640 = 0.1e1_dp / t639
            t641 = t10 * t640
            t657 = t219 ** 2
            t667 = t84 ** 2
            t669 = 0.1e1_dp / t85 / t667
            t677 = -0.9125925923e2_dp * t6 * t624 * t17 - 0.5866666667e2_dp * t6 * t628  &
            * t90 - 0.3200000001e2_dp * t6 * t632 * t232 + 0.1600000000e2_dp * t6  &
            * t225 * t254 - 0.120e2_dp * t6 * t641 * t232 * t90 + 0.120e2_dp * t382  &
            * t383 * t254 - 0.20e1_dp * t6 * t77 * (-0.6222222223e2_dp * t11 / t1 /  &
            t219 * t12 - 0.2115555556e3_dp * t6 * t624 * t86 + 0.1315555556e3_dp *  &
            t245 / t1 / t657 * t250 - 0.4266666668e2_dp * beta * t244 * t5 / t657  &
            / t69 * t669)
            t687 = t28 * t539
            t716 = t264 ** 2
            t722 = omega * t270 * t27
            t735 = t79 * t55
            t739 = t3 * t151
            t746 = t37 * t354
            t769 = t40 * t573
            t791 = 0.5400000000e2_dp * t769 * t43 * t606 - 0.3600000000e2_dp * t313 *  &
            t132 * t212 - 0.5400000000e2_dp * t451 * t452 * t258 - 0.6000000000e1_dp  &
            * t128 * t323 * t94 + 0.1800000000e2_dp * t128 * t132 * t258 + 0.8999999998e1_dp  &
            * t128 * t43 * t677 - 0.2666666667e1_dp * t42 * pi * t79
            t793 = t328 * t135
            t838 = REAL(3 * t326 * t135 * t46,KIND=dp) + REAL(t791 * t46,KIND=dp) + REAL(t793 *  &
            t46,KIND=dp) - 0.5555555558e-1_dp * t39 * t677 * t52 + 0.1111111112e0_dp * t330  &
            * t143 - 0.1666666668e0_dp * t330 * t48 * t148 - 0.1851851853e0_dp * t137  &
            * t339 + 0.2222222223e0_dp * t335 * t343 - 0.1666666668e0_dp * t335 *  &
            t347 - 0.1666666668e0_dp * t335 * t351 + 0.1646090535e0_dp * t47 * t48  &
            * t71 * t51 - 0.1851851853e0_dp * REAL(t146,KIND=dp) * REAL(t10,KIND=dp) * REAL(t135,KIND=dp)  &
            * REAL(t46,KIND=dp) + 0.1111111112e0_dp * REAL(t146,KIND=dp) * REAL(t141,KIND=dp) * REAL(t326,KIND=dp)  &
            * REAL(t46,KIND=dp) + 0.1111111112e0_dp * REAL(t146,KIND=dp) * REAL(t141,KIND=dp) * REAL(t328,KIND=dp)  &
            * REAL(t46,KIND=dp) - 0.5555555558e-1_dp * REAL(t146,KIND=dp) * REAL(t49,KIND=dp) * REAL(t791,KIND=dp)  &
            * REAL(t46,KIND=dp) - 0.1666666668e0_dp * REAL(t146,KIND=dp) * REAL(t346,KIND=dp) * REAL(t136,KIND=dp)  &
            - 0.5555555558e-1_dp * REAL(t146,KIND=dp) * REAL(t49,KIND=dp) * REAL(t793,KIND=dp) *  &
            REAL(t46,KIND=dp)
            t842 = 0.2e1_dp * t101 * (-t316 + t319 + 0.9000000000e1_dp * t321 + t325)  &
            * t103 * t112 + 0.2e1_dp * t101 * t103 * (-0.5625000000e1_dp * t687 * t31  &
            * t606 + 0.2250000000e1_dp * t271 * t109 * t212 + 0.6750000000e1_dp *  &
            t417 * t418 * t258 + 0.1000000000e1_dp * t105 * t281 * t94 - 0.1500000000e1_dp  &
            * t105 * t109 * t258 - 0.1500000000e1_dp * t105 * t31 * t677  &
            + 0.1111111111e1_dp * t30 * t26 * t10) + 0.4e1_dp * t101 * t265 * t284 +  &
            0.2e1_dp * t101 * t716 * t103 * t112 + 0.1250000000e0_dp * t722 * t118  &
            * t606 + 0.8333333333e-1_dp * t289 * t293 * t212 - 0.2500000000e0_dp * t289  &
            * t297 * t212 - 0.2500000000e0_dp * t289 * t118 * t613 + 0.2222222222e0_dp  &
            * t117 * t735 * t94 - 0.3333333333e0_dp * t117 * t739 * t94 -  &
            0.1666666667e0_dp * t117 * t293 * t258 + 0.5000000001e0_dp * t117 * t746  &
            * t94 + 0.5000000001e0_dp * t117 * t297 * t258 + 0.1666666667e0_dp * t117  &
            * t118 * t677 - 0.3456790122e0_dp * t36 * t27 * t237 * t55 + 0.4444444444e0_dp  &
            * t36 * t304 * t151 - 0.3333333333e0_dp * t36 * t122 * t354  &
            + 0.3333333334e0_dp * t36 * t38 * t838
            t846 = -0.5000000004e0_dp * t116 * t213 - 0.2000000001e1_dp * t36 * t216  &
            - 0.1000000001e1_dp * t36 * t259 - 0.6666666672e0_dp * t64 * t601 + 0.8333333340e-1_dp  &
            * t605 * t60 * t606 - 0.5000000004e0_dp * t211 * t207 * t212  &
            - 0.5000000004e0_dp * t211 * t60 * t613 - 0.1000000001e1_dp * t68 *  &
            t601 * t94 - 0.1000000001e1_dp * t68 * t207 * t258 - 0.3333333336e0_dp *  &
            t68 * t60 * t677 - 0.2222222224e0_dp * t24 * t98 * t842

            e_rho_rho_rho(ii) = e_rho_rho_rho(ii) +  t846 * sx

            t857 = t27 * t496
            t860 = t212 * t172
            t867 = t94 * t404
            t880 = t258 * t172
            t933 = 0.3911111110e2_dp * t11 * t222 - 0.1955555555e2_dp * t6 * t628 * t168  &
            + 0.2133333334e2_dp * t11 * t226 - 0.2133333334e2_dp * t6 * t71 * t384  &
            + 0.1066666667e2_dp * t6 * t225 * t400 + 0.80e1_dp * t11 * t233 - 0.120e2_dp  &
            * t382 * t640 * t232 * t168 + 0.80e1_dp * t382 * t383 * t400 - 0.40e1_dp  &
            * t11 * t255 + 0.40e1_dp * t382 * t230 * t254 * t168 - 0.20e1_dp *  &
            t6 * t77 * (0.1866666667e2_dp * beta * t237 * t12 + 0.9866666667e2_dp *  &
            t11 * t241 - 0.8266666668e2_dp * t393 * t251 + 0.3200000001e2_dp * beta  &
            * t244 * my_ndrho / t657 / t7 * t669)
            t961 = t687 * t26
            t1002 = t3 * t196
            t1009 = 0.2e1_dp * t101 * (-t455 + 0.9000000000e1_dp * t457 + 0.6000000000e1_dp  &
            * t460) * t103 * t112 + 0.1800000000e2_dp * t412 * t264 * t40 * t127  &
            * t8 * t172 * t103 * t112 + 0.2e1_dp * t101 * t265 * t428 + 0.1800000000e2_dp  &
            * t413 * t186 * t285 + 0.2e1_dp * t101 * t103 * (-0.5625000000e1_dp  &
            * t961 * t1 * t212 * t172 + 0.4500000000e1_dp * t417 * t418 * t404  &
            + 0.1500000000e1_dp * t417 * t49 * t94 * t172 - 0.1000000000e1_dp *  &
            t105 * t109 * t404 + 0.2250000000e1_dp * t417 * t1 * t258 * t172 - 0.1500000000e1_dp  &
            * t105 * t31 * t933 + 0.3333333334e0_dp * t105 * t281 *  &
            t172) + 0.1250000000e0_dp * t722 * t118 * t860 - 0.8333333335e-1_dp * t289  &
            * t435 * t212 - 0.1666666667e0_dp * t289 * t118 * t867 + 0.5555555555e-1_dp  &
            * t289 * t293 * t368 - 0.1111111111e0_dp * t117 * t1002 * t94  &
            - 0.1111111111e0_dp * t117 * t293 * t404
            t1013 = t37 * t492
            t1044 = t769 * pi
            t1069 = 0.5400000000e2_dp * t1044 * t8 * t212 * t172 - 0.3600000000e2_dp  &
            * t451 * t452 * t404 - 0.2400000000e2_dp * t451 * t37 * t94 * t172 +  &
            0.1200000000e2_dp * t128 * t132 * t404 - 0.1800000000e2_dp * t451 * t8 *  &
            t258 * t172 + 0.8999999998e1_dp * t128 * t43 * t933 - 0.2000000000e1_dp  &
            * t128 * t323 * t172
            t1091 = t127 * t172 * t46
            t1102 = t1069 * t46 + 0.8999999998e1_dp * t326 * t40 * t127 * t467 + REAL(2  &
            * t136 * t462,KIND=dp) + 0.8999999998e1_dp * t328 * t40 * t127 * t467 -  &
            0.5555555558e-1_dp * t39 * t933 * t52 - 0.5000000001e0_dp * t258 * t127  &
            * t466 + 0.7407407410e-1_dp * t470 * t143 + 0.6666666668e0_dp * t94 * t478  &
            * t1091 - 0.1111111112e0_dp * t470 * t48 * t148 - 0.1111111112e0_dp  &
            * t335 * t486 - 0.1000000001e1_dp * t94 * t135 * t1091
            t1136 = -0.6172839508e-1_dp * t190 * t339 - 0.5555555556e0_dp * t41 / t7  &
            * t466 + 0.7407407410e-1_dp * t482 * t343 + 0.7407407410e-1_dp * t146 *  &
            t141 * t462 * t46 + 0.6666666668e0_dp * t479 * t135 * t172 * t46 - 0.5555555558e-1_dp  &
            * t482 * t347 - 0.5555555558e-1_dp * t146 * t49 * t1069  &
            * t46 - 0.5000000001e0_dp * t41 * t326 * t466 - 0.5555555558e-1_dp * t482  &
            * t351 - 0.1111111112e0_dp * t146 * t147 * t463 - 0.5000000001e0_dp *  &
            t41 * t328 * t466
            t1141 = -0.1666666667e0_dp * t289 * t297 * t368 + 0.3333333334e0_dp * t117  &
            * t1013 * t94 + 0.3333333334e0_dp * t117 * t297 * t404 - 0.8333333335e-1_dp  &
            * t289 * t118 * t880 + 0.1666666667e0_dp * t117 * t435 * t258 +  &
            0.1666666667e0_dp * t117 * t118 * t933 + 0.7407407405e-1_dp * t117 * t735  &
            * t172 + 0.1481481481e0_dp * t36 * t304 * t196 - 0.1111111111e0_dp *  &
            t117 * t739 * t172 - 0.2222222222e0_dp * t36 * t122 * t492 + 0.1666666667e0_dp  &
            * t117 * t746 * t172 + 0.3333333334e0_dp * t36 * t38 * (t1102  &
            + t1136)
            t1146 = -0.3333333336e0_dp * t117 * t59 * t94 * t172 - 0.6666666672e0_dp  &
            * t36 * t372 - 0.6666666672e0_dp * t36 * t405 - 0.6666666672e0_dp * t36  &
            * t408 - 0.4444444448e0_dp * t64 * t857 + 0.8333333340e-1_dp * t605 * t60  &
            * t860 - 0.1666666668e0_dp * t211 * t365 * t212 - 0.3333333336e0_dp *  &
            t211 * t60 * t867 - 0.3333333336e0_dp * t211 * t207 * t368 - 0.6666666672e0_dp  &
            * t68 * t857 * t94 - 0.6666666672e0_dp * t68 * t207 * t404 - 0.1666666668e0_dp  &
            * t211 * t60 * t880 - 0.3333333336e0_dp * t68 * t365 *  &
            t258 - 0.3333333336e0_dp * t68 * t60 * t933 - 0.3333333336e0_dp * t68 *  &
            t601 * t172 - 0.2222222224e0_dp * t24 * t98 * (t1009 + t1141)

            e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii) +  t1146 * sx

            t1153 = t27 * t590
            t1156 = t94 * t501
            t1163 = t404 * t172
            t1167 = t94 * t529
            t1177 = beta * t71
            t1220 = -0.1066666667e2_dp * t1177 * t17 + 0.2133333334e2_dp * t11 * t377  &
            - 0.1066666667e2_dp * t6 * t632 * t513 + 0.5333333333e1_dp * t6 * t225  &
            * t525 - 0.40e1_dp * t508 * t76 * t90 + 0.160e2_dp * t11 * t10 * t384 -  &
            0.80e1_dp * t11 * t401 - 0.120e2_dp * t382 * t640 * t90 * t513 + 0.80e1_dp  &
            * t382 * t230 * t400 * t168 + 0.40e1_dp * t382 * t383 * t525 - 0.20e1_dp  &
            * t6 * t77 * (-0.3200000000e2_dp * t1177 * t86 + 0.4800000000e2_dp * t6  &
            * t397 - 0.2400000000e2_dp * t245 / t657 / my_rho * t669)
            t1284 = t37 * t586
            t1334 = 0.5400000000e2_dp * t1044 * t452 * t501 - 0.3600000000e2_dp * t451  &
            * t8 * t404 * t172 - 0.1800000000e2_dp * t451 * t452 * t529 + 0.8999999998e1_dp  &
            * t128 * t43 * t1220 - 0.1200000000e2_dp * t313 * t132 * t501  &
            + 0.5999999999e1_dp * t128 * t132 * t529
            t1341 = t501 * t46
            t1345 = t529 * t46
            t1370 = t40 * pi
            t1396 = t1334 * t46 + 0.1800000000e2_dp * t462 * t40 * t127 * t467 - 0.2250000000e2_dp  &
            * t464 * t312 * t43 * t1341 + 0.8999999998e1_dp * t465  &
            * t43 * t1345 - 0.1000000000e1_dp * t404 * t127 * t466 + 0.8099999996e2_dp  &
            * t135 * t571 * t573 * t533 * t2 * t1341 - 0.5555555558e-1_dp * t39  &
            * t1220 * t52 - 0.5000000001e0_dp * t473 * t1345 + 0.3333333334e0_dp *  &
            t479 * t1345 + 0.1000000000e1_dp * t94 * t312 * t1341 - 0.4500000000e1_dp  &
            * t94 * t573 * t501 * t1370 * t8 * t46 + 0.3703703705e-1_dp * t580  &
            * t143 + 0.3000000000e1_dp * t312 * t37 * t501 * t1370 * t46 - 0.5555555558e-1_dp  &
            * t580 * t48 * t148 - 0.1111111112e0_dp * t482 * t486 - 0.5555555558e-1_dp  &
            * t146 * t49 * t1334 * t46 - 0.1000000000e1_dp * t41 *  &
            t462 * t466 - 0.5000000001e0_dp * t489 * t1345
            t1400 = -0.3600000000e2_dp * t412 * t313 * t562 * t113 + 0.1800000000e2_dp  &
            * t413 * t566 * t113 + 0.3600000000e2_dp * t413 * t186 * t429 + 0.1620000000e3_dp  &
            * t26 * t533 * t100 * t574 * t576 * t113 + 0.2e1_dp * t101  &
            * t103 * (-0.5625000000e1_dp * t961 * t418 * t501 + 0.4500000000e1_dp  &
            * t417 * t1 * t404 * t172 + 0.2250000000e1_dp * t417 * t418 * t529 -  &
            0.1500000000e1_dp * t105 * t31 * t1220 + 0.7500000000e0_dp * t271 * t109  &
            * t501 - 0.5000000000e0_dp * t105 * t109 * t529) + 0.1250000000e0_dp *  &
            t722 * t118 * t1156 - 0.1666666667e0_dp * t289 * t435 * t368 - 0.1666666667e0_dp  &
            * t289 * t118 * t1163 - 0.8333333335e-1_dp * t289 * t118 * t1167  &
            + 0.1666666667e0_dp * t117 * t1284 * t94 + 0.3333333334e0_dp * t117  &
            * t435 * t404 + 0.1666666667e0_dp * t117 * t118 * t1220 + 0.2777777778e-1_dp  &
            * t289 * t293 * t501 - 0.1111111111e0_dp * t117 * t1002 * t172  &
            - 0.5555555555e-1_dp * t117 * t293 * t529 - 0.1111111111e0_dp * t36 * t122  &
            * t586 - 0.8333333335e-1_dp * t289 * t297 * t501 + 0.3333333334e0_dp  &
            * t117 * t1013 * t172 + 0.1666666667e0_dp * t117 * t297 * t529 + 0.3333333334e0_dp  &
            * t36 * t38 * t1396
            t1404 = -0.1666666668e0_dp * t116 * t502 - 0.6666666672e0_dp * t36 * t505  &
            - 0.3333333336e0_dp * t36 * t530 - 0.2222222224e0_dp * t64 * t1153 + 0.8333333340e-1_dp  &
            * t605 * t60 * t1156 - 0.3333333336e0_dp * t211 * t365  &
            * t368 - 0.3333333336e0_dp * t211 * t60 * t1163 - 0.1666666668e0_dp * t211  &
            * t60 * t1167 - 0.3333333336e0_dp * t68 * t1153 * t94 - 0.6666666672e0_dp  &
            * t68 * t365 * t404 - 0.3333333336e0_dp * t68 * t60 * t1220 - 0.1666666668e0_dp  &
            * t211 * t207 * t501 - 0.6666666672e0_dp * t68 * t857 *  &
            t172 - 0.3333333336e0_dp * t68 * t207 * t529 - 0.2222222224e0_dp * t24 *  &
            t98 * t1400

            e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii) +  t1404 * sx

            t1405 = t501 * t172
            t1412 = t172 * t529
            t1449 = -0.120e2_dp * t508 * t76 * t168 + 0.240e2_dp * t11 * t514 - 0.120e2_dp  &
            * t11 * t526 - 0.120e2_dp * t6 * t641 * t513 * t168 + 0.120e2_dp * t382  &
            * t230 * t168 * t525 - 0.20e1_dp * t6 * t77 * (-0.240e2_dp * beta * t521  &
            * t250 * my_ndrho + 0.180e2_dp * t393 / t657 * t669)
            t1456 = t1405 * t103
            t1467 = t533 * pi
            t1472 = t572 * t21
            t1553 = 0.1350000000e3_dp * t537 / t22 / t572 * my_rho * t1456 - 0.8100000000e2_dp  &
            * t534 * t536 * t539 * my_rho * t172 * t103 * t529 - 0.2430000000e3_dp  &
            * t1467 * t100 / t570 / omega / t22 / t1472 * t140 * t1456 -  &
            0.1125000000e2_dp * t177 * t687 * t1 * t1405 + 0.1350000000e2_dp * t176  &
            * t103 * t28 * t270 * t1 * t1412 - 0.3000000000e1_dp * t177 * t105 *  &
            t1 * t1449 + 0.1250000000e0_dp * t722 * t118 * t1405 - 0.2500000000e0_dp  &
            * t289 * t435 * t501 - 0.2500000000e0_dp * t289 * t118 * t1412 + 0.5000000001e0_dp  &
            * t117 * t1284 * t172 + 0.5000000001e0_dp * t117 * t435  &
            * t529 + 0.1666666667e0_dp * t117 * t118 * t1449 + 0.3333333334e0_dp * t36  &
            * t38 * (0.6750000000e2_dp * t1044 * t8 * t1405 * t46 - 0.6750000000e2_dp  &
            * t451 * t186 * t1345 - 0.5264999998e3_dp * t571 / t1472 * t533  &
            * t2 * t1405 * t46 + 0.8999999998e1_dp * t185 * t8 * t1449 * t46 + 0.2429999999e3_dp  &
            * t575 * t2 * t529 * t466 + 0.7289999995e3_dp / t570 / t39  &
            / t572 / t126 * t1467 * t7 * t1405 * t46 - 0.5555555558e-1_dp * t39  &
            * t1449 * t52 - 0.5000000001e0_dp * t41 * t1449 * t46)

            e_ndrho_ndrho_ndrho(ii) = e_ndrho_ndrho_ndrho(ii) +  (0.8333333340e-1_dp * t605 * t60 * t1405 - &
            0.5000000004e0_dp * t211 * t365 * t501 - 0.5000000004e0_dp * t211 * t60 * t1412 - 0.1000000001e1_dp  &
            * t68 * t1153 * t172 - 0.1000000001e1_dp * t68 * t365 * t529 - 0.3333333336e0_dp  &
            * t68 * t60 * t1449 - 0.2222222224e0_dp * t24 * t98 * t1553)  &
            * sx

          END IF
        END IF
      END DO
      !$omp end do
    END IF

  END SUBROUTINE xb88_lr_lda_calc

! *****************************************************************************
!> \brief evaluates the becke 88 longrange exchange functional for lsd
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param xb88_lr_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xb88_lr_lsd_eval(rho_set,deriv_set,grad_deriv,xb88_lr_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xb88_lr_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, npoints, &
                                                stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, omega, sx
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

    CALL section_vals_val_get(xb88_lr_params,"scale_x",r_val=sx,error=error)
    CALL section_vals_val_get(xb88_lr_params, "omega",r_val=omega, error=error)

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
          !$omp          shared(rho, norm_drho, e_0, e_rho, e_ndrho) &
          !$omp          shared(e_rho_rho, e_ndrho_rho, e_ndrho_ndrho) &
          !$omp          shared(e_rho_rho_rho, e_ndrho_rho_rho) &
          !$omp          shared(e_ndrho_ndrho_rho, e_ndrho_ndrho_ndrho) &
          !$omp          shared(grad_deriv, npoints, epsilon_rho, sx, omega) &
          !$omp          shared(error, ispin)

          CALL xb88_lr_lsd_calc(&
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
               sx=sx, omega=omega,&
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

  END SUBROUTINE xb88_lr_lsd_eval

! *****************************************************************************
!> \brief low level calculation of the becke 88 exchange functional for lsd
!> \param rho_spin alpha or beta spin density
!> \param norm_drho_spin || grad rho_spin ||
!> \param e_0 adds to it the local value of the functional
!> \param e_rho_spin e_*_spin derivative of the functional wrt. to the variables
!>        named where the * is. Everything wrt. to the spin of the arguments.
!> \param e_ndrho_spin ...
!> \param e_rho_rho_spin ...
!> \param e_ndrho_rho_spin ...
!> \param e_ndrho_ndrho_spin ...
!> \param e_rho_rho_rho_spin ...
!> \param e_ndrho_rho_rho_spin ...
!> \param e_ndrho_ndrho_rho_spin ...
!> \param e_ndrho_ndrho_ndrho_spin ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sx scaling-parameter for exchange
!> \param omega ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
!> \note
!>      - Derivatives higher than one not tested
! *****************************************************************************
  SUBROUTINE xb88_lr_lsd_calc(rho_spin, norm_drho_spin,e_0,&
       e_rho_spin,e_ndrho_spin,e_rho_rho_spin,e_ndrho_rho_spin,&
       e_ndrho_ndrho_spin, e_rho_rho_rho_spin, e_ndrho_rho_rho_spin,&
       e_ndrho_ndrho_rho_spin,&
       e_ndrho_ndrho_ndrho_spin,grad_deriv,npoints,epsilon_rho, sx,&
       omega, error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho_spin, norm_drho_spin
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_rho_spin, &
      e_ndrho_spin, e_rho_rho_spin, e_ndrho_rho_spin, e_ndrho_ndrho_spin, &
      e_rho_rho_rho_spin, e_ndrho_rho_rho_spin, e_ndrho_ndrho_rho_spin, &
      e_ndrho_ndrho_ndrho_spin
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, omega
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lr_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: Cx, epsilon_rho43, my_epsilon_rho, ndrho, rho, t1, t10, &
      t100, t1002, t1009, t101, t1013, t103, t104, t1044, t105, t1069, t109, &
      t1091, t11, t1102, t112, t113, t1136, t1141, t1146, t1153, t1156, t116, &
      t1163, t1167, t117, t1177, t118, t12, t122, t1220, t126, t127, t128, &
      t1284, t130, t132, t133, t1334, t1341, t1345, t135, t136, t137, t1370, &
      t1396, t140, t1400, t1404, t1405, t141, t1412, t143, t1449, t1456, &
      t146, t1467, t147, t1472, t148, t151, t155, t1553, t16, t168, t169, &
      t17, t172, t173, t176, t177, t18, t185, t186, t190, t196, t2, t200, &
      t207, t21, t211, t212, t213, t216, t219
    REAL(kind=dp) :: t22, t221, t222, t225, t226, t23, t230, t231, t232, &
      t233, t237, t24, t241, t244, t245, t250, t251, t254, t255, t258, t259, &
      t26, t264, t265, t27, t270, t271, t28, t281, t284, t285, t289, t29, &
      t293, t297, t3, t30, t304, t31, t311, t312, t313, t316, t319, t321, &
      t323, t325, t326, t328, t330, t335, t339, t34, t343, t346, t347, t351, &
      t354, t358, t36, t365, t368, t37, t372, t377, t38, t382, t383, t384, &
      t39, t393, t397, t40, t400, t401, t404, t405, t408, t41, t412, t413, &
      t417, t418, t42, t428, t429, t43, t435, t44, t451, t452, t455, t457, &
      t46, t460, t462, t463, t464, t465, t466
    REAL(kind=dp) :: t467, t47, t470, t473, t478, t479, t48, t482, t486, &
      t489, t49, t492, t496, t5, t501, t502, t505, t508, t51, t513, t514, &
      t519, t52, t521, t525, t526, t529, t530, t533, t534, t536, t537, t539, &
      t55, t562, t566, t570, t571, t572, t573, t574, t575, t576, t580, t586, &
      t59, t590, t6, t60, t601, t605, t606, t613, t624, t628, t632, t639, &
      t64, t640, t641, t657, t667, t669, t677, t68, t687, t69, t7, t71, t716, &
      t72, t722, t735, t739, t746, t75, t76, t769, t77, t79, t791, t793, t8, &
      t838, t84, t842, t846, t85, t857, t86, t860, t867, t87, t880, t90, t91, &
      t933, t94, t95, t961, t98, t99, xx

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
            t2 = t1 * rho
            t3 = 0.1e1_dp / t2
            xx = ndrho * MAX(t3,epsilon_rho43)
            t5 = ndrho ** 2
            t6 = beta * t5
            t7 = rho ** 2
            t8 = t1 ** 2
            t10 = 0.1e1_dp / t8 / t7
            t11 = beta * ndrho
            t12 = LOG(xx + SQRT(xx ** 0.2e1_dp + 0.1e1_dp))
            t16 = 0.10e1_dp + 0.60e1_dp * t11 * t3 * t12
            t17 = 0.1e1_dp / t16
            t18 = t10 * t17
            t21 = 0.20e1_dp * Cx + 0.20e1_dp * t6 * t18
            t22 = SQRT(t21)
            t23 = t22 * t21
            t24 = rho * t23
            t26 = rootpi
            t27 = 0.1e1_dp / t26
            t28 = 0.1e1_dp / omega
            t29 = 0.1e1_dp / t22
            t30 = t28 * t29
            t31 = t26 * t1
            t34 = erf(0.3000000000e1_dp * t30 * t31)
            t36 = omega * t22
            t37 = 0.1e1_dp / t1
            t38 = t27 * t37
            t39 = omega ** 2
            t40 = 0.1e1_dp / t39
            t41 = 0.1e1_dp / t21
            t42 = t40 * t41
            t43 = pi * t8
            t44 = t42 * t43
            t46 = EXP(-0.8999999998e1_dp * t44)
            t47 = t39 * t21
            t48 = 0.1e1_dp / pi
            t49 = 0.1e1_dp / t8
            t51 = t46 - 0.10e1_dp
            t52 = t48 * t49 * t51
            t55 = t46 - 0.15e1_dp - 0.5555555558e-1_dp * t47 * t52
            t59 = t26 * t34 + 0.3333333334e0_dp * t36 * t38 * t55
            t60 = t27 * t59

            e_0(ii) = e_0(ii) -0.2222222224e0_dp * t24 * omega * t60 * sx

          END IF
          IF( grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
            t64 = t23 * omega
            t68 = rho * t22 * omega
            t69 = t7 * rho
            t71 = 0.1e1_dp / t8 / t69
            t72 = t71 * t17
            t75 = t16 ** 2
            t76 = 0.1e1_dp / t75
            t77 = t10 * t76
            t79 = 0.1e1_dp / t1 / t7
            t84 = 1 + t5 * t10
            t85 = SQRT(t84)
            t86 = 0.1e1_dp / t85
            t87 = t71 * t86
            t90 = -0.8000000000e1_dp * t11 * t79 * t12 - 0.8000000000e1_dp * t6 * t87
            t91 = t77 * t90
            t94 = -0.5333333333e1_dp * t6 * t72 - 0.20e1_dp * t6 * t91
            t95 = t60 * t94
            t98 = omega * t27
            t99 = SQRT(0.3141592654e1_dp)
            t100 = 0.1e1_dp / t99
            t101 = t26 * t100
            t103 = EXP(-0.9000000000e1_dp * t44)
            t104 = 0.1e1_dp / t23
            t105 = t28 * t104
            t109 = t26 * t49
            t112 = -0.1500000000e1_dp * t105 * t31 * t94 + 0.1000000000e1_dp * t30 *  &
            t109
            t113 = t103 * t112
            t116 = omega * t29
            t117 = t116 * t27
            t118 = t37 * t55
            t122 = t27 * t3
            t126 = t21 ** 2
            t127 = 0.1e1_dp / t126
            t128 = t40 * t127
            t130 = t128 * t43 * t94
            t132 = pi * t37
            t133 = t42 * t132
            t135 = 0.8999999998e1_dp * t130 - 0.5999999999e1_dp * t133
            t136 = t135 * t46
            t137 = t39 * t94
            t140 = t8 * rho
            t141 = 0.1e1_dp / t140
            t143 = t48 * t141 * t51
            t146 = t47 * t48
            t147 = t49 * t135
            t148 = t147 * t46
            t151 = t136 - 0.5555555558e-1_dp * t137 * t52 + 0.3703703705e-1_dp * t47  &
            * t143 - 0.5555555558e-1_dp * t146 * t148
            t155 = REAL(2 * t101 * t113,KIND=dp) + 0.1666666667e0_dp * t117 * t118 * t94 -  &
            0.1111111111e0_dp * t36 * t122 * t55 + 0.3333333334e0_dp * t36 * t38 *  &
            t151

            e_rho_spin(ii) = e_rho_spin(ii) +  (-0.2222222224e0_dp * t64 * t60 - 0.3333333336e0_dp * t68 * t95 -  &
            0.2222222224e0_dp * t24 * t98 * t155) * sx

            t168 = 0.60e1_dp * beta * t3 * t12 + 0.60e1_dp * t11 * t10 * t86
            t169 = t77 * t168
            t172 = 0.40e1_dp * t11 * t18 - 0.20e1_dp * t6 * t169
            t173 = t60 * t172
            t176 = pi * t100
            t177 = t176 * t103
            t185 = t128 * pi
            t186 = t8 * t172
            t190 = t39 * t172
            t196 = 0.8999999998e1_dp * t185 * t186 * t46 - 0.5555555558e-1_dp * t190  &
            * t52 - 0.5000000001e0_dp * t41 * t172 * t46
            t200 = -0.3000000000e1_dp * t177 * t105 * t1 * t172 + 0.1666666667e0_dp *  &
            t117 * t118 * t172 + 0.3333333334e0_dp * t36 * t38 * t196

            e_ndrho_spin(ii) = e_ndrho_spin(ii) +  (-0.3333333336e0_dp * t68 * t173 - 0.2222222224e0_dp * t24 * t98  &
            * t200) * sx

          END IF
          IF( grad_deriv >=2 .OR. grad_deriv == -2 ) THEN
            t207 = t27 * t155
            t211 = rho * t29 * omega
            t212 = t94 ** 2
            t213 = t60 * t212
            t216 = t207 * t94
            t219 = t7 ** 2
            t221 = 0.1e1_dp / t8 / t219
            t222 = t221 * t17
            t225 = t71 * t76
            t226 = t225 * t90
            t230 = 0.1e1_dp / t75 / t16
            t231 = t10 * t230
            t232 = t90 ** 2
            t233 = t231 * t232
            t237 = 0.1e1_dp / t1 / t69
            t241 = t221 * t86
            t244 = t5 ** 2
            t245 = beta * t244
            t250 = 0.1e1_dp / t85 / t84
            t251 = 0.1e1_dp / t1 / t219 / t69 * t250
            t254 = 0.1866666667e2_dp * t11 * t237 * t12 + 0.4000000000e2_dp * t6 * t241  &
            - 0.1066666667e2_dp * t245 * t251
            t255 = t77 * t254
            t258 = 0.1955555555e2_dp * t6 * t222 + 0.1066666667e2_dp * t6 * t226 + 0.40e1_dp  &
            * t6 * t233 - 0.20e1_dp * t6 * t255
            t259 = t60 * t258
            t264 = 0.9000000000e1_dp * t130 - 0.6000000000e1_dp * t133
            t265 = t264 * t103
            t270 = 0.1e1_dp / t22 / t126
            t271 = t28 * t270
            t281 = t26 * t141
            t284 = 0.2250000000e1_dp * t271 * t31 * t212 - 0.1000000000e1_dp * t105 *  &
            t109 * t94 - 0.1500000000e1_dp * t105 * t31 * t258 - 0.6666666667e0_dp  &
            * t30 * t281
            t285 = t103 * t284
            t289 = omega * t104 * t27
            t293 = t3 * t55
            t297 = t37 * t151
            t304 = t27 * t79
            t311 = t126 * t21
            t312 = 0.1e1_dp / t311
            t313 = t40 * t312
            t316 = 0.1800000000e2_dp * t313 * t43 * t212
            t319 = 0.1200000000e2_dp * t128 * t132 * t94
            t321 = t128 * t43 * t258
            t323 = pi * t3
            t325 = 0.2000000000e1_dp * t42 * t323
            t326 = -t316 + t319 + 0.8999999998e1_dp * t321 + t325
            t328 = t135 ** 2
            t330 = t39 * t258
            t335 = t137 * t48
            t339 = t48 * t10 * t51
            t343 = t141 * t135 * t46
            t346 = t49 * t326
            t347 = t346 * t46
            t351 = t49 * t328 * t46
            t354 = t326 * t46 + t328 * t46 - 0.5555555558e-1_dp * t330 * t52 + 0.7407407410e-1_dp  &
            * t137 * t143 - 0.1111111112e0_dp * t335 * t148 - 0.6172839508e-1_dp  &
            * t47 * t339 + 0.7407407410e-1_dp * t146 * t343 - 0.5555555558e-1_dp  &
            * t146 * t347 - 0.5555555558e-1_dp * t146 * t351
            t358 = REAL(2 * t101 * t265 * t112,KIND=dp) + REAL(2 * t101 * t285,KIND=dp) - 0.8333333335e-1_dp  &
            * t289 * t118 * t212 - 0.1111111111e0_dp * t117 * t293 * t94  &
            + 0.3333333334e0_dp * t117 * t297 * t94 + 0.1666666667e0_dp * t117 *  &
            t118 * t258 + 0.1481481481e0_dp * t36 * t304 * t55 - 0.2222222222e0_dp *  &
            t36 * t122 * t151 + 0.3333333334e0_dp * t36 * t38 * t354

            e_rho_rho_spin(ii) = e_rho_rho_spin(ii) +  (-0.6666666672e0_dp * t36 * t95 - 0.4444444448e0_dp * t64 * t207 -  &
            0.1666666668e0_dp * t211 * t213 - 0.6666666672e0_dp * t68 * t216 - 0.3333333336e0_dp  &
            * t68 * t259 - 0.2222222224e0_dp * t24 * t98 * t358) * sx

            t365 = t27 * t200
            t368 = t94 * t172
            t372 = t365 * t94
            t377 = t225 * t168
            t382 = t6 * t10
            t383 = t230 * t90
            t384 = t383 * t168
            t393 = beta * t5 * ndrho
            t397 = 0.1e1_dp / t1 / t219 / t7 * t250
            t400 = -0.8000000000e1_dp * beta * t79 * t12 - 0.2400000000e2_dp * t11 *  &
            t87 + 0.8000000000e1_dp * t393 * t397
            t401 = t77 * t400
            t404 = -0.1066666667e2_dp * t11 * t72 + 0.5333333333e1_dp * t6 * t377 - 0.40e1_dp  &
            * t11 * t91 + 0.40e1_dp * t382 * t384 - 0.20e1_dp * t6 * t401
            t405 = t60 * t404
            t408 = t207 * t172
            t412 = t26 * pi * t100
            t413 = t412 * t128
            t417 = t271 * t26
            t418 = t1 * t94
            t428 = 0.2250000000e1_dp * t417 * t418 * t172 - 0.1500000000e1_dp * t105  &
            * t31 * t404 - 0.5000000000e0_dp * t105 * t109 * t172
            t429 = t103 * t428
            t435 = t37 * t196
            t451 = t313 * pi
            t452 = t8 * t94
            t455 = 0.1800000000e2_dp * t451 * t452 * t172
            t457 = t128 * t43 * t404
            t460 = t128 * t132 * t172
            t462 = -t455 + 0.8999999998e1_dp * t457 + 0.5999999999e1_dp * t460
            t463 = t462 * t46
            t464 = t135 * t40
            t465 = t464 * t127
            t466 = t172 * t46
            t467 = t43 * t466
            t470 = t39 * t404
            t473 = t94 * t127
            t478 = 0.1e1_dp / rho
            t479 = t41 * t478
            t482 = t190 * t48
            t486 = t49 * t462 * t46
            t489 = t41 * t135
            t492 = t463 + 0.8999999998e1_dp * t465 * t467 - 0.5555555558e-1_dp * t470  &
            * t52 - 0.5000000001e0_dp * t473 * t466 + 0.3703703705e-1_dp * t190 * t143  &
            + 0.3333333334e0_dp * t479 * t466 - 0.5555555558e-1_dp * t482 * t148  &
            - 0.5555555558e-1_dp * t146 * t486 - 0.5000000001e0_dp * t489 * t466
            t496 = 0.1800000000e2_dp * t413 * t186 * t113 + REAL(2 * t101 * t429,KIND=dp)  &
            - 0.8333333335e-1_dp * t289 * t118 * t368 + 0.1666666667e0_dp * t117 * t435  &
            * t94 + 0.1666666667e0_dp * t117 * t118 * t404 - 0.5555555555e-1_dp  &
            * t117 * t293 * t172 - 0.1111111111e0_dp * t36 * t122 * t196 + 0.1666666667e0_dp  &
            * t117 * t297 * t172 + 0.3333333334e0_dp * t36 * t38 * t492

            e_ndrho_rho_spin(ii) = e_ndrho_rho_spin(ii) +  (-0.3333333336e0_dp * t36 * t173 - 0.2222222224e0_dp * t64 * t365 -  &
            0.1666666668e0_dp * t211 * t60 * t368 - 0.3333333336e0_dp * t68 * t372  &
            - 0.3333333336e0_dp * t68 * t405 - 0.3333333336e0_dp * t68 * t408 - 0.2222222224e0_dp  &
            * t24 * t98 * t496) * sx

            t501 = t172 ** 2
            t502 = t60 * t501
            t505 = t365 * t172
            t508 = beta * t10
            t513 = t168 ** 2
            t514 = t231 * t513
            t519 = t219 * rho
            t521 = 0.1e1_dp / t1 / t519
            t525 = 0.120e2_dp * t508 * t86 - 0.60e1_dp * t6 * t521 * t250
            t526 = t77 * t525
            t529 = 0.40e1_dp * t508 * t17 - 0.80e1_dp * t11 * t169 + 0.40e1_dp * t6 * t514  &
            - 0.20e1_dp * t6 * t526
            t530 = t60 * t529
            t533 = pi ** 2
            t534 = t533 * t100
            t536 = 0.1e1_dp / t39 / omega
            t537 = t534 * t536
            t539 = 0.1e1_dp / t22 / t311
            t562 = t8 * t501
            t566 = t8 * t529
            t570 = t39 ** 2
            t571 = 0.1e1_dp / t570
            t572 = t126 ** 2
            t573 = 0.1e1_dp / t572
            t574 = t571 * t573
            t575 = t574 * t533
            t576 = t2 * t501
            t580 = t39 * t529
            t586 = -0.2250000000e2_dp * t451 * t562 * t46 + 0.8999999998e1_dp * t185  &
            * t566 * t46 + 0.8099999996e2_dp * t575 * t576 * t46 - 0.5555555558e-1_dp  &
            * t580 * t52 - 0.5000000001e0_dp * t41 * t529 * t46
            t590 = -0.2700000000e2_dp * t537 * t539 * rho * t501 * t103 + 0.4500000000e1_dp  &
            * t177 * t271 * t1 * t501 - 0.3000000000e1_dp * t177 * t105 *  &
            t1 * t529 - 0.8333333335e-1_dp * t289 * t118 * t501 + 0.3333333334e0_dp  &
            * t117 * t435 * t172 + 0.1666666667e0_dp * t117 * t118 * t529 + 0.3333333334e0_dp  &
            * t36 * t38 * t586

            e_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_spin(ii) +  (-0.1666666668e0_dp * t211 * t502 - 0.6666666672e0_dp * t68 * t505  &
            - 0.3333333336e0_dp * t68 * t530 - 0.2222222224e0_dp * t24 * t98 * t590)  &
            * sx

          END IF
          IF( grad_deriv >=3 .OR. grad_deriv == -3 ) THEN
            t601 = t27 * t358
            t605 = rho * t104 * omega
            t606 = t212 * t94
            t613 = t94 * t258
            t624 = 0.1e1_dp / t8 / t519
            t628 = t221 * t76
            t632 = t71 * t230
            t639 = t75 ** 2
            t640 = 0.1e1_dp / t639
            t641 = t10 * t640
            t657 = t219 ** 2
            t667 = t84 ** 2
            t669 = 0.1e1_dp / t85 / t667
            t677 = -0.9125925923e2_dp * t6 * t624 * t17 - 0.5866666667e2_dp * t6 * t628  &
            * t90 - 0.3200000001e2_dp * t6 * t632 * t232 + 0.1600000000e2_dp * t6  &
            * t225 * t254 - 0.120e2_dp * t6 * t641 * t232 * t90 + 0.120e2_dp * t382  &
            * t383 * t254 - 0.20e1_dp * t6 * t77 * (-0.6222222223e2_dp * t11 / t1 /  &
            t219 * t12 - 0.2115555556e3_dp * t6 * t624 * t86 + 0.1315555556e3_dp *  &
            t245 / t1 / t657 * t250 - 0.4266666668e2_dp * beta * t244 * t5 / t657  &
            / t69 * t669)
            t687 = t28 * t539
            t716 = t264 ** 2
            t722 = omega * t270 * t27
            t735 = t79 * t55
            t739 = t3 * t151
            t746 = t37 * t354
            t769 = t40 * t573
            t791 = 0.5400000000e2_dp * t769 * t43 * t606 - 0.3600000000e2_dp * t313 *  &
            t132 * t212 - 0.5400000000e2_dp * t451 * t452 * t258 - 0.6000000000e1_dp  &
            * t128 * t323 * t94 + 0.1800000000e2_dp * t128 * t132 * t258 + 0.8999999998e1_dp  &
            * t128 * t43 * t677 - 0.2666666667e1_dp * t42 * pi * t79
            t793 = t328 * t135
            t838 = REAL(3 * t326 * t135 * t46,KIND=dp) + REAL(t791 * t46,KIND=dp) + REAL(t793 *  &
            t46,KIND=dp) - 0.5555555558e-1_dp * t39 * t677 * t52 + 0.1111111112e0_dp * t330  &
            * t143 - 0.1666666668e0_dp * t330 * t48 * t148 - 0.1851851853e0_dp * t137  &
            * t339 + 0.2222222223e0_dp * t335 * t343 - 0.1666666668e0_dp * t335 *  &
            t347 - 0.1666666668e0_dp * t335 * t351 + 0.1646090535e0_dp * t47 * t48  &
            * t71 * t51 - 0.1851851853e0_dp * REAL(t146,KIND=dp) * REAL(t10,KIND=dp) * REAL(t135,KIND=dp)  &
            * REAL(t46,KIND=dp) + 0.1111111112e0_dp * REAL(t146,KIND=dp) * REAL(t141,KIND=dp) * REAL(t326,KIND=dp)  &
            * REAL(t46,KIND=dp) + 0.1111111112e0_dp * REAL(t146,KIND=dp) * REAL(t141,KIND=dp) * REAL(t328,KIND=dp)  &
            * REAL(t46,KIND=dp) - 0.5555555558e-1_dp * REAL(t146,KIND=dp) * REAL(t49,KIND=dp) * REAL(t791,KIND=dp)  &
            * REAL(t46,KIND=dp) - 0.1666666668e0_dp * REAL(t146,KIND=dp) * REAL(t346,KIND=dp) * REAL(t136,KIND=dp)  &
            - 0.5555555558e-1_dp * REAL(t146,KIND=dp) * REAL(t49,KIND=dp) * REAL(t793,KIND=dp) *  &
            REAL(t46,KIND=dp)
            t842 = 0.2e1_dp * t101 * (-t316 + t319 + 0.9000000000e1_dp * t321 + t325)  &
            * t103 * t112 + 0.2e1_dp * t101 * t103 * (-0.5625000000e1_dp * t687 * t31  &
            * t606 + 0.2250000000e1_dp * t271 * t109 * t212 + 0.6750000000e1_dp *  &
            t417 * t418 * t258 + 0.1000000000e1_dp * t105 * t281 * t94 - 0.1500000000e1_dp  &
            * t105 * t109 * t258 - 0.1500000000e1_dp * t105 * t31 * t677  &
            + 0.1111111111e1_dp * t30 * t26 * t10) + 0.4e1_dp * t101 * t265 * t284 +  &
            0.2e1_dp * t101 * t716 * t103 * t112 + 0.1250000000e0_dp * t722 * t118  &
            * t606 + 0.8333333333e-1_dp * t289 * t293 * t212 - 0.2500000000e0_dp * t289  &
            * t297 * t212 - 0.2500000000e0_dp * t289 * t118 * t613 + 0.2222222222e0_dp  &
            * t117 * t735 * t94 - 0.3333333333e0_dp * t117 * t739 * t94 -  &
            0.1666666667e0_dp * t117 * t293 * t258 + 0.5000000001e0_dp * t117 * t746  &
            * t94 + 0.5000000001e0_dp * t117 * t297 * t258 + 0.1666666667e0_dp * t117  &
            * t118 * t677 - 0.3456790122e0_dp * t36 * t27 * t237 * t55 + 0.4444444444e0_dp  &
            * t36 * t304 * t151 - 0.3333333333e0_dp * t36 * t122 * t354  &
            + 0.3333333334e0_dp * t36 * t38 * t838
            t846 = -0.5000000004e0_dp * t116 * t213 - 0.2000000001e1_dp * t36 * t216  &
            - 0.1000000001e1_dp * t36 * t259 - 0.6666666672e0_dp * t64 * t601 + 0.8333333340e-1_dp  &
            * t605 * t60 * t606 - 0.5000000004e0_dp * t211 * t207 * t212  &
            - 0.5000000004e0_dp * t211 * t60 * t613 - 0.1000000001e1_dp * t68 *  &
            t601 * t94 - 0.1000000001e1_dp * t68 * t207 * t258 - 0.3333333336e0_dp *  &
            t68 * t60 * t677 - 0.2222222224e0_dp * t24 * t98 * t842

            e_rho_rho_rho_spin(ii) = e_rho_rho_rho_spin(ii) +  t846 * sx

            t857 = t27 * t496
            t860 = t212 * t172
            t867 = t94 * t404
            t880 = t258 * t172
            t933 = 0.3911111110e2_dp * t11 * t222 - 0.1955555555e2_dp * t6 * t628 * t168  &
            + 0.2133333334e2_dp * t11 * t226 - 0.2133333334e2_dp * t6 * t71 * t384  &
            + 0.1066666667e2_dp * t6 * t225 * t400 + 0.80e1_dp * t11 * t233 - 0.120e2_dp  &
            * t382 * t640 * t232 * t168 + 0.80e1_dp * t382 * t383 * t400 - 0.40e1_dp  &
            * t11 * t255 + 0.40e1_dp * t382 * t230 * t254 * t168 - 0.20e1_dp *  &
            t6 * t77 * (0.1866666667e2_dp * beta * t237 * t12 + 0.9866666667e2_dp *  &
            t11 * t241 - 0.8266666668e2_dp * t393 * t251 + 0.3200000001e2_dp * beta  &
            * t244 * ndrho / t657 / t7 * t669)
            t961 = t687 * t26
            t1002 = t3 * t196
            t1009 = 0.2e1_dp * t101 * (-t455 + 0.9000000000e1_dp * t457 + 0.6000000000e1_dp  &
            * t460) * t103 * t112 + 0.1800000000e2_dp * t412 * t264 * t40 * t127  &
            * t8 * t172 * t103 * t112 + 0.2e1_dp * t101 * t265 * t428 + 0.1800000000e2_dp  &
            * t413 * t186 * t285 + 0.2e1_dp * t101 * t103 * (-0.5625000000e1_dp  &
            * t961 * t1 * t212 * t172 + 0.4500000000e1_dp * t417 * t418 * t404  &
            + 0.1500000000e1_dp * t417 * t49 * t94 * t172 - 0.1000000000e1_dp *  &
            t105 * t109 * t404 + 0.2250000000e1_dp * t417 * t1 * t258 * t172 - 0.1500000000e1_dp  &
            * t105 * t31 * t933 + 0.3333333334e0_dp * t105 * t281 *  &
            t172) + 0.1250000000e0_dp * t722 * t118 * t860 - 0.8333333335e-1_dp * t289  &
            * t435 * t212 - 0.1666666667e0_dp * t289 * t118 * t867 + 0.5555555555e-1_dp  &
            * t289 * t293 * t368 - 0.1111111111e0_dp * t117 * t1002 * t94  &
            - 0.1111111111e0_dp * t117 * t293 * t404
            t1013 = t37 * t492
            t1044 = t769 * pi
            t1069 = 0.5400000000e2_dp * t1044 * t8 * t212 * t172 - 0.3600000000e2_dp  &
            * t451 * t452 * t404 - 0.2400000000e2_dp * t451 * t37 * t94 * t172 +  &
            0.1200000000e2_dp * t128 * t132 * t404 - 0.1800000000e2_dp * t451 * t8 *  &
            t258 * t172 + 0.8999999998e1_dp * t128 * t43 * t933 - 0.2000000000e1_dp  &
            * t128 * t323 * t172
            t1091 = t127 * t172 * t46
            t1102 = t1069 * t46 + 0.8999999998e1_dp * t326 * t40 * t127 * t467 + REAL(2  &
            * t136 * t462,KIND=dp) + 0.8999999998e1_dp * t328 * t40 * t127 * t467 -  &
            0.5555555558e-1_dp * t39 * t933 * t52 - 0.5000000001e0_dp * t258 * t127  &
            * t466 + 0.7407407410e-1_dp * t470 * t143 + 0.6666666668e0_dp * t94 * t478  &
            * t1091 - 0.1111111112e0_dp * t470 * t48 * t148 - 0.1111111112e0_dp  &
            * t335 * t486 - 0.1000000001e1_dp * t94 * t135 * t1091
            t1136 = -0.6172839508e-1_dp * t190 * t339 - 0.5555555556e0_dp * t41 / t7  &
            * t466 + 0.7407407410e-1_dp * t482 * t343 + 0.7407407410e-1_dp * t146 *  &
            t141 * t462 * t46 + 0.6666666668e0_dp * t479 * t135 * t172 * t46 - 0.5555555558e-1_dp  &
            * t482 * t347 - 0.5555555558e-1_dp * t146 * t49 * t1069  &
            * t46 - 0.5000000001e0_dp * t41 * t326 * t466 - 0.5555555558e-1_dp * t482  &
            * t351 - 0.1111111112e0_dp * t146 * t147 * t463 - 0.5000000001e0_dp *  &
            t41 * t328 * t466
            t1141 = -0.1666666667e0_dp * t289 * t297 * t368 + 0.3333333334e0_dp * t117  &
            * t1013 * t94 + 0.3333333334e0_dp * t117 * t297 * t404 - 0.8333333335e-1_dp  &
            * t289 * t118 * t880 + 0.1666666667e0_dp * t117 * t435 * t258 +  &
            0.1666666667e0_dp * t117 * t118 * t933 + 0.7407407405e-1_dp * t117 * t735  &
            * t172 + 0.1481481481e0_dp * t36 * t304 * t196 - 0.1111111111e0_dp *  &
            t117 * t739 * t172 - 0.2222222222e0_dp * t36 * t122 * t492 + 0.1666666667e0_dp  &
            * t117 * t746 * t172 + 0.3333333334e0_dp * t36 * t38 * (t1102  &
            + t1136)
            t1146 = -0.3333333336e0_dp * t117 * t59 * t94 * t172 - 0.6666666672e0_dp  &
            * t36 * t372 - 0.6666666672e0_dp * t36 * t405 - 0.6666666672e0_dp * t36  &
            * t408 - 0.4444444448e0_dp * t64 * t857 + 0.8333333340e-1_dp * t605 * t60  &
            * t860 - 0.1666666668e0_dp * t211 * t365 * t212 - 0.3333333336e0_dp *  &
            t211 * t60 * t867 - 0.3333333336e0_dp * t211 * t207 * t368 - 0.6666666672e0_dp  &
            * t68 * t857 * t94 - 0.6666666672e0_dp * t68 * t207 * t404 - 0.1666666668e0_dp  &
            * t211 * t60 * t880 - 0.3333333336e0_dp * t68 * t365 *  &
            t258 - 0.3333333336e0_dp * t68 * t60 * t933 - 0.3333333336e0_dp * t68 *  &
            t601 * t172 - 0.2222222224e0_dp * t24 * t98 * (t1009 + t1141)

            e_ndrho_rho_rho_spin(ii) = e_ndrho_rho_rho_spin(ii) +  t1146 * sx

            t1153 = t27 * t590
            t1156 = t94 * t501
            t1163 = t404 * t172
            t1167 = t94 * t529
            t1177 = beta * t71
            t1220 = -0.1066666667e2_dp * t1177 * t17 + 0.2133333334e2_dp * t11 * t377  &
            - 0.1066666667e2_dp * t6 * t632 * t513 + 0.5333333333e1_dp * t6 * t225  &
            * t525 - 0.40e1_dp * t508 * t76 * t90 + 0.160e2_dp * t11 * t10 * t384 -  &
            0.80e1_dp * t11 * t401 - 0.120e2_dp * t382 * t640 * t90 * t513 + 0.80e1_dp  &
            * t382 * t230 * t400 * t168 + 0.40e1_dp * t382 * t383 * t525 - 0.20e1_dp  &
            * t6 * t77 * (-0.3200000000e2_dp * t1177 * t86 + 0.4800000000e2_dp * t6  &
            * t397 - 0.2400000000e2_dp * t245 / t657 / rho * t669)
            t1284 = t37 * t586
            t1334 = 0.5400000000e2_dp * t1044 * t452 * t501 - 0.3600000000e2_dp * t451  &
            * t8 * t404 * t172 - 0.1800000000e2_dp * t451 * t452 * t529 + 0.8999999998e1_dp  &
            * t128 * t43 * t1220 - 0.1200000000e2_dp * t313 * t132 * t501  &
            + 0.5999999999e1_dp * t128 * t132 * t529
            t1341 = t501 * t46
            t1345 = t529 * t46
            t1370 = t40 * pi
            t1396 = t1334 * t46 + 0.1800000000e2_dp * t462 * t40 * t127 * t467 - 0.2250000000e2_dp  &
            * t464 * t312 * t43 * t1341 + 0.8999999998e1_dp * t465  &
            * t43 * t1345 - 0.1000000000e1_dp * t404 * t127 * t466 + 0.8099999996e2_dp  &
            * t135 * t571 * t573 * t533 * t2 * t1341 - 0.5555555558e-1_dp * t39  &
            * t1220 * t52 - 0.5000000001e0_dp * t473 * t1345 + 0.3333333334e0_dp *  &
            t479 * t1345 + 0.1000000000e1_dp * t94 * t312 * t1341 - 0.4500000000e1_dp  &
            * t94 * t573 * t501 * t1370 * t8 * t46 + 0.3703703705e-1_dp * t580  &
            * t143 + 0.3000000000e1_dp * t312 * t37 * t501 * t1370 * t46 - 0.5555555558e-1_dp  &
            * t580 * t48 * t148 - 0.1111111112e0_dp * t482 * t486 - 0.5555555558e-1_dp  &
            * t146 * t49 * t1334 * t46 - 0.1000000000e1_dp * t41 *  &
            t462 * t466 - 0.5000000001e0_dp * t489 * t1345
            t1400 = -0.3600000000e2_dp * t412 * t313 * t562 * t113 + 0.1800000000e2_dp  &
            * t413 * t566 * t113 + 0.3600000000e2_dp * t413 * t186 * t429 + 0.1620000000e3_dp  &
            * t26 * t533 * t100 * t574 * t576 * t113 + 0.2e1_dp * t101  &
            * t103 * (-0.5625000000e1_dp * t961 * t418 * t501 + 0.4500000000e1_dp  &
            * t417 * t1 * t404 * t172 + 0.2250000000e1_dp * t417 * t418 * t529 -  &
            0.1500000000e1_dp * t105 * t31 * t1220 + 0.7500000000e0_dp * t271 * t109  &
            * t501 - 0.5000000000e0_dp * t105 * t109 * t529) + 0.1250000000e0_dp *  &
            t722 * t118 * t1156 - 0.1666666667e0_dp * t289 * t435 * t368 - 0.1666666667e0_dp  &
            * t289 * t118 * t1163 - 0.8333333335e-1_dp * t289 * t118 * t1167  &
            + 0.1666666667e0_dp * t117 * t1284 * t94 + 0.3333333334e0_dp * t117  &
            * t435 * t404 + 0.1666666667e0_dp * t117 * t118 * t1220 + 0.2777777778e-1_dp  &
            * t289 * t293 * t501 - 0.1111111111e0_dp * t117 * t1002 * t172  &
            - 0.5555555555e-1_dp * t117 * t293 * t529 - 0.1111111111e0_dp * t36 * t122  &
            * t586 - 0.8333333335e-1_dp * t289 * t297 * t501 + 0.3333333334e0_dp  &
            * t117 * t1013 * t172 + 0.1666666667e0_dp * t117 * t297 * t529 + 0.3333333334e0_dp  &
            * t36 * t38 * t1396
            t1404 = -0.1666666668e0_dp * t116 * t502 - 0.6666666672e0_dp * t36 * t505  &
            - 0.3333333336e0_dp * t36 * t530 - 0.2222222224e0_dp * t64 * t1153 + 0.8333333340e-1_dp  &
            * t605 * t60 * t1156 - 0.3333333336e0_dp * t211 * t365  &
            * t368 - 0.3333333336e0_dp * t211 * t60 * t1163 - 0.1666666668e0_dp * t211  &
            * t60 * t1167 - 0.3333333336e0_dp * t68 * t1153 * t94 - 0.6666666672e0_dp  &
            * t68 * t365 * t404 - 0.3333333336e0_dp * t68 * t60 * t1220 - 0.1666666668e0_dp  &
            * t211 * t207 * t501 - 0.6666666672e0_dp * t68 * t857 *  &
            t172 - 0.3333333336e0_dp * t68 * t207 * t529 - 0.2222222224e0_dp * t24 *  &
            t98 * t1400

            e_ndrho_ndrho_rho_spin(ii) = e_ndrho_ndrho_rho_spin(ii) +  t1404 * sx

            t1405 = t501 * t172
            t1412 = t172 * t529
            t1449 = -0.120e2_dp * t508 * t76 * t168 + 0.240e2_dp * t11 * t514 - 0.120e2_dp  &
            * t11 * t526 - 0.120e2_dp * t6 * t641 * t513 * t168 + 0.120e2_dp * t382  &
            * t230 * t168 * t525 - 0.20e1_dp * t6 * t77 * (-0.240e2_dp * beta * t521  &
            * t250 * ndrho + 0.180e2_dp * t393 / t657 * t669)
            t1456 = t1405 * t103
            t1467 = t533 * pi
            t1472 = t572 * t21
            t1553 = 0.1350000000e3_dp * t537 / t22 / t572 * rho * t1456 - 0.8100000000e2_dp  &
            * t534 * t536 * t539 * rho * t172 * t103 * t529 - 0.2430000000e3_dp  &
            * t1467 * t100 / t570 / omega / t22 / t1472 * t140 * t1456 -  &
            0.1125000000e2_dp * t177 * t687 * t1 * t1405 + 0.1350000000e2_dp * t176  &
            * t103 * t28 * t270 * t1 * t1412 - 0.3000000000e1_dp * t177 * t105 *  &
            t1 * t1449 + 0.1250000000e0_dp * t722 * t118 * t1405 - 0.2500000000e0_dp  &
            * t289 * t435 * t501 - 0.2500000000e0_dp * t289 * t118 * t1412 + 0.5000000001e0_dp  &
            * t117 * t1284 * t172 + 0.5000000001e0_dp * t117 * t435  &
            * t529 + 0.1666666667e0_dp * t117 * t118 * t1449 + 0.3333333334e0_dp * t36  &
            * t38 * (0.6750000000e2_dp * t1044 * t8 * t1405 * t46 - 0.6750000000e2_dp  &
            * t451 * t186 * t1345 - 0.5264999998e3_dp * t571 / t1472 * t533  &
            * t2 * t1405 * t46 + 0.8999999998e1_dp * t185 * t8 * t1449 * t46 + 0.2429999999e3_dp  &
            * t575 * t2 * t529 * t466 + 0.7289999995e3_dp / t570 / t39  &
            / t572 / t126 * t1467 * t7 * t1405 * t46 - 0.5555555558e-1_dp * t39  &
            * t1449 * t52 - 0.5000000001e0_dp * t41 * t1449 * t46)

            e_ndrho_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_ndrho_spin(ii) +  (0.8333333340e-1_dp * t605 * t60 * t1405 - &
            0.5000000004e0_dp * t211 * t365 * t501 - 0.5000000004e0_dp * t211 * t60 * t1412 - 0.1000000001e1_dp  &
            * t68 * t1153 * t172 - 0.1000000001e1_dp * t68 * t365 * t529 - 0.3333333336e0_dp  &
            * t68 * t60 * t1449 - 0.2222222224e0_dp * t24 * t98 * t1553)  &
            * sx

          END IF
        END IF
      END DO
      !$omp end do
    END IF

  END SUBROUTINE xb88_lr_lsd_calc

END MODULE xc_xbecke88_long_range
