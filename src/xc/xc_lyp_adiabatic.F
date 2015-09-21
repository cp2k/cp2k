!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the density_scaled Lyp functional when used in adiabatic hybrids.
!>      The energy is given as
!>
!>        Ec = 2*lambda*Ec(rho/lambda) + lambda^2*d/dlambda(Ec(rho/lambda)),
!>
!>      where rho/lambda is the scaled density
!> \par History
!>      1.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE xc_lyp_adiabatic
  USE bibliography,                    ONLY: Lee1988,&
                                             cite_reference
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
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
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_lyp_adiabatic'
  REAL(kind=dp), PARAMETER, PRIVATE :: a=0.04918_dp, b=0.132_dp,&
       c=0.2533_dp,d=0.349_dp

  PUBLIC :: lyp_adiabatic_lda_info, lyp_adiabatic_lsd_info, lyp_adiabatic_lda_eval, lyp_adiabatic_lsd_eval

CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE lyp_adiabatic_lda_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lda_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "C. Lee, W. Yang, R.G. Parr, Phys. Rev. B, 37, 785 (1988) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Lee-Yang-Parr correlation energy functional (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%rho_1_3=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE lyp_adiabatic_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE lyp_adiabatic_lsd_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lsd_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
      reference = "C. Lee, W. Yang, R.G. Parr, Phys. Rev. B, 37, 785 (1988) {LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Lee-Yang-Parr correlation energy functional (LSD)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%norm_drho_spin=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1
  END SUBROUTINE lyp_adiabatic_lsd_info

! *****************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param lyp_adiabatic_params ...
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE lyp_adiabatic_lda_eval(rho_set,deriv_set,grad_deriv,lyp_adiabatic_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: lyp_adiabatic_params

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, lambda
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrho, e_rho, &
                                                norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(lyp_adiabatic_params,"LAMBDA",r_val=lambda)
    CALL cite_reference(Lee1988)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
         norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
         drho_cutoff=epsilon_norm_drho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
    ELSE
      dummy=> rho
    END IF

    e_0 => dummy
    e_rho => dummy
    e_ndrho => dummy

    IF (grad_deriv>=0) THEN
      deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rho)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
    END IF
    IF (grad_deriv>1.OR.grad_deriv<-1) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
           message="derivatives bigger than 1 not implemented", &
           error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp          shared(rho, norm_drho, e_0, e_rho, e_ndrho) &
    !$omp          shared(grad_deriv, npoints) &
    !$omp          shared(epsilon_rho, lambda)

    CALL lyp_adiabatic_lda_calc(rho=rho, norm_drho=norm_drho,&
         e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
         grad_deriv=grad_deriv,&
         npoints=npoints,epsilon_rho=epsilon_rho,lambda=lambda)

    !$omp end parallel

    IF (cp_debug) THEN
      DEALLOCATE(dummy)
    ELSE
      NULLIFY(dummy)
    END IF
    CALL timestop(handle)
  END SUBROUTINE lyp_adiabatic_lda_eval

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param lambda ...
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE lyp_adiabatic_lda_calc(rho, norm_drho,&
     e_0,e_rho,e_ndrho, &
     grad_deriv,npoints,epsilon_rho, lambda)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, lambda

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    REAL(kind=dp) :: cf, my_ndrho, my_rho, t10, t107, t11, t117, t12, t122, &
      t125, t13, t14, t15, t153, t16, t17, t180, t189, t19, t195, t2, t20, &
      t25, t28, t29, t3, t34, t36, t37, t38, t4, t40, t41, t42, t43, t45, &
      t46, t47, t50, t51, t52, t57, t58, t59, t6, t63, t65, t7, t71, t77, &
      t78, t87, t9, t94

    cf=0.3_dp*(3._dp*pi*pi)**(2._dp/3._dp)

    !$omp do

    DO ii=1,npoints
      my_rho = rho(ii)
      IF (my_rho>epsilon_rho) THEN
        IF( grad_deriv >=0 )THEN
          my_ndrho = norm_drho(ii)
          t2 = d * lambda
          t3 = my_rho ** (0.1e1_dp / 0.3e1_dp)
          t4 = 0.1e1_dp / t3
          t6 = 0.10e1_dp + t2 * t4
          t7 = 0.1e1_dp / t6
          t9 = a * b
          t10 = t9 * my_rho
          t11 = c * lambda
          t12 = t11 * t4
          t13 = EXP(-t12)
          t14 = t13 * t7
          t15 = my_ndrho ** 2
          t16 = my_rho ** 2
          t17 = t3 ** 2
          t19 = 0.1e1_dp / t17 / t16
          t20 = t15 * t19
          t25 = 0.30e1_dp + 0.70e1_dp * t12 + 0.70e1_dp * t2 * t4 * t7
          t28 = Cf - 0.1388888889e-1_dp * t20 * t25
          t29 = t14 * t28
          t34 = lambda ** 2
          t36 = t6 ** 2
          t37 = 0.1e1_dp / t36
          t38 = t37 * d
          t40 = t9 * t17
          t41 = c * t13
          t42 = t7 * t28
          t43 = t41 * t42
          t45 = t13 * t37
          t46 = t28 * d
          t47 = t45 * t46
          t50 = 0.1e1_dp / t17 / my_rho
          t51 = t9 * t50
          t52 = c * t4
          t57 = d ** 2
          t58 = t57 * lambda
          t59 = 0.1e1_dp / t17
          t63 = 0.70e1_dp * t52 + 0.70e1_dp * d * t4 * t7 - 0.70e1_dp * t58 * t59 * t37
          t65 = t14 * t15 * t63

          e_0(ii) = e_0(ii) +  0.20e1_dp * lambda * (-a * my_rho * t7 - t10 * t29) + t34 * (a * t17  &
          * t38 + t40 * t43 + t40 * t47 + 0.13888888888888888889e-1_dp * t51 *  &
          t65)

        END IF
        IF(grad_deriv >=1 ) THEN
          t71 = a * t4
          t77 = lambda * t13
          t78 = t77 * t42
          t87 = t16 * my_rho
          t94 = 0.1e1_dp / t3 / my_rho
          t107 = 0.37037037037037037037e-1_dp * t15 / t17 / t87 * t25 - 0.1388888889e-1_dp  &
          * t20 * (-0.2333333333e1_dp * t11 * t94 - 0.2333333333e1_dp * t2  &
          * t94 * t7 + 0.23333333333333333333e1_dp * t57 * t34 * t50 * t37)
          t117 = 0.1e1_dp / t36 / t6
          t122 = t9 * t4
          t125 = c ** 2
          t153 = 0.1e1_dp / t87
          t180 = 0.1e1_dp / t16
          t189 = 0.2e1_dp / 0.3e1_dp * t71 * t38 + 0.2e1_dp / 0.3e1_dp * a * t59 * t117 *  &
          t57 * lambda + 0.2e1_dp / 0.3e1_dp * t122 * t43 + t9 * t59 * t125 * t78  &
          / 0.3e1_dp + 0.2e1_dp / 0.3e1_dp * t9 * t59 * c * t45 * t46 * lambda + t40  &
          * t41 * t7 * t107 + 0.2e1_dp / 0.3e1_dp * t122 * t47 + 0.2e1_dp / 0.3e1_dp *  &
          t9 * t59 * t13 * t117 * t28 * t58 + t40 * t45 * t107 * d - 0.2314814815e-1_dp  &
          * t9 * t19 * t65 + 0.46296296296296296297e-2_dp * t9 * t153  &
          * c * t77 * t7 * t15 * t63 + 0.46296296296296296297e-2_dp * t9 * t153  &
          * t13 * t37 * t15 * t63 * d * lambda + 0.13888888888888888889e-1_dp  &
          * t51 * t14 * t15 * (-0.2333333333e1_dp * c * t94 - 0.2333333333e1_dp *  &
          d * t94 * t7 + 0.70000000000000000000e1_dp * t57 * t50 * t37 * lambda  &
          - 0.4666666667e1_dp * t57 * d * t34 * t180 * t117)

          e_rho(ii) = e_rho(ii) +  0.20e1_dp * lambda * (-a * t7 - t71 * t38 * lambda / 0.3e1_dp - t9 *  &
          t29 - t9 * t52 * t78 / 0.3e1_dp - t9 * t4 * t13 * t37 * t28 * t2 / 0.3e1_dp  &
          - t10 * t14 * t107) + t34 * t189
          t195 = t14 * my_ndrho * t25

          e_ndrho(ii) = e_ndrho(ii) +  0.55555555555555555556e-1_dp * lambda * a * b * t50 * t195 + t34  &
          * (-0.2777777778e-1_dp * t9 * t180 * c * t195 - 0.2777777778e-1_dp * t9  &
          * t180 * t13 * t37 * my_ndrho * t25 * d + 0.27777777777777777778e-1_dp *  &
          t51 * t14 * my_ndrho * t63)

        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE lyp_adiabatic_lda_calc

! *****************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param lyp_adiabatic_params ...
!> \par History
!>      01.2008 created [fawzi]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE lyp_adiabatic_lsd_eval(rho_set,deriv_set,grad_deriv,lyp_adiabatic_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: lyp_adiabatic_params

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, lambda
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndr, &
      e_ndr_ndr, e_ndr_ra, e_ndr_rb, e_ndra, e_ndra_ndra, e_ndra_ra, &
      e_ndra_rb, e_ndrb, e_ndrb_ndrb, e_ndrb_ra, e_ndrb_rb, e_ra, e_ra_ra, &
      e_ra_rb, e_rb, e_rb_rb, norm_drho, norm_drhoa, norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(deriv, bo)

    CALL section_vals_val_get(lyp_adiabatic_params,"LAMBDA",r_val=lambda)
    CALL cite_reference(Lee1988)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    CALL xc_rho_set_get(rho_set,&
         rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
         norm_drhob=norm_drhob, norm_drho=norm_drho, &
         rho_cutoff=epsilon_rho,&
         local_bounds=bo)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
    ELSE
      dummy=> rhoa
    END IF
    e_0 => dummy
    e_ra => dummy
    e_rb => dummy
    e_ndra_ra => dummy
    e_ndra_rb => dummy
    e_ndrb_ra => dummy
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
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_0)
    END IF
    IF (grad_deriv==1.OR.grad_deriv==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ra)
      deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rb)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndr)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndra)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrb)
    END IF
    IF (grad_deriv>1.OR.grad_deriv<-1) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
           message="derivatives bigger than 1 not implemented", &
           error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp          shared(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob) &
    !$omp          shared(e_0, e_ra, e_rb, e_ndr, e_ndra, e_ndrb) &
    !$omp          shared(grad_deriv, npoints) &
    !$omp          shared(epsilon_rho, lambda)

    CALL lyp_adiabatic_lsd_calc(&
         rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, norm_drhoa=norm_drhoa,&
         norm_drhob=norm_drhob, e_0=e_0, e_ra=e_ra, e_rb=e_rb,&
         e_ndr=e_ndr,&
         e_ndra=e_ndra, e_ndrb=e_ndrb,&
         grad_deriv=grad_deriv, npoints=npoints, &
         epsilon_rho=epsilon_rho,lambda=lambda)

    !$omp end parallel

    IF (cp_debug) THEN
      DEALLOCATE(dummy)
    ELSE
      NULLIFY(dummy)
    END IF
    CALL timestop(handle)
  END SUBROUTINE lyp_adiabatic_lsd_eval

! *****************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param norm_drho ...
!> \param norm_drhoa ...
!> \param norm_drhob ...
!> \param e_0 ...
!> \param e_ra ...
!> \param e_rb ...
!> \param e_ndr ...
!> \param e_ndra ...
!> \param e_ndrb ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param lambda ...
!> \par History
!>      08.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE lyp_adiabatic_lsd_calc(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
     e_0, e_ra, e_rb,&
     e_ndr,&
     e_ndra, e_ndrb, &
     grad_deriv,npoints,epsilon_rho,lambda)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drho, &
                                                norm_drhoa, norm_drhob
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_ra, e_rb, e_ndr, &
                                                e_ndra, e_ndrb
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, lambda

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    REAL(KIND=dp) :: cf, my_ndrho, my_ndrhoa, my_ndrhob, my_rhoa, my_rhob, &
      t1, t10, t100, t102, t103, t106, t108, t113, t115, t118, t119, t124, &
      t125, t128, t129, t132, t135, t138, t14, t140, t141, t143, t145, t146, &
      t15, t151, t153, t157, t16, t162, t165, t169, t17, t171, t174, t179, &
      t18, t183, t186, t187, t188, t19, t194, t196, t199, t2, t200, t202, &
      t21, t212, t216, t220, t222, t223, t225, t23, t231, t237, t24, t246, &
      t25, t250, t259, t26, t266, t27, t270, t273, t276, t28, t280, t285, &
      t288, t294, t3, t30, t300, t307, t31, t316, t32, t325, t348, t351, &
      t355, t362, t387, t39, t394, t4, t41, t42
    REAL(KIND=dp) :: t421, t46, t47, t48, t49, t5, t51, t55, t58, t6, t62, &
      t63, t65, t67, t7, t73, t74, t76, t77, t78, t80, t83, t84, t85, t86, &
      t87, t9, t90, t91, t94, t95, t96, t97

    cf=0.3_dp*(3._dp*pi*pi)**(2._dp/3._dp)

    !$omp do

    DO ii=1,npoints
      my_rhoa=MAX(rhoa(ii),0.0_dp)
      my_rhob=MAX(rhob(ii),0.0_dp)
      IF (my_rhoa+my_rhob>epsilon_rho) THEN
        my_ndrhoa=norm_drhoa(ii)
        my_ndrhob=norm_drhob(ii)
        my_ndrho = norm_drho(ii)
        IF( grad_deriv >= 0 ) THEN
          t1 = a * my_rhoa
          t2 = my_rhoa + my_rhob
          t3 = 0.1e1_dp / t2
          t4 = my_rhob * t3
          t5 = d * lambda
          t6 = t2 ** (0.1e1_dp / 0.3e1_dp)
          t7 = 0.1e1_dp / t6
          t9 = 0.10e1_dp + t5 * t7
          t10 = 0.1e1_dp / t9
          t14 = a * b
          t15 = c * lambda
          t16 = t15 * t7
          t17 = EXP(-t16)
          t18 = t14 * t17
          t19 = t2 ** 2
          t21 = t6 ** 2
          t23 = 0.1e1_dp / t21 / t19 / t2
          t24 = t10 * t23
          t25 = my_rhoa * my_rhob
          t26 = my_rhoa ** 2
          t27 = my_rhoa ** (0.1e1_dp / 0.3e1_dp)
          t28 = t27 ** 2
          t30 = my_rhob ** 2
          t31 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
          t32 = t31 ** 2
          t39 = t5 * t7 * t10
          t41 = 0.26111111111111111111e1_dp - 0.3888888889e0_dp * t16 - 0.3888888889e0_dp  &
          * t39
          t42 = my_ndrho ** 2
          t46 = 0.25000000000000000000e1_dp - 0.5555555556e-1_dp * t16 - 0.5555555556e-1_dp  &
          * t39
          t47 = my_ndrhoa ** 2
          t48 = my_ndrhob ** 2
          t49 = t47 + t48
          t51 = t16 + t39 - 0.110e2_dp
          t55 = my_rhoa * t3 * t47 + t4 * t48
          t58 = 0.12699208415745595798e2_dp * Cf * (t28 * t26 + t32 * t30) + t41  &
          * t42 - t46 * t49 - 0.1111111111e0_dp * t51 * t55
          t62 = 0.66666666666666666667e0_dp * t19
          t63 = t62 - t26
          t65 = t62 - t30
          t67 = t25 * t58 - 0.6666666667e0_dp * t19 * t42 + t63 * t48 + t65 * t47
          t73 = lambda ** 2
          t74 = t1 * my_rhob
          t76 = 0.1e1_dp / t6 / t2
          t77 = t9 ** 2
          t78 = 0.1e1_dp / t77
          t80 = t76 * t78 * d
          t83 = t14 * c
          t84 = t19 ** 2
          t85 = 0.1e1_dp / t84
          t86 = t85 * t17
          t87 = t10 * t67
          t90 = t78 * t85
          t91 = t67 * d
          t94 = t17 * t10
          t95 = t14 * t94
          t96 = t23 * my_rhoa
          t97 = c * t7
          t100 = d * t7 * t10
          t102 = d ** 2
          t103 = t102 * lambda
          t106 = t103 / t21 * t78
          t108 = -0.3888888889e0_dp * t97 - 0.3888888889e0_dp * t100 + 0.38888888888888888889e0_dp  &
          * t106
          t113 = -0.5555555556e-1_dp * t97 - 0.5555555556e-1_dp * t100 + 0.55555555555555555556e-1_dp  &
          * t106
          t115 = t97 + t100 - t106
          t118 = t108 * t42 - t113 * t49 - 0.1111111111e0_dp * t115 * t55
          t119 = my_rhob * t118

          e_0(ii) = e_0(ii) +  0.20e1_dp * lambda * (-0.40e1_dp * t1 * t4 * t10 - t18 * t24 * t67)  &
          + t73 * (0.40e1_dp * t74 * t80 + t83 * t86 * t87 + t18 * t90 * t91 -  &
          t95 * t96 * t119)

        END IF
        IF( grad_deriv == 1 .OR. grad_deriv == -1 ) THEN
          t124 = a * my_rhob
          t125 = t3 * t10
          t128 = 0.1e1_dp / t19
          t129 = my_rhob * t128
          t132 = 0.40e1_dp * t1 * t129 * t10
          t135 = 0.1e1_dp / t6 / t19 * t78
          t138 = 0.1333333333e1_dp * t74 * t135 * t5
          t140 = t84 * t2
          t141 = 0.1e1_dp / t140
          t143 = t141 * t17 * t87
          t145 = t14 * t15 * t143 / 0.3e1_dp
          t146 = t17 * t78
          t151 = t14 * t146 * t141 * t67 * t5 / 0.3e1_dp
          t153 = 0.1e1_dp / t21 / t84
          t157 = 0.11e2_dp / 0.3e1_dp * t18 * t10 * t153 * t67
          t162 = t15 * t76
          t165 = t5 * t76 * t10
          t169 = 0.1e1_dp / t21 / t2
          t171 = t102 * t73 * t169 * t78
          t174 = (0.12962962962962962963e0_dp * t162 + 0.12962962962962962963e0_dp  &
          * t165 - 0.1296296296e0_dp * t171) * t42
          t179 = (0.18518518518518518519e-1_dp * t162 + 0.18518518518518518519e-1_dp  &
          * t165 - 0.1851851852e-1_dp * t171) * t49
          t183 = 0.1111111111e0_dp * (-t162 / 0.3e1_dp - t165 / 0.3e1_dp + t171 / 0.3e1_dp)  &
          * t55
          t186 = my_rhoa * t128 * t47
          t187 = t129 * t48
          t188 = t3 * t47 - t186 - t187
          t194 = 0.1333333333e1_dp * t2 * t42
          t196 = 0.13333333333333333333e1_dp * my_rhob
          t199 = 0.13333333333333333333e1_dp * my_rhoa
          t200 = t199 + t196
          t202 = my_rhob * t58 + t25 * (0.33864555775321588795e2_dp * Cf * t28 * my_rhoa  &
          + t174 - t179 - t183 - 0.1111111111e0_dp * t51 * t188) - t194 + (-0.6666666667e0_dp  &
          * my_rhoa + t196) * t48 + t200 * t47
          t212 = 0.5333333333e1_dp * t74 * t135 * d
          t216 = 0.1e1_dp / t77 / t9
          t220 = 0.26666666666666666667e1_dp * t74 / t21 / t19 * t216 * t103
          t222 = 4 * t83 * t143
          t223 = c ** 2
          t225 = 0.1e1_dp / t6 / t140
          t231 = t14 * t223 * t225 * lambda * t17 * t87 / 0.3e1_dp
          t237 = 0.2e1_dp / 0.3e1_dp * t14 * c * t225 * t146 * t91 * lambda
          t246 = 0.2e1_dp / 0.3e1_dp * t14 * t17 * t216 * t225 * t67 * t103
          t250 = 4 * t18 * t78 * t141 * t91
          t259 = t14 * t15 * t141 * t94 * t25 * t118 / 0.3e1_dp
          t266 = t14 * t146 * t141 * t25 * t118 * d * lambda / 0.3e1_dp
          t270 = 0.11e2_dp / 0.3e1_dp * t95 * t153 * my_rhoa * t119
          t273 = c * t76
          t276 = d * t76 * t10
          t280 = t102 * t169 * t78 * lambda
          t285 = t102 * d * t73 * t128 * t216
          t288 = (0.12962962962962962963e0_dp * t273 + 0.12962962962962962963e0_dp  &
          * t276 - 0.3888888889e0_dp * t280 + 0.25925925925925925926e0_dp * t285)  &
          * t42
          t294 = (0.18518518518518518519e-1_dp * t273 + 0.18518518518518518519e-1_dp  &
          * t276 - 0.5555555556e-1_dp * t280 + 0.37037037037037037037e-1_dp * t285)  &
          * t49
          t300 = 0.1111111111e0_dp * (-t273 / 0.3e1_dp - t276 / 0.3e1_dp + t280 - 0.2e1_dp  &
          / 0.3e1_dp * t285) * t55
          t307 = 0.40e1_dp * t124 * t80 - t212 + t220 - t222 + t231 + t237 + t83  &
          * t86 * t10 * t202 + t246 - t250 + t18 * t90 * t202 * d - t259 -  &
          t266 + t270 - t18 * t24 * t119 - t95 * t96 * my_rhob * (t288 - t294 -  &
          t300 - 0.1111111111e0_dp * t115 * t188)

          e_ra(ii) = e_ra(ii) +  0.20e1_dp * lambda * (-0.40e1_dp * t124 * t125 + t132 - t138 - t145  &
          - t151 + t157 - t18 * t24 * t202) + t73 * t307

          t316 = -t186 + t3 * t48 - t187
          t325 = my_rhoa * t58 + t25 * (0.33864555775321588795e2_dp * Cf * t32 * my_rhob  &
          + t174 - t179 - t183 - 0.1111111111e0_dp * t51 * t316) - t194 + t200  &
          * t48 + (t199 - 0.6666666667e0_dp * my_rhob) * t47
          t348 = 0.40e1_dp * t1 * t80 - t212 + t220 - t222 + t231 + t237 + t83 *  &
          t86 * t10 * t325 + t246 - t250 + t18 * t90 * t325 * d - t259 - t266  &
          + t270 - t18 * t24 * my_rhoa * t118 - t95 * t96 * my_rhob * (t288 - t294  &
          - t300 - 0.1111111111e0_dp * t115 * t316)

          e_rb(ii) = e_rb(ii) +  0.20e1_dp * lambda * (-0.40e1_dp * t1 * t125 + t132 - t138 - t145 -  &
          t151 + t157 - t18 * t24 * t325) + t73 * t348

          t351 = lambda * a * b
          t355 = t3 * my_ndrhoa
          t362 = t25 * (-REAL(2 * t46 * my_ndrhoa,dp) - 0.2222222222e0_dp * t51 * my_rhoa  &
          * t355) + REAL(2 * t65 * my_ndrhoa,dp)

          e_ndra(ii) = e_ndra(ii)   -0.20e1_dp * t351 * t94 * t23 * t362 + t73 * (t83 * t86 * t10 *  &
          t362 + t18 * t90 * t362 * d - t95 * t96 * my_rhob * (-REAL(2 * t113 *  &
          my_ndrhoa,dp) - 0.2222222222e0_dp * t115 * my_rhoa * t355))

          t387 = t3 * my_ndrhob
          t394 = t25 * (-REAL(2 * t46 * my_ndrhob,dp) - 0.2222222222e0_dp * t51 * my_rhob  &
          * t387) + REAL(2 * t63 * my_ndrhob,dp)

          e_ndrb(ii) = e_ndrb(ii)   -0.20e1_dp * t351 * t94 * t23 * t394 + t73 * (t83 * t86 * t10 *  &
          t394 + t18 * t90 * t394 * d - t95 * t96 * my_rhob * (-REAL(2 * t113 *  &
          my_ndrhob,dp) - 0.2222222222e0_dp * t115 * my_rhob * t387))

          t421 = REAL(2 * t25 * t41 * my_ndrho,dp) - 0.1333333333e1_dp * REAL(t19,dp) * REAL(my_ndrho,dp)

          e_ndr(ii) = e_ndr(ii)   -0.20e1_dp * t351 * t94 * t23 * t421 + t73 * (t83 * t86 * t10 * t421  &
          + t18 * t90 * t421 * d - REAL(2 * t95 * t96 * my_rhob * t108 * my_ndrho,dp))

        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE lyp_adiabatic_lsd_calc

END MODULE xc_lyp_adiabatic
