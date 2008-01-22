!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/xc_lyp_adiabatic [1.0] *
!!
!!   NAME
!!     xc_lyp
!!
!!   FUNCTION
!!     calculates the lyp correlation functional
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_lyp_adiabatic
  USE bibliography,                    ONLY: Lee1988,&
                                             cite_reference
  USE f77_blas
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
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_lyp'
  REAL(kind=dp), PARAMETER, PRIVATE :: a=0.04918_dp, b=0.132_dp,&
       c=0.2533_dp,d=0.349_dp

  PUBLIC :: lyp_adiabatic_lda_info, lyp_adiabatic_lsd_info, lyp_adiabatic_lda_eval, lyp_adiabatic_lsd_eval
!!***
!****************************************************************************

CONTAINS

!!****f* xc_lyp_adiabatic/lyp_adiabatic_lda_info [1.0] *
!!
!!   NAME
!!     lyp_adiabatic_lda_info
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [mguidon]
!!
!!*** **********************************************************************
  SUBROUTINE lyp_adiabatic_lda_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

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
!***************************************************************************

!!****f* xc_lyp_adiabatic/lyp_adiabatic_lsd_info [1.0] *
!!
!!   NAME
!!     lyp_adiabatic_lsd_info
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [mguidon]
!!
!!*** **********************************************************************
  SUBROUTINE lyp_adiabatic_lsd_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

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
!***************************************************************************

!!****f* xc_lyp_adiabatic/lyp_adiabatic_lda_eval [1.0] *
!!
!!   NAME
!!     lyp_adiabatic_lda_eval
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [mguidon]
!!
!!*************************************************************************
  SUBROUTINE lyp_adiabatic_lda_eval(rho_set,deriv_set,grad_deriv,lyp_adiabatic_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: lyp_adiabatic_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, lambda
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, e_ndrho_rho_rho, e_rho, &
      e_rho_rho, e_rho_rho_rho, norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(lyp_adiabatic_params,"LAMBDA",r_val=lambda,error=error)
    CALL cite_reference(Lee1988)

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
      IF (grad_deriv>1.OR.grad_deriv<-1) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 1 not implemented", &
             error=error, error_level=cp_failure_level)
      END IF

      CALL lyp_adiabatic_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
           e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
           e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
           e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
           e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,&
           grad_deriv=grad_deriv,&
           npoints=npoints,epsilon_rho=epsilon_rho,epsilon_norm_drho=epsilon_norm_drho,lambda=lambda,&
           error=error)
      IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ELSE
        NULLIFY(dummy)
      END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE lyp_adiabatic_lda_eval
!***************************************************************************

!!****f* xc_lyp_adiabatic/lyp_adiabatic_lda_calc [1.0] *
!!
!!   NAME
!!     lyp_adiabatic_lda_eval
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [mguidon]
!!
!!*** **********************************************************************
  SUBROUTINE lyp_adiabatic_lda_calc(rho, rho_1_3, norm_drho,&
     e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
     e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
     grad_deriv,npoints,epsilon_rho, epsilon_norm_drho,lambda,error)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: e_ndrho_ndrho_rho, &
      e_ndrho_rho_rho, e_rho_rho_rho, e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, &
      e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho_1_3, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, lambda
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lda_calc', &
      routineP = moduleN//':'//routineN

  END SUBROUTINE lyp_adiabatic_lda_calc
!***************************************************************************

!!****f* xc_lyp_adiabatic/lyp_adiabatic_lsd_eval [1.0] *
!!
!!   NAME
!!     lyp_adiabatic_lsd_eval
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE lyp_adiabatic_lsd_eval(rho_set,deriv_set,grad_deriv,lyp_adiabatic_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: lyp_adiabatic_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho, lambda
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndr, &
      e_ndr_ndr, e_ndr_ra, e_ndr_rb, e_ndra, e_ndra_ndra, e_ndra_ra, &
      e_ndra_rb, e_ndrb, e_ndrb_ndrb, e_ndrb_ra, e_ndrb_rb, e_ra, e_ra_ra, &
      e_ra_rb, e_rb, e_rb_rb, norm_drho, norm_drhoa, norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(deriv, bo)

    CALL section_vals_val_get(lyp_adiabatic_params,"LAMBDA",r_val=lambda,error=error)
    CALL cite_reference(Lee1988)

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
              allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)
      END IF
      IF (grad_deriv==1.OR.grad_deriv==-1) THEN
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
      IF (grad_deriv>1.OR.grad_deriv<-1) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 1 not implemented", &
             error=error, error_level=cp_failure_level)
      END IF

      CALL lyp_adiabatic_lsd_calc(&
           rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, norm_drhoa=norm_drhoa,&
           norm_drhob=norm_drhob, e_0=e_0, e_ra=e_ra, e_rb=e_rb,&
           & e_ndra_ra=e_ndra_ra, e_ndra_rb=e_ndra_rb, e_ndrb_ra&
           &=e_ndrb_ra, e_ndrb_rb=e_ndrb_rb, e_ndr_ndr=e_ndr_ndr,&
           e_ndra_ndra=e_ndra_ndra, e_ndrb_ndrb=e_ndrb_ndrb, e_ndr=e_ndr,&
           e_ndra=e_ndra, e_ndrb=e_ndrb, e_ra_ra=e_ra_ra, &
           e_ra_rb=e_ra_rb, e_rb_rb=e_rb_rb, e_ndr_ra=e_ndr_ra,&
           e_ndr_rb=e_ndr_rb,&
           grad_deriv=grad_deriv, npoints=npoints, &
           epsilon_rho=epsilon_rho,epsilon_drho=epsilon_drho,lambda=lambda,error=error)
 
      IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ELSE
        NULLIFY(dummy)
      END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE lyp_adiabatic_lsd_eval
!***************************************************************************

!!****f* xc_lyp_adiabatic/lyp_adiabatic_lsd_calc [1.0] *
!!
!!   NAME
!!     lyp_adiabatic_lsd_calc
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     08.2008 created [mguidon]
!!
!!*** **********************************************************************
  SUBROUTINE lyp_adiabatic_lsd_calc(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
     e_0, e_ra, e_rb, e_ndra_ra, e_ndra_rb, e_ndrb_ra, e_ndrb_rb, e_ndr_ndr,&
     e_ndra_ndra, e_ndrb_ndrb, e_ndr,&
     e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ndr_ra, e_ndr_rb,&
     grad_deriv,npoints,epsilon_rho,epsilon_drho,lambda,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drho, &
                                                norm_drhoa, norm_drhob
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_ra, e_rb, e_ndra_ra, &
      e_ndra_rb, e_ndrb_ra, e_ndrb_rb, e_ndr_ndr, e_ndra_ndra, e_ndrb_ndrb, &
      e_ndr, e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ndr_ra, e_ndr_rb
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho, lambda
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_adiabatic_lsd_calc', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER :: small = 1.0e-20, t_1_9 = 1.0_dp/9.0_dp, &
      t_2_3 = 2.0_dp/3.0_dp, t_4_27 = 4.0_dp/27.0_dp, t_4_3 = 4.0_dp/3.0_dp, &
      t_4_9 = 4.0_dp/9.0_dp, t_8_27 = 8.0_dp/27.0_dp

    REAL(KIND=dp) ::  cf, my_rhoa, my_rhob, my_ndrhoa, my_ndrhob, my_ndrho


    REAL(KIND=dp) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t13,t14,t15,t16,t17,t19,t22,&
        t23,t24,t25,t26,t28,t29,t30,t36,t38,t39,t43,t44,t45,t46,t48,t52,t55,&
        t59,t60,t62,t64,t68,t71,t72,t75,t79,t80,t84,t86,t88,t93,t98,t104,t110,&
        t111,t114,t116,t120,t123,t128,t132,t135,t136,t143,t145,t148,t149


    INTEGER                                  :: ii
    LOGICAL                                  :: failure

    cf=0.3_dp*(3._dp*pi*pi)**(2._dp/3._dp)

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
          t5 = t2 ** (0.1e1_dp / 0.3e1_dp)
          t6 = 0.1e1_dp / t5
          t7 = d * t6
          t8 = 0.10e1_dp + t7
          t9 = 0.1e1_dp / t8
          t13 = a * b
          t14 = c * t6
          t15 = exp(-t14)
          t16 = t13 * t15
          t17 = t2 ** 2
          t19 = t5 ** 2
          t22 = t9 / t19 / t17 / t2
          t23 = my_rhoa * my_rhob
          t24 = my_rhoa ** 2
          t25 = my_rhoa ** (0.1e1_dp / 0.3e1_dp)
          t26 = t25 ** 2
          t28 = my_rhob ** 2
          t29 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
          t30 = t29 ** 2
          t36 = t7 * t9
          t38 = 0.26111111111111111111e1_dp - 0.3888888889e0_dp * t14 - 0.3888888889e0_dp  &
          * t36
          t39 = my_ndrho ** 2
          t43 = 0.25000000000000000000e1_dp - 0.5555555556e-1_dp * t14 - 0.5555555556e-1_dp  &
          * t36
          t44 = my_ndrhoa ** 2
          t45 = my_ndrhob ** 2
          t46 = t44 + t45
          t48 = t14 + t36 - 0.110e2_dp
          t52 = my_rhoa * t3 * t44 + t4 * t45
          t55 = 0.12699208415745595798e2_dp * Cf * (t26 * t24 + t30 * t28) + t38  &
          * t39 - t43 * t46 - 0.1111111111e0_dp * t48 * t52
          t59 = 0.66666666666666666667e0_dp * t17
          t60 = t59 - t24
          t62 = t59 - t28
          t64 = t23 * t55 - 0.6666666667e0_dp * t17 * t39 + t60 * t45 + t62 * t44

          e_0(ii) = e_0(ii)  -0.40e1_dp * t1 * t4 * t9 - t16 * t22 * t64
        END IF
        IF(grad_deriv>=1 .OR. grad_deriv ==-1) THEN
          t68 = t3 * t9
          t71 = 0.1e1_dp / t17
          t72 = my_rhob * t71
          t75 = 0.40e1_dp * t1 * t72 * t9
          t79 = t8 ** 2
          t80 = 0.1e1_dp / t79
          t84 = 0.1333333333e1_dp * t1 * my_rhob / t5 / t17 * t80 * d
          t86 = t17 ** 2
          t88 = 0.1e1_dp / t86 / t2
          t93 = t13 * c * t88 * t15 * t9 * t64 / 0.3e1_dp
          t98 = t16 * t80 * t88 * t64 * d / 0.3e1_dp
          t104 = 0.11e2_dp / 0.3e1_dp * t16 * t9 / t19 / t86 * t64
          t110 = 0.1e1_dp / t5 / t2
          t111 = c * t110
          t114 = d * t110 * t9
          t116 = d ** 2
          t120 = t116 / t19 / t2 * t80
          t123 = (0.12962962962962962963e0_dp * t111 + 0.12962962962962962963e0_dp  &
          * t114 - 0.1296296296e0_dp * t120) * t39
          t128 = (0.18518518518518518519e-1_dp * t111 + 0.18518518518518518519e-1_dp  &
          * t114 - 0.1851851852e-1_dp * t120) * t46
          t132 = 0.1111111111e0_dp * (-t111 / 0.3e1_dp - t114 / 0.3e1_dp + t120 / 0.3e1_dp)  &
          * t52
          t135 = my_rhoa * t71 * t44
          t136 = t72 * t45
          t143 = 0.1333333333e1_dp * t2 * t39
          t145 = 0.13333333333333333333e1_dp * my_rhob
          t148 = 0.13333333333333333333e1_dp * my_rhoa
          t149 = t148 + t145
          e_ra(ii) = e_ra(ii)   -0.40e1_dp * a * my_rhob * t68 + t75 - t84 - t93 - t98 + t104 - t16  &
          * t22 * (my_rhob * t55 + t23 * (0.33864555775321588795e2_dp * Cf * t26 *  &
          my_rhoa + t123 - t128 - t132 - 0.1111111111e0_dp * t48 * (t3 * t44 - t135  &
          - t136)) - t143 + (-0.6666666667e0_dp * my_rhoa + t145) * t45 + t149  &
          * t44)
          e_rb(ii) = e_rb(ii)   -0.40e1_dp * t1 * t68 + t75 - t84 - t93 - t98 + t104 - t16 * t22  &
          * (my_rhoa * t55 + t23 * (0.33864555775321588795e2_dp * Cf * t30 * my_rhob  &
          + t123 - t128 - t132 - 0.1111111111e0_dp * t48 * (-t135 + t3 * t45 -  &
          t136)) - t143 + t149 * t45 + (t148 - 0.6666666667e0_dp * my_rhob) * t44)
          e_ndra(ii) = e_ndra(ii)   -t16 * t22 * (t23 * (-dble(2 * t43 * my_ndrhoa) - 0.2222222222e0_dp  &
          * dble(t48) * dble(my_rhoa) * dble(t3) * dble(my_ndrhoa)) + dble(2 * t62  &
          * my_ndrhoa))
          e_ndrb(ii) = e_ndrb(ii)   -t16 * t22 * (t23 * (-dble(2 * t43 * my_ndrhob) - 0.2222222222e0_dp  &
          * dble(t48) * dble(my_rhob) * dble(t3) * dble(my_ndrhob)) + dble(2 * t60  &
          * my_ndrhob))
          e_ndr(ii) = e_ndr(ii)   -t16 * t22 * (dble(2 * t23 * t38 * my_ndrho) - 0.1333333333e1_dp * dble(t17)  &
          * dble(my_ndrho))
        END IF
      END IF
    END DO
  END SUBROUTINE lyp_adiabatic_lsd_calc
!***************************************************************************

END MODULE xc_lyp_adiabatic
