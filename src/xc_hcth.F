!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_hcth [1.0] *
!!
!!   NAME
!!     xc_hcth
!!
!!   FUNCTION
!!     calculate the Hamprecht, Cohen, Tozer, and Handy (HCTH) exchange 
!!     functional
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!****************************************************************************
MODULE xc_hcth
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
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
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_hcth'

  PUBLIC :: hcth_lda_info, hcth_lda_eval
!!***
  !****************************************************************************
CONTAINS

!!****f* xc_hcth/hcth93_lda_info [1.0] *
!!
!!   NAME
!!     hcth93_lda_info
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
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE hcth_lda_info(iparset,reference,shortform, needs, max_deriv, error)
    INTEGER, INTENT(in)                      :: iparset
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hcth93_lda_info', &
      routineP = moduleN//':'//routineN

    SELECT CASE(iparset)
    CASE(93)
       IF ( PRESENT ( reference ) ) THEN
          reference = "F. A. Hamprecht, A. J. Cohen, D. J. Tozer, and N. C. Handy, J. Chem. Phys. 109, 6264 (1998);"//&
               " HCTH/93 xc functional {LDA version}"
       END IF
       IF ( PRESENT ( shortform ) ) THEN
          shortform = "HCTH/93 xc energy functional (LDA)"
       END IF
    CASE(-120)
       IF ( PRESENT ( reference ) ) THEN
          reference = "A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik, J. Chem. Phys. 112, 1670 (2000);"//&
               " HCTH/120 xc functional {LDA CPMD version}"
       END IF
       IF ( PRESENT ( shortform ) ) THEN
          shortform = "HCTH/120 xc energy functional (LDA CPMD)"
       END IF
    CASE(120)
       IF ( PRESENT ( reference ) ) THEN
          reference = "A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik, J. Chem. Phys. 112, 1670 (2000);"//&
               " HCTH/120 xc functional {LDA version}"
       END IF
       IF ( PRESENT ( shortform ) ) THEN
          shortform = "HCTH/120 xc energy functional (LDA)"
       END IF
    CASE(147)
       IF ( PRESENT ( reference ) ) THEN
          reference = "A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik, J. Chem. Phys. 112, 1670 (2000);"//&
               " HCTH/147 xc functional {LDA Version}"
       END IF
       IF ( PRESENT ( shortform ) ) THEN
          shortform = "HCTH/147 xc energy functional (LDA)"
       END IF
    CASE(407)
       IF ( PRESENT ( reference ) ) THEN
          reference = "A. D. Boese and N. C. Handy, J. Chem. Phys. 114, 5497 (2001);"//&
               "HCTH/407 xc functional {LDA version}"
       END IF
       IF ( PRESENT ( shortform ) ) THEN
          shortform = "HCTH/407 xc energy functional (LDA)"
       END IF
    CASE(-407)
       IF ( PRESENT ( reference ) ) THEN
          reference = "A. D. Boese and N. C. Handy, J. Chem. Phys. 114, 5497 (2001);"//&
               "HCTH/407 xc functional {LDA version}"
       END IF
       IF ( PRESENT ( shortform ) ) THEN
          shortform = "HCTH/407 xc energy functional (LDA CPMD)"
       END IF
    CASE default
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Invalid HCTH parameter set requested ("//cp_to_string(iparset)//")",&
            error)
    END SELECT
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE hcth_lda_info
  !***************************************************************************

!!****f* xc_hcth/hcth_lda_eval [1.0] *
!!
!!   NAME
!!     hcth_lda_eval
!!
!!   SYNOPSIS
!!     Subroutine hcth_lda_eval(iparset, rho_set, deriv_set, grad_deriv,&
!!         error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: grad_deriv, iparset
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine hcth_lda_eval
!!
!!   FUNCTION
!!     evaluates the hcth functional for lda
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
!!     - iparset: the parameter set that should be used (93,120,147,407)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE hcth_lda_eval(iparset,rho_set,deriv_set,grad_deriv,error)
    INTEGER, INTENT(in)                      :: iparset
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hcth_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_ndrho, e_rho, &
                                                norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    failure=.FALSE.
    NULLIFY(bo,e_0, e_ndrho, e_rho, norm_drho, rho)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
            drho_cutoff=epsilon_drho,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)


       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       IF (iparset==-120) THEN
          CALL hcth120_lda_calc(rho=rho, norm_drho=norm_drho,&
               e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
               npoints=npoints,epsilon_rho=epsilon_rho,&
               epsilon_drho=epsilon_drho, error=error)
       ELSEIF (iparset==-407) THEN
          CALL hcth407_lda_calc(rho=rho, norm_drho=norm_drho,&
               e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
               npoints=npoints,epsilon_rho=epsilon_rho,&
               epsilon_drho=epsilon_drho, error=error)
       ELSE
          CALL hcth_lda_calc(iparset=iparset,rho=rho, norm_drho=norm_drho,&
               e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
               npoints=npoints,epsilon_rho=epsilon_rho,&
               epsilon_drho=epsilon_drho, error=error)
       END IF
    END IF
  END SUBROUTINE hcth_lda_eval
  !***************************************************************************

!!****f* xc_hcth/hcth_lda_calc *
!!
!!   NAME
!!     hcth_lda_calc
!!
!!   FUNCTION
!!     Calculate the gradient-corrected xc energy and potential
!!     of Hamprecht, Cohen, Tozer, and Handy (HCTH) for a closed shell
!!     density.
!!
!!   NOTES
!!    Literature:- F. A. Hamprecht, A. J. Cohen, D. J. Tozer, and N. C. Handy,
!!                 J. Chem. Phys. 109, 6264 (1998) -> HCTH/93
!!               - A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik,
!!                 J. Chem. Phys. 112, 1670 (2000) -> HCTH/120 and HCTH/147
!!               - A. D. Boese and N. C. Handy,
!!                 J. Chem. Phys. 114, 5497 (2001) -> HCTH/407
!!               - J. P. Perdew and Y. Wang,
!!                 Phys. Rev. B 45, 13244 (1992) -> PW92
!!     (routine shamelessly copied from the one that Matthias Krack had 
!!     taken from CPMD)
!!
!!   ARGUMENTS
!!     - iparset: the parameter set that should be used (93,120,147,407)
!!     - rho: the density
!!     - norm_drho: the norm of the gradient of the density
!!     - e_0: the value of the functional in that point
!!     - e_rho: the derivative of the functional wrt. rho
!!     - e_ndrho: the derivative of the functional wrt. norm_drho
!!     - epsilon_rho: the cutoff on rho
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 rewamped [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE hcth_lda_calc(iparset,rho,norm_drho,e_0,e_rho,e_ndrho,&
       epsilon_rho,epsilon_drho,npoints,error)
    INTEGER, INTENT(IN)                      :: iparset
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0, e_rho, e_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(IN)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "hcth_lda_calc", &
      routineP = moduleN//":"//routineN
    REAL(KIND=dp), DIMENSION(4), PARAMETER :: &
      beta0 = (/ 7.59570_dp,3.58760_dp,1.63820_dp,0.49294_dp/), &
      beta1 = (/14.11890_dp,6.19770_dp,3.36620_dp,0.62517_dp/)
    REAL(KIND=dp), PARAMETER :: a0 = 0.031091_dp, a1 = 0.015545_dp, &
      alpha0 = 0.21370_dp, alpha1 = 0.20548_dp, f13 = 1.0_dp/3.0_dp, &
      f43 = 4.0_dp*f13, f83 = 8.0_dp*f13, gamma_cab = 0.006_dp, &
      gamma_css = 0.200_dp, gamma_xss = 0.004_dp

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: cx_vwn_e, cx_vwn_v, dgcabddrho, dgcabdrho, dgcabds, &
      dgcssddrho, dgcssdrho, dgcssds, dgdrs, dgxssddrho, dgxssdrho, dgxssds, &
      drho, drhos, drsdrho, ecab, ecss, exss, g, gcab, gcss, gs2, gxss, &
      my_rho, p, q, rho13, rho43, rhos, rhos13, rhos43, rs, rs12, rsfac, s, &
      s2, two13, u, vcab, vcss, vxss, x, y
    REAL(KIND=dp), DIMENSION(0:4)            :: ccab, ccss, cxss

    failure=.FALSE.
    cx_vwn_e = -0.75_dp*(3.0_dp/pi)**f13
    cx_vwn_v = f43*cx_vwn_e
    rsfac = (f43*pi)**(-f13)
    two13 = 2.0_dp**f13

    !   *** LSDA correlation parametrisation (PW92) ***
    !   *** GGA parametrisation (HCTH/iparset) ***
    !   ---------------------------------------------------------------------------

    !     *** Load the HCTH parameter set HCTH/iparset ***

    SELECT CASE (iparset)
    CASE (93)
       cxss(0) =  0.109320E+01_dp
       ccss(0) =  0.222601E+00_dp
       ccab(0) =  0.729974E+00_dp
       cxss(1) = -0.744056E+00_dp
       ccss(1) = -0.338622E-01_dp
       ccab(1) =  0.335287E+01_dp
       cxss(2) =  0.559920E+01_dp
       ccss(2) = -0.125170E-01_dp
       ccab(2) = -0.115430E+02_dp
       cxss(3) = -0.678549E+01_dp
       ccss(3) = -0.802496E+00_dp
       ccab(3) =  0.808564E+01_dp
       cxss(4) =  0.449357E+01_dp
       ccss(4) =  0.155396E+01_dp
       ccab(4) = -0.447857E+01_dp
    CASE (120)
       cxss(0) =  0.109163E+01_dp
       ccss(0) =  0.489508E+00_dp
       ccab(0) =  0.514730E+00_dp
       cxss(1) = -0.747215E+00_dp
       ccss(1) = -0.260699E+00_dp
       ccab(1) =  0.692982E+01_dp
       cxss(2) =  0.507833E+01_dp
       ccss(2) =  0.432917E+00_dp
       ccab(2) = -0.247073E+02_dp
       cxss(3) = -0.410746E+01_dp
       ccss(3) = -0.199247E+01_dp
       ccab(3) =  0.231098E+02_dp
       cxss(4) =  0.117173E+01_dp
       ccss(4) =  0.248531E+01_dp
       ccab(4) = -0.113234E+02_dp
    CASE (147)
       cxss(0) =  0.109025E+01_dp
       ccss(0) =  0.562576E+00_dp
       ccab(0) =  0.542352E+00_dp
       cxss(1) = -0.799194E+00_dp
       ccss(1) =  0.171436E-01_dp
       ccab(1) =  0.701464E+01_dp
       cxss(2) =  0.557212E+01_dp
       ccss(2) = -0.130636E+01_dp
       ccab(2) = -0.283822E+02_dp
       cxss(3) = -0.586760E+01_dp
       ccss(3) =  0.105747E+01_dp
       ccab(3) =  0.350329E+02_dp
       cxss(4) =  0.304544E+01_dp
       ccss(4) =  0.885429E+00_dp
       ccab(4) = -0.204284E+02_dp
    CASE (407)
       cxss(0) =  0.108184E+01_dp
       ccss(0) =  0.118777E+01_dp
       ccab(0) =  0.589076E+00_dp
       cxss(1) = -0.518339E+00_dp
       ccss(1) = -0.240292E+01_dp
       ccab(1) =  0.442374E+01_dp
       cxss(2) =  0.342562E+01_dp
       ccss(2) =  0.561741E+01_dp
       ccab(2) = -0.192218E+02_dp
       cxss(3) = -0.262901E+01_dp
       ccss(3) = -0.917923E+01_dp
       ccab(3) =  0.425721E+02_dp
       cxss(4) =  0.228855E+01_dp
       ccss(4) =  0.624798E+01_dp
       ccab(4) = -0.420052E+02_dp
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Invalid HCTH parameter set requested ("//cp_to_string(iparset)//")",&
            error)
    END SELECT

    !$omp parallel do default(none) shared(rho,norm_drho,cxss,ccss,&
    !$omp         ccab,cx_vwn_e, cx_vwn_v, rsfac, two13,epsilon_rho) &
    !$omp       private(ii,dgcabddrho, dgcabdrho, dgcabds, dgcssddrho, &
    !$omp         dgcssdrho, dgcssds, dgdrs, dgxssddrho, dgxssdrho, dgxssds,&
    !$omp         drhos, drsdrho, ecab, ecss, exss, g, gcab, gcss, gs2, &
    !$omp         gxss, p, q, rho13, rho43, rhos, rhos13, rhos43, rs, rs12,&
    !$omp         s, s2, u, vcab, vcss, vxss, x, y, my_rho, drho)
    DO ii=1,npoints
       !     *** rho_sigma = rho/2 = rho_alpha = rho_beta (same for |nabla rho|) ***

       IF (rho(ii)>epsilon_rho) THEN
          my_rho=MAX(rho(ii),epsilon_rho)
          drho=norm_drho(ii)
          rhos = 0.5_dp*my_rho
          drhos = 0.5_dp*drho

          rhos13 = rhos**f13
          rhos43 = rhos13*rhos

          rho13 = two13*rhos13
          rho43 = rho13*my_rho

          !     *** LSDA exchange part (VWN) ***

          exss = cx_vwn_e*rho43
          vxss = cx_vwn_v*rho13

          !     *** LSDA correlation part (PW92) ***

          !     *** G(rho_sigma,0) => spin polarisation zeta = 1 ***

          rs = rsfac/rhos13
          rs12 = SQRT(rs)
          q = 2.0_dp*a1*(beta1(1) + (beta1(2) + (beta1(3) +&
               beta1(4)*rs12)*rs12)*rs12)*rs12
          p = 1.0_dp + 1.0_dp/q
          x = -2.0_dp*a1*(1.0_dp + alpha1*rs)
          y = LOG(p)
          g = x*y
          dgdrs = -2.0_dp*a1*alpha1*y -&
               x*a1*(beta1(1)/rs12 + 2.0_dp*beta1(2) +&
               3.0_dp*beta1(3)*rs12 + 4.0_dp*beta1(4)*rs)/(p*q*q)
          drsdrho = -f13*rs/my_rho
          ecss = my_rho*g
          vcss = g + my_rho*dgdrs*drsdrho

          !     *** G(rho_alpha,rho_beta) => spin polarisation zeta = 0 ***

          rs = rsfac/rho13
          rs12 = SQRT(rs)
          q = 2.0_dp*a0*(beta0(1) + (beta0(2) + (beta0(3) +&
               beta0(4)*rs12)*rs12)*rs12)*rs12
          p = 1.0_dp + 1.0_dp/q
          x = -2.0_dp*a0*(1.0_dp + alpha0*rs)
          y = LOG(p)
          g = x*y
          dgdrs = -2.0_dp*a0*alpha0*y -&
               x*a0*(beta0(1)/rs12 + 2.0_dp*beta0(2) +&
               3.0_dp*beta0(3)*rs12 + 4.0_dp*beta0(4)*rs)/(p*q*q)
          drsdrho = -f13*rs/my_rho
          ecab = my_rho*g - ecss
          vcab = g + my_rho*dgdrs*drsdrho - vcss

          !     *** GGA part (HCTH) ***

          s = drhos/rhos43
          s2 = s*s
          x = -f83/my_rho
          y = 2.0_dp/(drho*drho)

          !     *** g_x(rho_sigma,rho_sigma) ***

          gs2 = gamma_xss*s2
          q = 1.0_dp/(1.0_dp + gs2)
          u = gs2*q
          gxss = cxss(0) + (cxss(1) + (cxss(2) + (cxss(3) + cxss(4)*u)*u)*u)*u
          dgxssds = q*(cxss(1) + (2.0_dp*cxss(2) + (3.0_dp*cxss(3) +&
               4.0_dp*cxss(4)*u)*u)*u)*u
          dgxssdrho = x*dgxssds
          dgxssddrho = y*dgxssds

          !     *** g_c(rho_sigma,rho_sigma) ***

          gs2 = gamma_css*s2
          q = 1.0_dp/(1.0_dp + gs2)
          u = gs2*q
          gcss = ccss(0) + (ccss(1) + (ccss(2) + (ccss(3) + ccss(4)*u)*u)*u)*u
          dgcssds = q*(ccss(1) + (2.0_dp*ccss(2) + (3.0_dp*ccss(3) +&
               4.0_dp*ccss(4)*u)*u)*u)*u
          dgcssdrho = x*dgcssds
          dgcssddrho = y*dgcssds

          !     *** g_c(rho_alpha,rho_beta) ***

          gs2 = gamma_cab*s2
          q = 1.0_dp/(1.0_dp + gs2)
          u = gs2*q
          gcab = ccab(0) + (ccab(1) + (ccab(2) + (ccab(3) + ccab(4)*u)*u)*u)*u
          dgcabds = q*(ccab(1) + (2.0_dp*ccab(2) + (3.0_dp*ccab(3) +&
               4.0_dp*ccab(4)*u)*u)*u)*u
          dgcabdrho = x*dgcabds
          dgcabddrho = y*dgcabds

          !     *** Finally collect all contributions ***

          e_0(ii) = e_0(ii)+exss*gxss + ecss*gcss + ecab*gcab
          e_rho(ii) = e_rho(ii)+vxss*gxss + exss*dgxssdrho +&
               vcss*gcss + ecss*dgcssdrho +&
               vcab*gcab + ecab*dgcabdrho
          e_ndrho(ii) = e_ndrho(ii)+ (exss*dgxssddrho + ecss*dgcssddrho + ecab*dgcabddrho)*drho
       END IF
    END DO

  END SUBROUTINE hcth_lda_calc
  !***************************************************************************

!!****f* xc_hcth/hcth120_lda_cpmd_calc *
!!
!!   NAME
!!     hcth_lda_calc
!!
!!   FUNCTION
!!     Calculate the gradient-corrected xc energy and potential
!!     of Hamprecht, Cohen, Tozer, and Handy (HCTH) for a closed shell
!!     density.
!!
!!   NOTES
!!    Literature:- F. A. Hamprecht, A. J. Cohen, D. J. Tozer, and N. C. Handy,
!!                 J. Chem. Phys. 109, 6264 (1998) -> HCTH/93
!!               - A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik,
!!                 J. Chem. Phys. 112, 1670 (2000) -> HCTH/120 and HCTH/147
!!               - A. D. Boese and N. C. Handy,
!!                 J. Chem. Phys. 114, 5497 (2001) -> HCTH/407
!!               - J. P. Perdew and Y. Wang,
!!                 Phys. Rev. B 45, 13244 (1992) -> PW92
!!     (routine shamelessly copied from CPMD, for testing, should be removed soon)
!!
!!   ARGUMENTS
!!     - iparset: the parameter set that should be used (93,120,147,407)
!!     - rho: the density
!!     - norm_drho: the norm of the gradient of the density
!!     - e_0: the value of the functional in that point
!!     - e_rho: the derivative of the functional wrt. rho
!!     - e_ndrho: the derivative of the functional wrt. norm_drho
!!     - epsilon_rho: the cutoff on rho
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 rewamped [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE hcth120_lda_calc(rho,norm_drho,e_0,e_rho,e_ndrho,&
       epsilon_rho,epsilon_drho,npoints,error)
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0, e_rho, e_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    REAL(kind=dp), PARAMETER :: fr83 = 8._dp/3._dp, o3 = 1.0e0_dp/3.0e0_dp, &
      o34 = 4.0e0_dp/3.0e0_dp , small = 1.e-24_dp

    INTEGER                                  :: i, ii
    REAL(kind=dp) :: BYGR, DENAA, DENAB, DENX, DERA1_DRA, DERAB0_DRAB, &
      DEX_DRHO, DFFAA_DRHO, DFFAB_DRHO, DG, DGAA_DGR, DGAA_DRHO, DGAB_DGR, &
      DGAB_DRHO, DGX_DGR, DGX_DRHO, DRA_DRHO, DRAB_DRHO, ERA1, ERAB0, EX, &
      F83RHO, FFAA, FFAB, G, GAA, GAB, GR, GX, myrho, R3PI, R3Q2, RA, RAB, &
      TAA, TAB, TXX, UAA, UAB, UX, XA, XA2
    REAL(kind=dp), DIMENSION(6)              :: caa, cab, cg0, cg1, cx

    r3q2=2.e0_dp**(-o3)
    r3pi=(3.e0_dp/pi)**o3
    !.....coefficients for PW correlation......................................
    cg0(1)= 0.031091e0_dp
    cg0(2)= 0.213700e0_dp
    cg0(3)= 7.595700e0_dp
    cg0(4)= 3.587600e0_dp
    cg0(5)= 1.638200e0_dp
    cg0(6)= 0.492940e0_dp
    cg1(1)= 0.015545e0_dp
    cg1(2)= 0.205480e0_dp
    cg1(3)=14.118900e0_dp
    cg1(4)= 6.197700e0_dp
    cg1(5)= 3.366200e0_dp
    cg1(6)= 0.625170e0_dp
    !......HCTH-19-4.....................................
    caa(1)=  0.489508e+00_dp
    caa(2)= -0.260699e+00_dp
    caa(3)=  0.432917e+00_dp
    caa(4)= -0.199247e+01_dp
    caa(5)=  0.248531e+01_dp
    caa(6)=  0.200000e+00_dp
    cab(1)=  0.514730e+00_dp
    cab(2)=  0.692982e+01_dp
    cab(3)= -0.247073e+02_dp
    cab(4)=  0.231098e+02_dp
    cab(5)= -0.113234e+02_dp
    cab(6)=  0.006000e+00_dp
    cx(1) =  0.109163e+01_dp
    cx(2) = -0.747215e+00_dp
    cx(3) =  0.507833e+01_dp
    cx(4) = -0.410746e+01_dp
    cx(5) =  0.117173e+01_dp
    cx(6)=   0.004000e+00_dp
    !...........................................................................

    DO ii=1,npoints

       IF (rho(ii)> epsilon_rho .AND. norm_drho(ii)> epsilon_drho) THEN
          gr=norm_drho(ii)
          myrho=rho(ii)
          xa=2.0e0_dp**o3*gr/myrho**(4.0e0_dp*o3)
          xa2=xa*xa
          ra=(3.0e0_dp/(2.0e0_dp*pi*myrho))**o3
          rab=r3q2*ra
          dra_drho=-(18d0*pi*myrho**4.0e0_dp)**(-o3)
          drab_drho=r3q2*dra_drho
          CALL pwcorr(ra,cg1,g,dg)
          era1=g
          dera1_dra=dg
          CALL pwcorr(rab,cg0,g,dg)
          erab0=g
          derab0_drab=dg
          ex=-0.75e0_dp*r3pi*myrho**o34
          dex_drho=-r3pi*myrho**o3
          uaa=caa(6)*xa2
          uaa=uaa/(1.0e0_dp+uaa)
          uab=cab(6)*xa2
          uab=uab/(1.0e0_dp+uab)
          ux=cx(6)*xa2
          ux=ux/(1.0e0_dp+ux)
          ffaa=myrho*era1
          ffab=myrho*erab0-ffaa
          dffaa_drho=era1+myrho*dera1_dra*dra_drho
          dffab_drho=erab0+myrho*derab0_drab*drab_drho-dffaa_drho
          gaa=0d0
          gab=0d0
          gx=0d0
          dgaa_drho=0d0
          dgab_drho=0d0
          dgx_drho=0d0
          dgaa_dgr=0d0
          dgab_dgr=0d0
          dgx_dgr=0d0
          !mb-> non i-dependent stuff
          denaa=1.e0_dp/(1.0e0_dp+caa(6)*xa2)
          denab=1.e0_dp/(1.0e0_dp+cab(6)*xa2)
          denx =1.e0_dp/(1.0e0_dp+cx(6)*xa2)
          f83rho=fr83/myrho
          bygr=2.0e0_dp/MAX(gr,small)
          DO i=0,4
             taa=caa(i+1)*uaa**i
             tab=cab(i+1)*uab**i
             txx=cx(i+1)*ux**i
             gaa=gaa+taa
             gab=gab+tab
             gx=gx+txx
             dgaa_drho=dgaa_drho-f83rho*denaa*taa*REAL(i,dp)
             dgab_drho=dgab_drho-f83rho*denab*tab*REAL(i,dp)
             dgx_drho=dgx_drho-f83rho*denx*txx*REAL(i,dp)
             dgaa_dgr=dgaa_dgr+bygr*denaa*taa*REAL(i,dp)
             dgab_dgr=dgab_dgr+bygr*denab*tab*REAL(i,dp)
             dgx_dgr=dgx_dgr+bygr*denx*txx*REAL(i,dp)
          ENDDO
          e_0(ii)=e_0(ii)&
               +ex*gx+ffaa*gaa+ffab*gab
          e_rho(ii)=e_rho(ii)&
               +dex_drho*gx+ex*dgx_drho &
               +dffaa_drho*gaa+ffaa*dgaa_drho &
               +dffab_drho*gab+ffab*dgab_drho
          e_ndrho(ii)=e_ndrho(ii)&
               +(ex*dgx_dgr+ffaa*dgaa_dgr+ffab*dgab_dgr)
       END IF
    END DO

  END SUBROUTINE hcth120_lda_calc

  SUBROUTINE pwcorr(r,c,g,dg)
    REAL(kind=dp), INTENT(in)                :: r
    REAL(kind=dp), DIMENSION(6), INTENT(in)  :: c
    REAL(kind=dp), INTENT(out)               :: g, dg

    REAL(kind=dp)                            :: drb, r12, r2, r32, rb, sb

    r12=dsqrt(r)
    r32=r*r12
    r2=r*r
    rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
    sb=1.0e0_dp+1.0e0_dp/(2.0e0_dp*c(1)*rb)
    g=-2.0e0_dp*c(1)*(1.0e0_dp+c(2)*r)*dlog(sb)
    drb=c(3)/(2.0e0_dp*r12)+c(4)+1.5e0_dp*c(5)*r12+2.0e0_dp*c(6)*r
    dg=(1.0e0_dp+c(2)*r)*drb/(rb**2.0e0_dp*sb)-2.0e0_dp*c(1)*c(2)*dlog(sb)
  END SUBROUTINE pwcorr

!!****f* xc_hcth/hcth120_lda_cpmd_calc *
!!
!!   NAME
!!     hcth_lda_calc
!!
!!   FUNCTION
!!     Calculate the gradient-corrected xc energy and potential
!!     of Hamprecht, Cohen, Tozer, and Handy (HCTH) for a closed shell
!!     density.
!!
!!   NOTES
!!     HCTH/407, Handy et al. JCP 114, p. 5497 (2001)
!!     Present release: Tsukuba, 15/5/2002 Mauro Boero
!!     (routine shamelessly copied from CPMD, for testing, shoud be removed soon)
!!
!!   ARGUMENTS
!!     - iparset: the parameter set that should be used (93,120,147,407)
!!     - rho: the density
!!     - norm_drho: the norm of the gradient of the density
!!     - e_0: the value of the functional in that point
!!     - e_rho: the derivative of the functional wrt. rho
!!     - e_ndrho: the derivative of the functional wrt. norm_drho
!!     - epsilon_rho: the cutoff on rho
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** *********************************************************************
  SUBROUTINE hcth407_lda_calc(rho,norm_drho,e_0,e_rho,e_ndrho,&
       epsilon_rho,epsilon_drho,npoints,error)
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0, e_rho, e_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    REAL(kind=dp), PARAMETER :: fr83 = 8._dp/3._dp, o3 = 1.0e0_dp/3.0e0_dp, &
      o34 = 4.0e0_dp/3.0e0_dp, small = 1.e-24_dp

    INTEGER                                  :: i, ii
    REAL(kind=dp) :: bygr, denaa, denab, denx, dera1_dra, derab0_drab, &
      dex_drho, dffaa_drho, dffab_drho, dg, dgaa_dgr, dgaa_drho, dgab_dgr, &
      dgab_drho, dgx_dgr, dgx_drho, dra_drho, drab_drho, era1, erab0, ex, &
      f83rho, ffaa, ffab, g, gaa, gab, gr, gx, my_rho, r3pi, r3q2, ra, rab, &
      taa, tab, txx, uaa, uab, ux, xa, xa2
    REAL(kind=dp), DIMENSION(6)              :: caa, cab, cg0, cg1, cx

    r3q2=2.e0_dp**(-o3)
    r3pi=(3.e0_dp/pi)**o3
    !.....coefficients for PW correlation...............
    cg0(1)= 0.031091e0_dp
    cg0(2)= 0.213700e0_dp
    cg0(3)= 7.595700e0_dp
    cg0(4)= 3.587600e0_dp
    cg0(5)= 1.638200e0_dp
    cg0(6)= 0.492940e0_dp
    cg1(1)= 0.015545e0_dp
    cg1(2)= 0.205480e0_dp
    cg1(3)=14.118900e0_dp
    cg1(4)= 6.197700e0_dp
    cg1(5)= 3.366200e0_dp
    cg1(6)= 0.625170e0_dp
    !......HCTH/407.....................................
    caa(1)=  1.18777e+00_dp
    caa(2)= -2.40290e+00_dp
    caa(3)=  5.61740e+00_dp
    caa(4)= -9.17920e+00_dp
    caa(5)=  6.24800e+00_dp
    caa(6)=  0.20000e+00_dp
    cab(1)=  0.58908e+00_dp
    cab(2)=  4.42370e+00_dp
    cab(3)=-19.22200e+00_dp
    cab(4)= 42.57200e+00_dp
    cab(5)=-42.00500e+00_dp
    cab(6)=  0.00600e+00_dp
    cx(1) =  1.08184e+00_dp
    cx(2) = -0.51830e+00_dp
    cx(3) =  3.42560e+00_dp
    cx(4) = -2.62900e+00_dp
    cx(5) =  2.28860e+00_dp
    cx(6)=   0.00400e+00_dp
    !...........................................................................

    DO ii=1,npoints

       IF (rho(ii)> epsilon_rho) THEN
          my_rho=rho(ii)
          gr=norm_drho(ii)
          xa=2.0e0_dp**o3*gr/my_rho**(4.0e0_dp*o3)
          xa2=xa*xa
          ra=(3.0e0_dp/(2.0e0_dp*pi*my_rho))**o3
          rab=r3q2*ra
          dra_drho=-(18d0*pi*my_rho**4.0e0_dp)**(-o3)
          drab_drho=r3q2*dra_drho
          CALL pwc407(ra,cg1,g,dg)
          era1=g
          dera1_dra=dg
          CALL pwc407(rab,cg0,g,dg)
          erab0=g
          derab0_drab=dg
          ex=-0.75e0_dp*r3pi*my_rho**o34
          dex_drho=-r3pi*my_rho**o3
          uaa=caa(6)*xa2
          uaa=uaa/(1.0e0_dp+uaa)
          uab=cab(6)*xa2
          uab=uab/(1.0e0_dp+uab)
          ux=cx(6)*xa2
          ux=ux/(1.0e0_dp+ux)
          ffaa=my_rho*era1
          ffab=my_rho*erab0-ffaa
          dffaa_drho=era1+my_rho*dera1_dra*dra_drho
          dffab_drho=erab0+my_rho*derab0_drab*drab_drho-dffaa_drho
          gaa=0._dp
          gab=0._dp
          gx=0._dp
          dgaa_drho=0._dp
          dgab_drho=0._dp
          dgx_drho=0._dp
          dgaa_dgr=0._dp
          dgab_dgr=0._dp
          dgx_dgr=0._dp
          !mb-> non i-dependent stuff
          denaa=1.e0_dp/(1.0e0_dp+caa(6)*xa2)
          denab=1.e0_dp/(1.0e0_dp+cab(6)*xa2)
          denx =1.e0_dp/(1.0e0_dp+cx(6)*xa2)
          f83rho=fr83/my_rho
          bygr=2.0e0_dp/MAX(gr,small)
          DO i=0,4
             taa=caa(i+1)*uaa**i
             tab=cab(i+1)*uab**i
             txx=cx(i+1)*ux**i
             gaa=gaa+taa
             gab=gab+tab
             gx=gx+txx
             dgaa_drho=dgaa_drho-f83rho*denaa*taa*REAL(i,dp)
             dgab_drho=dgab_drho-f83rho*denab*tab*REAL(i,dp)
             dgx_drho=dgx_drho-f83rho*denx*txx*REAL(i,dp)
             dgaa_dgr=dgaa_dgr+bygr*denaa*taa*REAL(i,dp)
             dgab_dgr=dgab_dgr+bygr*denab*tab*REAL(i,dp)
             dgx_dgr=dgx_dgr+bygr*denx*txx*REAL(i,dp)
          ENDDO
          e_0(ii)=e_0(ii)&
               +ex*gx+ffaa*gaa+ffab*gab
          e_rho(ii)=e_rho(ii)&
               +dex_drho*gx+ex*dgx_drho &
               +dffaa_drho*gaa+ffaa*dgaa_drho &
               +dffab_drho*gab+ffab*dgab_drho
          e_ndrho(ii)=e_ndrho(ii)&
               +ex*dgx_dgr+ffaa*dgaa_dgr+ffab*dgab_dgr
       END IF
    END DO
  END SUBROUTINE hcth407_lda_calc
! =-------------------------------------------------------------------=
  SUBROUTINE pwc407(r,c,g,dg)
    REAL(kind=dp), INTENT(in)                :: r
    REAL(kind=dp), DIMENSION(6), INTENT(in)  :: c
    REAL(kind=dp), INTENT(out)               :: g, dg

    REAL(kind=dp)                            :: drb, r12, r2, r32, rb, sb

    r12=dsqrt(r)
    r32=r*r12
    r2=r*r
    rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
    sb=1.0e0_dp+1.0e0_dp/(2.0e0_dp*c(1)*rb)
    g=-2.0e0_dp*c(1)*(1.0e0_dp+c(2)*r)*dlog(sb)
    drb=c(3)/(2.0e0_dp*r12)+c(4)+1.5e0_dp*c(5)*r12+2.0e0_dp*c(6)*r
    dg=(1.0e0_dp+c(2)*r)*drb/(rb**2.0e0_dp*sb)-2.0e0_dp*c(1)*c(2)*dlog(sb)
  END SUBROUTINE pwc407
!-----------------------------------------------------------------------------

END MODULE xc_hcth
