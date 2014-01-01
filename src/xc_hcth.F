!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculate the Hamprecht, Cohen, Tozer, and Handy (HCTH) exchange
!>      functional
!> \author fawzi
! *****************************************************************************
MODULE xc_hcth
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_hcth'

  PUBLIC :: hcth_lda_info, hcth_lda_eval
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE hcth_lda_info(iparset,reference,shortform, needs, max_deriv, error)
    INTEGER, INTENT(in)                      :: iparset
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hcth_lda_info', &
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

! *****************************************************************************
!> \brief evaluates the hcth functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param iparset the parameter set that should be used (93,120,147,407)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE hcth_lda_eval(iparset,rho_set,deriv_set,grad_deriv,error)
    INTEGER, INTENT(in)                      :: iparset
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

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

       CALL hcth_lda_calc(iparset=iparset,rho=rho, norm_drho=norm_drho,&
            e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
            npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_drho=epsilon_drho, error=error)
    END IF
  END SUBROUTINE hcth_lda_eval

! *****************************************************************************
!> \brief Calculate the gradient-corrected xc energy and potential
!>      of Hamprecht, Cohen, Tozer, and Handy (HCTH) for a closed shell
!>      density.
!> \param iparset the parameter set that should be used (93,120,147,407)
!> \param rho the density
!> \param norm_drho the norm of the gradient of the density
!> \param e_ 0: the value of the functional in that point
!> \param e_rho the derivative of the functional wrt. rho
!> \param e_ndrho the derivative of the functional wrt. norm_drho
!> \param epsilon_rho the cutoff on rho
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>     Literature:- F. A. Hamprecht, A. J. Cohen, D. J. Tozer, and N. C. Handy,
!>                  J. Chem. Phys. 109, 6264 (1998) -> HCTH/93
!>                - A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik,
!>                  J. Chem. Phys. 112, 1670 (2000) -> HCTH/120 and HCTH/147
!>                - A. D. Boese and N. C. Handy,
!>                  J. Chem. Phys. 114, 5497 (2001) -> HCTH/407
!>                - J. P. Perdew and Y. Wang,
!>                  Phys. Rev. B 45, 13244 (1992) -> PW92
!> \author fawzi (copying from the routine of Matthias Krack in functionals.F)
! *****************************************************************************
  SUBROUTINE hcth_lda_calc(iparset,rho,norm_drho,e_0,e_rho,e_ndrho,&
       epsilon_rho,epsilon_drho,npoints,error)
    INTEGER, INTENT(IN)                      :: iparset
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0, e_rho, e_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(IN)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hcth_lda_calc', &
      routineP = moduleN//':'//routineN
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
    !$omp         ccab,cx_vwn_e, cx_vwn_v, rsfac, two13,epsilon_rho,npoints, &
    !$omp         e_0,e_rho,e_ndrho)&
    !$omp       private(ii, dgcabddrho, dgcabdrho, dgcabds, dgcssddrho, &
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

END MODULE xc_hcth
