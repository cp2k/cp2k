!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief This functional is a combination of three different exchange hole
!>        models. The ingredients are:
!>
!>          1. Becke Roussel exchange hole
!>          2. PBE exchange hole
!>          3. LDA exchange hole
!>
!>        The full functionals is given as follows
!>
!>        Fx    = eps_lr_lda/eps_lr_br
!>        Fcorr = alpha/( exp( (Fx-mu)/N ) + 1)
!>        rhox  = Fcorr * eps_lr_pbe + (1-Fcorr)*eps_lr_br
!>        eps   = int_{R}^{\infty} rhox*s*ds
!>
!>        with alpha, mu and N fitting parameters
!> \par History
!>      01.2009 created [mguidon]
!> \author mguidon
! *****************************************************************************

MODULE xc_xbr_pbe_lda_hole_t_c_lr
  
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
  USE xc_xbecke_roussel,               ONLY: x_br_lsd_y_gt_0,&
                                             x_br_lsd_y_gt_0_cutoff,&
                                             x_br_lsd_y_lte_0,&
                                             x_br_lsd_y_lte_0_cutoff
  USE xc_xlda_hole_t_c_lr,             ONLY: xlda_hole_t_c_lr_lda_calc_0
  USE xc_xpbe_hole_t_c_lr,             ONLY: xpbe_hole_t_c_lr_lda_calc_1,&
                                             xpbe_hole_t_c_lr_lda_calc_2
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xbr_pbe_lda_hole_t_c_lr'

  REAL(dp), PARAMETER, PRIVATE  :: br_a1 =  1.5255251812009530_dp,&
                                   br_a2 =  0.4576575543602858_dp,&
                                   br_a3 =  0.4292036732051034_dp,&
                                   br_c0 =  0.7566445420735584_dp,&
                                   br_c1 = -2.6363977871370960_dp,&
                                   br_c2 =  5.4745159964232880_dp,&
                                   br_c3 = -12.657308127108290_dp,&
                                   br_c4 =  4.1250584725121360_dp,&
                                   br_c5 = -30.425133957163840_dp,&
                                   br_b0 =  0.4771976183772063_dp,&
                                   br_b1 = -1.7799813494556270_dp,&
                                   br_b2 =  3.8433841862302150_dp,&
                                   br_b3 = -9.5912050880518490_dp,&
                                   br_b4 =  2.1730180285916720_dp,&
                                   br_b5 = -30.425133851603660_dp,&
                                   br_d0 =  0.00004435009886795587_dp,&
                                   br_d1 =  0.58128653604457910_dp,&
                                   br_d2 =  66.742764515940610_dp,&
                                   br_d3 =  434.26780897229770_dp,&
                                   br_d4 =  824.7765766052239000_dp,&
                                   br_d5 =  1657.9652731582120_dp,&
                                   br_e0 =  0.00003347285060926091_dp,&
                                   br_e1 =  0.47917931023971350_dp,&
                                   br_e2 =  62.392268338574240_dp,&
                                   br_e3 =  463.14816427938120_dp,&
                                   br_e4 =  785.2360350104029000_dp,&
                                   br_e5 =  1657.962968223273000000_dp,&
                                   br_BB = 2.085749716493756_dp

  REAL(dp), PARAMETER, PRIVATE  :: smax = 8.572844_dp,&
                                   scutoff = 8.3_dp,&
                                   sconst = 18.79622316_dp,&
                                   gcutoff = 0.08_dp

  REAL(dp), PARAMETER, PRIVATE  :: alpha = 0.3956891_dp,&
                                   N = -0.0009800242_dp,&
                                   mu = 0.00118684_dp

  PUBLIC :: xbr_pbe_lda_hole_tc_lr_lda_info, &
            xbr_pbe_lda_hole_tc_lr_lsd_info, &
            xbr_pbe_lda_hole_tc_lr_lda_eval, &
            xbr_pbe_lda_hole_tc_lr_lsd_eval
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
!> \author mguidon (01.2009)
! *****************************************************************************
  SUBROUTINE xbr_pbe_lda_hole_tc_lr_lda_info(reference,shortform, needs, max_deriv,&
       error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'xbr_pbe_lda_hole_tc_lr_lda_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF ( PRESENT ( reference ) ) THEN
      reference = "{LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "{LDA}"
    END IF

    IF (PRESENT(needs)) THEN
      needs%rho=.TRUE.
      needs%norm_drho=.TRUE.
      needs%tau = .TRUE.
      needs%laplace_rho = .TRUE.
    END IF

    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xbr_pbe_lda_hole_tc_lr_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author mguidon (01.2009)
! *****************************************************************************
  SUBROUTINE xbr_pbe_lda_hole_tc_lr_lsd_info(reference,shortform, needs, max_deriv,&
     error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'xbr_pbe_lda_hole_tc_lr_lsd_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF ( PRESENT ( reference ) ) THEN
      reference = "{LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "{LDA}"
    END IF

    IF (PRESENT(needs)) THEN
      needs%rho_spin =.TRUE.
      needs%norm_drho_spin =.TRUE.
      needs%tau_spin = .TRUE.
      needs%laplace_rho_spin = .TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xbr_pbe_lda_hole_tc_lr_lsd_info

! *****************************************************************************
!> \brief Intermediate routine that gets grids, derivatives and some params
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param params parameters for functional
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author mguidon (01.2009)
! *****************************************************************************
  SUBROUTINE xbr_pbe_lda_hole_tc_lr_lda_eval(rho_set,deriv_set,grad_deriv,params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'xbr_pbe_lda_hole_tc_lr_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: gamma, R, sx
    REAL(kind=dp)                            :: epsilon_norm_drho, epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_laplace_rho, &
                                                e_ndrho, e_rho, e_tau, &
                                                laplace_rho, norm_drho, rho, &
                                                tau
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL xc_rho_set_get(rho_set, rho=rho, norm_drho=norm_drho, &
           tau=tau,laplace_rho=laplace_rho,local_bounds=bo,&
           rho_cutoff=epsilon_rho,drho_cutoff=epsilon_norm_drho,error=error)
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
      e_tau => dummy
      e_laplace_rho => dummy

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
        deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_tau,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho,error=error)
      END IF
      IF (grad_deriv>1.OR.grad_deriv<-1) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 1 not implemented", &
             error=error, error_level=cp_failure_level)
      END IF

      CALL section_vals_val_get(params,"scale_x",r_val=sx,error=error)
      CALL section_vals_val_get(params,"CUTOFF_RADIUS",r_val=R,error=error)
      CALL section_vals_val_get(params,"GAMMA",r_val=gamma,error=error)

      IF ( R == 0.0_dp ) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
                       message="Cutoff_Radius 0.0 not implemented", &
                       error=error, error_level=cp_failure_level)
      END IF

      !$omp parallel default(none) &
      !$omp          shared(rho, norm_drho, laplace_rho, tau, e_0, e_rho) &
      !$omp          shared(e_ndrho, e_tau, e_laplace_rho, grad_deriv) &
      !$omp          shared(npoints, epsilon_rho) &
      !$omp          shared(epsilon_norm_drho, sx, r, gamma, error)

      CALL xbr_pbe_lda_hole_tc_lr_lda_calc(rho=rho, norm_drho=norm_drho,&
          laplace_rho=laplace_rho,tau=tau,e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
          e_tau=e_tau,e_laplace_rho=e_laplace_rho,grad_deriv=grad_deriv, &
          npoints=npoints,epsilon_rho=epsilon_rho,&
          epsilon_norm_drho=epsilon_norm_drho,sx=sx,R=R,gamma=gamma,error=error)

      !$omp end parallel

      IF (cp_debug) THEN
         DEALLOCATE(dummy,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ELSE
         NULLIFY(dummy)
      END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE xbr_pbe_lda_hole_tc_lr_lda_eval

! *****************************************************************************
!> \brief Low level routine that calls the three involved holes and puts them
!>        together
!> \param rho values on the grid
!> \param norm_drho values on the grid
!> \param laplace_rho values on the grid
!> \param tau values on the grid
!> \param e_0 derivatives on the grid
!> \param e_rho derivatives on the grid
!> \param e_ndrho derivatives on the grid
!> \param e_tau derivatives on the grid
!> \param e_laplace_rho derivatives on the grid
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of gridpoints
!> \param epsilon_rho cutoffs
!> \param epsilon_norm_drho cutoffs
!> \param sx parameters for  functional
!> \param R parameters for  functional
!> \param gamma parameters for  functional
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author mguidon (01.2009)
! *****************************************************************************
  SUBROUTINE xbr_pbe_lda_hole_tc_lr_lda_calc(rho,norm_drho,laplace_rho,tau,e_0,e_rho,&
                                     e_ndrho,e_tau,e_laplace_rho,grad_deriv,npoints,&
                                     epsilon_rho,epsilon_norm_drho,sx,R,gamma,error)

    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: e_laplace_rho, e_tau, &
                                                e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: tau, laplace_rho, norm_drho, &
                                                rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, sx, R, &
                                                gamma
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'xbr_pbe_lda_hole_tc_lr_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp) :: dFermi_dlaplace_rho, dFermi_dndrho, dFermi_drho, dFermi_dtau, &
      e_0_br, e_0_lda, e_0_pbe, e_laplace_rho_br, e_ndrho_br, e_ndrho_pbe, &
      e_rho_br, e_rho_lda, e_rho_pbe, e_tau_br, Fermi, Fx, my_laplace_rho, &
      my_ndrho, my_rho, my_tau, ss, ss2, sscale, t1, t15, t16, t2, t3, t4, &
      t5, t6, t7, t8, t9, yval

    !$omp do

    DO ip = 1,npoints
      my_rho = 0.5_dp*MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho = 0.5_dp*MAX(norm_drho(ip),EPSILON(0.0_dp)*1.e4_dp)
        my_tau = 0.5_dp*MAX(EPSILON(0.0_dp)*1.e4_dp,tau(ip))
        my_laplace_rho = 0.5_dp*laplace_rho(ip)

        ! ** We calculate first the Becke-Roussel part, saving everything in local variables
        t1 = pi ** (0.1e1_dp / 0.3e1_dp)
        t2 = t1 ** 2
        t3 = my_rho ** (0.1e1_dp / 0.3e1_dp)
        t4 = t3 ** 2
        t5 = t4 * my_rho
        t8 = my_ndrho ** 2
        t9 = 0.1e1_dp / my_rho
        t15 = my_laplace_rho / 0.6e1_dp - gamma * (2.0_dp*my_tau - t8 * t9 / 0.4e1_dp) / 0.3e1_dp
        t16 = 0.1e1_dp / t15
        yval = 0.2e1_dp / 0.3e1_dp * t2 * t5 * t16

        e_0_br = 0.0_dp
        e_rho_br = 0.0_dp
        e_ndrho_br = 0.0_dp
        e_tau_br = 0.0_dp
        e_laplace_rho_br = 0.0_dp

        IF( R == 0.0_dp ) THEN
          IF( yval <= 0.0_dp) THEN
            CALL x_br_lsd_y_lte_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                  e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                  sx, R, gamma, grad_deriv, error)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_0_br = 2.0_dp * e_0_br
          ELSE
            CALL x_br_lsd_y_gt_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                 e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                 sx, R, gamma, grad_deriv, error)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_0_br = 2.0_dp * e_0_br
          END IF
        ELSE
          IF( yval <= 0.0_dp) THEN
            CALL x_br_lsd_y_lte_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                         e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                         sx, R, gamma, grad_deriv, error)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_0_br = 2.0_dp * e_0_br
          ELSE
            CALL x_br_lsd_y_gt_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                        e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                        sx, R, gamma, grad_deriv, error)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_0_br = 2.0_dp * e_0_br
          END IF
        END IF


        ! ** Now we calculate the pbe cutoff part
        ! ** Attention we need to scale rho, ndrho first
        my_rho = my_rho * 2.0_dp
        my_ndrho = my_ndrho * 2.0_dp

        ! ** Do some precalculation in order to catch the correct branch afterwards
        sscale = 1.0_dp
        t1 = pi ** 2
        t2 = t1 * my_rho
        t3 = t2 ** (0.1e1_dp / 0.3e1_dp)
        t4 = 0.1e1_dp / t3
        t6 = my_ndrho * t4
        t7 = 0.1e1_dp / my_rho
        t8 = t7 * sscale
        ss = 0.3466806371753173524216762e0_dp * t6 * t8
        IF( ss > scutoff) THEN
          ss2 = ss*ss
          sscale = (smax*ss2-sconst)/(ss2*ss)
        END IF
        e_0_pbe = 0.0_dp
        e_rho_pbe = 0.0_dp
        e_ndrho_pbe = 0.0_dp
        IF(ss*sscale>gcutoff) THEN
          CALL xpbe_hole_t_c_lr_lda_calc_1(e_0_pbe, e_rho_pbe, e_ndrho_pbe,&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, grad_deriv)
        ELSE
          CALL xpbe_hole_t_c_lr_lda_calc_2(e_0_pbe, e_rho_pbe, e_ndrho_pbe,&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, grad_deriv)
        END IF

        ! ** Finally we get the LDA part

        e_0_lda = 0.0_dp
        e_rho_lda = 0.0_dp
        CALL xlda_hole_t_c_lr_lda_calc_0(grad_deriv, my_rho, e_0_lda, e_rho_lda,&
                                         sx, R, error)


        Fx = e_0_br/e_0_lda

        Fermi = alpha/(EXP( (Fx-mu)/N ) + 1.0_dp)

        dFermi_drho = -Fermi**2/alpha/N*(e_rho_br/e_0_lda-e_0_br*e_rho_lda/e_0_lda**2)*EXP((Fx-mu)/N)
        dFermi_dndrho = -Fermi**2/alpha/N*(e_ndrho_br/e_0_lda)*EXP((Fx-mu)/N)
        dFermi_dtau = -Fermi**2/alpha/N*(e_tau_br/e_0_lda)*EXP((Fx-mu)/N)
        dFermi_dlaplace_rho = -Fermi**2/alpha/N*(e_laplace_rho_br/e_0_lda)*EXP((Fx-mu)/N)


        e_0(ip) = e_0(ip) + ( Fermi * e_0_pbe + (1.0_dp-Fermi) * e_0_br ) * sx

        IF( grad_deriv >=1 .OR. grad_deriv == -1) THEN

          e_rho(ip) = e_rho(ip) + ( Fermi * e_rho_pbe + dFermi_drho * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_rho_br - dFermi_drho*e_0_br ) * sx

          e_ndrho(ip) = e_ndrho(ip) + ( Fermi * e_ndrho_pbe + dFermi_dndrho * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_ndrho_br - dFermi_dndrho*e_0_br ) * sx

          e_tau(ip) = e_tau(ip) + ( dFermi_dtau * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_tau_br - dFermi_dtau*e_0_br ) * sx

          e_laplace_rho(ip) = e_laplace_rho(ip) + ( dFermi_dlaplace_rho * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_laplace_rho_br - dFermi_dlaplace_rho*e_0_br ) * sx
        END IF


      END IF
    END DO

    !$omp end do

  END SUBROUTINE xbr_pbe_lda_hole_tc_lr_lda_calc

! *****************************************************************************
!> \brief Intermediate routine that gets grids, derivatives and some params
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param params parameters for functional
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author mguidon (01.2009)
! *****************************************************************************
  SUBROUTINE xbr_pbe_lda_hole_tc_lr_lsd_eval(rho_set,deriv_set,grad_deriv,params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'xbr_pbe_lda_hole_tc_lr_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: gamma, R, sx
    REAL(kind=dp)                            :: epsilon_norm_drho, epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rhoa, &
      e_laplace_rhob, e_ndrhoa, e_ndrhob, e_rhoa, e_rhob, e_tau_a, e_tau_b, &
      laplace_rhoa, laplace_rhob, norm_drhoa, norm_drhob, rhoa, rhob, tau_a, &
      tau_b
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL xc_rho_set_get(rho_set, rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
           norm_drhob=norm_drhob, tau_a=tau_a, tau_b=tau_b, laplace_rhoa=laplace_rhoa,&
           laplace_rhob=laplace_rhob,local_bounds=bo,&
           rho_cutoff=epsilon_rho,drho_cutoff=epsilon_norm_drho,error=error)
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
      e_rhoa => dummy
      e_rhob => dummy
      e_ndrhoa => dummy
      e_ndrhob => dummy
      e_tau_a => dummy
      e_tau_b => dummy
      e_laplace_rhoa => dummy
      e_laplace_rhob => dummy

      IF (grad_deriv>=0) THEN
        deriv => xc_dset_get_derivative(deriv_set,"",&
             allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
      END IF
      IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rhoa,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rhob,error=error)

        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrhob,error=error)

        deriv => xc_dset_get_derivative(deriv_set,"(tau_a)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_tau_a,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(tau_b)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_tau_b,error=error)

        deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob,error=error)
      END IF
      IF (grad_deriv>1.OR.grad_deriv<-1) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 1 not implemented", &
             error=error, error_level=cp_failure_level)
      END IF

      CALL section_vals_val_get(params,"scale_x",r_val=sx,error=error)
      CALL section_vals_val_get(params,"CUTOFF_RADIUS",r_val=R,error=error)
      CALL section_vals_val_get(params,"GAMMA",r_val=gamma,error=error)

      IF ( R == 0.0_dp ) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
                       message="Cutoff_Radius 0.0 not implemented", &
                       error=error, error_level=cp_failure_level)
      END IF

      !$omp parallel default(none) &
      !$omp          shared(rhoa, norm_drhoa, laplace_rhoa, tau_a, e_0) &
      !$omp          shared(e_rhoa, e_ndrhoa, e_tau_a, e_laplace_rhoa) &
      !$omp          shared(grad_deriv, npoints, epsilon_rho) &
      !$omp          shared(epsilon_norm_drho, sx, r, gamma, error) &
      !$omp          shared(rhob, norm_drhob, laplace_rhob, tau_b, e_rhob) &
      !$omp          shared(e_ndrhob, e_tau_b, e_laplace_rhob)

      CALL xbr_pbe_lda_hole_tc_lr_lsd_calc(rho=rhoa, norm_drho=norm_drhoa,&
          laplace_rho=laplace_rhoa,tau=tau_a,e_0=e_0,e_rho=e_rhoa,e_ndrho=e_ndrhoa,&
          e_tau=e_tau_a,e_laplace_rho=e_laplace_rhoa,grad_deriv=grad_deriv, &
          npoints=npoints,epsilon_rho=epsilon_rho,&
          epsilon_norm_drho=epsilon_norm_drho,sx=sx,R=R,gamma=gamma,error=error)

      CALL xbr_pbe_lda_hole_tc_lr_lsd_calc(rho=rhob, norm_drho=norm_drhob,&
          laplace_rho=laplace_rhob,tau=tau_b,e_0=e_0,e_rho=e_rhob,e_ndrho=e_ndrhob,&
          e_tau=e_tau_b,e_laplace_rho=e_laplace_rhob,grad_deriv=grad_deriv, &
          npoints=npoints,epsilon_rho=epsilon_rho,&
          epsilon_norm_drho=epsilon_norm_drho,sx=sx,R=R,gamma=gamma,error=error)

      !$omp end parallel

      IF (cp_debug) THEN
         DEALLOCATE(dummy,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ELSE
         NULLIFY(dummy)
      END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE xbr_pbe_lda_hole_tc_lr_lsd_eval
! *****************************************************************************
!> \brief Low level routine that calls the three involved holes and puts them
!>        together
!> \param rho values on the grid
!> \param norm_drho values on the grid
!> \param laplace_rho values on the grid
!> \param tau values on the grid
!> \param e_0 derivatives on the grid
!> \param e_rho derivatives on the grid
!> \param e_ndrho derivatives on the grid
!> \param e_tau derivatives on the grid
!> \param e_laplace_rho derivatives on the grid
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of gridpoints
!> \param epsilon_rho cutoffs
!> \param epsilon_norm_drho cutoffs
!> \param sx parameters for  functional 
!> \param R parameters for  functional 
!> \param gamma parameters for  functional 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author mguidon (01.2009)
! *****************************************************************************
  SUBROUTINE xbr_pbe_lda_hole_tc_lr_lsd_calc(rho,norm_drho,laplace_rho,tau,e_0,e_rho,&
                                     e_ndrho,e_tau,e_laplace_rho,grad_deriv,npoints,&
                                     epsilon_rho,epsilon_norm_drho,sx,R,gamma,error)

    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: e_laplace_rho, e_tau, &
                                                e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: tau, laplace_rho, norm_drho, &
                                                rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, sx, R, &
                                                gamma
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'xbr_pbe_lda_hole_tc_lr_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp) :: dFermi_dlaplace_rho, dFermi_dndrho, dFermi_drho, dFermi_dtau, &
      e_0_br, e_0_lda, e_0_pbe, e_laplace_rho_br, e_ndrho_br, e_ndrho_pbe, &
      e_rho_br, e_rho_lda, e_rho_pbe, e_tau_br, Fermi, Fx, my_laplace_rho, &
      my_ndrho, my_rho, my_tau, ss, ss2, sscale, t1, t15, t16, t2, t3, t4, &
      t5, t6, t7, t8, t9, yval

    !$omp do

    DO ip = 1,npoints
      my_rho = MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho = MAX(norm_drho(ip),EPSILON(0.0_dp)*1.e4_dp)
        my_tau=1.0_dp*MAX(EPSILON(0.0_dp)*1.e4_dp,tau(ip))
        my_laplace_rho = 1.0_dp*laplace_rho(ip)

        t1 = pi ** (0.1e1_dp / 0.3e1_dp)
        t2 = t1 ** 2
        t3 = my_rho ** (0.1e1_dp / 0.3e1_dp)
        t4 = t3 ** 2
        t5 = t4 * my_rho
        t8 = my_ndrho ** 2
        t9 = 0.1e1_dp / my_rho
        t15 = my_laplace_rho / 0.6e1_dp - gamma * (2.0_dp*my_tau - t8 * t9 / 0.4e1_dp) / 0.3e1_dp
        t16 = 0.1e1_dp / t15
        yval = 0.2e1_dp / 0.3e1_dp * t2 * t5 * t16

        e_0_br = 0.0_dp
        e_rho_br = 0.0_dp
        e_ndrho_br = 0.0_dp
        e_tau_br = 0.0_dp
        e_laplace_rho_br = 0.0_dp

        IF( R == 0.0_dp ) THEN
          IF( yval <= 0.0_dp) THEN
            CALL x_br_lsd_y_lte_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                  e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                  sx, R, gamma, grad_deriv, error)
          ELSE
            CALL x_br_lsd_y_gt_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                 e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                 sx, R, gamma, grad_deriv, error)
          END IF
        ELSE
          IF( yval <= 0.0_dp) THEN
            CALL x_br_lsd_y_lte_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                         e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                         sx, R, gamma, grad_deriv, error)
          ELSE
            CALL x_br_lsd_y_gt_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0_br, &
                                        e_rho_br, e_ndrho_br, e_tau_br, e_laplace_rho_br,&
                                        sx, R, gamma, grad_deriv, error)
          END IF
        END IF

        ! ** Now we calculate the pbe cutoff part
        ! ** Attention we need to scale rho, ndrho first
        my_rho = my_rho * 2.0_dp
        my_ndrho = my_ndrho * 2.0_dp

        ! ** Do some precalculation in order to catch the correct branch afterwards
        sscale = 1.0_dp
        t1 = pi ** 2
        t2 = t1 * my_rho
        t3 = t2 ** (0.1e1_dp / 0.3e1_dp)
        t4 = 0.1e1_dp / t3
        t6 = my_ndrho * t4
        t7 = 0.1e1_dp / my_rho
        t8 = t7 * sscale
        ss = 0.3466806371753173524216762e0_dp * t6 * t8
        IF( ss > scutoff) THEN
          ss2 = ss*ss
          sscale = (smax*ss2-sconst)/(ss2*ss)
        END IF
        e_0_pbe = 0.0_dp
        e_rho_pbe = 0.0_dp
        e_ndrho_pbe = 0.0_dp
        IF(ss*sscale>gcutoff) THEN
          CALL xpbe_hole_t_c_lr_lda_calc_1(e_0_pbe, e_rho_pbe, e_ndrho_pbe,&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, grad_deriv)
        ELSE
          CALL xpbe_hole_t_c_lr_lda_calc_2(e_0_pbe, e_rho_pbe, e_ndrho_pbe,&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, grad_deriv)
        END IF

        e_0_pbe = 0.5_dp * e_0_pbe


        ! ** Finally we get the LDA part

        e_0_lda = 0.0_dp
        e_rho_lda = 0.0_dp
        CALL xlda_hole_t_c_lr_lda_calc_0(grad_deriv, my_rho, e_0_lda, e_rho_lda,&
                                         sx, R, error)
        e_0_lda = 0.5_dp * e_0_lda

        Fx = e_0_br/e_0_lda

        Fermi = alpha/(EXP( (Fx-mu)/N ) + 1.0_dp)

        dFermi_drho = -Fermi**2/alpha/N*(e_rho_br/e_0_lda-e_0_br*e_rho_lda/e_0_lda**2)*EXP((Fx-mu)/N)
        dFermi_dndrho = -Fermi**2/alpha/N*(e_ndrho_br/e_0_lda)*EXP((Fx-mu)/N)
        dFermi_dtau = -Fermi**2/alpha/N*(e_tau_br/e_0_lda)*EXP((Fx-mu)/N)
        dFermi_dlaplace_rho = -Fermi**2/alpha/N*(e_laplace_rho_br/e_0_lda)*EXP((Fx-mu)/N)


        e_0(ip) = e_0(ip) + ( Fermi * e_0_pbe + (1.0_dp-Fermi) * e_0_br ) * sx

        IF( grad_deriv >=1 .OR. grad_deriv == -1) THEN

          e_rho(ip) = e_rho(ip) + ( Fermi * e_rho_pbe + dFermi_drho * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_rho_br - dFermi_drho*e_0_br ) * sx

          e_ndrho(ip) = e_ndrho(ip) + ( Fermi * e_ndrho_pbe + dFermi_dndrho * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_ndrho_br - dFermi_dndrho*e_0_br ) * sx

          e_tau(ip) = e_tau(ip) + ( dFermi_dtau * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_tau_br - dFermi_dtau*e_0_br ) * sx

          e_laplace_rho(ip) = e_laplace_rho(ip) + ( dFermi_dlaplace_rho * e_0_pbe + &
                                    (1.0_dp - Fermi) * e_laplace_rho_br - dFermi_dlaplace_rho*e_0_br ) * sx
       END IF

      END IF
    END DO

    !$omp end do

  END SUBROUTINE xbr_pbe_lda_hole_tc_lr_lsd_calc

END MODULE xc_xbr_pbe_lda_hole_t_c_lr
