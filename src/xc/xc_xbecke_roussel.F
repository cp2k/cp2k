!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the exchange energy based on the Becke-Roussel exchange
!>        hole. Takes advantage of an analytical representation of the hole
!>        in order to avoid solving a non-linear equation by means of Newton-
!>        Raphson algorithm
!> \note
!> \par History
!>      12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
MODULE xc_xbecke_roussel
  USE bibliography,                    ONLY: BeckeRoussel1989,&
                                             Proynov2007,&
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xbecke_roussel'

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


  PUBLIC :: xbecke_roussel_lda_info, xbecke_roussel_lsd_info, xbecke_roussel_lda_eval, &
            xbecke_roussel_lsd_eval, x_br_lsd_y_lte_0, x_br_lsd_y_gt_0,x_br_lsd_y_lte_0_cutoff, &
            x_br_lsd_y_gt_0_cutoff
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv controls the number of derivatives
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xbecke_roussel_lda_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'xbecke_roussel_lda_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL cite_reference(BeckeRoussel1989)
    CALL cite_reference(Proynov2007)
    IF ( PRESENT ( reference ) ) THEN
      reference = "A.D. Becke, M.R. Roussel, "//&
                  "Phys. Rev. A, vol. 39, n 8, pp. 3761-3767, (1989)"//&
                  "{spin unpolarized}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Becke-Roussel exchange hole (spin unpolarized)"
    END IF
   IF (PRESENT(needs)) THEN
      needs%rho=.TRUE.
      needs%norm_drho=.TRUE.
      needs%tau = .TRUE.
      needs%laplace_rho = .TRUE.
    END IF

    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xbecke_roussel_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xbecke_roussel_lsd_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'xbecke_roussel_lsd_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL cite_reference(BeckeRoussel1989)
    CALL cite_reference(Proynov2007)
    IF ( PRESENT ( reference ) ) THEN
      reference = "A.D. Becke, M.R. Roussel, "//&
                  "Phys. Rev. A, vol. 39, n 8, pp. 3761-3767, (1989)"//&
                  "{spin polarized}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Becke-Roussel exchange hole (spin polarized)"
    END IF
    IF (PRESENT(needs)) THEN
      needs%rho_spin =.TRUE.
      needs%norm_drho_spin =.TRUE.
      needs%tau_spin = .TRUE.
      needs%laplace_rho_spin = .TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xbecke_roussel_lsd_info

! *****************************************************************************
!> \brief evaluates the Becke Roussel exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param br_params parameters for the becke roussel functional
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xbecke_roussel_lda_eval(rho_set,deriv_set,grad_deriv,br_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: br_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xbecke_roussel_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: gamma, R, sx
    REAL(kind=dp)                            :: epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_laplace_rho, &
                                                e_ndrho, e_rho, e_tau, &
                                                laplace_rho, norm_drho, rho, &
                                                tau
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,failure)
    CALL xc_rho_set_get(rho_set, rho=rho, norm_drho=norm_drho, &
         tau=tau,laplace_rho=laplace_rho,local_bounds=bo,&
         rho_cutoff=epsilon_rho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)
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
      deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_tau)
      deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho)
    END IF
    IF (grad_deriv>1.OR.grad_deriv<-1) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
           message="derivatives bigger than 1 not implemented", &
           error_level=cp_failure_level)
    END IF

    CALL section_vals_val_get(br_params,"scale_x",r_val=sx)
    CALL section_vals_val_get(br_params,"CUTOFF_RADIUS",r_val=R)
    CALL section_vals_val_get(br_params,"GAMMA",r_val=gamma)

    !$omp parallel default(none) &
    !$omp          shared(rho, norm_drho, laplace_rho, tau, e_0, e_rho) &
    !$omp          shared(e_ndrho, e_tau, e_laplace_rho, grad_deriv) &
    !$omp          shared(npoints, epsilon_rho) &
    !$omp          shared(sx, r, gamma)

    CALL xbecke_roussel_lda_calc(rho=rho, norm_drho=norm_drho,&
        laplace_rho=laplace_rho,tau=tau,e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
        e_tau=e_tau,e_laplace_rho=e_laplace_rho,grad_deriv=grad_deriv, &
        npoints=npoints,epsilon_rho=epsilon_rho,&
        sx=sx,R=R,gamma=gamma)

    !$omp end parallel

    IF (cp_debug) THEN
       DEALLOCATE(dummy,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
    ELSE
       NULLIFY(dummy)
    END IF
    CALL timestop(handle)
  END SUBROUTINE xbecke_roussel_lda_eval

! *****************************************************************************
!> \brief Precalculates which branch of the code has to be taken
!>        There are two main branches, one for a truncated potential and a cutoff
!>        radius, the other for the full coulomb interaction. In the end, the code
!>        for the lsd part will be called and the lda part is obtained via spin
!>        scaling relations
!> \param rho grid values
!> \param norm_drho grid values
!> \param laplace_rho grid values
!> \param tau grid values
!> \param e_0 energies and derivatives
!> \param e_rho energies and derivatives
!> \param e_ndrho energies and derivatives
!> \param e_tau energies and derivatives
!> \param e_laplace_rho energies and derivatives
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints size of the grids
!> \param epsilon_rho cutoffs
!> \param sx scales the exchange potential and energies
!> \param R cutoff Radius for truncated case
!> \param gamma paramter from original publication, usually set to 1
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xbecke_roussel_lda_calc(rho,norm_drho,laplace_rho,tau,e_0,e_rho,&
                                     e_ndrho,e_tau,e_laplace_rho,grad_deriv,npoints,&
                                     epsilon_rho,sx,R,gamma)

    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: e_laplace_rho, e_tau, &
                                                e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: tau, laplace_rho, norm_drho, &
                                                rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, R, gamma

    CHARACTER(len=*), PARAMETER :: routineN = 'xbecke_roussel_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: e_diff, e_old, &
                                                my_laplace_rho, my_ndrho, &
                                                my_rho, my_tau, t1, t15, t16, &
                                                t2, t3, t4, t5, t8, t9, yval

! Precalculate y, in order to chose the correct branch afterwards

    !$omp do

    DO ip = 1,npoints
      my_rho = 0.5_dp*MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho = 0.5_dp*MAX(norm_drho(ip),EPSILON(0.0_dp)*1.e4_dp)
        my_tau = 0.5_dp*MAX(EPSILON(0.0_dp)*1.e4_dp,tau(ip))
        my_laplace_rho = 0.5_dp*laplace_rho(ip)

        t1 = pi ** (0.1e1_dp / 0.3e1_dp)
        t2 = t1 ** 2
        t3 = my_rho ** (0.1e1_dp / 0.3e1_dp)
        t4 = t3 ** 2
        t5 = t4 * my_rho
        t8 = my_ndrho ** 2
        t9 = 0.1e1_dp / my_rho
        ! *** CP2K defines tau in a different way as compared to Becke !!!
        t15 = my_laplace_rho / 0.6e1_dp - gamma * (2.0_dp*my_tau - t8 * t9 / 0.4e1_dp) / 0.3e1_dp
        t16 = 0.1e1_dp / t15
        yval = 0.2e1_dp / 0.3e1_dp * t2 * t5 * t16
        IF( R == 0.0_dp ) THEN
          IF( yval <= 0.0_dp) THEN
            e_old = e_0(ip)
            CALL x_br_lsd_y_lte_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                  e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                  sx, gamma, grad_deriv)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_diff = e_0(ip) - e_old
            e_0(ip) = e_0(ip) + e_diff
          ELSE
            e_old = e_0(ip)
            CALL x_br_lsd_y_gt_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                 e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                 sx, gamma, grad_deriv)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_diff = e_0(ip) - e_old
            e_0(ip) = e_0(ip) + e_diff
          END IF
        ELSE
          IF( yval <= 0.0_dp) THEN
            e_old = e_0(ip)
            CALL x_br_lsd_y_lte_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                         e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                         sx, R, gamma, grad_deriv)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_diff = e_0(ip) - e_old
            e_0(ip) = e_0(ip) + e_diff
          ELSE
            e_old = e_0(ip)
            CALL x_br_lsd_y_gt_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                        e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                        sx, R, gamma, grad_deriv)
            ! VERY UGLY HACK e_0 has to multiplied by the factor 2
            e_diff = e_0(ip) - e_old
            e_0(ip) = e_0(ip) + e_diff
          END IF
        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE xbecke_roussel_lda_calc

! *****************************************************************************
!> \brief evaluates the Becke Roussel exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param br_params parameters for the becke roussel functional
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xbecke_roussel_lsd_eval(rho_set,deriv_set,grad_deriv,br_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: br_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xbecke_roussel_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: gamma, R, sx
    REAL(kind=dp)                            :: epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rhoa, &
      e_laplace_rhob, e_ndrhoa, e_ndrhob, e_rhoa, e_rhob, e_tau_a, e_tau_b, &
      laplace_rhoa, laplace_rhob, norm_drhoa, norm_drhob, rhoa, rhob, tau_a, &
      tau_b
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,failure)
    CALL xc_rho_set_get(rho_set, rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
         norm_drhob=norm_drhob, tau_a=tau_a, tau_b=tau_b, laplace_rhoa=laplace_rhoa,&
         laplace_rhob=laplace_rhob,local_bounds=bo,&
         rho_cutoff=epsilon_rho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)
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
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
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

      deriv => xc_dset_get_derivative(deriv_set,"(tau_a)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_tau_a)
      deriv => xc_dset_get_derivative(deriv_set,"(tau_b)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_tau_b)

      deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa)
      deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)",&
           allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob)
    END IF
    IF (grad_deriv>1.OR.grad_deriv<-1) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
           message="derivatives bigger than 1 not implemented", &
           error_level=cp_failure_level)
    END IF

    CALL section_vals_val_get(br_params,"scale_x",r_val=sx)
    CALL section_vals_val_get(br_params,"CUTOFF_RADIUS",r_val=R)
    CALL section_vals_val_get(br_params,"GAMMA",r_val=gamma)

    !$omp parallel default (none) &
    !$omp          shared(rhoa, norm_drhoa, laplace_rhoa, tau_a, e_0) &
    !$omp          shared(e_rhoa, e_ndrhoa, e_tau_a, e_laplace_rhoa) &
    !$omp          shared(grad_deriv, npoints, epsilon_rho) &
    !$omp          shared(sx, r, gamma) &
    !$omp          shared(rhob, norm_drhob, laplace_rhob, tau_b, e_rhob) &
    !$omp          shared(e_ndrhob, e_tau_b, e_laplace_rhob)

    CALL xbecke_roussel_lsd_calc(rho=rhoa, norm_drho=norm_drhoa,&
        laplace_rho=laplace_rhoa,tau=tau_a,e_0=e_0,e_rho=e_rhoa,e_ndrho=e_ndrhoa,&
        e_tau=e_tau_a,e_laplace_rho=e_laplace_rhoa,grad_deriv=grad_deriv, &
        npoints=npoints,epsilon_rho=epsilon_rho,&
        sx=sx,R=R,gamma=gamma)

    CALL xbecke_roussel_lsd_calc(rho=rhob, norm_drho=norm_drhob,&
        laplace_rho=laplace_rhob,tau=tau_b,e_0=e_0,e_rho=e_rhob,e_ndrho=e_ndrhob,&
        e_tau=e_tau_b,e_laplace_rho=e_laplace_rhob,grad_deriv=grad_deriv, &
        npoints=npoints,epsilon_rho=epsilon_rho,&
        sx=sx,R=R,gamma=gamma)

    !$omp end parallel

    IF (cp_debug) THEN
       DEALLOCATE(dummy,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
    ELSE
       NULLIFY(dummy)
    END IF
    CALL timestop(handle)
  END SUBROUTINE xbecke_roussel_lsd_eval

! *****************************************************************************
!> \brief Precalculates which branch of the code has to be taken
!>        There are two main branches, one for a truncated potential and a cutoff
!>        radius, the other for the full coulomb interaction
!> \param rho grid values
!> \param norm_drho grid values
!> \param laplace_rho grid values
!> \param tau grid values
!> \param e_0 energies and derivatives
!> \param e_rho energies and derivatives
!> \param e_ndrho energies and derivatives
!> \param e_tau energies and derivatives
!> \param e_laplace_rho energies and derivatives
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints size of the grids
!> \param epsilon_rho cutoffs
!> \param sx scales the exchange potential and energies
!> \param R cutoff Radius for truncated case
!> \param gamma paramter from original publication, usually set to 1
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xbecke_roussel_lsd_calc(rho,norm_drho,laplace_rho,tau,e_0,e_rho,&
                                     e_ndrho,e_tau,e_laplace_rho,grad_deriv,npoints,&
                                     epsilon_rho,sx,R,gamma)

    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: e_laplace_rho, e_tau, &
                                                e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: tau, laplace_rho, norm_drho, &
                                                rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, R, gamma

    CHARACTER(len=*), PARAMETER :: routineN = 'xbecke_roussel_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: my_laplace_rho, my_ndrho, &
                                                my_rho, my_tau, t1, t15, t16, &
                                                t2, t3, t4, t5, t8, t9, yval

! Precalculate y, in order to chose the correct branch afterwards

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
        ! *** CP2K defines tau in a different way as compared to Becke !!!
        t15 = my_laplace_rho / 0.6e1_dp - gamma * (2.0_dp*my_tau - t8 * t9 / 0.4e1_dp) / 0.3e1_dp
        t16 = 0.1e1_dp / t15
        yval = 0.2e1_dp / 0.3e1_dp * t2 * t5 * t16
        IF( R == 0.0_dp ) THEN
          IF( yval <= 0.0_dp) THEN
            CALL x_br_lsd_y_lte_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                  e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                  sx, gamma, grad_deriv)
          ELSE
            CALL x_br_lsd_y_gt_0(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                 e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                 sx, gamma, grad_deriv)
          END IF
        ELSE
          IF( yval <= 0.0_dp) THEN
            CALL x_br_lsd_y_lte_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                         e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                         sx, R, gamma, grad_deriv)
          ELSE
            CALL x_br_lsd_y_gt_0_cutoff(my_rho, my_ndrho, my_tau, my_laplace_rho, e_0(ip), &
                                        e_rho(ip), e_ndrho(ip), e_tau(ip), e_laplace_rho(ip),&
                                        sx, R, gamma, grad_deriv)
          END IF
        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE xbecke_roussel_lsd_calc

! *****************************************************************************
!> \brief full range evaluation for y <= 0
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_lte_0(rho, ndrho, tau, laplace_rho, e_0, &
                              e_rho, e_ndrho, e_tau, e_laplace_rho,&
                              sx, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'x_br_lsd_y_lte_0', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t100, t101, t102, t108, t111, t113, t114, t117, &
      t118, t120, t121, t122, t129, t130, t134, t135, t138, t142, t143, t146, &
      t147, t150, t152, t153, t157, t158, t16, t161, t164, t165, t166, t169, &
      t17, t170, t172, t173, t19, t199, t2, t20, t202, t207, t208, t209, t21, &
      t217, t218, t22, t220, t227, t229, t23, t234, t246, t249, t252, t255, &
      t259, t26, t263, t267, t27, t271, t274, t275, t276, t28, t29, t293, &
      t295, t3, t304, t307, t308, t32, t320, t33, t34, t340, t341, t342, &
      t344, t346, t349, t35, t350, t353, t354, t357, t358, t36, t361, t362, &
      t365, t366, t367, t37, t379, t38, t381
    REAL(KIND=dp) :: t387, t39, t4, t401, t42, t422, t43, t434, t435, t436, &
      t44, t448, t45, t450, t47, t471, t48, t5, t51, t52, t53, t54, t55, t56, &
      t57, t61, t62, t63, t64, t66, t67, t70, t71, t72, t75, t78, t81, t84, &
      t87, t88, t89, t9, t91, t92, t93, t94, t95, t96, t97, yval

    IF( grad_deriv >= 0 ) THEN
      t1 = pi ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = rho ** (0.1e1_dp / 0.3e1_dp)
      t4 = t3 ** 2
      t5 = t4 * rho
      t9 = ndrho ** 2
      t16 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t9 / rho / 0.4e1_dp) / 0.3e1_dp
      t17 = 0.1e1_dp / t16
      yval = 0.2e1_dp / 0.3e1_dp * t2 * t5 * t17
      t19 = t3 * rho
      t20 = 0.3141592654e1_dp ** (0.1e1_dp / 0.3e1_dp)
      t21 = t19 * t20
      t22 = br_a1 * t2
      t23 = t5 * t17
      t26 = 0.2e1_dp / 0.3e1_dp * t22 * t23 + br_a2
      t27 = ATAN(t26)
      t28 = -t27 + br_a3
      t29 = br_c1 * t2
      t32 = t1 * pi
      t33 = br_c2 * t32
      t34 = rho ** 2
      t35 = t34 * rho
      t36 = t3 * t35
      t37 = t16 ** 2
      t38 = 0.1e1_dp / t37
      t39 = t36 * t38
      t42 = pi ** 2
      t43 = br_c3 * t42
      t44 = t34 ** 2
      t45 = t44 * rho
      t47 = 0.1e1_dp / t37 / t16
      t48 = t45 * t47
      t51 = t2 * t42
      t52 = br_c4 * t51
      t53 = t44 * t34
      t54 = t4 * t53
      t55 = t37 ** 2
      t56 = 0.1e1_dp / t55
      t57 = t54 * t56
      t61 = t1 * t42 * pi
      t62 = br_c5 * t61
      t63 = t44 ** 2
      t64 = t3 * t63
      t66 = 0.1e1_dp / t55 / t16
      t67 = t64 * t66
      t70 = br_c0 + 0.2e1_dp / 0.3e1_dp * t29 * t23 + 0.4e1_dp / 0.9e1_dp * t33 * t39  &
      + 0.8e1_dp / 0.27e2_dp * t43 * t48 + 0.16e2_dp / 0.81e2_dp * t52 * t57 + 0.32e2_dp  &
      / 0.243e3_dp * t62 * t67
      t71 = t28 * t70
      t72 = br_b1 * t2
      t75 = br_b2 * t32
      t78 = br_b3 * t42
      t81 = br_b4 * t51
      t84 = br_b5 * t61
      t87 = br_b0 + 0.2e1_dp / 0.3e1_dp * t72 * t23 + 0.4e1_dp / 0.9e1_dp * t75 * t39  &
      + 0.8e1_dp / 0.27e2_dp * t78 * t48 + 0.16e2_dp / 0.81e2_dp * t81 * t57 + 0.32e2_dp  &
      / 0.243e3_dp * t84 * t67
      t88 = 0.1e1_dp / t87
      t89 = t71 * t88
      t91 = EXP(t89 / 0.3e1_dp)
      t92 = t21 * t91
      t93 = 0.1e1_dp / t28
      t94 = 0.1e1_dp / t70
      t95 = t93 * t94
      t96 = EXP(-t89)
      t97 = t88 * t96
      t100 = 0.1e1_dp - t96 - t71 * t97 / 0.2e1_dp
      t101 = t87 * t100
      t102 = t95 * t101
      e_0 = e_0 + (-t92 * t102)* sx
    END IF
    IF(grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
      t108 = t4 * t17
      t111 = 0.1e1_dp / t3
      t113 = t38 * gamma
      t114 = t113 * t9
      t117 = 0.10e2_dp / 0.9e1_dp * t22 * t108 + t22 * t111 * t114 / 0.18e2_dp
      t118 = t26 ** 2
      t120 = 0.1e1_dp / (0.1e1_dp + t118)
      t121 = t117 * t120
      t122 = t70 * t88
      t129 = t3 * t34
      t130 = t129 * t38
      t134 = t47 * gamma
      t135 = t134 * t9
      t138 = t44 * t47
      t142 = t56 * gamma
      t143 = t142 * t9
      t146 = t4 * t45
      t147 = t146 * t56
      t150 = t4 * t44
      t152 = t66 * gamma
      t153 = t152 * t9
      t157 = t3 * t44 * t35
      t158 = t157 * t66
      t161 = t3 * t53
      t164 = 0.1e1_dp / t55 / t37
      t165 = t164 * gamma
      t166 = t165 * t9
      t169 = 0.10e2_dp / 0.9e1_dp * t29 * t108 + t29 * t111 * t114 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t33 * t130 + 0.2e1_dp / 0.27e2_dp * t33 * t19 * t135 +  &
      0.40e2_dp / 0.27e2_dp * t43 * t138 + 0.2e1_dp / 0.27e2_dp * t43 * t35 * t143 +  &
      0.320e3_dp / 0.243e3_dp * t52 * t147 + 0.16e2_dp / 0.243e3_dp * t52 * t150 *  &
      t153 + 0.800e3_dp / 0.729e3_dp * t62 * t158 + 0.40e2_dp / 0.729e3_dp * t62 * t161  &
      * t166
      t170 = t28 * t169
      t172 = t87 ** 2
      t173 = 0.1e1_dp / t172
      t199 = 0.10e2_dp / 0.9e1_dp * t72 * t108 + t72 * t111 * t114 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t75 * t130 + 0.2e1_dp / 0.27e2_dp * t75 * t19 * t135 +  &
      0.40e2_dp / 0.27e2_dp * t78 * t138 + 0.2e1_dp / 0.27e2_dp * t78 * t35 * t143 +  &
      0.320e3_dp / 0.243e3_dp * t81 * t147 + 0.16e2_dp / 0.243e3_dp * t81 * t150 *  &
      t153 + 0.800e3_dp / 0.729e3_dp * t84 * t158 + 0.40e2_dp / 0.729e3_dp * t84 * t161  &
      * t166
      t202 = -t121 * t122 + t170 * t88 - t71 * t173 * t199
      t207 = t28 ** 2
      t208 = 0.1e1_dp / t207
      t209 = t91 * t208
      t217 = t21 * t91 * t93
      t218 = t70 ** 2
      t220 = 0.1e1_dp / t218 * t87
      t227 = -t202
      t229 = t122 * t96
      t234 = t173 * t96
      e_rho = e_rho + (-0.4e1_dp / 0.3e1_dp * t3 * t20 * t91 * t102 - t21 * t202 * t91 *  &
      t102 / 0.3e1_dp - t21 * t209 * t94 * t87 * t100 * t117 * t120 + t217  &
      * t220 * t100 * t169 - t92 * t95 * t199 * t100 - t92 * t95 * t87 *  &
      (-t227 * t96 + t121 * t229 / 0.2e1_dp - t170 * t97 / 0.2e1_dp + t71 * t234  &
      * t199 / 0.2e1_dp - t71 * t88 * t227 * t96 / 0.2e1_dp))*sx
      t246 = t4 * t38
      t249 = t120 * t70
      t252 = t22 * t246 * gamma * ndrho * t249 * t88
      t255 = t113 * ndrho
      t259 = t134 * ndrho
      t263 = t142 * ndrho
      t267 = t152 * ndrho
      t271 = t165 * ndrho
      t274 = -t29 * t4 * t255 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t33 * t129 * t259  &
      - 0.4e1_dp / 0.27e2_dp * t43 * t44 * t263 - 0.32e2_dp / 0.243e3_dp * t52 * t146  &
      * t267 - 0.80e2_dp / 0.729e3_dp * t62 * t157 * t271
      t275 = t28 * t274
      t276 = t275 * t88
      t293 = -t72 * t4 * t255 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t75 * t129 * t259  &
      - 0.4e1_dp / 0.27e2_dp * t78 * t44 * t263 - 0.32e2_dp / 0.243e3_dp * t81 * t146  &
      * t267 - 0.80e2_dp / 0.729e3_dp * t84 * t157 * t271
      t295 = t71 * t173 * t293
      t304 = t208 * t94 * t87
      t307 = t100 * br_a1 * t2
      t308 = ndrho * t120
      t320 = -t252 / 0.9e1_dp - t276 + t295
      e_ndrho = e_ndrho + (-t21 * (t252 / 0.27e2_dp + t276 / 0.3e1_dp - t295 / 0.3e1_dp) * t91  &
      * t102 + t34 * t20 * t91 * t304 * t307 * t113 * t308 / 0.9e1_dp + t217  &
      * t220 * t100 * t274 - t92 * t95 * t293 * t100 - t92 * t95 * t87  &
      * (-t320 * t96 - t22 * t246 * gamma * t308 * t229 / 0.18e2_dp - t275  &
      * t97 / 0.2e1_dp + t71 * t234 * t293 / 0.2e1_dp - t71 * t88 * t320 * t96  &
      / 0.2e1_dp)) * sx
      t340 = t5 * t38
      t341 = t22 * t340
      t342 = gamma * t120
      t344 = t341 * t342 * t122
      t346 = t340 * gamma
      t349 = t36 * t47
      t350 = t349 * gamma
      t353 = t45 * t56
      t354 = t353 * gamma
      t357 = t54 * t66
      t358 = t357 * gamma
      t361 = t64 * t164
      t362 = t361 * gamma
      t365 = 0.4e1_dp / 0.9e1_dp * t29 * t346 + 0.16e2_dp / 0.27e2_dp * t33 * t350 +  &
      0.16e2_dp / 0.27e2_dp * t43 * t354 + 0.128e3_dp / 0.243e3_dp * t52 * t358 + 0.320e3_dp  &
      / 0.729e3_dp * t62 * t362
      t366 = t28 * t365
      t367 = t366 * t88
      t379 = 0.4e1_dp / 0.9e1_dp * t72 * t346 + 0.16e2_dp / 0.27e2_dp * t75 * t350 +  &
      0.16e2_dp / 0.27e2_dp * t78 * t354 + 0.128e3_dp / 0.243e3_dp * t81 * t358 + 0.320e3_dp  &
      / 0.729e3_dp * t84 * t362
      t381 = t71 * t173 * t379
      t387 = t35 * t20
      t401 = 0.4e1_dp / 0.9e1_dp * t344 - t367 + t381
      e_tau = e_tau + (-t21 * (-0.4e1_dp / 0.27e2_dp * t344 + t367 / 0.3e1_dp - t381 / 0.3e1_dp)  &
      * t91 * t102 - 0.4e1_dp / 0.9e1_dp * t387 * t91 * t304 * t307 * t113 *  &
      t120 + t217 * t220 * t100 * t365 - t92 * t95 * t379 * t100 - t92  &
      * t95 * t87 * (-t401 * t96 + 0.2e1_dp / 0.9e1_dp * t341 * t342 * t229 -  &
      t366 * t97 / 0.2e1_dp + t71 * t234 * t379 / 0.2e1_dp - t71 * t88 * t401  &
      * t96 / 0.2e1_dp)) * sx
      t422 = t22 * t5 * t38 * t120 * t122
      t434 = -t29 * t340 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t33 * t349 - 0.4e1_dp /  &
      0.27e2_dp * t43 * t353 - 0.32e2_dp / 0.243e3_dp * t52 * t357 - 0.80e2_dp / 0.729e3_dp  &
      * t62 * t361
      t435 = t28 * t434
      t436 = t435 * t88
      t448 = -t72 * t340 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t75 * t349 - 0.4e1_dp /  &
      0.27e2_dp * t78 * t353 - 0.32e2_dp / 0.243e3_dp * t81 * t357 - 0.80e2_dp / 0.729e3_dp  &
      * t84 * t361
      t450 = t71 * t173 * t448
      t471 = -t422 / 0.9e1_dp - t436 + t450
      e_laplace_rho = e_laplace_rho + (-t21 * (t422 / 0.27e2_dp + t436 / 0.3e1_dp - t450 / 0.3e1_dp) * t91 *  &
      t102 + t387 * t209 * t94 * t101 * br_a1 * t2 * t38 * t120 / 0.9e1_dp  &
      + t217 * t220 * t100 * t434 - t92 * t95 * t448 * t100 - t92 * t95  &
      * t87 * (-t471 * t96 - t341 * t249 * t97 / 0.18e2_dp - t435 * t97 / 0.2e1_dp  &
      + t71 * t234 * t448 / 0.2e1_dp - t71 * t88 * t471 * t96 / 0.2e1_dp)) *sx
    END IF
  END SUBROUTINE x_br_lsd_y_lte_0

! *****************************************************************************
!> \brief Full range evaluation for y > 0
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_gt_0(rho, ndrho, tau, laplace_rho, e_0, &
                             e_rho, e_ndrho, e_tau, e_laplace_rho,&
                             sx, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'x_br_lsd_y_gt_0', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t102, t103, t104, t106, t107, t108, t109, t110, &
      t111, t112, t115, t117, t124, t131, t134, t137, t138, t142, t143, t154, &
      t157, t158, t159, t16, t160, t162, t165, t167, t168, t17, t176, t180, &
      t181, t184, t185, t188, t19, t190, t191, t195, t196, t199, t2, t20, &
      t202, t203, t204, t207, t208, t21, t210, t211, t22, t23, t237, t24, &
      t240, t245, t248, t249, t25, t255, t256, t258, t26, t265, t267, t272, &
      t285, t288, t289, t29, t297, t298, t3, t30, t301, t305, t309, t31, &
      t313, t317, t32, t320, t321, t33, t338, t34, t341, t35, t356, t36, t37, &
      t376, t377, t382, t383, t387, t388, t391
    REAL(KIND=dp) :: t392, t395, t396, t399, t4, t400, t403, t404, t41, t416, &
      t419, t42, t43, t434, t458, t459, t47, t471, t472, t48, t484, t487, &
      t49, t5, t50, t502, t51, t54, t57, t58, t59, t6, t60, t62, t63, t66, &
      t67, t68, t69, t70, t71, t72, t76, t77, t78, t79, t81, t82, t85, t86, &
      t87, t9, t90, t93, t96, t99, yval

    IF( grad_deriv >= 0 ) THEN
      t1 = pi ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = rho ** (0.1e1_dp / 0.3e1_dp)
      t4 = t3 ** 2
      t5 = t4 * rho
      t6 = t2 * t5
      t9 = ndrho ** 2
      t16 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t9 / rho / 0.4e1_dp) / 0.3e1_dp
      t17 = 0.1e1_dp / t16
      yval = 0.2e1_dp / 0.3e1_dp * t6 * t17
      t19 = t3 * rho
      t20 = 0.3141592654e1_dp ** (0.1e1_dp / 0.3e1_dp)
      t21 = t19 * t20
      t22 = 0.1e1_dp / br_BB
      t23 = 0.1e1_dp / t2
      t24 = t22 * t23
      t25 = 0.1e1_dp / t5
      t26 = t25 * t16
      t29 = br_BB ** 2
      t30 = t1 * pi
      t31 = t29 * t30
      t32 = rho ** 2
      t33 = t32 * rho
      t34 = t3 * t33
      t35 = t16 ** 2
      t36 = 0.1e1_dp / t35
      t37 = t34 * t36
      t41 = SQRT(0.10e1_dp + 0.4e1_dp / 0.9e1_dp * t31 * t37)
      t42 = t41 * t22
      t43 = t23 * t25
      t47 = 0.1500000000000000e1_dp * t24 * t26 + 0.3e1_dp / 0.2e1_dp * t42 * t43  &
      * t16
      t48 = LOG(t47)
      t49 = t48 + 0.2e1_dp
      t50 = br_d1 * t2
      t51 = t5 * t17
      t54 = br_d2 * t30
      t57 = pi ** 2
      t58 = br_d3 * t57
      t59 = t32 ** 2
      t60 = t59 * rho
      t62 = 0.1e1_dp / t35 / t16
      t63 = t60 * t62
      t66 = t2 * t57
      t67 = br_d4 * t66
      t68 = t59 * t32
      t69 = t4 * t68
      t70 = t35 ** 2
      t71 = 0.1e1_dp / t70
      t72 = t69 * t71
      t76 = t1 * t57 * pi
      t77 = br_d5 * t76
      t78 = t59 ** 2
      t79 = t3 * t78
      t81 = 0.1e1_dp / t70 / t16
      t82 = t79 * t81
      t85 = br_d0 + 0.2e1_dp / 0.3e1_dp * t50 * t51 + 0.4e1_dp / 0.9e1_dp * t54 * t37  &
      + 0.8e1_dp / 0.27e2_dp * t58 * t63 + 0.16e2_dp / 0.81e2_dp * t67 * t72 + 0.32e2_dp  &
      / 0.243e3_dp * t77 * t82
      t86 = t49 * t85
      t87 = br_e1 * t2
      t90 = br_e2 * t30
      t93 = br_e3 * t57
      t96 = br_e4 * t66
      t99 = br_e5 * t76
      t102 = br_e0 + 0.2e1_dp / 0.3e1_dp * t87 * t51 + 0.4e1_dp / 0.9e1_dp * t90 * t37  &
      + 0.8e1_dp / 0.27e2_dp * t93 * t63 + 0.16e2_dp / 0.81e2_dp * t96 * t72 + 0.32e2_dp  &
      / 0.243e3_dp * t99 * t82
      t103 = 0.1e1_dp / t102
      t104 = t86 * t103
      t106 = EXP(t104 / 0.3e1_dp)
      t107 = t21 * t106
      t108 = 0.1e1_dp / t49
      t109 = 0.1e1_dp / t85
      t110 = t108 * t109
      t111 = EXP(-t104)
      t112 = t103 * t111
      t115 = 0.1e1_dp - t111 - t86 * t112 / 0.2e1_dp
      t117 = t110 * t102 * t115
      e_0 = e_0 +(-t107 * t117) * sx
    END IF
    IF(grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
      t124 = 0.1e1_dp / t4 / t32
      t131 = 0.1e1_dp / t4 / t33 * gamma * t9
      t134 = 0.1e1_dp / t41
      t137 = t3 * t32
      t138 = t137 * t36
      t142 = t62 * gamma
      t143 = t142 * t9
      t154 = t42 * t23
      t157 = -0.2500000000e1_dp * t24 * t124 * t16 - 0.1250000000e0_dp * t24 *  &
      t131 + 0.3e1_dp / 0.4e1_dp * t134 * t22 * t23 * t26 * (0.40e2_dp / 0.27e2_dp *  &
      t31 * t138 + 0.2e1_dp / 0.27e2_dp * t31 * t19 * t143) - 0.5e1_dp / 0.2e1_dp *  &
      t42 * t23 * t124 * t16 - t154 * t131 / 0.8e1_dp
      t158 = 0.1e1_dp / t47
      t159 = t157 * t158
      t160 = t85 * t103
      t162 = t4 * t17
      t165 = 0.1e1_dp / t3
      t167 = t36 * gamma
      t168 = t167 * t9
      t176 = t59 * t62
      t180 = t71 * gamma
      t181 = t180 * t9
      t184 = t4 * t60
      t185 = t184 * t71
      t188 = t4 * t59
      t190 = t81 * gamma
      t191 = t190 * t9
      t195 = t3 * t59 * t33
      t196 = t195 * t81
      t199 = t3 * t68
      t202 = 0.1e1_dp / t70 / t35
      t203 = t202 * gamma
      t204 = t203 * t9
      t207 = 0.10e2_dp / 0.9e1_dp * t50 * t162 + t50 * t165 * t168 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t54 * t138 + 0.2e1_dp / 0.27e2_dp * t54 * t19 * t143 +  &
      0.40e2_dp / 0.27e2_dp * t58 * t176 + 0.2e1_dp / 0.27e2_dp * t58 * t33 * t181 +  &
      0.320e3_dp / 0.243e3_dp * t67 * t185 + 0.16e2_dp / 0.243e3_dp * t67 * t188 *  &
      t191 + 0.800e3_dp / 0.729e3_dp * t77 * t196 + 0.40e2_dp / 0.729e3_dp * t77 * t199  &
      * t204
      t208 = t49 * t207
      t210 = t102 ** 2
      t211 = 0.1e1_dp / t210
      t237 = 0.10e2_dp / 0.9e1_dp * t87 * t162 + t87 * t165 * t168 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t90 * t138 + 0.2e1_dp / 0.27e2_dp * t90 * t19 * t143 +  &
      0.40e2_dp / 0.27e2_dp * t93 * t176 + 0.2e1_dp / 0.27e2_dp * t93 * t33 * t181 +  &
      0.320e3_dp / 0.243e3_dp * t96 * t185 + 0.16e2_dp / 0.243e3_dp * t96 * t188 *  &
      t191 + 0.800e3_dp / 0.729e3_dp * t99 * t196 + 0.40e2_dp / 0.729e3_dp * t99 * t199  &
      * t204
      t240 = t159 * t160 + t208 * t103 - t86 * t211 * t237
      t245 = t49 ** 2
      t248 = t21 * t106 / t245
      t249 = t109 * t102
      t255 = t21 * t106 * t108
      t256 = t85 ** 2
      t258 = 0.1e1_dp / t256 * t102
      t265 = -t240
      t267 = t160 * t111
      t272 = t211 * t111
      e_rho = e_rho + (-0.4e1_dp / 0.3e1_dp * t3 * t20 * t106 * t117 - t21 * t240 * t106  &
      * t117 / 0.3e1_dp + t248 * t249 * t115 * t157 * t158 + t255 * t258 *  &
      t115 * t207 - t107 * t110 * t237 * t115 - t107 * t110 * t102 * (-t265  &
      * t111 - t159 * t267 / 0.2e1_dp - t208 * t112 / 0.2e1_dp + t86 * t272  &
      * t237 / 0.2e1_dp - t86 * t103 * t265 * t111 / 0.2e1_dp)) * sx
      t285 = t124 * gamma * ndrho
      t288 = t134 * br_BB
      t289 = t288 * t2
      t297 = 0.2500000000000000e0_dp * t24 * t285 - t289 * t4 * t36 * gamma *  &
      ndrho / 0.9e1_dp + t154 * t285 / 0.4e1_dp
      t298 = t297 * t158
      t301 = t167 * ndrho
      t305 = t142 * ndrho
      t309 = t180 * ndrho
      t313 = t190 * ndrho
      t317 = t203 * ndrho
      t320 = -t50 * t4 * t301 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t54 * t137 * t305  &
      - 0.4e1_dp / 0.27e2_dp * t58 * t59 * t309 - 0.32e2_dp / 0.243e3_dp * t67 * t184  &
      * t313 - 0.80e2_dp / 0.729e3_dp * t77 * t195 * t317
      t321 = t49 * t320
      t338 = -t87 * t4 * t301 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t90 * t137 * t305  &
      - 0.4e1_dp / 0.27e2_dp * t93 * t59 * t309 - 0.32e2_dp / 0.243e3_dp * t96 * t184  &
      * t313 - 0.80e2_dp / 0.729e3_dp * t99 * t195 * t317
      t341 = t298 * t160 + t321 * t103 - t86 * t211 * t338
      t356 = -t341
      e_ndrho = e_ndrho + (-t21 * t341 * t106 * t117 / 0.3e1_dp + t248 * t249 * t115 * t297  &
      * t158 + t255 * t258 * t115 * t320 - t107 * t110 * t338 * t115  &
      - t107 * t110 * t102 * (-t356 * t111 - t298 * t267 / 0.2e1_dp - t321  &
      * t112 / 0.2e1_dp + t86 * t272 * t338 / 0.2e1_dp - t86 * t103 * t356 *  &
      t111 / 0.2e1_dp)) * sx
      t376 = t5 * t36
      t377 = t376 * gamma
      t382 = -0.1000000000e1_dp * t24 * t25 * gamma + 0.4e1_dp / 0.9e1_dp * t289 * t377  &
      - t42 * t43 * gamma
      t383 = t382 * t158
      t387 = t34 * t62
      t388 = t387 * gamma
      t391 = t60 * t71
      t392 = t391 * gamma
      t395 = t69 * t81
      t396 = t395 * gamma
      t399 = t79 * t202
      t400 = t399 * gamma
      t403 = 0.4e1_dp / 0.9e1_dp * t50 * t377 + 0.16e2_dp / 0.27e2_dp * t54 * t388 +  &
      0.16e2_dp / 0.27e2_dp * t58 * t392 + 0.128e3_dp / 0.243e3_dp * t67 * t396 + 0.320e3_dp  &
      / 0.729e3_dp * t77 * t400
      t404 = t49 * t403
      t416 = 0.4e1_dp / 0.9e1_dp * t87 * t377 + 0.16e2_dp / 0.27e2_dp * t90 * t388 +  &
      0.16e2_dp / 0.27e2_dp * t93 * t392 + 0.128e3_dp / 0.243e3_dp * t96 * t396 + 0.320e3_dp  &
      / 0.729e3_dp * t99 * t400
      t419 = t383 * t160 + t404 * t103 - t86 * t211 * t416
      t434 = -t419
      e_tau = e_tau + (-t21 * t419 * t106 * t117 / 0.3e1_dp + t248 * t249 * t115 * t382  &
      * t158 + t255 * t258 * t115 * t403 - t107 * t110 * t416 * t115 -  &
      t107 * t110 * t102 * (-t434 * t111 - t383 * t267 / 0.2e1_dp - t404 *  &
      t112 / 0.2e1_dp + t86 * t272 * t416 / 0.2e1_dp - t86 * t103 * t434 * t111  &
      / 0.2e1_dp))* sx
      t458 = 0.2500000000000000e0_dp * t24 * t25 - t288 * t6 * t36 / 0.9e1_dp +  &
      t42 * t43 / 0.4e1_dp
      t459 = t458 * t158
      t471 = -t50 * t376 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t54 * t387 - 0.4e1_dp /  &
      0.27e2_dp * t58 * t391 - 0.32e2_dp / 0.243e3_dp * t67 * t395 - 0.80e2_dp / 0.729e3_dp  &
      * t77 * t399
      t472 = t49 * t471
      t484 = -t87 * t376 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t90 * t387 - 0.4e1_dp /  &
      0.27e2_dp * t93 * t391 - 0.32e2_dp / 0.243e3_dp * t96 * t395 - 0.80e2_dp / 0.729e3_dp  &
      * t99 * t399
      t487 = t459 * t160 + t472 * t103 - t86 * t211 * t484
      t502 = -t487
      e_laplace_rho = e_laplace_rho+(-t21 * t487 * t106 * t117 / 0.3e1_dp + t248 * t249 * t115 * t458  &
      * t158 + t255 * t258 * t115 * t471 - t107 * t110 * t484 * t115 -  &
      t107 * t110 * t102 * (-t502 * t111 - t459 * t267 / 0.2e1_dp - t472 *  &
      t112 / 0.2e1_dp + t86 * t272 * t484 / 0.2e1_dp - t86 * t103 * t502 * t111  &
      / 0.2e1_dp))* sx
    END IF

  END SUBROUTINE x_br_lsd_y_gt_0

! *****************************************************************************
!> \brief Truncated long range part for y <= 0
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param R ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_lte_0_cutoff(rho, ndrho, tau, laplace_rho, e_0, &
                                     e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                     sx, R, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, R, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'x_br_lsd_y_lte_0_cutoff', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: brval, t1, t101, t11, t12, t18, t2, t20, t24, t25, t26, &
      t3, t31, t33, t36, t38, t4, t41, t43, t47, t50, t54, t56, t6, t60, t62, &
      t66, t69, t7, t70, t88, t89, t96

      t1 = 8 ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = pi ** (0.1e1_dp / 0.3e1_dp)
      t4 = t3 ** 2
      t6 = rho ** (0.1e1_dp / 0.3e1_dp)
      t7 = t6 ** 2
      t11 = ndrho ** 2
      t12 = 0.1e1_dp / rho
      t18 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t11 * t12 / 0.4e1_dp) / 0.3e1_dp
      t20 = t7 * rho / t18
      t24 = ATAN(0.2e1_dp / 0.3e1_dp * br_a1 * t4 * t20 + br_a2)
      t25 = -t24 + br_a3
      t26 = t25 ** 2
      t31 = t3 * pi
      t33 = rho ** 2
      t36 = t18 ** 2
      t38 = t6 * t33 * rho / t36
      t41 = pi ** 2
      t43 = t33 ** 2
      t47 = t43 * rho / t36 / t18
      t50 = t4 * t41
      t54 = t36 ** 2
      t56 = t7 * t43 * t33 / t54
      t60 = t3 * t41 * pi
      t62 = t43 ** 2
      t66 = t6 * t62 / t54 / t18
      t69 = br_c0 + 0.2e1_dp / 0.3e1_dp * br_c1 * t4 * t20 + 0.4e1_dp / 0.9e1_dp * br_c2  &
      * t31 * t38 + 0.8e1_dp / 0.27e2_dp * br_c3 * t41 * t47 + 0.16e2_dp / 0.81e2_dp  &
      * br_c4 * t50 * t56 + 0.32e2_dp / 0.243e3_dp * br_c5 * t60 * t66
      t70 = t69 ** 2
      t88 = br_b0 + 0.2e1_dp / 0.3e1_dp * br_b1 * t4 * t20 + 0.4e1_dp / 0.9e1_dp * br_b2  &
      * t31 * t38 + 0.8e1_dp / 0.27e2_dp * br_b3 * t41 * t47 + 0.16e2_dp / 0.81e2_dp  &
      * br_b4 * t50 * t56 + 0.32e2_dp / 0.243e3_dp * br_b5 * t60 * t66
      t89 = t88 ** 2
      t96 = EXP(-t25 * t69 / t88)
      t101 = (t26 * t25 * t70 * t69 / t89 / t88 * t96 / 0.3141592654e1_dp *  &
      t12) ** (0.1e1_dp / 0.3e1_dp)
      brval = REAL(t2,KIND=dp) * t101 / 0.8e1_dp

    IF(R>brval) THEN
      CALL x_br_lsd_y_lte_0_cutoff_R_gt_b(rho, ndrho, tau, laplace_rho, e_0, &
                                           e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                           sx, R, gamma, grad_deriv)
    ELSE
      CALL x_br_lsd_y_lte_0_cutoff_R_lte_b(rho, ndrho, tau, laplace_rho, e_0, &
                                          e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                          sx, R, gamma, grad_deriv)
    END IF

  END SUBROUTINE x_br_lsd_y_lte_0_cutoff

! *****************************************************************************
!> \brief Truncated long range part for y <= 0 and R > b
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param R ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_lte_0_cutoff_R_gt_b(rho, ndrho, tau, laplace_rho, e_0, &
                                             e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                             sx, R, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, R, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'x_br_lsd_y_lte_0_cutoff_R_gt_b', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t10, t100, t101, t102, t103, t104, t106, t108, t109, &
      t110, t114, t117, t119, t120, t123, t124, t129, t132, t134, t135, t138, &
      t139, t141, t142, t143, t149, t150, t153, t155, t156, t159, t16, t163, &
      t164, t167, t168, t17, t171, t173, t174, t178, t179, t18, t182, t185, &
      t186, t187, t190, t192, t193, t2, t21, t219, t22, t221, t223, t224, &
      t225, t227, t228, t23, t231, t232, t233, t235, t24, t240, t245, t247, &
      t259, t261, t263, t265, t268, t269, t27, t270, t272, t274, t275, t28, &
      t283, t288, t29, t291, t3, t30, t301, t305, t31, t312, t314, t315, &
      t317, t319, t32, t323, t327, t33
    REAL(KIND=dp) :: t331, t335, t338, t34, t340, t356, t358, t359, t361, &
      t363, t364, t366, t367, t37, t38, t388, t39, t390, t392, t394, t397, &
      t399, t4, t40, t403, t405, t409, t410, t414, t42, t423, t426, t43, &
      t434, t443, t450, t455, t456, t459, t46, t460, t463, t464, t467, t468, &
      t47, t471, t472, t475, t477, t48, t488, t49, t490, t493, t494, t496, &
      t497, t499, t5, t50, t51, t516, t518, t52, t520, t522, t525, t526, &
      t527, t530, t532, t545, t548, t56, t563, t57, t572, t574, t58, t585, &
      t587, t59, t598, t6, t600, t605, t606, t608, t609, t61, t62, t628, &
      t630, t632, t634, t639, t641, t645, t65, t655
    REAL(KIND=dp) :: t66, t67, t672, t70, t73, t76, t79, t82, t83, t84, t85, &
      t86, t87, t88, t89, t9, t90, t91, t93, t94, t95, t96, t97, t99

    IF( grad_deriv >= 0 ) THEN
      t1 = pi ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = br_a1 * t2
      t4 = rho ** (0.1e1_dp / 0.3e1_dp)
      t5 = t4 ** 2
      t6 = t5 * rho
      t9 = ndrho ** 2
      t10 = 0.1e1_dp / rho
      t16 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t9 * t10 / 0.4e1_dp) / 0.3e1_dp
      t17 = 0.1e1_dp / t16
      t18 = t6 * t17
      t21 = 0.2e1_dp / 0.3e1_dp * t3 * t18 + br_a2
      t22 = ATAN(t21)
      t23 = -t22 + br_a3
      t24 = br_c1 * t2
      t27 = t1 * pi
      t28 = br_c2 * t27
      t29 = rho ** 2
      t30 = t29 * rho
      t31 = t4 * t30
      t32 = t16 ** 2
      t33 = 0.1e1_dp / t32
      t34 = t31 * t33
      t37 = pi ** 2
      t38 = br_c3 * t37
      t39 = t29 ** 2
      t40 = t39 * rho
      t42 = 0.1e1_dp / t32 / t16
      t43 = t40 * t42
      t46 = t2 * t37
      t47 = br_c4 * t46
      t48 = t39 * t29
      t49 = t5 * t48
      t50 = t32 ** 2
      t51 = 0.1e1_dp / t50
      t52 = t49 * t51
      t56 = t1 * t37 * pi
      t57 = br_c5 * t56
      t58 = t39 ** 2
      t59 = t4 * t58
      t61 = 0.1e1_dp / t50 / t16
      t62 = t59 * t61
      t65 = br_c0 + 0.2e1_dp / 0.3e1_dp * t24 * t18 + 0.4e1_dp / 0.9e1_dp * t28 * t34  &
      + 0.8e1_dp / 0.27e2_dp * t38 * t43 + 0.16e2_dp / 0.81e2_dp * t47 * t52 + 0.32e2_dp  &
      / 0.243e3_dp * t57 * t62
      t66 = t23 * t65
      t67 = br_b1 * t2
      t70 = br_b2 * t27
      t73 = br_b3 * t37
      t76 = br_b4 * t46
      t79 = br_b5 * t56
      t82 = br_b0 + 0.2e1_dp / 0.3e1_dp * t67 * t18 + 0.4e1_dp / 0.9e1_dp * t70 * t34  &
      + 0.8e1_dp / 0.27e2_dp * t73 * t43 + 0.16e2_dp / 0.81e2_dp * t76 * t52 + 0.32e2_dp  &
      / 0.243e3_dp * t79 * t62
      t83 = 0.1e1_dp / t82
      t84 = t66 * t83
      t85 = 8 ** (0.1e1_dp / 0.3e1_dp)
      t86 = t23 ** 2
      t87 = t86 * t23
      t88 = t65 ** 2
      t89 = t88 * t65
      t90 = t87 * t89
      t91 = t82 ** 2
      t93 = 0.1e1_dp / t91 / t82
      t94 = t90 * t93
      t95 = EXP(-t84)
      t96 = 0.1e1_dp / 0.3141592654e1_dp
      t97 = t95 * t96
      t99 = t94 * t97 * t10
      t100 = t99 ** (0.1e1_dp / 0.3e1_dp)
      t101 = 0.1e1_dp / t100
      t102 = REAL(t85,KIND=dp) * t101
      t103 = t102 * R
      t104 = t84 * t103
      t106 = EXP(t84 - t104)
      t108 = t106 * t23
      t109 = t65 * t83
      t110 = t108 * t109
      t114 = t83 * REAL(t85,KIND=dp) * t101 * R
      t117 = EXP(-t84 - t104)
      t119 = t117 * t23
      t120 = t119 * t109
      t123 = -0.2e1_dp * t106 + t110 - t108 * t65 * t114 + 0.2e1_dp * t117 + t120  &
      + t119 * t65 * t114
      t124 = rho * t123
      e_0 = e_0 + (t124 * t102 / 0.8e1_dp) * sx
    END IF
    IF(grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
      t129 = t5 * t17
      t132 = 0.1e1_dp / t4
      t134 = t33 * gamma
      t135 = t134 * t9
      t138 = 0.10e2_dp / 0.9e1_dp * t3 * t129 + t3 * t132 * t135 / 0.18e2_dp
      t139 = t21 ** 2
      t141 = 0.1e1_dp / (0.1e1_dp + t139)
      t142 = t138 * t141
      t143 = t142 * t109
      t149 = t4 * t29
      t150 = t149 * t33
      t153 = t4 * rho
      t155 = t42 * gamma
      t156 = t155 * t9
      t159 = t39 * t42
      t163 = t51 * gamma
      t164 = t163 * t9
      t167 = t5 * t40
      t168 = t167 * t51
      t171 = t5 * t39
      t173 = t61 * gamma
      t174 = t173 * t9
      t178 = t4 * t39 * t30
      t179 = t178 * t61
      t182 = t4 * t48
      t185 = 0.1e1_dp / t50 / t32
      t186 = t185 * gamma
      t187 = t186 * t9
      t190 = 0.10e2_dp / 0.9e1_dp * t24 * t129 + t24 * t132 * t135 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t28 * t150 + 0.2e1_dp / 0.27e2_dp * t28 * t153 * t156 +  &
      0.40e2_dp / 0.27e2_dp * t38 * t159 + 0.2e1_dp / 0.27e2_dp * t38 * t30 * t164  &
      + 0.320e3_dp / 0.243e3_dp * t47 * t168 + 0.16e2_dp / 0.243e3_dp * t47 * t171 *  &
      t174 + 0.800e3_dp / 0.729e3_dp * t57 * t179 + 0.40e2_dp / 0.729e3_dp * t57 *  &
      t182 * t187
      t192 = t23 * t190 * t83
      t193 = 0.1e1_dp / t91
      t219 = 0.10e2_dp / 0.9e1_dp * t67 * t129 + t67 * t132 * t135 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t70 * t150 + 0.2e1_dp / 0.27e2_dp * t70 * t153 * t156 +  &
      0.40e2_dp / 0.27e2_dp * t73 * t159 + 0.2e1_dp / 0.27e2_dp * t73 * t30 * t164  &
      + 0.320e3_dp / 0.243e3_dp * t76 * t168 + 0.16e2_dp / 0.243e3_dp * t76 * t171 *  &
      t174 + 0.800e3_dp / 0.729e3_dp * t79 * t179 + 0.40e2_dp / 0.729e3_dp * t79 *  &
      t182 * t187
      t221 = t66 * t193 * t219
      t223 = t142 * t65 * t114
      t224 = t192 * t103
      t225 = t66 * t193
      t227 = t102 * R * t219
      t228 = t225 * t227
      t231 = REAL(t85,KIND=dp) / t100 / t99
      t232 = t86 * t89
      t233 = t93 * t95
      t235 = t96 * t10
      t240 = t87 * t88 * t93
      t245 = t91 ** 2
      t247 = t90 / t245
      t259 = -0.3e1_dp * t232 * t233 * t235 * t142 + 0.3e1_dp * t240 * t97 * t10  &
      * t190 - 0.3e1_dp * t247 * t97 * t10 * t219 + t94 * (t143 - t192 +  &
      t221) * t95 * t235 - t94 * t97 / t29
      t261 = t231 * R * t259
      t263 = t84 * t261 / 0.3e1_dp
      t265 = (-t143 + t192 - t221 + t223 - t224 + t228 + t263) * t106
      t268 = t106 * t138
      t269 = t141 * t65
      t270 = t269 * t83
      t272 = t190 * t83
      t274 = t65 * t193
      t275 = t274 * t219
      t283 = t108 * t274
      t288 = (t143 - t192 + t221 + t223 - t224 + t228 + t263) * t117
      t291 = t117 * t138
      t301 = t119 * t274
      t305 = -0.2e1_dp * t265 + t265 * t84 - t268 * t270 + t108 * t272 - t108  &
      * t275 - t265 * t66 * t114 + t268 * t269 * t114 - t108 * t190 *  &
      t114 + t283 * t227 + t110 * t261 / 0.3e1_dp + 0.2e1_dp * t288 + t288 * t84  &
      - t291 * t270 + t119 * t272 - t119 * t275 + t288 * t66 * t114 -  &
      t291 * t269 * t114 + t119 * t190 * t114 - t301 * t227 - t120 * t261  &
      / 0.3e1_dp
      e_rho = e_rho + (t123 * REAL(t85,KIND=dp) * t101 / 0.8e1_dp + rho * t305 * t102 / 0.8e1_dp  &
      - t124 * t231 * t259 / 0.24e2_dp ) * sx
      t312 = t5 * t33
      t314 = gamma * ndrho
      t315 = t314 * t270
      t317 = t3 * t312 * t315 / 0.9e1_dp
      t319 = t134 * ndrho
      t323 = t155 * ndrho
      t327 = t163 * ndrho
      t331 = t173 * ndrho
      t335 = t186 * ndrho
      t338 = -t24 * t5 * t319 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t28 * t149 * t323  &
      - 0.4e1_dp / 0.27e2_dp * t38 * t39 * t327 - 0.32e2_dp / 0.243e3_dp * t47 * t167  &
      * t331 - 0.80e2_dp / 0.729e3_dp * t57 * t178 * t335
      t340 = t23 * t338 * t83
      t356 = -t67 * t5 * t319 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t70 * t149 * t323  &
      - 0.4e1_dp / 0.27e2_dp * t73 * t39 * t327 - 0.32e2_dp / 0.243e3_dp * t76 * t167  &
      * t331 - 0.80e2_dp / 0.729e3_dp * t79 * t178 * t335
      t358 = t66 * t193 * t356
      t359 = t3 * t5
      t361 = t270 * t103
      t363 = t359 * t319 * t361 / 0.9e1_dp
      t364 = t340 * t103
      t366 = t102 * R * t356
      t367 = t225 * t366
      t388 = t232 * t93 * t97 * t132 * t3 * t33 * t314 * t141 / 0.3e1_dp + 0.3e1_dp  &
      * t240 * t97 * t10 * t338 - 0.3e1_dp * t247 * t97 * t10 * t356 +  &
      t94 * (-t317 - t340 + t358) * t95 * t235
      t390 = t231 * R * t388
      t392 = t84 * t390 / 0.3e1_dp
      t394 = (t317 + t340 - t358 - t363 - t364 + t367 + t392) * t106
      t397 = t106 * br_a1
      t399 = t2 * t5 * t33
      t403 = t338 * t83
      t405 = t274 * t356
      t409 = t397 * t2
      t410 = t312 * gamma
      t414 = ndrho * t141 * t65 * t114
      t423 = (-t317 - t340 + t358 - t363 - t364 + t367 + t392) * t117
      t426 = t117 * br_a1
      t434 = t426 * t2
      t443 = -0.2e1_dp * t394 + t394 * t84 + t397 * t399 * t315 / 0.9e1_dp + t108  &
      * t403 - t108 * t405 - t394 * t66 * t114 - t409 * t410 * t414 /  &
      0.9e1_dp - t108 * t338 * t114 + t283 * t366 + t110 * t390 / 0.3e1_dp +  &
      0.2e1_dp * t423 + t423 * t84 + t426 * t399 * t315 / 0.9e1_dp + t119 * t403  &
      - t119 * t405 + t423 * t66 * t114 + t434 * t410 * t414 / 0.9e1_dp  &
      + t119 * t338 * t114 - t301 * t366 - t120 * t390 / 0.3e1_dp
      e_ndrho = e_ndrho + (rho * t443 * t102 / 0.8e1_dp - t124 * t231 * t388 / 0.24e2_dp ) * sx
      t450 = t6 * t33
      t455 = 0.4e1_dp / 0.9e1_dp * t3 * t450 * gamma * t141 * t109
      t456 = t450 * gamma
      t459 = t31 * t42
      t460 = t459 * gamma
      t463 = t40 * t51
      t464 = t463 * gamma
      t467 = t49 * t61
      t468 = t467 * gamma
      t471 = t59 * t185
      t472 = t471 * gamma
      t475 = 0.4e1_dp / 0.9e1_dp * t24 * t456 + 0.16e2_dp / 0.27e2_dp * t28 * t460 +  &
      0.16e2_dp / 0.27e2_dp * t38 * t464 + 0.128e3_dp / 0.243e3_dp * t47 * t468 + 0.320e3_dp  &
      / 0.729e3_dp * t57 * t472
      t477 = t23 * t475 * t83
      t488 = 0.4e1_dp / 0.9e1_dp * t67 * t456 + 0.16e2_dp / 0.27e2_dp * t70 * t460 +  &
      0.16e2_dp / 0.27e2_dp * t73 * t464 + 0.128e3_dp / 0.243e3_dp * t76 * t468 + 0.320e3_dp  &
      / 0.729e3_dp * t79 * t472
      t490 = t66 * t193 * t488
      t493 = 0.4e1_dp / 0.9e1_dp * t3 * t456 * t361
      t494 = t477 * t103
      t496 = t102 * R * t488
      t497 = t225 * t496
      t499 = t232 * t233 * t96
      t516 = -0.4e1_dp / 0.3e1_dp * t499 * t359 * t134 * t141 + 0.3e1_dp * t240 *  &
      t97 * t10 * t475 - 0.3e1_dp * t247 * t97 * t10 * t488 + t94 * (t455 -  &
      t477 + t490) * t95 * t235
      t518 = t231 * R * t516
      t520 = t84 * t518 / 0.3e1_dp
      t522 = (-t455 + t477 - t490 + t493 - t494 + t497 + t520) * t106
      t525 = t2 * t6
      t526 = t397 * t525
      t527 = t134 * t270
      t530 = t475 * t83
      t532 = t274 * t488
      t545 = (t455 - t477 + t490 + t493 - t494 + t497 + t520) * t117
      t548 = t426 * t525
      t563 = -0.2e1_dp * t522 + t522 * t84 - 0.4e1_dp / 0.9e1_dp * t526 * t527 + t108  &
      * t530 - t108 * t532 - t522 * t66 * t114 + 0.4e1_dp / 0.9e1_dp * t409  &
      * t456 * t361 - t108 * t475 * t114 + t283 * t496 + t110 * t518 /  &
      0.3e1_dp + 0.2e1_dp * t545 + t545 * t84 - 0.4e1_dp / 0.9e1_dp * t548 * t527 +  &
      t119 * t530 - t119 * t532 + t545 * t66 * t114 - 0.4e1_dp / 0.9e1_dp * t434  &
      * t456 * t361 + t119 * t475 * t114 - t301 * t496 - t120 * t518  &
      / 0.3e1_dp
      e_tau = e_tau + (rho * t563 * t102 / 0.8e1_dp - t124 * t231 * t516 / 0.24e2_dp ) * sx
      t572 = t33 * t141 * t109
      t574 = t3 * t6 * t572 / 0.9e1_dp
      t585 = -t24 * t450 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t28 * t459 - 0.4e1_dp /  &
      0.27e2_dp * t38 * t463 - 0.32e2_dp / 0.243e3_dp * t47 * t467 - 0.80e2_dp / 0.729e3_dp  &
      * t57 * t471
      t587 = t23 * t585 * t83
      t598 = -t67 * t450 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t70 * t459 - 0.4e1_dp /  &
      0.27e2_dp * t73 * t463 - 0.32e2_dp / 0.243e3_dp * t76 * t467 - 0.80e2_dp / 0.729e3_dp  &
      * t79 * t471
      t600 = t66 * t193 * t598
      t605 = t3 * t450 * t141 * t109 * t103 / 0.9e1_dp
      t606 = t587 * t103
      t608 = t102 * R * t598
      t609 = t225 * t608
      t628 = t499 * t5 * br_a1 * t2 * t33 * t141 / 0.3e1_dp + 0.3e1_dp * t240 *  &
      t97 * t10 * t585 - 0.3e1_dp * t247 * t97 * t10 * t598 + t94 * (-t574  &
      - t587 + t600) * t95 * t235
      t630 = t231 * R * t628
      t632 = t84 * t630 / 0.3e1_dp
      t634 = (t574 + t587 - t600 - t605 - t606 + t609 + t632) * t106
      t639 = t585 * t83
      t641 = t274 * t598
      t645 = t525 * t33
      t655 = (-t574 - t587 + t600 - t605 - t606 + t609 + t632) * t117
      t672 = -0.2e1_dp * t634 + t634 * t84 + t526 * t572 / 0.9e1_dp + t108 * t639  &
      - t108 * t641 - t634 * t66 * t114 - t397 * t645 * t361 / 0.9e1_dp  &
      - t108 * t585 * t114 + t283 * t608 + t110 * t630 / 0.3e1_dp + 0.2e1_dp *  &
      t655 + t655 * t84 + t548 * t572 / 0.9e1_dp + t119 * t639 - t119 * t641  &
      + t655 * t66 * t114 + t426 * t645 * t361 / 0.9e1_dp + t119 * t585  &
      * t114 - t301 * t608 - t120 * t630 / 0.3e1_dp
      e_laplace_rho = e_laplace_rho + (rho * t672 * t102 / 0.8e1_dp - t124 * t231 * t628 / 0.24e2_dp) * sx
    END IF

  END SUBROUTINE x_br_lsd_y_lte_0_cutoff_R_gt_b


! *****************************************************************************
!> \brief Truncated long range part for y <= 0 and R <= b
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param R ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_lte_0_cutoff_R_lte_b(rho, ndrho, tau, laplace_rho, e_0, &
                                             e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                             sx, R, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, R, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'x_br_lsd_y_lte_0_cutoff_R_lte_b', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t10, t100, t101, t102, t103, t104, t106, t108, t109, &
      t110, t111, t113, t115, t116, t118, t119, t121, t123, t128, t131, t133, &
      t134, t137, t138, t140, t141, t143, t146, t153, t154, t157, t159, t16, &
      t160, t163, t167, t168, t17, t171, t172, t175, t177, t178, t18, t182, &
      t183, t186, t189, t190, t191, t194, t196, t197, t199, t2, t200, t21, &
      t22, t226, t229, t23, t232, t233, t234, t235, t237, t24, t242, t247, &
      t249, t254, t256, t264, t267, t27, t270, t272, t277, t28, t280, t281, &
      t29, t292, t295, t298, t299, t3, t30, t302, t31, t310, t314, t315, &
      t317, t318, t32, t324, t328, t33
    REAL(KIND=dp) :: t332, t336, t339, t34, t341, t342, t359, t362, t368, &
      t369, t37, t38, t383, t385, t387, t39, t392, t395, t398, t4, t40, t402, &
      t404, t409, t42, t420, t427, t428, t429, t43, t432, t443, t444, t446, &
      t450, t451, t454, t455, t458, t459, t46, t462, t463, t466, t468, t469, &
      t47, t48, t481, t484, t487, t49, t5, t50, t501, t504, t506, t51, t511, &
      t514, t517, t52, t521, t525, t529, t540, t547, t548, t549, t552, t56, &
      t566, t57, t578, t58, t580, t581, t59, t593, t596, t6, t61, t612, t613, &
      t614, t616, t618, t62, t623, t626, t629, t638, t65, t654, t655, t656, &
      t659, t66, t67, t70, t73, t76
    REAL(KIND=dp) :: t79, t82, t83, t84, t85, t86, t87, t88, t89, t9, t90, &
      t91, t93, t94, t95, t96, t97, t99

    IF( grad_deriv >= 0 ) THEN
      t1 = pi ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = br_a1 * t2
      t4 = rho ** (0.1e1_dp / 0.3e1_dp)
      t5 = t4 ** 2
      t6 = t5 * rho
      t9 = ndrho ** 2
      t10 = 0.1e1_dp / rho
      t16 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t9 * t10 / 0.4e1_dp) / 0.3e1_dp
      t17 = 0.1e1_dp / t16
      t18 = t6 * t17
      t21 = 0.2e1_dp / 0.3e1_dp * t3 * t18 + br_a2
      t22 = ATAN(t21)
      t23 = -t22 + br_a3
      t24 = br_c1 * t2
      t27 = t1 * pi
      t28 = br_c2 * t27
      t29 = rho ** 2
      t30 = t29 * rho
      t31 = t4 * t30
      t32 = t16 ** 2
      t33 = 0.1e1_dp / t32
      t34 = t31 * t33
      t37 = pi ** 2
      t38 = br_c3 * t37
      t39 = t29 ** 2
      t40 = t39 * rho
      t42 = 0.1e1_dp / t32 / t16
      t43 = t40 * t42
      t46 = t2 * t37
      t47 = br_c4 * t46
      t48 = t39 * t29
      t49 = t5 * t48
      t50 = t32 ** 2
      t51 = 0.1e1_dp / t50
      t52 = t49 * t51
      t56 = t1 * t37 * pi
      t57 = br_c5 * t56
      t58 = t39 ** 2
      t59 = t4 * t58
      t61 = 0.1e1_dp / t50 / t16
      t62 = t59 * t61
      t65 = br_c0 + 0.2e1_dp / 0.3e1_dp * t24 * t18 + 0.4e1_dp / 0.9e1_dp * t28 * t34  &
      + 0.8e1_dp / 0.27e2_dp * t38 * t43 + 0.16e2_dp / 0.81e2_dp * t47 * t52 + 0.32e2_dp  &
      / 0.243e3_dp * t57 * t62
      t66 = t23 * t65
      t67 = br_b1 * t2
      t70 = br_b2 * t27
      t73 = br_b3 * t37
      t76 = br_b4 * t46
      t79 = br_b5 * t56
      t82 = br_b0 + 0.2e1_dp / 0.3e1_dp * t67 * t18 + 0.4e1_dp / 0.9e1_dp * t70 * t34  &
      + 0.8e1_dp / 0.27e2_dp * t73 * t43 + 0.16e2_dp / 0.81e2_dp * t76 * t52 + 0.32e2_dp  &
      / 0.243e3_dp * t79 * t62
      t83 = 0.1e1_dp / t82
      t84 = t66 * t83
      t85 = 8 ** (0.1e1_dp / 0.3e1_dp)
      t86 = t23 ** 2
      t87 = t86 * t23
      t88 = t65 ** 2
      t89 = t88 * t65
      t90 = t87 * t89
      t91 = t82 ** 2
      t93 = 0.1e1_dp / t91 / t82
      t94 = t90 * t93
      t95 = EXP(-t84)
      t96 = 0.1e1_dp / 0.3141592654e1_dp
      t97 = t95 * t96
      t99 = t94 * t97 * t10
      t100 = t99 ** (0.1e1_dp / 0.3e1_dp)
      t101 = 0.1e1_dp / t100
      t102 = REAL(t85,KIND=dp) * t101
      t103 = t102 * R
      t104 = t84 * t103
      t106 = EXP(0.2e1_dp * t104)
      t108 = t106 * R
      t109 = t108 * t23
      t110 = t65 * t83
      t111 = t110 * t102
      t113 = t106 * t23
      t115 = t84 + t104
      t116 = EXP(t115)
      t118 = 0.2e1_dp * t106 - t109 * t111 + t113 * t110 + 0.2e1_dp + t84 + t104  &
      - 0.4e1_dp * t116
      t119 = rho * t118
      t121 = EXP(-t115)
      t123 = t121 * REAL(t85,KIND=dp) * t101
      e_0 = e_0 + (t119 * t123 / 0.8e1_dp) * sx
    END IF
    IF(grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
      t128 = t5 * t17
      t131 = 0.1e1_dp / t4
      t133 = t33 * gamma
      t134 = t133 * t9
      t137 = 0.10e2_dp / 0.9e1_dp * t3 * t128 + t3 * t131 * t134 / 0.18e2_dp
      t138 = t21 ** 2
      t140 = 0.1e1_dp / (0.1e1_dp + t138)
      t141 = t137 * t140
      t143 = t83 * REAL(t85,KIND=dp)
      t146 = t141 * t65 * t143 * t101 * R
      t153 = t4 * t29
      t154 = t153 * t33
      t157 = t4 * rho
      t159 = t42 * gamma
      t160 = t159 * t9
      t163 = t39 * t42
      t167 = t51 * gamma
      t168 = t167 * t9
      t171 = t5 * t40
      t172 = t171 * t51
      t175 = t5 * t39
      t177 = t61 * gamma
      t178 = t177 * t9
      t182 = t4 * t39 * t30
      t183 = t182 * t61
      t186 = t4 * t48
      t189 = 0.1e1_dp / t50 / t32
      t190 = t189 * gamma
      t191 = t190 * t9
      t194 = 0.10e2_dp / 0.9e1_dp * t24 * t128 + t24 * t131 * t134 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t28 * t154 + 0.2e1_dp / 0.27e2_dp * t28 * t157 * t160 +  &
      0.40e2_dp / 0.27e2_dp * t38 * t163 + 0.2e1_dp / 0.27e2_dp * t38 * t30 * t168  &
      + 0.320e3_dp / 0.243e3_dp * t47 * t172 + 0.16e2_dp / 0.243e3_dp * t47 * t175 *  &
      t178 + 0.800e3_dp / 0.729e3_dp * t57 * t183 + 0.40e2_dp / 0.729e3_dp * t57 *  &
      t186 * t191
      t196 = t23 * t194 * t83
      t197 = t196 * t103
      t199 = 0.1e1_dp / t91
      t200 = t66 * t199
      t226 = 0.10e2_dp / 0.9e1_dp * t67 * t128 + t67 * t131 * t134 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t70 * t154 + 0.2e1_dp / 0.27e2_dp * t70 * t157 * t160 +  &
      0.40e2_dp / 0.27e2_dp * t73 * t163 + 0.2e1_dp / 0.27e2_dp * t73 * t30 * t168  &
      + 0.320e3_dp / 0.243e3_dp * t76 * t172 + 0.16e2_dp / 0.243e3_dp * t76 * t175 *  &
      t178 + 0.800e3_dp / 0.729e3_dp * t79 * t183 + 0.40e2_dp / 0.729e3_dp * t79 *  &
      t186 * t191
      t229 = t200 * t102 * R * t226
      t232 = 0.1e1_dp / t100 / t99
      t233 = REAL(t85,KIND=dp) * t232
      t234 = t86 * t89
      t235 = t93 * t95
      t237 = t96 * t10
      t242 = t87 * t88 * t93
      t247 = t91 ** 2
      t249 = t90 / t247
      t254 = t141 * t110
      t256 = t66 * t199 * t226
      t264 = -0.3e1_dp * t234 * t235 * t237 * t141 + 0.3e1_dp * t242 * t97 * t10  &
      * t194 - 0.3e1_dp * t249 * t97 * t10 * t226 + t94 * (t254 - t196 +  &
      t256) * t95 * t237 - t94 * t97 / t29
      t267 = t84 * t233 * R * t264
      t270 = (-0.2e1_dp * t146 + 0.2e1_dp * t197 - 0.2e1_dp * t229 - 0.2e1_dp / 0.3e1_dp  &
      * t267) * t106
      t272 = R * t23
      t277 = t194 * t83
      t280 = t108 * t66
      t281 = t199 * REAL(t85,KIND=dp)
      t292 = t140 * t65 * t83
      t295 = t65 * t199
      t298 = t267 / 0.3e1_dp
      t299 = -t254 + t196 - t256 - t146 + t197 - t229 - t298
      t302 = 0.2e1_dp * t270 - t270 * t272 * t111 + t108 * t141 * t111 - t109  &
      * t277 * t102 + t280 * t281 * t101 * t226 + t280 * t143 * t232 *  &
      t264 / 0.3e1_dp + t270 * t84 - t106 * t137 * t292 + t113 * t277 - t113  &
      * t295 * t226 - t254 + t196 - t256 - t146 + t197 - t229 - t298  &
      - 0.4e1_dp * t299 * t116
      t310 = t119 * t121
      e_rho = e_rho + (t118 * t121 * t102 / 0.8e1_dp + rho * t302 * t123 / 0.8e1_dp - t119  &
      * t299 * t123 / 0.8e1_dp - t310 * t233 * t264 / 0.24e2_dp) * sx
      t314 = t3 * t5
      t315 = t133 * ndrho
      t317 = t292 * t103
      t318 = t314 * t315 * t317
      t324 = t159 * ndrho
      t328 = t167 * ndrho
      t332 = t177 * ndrho
      t336 = t190 * ndrho
      t339 = -t24 * t5 * t315 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t28 * t153 * t324  &
      - 0.4e1_dp / 0.27e2_dp * t38 * t39 * t328 - 0.32e2_dp / 0.243e3_dp * t47 * t171  &
      * t332 - 0.80e2_dp / 0.729e3_dp * t57 * t182 * t336
      t341 = t23 * t339 * t83
      t342 = t341 * t103
      t359 = -t67 * t5 * t315 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t70 * t153 * t324  &
      - 0.4e1_dp / 0.27e2_dp * t73 * t39 * t328 - 0.32e2_dp / 0.243e3_dp * t76 * t171  &
      * t332 - 0.80e2_dp / 0.729e3_dp * t79 * t182 * t336
      t362 = t200 * t102 * R * t359
      t368 = gamma * ndrho
      t369 = t368 * t140
      t383 = t368 * t292
      t385 = t3 * t5 * t33 * t383 / 0.9e1_dp
      t387 = t66 * t199 * t359
      t392 = t234 * t93 * t97 * t131 * t3 * t33 * t369 / 0.3e1_dp + 0.3e1_dp *  &
      t242 * t97 * t10 * t339 - 0.3e1_dp * t249 * t97 * t10 * t359 + t94 *  &
      (-t385 - t341 + t387) * t95 * t237
      t395 = t84 * t233 * R * t392
      t398 = (0.2e1_dp / 0.9e1_dp * t318 + 0.2e1_dp * t342 - 0.2e1_dp * t362 - 0.2e1_dp  &
      / 0.3e1_dp * t395) * t106
      t402 = t108 * br_a1
      t404 = t2 * t5 * t33
      t409 = t339 * t83
      t420 = t106 * br_a1
      t427 = t318 / 0.9e1_dp
      t428 = t395 / 0.3e1_dp
      t429 = t385 + t341 - t387 + t427 + t342 - t362 - t428
      t432 = 0.2e1_dp * t398 - t398 * t272 * t111 - t402 * t404 * t369 * t111  &
      / 0.9e1_dp - t109 * t409 * t102 + t280 * t281 * t101 * t359 + t280  &
      * t143 * t232 * t392 / 0.3e1_dp + t398 * t84 + t420 * t404 * t383 / 0.9e1_dp  &
      + t113 * t409 - t113 * t295 * t359 + t385 + t341 - t387 + t427  &
      + t342 - t362 - t428 - 0.4e1_dp * t429 * t116
      e_ndrho = e_ndrho + (rho * t432 * t123 / 0.8e1_dp - t119 * t429 * t123 / 0.8e1_dp - t310  &
      * t233 * t392 / 0.24e2_dp ) * sx
      t443 = t6 * t33
      t444 = t443 * gamma
      t446 = t3 * t444 * t317
      t450 = t31 * t42
      t451 = t450 * gamma
      t454 = t40 * t51
      t455 = t454 * gamma
      t458 = t49 * t61
      t459 = t458 * gamma
      t462 = t59 * t189
      t463 = t462 * gamma
      t466 = 0.4e1_dp / 0.9e1_dp * t24 * t444 + 0.16e2_dp / 0.27e2_dp * t28 * t451 +  &
      0.16e2_dp / 0.27e2_dp * t38 * t455 + 0.128e3_dp / 0.243e3_dp * t47 * t459 + 0.320e3_dp  &
      / 0.729e3_dp * t57 * t463
      t468 = t23 * t466 * t83
      t469 = t468 * t103
      t481 = 0.4e1_dp / 0.9e1_dp * t67 * t444 + 0.16e2_dp / 0.27e2_dp * t70 * t451 +  &
      0.16e2_dp / 0.27e2_dp * t73 * t455 + 0.128e3_dp / 0.243e3_dp * t76 * t459 + 0.320e3_dp  &
      / 0.729e3_dp * t79 * t463
      t484 = t200 * t102 * R * t481
      t487 = t234 * t235 * t96
      t501 = gamma * t140
      t504 = 0.4e1_dp / 0.9e1_dp * t3 * t443 * t501 * t110
      t506 = t66 * t199 * t481
      t511 = -0.4e1_dp / 0.3e1_dp * t487 * t314 * t133 * t140 + 0.3e1_dp * t242 *  &
      t97 * t10 * t466 - 0.3e1_dp * t249 * t97 * t10 * t481 + t94 * (t504 -  &
      t468 + t506) * t95 * t237
      t514 = t84 * t233 * R * t511
      t517 = (-0.8e1_dp / 0.9e1_dp * t446 + 0.2e1_dp * t469 - 0.2e1_dp * t484 - 0.2e1_dp  &
      / 0.3e1_dp * t514) * t106
      t521 = t2 * t6
      t525 = t143 * t101
      t529 = t466 * t83
      t540 = t420 * t521
      t547 = 0.4e1_dp / 0.9e1_dp * t446
      t548 = t514 / 0.3e1_dp
      t549 = -t504 + t468 - t506 - t547 + t469 - t484 - t548
      t552 = 0.2e1_dp * t517 - t517 * t272 * t111 + 0.4e1_dp / 0.9e1_dp * t402 * t521  &
      * t33 * t501 * t65 * t525 - t109 * t529 * t102 + t280 * t281 *  &
      t101 * t481 + t280 * t143 * t232 * t511 / 0.3e1_dp + t517 * t84 - 0.4e1_dp  &
      / 0.9e1_dp * t540 * t133 * t292 + t113 * t529 - t113 * t295 * t481  &
      - t504 + t468 - t506 - t547 + t469 - t484 - t548 - 0.4e1_dp * t549  &
      * t116
      e_tau = e_tau + (rho * t552 * t123 / 0.8e1_dp - t119 * t549 * t123 / 0.8e1_dp - t310  &
      * t233 * t511 / 0.24e2_dp ) * sx
      t566 = t3 * t443 * t140 * t110 * t103
      t578 = -t24 * t443 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t28 * t450 - 0.4e1_dp /  &
      0.27e2_dp * t38 * t454 - 0.32e2_dp / 0.243e3_dp * t47 * t458 - 0.80e2_dp / 0.729e3_dp  &
      * t57 * t462
      t580 = t23 * t578 * t83
      t581 = t580 * t103
      t593 = -t67 * t443 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t70 * t450 - 0.4e1_dp /  &
      0.27e2_dp * t73 * t454 - 0.32e2_dp / 0.243e3_dp * t76 * t458 - 0.80e2_dp / 0.729e3_dp  &
      * t79 * t462
      t596 = t200 * t102 * R * t593
      t612 = t3 * t6
      t613 = t33 * t140
      t614 = t613 * t110
      t616 = t612 * t614 / 0.9e1_dp
      t618 = t66 * t199 * t593
      t623 = t487 * t5 * br_a1 * t2 * t33 * t140 / 0.3e1_dp + 0.3e1_dp * t242 *  &
      t97 * t10 * t578 - 0.3e1_dp * t249 * t97 * t10 * t593 + t94 * (-t616  &
      - t580 + t618) * t95 * t237
      t626 = t84 * t233 * R * t623
      t629 = (0.2e1_dp / 0.9e1_dp * t566 + 0.2e1_dp * t581 - 0.2e1_dp * t596 - 0.2e1_dp  &
      / 0.3e1_dp * t626) * t106
      t638 = t578 * t83
      t654 = t566 / 0.9e1_dp
      t655 = t626 / 0.3e1_dp
      t656 = t616 + t580 - t618 + t654 + t581 - t596 - t655
      t659 = 0.2e1_dp * t629 - t629 * t272 * t111 - t108 * t612 * t613 * t65  &
      * t525 / 0.9e1_dp - t109 * t638 * t102 + t280 * t281 * t101 * t593 +  &
      t280 * t143 * t232 * t623 / 0.3e1_dp + t629 * t84 + t540 * t614 / 0.9e1_dp  &
      + t113 * t638 - t113 * t295 * t593 + t616 + t580 - t618 + t654  &
      + t581 - t596 - t655 - 0.4e1_dp * t656 * t116
      e_laplace_rho = e_laplace_rho + (rho * t659 * t123 / 0.8e1_dp - t119 * t656 * t123 / 0.8e1_dp - t310  &
      * t233 * t623 / 0.24e2_dp ) * sx
    END IF

  END SUBROUTINE x_br_lsd_y_lte_0_cutoff_R_lte_b


! *****************************************************************************
!> \brief Truncated long range part for y > 0
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param R ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_gt_0_cutoff(rho, ndrho, tau, laplace_rho, e_0, &
                                    e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                    sx, R, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, R, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'x_br_lsd_y_gt_0_cutoff', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: brval, t1, t10, t103, t104, t11, t111, t116, t14, t15, &
      t2, t21, t25, t26, t28, t3, t31, t33, t37, t4, t44, t45, t46, t5, t50, &
      t56, t58, t6, t62, t65, t69, t71, t75, t77, t8, t81, t84, t85, t9

      t1 = 8 ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = 0.1e1_dp / br_BB
      t4 = pi ** (0.1e1_dp / 0.3e1_dp)
      t5 = t4 ** 2
      t6 = 0.1e1_dp / t5
      t8 = rho ** (0.1e1_dp / 0.3e1_dp)
      t9 = t8 ** 2
      t10 = t9 * rho
      t11 = 0.1e1_dp / t10
      t14 = ndrho ** 2
      t15 = 0.1e1_dp / rho
      t21 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t14 * t15 / 0.4e1_dp) / 0.3e1_dp
      t25 = br_BB ** 2
      t26 = t4 * pi
      t28 = rho ** 2
      t31 = t21 ** 2
      t33 = t8 * t28 * rho / t31
      t37 = SQRT(0.10e1_dp + 0.4e1_dp / 0.9e1_dp * t25 * t26 * t33)
      t44 = LOG(0.1500000000000000e1_dp * t3 * t6 * t11 * t21 + 0.3e1_dp / 0.2e1_dp  &
      * t37 * t3 * t6 * t11 * t21)
      t45 = t44 + 0.2e1_dp
      t46 = t45 ** 2
      t50 = t10 / t21
      t56 = pi ** 2
      t58 = t28 ** 2
      t62 = t58 * rho / t31 / t21
      t65 = t5 * t56
      t69 = t31 ** 2
      t71 = t9 * t58 * t28 / t69
      t75 = t4 * t56 * pi
      t77 = t58 ** 2
      t81 = t8 * t77 / t69 / t21
      t84 = br_d0 + 0.2e1_dp / 0.3e1_dp * br_d1 * t5 * t50 + 0.4e1_dp / 0.9e1_dp * br_d2  &
      * t26 * t33 + 0.8e1_dp / 0.27e2_dp * br_d3 * t56 * t62 + 0.16e2_dp / 0.81e2_dp  &
      * br_d4 * t65 * t71 + 0.32e2_dp / 0.243e3_dp * br_d5 * t75 * t81
      t85 = t84 ** 2
      t103 = br_e0 + 0.2e1_dp / 0.3e1_dp * br_e1 * t5 * t50 + 0.4e1_dp / 0.9e1_dp * br_e2  &
      * t26 * t33 + 0.8e1_dp / 0.27e2_dp * br_e3 * t56 * t62 + 0.16e2_dp / 0.81e2_dp  &
      * br_e4 * t65 * t71 + 0.32e2_dp / 0.243e3_dp * br_e5 * t75 * t81
      t104 = t103 ** 2
      t111 = EXP(-t45 * t84 / t103)
      t116 = (t46 * t45 * t85 * t84 / t104 / t103 * t111 / 0.3141592654e1_dp  &
      * t15) ** (0.1e1_dp / 0.3e1_dp)
      brval = REAL(t2,KIND=dp) * t116 / 0.8e1_dp

    IF(R>brval) THEN
      CALL x_br_lsd_y_gt_0_cutoff_R_gt_b(rho, ndrho, tau, laplace_rho, e_0, &
                                           e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                           sx, R, gamma, grad_deriv)
    ELSE
      CALL x_br_lsd_y_gt_0_cutoff_R_lte_b(rho, ndrho, tau, laplace_rho, e_0, &
                                          e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                          sx, R, gamma, grad_deriv)
    END IF

  END SUBROUTINE x_br_lsd_y_gt_0_cutoff

! *****************************************************************************
!> \brief Truncated long range part for y > 0 and R > b
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param R ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_gt_0_cutoff_R_gt_b(rho, ndrho, tau, laplace_rho, e_0, &
                                            e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                            sx, R, gamma, grad_deriv)

    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, R, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'x_br_lsd_y_gt_0_cutoff_R_gt_b', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t100, t101, t102, t103, t104, t105, t106, t108, &
      t109, t110, t111, t112, t114, t115, t116, t117, t118, t119, t12, t121, &
      t123, t124, t125, t129, t13, t132, t134, t135, t138, t139, t145, t152, &
      t155, t158, t159, t162, t164, t165, t176, t179, t180, t181, t182, t183, &
      t186, t188, t189, t19, t197, t2, t20, t201, t202, t205, t206, t209, &
      t211, t212, t216, t217, t220, t223, t224, t225, t228, t23, t230, t231, &
      t24, t25, t257, t259, t26, t261, t262, t263, t265, t266, t269, t27, &
      t272, t273, t278, t28, t283, t285, t29, t297, t299, t3, t30, t301, &
      t303, t306, t307, t308, t31, t310, t312
    REAL(KIND=dp) :: t313, t321, t326, t329, t339, t343, t35, t351, t354, &
      t355, t36, t363, t364, t365, t367, t37, t371, t375, t379, t383, t386, &
      t388, t4, t404, t406, t408, t409, t41, t411, t412, t42, t428, t43, &
      t430, t432, t434, t437, t439, t44, t441, t45, t453, t456, t46, t469, &
      t479, t480, t485, t486, t487, t49, t490, t491, t494, t495, t498, t499, &
      t5, t502, t503, t506, t508, t519, t52, t521, t523, t524, t526, t527, &
      t53, t54, t543, t545, t547, t549, t55, t552, t554, t556, t568, t57, &
      t571, t58, t584, t599, t6, t600, t601, t61, t612, t614, t62, t625, &
      t627, t629, t63, t630, t632, t633, t64, t649
    REAL(KIND=dp) :: t65, t651, t653, t655, t658, t66, t660, t662, t67, t674, &
      t677, t690, t7, t71, t72, t73, t74, t76, t77, t8, t80, t81, t82, t85, &
      t88, t9, t91, t94, t97, t98, t99

    IF( grad_deriv >= 0 ) THEN
      t1 = 0.1e1_dp / br_BB
      t2 = pi ** (0.1e1_dp / 0.3e1_dp)
      t3 = t2 ** 2
      t4 = 0.1e1_dp / t3
      t5 = t1 * t4
      t6 = rho ** (0.1e1_dp / 0.3e1_dp)
      t7 = t6 ** 2
      t8 = t7 * rho
      t9 = 0.1e1_dp / t8
      t12 = ndrho ** 2
      t13 = 0.1e1_dp / rho
      t19 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t12 * t13 / 0.4e1_dp) / 0.3e1_dp
      t20 = t9 * t19
      t23 = br_BB ** 2
      t24 = t2 * pi
      t25 = t23 * t24
      t26 = rho ** 2
      t27 = t26 * rho
      t28 = t6 * t27
      t29 = t19 ** 2
      t30 = 0.1e1_dp / t29
      t31 = t28 * t30
      t35 = SQRT(0.10e1_dp + 0.4e1_dp / 0.9e1_dp * t25 * t31)
      t36 = t35 * t1
      t37 = t4 * t9
      t41 = 0.1500000000000000e1_dp * t5 * t20 + 0.3e1_dp / 0.2e1_dp * t36 * t37 *  &
      t19
      t42 = LOG(t41)
      t43 = t42 + 0.2e1_dp
      t44 = br_d1 * t3
      t45 = 0.1e1_dp / t19
      t46 = t8 * t45
      t49 = br_d2 * t24
      t52 = pi ** 2
      t53 = br_d3 * t52
      t54 = t26 ** 2
      t55 = t54 * rho
      t57 = 0.1e1_dp / t29 / t19
      t58 = t55 * t57
      t61 = t3 * t52
      t62 = br_d4 * t61
      t63 = t54 * t26
      t64 = t7 * t63
      t65 = t29 ** 2
      t66 = 0.1e1_dp / t65
      t67 = t64 * t66
      t71 = t2 * t52 * pi
      t72 = br_d5 * t71
      t73 = t54 ** 2
      t74 = t6 * t73
      t76 = 0.1e1_dp / t65 / t19
      t77 = t74 * t76
      t80 = br_d0 + 0.2e1_dp / 0.3e1_dp * t44 * t46 + 0.4e1_dp / 0.9e1_dp * t49 * t31  &
      + 0.8e1_dp / 0.27e2_dp * t53 * t58 + 0.16e2_dp / 0.81e2_dp * t62 * t67 + 0.32e2_dp  &
      / 0.243e3_dp * t72 * t77
      t81 = t43 * t80
      t82 = br_e1 * t3
      t85 = br_e2 * t24
      t88 = br_e3 * t52
      t91 = br_e4 * t61
      t94 = br_e5 * t71
      t97 = br_e0 + 0.2e1_dp / 0.3e1_dp * t82 * t46 + 0.4e1_dp / 0.9e1_dp * t85 * t31  &
      + 0.8e1_dp / 0.27e2_dp * t88 * t58 + 0.16e2_dp / 0.81e2_dp * t91 * t67 + 0.32e2_dp  &
      / 0.243e3_dp * t94 * t77
      t98 = 0.1e1_dp / t97
      t99 = t81 * t98
      t100 = 8 ** (0.1e1_dp / 0.3e1_dp)
      t101 = t43 ** 2
      t102 = t101 * t43
      t103 = t80 ** 2
      t104 = t103 * t80
      t105 = t102 * t104
      t106 = t97 ** 2
      t108 = 0.1e1_dp / t106 / t97
      t109 = t105 * t108
      t110 = EXP(-t99)
      t111 = 0.1e1_dp / 0.3141592654e1_dp
      t112 = t110 * t111
      t114 = t109 * t112 * t13
      t115 = t114 ** (0.1e1_dp / 0.3e1_dp)
      t116 = 0.1e1_dp / t115
      t117 = REAL(t100,KIND=dp) * t116
      t118 = t117 * R
      t119 = t99 * t118
      t121 = EXP(t99 - t119)
      t123 = t121 * t43
      t124 = t80 * t98
      t125 = t123 * t124
      t129 = t98 * REAL(t100,KIND=dp) * t116 * R
      t132 = EXP(-t99 - t119)
      t134 = t132 * t43
      t135 = t134 * t124
      t138 = -0.2e1_dp * t121 + t125 - t123 * t80 * t129 + 0.2e1_dp * t132 + t135  &
      + t134 * t80 * t129
      t139 = rho * t138
      e_0 = e_0 + (t139 * t117 / 0.8e1_dp) * sx
    END IF
    IF(grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
      t145 = 0.1e1_dp / t7 / t26
      t152 = 0.1e1_dp / t7 / t27 * gamma * t12
      t155 = 0.1e1_dp / t35
      t158 = t6 * t26
      t159 = t158 * t30
      t162 = t6 * rho
      t164 = t57 * gamma
      t165 = t164 * t12
      t176 = t36 * t4
      t179 = -0.2500000000e1_dp * t5 * t145 * t19 - 0.1250000000e0_dp * t5 * t152  &
      + 0.3e1_dp / 0.4e1_dp * t155 * t1 * t4 * t20 * (0.40e2_dp / 0.27e2_dp * t25  &
      * t159 + 0.2e1_dp / 0.27e2_dp * t25 * t162 * t165) - 0.5e1_dp / 0.2e1_dp * t36  &
      * t4 * t145 * t19 - t176 * t152 / 0.8e1_dp
      t180 = 0.1e1_dp / t41
      t181 = t179 * t180
      t182 = t181 * t124
      t183 = t7 * t45
      t186 = 0.1e1_dp / t6
      t188 = t30 * gamma
      t189 = t188 * t12
      t197 = t54 * t57
      t201 = t66 * gamma
      t202 = t201 * t12
      t205 = t7 * t55
      t206 = t205 * t66
      t209 = t7 * t54
      t211 = t76 * gamma
      t212 = t211 * t12
      t216 = t6 * t54 * t27
      t217 = t216 * t76
      t220 = t6 * t63
      t223 = 0.1e1_dp / t65 / t29
      t224 = t223 * gamma
      t225 = t224 * t12
      t228 = 0.10e2_dp / 0.9e1_dp * t44 * t183 + t44 * t186 * t189 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t49 * t159 + 0.2e1_dp / 0.27e2_dp * t49 * t162 * t165 +  &
      0.40e2_dp / 0.27e2_dp * t53 * t197 + 0.2e1_dp / 0.27e2_dp * t53 * t27 * t202  &
      + 0.320e3_dp / 0.243e3_dp * t62 * t206 + 0.16e2_dp / 0.243e3_dp * t62 * t209 *  &
      t212 + 0.800e3_dp / 0.729e3_dp * t72 * t217 + 0.40e2_dp / 0.729e3_dp * t72 *  &
      t220 * t225
      t230 = t43 * t228 * t98
      t231 = 0.1e1_dp / t106
      t257 = 0.10e2_dp / 0.9e1_dp * t82 * t183 + t82 * t186 * t189 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t85 * t159 + 0.2e1_dp / 0.27e2_dp * t85 * t162 * t165 +  &
      0.40e2_dp / 0.27e2_dp * t88 * t197 + 0.2e1_dp / 0.27e2_dp * t88 * t27 * t202  &
      + 0.320e3_dp / 0.243e3_dp * t91 * t206 + 0.16e2_dp / 0.243e3_dp * t91 * t209 *  &
      t212 + 0.800e3_dp / 0.729e3_dp * t94 * t217 + 0.40e2_dp / 0.729e3_dp * t94 *  &
      t220 * t225
      t259 = t81 * t231 * t257
      t261 = t181 * t80 * t129
      t262 = t230 * t118
      t263 = t81 * t231
      t265 = t117 * R * t257
      t266 = t263 * t265
      t269 = REAL(t100,KIND=dp) / t115 / t114
      t272 = t101 * t104 * t108 * t110
      t273 = t111 * t13
      t278 = t102 * t103 * t108
      t283 = t106 ** 2
      t285 = t105 / t283
      t297 = 0.3e1_dp * t272 * t273 * t181 + 0.3e1_dp * t278 * t112 * t13 * t228  &
      - 0.3e1_dp * t285 * t112 * t13 * t257 + t109 * (-t182 - t230 + t259)  &
      * t110 * t273 - t109 * t112 / t26
      t299 = t269 * R * t297
      t301 = t99 * t299 / 0.3e1_dp
      t303 = (t182 + t230 - t259 - t261 - t262 + t266 + t301) * t121
      t306 = t121 * t179
      t307 = t180 * t80
      t308 = t307 * t98
      t310 = t228 * t98
      t312 = t80 * t231
      t313 = t312 * t257
      t321 = t123 * t312
      t326 = (-t182 - t230 + t259 - t261 - t262 + t266 + t301) * t132
      t329 = t132 * t179
      t339 = t134 * t312
      t343 = -0.2e1_dp * t303 + t303 * t99 + t306 * t308 + t123 * t310 - t123  &
      * t313 - t303 * t81 * t129 - t306 * t307 * t129 - t123 * t228 *  &
      t129 + t321 * t265 + t125 * t299 / 0.3e1_dp + 0.2e1_dp * t326 + t326 * t99  &
      + t329 * t308 + t134 * t310 - t134 * t313 + t326 * t81 * t129 +  &
      t329 * t307 * t129 + t134 * t228 * t129 - t339 * t265 - t135 * t299  &
      / 0.3e1_dp
      e_rho = e_rho + (t138 * REAL(t100,KIND=dp) * t116 / 0.8e1_dp + rho * t343 * t117 / 0.8e1_dp  &
      - t139 * t269 * t297 / 0.24e2_dp ) * sx
      t351 = t145 * gamma * ndrho
      t354 = t155 * br_BB
      t355 = t354 * t3
      t363 = 0.2500000000000000e0_dp * t5 * t351 - t355 * t7 * t30 * gamma * ndrho  &
      / 0.9e1_dp + t176 * t351 / 0.4e1_dp
      t364 = t363 * t180
      t365 = t364 * t124
      t367 = t188 * ndrho
      t371 = t164 * ndrho
      t375 = t201 * ndrho
      t379 = t211 * ndrho
      t383 = t224 * ndrho
      t386 = -t44 * t7 * t367 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t49 * t158 * t371  &
      - 0.4e1_dp / 0.27e2_dp * t53 * t54 * t375 - 0.32e2_dp / 0.243e3_dp * t62 * t205  &
      * t379 - 0.80e2_dp / 0.729e3_dp * t72 * t216 * t383
      t388 = t43 * t386 * t98
      t404 = -t82 * t7 * t367 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t85 * t158 * t371  &
      - 0.4e1_dp / 0.27e2_dp * t88 * t54 * t375 - 0.32e2_dp / 0.243e3_dp * t91 * t205  &
      * t379 - 0.80e2_dp / 0.729e3_dp * t94 * t216 * t383
      t406 = t81 * t231 * t404
      t408 = t364 * t80 * t129
      t409 = t388 * t118
      t411 = t117 * R * t404
      t412 = t263 * t411
      t428 = 0.3e1_dp * t272 * t273 * t364 + 0.3e1_dp * t278 * t112 * t13 * t386  &
      - 0.3e1_dp * t285 * t112 * t13 * t404 + t109 * (-t365 - t388 + t406)  &
      * t110 * t273
      t430 = t269 * R * t428
      t432 = t99 * t430 / 0.3e1_dp
      t434 = (t365 + t388 - t406 - t408 - t409 + t412 + t432) * t121
      t437 = t121 * t363
      t439 = t386 * t98
      t441 = t312 * t404
      t453 = (-t365 - t388 + t406 - t408 - t409 + t412 + t432) * t132
      t456 = t132 * t363
      t469 = -0.2e1_dp * t434 + t434 * t99 + t437 * t308 + t123 * t439 - t123  &
      * t441 - t434 * t81 * t129 - t437 * t307 * t129 - t123 * t386 *  &
      t129 + t321 * t411 + t125 * t430 / 0.3e1_dp + 0.2e1_dp * t453 + t453 * t99  &
      + t456 * t308 + t134 * t439 - t134 * t441 + t453 * t81 * t129 +  &
      t456 * t307 * t129 + t134 * t386 * t129 - t339 * t411 - t135 * t430  &
      / 0.3e1_dp
      e_ndrho = e_ndrho+ (rho * t469 * t117 / 0.8e1_dp - t139 * t269 * t428 / 0.24e2_dp ) * sx
      t479 = t8 * t30
      t480 = t479 * gamma
      t485 = -0.1000000000e1_dp * t5 * t9 * gamma + 0.4e1_dp / 0.9e1_dp * t355 * t480  &
      - t36 * t37 * gamma
      t486 = t485 * t180
      t487 = t486 * t124
      t490 = t28 * t57
      t491 = t490 * gamma
      t494 = t55 * t66
      t495 = t494 * gamma
      t498 = t64 * t76
      t499 = t498 * gamma
      t502 = t74 * t223
      t503 = t502 * gamma
      t506 = 0.4e1_dp / 0.9e1_dp * t44 * t480 + 0.16e2_dp / 0.27e2_dp * t49 * t491 +  &
      0.16e2_dp / 0.27e2_dp * t53 * t495 + 0.128e3_dp / 0.243e3_dp * t62 * t499 + 0.320e3_dp  &
      / 0.729e3_dp * t72 * t503
      t508 = t43 * t506 * t98
      t519 = 0.4e1_dp / 0.9e1_dp * t82 * t480 + 0.16e2_dp / 0.27e2_dp * t85 * t491 +  &
      0.16e2_dp / 0.27e2_dp * t88 * t495 + 0.128e3_dp / 0.243e3_dp * t91 * t499 + 0.320e3_dp  &
      / 0.729e3_dp * t94 * t503
      t521 = t81 * t231 * t519
      t523 = t486 * t80 * t129
      t524 = t508 * t118
      t526 = t117 * R * t519
      t527 = t263 * t526
      t543 = 0.3e1_dp * t272 * t273 * t486 + 0.3e1_dp * t278 * t112 * t13 * t506  &
      - 0.3e1_dp * t285 * t112 * t13 * t519 + t109 * (-t487 - t508 + t521)  &
      * t110 * t273
      t545 = t269 * R * t543
      t547 = t99 * t545 / 0.3e1_dp
      t549 = (t487 + t508 - t521 - t523 - t524 + t527 + t547) * t121
      t552 = t121 * t485
      t554 = t506 * t98
      t556 = t312 * t519
      t568 = (-t487 - t508 + t521 - t523 - t524 + t527 + t547) * t132
      t571 = t132 * t485
      t584 = -0.2e1_dp * t549 + t549 * t99 + t552 * t308 + t123 * t554 - t123  &
      * t556 - t549 * t81 * t129 - t552 * t307 * t129 - t123 * t506 *  &
      t129 + t321 * t526 + t125 * t545 / 0.3e1_dp + 0.2e1_dp * t568 + t568 * t99  &
      + t571 * t308 + t134 * t554 - t134 * t556 + t568 * t81 * t129 +  &
      t571 * t307 * t129 + t134 * t506 * t129 - t339 * t526 - t135 * t545  &
      / 0.3e1_dp
      e_tau = e_tau + (rho * t584 * t117 / 0.8e1_dp - t139 * t269 * t543 / 0.24e2_dp) * sx
      t599 = 0.2500000000000000e0_dp * t5 * t9 - t354 * t3 * t8 * t30 / 0.9e1_dp  &
      + t36 * t37 / 0.4e1_dp
      t600 = t599 * t180
      t601 = t600 * t124
      t612 = -t44 * t479 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t49 * t490 - 0.4e1_dp /  &
      0.27e2_dp * t53 * t494 - 0.32e2_dp / 0.243e3_dp * t62 * t498 - 0.80e2_dp / 0.729e3_dp  &
      * t72 * t502
      t614 = t43 * t612 * t98
      t625 = -t82 * t479 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t85 * t490 - 0.4e1_dp /  &
      0.27e2_dp * t88 * t494 - 0.32e2_dp / 0.243e3_dp * t91 * t498 - 0.80e2_dp / 0.729e3_dp  &
      * t94 * t502
      t627 = t81 * t231 * t625
      t629 = t600 * t80 * t129
      t630 = t614 * t118
      t632 = t117 * R * t625
      t633 = t263 * t632
      t649 = 0.3e1_dp * t272 * t273 * t600 + 0.3e1_dp * t278 * t112 * t13 * t612  &
      - 0.3e1_dp * t285 * t112 * t13 * t625 + t109 * (-t601 - t614 + t627)  &
      * t110 * t273
      t651 = t269 * R * t649
      t653 = t99 * t651 / 0.3e1_dp
      t655 = (t601 + t614 - t627 - t629 - t630 + t633 + t653) * t121
      t658 = t121 * t599
      t660 = t612 * t98
      t662 = t312 * t625
      t674 = (-t601 - t614 + t627 - t629 - t630 + t633 + t653) * t132
      t677 = t132 * t599
      t690 = -0.2e1_dp * t655 + t655 * t99 + t658 * t308 + t123 * t660 - t123  &
      * t662 - t655 * t81 * t129 - t658 * t307 * t129 - t123 * t612 *  &
      t129 + t321 * t632 + t125 * t651 / 0.3e1_dp + 0.2e1_dp * t674 + t674 * t99  &
      + t677 * t308 + t134 * t660 - t134 * t662 + t674 * t81 * t129 +  &
      t677 * t307 * t129 + t134 * t612 * t129 - t339 * t632 - t135 * t651  &
      / 0.3e1_dp
      e_laplace_rho = e_laplace_rho + (rho * t690 * t117 / 0.8e1_dp - t139 * t269 * t649 / 0.24e2_dp) * sx
    END IF

  END SUBROUTINE x_br_lsd_y_gt_0_cutoff_R_gt_b


! *****************************************************************************
!> \brief Truncated long range part for y > 0 and R <= b
!> \param rho ...
!> \param ndrho ...
!> \param tau ...
!> \param laplace_rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_tau ...
!> \param e_laplace_rho ...
!> \param sx ...
!> \param R ...
!> \param gamma ...
!> \param grad_deriv ...
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE x_br_lsd_y_gt_0_cutoff_R_lte_b(rho, ndrho, tau, laplace_rho, e_0, &
                                            e_rho, e_ndrho, e_tau, e_laplace_rho,&
                                            sx, R, gamma, grad_deriv)
    REAL(dp), INTENT(IN)                     :: rho, ndrho, tau, laplace_rho
    REAL(dp), INTENT(INOUT)                  :: e_0, e_rho, e_ndrho, e_tau, &
                                                e_laplace_rho
    REAL(dp), INTENT(IN)                     :: sx, R, gamma
    INTEGER, INTENT(IN)                      :: grad_deriv

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'x_br_lsd_y_gt_0_cutoff_R_lte_b', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t100, t101, t102, t103, t104, t105, t106, t108, &
      t109, t110, t111, t112, t114, t115, t116, t117, t118, t119, t12, t121, &
      t123, t124, t125, t126, t128, t13, t130, t131, t133, t134, t136, t138, &
      t144, t151, t154, t157, t158, t161, t163, t164, t175, t178, t179, t180, &
      t182, t184, t185, t187, t19, t190, t192, t193, t2, t20, t201, t205, &
      t206, t209, t210, t213, t215, t216, t220, t221, t224, t227, t228, t229, &
      t23, t232, t234, t235, t237, t238, t24, t25, t26, t264, t267, t27, &
      t270, t271, t274, t275, t28, t280, t285, t287, t29, t292, t294, t3, &
      t30, t302, t305, t308, t31, t310, t315
    REAL(KIND=dp) :: t318, t319, t330, t333, t336, t337, t340, t348, t35, &
      t353, t356, t357, t36, t365, t366, t368, t37, t371, t375, t379, t383, &
      t387, t390, t392, t393, t4, t41, t410, t413, t42, t426, t428, t43, &
      t433, t436, t439, t44, t445, t45, t46, t461, t462, t465, t479, t480, &
      t485, t486, t488, t49, t492, t493, t496, t497, t5, t500, t501, t504, &
      t505, t508, t510, t511, t52, t523, t526, t53, t539, t54, t541, t546, &
      t549, t55, t552, t558, t57, t574, t575, t578, t58, t597, t598, t6, &
      t600, t61, t612, t614, t615, t62, t627, t63, t630, t64, t643, t645, &
      t65, t650, t653, t656, t66, t662, t67, t678, t679
    REAL(KIND=dp) :: t682, t7, t71, t72, t73, t74, t76, t77, t8, t80, t81, &
      t82, t85, t88, t9, t91, t94, t97, t98, t99

    IF( grad_deriv >= 0 ) THEN
      t1 = 0.1e1_dp / br_BB
      t2 = pi ** (0.1e1_dp / 0.3e1_dp)
      t3 = t2 ** 2
      t4 = 0.1e1_dp / t3
      t5 = t1 * t4
      t6 = rho ** (0.1e1_dp / 0.3e1_dp)
      t7 = t6 ** 2
      t8 = t7 * rho
      t9 = 0.1e1_dp / t8
      t12 = ndrho ** 2
      t13 = 0.1e1_dp / rho
      t19 = laplace_rho / 0.6e1_dp - gamma * (REAL(2 * tau,KIND=dp) - t12 * t13 / 0.4e1_dp) / 0.3e1_dp
      t20 = t9 * t19
      t23 = br_BB ** 2
      t24 = t2 * pi
      t25 = t23 * t24
      t26 = rho ** 2
      t27 = t26 * rho
      t28 = t6 * t27
      t29 = t19 ** 2
      t30 = 0.1e1_dp / t29
      t31 = t28 * t30
      t35 = SQRT(0.10e1_dp + 0.4e1_dp / 0.9e1_dp * t25 * t31)
      t36 = t35 * t1
      t37 = t4 * t9
      t41 = 0.1500000000000000e1_dp * t5 * t20 + 0.3e1_dp / 0.2e1_dp * t36 * t37 *  &
      t19
      t42 = LOG(t41)
      t43 = t42 + 0.2e1_dp
      t44 = br_d1 * t3
      t45 = 0.1e1_dp / t19
      t46 = t8 * t45
      t49 = br_d2 * t24
      t52 = pi ** 2
      t53 = br_d3 * t52
      t54 = t26 ** 2
      t55 = t54 * rho
      t57 = 0.1e1_dp / t29 / t19
      t58 = t55 * t57
      t61 = t3 * t52
      t62 = br_d4 * t61
      t63 = t54 * t26
      t64 = t7 * t63
      t65 = t29 ** 2
      t66 = 0.1e1_dp / t65
      t67 = t64 * t66
      t71 = t2 * t52 * pi
      t72 = br_d5 * t71
      t73 = t54 ** 2
      t74 = t6 * t73
      t76 = 0.1e1_dp / t65 / t19
      t77 = t74 * t76
      t80 = br_d0 + 0.2e1_dp / 0.3e1_dp * t44 * t46 + 0.4e1_dp / 0.9e1_dp * t49 * t31  &
      + 0.8e1_dp / 0.27e2_dp * t53 * t58 + 0.16e2_dp / 0.81e2_dp * t62 * t67 + 0.32e2_dp  &
      / 0.243e3_dp * t72 * t77
      t81 = t43 * t80
      t82 = br_e1 * t3
      t85 = br_e2 * t24
      t88 = br_e3 * t52
      t91 = br_e4 * t61
      t94 = br_e5 * t71
      t97 = br_e0 + 0.2e1_dp / 0.3e1_dp * t82 * t46 + 0.4e1_dp / 0.9e1_dp * t85 * t31  &
      + 0.8e1_dp / 0.27e2_dp * t88 * t58 + 0.16e2_dp / 0.81e2_dp * t91 * t67 + 0.32e2_dp  &
      / 0.243e3_dp * t94 * t77
      t98 = 0.1e1_dp / t97
      t99 = t81 * t98
      t100 = 8 ** (0.1e1_dp / 0.3e1_dp)
      t101 = t43 ** 2
      t102 = t101 * t43
      t103 = t80 ** 2
      t104 = t103 * t80
      t105 = t102 * t104
      t106 = t97 ** 2
      t108 = 0.1e1_dp / t106 / t97
      t109 = t105 * t108
      t110 = EXP(-t99)
      t111 = 0.1e1_dp / 0.3141592654e1_dp
      t112 = t110 * t111
      t114 = t109 * t112 * t13
      t115 = t114 ** (0.1e1_dp / 0.3e1_dp)
      t116 = 0.1e1_dp / t115
      t117 = REAL(t100,KIND=dp) * t116
      t118 = t117 * R
      t119 = t99 * t118
      t121 = EXP(0.2e1_dp * t119)
      t123 = t121 * R
      t124 = t123 * t43
      t125 = t80 * t98
      t126 = t125 * t117
      t128 = t121 * t43
      t130 = t99 + t119
      t131 = EXP(t130)
      t133 = 0.2e1_dp * t121 - t124 * t126 + t128 * t125 + 0.2e1_dp + t99 + t119  &
      - 0.4e1_dp * t131
      t134 = rho * t133
      t136 = EXP(-t130)
      t138 = t136 * REAL(t100,KIND=dp) * t116
      e_0 = e_0 + (t134 * t138 / 0.8e1_dp) * sx
    END IF
    IF(grad_deriv >= 1 .OR. grad_deriv == -1 ) THEN
      t144 = 0.1e1_dp / t7 / t26
      t151 = 0.1e1_dp / t7 / t27 * gamma * t12
      t154 = 0.1e1_dp / t35
      t157 = t6 * t26
      t158 = t157 * t30
      t161 = t6 * rho
      t163 = t57 * gamma
      t164 = t163 * t12
      t175 = t36 * t4
      t178 = -0.2500000000e1_dp * t5 * t144 * t19 - 0.1250000000e0_dp * t5 * t151  &
      + 0.3e1_dp / 0.4e1_dp * t154 * t1 * t4 * t20 * (0.40e2_dp / 0.27e2_dp * t25  &
      * t158 + 0.2e1_dp / 0.27e2_dp * t25 * t161 * t164) - 0.5e1_dp / 0.2e1_dp * t36  &
      * t4 * t144 * t19 - t175 * t151 / 0.8e1_dp
      t179 = 0.1e1_dp / t41
      t180 = t178 * t179
      t182 = t98 * REAL(t100,KIND=dp)
      t184 = t182 * t116 * R
      t185 = t180 * t80 * t184
      t187 = t7 * t45
      t190 = 0.1e1_dp / t6
      t192 = t30 * gamma
      t193 = t192 * t12
      t201 = t54 * t57
      t205 = t66 * gamma
      t206 = t205 * t12
      t209 = t7 * t55
      t210 = t209 * t66
      t213 = t7 * t54
      t215 = t76 * gamma
      t216 = t215 * t12
      t220 = t6 * t54 * t27
      t221 = t220 * t76
      t224 = t6 * t63
      t227 = 0.1e1_dp / t65 / t29
      t228 = t227 * gamma
      t229 = t228 * t12
      t232 = 0.10e2_dp / 0.9e1_dp * t44 * t187 + t44 * t190 * t193 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t49 * t158 + 0.2e1_dp / 0.27e2_dp * t49 * t161 * t164 +  &
      0.40e2_dp / 0.27e2_dp * t53 * t201 + 0.2e1_dp / 0.27e2_dp * t53 * t27 * t206  &
      + 0.320e3_dp / 0.243e3_dp * t62 * t210 + 0.16e2_dp / 0.243e3_dp * t62 * t213 *  &
      t216 + 0.800e3_dp / 0.729e3_dp * t72 * t221 + 0.40e2_dp / 0.729e3_dp * t72 *  &
      t224 * t229
      t234 = t43 * t232 * t98
      t235 = t234 * t118
      t237 = 0.1e1_dp / t106
      t238 = t81 * t237
      t264 = 0.10e2_dp / 0.9e1_dp * t82 * t187 + t82 * t190 * t193 / 0.18e2_dp + 0.40e2_dp  &
      / 0.27e2_dp * t85 * t158 + 0.2e1_dp / 0.27e2_dp * t85 * t161 * t164 +  &
      0.40e2_dp / 0.27e2_dp * t88 * t201 + 0.2e1_dp / 0.27e2_dp * t88 * t27 * t206  &
      + 0.320e3_dp / 0.243e3_dp * t91 * t210 + 0.16e2_dp / 0.243e3_dp * t91 * t213 *  &
      t216 + 0.800e3_dp / 0.729e3_dp * t94 * t221 + 0.40e2_dp / 0.729e3_dp * t94 *  &
      t224 * t229
      t267 = t238 * t117 * R * t264
      t270 = 0.1e1_dp / t115 / t114
      t271 = REAL(t100,KIND=dp) * t270
      t274 = t101 * t104 * t108 * t110
      t275 = t111 * t13
      t280 = t102 * t103 * t108
      t285 = t106 ** 2
      t287 = t105 / t285
      t292 = t180 * t125
      t294 = t81 * t237 * t264
      t302 = 0.3e1_dp * t274 * t275 * t180 + 0.3e1_dp * t280 * t112 * t13 * t232  &
      - 0.3e1_dp * t287 * t112 * t13 * t264 + t109 * (-t292 - t234 + t294)  &
      * t110 * t275 - t109 * t112 / t26
      t305 = t99 * t271 * R * t302
      t308 = (0.2e1_dp * t185 + 0.2e1_dp * t235 - 0.2e1_dp * t267 - 0.2e1_dp / 0.3e1_dp  &
      * t305) * t121
      t310 = R * t43
      t315 = t232 * t98
      t318 = t123 * t81
      t319 = t237 * REAL(t100,KIND=dp)
      t330 = t179 * t80 * t98
      t333 = t80 * t237
      t336 = t305 / 0.3e1_dp
      t337 = t292 + t234 - t294 + t185 + t235 - t267 - t336
      t340 = 0.2e1_dp * t308 - t308 * t310 * t126 - t123 * t180 * t126 - t124  &
      * t315 * t117 + t318 * t319 * t116 * t264 + t318 * t182 * t270 *  &
      t302 / 0.3e1_dp + t308 * t99 + t121 * t178 * t330 + t128 * t315 - t128  &
      * t333 * t264 + t292 + t234 - t294 + t185 + t235 - t267 - t336  &
      - 0.4e1_dp * t337 * t131
      t348 = t134 * t136
      e_rho = e_rho + (t133 * t136 * t117 / 0.8e1_dp + rho * t340 * t138 / 0.8e1_dp - t134  &
      * t337 * t138 / 0.8e1_dp - t348 * t271 * t302 / 0.24e2_dp ) * sx
      t353 = t144 * gamma * ndrho
      t356 = t154 * br_BB
      t357 = t356 * t3
      t365 = 0.2500000000000000e0_dp * t5 * t353 - t357 * t7 * t30 * gamma * ndrho  &
      / 0.9e1_dp + t175 * t353 / 0.4e1_dp
      t366 = t365 * t179
      t368 = t366 * t80 * t184
      t371 = t192 * ndrho
      t375 = t163 * ndrho
      t379 = t205 * ndrho
      t383 = t215 * ndrho
      t387 = t228 * ndrho
      t390 = -t44 * t7 * t371 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t49 * t157 * t375  &
      - 0.4e1_dp / 0.27e2_dp * t53 * t54 * t379 - 0.32e2_dp / 0.243e3_dp * t62 * t209  &
      * t383 - 0.80e2_dp / 0.729e3_dp * t72 * t220 * t387
      t392 = t43 * t390 * t98
      t393 = t392 * t118
      t410 = -t82 * t7 * t371 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t85 * t157 * t375  &
      - 0.4e1_dp / 0.27e2_dp * t88 * t54 * t379 - 0.32e2_dp / 0.243e3_dp * t91 * t209  &
      * t383 - 0.80e2_dp / 0.729e3_dp * t94 * t220 * t387
      t413 = t238 * t117 * R * t410
      t426 = t366 * t125
      t428 = t81 * t237 * t410
      t433 = 0.3e1_dp * t274 * t275 * t366 + 0.3e1_dp * t280 * t112 * t13 * t390  &
      - 0.3e1_dp * t287 * t112 * t13 * t410 + t109 * (-t426 - t392 + t428)  &
      * t110 * t275
      t436 = t99 * t271 * R * t433
      t439 = (0.2e1_dp * t368 + 0.2e1_dp * t393 - 0.2e1_dp * t413 - 0.2e1_dp / 0.3e1_dp  &
      * t436) * t121
      t445 = t390 * t98
      t461 = t436 / 0.3e1_dp
      t462 = t426 + t392 - t428 + t368 + t393 - t413 - t461
      t465 = 0.2e1_dp * t439 - t439 * t310 * t126 - t123 * t366 * t126 - t124  &
      * t445 * t117 + t318 * t319 * t116 * t410 + t318 * t182 * t270 *  &
      t433 / 0.3e1_dp + t439 * t99 + t121 * t365 * t330 + t128 * t445 - t128  &
      * t333 * t410 + t426 + t392 - t428 + t368 + t393 - t413 - t461  &
      - 0.4e1_dp * t462 * t131
      e_ndrho = e_ndrho + (rho * t465 * t138 / 0.8e1_dp - t134 * t462 * t138 / 0.8e1_dp - t348  &
      * t271 * t433 / 0.24e2_dp ) * sx
      t479 = t8 * t30
      t480 = t479 * gamma
      t485 = -0.1000000000e1_dp * t5 * t9 * gamma + 0.4e1_dp / 0.9e1_dp * t357 * t480  &
      - t36 * t37 * gamma
      t486 = t485 * t179
      t488 = t486 * t80 * t184
      t492 = t28 * t57
      t493 = t492 * gamma
      t496 = t55 * t66
      t497 = t496 * gamma
      t500 = t64 * t76
      t501 = t500 * gamma
      t504 = t74 * t227
      t505 = t504 * gamma
      t508 = 0.4e1_dp / 0.9e1_dp * t44 * t480 + 0.16e2_dp / 0.27e2_dp * t49 * t493 +  &
      0.16e2_dp / 0.27e2_dp * t53 * t497 + 0.128e3_dp / 0.243e3_dp * t62 * t501 + 0.320e3_dp  &
      / 0.729e3_dp * t72 * t505
      t510 = t43 * t508 * t98
      t511 = t510 * t118
      t523 = 0.4e1_dp / 0.9e1_dp * t82 * t480 + 0.16e2_dp / 0.27e2_dp * t85 * t493 +  &
      0.16e2_dp / 0.27e2_dp * t88 * t497 + 0.128e3_dp / 0.243e3_dp * t91 * t501 + 0.320e3_dp  &
      / 0.729e3_dp * t94 * t505
      t526 = t238 * t117 * R * t523
      t539 = t486 * t125
      t541 = t81 * t237 * t523
      t546 = 0.3e1_dp * t274 * t275 * t486 + 0.3e1_dp * t280 * t112 * t13 * t508  &
      - 0.3e1_dp * t287 * t112 * t13 * t523 + t109 * (-t539 - t510 + t541)  &
      * t110 * t275
      t549 = t99 * t271 * R * t546
      t552 = (0.2e1_dp * t488 + 0.2e1_dp * t511 - 0.2e1_dp * t526 - 0.2e1_dp / 0.3e1_dp  &
      * t549) * t121
      t558 = t508 * t98
      t574 = t549 / 0.3e1_dp
      t575 = t539 + t510 - t541 + t488 + t511 - t526 - t574
      t578 = 0.2e1_dp * t552 - t552 * t310 * t126 - t123 * t486 * t126 - t124  &
      * t558 * t117 + t318 * t319 * t116 * t523 + t318 * t182 * t270 *  &
      t546 / 0.3e1_dp + t552 * t99 + t121 * t485 * t330 + t128 * t558 - t128  &
      * t333 * t523 + t539 + t510 - t541 + t488 + t511 - t526 - t574  &
      - 0.4e1_dp * t575 * t131
      e_tau = e_tau + ( rho * t578 * t138 / 0.8e1_dp - t134 * t575 * t138 / 0.8e1_dp - t348  &
      * t271 * t546 / 0.24e2_dp ) * sx
      t597 = 0.2500000000000000e0_dp * t5 * t9 - t356 * t3 * t8 * t30 / 0.9e1_dp  &
      + t36 * t37 / 0.4e1_dp
      t598 = t597 * t179
      t600 = t598 * t80 * t184
      t612 = -t44 * t479 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t49 * t492 - 0.4e1_dp /  &
      0.27e2_dp * t53 * t496 - 0.32e2_dp / 0.243e3_dp * t62 * t500 - 0.80e2_dp / 0.729e3_dp  &
      * t72 * t504
      t614 = t43 * t612 * t98
      t615 = t614 * t118
      t627 = -t82 * t479 / 0.9e1_dp - 0.4e1_dp / 0.27e2_dp * t85 * t492 - 0.4e1_dp /  &
      0.27e2_dp * t88 * t496 - 0.32e2_dp / 0.243e3_dp * t91 * t500 - 0.80e2_dp / 0.729e3_dp  &
      * t94 * t504
      t630 = t238 * t117 * R * t627
      t643 = t598 * t125
      t645 = t81 * t237 * t627
      t650 = 0.3e1_dp * t274 * t275 * t598 + 0.3e1_dp * t280 * t112 * t13 * t612  &
      - 0.3e1_dp * t287 * t112 * t13 * t627 + t109 * (-t643 - t614 + t645)  &
      * t110 * t275
      t653 = t99 * t271 * R * t650
      t656 = (0.2e1_dp * t600 + 0.2e1_dp * t615 - 0.2e1_dp * t630 - 0.2e1_dp / 0.3e1_dp  &
      * t653) * t121
      t662 = t612 * t98
      t678 = t653 / 0.3e1_dp
      t679 = t643 + t614 - t645 + t600 + t615 - t630 - t678
      t682 = 0.2e1_dp * t656 - t656 * t310 * t126 - t123 * t598 * t126 - t124  &
      * t662 * t117 + t318 * t319 * t116 * t627 + t318 * t182 * t270 *  &
      t650 / 0.3e1_dp + t656 * t99 + t121 * t597 * t330 + t128 * t662 - t128  &
      * t333 * t627 + t643 + t614 - t645 + t600 + t615 - t630 - t678  &
      - 0.4e1_dp * t679 * t131
      e_laplace_rho = e_laplace_rho + (rho * t682 * t138 / 0.8e1_dp - t134 * t679 * t138 / 0.8e1_dp - t348  &
      * t271 * t650 / 0.24e2_dp) * sx
    END IF

  END SUBROUTINE x_br_lsd_y_gt_0_cutoff_R_lte_b

END MODULE xc_xbecke_roussel
