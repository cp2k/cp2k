!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the Exchange contribution in the BEEF-vdW functional
!         (Norskov Group)
!> \par History
!>      02.2014 created based on xc_xbecke88.F [rkoitz]
!> \author rkoitz
! *****************************************************************************
MODULE xc_xbeef

  USE bibliography,                    ONLY: Wellendorff2012,&
                                             cite_reference
  USE cp_array_utils_r,                ONLY: cp_3d_r_p_type
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xbeef'

  PUBLIC :: xbeef_lda_info, xbeef_lsd_info, xbeef_lda_eval, xbeef_lsd_eval

    !for beef functional pulled out of GPAW source, Feb 2014
  REAL (kind=dp), PARAMETER :: a(0:29) = (/ 1.516501714304992365356,  0.441353209874497942611, -0.091821352411060291887, &
                                           -0.023527543314744041314,  0.034188284548603550816,  0.002411870075717384172, &
                                           -0.014163813515916020766,  0.000697589558149178113,  0.009859205136982565273, &
                                           -0.006737855050935187551, -0.001573330824338589097,  0.005036146253345903309, &
                                           -0.002569472452841069059, -0.000987495397608761146,  0.002033722894696920677, &
                                           -0.000801871884834044583, -0.000668807872347525591,  0.001030936331268264214, &
                                           -0.000367383865990214423, -0.000421363539352619543,  0.000576160799160517858, &
                                           -0.000083465037349510408, -0.000445844758523195788,  0.000460129009232047457, &
                                           -0.000005231775398304339, -0.000423957047149510404,  0.000375019067938866537, &
                                            0.000021149381251344578, -0.000190491156503997170,  0.000073843624209823442 /)
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      02.2014 created 
!> \author rkoitz
! *****************************************************************************

  SUBROUTINE xbeef_lda_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'xbeef_lda_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "Wellendorff, J. et al., Phys. Rev. B 85, 235149 (2012) {LDA}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Exchange Contribution to BEEF-vdW Functional (Wellendorff, 2012) {LDA}"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%rho_1_3=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xbeef_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      02.2014 created 
!> \author rkoitz
! *****************************************************************************
  SUBROUTINE xbeef_lsd_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'xbeef_lsd_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "Wellendorff, J. et al., Phys. Rev. B 85, 235149 (2012) {LSD}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Exchange Contribution to BEEF-vdW Functional (Wellendorff, 2012) {LSD}"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%rho_spin_1_3=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xbeef_lsd_info

! *****************************************************************************
!> \brief evaluates the beef exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param xbeef_params input parameters (scaling)
!> \par History
!>      02.2014 created
!> \author rkoitz
! *****************************************************************************
  SUBROUTINE xbeef_lda_eval(rho_set,deriv_set,grad_deriv,xbeef_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xbeef_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xbeef_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrho, e_rho, &
                                                norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(xbeef_params,"scale_x",r_val=sx)

    CALL cite_reference(Wellendorff2012)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
         norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho)
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
            message="derivatives greater than 1 not implemented", &
            error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp          shared(rho, rho_1_3, norm_drho, e_0, e_rho) &
    !$omp          shared(e_ndrho) &
    !$omp          shared( grad_deriv, npoints) &
    !$omp          shared(epsilon_rho,sx)
    CALL xbeef_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
         e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
         grad_deriv=grad_deriv,&
         npoints=npoints,epsilon_rho=epsilon_rho,sx=sx)
    !$omp end parallel

    IF (cp_debug) THEN
       DEALLOCATE(dummy)
    ELSE
       NULLIFY(dummy)
    END IF
    CALL timestop(handle)
  END SUBROUTINE xbeef_lda_eval

! *****************************************************************************
!> \brief evaluates the beef exchange functional for lda
!> \param rho the density where you want to evaluate the functional
!> \param rho_1_3 ...
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sx scaling-parameter for exchange
!> \par History 
!>      02.2014 created based on xc_xbecke88
!> \author rkoitz
! *****************************************************************************
SUBROUTINE xbeef_lda_calc(rho, rho_1_3, norm_drho,&
       e_0,e_rho,e_ndrho,&
       grad_deriv,npoints,epsilon_rho,sx)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho_1_3, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx

    CHARACTER(len=*), PARAMETER :: routineN = 'xbeef_lda_calc', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 30

    INTEGER                                  :: i, ii
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ds_ndrho, ds_rho, dt, e_ueg, &
                                                e_ueg_drho, epsilon_rho43, &
                                                k_lda, kf, my_rho, &
                                                my_rho_1_3, s, s2, t, t3
    REAL(kind=dp), DIMENSION(0:29)           :: de_leg, e_leg

!energies, first and second derivatives from legendre pol

    failure=.FALSE.

    kf = (3.0_dp * pi**2)**(1.0_dp / 3.0_dp) !only constants, without n^1/3
    k_lda = -(2.0_dp / (2.0_dp**(4._dp/3._dp))) * (3.0_dp/2.0_dp) * (3.0_dp / (4.0_dp * pi))**(1.0_dp/3.0_dp) 
        !exchange energy density of the uniform electron gas, non spin-polarized

    epsilon_rho43=epsilon_rho**(4._dp/3._dp)

    !$omp do
    DO ii=1,npoints

      my_rho = rho(ii)

      IF (my_rho>epsilon_rho) THEN
         my_rho_1_3 = rho_1_3(ii)

         e_ueg =  k_lda * my_rho * my_rho_1_3 
         e_ueg_drho = (4.0_dp/3.0_dp) * k_lda * my_rho_1_3 

         t3 = my_rho_1_3 * my_rho * 2 * kf        !reduced gradient, denominator

         s = norm_drho(ii)/MAX(t3,epsilon_rho43)  !reduced gradient finally
         s2 = s**2
         t = 2.0_dp*s2 / (4.0_dp + s2) - 1.0_dp             

         IF (grad_deriv>=0) THEN       !asking for pure e evaluation or also derivatives
             e_leg(0) = 1              !first legendre pol
             e_leg(1) = t              !second legendre pol
         END IF
         
         IF ((grad_deriv>=1) .OR. (grad_deriv==-1)) THEN  !asking for first derivative or higher
             de_leg(0) = 0       
             de_leg(1) = 1       
             dt = 4.0_dp*s / (4.0_dp + s2) - 4.0_dp * s * s2 / (4.0_dp + s2)**2
             ds_rho = -(4.0_dp * s) / (3.0_dp * MAX(my_rho,epsilon_rho))
             ds_ndrho = 1.0_dp/(MAX(t3,epsilon_rho43))
         END IF
         
         DO i=2,m-1    !LEGENDRE PART
            e_leg(i) = 2.*(t)*e_leg(i-1) - e_leg(i-2) - ((t)*e_leg(i-1) - e_leg(i-2))/(FLOAT(i))
             !taken from quantum espresso beef library.
         
             IF (ABS(grad_deriv)>=1) THEN   !first derivative 
                 !the zero-derivatives need to be available for the first deriv.
                 de_leg(i) = e_leg(i-1)*i + de_leg(i-1)*(t) 
             END IF
         END DO

        !NO DERIVATIVE
         IF (grad_deriv >= 0) THEN
            !add the scaled legendre linear combination to e_0
            e_0(ii) = e_0(ii) + SUM(e_leg*a) * e_ueg  * sx
        END IF

        !FIRST DERIVATIVE
        IF ((grad_deriv>=1) .OR. (grad_deriv==-1)) THEN  !asking for first derivative or higher
            e_rho(ii)   = e_rho(ii)   + (SUM(e_leg*a) * e_ueg_drho + SUM(de_leg*a) * dt * ds_rho * e_ueg ) * sx 
            e_ndrho(ii) = e_ndrho(ii) + (SUM(de_leg*a) * dt * ds_ndrho * e_ueg  ) * sx 
        END IF

      END IF
   END DO

   !$omp end do

  END SUBROUTINE xbeef_lda_calc

! *****************************************************************************
!> \brief evaluates the beef 88 exchange functional for lsd
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param xbeef_params ...
!> \par History
!>         2/2014 rkoitz [created based on Becke 88]
!> \author rkoitz
! *****************************************************************************
  SUBROUTINE xbeef_lsd_eval(rho_set,deriv_set,grad_deriv,xbeef_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xbeef_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xbeef_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0
    TYPE(cp_3d_r_p_type), DIMENSION(2)       :: e_ndrho, e_rho, norm_drho, &
                                                rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    CALL cite_reference(Wellendorff2012)

    failure=.FALSE.
    NULLIFY(deriv, bo)
    DO i=1,2
       NULLIFY(norm_drho(i)%array, rho(i)%array, rho_1_3(i)%array)
    END DO

    CALL section_vals_val_get(xbeef_params,"scale_x",r_val=sx)
    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    CALL xc_rho_set_get(rho_set,rhoa_1_3=rho_1_3(1)%array,&
         rhob_1_3=rho_1_3(2)%array,rhoa=rho(1)%array,&
         rhob=rho(2)%array,norm_drhoa=norm_drho(1)%array, &
         norm_drhob=norm_drho(2)%array,rho_cutoff=epsilon_rho,&
         local_bounds=bo)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
       ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
    ELSE
       dummy=> rho(1)%array
    END IF

    e_0=>dummy
    DO i=1,2
       e_rho(i)%array => dummy
       e_ndrho(i)%array => dummy
    END DO

   IF (grad_deriv>=0) THEN
       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv, deriv_data=e_0)
    END IF
    IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho(2)%array)
    END IF
    IF (grad_deriv>1.OR.grad_deriv<-1) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives greater than 1 not implemented", &
            error_level=cp_failure_level)
    END IF

    DO ispin=1,2

       !$omp parallel default(none) &
       !$omp          shared(rho, ispin, rho_1_3, norm_drho, e_0) &
       !$omp          shared(e_rho, e_ndrho) &
       !$omp          shared(grad_deriv, npoints) &
       !$omp          shared(epsilon_rho, sx)

       CALL xbeef_lsd_calc(&
            rho_spin=rho(ispin)%array,&
            rho_1_3_spin=rho_1_3(ispin)%array,&
            norm_drho_spin=norm_drho(ispin)%array,&
            e_0=e_0,e_rho_spin=e_rho(ispin)%array,&
            e_ndrho_spin=e_ndrho(ispin)%array,&
            grad_deriv=grad_deriv, npoints=npoints,&
            epsilon_rho=epsilon_rho,sx=sx)

       !$omp end parallel

    END DO

    IF (cp_debug) THEN
       DEALLOCATE(dummy)
    ELSE
       NULLIFY(dummy)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xbeef_lsd_eval
! *****************************************************************************
!> \brief low level calculation of the beef exchange functional for lsd
!> \param rho_spin alpha or beta spin density
!> \param rho_1_3_spin rho_spin**(1./3.)
!> \param norm_drho_spin || grad rho_spin ||
!> \param e_0 adds to it the local value of the functional
!> \param e_rho_spin e_*_spin: derivative of the functional wrt. to the variables
!>        named where the * is. Everything wrt. to the spin of the arguments.
!> \param e_ndrho_spin ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sx scaling-parameter for exchange
!> \par History
!>      02.2014 created based on Becke88
!> \author rkoitz
! *****************************************************************************
SUBROUTINE xbeef_lsd_calc(rho_spin, rho_1_3_spin, norm_drho_spin,e_0,&
       e_rho_spin,e_ndrho_spin,grad_deriv,npoints,epsilon_rho,sx)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho_spin, rho_1_3_spin, &
                                                norm_drho_spin
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rho_spin, e_ndrho_spin
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx

    CHARACTER(len=*), PARAMETER :: routineN = 'xbeef_lsd_calc', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 30

    INTEGER                                  :: i, ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: ds_ndrho, ds_rho, dt, e_ueg, e_ueg_drho, epsilon_rho43, &
      k_lsd, kf, my_epsilon_rho, my_rho, my_rho_1_3, s, s2, t, t3
    REAL(kind=dp), DIMENSION(0:29)           :: de_leg, e_leg

!energies and first derivatives from legendre pol

    failure=.FALSE.
    kf = (3.0_dp * pi**2)**(1.0_dp / 3.0_dp) !only constants, without n^1/3
    k_lsd = (3.0_dp/2.0_dp) * (3.0_dp / (4.0_dp * pi))**(1.0_dp/3.0_dp) 
        !exchange energy density of the uniform electron gas, spin-polarized

    my_epsilon_rho=0.5_dp*epsilon_rho
    epsilon_rho43=my_epsilon_rho**(4._dp/3._dp)


    !$omp do
    DO ii=1,npoints
      my_rho = rho_spin(ii)

      IF (my_rho>epsilon_rho) THEN
         my_rho_1_3 = rho_1_3_spin(ii)

         e_ueg =  k_lsd * my_rho * my_rho_1_3 
         e_ueg_drho = (4.0_dp/3.0_dp) * k_lsd * my_rho_1_3

         t3 = my_rho_1_3 * my_rho * 2 * kf        !reduced gradient, denominator

         s = norm_drho_spin(ii)/MAX(t3,epsilon_rho43)  !reduced gradient finally
         s2 = s**2
         t = 2.0_dp*s**2 / (4.0_dp + s**2) - 1.0_dp             

         IF (grad_deriv>=0) THEN       !asking for pure e evaluation or also derivatives
             e_leg(0) = 1              !first legendre pol
             e_leg(1) = t              !second legendre pol
         END IF
         
         IF ((grad_deriv>=1) .OR. (grad_deriv==-1)) THEN  !asking for first derivative or higher
             de_leg(0) = 0       
             de_leg(1) = 1       
             dt = 4.0_dp*s / (4.0_dp + s2) - 4.0_dp * s * s2 / (4.0_dp + s2)**2
             ds_rho = -(4.0_dp * s) / (3.0_dp * MAX(my_rho,epsilon_rho))
             ds_ndrho = 1.0_dp/(MAX(t3,epsilon_rho43))
         END IF
         
         DO i=2,m-1    !LEGENDRE PART
             e_leg(i) = 2.*(t)*e_leg(i-1) - e_leg(i-2) - ((t)*e_leg(i-1) - e_leg(i-2))/(FLOAT(i))
             !taken from quantum espresso beef library.

             IF (ABS(grad_deriv)>=1) THEN   !first derivative
                 !the zero-derivatives need to be available for the first deriv.
                de_leg(i) = e_leg(i-1)*i + de_leg(i-1)*(t) 
             END IF
         
         END DO

        !NO DERIVATIVE
         IF (grad_deriv >= 0) THEN
            !add the scaled legendre linear combination to e_0
            e_0(ii) = e_0(ii) + SUM(e_leg*a) * e_ueg  * sx
        END IF


        !FIRST DERIVATIVE
        IF ((grad_deriv>=1) .OR. (grad_deriv==-1)) THEN  !asking for first derivative or higher
            e_rho_spin(ii)   = e_rho_spin(ii)   + (SUM(e_leg*a) * e_ueg_drho + SUM(de_leg*a) * dt * ds_rho * e_ueg ) * sx 
            e_ndrho_spin(ii) = e_ndrho_spin(ii) + (SUM(de_leg*a) * dt * ds_ndrho * e_ueg ) * sx 
        END IF
      END IF
    END DO
  !$omp end do

  END SUBROUTINE xbeef_lsd_calc

END MODULE xc_xbeef

