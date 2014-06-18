!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the b97 correlation functional
!> \note
!>      This was generated with the help of a maple worksheet that you can
!>      find under doc/b97.mw .
!>      I did not add 3. derivatives in the polazied (lsd) case because it
!>      would have added lots of code. If you need them the worksheet
!>      is already prepared for them, and by uncommenting a couple of lines
!>      you should be able to generate them.
!>      Other parametrizations (B97-1 by FA Hamprecht, AJ Cohen, DJ Tozer, NC Handy)
!>      could be easily added, the maple file should be straightforward to extend
!>      to HCTH (and thus implement it also for unrestricted calculations).
!> \par History
!>      01.2009 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE xc_b97
  USE bibliography,                    ONLY: Becke1997,&
                                             Grimme2006,&
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
  USE xc_input_constants,              ONLY: xc_b97_grimme,&
                                             xc_b97_mardirossian,&
                                             xc_b97_orig
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_b97'

  PUBLIC :: b97_lda_info, b97_lsd_info, b97_lda_eval, b97_lsd_eval

  REAL(dp), DIMENSION(10) :: params_b97_orig=(/ 0.8094_dp, 0.5073_dp, 0.7481_dp, &
    0.9454_dp, 0.7471_dp, -4.5961_dp, 0.1737_dp, 2.3487_dp, -2.4868_dp , 1.0_dp-0.1943_dp /)
  REAL(dp), DIMENSION(10) :: params_b97_grimme=(/ 1.08662_dp, -0.52127_dp, 3.25429_dp, &
    0.69041_dp, 6.30270_dp, -14.9712_dp, 0.22340_dp, -1.56208_dp, 1.94293_dp, 1.0_dp /)
  REAL(dp), DIMENSION(10) :: params_b97_mardirossian=(/ 0.833_dp, 0.603_dp, 1.194_dp, &
    1.219_dp, -1.850_dp, 0.00_dp, 0.556_dp, -0.257_dp, 0.00_dp, 1.0_dp /)

CONTAINS

  SUBROUTINE b97_ref(param,lda,sx,sc,reference,shortform,error)
    INTEGER                                  :: param
    LOGICAL                                  :: lda
    REAL(dp), INTENT(in)                     :: sx, sc
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_ref', &
      routineP = moduleN//':'//routineN

    CHARACTER(20)                            :: pol

    IF (lda) THEN
       pol="(unpolarized)"
    ELSE
       pol="(polarized)"
    END IF
    SELECT CASE(param)
    CASE(xc_b97_orig)
       CALL cite_reference(Becke1997)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference = "A.D.Becke, "//&
                  "J. Chem. Phys, vol. 107, pp. 8554, (1997), needs exact x, "//&
                  pol
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "B97, Becke 1997 xc functional, needs exact x "//pol
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
             WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a,a)")&
                  "A.D.Becke, ",&
                  "J. Chem. Phys, vol. 107, pp. 8554, (1997)",&
                  sx,sc,", needs exact x ",pol
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,a,'sx=',f5.3,'sc=',f5.3)")&
                  "B97, Becke 1997 xc functional (unpolarized)",pol,sx,sc
          END IF
       END IF
    CASE(xc_b97_grimme)
       CALL cite_reference(Becke1997)
       CALL cite_reference(Grimme2006)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference = "B97-D, Grimme xc functional,"//&
                  " J Comput Chem 27: 1787-1799 (2006)"//&
                  ", needs C6 dispersion "//pol
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "B97-D, Grimme xc functional "//pol
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
             WRITE (reference,"(a,a,3x,' sx=',f6.3,' sc=',f6.3,1x,a)")&
                  "B97-D, Grimme xc functional,",&
                  " J Comput Chem 27: 1787-1799 (2006) ",&
                  sx,sc,pol
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,a,3x,' sx=',f6.3,' sc=',f6.3,' (LDA)')")&
                  "B97-D, Grimme xc functional ",pol,sx,sc
          END IF
       END IF
    CASE(xc_b97_mardirossian)
       CALL cite_reference(Becke1997)
!       CALL cite_reference(Mardirossian2014)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference = "wB97X-V, xc functional,"//&
                  " Mardirossian and Head-Gordon, PCCP DOI: 10.1039/c3cp54374a"//&
                  ", needs HFX exchange and VV10 dispersion (NOT TESTED)"//pol
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "wB97X-V, HFX+B97+VV10 functional (NOT TESTED)"//pol
          END IF
       ELSE
         CPAssertNoFail(.FALSE.,cp_failure_level,routineP,error)
       END IF
    CASE default
       CPAssertNoFail(.FALSE.,cp_failure_level,routineP,error)
    END SELECT
  END SUBROUTINE b97_ref

  ! *****************************************************************************
  !> \brief return various information on the functional
  !> \param b97_params section selecting the varius parameters for the functional
  !> \param reference string with the reference of the actual functional
  !> \param shortform string with the shortform of the functional name
  !> \param needs the components needed by this functional are set to
  !>        true (does not set the unneeded components to false)
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! *****************************************************************************
  SUBROUTINE b97_lda_info(b97_params,reference,shortform, needs, max_deriv,&
       error)
    TYPE(section_vals_type), POINTER         :: b97_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_lda_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: param
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: sc, sx

    failure=.FALSE.
    CALL section_vals_val_get(b97_params,"parametrization",i_val=param,error=error)
    CALL section_vals_val_get(b97_params,"scale_x",r_val=sx,error=error)
    CALL section_vals_val_get(b97_params,"scale_c",r_val=sc,error=error)

    CALL b97_ref(param,.TRUE.,sx,sc,reference,shortform,error)
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=2

  END SUBROUTINE b97_lda_info

  ! *****************************************************************************
  !> \brief return various information on the functional
  !> \param b97_params section selecting the varius parameters for the functional
  !> \param reference string with the reference of the actual functional
  !> \param shortform string with the shortform of the functional name
  !> \param needs the components needed by this functional are set to
  !>        true (does not set the unneeded components to false)
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! *****************************************************************************
  SUBROUTINE b97_lsd_info(b97_params,reference,shortform, needs, max_deriv,&
       error)
    TYPE(section_vals_type), POINTER         :: b97_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_lsd_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: param
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: sc, sx

    failure=.FALSE.
    CALL section_vals_val_get(b97_params,"parametrization",i_val=param,error=error)
    CALL section_vals_val_get(b97_params,"scale_x",r_val=sx,error=error)
    CALL section_vals_val_get(b97_params,"scale_c",r_val=sc,error=error)

    CALL b97_ref(param,.FALSE.,sx,sc,reference,shortform,error)
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=2

  END SUBROUTINE b97_lsd_info

  ! *****************************************************************************
  !> \brief evaluates the b97 correlation functional for lda
  !> \param rho_set the density where you want to evaluate the functional
  !> \param deriv_set place where to store the functional derivatives (they are
  !>        added to the derivatives)
  !> \param grad_deriv degree of the derivative that should be evalated,
  !>        if positive all the derivatives up to the given degree are evaluated,
  !>        if negative only the given degree is calculated
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! *****************************************************************************
  SUBROUTINE b97_lda_eval(rho_set,deriv_set,grad_deriv,b97_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: b97_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, param, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, scale_c, scale_x
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, &
      e_ndrho_rho_rho, e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
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
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)
       END IF
       IF (grad_deriv>3.OR.grad_deriv<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL section_vals_val_get(b97_params,"parametrization",i_val=param,error=error)
       CALL section_vals_val_get(b97_params,"scale_c",r_val=scale_c,error=error)
       CALL section_vals_val_get(b97_params,"scale_x",r_val=scale_x,error=error)

       !$omp parallel default(none) &
       !$omp          shared(rho, norm_drho, e_0, e_rho, e_ndrho, e_rho_rho) &
       !$omp          shared(e_ndrho_rho, e_ndrho_ndrho) &
       !$omp          shared(grad_deriv, npoints, epsilon_rho) &
       !$omp          shared(epsilon_norm_drho, param, scale_c, scale_x, error)

       CALL b97_lda_calc(rho_tot=rho, norm_drho=norm_drho,&
            e_0=e_0,e_r=e_rho,e_ndr=e_ndrho,e_r_r=e_rho_rho,&
            e_r_ndr=e_ndrho_rho, e_ndr_ndr=e_ndrho_ndrho, &
           ! e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
           ! e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,&
            grad_deriv=grad_deriv,&
            npoints=npoints,epsilon_rho=epsilon_rho,epsilon_drho=epsilon_norm_drho,&
            param=param,scale_c_in=scale_c,scale_x_in=scale_x,error=error)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE b97_lda_eval

  ! *****************************************************************************
  !> \brief evaluates the b 97 xc functional for lsd
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! *****************************************************************************
  SUBROUTINE b97_lsd_eval(rho_set,deriv_set,grad_deriv,b97_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: b97_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, param, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho, &
                                                scale_c, scale_x
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndra, &
      e_ndra_ndra, e_ndra_ndrb, e_ndra_ra, e_ndra_rb, e_ndrb, e_ndrb_ndrb, &
      e_ndrb_ra, e_ndrb_rb, e_ra, e_ra_ra, e_ra_rb, e_rb, e_rb_rb, &
      norm_drhoa, norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(deriv, bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,&
            rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
            norm_drhob=norm_drhob, &
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
       e_ndra => dummy
       e_ndrb => dummy

       e_ndra_ra => dummy
       e_ndra_rb => dummy
       e_ndrb_rb => dummy
       e_ndrb_ra => dummy

       e_ndra_ndra => dummy
       e_ndrb_ndrb => dummy
       e_ndra_ndrb => dummy

       e_ra_ra => dummy
       e_ra_rb => dummy
       e_rb_rb => dummy

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
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndra_ra,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndra_rb,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrb_rb,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrb_ra,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndra_ndra,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhoa)(norm_drhob)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndra_ndrb,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrb_ndrb,error=error)
       END IF
       IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
          ! to do
       END IF

       CALL section_vals_val_get(b97_params,"parametrization",i_val=param,error=error)
       CALL section_vals_val_get(b97_params,"scale_x",r_val=scale_x,error=error)
       CALL section_vals_val_get(b97_params,"scale_c",r_val=scale_c,error=error)

       !$omp parallel default (none) &
       !$omp          shared(rhoa, rhob, norm_drhoa, norm_drhob, e_0, e_ra) &
       !$omp          shared(e_rb, e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb) &
       !$omp          shared(e_ndra_ra, e_ndrb_ra, e_ndrb_rb, e_ndra_rb) &
       !$omp          shared(e_ndra_ndra, e_ndrb_ndrb, e_ndra_ndrb) &
       !$omp          shared(grad_deriv, npoints, param, scale_c, scale_x) &
       !$omp          shared(epsilon_rho, epsilon_drho, error)

       CALL b97_lsd_calc(&
            rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa,&
            norm_drhob=norm_drhob, e_0=e_0, &
            e_ra=e_ra, e_rb=e_rb, &
            e_ndra=e_ndra, e_ndrb=e_ndrb, &
            e_ra_ra=e_ra_ra, e_ra_rb=e_ra_rb, e_rb_rb=e_rb_rb,&
            e_ra_ndra=e_ndra_ra, e_ra_ndrb=e_ndrb_ra, &
            e_rb_ndrb=e_ndrb_rb, e_rb_ndra=e_ndra_rb,&
            e_ndra_ndra=e_ndra_ndra, e_ndrb_ndrb=e_ndrb_ndrb,&
            e_ndra_ndrb=e_ndra_ndrb,&
            grad_deriv=grad_deriv, npoints=npoints, &
            epsilon_rho=epsilon_rho,epsilon_drho=epsilon_drho,&
            param=param,scale_c_in=scale_c,scale_x_in=scale_x,error=error)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE b97_lsd_eval

  FUNCTION b97_coeffs(param, error) RESULT(res)
    INTEGER, INTENT(in)                      :: param
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(dp), DIMENSION(10)                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_coeffs', &
      routineP = moduleN//':'//routineN

    SELECT CASE(param)
    CASE(xc_b97_orig)
       res = params_b97_orig
    CASE(xc_b97_grimme)
       res = params_b97_grimme
    CASE(xc_b97_mardirossian)
       res = params_b97_mardirossian
    CASE default
       CPAssertNoFail(.FALSE.,cp_failure_level,routineP,error)
       res=0.0_dp
    END SELECT
  END FUNCTION b97_coeffs

  ! *****************************************************************************
  !> \brief low level calculation of the b97 exchange-correlation functional for lsd
  !> \param rhoa ,rhob: alpha or beta spin density
  !> \param norm_drhoa ,norm_drhob,norm_drho: || grad rhoa |||,| grad rhoa ||,
  !>        || grad (rhoa+rhob) ||
  !> \param e_ 0: adds to it the local value of the functional
  !> \param e_ *: derivative of the functional wrt. to the variables
  !>        named where the * is.
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! *****************************************************************************
  SUBROUTINE b97_lsd_calc(rhoa, rhob, norm_drhoa, norm_drhob,&
       e_0, e_ra, e_rb, e_ndra, e_ndrb, &
       e_ra_ndra,e_ra_ndrb, e_rb_ndra, e_rb_ndrb,&
       e_ndra_ndra, e_ndrb_ndrb, e_ndra_ndrb, &
       e_ra_ra, e_ra_rb, e_rb_rb,&
       grad_deriv,npoints,epsilon_rho,epsilon_drho, &
       param, scale_c_in, scale_x_in, error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drhoa, &
                                                norm_drhob
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_ra, e_rb, e_ndra, &
      e_ndrb, e_ra_ndra, e_ra_ndrb, e_rb_ndra, e_rb_ndrb, e_ndra_ndra, &
      e_ndrb_ndrb, e_ndra_ndrb, e_ra_ra, e_ra_rb, e_rb_rb
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(in)                      :: param
    REAL(kind=dp), INTENT(in)                :: scale_c_in, scale_x_in
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_lsd_calc', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: small = 1.0e-20

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: A_1, A_2, A_3, alpha_1_1, alpha_1_2, alpha_1_3, alpha_c, &
      alpha_c1rhoa, alpha_c1rhob, alpha_crhoa, alpha_crhob, beta_1_1, &
      beta_1_2, beta_1_3, beta_2_1, beta_2_2, beta_2_3, beta_3_1, beta_3_2, &
      beta_3_3, beta_4_1, beta_4_2, beta_4_3, c_cab_0, c_cab_1, c_cab_2, &
      c_css_0, c_css_1, c_css_2, c_x_0, c_x_1, c_x_2, chi, chirhoa, &
      chirhoarhoa, chirhoarhob, chirhob, chirhobrhob, e_c_u_0, e_c_u_01rhoa, &
      e_c_u_01rhob, e_c_u_0rhoa, e_c_u_0rhoarhoa, e_c_u_0rhoarhob, &
      e_c_u_0rhob, e_c_u_0rhobrhob, e_c_u_1rhoa, e_c_u_1rhob, e_lsda_c_a, &
      e_lsda_c_a1rhoa, e_lsda_c_ab, e_lsda_c_abrhoa
    REAL(kind=dp) :: e_lsda_c_abrhob, e_lsda_c_arhoa, e_lsda_c_arhoarhoa, &
      e_lsda_c_b, e_lsda_c_b1rhob, e_lsda_c_brhob, e_lsda_c_brhobrhob, &
      e_lsda_x_a, e_lsda_x_arhoa, e_lsda_x_b, e_lsda_x_brhob, epsilon_c_unif, &
      epsilon_c_unif1rhoa, epsilon_c_unif1rhob, epsilon_c_unif_a, &
      epsilon_c_unif_a1rhoa, epsilon_c_unif_arhoa, epsilon_c_unif_b, &
      epsilon_c_unif_b1rhob, epsilon_c_unif_brhob, epsilon_c_unifrhoa, &
      epsilon_c_unifrhob, exc, exc_norm_drhoa, exc_norm_drhoa_norm_drhoa, &
      exc_norm_drhoa_norm_drhob, exc_norm_drhob, exc_norm_drhob_norm_drhob, &
      exc_rhoa, exc_rhoa_norm_drhoa, exc_rhoa_norm_drhob, exc_rhoa_rhoa
    REAL(kind=dp) :: exc_rhoa_rhob, exc_rhob, exc_rhob_norm_drhoa, &
      exc_rhob_norm_drhob, exc_rhob_rhob, f, f1rhoa, f1rhob, frhoa, &
      frhoarhoa, frhoarhob, frhob, frhobrhob, gamma_c_ab, gamma_c_ss, &
      gamma_x, gc_a, gc_ab, gc_abnorm_drhoa, gc_abnorm_drhob, gc_abrhoa, &
      gc_abrhob, gc_anorm_drhoa, gc_arhoa, gc_b, gc_bnorm_drhob, gc_brhob, &
      gx_a, gx_anorm_drhoa, gx_arhoa, gx_b, gx_bnorm_drhob, gx_brhob, &
      my_norm_drhoa, my_norm_drhob, my_rhoa, my_rhob, p_1, p_2, p_3, rho, rs, &
      rs_a, rs_arhoa, rs_arhoarhoa, rs_b, rs_brhob, rs_brhobrhob, rsrhoa, &
      rsrhoarhoa, rsrhoarhob, rsrhob, rsrhobrhob, s_a, s_a_2
    REAL(kind=dp) :: s_a_21norm_drhoa, s_a_21rhoa, s_a_2norm_drhoa, &
      s_a_2norm_drhoanorm_drhoa, s_a_2rhoa, s_a_2rhoanorm_drhoa, &
      s_a_2rhoarhoa, s_anorm_drhoa, s_arhoa, s_arhoanorm_drhoa, s_arhoarhoa, &
      s_avg_2, s_avg_21norm_drhoa, s_avg_21norm_drhob, s_avg_21rhoa, &
      s_avg_21rhob, s_avg_2norm_drhoa, s_avg_2norm_drhoanorm_drhoa, &
      s_avg_2norm_drhob, s_avg_2norm_drhobnorm_drhob, s_avg_2rhoa, &
      s_avg_2rhoanorm_drhoa, s_avg_2rhoarhoa, s_avg_2rhob, &
      s_avg_2rhobnorm_drhob, s_avg_2rhobrhob, s_b, s_b_2, s_b_21norm_drhob, &
      s_b_21rhob, s_b_2norm_drhob, s_b_2norm_drhobnorm_drhob, s_b_2rhob, &
      s_b_2rhobnorm_drhob
    REAL(kind=dp) :: s_b_2rhobrhob, s_bnorm_drhob, s_brhob, &
      s_brhobnorm_drhob, s_brhobrhob, scale_c, scale_x, t1, t101, t1012, &
      t1014, t102, t1047, t1049, t105, t106, t107, t108, t110, t1107, t112, &
      t113, t1136, t1152, t1157, t116, t1161, t1165, t117, t1172, t1173, &
      t1175, t12, t120, t1201, t1205, t122, t1236, t125, t127, t1270, t128, &
      t129, t1299, t13, t1321, t1324, t133, t134, t1341, t1348, t1351, t1360, &
      t1368, t138, t1388, t139, t1394, t14, t1406, t1407, t1419, t142, t1422, &
      t1436, t1437, t144, t1440, t1451, t1452, t1455, t1457, t1458, t147, &
      t149, t15, t150, t151, t155, t156, t1571, t1589, t1590
    REAL(kind=dp) :: t1593, t16, t160, t1605, t161, t162, t163, t164, t165, &
      t166, t167, t168, t170, t1719, t173, t1738, t1753, t176, t18, t186, &
      t188, t191, t192, t194, t196, t197, t198, t199, t2, t20, t207, t208, &
      t209, t21, t210, t212, t216, t220, t221, t222, t223, t224, t228, t232, &
      t235, t236, t237, t239, t243, t244, t245, t246, t25, t250, t256, t257, &
      t258, t26, t260, t264, t265, t266, t267, t27, t271, t277, t278, t279, &
      t28, t285, t287, t289, t29, t290, t291, t293, t294, t295, t297, t299, &
      t3, t301, t302, t304, t31, t312, t313, t314, t315, t316, t320, t324, &
      t327, t328, t33, t336, t337, t338
    REAL(kind=dp) :: t339, t34, t344, t345, t346, t347, t35, t36, t365, t367, &
      t37, t370, t371, t374, t375, t376, t377, t396, t4, t40, t410, t42, &
      t424, t43, t431, t433, t435, t437, t438, t439, t441, t443, t445, t446, &
      t448, t456, t457, t458, t459, t46, t460, t464, t468, t471, t472, t48, &
      t480, t484, t485, t486, t49, t50, t51, t512, t516, t539, t543, t55, &
      t555, t56, t560, t564, t568, t575, t576, t577, t579, t6, t60, t600, &
      t601, t605, t606, t608, t619, t62, t621, t626, t627, t631, t632, t633, &
      t639, t644, t646, t647, t659, t66, t661, t662, t663, t67, t671, t673, &
      t678, t679, t68, t683, t689, t69, t694
    REAL(kind=dp) :: t7, t707, t709, t710, t711, t719, t721, t726, t727, t73, &
      t731, t737, t74, t742, t755, t757, t758, t759, t764, t765, t766, t771, &
      t772, t78, t790, t793, t796, t8, t80, t811, t818, t821, t830, t838, &
      t84, t85, t858, t86, t864, t87, t876, t877, t889, t892, t906, t907, &
      t91, t911, t913, t914, t92, t925, t926, t929, t930, t932, t933, t94, &
      t97, t974, t976, t98, t981, t99, t993, u_c_a, u_c_a1rhoa, u_c_ab, &
      u_c_ab1rhoa, u_c_ab1rhob, u_c_abnorm_drhoa, u_c_abnorm_drhoanorm_drhoa, &
      u_c_abnorm_drhoanorm_drhob, u_c_abnorm_drhob, &
      u_c_abnorm_drhobnorm_drhob, u_c_abrhoa, u_c_abrhoanorm_drhoa
    REAL(kind=dp) :: u_c_abrhoanorm_drhob, u_c_abrhoarhoa, u_c_abrhoarhob, &
      u_c_abrhob, u_c_abrhobnorm_drhoa, u_c_abrhobnorm_drhob, u_c_abrhobrhob, &
      u_c_anorm_drhoa, u_c_anorm_drhoanorm_drhoa, u_c_arhoa, &
      u_c_arhoanorm_drhoa, u_c_arhoarhoa, u_c_b, u_c_b1rhob, u_c_bnorm_drhob, &
      u_c_bnorm_drhobnorm_drhob, u_c_brhob, u_c_brhobnorm_drhob, &
      u_c_brhobrhob, u_x_a, u_x_a1rhoa, u_x_anorm_drhoa, &
      u_x_anorm_drhoanorm_drhoa, u_x_arhoa, u_x_arhoanorm_drhoa, &
      u_x_arhoarhoa, u_x_b, u_x_b1rhob, u_x_bnorm_drhob, &
      u_x_bnorm_drhobnorm_drhob, u_x_brhob, u_x_brhobnorm_drhob, u_x_brhobrhob
    REAL(kind=dp), DIMENSION(10)             :: coeffs

    failure=.FALSE.

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
    A_3 = 0.16887e-1_dp
    alpha_1_3 = 0.11125e0_dp
    beta_1_3 = 0.10357e2_dp
    beta_2_3 = 0.36231e1_dp
    beta_3_3 = 0.88026e0_dp
    beta_4_3 = 0.49671e0_dp

    coeffs=b97_coeffs(param,error)
    c_x_0=coeffs(1)
    c_x_1=coeffs(2)
    c_x_2=coeffs(3)
    c_cab_0=coeffs(4)
    c_cab_1=coeffs(5)
    c_cab_2=coeffs(6)
    c_css_0=coeffs(7)
    c_css_1=coeffs(8)
    c_css_2=coeffs(9)

    scale_c = scale_c_in
    scale_x = scale_x_in
    IF (scale_x==-1.0_dp) scale_x=coeffs(10)

    gamma_x    = 0.004_dp
    gamma_c_ss = 0.2_dp
    gamma_c_ab = 0.006_dp

    SELECT CASE(grad_deriv)
    CASE default
       t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
       t2 = 4 ** (0.1e1_dp / 0.3e1_dp)
       t3 = t2 ** 2
       t4 = t1 * t3
       t6 = (0.1e1_dp / pi) ** (0.1e1_dp / 0.3e1_dp)
       !$omp do
       DO ii=1,npoints
          my_rhoa=MAX(rhoa(ii),0.0_dp)
          my_rhob=MAX(rhob(ii),0.0_dp)
          rho=my_rhoa+my_rhob
          IF (rho>epsilon_rho) THEN
             my_rhoa=MAX(my_rhoa,0.5_dp*epsilon_rho)
             my_rhob=MAX(my_rhob,0.5_dp*epsilon_rho)
             rho=my_rhoa+my_rhob
             my_norm_drhoa = norm_drhoa(ii)
             my_norm_drhob = norm_drhob(ii)
             t7 = my_rhoa ** (0.1e1_dp / 0.3e1_dp)
             t8 = t7 * my_rhoa
             e_lsda_x_a = -0.3e1_dp / 0.8e1_dp * t4 * t6 * t8
             t12 = 0.1e1_dp / t8
             s_a = my_norm_drhoa * t12
             t13 = s_a ** 2
             t14 = gamma_x * t13
             t15 = 0.1e1_dp + t14
             t16 = 0.1e1_dp / t15
             u_x_a = t14 * t16
             t18 = c_x_1 + u_x_a * c_x_2
             gx_a = c_x_0 + u_x_a * t18
             t20 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
             t21 = t20 * my_rhob
             e_lsda_x_b = -0.3e1_dp / 0.8e1_dp * t4 * t6 * t21
             t25 = 0.1e1_dp / t21
             s_b = my_norm_drhob * t25
             t26 = s_b ** 2
             t27 = gamma_x * t26
             t28 = 0.1e1_dp + t27
             t29 = 0.1e1_dp / t28
             u_x_b = t27 * t29
             t31 = c_x_1 + u_x_b * c_x_2
             gx_b = c_x_0 + u_x_b * t31
             t33 = my_rhoa - my_rhob
             t34 = 0.1e1_dp / rho
             chi = t33 * t34
             t35 = 0.1e1_dp / pi
             t36 = t35 * t34
             t37 = t36 ** (0.1e1_dp / 0.3e1_dp)
             rs = t4 * t37 / 0.4e1_dp
             t40 = 0.1e1_dp + alpha_1_1 * rs
             t42 = 0.1e1_dp / A_1
             t43 = SQRT(rs)
             t46 = t43 * rs
             t48 = p_1 + 0.1e1_dp
             t49 = rs ** t48
             t50 = beta_4_1 * t49
             t51 = beta_1_1 * t43 + beta_2_1 * rs + beta_3_1 * t46 + t50
             t55 = 0.1e1_dp + t42 / t51 / 0.2e1_dp
             t56 = LOG(t55)
             e_c_u_0 = -0.2e1_dp * A_1 * t40 * t56
             t60 = 0.1e1_dp + alpha_1_2 * rs
             t62 = 0.1e1_dp / A_2
             t66 = p_2 + 0.1e1_dp
             t67 = rs ** t66
             t68 = beta_4_2 * t67
             t69 = beta_1_2 * t43 + beta_2_2 * rs + beta_3_2 * t46 + t68
             t73 = 0.1e1_dp + t62 / t69 / 0.2e1_dp
             t74 = LOG(t73)
             t78 = 0.1e1_dp + alpha_1_3 * rs
             t80 = 0.1e1_dp / A_3
             t84 = p_3 + 0.1e1_dp
             t85 = rs ** t84
             t86 = beta_4_3 * t85
             t87 = beta_1_3 * t43 + beta_2_3 * rs + beta_3_3 * t46 + t86
             t91 = 0.1e1_dp + t80 / t87 / 0.2e1_dp
             t92 = LOG(t91)
             alpha_c = 0.2e1_dp * A_3 * t78 * t92
             t94 = 2 ** (0.1e1_dp / 0.3e1_dp)
             t97 = 1 / (2 * t94 - 2)
             t98 = 0.1e1_dp + chi
             t99 = t98 ** (0.1e1_dp / 0.3e1_dp)
             t101 = 0.1e1_dp - chi
             t102 = t101 ** (0.1e1_dp / 0.3e1_dp)
             f = (t99 * t98 + t102 * t101 - 0.2e1_dp) * t97
             t105 = alpha_c * f
             t106 = 0.9e1_dp / 0.8e1_dp / t97
             t107 = chi ** 2
             t108 = t107 ** 2
             t110 = t106 * (0.1e1_dp - t108)
             t112 = -0.2e1_dp * A_2 * t60 * t74 - e_c_u_0
             t113 = t112 * f
             epsilon_c_unif = e_c_u_0 + t105 * t110 + t113 * t108
             t116 = t35 / my_rhoa
             t117 = t116 ** (0.1e1_dp / 0.3e1_dp)
             rs_a = t4 * t117 / 0.4e1_dp
             t120 = 0.1e1_dp + alpha_1_2 * rs_a
             t122 = SQRT(rs_a)
             t125 = t122 * rs_a
             t127 = rs_a ** t66
             t128 = beta_4_2 * t127
             t129 = beta_1_2 * t122 + beta_2_2 * rs_a + beta_3_2 * t125 + t128
             t133 = 0.1e1_dp + t62 / t129 / 0.2e1_dp
             t134 = LOG(t133)
             epsilon_c_unif_a = -0.2e1_dp * A_2 * t120 * t134
             t138 = t35 / my_rhob
             t139 = t138 ** (0.1e1_dp / 0.3e1_dp)
             rs_b = t4 * t139 / 0.4e1_dp
             t142 = 0.1e1_dp + alpha_1_2 * rs_b
             t144 = SQRT(rs_b)
             t147 = t144 * rs_b
             t149 = rs_b ** t66
             t150 = beta_4_2 * t149
             t151 = beta_1_2 * t144 + beta_2_2 * rs_b + beta_3_2 * t147 + t150
             t155 = 0.1e1_dp + t62 / t151 / 0.2e1_dp
             t156 = LOG(t155)
             epsilon_c_unif_b = -0.2e1_dp * A_2 * t142 * t156
             s_a_2 = t13
             s_b_2 = t26
             s_avg_2 = s_a_2 / 0.2e1_dp + s_b_2 / 0.2e1_dp
             e_lsda_c_a = epsilon_c_unif_a * my_rhoa
             e_lsda_c_b = epsilon_c_unif_b * my_rhob
             t160 = gamma_c_ab * s_avg_2
             t161 = 0.1e1_dp + t160
             t162 = 0.1e1_dp / t161
             u_c_ab = t160 * t162
             t163 = gamma_c_ss * s_a_2
             t164 = 0.1e1_dp + t163
             t165 = 0.1e1_dp / t164
             u_c_a = t163 * t165
             t166 = gamma_c_ss * s_b_2
             t167 = 0.1e1_dp + t166
             t168 = 0.1e1_dp / t167
             u_c_b = t166 * t168
             e_lsda_c_ab = epsilon_c_unif * rho - e_lsda_c_a - e_lsda_c_b
             t170 = c_cab_1 + u_c_ab * c_cab_2
             gc_ab = c_cab_0 + u_c_ab * t170
             t173 = c_css_1 + u_c_a * c_css_2
             gc_a = c_css_0 + u_c_a * t173
             t176 = c_css_1 + u_c_b * c_css_2
             gc_b = c_css_0 + u_c_b * t176

             IF (grad_deriv>=0) THEN
                exc = scale_x * (e_lsda_x_a * gx_a + e_lsda_x_b * gx_b) + scale_c &
                     * (e_lsda_c_ab * gc_ab + e_lsda_c_a * gc_a + e_lsda_c_b * gc_b)
                e_0(ii)=e_0(ii)+exc
             END IF

             IF (grad_deriv/=0) THEN

                e_lsda_x_arhoa = -t4 * t6 * t7 / 0.2e1_dp
                t186 = my_rhoa ** 2
                t188 = 0.1e1_dp / t7 / t186
                s_arhoa = -0.4e1_dp / 0.3e1_dp * my_norm_drhoa * t188
                t191 = gamma_x * s_a
                t192 = t16 * s_arhoa
                t194 = gamma_x ** 2
                t196 = t194 * s_a_2 * s_a
                t197 = t15 ** 2
                t198 = 0.1e1_dp / t197
                t199 = t198 * s_arhoa
                u_x_arhoa = 0.2e1_dp * t191 * t192 - 0.2e1_dp * t196 * t199
                gx_arhoa = u_x_arhoa * t18 + u_x_a * u_x_arhoa * c_x_2
                t207 = rho ** 2
                t208 = 0.1e1_dp / t207
                t209 = t33 * t208
                chirhoa = t34 - t209
                t210 = t37 ** 2
                t212 = 0.1e1_dp / t210 * t35
                rsrhoa = -t4 * t212 * t208 / 0.12e2_dp
                t216 = A_1 * alpha_1_1
                t220 = t51 ** 2
                t221 = 0.1e1_dp / t220
                t222 = t40 * t221
                t223 = 0.1e1_dp / t43
                t224 = beta_1_1 * t223
                t228 = beta_3_1 * t43
                t232 = 0.1e1_dp / rs
                t235 = t224 * rsrhoa / 0.2e1_dp + beta_2_1 * rsrhoa + 0.3e1_dp / &
                     0.2e1_dp * t228 * rsrhoa + t50 * t48 * rsrhoa * t232
                t236 = 0.1e1_dp / t55
                t237 = t235 * t236
                e_c_u_0rhoa = -0.2e1_dp * t216 * rsrhoa * t56 + t222 * t237
                t239 = A_2 * alpha_1_2
                t243 = t69 ** 2
                t244 = 0.1e1_dp / t243
                t245 = t60 * t244
                t246 = beta_1_2 * t223
                t250 = beta_3_2 * t43
                t256 = t246 * rsrhoa / 0.2e1_dp + beta_2_2 * rsrhoa + 0.3e1_dp / &
                     0.2e1_dp * t250 * rsrhoa + t68 * t66 * rsrhoa * t232
                t257 = 0.1e1_dp / t73
                t258 = t256 * t257
                e_c_u_1rhoa = -0.2e1_dp * t239 * rsrhoa * t74 + t245 * t258
                t260 = A_3 * alpha_1_3
                t264 = t87 ** 2
                t265 = 0.1e1_dp / t264
                t266 = t78 * t265
                t267 = beta_1_3 * t223
                t271 = beta_3_3 * t43
                t277 = t267 * rsrhoa / 0.2e1_dp + beta_2_3 * rsrhoa + 0.3e1_dp / &
                     0.2e1_dp * t271 * rsrhoa + t86 * t84 * rsrhoa * t232
                t278 = 0.1e1_dp / t91
                t279 = t277 * t278
                alpha_crhoa = 0.2e1_dp * t260 * rsrhoa * t92 - t266 * t279
                frhoa = (0.4e1_dp / 0.3e1_dp * t99 * chirhoa - 0.4e1_dp / 0.3e1_dp&
                     * t102 * chirhoa) * t97
                t285 = alpha_crhoa * f
                t287 = alpha_c * frhoa
                t289 = t107 * chi
                t290 = t106 * t289
                t291 = t290 * chirhoa
                t293 = 0.4e1_dp * t105 * t291
                t294 = e_c_u_1rhoa - e_c_u_0rhoa
                t295 = t294 * f
                t297 = t112 * frhoa
                t299 = t289 * chirhoa
                t301 = 0.4e1_dp * t113 * t299
                epsilon_c_unifrhoa = e_c_u_0rhoa + t285 * t110 + t287 * t110 - &
                     t293 + t295 * t108 + t297 * t108 + t301
                t302 = t117 ** 2
                t304 = 0.1e1_dp / t302 * t35
                rs_arhoa = -t4 * t304 / t186 / 0.12e2_dp
                t312 = t129 ** 2
                t313 = 0.1e1_dp / t312
                t314 = t120 * t313
                t315 = 0.1e1_dp / t122
                t316 = beta_1_2 * t315
                t320 = beta_3_2 * t122
                t324 = 0.1e1_dp / rs_a
                t327 = t316 * rs_arhoa / 0.2e1_dp + beta_2_2 * rs_arhoa + 0.3e1_dp&
                     / 0.2e1_dp * t320 * rs_arhoa + t128 * t66 * rs_arhoa * t324
                t328 = 0.1e1_dp / t133
                epsilon_c_unif_arhoa = -0.2e1_dp * t239 * rs_arhoa * t134 + t314 *&
                     t327 * t328
                s_a_2rhoa = 0.2e1_dp * s_a * s_arhoa
                s_avg_2rhoa = s_a_2rhoa / 0.2e1_dp
                e_lsda_c_arhoa = epsilon_c_unif_arhoa * my_rhoa + epsilon_c_unif_a
                t336 = gamma_c_ab ** 2
                t337 = t336 * s_avg_2
                t338 = t161 ** 2
                t339 = 0.1e1_dp / t338
                u_c_abrhoa = gamma_c_ab * s_avg_2rhoa * t162 - t337 * t339 * s_avg_2rhoa
                t344 = gamma_c_ss ** 2
                t345 = t344 * s_a_2
                t346 = t164 ** 2
                t347 = 0.1e1_dp / t346
                u_c_arhoa = gamma_c_ss * s_a_2rhoa * t165 - t345 * t347 * s_a_2rhoa
                e_lsda_c_abrhoa = epsilon_c_unifrhoa * rho + epsilon_c_unif - &
                     e_lsda_c_arhoa
                gc_abrhoa = u_c_abrhoa * t170 + u_c_ab * u_c_abrhoa * c_cab_2
                gc_arhoa = u_c_arhoa * t173 + u_c_a * u_c_arhoa * c_css_2

                IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                   exc_rhoa = scale_x * (e_lsda_x_arhoa * gx_a + e_lsda_x_a * &
                        gx_arhoa) + scale_c * (e_lsda_c_abrhoa * gc_ab + e_lsda_c_ab * &
                        gc_abrhoa + e_lsda_c_arhoa * gc_a + e_lsda_c_a * gc_arhoa)
                   e_ra(ii)=e_ra(ii)+exc_rhoa
                END IF

                e_lsda_x_brhob = -t4 * t6 * t20 / 0.2e1_dp
                t365 = my_rhob ** 2
                t367 = 0.1e1_dp / t20 / t365
                s_brhob = -0.4e1_dp / 0.3e1_dp * my_norm_drhob * t367
                t370 = gamma_x * s_b
                t371 = t29 * s_brhob
                t374 = t194 * s_b_2 * s_b
                t375 = t28 ** 2
                t376 = 0.1e1_dp / t375
                t377 = t376 * s_brhob
                u_x_brhob = 0.2e1_dp * t370 * t371 - 0.2e1_dp * t374 * t377
                gx_brhob = u_x_brhob * t31 + u_x_b * u_x_brhob * c_x_2
                chirhob = -t34 - t209
                rsrhob = rsrhoa
                t396 = t224 * rsrhob / 0.2e1_dp + beta_2_1 * rsrhob + 0.3e1_dp / &
                     0.2e1_dp * t228 * rsrhob + t50 * t48 * rsrhob * t232
                e_c_u_0rhob = -0.2e1_dp * t216 * rsrhob * t56 + t222 * t396 * t236
                t410 = t246 * rsrhob / 0.2e1_dp + beta_2_2 * rsrhob + 0.3e1_dp / &
                     0.2e1_dp * t250 * rsrhob + t68 * t66 * rsrhob * t232
                e_c_u_1rhob = -0.2e1_dp * t239 * rsrhob * t74 + t245 * t410 * t257
                t424 = t267 * rsrhob / 0.2e1_dp + beta_2_3 * rsrhob + 0.3e1_dp / &
                     0.2e1_dp * t271 * rsrhob + t86 * t84 * rsrhob * t232
                alpha_crhob = 0.2e1_dp * t260 * rsrhob * t92 - t266 * t424 * t278
                frhob = (0.4e1_dp / 0.3e1_dp * t99 * chirhob - 0.4e1_dp / 0.3e1_dp&
                     * t102 * chirhob) * t97
                t431 = alpha_crhob * f
                t433 = alpha_c * frhob
                t435 = t290 * chirhob
                t437 = 0.4e1_dp * t105 * t435
                t438 = e_c_u_1rhob - e_c_u_0rhob
                t439 = t438 * f
                t441 = t112 * frhob
                t443 = t289 * chirhob
                t445 = 0.4e1_dp * t113 * t443
                epsilon_c_unifrhob = e_c_u_0rhob + t431 * t110 + t433 * t110 - &
                     t437 + t439 * t108 + t441 * t108 + t445
                t446 = t139 ** 2
                t448 = 0.1e1_dp / t446 * t35
                rs_brhob = -t4 * t448 / t365 / 0.12e2_dp
                t456 = t151 ** 2
                t457 = 0.1e1_dp / t456
                t458 = t142 * t457
                t459 = 0.1e1_dp / t144
                t460 = beta_1_2 * t459
                t464 = beta_3_2 * t144
                t468 = 0.1e1_dp / rs_b
                t471 = t460 * rs_brhob / 0.2e1_dp + beta_2_2 * rs_brhob + 0.3e1_dp&
                     / 0.2e1_dp * rs_brhob * t464 + t150 * t66 * rs_brhob * t468
                t472 = 0.1e1_dp / t155
                epsilon_c_unif_brhob = -0.2e1_dp * t239 * rs_brhob * t156 + t458 *&
                     t471 * t472
                s_b_2rhob = 0.2e1_dp * s_b * s_brhob
                s_avg_2rhob = s_b_2rhob / 0.2e1_dp
                e_lsda_c_brhob = epsilon_c_unif_brhob * my_rhob + epsilon_c_unif_b
                t480 = t339 * s_avg_2rhob
                u_c_abrhob = gamma_c_ab * s_avg_2rhob * t162 - t337 * t480
                t484 = t344 * s_b_2
                t485 = t167 ** 2
                t486 = 0.1e1_dp / t485
                u_c_brhob = gamma_c_ss * s_b_2rhob * t168 - t484 * t486 * s_b_2rhob
                e_lsda_c_abrhob = epsilon_c_unifrhob * rho + epsilon_c_unif - &
                     e_lsda_c_brhob
                gc_abrhob = u_c_abrhob * t170 + u_c_ab * u_c_abrhob * c_cab_2
                gc_brhob = u_c_brhob * t176 + u_c_b * u_c_brhob * c_css_2

                IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                   exc_rhob = scale_x * (e_lsda_x_brhob * gx_b + e_lsda_x_b * &
                        gx_brhob) + scale_c * (e_lsda_c_abrhob * gc_ab + e_lsda_c_ab * &
                        gc_abrhob + e_lsda_c_brhob * gc_b + e_lsda_c_b * gc_brhob)
                   e_rb(ii)=e_rb(ii)+exc_rhob
                END IF

                s_anorm_drhoa = t12
                u_x_anorm_drhoa = 0.2e1_dp * t191 * t16 * s_anorm_drhoa - 0.2e1_dp&
                     * t196 * t198 * s_anorm_drhoa
                gx_anorm_drhoa = u_x_anorm_drhoa * t18 + u_x_a * u_x_anorm_drhoa * c_x_2
                s_a_2norm_drhoa = 0.2e1_dp * s_a * s_anorm_drhoa
                s_avg_2norm_drhoa = s_a_2norm_drhoa / 0.2e1_dp
                t512 = t339 * s_avg_2norm_drhoa
                u_c_abnorm_drhoa = gamma_c_ab * s_avg_2norm_drhoa * t162 - t337 * t512
                t516 = t347 * s_a_2norm_drhoa
                u_c_anorm_drhoa = gamma_c_ss * s_a_2norm_drhoa * t165 - t345 * t516
                gc_abnorm_drhoa = u_c_abnorm_drhoa * t170 + u_c_ab * &
                     u_c_abnorm_drhoa * c_cab_2
                gc_anorm_drhoa = u_c_anorm_drhoa * t173 + u_c_a * u_c_anorm_drhoa &
                     * c_css_2

                IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                   exc_norm_drhoa = scale_x * e_lsda_x_a * gx_anorm_drhoa + scale_c *&
                        (e_lsda_c_ab * gc_abnorm_drhoa + e_lsda_c_a * gc_anorm_drhoa)
                   e_ndra(ii)=e_ndra(ii)+exc_norm_drhoa
                END IF

                s_bnorm_drhob = t25
                u_x_bnorm_drhob = 0.2e1_dp * t370 * t29 * s_bnorm_drhob - 0.2e1_dp&
                     * t374 * t376 * s_bnorm_drhob
                gx_bnorm_drhob = u_x_bnorm_drhob * t31 + u_x_b * u_x_bnorm_drhob * c_x_2
                s_b_2norm_drhob = 0.2e1_dp * s_b * s_bnorm_drhob
                s_avg_2norm_drhob = s_b_2norm_drhob / 0.2e1_dp
                t539 = t339 * s_avg_2norm_drhob
                u_c_abnorm_drhob = gamma_c_ab * s_avg_2norm_drhob * t162 - t337 * t539
                t543 = t486 * s_b_2norm_drhob
                u_c_bnorm_drhob = gamma_c_ss * s_b_2norm_drhob * t168 - t484 * t543
                gc_abnorm_drhob = u_c_abnorm_drhob * t170 + u_c_ab * &
                     u_c_abnorm_drhob * c_cab_2
                gc_bnorm_drhob = u_c_bnorm_drhob * t176 + u_c_b * u_c_bnorm_drhob &
                     * c_css_2

                IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                   exc_norm_drhob = scale_x * e_lsda_x_b * gx_bnorm_drhob + scale_c *&
                        (e_lsda_c_ab * gc_abnorm_drhob + e_lsda_c_b * gc_bnorm_drhob)
                   e_ndrb(ii)=e_ndrb(ii)+exc_norm_drhob
                END IF

                IF (grad_deriv>1 .OR. grad_deriv<-1) THEN
                   t555 = t7 ** 2
                   t560 = t186 * my_rhoa
                   s_arhoarhoa = 0.28e2_dp / 0.9e1_dp * my_norm_drhoa / t7 / t560
                   t564 = s_arhoa ** 2
                   t568 = t194 * s_a_2
                   t575 = t194 * gamma_x
                   t576 = s_a_2 ** 2
                   t577 = t575 * t576
                   t579 = 0.1e1_dp / t197 / t15
                   u_x_arhoarhoa = 0.2e1_dp * gamma_x * t564 * t16 - 0.10e2_dp * t568&
                        * t198 * t564 + 0.2e1_dp * t191 * t16 * s_arhoarhoa + 0.8e1_dp * &
                        t577 * t579 * t564 - 0.2e1_dp * t196 * t198 * s_arhoarhoa
                   u_x_a1rhoa = u_x_arhoa
                   t600 = 0.1e1_dp / t207 / rho
                   t601 = t33 * t600
                   chirhoarhoa = -0.2e1_dp * t208 + 0.2e1_dp * t601
                   t605 = 0.3141592654e1_dp ** 2
                   t606 = 0.1e1_dp / t605
                   t608 = t207 ** 2
                   rsrhoarhoa = -t4 / t210 / t36 * t606 / t608 / 0.18e2_dp + &
                        t4 * t212 * t600 / 0.6e1_dp
                   t619 = alpha_1_1 * rsrhoa
                   t621 = t221 * t235 * t236
                   t626 = t40 / t220 / t51
                   t627 = t235 ** 2
                   t631 = 0.1e1_dp / t46
                   t632 = beta_1_1 * t631
                   t633 = rsrhoa ** 2
                   t639 = beta_3_1 * t223
                   t644 = t48 ** 2
                   t646 = rs ** 2
                   t647 = 0.1e1_dp / t646
                   t659 = t220 ** 2
                   t661 = t40 / t659
                   t662 = t55 ** 2
                   t663 = 0.1e1_dp / t662
                   e_c_u_0rhoarhoa = -0.2e1_dp * t216 * rsrhoarhoa * t56 + 0.2e1_dp *&
                        t619 * t621 - 0.2e1_dp * t626 * t627 * t236 + t222 * (-t632 * t633 &
                        / 0.4e1_dp + t224 * rsrhoarhoa / 0.2e1_dp + beta_2_1 * rsrhoarhoa + &
                        0.3e1_dp / 0.4e1_dp * t639 * t633 + 0.3e1_dp / 0.2e1_dp * t228 * &
                        rsrhoarhoa + t50 * t644 * t633 * t647 + t50 * t48 * rsrhoarhoa * &
                        t232 - t50 * t48 * t633 * t647) * t236 + t661 * t627 * t663 * t42 / &
                        0.2e1_dp
                   e_c_u_01rhoa = e_c_u_0rhoa
                   t671 = alpha_1_2 * rsrhoa
                   t673 = t244 * t256 * t257
                   t678 = t60 / t243 / t69
                   t679 = t256 ** 2
                   t683 = beta_1_2 * t631
                   t689 = beta_3_2 * t223
                   t694 = t66 ** 2
                   t707 = t243 ** 2
                   t709 = t60 / t707
                   t710 = t73 ** 2
                   t711 = 0.1e1_dp / t710
                   t719 = alpha_1_3 * rsrhoa
                   t721 = t265 * t277 * t278
                   t726 = t78 / t264 / t87
                   t727 = t277 ** 2
                   t731 = beta_1_3 * t631
                   t737 = beta_3_3 * t223
                   t742 = t84 ** 2
                   t755 = t264 ** 2
                   t757 = t78 / t755
                   t758 = t91 ** 2
                   t759 = 0.1e1_dp / t758
                   alpha_c1rhoa = alpha_crhoa
                   t764 = t99 ** 2
                   t765 = 0.1e1_dp / t764
                   t766 = chirhoa ** 2
                   t771 = t102 ** 2
                   t772 = 0.1e1_dp / t771
                   frhoarhoa = (0.4e1_dp / 0.9e1_dp * t765 * t766 + 0.4e1_dp / &
                        0.3e1_dp * t99 * chirhoarhoa + 0.4e1_dp / 0.9e1_dp * t772 * t766 - &
                        0.4e1_dp / 0.3e1_dp * t102 * chirhoarhoa) * t97
                   f1rhoa = frhoa
                   t790 = alpha_c1rhoa * f
                   t793 = alpha_c * f1rhoa
                   t796 = t106 * t107
                   t811 = e_c_u_1rhoa - e_c_u_01rhoa
                   t818 = t811 * f
                   t821 = t112 * f1rhoa
                   t830 = -0.4e1_dp * t105 * t290 * chirhoarhoa + (-0.2e1_dp * t239 *&
                        rsrhoarhoa * t74 + 0.2e1_dp * t671 * t673 - 0.2e1_dp * t678 * t679 &
                        * t257 + t245 * (-t683 * t633 / 0.4e1_dp + t246 * rsrhoarhoa / &
                        0.2e1_dp + beta_2_2 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t689 * t633&
                        + 0.3e1_dp / 0.2e1_dp * t250 * rsrhoarhoa + t68 * t694 * t633 * &
                        t647 + t68 * t66 * rsrhoarhoa * t232 - t68 * t66 * t633 * t647) * &
                        t257 + t709 * t679 * t711 * t62 / 0.2e1_dp - e_c_u_0rhoarhoa) * f * &
                        t108 + t294 * f1rhoa * t108 + 0.4e1_dp * t295 * t299 + t811 * frhoa &
                        * t108 + t112 * frhoarhoa * t108 + 0.4e1_dp * t297 * t299 + 0.4e1_dp&
                        * t818 * t299 + 0.4e1_dp * t821 * t299 + 0.12e2_dp * t113 * t107 * &
                        t766 + 0.4e1_dp * t113 * t289 * chirhoarhoa
                   epsilon_c_unif1rhoa = e_c_u_01rhoa + t790 * t110 + t793 * t110 - &
                        t293 + t818 * t108 + t821 * t108 + t301
                   t838 = t186 ** 2
                   rs_arhoarhoa = -t4 / t302 / t116 * t606 / t838 / 0.18e2_dp +&
                        t4 * t304 / t560 / 0.6e1_dp
                   t858 = t327 ** 2
                   t864 = rs_arhoa ** 2
                   t876 = rs_a ** 2
                   t877 = 0.1e1_dp / t876
                   t889 = t312 ** 2
                   t892 = t133 ** 2
                   epsilon_c_unif_a1rhoa = epsilon_c_unif_arhoa
                   s_a_2rhoarhoa = 0.2e1_dp * t564 + 0.2e1_dp * s_a * s_arhoarhoa
                   s_a_21rhoa = s_a_2rhoa
                   s_avg_2rhoarhoa = s_a_2rhoarhoa / 0.2e1_dp
                   s_avg_21rhoa = s_a_21rhoa / 0.2e1_dp
                   e_lsda_c_arhoarhoa = (-0.2e1_dp * t239 * rs_arhoarhoa * t134 + &
                        0.2e1_dp * alpha_1_2 * rs_arhoa * t313 * t327 * t328 - 0.2e1_dp * &
                        t120 / t312 / t129 * t858 * t328 + t314 * (-beta_1_2 / t125 * t864 /&
                        0.4e1_dp + t316 * rs_arhoarhoa / 0.2e1_dp + beta_2_2 * rs_arhoarhoa&
                        + 0.3e1_dp / 0.4e1_dp * beta_3_2 * t315 * t864 + 0.3e1_dp / &
                        0.2e1_dp * t320 * rs_arhoarhoa + t128 * t694 * t864 * t877 + t128 * &
                        t66 * rs_arhoarhoa * t324 - t128 * t66 * t864 * t877) * t328 + t120 &
                        / t889 * t858 / t892 * t62 / 0.2e1_dp) * my_rhoa + epsilon_c_unif_arhoa&
                        + epsilon_c_unif_a1rhoa
                   e_lsda_c_a1rhoa = epsilon_c_unif_a1rhoa * my_rhoa + epsilon_c_unif_a
                   t906 = t336 * s_avg_2rhoa
                   t907 = t339 * s_avg_21rhoa
                   t911 = t336 * gamma_c_ab * s_avg_2
                   t913 = 0.1e1_dp / t338 / t161
                   t914 = t913 * s_avg_2rhoa
                   u_c_abrhoarhoa = gamma_c_ab * s_avg_2rhoarhoa * t162 - 0.2e1_dp * &
                        t906 * t907 + 0.2e1_dp * t911 * t914 * s_avg_21rhoa - t337 * t339 * &
                        s_avg_2rhoarhoa
                   u_c_ab1rhoa = gamma_c_ab * s_avg_21rhoa * t162 - t337 * t907
                   t925 = t344 * s_a_2rhoa
                   t926 = t347 * s_a_21rhoa
                   t929 = t344 * gamma_c_ss
                   t930 = t929 * s_a_2
                   t932 = 0.1e1_dp / t346 / t164
                   t933 = t932 * s_a_2rhoa
                   u_c_arhoarhoa = gamma_c_ss * s_a_2rhoarhoa * t165 - 0.2e1_dp * &
                        t925 * t926 + 0.2e1_dp * t930 * t933 * s_a_21rhoa - t345 * t347 * &
                        s_a_2rhoarhoa
                   u_c_a1rhoa = gamma_c_ss * s_a_21rhoa * t165 - t345 * t926
                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhoa_rhoa = scale_x * (-t4 * t6 / t555 * gx_a / 0.6e1_dp&
                            + e_lsda_x_arhoa * (u_x_a1rhoa * t18 + u_x_a * u_x_a1rhoa * c_x_2) &
                            + e_lsda_x_arhoa * gx_arhoa + e_lsda_x_a * (u_x_arhoarhoa * t18 + &
                            0.2e1_dp * u_x_arhoa * u_x_a1rhoa * c_x_2 + u_x_a * u_x_arhoarhoa * &
                            c_x_2)) + scale_c * (((e_c_u_0rhoarhoa + (0.2e1_dp * t260 * &
                            rsrhoarhoa * t92 - 0.2e1_dp * t719 * t721 + 0.2e1_dp * t726 * t727 *&
                            t278 - t266 * (-t731 * t633 / 0.4e1_dp + t267 * rsrhoarhoa / &
                            0.2e1_dp + beta_2_3 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t737 * t633&
                            + 0.3e1_dp / 0.2e1_dp * t271 * rsrhoarhoa + t86 * t742 * t633 * &
                            t647 + t86 * t84 * rsrhoarhoa * t232 - t86 * t84 * t633 * t647) * &
                            t278 - t757 * t727 * t759 * t80 / 0.2e1_dp) * f * t110 + alpha_crhoa&
                            * f1rhoa * t110 - 0.4e1_dp * t285 * t291 + alpha_c1rhoa * frhoa * &
                            t110 + alpha_c * frhoarhoa * t110 - 0.4e1_dp * t287 * t291 - &
                            0.4e1_dp * t790 * t291 - 0.4e1_dp * t793 * t291 - 0.12e2_dp * t105 *&
                            t796 * t766 + t830) * rho + epsilon_c_unifrhoa + &
                            epsilon_c_unif1rhoa - e_lsda_c_arhoarhoa) * gc_ab + e_lsda_c_abrhoa &
                            * (u_c_ab1rhoa * t170 + u_c_ab * u_c_ab1rhoa * c_cab_2) + (&
                            epsilon_c_unif1rhoa * rho + epsilon_c_unif - e_lsda_c_a1rhoa) * &
                            gc_abrhoa + e_lsda_c_ab * (u_c_abrhoarhoa * t170 + 0.2e1_dp * &
                            u_c_abrhoa * u_c_ab1rhoa * c_cab_2 + u_c_ab * u_c_abrhoarhoa * &
                            c_cab_2) + e_lsda_c_arhoarhoa * gc_a + e_lsda_c_arhoa * (u_c_a1rhoa &
                            * t173 + u_c_a * u_c_a1rhoa * c_css_2) + e_lsda_c_a1rhoa * gc_arhoa &
                            + e_lsda_c_a * (u_c_arhoarhoa * t173 + 0.2e1_dp * u_c_arhoa * &
                            u_c_a1rhoa * c_css_2 + u_c_a * u_c_arhoarhoa * c_css_2))
                       e_ra_ra(ii)=e_ra_ra(ii)+exc_rhoa_rhoa
                   END IF

                   chirhoarhob = 0.2e1_dp * t601
                   rsrhoarhob = rsrhoarhoa
                   t974 = t221 * t396 * t236
                   t976 = alpha_1_1 * rsrhob
                   t981 = rsrhoa * rsrhob
                   t993 = rsrhob * t647 * rsrhoa
                   e_c_u_0rhoarhob = -0.2e1_dp * t216 * rsrhoarhob * t56 + t619 * &
                        t974 + t976 * t621 - 0.2e1_dp * t626 * t237 * t396 + t222 * (-t632 *&
                        t981 / 0.4e1_dp + t224 * rsrhoarhob / 0.2e1_dp + beta_2_1 * &
                        rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t639 * t981 + 0.3e1_dp / 0.2e1_dp&
                        * t228 * rsrhoarhob + t50 * t644 * t993 + t50 * t48 * rsrhoarhob * &
                        t232 - t50 * t48 * t993) * t236 + t661 * t235 * t663 * t42 * t396 / &
                        0.2e1_dp
                   t1012 = t244 * t410 * t257
                   t1014 = alpha_1_2 * rsrhob
                   t1047 = t265 * t424 * t278
                   t1049 = alpha_1_3 * rsrhob
                   frhoarhob = (0.4e1_dp / 0.9e1_dp * t765 * chirhoa * chirhob + &
                        0.4e1_dp / 0.3e1_dp * t99 * chirhoarhob + 0.4e1_dp / 0.9e1_dp * t772&
                        * chirhoa * chirhob - 0.4e1_dp / 0.3e1_dp * t102 * chirhoarhob) * &
                        t97
                   t1107 = t107 * chirhoa * chirhob
                   t1136 = -0.4e1_dp * t105 * t290 * chirhoarhob + (-0.2e1_dp * t239 &
                        * rsrhoarhob * t74 + t671 * t1012 + t1014 * t673 - 0.2e1_dp * t678 *&
                        t258 * t410 + t245 * (-t683 * t981 / 0.4e1_dp + t246 * rsrhoarhob /&
                        0.2e1_dp + beta_2_2 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t689 * &
                        t981 + 0.3e1_dp / 0.2e1_dp * t250 * rsrhoarhob + t68 * t694 * t993 +&
                        t68 * t66 * rsrhoarhob * t232 - t68 * t66 * t993) * t257 + t709 * &
                        t256 * t711 * t62 * t410 / 0.2e1_dp - e_c_u_0rhoarhob) * f * t108 + &
                        t294 * frhob * t108 + 0.4e1_dp * t295 * t443 + t438 * frhoa * t108 +&
                        t112 * frhoarhob * t108 + 0.4e1_dp * t297 * t443 + 0.4e1_dp * t439 &
                        * t299 + 0.4e1_dp * t441 * t299 + 0.12e2_dp * t113 * t1107 + &
                        0.4e1_dp * t113 * t289 * chirhoarhob
                   u_c_abrhoarhob = -0.2e1_dp * t906 * t480 + 0.2e1_dp * t911 * t914 &
                        * s_avg_2rhob

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhoa_rhob = scale_c * (((e_c_u_0rhoarhob + (0.2e1_dp * t260 * &
                            rsrhoarhob * t92 - t719 * t1047 - t1049 * t721 + 0.2e1_dp * t726 * &
                            t279 * t424 - t266 * (-t731 * t981 / 0.4e1_dp + t267 * rsrhoarhob / &
                            0.2e1_dp + beta_2_3 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t737 * t981&
                            + 0.3e1_dp / 0.2e1_dp * t271 * rsrhoarhob + t86 * t742 * t993 + t86&
                            * t84 * rsrhoarhob * t232 - t86 * t84 * t993) * t278 - t757 * t277 &
                            * t759 * t80 * t424 / 0.2e1_dp) * f * t110 + alpha_crhoa * frhob * &
                            t110 - 0.4e1_dp * t285 * t435 + alpha_crhob * frhoa * t110 + alpha_c&
                            * frhoarhob * t110 - 0.4e1_dp * t287 * t435 - 0.4e1_dp * t431 * &
                            t291 - 0.4e1_dp * t433 * t291 - 0.12e2_dp * t105 * t106 * t1107 + &
                            t1136) * rho + epsilon_c_unifrhoa + epsilon_c_unifrhob) * gc_ab + &
                            e_lsda_c_abrhoa * gc_abrhob + e_lsda_c_abrhob * gc_abrhoa + &
                            e_lsda_c_ab * (u_c_abrhoarhob * t170 + 0.2e1_dp * u_c_abrhoa * &
                            u_c_abrhob * c_cab_2 + u_c_ab * u_c_abrhoarhob * c_cab_2))
                       e_ra_rb(ii)=e_ra_rb(ii)+exc_rhoa_rhob
                   END IF

                   t1152 = t20 ** 2
                   t1157 = t365 * my_rhob
                   s_brhobrhob = 0.28e2_dp / 0.9e1_dp * my_norm_drhob / t20 / t1157
                   t1161 = s_brhob ** 2
                   t1165 = t194 * s_b_2
                   t1172 = s_b_2 ** 2
                   t1173 = t575 * t1172
                   t1175 = 0.1e1_dp / t375 / t28
                   u_x_brhobrhob = 0.2e1_dp * gamma_x * t1161 * t29 - 0.10e2_dp * &
                        t1165 * t376 * t1161 + 0.2e1_dp * t370 * t29 * s_brhobrhob + &
                        0.8e1_dp * t1173 * t1175 * t1161 - 0.2e1_dp * t374 * t376 * &
                        s_brhobrhob
                   u_x_b1rhob = u_x_brhob
                   chirhobrhob = 0.2e1_dp * t208 + 0.2e1_dp * t601
                   rsrhobrhob = rsrhoarhob
                   t1201 = t396 ** 2
                   t1205 = rsrhob ** 2
                   e_c_u_0rhobrhob = -0.2e1_dp * t216 * rsrhobrhob * t56 + 0.2e1_dp *&
                        t976 * t974 - 0.2e1_dp * t626 * t1201 * t236 + t222 * (-t632 * &
                        t1205 / 0.4e1_dp + t224 * rsrhobrhob / 0.2e1_dp + beta_2_1 * &
                        rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t639 * t1205 + 0.3e1_dp / &
                        0.2e1_dp * t228 * rsrhobrhob + t50 * t644 * t1205 * t647 + t50 * t48&
                        * rsrhobrhob * t232 - t50 * t48 * t1205 * t647) * t236 + t661 * &
                        t1201 * t663 * t42 / 0.2e1_dp
                   e_c_u_01rhob = e_c_u_0rhob
                   t1236 = t410 ** 2
                   t1270 = t424 ** 2
                   alpha_c1rhob = alpha_crhob
                   t1299 = chirhob ** 2
                   frhobrhob = (0.4e1_dp / 0.9e1_dp * t765 * t1299 + 0.4e1_dp / &
                        0.3e1_dp * t99 * chirhobrhob + 0.4e1_dp / 0.9e1_dp * t772 * t1299 - &
                        0.4e1_dp / 0.3e1_dp * t102 * chirhobrhob) * t97
                   f1rhob = frhob
                   t1321 = alpha_c1rhob * f
                   t1324 = alpha_c * f1rhob
                   t1341 = e_c_u_1rhob - e_c_u_01rhob
                   t1348 = t1341 * f
                   t1351 = t112 * f1rhob
                   t1360 = -0.4e1_dp * t105 * t290 * chirhobrhob + (-0.2e1_dp * t239 &
                        * rsrhobrhob * t74 + 0.2e1_dp * t1014 * t1012 - 0.2e1_dp * t678 * &
                        t1236 * t257 + t245 * (-t683 * t1205 / 0.4e1_dp + t246 * rsrhobrhob &
                        / 0.2e1_dp + beta_2_2 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t689 * &
                        t1205 + 0.3e1_dp / 0.2e1_dp * t250 * rsrhobrhob + t68 * t694 * t1205&
                        * t647 + t68 * t66 * rsrhobrhob * t232 - t68 * t66 * t1205 * t647) &
                        * t257 + t709 * t1236 * t711 * t62 / 0.2e1_dp - e_c_u_0rhobrhob) * f&
                        * t108 + t438 * f1rhob * t108 + 0.4e1_dp * t439 * t443 + t1341 * &
                        frhob * t108 + t112 * frhobrhob * t108 + 0.4e1_dp * t441 * t443 + &
                        0.4e1_dp * t1348 * t443 + 0.4e1_dp * t1351 * t443 + 0.12e2_dp * t113&
                        * t107 * t1299 + 0.4e1_dp * t113 * t289 * chirhobrhob
                   epsilon_c_unif1rhob = e_c_u_01rhob + t1321 * t110 + t1324 * t110 -&
                        t437 + t1348 * t108 + t1351 * t108 + t445
                   t1368 = t365 ** 2
                   rs_brhobrhob = -t4 / t446 / t138 * t606 / t1368 / 0.18e2_dp &
                        + t4 * t448 / t1157 / 0.6e1_dp
                   t1388 = t471 ** 2
                   t1394 = rs_brhob ** 2
                   t1406 = rs_b ** 2
                   t1407 = 0.1e1_dp / t1406
                   t1419 = t456 ** 2
                   t1422 = t155 ** 2
                   epsilon_c_unif_b1rhob = epsilon_c_unif_brhob
                   s_b_2rhobrhob = 0.2e1_dp * t1161 + 0.2e1_dp * s_b * s_brhobrhob
                   s_b_21rhob = s_b_2rhob
                   s_avg_2rhobrhob = s_b_2rhobrhob / 0.2e1_dp
                   s_avg_21rhob = s_b_21rhob / 0.2e1_dp
                   e_lsda_c_brhobrhob = (-0.2e1_dp * t239 * rs_brhobrhob * t156 + &
                        0.2e1_dp * alpha_1_2 * rs_brhob * t457 * t471 * t472 - 0.2e1_dp * &
                        t142 / t456 / t151 * t1388 * t472 + t458 * (-beta_1_2 / t147 * t1394&
                        / 0.4e1_dp + t460 * rs_brhobrhob / 0.2e1_dp + beta_2_2 * &
                        rs_brhobrhob + 0.3e1_dp / 0.4e1_dp * beta_3_2 * t459 * t1394 + &
                        0.3e1_dp / 0.2e1_dp * t464 * rs_brhobrhob + t150 * t694 * t1394 * &
                        t1407 + t150 * t66 * rs_brhobrhob * t468 - t150 * t66 * t1394 * &
                        t1407) * t472 + t142 / t1419 * t1388 / t1422 * t62 / 0.2e1_dp) * &
                        my_rhob + epsilon_c_unif_brhob + epsilon_c_unif_b1rhob
                   e_lsda_c_b1rhob = epsilon_c_unif_b1rhob * my_rhob + epsilon_c_unif_b
                   t1436 = t336 * s_avg_2rhob
                   t1437 = t339 * s_avg_21rhob
                   t1440 = t913 * s_avg_2rhob
                   u_c_abrhobrhob = gamma_c_ab * s_avg_2rhobrhob * t162 - 0.2e1_dp * &
                        t1436 * t1437 + 0.2e1_dp * t911 * t1440 * s_avg_21rhob - t337 * t339&
                        * s_avg_2rhobrhob
                   u_c_ab1rhob = gamma_c_ab * s_avg_21rhob * t162 - t337 * t1437
                   t1451 = t344 * s_b_2rhob
                   t1452 = t486 * s_b_21rhob
                   t1455 = t929 * s_b_2
                   t1457 = 0.1e1_dp / t485 / t167
                   t1458 = t1457 * s_b_2rhob
                   u_c_brhobrhob = gamma_c_ss * s_b_2rhobrhob * t168 - 0.2e1_dp * &
                        t1451 * t1452 + 0.2e1_dp * t1455 * t1458 * s_b_21rhob - t484 * t486 &
                        * s_b_2rhobrhob
                   u_c_b1rhob = gamma_c_ss * s_b_21rhob * t168 - t484 * t1452

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhob_rhob = scale_x * (-t4 * t6 / t1152 * gx_b / &
                            0.6e1_dp + e_lsda_x_brhob * (u_x_b1rhob * t31 + u_x_b * u_x_b1rhob *&
                            c_x_2) + e_lsda_x_brhob * gx_brhob + e_lsda_x_b * (u_x_brhobrhob * &
                            t31 + 0.2e1_dp * u_x_brhob * u_x_b1rhob * c_x_2 + u_x_b * &
                            u_x_brhobrhob * c_x_2)) + scale_c * (((e_c_u_0rhobrhob + (0.2e1_dp *&
                            t260 * rsrhobrhob * t92 - 0.2e1_dp * t1049 * t1047 + 0.2e1_dp * &
                            t726 * t1270 * t278 - t266 * (-t731 * t1205 / 0.4e1_dp + t267 * &
                            rsrhobrhob / 0.2e1_dp + beta_2_3 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp &
                            * t737 * t1205 + 0.3e1_dp / 0.2e1_dp * t271 * rsrhobrhob + t86 * &
                            t742 * t1205 * t647 + t86 * t84 * rsrhobrhob * t232 - t86 * t84 * &
                            t1205 * t647) * t278 - t757 * t1270 * t759 * t80 / 0.2e1_dp) * f * &
                            t110 + alpha_crhob * f1rhob * t110 - 0.4e1_dp * t431 * t435 + &
                            alpha_c1rhob * frhob * t110 + alpha_c * frhobrhob * t110 - 0.4e1_dp &
                            * t433 * t435 - 0.4e1_dp * t1321 * t435 - 0.4e1_dp * t1324 * t435 - &
                            0.12e2_dp * t105 * t796 * t1299 + t1360) * rho + epsilon_c_unifrhob &
                            + epsilon_c_unif1rhob - e_lsda_c_brhobrhob) * gc_ab + &
                            e_lsda_c_abrhob * (u_c_ab1rhob * t170 + u_c_ab * u_c_ab1rhob * &
                            c_cab_2) + (epsilon_c_unif1rhob * rho + epsilon_c_unif - &
                            e_lsda_c_b1rhob) * gc_abrhob + e_lsda_c_ab * (u_c_abrhobrhob * t170 &
                            + 0.2e1_dp * u_c_abrhob * u_c_ab1rhob * c_cab_2 + u_c_ab * &
                            u_c_abrhobrhob * c_cab_2) + e_lsda_c_brhobrhob * gc_b + &
                            e_lsda_c_brhob * (u_c_b1rhob * t176 + u_c_b * u_c_b1rhob * c_css_2) &
                            + e_lsda_c_b1rhob * gc_brhob + e_lsda_c_b * (u_c_brhobrhob * t176 + &
                            0.2e1_dp * u_c_brhob * u_c_b1rhob * c_css_2 + u_c_b * u_c_brhobrhob &
                            * c_css_2))
                       e_rb_rb(ii)=e_rb_rb(ii)+exc_rhob_rhob
                   END IF

                   s_arhoanorm_drhoa = -0.4e1_dp / 0.3e1_dp * t188
                   u_x_arhoanorm_drhoa = 0.2e1_dp * gamma_x * s_anorm_drhoa * t192 - &
                        0.10e2_dp * t568 * t199 * s_anorm_drhoa + 0.2e1_dp * t191 * t16 * &
                        s_arhoanorm_drhoa + 0.8e1_dp * t577 * t579 * s_arhoa * s_anorm_drhoa&
                        - 0.2e1_dp * t196 * t198 * s_arhoanorm_drhoa
                   s_a_2rhoanorm_drhoa = 0.2e1_dp * s_anorm_drhoa * s_arhoa + &
                        0.2e1_dp * s_a * s_arhoanorm_drhoa
                   s_avg_2rhoanorm_drhoa = s_a_2rhoanorm_drhoa / 0.2e1_dp
                   u_c_abrhoanorm_drhoa = gamma_c_ab * s_avg_2rhoanorm_drhoa * t162 -&
                        0.2e1_dp * t906 * t512 + 0.2e1_dp * t911 * t914 * s_avg_2norm_drhoa&
                        - t337 * t339 * s_avg_2rhoanorm_drhoa
                   u_c_arhoanorm_drhoa = gamma_c_ss * s_a_2rhoanorm_drhoa * t165 - &
                        0.2e1_dp * t925 * t516 + 0.2e1_dp * t930 * t933 * s_a_2norm_drhoa - &
                        t345 * t347 * s_a_2rhoanorm_drhoa

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhoa_norm_drhoa = scale_x * (e_lsda_x_arhoa * gx_anorm_drhoa +&
                            e_lsda_x_a * (u_x_arhoanorm_drhoa * t18 + 0.2e1_dp * u_x_arhoa * &
                            u_x_anorm_drhoa * c_x_2 + u_x_a * u_x_arhoanorm_drhoa * c_x_2)) + &
                            scale_c * (e_lsda_c_abrhoa * gc_abnorm_drhoa + e_lsda_c_ab * (&
                            u_c_abrhoanorm_drhoa * t170 + 0.2e1_dp * u_c_abrhoa * &
                            u_c_abnorm_drhoa * c_cab_2 + u_c_ab * u_c_abrhoanorm_drhoa * c_cab_2&
                            ) + e_lsda_c_arhoa * gc_anorm_drhoa + e_lsda_c_a * (&
                            u_c_arhoanorm_drhoa * t173 + 0.2e1_dp * u_c_arhoa * u_c_anorm_drhoa &
                            * c_css_2 + u_c_a * u_c_arhoanorm_drhoa * c_css_2))
                       e_ra_ndra(ii)=e_ra_ndra(ii)+exc_rhoa_norm_drhoa
                   END IF

                   u_c_abrhobnorm_drhoa = -0.2e1_dp * t1436 * t512 + 0.2e1_dp * t911 &
                        * t1440 * s_avg_2norm_drhoa

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhob_norm_drhoa = scale_c * (e_lsda_c_abrhob * gc_abnorm_drhoa&
                            + e_lsda_c_ab * (u_c_abrhobnorm_drhoa * t170 + 0.2e1_dp * &
                            u_c_abrhob * u_c_abnorm_drhoa * c_cab_2 + u_c_ab * &
                            u_c_abrhobnorm_drhoa * c_cab_2))
                       e_rb_ndra(ii)=e_rb_ndra(ii)+exc_rhob_norm_drhoa
                   END IF

                   t1571 = s_anorm_drhoa ** 2
                   u_x_anorm_drhoanorm_drhoa = 0.2e1_dp * gamma_x * t1571 * t16 - &
                        0.10e2_dp * t568 * t198 * t1571 + 0.8e1_dp * t577 * t579 * t1571
                   s_a_2norm_drhoanorm_drhoa = 0.2e1_dp * t1571
                   s_a_21norm_drhoa = s_a_2norm_drhoa
                   s_avg_2norm_drhoanorm_drhoa = s_a_2norm_drhoanorm_drhoa / 0.2e1_dp
                   s_avg_21norm_drhoa = s_a_21norm_drhoa / 0.2e1_dp
                   t1589 = t336 * s_avg_2norm_drhoa
                   t1590 = t339 * s_avg_21norm_drhoa
                   t1593 = t913 * s_avg_2norm_drhoa
                   u_c_abnorm_drhoanorm_drhoa = gamma_c_ab * &
                        s_avg_2norm_drhoanorm_drhoa * t162 - 0.2e1_dp * t1589 * t1590 + &
                        0.2e1_dp * t911 * t1593 * s_avg_21norm_drhoa - t337 * t339 * &
                        s_avg_2norm_drhoanorm_drhoa
                   t1605 = t347 * s_a_21norm_drhoa
                   u_c_anorm_drhoanorm_drhoa = gamma_c_ss * s_a_2norm_drhoanorm_drhoa&
                        * t165 - 0.2e1_dp * t344 * s_a_2norm_drhoa * t1605 + 0.2e1_dp * &
                        t930 * t932 * s_a_2norm_drhoa * s_a_21norm_drhoa - t345 * t347 * &
                        s_a_2norm_drhoanorm_drhoa

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_norm_drhoa_norm_drhoa = scale_x * e_lsda_x_a * (&
                            u_x_anorm_drhoanorm_drhoa * t18 + 0.2e1_dp * u_x_anorm_drhoa ** 2 * &
                            c_x_2 + u_x_a * u_x_anorm_drhoanorm_drhoa * c_x_2) + scale_c * (&
                            e_lsda_c_ab * (u_c_abnorm_drhoanorm_drhoa * t170 + 0.2e1_dp * &
                            u_c_abnorm_drhoa * (gamma_c_ab * s_avg_21norm_drhoa * t162 - t337 * &
                            t1590) * c_cab_2 + u_c_ab * u_c_abnorm_drhoanorm_drhoa * c_cab_2) + &
                            e_lsda_c_a * (u_c_anorm_drhoanorm_drhoa * t173 + 0.2e1_dp * &
                            u_c_anorm_drhoa * (gamma_c_ss * s_a_21norm_drhoa * t165 - t345 * &
                            t1605) * c_css_2 + u_c_a * u_c_anorm_drhoanorm_drhoa * c_css_2))
                       e_ndra_ndra(ii)=e_ndra_ndra(ii)+exc_norm_drhoa_norm_drhoa
                   END IF

                   u_c_abrhoanorm_drhob = -0.2e1_dp * t906 * t539 + 0.2e1_dp * t911 *&
                        t914 * s_avg_2norm_drhob

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhoa_norm_drhob = scale_c * (e_lsda_c_abrhoa * gc_abnorm_drhob&
                            + e_lsda_c_ab * (u_c_abrhoanorm_drhob * t170 + 0.2e1_dp * &
                            u_c_abrhoa * u_c_abnorm_drhob * c_cab_2 + u_c_ab * &
                            u_c_abrhoanorm_drhob * c_cab_2))
                       e_ra_ndrb(ii)=e_ra_ndrb(ii)+exc_rhoa_norm_drhob
                   END IF

                   s_brhobnorm_drhob = -0.4e1_dp / 0.3e1_dp * t367
                   u_x_brhobnorm_drhob = 0.2e1_dp * gamma_x * s_bnorm_drhob * t371 - &
                        0.10e2_dp * t1165 * t377 * s_bnorm_drhob + 0.2e1_dp * t370 * t29 * &
                        s_brhobnorm_drhob + 0.8e1_dp * t1173 * t1175 * s_brhob * &
                        s_bnorm_drhob - 0.2e1_dp * t374 * t376 * s_brhobnorm_drhob
                   s_b_2rhobnorm_drhob = 0.2e1_dp * s_bnorm_drhob * s_brhob + &
                        0.2e1_dp * s_b * s_brhobnorm_drhob
                   s_avg_2rhobnorm_drhob = s_b_2rhobnorm_drhob / 0.2e1_dp
                   u_c_abrhobnorm_drhob = gamma_c_ab * s_avg_2rhobnorm_drhob * t162 -&
                        0.2e1_dp * t1436 * t539 + 0.2e1_dp * t911 * t1440 * &
                        s_avg_2norm_drhob - t337 * t339 * s_avg_2rhobnorm_drhob
                   u_c_brhobnorm_drhob = gamma_c_ss * s_b_2rhobnorm_drhob * t168 - &
                        0.2e1_dp * t1451 * t543 + 0.2e1_dp * t1455 * t1458 * s_b_2norm_drhob&
                        - t484 * t486 * s_b_2rhobnorm_drhob

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_rhob_norm_drhob = scale_x * (e_lsda_x_brhob * gx_bnorm_drhob +&
                            e_lsda_x_b * (u_x_brhobnorm_drhob * t31 + 0.2e1_dp * u_x_brhob * &
                            u_x_bnorm_drhob * c_x_2 + u_x_b * u_x_brhobnorm_drhob * c_x_2)) + &
                            scale_c * (e_lsda_c_abrhob * gc_abnorm_drhob + e_lsda_c_ab * (&
                            u_c_abrhobnorm_drhob * t170 + 0.2e1_dp * u_c_abrhob * &
                            u_c_abnorm_drhob * c_cab_2 + u_c_ab * u_c_abrhobnorm_drhob * c_cab_2&
                            ) + e_lsda_c_brhob * gc_bnorm_drhob + e_lsda_c_b * (&
                            u_c_brhobnorm_drhob * t176 + 0.2e1_dp * u_c_brhob * u_c_bnorm_drhob &
                            * c_css_2 + u_c_b * u_c_brhobnorm_drhob * c_css_2))
                       e_rb_ndrb(ii)=e_rb_ndrb(ii)+exc_rhob_norm_drhob
                   END IF

                   u_c_abnorm_drhoanorm_drhob = -0.2e1_dp * t1589 * t539 + 0.2e1_dp *&
                        t911 * t1593 * s_avg_2norm_drhob

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_norm_drhoa_norm_drhob = scale_c * e_lsda_c_ab * (&
                            u_c_abnorm_drhoanorm_drhob * t170 + 0.2e1_dp * u_c_abnorm_drhoa * &
                            u_c_abnorm_drhob * c_cab_2 + u_c_ab * u_c_abnorm_drhoanorm_drhob * &
                            c_cab_2)
                       e_ndra_ndrb(ii)=e_ndra_ndrb(ii)+exc_norm_drhoa_norm_drhob
                   END IF

                   t1719 = s_bnorm_drhob ** 2
                   u_x_bnorm_drhobnorm_drhob = 0.2e1_dp * gamma_x * t1719 * t29 - &
                        0.10e2_dp * t1165 * t376 * t1719 + 0.8e1_dp * t1173 * t1175 * t1719
                   s_b_2norm_drhobnorm_drhob = 0.2e1_dp * t1719
                   s_b_21norm_drhob = s_b_2norm_drhob
                   s_avg_2norm_drhobnorm_drhob = s_b_2norm_drhobnorm_drhob / 0.2e1_dp
                   s_avg_21norm_drhob = s_b_21norm_drhob / 0.2e1_dp
                   t1738 = t339 * s_avg_21norm_drhob
                   u_c_abnorm_drhobnorm_drhob = gamma_c_ab * &
                        s_avg_2norm_drhobnorm_drhob * t162 - 0.2e1_dp * t336 * &
                        s_avg_2norm_drhob * t1738 + 0.2e1_dp * t911 * t913 * &
                        s_avg_2norm_drhob * s_avg_21norm_drhob - t337 * t339 * &
                        s_avg_2norm_drhobnorm_drhob
                   t1753 = t486 * s_b_21norm_drhob
                   u_c_bnorm_drhobnorm_drhob = gamma_c_ss * s_b_2norm_drhobnorm_drhob&
                        * t168 - 0.2e1_dp * t344 * s_b_2norm_drhob * t1753 + 0.2e1_dp * &
                        t1455 * t1457 * s_b_2norm_drhob * s_b_21norm_drhob - t484 * t486 * &
                        s_b_2norm_drhobnorm_drhob

                   IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                       exc_norm_drhob_norm_drhob = scale_x * e_lsda_x_b * (&
                            u_x_bnorm_drhobnorm_drhob * t31 + 0.2e1_dp * u_x_bnorm_drhob ** 2 * &
                            c_x_2 + u_x_b * u_x_bnorm_drhobnorm_drhob * c_x_2) + scale_c * (&
                            e_lsda_c_ab * (u_c_abnorm_drhobnorm_drhob * t170 + 0.2e1_dp * &
                            u_c_abnorm_drhob * (gamma_c_ab * s_avg_21norm_drhob * t162 - t337 * &
                            t1738) * c_cab_2 + u_c_ab * u_c_abnorm_drhobnorm_drhob * c_cab_2) + &
                            e_lsda_c_b * (u_c_bnorm_drhobnorm_drhob * t176 + 0.2e1_dp * &
                            u_c_bnorm_drhob * (gamma_c_ss * s_b_21norm_drhob * t168 - t484 * &
                            t1753) * c_css_2 + u_c_b * u_c_bnorm_drhobnorm_drhob * c_css_2))
                       e_ndrb_ndrb(ii)=e_ndrb_ndrb(ii)+exc_norm_drhob_norm_drhob
                   END IF
                END IF ! <1 || >1
             END IF ! /=0
          END IF ! rho>epsilon_rho
       END DO
       !$omp end do
    END SELECT
  END SUBROUTINE b97_lsd_calc

! *****************************************************************************
!> \brief low level calculation of the b97 exchange-correlation functional for lda
!> \note slow version
!> \param rho_tot,norm_drho
!> \param norm_drhoa ,norm_drhob,norm_drho: || grad rhoa |||,| grad rhoa ||,
!>        || grad (rhoa+rhob) ||
!> \param e_ 0: adds to it the local value of the functional
!> \param e_ *: derivative of the functional wrt. to the variables
!>        named where the * is.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE b97_lda_calc(rho_tot, norm_drho,&
     e_0, e_r, e_r_r, e_ndr, e_r_ndr, e_ndr_ndr, &
     grad_deriv,npoints,epsilon_rho,epsilon_drho, &
     param, scale_c_in, scale_x_in, error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho_tot, norm_drho
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_r, e_r_r, e_ndr, &
                                                e_r_ndr, e_ndr_ndr
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(in)                      :: param
    REAL(kind=dp), INTENT(in)                :: scale_c_in, scale_x_in
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'b97_lda_calc', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: small = 1.0e-20

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: A_1, A_2, A_3, alpha_1_1, alpha_1_2, alpha_1_3, alpha_c, &
      alpha_c1rhoa, alpha_c1rhob, alpha_crhoa, alpha_crhob, beta_1_1, &
      beta_1_2, beta_1_3, beta_2_1, beta_2_2, beta_2_3, beta_3_1, beta_3_2, &
      beta_3_3, beta_4_1, beta_4_2, beta_4_3, c_cab_0, c_cab_1, c_cab_2, &
      c_css_0, c_css_1, c_css_2, c_x_0, c_x_1, c_x_2, chi, chirhoa, &
      chirhoarhoa, chirhoarhob, chirhob, chirhobrhob, e_c_u_0, e_c_u_01rhoa, &
      e_c_u_01rhob, e_c_u_0rhoa, e_c_u_0rhoarhoa, e_c_u_0rhoarhob, &
      e_c_u_0rhob, e_c_u_0rhobrhob, e_c_u_1rhoa, e_c_u_1rhob, e_lsda_c_a, &
      e_lsda_c_a1rhoa, e_lsda_c_ab, e_lsda_c_abrhoa
    REAL(kind=dp) :: e_lsda_c_abrhob, e_lsda_c_arhoa, e_lsda_c_arhoarhoa, &
      e_lsda_c_b, e_lsda_c_b1rhob, e_lsda_c_brhob, e_lsda_c_brhobrhob, &
      e_lsda_x_a, e_lsda_x_arhoa, e_lsda_x_b, e_lsda_x_brhob, epsilon_c_unif, &
      epsilon_c_unif1rhoa, epsilon_c_unif1rhob, epsilon_c_unif_a, &
      epsilon_c_unif_a1rhoa, epsilon_c_unif_arhoa, epsilon_c_unif_b, &
      epsilon_c_unif_b1rhob, epsilon_c_unif_brhob, epsilon_c_unifrhoa, &
      epsilon_c_unifrhob, exc, exc_norm_drhoa, exc_norm_drhoa_norm_drhoa, &
      exc_norm_drhoa_norm_drhob, exc_norm_drhob, exc_norm_drhob_norm_drhob, &
      exc_rhoa, exc_rhoa_norm_drhoa, exc_rhoa_norm_drhob, exc_rhoa_rhoa
    REAL(kind=dp) :: exc_rhoa_rhob, exc_rhob, exc_rhob_norm_drhoa, &
      exc_rhob_norm_drhob, exc_rhob_rhob, f, f1rhoa, f1rhob, frhoa, &
      frhoarhoa, frhoarhob, frhob, frhobrhob, gamma_c_ab, gamma_c_ss, &
      gamma_x, gc_a, gc_ab, gc_abnorm_drhoa, gc_abnorm_drhob, gc_abrhoa, &
      gc_abrhob, gc_anorm_drhoa, gc_arhoa, gc_b, gc_bnorm_drhob, gc_brhob, &
      gx_a, gx_anorm_drhoa, gx_arhoa, gx_b, gx_bnorm_drhob, gx_brhob, &
      my_norm_drhoa, my_norm_drhob, my_rhoa, my_rhob, p_1, p_2, p_3, rho, rs, &
      rs_a, rs_arhoa, rs_arhoarhoa, rs_b, rs_brhob, rs_brhobrhob, rsrhoa, &
      rsrhoarhoa, rsrhoarhob, rsrhob, rsrhobrhob, s_a, s_a_2
    REAL(kind=dp) :: s_a_21norm_drhoa, s_a_21rhoa, s_a_2norm_drhoa, &
      s_a_2norm_drhoanorm_drhoa, s_a_2rhoa, s_a_2rhoanorm_drhoa, &
      s_a_2rhoarhoa, s_anorm_drhoa, s_arhoa, s_arhoanorm_drhoa, s_arhoarhoa, &
      s_avg_2, s_avg_21norm_drhoa, s_avg_21norm_drhob, s_avg_21rhoa, &
      s_avg_21rhob, s_avg_2norm_drhoa, s_avg_2norm_drhoanorm_drhoa, &
      s_avg_2norm_drhob, s_avg_2norm_drhobnorm_drhob, s_avg_2rhoa, &
      s_avg_2rhoanorm_drhoa, s_avg_2rhoarhoa, s_avg_2rhob, &
      s_avg_2rhobnorm_drhob, s_avg_2rhobrhob, s_b, s_b_2, s_b_21norm_drhob, &
      s_b_21rhob, s_b_2norm_drhob, s_b_2norm_drhobnorm_drhob, s_b_2rhob, &
      s_b_2rhobnorm_drhob
    REAL(kind=dp) :: s_b_2rhobrhob, s_bnorm_drhob, s_brhob, &
      s_brhobnorm_drhob, s_brhobrhob, scale_c, scale_x, t1, t101, t1012, &
      t1014, t102, t1047, t1049, t105, t106, t107, t108, t110, t1107, t112, &
      t113, t1136, t1152, t1157, t116, t1161, t1165, t117, t1172, t1173, &
      t1175, t12, t120, t1201, t1205, t122, t1236, t125, t127, t1270, t128, &
      t129, t1299, t13, t1321, t1324, t133, t134, t1341, t1348, t1351, t1360, &
      t1368, t138, t1388, t139, t1394, t14, t1406, t1407, t1419, t142, t1422, &
      t1436, t1437, t144, t1440, t1451, t1452, t1455, t1457, t1458, t147, &
      t149, t15, t150, t151, t155, t156, t1571, t1589, t1590
    REAL(kind=dp) :: t1593, t16, t160, t1605, t161, t162, t163, t164, t165, &
      t166, t167, t168, t170, t1719, t173, t1738, t1753, t176, t18, t186, &
      t188, t191, t192, t194, t196, t197, t198, t199, t2, t20, t207, t208, &
      t209, t21, t210, t212, t216, t220, t221, t222, t223, t224, t228, t232, &
      t235, t236, t237, t239, t243, t244, t245, t246, t25, t250, t256, t257, &
      t258, t26, t260, t264, t265, t266, t267, t27, t271, t277, t278, t279, &
      t28, t285, t287, t289, t29, t290, t291, t293, t294, t295, t297, t299, &
      t3, t301, t302, t304, t31, t312, t313, t314, t315, t316, t320, t324, &
      t327, t328, t33, t336, t337, t338
    REAL(kind=dp) :: t339, t34, t344, t345, t346, t347, t35, t36, t365, t367, &
      t37, t370, t371, t374, t375, t376, t377, t396, t4, t40, t410, t42, &
      t424, t43, t431, t433, t435, t437, t438, t439, t441, t443, t445, t446, &
      t448, t456, t457, t458, t459, t46, t460, t464, t468, t471, t472, t48, &
      t480, t484, t485, t486, t49, t50, t51, t512, t516, t539, t543, t55, &
      t555, t56, t560, t564, t568, t575, t576, t577, t579, t6, t60, t600, &
      t601, t605, t606, t608, t619, t62, t621, t626, t627, t631, t632, t633, &
      t639, t644, t646, t647, t659, t66, t661, t662, t663, t67, t671, t673, &
      t678, t679, t68, t683, t689, t69, t694
    REAL(kind=dp) :: t7, t707, t709, t710, t711, t719, t721, t726, t727, t73, &
      t731, t737, t74, t742, t755, t757, t758, t759, t764, t765, t766, t771, &
      t772, t78, t790, t793, t796, t8, t80, t811, t818, t821, t830, t838, &
      t84, t85, t858, t86, t864, t87, t876, t877, t889, t892, t906, t907, &
      t91, t911, t913, t914, t92, t925, t926, t929, t930, t932, t933, t94, &
      t97, t974, t976, t98, t981, t99, t993, u_c_a, u_c_a1rhoa, u_c_ab, &
      u_c_ab1rhoa, u_c_ab1rhob, u_c_abnorm_drhoa, u_c_abnorm_drhoanorm_drhoa, &
      u_c_abnorm_drhoanorm_drhob, u_c_abnorm_drhob, &
      u_c_abnorm_drhobnorm_drhob, u_c_abrhoa, u_c_abrhoanorm_drhoa
    REAL(kind=dp) :: u_c_abrhoanorm_drhob, u_c_abrhoarhoa, u_c_abrhoarhob, &
      u_c_abrhob, u_c_abrhobnorm_drhoa, u_c_abrhobnorm_drhob, u_c_abrhobrhob, &
      u_c_anorm_drhoa, u_c_anorm_drhoanorm_drhoa, u_c_arhoa, &
      u_c_arhoanorm_drhoa, u_c_arhoarhoa, u_c_b, u_c_b1rhob, u_c_bnorm_drhob, &
      u_c_bnorm_drhobnorm_drhob, u_c_brhob, u_c_brhobnorm_drhob, &
      u_c_brhobrhob, u_x_a, u_x_a1rhoa, u_x_anorm_drhoa, &
      u_x_anorm_drhoanorm_drhoa, u_x_arhoa, u_x_arhoanorm_drhoa, &
      u_x_arhoarhoa, u_x_b, u_x_b1rhob, u_x_bnorm_drhob, &
      u_x_bnorm_drhobnorm_drhob, u_x_brhob, u_x_brhobnorm_drhob, u_x_brhobrhob
    REAL(kind=dp), DIMENSION(10)             :: coeffs

  failure=.FALSE.

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
  A_3 = 0.16887e-1_dp
  alpha_1_3 = 0.11125e0_dp
  beta_1_3 = 0.10357e2_dp
  beta_2_3 = 0.36231e1_dp
  beta_3_3 = 0.88026e0_dp
  beta_4_3 = 0.49671e0_dp

  coeffs=b97_coeffs(param,error)
  c_x_0=coeffs(1)
  c_x_1=coeffs(2)
  c_x_2=coeffs(3)
  c_cab_0=coeffs(4)
  c_cab_1=coeffs(5)
  c_cab_2=coeffs(6)
  c_css_0=coeffs(7)
  c_css_1=coeffs(8)
  c_css_2=coeffs(9)

  scale_c = scale_c_in
  scale_x = scale_x_in
  IF (scale_x==-1.0_dp) scale_x=coeffs(10)

  gamma_x    = 0.004_dp
  gamma_c_ss = 0.2_dp
  gamma_c_ab = 0.006_dp

  SELECT CASE(grad_deriv)
  CASE default
     t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
     t2 = 4 ** (0.1e1_dp / 0.3e1_dp)
     t3 = t2 ** 2
     t4 = t1 * t3
     t6 = (0.1e1_dp / pi) ** (0.1e1_dp / 0.3e1_dp)
     !$omp do
     DO ii=1,npoints
        my_rhoa=0.5_dp*MAX(rho_tot(ii),0.0_dp)
        my_rhob=my_rhoa
        rho=my_rhoa+my_rhob
        IF (rho>epsilon_rho) THEN
           my_rhoa=MAX(my_rhoa,0.5_dp*epsilon_rho)
           my_rhob=MAX(my_rhob,0.5_dp*epsilon_rho)
           rho=my_rhoa+my_rhob
           my_norm_drhoa = 0.5_dp*norm_drho(ii)
           my_norm_drhob = 0.5_dp*norm_drho(ii)
           t7 = my_rhoa ** (0.1e1_dp / 0.3e1_dp)
           t8 = t7 * my_rhoa
           e_lsda_x_a = -0.3e1_dp / 0.8e1_dp * t4 * t6 * t8
           t12 = 0.1e1_dp / t8
           s_a = my_norm_drhoa * t12
           t13 = s_a ** 2
           t14 = gamma_x * t13
           t15 = 0.1e1_dp + t14
           t16 = 0.1e1_dp / t15
           u_x_a = t14 * t16
           t18 = c_x_1 + u_x_a * c_x_2
           gx_a = c_x_0 + u_x_a * t18
           t20 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
           t21 = t20 * my_rhob
           e_lsda_x_b = -0.3e1_dp / 0.8e1_dp * t4 * t6 * t21
           t25 = 0.1e1_dp / t21
           s_b = my_norm_drhob * t25
           t26 = s_b ** 2
           t27 = gamma_x * t26
           t28 = 0.1e1_dp + t27
           t29 = 0.1e1_dp / t28
           u_x_b = t27 * t29
           t31 = c_x_1 + u_x_b * c_x_2
           gx_b = c_x_0 + u_x_b * t31
           t33 = my_rhoa - my_rhob
           t34 = 0.1e1_dp / rho
           chi = t33 * t34
           t35 = 0.1e1_dp / pi
           t36 = t35 * t34
           t37 = t36 ** (0.1e1_dp / 0.3e1_dp)
           rs = t4 * t37 / 0.4e1_dp
           t40 = 0.1e1_dp + alpha_1_1 * rs
           t42 = 0.1e1_dp / A_1
           t43 = SQRT(rs)
           t46 = t43 * rs
           t48 = p_1 + 0.1e1_dp
           t49 = rs ** t48
           t50 = beta_4_1 * t49
           t51 = beta_1_1 * t43 + beta_2_1 * rs + beta_3_1 * t46 + t50
           t55 = 0.1e1_dp + t42 / t51 / 0.2e1_dp
           t56 = LOG(t55)
           e_c_u_0 = -0.2e1_dp * A_1 * t40 * t56
           t60 = 0.1e1_dp + alpha_1_2 * rs
           t62 = 0.1e1_dp / A_2
           t66 = p_2 + 0.1e1_dp
           t67 = rs ** t66
           t68 = beta_4_2 * t67
           t69 = beta_1_2 * t43 + beta_2_2 * rs + beta_3_2 * t46 + t68
           t73 = 0.1e1_dp + t62 / t69 / 0.2e1_dp
           t74 = LOG(t73)
           t78 = 0.1e1_dp + alpha_1_3 * rs
           t80 = 0.1e1_dp / A_3
           t84 = p_3 + 0.1e1_dp
           t85 = rs ** t84
           t86 = beta_4_3 * t85
           t87 = beta_1_3 * t43 + beta_2_3 * rs + beta_3_3 * t46 + t86
           t91 = 0.1e1_dp + t80 / t87 / 0.2e1_dp
           t92 = LOG(t91)
           alpha_c = 0.2e1_dp * A_3 * t78 * t92
           t94 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t97 = 1 / (2 * t94 - 2)
           t98 = 0.1e1_dp + chi
           t99 = t98 ** (0.1e1_dp / 0.3e1_dp)
           t101 = 0.1e1_dp - chi
           t102 = t101 ** (0.1e1_dp / 0.3e1_dp)
           f = (t99 * t98 + t102 * t101 - 0.2e1_dp) * t97
           t105 = alpha_c * f
           t106 = 0.9e1_dp / 0.8e1_dp / t97
           t107 = chi ** 2
           t108 = t107 ** 2
           t110 = t106 * (0.1e1_dp - t108)
           t112 = -0.2e1_dp * A_2 * t60 * t74 - e_c_u_0
           t113 = t112 * f
           epsilon_c_unif = e_c_u_0 + t105 * t110 + t113 * t108
           t116 = t35 / my_rhoa
           t117 = t116 ** (0.1e1_dp / 0.3e1_dp)
           rs_a = t4 * t117 / 0.4e1_dp
           t120 = 0.1e1_dp + alpha_1_2 * rs_a
           t122 = SQRT(rs_a)
           t125 = t122 * rs_a
           t127 = rs_a ** t66
           t128 = beta_4_2 * t127
           t129 = beta_1_2 * t122 + beta_2_2 * rs_a + beta_3_2 * t125 + t128
           t133 = 0.1e1_dp + t62 / t129 / 0.2e1_dp
           t134 = LOG(t133)
           epsilon_c_unif_a = -0.2e1_dp * A_2 * t120 * t134
           t138 = t35 / my_rhob
           t139 = t138 ** (0.1e1_dp / 0.3e1_dp)
           rs_b = t4 * t139 / 0.4e1_dp
           t142 = 0.1e1_dp + alpha_1_2 * rs_b
           t144 = SQRT(rs_b)
           t147 = t144 * rs_b
           t149 = rs_b ** t66
           t150 = beta_4_2 * t149
           t151 = beta_1_2 * t144 + beta_2_2 * rs_b + beta_3_2 * t147 + t150
           t155 = 0.1e1_dp + t62 / t151 / 0.2e1_dp
           t156 = LOG(t155)
           epsilon_c_unif_b = -0.2e1_dp * A_2 * t142 * t156
           s_a_2 = t13
           s_b_2 = t26
           s_avg_2 = s_a_2 / 0.2e1_dp + s_b_2 / 0.2e1_dp
           e_lsda_c_a = epsilon_c_unif_a * my_rhoa
           e_lsda_c_b = epsilon_c_unif_b * my_rhob
           t160 = gamma_c_ab * s_avg_2
           t161 = 0.1e1_dp + t160
           t162 = 0.1e1_dp / t161
           u_c_ab = t160 * t162
           t163 = gamma_c_ss * s_a_2
           t164 = 0.1e1_dp + t163
           t165 = 0.1e1_dp / t164
           u_c_a = t163 * t165
           t166 = gamma_c_ss * s_b_2
           t167 = 0.1e1_dp + t166
           t168 = 0.1e1_dp / t167
           u_c_b = t166 * t168
           e_lsda_c_ab = epsilon_c_unif * rho - e_lsda_c_a - e_lsda_c_b
           t170 = c_cab_1 + u_c_ab * c_cab_2
           gc_ab = c_cab_0 + u_c_ab * t170
           t173 = c_css_1 + u_c_a * c_css_2
           gc_a = c_css_0 + u_c_a * t173
           t176 = c_css_1 + u_c_b * c_css_2
           gc_b = c_css_0 + u_c_b * t176

           IF (grad_deriv>=0) THEN
              exc = scale_x * (e_lsda_x_a * gx_a + e_lsda_x_b * gx_b) + scale_c &
                   * (e_lsda_c_ab * gc_ab + e_lsda_c_a * gc_a + e_lsda_c_b * gc_b)
              e_0(ii)=e_0(ii)+exc
           END IF

           IF (grad_deriv/=0) THEN

              e_lsda_x_arhoa = -t4 * t6 * t7 / 0.2e1_dp
              t186 = my_rhoa ** 2
              t188 = 0.1e1_dp / t7 / t186
              s_arhoa = -0.4e1_dp / 0.3e1_dp * my_norm_drhoa * t188
              t191 = gamma_x * s_a
              t192 = t16 * s_arhoa
              t194 = gamma_x ** 2
              t196 = t194 * s_a_2 * s_a
              t197 = t15 ** 2
              t198 = 0.1e1_dp / t197
              t199 = t198 * s_arhoa
              u_x_arhoa = 0.2e1_dp * t191 * t192 - 0.2e1_dp * t196 * t199
              gx_arhoa = u_x_arhoa * t18 + u_x_a * u_x_arhoa * c_x_2
              t207 = rho ** 2
              t208 = 0.1e1_dp / t207
              t209 = t33 * t208
              chirhoa = t34 - t209
              t210 = t37 ** 2
              t212 = 0.1e1_dp / t210 * t35
              rsrhoa = -t4 * t212 * t208 / 0.12e2_dp
              t216 = A_1 * alpha_1_1
              t220 = t51 ** 2
              t221 = 0.1e1_dp / t220
              t222 = t40 * t221
              t223 = 0.1e1_dp / t43
              t224 = beta_1_1 * t223
              t228 = beta_3_1 * t43
              t232 = 0.1e1_dp / rs
              t235 = t224 * rsrhoa / 0.2e1_dp + beta_2_1 * rsrhoa + 0.3e1_dp / &
                   0.2e1_dp * t228 * rsrhoa + t50 * t48 * rsrhoa * t232
              t236 = 0.1e1_dp / t55
              t237 = t235 * t236
              e_c_u_0rhoa = -0.2e1_dp * t216 * rsrhoa * t56 + t222 * t237
              t239 = A_2 * alpha_1_2
              t243 = t69 ** 2
              t244 = 0.1e1_dp / t243
              t245 = t60 * t244
              t246 = beta_1_2 * t223
              t250 = beta_3_2 * t43
              t256 = t246 * rsrhoa / 0.2e1_dp + beta_2_2 * rsrhoa + 0.3e1_dp / &
                   0.2e1_dp * t250 * rsrhoa + t68 * t66 * rsrhoa * t232
              t257 = 0.1e1_dp / t73
              t258 = t256 * t257
              e_c_u_1rhoa = -0.2e1_dp * t239 * rsrhoa * t74 + t245 * t258
              t260 = A_3 * alpha_1_3
              t264 = t87 ** 2
              t265 = 0.1e1_dp / t264
              t266 = t78 * t265
              t267 = beta_1_3 * t223
              t271 = beta_3_3 * t43
              t277 = t267 * rsrhoa / 0.2e1_dp + beta_2_3 * rsrhoa + 0.3e1_dp / &
                   0.2e1_dp * t271 * rsrhoa + t86 * t84 * rsrhoa * t232
              t278 = 0.1e1_dp / t91
              t279 = t277 * t278
              alpha_crhoa = 0.2e1_dp * t260 * rsrhoa * t92 - t266 * t279
              frhoa = (0.4e1_dp / 0.3e1_dp * t99 * chirhoa - 0.4e1_dp / 0.3e1_dp&
                   * t102 * chirhoa) * t97
              t285 = alpha_crhoa * f
              t287 = alpha_c * frhoa
              t289 = t107 * chi
              t290 = t106 * t289
              t291 = t290 * chirhoa
              t293 = 0.4e1_dp * t105 * t291
              t294 = e_c_u_1rhoa - e_c_u_0rhoa
              t295 = t294 * f
              t297 = t112 * frhoa
              t299 = t289 * chirhoa
              t301 = 0.4e1_dp * t113 * t299
              epsilon_c_unifrhoa = e_c_u_0rhoa + t285 * t110 + t287 * t110 - &
                   t293 + t295 * t108 + t297 * t108 + t301
              t302 = t117 ** 2
              t304 = 0.1e1_dp / t302 * t35
              rs_arhoa = -t4 * t304 / t186 / 0.12e2_dp
              t312 = t129 ** 2
              t313 = 0.1e1_dp / t312
              t314 = t120 * t313
              t315 = 0.1e1_dp / t122
              t316 = beta_1_2 * t315
              t320 = beta_3_2 * t122
              t324 = 0.1e1_dp / rs_a
              t327 = t316 * rs_arhoa / 0.2e1_dp + beta_2_2 * rs_arhoa + 0.3e1_dp&
                   / 0.2e1_dp * t320 * rs_arhoa + t128 * t66 * rs_arhoa * t324
              t328 = 0.1e1_dp / t133
              epsilon_c_unif_arhoa = -0.2e1_dp * t239 * rs_arhoa * t134 + t314 *&
                   t327 * t328
              s_a_2rhoa = 0.2e1_dp * s_a * s_arhoa
              s_avg_2rhoa = s_a_2rhoa / 0.2e1_dp
              e_lsda_c_arhoa = epsilon_c_unif_arhoa * my_rhoa + epsilon_c_unif_a
              t336 = gamma_c_ab ** 2
              t337 = t336 * s_avg_2
              t338 = t161 ** 2
              t339 = 0.1e1_dp / t338
              u_c_abrhoa = gamma_c_ab * s_avg_2rhoa * t162 - t337 * t339 * s_avg_2rhoa
              t344 = gamma_c_ss ** 2
              t345 = t344 * s_a_2
              t346 = t164 ** 2
              t347 = 0.1e1_dp / t346
              u_c_arhoa = gamma_c_ss * s_a_2rhoa * t165 - t345 * t347 * s_a_2rhoa
              e_lsda_c_abrhoa = epsilon_c_unifrhoa * rho + epsilon_c_unif - &
                   e_lsda_c_arhoa
              gc_abrhoa = u_c_abrhoa * t170 + u_c_ab * u_c_abrhoa * c_cab_2
              gc_arhoa = u_c_arhoa * t173 + u_c_a * u_c_arhoa * c_css_2

              IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                 exc_rhoa = scale_x * (e_lsda_x_arhoa * gx_a + e_lsda_x_a * &
                      gx_arhoa) + scale_c * (e_lsda_c_abrhoa * gc_ab + e_lsda_c_ab * &
                      gc_abrhoa + e_lsda_c_arhoa * gc_a + e_lsda_c_a * gc_arhoa)
                 e_r(ii)=e_r(ii)+0.5_dp*exc_rhoa
              END IF

              e_lsda_x_brhob = -t4 * t6 * t20 / 0.2e1_dp
              t365 = my_rhob ** 2
              t367 = 0.1e1_dp / t20 / t365
              s_brhob = -0.4e1_dp / 0.3e1_dp * my_norm_drhob * t367
              t370 = gamma_x * s_b
              t371 = t29 * s_brhob
              t374 = t194 * s_b_2 * s_b
              t375 = t28 ** 2
              t376 = 0.1e1_dp / t375
              t377 = t376 * s_brhob
              u_x_brhob = 0.2e1_dp * t370 * t371 - 0.2e1_dp * t374 * t377
              gx_brhob = u_x_brhob * t31 + u_x_b * u_x_brhob * c_x_2
              chirhob = -t34 - t209
              rsrhob = rsrhoa
              t396 = t224 * rsrhob / 0.2e1_dp + beta_2_1 * rsrhob + 0.3e1_dp / &
                   0.2e1_dp * t228 * rsrhob + t50 * t48 * rsrhob * t232
              e_c_u_0rhob = -0.2e1_dp * t216 * rsrhob * t56 + t222 * t396 * t236
              t410 = t246 * rsrhob / 0.2e1_dp + beta_2_2 * rsrhob + 0.3e1_dp / &
                   0.2e1_dp * t250 * rsrhob + t68 * t66 * rsrhob * t232
              e_c_u_1rhob = -0.2e1_dp * t239 * rsrhob * t74 + t245 * t410 * t257
              t424 = t267 * rsrhob / 0.2e1_dp + beta_2_3 * rsrhob + 0.3e1_dp / &
                   0.2e1_dp * t271 * rsrhob + t86 * t84 * rsrhob * t232
              alpha_crhob = 0.2e1_dp * t260 * rsrhob * t92 - t266 * t424 * t278
              frhob = (0.4e1_dp / 0.3e1_dp * t99 * chirhob - 0.4e1_dp / 0.3e1_dp&
                   * t102 * chirhob) * t97
              t431 = alpha_crhob * f
              t433 = alpha_c * frhob
              t435 = t290 * chirhob
              t437 = 0.4e1_dp * t105 * t435
              t438 = e_c_u_1rhob - e_c_u_0rhob
              t439 = t438 * f
              t441 = t112 * frhob
              t443 = t289 * chirhob
              t445 = 0.4e1_dp * t113 * t443
              epsilon_c_unifrhob = e_c_u_0rhob + t431 * t110 + t433 * t110 - &
                   t437 + t439 * t108 + t441 * t108 + t445
              t446 = t139 ** 2
              t448 = 0.1e1_dp / t446 * t35
              rs_brhob = -t4 * t448 / t365 / 0.12e2_dp
              t456 = t151 ** 2
              t457 = 0.1e1_dp / t456
              t458 = t142 * t457
              t459 = 0.1e1_dp / t144
              t460 = beta_1_2 * t459
              t464 = beta_3_2 * t144
              t468 = 0.1e1_dp / rs_b
              t471 = t460 * rs_brhob / 0.2e1_dp + beta_2_2 * rs_brhob + 0.3e1_dp&
                   / 0.2e1_dp * rs_brhob * t464 + t150 * t66 * rs_brhob * t468
              t472 = 0.1e1_dp / t155
              epsilon_c_unif_brhob = -0.2e1_dp * t239 * rs_brhob * t156 + t458 *&
                   t471 * t472
              s_b_2rhob = 0.2e1_dp * s_b * s_brhob
              s_avg_2rhob = s_b_2rhob / 0.2e1_dp
              e_lsda_c_brhob = epsilon_c_unif_brhob * my_rhob + epsilon_c_unif_b
              t480 = t339 * s_avg_2rhob
              u_c_abrhob = gamma_c_ab * s_avg_2rhob * t162 - t337 * t480
              t484 = t344 * s_b_2
              t485 = t167 ** 2
              t486 = 0.1e1_dp / t485
              u_c_brhob = gamma_c_ss * s_b_2rhob * t168 - t484 * t486 * s_b_2rhob
              e_lsda_c_abrhob = epsilon_c_unifrhob * rho + epsilon_c_unif - &
                   e_lsda_c_brhob
              gc_abrhob = u_c_abrhob * t170 + u_c_ab * u_c_abrhob * c_cab_2
              gc_brhob = u_c_brhob * t176 + u_c_b * u_c_brhob * c_css_2

              IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                 exc_rhob = scale_x * (e_lsda_x_brhob * gx_b + e_lsda_x_b * &
                      gx_brhob) + scale_c * (e_lsda_c_abrhob * gc_ab + e_lsda_c_ab * &
                      gc_abrhob + e_lsda_c_brhob * gc_b + e_lsda_c_b * gc_brhob)
                 e_r(ii)=e_r(ii)+0.5_dp*exc_rhob
              END IF

              s_anorm_drhoa = t12
              u_x_anorm_drhoa = 0.2e1_dp * t191 * t16 * s_anorm_drhoa - 0.2e1_dp&
                   * t196 * t198 * s_anorm_drhoa
              gx_anorm_drhoa = u_x_anorm_drhoa * t18 + u_x_a * u_x_anorm_drhoa * c_x_2
              s_a_2norm_drhoa = 0.2e1_dp * s_a * s_anorm_drhoa
              s_avg_2norm_drhoa = s_a_2norm_drhoa / 0.2e1_dp
              t512 = t339 * s_avg_2norm_drhoa
              u_c_abnorm_drhoa = gamma_c_ab * s_avg_2norm_drhoa * t162 - t337 * t512
              t516 = t347 * s_a_2norm_drhoa
              u_c_anorm_drhoa = gamma_c_ss * s_a_2norm_drhoa * t165 - t345 * t516
              gc_abnorm_drhoa = u_c_abnorm_drhoa * t170 + u_c_ab * &
                   u_c_abnorm_drhoa * c_cab_2
              gc_anorm_drhoa = u_c_anorm_drhoa * t173 + u_c_a * u_c_anorm_drhoa &
                   * c_css_2

              IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                 exc_norm_drhoa = scale_x * e_lsda_x_a * gx_anorm_drhoa + scale_c *&
                      (e_lsda_c_ab * gc_abnorm_drhoa + e_lsda_c_a * gc_anorm_drhoa)
                 e_ndr(ii)=e_ndr(ii)+0.5_dp*exc_norm_drhoa
              END IF

              s_bnorm_drhob = t25
              u_x_bnorm_drhob = 0.2e1_dp * t370 * t29 * s_bnorm_drhob - 0.2e1_dp&
                   * t374 * t376 * s_bnorm_drhob
              gx_bnorm_drhob = u_x_bnorm_drhob * t31 + u_x_b * u_x_bnorm_drhob * c_x_2
              s_b_2norm_drhob = 0.2e1_dp * s_b * s_bnorm_drhob
              s_avg_2norm_drhob = s_b_2norm_drhob / 0.2e1_dp
              t539 = t339 * s_avg_2norm_drhob
              u_c_abnorm_drhob = gamma_c_ab * s_avg_2norm_drhob * t162 - t337 * t539
              t543 = t486 * s_b_2norm_drhob
              u_c_bnorm_drhob = gamma_c_ss * s_b_2norm_drhob * t168 - t484 * t543
              gc_abnorm_drhob = u_c_abnorm_drhob * t170 + u_c_ab * &
                   u_c_abnorm_drhob * c_cab_2
              gc_bnorm_drhob = u_c_bnorm_drhob * t176 + u_c_b * u_c_bnorm_drhob &
                   * c_css_2

              IF (grad_deriv>0 .OR. grad_deriv==-1) THEN
                 exc_norm_drhob = scale_x * e_lsda_x_b * gx_bnorm_drhob + scale_c *&
                      (e_lsda_c_ab * gc_abnorm_drhob + e_lsda_c_b * gc_bnorm_drhob)
                 e_ndr(ii)=e_ndr(ii)+0.5_dp*exc_norm_drhob
              END IF

              IF (grad_deriv>1 .OR. grad_deriv<-1) THEN
                 t555 = t7 ** 2
                 t560 = t186 * my_rhoa
                 s_arhoarhoa = 0.28e2_dp / 0.9e1_dp * my_norm_drhoa / t7 / t560
                 t564 = s_arhoa ** 2
                 t568 = t194 * s_a_2
                 t575 = t194 * gamma_x
                 t576 = s_a_2 ** 2
                 t577 = t575 * t576
                 t579 = 0.1e1_dp / t197 / t15
                 u_x_arhoarhoa = 0.2e1_dp * gamma_x * t564 * t16 - 0.10e2_dp * t568&
                      * t198 * t564 + 0.2e1_dp * t191 * t16 * s_arhoarhoa + 0.8e1_dp * &
                      t577 * t579 * t564 - 0.2e1_dp * t196 * t198 * s_arhoarhoa
                 u_x_a1rhoa = u_x_arhoa
                 t600 = 0.1e1_dp / t207 / rho
                 t601 = t33 * t600
                 chirhoarhoa = -0.2e1_dp * t208 + 0.2e1_dp * t601
                 t605 = 0.3141592654e1_dp ** 2
                 t606 = 0.1e1_dp / t605
                 t608 = t207 ** 2
                 rsrhoarhoa = -t4 / t210 / t36 * t606 / t608 / 0.18e2_dp + &
                      t4 * t212 * t600 / 0.6e1_dp
                 t619 = alpha_1_1 * rsrhoa
                 t621 = t221 * t235 * t236
                 t626 = t40 / t220 / t51
                 t627 = t235 ** 2
                 t631 = 0.1e1_dp / t46
                 t632 = beta_1_1 * t631
                 t633 = rsrhoa ** 2
                 t639 = beta_3_1 * t223
                 t644 = t48 ** 2
                 t646 = rs ** 2
                 t647 = 0.1e1_dp / t646
                 t659 = t220 ** 2
                 t661 = t40 / t659
                 t662 = t55 ** 2
                 t663 = 0.1e1_dp / t662
                 e_c_u_0rhoarhoa = -0.2e1_dp * t216 * rsrhoarhoa * t56 + 0.2e1_dp *&
                      t619 * t621 - 0.2e1_dp * t626 * t627 * t236 + t222 * (-t632 * t633 &
                      / 0.4e1_dp + t224 * rsrhoarhoa / 0.2e1_dp + beta_2_1 * rsrhoarhoa + &
                      0.3e1_dp / 0.4e1_dp * t639 * t633 + 0.3e1_dp / 0.2e1_dp * t228 * &
                      rsrhoarhoa + t50 * t644 * t633 * t647 + t50 * t48 * rsrhoarhoa * &
                      t232 - t50 * t48 * t633 * t647) * t236 + t661 * t627 * t663 * t42 / &
                      0.2e1_dp
                 e_c_u_01rhoa = e_c_u_0rhoa
                 t671 = alpha_1_2 * rsrhoa
                 t673 = t244 * t256 * t257
                 t678 = t60 / t243 / t69
                 t679 = t256 ** 2
                 t683 = beta_1_2 * t631
                 t689 = beta_3_2 * t223
                 t694 = t66 ** 2
                 t707 = t243 ** 2
                 t709 = t60 / t707
                 t710 = t73 ** 2
                 t711 = 0.1e1_dp / t710
                 t719 = alpha_1_3 * rsrhoa
                 t721 = t265 * t277 * t278
                 t726 = t78 / t264 / t87
                 t727 = t277 ** 2
                 t731 = beta_1_3 * t631
                 t737 = beta_3_3 * t223
                 t742 = t84 ** 2
                 t755 = t264 ** 2
                 t757 = t78 / t755
                 t758 = t91 ** 2
                 t759 = 0.1e1_dp / t758
                 alpha_c1rhoa = alpha_crhoa
                 t764 = t99 ** 2
                 t765 = 0.1e1_dp / t764
                 t766 = chirhoa ** 2
                 t771 = t102 ** 2
                 t772 = 0.1e1_dp / t771
                 frhoarhoa = (0.4e1_dp / 0.9e1_dp * t765 * t766 + 0.4e1_dp / &
                      0.3e1_dp * t99 * chirhoarhoa + 0.4e1_dp / 0.9e1_dp * t772 * t766 - &
                      0.4e1_dp / 0.3e1_dp * t102 * chirhoarhoa) * t97
                 f1rhoa = frhoa
                 t790 = alpha_c1rhoa * f
                 t793 = alpha_c * f1rhoa
                 t796 = t106 * t107
                 t811 = e_c_u_1rhoa - e_c_u_01rhoa
                 t818 = t811 * f
                 t821 = t112 * f1rhoa
                 t830 = -0.4e1_dp * t105 * t290 * chirhoarhoa + (-0.2e1_dp * t239 *&
                      rsrhoarhoa * t74 + 0.2e1_dp * t671 * t673 - 0.2e1_dp * t678 * t679 &
                      * t257 + t245 * (-t683 * t633 / 0.4e1_dp + t246 * rsrhoarhoa / &
                      0.2e1_dp + beta_2_2 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t689 * t633&
                      + 0.3e1_dp / 0.2e1_dp * t250 * rsrhoarhoa + t68 * t694 * t633 * &
                      t647 + t68 * t66 * rsrhoarhoa * t232 - t68 * t66 * t633 * t647) * &
                      t257 + t709 * t679 * t711 * t62 / 0.2e1_dp - e_c_u_0rhoarhoa) * f * &
                      t108 + t294 * f1rhoa * t108 + 0.4e1_dp * t295 * t299 + t811 * frhoa &
                      * t108 + t112 * frhoarhoa * t108 + 0.4e1_dp * t297 * t299 + 0.4e1_dp&
                      * t818 * t299 + 0.4e1_dp * t821 * t299 + 0.12e2_dp * t113 * t107 * &
                      t766 + 0.4e1_dp * t113 * t289 * chirhoarhoa
                 epsilon_c_unif1rhoa = e_c_u_01rhoa + t790 * t110 + t793 * t110 - &
                      t293 + t818 * t108 + t821 * t108 + t301
                 t838 = t186 ** 2
                 rs_arhoarhoa = -t4 / t302 / t116 * t606 / t838 / 0.18e2_dp +&
                      t4 * t304 / t560 / 0.6e1_dp
                 t858 = t327 ** 2
                 t864 = rs_arhoa ** 2
                 t876 = rs_a ** 2
                 t877 = 0.1e1_dp / t876
                 t889 = t312 ** 2
                 t892 = t133 ** 2
                 epsilon_c_unif_a1rhoa = epsilon_c_unif_arhoa
                 s_a_2rhoarhoa = 0.2e1_dp * t564 + 0.2e1_dp * s_a * s_arhoarhoa
                 s_a_21rhoa = s_a_2rhoa
                 s_avg_2rhoarhoa = s_a_2rhoarhoa / 0.2e1_dp
                 s_avg_21rhoa = s_a_21rhoa / 0.2e1_dp
                 e_lsda_c_arhoarhoa = (-0.2e1_dp * t239 * rs_arhoarhoa * t134 + &
                      0.2e1_dp * alpha_1_2 * rs_arhoa * t313 * t327 * t328 - 0.2e1_dp * &
                      t120 / t312 / t129 * t858 * t328 + t314 * (-beta_1_2 / t125 * t864 /&
                      0.4e1_dp + t316 * rs_arhoarhoa / 0.2e1_dp + beta_2_2 * rs_arhoarhoa&
                      + 0.3e1_dp / 0.4e1_dp * beta_3_2 * t315 * t864 + 0.3e1_dp / &
                      0.2e1_dp * t320 * rs_arhoarhoa + t128 * t694 * t864 * t877 + t128 * &
                      t66 * rs_arhoarhoa * t324 - t128 * t66 * t864 * t877) * t328 + t120 &
                      / t889 * t858 / t892 * t62 / 0.2e1_dp) * my_rhoa + epsilon_c_unif_arhoa&
                      + epsilon_c_unif_a1rhoa
                 e_lsda_c_a1rhoa = epsilon_c_unif_a1rhoa * my_rhoa + epsilon_c_unif_a
                 t906 = t336 * s_avg_2rhoa
                 t907 = t339 * s_avg_21rhoa
                 t911 = t336 * gamma_c_ab * s_avg_2
                 t913 = 0.1e1_dp / t338 / t161
                 t914 = t913 * s_avg_2rhoa
                 u_c_abrhoarhoa = gamma_c_ab * s_avg_2rhoarhoa * t162 - 0.2e1_dp * &
                      t906 * t907 + 0.2e1_dp * t911 * t914 * s_avg_21rhoa - t337 * t339 * &
                      s_avg_2rhoarhoa
                 u_c_ab1rhoa = gamma_c_ab * s_avg_21rhoa * t162 - t337 * t907
                 t925 = t344 * s_a_2rhoa
                 t926 = t347 * s_a_21rhoa
                 t929 = t344 * gamma_c_ss
                 t930 = t929 * s_a_2
                 t932 = 0.1e1_dp / t346 / t164
                 t933 = t932 * s_a_2rhoa
                 u_c_arhoarhoa = gamma_c_ss * s_a_2rhoarhoa * t165 - 0.2e1_dp * &
                      t925 * t926 + 0.2e1_dp * t930 * t933 * s_a_21rhoa - t345 * t347 * &
                      s_a_2rhoarhoa
                 u_c_a1rhoa = gamma_c_ss * s_a_21rhoa * t165 - t345 * t926
                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhoa_rhoa = scale_x * (-t4 * t6 / t555 * gx_a / 0.6e1_dp&
                          + e_lsda_x_arhoa * (u_x_a1rhoa * t18 + u_x_a * u_x_a1rhoa * c_x_2) &
                          + e_lsda_x_arhoa * gx_arhoa + e_lsda_x_a * (u_x_arhoarhoa * t18 + &
                          0.2e1_dp * u_x_arhoa * u_x_a1rhoa * c_x_2 + u_x_a * u_x_arhoarhoa * &
                          c_x_2)) + scale_c * (((e_c_u_0rhoarhoa + (0.2e1_dp * t260 * &
                          rsrhoarhoa * t92 - 0.2e1_dp * t719 * t721 + 0.2e1_dp * t726 * t727 *&
                          t278 - t266 * (-t731 * t633 / 0.4e1_dp + t267 * rsrhoarhoa / &
                          0.2e1_dp + beta_2_3 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t737 * t633&
                          + 0.3e1_dp / 0.2e1_dp * t271 * rsrhoarhoa + t86 * t742 * t633 * &
                          t647 + t86 * t84 * rsrhoarhoa * t232 - t86 * t84 * t633 * t647) * &
                          t278 - t757 * t727 * t759 * t80 / 0.2e1_dp) * f * t110 + alpha_crhoa&
                          * f1rhoa * t110 - 0.4e1_dp * t285 * t291 + alpha_c1rhoa * frhoa * &
                          t110 + alpha_c * frhoarhoa * t110 - 0.4e1_dp * t287 * t291 - &
                          0.4e1_dp * t790 * t291 - 0.4e1_dp * t793 * t291 - 0.12e2_dp * t105 *&
                          t796 * t766 + t830) * rho + epsilon_c_unifrhoa + &
                          epsilon_c_unif1rhoa - e_lsda_c_arhoarhoa) * gc_ab + e_lsda_c_abrhoa &
                          * (u_c_ab1rhoa * t170 + u_c_ab * u_c_ab1rhoa * c_cab_2) + (&
                          epsilon_c_unif1rhoa * rho + epsilon_c_unif - e_lsda_c_a1rhoa) * &
                          gc_abrhoa + e_lsda_c_ab * (u_c_abrhoarhoa * t170 + 0.2e1_dp * &
                          u_c_abrhoa * u_c_ab1rhoa * c_cab_2 + u_c_ab * u_c_abrhoarhoa * &
                          c_cab_2) + e_lsda_c_arhoarhoa * gc_a + e_lsda_c_arhoa * (u_c_a1rhoa &
                          * t173 + u_c_a * u_c_a1rhoa * c_css_2) + e_lsda_c_a1rhoa * gc_arhoa &
                          + e_lsda_c_a * (u_c_arhoarhoa * t173 + 0.2e1_dp * u_c_arhoa * &
                          u_c_a1rhoa * c_css_2 + u_c_a * u_c_arhoarhoa * c_css_2))
                     e_r_r(ii)=e_r_r(ii)+0.5_dp*0.5_dp*exc_rhoa_rhoa
                 END IF

                 chirhoarhob = 0.2e1_dp * t601
                 rsrhoarhob = rsrhoarhoa
                 t974 = t221 * t396 * t236
                 t976 = alpha_1_1 * rsrhob
                 t981 = rsrhoa * rsrhob
                 t993 = rsrhob * t647 * rsrhoa
                 e_c_u_0rhoarhob = -0.2e1_dp * t216 * rsrhoarhob * t56 + t619 * &
                      t974 + t976 * t621 - 0.2e1_dp * t626 * t237 * t396 + t222 * (-t632 *&
                      t981 / 0.4e1_dp + t224 * rsrhoarhob / 0.2e1_dp + beta_2_1 * &
                      rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t639 * t981 + 0.3e1_dp / 0.2e1_dp&
                      * t228 * rsrhoarhob + t50 * t644 * t993 + t50 * t48 * rsrhoarhob * &
                      t232 - t50 * t48 * t993) * t236 + t661 * t235 * t663 * t42 * t396 / &
                      0.2e1_dp
                 t1012 = t244 * t410 * t257
                 t1014 = alpha_1_2 * rsrhob
                 t1047 = t265 * t424 * t278
                 t1049 = alpha_1_3 * rsrhob
                 frhoarhob = (0.4e1_dp / 0.9e1_dp * t765 * chirhoa * chirhob + &
                      0.4e1_dp / 0.3e1_dp * t99 * chirhoarhob + 0.4e1_dp / 0.9e1_dp * t772&
                      * chirhoa * chirhob - 0.4e1_dp / 0.3e1_dp * t102 * chirhoarhob) * &
                      t97
                 t1107 = t107 * chirhoa * chirhob
                 t1136 = -0.4e1_dp * t105 * t290 * chirhoarhob + (-0.2e1_dp * t239 &
                      * rsrhoarhob * t74 + t671 * t1012 + t1014 * t673 - 0.2e1_dp * t678 *&
                      t258 * t410 + t245 * (-t683 * t981 / 0.4e1_dp + t246 * rsrhoarhob /&
                      0.2e1_dp + beta_2_2 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t689 * &
                      t981 + 0.3e1_dp / 0.2e1_dp * t250 * rsrhoarhob + t68 * t694 * t993 +&
                      t68 * t66 * rsrhoarhob * t232 - t68 * t66 * t993) * t257 + t709 * &
                      t256 * t711 * t62 * t410 / 0.2e1_dp - e_c_u_0rhoarhob) * f * t108 + &
                      t294 * frhob * t108 + 0.4e1_dp * t295 * t443 + t438 * frhoa * t108 +&
                      t112 * frhoarhob * t108 + 0.4e1_dp * t297 * t443 + 0.4e1_dp * t439 &
                      * t299 + 0.4e1_dp * t441 * t299 + 0.12e2_dp * t113 * t1107 + &
                      0.4e1_dp * t113 * t289 * chirhoarhob
                 u_c_abrhoarhob = -0.2e1_dp * t906 * t480 + 0.2e1_dp * t911 * t914 &
                      * s_avg_2rhob

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhoa_rhob = scale_c * (((e_c_u_0rhoarhob + (0.2e1_dp * t260 * &
                          rsrhoarhob * t92 - t719 * t1047 - t1049 * t721 + 0.2e1_dp * t726 * &
                          t279 * t424 - t266 * (-t731 * t981 / 0.4e1_dp + t267 * rsrhoarhob / &
                          0.2e1_dp + beta_2_3 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t737 * t981&
                          + 0.3e1_dp / 0.2e1_dp * t271 * rsrhoarhob + t86 * t742 * t993 + t86&
                          * t84 * rsrhoarhob * t232 - t86 * t84 * t993) * t278 - t757 * t277 &
                          * t759 * t80 * t424 / 0.2e1_dp) * f * t110 + alpha_crhoa * frhob * &
                          t110 - 0.4e1_dp * t285 * t435 + alpha_crhob * frhoa * t110 + alpha_c&
                          * frhoarhob * t110 - 0.4e1_dp * t287 * t435 - 0.4e1_dp * t431 * &
                          t291 - 0.4e1_dp * t433 * t291 - 0.12e2_dp * t105 * t106 * t1107 + &
                          t1136) * rho + epsilon_c_unifrhoa + epsilon_c_unifrhob) * gc_ab + &
                          e_lsda_c_abrhoa * gc_abrhob + e_lsda_c_abrhob * gc_abrhoa + &
                          e_lsda_c_ab * (u_c_abrhoarhob * t170 + 0.2e1_dp * u_c_abrhoa * &
                          u_c_abrhob * c_cab_2 + u_c_ab * u_c_abrhoarhob * c_cab_2))
                     e_r_r(ii)=e_r_r(ii)+0.5_dp*exc_rhoa_rhob
                 END IF

                 t1152 = t20 ** 2
                 t1157 = t365 * my_rhob
                 s_brhobrhob = 0.28e2_dp / 0.9e1_dp * my_norm_drhob / t20 / t1157
                 t1161 = s_brhob ** 2
                 t1165 = t194 * s_b_2
                 t1172 = s_b_2 ** 2
                 t1173 = t575 * t1172
                 t1175 = 0.1e1_dp / t375 / t28
                 u_x_brhobrhob = 0.2e1_dp * gamma_x * t1161 * t29 - 0.10e2_dp * &
                      t1165 * t376 * t1161 + 0.2e1_dp * t370 * t29 * s_brhobrhob + &
                      0.8e1_dp * t1173 * t1175 * t1161 - 0.2e1_dp * t374 * t376 * &
                      s_brhobrhob
                 u_x_b1rhob = u_x_brhob
                 chirhobrhob = 0.2e1_dp * t208 + 0.2e1_dp * t601
                 rsrhobrhob = rsrhoarhob
                 t1201 = t396 ** 2
                 t1205 = rsrhob ** 2
                 e_c_u_0rhobrhob = -0.2e1_dp * t216 * rsrhobrhob * t56 + 0.2e1_dp *&
                      t976 * t974 - 0.2e1_dp * t626 * t1201 * t236 + t222 * (-t632 * &
                      t1205 / 0.4e1_dp + t224 * rsrhobrhob / 0.2e1_dp + beta_2_1 * &
                      rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t639 * t1205 + 0.3e1_dp / &
                      0.2e1_dp * t228 * rsrhobrhob + t50 * t644 * t1205 * t647 + t50 * t48&
                      * rsrhobrhob * t232 - t50 * t48 * t1205 * t647) * t236 + t661 * &
                      t1201 * t663 * t42 / 0.2e1_dp
                 e_c_u_01rhob = e_c_u_0rhob
                 t1236 = t410 ** 2
                 t1270 = t424 ** 2
                 alpha_c1rhob = alpha_crhob
                 t1299 = chirhob ** 2
                 frhobrhob = (0.4e1_dp / 0.9e1_dp * t765 * t1299 + 0.4e1_dp / &
                      0.3e1_dp * t99 * chirhobrhob + 0.4e1_dp / 0.9e1_dp * t772 * t1299 - &
                      0.4e1_dp / 0.3e1_dp * t102 * chirhobrhob) * t97
                 f1rhob = frhob
                 t1321 = alpha_c1rhob * f
                 t1324 = alpha_c * f1rhob
                 t1341 = e_c_u_1rhob - e_c_u_01rhob
                 t1348 = t1341 * f
                 t1351 = t112 * f1rhob
                 t1360 = -0.4e1_dp * t105 * t290 * chirhobrhob + (-0.2e1_dp * t239 &
                      * rsrhobrhob * t74 + 0.2e1_dp * t1014 * t1012 - 0.2e1_dp * t678 * &
                      t1236 * t257 + t245 * (-t683 * t1205 / 0.4e1_dp + t246 * rsrhobrhob &
                      / 0.2e1_dp + beta_2_2 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t689 * &
                      t1205 + 0.3e1_dp / 0.2e1_dp * t250 * rsrhobrhob + t68 * t694 * t1205&
                      * t647 + t68 * t66 * rsrhobrhob * t232 - t68 * t66 * t1205 * t647) &
                      * t257 + t709 * t1236 * t711 * t62 / 0.2e1_dp - e_c_u_0rhobrhob) * f&
                      * t108 + t438 * f1rhob * t108 + 0.4e1_dp * t439 * t443 + t1341 * &
                      frhob * t108 + t112 * frhobrhob * t108 + 0.4e1_dp * t441 * t443 + &
                      0.4e1_dp * t1348 * t443 + 0.4e1_dp * t1351 * t443 + 0.12e2_dp * t113&
                      * t107 * t1299 + 0.4e1_dp * t113 * t289 * chirhobrhob
                 epsilon_c_unif1rhob = e_c_u_01rhob + t1321 * t110 + t1324 * t110 -&
                      t437 + t1348 * t108 + t1351 * t108 + t445
                 t1368 = t365 ** 2
                 rs_brhobrhob = -t4 / t446 / t138 * t606 / t1368 / 0.18e2_dp &
                      + t4 * t448 / t1157 / 0.6e1_dp
                 t1388 = t471 ** 2
                 t1394 = rs_brhob ** 2
                 t1406 = rs_b ** 2
                 t1407 = 0.1e1_dp / t1406
                 t1419 = t456 ** 2
                 t1422 = t155 ** 2
                 epsilon_c_unif_b1rhob = epsilon_c_unif_brhob
                 s_b_2rhobrhob = 0.2e1_dp * t1161 + 0.2e1_dp * s_b * s_brhobrhob
                 s_b_21rhob = s_b_2rhob
                 s_avg_2rhobrhob = s_b_2rhobrhob / 0.2e1_dp
                 s_avg_21rhob = s_b_21rhob / 0.2e1_dp
                 e_lsda_c_brhobrhob = (-0.2e1_dp * t239 * rs_brhobrhob * t156 + &
                      0.2e1_dp * alpha_1_2 * rs_brhob * t457 * t471 * t472 - 0.2e1_dp * &
                      t142 / t456 / t151 * t1388 * t472 + t458 * (-beta_1_2 / t147 * t1394&
                      / 0.4e1_dp + t460 * rs_brhobrhob / 0.2e1_dp + beta_2_2 * &
                      rs_brhobrhob + 0.3e1_dp / 0.4e1_dp * beta_3_2 * t459 * t1394 + &
                      0.3e1_dp / 0.2e1_dp * t464 * rs_brhobrhob + t150 * t694 * t1394 * &
                      t1407 + t150 * t66 * rs_brhobrhob * t468 - t150 * t66 * t1394 * &
                      t1407) * t472 + t142 / t1419 * t1388 / t1422 * t62 / 0.2e1_dp) * &
                      my_rhob + epsilon_c_unif_brhob + epsilon_c_unif_b1rhob
                 e_lsda_c_b1rhob = epsilon_c_unif_b1rhob * my_rhob + epsilon_c_unif_b
                 t1436 = t336 * s_avg_2rhob
                 t1437 = t339 * s_avg_21rhob
                 t1440 = t913 * s_avg_2rhob
                 u_c_abrhobrhob = gamma_c_ab * s_avg_2rhobrhob * t162 - 0.2e1_dp * &
                      t1436 * t1437 + 0.2e1_dp * t911 * t1440 * s_avg_21rhob - t337 * t339&
                      * s_avg_2rhobrhob
                 u_c_ab1rhob = gamma_c_ab * s_avg_21rhob * t162 - t337 * t1437
                 t1451 = t344 * s_b_2rhob
                 t1452 = t486 * s_b_21rhob
                 t1455 = t929 * s_b_2
                 t1457 = 0.1e1_dp / t485 / t167
                 t1458 = t1457 * s_b_2rhob
                 u_c_brhobrhob = gamma_c_ss * s_b_2rhobrhob * t168 - 0.2e1_dp * &
                      t1451 * t1452 + 0.2e1_dp * t1455 * t1458 * s_b_21rhob - t484 * t486 &
                      * s_b_2rhobrhob
                 u_c_b1rhob = gamma_c_ss * s_b_21rhob * t168 - t484 * t1452

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhob_rhob = scale_x * (-t4 * t6 / t1152 * gx_b / &
                          0.6e1_dp + e_lsda_x_brhob * (u_x_b1rhob * t31 + u_x_b * u_x_b1rhob *&
                          c_x_2) + e_lsda_x_brhob * gx_brhob + e_lsda_x_b * (u_x_brhobrhob * &
                          t31 + 0.2e1_dp * u_x_brhob * u_x_b1rhob * c_x_2 + u_x_b * &
                          u_x_brhobrhob * c_x_2)) + scale_c * (((e_c_u_0rhobrhob + (0.2e1_dp *&
                          t260 * rsrhobrhob * t92 - 0.2e1_dp * t1049 * t1047 + 0.2e1_dp * &
                          t726 * t1270 * t278 - t266 * (-t731 * t1205 / 0.4e1_dp + t267 * &
                          rsrhobrhob / 0.2e1_dp + beta_2_3 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp &
                          * t737 * t1205 + 0.3e1_dp / 0.2e1_dp * t271 * rsrhobrhob + t86 * &
                          t742 * t1205 * t647 + t86 * t84 * rsrhobrhob * t232 - t86 * t84 * &
                          t1205 * t647) * t278 - t757 * t1270 * t759 * t80 / 0.2e1_dp) * f * &
                          t110 + alpha_crhob * f1rhob * t110 - 0.4e1_dp * t431 * t435 + &
                          alpha_c1rhob * frhob * t110 + alpha_c * frhobrhob * t110 - 0.4e1_dp &
                          * t433 * t435 - 0.4e1_dp * t1321 * t435 - 0.4e1_dp * t1324 * t435 - &
                          0.12e2_dp * t105 * t796 * t1299 + t1360) * rho + epsilon_c_unifrhob &
                          + epsilon_c_unif1rhob - e_lsda_c_brhobrhob) * gc_ab + &
                          e_lsda_c_abrhob * (u_c_ab1rhob * t170 + u_c_ab * u_c_ab1rhob * &
                          c_cab_2) + (epsilon_c_unif1rhob * rho + epsilon_c_unif - &
                          e_lsda_c_b1rhob) * gc_abrhob + e_lsda_c_ab * (u_c_abrhobrhob * t170 &
                          + 0.2e1_dp * u_c_abrhob * u_c_ab1rhob * c_cab_2 + u_c_ab * &
                          u_c_abrhobrhob * c_cab_2) + e_lsda_c_brhobrhob * gc_b + &
                          e_lsda_c_brhob * (u_c_b1rhob * t176 + u_c_b * u_c_b1rhob * c_css_2) &
                          + e_lsda_c_b1rhob * gc_brhob + e_lsda_c_b * (u_c_brhobrhob * t176 + &
                          0.2e1_dp * u_c_brhob * u_c_b1rhob * c_css_2 + u_c_b * u_c_brhobrhob &
                          * c_css_2))
                     e_r_r(ii)=e_r_r(ii)+0.5_dp*0.5_dp*exc_rhob_rhob
                 END IF

                 s_arhoanorm_drhoa = -0.4e1_dp / 0.3e1_dp * t188
                 u_x_arhoanorm_drhoa = 0.2e1_dp * gamma_x * s_anorm_drhoa * t192 - &
                      0.10e2_dp * t568 * t199 * s_anorm_drhoa + 0.2e1_dp * t191 * t16 * &
                      s_arhoanorm_drhoa + 0.8e1_dp * t577 * t579 * s_arhoa * s_anorm_drhoa&
                      - 0.2e1_dp * t196 * t198 * s_arhoanorm_drhoa
                 s_a_2rhoanorm_drhoa = 0.2e1_dp * s_anorm_drhoa * s_arhoa + &
                      0.2e1_dp * s_a * s_arhoanorm_drhoa
                 s_avg_2rhoanorm_drhoa = s_a_2rhoanorm_drhoa / 0.2e1_dp
                 u_c_abrhoanorm_drhoa = gamma_c_ab * s_avg_2rhoanorm_drhoa * t162 -&
                      0.2e1_dp * t906 * t512 + 0.2e1_dp * t911 * t914 * s_avg_2norm_drhoa&
                      - t337 * t339 * s_avg_2rhoanorm_drhoa
                 u_c_arhoanorm_drhoa = gamma_c_ss * s_a_2rhoanorm_drhoa * t165 - &
                      0.2e1_dp * t925 * t516 + 0.2e1_dp * t930 * t933 * s_a_2norm_drhoa - &
                      t345 * t347 * s_a_2rhoanorm_drhoa

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhoa_norm_drhoa = scale_x * (e_lsda_x_arhoa * gx_anorm_drhoa +&
                          e_lsda_x_a * (u_x_arhoanorm_drhoa * t18 + 0.2e1_dp * u_x_arhoa * &
                          u_x_anorm_drhoa * c_x_2 + u_x_a * u_x_arhoanorm_drhoa * c_x_2)) + &
                          scale_c * (e_lsda_c_abrhoa * gc_abnorm_drhoa + e_lsda_c_ab * (&
                          u_c_abrhoanorm_drhoa * t170 + 0.2e1_dp * u_c_abrhoa * &
                          u_c_abnorm_drhoa * c_cab_2 + u_c_ab * u_c_abrhoanorm_drhoa * c_cab_2&
                          ) + e_lsda_c_arhoa * gc_anorm_drhoa + e_lsda_c_a * (&
                          u_c_arhoanorm_drhoa * t173 + 0.2e1_dp * u_c_arhoa * u_c_anorm_drhoa &
                          * c_css_2 + u_c_a * u_c_arhoanorm_drhoa * c_css_2))
                     e_r_ndr(ii)=e_r_ndr(ii)+0.5_dp*0.5_dp*exc_rhoa_norm_drhoa
                 END IF

                 u_c_abrhobnorm_drhoa = -0.2e1_dp * t1436 * t512 + 0.2e1_dp * t911 &
                      * t1440 * s_avg_2norm_drhoa

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhob_norm_drhoa = scale_c * (e_lsda_c_abrhob * gc_abnorm_drhoa&
                          + e_lsda_c_ab * (u_c_abrhobnorm_drhoa * t170 + 0.2e1_dp * &
                          u_c_abrhob * u_c_abnorm_drhoa * c_cab_2 + u_c_ab * &
                          u_c_abrhobnorm_drhoa * c_cab_2))
                     e_r_ndr(ii)=e_r_ndr(ii)+0.5_dp*0.5_dp*exc_rhob_norm_drhoa
                 END IF

                 t1571 = s_anorm_drhoa ** 2
                 u_x_anorm_drhoanorm_drhoa = 0.2e1_dp * gamma_x * t1571 * t16 - &
                      0.10e2_dp * t568 * t198 * t1571 + 0.8e1_dp * t577 * t579 * t1571
                 s_a_2norm_drhoanorm_drhoa = 0.2e1_dp * t1571
                 s_a_21norm_drhoa = s_a_2norm_drhoa
                 s_avg_2norm_drhoanorm_drhoa = s_a_2norm_drhoanorm_drhoa / 0.2e1_dp
                 s_avg_21norm_drhoa = s_a_21norm_drhoa / 0.2e1_dp
                 t1589 = t336 * s_avg_2norm_drhoa
                 t1590 = t339 * s_avg_21norm_drhoa
                 t1593 = t913 * s_avg_2norm_drhoa
                 u_c_abnorm_drhoanorm_drhoa = gamma_c_ab * &
                      s_avg_2norm_drhoanorm_drhoa * t162 - 0.2e1_dp * t1589 * t1590 + &
                      0.2e1_dp * t911 * t1593 * s_avg_21norm_drhoa - t337 * t339 * &
                      s_avg_2norm_drhoanorm_drhoa
                 t1605 = t347 * s_a_21norm_drhoa
                 u_c_anorm_drhoanorm_drhoa = gamma_c_ss * s_a_2norm_drhoanorm_drhoa&
                      * t165 - 0.2e1_dp * t344 * s_a_2norm_drhoa * t1605 + 0.2e1_dp * &
                      t930 * t932 * s_a_2norm_drhoa * s_a_21norm_drhoa - t345 * t347 * &
                      s_a_2norm_drhoanorm_drhoa

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_norm_drhoa_norm_drhoa = scale_x * e_lsda_x_a * (&
                          u_x_anorm_drhoanorm_drhoa * t18 + 0.2e1_dp * u_x_anorm_drhoa ** 2 * &
                          c_x_2 + u_x_a * u_x_anorm_drhoanorm_drhoa * c_x_2) + scale_c * (&
                          e_lsda_c_ab * (u_c_abnorm_drhoanorm_drhoa * t170 + 0.2e1_dp * &
                          u_c_abnorm_drhoa * (gamma_c_ab * s_avg_21norm_drhoa * t162 - t337 * &
                          t1590) * c_cab_2 + u_c_ab * u_c_abnorm_drhoanorm_drhoa * c_cab_2) + &
                          e_lsda_c_a * (u_c_anorm_drhoanorm_drhoa * t173 + 0.2e1_dp * &
                          u_c_anorm_drhoa * (gamma_c_ss * s_a_21norm_drhoa * t165 - t345 * &
                          t1605) * c_css_2 + u_c_a * u_c_anorm_drhoanorm_drhoa * c_css_2))
                     e_ndr_ndr(ii)=e_ndr_ndr(ii)+0.5_dp*0.5_dp*exc_norm_drhoa_norm_drhoa
                 END IF

                 u_c_abrhoanorm_drhob = -0.2e1_dp * t906 * t539 + 0.2e1_dp * t911 *&
                      t914 * s_avg_2norm_drhob

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhoa_norm_drhob = scale_c * (e_lsda_c_abrhoa * gc_abnorm_drhob&
                          + e_lsda_c_ab * (u_c_abrhoanorm_drhob * t170 + 0.2e1_dp * &
                          u_c_abrhoa * u_c_abnorm_drhob * c_cab_2 + u_c_ab * &
                          u_c_abrhoanorm_drhob * c_cab_2))
                     e_r_ndr(ii)=e_r_ndr(ii)+0.5_dp*0.5_dp*exc_rhoa_norm_drhob
                 END IF

                 s_brhobnorm_drhob = -0.4e1_dp / 0.3e1_dp * t367
                 u_x_brhobnorm_drhob = 0.2e1_dp * gamma_x * s_bnorm_drhob * t371 - &
                      0.10e2_dp * t1165 * t377 * s_bnorm_drhob + 0.2e1_dp * t370 * t29 * &
                      s_brhobnorm_drhob + 0.8e1_dp * t1173 * t1175 * s_brhob * &
                      s_bnorm_drhob - 0.2e1_dp * t374 * t376 * s_brhobnorm_drhob
                 s_b_2rhobnorm_drhob = 0.2e1_dp * s_bnorm_drhob * s_brhob + &
                      0.2e1_dp * s_b * s_brhobnorm_drhob
                 s_avg_2rhobnorm_drhob = s_b_2rhobnorm_drhob / 0.2e1_dp
                 u_c_abrhobnorm_drhob = gamma_c_ab * s_avg_2rhobnorm_drhob * t162 -&
                      0.2e1_dp * t1436 * t539 + 0.2e1_dp * t911 * t1440 * &
                      s_avg_2norm_drhob - t337 * t339 * s_avg_2rhobnorm_drhob
                 u_c_brhobnorm_drhob = gamma_c_ss * s_b_2rhobnorm_drhob * t168 - &
                      0.2e1_dp * t1451 * t543 + 0.2e1_dp * t1455 * t1458 * s_b_2norm_drhob&
                      - t484 * t486 * s_b_2rhobnorm_drhob

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_rhob_norm_drhob = scale_x * (e_lsda_x_brhob * gx_bnorm_drhob +&
                          e_lsda_x_b * (u_x_brhobnorm_drhob * t31 + 0.2e1_dp * u_x_brhob * &
                          u_x_bnorm_drhob * c_x_2 + u_x_b * u_x_brhobnorm_drhob * c_x_2)) + &
                          scale_c * (e_lsda_c_abrhob * gc_abnorm_drhob + e_lsda_c_ab * (&
                          u_c_abrhobnorm_drhob * t170 + 0.2e1_dp * u_c_abrhob * &
                          u_c_abnorm_drhob * c_cab_2 + u_c_ab * u_c_abrhobnorm_drhob * c_cab_2&
                          ) + e_lsda_c_brhob * gc_bnorm_drhob + e_lsda_c_b * (&
                          u_c_brhobnorm_drhob * t176 + 0.2e1_dp * u_c_brhob * u_c_bnorm_drhob &
                          * c_css_2 + u_c_b * u_c_brhobnorm_drhob * c_css_2))
                     e_r_ndr(ii)=e_r_ndr(ii)+0.5_dp*0.5_dp*exc_rhob_norm_drhob
                 END IF

                 u_c_abnorm_drhoanorm_drhob = -0.2e1_dp * t1589 * t539 + 0.2e1_dp *&
                      t911 * t1593 * s_avg_2norm_drhob

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_norm_drhoa_norm_drhob = scale_c * e_lsda_c_ab * (&
                          u_c_abnorm_drhoanorm_drhob * t170 + 0.2e1_dp * u_c_abnorm_drhoa * &
                          u_c_abnorm_drhob * c_cab_2 + u_c_ab * u_c_abnorm_drhoanorm_drhob * &
                          c_cab_2)
                     e_ndr_ndr(ii)=e_ndr_ndr(ii)+0.5_dp*exc_norm_drhoa_norm_drhob
                 END IF

                 t1719 = s_bnorm_drhob ** 2
                 u_x_bnorm_drhobnorm_drhob = 0.2e1_dp * gamma_x * t1719 * t29 - &
                      0.10e2_dp * t1165 * t376 * t1719 + 0.8e1_dp * t1173 * t1175 * t1719
                 s_b_2norm_drhobnorm_drhob = 0.2e1_dp * t1719
                 s_b_21norm_drhob = s_b_2norm_drhob
                 s_avg_2norm_drhobnorm_drhob = s_b_2norm_drhobnorm_drhob / 0.2e1_dp
                 s_avg_21norm_drhob = s_b_21norm_drhob / 0.2e1_dp
                 t1738 = t339 * s_avg_21norm_drhob
                 u_c_abnorm_drhobnorm_drhob = gamma_c_ab * &
                      s_avg_2norm_drhobnorm_drhob * t162 - 0.2e1_dp * t336 * &
                      s_avg_2norm_drhob * t1738 + 0.2e1_dp * t911 * t913 * &
                      s_avg_2norm_drhob * s_avg_21norm_drhob - t337 * t339 * &
                      s_avg_2norm_drhobnorm_drhob
                 t1753 = t486 * s_b_21norm_drhob
                 u_c_bnorm_drhobnorm_drhob = gamma_c_ss * s_b_2norm_drhobnorm_drhob&
                      * t168 - 0.2e1_dp * t344 * s_b_2norm_drhob * t1753 + 0.2e1_dp * &
                      t1455 * t1457 * s_b_2norm_drhob * s_b_21norm_drhob - t484 * t486 * &
                      s_b_2norm_drhobnorm_drhob

                 IF (grad_deriv>1 .OR. grad_deriv==-2) THEN
                     exc_norm_drhob_norm_drhob = scale_x * e_lsda_x_b * (&
                          u_x_bnorm_drhobnorm_drhob * t31 + 0.2e1_dp * u_x_bnorm_drhob ** 2 * &
                          c_x_2 + u_x_b * u_x_bnorm_drhobnorm_drhob * c_x_2) + scale_c * (&
                          e_lsda_c_ab * (u_c_abnorm_drhobnorm_drhob * t170 + 0.2e1_dp * &
                          u_c_abnorm_drhob * (gamma_c_ab * s_avg_21norm_drhob * t162 - t337 * &
                          t1738) * c_cab_2 + u_c_ab * u_c_abnorm_drhobnorm_drhob * c_cab_2) + &
                          e_lsda_c_b * (u_c_bnorm_drhobnorm_drhob * t176 + 0.2e1_dp * &
                          u_c_bnorm_drhob * (gamma_c_ss * s_b_21norm_drhob * t168 - t484 * &
                          t1753) * c_css_2 + u_c_b * u_c_bnorm_drhobnorm_drhob * c_css_2))
                     e_ndr_ndr(ii)=e_ndr_ndr(ii)+0.5_dp*0.5_dp*exc_norm_drhob_norm_drhob
                 END IF
              END IF ! <1 || >1
           END IF ! /=0
        END IF ! rho>epsilon_rho
     END DO
     !$omp end do
  END SELECT
END SUBROUTINE b97_lda_calc

END MODULE xc_b97
