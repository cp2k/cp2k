!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates a functional from libxc and its derivatives
!> \note
!>      LibXC:
!>      (Marques, Oliveira, Burnus, CPC 183, 2272 (2012)).
!>
!>      For subsequent versions of libxc, the following should be updated if
!>      necessary:
!>      1) The list of functionals for which it is possible to provide input
!>         parameters (in 'xc_libxc_wrap_functional_set_params'). For more
!>         information on the parameters, see subroutines xc_f90_xxx_set_par
!>         in libxc.f90 of the libxc package or xc_f90_lib_m.F.
!>         only checked for functionals up to 2.0.1
!>      2) Reactivate the functionals which are working correctly
!>         (in 'xc_libxc_wrap_functional_buggy').
!>         only checked for functionals up to 3.x.x
!>
!>      WARNING: In the subroutine libxc_lsd_calc, it could be that the
!>      ordering for the 1st index of v2lapltau, v2rholapl, v2rhotau,
!>      v2sigmalapl and v2sigmatau is not correct. For the moment it does not
!>      matter since the calculation of the 2nd derivatives for meta-GGA
!>      functionals is not implemented in CP2K.
!>
!> \par History
!>      01.2013 created [F. Tran]
!>      07.2014 updates to versions 2.1 [JGH]
!>      08.2015 refactoring [A. Gloess (agloess)]
!> \author F. Tran
! *****************************************************************************
MODULE xc_libxc
#if defined(__LIBXC2) || defined(__LIBXC3)
#define __LIBXC
#endif
  USE bibliography,                    ONLY: Marques2012,&
                                             cite_reference
  
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
  USE termination,                     ONLY: stop_program
#if defined (__LIBXC)
  USE xc_libxc_wrap,                   ONLY: xc_f90_pointer_t,&
                                             xc_f90_func_init,&
                                             xc_f90_func_end,&
                                             xc_f90_info_family,&
                                             xc_f90_info_kind,&
                                             xc_f90_info_name,&
                                             xc_f90_gga_exc,&
                                             xc_f90_gga_exc_vxc,&
                                             xc_f90_gga_fxc,&
                                             xc_f90_gga_lb_set_par,&
                                             xc_f90_gga_vxc,&
                                             xc_f90_gga_x_hjs_set_par,&
                                             xc_f90_gga_x_wpbeh_set_par,&
                                             xc_f90_hyb_gga_xc_hse_set_par,&
                                             xc_f90_lda,&
                                             xc_f90_lda_c_1d_csc_set_par,&
                                             xc_f90_lda_c_2d_prm_set_par,&
                                             xc_f90_lda_c_xalpha_set_par,&
                                             xc_f90_lda_exc,&
                                             xc_f90_lda_exc_vxc,&
                                             xc_f90_lda_fxc,&
                                             xc_f90_lda_kxc,&
                                             xc_f90_lda_vxc,&
                                             xc_f90_lda_x_1d_set_par,&
                                             xc_f90_lda_x_set_par,&
                                             xc_f90_mgga,&
                                             xc_f90_mgga_exc,&
                                             xc_f90_mgga_exc_vxc,&
                                             xc_f90_mgga_fxc,&
                                             xc_f90_mgga_vxc,&
                                             xc_f90_mgga_x_tb09_set_par,&
                                             XC_POLARIZED,&
                                             XC_UNPOLARIZED,&
                                             XC_FAMILY_LDA,&
                                             XC_FAMILY_GGA,&
                                             XC_FAMILY_MGGA,&
                                             XC_FAMILY_HYB_GGA,&
                                             XC_CORRELATION,&
                                             XC_EXCHANGE,&
                                             XC_EXCHANGE_CORRELATION,&
                                             XC_KINETIC,&
                                             xc_libxc_wrap_info_refs,&
                                             xc_libxc_wrap_version,&
                                             xc_libxc_wrap_functional_get_number,&
                                             xc_libxc_wrap_needs_laplace,&
                                             xc_libxc_wrap_functional_set_params,&
                                             xc_libxc_wrap_functional_buggy
#endif

#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc'

  PUBLIC :: libxc_lda_info, libxc_lda_eval, libxc_lsd_info, libxc_lsd_eval, &
            libxc_version_info

CONTAINS

! *****************************************************************************
!> \brief info about the functional from libxc
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv maximum implemented derivative of the xc functional
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_info(libxc_params,reference,shortform,needs,max_deriv,ifunc_name)

    TYPE(section_vals_type), POINTER         :: libxc_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_info', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    CHARACTER(LEN=120)                       :: s1, s2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id
    REAL(KIND=dp)                            :: sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: scale
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name)
    CALL section_vals_val_get(libxc_params,"scale",r_vals=scale)

    CALL cite_reference(Marques2012)

    IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
      sc = 1.0_dp
    ELSE
      sc = SCALE(ifunc_name)
    ENDIF

    func_id = xc_libxc_wrap_functional_get_number(func_name(ifunc_name))
    CALL xc_libxc_wrap_functional_buggy(func_id)
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_UNPOLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

    CALL xc_f90_info_name(xc_info,s1)
    SELECT CASE(xc_f90_info_kind(xc_info))
    CASE (XC_EXCHANGE);             WRITE(s2,'(a)') "exchange"
    CASE (XC_CORRELATION);          WRITE(s2,'(a)') "correlation"
    CASE (XC_EXCHANGE_CORRELATION); WRITE(s2,'(a)') "exchange-correlation"
    CASE (XC_KINETIC);              WRITE(s2,'(a)') "kinetic"
    CASE default
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name(ifunc_name))//": this XC_KIND is currently not supported.")
    END SELECT
    IF ( PRESENT ( shortform ) ) THEN
       shortform = TRIM(s1)//' ('//TRIM(s2)//')'
    END IF
    IF ( PRESENT ( reference ) ) THEN
       CALL xc_libxc_wrap_info_refs(xc_info, XC_UNPOLARIZED, sc, reference)
    END IF
    IF (PRESENT(needs)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          needs%rho = .TRUE.
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          needs%rho = .TRUE.
          needs%norm_drho = .TRUE.
       CASE (XC_FAMILY_MGGA)
          needs%rho = .TRUE.
          needs%norm_drho = .TRUE.
          needs%tau = .TRUE.
          needs%laplace_rho = xc_libxc_wrap_needs_laplace(func_id)
       CASE default
          CALL stop_program(routineN,moduleN,__LINE__,&
             TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
       END SELECT
    END IF
    IF (PRESENT(max_deriv)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          max_deriv = 3
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          max_deriv = 2
       CASE (XC_FAMILY_MGGA)
          max_deriv = 1
       CASE default
          CALL stop_program(routineN,moduleN,__LINE__,&
             TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
       END SELECT
    END IF

    CALL xc_f90_func_end(xc_func)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif

  END SUBROUTINE libxc_lda_info

! *****************************************************************************
!> \brief info about the functional from libxc
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv maximum implemented derivative of the xc functional
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_info(libxc_params,reference,shortform,needs,max_deriv,ifunc_name)

    TYPE(section_vals_type), POINTER         :: libxc_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_info', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    CHARACTER(LEN=120)                       :: s1, s2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id
    REAL(KIND=dp)                            :: sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: scale
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name)
    CALL section_vals_val_get(libxc_params,"scale",r_vals=scale)

    CALL cite_reference(Marques2012)

    IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
      sc = 1.0_dp
    ELSE
      sc = SCALE(ifunc_name)
    ENDIF

    func_id = xc_libxc_wrap_functional_get_number(func_name(ifunc_name))
    CALL xc_libxc_wrap_functional_buggy(func_id)
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_POLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

    CALL xc_f90_info_name(xc_info,s1)
    SELECT CASE(xc_f90_info_kind(xc_info))
    CASE (XC_EXCHANGE);             WRITE(s2,'(a)') "exchange"
    CASE (XC_CORRELATION);          WRITE(s2,'(a)') "correlation"
    CASE (XC_EXCHANGE_CORRELATION); WRITE(s2,'(a)') "exchange-correlation"
    CASE (XC_KINETIC);              WRITE(s2,'(a)') "kinetic"
    CASE default
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name(ifunc_name))//": this XC_KIND is currently not supported.")
    END SELECT
    IF ( PRESENT ( shortform ) ) THEN
       shortform = TRIM(s1)//' ('//TRIM(s2)//')'
    END IF
    IF ( PRESENT ( reference ) ) THEN
       CALL xc_libxc_wrap_info_refs(xc_info, XC_POLARIZED, sc, reference)
    END IF
    IF (PRESENT(needs)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          needs%rho_spin = .TRUE.
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          needs%rho_spin = .TRUE.
          needs%norm_drho = .TRUE.
          needs%norm_drho_spin = .TRUE.
       CASE (XC_FAMILY_MGGA)
          needs%rho_spin = .TRUE.
          needs%norm_drho = .TRUE.
          needs%norm_drho_spin = .TRUE.
          needs%tau_spin = .TRUE.
          needs%laplace_rho_spin = xc_libxc_wrap_needs_laplace(func_id)
       CASE default
          CALL stop_program(routineN,moduleN,__LINE__,&
             TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
       END SELECT
    END IF
    IF (PRESENT(max_deriv)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          max_deriv = 3
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          max_deriv = 2
       CASE (XC_FAMILY_MGGA)
          max_deriv = 1
       CASE default
          CALL stop_program(routineN,moduleN,__LINE__,&
             TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
       END SELECT
    END IF

    CALL xc_f90_func_end(xc_func)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif

  END SUBROUTINE libxc_lsd_info

! *****************************************************************************
!> \brief info about the LibXC version
!> \param version ...
!> \author A. Gloess (agloess)
! *****************************************************************************
  SUBROUTINE libxc_version_info(version)
    CHARACTER(LEN=*), INTENT(OUT)      :: version ! the string that is output

    CHARACTER(LEN=*), PARAMETER :: routineN = 'libxc_version_info', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    CALL xc_libxc_wrap_version(version)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif

  END SUBROUTINE libxc_version_info

! *****************************************************************************
!> \brief evaluates the functional from libxc
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_eval(rho_set,deriv_set,grad_deriv,libxc_params,ifunc_name)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: libxc_params
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_eval', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id, handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure, has_laplace
    REAL(KIND=dp)                            :: epsilon_rho, epsilon_tau, sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: params, scale
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rho, &
      e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_ndrho, &
      e_ndrho_laplace_rho, e_ndrho_ndrho, e_ndrho_rho, e_ndrho_tau, e_rho, &
      e_rho_laplace_rho, e_rho_rho, e_rho_rho_rho, e_rho_tau, e_tau, &
      e_tau_tau, laplace_rho, norm_drho, rho, tau
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL timeset(routineN,handle)

    failure = .FALSE.
    has_laplace = .FALSE.
    NULLIFY(bo, dummy)
    NULLIFY(rho, norm_drho, laplace_rho, tau)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    IF (.NOT. failure) THEN

       CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name)
       CALL section_vals_val_get(libxc_params,"scale",r_vals=scale)
       CALL section_vals_val_get(libxc_params,"parameters",r_vals=params)

       IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
         sc = 1.0_dp
       ELSE
         sc = SCALE(ifunc_name)
       ENDIF

       func_id = xc_libxc_wrap_functional_get_number(func_name(ifunc_name))
       CALL xc_libxc_wrap_functional_buggy(func_id, grad_deriv=grad_deriv)
!$OMP CRITICAL(libxc_init)
       CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_UNPOLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

       CALL xc_rho_set_get(rho_set, can_return_null=.TRUE.,&
            rho=rho,norm_drho=norm_drho, laplace_rho=laplace_rho,&
            rho_cutoff=epsilon_rho, tau_cutoff=epsilon_tau,&
            tau=tau, local_bounds=bo)

       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
       ELSE
          dummy => rho
       END IF

       ! due to assumed shape array usage in next routine
       IF (.NOT.ASSOCIATED(norm_drho)) norm_drho => dummy
       IF (.NOT.ASSOCIATED(tau)) tau => dummy

       ! only some MGGA functionals really need the Laplacian,
       ! all others can work with rho (read-only) as dummy
       IF (ASSOCIATED(laplace_rho)) has_laplace = .TRUE.
       IF (.NOT. has_laplace) laplace_rho => dummy

       e_0 => dummy
       e_rho => dummy
       e_ndrho => dummy
       e_laplace_rho => dummy
       e_tau => dummy
       e_rho_rho => dummy
       e_ndrho_rho => dummy
       e_ndrho_ndrho => dummy
       e_rho_laplace_rho => dummy
       e_rho_tau => dummy
       e_ndrho_laplace_rho => dummy
       e_ndrho_tau => dummy
       e_laplace_rho_laplace_rho => dummy
       e_laplace_rho_tau => dummy
       e_tau_tau => dummy
       e_rho_rho_rho => dummy

       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE.)
          CALL xc_derivative_get(deriv,deriv_data=e_0)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rho)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rho)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
          CASE(XC_FAMILY_MGGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rho)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
             deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_tau)
             IF (has_laplace) THEN
                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)",&
                     allocate_deriv=.TRUE.)
                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho)
             END IF
          CASE default
             CALL stop_program(routineN,moduleN,__LINE__,&
                TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
          END SELECT
       END IF
       IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)
          CASE(XC_FAMILY_MGGA)
            ! not implemented ...
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives larger than 1 not implemented or checked",&
                error_level=cp_failure_level)

!             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)
!             deriv => xc_dset_get_derivative(deriv_set,"(rho)(tau)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rho_tau)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(tau)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_tau)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau)(tau)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_tau)
!             IF (has_laplace) THEN
!                deriv => xc_dset_get_derivative(deriv_set,"(rho)(laplace_rho)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_rho_laplace_rho)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_rho)(laplace_rho)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrho_laplace_rho)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)(laplace_rho)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho_laplace_rho)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)(tau)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho_tau)
!             END IF
          CASE default
             CALL stop_program(routineN,moduleN,__LINE__,&
                TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
          END SELECT
       END IF
       IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA)
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives larger than 2 not implemented",&
                error_level=cp_failure_level)
          CASE default
             CALL stop_program(routineN,moduleN,__LINE__,&
                TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
          END SELECT
       END IF
       IF (grad_deriv>=4.OR.grad_deriv<=-4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP,&
             message="derivatives larger than 3 not implemented",&
             error_level=cp_failure_level)
       END IF

       !$omp parallel default(none), &
       !$omp shared(rho,norm_drho,laplace_rho,tau,e_0,e_rho,e_ndrho,e_laplace_rho),&
       !$omp shared(e_tau,e_rho_rho,e_ndrho_rho,e_ndrho_ndrho,e_rho_laplace_rho),&
       !$omp shared(e_rho_tau,e_ndrho_laplace_rho,e_ndrho_tau,e_laplace_rho_laplace_rho),&
       !$omp shared(e_laplace_rho_tau,e_tau_tau,e_rho_rho_rho),&
       !$omp shared(grad_deriv,npoints),&
       !$omp shared(epsilon_rho,epsilon_tau),&
       !$omp shared(func_name,ifunc_name,sc,params)

       CALL libxc_lda_calc(rho=rho,norm_drho=norm_drho,&
          laplace_rho=laplace_rho,tau=tau,&
          e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_laplace_rho=e_laplace_rho,&
          e_tau=e_tau,e_rho_rho=e_rho_rho,e_ndrho_rho=e_ndrho_rho,&
          e_ndrho_ndrho=e_ndrho_ndrho,e_rho_laplace_rho=e_rho_laplace_rho,&
          e_rho_tau=e_rho_tau,e_ndrho_laplace_rho=e_ndrho_laplace_rho,&
          e_ndrho_tau=e_ndrho_tau,e_laplace_rho_laplace_rho=e_laplace_rho_laplace_rho,&
          e_laplace_rho_tau=e_laplace_rho_tau,e_tau_tau=e_tau_tau,&
          e_rho_rho_rho=e_rho_rho_rho,&
          grad_deriv=grad_deriv,npoints=npoints,&
          epsilon_rho=epsilon_rho,&
          epsilon_tau=epsilon_tau,func_name=func_name(ifunc_name),&
          sc=sc,params=params)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy)
       ELSE
          NULLIFY(dummy)
       END IF

       CALL xc_f90_func_end(xc_func)

    END IF

    CALL timestop(handle)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif
  END SUBROUTINE libxc_lda_eval

! *****************************************************************************
!> \brief evaluates the functional from libxc
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param ifunc_name the index of the functional as given in the input file
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_eval(rho_set,deriv_set,grad_deriv,libxc_params,ifunc_name)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: libxc_params
    INTEGER, INTENT(in)                      :: ifunc_name

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_eval', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id, handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure, has_laplace
    REAL(KIND=dp)                            :: epsilon_rho, epsilon_tau, sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: params, scale
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rhoa, &
      e_laplace_rhoa_laplace_rhoa, e_laplace_rhoa_laplace_rhob, &
      e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, e_laplace_rhob, &
      e_laplace_rhob_laplace_rhob, e_laplace_rhob_tau_a, &
      e_laplace_rhob_tau_b, e_ndrho, e_ndrho_laplace_rhoa, &
      e_ndrho_laplace_rhob, e_ndrho_ndrho, e_ndrho_ndrhoa, e_ndrho_ndrhob, &
      e_ndrho_rhoa, e_ndrho_rhob, e_ndrho_tau_a, e_ndrho_tau_b, e_ndrhoa, &
      e_ndrhoa_laplace_rhoa, e_ndrhoa_laplace_rhob, e_ndrhoa_ndrhoa, &
      e_ndrhoa_ndrhob, e_ndrhoa_rhoa, e_ndrhoa_rhob, e_ndrhoa_tau_a, &
      e_ndrhoa_tau_b, e_ndrhob
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: e_ndrhob_laplace_rhoa, &
      e_ndrhob_laplace_rhob, e_ndrhob_ndrhob, e_ndrhob_rhoa, e_ndrhob_rhob, &
      e_ndrhob_tau_a, e_ndrhob_tau_b, e_rhoa, e_rhoa_laplace_rhoa, &
      e_rhoa_laplace_rhob, e_rhoa_rhoa, e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, &
      e_rhoa_rhob, e_rhoa_rhob_rhob, e_rhoa_tau_a, e_rhoa_tau_b, e_rhob, &
      e_rhob_laplace_rhoa, e_rhob_laplace_rhob, e_rhob_rhob, &
      e_rhob_rhob_rhob, e_rhob_tau_a, e_rhob_tau_b, e_tau_a, e_tau_a_tau_a, &
      e_tau_a_tau_b, e_tau_b, e_tau_b_tau_b, laplace_rhoa, laplace_rhob, &
      norm_drho, norm_drhoa, norm_drhob, rhoa, rhob, tau_a, tau_b
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL timeset(routineN,handle)

    failure = .FALSE.
    has_laplace = .FALSE.
    NULLIFY(bo, dummy)
    NULLIFY(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob, laplace_rhoa, &
            laplace_rhob, tau_a, tau_b)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    IF (.NOT. failure) THEN

       CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name)
       CALL section_vals_val_get(libxc_params,"scale",r_vals=scale)
       CALL section_vals_val_get(libxc_params,"parameters",r_vals=params)

       IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
         sc = 1.0_dp
       ELSE
         sc = SCALE(ifunc_name)
       ENDIF

       func_id = xc_libxc_wrap_functional_get_number(func_name(ifunc_name))
       CALL xc_libxc_wrap_functional_buggy(func_id, grad_deriv=grad_deriv)
!$OMP CRITICAL(libxc_init)
       CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_POLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

       CALL xc_rho_set_get(rho_set, can_return_null=.TRUE.,&
            rhoa=rhoa, rhob=rhob, norm_drho=norm_drho,&
            norm_drhoa=norm_drhoa, norm_drhob=norm_drhob,&
            laplace_rhoa=laplace_rhoa, laplace_rhob=laplace_rhob,&
            rho_cutoff=epsilon_rho, tau_cutoff=epsilon_tau,&
            tau_a=tau_a, tau_b=tau_b, local_bounds=bo)

       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
       ELSE
          dummy=> rhoa
       END IF

       ! due to assumed shape array usage in next routine
       IF (.NOT.ASSOCIATED(norm_drho)) norm_drho => dummy
       IF (.NOT.ASSOCIATED(norm_drhoa)) norm_drhoa => dummy
       IF (.NOT.ASSOCIATED(norm_drhob)) norm_drhob => dummy
       IF (.NOT.ASSOCIATED(tau_a)) tau_a => dummy
       IF (.NOT.ASSOCIATED(tau_b)) tau_b => dummy

       ! only some MGGA functionals really need the Laplacian,
       ! all others can work with rhoa (read-only) as dummy
       IF (ASSOCIATED(laplace_rhoa) .AND. ASSOCIATED(laplace_rhob)) has_laplace = .TRUE. 
       IF (.NOT. has_laplace) laplace_rhoa => dummy 
       IF (.NOT. has_laplace) laplace_rhob => dummy

       e_0 => dummy
       e_rhoa => dummy
       e_rhob => dummy
       e_ndrho => dummy
       e_ndrhoa => dummy
       e_ndrhob => dummy
       e_laplace_rhoa => dummy
       e_laplace_rhob => dummy
       e_tau_a => dummy
       e_tau_b => dummy
       e_rhoa_rhoa => dummy
       e_rhoa_rhob => dummy
       e_rhob_rhob => dummy
       e_ndrho_rhoa => dummy
       e_ndrho_rhob => dummy
       e_ndrhoa_rhoa => dummy
       e_ndrhoa_rhob => dummy
       e_ndrhob_rhoa => dummy
       e_ndrhob_rhob => dummy
       e_ndrho_ndrho => dummy
       e_ndrho_ndrhoa => dummy
       e_ndrho_ndrhob => dummy
       e_ndrhoa_ndrhoa => dummy
       e_ndrhoa_ndrhob => dummy
       e_ndrhob_ndrhob => dummy
       e_rhoa_laplace_rhoa => dummy
       e_rhoa_laplace_rhob => dummy
       e_rhob_laplace_rhoa => dummy
       e_rhob_laplace_rhob => dummy
       e_rhoa_tau_a => dummy
       e_rhoa_tau_b => dummy
       e_rhob_tau_a => dummy
       e_rhob_tau_b => dummy
       e_ndrho_laplace_rhoa => dummy
       e_ndrho_laplace_rhob => dummy
       e_ndrhoa_laplace_rhoa => dummy
       e_ndrhoa_laplace_rhob => dummy
       e_ndrhob_laplace_rhoa => dummy
       e_ndrhob_laplace_rhob => dummy
       e_ndrho_tau_a => dummy
       e_ndrho_tau_b => dummy
       e_ndrhoa_tau_a => dummy
       e_ndrhoa_tau_b => dummy
       e_ndrhob_tau_a => dummy
       e_ndrhob_tau_b => dummy
       e_laplace_rhoa_laplace_rhoa => dummy
       e_laplace_rhoa_laplace_rhob => dummy
       e_laplace_rhob_laplace_rhob => dummy
       e_laplace_rhoa_tau_a => dummy
       e_laplace_rhoa_tau_b => dummy
       e_laplace_rhob_tau_a => dummy
       e_laplace_rhob_tau_b => dummy
       e_tau_a_tau_a => dummy
       e_tau_a_tau_b => dummy
       e_tau_b_tau_b => dummy
       e_rhoa_rhoa_rhoa => dummy
       e_rhoa_rhoa_rhob => dummy
       e_rhoa_rhob_rhob => dummy
       e_rhob_rhob_rhob => dummy

       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE.)
          CALL xc_derivative_get(deriv,deriv_data=e_0)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob)
          CASE(XC_FAMILY_MGGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
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
             IF (has_laplace) THEN
                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)",&
                     allocate_deriv=.TRUE.)
                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa)
                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)",&
                     allocate_deriv=.TRUE.)
                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob)
             END IF
          CASE default
             CALL stop_program(routineN,moduleN,__LINE__,&
                TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
          END SELECT
       END IF
       IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhob)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(norm_drhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_ndrhob)
          CASE(XC_FAMILY_MGGA)
            ! not implemented ...
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives larger than 1 not implemented or checked",&
                error_level=cp_failure_level)

!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhoa)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhoa)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhoa)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhoa)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhoa)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhoa)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhoa)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhoa)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(norm_drhob)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_ndrhob)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(tau_a)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_tau_a)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_tau_b)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(tau_a)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhob_tau_a)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhob_tau_b)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(tau_a)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_tau_a)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_tau_b)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(tau_a)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_tau_a)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_tau_b)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(tau_a)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_tau_a)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_tau_b)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau_a)(tau_a)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_a_tau_a)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau_a)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_a_tau_b)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau_b)(tau_b)",&
!                  allocate_deriv=.TRUE.)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_b_tau_b)
!            IF (has_laplace) THEN
!                deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhoa_laplace_rhoa)
!                deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhoa_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(rhob)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhob_laplace_rhoa)
!                deriv => xc_dset_get_derivative(deriv_set,"(rhob)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhob_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrho_laplace_rhoa)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrho_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_laplace_rhoa)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_laplace_rhoa)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_laplace_rhoa)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob_laplace_rhob)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(tau_a)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_tau_a)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(tau_b)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_tau_b)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)(tau_a)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob_tau_a)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)(tau_b)",&
!                     allocate_deriv=.TRUE.)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob_tau_b)
!             END IF
          CASE default
             CALL stop_program(routineN,moduleN,__LINE__,&
                TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
          END SELECT
       END IF
       IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa_rhoa)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob_rhob)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
                  allocate_deriv=.TRUE.)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob_rhob)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA)
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives larger than 2 not implemented",&
                error_level=cp_failure_level)
          CASE default
             CALL stop_program(routineN,moduleN,__LINE__,&
                TRIM(func_name(ifunc_name))//": this XC_FAMILY is currently not supported.")
          END SELECT
       END IF
       IF (grad_deriv>=4.OR.grad_deriv<=-4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP,&
             message="derivatives larger than 3 not implemented",&
             error_level=cp_failure_level)
       END IF

       !$omp parallel default(none), &
       !$omp shared(rhoa,rhob,norm_drho,norm_drhoa,norm_drhob),&
       !$omp shared(laplace_rhoa,laplace_rhob,tau_a,tau_b),&
       !$omp shared(e_0,e_rhoa,e_rhob,e_ndrho,e_ndrhoa,e_ndrhob),&
       !$omp shared(e_laplace_rhoa,e_laplace_rhob,e_tau_a,e_tau_b),&
       !$omp shared(e_rhoa_rhoa,e_rhoa_rhob,e_rhob_rhob),&
       !$omp shared(e_ndrho_rhoa,e_ndrho_rhob),&
       !$omp shared(e_ndrhoa_rhoa,e_ndrhoa_rhob,e_ndrhob_rhoa,e_ndrhob_rhob),&
       !$omp shared(e_ndrho_ndrho,e_ndrho_ndrhoa,e_ndrho_ndrhob),&
       !$omp shared(e_ndrhoa_ndrhoa,e_ndrhoa_ndrhob,e_ndrhob_ndrhob),&
       !$omp shared(e_rhoa_laplace_rhoa,e_rhoa_laplace_rhob,e_rhob_laplace_rhoa,e_rhob_laplace_rhob),&
       !$omp shared(e_rhoa_tau_a,e_rhoa_tau_b,e_rhob_tau_a,e_rhob_tau_b),&
       !$omp shared(e_ndrho_laplace_rhoa,e_ndrho_laplace_rhob),&
       !$omp shared(e_ndrhoa_laplace_rhoa,e_ndrhoa_laplace_rhob,e_ndrhob_laplace_rhoa,e_ndrhob_laplace_rhob),&
       !$omp shared(e_ndrho_tau_a,e_ndrho_tau_b),&
       !$omp shared(e_ndrhoa_tau_a,e_ndrhoa_tau_b,e_ndrhob_tau_a,e_ndrhob_tau_b),&
       !$omp shared(e_laplace_rhoa_laplace_rhoa,e_laplace_rhoa_laplace_rhob,e_laplace_rhob_laplace_rhob),&
       !$omp shared(e_laplace_rhoa_tau_a,e_laplace_rhoa_tau_b,e_laplace_rhob_tau_a,e_laplace_rhob_tau_b),&
       !$omp shared(e_tau_a_tau_a,e_tau_a_tau_b,e_tau_b_tau_b),&
       !$omp shared(e_rhoa_rhoa_rhoa,e_rhoa_rhoa_rhob,e_rhoa_rhob_rhob,e_rhob_rhob_rhob),&
       !$omp shared(grad_deriv,npoints),&
       !$omp shared(epsilon_rho,epsilon_tau),&
       !$omp shared(func_name,ifunc_name,sc,params)

       CALL libxc_lsd_calc(rhoa=rhoa,rhob=rhob,norm_drho=norm_drho,&
          norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,laplace_rhoa=laplace_rhoa,&
          laplace_rhob=laplace_rhob,tau_a=tau_a,tau_b=tau_b,&
          e_0=e_0,e_rhoa=e_rhoa,e_rhob=e_rhob,e_ndrho=e_ndrho,&
          e_ndrhoa=e_ndrhoa,e_ndrhob=e_ndrhob,e_laplace_rhoa=e_laplace_rhoa,&
          e_laplace_rhob=e_laplace_rhob,e_tau_a=e_tau_a,e_tau_b=e_tau_b,&
          e_rhoa_rhoa=e_rhoa_rhoa,e_rhoa_rhob=e_rhoa_rhob,e_rhob_rhob=e_rhob_rhob,&
          e_ndrho_rhoa=e_ndrho_rhoa,e_ndrho_rhob=e_ndrho_rhob,&
          e_ndrhoa_rhoa=e_ndrhoa_rhoa,e_ndrhoa_rhob=e_ndrhoa_rhob,&
          e_ndrhob_rhoa=e_ndrhob_rhoa,e_ndrhob_rhob=e_ndrhob_rhob,&
          e_ndrho_ndrho=e_ndrho_ndrho,e_ndrho_ndrhoa=e_ndrho_ndrhoa,&
          e_ndrho_ndrhob=e_ndrho_ndrhob,e_ndrhoa_ndrhoa=e_ndrhoa_ndrhoa,&
          e_ndrhoa_ndrhob=e_ndrhoa_ndrhob,e_ndrhob_ndrhob=e_ndrhob_ndrhob,&
          e_rhoa_laplace_rhoa=e_rhoa_laplace_rhoa,&
          e_rhoa_laplace_rhob=e_rhoa_laplace_rhob,&
          e_rhob_laplace_rhoa=e_rhob_laplace_rhoa,&
          e_rhob_laplace_rhob=e_rhob_laplace_rhob,&
          e_rhoa_tau_a=e_rhoa_tau_a,e_rhoa_tau_b=e_rhoa_tau_b,&
          e_rhob_tau_a=e_rhob_tau_a,e_rhob_tau_b=e_rhob_tau_b,&
          e_ndrho_laplace_rhoa=e_ndrho_laplace_rhoa,&
          e_ndrho_laplace_rhob=e_ndrho_laplace_rhob,&
          e_ndrhoa_laplace_rhoa=e_ndrhoa_laplace_rhoa,&
          e_ndrhoa_laplace_rhob=e_ndrhoa_laplace_rhob,&
          e_ndrhob_laplace_rhoa=e_ndrhob_laplace_rhoa,&
          e_ndrhob_laplace_rhob=e_ndrhob_laplace_rhob,&
          e_ndrho_tau_a=e_ndrho_tau_a,e_ndrho_tau_b=e_ndrho_tau_b,&
          e_ndrhoa_tau_a=e_ndrhoa_tau_a,e_ndrhoa_tau_b=e_ndrhoa_tau_b,&
          e_ndrhob_tau_a=e_ndrhob_tau_a,e_ndrhob_tau_b=e_ndrhob_tau_b,&
          e_laplace_rhoa_laplace_rhoa=e_laplace_rhoa_laplace_rhoa,&
          e_laplace_rhoa_laplace_rhob=e_laplace_rhoa_laplace_rhob,&
          e_laplace_rhob_laplace_rhob=e_laplace_rhob_laplace_rhob,&
          e_laplace_rhoa_tau_a=e_laplace_rhoa_tau_a,&
          e_laplace_rhoa_tau_b=e_laplace_rhoa_tau_b,&
          e_laplace_rhob_tau_a=e_laplace_rhob_tau_a,&
          e_laplace_rhob_tau_b=e_laplace_rhob_tau_b,&
          e_tau_a_tau_a=e_tau_a_tau_a,&
          e_tau_a_tau_b=e_tau_a_tau_b,&
          e_tau_b_tau_b=e_tau_b_tau_b,&
          e_rhoa_rhoa_rhoa=e_rhoa_rhoa_rhoa,&
          e_rhoa_rhoa_rhob=e_rhoa_rhoa_rhob,&
          e_rhoa_rhob_rhob=e_rhoa_rhob_rhob,&
          e_rhob_rhob_rhob=e_rhob_rhob_rhob,&
          grad_deriv=grad_deriv,npoints=npoints,&
          epsilon_rho=epsilon_rho,&
          epsilon_tau=epsilon_tau,func_name=func_name(ifunc_name),&
          sc=sc,params=params)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy)
       ELSE
          NULLIFY(dummy)
       END IF

       CALL xc_f90_func_end(xc_func)

    END IF

    CALL timestop(handle)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif
  END SUBROUTINE libxc_lsd_eval

! *****************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rho density
!> \param norm_drho norm of the gradient of the density
!> \param laplace_rho laplacian of the density
!> \param tau kinetic-energy density
!> \param e_0 energy density
!> \param e_rho derivative of the energy density with respect to rho
!> \param e_ndrho derivative of the energy density with respect to ndrho
!> \param e_laplace_rho derivative of the energy density with respect to laplace_rho
!> \param e_tau derivative of the energy density with respect to tau
!> \param e_rho_rho derivative of the energy density with respect to rho_rho
!> \param e_ndrho_rho derivative of the energy density with respect to ndrho_rho
!> \param e_ndrho_ndrho derivative of the energy density with respect to ndrho_ndrho
!> \param e_rho_laplace_rho derivative of the energy density with respect to rho_laplace_rho
!> \param e_rho_tau derivative of the energy density with respect to rho_tau
!> \param e_ndrho_laplace_rho derivative of the energy density with respect to ndrho_laplace_rho
!> \param e_ndrho_tau derivative of the energy density with respect to ndrho_tau
!> \param e_laplace_rho_laplace_rho derivative of the energy density with respect to laplace_rho_laplace_rho
!> \param e_laplace_rho_tau derivative of the energy density with respect to laplace_rho_tau
!> \param e_tau_tau derivative of the energy density with respect to tau_tau
!> \param e_rho_rho_rho derivative of the energy density with respect to rho_rho_rho
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of points on the grid
!> \param epsilon_rho ...
!> \param epsilon_tau ...
!> \param func_name name of the functional
!> \param sc scaling factor
!> \param params parameters of the functional
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_calc(rho,norm_drho,laplace_rho,tau,&
          e_0,e_rho,e_ndrho,e_laplace_rho,e_tau,e_rho_rho,e_ndrho_rho,&
          e_ndrho_ndrho,e_rho_laplace_rho,e_rho_tau,e_ndrho_laplace_rho,&
          e_ndrho_tau,e_laplace_rho_laplace_rho,e_laplace_rho_tau,&
          e_tau_tau,e_rho_rho_rho,&
          grad_deriv,npoints,epsilon_rho,&
          epsilon_tau,func_name,sc,params)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho, laplace_rho, &
                                                tau
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rho, e_ndrho, &
      e_laplace_rho, e_tau, e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, &
      e_rho_laplace_rho, e_rho_tau, e_ndrho_laplace_rho, e_ndrho_tau, &
      e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_tau_tau, e_rho_rho_rho
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(KIND=dp), INTENT(in)                :: epsilon_rho, epsilon_tau
    CHARACTER(LEN=80), INTENT(IN)            :: func_name
    REAL(KIND=dp), INTENT(in)                :: sc
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_calc', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    INTEGER                                  :: func_id, ii
    LOGICAL                                  :: no_exc
    REAL(KIND=dp) :: exc, my_tau, v2lapl2, v2lapltau, v2rho2, v2rholapl, &
      v2rhosigma, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2tau2, &
      v3rho3, vlapl, vrho, vsigma, vtau
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    func_id = xc_libxc_wrap_functional_get_number(func_name)
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_UNPOLARIZED)
    CALL xc_libxc_wrap_functional_set_params(xc_func, xc_info, func_id, params, no_exc)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER
    SELECT CASE (xc_f90_info_family(xc_info))
    CASE(XC_FAMILY_LDA)
       IF (grad_deriv==0) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda_exc(xc_func,1,rho(ii),exc)
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-1) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda_vxc(xc_func,1,rho(ii),vrho)
                e_rho(ii) = e_rho(ii) + sc*vrho
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==1) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda_exc_vxc(xc_func,1,rho(ii),exc,vrho)
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-2) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda_fxc(xc_func,1,rho(ii),v2rho2)
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==2) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda_exc_vxc(xc_func,1,rho(ii),exc,vrho)
                CALL xc_f90_lda_fxc(xc_func,1,rho(ii),v2rho2)
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-3) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda_kxc(xc_func,1,rho(ii),v3rho3)
                e_rho_rho_rho(ii) = e_rho_rho_rho(ii) + sc*v3rho3
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==3) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_lda(xc_func,1,rho(ii),exc,vrho,v2rho2,v3rho3)
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
                e_rho_rho_rho(ii) = e_rho_rho_rho(ii) + sc*v3rho3
             END IF
          END DO
          !$end do
       END IF
    CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
       IF (grad_deriv==0) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_gga_exc(xc_func,1,rho(ii),norm_drho(ii)**2,exc)
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-1) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                CALL xc_f90_gga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,vrho,vsigma)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*vsigma*norm_drho(ii)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==1) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                IF (no_exc) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,vrho,vsigma)
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_gga_exc_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     exc,vrho,vsigma)
                END IF
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*vsigma*norm_drho(ii)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-2) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                IF (no_exc) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,vrho,vsigma)
                  CALL xc_f90_gga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     v2rho2,v2rhosigma,v2sigma2)
                ELSE
                  CALL xc_f90_gga_exc_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     exc,vrho,vsigma)
                  CALL xc_f90_gga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     v2rho2,v2rhosigma,v2sigma2)
                END IF
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
                e_ndrho_rho(ii) = e_ndrho_rho(ii) + sc*2.0_dp*v2rhosigma*norm_drho(ii)
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*2.0_dp*(2.0_dp*norm_drho(ii)**2*v2sigma2 + vsigma)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==2) THEN
          !$omp do
          DO ii=1,npoints
             IF (rho(ii) > epsilon_rho) THEN
                IF (no_exc) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,vrho,vsigma)
                  CALL xc_f90_gga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     v2rho2,v2rhosigma,v2sigma2)
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_gga_exc_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     exc,vrho,vsigma)
                  CALL xc_f90_gga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     v2rho2,v2rhosigma,v2sigma2)
                END IF
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*vsigma*norm_drho(ii)
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
                e_ndrho_rho(ii) = e_ndrho_rho(ii) + sc*2.0_dp*v2rhosigma*norm_drho(ii)
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*2.0_dp*(2.0_dp*norm_drho(ii)**2*v2sigma2 + vsigma)
             END IF
          END DO
          !$end do
       END IF
    CASE(XC_FAMILY_MGGA)
       IF (grad_deriv==0) THEN
          !$omp do
          DO ii=1,npoints
             IF ((rho(ii) > epsilon_rho) .AND. (tau(ii) > epsilon_tau)) THEN
                my_tau = MAX(tau(ii),norm_drho(ii)**2/(8.0_dp*rho(ii)))
                CALL xc_f90_mgga_exc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                   laplace_rho(ii),my_tau,exc)
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-1) THEN
          !$omp do
          DO ii=1,npoints
             IF ((rho(ii) > epsilon_rho) .AND. (tau(ii) > epsilon_tau)) THEN
                my_tau = MAX(tau(ii),norm_drho(ii)**2/(8.0_dp*rho(ii)))
                CALL xc_f90_mgga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                   laplace_rho(ii),my_tau,vrho,vsigma,vlapl,vtau)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*vsigma*norm_drho(ii)
                e_laplace_rho(ii) = e_laplace_rho(ii) + sc*vlapl
                e_tau(ii) = e_tau(ii) + sc*vtau
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==1) THEN
          !$omp do
          DO ii=1,npoints
             IF ((rho(ii) > epsilon_rho) .AND. (tau(ii) > epsilon_tau)) THEN
                my_tau = MAX(tau(ii),norm_drho(ii)**2/(8.0_dp*rho(ii)))
                IF (no_exc) THEN
                  CALL xc_f90_mgga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,vrho,vsigma,vlapl,vtau)
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_mgga_exc_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,exc,vrho,vsigma,vlapl,vtau)
                END IF
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*vsigma*norm_drho(ii)
                e_laplace_rho(ii) = e_laplace_rho(ii) + sc*vlapl
                e_tau(ii) = e_tau(ii) + sc*vtau
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-2) THEN
          !$omp do
          DO ii=1,npoints
             IF ((rho(ii) > epsilon_rho) .AND. (tau(ii) > epsilon_tau)) THEN
                my_tau = MAX(tau(ii),norm_drho(ii)**2/(8.0_dp*rho(ii)))
                IF (no_exc) THEN
                  CALL xc_f90_mgga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,vrho,vsigma,vlapl,vtau)
                  CALL xc_f90_mgga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,&
                     v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,&
                     v2rhotau,v2sigmalapl,v2sigmatau,v2lapltau)
                ELSE
                  CALL xc_f90_mgga(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,exc,vrho,vsigma,vlapl,vtau,&
                     v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,&
                     v2rhotau,v2sigmalapl,v2sigmatau,v2lapltau)
                END IF
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
                e_ndrho_rho(ii) = e_ndrho_rho(ii) + sc*2.0_dp*v2rhosigma*norm_drho(ii)
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*2.0_dp*(2.0_dp*norm_drho(ii)**2*v2sigma2 + vsigma)
                e_rho_laplace_rho(ii) = e_rho_laplace_rho(ii) + sc*v2rholapl
                e_rho_tau(ii) = e_rho_tau(ii) + sc*v2rhotau
                e_ndrho_laplace_rho(ii) = e_ndrho_laplace_rho(ii) + &
                   sc*2.0_dp*v2sigmalapl*norm_drho(ii)
                e_ndrho_tau(ii) = e_ndrho_tau(ii) + sc*2.0_dp*v2sigmatau*norm_drho(ii)
                e_laplace_rho_laplace_rho(ii) = e_laplace_rho_laplace_rho(ii) + sc*v2lapl2
                e_laplace_rho_tau(ii) = e_laplace_rho_tau(ii) + sc*v2lapltau
                e_tau_tau(ii) = e_tau_tau(ii) + sc*v2tau2
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==2) THEN
          !$omp do
          DO ii=1,npoints
             IF ((rho(ii) > epsilon_rho) .AND. (tau(ii) > epsilon_tau)) THEN
                my_tau = MAX(tau(ii),norm_drho(ii)**2/(8.0_dp*rho(ii)))
                IF (no_exc) THEN
                  CALL xc_f90_mgga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,vrho,vsigma,vlapl,vtau)
                  CALL xc_f90_mgga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,&
                     v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,&
                     v2rhotau,v2sigmalapl,v2sigmatau,v2lapltau)
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_mgga(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     laplace_rho(ii),my_tau,exc,vrho,vsigma,vlapl,vtau,&
                     v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,&
                     v2rhotau,v2sigmalapl,v2sigmatau,v2lapltau)
                END IF
                e_0(ii) = e_0(ii) + sc*exc*rho(ii)
                e_rho(ii) = e_rho(ii) + sc*vrho
                e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*vsigma*norm_drho(ii)
                e_laplace_rho(ii) = e_laplace_rho(ii) + sc*vlapl
                e_tau(ii) = e_tau(ii) + sc*vtau
                e_rho_rho(ii) = e_rho_rho(ii) + sc*v2rho2
                e_ndrho_rho(ii) = e_ndrho_rho(ii) + sc*2.0_dp*v2rhosigma*norm_drho(ii)
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*2.0_dp*(2.0_dp*norm_drho(ii)**2*v2sigma2 + vsigma)
                e_rho_laplace_rho(ii) = e_rho_laplace_rho(ii) + sc*v2rholapl
                e_rho_tau(ii) = e_rho_tau(ii) + sc*v2rhotau
                e_ndrho_laplace_rho(ii) = e_ndrho_laplace_rho(ii) + &
                   sc*2.0_dp*v2sigmalapl*norm_drho(ii)
                e_ndrho_tau(ii) = e_ndrho_tau(ii) + sc*2.0_dp*v2sigmatau*norm_drho(ii)
                e_laplace_rho_laplace_rho(ii) = e_laplace_rho_laplace_rho(ii) + sc*v2lapl2
                e_laplace_rho_tau(ii) = e_laplace_rho_tau(ii) + sc*v2lapltau
                e_tau_tau(ii) = e_tau_tau(ii) + sc*v2tau2
             END IF
          END DO
          !$end do
       END IF
    CASE default
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name)//": this XC_FAMILY is currently not supported.")
    END SELECT

    CALL xc_f90_func_end(xc_func)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif

  END SUBROUTINE libxc_lda_calc

! *****************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rhoa alpha density
!> \param rhob beta density
!> \param norm_drho ...
!> \param norm_drhoa norm of the gradient of the alpha density
!> \param norm_drhob norm of the gradient of the beta density
!> \param laplace_rhoa laplacian of the alpha density
!> \param laplace_rhob laplacian of the beta density
!> \param tau_a alpha kinetic-energy density
!> \param tau_b beta kinetic-energy density
!> \param e_0 energy density
!> \param e_rhoa derivative of the energy density with respect to rhoa
!> \param e_rhob derivative of the energy density with respect to rhob
!> \param e_ndrho derivative of the energy density with respect to ndrho
!> \param e_ndrhoa derivative of the energy density with respect to ndrhoa
!> \param e_ndrhob derivative of the energy density with respect to ndrhob
!> \param e_laplace_rhoa derivative of the energy density with respect to laplace_rhoa
!> \param e_laplace_rhob derivative of the energy density with respect to laplace_rhob
!> \param e_tau_a derivative of the energy density with respect to tau_a
!> \param e_tau_b derivative of the energy density with respect to tau_b
!> \param e_rhoa_rhoa derivative of the energy density with respect to rhoa_rhoa
!> \param e_rhoa_rhob derivative of the energy density with respect to rhoa_rhob
!> \param e_rhob_rhob derivative of the energy density with respect to rhob_rhob
!> \param e_ndrho_rhoa derivative of the energy density with respect to ndrho_rhoa
!> \param e_ndrho_rhob derivative of the energy density with respect to ndrho_rhob
!> \param e_ndrhoa_rhoa derivative of the energy density with respect to ndrhoa_rhoa
!> \param e_ndrhoa_rhob derivative of the energy density with respect to ndrhoa_rhob
!> \param e_ndrhob_rhoa derivative of the energy density with respect to ndrhob_rhoa
!> \param e_ndrhob_rhob derivative of the energy density with respect to ndrhob_rhob
!> \param e_ndrho_ndrho derivative of the energy density with respect to ndrho_ndrho
!> \param e_ndrho_ndrhoa derivative of the energy density with respect to ndrho_ndrhoa
!> \param e_ndrho_ndrhob derivative of the energy density with respect to ndrho_ndrhob
!> \param e_ndrhoa_ndrhoa derivative of the energy density with respect to ndrhoa_ndrhoa
!> \param e_ndrhoa_ndrhob derivative of the energy density with respect to ndrhoa_ndrhob
!> \param e_ndrhob_ndrhob derivative of the energy density with respect to ndrhob_ndrhob
!> \param e_rhoa_laplace_rhoa derivative of the energy density with respect to rhoa_laplace_rhoa
!> \param e_rhoa_laplace_rhob derivative of the energy density with respect to rhoa_laplace_rhob
!> \param e_rhob_laplace_rhoa derivative of the energy density with respect to rhob_laplace_rhoa
!> \param e_rhob_laplace_rhob derivative of the energy density with respect to rhob_laplace_rhob
!> \param e_rhoa_tau_a derivative of the energy density with respect to rhoa_tau_a
!> \param e_rhoa_tau_b derivative of the energy density with respect to rhoa_tau_b
!> \param e_rhob_tau_a derivative of the energy density with respect to rhob_tau_a
!> \param e_rhob_tau_b derivative of the energy density with respect to rhob_tau_b
!> \param e_ndrho_laplace_rhoa derivative of the energy density with respect to ndrho_laplace_rhoa
!> \param e_ndrho_laplace_rhob derivative of the energy density with respect to ndrho_laplace_rhob
!> \param e_ndrhoa_laplace_rhoa derivative of the energy density with respect to ndrhoa_laplace_rhoa
!> \param e_ndrhoa_laplace_rhob derivative of the energy density with respect to ndrhoa_laplace_rhob
!> \param e_ndrhob_laplace_rhoa derivative of the energy density with respect to ndrhob_laplace_rhoa
!> \param e_ndrhob_laplace_rhob derivative of the energy density with respect to ndrhob_laplace_rhob
!> \param e_ndrho_tau_a derivative of the energy density with respect to ndrho_tau_a
!> \param e_ndrho_tau_b derivative of the energy density with respect to ndrho_tau_b
!> \param e_ndrhoa_tau_a derivative of the energy density with respect to ndrhoa_tau_a
!> \param e_ndrhoa_tau_b derivative of the energy density with respect to ndrhoa_tau_b
!> \param e_ndrhob_tau_a derivative of the energy density with respect to ndrhob_tau_a
!> \param e_ndrhob_tau_b derivative of the energy density with respect to ndrhob_tau_b
!> \param e_laplace_rhoa_laplace_rhoa derivative of the energy density with respect to laplace_rhoa_laplace_rhoa
!> \param e_laplace_rhoa_laplace_rhob derivative of the energy density with respect to laplace_rhoa_laplace_rhob
!> \param e_laplace_rhob_laplace_rhob derivative of the energy density with respect to laplace_rhob_laplace_rhob
!> \param e_laplace_rhoa_tau_a derivative of the energy density with respect to laplace_rhoa_tau_a
!> \param e_laplace_rhoa_tau_b derivative of the energy density with respect to laplace_rhoa_tau_b
!> \param e_laplace_rhob_tau_a derivative of the energy density with respect to laplace_rhob_tau_a
!> \param e_laplace_rhob_tau_b derivative of the energy density with respect to laplace_rhob_tau_b
!> \param e_tau_a_tau_a derivative of the energy density with respect to tau_a_tau_a
!> \param e_tau_a_tau_b derivative of the energy density with respect to tau_a_tau_b
!> \param e_tau_b_tau_b derivative of the energy density with respect to tau_b_tau_b
!> \param e_rhoa_rhoa_rhoa derivative of the energy density with respect to rhoa_rhoa_rhoa
!> \param e_rhoa_rhoa_rhob derivative of the energy density with respect to rhoa_rhoa_rhob
!> \param e_rhoa_rhob_rhob derivative of the energy density with respect to rhoa_rhob_rhob
!> \param e_rhob_rhob_rhob derivative of the energy density with respect to rhob_rhob_rhob
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of points on the grid
!> \param epsilon_rho ...
!> \param epsilon_tau ...
!> \param func_name name of the functional
!> \param sc scaling factor
!> \param params parameters of the functional
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_calc(rhoa,rhob,norm_drho,norm_drhoa,&
          norm_drhob,laplace_rhoa,laplace_rhob,tau_a,tau_b,&
          e_0,e_rhoa,e_rhob,e_ndrho,e_ndrhoa,e_ndrhob,&
          e_laplace_rhoa,e_laplace_rhob,e_tau_a,e_tau_b,&
          e_rhoa_rhoa,e_rhoa_rhob,e_rhob_rhob,&
          e_ndrho_rhoa,e_ndrho_rhob,e_ndrhoa_rhoa,&
          e_ndrhoa_rhob,e_ndrhob_rhoa,e_ndrhob_rhob,&
          e_ndrho_ndrho,e_ndrho_ndrhoa,e_ndrho_ndrhob,&
          e_ndrhoa_ndrhoa,e_ndrhoa_ndrhob,e_ndrhob_ndrhob,&
          e_rhoa_laplace_rhoa,e_rhoa_laplace_rhob,&
          e_rhob_laplace_rhoa,e_rhob_laplace_rhob,&
          e_rhoa_tau_a,e_rhoa_tau_b,e_rhob_tau_a,e_rhob_tau_b,&
          e_ndrho_laplace_rhoa,e_ndrho_laplace_rhob,&
          e_ndrhoa_laplace_rhoa,e_ndrhoa_laplace_rhob,&
          e_ndrhob_laplace_rhoa,e_ndrhob_laplace_rhob,&
          e_ndrho_tau_a,e_ndrho_tau_b,&
          e_ndrhoa_tau_a,e_ndrhoa_tau_b,&
          e_ndrhob_tau_a,e_ndrhob_tau_b,&
          e_laplace_rhoa_laplace_rhoa,&
          e_laplace_rhoa_laplace_rhob,&
          e_laplace_rhob_laplace_rhob,&
          e_laplace_rhoa_tau_a,e_laplace_rhoa_tau_b,&
          e_laplace_rhob_tau_a,e_laplace_rhob_tau_b,&
          e_tau_a_tau_a,e_tau_a_tau_b,e_tau_b_tau_b,&
          e_rhoa_rhoa_rhoa,e_rhoa_rhoa_rhob,&
          e_rhoa_rhob_rhob,e_rhob_rhob_rhob,&
          grad_deriv,npoints,epsilon_rho,&
          epsilon_tau,func_name,sc,params)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rhoa, rhob, norm_drho, &
                                                norm_drhoa, norm_drhob, &
                                                laplace_rhoa, laplace_rhob, &
                                                tau_a, tau_b
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rhoa, e_rhob, &
      e_ndrho, e_ndrhoa, e_ndrhob, e_laplace_rhoa, e_laplace_rhob, e_tau_a, &
      e_tau_b, e_rhoa_rhoa, e_rhoa_rhob, e_rhob_rhob, e_ndrho_rhoa, &
      e_ndrho_rhob, e_ndrhoa_rhoa, e_ndrhoa_rhob, e_ndrhob_rhoa, &
      e_ndrhob_rhob, e_ndrho_ndrho, e_ndrho_ndrhoa, e_ndrho_ndrhob, &
      e_ndrhoa_ndrhoa, e_ndrhoa_ndrhob, e_ndrhob_ndrhob, e_rhoa_laplace_rhoa, &
      e_rhoa_laplace_rhob, e_rhob_laplace_rhoa, e_rhob_laplace_rhob, &
      e_rhoa_tau_a, e_rhoa_tau_b, e_rhob_tau_a, e_rhob_tau_b, &
      e_ndrho_laplace_rhoa, e_ndrho_laplace_rhob, e_ndrhoa_laplace_rhoa
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_ndrhoa_laplace_rhob, &
      e_ndrhob_laplace_rhoa, e_ndrhob_laplace_rhob, e_ndrho_tau_a, &
      e_ndrho_tau_b, e_ndrhoa_tau_a, e_ndrhoa_tau_b, e_ndrhob_tau_a, &
      e_ndrhob_tau_b, e_laplace_rhoa_laplace_rhoa, &
      e_laplace_rhoa_laplace_rhob, e_laplace_rhob_laplace_rhob, &
      e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, e_laplace_rhob_tau_a, &
      e_laplace_rhob_tau_b, e_tau_a_tau_a, e_tau_a_tau_b, e_tau_b_tau_b, &
      e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, e_rhoa_rhob_rhob, e_rhob_rhob_rhob
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(KIND=dp), INTENT(in)                :: epsilon_rho, epsilon_tau
    CHARACTER(LEN=80), INTENT(IN)            :: func_name
    REAL(KIND=dp), INTENT(in)                :: sc
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_calc', &
      routineP = moduleN//':'//routineN

#if defined (__LIBXC)
    INTEGER                                  :: func_id, ii
    LOGICAL                                  :: no_exc
    REAL(KIND=dp)                            :: exc, my_norm_drho, &
                                                my_norm_drhoa, my_norm_drhob, &
                                                my_rhoa, my_rhob, my_tau_a, &
                                                my_tau_b
    REAL(KIND=dp), DIMENSION(2, 1)           :: laplace_rhov, rhov, tauv, &
                                                vlapl, vrho, vtau
    REAL(KIND=dp), DIMENSION(3, 1)           :: sigmav, v2lapl2, v2rho2, &
                                                v2tau2, vsigma
    REAL(KIND=dp), DIMENSION(4, 1)           :: v2lapltau, v2rholapl, &
                                                v2rhotau, v3rho3
    REAL(KIND=dp), DIMENSION(6, 1)           :: v2rhosigma, v2sigma2, &
                                                v2sigmalapl, v2sigmatau
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

! these are just dummy variables, you need to use the correct size if working

    func_id = xc_libxc_wrap_functional_get_number(func_name)
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_POLARIZED)
    CALL xc_libxc_wrap_functional_set_params(xc_func, xc_info, func_id, params, no_exc)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

    SELECT CASE (xc_f90_info_family(xc_info))
    CASE(XC_FAMILY_LDA)
       IF (grad_deriv==0) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda_exc(xc_func,1,rhov(1,1),exc)
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-1) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda_vxc(xc_func,1,rhov(1,1),vrho(1,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==1) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda_exc_vxc(xc_func,1,rhov(1,1),exc,vrho(1,1))
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-2) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda_fxc(xc_func,1,rhov(1,1),v2rho2(1,1))
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==2) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda_exc_vxc(xc_func,1,rhov(1,1),exc,vrho(1,1))
                CALL xc_f90_lda_fxc(xc_func,1,rhov(1,1),v2rho2(1,1))
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-3) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda_kxc(xc_func,1,rhov(1,1),v3rho3(1,1))
                e_rhoa_rhoa_rhoa(ii) = e_rhoa_rhoa_rhoa(ii) + sc*v3rho3(1,1)
                e_rhoa_rhoa_rhob(ii) = e_rhoa_rhoa_rhob(ii) + sc*v3rho3(2,1)
                e_rhoa_rhob_rhob(ii) = e_rhoa_rhob_rhob(ii) + sc*v3rho3(3,1)
                e_rhob_rhob_rhob(ii) = e_rhob_rhob_rhob(ii) + sc*v3rho3(4,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==3) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                CALL xc_f90_lda(xc_func,1,rhov(1,1),exc,vrho(1,1),v2rho2(1,1),v3rho3(1,1))
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
                e_rhoa_rhoa_rhoa(ii) = e_rhoa_rhoa_rhoa(ii) + sc*v3rho3(1,1)
                e_rhoa_rhoa_rhob(ii) = e_rhoa_rhoa_rhob(ii) + sc*v3rho3(2,1)
                e_rhoa_rhob_rhob(ii) = e_rhoa_rhob_rhob(ii) + sc*v3rho3(3,1)
                e_rhob_rhob_rhob(ii) = e_rhob_rhob_rhob(ii) + sc*v3rho3(4,1)
             END IF
          END DO
          !$end do
       END IF
    CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
       IF (grad_deriv==0) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                CALL xc_f90_gga_exc(xc_func,1,rhov(1,1),sigmav(1,1),exc)
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-1) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                CALL xc_f90_gga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),vrho(1,1),vsigma(1,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_ndrho(ii) = e_ndrho(ii) + sc*vsigma(2,1)*my_norm_drho
                e_ndrhoa(ii) = e_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1))*my_norm_drhoa
                e_ndrhob(ii) = e_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1))*my_norm_drhob
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==1) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                IF (no_exc) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),vrho(1,1),vsigma(1,1))
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_gga_exc_vxc(xc_func,1,rhov(1,1),sigmav(1,1),exc,vrho(1,1),vsigma(1,1))
                END IF
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_ndrho(ii) = e_ndrho(ii) + sc*vsigma(2,1)*my_norm_drho
                e_ndrhoa(ii) = e_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1))*my_norm_drhoa
                e_ndrhob(ii) = e_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1))*my_norm_drhob
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-2) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                IF (no_exc) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),vrho(1,1),vsigma(1,1))
                  CALL xc_f90_gga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
                ELSE
                  CALL xc_f90_gga_exc_vxc(xc_func,1,rhov(1,1),sigmav(1,1),exc,vrho(1,1),vsigma(1,1))
                  CALL xc_f90_gga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
                END IF
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
                e_ndrho_rhoa(ii) = e_ndrho_rhoa(ii) + sc*v2rhosigma(2,1)*my_norm_drho
                e_ndrho_rhob(ii) = e_ndrho_rhob(ii) + sc*v2rhosigma(5,1)*my_norm_drho
                e_ndrhoa_rhoa(ii) = e_ndrhoa_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(1,1) - v2rhosigma(2,1))*my_norm_drhoa
                e_ndrhoa_rhob(ii) = e_ndrhoa_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(4,1) - v2rhosigma(5,1))*my_norm_drhoa
                e_ndrhob_rhoa(ii) = e_ndrhob_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(3,1) - v2rhosigma(2,1))*my_norm_drhob
                e_ndrhob_rhob(ii) = e_ndrhob_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(6,1) - v2rhosigma(5,1))*my_norm_drhob
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*(vsigma(2,1) + my_norm_drho**2*v2sigma2(4,1))
                e_ndrho_ndrhoa(ii) = e_ndrho_ndrhoa(ii) + &
                   sc*(2.0_dp*v2sigma2(2,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhoa
                e_ndrho_ndrhob(ii) = e_ndrho_ndrhob(ii) + &
                   sc*(2.0_dp*v2sigma2(5,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhob
                e_ndrhoa_ndrhoa(ii) = e_ndrhoa_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1) + my_norm_drhoa**2*( &
                   4.0_dp*v2sigma2(1,1) - 4.0_dp*v2sigma2(2,1) + v2sigma2(4,1)))
                e_ndrhoa_ndrhob(ii) = e_ndrhoa_ndrhob(ii) + &
                   sc*(4.0_dp*v2sigma2(3,1) - 2.0_dp*v2sigma2(2,1) - &
                   2.0_dp*v2sigma2(5,1) + v2sigma2(4,1))*my_norm_drhoa*my_norm_drhob
                e_ndrhob_ndrhob(ii) = e_ndrhob_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1) + my_norm_drhob**2*( &
                   4.0_dp*v2sigma2(6,1) - 4.0_dp*v2sigma2(5,1) + v2sigma2(4,1)))
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==2) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             IF ((my_rhoa+my_rhob) > epsilon_rho) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                IF (no_exc) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),vrho(1,1),vsigma(1,1))
                  CALL xc_f90_gga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_gga_exc_vxc(xc_func,1,rhov(1,1),sigmav(1,1),exc,vrho(1,1),vsigma(1,1))
                  CALL xc_f90_gga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
                END IF
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_ndrho(ii) = e_ndrho(ii) + sc*vsigma(2,1)*my_norm_drho
                e_ndrhoa(ii) = e_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1))*my_norm_drhoa
                e_ndrhob(ii) = e_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1))*my_norm_drhob
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
                e_ndrho_rhoa(ii) = e_ndrho_rhoa(ii) + sc*v2rhosigma(2,1)*my_norm_drho
                e_ndrho_rhob(ii) = e_ndrho_rhob(ii) + sc*v2rhosigma(5,1)*my_norm_drho
                e_ndrhoa_rhoa(ii) = e_ndrhoa_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(1,1) - v2rhosigma(2,1))*my_norm_drhoa
                e_ndrhoa_rhob(ii) = e_ndrhoa_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(4,1) - v2rhosigma(5,1))*my_norm_drhoa
                e_ndrhob_rhoa(ii) = e_ndrhob_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(3,1) - v2rhosigma(2,1))*my_norm_drhob
                e_ndrhob_rhob(ii) = e_ndrhob_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(6,1) - v2rhosigma(5,1))*my_norm_drhob
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*(vsigma(2,1) + my_norm_drho**2*v2sigma2(4,1))
                e_ndrho_ndrhoa(ii) = e_ndrho_ndrhoa(ii) + &
                   sc*(2.0_dp*v2sigma2(2,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhoa
                e_ndrho_ndrhob(ii) = e_ndrho_ndrhob(ii) + &
                   sc*(2.0_dp*v2sigma2(5,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhob
                e_ndrhoa_ndrhoa(ii) = e_ndrhoa_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1) + my_norm_drhoa**2*( &
                   4.0_dp*v2sigma2(1,1) - 4.0_dp*v2sigma2(2,1) + v2sigma2(4,1)))
                e_ndrhoa_ndrhob(ii) = e_ndrhoa_ndrhob(ii) + &
                   sc*(4.0_dp*v2sigma2(3,1) - 2.0_dp*v2sigma2(2,1) - &
                   2.0_dp*v2sigma2(5,1) + v2sigma2(4,1))*my_norm_drhoa*my_norm_drhob
                e_ndrhob_ndrhob(ii) = e_ndrhob_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1) + my_norm_drhob**2*( &
                   4.0_dp*v2sigma2(6,1) - 4.0_dp*v2sigma2(5,1) + v2sigma2(4,1)))
             END IF
          END DO
          !$end do
       END IF
    CASE(XC_FAMILY_MGGA)
       IF (grad_deriv==0) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             my_tau_a = MAX(tau_a(ii),0.0_dp)
             my_tau_b = MAX(tau_b(ii),0.0_dp)
             IF (((my_rhoa+my_rhob) > epsilon_rho) .AND. ((my_tau_a+my_tau_b) > epsilon_tau)) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                laplace_rhov(1,1) = laplace_rhoa(ii)
                laplace_rhov(2,1) = laplace_rhob(ii)
                tauv(1,1) = MAX(my_tau_a,EPSILON(0.0_dp)*1.e4_dp)
                tauv(2,1) = MAX(my_tau_b,EPSILON(0.0_dp)*1.e4_dp)
                tauv(1,1) = MAX(tauv(1,1),sigmav(1,1)/(8.0_dp*rhov(1,1)))
                tauv(2,1) = MAX(tauv(2,1),sigmav(3,1)/(8.0_dp*rhov(2,1)))
                CALL xc_f90_mgga_exc(xc_func,1,rhov(1,1),sigmav(1,1),&
                   laplace_rhov(1,1),tauv(1,1),exc)
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-1) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             my_tau_a = MAX(tau_a(ii),0.0_dp)
             my_tau_b = MAX(tau_b(ii),0.0_dp)
             IF (((my_rhoa+my_rhob) > epsilon_rho) .AND. ((my_tau_a+my_tau_b) > epsilon_tau)) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                laplace_rhov(1,1) = laplace_rhoa(ii)
                laplace_rhov(2,1) = laplace_rhob(ii)
                tauv(1,1) = MAX(my_tau_a,EPSILON(0.0_dp)*1.e4_dp)
                tauv(2,1) = MAX(my_tau_b,EPSILON(0.0_dp)*1.e4_dp)
                tauv(1,1) = MAX(tauv(1,1),sigmav(1,1)/(8.0_dp*rhov(1,1)))
                tauv(2,1) = MAX(tauv(2,1),sigmav(3,1)/(8.0_dp*rhov(2,1)))
                CALL xc_f90_mgga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                   laplace_rhov(1,1),tauv(1,1),vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_ndrho(ii) = e_ndrho(ii) + sc*vsigma(2,1)*my_norm_drho
                e_ndrhoa(ii) = e_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1))*my_norm_drhoa
                e_ndrhob(ii) = e_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1))*my_norm_drhob
                e_laplace_rhoa(ii) = e_laplace_rhoa(ii) + sc*vlapl(1,1)
                e_laplace_rhob(ii) = e_laplace_rhob(ii) + sc*vlapl(2,1)
                e_tau_a(ii) = e_tau_a(ii) + sc*vtau(1,1)
                e_tau_b(ii) = e_tau_b(ii) + sc*vtau(2,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==1) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             my_tau_a = MAX(tau_a(ii),0.0_dp)
             my_tau_b = MAX(tau_b(ii),0.0_dp)
             IF (((my_rhoa+my_rhob) > epsilon_rho) .AND. ((my_tau_a+my_tau_b) > epsilon_tau)) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                laplace_rhov(1,1) = laplace_rhoa(ii)
                laplace_rhov(2,1) = laplace_rhob(ii)
                tauv(1,1) = MAX(my_tau_a,EPSILON(0.0_dp)*1.e4_dp)
                tauv(2,1) = MAX(my_tau_b,EPSILON(0.0_dp)*1.e4_dp)
                tauv(1,1) = MAX(tauv(1,1),sigmav(1,1)/(8.0_dp*rhov(1,1)))
                tauv(2,1) = MAX(tauv(2,1),sigmav(3,1)/(8.0_dp*rhov(2,1)))
                IF (no_exc) THEN
                  CALL xc_f90_mgga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),vrho(1,1),vsigma(1,1),&
                     vlapl(1,1),vtau(1,1))
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_mgga_exc_vxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),exc,&
                     vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
                END IF
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_ndrho(ii) = e_ndrho(ii) + sc*vsigma(2,1)*my_norm_drho
                e_ndrhoa(ii) = e_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1))*my_norm_drhoa
                e_ndrhob(ii) = e_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1))*my_norm_drhob
                e_laplace_rhoa(ii) = e_laplace_rhoa(ii) + sc*vlapl(1,1)
                e_laplace_rhob(ii) = e_laplace_rhob(ii) + sc*vlapl(2,1)
                e_tau_a(ii) = e_tau_a(ii) + sc*vtau(1,1)
                e_tau_b(ii) = e_tau_b(ii) + sc*vtau(2,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==-2) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             my_tau_a = MAX(tau_a(ii),0.0_dp)
             my_tau_b = MAX(tau_b(ii),0.0_dp)
             IF (((my_rhoa+my_rhob) > epsilon_rho) .AND. ((my_tau_a+my_tau_b) > epsilon_tau)) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                laplace_rhov(1,1) = laplace_rhoa(ii)
                laplace_rhov(2,1) = laplace_rhob(ii)
                tauv(1,1) = MAX(my_tau_a,EPSILON(0.0_dp)*1.e4_dp)
                tauv(2,1) = MAX(my_tau_b,EPSILON(0.0_dp)*1.e4_dp)
                tauv(1,1) = MAX(tauv(1,1),sigmav(1,1)/(8.0_dp*rhov(1,1)))
                tauv(2,1) = MAX(tauv(2,1),sigmav(3,1)/(8.0_dp*rhov(2,1)))
                IF (no_exc) THEN
                  CALL xc_f90_mgga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),vrho(1,1),vsigma(1,1),&
                     vlapl(1,1),vtau(1,1))
                  CALL xc_f90_mgga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),&
                     v2rho2(1,1),v2sigma2(1,1),v2lapl2(1,1),v2tau2(1,1),&
                     v2rhosigma(1,1),v2rholapl(1,1),v2rhotau(1,1),&
                     v2sigmalapl(1,1),v2sigmatau(1,1),v2lapltau(1,1))
                ELSE
                  CALL xc_f90_mgga(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),exc,vrho(1,1),vsigma(1,1),&
                     vlapl(1,1),vtau(1,1),v2rho2(1,1),v2sigma2(1,1),&
                     v2lapl2(1,1),v2tau2(1,1),v2rhosigma(1,1),v2rholapl(1,1),&
                     v2rhotau(1,1),v2sigmalapl(1,1),v2sigmatau(1,1),v2lapltau(1,1))
                END IF
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
                e_ndrho_rhoa(ii) = e_ndrho_rhoa(ii) + sc*v2rhosigma(2,1)*my_norm_drho
                e_ndrho_rhob(ii) = e_ndrho_rhob(ii) + sc*v2rhosigma(5,1)*my_norm_drho
                e_ndrhoa_rhoa(ii) = e_ndrhoa_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(1,1) - v2rhosigma(2,1))*my_norm_drhoa
                e_ndrhoa_rhob(ii) = e_ndrhoa_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(4,1) - v2rhosigma(5,1))*my_norm_drhoa
                e_ndrhob_rhoa(ii) = e_ndrhob_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(3,1) - v2rhosigma(2,1))*my_norm_drhob
                e_ndrhob_rhob(ii) = e_ndrhob_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(6,1) - v2rhosigma(5,1))*my_norm_drhob
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*(vsigma(2,1) + my_norm_drho**2*v2sigma2(4,1))
                e_ndrho_ndrhoa(ii) = e_ndrho_ndrhoa(ii) + &
                   sc*(2.0_dp*v2sigma2(2,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhoa
                e_ndrho_ndrhob(ii) = e_ndrho_ndrhob(ii) + &
                   sc*(2.0_dp*v2sigma2(5,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhob
                e_ndrhoa_ndrhoa(ii) = e_ndrhoa_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1) + my_norm_drhoa**2*( &
                   4.0_dp*v2sigma2(1,1) - 4.0_dp*v2sigma2(2,1) + v2sigma2(4,1)))
                e_ndrhoa_ndrhob(ii) = e_ndrhoa_ndrhob(ii) + &
                   sc*(4.0_dp*v2sigma2(3,1) - 2.0_dp*v2sigma2(2,1) - &
                   2.0_dp*v2sigma2(5,1) + v2sigma2(4,1))*my_norm_drhoa*my_norm_drhob
                e_ndrhob_ndrhob(ii) = e_ndrhob_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1) + my_norm_drhob**2*( &
                   4.0_dp*v2sigma2(6,1) - 4.0_dp*v2sigma2(5,1) + v2sigma2(4,1)))
                e_rhoa_laplace_rhoa(ii) = e_rhoa_laplace_rhoa(ii) + sc*v2rholapl(1,1)
                e_rhoa_laplace_rhob(ii) = e_rhoa_laplace_rhob(ii) + sc*v2rholapl(2,1)
                e_rhob_laplace_rhoa(ii) = e_rhob_laplace_rhoa(ii) + sc*v2rholapl(3,1)
                e_rhob_laplace_rhob(ii) = e_rhob_laplace_rhob(ii) + sc*v2rholapl(4,1)
                e_rhoa_tau_a(ii) = e_rhoa_tau_a(ii) + sc*v2rhotau(1,1)
                e_rhoa_tau_b(ii) = e_rhoa_tau_b(ii) + sc*v2rhotau(2,1)
                e_rhob_tau_a(ii) = e_rhob_tau_a(ii) + sc*v2rhotau(3,1)
                e_rhob_tau_b(ii) = e_rhob_tau_b(ii) + sc*v2rhotau(4,1)
                e_ndrho_laplace_rhoa(ii) = e_ndrho_laplace_rhoa(ii) + sc*v2sigmalapl(3,1)*my_norm_drho
                e_ndrho_laplace_rhob(ii) = e_ndrho_laplace_rhob(ii) + sc*v2sigmalapl(4,1)*my_norm_drho
                e_ndrhoa_laplace_rhoa(ii) = e_ndrhoa_laplace_rhoa(ii) + &
                   sc*(2.0_dp*v2sigmalapl(1,1) - v2sigmalapl(3,1))*my_norm_drhoa
                e_ndrhoa_laplace_rhob(ii) = e_ndrhoa_laplace_rhob(ii) + &
                   sc*(2.0_dp*v2sigmalapl(2,1) - v2sigmalapl(4,1))*my_norm_drhoa
                e_ndrhob_laplace_rhoa(ii) = e_ndrhob_laplace_rhoa(ii) + &
                   sc*(2.0_dp*v2sigmalapl(5,1) - v2sigmalapl(3,1))*my_norm_drhob
                e_ndrhob_laplace_rhob(ii) = e_ndrhob_laplace_rhob(ii) + &
                   sc*(2.0_dp*v2sigmalapl(6,1) - v2sigmalapl(4,1))*my_norm_drhob
                e_ndrho_tau_a(ii) = e_ndrho_tau_a(ii) + sc*v2sigmatau(3,1)*my_norm_drho
                e_ndrho_tau_b(ii) = e_ndrho_tau_b(ii) + sc*v2sigmatau(4,1)*my_norm_drho
                e_ndrhoa_tau_a(ii) = e_ndrhoa_tau_a(ii) + &
                   sc*(2.0_dp*v2sigmatau(1,1) - v2sigmatau(3,1))*my_norm_drhoa
                e_ndrhoa_tau_b(ii) = e_ndrhoa_tau_b(ii) + &
                   sc*(2.0_dp*v2sigmatau(2,1) - v2sigmatau(4,1))*my_norm_drhoa
                e_ndrhob_tau_a(ii) = e_ndrhob_tau_a(ii) + &
                   sc*(2.0_dp*v2sigmatau(5,1) - v2sigmatau(3,1))*my_norm_drhob
                e_ndrhob_tau_b(ii) = e_ndrhob_tau_b(ii) + &
                   sc*(2.0_dp*v2sigmatau(6,1) - v2sigmatau(4,1))*my_norm_drhob
                e_laplace_rhoa_laplace_rhoa(ii) = e_laplace_rhoa_laplace_rhoa(ii) + sc*v2lapl2(1,1)
                e_laplace_rhoa_laplace_rhob(ii) = e_laplace_rhoa_laplace_rhob(ii) + sc*v2lapl2(2,1)
                e_laplace_rhob_laplace_rhob(ii) = e_laplace_rhob_laplace_rhob(ii) + sc*v2lapl2(3,1)
                e_laplace_rhoa_tau_a(ii) = e_laplace_rhoa_tau_a(ii) + sc*v2lapltau(1,1)
                e_laplace_rhoa_tau_b(ii) = e_laplace_rhoa_tau_b(ii) + sc*v2lapltau(2,1)
                e_laplace_rhob_tau_a(ii) = e_laplace_rhob_tau_a(ii) + sc*v2lapltau(3,1)
                e_laplace_rhob_tau_b(ii) = e_laplace_rhob_tau_b(ii) + sc*v2lapltau(4,1)
                e_tau_a_tau_a(ii) = e_tau_a_tau_a(ii) + sc*v2tau2(1,1)
                e_tau_a_tau_b(ii) = e_tau_a_tau_b(ii) + sc*v2tau2(2,1)
                e_tau_b_tau_b(ii) = e_tau_b_tau_b(ii) + sc*v2tau2(3,1)
             END IF
          END DO
          !$end do
       ELSE IF (grad_deriv==2) THEN
          !$omp do
          DO ii=1,npoints
             my_rhoa = MAX(rhoa(ii),0.0_dp)
             my_rhob = MAX(rhob(ii),0.0_dp)
             my_tau_a = MAX(tau_a(ii),0.0_dp)
             my_tau_b = MAX(tau_b(ii),0.0_dp)
             IF (((my_rhoa+my_rhob) > epsilon_rho) .AND. ((my_tau_a+my_tau_b) > epsilon_tau)) THEN
                rhov(1,1) = MAX(my_rhoa,EPSILON(0.0_dp)*1.e4_dp)
                rhov(2,1) = MAX(my_rhob,EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhoa = MAX(norm_drhoa(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drhob = MAX(norm_drhob(ii),EPSILON(0.0_dp)*1.e4_dp)
                my_norm_drho = MAX(norm_drho(ii),EPSILON(0.0_dp)*1.e4_dp)
                sigmav(1,1) = my_norm_drhoa**2
                sigmav(3,1) = my_norm_drhob**2
                sigmav(2,1) = 0.5_dp*(my_norm_drho**2 - sigmav(1,1) - sigmav(3,1))
                laplace_rhov(1,1) = laplace_rhoa(ii)
                laplace_rhov(2,1) = laplace_rhob(ii)
                tauv(1,1) = MAX(my_tau_a,EPSILON(0.0_dp)*1.e4_dp)
                tauv(2,1) = MAX(my_tau_b,EPSILON(0.0_dp)*1.e4_dp)
                tauv(1,1) = MAX(tauv(1,1),sigmav(1,1)/(8.0_dp*rhov(1,1)))
                tauv(2,1) = MAX(tauv(2,1),sigmav(3,1)/(8.0_dp*rhov(2,1)))
                IF (no_exc) THEN
                  CALL xc_f90_mgga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),vrho(1,1),vsigma(1,1),&
                     vlapl(1,1),vtau(1,1))
                  CALL xc_f90_mgga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),&
                     v2rho2(1,1),v2sigma2(1,1),v2lapl2(1,1),v2tau2(1,1),&
                     v2rhosigma(1,1),v2rholapl(1,1),v2rhotau(1,1),&
                     v2sigmalapl(1,1),v2sigmatau(1,1),v2lapltau(1,1))
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_mgga(xc_func,1,rhov(1,1),sigmav(1,1),&
                     laplace_rhov(1,1),tauv(1,1),exc,vrho(1,1),vsigma(1,1),&
                     vlapl(1,1),vtau(1,1),v2rho2(1,1),v2sigma2(1,1),&
                     v2lapl2(1,1),v2tau2(1,1),v2rhosigma(1,1),v2rholapl(1,1),&
                     v2rhotau(1,1),v2sigmalapl(1,1),v2sigmatau(1,1),v2lapltau(1,1))
                END IF
                e_0(ii) = e_0(ii) + sc*exc*(rhov(1,1) + rhov(2,1))
                e_rhoa(ii) = e_rhoa(ii) + sc*vrho(1,1)
                e_rhob(ii) = e_rhob(ii) + sc*vrho(2,1)
                e_ndrho(ii) = e_ndrho(ii) + sc*vsigma(2,1)*my_norm_drho
                e_ndrhoa(ii) = e_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1))*my_norm_drhoa
                e_ndrhob(ii) = e_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1))*my_norm_drhob
                e_laplace_rhoa(ii) = e_laplace_rhoa(ii) + sc*vlapl(1,1)
                e_laplace_rhob(ii) = e_laplace_rhob(ii) + sc*vlapl(2,1)
                e_tau_a(ii) = e_tau_a(ii) + sc*vtau(1,1)
                e_tau_b(ii) = e_tau_b(ii) + sc*vtau(2,1)
                e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*v2rho2(1,1)
                e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*v2rho2(2,1)
                e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*v2rho2(3,1)
                e_ndrho_rhoa(ii) = e_ndrho_rhoa(ii) + sc*v2rhosigma(2,1)*my_norm_drho
                e_ndrho_rhob(ii) = e_ndrho_rhob(ii) + sc*v2rhosigma(5,1)*my_norm_drho
                e_ndrhoa_rhoa(ii) = e_ndrhoa_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(1,1) - v2rhosigma(2,1))*my_norm_drhoa
                e_ndrhoa_rhob(ii) = e_ndrhoa_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(4,1) - v2rhosigma(5,1))*my_norm_drhoa
                e_ndrhob_rhoa(ii) = e_ndrhob_rhoa(ii) + &
                   sc*(2.0_dp*v2rhosigma(3,1) - v2rhosigma(2,1))*my_norm_drhob
                e_ndrhob_rhob(ii) = e_ndrhob_rhob(ii) + &
                   sc*(2.0_dp*v2rhosigma(6,1) - v2rhosigma(5,1))*my_norm_drhob
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                   sc*(vsigma(2,1) + my_norm_drho**2*v2sigma2(4,1))
                e_ndrho_ndrhoa(ii) = e_ndrho_ndrhoa(ii) + &
                   sc*(2.0_dp*v2sigma2(2,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhoa
                e_ndrho_ndrhob(ii) = e_ndrho_ndrhob(ii) + &
                   sc*(2.0_dp*v2sigma2(5,1) - v2sigma2(4,1))*my_norm_drho*my_norm_drhob
                e_ndrhoa_ndrhoa(ii) = e_ndrhoa_ndrhoa(ii) + &
                   sc*(2.0_dp*vsigma(1,1) - vsigma(2,1) + my_norm_drhoa**2*( &
                   4.0_dp*v2sigma2(1,1) - 4.0_dp*v2sigma2(2,1) + v2sigma2(4,1)))
                e_ndrhoa_ndrhob(ii) = e_ndrhoa_ndrhob(ii) + &
                   sc*(4.0_dp*v2sigma2(3,1) - 2.0_dp*v2sigma2(2,1) - &
                   2.0_dp*v2sigma2(5,1) + v2sigma2(4,1))*my_norm_drhoa*my_norm_drhob
                e_ndrhob_ndrhob(ii) = e_ndrhob_ndrhob(ii) + &
                   sc*(2.0_dp*vsigma(3,1) - vsigma(2,1) + my_norm_drhob**2*( &
                   4.0_dp*v2sigma2(6,1) - 4.0_dp*v2sigma2(5,1) + v2sigma2(4,1)))
                e_rhoa_laplace_rhoa(ii) = e_rhoa_laplace_rhoa(ii) + sc*v2rholapl(1,1)
                e_rhoa_laplace_rhob(ii) = e_rhoa_laplace_rhob(ii) + sc*v2rholapl(2,1)
                e_rhob_laplace_rhoa(ii) = e_rhob_laplace_rhoa(ii) + sc*v2rholapl(3,1)
                e_rhob_laplace_rhob(ii) = e_rhob_laplace_rhob(ii) + sc*v2rholapl(4,1)
                e_rhoa_tau_a(ii) = e_rhoa_tau_a(ii) + sc*v2rhotau(1,1)
                e_rhoa_tau_b(ii) = e_rhoa_tau_b(ii) + sc*v2rhotau(2,1)
                e_rhob_tau_a(ii) = e_rhob_tau_a(ii) + sc*v2rhotau(3,1)
                e_rhob_tau_b(ii) = e_rhob_tau_b(ii) + sc*v2rhotau(4,1)
                e_ndrho_laplace_rhoa(ii) = e_ndrho_laplace_rhoa(ii) + sc*v2sigmalapl(3,1)*my_norm_drho
                e_ndrho_laplace_rhob(ii) = e_ndrho_laplace_rhob(ii) + sc*v2sigmalapl(4,1)*my_norm_drho
                e_ndrhoa_laplace_rhoa(ii) = e_ndrhoa_laplace_rhoa(ii) + &
                   sc*(2.0_dp*v2sigmalapl(1,1) - v2sigmalapl(3,1))*my_norm_drhoa
                e_ndrhoa_laplace_rhob(ii) = e_ndrhoa_laplace_rhob(ii) + &
                   sc*(2.0_dp*v2sigmalapl(2,1) - v2sigmalapl(4,1))*my_norm_drhoa
                e_ndrhob_laplace_rhoa(ii) = e_ndrhob_laplace_rhoa(ii) + &
                   sc*(2.0_dp*v2sigmalapl(5,1) - v2sigmalapl(3,1))*my_norm_drhob
                e_ndrhob_laplace_rhob(ii) = e_ndrhob_laplace_rhob(ii) + &
                   sc*(2.0_dp*v2sigmalapl(6,1) - v2sigmalapl(4,1))*my_norm_drhob
                e_ndrho_tau_a(ii) = e_ndrho_tau_a(ii) + sc*v2sigmatau(3,1)*my_norm_drho
                e_ndrho_tau_b(ii) = e_ndrho_tau_b(ii) + sc*v2sigmatau(4,1)*my_norm_drho
                e_ndrhoa_tau_a(ii) = e_ndrhoa_tau_a(ii) + &
                   sc*(2.0_dp*v2sigmatau(1,1) - v2sigmatau(3,1))*my_norm_drhoa
                e_ndrhoa_tau_b(ii) = e_ndrhoa_tau_b(ii) + &
                   sc*(2.0_dp*v2sigmatau(2,1) - v2sigmatau(4,1))*my_norm_drhoa
                e_ndrhob_tau_a(ii) = e_ndrhob_tau_a(ii) + &
                   sc*(2.0_dp*v2sigmatau(5,1) - v2sigmatau(3,1))*my_norm_drhob
                e_ndrhob_tau_b(ii) = e_ndrhob_tau_b(ii) + &
                   sc*(2.0_dp*v2sigmatau(6,1) - v2sigmatau(4,1))*my_norm_drhob
                e_laplace_rhoa_laplace_rhoa(ii) = e_laplace_rhoa_laplace_rhoa(ii) + sc*v2lapl2(1,1)
                e_laplace_rhoa_laplace_rhob(ii) = e_laplace_rhoa_laplace_rhob(ii) + sc*v2lapl2(2,1)
                e_laplace_rhob_laplace_rhob(ii) = e_laplace_rhob_laplace_rhob(ii) + sc*v2lapl2(3,1)
                e_laplace_rhoa_tau_a(ii) = e_laplace_rhoa_tau_a(ii) + sc*v2lapltau(1,1)
                e_laplace_rhoa_tau_b(ii) = e_laplace_rhoa_tau_b(ii) + sc*v2lapltau(2,1)
                e_laplace_rhob_tau_a(ii) = e_laplace_rhob_tau_a(ii) + sc*v2lapltau(3,1)
                e_laplace_rhob_tau_b(ii) = e_laplace_rhob_tau_b(ii) + sc*v2lapltau(4,1)
                e_tau_a_tau_a(ii) = e_tau_a_tau_a(ii) + sc*v2tau2(1,1)
                e_tau_a_tau_b(ii) = e_tau_a_tau_b(ii) + sc*v2tau2(2,1)
                e_tau_b_tau_b(ii) = e_tau_b_tau_b(ii) + sc*v2tau2(3,1)
             END IF
          END DO
          !$end do
       END IF
    CASE default
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name)//": this XC_FAMILY is currently not supported.")
    END SELECT

    CALL xc_f90_func_end(xc_func)
#else
    CALL stop_program(routineN,moduleN,__LINE__,&
       "In order to use libxc you need to download and install it")
#endif

  END SUBROUTINE libxc_lsd_calc

END MODULE xc_libxc
