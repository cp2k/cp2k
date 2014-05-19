!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates a functional from libxc and its derivatives
!> \note
!>      Only for version v2.0.1 (2013-01-21) of libxc
!>      (Marques, Oliveira, Burnus, CPC 183, 2272 (2012)).
!>
!>      For subsequent versions of libxc, the following should be updated if
!>      necessary:
!>      1) The function get_func_id.
!>      2) The files libxc_funcs_m.F and xc_f90_lib_m.F.
!>      3) The list of functionals using the laplacian of the density
!>         (in libxc_lda_info, libxc_lsd_info, libxc_lda_eval and
!>         libxc_lsd_eval).
!>      4) The list of functionals for which it is possible to provide input
!>         parameters (in libxc_lda_calc and libxc_lsd_calc). For more
!>         information on the parameters, see subroutines xc_f90_xxx_set_par
!>         in libxc.f90 of the libxc package or xc_f90_lib_m.F.
!>      5) The list of functionals for which exc is not calculated by the
!>         libxc subroutines (exc=0 is arbitrarily chosen) (in libxc_lda_calc
!>         and libxc_lsd_calc).
!>      6) Reactivate the functionals which are working correctly
!>         (in libxc_lda_info and libxc_lsd_info).
!>
!>      WARNING: In the subroutine libxc_lsd_calc, it could be that the
!>      ordering for the 1st index of v2lapltau, v2rholapl, v2rhotau,
!>      v2sigmalapl and v2sigmatau is not correct. For the moment it does not
!>      matter since the calculation of the 2nd derivatives for meta-GGA
!>      functionals is not implemented in CP2K.
!>
!> \par History
!>      01.2013 created [F. Tran]
!> \author F. Tran
! *****************************************************************************
MODULE xc_libxc
#if defined (__LIBXC)
This version of CP2K ONLY works with libxc version 2.0.1
#endif
#if defined (__LIBXC2)
  USE bibliography,                    ONLY: Marques2012,&
                                             cite_reference
  USE f77_blas
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
  USE xc_f90_types_m,                  ONLY: xc_f90_pointer_t
  USE xc_f90_lib_m

#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc'

  PUBLIC :: libxc_lda_info, libxc_lda_eval, libxc_lsd_info, libxc_lsd_eval
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_info(libxc_params,reference,shortform,needs,max_deriv,ifunc_name,error)

    TYPE(section_vals_type), POINTER         :: libxc_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    INTEGER, INTENT(in)                      :: ifunc_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=120)                       :: s1, s2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id, ii
    REAL(KIND=dp)                            :: sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: scale
    TYPE(xc_f90_pointer_t)                   :: str, xc_func, xc_info

    CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name,error=error)
    CALL section_vals_val_get(libxc_params,"scale",r_vals=scale,error=error)

    CALL cite_reference(Marques2012)

    IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
      sc = 1.0_dp
    ELSE
      sc = SCALE(ifunc_name)
    ENDIF

    func_id = get_func_id(func_name(ifunc_name))
    IF (func_id == -1) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name(ifunc_name))//": wrong functional name")
    END IF
    IF ((func_id == 207) .OR. (func_id == 208) .OR. (func_id == 209)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name(ifunc_name))//": desactivated since (probably) buggy in v2.0.1 of libxc")
    END IF
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
    END SELECT
    IF ( PRESENT ( shortform ) ) THEN
       shortform = TRIM(s1)//' ('//TRIM(s2)//')'
    END IF
    ii = 0
    CALL xc_f90_info_refs(xc_info, ii, str, s1)
    IF ( PRESENT ( reference ) ) THEN
       WRITE(reference,"(a,' {scale=',f5.3,',spin-unpolarized}')") TRIM(s1), sc
    END IF
    IF (PRESENT(needs)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          needs%rho=.TRUE.
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          needs%rho=.TRUE.
          needs%norm_drho=.TRUE.
       CASE (XC_FAMILY_MGGA)
          needs%rho=.TRUE.
          needs%norm_drho=.TRUE.
          needs%tau=.TRUE.
          IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
              (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
             needs%laplace_rho=.TRUE.
          END IF
       END SELECT
    END IF
    IF (PRESENT(max_deriv)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          max_deriv=3
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          max_deriv=2
       CASE (XC_FAMILY_MGGA)
          max_deriv=1
!          max_deriv=2
       END SELECT
    END IF

    CALL xc_f90_func_end(xc_func)

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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_info(libxc_params,reference,shortform,needs,max_deriv,ifunc_name,error)

    TYPE(section_vals_type), POINTER         :: libxc_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    INTEGER, INTENT(in)                      :: ifunc_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=120)                       :: s1, s2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id, ii
    REAL(KIND=dp)                            :: sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: scale
    TYPE(xc_f90_pointer_t)                   :: str, xc_func, xc_info

    CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name,error=error)
    CALL section_vals_val_get(libxc_params,"scale",r_vals=scale,error=error)

    CALL cite_reference(Marques2012)

    IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
      sc = 1.0_dp
    ELSE
      sc = SCALE(ifunc_name)
    ENDIF

    func_id = get_func_id(func_name(ifunc_name))
    IF (func_id == -1) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name(ifunc_name))//": wrong functional name")
    END IF
    IF ((func_id == 207) .OR. (func_id == 208) .OR. (func_id == 209)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name(ifunc_name))//": desactivated since (probably) buggy in v2.0.1 of libxc")
    END IF
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
    END SELECT
    IF ( PRESENT ( shortform ) ) THEN
       shortform = TRIM(s1)//' ('//TRIM(s2)//')'
    END IF
    ii = 0
    CALL xc_f90_info_refs(xc_info, ii, str, s1)
    IF ( PRESENT ( reference ) ) THEN
       WRITE(reference,"(a,' {scale=',f5.3,',spin-polarized}')") TRIM(s1), sc
    END IF
    IF (PRESENT(needs)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          needs%rho_spin=.TRUE.
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          needs%rho_spin=.TRUE.
          needs%norm_drho=.TRUE.
          needs%norm_drho_spin=.TRUE.
       CASE (XC_FAMILY_MGGA)
          needs%rho_spin=.TRUE.
          needs%norm_drho=.TRUE.
          needs%norm_drho_spin=.TRUE.
          needs%tau_spin=.TRUE.
          IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
              (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
             needs%laplace_rho_spin=.TRUE.
          END IF
       END SELECT
    END IF
    IF (PRESENT(max_deriv)) THEN
       SELECT CASE(xc_f90_info_family(xc_info))
       CASE (XC_FAMILY_LDA)
          max_deriv=3
       CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          max_deriv=2
       CASE (XC_FAMILY_MGGA)
          max_deriv=1
!          max_deriv=2
       END SELECT
    END IF

    CALL xc_f90_func_end(xc_func)

  END SUBROUTINE libxc_lsd_info

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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_eval(rho_set,deriv_set,grad_deriv,libxc_params,ifunc_name,error)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: libxc_params
    INTEGER, INTENT(in)                      :: ifunc_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lda_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id, handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, epsilon_tau, sc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: params, scale
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rho, &
      e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_ndrho, &
      e_ndrho_laplace_rho, e_ndrho_ndrho, e_ndrho_rho, e_ndrho_tau, e_rho, &
      e_rho_laplace_rho, e_rho_rho, e_rho_rho_rho, e_rho_tau, e_tau, &
      e_tau_tau, laplace_rho, norm_drho, rho, tau
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name,error=error)
       CALL section_vals_val_get(libxc_params,"scale",r_vals=scale,error=error)
       CALL section_vals_val_get(libxc_params,"parameters",r_vals=params,error=error)

       IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
         sc = 1.0_dp
       ELSE
         sc = SCALE(ifunc_name)
       ENDIF

       func_id = get_func_id(func_name(ifunc_name))
!$OMP CRITICAL(libxc_init)
       CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_UNPOLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

       SELECT CASE (xc_f90_info_family(xc_info))
       CASE(XC_FAMILY_LDA)
          CALL xc_rho_set_get(rho_set,rho=rho,local_bounds=bo,&
               rho_cutoff=epsilon_rho,error=error)
          norm_drho => rho
          laplace_rho => rho
          tau => rho
       CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          CALL xc_rho_set_get(rho_set,rho=rho,norm_drho=norm_drho,&
               local_bounds=bo,rho_cutoff=epsilon_rho,&
               drho_cutoff=epsilon_norm_drho,error=error)
          laplace_rho => rho
          tau => rho
       CASE(XC_FAMILY_MGGA)
          IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR.&
              (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
             CALL xc_rho_set_get(rho_set,rho=rho,norm_drho=norm_drho,&
                  laplace_rho=laplace_rho,tau=tau,local_bounds=bo,&
                  rho_cutoff=epsilon_rho,drho_cutoff=epsilon_norm_drho,&
                  tau_cutoff=epsilon_tau,error=error)
          ELSE
             CALL xc_rho_set_get(rho_set,rho=rho,norm_drho=norm_drho,&
                  tau=tau,local_bounds=bo,rho_cutoff=epsilon_rho,&
                  drho_cutoff=epsilon_norm_drho,tau_cutoff=epsilon_tau,error=error)
             laplace_rho => rho
          END IF
       END SELECT

       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> rho
       END IF

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
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
          CASE(XC_FAMILY_MGGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_tau,error=error)
             IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
                 (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)",&
                     allocate_deriv=.TRUE.,error=error)
                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho,error=error)
             END IF
          END SELECT
       END IF
       IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
          CASE(XC_FAMILY_MGGA)
!             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rho)(tau)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rho_tau,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(tau)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_tau,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau)(tau)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_tau,error=error)
!             IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
!                 (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
!                deriv => xc_dset_get_derivative(deriv_set,"(rho)(laplace_rho)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_rho_laplace_rho,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_rho)(laplace_rho)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrho_laplace_rho,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)(laplace_rho)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho_laplace_rho,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rho)(tau)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rho_tau,error=error)
!             END IF
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives bigger than 1 not implemented",&
                error=error, error_level=cp_failure_level)
          END SELECT
       END IF
       IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA)
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives bigger than 2 not implemented",&
                error=error, error_level=cp_failure_level)
          END SELECT
       END IF
       IF (grad_deriv>=4.OR.grad_deriv<=-4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP,&
             message="derivatives bigger than 3 not implemented",&
             error=error, error_level=cp_failure_level)
       END IF

       !$omp parallel default(none), &
       !$omp shared(rho,norm_drho,laplace_rho,tau,e_0,e_rho,e_ndrho,e_laplace_rho),&
       !$omp shared(e_tau,e_rho_rho,e_ndrho_rho,e_ndrho_ndrho,e_rho_laplace_rho),&
       !$omp shared(e_rho_tau,e_ndrho_laplace_rho,e_ndrho_tau,e_laplace_rho_laplace_rho),&
       !$omp shared(e_laplace_rho_tau,e_tau_tau,e_rho_rho_rho),&
       !$omp shared(grad_deriv,npoints),&
       !$omp shared(epsilon_rho,epsilon_norm_drho,epsilon_tau),&
       !$omp shared(func_name,ifunc_name,sc,params,error)

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
          epsilon_rho=epsilon_rho,epsilon_norm_drho=epsilon_norm_drho,&
          epsilon_tau=epsilon_tau,func_name=func_name(ifunc_name),&
          sc=sc,params=params,error=error)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF

       CALL xc_f90_func_end(xc_func)

    END IF

    CALL timestop(handle)
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lsd_eval(rho_set,deriv_set,grad_deriv,libxc_params,ifunc_name,error)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: libxc_params
    INTEGER, INTENT(in)                      :: ifunc_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'libxc_lsd_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: func_id, handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, epsilon_tau, sc
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
      e_ndrhoa_tau_b, e_ndrhob, e_ndrhob_laplace_rhoa, e_ndrhob_laplace_rhob, &
      e_ndrhob_ndrhob
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: e_ndrhob_rhoa, &
      e_ndrhob_rhob, e_ndrhob_tau_a, e_ndrhob_tau_b, e_rhoa, &
      e_rhoa_laplace_rhoa, e_rhoa_laplace_rhob, e_rhoa_rhoa, &
      e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, e_rhoa_rhob, e_rhoa_rhob_rhob, &
      e_rhoa_tau_a, e_rhoa_tau_b, e_rhob, e_rhob_laplace_rhoa, &
      e_rhob_laplace_rhob, e_rhob_rhob, e_rhob_rhob_rhob, e_rhob_tau_a, &
      e_rhob_tau_b, e_tau_a, e_tau_a_tau_a, e_tau_a_tau_b, e_tau_b, &
      e_tau_b_tau_b, laplace_rhoa, laplace_rhob, norm_drho, norm_drhoa, &
      norm_drhob, rhoa, rhob, tau_a, tau_b
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_vals_val_get(libxc_params,"functional",c_vals=func_name,error=error)
       CALL section_vals_val_get(libxc_params,"scale",r_vals=scale,error=error)
       CALL section_vals_val_get(libxc_params,"parameters",r_vals=params,error=error)

       IF ((SIZE(scale) == 1) .AND. (ABS(SCALE(1)-1.0_dp) < 1.0e-10_dp)) THEN
         sc = 1.0_dp
       ELSE
         sc = SCALE(ifunc_name)
       ENDIF

       func_id = get_func_id(func_name(ifunc_name))
!$OMP CRITICAL(libxc_init)
       CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_POLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

       SELECT CASE (xc_f90_info_family(xc_info))
       CASE(XC_FAMILY_LDA)
          CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
               local_bounds=bo,rho_cutoff=epsilon_rho,error=error)
          norm_drho => rhoa
          norm_drhoa => rhoa
          norm_drhob => rhoa
          laplace_rhoa => rhoa
          laplace_rhob => rhoa
          tau_a => rhoa
          tau_b => rhoa
       CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
               norm_drho=norm_drho,norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,&
               local_bounds=bo,rho_cutoff=epsilon_rho,&
               drho_cutoff=epsilon_norm_drho,error=error)
          laplace_rhoa => rhoa
          laplace_rhob => rhoa
          tau_a => rhoa
          tau_b => rhoa
       CASE(XC_FAMILY_MGGA)
          IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
              (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
             CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
                  norm_drho=norm_drho,norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,&
                  laplace_rhoa=laplace_rhoa,laplace_rhob=laplace_rhob,&
                  tau_a=tau_a,tau_b=tau_b,local_bounds=bo,&
                  rho_cutoff=epsilon_rho,drho_cutoff=epsilon_norm_drho,&
                  tau_cutoff=epsilon_tau,error=error)
          ELSE
             CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
                  norm_drho=norm_drho,norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,&
                  tau_a=tau_a,tau_b=tau_b,local_bounds=bo,&
                  rho_cutoff=epsilon_rho,drho_cutoff=epsilon_norm_drho,&
                  tau_cutoff=epsilon_tau,error=error)
             laplace_rhoa => rhoa
             laplace_rhob => rhoa
          END IF
       END SELECT

       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> rhoa
       END IF

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
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob,error=error)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob,error=error)
          CASE(XC_FAMILY_MGGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
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
             IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
                 (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)",&
                     allocate_deriv=.TRUE.,error=error)
                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa,error=error)
                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)",&
                     allocate_deriv=.TRUE.,error=error)
                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob,error=error)
             END IF
          END SELECT
       END IF
       IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob,error=error)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(norm_drhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_ndrhob,error=error)
          CASE(XC_FAMILY_MGGA)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhoa,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhoa,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_rhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhoa)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhoa,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_rhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drho)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhoa)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhoa,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(norm_drhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhoa)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhoa,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(norm_drhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_ndrhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(norm_drhob)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_ndrhob,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(tau_a)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_tau_a,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_tau_b,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(tau_a)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhob_tau_a,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_rhob_tau_b,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(tau_a)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_tau_a,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_tau_b,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(tau_a)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_tau_a,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_tau_b,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(tau_a)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_tau_a,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_tau_b,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau_a)(tau_a)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_a_tau_a,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau_a)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_a_tau_b,error=error)
!             deriv => xc_dset_get_derivative(deriv_set,"(tau_b)(tau_b)",&
!                  allocate_deriv=.TRUE.,error=error)
!             CALL xc_derivative_get(deriv,deriv_data=e_tau_b_tau_b,error=error)
!             IF ((func_id == 206) .OR. (func_id == 207) .OR. (func_id == 208) .OR. &
!                 (func_id == 209) .OR. (func_id == 210) .OR. (func_id == 211)) THEN
!                deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhoa_laplace_rhoa,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhoa_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(rhob)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhob_laplace_rhoa,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(rhob)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_rhob_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrho_laplace_rhoa,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrho_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_laplace_rhoa,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_laplace_rhoa,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_ndrhob_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(laplace_rhoa)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_laplace_rhoa,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)(laplace_rhob)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob_laplace_rhob,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(tau_a)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_tau_a,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhoa)(tau_b)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhoa_tau_b,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)(tau_a)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob_tau_a,error=error)
!                deriv => xc_dset_get_derivative(deriv_set,"(laplace_rhob)(tau_b)",&
!                     allocate_deriv=.TRUE.,error=error)
!                CALL xc_derivative_get(deriv,deriv_data=e_laplace_rhob_tau_b,error=error)
!             END IF
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives bigger than 1 not implemented",&
                error=error, error_level=cp_failure_level)
          END SELECT
       END IF
       IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
          SELECT CASE (xc_f90_info_family(xc_info))
          CASE(XC_FAMILY_LDA)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa_rhoa,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhoa_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhoa_rhob_rhob,error=error)
             deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rhob_rhob_rhob,error=error)
          CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA)
             CALL cp_unimplemented_error(fromWhere=routineP,&
                message="derivatives bigger than 2 not implemented",&
                error=error, error_level=cp_failure_level)
          END SELECT
       END IF
       IF (grad_deriv>=4.OR.grad_deriv<=-4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP,&
             message="derivatives bigger than 3 not implemented",&
             error=error, error_level=cp_failure_level)
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
       !$omp shared(epsilon_rho,epsilon_norm_drho,epsilon_tau),&
       !$omp shared(func_name,ifunc_name,sc,params,error)

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
          epsilon_rho=epsilon_rho,epsilon_norm_drho=epsilon_norm_drho,&
          epsilon_tau=epsilon_tau,func_name=func_name(ifunc_name),&
          sc=sc,params=params,error=error)

       !$omp end parallel

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF

       CALL xc_f90_func_end(xc_func)

    END IF

    CALL timestop(handle)
  END SUBROUTINE libxc_lsd_eval

! *****************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rho density
!> \param norm_drho norm of the gradient of the density
!> \param laplace_rho laplacian of the density
!> \param tau kinetic-energy density
!> \param e_0 energy density
!> \param e_* derivative of the energy density with respect to *
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of points on the grid
!> \param func_name name of the functional
!> \param sc scaling factor
!> \param params parameters of the functional
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author F. Tran
! *****************************************************************************
  SUBROUTINE libxc_lda_calc(rho,norm_drho,laplace_rho,tau,&
          e_0,e_rho,e_ndrho,e_laplace_rho,e_tau,e_rho_rho,e_ndrho_rho,&
          e_ndrho_ndrho,e_rho_laplace_rho,e_rho_tau,e_ndrho_laplace_rho,&
          e_ndrho_tau,e_laplace_rho_laplace_rho,e_laplace_rho_tau,&
          e_tau_tau,e_rho_rho_rho,&
          grad_deriv,npoints,epsilon_rho,epsilon_norm_drho,&
          epsilon_tau,func_name,sc,params,error)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, norm_drho, laplace_rho, &
                                                tau
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rho, e_ndrho, &
      e_laplace_rho, e_tau, e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, &
      e_rho_laplace_rho, e_rho_tau, e_ndrho_laplace_rho, e_ndrho_tau, &
      e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_tau_tau, e_rho_rho_rho
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(KIND=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, epsilon_tau
    CHARACTER(LEN=80), INTENT(IN)            :: func_name
    REAL(KIND=dp), INTENT(in)                :: sc
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: func_id, ii
    LOGICAL                                  :: no_exc
    REAL(KIND=dp) :: exc, my_tau, v2lapl2, v2lapltau, v2rho2, v2rholapl, &
      v2rhosigma, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2tau2, &
      v3rho3, vlapl, vrho, vsigma, vtau
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    func_id = get_func_id(func_name)
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_UNPOLARIZED)

    IF (params(1) < 0.99e20_dp) THEN
       IF (func_id == 1) THEN
          CALL xc_f90_lda_x_set_par(xc_func,params(1),NINT(params(2)),params(3))
       ELSE IF (func_id == 6) THEN
          CALL xc_f90_lda_c_xalpha_set_par(xc_func,params(1))
       ELSE IF (func_id == 16) THEN
          CALL xc_f90_lda_c_2d_prm_set_par(xc_func,params(1))
       ELSE IF (func_id == 18) THEN
          CALL xc_f90_lda_c_1d_csc_set_par(xc_func,NINT(params(1)),params(2))
       ELSE IF (func_id == 21) THEN
          CALL xc_f90_lda_x_1d_set_par(xc_func,NINT(params(1)),params(2))
       ELSE IF (func_id == 160) THEN
          CALL xc_f90_gga_lb_set_par(xc_func,NINT(params(1)),params(2),params(3),params(4))
       ELSE IF (func_id == 208) THEN
          CALL xc_f90_mgga_x_tb09_set_par(xc_func,params(1))
       ELSE IF ((func_id == 427) .OR. (func_id == 428)) THEN
          CALL xc_f90_hyb_gga_xc_hse_set_par(xc_func,params(1))
       ELSE IF ((func_id == 429) .OR. (func_id == 430) .OR. (func_id == 431) .OR. &
                (func_id == 432) .OR. (func_id == 525) .OR. (func_id == 526) .OR. &
                (func_id == 527) .OR. (func_id == 528)) THEN
          CALL xc_f90_gga_x_hjs_set_par(xc_func,params(1))
       ELSE IF (func_id == 524) THEN
          CALL xc_f90_gga_x_wpbeh_set_par(xc_func,params(1))
       END IF
    END IF
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

    no_exc = .FALSE.
    IF ((func_id == 160) .OR. (func_id == 182) .OR. (func_id == 207) .OR. &
        (func_id == 208) .OR. (func_id == 209)) no_exc = .TRUE.

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
                IF (no_exc .EQV. .TRUE.) THEN
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
                IF (no_exc .EQV. .TRUE.) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,vrho,vsigma)
                  CALL xc_f90_gga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     v2rho2,v2rhosigma,v2sigma2)
                ELSE
                  CALL xc_f90_gga(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     exc,vrho,vsigma,v2rho2,v2rhosigma,v2sigma2)
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
                IF (no_exc .EQV. .TRUE.) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rho(ii),norm_drho(ii)**2,vrho,vsigma)
                  CALL xc_f90_gga_fxc(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     v2rho2,v2rhosigma,v2sigma2)
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_gga(xc_func,1,rho(ii),norm_drho(ii)**2,&
                     exc,vrho,vsigma,v2rho2,v2rhosigma,v2sigma2)
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
                IF (no_exc .EQV. .TRUE.) THEN
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
                IF (no_exc .EQV. .TRUE.) THEN
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
                IF (no_exc .EQV. .TRUE.) THEN
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
    END SELECT

    CALL xc_f90_func_end(xc_func)

  END SUBROUTINE libxc_lda_calc

! *****************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rhoa, rhob: alpha and beta density
!> \param norm_drhoa, norm_drhob: norm of the gradient of the alpha and beta
!>        density
!> \param laplace_rhoa, laplace_rhob: laplacian of the alpha and beta density
!> \param tau_a, tau_b: alpha and beta kinetic-energy density
!> \param e_0 energy density
!> \param e_* derivative of the energy density with respect to *
!> \param grad_deriv degree of the derivative that should be evaluated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints number of points on the grid
!> \param func_name name of the functional
!> \param sc scaling factor
!> \param params parameters of the functional
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
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
          grad_deriv,npoints,epsilon_rho,epsilon_norm_drho,&
          epsilon_tau,func_name,sc,params,error)

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
      e_ndrho_laplace_rhoa, e_ndrho_laplace_rhob, e_ndrhoa_laplace_rhoa, &
      e_ndrhoa_laplace_rhob, e_ndrhob_laplace_rhoa, e_ndrhob_laplace_rhob
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_ndrho_tau_a, &
      e_ndrho_tau_b, e_ndrhoa_tau_a, e_ndrhoa_tau_b, e_ndrhob_tau_a, &
      e_ndrhob_tau_b, e_laplace_rhoa_laplace_rhoa, &
      e_laplace_rhoa_laplace_rhob, e_laplace_rhob_laplace_rhob, &
      e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, e_laplace_rhob_tau_a, &
      e_laplace_rhob_tau_b, e_tau_a_tau_a, e_tau_a_tau_b, e_tau_b_tau_b, &
      e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, e_rhoa_rhob_rhob, e_rhob_rhob_rhob
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(KIND=dp), INTENT(in)                :: epsilon_rho, &
                                                epsilon_norm_drho, epsilon_tau
    CHARACTER(LEN=80), INTENT(IN)            :: func_name
    REAL(KIND=dp), INTENT(in)                :: sc
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

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

    func_id = get_func_id(func_name)
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func,xc_info,func_id,XC_POLARIZED)

    IF (params(1) < 0.99e20_dp) THEN
       IF (func_id == 1) THEN
          CALL xc_f90_lda_x_set_par(xc_func,params(1),NINT(params(2)),params(3))
       ELSE IF (func_id == 6) THEN
          CALL xc_f90_lda_c_xalpha_set_par(xc_func,params(1))
       ELSE IF (func_id == 16) THEN
          CALL xc_f90_lda_c_2d_prm_set_par(xc_func,params(1))
       ELSE IF (func_id == 18) THEN
          CALL xc_f90_lda_c_1d_csc_set_par(xc_func,NINT(params(1)),params(2))
       ELSE IF (func_id == 21) THEN
          CALL xc_f90_lda_x_1d_set_par(xc_func,NINT(params(1)),params(2))
       ELSE IF (func_id == 160) THEN
          CALL xc_f90_gga_lb_set_par(xc_func,NINT(params(1)),params(2),params(3),params(4))
       ELSE IF (func_id == 208) THEN
          CALL xc_f90_mgga_x_tb09_set_par(xc_func,params(1))
       ELSE IF ((func_id == 427) .OR. (func_id == 428)) THEN
          CALL xc_f90_hyb_gga_xc_hse_set_par(xc_func,params(1))
       ELSE IF ((func_id == 429) .OR. (func_id == 430) .OR. (func_id == 431) .OR. &
                (func_id == 432) .OR. (func_id == 525) .OR. (func_id == 526) .OR. &
                (func_id == 527) .OR. (func_id == 528)) THEN
          CALL xc_f90_gga_x_hjs_set_par(xc_func,params(1))
       ELSE IF (func_id == 524) THEN
          CALL xc_f90_gga_x_wpbeh_set_par(xc_func,params(1))
       END IF
    END IF
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER

    no_exc = .FALSE.
    IF ((func_id == 160) .OR. (func_id == 182) .OR. (func_id == 207) .OR. &
        (func_id == 208) .OR. (func_id == 209)) no_exc = .TRUE.

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
                IF (no_exc .EQV. .TRUE.) THEN
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
                IF (no_exc .EQV. .TRUE.) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),vrho(1,1),vsigma(1,1))
                  CALL xc_f90_gga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
                ELSE
                  CALL xc_f90_gga(xc_func,1,rhov(1,1),sigmav(1,1),&
                     exc,vrho(1,1),vsigma(1,1),v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
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
                IF (no_exc .EQV. .TRUE.) THEN
                  CALL xc_f90_gga_vxc(xc_func,1,rhov(1,1),sigmav(1,1),vrho(1,1),vsigma(1,1))
                  CALL xc_f90_gga_fxc(xc_func,1,rhov(1,1),sigmav(1,1),&
                     v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
                  exc = 0.0_dp
                ELSE
                  CALL xc_f90_gga(xc_func,1,rhov(1,1),sigmav(1,1),&
                     exc,vrho(1,1),vsigma(1,1),v2rho2(1,1),v2rhosigma(1,1),v2sigma2(1,1))
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
                IF (no_exc .EQV. .TRUE.) THEN
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
                IF (no_exc .EQV. .TRUE.) THEN
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
                IF (no_exc .EQV. .TRUE.) THEN
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
    END SELECT

    CALL xc_f90_func_end(xc_func)

  END SUBROUTINE libxc_lsd_calc

! *****************************************************************************
!> \brief returns the functional identifier
!> \param func_name the name of the functional
!> \note
!>      The correspondence between the functional name and identifier must be
!>      the same as in the file libxc_funcs.f90 of v2.0.1 (2013-01-21) of libxc
!>      and libxc_funcs_m.F
!> \author F. Tran
! *****************************************************************************
  FUNCTION get_func_id(func_name)

    CHARACTER(LEN=*), INTENT(IN)             :: func_name
    INTEGER                                  :: get_func_id

    get_func_id = -1
    IF (TRIM(func_name) == 'XC_LDA_X')                      get_func_id =   1  !  Exchange                     
    IF (TRIM(func_name) == 'XC_LDA_C_WIGNER')               get_func_id =   2  !  Wigner parametrization       
    IF (TRIM(func_name) == 'XC_LDA_C_RPA')                  get_func_id =   3  !  Random Phase Approximation   
    IF (TRIM(func_name) == 'XC_LDA_C_HL')                   get_func_id =   4  !  Hedin & Lundqvist            
    IF (TRIM(func_name) == 'XC_LDA_C_GL')                   get_func_id =   5  !  Gunnarson & Lundqvist        
    IF (TRIM(func_name) == 'XC_LDA_C_XALPHA')               get_func_id =   6  !  Slater Xalpha                
    IF (TRIM(func_name) == 'XC_LDA_C_VWN')                  get_func_id =   7  !  Vosko, Wilk, & Nussair (5)   
    IF (TRIM(func_name) == 'XC_LDA_C_VWN_RPA')              get_func_id =   8  !  Vosko, Wilk, & Nussair (RPA) 
    IF (TRIM(func_name) == 'XC_LDA_C_PZ')                   get_func_id =   9  !  Perdew & Zunger              
    IF (TRIM(func_name) == 'XC_LDA_C_PZ_MOD')               get_func_id =  10  !  Perdew & Zunger (Modified)   
    IF (TRIM(func_name) == 'XC_LDA_C_OB_PZ')                get_func_id =  11  !  Ortiz & Ballone (PZ)         
    IF (TRIM(func_name) == 'XC_LDA_C_PW')                   get_func_id =  12  !  Perdew & Wang                
    IF (TRIM(func_name) == 'XC_LDA_C_PW_MOD')               get_func_id =  13  !  Perdew & Wang (Modified)     
    IF (TRIM(func_name) == 'XC_LDA_C_OB_PW')                get_func_id =  14  !  Ortiz & Ballone (PW)         
    IF (TRIM(func_name) == 'XC_LDA_C_2D_AMGB')              get_func_id =  15  !  Attacalite et al             
    IF (TRIM(func_name) == 'XC_LDA_C_2D_PRM')               get_func_id =  16  !  Pittalis, Rasanen & Marques correlation in 2D 
    IF (TRIM(func_name) == 'XC_LDA_C_vBH')                  get_func_id =  17  !  von Barth & Hedin            
    IF (TRIM(func_name) == 'XC_LDA_C_1D_CSC')               get_func_id =  18  !  Casula, Sorella, and Senatore 1D correlation     
    IF (TRIM(func_name) == 'XC_LDA_X_2D')                   get_func_id =  19  !  Exchange in 2D 
    IF (TRIM(func_name) == 'XC_LDA_XC_TETER93')             get_func_id =  20  !  Teter 93 parametrization                
    IF (TRIM(func_name) == 'XC_LDA_X_1D')                   get_func_id =  21  !  Exchange in 1D     
    IF (TRIM(func_name) == 'XC_LDA_C_ML1')                  get_func_id =  22  !  Modified LSD (version 1) of Proynov and Salahub 
    IF (TRIM(func_name) == 'XC_LDA_C_ML2')                  get_func_id =  23  !  Modified LSD (version 2) of Proynov and Salahub 
    IF (TRIM(func_name) == 'XC_LDA_C_GOMBAS')               get_func_id =  24  !  Gombas parametrization       
    IF (TRIM(func_name) == 'XC_LDA_C_PW_RPA')               get_func_id =  25  !  Perdew & Wang fit of the RPA 
    IF (TRIM(func_name) == 'XC_LDA_C_1D_LOOS')              get_func_id =  26  !  P-F Loos correlation LDA     
    IF (TRIM(func_name) == 'XC_LDA_C_RC04')                 get_func_id =  27  !  Ragot-Cortona 
    IF (TRIM(func_name) == 'XC_LDA_C_VWN_1')                get_func_id =  28  !  Vosko, Wilk, & Nussair (1)   
    IF (TRIM(func_name) == 'XC_LDA_C_VWN_2')                get_func_id =  29  !  Vosko, Wilk, & Nussair (2)   
    IF (TRIM(func_name) == 'XC_LDA_C_VWN_3')                get_func_id =  30  !  Vosko, Wilk, & Nussair (3)   
    IF (TRIM(func_name) == 'XC_LDA_C_VWN_4')                get_func_id =  31  !  Vosko, Wilk, & Nussair (4)   
    IF (TRIM(func_name) == 'XC_LDA_K_TF')                   get_func_id =  50  !  Thomas-Fermi kinetic energy functional 
    IF (TRIM(func_name) == 'XC_LDA_K_LP')                   get_func_id =  51  !  Lee and Parr Gaussian ansatz           
    IF (TRIM(func_name) == 'XC_GGA_C_OP_XALPHA')            get_func_id =  84  !  one-parameter progressive functional (G96 version)  
    IF (TRIM(func_name) == 'XC_GGA_C_OP_G96')               get_func_id =  85  !  one-parameter progressive functional (G96 version)  
    IF (TRIM(func_name) == 'XC_GGA_C_OP_PBE')               get_func_id =  86  !  one-parameter progressive functional (PBE version)  
    IF (TRIM(func_name) == 'XC_GGA_C_OP_B88')               get_func_id =  87  !  one-parameter progressive functional (B88 version) 
    IF (TRIM(func_name) == 'XC_GGA_C_FT97')                 get_func_id =  88  !  Filatov & Thiel correlation 
    IF (TRIM(func_name) == 'XC_GGA_C_SPBE')                 get_func_id =  89  !  PBE correlation to be used with the SSB exchange   
    IF (TRIM(func_name) == 'XC_GGA_X_SSB_SW')               get_func_id =  90  !  Swarta, Sola and Bickelhaupt correction to PBE  
    IF (TRIM(func_name) == 'XC_GGA_X_SSB')                  get_func_id =  91  !  Swarta, Sola and Bickelhaupt  
    IF (TRIM(func_name) == 'XC_GGA_X_SSB_D')                get_func_id =  92  !  Swarta, Sola and Bickelhaupt dispersion  
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_407P')           get_func_id =  93  !  HCTH/407+                                
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_P76')            get_func_id =  94  !  HCTH p=7/6                               
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_P14')            get_func_id =  95  !  HCTH p=1/4                               
    IF (TRIM(func_name) == 'XC_GGA_XC_B97_GGA1')            get_func_id =  96  !  Becke 97 GGA-1                           
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_A')              get_func_id =  97  !  HCTH-A                                   
    IF (TRIM(func_name) == 'XC_GGA_X_BPCCAC')               get_func_id =  98  !  BPCCAC (GRAC for the energy) 
    IF (TRIM(func_name) == 'XC_GGA_C_REVTCA')               get_func_id =  99  !  Tognetti, Cortona, Adamo (revised) 
    IF (TRIM(func_name) == 'XC_GGA_C_TCA')                  get_func_id = 100  !  Tognetti, Cortona, Adamo 
    IF (TRIM(func_name) == 'XC_GGA_X_PBE')                  get_func_id = 101  !  Perdew, Burke & Ernzerhof exchange             
    IF (TRIM(func_name) == 'XC_GGA_X_PBE_R')                get_func_id = 102  !  Perdew, Burke & Ernzerhof exchange (revised)   
    IF (TRIM(func_name) == 'XC_GGA_X_B86')                  get_func_id = 103  !  Becke 86 Xalfa,beta,gamma                      
    IF (TRIM(func_name) == 'XC_GGA_X_HERMAN')               get_func_id = 104  !  Herman et al original GGA                  
    IF (TRIM(func_name) == 'XC_GGA_X_B86_MGC')              get_func_id = 105  !  Becke 86 Xalfa,beta,gamma (with mod. grad. correction) 
    IF (TRIM(func_name) == 'XC_GGA_X_B88')                  get_func_id = 106  !  Becke 88 
    IF (TRIM(func_name) == 'XC_GGA_X_G96')                  get_func_id = 107  !  Gill 96                                        
    IF (TRIM(func_name) == 'XC_GGA_X_PW86')                 get_func_id = 108  !  Perdew & Wang 86 
    IF (TRIM(func_name) == 'XC_GGA_X_PW91')                 get_func_id = 109  !  Perdew & Wang 91 
    IF (TRIM(func_name) == 'XC_GGA_X_OPTX')                 get_func_id = 110  !  Handy & Cohen OPTX 01                          
    IF (TRIM(func_name) == 'XC_GGA_X_DK87_R1')              get_func_id = 111  !  dePristo & Kress 87 (version R1)               
    IF (TRIM(func_name) == 'XC_GGA_X_DK87_R2')              get_func_id = 112  !  dePristo & Kress 87 (version R2)               
    IF (TRIM(func_name) == 'XC_GGA_X_LG93')                 get_func_id = 113  !  Lacks & Gordon 93 
    IF (TRIM(func_name) == 'XC_GGA_X_FT97_A')               get_func_id = 114  !  Filatov & Thiel 97 (version A) 
    IF (TRIM(func_name) == 'XC_GGA_X_FT97_B')               get_func_id = 115  !  Filatov & Thiel 97 (version B) 
    IF (TRIM(func_name) == 'XC_GGA_X_PBE_SOL')              get_func_id = 116  !  Perdew, Burke & Ernzerhof exchange (solids)    
    IF (TRIM(func_name) == 'XC_GGA_X_RPBE')                 get_func_id = 117  !  Hammer, Hansen & Norskov (PBE-like) 
    IF (TRIM(func_name) == 'XC_GGA_X_WC')                   get_func_id = 118  !  Wu & Cohen 
    IF (TRIM(func_name) == 'XC_GGA_X_MPW91')                get_func_id = 119  !  Modified form of PW91 by Adamo & Barone 
    IF (TRIM(func_name) == 'XC_GGA_X_AM05')                 get_func_id = 120  !  Armiento & Mattsson 05 exchange                
    IF (TRIM(func_name) == 'XC_GGA_X_PBEA')                 get_func_id = 121  !  Madsen (PBE-like) 
    IF (TRIM(func_name) == 'XC_GGA_X_MPBE')                 get_func_id = 122  !  Adamo & Barone modification to PBE             
    IF (TRIM(func_name) == 'XC_GGA_X_XPBE')                 get_func_id = 123  !  xPBE reparametrization by Xu & Goddard         
    IF (TRIM(func_name) == 'XC_GGA_X_2D_B86_MGC')           get_func_id = 124  !  Becke 86 MGC for 2D systems 
    IF (TRIM(func_name) == 'XC_GGA_X_BAYESIAN')             get_func_id = 125  !  Bayesian best fit for the enhancement factor 
    IF (TRIM(func_name) == 'XC_GGA_X_PBE_JSJR')             get_func_id = 126  !  JSJR reparametrization by Pedroza, Silva & Capelle 
    IF (TRIM(func_name) == 'XC_GGA_X_2D_B88')               get_func_id = 127  !  Becke 88 in 2D 
    IF (TRIM(func_name) == 'XC_GGA_X_2D_B86')               get_func_id = 128  !  Becke 86 Xalfa,beta,gamma                      
    IF (TRIM(func_name) == 'XC_GGA_X_2D_PBE')               get_func_id = 129  !  Perdew, Burke & Ernzerhof exchange in 2D          
    IF (TRIM(func_name) == 'XC_GGA_C_PBE')                  get_func_id = 130  !  Perdew, Burke & Ernzerhof correlation              
    IF (TRIM(func_name) == 'XC_GGA_C_LYP')                  get_func_id = 131  !  Lee, Yang & Parr 
    IF (TRIM(func_name) == 'XC_GGA_C_P86')                  get_func_id = 132  !  Perdew 86 
    IF (TRIM(func_name) == 'XC_GGA_C_PBE_SOL')              get_func_id = 133  !  Perdew, Burke & Ernzerhof correlation SOL          
    IF (TRIM(func_name) == 'XC_GGA_C_PW91')                 get_func_id = 134  !  Perdew & Wang 91 
    IF (TRIM(func_name) == 'XC_GGA_C_AM05')                 get_func_id = 135  !  Armiento & Mattsson 05 correlation             
    IF (TRIM(func_name) == 'XC_GGA_C_XPBE')                 get_func_id = 136  !  xPBE reparametrization by Xu & Goddard             
    IF (TRIM(func_name) == 'XC_GGA_C_LM')                   get_func_id = 137  !  Langreth and Mehl correlation          
    IF (TRIM(func_name) == 'XC_GGA_C_PBE_JRGX')             get_func_id = 138  !  JRGX reparametrization by Pedroza, Silva & Capelle 
    IF (TRIM(func_name) == 'XC_GGA_X_OPTB88_VDW')           get_func_id = 139  !  Becke 88 reoptimized to be used with vdW functional of Dion et al
    IF (TRIM(func_name) == 'XC_GGA_X_PBEK1_VDW')            get_func_id = 140  !  PBE reparametrization for vdW 
    IF (TRIM(func_name) == 'XC_GGA_X_OPTPBE_VDW')           get_func_id = 141  !  PBE reparametrization for vdW 
    IF (TRIM(func_name) == 'XC_GGA_X_RGE2')                 get_func_id = 142  !  Regularized PBE 
    IF (TRIM(func_name) == 'XC_GGA_C_RGE2')                 get_func_id = 143  !  Regularized PBE                                    
    IF (TRIM(func_name) == 'XC_GGA_X_RPW86')                get_func_id = 144  !  refitted Perdew & Wang 86 
    IF (TRIM(func_name) == 'XC_GGA_X_KT1')                  get_func_id = 145  !  Keal and Tozer version 1             
    IF (TRIM(func_name) == 'XC_GGA_XC_KT2')                 get_func_id = 146  !  Keal and Tozer version 2             
    IF (TRIM(func_name) == 'XC_GGA_C_WL')                   get_func_id = 147  !  Wilson & Levy 
    IF (TRIM(func_name) == 'XC_GGA_C_WI')                   get_func_id = 148  !  Wilson & Ivanov 
    IF (TRIM(func_name) == 'XC_GGA_X_MB88')                 get_func_id = 149  !  Modified Becke 88 for proton transfer 
    IF (TRIM(func_name) == 'XC_GGA_X_SOGGA')                get_func_id = 150  !  Second-order generalized gradient approximation 
    IF (TRIM(func_name) == 'XC_GGA_X_SOGGA11')              get_func_id = 151  !  Second-order generalized gradient approximation 2011 
    IF (TRIM(func_name) == 'XC_GGA_C_SOGGA11')              get_func_id = 152  !  Second-order generalized gradient approximation 2011 
    IF (TRIM(func_name) == 'XC_GGA_C_WI0')                  get_func_id = 153  !  Wilson & Ivanov initial version 
    IF (TRIM(func_name) == 'XC_GGA_XC_TH1')                 get_func_id = 154  !  Tozer and Handy v. 1 
    IF (TRIM(func_name) == 'XC_GGA_XC_TH2')                 get_func_id = 155  !  Tozer and Handy v. 2 
    IF (TRIM(func_name) == 'XC_GGA_XC_TH3')                 get_func_id = 156  !  Tozer and Handy v. 3 
    IF (TRIM(func_name) == 'XC_GGA_XC_TH4')                 get_func_id = 157  !  Tozer and Handy v. 4 
    IF (TRIM(func_name) == 'XC_GGA_X_C09X')                 get_func_id = 158  !  C09x to be used with the VdW of Rutgers-Chalmers     
    IF (TRIM(func_name) == 'XC_GGA_C_SOGGA11_X')            get_func_id = 159  !  To be used with hyb_gga_x_SOGGA11-X  
    IF (TRIM(func_name) == 'XC_GGA_X_LB')                   get_func_id = 160  !  van Leeuwen & Baerends 
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_93')             get_func_id = 161  !  HCTH functional fitted to  93 molecules  
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_120')            get_func_id = 162  !  HCTH functional fitted to 120 molecules  
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_147')            get_func_id = 163  !  HCTH functional fitted to 147 molecules  
    IF (TRIM(func_name) == 'XC_GGA_XC_HCTH_407')            get_func_id = 164  !  HCTH functional fitted to 407 molecules  
    IF (TRIM(func_name) == 'XC_GGA_XC_EDF1')                get_func_id = 165  !  Empirical functionals from Adamson, Gill, and Pople 
    IF (TRIM(func_name) == 'XC_GGA_XC_XLYP')                get_func_id = 166  !  XLYP functional 
    IF (TRIM(func_name) == 'XC_GGA_XC_B97')                 get_func_id = 167  !  Becke 97                                 
    IF (TRIM(func_name) == 'XC_GGA_XC_B97_1')               get_func_id = 168  !  Becke 97-1                               
    IF (TRIM(func_name) == 'XC_GGA_XC_B97_2')               get_func_id = 169  !  Becke 97-2                               
    IF (TRIM(func_name) == 'XC_GGA_XC_B97_D')               get_func_id = 170  !  Grimme functional to be used with C6 vdW term 
    IF (TRIM(func_name) == 'XC_GGA_XC_B97_K')               get_func_id = 171  !  Boese-Martin for Kinetics                
    IF (TRIM(func_name) == 'XC_GGA_XC_B97_3')               get_func_id = 172  !  Becke 97-3                               
    IF (TRIM(func_name) == 'XC_GGA_XC_PBE1W')               get_func_id = 173  !  Functionals fitted for water 
    IF (TRIM(func_name) == 'XC_GGA_XC_MPWLYP1W')            get_func_id = 174  !  Functionals fitted for water 
    IF (TRIM(func_name) == 'XC_GGA_XC_PBELYP1W')            get_func_id = 175  !  Functionals fitted for water 
    IF (TRIM(func_name) == 'XC_GGA_XC_SB98_1a')             get_func_id = 176  !  Schmider-Becke 98 parameterization 1a    
    IF (TRIM(func_name) == 'XC_GGA_XC_SB98_1b')             get_func_id = 177  !  Schmider-Becke 98 parameterization 1b    
    IF (TRIM(func_name) == 'XC_GGA_XC_SB98_1c')             get_func_id = 178  !  Schmider-Becke 98 parameterization 1c    
    IF (TRIM(func_name) == 'XC_GGA_XC_SB98_2a')             get_func_id = 179  !  Schmider-Becke 98 parameterization 2a    
    IF (TRIM(func_name) == 'XC_GGA_XC_SB98_2b')             get_func_id = 180  !  Schmider-Becke 98 parameterization 2b    
    IF (TRIM(func_name) == 'XC_GGA_XC_SB98_2c')             get_func_id = 181  !  Schmider-Becke 98 parameterization 2c    
    IF (TRIM(func_name) == 'XC_GGA_X_LBM')                  get_func_id = 182  !  van Leeuwen & Baerends modified
    IF (TRIM(func_name) == 'XC_GGA_X_OL2')                  get_func_id = 183  !  Exchange form based on Ou-Yang and Levy v.2 
    IF (TRIM(func_name) == 'XC_GGA_X_APBE')                 get_func_id = 184  !  mu fixed from the semiclassical neutral atom   
    IF (TRIM(func_name) == 'XC_GGA_K_APBE')                 get_func_id = 185  !  mu fixed from the semiclassical neutral atom   
    IF (TRIM(func_name) == 'XC_GGA_C_APBE')                 get_func_id = 186  !  mu fixed from the semiclassical neutral atom       
    IF (TRIM(func_name) == 'XC_GGA_K_TW1')                  get_func_id = 187  !  Tran and Wesolowski set 1 (Table II)           
    IF (TRIM(func_name) == 'XC_GGA_K_TW2')                  get_func_id = 188  !  Tran and Wesolowski set 2 (Table II)           
    IF (TRIM(func_name) == 'XC_GGA_K_TW3')                  get_func_id = 189  !  Tran and Wesolowski set 3 (Table II)           
    IF (TRIM(func_name) == 'XC_GGA_K_TW4')                  get_func_id = 190  !  Tran and Wesolowski set 4 (Table II)           
    IF (TRIM(func_name) == 'XC_GGA_X_HTBS')                 get_func_id = 191  !  Haas, Tran, Blaha, and Schwarz  
    IF (TRIM(func_name) == 'XC_GGA_X_AIRY')                 get_func_id = 192  !  Constantin et al based on the Airy gas 
    IF (TRIM(func_name) == 'XC_GGA_X_LAG')                  get_func_id = 193  !  Local Airy Gas 
    IF (TRIM(func_name) == 'XC_GGA_XC_MOHLYP')              get_func_id = 194  !  Functional for organometallic chemistry 
    IF (TRIM(func_name) == 'XC_GGA_XC_MOHLYP2')             get_func_id = 195  !  Functional for barrier heights 
    IF (TRIM(func_name) == 'XC_GGA_XC_TH_FL')               get_func_id = 196  !  Tozer and Handy v. FL  
    IF (TRIM(func_name) == 'XC_GGA_XC_TH_FC')               get_func_id = 197  !  Tozer and Handy v. FC  
    IF (TRIM(func_name) == 'XC_GGA_XC_TH_FCFO')             get_func_id = 198  !  Tozer and Handy v. FCFO 
    IF (TRIM(func_name) == 'XC_GGA_XC_TH_FCO')              get_func_id = 199  !  Tozer and Handy v. FCO 
    IF (TRIM(func_name) == 'XC_GGA_C_OPTC')                 get_func_id = 200  !  Optimized correlation functional of Cohen and Handy 
    IF (TRIM(func_name) == 'XC_GGA_K_VW')                   get_func_id = 500  !  von Weiszaecker functional 
    IF (TRIM(func_name) == 'XC_GGA_K_GE2')                  get_func_id = 501  !  Second-order gradient expansion (l = 1/9) 
    IF (TRIM(func_name) == 'XC_GGA_K_GOLDEN')               get_func_id = 502  !  TF-lambda-vW form by Golden (l = 13/45) 
    IF (TRIM(func_name) == 'XC_GGA_K_YT65')                 get_func_id = 503  !  TF-lambda-vW form by Yonei and Tomishima (l = 1/5) 
    IF (TRIM(func_name) == 'XC_GGA_K_BALTIN')               get_func_id = 504  !  TF-lambda-vW form by Baltin (l = 5/9) 
    IF (TRIM(func_name) == 'XC_GGA_K_LIEB')                 get_func_id = 505  !  TF-lambda-vW form by Lieb (l = 0.185909191) 
    IF (TRIM(func_name) == 'XC_GGA_K_ABSR1')                get_func_id = 506  !  gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)] 
    IF (TRIM(func_name) == 'XC_GGA_K_ABSR2')                get_func_id = 507  !  gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)] 
    IF (TRIM(func_name) == 'XC_GGA_K_GR')                   get_func_id = 508  !  gamma-TFvW form by Gázquez and Robles 
    IF (TRIM(func_name) == 'XC_GGA_K_LUDENA')               get_func_id = 509  !  gamma-TFvW form by Ludeña 
    IF (TRIM(func_name) == 'XC_GGA_K_GP85')                 get_func_id = 510  !  gamma-TFvW form by Ghosh and Parr 
    IF (TRIM(func_name) == 'XC_GGA_K_PEARSON')              get_func_id = 511  !  Pearson 
    IF (TRIM(func_name) == 'XC_GGA_K_OL1')                  get_func_id = 512  !  Ou-Yang and Levy v.1 
    IF (TRIM(func_name) == 'XC_GGA_K_OL2')                  get_func_id = 513  !  Ou-Yang and Levy v.2 
    IF (TRIM(func_name) == 'XC_GGA_K_FR_B88')               get_func_id = 514  !  Fuentealba & Reyes (B88 version) 
    IF (TRIM(func_name) == 'XC_GGA_K_FR_PW86')              get_func_id = 515  !  Fuentealba & Reyes (PW86 version) 
    IF (TRIM(func_name) == 'XC_GGA_K_DK')                   get_func_id = 516  !  DePristo and Kress                    
    IF (TRIM(func_name) == 'XC_GGA_K_PERDEW')               get_func_id = 517  !  Perdew                                
    IF (TRIM(func_name) == 'XC_GGA_K_VSK')                  get_func_id = 518  !  Vitos, Skriver, and Kollar            
    IF (TRIM(func_name) == 'XC_GGA_K_VJKS')                 get_func_id = 519  !  Vitos, Johansson, Kollar, and Skriver 
    IF (TRIM(func_name) == 'XC_GGA_K_ERNZERHOF')            get_func_id = 520  !  Ernzerhof 
    IF (TRIM(func_name) == 'XC_GGA_K_LC94')                 get_func_id = 521  !  Lembarki & Chermette 
    IF (TRIM(func_name) == 'XC_GGA_K_LLP')                  get_func_id = 522  !  Lee, Lee & Parr 
    IF (TRIM(func_name) == 'XC_GGA_K_THAKKAR')              get_func_id = 523  !  Thakkar 1992 
    IF (TRIM(func_name) == 'XC_GGA_X_WPBEH')                get_func_id = 524  !  short-range version of the PBE 
    IF (TRIM(func_name) == 'XC_GGA_X_HJS_PBE')              get_func_id = 525  !  HJS screened exchange PBE version 
    IF (TRIM(func_name) == 'XC_GGA_X_HJS_PBE_SOL')          get_func_id = 526  !  HJS screened exchange PBE_SOL version 
    IF (TRIM(func_name) == 'XC_GGA_X_HJS_B88')              get_func_id = 527  !  HJS screened exchange B88 version 
    IF (TRIM(func_name) == 'XC_GGA_X_HJS_B97X')             get_func_id = 528  !  HJS screened exchange B97x version 
    IF (TRIM(func_name) == 'XC_GGA_X_ITYH')                 get_func_id = 529  !  short-range recipe for exchange GGA functionals 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B3PW91')          get_func_id = 401  !  The original (ACM) hybrid of Becke    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B3LYP')           get_func_id = 402  !  The (in)famous B3LYP                  
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B3P86')           get_func_id = 403  !  Perdew 86 hybrid similar to B3PW91    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_O3LYP')           get_func_id = 404  !  hybrid using the optx functional 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_mPW1K')           get_func_id = 405  !  mixture of mPW91 and PW91 optimized for kinetics 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_PBEH')            get_func_id = 406  !  aka PBE0 or PBE1PBE 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B97')             get_func_id = 407  !  Becke 97                                 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B97_1')           get_func_id = 408  !  Becke 97-1                               
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B97_2')           get_func_id = 410  !  Becke 97-2                               
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_X3LYP')           get_func_id = 411  !  maybe the best hybrid 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B1WC')            get_func_id = 412  !  Becke 1-parameter mixture of WC and PBE 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B97_K')           get_func_id = 413  !  Boese-Martin for Kinetics                
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B97_3')           get_func_id = 414  !  Becke 97-3                               
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_MPW3PW')          get_func_id = 415  !  mixture with the mPW functional       
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B1LYP')           get_func_id = 416  !  Becke 1-parameter mixture of B88 and LYP 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_B1PW91')          get_func_id = 417  !  Becke 1-parameter mixture of B88 and PW91 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_mPW1PW')          get_func_id = 418  !  Becke 1-parameter mixture of mPW91 and PW91 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_MPW3LYP')         get_func_id = 419  !  mixture of mPW and LYP                
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_SB98_1a')         get_func_id = 420  !  Schmider-Becke 98 parameterization 1a    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_SB98_1b')         get_func_id = 421  !  Schmider-Becke 98 parameterization 1b    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_SB98_1c')         get_func_id = 422  !  Schmider-Becke 98 parameterization 1c    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_SB98_2a')         get_func_id = 423  !  Schmider-Becke 98 parameterization 2a    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_SB98_2b')         get_func_id = 424  !  Schmider-Becke 98 parameterization 2b    
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_SB98_2c')         get_func_id = 425  !  Schmider-Becke 98 parameterization 2c    
    IF (TRIM(func_name) == 'XC_HYB_GGA_X_SOGGA11_X')        get_func_id = 426  !  Hybrid based on SOGGA11 form 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_HSE03')           get_func_id = 427  !  the 2003 version of the screened hybrid HSE 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_HSE06')           get_func_id = 428  !  the 2006 version of the screened hybrid HSE 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_HJS_PBE')         get_func_id = 429  !  HJS hybrid screened exchange PBE version 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_HJS_PBE_SOL')     get_func_id = 430  !  HJS hybrid screened exchange PBE_SOL version 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_HJS_B88')         get_func_id = 431  !  HJS hybrid screened exchange B88 version 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_HJS_B97X')        get_func_id = 432  !  HJS hybrid screened exchange B97x version 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_CAM_B3LYP')       get_func_id = 433  !  CAM version of B3LYP 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_TUNED_CAM_B3LYP') get_func_id = 434  !  CAM version of B3LYP tunes for excitations
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_BHANDH')          get_func_id = 435  !  Becke half-and-half 
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_BHANDHLYP')       get_func_id = 436  !  Becke half-and-half with B88 exchange
    IF (TRIM(func_name) == 'XC_HYB_GGA_XC_MB3LYP_RC04')     get_func_id = 437  !  B3LYP with RC04 LDA                   
    IF (TRIM(func_name) == 'XC_MGGA_X_LTA')                 get_func_id = 201  !  Local tau approximation of Ernzerhof & Scuseria 
    IF (TRIM(func_name) == 'XC_MGGA_X_TPSS')                get_func_id = 202  !  Perdew, Tao, Staroverov & Scuseria exchange 
    IF (TRIM(func_name) == 'XC_MGGA_X_M06_L')               get_func_id = 203  !  M06-Local functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_GVT4')                get_func_id = 204  !  GVT4 from Van Voorhis and Scuseria 
    IF (TRIM(func_name) == 'XC_MGGA_X_TAU_HCTH')            get_func_id = 205  !  tau-HCTH from Boese and Handy 
    IF (TRIM(func_name) == 'XC_MGGA_X_BR89')                get_func_id = 206  !  Becke-Roussel 89  
    IF (TRIM(func_name) == 'XC_MGGA_X_BJ06')                get_func_id = 207  !  Becke & Johnson correction to Becke-Roussel 89  
    IF (TRIM(func_name) == 'XC_MGGA_X_TB09')                get_func_id = 208  !  Tran & Blaha correction to Becke & Johnson  
    IF (TRIM(func_name) == 'XC_MGGA_X_RPP09')               get_func_id = 209  !  Rasanen, Pittalis, and Proetto correction to Becke & Johnson  
    IF (TRIM(func_name) == 'XC_MGGA_X_2D_PRHG07')           get_func_id = 210  !  Pittalis, Rasanen, Helbig, Gross Exchange Functional 
    IF (TRIM(func_name) == 'XC_MGGA_X_2D_PRHG07_PRP10')     get_func_id = 211  !  PRGH07 with PRP10 correction 
    IF (TRIM(func_name) == 'XC_MGGA_X_REVTPSS')             get_func_id = 212  !  revised Perdew, Tao, Staroverov & Scuseria exchange 
    IF (TRIM(func_name) == 'XC_MGGA_X_PKZB')                get_func_id = 213  !  Perdew, Kurth, Zupan, and Blaha 
    IF (TRIM(func_name) == 'XC_MGGA_X_M05')                 get_func_id = 214  !  M05 functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_M05_2X')              get_func_id = 215  !  M05-2X functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_M06_HF')              get_func_id = 216  !  M06-HF functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_M06')                 get_func_id = 217  !  M06 functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_M06_2X')              get_func_id = 218  !  M06-2X functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_M08_HX')              get_func_id = 219  !  M08-HX functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_X_M08_SO')              get_func_id = 220  !  M08-SO functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_TPSS')                get_func_id = 231  !  Perdew, Tao, Staroverov & Scuseria correlation 
    IF (TRIM(func_name) == 'XC_MGGA_C_VSXC')                get_func_id = 232  !  VSxc from Van Voorhis and Scuseria (correlation part) 
    IF (TRIM(func_name) == 'XC_MGGA_C_M06_L')               get_func_id = 233  !  M06-Local functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_M06_HF')              get_func_id = 234  !  M06-HF functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_M06')                 get_func_id = 235  !  M06 functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_M06_2X')              get_func_id = 236  !  M06-2X functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_M05')                 get_func_id = 237  !  M05 functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_M05_2X')              get_func_id = 238  !  M05-2X functional of Minnesota 
    IF (TRIM(func_name) == 'XC_MGGA_C_PKZB')                get_func_id = 239  !  Perdew, Kurth, Zupan, and Blaha 
    IF (TRIM(func_name) == 'XC_MGGA_C_BC95')                get_func_id = 240  !  Becke correlation 95 
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_M05')            get_func_id = 438  !  M05 functional of Minnesota                      
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_M05_2X')         get_func_id = 439  !  M05-2X functional of Minnesota                   
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_B88B95')         get_func_id = 440  !  Mixture of B88 with BC95 (B1B95)                 
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_B86B95')         get_func_id = 441  !  Mixture of B86 with BC95                         
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_PW86B95')        get_func_id = 442  !  Mixture of PW86 with BC95                        
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_BB1K')           get_func_id = 443  !  Mixture of B88 with BC95 from Zhao and Truhlar   
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_M06_HF')         get_func_id = 444  !  M06-HF functional of Minnesota                   
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_MPW1B95')        get_func_id = 445  !  Mixture of mPW91 with BC95 from Zhao and Truhlar 
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_MPWB1K')         get_func_id = 446  !  Mixture of mPW91 with BC95 for kinetics          
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_X1B95')          get_func_id = 447  !  Mixture of X with BC95                           
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_XB1K')           get_func_id = 448  !  Mixture of X with BC95 for kinetics              
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_M06')            get_func_id = 449  !  M06 functional of Minnesota                      
    IF (TRIM(func_name) == 'XC_HYB_MGGA_XC_M06_2X')         get_func_id = 450  !  M06-2X functional of Minnesota                   

  END FUNCTION get_func_id
#endif
END MODULE xc_libxc
