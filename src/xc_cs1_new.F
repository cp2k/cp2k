!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_cs1 [1.0] *
!!
!!   NAME
!!     xc_cs1
!!
!!   FUNCTION
!!     Calculate the CS1 Functional (Handy s improved LYP functional)
!!
!!   AUTHOR
!!     JGH (16.03.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!******************************************************************************

MODULE xc_cs1_new

! *****************************************************************************

  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
       xc_dset_get_derivative
  USE xc_functionals_utilities,        ONLY: calc_rho13,&
                                             set_util,&
                                             setup_calculation
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
       cp_assert,&
       cp_error_get_logger,&
       cp_error_message,&
       cp_error_type,&
       cp_unimplemented_error,&
       cp_assertion_failed,&
       cp_debug
  USE cp_log_handling,                 ONLY: cp_failure_level,&
       cp_fatal_level,&
       cp_note_level,&
       cp_to_string,&
       cp_warning_level
  USE timings,                         ONLY: timeset,&
       timestop
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
       xc_derivative_type
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
       xc_rho_set_get,&
       xc_rho_set_type


  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: cs1_lda_info, cs1_lda_eval, cs1_lsd_info, cs1_lsd_eval

  REAL(dbl) :: eps_rho
  LOGICAL :: debug_flag
  REAL(dbl) :: fsig

  REAL(dbl), PARAMETER :: a = 0.04918_dbl, &
                          b = 0.132_dbl, &
                          c = 0.2533_dbl, &
                          d = 0.349_dbl

  REAL(dbl), PARAMETER :: c1 = 0.018897_dbl, &
                          c2 = -0.155240_dbl, &
                          c3 = 0.159068_dbl, &
                          c4 = -0.007953_dbl

  character(len=*), private, parameter :: moduleN='xc_cs1_new'

!!***
! *****************************************************************************

CONTAINS

  !!****f* xc_cs1_new/cs1_lda_info [1.0] *
  !!
  !!   NAME
  !!     cs1_lda_info
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
  !!*** **********************************************************************
  SUBROUTINE cs1_lda_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
         INTENT(inout), optional                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cs1_lda_info', &
         routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
      reference = "N.C. Handy and A.J. Cohen, J. Chem. Phys., 116, 5411 (2002)"//&
           "{LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "CS1: Handy improved LYP correlation energy functional {LDA}"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.true.
       needs%rho_1_3=.true.
       needs%norm_drho=.true.
    END IF
    if (present(max_deriv)) max_deriv=3

  END SUBROUTINE cs1_lda_info
  !***************************************************************************

  !!****f* xc_cs1_new/cs1_lsd_info [1.0] *
  !!
  !!   NAME
  !!     cs1_lsd_info
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
  !!*** **********************************************************************
  SUBROUTINE cs1_lsd_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
         INTENT(inout), optional                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cs1_lsd_info', &
         routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
      reference = "N.C. Handy and A.J. Cohen, J. Chem. Phys., 116, 5411 (2002)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "CS1: Handy improved LYP correlation energy functional"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.true.
       needs%rho_spin_1_3=.true.
       needs%norm_drho_spin=.true.
    END IF
    if (present(max_deriv)) max_deriv=1

  END SUBROUTINE cs1_lsd_info
  !***************************************************************************

! *****************************************************************************

  SUBROUTINE cs1_init ( cutoff, debug )

    REAL(dbl), INTENT(IN)                    :: cutoff
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    eps_rho = cutoff
    CALL set_util ( cutoff )

    IF ( PRESENT ( debug ) ) THEN
       debug_flag = debug
    ELSE
       debug_flag = .FALSE.
    END IF

    fsig = 2._dbl**f13

  END SUBROUTINE cs1_init

! *****************************************************************************

  SUBROUTINE cs1_lda_eval (rho_set,deriv_set,order,error)
    TYPE(xc_rho_set_type), pointer :: rho_set
    TYPE(xc_derivative_set_type), pointer :: deriv_set
    INTEGER, INTENT(in) :: order
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cs1_lda_eval',&
         routineP=moduleN//':'//routineN
    INTEGER :: npoints, handle
    integer, DIMENSION(:,:), pointer :: bo
    REAL(kind=dbl),DIMENSION(:,:,:), POINTER :: rho, rho13, grho, e_0,e_rho,e_ndrho,&
         e_rho_rho, e_rho_ndrho, e_ndrho_ndrho,&
         e_rho_rho_rho, e_rho_rho_ndrho, e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho
    TYPE(xc_derivative_type), POINTER :: deriv
    REAL(kind=dbl) :: epsilon_rho
    INTEGER :: m

    call timeset(routineN,"I","",handle)
    failure=.false.
    NULLIFY(bo,rho, rho13, grho, e_0,e_rho,e_ndrho,&
         e_rho_rho, e_rho_ndrho, e_ndrho_ndrho,&
         e_rho_rho_rho, e_rho_rho_ndrho, e_rho_ndrho_ndrho, &
         e_ndrho_ndrho_ndrho)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho_1_3=rho13,rho=rho,&
            norm_drho=grho,local_bounds=bo,rho_cutoff=epsilon_rho,&
            error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       m=ABS(order)
       call cs1_init(epsilon_rho)

       IF ( order>=0 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          call xc_derivative_get(deriv,deriv_data=e_0,error=error)

          CALL cs1_u_0 ( rho, grho, rho13, e_0,npoints, error)
       END IF
       IF ( order>=1 .OR. order==-1 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)

          CALL cs1_u_1 ( rho, grho, rho13, e_rho, e_ndrho,&
               npoints,error)
       END IF
       IF ( order>=2.OR.order==-2 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(norm_drho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)

          CALL cs1_u_2 ( rho, grho, rho13, e_rho_rho, e_rho_ndrho, &
               e_ndrho_ndrho, npoints, error )
       END IF
       IF ( order>=3.OR.order==-3 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(rho)(rho)(norm_drho)",allocate_deriv=.true.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_rho_rho_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(rho)(norm_drho)(norm_drho)",allocate_deriv=.true.,error=error)
          call xc_derivative_get(deriv,deriv_data=e_rho_ndrho_ndrho,error=error)
          deriv => xc_dset_get_derivative(deriv_set,&
               "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,&
               error=error)
          call xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)

          CALL cs1_u_3 ( rho, grho, rho13, e_rho_rho_rho, e_rho_rho_ndrho,&
               e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints, error)
       END IF
       IF ( order>3.OR.order<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

    END IF
    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE cs1_lda_eval

! *****************************************************************************

  SUBROUTINE cs1_lsd_eval (rho_set,deriv_set,order,error)
    TYPE(xc_rho_set_type), pointer :: rho_set
    TYPE(xc_derivative_set_type), pointer :: deriv_set
    INTEGER, INTENT(in) :: order
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cs1_lsd_eval',&
         routineP=moduleN//':'//routineN
    INTEGER :: npoints, handle
    integer, DIMENSION(:,:), pointer :: bo
    REAL(kind=dbl),DIMENSION(:,:,:), POINTER :: rhoa,rhob,norm_drhoa,&
         rhoa_1_3, rhob_1_3,norm_drhob, e_0,e_rhoa,e_rhob,e_ndrhoa,e_ndrhob
    TYPE(xc_derivative_type), POINTER :: deriv
    REAL(kind=dbl) :: epsilon_rho

    call timeset(routineN,"I","",handle)
    failure=.false.
    nullify(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
            rhoa_1_3=rhoa_1_3,rhob_1_3=rhob_1_3,&
            norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,&
            local_bounds=bo,rho_cutoff=epsilon_rho, error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       call cs1_init(epsilon_rho)

       IF ( order>=0 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          call xc_derivative_get(deriv,deriv_data=e_0,error=error)

          CALL cs1_ss_0 ( rhoa, rhob, norm_drhoa, norm_drhob,&
               rhoa_1_3, rhob_1_3, e_0,npoints,error)
       END IF
       IF ( order>=1.OR.order==-1 ) THEN
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

          CALL cs1_ss_1 ( rhoa, rhob, norm_drhoa, norm_drhob,&
               rhoa_1_3, rhob_1_3, e_rhoa, e_rhob, e_ndrhoa, e_ndrhob,&
               npoints,error)
       END IF
       if ( order>1 .or. order<1 ) then
          call cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF
    END IF
    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE cs1_lsd_eval

! *****************************************************************************

  SUBROUTINE cs1_u_0 ( rho, grho, r13, e_0, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: e_0
    INTEGER, intent(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: c2p, c3p, c4p, cp, dp, F1, &
                                                F2, F3, F4, g, oc, ocp, od, &
                                                odp, r, r3, x

!   ---------------------------------------------------------------------------

    dp = fsig*d
    cp = c*fsig*fsig
    c2p = c2*fsig**4
    c3p = -c3*0.25_dbl
    c4p = -c4*0.25_dbl

!$omp parallel do private(ip,r,r3,g,x,odp,ocp,od,oc,f1,f2,f3,f4)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         r = rho(ip)
         r3 = r13(ip)
         g = grho(ip)
         x = g/(r*r3)
         odp = 1._dbl/(r3+dp)
         ocp = 1._dbl/(r*r*r3*r3 + cp*g*g)
         od = 1._dbl/(r3+d)
         oc = 1._dbl/(r*r*r3*r3 + c*g*g)
         F1 = c1*r*r3*odp
         F2 = c2p*g**4*r3*r*odp*ocp*ocp
         F3 = c3p*r*r3*od
         F4 = c4p*g**4*r3*r*od*oc*oc
         e_0(ip) = e_0(ip) + F1 + F2 + F3 + F4
      END IF

    END DO

  END SUBROUTINE cs1_u_0

! *****************************************************************************

  SUBROUTINE cs1_u_1 ( rho, grho, r13, e_rho, e_ndrho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)  :: e_rho, e_ndrho
    INTEGER, intent(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: c2p, c3p, c4p, cp, dF1, dF3, &
                                                dgF2, dgF4, dp, drF2, drF4, &
                                                g, oc, ocp, od, odp, r, r3

!   ---------------------------------------------------------------------------

    dp = fsig*d
    cp = c*fsig*fsig
    c2p = c2*fsig**4
    c3p = -c3*0.25_dbl
    c4p = -c4*0.25_dbl

!$omp parallel do private(ip,r,r3,g,odp,ocp,od,oc,df1,drf2,dgf2,df3,drf4,dgf4)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         r = rho(ip)
         r3 = r13(ip)
         g = grho(ip)
         odp = 1._dbl/(r3+dp)
         ocp = 1._dbl/(r*r*r3*r3 + cp*g*g)
         od = 1._dbl/(r3+d)
         oc = 1._dbl/(r*r*r3*r3 + c*g*g)
         dF1 = c1*f13*r3*(3*r3+4*dp)*odp*odp
         drF2 = -f13*c2p*g**4*r3*(13*r**3-3*r3*cp*g*g+12*r*r*r3*r3*dp-&
                 4*dp*cp*g*g)*odp**2*ocp**3
         dgF2 = 4*c2p*g**3*r**4*odp*ocp**3
         dF3 = f13*r3*c3p*(3*r3+4*d)*od*od
         drF4 = -f13*c4p*g**4*r3*(13*r**3-3*r3*c*g*g+12*r*r*r3*r3*d-&
                 4*d*c*g*g)*od**2*oc**3
         dgF4 = 4*c4p*g**3*r**4*od*oc**3
         e_rho(ip) = e_rho(ip) + dF1 + drF2 + dF3 + drF4
         e_ndrho(ip) = e_ndrho(ip) + dgF2 + dgF4
      END IF

    END DO

  END SUBROUTINE cs1_u_1

! *****************************************************************************

  SUBROUTINE cs1_u_2 ( rho, grho, r13, e_rho_rho, e_rho_ndrho, e_ndrho_ndrho,&
       npoints, error)

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)  :: e_rho_rho, e_rho_ndrho,&
         e_ndrho_ndrho
    INTEGER, intent(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: c2p, c3p, c4p, cp, d2F1, &
                                                d2F3, d2gF2, d2gF4, d2rF2, &
                                                d2rF4, dp, drgF2, drgF4, g, &
                                                oc, ocp, od, odp, r, r3

!   ---------------------------------------------------------------------------

    dp = fsig*d
    cp = c*fsig*fsig
    c2p = c2*fsig**4
    c3p = -c3*0.25_dbl
    c4p = -c4*0.25_dbl

!$omp parallel do private(ip,r,r3,g,odp,ocp,od,oc,d2f1,d2rf2,drgf2,d2f3) &
!$omp             private(d2rf4,drgf4,d2gf4)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         r = rho(ip)
         r3 = r13(ip)
         g = grho(ip)
         odp = 1._dbl/(r3+dp)
         ocp = 1._dbl/(r*r*r3*r3 + cp*g*g)
         od = 1._dbl/(r3+d)
         oc = 1._dbl/(r*r*r3*r3 + c*g*g)
         d2F1 = c1*f23*f13*dp*r3/r*(r3+2*dp)*odp**3
         d2rF2 = c2p*f13*f23*g**4*r3/r * (193*dp*r**5*r3*r3+90*dp*dp*r**5*r3 &
                 -88*g*g*cp*r**3*r3-100*dp*dp*cp*g*g*r*r*r3*r3 &
                 +2*dp*dp*cp*cp*g**4-190*g*g*r**3*cp*dp+g**4*r3*cp*cp*dp &
                 +104*r**6)*odp**3*ocp**4
         drgF2 = c2p*f43*g**3*r*r*r3*(-13*r**3*r3*r3+11*cp*r*g*g-12*dp*r**3*r3 &
                 +12*r3*r3*dp*cp*g*g)*odp*odp*ocp**4
         d2gF2 = -12*c2p*g*g*r**4*(cp*g*g-r*r*r3*r3)*odp*ocp**4
         d2F3 = f13*f23*c3p*d*r3/r*(r3+2*d)*od**3
         d2rF4 = c4p*f13*f23*g**4*r3/r * (193*d*r**5*r3*r3+90*d*d*r**5*r3 &
                 -88*g*g*c*r**3*r3-100*d*d*c*g*g*r*r*r3*r3 &
                 +2*d*d*c*c*g**4-190*g*g*r**3*c*d+g**4*r3*c*c*d &
                 +104*r**6)*od**3*oc**4
         drgF4 = c4p*f43*g**3*r*r*r3*(-13*r**3*r3*r3+11*c*r*g*g-12*d*r**3*r3 &
                 +12*r3*r3*d*c*g*g)*od*od*oc**4
         d2gF4 = -12*c4p*g*g*r**4*(c*g*g-r*r*r3*r3)*od*oc**4
         e_rho_rho(ip) = e_rho_rho(ip) + d2F1 + d2rF2 + d2F3 + d2rF4
         e_rho_ndrho(ip) = e_rho_ndrho(ip) + drgF2 + drgF4
         e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + d2gF2 + d2gF4
      END IF

    END DO

  END SUBROUTINE cs1_u_2

! *****************************************************************************

  SUBROUTINE cs1_u_3 ( rho, grho, r13, e_rho_rho_rho, e_rho_rho_ndrho,&
       e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints, error)

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)  :: e_rho_rho_rho, e_rho_rho_ndrho,&
       e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho
    INTEGER, intent(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error

    INTEGER                                  :: ip
    REAL(dbl) :: c2p, c3p, c4p, cp, d2rgF2, d2rgF4, d3F1, d3F3, d3gF2, d3gF4, &
      d3rF2, d3rF4, dp, dr2gF2, dr2gF4, g, oc, ocp, od, odp, r, r3, t1, t10, &
      t11, t13, t15, t16, t17, t19, t2, t22, t23, t26, t29, t3, t32, t33, &
      t37, t4, t44, t45, t50, t51, t52, t58, t6, t61, t7, t74, t76, t77, t8, &
      t80, t81, t82, t9

!   ---------------------------------------------------------------------------

    dp = fsig*d
    cp = c*fsig*fsig
    c2p = c2*fsig**4
    c3p = -c3*0.25_dbl
    c4p = -c4*0.25_dbl

!$omp parallel do private(ip,r,r3,g,odp,ocp,od,oc,d3f1,d3rF2,d2rgF2,dr2gF2,d3gF2,d3F3) &
!$omp             private(d3rF4,d2rgF4,dr2gF4,d3gF4) &
!$omp             private( t1, t2, t3, t4, t6, t7, t8, t9, t10, t11, t13, t15, t16) &
!$omp             private( t17, t19, t22, t23, t26, t29, t32, t33, t37, t44, t45, t50) &
!$omp             private( t51, t52, t58, t61, t74, t76, t77, t80, t81, t82)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         r = rho(ip)
         r3 = r13(ip)
         g = grho(ip)
         odp = 1._dbl/(r3+dp)
         ocp = 1._dbl/(r*r*r3*r3 + cp*g*g)
         od = 1._dbl/(r3+d)
         oc = 1._dbl/(r*r*r3*r3 + c*g*g)
         d3F1=-c1*f23*f13*f13*dp*r3/(r*r)*(11*dp*r3+4*dp*dp+4*r/r3)*odp**4
         t1 = g**2; t2 = t1**2; t3 = r3; t4 = t3**2; t8 = dp**2; t9 = t8*dp
         t10 = cp**2; t11 = t10*cp; t13 = t2*t1; t16 = r**2; t17 = t4*t16
         t19 = t10*t2; t22 = t16**2; t23 = t22**2; t32 = t22*t16; t37 = t16*r
         t58 = t22*r; t61 = cp*t1
         t74 = 4*t9*t11*t13+668*t17*t9*t19+5524*t4*t23*dp+5171*t3*t23*t8+ &
               1620*t23*t9-3728*t3*t32*cp*t1+440*t4*t37*t10*t2+1500*t2*t3*t37*dp*t10 &
               +4*t13*t4*dp*t11+1737*t37*t8*t19+11*t3*t8*t11*t13-3860*t3*t58*t9*t61+ &
               1976*t23*r-11535*t4*t58*t8*t61-11412*t1*t32*cp*dp
         t76 = (t3+dp)**2; t77 = t76**2; t80 = t17+t61; t81 = t80**2; t82 = t81**2
         d3rF2 = -f23*f13*f13*c2p*t2/t4/r*t74/t77/t82/t80
         t4 = t3*r; t6 = r**2; t7 = t6**2; t8 = t7*t6; t9 = dp**2; t15 = t1**2
         t17 = cp**2; t23 = t3**2; t26 = t6*r; t29 = cp*t1; t33 = t17*t15; t44 = t3+dp
         t45 = t44**2; t50 = t23*t6+t29; t51 = t50**2; t52 = t51**2
         d2rgF2 = c2p*f23*f43*t1*g*t4*(90*t8*t9+193*t3*t8*dp+44*t15*t4*t17-236*t1 &
                  *t7*cp+104*t23*t8-240*t3*t26*t9*t29+54*t23*t9*t33-478*t23*t26*dp*t29 &
                  +97*r*dp*t33)/t45/t44/t52/t50
         dr2gF2 = -4*c2p*g*g*r*r*r3*(-40*r**3*r3*dp*cp*g*g+12*r3*r3*dp*cp*cp*g**4 &
                  +13*r**6*r3-40*r**3*r3*r3*cp*g*g+11*r*cp*cp*g**4+12*r**6*dp) &
                  *odp*odp*ocp**5
         d3gF2 = c2p*24*g*r**3*r3*(r**6-5*cp*g*g*r**3*r3+2*cp*cp*g**4*r3*r3)*odp*ocp**5
         d3F3 = -f23*f13*f13*c3p*d*r3/(r*r)*(11*d*r3+4*d*d+4*r3*r3)*od**4
         t1 = g**2; t2 = t1**2; t3 = r3; t4 = t3**2; t8 = d**2; t9 = t8*d
         t10 = c**2; t11 = t10*c; t13 = t2*t1; t16 = r**2; t17 = t4*t16
         t19 = t10*t2; t22 = t16**2; t23 = t22**2; t32 = t22*t16; t37 = t16*r
         t58 = t22*r; t61 = c*t1
         t74 = 4*t9*t11*t13+668*t17*t9*t19+5524*t4*t23*d+5171*t3*t23*t8+ &
               1620*t23*t9-3728*t3*t32*c*t1+440*t4*t37*t10*t2+1500*t2*t3*t37*d*t10 &
               +4*t13*t4*d*t11+1737*t37*t8*t19+11*t3*t8*t11*t13-3860*t3*t58*t9*t61+ &
               1976*t23*r-11535*t4*t58*t8*t61-11412*t1*t32*c*d
         t76 = (t3+d)**2; t77 = t76**2; t80 = t17+t61; t81 = t80**2; t82 = t81**2
         d3rF4 = -f23*f13*f13*c4p*t2/t4/r*t74/t77/t82/t80
         t4 = t3*r; t6 = r**2; t7 = t6**2; t8 = t7*t6; t9 = d**2; t15 = t1**2
         t17 = c**2; t23 = t3**2; t26 = t6*r; t29 = c*t1; t33 = t17*t15; t44 = t3+d
         t45 = t44**2; t50 = t23*t6+t29; t51 = t50**2; t52 = t51**2
         d2rgF4 = c4p*f23*f43*t1*g*t4*(90*t8*t9+193*t3*t8*d+44*t15*t4*t17-236*t1 &
                  *t7*c+104*t23*t8-240*t3*t26*t9*t29+54*t23*t9*t33-478*t23*t26*d*t29 &
                  +97*r*d*t33)/t45/t44/t52/t50
         dr2gF4 = -4*c4p*g*g*r*r*r3*(-40*r**3*r3*d*c*g*g+12*r3*r3*d*c*c*g**4 &
                  +13*r**6*r3-40*r**3*r3*r3*c*g*g+11*r*c*c*g**4+12*r**6*d) &
                  *od*od*oc**5
         d3gF4 = c4p*24*g*r**3*r3*(r**6-5*c*g*g*r**3*r3+2*c*c*g**4*r3*r3)*od*oc**5
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip) + d3F1 + d3rF2 + d3F3 + d3rF4
         e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip) + d2rgF2 + d2rgF4
         e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip) + dr2gF2 + dr2gF4
         e_ndrho_ndrho_ndrho(ip) = e_ndrho_ndrho_ndrho(ip) + d3gF2 + d3gF4
      END IF

    END DO

  END SUBROUTINE cs1_u_3

! *****************************************************************************

  SUBROUTINE cs1_ss_0 ( rhoa, rhob, grhoa, grhob, r13a, r13b, e_0,&
       npoints,error)

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rhoa, rhob, grhoa, grhob, &
                                                r13a, r13b
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: e_0
    INTEGER, intent(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: F1a, F1b, F2a, F2b, ga, gb, &
                                                oca, ocb, oda, odb, r3a, r3b, &
                                                ra, rb, xa, xb
    CHARACTER(len=*), PARAMETER :: routineN="cs1_ss_0",&
         routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    !FM IF ( .NOT. debug_flag ) CALL stop_program("cs1_ss_0","Routine not tested")
    CALL cp_assert(.FALSE.,cp_warning_level, cp_assertion_failed,routineP,&
         "not tested!"//&
         CPSourceFileRef,&
         error)

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         F1a = 0._dbl
         F2a = 0._dbl
      ELSE
         ra = rhoa(ip)
         r3a = r13a(ip)
         ga = grhoa(ip)
         xa = ga/(ra*r3a)
         oda = 1._dbl/(r3a+d)
         oca = 1._dbl/(ra*ra*r3a*r3a + c*ga*ga)
         F1a = c1*ra*r3a*oda
         F2a = c2*ga**4*r3a*ra*oda*oca*oca
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         F1b = 0._dbl
         F2b = 0._dbl
      ELSE
         rb = rhob(ip)
         r3b = r13b(ip)
         gb = grhob(ip)
         xb = gb/(rb*r3b)
         odb = 1._dbl/(r3b+d)
         ocb = 1._dbl/(rb*rb*r3b*r3b + c*gb*gb)
         F1b = c1*rb*r3b*odb
         F2b = c2*gb**4*r3b*rb*odb*ocb*ocb
      END IF

      e_0(ip) = e_0(ip) + F1a + F1b + F2a + F2b

    END DO

  END SUBROUTINE cs1_ss_0

! *****************************************************************************

  SUBROUTINE cs1_ss_1 ( rhoa, rhob, grhoa, grhob, r13a, r13b, e_rhoa,&
       e_rhob, e_ndrhoa, e_ndrhob, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rhoa, rhob, grhoa, grhob, &
                                                r13a, r13b
    REAL(dbl), DIMENSION(*), INTENT(INOUT)  :: e_rhoa,&
       e_rhob, e_ndrhoa, e_ndrhob
    INTEGER, intent(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error


    INTEGER                                  :: ip
    REAL(dbl)                                :: dF1a, dF1b, dgF2a, dgF2b, &
                                                drF2a, drF2b, ga, gb, oca, &
                                                ocb, oda, odb, r3a, r3b, ra, &
                                                rb
    CHARACTER(len=*), PARAMETER :: routineN="cs1_ss_1",&
         routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert(.FALSE.,cp_warning_level, cp_assertion_failed,routineP,&
         "not tested!"//&
         CPSourceFileRef,&
         error)

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         dF1a = 0._dbl
         drF2a = 0._dbl
         dgF2a = 0._dbl
      ELSE
         ra = rhoa(ip)
         r3a = r13a(ip)
         ga = grhoa(ip)
         oda = 1._dbl/(r3a+d)
         oca = 1._dbl/(ra*ra*r3a*r3a + c*ga*ga)
         dF1a = c1*f13*r3a*(3*r3a+4*d)*oda*oda

         drF2a = -f13*c2*ga**4*r3a*(13*ra**3-3*r3a*c*ga*ga+12*ra*ra*r3a*r3a*d-&
                 4*d*c*ga*ga)*oda**2*oca**3

         dgF2a = 4*c2*ga**3*ra**4*oda*oca**3

      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         dF1b = 0.0_dbl
         drF2b = 0._dbl
         dgF2b = 0._dbl
      ELSE
         rb = rhob(ip)
         r3b = r13b(ip)
         gb = grhob(ip)
         odb = 1._dbl/(r3b+d)
         ocb = 1._dbl/(rb*rb*r3b*r3b + c*gb*gb)
         dF1b = c1*f13*r3b*(3*r3b+4*d)*odb*odb

         drF2b = -f13*c2*gb**4*r3b*(13*rb**3-3*r3b*c*gb*gb+12*rb*rb*r3b*r3b*d-&
                 4*d*c*gb*gb)*odb**2*ocb**3

         dgF2b = 4*c2*gb**3*rb**4*odb*ocb**3

      END IF

      e_rhoa(ip) = e_rhoa(ip) + dF1a + drF2a
      e_ndrhoa(ip) = e_ndrhoa(ip) + dgF2a 
      e_rhob(ip) = e_rhob(ip) + dF1b + drF2b
      e_ndrhob(ip) = e_ndrhob(ip) + dgF2b

    END DO

  END SUBROUTINE cs1_ss_1

! *****************************************************************************

END MODULE xc_cs1_new

! *****************************************************************************
