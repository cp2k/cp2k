!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_tpss_u [1.0] *
!!
!!   NAME
!!     xc_tpss_u
!!
!!   FUNCTION
!!     Calculate the TPSS Functional (spin unpolarized)
!!
!!   AUTHOR
!!     JGH (04.01.2004)
!!
!!   MODIFICATION HISTORY
!!
!******************************************************************************

MODULE xc_tpss_u

! *****************************************************************************

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type
  USE xc_tpss_util,                    ONLY: pz_u,&
                                             tpbec,&
                                             tpbeca,&
                                             tpssux

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="xc_tpss_u"
  REAL(dbl),PARAMETER  :: small = 1.e-14_dbl
  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f83 = 8._dbl*f13, &
                          f1081 = 10._dbl/81._dbl

  PUBLIC :: tpss_lda_info, tpss_lda_eval

!!***
! *****************************************************************************

CONTAINS


!!****f* xc_tpss_u/tpss_lda_info [1.0] *
!!
!!   NAME
!!     tpss_lda_info
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
SUBROUTINE tpss_lda_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_info', &
      routineP = moduleN//':'//routineN

  IF ( PRESENT ( reference ) ) THEN
     reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LDA version}"
  END IF
  IF ( PRESENT ( shortform ) ) THEN
     shortform = "TPSS meta-GGA functional (LDA)"
  END IF
  IF (PRESENT(needs)) THEN
     needs%rho=.TRUE.
     needs%norm_drho=.TRUE.
     needs%tau=.TRUE.
  END IF
  IF (PRESENT(max_deriv)) max_deriv=1

END SUBROUTINE tpss_lda_info
!***************************************************************************

!!****f* xc_tpss_u/tpss_lda_eval [1.0] *
!!
!!   NAME
!!     tpss_lda_eval
!!
!!   FUNCTION
!!     evaluates the tpss functional (spin unpolarized)
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE tpss_lda_eval(rho_set,deriv_set,grad_deriv,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dbl)                           :: epsilon_drho, epsilon_rho
    REAL(kind=dbl), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_ndrho, e_rho, e_tau, &
                                                norm_drho, rho, tau
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo,e_0, e_ndrho, e_rho, e_tau, norm_drho, rho, tau,deriv)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,tau=tau, norm_drho=norm_drho,&
            can_return_null=.FALSE.,local_bounds=bo,rho_cutoff=epsilon_rho,&
            drho_cutoff=epsilon_drho,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)

       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)

       deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)

       deriv => xc_dset_get_derivative(deriv_set,"(tau)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_tau,error=error)

       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL tpss_u_a(rho=rho, norm_drho=norm_drho,tau=tau,&
            e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_tau=e_tau,&
            npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_drho=epsilon_drho, error=error)
    END IF
    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE tpss_lda_eval
  !***************************************************************************

!!****f* xc_tpss_u/tpss_u_a [1.0] *
!!
!!   NAME
!!     tpss_u_a
!!
!!   FUNCTION
!!     calculates the tpss functional on 1d arrays
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - npoints: the size of the 1d arrays
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     The others should be pretty self explaining (otherwise see
!!     tpss_lda_eval and the module xc_derivative_desc)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE tpss_u_a(rho, norm_drho,tau,e_0,e_rho,e_ndrho,e_tau,&
       npoints,epsilon_rho, epsilon_drho, error)
    REAL(kind=dbl), DIMENSION(*), INTENT(in) :: rho, norm_drho, tau
    REAL(kind=dbl), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rho, e_ndrho, e_tau
    INTEGER, INTENT(in)                      :: npoints
    REAL(kind=dbl), INTENT(in)               :: epsilon_rho, epsilon_drho
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_u_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    REAL(dbl)                                :: sc, sx, v1c, v1x, v2c, v2x, &
                                                vttc, vttx

    !$omp parallel do default(none) shared(rho,norm_drho,&
    !$omp          tau,e_0, e_rho, e_ndrho, e_tau) &
    !$omp      private(ii,sx,v1x,v2x,vttx,sc,v1c,v2c,vttc)
    DO ii=1,npoints
       IF (rho(ii)>epsilon_rho) THEN
          CALL tpss_u(rho=rho(ii),grho=norm_drho(ii)**2,tau=tau(ii),&
               sx=sx,v1x=v1x,v2x=v2x,vttx=vttx,sc=sc,v1c=v1c,v2c=v2c,vttc=vttc)
          e_0(ii)=e_0(ii)+sx+sc
          e_rho(ii)=e_rho(ii)+v1x+v1c
          e_ndrho(ii)=e_ndrho(ii)+v2x+v2c ! * norm_drho(ii) probably
          e_tau(ii)=e_tau(ii)+vttx+vttc
       END IF
    END DO
  END SUBROUTINE tpss_u_a
!***************************************************************************

! *****************************************************************************

      SUBROUTINE tpss_u(rho,grho,tau,sx,v1x,v2x,vttx,sc,v1c,v2c,vttc)

    REAL(dbl), INTENT(IN)                    :: rho, grho, tau
    REAL(dbl), INTENT(OUT)                   :: sx, v1x, v2x, vttx, sc, v1c, &
                                                v2c, vttc

        IF ( rho > small ) THEN
!         Exchange part
          CALL tpssux(rho,grho,tau,sx,v1x,v2x,vttx)

!         Correlation
          CALL tpssuc(rho,grho,tau,sc,v1c,v2c,vttc)
        ELSE
          sx   = 0._dbl
          v1x  = 0._dbl
          v2x  = 0._dbl
          vttx = 0._dbl
          sc   = 0._dbl
          v1c  = 0._dbl
          v2c  = 0._dbl
          vttc = 0._dbl
        END IF

      END SUBROUTINE tpss_u

! *****************************************************************************

      SUBROUTINE tpssuc(rho,grho,tau,sc,v1c,v2c,vttc)
!       correlation part of TPSS functional 

    REAL(dbl), INTENT(IN)                    :: rho, grho, tau
    REAL(dbl), INTENT(OUT)                   :: sc, v1c, v2c, vttc

    REAL(dbl), PARAMETER                     :: d = 2.8e0_dbl

    INTEGER                                  :: iflg
    REAL(dbl)                                :: dedg, dedr, dedt, dedz, dzdg, &
                                                dzdr, dzdt, e, edgt, edrt, &
                                                op, ro, rs, z

! -----------------------------------------------------------------------------

        IF ( ABS(tau) > small .AND. grho > small ) THEN
          z    = 0.125_dbl*grho/rho/tau
          dzdr = -z/rho
          dzdg = 2._dbl*z/SQRT(grho)
          dzdt = -z/tau
          CALL revpkzb(rho,grho,z,e,dedr,dedg,dedz)
          dedt = dedz*dzdt
          edrt = dedr + dedz*dzdr
          edgt = dedg + dedz*dzdg
          op  = 1._dbl + d * e * z**3
          ro  = rho * e * d * z**2
          sc  = rho * e * op
          v1c = (e + rho*edrt)*op + ro * (edrt*z+3._dbl*e*dzdr)
          v2c = rho*edgt*op + ro * (edgt*z+3._dbl*e*dzdg)
          v2c = v2c/SQRT(grho)
          vttc = rho*dedt*op + ro * (dedt*z+3._dbl*e*dzdt)
        ELSE
          rs   = (3._dbl/(4._dbl*pi*rho))**(1._dbl/3._dbl)
          iflg = 2
          IF(rs.lt.1.0_dbl) iflg=1
          CALL pz_u(rs,e,v1c,iflg)
          sc   = rho * e
          v2c  = 0._dbl
          vttc = 0._dbl
        ENDIF

      END SUBROUTINE tpssuc

! *****************************************************************************

      SUBROUTINE revpkzb(rho,grho,z,e,dedr,dedg,dedz)
!       Revised PKZB tau correlation functional

    REAL(dbl), INTENT(IN)                    :: rho, grho, z
    REAL(dbl), INTENT(OUT)                   :: e, dedr, dedg, dedz

    REAL(dbl), PARAMETER                     :: c00 = 0.53_dbl, c0p = 1.53_dbl

    REAL(dbl)                                :: ecpbe, ecpbea, t1, tecpbe, &
                                                tv1c, tv2c, v1c, v1ca, v2c, &
                                                v2ca, z2

! -----------------------------------------------------------------------------

        CALL tpbec(rho,grho,ecpbe,v1c,v2c)
        CALL tpbeca(0.5_dbl*rho,0.25_dbl*grho,ecpbea,v1ca,v2ca)
        IF ( ecpbea > ecpbe ) THEN
          tecpbe = ecpbea
          tv1c   = 0.5_dbl*v1ca
          tv2c   = 0.5_dbl*v2ca
        ELSE
          tecpbe = ecpbe
          tv1c   = v1c
          tv2c   = v2c
        END IF
        z2 = z*z
        t1 = 1._dbl + c00*z2
        e    = ecpbe*t1 - c0p*z2*tecpbe
        dedr = v1c*t1 - c0p*z2*tv1c
        dedg = v2c*t1 - c0p*z2*tv2c
        dedz = 2._dbl*z * (ecpbe*c00 - c0p*tecpbe)

      END SUBROUTINE revpkzb

! *****************************************************************************

END MODULE xc_tpss_u

! *****************************************************************************
