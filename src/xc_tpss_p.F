!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_tpss_p [1.0] *
!!
!!   NAME
!!     xc_tpss_p
!!
!!   FUNCTION
!!     Calculate the TPSS Functional (spin polarized)
!!
!!   AUTHOR
!!     JGH (04.01.2004)
!!
!!   MODIFICATION HISTORY
!!
!******************************************************************************

MODULE xc_tpss_p

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
  USE xc_tpss_util,                    ONLY: ccfun,&
                                             pz_p,&
                                             tpbeca,&
                                             tpbecs,&
                                             tpssux

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="xc_tpss_p"
  REAL(dbl),PARAMETER  :: small = 1.e-14_dbl
  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f83 = 8._dbl*f13, &
                          f1081 = 10._dbl/81._dbl

  PUBLIC :: tpss_lsd_info, tpss_lsd_eval

!!***
! *****************************************************************************

CONTAINS

!!****f* xc_tpss_p/tpss_lsd_info [1.0] *
!!
!!   NAME
!!     tpss_lsd_info
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
SUBROUTINE tpss_lsd_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_info', &
      routineP = moduleN//':'//routineN

  IF ( PRESENT ( reference ) ) THEN
     reference = "J. Tao, J.P.Perdew, V.N.Staroverov, E.Scuseria PRL, 91, 146401 (2003) {LSD version}"
  END IF
  IF ( PRESENT ( shortform ) ) THEN
     shortform = "TPSS meta-GGA functional (LSD)"
  END IF
  IF (PRESENT(needs)) THEN
     needs%rho_spin=.TRUE.
     needs%norm_drho_spin=.TRUE.
     needs%drhoa_drhob=.TRUE.
     needs%tau_spin=.TRUE.
  END IF
  IF (PRESENT(max_deriv)) max_deriv=1

END SUBROUTINE tpss_lsd_info
!***************************************************************************

!!****f* xc_tpss_p/tpss_lsd_eval [1.0] *
!!
!!   NAME
!!     tpss_lsd_eval
!!
!!   FUNCTION
!!     evaluates the tpss functional (spin polarized)
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
  SUBROUTINE tpss_lsd_eval(rho_set,deriv_set,grad_deriv,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dbl)                           :: epsilon_drho, epsilon_rho
    REAL(kind=dbl), DIMENSION(:, :, :), POINTER :: drhoa_drhob, e_0, &
      e_drhoa_drhob, e_ndrhoa, e_ndrhob, e_rhoa, e_rhob, e_tau_a, e_tau_b, &
      norm_drhoa, norm_drhob, rhoa, rhob, tau_a, tau_b
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo,e_0, e_rhoa, e_rhob, e_ndrhoa,&
      e_ndrhob,e_tau_a,e_tau_b,e_drhoa_drhob,&
      rhoa,rhob,norm_drhoa,norm_drhob,drhoa_drhob, tau_a,tau_b,deriv)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,tau_a=tau_a,&
            tau_b=tau_b,norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,&
            drhoa_drhob=drhoa_drhob,&
            can_return_null=.FALSE.,local_bounds=bo,rho_cutoff=epsilon_rho,&
            drho_cutoff=epsilon_drho,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)

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
       deriv => xc_dset_get_derivative(deriv_set,"(drhoa_drhob)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_drhoa_drhob,error=error)

       deriv => xc_dset_get_derivative(deriv_set,"(tau_a)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_tau_a,error=error)
       deriv => xc_dset_get_derivative(deriv_set,"(tau_b)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_tau_b,error=error)

       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL tpss_p_a(rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa,&
            norm_drhob=norm_drhob,drhoa_drhob=drhoa_drhob,&
            tau_a=tau_a,tau_b=tau_b,&
            e_0=e_0,e_rhoa=e_rhoa,e_rhob=e_rhob,e_ndrhoa=e_ndrhoa,&
            e_ndrhob=e_ndrhob,e_drhoa_drhob=e_drhoa_drhob,&
            e_tau_a=e_tau_a, e_tau_b=e_tau_b,&
            npoints=npoints,epsilon_rho=epsilon_rho,&
            epsilon_drho=epsilon_drho, error=error)
    END IF
    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE tpss_lsd_eval
  !***************************************************************************

!!****f* xc_tpss_p/tpss_p_a [1.0] *
!!
!!   NAME
!!     tpss_p_a
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
!!     tpss_lsd_eval and the module xc_derivative_desc)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE tpss_p_a(rhoa, rhob, norm_drhoa, norm_drhob,drhoa_drhob, tau_a,&
       tau_b,e_0,e_rhoa,e_rhob,e_tau_a,e_tau_b,e_ndrhoa, e_ndrhob,&
       e_drhoa_drhob, npoints,epsilon_rho,epsilon_drho, error)
    REAL(kind=dbl), DIMENSION(*), INTENT(in) :: rhoa, rhob, norm_drhoa, &
                                                norm_drhob, drhoa_drhob, &
                                                tau_a, tau_b
    REAL(kind=dbl), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rhoa, e_rhob, e_tau_a, &
                                                e_tau_b, e_ndrhoa, e_ndrhob, &
                                                e_drhoa_drhob
    INTEGER, INTENT(in)                      :: npoints
    REAL(kind=dbl), INTENT(in)               :: epsilon_rho, epsilon_drho
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tpss_p_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    REAL(dbl)                                :: sc, sx, v1ca, v1cb, v1xa, &
                                                v1xb, v2ca, v2cab, v2cb, &
                                                v2xa, v2xab, v2xb, vttca, &
                                                vttcb, vttxa, vttxb

    !$omp parallel do default(none) shared(rhoa,rhob,norm_drhoa,&
    !$omp         norm_drhob,drhoa_drhob, tau_a,tau_b,&
    !$omp         e_0, e_rhoa, e_rhob, e_ndrhoa, e_ndrhob, e_drhoa_drhob,&
    !$omp         e_tau_a, e_tau_b)
    !$omp     private(ii,sx,v1xa,v2xa,v1xb,v2xb,vttxa,vttxb,sc,v1ca,v2ca, &
    !$omp         v1cb,v2cb,v2xab,v2cab,vttca,vttcb)
    DO ii=1,npoints
       IF (rhoa(ii)+rhob(ii)>epsilon_rho) THEN
          CALL tpss_p(rhoa=rhoa(ii),rhob=rhob(ii),grhoaa=norm_drhoa(ii)**2,&
               grhobb=norm_drhob(ii)**2,grhoab=drhoa_drhob(ii),&
               taua=tau_a(ii),taub=tau_b(ii), &
               sx=sx,v1xa=v1xa,v2xa=v2xa,v1xb=v1xb,v2xb=v2xb,vttxa=vttxa,vttxb=vttxb,&
               sc=sc,v1ca=v1ca,v2ca=v2ca,v1cb=v1cb,v2cb=v2cb,v2xab=v2xab,v2cab=v2cab,&
               vttca=vttca,vttcb=vttcb)  
          e_0(ii)=e_0(ii)+sx+sc
          e_rhoa(ii)=e_rhoa(ii)+v1xa+v1ca
          e_rhob(ii)=e_rhob(ii)+v1xb+v1cb
          e_ndrhoa(ii)=e_ndrhoa(ii)+v2xa+v2ca ! * norm_drhoa(ii) probably
          e_ndrhob(ii)=e_ndrhob(ii)+v2xb+v2cb ! * norm_drhob(ii) probably
          e_drhoa_drhob(ii)=e_drhoa_drhob(ii)+v2xab+v2cab
          e_tau_a(ii)=e_tau_a(ii)+vttxa+vttca
          e_tau_b(ii)=e_tau_b(ii)+vttxb+vttcb
       END IF
  END DO
END SUBROUTINE tpss_p_a
!***************************************************************************

! *****************************************************************************

      SUBROUTINE tpss_p(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub, &
                        sx,v1xa,v2xa,v1xb,v2xb,vttxa,vttxb,sc,v1ca,v2ca, &
                        v1cb,v2cb,v2xab,v2cab,vttca,vttcb)

    REAL(dbl), INTENT(IN)                    :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, taua, taub
    REAL(dbl), INTENT(OUT)                   :: sx, v1xa, v2xa, v1xb, v2xb, &
                                                vttxa, vttxb, sc, v1ca, v2ca, &
                                                v1cb, v2cb, v2xab, v2cab, &
                                                vttca, vttcb

    REAL(dbl)                                :: grho, rho, sxx, tau, v1x, &
                                                v2x, vtt

        IF(rhoa > small .OR. rhob > small) THEN
!         Exchange part
!         alpha spin
          rho=2._dbl*rhoa
          grho=4._dbl*grhoaa
          tau=2._dbl*taua
          IF(rho > small) THEN
            CALL tpssux(rho,grho,tau,sxx,v1x,v2x,vtt)
            sx=0.5_dbl*sxx
            v1xa=v1x
            v2xa=2._dbl*v2x
            vttxa=vtt
          ELSE
            sx=0._dbl
            v1xa=0._dbl
            v2xa=0._dbl
            vttxa=0._dbl
          ENDIF
!         beta spin
          rho=2._dbl*rhob
          grho=4._dbl*grhobb
          tau=2._dbl*taub
          IF(RHO > SMALL) THEN
            CALL tpssux(rho,grho,tau,sxx,v1x,v2x,vtt)
            sx=sx + 0.5_dbl*sxx
            v1xb=v1x
            v2xb=2._dbl*v2x
            vttxb=vtt
          ELSE
            v1xb=0._dbl
            v2xb=0._dbl
            vttxb=0._dbl
          ENDIF
          v2xab=0._dbl
!         Correlation
          CALL tpsspc(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
                      sc,v1ca,v2ca,v1cb,v2cb,v2cab,vttca,vttcb)
        ELSE
          sx=0._dbl
          sc=0._dbl
          v1xa=0._dbl
          v1xb=0._dbl
          v2xa=0._dbl
          v2xb=0._dbl
          v1ca=0._dbl
          v1cb=0._dbl
          v2ca=0._dbl
          v2cb=0._dbl
          v2xab=0._dbl
          v2cab=0._dbl
          vttxa=0._dbl
          vttxb=0._dbl
          vttca=0._dbl
          vttcb=0._dbl
        ENDIF

      END SUBROUTINE tpss_p

! *****************************************************************************

      SUBROUTINE tpsspc(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub, &
                        sc,v1ca,v2ca,v1cb,v2cb,v2cab,vtta,vttb)

    REAL(dbl), INTENT(IN)                    :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, taua, taub
    REAL(dbl), INTENT(OUT)                   :: sc, v1ca, v2ca, v1cb, v2cb, &
                                                v2cab, vtta, vttb

    REAL(dbl), PARAMETER                     :: d = 2.80_dbl

    INTEGER                                  :: iflg
    REAL(dbl) :: dedga, dedgab, dedgb, dedra, dedrb, dedz, dzdg, dzdr, e, ec, &
      edaz, edbz, eta, grho, op, rdzdt, rho, ro, rs, tau, to, vca, vcb, z

        rho = rhoa + rhob
        grho = grhoaa + 2._dbl*grhoab + grhobb
        tau = taua + taub
        IF( ABS(tau) > small .AND. grho > small) THEN
          z    = 0.125_dbl*grho/rho/tau
          dzdr = -z/rho
          rdzdt = -z*rho/tau
          dzdg = z/grho
          CALL revpkzbp(rhoa,rhob,grhoaa,grhobb,grhoab,z,e, &
                         dedra,dedrb,dedga,dedgb,dedgab,dedz)
          op   = 1._dbl + d * e * z**3
          sc   = rho * e * op
          edaz = dedra + dedz*dzdr
          v1ca = e*op + rho*edaz*op + rho*e*d*z*z*(edaz*z+3._dbl*e*dzdr)
          edbz = dedrb + dedz*dzdr
          v1cb = e*op + rho*edbz*op + rho*e*d*z*z*(edbz*z+3._dbl*e*dzdr)
          ro   = e*d*e*3._dbl*z*z
          to   = rho * ( 1._dbl + 2._dbl*d*e*z**3 )
          v2ca = 2._dbl*rho*ro*dzdg + 2._dbl*to*dedz*dzdg + to*dedga
          v2cb = 2._dbl*rho*ro*dzdg + 2._dbl*to*dedz*dzdg + to*dedgb
          v2cab= 2._dbl*rho*ro*dzdg + 2._dbl*to*dedz*dzdg + to*dedgab
          vtta = vtta + dedz*rdzdt*op + e*d*dedz*rdzdt*z**3 + ro*rdzdt
          vttb = vttb + dedz*rdzdt*op + e*d*dedz*rdzdt*z**3 + ro*rdzdt
      ELSE
        eta=(rhoa-rhob)/rho
        IF(ABS(eta) > 1._dbl) eta=SIGN(1._dbl,eta)
        rs = (0.75_dbl/(pi*rho))**f13
        iflg=2
        IF(rs < 1.0_dbl) iflg=1
        CALL pz_p(rs,eta,ec,vca,vcb,iflg)
        sc   = rho * ec
        v1ca = vca
        v2ca = 0._dbl
        v1cb = vcb
        v2cb = 0._dbl
        v2cab= 0._dbl
      ENDIF

      END SUBROUTINE tpsspc

! *****************************************************************************

      SUBROUTINE revpkzbp(rhoa,rhob,grhoaa,grhobb,grhoab,z,e, &
                         dedra,dedrb,dedga,dedgb,dedgab,dedz)

    REAL(dbl), INTENT(IN)                    :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, z
    REAL(dbl), INTENT(OUT)                   :: e, dedra, dedrb, dedga, &
                                                dedgb, dedgab, dedz

    REAL(dbl), PARAMETER                     :: C00 = 0.54_dbl

    REAL(dbl) :: cxe, dcdgaa, dcdgab, dcdgbb, dcdra, dcdrb, e2, eca, ecb, &
      ef1, ef2, rho, sc, sca, scb, t1a, t1b, t2a, t2aba, t2abb, t2b, tv1ca, &
      tv1cb, tv2ca, tv2cb, txaba, txabb, tyaba, tyabb, v1ca, v1cb, v2ca, &
      v2cab, v2cb

        CALL tpbecs(rhoa,rhob,grhoaa,grhoab,grhobb,sc,v1ca,v2ca,v1cb,v2cb,v2cab)
        CALL tpbeca(rhoa,grhoaa,sca,tv1ca,tv2ca)
        CALL tpbeca(rhob,grhobb,scb,tv1cb,tv2cb)

        IF(sca > sc) THEN
          eca = sca
          t1a = tv1ca
          IF(grhoaa > small) THEN
             t2a = tv2ca/SQRT(grhoaa)
          ELSE
             t2a = 0._dbl
          ENDIF
          t2aba = 0._dbl
          txaba = 0._dbl
          tyaba = 0._dbl
        ELSE
          eca = sc
          t1a = v1ca
          t2a = v2ca
          t2aba = v2cab
          txaba = v1cb
          tyaba = v2cB
        ENDIF
        IF(scb > sc) THEN
          ecb = scb
          t1b = tv1cb
          IF(grhobb > small) THEN
             t2b = tv2cb/SQRT(grhobb)
          ELSE
             t2b = 0._dbl
          ENDIF
          t2abb = 0._dbl
          txabb = 0._dbl
          tyabb = 0._dbl
        ELSE
          ecb = sc
          t1b = v1cb
          t2b = v2cb
          t2abb = v2cab
          txabb = v1ca
          tyabb = v2cA
        ENDIF
!
        CALL ccfun(rhoa,rhob,grhoaa,grhobb,grhoab,cxe,dcdra,dcdrb, &
                   dcdgaa,dcdgbb,dcdgab)
!
        rho    = rhoa + rhob
        ef1    = 1._dbl + cxe * z * z
        ef2    = ( 1._dbl + cxe ) * z * z
        e2     = rhoa/rho * eca + rhob/rho * ecb
        e      = sc * ef1 - ef2 * e2
        dedra  = v1ca * ef1 - (-e2+eca + rhoa*t1a + rhob*txabb)/rho * ef2 &
                 + dcdra*z*z * ( sc - e2 )
        dedrb  = v1cb * ef1 - (-e2+ecb + rhob*t1b + rhoa*txaba)/rho * ef2 &
                 + dcdrb*z*z * ( sc - e2 )
        dedga  = v2ca * ef1 - ef2 * ( rhoa/rho*t2a + rhob/rho*tyabb ) &
                 + 2._dbl*dcdgaa*z*z * ( sc - e2 )
        dedgb  = v2cb * ef1 - ef2 * ( rhob/rho*t2b + rhoa/rho*tyaba ) &
                 + 2._dbl*dcdgbb*z*z * ( sc - e2 )
        dedgab = v2cab * ef1 - ef2 * (rhoa/rho*t2aba + rhob/rho*t2abb) &
                 + 2._dbl*dcdgab*z*z * ( sc - e2 )
        dedz   = 2._dbl * sc * cxe * z - 2._dbl * ( 1._dbl + cxe ) * z * e2

      END SUBROUTINE revpkzbp

! *****************************************************************************

END MODULE xc_tpss_p

! *****************************************************************************
