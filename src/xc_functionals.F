!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 CP2K developers group                                  !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****s* cp2k/xc_functionals [1.0] *
!!
!!   NAME
!!     xc_functionals
!!
!!   FUNCTION
!!     Meta Module of Exchange and Correlation Functionals
!!
!!   AUTHOR
!!     JGH (29-May-2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!   SOURCE
!******************************************************************************

MODULE xc_functionals
!MK *** Just a link to the old routine to have PBE ***
!MK
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE functionals,                     ONLY: init_functionals,&
                                             pbe_c,&
                                             pbe_x
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_cs1,                          ONLY: cs1
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_exchange_gga,                 ONLY: b88,&
                                             optx,&
                                             pbex,&
                                             revpbex,&
                                             x_pw86,&
                                             x_pw91
  USE xc_functionals_utilities,        ONLY: setup_calculation
  USE xc_ke_gga,                       ONLY: ke_pw86,&
                                             ke_pw91,&
                                             lc,&
                                             llp,&
                                             ol1,&
                                             ol2,&
                                             t92
  USE xc_lyp,                          ONLY: lyp
  USE xc_pade,                         ONLY: pade
  USE xc_perdew86,                     ONLY: perdew86
  USE xc_perdew_wang,                  ONLY: perdew_wang
  USE xc_perdew_zunger,                ONLY: perdew_zunger
  USE xc_tfw,                          ONLY: tfw
  USE xc_thomas_fermi,                 ONLY: thomas_fermi
  USE xc_vwn,                          ONLY: vwn
  USE xc_xalpha,                       ONLY: xalpha
  LOGICAL, PARAMETER :: GRADIENT_FUNCTIONAL = .TRUE., CROSS_TERM = .TRUE.

  PRIVATE
  PUBLIC :: b88, x_pw86, x_pw91, pbex, revpbex, optx, &
            ol1, ol2, llp, ke_pw86, ke_pw91, lc, t92, &
            lyp, pade, xalpha, thomas_fermi, tfw, vwn, perdew_wang, &
            perdew_zunger, perdew86, cs1
  PUBLIC :: xc_get_functional, xc_print_info, xc_set, xc_lda, xc_lsd, &
            xc_get_pot_size, xc_create_2nd_deriv_info, xc_calculate_lsd

  REAL ( dbl ), PARAMETER :: default_cutoff = 1.0E-12_dbl
  REAL ( dbl ), PARAMETER :: default_g_cutoff = 1.0E-10_dbl
  REAL ( dbl ) :: cutoff, g_cutoff, xalpha_coefficient
  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='xc_functionals'

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE xc_set ( density_cutoff, gradient_cutoff, falpha )
     IMPLICIT NONE
     REAL ( dbl ), INTENT(IN), OPTIONAL :: density_cutoff
     REAL ( dbl ), INTENT(IN), OPTIONAL :: gradient_cutoff
     REAL ( dbl ), INTENT(IN), OPTIONAL :: falpha

     IF ( PRESENT ( density_cutoff ) ) THEN
        cutoff = density_cutoff
     ELSE
        cutoff = default_cutoff
     ENDIF
     IF ( PRESENT ( gradient_cutoff ) ) THEN
        g_cutoff = gradient_cutoff
     ELSE
        g_cutoff = default_g_cutoff
     ENDIF
     IF ( PRESENT ( falpha ) ) THEN
        xalpha_coefficient = falpha
     ELSE
        xalpha_coefficient = 2._dbl/3._dbl
     ENDIF

  END SUBROUTINE xc_set

!******************************************************************************

  FUNCTION xc_get_pot_size ( lsd, gc, ct, order ) RESULT ( pot_size )
     LOGICAL, INTENT ( IN ) :: lsd, gc, ct
     INTEGER, INTENT ( IN ) :: order

     INTEGER :: pot_size

     INTEGER, PARAMETER :: omax = 3
     INTEGER :: m(0:omax,1:100)
     LOGICAL :: calc(0:omax)
     INTEGER :: tag, o

     tag = 100
     IF ( lsd ) tag = tag + 100
     IF ( gc )  tag = tag + 10
     IF ( ct )  tag = tag + 1

     m = 0
     CALL setup_calculation(order,m,calc,tag)

     o = ABS ( order )
     pot_size = MAXVAL(m(0:o,2)) - MINVAL(m(0:o,1)) + 1

  END FUNCTION xc_get_pot_size

!******************************************************************************

   SUBROUTINE xc_get_order ( lsd, gc, ct, gc_all, ct_all, order, index )

     IMPLICIT NONE

     LOGICAL, INTENT ( IN ) :: lsd, gc, ct, gc_all, ct_all
     INTEGER, INTENT ( IN ) :: order
     INTEGER, DIMENSION (:), INTENT ( OUT ) :: index
    
     INTEGER :: i

     IF ( .NOT. LSD ) THEN
        IF ( gc_all .AND. .NOT. gc ) THEN
           INDEX ( 1 ) = 1
           IF ( order >= 1 ) INDEX ( 2 ) = 2
           IF ( order >= 2 ) INDEX ( 3 ) = 4
           IF ( order >= 3 ) INDEX ( 4 ) = 7
        ELSE
           DO i = 1, SIZE ( index )
              INDEX ( i ) = i
           END DO
        END IF
     ELSE
        IF ( order == 0 ) THEN
           INDEX ( 1 ) = 1
        ELSEIF ( order == 1 ) THEN
           CALL stop_program ( "xc_get_order", "not implemented" )
        ELSE
           CALL stop_program ( "xc_get_order", "not implemented" )
        END IF
     END IF

  END SUBROUTINE xc_get_order 

!******************************************************************************

  SUBROUTINE xc_lda ( derivative_set, functional, gradient_functionals, crossterms, &
                      rho, order, drho, error )

     IMPLICIT NONE

     TYPE(xc_derivative_set_type), POINTER          :: derivative_set
     CHARACTER (LEN=*), DIMENSION ( 3 ), INTENT(IN) :: functional
     LOGICAL, DIMENSION ( 3 ), INTENT(IN)           :: gradient_functionals, crossterms
     REAL(dbl), DIMENSION (:), INTENT(in)           :: rho
     REAL(dbl), DIMENSION (:), INTENT(in), OPTIONAL :: drho
     INTEGER, INTENT(IN)                            :: order
     TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error

     LOGICAL :: gradient_correction, cross_term, docalc ( 3 )
     REAL ( dbl ), DIMENSION (:,:), ALLOCATABLE :: locpot
     INTEGER :: ifun, i, ierr, lp, lpmax, ind ( 100 )
     INTEGER :: npotentials, ngridpoints

     CHARACTER(len=*), PARAMETER :: routineN = 'xc_lda',&
                                    routineP = moduleN//':'//routineN

     gradient_correction = ANY(gradient_functionals)
     ! no crossterms in unpolarised spin calculations
     cross_term = .FALSE.

     CALL dcopy(SIZE(derivative_set%data,1)*SIZE(derivative_set%data,2), &
                0.0_dbl, 0, derivative_set%data, 1)

     IF ( .NOT. gradient_correction ) THEN
        npotentials = xc_get_pot_size(.FALSE., gradient_correction, cross_term, order)
        IF (npotentials > SIZE ( derivative_set%data, 2 )) &
             CALL stop_program ( "xc_lda", "Size of array DERIVATIVE_SET%DATA too small" )
        CALL xc_calculate_lda ( derivative_set, functional=functional(1), rho=rho,&
             pot=derivative_set%data, order=order, error=error )
        CALL xc_calculate_lda ( derivative_set, functional=functional(2), rho=rho,&
             pot=derivative_set%data, order=order, error=error )
        CALL xc_calculate_lda ( derivative_set, functional=functional(3), rho=rho,&
             pot=derivative_set%data, order=order, error=error )

     ELSE
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        docalc ( : ) = ( INDEX ( functional(:), "NONE" ) == 0 )
        ngridpoints = SIZE ( derivative_set%data, 1 )
        npotentials = SIZE ( derivative_set%data, 2 )
        lpmax = xc_get_pot_size ( .FALSE., .TRUE., cross_term, order )
        IF ( lpmax > npotentials ) CALL stop_program ( "xc_lda", &
             "Size of array DERIVATIVE_SET%DATA too small" )
        ALLOCATE ( locpot ( ngridpoints, npotentials ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "xc_lda", "locpot", ngridpoints*npotentials )
        DO ifun = 1, 3
           IF ( docalc ( ifun ) ) THEN
              lp = xc_get_pot_size ( .FALSE., gradient_functionals(ifun), .FALSE., order )
              CALL dcopy(SIZE(locpot,1)*SIZE(locpot,2), 0.0_dbl, 0, locpot, 1)
              CALL xc_calculate_lda ( derivative_set, functional=functional(ifun),&
                   rho=rho, drho=drho, pot=locpot, order=order, error=error )
              CALL xc_get_order ( .FALSE., gradient_functionals(ifun), .FALSE., &
                   gradient_correction, .FALSE., order, ind(1:lp) )
              DO i = 1, lp
                 CALL daxpy ( ngridpoints, 1._dbl, locpot(1,i), 1, derivative_set%data(1,ind(i)-1), 1 )
              END DO
           END IF
        END DO
        DEALLOCATE ( locpot, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "xc_lda", "locpot" )
     END IF

  END SUBROUTINE xc_lda

!******************************************************************************

  SUBROUTINE xc_lsd ( derivative_set, functional, gradient_functionals, crossterms, &
                      rho_a, rho_b, order, &
                      drho_a, drho_b, error )

     IMPLICIT NONE

     ! arguments
     TYPE(xc_derivative_set_type), POINTER :: derivative_set
     CHARACTER (LEN=*), DIMENSION ( 3 ), INTENT(IN) :: functional
     LOGICAL, DIMENSION ( 3 ), INTENT(IN)           :: gradient_functionals, crossterms
     REAL ( dbl ), DIMENSION (:), INTENT(in)        :: rho_a, rho_b
     REAL ( dbl ), DIMENSION (:), INTENT(in), OPTIONAL :: drho_a, drho_b
     INTEGER, INTENT(IN) :: order
     TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

     ! locals
     LOGICAL :: gradient_correction, cross_term, docalc ( 3 )
     REAL ( dbl ), DIMENSION (:,:), ALLOCATABLE :: locpot
     INTEGER :: ifun, i, ierr, np, mp, lp, lpmax, ind ( 100 )

     CHARACTER(len=*), PARAMETER :: routineN='xc_lsd',&
                                    routineP=moduleN//':'//routineN

     ! ------------------------------------------------------------

     gradient_correction = ANY ( gradient_functionals )
     cross_term = ANY ( crossterms )

     CALL dcopy(SIZE(derivative_set%data,1)*SIZE(derivative_set%data,2), 0.0_dbl, 0, derivative_set%data, 1)

     IF ( .NOT. gradient_correction ) THEN
        
        IF (xc_get_pot_size(.TRUE., gradient_correction, cross_term, order) > SIZE(derivative_set%data,2)) &
             CALL stop_program ( routineP, "Size of array DERIVATIVE_SET%DATA too small" )
        
        CALL xc_calculate_lsd ( derivative_set, functional(1), rho_a, rho_b, derivative_set%data, order, error=error )
        CALL xc_calculate_lsd ( derivative_set, functional(2), rho_a, rho_b, derivative_set%data, order, error=error )
        CALL xc_calculate_lsd ( derivative_set, functional(3), rho_a, rho_b, derivative_set%data, order, error=error )

     ELSE
        CPPreconditionNoFail(PRESENT(drho_a),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drho_b),cp_failure_level,routineP,error)
        
        docalc ( : ) = ( INDEX ( functional(:), "NONE" ) == 0 )
        
        np = SIZE ( derivative_set%data, 1 )
        mp = SIZE ( derivative_set%data, 2 )
        
        lpmax = xc_get_pot_size ( .TRUE., gradient_correction, cross_term, order )
        
        IF ( lpmax > mp ) CALL stop_program ( routineP, "Size of array DERIVATIVE_SET%DATA too small" )
        
        ALLOCATE ( locpot ( np, mp ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( routineP, "locpot", np*mp )
        DO ifun = 1, 3
           IF ( docalc ( ifun ) ) THEN
              lp = xc_get_pot_size ( .TRUE., gradient_functionals(ifun), crossterms(ifun), order )
              locpot(:,1:lp) = 0._dbl
              CALL xc_calculate_lsd ( derivative_set, functional(ifun),&
                   rho_a, rho_b, locpot, order, drhoa=drho_a, drhob=drho_b, error=error )
              CALL xc_get_order ( .TRUE., gradient_functionals(ifun), crossterms(ifun), &
                   gradient_correction, cross_term, order, ind(1:lp) )
              DO i = 1, lp
                 CALL daxpy ( np, 1._dbl, locpot(1,i), 1, derivative_set%data(1,ind(i)-1), 1 )
              END DO
           END IF
        END DO
        DEALLOCATE ( locpot, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "xc_lda", "locpot" )
     END IF

  END SUBROUTINE xc_lsd

!******************************************************************************

  SUBROUTINE xc_calculate_lda ( derivative_set, functional, rho, pot, order, drho ,error )
     IMPLICIT NONE
     TYPE(xc_derivative_set_type), POINTER :: derivative_set
     CHARACTER (LEN=*), INTENT(IN) :: functional
     REAL ( dbl ), DIMENSION (:), INTENT(in) :: rho
     REAL ( dbl ), DIMENSION (:), INTENT(in), OPTIONAL :: drho
     REAL ( dbl ), DIMENSION (:,:), INTENT(INOUT) :: pot
     INTEGER, INTENT(IN) :: order
     TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

     INTEGER :: i,ierr
     CHARACTER(len=3),PARAMETER :: routineN = 'xc_calculate_lda',&
                                   routineP = moduleN//':'//routineN

     ierr = 0

     SELECT CASE ( functional )
     CASE DEFAULT
        ierr = 1
     CASE ( "NONE" )
        ! do nothing
!exchange functional
     CASE ( "SLATER" )
        CALL xalpha ( cutoff, xalpha_coefficient )
        CALL xalpha ( derivative_set, rho, pot, order )
     CASE ( "BECKE88" )
        CALL b88 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL b88 ( derivative_set, rho, drho, pot, order )
     CASE ( "X_PW86" )
        CALL x_pw86 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL x_pw86 ( derivative_set, rho, drho, pot, order )
     CASE ( "X_PW91" )
        CALL x_pw91 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL x_pw91 ( derivative_set, rho, drho, pot, order )
     CASE ( "PBEX" )
!MK     CALL pbex ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
!MK     CALL pbex ( derivative_set, rho, drho, pot, order )
!MK *** Just a link to the old routine to have PBE ***
        CALL init_functionals(rho_threshold=cutoff,&
                              drho_threshold=g_cutoff)
        DO i=1,SIZE(rho)
          CALL pbe_x(rho(i),drho(i),pot(i,1),pot(i,2),pot(i,3))
          pot(i,3) = pot(i,3)*drho(i)
        END DO
!MK
     CASE ( "REVPBEX" )
        CALL revpbex ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL revpbex ( derivative_set, rho, drho, pot, order )
     CASE ( "OPTX" )
        CALL optx ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL optx ( derivative_set, rho, drho, pot, order )
!correlation functional
     CASE ( "PZ" )
        CALL perdew_zunger ( "ORIG", cutoff )
        CALL perdew_zunger ( derivative_set, rho, pot, order )
     CASE ( "PW" )
        CALL perdew_wang ( "PWO", cutoff )
        CALL perdew_wang ( derivative_set, rho, pot, order )
     CASE ( "VWN" )
        CALL vwn ( cutoff )
        CALL vwn ( derivative_set, rho, pot, order )
     CASE ( "PADE" )
        CALL pade ( cutoff )
        CALL pade ( derivative_set, rho, pot, order )
     CASE ( "PBEC" )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
!MK *** Just a link to the old routine to have PBE ***
        CALL init_functionals(rho_threshold=cutoff,&
                              drho_threshold=g_cutoff)
        DO i=1,SIZE(rho)
          CALL pbe_c(rho(i),drho(i),pot(i,1),pot(i,2),pot(i,3))
          pot(i,3) = pot(i,3)*drho(i)
        END DO
!MK
     CASE ( "PERDEW86" )
        CALL perdew86 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL perdew86 ( derivative_set, rho, drho, pot, order )
     CASE ( "LYP" )
        CALL lyp ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL lyp ( derivative_set, rho, drho, pot, order )
     CASE ( "CS1" )
        CALL cs1 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL cs1 ( derivative_set, rho, drho, pot, order )
!kinetic energy functional
     CASE ( "TF" )
        CALL thomas_fermi ( cutoff )
        CALL thomas_fermi ( derivative_set, rho, pot, order )
     CASE ( "TFW" )
        CALL tfw ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL tfw ( derivative_set, rho, drho, pot, order )
     CASE ( "OL1" )
        CALL ol1 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL ol1 ( derivative_set, rho, drho, pot, order )
     CASE ( "OL2" )
        CALL ol2 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL ol2 ( derivative_set, rho, drho, pot, order )
     CASE ( "LLP" )
        CALL llp ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL llp ( derivative_set, rho, drho, pot, order )
     CASE ( "KE_PW86" )
        CALL ke_pw86 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL ke_pw86 ( derivative_set, rho, drho, pot, order )
     CASE ( "KE_PW91" )
        CALL ke_pw91 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL ke_pw91 ( derivative_set, rho, drho, pot, order )
     CASE ( "LC" )
        CALL lc ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL lc ( derivative_set, rho, drho, pot, order )
     CASE ( "T92" )
        CALL t92 ( cutoff )
        CPPreconditionNoFail(PRESENT(drho),cp_failure_level,routineP,error)
        CALL t92 ( derivative_set, rho, drho, pot, order )
     END SELECT

     IF ( ierr /= 0 ) CALL stop_program ( "xc_calculate_lda", &
         " Unknown functional: "//functional )

  END SUBROUTINE xc_calculate_lda

!******************************************************************************

  SUBROUTINE xc_calculate_lsd ( derivative_set, functional, rhoa, rhob, pot, &
       order, drhoa, drhob, error )

     IMPLICIT NONE

     TYPE(xc_derivative_set_type), POINTER :: derivative_set
     CHARACTER (LEN=*), INTENT(IN) :: functional
     REAL ( dbl ), DIMENSION (:), INTENT(in) :: rhoa, rhob
     REAL ( dbl ), DIMENSION (:), INTENT(in), OPTIONAL :: drhoa, drhob
     REAL ( dbl ), DIMENSION (:,:), INTENT(INOUT) :: pot
     INTEGER, INTENT(IN) :: order
     TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

     INTEGER :: ierr
     CHARACTER(len=3),PARAMETER :: routineN = 'xc_calculate_lsd',&
                                   routineP = moduleN//':'//routineN

     ierr = 0

     SELECT CASE ( functional )
     CASE DEFAULT
        ierr = 1
     CASE ( "NONE" )
        ! do nothing
!exchange functional
     CASE ( "SLATER" )
        CALL xalpha ( cutoff, xalpha_coefficient )
        CALL xalpha ( derivative_set, rhoa, rhob, pot, order )
     CASE ( "BECKE88" )
        CALL b88 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL b88 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "X_PW86" )
        CALL x_pw86 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL x_pw86 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "X_PW91" )
        CALL x_pw91 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL x_pw91 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "PBEX" )
        CALL pbex ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL pbex ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "REVPBEX" )
        CALL revpbex ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL revpbex ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "OPTX" )
        CALL optx ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL optx ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
!correlation functional
     CASE ( "PZ" )
        CALL perdew_zunger ( "DMC", cutoff )
        CALL perdew_zunger ( derivative_set, rhoa, rhob, pot, order )
     CASE ( "PW" )
        CALL perdew_wang ( "PWO", cutoff )
        CALL perdew_wang ( derivative_set, rhoa, rhob, pot, order )
     CASE ( "VWN" )
        CALL vwn ( cutoff )
        CALL vwn ( derivative_set, rhoa, rhob, pot, order )
     CASE ( "PADE" )
        CALL pade ( cutoff )
        CALL pade ( derivative_set, rhoa, rhob, pot, order )
     CASE ( "PERDEW86" )
        CALL perdew86 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL perdew86 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "LYP" )
        CALL lyp ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL lyp ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
!kinetic energy functional
     CASE ( "TF" )
        CALL thomas_fermi ( cutoff )
        CALL thomas_fermi ( derivative_set, rhoa, rhob, pot, order )
     CASE ( "TFW" )
        CALL tfw ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL tfw ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "OL1" )
        CALL ol1 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL ol1 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "OL2" )
        CALL ol2 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL ol2 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "LLP" )
        CALL llp ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL llp ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "KE_PW86" )
        CALL ke_pw86 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL ke_pw86 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "KE_PW91" )
        CALL ke_pw91 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL ke_pw91 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "LC" )
        CALL lc ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL lc ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     CASE ( "T92" )
        CALL t92 ( cutoff )
        CPPreconditionNoFail(PRESENT(drhoa),cp_failure_level,routineP,error)
        CPPreconditionNoFail(PRESENT(drhob),cp_failure_level,routineP,error)
        CALL t92 ( derivative_set, rhoa, rhob, drhoa, drhob, pot, order )
     END SELECT

     IF ( ierr /= 0 ) CALL stop_program ( "xc_calculate_lsd", &
         " Unknown functional: "//functional )

  END SUBROUTINE xc_calculate_lsd

!******************************************************************************

  SUBROUTINE xc_get_functional ( tag, string, functional, gradient_functionals, crossterms )

     IMPLICIT NONE
     INTEGER, INTENT ( in ) :: tag
     CHARACTER (LEN=*), INTENT ( IN ) :: string
     CHARACTER (LEN=*), DIMENSION ( 3 ), INTENT ( INOUT ) :: functional
     LOGICAL, DIMENSION ( 3 ), INTENT (OUT) :: gradient_functionals, crossterms

     INTEGER :: ierr

     ierr = 0
     SELECT CASE ( tag )

     CASE DEFAULT
       CALL stop_program ( 'xc_get_functional', 'illegal tag' )

     CASE ( 0 )
     ! combined functionals
        gradient_functionals(1:2) = .FALSE.
        crossterms(1:2) = .FALSE.

        SELECT CASE ( string )
        CASE DEFAULT
          ierr = 1
        CASE ( "NONE" )
           functional ( 1 ) = "NONE"
           functional ( 2 ) = "NONE"
        CASE ( "XALPHA", "SLATER" )
           functional ( 1 ) = "SLATER"
           functional ( 2 ) = "NONE"
        CASE ( "LDA", "PZ" )
           functional ( 1 ) = "SLATER"
           functional ( 2 ) = "PZ"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "PW", "PERDEW-WANG" )
           functional ( 1 ) = "SLATER"
           functional ( 2 ) = "PW"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "VWN" )
           functional ( 1 ) = "SLATER"
           functional ( 2 ) = "VWN"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "PADE" )
           functional ( 1 ) = "NONE"
           functional ( 2 ) = "PADE"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "BP", "BP86", "BECKE-PERDEW" )
           functional ( 1 ) = "BECKE88"
           functional ( 2 ) = "PERDEW86"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "BLYP", "BECKE-LEE-YANG-PARR" )
           functional ( 1 ) = "BECKE88"
           functional ( 2 ) = "LYP"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "PBE" )
           functional ( 1 ) = "PBEX"
           functional ( 2 ) = "PBEC"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "OLYP" )
           functional ( 1 ) = "OPTX"
           functional ( 2 ) = "LYP"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        END SELECT

     CASE ( 1 )
     ! exchange functionals
        gradient_functionals(1) = .FALSE.
        crossterms(1) = .FALSE.

        SELECT CASE ( string )
        CASE DEFAULT
          ierr = 1
        CASE ( "NONE" )
           functional ( 1 ) = "NONE"
        CASE ( "XALPHA", "SLATER" )
           functional ( 1 ) = "SLATER"
        CASE ( "BECKE", "BECKE88" )
           functional ( 1 ) = "BECKE88"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        CASE ( "PW86", "PERDEW-WANG86" )
           functional ( 1 ) = "X_PW86"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        CASE ( "PW91", "PERDEW-WANG91" )
           functional ( 1 ) = "X_PW91"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        CASE ( "PBE" )
           functional ( 1 ) = "PBEX"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        CASE ( "REVPBE" )
           functional ( 1 ) = "REVPBEX"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        CASE ( "OPTX" )
           functional ( 1 ) = "OPTX"
           gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        END SELECT

     CASE ( 2 )
     ! correlation functionals
        gradient_functionals(2) = .FALSE.
        crossterms(2) = .FALSE.

        SELECT CASE ( string )
        CASE DEFAULT
          ierr = 1
        CASE ( "NONE" )
           functional ( 2 ) = "NONE"
        CASE ( "LDA", "PZ" )
           functional ( 2 ) = "PZ"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "PW", "PERDEW-WANG" )
           functional ( 2 ) = "PW"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "VWN" )
           functional ( 2 ) = "VWN"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "PADE" )
           functional ( 1 ) = "NONE"
           functional ( 2 ) = "PADE"
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "PBE" )
           functional ( 2 ) = "PBEC"
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "CS1" )
           functional ( 2 ) = "CS1"
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "LYP", "LEE-YANG-PARR" )
           functional ( 2 ) = "LYP"
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        CASE ( "P86", "PERDEW86" )
           functional ( 2 ) = "PERDEW86"
           gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
           crossterms ( 2 ) = CROSS_TERM
        END SELECT

     CASE ( 3 )
     ! kinetic energy functionals
        gradient_functionals(3) = .FALSE.
        crossterms(3) = .FALSE.

        SELECT CASE ( string )
        CASE DEFAULT
          ierr = 1
        CASE ( "NONE" )
           functional ( 3 ) = "NONE"
        CASE ( "TF", "THOMAS-FERMI" )
           functional ( 3 ) = "TF"
        CASE ( "TFW", "THOMAS-FERMI-WEIZSAECKER", "WEIZSAECKER" )
           functional ( 3 ) = "TFW"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "PW86", "PERDEW-WANG86" )
           functional ( 3 ) = "KE_PW86"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "PW91", "PERDEW-WANG91" )
           functional ( 3 ) = "KE_PW91"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "OL1", "OU-YANG-LEVY-1" )
           functional ( 3 ) = "OL1"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "OL2", "OU-YANG-LEVY-2" )
           functional ( 3 ) = "OL2"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "LLP", "LEE-LEE-PARR" )
           functional ( 3 ) = "LLP"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "LC", "LEMBARKI-CHERMETTE" )
           functional ( 3 ) = "LC"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        CASE ( "T92", "THAKKAR" )
           functional ( 3 ) = "T92"
           gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        END SELECT

     END SELECT

     IF ( ierr /= 0 ) THEN
       CALL stop_program ( 'xc_get_functional', 'Unknown Functional :'//string )
     END IF

  END SUBROUTINE xc_get_functional

!******************************************************************************

  SUBROUTINE xc_print_info ( iunit, functional, lsd )

    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: iunit
    CHARACTER ( LEN=* ), DIMENSION ( 3 ), INTENT ( IN ) :: functional
    LOGICAL, INTENT ( IN ) :: lsd

    CHARACTER (LEN=80), DIMENSION ( 3 ) :: reference
    CHARACTER (LEN=80), DIMENSION ( 3 ) :: shortform
    CHARACTER (LEN=80) :: not_tested
    LOGICAL, DIMENSION ( 3 ) :: tested

    tested = .FALSE.
    not_tested = " *** WARNING : THIS FUNCTIONAL IS NOT FULLY TESTED ***"

    SELECT CASE (functional(1))
    CASE ( "NONE" )
      reference(1)="NONE" 
      shortform(1)="NONE"
    CASE ( "SLATER" )
      CALL xalpha ( cutoff, xalpha_coefficient )
      CALL xalpha ( reference(1), shortform(1) )
      IF ( .NOT. LSD ) tested(1) = .TRUE.
    CASE ( "BECKE88" )
      CALL b88 ( reference(1), shortform(1) )
      IF ( .NOT. LSD ) tested(1) = .TRUE.
    CASE ( "X_PW86" )
      CALL x_pw86 ( reference(1), shortform(1) )
    CASE ( "X_PW91" )
      CALL x_pw91 ( reference(1), shortform(1) )
    CASE ( "PBEX" )
      CALL pbex ( reference(1), shortform(1) )
      IF ( .NOT. LSD ) tested(1) = .TRUE.
    CASE ( "REVPBEX" )
      CALL revpbex ( reference(1), shortform(1) )
      IF ( .NOT. LSD ) tested(1) = .TRUE.
    CASE ( "OPTX" )
      CALL optx ( reference(1), shortform(1) )
      IF ( .NOT. LSD ) tested(1) = .TRUE.
    CASE DEFAULT
      reference(1)="UNKNOWN" 
      shortform(1)="UNKNOWN"
    END SELECT

    SELECT CASE (functional(2))
    CASE ( "NONE" )
      reference(2)="NONE" 
      shortform(2)="NONE"
    CASE ( "PADE" )
      CALL pade ( reference(2), shortform(2) )
      IF ( .NOT. LSD ) tested(2) = .TRUE.
!MK
    CASE ( "PBEC" )
      CALL pbex ( reference(2), shortform(2) )
      IF ( .NOT. LSD ) tested(2) = .TRUE.
!MK
    CASE ( "VWN" )
      CALL vwn ( reference(2), shortform(2) )
      IF ( .NOT. LSD ) tested(2) = .TRUE.
    CASE ( "PW" )
      CALL perdew_wang ( reference(2), shortform(2) )
    CASE ( "PZ" )
      CALL perdew_zunger ( reference(2), shortform(2) )
    CASE ( "CS1" )
      CALL cs1 ( reference(2), shortform(2) )
    CASE ( "LYP" )
      CALL lyp ( reference(2), shortform(2) )
      IF ( .NOT. LSD ) tested(2) = .TRUE.
    CASE ( "PERDEW86" )
      CALL perdew86 ( reference(2), shortform(2) )
      IF ( .NOT. LSD ) tested(2) = .TRUE.
    CASE DEFAULT
      reference(2)="UNKNOWN" 
      shortform(2)="UNKNOWN"
    END SELECT

    SELECT CASE (functional(3))
    CASE ( "NONE" )
      reference(3)="NONE" 
      shortform(3)="NONE"
    CASE ( "TF" )
      CALL thomas_fermi ( reference(3), shortform(3) )
    CASE ( "TFW" )
      CALL tfw ( reference(3), shortform(3) )
    CASE ( "OL1" )
      CALL ol1 ( reference(3), shortform(3) )
    CASE ( "OL2" )
      CALL ol2 ( reference(3), shortform(3) )
    CASE ( "LLP" )
      CALL llp ( reference(3), shortform(3) )
    CASE ( "KE_PW86" )
      CALL ke_pw86 ( reference(3), shortform(3) )
    CASE ( "KE_PW91" )
      CALL ke_pw91 ( reference(3), shortform(3) )
    CASE ( "LC" )
      CALL lc ( reference(3), shortform(3) )
    CASE ( "T92" )
      CALL t92 ( reference(3), shortform(3) )
    CASE DEFAULT
      reference(3)="UNKNOWN" 
      shortform(3)="UNKNOWN"
    END SELECT

    WRITE ( iunit, '( A, T41, A )' ) " FUNCTIONAL| Exchange Functional", &
          ADJUSTR ( functional ( 1 ) (1:40) )
    IF ( functional ( 1 ) (1:4) /= "NONE" ) THEN
       WRITE ( iunit, '( A80 )' ) ADJUSTR ( reference ( 1 ) )
       IF ( .NOT. tested(1) ) WRITE ( iunit, '( A80 )' ) ADJUSTR ( not_tested )
    END IF
    WRITE ( iunit, '( A, T41, A )' ) " FUNCTIONAL| Correlation Functional", &
          ADJUSTR ( functional ( 2 ) (1:40) )
    IF ( functional ( 2 ) (1:4) /= "NONE" ) THEN
       WRITE ( iunit, '( A80 )' ) ADJUSTR ( reference ( 2 ) )
       IF ( .NOT. tested(2) ) WRITE ( iunit, '( A80 )' ) ADJUSTR ( not_tested )
    END IF
    IF ( functional ( 3 ) (1:4) /= "NONE" ) THEN
       WRITE ( iunit, '( A, T41, A )' ) " FUNCTIONAL| Kinetic Energy Functional", &
          ADJUSTR ( functional ( 3 ) (1:40) )
       WRITE ( iunit, '( A80 )' ) ADJUSTR ( reference ( 3 ) )
       IF ( .NOT. tested(3) ) WRITE ( iunit, '( A80 )' ) ADJUSTR ( not_tested )
    END IF

  END SUBROUTINE xc_print_info

!******************************************************************************

!!****f* xc_functionals_utilities/xc_create_2nd_deriv_info [1.0] *
!!
!!   NAME
!!     xc_create_2nd_deriv_info
!!
!!   SYNOPSIS
!!     Subroutine xc_create_2nd_deriv_info(lsd, gradient, crossterms,&
!!         spin_pot, grad_pot, ndiag_term, error)
!!       Logical, Intent (IN):: lsd, gradient, crossterms
!!       Integer, Dimension(:,:), Pointer:: spin_pot
!!       Logical, Dimension(:,:), Pointer:: grad_pot
!!       Logical, Dimension(:), Pointer:: ndiag_term
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_create_2nd_deriv_info
!!
!!   FUNCTION
!!     allocates and fills out spins_pot cross_pot and grad_pot
!!     with informations about the second derivatives.
!!
!!   NOTES
!!     The derivatives are always in the following order (but depending 
!!     on the functional some may be missing):
!!     aa,bb,daa,dbb,dada,dbdb,ab,dab,dba,dadb,dca,dcb,dcda,dcdb,dcdc
!!     ( a: deriv. wrt. alpha (spin 1), b: deriv. wrt. beta (spin 2),
!!     da: deriv. wrt. ||grad spin 1||, db deriv. wrt. ||grad spin 2||,
!!     dc: deriv. wrt. || grad(spin1+spin2)||).
!!     You are responsible of deallocating the arrays.
!!
!!   INPUTS
!!     - lsd: if the functional is using lsd (two spins)
!!     - gradient: if the functiuonal is a gradient corrected functional
!!     - crossterms: if the functional has crossterms
!!     - spins_pot(:,n) says with respect to which spins was derived the 
!!       n-th potential. In LSD spin 0 means alpha+beta (spin 1+2) and
!!       makes sense only for the gradient (norm of the gradient of alpha+
!!       beta)
!!     - grad_pot(s,n) says if the derivative of the n-th potential
!!       with respect to the spin spins_pot(s,n) was respect to the density
!!       (false) or its gradient (true)
!!     - ndiag_term(n) is true if the n-th potential derivative isn't a 
!!       diagonal term
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_create_2nd_deriv_info(lsd,gradient,crossterms,&
       spin_pot,grad_pot,ndiag_term,error)
    LOGICAL, INTENT(in) :: lsd, gradient,crossterms
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_create_2nd_deriv_info',&
         routineP=moduleN//':'//routineN
    INTEGER :: npot, stat
    INTEGER, DIMENSION(2,15) :: spin_template
    failure=.FALSE.

    npot=xc_get_pot_size(lsd,gradient,crossterms,-2)

    IF (.NOT.failure) THEN
       ALLOCATE(spin_pot(2,npot), grad_pot(2,npot), ndiag_term(npot), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CPAssert(npot<=15,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       spin_pot=1
       grad_pot=.FALSE.
       ndiag_term=.FALSE.
       IF (lsd) THEN
          spin_template=RESHAPE((/1,1,2,2,1,1,2,2,1,1,2,2,1,2,1,2,2,1,1,2,&
               0,1,0,2,0,1,0,2,0,0/),(/2,15/))
          spin_pot=spin_template(:,1:npot)
          IF (.NOT.gradient .AND. crossterms) spin_pot(:,1:3)=RESHAPE((/1,1,1,2,2,2/),(/2,3/))
          
          IF (gradient) THEN
             grad_pot(1,3:)=.TRUE.
             grad_pot(2,5:6)=.TRUE.
             IF (crossterms) THEN
                grad_pot(1,7)=.FALSE.
                grad_pot(2,10)=.TRUE.
                ndiag_term(7:)=.TRUE.
                ndiag_term(15)=.FALSE.
             END IF
          ELSE
             IF (crossterms) THEN
                ndiag_term(2)=.TRUE.
             END IF
          END IF
       ELSE
          ! lda
          IF (gradient) THEN
             grad_pot(1,2)=.TRUE.
             grad_pot(:,3)=.TRUE.
          END IF
       END IF
    END IF
  END SUBROUTINE xc_create_2nd_deriv_info
!***************************************************************************

END MODULE xc_functionals

!******************************************************************************
