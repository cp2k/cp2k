!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_functionals_tests [1.0] *
!!
!!   NAME
!!     xc_functionals_tests
!!
!!   FUNCTION
!!     Test routines for XC functionals
!!
!!   AUTHOR
!!     JGH (15.02.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     List of functionals
!!      1    X-Alpha
!!      2    Pade-XC-LDA
!!      3    VWN-Correlation LDA
!!      4    Thomas-Fermi
!!      5    Perdew-Wang Correlation LDA 
!!      6    Perdew-Zunger Correlation LDA 
!!      7    Thomas-Fermi-von Weizsaecker 
!!      8    OL1 Kinetic energy functional
!!      9    OL2 Kinetic energy functional
!!     10    LLP Kinetic energy functional
!!     11    PW86 Kinetic energy functional
!!     12    PW91 Kinetic energy functional
!!     13    LC Kinetic energy functional
!!     14    T92 Kinetic energy functional
!!     15    Becke 1988 Exchange
!!     16    Perdew-Wang 86 Exchange
!!     17    Perdew-Wang 91 Exchange
!!     18    Perdew-Burke-Ernzerhof Correlation
!!     19    PBE Correlation (Zhang Parameters)
!!     20    OPTX Handy Exchange
!!     21    Lee-Yang-Parr
!!     22    Perdew 1986 Correlation
!!     23    CS1 Correlation
!!
!!*****
!******************************************************************************

MODULE xc_functionals_tests
! *****************************************************************************
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals,                  ONLY: b88,&
                                             cs1,&
                                             ke_pw86,&
                                             ke_pw91,&
                                             lc,&
                                             llp,&
                                             lyp,&
                                             ol1,&
                                             ol2,&
                                             optx,&
                                             pade,&
                                             pbex,&
                                             perdew86,&
                                             perdew_wang,&
                                             perdew_zunger,&
                                             revpbex,&
                                             t92,&
                                             tfw,&
                                             thomas_fermi,&
                                             vwn,&
                                             x_pw86,&
                                             x_pw91,&
                                             xalpha

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xc_functional_test

  INTEGER, PARAMETER :: n=20
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rhoref, rhorefa, rhorefb
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rho, rhoa, rhob, grho, grhoref
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: grhoa, grhob, grhorefa, grhorefb, &
                                          grhorefab, grhoab
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: pot, potp, potm, potx
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgp, potgm
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgpa, potgpb, potgma, potgmb
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgpab, potgmab
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potpa, potpb, potma, potmb
  REAL(dbl), PARAMETER :: rrho=1.E-5_dbl, epsilon=1.E-15_dbl
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: relerr, abserr

  INTEGER :: iw

  INTEGER, PARAMETER :: nfun = 23  ! number of implemented functionals

  INTERFACE funlda
    MODULE PROCEDURE fun_lda, fun_lsd, fun_init, fun_info
  END INTERFACE

  INTERFACE fungga
    MODULE PROCEDURE funggau, funggap, funggab, fun_init, fun_info
  END INTERFACE

! *****************************************************************************

CONTAINS

!******************************************************************************
!!****** xc_functionals_tests/xc_functional_test [1.0] *
!!
!!   NAME
!!     xc_functional_test
!!
!!   SYNOPSIS
!!     Subroutine xc_functional_test(globenv)
!!       Type(global_environment_type), Intent (IN):: globenv
!!     End Subroutine xc_functional_test
!!
!!   FUNCTION
!!     Tests XC-functional derivaties by finite derivatives
!!
!!   AUTHOR
!!     JGH  15-Feb-2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****
!******************************************************************************

SUBROUTINE xc_functional_test ( globenv )

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=80)                        :: fun_string
    INTEGER                                  :: funtype, i
    LOGICAL                                  :: crossterm
    REAL(dbl)                                :: r
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set

! *****************************************************************************

  iw =globenv % scr

  IF ( globenv % ionode ) THEN

     CALL allocate_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test LDA functionals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF LDA FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     !
     DO i=1,n
       r = 0.01_dbl + (i-1)*25._dbl/n
       rhoref(i)=3.e3 * EXP(-r)
     END DO
     !
     DO i = 1, nfun
       CALL funlda ( i, epsilon, funtype, crossterm )
       IF ( funtype /= 1 .AND. funtype /= 3 ) CYCLE
       CALL funlda ( i, reference=fun_string )
       WRITE (iw,'(A80)') ADJUSTR(fun_string)
       CALL funlda ( derivative_set, i, rhoref, pot, 3 )
       rho = rhoref + rrho*rhoref
       CALL funlda ( derivative_set, i, rho, potp, 3 )
       rho = rhoref - rrho*rhoref
       CALL funlda ( derivative_set, i, rho, potm, 3 )
       !
       CALL testder (" First Derivative ",pot(:,1),potp(:,0),potm(:,0),rhoref)
       CALL testder (" Second Derivative ",pot(:,2),potp(:,1),potm(:,1),rhoref)
       CALL testder (" Third Derivative ",pot(:,3),potp(:,2),potm(:,2),rhoref)
     END DO
     !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test GGA-LDA functionals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF GGA-LDA FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     !
     DO i=1,n
       r = 0.01_dbl + (i-1)*25._dbl/n
       rhoref(i)=3.e3 * EXP(-r)
       grhoref(i)=7.e3 * EXP(-r)
     END DO
     !
     DO i = 1, nfun
       CALL fungga ( i, epsilon, funtype, crossterm )
       IF ( funtype /= 11 .AND. funtype /= 13 ) CYCLE
       CALL fungga ( i, reference=fun_string )
       WRITE (iw,'(A80)') ADJUSTR(fun_string)
       CALL fungga ( derivative_set, i, rhoref, grhoref, pot, 3 )
       rho = rhoref + rrho*rhoref
       CALL fungga ( derivative_set, i, rho, grhoref, potp, 3 )
       rho = rhoref - rrho*rhoref
       CALL fungga ( derivative_set, i, rho, grhoref, potm, 3 )
       grho = grhoref + rrho*grhoref
       CALL fungga ( derivative_set, i, rhoref, grho, potgp, 3 )
       grho = grhoref - rrho*grhoref
       CALL fungga ( derivative_set, i, rhoref, grho, potgm, 3 )
       !
       CALL testder (" First Derivative R",pot(:,1),potp(:,0),potm(:,0),rhoref)
       CALL testder (" First Derivative G",pot(:,2),potgp(:,0),potgm(:,0),grhoref)
       CALL testder (" Second Derivative RR",pot(:,3),potp(:,1),potm(:,1),rhoref)
       CALL testder (" Second Derivative RG",pot(:,4),potgp(:,1),potgm(:,1),grhoref)
       CALL testder (" Second Derivative GG",pot(:,5),potgp(:,2),potgm(:,2),grhoref)
       CALL testder (" Third Derivative RRR",pot(:,6),potp(:,3),potm(:,3),rhoref)
       CALL testder (" Third Derivative RRG",pot(:,7),potgp(:,3),potgm(:,3),grhoref)
       CALL testder (" Third Derivative RGG",pot(:,8),potgp(:,4),potgm(:,4),grhoref)
       CALL testder (" Third Derivative GGG",pot(:,9),potgp(:,5),potgm(:,5),grhoref)
       !
     END DO
     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test LSD functionals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF LSD FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     !
     DO i=1,n
       r = 0.01_dbl + (i-1)*25._dbl/n
       rhorefa(i)=5.e3 * EXP(-r)
       rhorefb(i)=3.e3 * EXP(-r)*(1._dbl-0.1_dbl*SIN(r))
     END DO
     rhoref = 2._dbl*rhorefa
     !
     DO i = 1, nfun
       CALL funlda ( i, epsilon, funtype, crossterm )
       IF ( funtype /= 2 .AND. funtype /= 3 ) CYCLE
       CALL funlda ( i, reference=fun_string )
       WRITE (iw,'(A80)') ADJUSTR(fun_string)
       CALL funlda ( derivative_set, i, rhoref, potx(:,1:1), 0 )
       CALL funlda ( derivative_set, i, rhorefa, rhorefa, potx(:,2:2), 0 )
       CALL funlda ( derivative_set, i, rhorefa, rhorefb, pot, 3 )
       rhoa = rhorefa + rrho*rhorefa
       rhob = rhorefb
       CALL funlda ( derivative_set, i, rhoa, rhob, potpa, 3 )
       rhoa = rhorefa - rrho*rhorefa
       CALL funlda ( derivative_set, i, rhoa, rhob, potma, 3 )
       rhoa = rhorefa
       rhob = rhorefb + rrho*rhorefb
       CALL funlda ( derivative_set, i, rhoa, rhob, potpb, 3 )
       rhob = rhorefb - rrho*rhorefb
       CALL funlda ( derivative_set, i, rhoa, rhob, potmb, 3 )
       !
       IF ( crossterm ) THEN
          abserr(:) = ABS(potx(:,1)-potx(:,2))
          relerr(:) = abserr(:)/ABS(potx(:,1))
          CALL output ( iw, rhoref, potx(:,1), abserr, relerr, &
                      " Compare LDA/LSD Functionals " )
          !
          CALL testder (" First Derivative A",pot(:,1),potpa(:,0),potma(:,0),rhorefa)
          CALL testder (" First Derivative B",pot(:,2),potpb(:,0),potmb(:,0),rhorefb)
          CALL testder (" Second Derivative AA",pot(:,3),potpa(:,1),potma(:,1),rhorefa)
          CALL testder (" Second Derivative AB",pot(:,4),potpa(:,2),potma(:,2),rhorefa)
          CALL testder (" Second Derivative BB",pot(:,5),potpb(:,2),potmb(:,2),rhorefb)
          CALL testder (" Third Derivative AAA",pot(:,6),potpa(:,3),potma(:,3),rhorefa)
          CALL testder (" Third Derivative AAB",pot(:,7),potpb(:,3),potmb(:,3),rhorefb)
          CALL testder (" Third Derivative ABB",pot(:,8),potpa(:,5),potma(:,5),rhorefa)
          CALL testder (" Third Derivative BBB",pot(:,9),potpb(:,5),potmb(:,5),rhorefb)
       ELSE
          abserr(:) = ABS(potx(:,1)-potx(:,2))
          relerr(:) = abserr(:)/ABS(potx(:,1))
          CALL output ( iw, rhoref, potx(:,1), abserr, relerr, &
                      " Compare LDA/LSD Functionals " )
          !
          CALL testder (" First Derivative A",pot(:,1),potpa(:,0),potma(:,0),rhorefa)
          CALL testder (" First Derivative B",pot(:,2),potpb(:,0),potmb(:,0),rhorefb)
          CALL testder (" Second Derivative AA",pot(:,3),potpa(:,1),potma(:,1),rhorefa)
          CALL testder (" Second Derivative BB",pot(:,4),potpb(:,2),potmb(:,2),rhorefb)
          CALL testder (" Third Derivative AAA",pot(:,5),potpa(:,3),potma(:,3),rhorefa)
          CALL testder (" Third Derivative BBB",pot(:,6),potpb(:,4),potmb(:,4),rhorefb)
       END IF
     END DO
     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test GGA-LSD functionals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF GGA-LSD FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     !
     DO i=1,n
       r = 0.1_dbl + (i-1)*25._dbl/n
       rhorefa(i)=3.e3 * EXP(-r)
       grhorefa(i)=7.e3 * EXP(-r)
       rhorefb(i)=3.e3 * EXP(-r)*(1._dbl-0.1*SIN(r))
       grhorefb(i)=7.e3 * EXP(-r)*(1._dbl-0.1*SIN(r))
       grhorefab(i)=(grhorefa(i)+grhorefb(i))*0.5_dbl
     END DO
     rhoref = 2._dbl*rhorefa
     grhoref= 2._dbl*grhorefa
     !
     DO i = 1, nfun
       CALL fungga ( i, epsilon, funtype, crossterm )
       IF ( funtype /= 12 .AND. funtype /= 13 ) CYCLE
       CALL fungga ( i, reference=fun_string )
       WRITE (iw,'(A80)') ADJUSTR(fun_string)
       IF ( crossterm ) THEN
          CALL fungga ( derivative_set, i, rhoref, grhoref, potx(:,1:1), 0 )
          CALL fungga ( derivative_set, i, rhorefa, rhorefa, grhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhorefb, grhorefab, pot, 3 )
          rhoa = rhorefa + rrho*rhorefa
          rhob = rhorefb
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, grhorefab, potpa, 3 )
          rhoa = rhorefa - rrho*rhorefa
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, grhorefab, potma, 3 )
          rhoa = rhorefa
          rhob = rhorefb + rrho*rhorefb
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, grhorefab, potpb, 3 )
          rhob = rhorefb - rrho*rhorefb
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, grhorefab, potmb, 3 )
          grhoa = grhorefa + rrho*grhorefa
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhoa, grhorefb, grhorefab, potgpa, 3 )
          grhoa = grhorefa - rrho*grhorefa
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhoa, grhorefb, grhorefab, potgma, 3 )
          grhob = grhorefb + rrho*grhorefb
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhob, grhorefab, potgpb, 3 )
          grhob = grhorefb - rrho*grhorefb
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhob, grhorefab, potgmb, 3 )
          grhoab = grhorefab + rrho*grhorefab
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhorefb, grhoab, potgpab, 3 )
          grhoab = grhorefab - rrho*grhorefab
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhorefb, grhoab, potgmab, 3 )
       ELSE
          CALL fungga ( derivative_set, i, rhoref, grhoref, potx(:,1:1), 0 )
          CALL fungga ( derivative_set, i, rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
          rhoa = rhorefa + rrho*rhorefa
          rhob = rhorefb
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
          rhoa = rhorefa - rrho*rhorefa
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, potma, 3 )
          rhoa = rhorefa
          rhob = rhorefb + rrho*rhorefb
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
          rhob = rhorefb - rrho*rhorefb
          CALL fungga ( derivative_set, i, rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
          grhoa = grhorefa + rrho*grhorefa
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
          grhoa = grhorefa - rrho*grhorefa
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
          grhob = grhorefb + rrho*grhorefb
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
          grhob = grhorefb - rrho*grhorefb
          CALL fungga ( derivative_set, i, rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
       END IF
       !
       abserr(:) = ABS(potx(:,1)-potx(:,2))
       relerr(:) = abserr(:)/ABS(potx(:,1))
       CALL output ( iw, rhoref, potx(:,1), abserr, relerr, &
                   " Compare GGA-LDA/GGA-LSD Functionals " )
       !
       IF ( crossterm ) THEN
         CALL testder (" First Derivative Ra",pot(:,1),potpa(:,0),potma(:,0),rhorefa)
         CALL testder (" First Derivative Ga",pot(:,2),potgpa(:,0),potgma(:,0),grhorefa)
         CALL testder (" First Derivative Rb",pot(:,3),potpb(:,0),potmb(:,0),rhorefb)
         CALL testder (" First Derivative Gb",pot(:,4),potgpb(:,0),potgmb(:,0),grhorefb)
         CALL testder (" First Derivative Gab",pot(:,5),potgpab(:,0),potgmab(:,0),grhorefab)
       ELSE
         CALL testder (" First Derivative Ra",pot(:,1),potpa(:,0),potma(:,0),rhorefa)
         CALL testder (" First Derivative Ga",pot(:,2),potgpa(:,0),potgma(:,0),grhorefa)
         CALL testder (" First Derivative Rb",pot(:,3),potpb(:,0),potmb(:,0),rhorefb)
         CALL testder (" First Derivative Gb",pot(:,4),potgpb(:,0),potgmb(:,0),grhorefb)
         CALL testder (" Second Derivative RaRa",pot(:,5),potpa(:,1),potma(:,1),rhorefa)
         CALL testder (" Second Derivative RaGa",pot(:,6),potgpa(:,1),potgma(:,1),grhorefa)
         CALL testder (" Second Derivative GaGa",pot(:,7),potgpa(:,2),potgma(:,2),grhorefa)
         CALL testder (" Second Derivative RbRb",pot(:,8),potpb(:,3),potmb(:,3),rhorefb)
         CALL testder (" Second Derivative RbGb",pot(:,9),potgpb(:,3),potgmb(:,3),grhorefb)
         CALL testder (" Second Derivative GbGb",pot(:,10),potgpb(:,4),potgmb(:,4),grhorefb)
         CALL testder (" Third Derivative RaRaRa",pot(:,11),potpa(:,5),potma(:,5),rhorefa)
         CALL testder (" Third Derivative RaRaGa",pot(:,12),potgpa(:,5),potgma(:,5),grhorefa)
         CALL testder (" Third Derivative RaGaGa",pot(:,13),potgpa(:,6),potgma(:,6),grhorefa)
         CALL testder (" Third Derivative GaGaGa",pot(:,14),potgpa(:,7),potgma(:,7),grhorefa)
         CALL testder (" Third Derivative RbRbRb",pot(:,15),potpb(:,8),potmb(:,8),rhorefb)
         CALL testder (" Third Derivative RbRbGb",pot(:,16),potgpb(:,8),potgmb(:,8),grhorefb)
         CALL testder (" Third Derivative RbGbGb",pot(:,17),potgpb(:,9),potgmb(:,9),grhorefb)
         CALL testder (" Third Derivative RbGbGb",pot(:,18),potgpb(:,10),potgmb(:,10),grhorefb)
       END IF
       !
     END DO
     !
     CALL deallocate_arrays

  END IF

END SUBROUTINE xc_functional_test

! *****************************************************************************

SUBROUTINE testder(message,p,pp,pm,rho)
    CHARACTER(LEN=*), INTENT(IN)             :: message
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: p, pp, pm, rho

  abserr(:) = ABS(p(:)-(pp(:)-pm(:))/(2._dbl*rrho*rho))
  relerr(:) = abserr(:)/ABS(p(:)+1.e-30_dbl)
  CALL output ( iw, rho, p(:), abserr, relerr, message )
END SUBROUTINE testder

SUBROUTINE output(iw,rhoref,pot,abserr,relerr,message)
    INTEGER, INTENT(IN)                      :: iw
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoref, pot, abserr, relerr
    CHARACTER(LEN=*), INTENT(IN)             :: message

    INTEGER                                  :: ia, ir

  ia=MAXVAL(MAXLOC(abserr))
  ir=MAXVAL(MAXLOC(relerr))
  IF ( relerr(ir) > rrho .AND. abserr(ir) > 1.e-12_dbl ) THEN
    WRITE (iw,'(A)') message
    WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
               rhoref(ia),pot(ia),abserr(ia),relerr(ia)
    WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
               rhoref(ir),pot(ir),abserr(ir),relerr(ir)
  ENDIF

END SUBROUTINE output

! *****************************************************************************

SUBROUTINE allocate_arrays
    INTEGER                                  :: ierr, mm = 57

     ALLOCATE ( abserr ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "abserr", n )
     abserr = 0._dbl
     ALLOCATE ( relerr ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "relerr", n )
     relerr = 0._dbl
     !
     ALLOCATE ( rhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref", n )
     ALLOCATE ( rho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho", n )
     ALLOCATE ( pot ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot", mm*n )
     ALLOCATE ( potp ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potp", mm*n )
     ALLOCATE ( potm ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potm", mm*n )
     !
     ALLOCATE ( grhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoref", n )
     ALLOCATE ( grho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grho", n )
     ALLOCATE ( potgp ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgp", mm*n )
     ALLOCATE ( potgm ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgm", mm*n )
     !
     ALLOCATE ( rhorefa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefa", n )
     ALLOCATE ( rhorefb ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefb", n )
     ALLOCATE ( grhorefa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefa", n )
     ALLOCATE ( grhorefb ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefb", n )
     ALLOCATE ( grhorefab ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefab", n )
     ALLOCATE ( rhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa", n )
     ALLOCATE ( rhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob", n )
     ALLOCATE ( grhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoa", n )
     ALLOCATE ( grhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhob", n )
     ALLOCATE ( grhoab ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoab", n )
     ALLOCATE ( potx ( n, mm ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potx", 2*mm )
     ALLOCATE ( potpa ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpa", mm*n )
     ALLOCATE ( potpb ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpb", mm*n )
     ALLOCATE ( potma ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potma", mm*n )
     ALLOCATE ( potmb ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potmb", mm*n )
     ALLOCATE ( potgpa ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpa", mm*n )
     ALLOCATE ( potgpb ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpb", mm*n )
     ALLOCATE ( potgpab ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpab", mm*n )
     ALLOCATE ( potgma ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgma", mm*n )
     ALLOCATE ( potgmb ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgmb", mm*n )
     ALLOCATE ( potgmab ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgmab", mm*n )
     !
END SUBROUTINE allocate_arrays

SUBROUTINE deallocate_arrays
    INTEGER                                  :: ierr

     DEALLOCATE ( abserr, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "abserr" )
     DEALLOCATE ( relerr, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "relerr" )
     !
     DEALLOCATE ( rhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref" )
     DEALLOCATE ( rho, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho" )
     DEALLOCATE ( pot, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot" )
     DEALLOCATE ( potp, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potp" )
     DEALLOCATE ( potm, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potm" )
     !
     DEALLOCATE ( grhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoref" )
     DEALLOCATE ( grho, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grho" )
     DEALLOCATE ( potgp, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgp" )
     DEALLOCATE ( potgm, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgm" )
     !
     DEALLOCATE ( rhorefa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefa" )
     DEALLOCATE ( rhorefb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefb" )
     DEALLOCATE ( rhoa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa" )
     DEALLOCATE ( rhob, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob" )
     DEALLOCATE ( grhorefa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefa" )
     DEALLOCATE ( grhorefb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefb" )
     DEALLOCATE ( grhorefab, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefab" )
     DEALLOCATE ( grhoa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoa" )
     DEALLOCATE ( grhob, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhob" )
     DEALLOCATE ( potx, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potx" )
     DEALLOCATE ( potpa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpa" )
     DEALLOCATE ( potpb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpb" )
     DEALLOCATE ( potma, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potma" )
     DEALLOCATE ( potmb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potmb" )
     DEALLOCATE ( potgpa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpa" )
     DEALLOCATE ( potgpb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpb" )
     DEALLOCATE ( potgma, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgma" )
     DEALLOCATE ( potgmb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgmb" )
     !
END SUBROUTINE deallocate_arrays

! *****************************************************************************

SUBROUTINE fun_init ( fun, cutoff, funtype, crossterm )

    INTEGER, INTENT(IN)                      :: fun
    REAL(dbl), INTENT(IN)                    :: cutoff
    INTEGER, INTENT(OUT)                     :: funtype
    LOGICAL, INTENT(OUT)                     :: crossterm

!
!   funtype : 0 nothing implemented
!             1 LDA implemented
!             2 LSD implemented
!             3 LDA+LSD implemented
!            11 GGA/LDA implemented
!            12 GGA/LSD implemented
!            13 GGA/LDA+LSD implemented
!

    SELECT CASE (fun)
    CASE (1)
      CALL xalpha ( cutoff, 2._dbl/3._dbl )
      crossterm=.FALSE.
      funtype=3
    CASE (2)
      CALL pade ( cutoff, debug=.TRUE. )
      crossterm=.TRUE.
      funtype=3
    CASE (3)
      CALL vwn ( cutoff )
      crossterm=.TRUE.
      funtype=1
    CASE (4)
      CALL thomas_fermi ( cutoff )
      crossterm=.FALSE.
      funtype=3
    CASE (5)
      CALL perdew_wang ( "PWO", cutoff )
      crossterm=.TRUE.
      funtype=3
    CASE (6)
      CALL perdew_zunger ( "ORIG", cutoff )
      crossterm=.TRUE.
      funtype=3
    CASE (7)
      CALL tfw ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (8)
      CALL ol1 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (9)
      CALL ol2 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (10)
      CALL llp ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (11)
      CALL ke_pw86 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (12)
      CALL ke_pw91 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (13)
      CALL lc ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (14)
      CALL t92 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (15)
      CALL b88 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (16)
      CALL x_pw86 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (17)
      CALL x_pw91 ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (18)
      CALL pbex ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (19)
      CALL revpbex ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (20)
      CALL optx ( cutoff )
      crossterm=.FALSE.
      funtype=13
    CASE (21)
      CALL lyp ( cutoff, debug=.TRUE. )
      crossterm=.TRUE.
      funtype=11
    CASE (22)
      CALL perdew86 ( cutoff, debug=.TRUE. )
      crossterm=.TRUE.
      funtype=11
    CASE (23)
      CALL cs1 ( cutoff, debug=.TRUE. )
      crossterm=.TRUE.
      funtype=13
    END SELECT

END SUBROUTINE fun_init

SUBROUTINE fun_info ( fun, reference, shortform )

    INTEGER, INTENT(IN)                      :: fun
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    SELECT CASE (fun)
    CASE (1)
      CALL xalpha ( reference, shortform )
    CASE (2)
      CALL pade ( reference, shortform )
    CASE (3)
      CALL vwn ( reference, shortform )
    CASE (4)
      CALL thomas_fermi ( reference, shortform )
    CASE (5)
      CALL perdew_wang ( reference, shortform )
    CASE (6)
      CALL perdew_zunger ( reference, shortform )
    CASE (7)
      CALL tfw ( reference, shortform )
    CASE (8)
      CALL ol1 ( reference, shortform )
    CASE (9)
      CALL ol2 ( reference, shortform )
    CASE (10)
      CALL llp ( reference, shortform )
    CASE (11)
      CALL ke_pw86 ( reference, shortform )
    CASE (12)
      CALL ke_pw91 ( reference, shortform )
    CASE (13)
      CALL lc ( reference, shortform )
    CASE (14)
      CALL t92 ( reference, shortform )
    CASE (15)
      CALL b88 ( reference, shortform )
    CASE (16)
      CALL x_pw86 ( reference, shortform )
    CASE (17)
      CALL x_pw91 ( reference, shortform )
    CASE (18)
      CALL pbex ( reference, shortform )
    CASE (19)
      CALL revpbex ( reference, shortform )
    CASE (20)
      CALL optx ( reference, shortform )
    CASE (21)
      CALL lyp ( reference, shortform )
    CASE (22)
      CALL perdew86 ( reference, shortform )
    CASE (23)
      CALL cs1 ( reference, shortform )
    END SELECT

END SUBROUTINE fun_info

SUBROUTINE fun_lda ( derivative_set, fun, rho, pot, order )

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    SELECT CASE (fun)
    CASE (1)
      pot = 0._dbl
      CALL xalpha ( derivative_set, rho, pot, order )
    CASE (2)
      pot = 0._dbl
      CALL pade ( derivative_set, rho, pot, order )
    CASE (3)
      pot = 0._dbl
      CALL vwn ( derivative_set, rho, pot, order )
    CASE (4)
      pot = 0._dbl
      CALL thomas_fermi ( derivative_set, rho, pot, order )
    CASE (5)
      pot = 0._dbl
      CALL perdew_wang ( derivative_set, rho, pot, order )
    CASE (6)
      pot = 0._dbl
      CALL perdew_zunger ( derivative_set, rho, pot, order )
    CASE DEFAULT
      CALL stop_program("fun_lda","functional not found")
    END SELECT

END SUBROUTINE fun_lda

SUBROUTINE fun_lsd ( derivative_set, fun, rhoa, rhob, pot, order )

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    SELECT CASE (fun)
    CASE (1)
      pot = 0._dbl
      CALL xalpha ( derivative_set, rhoa, rhob, pot, order )
    CASE (2)
      pot = 0._dbl
      CALL pade ( derivative_set, rhoa, rhob, pot, order )
    CASE (4)
      pot = 0._dbl
      CALL thomas_fermi ( derivative_set, rhoa, rhob, pot, order )
    CASE (5)
      pot = 0._dbl
      CALL perdew_wang ( derivative_set, rhoa, rhob, pot, order )
    CASE (6)
      pot = 0._dbl
      CALL perdew_zunger ( derivative_set, rhoa, rhob, pot, order )
    CASE DEFAULT
      CALL stop_program("fun_lsd","functional not found")
    END SELECT

END SUBROUTINE fun_lsd

SUBROUTINE funggau ( derivative_set, fun, rho, grho, pot, order )

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    SELECT CASE (fun)
    CASE (7)
      pot = 0._dbl
      CALL tfw ( derivative_set, rho, grho, pot, order )
    CASE (8)
      pot = 0._dbl
      CALL ol1 ( derivative_set, rho, grho, pot, order )
    CASE (9)
      pot = 0._dbl
      CALL ol2 ( derivative_set, rho, grho, pot, order )
    CASE (10)
      pot = 0._dbl
      CALL llp ( derivative_set, rho, grho, pot, order )
    CASE (11)
      pot = 0._dbl
      CALL ke_pw86 ( derivative_set, rho, grho, pot, order )
    CASE (12)
      pot = 0._dbl
      CALL ke_pw91 ( derivative_set, rho, grho, pot, order )
    CASE (13)
      pot = 0._dbl
      CALL lc ( derivative_set, rho, grho, pot, order )
    CASE (14)
      pot = 0._dbl
      CALL t92 ( derivative_set, rho, grho, pot, order )
    CASE (15)
      pot = 0._dbl
      CALL b88 ( derivative_set, rho, grho, pot, order )
    CASE (16)
      pot = 0._dbl
      CALL x_pw86 ( derivative_set, rho, grho, pot, order )
    CASE (17)
      pot = 0._dbl
      CALL x_pw91 ( derivative_set, rho, grho, pot, order )
    CASE (18)
      pot = 0._dbl
      CALL pbex ( derivative_set, rho, grho, pot, order )
    CASE (19)
      pot = 0._dbl
      CALL revpbex ( derivative_set, rho, grho, pot, order )
    CASE (20)
      pot = 0._dbl
      CALL optx ( derivative_set, rho, grho, pot, order )
    CASE (21)
      pot = 0._dbl
      CALL lyp ( derivative_set, rho, grho, pot, order )
    CASE (22)
      pot = 0._dbl
      CALL perdew86 ( derivative_set, rho, grho, pot, order )
    CASE (23)
      pot = 0._dbl
      CALL cs1 ( derivative_set, rho, grho, pot, order )
    CASE DEFAULT
      CALL stop_program("funggau","functional not found")
    END SELECT

END SUBROUTINE funggau

SUBROUTINE funggap ( derivative_set, fun, rhoa, rhob, grhoa, grhob, pot, order )

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    SELECT CASE (fun)
    CASE (7)
      pot = 0._dbl
      CALL tfw ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (8)
      pot = 0._dbl
      CALL ol1 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (9)
      pot = 0._dbl
      CALL ol2 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (10)
      pot = 0._dbl
      CALL llp ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (11)
      pot = 0._dbl
      CALL ke_pw86 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (12)
      pot = 0._dbl
      CALL ke_pw91 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (13)
      pot = 0._dbl
      CALL lc ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (14)
      pot = 0._dbl
      CALL t92 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (15)
      pot = 0._dbl
      CALL b88 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (16)
      pot = 0._dbl
      CALL x_pw86 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (17)
      pot = 0._dbl
      CALL x_pw91 ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (18)
      pot = 0._dbl
      CALL pbex ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (19)
      pot = 0._dbl
      CALL revpbex ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE (20)
      pot = 0._dbl
      CALL optx ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    CASE DEFAULT
      CALL stop_program("funggap","functional not found")
    END SELECT

END SUBROUTINE funggap

SUBROUTINE funggab ( derivative_set, fun, rhoa, rhob, grhoa, grhob, grhoab, pot, order )

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob, &
                                                grhoab
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    SELECT CASE (fun)
    CASE (23)
      pot = 0._dbl
      CALL cs1 ( derivative_set, rhoa, rhob, grhoa, grhob, grhoab, pot, order )
    CASE DEFAULT
      CALL stop_program("funggap","functional not found")
    END SELECT

END SUBROUTINE funggab

! *****************************************************************************

END MODULE xc_functionals_tests

! *****************************************************************************
