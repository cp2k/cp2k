!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/functionals_tests [1.0] *
!!
!!   NAME
!!     functionals_tests
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

MODULE functionals_tests

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE global_types, ONLY : global_environment_type
  USE exchange_gga_functional, ONLY : b88, x_pw86, x_pw91, pbex, revpbex, optx
  USE ke_gga_functional, ONLY : ol1, ol2, llp, ke_pw86, ke_pw91, lc, t92
  USE lyp_functional, ONLY : lyp
  USE pade_functional, ONLY : pade
  USE xalpha_functional, ONLY : xalpha
  USE thomas_fermi_functional, ONLY : thomas_fermi
  USE tfw_functional, ONLY : tfw
  USE termination, ONLY : stop_memory, stop_program
  USE vwn_functional, ONLY : vwn
  USE perdew_wang_functional, ONLY : perdew_wang
  USE perdew_zunger_functional, ONLY : perdew_zunger
  USE perdew86_functional, ONLY : perdew86
  USE cs1_functional, ONLY : cs1

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: functional_test

  INTEGER, PARAMETER :: n=20
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rhoref, rhorefa, rhorefb
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rho, rhoa, rhob, grho, grhoref
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: grhoa, grhob, grhorefa, grhorefb
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: pot, potp, potm, potx
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgp, potgm
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgpa, potgpb, potgma, potgmb
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potpa, potpb, potma, potmb
  REAL(dbl), PARAMETER :: rrho=1.E-5_dbl, epsilon=1.E-15_dbl
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: relerr, abserr

  INTEGER :: iw

  INTEGER, PARAMETER :: nfun = 23  ! number of implemented functionals

  INTERFACE funlda
    MODULE PROCEDURE fun_lda, fun_lsd, fun_init, fun_info
  END INTERFACE

  INTERFACE fungga
    MODULE PROCEDURE funggau, funggap, fun_init, fun_info
  END INTERFACE

! *****************************************************************************

CONTAINS

!******************************************************************************
!!****** functionals_tests/functional_test [1.0] *
!!
!!   NAME
!!     functional_test
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

SUBROUTINE functional_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  INTEGER :: i, j, ia, ir, funtype
  REAL(dbl) :: drho, r
  CHARACTER (LEN=80) :: fun_string
  LOGICAL :: crossterm

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
       CALL funlda ( i, rhoref, pot, 3 )
       rho = rhoref + rrho*rhoref
       CALL funlda ( i, rho, potp, 3 )
       rho = rhoref - rrho*rhoref
       CALL funlda ( i, rho, potm, 3 )
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
       CALL fungga ( i, rhoref, grhoref, pot, 3 )
       rho = rhoref + rrho*rhoref
       CALL fungga ( i, rho, grhoref, potp, 3 )
       rho = rhoref - rrho*rhoref
       CALL fungga ( i, rho, grhoref, potm, 3 )
       grho = grhoref + rrho*grhoref
       CALL fungga ( i, rhoref, grho, potgp, 3 )
       grho = grhoref - rrho*grhoref
       CALL fungga ( i, rhoref, grho, potgm, 3 )
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
       rhorefb(i)=3.e3 * EXP(-r)*(1._dbl-0.1_dbl*sin(r))
     END DO
     rhoref = 2._dbl*rhorefa
     !
     DO i = 1, nfun
       CALL funlda ( i, epsilon, funtype, crossterm )
       IF ( funtype /= 2 .AND. funtype /= 3 ) CYCLE
       CALL funlda ( i, reference=fun_string )
       WRITE (iw,'(A80)') ADJUSTR(fun_string)
       CALL funlda ( i, rhoref, potx(:,1:1), 0 )
       CALL funlda ( i, rhorefa, rhorefa, potx(:,2:2), 0 )
       CALL funlda ( i, rhorefa, rhorefb, pot, 3 )
       rhoa = rhorefa + rrho*rhorefa
       rhob = rhorefb
       CALL funlda ( i, rhoa, rhob, potpa, 3 )
       rhoa = rhorefa - rrho*rhorefa
       CALL funlda ( i, rhoa, rhob, potma, 3 )
       rhoa = rhorefa
       rhob = rhorefb + rrho*rhorefb
       CALL funlda ( i, rhoa, rhob, potpb, 3 )
       rhob = rhorefb - rrho*rhorefb
       CALL funlda ( i, rhoa, rhob, potmb, 3 )
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
       rhorefb(i)=3.e3 * EXP(-r)*(1._dbl-0.1*sin(r))
       grhorefb(i)=7.e3 * EXP(-r)*(1._dbl-0.1*sin(r))
     END DO
     rhoref = 2._dbl*rhorefa
     grhoref= 2._dbl*grhorefa
     !
     DO i = 1, nfun
       CALL fungga ( i, epsilon, funtype, crossterm )
       IF ( funtype /= 12 .AND. funtype /= 13 ) CYCLE
       CALL fungga ( i, reference=fun_string )
       WRITE (iw,'(A80)') ADJUSTR(fun_string)
       CALL fungga ( i, rhoref, grhoref, potx(:,1:1), 0 )
       CALL fungga ( i, rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
       CALL fungga ( i, rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
       rhoa = rhorefa + rrho*rhorefa
       rhob = rhorefb
       CALL fungga ( i, rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
       rhoa = rhorefa - rrho*rhorefa
       CALL fungga ( i, rhoa, rhob, grhorefa, grhorefb, potma, 3 )
       rhoa = rhorefa
       rhob = rhorefb + rrho*rhorefb
       CALL fungga ( i, rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
       rhob = rhorefb - rrho*rhorefb
       CALL fungga ( i, rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
       grhoa = grhorefa + rrho*grhorefa
       CALL fungga ( i, rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
       grhoa = grhorefa - rrho*grhorefa
       CALL fungga ( i, rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
       grhob = grhorefb + rrho*grhorefb
       CALL fungga ( i, rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
       grhob = grhorefb - rrho*grhorefb
       CALL fungga ( i, rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
       !
       abserr(:) = ABS(potx(:,1)-potx(:,2))
       relerr(:) = abserr(:)/ABS(potx(:,1))
       CALL output ( iw, rhoref, potx(:,1), abserr, relerr, &
                   " Compare GGA-LDA/GGA-LSD Functionals " )
       !
       IF ( crossterm ) THEN
         CALL stop_program("test","not available")
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

END SUBROUTINE functional_test

! *****************************************************************************

SUBROUTINE testder(message,p,pp,pm,rho)
  CHARACTER (LEN=*), INTENT(IN) :: message
  REAL(dbl), DIMENSION(:), INTENT(IN) :: p,pp,pm,rho
  abserr(:) = ABS(p(:)-(pp(:)-pm(:))/(2._dbl*rrho*rho))
  relerr(:) = abserr(:)/ABS(p(:)+1.e-30_dbl)
  CALL output ( iw, rho, p(:), abserr, relerr, message )
END SUBROUTINE testder

SUBROUTINE output(iw,rhoref,pot,abserr,relerr,message)
  INTEGER, INTENT(IN) :: iw
  REAL(dbl), DIMENSION(:), INTENT(IN) :: rhoref, pot, abserr, relerr
  CHARACTER (LEN=*), INTENT(IN) :: message

  INTEGER :: ia,ir

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
     INTEGER :: ierr
     INTEGER :: mm=19

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
     ALLOCATE ( rhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa", n )
     ALLOCATE ( rhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob", n )
     ALLOCATE ( grhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoa", n )
     ALLOCATE ( grhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhob", n )
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
     ALLOCATE ( potgma ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgma", mm*n )
     ALLOCATE ( potgmb ( n, 0:mm-1 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgmb", mm*n )
     !
END SUBROUTINE allocate_arrays

SUBROUTINE deallocate_arrays
     INTEGER :: ierr

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

    INTEGER, INTENT(IN) :: fun
    REAL(dbl), INTENT(IN) :: cutoff
    INTEGER, INTENT(OUT) :: funtype
    LOGICAL, INTENT(OUT) :: crossterm

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
      CALL perdew_zunger ( "DMC", cutoff )
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
      funtype=11
    END SELECT

END SUBROUTINE fun_init

SUBROUTINE fun_info ( fun, reference, shortform )

    INTEGER, INTENT(IN) :: fun
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

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

SUBROUTINE fun_lda ( fun, rho, pot, order )

    INTEGER, INTENT(IN) :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    SELECT CASE (fun)
    CASE (1)
      CALL xalpha ( rho, pot, order )
    CASE (2)
      CALL pade ( rho, pot, order )
    CASE (3)
      CALL vwn ( rho, pot, order )
    CASE (4)
      CALL thomas_fermi ( rho, pot, order )
    CASE (5)
      CALL perdew_wang ( rho, pot, order )
    CASE (6)
      CALL perdew_zunger ( rho, pot, order )
    CASE DEFAULT
    END SELECT

END SUBROUTINE fun_lda

SUBROUTINE fun_lsd ( fun, rhoa, rhob, pot, order )

    INTEGER, INTENT(IN) :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    SELECT CASE (fun)
    CASE (1)
      CALL xalpha ( rhoa, rhob, pot, order )
    CASE (2)
      CALL pade ( rhoa, rhob, pot, order )
    CASE (4)
      CALL thomas_fermi ( rhoa, rhob, pot, order )
    CASE (5)
      CALL perdew_wang ( rhoa, rhob, pot, order )
    CASE DEFAULT
    END SELECT

END SUBROUTINE fun_lsd

SUBROUTINE funggau ( fun, rho, grho, pot, order )

    INTEGER, INTENT(IN) :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    SELECT CASE (fun)
    CASE (7)
      CALL tfw ( rho, grho, pot, order )
    CASE (8)
      CALL ol1 ( rho, grho, pot, order )
    CASE (9)
      CALL ol2 ( rho, grho, pot, order )
    CASE (10)
      CALL llp ( rho, grho, pot, order )
    CASE (11)
      CALL ke_pw86 ( rho, grho, pot, order )
    CASE (12)
      CALL ke_pw91 ( rho, grho, pot, order )
    CASE (13)
      CALL lc ( rho, grho, pot, order )
    CASE (14)
      CALL t92 ( rho, grho, pot, order )
    CASE (15)
      CALL b88 ( rho, grho, pot, order )
    CASE (16)
      CALL x_pw86 ( rho, grho, pot, order )
    CASE (17)
      CALL x_pw91 ( rho, grho, pot, order )
    CASE (18)
      CALL pbex ( rho, grho, pot, order )
    CASE (19)
      CALL revpbex ( rho, grho, pot, order )
    CASE (20)
      CALL optx ( rho, grho, pot, order )
    CASE (21)
      CALL lyp ( rho, grho, pot, order )
    CASE (22)
      CALL perdew86 ( rho, grho, pot, order )
    CASE (23)
      CALL cs1 ( rho, grho, pot, order )
    CASE DEFAULT
    END SELECT

END SUBROUTINE funggau

SUBROUTINE funggap ( fun, rhoa, rhob, grhoa, grhob, pot, order )

    INTEGER, INTENT(IN) :: fun
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    SELECT CASE (fun)
    CASE (7)
      CALL tfw ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (8)
      CALL ol1 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (9)
      CALL ol2 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (10)
      CALL llp ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (11)
      CALL ke_pw86 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (12)
      CALL ke_pw91 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (13)
      CALL lc ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (14)
      CALL t92 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (15)
      CALL b88 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (16)
      CALL x_pw86 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (17)
      CALL x_pw91 ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (18)
      CALL pbex ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (19)
      CALL revpbex ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE (20)
      CALL optx ( rhoa, rhob, grhoa, grhob, pot, order )
    CASE DEFAULT
    END SELECT

END SUBROUTINE funggap

! *****************************************************************************

END MODULE functionals_tests

! *****************************************************************************
