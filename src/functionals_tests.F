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
!!*****
!******************************************************************************

MODULE functionals_tests

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE global_types, ONLY : global_environment_type
  USE pade_functional, ONLY : pade
  USE xalpha_functional, ONLY : xalpha
  USE thomas_fermi_functional, ONLY : thomas_fermi
  USE tfw_functional, ONLY : tfw
  USE termination, ONLY : stop_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: functional_test

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

  INTEGER, PARAMETER :: n=20
  INTEGER :: iw, ierr, i, nfun, j, ia, ir
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rhoref, rhorefa, rhorefb
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rho, rhoa, rhob, grho, grhoref
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: grhoa, grhob, grhorefa, grhorefb
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: pot, potp, potm, potx
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgp, potgm
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potgpa, potgpb, potgma, potgmb
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: potpa, potpb, potma, potmb
  REAL(dbl), PARAMETER :: rrho=1.E-5_dbl, epsilon=5.E-13_dbl
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: relerr, abserr
  REAL(dbl) :: drho
  CHARACTER (LEN=80) :: fun_string
  LOGICAL :: crossterm

! *****************************************************************************

  iw =globenv % scr

  IF ( globenv % ionode ) THEN
     ALLOCATE ( abserr ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "abserr", n )
     abserr = 0._dbl
     ALLOCATE ( relerr ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "relerr", n )
     relerr = 0._dbl
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test LDA functionals
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF LDA FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     ALLOCATE ( rhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref", n )
     ALLOCATE ( rho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho", n )
     ALLOCATE ( pot ( n, 0:3 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot", 4*n )
     ALLOCATE ( potp ( n, 0:3 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potp", 4*n )
     ALLOCATE ( potm ( n, 0:3 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potm", 4*n )
     !
     rhoref(1)=1.E-12_dbl
     DO i=2,n
       rhoref(i)=rhoref(i-1)*6._dbl
     END DO
     !
     nfun = 3
     DO i = 1, nfun

       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** X-ALPHA FUNCTIONAL "
         CALL xalpha ( epsilon )
         CALL xalpha ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL xalpha ( rhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL xalpha ( rho, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL xalpha ( rho, potm, 3 )
       CASE (2)
         WRITE (iw,'(/,A)') " *** PADE FUNCTIONAL "
         CALL pade ( epsilon )
         CALL pade ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string) 
         CALL pade ( rhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL pade ( rho, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL pade ( rho, potm, 3 )
       CASE (3)
         WRITE (iw,'(/,A)') " *** THOMAS-FERMI FUNCTIONAL "
         CALL thomas_fermi ( epsilon )
         CALL thomas_fermi ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL thomas_fermi ( rhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL thomas_fermi ( rho, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL thomas_fermi ( rho, potm, 3 )
       END SELECT
       !
       WRITE (iw,'(/,A)') " First Derivative "
       abserr(:) = ABS(pot(:,1)-(potp(:,0)-potm(:,0))/(2._dbl*rrho*rhoref))
       relerr(:) = abserr(:)/ABS(pot(:,1))
       CALL output ( iw, rhoref, pot(:,1), abserr, relerr )
       !
       WRITE (iw,'(/,A)') " Second Derivative "
       abserr(:) = ABS(pot(:,2)-(potp(:,1)-potm(:,1))/(2._dbl*rrho*rhoref))
       relerr(:) = abserr(:)/ABS(pot(:,2))
       CALL output ( iw, rhoref, pot(:,2), abserr, relerr )
       !
       WRITE (iw,'(/,A)') " Third Derivative "
       abserr(:) = ABS(pot(:,3)-(potp(:,2)-potm(:,2))/(2._dbl*rrho*rhoref))
       relerr(:) = abserr(:)/ABS(pot(:,3))
       CALL output ( iw, rhoref, pot(:,3), abserr, relerr )
     END DO
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

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test GGA-LDA functionals
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF GGA-LDA FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     ALLOCATE ( rhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref", n )
     ALLOCATE ( grhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoref", n )
     ALLOCATE ( rho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho", n )
     ALLOCATE ( grho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grho", n )
     ALLOCATE ( pot ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot", 10*n )
     ALLOCATE ( potp ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potp", 10*n )
     ALLOCATE ( potm ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potm", 10*n )
     ALLOCATE ( potgp ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgp", 10*n )
     ALLOCATE ( potgm ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgm", 10*n )
     !
     rhoref(1)=1.E-12_dbl
     grhoref(1)=1.E-13_dbl
     DO i=2,n
       rhoref(i)=rhoref(i-1)*6._dbl
       grhoref(i)=grhoref(i-1)*8._dbl
     END DO
     !
     nfun = 1
     DO i = 1, nfun
       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** Thomas-Fermi-von Weizsaecker Functional"
         CALL tfw ( epsilon )
         CALL tfw ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL tfw ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL tfw ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL tfw ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL tfw ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL tfw ( rhoref, grho, potgm, 3 )
       END SELECT
     END DO
     !
     WRITE (iw,'(/,A)') " First Derivative R"
     abserr(:) = ABS(pot(:,1)-(potp(:,0)-potm(:,0))/(2._dbl*rrho*rhoref))
     relerr(:) = abserr(:)/ABS(pot(:,1)+1.e-30)
     CALL output ( iw, rhoref, pot(:,1), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " First Derivative G"
     abserr(:) = ABS(pot(:,2)-(potgp(:,0)-potgm(:,0))/(2._dbl*rrho*grhoref))
     relerr(:) = abserr(:)/ABS(pot(:,2)+1.e-30)
     CALL output ( iw, grhoref, pot(:,2), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Second Derivative RR"
     abserr(:) = ABS(pot(:,3)-(potp(:,1)-potm(:,1))/(2._dbl*rrho*rhoref))
     relerr(:) = abserr(:)/ABS(pot(:,3)+1.e-30)
     CALL output ( iw, rhoref, pot(:,3), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Second Derivative RG"
     abserr(:) = ABS(pot(:,4)-(potgp(:,1)-potgm(:,1))/(2._dbl*rrho*grhoref))
     relerr(:) = abserr(:)/ABS(pot(:,4)+1.e-30)
     CALL output ( iw, grhoref, pot(:,4), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Second Derivative GG"
     abserr(:) = ABS(pot(:,5)-(potgp(:,2)-potgm(:,2))/(2._dbl*rrho*grhoref))
     relerr(:) = abserr(:)/ABS(pot(:,5)+1.e-30)
     CALL output ( iw, grhoref, pot(:,5), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Third Derivative RRR"
     abserr(:) = ABS(pot(:,6)-(potp(:,3)-potm(:,3))/(2._dbl*rrho*rhoref))
     relerr(:) = abserr(:)/ABS(pot(:,6)+1.e-30)
     CALL output ( iw, rhoref, pot(:,6), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Third Derivative RRG"
     abserr(:) = ABS(pot(:,7)-(potgp(:,3)-potgm(:,3))/(2._dbl*rrho*grhoref))
     relerr(:) = abserr(:)/ABS(pot(:,7)+1.e-30)
     CALL output ( iw, grhoref, pot(:,7), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Third Derivative RGG"
     abserr(:) = ABS(pot(:,8)-(potgp(:,4)-potgm(:,4))/(2._dbl*rrho*grhoref))
     relerr(:) = abserr(:)/ABS(pot(:,8)+1.e-30)
     CALL output ( iw, grhoref, pot(:,8), abserr, relerr )
     !
     WRITE (iw,'(/,A)') " Third Derivative GGG"
     abserr(:) = ABS(pot(:,9)-(potgp(:,5)-potgm(:,5))/(2._dbl*rrho*grhoref))
     relerr(:) = abserr(:)/ABS(pot(:,9)+1.e-30)
     CALL output ( iw, grhoref, pot(:,9), abserr, relerr )
     !
     DEALLOCATE ( rhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref" )
     DEALLOCATE ( rho, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho" )
     DEALLOCATE ( grhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoref" )
     DEALLOCATE ( grho, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grho" )
     DEALLOCATE ( pot, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot" )
     DEALLOCATE ( potp, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potp" )
     DEALLOCATE ( potm, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potm" )
     DEALLOCATE ( potgp, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgp" )
     DEALLOCATE ( potgm, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgm" )

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test LSD functionals
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF LSD FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     ALLOCATE ( rhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref", n )
     ALLOCATE ( rhorefa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefa", n )
     ALLOCATE ( rhorefb ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefb", n )
     ALLOCATE ( rhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa", n )
     ALLOCATE ( rhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob", n )
     ALLOCATE ( pot ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot", 10*n )
     ALLOCATE ( potpa ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpa", 10*n )
     ALLOCATE ( potpb ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpb", 10*n )
     ALLOCATE ( potma ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potma", 10*n )
     ALLOCATE ( potmb ( n, 0:9 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potmb", 10*n )
     ALLOCATE ( potx ( n, 2 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potx", 2*n )
     !
     rhorefa(1)=1.E-12_dbl
     DO i=2,n
       rhorefa(i)=rhorefa(i-1)*6._dbl
     END DO
     rhorefb(1)=5.E-12_dbl
     DO i=2,n/2
       rhorefb(i)=rhorefb(i-1)*7._dbl
     END DO
     DO i=n/2+1,n
       rhorefb(i)=rhorefb(i-1)*5._dbl
     END DO
     rhoref = 2._dbl*rhorefa
     !
     nfun = 3
     DO i = 1, nfun
       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** X-ALPHA FUNCTIONAL "
         CALL xalpha ( epsilon )
         CALL xalpha ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL xalpha ( rhoref, potx(:,1:1), 0 )
         CALL xalpha ( rhorefa, rhorefa, potx(:,2:2), 0 )
         CALL xalpha ( rhorefa, rhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL xalpha ( rhoa, rhob, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL xalpha ( rhoa, rhob, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL xalpha ( rhoa, rhob, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL xalpha ( rhoa, rhob, potmb, 3 )
         crossterm = .FALSE.
       CASE (2)
         WRITE (iw,'(/,A)') " *** PADE FUNCTIONAL "
         CALL pade ( epsilon )
         CALL pade ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL pade ( rhoref, potx(:,1:1), 0 )
         CALL pade ( rhorefa, rhorefa, potx(:,2:2), 0 )
         CALL pade ( rhorefa, rhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL pade ( rhoa, rhob, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL pade ( rhoa, rhob, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL pade ( rhoa, rhob, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL pade ( rhoa, rhob, potmb, 3 )
         crossterm = .TRUE.
       CASE (3)
         WRITE (iw,'(/,A)') " *** THOMAS-FERMI FUNCTIONAL "
         CALL thomas_fermi ( epsilon )
         CALL thomas_fermi ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL thomas_fermi ( rhoref, potx(:,1:1), 0 )
         CALL thomas_fermi ( rhorefa, rhorefa, potx(:,2:2), 0 )
         CALL thomas_fermi ( rhorefa, rhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL thomas_fermi ( rhoa, rhob, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL thomas_fermi ( rhoa, rhob, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL thomas_fermi ( rhoa, rhob, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL thomas_fermi ( rhoa, rhob, potmb, 3 )
         crossterm = .FALSE.
       END SELECT
       WRITE (iw,'(/,A)') " Compare LDA/LSD Functionals "
       abserr(:) = ABS(potx(:,1)-potx(:,2))
       relerr(:) = abserr(:)/ABS(potx(:,1))
       CALL output ( iw, rhoref, potx(:,1), abserr, relerr )
       !
       WRITE (iw,'(/,A)') " First Derivative A"
       abserr(:) = ABS(pot(:,1)-(potpa(:,0)-potma(:,0))/(2._dbl*rrho*rhorefa))
       relerr(:) = abserr(:)/ABS(pot(:,1))
       CALL output ( iw, rhorefa, pot(:,1), abserr, relerr )
       !
       WRITE (iw,'(/,A)') " First Derivative B"
       abserr(:) = ABS(pot(:,2)-(potpb(:,0)-potmb(:,0))/(2._dbl*rrho*rhorefb))
       relerr(:) = abserr(:)/ABS(pot(:,2))
       CALL output ( iw, rhorefb, pot(:,2), abserr, relerr )
       !
       WRITE (iw,'(/,A)') " Second Derivative AA"
       abserr(:) = ABS(pot(:,3)-(potpa(:,1)-potma(:,1))/(2._dbl*rrho*rhorefa))
       relerr(:) = abserr(:)/ABS(pot(:,3))
       CALL output ( iw, rhorefa, pot(:,3), abserr, relerr )
       !
       IF ( crossterm ) THEN
          WRITE (iw,'(/,A)') " Second Derivative AB"
          abserr(:) = ABS(pot(:,4)-(potpa(:,2)-potma(:,2))/(2._dbl*rrho*rhorefa))
          relerr(:) = abserr(:)/ABS(pot(:,4))
          CALL output ( iw, rhorefa, pot(:,4), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Second Derivative BB"
          abserr(:) = ABS(pot(:,5)-(potpb(:,2)-potmb(:,2))/(2._dbl*rrho*rhorefb))
          relerr(:) = abserr(:)/ABS(pot(:,5))
          CALL output ( iw, rhorefb, pot(:,5), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Third Derivative AAA"
          abserr(:) = ABS(pot(:,6)-(potpa(:,3)-potma(:,3))/(2._dbl*rrho*rhorefa))
          relerr(:) = abserr(:)/ABS(pot(:,6))
          CALL output ( iw, rhorefa, pot(:,6), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Third Derivative AAB"
          abserr(:) = ABS(pot(:,7)-(potpb(:,3)-potmb(:,3))/(2._dbl*rrho*rhorefb))
          relerr(:) = abserr(:)/ABS(pot(:,7))
          CALL output ( iw, rhorefb, pot(:,7), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Third Derivative ABB"
          abserr(:) = ABS(pot(:,8)-(potpa(:,5)-potma(:,5))/(2._dbl*rrho*rhorefa))
          relerr(:) = abserr(:)/ABS(pot(:,8))
          CALL output ( iw, rhorefa, pot(:,8), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Third Derivative BBB"
          abserr(:) = ABS(pot(:,9)-(potpb(:,5)-potmb(:,5))/(2._dbl*rrho*rhorefb))
          relerr(:) = abserr(:)/ABS(pot(:,9))
          CALL output ( iw, rhorefa, pot(:,9), abserr, relerr )
       ELSE
          !
          WRITE (iw,'(/,A)') " Second Derivative BB"
          abserr(:) = ABS(pot(:,4)-(potpb(:,2)-potmb(:,2))/(2._dbl*rrho*rhorefb))
          relerr(:) = abserr(:)/ABS(pot(:,4))
          CALL output ( iw, rhorefa, pot(:,4), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Third Derivative AAA"
          abserr(:) = ABS(pot(:,5)-(potpa(:,3)-potma(:,3))/(2._dbl*rrho*rhorefa))
          relerr(:) = abserr(:)/ABS(pot(:,5))
          CALL output ( iw, rhorefa, pot(:,5), abserr, relerr )
          !
          WRITE (iw,'(/,A)') " Third Derivative BBB"
          abserr(:) = ABS(pot(:,6)-(potpb(:,4)-potmb(:,4))/(2._dbl*rrho*rhorefb))
          relerr(:) = abserr(:)/ABS(pot(:,6))
          CALL output ( iw, rhorefa, pot(:,6), abserr, relerr )
       END IF
     END DO
     !
     DEALLOCATE ( rhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref" )
     DEALLOCATE ( rhorefa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefa" )
     DEALLOCATE ( rhorefb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefb" )
     DEALLOCATE ( rhoa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa" )
     DEALLOCATE ( rhob, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob" )
     DEALLOCATE ( pot, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot" )
     DEALLOCATE ( potpa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpa" )
     DEALLOCATE ( potpb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpb" )
     DEALLOCATE ( potma, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potma" )
     DEALLOCATE ( potmb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potmb" )
     DEALLOCATE ( potx, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potx" )

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Test GGA-LSD functionals
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE (iw,'(/,A,/)') " *** TEST OF GGA-LSD FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
     ALLOCATE ( rhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref", n )
     ALLOCATE ( rhorefa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefa", n )
     ALLOCATE ( rhorefb ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefb", n )
     ALLOCATE ( grhoref ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoref", n )
     ALLOCATE ( grhorefa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefa", n )
     ALLOCATE ( grhorefb ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefb", n )
     ALLOCATE ( rho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho", n )
     ALLOCATE ( rhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa", n )
     ALLOCATE ( rhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob", n )
     ALLOCATE ( grho ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grho", n )
     ALLOCATE ( grhoa ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoa", n )
     ALLOCATE ( grhob ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhob", n )
     ALLOCATE ( pot ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot", 19*n )
     ALLOCATE ( potx ( n, 2 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potx", 2*n )
     ALLOCATE ( potpa ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpa", 19*n )
     ALLOCATE ( potpb ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potpb", 19*n )
     ALLOCATE ( potma ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potma", 19*n )
     ALLOCATE ( potmb ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potmb", 19*n )
     ALLOCATE ( potgpa ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpa", 19*n )
     ALLOCATE ( potgpb ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgpb", 19*n )
     ALLOCATE ( potgma ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgma", 19*n )
     ALLOCATE ( potgmb ( n, 0:18 ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "potgmb", 19*n )
     !
     rhorefa(1)=1.E-12_dbl
     grhorefa(1)=1.E-13_dbl
     grhorefb(1)=1.E-13_dbl
     DO i=2,n
       rhorefa(i)=rhorefa(i-1)*6._dbl
       grhorefa(i)=grhorefa(i-1)*8._dbl
       grhorefb(i)=grhorefb(i-1)*7.77_dbl
     END DO
     rhorefb(1)=5.E-12_dbl
     DO i=2,n/2
       rhorefb(i)=rhorefb(i-1)*7._dbl
     END DO
     DO i=n/2+1,n
       rhorefb(i)=rhorefb(i-1)*5._dbl
     END DO
     rhoref = 2._dbl*rhorefa
     grhoref= 2._dbl*grhorefa
     !
     nfun = 1
     DO i = 1, nfun
       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** Thomas-Fermi-von Weizsaecker Functional"
         CALL tfw ( epsilon )
         CALL tfw ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL tfw ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL tfw ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL tfw ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL tfw ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL tfw ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL tfw ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL tfw ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )

         grhoa = grhorefa + rrho*grhorefa
         CALL tfw ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL tfw ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL tfw ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL tfw ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       END SELECT
       !
       WRITE (iw,'(/,A)') " Compare GGA-LDA/GGA-LSD Functionals "
       abserr(:) = ABS(potx(:,1)-potx(:,2))
       relerr(:) = abserr(:)/ABS(potx(:,1))
       CALL output ( iw, rhoref, potx(:,1), abserr, relerr )
       !
       IF ( crossterm ) THEN
       ELSE
         WRITE (iw,'(/,A)') " First Derivative Ra"
         abserr(:) = ABS(pot(:,1)-(potpa(:,0)-potma(:,0))/(2._dbl*rrho*rhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,1)+1.e-30)
         CALL output ( iw, rhorefa, pot(:,1), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " First Derivative Ga"
         abserr(:) = ABS(pot(:,2)-(potgpa(:,0)-potgma(:,0))/(2._dbl*rrho*grhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,2)+1.e-30)
         CALL output ( iw, grhorefa, pot(:,2), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " First Derivative Rb"
         abserr(:) = ABS(pot(:,3)-(potpb(:,0)-potmb(:,0))/(2._dbl*rrho*rhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,3)+1.e-30)
         CALL output ( iw, rhorefb, pot(:,3), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " First Derivative Gb"
         abserr(:) = ABS(pot(:,4)-(potgpb(:,0)-potgmb(:,0))/(2._dbl*rrho*grhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,4)+1.e-30)
         CALL output ( iw, grhorefb, pot(:,4), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Second Derivative RaRa"
         abserr(:) = ABS(pot(:,5)-(potpa(:,1)-potma(:,1))/(2._dbl*rrho*rhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,5)+1.e-30)
         CALL output ( iw, rhorefa, pot(:,5), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Second Derivative RaGa"
         abserr(:) = ABS(pot(:,6)-(potgpa(:,1)-potgma(:,1))/(2._dbl*rrho*grhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,6)+1.e-30)
         CALL output ( iw, grhorefa, pot(:,6), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Second Derivative GaGa"
         abserr(:) = ABS(pot(:,7)-(potgpa(:,2)-potgma(:,2))/(2._dbl*rrho*grhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,7)+1.e-30)
         CALL output ( iw, grhorefa, pot(:,7), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Second Derivative RbRb"
         abserr(:) = ABS(pot(:,8)-(potpb(:,3)-potmb(:,3))/(2._dbl*rrho*rhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,8)+1.e-30)
         CALL output ( iw, rhorefb, pot(:,8), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Second Derivative RbGb"
         abserr(:) = ABS(pot(:,9)-(potgpb(:,3)-potgmb(:,3))/(2._dbl*rrho*grhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,9)+1.e-30)
         CALL output ( iw, grhorefb, pot(:,9), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Second Derivative GbGb"
         abserr(:) = ABS(pot(:,10)-(potgpb(:,4)-potgmb(:,4))/(2._dbl*rrho*grhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,10)+1.e-30)
         CALL output ( iw, grhorefb, pot(:,10), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative RaRaRa"
         abserr(:) = ABS(pot(:,11)-(potpa(:,5)-potma(:,5))/(2._dbl*rrho*rhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,11)+1.e-30)
         CALL output ( iw, rhorefa, pot(:,11), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative RaRaGa"
         abserr(:) = ABS(pot(:,12)-(potgpa(:,5)-potgma(:,5))/(2._dbl*rrho*grhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,12)+1.e-30)
         CALL output ( iw, grhorefa, pot(:,12), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative RaGaGa"
         abserr(:) = ABS(pot(:,13)-(potgpa(:,6)-potgma(:,6))/(2._dbl*rrho*grhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,13)+1.e-30)
         CALL output ( iw, grhorefa, pot(:,13), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative GaGaGa"
         abserr(:) = ABS(pot(:,14)-(potgpa(:,7)-potgma(:,7))/(2._dbl*rrho*grhorefa))
         relerr(:) = abserr(:)/ABS(pot(:,14)+1.e-30)
         CALL output ( iw, grhorefa, pot(:,14), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative RbRbRb"
         abserr(:) = ABS(pot(:,15)-(potpb(:,8)-potmb(:,8))/(2._dbl*rrho*rhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,15)+1.e-30)
         CALL output ( iw, rhorefb, pot(:,15), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative RbRbGb"
         abserr(:) = ABS(pot(:,16)-(potgpb(:,8)-potgmb(:,8))/(2._dbl*rrho*grhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,16)+1.e-30)
         CALL output ( iw, grhorefb, pot(:,16), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative RbGbGb"
         abserr(:) = ABS(pot(:,17)-(potgpb(:,9)-potgmb(:,9))/(2._dbl*rrho*grhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,17)+1.e-30)
         CALL output ( iw, grhorefb, pot(:,17), abserr, relerr )
         !
         WRITE (iw,'(/,A)') " Third Derivative GbGbGb"
         abserr(:) = ABS(pot(:,18)-(potgpb(:,10)-potgmb(:,10))/(2._dbl*rrho*grhorefb))
         relerr(:) = abserr(:)/ABS(pot(:,18)+1.e-30)
         CALL output ( iw, grhorefb, pot(:,18), abserr, relerr )
         !
       END IF
       !
     END DO
     !
     DEALLOCATE ( rhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoref" )
     DEALLOCATE ( rhorefa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefa" )
     DEALLOCATE ( rhorefb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhorefb" )
     DEALLOCATE ( rho, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rho" )
     DEALLOCATE ( rhoa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhoa" )
     DEALLOCATE ( rhob, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "rhob" )
     DEALLOCATE ( grhoref, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoref" )
     DEALLOCATE ( grhorefa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefa" )
     DEALLOCATE ( grhorefb, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhorefb" )
     DEALLOCATE ( grho, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grho" )
     DEALLOCATE ( grhoa, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhoa" )
     DEALLOCATE ( grhob, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "grhob" )
     DEALLOCATE ( pot, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "pot" )
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
     DEALLOCATE ( abserr, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "abserr" )
     DEALLOCATE ( relerr, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "relerr" )

  END IF

END SUBROUTINE functional_test

! *****************************************************************************
SUBROUTINE output(iw,rhoref,pot,abserr,relerr)
  INTEGER :: iw
  REAL(dbl), DIMENSION(:) :: rhoref, pot, abserr, relerr

  INTEGER :: ia,ir

  ia=MAXVAL(MAXLOC(abserr))
  ir=MAXVAL(MAXLOC(relerr))
  WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ia),pot(ia),abserr(ia),relerr(ia)
  WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ir),pot(ir),abserr(ir),relerr(ir)

END SUBROUTINE output

! *****************************************************************************

END MODULE functionals_tests

! *****************************************************************************
