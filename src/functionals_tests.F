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
!!   SOURCE
!******************************************************************************

MODULE functionals_tests

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE global_types, ONLY : global_environment_type
  USE pade_functional, ONLY : pade
  USE xalpha_functional, ONLY : xalpha
  USE termination, ONLY : stop_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: functional_test

!!*****
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
!!   SOURCE
!******************************************************************************

SUBROUTINE functional_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  INTEGER, PARAMETER :: n=20
  INTEGER :: iw, ierr, i, nfun, j, ia, ir
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rhoref, rhorefa, rhorefb
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: rho, rhoa, rhob
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: pot, potp, potm
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
     ALLOCATE ( relerr ( n ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "relerr", n )
     ! Test LDA functionals
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
     nfun = 2
     DO i = 1, nfun

       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** X-ALPHA FUNCTIONAL "
         CALL xalpha ( epsilon )
         CALL xalpha ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL xalpha ( rhoref, pot, 3 )
         ! standard tests
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
         ! standard tests
         rho = rhoref + rrho*rhoref
         CALL pade ( rho, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL pade ( rho, potm, 3 )
       END SELECT

       WRITE (iw,'(/,A,/)') " First Derivative "
       abserr(:) = ABS(pot(:,1)-(potp(:,0)-potm(:,0))/(2._dbl*rrho*rhoref))
       relerr(:) = abserr(:)/ABS(pot(:,1))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ia),pot(ia,1),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ir),pot(ir,1),abserr(ir),relerr(ir)
       !
       WRITE (iw,'(/,A,/)') " Second Derivative "
       abserr(:) = ABS(pot(:,2)-(potp(:,1)-potm(:,1))/(2._dbl*rrho*rhoref))
       relerr(:) = abserr(:)/ABS(pot(:,2))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ia),pot(ia,1),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ir),pot(ir,1),abserr(ir),relerr(ir)
       !
       WRITE (iw,'(/,A,/)') " Third Derivative "
       abserr(:) = ABS(pot(:,3)-(potp(:,2)-potm(:,2))/(2._dbl*rrho*rhoref))
       relerr(:) = abserr(:)/ABS(pot(:,3))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ia),pot(ia,1),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhoref(ir),pot(ir,1),abserr(ir),relerr(ir)
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

     ! Test LSD functionals
     WRITE (iw,'(/,A,/)') " *** TEST OF LSD FUNCTIONALS *** "
     WRITE (iw,'(A,E10.3)') " Relative step size for finite difference ",rrho
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
     !
     nfun = 2
     DO i = 1, nfun
       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** X-ALPHA FUNCTIONAL "
         CALL xalpha ( epsilon )
         CALL xalpha ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL xalpha ( rhorefa, rhorefb, pot, 3 )
         ! standard tests
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
         CALL pade ( rhorefa, rhorefb, pot, 3 )
         ! standard tests
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
       END SELECT
       WRITE (iw,'(/,A,/)') " First Derivative A"
       abserr(:) = ABS(pot(:,1)-(potpa(:,0)-potma(:,0))/(2._dbl*rrho*rhorefa))
       relerr(:) = abserr(:)/ABS(pot(:,1))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefa(ia),pot(ia,1),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefa(ir),pot(ir,1),abserr(ir),relerr(ir)
       !
       WRITE (iw,'(/,A,/)') " First Derivative B"
       abserr(:) = ABS(pot(:,2)-(potpb(:,0)-potmb(:,0))/(2._dbl*rrho*rhorefb))
       relerr(:) = abserr(:)/ABS(pot(:,2))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefb(ia),pot(ia,2),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefb(ir),pot(ir,2),abserr(ir),relerr(ir)
       !
       WRITE (iw,'(/,A,/)') " Second Derivative AA"
       abserr(:) = ABS(pot(:,3)-(potpa(:,1)-potma(:,1))/(2._dbl*rrho*rhorefa))
       relerr(:) = abserr(:)/ABS(pot(:,3))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefa(ia),pot(ia,3),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefa(ir),pot(ir,3),abserr(ir),relerr(ir)
       !
       IF ( crossterm ) THEN
          WRITE (iw,'(/,A,/)') " Second Derivative AB"
          abserr(:) = ABS(pot(:,4)-(potpa(:,2)-potma(:,2))/(2._dbl*rrho*rhorefa))
          relerr(:) = abserr(:)/ABS(pot(:,4))
          ia=MAXVAL(MAXLOC(abserr))
          ir=MAXVAL(MAXLOC(relerr))
          WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
                rhorefa(ia),pot(ia,4),abserr(ia),relerr(ia)
          WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
                rhorefa(ir),pot(ir,4),abserr(ir),relerr(ir)
       END IF
       !
       WRITE (iw,'(/,A,/)') " Second Derivative BB"
       abserr(:) = ABS(pot(:,5)-(potpb(:,2)-potmb(:,2))/(2._dbl*rrho*rhorefb))
       relerr(:) = abserr(:)/ABS(pot(:,5))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefb(ia),pot(ia,4),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefb(ir),pot(ir,4),abserr(ir),relerr(ir)
       !
       WRITE (iw,'(/,A,/)') " Third Derivative AAA"
       abserr(:) = ABS(pot(:,6)-(potpa(:,3)-potma(:,3))/(2._dbl*rrho*rhorefa))
       relerr(:) = abserr(:)/ABS(pot(:,6))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefa(ia),pot(ia,6),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefa(ir),pot(ir,6),abserr(ir),relerr(ir)
       !
       IF ( crossterm ) THEN
          WRITE (iw,'(/,A,/)') " Third Derivative AAB"
          abserr(:) = ABS(pot(:,7)-(potpb(:,3)-potmb(:,3))/(2._dbl*rrho*rhorefb))
          relerr(:) = abserr(:)/ABS(pot(:,7))
          ia=MAXVAL(MAXLOC(abserr))
          ir=MAXVAL(MAXLOC(relerr))
          WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
                rhorefb(ia),pot(ia,7),abserr(ia),relerr(ia)
          WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
                rhorefb(ir),pot(ir,7),abserr(ir),relerr(ir)
          !
          WRITE (iw,'(/,A,/)') " Third Derivative ABB"
          abserr(:) = ABS(pot(:,8)-(potpa(:,5)-potma(:,5))/(2._dbl*rrho*rhorefa))
          relerr(:) = abserr(:)/ABS(pot(:,8))
          ia=MAXVAL(MAXLOC(abserr))
          ir=MAXVAL(MAXLOC(relerr))
          WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
                rhorefa(ia),pot(ia,8),abserr(ia),relerr(ia)
          WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
                rhorefa(ir),pot(ir,8),abserr(ir),relerr(ir)
       END IF
       !
       WRITE (iw,'(/,A,/)') " Third Derivative BBB"
       abserr(:) = ABS(pot(:,9)-(potpb(:,5)-potmb(:,5))/(2._dbl*rrho*rhorefb))
       relerr(:) = abserr(:)/ABS(pot(:,9))
       ia=MAXVAL(MAXLOC(abserr))
       ir=MAXVAL(MAXLOC(relerr))
       WRITE (iw,'("AERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefb(ia),pot(ia,9),abserr(ia),relerr(ia)
       WRITE (iw,'("RERR:",G14.6,5X,G14.8,5X,G16.8,5X,G16.8)') &
             rhorefb(ir),pot(ir,9),abserr(ir),relerr(ir)
     END DO
     !
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

     DEALLOCATE ( abserr, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "abserr" )
     DEALLOCATE ( relerr, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "functional_test", "relerr" )

  END IF

END SUBROUTINE functional_test

!!*****
! *****************************************************************************

END MODULE functionals_tests

! *****************************************************************************
