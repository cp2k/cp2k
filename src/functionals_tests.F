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
  USE exchange_gga_functional, ONLY : b88, x_pw86, x_pw91, pbex, revpbex, optx
  USE ke_gga_functional, ONLY : ol1, ol2, llp, ke_pw86, ke_pw91, lc, t92
  USE lyp_functional, ONLY : lyp
  USE pade_functional, ONLY : pade
  USE xalpha_functional, ONLY : xalpha
  USE thomas_fermi_functional, ONLY : thomas_fermi
  USE tfw_functional, ONLY : tfw
  USE termination, ONLY : stop_memory, stop_program
  USE vwn_functional, ONLY : vwn

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
  REAL(dbl), PARAMETER :: rrho=1.E-5_dbl, epsilon=5.E-13_dbl
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: relerr, abserr

  INTEGER :: iw

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

  INTEGER :: i, nfun, j, ia, ir
  REAL(dbl) :: drho
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
     rhoref(1)=1.E-12_dbl
     DO i=2,n
       rhoref(i)=rhoref(i-1)*6._dbl
     END DO
     !
     nfun = 4
     DO i = 1, nfun

       SELECT CASE (i)
       CASE (1)
         WRITE (iw,'(/,A)') " *** X-ALPHA FUNCTIONAL "
         CALL xalpha ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL xalpha ( epsilon )
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
         WRITE (iw,'(/,A)') " *** VWN FUNCTIONAL "
         CALL vwn ( epsilon )
         CALL vwn ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL vwn ( rhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL vwn ( rho, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL vwn ( rho, potm, 3 )
       CASE (4)
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
     rhoref(1)=1.E-12_dbl
     grhoref(1)=1.E-13_dbl
     DO i=2,n
       rhoref(i)=rhoref(i-1)*6._dbl
       grhoref(i)=grhoref(i-1)*8._dbl
     END DO
     !
     nfun = 14
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
       CASE (2)
         WRITE (iw,'(/,A)') " *** Ou-Yang-Levi Functional 1"
         CALL ol1 ( epsilon )
         CALL ol1 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ol1 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL ol1 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL ol1 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL ol1 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL ol1 ( rhoref, grho, potgm, 3 )
       CASE (3)
         WRITE (iw,'(/,A)') " *** Ou-Yang-Levi Functional 2"
         CALL ol2 ( epsilon )
         CALL ol2 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ol2 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL ol2 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL ol2 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL ol2 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL ol2 ( rhoref, grho, potgm, 3 )
       CASE (4)
         WRITE (iw,'(/,A)') " *** Lee-Lee-Parr Functional"
         CALL llp ( epsilon )
         CALL llp ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL llp ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL llp ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL llp ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL llp ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL llp ( rhoref, grho, potgm, 3 )
       CASE (5)
         WRITE (iw,'(/,A)') " *** Kinetic energy Perdew-Wang 86 "
         CALL ke_pw86 ( epsilon )
         CALL ke_pw86 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ke_pw86 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL ke_pw86 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL ke_pw86 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL ke_pw86 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL ke_pw86 ( rhoref, grho, potgm, 3 )
       CASE (6)
         WRITE (iw,'(/,A)') " *** Kinetic energy Perdew-Wang 91 "
         CALL ke_pw91 ( epsilon )
         CALL ke_pw91 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ke_pw91 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL ke_pw91 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL ke_pw91 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL ke_pw91 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL ke_pw91 ( rhoref, grho, potgm, 3 )
       CASE (7)
         WRITE (iw,'(/,A)') " *** Lembarki-Chermette Functional"
         CALL lc ( epsilon )
         CALL lc ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL lc ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL lc ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL lc ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL lc ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL lc ( rhoref, grho, potgm, 3 )
       CASE (8)
         WRITE (iw,'(/,A)') " *** Becke 1988 Functional"
         CALL b88 ( epsilon )
         CALL b88 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL b88 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL b88 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL b88 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL b88 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL b88 ( rhoref, grho, potgm, 3 )
       CASE (9)
         WRITE (iw,'(/,A)') " *** Perdew-Wang 1986 Exchange Functional"
         CALL x_pw86 ( epsilon )
         CALL x_pw86 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL x_pw86 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL x_pw86 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL x_pw86 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL x_pw86 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL x_pw86 ( rhoref, grho, potgm, 3 )
       CASE (10)
         WRITE (iw,'(/,A)') " *** Perdew-Wang 1991 Exchange Functional"
         CALL x_pw91 ( epsilon )
         CALL x_pw91 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL x_pw91 ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL x_pw91 ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL x_pw91 ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL x_pw91 ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL x_pw91 ( rhoref, grho, potgm, 3 )
       CASE (11)
         WRITE (iw,'(/,A)') " *** Perdew-Burke-Ernzerhof Exchange Functional"
         CALL pbex ( epsilon )
         CALL pbex ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL pbex ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL pbex ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL pbex ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL pbex ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL pbex ( rhoref, grho, potgm, 3 )
       CASE (12)
         WRITE (iw,'(/,A)') " *** PBE Exchange Functional (Zang Parameter)"
         CALL revpbex ( epsilon )
         CALL revpbex ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL revpbex ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL revpbex ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL revpbex ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL revpbex ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL revpbex ( rhoref, grho, potgm, 3 )
       CASE (13)
         WRITE (iw,'(/,A)') " *** OPTX (Handy) Exchange Functional"
         CALL optx ( epsilon )
         CALL optx ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL optx ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL optx ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL optx ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL optx ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL optx ( rhoref, grho, potgm, 3 )
       CASE (14)
         WRITE (iw,'(/,A)') " *** Lee-Yang-Parr Correlation Functional"
         CALL lyp ( epsilon )
         CALL lyp ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL lyp ( rhoref, grhoref, pot, 3 )
         rho = rhoref + rrho*rhoref
         CALL lyp ( rho, grhoref, potp, 3 )
         rho = rhoref - rrho*rhoref
         CALL lyp ( rho, grhoref, potm, 3 )
         grho = grhoref + rrho*grhoref
         CALL lyp ( rhoref, grho, potgp, 3 )
         grho = grhoref - rrho*grhoref
         CALL lyp ( rhoref, grho, potgm, 3 )
       END SELECT
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
         CALL pade ( rhorefa, rhorefb, pot, 3 )
         CALL pade ( rhoref, potx(:,1:1), 0 )
         CALL pade ( rhorefa, rhorefa, potx(:,2:2), 0 )
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
     nfun = 13
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
       CASE (2)
         WRITE (iw,'(/,A)') " *** Ou-Yang-Levi Functional 1"
         CALL ol1 ( epsilon )
         CALL ol1 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ol1 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL ol1 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL ol1 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL ol1 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL ol1 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL ol1 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL ol1 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL ol1 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL ol1 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL ol1 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL ol1 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (3)
         WRITE (iw,'(/,A)') " *** Ou-Yang-Levi Functional 2"
         CALL ol2 ( epsilon )
         CALL ol2 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ol2 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL ol2 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL ol2 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL ol2 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL ol2 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL ol2 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL ol2 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL ol2 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL ol2 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL ol2 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL ol2 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (4)
         WRITE (iw,'(/,A)') " *** Lee-Lee-Parr Functional"
         CALL llp ( epsilon )
         CALL llp ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL llp ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL llp ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL llp ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL llp ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL llp ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL llp ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL llp ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL llp ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL llp ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL llp ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL llp ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (5)
         WRITE (iw,'(/,A)') " *** Kinetic energy Perdew-Wang 86 "
         CALL ke_pw86 ( epsilon )
         CALL ke_pw86 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ke_pw86 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL ke_pw86 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL ke_pw86 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL ke_pw86 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL ke_pw86 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL ke_pw86 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL ke_pw86 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL ke_pw86 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL ke_pw86 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL ke_pw86 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL ke_pw86 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (6)
         WRITE (iw,'(/,A)') " *** Kinetic energy Perdew-Wang 91 "
         CALL ke_pw91 ( epsilon )
         CALL ke_pw91 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL ke_pw91 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL ke_pw91 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL ke_pw91 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL ke_pw91 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL ke_pw91 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL ke_pw91 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL ke_pw91 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL ke_pw91 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL ke_pw91 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL ke_pw91 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL ke_pw91 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (7)
         WRITE (iw,'(/,A)') " *** Lembarki-Chermette Functional"
         CALL lc ( epsilon )
         CALL lc ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL lc ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL lc ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL lc ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL lc ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL lc ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL lc ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL lc ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL lc ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL lc ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL lc ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL lc ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (8)
         WRITE (iw,'(/,A)') " *** Becke 1988 Functional"
         CALL b88 ( epsilon )
         CALL b88 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL b88 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL b88 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL b88 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL b88 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL b88 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL b88 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL b88 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL b88 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL b88 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL b88 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL b88 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (9)
         WRITE (iw,'(/,A)') " *** Perdew-Wang 1986 Exchange Functional"
         CALL x_pw86 ( epsilon )
         CALL x_pw86 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL x_pw86 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL x_pw86 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL x_pw86 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL x_pw86 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL x_pw86 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL x_pw86 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL x_pw86 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL x_pw86 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL x_pw86 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL x_pw86 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL x_pw86 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (10)
         WRITE (iw,'(/,A)') " *** Perdew-Wang 1991 Exchange Functional"
         CALL x_pw91 ( epsilon )
         CALL x_pw91 ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL x_pw91 ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL x_pw91 ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL x_pw91 ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL x_pw91 ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL x_pw91 ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL x_pw91 ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL x_pw91 ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL x_pw91 ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL x_pw91 ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL x_pw91 ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL x_pw91 ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (11)
         WRITE (iw,'(/,A)') " *** PBE Exchange Functional"
         CALL pbex ( epsilon )
         CALL pbex ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL pbex ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL pbex ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL pbex ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL pbex ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL pbex ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL pbex ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL pbex ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL pbex ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL pbex ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL pbex ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL pbex ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (12)
         WRITE (iw,'(/,A)') " *** PBE Exchange Functional (Zang Parameter)"
         CALL revpbex ( epsilon )
         CALL revpbex ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL revpbex ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL revpbex ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL revpbex ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL revpbex ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL revpbex ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL revpbex ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL revpbex ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL revpbex ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL revpbex ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL revpbex ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL revpbex ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       CASE (13)
         WRITE (iw,'(/,A)') " *** OPTX (Handy) Exchange Functional"
         CALL optx ( epsilon )
         CALL optx ( reference=fun_string )
         WRITE (iw,'(A80)') ADJUSTR(fun_string)
         CALL optx ( rhoref, grhoref, potx(:,1:1), 0 )
         CALL optx ( rhorefa, rhorefa, grhorefa, grhorefa, potx(:,2:2), 0 )
         CALL optx ( rhorefa, rhorefb, grhorefa, grhorefb, pot, 3 )
         rhoa = rhorefa + rrho*rhorefa
         rhob = rhorefb
         CALL optx ( rhoa, rhob, grhorefa, grhorefb, potpa, 3 )
         rhoa = rhorefa - rrho*rhorefa
         CALL optx ( rhoa, rhob, grhorefa, grhorefb, potma, 3 )
         rhoa = rhorefa
         rhob = rhorefb + rrho*rhorefb
         CALL optx ( rhoa, rhob, grhorefa, grhorefb, potpb, 3 )
         rhob = rhorefb - rrho*rhorefb
         CALL optx ( rhoa, rhob, grhorefa, grhorefb, potmb, 3 )
         grhoa = grhorefa + rrho*grhorefa
         CALL optx ( rhorefa, rhorefb, grhoa, grhorefb, potgpa, 3 )
         grhoa = grhorefa - rrho*grhorefa
         CALL optx ( rhorefa, rhorefb, grhoa, grhorefb, potgma, 3 )
         grhob = grhorefb + rrho*grhorefb
         CALL optx ( rhorefa, rhorefb, grhorefa, grhob, potgpb, 3 )
         grhob = grhorefb - rrho*grhorefb
         CALL optx ( rhorefa, rhorefb, grhorefa, grhob, potgmb, 3 )
         crossterm = .FALSE.
       END SELECT
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
  IF ( relerr(ir) > 1.e-6_dbl ) THEN
    WRITE (iw,'(/,A)') message
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

END MODULE functionals_tests

! *****************************************************************************
