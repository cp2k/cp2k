!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/radial_poisson [1.0] *
!!
!!   NAME
!!     radial_poisson
!!
!!   FUNCTION
!!     Poisson solver for radial grids
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE radial_poisson

  USE kinds, ONLY : dbl

  USE radial_grids, ONLY : radial_grid_type
  USE radial_solver, ONLY : numerov, diff6p
  USE radial_util, ONLY : radial_integrate
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rg_poisson

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!****** radial_poisson/rg_poisson [1.0] *
!!
!!   NAME
!!     rg_poisson
!!
!!   FUNCTION
!!     Poisson solver for radial grids
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     TCH 09-OCT-2001: - changed the way f is passed to the solvers
!!                        f is no longer multiplied by factors of h^2
!!                      - changed the name of diff5p to diff6p
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE rg_poisson ( rho, vh, nu, rg, method )

!Arguments
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( IN ) :: rho
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( OUT ) :: vh
  INTEGER, INTENT ( IN ) :: nu
  TYPE ( radial_grid_type ), INTENT ( IN ) :: rg
  CHARACTER ( LEN=* ), OPTIONAL :: method

!Locals
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: f, g
  REAL ( dbl ) :: corr ( 3 )
  REAL ( dbl ) :: rn, h12, h60, res
  INTEGER :: n, ierr, i

!------------------------------------------------------------------------------

  n = rg % n
  ALLOCATE ( f ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "radial_poisson", "f", n )
  ALLOCATE ( g ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "radial_poisson", "g", n )

  corr = 0._dbl
  g = ( rg % r ** nu ) * rho / ( rg % w ** 2 )
  CALL radial_integrate ( n+2, g, rg % h, res, ierr )
  IF ( ierr /= 0 ) call stop_program ( "radial_poisson", "radial_integrate" )
  corr ( 3 ) = res * rg % w ( n+1 ) / rg % rmax ** nu
  vh ( 0 ) = 0._dbl
  vh ( n+1 ) = corr ( 3 ) / rg % w ( n+1 )

  rn = REAL ( nu, dbl )
  f ( 0 ) = 0._dbl
  g ( 0 ) = 0._dbl
  DO i = 1, n + 1
    f ( i ) = - rg % wdd ( i ) / rg % w ( i ) - rn * ( rn + 1._dbl ) / &
      ( rg % w ( i )**4 * rg % r ( i )**2 )
    g ( i ) = - ( 2._dbl * rn + 1._dbl ) / &
      ( rg % w ( i )**3 * rg % r ( i ) ) * rho ( i )
  END DO

  IF ( PRESENT ( method ) ) THEN
    SELECT CASE ( method )
    CASE DEFAULT
      CALL stop_program("radial_poisson","method <"//TRIM(method)//"> unknown")
    CASE ( "NUMEROV" )
      CALL numerov ( f, g, rg % w, corr )
    CASE ( "DIFF6P" )
      CALL diff6p ( f, g, rg % w )
    END SELECT
  ELSE
!..use Numerov
    CALL numerov ( f, g, rg % w, corr )
  END IF

  vh ( 1:n ) = g ( 1:n )

  DEALLOCATE ( f, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "radial_poisson", "f" )
  DEALLOCATE ( g, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "radial_poisson", "g" )

  END SUBROUTINE rg_poisson

!!*****
!******************************************************************************

END MODULE radial_poisson

!******************************************************************************
