!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/atom_module [1.0] *
!!
!!   NAME
!!     atom_module
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH, 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE atom_module

  USE global_types, ONLY : global_environment_type
  USE header, ONLY : atom_header
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero, pi
  USE timings, ONLY : timeset, timestop
!deb
  USE radial_grids, ONLY : radial_grid_type, init_radial_grid
  USE radial_operators, ONLY : rad_kinetic
  USE radial_poisson, ONLY : rg_poisson
  USE radial_util, ONLY : radial_integrate, radial_diff
!deb

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atomic_code

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!****** atom_module/atomic_code [1.0] *
!!
!!   NAME
!!     atomic_code
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH, 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE atomic_code ( globenv )

    IMPLICIT NONE

!Arguments
    TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

!Locals
    REAL ( dbl ) :: xx
    INTEGER :: handle

!deb
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: f,df,fp,p
  REAL ( dbl ) :: h, x, res, erf, z, b(3), rmax
  Integer :: n, ierr, i, l
  TYPE ( radial_grid_type ) :: rg
!deb
!------------------------------------------------------------------------------

    CALL timeset ( 'ATOMIC_CODE', 'I', ' ', handle )

    IF ( globenv % ionode ) CALL atom_header ( globenv % scr )

!debdebdebdebdeb
    n = 100
    z=1._dbl
    b(1)=1._dbl
    b(2)=4._dbl
    rmax=10._dbl
    CALL init_radial_grid ( rg, "LINEAR", n, z, b, rmax )
!   CALL init_radial_grid ( rg, "LOG", n, z, b, rmax )
!   CALL init_radial_grid ( rg, "POWER", n, z, b, rmax )
!   CALL init_radial_grid ( rg, "RATIONAL", n, z, b, rmax )
!   CALL init_radial_grid ( rg, "ALGEBRAIC SQUARE ROOT", n, z, b, rmax )
    allocate ( f(0:n+1) )
    allocate ( p(0:n+1) )
    allocate ( fp(0:n+1) )
    allocate ( df(0:n+1) )
    do i = 0, n+1
      x = rg % r ( i )
      f ( i ) = exp ( - x*x )
      p ( i ) = x * f ( i )
      fp ( i ) = -2._dbl*x*f(i)
    enddo
    l = 0
    call rad_kinetic ( p, p, l, rg, res )
    res = res * 4._dbl * pi
    write ( 6,*) ' KINETIC ',res
    write ( 6,*) ' KINETIC ',1.5_dbl * ( 0.5_dbl*pi )**1.5_dbl
    write ( 6,*) ' KINETIC ',res - 1.5_dbl * ( 0.5_dbl*pi )**1.5_dbl
    call radial_diff ( n+2, f(0:n+1),df(0:n+1),1,rg,ierr )
    call radial_integrate ( n+2,f(0:n+1),rg,res,ierr )
    write ( 6,*) ' INTEGRAL ',res
    write ( 6,*) ' INTEGRAL ',0.5_dbl * sqrt ( pi ) !* erf (rmax)
    write ( 6,*) ' INTEGRAL ',res - 0.5_dbl * sqrt ( pi ) !* erf (rmax)
    write ( 6,*) ' DIFF     ',sum ( abs ( df - fp ) ) / real ( n )
    write ( 6,*) ' DIFF     ',minval ( abs ( df - fp ) )
    write ( 6,*) ' DIFF     ',maxval ( abs ( df - fp ) )

    do i = 0, n+1
      x = rg % r ( i )
      f ( i ) = exp ( - x*x )
      fp ( i ) = (4._dbl*x*x-2._dbl)*f(i)
    enddo
    call radial_diff ( n+2, f(0:n+1),df(0:n+1),2,rg,ierr )
    write ( 6,*) ' d2DIFF     ',sum ( abs ( df - fp ) ) / real ( n )
    write ( 6,*) ' d2DIFF     ',minval ( abs ( df - fp ) )
    write ( 6,*) ' d2DIFF     ',maxval ( abs ( df - fp ) )
    do i = 0, n+1
      x = rg % r ( i )
      f ( i ) = -z* (pi)**(-1.5) * x*x*exp ( - x*x ) * 4._dbl*pi
      p ( i ) = -z*erf ( x )
    enddo
!   call rg_poisson ( f, fp, l, rg, "NUMEROV" )
    call rg_poisson ( f, fp, l, rg, "DIFF5P" )
    write ( 6,*) ' POISSON  ',sum ( abs ( fp - p ) ) / real ( n )
    write ( 6,*) ' POISSON  ',minval ( abs ( fp - p ) )
    write ( 6,*) ' POISSON  ',maxval ( abs ( fp - p ) )
    do i = 0, n+1
      write(99,*) fp(i),p(i)
    enddo

    deallocate ( f,df,fp,p )
!debdebdebdebdeb

    CALL timestop ( zero, handle )

  END SUBROUTINE atomic_code

!!*****
!******************************************************************************

END MODULE atom_module

!******************************************************************************
