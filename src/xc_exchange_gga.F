!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_exchange_gga [1.0] *
!!
!!   NAME
!!     xc_exchange_gga
!!
!!   FUNCTION
!!     Calculate several different exchange energy functionals
!!     with a GGA form 
!!
!!   AUTHOR
!!     JGH (27.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!*** *************************************************************************

MODULE xc_exchange_gga
! *****************************************************************************
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals_utilities,        ONLY: calc_rho13,&
                                             calc_wave_vector,&
                                             set_util,&
                                             setup_calculation

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: b88, x_pw86, x_pw91, pbex, revpbex, optx

  REAL(dbl) :: cx, flda, flsd, sfac, t13
  REAL(dbl) :: fact, tact
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE b88
    MODULE PROCEDURE b88_u, b88_p, xgga_init, b88_info
  END INTERFACE
  INTERFACE x_pw86
    MODULE PROCEDURE pw86_u, pw86_p, xgga_init, pw86_info
  END INTERFACE
  INTERFACE x_pw91
    MODULE PROCEDURE pw91_u, pw91_p, xgga_init, pw91_info
  END INTERFACE
  INTERFACE pbex
    MODULE PROCEDURE pbex_u, pbex_p, xgga_init, pbex_info
  END INTERFACE
  INTERFACE revpbex
    MODULE PROCEDURE revpbex_u, revpbex_p, xgga_init, revpbex_info
  END INTERFACE
  INTERFACE optx
    MODULE PROCEDURE optx_u, optx_p, xgga_init, optx_info
  END INTERFACE

! *****************************************************************************

CONTAINS



! *****************************************************************************

  SUBROUTINE xgga_init ( cutoff )

    REAL(dbl), INTENT(IN)                    :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cx = -0.75_dbl*(3._dbl/pi)**f13
    t13 = 2._dbl**f13
    flda = cx
    flsd = cx * t13

    sfac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)

  END SUBROUTINE xgga_init

! *****************************************************************************

  SUBROUTINE b88_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "A. Becke, Phys. Rev. A 38, 3098 (1988)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Becke 1988 Exchange Functional"
    END IF

  END SUBROUTINE b88_info

  SUBROUTINE pw86_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800 (1986)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew-Wang 1986 Functional (exchange energy)"
    END IF

  END SUBROUTINE pw86_info

  SUBROUTINE pw91_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew et al., Phys. Rev. B, 46, 6671 (1992)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew-Wang 1991 Functional (exchange energy)"
    END IF

  END SUBROUTINE pw91_info

  SUBROUTINE pbex_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew, K. Burke, M Ernzerhof, Phys. Rev. Lett, 77, 3865 (1996)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "PBE exchange energy functional"
    END IF

  END SUBROUTINE pbex_info

  SUBROUTINE revpbex_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "Y. Zang et al., PRL, 80, 890 (1998) (Revised PBEX)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Revised PBEX by Zang et al."
    END IF

  END SUBROUTINE revpbex_info

  SUBROUTINE optx_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "Wee-Meng Hoe, A.J. Cohen, N.C. Handy, CPL, 341, 319 (2001)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "OPTX exchange energy functional"
    END IF

  END SUBROUTINE optx_info

! *****************************************************************************

  SUBROUTINE b88_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_u ( derivative_set, 1, rho, grho, pot, order )
  END SUBROUTINE b88_u
  SUBROUTINE b88_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_p ( derivative_set, 1, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE b88_p

  SUBROUTINE pw86_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_u ( derivative_set, 2, rho, grho, pot, order )
  END SUBROUTINE pw86_u
  SUBROUTINE pw86_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_p ( derivative_set, 2, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pw86_p

  SUBROUTINE pw91_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_u ( derivative_set, 3, rho, grho, pot, order )
  END SUBROUTINE pw91_u
  SUBROUTINE pw91_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_p ( derivative_set, 3, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pw91_p

  SUBROUTINE pbex_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_u ( derivative_set, 4, rho, grho, pot, order )
  END SUBROUTINE pbex_u
  SUBROUTINE pbex_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_p ( derivative_set, 4, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pbex_p

  SUBROUTINE revpbex_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_u ( derivative_set, 5, rho, grho, pot, order )
  END SUBROUTINE revpbex_u
  SUBROUTINE revpbex_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_p ( derivative_set, 5, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE revpbex_p

  SUBROUTINE optx_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_u ( derivative_set, 6, rho, grho, pot, order )
  END SUBROUTINE optx_u
  SUBROUTINE optx_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL x_p ( derivative_set, 6, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE optx_p

! *****************************************************************************

  SUBROUTINE x_u ( derivative_set, functional, rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: functional
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: error, m, m1, m2, mp, mt, &
                                                mx(0:3,2), n
    LOGICAL                                  :: calc(0:4)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: r13, s
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: fs

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,110)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "x_u","Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_u", "r13", n )
    ALLOCATE ( s(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_u", "s", n )
    ALLOCATE ( fs(n,m+1), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_u", "fs", n*(m+1) )

    fact = flda
    tact = t13
    
    CALL calc_rho13 ( rho, r13 )
    CALL calc_wave_vector ( "u", rho, grho, s )

    SELECT CASE (functional)
    CASE (1)
      CALL efactor_b88 ( s, fs, m )
    CASE (2)
      CALL efactor_pw86 ( s, fs, m )
    CASE (3)
      CALL efactor_pw91 ( s, fs, m )
    CASE (4)
      CALL efactor_pbex ( s, fs, m, 1 )
    CASE (5)
      CALL efactor_pbex ( s, fs, m, 2 )
    CASE (6)
      CALL efactor_optx ( s, fs, m )
    CASE DEFAULT
       CALL stop_program ( "x_u", "Unknown functional" )
    END SELECT

    IF ( calc(0) ) THEN
       CALL x_p_0 ( rho, r13, fs, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL x_p_1 ( rho, grho, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL x_p_2 ( rho, grho, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL x_p_3 ( rho, grho, r13, s, fs, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_u", "r13" )
    DEALLOCATE ( s, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_u", "s" )
    DEALLOCATE ( fs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_u", "fs" )

  END SUBROUTINE x_u

! *****************************************************************************

  SUBROUTINE x_p ( derivative_set, functional, rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: functional
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: error, m, m1, m2, mp, mt, &
                                                mx(0:3,2), n
    LOGICAL                                  :: calc(0:4)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: r13, s
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: fs

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,210)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "x_p", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_p", "r13", n )
    ALLOCATE ( s(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_p", "s", n )
    ALLOCATE ( fs(n,m+1), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_p", "fs", n*(m+1) )

    fact = flsd
    tact = 1._dbl

    CALL calc_rho13 ( rhoa, r13 )
    CALL calc_wave_vector ( "p", rhoa, grhoa, s )
    SELECT CASE (functional)
    CASE (1)
      CALL efactor_b88 ( s, fs, m )
    CASE (2)
      CALL efactor_pw86 ( s, fs, m )
    CASE (3)
      CALL efactor_pw91 ( s, fs, m )
    CASE (4)
      tact = t13
      CALL efactor_pbex ( s, fs, m, 1 )
      tact = 1._dbl
    CASE (5)
      tact = t13
      CALL efactor_pbex ( s, fs, m, 2 )
      tact = 1._dbl
    CASE (6)
      CALL efactor_optx ( s, fs, m )
    CASE DEFAULT
       CALL stop_program ( "x_p", "Unknown functional" )
    END SELECT

    IF ( calc(0) ) THEN
       CALL x_p_0 ( rhoa, r13, fs, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL x_p_1 ( rhoa, grhoa, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL x_p_2 ( rhoa, grhoa, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL x_p_3 ( rhoa, grhoa, r13, s, fs, pot(:,m1:m2) )
    END IF

    CALL calc_rho13 ( rhob, r13 )
    CALL calc_wave_vector ( "p", rhob, grhob, s )
    SELECT CASE (functional)
    CASE (1)
      CALL efactor_b88 ( s, fs, m )
    CASE (2)
      CALL efactor_pw86 ( s, fs, m )
    CASE (3)
      CALL efactor_pw91 ( s, fs, m )
    CASE (4)
      tact = t13
      CALL efactor_pbex ( s, fs, m, 1 )
      tact = 1._dbl
    CASE (5)
      tact = t13
      CALL efactor_pbex ( s, fs, m, 2 )
      tact = 1._dbl
    CASE (6)
      CALL efactor_optx ( s, fs, m )
    CASE DEFAULT
       CALL stop_program ( "x_p", "Unknown functional" )
    END SELECT

    IF ( calc(0) ) THEN
       CALL x_p_0 ( rhob, r13, fs, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL x_p_1 ( rhob, grhob, r13, s, fs, pot(:,m1+2:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL x_p_2 ( rhob, grhob, r13, s, fs, pot(:,m1+3:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL x_p_3 ( rhob, grhob, r13, s, fs, pot(:,m1+4:m2) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_p", "r13" )
    DEALLOCATE ( s, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_p", "s" )
    DEALLOCATE ( fs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "x_p", "fs" )

  END SUBROUTINE x_p

! *****************************************************************************

  SUBROUTINE x_p_0 ( rho, r13, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip, npoints

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "x_p_0", "Inconsistent array sizes" )

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + fact*r13(ip)*rho(ip) * fs(ip,1)
      END IF

    END DO

  END SUBROUTINE x_p_0

! *****************************************************************************

  SUBROUTINE x_p_1 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: a0, a1, sx, sy

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "x_p_1", "Inconsistent array sizes" )

!$omp parallel do private(ip,a0,a1,sx,sy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         pot(ip,1) = pot(ip,1) + a1*fs(ip,1) + a0*fs(ip,2)*sx
         pot(ip,2) = pot(ip,2) + a0*fs(ip,2)*sy

      END IF

    END DO

  END SUBROUTINE x_p_1

! *****************************************************************************

  SUBROUTINE x_p_2 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: a0, a1, a2, sx, sxx, sxy, sy

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "x_p_2", "Inconsistent array sizes" )

!$omp parallel do private(ip,a0,a1,a2,sx,sy,sxx,sxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         a2 = f13*f43*fact/(r13(ip)*r13(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         pot(ip,1) = pot(ip,1) + a2*fs(ip,1) +2._dbl*a1*fs(ip,2)*sx + &
                     a0*fs(ip,3)*sx*sx + a0*fs(ip,2)*sxx
         pot(ip,2) = pot(ip,2) + a1*fs(ip,2)*sy + a0*fs(ip,3)*sx*sy + &
                     a0*fs(ip,2)*sxy
         pot(ip,3) = pot(ip,3) + a0*fs(ip,3)*sy*sy

      END IF

    END DO

  END SUBROUTINE x_p_2

! *****************************************************************************

  SUBROUTINE x_p_3 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: a0, a1, a2, a3, sx, sxx, &
                                                sxxx, sxxy, sxy, sy

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "x_u_3", "Inconsistent array sizes" )

!$omp parallel do private(ip,a0,a1,a2,a3,sx,sy,sxx,sxy,sxxx,sxxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         a2 = f13*f43*fact/(r13(ip)*r13(ip))
         a3 = -f23*f13*f43*fact/(r13(ip)*r13(ip)*rho(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         sxxx= -280._dbl/27._dbl*s(ip)/(rho(ip)*rho(ip)*rho(ip))
         sxxy= 28._dbl/9._dbl*sfac*tact/(r13(ip)*rho(ip)*rho(ip)*rho(ip))
         pot(ip,1) = pot(ip,1) + a3*fs(ip,1) + 3._dbl*a2*fs(ip,2)*sx + &
                     3._dbl*a1*fs(ip,3)*sx*sx + 3._dbl*a1*fs(ip,2)*sxx + &
                     a0*fs(ip,4)*sx*sx*sx + 3._dbl*a0*fs(ip,3)*sx*sxx + &
                     a0*fs(ip,2)*sxxx
         pot(ip,2) = pot(ip,2) + a2*fs(ip,2)*sy + 2._dbl*a1*fs(ip,3)*sx*sy + &
                     2._dbl*a1*fs(ip,2)*sxy + a0*fs(ip,4)*sx*sx*sy + &
                     2._dbl*a0*fs(ip,3)*sx*sxy + a0*fs(ip,3)*sxx*sy + &
                     a0*fs(ip,2)*sxxy
         pot(ip,3) = pot(ip,3) + a1*fs(ip,3)*sy*sy + a0*fs(ip,4)*sx*sy*sy + &
                     2._dbl*a0*fs(ip,3)*sxy*sy
         pot(ip,4) = pot(ip,4) + a0*fs(ip,4)*sy*sy*sy

      END IF

    END DO

  END SUBROUTINE x_p_3

! *****************************************************************************
! Enhancement Factors
! *****************************************************************************
  SUBROUTINE efactor_b88 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl) :: as, asp, beta, bs, f0, p, q, sas, sbs, sbs3, t1, t10, t13, &
      t15, t16, t19, t2, t22, t24, t25, t32, t34, t36, t39, t4, t40, t41, &
      t44, t48, t49, t5, t6, t65, t8, t87, t9, x, ys

    beta = 0.0042_dbl
    f0 = 1._dbl/sfac
    p = -beta/flsd
    q = 6._dbl*beta
!$omp parallel do private(ip,x,bs,sbs,as,sas,ys,asp,sbs3) &
!$omp   private(t1,t2,t4,t5,t6,t8,t9,t10,t13,t15,t16,t19,t22,t24,t25,t32,t34) &
!$omp   private(t36,t39,t40,t41,t44,t48,t49,t65,t87)
    DO ip=1,SIZE(s)
      x = s(ip) * f0
      bs = beta*x
      sbs = SQRT(x*x+1._dbl)
      as = LOG(x+sbs)
      sas = x*as
      ys =  1._dbl/(1._dbl+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + p*x*x*ys
      CASE (1)
         asp = as + x / sbs
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = ( 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
      CASE (2)
         asp = as + x / sbs
         sbs3 = 1._dbl/(sbs*sbs*sbs)
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = ( 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
         fs(ip,3) = -f0*f0*p*ys**3*sbs3 * ( q*x*x*x*x*(q*sas + 5._dbl &
                    - 2._dbl*q*sbs) + 2._dbl*( x*x*(q*q*sas & 
                    + 3._dbl*q - sbs) - sbs))
      CASE (3)
         asp = as + x / sbs
         sbs3 = 1._dbl/(sbs*sbs*sbs)
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = ( 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
         fs(ip,3) = -f0*f0*p*ys**3*sbs3 * ( q*x*x*x*x*(q*sas + 5._dbl &
                    - 2._dbl*q*sbs) + 2._dbl*( x*x*(q*q*sas & 
                    + 3._dbl*q - sbs) - sbs))
      t1 = q*x
      t2 = x**2
      t4 = SQRT(1+t2)
      t5 = x+t4
      t6 = LOG(t5)
      t8 = 1+t1*t6
      t9 = t8**2
      t10 = 1/t9
      t13 = 1/t4
      t15 = 1+t13*x
      t16 = 1/t5
      t19 = q*t6+t1*t15*t16
      t22 = p*x
      t24 = 1/t9/t8
      t25 = t19**2
      t32 = t4**2
      t34 = 1/t32/t4
      t36 = -t34*t2+t13
      t39 = t15**2
      t40 = t5**2
      t41 = 1/t40
      t44 = 2*q*t15*t16+t1*t36*t16-t1*t39*t41
      t48 = p*t2
      t49 = t9**2
      t65 = t32**2
      t87 = -6*p*t10*t19+12*t22*t24*t25-6*t22*t10*t44-6*t48/t49*t25*t19+ &
      6*t48*t24*t19*t44-t48*t10*(3*q*t36*t16-3*q*t39*t41+3*t1*(1/t65/t4* &
      t2*x-t34*x)*t16-3*t1*t36*t41*t15+2*t1*t39*t15/t40/t5)

         fs(ip,4) = t87
         fs(ip,4) = f0*f0*f0*fs(ip,4)

      CASE DEFAULT
         CALL stop_program("efactor_b88","Illegal order")
      END SELECT
    END DO


  END SUBROUTINE efactor_b88
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: f15, p0, p1, p15, s2, s4, s6, &
                                                t1, t10, t12, t13, t14, t19, &
                                                t2, t25, t3, t8, t9

    t1 = 1.296_dbl
    t2 = 14._dbl
    t3 = 0.2_dbl
    f15 = 1._dbl/15._dbl
!$omp parallel do private(ip,s2,s4,s6,p0,p1,p15,f15)&
!$omp  private( t1, t2, t3, t8, t9, t10, t12, t13, t14, t19, t25)
    DO ip=1,SIZE(s)
      s2 = s(ip)*s(ip)
      s4 = s2*s2
      s6 = s2*s4
      SELECT CASE (m)
      CASE (0)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         fs(ip,1) = p0**f15
      CASE (1)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
      CASE (2)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14._dbl/225._dbl/t13/p0*t25 + &
                     1._dbl/t13*(2._dbl*t1+12*t2*s2+30._dbl*t3*s4)/15._dbl
      CASE (3)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14._dbl/225._dbl/t13/p0*t25 + &
                     1._dbl/t13*(2._dbl*t1+12*t2*s2+30._dbl*t3*s4)/15._dbl
         t8 = p0**2; t9 = p0**f15; t14 = p0/t9; t19 = s2*s(ip)
         fs(ip,4) = 406._dbl/3375._dbl/t14/t8*p1*p1*p1-14._dbl/&
                    75._dbl/t14/p0*p1*(2*t1+12*t2*s2+30*t3*s4)+&
                    1/t14*(24*t2*s(ip)+120*t3*t19)*f15
      CASE DEFAULT
         CALL stop_program("efactor_ol1","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_pw86
! *****************************************************************************
  SUBROUTINE efactor_optx ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    REAL(dbl), PARAMETER                     :: a1 = 1.05151_dbl, &
                                                a2 = 1.43169_dbl, &
                                                gamma = 0.006_dbl

    INTEGER                                  :: ip
    REAL(dbl)                                :: a, b, f0, x, y

    f0 = 1._dbl/sfac
    b = -a2/flsd
!$omp parallel do private(ip,x,a,y)
    DO ip=1,SIZE(s)
      x = s(ip) * f0
      a = gamma*x*x
      y = 1._dbl / (1._dbl + a)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = a1 + b * a*a*y*y
      CASE (1)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4._dbl*b*f0*a*gamma*x * y*y*y
      CASE (2)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4._dbl*b*f0*a*gamma*x * y*y*y
         fs(ip,3) = -12._dbl * b*f0*f0 * gamma*a*(a-1._dbl) * y*y*y*y
      CASE (3)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4._dbl*b*f0*a*gamma*x * y*y*y
         fs(ip,3) = -12._dbl * b*f0*f0 * gamma*a*(a-1._dbl) * y*y*y*y
         fs(ip,4) = 24._dbl*b*f0*f0*f0* gamma*gamma*x * &
                    (1._dbl-5._dbl*a+2._dbl*a*a) * y*y*y*y*y
      CASE DEFAULT
         CALL stop_program("efactor_optx","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_optx
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: a, a1, a2, a3, a4, b1, o, s1, &
                                                s2, s3, x

    o=1._dbl
    a1 = 0.19645_dbl
    a2 = 0.2743_dbl
    a3 = 0.1508_dbl
    a4 = 100._dbl
    a  = 0.8145161_dbl
    b1 = 0.004_dbl
    IF ( m >= 0 ) THEN
!$omp parallel do private(ip,x)
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,1) = (o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)&
                   /(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+b1*x**4)
      END DO
    END IF
    IF ( m >= 1 ) THEN
!$omp parallel do private(ip,x)
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,2) =&
           (a1*LOG(a*x+SQRT(o+a**2*x**2))+a1*x*a/SQRT(o+a**2*x**2)+2*a3*       &
           a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(o+a1*x*LOG(a*x+SQRT(o  &
           +a**2*x**2))+b1*x**4)-(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+a2-a3*exp  &
           (-a4*x**2)*x**2)/(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+b1*x**4)**2*(a1 &
           *LOG(a*x+SQRT(o+a**2*x**2))+a1*x*a/SQRT(o+a**2*x**2)+4*b1*x**3)
      END DO
    END IF
    IF ( m >= 2 ) THEN
!$omp parallel do private(ip,x,s1)
      DO ip=1,SIZE(s)
        x = s(ip)
        s1 = (2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+  &
             10*a3*a4*x**2*EXP(-a4*x**2)-4*a3*a4**2*x**4*EXP(-a4*x**2)     &
             -2*a3*EXP(-a4*x**2))/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+      &
             b1*x**4)-2*(a1*LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2) &
             +2*a3*a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*    &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**2*(a1*LOG(a*x+SQRT(1+a**2*x**2))&
             +a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)
        fs(ip,3) =&
             s1+2*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)&
             /(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*LOG(a*x+     &
             SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**2-       &
             (1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)     &
             /(1+a1*x*LOG(a*x +SQRT(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/        &
             SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+12*b1*x**2)
      END DO
    END IF
    IF ( m >= 3 ) THEN
!$omp parallel do private(ip,x,s1,s2,s3)
      DO ip=1,SIZE(s)
        x = s(ip)
        s2 = (-4*a1*a**3/SQRT(1+a**2*x**2)**3*x+3*a1*x**3*a**5/         &
             SQRT(1+a**2*x**2)**5+24*a3*a4*x*EXP(-a4*x**2)-36*a3*a4**2* &
             x**3*EXP(-a4*x**2)+8*a3*a4**3*x**5*EXP(-a4*x**2))/(1+a1*x* &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)
        s3 = -3*(2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*     &
             x**2)**3+10*a3*a4*x**2*EXP(-a4*x**2)-4*a3*a4**2*x**4*      &
             EXP(-a4*x**2)-2*a3*EXP(-a4*x**2))/(1+a1*x*LOG(a*x+         &
             SQRT(1+a**2*x**2))+b1*x**4)**2*(a1*LOG(a*x+SQRT(1+a**2*    &
             x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)+6*(a1*LOG(a*x+  &
             SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+2*a3*a4*x**3*  &
             EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*LOG(a*x+       &
             SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*LOG(a*x+SQRT(1+a**2*    &
             x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**2
        s1 = s2+s3
        s2 = s1-3*(a1*LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2) &
             +2*a3*a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*   &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/SQRT(1+a**2*  &
             x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+12*b1*x**2)
        s3 = s2-6*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)* &
             x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**4*(a1*    &
             LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**3
        fs(ip,4) =&
             s3+6*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)  &
             *x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*   &
             LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*    &
             x**3)*(2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*    &
             x**2)**3+12*b1*x**2)-(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-  &
             a3*EXP(-a4*x**2)*x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+   &
             b1*x**4)**2*(-4*a1*a**3/SQRT(1+a**2*x**2)**3*x+3*a1*x**3*    &
             a**5/SQRT(1+a**2*x**2)**5+24*b1*x)
      END DO
    END IF

  END SUBROUTINE efactor_pw91
! *****************************************************************************
  SUBROUTINE efactor_pbex ( s, fs, m, pset )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m, pset

    REAL(dbl), PARAMETER                     :: kappa1 = 0.804_dbl, &
                                                kappa2 = 1.245_dbl, &
                                                mu = 0.2195149727645171_dbl

    INTEGER                                  :: ip
    REAL(dbl)                                :: f0, mk, x, x2, y

    IF( pset==1 ) mk=mu/kappa1
    IF( pset==2 ) mk=mu/kappa2

    f0 = 1._dbl/tact
!$omp parallel do private(ip,x,x2,y)
    DO ip=1,SIZE(s)
      x = s(ip)*f0
      x2 = x*x
      y = 1._dbl/(1._dbl+mk*x2)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + mu*x2*y
      CASE (1)
         fs(ip,1) = 1._dbl + mu*x2*y
         fs(ip,2) = 2._dbl*mu*x*y*y*f0
      CASE (2)
         fs(ip,1) = 1._dbl + mu*x2*y
         fs(ip,2) = 2._dbl*mu*x*y*y*f0
         fs(ip,3) = -2._dbl*mu*(3._dbl*mk*x2-1._dbl)*y*y*y*f0*f0
      CASE (3)
         fs(ip,1) = 1._dbl + mu*x2*y
         fs(ip,2) = 2._dbl*mu*x*y*y*f0
         fs(ip,3) = -2._dbl*mu*(3._dbl*mk*x2-1._dbl)*y*y*y*f0*f0
         fs(ip,4) = 24._dbl*mu*mk*x*(mk*x2-1._dbl)*y*y*y*y*f0*f0*f0
      CASE DEFAULT
         CALL stop_program("efactor_pbex","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_pbex

! *****************************************************************************

END MODULE xc_exchange_gga

! *****************************************************************************
