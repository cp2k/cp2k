!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_ke_gga [1.0] *
!!
!!   NAME
!!     xc_ke_gga
!!
!!   FUNCTION
!!     Calculate the several different kinetic energy functionals
!!     with a GGA form 
!!
!!   AUTHOR
!!     JGH (20.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!******************************************************************************

MODULE xc_ke_gga
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

  PUBLIC :: ol1, ol2, llp, ke_pw86, ke_pw91, t92, lc

  REAL(dbl) :: cf, b, flda, flsd, sfac, t13
  REAL(dbl) :: fact, tact
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE ol1
    MODULE PROCEDURE ol1_u, ol1_p, kegga_init, ol1_info
  END INTERFACE
  INTERFACE ol2
    MODULE PROCEDURE ol2_u, ol2_p, kegga_init, ol2_info
  END INTERFACE
  INTERFACE llp
    MODULE PROCEDURE llp_u, llp_p, kegga_init, llp_info
  END INTERFACE
  INTERFACE ke_pw86
    MODULE PROCEDURE pw86_u, pw86_p, kegga_init, pw86_info
  END INTERFACE
  INTERFACE ke_pw91
    MODULE PROCEDURE pw91_u, pw91_p, kegga_init, pw91_info
  END INTERFACE
  INTERFACE lc
    MODULE PROCEDURE lc_u, lc_p, kegga_init, lc_info
  END INTERFACE
  INTERFACE t92
    MODULE PROCEDURE t92_u, t92_p, kegga_init, t92_info
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE kegga_init ( cutoff )

    REAL(dbl), INTENT(IN)                    :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23
    flda = cf
    flsd = flda * 2._dbl**f23
    b = 2._dbl**f43 * (3._dbl*pi*pi)**(-f13)
    sfac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)
    t13 = 2._dbl**f13

  END SUBROUTINE kegga_init

! *****************************************************************************

  SUBROUTINE ol1_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "H. Ou-Yang and M. Levy, "//&
                  "Intl. J. Quant. Chem. 40, 379 (1991); Functional 1"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Ou-Yang-Levy Functional 1"
    END IF

  END SUBROUTINE ol1_info

  SUBROUTINE ol2_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "H. Ou-Yang and M. Levy, "//&
                  "Intl. J. Quant. Chem. 40, 379 (1991); Functional 2"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Ou-Yang-Levy Functional 2"
    END IF

  END SUBROUTINE ol2_info

  SUBROUTINE llp_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "H. Lee, C. Lee, R.G. Parr, Phys. Rev. A, 44, 768 (1991)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Lee-Lee-Parr Functional"
    END IF

  END SUBROUTINE llp_info

  SUBROUTINE pw86_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800 (1986)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew-Wang 1986 Functional (kinetic energy)"
    END IF

  END SUBROUTINE pw86_info

  SUBROUTINE pw91_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew and Y. Wang, Electronic Structure of Solids 91"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew-Wang 1991 Functional (kinetic energy)"
    END IF

  END SUBROUTINE pw91_info

  SUBROUTINE lc_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "A. Lembarki and H. Chermette, Phys. Rev. A, 50, 5328 (1994)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Lembarki-Chermette kinetic energy functional"
    END IF

  END SUBROUTINE lc_info

  SUBROUTINE t92_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "A.J. Thakkar, Phys. Rev. A, 46, 6920 (1992)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Thakkar 1992 Functional"
    END IF

  END SUBROUTINE t92_info

! *****************************************************************************

  SUBROUTINE ol1_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 1, rho, grho, pot, order )
  END SUBROUTINE ol1_u
  SUBROUTINE ol1_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 1, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE ol1_p

  SUBROUTINE ol2_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 2, rho, grho, pot, order )
  END SUBROUTINE ol2_u
  SUBROUTINE ol2_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 2, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE ol2_p

  SUBROUTINE llp_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 3, rho, grho, pot, order )
  END SUBROUTINE llp_u
  SUBROUTINE llp_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 3, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE llp_p

  SUBROUTINE pw86_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 4, rho, grho, pot, order )
  END SUBROUTINE pw86_u
  SUBROUTINE pw86_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 4, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pw86_p

  SUBROUTINE pw91_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 5, rho, grho, pot, order )
  END SUBROUTINE pw91_u
  SUBROUTINE pw91_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 5, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pw91_p

  SUBROUTINE lc_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 6, rho, grho, pot, order )
  END SUBROUTINE lc_u
  SUBROUTINE lc_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 6, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE lc_p

  SUBROUTINE t92_u ( derivative_set, rho, grho, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_u ( derivative_set, 7, rho, grho, pot, order )
  END SUBROUTINE t92_u
  SUBROUTINE t92_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CALL kex_p ( derivative_set, 7, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE t92_p

! *****************************************************************************

  SUBROUTINE kex_u ( derivative_set, functional, rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: functional
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
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
    IF ( mt > mp ) CALL stop_program ( "kex_u","Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_u", "r13", n )
    ALLOCATE ( s(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_u", "s", n )
    ALLOCATE ( fs(n,m+1), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_u", "fs", n*(m+1) )
    
    CALL calc_rho13 ( rho, r13 )
    CALL calc_wave_vector ( "u", rho, grho, s )
    SELECT CASE (functional)
    CASE (1)
      CALL efactor_ol1 ( s, fs, m )
    CASE (2)
      CALL efactor_ol2 ( s, fs, m )
    CASE (3)
      CALL efactor_llp ( s, fs, m )
    CASE (4)
      CALL efactor_pw86 ( s, fs, m )
    CASE (5)
      CALL efactor_pw91 ( s, fs, m, 1 )
    CASE (6)
      CALL efactor_pw91 ( s, fs, m, 2 )
    CASE (7)
      CALL efactor_t92 ( s, fs, m )
    CASE DEFAULT
       CALL stop_program ( "kex_u", "Unknown functional" )
    END SELECT

    fact = flda
    tact = t13

    IF ( calc(0) ) THEN
       CALL kex_p_0 ( rho, r13, fs, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL kex_p_1 ( rho, grho, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL kex_p_2 ( rho, grho, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL kex_p_3 ( rho, grho, r13, s, fs, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_u", "r13" )
    DEALLOCATE ( s, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_u", "s" )
    DEALLOCATE ( fs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_u", "fs" )

  END SUBROUTINE kex_u

! *****************************************************************************

  SUBROUTINE kex_p ( derivative_set, functional, rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    INTEGER, INTENT(IN)                      :: functional
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, grhoa, grhob
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot
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
    IF ( mt > mp ) CALL stop_program ( "kex_p", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_p", "r13", n )
    ALLOCATE ( s(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_p", "s", n )
    ALLOCATE ( fs(n,m+1), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_p", "fs", n*(m+1) )

    fact = flsd
    tact = 1._dbl

    CALL calc_rho13 ( rhoa, r13 )
    CALL calc_wave_vector ( "p", rhoa, grhoa, s )
    SELECT CASE (functional)
    CASE (1)
      CALL efactor_ol1 ( s, fs, m )
    CASE (2)
      CALL efactor_ol2 ( s, fs, m )
    CASE (3)
      CALL efactor_llp ( s, fs, m )
    CASE (4)
      CALL efactor_pw86 ( s, fs, m )
    CASE (5)
      CALL efactor_pw91 ( s, fs, m, 1 )
    CASE (6)
      CALL efactor_pw91 ( s, fs, m, 2 )
    CASE (7)
      CALL efactor_t92 ( s, fs, m )
    CASE DEFAULT
       CALL stop_program ( "kex_u", "Unknown functional" )
    END SELECT

    IF ( calc(0) ) THEN
       CALL kex_p_0 ( rhoa, r13, fs, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL kex_p_1 ( rhoa, grhoa, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL kex_p_2 ( rhoa, grhoa, r13, s, fs, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL kex_p_3 ( rhoa, grhoa, r13, s, fs, pot(:,m1:m2) )
    END IF

    CALL calc_rho13 ( rhob, r13 )
    CALL calc_wave_vector ( "p", rhob, grhob, s )
    SELECT CASE (functional)
    CASE (1)
      CALL efactor_ol1 ( s, fs, m )
    CASE (2)
      CALL efactor_ol2 ( s, fs, m )
    CASE (3)
      CALL efactor_llp ( s, fs, m )
    CASE (4)
      CALL efactor_pw86 ( s, fs, m )
    CASE (5)
      CALL efactor_pw91 ( s, fs, m, 1 )
    CASE (6)
      CALL efactor_pw91 ( s, fs, m, 2 )
    CASE (7)
      CALL efactor_t92 ( s, fs, m )
    CASE DEFAULT
       CALL stop_program ( "kex_u", "Unknown functional" )
    END SELECT

    IF ( calc(0) ) THEN
       CALL kex_p_0 ( rhob, r13, fs, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL kex_p_1 ( rhob, grhob, r13, s, fs, pot(:,m1+2:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL kex_p_2 ( rhob, grhob, r13, s, fs, pot(:,m1+3:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL kex_p_3 ( rhob, grhob, r13, s, fs, pot(:,m1+4:m2) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_p", "r13" )
    DEALLOCATE ( s, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_p", "s" )
    DEALLOCATE ( fs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "kex_p", "fs" )

  END SUBROUTINE kex_p

! *****************************************************************************

  SUBROUTINE kex_p_0 ( rho, r13, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: pot

    INTEGER                                  :: ip, npoints

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_p_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + fact*r13(ip)*r13(ip)*rho(ip) * fs(ip,1)
      END IF

    END DO

  END SUBROUTINE kex_p_0

! *****************************************************************************

  SUBROUTINE kex_p_1 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: a0, a1, sx, sy

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_p_1", "Inconsistent array sizes" )

!$omp parallel do private(ip,a0,a1,sx,sy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         pot(ip,1) = pot(ip,1) + a1*fs(ip,1) + a0*fs(ip,2)*sx
         pot(ip,2) = pot(ip,2) + a0*fs(ip,2)*sy

      END IF

    END DO

  END SUBROUTINE kex_p_1

! *****************************************************************************

  SUBROUTINE kex_p_2 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: a0, a1, a2, sx, sxx, sxy, sy

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_p_2", "Inconsistent array sizes" )

!$omp parallel do private(ip,a0,a1,a2,sx,sy,sxx,sxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
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

  END SUBROUTINE kex_p_2

! *****************************************************************************

  SUBROUTINE kex_p_3 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: a0, a1, a2, a3, sx, sxx, &
                                                sxxx, sxxy, sxy, sy

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_u_3", "Inconsistent array sizes" )

!$omp parallel do private(ip,a0,a1,a2,a3,sx,sy,sxx,sxy,sxxx,sxxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         a3 = -f13*f23*f53*fact/(r13(ip)*rho(ip))
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

  END SUBROUTINE kex_p_3

! *****************************************************************************
! Enhancement Factors
! *****************************************************************************
  SUBROUTINE efactor_ol1 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: t1, t2

    t1 = b*b/(72._dbl*cf)
    t2 = 0.001878_dbl*b
!$omp parallel do private(ip)
    DO ip=1,SIZE(s)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)
      CASE (1)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)
         fs(ip,2) = 2._dbl*t1*s(ip) + t2
      CASE (2:3)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)
         fs(ip,2) = 2._dbl*t1*s(ip) + t2
         fs(ip,3) = 2._dbl*t1
      CASE DEFAULT
         CALL stop_program("efactor_ol1","Illegal order")
      END SELECT
    END DO
    IF(m==3) fs(:,4)=0._dbl

  END SUBROUTINE efactor_ol1
! *****************************************************************************
  SUBROUTINE efactor_ol2 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: t1, t2, t3, y

    t1 = b*b/(72._dbl*cf)
    t2 = 0.0245_dbl*b
    t3 = 2._dbl**f53*b
!$omp parallel do private(ip,y)
    DO ip=1,SIZE(s)
      y = 1._dbl/(1._dbl+t3*s(ip))
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
      CASE (1)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2._dbl*t1*s(ip) + t2*y*y
      CASE (2)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2._dbl*t1*s(ip) + t2*y*y
         fs(ip,3) = 2._dbl*(t1 - t2*t3*y*y*y)
      CASE (3)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)*y
         fs(ip,2) = 2._dbl*t1*s(ip) + t2*y*y
         fs(ip,3) = 2._dbl*(t1 - t2*t3*y*y*y)
         fs(ip,4) = 6._dbl*t2*t3*t3*y*y*y*y
      CASE DEFAULT
         CALL stop_program("efactor_ol2","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_ol2
! *****************************************************************************
  SUBROUTINE efactor_llp ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl) :: as, asp, asp2, asp3, bs, p, q, sas, sbs, sbs3, sbs5, t1, &
      t11, t12, t13, t16, t19, t22, t24, t25, t3, t33, t34, t36, t39, t4, &
      t43, t44, t60, t62, w1, x, ys

    p = 0.0045_dbl*b*b
    q = 0.0253_dbl*b
!$omp parallel do private(ip,x,bs,sbs,as,sas,ys,asp,sbs3,asp2,asp3,sbs5,w1) &
!$omp             private(t1,t3,t4,t11,t12,t13,t16,t19,t22,t24,t25,t33,t34,t36,t39) &
!$omp             private(t43,t44,t60,t62)
    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = SQRT(bs*bs+1._dbl)
      as = LOG(bs+sbs)
      sas = x*as
      ys =  1._dbl/(1._dbl+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + p*x*x*ys
      CASE (1)
         asp = as + bs / sbs
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys
      CASE (2)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys
         fs(ip,3) = 2._dbl*p*ys - p*q*x*(4._dbl*asp+x*asp2)*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys
      CASE (3)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         sbs5 = sbs3*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         asp3 = -4._dbl*b*b*bs / sbs3 + 3._dbl*b*b*bs*bs*bs/ sbs5
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys
         w1 = (4._dbl*asp+x*asp2)
         fs(ip,3) = 2._dbl*p*ys - p*q*x*w1*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys
         t1 = q*x; t3 = b**2; t4 = x**2; t11 = 1+t1*as; t12 = t11**2; t13 = 1/t12
         t16 = 1/sbs; t19 = q*as+t1*b*t16; t22 = p*x; t24 = 1/t12/t11; t25 = t19**2
         t33 = t3*b; t34 = sbs**2; t36 = 1/t34/sbs; t39 = 2*q*b*t16-q*t4*t33*t36
         t43 = p*t4; t44 = t12**2; t60 = t3**2; t62 = t34**2
         fs(ip,4) = &
            -6*p*t13*t19+12*t22*t24*t25-6*t22*t13*t39-6*t43/t44*t25*t19+ &
            6*t43*t24*t19*t39-t43*t13*(-4*q*t33*t36*x+3*q*t4*x*t60*b/t62/sbs)

      CASE DEFAULT
         CALL stop_program("efactor_llp","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_llp
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: f15, p0, p1, p15, s2, s4, s6, &
                                                t1, t10, t12, t13, t14, t19, &
                                                t2, t25, t3, t8, t9

    t1 = 1.296_dbl
    t2 = 14._dbl
    t3 = 0.2_dbl
    f15 = 1._dbl/15._dbl
!$omp parallel do private(ip,s2,s4,s6,p0,p1,p15,t8,t9,t10,t12,t13,t14,t19,t25)
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
  SUBROUTINE efactor_t92 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: a1, a2, as, asp, asp2, asp3, &
                                                bs, p, q, s1, s2, sas, sbs, &
                                                sbs3, sbs5, t0, w1, x, ys

    p = 0.0055_dbl*b*b
    q = 0.0253_dbl*b
    a1 = 0.072_dbl * b
    a2 = 2._dbl**f53 * b
!$omp parallel do private(ip,x,as,sas,asp,asp2,asp3,bs,ys,sbs,sbs3,sbs5,w1,t0,s1,s2)
    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = SQRT(bs*bs+1._dbl)
      as = LOG(bs+sbs)
      sas = x*as
      ys =  1._dbl/(1._dbl+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
      CASE (1)
         asp = as + bs / sbs
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
      CASE (2)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
         fs(ip,3) = 2._dbl*p*ys - p*q*x*(4._dbl*asp+x*asp2)*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys + 2._dbl*a1*a2/(1+a2*x)**3
      CASE (3)
         asp = as + bs / sbs
         sbs3 = sbs*sbs*sbs
         sbs5 = sbs3*sbs*sbs
         asp2 = 2._dbl*b / sbs - b*bs*bs/ sbs3
         asp3 = -4._dbl*b*b*bs / sbs3 + 3._dbl*b*b*bs*bs*bs/ sbs5
         w1 = (4._dbl*asp+x*asp2)
         fs(ip,1) = 1._dbl + p*x*x*ys - a1*x/(1+a2*x)
         fs(ip,2) = 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys - a1/(1+a2*x)**2
         fs(ip,3) = 2._dbl*p*ys - p*q*x*w1*ys*ys + &
                    2._dbl*p*q*q*x*x*asp*asp*ys*ys*ys + 2._dbl*a1*a2/(1+a2*x)**3

      s2 = -6*p/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**2*(q*LOG(b*x+SQRT(1+b**2*x**2))+&
           q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))+12*p*x/&
           (1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**3*(q*LOG(b*x+SQRT(1+b**2*x**2))+&
           q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))**2
      s1 = s2-6*p*x/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**2*(2*q*(b+1/SQRT(1+b**2*x**2)*b**2*x)/&
           (b*x+SQRT(1+b**2*x**2))+q*x*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+1/SQRT(1+b**2*x**2)*b**2)/&
           (b*x+SQRT(1+b**2*x**2))-q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)**2/&
           (b*x+SQRT(1+b**2*x**2))**2)-6*p*x**2/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**4&
           *(q*LOG(b*x+SQRT(1+b**2*x**2))+q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))**3
      s2 = s1+6*p*x**2/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**3*(q*LOG(b*x+SQRT(1+b**2*x**2))+&
           q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)/(b*x+SQRT(1+b**2*x**2)))*(2*q*(b+1/SQRT(1+b**2*x**2)*b**2*x)&
           /(b*x+SQRT(1+b**2*x**2))+q*x*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+1/SQRT(1+b**2*x**2)*&
           b**2)/(b*x+SQRT(1+b**2*x**2))-q*x*(b+1/SQRT(1+b**2*x**2)*b**2*x)**2/(b*x+SQRT(1+b**2*x**2))**2)
      t0 = s2-p*x**2/(1+q*x*LOG(b*x+SQRT(1+b**2*x**2)))**2*(3*q*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+&
           1/SQRT(1+b**2*x**2)*b**2)/(b*x+SQRT(1+b**2*x**2))-3*q*(b+1/SQRT(1+b**2*x**2)*b**2*x)**2/&
           (b*x+SQRT(1+b**2*x**2))**2+q*x*(3/SQRT(1+b**2*x**2)**5*b**6*x**3-3/SQRT(1+b**2*x**2)**3*b**4*x)/&
           (b*x+SQRT(1+b**2*x**2))-3*q*x*(-1/SQRT(1+b**2*x**2)**3*b**4*x**2+1/SQRT(1+b**2*x**2)*b**2)/&
           (b*x+SQRT(1+b**2*x**2))**2*(b+1/SQRT(1+b**2*x**2)*b**2*x)+2*q*x*(b+1/SQRT(1+b**2*x**2)*&
           b**2*x)**3/(b*x+SQRT(1+b**2*x**2))**3)-6*a1/(1+a2*x)**3*a2**2+6*a1*x/(1+a2*x)**4*a2**3

         fs(ip,4) = t0
!deb     fs(ip,4) = -2._dbl*p*q*asp*ys*ys+p*q*x*(5._dbl*asp2-x*asp3)*ys*ys - &
!deb                p*q*w1*ys*ys + 2._dbl*p*q*q*x*asp*w1*ys*ys*ys - &
!deb                6._dbl*p*q*q*q*x*x*asp*asp*asp*ys*ys*ys*ys + &
!deb                4._dbl*p*q*q*ys*ys*ys*(x*asp*asp+x*x*asp*asp2) - &
!deb                6._dbl * a1*a2*a2/(1+a2*x)**4
      CASE DEFAULT
         CALL stop_program("efactor_t92","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_t92
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m, pset )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m, pset

    INTEGER                                  :: ip
    REAL(dbl)                                :: a, a1, a2, a3, a4, b1, o, &
                                                pa(6,2), s1, s2, s3, x

! parameter set 1: Perdew-Wang
! parameter set 2: Lembarki-Chermette

    pa(1:6,1)=(/0.19645_dbl,    0.2743_dbl,  &
                0.1508_dbl,   100._dbl,      &
                0.8145161_dbl,  0.004_dbl/)
    pa(1:6,2)=(/0.093907_dbl,   0.26608_dbl, &
                0.0809615_dbl,100._dbl,      &
               76.320_dbl,      0.57767e-4_dbl/)
    o=1._dbl
    a1 = pa(1,pset)
    a2 = pa(2,pset)
    a3 = pa(3,pset)
    a4 = pa(4,pset)
    a  = pa(5,pset)
    b1 = pa(6,pset)
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
    IF ( m >= 2 ) THEN
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

END MODULE xc_ke_gga

! *****************************************************************************
