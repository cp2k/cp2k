!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/ke_gga_functional [1.0] *
!!
!!   NAME
!!     ke_gga_functional
!!
!!   FUNCTION
!!     Calculate the several different kinetic energy functionals
!!     with a GGA form 
!!
!!   AUTHOR
!!     JGH (20.02.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************

MODULE ke_gga_functional

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE functionals_utilities, ONLY: set_util, calc_rho13, setup_calculation, &
                                   calc_wave_vector
  USE termination, ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: ol1, ol2, llp, pw86, pw91, t92, lc

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
  INTERFACE pw86
    MODULE PROCEDURE pw86_u, pw86_p, kegga_init, pw86_info
  END INTERFACE
  INTERFACE pw91
    MODULE PROCEDURE pw91_u, pw91_p, kegga_init, pw91_info
  END INTERFACE
  INTERFACE lc
    MODULE PROCEDURE lc_u, lc_p, kegga_init, lc_info
  END INTERFACE
  INTERFACE t92
    MODULE PROCEDURE t92_u, t92_p, kegga_init, t92_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE kegga_init ( cutoff )

    REAL(dbl), INTENT(IN) :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23
    flda = cf
    flsd = flda * 2._dbl**f23
    fvw  = 1._dbl/72._dbl
    b = 2._dbl**f43 * (3._dbl*pi*pi)**(-f13)
    sfac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)
    t13 = 2._dbl**f13

  END SUBROUTINE kegga_init

! *****************************************************************************

  SUBROUTINE ol1_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "H. Ou-Yang and M. Levy, ",//&
                  "Intl. J. Quant. Chem. 40, 379 (1991); Functional 1"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Ou-Yang-Levy Functional 1"
    END IF

  END SUBROUTINE ol1_info

  SUBROUTINE ol2_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "H. Ou-Yang and M. Levy, ",//&
                  "Intl. J. Quant. Chem. 40, 379 (1991); Functional 2"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Ou-Yang-Levy Functional 2"
    END IF

  END SUBROUTINE ol2_info

  SUBROUTINE llp_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "H. Lee, C. Lee, R.G. Parr, Phys. Rev. A, 44, 768 (1991)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Lee-Lee-Parr Functional"
    END IF

  END SUBROUTINE llp_info

  SUBROUTINE pw86_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800 (1986)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew-Wang 1986 Functional (kinetic energy)"
    END IF

  END SUBROUTINE pw86_info

  SUBROUTINE pw91_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J.P. Perdew and Y. Wang, Electronic Structure of Solids 91"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew-Wang 1991 Functional (kinetic energy)"
    END IF

  END SUBROUTINE pw91_info

  SUBROUTINE lc_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "A. Lembarki and H. Chermette, Phys. Rev. A, 50, 5328 (1994)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Lembarki-Chermette kinetic energy functional"
    END IF

  END SUBROUTINE lc_info

  SUBROUTINE t92_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "A.J. Thakkar, Phys. Rev. A, 46, 6920 (1992)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Thakkar 1992 Functional"
    END IF

  END SUBROUTINE t92_info

! *****************************************************************************

  SUBROUTINE ol1_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 1, rho, grho, pot, order )
  END SUBROUTINE ol1_u
  SUBROUTINE ol1_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 1, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE ol1_p ( rhoa, rhob, grhoa, grhob, pot, order )

  SUBROUTINE ol2_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 2, rho, grho, pot, order )
  END SUBROUTINE ol2_u
  SUBROUTINE ol2_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 2, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE ol2_p ( rhoa, rhob, grhoa, grhob, pot, order )

  SUBROUTINE llp_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 3, rho, grho, pot, order )
  END SUBROUTINE llp_u
  SUBROUTINE llp_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 3, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE llp_p ( rhoa, rhob, grhoa, grhob, pot, order )

  SUBROUTINE pw86_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 4, rho, grho, pot, order )
  END SUBROUTINE pw86_u
  SUBROUTINE pw86_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 4, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pw86_p ( rhoa, rhob, grhoa, grhob, pot, order )

  SUBROUTINE pw91_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 5, rho, grho, pot, order )
  END SUBROUTINE pw91_u
  SUBROUTINE pw91_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 5, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE pw91_p ( rhoa, rhob, grhoa, grhob, pot, order )

  SUBROUTINE lc_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 6, rho, grho, pot, order )
  END SUBROUTINE lc_u
  SUBROUTINE lc1_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 6, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE lc1_p ( rhoa, rhob, grhoa, grhob, pot, order )

  SUBROUTINE t92_u ( rho, grho, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_u ( 7, rho, grho, pot, order )
  END SUBROUTINE t92_u
  SUBROUTINE t92_p ( rhoa, rhob, grhoa, grhob, pot, order )
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

    CALL kex_p ( 7, rhoa, rhob, grhoa, grhob, pot, order )
  END SUBROUTINE t92_p ( rhoa, rhob, grhoa, grhob, pot, order )

! *****************************************************************************

  SUBROUTINE kex_u ( functional, rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: functional
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13, s
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: fs
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2

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

  SUBROUTINE kex_p ( functional, rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: functional
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13, s
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: fs
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2

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
       pot(:,1) = 0._dbl
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

  END SUBROUTINE tfw_p

! *****************************************************************************

  SUBROUTINE kex_p_0 ( rho, r13, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fs
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_p_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = fact*r13(ip)*r13(ip)*rho(ip) * fs(ip,1)

      END IF

    END DO

  END SUBROUTINE kex_p_0

! *****************************************************************************

  SUBROUTINE kex_p_1 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: fs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: a0, a1, sx, sy
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_p_1", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl

      ELSE

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         pot(ip,1) = a1*fs(ip,1) + a0*fs(ip,2)*sx
         pot(ip,2) = a0*fs(ip,2)*sy

      END IF

    END DO

  END SUBROUTINE kex_p_1

! *****************************************************************************

  SUBROUTINE kex_p_2 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: fs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: a0, a1, a2, sx, sxx, sy, sxy
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_p_2", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl

      ELSE

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         pot(ip,1) = a2*fs(ip,1) +2._dbl*a1*fs(ip,2)*sx + a0*fs(ip,3)*sx*sx +&
                     a0*fs(ip,2)*sxx
         pot(ip,2) = a1*fs(ip,2)*sy + a0*fs(ip,3)*sx*sy + a0*fs(ip,2)*sxy
         pot(ip,3) = a0*fs(ip,3)*sy*sy

      END IF

    END DO

  END SUBROUTINE kex_p_2

! *****************************************************************************

  SUBROUTINE kex_p_3 ( rho, grho, r13, s, fs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: fs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: a0, a1, a2, a3, sx, sxx, sxxx, sy, sxy, sxxy
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "kex_u_3", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl

      ELSE

         a0 = fact*r13(ip)*r13(ip)*rho(ip)
         a1 = f53*fact*r13(ip)*r13(ip)
         a2 = f23*f53*fact/r13(ip)
         a3 = -f13*f23*f53*fact/(r13(ip)*rho(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         sxxx= 280._dbl/27._dbl*s(ip)/(rho(ip)*rho(ip)*rho(ip))
         sxxy= 28._dbl/9._dbl*sfac*tact/(r13(ip)*rho(ip)*rho(ip)*rho(ip))
         pot(ip,1) = a3*fs(ip,1) + 3._dbl*a2*fs(ip,2)*sx + &
                     3._dbl*a1*fs(ip,3)*sx*sx + 3._dbl*a1*fs(ip,2)*sxx + &
                     a0*fs(ip,4)*sx*sx*sx + 3._dbl*a0*fs(ip,3)*sx*sxx + &
                     a0*fs(ip,2)*sxxx
         pot(ip,2) = a2*fs(ip,2)*sy + 2._dbl*a1*fs(ip,3)*sx*sy + &
                     2._dbl*a1*fs(ip,2)*sxy + a0*fs(ip,3)*sx*sx*sy + &
                     2._dbl*a0*fs(ip,3)*sx*sxy + a0*fs(ip,3)*sxx*sy + &
                     a0*fs(ip,2)*sxxy
         pot(ip,3) = a1*fs(ip,3)*sy*sy + a0*fs(ip,4)*sx*sy*sy + &
                     2._dbl*a0*fs(ip,3)*sxy*sy
         pot(ip,4) = a0*fs(ip,4)*sy*sy*sy

      END IF

    END DO

  END SUBROUTINE kex_p_3

! *****************************************************************************
! Enhancement Factors
! *****************************************************************************
  SUBROUTINE efactor_ol1 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN) :: s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fs
    INTEGER, INTENT(IN) :: m

    INTEGER :: ip
    REAL(dbl) :: t1, t2

    t1 = b*b/(72._dbl*cf)
    t2 = 0.001878_dbl*b
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
    REAL(dbl), DIMENSION(:), INTENT(IN) :: s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fs
    INTEGER, INTENT(IN) :: m

    INTEGER :: ip
    REAL(dbl) :: t1, t2, t3

    t1 = b*b/(72._dbl*cf)
    t2 = 0.0245_dbl*b
    t3 = 2._dbl**f53*b
    DO ip=1,SIZE(s)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)/(1._dbl+t3*s(ip))
      CASE (1)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)/(1._dbl+t3*s(ip))
         fs(ip,2) = 2._dbl*t1*s(ip) + t2/(1._dbl+t3*s(ip))**2
      CASE (2)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)/(1._dbl+t3*s(ip))
         fs(ip,2) = 2._dbl*t1*s(ip) + t2/(1._dbl+t3*s(ip))**2
         fs(ip,3) = 2._dbl*t1 - t2*t3/(1._dbl+t3*s(ip))**3
      CASE (3)
         fs(ip,1) = 1._dbl + t1*s(ip)*s(ip) + t2*s(ip)/(1._dbl+t3*s(ip))
         fs(ip,2) = 2._dbl*t1*s(ip) + t2/(1._dbl+t3*s(ip))**2
         fs(ip,3) = 2._dbl*t1 - 2._dbl*t2*t3/(1._dbl+t3*s(ip))**3
         fs(ip,4) = - 2._dbl*t2*t3*t3/(1._dbl+t3*s(ip))**4
      CASE DEFAULT
         CALL stop_program("efactor_ol2","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_ol2
! *****************************************************************************
  SUBROUTINE efactor_llp ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN) :: s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fs
    INTEGER, INTENT(IN) :: m

    INTEGER :: ip
    REAL(dbl) :: p, q, x, as, sas, asp, asp2, asp3, bs, ys, sbs, sbs3, sbs5
    REAL(dbl) :: w1, w2

    p = 0.0045_dbl*b*b
    q = 0.0253_dbl*b
    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = sqrt(bs*bs+1._dbl)
      as = log(bs+sbs)
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
         fs(ip,4) = -2._dbl*p*q*asp*ys*ys p*q*x*(5._dbl*asp2-x*asp3)*ys*ys - &
                    p*q*w1*ys*ys + 2._dbl*p*q*q*x*asp*w1*ys*ys*ys - &
                    6._dbl*p*q*q*q*x*x*asp*asp*asp*ys*ys*ys*ys +
                    4._dbl*p*q*q*ys*ys*ys*(x*asp*asp+x*x*asp*asp2)
      CASE DEFAULT
         CALL stop_program("efactor_llp","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_llp
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN) :: s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fs
    INTEGER, INTENT(IN) :: m

    INTEGER :: ip
    REAL(dbl) :: t1, t2, t3, f15, s2, s4, s6
    REAL(dbl) :: p0, p1, p2, p3, p15

    t1 = 1.296_dbl
    t2 = 14._dbl
    t3 = 0.2_dbl
    f15 = 1._dbl/15._dbl
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
         fs(ip,3) = f15 * p15/p0 * ( p2 - 14._dbl/15._dbl*p1*p1/p )
      CASE (3)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         fs(ip,3) = f15 * p15/p0 * ( p2 - 14._dbl/15._dbl*p1*p1/p0 )
         fs(ip,4) = -14._dbl*15._dbl*f15 * p15/(p0*p0) &
                    * ( p2 - 14._dbl/15._dbl*p1*p1/p0 ) + &
                    f15 * p15/p0 * ( p3 + 14._dbl/15._dbl*p1*p1/(p0*p0) - &
                    28._dbl/15._dbl*p1*p2/p0 )
      CASE DEFAULT
         CALL stop_program("efactor_ol1","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_pw86
! *****************************************************************************
  SUBROUTINE efactor_t92 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN) :: s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fs
    INTEGER, INTENT(IN) :: m

    INTEGER :: ip
    REAL(dbl) :: p, q, x, as, sas, asp, asp2, asp3, bs, ys, sbs, sbs3, sbs5
    REAL(dbl) :: w1, w2, a1, a2

    p = 0.0055_dbl*b*b
    q = 0.0253_dbl*b
    a1 = 0.072 * b
    a2 = 2._dbl**f53 * b
    DO ip=1,SIZE(s)
      x = s(ip)
      bs = b*x
      sbs = sqrt(bs*bs+1._dbl)
      as = log(bs+sbs)
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
         fs(ip,4) = -2._dbl*p*q*asp*ys*ys p*q*x*(5._dbl*asp2-x*asp3)*ys*ys - &
                    p*q*w1*ys*ys + 2._dbl*p*q*q*x*asp*w1*ys*ys*ys - &
                    6._dbl*p*q*q*q*x*x*asp*asp*asp*ys*ys*ys*ys +
                    4._dbl*p*q*q*ys*ys*ys*(x*asp*asp+x*x*asp*asp2) - &
                    6._dbl * a1*a2*a2/(1+a2*x)**4
      CASE DEFAULT
         CALL stop_program("efactor_t92","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_t92
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m, pset )
    REAL(dbl), DIMENSION(:), INTENT(IN) :: s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fs
    INTEGER, INTENT(IN) :: m

    ! parameter set 1: Perdew-Wang
    ! parameter set 2: Lembarki-Chermette
    REAL(dbl), PARAMETER :: &
     pa(6,2)=(/0.19645_dbl,0.2743_dbl,0.1508_dbl,100._dbl,&
               0.8145161_dbl,0.004_dbl,
               0.093907_dbl,0.26608_dbl,0.0809615_dbl,100._dbl,&
               76.320_dbl,0.57767e-4_dbl/)
    INTEGER :: ip
    REAL(dbl) :: a1,a2,a3,a4,a,b1,o,s1,s2,s3

    o=1._dbl
    a1 = pa(1,pset)
    a2 = pa(2,pset)
    a3 = pa(3,pset)
    a4 = pa(4,pset)
    a  = pa(5,pset)
    b1 = pa(6,pset)
    IF ( m >= 0 ) THEN
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,1) = (o+a1*x*log(a*x+sqrt(o+a**2*x**2))+a2-a3*exp(-a4*x**2)*x**2)&
                   /(o+a1*x*log(a*x+sqrt(o+a**2*x**2))+b1*x**4)
      END DO
    END IF
    IF ( m >= 1 ) THEN
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,2) =&
           (a1*log(a*x+sqrt(o+a**2*x**2))+a1*x*a/sqrt(o+a**2*x**2)+2*a3*       &
           a4*x**3*exp(-a4*x**2)-2*a3*exp(-a4*x**2)*x)/(o+a1*x*log(a*x+sqrt(o  &
           +a**2*x**2))+b1*x**4)-(o+a1*x*log(a*x+sqrt(o+a**2*x**2))+a2-a3*exp  &
           (-a4*x**2)*x**2)/(o+a1*x*log(a*x+sqrt(o+a**2*x**2))+b1*x**4)**2*(a1 &
           *log(a*x+sqrt(o+a**2*x**2))+a1*x*a/sqrt(o+a**2*x**2)+4*b1*x**3)
      END DO
    END IF
    IF ( m >= 2 ) THEN
      DO ip=1,SIZE(s)
        x = s(ip)
        s1 = (2*a1*a/sqrt(1+a**2*x**2)-a1*x**2*a**3/sqrt(1+a**2*x**2)**3+  &
             10*a3*a4*x**2*exp(-a4*x**2)-4*a3*a4**2*x**4*exp(-a4*x**2)     &
             -2*a3*exp(-a4*x**2))/(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+      &
             b1*x**4)-2*(a1*log(a*x+sqrt(1+a**2*x**2))+a1*x*a/sqrt(1+a**2*x**2) &
             +2*a3*a4*x**3*exp(-a4*x**2)-2*a3*exp(-a4*x**2)*x)/(1+a1*x*    &
             log(a*x+sqrt(1+a**2*x**2))+b1*x**4)**2*(a1*log(a*x+sqrt(1+a**2*x**2))&
             +a1*x*a/sqrt(1+a**2*x**2)+4*b1*x**3)
        fs(ip,3) =&
             s1+2*(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+a2-a3*exp(-a4*x**2)*x**2)&
             /(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+b1*x**4)**3*(a1*log(a*x+     &
             sqrt(1+a**2*x**2))+a1*x*a/sqrt(1+a**2*x**2)+4*b1*x**3)**2-       &
             (1+a1*x*log(a*x+sqrt(1+a**2*x**2))+a2-a3*exp(-a4*x**2)*x**2)     &
             /(1+a1*x*log(a*x +sqrt(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/        &
             sqrt(1+a**2*x**2)-a1*x**2*a**3/sqrt(1+a**2*x**2)**3+12*b1*x**2)
      END DO
    END IF
    IF ( m >= 2 ) THEN
      DO ip=1,SIZE(s)
        x = s(ip)
        s2 = (-4*a1*a**3/sqrt(1+a**2*x**2)**3*x+3*a1*x**3*a**5/         &
             sqrt(1+a**2*x**2)**5+24*a3*a4*x*exp(-a4*x**2)-36*a3*a4**2* &
             x**3*exp(-a4*x**2)+8*a3*a4**3*x**5*exp(-a4*x**2))/(1+a1*x* &
             log(a*x+sqrt(1+a**2*x**2))+b1*x**4)
        s3 = -3*(2*a1*a/sqrt(1+a**2*x**2)-a1*x**2*a**3/sqrt(1+a**2*     &
             x**2)**3+10*a3*a4*x**2*exp(-a4*x**2)-4*a3*a4**2*x**4*      &
             exp(-a4*x**2)-2*a3*exp(-a4*x**2))/(1+a1*x*log(a*x+         &
             sqrt(1+a**2*x**2))+b1*x**4)**2*(a1*log(a*x+sqrt(1+a**2*    &
             x**2))+a1*x*a/sqrt(1+a**2*x**2)+4*b1*x**3)+6*(a1*log(a*x+  &
             sqrt(1+a**2*x**2))+a1*x*a/sqrt(1+a**2*x**2)+2*a3*a4*x**3*  &
             exp(-a4*x**2)-2*a3*exp(-a4*x**2)*x)/(1+a1*x*log(a*x+       &
             sqrt(1+a**2*x**2))+b1*x**4)**3*(a1*log(a*x+sqrt(1+a**2*    &
             x**2))+a1*x*a/sqrt(1+a**2*x**2)+4*b1*x**3)**2
        s1 = s2+s3
        s2 = s1-3*(a1*log(a*x+sqrt(1+a**2*x**2))+a1*x*a/sqrt(1+a**2*x**2) &
             +2*a3*a4*x**3*exp(-a4*x**2)-2*a3*exp(-a4*x**2)*x)/(1+a1*x*   &
             log(a*x+sqrt(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/sqrt(1+a**2*  &
             x**2)-a1*x**2*a**3/sqrt(1+a**2*x**2)**3+12*b1*x**2)
        s3 = s2-6*(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+a2-a3*exp(-a4*x**2)* &
             x**2)/(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+b1*x**4)**4*(a1*    &
             log(a*x+sqrt(1+a**2*x**2))+a1*x*a/sqrt(1+a**2*x**2)+4*b1*x**3)**3
        fs(ip,4) =&
             s3+6*(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+a2-a3*exp(-a4*x**2)  &
             *x**2)/(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+b1*x**4)**3*(a1*   &
             log(a*x+sqrt(1+a**2*x**2))+a1*x*a/sqrt(1+a**2*x**2)+4*b1*    &
             x**3)*(2*a1*a/sqrt(1+a**2*x**2)-a1*x**2*a**3/sqrt(1+a**2*    &
             x**2)**3+12*b1*x**2)-(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+a2-  &
             a3*exp(-a4*x**2)*x**2)/(1+a1*x*log(a*x+sqrt(1+a**2*x**2))+   &
             b1*x**4)**2*(-4*a1*a**3/sqrt(1+a**2*x**2)**3*x+3*a1*x**3*    &
             a**5/sqrt(1+a**2*x**2)**5+24*b1*x)
      END DO
    END IF

  END SUBROUTINE efactor_pw91

! *****************************************************************************

END MODULE ke_gga_functional

! *****************************************************************************
