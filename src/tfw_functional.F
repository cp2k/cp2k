!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/tfw_functional [1.0] *
!!
!!   NAME
!!     tfw_functional
!!
!!   FUNCTION
!!     Calculate the Thomas-Fermi kinetic energy functional
!!     plus the von Weizsaecker term
!!
!!   AUTHOR
!!     JGH (18.02.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************

MODULE tfw_functional

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE functionals_utilities, ONLY: set_util, calc_rho13
  USE termination, ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: tfw

  REAL(dbl) :: cf, flda, flsd, fvw
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE tfw
    MODULE PROCEDURE tfw_u, tfw_p, tfw_init, tfw_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE tfw_init ( cutoff )

    REAL(dbl), INTENT(IN) :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23
    flda = cf
    flsd = flda * 2._dbl**f23
    fvw  = 1._dbl/72._dbl

  END SUBROUTINE tfw_init

! *****************************************************************************

  SUBROUTINE tfw_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "Thomas-Fermi-Weizsaecker kinetic energy functional"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "TF+vW kinetic energy functional"
    END IF

  END SUBROUTINE tfw_info

! *****************************************************************************

  SUBROUTINE tfw_u ( rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13, s
    INTEGER :: n, error, m, mp, mt, ip

!   ---------------------------------------------------------------------------

    calc = .FALSE.
    mp = SIZE ( pot, 2 )
    allcalc = .TRUE.
    IF ( PRESENT ( order ) ) THEN
      IF ( order >= 0 ) THEN
        calc(0:order) = .TRUE.
        SELECT CASE (order)
        CASE (0)
          mt = 1
        CASE (1)
          mt = 3
        CASE (2)
          mt = 6
        CASE (3)
          mt = 10
        CASE DEFAULT
          mt = mp+100
        END SELECT
      ELSE
        calc(-order) = .TRUE.
        allcalc = .FALSE.
        mt = -order+1
      END IF
    ELSE
       calc(0) = .TRUE.
       calc(1) = .TRUE.
       mt = 3
    END IF
    IF ( abs(order) > 3 ) CALL stop_program ( "tfw_u", "order larger than 3" )
    IF ( mt > mp ) CALL stop_program ( "tfw_u", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_u", "r13", n )
    ALLOCATE ( s(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_u", "s", n )
    
    m = 1
    CALL calc_rho13 ( rho, r13 )
    DO ip = 1, n
      IF ( rho(ip) < eps_rho ) THEN
         s(ip) = 0.0_dbl
      ELSE
         s(ip) = grho(ip)*grho(ip) / rho(ip)
      END IF
    END DO
    IF ( calc(0) ) THEN
       CALL tfw_u_0 ( rho, r13, s, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       IF ( allcalc ) m = 2
       CALL tfw_u_1 ( rho, grho, r13, s, pot(:,m:m+1) )
    END IF
    IF ( calc(2) ) THEN
       IF ( allcalc ) m = 4
       CALL tfw_u_2 ( rho, grho, r13, s, pot(:,m:m+2) )
    END IF
    IF ( calc(3) ) THEN
       IF ( allcalc ) m = 7
       CALL tfw_u_3 ( rho, grho, r13, s, pot(:,m:m+3) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_u", "r13" )
    DEALLOCATE ( s, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_u", "s" )

  END SUBROUTINE tfw_u

! *****************************************************************************

  SUBROUTINE tfw_p ( rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13a, r13b, sa, sb
    INTEGER :: n, error, m, mp, mt, ip

!   ---------------------------------------------------------------------------

    calc = .FALSE.
    mp = SIZE ( pot, 2 )
    allcalc = .TRUE.
    IF ( PRESENT ( order ) ) THEN
      IF ( order >= 0 ) THEN
        calc(0:order) = .TRUE.
        SELECT CASE (order)
        CASE (0)
          mt = 1
        CASE (1)
          mt = 5
        CASE (2)
          mt = 11
        CASE (3)
          mt = 19
        CASE DEFAULT
          mt = mp+100
        END SELECT
      ELSE
        calc(-order) = .TRUE.
        allcalc = .FALSE.
        mt = 2*(ABS(order)+1)
      END IF
    ELSE
       calc(0) = .TRUE.
       calc(1) = .TRUE.
       mt = 3
    END IF
    IF ( abs(order) > 3 ) CALL stop_program ( "tfw_u", "order larger than 3" )
    IF ( mt > mp ) CALL stop_program ( "tfw_p", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_p", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_p", "r13b", n )
    ALLOCATE ( sa(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_p", "sa", n )
    ALLOCATE ( sb(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_p", "sb", n )

    m = 1
    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    DO ip = 1, n
      IF ( rhoa(ip) < eps_rho ) THEN
         sa(ip) = 0.0_dbl
      ELSE
         sa(ip) = grhoa(ip)*grhoa(ip) / rhoa(ip)
      END IF
    END DO
    DO ip = 1, n
      IF ( rhob(ip) < eps_rho ) THEN
         sb(ip) = 0.0_dbl
      ELSE
         sb(ip) = grhob(ip)*grhob(ip) / rhob(ip)
      END IF
    END DO
    DO ip = 1, n
      IF ( rhoa(ip) < eps_rho ) THEN
         sa(ip) = 0.0_dbl
      ELSE
         sa(ip) = grhoa(ip)*grhoa(ip) / rhoa(ip)
      END IF
    END DO
    DO ip = 1, n
      IF ( rhob(ip) < eps_rho ) THEN
         sb(ip) = 0.0_dbl
      ELSE
         sb(ip) = grhob(ip)*grhob(ip) / rhob(ip)
      END IF
    END DO
    IF ( calc(0) ) THEN
       CALL tfw_p_0 ( rhoa, rhob, grhoa, grhob, &
                      r13a, r13b, sa, sb, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       IF ( allcalc ) m = 2
       CALL tfw_p_1 ( rhoa, rhob, grhoa, grhob, &
                      r13a, r13b, sa, sb, pot(:,m:m+3) )
    END IF
    IF ( calc(2) ) THEN
       IF ( allcalc ) m = 6
       CALL tfw_p_2 ( rhoa, rhob, grhoa, grhob, &
                      r13a, r13b, sa, sb, pot(:,m:m+5) )
    END IF
    IF ( calc(3) ) THEN
       IF ( allcalc ) m = 12
       CALL tfw_p_3 ( rhoa, rhob, grhoa, grhob, &
                      r13a, r13b, sa, sb, pot(:,m:m+7) )
    END IF

    DEALLOCATE ( r13a, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_p", "r13a" )
    DEALLOCATE ( r13b, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "tfw_p", "r13b" )

  END SUBROUTINE tfw_p

! *****************************************************************************

  SUBROUTINE tfw_u_0 ( rho, r13, s, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13, s
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "tfw_u_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = flda*r13(ip)*r13(ip)*rho(ip) + fvw * s(ip)

      END IF

    END DO

  END SUBROUTINE tfw_u_0

! *****************************************************************************

  SUBROUTINE tfw_u_1 ( rho, grho, r13, s, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "tfw_u_1", "Inconsistent array sizes" )

    f = f53 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl

      ELSE

         pot(ip,1) = f * r13(ip)*r13(ip) - fvw * s(ip)/rho(ip)
         pot(ip,2) = 2._dbl * fvw * grho(ip)/rho(ip)

      END IF

    END DO

  END SUBROUTINE tfw_u_1

! *****************************************************************************

  SUBROUTINE tfw_u_2 ( rho, grho, r13, s, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "tfw_u_2", "Inconsistent array sizes" )

    f = f23 * f53 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl

      ELSE

         pot(ip,1) = f / r13(ip) + 2._dbl * fvw * s(ip)/(rho(ip)*rho(ip))
         pot(ip,2) = -2._dbl * fvw * grho(ip)/(rho(ip)*rho(ip))
         pot(ip,3) = 2._dbl * fvw / rho(ip)

      END IF

    END DO

  END SUBROUTINE tfw_u_2

! *****************************************************************************

  SUBROUTINE tfw_u_3 ( rho, grho, r13, s, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "tfw_u_3", "Inconsistent array sizes" )

    f = -f13 * f23 * f53 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl

      ELSE

         pot(ip,1) = f / ( r13(ip) * rho(ip) ) &
                     - 6._dbl * fvw * s(ip)/(rho(ip)*rho(ip)*rho(ip))
         pot(ip,2) = 4._dbl * fvw * grho(ip)/(rho(ip)*rho(ip)*rho(ip))
         pot(ip,3) = - 2._dbl * fvw / (rho(ip)*rho(ip))
         pot(ip,4) = 0.0_dbl

      END IF

    END DO

  END SUBROUTINE tfw_u_3

! *****************************************************************************

  SUBROUTINE tfw_p_0 ( rhoa, rhob, grhoa, grhob, r13a, r13b, sa, sb, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b, sa, sb
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "tfw_p_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip) = 0.0_dbl
      ELSE
         pot(ip) = flsd * r13a(ip) * r13a(ip) * rhoa(ip) + fvw * sa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13b(ip) * r13b(ip) * rhob(ip) + fvw*sb(ip)
      END IF

    END DO

  END SUBROUTINE tfw_p_0

! *****************************************************************************

  SUBROUTINE tfw_p_1 ( rhoa, rhob, grhoa, grhob, r13a, r13b, sa, sb, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b, sa, sb
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "tfw_p_1", "Inconsistent array sizes" )

    f = f53 * flsd

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,1) = f * r13a(ip)*r13a(ip) - fvw * sa(ip)/rhoa(ip)
         pot(ip,2) = 2._dbl * fvw * grhoa(ip)/rhoa(ip)
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      ELSE
         pot(ip,3) = f * r13b(ip)*r13b(ip) - fvw * sb(ip)/rhob(ip)
         pot(ip,4) = 2._dbl * fvw * grhob(ip)/rhob(ip)
      END IF

    END DO

  END SUBROUTINE tfw_p_1

! *****************************************************************************

  SUBROUTINE tfw_p_2 ( rhoa, rhob, grhoa, grhob, r13a, r13b, sa, sb, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b, sa, sb
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "tfw_p_2", "Inconsistent array sizes" )

    f = f23 * f53 * flsd

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
      ELSE
         pot(ip,1) = f / r13a(ip) + 2._dbl * fvw * sa(ip)/(rhoa(ip)*rhoa(ip))
         pot(ip,2) = -2._dbl * fvw * grhoa(ip)/(rhoa(ip)*rhoa(ip))
         pot(ip,3) = 2._dbl * fvw / rhoa(ip)
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,4) = 0.0_dbl
         pot(ip,5) = 0.0_dbl
         pot(ip,6) = 0.0_dbl
      ELSE
         pot(ip,4) = f / r13b(ip) + 2._dbl * fvw * sb(ip)/(rhob(ip)*rhob(ip))
         pot(ip,5) = -2._dbl * fvw * grhob(ip)/(rhob(ip)*rhob(ip))
         pot(ip,6) = 2._dbl * fvw / rhob(ip)
      END IF

    END DO

  END SUBROUTINE tfw_p_2

! *****************************************************************************

  SUBROUTINE tfw_p_3 ( rhoa, rhob, grhoa, grhob, r13a, r13b, sa, sb, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b, sa, sb
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "tfw_p_3", "Inconsistent array sizes" )

    f = -f13 * f23 * f53 * flsd

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      ELSE
         pot(ip,1) = f / ( r13a(ip) * rhoa(ip) ) &
                     - 6._dbl * fvw * sa(ip)/(rhoa(ip)*rhoa(ip)*rhoa(ip))
         pot(ip,2) = 4._dbl * fvw * grhoa(ip)/(rhoa(ip)*rhoa(ip)*rhoa(ip))
         pot(ip,3) = - 2._dbl * fvw / (rhoa(ip)*rhoa(ip))
         pot(ip,4) = 0.0_dbl
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,5) = 0.0_dbl
         pot(ip,6) = 0.0_dbl
         pot(ip,7) = 0.0_dbl
         pot(ip,8) = 0.0_dbl
      ELSE
         pot(ip,5) = f / ( r13b(ip) * rhob(ip) ) &
                     - 6._dbl * fvw * sb(ip)/(rhob(ip)*rhob(ip)*rhob(ip))
         pot(ip,6) = 4._dbl * fvw * grhob(ip)/(rhob(ip)*rhob(ip)*rhob(ip))
         pot(ip,7) = - 2._dbl * fvw / (rhob(ip)*rhob(ip))
         pot(ip,8) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE tfw_p_3

! *****************************************************************************

END MODULE tfw_functional

! *****************************************************************************
