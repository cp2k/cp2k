!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/perdew86_functional [1.0] *
!!
!!   NAME
!!     perdew86_functional
!!
!!   FUNCTION
!!     Calculate the Perdew Correlation from 1986
!!
!!   AUTHOR
!!     JGH (03.03.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************

MODULE lyp_functional

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE functionals_utilities, ONLY: set_util, setup_calculation, calc_rs
  USE termination, ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: perdew86

  REAL(dbl) :: eps_rho

  REAL(dbl), PARAMETER :: a = 0.04918_dbl, &
                          b = 0.132_dbl, &
                          c = 0.2533_dbl, &
                          d = 0.349_dbl

! *****************************************************************************

  INTERFACE perdew86
    MODULE PROCEDURE p86_u, p86_p, p86_init, p86_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE p86_init ( cutoff )

    REAL(dbl), INTENT(IN) :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

  END SUBROUTINE p86_init

! *****************************************************************************

  SUBROUTINE p86_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "J. P. Perdew, Phys. Rev. B, 33, 8822 (1986)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew 1986 correlation energy functional"
    END IF

  END SUBROUTINE p86_info

! *****************************************************************************

  SUBROUTINE p86_u ( rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: rs
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,110)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "p86_u", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( rs(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_u", "rs", n )
    
    CALL calc_rs ( rho, rs )
    IF ( calc(0) ) THEN
       CALL p86_u_0 ( rho, rs, grho, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL p86_u_1 ( rho, grho, rs, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL p86_u_2 ( rho, grho, rs, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL p86_u_3 ( rho, grho, rs, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( rs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_u", "rs" )

  END SUBROUTINE p86_u

! *****************************************************************************

  SUBROUTINE p86_p ( rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13a, r13b
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

!   CALL setup_calculation(order,mx,calc,210)
!   mp = SIZE ( pot, 2 )
!   m = ABS(order)
!   mt = mx(m,2)
!   IF ( mt > mp ) CALL stop_program ( "lyp_p", &
!      "Size of POT array too small" )

!   n = SIZE ( rhoa )
!   ALLOCATE ( r13a(n), STAT=error )
!   IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13a", n )
!   ALLOCATE ( r13b(n), STAT=error )
!   IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13b", n )

!   CALL calc_rho13 ( rhoa, r13a )
!   CALL calc_rho13 ( rhob, r13b )
!   IF ( calc(0) ) THEN
!      CALL p86_p_0 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,1) )
!   END IF
!   IF ( calc(1) ) THEN
!      m1=mx(1,1)
!      m2=mx(1,2)
!      CALL p86_p_1 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,m1:m2) )
!   END IF
!   IF ( calc(2) ) THEN
!      m1=mx(2,1)
!      m2=mx(2,2)
!      CALL p86_p_2 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,m1:m2) )
!   END IF
!   IF ( calc(3) ) THEN
!      m1=mx(3,1)
!      m2=mx(3,2)
!      CALL p86_p_3 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,m1:m2) )
!   END IF

!   DEALLOCATE ( r13a, STAT=error )
!   IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13a" )
!   DEALLOCATE ( r13b, STAT=error )
!   IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13b" )

  END SUBROUTINE p86_p

! *****************************************************************************

  SUBROUTINE p86_u_0 ( rho, rs, grho, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, rs, grho
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_0", "Inconsistent array sizes" )


    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN
         pot(ip) = 0.0_dbl
      ELSE
         pot(ip) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_u_0

! *****************************************************************************

  SUBROUTINE p86_u_1 ( rho, grho, rs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, rs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_1", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_u_1

! *****************************************************************************

  SUBROUTINE p86_u_2 ( rho, grho, rs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, rs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_2", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
      ELSE
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_u_2

! *****************************************************************************

  SUBROUTINE p86_u_3 ( rho, grho, rs, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, rs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_3", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      ELSE
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_u_3

! *****************************************************************************

  SUBROUTINE p86_p_0 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip) = 0.0_dbl
      ELSE
         pot(ip) = 0.0_dbl
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_p_0

! *****************************************************************************

  SUBROUTINE p86_p_1 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_1", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      ELSE
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_p_1

! *****************************************************************************

  SUBROUTINE p86_p_2 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_2", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
      ELSE
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,4) = 0.0_dbl
         pot(ip,5) = 0.0_dbl
         pot(ip,6) = 0.0_dbl
      ELSE
         pot(ip,4) = 0.0_dbl
         pot(ip,5) = 0.0_dbl
         pot(ip,6) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_p_2

! *****************************************************************************

  SUBROUTINE p86_p_3 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_3", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      ELSE
         pot(ip,1) = 0.0_dbl
         pot(ip,2) = 0.0_dbl
         pot(ip,3) = 0.0_dbl
         pot(ip,4) = 0.0_dbl
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,5) = 0.0_dbl
         pot(ip,6) = 0.0_dbl
         pot(ip,7) = 0.0_dbl
         pot(ip,8) = 0.0_dbl
      ELSE
         pot(ip,5) = 0.0_dbl
         pot(ip,6) = 0.0_dbl
         pot(ip,7) = 0.0_dbl
         pot(ip,8) = 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_p_3

! *****************************************************************************

END MODULE perdew86_functional

! *****************************************************************************
