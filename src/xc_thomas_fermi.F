!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_thomas_fermi [1.0] *
!!
!!   NAME
!!     xc_thomas_fermi
!!
!!   FUNCTION
!!     Calculate the Thomas-Fermi kinetic energy functional
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3;
!!                              LSD 0; a  b; aa bb; aaa bbb;
!!
!!   AUTHOR
!!     JGH (18.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!*** *************************************************************************

MODULE xc_thomas_fermi
! *****************************************************************************
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals_utilities,        ONLY: calc_rho13,&
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

  PUBLIC :: thomas_fermi

  REAL(dbl) :: cf, flda, flsd
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE thomas_fermi
    MODULE PROCEDURE thomas_fermi_lda, thomas_fermi_lsd, &
                     thomas_fermi_init, thomas_fermi_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE thomas_fermi_init ( cutoff )

    REAL(dbl), INTENT(IN)                    :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23
    flda = cf
    flsd = flda * 2._dbl**f23

  END SUBROUTINE thomas_fermi_init

! *****************************************************************************

  SUBROUTINE thomas_fermi_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "Thomas-Fermi kinetic energy functional: see Parr and Yang"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Thomas-Fermi kinetic energy functional"
    END IF

  END SUBROUTINE thomas_fermi_info

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda ( derivative_set, rho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: error, m, mp, mt, mx(0:3,2), n
    LOGICAL                                  :: calc(0:4)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: r13

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,100)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "thomas_fermi_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lda", "r13", n )
    
    CALL calc_rho13 ( rho, r13 )
    IF ( calc(0) ) THEN
       CALL thomas_fermi_lda_0 ( rho, r13, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m=mx(1,1)
       CALL thomas_fermi_lda_1 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(2) ) THEN
       m=mx(2,1)
       CALL thomas_fermi_lda_2 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(3) ) THEN
       m=mx(3,1)
       CALL thomas_fermi_lda_3 ( rho, r13, pot(:,m) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lda", "rs" )

  END SUBROUTINE thomas_fermi_lda

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd ( derivative_set, rhoa, rhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: error, m, m1, m2, mp, mt, &
                                                mx(0:3,2), n
    LOGICAL                                  :: calc(0:4)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: r13a, r13b

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,200)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "thomas_fermi_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13b", n )

    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    IF ( calc(0) ) THEN
       CALL thomas_fermi_lsd_0 ( rhoa, rhob, r13a, r13b, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL thomas_fermi_lsd_1 ( rhoa, rhob, r13a, r13b, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL thomas_fermi_lsd_2 ( rhoa, rhob, r13a, r13b, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL thomas_fermi_lsd_3 ( rhoa, rhob, r13a, r13b, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( r13a, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13a" )
    DEALLOCATE ( r13b, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13b" )

  END SUBROUTINE thomas_fermi_lsd

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_0 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip, npoints

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_0", "Inconsistent array sizes" )

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + flda*r13(ip)*r13(ip)*rho(ip)

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_0

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_1 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_1", "Inconsistent array sizes" )

    f = f53 * flda

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + f * r13(ip)*r13(ip)

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_1

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_2 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_2", "Inconsistent array sizes" )

    f = f23 * f53 * flda

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + f / r13(ip)

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_2

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_3 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_3", "Inconsistent array sizes" )

    f = -f13 * f23 * f53 * flda

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + f / ( r13(ip) * rho(ip) )

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_3

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_0 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip, npoints

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_0", "Inconsistent array sizes" )

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13a(ip) * r13a(ip) * rhoa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13b(ip) * r13b(ip) * rhob(ip)
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_0

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_1 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_1", "Inconsistent array sizes" )

    f = f53 * flsd

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + f * r13a(ip) * r13a(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,2) = pot(ip,2) + f * r13b(ip) * r13b(ip)
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_1

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_2 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_2", "Inconsistent array sizes" )

    f = f23 * f53 * flsd

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + f / r13a(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,2) = pot(ip,2) + f / r13b(ip)
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_2

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_3 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip, npoints
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_3", "Inconsistent array sizes" )

    f = -f13 * f23 * f53 * flsd

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + f / ( r13a(ip)*rhoa(ip))
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,2) = pot(ip,2) + f / ( r13b(ip)*rhob(ip))
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_3

! *****************************************************************************

END MODULE xc_thomas_fermi

! *****************************************************************************
