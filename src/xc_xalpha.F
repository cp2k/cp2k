!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_xalpha [1.0] *
!!
!!   NAME
!!     xc_xalpha
!!
!!   FUNCTION
!!     Calculate the local exchange functional
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3; 
!!                              LSD 0; a  b; aa bb; aaa bbb;
!!
!!   AUTHOR
!!     JGH (17.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!*** *************************************************************************

MODULE xc_xalpha
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
                          f43 = 4._dbl*f13

  PUBLIC :: xalpha

  REAL(dbl) :: xparam, flda, flsd
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE xalpha
    MODULE PROCEDURE xalpha_lda_pw, xalpha_lsd_pw, xalpha_init, xalpha_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE xalpha_init ( cutoff, PARAMETER )

    REAL(dbl), INTENT(IN)                    :: cutoff
    REAL(dbl), INTENT(IN), OPTIONAL          :: PARAMETER

    eps_rho = cutoff
    CALL set_util ( cutoff )
    IF ( PRESENT ( PARAMETER ) ) THEN
      xparam = PARAMETER
    ELSE
      xparam = 2._dbl / 3._dbl
    END IF

    flda = -9._dbl/8._dbl * xparam * (3._dbl/pi)**f13
    flsd = flda * 2._dbl**f13

  END SUBROUTINE xalpha_init

! *****************************************************************************

  SUBROUTINE xalpha_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
       WRITE (reference,'(A,F8.4)') &
         "Dirac/Slater local exchange; parameter=",xparam
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       WRITE (shortform,'(A,F8.4)') "Dirac/Slater exchange",xparam
    END IF

  END SUBROUTINE xalpha_info

! *****************************************************************************

  SUBROUTINE xalpha_lda_pw(derivative_set, rho, pot, order)

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rho
    REAL(dbl), DIMENSION(:, :), INTENT(INOUT)  :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: error, m, mp, mt, mx(0:3,2), n
    LOGICAL                                  :: calc(0:4)
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: r13

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,100)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "xalpha_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lda", "r13", n )
    
    CALL calc_rho13 ( rho, r13 )
    IF ( calc(0) ) THEN
       CALL xalpha_lda_0(n, rho, r13, &
                        pot(:,1:1))
    END IF
    IF ( calc(1) ) THEN
       m = mx(1,1)
       CALL xalpha_lda_1(n, rho, r13, &
                         pot(:,m:m))
    END IF
    IF ( calc(2) ) THEN
       m = mx(2,1)
       CALL xalpha_lda_2(n, rho, r13, &
                         pot(:,m:m))
    END IF
    IF ( calc(3) ) THEN
       m = mx(3,1)
       CALL xalpha_lda_3(n, rho, r13, &
                         pot(:,m:m))
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lda", "rs" )

  END SUBROUTINE xalpha_lda_pw

! *****************************************************************************

  SUBROUTINE xalpha_lsd_pw ( derivative_set, rhoa, rhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: rhoa, rhob
    REAL(dbl), DIMENSION(:, :), INTENT(INOUT)  :: pot
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
    IF ( mt > mp ) CALL stop_program ( "xalpha_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13b", n )

    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    IF ( calc(0) ) THEN
       CALL xalpha_lsd_0(n, rhoa, rhob, r13a, r13b, &
                         pot(:,1))
    END IF
    IF ( calc(1) ) THEN
       m1 = mx(1,1)
       m2 = mx(1,2)
       CALL xalpha_lsd_1(n, rhoa, rhob, r13a, r13b, &
                         pot(:,m1:m1), pot(:,m2:m2))
    END IF
    IF ( calc(2) ) THEN
       m1 = mx(2,1)
       m2 = mx(2,2)
       CALL xalpha_lsd_2(n, rhoa, rhob, r13a, r13b, &
                         pot(:,m1:m1), pot(:,m2:m2))
    END IF
    IF ( calc(3) ) THEN
       m1 = mx(3,1)
       m2 = mx(3,2)
       CALL xalpha_lsd_3(n, rhoa, rhob, r13a, r13b, &
                         pot(:,m1:m1), pot(:,m2:m2))
    END IF

    DEALLOCATE ( r13a, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13a" )
    DEALLOCATE ( r13b, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13b" )

  END SUBROUTINE xalpha_lsd_pw

! *****************************************************************************

  SUBROUTINE xalpha_lda_0(n, rho, r13, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: pot

    INTEGER                                  :: ip

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flda*r13(ip)*rho(ip)
      END IF
    END DO

  END SUBROUTINE xalpha_lda_0

! *****************************************************************************

  SUBROUTINE xalpha_lda_1(n, rho, r13, pot)
    
    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: pot

    INTEGER                                  :: ip
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    f = f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + f * r13(ip)
      END IF
    END DO

  END SUBROUTINE xalpha_lda_1

! *****************************************************************************

  SUBROUTINE xalpha_lda_2(n, rho, r13, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: pot

    INTEGER                                  :: ip
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    f = f13 * f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + f * r13(ip) / rho(ip)
      END IF
    END DO

  END SUBROUTINE xalpha_lda_2

! *****************************************************************************

  SUBROUTINE xalpha_lda_3(n, rho, r13, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: pot

    INTEGER                                  :: ip
    REAL(dbl)                                :: f

!   ---------------------------------------------------------------------------

    f = -f23 * f13 * f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + f * r13(ip) / ( rho(ip) * rho(ip) )
      END IF
    END DO

  END SUBROUTINE xalpha_lda_3

! *****************************************************************************

  SUBROUTINE xalpha_lsd_0 ( n, rhoa, rhob, r13a, r13b, pot )

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: pot

    INTEGER                                  :: ip

! number of points in array 
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13a(ip) * rhoa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13b(ip) * rhob(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_0

! *****************************************************************************

  SUBROUTINE xalpha_lsd_1 ( n, rhoa, rhob, r13a, r13b, pota, potb )

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: pota, potb

    INTEGER                                  :: ip
    REAL(dbl)                                :: f

! number of points in array 
!   ---------------------------------------------------------------------------

    f = f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         pota(ip) = pota(ip) + f * r13a(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         potb(ip) = potb(ip) + f * r13b(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_1

! *****************************************************************************

  SUBROUTINE xalpha_lsd_2 ( n, rhoa, rhob, r13a, r13b, potaa, potbb )

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: potaa, potbb

    INTEGER                                  :: ip
    REAL(dbl)                                :: f

! number of points in array 
!   ---------------------------------------------------------------------------

    f = f13 * f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         potaa(ip) = potaa(ip) + f * r13a(ip)/rhoa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         potbb(ip) = potbb(ip) + f * r13b(ip)/rhob(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_2

! *****************************************************************************

  SUBROUTINE xalpha_lsd_3 ( n, rhoa, rhob, r13a, r13b, potaaa, potbbb )

    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(*), INTENT(INOUT)     :: potaaa, potbbb

    INTEGER                                  :: ip
    REAL(dbl)                                :: f

! number of points in array 
!   ---------------------------------------------------------------------------

    f = -f23 * f13 * f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         potaaa(ip) = potaaa(ip) + f * r13a(ip)/(rhoa(ip)*rhoa(ip))
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         potbbb(ip) = potbbb(ip) + f * r13b(ip)/(rhob(ip)*rhob(ip))
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_3

! *****************************************************************************

END MODULE xc_xalpha

! *****************************************************************************
