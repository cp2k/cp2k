!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_lyp [1.0] *
!!
!!   NAME
!!     xc_lyp
!!
!!   FUNCTION
!!     Calculate the Lee-Yang-Parr Functional
!!
!!   AUTHOR
!!     JGH (27.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!******************************************************************************

MODULE xc_lyp

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE xc_derivative_types, ONLY: xc_derivative_set_type
  USE xc_functionals_utilities, ONLY: set_util, calc_rho13, setup_calculation
  USE termination, ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: lyp

  REAL(dbl) :: cf, flda, flsd, fvw
  REAL(dbl) :: eps_rho
  LOGICAL :: debug_flag

  REAL(dbl), PARAMETER :: a = 0.04918_dbl, &
                          b = 0.132_dbl, &
                          c = 0.2533_dbl, &
                          d = 0.349_dbl

! *****************************************************************************

  INTERFACE lyp
    MODULE PROCEDURE lyp_u, lyp_p, lyp_init, lyp_info
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE lyp_init ( cutoff, debug )

    REAL(dbl), INTENT(IN) :: cutoff
    LOGICAL, OPTIONAL, INTENT(IN) :: debug

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23

    IF ( PRESENT ( debug ) ) THEN
       debug_flag = debug
    ELSE
       debug_flag = .FALSE.
    END IF

  END SUBROUTINE lyp_init

! *****************************************************************************

  SUBROUTINE lyp_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "C. Lee, W. Yang, R.G. Parr, Phys. Rev. B, 37, 785 (1988)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Lee-Yang-Parr correlation energy functional"
    END IF

  END SUBROUTINE lyp_info

! *****************************************************************************

  SUBROUTINE lyp_u ( derivative_set, rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,110)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "lyp_u", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "lyp_u", "r13", n )
    
    CALL calc_rho13 ( rho, r13 )
    IF ( calc(0) ) THEN
       CALL lyp_u_0 ( rho, grho, r13, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL lyp_u_1 ( rho, grho, r13, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL lyp_u_2 ( rho, grho, r13, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL lyp_u_3 ( rho, grho, r13, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "lyp_u", "r13" )

  END SUBROUTINE lyp_u

! *****************************************************************************

  SUBROUTINE lyp_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13a, r13b
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,211)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "lyp_p", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13b", n )

    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    IF ( calc(0) ) THEN
       CALL lyp_p_0 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL lyp_p_1 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL lyp_p_2 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL lyp_p_3 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( r13a, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13a" )
    DEALLOCATE ( r13b, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "lyp_p", "r13b" )

  END SUBROUTINE lyp_p

! *****************************************************************************

  SUBROUTINE lyp_u_0 ( rho, grho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13, grho
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: x, odx, ecx, lda, ax, bx, g, f24
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "lyp_u_0", "Inconsistent array sizes" )

    f24 = 1._dbl/24._dbl

!$omp parallel do private(ip,x,g,odx,ecx,lda,ax,bx)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         x = 1._dbl/r13(ip)
         g = grho(ip)
         odx = 1._dbl/(1._dbl+d*x)
         ecx = EXP(-c*x)
         lda = -a*odx * (1._dbl+b*cf*ecx)
         ax  = a*b*x**5*ecx*odx
         bx  = 1._dbl+7*f13*(c*x+d*x*odx)
         pot(ip) = pot(ip) + lda*rho(ip) + f24*g*g*ax*bx
      END IF

    END DO

  END SUBROUTINE lyp_u_0

! *****************************************************************************

  SUBROUTINE lyp_u_1 ( rho, grho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: x, g, odx, ecx, lda, dlda, ax, bx, dax, dbx, f24
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "lyp_u_1", "Inconsistent array sizes" )

    f24 = 1._dbl/24._dbl

!$omp parallel do private(ip,x,g,odx,ecx,lda,dlda,ax,dax,bx,dbx)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         x = 1._dbl/r13(ip)
         g = grho(ip)
         odx = 1._dbl/(1._dbl+d*x)
         ecx = EXP(-c*x)
         lda = -a*odx * (1._dbl+b*cf*ecx)
         dlda = a*b*cf*c*ecx*odx+a*(1+b*cf*ecx)*odx**2*d
         ax  = a*b*x**5*ecx*odx
         dax = a*b*(5._dbl*x**4*ecx*odx-x**5*c*ecx*odx-x**5*ecx*odx*odx*d)
         bx  = 1._dbl+7*f13*(c*x+d*x*odx)
         dbx = 7*f13*(c+d*odx-d**2*x*odx**2)
         pot(ip,1) = pot(ip,1) + lda - f13*dlda*x - &
                     f24*f13*x**4*g*g*(dax*bx+ax*dbx)
         pot(ip,2) = pot(ip,2) + 2*f24*g*ax*bx
      END IF

    END DO

  END SUBROUTINE lyp_u_1

! *****************************************************************************

  SUBROUTINE lyp_u_2 ( rho, grho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: odx, ecx, x, g, lda, dlda, ax, bx, dax, dbx, dab, d2ab
    REAL(dbl) :: d2lda, d2ax, d2bx, f24
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "lyp_u_2", "Inconsistent array sizes" )

    f24 = 1._dbl/24._dbl

!$omp parallel do private(ip,odx,ecx,x,g,lda,dlda,ax,bx,dax,dbx,dab,d2ab,d2lda,d2ax,d2bx)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         x=1._dbl/r13(ip)
         g=grho(ip)
         odx = 1._dbl/(1._dbl+d*x)
         ecx = EXP(-c*x)
         lda = -a*odx * (1._dbl+b*cf*ecx)
         dlda = a*odx*odx*(d+d*b*cf*ecx+b*cf*c*ecx+b*c*cf*d*x*ecx)
         d2lda = -2*a*odx*odx*odx*(1+b*cf*ecx)*d*d-2*a*odx*odx*b*cf*c*ecx*d-&
                 a*odx*b*cf*c*c*ecx
         ax  = a*b*x**5*ecx*odx
         dax = a*b*(5._dbl*x**4*ecx*odx-x**5*c*ecx*odx-x**5*ecx*odx*odx*d)
         d2ax = 20._dbl*a*b*x**3*ecx*odx-10._dbl*a*b*x**4*c*ecx*odx-&
                10._dbl*a*b*x**4*ecx*odx*odx*d+a*b*x**5*c*c*ecx*odx+&
                2._dbl*a*b*x**5*c*ecx*odx*odx*d+&
                2._dbl*a*b*x**5*ecx*odx*odx*odx*d*d
         bx  = 1._dbl+7._dbl/3._dbl*(c*x+d*x*odx)
         dbx = 7./3.*(c+d*odx-d**2*x*odx**2)
         d2bx = 14./3.*d*d*(-odx**2+d*x*odx**3)
         dab = dax*bx+ax*dbx
         d2ab = d2ax*bx+2*dax*dbx+ax*d2bx
         pot(ip,1) = pot(ip,1) + f13*x**4*(-f23*dlda+f13*x*d2lda+&
                     f24*f43*x**3*g*g*dab+f24*f13*x**4*g*g*d2ab)
         pot(ip,2) = pot(ip,2) - f23*f24*x**4*g*dab
         pot(ip,3) = pot(ip,3) + 2._dbl*f24*ax*bx
      END IF

    END DO

  END SUBROUTINE lyp_u_2

! *****************************************************************************

  SUBROUTINE lyp_u_3 ( rho, grho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, r13
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: odx, ecx, x, g, lda, dlda, ax, bx, dax, dbx, dab, d2ab, d3ab
    REAL(dbl) :: d2lda, d2ax, d2bx, d3lda, d3ax, d3bx, f24
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "lyp_u_3", "Inconsistent array sizes" )

    f24 = 1._dbl/24._dbl

!$omp parallel do private(ip,odx,ecx,x,g,lda,dlda,ax,bx,dax,dbx,dab,d2ab,d3ab) &
!$omp             private(d2lda,d2ax,d2bx,d3lda,d3ax,d3bx)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         x=1._dbl/r13(ip)
         g=grho(ip)
         odx = 1._dbl/(1._dbl+d*x)
         ecx = EXP(-c*x)
         lda = -a*odx * (1._dbl+b*cf*ecx)
         dlda = a*odx*odx*(d+d*b*cf*ecx+b*cf*c*ecx+b*c*cf*d*x*ecx)
         d2lda = -2*a*odx*odx*odx*(1+b*cf*ecx)*d*d-2*a*odx*odx*b*cf*c*ecx*d-&
                 a*odx*b*cf*c*c*ecx
         d3lda = 6*a*odx**4*(1+b*cf*ecx)*d**3+6*a*odx**3*b*cf*c*ecx*d**2+&
                 3*a*odx**2*b*cf*c**2*ecx*d+a*odx*b*cf*c**3*ecx
         ax  = a*b*x**5*ecx*odx
         dax = a*b*(5._dbl*x**4*ecx*odx-x**5*c*ecx*odx-x**5*ecx*odx*odx*d)
         d2ax = 20._dbl*a*b*x**3*ecx*odx-10._dbl*a*b*x**4*c*ecx*odx-&
                10._dbl*a*b*x**4*ecx*odx*odx*d+a*b*x**5*c*c*ecx*odx+&
                2._dbl*a*b*x**5*c*ecx*odx*odx*d+&
                2._dbl*a*b*x**5*ecx*odx*odx*odx*d*d
         d3ax = 60*a*b*x**2*ecx*odx-60*a*b*x**3*c*ecx*odx-60*a*b*x**3*ecx*odx**2*d+&
                15*a*b*x**4*c**2*ecx*odx+30*a*b*x**4*c*ecx*odx**2*d+&
                30*a*b*x**4*ecx*odx**3*d**2-a*b*x**5*c**3*ecx*odx-&
                3*a*b*x**5*c**2*ecx*odx**2*d-6*a*b*x**5*c*ecx*odx**3*d**2-&
                6*a*b*x**5*ecx*odx**4*d**3
         bx  = 1._dbl+7._dbl/3._dbl*(c*x+d*x*odx)
         dbx = 7./3.*(c+d*odx-d**2*x*odx**2)
         d2bx = 14./3.*d*d*(-odx**2+d*x*odx**3)
         d3bx = 14*d**3*(odx**3-d*x*odx**4)
         dab = dax*bx+ax*dbx
         d2ab = d2ax*bx+2*dax*dbx+ax*d2bx
         d3ab = d3ax*bx+3*d2ax*dbx+3*dax*d2bx+ax*d3bx
         pot(ip,1) = pot(ip,1) + x**7/27*(8*dlda - 3*x*d2lda - x*x*d3lda) &
                     + f24/27*g*g* ( -28*x**10*dab - 12*x**11*d2ab - x**12*d3ab )
         pot(ip,2) = pot(ip,2) + f24*f23*g*x**7*(f43*dab+f13*x*d2ab)
         pot(ip,3) = pot(ip,3) - f24*f23*x**4*dab
         pot(ip,4) = pot(ip,4) + 0.0_dbl
      END IF

    END DO

  END SUBROUTINE lyp_u_3

! *****************************************************************************

  SUBROUTINE lyp_p_0 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "lyp_p_0", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program("lyp_p_0","Routine not tested")

    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + 0.0_dbl
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + 0.0_dbl
      END IF

    END DO

  END SUBROUTINE lyp_p_0

! *****************************************************************************

  SUBROUTINE lyp_p_1 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "lyp_p_1", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program("lyp_p_1","Routine not tested")

    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,3) = pot(ip,3) + 0.0_dbl
         pot(ip,4) = pot(ip,4) + 0.0_dbl
      END IF

    END DO

  END SUBROUTINE lyp_p_1

! *****************************************************************************

  SUBROUTINE lyp_p_2 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "lyp_p_2", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program("lyp_p_2","Routine not tested")

    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl
         pot(ip,3) = pot(ip,3) + 0.0_dbl
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,4) = pot(ip,4) + 0.0_dbl
         pot(ip,5) = pot(ip,5) + 0.0_dbl
         pot(ip,6) = pot(ip,6) + 0.0_dbl
      END IF

    END DO

  END SUBROUTINE lyp_p_2

! *****************************************************************************

  SUBROUTINE lyp_p_3 ( rhoa, rhob, grhoa, grhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "lyp_p_3", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program("lyp_p_3","Routine not tested")

    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl
         pot(ip,3) = pot(ip,3) + 0.0_dbl
         pot(ip,4) = pot(ip,4) + 0.0_dbl
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,5) = pot(ip,5) + 0.0_dbl
         pot(ip,6) = pot(ip,6) + 0.0_dbl
         pot(ip,7) = pot(ip,7) + 0.0_dbl
         pot(ip,8) = pot(ip,8) + 0.0_dbl
      END IF

    END DO

  END SUBROUTINE lyp_p_3

! *****************************************************************************

END MODULE xc_lyp

! *****************************************************************************
