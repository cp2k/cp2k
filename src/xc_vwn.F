!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_vwn [1.0] *
!!
!!   NAME
!!     xc_vwn
!!
!!   FUNCTION
!!     Calculate the LDA functional according to Vosk, Wilk and Nusair
!!     Literature: S. H. Vosko, L. Wilk and M. Nusair,
!!                 Can. J. Phys. 58, 1200 (1980)
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3;
!!                              LSD 0; a  b; aa ab bb; aaa aab abb bbb;
!!
!!   AUTHOR
!!     JGH (26.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!*** *************************************************************************

MODULE xc_vwn
! *****************************************************************************
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals_utilities,        ONLY: calc_fx,&
                                             calc_srs,&
                                             set_util,&
                                             setup_calculation
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13

  REAL(dbl), PARAMETER :: a  = 0.0310907_dbl,&
                          b  = 3.72744_dbl,&
                          c  = 12.9352_dbl,&
                          x0 = -0.10498_dbl

  PUBLIC :: vwn

  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE vwn
    MODULE PROCEDURE vwn_lda, vwn_lsd, vwn_init, vwn_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE vwn_init ( cutoff )

    REAL(dbl), INTENT(IN) :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

  END SUBROUTINE vwn_init

! *****************************************************************************

  SUBROUTINE vwn_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
       reference = "S. H. Vosko, L. Wilk and M. Nusair,"// &
                   " Can. J. Phys. 58, 1200 (1980)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Vosko-Wilk-Nusair Functional"
    END IF

  END SUBROUTINE vwn_info

! *****************************************************************************

  SUBROUTINE vwn_lda ( derivative_set, rho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: x
    INTEGER :: n, error, m, mp, mt, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,100)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "vwn_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( x(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "vwn_lda", "x", n )
    
    CALL calc_srs ( rho, x )
    IF ( calc(0) .AND. calc(1) ) THEN
       CALL vwn_lda_01 ( rho, x, pot(:,1:2) )
    ELSE IF ( calc(0) ) THEN
       CALL vwn_lda_0 ( rho, x, pot(:,1) )
    ELSE IF ( calc(1) ) THEN
       m1=mx(1,1)
       CALL vwn_lda_1 ( rho, x, pot(:,m1) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       CALL vwn_lda_2 ( rho, x, pot(:,m1) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       CALL vwn_lda_3 ( rho, x, pot(:,m1) )
    END IF

    DEALLOCATE ( x, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "vwn_lda", "x" )

  END SUBROUTINE vwn_lda

! *****************************************************************************

  SUBROUTINE vwn_lsd ( derivative_set, rhoa, rhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: x
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: fx
    INTEGER :: n, error, m, mp, mt, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,201)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "vwn_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( x(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "vwn_lsd", "x", n )
    m = ABS ( order ) + 1
    ALLOCATE ( fx(n,m), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "vwn_lsd", "fx", n*m )

    fx(1:n,1) = rhoa(1:n) + rhob(1:n)
    CALL calc_srs ( fx(:,1), x )
    CALL calc_fx ( rhoa, rhob, fx, m-1 )

    IF ( calc(0) .AND. calc(1) ) THEN
       m1=mx(0,1)
       m2=mx(1,2)
       CALL vwn_lsd_01 ( rhoa, rhob, x, fx, pot(:,m1:m2) )
    ELSE IF ( calc(0) ) THEN
       CALL vwn_lsd_0 ( rhoa, rhob, x, fx, pot(:,1) )
    ELSE IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL vwn_lsd_1 ( rhoa, rhob, x, fx, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL vwn_lsd_2 ( rhoa, rhob, x, fx, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL vwn_lsd_3 ( rhoa, rhob, x, fx, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( fx, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "vwn_lsd", "fx" )
    DEALLOCATE ( x, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "vwn_lsd", "x" )

  END SUBROUTINE vwn_lsd

! *****************************************************************************

  SUBROUTINE vwn_lda_0 ( rho, x, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, x
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: px, q, at, ln1, ln2, xb, dpx
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "vwn_lda_0", "Inconsistent array sizes" )

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b

!$omp parallel do private(ip,px,dpx,at,ln1,ln2)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         ln1 = LOG(x(ip)*x(ip)/px)
         ln2 = LOG((x(ip)-x0)**2/px)
         pot(ip) = pot(ip) + a * ( ln1 + b*at - b*x0/px* ( ln2 + xb * at ) ) * rho(ip)
      END IF

    END DO

  END SUBROUTINE vwn_lda_0

! *****************************************************************************

  SUBROUTINE vwn_lda_1 ( rho, x, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, x
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: px, q, at, dat, ln1, ln2, xb, ex, pa, dpx, dex, dln1, dln2
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "vwn_lda_1", "Inconsistent array sizes" )

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b
!$omp parallel do private(ip,px,dpx,at,pa,dat,ln1,dln1,ln2,dln2,ex,dex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2._dbl*c)/(x(ip)*px)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (b*x(ip)+2._dbl*c+2._dbl*x0*x(ip)+x0*b)/((x(ip)-x0)*px)
         ex = a * ( ln1 + b*at - b*x0/px* ( ln2 + xb * at ) )
         dex = a * ( dln1 + b*dat + b*x0*dpx/(px*px) * ( ln2 + xb * at ) - &
               b*x0/px* ( dln2 + xb * dat ) )
         pot(ip) = pot(ip) + ex - x(ip)*dex/6._dbl
      END IF

    END DO

  END SUBROUTINE vwn_lda_1

! *****************************************************************************

  SUBROUTINE vwn_lda_01 ( rho, x, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, x
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: px, dpx, q, at, dat, ln1, ln2, xb, ex, pa, dex, dln1, dln2
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "vwn_lda_01", "Inconsistent array sizes" )

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b
!$omp parallel do private(ip,px,dpx,at,pa,dat,ln1,dln1,ln2,dln2,ex,dex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2._dbl*c)/(x(ip)*px)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (x(ip)*(b+2._dbl*x0)+2._dbl*c+x0*b)/((x(ip)-x0)*px)
         ex = a * ( ln1 + b*at - b*x0/px* ( ln2 + xb * at ) )
         dex = a * ( dln1 + b*dat + b*x0*dpx/(px*px) * ( ln2 + xb * at ) - &
               b*x0/px* ( dln2 + xb * dat ) )
         pot(ip,1) = pot(ip,1) + ex * rho(ip)
         pot(ip,2) = pot(ip,2) + ex - x(ip)*dex/6._dbl
      END IF

    END DO


  END SUBROUTINE vwn_lda_01

! *****************************************************************************

  SUBROUTINE vwn_lda_2 ( rho, x, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, x
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: px, dpx, q, at, dat, ln1, ln2, xb, dex, pa, d2ex, dln1, dln2
    REAL(dbl) :: d2ln1, d2ln2, d2at, fp, dfp, d2fp
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "vwn_lda_2", "Inconsistent array sizes" )

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b

!$omp parallel do private(ip,px,dpx,at,pa,dat,d2at,ln1,dln1,d2ln1) &
!$omp             private(ln2,dln2,d2ln2,fp,dfp,d2fp,dex,d2ex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         d2at= 16._dbl*dpx/(pa*pa)
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2._dbl*c)/(x(ip)*px)
         d2ln1 = b/(x(ip)*px) - (b*x(ip)+2._dbl*c)/(x(ip)*px)**2 * (px+x(ip)*dpx)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (x(ip)*xb+2._dbl*c+x0*b)/((x(ip)-x0)*px)
         d2ln2 = xb/((x(ip)-x0)*px)-(x(ip)*xb+2._dbl*c+x0*b)/((x(ip)-x0)*px)**2 &
                 * (px+(x(ip)-x0)*dpx)
         fp = -b*x0/px
         dfp = b*x0*dpx/(px*px)
         d2fp = b*x0*2._dbl*(1._dbl-dpx*dpx/px)/(px*px)
         dex = a * ( dln1 + b*dat + dfp * ( ln2 + xb * at ) + fp * ( dln2 + xb * dat ) )
         d2ex= a * ( d2ln1 + b*d2at + d2fp * ( ln2 + xb * at ) + &
               2._dbl*dfp * ( dln2 + xb * dat ) + fp * ( d2ln2 + xb * d2at ) )
         pot(ip) = pot(ip) + x(ip)/(36._dbl*rho(ip))*(x(ip)*d2ex-5._dbl*dex)
      END IF

    END DO

  END SUBROUTINE vwn_lda_2

! *****************************************************************************

  SUBROUTINE vwn_lda_3 ( rho, x, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, x
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: px, dpx, q, at, dat, ln1, ln2, xb, dex, pa, d2ex, dln1, dln2
    REAL(dbl) :: d2ln1, d2ln2, d2at, d3ex, d3at, d3ln1, d3ln2, fp, dfp, d2fp
    REAL(dbl) :: d3fp, ax, bx, dbx, d2bx, cx, dx, ddx, d2dx
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "vwn_lda_2", "Inconsistent array sizes" )

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b

!$omp parallel do private(ip,px,dpx,at,pa,dat,d2at,d3at,ln1,ax,bx) &
!$omp             private(dbx,d2bx,dln1,d2ln1,d3ln1,ln2,cx,dx,ddx,d2dx,dln2,d2ln2,d3ln2) &
!$omp             private(fp,dfp,d2fp,d3fp,dex,d2ex,d3ex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         d2at= 16._dbl*dpx/(pa*pa)
         d3at= 32._dbl/(pa*pa)*(1._dbl-4._dbl*dpx*dpx/pa)
         ln1 = LOG(x(ip)*x(ip)/px)
         ax = b*x(ip)+2._dbl*c
         bx = x(ip)*px
         dbx = px + x(ip)*dpx
         d2bx = 2._dbl*(dpx+x(ip))
         dln1 = ax/bx
         d2ln1 = (b*bx-ax*dbx)/(bx*bx)
         d3ln1 = -ax*d2bx/(bx*bx) - 2._dbl*d2ln1*dbx/bx
         ln2 = LOG((x(ip)-x0)**2/px)
         cx = x(ip)*xb+2._dbl*c+x0*b
         dx = (x(ip)-x0)*px
         ddx = px + (x(ip)-x0)*dpx
         d2dx = 2._dbl*(dpx + (x(ip)-x0))
         dln2 = cx/dx
         d2ln2 = (xb*dx-cx*ddx)/(dx*dx)
         d3ln2 = -cx*d2dx/(dx*dx) - 2._dbl*d2ln2*ddx/dx
         fp = -b*x0/px
         dfp = b*x0*dpx/(px*px)
         d2fp = b*x0*2._dbl*(1._dbl-dpx*dpx/px)/(px*px)
         d3fp = b*x0*6._dbl*dpx*(dpx*dpx/px-2._dbl)/(px*px*px)
         dex = a * ( dln1 + b*dat + dfp * ( ln2 + xb * at ) + fp * ( dln2 + xb * dat ) )
         d2ex= a * ( d2ln1 + b*d2at + d2fp * ( ln2 + xb * at ) + &
               2._dbl*dfp * ( dln2 + xb * dat ) + fp * ( d2ln2 + xb * d2at ) )
         d3ex= a * ( d3ln1 + b*d3at + d3fp * ( ln2 + xb * at ) + &
               3._dbl*d2fp * ( dln2 + xb * dat ) + 3._dbl*dfp * ( d2ln2 + xb * d2at ) + &
               fp * ( d3ln2 + xb * d3at ) )
         pot(ip) = pot(ip) - 7._dbl*x(ip)/(216._dbl*rho(ip)*rho(ip))*(x(ip)*d2ex-5._dbl*dex) - &
                   x(ip)*x(ip)/(216._dbl*rho(ip)*rho(ip))*(x(ip)*d3ex-4._dbl*d2ex)
      END IF

    END DO

  END SUBROUTINE vwn_lda_3

! *****************************************************************************

  SUBROUTINE vwn_lsd_0 ( rhoa, rhob, x, fx, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, x
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, rhoab
!   ---------------------------------------------------------------------------

    CALL stop_program ( "vwn_lsd_0", "Not implemented" )

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "vwn_lsd_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         pot(ip) = pot(ip) + 0.0_dbl

      END IF

    END DO

  END SUBROUTINE vwn_lsd_0

! *****************************************************************************

  SUBROUTINE vwn_lsd_1 ( rhoa, rhob, x, fx, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, x
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
!   ---------------------------------------------------------------------------

    CALL stop_program ( "vwn_lsd_1", "Not implemented" )

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "vwn_lsd_1", "Inconsistent array sizes" )

    DO ip = 1, npoints

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl

      END IF

    END DO

  END SUBROUTINE vwn_lsd_1

! *****************************************************************************

  SUBROUTINE vwn_lsd_01 ( rhoa, rhob, x, fx, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, x
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
!   ---------------------------------------------------------------------------

    CALL stop_program ( "vwn_lsd_01", "Not implemented" )

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "vwn_lsd_01", "Inconsistent array sizes" )

    DO ip = 1, npoints

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl
         pot(ip,3) = pot(ip,3) + 0.0_dbl

      END IF

    END DO

  END SUBROUTINE vwn_lsd_01

! *****************************************************************************

  SUBROUTINE vwn_lsd_2 ( rhoa, rhob, x, fx, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, x
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
    REAL(dbl) :: d2p, d2q, dxp, dxq, drr, dxr, dxx, xt, yt, or
!   ---------------------------------------------------------------------------

    CALL stop_program ( "vwn_lsd_2", "Not implemented" )

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "vwn_lsd_2", "Inconsistent array sizes" )

    DO ip = 1, npoints

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl
         pot(ip,3) = pot(ip,3) + 0.0_dbl

      END IF

    END DO

  END SUBROUTINE vwn_lsd_2

! *****************************************************************************

  SUBROUTINE vwn_lsd_3 ( rhoa, rhob, x, fx, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, x
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
    REAL(dbl) :: d2p, d2q, dxp, dxq, drr, dxr, dxx, xt, yt, or
    REAL(dbl) :: d3p, d3q, d2xp, d2xq
!   ---------------------------------------------------------------------------

    CALL stop_program ( "vwn_lsd_3", "Not implemented" )

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "vwn_lsd_3", "Inconsistent array sizes" )

    DO ip = 1, npoints

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         pot(ip,1) = pot(ip,1) + 0.0_dbl
         pot(ip,2) = pot(ip,2) + 0.0_dbl
         pot(ip,3) = pot(ip,3) + 0.0_dbl
         pot(ip,4) = pot(ip,4) + 0.0_dbl

      END IF

    END DO

  END SUBROUTINE vwn_lsd_3

! *****************************************************************************

END MODULE xc_vwn

! *****************************************************************************
