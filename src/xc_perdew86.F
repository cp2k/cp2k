!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_perdew86 [1.0] *
!!
!!   NAME
!!     xc_perdew86
!!
!!   FUNCTION
!!     Calculate the Perdew Correlation from 1986
!!
!!   AUTHOR
!!     JGH (03.03.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!******************************************************************************

MODULE xc_perdew86

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE xc_derivative_types, ONLY: xc_derivative_set_type
  USE xc_functionals_utilities, ONLY: set_util, setup_calculation, calc_rs
  USE xc_perdew_zunger, ONLY : perdew_zunger
  USE termination, ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13, &
                          f76 = 7._dbl/6._dbl, &
                          frs = 1.6119919540164696407_dbl, &
                          fpe = 0.19199566167376364_dbl

  PUBLIC :: perdew86

  REAL(dbl) :: eps_rho
  LOGICAL :: debug_flag

  REAL(dbl), PARAMETER :: a = 0.023266_dbl, &
                          b = 7.389e-6_dbl, &
                          c = 8.723_dbl, &
                          d = 0.472_dbl, &
                          pc1 = 0.001667_dbl, &
                          pc2 = 0.002568_dbl, &
                          pci = pc1+pc2

! *****************************************************************************

  INTERFACE perdew86
    MODULE PROCEDURE p86_u, p86_p, p86_init, p86_info
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE p86_init ( cutoff, debug )

    REAL(dbl), INTENT(IN) :: cutoff
    LOGICAL, OPTIONAL, INTENT(IN) :: debug

    eps_rho = cutoff
    CALL set_util ( cutoff )

    CALL perdew_zunger('ORIG',cutoff)

    IF ( PRESENT ( debug ) ) THEN
       debug_flag = debug
    ELSE
       debug_flag = .FALSE.
    ENDIF

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

  SUBROUTINE p86_u ( derivative_set, rho, grho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
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
    IF ( mt > mp ) CALL stop_program ( "p86_u", "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( rs(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_u", "rs", n )
    
    CALL calc_rs ( rho, rs )
    IF ( calc(0) ) THEN
       CALL perdew_zunger ( derivative_set, rho, pot(:,1:1), 0 )
       CALL p86_u_0 ( derivative_set, rho, rs, grho, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL perdew_zunger ( derivative_set, rho, pot(:,m1:m2), -1 )
       CALL p86_u_1 ( derivative_set, rho, grho, rs, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL perdew_zunger ( derivative_set, rho, pot(:,m1:m2), -2 )
       CALL p86_u_2 ( derivative_set, rho, grho, rs, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL perdew_zunger ( derivative_set, rho, pot(:,m1:m2), -3 )
       CALL p86_u_3 ( derivative_set, rho, grho, rs, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( rs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_u", "rs" )

  END SUBROUTINE p86_u

! *****************************************************************************

  SUBROUTINE p86_p ( derivative_set, rhoa, rhob, grhoa, grhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grhoa, grhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: rsa, rsb
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: locpot
    INTEGER :: n, error, m, mp, mt, ip, mx(0:3,2), m1, m2, k

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,211)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "p86_p", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( rsa(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_p", "rsa", n )
    ALLOCATE ( rsb(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_p", "rsb", n )

    CALL calc_rs ( rhoa, rsa )
    CALL calc_rs ( rhob, rsb )
    IF ( calc(0) ) THEN
       CALL p86_p_0 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL p86_p_1 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL p86_p_2 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL p86_p_3 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( rsa, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_p", "rsa" )
    DEALLOCATE ( rsb, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "p86_p", "rsb" )

  END SUBROUTINE p86_p

! *****************************************************************************

  SUBROUTINE p86_u_0 ( derivative_set, rho, rs, grho, pot )

     TYPE(xc_derivative_set_type), POINTER :: derivative_set
   REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, rs, grho
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: r, cr, phi, ep, g, or, x
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_0", "Inconsistent array sizes" )


!$omp parallel do private(ip,g,r,x,or,cr,phi,ep)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1._dbl/rho(ip)
         cr= pc1 + (pc2+a*r+b*r*r)/(1._dbl+c*r+d*r*r+1.e4_dbl*b*r*r*r)
         phi = fpe*pci/cr * g * SQRT(x)*or
         ep = EXP(-phi)
         pot(ip) = pot(ip) + x*or*g*g*cr*ep
      END IF

    END DO

  END SUBROUTINE p86_u_0

! *****************************************************************************

  SUBROUTINE p86_u_1 ( derivative_set, rho, grho, rs, pot )

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, rs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: x, r, cr, phi, ep, g, or, p, q, dp, dq, dphir, dcr, dphig, ff
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_1", "Inconsistent array sizes" )

!$omp parallel do private(ip,g,r,x,or,p,dp,q,dq,cr,dcr,dphig,phi,dphir,ep,ff)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1._dbl/rho(ip)
         p = pc2+a*r+b*r*r
         dp = a+2._dbl*b*r
         q = 1._dbl+c*r+d*r*r+1.e4_dbl*b*r*r*r
         dq = c+2._dbl*d*r+3.e4_dbl*b*r*r
         cr= pc1 + p/q
         dcr = (dp*q-p*dq)/(q*q) * (-f13*r*or)
         dphig = fpe*pci/cr * SQRT(x)*or
         phi = dphig * g
         dphir = -phi*(dcr/cr + f76*or)
         ep = EXP(-phi)
         ff = x*or*g*ep
         pot(ip,1) = pot(ip,1) + ff*g*dcr-ff*g*cr*dphir-ff*g*cr*f43*or
         pot(ip,2) = pot(ip,2) + ff*cr*(2._dbl-g*dphig)
      END IF

    END DO

  END SUBROUTINE p86_u_1

! *****************************************************************************

  SUBROUTINE p86_u_2 ( derivative_set, rho, grho, rs, pot )

     TYPE(xc_derivative_set_type), POINTER :: derivative_set
   REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, rs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: x, r, cr, phi, ep, g, or, p, q, dp, dq, dphir, dcr, dphig
    REAL(dbl) :: dphigr, d2phir, d2cr, d2p, d2q
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_2", "Inconsistent array sizes" )

!$omp parallel do private(ip,x,r,cr,phi,ep,g,or,p,q,dp,dq,dphir,dcr) &
!$omp             private(dphig,dphigr,d2phir,d2cr,d2p,d2q)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1._dbl/rho(ip)
         p = pc2+a*r+b*r*r
         dp = a+2._dbl*b*r
         d2p = 2._dbl*b
         q = 1._dbl+c*r+d*r*r+1.e4_dbl*b*r*r*r
         dq = c+2._dbl*d*r+3.e4_dbl*b*r*r
         d2q = 2._dbl*d+6.e4*b*r
         cr= pc1 + p/q
         dcr = (dp*q-p*dq)/(q*q) * (-f13*r*or)
         d2cr = (d2p*q*q-p*q*d2q-2*dp*dq*q+2*p*dq*dq)/(q*q*q)*(f13*r*or)**2 + &
                (dp*q-p*dq)/(q*q) * f13*f43*r*or*or
         dphig = fpe*pci/cr * SQRT(x)*or
         phi = dphig * g
         dphir = -phi*(dcr/cr + f76*or)
         d2phir = -dphir*(dcr/cr + f76*or)-&
                   phi*((d2cr*cr-dcr*dcr)/(cr*cr) - f76*or*or)
         dphigr = -dphig*(dcr/cr + f76*or)
         ep = EXP(-phi)
         pot(ip,1) = pot(ip,1) + x*or*ep*g*g*(-f43*or*dcr+d2cr-dcr*dphir + &
                     f43*or*cr*dphir-dcr*dphir-cr*d2phir+cr*dphir*dphir + &
                     f43*or*(7.*f13*or*cr-dcr+cr*dphir))
         pot(ip,2) = pot(ip,2) + x*or*ep*g*(-2*f43*cr*or+2*dcr-2*cr*dphir+f43*or*g*cr*dphig-&
                     g*dcr*dphig+g*cr*dphir*dphig-g*cr*dphigr)
         pot(ip,3) = pot(ip,3) + x*or*ep*cr*(2._dbl-4._dbl*g*dphig+g*g*dphig*dphig)
      END IF

    END DO

  END SUBROUTINE p86_u_2

! *****************************************************************************

  SUBROUTINE p86_u_3 ( derivative_set, rho, grho, rs, pot )

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, grho, rs
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: x, r, cr, phi, ep, g, or, p, q, dp, dq, dphir, dcr, dphig
    REAL(dbl) :: dphigr, d2phir, d3phir, d2cr, d3cr, d2p, d2q, d2phirg, d3q
    REAL(dbl) :: pq, dpq, d2pq, d3pq, oz, dz, d2z, d3z
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "p86_u_3", "Inconsistent array sizes" )

!$omp parallel do private(ip,x, r, cr, phi, ep, g, or, p, q, dp, dq, dphir, dcr, dphig) &
!$omp             private(dphigr, d2phir, d3phir, d2cr, d3cr, d2p, d2q, d2phirg, d3q) &
!$omp             private(pq, dpq, d2pq, d3pq, oz, dz, d2z, d3z)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1._dbl/rho(ip)
         p = pc2+a*r+b*r*r
         dp = a+2._dbl*b*r
         d2p = 2._dbl*b
         q = 1._dbl+c*r+d*r*r+1.e4_dbl*b*r*r*r
         dq = c+2._dbl*d*r+3.e4_dbl*b*r*r
         d2q = 2._dbl*d+6.e4*b*r
         d3q = 6.e4*b
         pq = p/q
         dpq = (dp*q-p*dq)/(q*q)
         d2pq = (d2p*q*q-2*dp*dq*q+2*p*dq*dq-p*d2q*q)/(q*q*q)
         d3pq = -(3*d2p*dq*q*q-6*dp*dq*dq*q+3*dp*d2q*q*q+6*p*dq*dq*dq-6*p*dq*d2q*q &
                  +p*d3q*q*q)/(q*q*q*q)
         cr= pc1 + pq
         dcr = dpq*(-f13*r*or)
         d2cr = d2pq*f13*f13*r*r*or*or + dpq*f13*f43*r*or*or
         d3cr = d3pq*(-f13*r*or)**3 + 3*d2pq*(-f13*f13*f43*r*r*or*or*or) + &
                dpq*(-f13*f43*f13*7*r*or*or*or)
         oz = SQRT(x)*or/cr
         dz = dcr/cr+f76*or
         d2z = d2cr/cr+2*f76*dcr/cr*or+f76/6.*or*or
         d3z = d3cr/cr+3*f76*d2cr/cr*or+3*f76/6.*dcr/cr*or*or-5*f76/36.*or*or*or
         dphig = fpe*pci*oz
         phi = dphig * g
         dphir = -phi*dz
         dphigr = -dphig*dz
         d2phir = -phi*(d2z-2*dz*dz)
         d3phir = -phi*(d3z-6*d2z*dz+6*dz*dz*dz)
         d2phirg = -dphigr*dz-&
                   dphig*((d2cr*cr-dcr*dcr)/(cr*cr) - f76*or*or)
         ep = EXP(-phi)
         pot(ip,1) = pot(ip,1) + g*g*x*or*ep*(-280./27.*or*or*or*cr+3*28./9.*or*or*dcr+&
                     3*28./9.*or*or*cr*(-dphir)-4*or*d2cr-8*or*dcr*(-dphir)-&
                     4*or*cr*(-d2phir+dphir*dphir)+d3cr+3*d2cr*(-dphir)+&
                     3*dcr*(-d2phir+dphir*dphir)+cr*(-d3phir+3*dphir*d2phir-&
                     dphir**3))
         pot(ip,2) = pot(ip,2) + 2.*x*or*ep*g*(-f43*or*dcr+d2cr-dcr*dphir + &
                     f43*or*cr*dphir-dcr*dphir-cr*d2phir+cr*dphir*dphir + &
                     f43*or*(7.*f13*or*cr-dcr+cr*dphir)) - &
                     dphig*x*or*ep*g*g*(-f43*or*dcr+d2cr-dcr*dphir + &
                     f43*or*cr*dphir-dcr*dphir-cr*d2phir+cr*dphir*dphir + &
                     f43*or*(7.*f13*or*cr-dcr+cr*dphir)) + &
                     x*or*ep*g*g*(-dcr*dphigr+f43*or*cr*dphigr-dcr*dphigr-cr*d2phirg+&
                     2.*cr*dphigr*dphir+f43*or*cr*dphigr)
         pot(ip,3) = pot(ip,3) + x*or*ep*(-2*f43*cr*or+2*dcr-2*cr*dphir+f43*or*g*cr*dphig-&
                     g*dcr*dphig+g*cr*dphir*dphig-g*cr*dphigr) + &
                     x*or*ep*g*(-2*cr*dphigr+f43*or*cr*dphig-&
                     dcr*dphig+cr*dphir*dphig+g*cr*dphigr*dphig-cr*dphigr)-&
                     x*or*ep*g*dphig*(-2*f43*cr*or+2*dcr-2*cr*dphir+f43*or*g*cr*dphig-&
                     g*dcr*dphig+g*cr*dphir*dphig-g*cr*dphigr)
         pot(ip,4) = pot(ip,4) + x*or*ep*cr*dphig*(-6._dbl+6._dbl*g*dphig-g*g*dphig*dphig)
      END IF

    END DO

  END SUBROUTINE p86_u_3

! *****************************************************************************

  SUBROUTINE p86_p_0 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot )

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            rsa, rsb
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_0", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program ( "p86_p_0", "Routine not tested" )

    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + 0.0_dbl
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + 0.0_dbl
      END IF

    END DO

  END SUBROUTINE p86_p_0

! *****************************************************************************

  SUBROUTINE p86_p_1 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot )

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            rsa, rsb
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_1", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program ( "p86_p_1", "Routine not tested" )

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

  END SUBROUTINE p86_p_1

! *****************************************************************************

  SUBROUTINE p86_p_2 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot )

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            rsa, rsb
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_2", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program ( "p86_p_2", "Routine not tested" )

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

  END SUBROUTINE p86_p_2

! *****************************************************************************

  SUBROUTINE p86_p_3 ( derivative_set, rhoa, rhob, grhoa, grhob, rsa, rsb, pot )

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, grhoa, grhob, &
                                            rsa, rsb
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "p86_p_3", "Inconsistent array sizes" )

    IF ( .NOT. debug_flag ) CALL stop_program ( "p86_p_3", "Routine not tested" )

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

  END SUBROUTINE p86_p_3

! *****************************************************************************

END MODULE xc_perdew86

! *****************************************************************************
