!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_pade [1.0] *
!!
!!   NAME
!!     xc_pade
!!
!!   FUNCTION
!!     Calculate the LDA functional in the Pade approximation
!!     Literature: S. Goedecker, M. Teter and J. Hutter,
!!                 Phys. Rev. B 54, 1703 (1996)
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3;
!!                              LSD 0; a  b; aa ab bb; aaa aab abb bbb;
!!
!!   AUTHOR
!!     JGH (15.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!******************************************************************************

MODULE xc_pade
! *****************************************************************************
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals_utilities,        ONLY: calc_fx,&
                                             calc_rs,&
                                             set_util,&
                                             setup_calculation
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13

  REAL(dbl), PARAMETER :: a0 = 0.4581652932831429E+0_dbl,&
                          a1 = 0.2217058676663745E+1_dbl,&
                          a2 = 0.7405551735357053E+0_dbl,&
                          a3 = 0.1968227878617998E-1_dbl,&
                          b1 = 1.0000000000000000E+0_dbl,&
                          b2 = 0.4504130959426697E+1_dbl,&
                          b3 = 0.1110667363742916E+1_dbl,&
                          b4 = 0.2359291751427506E-1_dbl

  REAL(dbl), PARAMETER :: da0 = 0.119086804055547E+0_dbl,&
                          da1 = 0.6157402568883345E+0_dbl,&
                          da2 = 0.1574201515892867E+0_dbl,&
                          da3 = 0.3532336663397157E-2_dbl,&
                          db1 = 0.0000000000000000E+0_dbl,&
                          db2 = 0.2673612973836267E+0_dbl,&
                          db3 = 0.2052004607777787E+0_dbl,&
                          db4 = 0.4200005045691381E-2_dbl

  PUBLIC :: pade

  REAL(dbl) :: eps_rho
  LOGICAL :: debug_flag

! *****************************************************************************

  INTERFACE pade
    MODULE PROCEDURE pade_lda_pw, pade_lsd_pw, pade_init, pade_info
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE pade_init ( cutoff, debug )

    REAL(dbl), INTENT(IN) :: cutoff
    LOGICAL, OPTIONAL, INTENT(IN) :: debug

    eps_rho = cutoff
    CALL set_util ( cutoff )

    IF ( PRESENT ( debug ) ) THEN
       debug_flag = debug
    ELSE
       debug_flag = .FALSE.
    END IF

  END SUBROUTINE pade_init

! *****************************************************************************

  SUBROUTINE pade_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
       reference = "S. Goedecker, M. Teter and J. Hutter," &
                   //" Phys. Rev. B 54, 1703 (1996)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "S. Goedecker et al., PRB 54, 1703 (1996)"
    END IF

  END SUBROUTINE pade_info

! *****************************************************************************

  SUBROUTINE pade_lda_pw(derivative_set, rho, pot, order)

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:,:), INTENT(INOUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: rs
    INTEGER :: n, error, m, mp, mt, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,100)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "pade_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( rs(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "pade_lda", "rs", n )
    
    CALL calc_rs ( rho, rs )
    IF ( calc(0) .AND. calc(1) ) THEN
       m1=mx(0,1)
       m2=mx(1,2)
       CALL pade_lda_01(n, rho, rs, pot(:,m1:m1), pot(:,m2:m2))
    ELSE IF ( calc(0) ) THEN
       CALL pade_lda_0(n, rho, rs, pot(:,1:1))
    ELSE IF ( calc(1) ) THEN
       m1=mx(1,1)
       CALL pade_lda_1(n, rho, rs, pot(:,m1:m1))
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       CALL pade_lda_2(n, rho, rs, pot(:,m1:m1))
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       CALL pade_lda_3(n, rho, rs, pot(:,m1:m1))
    END IF

    DEALLOCATE ( rs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "pade_lda", "rs" )

  END SUBROUTINE pade_lda_pw

! *****************************************************************************

  SUBROUTINE pade_lsd_pw ( derivative_set, rhoa, rhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(INOUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: rs
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: fx
    INTEGER :: n, error, m, mp, mt, mx(0:3,2), m1, m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,201)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "pade_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( rs(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "pade_lsd", "rs", n )
    m = ABS ( order ) + 1
    ALLOCATE ( fx(n,m), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "pade_lsd", "fx", n*m )

!TC    fx(1:n,1) = rhoa(1:n) + rhob(1:n)
    CALL dcopy(n,rhoa,1,fx,1)
    CALL daxpy(n,1.0_dbl,rhob,1,fx,1)

    CALL calc_rs ( fx(:,1), rs )
    CALL calc_fx ( n, rhoa, rhob, fx, m-1 )

    IF ( calc(0) .AND. calc(1) ) THEN
       m1=mx(0,1)
       m2=mx(1,2)
       CALL pade_lsd_01 ( n, rhoa, rhob, rs, fx, &
                          pot(:,m1:m1), pot(:,m1+1:m1+1), pot(:,m2:m2) )
    ELSE IF ( calc(0) ) THEN
       CALL pade_lsd_0 ( n, rhoa, rhob, rs, fx, pot(:,1) )
    ELSE IF ( calc(1) ) THEN
       m1=mx(1,1)
       m2=mx(1,2)
       CALL pade_lsd_1 ( n, rhoa, rhob, rs, fx, &
                         pot(:,m1:m1), pot(:,m2:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       m2=mx(2,2)
       CALL pade_lsd_2 ( n, rhoa, rhob, rs, fx, &
                         pot(:,m1:m1), pot(:,m1+1:m1+1), pot(:,m2:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       m2=mx(3,2)
       CALL pade_lsd_3 ( n, rhoa, rhob, rs, fx, &
                         pot(:,m1:m1), pot(:,m1+1:m1+1), pot(:,m1+2:m1+2), pot(:,m2:m2) )
    END IF

    DEALLOCATE ( fx, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "pade_lsd", "fx" )
    DEALLOCATE ( rs, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "pade_lsd", "rs" )

  END SUBROUTINE pade_lsd_pw

! *****************************************************************************

  SUBROUTINE pade_lda_0(n, rho, rs, pot)

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rho, rs
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: pot

    INTEGER :: ip
    REAL(dbl) :: p, q, epade
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,p,q,epade)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         epade = -p/q
         pot(ip) = pot(ip) + epade*rho(ip)
      END IF
    END DO

  END SUBROUTINE pade_lda_0

! *****************************************************************************

  SUBROUTINE pade_lda_1(n, rho, rs, pot)

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rho, rs
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: pot

    INTEGER :: ip
    REAL(dbl) :: p, q, dp, dq, epade, depade
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,p,q,epade,dp,dq,depade)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         epade = -p/q

         dp = a1 + (2.0_dbl*a2 + 3.0_dbl*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dbl*b2 + (3.0_dbl*b3 + 4.0_dbl*b4*rs(ip))*rs(ip))*rs(ip)
         depade = f13 * rs(ip) * (dp*q - p*dq) / (q*q)

         pot(ip) = pot(ip) + epade + depade

      END IF
    END DO

  END SUBROUTINE pade_lda_1

! *****************************************************************************

  SUBROUTINE pade_lda_01(n, rho, rs, pot0, pot1)

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rho, rs
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: pot0, pot1

    INTEGER :: ip
    REAL(dbl) :: p, q, dp, dq, epade, depade
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,p,q,epade,dp,dq,depade)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         epade = -p/q

         dp = a1 + (2.0_dbl*a2 + 3.0_dbl*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dbl*b2 + (3.0_dbl*b3 + 4.0_dbl*b4*rs(ip))*rs(ip))*rs(ip)
         depade = f13 * rs(ip) * (dp*q - p*dq) / (q*q)

         pot0(ip) = pot0(ip) + epade * rho(ip)
         pot1(ip) = pot1(ip) + epade + depade

      END IF
    END DO

  END SUBROUTINE pade_lda_01

! *****************************************************************************

  SUBROUTINE pade_lda_2(n, rho, rs, pot)

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rho, rs
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: pot

    INTEGER :: ip
    REAL(dbl) :: p, q, dp, dq, d2p, d2q, rsr, t1, t2, t3
 !   ---------------------------------------------------------------------------

!$omp parallel do private(ip,p,q,dp,dq,d2p,d2q,rsr,t1,t2,t3)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dp = a1 + (2.0_dbl*a2 + 3.0_dbl*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dbl*b2 + (3.0_dbl*b3 + 4.0_dbl*b4*rs(ip))*rs(ip))*rs(ip)

         d2p = 2.0_dbl*a2 + 6.0_dbl*a3*rs(ip)
         d2q = 2.0_dbl*b2 + (6.0_dbl*b3 + 12.0_dbl*b4*rs(ip))*rs(ip)

         rsr = rs(ip)/rho(ip)
         t1 = (p*dq - dp*q)/(q*q)
         t2 = (d2p*q - p*d2q)/(q*q)
         t3 = (p*dq*dq - dp*q*dq)/(q*q*q)

         pot(ip) = pot(ip) - f13 * ( f23*t1 + f13*t2*rs(ip) + f23*t3*rs(ip) ) * rsr

      END IF
    END DO

  END SUBROUTINE pade_lda_2

! *****************************************************************************

  SUBROUTINE pade_lda_3(n, rho, rs, pot)

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rho, rs
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: pot

    INTEGER :: ip
    REAL(dbl) :: ab1, ab2, ab3, rsr1, rsr2, rsr3
    REAL(dbl) :: p, q, dp, dq, d2p, d2q, d3p, d3q
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,p,q,dp,dq,d2p,d2q,d3p,d3q,ab1,ab2,ab3,rsr1,rsr2,rsr3)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dp = a1 + (2.0_dbl*a2 + 3.0_dbl*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dbl*b2 + (3.0_dbl*b3 + 4.0_dbl*b4*rs(ip))*rs(ip))*rs(ip)

         d2p = 2.0_dbl*a2 + 6.0_dbl*a3*rs(ip)
         d2q = 2.0_dbl*b2 + (6.0_dbl*b3 + 12.0_dbl*b4*rs(ip))*rs(ip)

         d3p = 6.0_dbl*a3
         d3q = 6.0_dbl*b3 + 24.0_dbl*b4*rs(ip)

         ab1 = (dp*q - p*dq)/(q*q)
         ab2 = (d2p*q*q - p*q*d2q - 2._dbl*dp*q*dq + 2._dbl*p*dq*dq)/(q*q*q)
         ab3 = (d3p*q*q - p*q*d3q - 3._dbl*dp*q*d2q + 3._dbl*p*dq*d2q)/(q*q*q)
         ab3 = ab3 - 3._dbl*ab2*dq/q
         rsr1 = rs(ip) / (rho(ip)*rho(ip))
         rsr2 = f13*f13*rs(ip) * rsr1
         rsr3 = f13*rs(ip) * rsr2
         rsr1 = -f23*f23*f23*rsr1
         pot(ip) = pot(ip) + rsr1*ab1 + rsr2*ab2 + rsr3*ab3

      END IF
    END DO

  END SUBROUTINE pade_lda_3

! *****************************************************************************

  SUBROUTINE pade_lsd_0(n, rhoa, rhob, rs, fx, pot0 )

    INTEGER, INTENT(IN)                  :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)  :: rhoa, rhob, rs
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(*), INTENT(INOUT)  :: pot0

    INTEGER :: ip
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, rhoab
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,rhoab,p,q) &
!$omp private(fa0,fa1,fa2,fa3,fb1,fb2,fb3,fb4)
    DO ip = 1, n

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(ip,1)*da0
         fa1 = a1 + fx(ip,1)*da1
         fa2 = a2 + fx(ip,1)*da2
         fa3 = a3 + fx(ip,1)*da3
         fb1 = b1 + fx(ip,1)*db1
         fb2 = b2 + fx(ip,1)*db2
         fb3 = b3 + fx(ip,1)*db3
         fb4 = b4 + fx(ip,1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs(ip))*rs(ip))*rs(ip)
         q = (fb1 + (fb2 + (fb3 + fb4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         pot0(ip) = pot0(ip) - p/q*rhoab

      END IF

    END DO

  END SUBROUTINE pade_lsd_0

! *****************************************************************************

  SUBROUTINE pade_lsd_1 ( n, rhoa, rhob, rs, fx, pota, potb )

    REAL(dbl), DIMENSION(*), INTENT(IN)  :: rhoa, rhob, rs
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(*), INTENT(INOUT)  :: pota, potb

    INTEGER :: ip, n
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,rhoab,p,q,dp,dq,xp,xq,dr,dx,dc) &
!$omp private(fa0,fa1,fa2,fa3,fb1,fb2,fb3,fb4)
    DO ip = 1, n

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(ip,1)*da0
         fa1 = a1 + fx(ip,1)*da1
         fa2 = a2 + fx(ip,1)*da2
         fa3 = a3 + fx(ip,1)*da3
         fb1 = b1 + fx(ip,1)*db1
         fb2 = b2 + fx(ip,1)*db2
         fb3 = b3 + fx(ip,1)*db3
         fb4 = b4 + fx(ip,1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs(ip))*rs(ip))*rs(ip)
         q = (fb1 + (fb2 + (fb3 + fb4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         dp = fa1 + (2.0_dbl*fa2 + 3.0_dbl*fa3*rs(ip))*rs(ip)
         dq = fb1 + (2.0_dbl*fb2 + (3.0_dbl*fb3 + &
                     4.0_dbl*fb4*rs(ip))*rs(ip))*rs(ip)
         xp = da0 + (da1 + (da2 + da3*rs(ip))*rs(ip))*rs(ip)
         xq = (db1 + (db2 + (db3 + db4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dr = (dp*q - p*dq)/(q*q)
         dx = 2._dbl*(xp*q - p*xq)/(q*q)*fx(ip,2)/rhoab
         dc = f13*rs(ip)*dr - p/q

         pota(ip) = pota(ip) + dc - dx*rhob(ip)
         potb(ip) = potb(ip) + dc + dx*rhoa(ip)

      END IF

    END DO

  END SUBROUTINE pade_lsd_1

! *****************************************************************************

  SUBROUTINE pade_lsd_01 ( n, rhoa, rhob, rs, fx, pot0, pota, potb )

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rhoa, rhob, rs
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: pot0, pota, potb

    INTEGER :: ip
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,rhoab,p,q,dp,dq,xp,xq,dr,dx,dc) &
!$omp private(fa0,fa1,fa2,fa3,fb1,fb2,fb3,fb4)
    DO ip = 1, n

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(ip,1)*da0
         fa1 = a1 + fx(ip,1)*da1
         fa2 = a2 + fx(ip,1)*da2
         fa3 = a3 + fx(ip,1)*da3
         fb1 = b1 + fx(ip,1)*db1
         fb2 = b2 + fx(ip,1)*db2
         fb3 = b3 + fx(ip,1)*db3
         fb4 = b4 + fx(ip,1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs(ip))*rs(ip))*rs(ip)
         q = (fb1 + (fb2 + (fb3 + fb4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         dp = fa1 + (2.0_dbl*fa2 + 3.0_dbl*fa3*rs(ip))*rs(ip)
         dq = fb1 + (2.0_dbl*fb2 + (3.0_dbl*fb3 + &
                     4.0_dbl*fb4*rs(ip))*rs(ip))*rs(ip)
         xp = da0 + (da1 + (da2 + da3*rs(ip))*rs(ip))*rs(ip)
         xq = (db1 + (db2 + (db3 + db4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dr = (dp*q - p*dq)/(q*q)
         dx = 2._dbl*(xp*q - p*xq)/(q*q)*fx(ip,2)/rhoab
         dc = f13*rs(ip)*dr - p/q

         pot0(ip) = pot0(ip) - p/q*rhoab
         pota(ip) = pota(ip) + dc - dx*rhob(ip)
         potb(ip) = potb(ip) + dc + dx*rhoa(ip)

      END IF

    END DO

  END SUBROUTINE pade_lsd_01

! *****************************************************************************

  SUBROUTINE pade_lsd_2 (n, rhoa, rhob, rs, fx, potaa, potab, potbb )

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rhoa, rhob, rs
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: potaa, potab, potbb

    INTEGER :: ip
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
    REAL(dbl) :: d2p, d2q, dxp, dxq, drr, dxr, dxx, xt, yt, or
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,rhoab,p,q,dp,dq,d2p,d2q,xp,xq,dxp,dxq,dr,drr) &
!$omp private(dx,dxx,dxr,or,xt,yt) &
!$omp private(fa0,fa1,fa2,fa3,fb1,fb2,fb3,fb4)
    DO ip = 1, n

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(ip,1)*da0
         fa1 = a1 + fx(ip,1)*da1
         fa2 = a2 + fx(ip,1)*da2
         fa3 = a3 + fx(ip,1)*da3
         fb1 = b1 + fx(ip,1)*db1
         fb2 = b2 + fx(ip,1)*db2
         fb3 = b3 + fx(ip,1)*db3
         fb4 = b4 + fx(ip,1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs(ip))*rs(ip))*rs(ip)
         q = (fb1 + (fb2 + (fb3 + fb4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dp = fa1 + (2.0_dbl*fa2 + 3.0_dbl*fa3*rs(ip))*rs(ip)
         dq = fb1 + (2.0_dbl*fb2 + (3.0_dbl*fb3 + &
                     4.0_dbl*fb4*rs(ip))*rs(ip))*rs(ip)

         d2p = 2.0_dbl*fa2 + 6.0_dbl*fa3*rs(ip)
         d2q = 2.0_dbl*fb2 + (6.0_dbl*fb3 + 12.0_dbl*fb4*rs(ip))*rs(ip)

         xp = da0 + (da1 + (da2 + da3*rs(ip))*rs(ip))*rs(ip)
         xq = (db1 + (db2 + (db3 + db4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dxp = da1 + (2.0_dbl*da2 + 3.0_dbl*da3*rs(ip))*rs(ip)
         dxq = db1 + (2.0_dbl*db2 + (3.0_dbl*db3 + &
                     4.0_dbl*db4*rs(ip))*rs(ip))*rs(ip)

         dr = (dp*q - p*dq)/(q*q)
         drr = (d2p*q*q-p*q*d2q-2._dbl*dp*q*dq+2._dbl*p*dq*dq)/(q*q*q)
         dx = (xp*q - p*xq)/(q*q)
         dxx = 2._dbl*xq*(p*xq - xp*q)/(q*q*q)
         dxr = (dxp*q*q+dp*xq*q-xp*dq*q-p*dxq*q-2._dbl*dp*q*xq+2._dbl*p*dq*xq)/(q*q*q)

         or = 1._dbl/rhoab
         yt = rhob(ip)*or
         xt = rhoa(ip)*or

         potaa(ip) = potaa(ip) + f23*f13*dr*rs(ip)*or - f13*f13*drr*rs(ip)*rs(ip)*or &
                     + f43*rs(ip)*fx(ip,2)*dxr*yt*or &
                     - 4._dbl*fx(ip,2)*fx(ip,2)*dxx*yt*yt*or &
                     - 4._dbl*dx*fx(ip,3)*yt*yt*or
         potab(ip) = potab(ip) + f23*f13*dr*rs(ip)*or - f13*f13*drr*rs(ip)*rs(ip)*or &
                     + f23*rs(ip)*fx(ip,2)*dxr*(yt-xt)*or &
                     + 4._dbl*fx(ip,2)*fx(ip,2)*dxx*xt*yt*or &
                     + 4._dbl*dx*fx(ip,3)*xt*yt*or
         potbb(ip) = potbb(ip) + f23*f13*dr*rs(ip)*or - f13*f13*drr*rs(ip)*rs(ip)*or &
                     - f43*rs(ip)*fx(ip,2)*dxr*xt*or &
                     - 4._dbl*fx(ip,2)*fx(ip,2)*dxx*xt*xt*or &
                     - 4._dbl*dx*fx(ip,3)*xt*xt*or

      END IF

    END DO

  END SUBROUTINE pade_lsd_2

! *****************************************************************************

  SUBROUTINE pade_lsd_3 ( n, rhoa, rhob, rs, fx, potaaa, potaab, potabb, potbbb )

    INTEGER, INTENT(IN)                    :: n
    REAL(dbl), DIMENSION(*), INTENT(IN)    :: rhoa, rhob, rs
    REAL(dbl), DIMENSION(:,:), INTENT(IN)  :: fx
    REAL(dbl), DIMENSION(*), INTENT(INOUT) :: potaaa, potaab, potabb, potbbb

    INTEGER :: ip
    REAL(dbl) :: fa0, fa1, fa2, fa3, fb1, fb2, fb3, fb4
    REAL(dbl) :: p, q, dp, dq, xp, xq, rhoab, dr, dx, dc
    REAL(dbl) :: d2p, d2q, dxp, dxq, drr, dxr, dxx, xt, yt, or
    REAL(dbl) :: d3p, d3q, d2xp, d2xq
    REAL(dbl) :: drrr, dxrr, dxxr, dxxx
!   ---------------------------------------------------------------------------

    IF ( .NOT. debug_flag ) CALL stop_program("pade_lsd_3","Routine not tested")

!$omp parallel do private(ip,fa0,fa1,fa2,fa3,fb1,fb2,fb3,fb4) &
!$omp  private(p,q,dp,dq,xp,xq,rhoab,dr,dx,dc,d2p,d2q,dxp,dxq,drr,dxr,dxx) &
!$omp  private(xt,yt,or,d3p,d3q,d2xp,d2xq,drrr,dxrr,dxxr,dxxx)
    DO ip = 1, n

      rhoab = rhoa(ip) + rhob(ip)

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(ip,1)*da0
         fa1 = a1 + fx(ip,1)*da1
         fa2 = a2 + fx(ip,1)*da2
         fa3 = a3 + fx(ip,1)*da3
         fb1 = b1 + fx(ip,1)*db1
         fb2 = b2 + fx(ip,1)*db2
         fb3 = b3 + fx(ip,1)*db3
         fb4 = b4 + fx(ip,1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs(ip))*rs(ip))*rs(ip)
         q = (fb1 + (fb2 + (fb3 + fb4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dp = fa1 + (2.0_dbl*fa2 + 3.0_dbl*fa3*rs(ip))*rs(ip)
         dq = fb1 + (2.0_dbl*fb2 + (3.0_dbl*fb3 + &
                     4.0_dbl*fb4*rs(ip))*rs(ip))*rs(ip)

         d2p = 2.0_dbl*fa2 + 6.0_dbl*fa3*rs(ip)
         d2q = 2.0_dbl*fb2 + (6.0_dbl*fb3 + 12.0_dbl*fb4*rs(ip))*rs(ip)

         d3p = 6.0_dbl*fa3
         d3q = 6.0_dbl*fb3 + 24.0_dbl*fb4*rs(ip)

         xp = da0 + (da1 + (da2 + da3*rs(ip))*rs(ip))*rs(ip)
         xq = (db1 + (db2 + (db3 + db4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dxp = da1 + (2.0_dbl*da2 + 3.0_dbl*da3*rs(ip))*rs(ip)
         dxq = db1 + (2.0_dbl*db2 + (3.0_dbl*db3 + &
                     4.0_dbl*db4*rs(ip))*rs(ip))*rs(ip)

         d2xp = 2.0_dbl*da2 + 6.0_dbl*da3*rs(ip)
         d2xq = 2.0_dbl*db2 + (6.0_dbl*db3 + 12.0_dbl*db4*rs(ip))*rs(ip)

         dr = (dp*q - p*dq)/(q*q)
         drr = (d2p*q*q-p*q*d2q-2._dbl*dp*q*dq+2._dbl*p*dq*dq)/(q*q*q)
         drrr = (d3p*q*q*q-3._dbl*d2p*dq*q*q+6._dbl*dp*dq*dq*q-3._dbl*dp*d2q*q*q - &
                6._dbl*p*dq*dq*dq+6._dbl*p*dq*d2q*q-p*d3q*q*q)/(q*q*q*q)
         dx = (xp*q - p*xq)/(q*q)
         dxx = 2._dbl*xq*(p*xq - xp*q)/(q*q*q)
         dxxx = 6._dbl*xq*(q*xp*xq-p*xq*xq)/(q*q*q*q)
         dxr = (dxp*q*q+dp*xq*q-xp*dq*q-p*dxq*q-2._dbl*dp*q*xq+2._dbl*p*dq*xq)/(q*q*q)
         dxxr = 2._dbl*(2._dbl*dxq*q*p*xq-dxq*q*q*xp+xq*xq*q*dp-xq*q*q*dxp + &
                2._dbl*xq*q*xp*dq-3._dbl*xq*xq*dq*p)/(q*q*q*q)
         dxrr = (q*q*q*d2xp-2._dbl*q*q*dxp*dq-q*q*xp*d2q-q*q*d2p*xq - &
                2._dbl*q*q*dp*dxq-q*q*p*d2xq+4._dbl*dq*q*dp*xq+4._dbl*dq*q*p*dxq + &
                2._dbl*dq*dq*q*xp-6._dbl*dq*dq*p*xq+2._dbl*d2q*q*p*xq)/(q*q*q*q)

         or = 1._dbl/rhoab
         yt = rhob(ip)*or
         xt = rhoa(ip)*or

         potaaa(ip) = potaaa(ip) + 8._dbl/27._dbl*dr*rs(ip)*or*or + &
                     1._dbl/9._dbl*drr*rs(ip)*rs(ip)*or*or + &
                     1._dbl/27._dbl*drrr*rs(ip)**3*or*or + &
                     dxr*or*or*yt*rs(ip)*(-8._dbl/3._dbl*fx(ip,2)+4._dbl*fx(ip,3)*yt) 
         potaab(ip) = potaab(ip) + 0._dbl
         potabb(ip) = potabb(ip) + 0._dbl
         potbbb(ip) = potbbb(ip) + 0._dbl

      END IF

    END DO

  END SUBROUTINE pade_lsd_3

! *****************************************************************************

END MODULE xc_pade

! *****************************************************************************
