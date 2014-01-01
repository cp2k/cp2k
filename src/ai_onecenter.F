!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!   Calculates atomic integrals over unnormalized spherical Gaussian functions
!-----------------------------------------------------------------------------!
!
!   phi(r) = r^l * exp[-p*r^2] Ylm
!
!-----------------------------------------------------------------------------!
!   Calculates atomic integrals over normalized Slater type functions
!-----------------------------------------------------------------------------!
!
!   phi(r) = N(nlm) r^(n-1) * exp[-p*r] Ylm
!   N(nlm) = [(2n)!]^(-1/2) (2p)^(n+1/2)
!
!-----------------------------------------------------------------------------!
!   Calculates atomic integrals over spherical numerical functions
!-----------------------------------------------------------------------------!
!
!   phi(r) = R(r) Ylm
!
!-----------------------------------------------------------------------------!
MODULE ai_onecenter

  USE erf_fn,                          ONLY: erf
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             gamma0,&
                                             gamma1,&
                                             pi
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sg_overlap, sg_kinetic, sg_nuclear, sg_erf, sg_gpot, &
            sg_proj_ol, sg_conf, sg_coulomb, sg_exchange, sg_kinnuc, &
            sg_erfc
  PUBLIC :: sto_overlap, sto_kinetic, sto_nuclear, sto_conf
  PUBLIC :: num_overlap, num_kinetic, num_nuclear, num_erf, num_gpot, &
            num_proj_ol, num_conf, num_kinnuc

CONTAINS

!------------------------------------------------------------------------------
!
!  S(l,pq) = pi^(1/2)*(2*l+1)!! / 2^(l+2) / (p+q)^(l+1.5)
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_overlap ( smat, l, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: smat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: el, spi

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( smat, 1 ) .OR. m > SIZE ( smat, 2 ) ) STOP

    spi = SQRT(pi)/2.0_dp**(l+2)*dfac(2*l+1)
    el  = REAL(l,dp) + 1.5_dp

    DO iq = 1, m
      DO ip = 1, n
        smat ( ip, iq ) = spi/(pa(ip)+pb(iq))**el
      END DO
    END DO

  END SUBROUTINE sg_overlap

!------------------------------------------------------------------------------
!
!  T(l,pq) = (2l+3)!! pi^(1/2)/2^(l+2) [pq/(p+q)^(l+2.5)]
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_kinetic ( kmat, l, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: kmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: spi

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( kmat, 1 ) .OR. m > SIZE ( kmat, 2 ) ) STOP

    spi = dfac(2*l+3)*SQRT(pi)/2.0_dp**(l+2)
    DO iq = 1, m
      DO ip = 1, n
        kmat ( ip, iq ) = spi*pa(ip)*pb(iq)/(pa(ip)+pb(iq))**(l+2.5_dp)
      END DO
    END DO

  END SUBROUTINE sg_kinetic

!------------------------------------------------------------------------------
!
!  U(l,pq) = l!/2 / (p+q)^(l+1)
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_nuclear ( umat, l, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: umat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: tld

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) STOP

    tld = 0.5_dp*fac(l)
    DO iq = 1, m
      DO ip = 1, n
        umat ( ip, iq ) = tld / (pa(ip)+pb(iq))**(l+1)
      END DO
    END DO

  END SUBROUTINE sg_nuclear

!------------------------------------------------------------------------------
!
!  U(l,pq) = l!/2 / (p+q)^l * [4/(p+q)^2 *pq*(l+1) + 1]
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_kinnuc ( umat, l, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: umat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: ppq, pq, tld

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) STOP

    IF ( l > 0 ) THEN
      tld = 0.5_dp*fac(l)
      DO iq = 1, m
        DO ip = 1, n
          ppq = pa(ip)+pb(iq)
          pq = pa(ip)*pb(iq)
          umat ( ip, iq ) = tld / ppq**l * (4.0_dp/ppq**2*pq*REAL(l+1,dp)+1.0_dp)
        END DO
      END DO
    ELSE
      DO iq = 1, m
        DO ip = 1, n
          ppq = pa(ip)+pb(iq)
          pq = pa(ip)*pb(iq)
          umat ( ip, iq ) = 2.0_dp*pq/ppq**2
        END DO
      END DO
    END IF

  END SUBROUTINE sg_kinnuc

!------------------------------------------------------------------------------
!
!  z = a/(p+q)
!
!  UP(l,pq,a) = Gamma(l+3/2)*a/SQRT(Pi)/(p+q)^(l+3/2)*
!                      Hypergeom([1/2, 3/2 + l], [3/2], -z)
!
!  UP(l,pq,a) = a/2^(l+1)/(p+q)^(l+3/2)/(1+z)^(l+1/2) * F(z,l)
!
!  F(z,0) = 1
!  F(z,1) = 3 + 2*z
!  F(z,2) = 15 + 20*z + 8*z^2
!  F(z,3) = 35 + 70*z + 56*z^2 + 16*z^3
!  F(z,4) = 315 + 840*z + 1008*z^2 + 576*z^3 + 128*z^4
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_erf ( upmat, l, a, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: upmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: handle, ip, iq, m, n
    REAL(KIND=dp)                            :: a2, fpol, pq, tld, z

    CALL timeset("sg_erf",handle)

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( upmat, 1 ) .OR. m > SIZE ( upmat, 2 ) ) STOP

    a2 = a*a
    tld = a / 2._dp**(l+1)
    DO iq = 1, m
      DO ip = 1, n
        pq = pa(ip) + pb(iq)
        z = a2/pq
        upmat ( ip, iq ) = tld / (1._dp+z)**(l+0.5_dp) / pq**(l+1.5_dp)
      END DO
    END DO

    DO iq = 1, m
      SELECT CASE ( l )
      CASE DEFAULT
        STOP
      CASE ( 0 )
        ! nothing left to do
      CASE ( 1 )
        DO ip = 1, n
          pq = pa(ip) + pb(iq)
          z = a2/pq
          fpol = 2.0_dp*z + 3.0_dp
          upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
        END DO
      CASE ( 2 )
        DO ip = 1, n
          pq = pa(ip) + pb(iq)
          z = a2/pq
          fpol = 8.0_dp*z*z + 20.0_dp*z + 15.0_dp
          upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
        END DO
      CASE ( 3 )
        DO ip = 1, n
          pq = pa(ip) + pb(iq)
          z = a2/pq
          fpol = 16.0_dp*z*z*z + 56.0_dp*z*z + 70.0_dp*z + 35.0_dp
          fpol = 3._dp*fpol
          upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
        END DO
      CASE ( 4 )
        DO ip = 1, n
          pq = pa(ip) + pb(iq)
          z = a2/pq
          fpol = 128.0_dp*z*z*z*z + 576.0_dp*z*z*z + 1008.0_dp*z*z + 840.0_dp*z + 315.0_dp
          fpol = 3._dp*fpol
          upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
        END DO
      END SELECT
    END DO

    CALL timestop(handle)

  END SUBROUTINE sg_erf

!------------------------------------------------------------------------------
!
!  Overlap with Projectors P(l,k,rc) for k=0,1,..
!
!  P(l,k,rc) = SQRT(2)/SQRT(Gamma[l+2k+1.5])/rc^(l+2k+1.5) r^(l+2k) exp[-0.5(r/rc)^2]
!
!  SP(l,k,p,rc) = 2^(l+k+1) / SQRT(gamma[l+2k+1.5]) / rc^(l+2k+1.5)
!                    * Gamma(l+k+1.5) / (2p+1/rc^2)^(l+k+1.5)
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_proj_ol ( spmat, l, p, k, rc )


    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: spmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: p
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                :: rc

    REAL(KIND=dp)                            :: orc, pf

    IF ( SIZE ( spmat ) < SIZE ( p ) ) STOP

    pf = 2._dp**(l+k+1) * gamma1(l+k+1) / rc**(l+2*k+1.5_dp) / SQRT(gamma1(l+2*k+1))
    orc = 1._dp/(rc*rc)
    spmat(:) = pf / ( 2._dp*p(:) + orc )**(l+k+1.5_dp)

  END SUBROUTINE sg_proj_ol

!------------------------------------------------------------------------------
!
!  Matrix elements for Gaussian potentials
!
!  V(k,rc) = (r/rc)^2k exp[-1/2(r/rc)^2]
!
!  VP(l,k,p+q,rc) = 2^(l+k+0.5) * rc^(2l+3) * Gamma(l+k+1.5) / (1+2rc^2(p+q))^(l+k+1.5)
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_gpot ( vpmat, k, rc, l, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: vpmat
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                :: rc
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: tld

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( vpmat, 1 ) .OR. m > SIZE ( vpmat, 2 ) ) STOP

    tld = gamma1(l+k+1) * rc**(2*l+3) * 2._dp**(l+k+0.5)

    DO iq = 1, m
      DO ip = 1, n
        vpmat ( ip, iq ) = tld / ( 1._dp + 2._dp*rc*rc*(pa(ip)+pb(iq)) )**(l+k+1.5_dp)
      END DO
    END DO

  END SUBROUTINE sg_gpot

!------------------------------------------------------------------------------
!
!  G(l,k,pq) = <a|[r/rc]^2k|b>
!            = 0.5*Gamma(l+k+1.5)/rc^(2k)/(p+q)^(l+k+1.5)
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_conf ( gmat, rc, k, l, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: gmat
    REAL(KIND=dp), INTENT(IN)                :: rc
    INTEGER, INTENT(IN)                      :: k, l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: tld

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( gmat, 1 ) .OR. m > SIZE ( gmat, 2 ) ) STOP

    tld = 0.5_dp/rc**(2*k)*gamma1(l+k+1)
    DO iq = 1, m
      DO ip = 1, n
        gmat ( ip, iq ) = tld/(pa(ip)+pb(iq))**(l+k+1.5_dp)
      END DO
    END DO

  END SUBROUTINE sg_conf

!------------------------------------------------------------------------------
!
!  (plql,rl'sl')
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_coulomb ( eri, nu, pa, lab, pc, lcd )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: eri
    INTEGER, INTENT(IN)                      :: nu
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa
    INTEGER, INTENT(IN)                      :: lab
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pc
    INTEGER, INTENT(IN)                      :: lcd

    INTEGER                                  :: handle, ia, ib, ic, id, jab, &
                                                jcd, na, nc
    REAL(KIND=dp)                            :: cc1, cc2, p, q, r, s, sab, &
                                                scd, slab, slcd, ss, vab1, &
                                                vab2, vab3, vcd1, vcd2, vcd3, &
                                                xab, xcd

    CALL timeset("sg_coulomb",handle)

    na = SIZE(pa)
    nc = SIZE(pc)
    ss   = 2.0_dp*SQRT(2.0_dp/pi)
    slab = SQRT(pi)*dfac(2*lab+1)/2.0_dp**(lab+2)
    slcd = SQRT(pi)*dfac(2*lcd+1)/2.0_dp**(lcd+2)
    jab = 0
    DO ia=1,na
      p = pa(ia)
      DO ib=ia,na
        jab = jab + 1
        q = pa(ib)
        xab = 0.5_dp*(p+q)
        sab = slab / (p+q)**(lab+1.5_dp)
        vab1= vgau(2*lab-nu+1,xab)
        vab2= vgau(2*lab+2,xab)
        vab3= vgau(2*lab+nu+2,xab)
        jcd = 0
        DO ic=1,nc
          r = pc(ic)
          DO id=ic,nc
            jcd = jcd + 1
            s = pc(id)
            xcd = 0.5_dp*(r+s)
            scd = slcd / (r+s)**(lcd+1.5_dp)
            vcd1= vgau(2*lcd+nu+2,xcd)
            vcd2= vgau(2*lcd+2,xcd)
            vcd3= vgau(2*lcd-nu+1,xcd)
            cc1 = cgau(2*lab-nu+1,2*lcd+nu+2,xab/xcd)
            cc2 = cgau(2*lcd-nu+1,2*lab+nu+2,xcd/xab)

            eri(jab,jcd) = ss * sab * scd / (vab2*vcd2) * (cc1*vab1*vcd1 + cc2*vab3*vcd3)

          END DO
        END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE sg_coulomb

!------------------------------------------------------------------------------
!
!  (plql',rlsl')
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_exchange ( eri, nu, pa, lac, pb, lbd )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: eri
    INTEGER, INTENT(IN)                      :: nu
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa
    INTEGER, INTENT(IN)                      :: lac
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pb
    INTEGER, INTENT(IN)                      :: lbd

    INTEGER                                  :: handle, ia, ib, ic, id, jac, &
                                                jbd, na, nb
    REAL(KIND=dp) :: cc1, cc2, cc3, cc4, ee, p, q, r, s, sac, sbd, slac, &
      slbd, ss, v1pr, v1ps, v1qr, v1qs, v2pr, v2ps, v2qr, v2qs, vpq, vrs, &
      xab, xac, xad, xbc, xbd, xcd

    CALL timeset("sg_exchange",handle)

    eri=0.0_dp
    na = SIZE(pa)
    nb = SIZE(pb)
    ss   = 2.0_dp*SQRT(2.0_dp/pi)
    slac = SQRT(pi)*dfac(2*lac+1)/2.0_dp**(lac+2)
    slbd = SQRT(pi)*dfac(2*lbd+1)/2.0_dp**(lbd+2)
    jac = 0
    DO ia=1,na
      p = pa(ia)
      DO ic=ia,na
        jac = jac + 1
        q = pa(ic)
        xac = 0.5_dp*(p+q)
        sac = slac / (p+q)**(lac+1.5_dp)
        vpq = vgau(2*lac+2,xac)
        jbd = 0
        DO ib=1,nb
          r = pb(ib)
          xab = 0.5_dp*(p+r)
          xbc = 0.5_dp*(q+r)
          DO id=ib,nb
            jbd = jbd + 1
            s = pb(id)
            xbd = 0.5_dp*(r+s)
            xcd = 0.5_dp*(q+s)
            xad = 0.5_dp*(p+s)
            sbd = slbd / (r+s)**(lbd+1.5_dp)
            vrs = vgau(2*lbd+2,xbd)
            v1pr= vgau(lac+lbd-nu+1,xab)
            v1qs= vgau(lac+lbd-nu+1,xcd)
            v1ps= vgau(lac+lbd-nu+1,xad)
            v1qr= vgau(lac+lbd-nu+1,xbc)
            v2qs= vgau(lac+lbd+nu+2,xcd)
            v2pr= vgau(lac+lbd+nu+2,xab)
            v2qr= vgau(lac+lbd+nu+2,xbc)
            v2ps= vgau(lac+lbd+nu+2,xad)
            cc1 = cgau(lac+lbd-nu+1,lac+lbd+nu+2,xab/xcd)
            cc2 = cgau(lac+lbd-nu+1,lac+lbd+nu+2,xcd/xab)
            cc3 = cgau(lac+lbd-nu+1,lac+lbd+nu+2,xad/xbc)
            cc4 = cgau(lac+lbd-nu+1,lac+lbd+nu+2,xbc/xad)

            ee = 0.5_dp * ss * sac * sbd / (vpq*vrs)
            eri(jac,jbd) = ee * (v1pr*v2qs*cc1 + v1qs*v2pr*cc2 + &
                                 v1ps*v2qr*cc3 + v1qr*v2ps*cc4 )

          END DO
        END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE sg_exchange

!------------------------------------------------------------------------------
!
!  erfc(a*x)/x = 1/x - erf(a*x)/x
!
!------------------------------------------------------------------------------
  SUBROUTINE sg_erfc ( umat, l, a, pa, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: umat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa, pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: tld

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) STOP

    CALL sg_erf ( umat, l, a, pa, pb )

    tld = 0.5_dp*fac(l)
    DO iq = 1, m
      DO ip = 1, n
        umat ( ip, iq ) = tld / (pa(ip)+pb(iq))**(l+1) - umat ( ip, iq )
      END DO
    END DO

  END SUBROUTINE sg_erfc

! ******************************************************************************

  FUNCTION vgau(n,x) RESULT(v)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: v

    v = dfac(n-1)/x**(0.5_dp*(n+1))

  END FUNCTION vgau

  FUNCTION cgau(a,b,t) RESULT(c)
    INTEGER, INTENT(IN)                      :: a, b
    REAL(KIND=dp), INTENT(IN)                :: t
    REAL(KIND=dp)                            :: c

    INTEGER                                  :: l

    c = 0.0_dp
    DO l=0,(a-1)/2
      c = c + (t/(1.0_dp+t))**l * dfac(2*l+b-1)/dfac(2*l)
    END DO
    c = c * (1.0_dp + t)**(-0.5_dp*(b+1))/dfac(b-1)

  END FUNCTION cgau

!------------------------------------------------------------------------------
!
!  S(l,pn,qm) = ( V[2n,p]*V[2m,q] )^(-1/2) * V[n+m,(p+q)/2]
!
!------------------------------------------------------------------------------
  SUBROUTINE sto_overlap ( smat,  na, pa, nb, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: smat
    INTEGER, DIMENSION(:), INTENT(IN)        :: na
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa
    INTEGER, DIMENSION(:), INTENT(IN)        :: nb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: vp, vpq, vq

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( smat, 1 ) .OR. m > SIZE ( smat, 2 ) ) STOP

    DO iq = 1, m
      vq = vsto(2*nb(iq),pb(iq))
      DO ip = 1, n
        vp = vsto(2*na(ip),pa(ip))
        vpq = vsto(na(ip)+nb(iq),0.5_dp*(pa(ip)+pb(iq)))
        smat(ip,iq) = vpq/SQRT(vp*vq)
      END DO
    END DO

  END SUBROUTINE sto_overlap

!------------------------------------------------------------------------------
!
!  T(l,pn,qm) = 0.5*p*q*( V[2n,p]*V[2m,q] )^(-1/2) * V[n+m,(p+q)/2]
!                -(W[l,n,p]+W[l,m,q]) * V[n+m-1,(p+q)/2]
!                +W[l,n,p]*W[l,m,q] * V[n+m-2,(p+q)/2]
!
!------------------------------------------------------------------------------
  SUBROUTINE sto_kinetic ( kmat, l, na, pa, nb, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: kmat
    INTEGER, INTENT(IN)                      :: l
    INTEGER, DIMENSION(:), INTENT(IN)        :: na
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa
    INTEGER, DIMENSION(:), INTENT(IN)        :: nb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: vp, vpq, vpq1, vpq2, vq, wp, &
                                                wq

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( kmat, 1 ) .OR. m > SIZE ( kmat, 2 ) ) STOP

    DO iq = 1, m
      vq = vsto(2*nb(iq),pb(iq))
      wq = wsto(l,nb(iq),pb(iq))
      DO ip = 1, n
        vp = vsto(2*na(ip),pa(ip))
        vpq = vsto(na(ip)+nb(iq),0.5_dp*(pa(ip)+pb(iq)))
        vpq1 = vsto(na(ip)+nb(iq)-1,0.5_dp*(pa(ip)+pb(iq)))
        vpq2 = vsto(na(ip)+nb(iq)-2,0.5_dp*(pa(ip)+pb(iq)))
        wp = wsto(l,na(ip),pa(ip))
        kmat ( ip, iq ) = 0.5_dp*pa(ip)*pb(iq)/SQRT(vp*vq)*&
            (vpq -(wp+wq)*vpq1 + wp*wq*vpq2)
      END DO
    END DO

  END SUBROUTINE sto_kinetic

!------------------------------------------------------------------------------
!
!  U(l,pq) = 2( V[2n,p]*V[2m,q] )^(-1/2) * V[n+m-1,(p+q)/2]
!
!------------------------------------------------------------------------------
  SUBROUTINE sto_nuclear ( umat, na, pa, nb, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: umat
    INTEGER, DIMENSION(:), INTENT(IN)        :: na
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa
    INTEGER, DIMENSION(:), INTENT(IN)        :: nb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pb

    INTEGER                                  :: ip, iq, m, n
    REAL(KIND=dp)                            :: vp, vpq1, vq

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) STOP

    DO iq = 1, m
      vq = vsto(2*nb(iq),pb(iq))
      DO ip = 1, n
        vp = vsto(2*na(ip),pa(ip))
        vpq1 = vsto(na(ip)+nb(iq)-1,0.5_dp*(pa(ip)+pb(iq)))
        umat ( ip, iq ) = 2._dp/SQRT(vp*vq)*vpq1
      END DO
    END DO

  END SUBROUTINE sto_nuclear

!------------------------------------------------------------------------------
!
!  G(l,k,pq) = <aln|[r/rc]^2k|blm>
!            = N(aln)*N(blm) (a+b)^(-(n+m+2k+1))/rc^2k * GAMMA(n+m+2k+1)
!
!------------------------------------------------------------------------------
  SUBROUTINE sto_conf ( gmat, rc, k, na, pa, nb, pb )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: gmat
    REAL(KIND=dp), INTENT(IN)                :: rc
    INTEGER, INTENT(IN)                      :: k
    INTEGER, DIMENSION(:), INTENT(IN)        :: na
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pa
    INTEGER, DIMENSION(:), INTENT(IN)        :: nb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pb

    INTEGER                                  :: ip, iq, m, n

    n = SIZE ( pa )
    m = SIZE ( pb )

    IF ( n > SIZE ( gmat, 1 ) .OR. m > SIZE ( gmat, 2 ) ) STOP

    DO iq = 1, m
      DO ip = 1, n
        gmat ( ip, iq ) = (2._dp*pa(ip))**(na(ip)+0.5_dp)/SQRT(fac(2*na(ip))) &
                        * (2._dp*pb(iq))**(nb(iq)+0.5_dp)/SQRT(fac(2*nb(iq))) &
                        / rc**(2*k) / (pa(ip)+pb(iq))**(na(ip)+nb(iq)+2*k+1)  &
                        * gamma0(na(ip)+nb(iq)+2*k+1)
      END DO
    END DO

  END SUBROUTINE sto_conf

! *****************************************************************************

  FUNCTION vsto(n,x) RESULT(v)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: v

    v = fac(n)/x**(n+1)

  END FUNCTION vsto

  FUNCTION wsto(n,m,x) RESULT(w)
    INTEGER, INTENT(IN)                      :: n, m
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: w

    w = 2._dp*REAL(m-n-1,dp)/x

  END FUNCTION wsto
!------------------------------------------------------------------------------
!
!  S(l,pq) = INT(u^2 Ra(u) Rb(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_overlap ( smat, ra, rb, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: smat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, rb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wr

    INTEGER                                  :: ip, iq, m, n

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( smat, 1 ) .OR. m > SIZE ( smat, 2 ) ) STOP

    DO iq = 1, m
      DO ip = 1, n
        smat ( ip, iq ) = SUM(ra(:,ip)*rb(:,iq)*wr(:))
      END DO
    END DO

  END SUBROUTINE num_overlap

!------------------------------------------------------------------------------
!
!  T(l,pq) = 0.5 INT( u^2 dRa(u) dRb(u) + l(l+1) Ra(u) Rb(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_kinetic ( kmat, l, ra, dra, rb, drb, r, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: kmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, dra, rb, drb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ip, iq, m, n

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( kmat, 1 ) .OR. m > SIZE ( kmat, 2 ) ) STOP

    DO iq = 1, m
      DO ip = 1, n
        kmat ( ip, iq ) = 0.5_dp*SUM(wr(:)*dra(:,ip)*drb(:,iq) &
                             +wr(:)*REAL(l*(l+1),dp)*ra(:,ip)*rb(:,iq)/r(:)**2)
      END DO
    END DO

  END SUBROUTINE num_kinetic

!------------------------------------------------------------------------------
!
!  U(l,pq) = INT(u Ra(u) Rb(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_nuclear ( umat, ra, rb, r, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: umat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, rb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ip, iq, m, n

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) STOP

    DO iq = 1, m
      DO ip = 1, n
        umat ( ip, iq ) = SUM(wr(:)*ra(:,ip)*rb(:,iq)/r(:))
      END DO
    END DO

  END SUBROUTINE num_nuclear

!------------------------------------------------------------------------------
!
!  U(l,pq) = INT(u dRa(u) dRb(u))du + l(l+1) * INT(Ra(u) Rb(u) / u)du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_kinnuc ( umat, l, ra, dra, rb, drb, r, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: umat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, dra, rb, drb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ip, iq, m, n

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) STOP

    DO iq = 1, m
      DO ip = 1, n
        umat ( ip, iq ) = SUM(wr(:)*dra(:,ip)*drb(:,iq)/r(:) &
                            +wr(:)*REAL(l*(l+1),dp)*ra(:,ip)*rb(:,iq)/r(:)**3)
      END DO
    END DO

  END SUBROUTINE num_kinnuc

!------------------------------------------------------------------------------
!
!  U(l,pq) = INT(u erf(a*u) Ra(u) Rb(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_erf ( upmat, a, ra, rb, r, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: upmat
    REAL(KIND=dp), INTENT(IN)                :: a
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, rb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ip, iq, k, m, n

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( upmat, 1 ) .OR. m > SIZE ( upmat, 2 ) ) STOP

    DO iq = 1, m
      DO ip = 1, n
        upmat ( ip, iq ) = 0._dp
        DO k=1,SIZE(r)
          upmat ( ip, iq ) = upmat ( ip, iq ) + &
                             (wr(k)*ra(k,ip)*rb(k,iq)*erf(a*r(k))/r(k))
        ENDDO
      END DO
    END DO

  END SUBROUTINE num_erf

!------------------------------------------------------------------------------
!
!  Overlap with Projectors P(l,k,rc) for k=0,1,..
!
!  P(l,k,rc) = SQRT(2)/SQRT(Gamma[l+2k+1.5])/rc^(l+2k+1.5) r^(l+2k) exp[-0.5(r/rc)^2]
!
!  SP(l,k,p,rc) = INT(u^2 R(u) P(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_proj_ol ( spmat, l, ra, k, rc, r, wr )

    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: spmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                :: rc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ierr, ip, n
    REAL(KIND=dp)                            :: pf
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pro

    n =  SIZE ( ra, 2 )
    IF ( SIZE ( spmat ) < n ) STOP

    ALLOCATE(pro(n),STAT=ierr)
    IF ( ierr /= 0 ) STOP

    pf = SQRT(2._dp)/SQRT(gamma1(l+2*k+1))/rc**(l+2*k+1.5_dp)
    pro(:) = pf * r(:)**(l+2*k) * EXP(-0.5_dp*(r(:)/rc)**2)

    DO ip = 1, n
      spmat(ip) = SUM(wr(:)*pro(:)*ra(:,ip))
    END DO

    DEALLOCATE(pro,STAT=ierr)
    IF ( ierr /= 0 ) STOP

  END SUBROUTINE num_proj_ol

!------------------------------------------------------------------------------
!
!  Matrix elements for Gaussian potentials
!
!  V(k,rc) = (r/rc)^2k exp[-1/2(r/rc)^2]
!
!  VP(l,k,p+q,rc) = INT(u^2 V(u) Ra(u) Rb(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_gpot ( vpmat, k, rc,  ra, rb, r, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: vpmat
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                :: rc
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, rb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ierr, ip, iq, m, n
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: op

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( vpmat, 1 ) .OR. m > SIZE ( vpmat, 2 ) ) STOP

    ALLOCATE(op(n),STAT=ierr)
    IF ( ierr /= 0 ) STOP

    op(:) = (r(:)/rc)**(2*k) * EXP(-0.5_dp*(r(:)/rc)**2)

    DO iq = 1, m
      DO ip = 1, n
        vpmat ( ip, iq ) = SUM(wr(:)*ra(:,ip)*rb(:,iq)*op(:))
      END DO
    END DO

    DEALLOCATE(op,STAT=ierr)
    IF ( ierr /= 0 ) STOP

  END SUBROUTINE num_gpot

!------------------------------------------------------------------------------
!
!  G(l,k,pq) = <a|[r/rc]^2k|b>
!            = INT(u^2 [u/rc]^2k Ra(u) Rb(u))du
!
!------------------------------------------------------------------------------
  SUBROUTINE num_conf ( gmat, rc, k, ra, rb, r, wr )

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: gmat
    REAL(KIND=dp), INTENT(IN)                :: rc
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ra, rb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r, wr

    INTEGER                                  :: ierr, ip, iq, m, n
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: op

    n = SIZE ( ra, 2 )
    m = SIZE ( rb, 2 )

    IF ( n > SIZE ( gmat, 1 ) .OR. m > SIZE ( gmat, 2 ) ) STOP

    ALLOCATE(op(n),STAT=ierr)
    IF ( ierr /= 0 ) STOP

    op(:) = (r(:)/rc)**(2*k)

    DO iq = 1, m
      DO ip = 1, n
        gmat ( ip, iq ) = SUM(wr(:)*ra(:,ip)*rb(:,iq)*op(:))
      END DO
    END DO

    DEALLOCATE(op,STAT=ierr)
    IF ( ierr /= 0 ) STOP

  END SUBROUTINE num_conf

END MODULE ai_onecenter
