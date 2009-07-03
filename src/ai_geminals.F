!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Coulomb integrals over Correlated Cartesian Gaussian 
!>        Geminals (electron repulsion integrals, ERIs).
!> \par Literature
!>      S.L. Saito and Y. Suzuki, J. Chem. Phys. 114 1115 (2001)
!> \par History
!>      none
!> \author Juerg Hutter (06.2009)
! *****************************************************************************
MODULE ai_geminals
  USE ai_geminals_utils,               ONLY: acoff,&
                                             matv
  USE f77_blas
  USE gamma,                           ONLY: fgamma
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "cp_common_uses.h"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_geminals'
  PRIVATE

  PUBLIC :: gemint2, g2gemint, g4int

  CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Correlated Cartesian Gaussian Geminals
!> \author  Juerg Hutter
!> \date    06.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE gemint2(iab,zeta,na,lra,lsa,a,zetb,nb,lrb,lsb,b,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: iab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zeta
    INTEGER, INTENT(IN)                      :: na, lra, lsa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zetb
    INTEGER, INTENT(IN)                      :: nb, lrb, lsb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemint2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, n, nr, ns
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: iab0

    n = na*nb
    nr = ncoset(lra+lrb)
    ns = ncoset(lsa+lsb)

    ALLOCATE (iab0(n,nr,ns),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! vertical recursion
    CALL vrr(iab0,zeta,na,lra,lsa,a,zetb,nb,lrb,lsb,b,error)

    ! horizontal recursion
    CALL hrr(iab,iab0,na,lra,lsa,a,nb,lrb,lsb,b,error)

    DEALLOCATE (iab0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE gemint2
! *****************************************************************************
  SUBROUTINE g2gemint(intabc,la_max,npgfa,zeta,la_min,a,lb_max,npgfb,zetb,lb_min,b,&
                      lr_max,ls_max,zetc,c,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: intabc
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: la_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    INTEGER, INTENT(IN)                      :: lr_max, ls_max
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zetc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'g2gemint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat, j, k, lrab, lsab, &
                                                na, nab, nb, nc, nrab, nrc, &
                                                nsab, nsc
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iabc
    REAL(KIND=dp), DIMENSION(6)              :: ab

    na = SIZE(zeta,1)
    nb = SIZE(zetb,1)
    nab = na*nb

    ALLOCATE (zab(2,2,nab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i=1,na
      DO j=1,nb
        k=(i-1)*nb+j
        zab(1,1,k) = zeta(i) + zetb(j)
        zab(2,2,k) = zeta(i) + zetb(j)
        zab(1,2,k) = zeta(i) - zetb(j)
        zab(2,1,k) = zeta(i) - zetb(j)
      END DO
    END DO

    ab(1:3) = 0.5_dp*(a(1:3)+b(1:3))
    ab(4:6) = 0.5_dp*(a(1:3)-b(1:3))

    nc = SIZE(zetc,3)
    lrab = la_max + lb_max
    lsab = la_max + lb_max

    nrab = ncoset(lrab)
    nsab = ncoset(lsab)
    nrc = ncoset(lr_max)
    nsc = ncoset(ls_max)

    ALLOCATE (iabc(nc*na*nb,nrab,nsab,nrc,nsc),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL gemint2(iabc,zab,nab,lrab,lsab,ab,zetc,nc,lr_max,ls_max,c,error)

    CALL gg_trans1(iabc,intabc,la_max,lb_max,error)

    DEALLOCATE (zab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (iabc,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE g2gemint
! *****************************************************************************
  SUBROUTINE g4int(intabcd,la_max,npgfa,zeta,la_min,a,lb_max,npgfb,zetb,lb_min,b,&
                   lc_max,npgfc,zetc,lc_min,c,ld_max,npgfd,zetd,ld_min,d,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: intabcd
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: la_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    INTEGER, INTENT(IN)                      :: lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: c
    INTEGER, INTENT(IN)                      :: ld_max, npgfd
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetd
    INTEGER, INTENT(IN)                      :: ld_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'g4int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat, j, k, lrab, lrcd, &
                                                lsab, lscd, ma, mb, na, nab, &
                                                nb, nc, ncd, nd, nrab, nrcd, &
                                                nsab, nscd
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zab, zcd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: gab, iab
    REAL(KIND=dp), DIMENSION(6)              :: ab, cd

    na = SIZE(zeta,1)
    nb = SIZE(zetb,1)
    nab = na*nb
    nc = SIZE(zetc,1)
    nd = SIZE(zetd,1)
    ncd = nc*nd

    ALLOCATE (zab(2,2,nab),zcd(2,2,ncd),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i=1,na
      DO j=1,nb
        k=(i-1)*nb+j
        zab(1,1,k) = zeta(i) + zetb(j)
        zab(2,2,k) = zeta(i) + zetb(j)
        zab(1,2,k) = zeta(i) - zetb(j)
        zab(2,1,k) = zeta(i) - zetb(j)
      END DO
    END DO
    DO i=1,nc
      DO j=1,nd
        k=(i-1)*nd+j
        zcd(1,1,k) = zetc(i) + zetd(j)
        zcd(2,2,k) = zetc(i) + zetd(j)
        zcd(1,2,k) = zetc(i) - zetd(j)
        zcd(2,1,k) = zetc(i) - zetd(j)
      END DO
    END DO

    ab(1:3) = 0.5_dp*(a(1:3)+b(1:3))
    ab(4:6) = 0.5_dp*(a(1:3)-b(1:3))
    cd(1:3) = 0.5_dp*(c(1:3)+d(1:3))
    cd(4:6) = 0.5_dp*(c(1:3)-d(1:3))

    lrab = la_max + lb_max
    lsab = la_max + lb_max
    lrcd = lc_max + ld_max
    lscd = lc_max + ld_max

    nrab = ncoset(lrab)
    nsab = ncoset(lsab)
    nrcd = ncoset(lrcd)
    nscd = ncoset(lscd)

    ALLOCATE (iab(na*nb*nc*nd,nrab,nsab,nrcd,nscd),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL gemint2(iab,zab,nab,lrab,lsab,ab,zcd,ncd,lrcd,lscd,cd,error)

    DEALLOCATE (zab,zcd,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ma = ncoset(la_max)
    mb = ncoset(lb_max)
    ALLOCATE (gab(na*nb*nc*nd,ma,mb,nrcd,nscd),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL gg_trans1(iab,gab,la_max,lb_max,error)

    DEALLOCATE (iab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL gg_trans2(gab,intabcd,lc_max,ld_max,error)

    DEALLOCATE (gab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE g4int
! *****************************************************************************
  SUBROUTINE vrr(iab,zeta,na,lra,lsa,a,zetb,nb,lrb,lsb,b,error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: iab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zeta
    INTEGER, INTENT(IN)                      :: na, lra, lsa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zetb
    INTEGER, INTENT(IN)                      :: nb, lrb, lsb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vrr', &
      routineP = moduleN//':'//routineN
    INTEGER, DIMENSION(3), PARAMETER :: ii1 = (/1,0,0/), ii2 = (/0,1,0/), &
      ii3 = (/0,0,1/), ii4 = (/1,0,0/), ii5 = (/0,1,0/), ii6 = (/0,0,1/)
    INTEGER, PARAMETER                       :: nmem = 1000000

    INTEGER :: i1, i2, ii, iint, ir, irx, iry, irz, is, istat, isx, isy, isz, &
      ix, ixx(1), j10, j20, k, k10, k20, ka, kb, kk, ll, llr, lls, lr, ls, m, &
      mm, n, ni, nr, ns
    INTEGER, DIMENSION(3)                    :: irm, irmm, irr, ism, ismm, iss
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: aai, abi, abj, ari, asi, bab, &
                                                d01m, d01p, d02m, d02p, d10m, &
                                                d10p, d11m, d11p, d20m, d20p, &
                                                eab, kab, rho, zdet
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pfab, t
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: f00
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: za, zb, zg, zgm
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: intm
    REAL(KIND=dp), DIMENSION(2)              :: v, w
    REAL(KIND=dp), DIMENSION(2, 2)           :: ee
    REAL(KIND=dp), DIMENSION(6)              :: ab, c, d, q

    lr = lra + lrb
    ls = lsa + lsb
    ll = lr + ls

    n = na*nb

    nr = ncoset(lr)
    ns = ncoset(ls)

    CPPrecondition(n<=SIZE(iab,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(nr<=SIZE(iab,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(ns<=SIZE(iab,3),cp_failure_level,routineP,error,failure)
    CPPrecondition(na<=SIZE(zeta,3),cp_failure_level,routineP,error,failure)
    CPPrecondition(nb<=SIZE(zetb,3),cp_failure_level,routineP,error,failure)

    ab = a - b

    ni = MAX(nmem/(nr*ns*(ll+1)),1)
    ni = MIN(ni,n)
    ALLOCATE (za(2,2,ni),zb(2,2,ni),zg(2,2,ni),zgm(2,2,ni),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (intm(ni,0:nr,0:ns,0:ll),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (t(ni),pfab(ni),f00(ni,0:ll),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    intm(:,0,0:,0:)=0._dp
    intm(:,0:,0,0:)=0._dp

    DO iint=1,n,ni
      i1 = iint
      i2 = MIN(n,i1+ni-1)
      ii = i2-i1+1
      t = 0._dp
      DO k=i1,i2
        kk = k-i1+1
        ka = (k-1)/nb+1
        kb = MOD(k-1,nb)+1
        za(:,:,kk) = zeta(:,:,ka)
        zb(:,:,kk) = zetb(:,:,kb)
        zg(:,:,kk) = za(:,:,kk) + zb(:,:,kk)
        zdet = zg(1,1,kk)*zg(2,2,kk) - zg(1,2,kk)*zg(2,1,kk)
        zgm(1,1,kk) = zg(2,2,kk)/zdet
        zgm(1,2,kk) = -zg(1,2,kk)/zdet
        zgm(2,1,kk) = -zg(2,1,kk)/zdet
        zgm(2,2,kk) = zg(1,1,kk)/zdet
        c = matv(za(:,:,kk),A)+matv(zb(:,:,kk),B)
        q = matv(zgm(:,:,kk),C)
        rho = zdet/zg(1,1,kk)
        t(kk) = rho*(q(4)*q(4) + q(5)*q(5) + q(6)*q(6))
        ee = MATMUL(za(:,:,kk),MATMUL(zgm(:,:,kk),zb(:,:,kk)))
        d = matv(ee,ab)
        eab = SUM(ab*d)
        kab = EXP(-eab)
        pfab(kk) = 8._dp*kab*(pi/zg(1,1,kk))**1.5_dp * pi/rho
      END DO

      CALL fgamma(ll,t,f00)

      DO m=0,ll
        intm(1:ii,1,1,m) = f00(1:ii,m)*pfab(1:ii)
      END DO

      ! vertical recursion [ar,as||0]
      DO llr = 0, lr
        DO lls = 0, MIN(llr,ls)
          IF(llr==0 .AND. lls==0) CYCLE
          DO irx=0,llr
            DO iry=0,llr-irx
              irz=llr-irx-iry
              ir=coset(irx,iry,irz)
              irr(1)=irx; irr(2)=iry; irr(3)=irz
              ixx = MAXLOC(irr)
              ix = ixx(1)
              CPPostcondition(irr(ix)>0,cp_failure_level,routineP,error,failure)
              irm = irr
              irm(ix) = irm(ix) - 1
              irmm = irr
              irmm(ix) = irmm(ix) - 2
              DO isx=0,lls
                DO isy=0,lls-isx
                  isz=lls-isx-isy
                  is=coset(isx,isy,isz)
                  iss(1)=isx; iss(2)=isy; iss(3)=isz
                  ism = iss
                  ism(ix) = ism(ix) - 1
                  j10=coset(irm(1),irm(2),irm(3))
                  j20=coset(irmm(1),irmm(2),irmm(3))
                  k10=coset(ism(1),ism(2),ism(3))
                  ari = REAL(irm(ix),dp)
                  asi = REAL(iss(ix),dp)
                  abi = ab(ix)
                  abj = ab(3+ix)
                  aai = a(3+ix)
                  DO mm=0,ll-llr-lls
                    ! increase r
                    DO k=1,ii
                      d10m = intm(k,j10,is,mm)
                      d10p = intm(k,j10,is,mm+1)
                      d20m = intm(k,j20,is,mm)
                      d20p = intm(k,j20,is,mm+1)
                      d11m = intm(k,j10,k10,mm)
                      d11p = intm(k,j10,k10,mm+1)
                      v(1) = abi*(d10p-d10m)
                      v(2) = abj*(d10p-d10m)
                      v = MATMUL(zb(:,:,k),v)
                      v(1) = v(1) + 0.5_dp*ari*(d20m-d20p)
                      v(2) = v(2) + 0.5_dp*asi*(d11m-d11p)
                      w = MATMUL(zgm(:,:,k),v)
                      bab = zb(1,1,k)*abi+zb(1,2,k)*abj-zg(1,2,k)*aai
                      w(1) = w(1) - bab*d10p/zg(1,1,k)
                      w(1) = w(1) + 0.5_dp*ari*d20p/zg(1,1,k)
                      intm(k,ir,is,mm) = w(1)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

        DO lls=llr+1,ls
          DO irx=0,llr
            DO iry=0,llr-irx
              irz=llr-irx-iry
              ir=coset(irx,iry,irz)
              irr(1)=irx; irr(2)=iry; irr(3)=irz
              DO isx=0,lls
                DO isy=0,lls-isx
                  isz=lls-isx-isy
                  is=coset(isx,isy,isz)
                  iss(1)=isx; iss(2)=isy; iss(3)=isz
                  ixx = MAXLOC(iss)
                  ix = ixx(1)
                  CPPostcondition(iss(ix)>0,cp_failure_level,routineP,error,failure)
                  irm = irr
                  irm(ix) = irm(ix) - 1
                  ism = iss
                  ism(ix) = ism(ix) - 1
                  ismm = iss
                  ismm(ix) = ismm(ix) - 2
                  j10=coset(irm(1),irm(2),irm(3))
                  k10=coset(ism(1),ism(2),ism(3))
                  k20=coset(ismm(1),ismm(2),ismm(3))
                  ari = REAL(irr(ix),dp)
                  asi = REAL(ism(ix),dp)
                  abi = ab(ix)
                  abj = ab(3+ix)
                  aai = a(3+ix)
                  DO mm=0,ll-llr-lls
                    ! increase s
                    DO k=1,ii
                      d01m = intm(k,ir,k10,mm)
                      d01p = intm(k,ir,k10,mm+1)
                      d02m = intm(k,ir,k20,mm)
                      d02p = intm(k,ir,k20,mm+1)
                      d11m = intm(k,j10,k10,mm)
                      d11p = intm(k,j10,k10,mm+1)
                      v(1) = abi*(d01p-d01m)
                      v(2) = abj*(d01p-d01m)
                      v = MATMUL(zb(:,:,k),v)
                      v(1) = v(1) + 0.5_dp*ari*(d11m-d11p)
                      v(2) = v(2) + 0.5_dp*asi*(d02m-d02p)
                      w = MATMUL(zgm(:,:,k),v)
                      w(2) = w(2) - aai*d01p
                      intm(k,ir,is,mm) = w(2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END DO

      DO is=1,ns
        DO ir=1,nr
          DO k=i1,i2
            kk=k-i1+1
            iab(k,ir,is) = intm(kk,ir,is,0)
          END DO
        END DO
      END DO

    END DO

    DEALLOCATE (za,zb,zg,zgm,t,pfab,f00,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (intm,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE vrr
! *****************************************************************************
  SUBROUTINE hrr(int,iab0,na,lra,lsa,a,nb,lrb,lsb,b,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: int
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: iab0
    INTEGER, INTENT(IN)                      :: na, lra, lsa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: nb, lrb, lsb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hrr', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii(6), iim(6), ira, irap, irax, iray, iraz, irb, irbm, irbx, &
      irby, irbz, isa, isap, isax, isay, isaz, isb, isbm, isbx, isby, isbz, &
      istat, ix, ixx(1), jj(6), jjm(6), mra, mrb, msa, msb, n, nra, nrab, &
      nrb, nsa, nsab, nsb
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: abx
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab

    nra = ncoset(lra)
    nsa = ncoset(lsa)
    nrb = ncoset(lrb)
    nsb = ncoset(lsb)
    nrab = ncoset(lra+lrb)
    nsab = ncoset(lsa+lsb)

    CPPrecondition(na*nb<=SIZE(int,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(nra<=SIZE(int,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsa<=SIZE(int,3),cp_failure_level,routineP,error,failure)
    CPPrecondition(nrb<=SIZE(int,4),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsb<=SIZE(int,5),cp_failure_level,routineP,error,failure)

    CPPrecondition(na*nb<=SIZE(iab0,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(nrab<=SIZE(iab0,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsab<=SIZE(iab0,3),cp_failure_level,routineP,error,failure)

    n = na*nb
    ALLOCATE (iab(n,nrab,nsab,nrb,nsb),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO mrb=0,lrb
      DO msb=0,lsb
        IF (mrb+msb == 0) THEN
          iab(1:n,1:nrab,1:nsab,1,1) = iab0(1:n,1:nrab,1:nsab)
        ELSE
          DO irbx=0,mrb
            DO irby=0,mrb-irbx
              irbz=mrb-irbx-irby
              irb=coset(irbx,irby,irbz)
              ii(1)=irbx; ii(2)=irby; ii(3)=irbz
              DO isbx=0,msb
                DO isby=0,msb-isbx
                  isbz=msb-isbx-isby
                  isb=coset(isbx,isby,isbz)
                  ii(4)=isbx; ii(5)=isby; ii(6)=isbz
                  ixx = MAXLOC(ii)
                  ix = ixx(1)
                  abx = a(ix) - b(ix)
                  iim = ii
                  iim(ix) = iim(ix) - 1
                  irbm=coset(iim(1),iim(2),iim(3))
                  isbm=coset(iim(4),iim(5),iim(6))

                  DO mra=0,lra+lrb-mrb
                    DO msa=0,lsa+lsb-msb
                      DO irax=0,mra
                        DO iray=0,mra-irax
                          iraz=mra-irax-iray
                          ira=coset(irax,iray,iraz)
                          jj(1)=irax; jj(2)=iray; jj(3)=iraz
                          DO isax=0,msa
                            DO isay=0,msa-isax
                              isaz=msa-isax-isay
                              isa=coset(isax,isay,isaz)
                              jj(4)=isax; jj(5)=isay; jj(6)=isaz
                              jjm = jj
                              jjm(ix) = jjm(ix) + 1
                              irap=coset(jjm(1),jjm(2),jjm(3))
                              isap=coset(jjm(4),jjm(5),jjm(6))
          
                              iab(1:n,ira,isa,irb,isb) = iab(1:n,irap,isap,irbm,isbm) + &
                                                         abx * iab(1:n,ira,isa,irbm,isbm)

                            END DO
                          END DO
                        END DO
                      END DO
                    END DO
                  END DO

                END DO
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO

    INT(1:n,1:nra,1:nsa,1:nrb,1:nsb) = iab(1:n,1:nra,1:nsa,1:nrb,1:nsb)

    DEALLOCATE (iab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE hrr
! *****************************************************************************
  SUBROUTINE gg_trans1(in,iout,la,lb,error)
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(IN)                             :: in
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(OUT)                            :: iout
    INTEGER, INTENT(IN)                      :: la, lb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gg_trans1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, i3, istat, iu, iv, &
                                                iw, j1, j2, j3, k1, k2, l1, &
                                                l2, m1, m2, mk, mr, ms, nk, &
                                                nmax, nr, ns
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ac
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: a123

    nmax=MAX(la,lb)

    ALLOCATE (ac(0:nmax,0:nmax,0:2*nmax),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL acoff(ac,nmax)

    mk = SIZE(in,1)
    mr = SIZE(in,4)
    ms = SIZE(in,5)

    iout=0._dp
    DO l1=0,la
      DO l2=0,lb

        DO i1=0,l1
          DO i2=0,l1-i1
            i3=l1-i1-i2
            k1=coset(i1,i2,i3)
            DO j1=0,l2
              DO j2=0,l2-j1
                j3=l2-j1-j2
                k2=coset(j1,j2,j3)
                !transform integrals (ij|.)
                DO iu=0,i1+j1
                  DO iv=0,i2+j2
                    DO iw=0,i3+j3
                      m1=coset(i1+j1-iu,i2+j2-iv,i3+j3-iw)
                      m2=coset(iu,iv,iw)
                      a123 = ac(i1,j1,iu)*ac(i2,j2,iv)*ac(i3,j3,iw)
                      DO ns=1,ms
                        DO nr=1,mr
                          DO nk=1,mk
                            iout(nk,k1,k2,nr,ns) = iout(nk,k1,k2,nr,ns) + a123*in(nk,m1,m2,nr,ns)
                          END DO
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END DO
    END DO

    DEALLOCATE (ac,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE gg_trans1
! *****************************************************************************
  SUBROUTINE gg_trans2(in,iout,la,lb,error)
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(IN)                             :: in
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(OUT)                            :: iout
    INTEGER, INTENT(IN)                      :: la, lb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gg_trans2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, i3, istat, iu, iv, &
                                                iw, j1, j2, j3, k1, k2, l1, &
                                                l2, m1, m2, mk, mr, ms, nk, &
                                                nmax, nr, ns
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ac
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: a123

    nmax=MAX(la,lb)

    ALLOCATE (ac(0:nmax,0:nmax,0:2*nmax),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL acoff(ac,nmax)

    mk = SIZE(in,1)
    mr = SIZE(in,2)
    ms = SIZE(in,3)

    iout=0._dp
    DO l1=0,la
      DO l2=0,lb

        DO i1=0,l1
          DO i2=0,l1-i1
            i3=l1-i1-i2
            k1=coset(i1,i2,i3)
            DO j1=0,l2
              DO j2=0,l2-j1
                j3=l2-j1-j2
                k2=coset(j1,j2,j3)
                !transform integrals (.|ij)
                DO iu=0,i1+j1
                  DO iv=0,i2+j2
                    DO iw=0,i3+j3
                      m1=coset(i1+j1-iu,i2+j2-iv,i3+j3-iw)
                      m2=coset(iu,iv,iw)
                      a123 = ac(i1,j1,iu)*ac(i2,j2,iv)*ac(i3,j3,iw)
                      DO ns=1,ms
                        DO nr=1,mr
                          DO nk=1,mk
                            iout(nk,nr,ns,k1,k2) = iout(nk,nr,ns,k1,k2) + a123*in(nk,nr,ns,m1,m2)
                          END DO
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END DO
    END DO

    DEALLOCATE (ac,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE gg_trans2
! *****************************************************************************

END MODULE ai_geminals
