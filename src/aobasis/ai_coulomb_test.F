!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Test of Electron Repulsion Routines (ERI)
!> \par History
!>      none
!> \author JGH (01.07.2009)
! *****************************************************************************
MODULE ai_coulomb_test

  USE ai_coulomb,                      ONLY: coulomb2
  USE ai_eri_debug,                    ONLY: init_os,&
                                             os
  USE ai_geminals,                     ONLY: g2gemint,&
                                             g4int,&
                                             gemint2,&
                                             gemint2_derivative
  USE ai_geminals_debug,               ONLY: grs,&
                                             init_grs
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE orbital_pointers,                ONLY: coset,&
                                             deallocate_orbital_pointers,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE physcon,                         ONLY: bohr
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_coulomb_test'

  REAL(KIND=dp), PARAMETER :: threshold=1.0E-6_dp

  PRIVATE

  PUBLIC :: eri_test
! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE eri_test (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_test', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: lmax = 6

    CHARACTER(LEN=11), DIMENSION(0:lmax)     :: i2g
    CHARACTER(LEN=5), DIMENSION(0:lmax)      :: i2c
    CHARACTER(LEN=7), DIMENSION(0:lmax)      :: i2e
    CHARACTER(LEN=9), DIMENSION(0:lmax)      :: i2f
    INTEGER                                  :: i, ii, istat, l, la_max, &
                                                la_min, lc_max, lc_min, ll, &
                                                n, na, nb, npgfa, npgfb, &
                                                npgfc, npgfd
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: perf, rac2, t, tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vac
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: v, za, zb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabc
    REAL(KIND=dp), DIMENSION(3)              :: ra, rb, rc, rd
    REAL(KIND=dp), DIMENSION(6)              :: a, b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rpgf, zeta, zetb, zetc, zetd

    failure = .FALSE.
    IF ( iw>0 ) WRITE(iw,'(/,A)') " Test of Electron Repulsion Integrals (ERI) "

    CALL init_orbital_pointers(lmax)

    i2c(0) = "(s|s)"
    i2c(1) = "(p|p)"
    i2c(2) = "(d|d)"
    i2c(3) = "(f|f)"
    i2c(4) = "(g|g)"
    i2c(5) = "(h|h)"
    i2c(6) = "(i|i)"

    i2g(0) = "[(ss)|(ss)]"
    i2g(1) = "[(pp)|(pp)]"
    i2g(2) = "[(dd)|(dd)]"
    i2g(3) = "[(ff)|(ff)]"
    i2g(4) = "[(gg)|(gg)]"
    i2g(5) = "[(hh)|(hh)]"
    i2g(6) = "[(ii)|(ii)]"

    i2f(0) = "[ss|(ss)]"
    i2f(1) = "[pp|(pp)]"
    i2f(2) = "[dd|(dd)]"
    i2f(3) = "[ff|(ff)]"
    i2f(4) = "[gg|(gg)]"
    i2f(5) = "[hh|(hh)]"
    i2f(6) = "[ii|(ii)]"

    i2e(0) = "(ss|ss)"
    i2e(1) = "(pp|pp)"
    i2e(2) = "(dd|dd)"
    i2e(3) = "(ff|ff)"
    i2e(4) = "(gg|gg)"
    i2e(5) = "(hh|hh)"
    i2e(6) = "(ii|ii)"

    npgfa = 4
    npgfb = 2
    npgfc = 4
    npgfd = 1
    n = MAX(npgfa,npgfb,npgfc,npgfd)

    ALLOCATE(zeta(npgfa),zetb(npgfb),zetc(npgfc),zetd(npgfd),rpgf(n),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    zeta(1:npgfa) = 0.5_dp
    zetb(1:npgfb) = 0.4_dp
    zetc(1:npgfc) = 0.3_dp
    zetd(1:npgfd) = 0.2_dp

    ra = (/0.0_dp,0.0_dp,0.0_dp/)
    rb = (/1.0_dp,0.0_dp,0.0_dp/)
    rc = (/0.0_dp,0.3_dp,0.3_dp/)
    rd = (/0.7_dp,0.2_dp,0.1_dp/)

    rac2 = SUM((ra-rc)**2)
    rpgf = 1.e10_dp

    ! Performance test of coulomb2 routine
    IF ( iw>0 ) THEN

       WRITE(iw,'(//,A,/)') " Test of 2-Electron-2-Center Integrals (coulomb2) "
       DO l=0,lmax
         la_max=l
         la_min=l
         lc_max=l
         lc_min=l
         ll = ncoset(l)
         ALLOCATE(f(0:2*l+2),v(npgfa*ll,npgfc*ll,2*l+1),vac(npgfa*ll,npgfc*ll),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         vac = 0._dp
         ii = MAX(100/(l+1)**2,1)
         tstart = m_walltime ( )
         DO i=1,ii
           CALL coulomb2(la_max,npgfa,zeta,rpgf,la_min,lc_max,npgfc,zetc,rpgf,lc_min,rc,rac2,vac,v,f)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = REAL(ii*nco(l)**2)*1.e-6_dp * REAL(npgfa*npgfc,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.3)') " Performance [Mintegrals/s] ",i2c(l),perf
         DEALLOCATE(f,v,vac,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO

    END IF

    ! Performance test of geminal routines
    IF ( iw > 0 ) THEN
       WRITE(iw,'(//,A,/)') " Test of 2-Electron-Geminal Integrals (gemint2) "
       a(1:3) = 0.5_dp * (ra+rc)
       a(4:6) = 0.5_dp * (ra-rc)
       b(1:3) = 0.5_dp * (rb+rd)
       b(4:6) = 0.5_dp * (rb-rd)
       na = npgfa*npgfc
       nb = npgfb*npgfd
       ALLOCATE(za(2,2,na),zb(2,2,nb),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       za(1,1,1:na) = 0.8_dp
       za(2,2,1:na) = 0.3_dp
       za(1,2,1:na) = 0.1_dp
       za(2,1,1:na) = 0.1_dp
       zb(1,1,1:nb) = 0.7_dp
       zb(2,2,1:nb) = 0.6_dp
       zb(1,2,1:nb) = 0.2_dp
       zb(2,1,1:nb) = 0.2_dp
       DO l=0,lmax/2
         ll = ncoset(l)
         ALLOCATE(iab(na*nb,ll,ll,ll,ll),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         iab=0._dp
         ii = MAX(100/(l+1)**4,1)
         tstart = m_walltime ( )
         DO i=1,ii
           CALL gemint2(iab,za,na,l,l,2*l,a,zb,nb,l,l,2*l,b,error)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = ii*REAL(ll,dp)**4*1.e-6_dp * REAL(na*nb,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.3)') " Performance [Mintegrals/s] ",i2g(l),perf
         DEALLOCATE(iab,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(za,zb,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! Performance test of geminal routines
    IF ( iw > 0 ) THEN
       WRITE(iw,'(//,A,/)') " Test of mixed Gaussian/Geminal 2-Electron Integrals (g2emint) "
       b(1:3) = 0.5_dp * (rb+rd)
       b(4:6) = 0.5_dp * (rb-rd)
       nb = npgfb*npgfd
       ALLOCATE(zb(2,2,nb),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       zb(1,1,1:nb) = 0.7_dp
       zb(2,2,1:nb) = 0.6_dp
       zb(1,2,1:nb) = 0.2_dp
       zb(2,1,1:nb) = 0.2_dp
       DO l=0,lmax/2
         tstart = m_walltime ( )
         ll = ncoset(l)
         ALLOCATE(iabc(npgfa*npgfc*nb,ll,ll,ll,ll,1),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         iabc=0._dp
         ii = MAX(100/(l+1)**4,1)
         DO i=1,ii
           CALL g2gemint(iabc,l/2,npgfa,zeta,ra,l/2,npgfc,zetc,rc,l,l,nb,zb,b,0,error)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = ii*REAL(ll,dp)**4*1.e-6 * REAL(npgfa*npgfc*nb,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.3)') " Performance [Mintegrals/s] ",i2f(l),perf
         DEALLOCATE(iabc,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(zb,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(zeta,zetb,zetc,zetd,rpgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL geminal_test1 (iw, error)

    CALL geminal_test2 (iw, error)

    CALL geminal_test3 (iw, error)

    CALL geminal_test4 (iw, error)

    CALL deallocate_orbital_pointers()

  END SUBROUTINE eri_test

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE geminal_test1 (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geminal_test1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia1, ia2, ia3, ia4, ia5, ia6, &
                                                iax, iay, ib1, ib2, ib3, ib4, &
                                                ib5, ib6, ibx, iby, istat, &
                                                la, lb, lc, ld, lla, llb, ma, &
                                                mb
    INTEGER, DIMENSION(6)                    :: an, bn
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d1, dmax, intab, xa, xb, xc, &
                                                xd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab
    REAL(KIND=dp), DIMENSION(2, 2, 1)        :: za, zb
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d
    REAL(KIND=dp), DIMENSION(6)              :: ra, rb

    failure = .FALSE.
    IF ( iw>0 ) WRITE(iw,'(//,A)') " Test of Geminal Repulsion Integrals (GRI) "

    xa = 0.783300000000
    xb = 1.239648746700
    xc = 0.548370000000
    xd = 0.111100000000

    a = (/0.329309000000,0.28408240000,0.28408240000/) * bohr
    b = (/0.983983000000,0.00453720000,0.00432740000/) * bohr
    c = (/0.032380000000,1.23470000000,0.11137400000/) * bohr
    d = (/0.000392438010,0.03491028400,0.87875300000/) * bohr

    za(1,1,1) = xa + xc
    za(1,2,1) = xa - xc
    za(2,1,1) = xa - xc
    za(2,2,1) = xa + xc
    zb(1,1,1) = xb + xd
    zb(1,2,1) = xb - xd
    zb(2,1,1) = xb - xd
    zb(2,2,1) = xb + xd
    ra(1:3) = 0.5_dp*(a+c)
    ra(4:6) = 0.5_dp*(a-c)
    rb(1:3) = 0.5_dp*(b+d)
    rb(4:6) = 0.5_dp*(b-d)

    CALL init_grs(za(:,:,1),zb(:,:,1),ra,rb)

    la = 2
    lb = 1
    lc = 0
    ld = 1

    lla = la+lc
    llb = lb+ld
    ALLOCATE(iab(1,ncoset(lla),ncoset(lla),ncoset(llb),ncoset(llb)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    iab = 0._dp

    CALL gemint2(iab,za,1,lla,lla,2*lla,ra,zb,1,llb,llb,2*llb,rb,error)

    dmax = 0._dp
    DO ma=0,la+lc
      DO mb=0,lb+ld

        DO ia1=0,ma
          DO ia2=0,ma-ia1
            DO ia3=0,ma-ia1-ia2
              DO ia4=0,ma-ia1-ia2-ia3
                DO ia5=0,ma-ia1-ia2-ia3-ia4
                  ia6=ma-ia1-ia2-ia3-ia4-ia5
                  an(1)=ia1
                  an(2)=ia2
                  an(3)=ia3
                  an(4)=ia4
                  an(5)=ia5
                  an(6)=ia6
                  iax=coset(ia1,ia2,ia3)
                  iay=coset(ia4,ia5,ia6)
                  DO ib1=0,mb
                    DO ib2=0,mb-ib1
                      DO ib3=0,mb-ib1-ib2
                        DO ib4=0,mb-ib1-ib2-ib3
                          DO ib5=0,mb-ib1-ib2-ib3-ib4
                            ib6=mb-ib1-ib2-ib3-ib4-ib5
                            bn(1)=ib1
                            bn(2)=ib2
                            bn(3)=ib3
                            bn(4)=ib4
                            bn(5)=ib5
                            bn(6)=ib6
                            ibx=coset(ib1,ib2,ib3)
                            iby=coset(ib4,ib5,ib6)
                            intab = grs(an,bn)
                            d1 = ABS(iab(1,iax,iay,ibx,iby)-intab)
                            dmax=MAX(d1,dmax)
                            IF(iw>0 .AND. d1 > 1.e-10) THEN
                              WRITE(*,"(4I4,3F20.12)") iax,iay,ibx,iby,iab(1,iax,iay,ibx,iby),intab,d1
                            END IF
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

    IF ( iw>0 ) WRITE(iw,'(/,A,T61,F20.12)') " Max. Error detected ",dmax

    DEALLOCATE(iab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE geminal_test1

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE geminal_test2 (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geminal_test2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia1, iax, iay, iaz, ib1, ibx, &
                                                iby, ibz, ic1, icx, icy, icz, &
                                                id1, idx, idy, idz, istat, &
                                                la, lb, lc, ld, ma, mb, mc, md
    INTEGER, DIMENSION(3)                    :: na, nb, nc, nd
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dmax, res1, res2, xa, xb, xc, &
                                                xd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iabcd
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d

    failure = .FALSE.
    IF ( iw>0 ) WRITE(iw,'(//,A)') " Test of ERI calculated from from Geminal Repulsion Integrals (GRI) "

    xa = 0.783300000000
    xb = 1.239648746700
    xc = 0.548370000000
    xd = 0.111100000000

    a = (/0.329309000000,0.28408240000,0.28408240000/) * bohr
    b = (/0.983983000000,0.00453720000,0.00432740000/) * bohr
    c = (/0.032380000000,1.23470000000,0.11137400000/) * bohr
    d = (/0.000392438010,0.03491028400,0.87875300000/) * bohr

    CALL init_os(xa,xb,xc,xd,a,b,c,d)

    la = 2
    lb = 2
    lc = 1
    ld = 0

    ALLOCATE(iabcd(1,ncoset(la),ncoset(lc),ncoset(lb),ncoset(ld)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    iabcd = 0._dp

    CALL g4int(iabcd,la,1,(/xa/),a,lc,1,(/xc/),c,lb,1,(/xb/),b,ld,1,(/xd/),d,error)

    dmax = 0._dp
    DO ma=0,la
      DO mb=0,lb
        DO mc=0,lc
          DO md=0,ld
            DO iax=0,ma
              DO iay=0,ma-iax
                iaz=ma-iax-iay
                na(1)=iax; na(2)=iay; na(3)=iaz
                ia1=coset(iax,iay,iaz)
                DO ibx=0,mb
                  DO iby=0,mb-ibx
                    ibz=mb-ibx-iby
                    nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                    ib1=coset(ibx,iby,ibz)
                    DO icx=0,mc
                      DO icy=0,mc-icx
                        icz=mc-icx-icy
                        nc(1)=icx; nc(2)=icy; nc(3)=icz
                        ic1=coset(icx,icy,icz)
                        DO idx=0,md
                          DO idy=0,md-idx
                            idz=md-idx-idy
                            nd(1)=idx; nd(2)=idy; nd(3)=idz
                            id1=coset(idx,idy,idz)
                            res1=os(na,nb,nc,nd)
                            res2=iabcd(1,ia1,ic1,ib1,id1)
                            dmax=MAX(dmax,ABS(res1-res2))
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


    IF ( iw>0 ) WRITE(iw,'(/,A,T61,F20.12)') " Max. Error detected ",dmax

    DEALLOCATE(iabcd,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE geminal_test2

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE geminal_test3 (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geminal_test3', &
      routineP = moduleN//':'//routineN

    INTEGER :: ia1, ia2, ia3, ia4, ia5, ia6, iax, iay, ib1, ib2, ib3, ib4, &
      ib5, ib6, ibx, iby, istat, k, la, lb, lc, ld, lla, llb, lmax, ma, mb
    INTEGER, DIMENSION(6)                    :: an, bn
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d1, delta, dmax, intab, intm, &
                                                intp, xa, xb, xc, xd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabd
    REAL(KIND=dp), DIMENSION(2, 2, 1)        :: za, zb
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d
    REAL(KIND=dp), DIMENSION(6)              :: ra, rb

    failure = .FALSE.
    IF ( iw>0 ) WRITE(iw,'(//,A)') " Test of Geminal Repulsion Integral Derivatives (GRI) "

    delta = 0.00001_dp

    xa = 0.783300000000
    xb = 1.239648746700
    xc = 0.548370000000
    xd = 0.111100000000

    a = (/0.329309000000,0.28408240000,0.28408240000/) * bohr
    b = (/0.983983000000,0.00453720000,0.00432740000/) * bohr
    c = (/0.032380000000,1.23470000000,0.11137400000/) * bohr
    d = (/0.000392438010,0.03491028400,0.87875300000/) * bohr

    za(1,1,1) = xa + xc
    za(1,2,1) = xa - xc
    za(2,1,1) = xa - xc
    za(2,2,1) = xa + xc
    zb(1,1,1) = xb + xd
    zb(1,2,1) = xb - xd
    zb(2,1,1) = xb - xd
    zb(2,2,1) = xb + xd
    ra(1:3) = 0.5_dp*(a+c)
    ra(4:6) = 0.5_dp*(a-c)
    rb(1:3) = 0.5_dp*(b+d)
    rb(4:6) = 0.5_dp*(b-d)

    la = 2
    lb = 1
    lc = 0
    ld = 1

    lla = la+lc
    llb = lb+ld

    lmax=lla+llb+1
    CALL init_orbital_pointers(lmax)

    ALLOCATE(iab(1,ncoset(lla+1),ncoset(lla+1),ncoset(llb),ncoset(llb)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iabd(1,ncoset(lla),ncoset(lla),ncoset(llb),ncoset(llb),7),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    iab=0._dp

    CALL gemint2(iab,za,1,lla+1,lla+1,2*lla+1,ra,zb,1,llb,llb,2*llb,rb,error)
    CALL gemint2_derivative(iab,iabd,za,1,lla,lla,1,llb,llb,1,.FALSE.,error)

    dmax = 0._dp
    DO ma=0,la+lc
      DO mb=0,lb+ld

        DO ia1=0,ma
          DO ia2=0,ma-ia1
            DO ia3=0,ma-ia1-ia2
              DO ia4=0,ma-ia1-ia2-ia3
                DO ia5=0,ma-ia1-ia2-ia3-ia4
                  ia6=ma-ia1-ia2-ia3-ia4-ia5
                  an(1)=ia1
                  an(2)=ia2
                  an(3)=ia3
                  an(4)=ia4
                  an(5)=ia5
                  an(6)=ia6
                  iax=coset(ia1,ia2,ia3)
                  iay=coset(ia4,ia5,ia6)
                  DO ib1=0,mb
                    DO ib2=0,mb-ib1
                      DO ib3=0,mb-ib1-ib2
                        DO ib4=0,mb-ib1-ib2-ib3
                          DO ib5=0,mb-ib1-ib2-ib3-ib4
                            ib6=mb-ib1-ib2-ib3-ib4-ib5
                            bn(1)=ib1
                            bn(2)=ib2
                            bn(3)=ib3
                            bn(4)=ib4
                            bn(5)=ib5
                            bn(6)=ib6
                            ibx=coset(ib1,ib2,ib3)
                            iby=coset(ib4,ib5,ib6)
                            CALL init_grs(za(:,:,1),zb(:,:,1),ra,rb)
                            intab = grs(an,bn)
                            d1 = ABS(iabd(1,iax,iay,ibx,iby,1)-intab)
                            dmax=MAX(d1,dmax)
                            IF(iw>0 .AND. d1 > 1.e-10) THEN
                              WRITE(*,"('I    ',4I4,3F20.12)") iax,iay,ibx,iby,iabd(1,iax,iay,ibx,iby,1),intab,d1
                            END IF
                            DO k=1,6
                              ra(k) = ra(k) + 0.5_dp*delta
                              CALL init_grs(za(:,:,1),zb(:,:,1),ra,rb)
                              intp = grs(an,bn)
                              ra(k) = ra(k) - delta
                              CALL init_grs(za(:,:,1),zb(:,:,1),ra,rb)
                              intm = grs(an,bn)
                              ra(k) = ra(k) + 0.5_dp*delta
                              intab = (intp-intm)/delta
                              d1 = ABS(iabd(1,iax,iay,ibx,iby,k+1)-intab)
                              dmax=MAX(d1,dmax)
                              IF(iw>0 .AND. d1 > 1.e-7) THEN
                                WRITE(*,"('dI/dR',I1,4x,4I4,3F18.10)") k,iax,iay,ibx,iby,&
                                      iabd(1,iax,iay,ibx,iby,k+1),intab,d1
                              END IF
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
    END DO

    IF ( iw>0 ) WRITE(iw,'(/,A,T61,F20.12)') " Max. Error detected ",dmax

    DEALLOCATE(iab,iabd,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE geminal_test3

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE geminal_test4 (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'geminal_test4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iax, iay, iaz, ic, icx, &
                                                icy, icz, ir, irx, iry, irz, &
                                                is, istat, isx, isy, isz, k, &
                                                la, lb, lc, ld, lla, llb, &
                                                lmax, ma, mb, mc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d1, da, db, dc, delta, dmax, &
                                                xa, xb, xc, xd, xr, xs
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabc1m, iabc1p, iabc2m, &
                                                iabc2p, iabc3m, iabc3p, iabcd
    REAL(KIND=dp), DIMENSION(2, 2, 1)        :: za, zb
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d
    REAL(KIND=dp), DIMENSION(6)              :: ra, rb

    failure = .FALSE.
    IF ( iw>0 ) WRITE(iw,'(//,A)') " Test of Geminal/Gaussian Repulsion Integral Derivatives (GRI) "

    delta = 0.00001_dp
    dmax  = 0.0_dp

    xa = 0.783300000000
    xb = 1.239648746700
    xc = 0.548370000000
    xd = 0.111100000000

    a = (/0.329309000000,0.28408240000,0.28408240000/) * bohr
    b = (/0.983983000000,0.00453720000,0.00432740000/) * bohr
    c = (/0.032380000000,1.23470000000,0.11137400000/) * bohr
    d = (/0.000392438010,0.03491028400,0.87875300000/) * bohr

    za(1,1,1) = xa + xc
    za(1,2,1) = xa - xc
    za(2,1,1) = xa - xc
    za(2,2,1) = xa + xc
    zb(1,1,1) = xb + xd
    zb(1,2,1) = xb - xd
    zb(2,1,1) = xb - xd
    zb(2,2,1) = xb + xd
    ra(1:3) = 0.5_dp*(a+c)
    ra(4:6) = 0.5_dp*(a-c)
    rb(1:3) = 0.5_dp*(b+d)
    rb(4:6) = 0.5_dp*(b-d)

    la = 2
    lb = 2
    lc = 1
    ld = 1

    lla = la+lc
    llb = lb+ld

    lmax=lla+llb+1
    CALL init_orbital_pointers(lmax)

    ALLOCATE(iabcd(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),7),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    iabcd = 0._dp

    CALL g2gemint(iabcd,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,1,error)

    ALLOCATE(iabc1p(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iabc1m(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iabc2p(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iabc2m(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iabc3p(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iabc3m(1,ncoset(la),ncoset(lc),ncoset(llb),ncoset(llb),1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO k=1,3

      iabc1p =0._dp
      iabc1m =0._dp
      iabc2p =0._dp
      iabc2m =0._dp
      iabc3p =0._dp
      iabc3m =0._dp

      a(k) = a(k) + 0.5_dp*delta
      CALL g2gemint(iabc1p,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,0,error)
      a(k) = a(k) - delta
      CALL g2gemint(iabc1m,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,0,error)
      a(k) = a(k) + 0.5_dp*delta
      iabc1p(:,:,:,:,:,:) = (iabc1p-iabc1m)/delta

      c(k) = c(k) + 0.5_dp*delta
      CALL g2gemint(iabc2p,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,0,error)
      c(k) = c(k) - delta
      CALL g2gemint(iabc2m,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,0,error)
      c(k) = c(k) + 0.5_dp*delta
      iabc2p(:,:,:,:,:,:) = (iabc2p-iabc2m)/delta

      rb(k) = rb(k) + 0.5_dp*delta
      CALL g2gemint(iabc3p,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,0,error)
      rb(k) = rb(k) - delta
      CALL g2gemint(iabc3m,la,1,(/xa/),a,lc,1,(/xb/),c,llb,llb,1,zb,rb,0,error)
      rb(k) = rb(k) + 0.5_dp*delta
      iabc3p(:,:,:,:,:,:) = (iabc3p-iabc3m)/delta

      DO ma=0,la
        DO mc=0,lc

          DO iax=0,ma
            DO iay=0,ma-iax
              iaz=ma-iax-iay
              ia=coset(iax,iay,iaz)
              DO icx=0,mc
                DO icy=0,mc-icx
                  icz=mc-icx-icy
                  ic=coset(icx,icy,icz)

                  DO mb=0,llb

                    DO irx=0,mb
                      DO iry=0,mb-irx
                        irz=mb-irx-iry
                        ir=coset(irx,iry,irz)
                        DO isx=0,mb
                          DO isy=0,mb-isx
                            isz=mb-isx-isy
                            is=coset(isx,isy,isz)

                            xr = iabcd(1,ia,ic,ir,is,k+1)
                            xs = iabcd(1,ia,ic,ir,is,k+4)

                            da = iabc1p(1,ia,ic,ir,is,1)
                            db = iabc3p(1,ia,ic,ir,is,1)
                            dc = iabc2p(1,ia,ic,ir,is,1)

                            d1 = ABS(da-0.5_dp*(xr+xs))
                            dmax=MAX(d1,dmax)
                            IF(iw>0 .AND. d1 > 1.e-7) THEN
                              WRITE(*,"('dI/da',I1,4x,4I4,3F18.10)") k,ia,ic,ir,is,0.5_dp*(xr+xs),da,d1
                            END IF
                            d1 = ABS(dc-0.5_dp*(xr-xs))
                            dmax=MAX(d1,dmax)
                            IF(iw>0 .AND. d1 > 1.e-7) THEN
                              WRITE(*,"('dI/dc',I1,4x,4I4,3F18.10)") k,ia,ic,ir,is,0.5_dp*(xr-xs),dc,d1
                            END IF
                            d1 = ABS(db-(-xr))
                            dmax=MAX(d1,dmax)
                            IF(iw>0 .AND. d1 > 1.e-7) THEN
                              WRITE(*,"('dI/db',I1,4x,4I4,3F18.10)") k,ia,ic,ir,is,-xr,db,d1
                            END IF

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

    IF ( iw>0 ) WRITE(iw,'(/,A,T61,F20.12)') " Max. Error detected ",dmax

    DEALLOCATE(iabcd,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(iabc1p,iabc1m,iabc2p,iabc2m,iabc3p,iabc3m,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE geminal_test4

! *****************************************************************************

END MODULE ai_coulomb_test

