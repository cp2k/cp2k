!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
                                             gemint2
  USE ai_geminals_debug,               ONLY: grs,&
                                             init_grs
  USE hfx_energy_potential,            ONLY: coulomb4
  USE hfx_libint_wrapper,              ONLY: initialize_libint,&
                                             lib_int
  USE hfx_types,                       ONLY: hfx_potential_type,&
                                             hfx_screen_coeff_type
  USE input_constants,                 ONLY: do_hfx_potential_coulomb
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE orbital_pointers,                ONLY: coset,&
                                             deallocate_orbital_pointers,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE physcon,                         ONLY: bohr
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_coulomb_test'

  REAL(KIND=dp), PARAMETER :: threshold=1.0E-8_dp

  PRIVATE

  PUBLIC :: eri_test
! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE eri_test (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_test', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: lmax = 6

    CHARACTER*11, DIMENSION(0:lmax)          :: i2g
    CHARACTER*5, DIMENSION(0:lmax)           :: i2c
    CHARACTER*7, DIMENSION(0:lmax)           :: i2e
    CHARACTER*9, DIMENSION(0:lmax)           :: i2f
    INTEGER :: i, ii, istat, l, la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, ll, lmlib, n, na, nb, ncoa, ncob, ncoc, ncod, &
      npgfa, npgfb, npgfc, npgfd
    INTEGER(int_8)                           :: neris
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp) :: cartesian_estimate, eps_schwarz, log10_eps_schwarz, &
      log10_pmax, max_contraction_val, max_val2_set, perf, rab2, rac2, rcd2, &
      t, tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vac
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: v, za, zb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: pint
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab, iabc
    REAL(KIND=dp), DIMENSION(3)              :: ab, cd, ra, rb, rc, rd
    REAL(KIND=dp), DIMENSION(6)              :: a, b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rpgf, zeta, zetb, zetc, zetd
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: pgf1, pgf2, R1_pgf, R2_pgf
    TYPE(lib_int)                            :: lib

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
    npgfd = 2
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
         perf = ii*nco(l)**2 * REAL(npgfa*npgfc,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.0)') " Performance [integrals/s] ",i2c(l),perf
         DEALLOCATE(f,v,vac,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO

    END IF

#if defined (__LIBINT)
    ! Performance test of libint library
    IF ( iw>0 ) THEN

       WRITE(iw,'(//,A,/)') " Test of 2-Electron-4-Center Integrals (libint) "

       lmlib = MIN(lmax,3)
       CALL initialize_libint(lib,lmlib,error)

       ALLOCATE(R1_pgf(n,n),R2_pgf(n,n),pgf1(n,n),pgf2(n,n),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       R1_pgf(:,:)%x(1)=1._dp
       R1_pgf(:,:)%x(2)=1._dp
       R2_pgf(:,:)%x(1)=1._dp
       R2_pgf(:,:)%x(2)=1._dp
       pgf1(:,:)%x(1)=1._dp
       pgf1(:,:)%x(2)=1._dp
       pgf2(:,:)%x(1)=1._dp
       pgf2(:,:)%x(2)=1._dp
       potential_parameter%potential_type=do_hfx_potential_coulomb
       potential_parameter%omega=0._dp
       potential_parameter%scale_coulomb=1._dp
       potential_parameter%scale_longrange=1._dp
       potential_parameter%scale_gaussian=1._dp
       potential_parameter%cutoff_radius=100._dp
       potential_parameter%filename=""

       DO l=0,lmlib
         la_max=l
         la_min=0
         lb_max=l
         lb_min=0
         lc_max=l
         lc_min=0
         ld_max=l
         ld_min=0
         ll = ncoset(l)
         ncoa = npgfa*ll
         ncob = npgfb*ll
         ncoc = npgfc*ll
         ncod = npgfd*ll
         eps_schwarz=1.e-20_dp
         max_contraction_val=1._dp
         cartesian_estimate=1._dp
         rab2 = SUM((ra-rb)**2)
         rcd2 = SUM((rc-rd)**2)
         ab = ra - rb
         cd = rc - rd
         max_val2_set=1._dp
         log10_eps_schwarz=LOG(eps_schwarz)
         log10_pmax=1._dp
         ALLOCATE(pint(ncoa,ncob,ncoc,ncod),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         n=MAX(npgfa,npgfb,npgfc,npgfd)

         ii = MAX(200/(l+1)**4,1)
         tstart = m_walltime ( )
         DO i=1,ii
           CALL coulomb4(lib, ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd,&
                      la_min, la_max, lb_min, lb_max,&
                      lc_min, lc_max, ld_min, ld_max,&
                      ncoa, ncob, ncoc, ncod,&
                      zeta, zetb, zetc, zetd,&
                      pint,&
                      potential_parameter, eps_schwarz, max_contraction_val, &
                      cartesian_estimate, neris, rab2, rcd2, max_val2_set,&
                      log10_eps_schwarz, log10_pmax, R1_pgf, R2_pgf, pgf1, pgf2, ab, cd)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = ii*REAL(ncoa*ncob*ncoc*ncod,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.0)') " Performance [integrals/s] ",i2e(l),perf
         DEALLOCATE(pint,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO

       DEALLOCATE(R1_pgf,R2_pgf,pgf1,pgf2,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF
#endif

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
         ii = MAX(100/(l+1)**4,1)
         tstart = m_walltime ( )
         DO i=1,ii
           CALL gemint2(iab,za,na,l,l,a,zb,nb,l,l,b,error)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = ii*REAL(ll,dp)**4 * REAL(na*nb,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.0)') " Performance [integrals/s] ",i2g(l),perf
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
         ALLOCATE(iabc(npgfa*npgfc*nb,ll,ll,ll,ll),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ii = MAX(100/(l+1)**4,1)
         DO i=1,ii
           CALL g2gemint(iabc,l/2,npgfa,zeta,0,ra,l/2,npgfc,zetc,0,rc,l,l,zb,b,error)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = ii*REAL(ll,dp)**4 * REAL(npgfa*npgfc*nb,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.0)') " Performance [integrals/s] ",i2f(l),perf
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

    CALL deallocate_orbital_pointers()

  END SUBROUTINE eri_test

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
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: d1, dmax, intab, xa, xb, xc, &
                                                xd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab
    REAL(KIND=dp), DIMENSION(2, 2, 1)        :: za, zb
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d
    REAL(KIND=dp), DIMENSION(6)              :: ra, rb

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
    lb = 2
    lc = 1
    ld = 1

    lla = la+lc
    llb = lb+ld
    ALLOCATE(iab(1,ncoset(lla),ncoset(lla),ncoset(llb),ncoset(llb)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL gemint2(iab,za,1,lla,lla,ra,zb,1,llb,llb,rb,error)

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
                              WRITE(6,"(4I4,3F20.12)") iax,iay,ibx,iby,iab(1,iax,iay,ibx,iby),intab,d1
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
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: dmax, res1, res2, xa, xb, xc, &
                                                xd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iabcd
    REAL(KIND=dp), DIMENSION(3)              :: a, b, c, d

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
    ld = 1

    ALLOCATE(iabcd(1,ncoset(la),ncoset(lc),ncoset(lb),ncoset(ld)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL g4int(iabcd,la,1,(/xa/),0,a,lc,1,(/xc/),0,c,lb,1,(/xb/),0,b,ld,1,(/xd/),0,d,error)

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

END MODULE ai_coulomb_test

