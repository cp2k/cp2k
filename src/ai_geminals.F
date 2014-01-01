!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE gamma,                           ONLY: fgamma=>fgamma_1
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_geminals'
  PRIVATE

  PUBLIC :: gemint2, g2gemint, g4int, gemint2_derivative
  PUBLIC :: gemnorm

  CONTAINS

! *****************************************************************************
!> \brief   Calculation of the norm of Correlated Cartesian Gaussian Geminals
!> \author  Juerg Hutter
!> \date    08.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE gemnorm(anorm,zeta,na,lra,lsa,error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: anorm
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zeta
    INTEGER, INTENT(IN)                      :: na, lra, lsa
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemnorm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ir, irm, irx, iry, &
                                                irz, is, ism, istat, isx, &
                                                isy, isz, lr, ls
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pi3
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zdet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zg

    failure = .FALSE.
    CALL timeset(routineN,handle)

    anorm = 0._dp

    ALLOCATE (zg(2,2,na),zdet(na),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    pi3 = pi**3
    DO i=1,na
      zdet(i) = zeta(1,1,i)*zeta(2,2,i) - zeta(1,2,i)*zeta(2,1,i)
      zg(1,1,i) = zeta(2,2,i)/zdet(i)
      zg(1,2,i) = -zeta(1,2,i)/zdet(i)
      zg(2,1,i) = -zeta(2,1,i)/zdet(i)
      zg(2,2,i) = zeta(1,1,i)/zdet(i)
      anorm(i,1,1) = pi3/SQRT(zdet(i))
    END DO

    DO lr=0,lra,2
      DO irx=0,lr,2
        DO iry=0,lr-irx,2
          irz=lr-irx-iry
          ir=coset(irx,iry,irz)
          DO ls=0,lsa,2
            IF (lr+ls == 0) CYCLE
            DO isx=0,ls,2
              DO isy=0,ls-isx,2
                isz=ls-isx-isy
                is=coset(isx,isy,isz)
                IF (lr >= ls) THEN
                  IF ( irx >= iry .AND. irx >= irz ) THEN
                    irm=coset(irx-2,iry,irz)
                    anorm(1:na,ir,is)=0.5_dp*zg(1,1,1:na)*anorm(1:na,irm,is)
                    IF (isx > 0) THEN
                      ism=coset(isx-2,isy,isz)
                      anorm(1:na,ir,is)=anorm(1:na,ir,is) + 0.5_dp*zg(1,2,1:na)*anorm(1:na,ir,ism)
                    END IF
                  ELSE IF ( iry >= irz ) THEN
                    irm=coset(irx,iry-2,irz)
                    anorm(1:na,ir,is)=0.5_dp*zg(1,1,1:na)*anorm(1:na,irm,is)
                    IF (isy > 0) THEN
                      ism=coset(isx,isy-2,isz)
                      anorm(1:na,ir,is)=anorm(1:na,ir,is) + 0.5_dp*zg(1,2,1:na)*anorm(1:na,ir,ism)
                    END IF
                  ELSE
                    irm=coset(irx,iry,irz-2)
                    anorm(1:na,ir,is)=0.5_dp*zg(1,1,1:na)*anorm(1:na,irm,is)
                    IF (isz > 0) THEN
                      ism=coset(isx,isy,isz-2)
                      anorm(1:na,ir,is)=anorm(1:na,ir,is) + 0.5_dp*zg(1,2,1:na)*anorm(1:na,ir,ism)
                    END IF
                  END IF
                ELSE
                  IF ( isx >= isy .AND. isx >= isz ) THEN
                    ism=coset(isx-2,isy,isz)
                    anorm(1:na,ir,is)=0.5_dp*zg(2,2,1:na)*anorm(1:na,ir,ism)
                    IF (irx > 0) THEN
                      irm=coset(irx-2,iry,irz)
                      anorm(1:na,ir,is)=anorm(1:na,ir,is) + 0.5_dp*zg(2,1,1:na)*anorm(1:na,irm,is)
                    END IF
                  ELSE IF ( isy >= isz ) THEN
                    ism=coset(isx,isy-2,isz)
                    anorm(1:na,ir,is)=0.5_dp*zg(2,2,1:na)*anorm(1:na,ir,ism)
                    IF (iry > 0) THEN
                      irm=coset(irx,iry-2,irz)
                      anorm(1:na,ir,is)=anorm(1:na,ir,is) + 0.5_dp*zg(2,1,1:na)*anorm(1:na,irm,is)
                    END IF
                  ELSE
                    ism=coset(isx,isy,isz-2)
                    anorm(1:na,ir,is)=0.5_dp*zg(2,2,1:na)*anorm(1:na,ir,ism)
                    IF (irz > 0) THEN
                      irm=coset(irx,iry,irz-2)
                      anorm(1:na,ir,is)=anorm(1:na,ir,is) + 0.5_dp*zg(2,1,1:na)*anorm(1:na,irm,is)
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (zg,zdet,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE gemnorm

! *****************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Correlated Cartesian Gaussian Geminals
!> \author  Juerg Hutter
!> \date    06.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE gemint2(iab,zeta,na,lra,lsa,lta,a,zetb,nb,lrb,lsb,ltb,b,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: iab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zeta
    INTEGER, INTENT(IN)                      :: na, lra, lsa, lta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zetb
    INTEGER, INTENT(IN)                      :: nb, lrb, lsb, ltb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemint2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, n, nr, ns
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: iab0

    failure = .FALSE.
    CALL timeset(routineN,handle)

    n = na*nb
    nr = ncoset(lra+lrb)
    ns = ncoset(lsa+lsb)

    ALLOCATE (iab0(n,nr,ns),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! vertical recursion
    CALL vrr(iab0,zeta,na,lra,lsa,lta,a,zetb,nb,lrb,lsb,ltb,b,error)

    ! horizontal recursion
    CALL hrr(iab,iab0,na,lra,lsa,lta,a,nb,lrb,lsb,ltb,b,error)

    DEALLOCATE (iab0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE gemint2

! *****************************************************************************

  SUBROUTINE gemint2_derivative(iab,iabd,zeta,na,lra,lsa,nb,lrb,lsb,nder,r_only,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: iab
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :, :), &
      INTENT(INOUT)                          :: iabd
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zeta
    INTEGER, INTENT(IN)                      :: na, lra, lsa, nb, lrb, lsb, &
                                                nder
    LOGICAL, INTENT(IN)                      :: r_only
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gemint2_derivative', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, ij, ira, iram, irap, irax, iray, iraz, irb, irbx, irby, &
      irbz, isa, isam, isap, isax, isay, isaz, isb, isbx, isby, isbz, k, mra, &
      mrb, msa, msb
    INTEGER, DIMENSION(6)                    :: ii, jj
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ari, asi, za11, za12, za22

    failure = .FALSE.
    CPPrecondition(nder>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nder<=1,cp_failure_level,routineP,error,failure)

    DO mrb=0,lrb
      DO msb=0,lsb
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

                DO mra=0,lra
                  DO msa=0,lsa
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

                            DO i=1,na*nb
                              iabd(i,ira,isa,irb,isb,1) = iab(i,ira,isa,irb,isb)
                            END DO
                            IF (nder > 0) THEN
                              DO k=1,3
                                jj(k)=jj(k)-1
                                iram=coset(jj(1),jj(2),jj(3))
                                jj(1)=irax; jj(2)=iray; jj(3)=iraz
                                jj(k)=jj(k)+1
                                irap=coset(jj(1),jj(2),jj(3))
                                jj(k+3)=jj(k+3)+1
                                isap=coset(jj(4),jj(5),jj(6))
                                jj(1)=irax; jj(2)=iray; jj(3)=iraz
                                jj(4)=isax; jj(5)=isay; jj(6)=isaz
                                ari = -REAL(jj(k),dp)
                                DO i=1,na*nb
                                  ij=MOD(i-1,na)+1
                                  za11=2._dp*zeta(1,1,ij)
                                  za12=2._dp*zeta(1,2,ij)
                                  iabd(i,ira,isa,irb,isb,k+1) = ari*iab(i,iram,isa,irb,isb) +&
                                    za11*iab(i,irap,isa,irb,isb) + za12*iab(i,ira,isap,irb,isb)
                                END DO
                              END DO
                              IF (.NOT. r_only) THEN
                                DO k=1,3
                                  jj(k+3)=jj(k+3)-1
                                  isam=coset(jj(4),jj(5),jj(6))
                                  jj(4)=isax; jj(5)=isay; jj(6)=isaz
                                  jj(k+3)=jj(k+3)+1
                                  isap=coset(jj(4),jj(5),jj(6))
                                  jj(4)=isax; jj(5)=isay; jj(6)=isaz
                                  jj(k)=jj(k)+1
                                  irap=coset(jj(1),jj(2),jj(3))
                                  jj(1)=irax; jj(2)=iray; jj(3)=iraz
                                  asi = -REAL(jj(k+3),dp)
                                  DO i=1,na*nb
                                    ij=MOD(i-1,na)+1
                                    za22=2._dp*zeta(2,2,ij)
                                    za12=2._dp*zeta(1,2,ij)
                                    iabd(i,ira,isa,irb,isb,k+4) = asi*iab(i,ira,isam,irb,isb) +&
                                      za22*iab(i,ira,isap,irb,isb) + za12*iab(i,irap,isa,irb,isb)
                                  END DO
                                END DO
                              END IF
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


  END SUBROUTINE gemint2_derivative

! *****************************************************************************
  SUBROUTINE g2gemint(intabc,la_max,npgfa,zeta,a,lb_max,npgfb,zetb,b,&
                      lr_max,ls_max,ngemc,zetc,c,nderivative,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :, :), &
      INTENT(INOUT)                          :: intabc
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, a
    INTEGER, INTENT(IN)                      :: lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, b
    INTEGER, INTENT(IN)                      :: lr_max, ls_max, ngemc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zetc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: c
    INTEGER, INTENT(IN)                      :: nderivative
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'g2gemint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, istat, j, k, lrab, &
                                                lsab, ltab, ltc, na, nab, nb, &
                                                nc, nrab, nrc, nsab, nsc
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iabc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabcd
    REAL(KIND=dp), DIMENSION(6)              :: ab

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CPPrecondition(nderivative>=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nderivative<=1,cp_failure_level,routineP,error,failure)

    na = npgfa
    nb = npgfb
    nab = na*nb

    ALLOCATE (zab(2,2,nab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO j=1,nb
      DO i=1,na
        k=(j-1)*na+i
        zab(1,1,k) = zeta(i) + zetb(j)
        zab(2,2,k) = zeta(i) + zetb(j)
        zab(1,2,k) = zeta(i) - zetb(j)
        zab(2,1,k) = zeta(i) - zetb(j)
      END DO
    END DO

    ab(1:3) = 0.5_dp*(a(1:3)+b(1:3))
    ab(4:6) = 0.5_dp*(a(1:3)-b(1:3))

    nc = ngemc
    lrab = la_max + lb_max
    lsab = la_max + lb_max
    ltab = la_max + lb_max
    ltc  = lr_max + ls_max

    IF (nderivative==0) THEN

      nrab = ncoset(lrab)
      nsab = ncoset(lsab)
      nrc = ncoset(lr_max)
      nsc = ncoset(ls_max)

      ALLOCATE (iabc(nc*na*nb,nrab,nsab,nrc,nsc),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      iabc = 0._dp

      CALL gemint2(iabc,zab,nab,lrab,lsab,ltab,ab,zetc,nc,lr_max,ls_max,ltc,c,error)

      CALL gg_trans1(iabc,intabc(:,:,:,:,:,1),la_max,lb_max,error)

    ELSEIF (nderivative==1) THEN

      nrab = ncoset(lrab+1)
      nsab = ncoset(lsab+1)
      nrc = ncoset(lr_max)
      nsc = ncoset(ls_max)

      ALLOCATE (iabc(nc*na*nb,nrab,nsab,nrc,nsc),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      iabc = 0._dp

      CALL gemint2(iabc,zab,nab,lrab+1,lsab+1,ltab+1,ab,zetc,nc,lr_max,ls_max,ltc,c,error)

      nrab = ncoset(lrab)
      nsab = ncoset(lsab)
      ALLOCATE (iabcd(nc*na*nb,nrab,nsab,nrc,nsc,7),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL gemint2_derivative(iabc,iabcd,zab,nab,lrab,lsab,nc,lr_max,ls_max,nderivative,.FALSE.,error)

      DO k=1,7
        CALL gg_trans1(iabcd(:,:,:,:,:,k),intabc(:,:,:,:,:,k),la_max,lb_max,error)
      END DO

      DEALLOCATE (iabcd,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ELSE

      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)

    END IF

    DEALLOCATE (zab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (iabc,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE g2gemint
! *****************************************************************************
  SUBROUTINE g4int(intabcd,la_max,npgfa,zeta,a,lb_max,npgfb,zetb,b,&
                   lc_max,npgfc,zetc,c,ld_max,npgfd,zetd,d,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: intabcd
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, a
    INTEGER, INTENT(IN)                      :: lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, b
    INTEGER, INTENT(IN)                      :: lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc, c
    INTEGER, INTENT(IN)                      :: ld_max, npgfd
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetd, d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'g4int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat, j, k, lrab, lrcd, &
                                                lsab, lscd, ltab, ltcd, ma, &
                                                mb, na, nab, nb, nc, ncd, nd, &
                                                nrab, nrcd, nsab, nscd
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zab, zcd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: gab, iab
    REAL(KIND=dp), DIMENSION(6)              :: ab, cd

    failure = .FALSE.
    na = npgfa
    nb = npgfb
    nab = na*nb
    nc = npgfc
    nd = npgfd
    ncd = nc*nd

    ALLOCATE (zab(2,2,nab),zcd(2,2,ncd),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO j=1,nb
      DO i=1,na
        k=(j-1)*na+i
        zab(1,1,k) = zeta(i) + zetb(j)
        zab(2,2,k) = zeta(i) + zetb(j)
        zab(1,2,k) = zeta(i) - zetb(j)
        zab(2,1,k) = zeta(i) - zetb(j)
      END DO
    END DO
    DO j=1,nd
      DO i=1,nc
        k=(j-1)*nc+i
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
    ltab = la_max + lb_max
    lrcd = lc_max + ld_max
    lscd = lc_max + ld_max
    ltcd = lc_max + ld_max

    nrab = ncoset(lrab)
    nsab = ncoset(lsab)
    nrcd = ncoset(lrcd)
    nscd = ncoset(lscd)

    ALLOCATE (iab(na*nb*nc*nd,nrab,nsab,nrcd,nscd),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    iab = 0._dp

    CALL gemint2(iab,zab,nab,lrab,lsab,ltab,ab,zcd,ncd,lrcd,lscd,ltcd,cd,error)

    DEALLOCATE (zab,zcd,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ma = ncoset(la_max)
    mb = ncoset(lb_max)
    ALLOCATE (gab(na*nb*nc*nd,ma,mb,nrcd,nscd),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    gab = 0._dp

    CALL gg_trans1(iab,gab,la_max,lb_max,error)

    DEALLOCATE (iab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL gg_trans2(gab,intabcd,lc_max,ld_max,error)

    DEALLOCATE (gab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE g4int
! *****************************************************************************
  SUBROUTINE vrr(iab,zeta,na,lra,lsa,lta,a,zetb,nb,lrb,lsb,ltb,b,error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: iab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zeta
    INTEGER, INTENT(IN)                      :: na, lra, lsa, lta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zetb
    INTEGER, INTENT(IN)                      :: nb, lrb, lsb, ltb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vrr', &
      routineP = moduleN//':'//routineN
    INTEGER, DIMENSION(3), PARAMETER :: ii1 = (/1,0,0/), ii2 = (/0,1,0/), &
      ii3 = (/0,0,1/), ii4 = (/1,0,0/), ii5 = (/0,1,0/), ii6 = (/0,0,1/)
    INTEGER, PARAMETER                       :: nmem = 1000000

    INTEGER :: handle, i1, i2, ii, iint, ir, irx, iry, irz, is, istat, isx, &
      isy, isz, ix, ixx(1), j10, j20, k, k10, k20, ka, kb, kk, ll, llr, lls, &
      lr, ls, lsmax, lt, m, mm, n, ni, nr, ns
    INTEGER, DIMENSION(3)                    :: irm, irmm, irr, ism, ismm, iss
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aai, abi, abj, ari, asi, ba, &
                                                bab, bb, d01m, d01p, d02m, &
                                                d02p, d10m, d10p, d11m, d11p, &
                                                d20m, d20p, eab, kab, rho, &
                                                zdet
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

    failure = .FALSE.
    CALL timeset(routineN,handle)

    lr = lra + lrb
    ls = lsa + lsb
    lt = lta + ltb
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
        lsmax = MIN(llr,ls)
        DO lls = 0, lsmax
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
                  ari = 0.5_dp*REAL(irm(ix),dp)
                  asi = 0.5_dp*REAL(iss(ix),dp)
                  abi = ab(ix)
                  abj = ab(3+ix)
                  aai = a(3+ix)
                  IF ( ii==1 ) THEN
                    ba   = (zb(1,1,1)*abi+zb(1,2,1)*abj)
                    bb   = (zb(2,1,1)*abi+zb(2,2,1)*abj)
                    bab  = ba-zg(1,2,1)*aai
                    eab  = 1._dp/zg(1,1,1)
                    DO mm=0,ll-llr-lls
                      d10m = intm(1,j10,is,mm)
                      d10p = intm(1,j10,is,mm+1)
                      d20m = intm(1,j20,is,mm)
                      d20p = intm(1,j20,is,mm+1)
                      d11m = intm(1,j10,k10,mm)
                      d11p = intm(1,j10,k10,mm+1)
                      v(1) = ba*(d10p-d10m) + ari*(d20m-d20p)
                      v(2) = bb*(d10p-d10m) + asi*(d11m-d11p)
                      intm(1,ir,is,mm) = zgm(1,1,1)*v(1)+zgm(1,2,1)*v(2) - (bab*d10p-ari*d20p)*eab
                    END DO
                  ELSE
                    DO mm=0,ll-llr-lls
                      ! increase r
                      DO k=1,ii
                        d10m = intm(k,j10,is,mm)
                        d10p = intm(k,j10,is,mm+1)
                        d20m = intm(k,j20,is,mm)
                        d20p = intm(k,j20,is,mm+1)
                        d11m = intm(k,j10,k10,mm)
                        d11p = intm(k,j10,k10,mm+1)
                        bab  = (zb(1,1,k)*abi+zb(1,2,k)*abj)-zg(1,2,k)*aai
                        v(1) = (zb(1,1,k)*abi+zb(1,2,k)*abj)*(d10p-d10m) + ari*(d20m-d20p)
                        v(2) = (zb(2,1,k)*abi+zb(2,2,k)*abj)*(d10p-d10m) + asi*(d11m-d11p)
                        intm(k,ir,is,mm) = zgm(1,1,k)*v(1)+zgm(1,2,k)*v(2) - (bab*d10p-ari*d20p)/zg(1,1,k)
                      END DO
                    END DO
                  END IF
                END DO
              END DO
            END DO
          END DO
        END DO

        lsmax = ls
        DO lls=llr+1,lsmax
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

    CALL timestop(handle)

  END SUBROUTINE vrr
! *****************************************************************************
  SUBROUTINE hrr(intab,iab0,na,lra,lsa,lta,a,nb,lrb,lsb,ltb,b,error)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: intab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: iab0
    INTEGER, INTENT(IN)                      :: na, lra, lsa, lta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: nb, lrb, lsb, ltb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hrr', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ii(6), iim(6), ira, irap, irax, iray, iraz, irb, irbm, &
      irbx, irby, irbz, isa, isap, isax, isay, isaz, isb, isbm, isbx, isby, &
      isbz, istat, ix, ixx(1), jj(6), jjm(6), mra, mrb, msa, msax, msb, msbx, &
      n, nra, nrab, nrb, nsa, nsab, nsb, nua, nub
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: abx
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab

    failure = .FALSE.
    CALL timeset(routineN,handle)

    nra = ncoset(lra)
    nsa = ncoset(lsa)
    nrb = ncoset(lrb)
    nsb = ncoset(lsb)
    nrab = ncoset(lra+lrb)
    nsab = ncoset(lsa+lsb)

    CPPrecondition(na*nb<=SIZE(intab,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(nra<=SIZE(intab,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsa<=SIZE(intab,3),cp_failure_level,routineP,error,failure)
    CPPrecondition(nrb<=SIZE(intab,4),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsb<=SIZE(intab,5),cp_failure_level,routineP,error,failure)

    CPPrecondition(na*nb<=SIZE(iab0,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(nrab<=SIZE(iab0,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsab<=SIZE(iab0,3),cp_failure_level,routineP,error,failure)

    n = na*nb
    ALLOCATE (iab(n,nrab,nsab,nrb,nsb),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO mrb=0,lrb
      msbx = lsb
      DO msb=0,msbx
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
                    msax = lsa+lsb-msb
                    DO msa=0,msax
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

                              IF ( n==1 ) THEN
                                iab(1,ira,isa,irb,isb) = iab(1,irap,isap,irbm,isbm) + &
                                                         abx * iab(1,ira,isa,irbm,isbm)
                              ELSE
                                iab(1:n,ira,isa,irb,isb) = iab(1:n,irap,isap,irbm,isbm) + &
                                                           abx * iab(1:n,ira,isa,irbm,isbm)
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
        END IF
      END DO
    END DO

    DO mrb=0,lrb
      nub = ncoset(mrb-1)+1
      nrb = ncoset(mrb)
      msbx = MIN(ltb-mrb,lsb)
      nsb = ncoset(msbx)
      DO mra=0,lra
        nua = ncoset(mra-1)+1
        nra = ncoset(mra)
        msax = MIN(lta-mra,lsa)
        nsa = ncoset(msax)
        intab(1:n,nua:nra,1:nsa,nub:nrb,1:nsb) = intab(1:n,nua:nra,1:nsa,nub:nrb,1:nsb) +&
            iab(1:n,nua:nra,1:nsa,nub:nrb,1:nsb)
      END DO
    END DO

    DEALLOCATE (iab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE hrr
! *****************************************************************************
  SUBROUTINE gg_trans1(in,iout,la,lb,error)
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(IN)                             :: in
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: iout
    INTEGER, INTENT(IN)                      :: la, lb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gg_trans1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i1, i2, i3, istat, &
                                                iu, iv, iw, j1, j2, j3, k1, &
                                                k2, l1, l2, m1, m2, mk, mr, &
                                                ms, nk, nmax, nr, ns
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ac
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a123

    failure = .FALSE.
    CALL timeset(routineN,handle)

    nmax=MAX(la,lb)

    ALLOCATE (ac(0:nmax,0:nmax,0:2*nmax),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL acoff(ac,nmax)

    mk = SIZE(in,1)
    mr = SIZE(in,4)
    ms = SIZE(in,5)
    CPPrecondition(mk<=SIZE(iout,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(mr<=SIZE(iout,4),cp_failure_level,routineP,error,failure)
    CPPrecondition(ms<=SIZE(iout,5),cp_failure_level,routineP,error,failure)

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

    CALL timestop(handle)

  END SUBROUTINE gg_trans1
! *****************************************************************************
  SUBROUTINE gg_trans2(in,iout,la,lb,error)
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(IN)                             :: in
    REAL(dp), DIMENSION(:, :, :, :, :), &
      INTENT(INOUT)                          :: iout
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
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a123

    failure = .FALSE.
    nmax=MAX(la,lb)

    ALLOCATE (ac(0:nmax,0:nmax,0:2*nmax),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL acoff(ac,nmax)

    mk = SIZE(in,1)
    mr = SIZE(in,2)
    ms = SIZE(in,3)

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
