!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Debugs the  integral matrices for the LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created [07.2014] 
!> \authors Dorothea Golze  
! *****************************************************************************
MODULE lri_debug_integrals

  USE ai_overlap,                      ONLY: overlap
  USE ai_overlap3,                     ONLY: overlap3
  USE ai_overlap3_debug,               ONLY: init_os_overlap3,&
                                             os_overlap3
  USE ai_overlap_aabb,                 ONLY: overlap_aabb
  USE ai_overlap_debug,                ONLY: init_os_overlap2,&
                                             os_overlap2
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_debug_integrals'

  PUBLIC :: overlap_ab_test, overlap_abc_test, overlap_aabb_test

! *****************************************************************************

CONTAINS

! ******************************************************************************
!> \brief recursive test routines for integral (a,b) for only two exponents
!> \param error ...
! *****************************************************************************
  SUBROUTINE overlap_ab_test_simple(error)
   
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_ab_test_simple', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia1, iax, iay, iaz, ib1, ibx, &
                                                iby, ibz, la_max, la_min, &
                                                lb_max, lb_min, lds, ma, &
                                                maxl, mb, stat
    INTEGER, DIMENSION(3)                    :: na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dmax, res1, xa, xb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: swork
    REAL(KIND=dp), DIMENSION(1)              :: rpgfa, rpgfb, xa_work, xb_work
    REAL(KIND=dp), DIMENSION(3)              :: A, B, rab

    xa = 0.783300000000   ! exponents
    xb = 1.239648746700
 
    A = (/0.329309000000,0.28408240000,0.28408240000/)!* bohr  !positions
    B = (/0.983983000000,0.00453720000,0.00432740000/)!* bohr
   
    la_min = 0
    lb_min = 0

    la_max = 3
    lb_max = 4
  

    !---------------------------------------
    ALLOCATE(sab(ncoset(la_max),ncoset(lb_max)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    maxl=MAX(la_max, lb_max)
    lds = ncoset(maxl)
    ALLOCATE(swork(lds,lds,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sab=0._dp
    rab(:)= B(:) - A(:)
    dab = SQRT(DOT_PRODUCT(rab,rab))
    xa_work(1)=xa
    xb_work(1)=xb
    rpgfa =20._dp
    rpgfb =20._dp
    CALL overlap(la_max_set=la_max,la_min_set=la_min,npgfa=1,rpgfa=rpgfa,zeta=xa_work,&
               lb_max_set=lb_max,lb_min_set=lb_min,npgfb=1,rpgfb=rpgfb,zetb=xb_work,&
               rab=rab,dab=dab,sab=sab,da_max_set=0,return_derivatives=.FALSE.,s=swork,lds=lds)
    !---------------------------------------

    CALL init_os_overlap2(xa,xb,A,B)

    dmax=0._dp
    DO ma=la_min,la_max
      DO mb=lb_min,lb_max
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
                  res1=os_overlap2(na,nb)
                 ! write(*,*) "la, lb,na, nb, res1", ma, mb, na, nb, res1
                 ! write(*,*) "sab ia1, ib1", ia1, ib1, sab(ia1,ib1)
                  dmax=MAX(dmax,ABS(res1-sab(ia1,ib1)))
                END DO
              END DO
          END DO
        END DO
      END DO
    END DO
   
    WRITE(*,*) "dmax", dmax

   DEALLOCATE(sab,swork,STAT=stat) 
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_ab_test_simple

! ******************************************************************************
!> \brief recursive test routines for integral (a,b)
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param zetb ...
!> \param ra ...
!> \param rb ...
!> \param sab ...
!> \param dmax ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE overlap_ab_test(la_max,la_min,npgfa,zeta,lb_max,lb_min,npgfb,zetb,&
                           ra,rb,sab,dmax,error)
   
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sab
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_ab_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: coa, cob, ia1, iax, iay, iaz, &
                                                ib1, ibx, iby, ibz, ipgf, &
                                                jpgf, ma, mb
    INTEGER, DIMENSION(3)                    :: na, nb
    REAL(KIND=dp)                            :: res1, res2, xa, xb
    REAL(KIND=dp), DIMENSION(3)              :: A, B

    coa=0
    DO ipgf = 1, npgfa
      cob=0
      DO jpgf = 1, npgfb
        xa = zeta(ipgf)   !exponents
        xb = zetb(jpgf)
        A = ra            !positions
        B = rb
        CALL init_os_overlap2(xa,xb,A,B)
        DO ma=la_min,la_max
          DO mb=lb_min,lb_max
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
                    res1=os_overlap2(na,nb)
                    res2=sab(coa+ia1,cob+ib1)
                    dmax=MAX(dmax,ABS(res1-res2))
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        cob = cob + ncoset(lb_max)
      END DO
      coa = coa + ncoset(la_max)
    END DO
    !WRITE(*,*) "dmax overlap_ab_test", dmax
    
  END SUBROUTINE overlap_ab_test

! ******************************************************************************
!> \brief recursive test routines for integral (a,b,c) for only three exponents
!> \param error ...
! *****************************************************************************
  SUBROUTINE overlap_abc_test_simple(error)
   
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_abc_test_simple', &
      routineP = moduleN//':'//routineN

    INTEGER :: ia1, iax, iay, iaz, ib1, ibx, iby, ibz, ic1, icx, icy, icz, &
      la_max, la_min, lb_max, lb_min, lc_max, lc_min, ma, mb, mc, stat
    INTEGER, DIMENSION(3)                    :: na, nb, nc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc, dmax, res1, &
                                                xa, xb, xc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: sabc
    REAL(KIND=dp), DIMENSION(1)              :: rpgfa, rpgfb, rpgfc, xa_work, &
                                                xb_work, xc_work
    REAL(KIND=dp), DIMENSION(3)              :: A, B, C, rab, rac, rbc

    xa = 0.783300000000   ! exponents
    xb = 1.239648746700
    xc = 0.548370000000
 
    A = (/0.329309000000,0.28408240000,0.28408240000/)!* bohr  !positions
    B = (/0.983983000000,0.00453720000,0.00432740000/)!* bohr
    C = (/0.032380000000,1.23470000000,0.11137400000/)!* bohr
   
    la_min = 0
    lb_min = 0
    lc_min = 0

    la_max = 0
    lb_max = 0
    lc_max = 1
  
    !---------------------------------------
    rab(:)= B(:) - A(:)
    dab = SQRT(DOT_PRODUCT(rab,rab))
    rac(:)= C(:) - A(:)
    dac = SQRT(DOT_PRODUCT(rac,rac))
    rbc(:)= C(:) - B(:)
    dbc = SQRT(DOT_PRODUCT(rbc,rbc))
    ALLOCATE(sabc(ncoset(la_max),ncoset(lb_max),ncoset(lc_max)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    xa_work(1)=xa
    xb_work(1)=xb
    xc_work(1)=xc
    rpgfa =20._dp
    rpgfb =20._dp
    rpgfc =20._dp
    sabc=0._dp

    CALL overlap3(la_max_set=la_max,npgfa=1,zeta=xa_work,rpgfa=rpgfa,la_min_set=la_min,&
                  lb_max_set=lb_max,npgfb=1,zetb=xb_work,rpgfb=rpgfb,lb_min_set=lb_min,&
                  lc_max_set=lc_max,npgfc=1,zetc=xc_work,rpgfc=rpgfc,lc_min_set=lc_min,&
                  rab=rab,dab=dab,rac=rac,dac=dac,rbc=rbc,dbc=dbc,sabc=sabc,error=error)
    !sabc(:,:,:)=0._dp 
   !CALL overlap3_old(la_max=la_max,npgfa=1,zeta=xa_work,rpgfa=rpgfa,la_min=la_min,&
   !              lb_max=lb_max,npgfb=1,zetb=xb_work,rpgfb=rpgfb,lb_min=lb_min,&
   !              lc_max=lc_max,npgfc=1,zetc=xc,rpgfc=20._dp,lc_min=lc_min,&
   !              rab=rab,dab=dab,rac=rac,dac=dac,rbc=rbc,dbc=dbc,sabc=sabc)
  
  !---------------------------------------
 
    CALL init_os_overlap3(xa,xb,xc,A,B,C)

    dmax=0._dp
    DO ma=la_min,la_max
      DO mc=lc_min,lc_max
        DO mb=lb_min,lb_max
          DO iax=0,ma
            DO iay=0,ma-iax
              iaz=ma-iax-iay
              na(1)=iax; na(2)=iay; na(3)=iaz
              ia1=coset(iax,iay,iaz)
              DO icx=0,mc
                DO icy=0,mc-icx
                  icz=mc-icx-icy
                  nc(1)=icx; nc(2)=icy; nc(3)=icz
                  ic1=coset(icx,icy,icz)
                  DO ibx=0,mb
                    DO iby=0,mb-ibx
                      ibz=mb-ibx-iby
                      nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                      ib1=coset(ibx,iby,ibz)
                      res1=os_overlap3(na,nc,nb)
                      !write(*,*) "la, lc, lb,na,nc, nb, res1", ma, mc, mb, na, nc, nb, res1
                     !write(*,*) "sabc ia1, ib1, ic1", ia1, ib1, ic1, sabc(ia1,ib1,ic1)
                      dmax=MAX(dmax,ABS(res1-sabc(ia1,ib1,ic1)))
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
   
    WRITE(*,*) "dmax", dmax

   DEALLOCATE(sabc,STAT=stat) 
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
  END SUBROUTINE overlap_abc_test_simple

! ******************************************************************************
!> \brief recursive test routines for integral (a,b,c)
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param lb_min ...
!> \param lc_max ...
!> \param npgfc ...
!> \param zetc ...
!> \param lc_min ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param sabc ...
!> \param dmax ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE overlap_abc_test(la_max,npgfa,zeta,la_min,&
                           lb_max,npgfb,zetb,lb_min,&
                           lc_max,npgfc,zetc,lc_min,&
                           ra,rb,rc,sabc,dmax,error)
   
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    INTEGER, INTENT(IN)                      :: lb_min, lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: sabc
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_abc_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: coa, cob, coc, ia1, iax, iay, &
                                                iaz, ib1, ibx, iby, ibz, ic1, &
                                                icx, icy, icz, ipgf, jpgf, &
                                                kpgf, ma, mb, mc
    INTEGER, DIMENSION(3)                    :: na, nb, nc
    REAL(KIND=dp)                            :: res1, res2, xa, xb, xc
    REAL(KIND=dp), DIMENSION(3)              :: A, B, C

    coa=0
    DO ipgf = 1, npgfa
      cob=0
      DO jpgf = 1, npgfb
        coc=0
        DO kpgf = 1, npgfc
 
          xa = zeta(ipgf)   ! exponents
          xb = zetb(jpgf)
          xc = zetc(kpgf)
 
          A = Ra  !positions
          B = Rb
          C = Rc
          
          CALL init_os_overlap3(xa,xb,xc,A,B,C)

          DO ma=la_min,la_max
            DO mc=lc_min,lc_max
              DO mb=lb_min,lb_max
                DO iax=0,ma
                  DO iay=0,ma-iax
                    iaz=ma-iax-iay
                    na(1)=iax; na(2)=iay; na(3)=iaz
                    ia1=coset(iax,iay,iaz)    
                    DO icx=0,mc
                      DO icy=0,mc-icx
                        icz=mc-icx-icy
                        nc(1)=icx; nc(2)=icy; nc(3)=icz
                        ic1=coset(icx,icy,icz)
                        DO ibx=0,mb
                          DO iby=0,mb-ibx
                            ibz=mb-ibx-iby
                            nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                            ib1=coset(ibx,iby,ibz)
                            res1=os_overlap3(na,nc,nb)
                            res2=sabc(coa+ia1,cob+ib1,coc+ic1)
                            dmax=MAX(dmax,ABS(res1-res2))
                            !IF(dmax > 1.E-10) WRITE(*,*) "dmax in loop", dmax
                          END DO
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
          coc = coc + ncoset(lc_max)
        END DO 
        cob = cob + ncoset(lb_max)
      END DO
      coa = coa + ncoset(la_max)
    END DO
    !WRITE(*,*) "dmax abc", dmax
    
  END SUBROUTINE overlap_abc_test

! ******************************************************************************
!> \brief recursive test routines for integral (aa,bb) for only four exponents
!> \param error ...
! *****************************************************************************
  SUBROUTINE overlap_aabb_test_simple(error)
   
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_aabb_test_simple', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, iax, iay, iaz, ibx, iby, ibz, j, k, l, la_max, la_max1, &
      la_max2, la_min, la_min1, la_min2, lb_max, lb_max1, lb_max2, lb_min, &
      lb_min1, lb_min2, lds, ma, maxl, mb, stat
    INTEGER, DIMENSION(3)                    :: na, naa, nb, nbb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dmax, res1, xa, xa1, &
                                                xa2, xb, xb1, xb2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: swork
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: saabb
    REAL(KIND=dp), DIMENSION(1)              :: rpgfa1, rpgfa2, rpgfb1, &
                                                rpgfb2, xa_work1, xa_work2, &
                                                xb_work1, xb_work2
    REAL(KIND=dp), DIMENSION(3)              :: A, B, rab

    xa1 = 0.783300000000   ! exponents
    xb1 = 1.239648746700
    xa2 = 0.3400000000   ! exponents
    xb2 = 2.76
 
    A = (/0.329309000000,0.28408240000,0.28408240000/)!* bohr  !positions
    B = (/0.983983000000,0.00453720000,0.00432740000/)!* bohr
   
    la_min1 = 0
    la_min2 = 0
    lb_min1 = 3
    lb_min2 = 1

    la_max1 = 1
    la_max2 = 2
    lb_max1 = 3
    lb_max2 = 4
  

    !---------------------------------------
    ALLOCATE(saabb(ncoset(la_max1),ncoset(la_max2),ncoset(lb_max1),ncoset(lb_max2)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    maxl=MAX(la_max1+la_max2, lb_max1+lb_max2)
    lds = ncoset(maxl)
    ALLOCATE(swork(lds,lds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    saabb=0._dp
    rab(:)= B(:) - A(:)
    dab = SQRT(DOT_PRODUCT(rab,rab))
    xa_work1(1)=xa1
    xa_work2(1)=xa2
    xb_work1(1)=xb1
    xb_work2(1)=xb2
    rpgfa1 =20._dp
    rpgfa2 =20._dp
    rpgfb1 =20._dp
    rpgfb2 =20._dp
    CALL overlap_aabb(la_max_set1=la_max1,la_min_set1=la_min1,npgfa1=1,rpgfa1=rpgfa1,zeta1=xa_work1,&
                      la_max_set2=la_max2,la_min_set2=la_min2,npgfa2=1,rpgfa2=rpgfa2,zeta2=xa_work2,&
                      lb_max_set1=lb_max1,lb_min_set1=lb_min1,npgfb1=1,rpgfb1=rpgfb1,zetb1=xb_work1,&
                      lb_max_set2=lb_max2,lb_min_set2=lb_min2,npgfb2=1,rpgfb2=rpgfb2,zetb2=xb_work2,&
                      asets_equal=.FALSE.,bsets_equal=.FALSE.,rab=rab,dab=dab,saabb=saabb,s=swork,lds=lds)
    !---------------------------------------
     
    xa = xa1 + xa2
    xb = xb1 + xb2
    la_min = la_min1 + la_min2
    la_max = la_max1 + la_max2
    lb_min = lb_min1 + lb_min2
    lb_max = lb_max1 + lb_max2

    CALL init_os_overlap2(xa,xb,A,B)

    dmax=0._dp
    DO ma=la_min,la_max
      DO mb=lb_min,lb_max
        DO iax=0,ma
          DO iay=0,ma-iax
            iaz=ma-iax-iay
            na(1)=iax; na(2)=iay; na(3)=iaz
            DO ibx=0,mb
              DO iby=0,mb-ibx
                ibz=mb-ibx-iby
                nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                res1=os_overlap2(na,nb)
                DO i=ncoset(la_min1-1)+1,ncoset(la_max1)
                  DO j=ncoset(la_min2-1)+1,ncoset(la_max2) 
                    naa = indco(1:3,i) + indco(1:3,j)
                    DO k=ncoset(lb_min1-1)+1,ncoset(lb_max1)
                      DO l=ncoset(lb_min2-1)+1,ncoset(lb_max2)
                        nbb = indco(1:3,k) + indco(1:3,l)
                        IF(ALL(na == naa).AND.ALL(nb == nbb)) THEN
                          dmax=MAX(dmax,ABS(res1-saabb(i,j,k,l)))
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
   
    WRITE(*,*) "dmax", dmax

   DEALLOCATE(saabb,swork,STAT=stat) 
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_aabb_test_simple

! ******************************************************************************
!> \brief recursive test routines for integral (aa,bb)
!> \param la_max1 ...
!> \param la_min1 ...
!> \param npgfa1 ...
!> \param zeta1 ...
!> \param la_max2 ...
!> \param la_min2 ...
!> \param npgfa2 ...
!> \param zeta2 ...
!> \param lb_max1 ...
!> \param lb_min1 ...
!> \param npgfb1 ...
!> \param zetb1 ...
!> \param lb_max2 ...
!> \param lb_min2 ...
!> \param npgfb2 ...
!> \param zetb2 ...
!> \param ra ...
!> \param rb ...
!> \param saabb ...
!> \param dmax ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE overlap_aabb_test(la_max1,la_min1,npgfa1,zeta1,&
                                 la_max2,la_min2,npgfa2,zeta2,&
                                 lb_max1,lb_min1,npgfb1,zetb1,&
                                 lb_max2,lb_min2,npgfb2,zetb2,&
                                 ra,rb,saabb,dmax,error)
   
    INTEGER, INTENT(IN)                      :: la_max1, la_min1, npgfa1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta1
    INTEGER, INTENT(IN)                      :: la_max2, la_min2, npgfa2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta2
    INTEGER, INTENT(IN)                      :: lb_max1, lb_min1, npgfb1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb1
    INTEGER, INTENT(IN)                      :: lb_max2, lb_min2, npgfb2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: saabb
    REAL(KIND=dp), INTENT(INOUT)             :: dmax
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_aabb_test', &
      routineP = moduleN//':'//routineN

    INTEGER :: coa1, coa2, cob1, cob2, i, iax, iay, iaz, ibx, iby, ibz, ipgf, &
      j, jpgf, k, kpgf, l, la_max, la_min, lb_max, lb_min, lpgf, ma, mb
    INTEGER, DIMENSION(3)                    :: na, naa, nb, nbb
    REAL(KIND=dp)                            :: res1, xa, xb
    REAL(KIND=dp), DIMENSION(3)              :: A, B

    coa1=0
    DO ipgf = 1, npgfa1
      coa2=0
      DO jpgf = 1, npgfa2
        cob1=0
        DO kpgf = 1, npgfb1
          cob2=0
          DO lpgf = 1, npgfb2

            xa = zeta1(ipgf) + zeta2(jpgf)   ! exponents
            xb = zetb1(kpgf) + zetb2(lpgf)   ! exponents
            la_max = la_max1 + la_max2            
            lb_max = lb_max1 + lb_max2            
            la_min = la_min1 + la_min2            
            lb_min = lb_min1 + lb_min2            

            A = ra  !positions
            B = rb

            CALL init_os_overlap2(xa,xb,A,B)
 
            DO ma=la_min,la_max
              DO mb=lb_min,lb_max
                DO iax=0,ma
                  DO iay=0,ma-iax
                    iaz=ma-iax-iay
                    na(1)=iax; na(2)=iay; na(3)=iaz
                    DO ibx=0,mb
                      DO iby=0,mb-ibx
                        ibz=mb-ibx-iby
                        nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                        res1=os_overlap2(na,nb)
                        DO i=ncoset(la_min1-1)+1,ncoset(la_max1)
                          DO j=ncoset(la_min2-1)+1,ncoset(la_max2) 
                            naa = indco(1:3,i) + indco(1:3,j)
                            DO k=ncoset(lb_min1-1)+1,ncoset(lb_max1)
                              DO l=ncoset(lb_min2-1)+1,ncoset(lb_max2)
                                nbb = indco(1:3,k) + indco(1:3,l)
                                IF(ALL(na == naa).AND.ALL(nb == nbb)) THEN
                                  dmax=MAX(dmax,ABS(res1-saabb(coa1+i,coa2+j,cob1+k,cob2+l)))
                                ENDIF
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
            cob2 = cob2 + ncoset(lb_max2)
          END DO
          cob1 = cob1 + ncoset(lb_max1)
        END DO
        coa2 = coa2 + ncoset(la_max2)
      END DO
      coa1 = coa1 + ncoset(la_max1)
    END DO
   
    !WRITE(*,*) "dmax aabb", dmax

  END SUBROUTINE overlap_aabb_test

END MODULE lri_debug_integrals
