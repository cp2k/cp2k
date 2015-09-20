!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of general three-center integrals over Cartesian
!>      Gaussian-type functions and a spherical operator centered at position C
!>
!>      <a|V(local)|b> = <a|F(|r-C|)|b>
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Based in part on code by MK
!> \par Parameters
!>      -  ax,ay,az   : Angular momentum index numbers of orbital a.
!>      -  bx,by,bz   : Angular momentum index numbers of orbital b.
!>      -  coset      : Cartesian orbital set pointer.
!>      -  dab        : Distance between the atomic centers a and b.
!>      -  dac        : Distance between the atomic centers a and c.
!>      -  dbc        : Distance between the atomic centers b and c.
!>      -  l{a,b}     : Angular momentum quantum number of shell a or b.
!>      -  l{a,b}_max : Maximum angular momentum quantum number of shell a or b.
!>      -  ncoset     : Number of Cartesian orbitals up to l.
!>      -  rab        : Distance vector between the atomic centers a and b.
!>      -  rac        : Distance vector between the atomic centers a and c.
!>      -  rbc        : Distance vector between the atomic centers b and c.
!>      -  rpgf{a,b,c}: Radius of the primitive Gaussian-type function a or b.
!>      -  zet{a,b}   : Exponents of the Gaussian-type functions a or b.
!> \author jhu (05.2011)
! *****************************************************************************
MODULE ai_oneelectron
  
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_oneelectron'

! *** Public subroutines ***

  PUBLIC :: os_3center

CONTAINS

! *****************************************************************************
!> \brief   Calculation of three-center integrals <a|c|b> over
!>           Cartesian Gaussian functions and a spherical potential
!>
!> \param la_max_set ...
!> \param la_min_set ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max_set ...
!> \param lb_min_set ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param auxint ...
!> \param rpgfc ...
!> \param rab ...
!> \param dab ...
!> \param rac ...
!> \param dac ...
!> \param rbc ...
!> \param dbc ...
!> \param vab ...
!> \param s ...
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
!> \param fs ...
!> \date    May 2011
!> \author  Juerg Hutter
!> \version 1.0
! *****************************************************************************
  SUBROUTINE os_3center(la_max_set,la_min_set,npgfa,rpgfa,zeta,&
       lb_max_set,lb_min_set,npgfb,rpgfb,zetb,auxint,rpgfc,&
       rab,dab,rac,dac,rbc,dbc,vab,s,pab,force_a,force_b,fs)
    INTEGER, INTENT(IN)                      :: la_max_set, la_min_set, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max_set, lb_min_set, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(0:, :), &
      INTENT(IN)                             :: auxint
    REAL(KIND=dp), INTENT(IN)                :: rpgfc
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: dac
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rbc
    REAL(KIND=dp), INTENT(IN)                :: dbc
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pab
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: force_a, force_b
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: fs

    CHARACTER(len=*), PARAMETER :: routineN = 'os_3center', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, cdb, cdbx, &
      cdby, cdbz, coa, coamx, coamy, coamz, coapx, coapy, coapz, cob, cobmx, &
      cobmy, cobmz, cobpx, cobpy, cobpz, da, da_max, dax, day, daz, db, &
      db_max, dbx, dby, dbz, i, ia, iap, iax, iay, iaz, ib, ibm, ibx, iby, &
      ibz, ii(3), iim(3), ij, ipgf, ir, ir1, ir2, irm(3), irr(3), irx, iry, &
      irz, istat, ix, ixx(1), j, jj(3), jjp(3), jpgf, la, la_max, la_min, lb, &
      lb_max, lb_min, llr, m, ma, mb, mmax, na, nb
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: iiap
    LOGICAL                                  :: calculate_force_a, &
                                                calculate_force_b, failure
    REAL(KIND=dp)                            :: aai, abx, fax, fay, faz, fbx, &
                                                fby, fbz, ftz, orho, rho, s1, &
                                                s2
    REAL(KIND=dp), DIMENSION(3)              :: pai, pbi, pci

    failure = .FALSE.
    IF (PRESENT(pab)) THEN
      CPASSERT(PRESENT(fs))
      IF (PRESENT(force_a)) THEN
        calculate_force_a = .TRUE.
      ELSE
        calculate_force_a = .FALSE.
      END IF
      IF (PRESENT(force_b)) THEN
        calculate_force_b = .TRUE.
      ELSE
        calculate_force_b = .FALSE.
      END IF
    ELSE
      calculate_force_a = .FALSE.
      calculate_force_b = .FALSE.
    END IF

    IF( calculate_force_a ) THEN
       da_max = 1
       force_a = 0.0_dp
    ELSE
       da_max = 0
    END IF

    IF( calculate_force_b ) THEN
       db_max = 1
       force_b = 0.0_dp
    ELSE
       db_max = 0
    END IF

    la_max = la_max_set + da_max
    la_min = MAX(0,la_min_set-da_max)

    lb_max = lb_max_set + db_max
    lb_min = MAX(0,lb_min_set-db_max)

    mmax = la_max + lb_max

    ! precalculate indices for horizontal recursion
    ALLOCATE (iiap(ncoset(mmax),3),STAT=istat)
    CPASSERT(istat==0)
    DO ma=0,mmax
       DO iax=0,ma
          DO iay=0,ma-iax
             iaz=ma-iax-iay
             ia=coset(iax,iay,iaz)
             jj(1)=iax; jj(2)=iay; jj(3)=iaz
             jjp = jj
             jjp(1) = jjp(1) + 1
             iap=coset(jjp(1),jjp(2),jjp(3))
             iiap(ia,1) = iap
             jjp = jj
             jjp(2) = jjp(2) + 1
             iap=coset(jjp(1),jjp(2),jjp(3))
             iiap(ia,2) = iap
             jjp = jj
             jjp(3) = jjp(3) + 1
             iap=coset(jjp(1),jjp(2),jjp(3))
             iiap(ia,3) = iap
          END DO
       END DO
    END DO

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

!     *** Screening ***
      IF (rpgfa(ipgf) + rpgfc < dac) THEN
        na = na + ncoset(la_max_set)
        CYCLE
      END IF

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***
        IF ((rpgfb(jpgf) + rpgfc < dbc).OR.&
           (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
           nb = nb + ncoset(lb_max_set)
           CYCLE
        END IF

!       *** Calculate some prefactors ***
        rho = zeta(ipgf) + zetb(jpgf)
        pai(:) = zetb(jpgf)/rho*rab(:)
        pbi(:) = -zeta(ipgf)/rho*rab(:)
        pci(:) = -(zeta(ipgf)*rac(:) + zetb(jpgf)*rbc(:))/rho
        orho = 0.5_dp/rho

        ij = (ipgf-1)*npgfb + jpgf
        s(1,1,1:mmax+1) = auxint(0:mmax,ij)

        IF (la_max > 0) THEN
!         *** Recurrence steps: [s|c|s] -> [a|c|s]               ***
!         *** [a|c|s](m) = (Pi - Ai)*[a-1i|c|s](m) -             ***
!         ***              (Pi - Ci)*[a-1i|c|s](m+1)) +          ***
!         ***              Ni(a-1i)/2(a+b)*[a-2i|c|s](m) -       ***
!         ***              Ni(a-1i)/2(a+b)*[a-2i|c|s](m+1)       ***
          DO llr = 1,mmax
            IF (llr == 1)  THEN
               DO m=0,mmax-llr
                  s1 = s(1,1,m+1)
                  s2 = s(1,1,m+2)
                  s(2,1,m+1) = pai(1)*s1 - pci(1)*s2 ! [px|o|s]
                  s(3,1,m+1) = pai(2)*s1 - pci(2)*s2 ! [py|o|s]
                  s(4,1,m+1) = pai(3)*s1 - pci(3)*s2 ! [pz|o|s]
               END DO
            ELSE IF( llr == 2) THEN
               DO m=0,mmax-llr
                  s1 = s(1,1,m+1)-s(1,1,m+2)
                  s( 5,1,m+1) = pai(1)*s(2,1,m+1) - pci(1)*s(2,1,m+2) + orho*s1 ! [dx2|o|s]
                  s( 6,1,m+1) = pai(1)*s(3,1,m+1) - pci(1)*s(3,1,m+2)           ! [dxy|o|s]
                  s( 7,1,m+1) = pai(1)*s(4,1,m+1) - pci(1)*s(4,1,m+2)           ! [dxz|o|s]
                  s( 8,1,m+1) = pai(2)*s(3,1,m+1) - pci(2)*s(3,1,m+2) + orho*s1 ! [dy2|o|s]
                  s( 9,1,m+1) = pai(2)*s(4,1,m+1) - pci(2)*s(4,1,m+2)           ! [dyz|o|s]
                  s(10,1,m+1) = pai(3)*s(4,1,m+1) - pci(3)*s(4,1,m+2) + orho*s1 ! [dz2|o|s]
               END DO
            ELSE IF( llr == 3) THEN
               DO m=0,mmax-llr
                  s(11,1,m+1) = pai(1)*s( 5,1,m+1) - pci(1)*s( 5,1,m+2)  &      ! [fx3 |o|s]
                                + 2._dp*orho*(s(2,1,m+1)-s(2,1,m+2))
                  s(12,1,m+1) = pai(1)*s( 6,1,m+1) - pci(1)*s( 6,1,m+2)  &      ! [fx2y|o|s]
                                + orho*(s(3,1,m+1)-s(3,1,m+2))
                  s(13,1,m+1) = pai(1)*s( 7,1,m+1) - pci(1)*s( 7,1,m+2)  &      ! [fx2z|o|s]
                                + orho*(s(4,1,m+1)-s(4,1,m+2))
                  s(14,1,m+1) = pai(2)*s( 6,1,m+1) - pci(2)*s( 6,1,m+2)  &      ! [fxy2|o|s]
                                + orho*(s(2,1,m+1)-s(2,1,m+2))
                  s(15,1,m+1) = pai(1)*s( 9,1,m+1) - pci(1)*s( 9,1,m+2)         ! [fxyz|o|s]
                  s(16,1,m+1) = pai(3)*s( 7,1,m+1) - pci(3)*s( 7,1,m+2)  &      ! [fxz2|o|s]
                                + orho*(s(2,1,m+1)-s(2,1,m+2))
                  s(17,1,m+1) = pai(2)*s( 8,1,m+1) - pci(2)*s( 8,1,m+2)  &      ! [fy3 |o|s]
                                + 2._dp*orho*(s(3,1,m+1)-s(3,1,m+2))
                  s(18,1,m+1) = pai(2)*s( 9,1,m+1) - pci(2)*s( 9,1,m+2)  &      ! [fy2z|o|s]
                                + orho*(s(4,1,m+1)-s(4,1,m+2))
                  s(19,1,m+1) = pai(3)*s( 9,1,m+1) - pci(3)*s( 9,1,m+2)  &      ! [fyz2|o|s]
                                + orho*(s(3,1,m+1)-s(3,1,m+2))
                  s(20,1,m+1) = pai(3)*s(10,1,m+1) - pci(3)*s(10,1,m+2)  &      ! [fz3 |o|s]
                                + 2._dp*orho*(s(4,1,m+1)-s(4,1,m+2))
               END DO
            ELSE IF( llr == 4) THEN
               DO m=0,mmax-llr
                  s(21,1,m+1) = pai(1)*s(11,1,m+1) - pci(1)*s(11,1,m+2)  & ! [gx4  |s|s]
                                + 3._dp*orho*(s(5,1,m+1)-s(5,1,m+2))
                  s(22,1,m+1) = pai(1)*s(12,1,m+1) - pci(1)*s(12,1,m+2)  & ! [gx3y |s|s]
                                + 2._dp*orho*(s(6,1,m+1)-s(6,1,m+2))
                  s(23,1,m+1) = pai(1)*s(13,1,m+1) - pci(1)*s(13,1,m+2)  & ! [gx3z |s|s]
                                + 2._dp*orho*(s(7,1,m+1)-s(7,1,m+2))
                  s(24,1,m+1) = pai(1)*s(14,1,m+1) - pci(1)*s(14,1,m+2)  & ! [gx2y2|s|s]
                                + orho*(s(8,1,m+1)-s(8,1,m+2))
                  s(25,1,m+1) = pai(1)*s(15,1,m+1) - pci(1)*s(15,1,m+2)  & ! [gx2yz|s|s]
                                + orho*(s(9,1,m+1)-s(9,1,m+2))
                  s(26,1,m+1) = pai(1)*s(16,1,m+1) - pci(1)*s(16,1,m+2)  & ! [gx2z2|s|s]
                                + orho*(s(10,1,m+1)-s(10,1,m+2))
                  s(27,1,m+1) = pai(1)*s(17,1,m+1) - pci(1)*s(17,1,m+2)    ! [gxy3 |s|s]
                  s(28,1,m+1) = pai(1)*s(18,1,m+1) - pci(1)*s(18,1,m+2)    ! [gxy2z|s|s]
                  s(29,1,m+1) = pai(1)*s(19,1,m+1) - pci(1)*s(19,1,m+2)    ! [gxyz2|s|s]
                  s(30,1,m+1) = pai(1)*s(20,1,m+1) - pci(1)*s(20,1,m+2)    ! [gxz3 |s|s]
                  s(31,1,m+1) = pai(2)*s(17,1,m+1) - pci(2)*s(17,1,m+2)  & ! [gy4  |s|s]
                                + 3._dp*orho*(s(8,1,m+1)-s(8,1,m+2))
                  s(32,1,m+1) = pai(2)*s(18,1,m+1) - pci(2)*s(18,1,m+2)  & ! [gy3z |s|s]
                                + 2._dp*orho*(s(9,1,m+1)-s(9,1,m+2))
                  s(33,1,m+1) = pai(2)*s(19,1,m+1) - pci(2)*s(19,1,m+2)  & ! [gy2z2|s|s]
                                + orho*(s(10,1,m+1)-s(10,1,m+2))
                  s(34,1,m+1) = pai(2)*s(20,1,m+1) - pci(2)*s(20,1,m+2)    ! [gyz3 |s|s]
                  s(35,1,m+1) = pai(3)*s(20,1,m+1) - pci(3)*s(20,1,m+2)  & ! [gz4  |s|s]
                                + 3._dp*orho*(s(10,1,m+1)-s(10,1,m+2))
               END DO
            ELSE
               DO irx=0,llr
                 DO iry=0,llr-irx
                   irz=llr-irx-iry
                   irr(1)=irx; irr(2)=iry; irr(3)=irz
                   ixx = MAXLOC(irr)
                   ix = ixx(1)
                   ir=coset(irx,iry,irz)
                   irm = irr
                   irm(ix) = irm(ix) - 1
                   aai = REAL(MAX(irm(ix),0),dp)*orho
                   ir1=coset(irm(1),irm(2),irm(3))
                   irm(ix) = irm(ix) - 1
                   ir2=coset(irm(1),irm(2),irm(3))
                   DO m=0,mmax-llr
                     s(ir,1,m+1) = pai(ix)*s(ir1,1,m+1) - pci(ix)*s(ir1,1,m+2) &
                                   + aai*(s(ir2,1,m+1)-s(ir2,1,m+2))
                   END DO
                 END DO
               END DO
            END IF
          END DO

!         *** Horizontal recurrence steps ***
!         *** [a|c|b+1i] = [a+1i|c|b] + (Ai - Bi)*[a|c|b] ***

          DO mb=1,lb_max
            DO ibx=0,mb
              DO iby=0,mb-ibx
                ibz=mb-ibx-iby
                ib=coset(ibx,iby,ibz)
                ii(1)=ibx; ii(2)=iby; ii(3)=ibz
                ixx = MAXLOC(ii)
                ix = ixx(1)
                abx = -rab(ix)
                iim = ii
                iim(ix) = iim(ix) - 1
                ibm=coset(iim(1),iim(2),iim(3))
                DO ia=1,ncoset(mmax-mb)
                  iap = iiap(ia,ix)
                  s(ia,ib,1) = s(iap,ibm,1) + abx * s(ia,ibm,1)
                END DO
              END DO
            END DO
          END DO

        ELSE IF (lb_max > 0) THEN

!         *** Recurrence steps: [s|c|s] -> [s|c|b]               ***
!         *** [s|c|b](m) = (Pi - Bi)*[s|c|b-1i](m) -             ***
!         ***              (Pi - Ci)*[s|c|b-1i](m+1)) +          ***
!         ***              Ni(b-1i)/2(a+b)*[s|c|b-2i](m) -       ***
!         ***              Ni(b-1i)/2(a+b)*[s|c|b-2i](m+1)       ***
          DO llr = 1,lb_max
            IF (llr == 1)  THEN
               DO m=0,lb_max-llr
                  s1 = s(1,1,m+1)
                  s2 = s(1,1,m+2)
                  s(1,2,m+1) = pbi(1)*s1 - pci(1)*s2 ! [px|o|s]
                  s(1,3,m+1) = pbi(2)*s1 - pci(2)*s2 ! [py|o|s]
                  s(1,4,m+1) = pbi(3)*s1 - pci(3)*s2 ! [pz|o|s]
               END DO
            ELSE IF( llr == 2) THEN
               DO m=0,lb_max-llr
                  s1 = s(1,1,m+1)-s(1,1,m+2)
                  s(1, 5,m+1) = pbi(1)*s(1,2,m+1) - pci(1)*s(1,2,m+2) + orho*s1 ! [dx2|o|s]
                  s(1, 6,m+1) = pbi(1)*s(1,3,m+1) - pci(1)*s(1,3,m+2)           ! [dxy|o|s]
                  s(1, 7,m+1) = pbi(1)*s(1,4,m+1) - pci(1)*s(1,4,m+2)           ! [dxz|o|s]
                  s(1, 8,m+1) = pbi(2)*s(1,3,m+1) - pci(2)*s(1,3,m+2) + orho*s1 ! [dy2|o|s]
                  s(1, 9,m+1) = pbi(2)*s(1,4,m+1) - pci(2)*s(1,4,m+2)           ! [dyz|o|s]
                  s(1,10,m+1) = pbi(3)*s(1,4,m+1) - pci(3)*s(1,4,m+2) + orho*s1 ! [dz2|o|s]
               END DO
            ELSE IF( llr == 3) THEN
               DO m=0,lb_max-llr
                  s(1,11,m+1) = pbi(1)*s(1, 5,m+1) - pci(1)*s(1, 5,m+2)  &      ! [fx3 |o|s]
                                + 2._dp*orho*(s(1,2,m+1)-s(1,2,m+2))
                  s(1,12,m+1) = pbi(1)*s(1, 6,m+1) - pci(1)*s(1, 6,m+2)  &      ! [fx2y|o|s]
                                + orho*(s(1,3,m+1)-s(1,3,m+2))
                  s(1,13,m+1) = pbi(1)*s(1, 7,m+1) - pci(1)*s(1, 7,m+2)  &      ! [fx2z|o|s]
                                + orho*(s(1,4,m+1)-s(1,4,m+2))
                  s(1,14,m+1) = pbi(2)*s(1, 6,m+1) - pci(2)*s(1, 6,m+2)  &      ! [fxy2|o|s]
                                + orho*(s(1,2,m+1)-s(1,2,m+2))
                  s(1,15,m+1) = pbi(1)*s(1, 9,m+1) - pci(1)*s(1, 9,m+2)         ! [fxyz|o|s]
                  s(1,16,m+1) = pbi(3)*s(1, 7,m+1) - pci(3)*s(1, 7,m+2)  &      ! [fxz2|o|s]
                                + orho*(s(1,2,m+1)-s(1,2,m+2))
                  s(1,17,m+1) = pbi(2)*s(1, 8,m+1) - pci(2)*s(1, 8,m+2)  &      ! [fy3 |o|s]
                                + 2._dp*orho*(s(1,3,m+1)-s(1,3,m+2))
                  s(1,18,m+1) = pbi(2)*s(1, 9,m+1) - pci(2)*s(1, 9,m+2)  &      ! [fy2z|o|s]
                                + orho*(s(1,4,m+1)-s(1,4,m+2))
                  s(1,19,m+1) = pbi(3)*s(1, 9,m+1) - pci(3)*s(1, 9,m+2)  &      ! [fyz2|o|s]
                                + orho*(s(1,3,m+1)-s(1,3,m+2))
                  s(1,20,m+1) = pbi(3)*s(1,10,m+1) - pci(3)*s(1,10,m+2)  &      ! [fz3 |o|s]
                                + 2._dp*orho*(s(1,4,m+1)-s(1,4,m+2))
               END DO
            ELSE
               DO irx=0,llr
                 DO iry=0,llr-irx
                   irz=llr-irx-iry
                   irr(1)=irx; irr(2)=iry; irr(3)=irz
                   ixx = MAXLOC(irr)
                   ix = ixx(1)
                   ir=coset(irx,iry,irz)
                   irm = irr
                   irm(ix) = irm(ix) - 1
                   aai = REAL(MAX(irm(ix),0),dp)
                   ir1=coset(irm(1),irm(2),irm(3))
                   irm(ix) = irm(ix) - 1
                   ir2=coset(irm(1),irm(2),irm(3))
                   DO m=0,lb_max-llr
                     s(1,ir,m+1) = pbi(ix)*s(1,ir1,m+1) - pci(ix)*s(1,ir1,m+2) &
                                   + aai*orho*(s(1,ir2,m+1)-s(1,ir2,m+2))
                   END DO
                 END DO
               END DO
            END IF
          END DO

        END IF

!       *** Store the primitive three-center overlap integrals ***
        DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
          DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
            vab(na+i,nb+j) = vab(na+i,nb+j) + s(i,j,1)
          END DO
        END DO

!       *** Calculate the requested derivatives with respect  ***
!       *** to the nuclear coordinates of the atomic center a ***

        DO da=0,da_max-1
          ftz = 2.0_dp*zeta(ipgf)
          DO dax=0,da
            DO day=0,da-dax
              daz = da - dax - day
              cda = coset(dax,day,daz)
              cdax = coset(dax+1,day,daz)
              cday = coset(dax,day+1,daz)
              cdaz = coset(dax,day,daz+1)

!             *** [da/dAi|c|b] = 2*zeta*[a+1i|c|b] - Ni(a)[a-1i|c|b] ***

              DO la=la_min_set,la_max-da-1
                DO ax=0,la
                  fax = REAL(ax,dp)
                  DO ay=0,la-ax
                    fay = REAL(ay,dp)
                    az = la - ax - ay
                    faz = REAL(az,dp)
                    coa = coset(ax,ay,az)
                    coamx = coset(ax-1,ay,az)
                    coamy = coset(ax,ay-1,az)
                    coamz = coset(ax,ay,az-1)
                    coapx = coset(ax+1,ay,az)
                    coapy = coset(ax,ay+1,az)
                    coapz = coset(ax,ay,az+1)
                    DO cob=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
                       fs(coa,cob,cdax) = ftz*s(coapx,cob,cda) - fax*s(coamx,cob,cda)
                       fs(coa,cob,cday) = ftz*s(coapy,cob,cda) - fay*s(coamy,cob,cda)
                       fs(coa,cob,cdaz) = ftz*s(coapz,cob,cda) - faz*s(coamz,cob,cda)
                    END DO
                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO

!       *** Calculate the force contribution for the atomic center a ***

        IF (calculate_force_a) THEN
            DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
              DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                force_a(1) = force_a(1) + pab(na+i,nb+j)*fs(i,j,2)
                force_a(2) = force_a(2) + pab(na+i,nb+j)*fs(i,j,3)
                force_a(3) = force_a(3) + pab(na+i,nb+j)*fs(i,j,4)
              END DO
            END DO
        END IF

!       *** Calculate the requested derivatives with respect  ***
!       *** to the nuclear coordinates of the atomic center b ***

        DO db=0,db_max-1
          ftz = 2.0_dp*zetb(jpgf)
          DO dbx=0,db
            DO dby=0,db-dbx
              dbz = db - dbx - dby
              cdb = coset(dbx,dby,dbz)
              cdbx = coset(dbx+1,dby,dbz)
              cdby = coset(dbx,dby+1,dbz)
              cdbz = coset(dbx,dby,dbz+1)

!             *** [a|c|db/dBi] = 2*zetb*[a|c|b+1i] - Ni(b)[a|c|b-1i] ***

              DO lb=lb_min_set,lb_max-db-1
                DO bx=0,lb
                  fbx = REAL(bx,dp)
                  DO by=0,lb-bx
                    fby = REAL(by,dp)
                    bz = lb - bx - by
                    fbz = REAL(bz,dp)
                    cob = coset(bx,by,bz)
                    cobmx = coset(bx-1,by,bz)
                    cobmy = coset(bx,by-1,bz)
                    cobmz = coset(bx,by,bz-1)
                    cobpx = coset(bx+1,by,bz)
                    cobpy = coset(bx,by+1,bz)
                    cobpz = coset(bx,by,bz+1)
                    DO coa=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                       fs(coa,cob,cdbx) = ftz*s(coa,cobpx,cdb) - fbx*s(coa,cobmx,cdb)
                       fs(coa,cob,cdby) = ftz*s(coa,cobpy,cdb) - fby*s(coa,cobmy,cdb)
                       fs(coa,cob,cdbz) = ftz*s(coa,cobpz,cdb) - fbz*s(coa,cobmz,cdb)
                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO

!       *** Calculate the force contribution for the atomic center b ***

        IF (calculate_force_b) THEN
          DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
            DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
              force_b(1) = force_b(1) + pab(na+i,nb+j)*fs(i,j,2)
              force_b(2) = force_b(2) + pab(na+i,nb+j)*fs(i,j,3)
              force_b(3) = force_b(3) + pab(na+i,nb+j)*fs(i,j,4)
            END DO
          END DO
        END IF

        nb = nb + ncoset(lb_max_set)

      END DO

      na = na + ncoset(la_max_set)

    END DO

    DEALLOCATE (iiap,STAT=istat)
    CPASSERT(istat==0)

  END SUBROUTINE os_3center
! *****************************************************************************

END MODULE ai_oneelectron
