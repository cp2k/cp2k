!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the overlap integrals over Cartesian Gaussian-type
!>      functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - dab       : Distance between the atomic centers a and b.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - sab       : Shell set of overlap integrals.
!>      - zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      - zetp      : Reciprocal of the sum of the exponents of orbital a and b.
! *****************************************************************************
MODULE ai_overlap_aabb
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap_aabb'

! *** Public subroutines ***
  PUBLIC :: overlap_aabb

CONTAINS

! *****************************************************************************
!> \brief   Purpose: Calculation of the two-center overlap integrals [aa|bb] 
!>          over Cartesian Gaussian-type functions.
!> \author  Dorothea Golze
!> \date    06.2014
! *****************************************************************************
  SUBROUTINE overlap_aabb(la_max_set1,la_min_set1,npgfa1,rpgfa1,zeta1,&
                          la_max_set2,la_min_set2,npgfa2,rpgfa2,zeta2,&
                          lb_max_set1,lb_min_set1,npgfb1,rpgfb1,zetb1,&
                          lb_max_set2,lb_min_set2,npgfb2,rpgfb2,zetb2,&
                          asets_equal,bsets_equal,rab,dab,saabb,s,lds)

    INTEGER, INTENT(IN)                      :: la_max_set1, la_min_set1, &
                                                npgfa1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa1, zeta1
    INTEGER, INTENT(IN)                      :: la_max_set2, la_min_set2, &
                                                npgfa2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa2, zeta2
    INTEGER, INTENT(IN)                      :: lb_max_set1, lb_min_set1, &
                                                npgfb1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb1, zetb1
    INTEGER, INTENT(IN)                      :: lb_max_set2, lb_min_set2, &
                                                npgfb2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb2, zetb2
    LOGICAL, INTENT(IN)                      :: asets_equal, bsets_equal
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT)                          :: saabb
    INTEGER, INTENT(IN)                      :: lds
    REAL(KIND=dp), DIMENSION(lds, lds), &
      INTENT(INOUT)                          :: s

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_aabb', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, coa, coamx, coamy, coamz, coapx, &
      coapy, coapz, cob, cobm2x, cobm2y, cobm2z, cobmx, cobmy, cobmz, handle, &
      i, ia, ib, ipgf, j, ja, jb, jpgf, jpgf_start, kpgf, la, la_max, la_min, &
      la_start, lb, lb_max, lb_min, lpgf, lpgf_start, ncoa1, ncoa2, ncob1, &
      ncob2
    INTEGER, DIMENSION(3)                    :: na, naa, nb, nbb, nia, nib, &
                                                nja, njb
    REAL(KIND=dp)                            :: f0, f1, f2, f3, f4, fax, fay, &
                                                faz, zeta, zetb, zetp
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

    CALL timeset(routineN,handle)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    ncoa1 = 0
    ncoa2 = 0
    ncob1 = 0
    ncob2 = 0

    DO ipgf=1,npgfa1

      ncoa2 = 0

      IF(asets_equal) THEN
        jpgf_start=ipgf
        DO i = 1, jpgf_start-1
          ncoa2 = ncoa2 + ncoset(la_max_set2)
        ENDDO
      ELSE
        jpgf_start = 1
       ENDIF

      DO jpgf=jpgf_start,npgfa2

        ncob1 = 0
        zeta = zeta1(ipgf) + zeta2(jpgf)
        la_max = la_max_set1 + la_max_set2
        la_min = la_min_set1 + la_min_set2

        DO kpgf=1,npgfb1

          ncob2 = 0

          IF(bsets_equal) THEN
            lpgf_start=kpgf
            DO i = 1, lpgf_start-1
              ncob2 = ncob2 + ncoset(lb_max_set2)
            ENDDO
          ELSE
            lpgf_start = 1
          ENDIF

          DO lpgf=lpgf_start,npgfb2
 
            ! *** Screening ***
            IF ((rpgfa1(ipgf) + rpgfb1(kpgf) < dab).OR.&
                (rpgfa2(jpgf) + rpgfb1(kpgf) < dab).OR.&
                (rpgfa1(ipgf) + rpgfb2(lpgf) < dab).OR.&
                (rpgfa2(jpgf) + rpgfb2(lpgf) < dab)) THEN
               DO jb=ncoset(lb_min_set2-1)+1,ncoset(lb_max_set2)
                 DO ib=ncoset(lb_min_set1-1)+1,ncoset(lb_max_set1)
                   DO ja=ncoset(la_min_set2-1)+1,ncoset(la_max_set2)
                     DO ia =ncoset(la_min_set1-1)+1,ncoset(la_max_set1)
                       saabb(ncoa1+ia,ncoa2+ja,ncob1+ib,ncob2+jb) = 0._dp 
                       IF(asets_equal) saabb(ncoa2+ja,ncoa1+ia,ncob1+ib,ncob2+jb) = 0._dp
                       IF(bsets_equal) saabb(ncoa1+ia,ncoa2+ja,ncob2+jb,ncob1+ib) = 0._dp
                       IF(asets_equal.AND.bsets_equal) THEN
                         saabb(ncoa2+ja,ncoa1+ia,ncob1+ib,ncob2+jb) = 0._dp
                         saabb(ncoa1+ia,ncoa2+ja,ncob2+jb,ncob1+ib) = 0._dp
                         saabb(ncoa2+ja,ncoa1+ia,ncob2+jb,ncob1+ib) = 0._dp
                       ENDIF
                     ENDDO  
                   ENDDO
                 ENDDO
               ENDDO
               ncob2 = ncob2 + ncoset(lb_max_set2)
               CYCLE
            END IF

            zetb = zetb1(kpgf) + zetb2(lpgf)
            lb_max = lb_max_set1 + lb_max_set2
            lb_min = lb_min_set1 + lb_min_set2

!           *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta + zetb)

            f0 = SQRT((pi*zetp)**3)
            f1 = zetb*zetp
            f2 = 0.5_dp*zetp

!           *** Calculate the basic two-center overlap integral [s|s] ***

            s(1,1) = f0*EXP(-zeta*f1*dab*dab)

!           *** Recurrence steps: [s|s] -> [a|b] ***

            IF (la_max > 0) THEN

!             *** Vertical recurrence steps: [s|s] -> [a|s] ***

              rap(:) = f1*rab(:)

!             *** [p|s] = (Pi - Ai)*[s|s]  (i = x,y,z) ***

              s(2,1) = rap(1)*s(1,1) ! [px|s]
              s(3,1) = rap(2)*s(1,1) ! [py|s]
              s(4,1) = rap(3)*s(1,1) ! [pz|s]

              IF (la_max > 1) THEN

!               *** [d|s] ***

                f3 = f2*s(1,1)

                s( 5,1) = rap(1)*s(2,1) + f3 ! [dx2|s]
                s( 6,1) = rap(1)*s(3,1)      ! [dxy|s]
                s( 7,1) = rap(1)*s(4,1)      ! [dxz|s]
                s( 8,1) = rap(2)*s(3,1) + f3 ! [dy2|s]
                s( 9,1) = rap(2)*s(4,1)      ! [dyz|s]
                s(10,1) = rap(3)*s(4,1) + f3 ! [dz2|s]

                IF (la_max > 2) THEN

!                 *** [f|s] ***

                  f3 = 2.0_dp*f2

                  s(11,1) = rap(1)*s( 5,1) + f3*s(2,1) ! [fx3 |s]
                  s(12,1) = rap(1)*s( 6,1) + f2*s(3,1) ! [fx2y|s]
                  s(13,1) = rap(1)*s( 7,1) + f2*s(4,1) ! [fx2z|s]
                  s(14,1) = rap(1)*s( 8,1)             ! [fxy2|s]
                  s(15,1) = rap(1)*s( 9,1)             ! [fxyz|s]
                  s(16,1) = rap(1)*s(10,1)             ! [fxz2|s]
                  s(17,1) = rap(2)*s( 8,1) + f3*s(3,1) ! [fy3 |s]
                  s(18,1) = rap(2)*s( 9,1) + f2*s(4,1) ! [fy2z|s]
                  s(19,1) = rap(2)*s(10,1)             ! [fyz2|s]
                  s(20,1) = rap(3)*s(10,1) + f3*s(4,1) ! [fz3 |s]

                  IF (la_max > 3) THEN

!                   *** [g|s] ***

                    f4 = 3.0_dp*f2

                    s(21,1) = rap(1)*s(11,1) + f4*s( 5,1) ! [gx4  |s]
                    s(22,1) = rap(1)*s(12,1) + f3*s( 6,1) ! [gx3y |s]
                    s(23,1) = rap(1)*s(13,1) + f3*s( 7,1) ! [gx3z |s]
                    s(24,1) = rap(1)*s(14,1) + f2*s( 8,1) ! [gx2y2|s]
                    s(25,1) = rap(1)*s(15,1) + f2*s( 9,1) ! [gx2yz|s]
                    s(26,1) = rap(1)*s(16,1) + f2*s(10,1) ! [gx2z2|s]
                    s(27,1) = rap(1)*s(17,1)              ! [gxy3 |s]
                    s(28,1) = rap(1)*s(18,1)              ! [gxy2z|s]
                    s(29,1) = rap(1)*s(19,1)              ! [gxyz2|s]
                    s(30,1) = rap(1)*s(20,1)              ! [gxz3 |s]
                    s(31,1) = rap(2)*s(17,1) + f4*s( 8,1) ! [gy4  |s]
                    s(32,1) = rap(2)*s(18,1) + f3*s( 9,1) ! [gy3z |s]
                    s(33,1) = rap(2)*s(19,1) + f2*s(10,1) ! [gy2z2|s]
                    s(34,1) = rap(2)*s(20,1)              ! [gyz3 |s]
                    s(35,1) = rap(3)*s(20,1) + f4*s(10,1) ! [gz4  |s]

!                   *** [a|s] = (Pi - Ai)*[a-1i|s] + f2*Ni(a-1i)*[a-2i|s] ***

                    DO la=5,la_max

!                     *** Increase the angular momentum component z of a ***

                      s(coset(0,0,la),1) =&
                        rap(3)*s(coset(0,0,la-1),1) +&
                        f2*REAL(la-1,dp)*s(coset(0,0,la-2),1)

!                     *** Increase the angular momentum component y of a ***

                      az = la - 1
                      s(coset(0,1,az),1) = rap(2)*s(coset(0,0,az),1)
                      DO ay=2,la
                        az = la - ay
                        s(coset(0,ay,az),1) =&
                        rap(2)*s(coset(0,ay-1,az),1) +&
                        f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1)
                      END DO

!                     *** Increase the angular momentum component x of a ***

                      DO ay=0,la-1
                        az = la - 1 - ay
                        s(coset(1,ay,az),1) = rap(1)*s(coset(0,ay,az),1)
                      END DO
                      DO ax=2,la
                        f3 = f2*REAL(ax-1,dp)
                        DO ay=0,la-ax
                          az = la - ax - ay
                          s(coset(ax,ay,az),1) =&
                          rap(1)*s(coset(ax-1,ay,az),1) +&
                          f3*s(coset(ax-2,ay,az),1)
                        END DO
                      END DO

                    END DO

                  END IF

                END IF

              END IF

!             *** Recurrence steps: [a|s] -> [a|b] ***

              IF (lb_max > 0) THEN

                DO j=2,ncoset(lb_max)
                  DO i=1,ncoset(la_min)
                    s(i,j) = 0.0_dp
                  END DO
                END DO

!               *** Horizontal recurrence steps ***

                rbp(:) = rap(:) - rab(:)

!               *** [a|p] = [a+1i|s] - (Bi - Ai)*[a|s] ***

                IF (lb_max == 1) THEN
                  la_start = la_min
                ELSE
                  la_start = MAX(0,la_min-1)
                END IF

                DO la=la_start,la_max-1
                  DO ax=0,la
                    DO ay=0,la-ax
                      az = la - ax - ay
                      coa = coset(ax,ay,az)
                      coapx = coset(ax+1,ay,az)
                      coapy = coset(ax,ay+1,az)
                      coapz = coset(ax,ay,az+1)
                      s(coa,2) = s(coapx,1) - rab(1)*s(coa,1)
                      s(coa,3) = s(coapy,1) - rab(2)*s(coa,1)
                      s(coa,4) = s(coapz,1) - rab(3)*s(coa,1)
                    END DO
                  END DO
                END DO

!               *** Vertical recurrence step ***

!               *** [a|p] = (Pi - Bi)*[a|s] + f2*Ni(a)*[a-1i|s] ***

                DO ax=0,la_max
                  fax = f2*REAL(ax,dp)
                  DO ay=0,la_max-ax
                    fay = f2*REAL(ay,dp)
                    az = la_max - ax - ay
                    faz = f2*REAL(az,dp)
                    coa = coset(ax,ay,az)
                    coamx = coset(ax-1,ay,az)
                    coamy = coset(ax,ay-1,az)
                    coamz = coset(ax,ay,az-1)
                    s(coa,2) = rbp(1)*s(coa,1) + fax*s(coamx,1)
                    s(coa,3) = rbp(2)*s(coa,1) + fay*s(coamy,1)
                    s(coa,4) = rbp(3)*s(coa,1) + faz*s(coamz,1)
                  END DO
                END DO

!               *** Recurrence steps: [a|p] -> [a|b] ***

                DO lb=2,lb_max

!                 *** Horizontal recurrence steps ***

!                 *** [a|b] = [a+1i|b-1i] - (Bi - Ai)*[a|b-1i] ***

                  IF (lb == lb_max) THEN
                    la_start = la_min
                  ELSE
                    la_start = MAX(0,la_min-1)
                  END IF

                  DO la=la_start,la_max-1
                    DO ax=0,la
                      DO ay=0,la-ax
                        az = la - ax - ay
                        coa = coset(ax,ay,az)
                        coapx = coset(ax+1,ay,az)
                        coapy = coset(ax,ay+1,az)
                        coapz = coset(ax,ay,az+1)

!                       *** Shift of angular momentum component z from a to b ***

                        cob = coset(0,0,lb)
                        cobmz = coset(0,0,lb-1)
                        s(coa,cob) = s(coapz,cobmz) - rab(3)*s(coa,cobmz)

!                       *** Shift of angular momentum component y from a to b ***

                        DO by=1,lb
                          bz = lb - by
                          cob = coset(0,by,bz)
                          cobmy = coset(0,by-1,bz)
                          s(coa,cob) = s(coapy,cobmy) - rab(2)*s(coa,cobmy)
                        END DO

!                       *** Shift of angular momentum component x from a to b ***

                        DO bx=1,lb
                          DO by=0,lb-bx
                            bz = lb - bx - by
                            cob = coset(bx,by,bz)
                            cobmx = coset(bx-1,by,bz)
                            s(coa,cob) = s(coapx,cobmx) - rab(1)*s(coa,cobmx)
                          END DO
                        END DO

                      END DO
                    END DO
                  END DO

!                 *** Vertical recurrence step ***

!                 *** [a|b] = (Pi - Bi)*[a|b-1i] + f2*Ni(a)*[a-1i|b-1i] + ***
!                 ***         f2*Ni(b-1i)*[a|b-2i]                        ***

                  DO ax=0,la_max
                    fax = f2*REAL(ax,dp)
                    DO ay=0,la_max-ax
                      fay = f2*REAL(ay,dp)
                      az = la_max - ax - ay
                      faz = f2*REAL(az,dp)
                      coa = coset(ax,ay,az)
                      coamx = coset(ax-1,ay,az)
                      coamy = coset(ax,ay-1,az)
                      coamz = coset(ax,ay,az-1)

!                     *** Increase the angular momentum component z of b ***

                      f3 = f2*REAL(lb-1,dp)
                      cob = coset(0,0,lb)
                      cobmz = coset(0,0,lb-1)
                      cobm2z = coset(0,0,lb-2)
                      s(coa,cob) = rbp(3)*s(coa,cobmz) +&
                                   faz*s(coamz,cobmz) +&
                                   f3*s(coa,cobm2z)

!                     *** Increase the angular momentum component y of b ***

                      bz = lb - 1
                      cob = coset(0,1,bz)
                      cobmy = coset(0,0,bz)
                      s(coa,cob) = rbp(2)*s(coa,cobmy) +&
                                     fay*s(coamy,cobmy)
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        cob = coset(0,by,bz)
                        cobmy = coset(0,by-1,bz)
                        cobm2y = coset(0,by-2,bz)
                        s(coa,cob) = rbp(2)*s(coa,cobmy) +&
                                     fay*s(coamy,cobmy) +&
                                     f3*s(coa,cobm2y)
                      END DO

!                     *** Increase the angular momentum component x of b ***

                      DO by=0,lb-1
                        bz = lb - 1 - by
                        cob = coset(1,by,bz)
                        cobmx = coset(0,by,bz)
                        s(coa,cob) = rbp(1)*s(coa,cobmx) +&
                                     fax*s(coamx,cobmx)
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          cob = coset(bx,by,bz)
                          cobmx = coset(bx-1,by,bz)
                          cobm2x = coset(bx-2,by,bz)
                          s(coa,cob) = rbp(1)*s(coa,cobmx) +&
                                       fax*s(coamx,cobmx) +&
                                       f3*s(coa,cobm2x)
                        END DO
                      END DO

                    END DO
                  END DO

                END DO

              END IF

            ELSE

              IF (lb_max > 0) THEN

!               *** Vertical recurrence steps: [s|s] -> [s|b] ***

                rbp(:) = (f1 - 1.0_dp)*rab(:)

!               *** [s|p] = (Pi - Bi)*[s|s] ***

                s(1,2) = rbp(1)*s(1,1) ! [s|px]
                s(1,3) = rbp(2)*s(1,1) ! [s|py]
                s(1,4) = rbp(3)*s(1,1) ! [s|pz]

                IF (lb_max > 1) THEN

!                 *** [s|d] ***

                  f3 = f2*s(1,1)

                  s(1, 5) = rbp(1)*s(1,2) + f3 ! [s|dx2]
                  s(1, 6) = rbp(1)*s(1,3)      ! [s|dxy]
                  s(1, 7) = rbp(1)*s(1,4)      ! [s|dxz]
                  s(1, 8) = rbp(2)*s(1,3) + f3 ! [s|dy2]
                  s(1, 9) = rbp(2)*s(1,4)      ! [s|dyz]
                  s(1,10) = rbp(3)*s(1,4) + f3 ! [s|dz2]

!                 *** [s|b] = (Pi - Bi)*[s|b-1i] + f2*Ni(b-1i)*[s|b-2i] ***

                  DO lb=3,lb_max

!                   *** Increase the angular momentum component z of b ***

                    s(1,coset(0,0,lb)) =&
                      rbp(3)*s(1,coset(0,0,lb-1)) +&
                      f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2))

!                   *** Increase the angular momentum component y of b ***

                    bz = lb - 1
                    s(1,coset(0,1,bz)) = rbp(2)*s(1,coset(0,0,bz))
                    DO by=2,lb
                      bz = lb - by
                      s(1,coset(0,by,bz)) =&
                        rbp(2)*s(1,coset(0,by-1,bz)) +&
                        f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz))
                    END DO

!                   *** Increase the angular momentum component x of b ***

                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(1,coset(1,by,bz)) = rbp(1)*s(1,coset(0,by,bz))
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(1,coset(bx,by,bz)) =&
                          rbp(1)*s(1,coset(bx-1,by,bz)) +&
                          f3*s(1,coset(bx-2,by,bz))
                      END DO
                    END DO

                  END DO

                END IF

              END IF

            END IF
               
!           *** Store the primitive overlap integrals ***
            DO jb=ncoset(lb_min_set2-1)+1,ncoset(lb_max_set2)
              njb(1:3) = indco(1:3,jb)
              DO ib=ncoset(lb_min_set1-1)+1,ncoset(lb_max_set1)
                nib(1:3) = indco(1:3,ib)
                nbb(1:3) = nib + njb
                DO ja=ncoset(la_min_set2-1)+1,ncoset(la_max_set2)
                  nja(1:3) = indco(1:3,ja)
                  DO ia =ncoset(la_min_set1-1)+1,ncoset(la_max_set1)
                    nia(1:3) = indco(1:3,ia)
                    naa(1:3) = nia + nja
                    ! now loop over all elements of s
                    DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
                       nb(1:3) = indco(1:3,j)
                      DO i=ncoset(la_min-1)+1,ncoset(la_max)
                        na(1:3) = indco(1:3,i)
                        IF(ALL(na == naa).AND.ALL(nb == nbb)) THEN
                          saabb(ncoa1+ia,ncoa2+ja,ncob1+ib,ncob2+jb) = s(i,j)
                          IF(asets_equal) saabb(ncoa2+ja,ncoa1+ia,ncob1+ib,ncob2+jb) = s(i,j)
                          IF(bsets_equal) saabb(ncoa1+ia,ncoa2+ja,ncob2+jb,ncob1+ib) = s(i,j)
                          IF(asets_equal.AND.bsets_equal) THEN
                            saabb(ncoa2+ja,ncoa1+ia,ncob1+ib,ncob2+jb) = s(i,j)
                            saabb(ncoa1+ia,ncoa2+ja,ncob2+jb,ncob1+ib) = s(i,j)
                            saabb(ncoa2+ja,ncoa1+ia,ncob2+jb,ncob1+ib) = s(i,j) 
                          ENDIF
                        ENDIF 
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO

            ncob2 = ncob2 + ncoset(lb_max_set2)

          END DO

          ncob1 = ncob1 + ncoset(lb_max_set1)

        END DO

        ncoa2 = ncoa2 + ncoset(la_max_set2)

      END DO

      ncoa1 = ncoa1 + ncoset(la_max_set1)

    END DO

    CALL timestop(handle)

  END SUBROUTINE overlap_aabb

END MODULE ai_overlap_aabb
