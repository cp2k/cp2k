!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the overlap integrals over Cartesian Gaussian-type
!>      functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Derivatives added (02.05.2002,MK)
!>      - New OS routine with simpler logic (11.07.2014, JGH)
!> \author Matthias Krack (08.10.1999)
! *****************************************************************************
MODULE ai_overlap
  USE ai_os_rr,                        ONLY: os_rr_ovlp
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap'

! *** Public subroutines ***
  PUBLIC :: overlap, overlap_ab

CONTAINS

! *****************************************************************************
!> \brief   Purpose: Calculation of the two-center overlap integrals [a|b] over
!>          Cartesian Gaussian-type functions.
!> \param la_max_set Max L on center A
!> \param la_min_set Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb_max_set Max L on center B
!> \param lb_min_set Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param dab        Distance A-B
!> \param sab        Final Integrals, basic and derivatives
!> \param da_max_set Some additional derivative information
!> \param return_derivatives   Return integral derivatives
!> \param s          Work space
!> \param lds        Leading dimension of s
!> \param sdab       Return additional derivative integrals
!> \param pab        Density matrix block, used to calculate forces
!> \param force_a    Force vector [da/dR|b]
!> \date    19.09.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE overlap(la_max_set,la_min_set,npgfa,rpgfa,zeta,&
                     lb_max_set,lb_min_set,npgfb,rpgfb,zetb,&
                     rab,dab,sab,da_max_set,return_derivatives,s,lds,&
                     sdab,pab,force_a)
    INTEGER, INTENT(IN)                      :: la_max_set, la_min_set, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max_set, lb_min_set, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: sab
    INTEGER, INTENT(IN)                      :: da_max_set
    LOGICAL, INTENT(IN)                      :: return_derivatives
    INTEGER, INTENT(IN)                      :: lds
    REAL(KIND=dp), DIMENSION(lds, lds, *), &
      INTENT(INOUT)                          :: s
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: sdab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pab
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: force_a

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coamx, &
      coamy, coamz, coapx, coapy, coapz, cob, cobm2x, cobm2y, cobm2z, cobmx, &
      cobmy, cobmz, da, da_max, dax, day, daz, i, ipgf, j, jk, jpgf, jstart, &
      k, la, la_max, la_min, la_start, lb, lb_max, lb_min, lb_start, na, nb, &
      nda
    LOGICAL                                  :: calculate_force_a
    REAL(KIND=dp)                            :: f0, f1, f2, f3, f4, fax, fay, &
                                                faz, ftz, zetp
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

    IF (PRESENT(pab).AND.PRESENT(force_a)) THEN
      calculate_force_a = .TRUE.
      force_a(:) = 0.0_dp
    ELSE
      calculate_force_a = .FALSE.
    END IF

    IF (PRESENT(sdab).OR.calculate_force_a) THEN
      IF (da_max_set == 0) THEN
        da_max = 1
        la_max = la_max_set + 1
        la_min = MAX(0,la_min_set-1)
      ELSE
        da_max = da_max_set
        la_max = la_max_set + da_max_set + 1
        la_min = MAX(0,la_min_set-da_max_set-1)
      END IF
    ELSE
      da_max = da_max_set
      la_max = la_max_set + da_max_set
      la_min = MAX(0,la_min_set-da_max_set)
    END IF

    lb_max = lb_max_set
    lb_min = lb_min_set

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nda = 0

    DO ipgf=1,npgfa

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb+1,nb+ncoset(lb_max_set)
            DO i=na+1,na+ncoset(la_max_set)
              sab(i,j) = 0.0_dp
            END DO
          END DO
          IF (return_derivatives) THEN
            DO k=2,ncoset(da_max_set)
              jstart = (k - 1)*SIZE(sab,1)
              DO j=jstart+nb+1,jstart+nb+ncoset(lb_max_set)
                DO i=na+1,na+ncoset(la_max_set)
                  sab(i,j) = 0.0_dp
                END DO
              END DO
            END DO
          END IF
          nb = nb + ncoset(lb_max_set)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = SQRT((pi*zetp)**3)
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp

!       *** Calculate the basic two-center overlap integral [s|s] ***

        s(1,1,1) = f0*EXP(-zeta(ipgf)*f1*dab*dab)

!       *** Recurrence steps: [s|s] -> [a|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|s] -> [a|s] ***

          rap(:) = f1*rab(:)

!         *** [p|s] = (Pi - Ai)*[s|s]  (i = x,y,z) ***

          s(2,1,1) = rap(1)*s(1,1,1) ! [px|s]
          s(3,1,1) = rap(2)*s(1,1,1) ! [py|s]
          s(4,1,1) = rap(3)*s(1,1,1) ! [pz|s]

          IF (la_max > 1) THEN

!           *** [d|s] ***

            f3 = f2*s(1,1,1)

            s( 5,1,1) = rap(1)*s(2,1,1) + f3 ! [dx2|s]
            s( 6,1,1) = rap(1)*s(3,1,1)      ! [dxy|s]
            s( 7,1,1) = rap(1)*s(4,1,1)      ! [dxz|s]
            s( 8,1,1) = rap(2)*s(3,1,1) + f3 ! [dy2|s]
            s( 9,1,1) = rap(2)*s(4,1,1)      ! [dyz|s]
            s(10,1,1) = rap(3)*s(4,1,1) + f3 ! [dz2|s]

            IF (la_max > 2) THEN

!             *** [f|s] ***

              f3 = 2.0_dp*f2

              s(11,1,1) = rap(1)*s( 5,1,1) + f3*s(2,1,1) ! [fx3 |s]
              s(12,1,1) = rap(1)*s( 6,1,1) + f2*s(3,1,1) ! [fx2y|s]
              s(13,1,1) = rap(1)*s( 7,1,1) + f2*s(4,1,1) ! [fx2z|s]
              s(14,1,1) = rap(1)*s( 8,1,1)               ! [fxy2|s]
              s(15,1,1) = rap(1)*s( 9,1,1)               ! [fxyz|s]
              s(16,1,1) = rap(1)*s(10,1,1)               ! [fxz2|s]
              s(17,1,1) = rap(2)*s( 8,1,1) + f3*s(3,1,1) ! [fy3 |s]
              s(18,1,1) = rap(2)*s( 9,1,1) + f2*s(4,1,1) ! [fy2z|s]
              s(19,1,1) = rap(2)*s(10,1,1)               ! [fyz2|s]
              s(20,1,1) = rap(3)*s(10,1,1) + f3*s(4,1,1) ! [fz3 |s]

              IF (la_max > 3) THEN

!               *** [g|s] ***

                f4 = 3.0_dp*f2

                s(21,1,1) = rap(1)*s(11,1,1) + f4*s( 5,1,1) ! [gx4  |s]
                s(22,1,1) = rap(1)*s(12,1,1) + f3*s( 6,1,1) ! [gx3y |s]
                s(23,1,1) = rap(1)*s(13,1,1) + f3*s( 7,1,1) ! [gx3z |s]
                s(24,1,1) = rap(1)*s(14,1,1) + f2*s( 8,1,1) ! [gx2y2|s]
                s(25,1,1) = rap(1)*s(15,1,1) + f2*s( 9,1,1) ! [gx2yz|s]
                s(26,1,1) = rap(1)*s(16,1,1) + f2*s(10,1,1) ! [gx2z2|s]
                s(27,1,1) = rap(1)*s(17,1,1)                ! [gxy3 |s]
                s(28,1,1) = rap(1)*s(18,1,1)                ! [gxy2z|s]
                s(29,1,1) = rap(1)*s(19,1,1)                ! [gxyz2|s]
                s(30,1,1) = rap(1)*s(20,1,1)                ! [gxz3 |s]
                s(31,1,1) = rap(2)*s(17,1,1) + f4*s( 8,1,1) ! [gy4  |s]
                s(32,1,1) = rap(2)*s(18,1,1) + f3*s( 9,1,1) ! [gy3z |s]
                s(33,1,1) = rap(2)*s(19,1,1) + f2*s(10,1,1) ! [gy2z2|s]
                s(34,1,1) = rap(2)*s(20,1,1)                ! [gyz3 |s]
                s(35,1,1) = rap(3)*s(20,1,1) + f4*s(10,1,1) ! [gz4  |s]

!               *** [a|s] = (Pi - Ai)*[a-1i|s] + f2*Ni(a-1i)*[a-2i|s] ***

                DO la=5,la_max

!                 *** Increase the angular momentum component z of a ***

                  s(coset(0,0,la),1,1) =&
                    rap(3)*s(coset(0,0,la-1),1,1) +&
                    f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,1)

!                 *** Increase the angular momentum component y of a ***

                  az = la - 1
                  s(coset(0,1,az),1,1) = rap(2)*s(coset(0,0,az),1,1)
                  DO ay=2,la
                    az = la - ay
                    s(coset(0,ay,az),1,1) =&
                    rap(2)*s(coset(0,ay-1,az),1,1) +&
                    f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,1)
                  END DO

!                 *** Increase the angular momentum component x of a ***

                  DO ay=0,la-1
                    az = la - 1 - ay
                    s(coset(1,ay,az),1,1) = rap(1)*s(coset(0,ay,az),1,1)
                  END DO
                  DO ax=2,la
                    f3 = f2*REAL(ax-1,dp)
                    DO ay=0,la-ax
                      az = la - ax - ay
                      s(coset(ax,ay,az),1,1) =&
                      rap(1)*s(coset(ax-1,ay,az),1,1) +&
                      f3*s(coset(ax-2,ay,az),1,1)
                    END DO
                  END DO

                END DO

              END IF

            END IF

          END IF

!         *** Recurrence steps: [a|s] -> [a|b] ***

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_min)
                s(i,j,1) = 0.0_dp
              END DO
            END DO

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a|p] = [a+1i|s] - (Bi - Ai)*[a|s] ***

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
                  s(coa,2,1) = s(coapx,1,1) - rab(1)*s(coa,1,1)
                  s(coa,3,1) = s(coapy,1,1) - rab(2)*s(coa,1,1)
                  s(coa,4,1) = s(coapz,1,1) - rab(3)*s(coa,1,1)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|p] = (Pi - Bi)*[a|s] + f2*Ni(a)*[a-1i|s] ***

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
                s(coa,2,1) = rbp(1)*s(coa,1,1) + fax*s(coamx,1,1)
                s(coa,3,1) = rbp(2)*s(coa,1,1) + fay*s(coamy,1,1)
                s(coa,4,1) = rbp(3)*s(coa,1,1) + faz*s(coamz,1,1)
              END DO
            END DO

!           *** Recurrence steps: [a|p] -> [a|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|b] = [a+1i|b-1i] - (Bi - Ai)*[a|b-1i] ***

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

!                   *** Shift of angular momentum component z from a to b ***

                    cob = coset(0,0,lb)
                    cobmz = coset(0,0,lb-1)
                    s(coa,cob,1) = s(coapz,cobmz,1) - rab(3)*s(coa,cobmz,1)

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      cob = coset(0,by,bz)
                      cobmy = coset(0,by-1,bz)
                      s(coa,cob,1) = s(coapy,cobmy,1) - rab(2)*s(coa,cobmy,1)
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        cob = coset(bx,by,bz)
                        cobmx = coset(bx-1,by,bz)
                        s(coa,cob,1) = s(coapx,cobmx,1) - rab(1)*s(coa,cobmx,1)
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|b] = (Pi - Bi)*[a|b-1i] + f2*Ni(a)*[a-1i|b-1i] + ***
!             ***         f2*Ni(b-1i)*[a|b-2i]                        ***

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

!                 *** Increase the angular momentum component z of b ***

                  f3 = f2*REAL(lb-1,dp)
                  cob = coset(0,0,lb)
                  cobmz = coset(0,0,lb-1)
                  cobm2z = coset(0,0,lb-2)
                  s(coa,cob,1) = rbp(3)*s(coa,cobmz,1) +&
                                 faz*s(coamz,cobmz,1) +&
                                 f3*s(coa,cobm2z,1)

!                 *** Increase the angular momentum component y of b ***

                  bz = lb - 1
                  cob = coset(0,1,bz)
                  cobmy = coset(0,0,bz)
                  s(coa,cob,1) = rbp(2)*s(coa,cobmy,1) +&
                                 fay*s(coamy,cobmy,1)
                  DO by=2,lb
                    bz = lb - by
                    f3 = f2*REAL(by-1,dp)
                    cob = coset(0,by,bz)
                    cobmy = coset(0,by-1,bz)
                    cobm2y = coset(0,by-2,bz)
                    s(coa,cob,1) = rbp(2)*s(coa,cobmy,1) +&
                                   fay*s(coamy,cobmy,1) +&
                                   f3*s(coa,cobm2y,1)
                  END DO

!                 *** Increase the angular momentum component x of b ***

                  DO by=0,lb-1
                    bz = lb - 1 - by
                    cob = coset(1,by,bz)
                    cobmx = coset(0,by,bz)
                    s(coa,cob,1) = rbp(1)*s(coa,cobmx,1) +&
                                   fax*s(coamx,cobmx,1)
                  END DO
                  DO bx=2,lb
                    f3 = f2*REAL(bx-1,dp)
                    DO by=0,lb-bx
                      bz = lb - bx - by
                      cob = coset(bx,by,bz)
                      cobmx = coset(bx-1,by,bz)
                      cobm2x = coset(bx-2,by,bz)
                      s(coa,cob,1) = rbp(1)*s(coa,cobmx,1) +&
                                     fax*s(coamx,cobmx,1) +&
                                     f3*s(coa,cobm2x,1)
                    END DO
                  END DO

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|s] -> [s|b] ***

            rbp(:) = (f1 - 1.0_dp)*rab(:)

!           *** [s|p] = (Pi - Bi)*[s|s] ***

            s(1,2,1) = rbp(1)*s(1,1,1) ! [s|px]
            s(1,3,1) = rbp(2)*s(1,1,1) ! [s|py]
            s(1,4,1) = rbp(3)*s(1,1,1) ! [s|pz]

            IF (lb_max > 1) THEN

!             *** [s|d] ***

              f3 = f2*s(1,1,1)

              s(1, 5,1) = rbp(1)*s(1,2,1) + f3 ! [s|dx2]
              s(1, 6,1) = rbp(1)*s(1,3,1)      ! [s|dxy]
              s(1, 7,1) = rbp(1)*s(1,4,1)      ! [s|dxz]
              s(1, 8,1) = rbp(2)*s(1,3,1) + f3 ! [s|dy2]
              s(1, 9,1) = rbp(2)*s(1,4,1)      ! [s|dyz]
              s(1,10,1) = rbp(3)*s(1,4,1) + f3 ! [s|dz2]

!             *** [s|b] = (Pi - Bi)*[s|b-1i] + f2*Ni(b-1i)*[s|b-2i] ***

              DO lb=3,lb_max

!               *** Increase the angular momentum component z of b ***

                s(1,coset(0,0,lb),1) =&
                  rbp(3)*s(1,coset(0,0,lb-1),1) +&
                  f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),1)

!               *** Increase the angular momentum component y of b ***

                bz = lb - 1
                s(1,coset(0,1,bz),1) = rbp(2)*s(1,coset(0,0,bz),1)
                DO by=2,lb
                  bz = lb - by
                  s(1,coset(0,by,bz),1) =&
                    rbp(2)*s(1,coset(0,by-1,bz),1) +&
                    f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),1)
                END DO

!               *** Increase the angular momentum component x of b ***

                DO by=0,lb-1
                  bz = lb - 1 - by
                  s(1,coset(1,by,bz),1) = rbp(1)*s(1,coset(0,by,bz),1)
                END DO
                DO bx=2,lb
                  f3 = f2*REAL(bx-1,dp)
                  DO by=0,lb-bx
                    bz = lb - bx - by
                    s(1,coset(bx,by,bz),1) =&
                      rbp(1)*s(1,coset(bx-1,by,bz),1) +&
                      f3*s(1,coset(bx-2,by,bz),1)
                  END DO
                END DO

              END DO

            END IF

          END IF

        END IF

!       *** Store the primitive overlap integrals ***

        DO j=1,ncoset(lb_max_set)
          DO i=1,ncoset(la_max_set)
            sab(na+i,nb+j) = s(i,j,1)
          END DO
        END DO

!       *** Calculate the requested derivatives with respect  ***
!       *** to the nuclear coordinates of the atomic center a ***

        IF (PRESENT(sdab).OR.return_derivatives) THEN
          la_start = 0
          lb_start = 0
        ELSE
          la_start = la_min_set
          lb_start = lb_min_set
        END IF

        DO da=0,da_max-1
          ftz = 2.0_dp*zeta(ipgf)
          DO dax=0,da
            DO day=0,da-dax
              daz = da - dax - day
              cda = coset(dax,day,daz)
              cdax = coset(dax+1,day,daz)
              cday = coset(dax,day+1,daz)
              cdaz = coset(dax,day,daz+1)

!             *** [da/dAi|b] = 2*zeta*[a+1i|b] - Ni(a)[a-1i|b] ***

              DO la=la_start,la_max-da-1
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
                    DO lb=lb_start,lb_max_set
                      DO bx=0,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          cob = coset(bx,by,bz)
                          s(coa,cob,cdax) = ftz*s(coapx,cob,cda) -&
                                            fax*s(coamx,cob,cda)
                          s(coa,cob,cday) = ftz*s(coapy,cob,cda) -&
                                            fay*s(coamy,cob,cda)
                          s(coa,cob,cdaz) = ftz*s(coapz,cob,cda) -&
                                            faz*s(coamz,cob,cda)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO

!       *** Return all the calculated derivatives of the ***
!       *** primitive overlap integrals, if requested    ***

        IF (return_derivatives) THEN
          DO k=2,ncoset(da_max_set)
            jstart = (k - 1)*SIZE(sab,1)
            DO j=1,ncoset(lb_max_set)
              jk = jstart + j
              DO i=1,ncoset(la_max_set)
                sab(na+i,nb+jk) = s(i,j,k)
              END DO
            END DO
          END DO
        END IF

!       *** Calculate the force contribution for the atomic center a ***

        IF (calculate_force_a) THEN
          DO k=1,3
            DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
              DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                force_a(k) = force_a(k) + pab(na+i,nb+j)*s(i,j,k+1)
              END DO
            END DO
          END DO
        END IF

!       *** Store the first derivatives of the primitive overlap integrals ***
!       *** which are used as auxiliary integrals for the calculation of   ***
!       *** the kinetic energy integrals if requested                      ***

        IF (PRESENT(sdab)) THEN
          sdab(nda+1,nb+1,1) = s(1,1,1)
          DO k=2,4
            DO j=1,ncoset(lb_max_set)
              DO i=1,ncoset(la_max-1)
                sdab(nda+i,nb+j,k) = s(i,j,k)
              END DO
            END DO
          END DO
        END IF

        nb = nb + ncoset(lb_max_set)

      END DO

      na = na + ncoset(la_max_set)
      nda = nda + ncoset(la_max-1)

    END DO

  END SUBROUTINE overlap

! *****************************************************************************
!> \brief   Calculation of the two-center overlap integrals [a|b] over
!>          Cartesian Gaussian-type functions. First and second derivatives
!> \param la_max     Max L on center A
!> \param la_min     Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb_max     Max L on center B
!> \param lb_min     Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param sab        Final overlap integrals
!> \param dab        First derivative overlap integrals
!> \param ddab       Second derivative overlap integrals
!> \param error      CP2K error reporting
!> \date    01.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE overlap_ab(la_max,la_min,npgfa,rpgfa,zeta,&
                        lb_max,lb_min,npgfb,rpgfb,zetb,&
                        rab,sab,dab,ddab,error)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: sab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: dab, ddab
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_ab', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, coa, &
                                                cob, ia, ib, ipgf, istat, &
                                                jpgf, la, lb, ldrr, lma, lmb, &
                                                ma, mb, na, nb, ofa, ofb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, ambm, ambp, apbm, apbp, b, &
                                                dumx, dumy, dumz, f0, rab2, &
                                                tab, xhi, zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(sab)) THEN
       lma = la_max
       lmb = lb_max
    END IF
    IF (PRESENT(dab)) THEN
       lma = la_max + 1
       lmb = lb_max
    END IF
    IF (PRESENT(ddab)) THEN
       lma = la_max + 1
       lmb = lb_max + 1
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa = ncoset(la_min-1)
    ofb = ncoset(lb_min-1)
    na = ncoset(la_max) - ofa
    nb = ncoset(lb_max) - ofb
    IF (PRESENT(sab)) THEN
       CPPrecondition((SIZE(sab,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(sab,2) >= nb*npgfb),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(dab)) THEN
       CPPrecondition((SIZE(dab,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dab,2) >= nb*npgfb),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dab,3) >= 3),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(ddab)) THEN
       CPPrecondition((SIZE(ddab,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(ddab,2) >= nb*npgfb),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(ddab,3) >= 6),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all pairs of primitive Gaussian-type functions
    ma = 0
    DO ipgf=1,npgfa
       mb = 0
       DO jpgf=1,npgfb
          ! Distance Screening
          IF(rpgfa(ipgf) + rpgfb(jpgf) < tab) THEN
             IF (PRESENT(sab)) sab(ma+1:ma+na,mb+1:mb+nb) = 0.0_dp
             IF (PRESENT(dab)) dab(ma+1:ma+na,mb+1:mb+nb,1:3) = 0.0_dp
             IF (PRESENT(ddab)) ddab(ma+1:ma+na,mb+1:mb+nb,1:6) = 0.0_dp
             mb = mb + nb
             CYCLE
          ENDIF

          ! Calculate some prefactors
          a = zeta(ipgf)
          b = zetb(jpgf)
          zet =  a+b
          xhi =  a*b/zet
          rap =  b*rab/zet
          rbp = -a*rab/zet

          ! [s|s] integral
          f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

          ! Calculate the reccurence relation
          CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

          DO lb = lb_min,lb_max
          DO bx = 0,lb
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob = coset(bx,by,bz) - ofb
             ib = mb + cob
             DO la = la_min,la_max
             DO ax = 0,la
             DO ay = 0,la-ax
                az = la - ax - ay
                coa = coset(ax,ay,az) - ofa
                ia = ma + coa
                ! integrals
                IF (PRESENT(sab)) THEN
                   sab(ia,ib) = f0 * rr(ax,bx,1) * rr(ay,by,2) * rr(az,bz,3)
                END IF
                ! first derivatives
                IF (PRESENT(dab)) THEN
                   ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                   ! dx
                   dumx = 2.0_dp * a * rr(ax+1,bx,1)
                   IF(ax > 0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                   dab(ia,ib,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                   ! dy
                   dumy = 2.0_dp * a * rr(ay+1,by,2)
                   IF(ay > 0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                   dab(ia,ib,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                   ! dz
                   dumz = 2.0_dp * a * rr(az+1,bz,3)
                   IF(az > 0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                   dab(ia,ib,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                END IF
                ! 2nd derivatives
                IF (PRESENT(ddab)) THEN
                   ! (dda|b) = -4*a*b*(a+1|b+1) + 2*a*N(b)*(a+1|b-1)
                   !           + 2*b*N(a)*(a-1|b+1) - N(a)*N(b)*(a-1|b-1)
                   ! dx dx
                   apbp = f0 * rr(ax+1,bx+1,1) * rr(ay,by,2) * rr(az,bz,3)
                   IF(bx > 0) THEN
                       apbm = f0 * rr(ax+1,bx-1,1) * rr(ay,by,2) * rr(az,bz,3)
                   ELSE
                       apbm = 0.0_dp
                   END IF
                   IF(ax > 0) THEN
                       ambp = f0 * rr(ax-1,bx+1,1) * rr(ay,by,2) * rr(az,bz,3)
                   ELSE
                       ambp = 0.0_dp
                   END IF
                   IF(ax > 0 .AND. bx > 0) THEN
                       ambm = f0 * rr(ax-1,bx-1,1) * rr(ay,by,2) * rr(az,bz,3)
                   ELSE
                       ambm = 0.0_dp
                   END IF
                   ddab(ia,ib,1) = -4.0_dp*a*b*apbm + 2.0_dp*a*REAL(bx,dp)*apbm &
                      + 2.0_dp*b*REAL(ax,dp)*ambp - REAL(ax,dp)*REAL(bx,dp)*ambm
                   ! dx dy
                   apbp = f0 * rr(ax+1,bx,1) * rr(ay,by+1,2) * rr(az,bz,3)
                   IF(by > 0) THEN
                       apbm = f0 * rr(ax+1,bx,1) * rr(ay,by-1,2) * rr(az,bz,3)
                   ELSE
                       apbm = 0.0_dp
                   END IF
                   IF(ax > 0) THEN
                       ambp = f0 * rr(ax-1,bx,1) * rr(ay,by+1,2) * rr(az,bz,3)
                   ELSE
                       ambp = 0.0_dp
                   END IF
                   IF(ax > 0 .AND. by > 0) THEN
                       ambm = f0 * rr(ax-1,bx,1) * rr(ay,by-1,2) * rr(az,bz,3)
                   ELSE
                       ambm = 0.0_dp
                   END IF
                   ddab(ia,ib,2) = -4.0_dp*a*b*apbm + 2.0_dp*a*REAL(by,dp)*apbm &
                      + 2.0_dp*b*REAL(ax,dp)*ambp - REAL(ax,dp)*REAL(by,dp)*ambm
                   ! dx dz
                   apbp = f0 * rr(ax+1,bx,1) * rr(ay,by,2) * rr(az,bz+1,3)
                   IF(bz > 0) THEN
                       apbm = f0 * rr(ax+1,bx,1) * rr(ay,by,2) * rr(az,bz-1,3)
                   ELSE
                       apbm = 0.0_dp
                   END IF
                   IF(ax > 0) THEN
                       ambp = f0 * rr(ax-1,bx,1) * rr(ay,by,2) * rr(az,bz+1,3)
                   ELSE
                       ambp = 0.0_dp
                   END IF
                   IF(ax > 0 .AND. bz > 0) THEN
                       ambm = f0 * rr(ax-1,bx,1) * rr(ay,by,2) * rr(az,bz-1,3)
                   ELSE
                       ambm = 0.0_dp
                   END IF
                   ddab(ia,ib,3) = -4.0_dp*a*b*apbm + 2.0_dp*a*REAL(bz,dp)*apbm &
                      + 2.0_dp*b*REAL(ax,dp)*ambp - REAL(ax,dp)*REAL(bz,dp)*ambm
                   ! dy dy
                   apbp = f0 * rr(ax,bx,1) * rr(ay+1,by+1,2) * rr(az,bz,3)
                   IF(by > 0) THEN
                       apbm = f0 * rr(ax,bx,1) * rr(ay+1,by-1,2) * rr(az,bz,3)
                   ELSE
                       apbm = 0.0_dp
                   END IF
                   IF(ay > 0) THEN
                       ambp = f0 * rr(ax,bx,1) * rr(ay-1,by+1,2) * rr(az,bz,3)
                   ELSE
                       ambp = 0.0_dp
                   END IF
                   IF(ay > 0 .AND. by > 0) THEN
                       ambm = f0 * rr(ax,bx,1) * rr(ay-1,by-1,2) * rr(az,bz,3)
                   ELSE
                       ambm = 0.0_dp
                   END IF
                   ddab(ia,ib,4) = -4.0_dp*a*b*apbm + 2.0_dp*a*REAL(by,dp)*apbm &
                      + 2.0_dp*b*REAL(ay,dp)*ambp - REAL(ay,dp)*REAL(by,dp)*ambm
                   ! dy dz
                   apbp = f0 * rr(ax,bx,1) * rr(ay+1,by,2) * rr(az,bz+1,3)
                   IF(bz > 0) THEN
                       apbm = f0 * rr(ax,bx,1) * rr(ay+1,by,2) * rr(az,bz-1,3)
                   ELSE
                       apbm = 0.0_dp
                   END IF
                   IF(ay > 0) THEN
                       ambp = f0 * rr(ax,bx,1) * rr(ay-1,by,2) * rr(az,bz+1,3)
                   ELSE
                       ambp = 0.0_dp
                   END IF
                   IF(ay > 0 .AND. bz > 0) THEN
                       ambm = f0 * rr(ax,bx,1) * rr(ay-1,by,2) * rr(az,bz-1,3)
                   ELSE
                       ambm = 0.0_dp
                   END IF
                   ddab(ia,ib,5) = -4.0_dp*a*b*apbm + 2.0_dp*a*REAL(bz,dp)*apbm &
                      + 2.0_dp*b*REAL(ay,dp)*ambp - REAL(ay,dp)*REAL(bz,dp)*ambm
                   ! dz dz
                   apbp = f0 * rr(ax,bx,1) * rr(ay,by,2) * rr(az+1,bz+1,3)
                   IF(bz > 0) THEN
                       apbm = f0 * rr(ax,bx,1) * rr(ay,by,2) * rr(az+1,bz-1,3)
                   ELSE
                       apbm = 0.0_dp
                   END IF
                   IF(az > 0) THEN
                       ambp = f0 * rr(ax,bx,1) * rr(ay,by,2) * rr(az-1,bz+1,3)
                   ELSE
                       ambp = 0.0_dp
                   END IF
                   IF(az > 0 .AND. bz > 0) THEN
                       ambm = f0 * rr(ax,bx,1) * rr(ay,by,2) * rr(az-1,bz-1,3)
                   ELSE
                       ambm = 0.0_dp
                   END IF
                   ddab(ia,ib,6) = -4.0_dp*a*b*apbm + 2.0_dp*a*REAL(bz,dp)*apbm &
                      + 2.0_dp*b*REAL(az,dp)*ambp - REAL(az,dp)*REAL(bz,dp)*ambm
                END IF
                !
             ENDDO
             ENDDO
             ENDDO !la
          ENDDO
          ENDDO
          ENDDO !lb

          mb = mb + nb
       END DO
       ma = ma + na
    END DO

    DEALLOCATE(rr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_ab

! *****************************************************************************
!> \brief   Calculation of the two-center overlap integrals [aa|b] over
!>          Cartesian Gaussian-type functions.
!> \param la1_max    Max L on center A (basis 1)
!> \param la1_min    Min L on center A (basis 1)
!> \param npgfa1     Number of primitives on center A (basis 1)
!> \param rpgfa1     Range of functions on A, used for screening (basis 1)
!> \param zeta1      Exponents on center A (basis 1)
!> \param la2_max    Max L on center A (basis 2)
!> \param la2_min    Min L on center A (basis 2)
!> \param npgfa2     Number of primitives on center A (basis 2)
!> \param rpgfa2     Range of functions on A, used for screening (basis 2)
!> \param zeta2      Exponents on center A (basis 2)
!> \param lb_max     Max L on center B
!> \param lb_min     Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param saab       Final overlap integrals
!> \param daab       First derivative overlap integrals
!> \param error      CP2K error reporting
!> \date    01.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE overlap_aab(la1_max,la1_min,npgfa1,rpgfa1,zeta1,&
                         la2_max,la2_min,npgfa2,rpgfa2,zeta2,&
                         lb_max,lb_min,npgfb,rpgfb,zetb,&
                         rab,saab,daab,error)
    INTEGER, INTENT(IN)                      :: la1_max, la1_min, npgfa1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa1, zeta1
    INTEGER, INTENT(IN)                      :: la2_max, la2_min, npgfa2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa2, zeta2
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: saab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: daab
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_aab', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ax1, ax2, ay, ay1, ay2, az, az1, az2, bx, by, bz, coa1, &
      coa2, cob, i1pgf, i2pgf, ia1, ia2, ib, istat, jpgf, la1, la2, lb, ldrr, &
      lma, lmb, ma1, ma2, mb, na1, na2, nb, ofa1, ofa2, ofb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, dumx, dumy, dumz, f0, &
                                                rab2, rpgfa, tab, xhi, zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(saab)) THEN
       lma = la1_max + la2_max
       lmb = lb_max
    END IF
    IF (PRESENT(daab)) THEN
       lma = la1_max + la2_max + 1
       lmb = lb_max
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa1 = ncoset(la1_min-1)
    ofa2 = ncoset(la2_min-1)
    ofb = ncoset(lb_min-1)
    na1 = ncoset(la1_max) - ofa1
    na2 = ncoset(la2_max) - ofa2
    nb = ncoset(lb_max) - ofb
    IF (PRESENT(saab)) THEN
       CPPrecondition((SIZE(saab,1) >= na1*npgfa1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saab,2) >= na2*npgfa2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saab,3) >= nb*npgfb),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(daab)) THEN
       CPPrecondition((SIZE(daab,1) >= na1*npgfa1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daab,2) >= na2*npgfa2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daab,3) >= nb*npgfb),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daab,4) >= 3),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all primitive Gaussian-type functions
    ma1 = 0
    DO i1pgf=1,npgfa1
       ma2 = 0
       DO i2pgf=1,npgfa2
          rpgfa = MIN(rpgfa1(i1pgf),rpgfa2(i2pgf))
          mb = 0
          DO jpgf=1,npgfb
             ! Distance Screening
             IF(rpgfa + rpgfb(jpgf) < tab) THEN
                IF (PRESENT(saab)) saab(ma1+1:ma1+na1,ma2+1:ma2+na2,mb+1:mb+nb) = 0.0_dp
                IF (PRESENT(daab)) daab(ma1+1:ma1+na1,ma2+1:ma2+na2,mb+1:mb+nb,1:3) = 0.0_dp
                mb = mb + nb
                CYCLE
             ENDIF

             ! Calculate some prefactors
             a = zeta1(i1pgf) + zeta2(i2pgf)
             b = zetb(jpgf)
             zet =  a+b
             xhi =  a*b/zet
             rap =  b*rab/zet
             rbp = -a*rab/zet

             ! [ss|s] integral
             f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

             ! Calculate the reccurence relation
             CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

             DO lb = lb_min,lb_max
             DO bx = 0,lb
             DO by = 0,lb-bx
                bz = lb - bx - by
                cob = coset(bx,by,bz) - ofb
                ib = mb + cob
                DO la2 = la2_min,la2_max
                DO ax2 = 0,la2
                DO ay2 = 0,la2-ax2
                   az2 = la2 - ax2 - ay2
                   coa2 = coset(ax2,ay2,az2) - ofa2
                   ia2 = ma2 + coa2
                   DO la1 = la1_min,la1_max
                   DO ax1 = 0,la1
                   DO ay1 = 0,la1-ax1
                      az1 = la1 - ax1 - ay1
                      coa1 = coset(ax1,ay1,az1) - ofa1
                      ia1 = ma1 + coa1
                      ! integrals
                      IF (PRESENT(saab)) THEN
                         saab(ia1,ia2,ib) = f0 * rr(ax1+ax2,bx,1) * rr(ay1+ay2,by,2) * rr(az1+az2,bz,3)
                      END IF
                      ! first derivatives
                      IF (PRESENT(daab)) THEN
                         ax = ax1 + ax2
                         ay = ay1 + ay2
                         az = az1 + az2
                         ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                         ! dx
                         dumx = 2.0_dp * a * rr(ax+1,bx,1)
                         IF(ax > 0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                         daab(ia1,ia2,ib,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                         ! dy
                         dumy = 2.0_dp * a * rr(ay+1,by,2)
                         IF(ay > 0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                         daab(ia1,ia2,ib,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                         ! dz
                         dumz = 2.0_dp * a * rr(az+1,bz,3)
                         IF(az > 0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                         daab(ia1,ia2,ib,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                      END IF
                      !
                   ENDDO
                   ENDDO
                   ENDDO !la1
                ENDDO
                ENDDO
                ENDDO !la2
             ENDDO
             ENDDO
             ENDDO !lb

             mb = mb + nb
          END DO
          ma2 = ma2 + na2
       END DO
       ma1 = ma1 + na1
    END DO

    DEALLOCATE(rr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_aab

! *****************************************************************************
!> \brief   Calculation of the two-center overlap integrals [a|bb] over
!>          Cartesian Gaussian-type functions. 
!> \param la_max     Max L on center A
!> \param la_min     Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb1_max    Max L on center B (basis 1)
!> \param lb1_min    Min L on center B (basis 1)
!> \param npgfb1     Number of primitives on center B (basis 1)
!> \param rpgfb1     Range of functions on B, used for screening (basis 1)
!> \param zetb1      Exponents on center B (basis 1)
!> \param lb2_max    Max L on center B (basis 2)
!> \param lb2_min    Min L on center B (basis 2)
!> \param npgfb2     Number of primitives on center B (basis 2)
!> \param rpgfb2     Range of functions on B, used for screening (basis 2)
!> \param zetb2      Exponents on center B (basis 2)
!> \param rab        Distance vector A-B
!> \param sabb       Final overlap integrals
!> \param dabb       First derivative overlap integrals
!> \param error      CP2K error reporting
!> \date    01.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE overlap_abb(la_max,la_min,npgfa,rpgfa,zeta,&
                         lb1_max,lb1_min,npgfb1,rpgfb1,zetb1,&
                         lb2_max,lb2_min,npgfb2,rpgfb2,zetb2,&
                         rab,sabb,dabb,error)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb1_max, lb1_min, npgfb1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb1, zetb1
    INTEGER, INTENT(IN)                      :: lb2_max, lb2_min, npgfb2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb2, zetb2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: sabb
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: dabb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_abb', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, bx1, bx2, by, by1, by2, bz, bz1, bz2, coa, &
      cob1, cob2, ia, ib1, ib2, ipgf, istat, j1pgf, j2pgf, la, lb1, lb2, &
      ldrr, lma, lmb, ma, mb1, mb2, na, nb1, nb2, ofa, ofb1, ofb2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, dumx, dumy, dumz, f0, &
                                                rab2, rpgfb, tab, xhi, zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(sabb)) THEN
       lma = la_max
       lmb = lb1_max + lb2_max
    END IF
    IF (PRESENT(dabb)) THEN
       lma = la_max + 1
       lmb = lb1_max + lb2_max
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa = ncoset(la_min-1)
    ofb1 = ncoset(lb1_min-1)
    ofb2 = ncoset(lb2_min-1)
    na = ncoset(la_max) - ofa
    nb1 = ncoset(lb1_max) - ofb1
    nb2 = ncoset(lb2_max) - ofb2
    IF (PRESENT(sabb)) THEN
       CPPrecondition((SIZE(sabb,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(sabb,2) >= nb1*npgfb1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(sabb,3) >= nb2*npgfb2),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(dabb)) THEN
       CPPrecondition((SIZE(dabb,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabb,2) >= nb1*npgfb1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabb,3) >= nb2*npgfb2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabb,4) >= 3),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all pairs of primitive Gaussian-type functions
    ma = 0
    DO ipgf=1,npgfa
       mb1 = 0
       DO j1pgf=1,npgfb1
          mb2 = 0
          DO j2pgf=1,npgfb2
             ! Distance Screening
             rpgfb = MIN(rpgfb1(j1pgf),rpgfb2(j2pgf))
             IF(rpgfa(ipgf) + rpgfb < tab) THEN
                IF (PRESENT(sabb)) sabb(ma+1:ma+na,mb1+1:mb1+nb1,mb2+1:mb2+nb2) = 0.0_dp
                IF (PRESENT(dabb)) dabb(ma+1:ma+na,mb1+1:mb1+nb1,mb2+1:mb2+nb2,1:3) = 0.0_dp
                mb2 = mb2 + nb2
                CYCLE
             ENDIF

             ! Calculate some prefactors
             a = zeta(ipgf)
             b = zetb1(j1pgf) + zetb2(j2pgf)
             zet =  a+b
             xhi =  a*b/zet
             rap =  b*rab/zet
             rbp = -a*rab/zet

             ! [s|s] integral
             f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

             ! Calculate the reccurence relation
             CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

          DO lb2 = lb2_min,lb2_max
          DO bx2 = 0,lb2
          DO by2 = 0,lb2-bx2
             bz2 = lb2 - bx2 - by2
             cob2 = coset(bx2,by2,bz2) - ofb2
             ib2 = mb2 + cob2
             DO lb1 = lb1_min,lb1_max
             DO bx1 = 0,lb1
             DO by1 = 0,lb1-bx1
                bz1 = lb1 - bx1 - by1
                cob1 = coset(bx1,by1,bz1) - ofb1
                ib1 = mb1 + cob1
                DO la = la_min,la_max
                DO ax = 0,la
                DO ay = 0,la-ax
                   az = la - ax - ay
                   coa = coset(ax,ay,az) - ofa
                   ia = ma + coa
                   ! integrals
                   IF (PRESENT(sabb)) THEN
                      sabb(ia,ib1,ib2) = f0 * rr(ax,bx1+bx2,1) * rr(ay,by1+by2,2) * rr(az,bz1+bz2,3)
                   END IF
                   ! first derivatives
                   IF (PRESENT(dabb)) THEN
                      bx = bx1 + bx2
                      by = by1 + by2
                      bz = bz1 + bz2
                      ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                      ! dx
                      dumx = 2.0_dp * a * rr(ax+1,bx,1)
                      IF(ax > 0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                      dabb(ia,ib1,ib2,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                      ! dy
                      dumy = 2.0_dp * a * rr(ay+1,by,2)
                      IF(ay > 0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                      dabb(ia,ib1,ib2,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                      ! dz
                      dumz = 2.0_dp * a * rr(az+1,bz,3)
                      IF(az > 0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                      dabb(ia,ib1,ib2,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                   END IF
                   !
                   ENDDO
                   ENDDO
                   ENDDO !la
                ENDDO
                ENDDO
                ENDDO !lb1
             ENDDO
             ENDDO
             ENDDO !lb2

             mb2 = mb2 + nb2
          END DO
          mb1 = mb1 + nb1
       END DO
       ma = ma + na
    END DO

    DEALLOCATE(rr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_abb

! *****************************************************************************
!> \brief   Calculation of the two-center overlap integrals [aaa|b] over
!>          Cartesian Gaussian-type functions.
!> \param la1_max    Max L on center A (basis 1)
!> \param la1_min    Min L on center A (basis 1)
!> \param npgfa1     Number of primitives on center A (basis 1)
!> \param rpgfa1     Range of functions on A, used for screening (basis 1)
!> \param zeta1      Exponents on center A (basis 1)
!> \param la2_max    Max L on center A (basis 2)
!> \param la2_min    Min L on center A (basis 2)
!> \param npgfa2     Number of primitives on center A (basis 2)
!> \param rpgfa2     Range of functions on A, used for screening (basis 2)
!> \param zeta2      Exponents on center A (basis 2)
!> \param la3_max    Max L on center A (basis 3)
!> \param la3_min    Min L on center A (basis 3)
!> \param npgfa3     Number of primitives on center A (basis 3)
!> \param rpgfa3     Range of functions on A, used for screening (basis 3)
!> \param zeta3      Exponents on center A (basis 3)
!> \param lb_max     Max L on center B
!> \param lb_min     Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param saaab      Final overlap integrals
!> \param daaab      First derivative overlap integrals
!> \param error      CP2K error reporting
!> \date    01.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE overlap_aaab(la1_max,la1_min,npgfa1,rpgfa1,zeta1,&
                          la2_max,la2_min,npgfa2,rpgfa2,zeta2,&
                          la3_max,la3_min,npgfa3,rpgfa3,zeta3,&
                          lb_max,lb_min,npgfb,rpgfb,zetb,&
                          rab,saaab,daaab,error)
    INTEGER, INTENT(IN)                      :: la1_max, la1_min, npgfa1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa1, zeta1
    INTEGER, INTENT(IN)                      :: la2_max, la2_min, npgfa2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa2, zeta2
    INTEGER, INTENT(IN)                      :: la3_max, la3_min, npgfa3
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa3, zeta3
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: saaab
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: daaab
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_aaab', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ax1, ax2, ax3, ay, ay1, ay2, ay3, az, az1, az2, az3, bx, &
      by, bz, coa1, coa2, coa3, cob, i1pgf, i2pgf, i3pgf, ia1, ia2, ia3, ib, &
      istat, jpgf, la1, la2, la3, lb, ldrr, lma, lmb, ma1, ma2, ma3, mb, na1, &
      na2, na3, nb, ofa1, ofa2, ofa3, ofb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, dumx, dumy, dumz, f0, &
                                                rab2, rpgfa, tab, xhi, zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(saaab)) THEN
       lma = la1_max + la2_max + la3_max
       lmb = lb_max
    END IF
    IF (PRESENT(daaab)) THEN
       lma = la1_max + la2_max + la3_max + 1
       lmb = lb_max
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa1 = ncoset(la1_min-1)
    ofa2 = ncoset(la2_min-1)
    ofa3 = ncoset(la3_min-1)
    ofb = ncoset(lb_min-1)
    na1 = ncoset(la1_max) - ofa1
    na2 = ncoset(la2_max) - ofa2
    na3 = ncoset(la3_max) - ofa3
    nb = ncoset(lb_max) - ofb
    IF (PRESENT(saaab)) THEN
       CPPrecondition((SIZE(saaab,1) >= na1*npgfa1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saaab,2) >= na2*npgfa2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saaab,3) >= na3*npgfa3),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saaab,4) >= nb*npgfb),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(daaab)) THEN
       CPPrecondition((SIZE(daaab,1) >= na1*npgfa1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daaab,2) >= na2*npgfa2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daaab,3) >= na3*npgfa3),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daaab,4) >= nb*npgfb),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daaab,5) >= 3),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all primitive Gaussian-type functions
    ma1 = 0
    DO i1pgf=1,npgfa1
       ma2 = 0
       DO i2pgf=1,npgfa2
          ma3 = 0
          DO i3pgf=1,npgfa3
             rpgfa = MIN(rpgfa1(i1pgf),rpgfa2(i2pgf),rpgfa3(i3pgf))
             mb = 0
             DO jpgf=1,npgfb
                ! Distance Screening
                IF(rpgfa + rpgfb(jpgf) < tab) THEN
                   IF (PRESENT(saaab)) saaab(ma1+1:ma1+na1,ma2+1:ma2+na2,ma3+1:ma3+na3,mb+1:mb+nb) = 0.0_dp
                   IF (PRESENT(daaab)) daaab(ma1+1:ma1+na1,ma2+1:ma2+na2,ma3+1:ma3+na3,mb+1:mb+nb,1:3) = 0.0_dp
                   mb = mb + nb
                   CYCLE
                ENDIF

                ! Calculate some prefactors
                a = zeta1(i1pgf) + zeta2(i2pgf) + zeta3(i3pgf)
                b = zetb(jpgf)
                zet =  a+b
                xhi =  a*b/zet
                rap =  b*rab/zet
                rbp = -a*rab/zet

                ! [sss|s] integral
                f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

                ! Calculate the reccurence relation
                CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

                DO lb = lb_min,lb_max
                DO bx = 0,lb
                DO by = 0,lb-bx
                   bz = lb - bx - by
                   cob = coset(bx,by,bz) - ofb
                   ib = mb + cob
                   DO la3 = la3_min,la3_max
                   DO ax3 = 0,la3
                   DO ay3 = 0,la3-ax3
                      az3 = la3 - ax3 - ay3
                      coa3 = coset(ax3,ay3,az3) - ofa3
                      ia3 = ma3 + coa3
                      DO la2 = la2_min,la2_max
                      DO ax2 = 0,la2
                      DO ay2 = 0,la2-ax2
                         az2 = la2 - ax2 - ay2
                         coa2 = coset(ax2,ay2,az2) - ofa2
                         ia2 = ma2 + coa2
                         DO la1 = la1_min,la1_max
                         DO ax1 = 0,la1
                         DO ay1 = 0,la1-ax1
                            az1 = la1 - ax1 - ay1
                            coa1 = coset(ax1,ay1,az1) - ofa1
                            ia1 = ma1 + coa1
                            ! integrals
                            IF (PRESENT(saaab)) THEN
                               saaab(ia1,ia2,ia3,ib) = f0 * rr(ax1+ax2+ax3,bx,1) * &
                                                            rr(ay1+ay2+ay3,by,2) * rr(az1+az2+az3,bz,3)
                            END IF
                            ! first derivatives
                            IF (PRESENT(daaab)) THEN
                               ax = ax1 + ax2 + ax3
                               ay = ay1 + ay2 + ay3
                               az = az1 + az2 + az3
                               ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                               ! dx
                               dumx = 2.0_dp * a * rr(ax+1,bx,1)
                               IF(ax > 0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                               daaab(ia1,ia2,ia3,ib,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                               ! dy
                               dumy = 2.0_dp * a * rr(ay+1,by,2)
                               IF(ay > 0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                               daaab(ia1,ia2,ia3,ib,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                               ! dz
                               dumz = 2.0_dp * a * rr(az+1,bz,3)
                               IF(az > 0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                               daaab(ia1,ia2,ia3,ib,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                            END IF
                            !
                         ENDDO
                         ENDDO
                         ENDDO !la1
                      ENDDO
                      ENDDO
                      ENDDO !la2
                   ENDDO
                   ENDDO
                   ENDDO !la3
                ENDDO
                ENDDO
                ENDDO !lb

                mb = mb + nb
             END DO
             ma3 = ma3 + na3
          END DO
          ma2 = ma2 + na2
       END DO
       ma1 = ma1 + na1
    END DO

    DEALLOCATE(rr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_aaab
! *****************************************************************************
!> \brief   Calculation of the two-center overlap integrals [aa|bb] over
!>          Cartesian Gaussian-type functions.
!> \param la1_max    Max L on center A (basis 1)
!> \param la1_min    Min L on center A (basis 1)
!> \param npgfa1     Number of primitives on center A (basis 1)
!> \param rpgfa1     Range of functions on A, used for screening (basis 1)
!> \param zeta1      Exponents on center A (basis 1)
!> \param la2_max    Max L on center A (basis 2)
!> \param la2_min    Min L on center A (basis 2)
!> \param npgfa2     Number of primitives on center A (basis 2)
!> \param rpgfa2     Range of functions on A, used for screening (basis 2)
!> \param zeta2      Exponents on center A (basis 2)
!> \param lb1_max    Max L on center B (basis 3)
!> \param lb1_min    Min L on center B (basis 3)
!> \param npgfb1     Number of primitives on center B (basis 3)
!> \param rpgfb1     Range of functions on B, used for screening (basis 3)
!> \param zetb1      Exponents on center B (basis 3)
!> \param lb2_max    Max L on center B (basis 4)
!> \param lb2_min    Min L on center B (basis 4)
!> \param npgfb2     Number of primitives on center B (basis 4)
!> \param rpgfb2     Range of functions on B, used for screening (basis 4)
!> \param zetb2      Exponents on center B (basis 4)
!> \param rab        Distance vector A-B
!> \param saabb      Final overlap integrals
!> \param daabb      First derivative overlap integrals
!> \param error      CP2K error reporting
!> \date    01.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE overlap_aabb(la1_max,la1_min,npgfa1,rpgfa1,zeta1,&
                          la2_max,la2_min,npgfa2,rpgfa2,zeta2,&
                          lb1_max,lb1_min,npgfb1,rpgfb1,zetb1,&
                          lb2_max,lb2_min,npgfb2,rpgfb2,zetb2,&
                          rab,saabb,daabb,error)
    INTEGER, INTENT(IN)                      :: la1_max, la1_min, npgfa1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa1, zeta1
    INTEGER, INTENT(IN)                      :: la2_max, la2_min, npgfa2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa2, zeta2
    INTEGER, INTENT(IN)                      :: lb1_max, lb1_min, npgfb1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb1, zetb1
    INTEGER, INTENT(IN)                      :: lb2_max, lb2_min, npgfb2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb2, zetb2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: saabb
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: daabb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_aabb', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ax1, ax2, ay, ay1, ay2, az, az1, az2, bx, bx1, bx2, by, &
      by1, by2, bz, bz1, bz2, coa1, coa2, cob1, cob2, i1pgf, i2pgf, ia1, ia2, &
      ib1, ib2, istat, j1pgf, j2pgf, la1, la2, lb1, lb2, ldrr, lma, lmb, ma1, &
      ma2, mb1, mb2, na1, na2, nb1, nb2, ofa1, ofa2, ofb1, ofb2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, dumx, dumy, dumz, f0, &
                                                rab2, rpgfa, rpgfb, tab, xhi, &
                                                zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(saabb)) THEN
       lma = la1_max + la2_max
       lmb = lb1_max + lb2_max
    END IF
    IF (PRESENT(daabb)) THEN
       lma = la1_max + la2_max + 1
       lmb = lb1_max + lb2_max
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa1 = ncoset(la1_min-1)
    ofa2 = ncoset(la2_min-1)
    ofb1 = ncoset(lb1_min-1)
    ofb2 = ncoset(lb2_min-1)
    na1 = ncoset(la1_max) - ofa1
    na2 = ncoset(la2_max) - ofa2
    nb1 = ncoset(lb1_max) - ofb1
    nb2 = ncoset(lb2_max) - ofb2
    IF (PRESENT(saabb)) THEN
       CPPrecondition((SIZE(saabb,1) >= na1*npgfa1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saabb,2) >= na2*npgfa2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saabb,3) >= nb1*npgfb1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(saabb,4) >= nb2*npgfb2),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(daabb)) THEN
       CPPrecondition((SIZE(daabb,1) >= na1*npgfa1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daabb,2) >= na2*npgfa2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daabb,3) >= nb1*npgfb1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daabb,4) >= nb2*npgfb2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(daabb,5) >= 3),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all primitive Gaussian-type functions
    ma1 = 0
    DO i1pgf=1,npgfa1
       ma2 = 0
       DO i2pgf=1,npgfa2
          mb1 = 0
          rpgfa = MIN(rpgfa1(i1pgf),rpgfa2(i2pgf))
          DO j1pgf=1,npgfb1
             mb2 = 0
             DO j2pgf=1,npgfb2
                rpgfb = MIN(rpgfb1(j1pgf),rpgfb2(j2pgf))
                ! Distance Screening
                IF(rpgfa + rpgfb < tab) THEN
                   IF (PRESENT(saabb)) saabb(ma1+1:ma1+na1,ma2+1:ma2+na2,mb1+1:mb1+nb1,mb2+1:mb2+nb2) = 0.0_dp
                   IF (PRESENT(daabb)) daabb(ma1+1:ma1+na1,ma2+1:ma2+na2,mb1+1:mb1+nb1,mb2+1:mb2+nb2,1:3) = 0.0_dp
                   mb2 = mb2 + nb2
                   CYCLE
                ENDIF

                ! Calculate some prefactors
                a = zeta1(i1pgf) + zeta2(i2pgf)
                b = zetb1(j1pgf) + zetb2(j2pgf)
                zet =  a+b
                xhi =  a*b/zet
                rap =  b*rab/zet
                rbp = -a*rab/zet

                ! [ss|ss] integral
                f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

                ! Calculate the reccurence relation
                CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

                DO lb2 = lb2_min,lb2_max
                DO bx2 = 0,lb2
                DO by2 = 0,lb2-bx2
                   bz2 = lb2 - bx2 - by2
                   cob2 = coset(bx2,by2,bz2) - ofb2
                   ib2 = mb2 + cob2
                   DO lb1 = lb1_min,lb1_max
                   DO bx1 = 0,lb1
                   DO by1 = 0,lb1-bx1
                      bz1 = lb1 - bx1 - by1
                      cob1 = coset(bx1,by1,bz1) - ofb1
                      ib1 = mb1 + cob1
                      DO la2 = la2_min,la2_max
                      DO ax2 = 0,la2
                      DO ay2 = 0,la2-ax2
                         az2 = la2 - ax2 - ay2
                         coa2 = coset(ax2,ay2,az2) - ofa2
                         ia2 = ma2 + coa2
                         DO la1 = la1_min,la1_max
                         DO ax1 = 0,la1
                         DO ay1 = 0,la1-ax1
                            az1 = la1 - ax1 - ay1
                            coa1 = coset(ax1,ay1,az1) - ofa1
                            ia1 = ma1 + coa1
                            ! integrals
                            IF (PRESENT(saabb)) THEN
                               saabb(ia1,ia2,ib1,ib2) = f0 * rr(ax1+ax2,bx1+bx2,1) * &
                                                       rr(ay1+ay2,by1+by2,2) * rr(az1+az2,bz1+bz2,3)
                            END IF
                            ! first derivatives
                            IF (PRESENT(daabb)) THEN
                               ax = ax1 + ax2
                               ay = ay1 + ay2
                               az = az1 + az2
                               bx = bx1 + bx2
                               by = by1 + by2
                               bz = bz1 + bz2
                               ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                               ! dx
                               dumx = 2.0_dp * a * rr(ax+1,bx,1)
                               IF(ax > 0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                               daabb(ia1,ia2,ib1,ib2,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                               ! dy
                               dumy = 2.0_dp * a * rr(ay+1,by,2)
                               IF(ay > 0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                               daabb(ia1,ia2,ib1,ib2,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                               ! dz
                               dumz = 2.0_dp * a * rr(az+1,bz,3)
                               IF(az > 0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                               daabb(ia1,ia2,ib1,ib2,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                            END IF
                            !
                         ENDDO
                         ENDDO
                         ENDDO !la1
                      ENDDO
                      ENDDO
                      ENDDO !la2
                   ENDDO
                   ENDDO
                   ENDDO !lb1
                ENDDO
                ENDDO
                ENDDO !lb2

                mb2 = mb2 + nb2
             END DO
             mb1 = mb1 + nb1
          END DO
          ma2 = ma2 + na2
       END DO
       ma1 = ma1 + na1
    END DO

    DEALLOCATE(rr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_aabb
! *****************************************************************************
!> \brief   Calculation of the two-center overlap integrals [a|bbb] over
!>          Cartesian Gaussian-type functions. 
!> \param la_max     Max L on center A
!> \param la_min     Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb1_max    Max L on center B (basis 1)
!> \param lb1_min    Min L on center B (basis 1)
!> \param npgfb1     Number of primitives on center B (basis 1)
!> \param rpgfb1     Range of functions on B, used for screening (basis 1)
!> \param zetb1      Exponents on center B (basis 1)
!> \param lb2_max    Max L on center B (basis 2)
!> \param lb2_min    Min L on center B (basis 2)
!> \param npgfb2     Number of primitives on center B (basis 2)
!> \param rpgfb2     Range of functions on B, used for screening (basis 2)
!> \param zetb2      Exponents on center B (basis 2)
!> \param lb3_max    Max L on center B (basis 3)
!> \param lb3_min    Min L on center B (basis 3)
!> \param npgfb3     Number of primitives on center B (basis 3)
!> \param rpgfb3     Range of functions on B, used for screening (basis 3)
!> \param zetb3      Exponents on center B (basis 3)
!> \param rab        Distance vector A-B
!> \param sabbb      Final overlap integrals
!> \param dabbb      First derivative overlap integrals
!> \param error      CP2K error reporting
!> \date    01.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE overlap_abbb(la_max,la_min,npgfa,rpgfa,zeta,&
                          lb1_max,lb1_min,npgfb1,rpgfb1,zetb1,&
                          lb2_max,lb2_min,npgfb2,rpgfb2,zetb2,&
                          lb3_max,lb3_min,npgfb3,rpgfb3,zetb3,&
                          rab,sabbb,dabbb,error)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb1_max, lb1_min, npgfb1
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb1, zetb1
    INTEGER, INTENT(IN)                      :: lb2_max, lb2_min, npgfb2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb2, zetb2
    INTEGER, INTENT(IN)                      :: lb3_max, lb3_min, npgfb3
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb3, zetb3
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: sabbb
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), &
      INTENT(INOUT), OPTIONAL                :: dabbb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_abbb', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, bx1, bx2, bx3, by, by1, by2, by3, bz, bz1, &
      bz2, bz3, coa, cob1, cob2, cob3, ia, ib1, ib2, ib3, ipgf, istat, j1pgf, &
      j2pgf, j3pgf, la, lb1, lb2, lb3, ldrr, lma, lmb, ma, mb1, mb2, mb3, na, &
      nb1, nb2, nb3, ofa, ofb1, ofb2, ofb3
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, dumx, dumy, dumz, f0, &
                                                rab2, rpgfb, tab, xhi, zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(sabbb)) THEN
       lma = la_max
       lmb = lb1_max + lb2_max + lb3_max
    END IF
    IF (PRESENT(dabbb)) THEN
       lma = la_max + 1
       lmb = lb1_max + lb2_max + lb3_max
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa = ncoset(la_min-1)
    ofb1 = ncoset(lb1_min-1)
    ofb2 = ncoset(lb2_min-1)
    ofb3 = ncoset(lb3_min-1)
    na = ncoset(la_max) - ofa
    nb1 = ncoset(lb1_max) - ofb1
    nb2 = ncoset(lb2_max) - ofb2
    nb3 = ncoset(lb3_max) - ofb3
    IF (PRESENT(sabbb)) THEN
       CPPrecondition((SIZE(sabbb,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(sabbb,2) >= nb1*npgfb1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(sabbb,3) >= nb2*npgfb2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(sabbb,4) >= nb3*npgfb3),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(dabbb)) THEN
       CPPrecondition((SIZE(dabbb,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabbb,2) >= nb1*npgfb1),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabbb,3) >= nb2*npgfb2),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabbb,4) >= nb3*npgfb3),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dabbb,5) >= 3),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all pairs of primitive Gaussian-type functions
    ma = 0
    DO ipgf=1,npgfa
       mb1 = 0
       DO j1pgf=1,npgfb1
          mb2 = 0
          DO j2pgf=1,npgfb2
             mb3 = 0
             DO j3pgf=1,npgfb3
                ! Distance Screening
                rpgfb = MIN(rpgfb1(j1pgf),rpgfb2(j2pgf),rpgfb3(j3pgf))
                IF(rpgfa(ipgf) + rpgfb < tab) THEN
                   IF (PRESENT(sabbb)) sabbb(ma+1:ma+na,mb1+1:mb1+nb1,mb2+1:mb2+nb2,mb3+1:mb3+nb3) = 0.0_dp
                   IF (PRESENT(dabbb)) dabbb(ma+1:ma+na,mb1+1:mb1+nb1,mb2+1:mb2+nb2,mb3+1:mb3+nb3,1:3) = 0.0_dp
                   mb3 = mb3 + nb3
                   CYCLE
                ENDIF

                ! Calculate some prefactors
                a = zeta(ipgf)
                b = zetb1(j1pgf) + zetb2(j2pgf) + zetb3(j3pgf)
                zet =  a+b
                xhi =  a*b/zet
                rap =  b*rab/zet
                rbp = -a*rab/zet

                ! [s|s] integral
                f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

                ! Calculate the reccurence relation
                CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

                DO lb3 = lb3_min,lb3_max
                DO bx3 = 0,lb3
                DO by3 = 0,lb3-bx3
                   bz3 = lb3 - bx3 - by3
                   cob3 = coset(bx3,by3,bz3) - ofb3
                   ib3 = mb3 + cob3
                DO lb2 = lb2_min,lb2_max
                DO bx2 = 0,lb2
                DO by2 = 0,lb2-bx2
                   bz2 = lb2 - bx2 - by2
                   cob2 = coset(bx2,by2,bz2) - ofb2
                   ib2 = mb2 + cob2
                   DO lb1 = lb1_min,lb1_max
                   DO bx1 = 0,lb1
                   DO by1 = 0,lb1-bx1
                      bz1 = lb1 - bx1 - by1
                      cob1 = coset(bx1,by1,bz1) - ofb1
                      ib1 = mb1 + cob1
                      DO la = la_min,la_max
                      DO ax = 0,la
                      DO ay = 0,la-ax
                         az = la - ax - ay
                         coa = coset(ax,ay,az) - ofa
                         ia = ma + coa
                         ! integrals
                         IF (PRESENT(sabbb)) THEN
                            sabbb(ia,ib1,ib2,ib3) = f0 * rr(ax,bx1+bx2+bx3,1) * &
                                                    rr(ay,by1+by2+by3,2) * rr(az,bz1+bz2+bz3,3)
                         END IF
                         ! first derivatives
                         IF (PRESENT(dabbb)) THEN
                            bx = bx1 + bx2 + bx3
                            by = by1 + by2 + by3
                            bz = bz1 + bz2 + bz3
                            ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                            ! dx
                            dumx = 2.0_dp * a * rr(ax+1,bx,1)
                            IF(ax > 0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                            dabbb(ia,ib1,ib2,ib3,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                            ! dy
                            dumy = 2.0_dp * a * rr(ay+1,by,2)
                            IF(ay > 0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                            dabbb(ia,ib1,ib2,ib3,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                            ! dz
                            dumz = 2.0_dp * a * rr(az+1,bz,3)
                            IF(az > 0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                            dabbb(ia,ib1,ib2,ib3,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                         END IF
                         !
                         ENDDO
                         ENDDO
                         ENDDO !la
                      ENDDO
                      ENDDO
                      ENDDO !lb1
                   ENDDO
                   ENDDO
                   ENDDO !lb2
                ENDDO
                ENDDO
                ENDDO !lb3

                mb3 = mb3 + nb3
             END DO
             mb2 = mb2 + nb2
          END DO
          mb1 = mb1 + nb1
       END DO
       ma = ma + na
    END DO

    DEALLOCATE(rr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap_abbb

END MODULE ai_overlap
