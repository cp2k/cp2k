!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_overlap_ppl [1.0] *
!!
!!   NAME
!!     ai_overlap_ppl
!!
!!   FUNCTION
!!     Calculation of three-center overlap integrals over Cartesian
!!     Gaussian-type functions for the second term V(ppl) of the local
!!     part of the Goedecker pseudopotential (GTH):
!!
!!     <a|V(local)|b> = <a|V(erf) + V(ppl)|b>
!!                    = <a|V(erf)|b> + <a|V(ppl)|b>
!!                    = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!!                      (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!!                       C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!!
!!   AUTHOR
!!     Matthias Krack (04.10.2000)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!     S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!!     C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!!
!!   MODIFICATION HISTORY
!!     - Derivatives added (17.05.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE ai_overlap_ppl

! *****************************************************************************

! ax,ay,az   : Angular momentum index numbers of orbital a.
! bx,by,bz   : Angular momentum index numbers of orbital b.
! coset      : Cartesian orbital set pointer.
! dab        : Distance between the atomic centers a and b.
! dac        : Distance between the atomic centers a and c.
! dbc        : Distance between the atomic centers b and c.
! l{a,b,c}   : Angular momentum quantum number of shell a, b or c.
! l{a,b}_max : Maximum angular momentum quantum number of shell a, b or c.
! ncoset     : Number of Cartesian orbitals up to l.
! rab        : Distance vector between the atomic centers a and b.
! rac        : Distance vector between the atomic centers a and c.
! rbc        : Distance vector between the atomic centers b and c.
! rpgf{a,b,c}: Radius of the primitive Gaussian-type function a or b.
! zet{a,b,c} : Exponents of the Gaussian-type functions a or b.
! zetg       : Reciprocal of the sum of the exponents of orbital a, b and c.
! zetp       : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "ai_overlap_ppl"

! *** Public subroutines ***

  PUBLIC :: overlap_ppl

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE overlap_ppl(la_max_set,la_min_set,npgfa,rpgfa,zeta,&
                         lb_max_set,lb_min_set,npgfb,rpgfb,zetb,&
                         cexp_ppl,zetc,rpgfc,&
                         rab,dab,rac,dac,rbc,dbc,&
                         vab,da_max,db_max,return_derivatives,s,&
                         pab,force_a,force_b)

!   Purpose: Calculation of three-center overlap integrals <a|c|b> over
!            Cartesian Gaussian functions for the local part of the Goedecker
!            pseudopotential (GTH). c is a primitive Gaussian-type function
!            with a set of even angular momentum indices.

!            <a|V(ppl)|b> = <a| (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!                                C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>

!   History: - Creation (04.10.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: la_max_set, la_min_set, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max_set, lb_min_set, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb, cexp_ppl
    REAL(KIND=dp), INTENT(IN)                :: zetc, rpgfc
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: dac
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rbc
    REAL(KIND=dp), INTENT(IN)                :: dbc
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vab
    INTEGER, INTENT(IN)                      :: da_max, db_max
    LOGICAL, INTENT(IN)                      :: return_derivatives
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT)                          :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pab
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: force_a, force_b

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "overlap_ppl"

    INTEGER :: ax, ay, az, bx, by, bz, c, cda, cdax, cday, cdaz, cdb, cdbx, &
      cdby, cdbz, coa, coamx, coamy, coamz, coapx, coapy, coapz, cob, cobmx, &
      cobmy, cobmz, cobpx, cobpy, cobpz, da, dax, day, daz, db, dbx, dby, &
      dbz, i, ipgf, j, jk, jpgf, jstart, k, la, la_max, la_min, la_start, lb, &
      lb_max, lb_min, lc, maxder, na, nb, nexp_ppl
    LOGICAL                                  :: calculate_force_a, &
                                                calculate_force_b
    REAL(KIND=dp)                            :: f0, f1, f2, f3, f4, fax, fay, &
                                                faz, fbx, fby, fbz, fc, ftz, &
                                                rcg2, rcp2, zetg, zetp
    REAL(KIND=dp), DIMENSION(3)              :: rag, rbg, rcg, rcp

!   ---------------------------------------------------------------------------

    IF (PRESENT(pab)) THEN
      IF (PRESENT(force_a)) THEN
        calculate_force_a = .TRUE.
        force_a(:) = 0.0_dp
      ELSE
        calculate_force_a = .FALSE.
      END IF
      IF (PRESENT(force_b)) THEN
        calculate_force_b = .TRUE.
        force_b(:) = 0.0_dp
      ELSE
        calculate_force_b = .FALSE.
      END IF
    ELSE
      calculate_force_a = .FALSE.
      calculate_force_b = .FALSE.
    END IF

    la_max = la_max_set + da_max
    la_min = MAX(0,la_min_set-da_max)

    lb_max = lb_max_set + db_max
    lb_min = MAX(0,lb_min_set-db_max)

    nexp_ppl = SIZE(cexp_ppl)
    maxder = ncoset(MAX(da_max,db_max))

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

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
        zetg = 1.0_dp/(zeta(ipgf) + zetb(jpgf) + zetc)

        f0 = (pi*zetg)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetg

        rcp(:) = f1*rab(:) - rac(:)

        rcp2 = rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3)

!       *** Calculate the basic three-center overlap integral [s|s|s] ***

        s(1,1,1,1) = f0*EXP(-(zeta(ipgf)*f1*dab*dab + zetc*zetg*rcp2/zetp))

!       *** Recurrence steps: [s|s|s] -> [a|s|s] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|s|s] -> [a|s|s] ***

          rag(:) = zetg*(zetb(jpgf)*rab(:) + zetc*rac(:))

!         *** [p|s|s] = (Gi - Ai)*[s|s|s]  (i = x,y,z) ***

          s(2,1,1,1) = rag(1)*s(1,1,1,1) ! [px|s|s]
          s(3,1,1,1) = rag(2)*s(1,1,1,1) ! [py|s|s]
          s(4,1,1,1) = rag(3)*s(1,1,1,1) ! [pz|s|s]

          IF (la_max > 1) THEN

!           *** [d|s|s] ***

            f3 = f2*s(1,1,1,1)

            s( 5,1,1,1) = rag(1)*s(2,1,1,1) + f3 ! [dx2|s|s]
            s( 6,1,1,1) = rag(1)*s(3,1,1,1)      ! [dxy|s|s]
            s( 7,1,1,1) = rag(1)*s(4,1,1,1)      ! [dxz|s|s]
            s( 8,1,1,1) = rag(2)*s(3,1,1,1) + f3 ! [dy2|s|s]
            s( 9,1,1,1) = rag(2)*s(4,1,1,1)      ! [dyz|s|s]
            s(10,1,1,1) = rag(3)*s(4,1,1,1) + f3 ! [dz2|s|s]

            IF (la_max > 2) THEN

!             *** [f|s|s] ***

              f3 = 2.0_dp*f2

              s(11,1,1,1) = rag(1)*s( 5,1,1,1) + f3*s(2,1,1,1) ! [fx3 |s|s]
              s(12,1,1,1) = rag(1)*s( 6,1,1,1) + f2*s(3,1,1,1) ! [fx2y|s|s]
              s(13,1,1,1) = rag(1)*s( 7,1,1,1) + f2*s(4,1,1,1) ! [fx2z|s|s]
              s(14,1,1,1) = rag(1)*s( 8,1,1,1)                 ! [fxy2|s|s]
              s(15,1,1,1) = rag(1)*s( 9,1,1,1)                 ! [fxyz|s|s]
              s(16,1,1,1) = rag(1)*s(10,1,1,1)                 ! [fxz2|s|s]
              s(17,1,1,1) = rag(2)*s( 8,1,1,1) + f3*s(3,1,1,1) ! [fy3 |s|s]
              s(18,1,1,1) = rag(2)*s( 9,1,1,1) + f2*s(4,1,1,1) ! [fy2z|s|s]
              s(19,1,1,1) = rag(2)*s(10,1,1,1)                 ! [fyz2|s|s]
              s(20,1,1,1) = rag(3)*s(10,1,1,1) + f3*s(4,1,1,1) ! [fz3 |s|s]

              IF (la_max > 3) THEN

!               *** [g|s|s] ***

                f4 = 3.0_dp*f2

                s(21,1,1,1) = rag(1)*s(11,1,1,1) + f4*s( 5,1,1,1) ! [gx4  |s|s]
                s(22,1,1,1) = rag(1)*s(12,1,1,1) + f3*s( 6,1,1,1) ! [gx3y |s|s]
                s(23,1,1,1) = rag(1)*s(13,1,1,1) + f3*s( 7,1,1,1) ! [gx3z |s|s]
                s(24,1,1,1) = rag(1)*s(14,1,1,1) + f2*s( 8,1,1,1) ! [gx2y2|s|s]
                s(25,1,1,1) = rag(1)*s(15,1,1,1) + f2*s( 9,1,1,1) ! [gx2yz|s|s]
                s(26,1,1,1) = rag(1)*s(16,1,1,1) + f2*s(10,1,1,1) ! [gx2z2|s|s]
                s(27,1,1,1) = rag(1)*s(17,1,1,1)                  ! [gxy3 |s|s]
                s(28,1,1,1) = rag(1)*s(18,1,1,1)                  ! [gxy2z|s|s]
                s(29,1,1,1) = rag(1)*s(19,1,1,1)                  ! [gxyz2|s|s]
                s(30,1,1,1) = rag(1)*s(20,1,1,1)                  ! [gxz3 |s|s]
                s(31,1,1,1) = rag(2)*s(17,1,1,1) + f4*s( 8,1,1,1) ! [gy4  |s|s]
                s(32,1,1,1) = rag(2)*s(18,1,1,1) + f3*s( 9,1,1,1) ! [gy3z |s|s]
                s(33,1,1,1) = rag(2)*s(19,1,1,1) + f2*s(10,1,1,1) ! [gy2z2|s|s]
                s(34,1,1,1) = rag(2)*s(20,1,1,1)                  ! [gyz3 |s|s]
                s(35,1,1,1) = rag(3)*s(20,1,1,1) + f4*s(10,1,1,1) ! [gz4  |s|s]

!               *** [a|s|s] = (Gi - Ai)*[a-1i|s|s] + f2*Ni(a-1i)*[a-2i|s|s] ***

                DO la=5,la_max

!                 *** Increase the angular momentum component z of function a ***

                  s(coset(0,0,la),1,1,1) =&
                    rag(3)*s(coset(0,0,la-1),1,1,1) +&
                    f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,1,1)

!                 *** Increase the angular momentum component y of function a ***

                  az = la - 1
                  s(coset(0,1,az),1,1,1) = rag(2)*s(coset(0,0,az),1,1,1)

                  DO ay=2,la
                    az = la - ay
                    s(coset(0,ay,az),1,1,1) =&
                      rag(2)*s(coset(0,ay-1,az),1,1,1) +&
                      f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,1,1)
                  END DO

!                 *** Increase the angular momentum component x of function a ***

                  DO ay=0,la-1
                    az = la - 1 - ay
                    s(coset(1,ay,az),1,1,1) = rag(1)*s(coset(0,ay,az),1,1,1)
                  END DO

                  DO ax=2,la
                    f3 = f2*REAL(ax-1,dp)
                    DO ay=0,la-ax
                      az = la - ax - ay
                      s(coset(ax,ay,az),1,1,1) =&
                        rag(1)*s(coset(ax-1,ay,az),1,1,1) +&
                        f3*s(coset(ax-2,ay,az),1,1,1)
                    END DO
                  END DO

                END DO

              END IF

            END IF

          END IF

!         *** Recurrence steps: [a|s|s] -> [a|s|b] ***

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_min)
                s(i,j,1,1) = 0.0_dp
              END DO
            END DO

!           *** Horizontal recurrence steps ***

            rbg(:) = rag(:) - rab(:)

!           *** [a|s|p] = [a+1i|s|s] - (Bi - Ai)*[a|s|s] ***

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
                  s(coa,2,1,1) = s(coapx,1,1,1) - rab(1)*s(coa,1,1,1)
                  s(coa,3,1,1) = s(coapy,1,1,1) - rab(2)*s(coa,1,1,1)
                  s(coa,4,1,1) = s(coapz,1,1,1) - rab(3)*s(coa,1,1,1)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|s|p] = (Gi - Bi)*[a|s|s] + f2*Ni(a)*[a-1i|s|s] ***

            DO ax=0,la_max
              fax = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fay = f2*REAL(ay,dp)
                az = la_max - ax - ay
                faz = f2*REAL(az,dp)
                coa = coset(ax,ay,az)
                IF (ax == 0) THEN
                  s(coa,2,1,1) = rbg(1)*s(coa,1,1,1)
                ELSE
                  coamx = coset(ax-1,ay,az)
                  s(coa,2,1,1) = rbg(1)*s(coa,1,1,1) + fax*s(coamx,1,1,1)
                END IF
                IF (ay == 0) THEN
                  s(coa,3,1,1) = rbg(2)*s(coa,1,1,1)
                ELSE
                  coamy = coset(ax,ay-1,az)
                  s(coa,3,1,1) = rbg(2)*s(coa,1,1,1) + fay*s(coamy,1,1,1)
                END IF
                IF (az == 0) THEN
                  s(coa,4,1,1) = rbg(3)*s(coa,1,1,1)
                ELSE
                  coamz = coset(ax,ay,az-1)
                  s(coa,4,1,1) = rbg(3)*s(coa,1,1,1) + faz*s(coamz,1,1,1)
                END IF
              END DO
            END DO

!           *** Recurrence steps: [a|s|p] -> [a|s|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|s|b] = [a+1i|s|b-1i] - (Bi - Ai)*[a|s|b-1i] ***

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay

!                   *** Shift of angular momentum component z from a to b ***

                    s(coset(ax,ay,az),coset(0,0,lb),1,1) =&
                      s(coset(ax,ay,az+1),coset(0,0,lb-1),1,1) -&
                      rab(3)*s(coset(ax,ay,az),coset(0,0,lb-1),1,1)

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      s(coset(ax,ay,az),coset(0,by,bz),1,1) =&
                        s(coset(ax,ay+1,az),coset(0,by-1,bz),1,1) -&
                        rab(2)*s(coset(ax,ay,az),coset(0,by-1,bz),1,1)
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz),1,1) =&
                          s(coset(ax+1,ay,az),coset(bx-1,by,bz),1,1) -&
                          rab(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),1,1)
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|s|b] = (Gi - Bi)*[a|s|b-1i] +   ***
!             ***           f2*Ni(a)*[a-1i|s|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|s|b-2i]   ***

              DO ax=0,la_max
                fax = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fay = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  faz = f2*REAL(az,dp)

!                 *** Shift of angular momentum component z from a to b ***

                  f3 = f2*REAL(lb-1,dp)

                  IF (az == 0) THEN
                    s(coset(ax,ay,az),coset(0,0,lb),1,1) =&
                      rbg(3)*s(coset(ax,ay,az),coset(0,0,lb-1),1,1) +&
                      f3*s(coset(ax,ay,az),coset(0,0,lb-2),1,1)
                  ELSE
                    s(coset(ax,ay,az),coset(0,0,lb),1,1) =&
                      rbg(3)*s(coset(ax,ay,az),coset(0,0,lb-1),1,1) +&
                      faz*s(coset(ax,ay,az-1),coset(0,0,lb-1),1,1) +&
                      f3*s(coset(ax,ay,az),coset(0,0,lb-2),1,1)
                  END IF

!                 *** Shift of angular momentum component y from a to b ***

                  IF (ay == 0) THEN
                    bz = lb - 1
                    s(coset(ax,ay,az),coset(0,1,bz),1,1) =&
                      rbg(2)*s(coset(ax,ay,az),coset(0,0,bz),1,1)
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coset(ax,ay,az),coset(0,by,bz),1,1) =&
                        rbg(2)*s(coset(ax,ay,az),coset(0,by-1,bz),1,1) +&
                        f3*s(coset(ax,ay,az),coset(0,by-2,bz),1,1)
                    END DO
                  ELSE
                    bz = lb - 1
                    s(coset(ax,ay,az),coset(0,1,bz),1,1) =&
                      rbg(2)*s(coset(ax,ay,az),coset(0,0,bz),1,1) +&
                      fay*s(coset(ax,ay-1,az),coset(0,0,bz),1,1)
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coset(ax,ay,az),coset(0,by,bz),1,1) =&
                        rbg(2)*s(coset(ax,ay,az),coset(0,by-1,bz),1,1) +&
                        fay*s(coset(ax,ay-1,az),coset(0,by-1,bz),1,1) +&
                        f3*s(coset(ax,ay,az),coset(0,by-2,bz),1,1)
                    END DO
                  END IF

!                 *** Shift of angular momentum component x from a to b ***

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coset(ax,ay,az),coset(1,by,bz),1,1) =&
                        rbg(1)*s(coset(ax,ay,az),coset(0,by,bz),1,1)
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz),1,1) =&
                          rbg(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),1,1) +&
                          f3*s(coset(ax,ay,az),coset(bx-2,by,bz),1,1)
                      END DO
                    END DO
                  ELSE
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coset(ax,ay,az),coset(1,by,bz),1,1) =&
                        rbg(1)*s(coset(ax,ay,az),coset(0,by,bz),1,1) +&
                        fax*s(coset(ax-1,ay,az),coset(0,by,bz),1,1)
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz),1,1) =&
                          rbg(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),1,1) +&
                          fax*s(coset(ax-1,ay,az),coset(bx-1,by,bz),1,1) +&
                          f3*s(coset(ax,ay,az),coset(bx-2,by,bz),1,1)
                      END DO
                    END DO
                  END IF

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|s|s] -> [s|s|b] ***

            rbg(:) = -zetg*(zeta(ipgf)*rab(:) - zetc*rbc(:))

!           *** [s|s|p] = (Gi - Bi)*[s|s|s] ***

            s(1,2,1,1) = rbg(1)*s(1,1,1,1)
            s(1,3,1,1) = rbg(2)*s(1,1,1,1)
            s(1,4,1,1) = rbg(3)*s(1,1,1,1)

!           *** [s|s|b] = (Gi - Bi)*[s|s|b-1i] + f2*Ni(b-1i)*[s|s|b-2i] ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of function b ***

              s(1,coset(0,0,lb),1,1) =&
                rbg(3)*s(1,coset(0,0,lb-1),1,1) +&
                f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),1,1)

!             *** Increase the angular momentum component y of function b ***

              bz = lb - 1
              s(1,coset(0,1,bz),1,1) = rbg(2)*s(1,coset(0,0,bz),1,1)

              DO by=2,lb
                bz = lb - by
                s(1,coset(0,by,bz),1,1) =&
                  rbg(2)*s(1,coset(0,by-1,bz),1,1) +&
                  f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),1,1)
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                s(1,coset(1,by,bz),1,1) = rbg(1)*s(1,coset(0,by,bz),1,1)
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  s(1,coset(bx,by,bz),1,1) =&
                    rbg(1)*s(1,coset(bx-1,by,bz),1,1) +&
                    f3*s(1,coset(bx-2,by,bz),1,1)
                END DO
              END DO

            END DO

          END IF

        END IF

!       *** Calculate the contributions for lc > 0 ***

        DO c=2,nexp_ppl

          lc = 2*(c - 1)
          fc = f2*REAL(lc,dp)

!         *** Calculate the basic three-center overlap integrals [s|c|s] ***

          SELECT CASE (c)
          CASE (2)
            rcg(:) = -zetg*(zeta(ipgf)*rac(:) + zetb(jpgf)*rbc(:))
            rcg2 = rcg(1)*rcg(1) + rcg(2)*rcg(2) + rcg(3)*rcg(3)
            s(1,1,2,1) = (1.5_dp*zetg + rcg2)*s(1,1,1,1)
          CASE (3)
            s(1,1,3,1) = (2.5_dp*zetg + rcg2)*s(1,1,2,1) + rcg2*zetg*s(1,1,1,1)
          CASE (4)
            s(1,1,4,1) = (3.5_dp*zetg + rcg2)*s(1,1,3,1) +&
                         2.0_dp*rcg2*zetg*(s(1,1,2,1) + zetg*s(1,1,1,1))
          END SELECT

!         *** Recurrence steps: [s|c|s] -> [a|c|s] ***

          IF (la_max > 0) THEN

!           *** Vertical recurrence steps: [s|c|s] -> [a|c|s] ***

!           *** [p|c|s] = (Gi - Ai)*[s|c|s] +                          ***
!           ***           f2*Ni(c)*([p|c-2i|s] - (Ci - Ai)*[s|c-2i|s]) ***

            s(2,1,c,1) = rag(1)*s(1,1,c,1) +&
                         fc*(s(2,1,c-1,1) - rac(1)*s(1,1,c-1,1))
            s(3,1,c,1) = rag(2)*s(1,1,c,1) +&
                         fc*(s(3,1,c-1,1) - rac(2)*s(1,1,c-1,1))
            s(4,1,c,1) = rag(3)*s(1,1,c,1) +&
                         fc*(s(4,1,c-1,1) - rac(3)*s(1,1,c-1,1))

            IF (la_max > 1) THEN

!             *** [d|c|s] ***

              f3 = f2*s(1,1,c,1)

              s( 5,1,c,1) = rag(1)*s(2,1,c,1) + f3 +&
                            fc*(s( 5,1,c-1,1) - rac(1)*s(2,1,c-1,1))
              s( 6,1,c,1) = rag(1)*s(3,1,c,1) +&
                            fc*(s( 6,1,c-1,1) - rac(1)*s(3,1,c-1,1))
              s( 7,1,c,1) = rag(1)*s(4,1,c,1) +&
                            fc*(s( 7,1,c-1,1) - rac(1)*s(4,1,c-1,1))
              s( 8,1,c,1) = rag(2)*s(3,1,c,1) + f3 +&
                            fc*(s( 8,1,c-1,1) - rac(2)*s(3,1,c-1,1))
              s( 9,1,c,1) = rag(2)*s(4,1,c,1) +&
                            fc*(s( 9,1,c-1,1) - rac(2)*s(4,1,c-1,1))
              s(10,1,c,1) = rag(3)*s(4,1,c,1) + f3 +&
                            fc*(s(10,1,c-1,1) - rac(3)*s(4,1,c-1,1))

              IF (la_max > 2) THEN

!               *** [f|c|s] ***

                f3 = 2.0_dp*f2

                s(11,1,c,1) = rag(1)*s( 5,1,c,1) + f3*s(2,1,c,1) +&
                              fc*(s(11,1,c-1,1) - rac(1)*s( 5,1,c-1,1))
                s(12,1,c,1) = rag(1)*s( 6,1,c,1) + f2*s(3,1,c,1) +&
                              fc*(s(12,1,c-1,1) - rac(1)*s( 6,1,c-1,1))
                s(13,1,c,1) = rag(1)*s( 7,1,c,1) + f2*s(4,1,c,1) +&
                              fc*(s(13,1,c-1,1) - rac(1)*s( 7,1,c-1,1))
                s(14,1,c,1) = rag(1)*s( 8,1,c,1) +&
                              fc*(s(14,1,c-1,1) - rac(1)*s( 8,1,c-1,1))
                s(15,1,c,1) = rag(1)*s( 9,1,c,1) +&
                              fc*(s(15,1,c-1,1) - rac(1)*s( 9,1,c-1,1))
                s(16,1,c,1) = rag(1)*s(10,1,c,1) +&
                              fc*(s(16,1,c-1,1) - rac(1)*s(10,1,c-1,1))
                s(17,1,c,1) = rag(2)*s( 8,1,c,1) + f3*s(3,1,c,1) +&
                              fc*(s(17,1,c-1,1) - rac(2)*s( 8,1,c-1,1))
                s(18,1,c,1) = rag(2)*s( 9,1,c,1) + f2*s(4,1,c,1) +&
                              fc*(s(18,1,c-1,1) - rac(2)*s( 9,1,c-1,1))
                s(19,1,c,1) = rag(2)*s(10,1,c,1) +&
                              fc*(s(19,1,c-1,1) - rac(2)*s(10,1,c-1,1))
                s(20,1,c,1) = rag(3)*s(10,1,c,1) + f3*s(4,1,c,1) +&
                              fc*(s(20,1,c-1,1) - rac(3)*s(10,1,c-1,1))

                IF (la_max > 3) THEN

!                 *** [g|c|s] ***

                  f4 = 3.0_dp*f2

                  s(21,1,c,1) = rag(1)*s(11,1,c,1) + f4*s( 5,1,c,1) +&
                                fc*(s(21,1,c-1,1) - rac(1)*s(11,1,c-1,1))
                  s(22,1,c,1) = rag(1)*s(12,1,c,1) + f3*s( 6,1,c,1) +&
                                fc*(s(22,1,c-1,1) - rac(1)*s(12,1,c-1,1))
                  s(23,1,c,1) = rag(1)*s(13,1,c,1) + f3*s( 7,1,c,1) +&
                                fc*(s(23,1,c-1,1) - rac(1)*s(13,1,c-1,1))
                  s(24,1,c,1) = rag(1)*s(14,1,c,1) + f2*s( 8,1,c,1) +&
                                fc*(s(24,1,c-1,1) - rac(1)*s(14,1,c-1,1))
                  s(25,1,c,1) = rag(1)*s(15,1,c,1) + f2*s( 9,1,c,1) +&
                                fc*(s(25,1,c-1,1) - rac(1)*s(15,1,c-1,1))
                  s(26,1,c,1) = rag(1)*s(16,1,c,1) + f2*s(10,1,c,1) +&
                                fc*(s(26,1,c-1,1) - rac(1)*s(16,1,c-1,1))
                  s(27,1,c,1) = rag(1)*s(17,1,c,1) +&
                                fc*(s(27,1,c-1,1) - rac(1)*s(17,1,c-1,1))
                  s(28,1,c,1) = rag(1)*s(18,1,c,1) +&
                                fc*(s(28,1,c-1,1) - rac(1)*s(18,1,c-1,1))
                  s(29,1,c,1) = rag(1)*s(19,1,c,1) +&
                                fc*(s(29,1,c-1,1) - rac(1)*s(19,1,c-1,1))
                  s(30,1,c,1) = rag(1)*s(20,1,c,1) +&
                                fc*(s(30,1,c-1,1) - rac(1)*s(20,1,c-1,1))
                  s(31,1,c,1) = rag(2)*s(17,1,c,1) + f4*s( 8,1,c,1) +&
                                fc*(s(31,1,c-1,1) - rac(2)*s(17,1,c-1,1))
                  s(32,1,c,1) = rag(2)*s(18,1,c,1) + f3*s( 9,1,c,1) +&
                                fc*(s(32,1,c-1,1) - rac(2)*s(18,1,c-1,1))
                  s(33,1,c,1) = rag(2)*s(19,1,c,1) + f2*s(10,1,c,1) +&
                                fc*(s(33,1,c-1,1) - rac(2)*s(19,1,c-1,1))
                  s(34,1,c,1) = rag(2)*s(20,1,c,1) +&
                                fc*(s(34,1,c-1,1) - rac(2)*s(20,1,c-1,1))
                  s(35,1,c,1) = rag(3)*s(20,1,c,1) + f4*s(10,1,c,1) +&
                                fc*(s(35,1,c-1,1) - rac(3)*s(20,1,c-1,1))

!                 *** [a|c|s] = (Gi - Ai)*[a-1i|c|s] +             ***
!                 ***           f2*Ni(a-1i)*[a-2i|c|s] +           ***
!                 ***           f2*Ni(c)*([a|c-2i|s] -             ***
!                 ***                     (Ci - Ai)*[a-1i|c-2i|s]) ***

                  DO la=5,la_max

!                   *** Increase the angular momentum component z of function a ***

                    s(coset(0,0,la),1,c,1) =&
                      rag(3)*s(coset(0,0,la-1),1,c,1) +&
                      f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,c,1) +&
                      fc*(s(coset(0,0,la),1,c-1,1) -&
                          rac(3)*s(coset(0,0,la-1),1,c-1,1))

!                   *** Increase the angular momentum component y of function a ***

                    az = la - 1
                    s(coset(0,1,az),1,c,1) =&
                      rag(2)*s(coset(0,0,az),1,c,1) +&
                      fc*(s(coset(0,1,az),1,c-1,1) -&
                          rac(2)*s(coset(0,0,az),1,c-1,1))

                    DO ay=2,la
                      az = la - ay
                      s(coset(0,ay,az),1,c,1) =&
                        rag(2)*s(coset(0,ay-1,az),1,c,1) +&
                        f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,c,1) +&
                        fc*(s(coset(0,ay,az),1,c-1,1) -&
                            rac(2)*s(coset(0,ay-1,az),1,c-1,1))
                    END DO

!                   *** Increase the angular momentum component x of function a ***

                    DO ay=0,la-1
                      az = la - 1 - ay
                      s(coset(1,ay,az),1,c,1) =&
                        rag(1)*s(coset(0,ay,az),1,c,1) +&
                        fc*(s(coset(1,ay,az),1,c-1,1) -&
                            rac(1)*s(coset(0,ay,az),1,c-1,1))
                    END DO

                    DO ax=2,la
                      f3 = f2*REAL(ax-1,dp)
                      DO ay=0,la-ax
                        az = la - ax - ay
                        s(coset(ax,ay,az),1,c,1) =&
                          rag(1)*s(coset(ax-1,ay,az),1,c,1) +&
                          f3*s(coset(ax-2,ay,az),1,c,1) +&
                          fc*(s(coset(ax,ay,az),1,c-1,1) -&
                              rac(1)*s(coset(ax-1,ay,az),1,c-1,1))
                      END DO
                    END DO

                  END DO

                END IF

              END IF

            END IF

!           *** Recurrence steps: [a|c|s] -> [a|c|b] ***

            IF (lb_max > 0) THEN

              DO j=2,ncoset(lb_max)
                DO i=1,ncoset(la_min)
                  s(i,j,c,1) = 0.0_dp
                END DO
              END DO

!             *** Horizontal recurrence steps ***

!             *** [a|c|p] = [a+1i|c|s] - (Bi - Ai)*[a|c|s] ***

              IF (lb_max == 1) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    s(coset(ax,ay,az),2,c,1) =&
                      s(coset(ax+1,ay,az),1,c,1) -&
                      rab(1)*s(coset(ax,ay,az),1,c,1)
                    s(coset(ax,ay,az),3,c,1) =&
                      s(coset(ax,ay+1,az),1,c,1) -&
                      rab(2)*s(coset(ax,ay,az),1,c,1)
                    s(coset(ax,ay,az),4,c,1) =&
                      s(coset(ax,ay,az+1),1,c,1) -&
                      rab(3)*s(coset(ax,ay,az),1,c,1)
                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|c|p] = (Gi - Bi)*[a|c|s] +                          ***
!             ***           f2*Ni(a)*[a-1i|c|s] +                        ***
!             ***           f2*Ni(c)*([a|c-2i|p] - (Ci - Bi)*[a|c-2i|s]) ***

              DO ax=0,la_max
                fax = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fay = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  faz = f2*REAL(az,dp)
                  IF (ax == 0) THEN
                    s(coset(ax,ay,az),2,c,1) =&
                      rbg(1)*s(coset(ax,ay,az),1,c,1) +&
                      fc*(s(coset(ax,ay,az),2,c-1,1) -&
                          rbc(1)*s(coset(ax,ay,az),1,c-1,1))
                  ELSE
                    s(coset(ax,ay,az),2,c,1) =&
                      rbg(1)*s(coset(ax,ay,az),1,c,1) +&
                      fax*s(coset(ax-1,ay,az),1,c,1) +&
                      fc*(s(coset(ax,ay,az),2,c-1,1) -&
                          rbc(1)*s(coset(ax,ay,az),1,c-1,1))
                  END IF
                  IF (ay == 0) THEN
                    s(coset(ax,ay,az),3,c,1) =&
                      rbg(2)*s(coset(ax,ay,az),1,c,1) +&
                      fc*(s(coset(ax,ay,az),3,c-1,1) -&
                          rbc(2)*s(coset(ax,ay,az),1,c-1,1))
                  ELSE
                    s(coset(ax,ay,az),3,c,1) =&
                      rbg(2)*s(coset(ax,ay,az),1,c,1) +&
                      fay*s(coset(ax,ay-1,az),1,c,1) +&
                      fc*(s(coset(ax,ay,az),3,c-1,1) -&
                          rbc(2)*s(coset(ax,ay,az),1,c-1,1))
                  END IF
                  IF (az == 0) THEN
                    s(coset(ax,ay,az),4,c,1) =&
                      rbg(3)*s(coset(ax,ay,az),1,c,1) +&
                      fc*(s(coset(ax,ay,az),4,c-1,1) -&
                          rbc(3)*s(coset(ax,ay,az),1,c-1,1))
                  ELSE
                    s(coset(ax,ay,az),4,c,1) =&
                      rbg(3)*s(coset(ax,ay,az),1,c,1) +&
                      faz*s(coset(ax,ay,az-1),1,c,1) +&
                      fc*(s(coset(ax,ay,az),4,c-1,1) -&
                          rbc(3)*s(coset(ax,ay,az),1,c-1,1))
                  END IF
                END DO
              END DO

!             *** Recurrence steps: [a|s|p] -> [a|s|b] ***

              DO lb=2,lb_max

!               *** Horizontal recurrence steps ***

!               *** [a|c|b] = [a+1i|c|b-1i] - (Bi - Ai)*[a|c|b-1i] ***

                IF (lb == lb_max) THEN
                  la_start = la_min
                ELSE
                  la_start = MAX(0,la_min-1)
                END IF

                DO la=la_start,la_max-1
                  DO ax=0,la
                    DO ay=0,la-ax
                      az = la - ax - ay

!                     *** Shift of angular momentum component x from a to b ***

                      s(coset(ax,ay,az),coset(0,0,lb),c,1) =&
                        s(coset(ax,ay,az+1),coset(0,0,lb-1),c,1) -&
                        rab(3)*s(coset(ax,ay,az),coset(0,0,lb-1),c,1)

!                     *** Shift of angular momentum component y from a to b ***

                      DO by=1,lb
                        bz = lb - by
                        s(coset(ax,ay,az),coset(0,by,bz),c,1) =&
                          s(coset(ax,ay+1,az),coset(0,by-1,bz),c,1) -&
                          rab(2)*s(coset(ax,ay,az),coset(0,by-1,bz),c,1)
                      END DO

!                     *** Shift of angular momentum component z from a to b ***

                      DO bx=1,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(coset(ax,ay,az),coset(bx,by,bz),c,1) =&
                            s(coset(ax+1,ay,az),coset(bx-1,by,bz),c,1) -&
                            rab(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),c,1)
                        END DO
                      END DO

                    END DO
                  END DO
                END DO

!               *** Vertical recurrence step ***

!               *** [a|c|b] = (Gi - Bi)*[a|c|b-1i] +             ***
!               ***           f2*Ni(a)*[a-1i|c|b-1i] +           ***
!               ***           f2*Ni(b-1i)*[a|c|b-2i] +           ***
!               ***           f2*Ni(c)*([a|c-2i|b] -             ***
!               ***                     (Ci - Bi)*[a|c-2i|b-1i]) ***

                DO ax=0,la_max
                  fax = f2*REAL(ax,dp)
                  DO ay=0,la_max-ax
                    fay = f2*REAL(ay,dp)
                    az = la_max - ax - ay
                    faz = f2*REAL(az,dp)

!                   *** Shift of angular momentum component z from a to b ***

                    f3 = f2*REAL(lb-1,dp)

                    IF (az == 0) THEN
                      s(coset(ax,ay,az),coset(0,0,lb),c,1) =&
                        rbg(3)*s(coset(ax,ay,az),coset(0,0,lb-1),c,1) +&
                        f3*s(coset(ax,ay,az),coset(0,0,lb-2),c,1) +&
                        fc*(s(coset(ax,ay,az),coset(0,0,lb),c-1,1) -&
                            rbc(3)*s(coset(ax,ay,az),coset(0,0,lb-1),c-1,1))
                    ELSE
                      s(coset(ax,ay,az),coset(0,0,lb),c,1) =&
                        rbg(3)*s(coset(ax,ay,az),coset(0,0,lb-1),c,1) +&
                        faz*s(coset(ax,ay,az-1),coset(0,0,lb-1),c,1) +&
                        f3*s(coset(ax,ay,az),coset(0,0,lb-2),c,1) +&
                        fc*(s(coset(ax,ay,az),coset(0,0,lb),c-1,1) -&
                            rbc(3)*s(coset(ax,ay,az),coset(0,0,lb-1),c-1,1))
                    END IF

!                   *** Shift of angular momentum component y from a to b ***

                    IF (ay == 0) THEN
                      bz = lb - 1
                      s(coset(ax,ay,az),coset(0,1,bz),c,1) =&
                        rbg(2)*s(coset(ax,ay,az),coset(0,0,bz),c,1) +&
                        fc*(s(coset(ax,ay,az),coset(0,1,bz),c-1,1) -&
                            rbc(2)*s(coset(ax,ay,az),coset(0,0,bz),c-1,1))
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        s(coset(ax,ay,az),coset(0,by,bz),c,1) =&
                          rbg(2)*s(coset(ax,ay,az),coset(0,by-1,bz),c,1) +&
                          f3*s(coset(ax,ay,az),coset(0,by-2,bz),c,1) +&
                          fc*(s(coset(ax,ay,az),coset(0,by,bz),c-1,1) -&
                              rbc(2)*s(coset(ax,ay,az),&
                                          coset(0,by-1,bz),c-1,1))
                      END DO
                    ELSE
                      bz = lb - 1
                      s(coset(ax,ay,az),coset(0,1,bz),c,1) =&
                        rbg(2)*s(coset(ax,ay,az),coset(0,0,bz),c,1) +&
                        fay*s(coset(ax,ay-1,az),coset(0,0,bz),c,1) +&
                        fc*(s(coset(ax,ay,az),coset(0,1,bz),c-1,1) -&
                            rbc(2)*s(coset(ax,ay,az),coset(0,0,bz),c-1,1))
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        s(coset(ax,ay,az),coset(0,by,bz),c,1) =&
                          rbg(2)*s(coset(ax,ay,az),coset(0,by-1,bz),c,1) +&
                          fay*s(coset(ax,ay-1,az),coset(0,by-1,bz),c,1) +&
                          f3*s(coset(ax,ay,az),coset(0,by-2,bz),c,1) +&
                          fc*(s(coset(ax,ay,az),coset(0,by,bz),c-1,1) -&
                              rbc(2)*s(coset(ax,ay,az),&
                                          coset(0,by-1,bz),c-1,1))
                      END DO
                    END IF

!                   *** Shift of angular momentum component x from a to b ***

                    IF (ax == 0) THEN
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        s(coset(ax,ay,az),coset(1,by,bz),c,1) =&
                          rbg(1)*s(coset(ax,ay,az),coset(0,by,bz),c,1) +&
                          fc*(s(coset(ax,ay,az),coset(1,by,bz),c-1,1) -&
                              rbc(1)*s(coset(ax,ay,az),coset(0,by,bz),c-1,1))
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(coset(ax,ay,az),coset(bx,by,bz),c,1) =&
                            rbg(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),c,1) +&
                            f3*s(coset(ax,ay,az),coset(bx-2,by,bz),c,1) +&
                            fc*(s(coset(ax,ay,az),coset(bx,by,bz),c-1,1) -&
                                rbc(1)*s(coset(ax,ay,az),&
                                            coset(bx-1,by,bz),c-1,1))
                        END DO
                      END DO
                    ELSE
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        s(coset(ax,ay,az),coset(1,by,bz),c,1) =&
                          rbg(1)*s(coset(ax,ay,az),coset(0,by,bz),c,1) +&
                          fax*s(coset(ax-1,ay,az),coset(0,by,bz),c,1) +&
                          fc*(s(coset(ax,ay,az),coset(1,by,bz),c-1,1) -&
                                rbc(1)*s(coset(ax,ay,az),&
                                            coset(0,by,bz),c-1,1))
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(coset(ax,ay,az),coset(bx,by,bz),c,1) =&
                            rbg(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),c,1) +&
                            fax*s(coset(ax-1,ay,az),coset(bx-1,by,bz),c,1) +&
                            f3*s(coset(ax,ay,az),coset(bx-2,by,bz),c,1) +&
                            fc*(s(coset(ax,ay,az),coset(bx,by,bz),c-1,1) -&
                                rbc(1)*s(coset(ax,ay,az),&
                                            coset(bx-1,by,bz),c-1,1))
                        END DO
                      END DO
                    END IF

                  END DO
                END DO

              END DO

            END IF

          ELSE

            IF (lb_max > 0) THEN

!             *** Vertical recurrence steps: [s|c|s] -> [s|c|b] ***

!             *** [s|c|p] = (Gi - Bi)*[s|c|s] +                          ***
!             ***           f2*Ni(c)*([s|c-2i|p] - (Ci - Bi)*[s|c-2i|s]) ***

              s(1,2,c,1) = rbg(1)*s(1,1,c,1) +&
                            fc*(s(1,2,c-1,1) - rbc(1)*s(1,1,c-1,1))
              s(1,3,c,1) = rbg(2)*s(1,1,c,1) +&
                            fc*(s(1,3,c-1,1) - rbc(2)*s(1,1,c-1,1))
              s(1,4,c,1) = rbg(3)*s(1,1,c,1) +&
                            fc*(s(1,4,c-1,1) - rbc(3)*s(1,1,c-1,1))

!             *** [s|c|b] = (Gi - Bi)*[s|c|b-1i] +                          ***
!             ***           f2*Ni(b-1i)*[s|c|b-2i] +                        ***
!             ***           f2*Ni(c)*([s|c-2i|b] - (Ci - Bi)*[s|c-2i|b-1i]) ***

              DO lb=2,lb_max

!               *** Increase the angular momentum component z of function b ***

                s(1,coset(0,0,lb),c,1) =&
                  rbg(3)*s(1,coset(0,0,lb-1),c,1) +&
                  f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),c,1) +&
                  fc*(s(1,coset(0,0,lb),c-1,1) -&
                      rbc(3)*s(1,coset(0,0,lb-1),c-1,1))

!               *** Increase the angular momentum component y of function b ***

                bz = lb - 1
                s(1,coset(0,1,bz),c,1) =&
                  rbg(2)*s(1,coset(0,0,bz),c,1) +&
                  fc*(s(1,coset(0,1,bz),c-1,1) -&
                      rbc(2)*s(1,coset(0,0,bz),c-1,1))

                DO by=2,lb
                  bz = lb - by
                  s(1,coset(0,by,bz),c,1) =&
                    rbg(2)*s(1,coset(0,by-1,bz),c,1) +&
                    f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),c,1) +&
                    fc*(s(1,coset(0,by,bz),c-1,1) -&
                        rbc(2)*s(1,coset(0,by-1,bz),c-1,1))
                END DO

!               *** Increase the angular momentum component x of function b ***

                DO by=0,lb-1
                  bz = lb - 1 - by
                  s(1,coset(1,by,bz),c,1) =&
                    rbg(1)*s(1,coset(0,by,bz),c,1) +&
                    fc*(s(1,coset(1,by,bz),c-1,1) -&
                        rbc(1)*s(1,coset(0,by,bz),c-1,1))
                END DO

                DO bx=2,lb
                  f3 = f2*REAL(bx-1,dp)
                  DO by=0,lb-bx
                    bz = lb - bx - by
                    s(1,coset(bx,by,bz),c,1) =&
                      rbg(1)*s(1,coset(bx-1,by,bz),c,1) +&
                      f3*s(1,coset(bx-2,by,bz),c,1) +&
                      fc*(s(1,coset(bx,by,bz),c-1,1) -&
                          rbc(1)*s(1,coset(bx-1,by,bz),c-1,1))
                  END DO
                END DO

              END DO

            END IF

          END IF

        END DO

!       *** Store the primitive three-center overlap integrals ***

        DO c=1,nexp_ppl
          DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
            DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
              vab(na+i,nb+j) = vab(na+i,nb+j) + cexp_ppl(c)*s(i,j,c,1)
            END DO
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

              DO c=1,nexp_ppl
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
                      DO lb=lb_min_set,lb_max_set
                        DO bx=0,lb
                          DO by=0,lb-bx
                            bz = lb - bx - by
                            cob = coset(bx,by,bz)
                            s(coa,cob,c,cdax) = ftz*s(coapx,cob,c,cda) -&
                                                fax*s(coamx,cob,c,cda)
                            s(coa,cob,c,cday) = ftz*s(coapy,cob,c,cda) -&
                                                fay*s(coamy,cob,c,cda)
                            s(coa,cob,c,cdaz) = ftz*s(coapz,cob,c,cda) -&
                                                faz*s(coamz,cob,c,cda)
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

!       *** Return all the calculated derivatives of the primitive ***
!       *** three-center overlap integrals w.r.t. atomic center a, ***
!       *** if requested                                           ***

        IF (return_derivatives) THEN
          DO k=2,ncoset(da_max)
            jstart = (k - 1)*SIZE(vab,1)
            DO c=1,nexp_ppl
              DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
                jk = jstart + j
                DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                  vab(na+i,nb+jk) = vab(na+i,nb+jk) + cexp_ppl(c)*s(i,j,c,k)
                END DO
              END DO
            END DO
          END DO
        END IF

!       *** Calculate the force contribution for the atomic center a ***

        IF (calculate_force_a) THEN
          DO k=1,3
            DO c=1,nexp_ppl
              DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
                DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                  force_a(k) = force_a(k) + pab(na+i,nb+j)*&
                                            cexp_ppl(c)*s(i,j,c,k+1)
                END DO
              END DO
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

              DO c=1,nexp_ppl
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
                      DO la=la_min_set,la_max_set
                        DO ax=0,la
                          DO ay=0,la-ax
                            az = la - ax - ay
                            coa = coset(ax,ay,az)
                            s(coa,cob,c,cdbx) = ftz*s(coa,cobpx,c,cdb) -&
                                                fbx*s(coa,cobmx,c,cdb)
                            s(coa,cob,c,cdby) = ftz*s(coa,cobpy,c,cdb) -&
                                                fby*s(coa,cobmy,c,cdb)
                            s(coa,cob,c,cdbz) = ftz*s(coa,cobpz,c,cdb) -&
                                                fbz*s(coa,cobmz,c,cdb)
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

!       *** Return all the calculated derivatives of the primitive ***
!       *** three-center overlap integrals w.r.t. atomic center b, ***
!       *** if requested                                           ***

        IF (return_derivatives) THEN
          DO k=2,ncoset(db_max)
            jstart = (ncoset(da_max) + k - 2)*SIZE(vab,1)
            DO c=1,nexp_ppl
              DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
                jk = jstart + j
                DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                  vab(na+i,nb+jk) = vab(na+i,nb+jk) + cexp_ppl(c)*s(i,j,c,k)
                END DO
              END DO
            END DO
          END DO
        END IF

!       *** Calculate the force contribution for the atomic center b ***

        IF (calculate_force_b) THEN
          DO k=1,3
            DO c=1,nexp_ppl
              DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
                DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                  force_b(k) = force_b(k) + pab(na+i,nb+j)*&
                                            cexp_ppl(c)*s(i,j,c,k+1)
                END DO
              END DO
            END DO
          END DO
        END IF

        nb = nb + ncoset(lb_max_set)

      END DO

      na = na + ncoset(la_max_set)

    END DO

  END SUBROUTINE overlap_ppl

! *****************************************************************************

END MODULE ai_overlap_ppl
